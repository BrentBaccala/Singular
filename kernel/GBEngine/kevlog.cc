/**
 * @file kevlog.cc
 * @brief Global event log implementation (task 325, refactored by
 *        task parallel-bba-event-log-pcopy).
 *
 * See kevlog.h for the record schema + semantics.  Summary:
 *   - Preallocate 1 GiB buffer at bba entry (when gated).
 *   - Preallocate 256 MiB aux arena.
 *   - Global atomic seq; each event reserves its slot via fetch_add.
 *   - Overflow => abort.
 *   - Dump to disk only on SINGULAR_CHECK_IDEAL_MEMBERSHIP violation.
 *   - Poly capture: every event-captured poly is `p_Copy`d into a
 *     source->copy map at emit time.  Dumps render LMs from the
 *     copies (guaranteed-stable content).  Shutdown p_Delete's every
 *     copy.
 */

#include "kernel/GBEngine/kevlog.h"
#include "kernel/GBEngine/kthread.h"  // kt_lm_str (used at dump time only)
#include "kernel/polys.h"
#include "polys/monomials/p_polys.h"
#include "polys/monomials/ring.h"
#include "omalloc/omalloc.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <pthread.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

/* ------------------------------------------------------------------ */
/*  Gate + allocation                                                  */
/* ------------------------------------------------------------------ */
bool g_event_log_enabled = false;

static constexpr size_t EVLOG_CAP_BYTES  = (size_t)1 << 30;   // 1 GiB
static constexpr size_t EVLOG_CAP_RECS   = EVLOG_CAP_BYTES / sizeof(KEvtRecord);
static constexpr size_t EVLOG_AUX_BYTES  = (size_t)256 << 20; // 256 MiB

static KEvtRecord       *g_buf = nullptr;
static unsigned char    *g_aux = nullptr;
static std::atomic<uint64_t> g_seq{0};
static std::atomic<uint64_t> g_aux_off{1};  // reserve offset 0 for "no aux"
static int                g_dispatch_id = -1;

/* Every kevlog_capture call p_Copy's a fresh snapshot of the source
 * poly AT THE MOMENT of capture and appends one entry to this vector.
 * Rationale: ksReducePoly mutates ap->P's exponent vector in place
 * (task parallel-bba-raw-exp-probe / 336 uncovered that an earlier
 * first-sight-caching design produced stale LMs for mutated polys,
 * which invalidated the checker's LM-divides invariant-7 logic).
 * Consumers that need "same source poly across events" use the
 * source_ptr column in -polys.txt (copy -> source lookup). */
struct CapturedPoly {
  poly copy;            // fresh p_Copy at the moment of capture
  const void *source;   // raw source pointer at capture time
                        // (may be stale/freed post-capture)
  ring lmRing;          // ring used by p_Delete
  ring tailRing;        // separate tailRing, or same as lmRing
};
static std::mutex g_captured_lock;
static std::vector<CapturedPoly> *g_captured = nullptr;

/* Set of captured COPY pointers marked as final-S candidates (poly_ptr_1
 * of an EVT_ENTERS event).  Full p_String(p) is dumped into
 * -full-polys.txt at failure time for each entry. */
static std::mutex g_enters_set_lock;
static std::unordered_set<const void *> *g_enters_set = nullptr;

/* ------------------------------------------------------------------ */
/*  Abort helper                                                       */
/* ------------------------------------------------------------------ */
static void kevlog_overflow_abort(const char *what, uint64_t where)
{
  fprintf(stderr,
    "*** kevlog OVERFLOW: %s at %lu (seq cap=%zu, aux cap=%zu)\n"
    "*** Event log buffer full — raise capacity or shorten the run.\n",
    what, (unsigned long)where, EVLOG_CAP_RECS, EVLOG_AUX_BYTES);
  fflush(stderr);
  abort();
}

/* ------------------------------------------------------------------ */
/*  init / shutdown                                                    */
/* ------------------------------------------------------------------ */
void kevlog_init(int disp_id)
{
  // Gate: SINGULAR_EVENT_LOG=1 (any nonempty value enables).
  const char *s = getenv("SINGULAR_EVENT_LOG");
  g_event_log_enabled = (s != NULL && s[0] != 0 && s[0] != '0');
  if (!g_event_log_enabled) return;

  // Idempotent: if already initialised, just bump disp id.
  g_dispatch_id = disp_id;
  if (g_buf != nullptr) return;

  g_buf = (KEvtRecord *)calloc(EVLOG_CAP_RECS, sizeof(KEvtRecord));
  if (g_buf == nullptr) {
    fprintf(stderr, "*** kevlog: failed to calloc %zu bytes for event log\n",
            EVLOG_CAP_BYTES);
    g_event_log_enabled = false;
    return;
  }
  g_aux = (unsigned char *)calloc(EVLOG_AUX_BYTES, 1);
  if (g_aux == nullptr) {
    fprintf(stderr, "*** kevlog: failed to calloc %zu bytes for aux arena\n",
            EVLOG_AUX_BYTES);
    free(g_buf); g_buf = nullptr;
    g_event_log_enabled = false;
    return;
  }
  g_seq.store(0, std::memory_order_release);
  g_aux_off.store(1, std::memory_order_release);

  {
    std::lock_guard<std::mutex> lk(g_captured_lock);
    if (g_captured == nullptr)
      g_captured = new std::vector<CapturedPoly>();
    else
      g_captured->clear();
  }
  {
    std::lock_guard<std::mutex> lk(g_enters_set_lock);
    if (g_enters_set == nullptr)
      g_enters_set = new std::unordered_set<const void *>();
    else
      g_enters_set->clear();
  }

  fprintf(stderr, "[kevlog] init disp=%d seq_cap=%zu aux_cap=%zu\n",
          disp_id, EVLOG_CAP_RECS, EVLOG_AUX_BYTES);
  fflush(stderr);
}

/* Delete every captured poly-copy.  Called from shutdown; dump must
 * have already run (it needs the copies live to format LMs). */
static void kevlog_delete_captured_polys()
{
  std::vector<CapturedPoly> *v = nullptr;
  {
    std::lock_guard<std::mutex> lk(g_captured_lock);
    v = g_captured;
    g_captured = nullptr;  // detach so concurrent captures can't use it
  }
  if (v == nullptr) return;
  for (auto &e : *v)
  {
    poly p = e.copy;
    if (p == nullptr) continue;
    if (e.tailRing != nullptr && e.tailRing != e.lmRing)
      p_Delete(&p, e.lmRing, e.tailRing);
    else if (e.lmRing != nullptr)
      p_Delete(&p, e.lmRing);
  }
  delete v;
}

void kevlog_shutdown()
{
  if (!g_event_log_enabled) return;

  // Delete every captured poly-copy before freeing the buffer.
  kevlog_delete_captured_polys();

  // Drop the enters-set.  The pointers in it are copies already freed
  // by kevlog_delete_captured_polys above, so we just delete the set
  // itself without dereferencing.
  {
    std::lock_guard<std::mutex> lk(g_enters_set_lock);
    if (g_enters_set != nullptr) { delete g_enters_set; g_enters_set = nullptr; }
  }

  if (g_buf) { free(g_buf); g_buf = nullptr; }
  if (g_aux) { free(g_aux); g_aux = nullptr; }

  g_seq.store(0, std::memory_order_release);
  g_aux_off.store(1, std::memory_order_release);
  g_event_log_enabled = false;
}

/* ------------------------------------------------------------------ */
/*  Event emission                                                     */
/* ------------------------------------------------------------------ */
void kevlog_emit(uint16_t type, uint16_t tid, uint16_t atT, uint16_t flags,
                 uint32_t arg_a, uint32_t arg_b, uint32_t arg_c,
                 uint32_t aux_off,
                 const void *poly_ptr_1, const void *poly_ptr_2)
{
  if (!g_event_log_enabled || g_buf == nullptr) return;
  uint64_t my = g_seq.fetch_add(1, std::memory_order_acq_rel);
  if (my >= EVLOG_CAP_RECS)
    kevlog_overflow_abort("seq", my);
  KEvtRecord r;
  r.seq         = my;
  r.type        = type;
  r.tid         = tid;
  r.atT         = atT;
  r.flags       = flags;
  r.arg_a       = arg_a;
  r.arg_b       = arg_b;
  r.arg_c       = arg_c;
  r.aux_off     = aux_off;
  r.poly_ptr_1  = (uint64_t)poly_ptr_1;
  r.poly_ptr_2  = (uint64_t)poly_ptr_2;
  memcpy(&g_buf[my], &r, sizeof(r));
}

uint32_t kevlog_aux_alloc(size_t n, void **out_ptr)
{
  if (!g_event_log_enabled || g_aux == nullptr) {
    if (out_ptr) *out_ptr = nullptr;
    return 0;
  }
  if (n == 0) {
    if (out_ptr) *out_ptr = nullptr;
    return 0;
  }
  // 8-byte align to keep reads sane.
  size_t aligned = (n + 7) & ~((size_t)7);
  uint64_t off = g_aux_off.fetch_add(aligned, std::memory_order_acq_rel);
  if (off + aligned > EVLOG_AUX_BYTES)
    kevlog_overflow_abort("aux", off);
  if (out_ptr) *out_ptr = g_aux + off;
  return (uint32_t)off;
}

/* ------------------------------------------------------------------ */
/*  Poly capture (p_Copy under mutex, append-only — no source caching) */
/* ------------------------------------------------------------------ */
const void *kevlog_capture(const void *src, struct ip_sring *r)
{
  if (!g_event_log_enabled || src == nullptr || r == nullptr) return nullptr;
  // Always p_Copy fresh — every call yields a snapshot of `src` AT THE
  // MOMENT of this call.  ksReducePoly mutates ap->P's exponents in
  // place, so any source-caching scheme misrepresents later events.
  // p_Copy bottoms out in malloc (OMALLOC_USES_MALLOC=1 in this build),
  // so it's thread-safe.  The source's tail chain may be mid-mutation
  // on another thread; this is the same hazard the existing tracers
  // face and is observed to be benign for the 0.2s reproducer.
  std::lock_guard<std::mutex> lk(g_captured_lock);
  if (g_captured == nullptr) return nullptr;
  poly copy = p_Copy((poly)src, (ring)r);
  if (copy == nullptr) return nullptr;
  CapturedPoly e;
  e.copy     = copy;
  e.source   = src;
  e.lmRing   = (ring)r;
  e.tailRing = (ring)r;
  g_captured->push_back(e);
  return (const void *)copy;
}

const void *kevlog_capture_with_tail(const void *src,
                                     struct ip_sring *lmRing,
                                     struct ip_sring *tailRing)
{
  if (!g_event_log_enabled || src == nullptr || lmRing == nullptr) return nullptr;
  // Fall back to single-ring copy when the tailRing is identical.
  if (tailRing == nullptr || tailRing == lmRing)
    return kevlog_capture(src, lmRing);
  std::lock_guard<std::mutex> lk(g_captured_lock);
  if (g_captured == nullptr) return nullptr;
  poly copy = p_Copy((poly)src, (ring)lmRing, (ring)tailRing);
  if (copy == nullptr) return nullptr;
  CapturedPoly e;
  e.copy     = copy;
  e.source   = src;
  e.lmRing   = (ring)lmRing;
  e.tailRing = (ring)tailRing;
  g_captured->push_back(e);
  return (const void *)copy;
}

/* ------------------------------------------------------------------ */
/*  Mark a captured copy as final-S candidate (ENTERS poly_ptr_1)      */
/* ------------------------------------------------------------------ */
void kevlog_mark_enters_poly(const void *copy_ptr)
{
  if (!g_event_log_enabled || copy_ptr == nullptr) return;
  std::lock_guard<std::mutex> lk(g_enters_set_lock);
  if (g_enters_set == nullptr) return;
  g_enters_set->insert(copy_ptr);
}

/* ------------------------------------------------------------------ */
/*  Dump on failure                                                    */
/* ------------------------------------------------------------------ */
void kevlog_dump_on_failure(int disp_id)
{
  if (!g_event_log_enabled || g_buf == nullptr) return;

  // Dirs exist; audit-run is assumed.
  (void)mkdir("/tmp/audit-run", 0755);

  time_t now = time(NULL);
  pid_t pid = getpid();
  char base[256];
  snprintf(base, sizeof(base),
           "/tmp/audit-run/event-log-%d-%d-%ld",
           disp_id, (int)pid, (long)now);

  uint64_t nrec  = g_seq.load(std::memory_order_acquire);
  if (nrec > EVLOG_CAP_RECS) nrec = EVLOG_CAP_RECS;
  uint64_t nauxb = g_aux_off.load(std::memory_order_acquire);
  if (nauxb > EVLOG_AUX_BYTES) nauxb = EVLOG_AUX_BYTES;

  // 1) Event log binary.
  char pbin[320];
  snprintf(pbin, sizeof(pbin), "%s.bin", base);
  FILE *fb = fopen(pbin, "wb");
  if (fb != nullptr) {
    // Header: 32-byte magic/version block.
    //   offset 0..7: magic "KEVLOG01"
    //   offset 8..15: nrec (u64)
    //   offset 16..23: sizeof(KEvtRecord) = 48 (u64)
    //   offset 24..31: disp_id (u32) + pid (u32)
    uint64_t hdr[4];
    memcpy(&hdr[0], "KEVLOG01", 8);
    hdr[1] = nrec;
    hdr[2] = sizeof(KEvtRecord);
    hdr[3] = ((uint64_t)(uint32_t)disp_id) | (((uint64_t)(uint32_t)pid) << 32);
    fwrite(hdr, 1, sizeof(hdr), fb);
    fwrite(g_buf, sizeof(KEvtRecord), (size_t)nrec, fb);
    fclose(fb);
  }

  // 2) Aux arena binary.  No header; offset 0 is unused (reserved
  //    for "no aux" sentinel).  Writing the whole used range keeps
  //    aux_off values in records directly indexable.
  char paux[320];
  snprintf(paux, sizeof(paux), "%s-aux.bin", base);
  FILE *fa = fopen(paux, "wb");
  if (fa != nullptr) {
    fwrite(g_aux, 1, (size_t)nauxb, fa);
    fclose(fa);
  }

  // 3) Polys text dump.  Walk the source->copy map and render each
  //    COPY's LM — it's guaranteed-stable content because no one
  //    mutates the copy after we create it.  Keys of the map are
  //    the source pointers (historical identity), values are the
  //    copies (safe to render).  The records in *.bin reference
  //    the COPY pointer, so we emit copy addresses here for join.
  char ppoly[320];
  snprintf(ppoly, sizeof(ppoly), "%s-polys.txt", base);
  FILE *fp = fopen(ppoly, "w");
  if (fp != nullptr) {
    std::vector<CapturedPoly> entries;
    {
      std::lock_guard<std::mutex> lk(g_captured_lock);
      if (g_captured != nullptr) {
        entries.reserve(g_captured->size());
        for (auto &e : *g_captured) entries.push_back(e);
      }
    }
    fprintf(fp, "# log-owned copies, safe to render\n");
    fprintf(fp, "# FORMAT: copy_addr <TAB> lm <TAB> source_addr\n");
    fprintf(fp, "# One row PER capture (no source-to-copy dedup).  Each\n");
    fprintf(fp, "# capture is a live p_Copy snapshot taken at event time,\n");
    fprintf(fp, "# so repeated captures of the same source (e.g. ap->P.p\n");
    fprintf(fp, "# mutated by successive ksReducePoly calls) appear as\n");
    fprintf(fp, "# distinct rows with different copy_addr and potentially\n");
    fprintf(fp, "# different LMs, sharing source_addr.  Use source_addr\n");
    fprintf(fp, "# for \"same source poly\" comparisons across events.\n");
    fprintf(fp, "# %zu captures\n", entries.size());
    // Flush after every line so a SEGV partway through still leaves
    // a partial file.
    for (auto &e : entries) {
      poly copy = e.copy;
      char *lm = copy ? kt_lm_str(copy) : NULL;
      fprintf(fp, "0x%lx\t%s\t0x%lx\n",
              (unsigned long)(uintptr_t)copy,
              lm ? lm : "",
              (unsigned long)(uintptr_t)e.source);
      fflush(fp);
      if (lm) omFree(lm);
    }
    fclose(fp);
  }

  // 4) Full-polys text dump: p_String for every ENTERS-marked copy.
  //    Only candidate final-S members (poly_ptr_1 of EVT_ENTERS).
  //    These are what the Buchberger-closure checker needs to verify
  //    the final S is a correct Groebner basis.  Skipping all the
  //    mid-reduction intermediates keeps the file small.
  char pfull[320];
  snprintf(pfull, sizeof(pfull), "%s-full-polys.txt", base);
  FILE *ff = fopen(pfull, "w");
  if (ff != nullptr) {
    // Snapshot the enters set + resolve copy->ring via the poly map.
    std::vector<const void *> enters_copies;
    {
      std::lock_guard<std::mutex> lk(g_enters_set_lock);
      if (g_enters_set != nullptr) {
        enters_copies.reserve(g_enters_set->size());
        for (const void *p : *g_enters_set) enters_copies.push_back(p);
      }
    }
    // Build copy->CapturedPoly lookup from g_captured.  Each entry
    // has .copy + lmRing/tailRing — what we need to call p_String.
    // Multiple captures of the same source share nothing (each has a
    // distinct copy pointer), so last-writer-wins collisions on the
    // map key are impossible.
    std::unordered_map<const void *, CapturedPoly> copy_to_entry;
    {
      std::lock_guard<std::mutex> lk(g_captured_lock);
      if (g_captured != nullptr) {
        for (auto &e : *g_captured)
          copy_to_entry[(const void *)e.copy] = e;
      }
    }
    fprintf(ff, "# full polynomials for ENTERS arrival_id events\n");
    fprintf(ff, "# addr<TAB>full_poly_string\n");
    fprintf(ff, "# %zu enters polys marked\n", enters_copies.size());
    for (const void *cp : enters_copies) {
      auto it = copy_to_entry.find(cp);
      if (it == copy_to_entry.end()) {
        fprintf(ff, "0x%lx\t?no-entry\n", (unsigned long)(uintptr_t)cp);
        fflush(ff);
        continue;
      }
      poly p = it->second.copy;
      ring lmR = it->second.lmRing;
      ring tailR = it->second.tailRing;
      char *s = NULL;
      if (p != NULL && lmR != NULL) {
        if (tailR != NULL && tailR != lmR)
          s = p_String(p, lmR, tailR);
        else
          s = p_String(p, lmR, lmR);
      }
      fprintf(ff, "0x%lx\t%s\n",
              (unsigned long)(uintptr_t)cp,
              s ? s : "");
      fflush(ff);
      if (s) omFree(s);
    }
    fclose(ff);
  }

  fprintf(stderr, "[kevlog] DUMPED: disp=%d nrec=%lu aux=%lu -> %s.{bin,aux.bin,polys.txt,full-polys.txt}\n",
          disp_id, (unsigned long)nrec, (unsigned long)nauxb, base);
  fflush(stderr);
}
