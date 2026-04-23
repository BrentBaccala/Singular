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

/* Source-poly -> copy-poly map.  Lookup-first semantics: first sight
 * of a source allocates via p_Copy, subsequent sights return the
 * cached copy.  Deleted at shutdown. */
struct PolyCopyEntry {
  poly copy;        // fresh pointer from p_Copy (or p_Head)
  ring lmRing;      // ring used by p_Delete
  ring tailRing;    // separate tailRing, or same as lmRing
};
static std::mutex g_poly_map_lock;
static std::unordered_map<const void *, PolyCopyEntry> *g_poly_map = nullptr;

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
    std::lock_guard<std::mutex> lk(g_poly_map_lock);
    if (g_poly_map == nullptr)
      g_poly_map = new std::unordered_map<const void *, PolyCopyEntry>();
    else
      g_poly_map->clear();
  }

  fprintf(stderr, "[kevlog] init disp=%d seq_cap=%zu aux_cap=%zu\n",
          disp_id, EVLOG_CAP_RECS, EVLOG_AUX_BYTES);
  fflush(stderr);
}

/* Delete every captured poly-copy.  Called from shutdown; dump must
 * have already run (it needs the copies live to format LMs). */
static void kevlog_delete_captured_polys()
{
  std::unordered_map<const void *, PolyCopyEntry> *m = nullptr;
  {
    std::lock_guard<std::mutex> lk(g_poly_map_lock);
    m = g_poly_map;
    g_poly_map = nullptr;  // detach so concurrent captures can't use it
  }
  if (m == nullptr) return;
  for (auto &kv : *m)
  {
    poly p = kv.second.copy;
    if (p == nullptr) continue;
    if (kv.second.tailRing != nullptr && kv.second.tailRing != kv.second.lmRing)
      p_Delete(&p, kv.second.lmRing, kv.second.tailRing);
    else if (kv.second.lmRing != nullptr)
      p_Delete(&p, kv.second.lmRing);
  }
  delete m;
}

void kevlog_shutdown()
{
  if (!g_event_log_enabled) return;

  // Delete every captured poly-copy before freeing the buffer.
  kevlog_delete_captured_polys();

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
/*  Poly capture (p_Copy under mutex, source->copy identity map)       */
/* ------------------------------------------------------------------ */
const void *kevlog_capture(const void *src, struct ip_sring *r)
{
  if (!g_event_log_enabled || src == nullptr || r == nullptr) return nullptr;
  std::lock_guard<std::mutex> lk(g_poly_map_lock);
  if (g_poly_map == nullptr) return nullptr;
  auto it = g_poly_map->find(src);
  if (it != g_poly_map->end()) return (const void *)it->second.copy;
  // First sight: copy under the lock.  p_Copy bottoms out in malloc
  // (OMALLOC_USES_MALLOC=1 in this build), so it's thread-safe.  The
  // source's tail chain may be mid-mutation on another thread; this
  // is the same hazard the existing tracers face and is observed to
  // be benign for the 0.2s reproducer (see task writeup).
  poly copy = p_Copy((poly)src, (ring)r);
  if (copy == nullptr) return nullptr;
  PolyCopyEntry e;
  e.copy     = copy;
  e.lmRing   = (ring)r;
  e.tailRing = (ring)r;
  (*g_poly_map)[src] = e;
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
  std::lock_guard<std::mutex> lk(g_poly_map_lock);
  if (g_poly_map == nullptr) return nullptr;
  auto it = g_poly_map->find(src);
  if (it != g_poly_map->end()) return (const void *)it->second.copy;
  poly copy = p_Copy((poly)src, (ring)lmRing, (ring)tailRing);
  if (copy == nullptr) return nullptr;
  PolyCopyEntry e;
  e.copy     = copy;
  e.lmRing   = (ring)lmRing;
  e.tailRing = (ring)tailRing;
  (*g_poly_map)[src] = e;
  return (const void *)copy;
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
    std::vector<std::pair<const void *, PolyCopyEntry>> entries;
    {
      std::lock_guard<std::mutex> lk(g_poly_map_lock);
      if (g_poly_map != nullptr) {
        entries.reserve(g_poly_map->size());
        for (auto &kv : *g_poly_map) entries.push_back(kv);
      }
    }
    fprintf(fp, "# log-owned copies, safe to render\n");
    fprintf(fp, "# addr <TAB> lm  (addr is the COPY pointer, stored in *.bin records)\n");
    fprintf(fp, "# %zu unique poly pointers\n", entries.size());
    // Flush after every line so a SEGV partway through still leaves
    // a partial file.
    for (auto &e : entries) {
      poly copy = e.second.copy;
      char *lm = copy ? kt_lm_str(copy) : NULL;
      fprintf(fp, "0x%lx\t%s\n",
              (unsigned long)(uintptr_t)copy,
              lm ? lm : "");
      fflush(fp);
      if (lm) omFree(lm);
    }
    fclose(fp);
  }

  fprintf(stderr, "[kevlog] DUMPED: disp=%d nrec=%lu aux=%lu -> %s.{bin,aux.bin,polys.txt}\n",
          disp_id, (unsigned long)nrec, (unsigned long)nauxb, base);
  fflush(stderr);
}
