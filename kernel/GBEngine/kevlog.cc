/**
 * @file kevlog.cc
 * @brief Global event log implementation (task 325).
 *
 * See kevlog.h for the record schema + semantics.  Summary:
 *   - Preallocate 1 GiB buffer at bba entry (when gated).
 *   - Preallocate 256 MiB aux arena.
 *   - Global atomic seq; each event reserves its slot via fetch_add.
 *   - Overflow => abort.
 *   - Dump to disk only on SINGULAR_CHECK_IDEAL_MEMBERSHIP violation.
 *   - "Leak-on-trace" keeps poly pointers live: the kt_pLmFree /
 *     kt_pDelete / kt_p_LmFree / kt_p_Delete wrappers (see
 *     kevlog_wrap.h) are no-ops when g_defer_frees is true.  Process
 *     exit reclaims memory — acceptable for the 0.2s reproducer.
 */

#include "kernel/GBEngine/kevlog.h"
#include "kernel/GBEngine/kthread.h"  // kt_lm_str (used at dump time only)
#include "kernel/polys.h"
#include "polys/monomials/p_polys.h"
#include "omalloc/omalloc.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <atomic>
#include <mutex>
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
bool g_defer_frees = false;

static constexpr size_t EVLOG_CAP_BYTES  = (size_t)1 << 30;   // 1 GiB
static constexpr size_t EVLOG_CAP_RECS   = EVLOG_CAP_BYTES / sizeof(KEvtRecord);
static constexpr size_t EVLOG_AUX_BYTES  = (size_t)256 << 20; // 256 MiB

static KEvtRecord       *g_buf = nullptr;
static unsigned char    *g_aux = nullptr;
static std::atomic<uint64_t> g_seq{0};
static std::atomic<uint64_t> g_aux_off{1};  // reserve offset 0 for "no aux"
static int                g_dispatch_id = -1;

/* Registry of distinct poly pointers for the polys.txt dump. */
static std::mutex               g_poly_reg_lock;
static std::unordered_set<const void *> *g_poly_reg = nullptr;

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
    std::lock_guard<std::mutex> lk(g_poly_reg_lock);
    if (g_poly_reg == nullptr)
      g_poly_reg = new std::unordered_set<const void *>();
    else
      g_poly_reg->clear();
  }

  // "Leak-on-trace" mode: the free wrappers become no-ops so all
  // poly pointers captured in events survive until dump.
  g_defer_frees = true;

  fprintf(stderr, "[kevlog] init disp=%d seq_cap=%zu aux_cap=%zu\n",
          disp_id, EVLOG_CAP_RECS, EVLOG_AUX_BYTES);
  fflush(stderr);
}

void kevlog_shutdown()
{
  if (!g_event_log_enabled) return;

  // Leak-on-trace: we deliberately do NOT reclaim the poly memory
  // that the free wrappers skipped.  Process exit handles it.
  g_defer_frees = false;

  if (g_buf) { free(g_buf); g_buf = nullptr; }
  if (g_aux) { free(g_aux); g_aux = nullptr; }

  {
    std::lock_guard<std::mutex> lk(g_poly_reg_lock);
    if (g_poly_reg != nullptr) { delete g_poly_reg; g_poly_reg = nullptr; }
  }

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

  // Register the pointers for the polys.txt dump.
  if (poly_ptr_1 != nullptr) kevlog_register_poly(poly_ptr_1);
  if (poly_ptr_2 != nullptr) kevlog_register_poly(poly_ptr_2);
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

void kevlog_register_poly(const void *p)
{
  if (!g_event_log_enabled || p == nullptr) return;
  std::lock_guard<std::mutex> lk(g_poly_reg_lock);
  if (g_poly_reg != nullptr) g_poly_reg->insert(p);
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

  // 3) Polys text dump.  For every pointer we registered, try
  //    p_String.  Because defer-frees is ON throughout bba, these
  //    pointers are still live.
  char ppoly[320];
  snprintf(ppoly, sizeof(ppoly), "%s-polys.txt", base);
  FILE *fp = fopen(ppoly, "w");
  if (fp != nullptr) {
    std::vector<const void *> ptrs;
    {
      std::lock_guard<std::mutex> lk(g_poly_reg_lock);
      if (g_poly_reg != nullptr)
        for (const void *p : *g_poly_reg) ptrs.push_back(p);
    }
    fprintf(fp, "# addr <TAB> lm\n");
    fprintf(fp, "# %zu unique poly pointers\n", ptrs.size());
    // We only write the LM (via kt_lm_str which uses direct exponent
    // iteration, no ring-bin allocations) — it's safe even if the
    // pointer's tail chain is gone.  p_String would walk pNext and
    // risk a SEGV on partially-freed polys, so we avoid it.
    // Flush after every line so a SEGV partway through still leaves
    // a partial file.
    for (const void *p : ptrs) {
      char *lm = p ? kt_lm_str((poly)p) : NULL;
      fprintf(fp, "0x%lx\t%s\n",
              (unsigned long)(uintptr_t)p,
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
