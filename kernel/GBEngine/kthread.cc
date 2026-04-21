/**
 * @file kthread.cc
 * @brief Parallel Groebner basis reduction — continuous refill design
 *        (milestone (d) of the continuous-cursor redesign, task 283).
 *
 * All threads cooperatively sweep T[0..tl] for active polynomials, but
 * unlike the earlier barrier-based design there is no round structure.
 * Main continuously drains survivors + refills empty slots; workers
 * continuously pull tiles off a single atomic cursor and sweep them.
 *
 *   main:
 *     loop:
 *       drain survivor FIFO (main-only)
 *       for each EMPTY slot: try to refill from L, publish tile batch
 *       for each slot marked needs_republish: publish a fresh batch
 *       if L empty AND survivor queue empty AND all slots empty: break
 *       else wait on slot_freed_cv (short)
 *
 *   worker:
 *     loop:
 *       claim tile id by CAS on (tile_cursor, tile_end)
 *       if no tile available: wait on tiles_avail_cv (with shutdown check)
 *       decode (slot, slice) from tile id via the batch ring
 *       sweep slot's slice of T[0..sl_snapshot]
 *       if fetch_sub on tiles_remaining == 1: close the slot
 *                                              (merge + reduce-one-step)
 *       if slot transitions to EMPTY or produces a survivor:
 *         signal slot_freed_cv so main can advance
 *
 * Thread safety relies on:
 *   - posInT0 (append at end, no memmove)
 *   - T is read-only during sweep+reduce; modified only during
 *     process_survivor (after all slots are done)
 *   - --disable-omalloc build (thread-safe malloc/free)
 *   - currRing and si_opt_1/si_opt_2 are __thread
 *   - OPT_REDTHROUGH forced (no lazy pushback to L from reduction)
 *   - tailRing pre-expanded before parallel phase
 *   - Per-thread SweepResults: no contention during sweep
 *   - ksCreateSpoly race fix (task 267) preserved
 */

#include "kernel/GBEngine/kthread.h"
#include "kernel/GBEngine/kstd1.h"
#include "kernel/GBEngine/kutil.h"
#include "kernel/polys.h"
#include "kernel/ideals.h"
#include "polys/monomials/p_polys.h"
#include "misc/options.h"

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cstring>

/* ------------------------------------------------------------------ */
/*  Instrumentation helpers (task 482)                                 */
/* ------------------------------------------------------------------ */
#ifdef KTHREAD_INSTRUMENT
#  define KT_STATS(ctx)       ((ctx)->stats_enabled)
#  define KT_TS(ctx, tid)     ((ctx)->tstats[tid])
#  define KT_TIME_START(var)  long var = kt_now_ns()
#  define KT_TIME_DELTA(var)  (kt_now_ns() - (var))

static inline void kt_L_lock(SweepContext *ctx, int thread_id)
{
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    pthread_mutex_lock(&ctx->L_lock);
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.L_lock_wait_ns += dt;
    ts.L_lock_count++;
  }
  else
  {
    pthread_mutex_lock(&ctx->L_lock);
  }
}

static inline void kt_surv_q_lock(SweepContext *ctx, int thread_id)
{
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    pthread_mutex_lock(&ctx->survivor_queue_mutex);
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.surv_q_wait_ns += dt;
    ts.surv_q_count++;
  }
  else
  {
    pthread_mutex_lock(&ctx->survivor_queue_mutex);
  }
}

static inline void kt_record_reduce(SweepContext *ctx, int thread_id, long start_ns, long duration_ns)
{
  if (!KT_STATS(ctx)) return;
  pthread_mutex_lock(&ctx->stats_lock);
  ReduceEvent ev = { thread_id, start_ns, duration_ns };
  ctx->reduces_this_round->push_back(ev);
  pthread_mutex_unlock(&ctx->stats_lock);
}

/* Phase-0 S exclusive lock: short critical section for enterS(h). */
static inline void kt_S_lock_exclusive(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    strat->S.lock_exclusive();
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.phase0_wait_ns += dt;
  }
  else
  {
    strat->S.lock_exclusive();
  }
}

/* Phase-1 S shared lock: concurrent iteration while no writers run. */
static inline void kt_S_lock_shared(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    strat->S.lock_shared();
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.phase1_wait_ns += dt;
  }
  else
  {
    strat->S.lock_shared();
  }
}

/* L exclusive lock held during phase 1 (chainCritNormal merges local B
 * into strat->L).  Records wait time in phase1_l_wait_ns (renamed from
 * phase2_wait_ns in task 508 — phase 2 no longer exists). */
static inline void kt_L_lock_phase2(SweepContext *ctx, int thread_id)
{
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    pthread_mutex_lock(&ctx->L_lock);
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.phase1_l_wait_ns += dt;
  }
  else
  {
    pthread_mutex_lock(&ctx->L_lock);
  }
}
#else
#  define KT_STATS(ctx)       (false)
#  define KT_TIME_START(var)  ((void)0)
#  define KT_TIME_DELTA(var)  (0L)
static inline void kt_L_lock(SweepContext *ctx, int /*tid*/) { pthread_mutex_lock(&ctx->L_lock); }
static inline void kt_surv_q_lock(SweepContext *ctx, int /*tid*/) { pthread_mutex_lock(&ctx->survivor_queue_mutex); }
static inline void kt_record_reduce(SweepContext*, int, long, long) {}
static inline void kt_S_lock_exclusive(SweepContext *ctx, int /*tid*/) { ctx->strat->S.lock_exclusive(); }
static inline void kt_S_lock_shared(SweepContext *ctx, int /*tid*/) { ctx->strat->S.lock_shared(); }
static inline void kt_L_lock_phase2(SweepContext *ctx, int /*tid*/) { pthread_mutex_lock(&ctx->L_lock); }
#endif

/* ------------------------------------------------------------------ */
/*  posInT override: always append at end during parallel mode         */
/* ------------------------------------------------------------------ */

static int posInT_appendEnd(const BlockArray<TObject> & /*T*/, const int tl, LObject & /*h*/)
{
  return tl + 1;
}

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

int get_singular_threads()
{
  const char *env = getenv("SINGULAR_THREADS");
  if (env == NULL) return 1;
  int n = atoi(env);
  if (n < 1) return 1;
  if (n > 64) return 64;
  return n;
}

static inline SweepResult& sweep_result(SweepContext *ctx, int thread_id, int slot)
{
  return ctx->sweep_results[thread_id * ctx->max_active + slot];
}

/* ------------------------------------------------------------------ */
/*  Task 512: T-integrity audit (debug).                               */
/*                                                                     */
/*  Walks strat->T and checks that each entry's stored pLength matches */
/*  the actual chain length.  First mismatch aborts with a diagnostic. */
/*  Expensive — O(|T| * avg chain length).  Controlled by env var      */
/*  SINGULAR_AUDIT_T.  Run from chokepoints (drain entry, refill entry,*/
/*  reduce entry, redtailBba post) to narrow the window where          */
/*  corruption first appears.                                          */
/*                                                                     */
/*  Additions for the stale-pLength investigation:                     */
/*   1. Per-thread "current operation" tag (g_thread_debug), set by    */
/*      an RAII OpTag helper at key entry points.  On audit-fire we    */
/*      dump every thread's tag so we can see what the OTHER threads   */
/*      were doing at the moment a T entry was observed truncated.     */
/*   2. Ring buffer of pLength observations (g_plen_ring).  Every      */
/*      audit call records (tsc, tid, Tidx, head, stored, actual) so   */
/*      we can find the same head observed with different lengths on   */
/*      different threads — direct evidence of transient truncation.   */
/* ------------------------------------------------------------------ */
// SINGULAR_DEBUG_RING=1  — event ring, OpTag breadcrumbs, dErrorBreak
//                          hook active.  Low overhead; safe to leave on.
// SINGULAR_AUDIT_T=1     — additionally, audit_T_pLength walks every
//                          published T entry at each chokepoint.
//                          Expensive; floods the ring with PLEN events.
// Using the ring WITHOUT the walks lets TAG breadcrumbs (enterT,
// ksReducePoly sub-steps, etc.) survive long enough to reach a
// dErrorBreak hook dump without getting drowned out.
static std::atomic<bool> g_debug_ring_enabled{false};
static std::atomic<bool> g_audit_T_enabled{false};
static std::atomic<int>  g_audit_T_checked{0};
static FILE *g_audit_log = NULL;  // opened when ring/audit is enabled

#define MAX_AUDIT_THREADS 64

struct ThreadDebugState
{
  const char *op;       // tag string (static literal); NULL = never tagged
  void       *arg_poly; // poly pointer relevant to the op (LObject p, T[i].p, …)
  int         arg_slot; // active-slot or T-index
  int         arg_int;  // generic extra (e.g. i_r)
  uint64_t    tsc;      // timestamp of last tag (rdtsc)
};

static ThreadDebugState g_thread_debug[MAX_AUDIT_THREADS];

static inline void tag_log(int tid, const char *op, void *poly,
                           int slot, int arg);  // forward decl
static void dump_ring_generic(const char *reason);  // forward decl

// Thread-local tid used by the kt_debug_tag() breadcrumb API (kthread.h)
// so callers in other translation units (e.g. kspoly.cc) don't have to
// pass tid through their signatures.  Default 0 means main thread in
// serial code paths — harmless: the audit is disabled there.
thread_local int kt_debug_tid = 0;

// Public breadcrumb entry point.  Thin wrapper around tag_log.
void kt_debug_tag(const char *op, void *poly, int slot, int arg)
{
  tag_log(kt_debug_tid, op, poly, slot, arg);
}

// Dedicated audit-log line for the enterT "slot-in-use" probe:
// called unconditionally from enterT; if atT points at a slot whose
// .p is already non-NULL, we've either an OVERWRITE (breaks append-
// only) or a SHIFT (moves existing entries up).  Both invalidate
// pair.i_r2 values that reference affected slots.  Emits directly
// to g_audit_log so the line survives ring-buffer wrap.  Silent
// when the debug ring isn't enabled.
void kt_debug_enterT_slot_probe(int atT, int pre_T_size,
                                void *pre_T_atT_p, void *new_p)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
  bool overwrite = (pre_T_atT_p != NULL && pre_T_atT_p != new_p);
  bool shift = (atT < pre_T_size);
  if (!overwrite && !shift) return;   // pure append — not interesting
  FILE *log = g_audit_log ? g_audit_log : stderr;
  fprintf(log,
          "=== enterT:SLOT_IN_USE atT=%d  pre-T.size=%d  "
          "pre-T[atT].p=%p  new p=%p  overwrite=%d shift=%d ===\n",
          atT, pre_T_size, pre_T_atT_p, new_p,
          overwrite ? 1 : 0, shift ? 1 : 0);
  fflush(log);
  tag_log(kt_debug_tid,
          overwrite ? "enterT:OVERWRITE_SLOT" : "enterT:SHIFT_ONLY",
          pre_T_atT_p, atT, pre_T_size);
}

// General-purpose audit-log emitter.  Writes a formatted line to
// g_audit_log (or stderr) iff the debug ring is enabled.  Use for
// diagnostics that need to survive ring-buffer wrap.
void kt_debug_audit_printf(const char *fmt, ...)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
  FILE *log = g_audit_log ? g_audit_log : stderr;
  va_list ap;
  va_start(ap, fmt);
  vfprintf(log, fmt, ap);
  va_end(ap);
  fflush(log);
}

// Scan strat->L and count pairs whose R[i_r2]->p doesn't match
// pair.p2 (and same for p1/i_r1).  Emits direct to audit log on
// any nonzero count.  Call at parallel-phase checkpoints to
// time-bisect when R-mismatch first appears.  `label` is a short
// site name; included in the audit line.
//
// WARNING: walks all of strat->L without taking L_lock.  Cheap, but
// may observe a torn pair if a concurrent insert is mid-flight.
// Acceptable for diagnostic narrowing — false-positive rate from
// torn reads is negligible compared to the bug rate.
int kt_debug_R_scan_L(void *strat_void, const char *label)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return 0;
  kStrategy strat = (kStrategy)strat_void;
  int inconsistent_R = 0;
  int total = 0;
  int first_bad_k = -1;
  void *first_bad_pair_p = NULL;
  void *first_bad_R_p = NULL;
  for (auto it = strat->L.begin(); it != strat->L.end(); ++it)
  {
    total++;
    if (it->i_r2 >= 0 && it->i_r2 < (int)strat->T.size())
    {
      TObject *r_entry = strat->R[it->i_r2];
      poly r_p = (r_entry != NULL) ? r_entry->p : NULL;
      if (it->p2 != NULL && r_p != NULL && it->p2 != r_p)
      {
        inconsistent_R++;
        if (first_bad_k < 0)
        {
          first_bad_k = it->i_r2;
          first_bad_pair_p = (void*)it->p2;
          first_bad_R_p = (void*)r_p;
        }
      }
    }
    if (it->i_r1 >= 0 && it->i_r1 < (int)strat->T.size())
    {
      TObject *r_entry = strat->R[it->i_r1];
      poly r_p = (r_entry != NULL) ? r_entry->p : NULL;
      if (it->p1 != NULL && r_p != NULL && it->p1 != r_p)
      {
        inconsistent_R++;
        if (first_bad_k < 0)
        {
          first_bad_k = it->i_r1;
          first_bad_pair_p = (void*)it->p1;
          first_bad_R_p = (void*)r_p;
        }
      }
    }
  }
  if (inconsistent_R > 0)
  {
    FILE *log = g_audit_log ? g_audit_log : stderr;
    fprintf(log,
            "=== R_SCAN_L (%s): |L|=%d  inconsistent_R=%d  "
            "first_bad: k=%d pair_p=%p R[k]->p=%p ===\n",
            label, total, inconsistent_R, first_bad_k,
            first_bad_pair_p, first_bad_R_p);
    fflush(log);
  }
  return inconsistent_R;
}

// R-write probe.  Called BEFORE each write to strat->R[k].  Compares
// the pre-write R[k] (and its ->p) to the target being written.
// Emits a direct audit-log line when the rewrite is "meaningful":
//   - old R[k] was non-NULL, and
//   - old R[k]->p is different from new target's ->p.
// A shift-loop rewrite that keeps R[k]->p the same (just points at
// a new T-address) is expected and stays silent.  The "fresh-slot
// with stale contents" rewrite — the suspected bug — fires loud.
//
// `site` is a short string naming the call site (free-form).
// `old_target` is the pre-write R[k] TObject* (may be NULL).
// `new_target` is the T-addr about to be written (should be non-NULL).
// `old_p` / `new_p` are the .p fields of the TObjects, captured by
//   the caller (passing them in avoids re-dereferencing possibly-
//   stale pointers here).
void kt_debug_R_write_probe(const char *site, int k,
                            void *old_target, void *new_target,
                            void *old_p, void *new_p)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
  if (old_target == NULL) return;              // first write, benign
  if (old_p == new_p && old_target == new_target) return;  // no-op
  if (old_p == new_p) {
    // R[k] moved to a new T-address but still points at same .p —
    // this is the shift-loop's benign update.  Silent.
    return;
  }
  // Meaningful rewrite — old and new .p differ.
  FILE *log = g_audit_log ? g_audit_log : stderr;
  fprintf(log,
          "=== R_WRITE:OVERWRITE site=%s  R[%d]  "
          "old_target=%p old_p=%p  new_target=%p new_p=%p ===\n",
          site, k, old_target, new_target, old_p, new_p);
  fflush(log);
  tag_log(kt_debug_tid, "R:OVERWRITE_DIFF_POLY",
          old_p, k, 0);
}

// ---------------------------------------------------------------------
// T-node registry: records every chain node address that has been
// enterT'd into a T entry, along with the T index it belongs to.
// A mutation of any of these addresses' pNext is a bug — T entry
// chains are supposed to be immutable post-enterT.
//
// Implementation: open-addressed hash table with atomic CAS insert
// and atomic-acquire lookup.  Sized generously; overflow is silent
// (we'd see missing registrations in the ring, not crashes).
// ---------------------------------------------------------------------
#define TNODE_REGISTRY_SIZE (256 * 1024)   // 2 MB, plenty for ~100k T nodes
struct TNodeSlot { std::atomic<void*> addr; int tidx; };
static TNodeSlot g_tnode_registry[TNODE_REGISTRY_SIZE];

static inline uint32_t tnode_hash(void *addr)
{
  uint64_t h = (uint64_t)(uintptr_t)addr;
  h ^= h >> 33;
  h *= 0xff51afd7ed558ccdULL;
  h ^= h >> 33;
  return (uint32_t)h & (TNODE_REGISTRY_SIZE - 1);
}

static void tnode_register(void *addr, int tidx)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
  if (addr == NULL) return;
  uint32_t idx = tnode_hash(addr);
  for (int probe = 0; probe < 32; probe++)
  {
    void *expected = NULL;
    if (g_tnode_registry[idx].addr.compare_exchange_strong(
          expected, addr, std::memory_order_acq_rel))
    {
      g_tnode_registry[idx].tidx = tidx;
      return;
    }
    if (expected == addr) return;  // already registered
    idx = (idx + 1) & (TNODE_REGISTRY_SIZE - 1);
  }
  // Table is packed — silently drop; we'll lose this node but not crash.
}

// Returns the T index if addr is a registered T-entry chain node,
// else -1.  Lock-free read.
static int tnode_lookup(void *addr)
{
  if (addr == NULL) return -1;
  uint32_t idx = tnode_hash(addr);
  for (int probe = 0; probe < 32; probe++)
  {
    void *v = g_tnode_registry[idx].addr.load(std::memory_order_acquire);
    if (v == addr) return g_tnode_registry[idx].tidx;
    if (v == NULL) return -1;
    idx = (idx + 1) & (TNODE_REGISTRY_SIZE - 1);
  }
  return -1;
}

// Public hook called from hot-path pNext-write sites (kbuckets.cc etc.)
// to detect a mutation of a registered T chain node.  Fires a TAG
// event and doesn't abort — the subsequent STALE_l2 probe will fire
// normally and we can grep the ring for the `TNODE:mutation` tags.
void kt_debug_check_tnode_write(const char *site, void *addr)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
  int tidx = tnode_lookup(addr);
  if (tidx < 0) return;
  tag_log(kt_debug_tid, site, addr, tidx, /*arg=*/-1);
}

// Snapshot of each T entry's head pointer at enterT time.  Lets us
// detect whether T[atT].p gets overwritten post-enterT (which it
// shouldn't, per the immutability invariant).
#define TNODE_HEAD_SNAPSHOT_SIZE 16384
static std::atomic<void*> g_tnode_head_snapshot[TNODE_HEAD_SNAPSHOT_SIZE];

void kt_debug_snapshot_T_head(int tidx, void *addr)
{
  if (tidx < 0 || tidx >= TNODE_HEAD_SNAPSHOT_SIZE) return;
  g_tnode_head_snapshot[tidx].store(addr, std::memory_order_release);
}

void *kt_debug_lookup_T_head(int tidx)
{
  if (tidx < 0 || tidx >= TNODE_HEAD_SNAPSHOT_SIZE) return NULL;
  return g_tnode_head_snapshot[tidx].load(std::memory_order_acquire);
}

// Register a chain node in the T-node set.  Called from enterT for
// every node of a freshly-published T entry's chain.  Before the
// register, check if the node is already registered to a DIFFERENT
// T entry — that's chain aliasing at the node-address level, which
// is the bug we're hunting.  Fire a TAG event with (new_atT, prior
// tidx, address) and dump via the existing dErrorBreak hook.
void kt_debug_register_tnode(void *addr, int tidx)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) {
    tnode_register(addr, tidx);
    return;
  }
  int prior = tnode_lookup(addr);
  if (prior >= 0 && prior != tidx)
  {
    // Aliasing detected.  Encode (new_atT, prior) in the TAG args.
    tag_log(kt_debug_tid, "enterT:ALIASING", addr, tidx, prior);
    // Dump once and abort if the user wants an immediate stop.
    static std::atomic<int> once{0};
    if (once.fetch_add(1) == 0)
    {
      dump_ring_generic("enterT:ALIASING");
      if (getenv("SINGULAR_ABORT_ON_DERROR") != NULL) abort();
    }
  }
  tnode_register(addr, tidx);
}

// RAII tag helper: saves the previous tag on construction, restores on
// destruction, so nested tags work correctly without the tag silently
// becoming stale when the inner scope returns.
struct OpTag
{
  int tid;
  ThreadDebugState saved;

  OpTag(int t, const char *op, void *poly = NULL, int slot = -1, int arg = 0)
    : tid(t)
  {
    if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
    if (tid < 0 || tid >= MAX_AUDIT_THREADS) return;
    ThreadDebugState &d = g_thread_debug[tid];
    saved = d;
    d.op       = op;
    d.arg_poly = poly;
    d.arg_slot = slot;
    d.arg_int  = arg;
    d.tsc      = (uint64_t)__builtin_ia32_rdtsc();
    tag_log(tid, op, poly, slot, arg);
  }
  ~OpTag()
  {
    if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
    if (tid < 0 || tid >= MAX_AUDIT_THREADS) return;
    g_thread_debug[tid] = saved;
    // Log the restoration so the ring shows "tid left scope X, returned to Y".
    tag_log(tid, saved.op ? saved.op : "<none>",
            saved.arg_poly, saved.arg_slot, saved.arg_int);
  }
};

// Forward decls for tag_log (defined after OpTag uses it).
// (Already defined above; OpTag is declared after the event machinery.)

// Unified ring buffer of events.
// Two event kinds share one ring (so ordering is preserved across both):
//   TAG   — a thread entered an OpTag scope (op/poly/slot/arg valid, stored=-1)
//   PLEN  — audit observed pLength (Tidx/head/stored/actual valid, op=NULL)
// On dump we walk the ring backwards to reconstruct, for each thread, its
// most recent TAG event at or before the mismatch.
enum EventKind { EV_TAG = 0, EV_PLEN = 1 };

struct Event
{
  std::atomic<uint64_t> seq;  // written last on publish (acts as valid-marker)
  uint64_t    tsc;
  int         tid;
  int         kind;       // EventKind
  const char *op;         // TAG: function name; PLEN: NULL
  void       *poly;       // TAG: arg_poly; PLEN: chain head
  int         slot;       // TAG: arg_slot; PLEN: Tidx
  int         arg;        // TAG: arg_int; PLEN: stored
  int         actual;     // TAG: unused; PLEN: observed chain length
};

#define PLEN_RING_SIZE 1048576
static Event g_plen_ring[PLEN_RING_SIZE];
static std::atomic<uint64_t> g_plen_seq{0};
// Set to true at the moment of first mismatch so other threads stop
// logging and the dump sees a stable ring.  Threads continue running
// (they just don't record events), but the ring is frozen.
static std::atomic<bool> g_plen_frozen{false};

static inline void plen_log(int tid, int Tidx, void *head,
                            int stored, int actual)
{
  if (g_plen_frozen.load(std::memory_order_relaxed)) return;
  uint64_t s = g_plen_seq.fetch_add(1, std::memory_order_relaxed);
  Event &e = g_plen_ring[s & (PLEN_RING_SIZE - 1)];
  e.tsc    = (uint64_t)__builtin_ia32_rdtsc();
  e.tid    = tid;
  e.kind   = EV_PLEN;
  e.op     = NULL;
  e.poly   = head;
  e.slot   = Tidx;
  e.arg    = stored;
  e.actual = actual;
  e.seq.store(s, std::memory_order_release);
}

static inline void tag_log(int tid, const char *op, void *poly,
                           int slot, int arg)
{
  if (g_plen_frozen.load(std::memory_order_relaxed)) return;
  uint64_t s = g_plen_seq.fetch_add(1, std::memory_order_relaxed);
  Event &e = g_plen_ring[s & (PLEN_RING_SIZE - 1)];
  e.tsc    = (uint64_t)__builtin_ia32_rdtsc();
  e.tid    = tid;
  e.kind   = EV_TAG;
  e.op     = op;
  e.poly   = poly;
  e.slot   = slot;
  e.arg    = arg;
  e.actual = 0;
  e.seq.store(s, std::memory_order_release);
}

// Dump per-thread tag history, pLength observations for the matching head,
// and recent mismatches.  Called once, on the first audit mismatch, just
// before abort().
static void dump_debug_state(const char *tag, int my_tid, uint64_t mismatch_seq,
                             int match_Tidx, void *match_head)
{
  fprintf(g_audit_log,
          "\n=== audit fired: %s tid=%d T[%d] head=%p seq=%lu ===\n",
          tag, my_tid, match_Tidx, match_head, (unsigned long)mismatch_seq);

  uint64_t cur = g_plen_seq.load(std::memory_order_acquire);
  uint64_t start = (cur > PLEN_RING_SIZE) ? cur - PLEN_RING_SIZE : 0;

  // Per-thread tag history: for each thread, find its most recent TAG
  // event at or before the mismatch seq (i.e. what it was doing AT the
  // moment of truncation, not what it moved to afterward).
  fprintf(g_audit_log,
          "\n-- per-thread tag at mismatch (most recent TAG "
          "at seq <= %lu) --\n", (unsigned long)mismatch_seq);
  for (int t = 0; t < MAX_AUDIT_THREADS; t++)
  {
    // Walk backwards from mismatch_seq to find this thread's latest tag.
    bool found = false;
    uint64_t scan_lo = (mismatch_seq > PLEN_RING_SIZE)
                        ? mismatch_seq - PLEN_RING_SIZE : 0;
    for (uint64_t idx = mismatch_seq; idx-- > scan_lo; )
    {
      Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
      uint64_t ee = e.seq.load(std::memory_order_acquire);
      if (ee != idx) continue;
      if (e.tid != t) continue;
      if (e.kind != EV_TAG) continue;
      fprintf(g_audit_log,
              "  tid=%d seq=%lu tsc=%lu op=%s poly=%p slot=%d arg=%d\n",
              t, (unsigned long)ee, (unsigned long)e.tsc,
              e.op ? e.op : "<none>", e.poly, e.slot, e.arg);
      found = true;
      break;
    }
    if (!found)
    {
      // Also check if this thread has ANY recent event at all (to
      // distinguish "never seen" from "tag evicted from ring").
      for (uint64_t idx = cur; idx-- > start; )
      {
        Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
        if (e.seq.load(std::memory_order_acquire) != idx) continue;
        if (e.tid != t) continue;
        // Thread has events, but no TAG at or before mismatch_seq in ring.
        fprintf(g_audit_log, "  tid=%d <no tag in window> (has later events)\n", t);
        found = true;
        break;
      }
    }
  }

  // All PrepareRed and refill events in the ring (regardless of poly).
  // Lets us spot whether the match_head or an aliasing poly was ever
  // PrepareRed'd during the run.
  fprintf(g_audit_log,
          "\n-- all PrepareRed / refill events (entire ring) --\n");
  int pr_events = 0;
  for (uint64_t idx = start; idx < cur; idx++)
  {
    Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
    uint64_t ee = e.seq.load(std::memory_order_acquire);
    if (ee != idx) continue;
    if (e.kind != EV_TAG) continue;
    if (e.op == NULL) continue;
    if (strncmp(e.op, "PrepareRed", 10) != 0 &&
        strncmp(e.op, "refill_and_publish", 18) != 0)
      continue;
    const char *match = (e.poly == match_head) ? " <-- MATCH" : "";
    fprintf(g_audit_log,
            "  seq=%lu tsc=%lu tid=%d op=%s poly=%p slot=%d arg=%d%s\n",
            (unsigned long)ee, (unsigned long)e.tsc,
            e.tid, e.op, e.poly, e.slot, e.arg, match);
    pr_events++;
    if (pr_events > 400) { fprintf(g_audit_log, "  ... (truncated after 400)\n"); break; }
  }
  if (pr_events == 0)
    fprintf(g_audit_log, "  (no PrepareRed/refill events in ring)\n");

  // All events (TAG + PLEN) in the ring whose poly matches the mismatch
  // head.  This catches the PrepareRed / ksReducePoly / other ops that
  // took this poly as an argument — even ones before the neighborhood.
  fprintf(g_audit_log,
          "\n-- all events for poly=%p (entire ring) --\n", match_head);
  int poly_events = 0;
  for (uint64_t idx = start; idx < cur; idx++)
  {
    Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
    uint64_t ee = e.seq.load(std::memory_order_acquire);
    if (ee != idx) continue;
    if (e.poly != match_head) continue;
    if (e.kind == EV_TAG)
      fprintf(g_audit_log,
              "  seq=%lu tsc=%lu tid=%d TAG op=%s slot=%d arg=%d\n",
              (unsigned long)ee, (unsigned long)e.tsc,
              e.tid, e.op ? e.op : "<none>", e.slot, e.arg);
    else
      fprintf(g_audit_log,
              "  %c seq=%lu tsc=%lu tid=%d PLEN T[%d] stored=%d actual=%d\n",
              (e.arg != e.actual) ? '*' : ' ',
              (unsigned long)ee, (unsigned long)e.tsc,
              e.tid, e.slot, e.arg, e.actual);
    poly_events++;
    if (poly_events > 200) { fprintf(g_audit_log, "  ... (truncated after 200)\n"); break; }
  }
  if (poly_events == 0)
    fprintf(g_audit_log, "  (no events with this poly in ring)\n");

  fprintf(g_audit_log,
          "\n-- pLength observations for T[%d] "
          "(stored != actual marked *) --\n", match_Tidx);
  int printed = 0;
  for (uint64_t idx = start; idx < cur; idx++)
  {
    Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
    uint64_t ee = e.seq.load(std::memory_order_acquire);
    if (ee != idx) continue;
    if (e.kind != EV_PLEN) continue;
    if (e.slot != match_Tidx) continue;
    fprintf(g_audit_log,
            "  %c seq=%lu tsc=%lu tid=%d T[%d] head=%p stored=%d actual=%d\n",
            (e.arg != e.actual) ? '*' : ' ',
            (unsigned long)ee, (unsigned long)e.tsc,
            e.tid, e.slot, e.poly, e.arg, e.actual);
    printed++;
    if (printed > 60) { // keep output bounded
      fprintf(g_audit_log, "  ... (truncated after 60)\n");
      break;
    }
  }
  if (printed == 0)
    fprintf(g_audit_log, "  (no observations of T[%d] in ring)\n", match_Tidx);

  // Neighborhood: all events from all threads within +/- 64 seq of the
  // mismatch. Re-load cur so we include events logged during the earlier
  // parts of this dump (threads are still running).
  uint64_t cur2 = g_plen_seq.load(std::memory_order_acquire);
  uint64_t nlo = (mismatch_seq > 128) ? mismatch_seq - 128 : 0;
  uint64_t nhi = mismatch_seq + 64;
  if (nhi > cur2) nhi = cur2;
  fprintf(g_audit_log,
          "\n-- event neighborhood (seq %lu..%lu, cur=%lu) --\n",
          (unsigned long)nlo, (unsigned long)nhi, (unsigned long)cur2);
  int emitted = 0;
  for (uint64_t idx = nlo; idx < nhi; idx++)
  {
    Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
    uint64_t ee = e.seq.load(std::memory_order_acquire);
    if (ee != idx) continue;
    const char *marker = (idx == mismatch_seq) ? ">>" : "  ";
    if (e.kind == EV_TAG)
      fprintf(g_audit_log,
              "%s seq=%lu tsc=%lu tid=%d TAG op=%s poly=%p slot=%d arg=%d\n",
              marker, (unsigned long)ee, (unsigned long)e.tsc,
              e.tid, e.op ? e.op : "<none>", e.poly, e.slot, e.arg);
    else
      fprintf(g_audit_log,
              "%s seq=%lu tsc=%lu tid=%d PLEN T[%d] head=%p stored=%d actual=%d%s\n",
              marker, (unsigned long)ee, (unsigned long)e.tsc,
              e.tid, e.slot, e.poly, e.arg, e.actual,
              (e.arg != e.actual) ? " *" : "");
    emitted++;
  }
  fprintf(g_audit_log,
          "-- end of dump: %d events emitted --\n", emitted);
  fflush(g_audit_log);
  // Don't fclose: other threads may still call fprintf(g_audit_log,...)
  // via audit_T_pLength's pre-dump error path (subsequent mismatches
  // race with our abort() call).  fflush above already forces the write
  // to disk; fclose would just invite NULL derefs.
}

// Dump a generic ring-buffer snapshot (no specific T[j] to filter on).
// Called from the dErrorBreak hook so bucket-length dErrors and other
// fault sites get the same treatment the T-entry audit got.
static void dump_ring_generic(const char *reason)
{
  if (!g_debug_ring_enabled.load(std::memory_order_relaxed)) return;
  if (g_audit_log == NULL) return;
  // Freeze to keep the ring stable while we dump.  Note: this is a
  // one-shot — subsequent dErrors land in a frozen ring (that's fine,
  // the first one is the one we care about).
  bool was_frozen = g_plen_frozen.exchange(true, std::memory_order_acq_rel);
  if (was_frozen) return;  // someone else got here first

  uint64_t cur = g_plen_seq.load(std::memory_order_acquire);
  fprintf(g_audit_log,
          "\n=== dErrorBreak hook: %s cur=%lu ===\n",
          reason, (unsigned long)cur);

  // Per-thread tag at the moment of fault: walk the ring backwards
  // from `cur` to find each thread's most recent TAG event.  This
  // is the technique that pinned the _p_LmTest race — knowing what
  // OTHER threads were doing at the fault site is usually more
  // revealing than the faulting thread's own tag.
  fprintf(g_audit_log,
          "\n-- per-thread tag at fault (most recent TAG at seq < %lu) --\n",
          (unsigned long)cur);
  for (int t = 0; t < MAX_AUDIT_THREADS; t++)
  {
    uint64_t scan_lo = (cur > PLEN_RING_SIZE) ? cur - PLEN_RING_SIZE : 0;
    for (uint64_t idx = cur; idx-- > scan_lo; )
    {
      Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
      uint64_t ee = e.seq.load(std::memory_order_acquire);
      if (ee != idx) continue;
      if (e.tid != t) continue;
      if (e.kind != EV_TAG) continue;
      fprintf(g_audit_log,
              "  tid=%d seq=%lu tsc=%lu op=%s poly=%p slot=%d arg=%d\n",
              t, (unsigned long)ee, (unsigned long)e.tsc,
              e.op ? e.op : "<none>", e.poly, e.slot, e.arg);
      break;
    }
  }

  // Last 2048 events (thread-interleaved), newest first would be harder
  // to read — print oldest first so the immediate pre-fault events are
  // at the bottom.
  uint64_t nlo = (cur > 2048) ? cur - 2048 : 0;
  for (uint64_t idx = nlo; idx < cur; idx++)
  {
    Event &e = g_plen_ring[idx & (PLEN_RING_SIZE - 1)];
    uint64_t ee = e.seq.load(std::memory_order_acquire);
    if (ee != idx) continue;
    if (e.kind == EV_TAG)
      fprintf(g_audit_log,
              "  seq=%lu tsc=%lu tid=%d TAG op=%s poly=%p slot=%d arg=%d\n",
              (unsigned long)ee, (unsigned long)e.tsc,
              e.tid, e.op ? e.op : "<none>", e.poly, e.slot, e.arg);
    else
      fprintf(g_audit_log,
              "  seq=%lu tsc=%lu tid=%d PLEN T[%d] head=%p stored=%d actual=%d%s\n",
              (unsigned long)ee, (unsigned long)e.tsc,
              e.tid, e.slot, e.poly, e.arg, e.actual,
              (e.arg != e.actual) ? " *" : "");
  }
  fprintf(g_audit_log, "-- end dErrorBreak dump --\n");
  fflush(g_audit_log);
}

// Hook installers — kthread.cc plugs these into the reporter /
// kbuckets layers at bba_parallel_loop startup so libpolys can
// trigger our dump / tag without a reverse dependency on the kernel.
extern "C" {
extern void (*dErrorBreak_hook)(const char *reason);
extern void (*kbucket_debug_tag)(const char *op, void *lm, int slot, int arg);
extern void (*kbucket_debug_check_tnode)(const char *site, void *addr);
}

// Adapter with the right signature for the kbuckets hook.  The
// ring-buffer tag_log takes 5 args; we reuse it.
static void kbucket_debug_tag_adapter(const char *op, void *lm,
                                      int slot, int arg)
{
  tag_log(kt_debug_tid, op, lm, slot, arg);
}

static std::atomic<int> g_audit_hit_count{0};
static void audit_T_pLength(SweepContext *ctx, const char *tag, int thread_id)
{
  if (!g_audit_T_enabled.load(std::memory_order_relaxed)) return;
  kStrategy strat = ctx->strat;
  int n = strat->T.size();
  for (int j = 0; j < n; j++)
  {
    if (!tobject_published_load(strat->T[j])) continue;
    poly p = strat->T[j].p;
    if (p == NULL) p = strat->T[j].t_p;
    if (p == NULL) continue;
    int stored = strat->T[j].pLength;
    if (stored <= 0) continue;
    int actual = pLength(p);

    // Always log the observation to the ring buffer — lets us correlate
    // "same head observed with different lengths" after the fact.
    // Capture the seq so the dump knows exactly where in the ring the
    // mismatch sits (important because other threads keep logging).
    uint64_t my_seq = g_plen_seq.load(std::memory_order_relaxed);
    plen_log(thread_id, j, (void *)p, stored, actual);

    if (stored != actual)
    {
      int hit = g_audit_hit_count.fetch_add(1, std::memory_order_relaxed);
      if (hit < 20)  // log only first 20 hits per run
      {
        fprintf(g_audit_log,
                "[audit_T %s tid=%d] T[%d] stored pLength=%d, actual=%d "
                "(p=%p t_p=%p i_r=%d); Tsize=%d\n",
                tag, thread_id, j, stored, actual,
                (void*)strat->T[j].p, (void*)strat->T[j].t_p,
                strat->T[j].i_r, n);
        fflush(g_audit_log);
      }
      if (hit == 0)
      {
        // First hit: freeze the ring (other threads stop logging),
        // then dump and abort.  Freezing makes the dump see a stable
        // state — without it the ring wraps many times while the dump
        // is walking it.
        g_plen_frozen.store(true, std::memory_order_release);
        dump_debug_state(tag, thread_id, my_seq, j, (void *)p);
        abort();
      }
    }
  }
  g_audit_T_checked.fetch_add(1, std::memory_order_relaxed);
}

/* ------------------------------------------------------------------ */
/*  Lifecycle                                                          */
/* ------------------------------------------------------------------ */

SweepContext *sweep_context_init(kStrategy strat, int nthreads)
{
  SweepContext *ctx = (SweepContext *)calloc(1, sizeof(SweepContext));

  ctx->strat = strat;
  ctx->num_threads = nthreads;
  ctx->r = currRing;
  ctx->saved_si_opt_1 = si_opt_1;
  ctx->saved_si_opt_2 = si_opt_2;

  ctx->num_workers = nthreads - 1;
  {
    const char *wenv = getenv("SINGULAR_THREADS_WORKERS");
    if (wenv != NULL)
    {
      ctx->num_workers = atoi(wenv);
      if (ctx->num_workers < 0) ctx->num_workers = 0;
    }
  }

  // Milestone (d) pipeline depth. Default 16 or nthreads (whichever is
  // larger). Configurable via SINGULAR_PIPELINE_DEPTH for milestone (e)
  // to sweep without recompiling.
  {
    const char *penv = getenv("SINGULAR_PIPELINE_DEPTH");
    int depth = 16;
    if (penv != NULL) { depth = atoi(penv); if (depth < 1) depth = 1; }
    if (depth < nthreads) depth = nthreads;
    ctx->pipeline_depth = depth;
  }
  ctx->max_active = ctx->pipeline_depth;
  if (ctx->max_active < 1) ctx->max_active = 1;

  ctx->active = (ActivePoly *)calloc(ctx->max_active, sizeof(ActivePoly));
  for (int i = 0; i < ctx->max_active; i++)
  {
    ctx->active[i].occupied = false;
    ctx->active[i].is_survivor = false;
    ctx->active[i].state.store(SLOT_EMPTY, std::memory_order_relaxed);
    ctx->active[i].gen.store(0, std::memory_order_relaxed);
    ctx->active[i].needs_republish.store(false, std::memory_order_relaxed);
    ctx->active[i].tiles_remaining.store(0, std::memory_order_relaxed);
  }

  int total_threads = ctx->num_workers + 1;
  ctx->sweep_results = (SweepResult *)calloc(
      total_threads * ctx->max_active, sizeof(SweepResult));

  ctx->sweep_cursor.store(0, std::memory_order_relaxed);
  ctx->slot_counter.store(0, std::memory_order_relaxed);
  ctx->tile_cursor.store(0, std::memory_order_relaxed);
  ctx->tile_end.store(0, std::memory_order_relaxed);
  ctx->tiles_K = ctx->num_workers + 1;  // K slices per slot
  if (ctx->tiles_K < 1) ctx->tiles_K = 1;

  // Tile-batch ring. Size oversubscribed to accommodate many
  // re-publications per slot (multi-pass reductions).
  ctx->batch_ring_size = ctx->max_active * 64;
  if (ctx->batch_ring_size < 64) ctx->batch_ring_size = 64;
  ctx->batches = (TileBatch *)calloc(ctx->batch_ring_size, sizeof(TileBatch));
  ctx->batch_count.store(0, std::memory_order_relaxed);
  pthread_mutex_init(&ctx->publish_lock, NULL);
  pthread_cond_init(&ctx->tiles_avail_cv, NULL);
  pthread_cond_init(&ctx->slot_freed_cv, NULL);

  int barrier_count = ctx->num_workers + 1;
  if (barrier_count < 1) barrier_count = 1;
  pthread_barrier_init(&ctx->startup_barrier, NULL, barrier_count);

  int alloc_n = ctx->num_workers > 0 ? ctx->num_workers : 1;
  ctx->threads = (pthread_t *)calloc(alloc_n, sizeof(pthread_t));
  ctx->thread_ids = (int *)calloc(alloc_n, sizeof(int));

  pthread_mutex_init(&ctx->L_lock, NULL);

  ctx->done.store(false, std::memory_order_relaxed);
  ctx->stat_reductions.store(0, std::memory_order_relaxed);
  ctx->stat_zeros.store(0, std::memory_order_relaxed);
  ctx->stat_survivors.store(0, std::memory_order_relaxed);
  ctx->stat_rounds.store(0, std::memory_order_relaxed);

  // Asynchronous survivor drain: queue + mutex + CV + counter
  pthread_mutex_init(&ctx->survivor_queue_mutex, NULL);
  ctx->survivor_queue = new std::deque<LObject>();
  pthread_cond_init(&ctx->pairs_available, NULL);
  ctx->enterpairs_active.store(0, std::memory_order_relaxed);
  ctx->stat_max_queue_depth.store(0, std::memory_order_relaxed);

#ifdef KTHREAD_INSTRUMENT
  ctx->stats_enabled = (getenv("SINGULAR_KTHREAD_STATS") != NULL);
  int total_tt = ctx->num_workers + 1;
  ctx->tstats = (ThreadStats *)calloc(total_tt, sizeof(ThreadStats));
  ctx->rounds = new std::vector<RoundRecord>();
  ctx->reduces_this_round = new std::vector<ReduceEvent>();
  pthread_mutex_init(&ctx->stats_lock, NULL);
  ctx->start_ns = 0;
  ctx->workload_tag = getenv("SINGULAR_KTHREAD_TAG");
  if (ctx->workload_tag == NULL) ctx->workload_tag = "unknown";
#endif

  return ctx;
}

void sweep_context_destroy(SweepContext *ctx)
{
  if (ctx == NULL) return;
  pthread_barrier_destroy(&ctx->startup_barrier);
  pthread_mutex_destroy(&ctx->L_lock);
  pthread_mutex_destroy(&ctx->survivor_queue_mutex);
  pthread_mutex_destroy(&ctx->publish_lock);
  pthread_cond_destroy(&ctx->pairs_available);
  pthread_cond_destroy(&ctx->tiles_avail_cv);
  pthread_cond_destroy(&ctx->slot_freed_cv);
  delete ctx->survivor_queue;
  free(ctx->active);
  free(ctx->sweep_results);
  free(ctx->threads);
  free(ctx->thread_ids);
  free(ctx->batches);
#ifdef KTHREAD_INSTRUMENT
  free(ctx->tstats);
  delete ctx->rounds;
  delete ctx->reduces_this_round;
  pthread_mutex_destroy(&ctx->stats_lock);
#endif
  free(ctx);
}

/* ------------------------------------------------------------------ */
/*  Pop one polynomial from L and prepare for reduction.               */
/*  Caller must hold L_lock. Returns true if slot filled.              */
/* ------------------------------------------------------------------ */

/*
 * pop_and_prepare — pop next LObject from strat->L and set up ap for
 * reduction.  Caller holds ctx->L_lock.  The `sl_snapshot` parameter
 * is the T-size bound captured BEFORE taking L_lock (task 512
 * worker-side-drain); this avoids taking S-shared while holding L_lock,
 * which would invert the drainer's S → L order and risk deadlock.
 */
static BOOLEAN pop_and_prepare(SweepContext *ctx, ActivePoly *ap,
                               int sl_snapshot_arg)
{
  kStrategy strat = ctx->strat;

  while (!strat->L.empty())
  {
    if (strat->L.size() == 1) strat->interpt = TRUE;

    ap->P = strat->L.top();
    strat->L.pop();

    // pop-time sanity check.  pair.i_r{1,2} is an R-SLOT INDEX (stable
    // across T shifts — enterT updates R via the persistent .i_r
    // field).  The authoritative comparison is R[ir]->p.  T[ir].p
    // also reported as a secondary observation to spot T-shift
    // activity, but a mismatch there alone is NOT a bug.
    auto check_side = [&](const char *which, poly pp, int ir) {
      if (ir < 0 || ir >= strat->T.size()) return false;
      TObject *r_entry = (ir < (int)strat->T.size()) ? strat->R[ir] : NULL;
      poly exp_r = (r_entry != NULL) ? r_entry->p : NULL;
      poly exp_t = strat->T[ir].p;
      bool r_mismatch = (pp != NULL && exp_r != NULL && pp != exp_r);
      bool t_mismatch = (pp != NULL && pp != exp_t);
      if (r_mismatch || t_mismatch)
      {
        poly snap = (poly)kt_debug_lookup_T_head(ir);
        FILE *log = g_audit_log ? g_audit_log : stderr;
        fprintf(log,
                "\n=== POP CHECK (%s): pp=%p i_r=%d  "
                "R[ir]->p=%p (mismatch=%d, REAL BUG if 1)  "
                "T[ir].p=%p (mismatch=%d, cosmetic from T-shift)  "
                "snap[ir]=%p  T[ir-1].p=%p T[ir+1].p=%p  "
                "match_prev=%d match_next=%d  snap_eq_pp=%d "
                "snap_eq_tp=%d ===\n",
                which, (void*)pp, ir,
                (void*)exp_r, r_mismatch ? 1 : 0,
                (void*)exp_t, t_mismatch ? 1 : 0,
                (void*)snap,
                ir > 0 ? (void*)strat->T[ir-1].p : NULL,
                ir+1 < (int)strat->T.size() ? (void*)strat->T[ir+1].p : NULL,
                (ir > 0 && pp == strat->T[ir-1].p) ? 1 : 0,
                (ir+1 < (int)strat->T.size() && pp == strat->T[ir+1].p) ? 1 : 0,
                snap == pp ? 1 : 0,
                snap == exp_t ? 1 : 0);
        fflush(log);
      }
      // Only R-mismatch is treated as a real bug for gating
      // downstream tags.
      return r_mismatch;
    };
    bool p2_bad = check_side("p2", ap->P.p2, ap->P.i_r2);
    bool p1_bad = check_side("p1", ap->P.p1, ap->P.i_r1);
    if (p2_bad || p1_bad)
    {
      kt_debug_tag("pop:Pair_R_INCONSISTENT",
                   (void*)ap->P.p2, ap->P.i_r2,
                   (p1_bad ? 1 : 0) | (p2_bad ? 2 : 0));
      // Early-abort test: if SINGULAR_SKIP_BAD_POP is set, discard
      // this pair rather than passing it to ksReducePoly.  Lets us
      // verify whether the downstream SEGV is a consequence of the
      // R-inconsistent pair.
      if (getenv("SINGULAR_SKIP_BAD_POP") != NULL)
      {
        kt_debug_audit_printf(
          "=== SKIP_BAD_POP: discarding pair p1=%p i_r1=%d p2=%p "
          "i_r2=%d p1_bad=%d p2_bad=%d ===\n",
          (void*)ap->P.p1, ap->P.i_r1,
          (void*)ap->P.p2, ap->P.i_r2,
          p1_bad ? 1 : 0, p2_bad ? 1 : 0);
        kDeleteLcm(&ap->P);
        ap->P.Clear();
        continue;    // next iteration of `while (!strat->L.empty())`
      }
    }

    // Fingerprint experiment: recompute fp from the popped pair's
    // current (p1, p2, i_r1, i_r2) and compare to what was stamped
    // at creation.  fp==0 means the pair was created on a path that
    // didn't stamp (e.g., Ring variants, strong-poly variants) or
    // pre-dates the parallel window — skip those.  Only the
    // enterOnePairNormal path stamps.
    if (ap->P.dbg_fp != 0)
    {
      unsigned long now_fp = kt_debug_pair_fp((void*)ap->P.p1,
                                              (void*)ap->P.p2,
                                              ap->P.i_r1, ap->P.i_r2);
      bool fp_match = (now_fp == ap->P.dbg_fp);
      if (!fp_match || p1_bad || p2_bad)
      {
        FILE *log = g_audit_log ? g_audit_log : stderr;
        fprintf(log,
                "=== POP FP: stamped=%016lx now=%016lx match=%d  "
                "p1=%p p2=%p i_r1=%d i_r2=%d  "
                "p1_bad=%d p2_bad=%d ===\n",
                ap->P.dbg_fp, now_fp, fp_match ? 1 : 0,
                (void*)ap->P.p1, (void*)ap->P.p2,
                ap->P.i_r1, ap->P.i_r2,
                p1_bad ? 1 : 0, p2_bad ? 2 : 0);
        fflush(log);
        kt_debug_tag(fp_match ? "pop:fp_match" : "pop:fp_MISMATCH",
                     (void*)ap->P.p2, ap->P.i_r2, (int)ap->P.i_r1);
      }
    }

    // Create spoly if needed
    if (pNext(ap->P.p) == strat->tail)
    {
      if (rField_is_Ring(currRing))
        pLmDelete(ap->P.p);
      else
        pLmFree(ap->P.p);
      ap->P.p = NULL;
      poly m1 = NULL, m2 = NULL;
      while (strat->tailRing != currRing &&
             !kCheckSpolyCreation(&(ap->P), strat, m1, m2))
      {
        assume(m1 == NULL && m2 == NULL);
        break;
      }
      // Thread safety: ksCreateSpoly temporarily modifies Pair->p1/p2
      // via p_SetCompP when their components differ (module computations).
      // Since p1/p2 point to shared T entry polynomials that other threads
      // may be reading concurrently in ksReducePoly/pp_Mult_mm, we must
      // make deep copies so that p_SetCompP modifies copies, not shared data.
      bool need_copy = (currRing->pCompIndex >= 0) &&
          (ap->P.p1 != NULL) && (ap->P.p2 != NULL) &&
          (__p_GetComp(ap->P.p1, currRing) != __p_GetComp(ap->P.p2, currRing));

      if (need_copy)
      {
        poly orig_p1 = ap->P.p1;
        poly orig_p2 = ap->P.p2;
        ap->P.p1 = pCopy(orig_p1);
        ap->P.p2 = pCopy(orig_p2);

        ksCreateSpoly(&(ap->P), NULL, strat->use_buckets,
                      strat->tailRing, m1, m2, &strat->R);

        // Free the copies (ksCreateSpoly doesn't consume p1/p2, and the
        // spoly data is all freshly allocated from pp_Mult_mm etc.)
        pDelete(&ap->P.p1);
        pDelete(&ap->P.p2);
        ap->P.p1 = orig_p1;
        ap->P.p2 = orig_p2;
      }
      else
      {
        ksCreateSpoly(&(ap->P), NULL, strat->use_buckets,
                      strat->tailRing, m1, m2, &strat->R);
      }
    }
    else if (ap->P.p1 == NULL)
    {
      poly _tp = (ap->P.t_p != NULL) ? ap->P.t_p : ap->P.p;
      poly _next = (_tp != NULL) ? pNext(_tp) : NULL;
      // Log TWO tag events: (head, pNext) so we can find the mutation
      // via either the first or second node pointer.
      OpTag _pr(0, "PrepareRed(ap->P,p1==NULL) head", (void *)_tp);
      OpTag _pr2(0, "PrepareRed(ap->P,p1==NULL) pNext", (void *)_next);
      ap->P.PrepareRed(strat->use_buckets);
    }

    if ((ap->P.p == NULL) && (ap->P.t_p == NULL))
      continue;  // zero, try next

    // Initialize
    ap->P.SetLmCurrRing();
    ap->P.SetShortExpVector();
    ap->not_sev = ~ap->P.sev;
    ap->pass = 0;
    ap->d = ap->P.GetpFDeg() + ap->P.ecart;
    ap->reddeg = ap->d;
    ap->occupied = true;
    ap->is_survivor = false;
    ap->best_reducer = -1;
    ap->best_good = -1;

    // Task 280 milestone (a): record current T bound for this slot.
    // Task 512 worker-side-drain: sl_snapshot_arg was captured in the
    // caller BEFORE L_lock was taken, under S-shared — ensures a
    // consistent view of T.size() without inverting the lock order
    // (phase-1 drainers take S-shared then L-exclusive; main acquires
    // S-shared briefly and releases before L_lock).  The captured
    // value may be slightly stale by the time we use it, but enterT
    // only grows T monotonically and any entries added post-capture
    // are picked up in the next refill pass.
    ap->sl_snapshot = sl_snapshot_arg;

    {
      poly _tp = (ap->P.t_p != NULL) ? ap->P.t_p : ap->P.p;
      poly _next = (_tp != NULL) ? pNext(_tp) : NULL;
      OpTag _pr(0, "PrepareRed(ap->P,final) head", (void *)_tp);
      OpTag _pr2(0, "PrepareRed(ap->P,final) pNext", (void *)_next);
      ap->P.PrepareRed(strat->use_buckets);
    }
    return TRUE;
  }

  return FALSE;
}

/* ------------------------------------------------------------------ */
/*  Sweep phase: tile cursor with closer-reduce (milestone c, task 282)*/
/*                                                                     */
/*  A tile is a pair (slot, slice). For each occupied slot, its        */
/*  sl_snapshot range [0..sl_snapshot] is split into K slices where    */
/*  slice `i` visits j with (j % K == i). Workers pull tile ids from   */
/*  an atomic cursor; decoding is slot = id / K, slice = id % K.       */
/*                                                                     */
/*  Per-slot SweepResult remains per-thread to avoid contention (a     */
/*  single slot's K slices may be claimed by K different workers).     */
/*  After sweeping a tile, the worker fetch_sub(1)s the slot's         */
/*  tiles_remaining counter. The worker that drops it to zero is the   */
/*  "closer": it merges the K per-thread SweepResults for the slot     */
/*  and immediately calls reduce_slot_from_sweep on it (closer-        */
/*  reduces). Survivors are pushed to the survivor FIFO via            */
/*  queue_survivor; the main thread drains them between B1 and the     */
/*  next B0 as in milestone (b).                                       */
/*                                                                     */
/*  Safety: reduce_slot_from_sweep writes only slot-local state        */
/*  (ap->P, ap->not_sev, ap->d, ap->pass, ...) plus idempotent-store   */
/*  of strat->overflow and atomic stat_* counters. It reads strat->T   */
/*  (stable while enterT is serialized on main between B1 and B0) and  */
/*  calls ksReducePoly (operates on ap->P and read-only T entry).      */
/*  No two workers can close the same slot (fetch_sub(1)==1 is         */
/*  single-winner). Workers closing *different* slots run in           */
/*  parallel; their only shared writes are the atomic counters.        */
/* ------------------------------------------------------------------ */

// Forward declarations: closer-reduce calls these from inside sweep_phase.
static void reduce_slot_from_sweep(SweepContext *ctx, int slot, int thread_id);
static void queue_survivor(SweepContext *ctx, ActivePoly *ap);

/*
 * Merge the per-thread SweepResults for one slot into ap->best_*.
 * Called by the single worker that closes the slot's tiles_remaining
 * counter, so no locking is required: the acq_rel fetch_sub that
 * selects the closer synchronizes-with every other worker's release
 * of its per-thread SweepResult for this slot.
 */
static void merge_slot_results(SweepContext *ctx, int s)
{
  int total_threads = ctx->num_workers + 1;
  int best_reducer = -1, best_good = -1, best_pLength = 0;
  void *best_reducer_p = NULL, *best_good_p = NULL;

  for (int t = 0; t < total_threads; t++)
  {
    SweepResult &sr = sweep_result(ctx, t, s);
    if (sr.best_reducer >= 0 && best_reducer < 0)
    {
      best_reducer = sr.best_reducer;
      best_reducer_p = sr.best_reducer_p;
    }
    if (sr.best_good >= 0)
    {
      if (best_good < 0 || sr.best_pLength < best_pLength)
      {
        best_good = sr.best_good;
        best_good_p = sr.best_good_p;
        best_pLength = sr.best_pLength;
      }
    }
  }

  ctx->active[s].best_reducer = best_reducer;
  ctx->active[s].best_good = best_good;
  ctx->active[s].best_pLength = best_pLength;
  ctx->active[s].best_reducer_p = best_reducer_p;
  ctx->active[s].best_good_p = best_good_p;
}

/*
 * Reset per-thread sweep results for ONE slot across all threads.
 * Called immediately before publishing a fresh batch for a slot.
 */
static void reset_sweep_results_one(SweepContext *ctx, int slot)
{
  int total_threads = ctx->num_workers + 1;
  for (int t = 0; t < total_threads; t++)
  {
    SweepResult &sr = sweep_result(ctx, t, slot);
    sr.best_reducer = -1;
    sr.best_good = -1;
    sr.best_pLength = 0;
    sr.best_reducer_p = NULL;
    sr.best_good_p = NULL;
  }
}

/*
 * Publish a fresh tile batch for slot `s`. Caller holds publish_lock.
 * Increments the slot's generation, resets per-slot sweep results,
 * resets tiles_remaining to K, writes the batch record, and advances
 * tile_end. On exit, tiles_avail_cv is broadcast so waiting workers
 * wake up.
 *
 * Pre-conditions:
 *   - slot state == SLOT_FILLED (data initialised by pop_and_prepare
 *     or kept occupied after reduce)
 *   - tiles_remaining is either 0 (fresh slot) or already 0 after
 *     the previous close_slot (closer was the one that dropped it)
 */
static void publish_slot_tiles_locked(SweepContext *ctx, int s)
{
  ActivePoly *ap = &ctx->active[s];
  int K = ctx->tiles_K;

  uint64_t new_gen = ap->gen.fetch_add(1, std::memory_order_relaxed) + 1;
  ap->needs_republish.store(false, std::memory_order_relaxed);

  reset_sweep_results_one(ctx, s);
  ap->tiles_remaining.store(K, std::memory_order_relaxed);

  // Slot data must be visible before the batch is observable: use
  // release below when advancing tile_end.
  uint64_t b = ctx->batch_count.fetch_add(1, std::memory_order_relaxed);
  TileBatch &tb = ctx->batches[b % ctx->batch_ring_size];
  tb.slot = s;
  tb.gen = new_gen;
  tb.tl_snapshot = ap->sl_snapshot;

  // Publish slot state as FILLED before advancing tile_end. Workers
  // that claim a tile from this batch must observe state=FILLED.
  ap->state.store(SLOT_FILLED, std::memory_order_release);

  // Release: readers acquiring tile_end see batch fields initialised,
  // slot data (ap->P, not_sev, sl_snapshot), tiles_remaining, state.
  ctx->tile_end.store((b + 1) * (uint64_t)K, std::memory_order_release);
}

/*
 * Acquire publish_lock and publish, then broadcast tiles_avail_cv
 * so workers wake.
 */
static void publish_slot_tiles(SweepContext *ctx, int s)
{
  pthread_mutex_lock(&ctx->publish_lock);
  publish_slot_tiles_locked(ctx, s);
  pthread_cond_broadcast(&ctx->tiles_avail_cv);
  pthread_mutex_unlock(&ctx->publish_lock);
}

/*
 * Closer-reduce: called by the worker that drops tiles_remaining for
 * slot `s` to zero. Merges results and applies one ksReducePoly step.
 * If the slot becomes a survivor, pushes onto the survivor FIFO.
 * If the slot is zero'd or turned into a survivor, transitions the
 * slot state to SLOT_EMPTY and signals slot_freed_cv. If the slot is
 * still occupied (needs another sweep pass), sets needs_republish
 * so main will republish tiles for it.
 */
static void close_slot(SweepContext *ctx, int s, int thread_id)
{
  OpTag _close_tag(thread_id, "close_slot", NULL, s);
  ActivePoly *ap = &ctx->active[s];

  // Unoccupied or already-survivor: should not happen in the continuous
  // design (batches only fire for FILLED slots) but keep as a safety.
  if (!ap->occupied || ap->is_survivor) return;

  // Time-bisection: check L for R-mismatch BEFORE the work in close_slot.
  kt_debug_R_scan_L(ctx->strat, "close_slot:entry");

  merge_slot_results(ctx, s);

#ifdef KTHREAD_INSTRUMENT
  long rstart = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
  reduce_slot_from_sweep(ctx, s, thread_id);
#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx))
  {
    long rdur = kt_now_ns() - rstart;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.reduce_ns += rdur;
    ts.reduce_count++;
    kt_record_reduce(ctx, thread_id, rstart - ctx->start_ns, rdur);
  }
#endif

  // After reduce_slot_from_sweep, three possible outcomes:
  //   (1) poly reduced to zero / overflow / syzComp-out → ap->occupied = false
  //   (2) survivor               → ap->is_survivor = true, still occupied
  //   (3) still occupied & not survivor → needs another sweep pass
  if (ap->occupied && ap->is_survivor)
  {
    // Queue survivor for main to drain; transfer ownership out of slot.
    queue_survivor(ctx, ap);
    // queue_survivor clears the slot. Publish state=EMPTY.
    ap->state.store(SLOT_EMPTY, std::memory_order_release);
    pthread_mutex_lock(&ctx->publish_lock);
    pthread_cond_broadcast(&ctx->slot_freed_cv);
    pthread_mutex_unlock(&ctx->publish_lock);
  }
  else if (!ap->occupied)
  {
    // Slot cleared (zero or overflow). State → EMPTY.
    ap->state.store(SLOT_EMPTY, std::memory_order_release);
    pthread_mutex_lock(&ctx->publish_lock);
    pthread_cond_broadcast(&ctx->slot_freed_cv);
    pthread_mutex_unlock(&ctx->publish_lock);
  }
  else
  {
    // Needs another sweep pass on the same slot. Mark for main to
    // republish; main will pick this up in its loop and call
    // publish_slot_tiles(s).
    ap->needs_republish.store(true, std::memory_order_release);
    pthread_mutex_lock(&ctx->publish_lock);
    pthread_cond_broadcast(&ctx->slot_freed_cv);
    pthread_mutex_unlock(&ctx->publish_lock);
  }

  // And at exit: if R-mismatch first appears here, close_slot's body
  // (merge_slot_results + reduce_slot_from_sweep + queue_survivor) is
  // the writer.
  kt_debug_R_scan_L(ctx->strat, "close_slot:exit");
}

/*
 * Sweep one tile (target slot + slice). Pure read-only access to
 * strat->T up to tl_snapshot; writes only to the worker's private
 * SweepResult.
 *
 * Task 511 t-iterator-migrate-hot-path: the K-stride access pattern
 * visits indices j, j+K, j+2K, ..., which does not fit a sequential
 * iterator walk.  Instead we keep the integer j loop and add an
 * explicit tobject_published_load gate before dereferencing T[j]
 * fields.  sevT[j] is written by enterT before tobject_publish(T[j]),
 * so sev pre-filtering is safe on an in-flight slot (worst case the
 * read observes the old sev=0 and we fall through to the published
 * gate).  The acquire-load on published synchronises-with enterT's
 * release-publish on any concurrent drainer (future worker-side drain
 * task), so T[j].p / .ecart / .pLength are guaranteed visible once
 * the gate returns true.
 */
static void sweep_one_tile(SweepContext *ctx, int thread_id,
                           int s, int slice, int tl_snapshot)
{
  OpTag _sw_tag(thread_id, "sweep_one_tile", NULL, s, slice);
  kStrategy strat = ctx->strat;
  int K = ctx->tiles_K;
  ActivePoly *ap = &ctx->active[s];

  // If the slot has since transitioned to EMPTY (race with closer on
  // another slot batch?), skip. Should not happen: batches are only
  // published while state=FILLED and the closer fires after all K
  // tiles execute, clearing state only afterwards. Defensive anyway.
  if (ap->state.load(std::memory_order_acquire) != SLOT_FILLED)
    return;
  if (!ap->occupied || ap->is_survivor) return;

  unsigned long not_sev_s = ap->not_sev;
  SweepResult &sr = sweep_result(ctx, thread_id, s);

  for (int j = slice; j <= tl_snapshot; j += K)
  {
    unsigned long sev_j = strat->sevT[j];
    if (sev_j & not_sev_s) continue;
    // Acquire-load gate: only after observing published=true are the
    // T[j] field reads below guaranteed to synchronise-with the
    // release-store in enterT.
    if (!tobject_published_load(strat->T[j])) continue;
    if (!p_LmDivisibleBy(strat->T[j].p, ap->P.p, currRing))
      continue;

    if (sr.best_reducer < 0)
    {
      sr.best_reducer = j;
      sr.best_reducer_p = (void*)strat->T[j].p;
    }

    int ecart_j = strat->T[j].ecart;
    if (ecart_j <= ap->P.ecart)
    {
      int pLen = strat->T[j].pLength;
      if (pLen <= 0) pLen = 3;
      if (sr.best_good < 0 || pLen < sr.best_pLength)
      {
        sr.best_good = j;
        sr.best_good_p = (void*)strat->T[j].p;
        sr.best_pLength = pLen;
      }
    }
  }
}

/*
 * Continuous tile-pull loop. Called by worker_thread and also by
 * main between fill/drain cycles (pump_tiles_main) so main helps
 * sweep when there is nothing else for it to do.
 *
 * Termination: loops until ctx->done is set AND no pending tiles
 * remain. Workers called with block=true wait on tiles_avail_cv when
 * the cursor catches up to end; main calls with block=false and
 * returns as soon as the cursor is caught up, so it can go back to
 * refilling / draining.
 */
// Forward declaration — workers in tile_pull_loop check this between
// tile batches.  Defined later in this file.
static void drain_survivor_queue(SweepContext *ctx, int thread_id);

static void tile_pull_loop(SweepContext *ctx, int thread_id, bool block)
{
  OpTag _tp_tag(thread_id, block ? "tile_pull_loop(block)" : "tile_pull_loop(nb)");
  int K = ctx->tiles_K;

  while (true)
  {
    uint64_t end = ctx->tile_end.load(std::memory_order_acquire);
    uint64_t cur = ctx->tile_cursor.load(std::memory_order_acquire);

    if (cur >= end)
    {
      if (!block) return;
      if (ctx->done.load(std::memory_order_acquire)) return;
      // Task 512 worker-side-drain: before going to sleep, check whether
      // the survivor queue has work.  If so, drain instead of cond_wait.
      // Only workers do this (block=true); main has its own drain pass
      // in bba_parallel_loop and calls tile_pull_loop with block=false.
      //
      // Gating: skip the drain hop if the queue is observably empty
      // (snapshot under its own mutex so we don't race with a concurrent
      // push).  The cost of the lock is negligible vs the wait.
      if (getenv("SKIP_WORKER_DRAIN") == NULL)
      {
        bool queue_has_work;
        kt_surv_q_lock(ctx, thread_id);
        queue_has_work = !ctx->survivor_queue->empty();
        pthread_mutex_unlock(&ctx->survivor_queue_mutex);
        if (queue_has_work)
        {
#ifdef KTHREAD_INSTRUMENT
          long idle_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
          drain_survivor_queue(ctx, thread_id);
#ifdef KTHREAD_INSTRUMENT
          if (KT_STATS(ctx))
          {
            ThreadStats &ts = KT_TS(ctx, thread_id);
            ts.worker_drain_idle_ns += kt_now_ns() - idle_t0;
            ts.worker_drain_idle_count++;
          }
#endif
          continue;
        }
      }
      // Wait for more tiles or shutdown.
      pthread_mutex_lock(&ctx->publish_lock);
      // Re-check under lock.
      uint64_t end2 = ctx->tile_end.load(std::memory_order_acquire);
      uint64_t cur2 = ctx->tile_cursor.load(std::memory_order_acquire);
      bool done = ctx->done.load(std::memory_order_acquire);
      if (cur2 >= end2 && !done)
        pthread_cond_wait(&ctx->tiles_avail_cv, &ctx->publish_lock);
      pthread_mutex_unlock(&ctx->publish_lock);
      continue;
    }

    // Try to claim id via CAS.
    if (!ctx->tile_cursor.compare_exchange_weak(
          cur, cur + 1,
          std::memory_order_acq_rel, std::memory_order_acquire))
      continue;  // another thread took it; retry
    uint64_t id = cur;

    uint64_t b = id / (uint64_t)K;
    int slice = (int)(id % (uint64_t)K);
    TileBatch &tb = ctx->batches[b % ctx->batch_ring_size];
    int s = tb.slot;
    int tl = tb.tl_snapshot;
    uint64_t expected_gen = tb.gen;

    ActivePoly *ap = &ctx->active[s];
    uint64_t cur_gen = ap->gen.load(std::memory_order_acquire);

    if (cur_gen == expected_gen)
    {
      sweep_one_tile(ctx, thread_id, s, slice, tl);
    }
    // else: stale batch (should not happen given pipeline constraints) —
    // skip the sweep but still decrement tiles_remaining. But
    // tiles_remaining belongs to the *current* gen; decrementing it
    // would break the closer protocol. This is why we sized the ring
    // generously and only publish a new batch after tiles_remaining==0.

    // Decrement tiles_remaining; closer if we hit zero.
    if (ap->tiles_remaining.fetch_sub(1, std::memory_order_acq_rel) == 1)
      close_slot(ctx, s, thread_id);
  }
}

/* ------------------------------------------------------------------ */
/*  Merge per-thread sweep results into per-slot best reducer.         */
/*  (merge_sweep_results removed in milestone c of task 282:            */
/*   merging is now done per-slot by the closer worker via              */
/*   merge_slot_results inside sweep_phase.)                             */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/*  Reset sweep results for all threads and all slots.                 */
/* ------------------------------------------------------------------ */

/* reset_sweep_results removed in milestone (d): per-slot reset in
 * publish_slot_tiles_locked via reset_sweep_results_one. */

/* ------------------------------------------------------------------ */
/*  Reduce one slot: apply ONE ksReducePoly step from sweep result.    */
/*  Slot stays active for the next cooperative sweep round if the      */
/*  polynomial is non-zero and not a survivor.                         */
/*  Called by worker threads in parallel (one slot per thread).        */
/* ------------------------------------------------------------------ */

static void reduce_slot_from_sweep(SweepContext *ctx, int slot, int thread_id)
{
  OpTag _reduce_tag(thread_id, "reduce_slot_from_sweep", NULL, slot);
  audit_T_pLength(ctx, "reduce-entry", thread_id);
  kStrategy strat = ctx->strat;
  ActivePoly *ap = &ctx->active[slot];

  int best = (ap->best_good >= 0) ? ap->best_good : ap->best_reducer;

  if (best < 0)
  {
    // No reducer found by cooperative sweep — survivor
    ap->is_survivor = true;
    return;
  }

  // Apply first reduction from cooperative sweep result
  int ei = strat->T[best].ecart;

  // Consistency guard: verify T[best].p still matches what the sweep
  // observed (captured at merge time).  If it doesn't, either:
  //  - a concurrent enterT shifted T[best] (shouldn't happen under
  //    posInT_appendEnd + S-exclusive lock during enterT), or
  //  - some other path mutated T[best].p in-place,
  //  - or `best` is out of the currently-valid T range and we're
  //    reading uninitialised memory.
  // All three are bugs.  This is the earliest point we can detect
  // before the ksReducePoly call does the wild deref.
  {
    void *seen_at_sweep =
      (ap->best_good >= 0) ? ap->best_good_p : ap->best_reducer_p;
    poly now_at_reduce =
      (best >= 0 && best < (int)strat->T.size()) ? strat->T[best].p
                                                 : (poly)NULL;
    if (seen_at_sweep != NULL && seen_at_sweep != (void*)now_at_reduce)
    {
      kt_debug_audit_printf(
        "=== REDUCE_GUARD: T[best=%d].p changed between sweep and "
        "reduce!  seen_at_sweep=%p  now_at_reduce=%p  "
        "slot=%d thread=%d  T.size=%d  "
        "best_good=%d best_reducer=%d ===\n",
        best, seen_at_sweep, (void*)now_at_reduce,
        slot, thread_id, (int)strat->T.size(),
        ap->best_good, ap->best_reducer);
      kt_debug_tag("REDUCE_GUARD:T[best].p_CHANGED",
                   seen_at_sweep, best, slot);
    }
    if (best < 0 || best >= (int)strat->T.size())
    {
      kt_debug_audit_printf(
        "=== REDUCE_GUARD: best=%d OUT OF BOUNDS  T.size=%d  "
        "slot=%d thread=%d  seen_at_sweep=%p ===\n",
        best, (int)strat->T.size(), slot, thread_id, seen_at_sweep);
    }
  }

  // Time-bisection: scan L for R-inconsistency BEFORE ksReducePoly.
  // If 0, the writer hasn't struck yet (or has struck and recovered).
  kt_debug_R_scan_L(strat, "reduce_slot:pre_ksReducePoly");

  {
    // Tag includes the T[best] reducer so the dump shows exactly which
    // T entry this worker is consuming when the audit fires on some
    // OTHER thread.
    OpTag _r(thread_id, "ksReducePoly(using T[best])",
             (void *)strat->T[best].p, slot, best);
    ksReducePoly(&ap->P, strat->T.addr(best),
                 strat->kNoetherTail(), NULL, NULL, strat);
  }

  // And AFTER ksReducePoly: did THIS reduce step introduce the
  // R-mismatch?  If 0 before and N>0 after, ksReducePoly is the
  // writer.
  kt_debug_R_scan_L(strat, "reduce_slot:post_ksReducePoly");

  ctx->stat_reductions.fetch_add(1, std::memory_order_relaxed);

  if (ap->P.IsNull())
  {
    kDeleteLcm(&ap->P);
    ap->P.Clear();
    ap->occupied = false;
    ctx->stat_zeros.fetch_add(1, std::memory_order_relaxed);
    return;
  }

  // IDLIFT check
  if (UNLIKELY(TEST_OPT_IDLIFT))
  {
    poly hp = ap->P.p ? ap->P.p : ap->P.t_p;
    if (hp && p_GetComp(hp, currRing) > strat->syzComp)
    {
      ap->P.Delete();
      ap->occupied = false;
      ctx->stat_zeros.fetch_add(1, std::memory_order_relaxed);
      return;
    }
  }
  else if (UNLIKELY((strat->syzComp > 0) && (!TEST_OPT_REDTAIL_SYZ)))
  {
    poly hp = ap->P.p ? ap->P.p : ap->P.t_p;
    if (hp && p_GetComp(hp, currRing) > strat->syzComp)
    {
      ap->is_survivor = true;
      return;
    }
  }

  // Update ecart after first reduction
  ap->P.SetShortExpVector();
  int h_d = ap->P.SetpFDeg();
  if (ei <= ap->P.ecart)
    ap->P.ecart = ap->d - h_d;
  else
    ap->P.ecart = ap->d - h_d + ei - ap->P.ecart;
  ap->pass++;
  ap->d = h_d + ap->P.ecart;

  if (UNLIKELY(ap->d > ap->reddeg))
  {
    if (UNLIKELY(ap->d >= (long)strat->tailRing->bitmask))
    {
      if (ap->P.pTotalDeg() + ap->P.ecart >= (long)strat->tailRing->bitmask)
      {
        strat->overflow = TRUE;
        ap->P.GetP();
        ap->P.Clear();
        ap->occupied = false;
        return;
      }
    }
    ap->reddeg = ap->d;
  }

  ap->P.SetLmCurrRing();
  ap->P.SetShortExpVector();
  ap->not_sev = ~ap->P.sev;

  // Slot stays occupied — will be swept again on the next round
}

/* ------------------------------------------------------------------ */
/*  (reduce_phase_parallel removed in milestone c of task 282:          */
/*   reduction now happens inline in sweep_phase via close_slot when    */
/*   a worker drops the slot's tiles_remaining counter to zero.)        */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/*  Refill empty slots from L, then publish tile batches for any     */
/*  slot that is SLOT_EMPTY (just filled) or has needs_republish     */
/*  (closer left it occupied and needing another sweep). Returns the */
/*  number of slots currently FILLED (in pipeline) after the pass.   */
/*  Main-thread only.                                                 */
/* ------------------------------------------------------------------ */

static int refill_and_publish(SweepContext *ctx)
{
  OpTag _refill_tag(0, "refill_and_publish");
  audit_T_pLength(ctx, "refill-entry", 0);
  int in_pipeline = 0;

  for (int s = 0; s < ctx->max_active; s++)
  {
    ActivePoly *ap = &ctx->active[s];
    int st = ap->state.load(std::memory_order_acquire);

    if (st == SLOT_EMPTY)
    {
      // Task 512 worker-side-drain: capture sl_snapshot under S-shared
      // BEFORE taking L_lock, to avoid inverting the worker drainer's
      // S → L lock order.
      kt_S_lock_shared(ctx, 0);
      int sl_snapshot = ctx->strat->T.size() - 1;
      ctx->strat->S.unlock_shared();

      // Try to fill from L.
      kt_L_lock(ctx, 0);
      BOOLEAN got = pop_and_prepare(ctx, ap, sl_snapshot);
      pthread_mutex_unlock(&ctx->L_lock);
      if (got)
      {
        ctx->stat_rounds.fetch_add(1, std::memory_order_relaxed);
        // pop_and_prepare set ap->occupied=true; now publish.
        publish_slot_tiles(ctx, s);
        in_pipeline++;
      }
      continue;
    }
    // state == SLOT_FILLED
    if (ap->needs_republish.load(std::memory_order_acquire))
    {
      // Closer-reduce left this slot needing another sweep pass.
      // Re-snapshot tl under S-shared; no L_lock held here, so no
      // inversion risk.
      kt_S_lock_shared(ctx, 0);
      ap->sl_snapshot = ctx->strat->T.size() - 1;
      ctx->strat->S.unlock_shared();
      publish_slot_tiles(ctx, s);
    }
    in_pipeline++;
  }

  return in_pipeline;
}

/* ------------------------------------------------------------------ */
/*  Process survivor LObject                                           */
/*                                                                     */
/*  Runs the serial enterT / enterpairs / enterS sequence for one      */
/*  survivor polynomial. The caller guarantees mutual exclusion (at    */
/*  most one thread runs this at a time) via ctx->enterpairs_mutex.    */
/*                                                                     */
/*  This function does NOT touch any ActivePoly slot. The caller       */
/*  moves/copies the LObject into the queue and clears the slot        */
/*  independently.                                                     */
/* ------------------------------------------------------------------ */

/*
 * Phase-0 / Phase-1 design (task 508 enterpairs-parallel-phase1).
 *
 * Task 506 landed infrastructure (rwlock on sBasisSet, atomic CAS on
 * tombstone flags, explicit phase 0/1/2 structure, per-phase
 * instrumentation) but kept phase 1 under an S-exclusive lock because
 * enterpairs writes into the strat-global strat->B.  This task enables
 * concurrent phase-1 drainers under a SHARED S-lock by introducing
 * thread-local overrides (see kutil.h):
 *
 *   t_local_B_override       — thread-local LSet* replacing strat->B
 *   t_local_my_arrival       — filter S iteration by arrival_id
 *   t_local_pairtest_hits    — per-drainer pairtest-hit vector
 *
 * Phase 2 (the old "L merge" step) no longer exists as a separate
 * phase: initenterpairs already merges B into L via chainCritNormal →
 * kMergeBintoL.  The old phase-2 block was just enterS, which is now
 * in phase 0.
 *
 * Lock discipline:
 *
 *   Phase 0 (S-exclusive, brief): GetP / initEcart / find_pos /
 *     redtailBba / pCleardenom / pNorm / SetShortExpVector / enterT /
 *     my_arrival = arrival_counter.load() / enterS(h).  enterS uses
 *     fetch_add(1) internally; because only one thread is in phase 0
 *     under S-exclusive at a time, this fetch_add returns
 *     my_arrival, stamping the new SElement with the id we loaded.
 *
 *   Phase 1 (S-shared + L-exclusive): set up thread-local B /
 *     pairtest_hits / my_arrival; call enterpairs (walks S under
 *     shared lock, tombstone-only erases for clearS, builds local B,
 *     chainCrit merges local B into strat->L under L-lock).  Release
 *     L then S-shared.
 *
 * SORDER_STANDARD correctness note: enterS inserts h at a sorted
 * position via find_pos; any subsequent enterS-inserted elements (by
 * peer drainers) go to THEIR own sorted positions — but they can
 * only happen AFTER my phase 1 releases S-shared (they need
 * S-exclusive, they wait).  So during my phase 1, S is not mutated
 * (except via atomic CAS tombstones, which don't shift indices), and
 * my pos_it iterator stays valid.
 *
 * Peer drainer concurrency during phase 1 is shared-lock + tombstone,
 * so reads are safe.  What I've added to S in phase 0 is visible to
 * peer drainers' phase-1 reads (they hold shared lock; writes happened
 * under exclusive before they acquired shared — hb ordered).
 */
static void process_survivor_lobject(SweepContext *ctx, LObject *P, int thread_id)
{
  OpTag _ps_tag(thread_id, "process_survivor_lobject", P ? P->p : NULL);
  kStrategy strat = ctx->strat;
  BOOLEAN withT = ctx->withT;

#ifdef KTHREAD_INSTRUMENT
  long ps_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
  long redtail_accum = 0;
  long enterT_accum = 0;
  long enterpairs_accum = 0;
  long enterS_accum = 0;
  long phase0_work = 0;
  long phase1_work = 0;
#endif

  // ------------------------------------------------------------------
  // Phase 0 — S-exclusive, brief: everything that must run under
  // exclusive semantics.  GetP / initEcart / find_pos / redtailBba /
  // pCleardenom / pNorm / SetShortExpVector / enterT / enterS.
  //
  // Capturing my_arrival: because we hold S-exclusive, strat->arrival_counter
  // reflects the id of the NEXT entry enterS will take.  Read it
  // here; enterS then fetch_add's it, stamping the new SElement.
  // ------------------------------------------------------------------
  kt_S_lock_exclusive(ctx, thread_id);
#ifdef KTHREAD_INSTRUMENT
  long p0_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif

  P->GetP(strat->lmBin);
  if (strat->homog) strat->initEcart(P);

  if (TEST_OPT_PROT) PrintS("s");

  auto pos_it = strat->S.find_pos(P->p, P->ecart);

  strat->redTailChange = FALSE;

  audit_T_pLength(ctx, "ps-phase0-pre-redtailBba", thread_id);

  // Task 512 bisect: disable redtailBba in drain via env var.
  bool skip_redtail = (getenv("SKIP_DRAIN_REDTAIL") != NULL);

  if (rField_is_Z(currRing) && !rHasLocalOrMixedOrdering(currRing))
    redtailBbaAlsoLC_Z(P, strat);

  if (TEST_OPT_INTSTRATEGY)
  {
    P->pCleardenom();
    if ((TEST_OPT_REDSB) || (TEST_OPT_REDTAIL))
    {
#ifdef KTHREAD_INSTRUMENT
      long rt0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
      if (!skip_redtail)
      {
        OpTag _r(thread_id, "redtailBba(Z)", P->p);
        P->p = redtailBba(P, pos_it, strat, withT,
                          !TEST_OPT_CONTENTSB);
      }
#ifdef KTHREAD_INSTRUMENT
      if (KT_STATS(ctx)) redtail_accum += kt_now_ns() - rt0;
#endif
      P->pCleardenom();
      if (strat->redTailChange) P->t_p = NULL;
    }
  }
  else
  {
    P->pNorm();
    if ((TEST_OPT_REDSB) || (TEST_OPT_REDTAIL))
    {
#ifdef KTHREAD_INSTRUMENT
      long rt0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
      if (!skip_redtail)
      {
        OpTag _r(thread_id, "redtailBba", P->p);
        P->p = redtailBba(P, pos_it, strat, withT);
      }
#ifdef KTHREAD_INSTRUMENT
      if (KT_STATS(ctx)) redtail_accum += kt_now_ns() - rt0;
#endif
      if (strat->redTailChange) P->t_p = NULL;
    }
  }

  audit_T_pLength(ctx, "ps-phase0-post-redtailBba", thread_id);

  // Enter h into T and S under exclusive lock.  my_arrival is the
  // arrival_id enterS will stamp on h (enter_bba fetch_add's the
  // counter; we're the only thread doing enterS under S-exclusive).
  uint64_t my_arrival = UINT64_MAX;
  bool did_enterS = false;
  if ((!TEST_OPT_IDLIFT) || (pGetComp(P->p) <= strat->syzComp))
  {
    P->SetShortExpVector();
#ifdef KTHREAD_INSTRUMENT
    long et0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    {
      // Tag with the target T index (T.size() BEFORE enterT increments
      // it — that's the atT for appendEnd).  Lets us correlate pair
      // construction with the T slot it will occupy.
      OpTag _r(thread_id, "enterT", P->p, strat->T.size(), 0);
      enterT(*P, strat);
    }
    audit_T_pLength(ctx, "ps-phase0-post-enterT", thread_id);

    // L-scan after enterT, under L-lock, to avoid races with
    // concurrent enterpairs writing to strat->L.
    //
    // pair.i_r2 is an R-SLOT INDEX.  Primary check is R[i_r2]->p.
    // T-array lookup is reported as a secondary, COSMETIC observation
    // (off-by-one on T side is normal when posInT shifted T entries).
    if (g_debug_ring_enabled.load(std::memory_order_relaxed))
    {
      kt_L_lock(ctx, thread_id);
      int inconsistent_R = 0;   // REAL bug: R[i_r2]->p != pair.p2
      int inconsistent_T = 0;   // cosmetic: T[i_r2].p != pair.p2
      int offby1_T = 0;         // cosmetic: T[i_r2-1].p == pair.p2
      int fp_mismatches = 0;
      int first_bad_r_i_r2 = -1;
      poly first_bad_r_p2 = NULL;
      poly first_bad_r_r_p = NULL;
      int total_pairs = 0;
      for (auto it = strat->L.begin(); it != strat->L.end(); ++it)
      {
        total_pairs++;
        if (it->i_r2 < 0 || it->i_r2 >= (int)strat->T.size()) continue;
        poly t_p = strat->T[it->i_r2].p;
        TObject *r_entry = strat->R[it->i_r2];
        poly r_p = (r_entry != NULL) ? r_entry->p : NULL;
        if (it->dbg_fp != 0)
        {
          unsigned long now_fp = kt_debug_pair_fp((void*)it->p1,
                                                  (void*)it->p2,
                                                  it->i_r1, it->i_r2);
          if (now_fp != it->dbg_fp) fp_mismatches++;
        }
        if (it->p2 != NULL && r_p != NULL && it->p2 != r_p)
        {
          inconsistent_R++;
          if (first_bad_r_i_r2 < 0)
          {
            first_bad_r_i_r2 = it->i_r2;
            first_bad_r_p2 = it->p2;
            first_bad_r_r_p = r_p;
          }
        }
        if (it->p2 != NULL && it->p2 != t_p)
        {
          inconsistent_T++;
          if (it->i_r2 > 0 && it->p2 == strat->T[it->i_r2 - 1].p)
            offby1_T++;
        }
      }
      pthread_mutex_unlock(&ctx->L_lock);
      // Fire a distinct first-occurrence log line for REAL bugs
      // (R-inconsistency) versus cosmetic T-shift observations.
      static std::atomic<int> first_R_seen{-1};
      if (inconsistent_R > 0 && first_R_seen.exchange(inconsistent_R) < 0)
      {
        FILE *log = g_audit_log ? g_audit_log : stderr;
        fprintf(log,
                "\n=== FIRST L-SCAN R-INCONSISTENCY (REAL BUG, "
                "L-lock held): at enterT atT=%d, |L|=%d, "
                "inconsistent_R=%d  inconsistent_T=%d (cosmetic) "
                "offby1_T=%d (cosmetic) fp_mismatches=%d, "
                "first_bad: i_r2=%d pair.p2=%p R[i_r2]->p=%p ===\n",
                (int)strat->T.size()-1, total_pairs,
                inconsistent_R, inconsistent_T, offby1_T, fp_mismatches,
                first_bad_r_i_r2, (void*)first_bad_r_p2,
                (void*)first_bad_r_r_p);
        fflush(log);
      }
      kt_debug_tag("L-scan:after-enterT", NULL, inconsistent_R,
                   inconsistent_T);
      if (fp_mismatches > 0)
        kt_debug_tag("L-scan:fp_MISMATCH", NULL,
                     fp_mismatches, inconsistent_R);
    }
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx)) enterT_accum += kt_now_ns() - et0;
    long es0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    my_arrival = strat->arrival_counter.load(std::memory_order_relaxed);
    strat->enterS(*P, strat, strat->T.size()-1, strat->S.end());
    did_enterS = true;
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx)) enterS_accum += kt_now_ns() - es0;
#endif
  }

#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx)) phase0_work = kt_now_ns() - p0_t0;
#endif

  // ------------------------------------------------------------------
  // Phase 1 — S-shared + L-exclusive: thread-local B, arrival_id
  // filter, concurrent peer drainers possible.  enterpairs walks S
  // under shared lock (tombstone-only writes, arrival-id filter skips
  // h and peer survivors); enterOnePair writes to thread-local B;
  // chainCritNormal merges thread-local B into strat->L under L-lock
  // (the L-lock is acquired here and released on scope exit, so
  // chainCritNormal's L mutations are serialised against peer
  // drainers' own chainCritNormal and against fill_active_slots).
  //
  // Downgrade S-exclusive → S-shared: no atomic path, so unlock then
  // rdlock.  Peer drainers' phase 0 (which needs exclusive) could
  // steal in between; that is fine — their phase 0 runs, they finish
  // enterS, they drop exclusive, we get shared.
  // ------------------------------------------------------------------
  strat->S.unlock_exclusive();

  if (did_enterS)
  {
    // Stack-allocate thread-local B and pairtest-hit vector.  LSet's
    // comparator indirects through strat->compareL (via its
    // CompareLObject::strat pointer), so we must set strat on the
    // local B exactly as skStrategy::skStrategy() does for strat->B.
    LSet local_B;
    local_B.key_comp().strat = strat;
    std::vector<SElement*> local_pairtest_hits;

    // Save any previous thread-local overrides (NULL in practice; the
    // drain does not recurse through enterpairs, but defence in depth).
    LSet* saved_B_override = t_local_B_override;
    uint64_t saved_my_arrival = t_local_my_arrival;
    std::vector<SElement*>* saved_pairtest_hits = t_local_pairtest_hits;

    t_local_B_override = &local_B;
    t_local_my_arrival = my_arrival;
    t_local_pairtest_hits = &local_pairtest_hits;

    kt_S_lock_shared(ctx, thread_id);
    kt_L_lock_phase2(ctx, thread_id);

#ifdef KTHREAD_INSTRUMENT
    long p1_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
    long p1_l_t0 = p1_t0;  // L-lock held from here until unlock below
    long ep0 = KT_STATS(ctx) ? kt_now_ns() : 0;
    // Sample how many drainers are concurrently in phase 1 right now.
    // enterpairs_active counts drain_survivor_queue callers; this is an
    // upper bound on phase-1 concurrency (each caller passes through
    // phase 0 then phase 1).  Not all of these are necessarily in
    // phase 1 at the same instant, but the sample gives a reasonable
    // "peer count" indicator.
    if (KT_STATS(ctx))
    {
      long n = ctx->enterpairs_active.load(std::memory_order_relaxed);
      ThreadStats &ts = KT_TS(ctx, thread_id);
      if (n > ts.phase1_concurrent_max) ts.phase1_concurrent_max = n;
    }
#endif

    // pos_it captured in phase 0 remains valid because (a) SORDER_STANDARD
    // insert only happens in phase 0 which needs S-exclusive, peer drainers
    // are blocked while I hold S-shared; (b) tombstone CAS erases don't
    // shift indices.  pos_it now points at h (or the entry at find_pos'd
    // index); the arrival_id filter skips h in the clearS walk.
    // Hypothesis: `strat->T.size()-1` can be stale by the time we
    // reach here because we released S-exclusive and reacquired
    // S-shared between our own enterT and this enterpairs call.
    // A peer drainer can enterT in that window, incrementing
    // T.size().  Log atR and our own i_r so we can detect the
    // divergence.
    //
    // FIX (rr-replay watchpoint, 20Apr2026): use P->i_r (set under
    // our exclusive lock by enterT at kutil.cc:9092) instead of the
    // racy strat->T.size()-1.  Watchpoint proof: pair.i_r2=12 but
    // matching T entry was at T[5] (.i_r=5) — off by 7 = #peer
    // enterTs in the unlock-shared window.
    int atR_for_pairs = P->i_r;
    kt_debug_tag("enterpairs:atR_vs_P.i_r",
                 (void*)P->p, atR_for_pairs, P->i_r);
    // ALSO check: does R[atR_for_pairs]->p equal P->p at this moment?
    // atR is an R-slot index; T-based check is cosmetic (T shifts).
    TObject *r_entry = (atR_for_pairs >= 0 && atR_for_pairs < (int)strat->T.size())
                       ? strat->R[atR_for_pairs] : NULL;
    poly R_atR_p = (r_entry != NULL) ? r_entry->p : NULL;
    poly T_atR_p = strat->T[atR_for_pairs].p;
    if (R_atR_p != NULL && R_atR_p != P->p)
    {
      kt_debug_tag("enterpairs:R[atR]->p!=P->p (REAL)",
                   (void*)P->p, atR_for_pairs, P->i_r);
    }
    if (T_atR_p != P->p && R_atR_p == P->p)
    {
      kt_debug_tag("enterpairs:T[atR]!=P->p (cosmetic)",
                   (void*)P->p, atR_for_pairs, P->i_r);
    }
    if (rField_is_Ring(currRing))
      superenterpairs(P->p, strat->S.size()-1, P->ecart, pos_it, strat, atR_for_pairs);
    else
      enterpairs(P->p, strat->S.size()-1, P->ecart, pos_it, strat, atR_for_pairs);

#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx))
    {
      long now = kt_now_ns();
      enterpairs_accum += now - ep0;
      phase1_work = now - p1_t0;
      long p1_l_dt = now - p1_l_t0;
      ThreadStats &ts = KT_TS(ctx, thread_id);
      ts.phase1_l_ns += p1_l_dt;
    }
#endif

    // Safety: clear any residual local B / pairtest entries before
    // the LSet / vector destructors run (normally they're empty, but
    // a control-flow short-circuit could leave stragglers).
    //
    // Restore previous overrides BEFORE the LSet goes out of scope,
    // so nothing can reference it via the thread-local slot.
    t_local_B_override = saved_B_override;
    t_local_my_arrival = saved_my_arrival;
    t_local_pairtest_hits = saved_pairtest_hits;

    pthread_mutex_unlock(&ctx->L_lock);
    strat->S.unlock_shared();
  }
  else
  {
    // No enterS happened (IDLIFT/syzComp gated it); nothing to do in
    // phase 1.  The exclusive lock was already released above.
  }

  kDeleteLcm(P);
  ctx->stat_survivors.fetch_add(1, std::memory_order_relaxed);

  P->Init();

#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx))
  {
    long ps_total = kt_now_ns() - ps_t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.ps_redtail_ns += redtail_accum;
    ts.ps_enterT_ns += enterT_accum;
    ts.ps_enterpairs_ns += enterpairs_accum;
    ts.ps_enterS_ns += enterS_accum;
    long accounted = redtail_accum + enterT_accum + enterpairs_accum + enterS_accum;
    ts.ps_other_ns += (ps_total - accounted);
    ts.drain_survivors++;

    ts.phase0_ns += phase0_work;
    ts.phase1_ns += phase1_work;
    // phase1_l_ns: time with L-lock held during phase 1 (formerly
    // "phase2").  Accumulated via the lambda below — see the
    // in-phase-1 block above.
    ts.phase_survivors++;
  }
#endif
}

/* ------------------------------------------------------------------ */
/*  Queue a survivor from an ActivePoly slot.                          */
/*                                                                     */
/*  Transfers ownership of ap->P into the survivor_queue and clears    */
/*  the slot. Updates max-depth statistic.                             */
/* ------------------------------------------------------------------ */

static void queue_survivor(SweepContext *ctx, ActivePoly *ap)
{
  kt_surv_q_lock(ctx, 0);
  ctx->survivor_queue->push_back(ap->P);
  long depth = (long)ctx->survivor_queue->size();
  long prev = ctx->stat_max_queue_depth.load(std::memory_order_relaxed);
  while (depth > prev &&
         !ctx->stat_max_queue_depth.compare_exchange_weak(
             prev, depth, std::memory_order_relaxed))
    /* retry */;
  pthread_mutex_unlock(&ctx->survivor_queue_mutex);

  // Clear slot — ownership has been transferred to the queue.
  ap->P.Init();
  ap->occupied = false;
  ap->is_survivor = false;
}

/* ------------------------------------------------------------------ */
/*  Drain the survivor FIFO: pop one survivor at a time and delegate   */
/*  the enterT/enterpairs/enterS sequence to process_survivor_lobject. */
/*                                                                     */
/*  Multi-drainer safe: any number of threads may call this function   */
/*  concurrently. Each call pops at most as many survivors as remain   */
/*  in the queue at the time of the pop.                               */
/*                                                                     */
/*  Task 506: locking moved into process_survivor_lobject, which runs  */
/*  explicit phase 0 / 1 / 2 blocks.  drain_survivor_queue no longer   */
/*  holds S.lock() or L_lock around the call.                          */
/*                                                                     */
/*  Each drain worker increments ctx->enterpairs_active on entry and   */
/*  decrements on exit. The counter is polled at termination time to   */
/*  tell when all drain workers are quiescent.                         */
/* ------------------------------------------------------------------ */

static void drain_survivor_queue(SweepContext *ctx, int thread_id)
{
  OpTag _drain_tag(thread_id, "drain_survivor_queue");
#ifdef KTHREAD_INSTRUMENT
  long drain_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
  long drained_this_call = 0;
#endif

  audit_T_pLength(ctx, "drain-entry", thread_id);
  ctx->enterpairs_active.fetch_add(1, std::memory_order_acq_rel);

  while (true)
  {
    LObject P;
    kt_surv_q_lock(ctx, thread_id);
    if (ctx->survivor_queue->empty())
    {
      pthread_mutex_unlock(&ctx->survivor_queue_mutex);
      break;
    }
    P = ctx->survivor_queue->front();
    ctx->survivor_queue->pop_front();
    pthread_mutex_unlock(&ctx->survivor_queue_mutex);

    // process_survivor_lobject now manages its own locking via the
    // phase-0/1/2 structure (task 506).  It takes S-exclusive in
    // phase 0 and releases before returning.
    process_survivor_lobject(ctx, &P, thread_id);

#ifdef KTHREAD_INSTRUMENT
    drained_this_call++;
#endif

    // Wake anyone blocked waiting for L to refill.
    kt_L_lock(ctx, thread_id);
    pthread_cond_broadcast(&ctx->pairs_available);
    pthread_mutex_unlock(&ctx->L_lock);
  }

  ctx->enterpairs_active.fetch_sub(1, std::memory_order_acq_rel);

#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx) && drained_this_call > 0)
  {
    long dt = kt_now_ns() - drain_t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.drain_ns += dt;
    ts.drain_count++;
  }
#endif
}

/* (batch_reduce removed — replaced by cooperative sweep_phase) */

/* ------------------------------------------------------------------ */
/*  Worker thread function                                             */
/* ------------------------------------------------------------------ */

struct WorkerArg
{
  SweepContext *ctx;
  int thread_id;
};

/**
 * Worker thread (milestone d, task 283). No barriers. The worker is
 * a continuous tile-pull loop. It waits on tiles_avail_cv when the
 * tile cursor catches up to tile_end, and exits when ctx->done is
 * set and no more tiles remain.
 */
static void *worker_thread(void *arg)
{
  WorkerArg *wa = (WorkerArg *)arg;
  SweepContext *ctx = wa->ctx;
  int thread_id = wa->thread_id;
  free(wa);

  kt_debug_tid = thread_id;
  currRing = ctx->r;
  si_opt_1 = ctx->saved_si_opt_1;
  si_opt_2 = ctx->saved_si_opt_2;

  // Wait for all threads to be created and main to be ready
  pthread_barrier_wait(&ctx->startup_barrier);

#ifdef KTHREAD_INSTRUMENT
  long sw_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
  tile_pull_loop(ctx, thread_id, /*block=*/true);
#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx))
  {
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.sweep_ns += kt_now_ns() - sw_t0;
    ts.sweep_count++;
  }
#endif
  return NULL;
}

/* ------------------------------------------------------------------ */
/*  Main parallel loop                                                 */
/* ------------------------------------------------------------------ */

void bba_parallel_loop(SweepContext *ctx)
{
  kt_debug_tid = 0;  // main thread
  kStrategy strat = ctx->strat;
  int nthreads = ctx->num_threads;

  // Two layers of debug enablement:
  //   SINGULAR_AUDIT_T=1     -> audit-walks on top of the ring
  //   SINGULAR_DEBUG_RING=1  -> event ring + OpTag breadcrumbs + dError hook
  // SINGULAR_AUDIT_T implies SINGULAR_DEBUG_RING (the audit needs the ring).
  bool want_audit = (getenv("SINGULAR_AUDIT_T") != NULL);
  bool want_ring  = want_audit || (getenv("SINGULAR_DEBUG_RING") != NULL);
  g_audit_T_enabled.store(want_audit, std::memory_order_relaxed);
  g_debug_ring_enabled.store(want_ring, std::memory_order_relaxed);
  g_audit_T_checked.store(0, std::memory_order_relaxed);
  g_audit_hit_count.store(0, std::memory_order_relaxed);
  g_plen_seq.store(0, std::memory_order_relaxed);
  g_plen_frozen.store(false, std::memory_order_relaxed);
  for (int t = 0; t < MAX_AUDIT_THREADS; t++)
    g_thread_debug[t] = ThreadDebugState{};
  for (int i = 0; i < PLEN_RING_SIZE; i++)
    g_plen_ring[i].seq.store(UINT64_MAX, std::memory_order_relaxed);
  for (int i = 0; i < TNODE_REGISTRY_SIZE; i++)
    g_tnode_registry[i].addr.store(NULL, std::memory_order_relaxed);
  for (int i = 0; i < TNODE_HEAD_SNAPSHOT_SIZE; i++)
    g_tnode_head_snapshot[i].store(NULL, std::memory_order_relaxed);
  if (want_ring && g_audit_log == NULL)
  {
    const char *path = getenv("SINGULAR_AUDIT_LOG");
    if (path == NULL) path = "/tmp/audit-run/audit-debug.log";
    g_audit_log = fopen(path, "w");
    if (g_audit_log != NULL) setvbuf(g_audit_log, NULL, _IONBF, 0);
    else g_audit_log = stderr;
  }
  // Install the dErrorBreak hook so dReportError sites (bucket
  // length, etc.) get a ring dump before aborting.  Controlled at
  // fire-time by SINGULAR_ABORT_ON_DERROR=1 (checked in dErrorBreak).
  // Also install the kBucketInit tag hook so every call site shows
  // up in the ring.
  if (want_ring)
  {
    dErrorBreak_hook = dump_ring_generic;
    kbucket_debug_tag = kbucket_debug_tag_adapter;
    kbucket_debug_check_tnode = kt_debug_check_tnode_write;
  }

  // Startup-introduced-inconsistency probe: before the parallel phase
  // begins, scan strat->L and check each pair's (p2, i_r2) against
  // the current T array.  If ANY pair is inconsistent at this point,
  // the bug originates in startup code (reorderT / updateT are the
  // live suspects).  Also write a summary directly to g_audit_log
  // since these events would otherwise be evicted from the ring
  // long before any fault dump.
  if (want_ring)
  {
    int inconsistent_T = 0, inconsistent_R = 0, offby1 = 0, checked = 0;
    int r_reports = 0;
    FILE *log = g_audit_log ? g_audit_log : stderr;
    fprintf(log, "\n=== STARTUP L-scan: T.size=%d, L.size=%d ===\n",
            (int)strat->T.size(), (int)strat->L.size());
    for (auto it = strat->L.begin(); it != strat->L.end(); ++it)
    {
      checked++;
      if (it->i_r2 < 0 || it->i_r2 >= strat->T.size()) continue;
      poly t_p = strat->T[it->i_r2].p;
      poly r_p = (strat->R[it->i_r2] != NULL) ? strat->R[it->i_r2]->p : NULL;
      bool t_mismatch = (it->p2 != NULL && it->p2 != t_p);
      bool r_mismatch = (it->p2 != NULL && r_p != NULL && it->p2 != r_p);
      if (t_mismatch) inconsistent_T++;
      if (r_mismatch) inconsistent_R++;
      if (t_mismatch && it->i_r2 > 0
          && it->p2 == strat->T[it->i_r2 - 1].p)
        offby1++;
      // Print per-pair detail only for R-mismatch (REAL BUG).  T-only
      // mismatches are cosmetic (posInT shifted T after pair creation;
      // pair.i_r2 is an R-slot index, stable across shifts).
      if (r_mismatch && r_reports < 20)
      {
        fprintf(log,
                "  [REAL BUG] pair L[?]: p2=%p i_r2=%d  "
                "R[i_r2]->p=%p  T[i_r2].p=%p\n",
                (void*)it->p2, it->i_r2, (void*)r_p, (void*)t_p);
        r_reports++;
      }
    }
    fprintf(log,
            "STARTUP L-scan: checked=%d  "
            "inconsistent_R=%d (REAL)  "
            "inconsistent_T=%d (cosmetic)  "
            "offby1_T=%d (cosmetic)\n\n",
            checked, inconsistent_R, inconsistent_T, offby1);
    fflush(log);
    kt_debug_tag("STARTUP:L-scan-complete",
                 NULL, checked, inconsistent_R);
  }

  ctx->saved_posInT = strat->posInT;
  strat->posInT = posInT_appendEnd;

  SI_SAVE_OPT1(ctx->saved_opt1);
  si_opt_1 |= Sy_bit(OPT_REDTHROUGH);

  // Pre-expand tailRing
  while (strat->tailRing != currRing)
  {
    if (!kStratChangeTailRing(strat))
      break;
  }

  ctx->withT = !strat->homog;
  ctx->saved_si_opt_1 = si_opt_1;
  ctx->saved_si_opt_2 = si_opt_2 & ~Sy_bit(OPT_PROT);

  ctx->done.store(false, std::memory_order_release);

  // Task 280 milestone (a): reserve T capacity up front.
  //
  // T is BlockArray<TObject>. Individual element addresses are stable
  // across push_back (existing blocks are never moved), but the
  // directory of block pointers (`blocks[]`) is grown by realloc — via
  // free()+calloc()+memcpy in BlockArray::ensure_capacity. A reader
  // dereferencing the directory pointer (`blocks[i>>BLOCK_SHIFT][...]`)
  // concurrently with that realloc would observe a freed pointer.
  //
  // Later milestones let workers scan T concurrently with enterT's
  // push_back, so the directory must not be reallocated during the
  // parallel window. Reserve directory capacity now for 4x the current
  // T length (with a floor) so enterT during the parallel phase only
  // appends into existing directory slots.
  //
  // The BlockArray block size is 1024 (BLOCK_SHIFT=10); reserving
  // capacity for N elements allocates ceil(N/1024) blocks and sizes
  // the directory accordingly.
  // Task 512 worker-side-drain: BlockArray::ensure_capacity aborts if
  // it needs to grow the directory after freeze_dir() has been called
  // (the directory-realloc path is unsafe under concurrent readers).
  // Pre-allocate a generous reserve here while still single-threaded,
  // then freeze.  Staging-redsb workloads with THREADS=4 grow T to a
  // few thousand entries; 1M covers orders of magnitude more than any
  // observed test case.  If exceeded, the abort() fires at the next
  // enterT with a clear diagnostic, at which point we bump the reserve.
  //
  // strat->R is pre-allocated too (written in enterT; same race class).
  // strat->sevT same.
  {
    int cur_T = strat->T.size();
    long reserve_n = (long)cur_T * 16;
    if (reserve_n < 1048576) reserve_n = 1048576;
    strat->T.ensure_capacity((int)reserve_n);
    strat->sevT.ensure_capacity((int)reserve_n);
    strat->R.ensure_capacity((int)reserve_n);
    strat->T.freeze_dir();
    strat->sevT.freeze_dir();
    strat->R.freeze_dir();
  }

#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx))
  {
    ctx->start_ns = kt_now_ns();
    int tt = ctx->num_workers + 1;
    for (int i = 0; i < tt; i++)
      memset(&ctx->tstats[i], 0, sizeof(ThreadStats));
    ctx->rounds->clear();
  }
#endif

  for (int t = 0; t < ctx->num_workers; t++)
  {
    WorkerArg *wa = (WorkerArg *)malloc(sizeof(WorkerArg));
    wa->ctx = ctx;
    wa->thread_id = t + 1;
    ctx->thread_ids[t] = t + 1;
    pthread_create(&ctx->threads[t], NULL, worker_thread, wa);
  }

  // Task 511 t-iterator-migrate-hot-path: the startup-time pLength
  // pre-population loop is gone.  Every T entry in strat->T was
  // created by enterT (task 510), which computes pLength inline
  // before release-publishing the slot.  So by the time
  // bba_parallel_loop runs, every T[j] already has pLength > 0 and
  // published = true — no pre-scan needed.

  // Wait for all workers to start before entering main loop
  if (ctx->num_workers > 0)
    pthread_barrier_wait(&ctx->startup_barrier);

  // Continuous refill loop (milestone d, task 283). No barriers.
  //
  // Each iteration:
  //   1. If siCntrlc: shutdown path.
  //   2. Drain any survivors from the FIFO.  Main still takes a drain
  //      pass here (important when workers are all busy on tiles); in
  //      addition, workers now drain from their tile-idle branch
  //      (task 512 worker-side-drain).
  //      enterT (invoked from drain_survivor_queue) computes pLength
  //      inline and release-publishes T slots, so no post-drain
  //      refresh is needed (task 511).
  //   3. Refill empty slots from L and (re)publish tiles for any slot
  //      that closer-reduce marked needs_republish.
  //   4. If nothing to do (all slots empty, L empty, queue empty):
  //        a. Help sweep any outstanding tiles (non-blocking) so we
  //           don't sit on main while workers still have work.
  //        b. If still idle after that, break out — we are done.
  //   5. Else, wait briefly on slot_freed_cv so closers can wake us.
  //
  // Task 512 worker-side-drain lifted the main-only-drain architectural
  // invariant.  Workers in tile_pull_loop, when tile_cursor catches up
  // to tile_end, peek at the survivor queue and call drain_survivor_queue
  // if non-empty.  process_survivor_lobject's phase 0 (S-exclusive)
  // serialises peer drainers; phase 1 (S-shared + L-exclusive +
  // thread-local B) allows concurrent peer drainers to run enterpairs
  // in parallel.  Main's own drain pass remains for responsiveness
  // when all workers are busy on tiles (or, at the end of a run, when
  // workers are parked on tiles_avail_cv waiting for tiles but the
  // survivor queue is non-empty).
  //
  // T/L access audit for main outside drain (task 512):
  //   - pop_and_prepare: reads strat->L.{top,pop,empty,size} under
  //     L_lock (refill_and_publish acquires it); safe vs worker
  //     chainCritNormal which takes L_lock inside phase 1.
  //   - pop_and_prepare also reads strat->T.size()-1 to snapshot
  //     ap->sl_snapshot; this is a plain int load, monotonic (enterT
  //     only grows it).  Worst case the snapshot misses T entries a
  //     peer drainer is publishing right now; those entries show up
  //     in the next refill pass.  No lock needed.
  //   - refill_and_publish's republish path re-snapshots T.size()-1;
  //     same reasoning.
  //   - Termination check reads strat->L.empty() without L_lock.
  //     Safe iff no worker is actively modifying L.  Workers modify L
  //     only inside chainCritNormal (under L_lock) reached from
  //     process_survivor_lobject's phase 1.  Phase 1 entry requires
  //     a survivor in the queue.  If our check finds queue_empty AND
  //     no slot FILLED, no closer is running, so no survivor will be
  //     queued, so no worker will enter drain.  The L.empty() read is
  //     safe under this precondition.  (The shutdown paths at
  //     siCntrlc / overflow ARE racy: handled separately below.)
  //   - strat->T[i] accesses on main outside drain: none.
  //     sweep_one_tile runs through acquire-loaded published gate
  //     (task 511) and is called from tile_pull_loop which main
  //     invokes with block=false — main is a tile reader at the same
  //     terms as workers.
  //   - strat->B: accessed via strat_B(strat) which returns the
  //     thread-local override inside phase 1, and strat->B elsewhere.
  //     Main outside drain never touches strat->B (only serial code
  //     before / after bba_parallel_loop does).
  // Error-path L-clear flags (set here, acted on after workers join).
  // Task 512 worker-side-drain: we cannot touch strat->L while workers
  // may be inside chainCritNormal (which takes L_lock from phase 1).
  // Defer L mutation until after ctx->done + pthread_join below.
  bool clear_L_after_join = false;
  bool clear_slots_after_join = false;

  while (true)
  {
    if (siCntrlc)
    {
      clear_L_after_join = true;
      strat->noClearS = TRUE;
      goto parallel_shutdown;
    }

    if (strat->overflow || errorreported)
    {
      clear_L_after_join = true;
      clear_slots_after_join = true;
      goto parallel_shutdown;
    }

    // Step 1: drain survivor FIFO on main thread.
    //
    // Task 511 t-iterator-migrate-hot-path: the post-drain pLength
    // refresh loop that ran here after every drain is gone — enterT
    // now computes pLength inline and release-publishes the slot in
    // one step (see kutil.cc enterT, task 510).  Every T entry
    // created during drain is therefore already pLength-filled and
    // published by the time drain_survivor_queue returns; no external
    // refresh is needed.
    {
      kt_surv_q_lock(ctx, 0);
      bool queue_empty = ctx->survivor_queue->empty();
      pthread_mutex_unlock(&ctx->survivor_queue_mutex);
      if (!queue_empty)
      {
        drain_survivor_queue(ctx, 0);
      }
    }

    // Step 2: refill empty slots and republish any slot needing another pass.
    int in_pipeline = refill_and_publish(ctx);

    // Step 3: termination check.
    if (in_pipeline == 0)
    {
      kt_surv_q_lock(ctx, 0);
      bool queue_empty = ctx->survivor_queue->empty();
      pthread_mutex_unlock(&ctx->survivor_queue_mutex);
      // Task 512 worker-side-drain: L.empty() must be read under L_lock
      // to avoid racing with worker drainers' chainCritNormal pushes.
      kt_L_lock(ctx, 0);
      bool L_empty = strat->L.empty();
      pthread_mutex_unlock(&ctx->L_lock);
      // FIX (rr+gdb MCP, 20 Apr 2026 post-2200 investigation): worker
      // drainers may be mid-flight between popping a survivor from
      // the queue (queue becomes empty) and calling enterpairs (L
      // gains new pairs).  Without checking enterpairs_active, main
      // sees queue_empty && L_empty and breaks prematurely, losing
      // the pairs the worker is about to add.  Repro: /tmp/module3
      // at THREADS=4 gave 1 instead of 3 ~40% of runs; with this
      // check, 0 wrong in initial tests.
      int drainers_active = ctx->enterpairs_active.load(std::memory_order_acquire);
      if (queue_empty && L_empty && drainers_active == 0)
        break;
      // Else there is still work (drain produced survivors, or L has
      // new entries, or a worker drainer hasn't finished adding
      // pairs) — loop back immediately.
      continue;
    }

    // Step 4: help sweep with tiles while waiting. Non-blocking: main
    // returns as soon as the tile cursor catches up with tile_end so
    // it can go back to refilling/draining.
    tile_pull_loop(ctx, 0, /*block=*/false);

    // Step 5: wait for a slot-freed event or new tiles available.
    // We use a short timed wait so we don't rely solely on signals —
    // if a closer races our check, we won't deadlock.
    pthread_mutex_lock(&ctx->publish_lock);
    // Cheap re-check: if any slot is empty or needs republish, or the
    // survivor queue got work, don't sleep.
    bool any_free = false;
    for (int s = 0; s < ctx->max_active; s++)
    {
      int st = ctx->active[s].state.load(std::memory_order_acquire);
      if (st == SLOT_EMPTY ||
          ctx->active[s].needs_republish.load(std::memory_order_acquire))
      {
        any_free = true;
        break;
      }
    }
    if (!any_free)
    {
      struct timespec ts;
      clock_gettime(CLOCK_REALTIME, &ts);
      // 1 ms timeout.
      ts.tv_nsec += 1000000;
      if (ts.tv_nsec >= 1000000000) { ts.tv_sec++; ts.tv_nsec -= 1000000000; }
      pthread_cond_timedwait(&ctx->slot_freed_cv, &ctx->publish_lock, &ts);
    }
    pthread_mutex_unlock(&ctx->publish_lock);
  }

parallel_shutdown:
  // Signal done and wake all workers so they exit their tile_pull_loop.
  ctx->done.store(true, std::memory_order_release);
  pthread_mutex_lock(&ctx->publish_lock);
  pthread_cond_broadcast(&ctx->tiles_avail_cv);
  pthread_mutex_unlock(&ctx->publish_lock);

  for (int t = 0; t < ctx->num_workers; t++)
    pthread_join(ctx->threads[t], NULL);

  // Workers joined — unfreeze the BlockArray directories so downstream
  // serial code can grow T / sevT / R normally.
  strat->T.unfreeze_dir();
  strat->sevT.unfreeze_dir();
  strat->R.unfreeze_dir();

  // Now workers are joined — safe to mutate strat->L and slots.
  // (Task 512 worker-side-drain: chainCritNormal from worker drain is
  // the only writer of L outside main; joined workers cannot be inside
  // that critical section any longer.)
  if (clear_L_after_join)
  {
    while (!strat->L.empty()) strat->L.pop_and_erase();
  }
  if (clear_slots_after_join)
  {
    for (int i = 0; i < ctx->max_active; i++)
      ctx->active[i].occupied = false;
  }

  strat->posInT = ctx->saved_posInT;
  SI_RESTORE_OPT1(ctx->saved_opt1);

  if (TEST_OPT_PROT)
  {
    Print("\n// parallel bba: polys=%ld, reductions=%ld, zeros=%ld, "
          "survivors=%ld, max_queue_depth=%ld, threads=%d (workers=%d)\n",
          ctx->stat_rounds.load(),
          ctx->stat_reductions.load(),
          ctx->stat_zeros.load(),
          ctx->stat_survivors.load(),
          ctx->stat_max_queue_depth.load(),
          nthreads, ctx->num_workers);
  }

  // Emit diagnostic queue depth to stderr if requested.
  if (getenv("SINGULAR_DRAIN_DEPTH_LOG") != NULL)
  {
    fprintf(stderr,
            "[parallel-bba] survivors=%ld max_queue_depth=%ld rounds=%ld\n",
            ctx->stat_survivors.load(),
            ctx->stat_max_queue_depth.load(),
            ctx->stat_rounds.load());
  }

#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx))
  {
    fflush(stdout);
    long total_ns = kt_now_ns() - ctx->start_ns;
    int tt = ctx->num_workers + 1;
    fprintf(stderr, "\n============================================================\n");
    fprintf(stderr, "[kthread-stats] tag=%s threads=%d workers=%d wall_ns=%ld (%.3fs)\n",
            ctx->workload_tag, ctx->num_threads, ctx->num_workers,
            total_ns, total_ns / 1e9);
    fprintf(stderr, "[kthread-stats] reductions=%ld survivors=%ld rounds=%ld max_qd=%ld\n",
            ctx->stat_reductions.load(), ctx->stat_survivors.load(),
            ctx->stat_rounds.load(), ctx->stat_max_queue_depth.load());
    fprintf(stderr, "[kthread-stats] %-4s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
            "tid", "sweep_ns", "reduce_ns", "drain_ns", "B0_wait_ns", "B1_wait_ns",
            "L_wait_ns", "sweeps", "reduces", "drains", "survivors");
    for (int i = 0; i < tt; i++)
    {
      ThreadStats &ts = ctx->tstats[i];
      fprintf(stderr, "[kthread-stats] %-4d %12ld %12ld %12ld %12ld %12ld %12ld %12ld %12ld %12ld %12ld\n",
              i, ts.sweep_ns, ts.reduce_ns, ts.drain_ns,
              ts.wait_B0_ns, ts.wait_B1_ns, ts.L_lock_wait_ns,
              ts.sweep_count, ts.reduce_count, ts.drain_count,
              ts.drain_survivors);
    }
    // process_survivor breakdown (main + drain workers)
    long tot_rt = 0, tot_et = 0, tot_ep = 0, tot_es = 0, tot_oth = 0, tot_drain = 0;
    for (int i = 0; i < tt; i++)
    {
      ThreadStats &ts = ctx->tstats[i];
      tot_rt += ts.ps_redtail_ns;
      tot_et += ts.ps_enterT_ns;
      tot_ep += ts.ps_enterpairs_ns;
      tot_es += ts.ps_enterS_ns;
      tot_oth += ts.ps_other_ns;
      tot_drain += ts.drain_ns;
    }
    fprintf(stderr, "[kthread-stats] process_survivor: redtail=%ld enterT=%ld enterpairs=%ld enterS=%ld other=%ld (sum_drain=%ld)\n",
            tot_rt, tot_et, tot_ep, tot_es, tot_oth, tot_drain);
    fprintf(stderr, "[kthread-stats] trylock: ");
    for (int i = 0; i < tt; i++)
      fprintf(stderr, "t%d=%ld/%ld ", i,
              ctx->tstats[i].enterpairs_trylock_fail,
              ctx->tstats[i].enterpairs_trylock_count);
    fprintf(stderr, "\n");

    // Phase breakdown (task 507 enterpairs-parallel-measure +
    // task 508 enterpairs-parallel-phase1).  phase2_* renamed to
    // phase1_l_* (L-lock wait/held during phase 1, not a separate
    // phase).  ph1_cmax is the peak value of ctx->enterpairs_active
    // sampled by this thread on phase-1 entry (how many drainers
    // were concurrently in phase 1).
    fprintf(stderr, "[kthread-stats] %-4s %14s %14s %14s %14s %14s %14s %10s %10s %10s %10s\n",
            "tid", "phase0_wait", "phase0_ns", "phase1_wait", "phase1_ns",
            "phase1_l_wait", "phase1_l_ns", "ph_surv", "ph_s_cas", "ph_l_cas",
            "ph1_cmax");
    for (int i = 0; i < tt; i++)
    {
      ThreadStats &ts = ctx->tstats[i];
      fprintf(stderr, "[kthread-stats] p%-3d %14ld %14ld %14ld %14ld %14ld %14ld %10ld %10ld %10ld %10ld\n",
              i, ts.phase0_wait_ns, ts.phase0_ns, ts.phase1_wait_ns, ts.phase1_ns,
              ts.phase1_l_wait_ns, ts.phase1_l_ns, ts.phase_survivors,
              ts.phase_s_cas_fail, ts.phase_l_cas_fail,
              ts.phase1_concurrent_max);
    }

    // Worker-side drain participation (task 512 worker-side-drain).
    //   drain_survivors : # of survivors this thread handled
    //                     (was 0 for tid>0 before task 512; post-task,
    //                     non-zero on workloads with survivor-queue
    //                     backpressure).
    //   idle_drain_ns / idle_drain_count : time & count of drain hops
    //                     taken from the tile-idle branch of
    //                     tile_pull_loop.
    fprintf(stderr, "[kthread-stats] %-4s %14s %14s %14s\n",
            "tid", "drain_survivors", "wrk_drain_ns", "wrk_drain_cnt");
    for (int i = 0; i < tt; i++)
    {
      ThreadStats &ts = ctx->tstats[i];
      fprintf(stderr, "[kthread-stats] d%-3d %14ld %14ld %14ld\n",
              i, ts.drain_survivors,
              ts.worker_drain_idle_ns,
              ts.worker_drain_idle_count);
    }

    // CSV output for per-round records
    const char *csv_path = getenv("SINGULAR_KTHREAD_CSV");
    if (csv_path != NULL)
    {
      FILE *f = fopen(csv_path, "w");
      if (f != NULL)
      {
        fprintf(f, "round_id,ts_ns,qd_start,qd_end,active_slots,reductions,"
                   "min_reduce_ns,max_reduce_ns,sum_reduce_ns,sweep_ns_main,round_total_ns\n");
        for (auto &rr : *ctx->rounds)
        {
          fprintf(f, "%ld,%ld,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%ld\n",
                  rr.round_id, rr.timestamp_ns, rr.queue_depth_start,
                  rr.queue_depth_end, rr.active_slots, rr.reductions,
                  rr.min_reduce_ns, rr.max_reduce_ns, rr.sum_reduce_ns,
                  rr.sweep_ns_main, rr.round_total_ns);
        }
        fclose(f);
        fprintf(stderr, "[kthread-stats] wrote %zu round records to %s\n",
                ctx->rounds->size(), csv_path);
      }
      else
      {
        fprintf(stderr, "[kthread-stats] could not open CSV file %s\n", csv_path);
      }
    }
    fprintf(stderr, "============================================================\n\n");
    fflush(stderr);
  }
#endif
}
