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
#include <cstring>
#include <climits>

// Global ctx pointer used as the parallel-mode runtime signal.
// NULL  ⟺ outside bba_parallel_loop (serial mode); non-NULL ⟺ inside
// the parallel loop.  Read by LSet::erase to dispatch between
// physical (serial) and tombstone-only (parallel rdlock) semantics
// (task 363; commit f392b3c33), so this must be defined and
// maintained regardless of KTHREAD_INSTRUMENT.  Originally introduced
// (task 358; run 571) for gdb-driven kt_dump_stats calls.  Set at
// bba_parallel_loop entry, cleared at parallel_shutdown.
SweepContext *kt_current_ctx = NULL;

/* ------------------------------------------------------------------ */
/*  Instrumentation helpers (task 277; run 482)                       */
/* ------------------------------------------------------------------ */
#ifdef KTHREAD_INSTRUMENT
#  define KT_STATS(ctx)       ((ctx)->stats_enabled)
#  define KT_TS(ctx, tid)     ((ctx)->tstats[tid])
#  define KT_TIME_START(var)  long var = kt_now_ns()
#  define KT_TIME_DELTA(var)  (kt_now_ns() - (var))

// Thread-local thread_id so KTHREAD_INSTRUMENT-guarded sites in other
// translation units (kutil.cc chainCritNormal) can find their
// ThreadStats slot without plumbing a parameter through.  Set on
// process_survivor_lobject entry and main's drain entry; -1 means
// "outside the parallel phase, skip instrumentation".
__thread int kt_my_thread_id = -1;

// L-rwlock acquire helpers, split by call site (task 357 / 360 step 5).
//
// Pre-step-5 (task 360 step 4 and earlier) these took a single
// pthread_mutex_t ctx->L_lock that gave full mutual exclusion on every
// L access.  Step 5 (this commit) replaced the mutex with a
// reader/writer lock that lives on strat->L itself (LSet::rwlock_).
// In step 5, every site that took L_lock now takes the **writer** lock
// — full mutual exclusion is preserved because the underlying chunk is
// still a single, non-thread-safe std::multiset / std::vector pair.
// Step 6 (this commit, task 360) introduces the multi-chunk linked
// list.  The hot-path sites take rdlock instead of wrlock:
//
//   - kt_L_lock_drain  : pthread_cond_broadcast hints, no L mutation;
//                        rdlock suffices.
//   - kt_L_lock_term   : main's L.empty() probe — read-only; rdlock.
//   - kt_L_lock_refill : main's pop_and_prepare (top + pop + chunk
//                        traversal); pop is per-chunk and the chain
//                        walk is read-only against compact, so rdlock.
//   - kt_L_lock_phase2 : worker phase-1 chainCritNormal scan +
//                        erases + kMergeBintoL — chainCritNormal is
//                        a read-only scan + tombstone erases (which
//                        are per-element CAS, no chain mutation),
//                        kMergeBintoL is a CAS-append on the tail
//                        chunk's `next`.  All three of these admit
//                        concurrent readers, so rdlock.
//
// Compact (LSet::compact, called from refill_and_publish) is the
// only writer-lock holder; it takes wrlock via a drop-and-reacquire
// pattern: rdlock-pop-unlock, then wrlock-compact-unlock.
//
// The helper functions retain their existing names and bucketing
// (L_lock_*_wait_ns) so the kthread-stats dump format is preserved.
// What they record is rdlock acquire wait time — typically near-zero
// because rdlock has no contention with other rdlock holders;
// non-zero only when compact (wrlock) is in flight.
static inline void kt_L_lock_drain(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    strat->L.rdlock();
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.L_lock_drain_wait_ns += dt;
    ts.L_lock_drain_count++;
  }
  else
  {
    strat->L.rdlock();
  }
}

static inline void kt_L_lock_term(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    strat->L.rdlock();
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.L_lock_term_wait_ns += dt;
    ts.L_lock_term_count++;
  }
  else
  {
    strat->L.rdlock();
  }
}

static inline void kt_L_lock_refill(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    strat->L.rdlock();
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.L_lock_refill_wait_ns += dt;
    ts.L_lock_refill_count++;
  }
  else
  {
    strat->L.rdlock();
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

/* Phase-0 S shared lock: shared lock for redtailBba/pCleardenom/etc.
 * Same lock as phase-1 shared, but recorded in a separate wait bucket
 * (phase0_shared_wait_ns) so the attribution table can distinguish the
 * pre-enterS shared region from the post-enterS phase-1 shared region. */
static inline void kt_S_lock_shared_phase0(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    strat->S.lock_shared();
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.phase0_shared_wait_ns += dt;
  }
  else
  {
    strat->S.lock_shared();
  }
}

/* L rwlock READER acquire used during worker phase 1 (chainCritNormal
 * scan + erases + kMergeBintoL).
 *
 * Step 6 (task 360): downgraded from wrlock to rdlock.  Concurrent
 * worker phase-1 blocks all run as readers; their B-merge is a
 * CAS-append on the tail chunk's `next` pointer; their tombstone
 * erases are per-element CAS on the LObject->deleted flag.  Compact
 * is the only writer; it runs from main's refill path on a
 * drop-and-reacquire dance (rdlock-pop-release, wrlock-compact-release).
 *
 * Records wait time in phase1_l_wait_ns — typically near-zero now
 * since reader-vs-reader is uncontended; non-zero only when compact
 * (wrlock) is in flight, blocking all readers briefly.
 */
static inline void kt_L_lock_phase2(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    strat->L.rdlock();
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.phase1_l_wait_ns += dt;
  }
  else
  {
    strat->L.rdlock();
  }
}
#else
#  define KT_STATS(ctx)       (false)
#  define KT_TIME_START(var)  ((void)0)
#  define KT_TIME_DELTA(var)  (0L)
static inline void kt_L_lock_drain(SweepContext *ctx, int /*tid*/) { ctx->strat->L.rdlock(); }
static inline void kt_L_lock_term(SweepContext *ctx, int /*tid*/) { ctx->strat->L.rdlock(); }
static inline void kt_L_lock_refill(SweepContext *ctx, int /*tid*/) { ctx->strat->L.rdlock(); }
static inline void kt_surv_q_lock(SweepContext *ctx, int /*tid*/) { pthread_mutex_lock(&ctx->survivor_queue_mutex); }
static inline void kt_record_reduce(SweepContext*, int, long, long) {}
static inline void kt_S_lock_exclusive(SweepContext *ctx, int /*tid*/) { ctx->strat->S.lock_exclusive(); }
static inline void kt_S_lock_shared(SweepContext *ctx, int /*tid*/) { ctx->strat->S.lock_shared(); }
static inline void kt_S_lock_shared_phase0(SweepContext *ctx, int /*tid*/) { ctx->strat->S.lock_shared(); }
static inline void kt_L_lock_phase2(SweepContext *ctx, int /*tid*/) { ctx->strat->L.rdlock(); }
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

  // L_lock removed in step 5: strat->L's own rwlock now covers
  // L data protection (LSet::rwlock_, kutil.h).
  // redtail_lock removed in task 361 (run 581): see kthread.cc:1242
  // — redtailBba / redtailBbaAlsoLC_Z now take per-call out pointers
  // for redTailChange / completeReduce_retry, so there's no shared
  // strat write to serialize.

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

  // Phase-1 enterpairs ordering barrier (default on; opt out with
  // SINGULAR_DISABLE_ENTERPAIRS_BARRIER=1 for ablation).
  // Sync the barrier counter to strat's current arrival_counter so the
  // first survivor to enter S in this dispatch passes immediately.
  pthread_cond_init(&ctx->enterpairs_order_cv, NULL);
  pthread_mutex_init(&ctx->enterpairs_order_lock, NULL);
  ctx->next_enterpairs_arrival_id.store(
      strat->arrival_counter.load(std::memory_order_relaxed),
      std::memory_order_relaxed);
  ctx->serialize_enterpairs =
      (getenv("SINGULAR_DISABLE_ENTERPAIRS_BARRIER") == NULL);

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
  // L_lock removed in step 5; LSet's rwlock self-destructs.
  // redtail_lock removed in task 361 (run 581).
  pthread_mutex_destroy(&ctx->survivor_queue_mutex);
  pthread_mutex_destroy(&ctx->publish_lock);
  pthread_cond_destroy(&ctx->pairs_available);
  pthread_cond_destroy(&ctx->tiles_avail_cv);
  pthread_cond_destroy(&ctx->slot_freed_cv);
  pthread_cond_destroy(&ctx->enterpairs_order_cv);
  pthread_mutex_destroy(&ctx->enterpairs_order_lock);
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
/*  Caller must hold strat->L.rdlock(). Returns true if slot filled.   */
/* ------------------------------------------------------------------ */

/*
 * pop_and_prepare — pop next LObject from strat->L and set up ap for
 * reduction.  Caller holds strat->L's reader lock for the duration of
 * the call (acquired via kt_L_lock_refill or directly via rdlock()).
 * The `sl_snapshot` parameter is the T-size bound captured BEFORE
 * taking the L lock (task 305; run 512); this avoids taking S-shared
 * while holding the L lock, which would invert the drainer's S → L
 * order and risk deadlock.
 *
 * Step 5 (task 360) note: pop_and_prepare runs under the L rdlock,
 * which excludes only compact (the wrlock holder).  Concurrent
 * worker phase-1 chainCritNormal scans / B-merge appends are also
 * rdlock holders and run in parallel with this call.  Any tombstone
 * CAS that races with the begin()/top() filtered-iterator skip is
 * resolved under the chunk's own atomic deleted-flag dance — no
 * additional locking needed here.
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
    // Task 305 (run 512): sl_snapshot_arg was captured in the
    // caller BEFORE L_lock was taken, under S-shared — ensures a
    // consistent view of T.size() without inverting the lock order
    // (phase-1 drainers take S-shared then L-exclusive; main acquires
    // S-shared briefly and releases before L_lock).  The captured
    // value may be slightly stale by the time we use it, but enterT
    // only grows T monotonically and any entries added post-capture
    // are picked up in the next refill pass.
    ap->sl_snapshot = sl_snapshot_arg;

    ap->P.PrepareRed(strat->use_buckets);
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

  for (int t = 0; t < total_threads; t++)
  {
    SweepResult &sr = sweep_result(ctx, t, s);
    if (sr.best_reducer >= 0 && best_reducer < 0)
      best_reducer = sr.best_reducer;
    if (sr.best_good >= 0)
    {
      if (best_good < 0 || sr.best_pLength < best_pLength)
      {
        best_good = sr.best_good;
        best_pLength = sr.best_pLength;
      }
    }
  }

  ctx->active[s].best_reducer = best_reducer;
  ctx->active[s].best_good = best_good;
  ctx->active[s].best_pLength = best_pLength;
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

  // The closer-worker that last wrote ap->P / ap->not_sev finished its
  // writes before storing ap->needs_republish=true (release) at
  // close_slot's republish branch, or before storing ap->state=SLOT_EMPTY
  // (release) at the survivor / zero branches.  refill_and_publish
  // observes those release-stores via acquire-loads on needs_republish /
  // state before calling here.  An explicit acquire fence here makes the
  // cross-variable happens-before edge visible to TSan and to weakly-
  // ordered architectures (no-op on x86 TSO).
  std::atomic_thread_fence(std::memory_order_acquire);

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
  ActivePoly *ap = &ctx->active[s];

  // Unoccupied or already-survivor: should not happen in the continuous
  // design (batches only fire for FILLED slots) but keep as a safety.
  if (!ap->occupied || ap->is_survivor) return;

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
}

/*
 * Sweep one tile (target slot + slice). Pure read-only access to
 * strat->T up to tl_snapshot; writes only to the worker's private
 * SweepResult.
 *
 * Task 304 (run 511): the K-stride access pattern
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
      sr.best_reducer = j;

    int ecart_j = strat->T[j].ecart;
    if (ecart_j <= ap->P.ecart)
    {
      int pLen = strat->T[j].pLength;
      if (pLen <= 0) pLen = 3;
      if (sr.best_good < 0 || pLen < sr.best_pLength)
      {
        sr.best_good = j;
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
  int K = ctx->tiles_K;

  while (true)
  {
    uint64_t end = ctx->tile_end.load(std::memory_order_acquire);
    uint64_t cur = ctx->tile_cursor.load(std::memory_order_acquire);

    if (cur >= end)
    {
      if (!block) return;
      if (ctx->done.load(std::memory_order_acquire)) return;
      // Task 305 (run 512): before going to sleep, check whether
      // the survivor queue has work.  If so, drain instead of cond_wait.
      // Only workers do this (block=true); main has its own drain pass
      // in bba_parallel_loop and calls tile_pull_loop with block=false.
      //
      // Gating: skip the drain hop if the queue is observably empty
      // (snapshot under its own mutex so we don't race with a concurrent
      // push).  The cost of the lock is negligible vs the wait.
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
      {
#ifdef KTHREAD_INSTRUMENT
        long ti_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
        pthread_cond_wait(&ctx->tiles_avail_cv, &ctx->publish_lock);
#ifdef KTHREAD_INSTRUMENT
        if (KT_STATS(ctx))
        {
          ThreadStats &ts = KT_TS(ctx, thread_id);
          ts.tile_idle_ns += kt_now_ns() - ti_t0;
          ts.tile_idle_count++;
        }
#endif
      }
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
  kStrategy strat = ctx->strat;
  ActivePoly *ap = &ctx->active[slot];
  (void)thread_id;  // used only by instrumentation wrapper

  int best = (ap->best_good >= 0) ? ap->best_good : ap->best_reducer;

  if (best < 0)
  {
    // No reducer found by cooperative sweep — survivor
    ap->is_survivor = true;
    return;
  }

  // Apply first reduction from cooperative sweep result
  int ei = strat->T[best].ecart;

  ksReducePoly(&ap->P, strat->T.addr(best),
               strat->kNoetherTail(), NULL, NULL, strat);

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
  int in_pipeline = 0;

  for (int s = 0; s < ctx->max_active; s++)
  {
    ActivePoly *ap = &ctx->active[s];
    int st = ap->state.load(std::memory_order_acquire);

    if (st == SLOT_EMPTY)
    {
      // Task 305 (run 512): capture sl_snapshot under S-shared
      // BEFORE taking L_lock, to avoid inverting the worker drainer's
      // S → L lock order.
      kt_S_lock_shared(ctx, 0);
      int sl_snapshot = ctx->strat->T.size() - 1;
      ctx->strat->S.unlock_shared();

      // Try to fill from L.  Caller-side L-lock window:
      //   - pop_and_prepare physically removes the popped LObject's
      //     tree node (writable_set::erase deletes the LObject*); so
      //     it MUST run under wrlock to avoid UAFs in concurrent
      //     worker filtered_iterator scans (they dereference flat_
      //     pointers into the now-freed tree node).
      //   - After releasing wrlock, if needs_compact() flagged, take
      //     wrlock again for compact (we re-check under wrlock).
      // Step 6 (task 360): pop is the ONE write-lock holder besides
      // compact.  See ~/project/docs/parallel-bba-chunked-lsets.md.
      ctx->strat->L.wrlock();
      BOOLEAN got = pop_and_prepare(ctx, ap, sl_snapshot);
      bool needs = ctx->strat->L.needs_compact();
      if (needs) ctx->strat->L.compact();
      ctx->strat->L.unlock();
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
 * Phase-0 / Phase-1 design (task 301; run 508).
 *
 * Task 299 (run 506) landed infrastructure (rwlock on sBasisSet, atomic CAS on
 * tombstone flags, explicit phase 0/1/2 structure, per-phase
 * instrumentation) but kept phase 1 under an S-exclusive lock because
 * enterpairs writes into the strat-global strat->B.  This task enables
 * concurrent phase-1 drainers under a SHARED S-lock by introducing
 * thread-local overrides (see kutil.h):
 *
 *   t_local_B_override       — thread-local LSetChunk* replacing strat->B
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
  kStrategy strat = ctx->strat;
  BOOLEAN withT = ctx->withT;
  (void)thread_id;

#ifdef KTHREAD_INSTRUMENT
  // Bind the thread-local thread_id so chainCritNormal (kutil.cc) can
  // find its ThreadStats slot.  Restore at function exit so any
  // surrounding code that bumped this is preserved.  See task 358 (run 571).
  int saved_kt_my_thread_id = kt_my_thread_id;
  kt_my_thread_id = thread_id;
  long ps_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
  long redtail_accum = 0;
  long enterT_accum = 0;
  long enterpairs_accum = 0;
  long enterS_accum = 0;
  long phase0_work = 0;        // S-exclusive hold time (enterT+enterS)
  long phase0_shared_work = 0; // S-shared hold time (redtailBba etc.)
  long phase1_work = 0;
#endif

  // ------------------------------------------------------------------
  // Phase 0 (task 359; run 572; task 361 run 581 removed redtail_lock) —
  // split between three regions to shrink the S-exclusive critical
  // section to enterT + enterS only.
  //
  //   (a) Unlocked: GetP / initEcart / PrintS.  Pure LObject mutations.
  //
  //   (b) S-shared: redtailBbaAlsoLC_Z / pCleardenom (or pNorm) /
  //       redtailBba / pCleardenom (INTSTRATEGY post-redtail) /
  //       SetShortExpVector.  These read T/S but never mutate them.
  //       Multiple drainers can run this region concurrently.
  //       redtailBba is passed strat->S.end() instead of a precomputed
  //       find_pos result: in the common withT=true (non-homogeneous)
  //       case redtailBba doesn't read S at all (it uses
  //       kFindDivisibleByInT); in the withT=false (homogeneous) case
  //       kFindDivisibleByInS_T iterates S to end(), with the sev/
  //       pLmCmp prefilter making the divisor search cheap.
  //
  //       Task 361 (run 581) decoupled the two output flags
  //       strat->redTailChange and strat->completeReduce_retry from
  //       the redtailBba* call sites by adding default-nullptr out
  //       parameters; the parallel callers below pass thread-local
  //       pointers, eliminating the cross-thread race that previously
  //       required redtail_lock.  No mutex needed any more.
  //
  //   (c) S-exclusive: enterT + capture my_arrival + enterS.  enterS
  //       returns the iterator at the inserted position (h's index);
  //       phase 1 uses that as pos_it.  my_arrival is the arrival_id
  //       enterS will stamp on h (enter_bba fetch_add's the counter;
  //       we hold S-exclusive so the counter is stable between the
  //       load and the fetch_add).
  //
  // Motivation: staging-9454 attribution showed main's phase1_wait
  // (S-shared at refill) dominated by workers stuck in phase-0
  // S-exclusive doing redtailBba (~17 % of worker wall).  Splitting
  // redtailBba off S-exclusive directly attacks that contention.
  // ------------------------------------------------------------------

  // (a) Unlocked.
  P->GetP(strat->lmBin);
  if (strat->homog) strat->initEcart(P);

  if (TEST_OPT_PROT) PrintS("s");

  // (b) S-shared.  Per-call output flags for redtailBba / redtailBbaAlsoLC_Z;
  // see the long comment above (task 361 run 581) for the rationale.
  bool rt_change = false;
  bool rt_retry = false;
  kt_S_lock_shared_phase0(ctx, thread_id);
#ifdef KTHREAD_INSTRUMENT
  long p0s_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif

  if (rField_is_Z(currRing) && !rHasLocalOrMixedOrdering(currRing))
    redtailBbaAlsoLC_Z(P, strat, &rt_change, &rt_retry);

  if (TEST_OPT_INTSTRATEGY)
  {
    P->pCleardenom();
    if ((TEST_OPT_REDSB) || (TEST_OPT_REDTAIL))
    {
#ifdef KTHREAD_INSTRUMENT
      long rt0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
      P->p = redtailBba(P, strat->S.end(), strat, withT,
                        !TEST_OPT_CONTENTSB,
                        &rt_change, &rt_retry);
#ifdef KTHREAD_INSTRUMENT
      if (KT_STATS(ctx)) redtail_accum += kt_now_ns() - rt0;
#endif
      P->pCleardenom();
      if (rt_change) P->t_p = NULL;
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
      P->p = redtailBba(P, strat->S.end(), strat, withT, FALSE,
                        &rt_change, &rt_retry);
#ifdef KTHREAD_INSTRUMENT
      if (KT_STATS(ctx)) redtail_accum += kt_now_ns() - rt0;
#endif
      if (rt_change) P->t_p = NULL;
    }
  }
  // rt_retry is intentionally not consumed here: the parallel drain
  // doesn't act on it directly — the existing serial bba path picks up
  // strat->completeReduce_retry (still OR-merged by redtailBba's
  // write_back) and triggers the retry from the outer loop.

  // SetShortExpVector is a pure LObject mutation but is needed before
  // enterT (which reads p.sev to populate sevT[atT]).  Compute it here
  // so we don't redo it under exclusive.
  bool will_enterS = ((!TEST_OPT_IDLIFT) || (pGetComp(P->p) <= strat->syzComp));
  if (will_enterS)
    P->SetShortExpVector();

#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx)) phase0_shared_work = kt_now_ns() - p0s_t0;
#endif
  strat->S.unlock_shared();

  // (c) S-exclusive: enterT + enterS.  Capture pos_it from enterS's
  // return (the iterator at h's inserted position) for phase 1.
  // pos_it must be default-constructible / assignable across the
  // shared/exclusive boundary; sBasisSet::iterator is a thin wrapper.
  sBasisSet::iterator pos_it = strat->S.end();
  uint64_t my_arrival = UINT64_MAX;
  bool did_enterS = false;
  if (will_enterS)
  {
    kt_S_lock_exclusive(ctx, thread_id);
#ifdef KTHREAD_INSTRUMENT
    long p0_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
    long et0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    enterT(*P, strat);
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx)) enterT_accum += kt_now_ns() - et0;
    long es0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    my_arrival = strat->arrival_counter.load(std::memory_order_relaxed);
    pos_it = strat->enterS(*P, strat, strat->T.size()-1, strat->S.end());
    did_enterS = true;
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx)) enterS_accum += kt_now_ns() - es0;
    if (KT_STATS(ctx)) phase0_work = kt_now_ns() - p0_t0;
#endif
    strat->S.unlock_exclusive();
  }
  // else: nothing entered T/S; phase 1 is skipped (did_enterS == false)

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
  // S-exclusive was already released at the end of phase 0(c).  Peer
  // drainers' phase 0 (which needs exclusive) can run between (c) and
  // the kt_S_lock_shared below; that is fine — their phase 0 runs,
  // they finish enterS, they drop exclusive, we get shared.
  // ------------------------------------------------------------------

  if (did_enterS)
  {
    // Stack-allocate thread-local B and pairtest-hit vector.  LSetChunk's
    // comparator indirects through strat->compareL (via its
    // CompareLObject::strat pointer), so we must set strat on the
    // local B exactly as skStrategy::skStrategy() does for strat->B.
    LSetChunk local_B;
    local_B.key_comp().strat = strat;
    std::vector<SElement*> local_pairtest_hits;

    // Save any previous thread-local overrides (NULL in practice; the
    // drain does not recurse through enterpairs, but defence in depth).
    LSetChunk* saved_B_override = t_local_B_override;
    uint64_t saved_my_arrival = t_local_my_arrival;
    std::vector<SElement*>* saved_pairtest_hits = t_local_pairtest_hits;

    t_local_B_override = &local_B;
    t_local_my_arrival = my_arrival;
    t_local_pairtest_hits = &local_pairtest_hits;

    kt_S_lock_shared(ctx, thread_id);

    // Phase-1 enterpairs ordering barrier.  Block until predecessors
    // (arrival_id < my_arrival) have completed their enterpairs,
    // restoring the serial-equivalent invariant: enterpairs(h_i) sees
    // an S in which no entry has been tombstoned by anyone with
    // arrival_id > i.  Wait BEFORE acquiring L_lock so we don't hold
    // L_lock during the wait — predecessors need L_lock for their
    // enterpairs work, so blocking on it would self-deadlock.  See
    // ~/project/docs/parallel-bba-deferred-enterpairs-clearS-violation.md.
    //
    // Task 360: the cond_wait used to share L_lock as its mutex, which
    // forced the wait to happen *after* L-lock acquisition — making
    // L_lock load-bearing for ordering rather than just data
    // protection.  The chunked-LSet refactor decouples the two:
    // enterpairs_order_lock is dedicated to the barrier, held only
    // for the wait/broadcast pair, and uncontended outside that
    // window.
    if (ctx->serialize_enterpairs) {
      pthread_mutex_lock(&ctx->enterpairs_order_lock);
      while (ctx->next_enterpairs_arrival_id.load(std::memory_order_relaxed)
             != my_arrival) {
        pthread_cond_wait(&ctx->enterpairs_order_cv, &ctx->enterpairs_order_lock);
      }
      pthread_mutex_unlock(&ctx->enterpairs_order_lock);
    }

    kt_L_lock_phase2(ctx, thread_id);

#ifdef KTHREAD_INSTRUMENT
    // Task 357 (run 570) option F: sample timestamp at L-lock
    // acquisition (the point kt_L_lock_phase2 returned) so phase1_l_ns
    // measures genuine L-hold time, not phase-1 work time.  Previously
    // p1_l_t0 was initialised to p1_t0 below, making phase1_l_ns
    // algebraically equal to phase1_ns — useless.
    long p1_l_acquired_ns = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif

#ifdef KTHREAD_INSTRUMENT
    long p1_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
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
    //
    // atR: use P->i_r (set under our exclusive lock by enterT at
    // kutil.cc:9092) instead of the racy strat->T.size()-1.  After we
    // released S-exclusive and reacquired S-shared above, peer drainers
    // can run their own phase-0 enterT, growing T.size() — so
    // strat->T.size()-1 now points at a peer's R-slot, not ours.
    // P->i_r is the persistent R-slot for our T entry.
    int atR_for_pairs = P->i_r;
    if (rField_is_Ring(currRing))
      superenterpairs(P->p, strat->S.size()-1, P->ecart, pos_it, strat, atR_for_pairs);
    else
      enterpairs(P->p, strat->S.size()-1, P->ecart, pos_it, strat, atR_for_pairs);

#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx))
    {
      long now = kt_now_ns();
      long ep_dt = now - ep0;
      enterpairs_accum += ep_dt;
      phase1_work = now - p1_t0;
      // phase1_l_ns is now accumulated at the actual L-lock unlock site
      // below (option F of task 357; run 570) using
      // p1_l_acquired_ns sampled right after kt_L_lock_phase2 returned.

      // Per-call enterpairs distribution (task 358; run 571).  my_arrival
      // captured under S-exclusive in phase 0 (line ~1207 above) is
      // still in scope here — used as the "id of the slowest call"
      // tag so we can locate it later if useful.
      ThreadStats &ts = KT_TS(ctx, thread_id);
      ts.enterpairs_count++;
      if (ep_dt > ts.enterpairs_max_ns) {
        ts.enterpairs_max_ns = ep_dt;
        ts.enterpairs_max_arrival = (long)my_arrival;
      }
      if (ep_dt < ts.enterpairs_min_ns)
        ts.enterpairs_min_ns = ep_dt;
    }
#endif

    // Safety: clear any residual local B / pairtest entries before
    // the LSetChunk / vector destructors run (normally they're empty, but
    // a control-flow short-circuit could leave stragglers).
    //
    // Restore previous overrides BEFORE the LSetChunk goes out of scope,
    // so nothing can reference it via the thread-local slot.
    t_local_B_override = saved_B_override;
    t_local_my_arrival = saved_my_arrival;
    t_local_pairtest_hits = saved_pairtest_hits;

    // Phase-1 enterpairs ordering barrier: advance the counter and
    // wake successors waiting for arrival_id == my_arrival + 1.
    // Task 360: broadcast under enterpairs_order_lock (not L_lock) so
    // waiters that woke spuriously and re-checked the predicate before
    // we updated it cannot miss the wakeup.  release-store on the
    // counter pairs with the cond_wait's lock-acquire fence to ensure
    // the new value is visible to a freshly-woken waiter.
    if (ctx->serialize_enterpairs) {
      pthread_mutex_lock(&ctx->enterpairs_order_lock);
      ctx->next_enterpairs_arrival_id.store(my_arrival + 1,
                                            std::memory_order_release);
      pthread_cond_broadcast(&ctx->enterpairs_order_cv);
      pthread_mutex_unlock(&ctx->enterpairs_order_lock);
    }

#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx))
    {
      // Sample L-hold delta just before the unlock — covers the entire
      // span L-lock-acquired → L-lock-released.  The barrier wait now
      // happens BEFORE L_lock acquisition (task 360), so phase1_l_ns
      // no longer includes barrier-wait time — that's a separate
      // bucket on the predecessor's enterpairs.  See task 357 (run 570)
      // option F for the original L-lock instrumentation choice.
      KT_TS(ctx, thread_id).phase1_l_ns +=
          kt_now_ns() - p1_l_acquired_ns;
    }
#endif
    strat->L.unlock();
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
    ts.phase0_shared_ns += phase0_shared_work;
    ts.phase1_ns += phase1_work;
    // phase1_l_ns is accumulated inside the phase-1 block at the L-lock
    // unlock site (task 357; run 570, option F).  Strictly contains phase1_ns.
    ts.phase_survivors++;
  }
  kt_my_thread_id = saved_kt_my_thread_id;
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
/*  Task 299 (run 506): locking moved into process_survivor_lobject, which runs  */
/*  explicit phase 0 / 1 / 2 blocks.  drain_survivor_queue no longer   */
/*  holds S.lock() or L_lock around the call.                          */
/*                                                                     */
/*  Each drain worker increments ctx->enterpairs_active on entry and   */
/*  decrements on exit. The counter is polled at termination time to   */
/*  tell when all drain workers are quiescent.                         */
/* ------------------------------------------------------------------ */

static void drain_survivor_queue(SweepContext *ctx, int thread_id)
{
#ifdef KTHREAD_INSTRUMENT
  long drain_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
  long drained_this_call = 0;
#endif

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
    // phase-0/1/2 structure (task 299; run 506).  It takes S-exclusive in
    // phase 0 and releases before returning.
    process_survivor_lobject(ctx, &P, thread_id);

#ifdef KTHREAD_INSTRUMENT
    drained_this_call++;
#endif

    // Wake anyone blocked waiting for L to refill.  pairs_available
    // is currently never waited on (vestigial; left in place for
    // future use).  The rdlock acquire here is just exclusion for the
    // cond_broadcast call — broadcast under any L lock is fine since
    // the cv has no waiters today.
    kt_L_lock_drain(ctx, thread_id);
    pthread_cond_broadcast(&ctx->pairs_available);
    ctx->strat->L.unlock();
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

  currRing = ctx->r;
  si_opt_1 = ctx->saved_si_opt_1;
  si_opt_2 = ctx->saved_si_opt_2;

#ifdef KTHREAD_INSTRUMENT
  // Bind thread-local thread_id for chainCritNormal instrumentation
  // in kutil.cc (task 358; run 571).  Note: process_survivor_lobject also
  // sets this — that path is taken from the survivor-drain which
  // workers may execute (task 305; run 512), so the saved/
  // restore there preserves this value.
  kt_my_thread_id = thread_id;
#endif

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
/*  kt_dump_stats — extracted from end of bba_parallel_loop.            */
/*                                                                     */
/*  Task 358 (run 571): callable mid-run from gdb    */
/*  via `call kt_dump_stats(kt_current_ctx)` so we can capture stats   */
/*  on workloads where SIGINT-driven shutdown is unreachable (main    */
/*  blocked in pthread_mutex_lock, etc).                               */
/*                                                                     */
/*  Behaviour is byte-identical to the previous in-line dump when      */
/*  called from parallel_shutdown: — same fields in the same order,    */
/*  same fprintf format strings.  When called mid-run, wall_ns         */
/*  reflects "time when this dump was taken" (kt_now_ns() -            */
/*  ctx->start_ns), not the eventual final wall.                      */
/*                                                                     */
/*  Safety: holds no locks, writes only to stderr / the optional CSV   */
/*  file, reads atomic counters.  Per-thread fields are written by     */
/*  exactly one thread each (the owning thread); a mid-run dump may   */
/*  read fields concurrently with their owner's writes — those are    */
/*  long-store races, fine for diagnostics (we may see a value         */
/*  that's slightly stale or torn; never crash).                      */
/* ------------------------------------------------------------------ */
#ifdef KTHREAD_INSTRUMENT
void kt_dump_stats(SweepContext *ctx)
{
  if (ctx == NULL) {
    fprintf(stderr, "[kthread-stats] kt_dump_stats: ctx is NULL "
                    "(no parallel-bba run in flight)\n");
    fflush(stderr);
    return;
  }
  if (!KT_STATS(ctx)) {
    fprintf(stderr, "[kthread-stats] kt_dump_stats: stats_enabled is false "
                    "(SINGULAR_KTHREAD_STATS unset)\n");
    fflush(stderr);
    return;
  }

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
  // Per-thread bucket totals.  Task 357 (run 570) reshaped this.
  //   sweep_ns         : workers' tile_pull_loop wall (umbrella).
  //   reduce_ns        : reduce_slot_from_sweep wall.
  //   drain_ns         : drain_survivor_queue wall.
  //   L_drain_wait_ns  : kt_L_lock_drain wait (post-survivor broadcast).
  //   L_term_wait_ns   : kt_L_lock_term wait (main termination probe).
  //   surv_q_wait_ns   : survivor_queue_mutex wait.
  //   tile_idle_ns     : workers' pthread_cond_wait on tiles_avail_cv
  //                      (closes the previously-unmeasured idle wedge).
  fprintf(stderr, "[kthread-stats] %-4s %12s %12s %12s %12s %12s %12s %12s %12s %10s %10s %10s %10s\n",
          "tid", "sweep_ns", "reduce_ns", "drain_ns",
          "L_drain_wait", "L_term_wait", "L_refill_wait",
          "surv_q_wait", "tile_idle_ns",
          "sweeps", "reduces", "drains", "survivors");
  for (int i = 0; i < tt; i++)
  {
    ThreadStats &ts = ctx->tstats[i];
    fprintf(stderr, "[kthread-stats] %-4d %12ld %12ld %12ld %12ld %12ld %12ld %12ld %12ld %10ld %10ld %10ld %10ld\n",
            i, ts.sweep_ns, ts.reduce_ns, ts.drain_ns,
            ts.L_lock_drain_wait_ns, ts.L_lock_term_wait_ns,
            ts.L_lock_refill_wait_ns,
            ts.surv_q_wait_ns, ts.tile_idle_ns,
            ts.sweep_count, ts.reduce_count, ts.drain_count,
            ts.drain_survivors);
  }

  // Main-thread umbrella + sub-buckets (task 357; run 570, option C).
  //   main_loop_ns     : umbrella (the whole bba_parallel_loop while(true)).
  //   refill_ns        : refill_and_publish wall.
  //   tile_help_ns     : tile_pull_loop(block=false) wall.
  //   publish_wait_ns  : pthread_cond_timedwait on slot_freed_cv wall.
  //   drain_ns         : (already shown above) drain_survivor_queue wall.
  //   reduce_ns        : (already shown above) typically zero on main.
  // Plus: L_term_wait_ns + L_drain_wait_ns from the per-thread block,
  // and small loop overhead — sum should approximate main_loop_ns.
  {
    ThreadStats &ts0 = ctx->tstats[0];
    fprintf(stderr, "[kthread-stats] main_pie: main_loop_ns=%ld refill_ns=%ld(n=%ld) "
                    "tile_help_ns=%ld publish_wait_ns=%ld(n=%ld)\n",
            ts0.main_loop_ns, ts0.refill_ns, ts0.refill_count,
            ts0.tile_help_ns, ts0.publish_wait_ns, ts0.publish_wait_count);
  }

  // Per-thread process_survivor breakdown (task 357; run 570, option B).
  // Previously only the run-total summary line existed; per-thread
  // attribution needed scaling by drain-share.  Now both are emitted:
  // per-thread first, then the run-total for back-compat with older
  // analyze.py.
  fprintf(stderr, "[kthread-stats] %-4s %14s %14s %14s %14s %14s %14s\n",
          "tid", "ps_redtail", "ps_enterT", "ps_enterpairs",
          "ps_enterS", "ps_other", "ps_total");
  long tot_rt = 0, tot_et = 0, tot_ep = 0, tot_es = 0, tot_oth = 0, tot_drain = 0;
  for (int i = 0; i < tt; i++)
  {
    ThreadStats &ts = ctx->tstats[i];
    long ps_total =
        ts.ps_redtail_ns + ts.ps_enterT_ns + ts.ps_enterpairs_ns +
        ts.ps_enterS_ns + ts.ps_other_ns;
    fprintf(stderr, "[kthread-stats] s%-3d %14ld %14ld %14ld %14ld %14ld %14ld\n",
            i, ts.ps_redtail_ns, ts.ps_enterT_ns, ts.ps_enterpairs_ns,
            ts.ps_enterS_ns, ts.ps_other_ns, ps_total);
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

  // Phase breakdown (task 300 / run 507 +
  // task 301 / run 508).  phase2_* renamed to
  // phase1_l_* (L-lock wait/held during phase 1, not a separate
  // phase).  ph1_cmax is the peak value of ctx->enterpairs_active
  // sampled by this thread on phase-1 entry (how many drainers
  // were concurrently in phase 1).
  //
  // task 359 (run 572): phase 0 was split into a pre-enterS S-shared
  // region (phase0_shared_*) and the original S-exclusive enterT+enterS
  // region (still phase0_*).  Old phase0_ns ≈ new phase0_ns +
  // new phase0_shared_ns.
  fprintf(stderr, "[kthread-stats] %-4s %14s %14s %14s %14s %14s %14s %14s %14s %10s %10s %10s %10s\n",
          "tid", "phase0_wait", "phase0_ns", "p0_sh_wait", "p0_sh_ns",
          "phase1_wait", "phase1_ns",
          "phase1_l_wait", "phase1_l_ns", "ph_surv", "ph_s_cas", "ph_l_cas",
          "ph1_cmax");
  for (int i = 0; i < tt; i++)
  {
    ThreadStats &ts = ctx->tstats[i];
    fprintf(stderr, "[kthread-stats] p%-3d %14ld %14ld %14ld %14ld %14ld %14ld %14ld %14ld %10ld %10ld %10ld %10ld\n",
            i, ts.phase0_wait_ns, ts.phase0_ns,
            ts.phase0_shared_wait_ns, ts.phase0_shared_ns,
            ts.phase1_wait_ns, ts.phase1_ns,
            ts.phase1_l_wait_ns, ts.phase1_l_ns, ts.phase_survivors,
            ts.phase_s_cas_fail, ts.phase_l_cas_fail,
            ts.phase1_concurrent_max);
  }

  // Worker-side drain participation (task 305; run 512).
  //   drain_survivors : # of survivors this thread handled
  //                     (was 0 for tid>0 before task 305 / run 512; post-task,
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

  // Per-call enterpairs / chainCritNormal distribution (task 358;
  // run 571).  Cumulative ps_enterpairs_ns lives
  // in the s-block above; this block exposes count / min / max so
  // single-call tail latency is visible.  *_min printed as 0 if no
  // call has been made yet (initial sentinel value LONG_MAX).
  fprintf(stderr, "[kthread-stats] %-4s %10s %14s %14s %14s %10s %14s %14s %14s %14s\n",
          "tid", "ep_count", "ep_min", "ep_max", "ep_max@",
          "cc_count", "cc_min", "cc_max", "cc_max@", "cc_total");
  for (int i = 0; i < tt; i++)
  {
    ThreadStats &ts = ctx->tstats[i];
    long ep_min = (ts.enterpairs_min_ns == LONG_MAX) ? 0 : ts.enterpairs_min_ns;
    long cc_min = (ts.chaincrit_min_ns  == LONG_MAX) ? 0 : ts.chaincrit_min_ns;
    fprintf(stderr, "[kthread-stats] e%-3d %10ld %14ld %14ld %14ld %10ld %14ld %14ld %14ld %14ld\n",
            i, ts.enterpairs_count, ep_min, ts.enterpairs_max_ns,
            ts.enterpairs_max_arrival,
            ts.chaincrit_count, cc_min, ts.chaincrit_max_ns,
            ts.chaincrit_max_arrival, ts.chaincrit_total_ns);
  }

  // L-set tombstone accumulation hypothesis (task 358; run 571):
  // chainCritNormal's slowness is dominated by
  // LSetChunk::filtered_iterator::advance() linear-scanning sev_flat_ from
  // 0 to flat_size, which includes tombstones.  LSetChunk::compact() runs
  // exactly once per GB run (at exitBuchMora), never during the loop.
  // Per-call duration should scale with L_flat (not L_live).
  //
  // L_live  = strat->L.physical_size()   (multiset node count)
  // L_flat  = strat->L.sev_flat_size()   (= sev_flat_.size(), the array
  //                                        the iterator walks, parallel
  //                                        to writable_set::flat_)
  // tombstone_ratio = (L_flat - L_live) / L_flat.
  //
  // Step 6 additions: chunk_count (k) — bounded by COMPACT_CHUNK_COUNT_
  // THRESHOLD (256), and wrapper-level compact stats.  pair_index_lookup
  // is shown in the per-thread block below.
  {
    long L_live = (long)ctx->strat->L.physical_size();
    long L_flat = (long)ctx->strat->L.sev_flat_size();
    long L_size = (long)ctx->strat->L.size();
    long L_chunks = (long)ctx->strat->L.chunk_count();
    double tombstone_ratio = (L_flat > 0)
        ? (double)(L_flat - L_live) / (double)L_flat : 0.0;
    fprintf(stderr,
            "[kthread-stats] L_live=%ld L_flat=%ld tombstone_ratio=%.4f "
            "L_size=%ld chunks=%ld\n",
            L_live, L_flat, tombstone_ratio, L_size, L_chunks);
    fprintf(stderr,
            "[kthread-stats] compact: calls=%ld total_ms=%.3f max_ms=%.3f "
            "last_k=%zu\n",
            ctx->strat->L.compact_call_count(),
            ctx->strat->L.compact_total_ns() / 1e6,
            ctx->strat->L.compact_max_ns() / 1e6,
            ctx->strat->L.last_compact_k());
  }
  // pair_index lookup wall (per-thread).  Bumped by isInPairsetL.
  fprintf(stderr, "[kthread-stats] %-4s %14s %10s\n",
          "tid", "piX_ns", "piX_count");
  for (int i = 0; i < tt; i++)
  {
    ThreadStats &ts = ctx->tstats[i];
    fprintf(stderr, "[kthread-stats] x%-3d %14ld %10ld\n",
            i, ts.pair_index_lookup_ns, ts.pair_index_lookup_count);
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

/* ------------------------------------------------------------------ */
/*  Main parallel loop                                                 */
/* ------------------------------------------------------------------ */

void bba_parallel_loop(SweepContext *ctx)
{
  kStrategy strat = ctx->strat;
  int nthreads = ctx->num_threads;

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
  {
    int cur_T = strat->T.size();
    long reserve_n = (long)cur_T * 4;
    if (reserve_n < 4096) reserve_n = 4096;
    strat->T.ensure_capacity((int)reserve_n);
    strat->sevT.ensure_capacity((int)reserve_n);
  }

  // Expose ctx as the parallel-mode runtime signal (read by
  // LSet::erase dispatch) and for gdb-driven kt_dump_stats calls.
  // Set unconditionally — required for correct erase semantics
  // regardless of KTHREAD_INSTRUMENT.
  kt_current_ctx = ctx;
#ifdef KTHREAD_INSTRUMENT
  // Main thread runs as thread_id 0.
  kt_my_thread_id = 0;
  if (KT_STATS(ctx))
  {
    ctx->start_ns = kt_now_ns();
    int tt = ctx->num_workers + 1;
    for (int i = 0; i < tt; i++)
    {
      memset(&ctx->tstats[i], 0, sizeof(ThreadStats));
      // Per-call min trackers (task 358; run 571): initial value is "no calls
      // yet, treat as +inf"; the dump prints 0 if we exit before any
      // value has overwritten this.
      ctx->tstats[i].enterpairs_min_ns = LONG_MAX;
      ctx->tstats[i].chaincrit_min_ns  = LONG_MAX;
    }
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

  // Task 304 (run 511): the startup-time pLength
  // pre-population loop is gone.  Every T entry in strat->T was
  // created by enterT (task 303; run 510), which computes pLength inline
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
  //      (task 305; run 512).
  //      enterT (invoked from drain_survivor_queue) computes pLength
  //      inline and release-publishes T slots, so no post-drain
  //      refresh is needed (task 304; run 511).
  //   3. Refill empty slots from L and (re)publish tiles for any slot
  //      that closer-reduce marked needs_republish.
  //   4. If nothing to do (all slots empty, L empty, queue empty):
  //        a. Help sweep any outstanding tiles (non-blocking) so we
  //           don't sit on main while workers still have work.
  //        b. If still idle after that, break out — we are done.
  //   5. Else, wait briefly on slot_freed_cv so closers can wake us.
  //
  // Task 305 (run 512) lifted the main-only-drain architectural
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
  // T/L access audit for main outside drain (task 305; run 512):
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
  //     (task 304; run 511) and is called from tile_pull_loop which main
  //     invokes with block=false — main is a tile reader at the same
  //     terms as workers.
  //   - strat->B: accessed via strat_B(strat) which returns the
  //     thread-local override inside phase 1, and strat->B elsewhere.
  //     Main outside drain never touches strat->B (only serial code
  //     before / after bba_parallel_loop does).
  // Error-path L-clear flags (set here, acted on after workers join).
  // Task 305 (run 512): we cannot touch strat->L while workers
  // may be inside chainCritNormal (which takes L_lock from phase 1).
  // Defer L mutation until after ctx->done + pthread_join below.
  bool clear_L_after_join = false;
  bool clear_slots_after_join = false;

#ifdef KTHREAD_INSTRUMENT
  // main_loop umbrella — bracket the whole while(true) so the per-thread
  // pie for tid 0 has a closed total to attribute against (task 357;
  // run 570, option C).
  long main_loop_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif

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
    // Task 304 (run 511): the post-drain pLength
    // refresh loop that ran here after every drain is gone — enterT
    // now computes pLength inline and release-publishes the slot in
    // one step (see kutil.cc enterT, task 303 / run 510).  Every T entry
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
#ifdef KTHREAD_INSTRUMENT
    long refill_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    int in_pipeline = refill_and_publish(ctx);
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx))
    {
      ThreadStats &ts = KT_TS(ctx, 0);
      ts.refill_ns += kt_now_ns() - refill_t0;
      ts.refill_count++;
    }
#endif

    // Step 3: termination check.
    if (in_pipeline == 0)
    {
      kt_surv_q_lock(ctx, 0);
      bool queue_empty = ctx->survivor_queue->empty();
      pthread_mutex_unlock(&ctx->survivor_queue_mutex);
      // Task 305 (run 512): L.empty() must be read under L's lock
      // to avoid racing with worker drainers' chainCritNormal pushes.
      // Step 5 (task 360): rdlock instead of the old L_lock mutex.
      kt_L_lock_term(ctx, 0);
      bool L_empty = strat->L.empty();
      strat->L.unlock();
      // Worker drainers may be mid-flight between popping a survivor
      // (queue becomes empty) and calling enterpairs (L gains new
      // pairs).  Without checking enterpairs_active, main sees
      // queue_empty && L_empty and breaks prematurely, losing the
      // pairs the worker is about to add.  rr+gdb MCP, 20 Apr 2026.
      int drainers_active = ctx->enterpairs_active.load(std::memory_order_acquire);
      if (queue_empty && L_empty && drainers_active == 0)
        break;
      // Else there is still work (drain produced survivors, L has
      // new entries, or a worker drainer hasn't finished adding
      // pairs) — loop back immediately.
      continue;
    }

    // Step 4: help sweep with tiles while waiting. Non-blocking: main
    // returns as soon as the tile cursor catches up with tile_end so
    // it can go back to refilling/draining.
#ifdef KTHREAD_INSTRUMENT
    long th_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    tile_pull_loop(ctx, 0, /*block=*/false);
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx))
      KT_TS(ctx, 0).tile_help_ns += kt_now_ns() - th_t0;
#endif

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
#ifdef KTHREAD_INSTRUMENT
      long pw_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
      pthread_cond_timedwait(&ctx->slot_freed_cv, &ctx->publish_lock, &ts);
#ifdef KTHREAD_INSTRUMENT
      if (KT_STATS(ctx))
      {
        ThreadStats &mts = KT_TS(ctx, 0);
        mts.publish_wait_ns += kt_now_ns() - pw_t0;
        mts.publish_wait_count++;
      }
#endif
    }
    pthread_mutex_unlock(&ctx->publish_lock);
  }

parallel_shutdown:
#ifdef KTHREAD_INSTRUMENT
  if (KT_STATS(ctx))
    KT_TS(ctx, 0).main_loop_ns += kt_now_ns() - main_loop_t0;
#endif
  // Signal done and wake all workers so they exit their tile_pull_loop.
  ctx->done.store(true, std::memory_order_release);
  pthread_mutex_lock(&ctx->publish_lock);
  pthread_cond_broadcast(&ctx->tiles_avail_cv);
  pthread_mutex_unlock(&ctx->publish_lock);

  for (int t = 0; t < ctx->num_workers; t++)
    pthread_join(ctx->threads[t], NULL);

  // Now workers are joined — safe to mutate strat->L and slots.
  // (Task 305 (run 512): chainCritNormal from worker drain is
  // the only writer of L outside main; joined workers cannot be inside
  // that critical section any longer.)
  if (clear_L_after_join)
  {
    strat->L.clear_and_erase();
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
  // Only emit the dump if stats are enabled — kt_dump_stats's own
  // "stats_enabled is false" early-return is a diagnostic for gdb
  // callers (who explicitly asked for the dump and expect feedback);
  // we don't want it cluttering normal runs' stderr (which would
  // also break the regress.cmd diff tests).
  if (KT_STATS(ctx)) kt_dump_stats(ctx);
  kt_my_thread_id = -1;
#endif
  // Clear the parallel-mode signal so post-loop serial paths see
  // NULL (LSet::erase dispatches to physical_erase) and a stray gdb
  // call to kt_dump_stats(kt_current_ctx) is a NULL no-op.
  kt_current_ctx = NULL;
}
