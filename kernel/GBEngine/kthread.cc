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
#else
#  define KT_STATS(ctx)       (false)
#  define KT_TIME_START(var)  ((void)0)
#  define KT_TIME_DELTA(var)  (0L)
static inline void kt_L_lock(SweepContext *ctx, int /*tid*/) { pthread_mutex_lock(&ctx->L_lock); }
static inline void kt_surv_q_lock(SweepContext *ctx, int /*tid*/) { pthread_mutex_lock(&ctx->survivor_queue_mutex); }
static inline void kt_record_reduce(SweepContext*, int, long, long) {}
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

static BOOLEAN pop_and_prepare(SweepContext *ctx, ActivePoly *ap)
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
    // At this stage the sweep still reads the global strat->T.size()-1,
    // but we capture per-slot sl_snapshot so later milestones can tile
    // against it without reading a moving global bound.
    ap->sl_snapshot = strat->T.size() - 1;

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
      // Try to fill from L.
      kt_L_lock(ctx, 0);
      BOOLEAN got = pop_and_prepare(ctx, ap);
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
      // Re-snapshot tl here, in case T has grown since last pass.
      ap->sl_snapshot = ctx->strat->T.size() - 1;
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

static void process_survivor_lobject(SweepContext *ctx, LObject *P, int thread_id)
{
  kStrategy strat = ctx->strat;
  BOOLEAN withT = ctx->withT;
  (void)thread_id;

#ifdef KTHREAD_INSTRUMENT
  long ps_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
  long redtail_accum = 0;
  long enterT_accum = 0;
  long enterpairs_accum = 0;
  long enterS_accum = 0;
#endif

  P->GetP(strat->lmBin);
  if (strat->homog) strat->initEcart(P);

  if (TEST_OPT_PROT) PrintS("s");

  int pos = posInS(strat, strat->S.size()-1, P->p, P->ecart);

  strat->redTailChange = FALSE;

  if (rField_is_Z(currRing) && !rHasLocalOrMixedOrdering(currRing))
    redtailBbaAlsoLC_Z(P, strat->T.size()-1, strat);

  if (TEST_OPT_INTSTRATEGY)
  {
    P->pCleardenom();
    if ((TEST_OPT_REDSB) || (TEST_OPT_REDTAIL))
    {
#ifdef KTHREAD_INSTRUMENT
      long rt0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
      P->p = redtailBba(P, pos - 1, strat, withT,
                        !TEST_OPT_CONTENTSB);
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
      P->p = redtailBba(P, pos - 1, strat, withT);
#ifdef KTHREAD_INSTRUMENT
      if (KT_STATS(ctx)) redtail_accum += kt_now_ns() - rt0;
#endif
      if (strat->redTailChange) P->t_p = NULL;
    }
  }

  if ((!TEST_OPT_IDLIFT) || (pGetComp(P->p) <= strat->syzComp))
  {
    P->SetShortExpVector();
#ifdef KTHREAD_INSTRUMENT
    long et0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    enterT(*P, strat);
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx)) enterT_accum += kt_now_ns() - et0;
    long ep0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif

    if (rField_is_Ring(currRing))
      superenterpairs(P->p, strat->S.size()-1, P->ecart, pos, strat, strat->T.size()-1);
    else
      enterpairs(P->p, strat->S.size()-1, P->ecart, pos, strat, strat->T.size()-1);

#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx)) enterpairs_accum += kt_now_ns() - ep0;
    long es0 = KT_STATS(ctx) ? kt_now_ns() : 0;
#endif
    strat->enterS(*P, strat, strat->T.size()-1, -1);
#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx)) enterS_accum += kt_now_ns() - es0;
#endif
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
/*  Drain the survivor FIFO: pop one survivor at a time and process    */
/*  it under the S+L lock pair.                                        */
/*                                                                     */
/*  Multi-drainer safe: any number of threads may call this function   */
/*  concurrently. Each call pops at most as many survivors as remain   */
/*  in the queue at the time of the pop, processing each one while     */
/*  holding both strat->S.lock() and ctx->L_lock. The lock pair is     */
/*  released between survivors, so fill_active_slots (which only       */
/*  needs L_lock) can interleave with the drain, and two drain         */
/*  workers can execute process_survivor_lobject in pipelined          */
/*  fashion (one draining while another is blocked waiting on the      */
/*  locks for the next survivor).                                      */
/*                                                                     */
/*  Lock order: S.lock(), then L_lock. Always release L_lock first,    */
/*  then S.lock(). This avoids deadlock with any future code that      */
/*  takes S.lock() while holding L_lock. Currently no such code        */
/*  exists: fill_active_slots takes only L_lock.                       */
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
  kStrategy strat = ctx->strat;

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

    // Acquire S.lock() and L_lock in that order. S.lock() serializes
    // iteration and mutation of S across drain workers; L_lock
    // serializes mutation of L against itself and fill_active_slots.
    // The pair is held only for this one survivor, not across the
    // whole drain — so between survivors another drain worker can
    // take over and the main thread can fill from L.
    strat->S.lock();
    kt_L_lock(ctx, thread_id);

    process_survivor_lobject(ctx, &P, thread_id);

    pthread_mutex_unlock(&ctx->L_lock);
    strat->S.unlock();

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

  // Pre-compute pLength for all T entries to avoid lazy init during parallel phase
  // (must be done BEFORE startup barrier so workers don't race ahead)
  for (int j = 0; j < strat->T.size(); j++)
  {
    if (strat->T[j].pLength <= 0)
      strat->T[j].pLength = pLength(strat->T[j].p ? strat->T[j].p : strat->T[j].t_p);
  }

  // Wait for all workers to start before entering main loop
  if (ctx->num_workers > 0)
    pthread_barrier_wait(&ctx->startup_barrier);

  // Continuous refill loop (milestone d, task 283). No barriers.
  //
  // Each iteration:
  //   1. If siCntrlc: shutdown path.
  //   2. Drain any survivors from the FIFO (serially, on main).
  //   3. Refresh pLength for any new T entries added by drain.
  //   4. Refill empty slots from L and (re)publish tiles for any slot
  //      that closer-reduce marked needs_republish.
  //   5. If nothing to do (all slots empty, L empty, queue empty):
  //        a. Help sweep any outstanding tiles (non-blocking) so we
  //           don't sit on main while workers still have work.
  //        b. If still idle after that, break out — we are done.
  //   6. Else, wait briefly on slot_freed_cv so closers can wake us.
  //
  // Main only runs the drain; workers never touch strat->T / strat->S
  // / strat->L (except the L_lock for pop_and_prepare, which happens
  // only on the main thread via refill_and_publish).
  while (true)
  {
    if (siCntrlc)
    {
      while (!strat->L.empty()) strat->L.pop_and_erase();
      strat->noClearS = TRUE;
      goto parallel_shutdown;
    }

    if (strat->overflow || errorreported)
    {
      while (!strat->L.empty()) strat->L.pop_and_erase();
      for (int i = 0; i < ctx->max_active; i++)
        ctx->active[i].occupied = false;
      goto parallel_shutdown;
    }

    // Step 1: drain survivor FIFO on main thread.
    {
      kt_surv_q_lock(ctx, 0);
      bool queue_empty = ctx->survivor_queue->empty();
      pthread_mutex_unlock(&ctx->survivor_queue_mutex);
      if (!queue_empty)
      {
        drain_survivor_queue(ctx, 0);
        // Refresh pLength for any T entries added by process_survivor.
        for (int j = 0; j < strat->T.size(); j++)
        {
          if (strat->T[j].pLength <= 0)
            strat->T[j].pLength = pLength(strat->T[j].p ? strat->T[j].p : strat->T[j].t_p);
        }
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
      if (queue_empty && strat->L.empty())
        break;
      // Else there is still work (drain produced survivors, or L has
      // new entries) — loop back immediately.
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
