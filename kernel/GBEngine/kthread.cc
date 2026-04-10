/**
 * @file kthread.cc
 * @brief Parallel Groebner basis reduction — cooperative T-sweep design.
 *
 * All threads cooperatively sweep T[0..tl] for ALL active polynomials.
 * Each T entry is loaded into cache once and checked against N active
 * polynomials, instead of N threads each loading every T entry.
 *
 * Round structure:
 *   Fill:   Main thread fills active slots from L
 *   B0:     Barrier — start sweep
 *   Sweep:  All threads grab T indices via atomic sweep_cursor,
 *           check each T[j] against ALL active slots
 *   B1:     Barrier — sweep done
 *   Merge:  Main merges per-thread SweepResults into per-slot best
 *   Reduce: Main applies ksReducePoly for each slot with a reducer
 *   Update: Zero → mark empty; no reducer → survivor; still reducing → loop
 *   Process survivors, refill, next round
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

  ctx->max_active = nthreads;  // one slot per thread
  if (ctx->max_active < 1) ctx->max_active = 1;

  ctx->active = (ActivePoly *)calloc(ctx->max_active, sizeof(ActivePoly));
  for (int i = 0; i < ctx->max_active; i++)
  {
    ctx->active[i].occupied = false;
    ctx->active[i].is_survivor = false;
  }

  int total_threads = ctx->num_workers + 1;
  ctx->sweep_results = (SweepResult *)calloc(
      total_threads * ctx->max_active, sizeof(SweepResult));

  ctx->sweep_cursor.store(0, std::memory_order_relaxed);
  ctx->slot_counter.store(0, std::memory_order_relaxed);

  int barrier_count = ctx->num_workers + 1;
  if (barrier_count < 1) barrier_count = 1;
  pthread_barrier_init(&ctx->barrier_B0, NULL, barrier_count);
  pthread_barrier_init(&ctx->barrier_B1, NULL, barrier_count);
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

  // Asynchronous survivor drain: queue + mutex + CV + flag
  pthread_mutex_init(&ctx->survivor_queue_mutex, NULL);
  ctx->survivor_queue = new std::deque<LObject>();
  pthread_mutex_init(&ctx->enterpairs_mutex, NULL);
  pthread_cond_init(&ctx->pairs_available, NULL);
  ctx->enterpairs_active.store(false, std::memory_order_relaxed);
  ctx->stat_max_queue_depth.store(0, std::memory_order_relaxed);

  return ctx;
}

void sweep_context_destroy(SweepContext *ctx)
{
  if (ctx == NULL) return;
  pthread_barrier_destroy(&ctx->barrier_B0);
  pthread_barrier_destroy(&ctx->barrier_B1);
  pthread_barrier_destroy(&ctx->startup_barrier);
  pthread_mutex_destroy(&ctx->L_lock);
  pthread_mutex_destroy(&ctx->survivor_queue_mutex);
  pthread_mutex_destroy(&ctx->enterpairs_mutex);
  pthread_cond_destroy(&ctx->pairs_available);
  delete ctx->survivor_queue;
  free(ctx->active);
  free(ctx->sweep_results);
  free(ctx->threads);
  free(ctx->thread_ids);
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

    ap->P.PrepareRed(strat->use_buckets);
    return TRUE;
  }

  return FALSE;
}

/* ------------------------------------------------------------------ */
/*  Sweep phase: all threads cooperatively scan T for active polys.    */
/*  Each thread grabs T indices via atomic cursor and checks each T[j] */
/*  against ALL active slots, recording per-thread best reducers.      */
/* ------------------------------------------------------------------ */

static const int SWEEP_CHUNK = 64;

static void sweep_phase(SweepContext *ctx, int thread_id)
{
  kStrategy strat = ctx->strat;
  int tl = strat->T.size()-1;
  int max_active = ctx->max_active;

  while (true)
  {
    int start = ctx->sweep_cursor.fetch_add(SWEEP_CHUNK, std::memory_order_relaxed);
    if (start > tl) break;
    int end = start + SWEEP_CHUNK;
    if (end > tl + 1) end = tl + 1;

    for (int j = start; j < end; j++)
    {
      unsigned long sev_j = strat->sevT[j];

      for (int s = 0; s < max_active; s++)
      {
        if (!ctx->active[s].occupied) continue;
        if (ctx->active[s].is_survivor) continue;
        if (sev_j & ctx->active[s].not_sev) continue;
        if (!p_LmDivisibleBy(strat->T[j].p, ctx->active[s].P.p, currRing))
          continue;

        // Found a divisor for slot s
        SweepResult &sr = sweep_result(ctx, thread_id, s);
        if (sr.best_reducer < 0)
          sr.best_reducer = j;

        int ecart_j = strat->T[j].ecart;
        if (ecart_j <= ctx->active[s].P.ecart)
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
  }
}

/* ------------------------------------------------------------------ */
/*  Merge per-thread sweep results into per-slot best reducer.         */
/*  Called by main thread after sweep barrier.                          */
/* ------------------------------------------------------------------ */

static void merge_sweep_results(SweepContext *ctx)
{
  int total_threads = ctx->num_workers + 1;

  for (int s = 0; s < ctx->max_active; s++)
  {
    if (!ctx->active[s].occupied) continue;
    if (ctx->active[s].is_survivor) continue;

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
}

/* ------------------------------------------------------------------ */
/*  Reset sweep results for all threads and all slots.                 */
/* ------------------------------------------------------------------ */

static void reset_sweep_results(SweepContext *ctx)
{
  int total_threads = ctx->num_workers + 1;
  int total = total_threads * ctx->max_active;
  for (int i = 0; i < total; i++)
  {
    ctx->sweep_results[i].best_reducer = -1;
    ctx->sweep_results[i].best_good = -1;
    ctx->sweep_results[i].best_pLength = 0;
  }
}

/* ------------------------------------------------------------------ */
/*  Reduce one slot: apply ONE ksReducePoly step from sweep result.    */
/*  Slot stays active for the next cooperative sweep round if the      */
/*  polynomial is non-zero and not a survivor.                         */
/*  Called by worker threads in parallel (one slot per thread).        */
/* ------------------------------------------------------------------ */

static void reduce_slot_from_sweep(SweepContext *ctx, int slot)
{
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
/*  Parallel reduce phase: threads grab slots via atomic counter.      */
/* ------------------------------------------------------------------ */

static void reduce_phase_parallel(SweepContext *ctx)
{
  while (true)
  {
    int s = ctx->slot_counter.fetch_add(1, std::memory_order_relaxed);
    if (s >= ctx->max_active) break;

    ActivePoly *ap = &ctx->active[s];
    if (!ap->occupied) continue;
    if (ap->is_survivor) continue;

    reduce_slot_from_sweep(ctx, s);
  }
}

/* ------------------------------------------------------------------ */
/*  Fill empty active slots from L. Returns number of occupied slots.  */
/* ------------------------------------------------------------------ */

static int fill_active_slots(SweepContext *ctx)
{
  int occupied = 0;

  for (int s = 0; s < ctx->max_active; s++)
  {
    if (ctx->active[s].occupied)
    {
      if (!ctx->active[s].is_survivor)
        occupied++;
      continue;
    }

    // Try to fill this empty slot from L
    pthread_mutex_lock(&ctx->L_lock);
    BOOLEAN got = pop_and_prepare(ctx, &ctx->active[s]);
    pthread_mutex_unlock(&ctx->L_lock);

    if (got)
    {
      occupied++;
      ctx->stat_rounds.fetch_add(1, std::memory_order_relaxed);
    }
  }

  return occupied;
}

/* ------------------------------------------------------------------ */
/*  Process survivor (main thread only)                                */
/* ------------------------------------------------------------------ */

static void process_survivor(SweepContext *ctx, ActivePoly *ap)
{
  kStrategy strat = ctx->strat;
  BOOLEAN withT = ctx->withT;

  ap->P.GetP(strat->lmBin);
  if (strat->homog) strat->initEcart(&(ap->P));

  if (TEST_OPT_PROT) PrintS("s");

  int pos = posInS(strat, strat->S.size()-1, ap->P.p, ap->P.ecart);

  strat->redTailChange = FALSE;

  if (rField_is_Z(currRing) && !rHasLocalOrMixedOrdering(currRing))
    redtailBbaAlsoLC_Z(&(ap->P), strat->T.size()-1, strat);

  if (TEST_OPT_INTSTRATEGY)
  {
    ap->P.pCleardenom();
    if ((TEST_OPT_REDSB) || (TEST_OPT_REDTAIL))
    {
      ap->P.p = redtailBba(&(ap->P), pos - 1, strat, withT,
                            !TEST_OPT_CONTENTSB);
      ap->P.pCleardenom();
      if (strat->redTailChange) ap->P.t_p = NULL;
    }
  }
  else
  {
    ap->P.pNorm();
    if ((TEST_OPT_REDSB) || (TEST_OPT_REDTAIL))
    {
      ap->P.p = redtailBba(&(ap->P), pos - 1, strat, withT);
      if (strat->redTailChange) ap->P.t_p = NULL;
    }
  }

  if ((!TEST_OPT_IDLIFT) || (pGetComp(ap->P.p) <= strat->syzComp))
  {
    ap->P.SetShortExpVector();
    enterT(ap->P, strat);

    if (rField_is_Ring(currRing))
      superenterpairs(ap->P.p, strat->S.size()-1, ap->P.ecart, pos, strat, strat->T.size()-1);
    else
      enterpairs(ap->P.p, strat->S.size()-1, ap->P.ecart, pos, strat, strat->T.size()-1);

    strat->enterS(ap->P, strat, strat->T.size()-1, -1);
  }

  kDeleteLcm(&ap->P);
  ctx->stat_survivors.fetch_add(1, std::memory_order_relaxed);

  ap->P.Init();
  ap->occupied = false;
  ap->is_survivor = false;
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
 * Worker thread. Each round:
 *   1. B0: barrier_B0 (wait for main to signal start of round)
 *   2. sweep_phase: cooperatively scan T for all active polynomials
 *   3. B1: barrier_B1 (sweep done, main merges results)
 *   4. B0: barrier_B0 (wait for main to prepare reduce phase)
 *   5. reduce_phase_parallel: grab a slot, apply ONE ksReducePoly step
 *   6. B1: barrier_B1 (reduce done, main processes survivors + refills)
 * Workers check ctx->done after each B0 and exit if true.
 *
 * The 4-barrier protocol per round:
 *   B0: start sweep       B1: sweep done
 *   B0: start reduce      B1: reduce done
 * Between B1(reduce) and B0(sweep), main processes survivors and
 * refills empty slots. Multiple rounds may be needed per polynomial.
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

  while (true)
  {
    // Wait for main to signal start of batch (or done)
    pthread_barrier_wait(&ctx->barrier_B0);
    if (ctx->done.load(std::memory_order_acquire)) break;

    // Phase 1: cooperative sweep
    sweep_phase(ctx, thread_id);
    pthread_barrier_wait(&ctx->barrier_B1);

    // Phase 2: parallel reduction (main merges first, then signals)
    // Workers wait at B0 for main to finish merge, then reduce
    pthread_barrier_wait(&ctx->barrier_B0);
    if (ctx->done.load(std::memory_order_acquire)) break;

    reduce_phase_parallel(ctx);
    pthread_barrier_wait(&ctx->barrier_B1);

    // Main processes survivors between B1 and next B0
  }

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

  // Main loop: cooperative sweep rounds
  //
  // Each round:
  //   Fill:    main fills empty active slots from L
  //   B0:      barrier — start sweep
  //   Sweep:   all threads cooperatively scan T (atomic cursor)
  //   B1:      barrier — sweep done
  //   Merge:   main merges per-thread SweepResults
  //   Reduce:  main applies ksReducePoly for slots with reducers
  //   Process: survivors → enterT/enterpairs/enterS (serialized)
  //   Loop if any slots still need reduction
  //
  while (true)
  {
    if (siCntrlc)
    {
      while (!strat->L.empty()) strat->L.pop_and_erase();
      strat->noClearS = TRUE;
      if (ctx->num_workers > 0)
      {
        // Workers are waiting at B0; signal done and release them
        ctx->done.store(true, std::memory_order_release);
        pthread_barrier_wait(&ctx->barrier_B0);
      }
      goto parallel_cleanup;
    }

    // Fill phase: fill empty active slots from L
    int occupied = fill_active_slots(ctx);

    // If no occupied slots (all done or L empty), check for remaining work
    if (occupied == 0)
    {
      // Process any remaining survivors
      bool had_survivor = false;
      for (int i = 0; i < ctx->max_active; i++)
      {
        if (ctx->active[i].occupied && ctx->active[i].is_survivor)
        {
          process_survivor(ctx, &ctx->active[i]);
          had_survivor = true;
        }
      }

      if (had_survivor)
      {
        // Re-compute pLength for new T entries
        for (int j = 0; j < strat->T.size(); j++)
        {
          if (strat->T[j].pLength <= 0)
            strat->T[j].pLength = pLength(strat->T[j].p ? strat->T[j].p : strat->T[j].t_p);
        }
        // Try again — process_survivor may have added to L via enterpairs
        continue;
      }

      // Nothing left to do
      if (strat->L.empty())
        break;
      // L has new entries from enterpairs, loop back to fill
      continue;
    }

    // Prepare for sweep round
    reset_sweep_results(ctx);
    ctx->sweep_cursor.store(0, std::memory_order_relaxed);

    // B0: signal start of sweep (workers are waiting here)
    if (ctx->num_workers > 0)
      pthread_barrier_wait(&ctx->barrier_B0);

    // Phase 1: cooperative sweep — all threads scan T for active polys
    sweep_phase(ctx, 0);

    // B1: sweep done
    pthread_barrier_wait(&ctx->barrier_B1);

    // Merge per-thread sweep results into per-slot best reducer
    merge_sweep_results(ctx);

    // Prepare slot counter for parallel reduction
    ctx->slot_counter.store(0, std::memory_order_relaxed);

    // B0: signal start of reduction phase (workers waiting)
    if (ctx->num_workers > 0)
      pthread_barrier_wait(&ctx->barrier_B0);

    // Phase 2: parallel reduction — each thread grabs a slot and applies
    // ONE ksReducePoly step. Slots that are still non-zero and non-survivor
    // stay active and will be swept again on the next round.
    reduce_phase_parallel(ctx);

    // B1: reduction done
    pthread_barrier_wait(&ctx->barrier_B1);

    if (strat->overflow || errorreported)
    {
      while (!strat->L.empty()) strat->L.pop_and_erase();
      for (int i = 0; i < ctx->max_active; i++)
        ctx->active[i].occupied = false;
      break;
    }

    // Process all survivors (main thread only, serialized)
    bool had_survivor = false;
    for (int i = 0; i < ctx->max_active; i++)
    {
      if (ctx->active[i].occupied && ctx->active[i].is_survivor)
      {
        process_survivor(ctx, &ctx->active[i]);
        had_survivor = true;
      }
    }

    // Re-compute pLength for any new T entries added by process_survivor
    if (had_survivor)
    {
      for (int j = 0; j < strat->T.size(); j++)
      {
        if (strat->T[j].pLength <= 0)
          strat->T[j].pLength = pLength(strat->T[j].p ? strat->T[j].p : strat->T[j].t_p);
      }
    }

    // Slots may still be occupied and need more reduction rounds.
    // Loop back: fill empty slots from L, sweep again for all active slots.
  }

  // Signal workers to exit: hit B0 with ctx->done=true so they break
  if (ctx->num_workers > 0)
  {
    ctx->done.store(true, std::memory_order_release);
    pthread_barrier_wait(&ctx->barrier_B0);
  }

parallel_cleanup:
  for (int t = 0; t < ctx->num_workers; t++)
    pthread_join(ctx->threads[t], NULL);

  strat->posInT = ctx->saved_posInT;
  SI_RESTORE_OPT1(ctx->saved_opt1);

  if (TEST_OPT_PROT)
  {
    Print("\n// parallel bba: polys=%ld, reductions=%ld, zeros=%ld, "
          "survivors=%ld, threads=%d (workers=%d)\n",
          ctx->stat_rounds.load(),
          ctx->stat_reductions.load(),
          ctx->stat_zeros.load(),
          ctx->stat_survivors.load(),
          nthreads, ctx->num_workers);
  }
}
