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
#include <sched.h>

/* ------------------------------------------------------------------ */
/*  Instrumentation helpers (task 482)                                 */
/* ------------------------------------------------------------------ */
#ifdef KTHREAD_INSTRUMENT
#  define KT_STATS(ctx)       ((ctx)->stats_enabled)
#  define KT_TS(ctx, tid)     ((ctx)->tstats[tid])
#  define KT_TIME_START(var)  long var = kt_now_ns()
#  define KT_TIME_DELTA(var)  (kt_now_ns() - (var))

static inline void kt_barrier_wait_B0(SweepContext *ctx, int thread_id)
{
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    pthread_barrier_wait(&ctx->barrier_B0);
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.wait_B0_ns += dt;
    ts.wait_B0_count++;
    if (dt > ts.wait_B0_max_ns) ts.wait_B0_max_ns = dt;
  }
  else
  {
    pthread_barrier_wait(&ctx->barrier_B0);
  }
}

static inline void kt_barrier_wait_B1(SweepContext *ctx, int thread_id)
{
  if (KT_STATS(ctx))
  {
    long t0 = kt_now_ns();
    pthread_barrier_wait(&ctx->barrier_B1);
    long dt = kt_now_ns() - t0;
    ThreadStats &ts = KT_TS(ctx, thread_id);
    ts.wait_B1_ns += dt;
    ts.wait_B1_count++;
    if (dt > ts.wait_B1_max_ns) ts.wait_B1_max_ns = dt;
  }
  else
  {
    pthread_barrier_wait(&ctx->barrier_B1);
  }
}

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
static inline void kt_barrier_wait_B0(SweepContext *ctx, int /*tid*/) { pthread_barrier_wait(&ctx->barrier_B0); }
static inline void kt_barrier_wait_B1(SweepContext *ctx, int /*tid*/) { pthread_barrier_wait(&ctx->barrier_B1); }
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
  pthread_barrier_destroy(&ctx->barrier_B0);
  pthread_barrier_destroy(&ctx->barrier_B1);
  pthread_barrier_destroy(&ctx->startup_barrier);
  pthread_mutex_destroy(&ctx->L_lock);
  pthread_mutex_destroy(&ctx->survivor_queue_mutex);
  pthread_cond_destroy(&ctx->pairs_available);
  delete ctx->survivor_queue;
  free(ctx->active);
  free(ctx->sweep_results);
  free(ctx->threads);
  free(ctx->thread_ids);
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
/*  Parallel reduce phase: threads grab slots via atomic counter.      */
/* ------------------------------------------------------------------ */

static void reduce_phase_parallel(SweepContext *ctx, int thread_id)
{
  while (true)
  {
    int s = ctx->slot_counter.fetch_add(1, std::memory_order_relaxed);
    if (s >= ctx->max_active) break;

    ActivePoly *ap = &ctx->active[s];
    if (!ap->occupied) continue;
    if (ap->is_survivor) continue;

#ifdef KTHREAD_INSTRUMENT
    long rstart = kt_now_ns();
    reduce_slot_from_sweep(ctx, s, thread_id);
    long rdur = kt_now_ns() - rstart;
    if (KT_STATS(ctx))
    {
      ThreadStats &ts = KT_TS(ctx, thread_id);
      ts.reduce_ns += rdur;
      ts.reduce_count++;
      kt_record_reduce(ctx, thread_id, rstart - ctx->start_ns, rdur);
    }
#else
    reduce_slot_from_sweep(ctx, s, thread_id);
#endif
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
    kt_L_lock(ctx, 0);
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
    kt_barrier_wait_B0(ctx, thread_id);
    if (ctx->done.load(std::memory_order_acquire)) break;

    // Phase 1: cooperative sweep
    //
    // Task 483: workers do NOT drain during sweep. The drain runs
    // only on the main thread, between reduce-B1 and the next
    // round's fill_active_slots. This avoids racing drain-side
    // mutations of strat->T / strat->S / strat->L against main's
    // pLength refresh and sweep read of strat->T.
    //
    // The multi-drainer capability built into drain_survivor_queue
    // (S+L lock pair) is retained for future use — e.g., if a
    // follow-up task adds a dedicated drain barrier or continuous
    // design. On the measured workloads the survivor queue is empty
    // 99% of the time (parallel-bba-bottleneck-analysis.md), so
    // single-drainer on main is not a performance loss.
#ifdef KTHREAD_INSTRUMENT
    long sw_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
    sweep_phase(ctx, thread_id);
    if (KT_STATS(ctx))
    {
      ThreadStats &ts = KT_TS(ctx, thread_id);
      ts.sweep_ns += kt_now_ns() - sw_t0;
      ts.sweep_count++;
    }
#else
    sweep_phase(ctx, thread_id);
#endif
    kt_barrier_wait_B1(ctx, thread_id);

    // Phase 2: parallel reduction (main merges first, then signals)
    // Workers wait at B0 for main to finish merge, then reduce
    kt_barrier_wait_B0(ctx, thread_id);
    if (ctx->done.load(std::memory_order_acquire)) break;

    reduce_phase_parallel(ctx, thread_id);
    kt_barrier_wait_B1(ctx, thread_id);
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
      // Queue any remaining active-slot survivors.
      for (int i = 0; i < ctx->max_active; i++)
      {
        if (ctx->active[i].occupied && ctx->active[i].is_survivor)
          queue_survivor(ctx, &ctx->active[i]);
      }

      // Idle-drain: main thread drains any remaining survivors when
      // it has nothing else to do. With the task-483 multi-drainer
      // model, worker threads can also drain between rounds, but
      // at termination there are no more rounds, so the main thread
      // must do a final sweep to ensure the queue is empty. It may
      // run concurrently with workers that are still inside a drain
      // call — the S+L lock pair in drain_survivor_queue serializes
      // the inner work.
      bool drained_something = false;
      {
        kt_surv_q_lock(ctx, 0);
        bool queue_empty = ctx->survivor_queue->empty();
        pthread_mutex_unlock(&ctx->survivor_queue_mutex);
        if (!queue_empty)
        {
          drain_survivor_queue(ctx, 0);
          drained_something = true;
        }
      }

      // Wait for any in-flight worker drains to finish before we
      // decide the computation is done. A worker may have been mid-
      // drain when we polled strat->L and found it empty; if so, the
      // worker could add new entries to L before we terminate.
      // Spin briefly on the counter; on the measured workloads
      // drains complete in microseconds.
      while (ctx->enterpairs_active.load(std::memory_order_acquire) > 0)
      {
        // Yield to let drain workers run. This loop runs only at
        // termination, not in the hot path.
        sched_yield();
      }

      if (drained_something)
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

      // Nothing queued and no survivors. Either L is now populated
      // (loop back and fill) or we are done.
      if (strat->L.empty())
        break;
      // L has new entries, loop back to fill
      continue;
    }

    // Prepare for sweep round
    reset_sweep_results(ctx);
    ctx->sweep_cursor.store(0, std::memory_order_relaxed);

#ifdef KTHREAD_INSTRUMENT
    long round_t0 = 0;
    int depth_start = 0;
    if (KT_STATS(ctx))
    {
      round_t0 = kt_now_ns();
      pthread_mutex_lock(&ctx->survivor_queue_mutex);
      depth_start = (int)ctx->survivor_queue->size();
      pthread_mutex_unlock(&ctx->survivor_queue_mutex);
      pthread_mutex_lock(&ctx->stats_lock);
      ctx->reduces_this_round->clear();
      pthread_mutex_unlock(&ctx->stats_lock);
    }
#endif

    // B0: signal start of sweep (workers are waiting here)
    if (ctx->num_workers > 0)
      kt_barrier_wait_B0(ctx, 0);

    // Phase 1: cooperative sweep — all threads scan T for active polys
#ifdef KTHREAD_INSTRUMENT
    long sw_t0 = KT_STATS(ctx) ? kt_now_ns() : 0;
    sweep_phase(ctx, 0);
    long main_sweep_ns = 0;
    if (KT_STATS(ctx))
    {
      main_sweep_ns = kt_now_ns() - sw_t0;
      KT_TS(ctx, 0).sweep_ns += main_sweep_ns;
      KT_TS(ctx, 0).sweep_count++;
    }
#else
    sweep_phase(ctx, 0);
#endif

    // B1: sweep done
    kt_barrier_wait_B1(ctx, 0);

    // Merge per-thread sweep results into per-slot best reducer
    merge_sweep_results(ctx);

    // Prepare slot counter for parallel reduction
    ctx->slot_counter.store(0, std::memory_order_relaxed);

    // B0: signal start of reduction phase (workers waiting)
    if (ctx->num_workers > 0)
      kt_barrier_wait_B0(ctx, 0);

    // Phase 2: parallel reduction — each thread grabs a slot and applies
    // ONE ksReducePoly step. Slots that are still non-zero and non-survivor
    // stay active and will be swept again on the next round.
    reduce_phase_parallel(ctx, 0);

    // B1: reduction done
    kt_barrier_wait_B1(ctx, 0);

    if (strat->overflow || errorreported)
    {
      while (!strat->L.empty()) strat->L.pop_and_erase();
      for (int i = 0; i < ctx->max_active; i++)
        ctx->active[i].occupied = false;
      break;
    }

    // Queue all survivors into the FIFO, then drain them on the
    // main thread. Task 483 moved the drain onto main so it runs
    // single-threaded, avoiding races between worker enterT/enterS
    // mutations of strat->T / strat->S and main's pLength refresh
    // and sweep-time read of strat->T. See the worker_thread
    // comment for rationale.
    for (int i = 0; i < ctx->max_active; i++)
    {
      if (ctx->active[i].occupied && ctx->active[i].is_survivor)
        queue_survivor(ctx, &ctx->active[i]);
    }

    // Drain the queue on main while workers wait at the next B0.
    if (!ctx->survivor_queue->empty())
      drain_survivor_queue(ctx, 0);

#ifdef KTHREAD_INSTRUMENT
    if (KT_STATS(ctx))
    {
      RoundRecord rr;
      rr.round_id = (long)ctx->rounds->size();
      rr.timestamp_ns = round_t0 - ctx->start_ns;
      rr.queue_depth_start = depth_start;
      pthread_mutex_lock(&ctx->survivor_queue_mutex);
      rr.queue_depth_end = (int)ctx->survivor_queue->size();
      pthread_mutex_unlock(&ctx->survivor_queue_mutex);
      rr.active_slots = occupied;
      rr.sweep_ns_main = main_sweep_ns;
      rr.round_total_ns = kt_now_ns() - round_t0;
      pthread_mutex_lock(&ctx->stats_lock);
      rr.reductions = (int)ctx->reduces_this_round->size();
      long mn = -1, mx = 0, sm = 0;
      for (auto &ev : *ctx->reduces_this_round)
      {
        if (mn < 0 || ev.duration_ns < mn) mn = ev.duration_ns;
        if (ev.duration_ns > mx) mx = ev.duration_ns;
        sm += ev.duration_ns;
      }
      rr.min_reduce_ns = (mn < 0) ? 0 : mn;
      rr.max_reduce_ns = mx;
      rr.sum_reduce_ns = sm;
      pthread_mutex_unlock(&ctx->stats_lock);
      ctx->rounds->push_back(rr);
    }
#endif

    // Re-compute pLength for any new T entries. A worker may have
    // drained the queue during the sweep phase, adding new T entries
    // via enterT, so we unconditionally refresh pLength here.
    for (int j = 0; j < strat->T.size(); j++)
    {
      if (strat->T[j].pLength <= 0)
        strat->T[j].pLength = pLength(strat->T[j].p ? strat->T[j].p : strat->T[j].t_p);
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
