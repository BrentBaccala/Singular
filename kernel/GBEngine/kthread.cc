/**
 * @file kthread.cc
 * @brief Parallel Groebner basis reduction — poly-per-thread design.
 *
 * Each thread independently reduces a different polynomial from L,
 * sharing T[] and sevT[] as read-only. When a polynomial survives
 * (no divisor found), it enters a survivor queue processed by the
 * main thread (enterT/enterpairs/enterS). When a polynomial reduces
 * to zero, the thread immediately grabs the next polynomial from L.
 *
 * Synchronization:
 *   - L pops are under L_lock (cheap, since L.pop() is O(1))
 *   - Survivor processing is serialized on main thread
 *   - T/sevT/S are read-only during parallel reduction; modified
 *     only during survivor processing (serialized)
 *   - Barrier between batch reduce and batch survivor processing
 *
 * Thread safety relies on:
 *   - posInT0 (append at end, no memmove)
 *   - sevT/R/T written BEFORE tl++ (compiler barrier in enterT)
 *   - --disable-omalloc build (thread-safe malloc/free)
 *   - currRing and si_opt_1/si_opt_2 are __thread
 *   - OPT_REDTHROUGH forced (no lazy pushback to L from reduction)
 *   - tailRing pre-expanded before parallel phase
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
  pthread_barrier_init(&ctx->barrier, NULL, barrier_count);
  pthread_barrier_init(&ctx->startup_barrier, NULL, barrier_count);

  int alloc_n = ctx->num_workers > 0 ? ctx->num_workers : 1;
  ctx->threads = (pthread_t *)calloc(alloc_n, sizeof(pthread_t));
  ctx->thread_ids = (int *)calloc(alloc_n, sizeof(int));

  pthread_mutex_init(&ctx->L_lock, NULL);

  ctx->stat_reductions.store(0, std::memory_order_relaxed);
  ctx->stat_zeros.store(0, std::memory_order_relaxed);
  ctx->stat_survivors.store(0, std::memory_order_relaxed);
  ctx->stat_rounds.store(0, std::memory_order_relaxed);

  return ctx;
}

void sweep_context_destroy(SweepContext *ctx)
{
  if (ctx == NULL) return;
  pthread_barrier_destroy(&ctx->barrier);
  pthread_barrier_destroy(&ctx->startup_barrier);
  pthread_mutex_destroy(&ctx->L_lock);
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
      ksCreateSpoly(&(ap->P), NULL, strat->use_buckets,
                    strat->tailRing, m1, m2, &strat->R);
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
/*  Find best reducer for a polynomial by scanning T (sequential).     */
/*  This is equivalent to kFindDivisibleByInT_ecart.                   */
/* ------------------------------------------------------------------ */

static int find_reducer(SweepContext *ctx, ActivePoly *ap)
{
  kStrategy strat = ctx->strat;
  int tl = strat->tl;
  unsigned long not_sev = ap->not_sev;
  poly p_lm = ap->P.p;
  int p_ecart = ap->P.ecart;

  int best_reducer = -1;    // any divisor
  int best_good = -1;       // good ecart, shortest pLength
  int best_pLength = 0;

  for (int j = 0; j <= tl; j++)
  {
    unsigned long sev_j = strat->sevT[j];
    if (sev_j & not_sev) continue;
    if (!p_LmDivisibleBy(strat->T[j].p, p_lm, currRing)) continue;

    if (best_reducer < 0)
      best_reducer = j;

    int ecart_j = strat->T[j].ecart;
    if (ecart_j <= p_ecart)
    {
      int pLen = strat->T[j].pLength;
      if (pLen <= 0)
        pLen = 3;  // avoid thread-unsafe T[j].GetpLength()

      if (best_good < 0 || pLen < best_pLength)
      {
        best_good = j;
        best_pLength = pLen;
      }
      if (best_pLength <= 2) break;  // can't do better
    }
  }

  return (best_good >= 0) ? best_good : best_reducer;
}

/* ------------------------------------------------------------------ */
/*  Reduce one polynomial fully (like redHoney). Returns:              */
/*    0 = reduced to zero                                              */
/*    1 = survived (no reducer found)                                  */
/* ------------------------------------------------------------------ */

static int reduce_fully(SweepContext *ctx, ActivePoly *ap)
{
  kStrategy strat = ctx->strat;

  while (true)
  {
    int best = find_reducer(ctx, ap);

    if (best < 0)
    {
      // No reducer found — survivor
      ap->is_survivor = true;
      return 1;
    }

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
      return 0;
    }

    // Update ecart (redHoney lines 2421-2436)
    ap->P.SetShortExpVector();
    int h_d = ap->P.SetpFDeg();
    if (ei <= ap->P.ecart)
      ap->P.ecart = ap->d - h_d;
    else
      ap->P.ecart = ap->d - h_d + ei - ap->P.ecart;
    ap->pass++;
    ap->d = h_d + ap->P.ecart;

    // Update for next reducer search
    ap->P.SetLmCurrRing();
    ap->P.SetShortExpVector();
    ap->not_sev = ~ap->P.sev;
  }
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

  int pos = posInS(strat, strat->sl, ap->P.p, ap->P.ecart);

  strat->redTailChange = FALSE;

  if (rField_is_Z(currRing) && !rHasLocalOrMixedOrdering(currRing))
    redtailBbaAlsoLC_Z(&(ap->P), strat->tl, strat);

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
      superenterpairs(ap->P.p, strat->sl, ap->P.ecart, pos, strat, strat->tl);
    else
      enterpairs(ap->P.p, strat->sl, ap->P.ecart, pos, strat, strat->tl);

    strat->enterS(ap->P, pos, strat, strat->tl);
  }

  kDeleteLcm(&ap->P);
  ctx->stat_survivors.fetch_add(1, std::memory_order_relaxed);

  ap->P.Init();
  ap->occupied = false;
  ap->is_survivor = false;
}

/* ------------------------------------------------------------------ */
/*  Batch reduce function: each thread grabs and reduces polys from L  */
/* ------------------------------------------------------------------ */

static void batch_reduce(SweepContext *ctx, int thread_id)
{
  ActivePoly *ap = &ctx->active[thread_id];

  while (true)
  {
    // Grab a polynomial from L (under lock)
    pthread_mutex_lock(&ctx->L_lock);
    BOOLEAN got = pop_and_prepare(ctx, ap);
    pthread_mutex_unlock(&ctx->L_lock);

    if (!got) break;

    // Reduce it fully
    reduce_fully(ctx, ap);
    ctx->stat_rounds.fetch_add(1, std::memory_order_relaxed);

    // If survivor, stop — can't overwrite our slot
    if (ap->is_survivor) break;
  }
}

/* ------------------------------------------------------------------ */
/*  Worker thread function                                             */
/* ------------------------------------------------------------------ */

struct WorkerArg
{
  SweepContext *ctx;
  int thread_id;
};

static std::atomic<bool> g_done(false);

/**
 * Worker thread. Each batch:
 *   1. B0: barrier (wait for main to signal start of batch)
 *   2. batch_reduce: grab polys from L, reduce them
 *   3. B1: barrier (all done reducing, main processes survivors)
 * Workers check g_done after B0 and exit if true.
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
    pthread_barrier_wait(&ctx->barrier);  // B0: start of batch
    if (g_done.load(std::memory_order_acquire)) break;

    batch_reduce(ctx, thread_id);
    pthread_barrier_wait(&ctx->barrier);  // B1: all done reducing

    // Main processes survivors between B1 and B2
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

  g_done.store(false, std::memory_order_release);

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
  for (int j = 0; j <= strat->tl; j++)
  {
    if (strat->T[j].pLength <= 0)
      strat->T[j].pLength = pLength(strat->T[j].p ? strat->T[j].p : strat->T[j].t_p);
  }

  // Wait for all workers to start before entering main loop
  if (ctx->num_workers > 0)
    pthread_barrier_wait(&ctx->startup_barrier);

  // Main loop: batch reduce + process survivors
  //
  // Barrier protocol (3-phase per batch):
  //   B0: main signals start of batch (workers waiting here)
  //   batch_reduce: all threads reduce in parallel
  //   B1: all done reducing
  //   main processes survivors (workers idle)
  //   back to B0 for next batch
  //
  // Workers: B0 -> batch_reduce -> B1 -> (loop back to B0)
  // Main:    B0 -> batch_reduce -> B1 -> process_survivors -> (loop)
  //
  while (!strat->L.empty())
  {
    if (siCntrlc)
    {
      while (!strat->L.empty()) strat->L.pop_and_erase();
      strat->noClearS = TRUE;
      // Workers are waiting at B0; signal done and release them
      if (ctx->num_workers > 0)
      {
        g_done.store(true, std::memory_order_release);
        pthread_barrier_wait(&ctx->barrier);  // B0: workers check g_done and exit
      }
      goto parallel_cleanup;
    }

    // B0: signal start of batch (workers are waiting here)
    if (ctx->num_workers > 0)
      pthread_barrier_wait(&ctx->barrier);  // B0

    // Reduce: all threads participate in parallel
    batch_reduce(ctx, 0);
    pthread_barrier_wait(&ctx->barrier);  // B1: all done reducing

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
      for (int j = 0; j <= strat->tl; j++)
      {
        if (strat->T[j].pLength <= 0)
          strat->T[j].pLength = pLength(strat->T[j].p ? strat->T[j].p : strat->T[j].t_p);
      }
    }

    if (strat->overflow || errorreported)
    {
      // Drain L so workers find nothing in next batch
      while (!strat->L.empty()) strat->L.pop_and_erase();
      break;
    }
  }

  // Signal workers to exit: hit B0 with g_done=true so they break
  if (ctx->num_workers > 0)
  {
    g_done.store(true, std::memory_order_release);
    pthread_barrier_wait(&ctx->barrier);  // B0: workers check g_done and exit
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
