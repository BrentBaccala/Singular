/**
 * @file kthread.h
 * @brief Parallel Groebner basis reduction — parallel sweep design.
 *
 * For each polynomial from L, all threads cooperatively scan T[0..tl]
 * to find the best reducer (replacing kFindDivisibleByInT_ecart).
 * The main thread then performs ksReducePoly serially. Per-thread
 * sweep results avoid contention; a shared atomic cursor distributes
 * T positions across threads.
 *
 * See ~/project/docs/multithreading-plan.md for full design.
 */

#ifndef KTHREAD_H
#define KTHREAD_H

#include "kernel/GBEngine/kutil.h"
#include <pthread.h>
#include <atomic>

/**
 * Per-thread, per-slot sweep result — avoids contention on shared fields.
 */
struct SweepResult
{
  int best_reducer;   // first divisor found (any ecart), -1 = none
  int best_good;      // best divisor with ecart <= P.ecart, -1 = none
  int best_pLength;   // pLength of best_good
};

/**
 * An active polynomial slot in the parallel sweep.
 */
struct ActivePoly
{
  LObject P;                // the polynomial being reduced
  bool occupied;            // is this slot in use?
  unsigned long not_sev;    // ~P.sev for the sweep filter
  int best_reducer;         // merged: best T[j] found (any divisor), -1 = none
  int best_good;            // merged: best T[j] with ecart <= P.ecart, -1 = none
  int best_pLength;         // pLength of best_good reducer
  bool is_survivor;         // true if sweep found no divisor (true survivor)
  long d;                   // ecart tracking: h_d + h->ecart
  long reddeg;              // degree tracking from redHoney
  int pass;                 // reduction pass count
};

/**
 * Context for the parallel sweep.
 */
struct SweepContext
{
  ActivePoly *active;       // array of active polynomial slots
  int max_active;           // size of active array (2 * num_threads)
  int num_threads;          // total thread count (including main)
  int num_workers;          // worker threads (num_threads - 1, or from env)
  kStrategy strat;          // the Buchberger strategy
  ring r;                   // the ring (for setting currRing in workers)
  unsigned saved_si_opt_1;  // si_opt_1 for worker threads
  unsigned saved_si_opt_2;  // si_opt_2 for worker threads (PROT suppressed)

  // Per-thread sweep results: sweep_results[thread_id][slot]
  SweepResult *sweep_results;  // flat array [num_threads * max_active]

  // Shared atomic cursor for Phase 2 sweep
  std::atomic<int> sweep_cursor;

  // Shared atomic slot counter for Phase 1 and Phase 3
  std::atomic<int> slot_counter;

  // Barriers for phase synchronization (separate to prevent reuse races)
  pthread_barrier_t barrier_B0;   // start of batch
  pthread_barrier_t barrier_B1;   // end of batch (all done reducing)

  // Startup barrier: ensures all workers are running before main loop
  pthread_barrier_t startup_barrier;

  // Thread handles and IDs
  pthread_t *threads;
  int *thread_ids;          // thread_id for each worker

  // Mutex for L access (multiple threads pop from L)
  pthread_mutex_t L_lock;

  // Shutdown flag (per-context, not global static)
  std::atomic<bool> done;

  // Saved state for restoration after parallel phase
  int (*saved_posInT)(const BlockArray<TObject> &T, const int tl, LObject &h);
  BITSET saved_opt1;
  BOOLEAN withT;            // for redtailBba

  // Statistics
  std::atomic<long> stat_reductions;
  std::atomic<long> stat_zeros;
  std::atomic<long> stat_survivors;
  std::atomic<long> stat_rounds;
};

/* Public API */
int get_singular_threads();
SweepContext *sweep_context_init(kStrategy strat, int nthreads);
void sweep_context_destroy(SweepContext *ctx);
void bba_parallel_loop(SweepContext *ctx);

#endif /* KTHREAD_H */
