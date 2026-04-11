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
#include <deque>
#include <vector>

/* ------------------------------------------------------------------ */
/*  Instrumentation (task 482)                                         */
/*  Compile-time flag KTHREAD_INSTRUMENT enables measurement code.     */
/*  At runtime, SINGULAR_KTHREAD_STATS=1 must be set to actually       */
/*  emit output. When the flag is off, all instrumentation is stubbed  */
/*  out at zero overhead.                                              */
/* ------------------------------------------------------------------ */
#ifdef KTHREAD_INSTRUMENT
#include <time.h>
static inline long kt_now_ns()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (long)ts.tv_sec * 1000000000L + (long)ts.tv_nsec;
}

/**
 * Per-thread instrumentation accumulators. One per thread, zero
 * locking needed since each thread writes only its own slot.
 */
struct ThreadStats
{
  long wait_B0_ns;          // cumulative ns waiting at B0 barrier
  long wait_B1_ns;          // cumulative ns waiting at B1 barrier
  long wait_B0_count;       // number of B0 waits
  long wait_B1_count;       // number of B1 waits
  long wait_B0_max_ns;      // single longest B0 wait
  long wait_B1_max_ns;      // single longest B1 wait

  long sweep_ns;            // cumulative ns inside sweep_phase
  long sweep_count;         // number of sweep_phase calls
  long reduce_ns;           // cumulative ns inside reduce_slot_from_sweep
  long reduce_count;        // number of reduce_slot_from_sweep calls

  long drain_ns;            // cumulative ns inside drain_survivor_queue_locked
  long drain_count;         // number of drain calls (not survivors)
  long drain_survivors;     // total survivors processed on this thread

  long L_lock_wait_ns;      // cumulative ns waiting for L_lock
  long L_lock_count;        // number of L_lock acquisitions

  long surv_q_wait_ns;      // cumulative ns waiting for survivor_queue_mutex
  long surv_q_count;        // survivor queue mutex acquisitions

  long enterpairs_trylock_count;   // number of trylock attempts
  long enterpairs_trylock_fail;    // number of failed trylocks

  // process_survivor internal breakdown
  long ps_redtail_ns;       // ns inside redtailBba
  long ps_enterT_ns;        // ns inside enterT
  long ps_enterpairs_ns;    // ns inside enterpairs/superenterpairs
  long ps_enterS_ns;        // ns inside strat->enterS
  long ps_other_ns;         // remainder of process_survivor

  long round_start_ns;      // timestamp at start of current round
};

/**
 * Per-round event record.
 */
struct RoundRecord
{
  long round_id;
  long timestamp_ns;
  int queue_depth_start;    // depth when round began
  int queue_depth_end;      // depth after survivors added at round end
  int active_slots;         // active slots at round start
  int reductions;           // reductions applied this round
  long min_reduce_ns;       // min single-reduce duration this round
  long max_reduce_ns;       // max single-reduce duration this round
  long sum_reduce_ns;       // sum of all reduce durations
  long sweep_ns_main;       // main thread sweep duration this round
  long round_total_ns;      // wall-clock duration of this round
};

/**
 * Per-reduce event (transient, only collected per-round and summarized).
 */
struct ReduceEvent
{
  int thread_id;
  long start_ns;
  long duration_ns;
};
#endif  // KTHREAD_INSTRUMENT

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

  // ---- Asynchronous survivor enterpairs drain ---------------------
  // FIFO of survivor LObjects waiting for enterT/enterpairs/enterS
  // processing. Pushed by main thread at end of round; drained by a
  // single worker thread holding enterpairs_mutex.
  pthread_mutex_t survivor_queue_mutex;
  std::deque<LObject> *survivor_queue;

  // Serializes any drainer of survivor_queue (at most one thread at
  // a time runs process_survivor). Acquired by workers via trylock.
  // Main thread never acquires this (see kthread.cc comment).
  pthread_mutex_t enterpairs_mutex;

  // Signaled when enterpairs has added new entries to L, or when
  // the drain worker finishes draining. Used by the main thread to
  // wait for work while enterpairs is in progress. Protected by
  // L_lock.
  pthread_cond_t pairs_available;

  // True while a worker holds enterpairs_mutex and is draining.
  std::atomic<bool> enterpairs_active;

  // Max survivor queue depth observed (for diagnostics).
  std::atomic<long> stat_max_queue_depth;

#ifdef KTHREAD_INSTRUMENT
  // ---- Instrumentation (task 482) ---------------------------------
  bool stats_enabled;           // runtime toggle (SINGULAR_KTHREAD_STATS)
  ThreadStats *tstats;          // [num_workers+1]
  std::vector<RoundRecord> *rounds;
  std::vector<ReduceEvent> *reduces_this_round;
  pthread_mutex_t stats_lock;   // protects rounds / reduces_this_round
  long start_ns;                // t0 of bba_parallel_loop
  const char *workload_tag;     // optional label
#endif
};

/* Public API */
int get_singular_threads();
SweepContext *sweep_context_init(kStrategy strat, int nthreads);
void sweep_context_destroy(SweepContext *ctx);
void bba_parallel_loop(SweepContext *ctx);

#endif /* KTHREAD_H */
