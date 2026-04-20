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
#include <cstdint>
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

  // Phase-split instrumentation (task 506 enterpairs-parallel /
  // task 508 enterpairs-parallel-phase1).
  //   phase0: S-exclusive for setup + enterT + enterS
  //   phase1: S-shared + L-exclusive for enterpairs iteration +
  //           chainCrit (merges local B into strat->L) + clearS
  // As of task 508, phase 2 no longer exists (chainCritNormal does the
  // B-into-L merge inline under L-lock).  The phase2_* fields are
  // retained for ABI compatibility with the post-506 dump format but
  // renamed semantically.
  long phase0_wait_ns;      // blocked on S exclusive lock (phase 0)
  long phase0_ns;           // work inside phase 0
  long phase1_wait_ns;      // blocked on S shared lock (phase 1)
  long phase1_ns;           // iteration + chainCrit + clearS + B construction
  long phase1_l_wait_ns;    // blocked on L-exclusive during phase 1
                            // (formerly phase2_wait_ns)
  long phase1_l_ns;         // time under L-exclusive during phase 1
                            // (formerly phase2_ns)
  long phase_survivors;     // survivors passed through the phased path
  long phase_s_cas_fail;    // S tombstone CAS failures (peer drainer won)
  long phase_l_cas_fail;    // L tombstone CAS failures (peer drainer won)
  long phase1_concurrent_max;  // peak observed value of
                               // ctx->enterpairs_active during this
                               // thread's phase 1 (sampled at entry).

  // Worker-side drain participation (task 512 worker-side-drain).
  //
  // drain_survivors (above) already records the number of survivors this
  // thread processed through its drain call — tid > 0 values become
  // non-zero when worker-side drain is enabled.  The two fields below
  // capture how much time workers spend draining (as opposed to sweeping)
  // and how often they hop from tile-idle into drain.
  //
  //   worker_drain_idle_ns    : time spent in drain_survivor_queue calls
  //                             initiated from the tile-idle branch of
  //                             tile_pull_loop.  Includes the S-exclusive
  //                             wait in phase 0 if a peer drainer is
  //                             ahead.
  //   worker_drain_idle_count : number of such drain hops.
  long worker_drain_idle_ns;
  long worker_drain_idle_count;

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
 *
 * Milestone (d) of the continuous-cursor redesign (task 283): slots
 * have explicit state and generation counters so workers can tell
 * whether a slot they dispatched to via a tile batch is still live.
 */
struct ActivePoly
{
  LObject P;                // the polynomial being reduced
  bool occupied;            // is this slot in use? (main-only view)
  unsigned long not_sev;    // ~P.sev for the sweep filter
  int best_reducer;         // merged: best T[j] found (any divisor), -1 = none
  int best_good;            // merged: best T[j] with ecart <= P.ecart, -1 = none
  int best_pLength;         // pLength of best_good reducer
  bool is_survivor;         // true if sweep found no divisor (true survivor)
  long d;                   // ecart tracking: h_d + h->ecart
  long reddeg;              // degree tracking from redHoney
  int pass;                 // reduction pass count

  // T-snapshot bound (task 280). Captured on refill; used as the
  // per-slot upper bound on T indices in the sweep.
  int sl_snapshot;

  // Tile cursor infrastructure (task 281/282). tiles_remaining counts
  // down as worker-tiles sweep the slot; fetch_sub==1 is the closer
  // and runs reduce_slot_from_sweep.
  std::atomic<int> tiles_remaining;

  // Milestone (d) — slot state and generation.
  //
  //   state     : SLOT_EMPTY (main may refill) / SLOT_FILLED (worker
  //               may sweep). Transitions are release/acquire ordered.
  //   gen       : incremented every time the slot is re-published
  //               (fresh refill from L, or re-sweep after a reduce
  //               that kept the slot occupied). Workers read this
  //               against the tile batch's `gen` to discard tiles
  //               from stale publications (though the pipeline design
  //               avoids stale batches in practice).
  //   needs_republish: closer-reduce left the slot occupied-and-not-
  //               survivor — it needs another sweep pass. Main
  //               notices and publishes a fresh tile batch.
  std::atomic<int> state;
  std::atomic<uint64_t> gen;
  std::atomic<bool> needs_republish;
};

/**
 * Slot-state values.
 */
enum {
  SLOT_EMPTY  = 0,
  SLOT_FILLED = 1,
};

/**
 * A tile batch record: K tiles all targeting one (slot, gen). Main
 * stores one of these into the batches[] ring before publishing
 * tile_end; workers read it to decode the tile destination.
 */
struct TileBatch
{
  int slot;           // which ActivePoly slot
  uint64_t gen;       // expected generation of that slot
  int tl_snapshot;    // T upper bound to sweep to
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

  // Milestone (b) tile cursor (task 281).
  //
  // A tile is identified by an integer id in [0, tile_end). The
  // encoding is:
  //     slot  = id / tiles_K
  //     slice = id % tiles_K
  // so slices 0..K-1 of one slot are contiguous tile ids. This keeps
  // adjacent ids likely to hit the same slot (same LObject), which is
  // cache-friendly on the reducer side; slices within a slot step
  // through T with stride K (slice `i` visits j with j%K==i) which
  // still gives every thread a mix of small and large j values.
  //
  // Usage (this milestone): main thread sets tiles_K (= num_threads),
  // initializes each occupied slot's tiles_remaining to tiles_K, stores
  // 0 into tile_cursor, then publishes tile_end = max_active * tiles_K
  // with release semantics. Workers fetch_add on tile_cursor to claim
  // tiles. Each worker, after sweeping one tile, fetch_sub(1) on the
  // slot's tiles_remaining; the one that decrements to zero sets the
  // slot's `ready` flag. All this is still bracketed by B0/B1, so the
  // main thread can observe `ready` after B1 without further sync.
  std::atomic<uint64_t> tile_cursor;
  std::atomic<uint64_t> tile_end;
  int tiles_K;

  // Milestone (d) — tile-batch ring.
  //
  // Each tile id in [batch*K, (batch+1)*K) targets
  //   slot  = batches[batch % batch_ring_size].slot
  //   slice = id % K
  // batch_ring_size must be large enough that a batch's data is
  // never overwritten while any worker still holds a tile from it.
  // In practice the closer-reduce fires before any new batch for
  // the same slot is published, so a ring of size ~64 * pipeline
  // depth is plenty. publish_lock serialises main+closer writes
  // into this ring and the tile_end advance.
  TileBatch *batches;
  int batch_ring_size;
  std::atomic<uint64_t> batch_count;   // total batches published
  pthread_mutex_t publish_lock;
  pthread_cond_t  tiles_avail_cv;      // workers wait here for new tiles
  pthread_cond_t  slot_freed_cv;       // main waits here for slot/queue event

  // Shared atomic slot counter for Phase 1 and Phase 3
  std::atomic<int> slot_counter;

  // Barriers remain only for startup synchronisation. B0/B1 removed
  // in milestone (d) — main and workers now run continuously.
  pthread_barrier_t startup_barrier;

  // Pipeline depth: how many slots we allow to be filled at once.
  // Replaces max_active as the "target concurrency" knob.
  // Configurable via SINGULAR_PIPELINE_DEPTH (default 16 or num_threads,
  // whichever is larger). Note: max_active (number of slot records)
  // is set to pipeline_depth.
  int pipeline_depth;

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
  // processing. Pushed by main thread at end of round; drained by
  // any worker thread that observes the queue non-empty.
  //
  // Concurrency model (task 483, replacing task 275):
  //   - Queue push/pop is under survivor_queue_mutex (brief).
  //   - Each drain worker pops one survivor, then runs
  //     process_survivor_lobject while holding BOTH strat->S's mutex
  //     (via strat->S.lock()) and ctx->L_lock. This gives per-survivor
  //     mutual exclusion on S and L, not per-drain — task 275's
  //     "L_lock held across entire drain" is gone.
  //   - Multiple drain workers may be simultaneously popping from
  //     the queue and running process_survivor. The S+L lock pair
  //     serializes their inner work, but queue management and
  //     setup run in parallel, and the locks are released between
  //     survivors so fill_active_slots (which needs L_lock) can
  //     interleave.
  //   - Termination uses an atomic counter of in-flight drain
  //     workers (enterpairs_active) so the main thread can tell
  //     when all drainers have finished.
  pthread_mutex_t survivor_queue_mutex;
  std::deque<LObject> *survivor_queue;

  // Counter of workers currently inside process_survivor_lobject.
  // Incremented on drain entry, decremented on drain exit. Replaces
  // the old task-275 enterpairs_active boolean. Used by the main
  // thread to tell when the drain is quiescent at termination.
  std::atomic<int> enterpairs_active;

  // Signaled when a drain worker finishes a survivor (L may have
  // new entries). Kept for future CV-based wait loops; currently
  // only used for main-thread idle-drain termination check.
  pthread_cond_t pairs_available;

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

/* ------------------------------------------------------------------ */
/*  Debug breadcrumb API (SINGULAR_AUDIT_T).                           */
/*                                                                     */
/*  kt_debug_tid is thread-local and set by the parallel loop (main    */
/*  thread 0) and worker_thread at startup.  Files outside kthread.cc  */
/*  (e.g. kspoly.cc) can drop breadcrumbs into the shared event ring   */
/*  by calling kt_debug_tag().  A no-op when the audit is disabled.    */
/* ------------------------------------------------------------------ */
extern thread_local int kt_debug_tid;
void kt_debug_tag(const char *op, void *poly, int slot, int arg);

#endif /* KTHREAD_H */
