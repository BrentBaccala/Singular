/**
 * @file kevlog.h
 * @brief Global event log for parallel bba correctness auditing
 *        (task 325 parallel-bba-event-log, refactored by task
 *        parallel-bba-event-log-pcopy).
 *
 * Records every algorithmic decision (pop, reduce, enterT, enterS,
 * enterpairs start/end, LINSERT, KILL, LKILL, etc.) with poly
 * identity + a global seq number into a preallocated, lock-free
 * append-only buffer.
 *
 * Gate: SINGULAR_EVENT_LOG=1.  Without it, allocation and hot-path
 * cost are zero.  Combined with SINGULAR_CHECK_IDEAL_MEMBERSHIP=1,
 * only failed runs dump the log to disk.
 *
 * Record: fixed 48-byte packed struct, little-endian native layout.
 * Capacity: 1 GiB / 48 B = ~22 million events.
 * On overflow: abort().  Do NOT silently drop or wrap — the checker
 * relies on a complete log.
 *
 * Aux arena: 256 MiB bump allocator for variadic payloads (T
 * snapshots, etc.).  Records reference it by offset.
 *
 * --- Poly capture semantics (post-refactor) --------------------
 *
 * Instead of storing raw in-flight poly pointers, every captured
 * poly is `p_Copy`d at emit time and the *copy's* pointer is stored
 * in the event.  To preserve pointer identity across events (so
 * "same source poly" still compares equal), a thread-safe map from
 * source-pointer -> copy-pointer is consulted: the first sight of
 * a source yields a fresh copy, subsequent sights return the same
 * copy pointer.
 *
 * At shutdown, every captured (copy) pointer is `p_Delete`d and the
 * map is cleared.  The dump's `polys.txt` renders LMs from the
 * copies — they were not mutated after creation, so the LMs are
 * authoritative.
 *
 * This works in this build because omalloc is configured to be a
 * thin wrapper over glibc malloc/free (see `build/omalloc/_config.h`
 * OMALLOC_USES_MALLOC=1 and `omalloc/xalloc.h` with XALLOC_BIN
 * commented out for thread safety).  All p_Copy / p_Delete calls
 * therefore bottom out in `malloc` / `free`, which are thread-safe.
 */

#ifndef KEVLOG_H
#define KEVLOG_H

#include <cstdint>
#include <cstddef>
#include <atomic>

/* Forward-declare the poly/ring C structs used by the capture API
 * without dragging in the large polys.h header. */
struct spolyrec;
struct ip_sring;

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/*  Gate state (plain-bool load in the hot path when disabled).        */
/* ------------------------------------------------------------------ */
extern bool g_event_log_enabled;

/* ------------------------------------------------------------------ */
/*  Event type enum                                                    */
/* ------------------------------------------------------------------ */
enum kt_evt_type : uint16_t {
  EVT_NONE           = 0,
  EVT_BPL_START      = 1,
  EVT_BPL_END        = 2,
  EVT_POP            = 3,
  EVT_REDUCE_START   = 4,
  EVT_REDUCE_STEP    = 5,
  EVT_REDUCE_END     = 6,
  EVT_REDTAIL        = 7,
  EVT_ENTERT         = 8,
  EVT_ENTERS         = 9,
  EVT_ENTERPAIRS_START = 10,
  EVT_ENTERPAIRS_END = 11,
  EVT_ENTERPAIR      = 12,
  EVT_KEEP           = 13,
  EVT_KILL           = 14,
  EVT_LINSERT        = 15,
  EVT_LKILL          = 16,
  EVT_BKILL          = 17,
  EVT_CHAINKILL      = 18,
  EVT_MAIN_BREAK     = 19,
  EVT_CLEARS_TOMBSTONE = 20,
  EVT_ZERO_REDUCE    = 21,  // pair reduced to 0, not a survivor
  EVT_GMKILL         = 22,  // Gebauer-Moller dedup kill
  EVT_SWEEP_RESULT   = 23,  // per-slot cooperative-sweep merged result +
                            // per-T-entry accept/reject reasons in aux
                            // (V1: 16-byte aux rows).
  EVT_SWEEP_RESULT_V2 = 24, // Same shape, but 32-byte aux rows that
                            // additionally carry the sweep-time
                            // captured T[j].p pointer (raw, no
                            // p_Copy) and sevT[j] value.  Task
                            // parallel-bba-sweep-atomic-capture.
  EVT_SWEEP_RESULT_V3 = 25, // Same shape as V2, but 40-byte aux rows
                            // with an extra 8-byte tail column
                            // sev_computed_at_sweep =
                            // pGetShortExpVector(T[j].p) computed at
                            // sweep time.  Enables H1a (sev_sweep ==
                            // sev_computed) vs H1b (!=) sub-
                            // classification of stale-sev misses.
                            // Task parallel-bba-enterT-sev-guard.
  EVT_SWEEP_RESULT_V4 = 26, // Same V3 40-byte aux rows, but the aux
                            // payload now starts with a 32-byte
                            // header carrying the active-slot P-side
                            // state captured at sweep time:
                            //   u32 n_examined; u32 pad;
                            //   u64 apP_sev_stored;
                            //   u64 apP_sev_computed;
                            //   u64 apP_not_sev_stored;
                            // followed by n_examined * 40-byte V3
                            // rows.  Lets the checker classify H1a
                            // misses into H1a.1 (P.sev stale vs P.p),
                            // H1a.2 (not_sev stale vs P.sev), or
                            // H1a.3 (both consistent — MYSTERY).
                            // Task parallel-bba-apP-sev-probe.
};

/* ------------------------------------------------------------------ */
/*  kill_reason enum for EVT_KILL (arg_b).                             */
/*  ------------------------------------------------------------------ */
/*  The KILL event is emitted inside enterOnePairNormal on the three
 *  pair-entry-kill paths.  The KR_* names are the historical spelling
 *  (task 325).  The KILL_* aliases are the task-prompt spelling — same
 *  values; either name is accepted.  Keep values stable — checker
 *  reads them.
 */
enum kt_kill_reason : uint32_t {
  KR_NONE          = 0,
  KILL_UNKNOWN     = 0,
  KR_PROD_CRIT     = 1,   // lm(p)*lm(q) == lcm => product criterion
  KILL_PROD_CRIT   = 1,
  KR_DOM_BY_B      = 2,   // c3++ path (L-or-B dominated by existing B)
  KILL_DOM_BY_B    = 2,
  KR_FROM_T_ECART  = 3,   // fromT + ecart-too-big
  KILL_FROM_T_ECART= 3,
  KR_PCMP_CHAIN_EQ = 4,   // legacy (unused)
  KR_PCMP_CHAIN_LT = 5,
  KR_PCMP_CHAIN_GT = 6,
  KR_BVEC_TRI_BVEC = 7,
  KR_BVEC_TRI_L    = 8,
  KR_LOCAL_HITS    = 9,
  KR_S_PAIRTEST    = 10,
  KR_GM            = 11,
  KR_OTHER         = 99,
};

/* ------------------------------------------------------------------ */
/*  kill_src enum for EVT_LKILL / EVT_BKILL / EVT_CHAINKILL /          */
/*  EVT_GMKILL (arg_a).                                                */
/* ------------------------------------------------------------------ */
/*  Mirrors the textual `src=...` field the SINGULAR_TRACE_PAIRCRIT
 *  text log writes for each kill path.  Keep values stable — checker
 *  reads them.  Checker / inspector decodes arg_a through this enum.
 */
enum kt_kill_src : uint32_t {
  LKILL_UNKNOWN                 = 0,
  LKILL_BVEC_TRIANGLE_BVEC_ERASE= 1, // bvec dedup: erase bvec[ii]
  LKILL_BVEC_TRIANGLE_L_ERASE   = 2, // bvec dedup: erase lt in L
  LKILL_PCOMPARE_CHAIN_SUGAR_GM = 3, // Gebauer+sugarCrit pCompareChain
  LKILL_PCOMPARE_CHAIN_GM       = 4, // Gebauer !sugarCrit pCompareChain
  LKILL_PCOMPARE_CHAIN_NONGEBAUER = 5, // non-Gebauer pCompareChain
  BKILL_PAIRTEST                = 10, // BKILL: dom-by-B L-side (sugar + nonsugar)
  BKILL_DOMBYB                  = 11, // BKILL: dom-by-B sugar, lp_lcm<it
  CHAINKILL_LOCAL_HITS          = 20, // local_hits chain-crit path
  CHAINKILL_S_PAIRTEST          = 21, // legacy S_pairtest chain-crit path
  CHAINKILL_LOCAL_HITS_LP       = 22,
  CHAINKILL_S_PAIRTEST_LP       = 23,
  GMKILL_B_GEBAUER              = 30, // the Gebauer-classic branch (if any)
  GMKILL_LFIRSTGM_SUGAR         = 31, // sugar_GM in chainCritNormal
  GMKILL_LFIRSTGM_NONSUGAR      = 32, // GM_tiebreak in chainCritNormal
  GMKILL_SUGAR_GM_ECART         = 33, // ecart-tiebreak in sugarCrit GM
};

enum kt_reduce_outcome : uint32_t {
  RO_NONE     = 0,
  RO_SURVIVOR = 1,
  RO_ZERO     = 2,
  RO_CONTINUE = 3,   // re-swept for another reduce step
};

/* ------------------------------------------------------------------ */
/*  sweep_reject_reason enum for EVT_SWEEP_RESULT aux payload          */
/*  (task parallel-bba-event-log-sweep).                               */
/*  ------------------------------------------------------------------ */
/*  For every T[j] entry the cooperative sweep LOOKED AT (via the      */
/*  K-strided per-tile loop in sweep_one_tile), the aux payload        */
/*  records one row { T_idx, reject_reason, entry_poly_ptr }.          */
/*  SWEEP_ACCEPTED marks the row(s) that contributed to best_reducer   */
/*  or best_good for this slot (after merge).                          */
/*                                                                     */
/*  Keep values stable — checker reads them.                            */
enum kt_sweep_reject_reason : uint32_t {
  SWEEP_ACCEPTED          = 0,
  SWEEP_REJECT_SEV_FILTER = 1,  // sevT[j] & not_sev_s != 0
  SWEEP_REJECT_NOT_PUBLISHED = 2, // tobject_published_load was false
  SWEEP_REJECT_NOT_DIVISIBLE = 3, // p_LmDivisibleBy returned false
  SWEEP_REJECT_ECART      = 4,  // divided, but ecart > P.ecart (not a
                                // "good" reducer — still makes it a
                                // best_reducer fallback)
  SWEEP_REJECT_WORSE_PLEN = 5,  // divided with acceptable ecart but a
                                // better (shorter) best_good already
                                // recorded
  SWEEP_REJECT_TOMBSTONED = 6,  // reserved — T entry marked dead
                                // (current sweep has no explicit
                                // tombstone gate; left for future use)
  SWEEP_REJECT_SELF       = 7,  // reserved — T[j] is the slot's own P
                                // (current sweep doesn't special-case
                                // this; left for future use)
  SWEEP_REJECT_FROM_T_RULE = 8, // reserved — fromT path rejection
  // Atomic-capture distinct outcomes (observed at sweep time, before
  // merge; written by sweep_one_tile into the per-thread capture
  // buffer; drained at close_slot into V2 aux rows).  Used to
  // distinguish sweep-time outcome from rescan-time outcome in the
  // check-event-log.py H1/H2 classifier.
  SWEEP_ACCEPTED_ATOMIC   = 10, // Passed sev + published + divides
                                // + ecart <= P.ecart at sweep time
                                // (before merge picked a winner)
  SWEEP_ATOMIC_REJECT_SEV = 11, // sev_j & not_sev_s != 0
  SWEEP_ATOMIC_REJECT_NOT_PUBLISHED = 12,
  SWEEP_ATOMIC_REJECT_NOT_DIVISIBLE = 13,
  SWEEP_ATOMIC_REJECT_ECART = 14,
  SWEEP_REJECT_OTHER      = 99, // catch-all
};

/* ------------------------------------------------------------------ */
/*  Record schema — 48 bytes, packed, native little-endian             */
/* ------------------------------------------------------------------ */
#pragma pack(push, 1)
struct KEvtRecord {
  uint64_t seq;          //  0 : redundant with position (sanity)
  uint16_t type;         //  8 : kt_evt_type
  uint16_t tid;          // 10 : thread id at event time
  uint16_t atT;          // 12 : strat->T.size() at event time
  uint16_t flags;        // 14 : reserved flag bits
  uint32_t arg_a;        // 16 : event-specific (i_r1, T_idx, arrival_id, ...)
  uint32_t arg_b;        // 20 : event-specific (i_r2, killer_idx, ...)
  uint32_t arg_c;        // 24 : event-specific (|L|, ecart, ...)
  uint32_t aux_off;      // 28 : offset into aux bump buffer (0 if none)
  uint64_t poly_ptr_1;   // 32 : captured (p_Copy'd) poly pointer
  uint64_t poly_ptr_2;   // 40 : captured (p_Copy'd) poly pointer
};
#pragma pack(pop)

static_assert(sizeof(KEvtRecord) == 48, "KEvtRecord must be 48 bytes");

/* ------------------------------------------------------------------ */
/*  Lifecycle                                                          */
/* ------------------------------------------------------------------ */
/** Called by bba() at entry.  Checks SINGULAR_EVENT_LOG env.  If set,
 *  allocates the 1 GiB buffer + 256 MiB aux arena.  Idempotent. */
void kevlog_init(int disp_id);

/** Called by bba() at exit (whether the run was good or bad).
 *  Frees the buffer and aux arena AND deletes every captured (copied)
 *  poly.  Safe to call if init failed. */
void kevlog_shutdown();

/** Dump the buffer to disk.  Called on bad runs from the
 *  SINGULAR_CHECK_IDEAL_MEMBERSHIP violation handler.  disp_id, pid,
 *  and timestamp distinguish runs.  Writes:
 *    /tmp/audit-run/event-log-<disp>-<pid>-<ts>.bin
 *    /tmp/audit-run/event-log-<disp>-<pid>-<ts>-aux.bin
 *    /tmp/audit-run/event-log-<disp>-<pid>-<ts>-polys.txt
 *
 *  Must be called BEFORE kevlog_shutdown (the dump needs the copies
 *  alive to format LMs).
 */
void kevlog_dump_on_failure(int disp_id);

/* ------------------------------------------------------------------ */
/*  Event emission                                                     */
/* ------------------------------------------------------------------ */
/** Reserve a slot (atomic fetch_add on the global seq), memcpy the
 *  record.  On overflow: abort with a clear stderr message.  Callers
 *  should gate on g_event_log_enabled to avoid argument evaluation
 *  cost when disabled.
 *
 *  poly_ptr_1 / poly_ptr_2 are the *captured* pointers (usually
 *  obtained via kevlog_capture / kevlog_capture_with_tail).  Passing
 *  an uncaptured raw pointer is legal but its LM won't be rendered
 *  in the dump (nothing to p_Delete either).
 *
 *  Thread-safety: lock-free; safe from any worker thread.
 */
void kevlog_emit(uint16_t type, uint16_t tid, uint16_t atT, uint16_t flags,
                 uint32_t arg_a, uint32_t arg_b, uint32_t arg_c,
                 uint32_t aux_off,
                 const void *poly_ptr_1, const void *poly_ptr_2);

/** Allocate `n` bytes in the aux arena, return the offset.  On
 *  overflow: abort.  Thread-safe (atomic bump). */
uint32_t kevlog_aux_alloc(size_t n, void **out_ptr);

/** Capture `src` by p_Copy'ing it into the log's poly registry.
 *  First sight of a source pointer triggers `p_Copy(src, r)`; later
 *  sights of the same source return the cached copy.  The returned
 *  pointer is stable for the duration of the traced bba run and is
 *  p_Delete'd at shutdown.
 *
 *  Returns NULL if `src` is NULL, if the log is disabled, or if
 *  p_Copy fails (shouldn't happen, but defensive).
 *
 *  Thread-safe: uses an internal mutex.  Low contention — only the
 *  first sight of each source allocates, subsequent sights are a
 *  lookup.
 */
const void *kevlog_capture(const void *src, struct ip_sring *r);

/** Same, but for polys that carry a separate tailRing (LObjects).
 *  Uses the `p_Copy(poly, lmRing, tailRing)` overload. */
const void *kevlog_capture_with_tail(const void *src,
                                     struct ip_sring *lmRing,
                                     struct ip_sring *tailRing);

/** Mark a (captured) poly as a final-S candidate — i.e. the poly_ptr_1
 *  of an EVT_ENTERS event.  The dump writes full `p_String(p)` for
 *  every poly so marked into the companion `-full-polys.txt` file.
 *  Thread-safe.  Idempotent (marking the same copy twice is cheap). */
void kevlog_mark_enters_poly(const void *copy_ptr);

#ifdef __cplusplus
}
#endif

#endif /* KEVLOG_H */
