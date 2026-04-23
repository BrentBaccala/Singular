/**
 * @file kevlog.h
 * @brief Global event log for parallel bba correctness auditing
 *        (task 325 parallel-bba-event-log).
 *
 * Records every algorithmic decision (pop, reduce, enterT, enterS,
 * enterpairs start/end, LINSERT, KILL, LKILL, etc.) with full pointers
 * + a global seq number into a preallocated, lock-free append-only
 * buffer.  A post-hoc Python checker replays the log and identifies
 * the first illegal decision.
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
 * Aux arena: 256 MiB bump allocator for variadic payloads (LM
 * strings, T snapshots, etc.).  Records reference it by offset.
 */

#ifndef KEVLOG_H
#define KEVLOG_H

#include <cstdint>
#include <cstddef>
#include <atomic>

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
};

/* kill_reason / kill_src enums.  Small ints for the arg_b field of
 * KILL / LKILL / BKILL events.  Keep values stable — checker reads
 * them. */
enum kt_kill_reason : uint32_t {
  KR_NONE          = 0,
  KR_PROD_CRIT     = 1,   // lm(p)*lm(q) == lcm => product criterion
  KR_DOM_BY_B      = 2,   // LCM dominated by existing B entry
  KR_FROM_T_ECART  = 3,   // fromT + ecart rule
  KR_PCMP_CHAIN_EQ = 4,   // pCompareChain equal ecart/idx kill
  KR_PCMP_CHAIN_LT = 5,
  KR_PCMP_CHAIN_GT = 6,
  KR_BVEC_TRI_BVEC = 7,   // bvec_triangle B-side erase
  KR_BVEC_TRI_L    = 8,   // bvec_triangle L-side erase
  KR_LOCAL_HITS    = 9,   // chainCrit via local_hits
  KR_S_PAIRTEST    = 10,  // chainCrit via S_pairtest
  KR_GM            = 11,  // Gebauer-Moller
  KR_OTHER         = 99,
};

enum kt_reduce_outcome : uint32_t {
  RO_NONE     = 0,
  RO_SURVIVOR = 1,
  RO_ZERO     = 2,
  RO_CONTINUE = 3,   // re-swept for another reduce step
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
  uint64_t poly_ptr_1;   // 32 : raw poly pointer (identity across events)
  uint64_t poly_ptr_2;   // 40 : raw poly pointer
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
 *  Frees the buffer and aux arena.  Safe to call if init failed. */
void kevlog_shutdown();

/** Dump the buffer to disk.  Called on bad runs from the
 *  SINGULAR_CHECK_IDEAL_MEMBERSHIP violation handler.  disp_id, pid,
 *  and timestamp distinguish runs.  Writes:
 *    /tmp/audit-run/event-log-<disp>-<pid>-<ts>.bin
 *    /tmp/audit-run/event-log-<disp>-<pid>-<ts>-aux.bin
 *    /tmp/audit-run/event-log-<disp>-<pid>-<ts>-polys.txt
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
 *  Thread-safety: lock-free; safe from any worker thread.
 */
void kevlog_emit(uint16_t type, uint16_t tid, uint16_t atT, uint16_t flags,
                 uint32_t arg_a, uint32_t arg_b, uint32_t arg_c,
                 uint32_t aux_off,
                 const void *poly_ptr_1, const void *poly_ptr_2);

/** Allocate `n` bytes in the aux arena, return the offset.  On
 *  overflow: abort.  Thread-safe (atomic bump). */
uint32_t kevlog_aux_alloc(size_t n, void **out_ptr);

/** Register a poly pointer so that it ends up in the polys.txt dump
 *  at dump time.  Thread-safe.  Idempotent per pointer (deduped by
 *  a mutex-protected set).  It's fine to register pointers that are
 *  later freed — we only attempt p_String on pointers that survive
 *  (see defer-frees in kevlog.cc). */
void kevlog_register_poly(const void *p);

/* ------------------------------------------------------------------ */
/*  Leak-on-trace — simplest possible scheme for keeping poly pointers */
/*  captured in events valid at dump time.  The kt_pLmFree / kt_pDelete */
/*  / kt_p_LmFree / kt_p_Delete wrappers turn into no-ops when         */
/*  g_defer_frees is true (set by kevlog_init).  Process exit reclaims */
/*  leaked memory — fine for the 0.2s reproducer this is scoped to.    */
/* ------------------------------------------------------------------ */
extern bool g_defer_frees;

#ifdef __cplusplus
}
#endif

#endif /* KEVLOG_H */
