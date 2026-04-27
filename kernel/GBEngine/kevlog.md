# kevlog — Parallel bba Event Log

A lock-free, pid-local, dump-on-failure trace of every algorithmic
decision inside a parallel-bba run. Built across tasks 325-335 to
give the checker a complete ground truth for diagnosing
non-deterministic failures that `rr` can't reproduce.

Source: `kevlog.h` (schema + public API), `kevlog.cc` (emit / dump /
capture implementation). Event sites scattered across
`kthread.cc`, `kstd2.cc`, `kutil.cc`.

## What kevlog is (and isn't)

- **Is**: a single in-memory record buffer that every worker thread
  appends to lock-free during bba. Complete on both good and bad
  runs; written to disk only on bad.
- **Isn't**: a ring buffer — overflow aborts rather than wraps.
  A streaming log — disk I/O is failure-only. Cross-process —
  each Singular invocation owns its own buffer.

## Gating (two env vars)

- `SINGULAR_EVENT_LOG=1` — `kevlog_init(disp_id)` calloc's a
  1 GiB record buffer (~22M slots) plus a 256 MiB aux arena.
  `bool g_event_log_enabled` is the hot-path gate; when false,
  `kevlog_emit` is a single load + branch and all capture helpers
  return NULL.
- `SINGULAR_CHECK_IDEAL_MEMBERSHIP=1` — after bba returns, verify
  the result is correct; on mismatch call
  `kevlog_dump_on_failure(disp)` *before* `kevlog_shutdown()`.

Good runs emit to memory, then `kevlog_shutdown()` frees everything.
No disk I/O. Bad runs dump four files (see below), then shutdown.

## Record schema (48 bytes, packed, little-endian)

```
 0..7   seq (u64)        — redundant with position, sanity check
 8..9   type (u16)       — kt_evt_type
10..11  tid (u16)        — emitter thread id
12..13  atT (u16)        — strat->T.size() at event time
14..15  flags (u16)
16..19  arg_a (u32)      — event-specific
20..23  arg_b (u32)      — event-specific
24..27  arg_c (u32)      — event-specific
28..31  aux_off (u32)    — offset into aux arena; 0 = no payload
32..39  poly_ptr_1 (u64) — captured copy pointer
40..47  poly_ptr_2 (u64) — captured copy pointer
```

## Lock-free emission

```cpp
uint64_t my = g_seq.fetch_add(1);    // reserve a slot
if (my >= EVLOG_CAP_RECS) abort();    // never wrap
memcpy(&g_buf[my], &r, 48);           // fill in place
```

One CAS per event. Overflow aborts with a clear stderr message —
the checker relies on a complete record, so silent drop is worse
than crash.

## Aux arena (variadic payloads)

256 MiB bump allocator via `kevlog_aux_alloc(n, &out_ptr)` — atomic
fetch_add on `g_aux_off`, 8-byte aligned. Used for payloads that
don't fit the 48-byte record (T-snapshot at REDUCE_START, per-T-entry
rows at SWEEP_RESULT). `aux_off = 0` is reserved for "no payload";
offsets are otherwise valid indices into the arena.

## Poly capture (the subtle part)

Raw poly pointers are not content-stable — `ksReducePoly` mutates
exponent vectors in place (task 325/326 learned this the hard way
and spent a full task redesigning around it). Every poly that lands
in a `poly_ptr_{1,2}` or inside an aux payload goes through
`kevlog_capture(src, ring)`:

- Mutex-guarded append-only `vector<{copy, source, lmRing,
  tailRing}>`.
- Every call → fresh `p_Copy(src, ring)` → append one entry.
  Each capture is a **live snapshot** of the poly AT THE MOMENT
  of the call.
- At `kevlog_shutdown()`, every copy gets `p_Delete`d.

### Why no source→copy caching

Task 336 discovered that a first-sight-caching design (keyed by
source pointer, returning the same cached copy for all later
sights) produces stale LMs for any poly whose source is mutated
in place later. The concrete offender is `ap->P.p`: the same raw
pointer survives across successive `ksReducePoly` rounds carrying
progressively-reduced content. Cached captures froze the
*pre-mutation* snapshot forever, and the checker's LM-divides
invariant-7 logic mistakenly compared sweep-time live sev against
stale cached LM — producing a cascade of phantom sev-filter
false-rejects that drove tasks 328-336. Once captures became live
snapshots, every sev-filter "false-reject" disappeared.

### Cross-event "same source" comparisons

Losing the cache meant losing the free `copy_ptr == copy_ptr`
identity check between events. Consumers now look up each event's
`copy_ptr` in the 3-column `-polys.txt` to get its `source_ptr`,
and compare source pointers. The checker's `copy_to_src[copy] =
source` map makes this one line. See
`check_survivor_bookkeeping` for the pattern.

Thread-safety rests on this build's `OMALLOC_USES_MALLOC=1` and
`XALLOC_BIN` being commented out in `omalloc/xalloc.h` — so
`p_Copy` bottoms out in plain `malloc`/`free`, which is
thread-safe. Without that configuration, kevlog's p_Copy mutex
would be insufficient.

`kevlog_capture_with_tail(src, lmRing, tailRing)` is the LObject
variant for polys with a separate tailRing.

`kevlog_mark_enters_poly(copy_ptr)` flags a captured copy as a
final-S candidate; the dump writes its full `p_String` into
`-full-polys.txt` for Buchberger closure checking by the checker.

## Dump on failure (four files)

Path stem: `/tmp/audit-run/event-log-<disp>-<pid>-<ts>`.

1. **`.bin`** — 32-byte header (`KEVLOG01` magic + nrec + sizeof(rec)
   + disp/pid), followed by the raw record array.
2. **`-aux.bin`** — the aux arena contents. Offset 0 is the
   "no-aux" sentinel. Records reference this by `aux_off`.
3. **`-polys.txt`** — three columns:
   `copy_addr<TAB>lm<TAB>source_addr`. One row per capture (no
   source-to-copy dedup after the drop-the-cache change): repeated
   captures of the same source produce distinct `copy_addr` rows
   with potentially distinct LMs (e.g., after each in-place
   mutation), all sharing `source_addr`. Use `source_addr` for
   "same source poly across events" comparisons; use `copy_addr`
   for LM rendering of a specific event's snapshot. The checker
   loads both `polys_by_copy[copy] = lm` and `copy_to_src[copy]
   = source`.
4. **`-full-polys.txt`** — full `p_String(p)` for every poly
   marked via `kevlog_mark_enters_poly` (called at each
   `EVT_ENTERS` site). These are the final-S candidates needed
   by the Buchberger closure check in invariant 8.

## Events that carry aux payloads

Only two event types actually use the aux arena. All other events
emit with `aux_off = 0` and carry their state in the fixed
48-byte record.

### `EVT_REDUCE_START` (type 4) — T snapshot

Layout (variable length):

```
 0..3    u32 T_size
 4..     u64 poly_ptr [T_size]    — captured copies of strat->T[j].p
```

Records what `strat->T` looked like at the moment the reducer was
selected for this slot, so the checker can verify whether a sweep
that returned "no reducer" was really surveying a T that lacked a
divisor. The ptrs are `kevlog_capture`'d; LMs render via
`-polys.txt`. One payload per `EVT_REDUCE_START` event.

### `EVT_SWEEP_RESULT` (type 23) and its versioned descendants

One event per slot after the cooperative sweep merges its per-tile
results. The payload records per-T-entry attribution for the sweep
decision.

**Emit gate**: aux is filled only when `merged_best == -1` — i.e.,
the sweep produced a survivor. Accepted-reducer cases emit a bare
record with `aux_off = 0` since the attribution isn't interesting.
This gate, added in task 329, is the Heisenberg mitigation that
kept the bug rate in the 0.5-3% band despite heavy instrumentation.

Four versions have shipped, distinguished by the event type tag.
Older dumps remain parseable because each version uses a distinct
tag; the inspector/checker dispatch on type.

#### V1 (type 23) — task 329 (`parallel-bba-event-log-sweep`)

```
prelude:
 0..3    u32 n_examined
 4..7    u32 pad

row (16 bytes):
 0..3    u32 T_idx
 4..7    u32 reason         — kt_sweep_reject_reason
 8..15   u64 entry_ptr      — p_Copy of strat->T[j].p at emit time
```

Installed the per-T attribution framework. Introduced the reason
enum (`SWEEP_ACCEPTED`, `SWEEP_REJECT_SEV_FILTER`,
`SWEEP_REJECT_NOT_DIVISIBLE`, `SWEEP_REJECT_ECART`,
`SWEEP_REJECT_WORSE_PLEN`, …).

#### V2 (type 24) — task 332 (`parallel-bba-sweep-atomic-capture`)

Row widens to 32 bytes; prelude unchanged.

```
row (32 bytes):
 0..3    u32 T_idx
 4..7    u32 reason         — kt_sweep_reject_reason (now includes
                              SWEEP_ATOMIC_* variants — sweep-time
                              outcomes distinct from rescan-time)
 8..15   u64 entry_ptr      — p_Copy (rescan)
16..23   u64 sweep_ptr      — raw T[j].p read at sweep time
24..31   u64 sev_sweep      — sevT[j] read at sweep time
```

Added per-thread × per-slot sweep-time capture (no p_Copy in hot
path; raw pointer only). Enabled the H1 (stale sev) vs H2 (T[j]
replaced between sweep and rescan) classification in the checker.
Result: 6/6 H1 with hard pointer identity, H2 ruled out.

#### V3 (type 25) — task 334 (`parallel-bba-enterT-sev-guard`)

Row widens to 40 bytes; prelude unchanged.

```
row (40 bytes):
 0..3    u32 T_idx
 4..7    u32 reason
 8..15   u64 entry_ptr
16..23   u64 sweep_ptr
24..31   u64 sev_sweep
32..39   u64 sev_computed   — pGetShortExpVector(T[j].p) at
                              sweep time, fresh
```

Added a companion enterT publish guard (`kt_entert_sev_guard_hit`
→ `/tmp/audit-run/entert-sev-guard.log`) that checks
`sevT[atT] == pGetShortExpVector(T[atT].p)` immediately after
`tobject_publish`. The V3 row supports H1a (sev_sweep == sev_computed
→ sev consistent with T[j].p at sweep time) vs H1b (differ → torn
read or drift) sub-classification. Guard silent in 400 iters; V3
showed 14/14 H1a — T[j] side ruled out.

#### V4 (type 26) — task 335 (`parallel-bba-apP-sev-probe`) — current

Prelude extends to 32 bytes (replaces V1-V3's 8-byte prelude) by
carrying per-slot state for the active poly; rows unchanged from V3.

```
prelude (32 bytes):
 0..3    u32 n_examined
 4..7    u32 pad
 8..15   u64 apP_sev_stored        — ap->P.sev at sweep time
16..23   u64 apP_sev_computed      — pGetShortExpVector(ap->P.p)
                                     at sweep time, fresh
24..31   u64 apP_not_sev_stored    — ap->not_sev at sweep time

row (40 bytes, V3 layout):
 0..3    u32 T_idx
 4..7    u32 reason
 8..15   u64 entry_ptr
16..23   u64 sweep_ptr
24..31   u64 sev_sweep
32..39   u64 sev_computed
```

Added a companion apP publish guard (`kt_app_sev_guard_hit` →
`/tmp/audit-run/app-sev-guard.log`) that checks
`ap->P.sev == pGetShortExpVector(ap->P.p)` AND
`ap->not_sev == ~ap->P.sev` immediately after the release-store
in `publish_slot_tiles_locked`. The V4 slot block supports
H1a.1 / H1a.2 / H1a.3 sub-classification of H1a stale-sev misses.
Guard silent in 400 iters; V4 showed 4/4 H1a.3 (all invariants
self-consistent, yet the sev test still fires on an identity
divisor) — MYSTERY verdict that current work is discriminating.

## arg_a / arg_b / arg_c per event type

This section documents what each arg_a, arg_b, arg_c field carries
for every event type. Readers can use this to trace a pair's lifecycle
from POP through reduction and entry to the final basis, without
rereading the source code.

### EVT_BPL_START (type 1) — `kstd2.cc:2839`

Emitted once at bba() entry to mark the start of a parallel-bba run.

- `arg_a` = `get_singular_threads()` — number of worker threads configured
- `arg_b` = 0 (unused)
- `arg_c` = 0 (unused)
- `poly_ptr_1` = NULL
- `poly_ptr_2` = NULL

### EVT_BPL_END (type 2) — `kstd2.cc:3498`, `3511`, `3558`

Emitted at bba() exit or on violation (before dump). Multiple variants:
  - Normal failure detection (3498)
  - Force-dump on passing run (3511)
  - End of bba proper (3558)

- `arg_a` = 0 (unused)
- `arg_b` = 0 (unused)
- `arg_c` = 0 (unused)
- `poly_ptr_1` = NULL
- `poly_ptr_2` = NULL

### EVT_POP (type 3) — `kthread.cc:1673`

Emitted when a pair is popped from the L-queue for reduction.

- `arg_a` = `i_r1` — R-slot index of parent polynomial 1 (u32(-1) = sentinel for no T-slot; used for initial empty-L pre-pops that pull F polynomials directly)
- `arg_b` = `i_r2` — R-slot index of parent polynomial 2 (u32(-1) sentinel as for arg_a)
- `arg_c` = `strat->L.size() + 1` — L-queue size *before* this pop
- `poly_ptr_1` = `p_Copy(T[ir1].p)` — captured snapshot of parent 1 at pop time (or NULL if ir1=-1)
- `poly_ptr_2` = `p_Copy(T[ir2].p)` — captured snapshot of parent 2 at pop time (or NULL if ir2=-1)

### EVT_REDUCE_START (type 4) — `kthread.cc:2816`

Emitted when a slot's reducer is selected by cooperative sweep, before
reduction begins. Carries a T-snapshot in the aux arena showing which T
entries the sweep surveyed.

- `arg_a` = `slot` — active slot index
- `arg_b` = `best` — selected reducer T-slot index (or u32(-1) if no reducer found)
- `arg_c` = `Tsz` — number of T entries captured in aux
- `aux_off` ≠ 0 — offset into aux arena; payload is: `u32 T_size; u64 poly_ptr[T_size]`
- `poly_ptr_1` = `p_Copy(ap->P.p)` — captured snapshot of active poly before reduction
- `poly_ptr_2` = NULL

### EVT_REDUCE_STEP (type 5) — `kthread.cc:2913`

Emitted before each call to `ksReducePoly` in a reduction loop (may be
called multiple times per slot if the poly survives and needs re-sweeping).

- `arg_a` = `slot` — active slot index
- `arg_b` = `best` — reducer T-slot index
- `arg_c` = 0 (unused)
- `poly_ptr_1` = `p_Copy(input_poly)` — captured snapshot of poly *before* this reduction step
- `poly_ptr_2` = `p_Copy(reducer_poly)` — captured snapshot of the T-entry used to reduce

### EVT_REDUCE_END (type 6) — `kthread.cc:2834`, `2959`, `3067`

Emitted when reduction of a slot terminates. Multiple variants indicate
the outcome:
  - `2834`: survivor path (no reducer found, poly becomes final-S candidate)
  - `2959`: zero path (poly reduced to 0, not a survivor)
  - `3067`: continue path (poly survived this round, needs re-sweep)

- `arg_a` = `slot` — active slot index
- `arg_b` = outcome code (u32):
  - `RO_SURVIVOR = 1` — poly survived reduction
  - `RO_ZERO = 2` — poly reduced to zero
  - `RO_CONTINUE = 3` — poly reduced to nonzero but not survivor; will be re-swept
- `arg_c` = 0 (unused)
- `poly_ptr_1` = `p_Copy(ap->P.p)` — final poly after all reduction (or NULL if RO_ZERO)
- `poly_ptr_2` = NULL

### EVT_REDTAIL (type 7) — `kthread.cc:3274`, `3316`

Emitted when a final polynomial undergoes tail reduction (either Z-module
mode with `redTailZ` or standard mode with `redtailBba`).

- `arg_a` = mode code:
  - 0 = Z-module mode (`redTailZ`)
  - 1 = standard mode (`redtailBba`)
- `arg_b` = 0 (unused)
- `arg_c` = 0 (unused)
- `flags` carries `strat->redTailChange` in bit 0:
  - 0 = tail reduction made no change
  - 1 = tail reduction changed the poly
- `poly_ptr_1` = `p_Copy(rt_in)` — captured snapshot of poly *before* tail reduction
- `poly_ptr_2` = `p_Copy(P->p)` — captured snapshot of poly *after* tail reduction

### EVT_ENTERT (type 8) — `kthread.cc:3366`

Emitted when a final polynomial is entered into T (the reduction table).

- `arg_a` = new T-index (same as `strat->T.size()` at emit time)
- `arg_b` = `P->ecart` — Buchberger criterion parameter
- `arg_c` = 0 (unused)
- `poly_ptr_1` = `p_Copy(P->p)` — captured snapshot of the polynomial being added to T
- `poly_ptr_2` = NULL

### EVT_ENTERS (type 9) — `kthread.cc:3472`

Emitted when a final polynomial is entered into S (the Gröbner basis).
This marks a completed reduction chain from an initial basis element
through all reductions to final form.

- `arg_a` = `arrival_id` — global sequence number tracking when this poly was first created/inserted into L
- `arg_b` = `strat->S.size()` — S size *before* this entry
- `arg_c` = `P->ecart` — Buchberger criterion parameter
- `poly_ptr_1` = `p_Copy(P->p)` — captured snapshot of final polynomial (marked as "full-poly" candidate for Buchberger closure check)
- `poly_ptr_2` = NULL

**Note**: The captured poly is marked via `kevlog_mark_enters_poly` so the dump writes its full `p_String` into `-full-polys.txt` for the final checker's Buchberger closure invariant.

### EVT_ENTERPAIRS_START (type 10) — `kthread.cc:3627`

Bracket start: marks when a basis element enters pairing generation.
Called before `enterpairs` or `superenterpairs`.

- `arg_a` = `my_arrival` — arrival_id of the S entry being paired with all earlier entries
- `arg_b` = `strat->S.size()` — S size at entry
- `arg_c` = 0 (unused)
- `poly_ptr_1` = `p_Copy(P->p)` — captured snapshot of the pairing element
- `poly_ptr_2` = NULL

### EVT_ENTERPAIRS_END (type 11) — `kthread.cc:3640`

Bracket end: marks when pairing generation completes and all new pairs
have been merged into L by `chainCritNormal`.

- `arg_a` = `my_arrival` — arrival_id (same as START)
- `arg_b` = number of new pairs added to L (= L size after - L size before)
- `arg_c` = `strat->S.size()` — S size at completion
- `poly_ptr_1` = `p_Copy(P->p)` — captured snapshot of the pairing element (redactantly repeated)
- `poly_ptr_2` = NULL

### EVT_ENTERPAIR (type 12) — `kutil.cc:2285`

Emitted when a single pair (h, s[i]) is created and about to be tested
against existing kills.

- `arg_a` = `si.arrival_id` — arrival_id of the S entry s[i] being paired with h
- `arg_b` = 0 (unused)
- `arg_c` = `ecart` — Buchberger criterion parameter (initial ecart before spoly)
- `poly_ptr_1` = `p_Copy(h)` — captured snapshot of first pair element (from enterS)
- `poly_ptr_2` = `p_Copy(si.p)` — captured snapshot of second pair element (from S)

### EVT_KEEP (type 13) — `kutil.cc:2788`

Emitted when a pair survives all pre-insertion criteria and is about to
be added to B (the batch queue before L-merge).

- `arg_a` = `si.arrival_id` — arrival_id of the S entry
- `arg_b` = 0 (unused)
- `arg_c` = `Lp.ecart` — pair's spoly ecart after LCM computation
- `poly_ptr_1` = `p_Copy(p)` — captured snapshot of first pair element
- `poly_ptr_2` = `p_Copy(si.p)` — captured snapshot of second pair element

### EVT_KILL (type 14) — `kutil.cc:2309`, `2349`, `2376`, `2416`, `2491`, `2571`

Emitted when a pair is eliminated by one of three kill criteria during
pair creation. Multiple call sites with different kill reasons:
  - `2309`: `KR_FROM_T_ECART` (sugar mode, fromT path, ecart too big)
  - `2349`: `KR_PROD_CRIT` (sugar mode, product criterion)
  - `2376`: `KR_FROM_T_ECART` (nonsugar mode, fromT path)
  - `2416`: `KR_PROD_CRIT` (nonsugar mode, product criterion)
  - `2491`: `KR_DOM_BY_B` (pair dominated by existing B entry)

- `arg_a` = `si.arrival_id` — arrival_id of the killed S entry
- `arg_b` = `kill_reason` — code from `kt_kill_reason` enum:
  - `KR_PROD_CRIT = 1` — product criterion (lcm(p,q) == lm(p)*lm(q))
  - `KR_DOM_BY_B = 2` — dominated by existing B entry
  - `KR_FROM_T_ECART = 3` — fromT with ecart > current ecart
  - (others: `KR_PCMP_CHAIN_*`, `KR_BVEC_TRI_*`, `KR_LOCAL_HITS`, `KR_S_PAIRTEST`, `KR_GM`)
- `arg_c` = context-dependent (ecart or 0)
- `poly_ptr_1` = `p_Copy(p)` — captured snapshot of first pair element
- `poly_ptr_2` = `p_Copy(si.p)` — captured snapshot of second pair element

### EVT_LINSERT (type 15) — `kutil.cc:3867`, `3909`

Emitted when a pair is merged from B into L. Two variants:
  - `3867`: normal insert from B into L
  - `3909`: insert returned via iterator (flags=1 indicates `_ret` variant)

- `arg_a` = `Lobj.i_r1` — R-slot index of first parent
- `arg_b` = `Lobj.i_r2` — R-slot index of second parent
- `arg_c` = `Lobj.ecart` — pair's ecart
- `poly_ptr_1` = `p_Copy(Lobj.p1)` — captured snapshot of first parent
- `poly_ptr_2` = `p_Copy(Lobj.p2)` — captured snapshot of second parent

### EVT_LKILL (type 16) — `kutil.cc:4170`, `4339`, `4465`, `4524`, `4561`

Emitted when a pair in L is eliminated by chain criterion. Multiple
kill-src codes indicate the specific chain-criterion path:

- `arg_a` = `kill_src` — code from `kt_kill_src` enum:
  - `LKILL_PCOMPARE_CHAIN_SUGAR_GM = 3` — pCompareChain + sugar + Gebauer-Moller
  - `LKILL_PCOMPARE_CHAIN_GM = 4` — pCompareChain + nonsugar Gebauer-Moller
  - `LKILL_PCOMPARE_CHAIN_NONGEBAUER = 5` — non-Gebauer pCompareChain
  - (see `kevlog.h:209-226` for full enum)
- `arg_b` = u32(-1) (reserved, not used for LKILL)
- `arg_c` = `ecart` — pair's ecart at kill time
- `poly_ptr_1` = `p_Copy(killer)` — captured snapshot of the dominating element
- `poly_ptr_2` = `p_Copy(victim_lcm)` — captured snapshot of the L pair's lcm

### EVT_BKILL (type 17) — `kutil.cc:2519`, `2598`

Emitted when a pair in B is eliminated by comparison with another B
entry. Two variants indicate sugar vs. non-sugar mode:

- `arg_a` = `kill_src` — code from `kt_kill_src` enum:
  - `BKILL_PAIRTEST = 10` — dom-by-B via pairtest
  - `BKILL_DOMBYB = 11` — dom-by-B with sugar, lp_lcm < it
- `arg_b` = u32(-1) (reserved, not used for BKILL)
- `arg_c` = `Lp.ecart` — killer pair's ecart
- `poly_ptr_1` = `p_Copy(Lp.lcm)` — captured snapshot of killer pair's lcm
- `poly_ptr_2` = `p_Copy(it->lcm)` — captured snapshot of victim pair's lcm

### EVT_CHAINKILL (type 18) — `kutil.cc:3995`, `4033`, `4078`, `4122`

Emitted when a pair in B is eliminated by a chain-criterion element.
Multiple call sites indicate whether the pairtest used local_hits or
legacy S_pairtest, and sugar vs. non-sugar:

- `arg_a` = `kill_src` — code from `kt_kill_src` enum:
  - `CHAINKILL_LOCAL_HITS = 20` — local_hits chain-crit (sugar mode, short expo variant)
  - `CHAINKILL_S_PAIRTEST = 21` — S_pairtest legacy path
  - `CHAINKILL_LOCAL_HITS_LP = 22` — local_hits chain-crit (sugar mode, long poly variant)
  - `CHAINKILL_S_PAIRTEST_LP = 23` — S_pairtest legacy long-poly variant
- `arg_b` = `sit->arrival_id` — arrival_id of the killer S entry
- `arg_c` = 0 (unused)
- `poly_ptr_1` = `p_Copy(sit->p)` — captured snapshot of killer element (from S)
- `poly_ptr_2` = `p_Copy(it->lcm)` — captured snapshot of victim pair's lcm

### EVT_MAIN_BREAK (type 19) — `kthread.cc:4206`

Emitted once by the main thread when it detects that all work is
complete (L empty, active slots empty, all drainers idle).

- `arg_a` = `strat->S.size()` — final S size
- `arg_b` = `strat->L.size()` — final L size (should be 0)
- `arg_c` = 0 (unused)
- `poly_ptr_1` = NULL
- `poly_ptr_2` = NULL

### EVT_CLEARS_TOMBSTONE (type 20) — `kutil.cc:1260`

Emitted when a tombstoned S entry is physically erased from the basis set.
Tombstoning occurs during `enterpairs` walk to mark entries < arrival_id
as logically deleted without disrupting the iterator.

- `arg_a` = `i` — S array index of the tombstoned entry being erased
- `arg_b` = `victim.arrival_id` — arrival_id of the tombstoned entry
- `arg_c` = 0 (unused)
- `poly_ptr_1` = `p_Copy(victim)` — captured snapshot of the tombstoned (erased) entry
- `poly_ptr_2` = `p_Copy(p)` — captured snapshot of the killer/replacement poly

### EVT_ZERO_REDUCE (type 21) — `kthread.cc:2963`

Emitted when a pair reduces to zero (becomes a syzygy witness). Complements
the `EVT_REDUCE_END` with `RO_ZERO` by marking this outcome explicitly.

- `arg_a` = `slot` — active slot index
- `arg_b` = 0 (unused)
- `arg_c` = 0 (unused)
- `poly_ptr_1` = NULL
- `poly_ptr_2` = NULL

### EVT_GMKILL (type 22) — `kutil.cc:4220`, `4245`, `4273`, `4297`, `4389`, `4413`

Emitted when a pair in B is eliminated by Gebauer-Moller deduplication.
Multiple call sites indicate the ecart-comparison or tiebreak outcome:
  - `4220`, `4245`: `GMKILL_SUGAR_GM_ECART` (ecart-based selection)
  - `4273`, `4297`: `GMKILL_LFIRSTGM_SUGAR` (tiebreak selection in sugar mode)
  - `4389`, `4413`: `GMKILL_LFIRSTGM_NONSUGAR` (tiebreak selection in nonsugar mode)

- `arg_a` = `kill_src` — code from `kt_kill_src` enum:
  - `GMKILL_SUGAR_GM_ECART = 33` — ecart tiebreak
  - `GMKILL_LFIRSTGM_SUGAR = 31` — sugar mode tiebreak
  - `GMKILL_LFIRSTGM_NONSUGAR = 32` — nonsugar mode tiebreak
  - `GMKILL_B_GEBAUER = 30` — (legacy, used in some code paths)
- `arg_b` = u32(-1) (reserved, not used for GMKILL)
- `arg_c` = victim `ecart` (killer ecart when ecart-based selection, 0 otherwise)
- `poly_ptr_1` = `p_Copy(killer_lcm)` — captured snapshot of killer pair's lcm
- `poly_ptr_2` = `p_Copy(victim_lcm)` — captured snapshot of victim pair's lcm

### EVT_SWEEP_RESULT (type 23) — `kthread.cc:1995` and V1–V5 variants (types 23–27)

**For lightweight case (best >= 0):**
Emitted when cooperative sweep completes and either finds an accepted
reducer (lightweight emit, no aux) or finds no reducer (full emit with aux).

- `arg_a` = `s` — slot index
- `arg_b` = `best` — selected reducer T-index (or u32(-1) if no reducer found)
- `arg_c` = `n_examined` — number of T entries examined by sweep
- `aux_off` = 0 (lightweight case, no attribution payload)
- `poly_ptr_1` = `p_Copy(ap->P.p)` — captured snapshot of active poly at sweep completion
- `poly_ptr_2` = `best >= 0 ? p_Copy(strat->T[best].p) : NULL` — captured snapshot of selected reducer (if any)

**For survivor case (best < 0) and versioned variants (V1–V5):**
See "Events that carry aux payloads" section (lines 145–275) for detailed
layout of the 40–48 byte aux rows and the prelude.

- `arg_a` = `s` — slot index
- `arg_b` = `best` — merged best (u32(-1) for survivor, else best_good or best_reducer)
- `arg_c` = `n_examined` — number of T entries examined
- `aux_off` ≠ 0 (when best < 0) — offset to aux payload with per-T attribution rows
- `poly_ptr_1` = `p_Copy(ap->P.p)` — captured snapshot of active poly
- `poly_ptr_2` = per-version (NULL in V1, T[best].p in V2+)

**Versions:**
- **V1** (type 23): 16-byte aux rows `{T_idx, reason, entry_ptr}` (task 329)
- **V2** (type 24): 32-byte rows add `sweep_ptr, sev_sweep` (task 332)
- **V3** (type 25): 40-byte rows add `sev_computed` (task 334)
- **V4** (type 26): 32-byte prelude `{n_examined, pad, apP_sev_stored, apP_sev_computed, apP_not_sev_stored}` + V3 rows (task 335)
- **V5** (type 27): 64-byte prelude adds `{apP_p_ptr, apP_exps_packed, apP_currRing_ptr}` + 48-byte rows add `T_exps_packed` (task parallel-bba-raw-exp-probe)

For details, see kevlog.h:103–170 and the "Events that carry aux payloads" section.

## Checker consumers

All three live at `~/project/tools/`:

- `check-event-log.py` — the analyzer. Invariants 1-6 verify
  structural correctness (ENTERPAIRS completeness, L-drain,
  survivor bookkeeping, kill-reason validity, …). Invariant 7
  scans SWEEP_RESULT aux for SEV-filter false-rejects and does
  the H1/H2 / H1a/H1b / H1a.1-3 classification. Invariant 8
  runs Buchberger closure via a Singular subprocess on the
  `-full-polys.txt` generators plus a redundancy analysis of
  `std(S)`.
- `inspect-event-log.py` — the ad-hoc decoder. `--verbose`,
  `--filter-type`, `--dump N` flags. Handles all four SWEEP_RESULT
  versions.
- `sweep.sh` — the 200/400-iteration reproducer harness.

## Design decisions worth knowing

**Why 1 GiB, not a ring buffer**: the checker needs a complete log.
Dropping events to wrap would silently break invariants. Machine
has 29 GiB free; 1 GiB is comfortable. (Task 325 considered a ring
buffer and the user explicitly vetoed it.)

**Why only dump on failure**: dumping adds I/O, and the typical
pattern is `sweep.sh` running hundreds of iterations where almost
all pass. Lazy dump keeps the good-run overhead to the emit path
(a few ns per event).

**Why versioned SWEEP_RESULT events rather than in-place schema
changes**: each investigation task adds a probe. Old dumps from
previous tasks are still useful evidence for the checker to
classify. Versioning by type tag lets the checker handle
mixed-version archives and keeps the decoder dispatch clean.

**Why the `merged_best == -1` emit gate**: without it, every sweep
writes per-T aux for every slot — enough overhead to suppress the
bug (Heisenberg). The survivor case is the only one the checker
cares about for sev-filter analysis, so gating on it preserves
signal at a tiny fraction of the cost.
