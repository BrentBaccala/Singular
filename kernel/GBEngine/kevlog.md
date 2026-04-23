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

- Mutex-guarded `unordered_map<source_ptr, {copy, lmRing, tailRing}>`.
- First sight of a source → `p_Copy(src, ring)` → store the fresh
  copy in the map.
- Later sights → return the cached copy (preserves pointer
  identity across events: "same source poly" compares equal in
  the log).
- At `kevlog_shutdown()`, every copy gets `p_Delete`d.

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
   `copy_addr<TAB>lm<TAB>source_addr`. The copy column is what
   `*.bin` references; the source column (added in task 332) is
   the raw pre-copy pointer, used to match `SWEEP_RESULT_V2+`
   `sweep_ptr` entries against `entry_ptr` in the checker.
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
