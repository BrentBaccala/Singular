# Singular benchmarking framework — `~/Singular-benchmark/Tst/`

Tooling that runs Singular `.tst` files as performance benchmarks. The
core problem: many tests complete in <1 s, so single-run timing is
dominated by Singular's startup/parser overhead instead of the
computation. This framework wraps tests in loops (or chooses a
no-loop variant when looping doesn't work) so each measurement
runs for ~20 s of computation, then uses Singular's internal
`timer` / `rtimer` to isolate the work from process-launch cost.

Designed and built by **task 203 (`benchmark-test-suite`)** —
`~/project/prompts/benchmark-test-suite` reads like the spec. Reach
for that prompt before changing the framework's shape; this CLAUDE.md
is the day-to-day reference.

## Workflow

1. **Classify** each `.tst` once (writes a per-suite
   `<lstname>_classification.csv` next to the `.lst` file).
2. **Suite-run** with one or more builds (builds the wrapped input
   once per test using the saved class + iteration count, runs each
   build × test × N runs through a single Singular invocation, emits
   one CSV row per run).
3. **Analyze / plot** with the Python tools.

Classification is build-independent (the test's loop-tolerance
shape doesn't change with build), so the same classification CSV
serves every build, and the calibrated iteration count keeps each
benchmark cell at roughly the target wall time.

## Scripts

| Script | Role |
|---|---|
| `bench-wrapper.sh` | takes `--class X --iterations N TESTFILE`, emits wrapped Singular input on stdout. The unit primitive — pipe to `Singular -q`. |
| `bench-classify.sh` | tries Class C first to get a single-run baseline, then Class A with a tiny iteration count to see if looping works. Falls back to B / D as needed. Emits one CSV row: `testname.tst,CLASS,ITERATIONS,SINGLE_RUN_US`. The 4th column is microseconds (matching the wrapper's `BENCH_WALL:` units); not consumed by `bench-suite.sh`. |
| `bench-suite.sh` | the orchestrator. Takes builds + `.lst` files (or individual `.tst` files), looks up each test's classification, builds the wrapped input, runs under `setarch -R` (ASLR off) + optional `numactl` / `taskset`, emits a wide CSV. |
| `bench-launch.sh` | nohup-friendly background-run skeleton. scp to the host, `bash ~/bench-launch.sh`. |
| `bench-analyze.py` | per-test mean / stddev / CV from a benchmark CSV. `--during START END` flags rows during a busy time window. |
| `bench-stripplot.py` | strip plot of timings, sorted shortest-first, dot-per-run colour-by-build. |
| `bench-estimate.py` | wall-clock estimate for an upcoming run at a given `--scale`. |
| `bench-check-timing.py` | sanity-checks measured times against classification estimates. |

## Wrapper classes (the body shape `bench-wrapper.sh` emits)

| Class | Shape | Use when |
|---|---|---|
| **A** | one `for` loop running the body N times | most short tests; default first attempt |
| **B** | body wrapped in a `proc`, called N times | Class A redefines a name on iter 2 (rings, idents) — the proc gives a local scope that gets torn down each call |
| **C** | single run, no warmup | test already takes >10 s, or Class A / B would change semantics |
| **D** | textually concatenate the body N times in the wrapped file | neither loop nor proc works (e.g. test does its own `kill` / state mutation) — heavier on parse cost but still one Singular launch |
| **E** | same body shape as C, but Singular's internal `timer` reads 0 | tests that kill the basering or otherwise destroy `timer`; only `ext_wall_ns` (`/usr/bin/time`-measured) is meaningful |
| **X** | excluded — emit no wrapper, suite skips the test | tests that reset `--ticks-per-sec`, dump files to disk, or otherwise can't be benchmarked at all |

The wrapped input emits these stdout markers, parsed by
`bench-suite.sh` to fill the CSV columns:

```
BENCH_WARMUP_CPU: <us>
BENCH_WARMUP_WALL: <us>
BENCH_CPU: <us>
BENCH_WALL: <us>
```

Classes A and B include a single warmup iteration outside the timed
loop so library load + cache warm don't pollute the measurement.
Classes C, D, E omit the warmup.

## Classification CSV format

Filename: `<listname>_classification.csv` next to the corresponding
`.lst` (e.g. `Short/ok_s.lst` ↔ `Short/ok_s_classification.csv`).
Keys are filenames (basename), not paths — so the same `.tst`
referenced from multiple lists shares a single classification row.

Columns: `testname.tst,CLASS,ITERATIONS,TIMEOUT_SECS`.

Example rows:
```
abusalem.tst,A,8,60
arnold.tst,C,1,300
arr.tst,X,0,0
bug_732.tst,A,105287,60
```

The `TIMEOUT_SECS` column was added later (commit `ab279063a`); the
suite defaults are 60 s for looped classes and 300 s for Class C.
Bump per test if a slower build needs more headroom.

## Suites classified today

| Classification CSV | Tests covered | Source list |
|---|---|---|
| `Tst/Short/ok_s_classification.csv` | classic short suite | `Tst/Short/ok_s.lst` |
| `Tst/Long/ok_l_classification.csv` | long suite | `Tst/Long/ok_l.lst` |
| `Tst/Buch/buch_classification.csv` | Buchberger-focused tests | `Tst/Buch/buch.lst` |
| `Tst/Old/universal_classification.csv` | legacy universal tests | `Tst/Old/universal.lst` |
| `~/test-cases/singular/staging_classification.csv` | helium staging tests | `~/test-cases/singular/staging.lst` |

**Not classified** as of 2026-04-25: the `gb_invariant.lst` and
`gb_portable/gb_portable2.lst` lists in `~/Singular-gb-testsuite`.
These are the lists the rustgb invariant-validation harness drives
through (`~/project/run-rustgb-invariant-validation.sh`), and they
need classification before any of those tests can be benchmarked
under bench-suite.sh — `libehv_s` in particular runs in ~0.5 s so
single-run timing is mostly Singular startup.

## Adding a new suite

```
bench-suite.sh --save-classify <build-spec> -- <path/to/new.lst>
```

Without a saved classification, `bench-suite.sh` runs
`bench-classify.sh` on each test on the fly — works once but takes
~5–30 s per test for the calibration probes. With `--save-classify`
the result is written next to the `.lst` and reused on every later
run. After this completes, commit the new
`<listname>_classification.csv` so future runs (and other branches)
inherit it.

To classify a single test manually:

```
bench-classify.sh --singular <path/to/Singular> --target-time 20 \
  /path/to/test.tst > >(tee -a path/to/<list>_classification.csv)
```

## Loading dispatch shims / extra LIBs (`--preamble-file`)

For builds that need a Singular-level preamble before each test (e.g.
the rustgb dispatch shim's `LIB "singrust.so"; LIB "rustgb-dispatch.lib";`),
pass `--preamble-file FILE` to `bench-suite.sh`. The file's contents
are prepended to the wrapped body **before the timer-start line**,
so the LIB-load cost is paid once per Singular invocation but isn't
charged against the per-iteration timing.

A single preamble file works for all rustgb tests — the dispatch
shim's filter (Z/p, degrevlex, ≤31 vars) is runtime, so the same
preamble routes std() correctly across every workload.

```
# ~/rustgb-preamble.sing
LIB "singrust.so";
LIB "rustgb-dispatch.lib";
```

```bash
bench-suite.sh -n 3 --taskset "-c 11" --save-classify \
  --preamble-file ~/rustgb-preamble.sing \
  "rustgb=LD_LIBRARY_PATH=$HOME/rustgb/target/release \
         SINGULARPATH=$HOME/Singular-rustgb/install/lib/singular/MOD \
         $HOME/Singular-rustgb/install/bin/Singular" \
  -- ~/Singular-gb-testsuite/Tst/Short/libehv_s.tst
```

**Caveat:** classification (`bench-classify.sh`) does not yet honour
the preamble. Iteration counts get calibrated against the
no-preamble Singular, so any per-call overhead the preamble adds
(or savings it produces, in the rustgb-dispatch case) shifts the
post-preamble measurement away from the 20 s target. Acceptable
when the shift is small; if a future preamble adds significant
per-call cost, extend bench-classify.sh similarly.

## Standard invocation patterns

c200-1 isolated benchmark host (CPU 11, NUMA-local memory on socket 1):

```
bench-suite.sh -n 5 \
  --taskset "-c 11" --numactl "--membind=1" \
  -o results.csv \
  spielwiese=$HOME/Singular/install/bin/Singular \
  LSet=$HOME/Singular-LSet/install/bin/Singular \
  -- Short/ok_s.lst Buch/buch.lst Long/ok_l.lst
```

Build specs accept arbitrary command strings, so env-var-prefixed
invocations work directly:

```
"mimalloc=LD_PRELOAD=/usr/lib/libmimalloc.so /path/to/Singular"
```

`bench-launch.sh` is the canonical "kick off a long run via SSH"
template — edit the build list + lst file list, scp it to the
host, run with `ssh host 'bash ~/bench-launch.sh'`.

## Output CSV schema

```
build,test_name,class,iterations,run,
warmup_cpu_us,warmup_wall_us,cpu_us,wall_us,ext_wall_ns,
clean_exit,list_file,
instructions,cycles,cache_misses,branch_misses
```

`cpu_us` / `wall_us` come from Singular's internal `timer` / `rtimer`
(in microseconds, with `system("--ticks-per-sec", 1000000)` set in
the wrapper preamble). `ext_wall_ns` is `/usr/bin/time` measured
externally — useful when Class E's internal timer is unreliable.
The `instructions` / `cycles` / `cache_misses` / `branch_misses`
columns come from `perf stat` (commit `bc13bd031`); zero on hosts
where perf isn't usable.

## Common gotchas

- **Default classification is build-independent**, so calibration
  is done with the first build in the spec list. If that build is
  pathologically slow on a test, the iteration count gets
  under-calibrated for the others. The fix is `--scale`: e.g.
  `--scale 0.25` quarters every iteration count for a quick
  smoke-run, `--scale 2` doubles for a longer / lower-noise run.
- **ASLR is force-disabled** via `setarch -R` (commit `0a082378b`)
  for deterministic perf-stat instruction counts. Don't add
  `setarch` redundantly.
- **`bench-suite.sh` accepts individual `.tst` files** as well as
  `.lst` files (commit `dd5b62dca`) — useful for one-off cells.
- **Class A redefinition errors silently invalidate timing.**
  `bench-classify.sh` tries to detect this, but a test that
  succeeds-with-warnings on iter 2 can pass calibration and still
  be wrong. Sanity-check `cpu_us / iterations` against
  `bench-classify.sh`'s `SINGLE_RUN_US` for any new test.
- **`--no-warn`** is **not** set by the framework — Singular
  warnings (e.g. "redefining ring") go to stderr and are visible
  in the suite's log. Useful for catching silently-broken loop
  wraps.

## References

- Spec: `~/project/prompts/benchmark-test-suite` (task 203 prompt).
- First-run report: `~/project/docs/c200-1-benchmark-report.md`
  (task 224).
- Calibration data: `~/project/docs/benchmark-calibration-data.md`.
- Reviews / results: `benchmark-review.md`, `benchmark-results.md`,
  `sev-prefilter-benchmark-report.md` in `~/project/docs/`.
- Feature evolution: `git -C ~/Singular-benchmark log --
  Tst/bench-*.sh Tst/*_classification.csv`.
