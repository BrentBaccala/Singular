# Singular-rustgb — per-worktree notes

This is the `rustgb-integration` worktree of `~/Singular`. It hosts the
`singrust` dyn_module (links against `librustgb.so` from `~/rustgb`) and
the `rustgb-dispatch.lib` Singular shim. General rustgb docs live in
`~/CLAUDE.md` under "Rust Port (rustgb)".

## Staging validation under the subagent stream watchdog

The staging validation suite (`~/project/run-rustgb-staging-validation.sh`)
runs three large redsb tests from `~/test-cases/singular/`. Per-test
walls:

| test              | Vec     | List (pre-splice) | List (post-splice) |
|-------------------|--------:|------------------:|-------------------:|
| staging-5101449   |  ~190 s |    ~730 s         |   ~400–500 s (expected) |
| staging-5104053   |  ~4 min |   ~15 min         |   ~8–10 min (projected) |
| staging-5106746   |  ~6 min |   ~23 min         |   ~12–15 min (projected) |

Any single run of the full script exceeds the 600 s stream-watchdog cap
that `~/CLAUDE.md` describes. In subagents, a blocking MCP bash call
with `timeout=3600` returns no stream events while it waits, so the
watchdog kills the agent at 600 s — even though the command itself
would have finished correctly.

### The pattern: fire-and-poll

Use `~/project/run-rustgb-one-staging.sh <tag>` instead of the full
script. It:

1. Wraps one staging `.tst` with the `tst_init()` + `LIB` preamble.
2. Fires the Singular run in the background (nohup + setsid, output
   tee'd to a per-run `LOG` file).
3. Returns **immediately** to the caller with `PID=...`, `WORKDIR=...`,
   `LOG=...`, `RESULT=...`.

Typical flow:

```bash
# Kick off (returns in < 1 s):
~/project/run-rustgb-one-staging.sh staging-5101449 --build-flag vec
# → PID=12345, WORKDIR=..., LOG=..., RESULT=...

# Poll from later turns (quick tool call, no blocking):
ps -p <PID> > /dev/null && echo running || echo done

# When done, read the summary:
cat <RESULT>
```

Between polls, do other work (edit docs, run unit tests on a different
backend, etc.) so the stream keeps getting events. The watchdog only
cares about gaps in the JSONL stream, not about what's in them.

The `--build-flag` arg is a label used in the workdir path (not a
build-system flag). Pass `vec` / `list` / `list-splice` etc. to keep
concurrent and back-to-back runs' workdirs distinct.

### Do not rebuild `librustgb.so` while a run is in flight

The staging run `dlopen`s `librustgb.so` once at Singular startup, so a
rebuild partway through does not affect an already-running test. But
starting a new test after a rebuild picks up the new library — which
is exactly what you want when comparing Vec vs List. Sequence:

```
cargo build --release                             # Vec
./run-rustgb-one-staging.sh 5101449 --build-flag vec
# ... wait for it to finish ...
cargo build --release --features linked_list_poly # List
./run-rustgb-one-staging.sh 5101449 --build-flag list
# ...etc
```

### Commit between runs

Each poll-and-wait cycle is a natural commit point. If the subagent is
killed by the watchdog during a long idle spell, the work done before
the kill is already in git; only the in-flight staging run is lost
(and that's recoverable by re-running a single tag).

### When the full script is fine

From the top-level interactive session (not inside a subagent) the
watchdog does not apply, and `~/project/run-rustgb-staging-validation.sh`
works end-to-end. Only split across per-tag runs when you're running
under the task runner.

## Build

See the Rust Port section of `~/CLAUDE.md`. Short version:

```bash
cd ~/rustgb && ~/.cargo/bin/cargo build --release
cd ~/Singular-rustgb && make -j"$(nproc)" && make install
```

## Tests

`cargo test --release` (Vec) and
`cargo test --release --features linked_list_poly` (List) both in
`~/rustgb`; the Singular-side integration is exercised only through
`run-rustgb-staging-validation.sh` and the gb_invariant suite.
