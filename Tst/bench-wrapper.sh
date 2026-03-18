#!/bin/bash
# bench-wrapper.sh - Wraps a Singular test file for benchmarking
#
# Produces wrapped Singular input on stdout. Pipe to Singular -q.
#
# Usage:
#   bench-wrapper.sh --class A --iterations 10 Short/facstd.tst | Singular -q
#   bench-wrapper.sh --class C Short/facstd.tst | Singular -q

set -e

CLASS="A"
ITERATIONS=10
TESTFILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --class) CLASS="$2"; shift 2 ;;
        --iterations|-n) ITERATIONS="$2"; shift 2 ;;
        *) TESTFILE="$1"; shift ;;
    esac
done

if [[ -z "$TESTFILE" ]]; then
    echo "Usage: bench-wrapper.sh [--class A|B|C|D] [--iterations N] TESTFILE" >&2
    exit 1
fi

if [[ ! -f "$TESTFILE" ]]; then
    echo "Error: test file '$TESTFILE' not found" >&2
    exit 1
fi

# Strip the test body: remove tst_init/tst_status/tst.lib lines,
# remove trailing 'exit;', and get the core computation.
# Keep LIB loads except tst.lib.
strip_test_body() {
    sed -E \
        -e '/^[[:space:]]*LIB[[:space:]]+"tst\.lib"/d' \
        -e '/^[[:space:]]*tst_init\b/d' \
        -e '/^[[:space:]]*tst_status\b/d' \
        -e '/^[[:space:]]*exit[[:space:]]*;/d' \
        -e '/^[[:space:]]*\$[[:space:]]*$/d' \
        -e '/^[[:space:]]*timer[[:space:]]*=[[:space:]]*[0-9]/d' \
        -e '/^[[:space:]]*int[[:space:]]+t[[:space:]]*=[[:space:]]*timer/d' \
        -e '/^[[:space:]]*int[[:space:]]+elapsed[[:space:]]*=[[:space:]]*timer/d' \
        -e '/^[[:space:]]*printf\(.*elapsed/d' \
        -e '/^[[:space:]]*printf\(.*Basis has/d' \
        "$1"
}

case "$CLASS" in
    A)
        # Class A: Loop wrapper - run test body N times in a loop
        # Uses "execute" approach: load body into string, execute in loop
        # This avoids redefinition issues for some cases
        cat <<'HEADER'
system("--ticks-per-sec", 1000);
HEADER
        # Emit the test body once outside the loop for definitions,
        # then time repeated executions. Actually, for Class A we just
        # wrap in a simple loop. Variable redefinitions will be caught
        # by the classifier.
        echo "int benchX_N = $ITERATIONS;"
        # Warmup: one untimed iteration for library loading and cache priming
        strip_test_body "$TESTFILE"
        echo 'int benchX_start = timer;'
        echo 'for (int benchX_i = 1; benchX_i <= benchX_N; benchX_i++) {'
        strip_test_body "$TESTFILE"
        echo '}'
        echo 'int benchX_end = timer;'
        echo '"BENCH_TIME:", benchX_end - benchX_start;'
        echo 'exit;'
        ;;
    B)
        # Class B: Proc wrapper - wrap body in a proc for local scope
        cat <<'HEADER'
system("--ticks-per-sec", 1000);
HEADER
        # Extract LIB lines and put them outside the proc
        grep -E '^[[:space:]]*LIB[[:space:]]+' "$TESTFILE" | grep -v 'tst\.lib' || true
        echo 'proc benchX_body() {'
        # Strip LIB lines from the body since they're outside
        strip_test_body "$TESTFILE" | sed -E '/^[[:space:]]*LIB[[:space:]]+/d'
        echo '  return();'
        echo '}'
        echo "int benchX_N = $ITERATIONS;"
        # Warmup: one untimed iteration for library loading and cache priming
        echo '  benchX_body();'
        echo 'int benchX_start = timer;'
        echo 'for (int benchX_i = 1; benchX_i <= benchX_N; benchX_i++) {'
        echo '  benchX_body();'
        echo '}'
        echo 'int benchX_end = timer;'
        echo '"BENCH_TIME:", benchX_end - benchX_start;'
        echo 'exit;'
        ;;
    C)
        # Class C: Single run - just time one execution (no warmup)
        cat <<'HEADER'
system("--ticks-per-sec", 1000);
int benchX_start = timer;
HEADER
        strip_test_body "$TESTFILE"
        echo 'int benchX_end = timer;'
        echo '"BENCH_TIME:", benchX_end - benchX_start;'
        echo 'exit;'
        ;;
    D)
        # Class D: Concatenate file N times
        cat <<'HEADER'
system("--ticks-per-sec", 1000);
HEADER
        # Warmup: one untimed iteration
        strip_test_body "$TESTFILE"
        echo 'int benchX_start = timer;'
        for ((i = 1; i <= ITERATIONS; i++)); do
            strip_test_body "$TESTFILE"
        done
        echo 'int benchX_end = timer;'
        echo '"BENCH_TIME:", benchX_end - benchX_start;'
        echo 'exit;'
        ;;
    *)
        echo "Error: unknown class '$CLASS'" >&2
        exit 1
        ;;
esac
