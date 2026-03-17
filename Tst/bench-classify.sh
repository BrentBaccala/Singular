#!/bin/bash
# bench-classify.sh - Classify a test and determine iteration count
#
# Tries Class C (single run) to get baseline timing, then tests Class A
# with a small iteration count (3) to see if looping works.
# Reports the class and calculated iteration count for ~TARGET_TIME seconds.
#
# Usage:
#   bench-classify.sh [--singular PATH] [--target-time 20] TESTFILE
#
# Output (tab-separated):
#   TESTFILE  CLASS  ITERATIONS  SINGLE_RUN_MS

set -e

SINGULAR="${SINGULAR:-Singular}"
TARGET_TIME=20  # seconds
MAX_ITERATIONS=1000
TEST_ITERATIONS=3  # small count to test if class works
TIMEOUT=120  # seconds per attempt
TESTFILE=""
WRAPPER="$(dirname "$0")/bench-wrapper.sh"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --singular) SINGULAR="$2"; shift 2 ;;
        --target-time) TARGET_TIME="$2"; shift 2 ;;
        --timeout) TIMEOUT="$2"; shift 2 ;;
        *) TESTFILE="$1"; shift ;;
    esac
done

if [[ -z "$TESTFILE" ]]; then
    echo "Usage: bench-classify.sh [--singular PATH] TESTFILE" >&2
    exit 1
fi

# Try a class: returns 0 if it works, 1 if not
# Sets BENCH_TIME_MS on success
try_class() {
    local class="$1"
    local iters="$2"
    local output

    output=$(timeout "$TIMEOUT" bash -c \
        "'$WRAPPER' --class '$class' --iterations '$iters' '$TESTFILE' | '$SINGULAR' -q 2>&1") || return 1

    # Check for BENCH_TIME in output
    local bench_line
    bench_line=$(echo "$output" | grep '^BENCH_TIME:' | tail -1) || return 1

    if [[ -z "$bench_line" ]]; then
        return 1
    fi

    # Extract the time value
    BENCH_TIME_MS=$(echo "$bench_line" | sed 's/BENCH_TIME:[[:space:]]*//' | tr -d '[:space:]')

    # Check it's a valid number (may be negative if timer wraps)
    if ! [[ "$BENCH_TIME_MS" =~ ^-?[0-9]+$ ]]; then
        return 1
    fi

    return 0
}

# Phase 1: Get baseline timing with Class C (single run)
if ! try_class C 1; then
    # Test doesn't even run - skip it
    echo -e "${TESTFILE}\tFAIL\t0\t0"
    exit 1
fi

SINGLE_TIME_MS=$BENCH_TIME_MS

# If single run takes > 10s, use Class C with 1 iteration
if [[ $SINGLE_TIME_MS -ge 10000 ]]; then
    echo -e "${TESTFILE}\tC\t1\t${SINGLE_TIME_MS}"
    exit 0
fi

# Calculate target iterations based on single-run timing
if [[ $SINGLE_TIME_MS -le 0 ]]; then
    TARGET_N=$MAX_ITERATIONS
else
    TARGET_N=$(( (TARGET_TIME * 1000 + SINGLE_TIME_MS - 1) / SINGLE_TIME_MS ))
fi

# Cap iterations
if [[ $TARGET_N -gt $MAX_ITERATIONS ]]; then
    TARGET_N=$MAX_ITERATIONS
fi
if [[ $TARGET_N -lt 1 ]]; then
    TARGET_N=1
fi

# If target is 1, no need for looping - use Class C
if [[ $TARGET_N -le 1 ]]; then
    echo -e "${TESTFILE}\tC\t1\t${SINGLE_TIME_MS}"
    exit 0
fi

# Phase 2: Test if Class A works with a small iteration count
if try_class A "$TEST_ITERATIONS"; then
    echo -e "${TESTFILE}\tA\t${TARGET_N}\t${SINGLE_TIME_MS}"
    exit 0
fi

# Phase 3: Try Class B
if try_class B "$TEST_ITERATIONS"; then
    echo -e "${TESTFILE}\tB\t${TARGET_N}\t${SINGLE_TIME_MS}"
    exit 0
fi

# Phase 4: Try Class D
if try_class D "$TEST_ITERATIONS"; then
    echo -e "${TESTFILE}\tD\t${TARGET_N}\t${SINGLE_TIME_MS}"
    exit 0
fi

# Fallback: Class C single run
echo -e "${TESTFILE}\tC\t1\t${SINGLE_TIME_MS}"
exit 0
