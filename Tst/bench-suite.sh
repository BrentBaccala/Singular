#!/bin/bash
# bench-suite.sh - Benchmark all tests from a test list file
#
# Usage:
#   bench-suite.sh [OPTIONS] LISTFILE [LISTFILE2 ...]
#
# Options:
#   --singular PATH     Path to Singular binary
#   --target-time SECS  Target time per test (default: 20)
#   --timeout SECS      Timeout per individual run (default: 300)
#   --output FILE       CSV output file (default: stdout)
#   --classify-file F   Save/load classification from file
#   --numactl ARGS      numactl arguments (e.g., "--membind=1")
#   --taskset ARGS      taskset arguments (e.g., "-c 11")
#   --skip-classify     Skip classification, load from --classify-file
#   --skip-to TEST      Skip tests until TEST is reached (resume)

set -e

SINGULAR="${SINGULAR:-Singular}"
TARGET_TIME=20
TIMEOUT=300
OUTPUT=""
CLASSIFY_FILE=""
NUMACTL_ARGS=""
TASKSET_ARGS=""
SKIP_CLASSIFY=0
SKIP_TO=""
APPEND=0
LISTFILES=()

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WRAPPER="$SCRIPT_DIR/bench-wrapper.sh"
CLASSIFIER="$SCRIPT_DIR/bench-classify.sh"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --singular) SINGULAR="$2"; shift 2 ;;
        --target-time) TARGET_TIME="$2"; shift 2 ;;
        --timeout) TIMEOUT="$2"; shift 2 ;;
        --output|-o) OUTPUT="$2"; shift 2 ;;
        --classify-file) CLASSIFY_FILE="$2"; shift 2 ;;
        --numactl) NUMACTL_ARGS="$2"; shift 2 ;;
        --taskset) TASKSET_ARGS="$2"; shift 2 ;;
        --skip-classify) SKIP_CLASSIFY=1; shift ;;
        --skip-to) SKIP_TO="$2"; shift 2 ;;
        --append) APPEND=1; shift ;;
        *) LISTFILES+=("$1"); shift ;;
    esac
done

if [[ ${#LISTFILES[@]} -eq 0 ]]; then
    echo "Usage: bench-suite.sh [OPTIONS] LISTFILE [LISTFILE2 ...]" >&2
    exit 1
fi

# Build the Singular command with optional numactl/taskset
build_singular_cmd() {
    local cmd=""
    if [[ -n "$NUMACTL_ARGS" ]]; then
        cmd="numactl $NUMACTL_ARGS "
    fi
    if [[ -n "$TASKSET_ARGS" ]]; then
        cmd="${cmd}taskset $TASKSET_ARGS "
    fi
    cmd="${cmd}${SINGULAR} -q"
    echo "$cmd"
}

SINGULAR_CMD=$(build_singular_cmd)

# CSV header
csv_header() {
    echo "test_name,class,iterations,total_time_ms,per_iter_ms,clean_exit,list_file"
}

# Collect tests from all list files
collect_tests() {
    for listfile in "${LISTFILES[@]}"; do
        local dir
        dir=$(dirname "$listfile")
        while IFS= read -r line; do
            # Skip empty lines, comments, and disabled tests (;prefix)
            [[ -z "$line" ]] && continue
            [[ "$line" =~ ^# ]] && continue
            [[ "$line" =~ ^';' ]] && continue
            # Strip trailing whitespace
            line=$(echo "$line" | sed 's/[[:space:]]*$//')
            local testpath="${dir}/${line}.tst"
            if [[ -f "$testpath" ]]; then
                echo "$testpath|$listfile"
            else
                echo "Warning: test file not found: $testpath" >&2
            fi
        done < "$listfile"
    done
}

# Classification phase
declare -A TEST_CLASS
declare -A TEST_ITERS

classify_test() {
    local testfile="$1"
    local result

    # Export SINGULAR for the classifier
    export SINGULAR

    result=$(bash "$CLASSIFIER" --singular "$SINGULAR" --target-time "$TARGET_TIME" --timeout "$TIMEOUT" "$testfile" 2>/dev/null) || {
        TEST_CLASS["$testfile"]="FAIL"
        TEST_ITERS["$testfile"]=0
        return 1
    }

    local class iters
    class=$(echo "$result" | cut -f2)
    iters=$(echo "$result" | cut -f3)

    TEST_CLASS["$testfile"]="$class"
    TEST_ITERS["$testfile"]="$iters"
    return 0
}

# Benchmark a single test
benchmark_test() {
    local testfile="$1"
    local class="${TEST_CLASS[$testfile]}"
    local iters="${TEST_ITERS[$testfile]}"
    local listfile="$2"

    if [[ "$class" == "FAIL" ]]; then
        echo "${testfile},FAIL,0,0,0,no,${listfile}"
        return
    fi

    # Run the actual benchmark
    local output exit_code=0
    output=$(timeout "$TIMEOUT" bash -c \
        "'$WRAPPER' --class '$class' --iterations '$iters' '$testfile' | $SINGULAR_CMD 2>&1") || exit_code=$?

    local bench_time=0
    local clean="yes"

    if [[ $exit_code -ne 0 ]]; then
        clean="no"
    fi

    # Extract BENCH_TIME
    local bench_line
    bench_line=$(echo "$output" | grep '^BENCH_TIME:' | tail -1) || true

    if [[ -n "$bench_line" ]]; then
        bench_time=$(echo "$bench_line" | sed 's/BENCH_TIME:[[:space:]]*//' | tr -d '[:space:]')
        if ! [[ "$bench_time" =~ ^-?[0-9]+$ ]]; then
            bench_time=0
            clean="no"
        fi
    else
        clean="no"
    fi

    # Check for errors in output (but not warnings)
    if echo "$output" | grep -qi '^ *\? error\|^   \? Segment\|SIGSEGV\|Abort'; then
        clean="no"
    fi

    # Calculate per-iteration time
    local per_iter=0
    if [[ $iters -gt 0 && $bench_time -gt 0 ]]; then
        per_iter=$(( bench_time / iters ))
    fi

    echo "${testfile},${class},${iters},${bench_time},${per_iter},${clean},${listfile}"
}

# Load classification from file
load_classification() {
    local file="$1"
    while IFS=$'\t' read -r testfile class iters _; do
        TEST_CLASS["$testfile"]="$class"
        TEST_ITERS["$testfile"]="$iters"
    done < "$file"
}

# Save classification to file
save_classification() {
    local file="$1"
    for testfile in "${!TEST_CLASS[@]}"; do
        echo -e "${testfile}\t${TEST_CLASS[$testfile]}\t${TEST_ITERS[$testfile]}"
    done | sort > "$file"
}

# Main execution
main() {
    local tests
    tests=$(collect_tests)
    local total
    total=$(echo "$tests" | wc -l)

    echo "Benchmarking $total tests with Singular: $SINGULAR" >&2
    echo "Target time: ${TARGET_TIME}s, Timeout: ${TIMEOUT}s" >&2
    echo "Command: $SINGULAR_CMD" >&2
    echo "" >&2

    # Load existing classification if requested
    if [[ $SKIP_CLASSIFY -eq 1 && -n "$CLASSIFY_FILE" && -f "$CLASSIFY_FILE" ]]; then
        echo "Loading classification from $CLASSIFY_FILE" >&2
        load_classification "$CLASSIFY_FILE"
    fi

    # Phase 1: Classification (if not skipping)
    if [[ $SKIP_CLASSIFY -eq 0 ]]; then
        echo "=== Phase 1: Classifying tests ===" >&2
        local count=0
        local skipping=1
        [[ -z "$SKIP_TO" ]] && skipping=0

        while IFS='|' read -r testfile listfile; do
            if [[ $skipping -eq 1 ]]; then
                if [[ "$testfile" == *"$SKIP_TO"* ]]; then
                    skipping=0
                else
                    continue
                fi
            fi

            count=$((count + 1))
            echo -ne "\r  [$count/$total] Classifying: $testfile ... " >&2

            if [[ -n "${TEST_CLASS[$testfile]+x}" ]]; then
                echo "cached (${TEST_CLASS[$testfile]})" >&2
                continue
            fi

            if classify_test "$testfile"; then
                echo "${TEST_CLASS[$testfile]} (N=${TEST_ITERS[$testfile]})" >&2
            else
                echo "FAIL" >&2
            fi
        done <<< "$tests"

        echo "" >&2
        echo "Classification complete." >&2

        # Save classification
        if [[ -n "$CLASSIFY_FILE" ]]; then
            save_classification "$CLASSIFY_FILE"
            echo "Saved classification to $CLASSIFY_FILE" >&2
        fi
    fi

    # Phase 2: Benchmarking
    echo "=== Phase 2: Benchmarking ===" >&2

    # Output CSV
    local csv_out
    if [[ -n "$OUTPUT" ]]; then
        if [[ $APPEND -eq 1 && -f "$OUTPUT" ]]; then
            csv_out="$OUTPUT"
        else
            csv_header > "$OUTPUT"
            csv_out="$OUTPUT"
        fi
    else
        csv_header
    fi

    local count=0
    local success=0
    local failed=0
    local skipping=1
    [[ -z "$SKIP_TO" ]] && skipping=0

    while IFS='|' read -r testfile listfile; do
        if [[ $skipping -eq 1 ]]; then
            if [[ "$testfile" == *"$SKIP_TO"* ]]; then
                skipping=0
            else
                continue
            fi
        fi

        count=$((count + 1))

        # Classify if not already done
        if [[ -z "${TEST_CLASS[$testfile]+x}" ]]; then
            classify_test "$testfile" || true
        fi

        local class="${TEST_CLASS[$testfile]}"
        echo -ne "\r  [$count/$total] Benchmarking ($class): $testfile ... " >&2

        local result
        result=$(benchmark_test "$testfile" "$listfile")

        if [[ -n "$csv_out" ]]; then
            echo "$result" >> "$csv_out"
        else
            echo "$result"
        fi

        local clean
        clean=$(echo "$result" | cut -d, -f6)
        if [[ "$clean" == "yes" ]]; then
            success=$((success + 1))
            echo "done" >&2
        else
            failed=$((failed + 1))
            echo "ISSUES" >&2
        fi
    done <<< "$tests"

    echo "" >&2
    echo "=== Summary ===" >&2
    echo "Total: $count, Success: $success, Issues: $failed" >&2

    # Class breakdown
    local class_a=0 class_b=0 class_c=0 class_d=0 class_fail=0
    for class in "${TEST_CLASS[@]}"; do
        case "$class" in
            A) class_a=$((class_a + 1)) ;;
            B) class_b=$((class_b + 1)) ;;
            C) class_c=$((class_c + 1)) ;;
            D) class_d=$((class_d + 1)) ;;
            FAIL) class_fail=$((class_fail + 1)) ;;
        esac
    done
    echo "Classes: A=$class_a, B=$class_b, C=$class_c, D=$class_d, FAIL=$class_fail" >&2

    if [[ -n "$csv_out" ]]; then
        echo "Results written to: $csv_out" >&2
    fi
}

main
