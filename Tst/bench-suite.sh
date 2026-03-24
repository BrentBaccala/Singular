#!/bin/bash
# bench-suite.sh - Benchmark Singular test suite across builds
#
# Usage:
#   bench-suite.sh [OPTIONS] name=command [name=command ...] -- LISTFILE [LISTFILE2 ...]
#
# Build specifications:
#   name=command where command is the full Singular invocation, e.g.:
#     spielwiese=/home/claude/Singular/install/bin/Singular
#     spielwiese-mimalloc="LD_PRELOAD=/usr/lib/libmimalloc.so /home/claude/Singular/install/bin/Singular"
#
# Classification:
#   For each LISTFILE, looks for a classification file in the same directory
#   named by replacing .lst with _classification.csv (e.g., ok_s.lst ->
#   ok_s_classification.csv). Auto-loads if found. Tests not in any
#   classification file are classified on the fly.
#
#   Use --save-classify to write/update classification files after classifying.
#   Without this flag, classification files are never modified.
#
# Options:
#   -n RUNS             Number of runs per test (default: 1)
#   --target-time SECS  Target time per test for calibration (default: 20)
#   --timeout SECS      Timeout per individual run (default: 300)
#   --output FILE       CSV output file (default: stdout)
#   --numactl ARGS      numactl arguments (e.g., "--membind=1")
#   --taskset ARGS      taskset arguments (e.g., "-c 11")
#   --scale FACTOR      Scale looped iteration counts (e.g., 0.5 halves, 2.0 doubles)
#   --save-classify     Write classification files after classifying new tests
#   --classify-from F   Load classifications from file F (read-only)
#   --skip-to TEST      Skip tests until TEST is reached (resume)
#   --append            Append to output file instead of overwriting

set -e

RUNS=1
TARGET_TIME=20
TIMEOUT=300
SCALE=""
OUTPUT=""
NUMACTL_ARGS=""
TASKSET_ARGS=""
SAVE_CLASSIFY=0
CLASSIFY_FROM=""
SKIP_TO=""
APPEND=0
LISTFILES=()

# Build specs: parallel arrays
declare -a BUILD_NAMES
declare -a BUILD_COMMANDS

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WRAPPER="$SCRIPT_DIR/bench-wrapper.sh"
CLASSIFIER="$SCRIPT_DIR/bench-classify.sh"

# Parse arguments
parsing_builds=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n) RUNS="$2"; shift 2 ;;
        --target-time) TARGET_TIME="$2"; shift 2 ;;
        --timeout) TIMEOUT="$2"; shift 2 ;;
        --output|-o) OUTPUT="$2"; shift 2 ;;
        --numactl) NUMACTL_ARGS="$2"; shift 2 ;;
        --taskset) TASKSET_ARGS="$2"; shift 2 ;;
        --scale) SCALE="$2"; shift 2 ;;
        --save-classify) SAVE_CLASSIFY=1; shift ;;
        --classify-from) CLASSIFY_FROM="$2"; shift 2 ;;
        --skip-to) SKIP_TO="$2"; shift 2 ;;
        --append) APPEND=1; shift ;;
        --) parsing_builds=0; shift ;;
        *=*)
            if [[ $parsing_builds -eq 1 ]]; then
                local_name="${1%%=*}"
                local_cmd="${1#*=}"
                BUILD_NAMES+=("$local_name")
                BUILD_COMMANDS+=("$local_cmd")
                shift
            else
                LISTFILES+=("$1"); shift
            fi
            ;;
        *)
            LISTFILES+=("$1"); shift
            ;;
    esac
done

# Backwards compatibility: if no builds specified, check for SINGULAR env var
if [[ ${#BUILD_NAMES[@]} -eq 0 ]]; then
    if [[ -n "${SINGULAR:-}" ]]; then
        BUILD_NAMES+=("default")
        BUILD_COMMANDS+=("$SINGULAR")
    else
        echo "Usage: bench-suite.sh [OPTIONS] name=command [...] -- LISTFILE [...]" >&2
        echo "" >&2
        echo "Examples:" >&2
        echo "  bench-suite.sh spielwiese=/path/to/Singular -- Short/ok_s.lst" >&2
        echo "  bench-suite.sh -n 5 spiel=/path/to/Singular LSet=/other/Singular -- Short/ok_s.lst" >&2
        echo '  bench-suite.sh "mimalloc=LD_PRELOAD=/usr/lib/libmimalloc.so /path/to/Singular" -- Short/ok_s.lst' >&2
        exit 1
    fi
fi

if [[ ${#LISTFILES[@]} -eq 0 ]]; then
    echo "Error: no list files specified. Use -- to separate builds from list files." >&2
    exit 1
fi

# Derive classification filename from list filename
# ok_s.lst -> ok_s_classification.csv
classify_file_for() {
    local lstfile="$1"
    local base
    base=$(basename "$lstfile" .lst)
    local dir
    dir=$(dirname "$lstfile")
    echo "${dir}/${base}_classification.csv"
}

# Build a command string with optional numactl/taskset prefix
build_full_cmd() {
    local base_cmd="$1"
    local cmd=""
    if [[ -n "$NUMACTL_ARGS" ]]; then
        cmd="numactl $NUMACTL_ARGS "
    fi
    if [[ -n "$TASKSET_ARGS" ]]; then
        cmd="${cmd}taskset $TASKSET_ARGS "
    fi
    cmd="${cmd}${base_cmd} -q"
    echo "$cmd"
}

# Use the first build for classification (classification is build-independent)
CLASSIFY_SINGULAR="${BUILD_COMMANDS[0]}"

# CSV header
csv_header() {
    echo "build,test_name,class,iterations,run,warmup_cpu_us,warmup_wall_us,cpu_us,wall_us,ext_wall_ns,clean_exit,list_file,instructions,cycles,cache_misses,branch_misses"
}

# Collect tests from list files and/or individual .tst files
collect_tests() {
    for listfile in "${LISTFILES[@]}"; do
        if [[ "$listfile" == *.tst ]]; then
            # Individual test file
            if [[ -f "$listfile" ]]; then
                echo "$listfile|$listfile"
            else
                echo "Warning: test file not found: $listfile" >&2
            fi
        else
            # List file
            local dir
            dir=$(dirname "$listfile")
            while IFS= read -r line; do
                [[ -z "$line" ]] && continue
                [[ "$line" =~ ^# ]] && continue
                [[ "$line" =~ ^';' ]] && continue
                line=$(echo "$line" | sed 's/[[:space:]]*$//')
                local testpath="${dir}/${line}.tst"
                if [[ -f "$testpath" ]]; then
                    echo "$testpath|$listfile"
                else
                    echo "Warning: test file not found: $testpath" >&2
                fi
            done < "$listfile"
        fi
    done
}

# Classification — keyed by filename (basename), not full path
declare -A TEST_CLASS
declare -A TEST_ITERS
declare -A TEST_TIMEOUT

# Track which list files had new classifications (for --save-classify)
declare -A CLASSIFY_DIRTY

# Get the classification key (basename) for a test file path
classify_key() {
    basename "$1"
}

classify_test() {
    local testfile="$1"
    local listfile="$2"
    local key
    key=$(classify_key "$testfile")
    local result

    result=$(bash "$CLASSIFIER" --singular "$CLASSIFY_SINGULAR" \
        --target-time "$TARGET_TIME" --timeout "$TIMEOUT" "$testfile" 2>/dev/null) || {
        TEST_CLASS["$key"]="FAIL"
        TEST_ITERS["$key"]=0
        TEST_TIMEOUT["$key"]=$TIMEOUT
        CLASSIFY_DIRTY["$listfile"]=1
        return 1
    }

    local class iters
    class=$(echo "$result" | cut -d, -f2)
    iters=$(echo "$result" | cut -d, -f3)

    TEST_CLASS["$key"]="$class"
    TEST_ITERS["$key"]="$iters"
    # Default timeout: 60 for looped classes, 300 for single run
    if [[ "$class" == "C" ]]; then
        TEST_TIMEOUT["$key"]=300
    else
        TEST_TIMEOUT["$key"]=60
    fi
    CLASSIFY_DIRTY["$listfile"]=1
    return 0
}

# Benchmark a single test with a specific build
benchmark_test() {
    local build_name="$1"
    local build_cmd="$2"
    local testfile="$3"
    local run_num="$4"
    local listfile="$5"
    local key
    key=$(classify_key "$testfile")
    local class="${TEST_CLASS[$key]}"
    local iters="${TEST_ITERS[$key]}"

    # Apply scaling to looped tests
    if [[ -n "$SCALE" && "$class" != "C" && "$class" != "E" && "$class" != "X" && "$class" != "FAIL" ]]; then
        iters=$(python3 -c "print(max(1, round($iters * $SCALE)))")
    fi

    if [[ "$class" == "FAIL" ]]; then
        echo "${build_name},${testfile},FAIL,0,${run_num},0,0,0,0,0,no,${listfile},0,0,0,0"
        return
    fi

    if [[ "$class" == "X" ]]; then
        return  # excluded from benchmarking, no output
    fi

    local full_cmd
    full_cmd=$(build_full_cmd "$build_cmd")

    # Run the actual benchmark
    local test_timeout="${TEST_TIMEOUT[$key]:-$TIMEOUT}"
    local output exit_code=0
    local start_ns end_ns
    local perf_output_file
    perf_output_file=$(mktemp /tmp/bench-perf.XXXXXX)

    start_ns=$(date +%s%N)
    output=$(timeout "$test_timeout" perf stat -e instructions,cycles,cache-misses,branch-misses \
        -o "$perf_output_file" \
        bash -c "'$WRAPPER' --class '$class' --iterations '$iters' '$testfile' | $full_cmd 2>&1") || exit_code=$?
    end_ns=$(date +%s%N)

    local warmup_cpu_us=0
    local warmup_wall_us=0
    local cpu_us=0
    local wall_us=0
    local ext_wall_ns=$(( end_ns - start_ns ))
    local clean="yes"
    local instructions=0
    local cycles=0
    local cache_misses=0
    local branch_misses=0

    if [[ $exit_code -ne 0 ]]; then
        clean="no"
    fi

    # Extract all four Singular timers (microseconds)
    extract_timer() {
        local label="$1"
        local line
        line=$(echo "$output" | grep "^${label}:" | tail -1) || true
        if [[ -n "$line" ]]; then
            local val
            val=$(echo "$line" | sed "s/${label}:[[:space:]]*//" | tr -d '[:space:]')
            if [[ "$val" =~ ^-?[0-9]+$ ]]; then
                echo "$val"
                return
            fi
        fi
        echo "0"
    }

    warmup_cpu_us=$(extract_timer "BENCH_WARMUP_CPU")
    warmup_wall_us=$(extract_timer "BENCH_WARMUP_WALL")
    cpu_us=$(extract_timer "BENCH_CPU")
    wall_us=$(extract_timer "BENCH_WALL")

    # Extract perf stat counters
    extract_perf_counter() {
        local counter="$1"
        local val
        val=$(grep "$counter" "$perf_output_file" | head -1 | \
            sed 's/^[[:space:]]*//' | cut -d' ' -f1 | tr -d ',') || true
        if [[ -n "$val" && "$val" =~ ^[0-9]+$ ]]; then
            echo "$val"
        else
            echo "0"
        fi
    }

    instructions=$(extract_perf_counter "instructions")
    cycles=$(extract_perf_counter "cycles")
    cache_misses=$(extract_perf_counter "cache-misses")
    branch_misses=$(extract_perf_counter "branch-misses")

    rm -f "$perf_output_file"

    # Check for crashes only (not Singular runtime errors, which some tests produce deliberately)
    if echo "$output" | grep -qi 'SIGSEGV\|Segmentation fault\|Abort'; then
        clean="no"
    fi

    echo "${build_name},${testfile},${class},${iters},${run_num},${warmup_cpu_us},${warmup_wall_us},${cpu_us},${wall_us},${ext_wall_ns},${clean},${listfile},${instructions},${cycles},${cache_misses},${branch_misses}"
}

# Load classification from file (keys are filenames, not paths)
load_classification() {
    local file="$1"
    while IFS=',' read -r testname class iters tout; do
        # Strip any path prefix — classification files use filenames only
        local key
        key=$(basename "$testname")
        TEST_CLASS["$key"]="$class"
        TEST_ITERS["$key"]="$iters"
        if [[ -n "$tout" ]]; then
            TEST_TIMEOUT["$key"]="$tout"
        elif [[ "$class" == "C" ]]; then
            TEST_TIMEOUT["$key"]=300
        else
            TEST_TIMEOUT["$key"]=60
        fi
    done < "$file"
}

# Save classification for a specific list file (filenames only)
save_classification_for() {
    local listfile="$1"
    local cfile
    cfile=$(classify_file_for "$listfile")

    # Collect all tests belonging to this list file, write filename only
    {
        while IFS='|' read -r testfile tlistfile; do
            local key
            key=$(classify_key "$testfile")
            if [[ "$tlistfile" == "$listfile" && -n "${TEST_CLASS[$key]+x}" ]]; then
                echo "${key},${TEST_CLASS[$key]},${TEST_ITERS[$key]},${TEST_TIMEOUT[$key]:-300}"
            fi
        done <<< "$(collect_tests)"
    } | sort > "$cfile"

    echo "Saved classification to $cfile" >&2
}

# Main execution
main() {
    local tests
    tests=$(collect_tests)
    local total
    total=$(echo "$tests" | wc -l)

    echo "Benchmarking $total tests × ${#BUILD_NAMES[@]} builds × $RUNS runs" >&2
    for i in "${!BUILD_NAMES[@]}"; do
        echo "  ${BUILD_NAMES[$i]}: ${BUILD_COMMANDS[$i]}" >&2
    done
    echo "Target time: ${TARGET_TIME}s, Timeout: ${TIMEOUT}s" >&2
    echo "" >&2

    # Load explicit classification file if given
    if [[ -n "$CLASSIFY_FROM" && -f "$CLASSIFY_FROM" ]]; then
        echo "Loading classification from $CLASSIFY_FROM" >&2
        load_classification "$CLASSIFY_FROM"
    fi

    # Auto-load classification files for each list file
    for listfile in "${LISTFILES[@]}"; do
        local cfile
        cfile=$(classify_file_for "$listfile")
        if [[ -f "$cfile" ]]; then
            echo "Loading classification from $cfile" >&2
            load_classification "$cfile"
        fi
    done

    # Phase 1: Classify any tests not yet classified
    local need_classify=0
    while IFS='|' read -r testfile listfile; do
        local key
        key=$(classify_key "$testfile")
        if [[ -z "${TEST_CLASS[$key]+x}" ]]; then
            need_classify=1
            break
        fi
    done <<< "$tests"

    if [[ $need_classify -eq 1 ]]; then
        echo "=== Phase 1: Classifying unclassified tests (using ${BUILD_NAMES[0]}) ===" >&2
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

            local key
            key=$(classify_key "$testfile")
            if [[ -n "${TEST_CLASS[$key]+x}" ]]; then
                continue
            fi

            echo -ne "\r  [$count/$total] Classifying: $testfile ... " >&2

            if classify_test "$testfile" "$listfile"; then
                echo "${TEST_CLASS[$testfile]} (N=${TEST_ITERS[$testfile]})" >&2
            else
                echo "FAIL" >&2
            fi
        done <<< "$tests"

        echo "" >&2
        echo "Classification complete." >&2

        # Save classification files if requested
        if [[ $SAVE_CLASSIFY -eq 1 ]]; then
            for listfile in "${!CLASSIFY_DIRTY[@]}"; do
                save_classification_for "$listfile"
            done
        fi
    fi

    # Phase 2: Benchmarking
    echo "=== Phase 2: Benchmarking ===" >&2

    # Output CSV
    local csv_out=""
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

    local total_runs=$(( total * ${#BUILD_NAMES[@]} * RUNS ))
    local run_count=0
    local success=0
    local failed=0

    for run in $(seq 1 $RUNS); do
        for bi in "${!BUILD_NAMES[@]}"; do
            local build_name="${BUILD_NAMES[$bi]}"
            local build_cmd="${BUILD_COMMANDS[$bi]}"

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

                # Classify if not already done
                local key
                key=$(classify_key "$testfile")
                if [[ -z "${TEST_CLASS[$key]+x}" ]]; then
                    classify_test "$testfile" "$listfile" || true
                fi

                local class="${TEST_CLASS[$key]}"

                # Skip excluded tests
                if [[ "$class" == "X" ]]; then
                    continue
                fi

                run_count=$((run_count + 1))
                echo -ne "\r  [$run_count/$total_runs] Run $run, $build_name ($class): $testfile ... " >&2

                local result
                result=$(benchmark_test "$build_name" "$build_cmd" "$testfile" "$run" "$listfile")

                if [[ -n "$csv_out" ]]; then
                    echo "$result" >> "$csv_out"
                else
                    echo "$result"
                fi

                local clean
                clean=$(echo "$result" | cut -d, -f11)
                if [[ "$clean" == "yes" ]]; then
                    success=$((success + 1))
                    echo "done" >&2
                else
                    failed=$((failed + 1))
                    echo "ISSUES" >&2
                fi
            done <<< "$tests"
        done
    done

    echo "" >&2
    echo "=== Summary ===" >&2
    echo "Total runs: $run_count, Success: $success, Issues: $failed" >&2

    # Class breakdown
    local class_a=0 class_b=0 class_c=0 class_d=0 class_e=0 class_x=0 class_fail=0
    for class in "${TEST_CLASS[@]}"; do
        case "$class" in
            A) class_a=$((class_a + 1)) ;;
            B) class_b=$((class_b + 1)) ;;
            C) class_c=$((class_c + 1)) ;;
            D) class_d=$((class_d + 1)) ;;
            E) class_e=$((class_e + 1)) ;;
            X) class_x=$((class_x + 1)) ;;
            FAIL) class_fail=$((class_fail + 1)) ;;
        esac
    done
    echo "Classes: A=$class_a, B=$class_b, C=$class_c, D=$class_d, E=$class_e, X=$class_x, FAIL=$class_fail" >&2

    if [[ -n "$csv_out" ]]; then
        echo "Results written to: $csv_out" >&2
    fi
}

main
