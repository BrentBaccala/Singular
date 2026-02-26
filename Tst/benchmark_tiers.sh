#!/bin/bash
# Tiered benchmark wrapper for comparing Singular builds.
#
# Calls benchmark_singular.sh with appropriate parameters for each tier.
#
# Usage:
#   ./benchmark_tiers.sh --short -s BUILD1 -s BUILD2 ...
#   ./benchmark_tiers.sh --medium -s BUILD1 -s BUILD2 ...
#   ./benchmark_tiers.sh --long -s BUILD1 -s BUILD2 ...
#
# Calibration data (single core, spielwiese, Feb 2026):
#
#   Short tier tests (~30s per build):
#     katsura(9)  Fp std:    ~1s    slimgb: ~1s
#     katsura(8)  QQ std:    ~3s
#     cyclic(7)   Fp std:    ~2s    slimgb: ~2s
#     hydrogen-5  QQ std:    <1s
#
#   Medium tier adds (~5 min per build):
#     katsura(10) Fp std:    ~8s    slimgb: ~5s
#     katsura(9)  QQ std:    ~17s   slimgb: ~31s
#
#   Long tier adds (~35 min per build):
#     katsura(11) Fp std:    ~69s   slimgb: ~38s
#     cyclic(8)   Fp std:    ~41s   slimgb: ~43s
#     katsura(10) QQ std:    ~233s
#     katsura(9)  QQ modstd: ~150s

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_SCRIPT="$SCRIPT_DIR/benchmark_singular.sh"

if [ ! -x "$BENCHMARK_SCRIPT" ]; then
    echo "Error: benchmark_singular.sh not found at $BENCHMARK_SCRIPT"
    exit 1
fi

# Defaults
TIER=""
SINGULAR_ARGS=()
NUM_RUNS=""
EXTRA_ARGS=()

usage() {
    cat << 'EOF'
Usage: benchmark_tiers.sh TIER [OPTIONS] -s BUILD1 [-s BUILD2 ...]

TIERS (exactly one required):
  --short      Quick regression check (~2 min for 4 builds)
  --medium     Meaningful comparison (~20 min for 4 builds)
  --long       Comprehensive overnight benchmark (~2.5 hours for 4 builds)

OPTIONS:
  -s, --singular PATH   Singular build path (can specify multiple)
  -n, --num-runs N      Override default number of runs per tier
  --markdown             Output in markdown format
  --show-output          Show Singular output during runs
  --prot                 Enable protocol output
  --extra "ARGS"         Extra args passed to benchmark_singular.sh
  -h, --help             Show this help

EXAMPLES:
  # Quick regression check with all four builds
  ./benchmark_tiers.sh --short \
    -s ~/Singular/build-spielwiese/Singular/.libs/Singular \
    -s ~/Singular-LSet/build/Singular/.libs/Singular \
    -s ~/Singular-LSet2/build/Singular/.libs/Singular \
    -s ~/Singular-AVX2/build/Singular/.libs/Singular

  # Medium benchmark, 5 runs each
  ./benchmark_tiers.sh --medium -n 5 \
    -s ~/Singular/build-spielwiese/Singular/.libs/Singular \
    -s ~/Singular-LSet2/build/Singular/.libs/Singular

  # Long overnight benchmark
  ./benchmark_tiers.sh --long \
    -s ~/Singular/build-spielwiese/Singular/.libs/Singular \
    -s ~/Singular-LSet/build/Singular/.libs/Singular \
    -s ~/Singular-LSet2/build/Singular/.libs/Singular \
    -s ~/Singular-AVX2/build/Singular/.libs/Singular

TIER DETAILS:
  Short (~30s per build, ~2 min for 4 builds):
    - katsura(9)  Fp(32003) dp, std+slimgb, 3 runs
    - katsura(8)  QQ dp, std, 3 runs
    - cyclic(7)   Fp(32003) dp, std+slimgb, 3 runs
    - hydrogen-5.ssi, std, 3 runs
    Per-test timeout: 30s.  Overall timeout: 5 min.

  Medium (~5 min per build, ~20 min for 4 builds):
    - Everything in short, plus:
    - katsura(10) Fp(32003) dp, std+slimgb, 3 runs
    - katsura(9)  QQ dp, std+slimgb, 3 runs
    Per-test timeout: 120s.  Overall timeout: 30 min.

  Long (~35 min per build, ~2.5 hours for 4 builds):
    - Everything in medium, plus:
    - katsura(11) Fp(32003) dp, std+slimgb, 3 runs
    - cyclic(8)   Fp(32003) dp, std+slimgb, 3 runs
    - katsura(10) QQ dp, std, 3 runs
    - katsura(9)  QQ dp, modstd, 3 runs
    Per-test timeout: 600s.  Overall timeout: 3 hours.
    Note: modstd is multi-threaded; builds using modstd run sequentially.

EOF
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --short)  TIER="short";  shift ;;
        --medium) TIER="medium"; shift ;;
        --long)   TIER="long";   shift ;;
        -s|--singular)
            SINGULAR_ARGS+=("-s" "$2")
            shift 2 ;;
        -n|--num-runs)
            NUM_RUNS="$2"
            shift 2 ;;
        --markdown|--show-output|--prot|--warmup|--memory)
            EXTRA_ARGS+=("$1")
            shift ;;
        --extra)
            # shellcheck disable=SC2206
            EXTRA_ARGS+=($2)
            shift 2 ;;
        -h|--help) usage ;;
        *)
            echo "Unknown argument: $1"
            echo "Use --help for usage information."
            exit 1 ;;
    esac
done

if [ -z "$TIER" ]; then
    echo "Error: Must specify a tier (--short, --medium, or --long)"
    echo "Use --help for usage information."
    exit 1
fi

if [ ${#SINGULAR_ARGS[@]} -eq 0 ]; then
    echo "Error: Must specify at least one Singular build with -s"
    exit 1
fi

# Tier-specific defaults
case $TIER in
    short)
        DEFAULT_RUNS=3
        TEST_TIMEOUT=30
        OVERALL_TIMEOUT=300      # 5 minutes
        ;;
    medium)
        DEFAULT_RUNS=3
        TEST_TIMEOUT=120
        OVERALL_TIMEOUT=1800     # 30 minutes
        ;;
    long)
        DEFAULT_RUNS=3
        TEST_TIMEOUT=600
        OVERALL_TIMEOUT=10800    # 3 hours
        ;;
esac

RUNS="${NUM_RUNS:-$DEFAULT_RUNS}"

# Header
echo "========================================"
echo "  Singular Benchmark - ${TIER^^} tier"
echo "========================================"
echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Machine: $(hostname)"
echo "Runs per test: $RUNS"
echo "Per-test timeout: ${TEST_TIMEOUT}s"
echo "Overall timeout: ${OVERALL_TIMEOUT}s ($(echo "$OVERALL_TIMEOUT / 60" | bc) min)"
echo ""

# Function to run a single benchmark_singular.sh invocation with timeout
run_benchmark() {
    local desc="$1"
    shift
    local args=("$@")

    echo ""
    echo "----------------------------------------"
    echo "  $desc"
    echo "----------------------------------------"

    timeout "$OVERALL_TIMEOUT" "$BENCHMARK_SCRIPT" \
        "${SINGULAR_ARGS[@]}" \
        -n "$RUNS" \
        --timeout "$TEST_TIMEOUT" \
        "${EXTRA_ARGS[@]}" \
        "${args[@]}"
    local exit_code=$?

    if [ $exit_code -eq 124 ]; then
        echo ""
        echo "*** OVERALL TIMEOUT ($OVERALL_TIMEOUT s) reached for: $desc ***"
    fi

    return $exit_code
}

# Track overall start time
OVERALL_START=$(date +%s)

# Run the appropriate tier
case $TIER in
    short)
        # Short tier: fast regression check
        # Test 1: katsura(9) Fp std+slimgb
        run_benchmark "katsura(9) Fp(32003) dp, std+slimgb" \
            --katsura-n 9 --katsura-fp --algorithm std,slimgb || true

        # Test 2: katsura(8) QQ std
        run_benchmark "katsura(8) QQ dp, std" \
            --katsura-n 8 --katsura-qq --algorithm std || true

        # Test 3: cyclic(7) Fp std+slimgb
        run_benchmark "cyclic(7) Fp(32003) dp, std+slimgb" \
            --cyclic-n 7 --cyclic-fp-dp --algorithm std,slimgb || true

        # Test 4: hydrogen-5.ssi std
        if [ -f ~/helium/hydrogen-5.ssi ]; then
            run_benchmark "hydrogen-5.ssi (QQ, 11 vars) std" \
                --ssi ~/helium/hydrogen-5.ssi --algorithm std || true
        fi
        ;;

    medium)
        # Medium tier: meaningful comparison
        # Includes all short tests plus bigger problems
        # Uses warmup runs for more stable timings

        # --- Short tests (repeated for completeness) ---
        run_benchmark "katsura(9) Fp(32003) dp, std+slimgb" \
            --katsura-n 9 --katsura-fp --algorithm std,slimgb --warmup || true

        run_benchmark "katsura(8) QQ dp, std" \
            --katsura-n 8 --katsura-qq --algorithm std --warmup || true

        run_benchmark "cyclic(7) Fp(32003) dp, std+slimgb" \
            --cyclic-n 7 --cyclic-fp-dp --algorithm std,slimgb --warmup || true

        if [ -f ~/helium/hydrogen-5.ssi ]; then
            run_benchmark "hydrogen-5.ssi (QQ, 11 vars) std" \
                --ssi ~/helium/hydrogen-5.ssi --algorithm std --warmup || true
        fi

        # --- Medium-only tests ---
        run_benchmark "katsura(10) Fp(32003) dp, std+slimgb" \
            --katsura-n 10 --katsura-fp --algorithm std,slimgb --warmup || true

        run_benchmark "katsura(9) QQ dp, std+slimgb" \
            --katsura-n 9 --katsura-qq --algorithm std,slimgb --warmup || true
        ;;

    long)
        # Long tier: comprehensive benchmark
        # Includes all medium tests plus heavy problems
        # Uses warmup runs for more stable timings

        # --- Short tests ---
        run_benchmark "katsura(9) Fp(32003) dp, std+slimgb" \
            --katsura-n 9 --katsura-fp --algorithm std,slimgb --warmup || true

        run_benchmark "katsura(8) QQ dp, std" \
            --katsura-n 8 --katsura-qq --algorithm std --warmup || true

        run_benchmark "cyclic(7) Fp(32003) dp, std+slimgb" \
            --cyclic-n 7 --cyclic-fp-dp --algorithm std,slimgb --warmup || true

        if [ -f ~/helium/hydrogen-5.ssi ]; then
            run_benchmark "hydrogen-5.ssi (QQ, 11 vars) std" \
                --ssi ~/helium/hydrogen-5.ssi --algorithm std --warmup || true
        fi

        # --- Medium tests ---
        run_benchmark "katsura(10) Fp(32003) dp, std+slimgb" \
            --katsura-n 10 --katsura-fp --algorithm std,slimgb --warmup || true

        run_benchmark "katsura(9) QQ dp, std+slimgb" \
            --katsura-n 9 --katsura-qq --algorithm std,slimgb --warmup || true

        # --- Long-only tests ---
        run_benchmark "katsura(11) Fp(32003) dp, std+slimgb" \
            --katsura-n 11 --katsura-fp --algorithm std,slimgb --warmup || true

        run_benchmark "cyclic(8) Fp(32003) dp, std+slimgb" \
            --cyclic-n 8 --cyclic-fp-dp --algorithm std,slimgb --warmup || true

        run_benchmark "katsura(10) QQ dp, std" \
            --katsura-n 10 --katsura-qq --algorithm std --warmup || true

        # modstd is multi-threaded — run separately with note
        echo ""
        echo "Note: modstd is multi-threaded. Running sequentially."
        run_benchmark "katsura(9) QQ dp, modstd" \
            --katsura-n 9 --katsura-qq --algorithm modstd --warmup || true
        ;;
esac

# Summary
OVERALL_END=$(date +%s)
OVERALL_ELAPSED=$((OVERALL_END - OVERALL_START))
echo ""
echo "========================================"
echo "  ${TIER^^} tier complete"
echo "  Total wall time: ${OVERALL_ELAPSED}s ($(echo "$OVERALL_ELAPSED / 60" | bc) min)"
echo "========================================"
