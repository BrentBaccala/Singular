#!/bin/bash
# Benchmark script created by Claude Sonnet 4.5 (Anthropic)

# Default Configuration
NUM_RUNS=5
CYCLIC_N=5
KATSURA_N=5
GB_ALGORITHM="std"
V2_NAME="V2-build"
WARMUP_RUN=0
SKIP_V1=0
SHOW_INPUT=0
SHOW_OUTPUT=0
USE_PROT=0
NEWELLP1_FILE=~/Downloads/newellp1
TEMP_DIR=$(mktemp -d)

# Test flags (all enabled by default)
RUN_NEWELLP1=0
RUN_CYCLIC_QQ_DP=0
RUN_CYCLIC_ZZ_DP=0
RUN_CYCLIC_QQ_LP=0
RUN_CYCLIC_ZZ_LP=0
RUN_CYCLIC_HOM_QQ=0
RUN_CYCLIC_HOM_ZZ=0
RUN_KATSURA_QQ=0
RUN_KATSURA_ZZ=0

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS] [TESTS]

OPTIONS:
  -n, --num-runs N          Number of runs per test (default: 5)
  --cyclic-n N              Number of variables for cyclic tests (default: 5)
  --katsura-n N             Number of variables for katsura tests (default: 5)
  --algorithm ALG           Groebner basis algorithm to use (default: std)
                            Options: std, modstd, groebner, slimgb, hilb, fglm
  --v2-name NAME            Name for V2 version in output (default: V2-build)
  --warmup                  Perform an untimed warmup run before timed runs
  --skip-v1                 Only run Singular-build, skip Singular-spielwiese
  --show-input              Display the input file content before each test
  --show-output             Display the output from each run
  --prot                    Enable option(prot) for protocol output during computation
  -h, --help                Show this help message
  -a, --all                 Run all available tests (default if no tests specified)

TESTS (can specify multiple):
  --newellp1               Run newellp1 test
  --cyclic-qq-dp           Run cyclic(5) with QQ coefficients and dp ordering
  --cyclic-zz-dp           Run cyclic(5) with ZZ/32003 and dp ordering
  --cyclic-qq-lp           Run cyclic(5) with QQ and lp ordering
  --cyclic-zz-lp           Run cyclic(5) with ZZ/32003 and lp ordering
  --cyclic-hom-qq          Run homogenized cyclic(5) with QQ and dp
  --cyclic-hom-zz          Run homogenized cyclic(5) with ZZ/32003 and dp
  --katsura-qq             Run katsura(5) with QQ and dp
  --katsura-zz             Run katsura(5) with ZZ/32003 and dp
  --cyclic                 Run all cyclic tests
  --katsura                Run all katsura tests

EXAMPLES:
  $0 --num-runs 10 --cyclic-qq-dp --katsura-qq
  $0 -n 3 --cyclic --katsura
  $0 --all
  $0 --newellp1 --cyclic-qq-dp
  $0 --cyclic-n 6 --cyclic
  $0 --cyclic-n 7 --katsura-n 6 --all
  $0 --algorithm modstd --cyclic
  $0 --algorithm slimgb --cyclic-n 7 --cyclic-qq-dp
  $0 --v2-name V2-pr1301 --cyclic
  $0 --warmup --cyclic-qq-dp
  $0 --skip-v1 --v2-name V2-pr1301 --cyclic
  $0 --show-input --cyclic-qq-dp
  $0 --show-output --cyclic-qq-dp
  $0 --prot --show-output --cyclic-qq-dp

EOF
    exit 0
}

# Parse command line arguments
TESTS_SPECIFIED=0
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-runs)
            NUM_RUNS="$2"
            shift 2
            ;;
        --cyclic-n)
            CYCLIC_N="$2"
            shift 2
            ;;
        --katsura-n)
            KATSURA_N="$2"
            shift 2
            ;;
        --algorithm)
            GB_ALGORITHM="$2"
            case $GB_ALGORITHM in
                std|modstd|groebner|slimgb|hilb|fglm)
                    ;;
                *)
                    echo "Error: Unknown algorithm '$GB_ALGORITHM'"
                    echo "Valid options: std, modstd, groebner, slimgb, hilb, fglm"
                    exit 1
                    ;;
            esac
            # Map modstd to modStd for Singular
            if [ "$GB_ALGORITHM" == "modstd" ]; then
                GB_ALGORITHM_CALL="modStd"
            else
                GB_ALGORITHM_CALL="$GB_ALGORITHM"
            fi
            shift 2
            ;;
        --v2-name)
            V2_NAME="$2"
            shift 2
            ;;
        --warmup)
            WARMUP_RUN=1
            shift
            ;;
        --skip-v1)
            SKIP_V1=1
            shift
            ;;
        --show-input)
            SHOW_INPUT=1
            shift
            ;;
        --show-output)
            SHOW_OUTPUT=1
            shift
            ;;
        --prot)
            USE_PROT=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        -a|--all)
            RUN_NEWELLP1=1
            RUN_CYCLIC_QQ_DP=1
            RUN_CYCLIC_ZZ_DP=1
            RUN_CYCLIC_QQ_LP=1
            RUN_CYCLIC_ZZ_LP=1
            RUN_CYCLIC_HOM_QQ=1
            RUN_CYCLIC_HOM_ZZ=1
            RUN_KATSURA_QQ=1
            RUN_KATSURA_ZZ=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --newellp1)
            RUN_NEWELLP1=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --cyclic-qq-dp)
            RUN_CYCLIC_QQ_DP=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --cyclic-zz-dp)
            RUN_CYCLIC_ZZ_DP=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --cyclic-qq-lp)
            RUN_CYCLIC_QQ_LP=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --cyclic-zz-lp)
            RUN_CYCLIC_ZZ_LP=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --cyclic-hom-qq)
            RUN_CYCLIC_HOM_QQ=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --cyclic-hom-zz)
            RUN_CYCLIC_HOM_ZZ=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --katsura-qq)
            RUN_KATSURA_QQ=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --katsura-zz)
            RUN_KATSURA_ZZ=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --cyclic)
            RUN_CYCLIC_QQ_DP=1
            RUN_CYCLIC_ZZ_DP=1
            RUN_CYCLIC_QQ_LP=1
            RUN_CYCLIC_ZZ_LP=1
            RUN_CYCLIC_HOM_QQ=1
            RUN_CYCLIC_HOM_ZZ=1
            TESTS_SPECIFIED=1
            shift
            ;;
        --katsura)
            RUN_KATSURA_QQ=1
            RUN_KATSURA_ZZ=1
            TESTS_SPECIFIED=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# If no tests specified, run all
if [ $TESTS_SPECIFIED -eq 0 ]; then
    RUN_NEWELLP1=1
    RUN_CYCLIC_QQ_DP=1
    RUN_CYCLIC_ZZ_DP=1
    RUN_CYCLIC_QQ_LP=1
    RUN_CYCLIC_ZZ_LP=1
    RUN_CYCLIC_HOM_QQ=1
    RUN_CYCLIC_HOM_ZZ=1
    RUN_KATSURA_QQ=1
    RUN_KATSURA_ZZ=1
fi

echo "Comprehensive Singular Benchmark"
echo "================================"
echo "Number of runs per test: $NUM_RUNS"
echo "Cyclic n: $CYCLIC_N"
echo "Katsura n: $KATSURA_N"
echo "Algorithm: $GB_ALGORITHM"
if [ "$GB_ALGORITHM" == "modstd" ]; then
    echo "  (calling as: modStd)"
fi
echo "V2 name: $V2_NAME"
echo "Warmup run: $([ $WARMUP_RUN -eq 1 ] && echo 'enabled' || echo 'disabled')"
echo "Skip V1: $([ $SKIP_V1 -eq 1 ] && echo 'yes' || echo 'no')"
echo "Show input: $([ $SHOW_INPUT -eq 1 ] && echo 'yes' || echo 'no')"
echo "Show output: $([ $SHOW_OUTPUT -eq 1 ] && echo 'yes' || echo 'no')"
echo "Protocol output: $([ $USE_PROT -eq 1 ] && echo 'enabled' || echo 'disabled')"
echo "Temporary directory: $TEMP_DIR"
echo ""

# Cleanup on exit
trap "rm -rf $TEMP_DIR" EXIT

# Function to create test input files
create_test_file() {
    local test_name=$1
    local filename="$TEMP_DIR/${test_name}.sing"
    
    # Generate variable list for cyclic
    local cyclic_vars=""
    for ((i=1; i<=$CYCLIC_N; i++)); do
        if [ $i -eq 1 ]; then
            cyclic_vars="x$i"
        else
            cyclic_vars="${cyclic_vars},x$i"
        fi
    done
    
    # Generate variable list for katsura (0-indexed)
    local katsura_vars=""
    for ((i=0; i<=$KATSURA_N; i++)); do
        if [ $i -eq 0 ]; then
            katsura_vars="x$i"
        else
            katsura_vars="${katsura_vars},x$i"
        fi
    done
    
    case $test_name in
        cyclic_qq_dp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 0,($cyclic_vars),dp;
ideal i = cyclic($CYCLIC_N);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
        cyclic_zz_dp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 32003,($cyclic_vars),dp;
ideal i = cyclic($CYCLIC_N);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
        cyclic_qq_lp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 0,($cyclic_vars),lp;
ideal i = cyclic($CYCLIC_N);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
        cyclic_zz_lp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 32003,($cyclic_vars),lp;
ideal i = cyclic($CYCLIC_N);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
        cyclic_hom_qq_dp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 0,($cyclic_vars,h),dp;
ideal i = homog(cyclic($CYCLIC_N),h);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
        cyclic_hom_zz_dp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 32003,($cyclic_vars,h),dp;
ideal i = homog(cyclic($CYCLIC_N),h);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
        katsura_qq_dp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 0,($katsura_vars),dp;
ideal i = katsura($KATSURA_N);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
        katsura_zz_dp)
            cat > "$filename" << EOF
LIB "polylib.lib";
LIB "modstd.lib";
ring r = 32003,($katsura_vars),dp;
ideal i = katsura($KATSURA_N);
$([ $USE_PROT -eq 1 ] && echo "option(prot);")
option(redSB);
ideal j = $GB_ALGORITHM_CALL(i);
quit;
EOF
            ;;
    esac
    
    echo "$filename"
}

# Function to run benchmark
run_benchmark() {
    local version_name=$1
    local ld_library_path=$2
    local singular_path=$3
    local executable=$4
    local test_name=$5
    local input_file=$6
    
    echo -e "${BLUE}Testing $version_name - $test_name${NC}"
    echo "----------------------------------------"
    
    # Show input file if requested
    if [ $SHOW_INPUT -eq 1 ]; then
        echo -e "${YELLOW}Input file ($input_file):${NC}"
        echo "- - - - - - - - - - - - - - - - - - - -"
        cat "$input_file"
        echo "- - - - - - - - - - - - - - - - - - - -"
        echo ""
    fi
    
    # Perform warmup run if requested
    if [ $WARMUP_RUN -eq 1 ]; then
        echo -n "Warmup run... "
        if [ $SHOW_OUTPUT -eq 1 ]; then
            echo ""
            echo -e "${YELLOW}Warmup output:${NC}"
            echo "- - - - - - - - - - - - - - - - - - - -"
            LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" 2>&1
            echo "- - - - - - - - - - - - - - - - - - - -"
        else
            LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" > /dev/null 2>&1
            echo "done"
        fi
    fi
    
    local times=()
    local total=0
    
    for i in $(seq 1 $NUM_RUNS); do
        echo -n "Run $i/$NUM_RUNS... "
        
        # Run with time measurement
        local start=$(date +%s.%N)
        if [ $SHOW_OUTPUT -eq 1 ]; then
            # Save output to temp file for timing purposes
            local output_file="$TEMP_DIR/output_${version_name}_${test_name}_run${i}.txt"
            LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" | tee "$output_file" 2>&1
        else
            LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" > /dev/null 2>&1
        fi
        local end=$(date +%s.%N)
        
        local runtime=$(echo "$end - $start" | bc)
        times+=($runtime)
        total=$(echo "$total + $runtime" | bc)
        
        echo "${runtime}s"
        
        # Show output if requested
        if [ $SHOW_OUTPUT -eq 1 ]; then
            echo -e "${YELLOW}Output from run $i:${NC}"
            echo "- - - - - - - - - - - - - - - - - - - -"
            cat "$output_file"
            echo "- - - - - - - - - - - - - - - - - - - -"
            echo ""
        fi
    done
    
    # Calculate statistics
    local avg=$(echo "scale=4; $total / $NUM_RUNS" | bc)
    
    # Find min and max
    local min=${times[0]}
    local max=${times[0]}
    for time in "${times[@]}"; do
        if (( $(echo "$time < $min" | bc -l) )); then
            min=$time
        fi
        if (( $(echo "$time > $max" | bc -l) )); then
            max=$time
        fi
    done
    
    # Calculate standard deviation
    local sum_sq_diff=0
    for time in "${times[@]}"; do
        local diff=$(echo "$time - $avg" | bc)
        local sq_diff=$(echo "$diff * $diff" | bc)
        sum_sq_diff=$(echo "$sum_sq_diff + $sq_diff" | bc)
    done
    local variance=$(echo "scale=6; $sum_sq_diff / $NUM_RUNS" | bc)
    local stddev=$(echo "scale=4; sqrt($variance)" | bc)
    
    echo ""
    echo -e "${GREEN}Results:${NC}"
    echo "  Average: ${avg}s"
    echo "  Std Dev: ${stddev}s"
    echo "  Min:     ${min}s"
    echo "  Max:     ${max}s"
    echo "  Total:   ${total}s"
    echo ""
    
    # Store results for summary
    echo "$version_name|$test_name|$avg|$stddev|$min|$max|$total" >> "$TEMP_DIR/results.txt"
}

# Version 1: Singular-spielwiese
LD_PATH_V1=$(find ~/src/Singular-spielwiese/ -name .libs | tr '\n' ' ' | sed 's/[[:space:]]/:/g')
SING_PATH_V1=$(echo ~/src/Singular-spielwiese/Singular/dyn_modules/*/.libs | sed 's/ /:/g')
EXEC_V1=~/src/Singular-spielwiese/Singular/.libs/Singular

# Version 2: Singular-build
LD_PATH_V2=$(find ~/src/Singular-build -name .libs | tr '\n' ' ' | sed 's/[[:space:]]/:/g')
SING_PATH_V2=$(echo ~/src/Singular-build/Singular/dyn_modules/*/.libs | sed 's/ /:/g')
EXEC_V2=~/src/Singular-build/Singular/.libs/Singular

# Initialize results file
echo "Version|Test|Average|StdDev|Min|Max|Total" > "$TEMP_DIR/results.txt"

# Build list of tests to run
declare -a TESTS_TO_RUN

if [ $RUN_NEWELLP1 -eq 1 ] && [ -f "$NEWELLP1_FILE" ]; then
    TESTS_TO_RUN+=("newellp1:$NEWELLP1_FILE")
fi

if [ $RUN_CYCLIC_QQ_DP -eq 1 ]; then
    TESTS_TO_RUN+=("cyclic${CYCLIC_N}-QQ-dp:$(create_test_file cyclic_qq_dp)")
fi

if [ $RUN_CYCLIC_ZZ_DP -eq 1 ]; then
    TESTS_TO_RUN+=("cyclic${CYCLIC_N}-ZZ-dp:$(create_test_file cyclic_zz_dp)")
fi

if [ $RUN_CYCLIC_QQ_LP -eq 1 ]; then
    TESTS_TO_RUN+=("cyclic${CYCLIC_N}-QQ-lp:$(create_test_file cyclic_qq_lp)")
fi

if [ $RUN_CYCLIC_ZZ_LP -eq 1 ]; then
    TESTS_TO_RUN+=("cyclic${CYCLIC_N}-ZZ-lp:$(create_test_file cyclic_zz_lp)")
fi

if [ $RUN_CYCLIC_HOM_QQ -eq 1 ]; then
    TESTS_TO_RUN+=("cyclic${CYCLIC_N}-hom-QQ-dp:$(create_test_file cyclic_hom_qq_dp)")
fi

if [ $RUN_CYCLIC_HOM_ZZ -eq 1 ]; then
    TESTS_TO_RUN+=("cyclic${CYCLIC_N}-hom-ZZ-dp:$(create_test_file cyclic_hom_zz_dp)")
fi

if [ $RUN_KATSURA_QQ -eq 1 ]; then
    TESTS_TO_RUN+=("katsura${KATSURA_N}-QQ-dp:$(create_test_file katsura_qq_dp)")
fi

if [ $RUN_KATSURA_ZZ -eq 1 ]; then
    TESTS_TO_RUN+=("katsura${KATSURA_N}-ZZ-dp:$(create_test_file katsura_zz_dp)")
fi

if [ ${#TESTS_TO_RUN[@]} -eq 0 ]; then
    echo -e "${RED}No tests selected to run!${NC}"
    echo "Use --help to see available options"
    exit 1
fi

echo -e "${YELLOW}Tests to run: ${#TESTS_TO_RUN[@]}${NC}"
echo ""

# Run all benchmarks for Version 1
if [ $SKIP_V1 -eq 0 ]; then
    echo -e "${YELLOW}=== Running benchmarks for Version 1 (spielwiese) ===${NC}"
    echo ""

    for test_spec in "${TESTS_TO_RUN[@]}"; do
        IFS=':' read -r test_name test_file <<< "$test_spec"
        run_benchmark "V1-spielwiese" "$LD_PATH_V1" "$SING_PATH_V1" "$EXEC_V1" "$test_name" "$test_file"
    done
else
    echo -e "${YELLOW}=== Skipping Version 1 (spielwiese) ===${NC}"
    echo ""
fi

# Run all benchmarks for Version 2
echo -e "${YELLOW}=== Running benchmarks for $V2_NAME ===${NC}"
echo ""

for test_spec in "${TESTS_TO_RUN[@]}"; do
    IFS=':' read -r test_name test_file <<< "$test_spec"
    run_benchmark "$V2_NAME" "$LD_PATH_V2" "$SING_PATH_V2" "$EXEC_V2" "$test_name" "$test_file"
done

# Print summary
echo "================================"
echo -e "${YELLOW}SUMMARY OF ALL BENCHMARKS${NC}"
echo "================================"
echo ""

# Parse results and create comparison
declare -A v1_times
declare -A v2_times
declare -A v1_stddev
declare -A v2_stddev

while IFS='|' read -r version test avg stddev min max total; do
    if [ "$version" != "Version" ]; then
        if [ "$version" == "V1-spielwiese" ]; then
            v1_times["$test"]="$avg"
            v1_stddev["$test"]="$stddev"
        elif [ "$version" == "$V2_NAME" ]; then
            v2_times["$test"]="$avg"
            v2_stddev["$test"]="$stddev"
        fi
    fi
done < "$TEMP_DIR/results.txt"

# Print comparison table (text format)
if [ $SKIP_V1 -eq 1 ]; then
    # Only V2 results
    echo -e "${BLUE}Results for $V2_NAME:${NC}"
    printf "%-30s %-20s\n" "Test" "Time"
    printf "%-30s %-20s\n" "----" "----"
    
    for test in "${!v2_times[@]}"; do
        v2="${v2_times[$test]}"
        v2_sd="${v2_stddev[$test]}"
        printf "%-30s %-20s\n" "$test" "${v2}s (±${v2_sd})"
    done
else
    # Comparison table
    echo -e "${BLUE}Average Times Comparison:${NC}"
    printf "%-30s %-25s %-25s\n" "Test" "V1-spielwiese" "$V2_NAME"
    printf "%-30s %-25s %-25s\n" "----" "-------------" "--------"
    
    for test in "${!v1_times[@]}"; do
        v1="${v1_times[$test]}"
        v2="${v2_times[$test]}"
        v1_sd="${v1_stddev[$test]}"
        v2_sd="${v2_stddev[$test]}"
        
        printf "%-30s %-25s %-25s" "$test" "${v1}s (±${v1_sd})" "${v2}s (±${v2_sd})"
        
        # Calculate and show speedup
        if [ -n "$v1" ] && [ -n "$v2" ]; then
            speedup=$(echo "scale=2; $v1 / $v2" | bc)
            if (( $(echo "$speedup > 1.05" | bc -l) )); then
                echo -e " ${GREEN}(V2 ${speedup}x faster)${NC}"
            elif (( $(echo "$speedup < 0.95" | bc -l) )); then
                speedup=$(echo "scale=2; $v2 / $v1" | bc)
                echo -e " ${RED}(V1 ${speedup}x faster)${NC}"
            else
                echo " (similar)"
            fi
        else
            echo ""
        fi
    done
fi

echo ""
echo "================================"

# Now create markdown format
MARKDOWN_FILE="$TEMP_DIR/results.md"

if [ $SKIP_V1 -eq 1 ]; then
    # Only V2 markdown
    cat > "$MARKDOWN_FILE" << MDEOF
## Results for $V2_NAME

| Test | Time |
|------|------|
MDEOF

    # Sort tests for consistent output
    readarray -t sorted_tests < <(printf '%s\n' "${!v2_times[@]}" | sort)
    
    for test in "${sorted_tests[@]}"; do
        v2="${v2_times[$test]}"
        v2_sd="${v2_stddev[$test]}"
        
        # Format numbers with leading zeros
        v2_fmt=$(printf "%.4f" "$v2")
        v2_sd_fmt=$(printf "%.4f" "$v2_sd")
        
        echo "| $test | ${v2_fmt}s (±${v2_sd_fmt}) |" >> "$MARKDOWN_FILE"
    done
else
    # Comparison markdown
    cat > "$MARKDOWN_FILE" << MDEOF
## Average Times Comparison

| Test | V1-spielwiese | $V2_NAME | Comparison |
|------|---------------|----------|------------|
MDEOF

    # Sort tests for consistent output
    readarray -t sorted_tests < <(printf '%s\n' "${!v1_times[@]}" | sort)
    
    for test in "${sorted_tests[@]}"; do
        v1="${v1_times[$test]}"
        v2="${v2_times[$test]}"
        v1_sd="${v1_stddev[$test]}"
        v2_sd="${v2_stddev[$test]}"
        
        # Format numbers with leading zeros
        v1_fmt=$(printf "%.4f" "$v1")
        v2_fmt=$(printf "%.4f" "$v2")
        v1_sd_fmt=$(printf "%.4f" "$v1_sd")
        v2_sd_fmt=$(printf "%.4f" "$v2_sd")
        
        # Calculate comparison
        comparison=""
        if [ -n "$v1" ] && [ -n "$v2" ]; then
            speedup=$(echo "scale=2; $v1 / $v2" | bc)
            if (( $(echo "$speedup > 1.05" | bc -l) )); then
                comparison="V2 ${speedup}x faster"
            elif (( $(echo "$speedup < 0.95" | bc -l) )); then
                speedup=$(echo "scale=2; $v2 / $v1" | bc)
                comparison="V1 ${speedup}x faster"
            else
                comparison="similar"
            fi
        fi
        
        echo "| $test | ${v1_fmt}s (±${v1_sd_fmt}) | ${v2_fmt}s (±${v2_sd_fmt}) | $comparison |" >> "$MARKDOWN_FILE"
    done
fi

echo ""
echo -e "${YELLOW}Markdown Output:${NC}"
echo "================================"
cat "$MARKDOWN_FILE"
echo "================================"
echo ""
echo "Detailed results saved to: $TEMP_DIR/results.txt"
echo "Markdown table saved to: $MARKDOWN_FILE"
echo "Benchmark Complete!"
