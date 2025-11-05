#!/bin/bash
# Benchmark script created by Claude Sonnet 4.5 (Anthropic)

# Default Configuration
NUM_RUNS=5
CYCLIC_N=5
KATSURA_N=5
GB_ALGORITHMS=("std")
WARMUP_RUN=0
SHOW_INPUT=0
SHOW_OUTPUT=0
USE_PROT=0
SHOW_STRATEGY=0
TEMP_DIR=$(mktemp -d)

# Arrays to store Singular executables and their info
declare -a SINGULAR_EXECS
declare -a SINGULAR_NAMES
declare -a SINGULAR_LD_PATHS
declare -a SINGULAR_PATHS

# Test flags (all disabled by default)
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
  --algorithm ALG[,ALG2,...]  Groebner basis algorithm(s) to use (default: std)
                            Options: std, modstd, groebner, slimgb, sba, all
                            'all' runs: std, modstd, groebner, slimgb, sba
                            Can specify multiple separated by commas
  -s, --singular PATH       Path to Singular executable (can be specified multiple times)
                            If not specified, uses 'Singular' from PATH
  --warmup                  Perform an untimed warmup run before timed runs
  --show-input              Display the input file content before each test
  --show-output             Display the output from each run (shows live output)
  --prot                    Enable option(prot) for protocol output during computation
  --show-strategy           Enable option bit 23 to display strategy information
  -h, --help                Show this help message
  -a, --all                 Run all available tests (default if no tests specified)

TESTS (can specify multiple):
  --newellp1               Run newellp1 test (embedded in script)
  --cyclic-qq-dp           Run cyclic(n) with QQ coefficients and dp ordering
  --cyclic-zz-dp           Run cyclic(n) with ZZ/32003 and dp ordering
  --cyclic-qq-lp           Run cyclic(n) with QQ and lp ordering
  --cyclic-zz-lp           Run cyclic(n) with ZZ/32003 and lp ordering
  --cyclic-hom-qq          Run homogenized cyclic(n) with QQ and dp
  --cyclic-hom-zz          Run homogenized cyclic(n) with ZZ/32003 and dp
  --katsura-qq             Run katsura(n) with QQ and dp
  --katsura-zz             Run katsura(n) with ZZ/32003 and dp
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
  $0 --warmup --cyclic-qq-dp
  $0 --show-output --cyclic-qq-dp
  $0 --prot --show-output --cyclic-qq-dp
  $0 --show-strategy --show-output --cyclic-qq-dp
  $0 --algorithm std,modstd,sba --cyclic-qq-dp
  $0 --algorithm all --cyclic-qq-dp
  $0 -s ~/src/Singular-build/Singular/.libs/Singular --cyclic
  $0 -s ~/build1/Singular -s ~/build2/Singular --algorithm std,sba --cyclic

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
            # Handle 'all' keyword
            if [ "$2" == "all" ]; then
                GB_ALGORITHMS=("std" "modstd" "groebner" "slimgb" "sba")
            else
                IFS=',' read -ra GB_ALGORITHMS <<< "$2"
                # Validate each algorithm
                for alg in "${GB_ALGORITHMS[@]}"; do
                    case $alg in
                        std|modstd|groebner|slimgb|sba)
                            ;;
                        all)
                            # Expand 'all' in a comma-separated list
                            echo "Error: 'all' should be used alone, not in a comma-separated list"
                            exit 1
                            ;;
                        *)
                            echo "Error: Unknown algorithm '$alg'"
                            echo "Valid options: std, modstd, groebner, slimgb, sba, all"
                            exit 1
                            ;;
                    esac
                done
            fi
            shift 2
            ;;
        -s|--singular)
            # Check if it's a directory or file
            if [ -d "$2" ]; then
                # It's a directory - look for Singular executable
                if [ -f "$2/Singular/.libs/Singular" ]; then
                    SINGULAR_EXECS+=("$2/Singular/.libs/Singular")
                elif [ -f "$2/Singular" ]; then
                    SINGULAR_EXECS+=("$2/Singular")
                else
                    echo "Error: Could not find Singular executable in directory: $2"
                    echo "Looked for: $2/Singular/.libs/Singular or $2/Singular"
                    exit 1
                fi
            elif [ -f "$2" ]; then
                # It's a file
                SINGULAR_EXECS+=("$2")
            else
                echo "Error: Path does not exist: $2"
                exit 1
            fi
            shift 2
            ;;
        --warmup)
            WARMUP_RUN=1
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
        --show-strategy)
            SHOW_STRATEGY=1
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

# If no Singular executables specified, use default from PATH
if [ ${#SINGULAR_EXECS[@]} -eq 0 ]; then
    SINGULAR_EXECS=("Singular")
fi

# Process each Singular executable to determine paths and names
for exec_path in "${SINGULAR_EXECS[@]}"; do
    # Get absolute path if it exists
    if [ -f "$exec_path" ]; then
        exec_path=$(realpath "$exec_path")
    fi
    
    # Determine name for this Singular version
    if [ "$exec_path" == "Singular" ]; then
        # Using Singular from PATH
        name="Singular"
        ld_path=""
        sing_path=""
    else
        # Extract a meaningful name from the path
        # Check if it's in a .libs directory (built but not installed)
        if [[ "$exec_path" == *"/.libs/"* ]]; then
            # It's a built-but-not-installed version
            # Go up from .libs/Singular to the Singular directory, then up to the build directory
            build_dir=$(dirname $(dirname $(dirname "$exec_path")))
            name=$(basename "$build_dir")
            
            # Find all .libs directories
            ld_path=$(find "$build_dir" -name .libs -type d 2>/dev/null | tr '\n' ':' | sed 's/:$//')
            
            # Find dyn_modules/.libs
            sing_path=$(find "$build_dir" -path "*/Singular/dyn_modules/*/.libs" -type d 2>/dev/null | tr '\n' ':' | sed 's/:$//')
        else
            # Assume it's an installed version or custom location
            name=$(basename $(dirname "$exec_path"))
            ld_path=""
            sing_path=""
        fi
    fi
    
    SINGULAR_NAMES+=("$name")
    SINGULAR_LD_PATHS+=("$ld_path")
    SINGULAR_PATHS+=("$sing_path")
done

echo "Comprehensive Singular Benchmark"
echo "================================"
echo "Number of runs per test: $NUM_RUNS"
echo "Cyclic n: $CYCLIC_N"
echo "Katsura n: $KATSURA_N"
echo "Algorithms: ${GB_ALGORITHMS[*]}"
echo "Singular versions: ${#SINGULAR_EXECS[@]}"
for i in "${!SINGULAR_EXECS[@]}"; do
    echo "  [$((i+1))] ${SINGULAR_NAMES[$i]}: ${SINGULAR_EXECS[$i]}"
done
echo "Warmup run: $([ $WARMUP_RUN -eq 1 ] && echo 'enabled' || echo 'disabled')"
echo "Show input: $([ $SHOW_INPUT -eq 1 ] && echo 'yes' || echo 'no')"
echo "Show output: $([ $SHOW_OUTPUT -eq 1 ] && echo 'yes' || echo 'no')"
echo "Protocol output: $([ $USE_PROT -eq 1 ] && echo 'enabled' || echo 'disabled')"
echo "Show strategy: $([ $SHOW_STRATEGY -eq 1 ] && echo 'enabled' || echo 'disabled')"
echo "Temporary directory: $TEMP_DIR"
echo ""

# Cleanup on exit
trap "rm -rf $TEMP_DIR" EXIT

# Function to map algorithm name to Singular function call
get_algorithm_call() {
    local alg=$1
    if [ "$alg" == "modstd" ]; then
        echo "modStd"
    else
        echo "$alg"
    fi
}

# Function to create test input files
create_test_file() {
    local test_name=$1
    local algorithm=$2
    local filename="$TEMP_DIR/${test_name}_${algorithm}.sing"
    local GB_ALGORITHM_CALL=$(get_algorithm_call "$algorithm")
    
    # Handle newellp1 specially
    if [ "$test_name" == "newellp1" ]; then
        cat > "$filename" << 'EOF'
LIB "modstd.lib";
ring R=QQ, (x,y,z,u,v),M(0,0,0,1,1, 1,1,1,0,0, 1,1,0,0,0, 1,0,0,0,0, 0,0,0,1,0);
ideal I=  -x + 7/5 - 231/125 * v^2 + 39/80 * u^2 - 1/5 * u^3 + 99/400 * u * v^2 - 1287/2000 * u^2 * v^2 + 33/125 * u^3 * v^2 - 3/16 * u + 56/125 * v^3 - 3/50 * u * v^3 + 39/250 * u^2 * v^3 - 8/125 * u^3 * v^3,
-y + 63/125 * v^2 - 294/125 * v + 56/125 * v^3 - 819/1000 * u^2 * v + 42/125 * u^3 * v - 3/50 * u * v^3 + 351/2000 * u^2 * v^2 + 39/250 * u^2 * v^3 - 9/125 * u^3 * v^2 - 8/125 * u^3 * v^3,
-z + 12/5 - 63/160 * u^2 + 63/160 * u;
EOF
        # Add optional settings
        [ $SHOW_STRATEGY -eq 1 ] && cat >> "$filename" << 'EOF'
intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);
EOF
        [ $USE_PROT -eq 1 ] && echo "option(prot);" >> "$filename"
        
        # Add computation
        cat >> "$filename" << EOF
int t=timer;
ideal J=$GB_ALGORITHM_CALL(I);
timer-t;
quit;
EOF
        echo "$filename"
        return
    fi
    
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
$([ $SHOW_STRATEGY -eq 1 ] && echo 'intvec options = option(get);
options[2] = options[2] + 2^23;
option(set, options);')
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
            if [ -n "$ld_library_path" ]; then
                LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" 2>&1 | tee "$TEMP_DIR/warmup_output.txt"
            else
                "$executable" < "$input_file" 2>&1 | tee "$TEMP_DIR/warmup_output.txt"
            fi
            echo "- - - - - - - - - - - - - - - - - - - -"
        else
            if [ -n "$ld_library_path" ]; then
                LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" > /dev/null 2>&1
            else
                "$executable" < "$input_file" > /dev/null 2>&1
            fi
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
            # Save output to temp file and display with tee
            local output_file="$TEMP_DIR/output_${version_name}_${test_name}_run${i}.txt"
            echo ""
            echo -e "${YELLOW}Output from run $i:${NC}"
            echo "- - - - - - - - - - - - - - - - - - - -"
            if [ -n "$ld_library_path" ]; then
                LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" 2>&1 | tee "$output_file"
            else
                "$executable" < "$input_file" 2>&1 | tee "$output_file"
            fi
            echo "- - - - - - - - - - - - - - - - - - - -"
        else
            if [ -n "$ld_library_path" ]; then
                LD_LIBRARY_PATH="$ld_library_path" SINGULARPATH="$singular_path" "$executable" < "$input_file" > /dev/null 2>&1
            else
                "$executable" < "$input_file" > /dev/null 2>&1
            fi
        fi
        local end=$(date +%s.%N)
        
        local runtime=$(echo "$end - $start" | bc)
        times+=($runtime)
        total=$(echo "$total + $runtime" | bc)
        
        echo "Time: ${runtime}s"
        echo ""
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
    
    if [ $NUM_RUNS -gt 1 ]; then
        echo ""
        echo -e "${GREEN}Results:${NC}"
        echo "  Average: ${avg}s"
        echo "  Std Dev: ${stddev}s"
        echo "  Min:     ${min}s"
        echo "  Max:     ${max}s"
        echo "  Total:   ${total}s"
        echo ""
    fi
    
    # Store results for summary
    echo "$version_name|$test_name|$avg|$stddev|$min|$max|$total" >> "$TEMP_DIR/results.txt"
}

# Initialize results file
echo "Version|Test|Average|StdDev|Min|Max|Total" > "$TEMP_DIR/results.txt"

# Build list of tests to run
declare -a TESTS_TO_RUN

# For each algorithm, create test files
for algorithm in "${GB_ALGORITHMS[@]}"; do
    if [ $RUN_NEWELLP1 -eq 1 ]; then
        TESTS_TO_RUN+=("newellp1-${algorithm}:$(create_test_file newellp1 $algorithm)")
    fi

    if [ $RUN_CYCLIC_QQ_DP -eq 1 ]; then
        TESTS_TO_RUN+=("cyclic${CYCLIC_N}-QQ-dp-${algorithm}:$(create_test_file cyclic_qq_dp $algorithm)")
    fi

    if [ $RUN_CYCLIC_ZZ_DP -eq 1 ]; then
        TESTS_TO_RUN+=("cyclic${CYCLIC_N}-ZZ-dp-${algorithm}:$(create_test_file cyclic_zz_dp $algorithm)")
    fi

    if [ $RUN_CYCLIC_QQ_LP -eq 1 ]; then
        TESTS_TO_RUN+=("cyclic${CYCLIC_N}-QQ-lp-${algorithm}:$(create_test_file cyclic_qq_lp $algorithm)")
    fi

    if [ $RUN_CYCLIC_ZZ_LP -eq 1 ]; then
        TESTS_TO_RUN+=("cyclic${CYCLIC_N}-ZZ-lp-${algorithm}:$(create_test_file cyclic_zz_lp $algorithm)")
    fi

    if [ $RUN_CYCLIC_HOM_QQ -eq 1 ]; then
        TESTS_TO_RUN+=("cyclic${CYCLIC_N}-hom-QQ-dp-${algorithm}:$(create_test_file cyclic_hom_qq_dp $algorithm)")
    fi

    if [ $RUN_CYCLIC_HOM_ZZ -eq 1 ]; then
        TESTS_TO_RUN+=("cyclic${CYCLIC_N}-hom-ZZ-dp-${algorithm}:$(create_test_file cyclic_hom_zz_dp $algorithm)")
    fi

    if [ $RUN_KATSURA_QQ -eq 1 ]; then
        TESTS_TO_RUN+=("katsura${KATSURA_N}-QQ-dp-${algorithm}:$(create_test_file katsura_qq_dp $algorithm)")
    fi

    if [ $RUN_KATSURA_ZZ -eq 1 ]; then
        TESTS_TO_RUN+=("katsura${KATSURA_N}-ZZ-dp-${algorithm}:$(create_test_file katsura_zz_dp $algorithm)")
    fi
done

if [ ${#TESTS_TO_RUN[@]} -eq 0 ]; then
    echo -e "${RED}No tests selected to run!${NC}"
    echo "Use --help to see available options"
    exit 1
fi

echo -e "${YELLOW}Tests to run: ${#TESTS_TO_RUN[@]}${NC}"
echo ""

# Run all benchmarks for each Singular version
for idx in "${!SINGULAR_EXECS[@]}"; do
    version_name="${SINGULAR_NAMES[$idx]}"
    executable="${SINGULAR_EXECS[$idx]}"
    ld_path="${SINGULAR_LD_PATHS[$idx]}"
    sing_path="${SINGULAR_PATHS[$idx]}"
    
    echo -e "${YELLOW}=== Running benchmarks for $version_name ===${NC}"
    echo ""

    for test_spec in "${TESTS_TO_RUN[@]}"; do
        IFS=':' read -r test_name test_file <<< "$test_spec"
        run_benchmark "$version_name" "$ld_path" "$sing_path" "$executable" "$test_name" "$test_file"
    done
done

# Print summary
echo "================================"
echo -e "${YELLOW}SUMMARY OF ALL BENCHMARKS${NC}"
echo "================================"
echo ""

# Parse results and create comparison
declare -A version_times
declare -A version_stddev

while IFS='|' read -r version test avg stddev min max total; do
    if [ "$version" != "Version" ]; then
        version_times["$version|$test"]="$avg"
        version_stddev["$version|$test"]="$stddev"
    fi
done < "$TEMP_DIR/results.txt"

# Get unique test names
declare -A test_set
for key in "${!version_times[@]}"; do
    test=$(echo "$key" | cut -d'|' -f2)
    test_set["$test"]=1
done

# Print comparison table (text format)
if [ ${#SINGULAR_EXECS[@]} -eq 1 ]; then
    # Only one version - simple results table
    echo -e "${BLUE}Results for ${SINGULAR_NAMES[0]}:${NC}"
    printf "%-40s %-20s\n" "Test" "Time"
    printf "%-40s %-20s\n" "----" "----"
    
    for test in "${!test_set[@]}"; do
        key="${SINGULAR_NAMES[0]}|$test"
        time="${version_times[$key]}"
        stddev="${version_stddev[$key]}"
        if [ $NUM_RUNS -gt 1 ]; then
            printf "%-40s %-20s\n" "$test" "${time}s (±${stddev})"
        else
            printf "%-40s %-20s\n" "$test" "${time}s"
        fi
    done
else
    # Multiple versions - comparison table
    echo -e "${BLUE}Average Times Comparison:${NC}"
    
    # Dynamic header based on number of versions
    printf "%-40s" "Test"
    for name in "${SINGULAR_NAMES[@]}"; do
        printf " %-25s" "$name"
    done
    echo ""
    
    printf "%-40s" "----"
    for name in "${SINGULAR_NAMES[@]}"; do
        printf " %-25s" "--------"
    done
    echo ""
    
    for test in "${!test_set[@]}"; do
        printf "%-40s" "$test"
        
        declare -a test_times
        for name in "${SINGULAR_NAMES[@]}"; do
            key="$name|$test"
            time="${version_times[$key]}"
            stddev="${version_stddev[$key]}"
            if [ $NUM_RUNS -gt 1 ]; then
                printf " %-25s" "${time}s (±${stddev})"
            else
                printf " %-25s" "${time}s"
            fi
            test_times+=("$time")
        done
        
        # Find best time and show comparison if multiple versions
        if [ ${#SINGULAR_NAMES[@]} -eq 2 ]; then
            t1="${test_times[0]}"
            t2="${test_times[1]}"
            if [ -n "$t1" ] && [ -n "$t2" ]; then
                speedup=$(echo "scale=2; $t1 / $t2" | bc)
                if (( $(echo "$speedup > 1.05" | bc -l) )); then
                    echo -e " ${GREEN}(${SINGULAR_NAMES[1]} ${speedup}x faster)${NC}"
                elif (( $(echo "$speedup < 0.95" | bc -l) )); then
                    speedup=$(echo "scale=2; $t2 / $t1" | bc)
                    echo -e " ${RED}(${SINGULAR_NAMES[0]} ${speedup}x faster)${NC}"
                else
                    echo " (similar)"
                fi
            else
                echo ""
            fi
        else
            echo ""
        fi
        
        unset test_times
    done
fi

echo ""
echo "================================"

# Now create markdown format
MARKDOWN_FILE="$TEMP_DIR/results.md"

if [ ${#SINGULAR_EXECS[@]} -eq 1 ]; then
    # Only one version markdown
    cat > "$MARKDOWN_FILE" << MDEOF
## Results for ${SINGULAR_NAMES[0]}

| Test | Time |
|------|------|
MDEOF

    # Sort tests for consistent output
    readarray -t sorted_tests < <(printf '%s\n' "${!test_set[@]}" | sort)
    
    for test in "${sorted_tests[@]}"; do
        key="${SINGULAR_NAMES[0]}|$test"
        time="${version_times[$key]}"
        stddev="${version_stddev[$key]}"
        
        # Format numbers with leading zeros
        time_fmt=$(printf "%.4f" "$time")
        stddev_fmt=$(printf "%.4f" "$stddev")
        
        if [ $NUM_RUNS -gt 1 ]; then
            echo "| $test | ${time_fmt}s (±${stddev_fmt}) |" >> "$MARKDOWN_FILE"
        else
            echo "| $test | ${time_fmt}s |" >> "$MARKDOWN_FILE"
        fi
    done
else
    # Comparison markdown
    cat > "$MARKDOWN_FILE" << MDEOF
## Average Times Comparison

| Test |$(for name in "${SINGULAR_NAMES[@]}"; do echo -n " $name |"; done)
|------|$(for name in "${SINGULAR_NAMES[@]}"; do echo -n "----------|"; done)
MDEOF

    # Sort tests for consistent output
    readarray -t sorted_tests < <(printf '%s\n' "${!test_set[@]}" | sort)
    
    for test in "${sorted_tests[@]}"; do
        echo -n "| $test |" >> "$MARKDOWN_FILE"
        
        for name in "${SINGULAR_NAMES[@]}"; do
            key="$name|$test"
            time="${version_times[$key]}"
            stddev="${version_stddev[$key]}"
            
            # Format numbers with leading zeros
            time_fmt=$(printf "%.4f" "$time")
            stddev_fmt=$(printf "%.4f" "$stddev")
            
            if [ $NUM_RUNS -gt 1 ]; then
                echo -n " ${time_fmt}s (±${stddev_fmt}) |" >> "$MARKDOWN_FILE"
            else
                echo -n " ${time_fmt}s |" >> "$MARKDOWN_FILE"
            fi
        done
        echo "" >> "$MARKDOWN_FILE"
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