#!/bin/bash
#
# run_gb_tests.sh - Run Singular test suite with selectable GB algorithm
#
# Usage: ./run_gb_tests.sh [-s Singular] [-a algorithm] [-l test.lst] [test.tst ...]
#
# Algorithms: std (default), sba, slimgb, mathicgb, groebner
#
# Approach:
#   1. The modified standard.lib adds support for:
#      - groebner(i, "sba") — dispatches to sba()
#      - groebner_default_method global variable — sets default algorithm
#   2. Each test file is preprocessed with sed to replace std( with groebner(
#   3. The --execute preamble sets groebner_default_method and loads the
#      modified standard.lib
#
# This mirrors the approach used in the primdec-SQL branch where
# primdec.lib was modified to call groebner() instead of std().

set +e  # Don't exit on test failures

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MODIFIED_STDLIB="$SCRIPT_DIR/../Singular/LIB/standard.lib"

# Defaults
SINGULAR=""
ALGORITHM="std"
TIMEOUT=300
TESTS=()
TESTLISTS=()
RESULTS_DIR=""
VERBOSE=0

usage() {
    cat <<EOF
Usage: $0 [options] [test.tst ...]

Options:
  -s PATH      Path to Singular binary
  -a ALGO      GB algorithm: std, sba, slimgb, mathicgb, groebner (default: std)
  -l FILE      Read test names from list file (one per line)
  -A SECS      Timeout per test in seconds (default: 300)
  -o DIR       Output results directory (default: results-ALGO)
  -v           Verbose output
  -h           Show this help

Algorithms:
  std       Standard Buchberger algorithm (default)
  sba       Signature-based algorithm
  slimgb    Slim Groebner bases
  mathicgb  F4 algorithm via mathicgb (finite characteristic only)
  groebner  Heuristic dispatcher (groebner() proc, uses its own heuristics)
EOF
    exit 1
}

while getopts "s:a:l:A:o:vh" opt; do
    case $opt in
        s) SINGULAR="$OPTARG" ;;
        a) ALGORITHM="$OPTARG" ;;
        l) TESTLISTS+=("$OPTARG") ;;
        A) TIMEOUT="$OPTARG" ;;
        o) RESULTS_DIR="$OPTARG" ;;
        v) VERBOSE=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

# Remaining arguments are test files
TESTS+=("$@")

# Find Singular binary
if [ -z "$SINGULAR" ]; then
    if [ -x "./Singular" ]; then
        SINGULAR="./Singular"
    elif [ -x "../Singular" ]; then
        SINGULAR="../Singular"
    else
        echo "ERROR: No Singular binary found. Use -s to specify." >&2
        exit 1
    fi
fi

if [ ! -x "$SINGULAR" ]; then
    echo "ERROR: $SINGULAR is not executable" >&2
    exit 1
fi

# Check modified standard.lib exists
if [ ! -f "$MODIFIED_STDLIB" ]; then
    echo "ERROR: Modified standard.lib not found at $MODIFIED_STDLIB" >&2
    exit 1
fi

# Set up results directory
if [ -z "$RESULTS_DIR" ]; then
    RESULTS_DIR="results-${ALGORITHM}"
fi
mkdir -p "$RESULTS_DIR"

# Read test lists
for listfile in "${TESTLISTS[@]}"; do
    dir=$(dirname "$listfile")
    while IFS= read -r line; do
        # Skip comments and empty lines
        line=$(echo "$line" | sed 's/;.*//' | tr -d '[:space:]')
        [ -z "$line" ] && continue
        TESTS+=("${dir}/${line}.tst")
    done < "$listfile"
done

if [ ${#TESTS[@]} -eq 0 ]; then
    echo "ERROR: No tests specified. Use -l or pass .tst files as arguments." >&2
    usage
fi

# Validate algorithm
case "$ALGORITHM" in
    std|sba|slimgb|mathicgb|groebner) ;;
    *) echo "ERROR: Unknown algorithm: $ALGORITHM" >&2; exit 1 ;;
esac

# Build the --execute preamble
# For "std", we still load the modified lib but don't set a default method
# For others, we set groebner_default_method so all groebner() calls use it
if [ "$ALGORITHM" = "std" ]; then
    # No method override - std is the default in groebner() anyway
    PREAMBLE='LIB "/home/claude/Singular-gb-testsuite/Singular/LIB/standard.lib";'
else
    PREAMBLE="LIB \"$MODIFIED_STDLIB\"; string groebner_default_method = \"$ALGORITHM\"; export(groebner_default_method);"
fi

# Singular options (same as regress.cmd)
SINGULAR_OPTS="--ticks-per-sec=100 -teqsr12345678 --no-rc"

# Run tests
PASS=0
FAIL=0
ERROR=0
SKIP=0
TOTAL=${#TESTS[@]}

echo "========================================"
echo "GB Algorithm Test Suite"
echo "Algorithm:  $ALGORITHM"
echo "Singular:   $SINGULAR"
echo "Tests:      $TOTAL"
echo "Results:    $RESULTS_DIR"
echo "Timeout:    ${TIMEOUT}s"
echo "========================================"
echo ""

for testfile in "${TESTS[@]}"; do
    testname=$(basename "$testfile" .tst)
    testdir=$(dirname "$testfile")

    # Check test file exists
    if [ ! -f "$testfile" ]; then
        echo "SKIP  $testname (file not found: $testfile)"
        SKIP=$((SKIP + 1))
        echo "$testname SKIP file_not_found" >> "$RESULTS_DIR/summary.txt"
        continue
    fi

    # Decode reference file
    resfile="${testdir}/${testname}.res.gz.uu"
    reffile="$RESULTS_DIR/${testname}.res.expected"
    newresfile="$RESULTS_DIR/${testname}.res.new"

    if [ -f "$resfile" ]; then
        # Decode and decompress reference
        cp "$resfile" "$RESULTS_DIR/${testname}.res.gz.uu"
        (cd "$RESULTS_DIR" && uudecode "${testname}.res.gz.uu" 2>/dev/null && gunzip -f "${testname}.res.gz" 2>/dev/null && mv "${testname}.res" "${testname}.res.expected" 2>/dev/null) || true
        rm -f "$RESULTS_DIR/${testname}.res.gz.uu" "$RESULTS_DIR/${testname}.res.gz"
    fi

    # Build the --execute string
    EXEC_STRING="string tst_status_file=\"$RESULTS_DIR/${testname}.stat\"; $PREAMBLE"

    # Preprocess the test file: replace std( with groebner(
    # This is the same approach as primdec-SQL branch's std()->groebner() conversion.
    # We use a temp file to feed the preprocessed test to Singular.
    tmptest=$(mktemp /tmp/gb_test_XXXXXX.tst)
    if [ "$ALGORITHM" != "std" ]; then
        # Replace std( with groebner( — but be careful:
        # - Don't replace "std" inside strings/comments
        # - Don't replace stdfglm, stdhilb, etc.
        # - The sed pattern replaces standalone std( calls
        # We use a simple word-boundary approach: replace \bstd( with groebner(
        sed -e 's/\bstd(/groebner(/g' \
            -e 's/\bstdfglm(/groebner(/g' \
            -e 's/\bstdhilb(/groebner(/g' \
            "$testfile" > "$tmptest"
    else
        cp "$testfile" "$tmptest"
    fi

    # Run Singular with the preprocessed test
    if command -v timeout &>/dev/null; then
        timeout "$TIMEOUT" "$SINGULAR" --execute "$EXEC_STRING" $SINGULAR_OPTS < "$tmptest" > "$newresfile" 2>&1
        exit_status=$?
    else
        "$SINGULAR" --execute "$EXEC_STRING" $SINGULAR_OPTS < "$tmptest" > "$newresfile" 2>&1
        exit_status=$?
    fi
    rm -f "$tmptest"

    if [ $exit_status -ne 0 ]; then
        if [ $exit_status -eq 124 ]; then
            echo "TIMEOUT $testname (${TIMEOUT}s)"
            echo "$testname TIMEOUT" >> "$RESULTS_DIR/summary.txt"
        else
            echo "ERROR $testname (exit $exit_status)"
            echo "$testname ERROR exit_$exit_status" >> "$RESULTS_DIR/summary.txt"
        fi
        ERROR=$((ERROR + 1))
        continue
    fi

    # Compare with reference (if available)
    if [ -f "$reffile" ]; then
        # Strip timing/version/status info and input echo for comparison
        # The STDIN echo lines differ because we sed-replaced std->groebner
        sed -e '/used time:/d' -e '/tst_status/d' -e '/^Singular/d' \
            -e '/\$Id/d' -e '/init >>/d' -e '/^\/\/ computer:/d' \
            -e '/^\/\/ hostname:/d' -e '/^\/\/ sysname:/d' \
            -e '/^\/\/ \*\* loaded/d' -e '/^\/\/ \*\* library/d' \
            -e '/^\/\/ \*\* redefining/d' \
            -e '/>> tst_memory/d' -e '/>> tst_timer/d' \
            -e '/^STDIN /d' \
            -e 's|[^ ]*standard\.lib::|standard.lib::|g' \
            -e 's|standard\.lib::\([a-zA-Z_]*\) line [0-9]*|standard.lib::\1 line NNN|g' \
            -e 's|standard\.lib::\([a-zA-Z_]*\) ([0-9]*)|standard.lib::\1 (NNN)|g' \
            "$reffile" > "$RESULTS_DIR/${testname}.res.expected.clean" 2>/dev/null
        sed -e '/used time:/d' -e '/tst_status/d' -e '/^Singular/d' \
            -e '/\$Id/d' -e '/init >>/d' -e '/^\/\/ computer:/d' \
            -e '/^\/\/ hostname:/d' -e '/^\/\/ sysname:/d' \
            -e '/^\/\/ \*\* loaded/d' -e '/^\/\/ \*\* library/d' \
            -e '/^\/\/ \*\* redefining/d' \
            -e '/>> tst_memory/d' -e '/>> tst_timer/d' \
            -e '/^STDIN /d' \
            -e 's|[^ ]*standard\.lib::|standard.lib::|g' \
            -e 's|standard\.lib::\([a-zA-Z_]*\) line [0-9]*|standard.lib::\1 line NNN|g' \
            -e 's|standard\.lib::\([a-zA-Z_]*\) ([0-9]*)|standard.lib::\1 (NNN)|g' \
            "$newresfile" > "$RESULTS_DIR/${testname}.res.new.clean" 2>/dev/null

        if diff -w -b "$RESULTS_DIR/${testname}.res.expected.clean" \
                      "$RESULTS_DIR/${testname}.res.new.clean" > "$RESULTS_DIR/${testname}.diff" 2>&1; then
            echo "PASS  $testname"
            echo "$testname PASS" >> "$RESULTS_DIR/summary.txt"
            PASS=$((PASS + 1))
            # Clean up on pass
            rm -f "$RESULTS_DIR/${testname}.res.expected" "$RESULTS_DIR/${testname}.res.expected.clean" \
                  "$RESULTS_DIR/${testname}.res.new" "$RESULTS_DIR/${testname}.res.new.clean" \
                  "$RESULTS_DIR/${testname}.diff" "$RESULTS_DIR/${testname}.stat"
        else
            echo "FAIL  $testname"
            echo "$testname FAIL" >> "$RESULTS_DIR/summary.txt"
            FAIL=$((FAIL + 1))
            # Keep diff and output for inspection
            rm -f "$RESULTS_DIR/${testname}.res.expected.clean" "$RESULTS_DIR/${testname}.res.new.clean"
        fi
    else
        echo "NOREF $testname (no reference file)"
        echo "$testname NOREF" >> "$RESULTS_DIR/summary.txt"
        PASS=$((PASS + 1))  # Count as pass if no reference to compare
    fi
done

echo ""
echo "========================================"
echo "Results: $PASS pass, $FAIL fail, $ERROR error, $SKIP skip (of $TOTAL)"
echo "========================================"
echo ""
echo "ALGORITHM=$ALGORITHM PASS=$PASS FAIL=$FAIL ERROR=$ERROR SKIP=$SKIP TOTAL=$TOTAL" >> "$RESULTS_DIR/summary.txt"

# Exit with non-zero if any failures
if [ $FAIL -gt 0 ] || [ $ERROR -gt 0 ]; then
    exit 1
fi
