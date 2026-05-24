#!/bin/bash
# run-tests.sh - Run regress.cmd in the background via nohup.
#
# Avoids SSH/nohup quoting problems and Bash tool timeouts.
# Launch with: bash Tst/run-tests.sh
# Monitor with: tail -f ~/project/test-logs/<logfile>
#
# Edit the variables below to configure, or pass arguments:
#   bash run-tests.sh [SINGULAR_BINARY] [LIST_FILES...]
#
# Defaults use the install binary from this worktree and
# the comprehensive test suite (all five lists).

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKTREE="$(cd "$SCRIPT_DIR/.." && pwd)"
BRANCH=$(cd "$WORKTREE" && git branch --show-current 2>/dev/null || basename "$WORKTREE")
DATE=$(date +%d%b%Y-%H%M)
LOGDIR=~/project/test-logs
LOGFILE="${LOGDIR}/${BRANCH}-${DATE}.log"

# Default binary: install/bin/Singular from this worktree
SINGULAR="${1:-${WORKTREE}/install/bin/Singular}"
shift 2>/dev/null

# Default lists: comprehensive suite
if [ $# -eq 0 ]; then
    LISTS="Old/universal.lst Buch/buch.lst Plural/short.lst Short/ok_s.lst Long/ok_l.lst"
else
    LISTS="$*"
fi

# Verify binary exists
if [ ! -x "$SINGULAR" ]; then
    echo "Error: Singular binary not found: $SINGULAR" >&2
    exit 1
fi

# Create log directory
mkdir -p "$LOGDIR"

# Per-test execution timeout (seconds), passed to regress.cmd via -A.
# Without it regress.cmd's default is 0 = no timeout, so a single hung
# Singular (e.g. heap corruption under SINGULAR_THREADS>1 that lands the
# process in its SIGSEGV-handler restart loop) blocks the whole run
# indefinitely.  REGRESS_TIMEOUT=0 disables the cap (legacy behaviour).
REGRESS_TIMEOUT="${REGRESS_TIMEOUT:-300}"
TIMEOUT_ARGS=""
if [ "$REGRESS_TIMEOUT" -gt 0 ] 2>/dev/null; then
    TIMEOUT_ARGS="-A $REGRESS_TIMEOUT"
fi

echo "Branch:    $BRANCH"
echo "Binary:    $SINGULAR"
echo "Lists:     $LISTS"
echo "Log:       $LOGFILE"
echo "Per-test timeout: ${REGRESS_TIMEOUT}s"
echo ""

cd "$SCRIPT_DIR"

nohup ./regress.cmd $TIMEOUT_ARGS -s "$SINGULAR" $LISTS > "$LOGFILE" 2>&1 &
PID=$!

echo "PID: $PID"
echo "Monitor: tail -f $LOGFILE"
echo "Results: grep -E 'Summary:|FAIL|!!!' $LOGFILE"
