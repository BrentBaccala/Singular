#!/bin/bash
# Run bug_532_s repeatedly until it hangs, then wait for GDB attach.
# Usage: bash hang-test.sh [binary] [count]

SINGULAR=${1:-~/Singular-parallel-bba/install/bin/Singular}
COUNT=${2:-100}

cd ~/Singular-gb-testsuite/Tst

for i in $(seq 1 $COUNT); do
  echo -n "Run $i... "
  SINGULAR_THREADS=4 "$SINGULAR" -t Short/bug_532_s.tst > /dev/null 2>&1 &
  PID=$!
  # Wait up to 10 seconds
  for j in $(seq 1 100); do
    if ! kill -0 $PID 2>/dev/null; then break; fi
    sleep 0.1
  done
  if kill -0 $PID 2>/dev/null; then
    echo "HUNG on run $i (PID $PID) — attach GDB:"
    echo "  gdb -p $PID -batch -ex 'thread apply all bt'"
    echo "Press Enter to kill it and continue, or Ctrl-C to stop."
    read
    kill -9 $PID 2>/dev/null
    wait $PID 2>/dev/null
  else
    wait $PID
    echo "ok"
  fi
done
echo "Completed $COUNT runs without hanging."
