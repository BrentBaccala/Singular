#!/bin/bash
# bench-launch.sh - Launch a benchmark run in the background via nohup.
#
# Write this script to avoid SSH/nohup quoting problems. scp it to the
# remote machine and run it with: ssh host 'bash ~/bench-launch.sh'
#
# Edit the variables below to configure the run, then scp and execute.

cd ~/Singular-benchmark/Tst

nohup ./bench-suite.sh \
    -n 5 --taskset "-c 11" --numactl "--membind=1" \
    -o ~/bench-next-opt-staging.csv \
    "next-opt=~/Singular-next-opt/install/bin/Singular --random=42" \
    -- ~/test-cases/singular/staging.lst \
    > ~/bench-next-opt-staging.log 2>&1 &

echo "PID: $!"
