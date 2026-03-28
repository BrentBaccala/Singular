#!/bin/bash
# bench-launch.sh - Launch a benchmark run in the background via nohup.
#
# Write this script to avoid SSH/nohup quoting problems. scp it to the
# remote machine and run it with: ssh host 'bash ~/bench-launch.sh'
#
# Edit the variables below to configure the run, then scp and execute.

cd ~/Singular-benchmark/Tst

nohup ./bench-suite.sh \
    -n 5 --scale 0.25 --taskset "-c 11" --numactl "--membind=1" \
    -o ~/bench-spielwiese-full.csv \
    "spielwiese=~/Singular/install/bin/Singular --random=42" \
    -- Short/ok_s.lst Buch/buch.lst Long/ok_l.lst Old/universal.lst \
       ~/test-cases/singular/staging.lst \
    > ~/bench-spielwiese-full.log 2>&1 &

echo "PID: $!"
