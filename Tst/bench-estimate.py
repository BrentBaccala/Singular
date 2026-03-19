#!/usr/bin/env python3
"""Estimate wall clock time for a benchmark run at a given scale factor.

Uses measured timing data to predict how long a run will take. For looped
tests (A/B/D), scales the iteration count and estimates from per-iteration
wall time. For Class C tests, uses the measured time as-is.

Usage:
    bench-estimate.py [--scale FACTOR] FILE [FILE ...]
"""

import csv
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description='Estimate benchmark run time')
    parser.add_argument('files', nargs='+', help='CSV benchmark files')
    parser.add_argument('--scale', type=float, default=1.0,
                        help='Scale factor for looped iteration counts (default: 1.0)')
    args = parser.parse_args()

    total_s = 0
    tests = []

    for fname in args.files:
        with open(fname) as f:
            for row in csv.DictReader(f):
                cls = row['class']
                if cls in ('X', 'FAIL'):
                    continue

                iters = int(row['iterations'])
                warmup_wall_us = int(row['warmup_wall_us'])
                wall_us = int(row['wall_us'])
                ext_wall_ns = int(row['ext_wall_ns'])

                if cls == 'C':
                    # Single run, no scaling
                    est_s = ext_wall_ns / 1e9
                else:
                    # Compute scaled iterations (same math as bench-suite.sh)
                    scaled_iters = max(1, round(iters * args.scale))

                    # Per-iteration wall time
                    if iters > 0 and wall_us > 0:
                        per_iter_us = wall_us / iters
                    else:
                        per_iter_us = 0

                    # Overhead = ext_wall - (warmup_wall + loop_wall), in seconds
                    singular_total_us = warmup_wall_us + wall_us
                    overhead_s = max(0, ext_wall_ns / 1e9 - singular_total_us / 1e6)

                    # Estimate: warmup + scaled loop + overhead
                    est_s = warmup_wall_us / 1e6 + scaled_iters * per_iter_us / 1e6 + overhead_s

                tests.append({
                    'test': row['test_name'],
                    'class': cls,
                    'iters': iters,
                    'scaled_iters': scaled_iters if cls != 'C' else 1,
                    'est_s': est_s,
                })
                total_s += est_s

    # Print by descending estimated time
    print(f"Scale factor: {args.scale}")
    print()
    print(f"{'Test':<50} {'Class':>5} {'Iters':>7} {'Est(s)':>8}")
    print("-" * 73)
    for t in sorted(tests, key=lambda x: -x['est_s'])[:30]:
        print(f"{t['test']:<50} {t['class']:>5} {t['scaled_iters']:>7} {t['est_s']:>8.1f}")

    if len(tests) > 30:
        remaining = sum(t['est_s'] for t in sorted(tests, key=lambda x: -x['est_s'])[30:])
        print(f"{'... and ' + str(len(tests) - 30) + ' more':<50} {'':>5} {'':>7} {remaining:>8.1f}")

    print("-" * 73)
    print(f"{'TOTAL (' + str(len(tests)) + ' tests)':<50} {'':>5} {'':>7} {total_s:>8.1f}")
    print(f"{'':50} {'':>5} {'':>7} {total_s/60:>7.1f}m")
    print(f"{'':50} {'':>5} {'':>7} {total_s/3600:>7.2f}h")


if __name__ == '__main__':
    main()
