#!/usr/bin/env python3
"""Estimate wall clock time for a benchmark run at a given scale factor.

Accepts a mix of benchmark CSV files and classification CSV files.
Uses classification files for iteration counts (what will run next time)
and benchmark files for per-iteration timing (measured data).

Automatically distinguishes file types by their CSV headers.

Usage:
    bench-estimate.py [--scale FACTOR] FILE [FILE ...]

Examples:
    # Estimate from benchmark data alone (uses measured iteration counts)
    bench-estimate.py bench-short-run1.csv bench-buch-run1.csv

    # Estimate with recalibrated iteration counts from classification files
    bench-estimate.py bench-short-run1.csv Short/ok_s_classification.csv

    # Estimate at half scale
    bench-estimate.py --scale 0.5 bench-*.csv *_classification.csv
"""

import csv
import sys
import os
import argparse


def is_benchmark_file(fname):
    """Check if a CSV file is a benchmark results file (vs classification)."""
    with open(fname) as f:
        header = f.readline().strip()
    return header.startswith('build,')


def load_benchmark_data(fname):
    """Load per-test timing from a benchmark CSV. Returns dict keyed by test basename."""
    data = {}
    with open(fname) as f:
        for row in csv.DictReader(f):
            key = os.path.basename(row['test_name'])
            cls = row['class']
            if cls in ('X', 'FAIL'):
                continue

            iters = int(row['iterations'])
            warmup_wall_us = int(row['warmup_wall_us'])
            wall_us = int(row['wall_us'])
            ext_wall_ns = int(row['ext_wall_ns'])

            # Per-iteration wall time
            if iters > 0 and wall_us > 0:
                per_iter_us = wall_us / iters
            else:
                per_iter_us = 0

            # Overhead = ext_wall - (warmup_wall + loop_wall)
            singular_total_us = warmup_wall_us + wall_us
            overhead_s = max(0, ext_wall_ns / 1e9 - singular_total_us / 1e6)

            data[key] = {
                'class': cls,
                'measured_iters': iters,
                'per_iter_us': per_iter_us,
                'warmup_wall_s': warmup_wall_us / 1e6,
                'overhead_s': overhead_s,
                'ext_wall_s': ext_wall_ns / 1e9,
            }
    return data


def load_classification(fname):
    """Load classification CSV. Returns dict keyed by test basename."""
    classifications = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            key = os.path.basename(parts[0])
            cls = parts[1]
            iters = int(parts[2])
            classifications[key] = {
                'class': cls,
                'iters': iters,
            }
    return classifications


def main():
    parser = argparse.ArgumentParser(description='Estimate benchmark run time')
    parser.add_argument('files', nargs='+', help='Benchmark CSV and/or classification CSV files')
    parser.add_argument('--scale', type=float, default=1.0,
                        help='Scale factor for looped iteration counts (default: 1.0)')
    args = parser.parse_args()

    # Separate file types
    bench_data = {}  # test basename -> timing data
    classifications = {}  # test basename -> class + iters

    for fname in args.files:
        if is_benchmark_file(fname):
            bench_data.update(load_benchmark_data(fname))
        else:
            classifications.update(load_classification(fname))

    # For tests in benchmark data but not in classification, use measured iters
    # For tests in classification, use classification iters
    # Only estimate tests that have benchmark data (we need timing)

    total_s = 0
    tests = []

    for key, timing in bench_data.items():
        cls = timing['class']

        if cls in ('X', 'FAIL'):
            continue

        # Get iteration count: prefer classification, fall back to measured
        if key in classifications:
            cal = classifications[key]
            cal_cls = cal['class']
            iters = cal['iters']
            if cal_cls in ('X', 'FAIL'):
                continue
            if cal_cls == 'C':
                cls = 'C'  # classification overrides
        else:
            iters = timing['measured_iters']

        if cls == 'C':
            est_s = timing['ext_wall_s']
            scaled_iters = 1
        else:
            scaled_iters = max(1, round(iters * args.scale))
            est_s = (timing['warmup_wall_s']
                     + scaled_iters * timing['per_iter_us'] / 1e6
                     + timing['overhead_s'])

        tests.append({
            'test': key,
            'class': cls,
            'iters': iters,
            'scaled_iters': scaled_iters,
            'est_s': est_s,
        })
        total_s += est_s

    # Warn about classification entries with no benchmark data
    missing = set(classifications.keys()) - set(bench_data.keys())
    missing = {k for k in missing if classifications[k]['class'] not in ('X', 'FAIL')}
    if missing:
        print(f"Warning: {len(missing)} classified tests have no benchmark data:", file=sys.stderr)
        for k in sorted(missing)[:10]:
            print(f"  {k}", file=sys.stderr)
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more", file=sys.stderr)

    # Print by descending estimated time
    print(f"Scale factor: {args.scale}")
    if classifications:
        print(f"Using calibrated iteration counts from {sum(1 for f in args.files if not is_benchmark_file(f))} classification file(s)")
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
