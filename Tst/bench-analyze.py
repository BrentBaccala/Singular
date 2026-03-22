#!/usr/bin/env python3
"""Analyze benchmark results: mean, stddev, CV% per test.

Optionally flags tests that overlap with a given time window (e.g., a
concurrent CPU-heavy task) to check for L3 cache or scheduling interference.

Usage:
    bench-analyze.py FILE [FILE ...]
    bench-analyze.py --during START END FILE    # flag tests during time window
    bench-analyze.py --log-end TIMESTAMP FILE   # reconstruct absolute times

Time window analysis (--during):
    Requires --log-end to set the absolute end time of the benchmark log.
    Reconstructs when each test ran by cumulating ext_wall_ns backwards from
    the log end time, then splits statistics into "during" and "outside" the
    given window.

    TIMESTAMP format: "2026-03-21 03:18:22" (local time)
    START/END format:  "2026-03-20 18:54" (local time)

Examples:
    # Basic stats
    bench-analyze.py bench-10runs.csv

    # Check if task 209 (18:54-21:31) affected benchmarks
    bench-analyze.py --log-end "2026-03-21 03:18:22" \\
        --during "2026-03-20 18:54" "2026-03-20 21:31" \\
        bench-10runs.csv

    # Sort by CV% descending
    bench-analyze.py --sort cv bench-10runs.csv

    # Only show tests with CV > 5%
    bench-analyze.py --min-cv 5 bench-10runs.csv
"""

import csv
import sys
import argparse
import math
from datetime import datetime, timedelta
from collections import defaultdict


def parse_timestamp(s):
    """Parse a timestamp string."""
    for fmt in ('%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M'):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(f"Cannot parse timestamp: {s}")


def compute_stats(values):
    """Return (mean, stddev, cv_pct) for a list of numeric values."""
    n = len(values)
    if n == 0:
        return (0, 0, 0)
    mean = sum(values) / n
    if n == 1:
        return (mean, 0, 0)
    variance = sum((x - mean) ** 2 for x in values) / (n - 1)
    stddev = math.sqrt(variance)
    cv = (stddev / mean * 100) if mean > 0 else 0
    return (mean, stddev, cv)


def main():
    parser = argparse.ArgumentParser(
        description='Analyze benchmark results: mean, stddev, CV%')
    parser.add_argument('files', nargs='+', help='CSV benchmark files')
    parser.add_argument('--sort', choices=['name', 'cv', 'mean', 'stddev'],
                        default='name', help='Sort order (default: name)')
    parser.add_argument('--min-cv', type=float, default=0,
                        help='Only show tests with CV%% above this threshold')
    parser.add_argument('--min-runs', type=int, default=2,
                        help='Minimum runs needed for stats (default: 2)')
    parser.add_argument('--metric', choices=['cpu_us', 'wall_us', 'ext_wall_ns'],
                        default='wall_us',
                        help='Which timing metric to analyze (default: wall_us)')
    parser.add_argument('--log-end', type=str, default=None,
                        help='Log file end timestamp for absolute time reconstruction')
    parser.add_argument('--during', nargs=2, metavar=('START', 'END'),
                        help='Flag tests overlapping this time window')
    parser.add_argument('--summary', action='store_true',
                        help='Print only summary statistics')
    args = parser.parse_args()

    if args.during and not args.log_end:
        parser.error('--during requires --log-end')

    # Read all rows
    rows = []
    for fname in args.files:
        with open(fname) as f:
            for row in csv.DictReader(f):
                if row['class'] in ('X', 'FAIL'):
                    continue
                rows.append(row)

    if not rows:
        print("No data rows found.")
        return

    # Reconstruct absolute timestamps if needed
    row_times = {}
    if args.log_end:
        log_end = parse_timestamp(args.log_end)
        # Walk backwards from end, cumulating ext_wall_ns
        cumulative_ns = 0
        for i in range(len(rows) - 1, -1, -1):
            ext_ns = int(rows[i]['ext_wall_ns'])
            cumulative_ns += ext_ns
            test_start = log_end - timedelta(seconds=cumulative_ns / 1e9)
            test_end = test_start + timedelta(seconds=ext_ns / 1e9)
            row_times[i] = (test_start, test_end)

    # Parse the "during" window
    during_start = during_end = None
    if args.during:
        during_start = parse_timestamp(args.during[0])
        during_end = parse_timestamp(args.during[1])

    # Group by test name
    # Key: test_name, Value: list of (metric_value, is_during_window)
    tests = defaultdict(list)
    for i, row in enumerate(rows):
        metric = args.metric
        if metric == 'ext_wall_ns':
            val = int(row[metric]) / 1000  # convert to us for display
        else:
            val = int(row[metric])

        is_during = False
        if during_start and i in row_times:
            t_start, t_end = row_times[i]
            # Overlap check
            if t_start < during_end and t_end > during_start:
                is_during = True

        tests[row['test_name']].append((val, is_during))

    # Compute stats per test
    results = []
    for name, entries in tests.items():
        values = [v for v, _ in entries]
        if len(values) < args.min_runs:
            continue

        mean, stddev, cv = compute_stats(values)

        # Split during/outside if requested
        during_vals = [v for v, d in entries if d]
        outside_vals = [v for v, d in entries if not d]

        result = {
            'name': name,
            'n': len(values),
            'mean': mean,
            'stddev': stddev,
            'cv': cv,
        }

        if args.during:
            result['n_during'] = len(during_vals)
            result['n_outside'] = len(outside_vals)
            if during_vals and outside_vals:
                d_mean, d_std, d_cv = compute_stats(during_vals)
                o_mean, o_std, o_cv = compute_stats(outside_vals)
                result['during_mean'] = d_mean
                result['outside_mean'] = o_mean
                # Percent difference: (during - outside) / outside * 100
                if o_mean > 0:
                    result['pct_diff'] = (d_mean - o_mean) / o_mean * 100
                else:
                    result['pct_diff'] = 0
            else:
                result['during_mean'] = compute_stats(during_vals)[0] if during_vals else None
                result['outside_mean'] = compute_stats(outside_vals)[0] if outside_vals else None
                result['pct_diff'] = None

        results.append(result)

    if not results:
        print("No tests with enough runs.")
        return

    # Filter by min CV
    if args.min_cv > 0:
        results = [r for r in results if r['cv'] >= args.min_cv]

    # Sort
    sort_keys = {
        'name': lambda r: r['name'],
        'cv': lambda r: -r['cv'],
        'mean': lambda r: -r['mean'],
        'stddev': lambda r: -r['stddev'],
    }
    results.sort(key=sort_keys[args.sort])

    # Summary stats
    all_cvs = [r['cv'] for r in results]
    median_cv = sorted(all_cvs)[len(all_cvs) // 2] if all_cvs else 0
    mean_cv = sum(all_cvs) / len(all_cvs) if all_cvs else 0
    high_cv = sum(1 for c in all_cvs if c > 5)

    if not args.summary:
        # Print per-test table
        if args.during:
            print(f"{'Test':<42} {'N':>3} {'Mean(us)':>12} {'Stddev':>10} {'CV%':>6}"
                  f"  {'Dur':>3} {'Out':>3} {'Dur_mean':>12} {'Out_mean':>12} {'Diff%':>7}")
            print("-" * 120)
            for r in results:
                dur_str = f"{r['during_mean']:>12.1f}" if r.get('during_mean') is not None else f"{'—':>12}"
                out_str = f"{r['outside_mean']:>12.1f}" if r.get('outside_mean') is not None else f"{'—':>12}"
                diff_str = f"{r['pct_diff']:>6.2f}%" if r.get('pct_diff') is not None else f"{'—':>7}"
                print(f"{r['name']:<42} {r['n']:>3} {r['mean']:>12.1f} {r['stddev']:>10.1f} {r['cv']:>5.2f}%"
                      f"  {r.get('n_during', 0):>3} {r.get('n_outside', 0):>3} {dur_str} {out_str} {diff_str}")
        else:
            print(f"{'Test':<42} {'N':>3} {'Mean(us)':>12} {'Stddev':>10} {'CV%':>6}")
            print("-" * 78)
            for r in results:
                print(f"{r['name']:<42} {r['n']:>3} {r['mean']:>12.1f} {r['stddev']:>10.1f} {r['cv']:>5.2f}%")

        print()

    # Print summary
    print(f"Tests analyzed: {len(results)}")
    print(f"CV%: mean={mean_cv:.2f}%, median={median_cv:.2f}%, max={max(all_cvs):.2f}%")
    print(f"Tests with CV > 5%: {high_cv}/{len(results)}")
    print(f"Tests with CV > 2%: {sum(1 for c in all_cvs if c > 2)}/{len(results)}")
    print(f"Tests with CV > 1%: {sum(1 for c in all_cvs if c > 1)}/{len(results)}")

    if args.during:
        # Summary of during vs outside
        diffs = [r['pct_diff'] for r in results if r.get('pct_diff') is not None]
        if diffs:
            mean_diff = sum(diffs) / len(diffs)
            median_diff = sorted(diffs)[len(diffs) // 2]
            print(f"\nDuring-window impact ({args.during[0]} to {args.during[1]}):")
            print(f"  Tests with both during+outside data: {len(diffs)}")
            print(f"  Mean difference: {mean_diff:+.2f}%")
            print(f"  Median difference: {median_diff:+.2f}%")
            print(f"  Tests >2% slower during: {sum(1 for d in diffs if d > 2)}")
            print(f"  Tests >2% faster during: {sum(1 for d in diffs if d < -2)}")


if __name__ == '__main__':
    main()
