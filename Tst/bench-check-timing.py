#!/usr/bin/env python3
"""Check that warmup + timed loop wall time ≈ external wall clock.

Reports tests where the discrepancy exceeds a threshold, which may
indicate timer corruption, timeouts, or other issues.

Usage:
    bench-check-timing.py [--threshold PERCENT] FILE [FILE ...]
"""

import csv
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Check benchmark timing consistency')
    parser.add_argument('files', nargs='+', help='CSV benchmark files')
    parser.add_argument('--threshold', type=float, default=10.0,
                        help='Report if discrepancy exceeds this percent (default: 10)')
    args = parser.parse_args()

    problems = []
    total = 0

    for fname in args.files:
        with open(fname) as f:
            for row in csv.DictReader(f):
                if row['class'] in ('X', 'FAIL'):
                    continue

                total += 1
                warmup_us = int(row['warmup_cpu_us'])
                warmup_wall_us = int(row['warmup_wall_us'])
                cpu_us = int(row['cpu_us'])
                wall_us = int(row['wall_us'])
                ext_wall_ns = int(row['ext_wall_ns'])

                # Singular wall times for warmup + timed loop (microseconds)
                singular_total_us = warmup_wall_us + wall_us

                # External wall clock (convert ns to us)
                ext_us = ext_wall_ns / 1000

                if ext_us == 0:
                    continue

                discrepancy_us = ext_us - singular_total_us
                discrepancy_pct = discrepancy_us / ext_us * 100

                if abs(discrepancy_pct) > args.threshold:
                    problems.append({
                        'file': fname,
                        'test': row['test_name'],
                        'build': row['build'],
                        'run': row['run'],
                        'warmup_wall_s': warmup_wall_us / 1e6,
                        'loop_wall_s': wall_us / 1e6,
                        'singular_total_s': singular_total_us / 1e6,
                        'ext_s': ext_us / 1e6,
                        'discrepancy_s': discrepancy_us / 1e6,
                        'discrepancy_pct': discrepancy_pct,
                    })

    if problems:
        print(f"{'Test':<45} {'Warmup':>7} {'Loop':>7} {'Sing':>7} {'Ext':>7} {'Diff':>7} {'Diff%':>6}")
        print("-" * 90)
        for p in sorted(problems, key=lambda x: -abs(x['discrepancy_pct'])):
            print(f"{p['test']:<45} {p['warmup_wall_s']:>7.2f} {p['loop_wall_s']:>7.2f} "
                  f"{p['singular_total_s']:>7.2f} {p['ext_s']:>7.2f} {p['discrepancy_s']:>7.2f} {p['discrepancy_pct']:>5.1f}%")
        print()

    ok = total - len(problems)
    print(f"{ok}/{total} tests within {args.threshold}% threshold, {len(problems)} outside")

if __name__ == '__main__':
    main()
