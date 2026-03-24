#!/usr/bin/env python3
"""Generate strip plots of benchmark timings, sorted by mean time.

Produces a multi-page PDF with one horizontal line per test, dots for
each run, colored by build. Tests are sorted shortest-first, split
across pages (~45 per page) with appropriate x-axis scaling.

Usage:
    bench-stripplot.py [-o OUTPUT] FILE [FILE ...]
    bench-stripplot.py -o wall.pdf --metric wall_us FILE
    bench-stripplot.py -o cpu.pdf --metric cpu_us FILE

Examples:
    # Generate both wall and CPU time plots
    bench-stripplot.py -o timings-wall.pdf --metric wall_us data.csv
    bench-stripplot.py -o timings-cpu.pdf --metric cpu_us data.csv

    # Multiple CSV files (e.g., different builds)
    bench-stripplot.py -o compare.pdf build1.csv build2.csv
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import csv
import argparse
import math
from collections import defaultdict


def format_time(us, _pos=None):
    """Format microseconds into a human-readable string."""
    if us >= 1e6:
        return f'{us/1e6:.1f}s'
    elif us >= 1e3:
        return f'{us/1e3:.1f}ms'
    else:
        return f'{us:.0f}µs'


def main():
    parser = argparse.ArgumentParser(
        description='Strip plot of benchmark timings')
    parser.add_argument('files', nargs='+', help='CSV benchmark files')
    parser.add_argument('-o', '--output', default='bench-timings.pdf',
                        help='Output PDF file (default: bench-timings.pdf)')
    parser.add_argument('--metric', choices=['wall_us', 'cpu_us'],
                        default='wall_us',
                        help='Timing metric (default: wall_us)')
    parser.add_argument('--per-page', type=int, default=45,
                        help='Tests per page (default: 45)')
    parser.add_argument('--title', type=str, default=None,
                        help='Plot title (default: auto from metric)')
    args = parser.parse_args()

    metric = args.metric
    metric_label = 'Wall Time' if metric == 'wall_us' else 'CPU Time'
    title = args.title or f'Benchmark {metric_label} by Test (sorted by mean)'

    # Read all data
    # Key: (build, test_name) -> list of metric values
    data = defaultdict(list)
    for fname in args.files:
        with open(fname) as f:
            for row in csv.DictReader(f):
                if row['class'] in ('X', 'FAIL'):
                    continue
                build = row['build']
                test = row['test_name']
                val = int(row[metric])
                data[(build, test)].append(val)

    # Get unique builds and tests
    builds = sorted(set(b for b, t in data.keys()))
    tests = sorted(set(t for b, t in data.keys()))

    # Compute mean per test (across all builds) for sorting
    test_means = {}
    for test in tests:
        all_vals = []
        for build in builds:
            all_vals.extend(data.get((build, test), []))
        test_means[test] = sum(all_vals) / len(all_vals) if all_vals else 0

    # Sort tests by mean time (shortest first)
    tests_sorted = sorted(tests, key=lambda t: test_means[t])

    # Build colors
    build_colors = {b: f'C{i}' for i, b in enumerate(builds)}

    # Split into pages, balanced so last page is similar size to others
    max_per_page = args.per_page
    n_tests_total = len(tests_sorted)
    n_pages = max(1, math.ceil(n_tests_total / max_per_page))
    base = n_tests_total // n_pages
    extra = n_tests_total % n_pages  # first 'extra' pages get base+1

    with PdfPages(args.output) as pdf:
        offset = 0
        for page in range(n_pages):
            count = base + (1 if page < extra else 0)
            page_tests = tests_sorted[offset:offset + count]
            offset += count
            n_tests = len(page_tests)

            fig_height = max(6, n_tests * 0.28 + 1.5)
            fig, ax = plt.subplots(figsize=(14, fig_height))

            # Draw thin gray horizontal lines for each test
            for yi in range(n_tests):
                ax.axhline(y=yi, color='gray', linewidth=0.3, zorder=1)

            # Offset builds vertically so dots sit above/below the guide line
            # First build: shifted up (bottom touches line)
            # Second build: shifted down (top touches line)
            n_builds = len(builds)
            if n_builds == 1:
                offsets = {builds[0]: 0}
            else:
                offsets = {builds[i]: -0.12 + i * 0.24 / (n_builds - 1)
                           for i in range(n_builds)}

            for yi, test in enumerate(page_tests):
                for build in builds:
                    vals = data.get((build, test), [])
                    if not vals:
                        continue
                    y_pos = [yi + offsets[build]] * len(vals)
                    ax.scatter(vals, y_pos,
                               s=16, color=build_colors[build],
                               label=build if yi == 0 else None,
                               alpha=0.8, zorder=5, edgecolors='none')

            ax.set_ylim(-0.7, n_tests - 1 + 0.7)
            ax.set_yticks(range(n_tests))
            # Strip common path prefixes for readability
            labels = []
            for t in page_tests:
                # Strip leading path components that are in list_file
                label = t.split('/')[-1] if '/' in t else t
                labels.append(label)
            ax.set_yticklabels(labels, fontsize=7, family='monospace')
            ax.invert_yaxis()  # shortest at top

            ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_time))
            ax.set_xlabel(metric_label, fontsize=10)
            ax.grid(axis='x', alpha=0.3)

            # x-axis: auto-scale per page, never show negative time
            ax.margins(x=0.05)
            ax.set_xlim(left=max(0, ax.get_xlim()[0]))

            page_title = f'{title}  (page {page+1}/{n_pages})'
            ax.set_title(page_title, fontsize=11, fontweight='bold')

            if len(builds) > 1:
                ax.legend(loc='upper right', markerscale=2, fontsize=8)

            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

    print(f'Wrote {n_pages} pages to {args.output}')
    print(f'{len(tests_sorted)} tests, {len(builds)} build(s): {", ".join(builds)}')


if __name__ == '__main__':
    main()
