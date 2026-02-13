#!/usr/bin/env python3
"""
SSI Trace Viewer - Pretty print SSI trace files

This script reads SSI trace files and displays them in a human-readable format,
showing ring details and summarizing polynomials/ideals.

Usage:
    python3 ssi_trace_viewer.py <trace_file.ssi>
    python3 ssi_trace_viewer.py <trace_file.ssi> --verbose  # Show polynomial details
"""

import sys
import os

# Add current directory to path for importing ssi_reader
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ssi_reader import (
    SSIReader, SSIRing, SSIPoly, SSIIdeal, SSIModule, SSIMatrix,
    SSICommand, SSIVersion, SSIQuit, SSINone, SSIProc, SSIDef,
    SSIBlackbox, SSIAttributed, SSITerm
)
from typing import Any


class SSITraceViewer:
    """Pretty printer for SSI traces."""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.message_count = 0
        self.current_ring = None

    def format_ordering_block(self, block):
        """Format an ordering block tuple."""
        ord_type, block0, block1, weights = block

        ord_names = {
            0: "lp", 1: "dp", 2: "Dp", 3: "wp", 4: "Wp",
            5: "ws", 6: "Ws", 7: "M", 8: "a", 9: "aa",
            2: "C", 4: "c"  # Component orderings
        }

        ord_name = ord_names.get(ord_type, f"ord{ord_type}")

        if weights:
            return f"{ord_name}({block0}..{block1}, weights={weights})"
        else:
            return f"{ord_name}({block0}..{block1})"

    def format_ring(self, ring: SSIRing, indent=""):
        """Format ring details."""
        lines = []

        # Coefficient field
        if ring.characteristic == 0:
            if ring.extension_type == 'transcendental':
                cf = "Q(transcendental extension)"
            elif ring.extension_type == 'algebraic':
                cf = "Q(algebraic extension)"
            else:
                cf = "Q"
        elif ring.characteristic > 0:
            cf = f"Z/{ring.characteristic}"
        elif ring.characteristic == -3:
            cf = "Named coefficient field"
        else:
            cf = f"Special char={ring.characteristic}"

        lines.append(f"{indent}Ring: {cf}")

        # Variables
        if ring.variables:
            vars_str = ", ".join(ring.variables)
            lines.append(f"{indent}  Variables ({len(ring.variables)}): {vars_str}")

        # Ordering
        if ring.ordering_blocks:
            lines.append(f"{indent}  Ordering ({len(ring.ordering_blocks)} blocks):")
            for i, block in enumerate(ring.ordering_blocks):
                lines.append(f"{indent}    [{i}] {self.format_ordering_block(block)}")

        # Extension ring
        if ring.extension_ring:
            lines.append(f"{indent}  Extension ring:")
            lines.extend(self.format_ring(ring.extension_ring, indent + "    "))

        # Quotient ideal
        if ring.quotient_ideal and ring.quotient_ideal.generators:
            lines.append(f"{indent}  Quotient ideal: {len(ring.quotient_ideal.generators)} generators")

        # Properties
        if ring.is_plural:
            lines.append(f"{indent}  Properties: Plural (non-commutative)")
        if ring.is_letterplace:
            lines.append(f"{indent}  Properties: Letterplace")
        if ring.bitmask is not None:
            lines.append(f"{indent}  Bitmask: 0x{ring.bitmask:x}")

        return lines

    def format_poly(self, poly: SSIPoly, indent=""):
        """Format polynomial (summary or detailed)."""
        if not poly.terms:
            return [f"{indent}Polynomial: 0"]

        lines = [f"{indent}Polynomial: {len(poly.terms)} terms"]

        if self.verbose and len(poly.terms) <= 10:
            for i, term in enumerate(poly.terms):
                coef = term.coefficient
                exp_str = "".join([f"*x{j}^{e}" if e > 0 else ""
                                   for j, e in enumerate(term.exponents) if e > 0])
                if not exp_str:
                    exp_str = "1"
                else:
                    exp_str = exp_str.lstrip("*")

                if term.component > 0:
                    lines.append(f"{indent}  [{i}] {coef}*{exp_str} * gen({term.component})")
                else:
                    lines.append(f"{indent}  [{i}] {coef}*{exp_str}")
        elif self.verbose:
            lines.append(f"{indent}  (too many terms to display)")

        return lines

    def format_ideal(self, ideal: SSIIdeal, indent=""):
        """Format ideal (summary)."""
        lines = [f"{indent}Ideal: {len(ideal.generators)} generators"]

        if self.verbose and len(ideal.generators) <= 5:
            for i, gen in enumerate(ideal.generators):
                lines.append(f"{indent}  Generator {i}:")
                lines.extend(self.format_poly(gen, indent + "    "))
        elif ideal.generators:
            # Show term counts
            term_counts = [len(g.terms) for g in ideal.generators]
            lines.append(f"{indent}  Terms per generator: {term_counts}")

        return lines

    def format_module(self, module: SSIModule, indent=""):
        """Format module (summary)."""
        lines = [f"{indent}Module: rank={module.rank}, {len(module.generators)} generators"]

        if self.verbose and len(module.generators) <= 3:
            for i, gen in enumerate(module.generators):
                lines.append(f"{indent}  Generator {i}:")
                lines.extend(self.format_poly(gen, indent + "    "))

        return lines

    def format_matrix(self, matrix: SSIMatrix, indent=""):
        """Format matrix (summary)."""
        lines = [f"{indent}Matrix: {matrix.rows}x{matrix.cols}"]

        if self.verbose and matrix.rows * matrix.cols <= 9:
            for i in range(matrix.rows):
                for j in range(matrix.cols):
                    entry = matrix[i, j]
                    lines.append(f"{indent}  [{i},{j}]: {len(entry.terms)} terms")

        return lines

    def format_command(self, cmd: SSICommand, indent=""):
        """Format command."""
        lines = [f"{indent}Command: opcode={cmd.opcode}, {len(cmd.arguments)} arguments"]

        if self.verbose and len(cmd.arguments) <= 3:
            for i, arg in enumerate(cmd.arguments):
                lines.append(f"{indent}  Arg {i}: {type(arg).__name__}")

        return lines

    def format_list(self, lst: list, indent=""):
        """Format list."""
        lines = [f"{indent}List: {len(lst)} elements"]

        if self.verbose and len(lst) <= 5:
            for i, item in enumerate(lst):
                lines.append(f"{indent}  [{i}] {type(item).__name__}")

        return lines

    def format_value(self, value: Any, indent=""):
        """Format any SSI value."""
        lines = []

        if isinstance(value, SSIVersion):
            lines.append(f"{indent}Handshake: version={value.version}, max_tok={value.max_tok}, opt1=0x{value.opt1:x}, opt2=0x{value.opt2:x}")

        elif isinstance(value, SSIRing):
            lines.extend(self.format_ring(value, indent))
            self.current_ring = value

        elif isinstance(value, SSIPoly):
            lines.extend(self.format_poly(value, indent))

        elif isinstance(value, SSIIdeal):
            lines.extend(self.format_ideal(value, indent))

        elif isinstance(value, SSIModule):
            lines.extend(self.format_module(value, indent))

        elif isinstance(value, SSIMatrix):
            lines.extend(self.format_matrix(value, indent))

        elif isinstance(value, SSICommand):
            lines.extend(self.format_command(value, indent))

        elif isinstance(value, list):
            lines.extend(self.format_list(value, indent))

        elif isinstance(value, SSIProc):
            preview = value.body[:50] + "..." if len(value.body) > 50 else value.body
            lines.append(f"{indent}Procedure: {repr(preview)}")

        elif isinstance(value, SSIDef):
            lines.append(f"{indent}Def: {value.name}")

        elif isinstance(value, SSINone):
            lines.append(f"{indent}None")

        elif isinstance(value, SSIQuit):
            lines.append(f"{indent}Quit")

        elif isinstance(value, SSIBlackbox):
            lines.append(f"{indent}Blackbox: {value.typename}")

        elif isinstance(value, SSIAttributed):
            lines.append(f"{indent}Attributed: flags=0x{value.flags:x}")
            lines.extend(self.format_value(value.value, indent + "  "))

        elif isinstance(value, int):
            lines.append(f"{indent}Int: {value}")

        elif isinstance(value, str):
            preview = value[:50] + "..." if len(value) > 50 else value
            lines.append(f"{indent}String: {repr(preview)}")

        elif isinstance(value, list) and all(isinstance(x, int) for x in value):
            if len(value) <= 10:
                lines.append(f"{indent}IntVec: {value}")
            else:
                lines.append(f"{indent}IntVec: {len(value)} elements")

        else:
            lines.append(f"{indent}{type(value).__name__}: {value}")

        return lines

    def print_message(self, value: Any):
        """Print a single SSI message."""
        self.message_count += 1

        print(f"\n{'='*70}")
        print(f"Message {self.message_count}: {type(value).__name__}")
        print('='*70)

        for line in self.format_value(value):
            print(line)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 ssi_trace_viewer.py <trace_file.ssi> [--verbose]")
        print()
        print("Displays SSI trace files in a human-readable format.")
        print("  --verbose: Show detailed polynomial/ideal contents")
        sys.exit(1)

    filename = sys.argv[1]
    verbose = "--verbose" in sys.argv or "-v" in sys.argv

    viewer = SSITraceViewer(verbose=verbose)

    print(f"Reading SSI trace file: {filename}")
    print(f"Verbose mode: {'ON' if verbose else 'OFF'}")

    try:
        with open(filename, 'rb') as f:
            reader = SSIReader(f)

            while True:
                try:
                    value = reader.read_message()
                    if isinstance(value, SSIQuit):
                        viewer.print_message(value)
                        break
                    viewer.print_message(value)
                except EOFError:
                    break

        print(f"\n{'='*70}")
        print(f"Total messages: {viewer.message_count}")
        print('='*70)

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
