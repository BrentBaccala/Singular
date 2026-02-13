#!/usr/bin/env python3
"""
SSI Trace Ring Viewer - Extract and display ring information from SSI traces
"""

import sys
import re

def extract_rings(filename):
    """Extract ring definitions from SSI trace file."""

    with open(filename, 'r') as f:
        content = f.read()

    # Split into messages (newline-separated, but some messages span multiple lines)
    # We'll process it as a token stream
    tokens = content.split()

    rings = []
    i = 0
    message_num = 0

    while i < len(tokens):
        try:
            type_code = int(tokens[i])
        except:
            i += 1
            continue

        if type_code == 98:  # VERSION
            message_num += 1
            i += 5  # skip version message (type + 4 parameters)

        elif type_code == 15:  # RING_DATA
            message_num += 1
            i += 1  # skip type code

            ring_info = {"message": message_num}

            # Read characteristic
            ch = int(tokens[i])
            i += 1

            # Handle caching
            if ch == -6:
                cache_idx = int(tokens[i])
                ring_info["cache_index"] = cache_idx
                i += 1
                ch = int(tokens[i])
                i += 1
            elif ch == -5:
                ref_idx = int(tokens[i])
                ring_info["cache_ref"] = ref_idx
                i += 2
                rings.append(ring_info)
                continue

            # Characteristic
            if ch == 0:
                ring_info["coeff"] = "Q"
            elif ch > 0:
                ring_info["coeff"] = f"Z/{ch}"
            else:
                ring_info["coeff"] = f"special({ch})"

            # Number of variables
            n_vars = int(tokens[i])
            ring_info["nvars"] = n_vars
            i += 1

            # Variable names
            vars = []
            for _ in range(n_vars):
                length = int(tokens[i])
                i += 1
                varname = tokens[i]
                i += 1
                vars.append(varname)

            ring_info["variables"] = vars

            # Skip ordering blocks (we won't parse them in detail)
            n_blocks = int(tokens[i])
            ring_info["n_ordering_blocks"] = n_blocks
            i += 1

            # Heuristic: skip ahead past the ordering blocks
            # This is approximate but works for most cases
            # We'll just skip until we hit something that looks like data
            skip_count = 0
            max_skip = 50  # safety limit
            while skip_count < max_skip and i < len(tokens):
                try:
                    next_tok = int(tokens[i])
                    # If we see a type code that looks like data (7=IDEAL, 10=MODULE, etc.)
                    # and the next token is a reasonable generator count, we're done with ring
                    if next_tok in (7, 10, 14, 16) and i+1 < len(tokens):
                        try:
                            next_next = int(tokens[i+1])
                            if 0 <= next_next <= 1000:  # reasonable gen count
                                break
                        except:
                            pass
                    i += 1
                    skip_count += 1
                except:
                    i += 1
                    skip_count += 1

            rings.append(ring_info)

        elif type_code == 14:  # LIST
            message_num += 1
            i += 2  # skip type and length

        elif type_code in (1, 2, 7, 10, 16):
            message_num += 1
            i += 1  # skip for now

        else:
            i += 1

    return rings

def print_rings(rings):
    """Print extracted ring information."""
    print("="*70)
    print("Ring Definitions in SSI Trace")
    print("="*70)

    for i, ring in enumerate(rings, 1):
        print(f"\nRing {i} (Message {ring['message']}):")

        if "cache_ref" in ring:
            print(f"  Reference to cached ring: index {ring['cache_ref']}")
            continue

        if "cache_index" in ring:
            print(f"  Cached as: index {ring['cache_index']}")

        print(f"  Coefficient field: {ring['coeff']}")
        print(f"  Number of variables: {ring['nvars']}")
        print(f"  Variables: {', '.join(ring['variables'])}")
        print(f"  Ordering blocks: {ring['n_ordering_blocks']}")

    print(f"\n{'='*70}")
    print(f"Total rings found: {len(rings)}")
    print("="*70)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 ssi_trace_rings.py <trace_file.ssi>")
        sys.exit(1)

    filename = sys.argv[1]

    try:
        rings = extract_rings(filename)
        print_rings(rings)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
