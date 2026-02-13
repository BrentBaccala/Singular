#!/usr/bin/env python3
"""
SSI Trace Message Viewer - Shows message-by-message structure
"""

import sys

def parse_trace(filename):
    """Parse SSI trace file message by message."""

    with open(filename, 'r') as f:
        lines = f.readlines()

    messages = []

    for line_num, line in enumerate(lines, 1):
        tokens = line.strip().split()
        if not tokens:
            continue

        try:
            type_code = int(tokens[0])
        except:
            continue

        msg_info = {
            'line': line_num,
            'type_code': type_code,
            'tokens': tokens
        }

        if type_code == 98:  # VERSION
            msg_info['type'] = 'VERSION'
            if len(tokens) >= 5:
                msg_info['version'] = int(tokens[1])
                msg_info['max_tok'] = int(tokens[2])
                msg_info['opt1'] = int(tokens[3])
                msg_info['opt2'] = int(tokens[4])

        elif type_code == 15:  # RING_DATA
            msg_info['type'] = 'RING_DATA'
            # Parse ring definition
            idx = 1
            ch = int(tokens[idx])
            idx += 1

            if ch == -6:
                msg_info['cache_index'] = int(tokens[idx])
                idx += 1
                ch = int(tokens[idx])
                idx += 1

            msg_info['characteristic'] = ch

            if ch >= 0:  # Not a reference
                n_vars = int(tokens[idx])
                msg_info['nvars'] = n_vars
                idx += 1

                # Read variable names
                vars = []
                for _ in range(n_vars):
                    length = int(tokens[idx])
                    idx += 1
                    varname = tokens[idx]
                    idx += 1
                    vars.append(varname)
                msg_info['variables'] = vars

                # Parse ordering blocks
                if idx < len(tokens):
                    n_blocks = int(tokens[idx])
                    msg_info['n_blocks'] = n_blocks
                    idx += 1

                    orderings = []
                    for _ in range(n_blocks):
                        if idx >= len(tokens):
                            break
                        ord_type = int(tokens[idx])
                        idx += 1
                        block0 = int(tokens[idx])
                        idx += 1
                        block1 = int(tokens[idx])
                        idx += 1

                        # Read weights for weighted orderings
                        weights = []
                        if ord_type in (3, 4, 5, 6, 8, 9):  # wp, Wp, ws, Ws, a, aa
                            # NOTE: According to ssiLink.cc, n_weights = block1 - block0 + 1
                            # However, empirical SSI data shows this may not be accurate for
                            # ord_type 8 (a) and 9 (aa) orderings. May need special handling.
                            n_weights = block1 - block0 + 1
                            for _ in range(n_weights):
                                if idx < len(tokens):
                                    weights.append(int(tokens[idx]))
                                    idx += 1
                        elif ord_type == 7:  # M (matrix ordering)
                            n = block1 - block0 + 1
                            n_weights = n * n
                            for _ in range(n_weights):
                                if idx < len(tokens):
                                    weights.append(int(tokens[idx]))
                                    idx += 1

                        orderings.append({
                            'type': ord_type,
                            'block0': block0,
                            'block1': block1,
                            'weights': weights
                        })

                    msg_info['orderings'] = orderings

        elif type_code == 5:  # RING
            msg_info['type'] = 'RING'

        elif type_code == 7:  # IDEAL
            msg_info['type'] = 'IDEAL'
            if len(tokens) >= 2:
                msg_info['n_generators'] = int(tokens[1])

        elif type_code == 10:  # MODULE
            msg_info['type'] = 'MODULE'
            if len(tokens) >= 3:
                msg_info['rank'] = int(tokens[1])
                msg_info['n_generators'] = int(tokens[2])

        elif type_code == 14:  # LIST
            msg_info['type'] = 'LIST'
            if len(tokens) >= 2:
                msg_info['length'] = int(tokens[1])

        elif type_code == 1:  # INT
            msg_info['type'] = 'INT'
            if len(tokens) >= 2:
                msg_info['value'] = int(tokens[1])

        elif type_code == 2:  # STRING
            msg_info['type'] = 'STRING'

        elif type_code == 16:  # NONE
            msg_info['type'] = 'NONE'

        elif type_code == 99:  # QUIT
            msg_info['type'] = 'QUIT'

        else:
            msg_info['type'] = f'UNKNOWN({type_code})'

        messages.append(msg_info)

    return messages


def print_messages(messages, show_rings=True):
    """Print message summary."""

    print("="*70)
    print("SSI Trace Message Structure")
    print("="*70)

    for i, msg in enumerate(messages, 1):
        print(f"\n{'─'*70}")
        print(f"Message {i} (Line {msg['line']}): {msg['type']}")
        print(f"{'─'*70}")

        if msg['type'] == 'VERSION':
            print(f"  Protocol version: {msg.get('version', '?')}")
            print(f"  Max token: {msg.get('max_tok', '?')}")
            print(f"  Options: opt1=0x{msg.get('opt1', 0):x}, opt2=0x{msg.get('opt2', 0):x}")

        elif msg['type'] == 'RING_DATA':
            if 'cache_index' in msg:
                print(f"  Cache index: {msg['cache_index']}")
            if msg['characteristic'] == 0:
                print(f"  Coefficient field: Q")
            elif msg['characteristic'] > 0:
                print(f"  Coefficient field: Z/{msg['characteristic']}")
            else:
                print(f"  Coefficient field: special({msg['characteristic']})")

            if show_rings and 'variables' in msg:
                print(f"  Variables ({msg['nvars']}): {', '.join(msg['variables'])}")

                # Display orderings
                if 'orderings' in msg:
                    ord_names = {
                        0: 'lp', 1: 'dp', 2: 'Dp', 3: 'wp', 4: 'Wp',
                        5: 'ws', 6: 'Ws', 7: 'M', 8: 'a', 9: 'aa',
                        10: 'ls', 11: 'rs', 12: 'ds', 13: 'Ds',
                        # Component orderings
                        2: 'C', 4: 'c'  # Note: These may overlap with Dp(2) and Wp(4)
                    }

                    # For component orderings, check if block0/block1 suggest it's a component ordering
                    def get_ord_name(ord_type, block0, block1):
                        if ord_type in (2, 4) and block0 > 100:  # Heuristic for component ordering
                            return 'C' if ord_type == 2 else 'c'
                        return ord_names.get(ord_type, f'ord{ord_type}')

                    ord_strs = []
                    for ord_info in msg['orderings']:
                        ord_type = ord_info['type']
                        block0 = ord_info['block0']
                        block1 = ord_info['block1']
                        ord_name = get_ord_name(ord_type, block0, block1)

                        if ord_info['weights']:
                            if ord_type == 7:  # Matrix ordering
                                ord_strs.append(f"{ord_name}({block0}..{block1}, matrix)")
                            else:
                                weights_str = ','.join(map(str, ord_info['weights']))
                                ord_strs.append(f"{ord_name}({weights_str})")
                        else:
                            if block0 == block1 and block0 == 0:
                                ord_strs.append(ord_name)
                            elif ord_name in ('C', 'c'):
                                ord_strs.append(ord_name)
                            else:
                                ord_strs.append(f"{ord_name}({block0}..{block1})")

                    print(f"  Ordering: ({', '.join(ord_strs)})")
                else:
                    print(f"  Ordering blocks: {msg.get('n_blocks', '?')}")

        elif msg['type'] == 'IDEAL':
            n_gens = msg.get('n_generators', '?')
            n_tokens = len(msg['tokens'])
            print(f"  Generators: {n_gens}")
            print(f"  Total tokens: {n_tokens}")

        elif msg['type'] == 'MODULE':
            print(f"  Rank: {msg.get('rank', '?')}")
            print(f"  Generators: {msg.get('n_generators', '?')}")

        elif msg['type'] == 'LIST':
            print(f"  Length: {msg.get('length', '?')}")

        elif msg['type'] == 'INT':
            print(f"  Value: {msg.get('value', '?')}")

        elif msg['type'] == 'STRING':
            print(f"  (string value not parsed)")

    print(f"\n{'='*70}")
    print(f"Total messages: {len(messages)}")
    print("="*70)


def print_summary(messages):
    """Print compact summary."""

    print("\n" + "="*70)
    print("Message Summary")
    print("="*70)

    from collections import Counter
    type_counts = Counter(msg['type'] for msg in messages)

    print("\nMessage type counts:")
    for msg_type, count in sorted(type_counts.items()):
        print(f"  {msg_type}: {count}")

    print("\nMessage sequence pattern:")
    pattern = []
    for msg in messages[:20]:  # First 20 messages
        if msg['type'] == 'RING_DATA':
            nvars = msg.get('nvars', '?')
            pattern.append(f"RING[{nvars}vars]")
        elif msg['type'] == 'IDEAL':
            ngens = msg.get('n_generators', '?')
            pattern.append(f"IDEAL[{ngens}gen]")
        else:
            pattern.append(msg['type'])

    print("First 20: " + " → ".join(pattern))

    print("\nRing contexts used:")
    ring_num = 0

    ord_names = {
        0: 'lp', 1: 'dp', 2: 'Dp', 3: 'wp', 4: 'Wp',
        5: 'ws', 6: 'Ws', 7: 'M', 8: 'a', 9: 'aa',
        10: 'ls', 11: 'rs', 12: 'ds', 13: 'Ds'
    }

    def get_ord_name(ord_type, block0, block1):
        if ord_type in (2, 4) and block0 > 100:
            return 'C' if ord_type == 2 else 'c'
        return ord_names.get(ord_type, f'ord{ord_type}')

    for msg in messages:
        if msg['type'] == 'RING_DATA' and 'variables' in msg:
            ring_num += 1
            vars = msg.get('variables', [])

            # Format orderings
            ord_str = ""
            if 'orderings' in msg:
                ord_parts = []
                for ord_info in msg['orderings']:
                    ord_type = ord_info['type']
                    block0 = ord_info['block0']
                    block1 = ord_info['block1']
                    ord_name = get_ord_name(ord_type, block0, block1)

                    if ord_info['weights']:
                        weights_str = ','.join(map(str, ord_info['weights']))
                        ord_parts.append(f"{ord_name}({weights_str})")
                    elif ord_name in ('C', 'c'):
                        ord_parts.append(ord_name)
                    else:
                        ord_parts.append(ord_name)

                ord_str = f", ({', '.join(ord_parts)})"

            print(f"  Ring {ring_num}: Q[{', '.join(vars)}]{ord_str}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 ssi_trace_messages.py <trace_file.ssi> [--full]")
        print("  --full: Show full message details (default: summary only)")
        sys.exit(1)

    filename = sys.argv[1]
    full = '--full' in sys.argv

    try:
        messages = parse_trace(filename)

        if full:
            print_messages(messages, show_rings=True)

        print_summary(messages)

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
