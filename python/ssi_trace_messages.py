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

                # Skip ordering blocks
                if idx < len(tokens):
                    n_blocks = int(tokens[idx])
                    msg_info['n_blocks'] = n_blocks

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
    for msg in messages:
        if msg['type'] == 'RING_DATA' and 'variables' in msg:
            ring_num += 1
            vars = msg.get('variables', [])
            print(f"  Ring {ring_num}: Q[{', '.join(vars)}]")


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
