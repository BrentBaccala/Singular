#!/usr/bin/env python3
"""
Simple SSI Trace Viewer - Focuses on rings and message summaries
"""

import sys

class SimpleSSIReader:
    """Simple SSI reader that focuses on structure, not perfect parsing."""

    def __init__(self, filename):
        self.filename = filename
        self.tokens = []
        self.pos = 0
        self._read_all_tokens()

    def _read_all_tokens(self):
        """Read all tokens from file."""
        with open(self.filename, 'rb') as f:
            while True:
                token = bytearray()
                # Skip whitespace
                while True:
                    b = f.read(1)
                    if not b:
                        return
                    if b[0] not in (ord(' '), ord('\n'), ord('\r'), ord('\t')):
                        token.append(b[0])
                        break

                if not token:
                    return

                # Read until whitespace
                while True:
                    b = f.read(1)
                    if not b:
                        break
                    if b[0] in (ord(' '), ord('\n'), ord('\r'), ord('\t')):
                        break
                    token.append(b[0])

                if token:
                    self.tokens.append(token.decode('ascii', errors='replace'))

    def peek(self):
        """Peek at current token."""
        if self.pos < len(self.tokens):
            return self.tokens[self.pos]
        return None

    def read(self):
        """Read and advance."""
        if self.pos < len(self.tokens):
            tok = self.tokens[self.pos]
            self.pos += 1
            return tok
        return None

    def read_int(self):
        """Read an integer."""
        tok = self.read()
        try:
            return int(tok)
        except:
            return 0

    def read_string(self):
        """Read a length-prefixed string."""
        length = self.read_int()
        # The next token IS the string (it was already read as a token)
        # Note: this works because SSI strings don't contain spaces
        string_tok = self.read()
        return string_tok if string_tok else ""

    def skip(self, n=1):
        """Skip n tokens."""
        self.pos += n

    def format_ring(self, start_pos):
        """Try to parse and format a ring definition."""
        self.pos = start_pos
        lines = []

        ch = self.read_int()

        # Handle ring caching
        cache_idx = None
        if ch == -6:
            cache_idx = self.read_int()
            ch = self.read_int()
            lines.append(f"  Cached ring reference: index {cache_idx}")
        elif ch == -5:
            ref_idx = self.read_int()
            lines.append(f"  Ring reference: index {ref_idx}")
            return lines, self.pos

        # Characteristic
        if ch == 0:
            lines.append(f"  Coefficient field: Q (rationals)")
        elif ch > 0:
            lines.append(f"  Coefficient field: Z/{ch}")
        elif ch == -1:
            lines.append(f"  Coefficient field: Transcendental extension")
        elif ch == -2:
            lines.append(f"  Coefficient field: Algebraic extension")
        elif ch == -3:
            lines.append(f"  Coefficient field: Named field")
        elif ch == -4:
            lines.append(f"  NULL ring")
            return lines, self.pos

        # Number of variables
        n_vars = self.read_int()
        lines.append(f"  Variables ({n_vars}):")

        # Read variable names
        vars = []
        for i in range(n_vars):
            var = self.read_string()
            vars.append(var)
        lines.append(f"    {', '.join(vars)}")

        # Ordering blocks
        n_blocks = self.read_int()
        lines.append(f"  Ordering ({n_blocks} blocks):")

        ord_names = {
            0: "lp", 1: "dp", 2: "Dp", 3: "wp", 4: "Wp",
            5: "ws", 6: "Ws", 7: "M", 8: "a", 9: "aa"
        }

        for i in range(n_blocks):
            ord_type = self.read_int()
            block0 = self.read_int()
            block1 = self.read_int()

            ord_name = ord_names.get(ord_type, f"ord{ord_type}")

            # Skip weights for weighted orderings
            if ord_type in (3, 4, 5, 6, 8, 9):
                n_weights = block1 - block0 + 1
                weights = [self.read_int() for _ in range(n_weights)]
                lines.append(f"    [{i}] {ord_name}({block0}..{block1}, {n_weights} weights)")
            elif ord_type == 7:  # Matrix ordering
                n_weights = (block1 - block0 + 1) ** 2
                self.skip(n_weights)
                lines.append(f"    [{i}] M({block0}..{block1}, matrix)")
            else:
                lines.append(f"    [{i}] {ord_name}({block0}..{block1})")

        # Quotient ideal (we'll skip parsing it fully)
        n_gens = self.read_int()
        if n_gens > 0:
            lines.append(f"  Quotient ideal: {n_gens} generators (skipped)")
            # Skip the generators
            for _ in range(n_gens):
                n_terms = self.read_int()
                # Skip polynomial terms (rough estimate)
                self.skip(n_terms * (2 + 1 + n_vars))

        return lines, self.pos

    def summarize_trace(self):
        """Summarize the entire trace."""
        self.pos = 0
        message_num = 0

        print("="*70)
        print(f"SSI Trace Summary: {self.filename}")
        print("="*70)

        while self.pos < len(self.tokens):
            type_code_str = self.peek()
            if not type_code_str:
                break

            try:
                type_code = int(type_code_str)
            except:
                self.skip()
                continue

            message_num += 1
            print(f"\nMessage {message_num}:")

            if type_code == 98:  # VERSION
                self.skip()  # Skip type code
                version = self.read_int()
                max_tok = self.read_int()
                opt1 = self.read_int()
                opt2 = self.read_int()
                print(f"  Type: Handshake (VERSION)")
                print(f"  Protocol version: {version}")
                print(f"  Max token: {max_tok}")
                print(f"  Options: opt1=0x{opt1:x}, opt2=0x{opt2:x}")

            elif type_code == 15:  # RING_DATA
                self.skip()  # Skip type code
                print(f"  Type: RING_DATA (ring + data)")
                ring_start = self.pos
                try:
                    ring_lines, new_pos = self.format_ring(ring_start)
                    print("  Ring:")
                    for line in ring_lines:
                        print(line)
                    self.pos = new_pos

                    # Now read the data type (should be next token)
                    data_type_tok = self.peek()
                    if data_type_tok:
                        try:
                            data_type = int(data_type_tok)
                            self.skip()  # consume it
                            if data_type == 7:  # IDEAL
                                n_gens = self.read_int()
                                print(f"  Data: IDEAL with {n_gens} generators")
                                # Skip the entire ideal (approximate)
                            elif data_type == 10:  # MODULE
                                rank = self.read_int()
                                n_gens = self.read_int()
                                print(f"  Data: MODULE rank={rank}, {n_gens} generators")
                            elif data_type == 16:  # NONE
                                print(f"  Data: NONE (empty)")
                            else:
                                print(f"  Data: Type {data_type}")
                        except:
                            print(f"  Data: (couldn't parse data type)")
                except Exception as e:
                    print(f"  (Error parsing: {e})")
                    import traceback
                    traceback.print_exc()

            elif type_code == 5:  # RING
                self.skip()  # Skip type code
                print(f"  Type: RING")
                try:
                    ring_lines, new_pos = self.format_ring(self.pos)
                    for line in ring_lines:
                        print(line)
                    self.pos = new_pos
                except Exception as e:
                    print(f"  (Error parsing ring: {e})")

            elif type_code == 7:  # IDEAL
                self.skip()
                n_gens = self.read_int()
                print(f"  Type: IDEAL with {n_gens} generators")

            elif type_code == 14:  # LIST
                self.skip()
                length = self.read_int()
                print(f"  Type: LIST with {length} elements")

            elif type_code == 1:  # INT
                self.skip()
                val = self.read_int()
                print(f"  Type: INT = {val}")

            elif type_code == 2:  # STRING
                self.skip()
                s = self.read_string()
                print(f"  Type: STRING = {repr(s[:50])}")

            elif type_code == 16:  # NONE
                self.skip()
                print(f"  Type: NONE")

            elif type_code == 99:  # QUIT
                self.skip()
                print(f"  Type: QUIT")
                break

            else:
                self.skip()
                print(f"  Type: {type_code} (unsupported in summary)")

        print(f"\n{'='*70}")
        print(f"Total messages: {message_num}")
        print(f"Total tokens: {len(self.tokens)}")
        print("="*70)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 ssi_trace_simple.py <trace_file.ssi>")
        sys.exit(1)

    filename = sys.argv[1]
    reader = SimpleSSIReader(filename)
    reader.summarize_trace()


if __name__ == '__main__':
    main()
