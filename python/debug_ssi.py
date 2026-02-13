#!/usr/bin/env python3
"""Debug SSI file reading"""

import sys

with open(sys.argv[1], 'rb') as f:
    # Read all tokens
    tokens = []
    while True:
        token = bytearray()
        # Skip whitespace
        while True:
            b = f.read(1)
            if not b:
                break
            if b[0] not in (ord(' '), ord('\n'), ord('\r'), ord('\t')):
                token.append(b[0])
                break

        if not token:
            break

        # Read until whitespace
        while True:
            b = f.read(1)
            if not b:
                break
            if b[0] in (ord(' '), ord('\n'), ord('\r'), ord('\t')):
                break
            token.append(b[0])

        if token:
            tokens.append(token.decode('ascii'))

    # Print first 100 tokens
    for i, tok in enumerate(tokens[:100]):
        print(f"[{i}] {tok}")
