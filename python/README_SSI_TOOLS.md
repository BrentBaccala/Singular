# SSI Tools for Singular

Python tools for reading and analyzing Singular SSI (Singular Serialization Interface) protocol traces.

## Tools

### 1. ssi_reader.py
Complete SSI protocol reader library (original from Claude's documentation commit).

**Usage:**
```python
from ssi_reader import SSIReader, SSIConnection

# Read from file
with open('data.ssi', 'rb') as f:
    reader = SSIReader(f)
    value = reader.read_value()

# Connect to TCP server
conn = SSIConnection('localhost', 12345)
value = conn.read()
conn.close()
```

**Note:** Currently has parsing issues with some coefficient encodings.

### 2. ssi_trace_rings.py (RECOMMENDED)
Simple tool to extract and display ring definitions from SSI trace files.

**Usage:**
```bash
python3 ssi_trace_rings.py <trace_file.ssi>
```

**Output:**
- Ring number and message position
- Coefficient field (Q, Z/p, etc.)
- Number of variables
- Variable names
- Number of ordering blocks

**Example:**
```bash
$ python3 ssi_trace_rings.py hydrogen-5.ssi

Ring 1 (Message 2):
  Cached as: index 0
  Coefficient field: Q
  Number of variables: 11
  Variables: E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1
  Ordering blocks: 2
...
Total rings found: 13
```

### 3. debug_ssi.py
Debug utility that tokenizes an SSI file and displays raw tokens.

**Usage:**
```bash
python3 debug_ssi.py <trace_file.ssi>
```

## Creating SSI Traces

In Singular, use `set_groebner_ssi()` to create traces:

```singular
// Load the library
LIB "standard.lib";

// Set SSI trace target (append mode)
set_groebner_ssi("ssi:a /tmp/trace.ssi");

// Now all groebner() calls will be logged
ideal I = /* your ideal */;
ideal G = groebner(I);

// To stop tracing
set_groebner_ssi("");
```

## SSI Protocol Documentation

See `doc/SSI-PROTOCOL.md` for complete protocol specification (Version 15).

## Known Issues

- `ssi_reader.py` has issues with certain coefficient encodings in polynomial terms
- Ordering block parsing is approximate in the trace tools (doesn't affect ring variables/coefficients)

## Future Enhancements

- Fix coefficient parsing in ssi_reader.py
- Add polynomial pretty-printing
- Add ideal statistics (number of generators, term counts, degrees)
- Support for remote Singular computation via SSI
