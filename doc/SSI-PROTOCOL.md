# SSI Protocol Specification

## Singular Serialization Interface (SSI) Protocol

**Version:** 15
**Document Version:** 1.0
**Based on:** Singular source code analysis (ssiLink.cc, longrat.cc, modulop.cc)

---

## 1. Overview

SSI is a text-based serialization protocol used by the Singular computer algebra system
for inter-process communication and data persistence. It encodes Singular's algebraic
data types (polynomials, ideals, rings, etc.) in a human-readable, space-separated format.

### 1.1 Design Principles

- **Text-based:** All data encoded as ASCII text with space separators
- **Self-describing:** Each value is prefixed with a type code
- **Streaming:** Values can be read/written incrementally
- **Ring-aware:** Polynomials and related types include ring context

### 1.2 Basic Syntax

- All tokens are separated by spaces (ASCII 0x20)
- Messages are terminated by newline (ASCII 0x0A)
- Integers are encoded in decimal
- Large integers (mpz_t) are encoded in hexadecimal (base 16)
- Strings are length-prefixed

---

## 2. Transport Layer

### 2.1 Connection Modes

| Mode | Description |
|------|-------------|
| `ssi:r filename` | Read from file |
| `ssi:w filename` | Write to file (overwrite) |
| `ssi:a filename` | Write to file (append) |
| `ssi:tcp` | Server mode: listen on first available port (1025-50000) |
| `ssi:connect host:port` | Client mode: connect to server |
| `ssi:tcp host:path` | Launch mode: start remote Singular via SSH |
| `ssi:fork` | Fork mode: create child process |

### 2.2 TCP Connection Establishment

**Server mode (`ssi:tcp`):**
1. Create TCP socket (AF_INET, SOCK_STREAM)
2. Bind to first available port in range 1025-50000
3. Listen with backlog of 1
4. Accept single connection
5. Send handshake message (type 98)

**Client mode (`ssi:connect`):**
1. Parse `host:port` from link name
2. Create TCP socket
3. Resolve hostname via gethostbyname()
4. Connect to server
5. Receive and validate handshake

### 2.3 Communication Model

In fork/tcp modes, Singular uses a **master-worker** model:

- **Master (client):** Runs user code, sends commands, receives results
- **Worker (server):** Runs read-eval-write loop:
  ```
  loop {
      command = read(link)
      result = evaluate(command)
      write(link, result)
  }
  ```

For peer-to-peer operation, both endpoints can use `ssi:connect` mode
via a relay (e.g., socat), but coordination must be handled by the application.

---

## 3. Handshake

### 3.1 Handshake Message (Type 98)

Sent by server immediately after connection establishment:

```
98 <version> <max_tok> <opt1> <opt2>\n
```

| Field | Type | Description |
|-------|------|-------------|
| version | int | SSI protocol version (currently 15) |
| max_tok | int | Maximum token value for commands |
| opt1 | uint | si_opt_1 options bitmask |
| opt2 | uint | si_opt_2 options bitmask |

**Example:**
```
98 15 512 0 0
```

### 3.2 Version Compatibility

Client should verify version compatibility. Version history:
- 5→6: Changed newstruct representation
- 6→7: Added attributes
- 7→8: Added qring
- 8→9: Added module rank
- 9→10: Reorganized tokens in grammar.h/tok.h
- 10→11: Extended ring description for named coefficients
- 11→12: Added rank to ideal/module, added smatrix
- 12→13: Added NC (non-commutative) rings
- 13→14: Added ring references
- 14→15: Added bigintvec, prune_map, mres_map

---

## 4. Message Format

### 4.1 General Structure

Each message consists of:
```
[<attribute>] <type_code> <type_specific_data>\n
```

Optional attribute prefix (type 21) may precede the main type.

### 4.2 Type Codes

| Code | Type | Description |
|------|------|-------------|
| 1 | INT | Machine integer |
| 2 | STRING | Character string |
| 3 | NUMBER | Coefficient (ring-dependent) |
| 4 | BIGINT | Arbitrary precision integer |
| 5 | RING | Ring definition |
| 6 | POLY | Polynomial |
| 7 | IDEAL | Ideal (list of polynomials) |
| 8 | MATRIX | Matrix of polynomials |
| 9 | VECTOR | Vector (polynomial with module component) |
| 10 | MODULE | Module (ideal with rank) |
| 11 | COMMAND | Unevaluated command |
| 12 | DEF | Undefined/quoted name |
| 13 | PROC | Procedure definition |
| 14 | LIST | Heterogeneous list |
| 15 | RING+DATA | Ring context switch followed by data |
| 16 | NONE | Void/nothing |
| 17 | INTVEC | Integer vector |
| 18 | INTMAT | Integer matrix |
| 19 | BIGINTMAT | Big integer matrix |
| 20 | BLACKBOX | Custom/plugin type |
| 21 | ATTRIB | Attributes prefix |
| 22 | SMATRIX | Sparse matrix |
| 23 | RINGPROP | Ring properties (bitmask, LP, NC) |
| 24 | BIGINTVEC | Big integer vector |
| 98 | VERSION | Handshake (see section 3) |
| 99 | QUIT | Close connection |

---

## 5. Scalar Types

### 5.1 Integer (Type 1)

```
1 <int>
```

**Example:** `1 42` represents the integer 42.

### 5.2 String (Type 2)

```
2 <length> <characters>
```

| Field | Type | Description |
|-------|------|-------------|
| length | int | Number of characters |
| characters | chars | Exactly `length` characters (may contain spaces) |

**Example:** `2 5 hello` represents the string "hello".

### 5.3 Big Integer (Type 4)

Uses the rational number encoding (see section 6.1) with subtype constraints.

```
4 <rational_encoding>
```

The rational must have denominator 1 (i.e., subtypes 4 or 8 only).

---

## 6. Coefficient Encoding

Coefficients are encoded according to the ring's coefficient type.

### 6.1 Rational Numbers (Q)

Subtype code determines encoding:

| Subtype | Description | Format |
|---------|-------------|--------|
| 0 | Unnormalized fraction | `0 <mpz_num> <mpz_den>` |
| 1 | Normalized fraction | `1 <mpz_num> <mpz_den>` |
| 3 | Integer (old format) | `3 <mpz>` |
| 4 | Small integer | `4 <int>` |
| 5 | Raw unnormalized | `5 <hex_num> <hex_den>` |
| 6 | Raw normalized | `6 <hex_num> <hex_den>` |
| 8 | Raw integer | `8 <hex>` |

Where:
- `<int>` is a decimal machine integer
- `<mpz>` is a GMP integer in decimal
- `<hex>` is a GMP integer in hexadecimal (base 16)

**Examples:**
```
4 42           # Small integer 42
8 1a           # Integer 26 (hex)
6 3 4          # Fraction 3/4 (normalized, hex)
```

### 6.2 Modular Integers (Z/p)

For prime field Z/p, coefficient is a single integer in range [0, p-1]:

```
<int>
```

**Example:** In Z/7, the element 5 is encoded as just `5`.

### 6.3 Transcendental Extension

For Q(t1,...,tn), coefficient is encoded as fraction of polynomials:

```
<poly_numerator> <poly_denominator>
```

Each polynomial is encoded per section 7.2, over the extension ring.

### 6.4 Algebraic Extension

For Q[x]/(minpoly), coefficient is encoded as a single polynomial:

```
<poly>
```

---

## 7. Polynomial Types

### 7.1 Ring Context

Polynomial types require an active ring. If the current ring differs from
the data's ring, a ring context switch (type 15) precedes the data:

```
15 <ring_definition> <actual_data>
```

### 7.2 Polynomial (Type 6)

```
6 <nterms> [<term>]*
```

Where each term is:
```
<coefficient> <component> <exp1> <exp2> ... <expN>
```

| Field | Type | Description |
|-------|------|-------------|
| nterms | int | Number of terms |
| coefficient | number | See section 6 |
| component | int | Module component (0 for polynomials) |
| exp1..expN | int | Exponents for each variable |

**Example:** In ring Q[x,y], polynomial `3*x^2*y + 2`:
```
6 2 4 3 0 2 1 4 2 0 0 0
```
Breakdown:
- `6` - type code (POLY)
- `2` - two terms
- `4 3 0 2 1` - term 1: coeff=3 (small int), comp=0, x^2, y^1
- `4 2 0 0 0` - term 2: coeff=2, comp=0, x^0, y^0

### 7.3 Vector (Type 9)

Same as polynomial, but component > 0:

```
9 <nterms> [<term>]*
```

### 7.4 Ideal (Type 7)

```
7 <ngenerators> [<poly>]*
```

| Field | Type | Description |
|-------|------|-------------|
| ngenerators | int | Number of generators |
| poly | polynomial | Each generator (without type prefix) |

### 7.5 Module (Type 10)

```
10 <rank> <ngenerators> [<vector>]*
```

| Field | Type | Description |
|-------|------|-------------|
| rank | int | Free module rank |
| ngenerators | int | Number of generators |
| vector | vector | Each generator |

### 7.6 Matrix (Type 8)

```
8 <rows> <cols> [<poly>]*
```

Entries are in row-major order (row 1 col 1, row 1 col 2, ...).

### 7.7 Sparse Matrix (Type 22)

```
22 <rank> <ngenerators> [<vector>]*
```

Same format as module.

---

## 8. Ring Encoding (Type 5)

### 8.1 Ring Definition

```
5 <ch> <N> [<varname>]* <nblocks> [<ordering_block>]* [<extension>] <qideal>
```

| Field | Type | Description |
|-------|------|-------------|
| ch | int | Characteristic (see 8.2) |
| N | int | Number of variables |
| varname | string | Variable name (length-prefixed) |
| nblocks | int | Number of ordering blocks |
| ordering_block | block | See section 8.3 |
| extension | ring | For ch=-1,-2: coefficient ring |
| qideal | ideal | Quotient ideal (0 if none) |

### 8.2 Characteristic Codes

| Code | Meaning |
|------|---------|
| 0 | Rational numbers Q |
| p > 0 | Prime field Z/p |
| -1 | Transcendental extension (ring follows) |
| -2 | Algebraic extension (ring with qideal follows) |
| -3 | Named coefficient field (name follows) |
| -4 | NULL ring |
| -5 | Ring reference (cached) |
| -6 | New ring reference (cache and define) |

### 8.3 Ordering Block

```
<ord_type> <block0> <block1> [<weights>]*
```

| Field | Type | Description |
|-------|------|-------------|
| ord_type | int | Ordering type code |
| block0 | int | First variable index |
| block1 | int | Last variable index |
| weights | int[] | Weight vector (for weighted orderings) |

**Ordering Types:**

| Code | Ordering | Weights |
|------|----------|---------|
| 0 | lp (lex) | none |
| 1 | dp (degrevlex) | none |
| 2 | Dp (deglex) | none |
| 3 | wp (weighted degrevlex) | n weights |
| 4 | Wp (weighted deglex) | n weights |
| 5 | ws (weighted revlex, neg) | n weights |
| 6 | Ws (weighted lex, neg) | n weights |
| 7 | M (matrix) | n×n weights |
| 8 | a (extra weight) | n weights |
| ... | (others) | ... |

### 8.4 Ring Caching

To avoid retransmitting ring definitions, SSI maintains a cache of 20 rings.

**New reference:** `-6 <index> <ring_definition>`
Stores ring in cache slot `index`, then defines it.

**Cached reference:** `-5 <index>`
References previously cached ring.

### 8.5 Ring Properties (Type 23)

Additional ring properties after the main definition:

```
23 <subtype> <data>
```

| Subtype | Description | Data |
|---------|-------------|------|
| 0 | Custom bitmask | `<log2_bitmask>` |
| 1 | Letterplace ring | `<log2_bitmask> <isLPring>` |
| 2 | Plural (NC) ring | `<matrix_C> <matrix_D>` |

---

## 9. Compound Types

### 9.1 List (Type 14)

```
14 <length> [<element>]*
```

Elements are self-describing (include their type codes).

### 9.2 Integer Vector (Type 17)

```
17 <length> [<int>]*
```

### 9.3 Integer Matrix (Type 18)

```
18 <rows> <cols> [<int>]*
```

Row-major order.

### 9.4 Big Integer Matrix (Type 19)

```
19 <rows> <cols> [<bigint>]*
```

### 9.5 Big Integer Vector (Type 24)

```
24 <length> [<bigint>]*
```

---

## 10. Command Encoding (Type 11)

Commands represent unevaluated expressions (created by `quote()`).

```
11 <argc> <opcode> [<arg>]*
```

| Field | Type | Description |
|-------|------|-------------|
| argc | int | Argument count (1-3, or more) |
| opcode | int | Operation token number |
| arg | value | Self-describing argument values |

For argc ≤ 3, arguments are positional (arg1, arg2, arg3).
For argc > 3, arguments form a linked list from arg1.

### 10.1 Command Evaluation

When a COMMAND is read by the interpreter:
- In the worker's read-eval-write loop: automatically evaluated
- In peer-to-peer mode: evaluated when the read() result is used
- Can be explicitly evaluated with `eval()`

---

## 11. Attributes (Type 21)

Attributes prefix any other value:

```
21 <flags> <nattrs> [<attr_data>]* <value>
```

| Field | Type | Description |
|-------|------|-------------|
| flags | int | Bitset of attribute flags |
| nattrs | int | Number of named attributes |
| attr_data | varies | Named attribute data |
| value | value | The attributed value (any type) |

---

## 12. Special Types

### 12.1 None (Type 16)

Represents void/nothing:

```
16
```

### 12.2 Procedure (Type 13)

```
13 <body>
```

Where `body` is a length-prefixed string containing the procedure source.

### 12.3 Def (Type 12)

Undefined/quoted name:

```
12 <name>
```

Where `name` is a length-prefixed string.

### 12.4 Blackbox (Type 20)

Custom types registered by plugins:

```
20 <typename> <serialized_data>
```

The typename is a string; serialized_data depends on the blackbox implementation.

---

## 13. Control Messages

### 13.1 Quit (Type 99)

Signals connection close:

```
99
```

Upon receiving type 99, the recipient should close the connection and exit
(in batch/fork mode).

---

## 14. Examples

### 14.1 Simple Integer

```
1 42
```

### 14.2 Polynomial in Q[x,y,z]

Ring definition followed by polynomial x^2 + y*z:

```
5 0 3 1 x 1 y 1 z 2 1 1 3 2 4 6 0 0 0 6 2 4 1 0 2 0 0 4 1 0 0 1 1
```

Breakdown:
- `5` - RING type
- `0 3` - char=0 (Q), 3 variables
- `1 x 1 y 1 z` - variable names
- `2` - 2 ordering blocks
- `1 1 3` - dp ordering, vars 1-3
- `2 4 6` - C ordering (component), positions 4-6
- `0` - no quotient ideal
- `6 2` - POLY with 2 terms
- `4 1 0 2 0 0` - term: coeff=1, comp=0, x^2
- `4 1 0 0 1 1` - term: coeff=1, comp=0, y*z

### 14.3 Quoted Command

`quote(std(I))` where I is an ideal:

```
11 1 285 7 2 6 1 4 1 0 1 0 0 6 1 4 1 0 0 1 0
```

- `11 1 285` - COMMAND with 1 arg, opcode 285 (std)
- `7 2 ...` - IDEAL with 2 generators

---

## 15. Implementation Notes

### 15.1 Buffering

- Write: Standard FILE* buffering, flushed after each complete message
- Read: Custom s_buff buffered reader for efficient parsing

### 15.2 Error Handling

- Parse errors: Return NULL/error type
- Connection errors: Close link, set error flag
- Version mismatch: Warning, attempt to continue

### 15.3 Security Considerations

- No authentication mechanism
- No encryption
- Command execution: received COMMANDs are evaluated
- File access: limited to Singular's file I/O capabilities

---

## Appendix A: Grammar (BNF-like)

```
message     := [attrib] value NEWLINE
attrib      := '21' INT INT {attr_data}*
value       := int_val | string_val | number_val | bigint_val
             | ring_val | poly_val | ideal_val | matrix_val
             | vector_val | module_val | command_val | def_val
             | proc_val | list_val | none_val | intvec_val
             | intmat_val | bigintmat_val | smatrix_val
             | blackbox_val | version_val | quit_val

int_val     := '1' INT
string_val  := '2' INT CHARS
number_val  := '3' coeff_encoding
bigint_val  := '4' rational_encoding
ring_val    := '5' ring_def | '-4' | '-5' INT | '-6' INT ring_def
poly_val    := '6' INT {term}*
ideal_val   := '7' INT {poly_body}*
matrix_val  := '8' INT INT {poly_body}*
vector_val  := '9' INT {term}*
module_val  := '10' INT INT {vector_body}*
command_val := '11' INT INT {value}*
list_val    := '14' INT {value}*
none_val    := '16'
version_val := '98' INT INT INT INT
quit_val    := '99'

ring_def    := INT INT {varname}* INT {ord_block}* [ring_def] ideal_body
varname     := INT CHARS
ord_block   := INT INT INT {INT}*
term        := coeff_encoding INT {INT}+
poly_body   := INT {term}*

coeff_encoding := (depends on coefficient type)
rational_encoding := '4' INT | '8' HEX | '5' HEX HEX | '6' HEX HEX | ...

INT         := decimal integer
HEX         := hexadecimal integer (base 16)
CHARS       := sequence of characters
NEWLINE     := '\n' (ASCII 0x0A)
```

---

## Appendix B: Version History

| Version | Changes |
|---------|---------|
| 15 | bigintvec, prune_map, mres_map |
| 14 | Ring references |
| 13 | NC (non-commutative) rings |
| 12 | Rank in ideal/module, smatrix |
| 11 | Extended ring for named coefficients |
| 10 | Token reorganization |
| 9 | Module rank |
| 8 | Quotient rings |
| 7 | Attributes |
| 6 | Newstruct representation |

---

## Appendix C: References

- Singular source: `Singular/links/ssiLink.cc`
- Coefficient I/O: `libpolys/coeffs/longrat.cc`, `libpolys/coeffs/modulop.cc`
- Buffer handling: `libpolys/reporter/s_buff.h`
- Documentation: `doc/types.doc`, `doc/reference.doc`
