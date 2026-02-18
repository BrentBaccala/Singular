#!/usr/bin/env python3
"""
SSI Protocol Reader for Singular

A Python implementation of the Singular Serialization Interface (SSI) protocol
reader. This allows Python programs to read data serialized by Singular.

Based on SSI Protocol Version 15.

Author: Claude (Anthropic) - claude-opus-4-5-20250101
Date: December 2024
License: GNU General Public License (same as Singular)

Usage:
    from ssi_reader import SSIReader, SSIConnection

    # Read from file
    with open('data.ssi', 'rb') as f:
        reader = SSIReader(f)
        value = reader.read_value()

    # Connect to Singular TCP server
    conn = SSIConnection('localhost', 12345)
    value = conn.read()
    conn.write_int(42)
    conn.close()
"""

import socket
import io
from dataclasses import dataclass, field
from typing import Any, Optional, Union, List, Dict, BinaryIO
from fractions import Fraction


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class SSIRing:
    """Represents a Singular ring."""
    characteristic: int
    variables: List[str]
    ordering_blocks: List[tuple]
    quotient_ideal: Optional['SSIIdeal'] = None
    extension_ring: Optional['SSIRing'] = None
    # For algebraic/transcendental extensions
    extension_type: Optional[str] = None  # 'algebraic', 'transcendental', or None
    # Ring properties
    bitmask: Optional[int] = None
    is_letterplace: bool = False
    is_plural: bool = False
    nc_C: Optional['SSIMatrix'] = None  # For plural rings
    nc_D: Optional['SSIMatrix'] = None

    def __repr__(self):
        if self.characteristic == 0:
            cf = "Q"
        elif self.characteristic > 0:
            cf = f"Z/{self.characteristic}"
        else:
            cf = f"ext({self.extension_type})"
        vars_str = ",".join(self.variables)
        return f"SSIRing({cf}[{vars_str}])"


@dataclass
class SSITerm:
    """A single term in a polynomial."""
    coefficient: Any  # int, Fraction, or other coefficient type
    component: int    # Module component (0 for polynomials)
    exponents: List[int]

    def __repr__(self):
        return f"SSITerm(coef={self.coefficient}, comp={self.component}, exp={self.exponents})"


@dataclass
class SSIPoly:
    """A polynomial (list of terms)."""
    terms: List[SSITerm]
    ring: Optional[SSIRing] = None

    def __repr__(self):
        if not self.terms:
            return "SSIPoly(0)"
        return f"SSIPoly({len(self.terms)} terms)"


@dataclass
class SSIIdeal:
    """An ideal (list of polynomials)."""
    generators: List[SSIPoly]
    ring: Optional[SSIRing] = None

    def __repr__(self):
        return f"SSIIdeal({len(self.generators)} generators)"


@dataclass
class SSIModule:
    """A module (ideal with rank)."""
    rank: int
    generators: List[SSIPoly]
    ring: Optional[SSIRing] = None

    def __repr__(self):
        return f"SSIModule(rank={self.rank}, {len(self.generators)} generators)"


@dataclass
class SSIMatrix:
    """A matrix of polynomials."""
    rows: int
    cols: int
    entries: List[SSIPoly]  # Row-major order
    ring: Optional[SSIRing] = None

    def __getitem__(self, idx):
        i, j = idx
        return self.entries[i * self.cols + j]

    def __repr__(self):
        return f"SSIMatrix({self.rows}x{self.cols})"


@dataclass
class SSICommand:
    """An unevaluated command (from quote())."""
    opcode: int
    arguments: List[Any]

    def __repr__(self):
        return f"SSICommand(op={self.opcode}, argc={len(self.arguments)})"


@dataclass
class SSIProc:
    """A Singular procedure."""
    body: str

    def __repr__(self):
        preview = self.body[:50] + "..." if len(self.body) > 50 else self.body
        return f"SSIProc({preview!r})"


@dataclass
class SSIDef:
    """An undefined/quoted name."""
    name: str

    def __repr__(self):
        return f"SSIDef({self.name!r})"


@dataclass
class SSIVersion:
    """Handshake/version information."""
    version: int
    max_tok: int
    opt1: int
    opt2: int

    def __repr__(self):
        return f"SSIVersion(v={self.version}, max_tok={self.max_tok})"


@dataclass
class SSIQuit:
    """Quit/close signal."""
    pass


@dataclass
class SSINone:
    """Void/nothing value."""
    pass


@dataclass
class SSIBlackbox:
    """Custom/plugin type."""
    typename: str
    data: Any

    def __repr__(self):
        return f"SSIBlackbox({self.typename})"


@dataclass
class SSIAttributed:
    """A value with attributes."""
    flags: int
    attributes: Dict[str, Any]
    value: Any

    def __repr__(self):
        return f"SSIAttributed(flags={self.flags}, value={self.value})"


# =============================================================================
# Reader Implementation
# =============================================================================

class SSIReader:
    """
    Reads SSI-encoded data from a binary stream.

    The stream should be opened in binary mode ('rb').
    """

    # Type codes
    TYPE_INT = 1
    TYPE_STRING = 2
    TYPE_NUMBER = 3
    TYPE_BIGINT = 4
    TYPE_RING = 5
    TYPE_POLY = 6
    TYPE_IDEAL = 7
    TYPE_MATRIX = 8
    TYPE_VECTOR = 9
    TYPE_MODULE = 10
    TYPE_COMMAND = 11
    TYPE_DEF = 12
    TYPE_PROC = 13
    TYPE_LIST = 14
    TYPE_RING_DATA = 15
    TYPE_NONE = 16
    TYPE_INTVEC = 17
    TYPE_INTMAT = 18
    TYPE_BIGINTMAT = 19
    TYPE_BLACKBOX = 20
    TYPE_ATTRIB = 21
    TYPE_SMATRIX = 22
    TYPE_RINGPROP = 23
    TYPE_BIGINTVEC = 24
    TYPE_VERSION = 98
    TYPE_QUIT = 99

    def __init__(self, stream: BinaryIO):
        """
        Initialize reader with a binary stream.

        Args:
            stream: A file-like object opened in binary mode
        """
        self.stream = stream
        self.current_ring: Optional[SSIRing] = None
        self.ring_cache: Dict[int, SSIRing] = {}
        self._buffer = b''

    def _read_byte(self) -> int:
        """Read a single byte, returning its integer value."""
        b = self.stream.read(1)
        if not b:
            raise EOFError("Unexpected end of stream")
        return b[0]

    def _peek_byte(self) -> int:
        """Peek at the next byte without consuming it."""
        b = self.stream.read(1)
        if not b:
            raise EOFError("Unexpected end of stream")
        # Put it back by seeking, or use buffer
        if hasattr(self.stream, 'seek'):
            self.stream.seek(-1, 1)
        else:
            self._buffer = b + self._buffer
        return b[0]

    def _read_bytes(self, n: int) -> bytes:
        """Read exactly n bytes."""
        result = self.stream.read(n)
        if len(result) < n:
            raise EOFError(f"Expected {n} bytes, got {len(result)}")
        return result

    def _skip_whitespace(self):
        """Skip spaces (but not newlines)."""
        while True:
            b = self.stream.read(1)
            if not b:
                return
            if b[0] != ord(' '):
                # Put it back
                if hasattr(self.stream, 'seek'):
                    self.stream.seek(-1, 1)
                else:
                    self._buffer = b + self._buffer
                return

    def read_token(self) -> bytes:
        """Read a whitespace-delimited token."""
        # Skip leading spaces
        result = bytearray()

        # Skip spaces
        while True:
            b = self.stream.read(1)
            if not b:
                if result:
                    return bytes(result)
                raise EOFError("Unexpected end of stream")
            if b[0] not in (ord(' '), ord('\n'), ord('\r'), ord('\t')):
                result.append(b[0])
                break

        # Read until whitespace
        while True:
            b = self.stream.read(1)
            if not b:
                break
            if b[0] in (ord(' '), ord('\n'), ord('\r'), ord('\t')):
                break
            result.append(b[0])

        return bytes(result)

    def read_int(self) -> int:
        """Read a space-separated integer."""
        token = self.read_token()
        return int(token.decode('ascii'))

    def read_long(self) -> int:
        """Read a space-separated long integer."""
        return self.read_int()

    def read_string(self) -> str:
        """Read a length-prefixed string."""
        length = self.read_int()
        # read_int()'s read_token() already consumed the trailing space
        # Read exactly length bytes
        data = self._read_bytes(length)
        return data.decode('utf-8', errors='replace')

    def read_hex_mpz(self) -> int:
        """Read a hexadecimal-encoded GMP integer."""
        token = self.read_token()
        s = token.decode('ascii')
        if s.startswith('-'):
            return -int(s[1:], 16)
        return int(s, 16)

    def read_rational(self) -> Union[int, Fraction]:
        """
        Read a rational number (Q coefficient).

        Returns int for integers, Fraction for true rationals.
        """
        subtype = self.read_int()

        if subtype == 4:
            # Small integer
            return self.read_int()
        elif subtype == 8:
            # Raw mpz integer (hex)
            return self.read_hex_mpz()
        elif subtype == 3:
            # mpz integer (decimal) - old format
            return self.read_int()
        elif subtype in (0, 1):
            # Fraction with decimal mpz
            num = self.read_int()
            den = self.read_int()
            return Fraction(num, den)
        elif subtype in (5, 6):
            # Raw fraction (hex)
            num = self.read_hex_mpz()
            den = self.read_hex_mpz()
            return Fraction(num, den)
        else:
            raise ValueError(f"Unknown rational subtype: {subtype}")

    def read_number(self, ring: Optional[SSIRing] = None) -> Any:
        """
        Read a coefficient based on ring type.
        """
        ring = ring or self.current_ring

        if ring is None:
            # Assume Q
            return self.read_rational()

        if ring.characteristic == 0:
            # Q or extension
            if ring.extension_type == 'transcendental':
                # Fraction of polynomials
                num = self.read_poly_body(ring.extension_ring)
                den = self.read_poly_body(ring.extension_ring)
                return (num, den)  # Return as tuple
            elif ring.extension_type == 'algebraic':
                # Single polynomial
                return self.read_poly_body(ring.extension_ring)
            else:
                return self.read_rational()
        else:
            # Z/p - just an integer
            return self.read_int()

    def read_ring(self) -> Optional[SSIRing]:
        """Read a ring definition."""
        ch = self.read_int()

        # Handle special characteristic codes
        if ch == -4:
            # NULL ring
            return None

        if ch == -5:
            # Cached ring reference
            index = self.read_int()
            return self.ring_cache.get(index)

        new_cache_index = None
        if ch == -6:
            # New cached ring
            new_cache_index = self.read_int()
            ch = self.read_int()

        if ch == -5:
            # Reference after -6
            index = self.read_int()
            ring = self.ring_cache.get(index)
            if new_cache_index is not None and ring:
                self.ring_cache[new_cache_index] = ring
            return ring

        if ch == -4:
            return None

        # Read number of variables
        n_vars = self.read_int()

        # Handle named coefficient field
        cf_name = None
        if ch == -3:
            cf_name = self.read_string()

        # Read variable names
        variables = []
        for _ in range(n_vars):
            variables.append(self.read_string())

        # Read ordering blocks
        n_blocks = self.read_int()
        ordering_blocks = []
        for _ in range(n_blocks):
            ord_type = self.read_int()
            block0 = self.read_int()
            block1 = self.read_int()

            # Read weights for weighted orderings
            # Enum values from Singular's ring.h:
            #   1=a, 5=M, 12=wp, 13=Wp, 18=ws, 19=Ws, 22=aa
            weights = []
            if ord_type in (1, 12, 13, 18, 19, 22):  # a, wp, Wp, ws, Ws, aa
                n_weights = block1 - block0 + 1
                weights = [self.read_int() for _ in range(n_weights)]
            elif ord_type == 5:  # M (matrix ordering)
                n_weights = block1 - block0 + 1
                weights = [self.read_int() for _ in range(n_weights * n_weights)]

            ordering_blocks.append((ord_type, block0, block1, weights))

        # Read extension ring for transcendental/algebraic extensions
        extension_ring = None
        extension_type = None
        if ch == -1:
            extension_type = 'transcendental'
            extension_ring = self.read_ring()
        elif ch == -2:
            extension_type = 'algebraic'
            extension_ring = self.read_ring()

        # Read quotient ideal
        quotient_ideal = self.read_ideal_body(n_vars)
        if quotient_ideal and len(quotient_ideal.generators) == 0:
            quotient_ideal = None

        ring = SSIRing(
            characteristic=ch if ch >= 0 else 0,
            variables=variables,
            ordering_blocks=ordering_blocks,
            quotient_ideal=quotient_ideal,
            extension_ring=extension_ring,
            extension_type=extension_type
        )

        if new_cache_index is not None:
            self.ring_cache[new_cache_index] = ring

        return ring

    def read_term(self, ring: Optional[SSIRing] = None) -> SSITerm:
        """Read a single polynomial term."""
        ring = ring or self.current_ring
        n_vars = len(ring.variables) if ring else 1

        coef = self.read_number(ring)
        comp = self.read_int()
        exponents = [self.read_int() for _ in range(n_vars)]

        return SSITerm(coefficient=coef, component=comp, exponents=exponents)

    def read_poly_body(self, ring: Optional[SSIRing] = None) -> SSIPoly:
        """Read polynomial data (without type code)."""
        ring = ring or self.current_ring
        n_terms = self.read_int()
        terms = [self.read_term(ring) for _ in range(n_terms)]
        return SSIPoly(terms=terms, ring=ring)

    def read_ideal_body(self, n_vars: Optional[int] = None) -> SSIIdeal:
        """Read ideal data (without type code)."""
        n_gens = self.read_int()
        generators = [self.read_poly_body() for _ in range(n_gens)]
        return SSIIdeal(generators=generators, ring=self.current_ring)

    def read_matrix_body(self) -> SSIMatrix:
        """Read matrix data (without type code)."""
        rows = self.read_int()
        cols = self.read_int()
        entries = [self.read_poly_body() for _ in range(rows * cols)]
        return SSIMatrix(rows=rows, cols=cols, entries=entries, ring=self.current_ring)

    def read_module_body(self) -> SSIModule:
        """Read module data (without type code)."""
        rank = self.read_int()
        n_gens = self.read_int()
        generators = [self.read_poly_body() for _ in range(n_gens)]
        return SSIModule(rank=rank, generators=generators, ring=self.current_ring)

    def read_command(self) -> SSICommand:
        """Read a command (unevaluated expression)."""
        argc = self.read_int()
        opcode = self.read_int()
        arguments = [self.read_value() for _ in range(argc)]
        return SSICommand(opcode=opcode, arguments=arguments)

    def read_list(self) -> list:
        """Read a heterogeneous list."""
        length = self.read_int()
        return [self.read_value() for _ in range(length)]

    def read_intvec(self) -> List[int]:
        """Read an integer vector."""
        length = self.read_int()
        return [self.read_int() for _ in range(length)]

    def read_intmat(self) -> List[List[int]]:
        """Read an integer matrix."""
        rows = self.read_int()
        cols = self.read_int()
        data = [self.read_int() for _ in range(rows * cols)]
        # Convert to 2D list
        return [data[i*cols:(i+1)*cols] for i in range(rows)]

    def read_bigintmat(self) -> List[List[int]]:
        """Read a big integer matrix."""
        rows = self.read_int()
        cols = self.read_int()
        data = [self.read_rational() for _ in range(rows * cols)]
        return [data[i*cols:(i+1)*cols] for i in range(rows)]

    def read_bigintvec(self) -> List[int]:
        """Read a big integer vector."""
        length = self.read_int()
        return [self.read_rational() for _ in range(length)]

    def read_ring_properties(self):
        """Read and apply ring properties (type 23)."""
        subtype = self.read_int()

        if subtype == 0:
            # Custom bitmask
            log2_bitmask = self.read_int()
            if self.current_ring:
                self.current_ring.bitmask = (1 << log2_bitmask) - 1
        elif subtype == 1:
            # Letterplace ring
            log2_bitmask = self.read_int()
            is_lp = self.read_int()
            if self.current_ring:
                self.current_ring.bitmask = (1 << log2_bitmask) - 1
                self.current_ring.is_letterplace = bool(is_lp)
        elif subtype == 2:
            # Plural (NC) ring
            C = self.read_matrix_body()
            D = self.read_matrix_body()
            if self.current_ring:
                self.current_ring.is_plural = True
                self.current_ring.nc_C = C
                self.current_ring.nc_D = D

    def read_attributes(self) -> SSIAttributed:
        """Read attributed value (type 21)."""
        flags = self.read_int()
        n_attrs = self.read_int()

        attributes = {}
        # Note: attribute reading is incomplete in original code
        # For now, just skip them
        for _ in range(n_attrs):
            pass  # Would need to read attr name and value

        value = self.read_value()
        return SSIAttributed(flags=flags, attributes=attributes, value=value)

    def read_blackbox(self) -> SSIBlackbox:
        """Read a blackbox (custom type)."""
        # First read the typename
        typename = self.read_value()  # Usually a string
        # The rest depends on the blackbox type
        # We can't properly deserialize without knowing the type
        return SSIBlackbox(typename=str(typename), data=None)

    def read_value(self) -> Any:
        """
        Read a single SSI value.

        Returns the appropriate Python object for the SSI type.
        """
        type_code = self.read_int()

        if type_code == self.TYPE_INT:
            return self.read_int()

        elif type_code == self.TYPE_STRING:
            return self.read_string()

        elif type_code == self.TYPE_NUMBER:
            return self.read_number()

        elif type_code == self.TYPE_BIGINT:
            return self.read_rational()

        elif type_code == self.TYPE_RING:
            ring = self.read_ring()
            self.current_ring = ring
            return ring

        elif type_code == self.TYPE_POLY:
            return self.read_poly_body()

        elif type_code == self.TYPE_IDEAL:
            return self.read_ideal_body()

        elif type_code == self.TYPE_MATRIX:
            return self.read_matrix_body()

        elif type_code == self.TYPE_VECTOR:
            return self.read_poly_body()  # Same as poly but component > 0

        elif type_code == self.TYPE_MODULE:
            return self.read_module_body()

        elif type_code == self.TYPE_COMMAND:
            return self.read_command()

        elif type_code == self.TYPE_DEF:
            return SSIDef(name=self.read_string())

        elif type_code == self.TYPE_PROC:
            return SSIProc(body=self.read_string())

        elif type_code == self.TYPE_LIST:
            return self.read_list()

        elif type_code == self.TYPE_RING_DATA:
            # Ring context switch
            ring = self.read_ring()
            self.current_ring = ring
            # Now read the actual value
            return self.read_value()

        elif type_code == self.TYPE_NONE:
            return SSINone()

        elif type_code == self.TYPE_INTVEC:
            return self.read_intvec()

        elif type_code == self.TYPE_INTMAT:
            return self.read_intmat()

        elif type_code == self.TYPE_BIGINTMAT:
            return self.read_bigintmat()

        elif type_code == self.TYPE_BLACKBOX:
            return self.read_blackbox()

        elif type_code == self.TYPE_ATTRIB:
            return self.read_attributes()

        elif type_code == self.TYPE_SMATRIX:
            return self.read_module_body()  # Same format as module

        elif type_code == self.TYPE_RINGPROP:
            self.read_ring_properties()
            return self.read_value()  # Continue reading next value

        elif type_code == self.TYPE_BIGINTVEC:
            return self.read_bigintvec()

        elif type_code == self.TYPE_VERSION:
            version = self.read_int()
            max_tok = self.read_int()
            opt1 = self.read_int()
            opt2 = self.read_int()
            return SSIVersion(version=version, max_tok=max_tok, opt1=opt1, opt2=opt2)

        elif type_code == self.TYPE_QUIT:
            return SSIQuit()

        elif type_code == 0:
            # EOF or error
            return SSINone()

        else:
            raise ValueError(f"Unknown type code: {type_code}")

    def read_message(self) -> Any:
        """
        Read a complete SSI message (value + newline).

        This is the main entry point for reading SSI data.
        """
        value = self.read_value()
        # Consume trailing whitespace/newline
        try:
            while True:
                b = self.stream.read(1)
                if not b or b[0] not in (ord(' '), ord('\n'), ord('\r'), ord('\t')):
                    if b and hasattr(self.stream, 'seek'):
                        self.stream.seek(-1, 1)
                    break
        except:
            pass
        return value


# =============================================================================
# Writer Implementation
# =============================================================================

class SSIWriter:
    """
    Writes SSI-encoded data to a binary stream.
    """

    def __init__(self, stream: BinaryIO):
        self.stream = stream
        self.current_ring: Optional[SSIRing] = None
        self.ring_cache: Dict[int, SSIRing] = {}
        self._next_ring_id = 0

    def write(self, data: bytes):
        """Write raw bytes."""
        self.stream.write(data)

    def write_int(self, n: int):
        """Write a space-separated integer."""
        self.write(f"{n} ".encode('ascii'))

    def write_string(self, s: str):
        """Write a length-prefixed string."""
        encoded = s.encode('utf-8')
        self.write(f"{len(encoded)} ".encode('ascii'))
        self.write(encoded)
        self.write(b' ')

    def write_rational(self, n: Union[int, Fraction]):
        """Write a rational number."""
        if isinstance(n, int):
            if -2**27 <= n < 2**27:
                # Small integer
                self.write(f"4 {n} ".encode('ascii'))
            else:
                # Large integer - use hex
                if n < 0:
                    self.write(f"8 -{abs(n):x} ".encode('ascii'))
                else:
                    self.write(f"8 {n:x} ".encode('ascii'))
        elif isinstance(n, Fraction):
            num, den = n.numerator, n.denominator
            if num < 0:
                num_hex = f"-{abs(num):x}"
            else:
                num_hex = f"{num:x}"
            den_hex = f"{den:x}"
            self.write(f"6 {num_hex} {den_hex} ".encode('ascii'))
        else:
            # Try to convert to int
            self.write_rational(int(n))

    def write_poly(self, poly: SSIPoly):
        """Write a polynomial."""
        self.write_int(6)  # TYPE_POLY
        self.write_int(len(poly.terms))
        for term in poly.terms:
            self.write_rational(term.coefficient)
            self.write_int(term.component)
            for exp in term.exponents:
                self.write_int(exp)

    def write_ideal(self, ideal: SSIIdeal):
        """Write an ideal."""
        self.write_int(7)  # TYPE_IDEAL
        self.write_int(len(ideal.generators))
        for gen in ideal.generators:
            self.write_int(len(gen.terms))
            for term in gen.terms:
                self.write_rational(term.coefficient)
                self.write_int(term.component)
                for exp in term.exponents:
                    self.write_int(exp)

    def write_list(self, lst: list):
        """Write a heterogeneous list."""
        self.write_int(14)  # TYPE_LIST
        self.write_int(len(lst))
        for item in lst:
            self.write_value(item)

    def write_intvec(self, vec: List[int]):
        """Write an integer vector."""
        self.write_int(17)  # TYPE_INTVEC
        self.write_int(len(vec))
        for v in vec:
            self.write_int(v)

    def write_value(self, value: Any):
        """Write any supported value."""
        if isinstance(value, int):
            self.write_int(1)  # TYPE_INT
            self.write_int(value)
        elif isinstance(value, str):
            self.write_int(2)  # TYPE_STRING
            self.write_string(value)
        elif isinstance(value, Fraction):
            self.write_int(4)  # TYPE_BIGINT
            self.write_rational(value)
        elif isinstance(value, SSIPoly):
            self.write_poly(value)
        elif isinstance(value, SSIIdeal):
            self.write_ideal(value)
        elif isinstance(value, list):
            if all(isinstance(x, int) for x in value):
                self.write_intvec(value)
            else:
                self.write_list(value)
        elif isinstance(value, SSINone):
            self.write_int(16)  # TYPE_NONE
        elif value is None:
            self.write_int(16)  # TYPE_NONE
        else:
            raise TypeError(f"Cannot write type: {type(value)}")

    def write_message(self, value: Any):
        """Write a complete message with newline."""
        self.write_value(value)
        self.write(b'\n')
        self.stream.flush()

    def write_quit(self):
        """Write quit message."""
        self.write_int(99)
        self.write(b'\n')
        self.stream.flush()


# =============================================================================
# Connection Helper
# =============================================================================

class SSIConnection:
    """
    Manages an SSI TCP connection to Singular.

    Example:
        conn = SSIConnection('localhost', 12345)
        version = conn.handshake
        conn.write(42)
        result = conn.read()
        conn.close()
    """

    def __init__(self, host: str, port: int, timeout: float = 30.0):
        """
        Connect to an SSI server.

        Args:
            host: Server hostname
            port: Server port
            timeout: Connection timeout in seconds
        """
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.settimeout(timeout)
        self.socket.connect((host, port))

        # Create buffered stream
        self._rfile = self.socket.makefile('rb')
        self._wfile = self.socket.makefile('wb')

        self.reader = SSIReader(self._rfile)
        self.writer = SSIWriter(self._wfile)

        # Read handshake
        self.handshake: Optional[SSIVersion] = None
        try:
            first_value = self.reader.read_value()
            if isinstance(first_value, SSIVersion):
                self.handshake = first_value
        except:
            pass

    def read(self) -> Any:
        """Read a value from the connection."""
        return self.reader.read_message()

    def write(self, value: Any):
        """Write a value to the connection."""
        self.writer.write_message(value)

    def close(self):
        """Close the connection."""
        try:
            self.writer.write_quit()
        except:
            pass
        try:
            self._rfile.close()
            self._wfile.close()
            self.socket.close()
        except:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


# =============================================================================
# Utility Functions
# =============================================================================

def read_ssi_file(filename: str) -> List[Any]:
    """
    Read all values from an SSI file.

    Args:
        filename: Path to the SSI file

    Returns:
        List of values read from the file
    """
    values = []
    with open(filename, 'rb') as f:
        reader = SSIReader(f)
        while True:
            try:
                value = reader.read_message()
                if isinstance(value, SSIQuit):
                    break
                values.append(value)
            except EOFError:
                break
    return values


def write_ssi_file(filename: str, values: List[Any]):
    """
    Write values to an SSI file.

    Args:
        filename: Path to the SSI file
        values: List of values to write
    """
    with open(filename, 'wb') as f:
        writer = SSIWriter(f)
        for value in values:
            writer.write_message(value)


# =============================================================================
# Main / Demo
# =============================================================================

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        # Read and dump SSI file
        filename = sys.argv[1]
        print(f"Reading SSI file: {filename}")
        try:
            values = read_ssi_file(filename)
            for i, v in enumerate(values):
                print(f"[{i}] {type(v).__name__}: {v}")
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("Usage: python ssi_reader.py <filename.ssi>")
        print()
        print("Or use as a library:")
        print("  from ssi_reader import SSIReader, SSIConnection")
        print()
        print("  # Read from file")
        print("  with open('data.ssi', 'rb') as f:")
        print("      reader = SSIReader(f)")
        print("      value = reader.read_message()")
        print()
        print("  # Connect to Singular")
        print("  conn = SSIConnection('localhost', 12345)")
        print("  value = conn.read()")
