# SPDX-License-Identifier: Apache-2.0
"""qoptct.crypto -- circuit generators for cryptographic arithmetic on logical ANDs (pure functions, no optimizer
state). Model: Gidney 2018 (arXiv:1709.06648), AND = 4 T (+ 3 CX), erased by X-basis measurement with a Clifford
correction (0 T). Every generator checks its circuit at build time or offers an exact check.

Modules:
    kaliski  kaliski_inverse(p, n=None, reduction='generic', clean=False) -- Kaliski modular inversion,
             12n^2 + n ANDs forward (counted for every n tested, not proved for all n)
    mulsq    multiply(n, karatsuba_levels=0, ...), square(n, ...) -- integer multiplication / squaring on
             conditionally inverted additions, optional Karatsuba levels with garbage
    fold     fold.reduce_special(p, ...) -- reduction of a 2n-bit product modulo special primes (p = 2^n - c,
             NIST/Solinas P-256, Mersenne P-521); addition of a classical constant and its controlled version
    andcirc  small engine used by fold: wires, AND / measurement-based erase, Gidney adders (register or constant
             addend), exact bit-sliced simulation, OpenQASM 3 output (qoptct dialect)

The names exported here come from kaliski and mulsq; fold is used as a module
(``from qoptct.crypto import fold``).

    from qoptct.crypto import multiply, square
    c = multiply(256, karatsuba_levels=3)      # (x, y, 0) -> (x, y, x*y) + garbage
    c.counts()['AND']                          # 37 880 Toffoli (= logical AND, 4 T each)

See each module's docstring for the constructions, wire layouts and the comparison with published circuits.
"""
from ._errors import CheckError
from .kaliski import kaliski_inverse, kaliski_and_count, KaliskiCircuit
from .mulsq import multiply, square, MulSqCircuit, reference_qc, prove

__all__ = ['CheckError',
           'kaliski_inverse', 'kaliski_and_count', 'KaliskiCircuit',
           'multiply', 'square', 'MulSqCircuit', 'reference_qc', 'prove']
