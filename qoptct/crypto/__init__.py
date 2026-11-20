# SPDX-License-Identifier: Apache-2.0
"""qoptct.crypto -- circuit generators for cryptographic arithmetic on logical ANDs (pure functions, no optimizer
state). Model: Gidney 2018 (arXiv:1709.06648), AND = 4 T (+ 3 CX), erased by X-basis measurement with a Clifford
correction (0 T). Every generator checks its circuit at build time or offers an exact check.

Modules:
    kaliski  kaliski_inverse(p, n=None, reduction='generic', clean=False, legacy=False) -- Kaliski modular
             inversion, 12n^2 - 8n - 3 ANDs forward (legacy=True: 12n^2 + n; counted for every n tested, not
             proved for all n)
    mulsq    multiply(n, karatsuba_levels=0, ...), square(n, ...) -- integer multiplication / squaring on
             conditionally inverted additions, optional Karatsuba levels with garbage
    fold     fold.reduce_special(p, ...) -- reduction of a 2n-bit product modulo special primes (p = 2^n - c,
             NIST/Solinas P-256, Mersenne P-521); addition of a classical constant and its controlled version
    andcirc  small engine used by fold: wires, AND / measurement-based erase, Gidney adders (register or constant
             addend), exact bit-sliced simulation, OpenQASM 3 output (qoptct dialect)
    lookup   lookup.narrow_lookup(table, c, w, ...) -- narrow table lookups z ^= T[a] (bilinear split, QROAM s=1,
             ANF/Davio), measurement-based unlookup, lookup-accumulate acc += T[a]; fanout='gray' option for Davio
    gf2m     gf2m_multiply / gf2m_square / gf2m_inverse / gf2m_divide -- GF(2^m) arithmetic (binary fields of the
             NIST/SEC curves): CRT bilinear multiplier, in-place squaring, Itoh-Tsujii inversion and division
    ecpoint  ec_point_add / ec_point_add_window -- exact complete addition of a classical point to a quantum point on a
             prime-field curve (core in the private subpackage _ecp; the default construction is proved at n = 256 by
             qoptct.verify.compose.verify_ecadd); approximate schemes are cheaper
    residue  fold_accumulator(q, clen, nbits) -- accumulator compression modulo q (loop 2 of Gidney 2025),
             accumulator_bits
    modadd   mod_add / mod_sub / mod_double / compare_geq ... -- modular adders (4n - 2 ANDs)
    gf2ec    gf2_point_add(m, ...) -- point addition on binary curves over GF(2^m), exceptional cases not handled

The names exported here come from kaliski, mulsq and gf2m; fold, lookup, ecpoint, residue, modadd and gf2ec are used as
modules (``from qoptct.crypto import fold, lookup, ecpoint``), since lookup has its own ``reference_qc``.

    from qoptct.crypto import multiply, square
    c = multiply(256, karatsuba_levels=3)      # (x, y, 0) -> (x, y, x*y) + garbage
    c.counts()['AND']                          # 37 826 Toffoli (= logical AND, 4 T each); legacy=True: 37 880

See each module's docstring for the constructions, wire layouts and the comparison with published circuits.
"""
from ._errors import CheckError
from .kaliski import kaliski_inverse, kaliski_and_count, KaliskiCircuit
from .mulsq import multiply, square, MulSqCircuit, reference_qc, prove
from .gf2m import (default_modulus, gf2m_bilinear, gf2m_multiply, gf2m_mult_count, gf2m_square, gf2m_inverse,
                   gf2m_divide, gf2m_inverse_count, gf2m_divide_count, gf2m_itoh_gate_check, verify_multiplier)

__all__ = ['CheckError',
           'kaliski_inverse', 'kaliski_and_count', 'KaliskiCircuit',
           'multiply', 'square', 'MulSqCircuit', 'reference_qc', 'prove',
           'default_modulus', 'gf2m_bilinear', 'gf2m_multiply', 'gf2m_mult_count', 'gf2m_square', 'gf2m_inverse',
           'gf2m_divide', 'gf2m_inverse_count', 'gf2m_divide_count', 'gf2m_itoh_gate_check', 'verify_multiplier']
