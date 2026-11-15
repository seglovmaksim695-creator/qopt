# SPDX-License-Identifier: Apache-2.0
"""qoptct.crypto -- generators of arithmetic circuits for cryptanalytic resource estimates (logical-AND model).

    from qoptct.crypto import multiply, square
    c = multiply(256, karatsuba_levels=3)      # (x, y, 0) -> (x, y, x*y) + garbage
    c.counts()['AND']                          # 37 880 Toffoli (= logical AND, 4 T each)

See qoptct.crypto.mulsq for the constructions, the wire layout and the comparison with published circuits.
"""
from .mulsq import multiply, square, MulSqCircuit, CheckError, reference_qc, prove

__all__ = ['multiply', 'square', 'MulSqCircuit', 'CheckError', 'reference_qc', 'prove']
