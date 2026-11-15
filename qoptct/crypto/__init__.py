# SPDX-License-Identifier: Apache-2.0
"""qoptct.crypto: generators of classical reversible arithmetic for cryptanalytic circuits, built on logical ANDs
(Gidney 2018: 4 T + 3 CX per AND, uncomputation by X-basis measurement with a Clifford correction, 0 T).

Modules:
    andcirc -- small engine: wires, AND / measurement-based erase, Gidney adders (register or constant addend),
               exact bit-sliced simulation, OpenQASM 3 output (qoptct dialect).
    fold    -- reduction of a 2n-bit product modulo special primes (p = 2^n - c, NIST/Solinas P-256,
               Mersenne P-521); addition of a classical constant and its controlled version.
"""
