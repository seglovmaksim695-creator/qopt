# SPDX-License-Identifier: Apache-2.0
"""qoptct.verify.dyadic -- exact verification of circuits with rotations pi*p/2^k (QFT, Draper adders, QFT
multipliers) by a sum-over-paths whose phases are integers modulo 2^K (Amy, arXiv:1805.06908v2, generalised from Z8).

    from qoptct.verify.dyadic import DCircuit, Pi, prove_dyadic, prove_equal, prove_spec, QFTSpec, adder_spec

    c = DCircuit(3); c.add('h', (2,)); c.add('cp', (1, 2), 'pi/2'); ...
    prove_equal(c, reference)                    # circuit = exp(i g) reference ?
    prove_spec(c, QFTSpec([0, 1, 2]))            # circuit = exp(i g) QFT ?
    prove_spec(adder, adder_spec(a, b))          # circuit = (b += a mod 2^n) ?  (assumption: fourier-shift lemma)
    prove_dyadic(c, reference_or_spec)           # one entry point for the bridges
    -> {'verdict': 'PROVED'|'REFUTED'|'UNKNOWN'|'ERROR'|'INTERNAL', 'method', 'global_phase', 'counterexample', ...}

The existing Clifford+T prover qoptct.verify.pathsum (phases in Z8) is not used or changed.  See docs/DYADIC.md.
"""
from .circuit import DCircuit, Pi, NotDyadic, Unsupported, parse_qasm, from_qc, to_angle
from .spec import (QFTSpec, ArithSpec, PermutationSpec, ClassicalSpec, adder_spec, const_adder_spec,
                   multiplier_spec)
from .api import prove_dyadic, prove_equal, prove_spec, prove_qasm, VERDICTS

__all__ = ['DCircuit', 'Pi', 'NotDyadic', 'Unsupported', 'parse_qasm', 'from_qc', 'to_angle', 'QFTSpec',
           'ArithSpec', 'PermutationSpec', 'ClassicalSpec', 'adder_spec', 'const_adder_spec', 'multiplier_spec',
           'prove_equal', 'prove_spec', 'prove_qasm', 'prove_dyadic', 'VERDICTS']
