"""ProjectQ library blocks that fit into NCT / Clifford+T (no rotations), each with an independent specification.

Every block is a pair (program, specification) for a size ``n``; the program uses the block exactly as a ProjectQ user does
(``AddQuantum | (a, b)``, ``with Control(eng, c): ...``), and ProjectQ's own decomposition rules turn it into
X / CNOT / Toffoli.  The specification is plain Python arithmetic, written here, not taken from ProjectQ.

Blocks with rotations (QFT, Draper adders: AddConstant, AddConstantModN, MultiplyByConstantModN, ...) are not here: they
are the business of the fractional-phase prover.  ``rotation_blocks()`` lists them and why.
"""
from __future__ import annotations

from ._classical import ClassicalSpec

__all__ = ['BLOCKS', 'block', 'rotation_blocks']


def _prog(alloc, body, ins, outs):
    def program(eng):
        regs = {name: (eng.allocate_qubit() if k == 1 else eng.allocate_qureg(k)) for name, k in alloc}
        body(eng, regs)
        return dict(inputs={n: regs[n] for n in ins}, outputs={n: regs[n] for n in outs})
    return program


def add_quantum(n):
    from projectq.libs.math import AddQuantum
    prog = _prog([('a', n), ('b', n)], lambda e, r: AddQuantum | (r['a'], r['b']), ['a', 'b'], ['a', 'b'])
    return prog, ClassicalSpec(lambda v: dict(a=v['a'], b=(v['a'] + v['b']) % (1 << n)), name=f'add_quantum_n{n}')


def add_quantum_carry(n):
    from projectq.libs.math import AddQuantum
    prog = _prog([('a', n), ('b', n), ('c', 1)], lambda e, r: AddQuantum | (r['a'], r['b'], r['c']), ['a', 'b', 'c'],
                 ['a', 'b', 'c'])
    return prog, ClassicalSpec(lambda v: dict(a=v['a'], b=(v['a'] + v['b']) % (1 << n),
                                              c=v['c'] ^ ((v['a'] + v['b']) >> n & 1)), name=f'add_quantum_carry_n{n}')


def add_quantum_controlled(n):
    from projectq.libs.math import AddQuantum
    from projectq.meta import Control

    def body(eng, r):
        with Control(eng, r['k']):
            AddQuantum | (r['a'], r['b'])
    prog = _prog([('k', 1), ('a', n), ('b', n)], body, ['k', 'a', 'b'], ['k', 'a', 'b'])
    return prog, ClassicalSpec(lambda v: dict(k=v['k'], a=v['a'], b=(v['a'] + v['b']) % (1 << n) if v['k'] else v['b']),
                               name=f'add_quantum_controlled_n{n}')


def subtract_quantum(n):
    from projectq.libs.math import SubtractQuantum
    prog = _prog([('a', n), ('b', n)], lambda e, r: SubtractQuantum | (r['a'], r['b']), ['a', 'b'], ['a', 'b'])
    return prog, ClassicalSpec(lambda v: dict(a=v['a'], b=(v['b'] - v['a']) % (1 << n)), name=f'subtract_quantum_n{n}')


def comparator_quantum(n):
    from projectq.libs.math import ComparatorQuantum
    prog = _prog([('a', n), ('b', n), ('c', 1)], lambda e, r: ComparatorQuantum | (r['a'], r['b'], r['c']), ['a', 'b', 'c'],
                 ['a', 'b', 'c'])
    return prog, ClassicalSpec(lambda v: dict(a=v['a'], b=v['b'], c=v['c'] ^ int(v['a'] > v['b'])), name=f'comparator_quantum_n{n}')


def multiply_quantum(n):
    from projectq.libs.math import MultiplyQuantum
    prog = _prog([('a', n), ('b', n), ('c', 2 * n + 1)], lambda e, r: MultiplyQuantum | (r['a'], r['b'], r['c']), ['a', 'b', 'c'],
                 ['a', 'b', 'c'])
    return prog, ClassicalSpec(lambda v: dict(a=v['a'], b=v['b'], c=(v['c'] + v['a'] * v['b']) % (1 << 2 * n + 1)),
                               promise=lambda v: v['c'] == 0, name=f'multiply_quantum_n{n}')


def divide_quantum(n):
    from projectq.libs.math import DivideQuantum
    prog = _prog([('a', n), ('r', n), ('c', n)], lambda e, r: DivideQuantum | (r['a'], r['r'], r['c']), ['a', 'r', 'c'],
                 ['a', 'r', 'c'])
    return prog, ClassicalSpec(lambda v: dict(a=v['a'] % v['c'], r=v['a'] // v['c'], c=v['c']),
                               promise=lambda v: v['r'] == 0 and v['c'] > 0 and v['a'] >= v['c'], name=f'divide_quantum_n{n}')


def multi_controlled_x(n):
    from projectq.meta import Control
    from projectq.ops import X

    def body(eng, r):
        with Control(eng, list(r['k'])):
            X | r['t']
    prog = _prog([('k', n), ('t', 1)], body, ['k', 't'], ['k', 't'])
    return prog, ClassicalSpec(lambda v: dict(k=v['k'], t=v['t'] ^ int(v['k'] == (1 << n) - 1)), name=f'mcx_n{n}')


def controlled_swap(n):
    from projectq.meta import Control
    from projectq.ops import Swap

    def body(eng, r):
        with Control(eng, r['k']):
            Swap | (r['a'], r['b'])
    prog = _prog([('k', 1), ('a', 1), ('b', 1)], body, ['k', 'a', 'b'], ['k', 'a', 'b'])
    return prog, ClassicalSpec(lambda v: dict(k=v['k'], a=v['b'] if v['k'] else v['a'], b=v['a'] if v['k'] else v['b']),
                               name='fredkin')


def permutation_oracle(n):
    """``projectq.libs.revkit.PermutationOracle`` (needs the external ``revkit`` package) for a fixed pseudo-random permutation
    of 0..2^n-1 (seed 2026): |x> -> |perm[x]>."""
    import random
    from projectq.libs.revkit import PermutationOracle
    perm = list(range(1 << n))
    random.Random(2026 + n).shuffle(perm)
    prog = _prog([('q', n)], lambda e, r: PermutationOracle(perm) | r['q'], ['q'], ['q'])
    return prog, ClassicalSpec(lambda v: dict(q=perm[v['q']]), name=f'revkit_permutation_n{n}')


BLOCKS = {
    'revkit_permutation': permutation_oracle,
    'add_quantum': add_quantum,
    'add_quantum_carry': add_quantum_carry,
    'add_quantum_controlled': add_quantum_controlled,
    'subtract_quantum': subtract_quantum,
    'comparator_quantum': comparator_quantum,
    'multiply_quantum': multiply_quantum,
    'divide_quantum': divide_quantum,
    'multi_controlled_x': multi_controlled_x,
    'controlled_swap': controlled_swap,
}


def block(name, n):
    """-> (program, ClassicalSpec) of the library block ``name`` at size ``n``."""
    if name not in BLOCKS:
        raise KeyError(f'unknown block {name!r}; known: {sorted(BLOCKS)}')
    return BLOCKS[name](n)


def rotation_blocks():
    """ProjectQ blocks that need rotations (Q2, the fractional-phase prover): name -> reason.  Found by running every block
    through ``record`` in both gate sets; see docs/PROJECTQ.md for the table."""
    return {
        'AddConstant': 'Draper adder: QFT and controlled phase rotations R(pi/2^k)',
        'SubConstant': 'inverse of AddConstant',
        'AddConstantModN': 'Beauregard adder built from AddConstant',
        'SubConstantModN': 'inverse of AddConstantModN',
        'MultiplyByConstantModN': 'Beauregard multiplier built from AddConstantModN',
        'QFT': 'quantum Fourier transform: Hadamards and controlled phase rotations',
        'Rx/Ry/Rz/R/Ph': 'arbitrary rotations; only exact multiples of pi/4 (Z8) are atoms of the translation',
        'PhaseOracle/PermutationOracle (libs.revkit)': 'need the external revkit package',
    }
