# SPDX-License-Identifier: Apache-2.0
"""Qiskit arithmetic circuits without rotations: builders with their specifications.

Every family function returns ``Case(circuit, spec, clean, name, n)``; ``prove_case(case)`` runs check_qiskit.
Register layouts follow Qiskit (LSB first).  ``controlled(case)`` builds the controlled version from Qiskit's own
``.control(1)`` of every primitive (x->cx, cx->ccx, ccx->C3X, mcx->mcx, swap->cswap): Qiskit's generic
``Gate.control`` of a whole adder inserts fractional multi-controlled phases (see docs/QISKIT.md, list for Q2).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library import (AndGate, BitwiseXorGate, CDKMRippleCarryAdder, CCXGate, CXGate, FullAdderGate,
                                    HalfAdderGate, IntegerComparator, IntegerComparatorGate, MCXGate,
                                    ModularAdderGate, OrGate, SwapGate, VBERippleCarryAdder, XGate)

from .check import Spec, check_qiskit
from .translate import translate


@dataclass
class Case:
    name: str
    n: int
    circuit: QuantumCircuit
    spec: Spec
    clean: tuple
    note: str = ''


def _m(v, n):
    return v & ((1 << n) - 1)


def _bit(S, x):
    from qoptct.verify.ksym import SymInt
    return SymInt(S.B, [x])


# ------------------------------------------------------------------------------------------------ adders
def adder_spec(n, cin=False, cout=False, name='adder'):
    """b := (cin + a + b) mod 2^n;  cout ^= carry  (cout may be dirty)."""
    def classical(v):
        s = v['a'] + v['b'] + (v['cin'] if cin else 0)
        out = {'b': _m(s, n)}
        if cout:
            out['cout'] = v['cout'] ^ (s >> n)
        return out

    def symbolic(S, R):
        s = R['a'].add(R['b'], cin=R['cin'].bits[0]) if cin else R['a'].add(R['b'])
        out = {'b': s.resize(n)}
        if cout:
            out['cout'] = _bit(S, S.B.XOR(R['cout'].bits[0], s.bits[n]))
        return out
    return Spec(name, classical, symbolic, 'b = (a + b%s) mod 2^%d%s' % (' + cin' if cin else '', n,
                                                                     ', cout ^= carry' if cout else ''))


def _kind_flags(kind):
    return dict(half=(False, True), full=(True, True), fixed=(False, False))[kind]


def cdkm(n, kind='half'):
    cin, cout = _kind_flags(kind)
    return Case('CDKMRippleCarryAdder-' + kind, n, CDKMRippleCarryAdder(n, kind=kind),
                adder_spec(n, cin, cout, 'cdkm-' + kind), ('help',))


def vbe(n, kind='half'):
    cin, cout = _kind_flags(kind)
    return Case('VBERippleCarryAdder-' + kind, n, VBERippleCarryAdder(n, kind=kind),
                adder_spec(n, cin, cout, 'vbe-' + kind), ('helper',))


def _wrap(gate, layout):
    regs = [QuantumRegister(size, name) for name, size in layout]
    qc = QuantumCircuit(*regs, name=gate.name)
    qc.append(gate, [q for r in regs for q in r])
    return qc


def half_adder_gate(n):
    return Case('HalfAdderGate', n, _wrap(HalfAdderGate(n), [('a', n), ('b', n), ('cout', 1)]),
                adder_spec(n, False, True, 'half-adder-gate'), ('cout',),
                'cout is a clean output (|a>|b>|0> -> |a>|a+b>)')


def full_adder_gate(n):
    return Case('FullAdderGate', n, _wrap(FullAdderGate(n), [('cin', 1), ('a', n), ('b', n), ('cout', 1)]),
                adder_spec(n, True, True, 'full-adder-gate'), ('cout',))


def modular_adder_gate(n):
    return Case('ModularAdderGate', n, _wrap(ModularAdderGate(n), [('a', n), ('b', n)]),
                adder_spec(n, False, False, 'modular-adder-gate'), ())


# ------------------------------------------------------------------------------------------------ comparators
def comparator_spec(n, value, geq=True):
    def classical(v):
        r = (v['state'] >= value) if geq else (v['state'] < value)
        return {'compare': v['compare'] ^ int(r)}

    def symbolic(S, R):
        lt = R['state'].lt(value)
        r = lt ^ 1 if geq else lt
        return {'compare': _bit(S, S.B.XOR(R['compare'].bits[0], r))}
    return Spec('comparator', classical, symbolic, 'compare ^= [state %s %d]' % ('>=' if geq else '<', value), msb_first=True)


def integer_comparator(n, value, geq=True):
    qc = IntegerComparator(n, value, geq=geq)
    anc = tuple(r.name for r in qc.qregs if r.name not in ('state', 'compare'))
    return Case('IntegerComparator', n, qc, comparator_spec(n, value, geq), anc)


def integer_comparator_gate(n, value, geq=True):
    return Case('IntegerComparatorGate', n, _wrap(IntegerComparatorGate(n, value, geq=geq),
                                                  [('state', n), ('compare', 1)]),
                comparator_spec(n, value, geq), ())


# ------------------------------------------------------------------------------------------------ logic gates
def and_gate(n):
    def classical(v):
        return {'out': v['out'] ^ int(v['x'] == (1 << n) - 1)}

    def symbolic(S, R):
        r = R['x'].eq((1 << n) - 1)
        return {'out': _bit(S, S.B.XOR(R['out'].bits[0], r))}
    return Case('AndGate', n, _wrap(AndGate(n), [('x', n), ('out', 1)]), Spec('and', classical, symbolic), ())


def or_gate(n):
    def classical(v):
        return {'out': v['out'] ^ int(v['x'] != 0)}

    def symbolic(S, R):
        r = R['x'].eq(0) ^ 1
        return {'out': _bit(S, S.B.XOR(R['out'].bits[0], r))}
    return Case('OrGate', n, _wrap(OrGate(n), [('x', n), ('out', 1)]), Spec('or', classical, symbolic), ())


def bitwise_xor(n, amount):
    def classical(v):
        return {'x': v['x'] ^ _m(amount, n)}

    def symbolic(S, R):
        return {'x': R['x'] ^ S.const(_m(amount, n), n)}
    return Case('BitwiseXorGate', n, _wrap(BitwiseXorGate(n, amount), [('x', n)]), Spec('xor', classical, symbolic), ())


# ------------------------------------------------------------------------------------------------ controlled
def _leaf_gate(op):
    name, k = op[0], len(op) - 1
    return {'x': XGate, 'cx': CXGate, 'ccx': CCXGate, 'swap': SwapGate}.get(name, lambda: MCXGate(k - 1))()


def controlled(case):
    """Controlled version (control register 'ctrl' first), built from Qiskit's .control(1) of every primitive."""
    tr = translate(case.circuit)
    if not tr.is_classical():
        raise ValueError('only classical circuits are controlled here')
    regs = [QuantumRegister(1, 'ctrl')] + list(case.circuit.qregs)
    qc = QuantumCircuit(*regs, name='c' + case.circuit.name)
    pos = {i: qc.find_bit(q).index for i, q in enumerate(case.circuit.qubits)}
    # wire i of the original is wire i+1 in the new circuit (ctrl is register 0, size 1)
    for op in tr.ops:
        g = _leaf_gate(op).control(1)
        qc.append(g, [0] + [op_w + 1 for op_w in op[1:]])
    base = case.spec

    def classical(v):
        c = v['ctrl']
        out = base.classical({k: x for k, x in v.items() if k != 'ctrl'})
        return out if c else {}

    def symbolic(S, R):
        res = base.symbolic(S, {k: x for k, x in R.items() if k != 'ctrl'})
        return {k: S.mux(R['ctrl'].bits[0], v, R[k]) for k, v in res.items()}
    spec = Spec('controlled-' + base.name, classical, symbolic, 'ctrl ? ' + base.description + ' : identity',
                msb_first=base.msb_first)
    return Case('c' + case.name, case.n, qc, spec, case.clean, 'controlled by .control(1) of each primitive')


def prove_case(case, **kw):
    clean = tuple(n for n in case.clean if n in {r.name for r in case.circuit.qregs})
    r = check_qiskit(case.circuit, case.spec, clean=clean, **kw)
    r['case'] = case.name
    r['n'] = case.n
    return r


FAMILIES = {
    'cdkm-half': lambda n: cdkm(n, 'half'), 'cdkm-full': lambda n: cdkm(n, 'full'),
    'cdkm-fixed': lambda n: cdkm(n, 'fixed'),
    'vbe-half': lambda n: vbe(n, 'half'), 'vbe-full': lambda n: vbe(n, 'full'), 'vbe-fixed': lambda n: vbe(n, 'fixed'),
    'half-adder-gate': half_adder_gate, 'full-adder-gate': full_adder_gate, 'modular-adder-gate': modular_adder_gate,
    'comparator-geq': lambda n: integer_comparator(n, (1 << n) // 3 + 1, True),
    'comparator-lt': lambda n: integer_comparator(n, (1 << n) // 3 + 1, False),
    'comparator-gate': lambda n: integer_comparator_gate(n, (1 << n) // 3 + 1, True),
    'and': and_gate, 'or': or_gate, 'xor-const': lambda n: bitwise_xor(n, (1 << n) // 5 * 3 + 1),
}
for _k in ('cdkm-half', 'cdkm-full', 'cdkm-fixed', 'vbe-half', 'vbe-full', 'vbe-fixed'):
    FAMILIES['c-' + _k] = (lambda k: lambda n: controlled(FAMILIES[k](n)))(_k)
