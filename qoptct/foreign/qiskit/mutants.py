# SPDX-License-Identifier: Apache-2.0
"""Canary mutants of a Qiskit circuit: the circuit is flattened by the strict translator and rebuilt with one defect.

Every mutant must be REFUTED (or at least never PROVED) against the spec of the original.  The set covers seams
(first / last / middle gates), wrong wires, a duplicated or reordered gate, an unremoved ancilla, a phase-only defect
and a defect on a single rare input (an all-ones multi-controlled flip).
"""
from __future__ import annotations

from qiskit import QuantumCircuit
from qiskit.circuit.library import (CCXGate, CXGate, CZGate, HGate, MCXGate, SGate, SdgGate, SwapGate, TGate,
                                    TdgGate, XGate, ZGate)

from .translate import translate

_GATES = {'x': XGate, 'cx': CXGate, 'ccx': CCXGate, 'swap': SwapGate, 'h': HGate, 's': SGate, 'sdg': SdgGate,
          't': TGate, 'tdg': TdgGate, 'z': ZGate, 'cz': CZGate}


def rebuild(circuit, ops):
    qc = QuantumCircuit(*circuit.qregs, name='mutant')
    for op in ops:
        g = _GATES[op[0]]() if op[0] in _GATES else MCXGate(len(op) - 2)
        qc.append(g, list(op[1:]))
    return qc


def make_mutants(circuit, rare_wires=None, count=None):
    """-> list of (label, QuantumCircuit).  rare_wires: wires of one input register (for the rare-input mutant)."""
    tr = translate(circuit)
    ops = list(tr.ops)
    n = len(ops)
    out = []

    def add(label, new):
        out.append((label, rebuild(circuit, new)))
    pos = sorted({0, 1, n // 8, n // 4, (3 * n) // 8, n // 2 - 1, n // 2, n // 2 + 1, (5 * n) // 8, (3 * n) // 4,
                  (7 * n) // 8, n - 2, n - 1})
    for i in pos:
        if 0 <= i < n:
            add('drop[%d]' % i, ops[:i] + ops[i + 1:])
    for i in (n // 2, n - 1):
        add('dup-move[%d]' % i, ops[:i] + [ops[i], ops[i]] + ops[i + 1:] if False else ops[:i] + ops[i + 1:] + [ops[i]])
    # swap two adjacent gates at the middle seam
    m = n // 2
    if ops[m] != ops[m + 1]:
        add('swap-adjacent[%d]' % m, ops[:m] + [ops[m + 1], ops[m]] + ops[m + 2:])
    # wrong control / wrong target on the first, middle and last multi-qubit gates
    multi = [i for i, o in enumerate(ops) if o[0] in ('cx', 'ccx', 'mcx')]
    for i in sorted({multi[0], multi[len(multi) // 4], multi[len(multi) // 2], multi[(3 * len(multi)) // 4], multi[-1]}):
        o = ops[i]
        wires = set(range(tr.nq))
        alt = next(w for w in sorted(wires - set(o[1:])))
        add('wrong-control[%d]' % i, ops[:i] + [(o[0], alt) + tuple(o[2:])] + ops[i + 1:])
        add('wrong-target[%d]' % i, ops[:i] + [(o[0],) + tuple(o[1:-1]) + (alt,)] + ops[i + 1:])
    j = multi[len(multi) // 2]
    o = ops[j]
    if o[0] == 'ccx':
        add('ccx-to-cx[%d]' % j, ops[:j] + [('cx', o[1], o[3])] + ops[j + 1:])
    # extra gates: bit flip, extra cnot, left-over on an ancilla
    last = max(w for w in range(tr.nq))
    add('extra-x-last', ops + [('x', last)])
    add('extra-cx-first-last', ops + [('cx', 0, last)])
    add('extra-x-middle', ops[:m] + [('x', 1)] + ops[m:])
    # phase-only defects
    add('phase-z-end', ops + [('z', 0)])
    add('phase-t-middle', ops[:m] + [('t', 0)] + ops[m:])
    add('phase-s-first', [('s', last)] + ops)
    # a defect on a single rare input: flips wire 0's register bit only when the whole register is all ones
    if rare_wires:
        rw = list(rare_wires)
        other = next(w for w in range(tr.nq) if w not in rw)
        add('rare-allones', ops + [('mcx',) + tuple(rw) + (other,)] if len(rw) >= 3 else ops + [('ccx', rw[0], rw[1], other)])
    if count:
        out = out[:count]
    return out
