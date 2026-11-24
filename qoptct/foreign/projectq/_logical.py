"""Recording -> qoptct logical program (x, cx, and, erase, and the Clifford+T atoms) and its OpenQASM 3 text.

Rules (every operation has exactly one rule; no recognition of circuit patterns, no reordering):

* ``x``                      -> x
* ``x`` with 1 control       -> cx
* ``x`` with k >= 2 controls -> a chain of k-1 logical ANDs on fresh scratch wires (one per link), cx onto the target,
                                then the AND erasures in reverse order (Gidney 2018).  The scratch wires are new wires
                                after the program's own wires; every one of them is |0> before and after the gate.
* ``z`` -> s, s;  ``cz a b`` -> h b, cx a b, h b;  ``h s sdg t tdg`` -> the same atom.

The AND gadget needs a clean target and its erasure needs ``target = p AND q``; both are decided EXACTLY by the provers
(exact2: for all inputs of weight <= 2; the exhaustive simulator: for every input).  ``alloc`` / ``free`` produce no
operation: a freed wire is simply required to end in |0> (checked by the provers like every other clean wire).
"""
from __future__ import annotations

import hashlib
import json

from ._record import Recording, Unsupported

__all__ = ['Logical', 'to_logical', 'emit_qasm', 'AND10']

AND10 = (('h', 2), ('t', 2), ('cx', 1, 2), ('tdg', 2), ('cx', 0, 2),
         ('t', 2), ('cx', 1, 2), ('tdg', 2), ('h', 2), ('sdg', 2))
_ONE = ('h', 's', 'sdg', 't', 'tdg', 'x')


class Logical:
    """A logical program with its layout: ``inputs`` / ``outputs`` / ``garbage`` map names to wires."""

    def __init__(self, ops, nq, inputs, outputs, garbage, free_wires, scratch, mode, source=None):
        self.ops = ops
        self.nq = nq
        self.inputs = inputs
        self.outputs = outputs
        self.garbage = garbage
        self.free_wires = free_wires      # wires with a free initial value, in the fixed order of Recording.free_input_wires
        self.scratch = scratch
        self.mode = mode
        self.source = source              # SHA-256 of the recording this program was made from

    def copy(self, ops=None):
        return Logical(list(self.ops if ops is None else ops), self.nq, {k: list(v) for k, v in self.inputs.items()},
                       {k: list(v) for k, v in self.outputs.items()}, {k: list(v) for k, v in self.garbage.items()},
                       list(self.free_wires), list(self.scratch), self.mode, self.source)

    def qasm(self):
        return emit_qasm(self.nq, self.ops)

    def sha256(self):
        blob = json.dumps([self.mode, self.nq, self.ops, sorted(self.inputs.items()), sorted(self.outputs.items()),
                           sorted(self.garbage.items())], separators=(',', ':'), default=list).encode()
        return hashlib.sha256(blob).hexdigest()

    def is_classical(self):
        return all(op[0] in ('x', 'cx', 'and', 'erase') for op in self.ops)


def to_logical(rec: Recording) -> Logical:
    ops = []
    nq = rec.n_wires
    scratch = []

    def need(k):
        nonlocal nq
        while len(scratch) < k:
            scratch.append(nq)
            nq += 1

    for op in rec.ops:
        k = op[0]
        if k in ('alloc', 'free'):
            continue
        if k == 'x':
            ctl, t = op[1], op[2]
            if len(ctl) == 0:
                ops.append(('x', t))
            elif len(ctl) == 1:
                ops.append(('cx', ctl[0], t))
            else:
                need(len(ctl) - 1)
                s = scratch[:len(ctl) - 1]
                ops.append(('and', ctl[0], ctl[1], s[0]))
                for i in range(2, len(ctl)):
                    ops.append(('and', s[i - 2], ctl[i], s[i - 1]))
                ops.append(('cx', s[-1], t))
                for i in range(len(ctl) - 1, 1, -1):
                    ops.append(('erase', s[i - 1], (('cz', s[i - 2], ctl[i]),)))
                ops.append(('erase', s[0], (('cz', ctl[0], ctl[1]),)))
        elif k == 'z':
            ops += [('s', op[1]), ('s', op[1])]
        elif k == 'cz':
            ops += [('h', op[2]), ('cx', op[1], op[2]), ('h', op[2])]
        elif k in ('h', 's', 'sdg', 't', 'tdg'):
            ops.append((k, op[1]))
        else:
            raise Unsupported(f'no rule for {k}')
    return Logical(ops, nq, {k: list(v) for k, v in rec.inputs.items()}, {k: list(v) for k, v in rec.outputs.items()},
                   {k: list(v) for k, v in rec.garbage.items()}, rec.free_input_wires(), scratch, rec.mode, rec.sha256())


def emit_qasm(nq, ops):
    """Canonical OpenQASM 3 text (the same AND gadget and erasure as the Qualtran bridge)."""
    lines = ['OPENQASM 3;', 'include "stdgates.inc";', f'qubit[{nq}] q;',
             f'bit[{max(1, sum(op[0] == "erase" for op in ops))}] c;']
    k = 0

    def gate(g, wires):
        if len(set(wires)) != len(wires) or any(type(w) is not int or not 0 <= w < nq for w in wires):
            raise ValueError(f'Invalid wires {g} {wires}')
        return g + ' ' + ', '.join(f'q[{w}]' for w in wires) + ';'
    for op in ops:
        if op[0] in _ONE or op[0] == 'cx':
            lines.append(gate(op[0], op[1:]))
        elif op[0] == 'and':
            a, b, t = op[1:]
            gate('and', (a, b, t))
            for g, *w in AND10:
                lines.append(gate(g, [(a, b, t)[i] for i in w]))
        elif op[0] == 'erase':
            _, t, corr = op
            lines += [gate('h', (t,)), f'c[{k}] = measure q[{t}];']
            body = []
            for g, *w in corr:
                if g not in ('z', 'cz') or t in w:
                    raise ValueError('Unsupported erase correction')
                body.append(gate(g, w))
            body.append(gate('x', (t,)))
            lines.append(f'if (c[{k}]) ' + '{' + ' '.join(body) + '}')
            k += 1
        else:
            raise ValueError(f'Unknown operation {op}')
    return '\n'.join(lines) + '\n'
