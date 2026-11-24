# SPDX-License-Identifier: Apache-2.0
"""Strict translation of a Qiskit ``QuantumCircuit`` into a qoptct circuit (Clifford+T, X-type gates, measurement).

Nothing is optimised, recognised or re-synthesised.  A gate is either a *leaf* whose exact matrix was compared with
the rule table, or it is expanded from its OWN ``definition``.  A rotation that is not a multiple of pi/4 is not
silently rounded: the translation stops with :class:`NeedsFractionalPhases` (such circuits need fractional phases,
see ``docs/QISKIT.md``).  Anything else that is not understood raises :class:`QiskitTranslationError`.

The result is a :class:`Translation` with the leaf list (qoptct operations), OpenQASM 3 text accepted by
``qoptct.verify`` and a map from every register of the circuit to its wires (LSB first).

Leaf rules (Qiskit numbering, first argument = least significant bit of the matrix):
  x cx ccx swap h s sdg t tdg z cz   exact gate classes only, matrix checked against the table
  mcx (C3X, C4X, MCXGate)            classical k-controlled X; matrix checked up to 9 wires, exact class and
                                     attributes above; ``mcx='expand'`` expands it from Qiskit's own definition
  p(theta) / u1(theta)               theta = k*pi/4 (|k - round(k)| < 1e-12)  ->  t^k
  measure, if_else (one bit == 1, body of x/z/cx/cz)
Global phase is allowed and recorded in ``global_phase``: every verdict is "up to a global phase", the same contract
as the qoptct provers (the phase must not depend on the input or on the measurement outcomes).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import ControlledGate
from qiskit.circuit.library import (C3XGate, C4XGate, CCXGate, CXGate, CZGate, HGate, MCXGate, PhaseGate, SGate,
                                    SdgGate, SwapGate, TGate, TdgGate, U1Gate, XGate, ZGate)


class QiskitTranslationError(ValueError):
    """The circuit is outside the strict supported language."""


class NeedsFractionalPhases(QiskitTranslationError):
    """A rotation that is not a multiple of pi/4 was met: fractional phases are needed.

    ``gates`` is a list of (path, gate name, angle) for the first offending rotations."""

    def __init__(self, message, gates=()):
        super().__init__(message)
        self.gates = list(gates)


_K_TOL = 1e-12
_MAT_QUBITS = 9

_NAMES = {XGate: ('x', 1), CXGate: ('cx', 2), CCXGate: ('ccx', 3), SwapGate: ('swap', 2), HGate: ('h', 1),
          SGate: ('s', 1), SdgGate: ('sdg', 1), TGate: ('t', 1), TdgGate: ('tdg', 1), ZGate: ('z', 1),
          CZGate: ('cz', 2)}
_LEAF = {}
for _cls, (_n, _w) in _NAMES.items():
    for _typ in (_cls, type(_cls())):
        _LEAF[_typ] = (_cls, _n, _w)
_MCX_TYPES = (MCXGate, C3XGate, C4XGate)
_PHASE_TYPES = (PhaseGate, U1Gate)
_ROTATIONS = frozenset({'u', 'u2', 'u3', 'r', 'rx', 'ry', 'rz', 'rv'})  # arbitrary angles: fractional phases (Q2)
_FORBIDDEN = frozenset({'reset', 'store', 'while_loop', 'for_loop', 'switch_case', 'break_loop', 'continue_loop',
                        'delay', 'global_phase'})
_INVERSE = {'x': 'x', 'cx': 'cx', 'ccx': 'ccx', 'mcx': 'mcx', 'swap': 'swap', 'h': 'h', 'z': 'z', 'cz': 'cz',
            's': 'sdg', 'sdg': 's', 't': 'tdg', 'tdg': 't'}


def _expected(name, width):
    """Matrix of a rule in Qiskit numbering (first argument = least significant bit)."""
    if name == 'h':
        r = 1 / math.sqrt(2)
        return np.array([[r, r], [r, -r]], dtype=complex)
    w8 = complex(math.cos(math.pi / 4), math.sin(math.pi / 4))
    m = np.zeros((1 << width, 1 << width), dtype=complex)
    for col in range(1 << width):
        row, val = col, 1
        if name == 'x':
            row ^= 1
        elif name == 'cx' and col & 1:
            row ^= 2
        elif name == 'ccx' and col & 3 == 3:
            row ^= 4
        elif name == 'swap' and (col & 1) != ((col >> 1) & 1):
            row ^= 3
        elif name == 'z' and col & 1:
            val = -1
        elif name == 's' and col & 1:
            val = 1j
        elif name == 'sdg' and col & 1:
            val = -1j
        elif name == 't' and col & 1:
            val = w8
        elif name == 'tdg' and col & 1:
            val = w8.conjugate()
        elif name == 'cz' and col & 3 == 3:
            val = -1
        m[row, col] = val
    return m


_MATRICES = {n: _expected(n, w) for n, w in set(_NAMES.values())}


@dataclass
class Translation:
    nq: int
    ops: list                     # ('x',q) ('cx',c,t) ('ccx',a,b,t) ('mcx',c1..ck,t) ('swap',a,b) ('h',q) ...
    registers: dict               # register name -> [wire, ...] (LSB first)
    nbits: int = 0                # ('measure',q,bit) ('if',bit,[ops])
    global_phase: float = 0.0
    primitives: dict = field(default_factory=dict)
    leaf_paths: list = field(default_factory=list)

    def counts(self):
        out = {}
        for op in self.ops:
            out[op[0]] = out.get(op[0], 0) + 1
        return out

    def is_classical(self):
        return all(op[0] in ('x', 'cx', 'ccx', 'mcx', 'swap') for op in self.ops)

    def qasm(self, nq=None):
        return emit_qasm(nq or self.nq, self.ops, self.nbits)

    def inverse_ops(self):
        if any(op[0] in ('measure', 'if') for op in self.ops):
            raise QiskitTranslationError('a circuit with measurements has no inverse')
        return [(_INVERSE[op[0]],) + tuple(op[1:]) for op in reversed(self.ops)]


def emit_qasm(nq, ops, nbits=0):
    lines = ['OPENQASM 3;', 'include "stdgates.inc";', 'qubit[%d] q;' % nq, 'bit[%d] c;' % max(1, nbits)]

    def g(name, wires):
        if len(set(wires)) != len(wires) or any(type(w) is not int or not 0 <= w < nq for w in wires):
            raise QiskitTranslationError('invalid wires %s %s' % (name, wires))
        return name + ' ' + ', '.join('q[%d]' % w for w in wires) + ';'
    for op in ops:
        if op[0] == 'measure':
            lines.append('c[%d] = measure q[%d];' % (op[2], op[1]))
        elif op[0] == 'if':
            lines.append('if (c[%d]) {%s}' % (op[1], ' '.join(g(o[0], o[1:]) for o in op[2])))
        else:
            lines.append(g(op[0], op[1:]))
    return '\n'.join(lines) + '\n'


def _angle(op):
    if len(op.params) != 1:
        raise QiskitTranslationError('phase gate with %d parameters' % len(op.params))
    try:
        theta = float(op.params[0])
    except (TypeError, ValueError) as exc:
        raise QiskitTranslationError('unbound or symbolic angle %r' % (op.params[0],)) from exc
    if not math.isfinite(theta):
        raise QiskitTranslationError('non-finite angle')
    return theta


def _phase_ops(k, q):
    return {0: [], 1: [('t', q)], 2: [('s', q)], 3: [('s', q), ('t', q)], 4: [('z', q)], 5: [('z', q), ('t', q)],
            6: [('sdg', q)], 7: [('tdg', q)]}[k % 8]


def translate(circuit, *, mcx='primitive', max_depth=64, max_operations=2_000_000):
    """QuantumCircuit -> :class:`Translation`.

    mcx: 'primitive' (default; C3X, C4X, MCX with all-ones control state are leaves) or 'expand'."""
    if not isinstance(circuit, QuantumCircuit):
        raise QiskitTranslationError('a QuantumCircuit is expected')
    if mcx not in ('primitive', 'expand'):
        raise QiskitTranslationError("mcx must be 'primitive' or 'expand'")
    if getattr(circuit, 'num_vars', 0) or getattr(circuit, 'num_stretches', 0):
        raise QiskitTranslationError('real-time variables are not supported')
    if circuit.parameters:
        raise QiskitTranslationError('unbound circuit parameters: %r' % (tuple(circuit.parameters),))
    ops, prim, leaf_paths, fractional = [], {}, [], []
    phase = [0.0]
    count = [0]

    def walk(block, qmap, cmap, path, depth, out):
        if depth > max_depth:
            raise QiskitTranslationError('%s: decomposition deeper than %d' % (path, max_depth))
        if block.parameters:
            raise QiskitTranslationError('%s: unbound parameters' % path)
        try:
            phase[0] += float(block.global_phase)
        except (TypeError, ValueError) as exc:
            raise QiskitTranslationError('%s: unbound global phase' % path) from exc
        if len(qmap) != block.num_qubits or len(set(qmap)) != len(qmap):
            raise QiskitTranslationError('%s: bad wire layout' % path)
        for idx, ins in enumerate(block.data):
            op = ins.operation
            name = getattr(op, 'name', type(op).__name__)
            cur = '%s/%d:%s' % (path, idx, name)
            count[0] += 1
            if count[0] > max_operations:
                raise QiskitTranslationError('%s: more than %d operations' % (cur, max_operations))
            if name in _FORBIDDEN:
                raise QiskitTranslationError('%s: %r is not supported' % (cur, name))
            try:
                qs = [qmap[block.find_bit(q).index] for q in ins.qubits]
                cs = [cmap[block.find_bit(c).index] for c in ins.clbits]
            except (KeyError, IndexError) as exc:
                raise QiskitTranslationError('%s: wire mapping failed: %s' % (cur, exc)) from exc
            if len(set(qs)) != len(qs):
                raise QiskitTranslationError('%s: repeated wire' % cur)
            if name == 'barrier':
                continue
            if name == 'measure':
                if len(qs) != 1 or len(cs) != 1:
                    raise QiskitTranslationError('%s: bad measure' % cur)
                out.append(('measure', qs[0], cs[0]))
                continue
            if name == 'if_else':
                cond = op.condition
                if not (isinstance(cond, tuple) and len(cond) == 2 and cond[1] in (1, True)) or len(cs) != 1:
                    raise QiskitTranslationError('%s: only a one-bit "== 1" condition is supported' % cur)
                if len(op.blocks) != 1 or len(qs) != op.blocks[0].num_qubits:
                    raise QiskitTranslationError('%s: only if without else is supported' % cur)
                sub = []
                walk(op.blocks[0], qs, [], cur, depth + 1, sub)
                if any(o[0] not in ('x', 'z', 'cx', 'cz') for o in sub):
                    raise QiskitTranslationError('%s: if body must consist of x/z/cx/cz' % cur)
                out.append(('if', cs[0], sub))
                continue
            if getattr(op, 'condition', None) is not None or cs:
                raise QiskitTranslationError('%s: classical control is not supported here' % cur)
            if type(op) in _PHASE_TYPES:
                theta = _angle(op)
                k = theta / (math.pi / 4)
                if abs(k - round(k)) > _K_TOL * max(1.0, abs(k)):
                    fractional.append((cur, name, theta))
                    continue
                out.extend(_phase_ops(int(round(k)), qs[0]))
                leaf_paths.append(cur)
                continue
            info = _LEAF.get(type(op))
            open_ctrl = isinstance(op, ControlledGate) and op.ctrl_state != (1 << op.num_ctrl_qubits) - 1
            if info is not None and not open_ctrl:
                _, lname, width = info
                if op.num_qubits != width or op.name != lname:
                    raise QiskitTranslationError('%s: wrong attributes of a known gate' % cur)
                try:
                    mat = np.asarray(op.to_matrix(), dtype=complex)
                except Exception as exc:
                    raise QiskitTranslationError('%s: matrix unavailable: %s' % (cur, exc)) from exc
                if not np.allclose(mat, _MATRICES[lname], atol=1e-13, rtol=0):
                    raise QiskitTranslationError('%s: matrix differs from the rule table' % cur)
                out.append((lname,) + tuple(qs))
                leaf_paths.append(cur)
                continue
            if type(op) in _MCX_TYPES and not open_ctrl and mcx == 'primitive':
                k = op.num_qubits
                if k < 4:
                    raise QiskitTranslationError('%s: mcx with %d wires' % (cur, k))
                if op.num_ctrl_qubits != k - 1 or op.ctrl_state != (1 << (k - 1)) - 1 or op.name not in (
                        'mcx', 'c3x', 'c4x'):
                    raise QiskitTranslationError('%s: bad mcx attributes' % cur)
                if type(op.base_gate) is not type(XGate()) and not isinstance(op.base_gate, XGate):
                    raise QiskitTranslationError('%s: mcx base gate is not X' % cur)
                if k <= _MAT_QUBITS and type(op) is not MCXGate:
                    want = np.eye(1 << k, dtype=complex)
                    lo, hi = (1 << (k - 1)) - 1, (1 << k) - 1
                    want[[lo, hi]] = want[[hi, lo]]
                    if not np.allclose(np.asarray(op.to_matrix(), dtype=complex), want, atol=1e-13, rtol=0):
                        raise QiskitTranslationError('%s: mcx matrix differs from the rule' % cur)
                out.append(('mcx',) + tuple(qs))
                prim['mcx%d' % (k - 1)] = prim.get('mcx%d' % (k - 1), 0) + 1
                leaf_paths.append(cur)
                continue
            if name in _ROTATIONS:
                par = [x for x in op.params]
                try:
                    ang = float(par[0])
                except (TypeError, ValueError, IndexError):
                    raise QiskitTranslationError('%s: unbound rotation parameter' % cur)
                fractional.append((cur, name, ang))
                continue
            try:
                definition = op.definition
            except Exception as exc:
                raise QiskitTranslationError('%s: no definition (%s)' % (cur, exc)) from exc
            if not isinstance(definition, QuantumCircuit) or not definition.data:
                raise QiskitTranslationError('%s: %r has no usable definition' % (cur, name))
            if definition.num_qubits != len(qs):
                raise QiskitTranslationError('%s: definition width mismatch' % cur)
            if definition.num_clbits:
                raise QiskitTranslationError('%s: definition with classical bits' % cur)
            walk(definition, qs, [], cur, depth + 1, out)

    walk(circuit, list(range(circuit.num_qubits)), list(range(circuit.num_clbits)), 'root', 0, ops)
    if fractional:
        shown = ', '.join('%s=%.6g' % (n, t) for _, n, t in fractional[:4])
        raise NeedsFractionalPhases(
            'fractional phases needed (Q2): %d rotation(s) that are not multiples of pi/4, e.g. %s at %s'
            % (len(fractional), shown, fractional[0][0]), fractional[:20])
    regs = {}
    for reg in circuit.qregs:
        if reg.name in regs:
            raise QiskitTranslationError('ambiguous register name %r' % reg.name)
        regs[reg.name] = [circuit.find_bit(q).index for q in reg]
    return Translation(nq=circuit.num_qubits, ops=ops, registers=regs, nbits=circuit.num_clbits,
                       global_phase=phase[0], primitives=prim, leaf_paths=leaf_paths)
