#!/usr/bin/env python
"""Dump ProjectQ circuits with rotations pi/2^k as JSON gate lists for scripts/dyadic_runs.py (needs ProjectQ).

    python scripts/dyadic_projectq_dump.py [--out results/dyadic/projectq] [--sizes 2,4,8,16,32,64]

Circuits: AddConstant(c) on an n-qubit quint (ProjectQ's Draper/Beauregard constant adder: QFT, R rotations,
inverse QFT) and QFT on n qubits, each compiled by ProjectQ's own AutoReplacer with the default decomposition rules
(projectq.setups.decompositions + projectq.libs.math), with or without LocalOptimizer (--optimizer on merges
consecutive rotations: the merged angles pi*p/2^k have large p and ProjectQ keeps only 12 decimals of every angle),
down to one of two backends:
  native: H, X, CNOT, Toffoli, Swap, single-qubit R / Rz / Rx / Ry / Ph, controlled R (one control);
  basic:  single-qubit gates and CNOT only (controlled R is decomposed by ProjectQ into CNOT, R/Rz and Ph).
Gate names in the dump are those of qoptct.dyadic (controls first); angles are ProjectQ's floats (radians).
"""
import argparse
import json
import os
import sys

from projectq import MainEngine
from projectq.cengines import AutoReplacer, BasicEngine, DecompositionRuleSet, LocalOptimizer, TagRemover
from projectq.meta import get_control_count
from projectq.ops import (AllocateQubitGate, DeallocateQubitGate, FlushGate, QFT, R, Rx, Ry, Rz, Ph, HGate, XGate,
                          YGate, ZGate, SwapGate)
import projectq.libs.math as pqmath
import projectq.setups.decompositions as pqdecomp
from projectq.libs.math import AddConstant


ROT = {R: 'p', Rz: 'rz', Rx: 'rx', Ry: 'ry'}


class Recorder(BasicEngine):
    """Backend that records the commands it accepts."""

    def __init__(self, mode):
        super().__init__()
        self.mode = mode
        self.cmds = []

    def is_available(self, cmd):
        g = cmd.gate
        if isinstance(g, (AllocateQubitGate, DeallocateQubitGate, FlushGate)):
            return True
        k = get_control_count(cmd)
        if isinstance(g, XGate):
            return k <= (2 if self.mode == 'native' else 1)
        if isinstance(g, (HGate, YGate, ZGate)) or str(g) in ('S', 'Sdag', 'T', 'Tdag'):
            return k == 0
        if isinstance(g, SwapGate):
            return k == 0 and self.mode == 'native'
        if isinstance(g, R):
            return k == 0 or (k == 1 and self.mode == 'native')
        if isinstance(g, (Rz, Rx, Ry, Ph)):
            return k == 0
        return False

    def receive(self, command_list):
        for cmd in command_list:
            if not isinstance(cmd.gate, (FlushGate, AllocateQubitGate, DeallocateQubitGate)):
                self.cmds.append(cmd)


def to_op(cmd):
    g = cmd.gate
    ctrl = [q.id for q in cmd.control_qubits]
    tgt = [q.id for qr in cmd.qubits for q in qr]
    k = len(ctrl)
    if isinstance(g, XGate):
        return ['x' if k == 0 else 'cx' if k == 1 else 'ccx' if k == 2 else 'mcx', ctrl + tgt, []]
    if isinstance(g, HGate):
        return ['h', tgt, []]
    if isinstance(g, YGate):
        return ['y', tgt, []]
    if isinstance(g, ZGate):
        return ['z', tgt, []]
    s = str(g)
    if s in ('S', 'Sdag', 'T', 'Tdag'):
        return [{'S': 's', 'Sdag': 'sdg', 'T': 't', 'Tdag': 'tdg'}[s], tgt, []]
    if isinstance(g, SwapGate):
        return ['swap', tgt, []]
    for cls, nm in ROT.items():
        if isinstance(g, cls):
            if k == 0:
                return [nm, tgt, [g.angle]]
            if cls is R and k == 1:
                return ['cp', ctrl + tgt, [g.angle]]
    if isinstance(g, Ph) and k == 0:
        return ['gphase', [], [g.angle]]
    raise ValueError('unexpected command %s' % cmd)


def compile_circuit(mode, build, optimizer=False):
    rec = Recorder(mode)
    rules = DecompositionRuleSet(modules=[pqmath, pqdecomp])
    engines = [TagRemover(), AutoReplacer(rules), TagRemover()]
    if optimizer:
        engines = [TagRemover(), LocalOptimizer(10), AutoReplacer(rules), TagRemover(), LocalOptimizer(10)]
    eng = MainEngine(backend=rec, engine_list=engines, verbose=False)
    meta = build(eng)
    eng.flush()
    ops, gph = [], 0.0
    nq = 0
    for cmd in rec.cmds:
        op = to_op(cmd)
        if op[0] == 'gphase':
            gph += op[2][0]
            continue
        ops.append(op)
        nq = max([nq] + [q + 1 for q in op[1]])
    meta.update(n=max(nq, meta.get('n', 0)), ops=ops, global_phase=gph, engine=mode, optimizer=optimizer)
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                                  'results', 'dyadic', 'projectq'))
    ap.add_argument('--sizes', default='2,4,8,16,32,64')
    ap.add_argument('--mode', default='native', choices=['native', 'basic'])
    ap.add_argument('--optimizer', default='off', choices=['on', 'off'])
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    try:
        from importlib.metadata import version
        pqv = version('projectq')
    except Exception:
        pqv = 'unknown'
    for n in [int(s) for s in a.sizes.split(',')]:
        c = (0xB5C3A1F7D94E62B8 * 0x9E3779B97F4A7C15 + n) % (1 << n) | 1    # mixed bits, odd

        def build_add(eng, n=n, c=c):
            q = eng.allocate_qureg(n)
            AddConstant(c) | q
            return dict(constant=c, b=[x.id for x in q], n=n)

        def build_qft(eng, n=n):
            q = eng.allocate_qureg(n)
            QFT | q
            return dict(wires=[x.id for x in q], swaps=False, n=n)
        for fam, build in (('pq_addconst', build_add), ('pq_qft', build_qft)):
            d = compile_circuit(a.mode, build, a.optimizer == 'on')
            d.update(generator='scripts/dyadic_projectq_dump.py', projectq_version=pqv, family=fam)
            name = (fam + ('_basic' if a.mode == 'basic' else '') + ('_opt' if a.optimizer == 'on' else '')
                    + '_%d.json' % n)
            with open(os.path.join(a.out, name), 'w') as f:
                json.dump(d, f)
            print(fam, n, a.mode, len(d['ops']), 'ops', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
