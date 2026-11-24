# SPDX-License-Identifier: Apache-2.0
"""check_qiskit(circuit, spec | reference) -> PROVED / REFUTED (with witness) / UNKNOWN.

spec is a :class:`Spec` (classical function + symbolic word-level form) or a plain classical function
``f(vals) -> dict`` (``vals``: register name -> int; the result lists only the registers that change; inputs that are
not listed stay as they are, clean registers must be 0 again).  A reference is another Qiskit circuit that is
compared through the miter ``circuit * reference^-1`` (all wires, clean wires are 0 at the start and at the end).

PROVED comes only from the exact prover (kprove, BDD sum-over-paths with symbolic outcomes) or from an exhaustive
enumeration of ALL inputs of a classical circuit.  REFUTED always carries a witness that was replayed by an
independent classical simulator.  Everything else (fractional phases, unsupported gate, budget, no symbolic form) is
UNKNOWN; an unproved sample never becomes PROVED.
"""
from __future__ import annotations

import os
import random
import tempfile
import time
from dataclasses import dataclass
from typing import Callable, Optional

from .translate import (NeedsFractionalPhases, QiskitTranslationError, Translation, emit_qasm, translate,
                        _INVERSE)

EXHAUSTIVE_BITS = 18


@dataclass
class Spec:
    """classical(vals)->{name: int}; symbolic(S, R)->{name: SymInt} (R: name -> SymInt of the input registers)."""
    name: str
    classical: Callable
    symbolic: Optional[Callable] = None
    description: str = ''
    msb_first: bool = False      # variable order hint for the BDD prover (comparators)


def _regs_info(tr, clean):
    clean = set(clean)
    unknown = clean - set(tr.registers)
    if unknown:
        raise QiskitTranslationError('unknown clean registers: %s' % sorted(unknown))
    inputs = [n for n in tr.registers if n not in clean]
    return inputs, sorted(clean)


def run_classical(ops, state):
    """Independent classical execution of X/CX/CCX/MCX/SWAP leaves on an integer state (wire 0 = LSB)."""
    for op in ops:
        g = op[0]
        if g == 'x':
            state ^= 1 << op[1]
        elif g in ('cx', 'ccx', 'mcx'):
            *cs, t = op[1:]
            if all((state >> c) & 1 for c in cs):
                state ^= 1 << t
        elif g == 'swap':
            a, b = op[1], op[2]
            if ((state >> a) ^ (state >> b)) & 1:
                state ^= (1 << a) | (1 << b)
        else:
            raise QiskitTranslationError('not a classical gate: %s' % g)
    return state


def _pack(tr, vals):
    s = 0
    for name, v in vals.items():
        for i, w in enumerate(tr.registers[name]):
            s |= ((v >> i) & 1) << w
    return s


def _unpack(tr, state):
    return {name: sum(((state >> w) & 1) << i for i, w in enumerate(ws)) for name, ws in tr.registers.items()}


def _expected(tr, spec, inputs, clean, vals):
    out = {n: (vals.get(n, 0) if n in inputs else 0) for n in tr.registers}
    vv = {n: 0 for n in tr.registers}
    vv.update(vals)
    res = spec.classical(vv) if isinstance(spec, Spec) else spec(vv)
    for n, v in res.items():
        if n not in tr.registers:
            raise QiskitTranslationError('spec names unknown register %r' % n)
        out[n] = v & ((1 << len(tr.registers[n])) - 1)
    return out


def _sample_inputs(tr, inputs, rng, count):
    edges = []
    for pat in (lambda w: 0, lambda w: (1 << w) - 1, lambda w: int('10' * w, 2) & ((1 << w) - 1),
                lambda w: int('01' * w, 2) & ((1 << w) - 1), lambda w: 1, lambda w: 1 << (w - 1)):
        edges.append({n: pat(len(tr.registers[n])) for n in inputs})
    for k in range(count):
        edges.append({n: rng.getrandbits(len(tr.registers[n])) if rng.random() < 0.7 else rng.choice(
            [0, (1 << len(tr.registers[n])) - 1, 1]) for n in inputs})
    return edges


def classical_validate(tr, spec, clean, samples=64, seed=1, exhaustive_bits=EXHAUSTIVE_BITS):
    """Independent replay of the translated leaves against the classical spec.
    -> (verdict, info): verdict 'ok' | 'mismatch' (witness in info); info['exhaustive'] tells if ALL inputs ran."""
    inputs, _ = _regs_info(tr, clean)
    width = sum(len(tr.registers[n]) for n in inputs)
    if width <= exhaustive_bits:
        cases = []
        for x in range(1 << width):
            vals, k = {}, 0
            for n in inputs:
                w = len(tr.registers[n])
                vals[n] = (x >> k) & ((1 << w) - 1)
                k += w
            cases.append(vals)
        exhaustive = True
    else:
        cases = _sample_inputs(tr, inputs, random.Random(seed), samples)
        exhaustive = False
    for vals in cases:
        got = _unpack(tr, run_classical(tr.ops, _pack(tr, vals)))
        want = _expected(tr, spec, inputs, clean, vals)
        if got != want:
            return 'mismatch', dict(witness=dict(vals), got=got, expected=want, exhaustive=exhaustive, cases=len(cases))
    return 'ok', dict(exhaustive=exhaustive, cases=len(cases))


def _order(tr, inputs, msb_first=False):
    """Interleave the input registers bit by bit (the order ripple carries need)."""
    words = [tr.registers[n][::-1] if msb_first else tr.registers[n] for n in inputs]
    m = max(map(len, words)) if words else 0
    return [w[i] for i in range(m) for w in words if i < len(w)]


def _prove_symbolic(tr, spec, clean, budget_sec, max_mb):
    from qoptct.verify.kprove import kprove
    from qoptct.verify.ksym import SymInt
    inputs, _ = _regs_info(tr, clean)
    clean_wires = sorted(w for n in clean for w in tr.registers[n])

    def full_spec(S):
        R = {n: S.reg(tr.registers[n]) for n in tr.registers}
        res = spec.symbolic(S, R)
        out = {}
        for n, v in res.items():
            if n not in tr.registers:
                raise ValueError('spec names unknown register %r' % n)
            if not isinstance(v, SymInt):
                raise ValueError('spec value of %r is not a SymInt' % n)
            out.update(S.out(tr.registers[n], v))
        return out
    with tempfile.TemporaryDirectory(prefix='qiskit_proof_') as d:
        path = os.path.join(d, 'c.qasm')
        with open(path, 'w') as f:
            f.write(tr.qasm())
        return kprove(path, spec=full_spec, n=tr.nq, clean_wires=clean_wires, order=_order(tr, inputs, spec.msb_first),
                      budget_sec=budget_sec, max_mb=max_mb, sample=False, confirm=True)


def _finish(out, t0):
    out['time'] = round(time.time() - t0, 2)
    return out


def _info(tr):
    return dict(nq=tr.nq, counts=tr.counts(), primitives=tr.primitives, global_phase=tr.global_phase,
                registers={k: len(v) for k, v in tr.registers.items()})


def check_qiskit(circuit, spec=None, reference=None, *, clean=(), ref_clean=(), budget_sec=300, max_mb=900,
                 mcx='primitive', samples=64, seed=1, exhaustive_bits=EXHAUSTIVE_BITS):
    """See the module docstring.  Never raises: bad input and internal errors give UNKNOWN."""
    t0 = time.time()
    out = dict(status='UNKNOWN', method='none', detail='')
    try:
        if (spec is None) == (reference is None):
            raise QiskitTranslationError('give exactly one of spec and reference')
        tr = translate(circuit, mcx=mcx)
        out['translation'] = _info(tr)
        if reference is not None:
            return _finish(_check_reference(out, tr, reference, clean, ref_clean, budget_sec, max_mb, mcx,
                                            exhaustive_bits), t0)
        if tr.nbits or not tr.is_classical():
            pass
        if tr.is_classical():
            verdict, info = classical_validate(tr, spec, clean, samples, seed, exhaustive_bits)
            out['independent_check'] = info
            if verdict == 'mismatch':
                out.update(status='REFUTED', method='replay', detail='the classical spec differs from the circuit',
                           witness=info['witness'])
                return _finish(out, t0)
            if info['exhaustive']:
                out.update(status='PROVED', method='exhaustive',
                           detail='all %d basis inputs replayed; classical circuit' % info['cases'])
                return _finish(out, t0)
        if not isinstance(spec, Spec) or spec.symbolic is None:
            out.update(method='classical-sample',
                       detail='no symbolic form of the spec: only a sampled check (not a proof)')
            return _finish(out, t0)
        r = _prove_symbolic(tr, spec, clean, budget_sec, max_mb)
        st = r.get('status')
        out.update(method='kprove-bdd', detail=str(r.get('detail', ''))[:400], prover=dict(
            time=r.get('time'), mem=r.get('mem'), counts=r.get('counts')))
        if st in ('PROVED', 'REFUTED'):
            out['status'] = st
            for k in ('witness', 'counterexample'):
                if r.get(k):
                    out[k] = r[k]
        else:
            out['detail'] = 'prover: %s %s' % (st, out['detail'])
    except NeedsFractionalPhases as e:
        out.update(method='translate', needs='fractional_phases', detail=str(e), rotations=e.gates)
    except QiskitTranslationError as e:
        out.update(method='translate', detail=str(e))
    except Exception as e:
        out.update(method='error', detail='%s: %s' % (type(e).__name__, e))
    return _finish(out, t0)


def _map_reference(tr, rtr, clean, ref_clean):
    """Wire map of the reference into the circuit: registers with the same name and size are shared, the other
    reference wires are fresh clean wires appended after the circuit."""
    wmap, nxt = {}, tr.nq
    for name, ws in rtr.registers.items():
        if name in tr.registers and len(tr.registers[name]) == len(ws) and name not in ref_clean:
            for a, b in zip(ws, tr.registers[name]):
                wmap[a] = b
    for w in range(rtr.nq):
        if w not in wmap:
            wmap[w] = nxt
            nxt += 1
    return wmap, nxt


def _check_reference(out, tr, reference, clean, ref_clean, budget_sec, max_mb, mcx, exhaustive_bits):
    from qoptct.verify.kprove import kprove
    rtr = translate(reference, mcx=mcx)
    out['reference'] = _info(rtr)
    if tr.nbits or rtr.nbits or any(o[0] in ('measure', 'if') for o in tr.ops + rtr.ops):
        raise QiskitTranslationError('reference mode needs unitary circuits')
    wmap, nq = _map_reference(tr, rtr, clean, ref_clean)
    inv = [(o[0],) + tuple(wmap[w] for w in o[1:]) for o in rtr.inverse_ops()]
    ops = list(tr.ops) + inv
    clean_wires = sorted({w for n in clean for w in tr.registers[n]} | set(range(tr.nq, nq)) | {
        wmap[w] for n in ref_clean for w in rtr.registers[n]})
    shared_inputs = [w for w in range(tr.nq) if w not in clean_wires]
    out['miter'] = dict(nq=nq, clean=len(clean_wires))
    if all(o[0] in ('x', 'cx', 'ccx', 'mcx', 'swap') for o in ops) and len(shared_inputs) <= exhaustive_bits:
        for x in range(1 << len(shared_inputs)):
            s = sum(((x >> i) & 1) << w for i, w in enumerate(shared_inputs))
            if run_classical(ops, s) != s:
                out.update(status='REFUTED', method='replay', witness=dict(input=s),
                           detail='circuit * reference^-1 is not the identity')
                return out
        out.update(status='PROVED', method='exhaustive', detail='miter replayed on all %d inputs' % (1 << len(shared_inputs)))
        return out
    with tempfile.TemporaryDirectory(prefix='qiskit_miter_') as d:
        path = os.path.join(d, 'm.qasm')
        with open(path, 'w') as f:
            f.write(emit_qasm(nq, ops))
        order = [w for w in range(nq) if w not in clean_wires]
        r = kprove(path, spec=lambda S: {}, n=nq, clean_wires=clean_wires, order=order, budget_sec=budget_sec,
                   max_mb=max_mb, sample=False, confirm=True)
    st = r.get('status')
    out.update(method='kprove-miter', detail=str(r.get('detail', ''))[:400])
    if st in ('PROVED', 'REFUTED'):
        out['status'] = st
        for k in ('witness', 'counterexample'):
            if r.get(k):
                out[k] = r[k]
    else:
        out['detail'] = 'prover: %s %s' % (st, out['detail'])
    return out
