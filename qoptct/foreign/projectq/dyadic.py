# SPDX-License-Identifier: Apache-2.0
"""ProjectQ programs with rotations pi*p/2^k: record them and hand them to the dyadic prover (``qoptct.dyadic``).

    from qoptct.foreign.projectq.dyadic import check_projectq_dyadic
    from qoptct.dyadic import QFTSpec

    def program(eng):
        q = eng.allocate_qureg(4)
        QFT | q
        return dict(registers=dict(q=q))

    r = check_projectq_dyadic(program, lambda regs: QFTSpec(regs['q'], swaps=False))
    r['status']   # 'PROVED' | 'REFUTED' | 'UNKNOWN'

The program is run on a ProjectQ engine whose own decomposition rules (``projectq.setups.decompositions``,
``projectq.libs.math``) bring it to X with up to two controls, H, Y, Z, S, T, swap, R (up to one control), Rz, Rx, Ry and
Ph.  No LocalOptimizer is used (it merges angles and rounds the sum again).  The target is a ``DCircuit`` reference, a
dyadic specification, or a function of the registers (name -> list of wire numbers) that returns one of them.

ANGLES.  ProjectQ stores every angle as a float rounded to 12 decimals (``ANGLE_PRECISION``) and drops what rounds to 0:
a rotation by pi/2^43 or less is a gate with angle 0.  An angle is therefore accepted only when it is EXACTLY what
ProjectQ stores for a rotation by p*pi/2^k, i.e. ``stored == round((p*pi/2^k) mod period, 12)`` (or the same value for the inverse gate, which
ProjectQ builds from the rounded angle of the original) for some p and some k <= ``MAX_K`` (38: then the multiples of pi/2^k are further apart than the rounding step, so the reading is unique).
An angle that is 0, within 1e-11 of 0 or of the period, or that matches no such value (for example the sum of two
rounded angles) makes the answer UNKNOWN: it could be several different rotations, so nothing is guessed.  The
accepted number of angles and the largest k are in ``result['angles']``.

Assumptions of the answer: those of ``qoptct.dyadic`` (for ``ArithSpec`` the fourier-shift lemma, in
``result['assumptions']``) and the exact angle reading above.  Allocated and deallocated wires of the recording are clean
ancillas (start in 0, must end in 0), as well as the wires of ``program``'s ``ancillas`` entry.  The verdict is up to a
global phase.
"""
from __future__ import annotations

import math
from fractions import Fraction

from ._record import ProjectQUnavailable, Unsupported, need_projectq

__all__ = ['record_dyadic', 'check_projectq_dyadic', 'exact_angle', 'MAX_K', 'ANGLE_DECIMALS']

ANGLE_DECIMALS = 12          # projectq.ops._basics.ANGLE_PRECISION
MAX_K = 38                   # pi / 2^38 = 1.1e-11 > 2 x 5e-13 rounding step: the reading of the angle is unique
NEAR_ZERO = 1e-11


def _stored_values(val, period):
    """What ProjectQ stores for a rotation by ``val``: the angle reduced and rounded once (a gate built from ``val``), or
    the inverse of the gate for ``-val`` (``get_inverse`` negates the already rounded angle and rounds again)."""
    direct = round(val % period, ANGLE_DECIMALS)
    s_neg = round((-val) % period, ANGLE_DECIMALS)
    return {direct, round((-s_neg) % period, ANGLE_DECIMALS)}


def exact_angle(stored, period):
    """The Fraction r with ``stored`` = what ProjectQ stores for a rotation by r*pi (see ``_stored_values``) and r = p/2^k, k <= MAX_K (smallest k), or
    raises Unsupported when the angle is 0, close to 0 or to the period, or is not such a value."""
    stored = float(stored)
    if not 0.0 <= stored < period:
        raise Unsupported('ProjectQ angle %r outside [0, %r)' % (stored, period))
    if stored < NEAR_ZERO or period - stored < NEAR_ZERO:
        raise Unsupported('ProjectQ angle %r is 0 or within 1e-11 of the period: it may be a rotation that ProjectQ '
                          'rounded away (ambiguous)' % stored)
    x = stored / math.pi
    for k in range(0, MAX_K + 1):
        p = round(x * (1 << k))
        if p == 0 or p % 2 == 0 and k > 0:
            continue                 # an even p belongs to a smaller k (already tried)
        for val in (p * math.pi / (1 << k), math.pi * p / (1 << k), math.ldexp(math.pi, -k) * p):
            if stored in _stored_values(val, period):
                return Fraction(p, 1 << k)
    raise Unsupported('ProjectQ angle %r is not exactly p*pi/2^k with k <= %d (its rounding to 12 decimals matches no '
                      'such value)' % (stored, MAX_K))


def _make_recorder():
    from projectq.cengines import BasicEngine
    from projectq.meta import get_control_count
    from projectq.ops import (AllocateQubitGate, DeallocateQubitGate, FlushGate, HGate, Ph, R, Rx, Ry, Rz, SwapGate,
                              XGate, YGate, ZGate)

    class Recorder(BasicEngine):
        def __init__(self):
            super().__init__()
            self.cmds = []

        def is_available(self, cmd):
            g = cmd.gate
            if isinstance(g, (AllocateQubitGate, DeallocateQubitGate, FlushGate)):
                return True
            k = get_control_count(cmd)
            if isinstance(g, XGate):
                return k <= 2
            if isinstance(g, (HGate, YGate, ZGate)) or str(g) in ('S', 'Sdag', 'T', 'Tdag'):
                return k == 0
            if isinstance(g, SwapGate):
                return k == 0
            if isinstance(g, R):
                return k <= 1
            if isinstance(g, (Rz, Rx, Ry, Ph)):
                return k == 0
            return False

        def receive(self, command_list):
            self.cmds.extend(command_list)
    return Recorder()


def _qureg_ids(q):
    return [x.id for x in q]


def record_dyadic(program):
    """Run ``program(eng)`` and return (DCircuit, info).  Raises Unsupported when something is not understood."""
    need_projectq()
    from projectq import MainEngine
    from projectq.cengines import AutoReplacer, DecompositionRuleSet, TagRemover
    from projectq.meta import get_control_count
    from projectq.ops import (AllocateQubitGate, DeallocateQubitGate, FlushGate, HGate, Ph, R, Rx, Ry, Rz, SwapGate,
                              XGate, YGate, ZGate)
    import projectq.libs.math as pqmath
    import projectq.setups.decompositions as pqdecomp
    from qoptct.dyadic import DCircuit, Pi

    rec = _make_recorder()
    rules = DecompositionRuleSet(modules=[pqmath, pqdecomp])
    eng = MainEngine(backend=rec, engine_list=[TagRemover(), AutoReplacer(rules), TagRemover()], verbose=False)
    meta = program(eng)
    eng.flush()
    if not isinstance(meta, dict):
        raise Unsupported('the program must return a dict (registers=..., ancillas=...)')
    allocated, freed, ops, nq = [], set(), [], 0
    ks, nang = [], 0
    for cmd in rec.cmds:
        g = cmd.gate
        if isinstance(g, FlushGate):
            continue
        if isinstance(g, AllocateQubitGate):
            for qr in cmd.qubits:
                allocated += [q.id for q in qr]
            continue
        if isinstance(g, DeallocateQubitGate):
            for qr in cmd.qubits:
                freed |= {q.id for q in qr}
            continue
        ctrl = [q.id for q in cmd.control_qubits]
        tgt = [q.id for qr in cmd.qubits for q in qr]
        k = len(ctrl)
        if isinstance(g, XGate):
            ops.append(('x' if k == 0 else 'cx' if k == 1 else 'ccx', ctrl + tgt, ()))
        elif isinstance(g, HGate):
            ops.append(('h', tgt, ()))
        elif isinstance(g, YGate):
            ops.append(('y', tgt, ()))
        elif isinstance(g, ZGate):
            ops.append(('z', tgt, ()))
        elif str(g) in ('S', 'Sdag', 'T', 'Tdag'):
            ops.append(({'S': 's', 'Sdag': 'sdg', 'T': 't', 'Tdag': 'tdg'}[str(g)], tgt, ()))
        elif isinstance(g, SwapGate):
            ops.append(('swap', tgt, ()))
        elif isinstance(g, (R, Rz, Rx, Ry, Ph)):
            period = 2 * math.pi if isinstance(g, (R, Ph)) else 4 * math.pi
            r = exact_angle(g.angle, period)
            nang += 1
            ks.append(r.denominator.bit_length() - 1)
            name = {R: 'p', Rz: 'rz', Rx: 'rx', Ry: 'ry'}.get(type(g)) if not isinstance(g, Ph) else 'gphase'
            if isinstance(g, Ph):
                ops.append(('gphase', [], (Pi(r.numerator, r.denominator),)))
            elif isinstance(g, R) and k == 1:
                ops.append(('cp', ctrl + tgt, (Pi(r.numerator, r.denominator),)))
            else:
                ops.append((name, tgt, (Pi(r.numerator, r.denominator),)))
        else:
            raise Unsupported('unexpected command %s' % cmd)
        nq = max([nq] + [q + 1 for q in ctrl + tgt])
    nq = max([nq, (max(allocated) + 1) if allocated else 0])
    c = DCircuit(nq, name='projectq')
    for name, qs, angles in ops:
        if name == 'gphase':
            c.add_global_phase(angles[0])
        else:
            c.add(name, qs, *angles)
    regs = {k: _qureg_ids(v) for k, v in (meta.get('registers') or {}).items()}
    anc = set(freed) | set(_qureg_ids(meta.get('ancillas') or []))
    info = dict(wires=nq, ops=len(ops), registers=regs, ancillas=sorted(anc),
                angles=dict(count=nang, max_k=max(ks) if ks else None), projectq=_pq_version())
    return c, info


def _pq_version():
    try:
        from importlib.metadata import version
        return version('projectq')
    except Exception:       # noqa: BLE001
        return None


def check_projectq_dyadic(program, target, *, phase='global', budget_sec=300.0, **kw):
    """Prove ``program`` (a function ``program(eng)`` that returns ``dict(registers={name: qureg}, ancillas=qureg)``)
    against ``target`` (see the module docstring).  Never raises, except ProjectQUnavailable.
    -> dict(status PROVED | REFUTED | UNKNOWN, method, detail, assumptions, angles, dyadic=<the answer of the prover>)."""
    res = dict(status='UNKNOWN', method='projectq-dyadic', detail='', assumptions=[])
    try:
        circ, info = record_dyadic(program)
    except ProjectQUnavailable:
        raise
    except Unsupported as e:
        return dict(res, detail='unsupported: %s' % e)
    except Exception as e:      # noqa: BLE001 - UNKNOWN is never a proof
        return dict(res, detail='%s: %s' % (type(e).__name__, e))
    res['circuit'] = info
    res['angles'] = info['angles']
    try:
        from qoptct.dyadic import DCircuit, prove_dyadic
        tgt = target(info['registers']) if callable(target) and not isinstance(target, DCircuit) else target
        anc = set(info['ancillas'])
        if isinstance(tgt, DCircuit):
            anc |= set(range(tgt.n, circ.n))
        r = prove_dyadic(circ, tgt, phase=phase, ancillas=sorted(anc), budget_sec=budget_sec, **kw)
    except Exception as e:      # noqa: BLE001
        return dict(res, detail='dyadic prover not run: %s: %s' % (type(e).__name__, e))
    v = r.get('verdict')
    res['dyadic'] = r
    res['assumptions'] = list(r.get('assumptions') or []) + [
        'ProjectQ angles are accepted only as exactly round(p*pi/2^k mod period, 12) with k <= %d (unique reading); '
        'a rotation that ProjectQ rounded to 0 cannot be told from no rotation, such a program is UNKNOWN' % MAX_K]
    if v in ('PROVED', 'REFUTED'):
        res['status'] = v
        res['detail'] = 'dyadic prover: %s (%s)' % (v, r.get('method'))
        if r.get('counterexample') is not None:
            res['witness'] = r['counterexample']
    else:
        res['detail'] = 'dyadic prover: %s %s' % (v, str(r.get('reason', ''))[:300])
        if v == 'INTERNAL':
            res['internal_error'] = True
    return res
