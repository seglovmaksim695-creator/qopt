# SPDX-License-Identifier: Apache-2.0
"""Random-circuit harness for qoptct.verify.dyadic (used by tests/test_dyadic.py and scripts/dyadic_fuzz.py).

Pairs (circuit, reference) are built from a random circuit by exact rewrites (equivalent pairs) and by mutations
(usually different).  The ground truth is the independent dense simulator qoptct.verify.dyadic.sim; a verdict is
consistent if PROVED => dense-equal with the reported global phase, REFUTED => dense-different (or, with
phase='exact', a global phase != 1).  Nothing here is used by the prover itself.
"""
import math
import random
from fractions import Fraction

import numpy as np

from .circuit import DCircuit, GATES, Pi
from . import sim

__all__ = ['rand_circuit', 'rewrite', 'mutate', 'dense_equal', 'consistent', 'run_fuzz', 'qft_circuit']


def qft_circuit(n, swaps=True, inverse=False, wires=None, nq=None):
    """Textbook QFT in Qiskit order (H on the top qubit first, then CP(pi/2^(j-k)) to the lower ones)."""
    wires = list(range(n)) if wires is None else list(wires)
    ops = []
    for j in reversed(range(n)):
        ops.append(('h', (wires[j],), ()))
        for k in reversed(range(j)):
            ops.append(('cp', (wires[j], wires[k]), (Fraction(1, 2 ** (j - k)),)))
    if swaps:
        for i in range(n // 2):
            ops.append(('swap', (wires[i], wires[n - 1 - i]), ()))
    if inverse:
        ops = [(o[0], o[1], tuple(-a for a in o[2])) for o in reversed(ops)]
    c = DCircuit(nq or max(wires) + 1)
    c.ops = ops
    return c


def rand_angle(rng, kmax=5):
    """pi * p / 2^k with |p / 2^k| < 2 (angles in (-2 pi, 2 pi): the 4 pi period of the half-angle gates matters)."""
    k = rng.randint(0, kmax)
    return Fraction(rng.randrange(-(2 << k), 2 << k) or 1, 1 << k)


def rand_circuit(rng, n, depth, kmax=5, names=None):
    names = names or [g for g in GATES if g != 'gphase']
    c = DCircuit(n)
    for _ in range(depth):
        while True:
            g = rng.choice(names)
            npar, ar = GATES[g]
            ar = ar if ar is not None else rng.randint(1, min(n, 4))
            if ar <= n:
                break
        qs = rng.sample(range(n), ar)
        c.add(g, qs, *[Pi(rand_angle(rng, kmax)) for _ in range(npar)])
    if rng.random() < 0.3:
        c.add_global_phase(Pi(rand_angle(rng, kmax)))
    return c


_INV = {'s': 'sdg', 'sdg': 's', 't': 'tdg', 'tdg': 't', 'sx': 'sxdg', 'sxdg': 'sx', 'cs': 'csdg', 'csdg': 'cs'}


def inverse_op(op):
    nm, qs, rs = op
    if nm in _INV:
        return (_INV[nm], qs, ())
    if nm in ('u', 'u3'):
        return (nm, qs, (-rs[0], -rs[2], -rs[1]))
    if nm == 'u2':
        return ('u', qs, (Fraction(-1, 2), -rs[1], -rs[0]))
    if nm == 'dcx':
        return ('dcx', (qs[1], qs[0]), ())
    if nm == 'csx':
        return ('h', (qs[1],), ()), ('csdg', qs, ()), ('h', (qs[1],), ())
    return (nm, qs, tuple(-r for r in rs))


_TOF = [('h', 2), ('cx', 1, 2), ('tdg', 2), ('cx', 0, 2), ('t', 2), ('cx', 1, 2), ('tdg', 2), ('cx', 0, 2),
        ('t', 1), ('t', 2), ('h', 2), ('cx', 0, 1), ('t', 0), ('tdg', 1), ('cx', 0, 1)]


def rewrite(rng, c, rounds=None):
    """An equivalent circuit, by exact identities (the dense simulator is the judge anyway)."""
    d = c.copy()
    ops = list(d.ops)
    for _ in range(rounds or rng.randint(1, 5)):
        kind = rng.randrange(11)
        i = rng.randrange(len(ops) + 1)
        js = {}
        for j, o in enumerate(ops):
            js.setdefault(o[0], []).append(j)
        if kind == 0 and ops:                                   # g g^dagger
            g = ops[rng.randrange(len(ops))]
            gi = inverse_op(g)
            ops[i:i] = [g] + (list(gi) if isinstance(gi[0], tuple) else [gi])
        elif kind == 1 and js.get('cp'):                        # CP(t) = P(t/2) P(t/2) CX P(-t/2) CX
            j = rng.choice(js['cp'])
            _, (a, b), (r,) = ops[j]
            ops[j:j + 1] = [('p', (a,), (r / 2,)), ('p', (b,), (r / 2,)), ('cx', (a, b), ()), ('p', (b,), (-r / 2,)),
                            ('cx', (a, b), ())]
        elif kind == 2 and js.get('rz'):                        # RZ(t) = exp(-i t/2) P(t)
            j = rng.choice(js['rz'])
            _, qs, (r,) = ops[j]
            ops[j] = ('p', qs, (r,))
            d.gphase = (d.gphase - r / 2) % 2
        elif kind == 3 and js.get('swap'):                      # SWAP = 3 CX
            j = rng.choice(js['swap'])
            _, (a, b), _ = ops[j]
            ops[j:j + 1] = [('cx', (a, b), ()), ('cx', (b, a), ()), ('cx', (a, b), ())]
        elif kind == 4 and len(ops) > 1:                        # commute gates on disjoint wires
            j = rng.randrange(len(ops) - 1)
            if not set(ops[j][1]) & set(ops[j + 1][1]):
                ops[j], ops[j + 1] = ops[j + 1], ops[j]
        elif kind == 5 and js.get('h'):                         # H = exp(-i pi/4) S SX S
            j = rng.choice(js['h'])
            qs = ops[j][1]
            ops[j:j + 1] = [('s', qs, ()), ('sx', qs, ()), ('s', qs, ())]
            d.gphase = (d.gphase - Fraction(1, 4)) % 2
        elif kind == 6 and (js.get('cx') or js.get('cz')):     # CX = H_t CZ H_t and back
            j = rng.choice(js.get('cx', []) + js.get('cz', []))
            nm, (a, b), _ = ops[j]
            ops[j:j + 1] = [('h', (b,), ()), ('cz' if nm == 'cx' else 'cx', (a, b), ()), ('h', (b,), ())]
        elif kind == 7:                                         # (S H)^3 = exp(i pi/4): needs the omega rule
            q = rng.randrange(c.n)
            ops[i:i] = [('s', (q,), ()), ('h', (q,), ())] * 3
            d.gphase = (d.gphase - Fraction(1, 4)) % 2
        elif kind == 8 and (js.get('ccx') or js.get('toffoli')):   # Toffoli = 7-T Clifford+T circuit
            j = rng.choice(js.get('ccx', []) + js.get('toffoli', []))
            qs = ops[j][1]
            ops[j:j + 1] = [(g[0], tuple(qs[k] for k in g[1:]), ()) for g in _TOF]
        elif kind == 9 and (js.get('mcp') or js.get('mcphase')):  # CCP(t) from CP(t/2) and CX
            j = rng.choice(js.get('mcp', []) + js.get('mcphase', []))
            nm, qs, (r,) = ops[j]
            if len(qs) == 3:
                c0, c1, t = qs
                ops[j:j + 1] = [('cp', (c1, t), (r / 2,)), ('cx', (c0, c1), ()), ('cp', (c1, t), (-r / 2,)),
                                ('cx', (c0, c1), ()), ('cp', (c0, t), (r / 2,))]
        elif kind == 10 and (js.get('rx') or js.get('ry')):     # RX = H RZ H, RY = S RX Sdg
            j = rng.choice(js.get('rx', []) + js.get('ry', []))
            nm, qs, (r,) = ops[j]
            mid = [('h', qs, ()), ('rz', qs, (r,)), ('h', qs, ())]
            ops[j:j + 1] = mid if nm == 'rx' else [('sdg', qs, ())] + mid + [('s', qs, ())]
    d.ops = ops
    return d


def mutate(rng, c):
    """A (usually) different circuit: shifted angle, lost gate, wrong qubit, global phase, extra gate."""
    d = c.copy()
    ops = list(d.ops)
    kind = rng.randrange(5)
    angled = [j for j, o in enumerate(ops) if o[2]]
    if kind == 0 and angled:
        j = rng.choice(angled)
        nm, qs, rs = ops[j]
        rs = list(rs)
        rs[rng.randrange(len(rs))] += Fraction(1, 2 ** rng.randint(0, 7))
        ops[j] = (nm, qs, tuple(rs))
    elif kind == 1 and ops:
        del ops[rng.randrange(len(ops))]
    elif kind == 2 and ops:
        j = rng.randrange(len(ops))
        nm, qs, rs = ops[j]
        free = [q for q in range(c.n) if q not in qs]
        if free:
            qs = list(qs)
            qs[rng.randrange(len(qs))] = rng.choice(free)
            ops[j] = (nm, tuple(qs), rs)
    elif kind == 3:
        d.gphase = (d.gphase + Fraction(1, 2 ** rng.randint(0, 7))) % 2
    else:
        e = rand_circuit(rng, c.n, 1)
        ops.insert(rng.randrange(len(ops) + 1), e.ops[0])
    d.ops = ops
    return d


def _widen(c, n):
    w = DCircuit(n)
    w.ops = list(c.ops)
    w.gphase = c.gphase
    return w


def dense_equal(c, ref=None, spec=None, anc=()):
    """(equal up to a global phase on the inputs with ancillas 0, lambda) by the dense simulator."""
    n = c.n
    U = sim.unitary(c)
    if ref is not None:
        V = sim.unitary(ref if ref.n == n else _widen(ref, n))
    else:
        V = sim.spec_unitary(spec, n)
    keep = [b for b in range(1 << n) if not any((b >> a) & 1 for a in anc)]
    Uk, Vk = U[:, keep], V[:, keep]
    lam = complex(np.vdot(Vk[:, 0], Uk[:, 0]))
    ok = abs(abs(lam) - 1) < 1e-7 and float(np.max(np.abs(Uk - lam * Vk))) < 1e-7
    return ok, lam


def consistent(r, ok, lam, exact=False):
    v = r['verdict']
    if v == 'PROVED':
        if not ok:
            return False
        g = r.get('global_phase_pi')
        if g is not None and abs(complex(math.cos(math.pi * float(g)), math.sin(math.pi * float(g))) - lam) > 1e-7:
            return False
        return not exact or abs(lam - 1) < 1e-7
    if v == 'REFUTED':
        return (not ok) or (exact and abs(lam - 1) > 1e-7)
    return v == 'UNKNOWN'


def run_fuzz(seed, cases, nmax=6, depth=(3, 14), selfcheck=False, budget_sec=20, kmax=5, log=None):
    """-> dict(counts, bad=[...]).  Half of the pairs are rewrites, half rewrite+mutation; 25% use phase='exact'."""
    from .api import prove_equal
    rng = random.Random(seed)
    st = dict(PROVED=0, REFUTED=0, UNKNOWN=0, equal=0, different=0, methods={}, bad=[], cases=cases)
    for t in range(cases):
        n = rng.randint(1, nmax)
        c = rand_circuit(rng, n, rng.randint(*depth), kmax=kmax)
        other = rewrite(rng, c) if rng.random() < 0.5 else mutate(rng, rewrite(rng, c))
        exact = rng.random() < 0.25
        r = prove_equal(other, c, phase='exact' if exact else 'global', budget_sec=budget_sec, selfcheck=selfcheck)
        ok, lam = dense_equal(other, ref=c)
        st['equal' if ok else 'different'] += 1
        v = r['verdict']
        st[v if v in ('PROVED', 'REFUTED', 'UNKNOWN') else 'UNKNOWN'] += 1
        key = '%s/%s' % (v, r.get('method'))
        st['methods'][key] = st['methods'].get(key, 0) + 1
        if v in ('ERROR', 'INTERNAL') or not consistent(r, ok, lam, exact):
            st['bad'].append(dict(case=t, verdict=v, dense_equal=ok, reason=str(r.get('reason'))[:200],
                                  circuit=other.to_qasm(), reference=c.to_qasm()))
        if log is not None:
            log(t, r, ok)
    return st


def _draper(n, kind='fixed'):
    """Draper adder written out as in Qiskit (QFT on the sum register, CP(pi/2^k) from a, inverse QFT)."""
    m = n if kind == 'fixed' else n + 1
    nq = 2 * n + (0 if kind == 'fixed' else 1)
    sumw = list(range(n, 2 * n)) + ([2 * n] if kind == 'half' else [])
    c = qft_circuit(m, wires=sumw, nq=nq)
    for j in range(n):
        for k in range(m - j):
            c.ops.append(('cp', (j, sumw[m - 1 - (j + k)]), (Fraction(1, 2 ** k),)))
    c.ops += qft_circuit(m, wires=sumw, nq=nq, inverse=True).ops
    return c, list(range(n)), list(range(n, 2 * n)), (2 * n if kind == 'half' else None)


def _const_adder(m, const, ctrl):
    """ProjectQ-style constant adder: QFT (no swaps), P(pi/2^(i-j)) per set bit, inverse QFT."""
    n = m + (1 if ctrl else 0)
    c = qft_circuit(m, swaps=False, nq=n)
    for i in range(m):
        for j in range(i, -1, -1):
            if (const >> j) & 1:
                q = m - 1 - i
                if ctrl:
                    c.ops.append(('cp', (m, q), (Fraction(1, 1 << (i - j)),)))
                else:
                    c.ops.append(('p', (q,), (Fraction(1, 1 << (i - j)),)))
    c.ops += qft_circuit(m, swaps=False, nq=n, inverse=True).ops
    return c


def _rgqft(n):
    """RGQFT multiplier written out as in Qiskit (QFT on out, CCP(2 pi 2^(2n-i-j-k)), inverse QFT)."""
    m = 2 * n
    out = list(range(2 * n, 4 * n))
    c = qft_circuit(m, wires=out, nq=4 * n)
    for j in range(1, n + 1):
        for i in range(1, n + 1):
            for k in range(1, m + 1):
                ang = Fraction(2) * Fraction(2) ** (2 * n - i - j - k)
                c.ops.append(('mcp', (n - j, 2 * n - i, out[m - k]), (ang,)))
    c.ops += qft_circuit(m, wires=out, nq=4 * n, inverse=True).ops
    return c


def run_spec_fuzz(seed, cases, budget_sec=20, log=None):
    """Specifications against the dense simulator: QFT variants, Draper adders, constant adders (controlled or not),
    RGQFT multipliers, each correct or mutated, up to 10 qubits."""
    from .api import prove_spec
    from .spec import QFTSpec, adder_spec, const_adder_spec, multiplier_spec
    rng = random.Random(seed)
    st = dict(PROVED=0, REFUTED=0, UNKNOWN=0, equal=0, different=0, methods={}, bad=[], cases=cases, kinds={})
    for t in range(cases):
        kind = rng.choice(['qft', 'qft', 'draper', 'const', 'rgqft'])
        if kind == 'qft':
            n = rng.randint(1, 10)
            m = rng.randint(1, n)
            wires = rng.sample(range(n), m)
            sw, inv = rng.random() < 0.5, rng.random() < 0.5
            c = qft_circuit(m, swaps=sw, inverse=inv, wires=wires, nq=n)
            spec = QFTSpec(wires, swaps=sw, inverse=inv)
        elif kind == 'draper':
            n = rng.randint(1, 4)
            c, a, b, cout = _draper(n, rng.choice(['fixed', 'half']))
            spec = adder_spec(a, b, cout)
        elif kind == 'const':
            m = rng.randint(1, 8)
            const = rng.randrange(1 << m)
            ctrl = rng.random() < 0.5
            c = _const_adder(m, const, ctrl)
            spec = const_adder_spec(const, list(range(m)), controls=(m,) if ctrl else ())
        else:
            n = rng.randint(1, 2)
            c = _rgqft(n)
            spec = multiplier_spec(range(n), range(n, 2 * n), range(2 * n, 4 * n))
        if rng.random() < 0.5:
            c = mutate(rng, c)
        elif rng.random() < 0.5:
            c = rewrite(rng, c)
        r = prove_spec(c, spec, budget_sec=budget_sec, selfcheck=False)
        ok, lam = dense_equal(c, spec=spec)
        st['equal' if ok else 'different'] += 1
        v = r['verdict']
        st[v if v in ('PROVED', 'REFUTED', 'UNKNOWN') else 'UNKNOWN'] += 1
        st['kinds'][kind] = st['kinds'].get(kind, 0) + 1
        key = '%s/%s' % (v, r.get('method'))
        st['methods'][key] = st['methods'].get(key, 0) + 1
        if v in ('ERROR', 'INTERNAL') or not consistent(r, ok, lam):
            st['bad'].append(dict(case=t, kind=kind, verdict=v, dense_equal=ok, reason=str(r.get('reason'))[:200],
                                  circuit=c.to_qasm()))
        if log is not None:
            log(t, r, ok)
    return st
