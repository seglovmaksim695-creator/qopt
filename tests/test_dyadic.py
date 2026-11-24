"""Tests of qoptct.verify.dyadic (exact path sums with phases in Z_{2^K}).

Correctness is checked against things that share no code with the path-sum engine:
  * the lift and every rewrite rule against an exact enumeration of the operator in Z[zeta] (random path sums);
  * random circuits with rotations pi*p/2^k (equivalent pairs and mutants) against the independent dense simulator
    (qoptct.verify.dyadic.sim): a PROVED must be dense-equal with the same global phase, a REFUTED dense-different;
  * canaries on QFT and Draper adders (shifted angle, lost rotation, wrong qubit, global phase): never PROVED;
  * the cases of tests/data/kprove and the audit circuits: a wrong circuit is never PROVED.
Fast by default; the large runs (n = 32 / 64) are marked slow.
"""
import json
import math
import os
import random
from fractions import Fraction

import pytest

from qoptct.verify.dyadic import (DCircuit, Pi, prove_equal, prove_spec, prove_qasm, parse_qasm, QFTSpec, ArithSpec,
                                  PermutationSpec, adder_spec, const_adder_spec, multiplier_spec, NotDyadic)
from qoptct.verify.dyadic import sim
from qoptct.verify.dyadic.core import DPathSum, ZERO, ONE, peval
from qoptct.verify.dyadic.api import _canon, _lam_equal
from qoptct.verify.dyadic.circuit import recognise_float, to_angle, GATES

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')


# ----------------------------------------------------------------------------- helpers
from qoptct.verify.dyadic.fuzz import (qft_circuit, dense_equal, consistent, rand_circuit, rewrite, mutate,
                                       run_fuzz)


# ----------------------------------------------------------------------------- angles and parsing
def test_angle_recognition():
    assert recognise_float(math.pi / 4) == Fraction(1, 4)
    assert recognise_float(-math.pi / 2 ** 40) == 4 - Fraction(1, 2 ** 40)
    assert recognise_float(math.pi / 2 ** 63) == Fraction(1, 2 ** 63)
    assert recognise_float(2 * math.pi * 2 ** 40) == 0
    assert recognise_float(2 * math.pi - math.pi / 2 ** 45) == 2 - Fraction(1, 2 ** 45)
    assert recognise_float(4 * math.pi - math.pi / 2 ** 45) == 4 - Fraction(1, 2 ** 45)
    assert recognise_float(round(math.pi / 2 ** 30, 12), abs_tol=2e-11) == Fraction(1, 2 ** 30)   # ProjectQ rounding
    with pytest.raises(NotDyadic):
        recognise_float(round(math.pi / 2 ** 30, 12))
    for tiny in (round(math.pi / 2 ** 36, 12), 1e-12, 4 * math.pi - 1e-11):
        with pytest.raises(NotDyadic):                 # below 64 x the tolerance: ambiguous
            recognise_float(tiny, abs_tol=2e-11)
    assert recognise_float(0.0, abs_tol=2e-11) == 0
    assert recognise_float(3 * math.pi / 8) == Fraction(3, 8)
    for bad in (0.1, 1.0, math.e, math.pi / 3, math.pi * (2 ** 21 + 1) / 2 ** 22):
        with pytest.raises(NotDyadic):
            recognise_float(bad)
    assert to_angle('pi/4') == Fraction(1, 4)
    assert to_angle('-3*pi/8') == Fraction(29, 8)
    assert to_angle('pi*2**-5') == Fraction(1, 32)
    assert to_angle('2*pi/2^5') == Fraction(1, 16)
    assert to_angle(Pi(-1, 2)) == Fraction(7, 2)
    with pytest.raises(NotDyadic):
        to_angle('pi/3')
    with pytest.raises(TypeError):
        to_angle(Fraction(1, 4))           # ambiguous: refused
    with pytest.raises(NotDyadic):
        to_angle('__import__("os")')


def test_half_angle_gates_have_period_4pi():
    """RZ(t + 2 pi) = -RZ(t); CRZ(t + 2 pi) has a relative phase -1 on the control: angles are kept modulo 4 pi."""
    for g, qs in (('rz', (0,)), ('rx', (0,)), ('ry', (0,)), ('u', (0,)), ('crz', (0, 1)), ('rzz', (0, 1))):
        args = (Pi(1, 4), Pi(1, 2), Pi(1, 8)) if g == 'u' else (Pi(1, 4),)
        args2 = (Pi(9, 4),) + args[1:]
        c1 = DCircuit(2).add(g, qs, *args)
        c2 = DCircuit(2).add(g, qs, *args2)
        r_exact = prove_equal(c1, c2, phase='exact')
        r_glob = prove_equal(c1, c2)
        assert r_exact['verdict'] == 'REFUTED', (g, r_exact)
        if g == 'crz':
            assert r_glob['verdict'] == 'REFUTED', r_glob
        else:
            assert r_glob['verdict'] == 'PROVED' and r_glob['global_phase'] == 'pi*1/1', (g, r_glob)
        c3 = DCircuit(2).add(g, qs, *((Pi(17, 4),) + args[1:]))      # t + 4 pi: the same gate
        assert prove_equal(c1, c3, phase='exact')['verdict'] == 'PROVED'


def test_qasm_parser_and_unsupported():
    txt = '''OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
gate mygate(t) x, y { cp(t/2) x, y; h y; }
mygate(pi/2) a[0], b[0];
cu1(pi/8) a[1], b[0];
u3(pi/2, 0, pi) a[0];
barrier a[0], b[0];
'''
    c = parse_qasm(txt)
    assert c.n == 3 and [o[0] for o in c.ops] == ['cp', 'h', 'cu1', 'u3']
    assert c.ops[0][2] == (Fraction(1, 4),)
    r = prove_qasm(txt, txt)
    assert r['verdict'] == 'PROVED'
    for bad in ('OPENQASM 3.0;\nqubit[1] q;\nbit[1] c;\nc[0] = measure q[0];\n',
                'OPENQASM 3.0;\nqubit[1] q;\nreset q[0];\n',
                'OPENQASM 3.0;\nqubit[1] q;\nrz(0.1) q[0];\n',
                'OPENQASM 3.0;\nqubit[2] q;\nch q[0], q[1];\n'):
        r = prove_qasm(bad, 'OPENQASM 3.0;\nqubit[2] q;\n')
        assert r['verdict'] == 'UNKNOWN', (bad, r)
    r = prove_qasm('OPENQASM 3.0;\nqubit[1] q;\nh q[3];\n', 'OPENQASM 3.0;\nqubit[1] q;\n')
    assert r['verdict'] == 'ERROR'


# ----------------------------------------------------------------------------- lift and rules vs exact enumeration
def _rand_poly(rng, vars_, maxterms):
    P = set()
    for _ in range(rng.randint(0, maxterms)):
        mo = 0
        for v in vars_:
            if rng.random() < 0.4:
                mo |= 1 << v
        P ^= {mo}
    return frozenset(P)


def test_lift_exhaustive():
    rng = random.Random(3)
    for _ in range(300):
        K = rng.randint(3, 9)
        ps = DPathSum(K)
        vs = [ps.var('x') for _ in range(rng.randint(1, 5))]
        P = _rand_poly(rng, vs, 6)
        c = rng.randrange(1, 1 << K)
        ps.add_lift(P, c)
        for val in range(1 << len(vs)):
            f = sum(cc for mo, cc in ps.ph.items() if mo & val == mo) % (1 << K)
            assert f == (c * peval(P, val)) % (1 << K)


def _operator(ps, xs):
    """Exact operator: {(x assignment, output): (canonical element, e2)}."""
    ys = sorted(ps.ys)
    out = {}
    for xa in range(1 << len(xs)):
        xval = sum(1 << v for i, v in enumerate(xs) if (xa >> i) & 1)
        acc = {}
        for ya in range(1 << len(ys)):
            val = xval | sum(1 << v for i, v in enumerate(ys) if (ya >> i) & 1)
            o, f = ps.evaluate(val)
            d = acc.setdefault(o, {})
            d[f] = d.get(f, 0) + 1
        for o, d in acc.items():
            d = _canon(d, ps.N)
            if d:
                out[(xa, o)] = d
    return out, ps.e2


def test_rules_sound_on_random_path_sums():
    """Random path sums built to trigger Elim / HH / omega: the exact operator before and after the reduction."""
    rng = random.Random(11)
    fired = dict(elim=0, hh=0, omega=0)
    for t in range(600):
        K = rng.randint(3, 7)
        N = 1 << K
        ps = DPathSum(K)
        xs = [ps.var('x') for _ in range(rng.randint(0, 3))]
        ys = [ps.var('y') for _ in range(rng.randint(1, 5))]
        allv = xs + ys
        for q in range(rng.randint(1, 3)):
            ps.add_wire(_rand_poly(rng, [v for v in allv if rng.random() < 0.5], 2))
        for _ in range(rng.randint(1, 8)):
            mo = 0
            for v in allv:
                if rng.random() < 0.35:
                    mo |= 1 << v
            c = rng.choice([N >> 1, N >> 1, N >> 1, N >> 2, 3 * (N >> 2), rng.randrange(1, N)])
            ps.addc(mo, c)
        before, e_b = _operator(ps, xs)
        ps.dirty |= set(ps.ys)
        ps.reduce()
        for k in fired:
            fired[k] += ps.stats[k]
        if ps.zero:
            assert not before, 'zero flag on a nonzero operator'
            continue
        after, e_a = _operator(ps, xs)
        assert set(before) == set(after), t
        for key in before:
            assert _lam_equal((before[key], e_b), (after[key], e_a), N), (t, key)
    assert all(v > 20 for v in fired.values()), fired


# ----------------------------------------------------------------------------- random circuits vs dense
def test_fuzz_against_dense_fast():
    st = run_fuzz(7, 200)
    assert not st['bad'], st['bad'][:3]
    assert st['PROVED'] >= 60 and st['REFUTED'] >= 40, st


@pytest.mark.slow
def test_fuzz_against_dense_large():
    st = run_fuzz(2026, 1500, nmax=10, depth=(3, 30))
    assert not st['bad'], st['bad'][:3]


def test_specs_against_dense():
    rng = random.Random(5)
    for t in range(60):
        n = rng.randint(1, 6)
        m = rng.randint(1, n)
        wires = rng.sample(range(n), m)
        swaps, inv = rng.random() < 0.5, rng.random() < 0.5
        c = qft_circuit(m, swaps=swaps, inverse=inv, wires=wires, nq=n)
        if rng.random() < 0.5:
            c = mutate(rng, c)
        spec = QFTSpec(wires, swaps=swaps, inverse=inv)
        r = prove_spec(c, spec, selfcheck=False, budget_sec=20)
        ok, lam = dense_equal(c, spec=spec)
        assert consistent(r, ok, lam), (t, r)
    # random small arithmetic: Draper-style constant adders, controlled or not
    for t in range(40):
        m = rng.randint(1, 5)
        const = rng.randrange(1 << m)
        ctrl = rng.random() < 0.5
        n = m + (1 if ctrl else 0)
        b = list(range(m))
        c = qft_circuit(m, swaps=False, nq=n)
        for i in range(m):
            for j in range(i, -1, -1):
                if (const >> j) & 1:
                    ang = Pi(1, 1 << (i - j))
                    if ctrl:
                        c.add('cp', (m, b[m - 1 - i]), ang)
                    else:
                        c.add('p', (b[m - 1 - i],), ang)
        inv = qft_circuit(m, swaps=False, inverse=True, nq=n)
        for o in inv.ops:
            c.add(o[0], o[1], *[Pi(x.numerator, x.denominator) for x in o[2]])
        if rng.random() < 0.5:
            c = mutate(rng, c)
        spec = const_adder_spec(const, b, controls=(m,) if ctrl else ())
        r = prove_spec(c, spec, selfcheck=False, budget_sec=20)
        ok, lam = dense_equal(c, spec=spec)
        assert consistent(r, ok, lam), (t, r, ok)


def test_permutation_spec():
    rng = random.Random(9)
    for t in range(60):
        n = rng.randint(1, 5)
        c = rand_circuit(rng, n, rng.randint(1, 10), names=['x', 'cx', 'ccx', 'swap', 'cswap', 'mcx'])
        if rng.random() < 0.5:
            c.add('h', (0,))
            c.add('z', (0,))
            c.add('h', (0,))                       # = X: still a permutation
        func = _perm_of(c)
        spec = PermutationSpec(n, func if rng.random() < 0.7 else (lambda b, f=func: f(b) ^ 1))
        r = prove_spec(c, spec, selfcheck=False)
        ok, lam = dense_equal(c, spec=spec)
        assert consistent(r, ok, lam), (t, r)


def _perm_of(c):
    U = sim.unitary(c)
    import numpy as np
    tab = [int(np.argmax(np.abs(U[:, b]))) for b in range(1 << c.n)]
    return lambda b: tab[b]


# ----------------------------------------------------------------------------- QFT, Draper, canaries
@pytest.mark.parametrize('n', [1, 2, 3, 5, 8, 16])
def test_qft_proved(n):
    assert prove_spec(qft_circuit(n), QFTSpec(range(n)))['verdict'] == 'PROVED'
    assert prove_spec(qft_circuit(n, swaps=False), QFTSpec(range(n), swaps=False))['verdict'] == 'PROVED'
    assert prove_spec(qft_circuit(n, inverse=True), QFTSpec(range(n), inverse=True))['verdict'] == 'PROVED'
    c = qft_circuit(n)
    for o in qft_circuit(n, inverse=True).ops:
        c.add(o[0], o[1], *[Pi(x.numerator, x.denominator) for x in o[2]])
    assert prove_equal(c, DCircuit(n))['verdict'] == 'PROVED'           # QFT QFT^dagger = I


def canaries_qft(n, rng):
    """Targeted mutants of the n-qubit QFT: every one is wrong."""
    base = qft_circuit(n)
    cps = [i for i, o in enumerate(base.ops) if o[0] == 'cp']
    small = min(cps, key=lambda i: base.ops[i][2][0])           # the smallest rotation
    out = []
    for i in [small] + rng.sample(cps, min(5, len(cps))):
        for k in sorted({1, n - 1, n, n + 1}):
            c = base.copy()
            nm, qs, rs = c.ops[i]
            c.ops[i] = (nm, qs, (rs[0] + Fraction(1, 2 ** k),))
            out.append(('shift_pi/2^%d' % k, c))
    for i in [small] + rng.sample(cps, min(4, len(cps))):
        c = base.copy()
        del c.ops[i]
        out.append(('lost_rotation', c))
    for i in rng.sample(cps, min(5, len(cps))):
        c = base.copy()
        nm, (a, b), rs = c.ops[i]
        free = [q for q in range(n) if q not in (a, b)]
        if free:
            c.ops[i] = (nm, (a, rng.choice(free)), rs)
            out.append(('wrong_qubit', c))
    if n >= 2:
        c = base.copy()
        sw = [i for i, o in enumerate(c.ops) if o[0] == 'swap']
        if sw:
            del c.ops[sw[0]]
            out.append(('lost_swap', c))
        c = base.copy()
        hs = [i for i, o in enumerate(c.ops) if o[0] == 'h']
        q = c.ops[hs[1]][1][0]
        j = min(i for i, o in enumerate(c.ops) if o[0] == 'cp' and q in o[1])
        h = c.ops.pop(hs[1])
        c.ops.insert(j, h)                          # H moved before a CP on its own qubit (does not commute)
        out.append(('order', c))
    c = base.copy()
    c.gphase = Fraction(1, 2 ** n)
    out.append(('global_phase', c))
    return out


@pytest.mark.parametrize('n', [8, 16])
def test_qft_canaries(n):
    rng = random.Random(n)
    cans = canaries_qft(n, rng)
    assert len(cans) >= 30
    for name, c in cans:
        r = prove_spec(c, QFTSpec(range(n)), phase='exact' if name == 'global_phase' else 'global', budget_sec=60)
        assert r['verdict'] == 'REFUTED', (name, r)


def draper(n, kind='fixed'):
    """DraperQFTAdder written out as in Qiskit (QFT on the sum register, CP(pi/2^k) from a, inverse QFT)."""
    m = n if kind == 'fixed' else n + 1
    nq = 2 * n + (0 if kind == 'fixed' else 1)
    sumw = list(range(n, 2 * n)) + ([2 * n] if kind == 'half' else [])
    c = qft_circuit(m, wires=sumw, nq=nq)
    for j in range(n):
        for k in range(m - j):
            c.add('cp', (j, sumw[m - 1 - (j + k)]), Pi(1, 2 ** k))
    for o in qft_circuit(m, wires=sumw, nq=nq, inverse=True).ops:
        c.add(o[0], o[1], *[Pi(x.numerator, x.denominator) for x in o[2]])
    return c, list(range(n)), list(range(n, 2 * n)), (2 * n if kind == 'half' else None)


@pytest.mark.parametrize('n', [1, 2, 3, 4, 8])
@pytest.mark.parametrize('kind', ['fixed', 'half'])
def test_draper_proved_and_dense(n, kind):
    c, a, b, cout = draper(n, kind)
    spec = adder_spec(a, b, cout)
    r = prove_spec(c, spec)
    assert r['verdict'] == 'PROVED' and r['assumptions'], r
    if c.n <= 9:
        ok, lam = dense_equal(c, spec=spec)
        assert ok and consistent(r, ok, lam)


def test_draper_small_without_lemma():
    """The same adders against the classical permutation: no Fourier lemma involved (direct normal form)."""
    for n in (1, 2, 3, 4):
        c, a, b, _ = draper(n)
        spec = PermutationSpec(2 * n, lambda v, n=n: (v & ((1 << n) - 1)) | (((((v >> n) + (v & ((1 << n) - 1)))
                                                                              % (1 << n))) << n))
        r = prove_spec(c, spec)
        assert r['verdict'] == 'PROVED' and r['method'] == 'normal-form' and 'assumptions' not in r, r


def test_draper_canaries():
    rng = random.Random(1)
    n = 6
    c0, a, b, _ = draper(n)
    spec = adder_spec(a, b)
    cps = [i for i, o in enumerate(c0.ops) if o[0] == 'cp']
    cnt = 0
    for i in rng.sample(cps, 12):
        for k in (0, 3, n):
            c = c0.copy()
            nm, qs, rs = c.ops[i]
            c.ops[i] = (nm, qs, (rs[0] + Fraction(1, 2 ** k),))
            assert prove_spec(c, spec, budget_sec=30)['verdict'] == 'REFUTED'
            cnt += 1
    for i in rng.sample(cps, 6):
        c = c0.copy()
        del c.ops[i]
        assert prove_spec(c, spec, budget_sec=30)['verdict'] == 'REFUTED'
        cnt += 1
    assert prove_spec(c0, adder_spec(a[::-1], b))['verdict'] == 'REFUTED'
    assert prove_spec(c0, adder_spec(a, b[::-1]))['verdict'] == 'REFUTED'
    c = c0.copy()
    c.gphase = Fraction(1, 64)
    assert prove_spec(c, spec, phase='exact')['verdict'] == 'REFUTED'
    assert prove_spec(c, spec)['verdict'] == 'PROVED'
    assert cnt >= 30


def test_spec_domain_checked_first():
    c, a, b, _ = draper(3)
    for bad in (ArithSpec([3, 3, 4], [(1, (0,))]), ArithSpec([3, 4, 5], [(1, (3,))]), QFTSpec([0, 9]),
                QFTSpec([0, 0]), ArithSpec([3, 4, 9], [(1, (0,))])):
        assert prove_spec(c, bad)['verdict'] == 'ERROR'
    with pytest.raises(ValueError):
        ArithSpec([3, 4, 5], [(Fraction(1, 2), (0,))])


def test_never_proved_on_errors():
    assert prove_equal('not a circuit', DCircuit(1))['verdict'] == 'ERROR'
    c = DCircuit(2)
    c.add('h', (0,))
    r = prove_equal(c, DCircuit(3))
    assert r['verdict'] == 'ERROR'


# ----------------------------------------------------------------------------- old cases: wrong circuits never PROVED
def test_kprove_cases_never_falsely_proved():
    cases = os.path.join(DATA, 'kprove', 'cases.json')
    if not os.path.exists(cases):
        pytest.skip('kprove cases not present')
    seen = 0
    for case in json.load(open(cases)):
        q = os.path.join(DATA, 'kprove', case['qasm'])
        qc = os.path.join(DATA, 'kprove', case['qc'])
        if not (os.path.exists(q) and os.path.exists(qc)):
            continue
        r = prove_qasm(q, qc, clean=case.get('clean', False), budget_sec=20)
        seen += 1
        if case.get('truth') != 'ok':
            assert r['verdict'] != 'PROVED', case
        assert r['verdict'] != 'INTERNAL', (case, r)
    assert seen > 0


def test_audit_circuits_not_proved():
    audit = os.path.join(DATA, 'audit')
    if not os.path.isdir(audit):
        pytest.skip('audit data not present')
    for name in ('hole1_mirror', 'hole2_virtual_carry'):
        r = prove_qasm(os.path.join(audit, name + '.qasm'), 'OPENQASM 3.0;\nqubit[1] q;\n')
        assert r['verdict'] != 'PROVED'


# ----------------------------------------------------------------------------- translation from Qiskit vs Qiskit itself
def test_qiskit_translation_matches_qiskit_operator():
    """Random Qiskit circuits (angles in (-4 pi, 4 pi), library blocks): Qiskit's own Operator of the original circuit
    equals the dense unitary of the translated DCircuit, global phase included (external check of qiskit_io, of the
    angle reading and of sim.py)."""
    qiskit = pytest.importorskip('qiskit')
    import numpy as np
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    from qiskit.circuit.library import QFTGate, MCPhaseGate, MCXGate
    from qoptct.verify.dyadic.qiskit_io import from_qiskit
    rng = random.Random(4)
    one = ['h', 'x', 'y', 'z', 's', 'sdg', 't', 'tdg', 'sx', 'sxdg']
    par1 = ['p', 'rz', 'rx', 'ry']
    for t in range(80):
        n = rng.randint(1, 4)
        qc = QuantumCircuit(n, global_phase=rng.choice([0, math.pi / 8, -math.pi / 2]))
        for _ in range(rng.randint(1, 12)):
            ang = math.pi * rng.randrange(-31, 32) / 8
            k = rng.random()
            if k < 0.3:
                getattr(qc, rng.choice(one))(rng.randrange(n))
            elif k < 0.6:
                getattr(qc, rng.choice(par1))(ang, rng.randrange(n))
            elif k < 0.7:
                qc.u(ang, ang / 2, -ang / 4, rng.randrange(n))
            elif n >= 2 and k < 0.85:
                a, b = rng.sample(range(n), 2)
                g = rng.choice(['cx', 'cz', 'cy', 'swap', 'cp', 'crz', 'rzz', 'cs', 'csx'])
                if g in ('cp', 'crz', 'rzz'):
                    getattr(qc, g)(ang, a, b)
                else:
                    getattr(qc, g)(a, b)
            elif n >= 3 and k < 0.92:
                qs = rng.sample(range(n), 3)
                rng.choice([lambda: qc.ccx(*qs), lambda: qc.cswap(*qs), lambda: qc.ccz(*qs),
                            lambda: qc.append(MCPhaseGate(ang, 2), qs),
                            lambda: qc.append(MCXGate(2, ctrl_state=rng.randrange(4)), qs)])()
            else:
                m = rng.randint(1, n)
                qs = rng.sample(range(n), m)
                qc.append(QFTGate(m) if rng.random() < 0.5 else QFTGate(m).inverse(), qs)
        c = from_qiskit(qc)
        U = Operator(qc).data
        V = sim.unitary(c)
        assert np.max(np.abs(U - V)) < 1e-9, (t, qc)
