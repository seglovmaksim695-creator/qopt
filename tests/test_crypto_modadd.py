"""Tests of qoptct.crypto.modadd: modular adders / subtractors, doubling, comparators.

Levels: (1) exact simulation of the gate list (andcirc.simulate) on EVERY valid input for all small moduli; (2) the
emitted OpenQASM through qoptct.verify.arith.check_sim (independent simulator, random measurement branches); (3) kprove (all
inputs of the domain, all outcomes) against word-level specifications; (4) canaries; (5) AND counts pinned to the lab
(agent a3: 4n-2, 3n-2, 3n-1, 2n-2 ...)."""
import os
import random
import re
import sys

import numpy as np
import pytest

from qoptct.crypto import modadd as M


# ------------------------------------------------------------------ API
def test_api_errors():
    with pytest.raises(ValueError):
        M.mod_add(1)
    with pytest.raises(ValueError):
        M.mod_add(20, n=4)                       # p > 2^n
    with pytest.raises(ValueError):
        M.mod_add_out(8)
    with pytest.raises(ValueError):
        M.mod_double(8)
    with pytest.raises(ValueError):
        M.mod_add_const(13, 0)
    with pytest.raises(ValueError):
        M.mod_add_const(13, 13)
    with pytest.raises(ValueError):
        M.mod_add_mersenne(1)
    with pytest.raises(ValueError):
        M.compare_geq(0)
    with pytest.raises(ValueError):
        M.compare_const_geq(3, 9)


def test_pure_deterministic():
    assert M.mod_add(13).qasm() == M.mod_add(13).qasm()
    assert M.mod_sub(13).qasm() == M.mod_sub(13).qasm()


# ------------------------------------------------------------------ exact simulation, every valid input
def test_mod_add_sub_all_moduli():
    for n in (2, 3, 4, 5):
        for p in range(2, (1 << n) + 1):
            if (p - 1).bit_length() != n and p > 2:
                continue                           # exactly n bits, plus p = 2^n handled below
            for build in (M.mod_add, M.mod_sub):
                A = build(p, n)
                r = A.check(A.domain())
                assert r['ok'], (build.__name__, n, p, r)


def test_wider_register_than_needed():
    """p < 2^(n-1): the registers may be wider than needed."""
    for p, n in [(5, 5), (11, 6), (3, 4)]:
        A = M.mod_add(p, n)
        assert A.check(A.domain())['ok']
        A = M.mod_sub(p, n)
        assert A.check(A.domain())['ok']


def test_mod_add_out_all_odd_moduli():
    for n in (2, 3, 4, 5):
        for p in range(3, (1 << n) + 1, 2):
            if (p - 1).bit_length() != n:
                continue
            A = M.mod_add_out(p, n)
            r = A.check(A.domain())
            assert r['ok'], (n, p, r)


def test_mersenne_all():
    for n in range(2, 8):
        A = M.mod_add_mersenne(n)
        assert A.check(A.domain())['ok'], n


def test_const_all_a_small():
    for n in (2, 3, 4):
        for p in range(2, (1 << n) + 1):
            if (p - 1).bit_length() != n and p > 2:
                continue
            for a in range(1, p):
                try:
                    A = M.mod_add_const(p, a, n)
                except ValueError:
                    continue
                assert A.check(A.domain())['ok'], (n, p, a)
                B = M.mod_sub_const(p, a, n)
                assert B.check(B.domain())['ok'], (n, p, a)
    # p = 5, n = 3: every a
    for a in range(1, 5):
        assert M.mod_add_const(5, a, 3).check(M.mod_add_const(5, a, 3).domain())['ok']


def test_double_all_odd_moduli():
    for n in (2, 3, 4, 5):
        for p in range(3, (1 << n) + 1, 2):
            if (p - 1).bit_length() != n:
                continue
            A = M.mod_double(p, n)
            r = A.check(A.domain())
            assert r['ok'], (n, p, r)


def test_comparators_all():
    for n in (1, 2, 3, 4):
        A = M.compare_geq(n)
        assert A.check(A.domain())['ok'], n
        for c in range(0, (1 << n) + 1):
            B = M.compare_const_geq(n, c)
            assert B.check(B.domain())['ok'], (n, c)


def test_full_size_random():
    """n = 64, 256: random and edge inputs."""
    rng = random.Random(3)
    for n in (64, 256):
        p = (1 << n) - rng.choice([189, 237, 357]) if n == 256 else (1 << 64) - 59
        pts = [dict(x=a, y=b) for a, b in [(0, 0), (p - 1, p - 1), (p - 1, 1), (1, p - 1), (0, p - 1), (p // 2, p // 2 + 1)]]
        pts += [dict(x=rng.randrange(p), y=rng.randrange(p)) for _ in range(150)]
        for build in (M.mod_add, M.mod_sub, M.mod_add_out):
            A = build(p, n)
            r = A.check(pts)
            assert r['ok'], (build.__name__, n, r)
        a = rng.randrange(1, p)
        pts1 = [dict(y=pt['y']) for pt in pts]
        for build in (M.mod_add_const, M.mod_sub_const):
            r = build(p, a, n).check(pts1)
            assert r['ok'], (build.__name__, n, r)
        r = M.mod_double(p, n).check(pts1)
        assert r['ok'], ('double', n, r)
    n = 255
    A = M.mod_add_mersenne(n)
    p = (1 << n) - 1
    pts = [dict(x=rng.randrange(p), y=rng.randrange(p)) for _ in range(100)] + [dict(x=p - 1, y=p - 1), dict(x=0, y=0)]
    assert A.check(pts)['ok']


# ------------------------------------------------------------------ AND counts (lab a3)
def test_and_counts_pinned():
    for n, k in [(4, 14), (5, 18), (6, 22), (8, 30)]:
        p = (1 << n) - 3
        assert M.mod_add(p, n).n_and == M.mod_sub(p, n).n_and == k == 4 * n - 2
    for n in (4, 5, 6, 8):
        p = (1 << n) - 3
        assert M.mod_add_out(p, n).n_and == 3 * n - 2
        assert M.mod_add_mersenne(n).n_and == 3 * n - 1
        assert M.mod_double(p, n).n_and == 2 * n - 2
    assert M.mod_add_const(29, 7, 5).n_and == 11 and M.mod_add_const(61, 10, 6).n_and == 14
    assert M.mod_add_const(251, 100, 8).n_and == 19
    assert M.mod_sub_const(251, 100, 8).n_and == 19
    # large
    p = (1 << 256) - 189
    assert M.mod_add(p, 256).n_and == 4 * 256 - 2 and M.mod_add_out(p, 256).n_and == 3 * 256 - 2
    for n in (1, 2, 8, 64):
        assert M.compare_geq(n).n_and == n
    assert M.compare_const_geq(4, 1).n_and == 3 and M.compare_const_geq(4, 8).n_and == 0
    c = M.mod_add(29, 5).counts()
    assert c['AND'] == 18 and c['T'] == 72


def test_layout():
    A = M.mod_add(13, 4)
    assert A.layout['inputs'] == dict(x=[0, 1, 2, 3], y=[4, 5, 6, 7]) and A.layout['outputs'] == dict(y=[4, 5, 6, 7])
    Z = M.mod_add_out(13, 4)
    assert Z.layout['outputs']['z'] == [8, 9, 10, 11] and Z.nref == 12
    D = M.mod_double(13, 4)
    assert D.layout['outputs']['z'][1:] == [0, 1, 2] and D.nref == 5


# ------------------------------------------------------------------ gate level (emitted QASM)
def _check_sim(A, pts, tmp_path, txt=None, m_trials=4):
    from qoptct.verify import arith
    nref = A.nref
    path = tmp_path / ('ma_%d.qasm' % random.randrange(1 << 30))
    path.write_text(txt if txt is not None else A.qasm())
    Lin, Lout = A.layout['inputs'], A.layout['outputs']
    X = np.zeros((len(pts), nref), np.uint8)
    for i, pt in enumerate(pts):
        for name, ws in Lin.items():
            for j, w in enumerate(ws):
                X[i, w] = (pt[name] >> j) & 1

    def f(Xb):
        D = np.zeros((len(Xb), nref), np.uint8)
        for row in range(len(Xb)):
            pt = {name: sum(int(Xb[row, w]) << j for j, w in enumerate(ws)) for name, ws in Lin.items()}
            ref = A.reference(pt)
            for name, ws in Lin.items():
                if name not in Lout and name not in A._consumed:
                    for j, w in enumerate(ws):
                        D[row, w] = (pt[name] >> j) & 1
            for name, ws in Lout.items():
                for j, w in enumerate(ws):
                    D[row, w] = (ref[name] >> j) & 1
        return D
    return arith.check_sim(str(path), (nref, f, sorted({w for ws in Lin.values() for w in ws})), X=X,
                           domain_complete=False, m_trials=m_trials)


BUILDS = {
    'add': lambda: M.mod_add(13, 4), 'sub': lambda: M.mod_sub(13, 4), 'out': lambda: M.mod_add_out(13, 4),
    'mers': lambda: M.mod_add_mersenne(4), 'const': lambda: M.mod_add_const(13, 5, 4),
    'subc': lambda: M.mod_sub_const(13, 5, 4), 'dbl': lambda: M.mod_double(13, 4),
    'geq': lambda: M.compare_geq(3), 'cgeq': lambda: M.compare_const_geq(4, 6),
}


@pytest.mark.parametrize('name', list(BUILDS))
def test_qasm_check_sim_small(name, tmp_path):
    A = BUILDS[name]()
    r = _check_sim(A, A.domain(), tmp_path)
    assert r['ok'], (name, r.get('error'), r.get('worst'))


def test_qasm_check_sim_full_size(tmp_path):
    rng = random.Random(8)
    n = 64
    p = (1 << 64) - 59
    pts = [dict(x=rng.randrange(p), y=rng.randrange(p)) for _ in range(24)] + [dict(x=p - 1, y=p - 1)]
    for build in (M.mod_add, M.mod_sub):
        A = build(p, n)
        r = _check_sim(A, pts, tmp_path, m_trials=2)
        assert r['ok'], (build.__name__, r.get('error'))


# ------------------------------------------------------------------ canaries
def _mutants(txt, rng):
    L = txt.split('\n'); out = []
    ifs = [i for i, l in enumerate(L) if l.startswith('if') and 'cz' in l]
    ts = [i for i, l in enumerate(L) if l.startswith('t ')]
    cxs = [i for i, l in enumerate(L) if l.startswith('cx ')]
    i = rng.choice(ts); M_ = list(L); M_[i] = 'tdg' + M_[i][1:]; out.append(('t_to_tdg', '\n'.join(M_)))
    i = rng.choice(cxs); M_ = list(L); del M_[i]; out.append(('drop_cx', '\n'.join(M_)))
    i = ifs[0]; M_ = list(L); M_[i] = re.sub(r'cz q\[\d+\], q\[\d+\]; ?', '', M_[i], count=1)
    out.append(('drop_cz', '\n'.join(M_)))
    return out


@pytest.mark.parametrize('name', ['add', 'sub', 'dbl'])
def test_canaries_gate_level(name, tmp_path):
    A = BUILDS[name]()
    rng = random.Random(5)
    for kind, mtxt in _mutants(A.qasm(), rng):
        r = _check_sim(A, A.domain(), tmp_path, txt=mtxt, m_trials=8)
        assert not r['ok'], (name, kind)


def test_canaries_op_level():
    A = M.mod_add(13, 4)
    pts = A.domain()
    assert A.check(pts)['ok']
    ops = A.circ.ops
    kinds = {'corr': [j for j, o in enumerate(ops) if o[0] == 'MZ' and o[2]],
             'cx': [j for j, o in enumerate(ops) if o[0] == 'cx']}
    rng = random.Random(2)
    for kind, idx in kinds.items():
        caught = 0
        for j in rng.sample(idx, min(8, len(idx))):
            m = list(ops)
            if kind == 'corr':
                m[j] = ('MZ', ops[j][1], ())
            else:
                del m[j]
            A.circ.ops = m
            try:
                caught += not A.check(pts)['ok']
            finally:
                A.circ.ops = ops
        assert caught >= 1, kind


# ------------------------------------------------------------------ kprove (all valid inputs, all outcomes)
def _kprove():
    try:
        from qoptct.verify.kprove import kprove, reg_value_lt
        return kprove, reg_value_lt
    except ImportError:
        pass
    d = os.environ.get('QOPT_KPROVE_DIR')
    if d and os.path.exists(os.path.join(d, 'kprove.py')):
        if d not in sys.path:
            sys.path.insert(0, d)
        from kprove import kprove, reg_value_lt
        return kprove, reg_value_lt
    pytest.skip('kprove not available (qoptct.verify.kprove or $QOPT_KPROVE_DIR)')


def _pre(A, reg_value_lt, regs):
    """promise: every register in `regs` is < p."""
    names = [''] * A.nref
    for name, ws in A.layout['inputs'].items():
        for j, w in enumerate(ws):
            names[w] = '%s_%d' % (name, j)
    parts = [reg_value_lt(names, name + '_', A.p) for name in regs]

    def build(ps, xedge):
        r = None
        for b, _ in parts:
            e = b(ps, xedge)
            r = e if r is None else ps.B.AND(r, e)
        return r

    def sample(x, rng):
        for _, s in parts:
            x = s(x, rng)
        return x
    return build, sample


def _spec(A):
    """independent word-level specifications."""
    p, n, Li, Lo = A.p, A.n, A.layout['inputs'], A.layout['outputs']
    k = A.kind
    from qoptct.verify.bdd import ZERO

    def spec(S):
        out = {}
        if k in ('mod_add', 'mod_sub', 'mod_add_out', 'mod_add_mersenne'):
            X, Y = S.reg(Li['x']), S.reg(Li['y'])
            if k == 'mod_sub':
                d = Y.sub(X)                                    # n bits + borrow
                r = S.mux(d.bits[-1], (d.resize(n) + p).resize(n), d.resize(n))
            else:
                r = (X + Y).mod_reduce(p, n)
            out.update(S.out(Li['x'], X))
            if k == 'mod_add_out':
                out.update(S.out(Li['y'], Y)); out.update(S.out(Lo['z'], r))
            else:
                out.update(S.out(Lo['y'], r))
        elif k in ('mod_add_const', 'mod_sub_const'):
            a = A.a
            Y = S.reg(Li['y'])
            if k == 'mod_add_const':
                r = (Y + a).mod_reduce(p, n)
            else:
                d = Y.sub(a)
                r = S.mux(d.bits[-1], (d.resize(n) + p).resize(n), d.resize(n))
            out.update(S.out(Lo['y'], r))
        elif k == 'mod_double':
            Y = S.reg(Li['y'])
            r = Y.shl(1).mod_reduce(p, n)
            for w in Li['y']:
                out[w] = ZERO
            out.update(S.out(Lo['z'], r))
        elif k == 'compare_geq':
            X, Y = S.reg(Li['x']), S.reg(Li['y'])
            out.update(S.out(Li['x'], X)); out.update(S.out(Li['y'], Y))
            out[Li['t'][0]] = S.B.XOR(S.x(Li['t'][0]), X >= Y)
        for w in range(A.circ.N):          # every other wire (ancillas) ends at 0
            out.setdefault(w, ZERO)
        return out
    return spec


def _kp(A, tmp_path, kprove, reg_value_lt, txt=None, budget=300, max_mb=None):
    path = tmp_path / 'm.qasm'; path.write_text(txt if txt is not None else A.qasm())
    Li, Lo = A.layout['inputs'], A.layout['outputs']
    # ALL wires are listed as reference wires (ancillas too): kprove checks the ancillas of a circuit over every x, not only
    # over the promise; as data wires they are compared to 0 on the promise only
    clean = sorted(set(range(A.circ.N)) - {w for ws in Li.values() for w in ws})
    regs = [r for r in Li if r in ('x', 'y')]
    kw = {} if max_mb is None else dict(max_mb=max_mb)
    pre = None if A.kind == 'compare_geq' else _pre(A, reg_value_lt, regs)
    return kprove(str(path), spec=_spec(A), n=A.circ.N, clean_wires=clean, pre=pre, budget_sec=budget, **kw)


KP_CASES = [
    ('add', lambda: M.mod_add(13, 4)), ('add2', lambda: M.mod_add(29, 5)), ('sub', lambda: M.mod_sub(13, 4)),
    ('out', lambda: M.mod_add_out(13, 4)), ('mers', lambda: M.mod_add_mersenne(5)),
    ('const', lambda: M.mod_add_const(13, 5, 4)), ('subc', lambda: M.mod_sub_const(13, 5, 4)),
    ('dbl', lambda: M.mod_double(13, 4)), ('geq', lambda: M.compare_geq(4)),
]


@pytest.mark.parametrize('name,build', KP_CASES, ids=[c[0] for c in KP_CASES])
def test_kprove_small(name, build, tmp_path):
    kprove, rvl = _kprove()
    A = build()
    r = _kp(A, tmp_path, kprove, rvl)
    assert r['status'] == 'PROVED', (name, r.get('detail'))


@pytest.mark.slow
@pytest.mark.parametrize('name,build', [('add8', lambda: M.mod_add(251, 8)), ('sub8', lambda: M.mod_sub(251, 8)),
                                        ('out8', lambda: M.mod_add_out(251, 8)),
                                        ('add16', lambda: M.mod_add((1 << 16) - 15, 16)),
                                        ('add32', lambda: M.mod_add((1 << 32) - 5, 32)),
                                        ('sub32', lambda: M.mod_sub((1 << 32) - 5, 32)),
                                        ('add64', lambda: M.mod_add((1 << 64) - 59, 64)),
                                        ('mers64', lambda: M.mod_add_mersenne(64)),
                                        ('const32', lambda: M.mod_add_const((1 << 32) - 5, 12345, 32)),
                                        ('dbl32', lambda: M.mod_double((1 << 32) - 5, 32))],
                         ids=['add8', 'sub8', 'out8', 'add16', 'add32', 'sub32', 'add64', 'mers64', 'const32', 'dbl32'])
def test_kprove_bigger(name, build, tmp_path):
    kprove, rvl = _kprove()
    r = _kp(build(), tmp_path, kprove, rvl, max_mb=1400)
    assert r['status'] == 'PROVED', (name, r.get('detail'), r.get('mem'))


def test_kprove_mutants(tmp_path):
    kprove, rvl = _kprove()
    A = M.mod_add(13, 4)
    assert _kp(A, tmp_path, kprove, rvl)['status'] == 'PROVED'
    rng = random.Random(4)
    for kind, mtxt in _mutants(A.qasm(), rng):
        r = _kp(A, tmp_path, kprove, rvl, txt=mtxt, budget=120)
        assert r['status'] == 'REFUTED', (kind, r['status'])
