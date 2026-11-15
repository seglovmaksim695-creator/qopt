# SPDX-License-Identifier: Apache-2.0
"""Tests of qoptct.crypto.fold: reduction modulo special primes and constant addition.

Levels: (1) exact simulation of the gate list (andcirc.simulate) on ALL inputs for small analogues p = 2^n - c and
on edge/random inputs at full size; (2) the emitted OpenQASM through qoptct.verify.arith.check_sim (independent
simulator, random measurement branches); (3) kprove (all inputs, all outcomes) for small n when it is available
(qoptct.verify.kprove from v0.5, or $QOPT_KPROVE_DIR with kprove.py); (4) canaries; (5) AND counts at full size
without simulation (pinned to the lab claims)."""
import os
import random
import re
import sys

import numpy as np
import pytest

from qoptct.crypto import fold as F
from qoptct.crypto.andcirc import Circ, text_counts


def small_cases(nmax):
    for n in range(3, nmax + 1):
        for c in range(1, 1 << n):
            if c * c + 2 * c + 2 < (1 << n):
                yield n, c


# ------------------------------------------------------------------ API
def test_api_errors():
    with pytest.raises(ValueError):
        F.reduce_special((1 << 8) - 21)               # c^2 + 2c + 2 >= 2^n: one fold is not enough
    with pytest.raises(ValueError):
        F.reduce_special(F.P25519, method='solinas')
    with pytest.raises(ValueError):
        F.reduce_special(F.P25519, method='mersenne')
    with pytest.raises(ValueError):
        F.reduce_special(F.P25519, method='montgomery')
    with pytest.raises(ValueError):
        F.reduce_special(F.P25519, step2='x')
    with pytest.raises(ValueError):
        F.add_const(4, 16)


def test_auto_method_and_layout():
    exp = {'curve25519': ('fold', 'M', 255), 'secp256k1': ('fold', 'M', 256), 'p256': ('solinas', 'direct', 256),
           'p384': ('fold', 'M', 384), 'p521': ('mersenne', 'direct', 521)}
    for name, p in F.CURVES.items():
        R = F.reduce_special(p)
        assert (R.method, R.step2, R.n) == exp[name]
        n = R.n
        assert R.layout['r'] == list(range(n)) and R.layout['H'] == list(range(n, 2 * n))
        regs = [q for k in R.layout if k not in ('P', 'r', 'H') for q in R.layout[k]]
        assert regs == list(range(2 * n, R.nref))


def test_signed_digits():
    for c in range(1, 3000):
        d = F.signed_digits(c)
        assert sum(sg << s for s, sg in d) == c
        assert len(d) <= bin(c).count('1')


def test_p256_formula():
    rng = random.Random(5)
    for _ in range(300):
        P = rng.getrandbits(512)
        t = F.p256_formula(P)
        assert (t - P) % F.P256 == 0 and -4 * F.P256 <= t < 5 * F.P256


# ------------------------------------------------------------------ exact simulation, all inputs (small n)
@pytest.mark.parametrize('n', [3, 4, 5, 6, 7, 8])
def test_fold_all_inputs_all_c(n):
    """every c with c^2 + 2c + 2 < 2^n, both step2 variants, ALL P < 2^(2n)."""
    allP = np.arange(1 << (2 * n), dtype=np.int64)
    k = 0
    for nn, c in small_cases(n):
        if nn != n:
            continue
        for st in ('M', 'direct'):
            R = F.reduce_special((1 << n) - c, n=n, step2=st)
            r = R.check(allP)
            assert r['ok'], (n, c, st, r)
            k += 1
    assert k > 0


@pytest.mark.parametrize('n,c', [(9, 1), (9, 3), (9, 19), (9, 21), (10, 1), (10, 3), (10, 19), (10, 29)])
def test_fold_all_inputs_n9_n10(n, c):
    allP = np.arange(1 << (2 * n), dtype=np.int64)
    for st in ('M', 'direct'):
        R = F.reduce_special((1 << n) - c, n=n, step2=st)
        assert R.check(allP)['ok'], (n, c, st)


def test_fold_edge_inputs_mid_size():
    for n, c in [(16, 19), (20, 19), (24, 5), (32, 977), (48, (1 << 20) + 977)]:
        p = (1 << n) - c
        R = F.reduce_special(p, n=n)
        assert R.check(F.edge_inputs(p, n, 400))['ok'], (n, c)


# ------------------------------------------------------------------ full size
def p256_extreme_inputs():
    """inputs whose words are 0 or 2^32-1 that maximise / minimise the NIST formula (it is linear in the words)."""
    S, coef = F.p256_terms(); wgt = [0] * 16
    for k, words in S.items():
        for pos, idx in enumerate(reversed(words)):
            if idx is not None:
                wgt[idx] += coef[k] << (32 * pos)
    hi = sum(0xffffffff << (32 * i) for i in range(16) if wgt[i] > 0)
    lo = sum(0xffffffff << (32 * i) for i in range(16) if wgt[i] < 0)
    return [hi, lo, hi ^ 1, lo ^ (1 << 300)]


@pytest.mark.parametrize('name', list(F.CURVES))
def test_full_size_edge_inputs(name):
    p = F.CURVES[name]
    R = F.reduce_special(p)
    vals = F.edge_inputs(p, R.n, 96, seed=7)
    if name == 'p256':
        ex = p256_extreme_inputs()
        assert max(F.p256_formula(x) for x in ex) > 4 * p and min(F.p256_formula(x) for x in ex) < -3 * p
        vals += ex
    r = R.check(vals)
    assert r['ok'], r


# ------------------------------------------------------------------ AND counts (no simulation)
FULL = {'curve25519': 1304, 'secp256k1': 2175, 'p256': 2744, 'p384': 3333, 'p521': 1562}


@pytest.mark.parametrize('name', list(F.CURVES))
def test_full_size_and_counts(name):
    R = F.reduce_special(F.CURVES[name])
    assert R.n_and == FULL[name]
    tc = R.text_counts()
    assert tc['T'] == 4 * R.n_and and tc['AND'] == R.n_and
    # every AND is erased by one measurement, except the carry-out AND that writes g (step2 'M'; g is garbage)
    assert tc['measurements'] == R.n_and - (1 if R.step2 == 'M' else 0)
    cm = R.counts()
    assert cm['T'] == 4 * R.n_and


def test_and_count_formulas():
    for n in list(range(100, 522, 23)) + [255, 256, 521]:
        assert F.fold_and_count((1 << n) - 19) == 5 * n + 29, n
        assert F.fold_and_count((1 << n) - (1 << 32) - 977) == 7 * n + 383, n
        assert F.fold_and_count((1 << n) - 1) == 3 * n - 1, n


# ------------------------------------------------------------------ gate level (emitted QASM)
def _check_sim(R, vals, tmp_path, txt=None, m_trials=3, complete=False):
    from qoptct.verify import arith
    n = R.n; nref = R.nref
    path = tmp_path / ('red_%d.qasm' % random.randrange(1 << 30))
    path.write_text(txt if txt is not None else R.qasm())
    X = np.zeros((len(vals), nref), np.uint8)
    for i, P in enumerate(vals):
        for j in range(2 * n):
            X[i, j] = (int(P) >> j) & 1

    def f(Xb):
        D = np.zeros((len(Xb), nref), np.uint8)
        for row in range(len(Xb)):
            P = 0
            for j in range(2 * n):
                P |= int(Xb[row, j]) << j
            ref = R.reference(P)
            for reg, ws in R.layout.items():
                if reg == 'P':
                    continue
                val = ref['r' if reg == 'r' else reg]
                for j, q in enumerate(ws):
                    D[row, q] = (val >> j) & 1
        return D
    return arith.check_sim(str(path), (nref, f, list(range(2 * n))), X=X, domain_complete=complete, m_trials=m_trials)


def test_qasm_check_sim_small_all_inputs(tmp_path):
    for n, c, st in [(5, 3, 'M'), (6, 3, 'direct'), (6, 5, 'M')]:
        R = F.reduce_special((1 << n) - c, n=n, step2=st)
        r = _check_sim(R, list(range(1 << (2 * n))), tmp_path, complete=True)
        assert r['ok'] and r['proved_all_x'], (n, c, st, r.get('error'), r.get('worst'))


@pytest.mark.parametrize('name', list(F.CURVES))
def test_qasm_check_sim_full_size(name, tmp_path):
    p = F.CURVES[name]
    R = F.reduce_special(p)
    vals = F.edge_inputs(p, R.n, 64, seed=3) + (p256_extreme_inputs() if name == 'p256' else [])
    r = _check_sim(R, vals, tmp_path, m_trials=3)
    assert r['ok'], (r.get('error'), r.get('worst'))


# ------------------------------------------------------------------ canaries
def _mutants(txt, rng):
    L = txt.split('\n'); out = []
    ifs = [i for i, l in enumerate(L) if l.startswith('if') and 'cz' in l]
    ts = [i for i, l in enumerate(L) if l.startswith('t ')]
    cxs = [i for i, l in enumerate(L) if l.startswith('cx ')]
    i = rng.choice(ts); M = list(L); M[i] = 'tdg' + M[i][1:]; out.append(('t_to_tdg', '\n'.join(M)))
    i = rng.choice(cxs); M = list(L); del M[i]; out.append(('drop_cx', '\n'.join(M)))
    i = rng.choice(ifs); M = list(L); M[i] = re.sub(r'cz q\[\d+\], q\[\d+\]; ?', '', M[i], count=1)
    out.append(('drop_corr_cz', '\n'.join(M)))
    return out


def test_canaries_gate_level(tmp_path):
    n, c = 5, 3
    R = F.reduce_special((1 << n) - c, n=n)
    rng = random.Random(11)
    allx = list(range(1 << (2 * n)))
    for kind, mtxt in _mutants(R.qasm(), rng):
        r = _check_sim(R, allx, tmp_path, txt=mtxt, m_trials=8, complete=True)
        assert not r['ok'], kind


def test_canaries_op_level():
    """mutated gate lists must fail the exact simulation on all inputs."""
    n, c = 6, 5
    allP = np.arange(1 << (2 * n), dtype=np.int64)
    R = F.reduce_special((1 << n) - c, n=n)
    assert R.check(allP)['ok']
    ops = R.circ.ops
    kinds = {'cz': [j for j, o in enumerate(ops) if o[0] == 'MZ' and o[2]],
             'cx': [j for j, o in enumerate(ops) if o[0] == 'cx'],
             'x': [j for j, o in enumerate(ops) if o[0] == 'x']}
    rng = random.Random(2)
    for kind, idx in kinds.items():
        caught = 0
        for j in rng.sample(idx, min(6, len(idx))):
            M = list(ops)
            if kind == 'cz':
                M[j] = ('MZ', ops[j][1], ())          # drop the correction of an erased AND
            else:
                del M[j]
            R.circ.ops = M
            try:
                caught += not R.check(allP)['ok']
            finally:
                R.circ.ops = ops
        assert caught >= 1, kind


# ------------------------------------------------------------------ kprove (all inputs, all outcomes)
def _kprove():
    try:
        from qoptct.verify.kprove import kprove
        return kprove
    except ImportError:
        pass
    d = os.environ.get('QOPT_KPROVE_DIR')
    if d and os.path.exists(os.path.join(d, 'kprove.py')):
        if d not in sys.path:
            sys.path.insert(0, d)
        from kprove import kprove
        return kprove
    pytest.skip('kprove not available (qoptct.verify.kprove or $QOPT_KPROVE_DIR)')


def spec_fold(n, p):
    """independent word-level spec: r = P mod p by two folds with k = 2^n mod p and one conditional subtraction;
    H unchanged; garbage wires are free."""
    k = (1 << n) % p

    def spec(S):
        P = S.reg(range(2 * n)); L, H = P.resize(n), P.shr(n)
        wt = (((1 << n) - 1) * (k + 1)).bit_length() + 1
        t = (L.resize(wt) + H.mul_const(k, wt)).resize(wt)
        t2 = (t.resize(n).resize(n + 2) + t.shr(n).mul_const(k, n + 2)).resize(n + 2)
        r = t2.mod_reduce(p, n + 1)
        return {**S.out(range(n), r), **S.out(range(n, 2 * n), H)}
    return spec


@pytest.mark.parametrize('n,c', [(6, 3), (8, 5), (11, 19)])
def test_kprove_fold(n, c, tmp_path):
    kprove = _kprove()
    p = (1 << n) - c
    R = F.reduce_special(p, n=n)
    path = tmp_path / 'red.qasm'; path.write_text(R.qasm())
    garb = range(2 * n, R.nref)
    r = kprove(str(path), spec=spec_fold(n, p), n=R.nref, clean_wires=garb, free_wires=garb, budget_sec=300)
    assert r['status'] == 'PROVED', r.get('detail')
    if n == 8:
        rng = random.Random(4)
        for kind, mtxt in _mutants(R.qasm(), rng):
            mp = tmp_path / ('mut_%s.qasm' % kind); mp.write_text(mtxt)
            rm = kprove(str(mp), spec=spec_fold(n, p), n=R.nref, clean_wires=garb, free_wires=garb, budget_sec=120)
            assert rm['status'] != 'PROVED', kind


@pytest.mark.parametrize('n', [16, 64])
def test_kprove_const_add(n, tmp_path):
    kprove = _kprove()
    rng = random.Random(5)
    for K in (rng.getrandbits(n) | 1, (rng.getrandbits(n) << 3) % (1 << n) | 8):
        for ctl in (False, True):
            A = F.controlled_add_const(n, K) if ctl else F.add_const(n, K)
            path = tmp_path / 'ca.qasm'; path.write_text(A.qasm())

            def spec(S, K=K, ctl=ctl):
                X = S.reg(range(n)); s = (X + K).resize(n)
                return S.out(range(n), S.mux(S.x(n), s, X) if ctl else s)
            r = kprove(str(path), spec=spec, n=A.nref, budget_sec=120)
            assert r['status'] == 'PROVED', (n, K, ctl, r.get('detail'))


# ------------------------------------------------------------------ constant addition
@pytest.mark.parametrize('n', [1, 2, 3, 4, 5, 6])
def test_const_add_all(n):
    xs = list(range(1 << n))
    for K in range(1 << n):
        A = F.add_const(n, K)
        assert A.n_and == A.expected_and() == (0 if K == 0 else max(0, n - 2 - F.tz(K)))
        assert A.check(xs)['ok'], (n, K)
        B = F.controlled_add_const(n, K)
        assert B.n_and == B.expected_and() == (0 if K == 0 else max(0, n - 1 - F.tz(K)))
        assert B.check(xs + xs, [0] * len(xs) + [1] * len(xs))['ok'], (n, K)


def test_const_add_counts_large():
    rng = random.Random(9)
    for n in (8, 64, 256, 521):
        for K in [rng.getrandbits(n) | 1, 1 << (n - 1), (rng.getrandbits(n) << 3) % (1 << n) | 8]:
            A, B = F.add_const(n, K), F.controlled_add_const(n, K)
            assert A.n_and == max(0, n - 2 - F.tz(K)) and B.n_and == max(0, n - 1 - F.tz(K))
            assert text_counts(A.qasm())['AND'] == A.n_and
        A = F.add_const(n, K)
        xs = [0, 1, (1 << n) - 1, (1 << n) - K] + [rng.getrandbits(n) for _ in range(20)]
        assert A.check([x % (1 << n) for x in xs])['ok']
