# SPDX-License-Identifier: Apache-2.0
"""Tests of qoptct.crypto.lookup (narrow table lookups z ^= T[a] on logical ANDs).

Every builder already checks its circuit exactly on all addresses (internal batch simulation).  The tests add an
INDEPENDENT check of the emitted OpenQASM text against the Toffoli reference of reference_qc:
qoptct.verify.arith.check_sim (other simulator, all addresses x several z, random outcome branches) and, for small
sizes (c <= 6), the exact prover qoptct.verify.prove (and qoptct.verify.kprove when present, qoptct >= 0.5);
mutated circuits (canaries) must be rejected.  Slow test (-m slow): the AND-count distribution on 40 random tables."""
import random

import numpy as np
import pytest

from qoptct.crypto import lookup as L


def table(c, w, seed):
    rng = random.Random(seed)
    return [rng.getrandbits(w) for _ in range(1 << c)]


def n_ref_wires(qc_text):
    return len(qc_text.split('\n', 1)[0].split()) - 1


def sim_check(r, tmp_path, zs=(0, 1, 5), m_trials=4, qasm=None, tag='c'):
    """independent simulation of the emitted text: all addresses x z values (z / acc wires), other wires 0."""
    from qoptct.verify import arith as VA
    qp = tmp_path / (tag + '.qasm')
    qcp = tmp_path / (tag + '.qc')
    r.write(str(qp), str(qcp))
    if qasm is not None:
        qp.write_text(qasm)
    n = n_ref_wires(qcp.read_text())
    full = (1 << r.w) - 1
    zv = sorted({z & full for z in zs} | {full})
    X = np.zeros(((1 << r.c) * len(zv), n), np.uint8)
    row = 0
    for a in range(1 << r.c):
        for z in zv:
            for i in range(r.c):
                X[row, i] = (a >> i) & 1
            for i in range(r.w):
                X[row, r.c + i] = (z >> i) & 1
            row += 1
    res = VA.check_sim(str(qp), str(qcp), X=X, m_trials=m_trials, anc='zero_clean', seed=3)
    return bool(res.get('ok')), res


def prove(r, tmp_path, qasm=None, tag='p', budget=300):
    from qoptct import verify
    qp = tmp_path / (tag + '.qasm')
    qcp = tmp_path / (tag + '.qc')
    r.write(str(qp), str(qcp))
    if qasm is not None:
        qp.write_text(qasm)
    return verify.prove(str(qp), str(qcp), clean=True, budget_sec=budget)['status']


def canaries(qasm):
    """mutants: t -> tdg, drop the CZ of one erasure fix, drop one CX, drop one parity-conditioned block."""
    lines = qasm.splitlines()
    out = {}
    i = next((i for i, s in enumerate(lines) if s.startswith('t ')), None)
    if i is not None:
        out['t2tdg'] = lines[:i] + ['tdg ' + lines[i][2:]] + lines[i + 1:]
    i = next((i for i, s in enumerate(lines) if s.startswith('if (c[') and 'cz ' in s and '^' not in s), None)
    if i is not None:
        s = lines[i]
        j = s.index('cz ')
        out['dropCZfix'] = lines[:i] + [s[:j] + s[s.index(';', j) + 1:].lstrip()] + lines[i + 1:]
    idx = [i for i, s in enumerate(lines) if s.startswith('cx ')]
    if idx:
        i = idx[len(idx) // 2]
        out['dropCX'] = lines[:i] + lines[i + 1:]
    idx = [i for i, s in enumerate(lines) if s.startswith('if (') and '^' in s]
    if idx:
        i = idx[len(idx) // 2]
        out['dropParity'] = lines[:i] + lines[i + 1:]
    return {k: '\n'.join(v) + '\n' for k, v in out.items()}


# ------------------------------------------------------------------ AND counts + independent simulation
@pytest.mark.parametrize('c,w,seed,ands,h', [
    (4, 4, 1, 6, 2),        # base 2 + R 4          (Davio 11)
    (4, 8, 2, 11, None),    # core dim 9 -> 2 + 9 = 11 = Davio: falls back to Davio
    (5, 8, 3, None, 2),     # base 5 + R (R >= 8)    (Davio 26)
    (5, 16, 2, 21, 2),      # 5 + 16                 (Davio 26)
    (6, 21, 5, 33, 2),      # 12 + 21                (Davio 57; h=3 gives 36..40)
    (6, 32, 1, 40, 3),      # 8 + 32                 (Davio 57)
])
def test_narrow_counts_and_independent_sim(c, w, seed, ands, h, tmp_path):
    T = table(c, w, seed)
    r = L.narrow_lookup(T, c, w, seed=seed)
    assert r.checked.startswith('exact on all %d addresses' % (1 << c))
    assert r.n_and <= L.davio_ands(c)
    assert r.t_count == 4 * r.n_and
    if ands is not None:
        assert r.n_and == ands
    if h is not None:
        assert r.method == 'bilinear' and r.plan.h == h
        assert r.n_and == r.plan.base + r.plan.R
        assert r.n_and >= r.plan.base + r.plan.dim_core          # lower bound of the construction
    else:
        assert r.method == 'davio'
    ok, res = sim_check(r, tmp_path)
    assert ok, res


def test_narrow_davio_fallback_and_degenerate_tables(tmp_path):
    # wide output: the core has full rank 9 -> no gain over Davio (11), Davio is returned
    r = L.narrow_lookup(table(4, 16, 7), 4, 16)
    assert r.method == 'davio' and r.n_and == 11
    # zero / affine tables: no AND at all
    assert L.narrow_lookup([0] * 32, 5, 3).n_and == 0
    aff = [(a & 1) * 5 ^ ((a >> 2) & 1) * 3 ^ 6 for a in range(16)]
    r = L.narrow_lookup(aff, 4, 3)
    assert r.n_and == 0
    assert sim_check(r, tmp_path)[0]
    # one address bit
    r = L.narrow_lookup([2, 1], 1, 2)
    assert r.n_and == 0 and r.method == 'davio'
    assert sim_check(r, tmp_path, tag='one')[0]
    # a pure-hi table: Davio needs only the hi monomials, the split would build all of them
    hi = [((a >> 3) & (a >> 4) & (a >> 5) & 1) * 7 for a in range(64)]
    r = L.narrow_lookup(hi, 6, 3)
    assert r.n_and <= 2


def test_lower_bound_and_plan_fields():
    T = table(5, 6, 11)
    p = L.plan_narrow(T, 5, 6, seed=11)
    assert p is not None and p.c == 5 and p.h + p.l == 5
    assert p.n_and == L.base_ands(p.h, p.l) + p.R and p.R >= p.dim_core
    assert p.seconds >= 0
    assert L.plan_narrow(T, 5, 6, splits=[1]).h == 1


# ------------------------------------------------------------------ other constructions
@pytest.mark.parametrize('c,w', [(3, 2), (4, 3), (5, 8), (6, 21)])
def test_qroam1_formula_and_sim(c, w, tmp_path):
    T = table(c, w, 40 + c)
    r = L.qroam1_lookup(T, c, w)
    assert r.method == 'qroam1'
    assert r.n_and == L.qroam1_ands(c, w)
    if (c, w) == (6, 21):
        assert r.n_and == 52
    assert sim_check(r, tmp_path)[0]


def test_unlookup_counts():
    assert [L.unlookup_ands(c) for c in (2, 3, 4, 5, 6, 7)] == [0, 1, 2, 5, 8, 15]
    assert L.davio_ands(6) == 57 and L.base_ands(3, 3) == 8 and L.qroam1_ands(6, 21) == 52


@pytest.mark.parametrize('mode', ['narrow', 'davio', 'qroam1'])
def test_pair_small_sim(mode, tmp_path):
    T = table(4, 4, 21)
    r = L.lookup_pair(T, 4, 4, mode=mode)
    assert r.kind == 'pair'
    assert sim_check(r, tmp_path)[0]


def test_pair_c6_counts(tmp_path):
    T = table(6, 21, 13)
    q = L.lookup_pair(T, 6, 21, mode='qroam1')
    d = L.lookup_pair(T, 6, 21, mode='davio')
    n = L.lookup_pair(T, 6, 21, mode='narrow')
    assert d.n_and == 57 + 8
    assert q.n_and == 26 + 21 + 8 == 55
    assert n.method == 'bilinear' and n.n_and == n.plan.n_and + 8 and n.n_and < q.n_and
    assert sim_check(n, tmp_path, m_trials=3)[0]


def test_lookup_add_small(tmp_path):
    T = table(5, 8, 3)
    r = L.lookup_add(T, 5, 8, seed=3)
    assert r.kind == 'add' and r.ref_kind == 'add' and r.method == 'bilinear'
    assert r.n_and == r.plan.n_and + 7 + L.unlookup_ands(5)
    assert sim_check(r, tmp_path, zs=(0, 1, 200, 255))[0]
    # Davio fallback inside lookup_add (c=4, w=8: core rank 9 -> 2 + 9 = 11 = Davio)
    d = L.lookup_add(table(4, 8, 2), 4, 8, seed=2)
    assert d.method == 'davio' and d.n_and == 11 + 7 + 2


def test_lookup_add_loop1_step(tmp_path):
    """acc(32) += T[a(6)]: 40 + 31 + 8 = 79 ANDs (Gidney 2025 counting: 57 + 31 + phaseup)."""
    T = table(6, 32, 4)
    r = L.lookup_add(T, 6, 32, seed=4)
    assert r.method == 'bilinear' and r.n_and == 79
    assert sim_check(r, tmp_path, zs=(0, 1, (1 << 31), (1 << 32) - 2), m_trials=3)[0]


# ------------------------------------------------------------------ exact proofs and canaries
@pytest.mark.parametrize('c,w,seed', [(4, 4, 1), (5, 8, 3)])
def test_prove_narrow_small(c, w, seed, tmp_path):
    r = L.narrow_lookup(table(c, w, seed), c, w, seed=seed)
    assert r.method == 'bilinear'
    assert prove(r, tmp_path) == 'PROVED'


def test_prove_qroam1_pair_add_small(tmp_path):
    assert prove(L.qroam1_lookup(table(4, 3, 5), 4, 3), tmp_path, tag='q') == 'PROVED'
    assert prove(L.lookup_pair(table(4, 4, 6), 4, 4, mode='qroam1'), tmp_path, tag='pq') == 'PROVED'
    assert prove(L.lookup_pair(table(4, 4, 6), 4, 4, mode='narrow'), tmp_path, tag='pn') == 'PROVED'
    assert prove(L.lookup_add(table(4, 8, 2), 4, 8, seed=2), tmp_path, tag='a') == 'PROVED'
    assert prove(L.lookup_add(table(5, 8, 3), 5, 8, seed=3), tmp_path, tag='a5') == 'PROVED'


def test_canaries_rejected(tmp_path):
    r = L.lookup_pair(table(4, 4, 6), 4, 4, mode='narrow')
    assert prove(r, tmp_path, tag='orig') == 'PROVED'
    bad = canaries(r.qasm)
    assert set(bad) == {'t2tdg', 'dropCZfix', 'dropCX', 'dropParity'}
    for name, txt in bad.items():
        st = prove(r, tmp_path, qasm=txt, tag=name, budget=120)
        if st != 'REFUTED':        # the prover may give up on a mutant; the independent simulator must catch it
            assert st != 'PROVED', name
            assert not sim_check(r, tmp_path, qasm=txt, tag=name + '_s', m_trials=6)[0], name


@pytest.mark.parametrize('w,seed,ands', [(21, 5, 33), (32, 1, 40)])
def test_prove_narrow_c6(w, seed, ands, tmp_path):
    r = L.narrow_lookup(table(6, w, seed), 6, w, seed=seed)
    assert r.n_and == ands
    assert prove(r, tmp_path, budget=600) == 'PROVED'


def test_canaries_rejected_c6_w32(tmp_path):
    r = L.narrow_lookup(table(6, 32, 1), 6, 32, seed=1)
    for name, txt in canaries(r.qasm).items():
        st = prove(r, tmp_path, qasm=txt, tag=name, budget=120)
        if st != 'REFUTED':
            assert st != 'PROVED', name
            assert not sim_check(r, tmp_path, qasm=txt, tag=name + '_s', m_trials=6)[0], name


@pytest.mark.parametrize('what', ['narrow_6_32', 'add_6_32'])
def test_kprove_c6(what, tmp_path):
    kp = pytest.importorskip('qoptct.verify.kprove')       # qoptct >= 0.5 (BDD path-sum prover)
    if what == 'narrow_6_32':
        r = L.narrow_lookup(table(6, 32, 1), 6, 32, seed=1)
        assert r.n_and == 40
    else:
        r = L.lookup_add(table(6, 32, 4), 6, 32, seed=4)
        assert r.n_and == 79
    qp, qcp = tmp_path / 'k.qasm', tmp_path / 'k.qc'
    r.write(str(qp), str(qcp))
    assert kp.kprove(str(qp), str(qcp), clean=True, budget_sec=600)['status'] == 'PROVED'


@pytest.mark.slow
def test_and_distribution_c6_w32():
    """k1 measurement (Random(10000+t), 40 tables): 40 ANDs in 38/40, 41 in 2/40; here: every table <= 41."""
    cnt = {}
    for t in range(40):
        p = L.plan_narrow(table(6, 32, 10000 + t), 6, 32, seed=t, splits=[3])
        cnt[p.n_and] = cnt.get(p.n_and, 0) + 1
    assert set(cnt) <= {40, 41} and cnt.get(40, 0) >= 35, cnt


# ------------------------------------------------------------------ purity and input validation
def test_pure_and_deterministic():
    T = table(5, 6, 8)
    T0 = list(T)
    a = L.narrow_lookup(T, 5, 6, seed=2)
    b = L.narrow_lookup(T, 5, 6, seed=2)
    assert T == T0 and a.qasm == b.qasm and a.table == tuple(T0)


@pytest.mark.parametrize('args', [([1, 2, 3], 2, 2), ([0, 1, 2, 4], 2, 2), ([0, -1, 0, 0], 2, 2), ([0, 0], 0, 1),
                                  ([0, 0], 1, 0)])
def test_bad_input(args):
    with pytest.raises(ValueError):
        L.narrow_lookup(*args)
