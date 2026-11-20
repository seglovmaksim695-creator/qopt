# SPDX-License-Identifier: Apache-2.0
"""Tests of qoptct.crypto.ecpoint (exact / complete mixed elliptic-curve point addition over a prime field).

Function checks are classical and exact (bit-sliced run; every AND target must be 0, every erase condition is checked,
all ancillas must end at 0). Reference = the textbook affine addition of qoptct.crypto._ecp.ec. Phases: sparse quantum
simulation of the emitted OpenQASM on a superposition, plus mutants of the circuit that must be rejected."""
import random

import pytest

from qoptct.crypto import ecpoint as E
from qoptct.crypto import CheckError
from qoptct.crypto._ecp.ec import small_curve, is_prime


def curves(n, k=2, minpts=6):
    """k curves of ODD prime group order with p = 2^n - small (special prime: the reduction of qoptct.crypto.fold)."""
    rng = random.Random(1000 + n)
    out = []
    while len(out) < k:
        cv = small_curve(n, rng, prime_order=True, min_points=minpts,
                         pfilter=lambda q: (1 << n) - q < (1 << (n // 2 + 1)))
        if all(cv['p'] != o['p'] or cv['a'] != o['a'] or cv['b'] != o['b'] for o in out):
            out.append(cv)
    return out


def rand_curve(n, rng):
    """random curve over the largest prime below 2^n with odd group order (no 2-torsion)."""
    p = next(q for q in range((1 << n) - 1, 1 << (n - 1), -2) if is_prime(q))
    while True:
        a, b = rng.randrange(p), rng.randrange(p)
        if (4 * a ** 3 + 27 * b * b) % p and not E.has_two_torsion(p, a, b):
            return p, a, b


def special_points(C):
    return C.exceptional_points()


# ------------------------------------------------------------------ small curves: ALL P x ALL Q
@pytest.mark.parametrize('n,idx', [(n, i) for n in range(4, 9) for i in (0, 1)])
@pytest.mark.parametrize('mult', ['sd', 'school'])
def test_all_P_all_Q(n, idx, mult):
    cv = curves(n)[idx]
    p, a, b = cv['p'], cv['a'], cv['b']
    pts = [None] + cv['points']
    qs = cv['points']                          # ALL Q
    if mult == 'school':
        qs = qs[::4]
    ands = set()
    for Q in qs:
        C = E.ec_point_add(p, a, b, Q, karatsuba_levels=1, multiplier=mult)
        assert C.check(pts) == 0, (cv, Q)
        ands.add(C.counts()['AND'])
    assert len(ands) <= 4                      # the count depends on Q only through the constants (a few ANDs)


@pytest.mark.parametrize('n', [4, 5])
def test_all_P_all_Q_exhaustive_small(n):
    """n = 4, 5: every Q of every curve, both Karatsuba depths 0 and 2 (small bases are school/sd anyway)."""
    for cv in curves(n):
        p, a, b = cv['p'], cv['a'], cv['b']
        pts = [None] + cv['points']
        for L in (0, 2):
            for Q in cv['points']:
                C = E.ec_point_add(p, a, b, Q, karatsuba_levels=L)
                assert C.check(pts) == 0


@pytest.mark.parametrize('n', [16, 24, 32, 48, 64])
def test_random_and_exceptional_points(n):
    rng = random.Random(n)
    p, a, b = rand_curve(n, rng)
    Q = E.find_point(p, a, b, rng)
    C = E.ec_point_add(p, a, b, Q, karatsuba_levels=2)
    pts = C.exceptional_points() + [E.find_point(p, a, b, rng) for _ in range(40)]
    assert C.check(pts) == 0
    # also the images of the exceptional set under the safe point: D + Q and -(D + Q)
    DQ = E.ec_add(C.D, Q, p, a)
    assert C.check([DQ, E.ec_neg(DQ, p), E.ec_add(Q, E.ec_add(Q, Q, p, a), p, a)]) == 0


@pytest.mark.parametrize('n', [5, 6, 8])
def test_reduction_solinas_option(n):
    """the package's own Solinas search (reduction='solinas') gives the same function as the fold reduction."""
    cv = curves(n)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][1]
    C = E.ec_point_add(p, a, b, Q, karatsuba_levels=1, reduction='solinas')
    assert C.check([None] + cv['points']) == 0


def test_counts_structure():
    cv = curves(6)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][3]
    c3 = E.ec_point_add(p, a, b, Q, karatsuba_levels=3).counts()
    c0 = E.ec_point_add(p, a, b, Q, karatsuba_levels=0).counts()
    assert c3['T'] == 4 * c3['AND'] and c3['CCX'] == 0
    assert c3['CX'] > 3 * c3['AND'] and c3['qubits'] > 0 and c3['n'] == 6
    assert c0['AND'] > 0 and c0['qubits'] <= c3['qubits']


# ------------------------------------------------------------------ errors
def test_argument_checks():
    cv = curves(5)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][0]
    with pytest.raises(ValueError):
        E.ec_point_add(p, a, b, None)
    with pytest.raises(ValueError):
        E.ec_point_add(p, a, b, (Q[0], (Q[1] + 1) % p))          # not on the curve
    with pytest.raises(ValueError):
        E.ec_point_add(p, a, b, Q, karatsuba_levels=9)
    with pytest.raises(ValueError):
        E.ec_point_add(p, a, b, Q, multiplier='x')
    with pytest.raises(ValueError):
        E.ec_point_add(p, a, b, Q, n=p.bit_length() + 1)
    with pytest.raises(ValueError):
        E.ec_point_add(15, a, b, Q)                               # not prime
    C = E.ec_point_add(p, a, b, Q)
    with pytest.raises(ValueError):
        C.simulate([(1, 1)])                                      # not on the curve
    with pytest.raises(ValueError):
        C.simulate([None], ks=[0])


def test_even_order_curve_rejected():
    rng = random.Random(3)
    n = 7
    p = next(q for q in range((1 << n) - 1, 1 << (n - 1), -2) if is_prime(q))
    seen = [0, 0]
    for _ in range(40):
        a, b = rng.randrange(p), rng.randrange(p)
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        brute = any((x ** 3 + a * x + b) % p == 0 for x in range(p))
        assert E.has_two_torsion(p, a, b) == brute
        seen[brute] += 1
    assert seen[0] and seen[1]
    # a curve with a 2-torsion point is refused by the builder
    a, b = next((a, b) for a in range(p) for b in range(p)
                if (4 * a ** 3 + 27 * b * b) % p and E.has_two_torsion(p, a, b))
    with pytest.raises(ValueError, match='even order'):
        E.ec_point_add(p, a, b, (0, 0))


def test_order_greater_than_4_for_Q():
    """curves of prime order 5 or 7 over tiny p: a Q of order 5 is accepted and handled (D is searched)."""
    for p in (11, 13, 17, 19, 23):
        for a in range(p):
            for b in range(p):
                if (4 * a ** 3 + 27 * b * b) % p == 0 or E.has_two_torsion(p, a, b):
                    continue
                pts = [(x, y) for x in range(p) for y in range(p) if (y * y - x ** 3 - a * x - b) % p == 0]
                if len(pts) + 1 == 7:
                    C = E.ec_point_add(p, a, b, pts[0], karatsuba_levels=0)
                    assert C.check([None] + pts) == 0
                    return
    pytest.skip('no curve of order 7 found')


# ------------------------------------------------------------------ window variant
def make_table(cv, w, rng):
    p, a = cv['p'], cv['a']
    pts = cv['points']
    while True:
        T = rng.sample(pts, 1 << w)
        if all(E.ec_add(t, t, p, a) is not None for t in T):
            break
    return T


@pytest.mark.parametrize('n,w', [(5, 2), (5, 3), (6, 2), (6, 3), (6, 4), (7, 3), (8, 4)])
@pytest.mark.parametrize('mult', ['sd', 'school'])
def test_window_all_pairs(n, w, mult):
    cv = curves(n, k=1, minpts=4 * (1 << w) + 8)[0]
    rng = random.Random(n * 100 + w)
    p, a, b = cv['p'], cv['a'], cv['b']
    T = make_table(cv, w, rng)
    C = E.ec_point_add_window(p, a, b, T, w, karatsuba_levels=1, multiplier=mult)
    pairs = [(P, k) for P in [None] + cv['points'] for k in range(1 << w)]
    assert C.check([P for P, _ in pairs], [k for _, k in pairs]) == 0
    cnt = C.counts()
    assert cnt['AND'] > 0 and cnt['qubits'] > 0


def test_window_errors():
    cv = curves(6, k=1, minpts=40)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    T = make_table(cv, 2, random.Random(1))
    with pytest.raises(ValueError):
        E.ec_point_add_window(p, a, b, T[:3], 2)
    with pytest.raises(ValueError):
        E.ec_point_add_window(p, a, b, T, 2, D=T[0])                # D in the exceptional set
    C = E.ec_point_add_window(p, a, b, T, 2)
    with pytest.raises(ValueError):
        C.simulate([None])                                          # needs ks


def test_window_medium():
    rng = random.Random(77)
    n, w = 32, 3
    p, a, b = rand_curve(n, rng)
    T = []
    while len(T) < 1 << w:
        t = E.find_point(p, a, b, rng)
        if t not in T:
            T.append(t)
    C = E.ec_point_add_window(p, a, b, T, w, karatsuba_levels=2)
    pairs = []
    for k in range(1 << w):
        Q = T[k]
        Q2 = E.ec_add(Q, Q, p, a)
        for P in (None, Q, E.ec_neg(Q, p), Q2, E.ec_neg(Q2, p), C.D, E.ec_neg(C.D, p)):
            pairs.append((P, k))
    pairs += [(E.find_point(p, a, b, rng), rng.randrange(1 << w)) for _ in range(20)]
    assert C.check([P for P, _ in pairs], [k for _, k in pairs]) == 0


# ------------------------------------------------------------------ mutants (classical) and phases (sparse quantum)
def _mutants(C, kinds, count, seed):
    """circuits with ONE operation deleted (kind in kinds), as lists of replacement op arrays."""
    ops = C.circ.ops
    a = ops.a
    idx = [i for i in range(len(ops)) if a[5 * i] in kinds]
    r = random.Random(seed)
    return r.sample(idx, count)


def _run_without(C, i, pts):
    ops = C.circ.ops
    saved = ops.a[5 * i:5 * i + 5]
    del ops.a[5 * i:5 * i + 5]
    try:
        return C.check(pts)
    except CheckError:
        return -1
    finally:
        ops.a[5 * i:5 * i] = saved


def test_classical_mutants_rejected():
    cv = curves(4)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    pts = [None] + cv['points']
    Q = cv['points'][2]
    C = E.ec_point_add(p, a, b, Q, karatsuba_levels=1)
    assert C.check(pts) == 0
    # n = 4 is tiny: many CX / AND are redundant on the valid data (constant products), so only rates are asserted;
    # every deleted X gate must be noticed.
    rate = {}
    for kind, cnt in ((0, 40), (1, 80), (3, 40)):
        hit = [_run_without(C, i, pts) != 0 for i in _mutants(C, (kind,), cnt, kind + 1)]
        rate[kind] = sum(hit) / len(hit)
    assert rate[0] == 1.0, rate
    assert rate[1] >= 0.55 and rate[3] >= 0.4, rate
    assert C.check(pts) == 0                                             # the circuit is restored


def test_phase_sparse_superposition_and_canaries():
    import re
    cv = curves(4)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][0]
    C = E.ec_point_add(p, a, b, Q, karatsuba_levels=1)
    pts = [None] + cv['points']
    txt = C.qasm()
    assert C.phase_check(pts, K=6, trials=3, seed=1, qasm=txt) > 1 - 1e-6
    exc = C.exceptional_points()[:5]                    # superposition of {O, Q, -Q, 2Q, -2Q}
    assert C.phase_check(exc, K=5, trials=2, seed=2, qasm=txt) > 1 - 1e-6
    # canaries: drop the cz of one correction / T -> T^dagger
    L = txt.splitlines()
    r = random.Random(0)
    worst_ok = 0
    for kind in range(4):
        L2 = list(L)
        if kind < 2:
            idx = [i for i, l in enumerate(L2) if l.startswith('if') and 'cz' in l]
            i = r.choice(idx)
            L2[i] = re.sub(r'cz q\[\d+\], q\[\d+\]; ', '', L2[i])
        else:
            idx = [i for i, l in enumerate(L2) if l.startswith('t q')]
            i = r.choice(idx)
            L2[i] = L2[i].replace('t q', 'tdg q', 1)
        f = C.phase_check(pts, K=6, trials=2, seed=3, qasm='\n'.join(L2) + '\n')
        if f < 1 - 1e-3:
            worst_ok += 1
    assert worst_ok >= 3, worst_ok                       # the phase check rejects the mutants


def test_qasm_counts_match():
    import re
    cv = curves(4)[0]
    C = E.ec_point_add(cv['p'], cv['a'], cv['b'], cv['points'][1], karatsuba_levels=1)
    txt = C.qasm()
    T = len(re.findall(r'\b(t|tdg) q', txt))
    assert T == C.counts()['T']


# ------------------------------------------------------------------ n = 256 (slow: minutes, about 0.5 GB)
BIG = {'P-256': dict(L3=3041421, L2=3082199, L0=3293519), 'secp256k1': dict(L3=3033954, L2=3074732, L0=3286052)}


@pytest.mark.slow
@pytest.mark.parametrize('name', ['P-256', 'secp256k1'])
def test_n256_points_and_counts(name):
    cv = E.CURVES[name]
    p, a, b = cv['p'], cv['a'], cv['b']
    rng = random.Random(5)
    Q = E.find_point(p, a, b, rng)
    C = E.ec_point_add(p, a, b, Q, karatsuba_levels=3)
    assert C.counts()['AND'] == BIG[name]['L3']                  # the lab number (K_CLAIMS, rows k4)
    pts = C.exceptional_points() + [E.ec_neg(E.ec_add(C.D, Q, p, a), p)] + [E.find_point(p, a, b, rng) for _ in range(12)]
    assert len(pts) >= 16 and None in pts and Q in pts
    assert C.check(pts) == 0


@pytest.mark.slow
def test_n256_count_only_levels_p256():
    cv = E.CURVES['P-256']
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = E.find_point(p, a, b, random.Random(5))
    for L in (0, 2):
        C = E.ec_point_add(p, a, b, Q, karatsuba_levels=L)
        assert C.counts()['AND'] == BIG['P-256']['L%d' % L]
