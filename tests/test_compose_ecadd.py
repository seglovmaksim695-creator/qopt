"""Tests of qoptct.verify.compose.ecadd / ecmath: the compositional proof of the exact elliptic-curve point addition
(qoptct.crypto.ecpoint.ec_point_add, default construction).

Fast part: small curves (n = 4..8): PROVED through the lifted ops and through the OpenQASM text; the word-level
program W_ec on all points; Theorem T by enumeration; the polynomial identities; mutants compared with the truth
(exhaustive classical run over all points): a wrong circuit is never PROVED; bad parameters; doctored certificates.
Slow part (pytest -m slow): n = 256 for secp256k1 and P-256 (about 2 minutes and 1.7 GB each), full-size canaries."""
import importlib.util
import json
import os
import random

import pytest

from qoptct.crypto import ecpoint as E
from qoptct.crypto._ecp.ec import small_curve
from qoptct.verify.compose import verify_ecadd, wec_test, lift_circuit, lift, verify_ecadd_window, wec_test_window
from qoptct.verify.compose import ecadd, ecmath
from qoptct.verify.compose.lift import Lifted
from qoptct.verify.compose.__main__ import main as compose_main


def _canaries():
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts', 'ecadd_canaries.py')
    if not os.path.exists(path):
        pytest.skip('scripts/ecadd_canaries.py not available')
    spec = importlib.util.spec_from_file_location('ecadd_canaries', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def curves(n, k=2):
    rng = random.Random(1000 + n)
    out = []
    while len(out) < k:
        cv = small_curve(n, rng, prime_order=True, min_points=6, pfilter=lambda q: (1 << n) - q < (1 << (n // 2 + 1)))
        if all((cv['p'], cv['a'], cv['b']) != (o['p'], o['a'], o['b']) for o in out):
            out.append(cv)
    return out


def build(cv, Q, L=1, **kw):
    C = E.ec_point_add(cv['p'], cv['a'], cv['b'], Q, karatsuba_levels=L, **kw)
    return C, json.loads(json.dumps(C.hints()))


# ------------------------------------------------------------------ mathematics (Theorem T)
def test_polynomial_identities():
    r = ecmath.check_identities()
    assert len(r) == 4 and all(r.values())
    # the checker is not vacuous: a wrong identity is not accepted
    a, x, y, xq, yq = (ecmath.Pol.var(v) for v in ecmath.Pol.VARS)
    assert not ((x + xq) ** 2 - x ** 2 - xq ** 2).is_zero()
    assert ((x + xq) ** 2 - x ** 2 - 2 * x * xq - xq ** 2).is_zero()
    assert (y ** 3 - y * (x + 1)).reduce_y(x + 1).is_zero() and not (y ** 2 - x).reduce_y(x + 1).is_zero()


def test_theorem_T_enumeration():
    r = ecmath.exhaustive_T([5, 7, 11], all_D=True)
    assert r['curves'] == 64 and r['cases'] > 4000 and not r['violations']
    r = ecmath.exhaustive_T([13, 17, 19], max_curves=40)
    assert r['cases'] > 500 and not r['violations']


def test_theorem_T_model_detects_violations():
    """the field-level model is not vacuous: an unsafe D or a Q of order 3 gives violations."""
    cv = curves(5)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][0]
    assert ecmath.theorem_T_violations(p, a, b, Q, Q)                       # D = Q: the first denominator is 0
    m2Q = ecadd.ec_neg(ecadd.ec_add(Q, Q, p, a), p)
    assert ecmath.theorem_T_violations(p, a, b, Q, m2Q)                     # D = -2Q: the second denominator is 0
    assert ecadd.check_parameters(p, a, b, Q, Q) and ecadd.check_parameters(p, a, b, Q, m2Q)


@pytest.mark.slow
def test_theorem_T_enumeration_slow():
    r = ecmath.exhaustive_T([5, 7, 11, 13], all_D=True)
    assert r['curves'] == 120 and not r['violations']
    r = ecmath.exhaustive_T([17, 19, 23, 29, 31])
    assert r['curves'] == 992 and r['cases'] == 24020 and not r['violations']


def test_two_torsion_check():
    rng = random.Random(3)
    for p in (5, 7, 11, 13, 31, 61, 127, 251):
        for _ in range(30):
            a, b = rng.randrange(p), rng.randrange(p)
            brute = any((x ** 3 + a * x + b) % p == 0 for x in range(p))
            assert ecadd._two_torsion(p, a, b) == brute
    for p in (2003, 4099, 65537, 2 ** 31 - 1):                              # the gcd path against the generator's
        for _ in range(12):
            a, b = rng.randrange(p), rng.randrange(p)
            if (4 * a ** 3 + 27 * b * b) % p:
                assert ecadd._two_torsion(p, a, b) == E.has_two_torsion(p, a, b)
    for p in (2003, 4099):
        for _ in range(6):
            a, b = rng.randrange(p), rng.randrange(p)
            assert ecadd._two_torsion(p, a, b) == any((x ** 3 + a * x + b) % p == 0 for x in range(p))


# ------------------------------------------------------------------ PROVED on small curves
@pytest.mark.parametrize('n,L', [(4, 1), (5, 0), (5, 1), (6, 1), (6, 3), (7, 1), (8, 2)])
def test_proved_small(n, L):
    for cv in curves(n):
        p, a, b = cv['p'], cv['a'], cv['b']
        for Q in cv['points'][:2]:
            C, h = build(cv, Q, L)
            r = verify_ecadd(None, p, a, b, Q, hints=h, lifted=lift_circuit(C.circ))
            assert r['status'] == 'PROVED', r
            assert r['mirrors'] == 5 and r['wec_test']['ok'] and r['wec_test']['circuit_agrees']
            assert r['wec_test']['inputs'] == len(cv['points']) + 1               # all points and O
            for k in ('D-input', 'G1-denominator-nonzero', 'G2-lambda-recomputed', 'G3-flagged-result',
                      'G4-generic-result', 'M1-width'):
                assert r['assumptions'].get(k), (k, r['assumptions'])
            assert r['assumptions']['G1-denominator-nonzero'] == 2
            assert all(x['status'] == 'PROVED' for x in r['algebra'])


def test_proved_through_the_text(tmp_path):
    cv = curves(5)[1]                                        # p = 29: fold reduction (algebra sub-goals)
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][4]
    C, h = build(cv, Q)
    path = str(tmp_path / 'c.qasm')
    C.write_qasm(path)
    assert open(path).read() == C.qasm()                     # the streamed writer writes the same text
    L1, L2 = lift(path), lift_circuit(C.circ)
    assert L1.nq == L2.nq and L1.ops == L2.ops
    r = verify_ecadd(path, p, a, b, Q, hints=h)
    assert r['status'] == 'PROVED' and len(r['algebra']) == 14
    r = verify_ecadd(C.qasm(), p, a, b, Q, karatsuba_levels=1)   # the certificate regenerated from the library
    assert r['status'] == 'PROVED'
    # a text mutant of a gadget is rejected by the lifter
    txt = C.qasm().replace('tdg q[', 't q[', 1)
    r = verify_ecadd(txt, p, a, b, Q, hints=h)
    assert r['status'] == 'UNKNOWN' and r['stage'] == 'lift'


def test_special_reduction_path_has_no_algebra_folds():
    cv = curves(4)[0]                                        # p = 11: no one-fold reduction
    Q = cv['points'][3]
    C, h = build(cv, Q)
    assert not any(m[1] == 'fold' for m in h['marks']) and any(m[1] == 'sp_q' for m in h['marks'])
    r = verify_ecadd(None, cv['p'], cv['a'], cv['b'], Q, hints=h, lifted=lift_circuit(C.circ))
    assert r['status'] == 'PROVED' and {x['label'] for x in r['algebra']} == {'mul', 'sq'}
    assert 'product bound: out == x*y with x, y <= p - 1' in r['lemma_facts']


def test_all_Q_of_a_curve():
    cv = curves(5)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    for Q in cv['points']:
        C, h = build(cv, Q, 0)
        r = verify_ecadd(None, p, a, b, Q, hints=h, lifted=lift_circuit(C.circ))
        assert r['status'] == 'PROVED', (Q, r)


def test_wec_program_all_points_and_violations():
    cv = curves(6)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][5]
    C, h = build(cv, Q)
    t = wec_test(p, a, b, Q, hints=h, lifted=lift_circuit(C.circ))
    assert t['ok'] and t['inputs'] == len(cv['points']) + 1 and t['circuit_agrees']
    # a point outside the curve violates the domain assumption (and usually more)
    bad = next((x, y) for x in range(1, p) for y in range(1, p)
               if (y * y - x ** 3 - a * x - b) % p and x == Q[0])
    t = wec_test(p, a, b, Q, hints=h, points=[None, Q, bad])
    assert not t['ok'] and any(v['what'] == 'D-input' for v in t['violations'])


# ------------------------------------------------------------------ never PROVED wrongly
def test_mutants_never_proved_wrongly():
    """single-op mutants in all segment kinds (and doctored certificates); verdict of the chain alone (no sampled
    simulation) against the exhaustive truth over all points of the curve."""
    can = _canaries()
    st = can.run_small(count=45, seed=11, ns=[4, 5], levels=1, progress=False)
    assert st['base_proved'] == st['base'] == 8
    assert st['mutants'] > 200 and st['wrong'] > 120
    assert st['wrong_proved'] == 0, st['false_proved']
    assert st['correct_proved'] >= 1                         # equivalent mutants are still proved (semantic check)


@pytest.mark.slow
def test_mutants_never_proved_wrongly_slow():
    can = _canaries()
    st = can.run_small(count=250, seed=12, ns=[4, 5, 6], levels=1, progress=False)
    assert st['mutants'] > 2000 and st['wrong_proved'] == 0, st['false_proved']
    st = can.run_small(count=120, seed=13, ns=[5, 6], levels=3, progress=False)
    assert st['wrong_proved'] == 0, st['false_proved']


def test_mirror_frame_condition():
    """the apply part must not change a wire of the compute part (hole 1 of the v0.6 audit, here by construction):
    an x inserted on a wire of the product register inside the apply part is refuted or unknown, never PROVED."""
    can = _canaries()
    cv = curves(5)[1]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][2]
    C, h = build(cv, Q)
    L0 = lift_circuit(C.circ)
    marks = h['marks']
    j = next(i for i, m in enumerate(marks) if m[1] == 'apply')
    gar = next(m for m in marks[:j] if m[1] == 'mul')[2]['out'][7]          # a wire of the product register
    at = marks[j][0] - 1                                                    # inside the apply part
    ops = list(L0.ops)
    ops[5 * at:5 * at] = [0, gar, 0, 0, -1]
    L = Lifted(L0.nq, __import__('array').array('i', ops), 0)
    h2 = can.shift_marks(h, at, 1)
    r = verify_ecadd(None, p, a, b, Q, hints=h2, lifted=L)
    assert r['status'] in ('REFUTED', 'UNKNOWN')
    ok, why = can.chain_verdict(p, a, b, Q, h2, L)
    assert not ok


def test_wrong_parameters_and_certificates():
    cv = curves(5)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q, Q1 = cv['points'][0], cv['points'][7]
    C, h = build(cv, Q)
    L = lift_circuit(C.circ)
    assert verify_ecadd(None, p, a, b, Q1, hints=h, lifted=L)['status'] != 'PROVED'                 # another Q
    assert verify_ecadd(None, p, a, (b + 1) % p, Q, hints=h, lifted=L)['status'] == 'ERROR'         # Q not on E
    assert verify_ecadd(None, 33, 1, 1, (0, 1), hints=h, lifted=L)['status'] == 'ERROR'             # composite p
    h2 = dict(h, D=list(Q))
    assert verify_ecadd(None, p, a, b, Q, hints=h2, lifted=L)['status'] == 'ERROR'                  # unsafe D
    other = next(P for P in cv['points'] if ecadd.check_parameters(p, a, b, Q, P) is None and list(P) != h['D'])
    assert verify_ecadd(None, p, a, b, Q, hints=dict(h, D=list(other)), lifted=L)['status'] != 'PROVED'
    assert verify_ecadd(None, p, a, b, Q, hints=dict(h, marks=h['marks'][:-1]), lifted=L)['status'] == 'UNKNOWN'
    assert verify_ecadd(None, p, a, b, Q, hints=dict(h, marks='x'), lifted=L)['status'] == 'UNKNOWN'
    # a curve with a point of order 2, and a Q of order 3
    assert 'y = 0' in ecadd.check_parameters(7, 0, 1, (0, 1), (2, 3))                               # x^3 + 1: root 3
    rng = random.Random(1)
    while True:
        q = rng.choice([13, 19, 31, 37, 43])
        aa, bb = rng.randrange(q), rng.randrange(q)
        if (4 * aa ** 3 + 27 * bb * bb) % q == 0 or ecadd._two_torsion(q, aa, bb):
            continue
        t3 = [P for P in ecadd.curve_points(q, aa, bb) if ecadd.ec_add(P, P, q, aa)[0] == P[0]]
        if t3:
            assert 'order of Q' in ecadd.check_parameters(q, aa, bb, t3[0], t3[0])
            break


def test_legacy_construction_is_not_proved():
    cv = curves(5)[0]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][0]
    C, h = build(cv, Q, legacy=True)
    r = verify_ecadd(None, p, a, b, Q, hints=h, lifted=lift_circuit(C.circ))
    assert r['status'] == 'UNKNOWN'
    assert C.check([None] + cv['points']) == 0


def test_cli(tmp_path, capsys):
    cv = curves(5)[1]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][1]
    curve = '%d,%d,%d' % (p, a, b)
    q = ['--qx', str(Q[0]), '--qy', str(Q[1])]
    assert compose_main(['ecgen', '--curve', curve] + q + ['--levels', '1', '-o', str(tmp_path)]) == 0
    base = str(tmp_path / ('ecadd_p%d_L1' % p))
    assert compose_main(['ecverify', '--curve', curve, base + '.qasm', '--hints', base + '.hints.json']) == 0
    out = capsys.readouterr().out
    assert 'PROVED' in out and 'G4-generic-result' in out
    assert compose_main(['ecverify', '--curve', curve] + q + [base + '.qasm', '--levels', '1', '--json']) == 0
    assert json.loads(capsys.readouterr().out.strip().splitlines()[-1])['status'] == 'PROVED'
    assert compose_main(['ecwec', '--curve', curve] + q + ['--levels', '1']) == 0
    txt = open(base + '.qasm').read().replace('cx q[', 'cx q[1', 1)
    open(base + '.bad.qasm', 'w').write(txt)
    assert compose_main(['ecverify', '--curve', curve, base + '.bad.qasm', '--hints', base + '.hints.json',
                         '--quiet']) in (1, 2)


# ------------------------------------------------------------------ window variant (P, k) -> (P + T[k], k)
def window_case(n, w, L=1, seed=0):
    rng = random.Random(1000 + n)
    cv = small_curve(n, rng, prime_order=True, min_points=4 * (1 << w) + 8,
                     pfilter=lambda q: (1 << n) - q < (1 << (n // 2 + 1)))
    r2 = random.Random(n * 100 + w + seed)
    T = r2.sample(cv['points'], 1 << w)
    C = E.ec_point_add_window(cv['p'], cv['a'], cv['b'], T, w, karatsuba_levels=L)
    return cv, T, C, json.loads(json.dumps(C.hints()))


@pytest.mark.parametrize('n,w', [(5, 2), (6, 2), (6, 3), (7, 3), (8, 4)])
def test_window_proved_small(n, w):
    cv, T, C, h = window_case(n, w)
    p, a, b = cv['p'], cv['a'], cv['b']
    r = verify_ecadd_window(None, p, a, b, T, w, hints=h, lifted=lift_circuit(C.circ))
    assert r['status'] == 'PROVED', r
    assert r['mirrors'] == 5 and r['wec_test']['ok'] and r['wec_test']['circuit_agrees']
    assert r['wec_test']['inputs'] == (len(cv['points']) + 1) << w                 # all pairs (P, k)
    assert r['assumptions']['G1-denominator-nonzero'] == 2 and r['assumptions']['D-input'] == 1


def test_window_text_cli_and_errors(tmp_path, capsys):
    cv, T, C, h = window_case(6, 2)
    p, a, b = cv['p'], cv['a'], cv['b']
    path = str(tmp_path / 'w.qasm')
    C.write_qasm(path)
    r = verify_ecadd_window(path, p, a, b, T, 2, hints=h)
    assert r['status'] == 'PROVED'
    L = lift_circuit(C.circ)
    # another table (two entries exchanged): the circuit does not compute it
    T2 = [T[1], T[0]] + T[2:]
    assert verify_ecadd_window(None, p, a, b, T2, 2, hints=h, lifted=L)['status'] != 'PROVED'
    assert verify_ecadd_window(None, p, a, b, T[:3], 2, hints=h, lifted=L)['status'] == 'ERROR'
    assert verify_ecadd_window(None, p, a, b, T, 2, hints=dict(h, D=list(T[0])), lifted=L)['status'] == 'ERROR'
    t = wec_test_window(p, a, b, T, 2, hints=h, lifted=L)
    assert t['ok'] and t['inputs'] == 4 * (len(cv['points']) + 1)
    # CLI: the table (K + k) * B
    B = cv['points'][3]
    curve = '%d,%d,%d' % (p, a, b)
    args = ['--curve', curve, '--qx', str(B[0]), '--qy', str(B[1]), '--window', '2', '--levels', '1']
    table = [ecadd.ec_mul(1 + k, B, p, a) for k in range(4)]
    try:
        E.ec_point_add_window(p, a, b, table, 2, karatsuba_levels=1)
    except ValueError:
        pytest.skip('no safe point for this table')
    assert compose_main(['ecgen'] + args + ['-o', str(tmp_path)]) == 0
    base = str(tmp_path / ('ecadd_p%d_L1_w2' % p))
    assert compose_main(['ecverify', '--curve', curve, '--window', '2', base + '.qasm', '--hints',
                         base + '.hints.json', '--quiet']) == 0
    assert 'PROVED' in capsys.readouterr().out
    assert compose_main(['ecwec'] + args) == 0


def test_window_mutants_never_proved_wrongly():
    can = _canaries()
    st = can.run_small(count=40, seed=21, ns=[5, 6], levels=1, progress=False, window=2)
    assert st['base'] >= 4 and st['base_proved'] == st['base']
    assert st['wrong'] > 100 and st['wrong_proved'] == 0, st['false_proved']


@pytest.mark.slow
def test_window_mutants_never_proved_wrongly_slow():
    can = _canaries()
    st = can.run_small(count=200, seed=22, ns=[5, 6], levels=1, progress=False, window=2)
    assert st['wrong'] > 800 and st['wrong_proved'] == 0, st['false_proved']
    st = can.run_small(count=120, seed=23, ns=[6], levels=1, progress=False, window=3)
    assert st['wrong_proved'] == 0, st['false_proved']


# ------------------------------------------------------------------ n = 256
NEW = {'secp256k1': (3023913, 6916), 'P-256': (3031380, 6724)}
NEWW = {'secp256k1': (3027017, 9225), 'P-256': (3034484, 9033)}          # window, w = 4


@pytest.mark.slow
@pytest.mark.parametrize('name', ['secp256k1', 'P-256'])
def test_n256_window_proved(name):
    """w = 4; about 2.5 minutes and 2.1 GB per curve."""
    import gc
    cv = E.CURVES[name]
    p, a, b = cv['p'], cv['a'], cv['b']
    T = [E.ec_mul(123456789 + k, (cv['gx'], cv['gy']), p, a) for k in range(16)]
    C = E.ec_point_add_window(p, a, b, T, 4)
    c = C.counts()
    assert (c['AND'], c['qubits']) == NEWW[name]
    h = json.loads(json.dumps(C.hints()))
    L = lift_circuit(C.circ)
    del C
    gc.collect()
    r = verify_ecadd_window(None, p, a, b, T, 4, hints=h, lifted=L)
    assert r['status'] == 'PROVED', {k: v for k, v in r.items() if k != 'algebra'}
    assert len(r['algebra']) == 14 and r['mirrors'] == 5 and r['wec_test']['ok']


@pytest.mark.slow
@pytest.mark.parametrize('name', ['secp256k1', 'P-256'])
def test_n256_proved(name):
    """about 2 minutes and 1.7 GB per curve (lifted directly from the ops; the proofs of record lift the 1 GB text:
    `qoptct compose ecgen / ecverify`)."""
    import gc
    cv = E.CURVES[name]
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = E.ec_mul(123456789, (cv['gx'], cv['gy']), p, a)
    C = E.ec_point_add(p, a, b, Q)
    c = C.counts()
    assert (c['AND'], c['qubits']) == NEW[name]
    h = json.loads(json.dumps(C.hints()))
    L = lift_circuit(C.circ)
    del C
    gc.collect()
    r = verify_ecadd(None, p, a, b, Q, hints=h, lifted=L)
    assert r['status'] == 'PROVED', {k: v for k, v in r.items() if k != 'algebra'}
    assert len(r['algebra']) == 14 and r['mirrors'] == 5 and r['wec_test']['ok']
    # two full-size canaries in the same process: a flipped cx in the wrapper and a changed AND operand in a product
    can = _canaries()
    marks = h['marks']
    rnd = random.Random(1)
    for lab, kind in (('w_post', 'cx_flip'), ('mul', 'and_ctrl')):
        j = next(i for i, m in enumerate(marks) if m[1] == lab)
        lo, hi = marks[j - 1][0], marks[j][0]
        i = next(i for i in range((lo + hi) // 2, hi) if L.ops[5 * i] == (1 if kind == 'cx_flip' else 2))
        saved = L.ops[5 * i:5 * i + 5]
        m = can.mutate(Lifted(L.nq, saved, 0), kind, 0, rnd)
        L.ops[5 * i:5 * i + 5] = m[0]
        r = verify_ecadd(None, p, a, b, Q, hints=h, lifted=L, test_samples=0)
        L.ops[5 * i:5 * i + 5] = saved
        assert r['status'] in ('REFUTED', 'UNKNOWN'), (lab, r['status'])
