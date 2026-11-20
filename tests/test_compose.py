"""qoptct.verify.compose: compositional proof of the Kaliski inversion (lifting, contracts on ROBDDs, chain)."""
import copy
import json
import random

import pytest

from qoptct.crypto import kaliski as K
from qoptct.crypto.kaliski import kaliski_inverse
from qoptct.verify.bdd import BDD
from qoptct.verify.compose import (verify_kaliski, pstd_test, check_gadget_lemmas, lift, LiftError, kaliski_hints)
from qoptct.verify.compose import bv
from qoptct.verify.compose.__main__ import main as cli


def _primes(lo, hi):
    return [p for p in range(max(lo, 3), hi, 2) if all(p % q for q in range(3, int(p ** 0.5) + 1, 2))]


def test_gadget_lemmas_exact():
    assert check_gadget_lemmas() is True


def test_bv_definitions_match_integers():
    rnd = random.Random(5)
    for w in (1, 2, 3, 7, 16, 64, 256):
        S = 64
        al = bv.IntAlg(S)
        xs = [rnd.randrange(1 << w) for _ in range(S)]
        ys = [rnd.choice([x, rnd.randrange(1 << w)]) for x in xs]
        kc = rnd.randrange(1 << (w + 1))

        def pack(vals):
            return [sum(((v >> i) & 1) << j for j, v in enumerate(vals)) for i in range(w)]

        def unpack(bits):
            return [sum(((b >> j) & 1) << i for i, b in enumerate(bits)) for j in range(S)]
        X, Y = pack(xs), pack(ys)
        s, c = bv.add(al, X, Y)
        assert unpack(s) == [(x + y) % (1 << w) for x, y in zip(xs, ys)]
        assert [(c >> j) & 1 for j in range(S)] == [(x + y) >> w for x, y in zip(xs, ys)]
        d, b = bv.sub(al, X, Y)
        assert unpack(d) == [(x - y) % (1 << w) for x, y in zip(xs, ys)]
        assert [(b >> j) & 1 for j in range(S)] == [int(x < y) for x, y in zip(xs, ys)]
        assert [(bv.le(al, X, Y) >> j) & 1 for j in range(S)] == [int(x <= y) for x, y in zip(xs, ys)]
        assert [(bv.lt_const(al, X, kc) >> j) & 1 for j in range(S)] == [int(x < kc) for x in xs]
        assert unpack(bv.const(al, kc % (1 << w), w)) == [kc % (1 << w)] * S
        assert [(bv.allzero(al, X) >> j) & 1 for j in range(S)] == [int(x == 0) for x in xs]


def test_bv_bdd_algebra_exhaustive():
    w = 4
    B = BDD()
    al = bv.BddAlg(B)
    xv = [B.new_var() for _ in range(w)]
    yv = [B.new_var() for _ in range(w)]
    X, Y = [e for _, e in xv], [e for _, e in yv]
    s, c = bv.add(al, X, Y)
    d, b = bv.sub(al, X, Y)
    le = bv.le(al, X, Y)
    lt9 = bv.lt_const(al, X, 9)
    for x in range(1 << w):
        for y in range(1 << w):
            asg = {xv[i][0]: (x >> i) & 1 for i in range(w)}
            asg.update({yv[i][0]: (y >> i) & 1 for i in range(w)})
            ev = lambda bits: sum(int(B.eval(e, asg)) << i for i, e in enumerate(bits))
            assert ev(s) == (x + y) % 16 and B.eval(c, asg) == (x + y >= 16)
            assert ev(d) == (x - y) % 16 and B.eval(b, asg) == (x < y)
            assert B.eval(le, asg) == (x <= y) and B.eval(lt9, asg) == (x < 9)


def test_lift_roundtrip_and_streamed_text(tmp_path):
    C = kaliski_inverse(31, clean=True)
    txt = C.to_qasm()
    C.write_qasm(tmp_path / 'c.qasm')
    assert (tmp_path / 'c.qasm').read_text() == txt
    L = lift(str(tmp_path / 'c.qasm'))
    assert L.nq == C.nqubits and L.n == len(C.ops)
    from qoptct.verify.compose.lift import write_qasm
    write_qasm(L, tmp_path / 'back.qasm', chunk=1000)
    assert (tmp_path / 'back.qasm').read_text() == txt
    for i, o in enumerate(C.ops):
        k, a, b, t, l = L.op(i)
        if o[0] == 'x':
            assert (k, a) == (0, o[1])
        elif o[0] == 'cx':
            assert (k, a, b) == (1, o[1], o[2])
        elif o[0] == 'and':
            assert (k, a, b, t) == (2, o[1], o[2], o[3])
        else:
            assert (k, t, a, b, l) == (3, o[1], o[2], o[3], o[4][0] if o[4] else -1)


def test_lift_is_strict():
    C = kaliski_inverse(13)
    L = C.to_qasm().split('\n')
    it = next(i for i, s in enumerate(L) if s.startswith('tdg '))
    im = next(i for i, s in enumerate(L) if 'measure' in s)
    ii = next(i for i, s in enumerate(L) if s.startswith('if '))
    muts = [L[:it] + L[it + 1:],                                  # a T gate removed from an AND gadget
            L[:it] + ['t' + L[it][3:]] + L[it + 1:],              # tdg -> t
            L[:im] + L[im + 1:],                                  # measurement removed
            L[:ii] + [L[ii].replace('cz ', 'cx ')] + L[ii + 1:],  # wrong correction
            L[:ii] + [L[ii].replace(' x q', ' z q')] + L[ii + 1:],
            L[:it] + ['s q[0];'] + L[it:],                        # foreign gate
            L[:40] + ['h q[1];'] + L[40:]]
    for m in muts:
        with pytest.raises(LiftError):
            lift('\n'.join(m))
        r = verify_kaliski('\n'.join(m), 13)
        assert r['status'] == 'UNKNOWN' and r['stage'] == 'lift'


@pytest.mark.parametrize('p', [5, 7, 13, 31, 61, 127, 251])
@pytest.mark.parametrize('reduction', ['generic', 'special'])
@pytest.mark.parametrize('clean', [False, True])
def test_kaliski_proved(p, reduction, clean):
    C = kaliski_inverse(p, reduction=reduction, clean=clean)
    r = verify_kaliski(C.to_qasm(), p, reduction, clean)
    assert r['status'] == 'PROVED', r
    assert r['pstd_test']['ok'] and r['segments'] > 8 * C.n
    assert r['out'] == list(C.out)
    assert set(r['assumptions']) <= {'M1-width', 'M2-P3', 'M3-vshort', 'M4-rfit', 'M4-rs', 'M4-sfit', 'M5-last-even',
                                     'M5-last-w', 'M5-w', 'M6-rtop'}


@pytest.mark.parametrize('p', [pytest.param(2 ** 31 - 1, marks=pytest.mark.slow),
                               pytest.param(2 ** 61 - 1, marks=pytest.mark.slow)])
def test_kaliski_proved_wider(p):
    for reduction in ('generic', 'special'):
        C = kaliski_inverse(p, reduction=reduction, clean=True)
        r = verify_kaliski(C.to_qasm(), p, reduction, True)
        assert r['status'] == 'PROVED', r


def test_pstd_assumptions_exhaustive():
    """the mathematical assumptions M1..M6 and the result of P_std: every odd prime p < 2^7, every x in [1, p)."""
    for p in _primes(5, 128):
        for reduction in ('generic', 'special'):
            t = pstd_test(p, reduction, xs=range(1, p))
            assert t['ok'] and t['inputs'] == p - 1, (p, reduction, t)


@pytest.mark.slow
def test_pstd_assumptions_exhaustive_slow():
    for p in _primes(128, 1024):
        t = pstd_test(p, 'generic', xs=range(1, p))
        assert t['ok'], (p, t)


def test_pstd_detects_a_false_assumption():
    """outside the domain (x = 0, or x >= p) the assumptions do not hold: the test must say so."""
    assert not pstd_test(13, xs=[0])['ok']
    assert not pstd_test(15, xs=[3, 5, 6])['ok']               # p not prime: gcd(x, p) > 1


def _truth(C, ops):
    p, nw = C.p, C.nqubits
    xs = list(range(1, p))
    S = len(xs)
    try:
        v = K._simulate(ops, nw, K._pack(xs, C.x), S)
    except Exception:
        return False
    if K._unpack(v, C.out, S) != [C.expected(x) for x in xs]:
        return False
    if C.clean:
        keep = set(C.x) | set(C.out)
        return K._unpack(v, C.x, S) == xs and all(v[w] == 0 for w in range(nw) if w not in keep)
    return K._unpack(v, C.regs['u'], S) == [1] * S and K._unpack(v, C.regs['v'], S) == [0] * S


def _mutate(rnd, ops, nw):
    ops = list(ops)
    i = rnd.randrange(len(ops))
    o = ops[i]
    kind = rnd.choice(['del', 'dup', 'wire', 'swapadj', 'insx', 'inscx', 'flipcx', 'droplin', 'kindswap'])
    if kind == 'del':
        del ops[i]
    elif kind == 'dup':
        ops.insert(i, o)
    elif kind == 'wire':
        o = list(o)
        j = rnd.choice([j for j in range(1, len(o)) if isinstance(o[j], int)])
        o[j] = rnd.randrange(nw)
        ops[i] = tuple(o)
    elif kind == 'swapadj' and i + 1 < len(ops):
        ops[i], ops[i + 1] = ops[i + 1], ops[i]
    elif kind == 'insx':
        ops.insert(i, ('x', rnd.randrange(nw)))
    elif kind == 'inscx':
        ops.insert(i, ('cx',) + tuple(rnd.sample(range(nw), 2)))
    elif kind == 'flipcx' and o[0] == 'cx':
        ops[i] = ('cx', o[2], o[1])
    elif kind == 'droplin' and o[0] == 'erase' and o[4]:
        ops[i] = o[:4] + ((),)
    elif kind == 'kindswap' and o[0] == 'and':
        ops[i] = ('erase', o[3], o[1], o[2], ())
    else:
        return None
    for o in ops:
        if (o[0] == 'cx' and o[1] == o[2]) or (o[0] == 'and' and len({o[1], o[2], o[3]}) != 3) or \
                (o[0] == 'erase' and len({o[1], o[2], o[3], *o[4]}) != 3 + len(o[4])):
            return None
    return ops


@pytest.mark.parametrize('reduction', ['generic', 'special'])
@pytest.mark.parametrize('clean', [False, True])
def test_mutants_never_proved_wrongly(reduction, clean):
    """random single-operation mutants at p = 13: the verdict is compared with the ground truth (exhaustive classical
    run with all gadget conditions): PROVED only for mutants that are really correct, REFUTED only for wrong ones."""
    p = 13
    rnd = random.Random(2026)
    C = kaliski_inverse(p, reduction=reduction, clean=clean)
    h = C.hints()
    seen = {'PROVED': 0, 'REFUTED': 0, 'UNKNOWN': 0}
    n = 0
    while n < 60:
        ops = _mutate(rnd, C.ops, C.nqubits)
        if ops is None:
            continue
        n += 1
        r = verify_kaliski(K._to_qasm(ops, C.nqubits), p, reduction, clean, hints=h, test_samples=0)
        tr = _truth(C, ops)
        seen[r['status']] += 1
        assert not (r['status'] == 'PROVED' and not tr), r
        assert not (r['status'] == 'REFUTED' and tr), r
    assert seen['REFUTED'] > 20


def test_bad_certificate_cannot_prove():
    p = 31
    C = kaliski_inverse(p)
    q = C.to_qasm()
    h0 = json.loads(json.dumps(C.hints()))
    assert verify_kaliski(q, p, hints=h0, test_samples=0)['status'] == 'PROVED'
    rnd = random.Random(3)
    for _ in range(40):
        h = copy.deepcopy(h0)
        m = h['marks'][rnd.randrange(len(h['marks']))]
        what = rnd.choice(['pos', 'wire', 'drop', 'label'])
        if what == 'pos':
            m[0] += rnd.choice([-3, -1, 1, 2, 7])
        elif what == 'wire':
            keys = [k for k, v in m[2].items() if isinstance(v, int) and k not in ('k', 'j')]
            if not keys:
                continue
            m[2][rnd.choice(keys)] = rnd.randrange(C.nqubits)
        elif what == 'drop':
            h['marks'].remove(m)
        else:
            m[1] = 'swap' if m[1] != 'swap' else 'sub'
        r = verify_kaliski(q, p, hints=h, test_samples=0)
        # the circuit is correct, so REFUTED is impossible; a damaged certificate may only lose the proof
        assert r['status'] in ('UNKNOWN', 'PROVED'), r
        if what in ('drop', 'label'):
            assert r['status'] == 'UNKNOWN'


def test_wrong_parameters_not_proved():
    C = kaliski_inverse(61, clean=True)
    q = C.to_qasm()
    assert verify_kaliski(q, 59, clean=True)['status'] != 'PROVED'            # another prime of the same width
    assert verify_kaliski(q, 61, clean=False)['status'] != 'PROVED'
    assert verify_kaliski(q, 61, 'special', clean=True)['status'] != 'PROVED'
    assert verify_kaliski(kaliski_inverse(61, legacy=True).to_qasm(), 61)['status'] != 'PROVED'   # not supported
    assert verify_kaliski('no such file.qasm', 61)['status'] == 'UNKNOWN'


def test_cli(tmp_path, capsys):
    assert cli(['gen', '--prime', '31', '--reduction', 'special', '--clean', '-o', str(tmp_path)]) == 0
    base = tmp_path / 'kaliski_p31_special_clean'
    assert cli(['verify', '--prime', '31', '--reduction', 'special', '--clean', str(base) + '.qasm',
                '--hints', str(base) + '.hints.json', '--quiet']) == 0
    assert 'PROVED' in capsys.readouterr().out
    txt = (tmp_path / 'kaliski_p31_special_clean.qasm').read_text().split('\n')
    i = next(j for j, s in enumerate(txt) if s.startswith('cx ') and j > 300)
    (tmp_path / 'bad.qasm').write_text('\n'.join(txt[:i] + txt[i + 1:]))
    assert cli(['verify', '--prime', '31', '--reduction', 'special', '--clean', str(tmp_path / 'bad.qasm'),
                '--hints', str(base) + '.hints.json', '--json']) == 1
    assert json.loads(capsys.readouterr().out)['status'] == 'REFUTED'
    assert cli(['pstd', '--prime', '31', '--all']) == 0
    assert kaliski_hints(31)['marks'][0][1] == 'init'


def test_lift_andcirc_dialect(tmp_path):
    """the 13-line AND gadget of qoptct.crypto.andcirc is lifted too (gadget lemma checked for both templates)."""
    from qoptct.crypto.andcirc import Circ
    from qoptct.verify.compose.lift import write_qasm
    from qoptct.verify.compose.walker import simulate_flat
    c = Circ()
    a, b = c.new(5), c.new(5)
    co = c.new()[0]
    c.add(b, a, cout=co)
    txt = c.qasm()
    L = lift(txt)
    assert L.style == 2 and sum(1 for i in range(L.n) if L.op(i)[0] == 2) == c.nand
    write_qasm(L, tmp_path / 'back.qasm')
    assert (tmp_path / 'back.qasm').read_text().split('\n', 1)[1] == txt.split('\n', 1)[1]
    xs = [(x, y) for x in range(32) for y in range(32)]
    init = {}
    for i in range(5):
        init[a[i]] = sum(((x >> i) & 1) << j for j, (x, y) in enumerate(xs))
        init[b[i]] = sum(((y >> i) & 1) << j for j, (x, y) in enumerate(xs))
    v, bad = simulate_flat(L, init, len(xs))
    assert bad is None
    for j, (x, y) in enumerate(xs):
        got = sum(((v[w] >> j) & 1) << i for i, w in enumerate(b + [co]))
        assert got == x + y


def _z3():
    from qoptct.verify.compose import lemmas_z3
    if not lemmas_z3.z3_available():
        pytest.skip('z3 binary not installed')
    return lemmas_z3


def test_lemmas_z3_quick():
    """a part of the Z3 proof of M1..M6: the arithmetic lemmas, INIT, FINAL, the invariant on real runs."""
    L = _z3()
    for only in ('LEMMA', 'DIV', 'FOLD', 'QUOT', 'INIT', 'FINAL', 'CONCRETE: Inv holds along all runs for p = 13', 'M5: k = 2n-1'):
        r = L.run(timeout=60, only=only)
        assert r['ok'] and r['total'] >= 1, (only, r)
    # the trace used by the concrete queries is Kaliski's algorithm: it ends with u = 1, v = 0, r = K
    for p in (13, 31, 251):
        n = p.bit_length()
        for x in range(1, p):
            k, wp, A, B, R, S, a, b, cR, cS = L.kaliski_trace(p, x)[-1]
            assert (k, wp, A, B) == (2 * n, 1, 0, 1) and S % p == (-pow(x, -1, p) * pow(2, 2 * n, p)) % p


@pytest.mark.slow
def test_lemmas_z3_all():
    """M1..M6 for every n, every odd prime p, every x: all implications unsat, no negative control unsat."""
    r = _z3().run(timeout=120)
    assert r['ok'] and r['total'] >= 290 and not r['bad_controls'], r
