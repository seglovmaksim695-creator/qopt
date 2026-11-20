"""Tests for qoptct.verify.algebra (word-level algebraic prover).

* correct circuits of qoptct.crypto (multiply / square with Karatsuba levels, fold reductions, modular multiplication)
  must be PROVED;
* soundness fuzz: single-gate mutants of small circuits are compared with an exhaustive classical check written here
  (all inputs; AND targets |0>, erasure conditions, claims, clean ancillas) -- a mutant that is wrong must never be
  PROVED;
* strict parser, CLI, kprove(engine='algebra').
"""
import json
import os
import random
import zlib

import pytest

from qoptct.crypto import fold, mulsq
from qoptct.verify import algebra as A


# ----------------------------------------------------------------------------- exhaustive truth (independent)
def exhaustive_ok(N, ops, spec):
    """True iff the op list satisfies the spec for ALL inputs in the gate model of the prover's docstring."""
    sp = A.normalize_spec(spec, N)
    inw = sp['in_wires']
    k = len(inw)
    assert k <= 16
    S = 1 << k
    M = (1 << S) - 1
    v = [0] * N
    for i, w in enumerate(inw):
        blk = ((1 << (1 << i)) - 1) << (1 << i)
        t = 0
        for off in range(0, S, 2 << i):
            t |= blk << off
        v[w] = t
    init = list(v)
    for o in ops:
        if o[0] == 'cx':
            v[o[2]] ^= v[o[1]]
        elif o[0] == 'x':
            v[o[1]] ^= M
        elif o[0] == 'and':
            if v[o[3]]:
                return False
            v[o[3]] = v[o[1]] & v[o[2]]
        else:
            d = v[o[1]]
            for c in o[2]:
                d ^= (v[c[1]] & v[c[2]]) if c[0] == 'cz' else v[c[1]]
            if d not in (0, M):
                return False
            v[o[1]] = 0
    outw = set(w for ws in sp['outputs'].values() for w in ws)
    for w in range(N):
        if w in sp['keep']:
            if v[w] != init[w]:
                return False
        elif w not in outw and w not in sp['garbage'] and v[w]:
            return False
    for s in range(S):
        env = {name: sum(((init[w] >> s) & 1) << i for i, w in enumerate(ws)) for name, ws in sp['inputs'].items()}
        for name, ws in sp['outputs'].items():
            env[name] = sum(((v[w] >> s) & 1) << i for i, w in enumerate(ws))
        for text, op, l, r in sp['claims']:
            a, b = A._ev(l, env), A._ev(r, env)
            if not {'==': a == b, '<=': a <= b, '<': a < b, '>=': a >= b, '>': a > b}[op]:
                return False
    return True


def mutants(N, ops, rng, count):
    out = []
    for _ in range(count):
        i = rng.randrange(len(ops))
        o = ops[i]
        kind = rng.randrange(4)
        new = list(ops)
        if o[0] in ('cx', 'x') and kind == 0:
            del new[i]
        elif o[0] == 'cx':
            w = rng.randrange(N)
            c, t = (w, o[2]) if kind == 1 else (o[1], w)
            if c == t:
                continue
            new[i] = ('cx', c, t)
        elif o[0] == 'x':
            new[i] = ('x', rng.randrange(N))
        elif o[0] == 'and':
            w = rng.randrange(N)
            a, b, t = [(w, o[2], o[3]), (o[1], w, o[3]), (o[1], o[2], w), (o[2], o[1], o[3])][kind]
            if len({a, b, t}) != 3:
                continue
            new[i] = ('and', a, b, t)
        else:
            corr = list(o[2])
            if corr and kind < 2:
                del corr[rng.randrange(len(corr))]
            elif kind == 2:
                w = rng.randrange(N)
                if w == o[1]:
                    continue
                corr.append(('z', w))
            else:
                del new[i]
                out.append(new)
                continue
            new[i] = ('erase', o[1], tuple(corr))
        out.append(new)
    return out


def _mul(n, L=0, kind='mul', **kw):
    c = (mulsq.multiply if kind == 'mul' else mulsq.square)(n, karatsuba_levels=L, **kw)
    N, ops = A.ops_of_mulsq(c)
    return N, ops, A.spec_mulsq(c)


def _fold(n, c):
    R = fold.reduce_special((1 << n) - c, n)
    N, ops = A.ops_of_andcirc(R.circ, R.order)
    return N, ops, A.spec_fold(R)


def _modmul(n, c, L=0, kind='mul'):
    R = fold.reduce_special((1 << n) - c, n)
    m = (mulsq.multiply if kind == 'mul' else mulsq.square)(n, karatsuba_levels=L)
    return A.compose_modmul(m, R)


# ----------------------------------------------------------------------------- correct circuits are proved
@pytest.mark.parametrize('kind', ['mul', 'sq'])
@pytest.mark.parametrize('n,L', [(1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (8, 0), (13, 0), (6, 1), (8, 1), (12, 1),
                                 (16, 1), (8, 2), (16, 2), (32, 2), (32, 3)])
def test_mulsq_proved(kind, n, L):
    N, ops, spec = _mul(n, L, kind)
    r = A.prove_ops(N, ops, spec)
    assert r['status'] == 'PROVED', r['detail']


def test_mul_legacy_and_text_roundtrip():
    c = mulsq.multiply(12, karatsuba_levels=1, legacy=True)
    r = A.aprove(c.to_qasm(), A.spec_mulsq(c))
    assert r['status'] == 'PROVED', r['detail']
    N, ops = A.ops_of_mulsq(c)
    assert A.parse_strict(A.ops_to_qasm(N, ops)) == (N, ops) == A.parse_strict(c.to_qasm())


@pytest.mark.parametrize('kind,n,L', [('mul', 4, 0), ('mul', 8, 1), ('sq', 8, 1), ('mul', 16, 2), ('sq', 16, 2)])
def test_clean_bennett_proved(kind, n, L):
    N, ops, spec = _mul(n, L, kind, clean=True)
    assert not spec['garbage']
    r = A.prove_ops(N, ops, spec)
    assert r['status'] == 'PROVED', r['detail']
    assert r['mirror']['uncompute'] > 0
    k = A.mirror_split(ops)[1]
    broken = ops[:k + 3] + ops[k + 4:]                 # the mirror is no longer exact: no shortcut, no proof
    assert A.prove_ops(N, broken, spec)['status'] != 'PROVED'


@pytest.mark.parametrize('n,c', [(8, 3), (12, 5), (16, 19), (24, 977), (20, 1), (48, 5073)])
def test_fold_proved(n, c):
    N, ops, spec = _fold(n, c)
    r = A.prove_ops(N, ops, spec)
    assert r['status'] == 'PROVED', r['detail']


@pytest.mark.parametrize('name', ['curve25519', 'secp256k1', 'p256'])
def test_fold_curves_proved(name):
    R = fold.reduce_special(fold.CURVES[name])
    r = A.aprove(R.qasm(), A.spec_fold(R))
    assert r['status'] == 'PROVED', r['detail']
    assert r['counts']['AND'] == R.n_and


@pytest.mark.parametrize('n,c,L,kind', [(8, 3, 0, 'mul'), (12, 5, 1, 'mul'), (16, 19, 0, 'sq'), (24, 977, 2, 'mul')])
def test_modmul_proved(n, c, L, kind):
    N, ops, spec = _modmul(n, c, L, kind)
    r = A.prove_ops(N, ops, spec)
    assert r['status'] == 'PROVED', r['detail']


def test_exhaustive_agrees_on_correct():
    for N, ops, spec in (_mul(3), _mul(4, 1), _mul(4, 0, 'sq'), _fold(5, 3), _modmul(5, 3)):
        assert exhaustive_ok(N, ops, spec)


# ----------------------------------------------------------------------------- soundness fuzz
@pytest.mark.parametrize('case,count', [(('mul', 3, 0), 120), (('mul', 4, 1), 120), (('sq', 4, 0), 120),
                                        (('sq', 6, 1), 80), (('fold', 5, 3), 120), (('fold', 6, 5), 80),
                                        (('modmul', 5, 3), 80), (('clean', 3, 0), 150), (('clean', 4, 1), 150)])
def test_mutants_never_proved_when_wrong(case, count):
    if case[0] in ('mul', 'sq'):
        N, ops, spec = _mul(case[1], case[2], case[0])
    elif case[0] == 'clean':
        N, ops, spec = _mul(case[1], case[2], 'mul', clean=True)
    elif case[0] == 'fold':
        N, ops, spec = _fold(case[1], case[2])
    else:
        N, ops, spec = _modmul(case[1], case[2])
    rng = random.Random(zlib.crc32(repr(case).encode()))
    nbad = 0
    for new in mutants(N, ops, rng, count):
        truth = exhaustive_ok(N, new, spec)
        r = A.prove_ops(N, new, spec, budget_sec=30)
        assert r['status'] in ('PROVED', 'REFUTED', 'UNKNOWN'), r
        if not truth:
            nbad += 1
            assert r['status'] != 'PROVED', (case, r['detail'])
        if r['status'] == 'REFUTED':
            assert not truth, (case, r['detail'])
    assert nbad > count // 3


def test_phase_canary():
    """dropping the Z correction of one erasure leaves the classical function intact: only the phase is wrong."""
    N, ops, spec = _mul(6)
    idx = [i for i, o in enumerate(ops) if o[0] == 'erase' and any(c[0] == 'z' for c in o[2])]
    assert idx
    for i in idx[::7]:
        o = ops[i]
        new = list(ops)
        new[i] = ('erase', o[1], tuple(c for c in o[2] if c[0] != 'z'))
        assert A.prove_ops(N, new, spec)['status'] in ('REFUTED', 'UNKNOWN')
        new[i] = ('erase', o[1], tuple(c for c in o[2] if c[0] != 'cz'))
        assert A.prove_ops(N, new, spec)['status'] in ('REFUTED', 'UNKNOWN')


def test_wrong_claims():
    N, ops, spec = _mul(8, 1)
    for claim in ('out == x*y + 1', 'out == x*x', 'out <= 65000', 'out + x == x*y'):
        s2 = dict(spec, claims=[claim])
        assert A.prove_ops(N, ops, s2)['status'] in ('REFUTED', 'UNKNOWN'), claim
    assert A.prove_ops(N, ops, dict(spec, claims=['out <= 65025', 'out >= 0']))['status'] == 'PROVED'
    assert A.prove_ops(N, ops, dict(spec, garbage=[]))['status'] != 'PROVED'       # Karatsuba garbage is not |0>


# ----------------------------------------------------------------------------- parser, CLI
def test_strict_parser_rejects():
    c = mulsq.multiply(3)
    txt = c.to_qasm()
    for bad in (txt.replace('cx q[', 'cy q[', 1), txt.replace('tdg q[', 't q[', 1), txt + 'z q[0];\n',
                txt + 'h q[0];\n', txt.replace('OPENQASM 3.0;', 'OPENQASM 2.0;'), txt + 'ccx q[0], q[1], q[2];\n',
                txt.replace('x q[', 'x q[9', 1)):
        with pytest.raises(A.ParseError):
            A.parse_strict(bad)
    r = A.aprove(txt + 'z q[0];\n', A.spec_mulsq(c))
    assert r['status'] == 'ERROR'
    assert A.aprove(txt, dict(inputs={'x': [0, 1, 2]}, claims=[]))['status'] == 'ERROR'
    assert A.aprove(txt, dict(inputs={'x': [0, 1, 99]}, outputs={'o': [3]}, claims=['o == x']))['status'] == 'ERROR'


def test_gadget_checks():
    assert A.gadget_check(A.AND13) < 1e-9 and A.gadget_check(A.AND10) < 1e-9
    broken = A.AND13[:6] + (('t', 'a'),) + A.AND13[7:]
    assert A.gadget_check(broken) > 1e-3


def test_cli_and_kprove_engine(tmp_path, capsys):
    from qoptct.verify.kprove import kprove, main as kprove_main
    c = mulsq.square(10, karatsuba_levels=1)
    q = tmp_path / 'c.qasm'
    s = tmp_path / 'spec.json'
    q.write_text(c.to_qasm())
    s.write_text(json.dumps(A.spec_mulsq(c)))
    assert A.main([str(q), str(s)]) == 0
    assert 'PROVED' in capsys.readouterr().out
    assert kprove(str(q), str(s), engine='algebra')['status'] == 'PROVED'
    assert kprove_main([str(q), str(s), '--engine', 'algebra']) == 0
    bad = tmp_path / 'bad.json'
    bad.write_text(json.dumps(dict(A.spec_mulsq(c), claims=['out == x*x + 1'])))
    assert A.main([str(q), str(bad), '--json']) == 1
    assert json.loads(capsys.readouterr().out.strip().splitlines()[-1])['status'] == 'REFUTED'
    assert A.main([str(q), str(tmp_path / 'missing.json')]) == 3


@pytest.mark.skipif(not os.environ.get('QOPTCT_SLOW'), reason='n = 256 proofs: set QOPTCT_SLOW=1 (about 1 min, 1 GB)')
@pytest.mark.parametrize('kind,L', [('mul', 0), ('mul', 3), ('sq', 0), ('sq', 3)])
def test_n256(kind, L):
    N, ops, spec = _mul(256, L, kind)
    r = A.prove_ops(N, ops, spec, budget_sec=1200)
    assert r['status'] == 'PROVED', r['detail']
