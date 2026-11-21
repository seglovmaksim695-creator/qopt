"""Regression tests for the 'erasure with constant 1' hole (found by an independent audit, 20.11.2026).

An erasure whose condition (wire xor correction) is the CONSTANT 1 used to be accepted by the algebra prover
(qoptct.verify.algebra).  The mirror lemma of the point-addition chain and of `mirror_split` assume that every
erasure is exact (condition 0), so a circuit with such an erasure in the compute half got a false PROVED:
  * from the chain alone (no final run of the circuit: `test_samples=0`, `--no-test`);
  * from the public algebra prover `aprove` / `prove_ops` (Bennett form F, copy, F^-1).
Now: the algebra prover accepts an erasure only when the condition is proved to be 0; verify_ecadd without the final
run never answers PROVED (status UNKNOWN, chain_status shows the chain verdict).
"""
import array
import json
import random

import pytest

from qoptct.crypto import ecpoint as E
from qoptct.crypto._ecp.ec import small_curve
from qoptct.verify import algebra as A
from qoptct.verify.compose import ecadd, lift_circuit, verify_ecadd, verify_ecadd_window
from qoptct.verify.compose.lift import Lifted
from qoptct.verify.compose.walker import ProofFailure

KX, KCX, KAND, KER = 0, 1, 2, 3


# ---------------------------------------------------------------- the algebra prover (hole 4)
def _hole4():
    """Bennett form: F has an erasure with the constant 1 (t = x0*x1 xor 1), copy, exact mirror F^-1."""
    k = 12
    pattern = 0b101101110101
    x0, x1 = k, k + 1
    t = k + 2
    base = k + 3
    F = [('x', i) for i in range(k) if not (pattern >> i) & 1]
    F.append(('and', 0, 1, base))
    for i in range(2, k):
        F.append(('and', base + i - 2, i, base + i - 1))
    C = base + k - 2
    u = C + 1
    o = u + 1
    N = o + 1
    F += [('and', x0, x1, t), ('and', t, C, u), ('x', t), ('erase', t, (('cz', x0, x1),))]
    Finv = []
    for g in reversed(F):
        Finv += A._mirror(g)
    ops = F + [('cx', C, o)] + Finv
    spec = dict(inputs=dict(a=list(range(k)), x=[x0, x1]), outputs=dict(o=[o]), claims=['o <= 1'])
    return N, ops, spec


def _exhaustive_phase_problem(N, ops, spec):
    """exhaustive bit-sliced run: (index of an erasure with an input-dependent phase or None, constant-1 erasures)."""
    ni = sum(len(w) for w in spec['inputs'].values())
    S = 1 << ni
    M = (1 << S) - 1
    v = [0] * N
    for i in range(ni):
        m = 0
        blk = (1 << (1 << i)) - 1
        for s in range(1 << i, S, 2 << i):
            m |= blk << s
        v[i] = m
    dep, const1 = None, []
    for i, g in enumerate(ops):
        if g[0] == 'cx':
            v[g[2]] ^= v[g[1]]
        elif g[0] == 'x':
            v[g[1]] ^= M
        elif g[0] == 'and':
            assert not v[g[3]]
            v[g[3]] = v[g[1]] & v[g[2]]
        else:
            d = v[g[1]]
            for c in g[2]:
                d ^= (v[c[1]] & v[c[2]]) if c[0] == 'cz' else v[c[1]]
            if d == M:
                const1.append(i)
            elif d and dep is None:
                dep = i
            v[g[1]] = 0
    return dep, const1


def test_hole4_bennett_const1_is_not_proved(tmp_path):
    N, ops, spec = _hole4()
    dep, const1 = _exhaustive_phase_problem(N, ops, spec)
    assert const1 and dep is not None            # the circuit really is wrong: an input-dependent phase appears
    r = A.prove_ops(N, ops, spec)
    assert r['status'] != 'PROVED', r
    q, s = tmp_path / 'h.qasm', tmp_path / 'h.spec.json'
    q.write_text(A.ops_to_qasm(N, ops))
    s.write_text(json.dumps(spec))
    r = A.aprove(str(q), str(s))                 # the public entry point, all defaults
    assert r['status'] != 'PROVED', r


def test_lone_and_bennett_const1_erasure_is_not_proved():
    ops = [('and', 0, 1, 2), ('cx', 2, 3), ('x', 2), ('erase', 2, (('cz', 0, 1),))]
    spec = dict(inputs=dict(a=[0], b=[1]), outputs=dict(o=[3]), claims=['o == a*b'])
    assert A.prove_ops(4, ops, spec)['status'] != 'PROVED'
    ops = [('and', 0, 1, 2), ('cx', 2, 3), ('x', 2), ('erase', 2, (('cz', 0, 1),)),
           ('cx', 3, 4),
           ('and', 0, 1, 2), ('x', 2), ('cx', 2, 3), ('erase', 2, (('cz', 0, 1),))]
    spec = dict(inputs=dict(a=[0], b=[1]), outputs=dict(o=[4]), claims=['o == a*b'])
    assert A.prove_ops(5, ops, spec)['status'] != 'PROVED'


def test_exact_erasure_still_proved():
    """the fix does not make the prover useless: an exact AND / erase pair is proved."""
    ops = [('and', 0, 1, 2), ('cx', 2, 3), ('erase', 2, (('cz', 0, 1),))]
    spec = dict(inputs=dict(a=[0], b=[1]), outputs=dict(o=[3]), claims=['o == a*b'])
    assert A.prove_ops(4, ops, spec)['status'] == 'PROVED'


# ---------------------------------------------------------------- the point addition chain (the 'paired mutant')
def _curve(n):
    rng = random.Random(1000 + n)
    return small_curve(n, rng, prime_order=True, min_points=6, pfilter=lambda q: (1 << n) - q < (1 << (n // 2 + 1)))


def _inv_ops(comp):
    out = []
    for kd, a, b, t, l in reversed(comp):
        if kd in (KX, KCX):
            out.append([kd, a, b, t, l])
        elif kd == KAND:
            out.append([KER, a, b, t, -1])
        else:
            out.append([KAND, a, b, t, -1])
            if l >= 0:
                out.append([KCX, l, t, 0, -1])
    return out


def _flat(nq, ops):
    o = array.array('i')
    for x in ops:
        o.extend(x)
    return Lifted(nq, o, 0)


def _pair_mutants(L, h, chain_cls, run_chain, labels=('mul', 'sq', 'fold'), per_label=1):
    """'X on the erased wire just before an erasure of a compute half, mirror half rebuilt as the exact inverse'.
    Yields (Lifted, hints).  The mirror check passes by construction: only the proofs of the segments stand between
    the mutant and PROVED."""
    blocks = []
    orig = chain_cls._mirror

    def rec(self, snap, mid, pos, written):
        if self.prove:
            blocks.append((snap['pos'], mid, self.pos, pos))
        return orig(self, snap, mid, pos, written)
    chain_cls._mirror = rec
    try:
        run_chain(L, h)
    finally:
        chain_cls._mirror = orig
    ops = [list(L.ops[5 * i:5 * i + 5]) for i in range(L.n)]
    marks = h['marks']
    out = []
    seen = {}
    for (s, e, lo, pos) in blocks:
        lab = next((m[1] for m in marks if s < m[0] <= e and m[1] in labels), None)
        if lab is None:
            continue
        cand = [i for i in range(s, e) if ops[i][0] == KER]
        if not cand or seen.get(lab, 0) >= per_label:
            continue
        seen[lab] = seen.get(lab, 0) + 1
        i = cand[len(cand) // 2]
        t = ops[i][3]
        comp = ops[s:i] + [[KX, t, 0, 0, -1]] + ops[i:e]
        new = ops[:s] + comp + ops[e:lo] + _inv_ops(comp) + ops[pos:]
        d1, d2 = 1, len(_inv_ops(comp)) - (pos - lo)          # +1 inside the compute half; marks after the mirror also shift by its growth
        g = dict(h)
        g['marks'] = [[m[0] + (d1 if m[0] > i else 0) + (d2 if m[0] >= pos else 0), m[1], m[2]] for m in marks]
        out.append((lab, i, _flat(L.nq, new), g))
    return out


def _chain_only(p, a, b, Q, h, L):
    try:
        ecadd.ECChain(p, a, b, tuple(Q), tuple(h['D']), h, lifted=L, budget=120.0).run()
        return True
    except ProofFailure:
        return False


@pytest.mark.parametrize('n', [5, 6])
def test_paired_const1_mutants_not_proved(n):
    cv = _curve(n)
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][1]
    C = E.ec_point_add(p, a, b, Q, karatsuba_levels=1)
    h = json.loads(json.dumps(C.hints()))
    L = lift_circuit(C.circ)
    assert _chain_only(p, a, b, Q, h, L)
    base = verify_ecadd(None, p, a, b, Q, hints=h, lifted=L)
    assert base['status'] == 'PROVED'
    muts = _pair_mutants(L, h, ecadd.ECChain, lambda L_, h_: _chain_only(p, a, b, Q, h_, L_))
    assert len(muts) >= 2
    for lab, i, Lm, hm in muts:
        assert not _chain_only(p, a, b, Q, hm, Lm), (lab, i)                       # the chain alone
        r = verify_ecadd(None, p, a, b, Q, hints=hm, lifted=Lm, test_samples=0)    # == --no-test
        assert r['status'] != 'PROVED' and r.get('chain_status') != 'PROVED', (lab, i, r['status'])
        r = verify_ecadd(None, p, a, b, Q, hints=hm, lifted=Lm)                    # default: full check
        assert r['status'] != 'PROVED', (lab, i, r['status'])


def test_paired_const1_mutant_window_not_proved():
    n, w = 5, 2
    cv = small_curve(n, random.Random(7000 + n), prime_order=True, min_points=4 * (1 << w) + 8,
                     pfilter=lambda q: (1 << n) - q < (1 << (n // 2 + 1)))
    p, a, b = cv['p'], cv['a'], cv['b']
    table = random.Random(1).sample(cv['points'], 1 << w)
    C = E.ec_point_add_window(p, a, b, table, w, karatsuba_levels=1)
    h = json.loads(json.dumps(C.hints()))
    L = lift_circuit(C.circ)

    def chain(L_, h_):
        try:
            ecadd.ECWinChain(p, a, b, table, w, tuple(h_['D']), h_, lifted=L_, budget=120.0).run()
            return True
        except ProofFailure:
            return False
    assert chain(L, h)
    muts = _pair_mutants(L, h, ecadd.ECWinChain, chain, per_label=1)
    assert muts
    for lab, i, Lm, hm in muts:
        assert not chain(Lm, hm), (lab, i)
        r = verify_ecadd_window(None, p, a, b, table, w, hints=hm, lifted=Lm, test_samples=0)
        assert r['status'] != 'PROVED' and r.get('chain_status') != 'PROVED', (lab, i, r['status'])
        r = verify_ecadd_window(None, p, a, b, table, w, hints=hm, lifted=Lm)
        assert r['status'] != 'PROVED', (lab, i, r['status'])


def test_no_test_never_proved():
    """test_samples=0 (the chain without the final strict run of the circuit) is never reported as PROVED."""
    cv = _curve(5)
    p, a, b = cv['p'], cv['a'], cv['b']
    Q = cv['points'][1]
    C = E.ec_point_add(p, a, b, Q, karatsuba_levels=1)
    h = json.loads(json.dumps(C.hints()))
    L = lift_circuit(C.circ)
    r = verify_ecadd(None, p, a, b, Q, hints=h, lifted=L, test_samples=0)
    assert r['status'] == 'UNKNOWN' and r['chain_status'] == 'PROVED' and r['stage'] == 'chain_only', r
    assert verify_ecadd(None, p, a, b, Q, hints=h, lifted=L)['status'] == 'PROVED'
