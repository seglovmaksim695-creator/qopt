# SPDX-License-Identifier: Apache-2.0
"""Tests of qoptct.crypto.mulsq (multiplication / squaring on conditionally inverted additions + Karatsuba).

Exhaustive classical checks (every AND / erase precondition, inputs restored, ancillas |0>) on all inputs up to
6 bits, AND counts against the closed formulas and the published n = 256 numbers (counting only, no simulation),
exact proofs (qoptct.verify.prove, clean model) against an independent schoolbook Toffoli reference, mutants."""
import itertools
import random
import re

import pytest

from qoptct.crypto import mulsq as M
from qoptct.crypto import multiply, square, CheckError


def run_all(c):
    m = c.n
    if c.kind == 'mul':
        P = list(itertools.product(range(1 << m), repeat=2))
        got = c.simulate([a for a, _ in P], [b for _, b in P])
        return sum(g != a * b for g, (a, b) in zip(got, P))
    xs = list(range(1 << m))
    return sum(g != a * a for g, a in zip(c.simulate(xs), xs))


# ------------------------------------------------------------------------------------- AND counts

@pytest.mark.parametrize('n', range(1, 17))
def test_base_formulas(n):
    c = multiply(n)
    k = c.counts()
    assert k['AND'] == n * n + 4 * n
    assert k['ERASE'] == k['AND'] and k['garbage'] == 0 and c.garbage == []
    s = square(n)
    k = s.counts()
    assert 2 * k['AND'] == n * n + 7 * n - 6
    assert k['ERASE'] == k['AND'] and s.garbage == []


# small widths, base = 1 (recursion down to 1 bit): (m, L) -> AND
KARA_MUL = {(2, 1): 32, (3, 1): 56, (4, 1): 64, (5, 1): 92, (6, 1): 102,
            (2, 2): 52, (3, 2): 111, (4, 2): 139, (5, 2): 179, (6, 2): 204}
KARA_SQ = {(2, 1): 17, (3, 1): 35, (4, 1): 41, (5, 1): 61, (6, 1): 68,
           (2, 2): 28, (3, 2): 69, (4, 2): 86, (5, 2): 117, (6, 2): 136}


@pytest.mark.parametrize('m,L', sorted(KARA_MUL))
def test_karatsuba_table_small(m, L):
    assert multiply(m, L, base=1).counts()['AND'] == KARA_MUL[m, L]
    assert square(m, L, base=1).counts()['AND'] == KARA_SQ[m, L]
    assert multiply(m, L, base=1, clean=True).counts()['AND'] == 2 * KARA_MUL[m, L]


# n = 256, base = 4: K_CLAIMS (k4) and the x5 audit (qubits of multiply L2 / L3)
N256 = [('mul', 0, 66560, 1537, 0), ('mul', 2, 42727, 2707, 1300), ('mul', 3, 37880, 3911, 2504),
        ('sq', 0, 33661, 1279, 0), ('sq', 2, 23333, 2126, 975), ('sq', 3, 22332, 3029, 1878)]


@pytest.mark.parametrize('kind,L,AND,qubits,garbage', N256)
def test_counts_n256(kind, L, AND, qubits, garbage):
    c = (multiply if kind == 'mul' else square)(256, karatsuba_levels=L)
    k = c.counts()
    assert (k['AND'], k['qubits'], k['garbage']) == (AND, qubits, garbage)
    assert k['T'] == 4 * AND


def test_karatsuba_saves_at_256():
    assert square(256, 3).counts()['AND'] < 0.68 * (256 * 259 // 2 - 1)     # -33 % vs n(n+3)/2 - 1 (2609.05625)
    assert multiply(256, 3).counts()['AND'] < 0.58 * 256 ** 2


# ------------------------------------------------------------------------------- function (classical)

@pytest.mark.parametrize('kind', ['mul', 'sq'])
@pytest.mark.parametrize('m', range(1, 7))
@pytest.mark.parametrize('L', [0, 1, 2])
def test_all_inputs(kind, m, L):
    f = multiply if kind == 'mul' else square
    c = f(m, L, base=1)
    assert run_all(c) == 0
    if L == 0:
        assert c.garbage == []
    elif m >= 2:
        assert c.garbage
    cc = f(m, L, base=1, clean=True)
    assert cc.garbage == [] and run_all(cc) == 0


@pytest.mark.parametrize('kind,n,L', [('mul', 16, 2), ('mul', 64, 3), ('sq', 64, 3), ('mul', 128, 2),
                                      ('sq', 128, 4)])
def test_random_and_edge_inputs(kind, n, L):
    rng = random.Random(n * 10 + L)
    top = (1 << n) - 1
    xs = [0, 1, top, top, 0, 1 << (n - 1)] + [rng.getrandbits(n) for _ in range(40)]
    ys = [top, top, top, 0, 0, top] + [rng.getrandbits(n) for _ in range(40)]
    if kind == 'mul':
        got = multiply(n, L).simulate(xs, ys)
        assert got == [a * b for a, b in zip(xs, ys)]
    else:
        assert square(n, L).simulate(xs) == [a * a for a in xs]


def test_layout_and_qasm():
    n = 5
    c = multiply(n, 1, base=2)
    assert c.x == list(range(n)) and c.y == list(range(n, 2 * n)) and c.out == list(range(2 * n, 4 * n))
    assert c.garbage == list(range(4 * n, 4 * n + len(c.garbage)))
    assert sorted(c.x + c.y + c.out + c.garbage + c.ancillas) == list(range(c.num_qubits))
    s = square(n)
    assert s.x == list(range(n)) and s.out == list(range(n, 3 * n)) and s.y == []
    txt = s.to_qasm()
    k = s.counts()
    assert txt.startswith('OPENQASM 3.0;') and 'qubit[%d] q;' % s.num_qubits in txt
    assert len(re.findall(r'^(t|tdg) q', txt, re.M)) == 4 * k['AND']
    assert txt.count('measure') == k['ERASE']


def test_errors():
    with pytest.raises(ValueError):
        multiply(0)
    with pytest.raises(ValueError):
        square(4, karatsuba_levels=-1)
    with pytest.raises(ValueError):
        M.prove(multiply(4, 1, base=2))          # garbage: not provable as a clean function


# ------------------------------------------------------------------------- mutants (classical)

def test_classical_mutants_rejected():
    c = multiply(4)
    # drop the last CX that writes into the product register: wrong product on some input
    i = max(j for j, o in enumerate(c.ops) if o[0] == 'cx' and o[2] in c.out)
    bad = M.MulSqCircuit('mul', 4, 0, 4, False, c.ops[:i] + c.ops[i + 1:], c.num_qubits, c.x, c.y, c.out, [])
    assert _broken(bad)
    # AND -> CX on the first AND: precondition / value broken
    j = next(j for j, o in enumerate(c.ops) if o[0] == 'and')
    a, b, t = c.ops[j][1:]
    ops = list(c.ops)
    ops[j] = ('cx', a, t)
    bad = M.MulSqCircuit('mul', 4, 0, 4, False, ops, c.num_qubits, c.x, c.y, c.out, [])
    assert _broken(bad)
    # same for the square: drop the first X
    s = square(4)
    j = next(j for j, o in enumerate(s.ops) if o[0] == 'x')
    bad = M.MulSqCircuit('sq', 4, 0, 4, False, s.ops[:j] + s.ops[j + 1:], s.num_qubits, s.x, [], s.out, [])
    assert _broken(bad)


def _broken(c):
    try:
        return run_all(c) > 0
    except CheckError:
        return True


# ------------------------------------------------------------------------ independent reference

def _qc_eval(text, inputs):
    """tiny classical evaluator of a Toffoli-only .qc (independent of qoptct): inputs {name: bit}."""
    names = None
    v = {}
    for line in text.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == '.v':
            names = p[1:]
            v = {q: 0 for q in names}
            v.update(inputs)
        elif p[0] == 'tof':
            *cs, t = p[1:]
            if all(v[q] for q in cs):
                v[t] ^= 1
    return names, v


@pytest.mark.parametrize('kind,n', [('mul', 1), ('mul', 3), ('sq', 1), ('sq', 3), ('sq', 4)])
def test_reference_qc_is_the_product(kind, n):
    txt = M.reference_qc(kind, n)
    o0 = 2 * n if kind == 'mul' else n
    for x in range(1 << n):
        for y in (range(1 << n) if kind == 'mul' else [x]):
            inp = {'q%d' % i: x >> i & 1 for i in range(n)}
            if kind == 'mul':
                inp.update({'q%d' % (n + i): y >> i & 1 for i in range(n)})
            names, v = _qc_eval(txt, inp)
            out = sum(v['q%d' % (o0 + i)] << i for i in range(2 * n))
            assert out == x * y
            rest = [q for q in names if int(q[1:]) >= o0 + 2 * n]
            assert all(v[q] == 0 for q in rest)                       # its ancillas are clean again
            assert all(v[q] == b for q, b in inp.items())             # inputs restored


# ------------------------------------------------------------------------------ exact proofs

# qoptct.verify.prove (path sum): about 0.1-1 s each, multiply n=5 about 15 s; multiply n=6 is out of its reach
# (UNKNOWN after 1000 s) -> the BDD prover kprove below (lab k9 / qoptct >= 0.5): n=6..8 in seconds.
PROVE = [('mul', 4, 0, False), ('sq', 4, 0, False), ('sq', 5, 0, False), ('sq', 6, 0, False),
         ('mul', 3, 1, True), ('mul', 4, 1, True), ('sq', 4, 1, True),
         pytest.param('mul', 5, 0, False, marks=pytest.mark.slow)]


@pytest.mark.parametrize('kind,n,L,clean', PROVE)
def test_prove_against_schoolbook(kind, n, L, clean):
    c = (multiply if kind == 'mul' else square)(n, L, base=2, clean=clean)
    r = M.prove(c, budget_sec=600)
    assert r['status'] == 'PROVED', r


KPROVE = [('mul', 6, 0, False), ('sq', 6, 0, False), ('mul', 5, 1, True), ('sq', 6, 1, True),
          ('mul', 6, 2, True), ('sq', 6, 2, True)]


@pytest.mark.parametrize('kind,n,L,clean', KPROVE)
def test_kprove_against_schoolbook(kind, n, L, clean):
    pytest.importorskip('qoptct.verify.kprove')
    c = (multiply if kind == 'mul' else square)(n, L, base=2, clean=clean)
    r = M.prove(c, budget_sec=600, engine='kprove')
    assert r['status'] == 'PROVED', r


def _prove_text(c, txt):
    import os
    import tempfile
    from qoptct import verify
    d = tempfile.mkdtemp()
    qa, qc = os.path.join(d, 'm.qasm'), os.path.join(d, 'r.qc')
    with open(qa, 'w') as fh:
        fh.write(txt)
    with open(qc, 'w') as fh:
        fh.write(c.reference_qc())
    return verify.prove(qa, qc, clean=True, budget_sec=600)


def test_prove_mutants_refuted():
    c = multiply(4)
    txt = c.to_qasm(num_qubits=c.reference_width())
    i = txt.index('\ntdg q[')
    mut = txt[:i] + '\nt q[' + txt[i + len('\ntdg q['):]                  # T-dagger -> T in the first AND
    assert _prove_text(c, mut)['status'] == 'REFUTED'
    L = txt.split('\n')
    j = max(k for k, s in enumerate(L) if s.startswith('cx ') and int(s.split('q[')[2].split(']')[0]) in c.out)
    assert _prove_text(c, '\n'.join(L[:j] + L[j + 1:]))['status'] == 'REFUTED'   # last CX into the product
    assert _prove_text(c, txt)['status'] == 'PROVED'
