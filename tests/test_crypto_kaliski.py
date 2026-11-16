# SPDX-License-Identifier: Apache-2.0
"""Tests of qoptct.crypto.kaliski (Kaliski modular inversion, 12n^2 - 8n - 3 logical ANDs; legacy=True: 12n^2 + n).

Function checks are classical and exact (bit-sliced run over ALL x in [1, p), with a check of every AND target and
erase condition); the emitted OpenQASM text is checked independently (qoptct.verify: sparse simulation and the exact
prover) against a truth-table reference .qc built here from the definition, not from the circuit."""
import random
import re

import pytest

from qoptct.crypto import kaliski_inverse, kaliski_and_count, CheckError
from qoptct.crypto import kaliski as KM


def is_prime(q):
    return q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))


def primes(n, k=3):
    """smallest, a middle one and the largest odd prime with bit length n (k of them, distinct)."""
    ps = [q for q in range(2 ** (n - 1) + 1, 2 ** n, 2) if is_prime(q)]
    pick = [ps[0], ps[len(ps) // 2], ps[-1]][:k]
    return sorted(set(pick))


def pseudo_mersenne(n):
    return max(q for q in range(2 ** (n - 1) + 1, 2 ** n, 2) if is_prime(q))


def expect(p, n, x):
    return (-pow(x, -1, p) * pow(2, 2 * n, p)) % p


CASES = [(n, p) for n in range(3, 9) for p in primes(n)]


# ------------------------------------------------------------------ function, all inputs
@pytest.mark.parametrize('legacy', [False, True])
@pytest.mark.parametrize('n,p', CASES)
def test_forward_all_inputs_and_count(n, p, legacy):
    C = kaliski_inverse(p, legacy=legacy)
    assert C.n == n and C.x == list(range(n)) and C.legacy == legacy
    xs = list(range(1, p))
    got, v = C.evaluate(xs, return_wires=True)
    assert got == [expect(p, n, x) for x in xs]
    assert C.read(v, C.regs['u'], len(xs)) == [1] * len(xs)
    assert C.read(v, C.regs['v'], len(xs)) == [0] * len(xs)
    c = C.counts()
    if legacy:
        assert c['AND'] == kaliski_and_count(n, legacy=True) == 12 * n * n + n
        assert c['qubits'] == 11 * n
    else:
        assert c['AND'] == kaliski_and_count(n) == 12 * n * n - 8 * n - 3
        assert c['qubits'] == 11 * n - 2
        L = kaliski_inverse(p, legacy=True).counts()
        assert c['CX'] < L['CX'] and c['ERASE'] == L['ERASE'] - 9 * n
    assert c['AND'] < 26 * n * n + 2 * n          # Litinski 2023


def test_smallest_prime():
    C = kaliski_inverse(3)
    assert C.evaluate([1, 2]) == [expect(3, 2, 1), expect(3, 2, 2)]
    assert C.counts()['AND'] == kaliski_and_count(2) == 29


def _constant_ands(C, xs):
    """indices of the ANDs whose product a*b is the same on every input (found by a bit-sliced run)."""
    S = len(xs)
    M = (1 << S) - 1
    v = KM._simulate([], C.nqubits, KM._pack(xs, C.x), S)
    out = []
    for i, o in enumerate(C.ops):
        if o[0] == 'and' and (v[o[1]] & v[o[2]]) in (0, M):
            out.append(i)
        v = KM._simulate([o], C.nqubits, dict(enumerate(v)), S)
    return out


@pytest.mark.parametrize('n,p', [(n, p) for n in range(3, 8) for p in primes(n, 2)])
def test_removed_ands_were_constant(n, p):
    """of the 9n + 3 ANDs removed by default (P1-P5 in the module docstring), 9n + 1 are ANDs of identically constant
    product in the legacy circuit; the other 2 (w_{k-1}(u_0 ^ v_0) and the comparison carry of step 2n-2, P5) have a
    varying product but only serve to compute the constant w_{2n-2}. The remaining constant-product ANDs stay."""
    xs = list(range(1, p))
    old = _constant_ands(kaliski_inverse(p, legacy=True), xs)
    new = _constant_ands(kaliski_inverse(p), xs)
    assert len(old) - len(new) == 9 * n + 1


@pytest.mark.parametrize('legacy', [False, True])
@pytest.mark.parametrize('n,p', [(n, p) for n in range(3, 7) for p in primes(n, 2)])
def test_clean_all_inputs(n, p, legacy):
    C = kaliski_inverse(p, clean=True, legacy=legacy)
    assert C.x == list(range(n)) and C.out == list(range(n, 2 * n)) and C.garbage == []
    xs = list(range(1, p))
    got, v = C.evaluate(xs, return_wires=True)
    assert got == [expect(p, n, x) for x in xs]
    assert C.read(v, C.x, len(xs)) == xs
    assert all(v[w] == 0 for w in range(C.nqubits) if w not in C.x + C.out)
    F = kaliski_inverse(p, legacy=legacy).counts()
    c = C.counts()
    assert c['AND'] == F['AND'] + F['ERASE']
    assert c['qubits'] == 12 * n - (0 if legacy else 2)


@pytest.mark.parametrize('legacy', [False, True])
@pytest.mark.parametrize('n', range(4, 9))
def test_special_reduction(n, legacy):
    p = pseudo_mersenne(n)
    C = kaliski_inverse(p, reduction='special', legacy=legacy)
    xs = list(range(1, p))
    assert C.evaluate(xs) == [expect(p, n, x) for x in xs]


@pytest.mark.parametrize('legacy', [False, True])
def test_special_reduction_any_prime_clean(legacy):
    for p in (37, 41, 53):
        C = kaliski_inverse(p, reduction='special', clean=True, legacy=legacy)
        xs = list(range(1, p))
        got, v = C.evaluate(xs, return_wires=True)
        assert got == [expect(p, C.n, x) for x in xs]
        assert all(v[w] == 0 for w in range(C.nqubits) if w not in C.x + C.out)


def test_register_width_is_bit_length():
    """n must be p.bit_length(): the final long division has n+1 trial subtractions, enough only for p >= 2^(n-1)
    (with a wider register the result would not be fully reduced)."""
    assert kaliski_inverse(11, n=4).n == 4
    for n in (3, 5, 6):
        with pytest.raises(ValueError):
            kaliski_inverse(11, n=n)


@pytest.mark.parametrize('n', [12, 16, 24, 32])
def test_and_count_formula_count_only(n):
    """count only (no simulation) at larger n; one random sample of inputs at n = 12."""
    rng = random.Random(n)
    while True:
        p = rng.randrange(2 ** (n - 1) + 1, 2 ** n, 2)
        if is_prime(p):
            break
    C = kaliski_inverse(p)
    assert C.counts()['AND'] == 12 * n * n - 8 * n - 3
    assert C.counts()['qubits'] == 11 * n - 2
    L = kaliski_inverse(p, legacy=True).counts()
    assert L['AND'] == 12 * n * n + n and L['qubits'] == 11 * n
    if n == 12:
        xs = [1, 2, p - 1, (p + 1) // 2] + [rng.randrange(1, p) for _ in range(60)]
        assert C.evaluate(xs) == [expect(p, n, x) for x in xs]


def test_bad_arguments():
    for bad in (4, 1, 2, -7, 7.0):
        with pytest.raises(ValueError):
            kaliski_inverse(bad)
    with pytest.raises(ValueError):
        kaliski_inverse(11, reduction='montgomery')


# ------------------------------------------------------------------ mutants of the logical circuit are rejected
def _mutants(C, rng, k):
    ops = C.ops
    for _ in range(k):
        L = list(ops)
        kind = rng.randrange(4)
        if kind == 0:      # drop a CX
            i = rng.choice([i for i, o in enumerate(L) if o[0] == 'cx'])
            del L[i]
        elif kind == 1:    # AND -> CX from its first control
            i = rng.choice([i for i, o in enumerate(L) if o[0] == 'and'])
            L[i] = ('cx', L[i][1], L[i][3])
        elif kind == 2:    # drop an X
            i = rng.choice([i for i, o in enumerate(L) if o[0] == 'x'])
            del L[i]
        else:              # swap control and target of a CX
            i = rng.choice([i for i, o in enumerate(L) if o[0] == 'cx'])
            L[i] = ('cx', L[i][2], L[i][1])
        yield L


@pytest.mark.parametrize('n,p', [(5, 29), (6, 61)])
def test_logical_mutants_rejected(n, p):
    """random single-gate mutants of the clean circuit: the exact classical check (all x) must reject them; a few
    mutants may be neutral (a gate whose effect is cancelled or acts on a constant wire) and are tolerated."""
    C = kaliski_inverse(p, clean=True)
    xs = list(range(1, p))
    want = [expect(p, n, x) for x in xs]
    rng = random.Random(p)
    rejected = 0
    total = 16
    for L in _mutants(C, rng, total):
        try:
            v = KM._simulate(L, C.nqubits, KM._pack(xs, C.x), len(xs))
        except CheckError:
            rejected += 1
            continue
        got = KM._unpack(v, C.out, len(xs))
        dirty = any(v[w] for w in range(C.nqubits) if w not in C.x + C.out)
        kept = KM._unpack(v, C.x, len(xs)) == xs
        rejected += int(got != want or dirty or not kept)
    assert rejected >= total - 3, rejected


# ------------------------------------------------------------------ the OpenQASM text, independent checks
def reference_qc(p, n, path):
    """truth-table reference built from the definition: wires x0..x{n-1} (inputs) and o0..o{n-1} (clean, result);
    for every n-bit x: o ^= table(p, n)[x] under the n-controlled condition [x wires == x] (negative controls by X
    conjugation)."""
    X = ['x%d' % i for i in range(n)]
    O = ['o%d' % i for i in range(n)]
    L = ['.v ' + ' '.join(X + O), '.i ' + ' '.join(X), '', 'BEGIN']
    T = table(p, n)
    for x in range(1 << n):
        k = T[x]
        if not k:
            continue
        neg = [X[i] for i in range(n) if not (x >> i) & 1]
        L += ['tof %s' % w for w in neg]
        for j in range(n):
            if (k >> j) & 1:
                L.append('tof %s %s' % (' '.join(X), O[j]))
        L += ['tof %s' % w for w in neg]
    L.append('END')
    path.write_text('\n'.join(L) + '\n')
    return str(path)


def table(p, n):
    """the reference on all n-bit x: K(x mod p), and 0 for x == 0 mod p. On [1, p) this is the specification; the
    extension outside [1, p) is what the clean circuit does at these sizes (test_clean_out_of_domain)."""
    return {x: (expect(p, n, x % p) if x % p else 0) for x in range(1 << n)}


@pytest.mark.parametrize('legacy', [False, True])
@pytest.mark.parametrize('p', [5, 7, 11, 13])
def test_clean_out_of_domain(p, legacy):
    """outside [1, p) the clean circuit is still clean and gives K(x mod p) (0 for x = 0 mod p) at these sizes."""
    C = kaliski_inverse(p, clean=True, legacy=legacy)
    xs = list(range(1 << C.n))
    got, v = C.evaluate(xs, return_wires=True)
    T = table(p, C.n)
    assert got == [T[x] for x in xs]
    assert C.read(v, C.x, len(xs)) == xs
    assert all(v[w] == 0 for w in range(C.nqubits) if w not in C.x + C.out)


def _qasm_counts(txt):
    T = len(re.findall(r'^(t|tdg) q', txt, re.M))
    meas = len(re.findall(r'= measure', txt))
    return T, meas


def test_qasm_text_counts():
    for p in (7, 13):
        for clean in (False, True):
            C = kaliski_inverse(p, clean=clean)
            T, meas = _qasm_counts(C.to_qasm())
            c = C.counts()
            assert T == 4 * c['AND'] and meas == c['ERASE']
            assert ('qubit[%d] q;' % c['qubits']) in C.to_qasm()


@pytest.mark.parametrize('legacy', [False, True])
@pytest.mark.parametrize('p', [5, 7, 13, 31])
def test_qasm_sparse_and_prove(p, legacy, tmp_path):
    """the emitted OpenQASM text against the truth-table reference: sampled sparse simulation, then the exact
    path-sum prover (all inputs, all measurement outcomes, clean model)."""
    from qoptct import verify
    C = kaliski_inverse(p, clean=True, legacy=legacy)
    ref = reference_qc(p, C.n, tmp_path / 'ref.qc')
    q = C.to_qasm()
    s = verify.check_sparse(q, ref, clean=True, trials=2)
    assert s['ok'], s
    r = verify.prove(q, ref, clean=True, budget_sec=240)
    assert r['status'] == 'PROVED', r


@pytest.mark.slow
@pytest.mark.parametrize('legacy', [False, True])
@pytest.mark.parametrize('p', [61, 127])
def test_qasm_prove_slow(p, legacy, tmp_path):
    from qoptct import verify
    C = kaliski_inverse(p, clean=True, legacy=legacy)
    r = verify.prove(C.to_qasm(), reference_qc(p, C.n, tmp_path / 'ref.qc'), clean=True, budget_sec=600)
    assert r['status'] == 'PROVED', r


def _kprove(C, ref, tmp_path, txt=None, budget=300):
    from qoptct.verify.kprove import kprove
    qa = tmp_path / 'c.qasm'
    qa.write_text(C.to_qasm() if txt is None else txt)
    return kprove(str(qa), ref, clean=True, budget_sec=budget, engine='bdd', sample=False)


@pytest.mark.parametrize('p', [5, 7, 13, 31, pytest.param(61, marks=pytest.mark.slow),
                               pytest.param(127, marks=pytest.mark.slow)])
def test_qasm_kprove(p, tmp_path):
    """the second exact prover (BDD sum over paths, qoptct.verify.kprove) on the default clean circuit."""
    C = kaliski_inverse(p, clean=True)
    r = _kprove(C, reference_qc(p, C.n, tmp_path / 'ref.qc'), tmp_path)
    assert r['status'] == 'PROVED', r


_LINES = {'x': 1, 'cx': 1, 'and': 10, 'erase': 3}
_HEADER = 4


def _qasm_line(C, op_index):
    """first OpenQASM line of C.ops[op_index] in C.to_qasm() (x/cx: 1 line, AND: 10, erase: 3)."""
    return _HEADER + sum(_LINES[o[0]] for o in C.ops[:op_index])


def _live_erases(C, xs):
    """indices of erase ops whose CZ correction matters: a*b is not constant over the inputs at that point (a
    correction on a constant product is a global phase or nothing, and dropping it is a neutral mutant)."""
    out = []
    S = len(xs)
    v = KM._simulate([], C.nqubits, KM._pack(xs, C.x), S)
    for i, o in enumerate(C.ops):
        if o[0] == 'erase':
            ab = v[o[2]] & v[o[3]]
            if ab not in (0, (1 << S) - 1):
                out.append(i)
        v = KM._simulate([o], C.nqubits, dict(enumerate(v)), S)
    return out


@pytest.mark.parametrize('legacy', [False, True])
def test_qasm_mutants_refuted(legacy, tmp_path):
    """canaries on the OpenQASM text: a T -> T^dagger in an AND, a live CZ correction dropped, the copy CX of the
    result dropped, a CX inside the forward part dropped. The exact prover must refute each; the sampled sparse check
    (2 trials) must reject the three non-phase mutants (a dropped CZ is a phase error that 2 trials can miss)."""
    from qoptct import verify
    p = 13
    C = kaliski_inverse(p, clean=True, legacy=legacy)
    ref = reference_qc(p, C.n, tmp_path / 'ref.qc')
    L = C.to_qasm().splitlines()
    rng = random.Random(5)
    xs = list(range(1 << C.n))
    muts = []
    ia = rng.choice([i for i, o in enumerate(C.ops) if o[0] == 'and'])
    k = _qasm_line(C, ia) + 1
    assert L[k].startswith('t q')
    muts.append(L[:k] + [L[k].replace('t q', 'tdg q', 1)] + L[k + 1:])
    ie = rng.choice(_live_erases(C, xs))
    k = _qasm_line(C, ie) + 2
    assert L[k].startswith('if') and 'cz' in L[k]
    muts.append(L[:k] + [re.sub(r'cz q\[\d+\], q\[\d+\]; ', '', L[k], count=1)] + L[k + 1:])
    k = next(i for i, l in enumerate(L) if re.fullmatch(r'cx q\[\d+\], q\[%d\];' % C.out[0], l))
    muts.append(L[:k] + L[k + 1:])
    k = _qasm_line(C, next(i for i, o in enumerate(C.ops) if o[0] == 'cx' and i > 200))
    assert L[k].startswith('cx q')
    muts.append(L[:k] + L[k + 1:])
    for j, m in enumerate(muts):
        r = verify.prove('\n'.join(m) + '\n', ref, clean=True, budget_sec=120)
        assert r['status'] == 'REFUTED', r
        if j != 1:
            assert not verify.check_sparse('\n'.join(m) + '\n', ref, clean=True, trials=2)['ok']
        if j >= 1 and not legacy:      # kprove on T -> T^dagger takes minutes and ~0.9 GB: checked outside the tests
            assert _kprove(C, ref, tmp_path, txt='\n'.join(m) + '\n')['status'] == 'REFUTED'


@pytest.mark.parametrize('legacy', [False, True])
def test_reference_rejects_wrong_function(legacy, tmp_path):
    """the truth-table reference is not vacuous: the circuit for p = 13 does not prove against the table of p = 11."""
    from qoptct import verify
    C = kaliski_inverse(13, clean=True, legacy=legacy)
    ref = reference_qc(11, 4, tmp_path / 'ref.qc')
    r = verify.prove(C.to_qasm(), ref, clean=True, budget_sec=120)
    assert r['status'] == 'REFUTED', r
    assert _kprove(C, ref, tmp_path)['status'] == 'REFUTED'


@pytest.mark.slow
@pytest.mark.parametrize('n', [48, 64])
def test_and_count_formula_count_only_slow(n):
    rng = random.Random(n)
    while True:
        p = rng.randrange(2 ** (n - 1) + 1, 2 ** n, 2)
        if pow(2, p - 1, p) == 1 and all(pow(a, p - 1, p) == 1 for a in (3, 5, 7, 11, 13)):
            break
    assert kaliski_inverse(p).counts()['AND'] == 12 * n * n - 8 * n - 3
    assert kaliski_inverse(p, legacy=True).counts()['AND'] == 12 * n * n + n


def test_qasm_regression_hashes():
    """legacy=True: the emitted clean circuits are byte-identical to the lab circuits they were ported from (k4,
    circ/kal_k4_div_n{4,5,6}.qasm, already checked there by sparse simulation). Default: the hashes of the circuits
    proved in this version (P1-P5). A change of the generator shows here."""
    import hashlib
    want = {13: '17b6191c88329c6a', 31: '9cf5ec8ab0dc387f', 61: 'f5054b32241ae363'}
    for p, h in want.items():
        assert hashlib.sha256(kaliski_inverse(p, clean=True, legacy=True).to_qasm().encode()).hexdigest()[:16] == h
    want = {13: '23165abb8e6a8074', 31: 'e99bb8bcbc876500', 61: 'd4ddf16031cd4fe5'}
    for p, h in want.items():
        assert hashlib.sha256(kaliski_inverse(p, clean=True).to_qasm().encode()).hexdigest()[:16] == h
