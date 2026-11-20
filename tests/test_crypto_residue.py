"""Tests of qoptct.crypto.residue: compression of an accumulator modulo q ("loop2" of Gidney 2025) and the width of
the accumulator.

Levels: (1) exact simulation of the gate list (andcirc.simulate) on ALL inputs for small analogues (clen = 5..7) and on
edge/random inputs at the size of the paper (clen = 21, 32 / 29 bits); (2) the emitted OpenQASM through
qoptct.verify.arith.check_sim; (3) kprove (all inputs, all outcomes) for small analogues when available;
(4) canaries (mutants must be caught); (5) AND counts pinned to the lab claims of agent k2 (148 / 138 / 106)."""
import os
import random
import re
import sys

import numpy as np
import pytest

from qoptct.crypto import residue as R

Q21 = 1574500            # p - 1 for a 21-bit prime (the lab's RSA-2048 example)
PRIMES21 = [1574500, 1185822, 2023528]       # three values of p - 1, as in the lab's runs (clen = 21)


# ------------------------------------------------------------------ API and accumulator width
def test_api_errors():
    with pytest.raises(ValueError):
        R.fold_accumulator(100, 5, 9)             # q > 2^clen
    with pytest.raises(ValueError):
        R.fold_accumulator(10, 5, 5)              # nbits must exceed clen
    with pytest.raises(ValueError):
        R.fold_accumulator(0, 5, 9)
    with pytest.raises(ValueError):
        R.fold_accumulator(10, 5, 9, w=0)
    with pytest.raises(ValueError):
        R.fold_accumulator_ands(10, 5, 9, None)
    with pytest.raises(ValueError):
        R.accumulator_bits(-1, 3)
    with pytest.raises(ValueError):
        R.fold_accumulator(10, 5, 9).reference(1 << 9)


def test_accumulator_bits_paper_case():
    assert R.accumulator_bound(214, 21) == 214 * (2 ** 21 - 1) == 448790314
    assert R.accumulator_bits(214, 21) == 29
    assert 2 ** 28 < 448790314 < 2 ** 29
    assert R.accumulator_bits(213, 21) == 29 and R.accumulator_bits(128, 21) == 28
    assert R.accumulator_bits(0, 21) == 0 and R.accumulator_bits(1, 21) == 21


def test_accumulator_bits_is_least_width_bruteforce():
    """for small parameters: every sum fits in b bits and some sum needs exactly b bits."""
    import itertools
    for value_bits in (1, 2, 3):
        for terms in range(1, 6):
            sums = {sum(t) for t in itertools.product(range(1 << value_bits), repeat=terms)}
            b = R.accumulator_bits(terms, value_bits)
            assert max(sums) < (1 << b) and max(sums) >= (1 << (b - 1))


# ------------------------------------------------------------------ AND counts (lab claims of agent k2)
def test_and_counts_pinned():
    assert R.fold_accumulator(Q21, 21, 32, w=3).n_and == 148
    assert R.fold_accumulator(Q21, 21, 32, w=4).n_and == 138
    assert R.fold_accumulator(Q21, 21, 29, w=4).n_and == 106
    assert R.fold_accumulator_ands(Q21, 21, 32, 3) == 148 and R.fold_accumulator_ands(Q21, 21, 29, 4) == 106
    F = R.fold_accumulator(Q21, 21, 29)                       # auto window: never worse than the pinned ones
    assert F.n_and <= 106 and F.w in range(2, 7)
    c = F.counts()
    assert c['AND'] == F.n_and and c['T'] == 4 * F.n_and


def test_layout_and_text_counts():
    F = R.fold_accumulator(Q21, 21, 32, w=3)
    L = F.layout
    assert L['r'] == list(range(21)) and L['H'] == list(range(21, 32))
    assert L['top'] + L['g'] == list(range(32, F.nref)) and len(L['g']) == 1
    tc = F.text_counts()
    assert tc['AND'] == 148 and tc['T'] == 4 * 148 and tc['qubits'] == F.circ.N


# ------------------------------------------------------------------ exact simulation, all inputs (small analogues)
@pytest.mark.parametrize('clen,nbits', [(5, 9), (6, 10), (7, 11)])
def test_all_inputs_small(clen, nbits):
    """ALL S < 2^nbits, for many q <= 2^clen and windows 1..4: end property (r < 2^clen, r = S mod q), registers equal
    the reference, ancillas clean, every measurement erase consistent."""
    rng = random.Random(clen)
    qs = list(range(1, (1 << clen) + 1)) if clen == 5 else sorted({1 << clen, (1 << clen) - 1, (1 << clen) - 3,
                                                                  (1 << (clen - 1)) + 1, *(rng.randrange(1, 1 << clen)
                                                                                           for _ in range(6))})
    allS = list(range(1 << nbits))
    for q in qs:
        for w in ((1, 2, 3, 4) if clen == 5 and q % 4 == 0 else (2, 3)):
            F = R.fold_accumulator(q, clen, nbits, w)
            r = F.check(allS)
            assert r['ok'], (clen, nbits, q, w, r)


@pytest.mark.parametrize('nbits,w', [(32, 3), (32, 4), (29, 4)])
@pytest.mark.parametrize('q', PRIMES21)
def test_full_size_random(q, nbits, w):
    F = R.fold_accumulator(q, 21, nbits, w)
    dom = 1 << nbits
    rng = random.Random(q + nbits)
    vals = R.edge_inputs(q, 21, nbits, k=64, seed=q) + [rng.randrange(dom) for _ in range(300)]
    r = F.check(vals)
    assert r['ok'], r


def test_reference_end_property_bruteforce():
    """the reference alone (no circuit): r < 2^clen and r = S mod q for all S, many parameters."""
    for clen in (4, 5, 6):
        for q in range(1, (1 << clen) + 1):
            for w in (1, 2, 3):
                F = R.fold_accumulator(q, clen, clen + 5, w)
                for S in range(1 << (clen + 5)):
                    r = F.reference(S)['r']
                    assert r < (1 << clen) and (r - S) % q == 0


def test_not_fully_reduced_is_possible():
    """r < 2^clen only (not r < q): witness input, so the docstring caveat is real."""
    q, clen = 22, 5
    F = R.fold_accumulator(q, clen, 9, 2)
    assert any(F.reference(S)['r'] >= q for S in range(1 << 9))


# ------------------------------------------------------------------ gate level (emitted QASM)
def _check_sim(F, vals, tmp_path, txt=None, m_trials=3, complete=False):
    from qoptct.verify import arith
    nb = F.nbits; nref = F.nref
    path = tmp_path / ('res_%d.qasm' % random.randrange(1 << 30))
    path.write_text(txt if txt is not None else F.qasm())
    X = np.zeros((len(vals), nref), np.uint8)
    for i, S in enumerate(vals):
        for j in range(nb):
            X[i, j] = (int(S) >> j) & 1

    def f(Xb):
        D = np.zeros((len(Xb), nref), np.uint8)
        for row in range(len(Xb)):
            S = 0
            for j in range(nb):
                S |= int(Xb[row, j]) << j
            ref = F.reference(S)
            for reg in ('r', 'H', 'top', 'g'):
                for j, q in enumerate(F.layout[reg]):
                    D[row, q] = (ref[reg] >> j) & 1
        return D
    return arith.check_sim(str(path), (nref, f, list(range(nb))), X=X, domain_complete=complete, m_trials=m_trials)


@pytest.mark.parametrize('q,clen,nbits,w', [(22, 5, 9, 2), (18, 5, 9, 2), (66, 7, 11, 2), (45, 6, 10, 3)])
def test_qasm_check_sim_small_all_inputs(q, clen, nbits, w, tmp_path):
    F = R.fold_accumulator(q, clen, nbits, w)
    r = _check_sim(F, list(range(1 << nbits)), tmp_path, complete=True)
    assert r['ok'] and r['proved_all_x'], (r.get('error'), r.get('worst'))


@pytest.mark.parametrize('nbits,w', [(32, 3), (29, 4)])
def test_qasm_check_sim_full_size(nbits, w, tmp_path):
    F = R.fold_accumulator(Q21, 21, nbits, w)
    vals = R.edge_inputs(Q21, 21, nbits, k=32, seed=3)
    r = _check_sim(F, vals, tmp_path, m_trials=3)
    assert r['ok'], (r.get('error'), r.get('worst'))


# ------------------------------------------------------------------ canaries
def _mutants(txt, rng):
    L = txt.split('\n'); out = []
    ifs = [i for i, l in enumerate(L) if l.startswith('if') and ('cz' in l or ' z ' in l)]
    ts = [i for i, l in enumerate(L) if l.startswith('t ')]
    cxs = [i for i, l in enumerate(L) if l.startswith('cx ')]
    i = rng.choice(ts); M = list(L); M[i] = 'tdg' + M[i][1:]; out.append(('t_to_tdg', '\n'.join(M)))
    i = rng.choice(cxs); M = list(L); del M[i]; out.append(('drop_cx', '\n'.join(M)))
    # the first erase correction (a CZ of the first lookup's monomial or adder carry); some later corrections are
    # redundant on the reachable domain (e.g. top0 AND top1 is identically 0), so a random choice would be flaky
    i = ifs[0]; M = list(L); M[i] = re.sub(r'(cz q\[\d+\], q\[\d+\]|z q\[\d+\]); ?', '', M[i], count=1)
    out.append(('drop_corr', '\n'.join(M)))
    return out


def test_canaries_gate_level(tmp_path):
    F = R.fold_accumulator(22, 5, 9, 2)
    rng = random.Random(11)
    allx = list(range(1 << 9))
    for kind, mtxt in _mutants(F.qasm(), rng):
        r = _check_sim(F, allx, tmp_path, txt=mtxt, m_trials=8, complete=True)
        assert not r['ok'], kind


def test_canaries_op_level():
    q, clen, nbits = 45, 6, 10
    F = R.fold_accumulator(q, clen, nbits, 3)
    allS = list(range(1 << nbits))
    assert F.check(allS)['ok']
    ops = F.circ.ops
    kinds = {'corr': [j for j, o in enumerate(ops) if o[0] == 'MZ' and o[2]],
             'cx': [j for j, o in enumerate(ops) if o[0] == 'cx']}
    rng = random.Random(2)
    for kind, idx in kinds.items():
        caught = 0
        for j in rng.sample(idx, min(6, len(idx))):
            M = list(ops)
            if kind == 'corr':
                M[j] = ('MZ', ops[j][1], ())
            else:
                del M[j]
            F.circ.ops = M
            try:
                caught += not F.check(allS)['ok']
            finally:
                F.circ.ops = ops
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


def spec_residue(F):
    """independent word-level spec of the construction: windowed table sums, a second fold of the carry bits, one
    conditional addition of 2^clen - q."""
    q, clen, nbits, w, e, J = F.q, F.clen, F.nbits, F.w, F.e, F.J
    hb = nbits - clen

    def spec(S):
        X = S.reg(range(nbits)); L = X.resize(clen)
        v = L.resize(clen + e + 2)
        for j in range(J):
            ww = min(w, hb - w * j)
            addr = S.reg(range(clen + w * j, clen + w * j + ww))
            tab = [(a << (clen + w * j)) % q for a in range(1 << ww)]
            v = (v + S.lookup(addr, tab, clen)).resize(clen + e + 2)
        top = v.shr(clen).resize(e); lo = v.resize(clen)
        u = (lo.resize(clen + 1) + S.lookup(top, [(a << clen) % q for a in range(1 << e)], clen)).resize(clen + 1)
        g = u.bits[clen]; lo2 = u.resize(clen)
        r = S.mux(g, (lo2 + ((1 << clen) - q)).resize(clen), lo2)
        out = {**S.out(range(clen), r), **S.out(range(nbits, nbits + e), top)}
        out[nbits + e] = g
        return out
    return spec


SLOW_MAX_MB = 1400


def _kprove_residue(kprove, q, clen, nbits, w, tmp_path, txt=None, max_mb=None):
    F = R.fold_accumulator(q, clen, nbits, w)
    path = tmp_path / 'res.qasm'; path.write_text(txt if txt is not None else F.qasm())
    kw = {} if max_mb is None else dict(max_mb=max_mb)
    r = kprove(str(path), spec=spec_residue(F), n=F.nref, clean_wires=list(range(nbits, F.nref)), budget_sec=300, **kw)
    return F, r


@pytest.mark.parametrize('q,clen,nbits,w', [(22, 5, 9, 2), (45, 6, 10, 3), (66, 7, 11, 2), (100, 7, 12, 3)])
def test_kprove_small(q, clen, nbits, w, tmp_path):
    kprove = _kprove()
    F, r = _kprove_residue(kprove, q, clen, nbits, w, tmp_path)
    assert r['status'] == 'PROVED', r.get('detail')


@pytest.mark.slow
@pytest.mark.parametrize('q,clen,nbits,w', [(1018, 10, 14, 3)])
def test_kprove_big(q, clen, nbits, w, tmp_path):
    kprove = _kprove()
    F, r = _kprove_residue(kprove, q, clen, nbits, w, tmp_path, max_mb=SLOW_MAX_MB)
    assert r['status'] == 'PROVED', (r.get('detail'), r.get('mem'))


def test_kprove_mutants(tmp_path):
    kprove = _kprove()
    q, clen, nbits, w = 22, 5, 9, 2
    F, r = _kprove_residue(kprove, q, clen, nbits, w, tmp_path)
    assert r['status'] == 'PROVED', r.get('detail')
    rng = random.Random(4)
    for kind, mtxt in _mutants(F.qasm(), rng):
        _, rm = _kprove_residue(kprove, q, clen, nbits, w, tmp_path, txt=mtxt)
        assert rm['status'] == 'REFUTED', (kind, rm['status'])
