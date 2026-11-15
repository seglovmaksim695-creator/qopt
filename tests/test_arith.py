# SPDX-License-Identifier: Apache-2.0
"""Tests of qoptct.contrib.arith (recognition of reversible arithmetic + synthesis on logical ANDs).

Benchmarks are read from $QOPT_BENCH (see conftest.py); tests that need a missing file are skipped.
The hidden-weighted-bit tests build their own reference .qc from the definition (no licensed file needed)."""
import math
import os
import random
import re
from pathlib import Path

import pytest

from qoptct.contrib import arith
from qoptct.contrib.arith import ir, tt

BENCH = Path(os.environ["QOPT_BENCH"])


def bench(rel):
    p = BENCH / rel
    if not p.exists():
        pytest.skip('benchmark %s not available' % rel)
    return str(p)


def dense_ok(qasm, qc, clean, tmp_path):
    """independent dense check (qoptct.verify_qasm) for small circuits."""
    from qoptct import verify_qasm
    p = tmp_path / ('c%d.qasm' % random.randrange(1 << 30))
    p.write_text(qasm)
    return verify_qasm.check(str(p), qc, trials=3, clean=clean) > 1 - 1e-9


def sparse_ok(qasm, qc, clean):
    from qoptct import verify
    return verify.check_sparse(qasm, qc, clean=clean, trials=3)['ok']


def counts(r):
    return {(c['T'], c['CX'], c['extra']) for c in r}


# ------------------------------------------------------------------ gadgets
@pytest.mark.parametrize('pol', ['ab', 'a!b', '!ab', 'or'])
def test_and_gadget_exact(pol):
    f = {'ab': lambda a, b: a & b, 'a!b': lambda a, b: a & ~b & 1, '!ab': lambda a, b: ~a & b & 1,
         'or': lambda a, b: a | b}[pol]
    L, nb = ir.lower_ops([('and', (0, 1, 2, pol))])
    prog = ir.parse(L)
    amps = set()
    for x in range(4):
        st = ir.run_forced(prog, x, [])
        want = x | (f(x & 1, x >> 1) << 2)
        assert list(st) == [want]
        amps.add(complex(round(st[want].real, 9), round(st[want].imag, 9)))
    assert amps == {1 + 0j}
    L, nb = ir.lower_ops([('and', (0, 1, 2, pol)), ('erase', (2, tuple(ir.and_phase(0, 1, pol))))])
    prog = ir.parse(L)
    for m in (0, 1):
        outs = [ir.run_forced(prog, x, [m]) for x in range(4)]
        assert all(list(o) == [x] for x, o in enumerate(outs))
        assert len({complex(round(o[x].real, 9), round(o[x].imag, 9)) for x, o in enumerate(outs)}) == 1


def test_toffoli_lowering():
    L, nb = ir.lower_ops([('ccx', (0, 1, 2))])
    prog = ir.parse(L)
    for x in range(8):
        st = ir.run_forced(prog, x, [])
        want = x ^ (4 if x & 3 == 3 else 0)
        assert list(st) == [want] and abs(st[want] - 1) < 1e-9


# ------------------------------------------------------------------ recognition + synthesis on benchmarks
def test_gidney_adder_n4_candidates(tmp_path):
    qc = bench('crypto/qc/add_gidney_n4.qc')
    r = arith.candidates(qc, clean=True, budget_sec=60)
    assert r and all('arith/add' in c['note'] for c in r)
    Ts = {ir.count(c['qasm'])[0] for c in r}
    assert 16 in Ts                     # 4 logical ANDs (n ANDs with carry-out)
    for c in r:
        assert dense_ok(c['qasm'], qc, True, tmp_path)


def test_rc_adder_6_clean(tmp_path):
    qc = bench('std/rc_adder_6.qc')
    r = arith.build_all(qc, True, 60, external=False)
    assert {(24, 51, 5), (43, 64, 0), (39, 61, 1)} <= counts(r)   # Gidney and Gidney/Cuccaro hybrids
    for c in r:
        if c['N'] <= 19:
            assert dense_ok(c['qasm'], qc, True, tmp_path)


def test_adder_8_clean_and_strict():
    qc = bench('std/adder_8.qc')
    r = arith.build_all(qc, True, 60, external=False)
    assert {(32, 69, 1), (35, 70, 0)} <= counts(r)
    for c in r:
        assert sparse_ok(c['qasm'], qc, True)
    # with dirty ancillas the reference is not an adder: no strict candidate
    assert arith.build_all(qc, False, 60, external=False) == []


def test_qcla_com_7_comparator():
    qc = bench('std/qcla_com_7.qc')
    r = arith.build_all(qc, True, 60, external=False)
    assert (32, 62, 0) in counts(r)
    c = [c for c in r if c['T'] == 32][0]
    assert c['qasm'].startswith('OPENQASM 3;') and 'measure' in c['qasm']
    assert sparse_ok(c['qasm'], qc, True)
    assert arith.build_all(qc, False, 60, external=False) == []


def test_qcla_adder_10_and_mod_7():
    r = arith.build_all(bench('large/qcla_adder_10.qc'), True, 60, external=False)
    assert (40, 95, 0) in counts(r)
    r = arith.build_all(bench('std/qcla_mod_7.qc'), True, 60, external=False)
    assert any(c['T'] == 100 and c['extra'] == 0 for c in r)


def test_mod_adder_1024():
    qc = bench('large/mod_adder_1024.qc')
    r = arith.build_all(qc, True, 60, external=False)
    assert {(36, 79, 1), (39, 80, 0)} <= counts(r)


def test_not_arithmetic():
    for rel in ('std/tof_3.qc', 'std/gf2^4_mult.qc', 'std/qft_4.qc'):
        assert arith.candidates(bench(rel), clean=True, budget_sec=20) == []
    assert arith.candidates('/nonexistent/x.qc', clean=True) == []


# ------------------------------------------------------------------ canaries
def test_broken_reference_is_not_recognised(tmp_path):
    src = Path(bench('crypto/qc/add_gidney_n4.qc')).read_text().splitlines()
    idx = [i for i, l in enumerate(src) if l.split()[:1] == ['tof'] and len(l.split()) == 4]
    for k in (0, len(idx) // 2, len(idx) - 1):
        lines = list(src)
        del lines[idx[k]]
        p = tmp_path / ('broken%d.qc' % k)
        p.write_text('\n'.join(lines) + '\n')
        for c in arith.build_all(str(p), True, 30, external=False):
            assert dense_ok(c['qasm'], str(p), True, tmp_path)   # anything returned must match the broken reference


def _mutants(qasm):
    L = qasm.splitlines()
    out = []
    i = next(i for i, l in enumerate(L) if l.startswith('t '))
    out.append(L[:i] + ['tdg ' + L[i][2:]] + L[i + 1:])
    i = next(i for i, l in enumerate(L) if l.startswith('if (') and 'cz' in l)
    out.append(L[:i] + [re.sub(r'cz q\[\d+\], q\[\d+\];', '', L[i], count=1)] + L[i + 1:])
    i = next(i for i, l in enumerate(L) if l.startswith('cx '))
    out.append(L[:i] + L[i + 1:])
    return ['\n'.join(m) + '\n' for m in out]


def test_canaries_rejected(tmp_path):
    qc = bench('std/rc_adder_6.qc')
    T = tt.load(qc, True)
    r = [c for c in arith.build_all(qc, True, 60, external=False) if c['N'] <= 19]
    assert r
    for bad in _mutants(r[0]['qasm']):
        body = [l for l in bad.splitlines() if not l.startswith(('OPENQASM', 'include', 'qubit', 'bit'))]
        nb = int(re.search(r'bit\[(\d+)\]', bad).group(1))
        ok, _ = ir.check_text(body, nb, T)
        assert not ok
        assert not dense_ok(bad, qc, True, tmp_path)


# ------------------------------------------------------------------ hidden weighted bit (self-made reference)
def hwb_def(x, n):
    w = bin(x).count('1')
    y = 0
    for i in range(n):
        if x >> i & 1:
            y |= 1 << ((i + w) % n)
    return y


def write_hwb_qc(path, n):
    """naive reference: weight counter (Toffolis), controlled rotations (Fredkin), counter uncomputed."""
    nb = max(1, (n - 1).bit_length())
    assert nb <= 3
    xs = ['x%d' % i for i in range(n)]
    ws = ['w%d' % k for k in range(nb)]
    names = xs + ws + ['t']
    inc = []
    for x in xs:
        if nb >= 3:
            inc += [('t', x, 'w0'), ('w2', 't', 'w1'), ('t', x, 'w0')]
        if nb >= 2:
            inc += [('w1', x, 'w0')]
        inc += [('w0', x)]
    rot = []
    for k in range(nb):
        s = (1 << k) % n
        if s == 0:
            continue
        g = math.gcd(n, s)
        for r0 in range(g):
            cyc = [(r0 + j * s) % n for j in range(n // g)]
            for j in cyc[1:]:
                a, b = xs[cyc[0]], xs[j]
                rot += [(a, b), (b, ws[k], a), (a, b)]
    gates = inc + rot + inc[::-1]

    def tof(g):
        tgt, ctl = g[0], g[1:]
        return 'tof ' + ' '.join(list(ctl) + [tgt])
    txt = ['.v ' + ' '.join(names), '.i ' + ' '.join(xs), '', 'BEGIN'] + [tof(g) for g in gates] + ['END']
    Path(path).write_text('\n'.join(txt) + '\n')


@pytest.mark.parametrize('n', [5, 6])
def test_hwb_from_definition(n, tmp_path):
    qc = str(tmp_path / ('hwb%d_ref.qc' % n))
    write_hwb_qc(qc, n)
    T = tt.load(qc, True)
    for x in range(1 << n):          # the self-made reference implements the definition
        assert sum(T.at(w, x) << i for i, w in enumerate(T.vars)) == hwb_def(x, n)
    r = arith.build_all(qc, True, 60, external=False)
    assert r and all(c['family'] == 'hwb' for c in r)
    for c in r:
        prog = ir.parse(c['qasm'].splitlines())
        rng = random.Random(n)
        amp = None
        outc = [rng.getrandbits(1) for _ in range(c['qasm'].count('measure'))]
        for x in range(1 << n):       # definition check on every input, independent of the truth table code
            st = ir.run_forced(prog, x, outc)
            assert len(st) == 1
            (k, a), = st.items()
            assert k == hwb_def(x, n)
            amp = a if amp is None else amp
            assert abs(a - amp) < 1e-9
        if c['N'] <= 18:
            assert dense_ok(c['qasm'], qc, True, tmp_path)


def test_hwb8_bench_if_present():
    p = BENCH / 'extra' / 'hwb8.qc'
    if not p.exists():
        pytest.skip('hwb8.qc is not distributed')
    r = arith.build_all(str(p), True, 60, external=False)
    assert any(c['T'] == 101 for c in r)
