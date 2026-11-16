"""Tests of qoptct.crypto.gf2m (v0.5): CRT multiplier, in-place squaring, Itoh-Tsujii inversion / division.

Fast tests: small m with all inputs (bitsliced simulation, every AND / erasure precondition enforced), the exact
degree-2 argument, the exact prover on m = 3, 4, canaries, and the AND counts of the NIST sizes by counting (no
simulation).  @slow: gate-level circuits and exact checks at m = 163, counts at m = 409, 571."""
import os
import random
import tempfile

import pytest

from qoptct.crypto import gf2m as G
from qoptct.contrib.gf2m import crt, itoh
from qoptct.contrib.gf2m.field import pmod, pmul, is_irreducible
from qoptct.contrib.gf2m.qcirc import GF2Circuit, PreconditionError

NIST_MULT_AND = {163: 908, 233: 1343, 283: 1669, 409: 2501, 571: 3569}


def mulmod(a, b, f):
    return pmod(pmul(a, b), f)


def inv(x, m, f):
    return itoh.Field(m, f).inv(x)


def write_tmp(text, suffix):
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


def schoolbook_qc(m, f):
    """independent reference: c ^= a*b mod f with Toffolis (wires a, b, c in this order)."""
    names = ["a%d" % i for i in range(m)] + ["b%d" % i for i in range(m)] + ["c%d" % i for i in range(m)]
    lines = [".v " + " ".join(names), ".i " + " ".join(names[:2 * m]), "BEGIN"]
    for i in range(m):
        for j in range(m):
            v = pmod(1 << (i + j), f)
            for k in range(m):
                if v >> k & 1:
                    lines.append("tof a%d b%d c%d" % (i, j, k))
    lines.append("END")
    return "\n".join(lines) + "\n"


def all_mult_inputs(m, clean):
    A, B, Cv = [], [], []
    for a in range(1 << m):
        for b in range(1 << m):
            for c in ([0] if clean else range(1 << m)):
                A.append(a)
                B.append(b)
                Cv.append(c)
    return A, B, Cv


def mutate(C, kind, k=0):
    """canary circuits: drop the k-th CX / AND, or swap the sides of the k-th CX."""
    ops = list(C.ops)
    idx = [i for i, o in enumerate(ops) if o[0] == ('a' if kind == 'drop_and' else 'c')]
    i = idx[k % len(idx)]
    if kind == 'drop_cx' or kind == 'drop_and':
        del ops[i]
    elif kind == 'swap_cx':
        ops[i] = ('c', ops[i][2], ops[i][1])
    D = GF2Circuit(C.nq, ops, dict(C.regs), dict(C.meta))
    return D


# ---------------------------------------------------------------------------------------------- data, plans

def test_moduli_and_plans():
    for m in range(2, 40):
        assert is_irreducible(G.default_modulus(m))
    assert G.default_modulus(163) == (1 << 163) | (1 << 7) | (1 << 6) | (1 << 3) | 1
    for m, r in NIST_MULT_AND.items():
        assert crt.plan(m)[0] == r
    with pytest.raises(ValueError):
        G.gf2m_multiply(4, f=0b10101)          # reducible


def test_ring_table_is_valid():
    for h in crt.ranks():
        crt.ring_dec(h)                         # raises on a corrupt entry


@pytest.mark.parametrize("m", [163, 233, 283])
def test_mult_count_nist(m):
    c = G.gf2m_mult_count(m)
    assert c["AND"] == NIST_MULT_AND[m] and c["T"] == 4 * NIST_MULT_AND[m] and c["qubits"] == 3 * m + 1


@pytest.mark.slow
@pytest.mark.parametrize("m", [409, 571])
def test_mult_count_nist_big(m):
    assert G.gf2m_mult_count(m)["AND"] == NIST_MULT_AND[m]


def test_bilinear_auto_picks_shorter():
    for m in range(2, 11):
        dec, info = G.gf2m_bilinear(m)
        assert crt.check_algorithm(dec, m, G.default_modulus(m))
        assert len(dec) <= len(G.gf2m_bilinear(m, algorithm="crt")[0])
    assert len(G.gf2m_bilinear(4)[0]) == 9 and len(G.gf2m_bilinear(8)[0]) == 24


# ---------------------------------------------------------------------------------------------- multiplier

@pytest.mark.parametrize("mode", ["strict", "clean"])
@pytest.mark.parametrize("alg", ["auto", "crt"])
def test_multiply_small_all_inputs(mode, alg):
    clean = mode == "clean"
    for m in range(2, 7 if clean else 6):
        f = G.default_modulus(m)
        C = G.gf2m_multiply(m, mode=mode, algorithm=alg)
        A, B, Cv = all_mult_inputs(m, clean)
        out = C.simulate(dict(a=A, b=B, c=Cv), len(A))
        assert out["a"] == A and out["b"] == B
        assert out["c"] == [c ^ mulmod(a, b, f) for a, b, c in zip(A, B, Cv)]
        cnt = C.counts()
        assert cnt["AND"] == C.meta["rank"] and cnt["T"] == 4 * cnt["AND"] and cnt["qubits"] == 3 * m + 1


@pytest.mark.parametrize("mode", ["strict", "clean"])
def test_multiply_exact_degree2_small(mode):
    for m in range(3, 11):
        for alg in ("auto", "crt"):
            C = G.gf2m_multiply(m, mode=mode, algorithm=alg)
            ok, why = G.verify_multiplier(C, exact=True)
            assert ok, (m, alg, why)


def test_multiply_other_modulus_and_seed():
    f = (1 << 8) | (1 << 4) | (1 << 3) | (1 << 2) | 1          # another irreducible of degree 8
    assert is_irreducible(f)
    for seed in (0, 3):
        C = G.gf2m_multiply(8, f=f, algorithm="crt", seed=seed, n_anc=2)
        assert C.nq == 26
        assert G.verify_multiplier(C, exact=True)[0]
    assert G.gf2m_multiply(8, f=f).qasm() == G.gf2m_multiply(8, f=f).qasm()      # pure / deterministic


def test_multiply_canaries():
    for mode in ("strict", "clean"):
        C = G.gf2m_multiply(6, mode=mode, algorithm="crt")
        assert G.verify_multiplier(C, exact=True)[0]
        for kind, k in (("drop_cx", 0), ("drop_cx", 7), ("swap_cx", 3), ("drop_and", 0), ("drop_and", 5)):
            D = mutate(C, kind, k)
            assert not G.verify_multiplier(D, exact=True)[0], (mode, kind, k)


def test_multiply_prove_small():
    from qoptct.verify import prove
    for m in (3, 4):
        f = G.default_modulus(m)
        qc = write_tmp(schoolbook_qc(m, f), ".qc")
        for mode in ("strict", "clean"):
            C = G.gf2m_multiply(m, mode=mode, algorithm="crt")
            r = prove(C.qasm(), qc, clean=(mode == "clean"), budget_sec=120)
            assert r["status"] == "PROVED", (m, mode, r)
        D = mutate(G.gf2m_multiply(m, mode="strict", algorithm="crt"), "drop_and", 1)
        assert prove(D.qasm(), qc, clean=False, budget_sec=120)["status"] == "REFUTED"


@pytest.mark.slow
@pytest.mark.parametrize("mode", ["strict", "clean"])
def test_multiply_163_exact(mode):
    C = G.gf2m_multiply(163, mode=mode)
    cnt = C.counts()
    assert cnt["AND"] == 908 and cnt["T"] == 3632 and cnt["qubits"] == 490
    ok, why = G.verify_multiplier(C, lanes=64)
    assert ok, why
    ok, why = G.verify_multiplier(C, exact=True)          # all inputs (degree <= 2 argument)
    assert ok, why
    from qoptct.contrib.gf2m.qcirc import parse_qasm
    nq, ops = parse_qasm(C.qasm())                         # the emitted text is exactly the checked circuit
    assert nq == C.nq and ops == C.ops
    D = mutate(C, "swap_cx", 11)
    assert not G.verify_multiplier(D, lanes=16)[0]


# ---------------------------------------------------------------------------------------------- squaring

def test_square_small_exhaustive():
    for m in range(2, 12):
        f = G.default_modulus(m)
        for relabel in (False, True):
            C = G.gf2m_square(m, relabel=relabel)
            assert C.counts()["AND"] == 0 and C.nq == m
            xs = list(range(1 << m))
            st = [0] * m
            for i in range(m):
                st[i] = sum(1 << k for k, x in enumerate(xs) if x >> i & 1)
            C.run_lanes(st, (1 << len(xs)) - 1)
            ow = C.regs["out"]
            for k, x in enumerate(xs):
                y = sum((st[ow[i]] >> k & 1) << i for i in range(m))
                assert y == mulmod(x, x, f)
            if not relabel:
                assert ow == list(range(m))


def test_square_canary():
    from qoptct.contrib.gf2m import square
    C = G.gf2m_square(11)
    ops = [(o[1], o[2]) for o in C.ops]
    assert square.check_linear(ops, 11, G.default_modulus(11))
    assert not square.check_linear(ops[1:], 11, G.default_modulus(11))


@pytest.mark.slow
def test_square_nist_counts():
    exact = {163: 798, 233: 825, 283: 1386, 409: 1418, 571: 3239}
    relab = {163: 408, 233: 153, 283: 713, 409: 204, 571: 1538}
    for m in exact:
        assert G.gf2m_square(m).counts()["CX"] == exact[m]
        assert G.gf2m_square(m, relabel=True).counts()["CX"] == relab[m]


# ---------------------------------------------------------------------------------------------- Itoh-Tsujii

def test_chains():
    for m in range(2, 80):
        for ch in (itoh.it_chain(m - 1), itoh.short_chain(m - 1)):
            assert itoh.make_chain(m, ch)[-1][0] == m - 1
    assert len(itoh.it_chain(162)) - 1 == 9 and len(itoh.short_chain(162)) - 1 == 9
    assert len(itoh.short_chain(408)) - 1 == 10 and len(itoh.short_chain(570)) - 1 == 12


# lab audit x5 (k6 circuits, all pairs, own Euclid).  m = 4, 5, 6, 8 use the same decomposition here and must
# agree exactly; for m = 7 'auto' picks the CRT algorithm (same rank 22, other terms), fewer degenerate terms.
DIV_AND = {4: 27, 5: 38, 6: 56, 7: 88, 8: 120}
INV_AND = {4: 18, 5: 25, 6: 41, 7: 66, 8: 96}
SAME_DEC = (4, 5, 6, 8)


def test_divide_small_all_pairs():
    for m in range(3, 9):
        f = G.default_modulus(m)
        C = G.gf2m_divide(m)
        X, Y, T = [], [], []
        for x in range(1 << m):
            for y in range(1 << m):
                X.append(x)
                Y.append(y)
                T.append((x * 7 + y) % (1 << m))
        out = C.simulate(dict(x=X, y=Y, t=T), len(X))
        assert out["x"] == X and out["y"] == Y
        assert out["t"] == [t ^ (mulmod(y, inv(x, m, f), f) if x else 0) for x, y, t in zip(X, Y, T)]
        cost = C.meta["cost"]
        assert cost["mult"] == cost["chain_length"] + 1
        assert C.counts()["AND"] == cost["AND"] <= cost["formula_AND"]
        if m in SAME_DEC:
            assert cost["AND"] == DIV_AND[m]
        elif m in DIV_AND:
            assert cost["AND"] <= DIV_AND[m]


def test_itoh_lowerings_all_inputs():
    """CRT lowering (in-place compiler, copies for x^(2^k) * x) and parity lowering: same function, same ANDs."""
    for m in range(3, 9):
        f = G.default_modulus(m)
        X, Y, T = [], [], []
        for x in range(1 << m):
            for y in range(1 << m):
                X.append(x)
                Y.append(y)
                T.append((5 * x + 3 * y) % (1 << m))
        want = [t ^ (mulmod(y, inv(x, m, f), f) if x else 0) for x, y, t in zip(X, Y, T)]
        ands = set()
        for low in ("crt", "parity"):
            C = G.gf2m_divide(m, algorithm="crt", lowering=low)
            out = C.simulate(dict(x=X, y=Y, t=T), len(X))
            assert out["x"] == X and out["y"] == Y and out["t"] == want, (m, low)
            ands.add(C.counts()["AND"])
            I = G.gf2m_inverse(m, algorithm="crt", lowering=low)
            o2 = I.simulate(dict(x=list(range(1 << m)), t=[0] * (1 << m)), 1 << m)
            assert o2["t"] == [inv(x, m, f) if x else 0 for x in range(1 << m)], (m, low)
        assert len(ands) == 1


def test_inverse_small_all_inputs():
    for m in range(3, 9):
        f = G.default_modulus(m)
        C = G.gf2m_inverse(m)
        X = list(range(1 << m))
        T = [0] * len(X)
        out = C.simulate(dict(x=X, t=T), len(X))
        assert out["x"] == X and out["t"] == [inv(x, m, f) if x else 0 for x in X]
        assert C.meta["cost"]["mult"] == C.meta["cost"]["chain_length"]
        if m in SAME_DEC:
            assert C.counts()["AND"] == INV_AND[m]
        elif m in INV_AND:
            assert C.counts()["AND"] <= INV_AND[m]


def test_divide_fewer_registers():
    m = 12
    f = G.default_modulus(m)
    full = G.gf2m_divide_count(m)
    L = full["chain_length"]
    assert full["mult"] == L + 1
    with pytest.raises(ValueError):
        G.gf2m_divide(m, p=1)
    C = G.gf2m_divide(m, p=L - 2)                # fewer chain registers -> recomputation
    assert C.meta["cost"]["mult"] > full["mult"]
    assert C.meta["cost"]["peak_registers"] < full["peak_registers"]
    rng = random.Random(4)
    X = [rng.getrandbits(m) for _ in range(64)]
    Y = [rng.getrandbits(m) for _ in range(64)]
    out = C.simulate(dict(x=X, y=Y, t=[0] * 64), 64)
    assert out["t"] == [mulmod(y, inv(x, m, f), f) if x else 0 for x, y in zip(X, Y)]


def test_itoh_canaries():
    m = 5
    C = G.gf2m_divide(m)
    X, Y = [], []
    for x in range(1 << m):
        for y in range(1 << m):
            X.append(x)
            Y.append(y)
    f = G.default_modulus(m)
    want = [mulmod(y, inv(x, m, f), f) if x else 0 for x, y in zip(X, Y)]
    caught = 0
    kinds = (("drop_and", 0), ("drop_and", 9), ("drop_cx", 0), ("drop_cx", 40), ("swap_cx", 5), ("swap_cx", 77))
    for kind, k in kinds:
        D = mutate(C, kind, k)
        try:
            out = D.simulate(dict(x=X, y=Y, t=[0] * len(X)), len(X))
            caught += out["t"] != want or out["x"] != X or out["y"] != Y
        except PreconditionError:
            caught += 1
    assert caught == len(kinds)
    # a register erasure with a wrong phase polynomial is invisible to the classical simulation; the prover
    # below catches it


def _drop_phase(C):
    """canary: remove one gate of the first register-erasure phase correction."""
    ops = list(C.ops)
    i = next(i for i, o in enumerate(ops) if o[0] == 'p' or (o[0] == 'r' and o[2]))
    if ops[i][0] == 'p':
        del ops[i]
    else:
        a, poly, const = ops[i][1:]
        ops[i] = ('r', a, poly[1:], const)
    return GF2Circuit(C.nq, ops, C.regs, C.meta)


@pytest.mark.parametrize("erase,lowering,alg", [("terms", "crt", "crt"), ("terms", "parity", "auto"),
                                                ("monomials", "parity", "auto")])
def test_itoh_prove_small(erase, lowering, alg):
    """exact prover: the measurement-based circuit (all outcomes, phases) equals its Toffoli-level reference;
    the reference's classical map is the field division / inversion (checked on all inputs above)."""
    from qoptct.verify import prove
    for m, kind in ((3, "divide"), (3, "inverse"), (4, "inverse")):
        kw = dict(erase=erase, lowering=lowering, algorithm=alg)
        C = (G.gf2m_divide(m, **kw) if kind == "divide" else G.gf2m_inverse(m, **kw))
        qc = write_tmp(C.reference_qc(), ".qc")
        r = prove(C.qasm(), qc, clean=True, budget_sec=150)
        assert r["status"] == "PROVED", (m, kind, r)
        assert prove(_drop_phase(C).qasm(), qc, clean=True, budget_sec=150)["status"] != "PROVED"


# ------------------------------------------------- register erasure by terms: exact proofs at larger m (0.5.0)

def table_ref_qc(C, kind, m, f):
    """INDEPENDENT reference from the field definition: t ^= x^-1 (inverse) or t ^= y * x^-1 (divide), one
    multi-controlled X per (value v of x, output bit) with x^-1 = x^(2^m - 2) by square-and-multiply; wire names
    q<i> of the circuit, inputs = the registers x, (y,) t."""
    F = itoh.Field(m, f)
    nm = ["q%d" % i for i in range(C.nq)]
    X, T, Y = C.regs["x"], C.regs["t"], C.regs.get("y")
    ins = sorted(w for ws in C.regs.values() for w in ws)
    L = [".v " + " ".join(nm), ".i " + " ".join(nm[i] for i in ins), "BEGIN"]
    ctl = " ".join(nm[w] for w in X)
    for v in range(1, 1 << m):
        zero = [X[i] for i in range(m) if not v >> i & 1]
        L += ["tof %s" % nm[w] for w in zero]
        iv = F.inv(v)
        for j in range(m if kind == "divide" else 1):
            col = F.mul(iv, 1 << j) if kind == "divide" else iv
            for k in range(m):
                if col >> k & 1:
                    L.append("tof %s %s%s" % (ctl, (nm[Y[j]] + " ") if kind == "divide" else "", nm[T[k]]))
        L += ["tof %s" % nm[w] for w in zero]
    L.append("END")
    return "\n".join(L) + "\n"


def _drop_cphase(C, k):
    ops = list(C.ops)
    idx = [i for i, o in enumerate(ops) if o[0] == 'p']
    del ops[idx[k % len(idx)]]
    return GF2Circuit(C.nq, ops, C.regs, C.meta)


def _kprove_itoh(kind, m, tmp_path, max_mb=1400, canary=True):
    """kprove (BDD engine, all inputs and all measurement outcomes) of the default circuit (erase='terms', CRT
    lowering) against the independent table reference; canaries: one conditional phase removed -> REFUTED."""
    from qoptct.verify.kprove import kprove
    f = G.default_modulus(m)
    C = G.gf2m_divide(m, algorithm="crt") if kind == "divide" else G.gf2m_inverse(m, algorithm="crt")
    assert C.counts()["CZ"] + C.counts()["Z"] > 0
    q, r = tmp_path / "c.qasm", tmp_path / "r.qc"
    q.write_text(C.qasm())
    r.write_text(table_ref_qc(C, kind, m, f))
    order = list(C.regs["t"]) + list(C.regs.get("y", [])) + list(C.regs["x"])      # x on top of the BDD order
    kw = dict(clean=True, engine="bdd", sample=False, budget_sec=600, max_mb=max_mb, order=order)
    res = kprove(str(q), str(r), **kw)
    assert res["status"] == "PROVED", (kind, m, res.get("detail"), res.get("mem"))
    if canary:
        for k in (0, 7):
            d = tmp_path / "d.qasm"
            d.write_text(_drop_cphase(C, k).qasm())
            assert kprove(str(d), str(r), **kw)["status"] == "REFUTED", (kind, m, k)
    return res


@pytest.mark.parametrize("kind,m", [("inverse", 3), ("inverse", 4), ("inverse", 5), ("inverse", 6), ("inverse", 7),
                                    ("inverse", 8), ("divide", 3), ("divide", 4), ("divide", 5), ("divide", 6)])
def test_itoh_kprove_table_reference(kind, m, tmp_path):
    _kprove_itoh(kind, m, tmp_path)


@pytest.mark.slow
@pytest.mark.parametrize("kind,m", [("inverse", 10), ("inverse", 12), ("divide", 7), ("divide", 8), ("divide", 9),
                                    ("divide", 10)])
def test_itoh_kprove_table_reference_big(kind, m, tmp_path):
    """0.1-0.6 GB of BDD (growth budget 1.4 GB); larger proofs (inverse m = 14, divide m = 11: ~1.4 GB / ~0.75 GB)
    are in the 0.5.0 lab report, inverse m = 16 and divide m = 12 do not fit into 1.4 GB."""
    _kprove_itoh(kind, m, tmp_path)


def _phase_checked(C, vals, lanes):
    from qoptct.contrib.gf2m.qcirc import ErasurePhaseChecker
    ch = ErasurePhaseChecker()
    for o in C.ops:
        ch.feed(o)
    ok, why = ch.result()
    out = C.simulate(vals, lanes, zero_at=ch.assumed_zero)      # the checker's |0> assumptions hold
    return ok, why, out, ch


@pytest.mark.parametrize("lowering", ["crt", "parity"])
def test_erasure_phase_check_any_m(lowering):
    """symbolic exact check of every register-erasure phase correction (qcirc.ErasurePhaseChecker); canaries."""
    from qoptct.contrib.gf2m.qcirc import erasure_phase_check
    for m in (3, 4, 5, 6, 8, 11, 16, 23):
        f = G.default_modulus(m)
        F = itoh.Field(m, f)
        for kind in ("inverse", "divide"):
            fn = G.gf2m_divide if kind == "divide" else G.gf2m_inverse
            C = fn(m, lowering=lowering, algorithm="crt" if lowering == "crt" else "auto")
            rng = random.Random(m)
            X = list(range(1 << m)) if m <= 8 else [0, 1] + [rng.getrandbits(m) for _ in range(30)]
            vals = dict(x=X, t=[rng.getrandbits(m) for _ in X])
            if kind == "divide":
                vals["y"] = [rng.getrandbits(m) for _ in X]
            ok, why, out, ch = _phase_checked(C, vals, len(X))
            assert ok and ch.n_groups == C.meta["cost"]["register_erasures"], (m, kind, why)
            want = [t ^ (F.mul(vals["y"][k] if kind == "divide" else 1, F.inv(x)) if x else 0)
                    for k, (x, t) in enumerate(zip(X, vals["t"]))]
            assert out["t"] == want
            ps = [i for i, o in enumerate(C.ops) if o[0] == 'p']
            for j in (0, len(ps) // 2, len(ps) - 1):              # one conditional phase removed
                ops = list(C.ops)
                del ops[ps[j]]
                assert not erasure_phase_check(ops)[0], (m, kind, j)
            ops = list(C.ops)                                      # one CZ moved to another wire
            i2 = next(i for i in ps if len(ops[i][2]) == 2)
            mid, ws = ops[i2][1:]
            other = next(w for w in C.regs["x"] + C.regs["t"] if w not in ws)
            ops[i2] = ('p', mid, (ws[0], other))
            assert not erasure_phase_check(ops)[0], (m, kind)


@pytest.mark.slow
@pytest.mark.parametrize("kind", ["divide", "inverse"])
def test_erasure_phase_check_163(kind):
    """the full gate-level circuit at m = 163 (~1.8M ops, ~0.25 GB): every one of the 668 520 (divide) conditional
    phases of the register erasures is checked symbolically; classical part on 16 inputs."""
    m = 163
    f = G.default_modulus(m)
    F = itoh.Field(m, f)
    C = G.gf2m_divide(m) if kind == "divide" else G.gf2m_inverse(m)
    rng = random.Random(3)
    X = [0, 1, (1 << m) - 1] + [rng.getrandbits(m) for _ in range(13)]
    vals = dict(x=X, t=[rng.getrandbits(m) for _ in X])
    if kind == "divide":
        vals["y"] = [rng.getrandbits(m) for _ in X]
    ok, why, out, ch = _phase_checked(C, vals, len(X))
    assert ok and ch.n_groups == 9 and ch.n_phase == C.meta["cost"]["CZ_register_erasures"], why
    want = [t ^ (F.mul(vals["y"][k] if kind == "divide" else 1, F.inv(x)) if x else 0)
            for k, (x, t) in enumerate(zip(X, vals["t"]))]
    assert out["t"] == want


# ------------------------------------------------- exact2 as a prover in qoptct.verify (0.5.0)

def test_exact2_prover_multiplier(tmp_path):
    from qoptct.verify.exact2 import prove_degree2, gadget_check
    from qoptct.verify.kprove import kprove
    assert gadget_check() < 1e-12
    for m in (3, 5, 8):
        f = G.default_modulus(m)
        for mode in ("strict", "clean"):
            C = G.gf2m_multiply(m, mode=mode, algorithm="crt")
            qc = tmp_path / ("s%d%s.qc" % (m, mode))
            txt = schoolbook_qc(m, f)
            if mode == "strict":          # c is an input as well
                txt = txt.replace(".i " + " ".join(["a%d" % i for i in range(m)] + ["b%d" % i for i in range(m)]),
                                  ".i " + " ".join(["a%d" % i for i in range(m)] + ["b%d" % i for i in range(m)]
                                                   + ["c%d" % i for i in range(m)]))
            qc.write_text(txt)
            q = tmp_path / "m.qasm"
            q.write_text(C.qasm())
            r = prove_degree2(str(q), str(qc), clean=(mode == "clean"))
            assert r["status"] == "PROVED", (m, mode, r)
            assert kprove(str(q), str(qc), clean=(mode == "clean"))["status"] == "PROVED"
            for kind, k in (("drop_cx", 0), ("drop_cx", 9), ("swap_cx", 3)):
                q.write_text(mutate(C, kind, k).qasm())
                assert prove_degree2(str(q), str(qc), clean=(mode == "clean"))["status"] != "PROVED", (m, mode, kind)
    # not applicable: the division (degree > 2) -> no verdict
    D = G.gf2m_divide(3)
    q = tmp_path / "d.qasm"
    q.write_text(D.qasm())
    qc = tmp_path / "d.qc"
    qc.write_text(D.reference_qc())
    assert prove_degree2(str(q), str(qc), clean=True)["status"] in ("ERROR", "UNKNOWN")


@pytest.mark.slow
@pytest.mark.parametrize("mode", ["strict", "clean"])
def test_exact2_prover_163(mode, tmp_path):
    from qoptct.verify.exact2 import prove_degree2
    m = 163
    f = G.default_modulus(m)
    C = G.gf2m_multiply(m, mode=mode)
    names = ["a%d" % i for i in range(m)] + ["b%d" % i for i in range(m)] + ["c%d" % i for i in range(m)]
    txt = schoolbook_qc(m, f)
    if mode == "strict":
        txt = txt.replace(".i " + " ".join(names[:2 * m]), ".i " + " ".join(names))
    qc = tmp_path / "s.qc"
    qc.write_text(txt)
    q = tmp_path / "m.qasm"
    q.write_text(C.qasm())
    r = prove_degree2(str(q), str(qc), clean=(mode == "clean"))
    assert r["status"] == "PROVED", r
    n = 3 * m if mode == "strict" else 2 * m           # free input bits
    assert r["lanes"] == 1 + n * (n + 1) // 2
    q.write_text(mutate(C, "drop_cx", 1234).qasm())
    assert prove_degree2(str(q), str(qc), clean=(mode == "clean"))["status"] != "PROVED"


def test_itoh_erase_modes_agree():
    for m in (5, 7):
        a = G.gf2m_divide(m, erase="terms")
        b = G.gf2m_divide(m, erase="monomials")
        assert a.counts()["AND"] == b.counts()["AND"]
        X = list(range(1 << m)) * 4
        Y = [(7 * i + 3) % (1 << m) for i in range(len(X))]
        oa = a.simulate(dict(x=X, y=Y, t=[0] * len(X)), len(X))
        ob = b.simulate(dict(x=X, y=Y, t=[0] * len(X)), len(X))
        assert oa == ob


@pytest.mark.parametrize("m,chain,mult,and_", [(163, "it", 10, 9080), (233, "it", 11, 14773),
                                               (283, "it", 12, 20028)])
def test_divide_counts_nist(m, chain, mult, and_):
    d = G.gf2m_divide_count(m, chain=chain)
    assert d["mult"] == mult and d["formula_AND"] == and_ and d["AND"] <= and_
    i = G.gf2m_inverse_count(m, chain=chain)
    assert i["mult"] == mult - 1 and i["formula_AND"] == and_ - NIST_MULT_AND[m]


def test_stream_gate_check_matches_circuit():
    for kind in ("divide", "inverse"):
        C = G.gf2m_divide(9) if kind == "divide" else G.gf2m_inverse(9)
        s = G.gf2m_itoh_gate_check(kind, 9, lanes=40)
        cnt = C.counts()
        assert s["gates"]["AND"] == cnt["AND"] == s["AND"]
        assert s["gates"]["CX_lin"] == cnt["CX_lin"]


@pytest.mark.slow
@pytest.mark.parametrize("kind", ["divide", "inverse"])
@pytest.mark.parametrize("lowering", ["crt", "parity"])
def test_itoh_163_gate_level(kind, lowering):
    s = G.gf2m_itoh_gate_check(kind, 163, lanes=24, lowering=lowering)
    assert s["gates"]["AND"] == s["AND"] == (9080 if kind == "divide" else 8172)
    if lowering == "crt":
        assert s["gates"]["CX"] < 1.3e6


@pytest.mark.slow
@pytest.mark.parametrize("m,mult", [(409, 11), (571, 13)])
def test_divide_counts_nist_short_chain(m, mult):
    d = G.gf2m_divide_count(m, chain="short")
    assert d["mult"] == mult and d["formula_AND"] == mult * NIST_MULT_AND[m]
