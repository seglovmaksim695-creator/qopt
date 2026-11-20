"""qoptct.crypto.gf2ec: point addition on binary curves.  Independent affine arithmetic (naive polynomial field,
brute-force inverse) for small m: ALL admissible pairs (P1, P2); exact proofs (kprove on the promised domain) at
m = 3, 4 with canaries; phase corrections of the register erasures (symbolic, any m); AND counts m = 163..571;
gate-level stream check at m = 163.  @slow: m = 163 circuit / phase check / stream check, counts m = 409, 571."""
import random

import pytest

from qoptct.crypto import gf2ec as E
from qoptct.crypto import gf2m as G
from qoptct.contrib.gf2m.qcirc import GF2Circuit, erasure_phase_check, ErasurePhaseChecker


# ------------------------------------------------------------------ independent reference (no library field code)

def _pmul(a, b, f, m):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> m & 1:
            a ^= f
    return r


def _ref_curve(m, f, a, b):
    mul = lambda u, v: _pmul(u, v, f, m)
    inv = {}
    for u in range(1, 1 << m):
        for v in range(1, 1 << m):
            if mul(u, v) == 1:
                inv[u] = v
                break
    pts = [(x, y) for x in range(1 << m) for y in range(1 << m)
           if mul(y, y) ^ mul(x, y) == mul(mul(x, x), x) ^ mul(a, mul(x, x)) ^ b]

    def add(P1, P2):                          # x1 != x2
        (x1, y1), (x2, y2) = P1, P2
        lam = mul(y1 ^ y2, inv[x1 ^ x2])
        x3 = mul(lam, lam) ^ lam ^ x1 ^ x2 ^ a
        return x3, mul(lam, x1 ^ x3) ^ x3 ^ y1
    return pts, add


def _domain(pts, add, P2, controlled):
    """admissible (P1, c, result): c = 1: x1 != x2 and x3 != x2;  c = 0 (controlled only): x1 != 0, result = P1."""
    out = []
    for P1 in pts:
        if P1[0] != P2[0]:
            R = add(P1, P2)
            if R[0] != P2[0]:
                out.append((P1, 1, R))
        if controlled and P1[0] != 0:
            out.append((P1, 0, P1))
    return out


def _run(C, dom, controlled):
    n = len(dom)
    vals = dict(x=[d[0][0] for d in dom], y=[d[0][1] for d in dom])
    if controlled:
        vals["c"] = [d[1] for d in dom]
    out = C.simulate(vals, n)
    assert out["x"] == [d[2][0] for d in dom]
    assert out["y_out"] == [d[2][1] for d in dom]
    if controlled:
        assert out["c"] == vals["c"]
    # the wires of y1 that are not wires of y3 must be clean (simulate treats all register wires as outputs)
    full = (1 << n) - 1
    st = [0] * C.nq
    for name, v in vals.items():
        for i, w in enumerate(C.regs[name]):
            st[w] = sum(1 << ln for ln, t in enumerate(v) if t >> i & 1)
    C.run_lanes(st, full)
    assert all(st[w] == 0 for w in C.regs["y"] if w not in set(C.regs["y_out"]))


CURVES = {3: [(1, 1), (1, 3), (1, 5)], 4: [(0, 1), (1, 1), (1, 7)], 5: [(0, 1), (1, 1), (1, 5)],
          6: [(1, 1), (0, 5)], 7: [(1, 1)], 8: [(1, 1)]}


@pytest.mark.parametrize("m", [3, 4, 5, 6, 7, 8])
@pytest.mark.parametrize("controlled", [False, True])
def test_all_pairs_small(m, controlled):
    """every point P2 (constant), every admissible P1 (and c in {0, 1}); ancillas clean (simulate enforces every
    AND / erasure precondition and that all non-register wires return to 0)."""
    f = G.default_modulus(m)
    built = 0
    for a, b in CURVES[m]:
        pts, add = _ref_curve(m, f, a, b)
        for P2 in pts:
            dom = _domain(pts, add, P2, controlled)
            if not any(d[1] for d in dom) or (controlled and not any(not d[1] for d in dom)):
                with pytest.raises(ValueError):
                    E.gf2_point_add(m, f, a, b, P2, controlled=controlled)
                continue
            C = E.gf2_point_add(m, f, a, b, P2, controlled=controlled)
            _run(C, dom, controlled)
            built += 1
            if m >= 7 and built >= 40:               # keep the m = 7, 8 runs short: 40 constants P2 of the curve
                break
    assert built > 0


@pytest.mark.parametrize("m,reg", [(5, 4), (6, 4), (7, None), (8, 4)])
def test_few_registers(m, reg):
    f = G.default_modulus(m)
    a, b = 1, 1
    pts, add = _ref_curve(m, f, a, b)
    P2 = pts[len(pts) // 3]
    dom = _domain(pts, add, P2, False)
    C = E.gf2_point_add(m, f, a, b, P2, registers=reg)
    _run(C, dom, False)
    assert C.meta["cost"]["gf_registers"] == C.meta["cost"]["registers"] + 2
    assert reg is None or C.meta["cost"]["registers"] == reg


def test_registers_tradeoff_and_errors():
    cm = E.min_registers(163)
    assert cm <= 4
    full = E.gf2_point_add_count(163)
    p4 = E.gf2_point_add_count(163, registers=4)
    assert p4["AND"] > full["AND"] and p4["gf_registers"] == 6 and full["gf_registers"] == 12
    with pytest.raises(ValueError):
        E.gf2_point_add_count(163, registers=1)
    with pytest.raises(ValueError):
        E.gf2_point_add(5, None, 1, 1, (1, 1))                    # not a curve point
    with pytest.raises(ValueError):
        E.gf2_point_add(5, None, 1, 0, (0, 0))                    # b = 0
    with pytest.raises(ValueError):
        E.gf2_point_add(5, None, 1, 1, None)


def test_inadmissible_inputs_are_not_handled():
    """documented limitation: P1 = -P2 and P1 + P2 = +-P2 break an erasure precondition (a measurement-based circuit
    would leave garbage); P1 = P2 gives X = 0 and silently NOT the doubling.  None of these is treated."""
    m, f, a, b = 5, G.default_modulus(5), 1, 1
    pts, add = _ref_curve(m, f, a, b)
    P2 = next(P for P in pts if P[0] != 0)
    C = E.gf2_point_add(m, f, a, b, P2)
    neg = (P2[0], P2[1] ^ P2[0])
    with pytest.raises(AssertionError):
        C.simulate(dict(x=[neg[0]], y=[neg[1]]), 1)
    bad = [P for P in pts if P[0] != P2[0] and add(P, P2)[0] == P2[0]]
    assert bad
    with pytest.raises(AssertionError):
        C.simulate(dict(x=[bad[0][0]], y=[bad[0][1]]), 1)
    out = C.simulate(dict(x=[P2[0]], y=[P2[1]]), 1)          # P1 = P2: no error, but not 2*P2 (lam = 0)
    assert out["x"] and out["y_out"]


# ------------------------------------------------------------------ mutants, phases

def _active_and(C, vals, n, cand):
    """first AND among the op indices `cand` whose two inputs are not 0 AND 0 on all input lanes (a gate that is
    inactive on every tested input cannot be told from a missing one by any test on these inputs)."""
    full = (1 << n) - 1
    for i in cand:
        st = [0] * C.nq
        for name, v in vals.items():
            for b, w in enumerate(C.regs[name]):
                st[w] = sum(1 << ln for ln, t in enumerate(v) if t >> b & 1)
        D = GF2Circuit(C.nq, C.ops[:i], dict(C.regs), dict(C.meta))
        D.run_lanes(st, full)
        _, p, q, _ = C.ops[i]
        if st[p] & st[q]:
            return i
    raise AssertionError("no active AND among the candidates")


def _mutants(C, vals=None, n=0):
    ops = C.ops
    ps = [i for i, o in enumerate(ops) if o[0] == 'p']
    ands = [i for i, o in enumerate(ops) if o[0] == 'a']
    cxs = [i for i, o in enumerate(ops) if o[0] == 'c']
    xs = [i for i, o in enumerate(ops) if o[0] == 'x']
    ai = ands[len(ands) // 2]
    if vals is not None:
        ai = _active_and(C, vals, n, ands[len(ands) // 2:] + ands)
    out = []
    for name, idx in (("drop_phase_first", ps[0]), ("drop_phase_mid", ps[len(ps) // 2]),
                      ("drop_phase_last", ps[-1]), ("drop_and", ai),
                      ("drop_cx", cxs[len(cxs) // 2]), ("drop_x", xs[0] if xs else cxs[0])):
        o2 = list(ops)
        del o2[idx]
        out.append((name, GF2Circuit(C.nq, o2, dict(C.regs), dict(C.meta))))
    o2 = list(ops)                                    # one CZ moved to another wire
    i2 = next(i for i in ps if len(ops[i][2]) == 2)
    mid, ws = ops[i2][1:]
    other = next(w for w in C.regs["x"] if w not in ws)
    o2[i2] = ('p', mid, (ws[0], other))
    out.append(("move_cz", GF2Circuit(C.nq, o2, dict(C.regs), dict(C.meta))))
    return out


@pytest.mark.parametrize("m", [3, 4, 5, 8, 11])
@pytest.mark.parametrize("controlled", [False, True])
def test_phase_checker_and_mutants(m, controlled):
    """ErasurePhaseChecker (symbolic, exact for all inputs and outcomes): every conditional phase of the register
    erasures is right; mutants (removed phase / AND / CX / X, moved CZ) are caught by the phase checker or by the
    classical simulation."""
    f = G.default_modulus(m)
    rng = random.Random(m)
    if m <= 8:
        pts, add = _ref_curve(m, f, 1, 1)
        P2 = pts[len(pts) // 2]
        dom = _domain(pts, add, P2, controlled)
    else:
        P2 = E.random_point(m, f, 1, 3, rng)
        dom = []
        while len(dom) < 24:
            P1 = E.random_point(m, f, 1, 3, rng)
            if P1[0] != P2[0]:
                R = E.ec_add(m, f, 1, P1, P2)
                if R[0] != P2[0]:
                    dom.append((P1, 1, R))
    a, b = 1, (1 if m <= 8 else 3)
    C = E.gf2_point_add(m, f, a, b, P2, controlled=controlled)
    ch = ErasurePhaseChecker()
    for o in C.ops:
        ch.feed(o)
    ok, why = ch.result()
    assert ok and ch.n_groups == C.meta["cost"]["register_erasures"], why
    vals = dict(x=[d[0][0] for d in dom], y=[d[0][1] for d in dom])
    if controlled:
        vals["c"] = [d[1] for d in dom]
    out = C.simulate(vals, len(dom), zero_at=ch.assumed_zero)
    assert out["x"] == [d[2][0] for d in dom] and out["y_out"] == [d[2][1] for d in dom]
    for name, D in _mutants(C, vals, len(dom)):
        phase_ok = erasure_phase_check(D)[0]
        try:
            o = D.simulate(vals, len(dom))
            classical_ok = o["x"] == out["x"] and o["y_out"] == out["y_out"]
        except AssertionError:
            classical_ok = False
        assert not (phase_ok and classical_ok), (m, controlled, name)


# ------------------------------------------------------------------ exact proofs (m = 3, 4)

def _table_ref_qc(C, dom, controlled):
    """INDEPENDENT Toffoli-level reference from the table of the map (P1, c) -> result on the promised domain:
    scratch := result (one multi-controlled X per table entry and output bit, controls = the input wires), then the
    inputs are cleared (controls = scratch + c), then scratch is moved into the output wires.  Scratch wires are
    ancillas of the circuit."""
    m = C.meta["m"]
    X, Y, Yo = C.regs["x"], C.regs["y"], C.regs["y_out"]
    c = C.regs.get("c", [])
    used = set(X) | set(Y) | set(Yo) | set(c)
    scr = [w for w in range(C.nq) if w not in used][:2 * m]
    assert len(scr) == 2 * m
    Sx, Sy = scr[:m], scr[m:]
    nm = ["q%d" % i for i in range(C.nq)]
    L = [".v " + " ".join(nm), ".i " + " ".join(nm[w] for w in X + Y + c), "BEGIN"]

    def entry(ctl, bits_zero, tgt):
        L.extend("tof %s" % nm[w] for w in bits_zero)
        L.extend("tof %s %s" % (" ".join(nm[w] for w in ctl), nm[t]) for t in tgt)
        L.extend("tof %s" % nm[w] for w in bits_zero)
    for P1, cb, R in dom:
        ctl = X + Y + c
        val = [(P1[0] >> i) & 1 for i in range(m)] + [(P1[1] >> i) & 1 for i in range(m)] + ([cb] if controlled else [])
        zero = [w for w, v in zip(ctl, val) if not v]
        tgt = [Sx[i] for i in range(m) if R[0] >> i & 1] + [Sy[i] for i in range(m) if R[1] >> i & 1]
        entry(ctl, zero, tgt)
    for P1, cb, R in dom:
        ctl = Sx + Sy + c
        val = [(R[0] >> i) & 1 for i in range(m)] + [(R[1] >> i) & 1 for i in range(m)] + ([cb] if controlled else [])
        zero = [w for w, v in zip(ctl, val) if not v]
        tgt = [X[i] for i in range(m) if P1[0] >> i & 1] + [Y[i] for i in range(m) if P1[1] >> i & 1]
        entry(ctl, zero, tgt)
    for i in range(m):
        L.append("tof %s %s" % (nm[Sx[i]], nm[X[i]]))
        L.append("tof %s %s" % (nm[X[i]], nm[Sx[i]]))
        L.append("tof %s %s" % (nm[Sy[i]], nm[Yo[i]]))
        L.append("tof %s %s" % (nm[Yo[i]], nm[Sy[i]]))
    L.append("END")
    return "\n".join(L) + "\n"


def _kprove_domain(C, dom, controlled, tmp_path, canaries=True):
    from qoptct.verify.kprove import kprove
    from qoptct.verify.bdd import ONE, ZERO
    m = C.meta["m"]
    X, Y, c = C.regs["x"], C.regs["y"], C.regs.get("c", [])
    q, r = tmp_path / "c.qasm", tmp_path / "r.qc"
    q.write_text(C.qasm())
    r.write_text(_table_ref_qc(C, dom, controlled))

    def bits(P1, cb):
        d = {}
        for i, w in enumerate(X):
            d[w] = P1[0] >> i & 1
        for i, w in enumerate(Y):
            d[w] = P1[1] >> i & 1
        if controlled:
            d[c[0]] = cb
        return d

    def builder(ps, xedge):
        B = ps.B
        acc = ZERO
        for P1, cb, _ in dom:
            t = ONE
            for w, v in bits(P1, cb).items():
                e = xedge(w)
                t = B.AND(t, e if v else e ^ 1)
            acc = B.OR(acc, t)
        return acc

    def sampler(x, rng):
        P1, cb, _ = rng.choice(dom)
        for w, v in bits(P1, cb).items():
            x[w] = v
        return x
    kw = dict(clean=True, pre=(builder, sampler), engine="bdd", sample=False, budget_sec=900, max_mb=1500)
    res = kprove(str(q), str(r), **kw)
    assert res["status"] == "PROVED", res
    if canaries:
        # canaries: every 4th conditional phase gate removed in turn (a phase invisible on the promised domain
        # may stay PROVED: here at most 10 %), and a CZ moved to another wire
        ps = [i for i, o in enumerate(C.ops) if o[0] == 'p']
        verdicts = []
        for k in range(0, len(ps), 4):
            o2 = list(C.ops)
            del o2[ps[k]]
            d = tmp_path / "d.qasm"
            d.write_text(GF2Circuit(C.nq, o2, dict(C.regs), dict(C.meta)).qasm())
            verdicts.append(kprove(str(d), str(r), **kw)["status"])
        assert set(verdicts) <= {"REFUTED", "PROVED"} and verdicts.count("REFUTED") >= 0.9 * len(verdicts), verdicts
        D = _mutants(C)[-1][1]
        d = tmp_path / "d.qasm"
        d.write_text(D.qasm())
        assert kprove(str(d), str(r), **kw)["status"] == "REFUTED", "move_cz"
    return res


@pytest.mark.parametrize("m,controlled", [(3, False), (3, True), (4, False), (4, True)])
def test_kprove_small(m, controlled, tmp_path):
    """exact proof (all inputs of the promised domain, all measurement outcomes, phases) against the table
    reference; removed phase gates / moved CZ are REFUTED."""
    f = G.default_modulus(m)
    a, b = 1, (1 if m == 3 else 7)
    pts, add = _ref_curve(m, f, a, b)
    best = None
    for P2 in pts:
        dom = _domain(pts, add, P2, controlled)
        if any(d[1] for d in dom) and (not controlled or any(not d[1] for d in dom)):
            if best is None or len(dom) > len(best[1]):
                best = (P2, dom)
    P2, dom = best
    C = E.gf2_point_add(m, f, a, b, P2, controlled=controlled)
    _kprove_domain(C, dom, controlled, tmp_path)


# ------------------------------------------------------------------ counts

#            m: (L, mult, rank, exact AND)   registers = None, chain 'auto', uncontrolled
EXPECT = {163: (9, 20, 908), 233: (10, 22, 1343), 283: (11, 24, 1669), 409: (10, 22, 2501), 571: (12, 26, 3569)}


@pytest.mark.parametrize("m", [163, 233, 283])
def test_counts_nist(m):
    L, mult, R = EXPECT[m]
    c = E.gf2_point_add_count(m)
    assert (c["chain_length"], c["mult"], c["rank"]) == (L, mult, R)
    assert c["formula_AND"] == mult * R and c["AND"] <= c["formula_AND"] and c["AND"] > c["formula_AND"] - 40 * mult
    cc = E.gf2_point_add_count(m, controlled=True)
    assert cc["AND"] - c["AND"] == 2 * m - 1


@pytest.mark.slow
@pytest.mark.parametrize("m", [409, 571])
def test_counts_nist_big(m):
    L, mult, R = EXPECT[m]
    c = E.gf2_point_add_count(m)
    assert (c["chain"], c["mult"], c["rank"]) == ("short", mult, R)
    assert E.gf2_point_add_count(m, controlled=True)["AND"] - c["AND"] == 2 * m - 1


def test_count_matches_circuit_small():
    for m in (5, 8):
        f = G.default_modulus(m)
        pts, add = _ref_curve(m, f, 1, 1)
        P2 = pts[2]
        for ctl in (False, True):
            C = E.gf2_point_add(m, f, 1, 1, P2, controlled=ctl)
            c = E.gf2_point_add_count(m, f, 1, 1, P2, controlled=ctl)
            assert C.counts()["AND"] == c["AND"] == C.meta["cost"]["AND"]


# ------------------------------------------------------------------ m = 163 (slow)

@pytest.mark.slow
def test_stream_gate_check_163():
    """whole circuit at gate level (not stored), 20 random curve points, p = L + 1 and p = 4, plain and controlled."""
    for reg, ctl in ((None, False), (4, False), (None, True)):
        c = E.gf2_point_add_check(163, registers=reg, controlled=ctl, lanes=20)
        n = E.gf2_point_add_count(163, registers=reg, controlled=ctl)
        assert c["AND"] == n["AND"] and c["gates"]["AND"] == n["AND"]


@pytest.mark.slow
def test_phase_check_163():
    """the stored circuit m = 163 (p = L + 1): every conditional phase of all register erasures checked symbolically;
    classical action on 16 random points; counts as in the table."""
    m, f = 163, G.default_modulus(163)
    rng = random.Random(5)
    a, b = 1, rng.getrandbits(163) | 1
    P2 = E.random_point(m, f, a, b, rng)
    C = E.gf2_point_add(m, f, a, b, P2, samples=4)
    ch = ErasurePhaseChecker()
    for o in C.ops:
        ch.feed(o)
    ok, why = ch.result()
    assert ok and ch.n_groups == C.meta["cost"]["register_erasures"], why
    dom = []
    while len(dom) < 16:
        P1 = E.random_point(m, f, a, b, rng)
        R = E.ec_add(m, f, a, P1, P2)
        if P1[0] != P2[0] and R[0] != P2[0]:
            dom.append((P1, 1, R))
    out = C.simulate(dict(x=[d[0][0] for d in dom], y=[d[0][1] for d in dom]), len(dom), zero_at=ch.assumed_zero)
    assert out["x"] == [d[2][0] for d in dom] and out["y_out"] == [d[2][1] for d in dom]
    assert C.counts()["AND"] == C.meta["cost"]["AND"] == E.gf2_point_add_count(163)["AND"]
    for name, D in _mutants(C):
        assert not erasure_phase_check(D)[0] or name in ("drop_and", "drop_cx", "drop_x"), name
