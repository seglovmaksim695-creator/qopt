# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""qoptct.crypto.gf2ec -- point addition on binary elliptic curves over GF(2^m) (affine, second point a constant)
[qoptct 0.6; lab k6 gf2ec.py: point_add2 / point_add2_ctrl, rebuilt on qoptct.crypto.gf2m / contrib.gf2m].

    gf2_point_add(m, f=None, a=0, b=1, P2=(x2, y2), registers=None, controlled=False, chain='auto')
        circuit (x1, y1) -> (x3, y3) = (x1, y1) + (x2, y2) on  y^2 + x*y = x^3 + a*x^2 + b  over GF(2^m) = F2[x]/(f)
    gf2_point_add_count(...)        the same costs without building the circuit (m = 163 ... 571 in seconds)
    gf2_point_add_check(...)        gate-level check of the whole circuit on random curve points WITHOUT storing it
                                    (every AND / erasure precondition enforced; for m = 163 ... 571)
    ec_add(m, f, a, P1, P2), on_curve(...), curve_points(...), random_point(...)   classical reference arithmetic

Method (affine formulas, P2 = (x2, y2) a classical constant):
    X = x1 + x2,  Y = y1 + y2,  lam = Y / X                       (division: Itoh-Tsujii, gf2m.py)
    x3 = lam^2 + lam + x1 + x2 + a,   y3 = lam * (x3 + x2) + x3 + y2
    lam is erased by a second division  lam = Y' / X'  with X' = x3 + x2,  Y' = lam * X' (= y3 + x3 + y2).
Two divisions and two multiplications.  The registers Y (after lam is computed) and lam (after y3 is computed) are
erased WITHOUT further ANDs, by an X-basis measurement with a classically controlled CZ/Z phase correction
(qoptct.contrib.gf2m.itoh, 'erase=terms'): a product register is removed for free while its two factors are alive.
Inside the Itoh-Tsujii chains the same free erasure is used (pebble game with free removal).

Cost (logical ANDs = Toffoli equivalents; one product = R ANDs, R = rank of the bilinear algorithm, M = R):
    registers = None (p = L + 1 chain registers, L = length of the addition chain):  2 (L + 1) multiplications,
        e.g. m = 163 (L = 9, R = 908):  2 * 10 * 908 = 18 160 ANDs minus the few product terms that are linear
        (see gf2_point_add_count for the exact numbers), (2 + p) GF(2^m) registers
    registers = p smaller: more multiplications (recomputation in the chains); p = 4 gives 32 multiplications
        (29 056 ANDs at m = 163) with 6 registers
    controlled = True: the addition is applied only if a control qubit c is 1 (c = 0: identity); + (2m - 1) ANDs
        (two controlled linear maps: rank(S + I) = m - 1 and the identity, m), the control wire is regs['c'].
Measured at m = 163 (gate level, tests / lab report): lowering='crt' (default): p = L + 1: 18 160 ANDs, 1 957
qubits, 2.33M CX; p = 4: 29 056 ANDs, 1 142 qubits (7m + 1: a product of a register with itself works on a
temporary copy), 3.72M CX.  lowering='parity' (term by term on parity wires): 1 959 / 981 qubits ((2 + p) m + 3),
13.2M / 21.1M CX.  CX counts include the 3 CX of every AND gadget and the CZ of the erasures.
The multiplication count of each variant is exactly the pebble-plan optimum of the chain; chain='auto' takes the
better of the Itoh-Tsujii chain and a shortest star chain (409: 10 instead of 11 multiplications per division, 571:
12 instead of 13).

ADMISSIBLE INPUTS -- WHAT IS NOT HANDLED (read this before comparing numbers).
    The circuit is correct only for inputs with  x1 != x2  and  x3 != x2  (i.e. P1 != +-P2 and P1 + P2 != +-P2);
    equivalently the two divisions are by nonzero numbers.  The point at infinity O, doubling P1 = P2, P1 = -P2,
    and P1 + P2 = +-P2 are NOT treated: no case distinction and no flag qubit.  On P1 = -P2 and on P1 + P2 = +-P2
    an erasure precondition fails (the classical simulation raises; a measurement-based circuit would leave a
    wrong phase / garbage); on P1 = P2 the first division is by 0 and the circuit silently returns something that
    is not 2*P2.  The controlled
    version additionally needs x1 != 0 when c = 0.  This is the same convention as
      * Banegas, Bernstein, van Hoof, Lange (BBHL), TCHES 2021 (ePrint 2020/1180), Algorithm 3,
      * Taguchi-Takayasu (TT), ePrint 2023/553 (m = 163: 40 169 Toffoli per addition),
    so the comparison with these two is like for like (all three count ANDs/Toffoli of a point addition with
    a constant second point and no exceptional cases).  It is NOT like for like with Garn-Kan, arXiv:2503.02984
    (Table IV: 68 200 Toffoli at m = 163), whose point addition is exact with respect to all exceptional cases
    (four inversions); our number is smaller partly because those cases are dropped.  In a complete Shor / ECDLP
    algorithm the exceptional cases are usually handled by the choice of the accumulator offsets or are dropped
    with a negligible failure probability; that analysis is NOT part of this module and no total Shor cost is
    claimed here.
    Prior art.  Erasing a GF(2^m) product register by X-basis measurement + CZ correction is published (Khattar,
    Shutty, Gidney et al., arXiv:2510.10967, App. C.3, Qualtran GF2MulMBUC); here it is applied to the whole
    point-addition (second division for lam, the chains, y).  The affine formulas are those of BBHL / TT.

Wires.  GF2Circuit (qoptct.contrib.gf2m.qcirc): regs['x'] = wires of x1 (and of x3: x is updated in place),
regs['y'] = wires of y1 (input), regs['y_out'] = wires of y3.  The register of y is released after lam is
computed and allocated again later, so y3 is in general NOT on the wires of y1: after the circuit the wires
regs['y'] are 0 (unless they coincide with regs['y_out']).  regs['c'] = [control wire] if controlled.  Everything
else is ancilla and returns to |0>.  bit i of a register is on wire regs[name][i] (coefficient of t^i).
`C.simulate(dict(x=[...], y=[...]), lanes)` runs the classical action on many inputs at once (every AND / erasure
precondition is enforced); `C.qasm()` / `C.reference_qc()` give OpenQASM 3 and the Toffoli-level reference.

Verification (tests/test_crypto_gf2ec.py): all admissible pairs (P1, P2) on small fields m = 3..8 (independent
affine formula, ancillas clean); exact proofs (qoptct.verify.prove, kprove on the promised domain) at m = 3, 4 with
an independent table reference and canaries (removed phase gate -> REFUTED); the phase corrections of all register
erasures are checked symbolically for any m (qcirc.ErasurePhaseChecker; m = 163 in the slow tests); gate-level
stream check at m = 163 on random points.
"""
import random

from qoptct.contrib.gf2m import itoh as _itoh
from qoptct.crypto import gf2m as _G
from qoptct.crypto._errors import CheckError

__all__ = ["gf2_point_add", "gf2_point_add_count", "gf2_point_add_check", "ec_add", "on_curve", "curve_points",
           "random_point", "min_registers"]


# ------------------------------------------------------------------------------------------- classical curve

def _field(m, f=None):
    return _itoh.Field(m, _G._modulus(m, f))


def on_curve(m, f, a, b, P):
    """(x, y) on  y^2 + x y = x^3 + a x^2 + b  over GF(2^m) = F2[t]/(f)."""
    F = _field(m, f)
    x, y = P
    return F.sq(y) ^ F.mul(x, y) == F.mul(F.sq(x), x) ^ F.mul(a, F.sq(x)) ^ b


def ec_add(m, f, a, P1, P2):
    """affine addition formula for P1 != +-P2 (x1 != x2):  lam = (y1+y2)/(x1+x2), x3 = lam^2 + lam + x1 + x2 + a,
    y3 = lam (x1 + x3) + x3 + y1."""
    F = _field(m, f)
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2:
        raise ValueError("x1 == x2: not an admissible pair")
    lam = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
    x3 = F.sq(lam) ^ lam ^ x1 ^ x2 ^ a
    y3 = F.mul(lam, x1 ^ x3) ^ x3 ^ y1
    return x3, y3


def curve_points(m, f, a, b):
    """all affine points (small m only: loops over 4^m pairs)."""
    if m > 12:
        raise ValueError("curve_points is for small m")
    F = _field(m, f)
    pts = []
    for x in range(1 << m):
        c = F.mul(F.sq(x), x) ^ F.mul(a, F.sq(x)) ^ b
        for y in range(1 << m):
            if F.sq(y) ^ F.mul(x, y) == c:
                pts.append((x, y))
    return pts


def random_point(m, f, a, b, rng):
    """a random affine point with x != 0 (odd m: half-trace; even m <= 16: search)."""
    F = _field(m, f)
    for _ in range(10000):
        x = rng.getrandbits(m)
        if x == 0:
            continue
        c = F.mul(F.inv(F.sq(x)), b) ^ x ^ a          # z = y / x solves z^2 + z = c
        if m % 2:
            h = t = c
            for _ in range((m - 1) // 2):
                t = F.sq(F.sq(t))
                h ^= t
            if F.sq(h) ^ h == c:
                return x, F.mul(x, h ^ rng.getrandbits(1))
        else:
            if m > 16:
                raise ValueError("even m > 16 not supported by random_point")
            zs = [z for z in range(1 << m) if F.sq(z) ^ z == c]
            if zs:
                return x, F.mul(x, rng.choice(zs))
    raise ValueError("no point found")


# ------------------------------------------------------------------------------------------- chain choice

def _make_chain(m, registers, chain):
    """-> itoh.Chain with p = registers chain registers (default L + 1: all chain registers and the sink alive) and
    use_p = p - 1 (the sink state keeps one register free for lam / y).  chain 'auto': the better of 'it' / 'short'."""
    names = ("it", "short") if chain == "auto" else (chain,)
    best = None
    for nm in names:
        ch = _itoh.make_chain(m, nm)
        L = len(ch) - 1
        p = L + 1 if registers is None else registers
        r = _itoh.pebble_plan(ch, p, p - 1) if p >= 2 else None
        if r is None:
            continue
        cost = r[0] + r[1]
        if best is None or cost < best[0]:
            best = (cost, nm, ch, p)
    if best is None:
        raise ValueError("no pebbling strategy with registers=%r chain registers" % (registers,))
    _, nm, ch, p = best
    C = _itoh.Chain(m, ch, p, p - 1)
    C.name = nm
    return C


def min_registers(m, chain="auto"):
    """smallest number of chain registers p for which a strategy exists (use_p = p - 1)."""
    for p in range(2, 40):
        try:
            _make_chain(m, p, chain)
            return p
        except ValueError:
            continue
    raise ValueError("no strategy")


# ------------------------------------------------------------------------------------------- program primitives

def _addc(P, R, const):
    """R ^= const (X gates)."""
    P.val[R] = [v ^ const for v in P.val[R]]
    P._fresh.discard(R)
    if P.gates:
        for i in range(P.m):
            if const >> i & 1:
                P.C.x(P.wires[R][i])


def _caddc(P, R, c, const):
    """R ^= c * const (CX from the control wire)."""
    P.val[R] = [v ^ (const if b else 0) for v, b in zip(P.val[R], P.val[c])]
    P._fresh.discard(R)
    if P.gates:
        for i in range(P.m):
            if const >> i & 1:
                P.C.cx(P.wires[c][0], P.wires[R][i])


def _clin(P, T, A, ks, c):
    """T ^= c * (sum_k A^(2^k)); ANDs = rank of the linear map (one AND of c with a parity per basis row)."""
    F, m = P.F, P.m

    def Lv(v):
        r = 0
        for k in ks:
            r ^= F.frob(v, k)
        return r
    P.val[T] = [t ^ (Lv(a) if b else 0) for t, a, b in zip(P.val[T], P.val[A], P.val[c])]
    P._fresh.discard(T)
    rows = [0] * m
    for j in range(m):
        col = Lv(1 << j)
        for i in range(m):
            if col >> i & 1:
                rows[i] |= 1 << j
    basis, expr = [], []
    for r in rows:
        e, x = 0, r
        for bi, (pv, bv) in enumerate(basis):
            if x >> pv & 1:
                x ^= bv
                e ^= 1 << bi
        if x:
            basis.append((x.bit_length() - 1, x))
            e ^= 1 << (len(basis) - 1)
        expr.append(e)
    P.and_count += len(basis)
    P.clin_and = getattr(P, "clin_and", 0) + len(basis)
    if P.gates:
        C = P.C
        wa, wt, cw = P.wires[A], P.wires[T], P.wires[c][0]
        for bi, (pv, bv) in enumerate(basis):
            p, idx = P._par(wa, bv)
            y = P._new()
            C.AND(cw, p, y)
            for i in range(m):
                if expr[i] >> bi & 1:
                    C.cx(y, wt[i])
            C.erase(cw, p, y)
            P._release(y)
            P._unpar(p, idx, wa)


def _program(P, ch, x2, y2, a, controlled):
    """the register program of the addition (k6 point_add2 / point_add2_ctrl)."""
    X, Y = "x", "y"
    if not controlled:
        _addc(P, X, x2)
        _addc(P, Y, y2)                                  # X = x1 + x2, Y = y1 + y2
        ch.run(P, X, ch.up, "c")
        P.alloc("lam")
        P.mul("lam", ch.sink("c"), 1, Y, 0)              # lam = Y / X
        P.erase(Y, "lam", 0, X, 0)                       # Y == lam * X: erased, register released
        ch.run(P, X, ch.down, "c")
        P.lin(X, "lam", 1)
        P.lin(X, "lam", 0)                               # X = x1 + x2 + lam^2 + lam
        _addc(P, X, a ^ x2)                              # X = x3 + x2
        ch.run(P, X, ch.up, "d")
        P.alloc(Y)
        P.mul(Y, "lam", 0, X, 0)                         # Y = lam * (x3 + x2)
        P.erase("lam", ch.sink("d"), 1, Y, 0)            # lam == Y / X: erased
        ch.run(P, X, ch.down, "d")
        _addc(P, X, x2)                                  # X = x3
        P.lin(Y, X, 0)
        _addc(P, Y, y2)                                  # Y = lam (x3 + x2) + x3 + y2 = y3
    else:
        c = "ctl"
        _caddc(P, X, c, x2)
        _caddc(P, Y, c, y2)
        ch.run(P, X, ch.up, "c")
        P.alloc("lam")
        P.mul("lam", ch.sink("c"), 1, Y, 0)
        P.erase(Y, "lam", 0, X, 0)
        ch.run(P, X, ch.down, "c")
        _clin(P, X, "lam", (1, 0), c)                    # X ^= c (lam^2 + lam)
        _caddc(P, X, c, a ^ x2)                          # c = 1: x3 + x2;  c = 0: x1
        ch.run(P, X, ch.up, "d")
        P.alloc(Y)
        P.mul(Y, "lam", 0, X, 0)                         # c = 0: lam * x1 = y1
        P.erase("lam", ch.sink("d"), 1, Y, 0)
        ch.run(P, X, ch.down, "d")
        _caddc(P, X, c, x2)                              # x3 / x1
        _clin(P, Y, X, (0,), c)
        _caddc(P, Y, c, y2)                              # y3 / y1


# ------------------------------------------------------------------------------------------- samples, building

def _samples(m, f, a, b, P2, controlled, n, seed):
    """n admissible inputs (x1, y1, c, x3, y3), deterministic in seed; small m: all admissible points first."""
    rng = random.Random(seed)
    F = _field(m, f)
    out = []
    seen = set()

    def cands():
        if m <= 8:
            pts = curve_points(m, f, a, b)
            rng.shuffle(pts)
            for P1 in pts:
                yield P1
        else:
            while True:
                yield random_point(m, f, a, b, rng)
    tries = 0
    for P1 in cands():
        tries += 1
        if tries > 5000 + 10 * n:
            break
        if P1 in seen:
            continue
        if controlled and len(out) % 2 == 1:
            if P1[0] == 0:
                continue
            out.append((P1[0], P1[1], 0, P1[0], P1[1]))
            seen.add(P1)
        else:
            if P1[0] == P2[0]:
                continue
            R = ec_add(m, f, a, P1, P2)
            if R[0] == P2[0]:
                continue
            out.append((P1[0], P1[1], 1, R[0], R[1]))
            seen.add(P1)
        if len(out) >= n:
            break
    if not out or (controlled and len(out) < 2):
        raise ValueError("no admissible input points on this curve for P2 (too few points)")
    return out


def _build(m, f, a, b, P2, registers, controlled, chain, gates, nsamples, seed, algorithm, lowering):
    f = _G._modulus(m, f)
    if not 0 <= a < 1 << m:
        raise ValueError("a must be a field element")
    if not 0 < b < 1 << m:
        raise ValueError("b must be a nonzero field element (b = 0 is singular)")
    if P2 is None:
        raise ValueError("P2 = (x2, y2), a point of the curve, is required")
    x2, y2 = P2
    if not (0 <= x2 < 1 << m and 0 <= y2 < 1 << m) or not on_curve(m, f, a, b, (x2, y2)):
        raise ValueError("P2 is not a point of the curve")
    if lowering not in ("auto", "crt", "parity"):
        raise ValueError("lowering must be 'auto', 'crt' or 'parity'")
    dec, info = _G.gf2m_bilinear(m, f, algorithm)
    if lowering == "crt" and info["algorithm"] != "crt":
        raise ValueError("lowering='crt' needs the CRT algorithm")
    use_crt = bool(gates) and info["algorithm"] == "crt" and lowering in ("auto", "crt")
    ch = _make_chain(m, registers, chain)
    smp = _samples(m, f, a, b, (x2, y2), controlled, nsamples, seed)
    P = _itoh.Program(m, f, dec, dict(x=[s[0] for s in smp], y=[s[1] for s in smp]), gates=gates, count_cz=False,
                      erase_mode="terms",
                      crt=_G._crt_data(m, f, info["places"]) if use_crt else None)
    y_in = list(P.wires["y"])
    if controlled:
        w = P._new()
        P.val["ctl"] = [s[2] for s in smp]
        P.wires["ctl"] = [w]
        if P.stream:
            P.C.write([w], P.val["ctl"])
    _program(P, ch, x2, y2, a, controlled)
    want_live = {"x", "y"} | ({"ctl"} if controlled else set())
    if set(P.live_regs) != {"x", "y"} or set(P.val) != want_live:
        raise AssertionError("chain registers left alive")
    if P.val["x"] != [s[3] for s in smp] or P.val["y"] != [s[4] for s in smp]:
        raise AssertionError("register program gives a wrong result")
    P.finish()
    R = len(dec)
    cost = dict(mult=P.n_mul, AND=P.and_count, T=4 * P.and_count, rank=R, formula_AND=P.n_mul * R + getattr(P, "clin_and", 0),
                chain_length=ch.L, chain=ch.name, registers=ch.p, gf_registers=ch.p + 2,
                register_erasures=P.n_erase, qubits_registers=(ch.p + 2) * m + (1 if controlled else 0),
                peak_wires=P.peak_wires, controlled=bool(controlled), algorithm=info["algorithm"], samples=len(smp))
    return P, cost, y_in, smp


def gf2_point_add(m, f=None, a=0, b=1, P2=None, registers=None, controlled=False, chain="auto", samples=8, seed=1,
                  algorithm="auto", lowering="auto"):
    """circuit (x1, y1) -> (x1, y1) + P2 on  y^2 + x y = x^3 + a x^2 + b  over GF(2^m) = F2[t]/(f)  (f: int, bit i =
    coefficient of t^i; default gf2m.default_modulus(m)); P2 = (x2, y2) a point of the curve (constant).
    registers: chain registers p (default None = L + 1, fewest multiplications; at least min_registers(m));
    total (2 + p) GF(2^m) registers.  controlled: add a control wire regs['c'] (c = 0: identity).
    chain: 'auto' (better of 'it' and 'short'), 'it', 'short' or an explicit chain.
    samples / seed: number of random admissible inputs (curve points) on which the register program is checked
    while building (the circuit itself does not depend on them).
    -> qoptct.contrib.gf2m.qcirc.GF2Circuit with regs 'x', 'y', 'y_out' (and 'c'), meta['cost'].
    ONLY ADMISSIBLE INPUTS x1 != x2, x3 != x2 are handled (see the module docstring)."""
    P, cost, y_in, smp = _build(m, f, a, b, P2, registers, controlled, chain, True, samples, seed, algorithm, lowering)
    C = P.C
    C.nq = P._nw
    C.regs = dict(x=list(P.wires["x"]), y=y_in, y_out=list(P.wires["y"]))
    if controlled:
        C.regs["c"] = list(P.wires["ctl"])
    cost["gates"] = C.counts()
    C.meta.update(kind="gf2_point_add", m=m, f=_G._modulus(m, f), a=a, b=b, P2=tuple(P2), controlled=bool(controlled),
                  cost=cost, algorithm=cost["algorithm"])
    return C


def gf2_point_add_count(m, f=None, a=1, b=None, P2=None, registers=None, controlled=False, chain="auto", samples=4,
                        seed=1, algorithm="auto"):
    """costs of gf2_point_add without building gates (m = 163 ... 571 in seconds): exact AND count (nonzero bilinear
    terms of every product + the controlled linear maps), multiplications, registers; the register program is run
    on `samples` random curve points (every erasure condition and the result are checked).  b, P2 default to a
    random curve (b != 0) and a point on it."""
    f = _G._modulus(m, f)
    rng = random.Random(seed)
    if b is None:
        b = rng.getrandbits(m) or 1
    if P2 is None:
        P2 = random_point(m, f, a, b, rng)
    P, cost, _, _ = _build(m, f, a, b, P2, registers, controlled, chain, False, samples, seed, algorithm, "auto")
    cost.pop("peak_wires")
    return cost


def gf2_point_add_check(m, f=None, a=1, b=None, P2=None, registers=None, controlled=False, chain="auto", lanes=16,
                        seed=1, algorithm="auto", lowering="auto"):
    """gate-level check for large m without storing the circuit: every gate of every product and of the parity wires
    is executed bitsliced on `lanes` random admissible inputs (curve points), every AND / erasure precondition
    enforced (register erasures lazily, their phases: see ErasurePhaseChecker); the output is compared with the
    classical formula and all ancillas must be 0.  -> cost dict with 'gates' counts and 'gate_check'; raises
    on any failure."""
    f = _G._modulus(m, f)
    rng = random.Random(seed)
    if b is None:
        b = rng.getrandbits(m) or 1
    if P2 is None:
        P2 = random_point(m, f, a, b, rng)
    P, cost, y_in, smp = _build(m, f, a, b, P2, registers, controlled, chain, "stream", lanes, seed, algorithm, lowering)
    S = P.C
    if S.read(P.wires["x"]) != [s[3] for s in smp] or S.read(P.wires["y"]) != [s[4] for s in smp]:
        raise AssertionError("gate-level stream: wrong result")
    keep = set(P.wires["x"]) | set(P.wires["y"]) | (set(P.wires["ctl"]) if controlled else set())
    if any(v for w, v in S.st.items() if w not in keep):
        raise AssertionError("gate-level stream: ancilla not clean")
    cost["gates"] = S.counts()
    cost["gate_check"] = "%d random admissible inputs at gate level (every AND / erasure precondition)" % len(smp)
    return cost
