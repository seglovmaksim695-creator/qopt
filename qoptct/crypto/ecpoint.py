# SPDX-License-Identifier: Apache-2.0
"""ecpoint: EXACT (complete) mixed addition of an elliptic-curve point over a prime field, as a reversible circuit on
logical ANDs (Gidney 2018: AND = 4 T + 3 CX, erased by an X-basis measurement with a Clifford correction, 0 T).
Ported from the lab agent k4 (work/k4: padd_full / padd_k4 / padd_win), quantum-arithmetic building block for
resource estimation and verification.

    from qoptct.crypto import ecpoint
    C = ecpoint.ec_point_add(p, a, b, Q=(xq, yq), karatsuba_levels=3)     # (P) -> (P + Q), Q classical
    C.counts()                    # {'AND': ..., 'CX': ..., 'T': ..., 'qubits': ..., ...}
    C.simulate([None, Q, (x, y)]) # exact classical bit-sliced run, every AND / erase checked, ancillas must end at 0
    W = ecpoint.ec_point_add_window(p, a, b, table, w)                      # (P, k) -> (P + T[k], k), k = w bits
    W.simulate([(P, k), ...])

What the circuit computes. The register (x, y) holds a point of  y^2 = x^3 + a x + b  over F_p, or the point at
infinity O encoded as (0, 0) (needs a curve of ODD group order: no point has y = 0; checked at build time). The
circuit maps  P -> P + Q  for EVERY P (O included) and a classical Q != O of order > 4: no exceptional case is
dropped (P = O, P = Q, P = -Q, P = -2Q are handled; the others are the generic chord formula). Because Q is classical,
the answers on the exceptional set E = {O, Q, -Q, -2Q} are constants (Q, 2Q, O, -Q): flags are computed (about 6n
ANDs), P is replaced by a safe point D = 3Q (CX only), the generic addition runs, the image is corrected (CX only) and
the flags are erased from the result. The generic addition has two divisions (Kaliski inversion + multiplication,
computed, applied, and uncomputed), 8 multiplications and 2 squarings in "compute-apply-uncompute" form, in the
ordinary (non-Montgomery) representation; the divisions use x*2^(2n) mod p (sparse constant) and the Kaliski inverse
with the "halfB" saving; products use conditionally inverted additions (signed digits 2y_i - 1; n^2 + 4n ANDs) with
Karatsuba levels (garbage uncomputed by the reversal of the whole block) and the reduction of qoptct.crypto.fold
(2 744 ANDs P-256, 2 175 secp256k1) for special primes.

Costs (n = 256; ANDs = Toffoli-equivalents, worst case, not average; CX = explicit CX + 3 per AND + 1 per erase):

    curve        karatsuba_levels   ANDs        ANDs/n^2   qubits (peak wires)
    secp256k1    3                  3 033 954   46.3       6 934
    secp256k1    2                  3 074 732   46.9       5 730
    secp256k1    0                  3 286 052   50.1       4 560
    P-256        3                  3 041 421   46.4       6 742   (T = 12.17 M, CX = 27.5 M)
    P-256        2 / 0              3 082 199 / 3 293 519         5 538 / 4 368
    multiplier='school' (partial-product ANDs) instead of 'sd': P-256 L3 3 261 829 (49.8 n^2)
    window (table lookup of Q, w = 4, field-lookups of 9 values + 4 flags) adds about 18 n ANDs and 9n qubits;
    the lookup itself (2 * (2^w - 2) ANDs by minterms) is counted in the circuit.

Honest comparison. Same cost model (logical AND = Toffoli equivalent, EXACT circuits, worst case): Litinski 2023
(arXiv:2306.08585) / Qualtran (Papa 2025, arXiv:2506.03318) need about 126.5 n^2 (about 8.3 M ANDs) per exact
addition. Our exact circuit is about 63 % below that. Approximate schemes of 2026 are cheaper and use fewer qubits and
this circuit is NOT a record against them: Haener et al. arXiv:2609.05625 about 1.2 M ANDs (failure probability > 0,
exceptional cases not handled), Schrottenloher arXiv:2606.02235 1.86 M ANDs at 1446 qubits, Initio 2026. We use
4.4 - 6.9 thousand qubits (depending on karatsuba_levels). Approximate and exact schemes must not be mixed in one
table: p_fail > 0 versus an exact (error-free) map.

Choice of construction (honest). The Kaliski inverse and the multipliers of this package are NOT the blocks of
qoptct.crypto.kaliski / mulsq: the lab k4 circuit needs a "compute / apply / uncompute" builder with shared garbage
wires and a specific halfB variant of the Kaliski step, which the stand-alone (own-engine) generators of the library do
not expose as embeddable sub-circuits. The k4 core was therefore ported as is (subpackage _ecp, own engine from the lab
modules arith / circk), and only the reduction of special primes is taken from the library (qoptct.crypto.fold, through
an adapter). The library kaliski_inverse has 12 n^2 - 8 n - 3 ANDs against 12 n^2 + n here (4 passes: a saving of
about 36 n = 9 thousand ANDs, 0.3 %, not applied). The numbers above are those of the lab (K_CLAIMS rows k4),
reproduced exactly.

Checks (tests/test_crypto_ecpoint.py): all P (O included) x all Q on curves of prime order n = 4..8 (two curves per
n), every AND target and erase condition, ancillas clean; random and exceptional points n = 16..64; n = 256
(secp256k1, P-256) 16+ points including O, Q, -Q, 2Q, -2Q; the window variant on all (P, k); phases: sparse
quantum simulation of the emitted QASM on a superposition (n = 4) and mutants rejected. An exact proof of the whole
circuit by qoptct.verify.kprove is out of reach (no independent reference of the whole point addition; the BDD of a
multiplier already grows exponentially at n = 10): sub-blocks of this kind are proved in qoptct.crypto.mulsq / kaliski.
"""
import contextlib
import random

from ._errors import CheckError
from ._ecp import engine as _E
from ._ecp import ec as _ec
from ._ecp import padd as _padd
from ._ecp import full as _full
from ._ecp import win as _win
from ._ecp import emit as _emit
from ._ecp.ec import ec_add, ec_neg, ec_mul, on_curve, is_prime, P256, SECP256K1

__all__ = ['ec_point_add', 'ec_point_add_window', 'ECPointAddCircuit', 'CURVES', 'ec_add', 'ec_neg', 'ec_mul',
           'has_two_torsion', 'find_point']

CURVES = {'P-256': P256, 'secp256k1': SECP256K1}


# ---------------------------------------------------------------- classical helpers
def _sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def find_point(p, a, b, rng=None):
    """a random point (x, y) with y != 0 on y^2 = x^3 + a x + b over F_p."""
    rng = rng or random.Random(1)
    while True:
        x = rng.randrange(p)
        y = _sqrt_mod(x ** 3 + a * x + b, p)
        if y:
            return (x, y if rng.random() < .5 else p - y)


def _pmul(f, g, h, p):
    """polynomials over F_p as coefficient lists (low first), product mod h (monic cubic)."""
    r = [0] * 5
    for i, u in enumerate(f):
        for j, v in enumerate(g):
            r[i + j] = (r[i + j] + u * v) % p
    for d in (4, 3):
        c = r[d]
        if c:
            for k in range(3):
                r[d - 3 + k] = (r[d - 3 + k] - c * h[k]) % p
            r[d] = 0
    return r[:3]


def _trim(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def _polymod(u, v, p):
    """u mod v over F_p (coefficient lists, low first, trimmed)."""
    u = _trim(list(u)); inv = pow(v[-1], -1, p)
    while len(u) >= len(v):
        c = u[-1] * inv % p; sft = len(u) - len(v)
        for i, t in enumerate(v):
            u[sft + i] = (u[sft + i] - c * t) % p
        _trim(u)
    return u


def has_two_torsion(p, a, b):
    """True iff x^3 + a x + b has a root in F_p, i.e. the group has a point with y = 0 (even order):
    gcd(x^p - x, x^3 + a x + b) != 1."""
    h = [b % p, a % p, 0]                      # x^3 = -(h0 + h1 x + h2 x^2) modulo f
    res, base, e = [1, 0, 0], [0, 1, 0], p
    while e:
        if e & 1:
            res = _pmul(res, base, h, p)
        base = _pmul(base, base, h, p); e >>= 1
    g = _trim([res[0], (res[1] - 1) % p, res[2]])      # x^p - x mod f
    u, v = [b % p, a % p, 0, 1], g
    while v:
        u, v = v, _polymod(u, v, p)
    return len(u) > 1                          # deg gcd >= 1 (g = 0: f divides x^p - x)


# ---------------------------------------------------------------- build parameters
_BASES = {'sd': 'sd', 'school': 'school'}
_REDS = {'fold': 'k2', 'solinas': 'sol'}


@contextlib.contextmanager
def _modes(multiplier, reduction):
    old = dict(_padd.MUL), dict(_padd.RED)
    _padd.MUL['base'] = _BASES[multiplier]
    _padd.RED['mode'] = _REDS[reduction]
    try:
        yield
    finally:
        _padd.MUL.clear(); _padd.MUL.update(old[0])
        _padd.RED.clear(); _padd.RED.update(old[1])


def _check_args(p, a, b, n, L, multiplier, reduction):
    if not (isinstance(p, int) and p > 3 and p % 2 == 1 and is_prime(p)):
        raise ValueError('p must be an odd prime')
    if n is None:
        n = p.bit_length()
    if n != p.bit_length():
        raise ValueError('n must be the bit length of p (%d)' % p.bit_length())
    if multiplier not in _BASES:
        raise ValueError("multiplier must be 'sd' or 'school'")
    if reduction not in _REDS:
        raise ValueError("reduction must be 'fold' or 'solinas'")
    if not (isinstance(L, int) and 0 <= L <= 6):
        raise ValueError('karatsuba_levels must be an integer in 0..6')
    if (4 * a ** 3 + 27 * b * b) % p == 0:
        raise ValueError('singular curve')
    if has_two_torsion(p, a, b):
        raise ValueError('the group has even order (a point with y = 0): O = (0, 0) is not representable; '
                         'use a curve of odd order (P-256, secp256k1, ...)')
    return n


def _check_point(P, p, a, b, what='point'):
    if P is None:
        return
    if not (0 <= P[0] < p and 0 <= P[1] < p) or (P[1] * P[1] - P[0] ** 3 - a * P[0] - b) % p:
        raise ValueError('%s is not on the curve' % what)


# ---------------------------------------------------------------- circuit objects
class ECPointAddCircuit:
    """A built circuit (P) -> (P + Q) [window: (P, k) -> (P + T[k], k)].

    Attributes: p, a, b, n, Q (classical) / table, w (window), D (safe replacement point), karatsuba_levels,
    multiplier, reduction, circ (engine circuit: ops, wires), X, Y (the n + n input/output wires, LSB first), K (the
    w index wires, window only), window (bool).
    """

    def __init__(s, p, a, b, n, Q, table, w, D, L, multiplier, reduction):
        s.p, s.a, s.b, s.n, s.Q, s.table, s.w, s.D = p, a, b, n, Q, table, w, D
        s.karatsuba_levels, s.multiplier, s.reduction = L, multiplier, reduction
        s.window = table is not None
        cv = dict(p=p, a=a, b=b)
        with _modes(multiplier, reduction):
            C = _E.CircK(); X = C.reg('x', n); Y = C.reg('y', n)
            if s.window:
                K = C.reg('k', w)
                _win.point_add_window(C, X, Y, K, table, cv, D, L)
                s.K = K
            else:
                generic = _padd.make(p, L)
                info = _full.point_add_full(C, X, Y, Q, cv, lambda v: v % p, generic, D=D)
                s.D = info['D']; s.K = []
        s.circ, s.X, s.Y = C, X, Y

    # -- counts
    def counts(s):
        """dict: AND (= Toffoli equivalents, 4 T each), T, CX (explicit + 3 per AND + 1 per erase), ERASE, qubits
        (all wires ever allocated by the builder = peak of simultaneously live wires), n."""
        k = s.circ.counts()
        return dict(AND=k['AND'], CCX=k['CCX'], ERASE=k['ERASE'], T=k['T'], CX=k['CX'], qubits=k['qubits'],
                    CXops=k['CXops'], n=s.n)

    # -- classical, exact
    def _code(s, P):
        return (0, 0) if P is None else (P[0], P[1])

    def simulate(s, points, ks=None, check_curve=True):
        """Exact bit-sliced run on a batch of inputs. points: list of P (None = O, else (x, y) on the curve); window:
        ks = list of indices (same length). Checks every AND target = 0 and every erase condition, and that all
        ancillas end at 0 (CheckError otherwise). Returns the list of results (None for O); window: the index
        register is also checked to be unchanged."""
        S = len(points)
        if s.window:
            if ks is None or len(ks) != S:
                raise ValueError('window circuit: pass ks, one table index per point')
        elif ks is not None:
            raise ValueError('ks is for the window circuit only')
        for P in points:
            if check_curve:
                _check_point(P, s.p, s.a, s.b)
            if P is not None and P[1] == 0:
                raise ValueError('a point with y = 0 cannot be encoded')
        init = _E.pack([s._code(P)[0] for P in points], s.X)
        init.update(_E.pack([s._code(P)[1] for P in points], s.Y))
        if s.window:
            init.update(_E.pack(list(ks), s.K))
        v = _E.simulate(s.circ, init, S)
        gx, gy = _E.unpack(v, s.X, S), _E.unpack(v, s.Y, S)
        io = set(s.X) | set(s.Y) | set(s.K)
        dirty = [w for w in range(s.circ.nw) if w not in io and v[w]]
        if dirty:
            raise CheckError('%d ancilla wires are not 0 at the end' % len(dirty))
        if s.window:
            if _E.unpack(v, s.K, S) != list(ks):
                raise CheckError('the index register changed')
        return [None if (x, y) == (0, 0) else (x, y) for x, y in zip(gx, gy)]

    def reference(s, P, k=None):
        """classical P + Q (window: P + T[k])."""
        Q = s.table[k] if s.window else s.Q
        return ec_add(P, Q, s.p, s.a)

    def check(s, points, ks=None):
        """simulate + compare with the classical sum. Returns the number of wrong results (0 = correct)."""
        out = s.simulate(points, ks)
        return sum(1 for i, (P, R) in enumerate(zip(points, out))
                   if R != s.reference(P, None if ks is None else ks[i]))

    def exceptional_points(s, k=None):
        """the points that exercise the special branches: O, Q, -Q, 2Q, -2Q, D, -D (for window: index k)."""
        Q = s.table[k] if s.window else s.Q
        Q2 = ec_add(Q, Q, s.p, s.a)
        return [None, Q, ec_neg(Q, s.p), Q2, ec_neg(Q2, s.p), s.D, ec_neg(s.D, s.p)]

    # -- text / phase
    def qasm(s):
        """OpenQASM 3 text of the circuit (AND gadget of Gidney 2018 fig. 3, X-basis measurement erase)."""
        return _emit.to_qasm(s.circ)

    def phase_check(s, points, K=6, trials=4, seed=1, qasm=None):
        """sparse quantum simulation (independent of the circuit object, reads only the QASM text): a superposition
        of K basis inputs with random complex amplitudes, random measurement outcomes; returns the worst overlap
        |<expected|got>| over the trials (1.0 = exact up to a global phase). Classical (no window) only; small n."""
        if s.window:
            raise NotImplementedError('phase_check: classical-Q circuit only')
        txt = qasm or s.qasm()
        pos = lambda P: sum(((s._code(P)[0] >> i) & 1) << w for i, w in enumerate(s.X)) + \
            sum(((s._code(P)[1] >> i) & 1) << w for i, w in enumerate(s.Y))
        fun = lambda P: pos(ec_add(P, s.Q, s.p, s.a))
        return _emit.sparse_check(txt, pos, fun, list(points), K=K, trials=trials, seed=seed)


def ec_point_add(p, a, b, Q, n=None, karatsuba_levels=3, reduction='fold', multiplier='sd', D=None):
    """Build the exact mixed addition P -> P + Q (Q classical, Q != O of order > 4) on y^2 = x^3 + a x + b over F_p.

    p odd prime (n = bit length), the curve of odd group order (no point with y = 0). karatsuba_levels L (0..6; the
    qubit/AND trade-off, see the module table). reduction: 'fold' (qoptct.crypto.fold for p = 2^n - c and P-256;
    other primes fall back to the Solinas search of this package) or 'solinas'. multiplier: 'sd' (conditionally
    inverted additions, default) or 'school' (partial products). D: the safe replacement point (default 3Q).
    Returns an ECPointAddCircuit."""
    n = _check_args(p, a, b, n, karatsuba_levels, multiplier, reduction)
    _check_point(Q, p, a, b, 'Q')
    if Q is None or Q[1] == 0:
        raise ValueError('Q must be a point other than O')
    Q2 = ec_add(Q, Q, p, a)
    if Q2 is None or ec_add(Q2, Q2, p, a) is None or ec_neg(Q2, p) in (Q, ec_neg(Q, p)):
        raise ValueError('the order of Q must be greater than 4')
    Q = tuple(Q)
    if D is None:
        D = ec_add(Q2, Q, p, a)                                  # 3Q
        bad = {Q, ec_neg(Q, p), ec_neg(Q2, p), None}
        okD = lambda Dc: Dc not in bad and ec_add(Dc, Q, p, a) not in (Q, ec_neg(Q, p), None)
        if not okD(D):                                           # order 5 etc.: look for another safe point
            rng = random.Random(1)
            for _ in range(2000):
                Dc = find_point(p, a, b, rng)
                if okD(Dc):
                    D = Dc
                    break
            else:
                raise ValueError('no safe point D found; pass D explicitly')
    return ECPointAddCircuit(p, a, b, n, Q, None, 0, D, karatsuba_levels, multiplier, reduction)


def ec_point_add_window(p, a, b, table, w, n=None, karatsuba_levels=3, reduction='fold', multiplier='sd', D=None,
                        seed=1):
    """Windowed exact addition (P, k) -> (P + table[k], k) for a w-bit quantum index k (Shor's algorithm window).

    table: the 2^w classical points T[k] != O (the offset T[k] = k G + R0 is accounted for classically by the
    caller). D: a classical point safe for ALL k (D not in {T, -T, -2T}, D + T not in {T, -T, O}); by default it is
    searched (random points; for tiny curves the table must leave enough points). Lookup by minterms (2^w - 2 ANDs each
    way); the unlookup is a repeated XOR of the same table. Returns an ECPointAddCircuit (window = True)."""
    n = _check_args(p, a, b, n, karatsuba_levels, multiplier, reduction)
    table = [tuple(T) for T in table]
    if len(table) != (1 << w) or w < 1:
        raise ValueError('table must have 2^w entries, w >= 1')
    for T in table:
        _check_point(T, p, a, b, 'a table point')
        if T[1] == 0:
            raise ValueError('table points must differ from O and have y != 0')
        Q2 = ec_add(T, T, p, a)
        if Q2 is None or ec_neg(Q2, p) in (T, ec_neg(T, p)):
            raise ValueError('every table point must have order > 4')
    bad = set()
    for T in table:
        bad |= {T, ec_neg(T, p), ec_neg(ec_add(T, T, p, a), p)}
    def safe(Dc):
        return (Dc is not None and Dc not in bad and
                all(ec_add(Dc, T, p, a) not in (T, ec_neg(T, p), None) for T in table))
    if D is None:
        rng = random.Random(seed)
        for _ in range(2000):
            Dc = find_point(p, a, b, rng)
            if safe(Dc):
                D = Dc
                break
        else:
            raise ValueError('no safe point D found; pass D explicitly')
    elif not safe(tuple(D)):
        raise ValueError('D is not safe for every table entry')
    return ECPointAddCircuit(p, a, b, n, None, table, w, tuple(D), karatsuba_levels, multiplier, reduction)
