# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""crt: CRT (Chinese-remainder) bilinear algorithms for multiplication in GF(2^m) = F2[x]/(f)  [qoptct v0.5, lab k5/x3].

The full product P = a*b of two polynomials with m coefficients (2m-1 coefficients) is recovered from its residues
modulo pairwise coprime "places" h = p^e (p irreducible) and from its top coefficients ("place at infinity" of
order k); every residue product is computed with a short bilinear algorithm of the ring F2[x]/(h) (stored table
data/crt_ranks.json, flip-graph searched in the lab) and the interpolation is folded, together with the reduction
modulo f, into the output map.  The places are chosen by a knapsack DP over the table (`plan`).
A bilinear algorithm is a list of (alpha, beta, gamma) int masks: c = sum_k (alpha_k . a)(beta_k . b) gamma_k.

Ranks (= logical ANDs of the circuit) of the stored plans: m=163: 908, 233: 1343, 283: 1669, 409: 2501, 571: 3569.
The table covers the places needed by the optimal plan for every 2 <= m <= 600; for larger m the DP still works on
the stored places (result possibly not optimal).
"""
import json
import os

from .field import pmod, pmul, is_irreducible

_HERE = os.path.dirname(os.path.abspath(__file__))
_RK = None
_CHECKED = set()

# NIST / SEC binary-field moduli and a few low-weight moduli used in the literature tables
NIST = {163: (163, 7, 6, 3, 0), 233: (233, 74, 0), 283: (283, 12, 7, 5, 0), 409: (409, 87, 0),
        571: (571, 10, 5, 2, 0)}


def poly(*exps):
    f = 0
    for e in exps:
        f |= 1 << e
    return f


def default_modulus(m):
    """NIST modulus for m in (163, 233, 283, 409, 571); otherwise the lowest-weight irreducible
    (first trinomial x^m + x^k + 1, else first pentanomial, in lexicographic order of the middle exponents)."""
    if m in NIST:
        return poly(*NIST[m])
    if m == 1:
        return 0b11
    for k in range(1, m):
        f = poly(m, k, 0)
        if is_irreducible(f):
            return f
    for k3 in range(3, m):
        for k2 in range(2, k3):
            for k1 in range(1, k2):
                f = poly(m, k3, k2, k1, 0)
                if is_irreducible(f):
                    return f
    raise ValueError("no irreducible trinomial/pentanomial of degree %d" % m)


def ranks():
    global _RK
    if _RK is None:
        with open(os.path.join(_HERE, "data", "crt_ranks.json")) as fh:
            d = json.load(fh)
        _RK = {int(k): v for k, v in d["ranks"].items()}
    return _RK


def ring_dec(h):
    """verified bilinear algorithm of F2[x]/(h) (h a key of the table)."""
    v = ranks()[h]
    dec = [tuple(t) for t in v["dec"]]
    if h not in _CHECKED:
        k = h.bit_length() - 1
        for i in range(k):
            for j in range(k):
                c = 0
                for al, be, ga in dec:
                    if al >> i & 1 and be >> j & 1:
                        c ^= ga
                if c != pmod(1 << (i + j), h):
                    raise ValueError("corrupt ring algorithm for h=%d" % h)
        _CHECKED.add(h)
    return dec


def _groups():
    byp = {}
    for h, v in ranks().items():
        byp.setdefault(v["p"], []).append((v["deg"], v["rank"], ("mod", h)))
    G = [sorted(o) for p, o in sorted(byp.items())]
    xk = {v["e"]: v["rank"] for h, v in ranks().items() if v["p"] == 2}
    G.append([(k, r, ("inf", k)) for k, r in sorted(xk.items())])
    return G


def plan(n):
    """-> (rank, places) for the product of two polynomials with n coefficients (2n-1 product coefficients).
    places: tuple of ('mod', h) / ('inf', k); at most one place per irreducible p (one power), one at infinity."""
    N = 2 * n - 1
    INF = 10 ** 9
    dp = [(INF, ())] * (N + 1)
    dp[0] = (0, ())
    for g in _groups():
        nd = list(dp)
        for s in range(N + 1):
            c0, pl = dp[s]
            if c0 >= INF:
                continue
            for deg, r, place in g:
                t = s + deg
                if t <= N and c0 + r < nd[t][0]:
                    nd[t] = (c0 + r, pl + (place,))
        dp = nd
    if dp[N][0] >= INF:
        raise ValueError("no CRT plan for n=%d" % n)
    return dp[N]


def _inv_rows(rows, N):
    A = [(rows[i], 1 << i) for i in range(N)]
    for c in range(N):
        p = next((r for r in range(c, N) if A[r][0] >> c & 1), None)
        if p is None:
            raise ValueError("singular CRT system")
        A[c], A[p] = A[p], A[c]
        pc = A[c]
        for r in range(N):
            if r != c and A[r][0] >> c & 1:
                A[r] = (A[r][0] ^ pc[0], A[r][1] ^ pc[1])
    return [A[i][1] for i in range(N)]


def _ring_key(place):
    kind, v = place
    return v if kind == "mod" else 1 << v


def place_degree(place):
    kind, v = place
    return v.bit_length() - 1 if kind == "mod" else v


def powers_mod(v, n):
    """[x^i mod v for i < n] (incrementally)."""
    d = v.bit_length() - 1
    out = []
    r = 1 if d > 0 else 0
    for _ in range(n):
        out.append(r)
        r <<= 1
        if r >> d & 1:
            r ^= v
    return out


def place_rows(m, place):
    """linear forms (masks over the m input coefficients) of the residue coordinates of a place."""
    kind, v = place
    if kind == "mod":
        D = v.bit_length() - 1
        red = powers_mod(v, m)
        return [sum(1 << i for i in range(m) if red[i] >> j & 1) for j in range(D)]
    return [1 << (m - 1 - j) if m - 1 - j >= 0 else 0 for j in range(v)]


def slot_columns(m, f, places):
    """{(place index, j): mask over the m output bits} -- interpolation column of every CRT slot, reduced mod f."""
    N = 2 * m - 1
    rows, slots = [], []
    for pi, (kind, v) in enumerate(places):
        if kind == "mod":
            k = v.bit_length() - 1
            red = powers_mod(v, N)
            for j in range(k):
                rows.append(sum(1 << i for i in range(N) if red[i] >> j & 1))
                slots.append((pi, j))
        else:
            for j in range(v):
                rows.append(1 << (N - 1 - j))
                slots.append((pi, j))
    if len(rows) != N:
        raise ValueError("places do not have total degree 2m-1")
    Minv = _inv_rows(rows, N)
    red_f = powers_mod(f, N) if f else [1 << i for i in range(N)]
    col = {s: 0 for s in slots}
    for i in range(N):
        x = Minv[i]
        while x:
            z = (x & -x).bit_length() - 1
            x &= x - 1
            col[slots[z]] ^= red_f[i]
    return col


def _comb(base, mask):
    r = 0
    j = 0
    while mask:
        if mask & 1:
            r ^= base[j]
        mask >>= 1
        j += 1
    return r


def place_terms(m, place):
    """bilinear terms of one place: list of (alpha over a, beta over b, gamma over the place's slots)."""
    base = place_rows(m, place)
    out = []
    for al, be, ga in ring_dec(_ring_key(place)):
        out.append((_comb(base, al), _comb(base, be), ga))
    return out


def field_algorithm(m, f, places=None):
    """-> (dec, places): bilinear algorithm of c = a*b mod f (polynomial basis) from the CRT plan."""
    if places is None:
        _, places = plan(m)
    cols = slot_columns(m, f, places)
    dec = []
    for pi, pl in enumerate(places):
        D = place_degree(pl)
        for al, be, ga in place_terms(m, pl):
            g = 0
            for j in range(D):
                if ga >> j & 1:
                    g ^= cols[(pi, j)]
            if al and be and g:
                dec.append((al, be, g))
    return dec, tuple(places)


def check_algorithm(dec, m, f):
    """exact check of a bilinear algorithm against c = a*b mod f (all m^2 coefficient pairs)."""
    acc = {}
    for al, be, ga in dec:
        x = al
        while x:
            i = (x & -x).bit_length() - 1
            x &= x - 1
            y = be
            while y:
                j = (y & -y).bit_length() - 1
                y &= y - 1
                acc[(i, j)] = acc.get((i, j), 0) ^ ga
    for i in range(m):
        for j in range(m):
            if acc.get((i, j), 0) != pmod(1 << (i + j), f):
                return False
    return True


def check_random(dec, m, f, trials=64, seed=1):
    """randomised check of a bilinear algorithm (edge cases + `trials` random pairs).  Used for large m, where the
    exact tensor check costs ~ rank * m^2; the CRT construction itself is exact (verified ring algorithms,
    interpolation by an inverted matrix)."""
    import random
    rng = random.Random(seed)
    cases = [((1 << m) - 1, (1 << m) - 1), (1, (1 << m) - 1), (1 << (m - 1), 1 << (m - 1))]
    cases += [(rng.getrandbits(m), rng.getrandbits(m)) for _ in range(trials)]
    for a, b in cases:
        c = 0
        for al, be, ga in dec:
            if bin(al & a).count("1") & 1 and bin(be & b).count("1") & 1:
                c ^= ga
        if c != pmod(pmul(a, b), f):
            return False
    return True


def mulmod(a, b, f):
    return pmod(pmul(a, b), f)
