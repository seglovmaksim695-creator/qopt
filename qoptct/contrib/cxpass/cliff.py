# SPDX-License-Identifier: Apache-2.0
"""Freedom of the even (Clifford) phase weights. Odd terms O keep their weights; everything else is free.
For a set V of visited parities solve over Z4:  w'_p = w0_p + 2 u_p  (p in V; w0 = 0 outside O).
 (1) sum_{p in V, p contains i} u_p = a_i (mod 4);  (2) sum_{p in V, p contains ij} 2 u_p = 2 b_ij (mod 4).
Exact formula: coefficient of monomial x_S (|S| <= 3) is C_S = (-2)^{|S|-1} * sum_{p contains S} w'_p (mod 8)."""
import itertools


def coeffs(terms):
    C = {}
    for m, w in terms.items():
        bits = [i for i in range(m.bit_length()) if m >> i & 1]
        for r in (1, 2, 3):
            f = w * (-2) ** (r - 1)
            for S in itertools.combinations(bits, r):
                C[S] = (C.get(S, 0) + f) % 8
    return {k: v for k, v in C.items() if v}


def residual(n, terms):
    """-> (O: dict mask->w0, a: list[n] mod 4, b: dict (i,j)->bit) or None if inconsistent."""
    O = {m: w for m, w in terms.items() if w & 1}
    C = coeffs(terms)
    CO = coeffs(O)
    a = [0] * n; b = {}
    for S in set(C) | set(CO):
        d = (C.get(S, 0) - CO.get(S, 0)) % 8
        if len(S) == 1:
            if d & 1: return None
            a[S[0]] = (d // 2) % 4
        elif len(S) == 2:
            if d % 4: return None
            # (-2)*2*u_sum = d  -> u_sum = -d/4 mod 2
            b[S] = ((-d) // 4) % 2
        else:
            if d: return None
    return O, a, b


def solve_z4(rows, rhs, ncols):
    """rows: list of lists mod 4 (length ncols); rhs: list. Returns a solution (list mod 4) or None."""
    M = [r[:] + [h] for r, h in zip(rows, rhs)]
    nr = len(M)
    used = [False] * nr
    pivs = []  # (row, col)
    for c in range(ncols):
        pr = None
        for r in range(nr):
            if not used[r] and M[r][c] & 1:
                pr = r; break
        if pr is None: continue
        used[pr] = True
        if M[pr][c] == 3:
            M[pr] = [(3 * v) % 4 for v in M[pr]]
        for r in range(nr):
            if r != pr and M[r][c]:
                f = M[r][c]
                M[r] = [(v - f * u) % 4 for v, u in zip(M[r], M[pr])]
        pivs.append((pr, c))
    pcols = {c for _, c in pivs}
    free = [c for c in range(ncols) if c not in pcols]
    left = []
    for r in range(nr):
        if used[r]: continue
        if M[r][ncols] & 1: return None
        bits = 0
        for k, c in enumerate(free):
            v = M[r][c]
            if v & 1: return None
            if v == 2: bits |= 1 << k
        left.append((bits, (M[r][ncols] // 2) & 1))
    piv = {}
    for bits, rh in left:
        for p, (pb, pr_) in piv.items():
            if bits >> p & 1: bits ^= pb; rh ^= pr_
        if bits == 0:
            if rh: return None
            continue
        p = bits.bit_length() - 1
        for q in list(piv):
            if piv[q][0] >> p & 1: piv[q] = (piv[q][0] ^ bits, piv[q][1] ^ rh)
        piv[p] = (bits, rh)
    fv = [0] * len(free)
    for p, (pb, rh) in piv.items():
        fv[p] = rh
    x = [0] * ncols
    for k, c in enumerate(free): x[c] = fv[k]
    for pr, c in pivs:
        s = M[pr][ncols]
        for k, cc in enumerate(free):
            s -= M[pr][cc] * x[cc]
        x[c] = s % 4
    return x


def solve_V(n, V, a, b):
    """V: list of masks. Returns u (list mod 4 over V) or None."""
    rows = []; rhs = []
    for i in range(n):
        rows.append([(p >> i) & 1 for p in V]); rhs.append(a[i] % 4)
    for i in range(n):
        for j in range(i + 1, n):
            r = [2 * ((p >> i) & (p >> j) & 1) for p in V]
            bj = b.get((i, j), 0)
            if any(r) or bj:
                rows.append(r); rhs.append(2 * bj)
    return solve_z4(rows, rhs, len(V))
