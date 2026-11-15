# SPDX-License-Identifier: Apache-2.0
"""Linear CX synthesis helpers: Gaussian elimination and a greedy path search in row coordinates
(in the spirit of Paar / de Brugiere et al. and GraySynth). A term p has coordinates x = S^-1 p; CX(c->t) acts as x_c ^= x_t.
A term is visited when its weight is 1; the final rows A are tracked in coordinates as well (target = identity)."""
import random


def lin_cost_gauss(rows, target):
    """CX operations (c, t) taking the invertible row set `rows` to `target` by plain Gaussian elimination
    (rows -> I -> target)."""
    n = len(rows)

    def reduce(M):
        M = list(M); ops = []
        for col in range(n):
            if not (M[col] >> col & 1):
                for r in range(col + 1, n):
                    if M[r] >> col & 1:
                        M[col] ^= M[r]; ops.append((r, col)); break
                else:
                    raise ValueError('singular')
            for r in range(n):
                if r != col and M[r] >> col & 1:
                    M[r] ^= M[col]; ops.append((col, r))
        return ops
    o1 = reduce(rows)
    o2 = reduce(target)
    return o1 + o2[::-1]


def greedy_path(n, tmasks, A, lam=1, bonus=2, noise=0.0, rng=None, start=None, maxsteps=100000):
    """tmasks: mandatory masks (odd terms); A: target rows. Returns the full ops list [(c, t)] including the final linear part."""
    rng = rng or random.Random()
    S = list(start) if start else [1 << i for i in range(n)]
    nt = len(tmasks)
    assert start is None
    cols = [0] * n
    for k, m in enumerate(tmasks):
        for i in range(n):
            if m >> i & 1: cols[i] |= 1 << k
    for k, m in enumerate(A):
        for i in range(n):
            if m >> i & 1: cols[i] |= 1 << (nt + k)
    Tm = (1 << nt) - 1
    Am = ((1 << n) - 1) << nt
    wt = [0] * (nt + n)
    for i in range(n):
        c = cols[i]
        while c:
            b = c & -c; wt[b.bit_length() - 1] += 1; c ^= b
    active = Tm
    for k in range(nt):
        if wt[k] == 1: active &= ~(1 << k)
    ops = []
    def A_done():
        for k in range(n):
            if wt[nt + k] != 1: return False
        for k in range(n):
            if not (cols[k] >> (nt + k) & 1): return False
        return True
    steps = 0
    while steps < maxsteps:
        if active == 0 and A_done(): break
        steps += 1
        popT = [(c & active).bit_count() for c in cols]
        popA = [(c & Am).bit_count() for c in cols]
        best = None; bl = []
        for j in range(n):
            cj = cols[j]
            if cj & (active | Am) == 0: continue
            tj = popT[j]; aj = popA[j]
            for i in range(n):
                if i == j: continue
                ci = cols[i]
                gT = 2 * (ci & cj & active).bit_count() - tj
                gA = 2 * (ci & cj & Am).bit_count() - aj
                g = gT + lam * gA
                if active == 0: g = gA
                if best is None or g > best:
                    best = g; bl = [(i, j)]
                elif g == best or (noise and g >= best - noise and g > 0 and rng.random() < 0.2):
                    bl.append((i, j))
        if best is None or best <= 0:
            # fallback
            if active:
                kk = None; kw = 10 ** 9
                a = active
                while a:
                    b = a & -a; k = b.bit_length() - 1; a ^= b
                    if wt[k] < kw: kw = wt[k]; kk = k
                sup = [i for i in range(n) if cols[i] >> kk & 1]
                piv = sup[0]
                for i in sup[1:]:
                    _apply((i, piv), cols, wt, nt, n)
                    ops.append((i, piv))
                    S[piv] ^= S[i]
                active &= ~(1 << kk)
                for k in range(nt):
                    if active >> k & 1 and wt[k] == 1: active &= ~(1 << k)
                continue
            else:
                tail = lin_cost_gauss(S, A)
                ops.extend(tail)
                for c, t in tail: S[t] ^= S[c]
                break
        c, t = rng.choice(bl)
        _apply((c, t), cols, wt, nt, n)
        ops.append((c, t)); S[t] ^= S[c]
        a = active
        while a:
            b = a & -a; k = b.bit_length() - 1; a ^= b
            if wt[k] == 1: active &= ~b
    else:
        return None
    return ops


def _apply(ct, cols, wt, nt, n):
    c, t = ct
    d = cols[t]
    cc = cols[c]
    x = d
    while x:
        b = x & -x; k = b.bit_length() - 1; x ^= b
        if cc >> k & 1: wt[k] -= 1
        else: wt[k] += 1
    cols[c] = cc ^ d
