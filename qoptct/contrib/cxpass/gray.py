# SPDX-License-Identifier: Apache-2.0
"""GraySynth (Amy, Azimzadeh, Mosca 2018, arXiv:1712.01859) in coordinate form + a linear-synthesis tail to the target rows A.
Operation (c,t): S_t ^= S_c; on coordinates x_c ^= x_t."""
import random
from .linsynth import lin_cost_gauss, greedy_path


def inv_matrix(rows):
    n = len(rows)
    M = list(rows); R = [1 << i for i in range(n)]
    for col in range(n):
        p = None
        for r in range(col, n):
            if M[r] >> col & 1: p = r; break
        if p is None: raise ValueError('singular')
        M[col], M[p] = M[p], M[col]; R[col], R[p] = R[p], R[col]
        for r in range(n):
            if r != col and M[r] >> col & 1:
                M[r] ^= M[col]; R[r] ^= R[col]
    return R


def coords(S, targets):
    R = inv_matrix(S)
    out = []
    for a in targets:
        x = 0; i = 0
        while a:
            if a & 1: x ^= R[i]
            a >>= 1; i += 1
        out.append(x)
    return out


def tail_ops(S, A, rng=None, tries=3):
    """CX operations taking rows S to A: best of Gaussian elimination and the greedy path."""
    n = len(S)
    if list(S) == list(A): return []
    best = lin_cost_gauss(S, A)
    E = coords(S, A)
    rng = rng or random.Random(0)
    for lam in (1,):
        for r in range(tries):
            o = greedy_path(n, [], E, lam=lam, noise=0 if r == 0 else 1.0, rng=rng)
            if o is not None and len(o) < len(best): best = o
    return best


def gray_path(n, tmasks, A, rng, mode='max', tail=True, tailtries=3):
    X = list(tmasks)
    S = [1 << i for i in range(n)]
    ops = []

    def op(c, t):
        bt = 1 << t; bc = 1 << c
        for k in range(len(X)):
            if X[k] & bt: X[k] ^= bc
        S[t] ^= S[c]; ops.append((c, t))

    def rec(Sidx, I, j):
        Sidx = [k for k in Sidx if X[k] & (X[k] - 1)]
        if not Sidx: return
        if not I: raise RuntimeError('I empty')
        cnt = {}
        for l in I:
            bl = 1 << l
            cnt[l] = sum(1 for k in Sidx if X[k] & bl)
        cand = [l for l in I if cnt[l] > 0]
        if not cand: raise RuntimeError('no col')
        if mode == 'max':
            m = max(cnt[l] for l in cand); pick = [l for l in cand if cnt[l] == m]
        elif mode == 'bal':
            half = len(Sidx) / 2
            m = min(abs(cnt[l] - half) for l in cand); pick = [l for l in cand if abs(cnt[l] - half) == m]
        else:
            pick = cand
        l = rng.choice(pick)
        bl = 1 << l
        S0 = [k for k in Sidx if not X[k] & bl]
        S1 = [k for k in Sidx if X[k] & bl]
        I2 = [x for x in I if x != l]
        if j is None:
            rec(S0, I2, None)
            rec(S1, I2, l)
        else:
            rec(S0, I2, j)
            if S1:
                op(l, j)
                rec(S1, I2, j)

    rec(list(range(len(X))), list(range(n)), None)
    if tail:
        t = tail_ops(S, A, rng, tailtries)
        ops.extend(t)
    return ops
