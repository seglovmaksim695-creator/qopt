# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""linsynth: CNOT synthesis of invertible GF(2) maps on a register (rows = int masks)  [qoptct v0.5, from lab k5].
op (c, t): rows[t] ^= rows[c].
  inverse(rows)                      -> inverse matrix (rows)
  matmul(A, B)                       -> A*B (rows)
  to_exact(cur, target)              -> ops with apply(cur) == target
  to_perm(cur, targets)              -> ops, perm ; apply(cur)[perm[i]] == targets[i]   (wire assignment free)
Methods: Gauss-Jordan with greedy pivots (random restarts), PMH (Patel-Markov-Hayes) sections 1..6, both on the
matrix and on its transpose (column ops = reversed transposed row ops); the shortest is returned."""
import random


def apply(rows, ops):
    r = list(rows)
    for c, t in ops: r[t] ^= r[c]
    return r


def inverse(rows):
    n = len(rows)
    A = [[rows[i], 1 << i] for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][0] >> c & 1), None)
        if p is None: raise ValueError('singular')
        A[c], A[p] = A[p], A[c]
        a0, a1 = A[c]
        for r in range(n):
            if r != c and A[r][0] >> c & 1:
                A[r][0] ^= a0; A[r][1] ^= a1
    return [A[i][1] for i in range(n)]


def matmul(A, B):
    out = []
    for a in A:
        r = 0; x = a
        while x:
            j = (x & -x).bit_length() - 1; x &= x - 1
            r ^= B[j]
        out.append(r)
    return out


def transpose(rows, n):
    out = [0] * n
    for i, r in enumerate(rows):
        x = r
        while x:
            j = (x & -x).bit_length() - 1; x &= x - 1
            out[j] |= 1 << i
    return out


def popc(x): return bin(x).count('1')


# ---------- reductions of X to the identity / a permutation (ops in application order) --------------------

def gauss_to_perm(X, rng=None, greedy=True):
    """row ops reducing X to a permutation matrix; pivot for column col = a not-yet-used row with that bit,
    chosen with minimum weight (ties random).  Returns ops, pos (pos[col] = row holding e_col)."""
    n = len(X); X = list(X); ops = []; used = [False] * n; pos = [0] * n
    cols = list(range(n))
    if rng: rng.shuffle(cols)
    for col in cols:
        cands = [r for r in range(n) if not used[r] and X[r] >> col & 1]
        if greedy:
            r0 = min(cands, key=lambda r: (popc(X[r]), rng.random() if rng else 0))
        else:
            r0 = cands[0]
        used[r0] = True; pos[col] = r0
        v = X[r0]
        for r in range(n):
            if r != r0 and X[r] >> col & 1:
                X[r] ^= v; ops.append((r0, r))
    return ops, pos


def gauss_to_id(X, rng=None):
    ops, pos = gauss_to_perm(X, rng)
    # finish the permutation with CX swaps (3 CX each) along cycles
    n = len(X); cur = apply(X, ops)
    where = {cur[r].bit_length() - 1: r for r in range(n)}  # column -> row
    for col in range(n):
        r = where[col]
        if r != col:
            # swap rows r and col
            ops += [(r, col), (col, r), (r, col)]
            cc = cur[col].bit_length() - 1
            cur[r], cur[col] = cur[col], cur[r]
            where[cc] = r; where[col] = col
    return ops


def _lower(rows, n, ss):
    rows = list(rows); ops = []
    for sec in range((n + ss - 1) // ss):
        lo, hi = sec * ss, min((sec + 1) * ss, n); msk = ((1 << (hi - lo)) - 1) << lo
        patt = {}
        for row in range(lo, n):
            sub = rows[row] & msk
            if sub == 0: continue
            if sub not in patt: patt[sub] = row
            else:
                rows[row] ^= rows[patt[sub]]; ops.append((patt[sub], row))
        for col in range(lo, hi):
            diag = (rows[col] >> col) & 1
            for row in range(col + 1, n):
                if (rows[row] >> col) & 1:
                    if not diag:
                        rows[col] ^= rows[row]; ops.append((row, col)); diag = 1
                    rows[row] ^= rows[col]; ops.append((col, row))
    return rows, ops


def pmh_to_id(X, ss):
    n = len(X)
    u1, p1 = _lower(X, n, ss)
    t = transpose(u1, n)
    _, p2 = _lower(t, n, ss)
    return list(p1) + [(tt, cc) for (cc, tt) in reversed(p2)]


def greedy_reduce(X):
    """greedy weight reduction: repeatedly apply the row op (i -> j) with the largest decrease of the total
    weight (gain = 2|X_i & X_j| - |X_i|); returns (ops, remaining matrix as int rows)."""
    import numpy as np
    n = len(X)
    M = np.zeros((n, n), np.int32)
    for i, r in enumerate(X):
        x = r
        while x:
            j = (x & -x).bit_length() - 1; x &= x - 1
            M[i, j] = 1
    O = M @ M.T
    w = np.diag(O).copy()
    ops = []
    while True:
        G = 2 * O - w[:, None]          # G[i, j] gain of rows[j] ^= rows[i]
        np.fill_diagonal(G, -10 ** 6)
        k = int(np.argmax(G)); i, j = divmod(k, n)
        if G[i, j] <= 0: break
        M[j] ^= M[i]
        ops.append((i, j))
        newrow = M @ M[j]
        O[j, :] = newrow; O[:, j] = newrow; w[j] = newrow[j]
    rows = [int(''.join('1' if b else '0' for b in M[i][::-1]), 2) for i in range(n)]
    return ops, rows


def reduce_best(X, perm, tries=4, seed=0, n_pmh=(2, 3, 4, 5, 6)):
    """shortest found op sequence reducing X to identity (perm=False) or to a permutation (perm=True)."""
    n = len(X); rng = random.Random(seed)
    ident = [1 << i for i in range(n)]
    cands = []
    # greedy then finish with Gauss/PMH
    try:
        g_ops, rest = greedy_reduce(X)
        fin = []
        for ss in (1, 2, 3, 4):
            fin.append(pmh_to_id(rest, ss) if not perm else None)
        if perm:
            fin = [gauss_to_perm(rest, None)[0]] + [gauss_to_perm(rest, random.Random(s))[0] for s in range(tries)]
        else:
            fin += [gauss_to_id(rest, None)]
        fin = [f_ for f_ in fin if f_ is not None]
        cands += [g_ops + f_ for f_ in fin]
    except Exception:
        pass
    for ss in n_pmh:
        if ss <= n:
            try: cands.append(pmh_to_id(X, ss))
            except Exception: pass
    if perm:
        for t in range(tries):
            cands.append(gauss_to_perm(X, rng if t else None)[0])
    else:
        for t in range(tries):
            cands.append(gauss_to_id(X, rng if t else None))
    best = None
    for c in cands:
        r = apply(X, c)
        ok = (r == ident) if not perm else all(v and v & (v - 1) == 0 for v in r) and len(set(r)) == n
        if ok and (best is None or len(c) < len(best)): best = c
    return best


def to_exact(cur, target, tries=4, seed=0):
    """ops O with O applied to state cur gives target:  reduce X = (target cur^-1)^-1 = cur target^-1 to I."""
    T = matmul(target, inverse(cur))          # target = T cur
    best = None
    for X, transposed in ((inverse(T), False), (transpose(T, len(T)), True)):
        ops = reduce_best(X, False, tries, seed)
        if transposed:
            # row ops reducing T^t to I  <=>  T^t = E1..Ek  ->  T = Ek^t..E1^t ; E^t of (c,t) is (t,c)
            ops = [(t, c) for (c, t) in ops]
        if apply(cur, ops) != list(target): continue
        if best is None or len(ops) < len(best): best = ops
    # reverse variants: reduce T (or T^-t) itself and apply the ops in reverse order
    n = len(T)
    for X, transposed in ((T, False), (transpose(inverse(T), n), True)):
        ops = reduce_best(X, False, tries, seed)
        if transposed: ops = [(t, c) for (c, t) in ops]
        ops = ops[::-1]
        if apply(cur, ops) != list(target): continue
        if best is None or len(ops) < len(best): best = ops
    assert best is not None
    return best


def to_perm(cur, targets, tries=6, seed=0):
    """ops and perm: apply(cur, ops)[perm[i]] == targets[i]."""
    n = len(cur)
    T = matmul(targets, inverse(cur))
    best = None
    for X, transposed in ((inverse(T), False), (transpose(T, n), True)):
        ops = reduce_best(X, True, tries, seed)
        if transposed: ops = [(t, c) for (c, t) in ops]
        res = apply(cur, ops)
        idx = {v: i for i, v in enumerate(res)}
        if not all(t in idx for t in targets): continue
        if best is None or len(ops) < len(best[0]): best = (ops, [idx[t] for t in targets])
    # reverse variant: E_k..E_1 T = Q (permutation)  ->  T = E_1..E_k Q ; virtual wire i = physical pi[i]
    try:
        ops = reduce_best(T, True, tries, seed)
        Q = apply(T, ops)
        pi = [Q[i].bit_length() - 1 for i in range(n)]
        ops2 = [(pi[c], pi[t]) for (c, t) in ops[::-1]]
        res = apply(cur, ops2)
        if all(res[pi[i]] == targets[i] for i in range(n)) and len(ops2) < len(best[0]):
            best = (ops2, list(pi))
    except Exception:
        pass
    assert best is not None
    return best
