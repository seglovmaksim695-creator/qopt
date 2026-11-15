# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""GF(2) linear algebra on int bitmask rows and CNOT-network synthesis.

Conventions: a register state is a list of row masks; a CNOT op (c, t) means rows[t] ^= rows[c].
"""
import random


def popc(x):
    return bin(x).count("1")


def express(alpha, rows):
    """alpha as XOR of a subset of `rows` -> set of indices, or None if alpha is not in the span."""
    basis = []  # (pivot bit, vector, index set)
    for i, r in enumerate(rows):
        s = {i}
        v = r
        for pb, bv, bs in basis:
            if v >> pb & 1:
                v ^= bv
                s ^= bs
        if v:
            pb = v.bit_length() - 1
            nb = []
            for qb, qv, qs in basis:
                if qv >> pb & 1:
                    qv ^= v
                    qs = qs ^ s
                nb.append((qb, qv, qs))
            basis = nb + [(pb, v, s)]
    v = alpha
    s = set()
    for pb, bv, bs in basis:
        if v >> pb & 1:
            v ^= bv
            s ^= bs
    return s if v == 0 else None


def rank(rows):
    piv = {}
    r = 0
    for v in rows:
        while v:
            b = v.bit_length() - 1
            if b in piv:
                v ^= piv[b]
            else:
                piv[b] = v
                r += 1
                break
    return r


def apply_cnots(rows, ops):
    r = list(rows)
    for c, t in ops:
        r[t] ^= r[c]
    return r


def gauss_synth(M):
    """CNOT circuit implementing x -> M x (rows of M = new value of wire i over old wires)."""
    n = len(M)
    A = list(M)
    ops = []
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r] >> col & 1:
                piv = r
                break
        if piv is None:
            raise ValueError("singular")
        if piv != col:
            A[col] ^= A[piv]
            ops.append((piv, col))
        for r in range(n):
            if r != col and A[r] >> col & 1:
                A[r] ^= A[col]
                ops.append((col, r))
    return list(reversed(ops))


def synth_cnot(rows, target, rng=None, tries=1):
    """Greedy CNOT circuit taking register `rows` to `target` (same span); None on failure."""
    best = None
    n = len(rows)
    for _ in range(tries):
        cur = list(rows)
        ops = []
        guard = 0
        while cur != target and guard < 10 * n * n:
            guard += 1
            cands = []
            for w in range(n):
                if cur[w] == target[w]:
                    continue
                S = express(target[w], cur)
                if S is None:
                    raise ValueError("not in span")
                cands.append((len(S - {w}) if w in S else 99, w, S))
            cands.sort(key=lambda x: x[0])
            ok = [c for c in cands if c[0] < 99]
            if ok:
                k = ok[0][0]
                pool = [c for c in ok if c[0] == k]
                _, w, S = pool[rng.randrange(len(pool))] if rng else pool[0]
                for i in sorted(S - {w}):
                    ops.append((i, w))
                    cur[w] ^= cur[i]
            else:
                w = cands[0][1]
                i = sorted(cands[0][2])[0]
                ops.append((w, i))
                cur[i] ^= cur[w]
        if cur != target:
            continue
        if best is None or len(ops) < len(best):
            best = ops
    return best


# ---- PMH-style reductions (Patel-Markov-Hayes) --------------------------------------------------------------

def _gauss_reduce(rows, n, rng=None):
    rows = list(rows)
    ops = []
    for col in range(n):
        if not (rows[col] >> col) & 1:
            cands = [r for r in range(col + 1, n) if (rows[r] >> col) & 1]
            r = rng.choice(cands) if rng else cands[0]
            rows[col] ^= rows[r]
            ops.append((r, col))
        for r in range(n):
            if r != col and (rows[r] >> col) & 1:
                rows[r] ^= rows[col]
                ops.append((col, r))
    return ops


def _lower(rows, n, ss):
    rows = list(rows)
    ops = []
    for sec in range((n + ss - 1) // ss):
        lo, hi = sec * ss, min((sec + 1) * ss, n)
        msk = ((1 << (hi - lo)) - 1) << lo
        patt = {}
        for row in range(lo, n):
            sub = rows[row] & msk
            if sub == 0:
                continue
            if sub not in patt:
                patt[sub] = row
            else:
                rows[row] ^= rows[patt[sub]]
                ops.append((patt[sub], row))
        for col in range(lo, hi):
            diag = (rows[col] >> col) & 1
            for row in range(col + 1, n):
                if (rows[row] >> col) & 1:
                    if not diag:
                        rows[col] ^= rows[row]
                        ops.append((row, col))
                        diag = 1
                    rows[row] ^= rows[col]
                    ops.append((col, row))
    return rows, ops


def _transpose(rows, n):
    return [sum(((rows[i] >> j) & 1) << i for i in range(n)) for j in range(n)]


def _pmh_reduce(rows, n, ss=2):
    u1, p1 = _lower(rows, n, ss)
    t = _transpose(u1, n)
    _, p2 = _lower(t, n, ss)
    return list(p1) + [(tt, cc) for (cc, tt) in reversed(p2)]


def reduce_ops(rows, n, fast=True, tries=30, seed=0):
    """CNOT ops taking register state `rows` (invertible) to the identity."""
    if fast:
        return _gauss_reduce(rows, n)
    rng = random.Random(seed)
    cands = [_gauss_reduce(rows, n)]
    for ss in (1, 2, 3, 4):
        try:
            cands.append(_pmh_reduce(rows, n, ss))
        except Exception:
            pass
    for _ in range(tries):
        try:
            cands.append(_gauss_reduce(rows, n, rng))
        except Exception:
            pass
    ident = [1 << i for i in range(n)]
    best = None
    for c in cands:
        if apply_cnots(rows, c) == ident and (best is None or len(c) < len(best)):
            best = c
    return best


def best_synth(M, rng, tries=30):
    """Shortest found CNOT circuit for x -> M x (identity -> rows M)."""
    n = len(M)
    best = gauss_synth(M)
    if n > 16:  # the greedy search is O(n^4) per try: PMH-style reductions only
        try:
            r = reduce_ops(list(M), n, fast=False, tries=2)
            if r is not None and len(r) < len(best):
                alt = list(reversed(r))
                if apply_cnots([1 << i for i in range(n)], alt) == list(M):
                    best = alt
        except Exception:
            pass
        return best
    try:
        g = synth_cnot([1 << i for i in range(n)], list(M), rng, tries)
        if g is not None and len(g) < len(best):
            best = g
    except Exception:
        pass
    return best


def relmat(cur, target):
    """matrix T (rows as masks over indices of cur) with target[i] = XOR_j T[i][j] cur[j]."""
    T = []
    for t in target:
        S = express(t, cur)
        if S is None:
            raise ValueError("not in span")
        T.append(sum(1 << j for j in S))
    return T


def transform(cur, target, rng=None, tries=20):
    """CNOT ops taking register state `cur` to `target`."""
    if list(cur) == list(target):
        return []
    T = relmat(cur, target)
    ops = best_synth(T, rng or random.Random(1), tries)
    if apply_cnots(cur, ops) != list(target):
        raise AssertionError("CNOT synthesis failed")
    return ops
