# SPDX-License-Identifier: Apache-2.0
"""Window optimisation of a CX path: a window of consecutive CX on <= wmax wires is replaced by the shortest (IDA*)
sequence with the same output rows that visits every mandatory (odd) term not visited outside the window.
Clifford weights / CZ are recomputed afterwards (big.solve_clifford)."""
import random, time
from . import big


class Limit(Exception): pass


def popc(x): return bin(x).count('1')


def exact_local(n, M, A, ub, node_limit=60000, timeout=None, rng=None):
    """Shortest CX sequence (< ub) from I to the rows A (tuple) visiting the masks M. Returns list | None | 'limit'."""
    rows0 = tuple(1 << i for i in range(n))
    Ml = sorted(set(M) - set(rows0))
    idx = {m: k for k, m in enumerate(Ml)}
    R0 = (1 << len(Ml)) - 1
    A = tuple(A)
    moves = [(c, t) for c in range(n) for t in range(n) if c != t]
    if rng: rng.shuffle(moves)
    nodes = [0]; t0 = time.time()

    def h(rows, R):
        pr = popc(R); ex = 0
        for i in range(n):
            if rows[i] != A[i]:
                ai = A[i]
                if ai not in idx or not (R >> idx[ai] & 1): ex += 1
        hv = pr + ex
        if pr and ex == 0:
            rs = set(rows); reach = False
            for k in range(len(Ml)):
                if R >> k & 1:
                    p = Ml[k]
                    for ri in rows:
                        if (p ^ ri) in rs: reach = True; break
                    if reach: break
            if not reach: hv += 1
        return hv

    def search(bound):
        seen = {}; path = []

        def dfs(rows, R, g):
            if g + h(rows, R) > bound: return False
            if R == 0 and rows == A: return True
            key = (rows, R)
            sg = seen.get(key)
            if sg is not None and sg <= g: return False
            seen[key] = g
            nodes[0] += 1
            if nodes[0] > node_limit or (timeout and time.time() - t0 > timeout): raise Limit()
            cand = []
            for c, t in moves:
                nr = rows[t] ^ rows[c]
                if nr == 0: continue
                R2 = R
                if nr in idx: R2 &= ~(1 << idx[nr])
                rows2 = rows[:t] + (nr,) + rows[t + 1:]
                f2 = g + 1 + h(rows2, R2)
                if f2 <= bound: cand.append((f2, -(popc(R) - popc(R2)), c, t, rows2, R2))
            cand.sort(key=lambda x: (x[0], x[1]))
            for f2, _, c, t, rows2, R2 in cand:
                path.append((c, t))
                if dfs(rows2, R2, g + 1): return True
                path.pop()
            return False
        return path if dfs(rows0, R0, 0) else None

    lb = h(rows0, R0)
    try:
        for bound in range(lb, ub):
            r = search(bound)
            if r is not None: return r
    except Limit:
        return 'limit'
    return None


def convex_window(n, ops, a, wmax, Lmax):
    """Convex sub-range of operations (commuting forward to position a). Returns (indices, wires, snapshots)."""
    W = []; blocked = set(); S = []
    loc = {}
    lrows = {}
    snap = {}
    i = a
    while i < len(ops) and len(S) < Lmax:
        c, t = ops[i]
        take = False
        if c not in blocked and t not in blocked:
            nw = set(W) | {c, t}
            if len(nw) <= wmax:
                take = True
                for w in (c, t):
                    if w not in loc:
                        loc[w] = len(W); W.append(w)
                        lrows[w] = None
        if take:
            S.append(i)
        else:
            for w in (c, t):
                if w in loc and w not in blocked:
                    snap[w] = ('blk', i)
                blocked.add(w)
        i += 1
    return S, W, snap


def window_pass(n, ops, O, terms, rng, wmax=7, Lmax=14, tries=200, node_limit=40000, verbose=False, cz_iters=1500, deadline=None, plateau=True):
    """Returns (ops, total_cost, sol)."""
    Oset = set(O)
    sol = big.solve_clifford(n, terms, ops, iters=cz_iters)
    cost = len(ops) + sol['ncz']
    for it in range(tries):
        if deadline and time.time() > deadline: break
        if len(ops) < 3: break
        rows = [1 << i for i in range(n)]
        created = []
        for c, t in ops:
            rows[t] ^= rows[c]; created.append(rows[t])
        occ = {}
        for i in range(n): occ[1 << i] = 1
        for m in created: occ[m] = occ.get(m, 0) + 1
        a = rng.randrange(len(ops))
        S, W, snap = convex_window(n, ops, a, wmax, Lmax)
        if len(S) < 3: continue
        rows = [1 << i for i in range(n)]
        for c, t in ops[:a]: rows[t] ^= rows[c]
        nl = len(W); loc = {w: k for k, w in enumerate(W)}
        lrows = [1 << k for k in range(nl)]
        grow = {w: rows[w] for w in W}
        Sset = set(S)
        mand = []
        own = {}
        for i in S: own[created[i]] = own.get(created[i], 0) + 1
        final = {}
        blk_at = {}
        for w, v in snap.items():
            blk_at.setdefault(v[1], []).append(w)
        for i in range(a, S[-1] + 1):
            if i in blk_at:
                for w in blk_at[i]: final[w] = lrows[loc[w]]
            if i in Sset:
                c, t = ops[i]
                lrows[loc[t]] ^= lrows[loc[c]]; grow[t] ^= grow[c]
                g = grow[t]
                if g in Oset and occ[g] - own[g] <= 0:
                    mand.append(lrows[loc[t]])
        for w in W:
            if w not in final: final[w] = lrows[loc[w]]
        A_loc = tuple(final[w] for w in W)
        res = exact_local(nl, mand, A_loc, len(S) + (1 if plateau else 0), node_limit=node_limit, rng=rng if plateau else None)
        if isinstance(res, list):
            rest = [ops[i] for i in range(a, len(ops)) if i not in Sset]
            new_ops = ops[:a] + [(W[c], W[t]) for c, t in res] + rest
            if len(new_ops) > len(ops) or (len(new_ops) == len(ops) and not plateau): continue
            if new_ops == ops: continue
            sol2 = big.solve_clifford(n, terms, new_ops, iters=cz_iters)
            if sol2 is None: continue
            c2 = len(new_ops) + sol2['ncz']
            if c2 < cost or (plateau and c2 == cost and rng.random() < 0.3):
                if verbose and c2 < cost: print('  convex a=%d |S|=%d wires=%d: %d -> %d (cx %d, cz %d)' % (a, len(S), nl, cost, c2, len(new_ops), sol2['ncz']), flush=True)
                ops, cost, sol = new_ops, c2, sol2
    return ops, cost, sol
