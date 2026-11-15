# SPDX-License-Identifier: Apache-2.0
"""Exact IDA* with free even (Clifford) weights. Odd terms O are mandatory."""
import time
from . import cliff
from .blocks import phase_gates


def popc(x): return bin(x).count('1')


class Limit(Exception): pass


def exact2(n, p_terms, ub, node_limit=500000, timeout=None, A=None, lo=None):
    """p_terms: dict mask->w (block polynomial), A: target rows. Searches a path shorter than ub.
    Returns (ops, Vlist, u) | None | 'limit'."""
    res = cliff.residual(n, p_terms)
    if res is None: return None
    O, a, b = res
    Ol = sorted(O)
    idx = {m: k for k, m in enumerate(Ol)}
    rows0 = tuple(1 << i for i in range(n))
    R0 = (1 << len(Ol)) - 1
    for r in rows0:
        if r in idx: R0 &= ~(1 << idx[r])
    A = tuple(A)
    moves_all = [(c, t) for c in range(n) for t in range(n) if c != t]
    nodes = [0]
    t0 = time.time()
    solcache = {}

    Aset_idx = idx

    def h(rows, R):
        pr = popc(R)
        ex = 0
        for i in range(n):
            ri = rows[i]
            if ri != A[i]:
                ai = A[i]
                if ai not in idx or not (R >> idx[ai] & 1): ex += 1
        hv = pr + ex
        if pr and ex == 0:
            rs = set(rows)
            reach = False
            for k in range(len(Ol)):
                if R >> k & 1:
                    p = Ol[k]
                    for ri in rows:
                        if (p ^ ri) in rs: reach = True; break
                    if reach: break
            if not reach: hv += 1
        return hv

    def search(bound):
        seen = {}
        path = []

        def dfs(rows, R, V, g):
            hh = h(rows, R)
            if g + hh > bound: return None
            if R == 0 and rows == A:
                key = frozenset(V)
                if key not in solcache:
                    solcache[key] = cliff.solve_V(n, sorted(V), a, b)
                u = solcache[key]
                if u is not None: return (list(path), sorted(V), u)
            key = (rows, R, frozenset(V))
            sg = seen.get(key)
            if sg is not None and sg <= g: return None
            seen[key] = g
            nodes[0] += 1
            if nodes[0] > node_limit or (timeout and time.time() - t0 > timeout): raise Limit()
            cand = []
            for c, t in moves_all:
                nr = rows[t] ^ rows[c]
                if nr == 0: continue
                R2 = R
                if nr in idx: R2 &= ~(1 << idx[nr])
                rows2 = rows[:t] + (nr,) + rows[t + 1:]
                f2 = g + 1 + h(rows2, R2)
                if f2 <= bound:
                    cand.append((f2, -(popc(R) - popc(R2)), nr in V, c, t, rows2, R2, nr))
            cand.sort(key=lambda x: (x[0], x[1], x[2]))
            for f2, _, inV, c, t, rows2, R2, nr in cand:
                path.append((c, t))
                V2 = V if inV else (V | {nr})
                r = dfs(rows2, R2, V2, g + 1)
                if r is not None: return r
                path.pop()
            return None
        return dfs(rows0, R0, frozenset(rows0), 0)

    lb = h(rows0, R0)
    if lo: lb = max(lb, lo)
    try:
        for bound in range(lb, ub):
            r = search(bound)
            if r is not None: return r
    except Limit:
        return 'limit'
    return None


def build_gates2(n, p_terms, off, ops, Vs, u):
    O, a, b = cliff.residual(n, p_terms)
    uw = {p: x for p, x in zip(Vs, u)}
    rows = [1 << i for i in range(n)]
    done = set(); g = []

    def visit(w):
        m = rows[w]
        if m in done: return
        done.add(m)
        wt = (O.get(m, 0) + 2 * uw.get(m, 0)) % 8
        g.extend(phase_gates(wt, w))
    for w in range(n): visit(w)
    for c, t in ops:
        rows[t] ^= rows[c]
        g.append(('cx', [c, t])); visit(t)
    for w in range(n):
        if off[w]: g.append(('x', [w]))
    return g
