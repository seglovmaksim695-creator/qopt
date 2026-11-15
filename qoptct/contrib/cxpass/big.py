# SPDX-License-Identifier: Apache-2.0
"""Block synthesis for a fixed CX path: even (Clifford) weights and CZ gates are recomputed.
Odd terms O are mandatory; S/Z on any visited row p adds p p^T to a symmetric F2 matrix (off-diagonal part);
CZ(p,q) on two rows present at the same time adds p q^T + q p^T. The number of CZ is minimised by Gaussian reduction of
the residual, weight descent in the pair basis and a greedy cover with dense CZ pairs over all moments of the path.
The diagonal and the mod-4 carry are closed with rows e_i (always visited)."""
import random
from . import cliff
from .blocks import phase_gates


def path_rows(n, ops):
    rows = [1 << i for i in range(n)]
    vis = list(rows); seen = set(rows)
    for c, t in ops:
        rows[t] ^= rows[c]
        if rows[t] not in seen: seen.add(rows[t]); vis.append(rows[t])
    return vis, rows


def _bits(m):
    out = []; i = 0
    while m:
        if m & 1: out.append(i)
        m >>= 1; i += 1
    return out


def solve_clifford(n, terms, ops, iters=3000, seed=0, A=None, cands=True, maxcand_pairs=400000):
    """Returns dict(czs=[(i, j, k)], wt={mask: weight mod 8}, ncz) or None."""
    res = cliff.residual(n, terms)
    if res is None: return None
    O, a, b = res
    vis, final = path_rows(n, ops)
    vs = set(vis)
    if any(m not in vs for m in O): return None
    pos = {}
    k_ = 0
    for i in range(n):
        for j in range(i + 1, n):
            pos[(i, j)] = k_; k_ += 1

    def vec_pair(p, q):
        v = 0
        bp = _bits(p); bq = _bits(q)
        for x in bp:
            for y in bq:
                if x != y:
                    v ^= 1 << pos[(x, y) if x < y else (y, x)]
        return v

    def vec(pm):
        bits = _bits(pm); v = 0
        for x in range(len(bits)):
            for y in range(x + 1, len(bits)):
                v |= 1 << pos[(bits[x], bits[y])]
        return v
    gens = []; gm = []
    for m in vis:
        if m & (m - 1):
            gens.append(vec(m)); gm.append(m)
    r = 0
    for (i, j), x in b.items():
        if x & 1: r |= 1 << pos[(i, j)]
    rng = random.Random(seed)
    piv = {}
    for gi, g in enumerate(gens):
        u_ = 1 << gi
        while g:
            h = g.bit_length() - 1
            if h in piv:
                g ^= piv[h][0]; u_ ^= piv[h][1]
            else:
                piv[h] = (g, u_); break
    use = 0
    for h in sorted(piv, reverse=True):
        if r >> h & 1:
            r ^= piv[h][0]; use ^= piv[h][1]
    G = len(gens); order = list(range(G))

    def descend(r, use):
        improved = True
        while improved:
            improved = False
            rng.shuffle(order)
            for gi in order:
                nr = r ^ gens[gi]
                if nr.bit_count() < r.bit_count():
                    r = nr; use ^= 1 << gi; improved = True
        return r, use
    r, use = descend(r, use)
    chosen = []   # (vec, k, i, j, p, q)
    best = (r.bit_count(), r, use, list(chosen))
    if cands and r.bit_count() > 1:
        rows = [1 << i for i in range(n)]
        cl = []; seenv = set()
        for k, (c, t) in enumerate(ops):
            rows[t] ^= rows[c]
            p = rows[t]
            if p & (p - 1) == 0: continue
            for j in range(n):
                if j == t: continue
                q = rows[j]
                v = vec_pair(p, q)
                if v and v not in seenv:
                    seenv.add(v); cl.append((v, k + 1, t, j, p, q))
            if len(cl) > maxcand_pairs: break
        while True:
            w0 = r.bit_count()
            if w0 <= 1: break
            bw = w0 - 1; bc = None
            for cnd in cl:
                nw = (r ^ cnd[0]).bit_count()
                if nw < bw: bw = nw; bc = cnd
            if bc is None: break
            r ^= bc[0]; chosen.append(bc)
            r, use = descend(r, use)
            if r.bit_count() + len(chosen) < best[0] + len(best[3]):
                best = (r.bit_count(), r, use, list(chosen))
        r, use, chosen = best[1], best[2], best[3]
    czs = [(i, j, 0) for (i, j), q in pos.items() if r >> q & 1]
    cdiag = [0] * n
    for v, k, i, j, p, q in chosen:
        czs.append((i, j, k))
        for x in _bits(p & q): cdiag[x] ^= 1
    u0 = {}
    for gi in range(G):
        if use >> gi & 1: u0[gm[gi]] = 1
    for i in range(n):
        s = sum(1 for m in u0 if m >> i & 1 and m != (1 << i))
        if (s & 1) != (a[i] & 1): u0[1 << i] = 1
    u = {}
    for i in range(n):
        s = sum(1 for m in u0 if m >> i & 1)
        d = (a[i] - 2 * cdiag[i] - s) % 4
        assert d % 2 == 0, 'odd difference'
        e = 1 << i
        u[e] = (u0.get(e, 0) + 2 * (d // 2)) % 4
    for m, x in u0.items():
        if m & (m - 1): u[m] = x % 4
    wt = {}
    for m in vis:
        w = (O.get(m, 0) + 2 * u.get(m, 0)) % 8
        if w: wt[m] = w
    return dict(czs=czs, wt=wt, ncz=len(czs))


def build_big(n, off, ops, sol):
    """CZ at their moments, phases at the first visit of a row, CX, X gates at the end."""
    g = []
    cz_at = {}
    for i, j, k in sol['czs']: cz_at.setdefault(k, []).append((i, j))
    rows = [1 << i for i in range(n)]
    wt = dict(sol['wt']); done = set()

    def visit(w):
        m = rows[w]
        if m in done: return
        done.add(m)
        g.extend(phase_gates(wt.get(m, 0), w))
    for i, j in cz_at.get(0, []): g.append(('cz', [i, j]))
    for w in range(n): visit(w)
    for k, (c, t) in enumerate(ops):
        rows[t] ^= rows[c]
        g.append(('cx', [c, t])); visit(t)
        for i, j in cz_at.get(k + 1, []): g.append(('cz', [i, j]))
    for w in range(n):
        if off[w]: g.append(('x', [w]))
    return g
