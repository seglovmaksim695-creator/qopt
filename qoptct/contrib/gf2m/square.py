# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""square: in-place squaring  a <- a^2 mod f  in GF(2^m) as a CNOT network (0 T)  [qoptct v0.5, lab k5 square.py,
square2.py].

Squaring is GF(2)-linear (Frobenius), so the circuit is a CNOT synthesis of the m x m matrix S with columns
x^(2j) mod f.  Two variants:
  * exact wire order (output bit i on input wire i): best of Gauss/PMH/greedy syntheses (linsynth.to_exact);
  * with free wire relabelling (relabel=True): output bit i ends on wire out_wires[i].  For odd m the structured
    construction of k5/square2 is used: wire j is relabelled to position 2j (2j < m) or 2j - m; the map becomes
    I + N with N sparse (reduction terms) and is realised with nnz(N) CNOTs in topological order of the dependency
    graph (strongly connected parts synthesised separately).  For even m: linsynth.to_perm.
Every returned network is checked exactly (images of all unit vectors; the map is linear).
Numbers of this implementation (NIST moduli, seeds=3, tests/test_crypto_gf2m.py): exact wire order
798 / 825 / 1386 / 1418 / 3239 CX (m = 163/233/283/409/571; lab k5 files: 875 / 825 / 1534 / 1418 / 3239 --
these older files are the proved squarers kept in results/crypto; the greedy candidate needs numpy, without it
the exact-order counts are 888 / 825 / 1653 / 1418 / 3827);
with relabelling 408 / 153 / 713 / 204 / 1538 CX (the x5 audit re-ran k5/square2 alone: 408/153/713/247/1538;
for m = 409 the linsynth.to_perm candidate, 204, is shorter than the structured one).
"""
import collections

from . import linsynth as LS
from .field import pmod
from .qcirc import GF2Circuit


def square_rows(m, f):
    cols = [pmod(1 << (2 * j), f) for j in range(m)]
    return LS.transpose(cols, m)


def exact_ops(m, f, seeds=3):
    rows = square_rows(m, f)
    best = None
    for s in range(seeds):
        ops = LS.to_exact([1 << i for i in range(m)], rows, tries=4, seed=s)
        if best is None or len(ops) < len(best):
            best = ops
    return best


def _structured(m, f):
    sig = [2 * j if 2 * j < m else 2 * j - m for j in range(m)]
    img = [pmod(1 << (2 * j), f) for j in range(m)]
    T = [sum(1 << j for j in range(m) if img[j] >> sig[w] & 1) for w in range(m)]
    if not all(T[w] >> w & 1 for w in range(m)):
        return None
    N = [T[w] ^ (1 << w) for w in range(m)]
    adj = collections.defaultdict(list)
    for u in range(m):
        x = N[u]
        while x:
            v = (x & -x).bit_length() - 1
            x &= x - 1
            adj[u].append(v)
    # strongly connected components (iterative Tarjan); edges u -> v mean u must be done before v changes
    index, low, onst, st, comps = {}, {}, set(), [], []
    cnt = 0
    for root in range(m):
        if root in index:
            continue
        work = [(root, 0)]
        while work:
            u, k = work.pop()
            if k == 0:
                index[u] = low[u] = cnt
                cnt += 1
                st.append(u)
                onst.add(u)
            recurse = False
            for i in range(k, len(adj[u])):
                v = adj[u][i]
                if v not in index:
                    work.append((u, i + 1))
                    work.append((v, 0))
                    recurse = True
                    break
                elif v in onst:
                    low[u] = min(low[u], index[v])
            if recurse:
                continue
            if low[u] == index[u]:
                comp = []
                while True:
                    w = st.pop()
                    onst.discard(w)
                    comp.append(w)
                    if w == u:
                        break
                comps.append(comp)
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[u])
    comps.reverse()
    ops = []
    for comp in comps:
        S = sorted(comp)
        Sset = set(S)
        if len(S) > 1:
            loc = {w: i for i, w in enumerate(S)}
            Tl = [sum(1 << loc[v] for v in S if T[u] >> v & 1) for u in S]
            lops = LS.to_exact([1 << i for i in range(len(S))], Tl, 4, 0)
            ops += [(S[c], S[t]) for c, t in lops]
        for u in S:
            x = N[u]
            while x:
                v = (x & -x).bit_length() - 1
                x &= x - 1
                if v not in Sset:
                    ops.append((v, u))
    if LS.apply([1 << i for i in range(m)], ops) != T:
        return None
    out = [0] * m
    for w in range(m):
        out[sig[w]] = w
    return ops, out


def relabel_ops(m, f, seeds=3):
    """-> (ops, out_wires): output bit i of a^2 on wire out_wires[i]."""
    cands = []
    if m % 2 == 1:
        r = _structured(m, f)
        if r is not None:
            cands.append(r)
    rows = square_rows(m, f)
    for s in range(seeds if not cands else 1):
        ops, perm = LS.to_perm([1 << i for i in range(m)], rows, tries=6, seed=s)
        cands.append((ops, list(perm)))
    return min(cands, key=lambda c: len(c[0]))


def check_linear(ops, m, f, out_wires=None):
    """exact: the CNOT network maps every unit vector a = x^j to x^(2j) mod f (on out_wires)."""
    st = LS.apply([1 << i for i in range(m)], ops)
    ow = list(range(m)) if out_wires is None else out_wires
    for j in range(m):
        v = sum((st[ow[i]] >> j & 1) << i for i in range(m))
        if v != pmod(1 << (2 * j), f):
            return False
    return True


def build(m, f, relabel=False, seeds=3):
    if relabel:
        ops, ow = relabel_ops(m, f, seeds)
    else:
        ops, ow = exact_ops(m, f, seeds), list(range(m))
    if not check_linear(ops, m, f, ow):
        raise AssertionError("squaring network failed its exact check")
    C = GF2Circuit(m, [('c', c, t) for c, t in ops], regs=dict(a=list(range(m))))
    C.regs["out"] = list(ow)
    return C
