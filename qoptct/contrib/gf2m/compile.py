# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""Compilers of  y <- M y + Lin(x) + c + sum_k (alpha_k . a)(beta_k . b) gamma_k  into Clifford+T with logical ANDs.

Every product costs one logical AND (4 T, 3 CX) into an ancilla, a fan-out of CX into the targets and a
measurement-based erasure (0 T, one conditional CZ).  Two strategies (both from the qoptct lab):
  * ``frame``  (a2): operand registers tracked as linear states, output side either through a fixed basis
    (N1, N2 with N2 N1 = M) or by tracking the column frame backwards in time; clean model may AND directly
    into an untouched target.
  * ``walk``   (a6): operand registers walked with a dual-basis tracker, sparsified output matrix (row ops on the
    targets before/after the products), strict model applies M first.
"""
import math
import random

from . import gf2lin as lin
from .gf2lin import popc


class Problem:
    """wire-level description: A, B, T wire lists; M rows (masks over T indices); lin {t: [sources, const]};
    zero (targets start at 0, clean); anc = ancilla wire; nq = number of qubits."""

    def __init__(self, A, B, T, M, linpart, zero, anc, nq):
        self.A, self.B, self.T, self.M = list(A), list(B), list(T), list(M)
        self.lin, self.zero, self.anc, self.nq = linpart, zero, anc, nq


class Emit:
    def __init__(self, nq):
        self.L = []
        self.nb = 0
        self.nq = nq

    def cx(self, c, t):
        self.L.append(f"cx q[{c}], q[{t}];")

    def x(self, t):
        self.L.append(f"x q[{t}];")

    def AND(self, p, q, a):
        self.L += [f"h q[{a}];", f"t q[{a}];", f"cx q[{q}], q[{a}];", f"tdg q[{a}];", f"cx q[{p}], q[{a}];",
                   f"t q[{a}];", f"cx q[{q}], q[{a}];", f"tdg q[{a}];", f"h q[{a}];", f"sdg q[{a}];"]

    def erase(self, p, q, a):
        k = self.nb
        self.nb += 1
        self.L += [f"h q[{a}];", f"c[{k}] = measure q[{a}];", f"if (c[{k}]) {{cz q[{p}], q[{q}]; x q[{a}];}}"]

    def text(self):
        return "\n".join(["OPENQASM 3;", 'include "stdgates.inc";', f"qubit[{self.nq}] q;",
                          f"bit[{max(self.nb, 1)}] c;"] + self.L) + "\n"


def _finish_linear(E, pb):
    for t in pb.T:
        for s in sorted(pb.lin[t][0]):
            E.cx(s, t)
        if pb.lin[t][1]:
            E.x(t)


# ============================== strategy "frame" (a2) =========================================================

class _Reg:
    def __init__(self, wires):
        self.w = list(wires)
        self.cur = [1 << i for i in range(len(wires))]

    def get(self, alpha, E, rng, future=()):
        if alpha in self.cur:
            return self.w[self.cur.index(alpha)]
        S = lin.express(alpha, self.cur)
        fut = list(future)

        def score(i):
            c2 = list(self.cur)
            for j in S - {i}:
                c2[i] ^= c2[j]
            v = 0
            if fut and fut[0] not in c2:
                v = len(lin.express(fut[0], c2)) - 1
            return (v, self.cur[i] in set(fut[1:]), rng.random())
        piv = min(S, key=score)
        for i in sorted(S - {piv}):
            E.cx(self.w[i], self.w[piv])
            self.cur[piv] ^= self.cur[i]
        return self.w[piv]

    def restore(self, E, rng):
        ops = lin.transform(self.cur, [1 << i for i in range(len(self.w))], rng)
        for c, t in ops:
            E.cx(self.w[c], self.w[t])
        self.cur = [1 << i for i in range(len(self.w))]


def out_basis(M, G, rng, zero, iters=400):
    """choose N2 (products XORed between N1 and N2, N2 N1 = M) minimising |N2^-1 G| + cost(N1) + cost(N2)."""
    n = len(M)
    I = [1 << i for i in range(n)]

    def cols_apply(ops, cols):
        out = []
        for g in cols:
            for c, t in ops:
                if g >> c & 1:
                    g ^= 1 << t
            out.append(g)
        return out

    def cost(n2inv_ops):
        Gp = cols_apply(n2inv_ops, G)
        N2inv = lin.apply_cnots(I, n2inv_ops)
        opsN2 = list(reversed(n2inv_ops))
        N1 = [0] * n
        for i in range(n):
            r = 0
            for j in range(n):
                if N2inv[i] >> j & 1:
                    r ^= M[j]
            N1[i] = r
        opsN1 = [] if zero else lin.best_synth(N1, rng, 5)
        return sum(popc(g) for g in Gp) + len(opsN1) + len(opsN2), Gp, opsN1, opsN2

    best = cost([])
    cur_ops = []
    curc = best
    for _ in range(iters):
        ops = list(cur_ops)
        if ops and rng.random() < 0.3:
            ops.pop(rng.randrange(len(ops)))
        else:
            c, t = rng.sample(range(n), 2)
            ops.insert(rng.randrange(len(ops) + 1), (c, t))
        cc = cost(ops)
        if cc[0] <= curc[0]:
            cur_ops, curc = ops, cc
            if cc[0] < best[0]:
                best = cc
    return best[1], best[2], best[3]


def _colops(cols, target, rng):
    ops = lin.transform(cols, target, rng)
    return [(dst, src) for src, dst in reversed(ops)]


def out_track(order, G, M, rng, zero, look=3):
    """backward column-frame tracking -> (ops0, sched); sched[pos] = (injection target indices, ops after)."""
    n = len(M)
    cols = [1 << j for j in range(n)]
    sched = [None] * len(order)
    for pos in range(len(order) - 1, -1, -1):
        g = G[order[pos]]
        S = lin.express(g, cols)
        ops_after = []
        if len(S) > 1:
            fut = [G[order[j]] for j in range(pos - 1, max(-1, pos - 1 - look), -1)]
            base = sum(len(lin.express(f, cols)) for f in fut)
            bestp = None
            for w in S:
                c2 = list(cols)
                for i in S - {w}:
                    c2[w] ^= c2[i]
                v = sum(len(lin.express(f, c2)) for f in fut)
                if v < base and (bestp is None or v < bestp[0]):
                    bestp = (v, w, c2)
            if bestp:
                _, w, c2 = bestp
                ops_after = [(w, i) for i in sorted(S - {w})]
                cols = c2
                S = {w}
        sched[pos] = (sorted(S), ops_after)
    Mcols = [sum(1 << i for i in range(n) if M[i] >> j & 1) for j in range(n)]
    ops0 = [] if zero else _colops(cols, Mcols, rng)
    return ops0, sched


def compile_frame(pb, dec, seed=0, order=None, basis=None, track=False):
    rng = random.Random(seed)
    E = Emit(pb.nq)
    G = [g for _, _, g in dec]
    idx = list(range(len(dec))) if order is None else list(order)
    if track:
        ops0, sched = out_track(idx, G, pb.M, rng, pb.zero)
        opsN1, opsN2 = ops0, []
    else:
        Gp, opsN1, opsN2 = basis if basis else out_basis(pb.M, G, rng, pb.zero)
        sched = [([i for i in range(len(pb.T)) if Gp[k] >> i & 1], []) for k in idx]
    T = pb.T
    for c, t in opsN1:
        E.cx(T[c], T[t])
    RA, RB = _Reg(pb.A), _Reg(pb.B)
    touched = set()
    for pos, k in enumerate(idx):
        al, be, _ = dec[k]
        tgi, after = sched[pos]
        futA = [dec[j][0] for j in idx[pos + 1:pos + 4]]
        futB = [dec[j][1] for j in idx[pos + 1:pos + 4]]
        p = RA.get(al, E, rng, futA)
        q = RB.get(be, E, rng, futB)
        tg = [T[i] for i in tgi]
        if pb.zero and len(tg) == 1 and tg[0] not in touched:
            E.AND(p, q, tg[0])
            touched.add(tg[0])
        else:
            E.AND(p, q, pb.anc)
            for t in tg:
                E.cx(pb.anc, t)
                touched.add(t)
            E.erase(p, q, pb.anc)
        for c, t in after:
            E.cx(T[c], T[t])
            if T[c] in touched:
                touched.add(T[t])
    RA.restore(E, rng)
    RB.restore(E, rng)
    for c, t in opsN2:
        E.cx(T[c], T[t])
    _finish_linear(E, pb)
    return E.text()


class _Null:
    def __init__(self):
        self.n = 0

    def cx(self, c, t):
        self.n += 1


def setup_cost(dec, nA, nB, order, seed=0):
    rng = random.Random(seed)
    E = _Null()
    RA, RB = _Reg(range(nA)), _Reg(range(nB))
    for pos, k in enumerate(order):
        RA.get(dec[k][0], E, rng, [dec[j][0] for j in order[pos + 1:pos + 4]])
        RB.get(dec[k][1], E, rng, [dec[j][1] for j in order[pos + 1:pos + 4]])
    for R in (RA, RB):
        E.n += len(lin.gauss_synth(lin.relmat(R.cur, [1 << i for i in range(len(R.w))])))
    return E.n


def anneal_order(dec, nA, nB, iters=1500, seed=0, deadline=None):
    import time
    rng = random.Random(seed)
    order = list(range(len(dec)))
    rng.shuffle(order)
    cur = setup_cost(dec, nA, nB, order)
    best = (cur, list(order))
    for it in range(iters):
        if deadline and time.time() > deadline:
            break
        Tm = 2.0 * (1 - it / iters) + 0.05
        o = list(order)
        i, j = sorted(rng.sample(range(len(o)), 2))
        if rng.random() < 0.5:
            o[i], o[j] = o[j], o[i]
        else:
            o[i:j + 1] = reversed(o[i:j + 1])
        c = setup_cost(dec, nA, nB, o)
        if c <= cur or rng.random() < math.exp((cur - c) / Tm):
            order, cur = o, c
            if c < best[0]:
                best = (c, list(o))
    return best


# ============================== strategy "walk" (a6) ==========================================================

class _Walk:
    def __init__(self, wires, E, rng=None):
        self.w = list(wires)
        n = len(wires)
        self.n = n
        self.E = E
        self.rng = rng
        self.A = [1 << i for i in range(n)]
        self.G = [1 << i for i in range(n)]

    def cx(self, c, t):
        self.A[t] ^= self.A[c]
        for i in range(self.n):
            if (self.G[i] >> t) & 1:
                self.G[i] ^= 1 << c
        self.E.cx(self.w[c], self.w[t])

    def coords(self, u):
        s = 0
        for i in range(self.n):
            if (u >> i) & 1:
                s ^= self.G[i]
        return s

    def bring(self, u, future):
        s = self.coords(u)
        if popc(s) == 1:
            return self.w[s.bit_length() - 1]
        ws = [j for j in range(self.n) if (s >> j) & 1]
        fut = set(future)
        cand = [j for j in ws if self.A[j] not in fut]
        pool = cand or ws
        t = self.rng.choice(pool) if self.rng else pool[0]
        for j in ws:
            if j != t:
                self.cx(j, t)
        return self.w[t]

    def finish(self, fast=True):
        ops = lin.reduce_ops(self.A, self.n, fast=fast, tries=40)
        for c, t in ops:
            self.cx(c, t)


def _sparsify(Wm, m, min_gain=2):
    W = list(Wm)
    ops = []
    while True:
        best = None
        for i in range(m):
            for j in range(m):
                if i == j:
                    continue
                g = popc(W[i]) - popc(W[i] ^ W[j])
                if g >= min_gain and (best is None or g > best[0]):
                    best = (g, i, j)
        if best is None:
            break
        _, i, j = best
        W[i] ^= W[j]
        ops.append((j, i))
    return W, ops


def compile_walk(pb, dec, order=None, gain=2, seed=None, fast=True):
    """a6 strategy: c-row sparsification + walking operand registers.  Products indexed by r, W rows by target."""
    rng = random.Random(seed) if seed is not None else None
    m = len(pb.T)
    R = len(dec)
    Wm = [sum(((dec[r][2] >> k) & 1) << r for r in range(R)) for k in range(m)]
    Wp, zops = _sparsify(Wm, m, gain + (0 if pb.zero else 1))
    order = list(order) if order is not None else list(range(R))
    E = Emit(pb.nq)
    T = pb.T
    if not pb.zero:
        # c <- Z M c in one CNOT network (Z = row operations of the sparsification); in the clean model the
        # targets are 0 and both steps are skipped
        for c, t in lin.best_synth(lin.apply_cnots(pb.M, zops), random.Random(seed or 0), 5):
            E.cx(T[c], T[t])
    wa = _Walk(pb.A, E, rng)
    wb = _Walk(pb.B, E, rng)
    um = [d[0] for d in dec]
    vm = [d[1] for d in dec]
    touched = set()
    for idx, r in enumerate(order):
        p = wa.bring(um[r], [um[x] for x in order[idx + 1:]])
        q = wb.bring(vm[r], [vm[x] for x in order[idx + 1:]])
        fan = [T[k] for k in range(m) if (Wp[k] >> r) & 1]
        if pb.zero and len(fan) == 1 and fan[0] not in touched:
            E.AND(p, q, fan[0])
            touched.add(fan[0])
            continue
        E.AND(p, q, pb.anc)
        for t in fan:
            E.cx(pb.anc, t)
            touched.add(t)
        E.erase(p, q, pb.anc)
    wa.finish(fast)
    wb.finish(fast)
    for c, t in reversed(zops):
        E.cx(T[c], T[t])
    _finish_linear(E, pb)
    return E.text()


def walk_setup_cost(dec, nA, nB, order):
    """cheap proxy: CX spent on the operand registers by the walk strategy (deterministic pivots)."""
    class _N:
        n = 0

        def cx(self, c, t):
            self.n += 1
    E = _N()
    wa = _Walk(range(nA), E)
    wb = _Walk(range(nB), E)
    um = [d[0] for d in dec]
    vm = [d[1] for d in dec]
    for idx, r in enumerate(order):
        wa.bring(um[r], [um[x] for x in order[idx + 1:]])
        wb.bring(vm[r], [vm[x] for x in order[idx + 1:]])
    wa.finish()
    wb.finish()
    return E.n
