"""c1/csyn.py: phase-polynomial body synthesis with FREE Clifford absorption.

Body = T-parity network (only odd parities are synthesised) + Clifford remainder R (even weights).
Key observation: R/2 is a Z4 quadratic form.  Its mod-2 'symmetric matrix' B (diag = singles/2 mod 2, offdiag = pairs/4)
can be produced for free by S gates on ANY parity l that is present on some wire at some moment of the network
(contributes l l^T), incl. the start rows Winit, the end rows Wend and every intermediate wire value.  CZ between two wires
holding a,b contributes a b^T + b a^T (cost 1).  After that only singles with weight 4 remain -> Z gates on start rows
(Winit is a basis, always solvable, cost 0).  So the Clifford remainder costs CZ only for the part of B outside
span{l l^T : l available}.
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, random
for p in ((_QPKG + "/contrib/tg3"), (_QPKG + "/contrib/eng1"), _QPKG,
          (_QPKG + "/portfolio")):
    if p not in sys.path: sys.path.insert(0, p)
from core import *            # tg3 core (pp, coeffs, canon, evo_lib, det2, ...)


def bitsof(m):
    return [k for k in range(m.bit_length()) if m >> k & 1]


def llT(l, N):
    b = bitsof(l); v = 0
    for x in range(len(b)):
        i = b[x]
        for y in range(x, len(b)):
            v |= 1 << (i * N + b[y])
    return v


def abT(a, b, N):
    return llT(a, N) ^ llT(b, N) ^ llT(a ^ b, N)


def target_B(R, N):
    B = 0
    for k, v in R.items():
        v %= 8
        if not v: continue
        if len(k) == 1:
            if v % 2: return None
            if (v // 2) % 2: B |= 1 << (k[0] * N + k[0])
        elif len(k) == 2:
            if v % 4: return None
            if v == 4: B |= 1 << (k[0] * N + k[1])
        else:
            return None
    return B


class Elim:
    def __init__(self):
        self.rows = {}

    def reduce(self, v, comb=0):
        while v:
            p = v.bit_length() - 1
            r = self.rows.get(p)
            if r is None: return v, comb
            v ^= r[0]; comb ^= r[1]
        return 0, comb

    def full_reduce(self, v, comb=0):
        """reduce every bit that has a pivot (not only the leading one) -> canonical representative."""
        out = 0
        while v:
            p = v.bit_length() - 1
            r = self.rows.get(p)
            if r is None: out |= 1 << p; v ^= 1 << p
            else: v ^= r[0]; comb ^= r[1]
        return out, comb

    def add(self, v, comb):
        v, comb = self.reduce(v, comb)
        if v: self.rows[v.bit_length() - 1] = (v, comb); return True
        return False


def poly_of(body, Winit, N):
    W = list(Winit); K = [0] * N; P = {}
    pp.block_poly([g for g in body if g[0] != "x"], W, K, P)
    return P, W, K


def sub_canon(A, B):
    C = dict(A)
    for k, v in B.items(): C[k] = (C.get(k, 0) - v) % 8
    return {k: v for k, v in C.items() if v}


EXACT = [False]


def shrink_cz(r, czs, red, E, N, exact=False):
    """local search: replace 2 CZs by 1 (and, exact=True, 3 by 2) when the reduced vectors allow it."""
    rep = lambda c: E.full_reduce(abT(c[3], c[4], N))[0]
    sol = [rep(c) for c in czs]
    imp = True
    pairsum = None
    if exact and len(red) <= 2500:
        keys = list(red); pairsum = {}
        for i in range(len(keys)):
            ki = keys[i]
            for j in range(i + 1, len(keys)):
                pairsum.setdefault(ki ^ keys[j], (ki, keys[j]))
    while imp:
        imp = False
        L = len(sol)
        for i in range(L):
            for j in range(i + 1, L):
                v = sol[i] ^ sol[j]
                if v == 0:
                    sol = [x for k, x in enumerate(sol) if k not in (i, j)]; imp = True; break
                if v in red:
                    sol = [x for k, x in enumerate(sol) if k not in (i, j)] + [v]; imp = True; break
            if imp: break
        if imp or pairsum is None: continue
        for i in range(L):
            for j in range(i + 1, L):
                for k in range(j + 1, L):
                    v = sol[i] ^ sol[j] ^ sol[k]
                    if v in pairsum:
                        a, b = pairsum[v]
                        sol = [x for t, x in enumerate(sol) if t not in (i, j, k)] + [a, b]; imp = True; break
                if imp: break
            if imp: break
    acc = 0
    for v in sol: acc ^= v
    if acc != r: return czs
    return [red[v] for v in sol]


def clifford_absorb(body0, Winit, Ptarget, N, cz_mode="any", canP=None):
    """body0: CX + T-type gates starting from wires Winit (implements the T part).  Inserts free S gates on available
    parities + minimal CZs + Z on start rows so that the body implements Ptarget exactly.  -> body or None."""
    P0, Wend, _ = poly_of(body0, Winit, N)
    R = sub_canon(canP if canP is not None else canon(Ptarget), canon(P0))
    B = target_B(R, N)
    if B is None: return None
    seen = {}
    W = list(Winit)
    for q in range(N):
        if W[q] not in seen: seen[W[q]] = (0, q)
    snaps = [(0, list(W))]
    for idx, (g, qs) in enumerate(body0):
        if g == "cx":
            W[qs[1]] ^= W[qs[0]]
            t = qs[1]
            if W[t] not in seen: seen[W[t]] = (idx + 1, t)
            snaps.append((idx + 1, list(W)))
    gens = list(seen.items())
    E = Elim()
    for gi, (l, _) in enumerate(gens):
        E.add(llT(l, N), 1 << gi)
    r, comb = E.full_reduce(B)
    czs = []
    if r:
        sn = snaps if cz_mode == "any" else [snaps[0], snaps[-1]]
        red = {}
        for pos, Ws in sn:
            for a in range(N):
                for b in range(a + 1, N):
                    v, _ = E.full_reduce(abT(Ws[a], Ws[b], N))
                    if v and v not in red: red[v] = (pos, a, b, Ws[a], Ws[b])
        if r in red:
            czs = [red[r]]
        else:
            for v in list(red):
                if (r ^ v) in red:
                    czs = [red[v], red[r ^ v]]; break
            if not czs:
                E2 = Elim(); lst = list(red.items())
                for j, (v, _) in enumerate(lst): E2.add(v, 1 << j)
                rr, cc = E2.reduce(r)
                if rr: return None
                czs = [lst[j][1] for j in bitsof(cc)]
        if len(czs) >= 2:
            czs = shrink_cz(r, czs, red, E, N, exact=EXACT[0])
        B2 = B
        for c in czs: B2 ^= abT(c[3], c[4], N)
        r2, comb = E.reduce(B2)
        if r2: return None
    ins = {}
    for gi in bitsof(comb):
        l, (pos, w) = gens[gi]
        ins.setdefault(pos, []).append(("s", (w,)))
    for (pos, a, b, _, _) in czs:
        ins.setdefault(pos, []).append(("cz", (a, b)))
    body = []
    for idx in range(len(body0) + 1):
        body += ins.get(idx, [])
        if idx < len(body0): body.append(body0[idx])
    P1, _, _ = poly_of(body, Winit, N)
    D = sub_canon(canP if canP is not None else canon(Ptarget), canon(P1))
    z = 0
    for k, v in D.items():
        if len(k) != 1 or v != 4: return None
        z |= 1 << k[0]
    pre = []
    if z:
        c = pp.solve(Winit, N, z)
        pre = [("z", (q,)) for q in bitsof(c)]
    return merge_phases(pre + body)


W8 = {"t": 1, "s": 2, "z": 4, "sdg": 6, "tdg": 7}
SEQ8 = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["z"], 5: ["z", "t"], 6: ["sdg"], 7: ["tdg"]}


def merge_phases(body):
    out = []; acc = {}
    def flush(q):
        v = acc.pop(q, 0) % 8
        out.extend((g, (q,)) for g in SEQ8[v])
    for g, qs in body:
        if g in W8: acc[qs[0]] = acc.get(qs[0], 0) + W8[g]
        elif g == "cz": out.append((g, qs))
        else:
            for q in qs:
                if q in acc: flush(q)
            out.append((g, qs))
    for q in list(acc): flush(q)
    return out


def cost(body):
    return sum(1 for g, _ in body if g in ("cx", "cz"))


def t_parities(P):
    return [m for m, w in P.items() if w % 2]


def net_variants(PT, Winit, Wend, N, rng, s):
    W = list(Winit)
    try:
        if s % 3 == 0: net = pp.graysynth(dict(PT), W, N, rng)
        elif s % 3 == 1: net = pp.greedy(dict(PT), W, N, rng, Wend, (0.0, 0.3, 0.6, 1.0, 1.5)[(s // 3) % 5])
        else: net = pp.greedy(dict(PT), W, N, rng, Wend, rng.uniform(0, 1.5), rng.uniform(0, 1.5), rng.uniform(0.1, 1.0))
    except Exception:
        net = None
    return net, W


def synth_abs(ctx, Winit, Wend, seeds=8, base_seed=0, P=None, lin_tries=6):
    N, n, K = ctx["N"], ctx["n"], ctx["K"]
    P = P if P is not None else ctx["P2"]
    PT = {m: 1 for m in t_parities(P)}
    canP = canon(P)
    best = None
    for s in range(seeds):
        rng = random.Random(base_seed + s)
        net, W = net_variants(PT, Winit, Wend, N, rng, s)
        if net is None: continue
        M = [pp.solve(Wend, N, W[r]) for r in range(N)]
        ops = pp.lin_best(M, N, range(N), rng, lin_tries)
        body = clifford_absorb(list(net) + [("cx", o) for o in ops], Winit, P, N, canP=canP)
        if body is None: continue
        body = body + [("x", (q,)) for q in range(n) if K[q]]
        c = cost(body)
        if best is None or c < best[0]: best = (c, body)
    return None if best is None else best[1]
