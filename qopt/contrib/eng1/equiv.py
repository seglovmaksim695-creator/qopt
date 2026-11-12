"""Fibre-preserving change of variables (the freedom that plain TODD does not use).
Amplitude  A(o|x) = sum_{y : L(x,y)=o} w^{f(x,y)}.  For any invertible linear phi on F_2^N with  x∘phi = x  and  L∘phi = L
we have  A_f = A_{f∘phi}  (substitution v = phi(u) inside the sum).  Elementary moves: transvections
  phi(v) = v + w (l.v),  w in y-space, L w = 0, l.w = 0   (acts on parity masks as  m -> m ^ (m.w) l ).
After a move the polynomial is re-reduced by TODD (pointwise-exact).  Whole chain is certified by replay (certify)."""
import random, time
import numpy as np
import gad
from gad import pp, split_odd, clifford_rest, canon, rank_ok
from fasttodd import todd_fast

par = lambda x: bin(x).count("1") & 1


def kernel_dirs(Lrows, n, N):
    """Basis (as N-bit masks supported on y-vars n..N-1) of {w : <L_i, w> = 0 for all output rows}."""
    h = N - n
    rows = [(Lm >> n) for Lm in Lrows]  # restricted to y
    # nullspace over GF(2) of matrix rows (each an h-bit int)
    piv = {}
    for r in rows:
        for p, pr in piv.items():
            if r >> p & 1: r ^= pr
        if r:
            p = r.bit_length() - 1
            for q in list(piv):
                if piv[q] >> p & 1: piv[q] ^= r
            piv[p] = r
    free = [j for j in range(h) if j not in piv]
    out = []
    for f in free:
        y = 1 << f
        for p, pr in piv.items():
            if pr >> f & 1: y |= 1 << p
        out.append(y << n)
    return out


def relabel(P, w, ell):
    out = {}
    for m, wt in P.items():
        m2 = m ^ ell if par(m & w) else m
        out[m2] = (out.get(m2, 0) + wt) % 8
    return {m: v for m, v in out.items() if v}


def certify_chain(P0, W0, n, N, steps):
    """Check every transvection (L w=0, w in y-space, l.w=0) and return the replayed canonical polynomial."""
    P = dict(P0)
    for w, ell in steps:
        assert w >> n << n == w and w, "w must live on y variables"
        assert all(par(Lm & w) == 0 for Lm in W0[:n]), "L w != 0"
        assert par(ell & w) == 0, "l.w != 0"
        P = relabel(P, w, ell)
    return canon(P)


def search(P, W, N, n, rng, total_time=60.0, todd_tl=3.0, tries=120, log=None, anneal=0.0):
    """Random walk over (transvection, TODD) moves; returns (P_best, steps_best)."""
    t0 = time.time()
    dirs = kernel_dirs(W[:n], n, N)
    odd, even = split_odd(P)
    cur = (dict(P), len(odd), []); best = cur
    if log: log(f"   equiv: kernel dirs {len(dirs)}")
    if not dirs: return best[0], best[2]
    while time.time() - t0 < total_time:
        Pc, mc, steps = cur
        w = 0
        for d in rng.sample(dirs, min(len(dirs), rng.choice((1, 1, 2)))): w ^= d
        if not w: continue
        sup = [i for i in range(N) if w >> i & 1]
        r = rng.random()
        if r < 0.5: ell = 1 << rng.randrange(N)
        elif r < 0.8: ell = (1 << rng.randrange(N)) ^ (1 << rng.randrange(N))
        else:
            cols = [m for m, wt in Pc.items() if wt % 2]
            ell = rng.choice(cols) ^ rng.choice(cols)
        if par(ell & w): ell ^= 1 << rng.choice(sup)
        if ell == 0 or par(ell & w): continue
        P2 = relabel(Pc, w, ell)
        odd2, even2 = split_odd(P2)
        new = todd_fast(odd2, N, rng, todd_tl, tries=tries, plateau_max=0) if False else todd_fast(odd2, N, rng, todd_tl, 0, tries)
        if len(new) <= mc + anneal * (rng.random() < 0.1):
            P3 = clifford_rest(odd2, new, even2)
            cur = (P3, len(new), steps + [(w, ell)])
            if len(new) < best[1] or (len(new) == best[1] and False): best = cur
            if log and len(new) < mc: log(f"   equiv: T {mc}->{len(new)} t={time.time()-t0:.0f}s")
    return best[0], best[2]
