# SPDX-License-Identifier: Apache-2.0
"""Freedom of |+> wires at a block input: the start row of such a wire is arbitrary (completed to a basis).
Change of variables x'' = S0 x: terms p -> coords(S0, p), final rows A' -> coords(S0, A'), path from I."""
from . import gray


def make_S0(n, P, Aprime):
    P = list(P)
    S0 = [1 << w for w in range(n)]
    Pmask = 0
    for t in P: Pmask |= 1 << t
    piv = {}
    chosen = {}
    def red(v):
        for h, b in piv.items():
            if v >> h & 1: v ^= b
        return v
    def add(v):
        v = red(v)
        if v == 0: return False
        h = v.bit_length() - 1
        for k in list(piv):
            if piv[k] >> h & 1: piv[k] ^= v
        piv[h] = v; return True
    for t in P:
        v = Aprime[t] & Pmask
        if v and add(v): chosen[t] = Aprime[t]
    for t in P:
        if t in chosen: continue
        if add(1 << t): chosen[t] = 1 << t
        else:
            for u in P:
                if add(1 << u): chosen[t] = (Aprime[t] & ~Pmask) | (1 << u) if False else (1 << u) | (1 << t) * 0; break
    for t in P:
        if t in chosen: S0[t] = chosen[t]
    try:
        gray.inv_matrix(S0)
    except ValueError:
        return [1 << w for w in range(n)]
    return S0


def transform(S0, O, Aprime):
    O2 = gray.coords(S0, O)
    A2 = gray.coords(S0, Aprime)
    return O2, A2
