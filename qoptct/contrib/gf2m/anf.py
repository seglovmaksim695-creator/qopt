# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""Classical function of a reversible Feynman .qc circuit as GF(2) polynomials (algebraic normal form) and
extraction of the bilinear structure  y <- M y + Lin(x) + c + B(a, b).

A polynomial is a Python set of monomials; a monomial is an int bitmask of variables (0 = constant 1).
strict model: variable k = wire k.  clean model: variable k = k-th wire of the .i line; other wires start at 0.
"""


class NotClassical(Exception):
    pass


def read_qc(path):
    """-> (names, gates, inputs) ; gates as (name, qubit tuple); inputs = wire indices of the .i line (or all)."""
    from qoptct.qopt_core import load_qc
    names, gates = load_qc(path)
    ins = None
    for line in open(path):
        p = line.split()
        if p and p[0] == ".i":
            ins = p[1:]
    inp = [names.index(x) for x in ins] if ins is not None else list(range(len(names)))
    return names, gates, inp


def logical(gates, n):
    """H-frames -> classical gates x, cx, ccx and phase gates z, cz, ccz. Raises NotClassical otherwise."""
    fr = [0] * n
    out = []
    for g, q in gates:
        q = tuple(q)
        if g == "h":
            fr[q[0]] ^= 1
            continue
        if g in ("ccz", "cz", "z"):
            F = [x for x in q if fr[x]]
            if not F:
                out.append((g, q))
            elif len(F) == 1:
                t = F[0]
                c = tuple(x for x in q if x != t)
                out.append(({"ccz": "ccx", "cz": "cx", "z": "x"}[g], c + (t,)))
            else:
                raise NotClassical("two H-frames in one phase gate")
        elif g == "x":
            out.append(("z" if fr[q[0]] else "x", q))
        elif g in ("cx", "ccx"):
            if any(fr[x] for x in q[:-1]):
                raise NotClassical("control inside an H-frame")
            out.append((({"cx": "cz", "ccx": "ccz"}[g]) if fr[q[-1]] else g, q))
        else:
            raise NotClassical(g)
    if any(fr):
        raise NotClassical("open H-frame")
    return out


def pmul(p, q):
    r = set()
    for a in p:
        for b in q:
            r ^= {a | b}
    return r


def initial(n, inp, clean):
    if clean:
        P = [set() for _ in range(n)]
        for k, w in enumerate(inp):
            P[w] = {1 << k}
        return P
    return [{1 << w} for w in range(n)]


MAX_TERMS = 20000   # per polynomial; beyond this the function is far from bilinear (and ANF may explode)
MAX_DEGREE = 16


def anf(gates, n, inp, clean):
    """-> (P, ph): final polynomial of each wire and the +-1 phase polynomial.
    Raises NotClassical for non-classical circuits or when the ANF grows beyond MAX_TERMS / MAX_DEGREE."""
    P = initial(n, inp, clean)
    ph = set()
    for g, q in logical(gates, n):
        if g == "x":
            P[q[0]] ^= {0}
        elif g == "cx":
            P[q[1]] ^= P[q[0]]
        elif g == "ccx":
            if len(P[q[0]]) * len(P[q[1]]) > MAX_TERMS:
                raise NotClassical("ANF too large")
            pr = pmul(P[q[0]], P[q[1]])
            if any(bin(mm).count("1") > MAX_DEGREE for mm in pr):
                raise NotClassical("degree too high")
            P[q[2]] ^= pr
            if len(P[q[2]]) > MAX_TERMS:
                raise NotClassical("ANF too large")
        elif g == "z":
            ph ^= P[q[0]]
        elif g == "cz":
            ph ^= pmul(P[q[0]], P[q[1]])
        elif g == "ccz":
            ph ^= pmul(pmul(P[q[0]], P[q[1]]), P[q[2]])
        else:
            raise NotClassical(g)
    return P, ph


def structure(P, n, inp, clean):
    """Bilinear structure of the function or None.

    Returns dict(T=target wires, M={t: set of target wires}, lin={t: [set of source wires, const]},
    quad=[(s1, s2, t)], zero=bool) where t_new = sum_{u in M[t]} u + sum lin + const + sum_{(s1,s2,t)} s1 s2.
    zero=True (clean only): every target starts at 0 and has no target-linear part."""
    init = initial(n, inp, clean)
    wire_of_var = list(inp) if clean else list(range(n))
    T = [w for w in range(n) if P[w] != init[w]]
    if not T:
        return None
    Tset = set(T)
    M, lin, quad = {}, {}, []
    for t in T:
        M[t] = set()
        lin[t] = [set(), 0]
        for m in P[t]:
            k = bin(m).count("1")
            if k == 0:
                lin[t][1] ^= 1
                continue
            if k == 1:
                w = wire_of_var[m.bit_length() - 1]
                if w in Tset:
                    M[t].add(w)
                else:
                    lin[t][0].add(w)
            elif k == 2:
                i = (m & -m).bit_length() - 1
                j = (m ^ (1 << i)).bit_length() - 1
                wi, wj = wire_of_var[i], wire_of_var[j]
                if wi in Tset or wj in Tset:
                    return None
                quad.append((min(wi, wj), max(wi, wj), t))
            else:
                return None
    if not quad:
        return None
    zero = clean and all(t not in inp for t in T) and all(not M[t] for t in T)
    if clean and not zero:
        # some targets are inputs (e.g. a wire listed in .i by the benchmark): compile in place; a target that
        # starts at 0 keeps its own (zero) value as the identity part, which makes M invertible when possible
        inps = set(inp)
        for t in T:
            if t not in inps:
                M[t].add(t)
    if not zero:
        # the target-linear part must be invertible (it is applied in place)
        from .gf2lin import rank
        rows = [sum(1 << T.index(u) for u in M[t]) for t in T]
        if rank(rows) != len(T):
            return None
    return dict(T=T, M=M, lin=lin, quad=quad, zero=zero, n=n)


def bipartition(st):
    """split the source wires of quadratic monomials into (A, B), each monomial with one wire in each."""
    adj = {}
    for i, j, _ in st["quad"]:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    col = {}
    for s0 in sorted(adj):
        if s0 in col:
            continue
        col[s0] = 0
        stack = [s0]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in col:
                    col[v] = 1 - col[u]
                    stack.append(v)
                elif col[v] == col[u]:
                    return None
    A = sorted(v for v in col if col[v] == 0)
    B = sorted(v for v in col if col[v] == 1)
    if len(A) > len(B) or (len(A) == len(B) and A > B):
        A, B = B, A
    return A, B


def tensor(st):
    """-> (A, B, T, terms) with terms[(i, j)] = mask over T of targets receiving a_i b_j."""
    AB = bipartition(st)
    if AB is None:
        return None
    A, B = AB
    T = st["T"]
    Ai = {w: k for k, w in enumerate(A)}
    Bi = {w: k for k, w in enumerate(B)}
    Ti = {w: k for k, w in enumerate(T)}
    terms = {}
    for i, j, t in st["quad"]:
        if i in Bi:
            i, j = j, i
        key = (Ai[i], Bi[j])
        terms[key] = terms.get(key, 0) ^ (1 << Ti[t])
    terms = {k: v for k, v in terms.items() if v}
    return A, B, T, terms


def tensor_of(dec):
    """tensor of a decomposition [(alpha, beta, gamma)] as {(i, j): gamma mask}."""
    tt = {}
    for a, b, c in dec:
        ia = [i for i in range(a.bit_length()) if a >> i & 1]
        jb = [j for j in range(b.bit_length()) if b >> j & 1]
        for i in ia:
            for j in jb:
                tt[(i, j)] = tt.get((i, j), 0) ^ c
    return {k: v for k, v in tt.items() if v}
