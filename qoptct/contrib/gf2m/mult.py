# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""mult: GF(2^m) multiplier circuits  c ^= a*b mod f  (strict)  /  c = a*b mod f  (clean, c starts at 0)
[qoptct v0.5; algorithm and compiler from the qoptct lab, k5 (crt.py, plan.py, rns.py) and x3].

CRT compiler ("residue number system" form).  The CRT places are split into two batches B1 (total degree m) and B2
(the rest, total degree m-1, may include infinity).  For a batch the a and b registers are transformed IN PLACE by a
CNOT network so that every wire holds one residue coordinate of the batch (a -> (a mod p)_{p in B1} is invertible on
m coefficients); the c register is put into the frame in which every wire is one CRT slot of the batch.  Every
product of a place is then one logical AND of short residue forms (short walks inside the place's wires), XORed
into the slot wires and erased by measurement.  One CNOT network per register between and after the batches.
Cost: #AND = rank of the CRT algorithm; T = 4 #AND; qubits 3m + n_anc.

For small m (and whenever it is shorter) the stored flip-graph/tower decompositions of data/decomps.json are used
with the v0.4 compilers of qoptct.contrib.gf2m.compile instead ("table" algorithm).
"""
import random

from . import crt as _crt
from . import linsynth as LS
from .gf2lin import express, rank as f2rank
from .qcirc import GF2Circuit, parse_qasm


def split_batches(m, places, seed=0):
    """B1: finite places with total degree exactly m (subset-sum DP); B2: the rest."""
    deg = [_crt.place_degree(p) for p in places]
    idx = [i for i, p in enumerate(places) if p[0] == 'mod']
    rng = random.Random(seed)
    rng.shuffle(idx)
    reach = {0: []}
    for i in idx:
        for s, L in list(reach.items()):
            t = s + deg[i]
            if t <= m and t not in reach:
                reach[t] = L + [i]
    if m not in reach:
        raise ValueError("no batch of total degree m")
    B1 = reach[m]
    B2 = [i for i in range(len(places)) if i not in B1]
    return B1, B2


class _Reg:
    def __init__(self, wires, C):
        self.w = list(wires)
        self.rows = [1 << i for i in range(len(wires))]
        self.C = C

    def run(self, ops):
        for c, t in ops:
            self.C.cx(self.w[c], self.w[t])
            self.rows[t] ^= self.rows[c]


def _express(form, rows):
    S = express(form, rows)
    if S is None:
        raise AssertionError("form not in span")
    return sorted(S)


def _bring(reg, local, form, future):
    rows = [reg.rows[i] for i in local]
    if form in rows:
        return reg.w[local[rows.index(form)]]
    S = [local[i] for i in _express(form, rows)]
    fut = set(future)
    cand = [i for i in S if reg.rows[i] not in fut] or S
    piv = min(cand, key=lambda i: (sum(1 for f_ in future[:3] if f_ == reg.rows[i]), i))
    for i in S:
        if i != piv:
            reg.run([(i, piv)])
    return reg.w[piv]


def _sim_walk(forms, D):
    rows = [1 << j for j in range(D)]
    cost = 0
    for k, form in enumerate(forms):
        if form in rows:
            continue
        S = _express(form, rows)
        fut = forms[k + 1:k + 4]
        cand = [i for i in S if rows[i] not in fut] or S
        piv = min(cand, key=lambda i: (sum(1 for f_ in fut if f_ == rows[i]), i))
        for i in S:
            if i != piv:
                rows[piv] ^= rows[i]
                cost += 1
    return cost


def _best_local_order(dec, D, seed=0, iters=300):
    rng = random.Random(seed * 7919 + len(dec))

    def cost(o):
        return _sim_walk([dec[i][0] for i in o], D) + _sim_walk([dec[i][1] for i in o], D)
    order = list(range(len(dec)))
    cur = cost(order)
    best = (cur, list(order))
    for _ in range(iters):
        o = list(order)
        i, j = rng.sample(range(len(o)), 2) if len(o) > 1 else (0, 0)
        if rng.random() < 0.5:
            o[i], o[j] = o[j], o[i]
        else:
            o.insert(j, o.pop(i))
        c = cost(o)
        if c <= cur or rng.random() < 0.05:
            order, cur = o, c
            if c < best[0]:
                best = (c, list(o))
    return [dec[i] for i in best[1]]


def _comb(base, mask):
    return _crt._comb(base, mask)


def _complete(vecs, m):
    if len(vecs) < m:
        for i in range(m):
            if f2rank(vecs + [1 << i]) == len(vecs) + 1:
                return vecs + [1 << i], True
    return vecs, False


def build_crt(m, f, clean=False, seed=0, tries=4, places=None, n_anc=1, order=True):
    """CRT multiplier circuit (k5/x3 'rns' compiler).  Wires: a = 0..m-1, b = m..2m-1, c = 2m..3m-1, ancillas
    3m..3m+n_anc-1."""
    if places is None:
        _, places = _crt.plan(m)
    places = list(places)
    cols = _crt.slot_columns(m, f, places)
    B1, B2 = split_batches(m, places, seed)
    nq = 3 * m + n_anc
    C = GF2Circuit(nq)
    RA = _Reg(range(m), C)
    RB = _Reg(range(m, 2 * m), C)
    RC = _Reg(range(2 * m, 3 * m), C)
    anc0 = 3 * m
    npl = 0
    first = True
    for batch in (B1, B2):
        trow, owner = [], []
        for pi in batch:
            for j, r in enumerate(_crt.place_rows(m, places[pi])):
                trow.append(r)
                owner.append((pi, j))
        trow, added = _complete(trow, m)
        if added:
            owner.append(None)
        if len(trow) != m:
            raise AssertionError("batch rows")
        ops, perm = LS.to_perm(RA.rows, trow, tries, seed)
        RA.run(ops)
        wire_of = {owner[i]: perm[i] for i in range(m) if owner[i] is not None}
        opsb, permb = LS.to_perm(RB.rows, trow, tries, seed)
        RB.run(opsb)
        wire_of_b = {owner[i]: permb[i] for i in range(m) if owner[i] is not None}
        ccols, cown = [], []
        for pi in batch:
            for j in range(_crt.place_degree(places[pi])):
                ccols.append(cols[(pi, j)])
                cown.append((pi, j))
        ccols, added = _complete(ccols, m)
        if added:
            cown.append(None)
        Kt = LS.inverse(LS.transpose(ccols, m))          # rows of K, K^-1 has the slot columns
        if clean and first:
            cperm = list(range(m))                        # c = 0: the frame change is a relabelling only
            RC.rows = list(Kt)
        else:
            cops, cperm = LS.to_perm(RC.rows, Kt, tries, seed)
            RC.run(cops)
        cwire = {cown[i]: cperm[i] for i in range(m) if cown[i] is not None}
        zero_ok = clean and first
        touched = set()
        first = False
        for pi in batch:
            pl = places[pi]
            anc = anc0 + (npl % n_anc)
            npl += 1
            D = _crt.place_degree(pl)
            la = [wire_of[(pi, j)] for j in range(D)]
            lb = [wire_of_b[(pi, j)] for j in range(D)]
            base = _crt.place_rows(m, pl)
            dec = _crt.ring_dec(_crt._ring_key(pl))
            if order:
                dec = _best_local_order(dec, D, seed)
            forms_a = [_comb(base, al) for al, be, ga in dec]
            forms_b = [_comb(base, be) for al, be, ga in dec]
            for k, (al, be, ga) in enumerate(dec):
                if not forms_a[k] or not forms_b[k]:
                    continue
                p = _bring(RA, la, forms_a[k], forms_a[k + 1:k + 4])
                q = _bring(RB, lb, forms_b[k], forms_b[k + 1:k + 4])
                tg = [RC.w[cwire[(pi, j)]] for j in range(D) if ga >> j & 1]
                if not tg:
                    continue
                if zero_ok and len(tg) == 1 and tg[0] not in touched:
                    C.AND(p, q, tg[0])
                    touched.add(tg[0])
                    continue
                C.AND(p, q, anc)
                for t in tg:
                    C.cx(anc, t)
                    touched.add(t)
                C.erase(p, q, anc)
    ident = [1 << i for i in range(m)]
    RA.run(LS.to_exact(RA.rows, ident, tries, seed))
    RB.run(LS.to_exact(RB.rows, ident, tries, seed))
    RC.run(LS.to_exact(RC.rows, ident, tries, seed))
    if RA.rows != ident or RB.rows != ident or RC.rows != ident:
        raise AssertionError("registers not restored")
    C.regs = dict(a=list(range(m)), b=list(range(m, 2 * m)), c=list(range(2 * m, 3 * m)))
    return C, tuple(places)


# ------------------------------------------------------------------------------ CRT products inside larger circuits

def _pull(mask, cols):
    """beta with mask.(L x) = beta.x for the linear map L with columns cols (None = identity)."""
    if cols is None:
        return mask
    b = 0
    for j, c in enumerate(cols):
        if bin(mask & c).count("1") & 1:
            b |= 1 << j
    return b


class CRTData:
    """per-field data of the CRT compiler (places, slot columns, batches, ordered ring algorithms)."""

    def __init__(self, m, f, places, seed=0, order=True):
        self.m, self.f, self.places = m, f, list(places)
        self.cols = _crt.slot_columns(m, f, self.places)
        self.batches = split_batches(m, self.places, seed)
        self.base = [_crt.place_rows(m, pl) for pl in self.places]
        self.decs = []
        for pl in self.places:
            dec = _crt.ring_dec(_crt._ring_key(pl))
            if order:
                dec = _best_local_order(dec, _crt.place_degree(pl), seed)
            self.decs.append(dec)
        self.rows = []
        for batch in self.batches:
            trow, owner = [], []
            for pi in batch:
                for j, r in enumerate(self.base[pi]):
                    trow.append(r)
                    owner.append((pi, j))
            trow, added = _complete(trow, m)
            if added:
                owner.append(None)
            ccols, cown = [], []
            for pi in batch:
                for j in range(_crt.place_degree(self.places[pi])):
                    ccols.append(self.cols[(pi, j)])
                    cown.append((pi, j))
            ccols, added = _complete(ccols, m)
            if added:
                cown.append(None)
            Kt = LS.inverse(LS.transpose(ccols, m))
            self.rows.append((trow, owner, Kt, cown))
        self._syn = {}

    def to_perm(self, cur, target, tries, seed):
        key = ('p', tuple(cur), tuple(target))
        if key not in self._syn:
            if len(self._syn) > 4096:
                self._syn.clear()
            self._syn[key] = LS.to_perm(cur, target, tries, seed)
        return self._syn[key]

    def to_exact(self, cur, target, tries, seed):
        key = ('e', tuple(cur), tuple(target))
        if key not in self._syn:
            if len(self._syn) > 4096:
                self._syn.clear()
            self._syn[key] = LS.to_exact(cur, target, tries, seed)
        return self._syn[key]


def emit_product(C, D, wa, wb, wt, anc, cols_a=None, cols_b=None, mode="mul", zero_target=False, mids=None,
                 same=False, tries=4, seed=0):
    """CRT product inside a circuit C (GF2Circuit or itoh.StreamCircuit), registers given by wire lists.
    mode 'mul':   wt ^= F_a(A) * F_b(B)  (F = linear maps with columns cols_*, e.g. Frobenius powers;
                  zero_target: wt starts at 0, the first frame change is free).
    mode 'erase': the phase part of a register erasure: wt has been measured (outcome ids `mids`, one per bit);
                  applies (-1)^(s . F_a(A) F_b(B)) as classically controlled CZ / Z (mids None: only the CX).
    same=True: wb is a copy of A (the caller made it); a term whose two pulled forms coincide is linear (no AND).
    A, B (and wt in 'mul') are transformed in place and restored.  -> (ANDs emitted, conditional phases)."""
    m = D.m
    RA, RB = _Reg(wa, C), _Reg(wb, C)
    RC = _Reg(wt, C) if mode == "mul" else None
    n_and = 0
    n_ph = 0
    first = True
    for bi, batch in enumerate(D.batches):
        trow, owner, Kt, cown = D.rows[bi]
        ops, perm = D.to_perm(RA.rows, [_pull(r, cols_a) for r in trow], tries, seed)
        RA.run(ops)
        wire_a = {owner[i]: perm[i] for i in range(m) if owner[i] is not None}
        ops, perm = D.to_perm(RB.rows, [_pull(r, cols_b) for r in trow], tries, seed)
        RB.run(ops)
        wire_b = {owner[i]: perm[i] for i in range(m) if owner[i] is not None}
        zero_ok = False
        if mode == "mul":
            if zero_target and first:
                cperm = list(range(m))
                RC.rows = list(Kt)
                zero_ok = True
            else:
                cops, cperm = D.to_perm(RC.rows, Kt, tries, seed)
                RC.run(cops)
            cwire = {cown[i]: cperm[i] for i in range(m) if cown[i] is not None}
        first = False
        touched = set()
        for pi in batch:
            Dg = _crt.place_degree(D.places[pi])
            la = [wire_a[(pi, j)] for j in range(Dg)]
            lb = [wire_b[(pi, j)] for j in range(Dg)]
            base_a = [_pull(r, cols_a) for r in D.base[pi]]
            base_b = [_pull(r, cols_b) for r in D.base[pi]]
            dec = D.decs[pi]
            fa = [_comb(base_a, al) for al, be, ga in dec]
            fb = [_comb(base_b, be) for al, be, ga in dec]
            for k, (al, be, ga) in enumerate(dec):
                if not fa[k] or not fb[k]:
                    continue
                p = _bring(RA, la, fa[k], fa[k + 1:k + 4])
                q = _bring(RB, lb, fb[k], fb[k + 1:k + 4])
                lin_term = same and fa[k] == fb[k]
                if mode == "mul":
                    tg = [RC.w[cwire[(pi, j)]] for j in range(Dg) if ga >> j & 1]
                    if not tg:
                        continue
                    if lin_term:
                        for t in tg:
                            C.cx(p, t)
                            touched.add(t)
                        continue
                    n_and += 1
                    if zero_ok and len(tg) == 1 and tg[0] not in touched:
                        C.AND(p, q, tg[0])
                        touched.add(tg[0])
                        continue
                    C.AND(p, q, anc)
                    for t in tg:
                        C.cx(anc, t)
                        touched.add(t)
                    C.erase(p, q, anc)
                else:
                    g = 0
                    for j in range(Dg):
                        if ga >> j & 1:
                            g ^= D.cols[(pi, j)]
                    n_ph += bin(g).count("1")
                    if mids is not None:
                        for i in range(m):
                            if g >> i & 1:
                                C.cphase(mids[i], (p,) if lin_term else (p, q))
    ident = [1 << i for i in range(m)]
    RA.run(D.to_exact(RA.rows, ident, tries, seed))
    RB.run(D.to_exact(RB.rows, ident, tries, seed))
    if RC is not None:
        RC.run(D.to_exact(RC.rows, ident, tries, seed))
    return n_and, n_ph


def table_decs(m, f):
    """stored decompositions (data/decomps.json) of F2[x]/(f), transported to f; [] if none."""
    from . import _data
    from . import field as _field
    table = _data("decomps.json")
    out = []
    if not table:
        return out
    for ent in table["fields"]:
        if ent["m"] != m:
            continue
        for d in ent["decs"]:
            try:
                dd = _field.change_basis([tuple(t) for t in d["dec"]], m, d["f"], f)
            except ValueError:
                continue
            if _crt.check_algorithm(dd, m, f):
                out.append((dd, d["src"]))
    return out


def build_table(m, f, dec, clean=False, seeds=(0, 1, 2)):
    """multiplier from an explicit bilinear algorithm with the v0.4 compilers (frame / walk); shortest CX wins."""
    from . import compile as _comp
    A, B, T = list(range(m)), list(range(m, 2 * m)), list(range(2 * m, 3 * m))
    M = [1 << i for i in range(m)]
    lin = {t: [[], 0] for t in T}
    pb = _comp.Problem(A, B, T, M, lin, clean, 3 * m, 3 * m + 1)
    best = None
    order = list(range(len(dec)))
    texts = []
    for s in seeds:
        texts.append(_comp.compile_walk(pb, dec, order=order, gain=2, seed=s))
        texts.append(_comp.compile_frame(pb, dec, seed=s, order=order, track=True))
        texts.append(_comp.compile_frame(pb, dec, seed=s, order=order))
    for txt in texts:
        nq, ops = parse_qasm(txt)
        C = GF2Circuit(nq, ops)
        cx = C.counts()["CX"]
        if best is None or cx < best[0]:
            best = (cx, C)
    C = best[1]
    C.regs = dict(a=A, b=B, c=T)
    return C
