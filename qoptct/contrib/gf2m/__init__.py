# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""gf2m: bilinear functions (z ^= B(x, y)) with few logical ANDs.

Recognition is by function: the reversible .qc is turned into its ANF (strict: all wires are variables;
clean: wires missing from the .i line start at 0).  If every changed wire t satisfies
    t_new = M t + Lin(sources) + const + sum of products a_i b_j        (a, b disjoint source registers)
the map is a bilinear form.  Multiplication in GF(2^m) (polynomial basis, any wire order, any irreducible
modulus) is detected and synthesised from stored short bilinear decompositions (flip graph, tower + CRT;
data/decomps.json) transported to the circuit's modulus; other bilinear forms get a slice decomposition
improved by a flip-graph search.  Each product becomes one logical AND (4 T) with measurement-based erasure.
A small table of verified constructions for specific functions (ham15, data/known.json) is matched by ANF.

API (qoptct v0.4 contrib interface):
    candidates(qc_path, clean, budget_sec=60.0) -> [{"qasm": str, "note": str}, ...]
Every returned circuit passed the exact symbolic check of check.py (all inputs, all measurement outcomes) and,
for small circuits, the state-vector check of qoptct.verify_qasm.  Unrecognised input -> [].
"""
import json
import math
import os
import random
import time

from . import anf as _anf
from . import field as _field
from . import compile as _comp
from . import check as _check
from . import known as _known
from . import flip as _flip

_HERE = os.path.dirname(os.path.abspath(__file__))
_CACHE = {}
_BASIS = {}
NUMERIC_MAX_QUBITS = 16


def _data(name):
    if name not in _CACHE:
        try:
            with open(os.path.join(_HERE, "data", name)) as fh:
                _CACHE[name] = json.load(fh)
        except (OSError, ValueError):
            _CACHE[name] = None
    return _CACHE[name]


# ------------------------------------------------------------------------------------------------- analysis

def analyse(qc_path, clean):
    """-> dict with the reference function and (if any) the bilinear structure; None if not classical."""
    names, gates, inp = _anf.read_qc(qc_path)
    n = len(names)
    try:
        P, ph = _anf.anf(gates, n, inp, clean)
    except _anf.NotClassical:
        return None
    info = dict(names=names, n=n, inp=inp, P=P, ph=ph, clean=clean, st=None)
    if ph:
        return info
    st = _anf.structure(P, n, inp, clean)
    if st is None:
        return info
    tn = _anf.tensor(st)
    if tn is None:
        return info
    info["st"] = st
    info["tensor"] = tn
    return info


def _problem(info):
    st = info["st"]
    A, B, T, _ = info["tensor"]
    n, inp, P = info["n"], info["inp"], info["P"]
    anc, nq = n, n + 1
    if info["clean"]:
        used = set(A) | set(B) | set(T)
        free = [w for w in range(n) if w not in inp and w not in used and not P[w]]
        if free:
            anc, nq = free[0], n
    M = [sum(1 << T.index(u) for u in st["M"][t]) for t in T]
    return _comp.Problem(A, B, T, M, st["lin"], st["zero"], anc, nq)


def _slice_decomposition(nA, nB, nT, terms):
    """rank decomposition of the slices along the cheapest axis."""
    best = None
    for axis in range(3):
        out = []
        size = (nA, nB, nT)[axis]
        for s in range(size):
            # matrix rows: (other index u) -> vector over the third index v
            rows = {}
            for (i, j), c in terms.items():
                if axis == 0 and i == s:
                    rows[j] = c                                   # b_j -> targets
                elif axis == 1 and j == s:
                    rows[i] = c                                   # a_i -> targets
                elif axis == 2 and c >> s & 1:
                    rows[i] = rows.get(i, 0) | (1 << j)           # a_i -> b mask
            if not rows:
                continue
            # basis of the row space, coefficients of every row
            basis = []
            for v in rows.values():
                if _field_express(v, basis) is None:
                    basis.append(v)
            coef = {u: _field_express(v, basis) for u, v in rows.items()}
            for k, v in enumerate(basis):
                um = sum(1 << u for u, S in coef.items() if k in S)
                if axis == 0:
                    out.append((1 << s, um, v))
                elif axis == 1:
                    out.append((um, 1 << s, v))
                else:
                    out.append((um, v, 1 << s))
        if best is None or len(out) < len(best):
            best = out
    return best


def _field_express(v, basis):
    from .gf2lin import express
    return express(v, basis)


def decompositions(info, deadline):
    """list of (dec over tensor indices, source note, field info or None)."""
    A, B, T, terms = info["tensor"]
    out = []
    fid = _field.identify(A, B, T, terms)
    if fid:
        m, f, pa, pb, pc = fid

        def to_tensor(dec):
            res = []
            for a, b, c in dec:
                res.append((sum(1 << pa[k] for k in range(m) if a >> k & 1),
                            sum(1 << pb[k] for k in range(m) if b >> k & 1),
                            sum(1 << pc[k] for k in range(m) if c >> k & 1)))
            return res
        table = _data("decomps.json")
        if _field.is_irreducible(f) and table:
            for ent in table["fields"]:
                if ent["m"] != m:
                    continue
                for k, d in enumerate(ent["decs"]):
                    try:
                        dd = _field.change_basis([tuple(t) for t in d["dec"]], m, d["f"], f)
                    except ValueError:
                        continue
                    out.append((to_tensor(dd), "table m=%d #%d (%s), rank %d" % (m, k, d["src"], len(dd)),
                                dict(m=m, f=f, pa=pa, pb=pb, pc=pc, field_dec=dd, key=k)))
        kd = _field.karatsuba_field(m, f)
        if not out or len(kd) < min(len(d) for d, _, _ in out):
            out.append((to_tensor(kd), "Karatsuba, rank %d" % len(kd),
                        dict(m=m, f=f, pa=pa, pb=pb, pc=pc, field_dec=kd, key=-1)))
    else:
        sd = _slice_decomposition(len(A), len(B), len(T), terms)
        budget = max(0.0, min(deadline - time.time(), 20.0)) * 0.5
        if budget > 0.2 and len(sd) > 1:
            sd = _flip.search(sd, time.time() + budget, seed=1)
        out.append((sd, "slice + flip graph, rank %d" % len(sd), None))
    good = [x for x in out if _anf.tensor_of(x[0]) == terms]
    return good


def _variant(dinfo, rng):
    """random isotopic variant (Frobenius, scalings) of a field decomposition, mapped to tensor indices."""
    m, f = dinfo["m"], dinfo["f"]
    k = rng.randrange(m)
    lam = rng.randrange(1, 1 << m)
    mu = rng.randrange(1, 1 << m)
    dd = _field.isotope(dinfo["field_dec"], m, f, k, lam, mu)
    pa, pb, pc = dinfo["pa"], dinfo["pb"], dinfo["pc"]
    res = [(sum(1 << pa[i] for i in range(m) if a >> i & 1), sum(1 << pb[i] for i in range(m) if b >> i & 1),
            sum(1 << pc[i] for i in range(m) if c >> i & 1)) for a, b, c in dd]
    return res, (k, lam, mu)


# ------------------------------------------------------------------------------------------------- search

def _run(pb, dec, p):
    """compile with parameter dict p -> qasm text."""
    if p["s"] == "frame":
        basis = None
        if not p.get("track"):
            G = [g for _, _, g in dec]
            key = (tuple(G), tuple(pb.M), pb.zero, p.get("bseed", 0))
            basis = _BASIS.get(key)
            if basis is None:
                if len(_BASIS) > 256:
                    _BASIS.clear()
                basis = _comp.out_basis(pb.M, G, random.Random(p.get("bseed", 0)), pb.zero,
                                        iters=p.get("biters", 400 if len(pb.T) <= 12 else 40))
                _BASIS[key] = basis
        return _comp.compile_frame(pb, dec, seed=p.get("seed", 0), order=p["order"], basis=basis,
                                   track=p.get("track", False))
    return _comp.compile_walk(pb, dec, order=p["order"], gain=p.get("gain", 2), seed=p.get("seed"),
                              fast=p.get("fast", True))


def _dec_for(decs, p):
    d, _, dinfo = decs[p["d"]]
    if p.get("iso") and dinfo:
        m, f = dinfo["m"], dinfo["f"]
        dd = _field.isotope(dinfo["field_dec"], m, f, *p["iso"])
        pa, pb_, pc = dinfo["pa"], dinfo["pb"], dinfo["pc"]
        d = [(sum(1 << pa[i] for i in range(m) if a >> i & 1), sum(1 << pb_[i] for i in range(m) if b >> i & 1),
              sum(1 << pc[i] for i in range(m) if c >> i & 1)) for a, b, c in dd]
    return d


def search(info, deadline, seed=1, recipes=(), keep=3):
    """compile many variants until `deadline`; -> list of (CX, T, qasm, params, note) best first."""
    pb = _problem(info)
    decs = decompositions(info, deadline)
    if not decs:
        return []
    rng = random.Random(seed)
    nA, nB = len(pb.A), len(pb.B)
    best = []
    seen = set()

    def consider(dec, p, note):
        txt = _run(pb, dec, p)
        T, CX = _check.counts(txt)
        if txt in seen:
            return
        seen.add(txt)
        best.append((CX, T, txt, dict(p), note))
        best.sort(key=lambda x: (x[1], x[0]))
        del best[keep:]

    min_rank = min(len(d) for d, _, _ in decs)
    useful = [i for i, (d, _, _) in enumerate(decs) if len(d) == min_rank]
    # 1) stored recipes (deterministic replays of the best runs found offline)
    for p in recipes:
        if p.get("d", 0) < len(decs) and time.time() < deadline:
            try:
                consider(_dec_for(decs, p), p, "recipe")
            except Exception:
                pass
    # 2) plain compilations of every decomposition
    for i in useful:
        d = decs[i][0]
        for p in (dict(s="frame", d=i, order=list(range(len(d))), seed=0),
                  dict(s="frame", d=i, order=list(range(len(d))), seed=0, track=True),
                  dict(s="walk", d=i, order=list(range(len(d))), gain=2)):
            if time.time() < deadline or not best:
                consider(d, p, decs[i][1])
    # 3) stage A: screen isotopic variants (Frobenius, scalings) of every decomposition with a fixed order
    t_start = time.time()
    stage_a = t_start + 0.35 * max(0.0, deadline - t_start)
    pool = []  # (cost, p)
    it = 0
    while time.time() < stage_a:
        it += 1
        i = useful[it % len(useful)]
        dinfo = decs[i][2]
        p = dict(d=i)
        if dinfo and it > len(useful):
            _, p["iso"] = _variant(dinfo, rng)
        d = _dec_for(decs, p)
        order = list(range(len(d)))
        qs = [dict(s="walk", gain=2), dict(s="frame", seed=0, track=True)]
        if it % 3 == 0:
            qs.append(dict(s="frame", seed=0, track=False, bseed=it))
        for q in qs:
            pp = dict(p, order=order, **q)
            try:
                txt = _run(pb, d, pp)
            except Exception:
                continue
            pool.append((_check.counts(txt)[1], pp))
        pool.sort(key=lambda x: x[0])
        del pool[12:]
    # 4) stage B: anneal the product order of the best variants with the true CX cost
    k = 0
    while time.time() < deadline and pool:
        _, p0 = pool[k % len(pool)]
        k += 1
        d = _dec_for(decs, p0)
        p = dict(p0)
        if p["s"] == "walk":
            p.update(gain=rng.choice((2, 3)), seed=rng.randrange(10 ** 6), fast=True)
        else:
            p["seed"] = rng.randrange(10 ** 6)
            if rng.random() < 0.3:
                p.update(track=rng.random() < 0.5, bseed=rng.randrange(10 ** 6))
        order = list(p["order"])
        if k > len(pool):
            rng.shuffle(order)
        p["order"] = order
        cur = _check.counts(_run(pb, d, p))[1]
        iters = 120 if len(d) <= 30 else 60
        temp0 = 3.0
        for step in range(iters):
            if time.time() > deadline:
                break
            o2 = list(p["order"])
            a, b = rng.sample(range(len(o2)), 2)
            if rng.random() < 0.5:
                o2[a], o2[b] = o2[b], o2[a]
            else:
                o2.insert(b, o2.pop(a))
            p2 = dict(p, order=o2)
            try:
                c2 = _check.counts(_run(pb, d, p2))[1]
            except Exception:
                continue
            tt = temp0 * (1 - step / iters) + 0.05
            if c2 <= cur or rng.random() < math.exp((cur - c2) / tt):
                p, cur = p2, c2
        try:
            consider(d, p, decs[p["d"]][1])
            if p["s"] == "walk":  # polish: slower register reduction
                consider(d, dict(p, fast=False), decs[p["d"]][1])
        except Exception:
            pass
    return best


# ------------------------------------------------------------------------------------------------- public API

def recipe_key(info):
    """m:f:model:sig, sig = 'z' (targets start at 0) or a digest of M written in the power basis."""
    if info.get("st") is None:
        return None
    A, B, T, terms = info["tensor"]
    fid = _field.identify(A, B, T, terms)
    if not fid:
        return None
    m, f, pa, pb, pc = fid
    if info["st"]["zero"]:
        sig = "z"
    else:
        import hashlib
        Tw = [T[pc[k]] for k in range(m)]
        rows = [sorted(Tw.index(u) for u in info["st"]["M"][t]) for t in Tw]
        sig = hashlib.sha1(json.dumps(rows).encode()).hexdigest()[:10]
    return "%d:%d:%s:%s" % (m, f, "clean" if info["clean"] else "strict", sig)


def _recipes_for(info):
    data = _data("recipes.json")
    key = recipe_key(info)
    if not data or key is None:
        return []
    return data.get("recipes", {}).get(key, [])


def _verify(qasm, info, qc_path, clean):
    ok, why = _check.structural_check(qasm, info["P"], info["n"], info["inp"], clean)
    if not ok:
        return False, why
    nq = int(qasm.split("qubit[", 1)[1].split("]", 1)[0])
    if nq <= NUMERIC_MAX_QUBITS:
        fid = _check.numeric_check(qasm, qc_path, clean, trials=2)
        if fid is None or fid < 1 - 1e-9:
            return False, "verify_qasm fidelity %r" % fid
        why += "; verify_qasm F=%.12f" % fid
    return True, why


def candidates(qc_path, clean, budget_sec=60.0):
    """verified measurement-based circuits for bilinear functions / known functions; [] if not recognised."""
    try:
        t0 = time.time()
        info = analyse(qc_path, clean)
        if info is None or info["ph"]:
            return []
        found = []
        # known constructions matched by function
        kn = _data("known.json")
        if kn:
            fp = _known.fingerprint(info["P"], info["n"], info["inp"], clean)
            for ent in kn.get("entries", []):
                if ent["fp"] == fp and bool(ent["clean"]) == bool(clean):
                    found.append((_known.decode(ent["nq"], ent["ops"]), "known construction: " + ent["src"]))
        if info["st"] is not None:
            deadline = t0 + 0.75 * budget_sec
            res = search(info, deadline, recipes=_recipes_for(info))
            for CX, T, txt, p, note in res:
                found.append((txt, "gf2m bilinear: %s; %d ANDs; strategy %s" % (note, T // 4, p.get("s"))))
        out = []
        for txt, note in found:
            try:
                ok, why = _verify(txt, info, qc_path, clean)
            except Exception as e:  # noqa: BLE001
                ok, why = False, str(e)
            if ok:
                out.append({"qasm": txt, "note": note + " | verified: " + why})
        out.sort(key=lambda d: _check.counts(d["qasm"]))
        return out
    except Exception:
        return []
