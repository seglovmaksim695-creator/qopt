"""tg3/core.py: post-optimisation of deterministic H-gadget circuits (eng1/det2 pipeline) by exploiting gauge freedoms.

Path sum:  state = sum_v w^{f(v)} |W v>,  v = (x, a_1..a_h).
Free (CX-cost-free) freedoms that keep the measured behaviour exact:
 (F1) initial ancilla wires are |+>:  any CX network R on (data, ancilla) wires with data rows fixed maps |x>|+..+> to itself, so the
      body may START from wires S0 = R*identity (ancilla rows = A a + B x), the physical circuit simply has no R.
 (F2) garbage wire i is X-measured: its measured content Y_i may be changed by adding any other rows (data rows give a
      conditional Z, free); corrections are re-derived in the new coordinates (det2.corrections / build_prog), feasibility =
      acyclicity of the dependency graph (det2.topo_order).
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, json, time, random, itertools, os
from pathlib import Path
for p in ((_QPKG + "/contrib/eng1"), (_QPKG + "/contrib/eng2"), _QPKG,
          (_QPKG + "/portfolio")):
    if p not in sys.path: sys.path.insert(0, p)
from gad import *          # pp, coeffs, build, split_odd, clifford_rest, canon, prep_gates, rank_ok, extract, ...
import det2, fasttodd, gidney
HERE = Path((_QPKG + "/contrib/tg3"))
(_QCACHE / "tg3").mkdir(parents=True, exist_ok=True)


def bits(m): return [k for k in range(m.bit_length()) if m >> k & 1]


def setup(path, start="raw", time_limit=20.0, tries=300, seed=100):
    names, hl = qoptct.load_qc(path); n = len(names); nm = Path(path).stem
    if start == "raw": gates = prep_gates(hl)
    else:
        _f = qoptct.ours_file(path, f"{_QPKG}/contrib/eng1/ours/{nm}.json")
        if _f is None: raise FileNotFoundError("no cached unitary start for " + nm)
        d = json.load(open(_f)); gates = [(g, tuple(q)) for g, q in d["best"]]
    N, P, W, K = build(gates, n); h = N - n
    assert rank_ok(W, N)
    odd, even = split_odd(P)
    cf = _QCACHE / "tg3" / f"{nm}_{start}_{seed}.json"
    if cf.exists(): new = json.load(open(cf))["new"]
    else:
        rng = random.Random(seed)
        new = fasttodd.todd_fast(odd, N, rng, time_limit, 300, tries) if odd else []
        json.dump(dict(new=new), open(cf, "w"))
    P2 = clifford_rest(odd, new, even)
    return dict(name=nm, path=path, n=n, N=N, h=h, P=P, P2=P2, W=W, K=K, hl=hl, start=start, new=new)


def mat_cx(rows, N):
    """CX list that turns identity wires into wires holding `rows` (masks) -- only for algebraic checks."""
    cur = [1 << i for i in range(N)]; ops = []
    # solve: find ops transforming cur->rows; do elimination on rows -> identity then reverse
    M = list(rows); inv = []
    for c in range(N):
        b = 1 << c
        if not M[c] & b:
            r = next(k for k in range(c + 1, N) if M[k] & b); M[c] ^= M[r]; inv.append((r, c))
        for k in range(N):
            if k != c and M[k] & b: M[k] ^= M[c]; inv.append((c, k))
    return [("cx", o) for o in reversed(inv)]


def synth2(ctx, Winit, Wend, seeds=8, base_seed=0, lam_set=(0.0, 0.3, 0.6, 1.0, 1.5), P2=None):
    N, n, K = ctx["N"], ctx["n"], ctx["K"]; P2 = P2 if P2 is not None else ctx["P2"]
    best = None
    for s in range(seeds):
        rng = random.Random(base_seed + s)
        W = list(Winit)
        try:
            if s % 3 == 0: out = pp.graysynth(dict(P2), W, N, rng)
            elif s % 3 == 1: out = pp.greedy(dict(P2), W, N, rng, Wend, lam_set[(s // 3) % len(lam_set)])
            else: out = pp.greedy(dict(P2), W, N, rng, Wend, rng.uniform(0, 1.5), rng.uniform(0, 1.5), rng.uniform(0.1, 1.0))
        except Exception:
            out = None
        if out is None: continue
        M = [pp.solve(Wend, N, W[r]) for r in range(N)]
        ops = pp.lin_best(M, N, range(N), rng, 6)
        out = out + [("cx", o) for o in ops] + [("x", (q,)) for q in range(n) if K[q]]
        out = evo_lib.fold_cancel(out, N)
        k = (evo_lib.tcount(out), evo_lib.cx(out))
        if best is None or k < best[0]: best = (k, out)
    return None if best is None else best[1]


def corr_info(ctx, Wend):
    n, N, h = ctx["n"], ctx["N"], ctx["h"]
    Cs = det2.corrections(ctx["P"], Wend, n, N); order = det2.topo_order(Cs, n, h)
    return Cs, order


def make_prog(ctx, body, Wend, Cs, order, seed=5):
    return det2.build_prog(body, Cs, order, Wend, ctx["n"], ctx["N"], ctx["P"], random.Random(seed))


def check_body(ctx, body, Winit, Wend):
    N, n = ctx["N"], ctx["n"]
    pre = mat_cx(Winit, N)
    c1, W1, K1 = extract(pre + list(body), N)
    return c1 == canon(ctx["P"]) and W1 == list(Wend) and K1[:n] == list(ctx["K"][:n])


def identity_init(ctx): return [1 << i for i in range(ctx["N"])]


def verify(ctx, prog, trials=4):
    N, n, h = ctx["N"], ctx["n"], ctx["h"]
    if N > 22: return None, None
    ref = evo_lib.lower_std(ctx["hl"])
    return det2.verify(prog, ref, n, h, trials=trials)[:2]


# ------------------------------------------------------------------ start from an existing deterministic qasm (same T!)
import re


def parse_det_qasm(path):
    """eng1/det2 qasm -> (N, prep qubits, body gates, rest lines).  body = unconditional gates after the |+> preps, before the first
    measurement section (which starts with an 'h')."""
    N = None; ops = []
    for L in open(path):
        L = L.strip()
        m = re.match(r"qubit\[(\d+)\]", L)
        if m: N = int(m.group(1)); continue
        if not L or L.startswith(("OPENQASM", "include", "bit")): continue
        ops.append(L.rstrip(";"))
    i = 0; prep = []
    while i < len(ops) and ops[i].startswith("h "):
        prep.append(int(re.findall(r"q\[(\d+)\]", ops[i])[0])); i += 1
    body = []
    while i < len(ops) and not ops[i].startswith("h ") and not ops[i].startswith(("if", "c[")):
        nm, args = ops[i].split(None, 1)
        body.append((nm, tuple(int(x) for x in re.findall(r"q\[(\d+)\]", args)))); i += 1
    return N, prep, body


def setup_qasm(qasm_path, qc_path):
    names, hl = qoptct.load_qc(qc_path); n = len(names); nm = Path(qc_path).stem
    N, prep, body = parse_det_qasm(qasm_path)
    assert prep == list(range(n, N)), (prep, n, N)
    W = [1 << i for i in range(N)]; K = [0] * N; P = {}
    pp.block_poly(body, W, K, P)
    assert rank_ok(W, N)
    P = {m: w for m, w in P.items() if w % 8}
    return dict(name=nm, path=qc_path, n=n, N=N, h=N - n, P=P, P2=P, W=W, K=K, hl=hl, start="qasm", body0=body,
                T=evo_lib.tcount(body), CX0=evo_lib.cx(body))


def cx_count_qasm(path):
    cx = 0; T = 0; meas = 0
    for L in open(path):
        L = L.strip()
        if "measure" in L: meas += 1
        for mm in re.finditer(r"\b(cx|cz) q", L): cx += 1
        for mm in re.finditer(r"\b(t|tdg) q", L): T += 1
    return T, cx, meas


# ---- free-CZ aware correction fragment choice (R-CZ in post.py): CZ touching a still-alive garbage wire is free
import post


def _gcount_free(frag, n=None):
    return sum(1 for o in frag if o[1] == "cx" or (o[1] == "cz" and not any(w >= _GN[0] for w in o[2])))


_GN = [0]


def set_n(n):
    _GN[0] = n


def final_prog(ctx, body, Wend, Cs, order):
    set_n(ctx["n"])
    prog = make_prog(ctx, body, Wend, Cs, order)
    return post.free_cz_pass(prog, ctx["n"])


def verify_file(qasm, ctx):
    import verify_qasm
    return verify_qasm.check(qasm, ctx["path"])


# ---------------------------------------------------------------- split variants (Clifford remainder as S/CZ at start / end)
def _canon_to_masks(Cq):
    masks = {}
    for k, v in Cq.items():
        for m, w in det2.gad_contrib(k, v): masks[m] = (masks.get(m, 0) + w) % 8
    return {m: w for m, w in masks.items() if w}


def synth_split2(ctx, Winit, Wend, place, seeds=4, base_seed=0):
    """T-columns by phase-network synthesis; Clifford remainder as S/Z/CZ gates at the start (physical variable basis)
    or at the end (coordinate basis)."""
    N, n, K, P2 = ctx["N"], ctx["n"], ctx["K"], ctx["P2"]
    if place == "end" and any(K): return None
    new = [m for m, w in P2.items() if w % 2]
    PT = {m: 1 for m in new}
    Cq = dict(canon(P2))
    for k, v in canon(PT).items(): Cq[k] = (Cq.get(k, 0) - v) % 8
    Cq = {k: v for k, v in Cq.items() if v}
    masks = _canon_to_masks(Cq)
    basis = Winit if place == "start" else Wend
    C2 = {}
    for m, w in masks.items():
        for k2, v2 in coeffs([(pp.solve(basis, N, m), w)]).items(): C2[k2] = (C2.get(k2, 0) + v2) % 8
    C2 = {k: v for k, v in C2.items() if v}
    try: qg = det2.clifford_gates(C2)
    except AssertionError: return None
    ctx2 = dict(ctx); ctx2["P2"] = PT
    if PT:
        body = synth2(ctx2, Winit, Wend, seeds=seeds, base_seed=base_seed)
        if body is None: return None
    else:
        body = []
    return (qg + body) if place == "start" else (body + qg)


def cost_body(b): return evo_lib.cx(b) + sum(1 for g, _ in b if g == "cz")


def synth3(ctx, Winit, Wend, seeds=8, base_seed=0, variants=("full", "start", "end")):
    best = None
    for v in variants:
        if v == "full": b = synth2(ctx, Winit, Wend, seeds=seeds, base_seed=base_seed)
        else: b = synth_split2(ctx, Winit, Wend, v, seeds=seeds, base_seed=base_seed)
        if b is None: continue
        k = (evo_lib.tcount(b), cost_body(b))
        if best is None or k < best[0]: best = (k, b, v)
    return best
