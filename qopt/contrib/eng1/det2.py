"""eng1/det2.py: DETERMINISTIC H-gadget compilation, protocol B (all ancillas alive, sequential measurement + feed-forward).

Path sum:  state = sum_v w^{f(v)} |W v>,  v=(x, a_1..a_h), wires after the body = coordinates y = W v  (n outputs L, h garbage Z_i).
Phase 1 (no classical control at all): ancillas |+>, phase polynomial f (T-count reduced by TODD, exact) + linear map to W.
Phase 2: measure garbage wires Z_i in X basis in an order sigma.  Outcome m_i == substitution a_i -> a_i xor m_i; the needed
correction is the FIXED Clifford polynomial  G_i(v) = f(v+e_i) - f(v) + 4 Z_i(v)  (applied if m_i = 1, in the frame of earlier
outcomes: every gate gets an extra frame correction when its wire parity contains measured a_j).  It is applicable iff G_i, written
in coordinates y, does not depend on garbage coordinates already measured (incl. Z_i):  edges i -> j for j in dep(G_i); sigma = a
topological order.  G_i depends only on the FUNCTION f, so TODD (function preserving) never breaks feasibility.
Finally outputs get X^{(L_r . eps)}.   Extra qubits = h.   T gates appear only in phase 1.
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, json, time, random, itertools
sys.path.insert(0, (_QPKG + "/contrib/eng1"))
from gad import *
import gad2, fasttodd, gidney
CORR = {"t": "sdg", "tdg": "s", "s": "z", "sdg": "z"}


def bits(m):
    return [k for k in range(m.bit_length()) if m >> k & 1]


def corrections(P, W, n, N):
    """-> list over i of canonical coefficient dict (keys: tuples of coordinate indices) of G_i in coordinates."""
    h = N - n; out = []
    for i in range(h):
        bit = 1 << (n + i); C = {}
        for M, w in P.items():
            if M & bit:
                Mp = pp.solve(W, N, M)
                for k, v in coeffs([(Mp, -2 * w)]).items(): C[k] = (C.get(k, 0) + v) % 8
        for k, v in coeffs([(1 << (n + i), 4)]).items(): C[k] = (C.get(k, 0) + v) % 8
        out.append({k: v for k, v in C.items() if v})
    return out


def topo_order(Cs, n, h):
    indeg = [0] * h; adj = [[] for _ in range(h)]
    for i, C in enumerate(Cs):
        d = {x - n for k in C for x in k if x >= n}
        if i in d: return None
        for j in d: adj[i].append(j); indeg[j] += 1
    q = [i for i in range(h) if indeg[i] == 0]; res = []
    while q:
        x = q.pop(0); res.append(x)
        for y in adj[x]:
            indeg[y] -= 1
            if indeg[y] == 0: q.append(y)
    return res if len(res) == h else None


def frag_pairs(i, C, W, n, measured, alive):
    """Correction block as S/CZ gates on coordinate wires (all conditional on c_i); frame corrections included."""
    out = []; c1 = 1 << i
    phi = lambda k: (W[k] >> n) & measured
    for key, v in sorted(C.items(), key=lambda kv: (len(kv[0]), kv[0])):
        assert all(x in alive for x in key), "correction touches a dead wire"
        if len(key) == 1:
            assert v % 2 == 0
            for g in pp.SEQ[v]:
                out.append(("g", g, (key[0],), (c1, 0)))
                if phi(key[0]) and g in CORR: out.append(("g", CORR[g], (key[0],), (c1, phi(key[0]))))
        elif len(key) == 2:
            assert v % 4 == 0, f"non-Clifford pair {v}"
            if v == 4:
                a_, b_ = key; out.append(("g", "cz", (a_, b_), (c1, 0)))
                if phi(b_): out.append(("g", "z", (a_,), (c1, phi(b_))))
                if phi(a_): out.append(("g", "z", (b_,), (c1, phi(a_))))
        else:
            assert v % 8 == 0, "triple coefficient in correction"
    return out


def frag_net(i, masks, W, n, measured, alive, rng, variants=6):
    """Correction block as a parity network (CX + conditional S/Z gates) over the alive coordinate wires, then restore."""
    al = sorted(alive); loc = {k: j for j, k in enumerate(al)}; nl = len(al); c1 = 1 << i
    lm = {}
    for m, w in masks.items():
        if not w % 8: continue
        if m >> n << n != m and False: pass
        mm = 0
        for k in bits(m):
            if k not in loc: return None
            mm |= 1 << loc[k]
        lm[mm] = (lm.get(mm, 0) + w) % 8
    lm = {m: w for m, w in lm.items() if w}
    best = None
    for v in range(variants):
        Wl = [1 << j for j in range(nl)]
        try:
            net = pp.graysynth(dict(lm), Wl, nl, rng) if v % 2 == 0 else pp.greedy(dict(lm), Wl, nl, rng, [1 << j for j in range(nl)], rng.uniform(0, 1.0))
        except Exception:
            net = None
        if net is None: continue
        M = list(Wl); ops = pp.lin_best(M, nl, range(nl), rng, 4)
        cost = sum(g == "cx" for g, _ in net) + len(ops)
        if best is None or cost < best[0]: best = (cost, net, ops)
    if best is None: return None
    _, net, ops = best
    cur = [1 << j for j in range(nl)]      # coordinate-mask (local) held by each local wire
    out = []
    def varmask(cm):
        M = 0
        for j in bits(cm): M ^= W[al[j]]
        return M
    for g, q in net:
        if g == "cx":
            c, t = q; cur[t] ^= cur[c]; out.append(("g", "cx", (al[c], al[t]), (c1, 0)))
        else:
            p = q[0]; out.append(("g", g, (al[p],), (c1, 0)))
            phi = (varmask(cur[p]) >> n) & measured
            if phi and g in CORR: out.append(("g", CORR[g], (al[p],), (c1, phi)))
    for c, t in ops: out.append(("g", "cx", (al[c], al[t]), (c1, 0)))
    return out


def gcount(frag): return sum(1 for o in frag if o[1] in ("cx", "cz"))


def build_prog(body, Cs, order, W, n, N, P=None, rng=None):
    h = N - n; prog = [("prep", q) for q in range(n, N)] + [("g", g, tuple(q), 0) for g, q in body]
    measured = 0; rng = rng or random.Random(5)
    for i in order:
        prog.append(("meas", n + i, i))
        alive = set(range(n)) | {n + j for j in range(h) if not (measured >> j & 1) and j != i}
        cands = [frag_pairs(i, Cs[i], W, n, measured, alive)]
        # canonical coefficients -> masks (coordinate basis)
        mk = {}
        for k, v in Cs[i].items():
            for m, w in gad_contrib(k, v): mk[m] = (mk.get(m, 0) + w) % 8
        fn = frag_net(i, mk, W, n, measured, alive, rng)
        if fn is not None: cands.append(fn)
        if P is not None:       # raw masks of the original polynomial (usually fewer than canonical ones)
            raw = {}
            bit = 1 << (n + i)
            for M, w in P.items():
                if M & bit:
                    Mp = pp.solve(W, N, M); raw[Mp] = (raw.get(Mp, 0) - 2 * w) % 8
            raw[1 << (n + i)] = (raw.get(1 << (n + i), 0) + 4) % 8
            fr = frag_net(i, raw, W, n, measured, alive, rng)
            if fr is not None and check_frag(fr, Cs[i], i): cands.append(fr)
        cands = [c for c in cands if check_frag(c, Cs[i], i)]
        frag = min(cands, key=lambda c: (gcount(c), len(c)))
        prog += frag
        measured |= 1 << i
    for r in range(n):
        c = W[r] >> n
        if c: prog.append(("g", "x", (r,), (c, 0)))
    return prog


def check_frag(frag, Ci, i):
    """Algebraic check: the (c_i,0)-conditioned gates of the fragment implement exactly G_i (canonical coefficients)."""
    wires = {}
    cur = {}
    C = {}
    Wc = {}
    # replay on coordinate masks
    state = {}
    def get(q): return state.get(q, 1 << q)
    for o in frag:
        if o[3] != (1 << i, 0): continue
        g, q = o[1], o[2]
        if g == "x": continue
        if g == "cx": state[q[1]] = get(q[1]) ^ get(q[0])
        elif g == "cz":
            a_, b_ = get(q[0]), get(q[1])
            for k, v in coeffs([(a_, 2), (b_, 2), (a_ ^ b_, -2)]).items(): C[k] = (C.get(k, 0) + v) % 8
        else:
            for k, v in coeffs([(get(q[0]), {"s": 2, "sdg": 6, "z": 4}[g])]).items(): C[k] = (C.get(k, 0) + v) % 8
    for q, m in state.items():
        if m != 1 << q: return False      # wires must be restored
    return {k: v for k, v in C.items() if v} == Ci


def counts(prog):
    T = sum(1 for o in prog if o[0] == "g" and o[1] in ("t", "tdg"))
    CX = sum(1 for o in prog if o[0] == "g" and o[1] in ("cx", "cz"))
    CXu = sum(1 for o in prog if o[0] == "g" and o[1] in ("cx", "cz") and not o[3])
    cond = sum(1 for o in prog if o[0] == "g" and o[3])
    return dict(T=T, CX=CX, CX_uncond=CXu, cond_cliffords=cond, measurements=sum(1 for o in prog if o[0] == "meas"))


def _cond_ok(cond, bits_):
    if not cond: return True
    c1, c2 = cond
    if sum(bits_.get(k, 0) for k in bits(c1)) % 2 == 0: return False
    if c2 and sum(bits_.get(k, 0) for k in bits(c2)) % 2 == 0: return False
    return True


def run_prog(prog, nq, psi, rng, forced=None):
    bits_ = {}
    for o in prog:
        if o[0] == "prep": psi = gidney.apply(psi, nq, "h", (o[1],))
        elif o[0] == "g":
            if _cond_ok(o[3], bits_): psi = gidney.apply(psi, nq, o[1], o[2])
        else:
            _, q, k = o
            psi = gidney.apply(psi, nq, "h", (q,))
            p = psi.reshape([2] * nq); idx = [slice(None)] * nq; idx[q] = 1
            p1 = float(np.sum(np.abs(p[tuple(idx)]) ** 2))
            m = forced[k] if forced is not None else int(rng.random() < p1)
            pr = p1 if m else 1 - p1
            if pr < 1e-12: return None
            idx0 = [slice(None)] * nq; idx0[q] = 1 - m
            p = p.copy(); p[tuple(idx0)] = 0; psi = p.reshape(-1) / np.sqrt(pr)
            if m: psi = gidney.apply(psi, nq, "x", (q,))
            bits_[k] = m
    return psi


def verify(prog, ref, n, h, trials=4, enum_max=7, seed=1, max_runs=200):
    """Random data states, ancillas |0> (prep makes |+>), REAL sampled outcomes (+ all outcomes for small h)."""
    nq = n + h; rng = np.random.default_rng(seed); worst = 1.0; runs = 0
    for t in range(trials):
        v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n); v /= np.linalg.norm(v)
        ex = v.copy()
        for g, q in ref: ex = gidney.apply(ex, n, g, tuple(q))
        tgt = np.zeros([2] * nq, complex); tgt[tuple([slice(None)] * n + [0] * h)] = ex.reshape([2] * n); tgt = tgt.reshape(-1)
        st = np.zeros([2] * nq, complex); st[tuple([slice(None)] * n + [0] * h)] = v.reshape([2] * n); st = st.reshape(-1)
        outs = [None] * 3
        if h <= enum_max and t == 0: outs = list(itertools.product([0, 1], repeat=h))
        for forced in outs:
            psi = run_prog(prog, nq, st.copy(), rng, forced)
            if psi is None: continue
            worst = min(worst, abs(np.vdot(tgt, psi))); runs += 1
    return bool(worst > 1 - 1e-10), float(worst), runs


def check_emission(prog, Cs, order, n, N):
    return all(check_frag([o for o in prog if o[0] == 'g' and o[3] and o[3][0] == (1 << i)], Cs[i], i) for i in range(N - n))


def check_emission_old(prog, Cs, order, n, N):
    """Algebraic check of the emitted correction code: gates conditioned on exactly (1<<i,0) re-extract to G_i."""
    for i in range(N - n):
        C = {}
        for o in prog:
            if o[0] != "g" or o[3] != (1 << i, 0): continue
            g, q = o[1], o[2]
            if g == "x": continue
            if g == "cz": add = coeffs([(1 << q[0], 2), (1 << q[1], 2), ((1 << q[0]) | (1 << q[1]), -2)])
            else: add = coeffs([(1 << q[0], {"s": 2, "sdg": 6, "z": 4}[g])])
            for k, v in add.items(): C[k] = (C.get(k, 0) + v) % 8
        if {k: v for k, v in C.items() if v} != Cs[i]: return False
    return True


def to_qasm3(prog, n, h, path=None):
    N = n + h
    L = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{N}] q;", f"bit[{max(h,1)}] c;"]
    ce = lambda m: " ^ ".join(f"c[{k}]" for k in bits(m))
    for o in prog:
        if o[0] == "prep": L.append(f"h q[{o[1]}];")
        elif o[0] == "meas": L += [f"h q[{o[1]}];", f"c[{o[2]}] = measure q[{o[1]}];"]
        else:
            s = f"{o[1]} " + ", ".join(f"q[{x}]" for x in o[2]) + ";"
            if o[3]:
                c1, c2 = o[3]; cond = f"({ce(c1)})" + (f" && ({ce(c2)})" if c2 else "")
                s = f"if ({cond}) {{ {s} }}"
            L.append(s)
    txt = "\n".join(L) + "\n"
    if path: Path(path).write_text(txt)
    return txt



def clifford_gates(C, place_masks=None):
    """canonical coefficient dict (singles even, pairs in {0,4}, no triples) -> gate list on wires = coordinate indices."""
    out = []
    for key, v in sorted(C.items(), key=lambda kv: (len(kv[0]), kv[0])):
        if len(key) == 1:
            assert v % 2 == 0
            out += [(g, (key[0],)) for g in pp.SEQ[v]]
        elif len(key) == 2:
            assert v % 4 == 0, "non-Clifford pair"
            if v == 4: out.append(("cz", key))
        else:
            assert v % 8 == 0, "triple in Clifford part"
    return out


def synth_split(P2, new, W, K, N, n, seeds, base_seed, place):
    """f = T-columns (weight 1) + Clifford remainder Q (exact).  Q as S/CZ gates at the start (variable basis) or at the end
    (coordinate basis), T-columns by phase-network synthesis.  Returns gate list (with cz) or None."""
    if place == "end" and any(K): return None
    PT = {m: 1 for m in new}
    Cq = {}
    for k, v in canon(P2).items(): Cq[k] = v
    for k, v in canon(PT).items(): Cq[k] = (Cq.get(k, 0) - v) % 8
    Cq = {k: v for k, v in Cq.items() if v}
    if place == "end":   # convert to coordinate basis
        C2 = {}
        masks = {}
        for k, v in Cq.items():
            for m, w in gad_contrib(k, v): masks[m] = (masks.get(m, 0) + w) % 8
        for m, w in masks.items():
            if w:
                for k2, v2 in coeffs([(pp.solve(W, N, m), w)]).items(): C2[k2] = (C2.get(k2, 0) + v2) % 8
        Cq = {k: v for k, v in C2.items() if v}
    try: qg = clifford_gates(Cq)
    except AssertionError: return None
    body = synth(PT, W, K, N, n, seeds=seeds, base_seed=base_seed) if PT else None
    if PT and body is None: return None
    if not PT:
        body = []   # no T: only the linear map remains
        rng = random.Random(0); M = [pp.solve(W, N, 1 << i) for i in range(N)]
        ops = pp.lin_best(list(M), N, range(N), rng, 4)
        body = [("cx", o) for o in ops] + [("x", (q,)) for q in range(n) if K[q]]
    return (qg + body) if place == "start" else (body + qg)


def gad_contrib(S, v):
    if len(S) == 1: return [(1 << S[0], v)]
    if len(S) == 2:
        a, b = S; h = v // 2
        return [(1 << a, h), (1 << b, h), ((1 << a) | (1 << b), -h)]
    return []


def det_gates(gates, hl, n, time_limit=20.0, restarts=2, tries=300, seeds=10, log=None, tag="", sim_max=22, use_eng2=False):
    t0 = time.time()
    N, P, W, K = build(gates, n); h = N - n
    assert rank_ok(W, N)
    Cs = corrections(P, W, n, N); order = topo_order(Cs, n, h)
    res = dict(n=n, h=h, N=N, start=tag)
    if order is None:
        res.update(feasible=False); return res, None
    odd, even = split_odd(P); res["t_poly"] = len(odd); best = None
    ref = evo_lib.lower_std(hl)
    for r in range(restarts):
        rng = random.Random(100 + r)
        new = fasttodd.todd_fast(odd, N, rng, time_limit, 300, tries) if odd else []
        P2 = clifford_rest(odd, new, even)
        variants = [("graysynth", synth(P2, W, K, N, n, seeds=seeds, base_seed=r * 100))]
        for place in ("start", "end"):
            variants.append(("split_" + place, synth_split(P2, new, W, K, N, n, seeds, r * 100, place)))
        for vn, body in variants:
            if body is None: continue
            c = (evo_lib.tcount(body), evo_lib.cx(body) + sum(g == "cz" for g, _ in body))
            if log: log(f"  [{tag}] r{r} {vn}: T={c[0]} CX(body)={c[1]} ({time.time()-t0:.0f}s)")
            if best is None or c < best[0]: best = (c, body)
    body = best[1]
    if use_eng2:
        sys.path.insert(0, (_QPKG + "/contrib/eng2"))
        import improve as e2
        b2, gained = e2.improve_gates(body, N, budget=2.0)
        if check_alg(b2, N, n, P, W, K) and evo_lib.tcount(b2) <= evo_lib.tcount(body) and evo_lib.cx(b2) < evo_lib.cx(body):
            body = b2
    assert check_alg(body, N, n, P, W, K), "body not equal to f"
    prog = build_prog(body, Cs, order, W, n, N, P, random.Random(5))
    res.update(counts(prog)); res.update(extra_qubits=h, feasible=True, order=order,
                                         alg_body=True, alg_corr=check_emission(prog, Cs, order, n, N))
    if N <= sim_max:
        ok, worst, runs = verify(prog, ref, n, h)
        res.update(verified=ok, min_fidelity=worst, sim_runs=runs)
    else:
        res.update(verified=None, min_fidelity=None, sim_runs=0)
    res["seconds"] = round(time.time() - t0, 1)
    return res, prog


def improve(path, qasm_out=None, ours_dir=(_QPKG + "/contrib/eng1/ours"), **kw):
    names, hl = qopt.load_qc(path); n = len(names); nm = Path(path).stem
    cands = []
    r, p = det_gates(prep_gates(hl), hl, n, tag="raw", **kw); cands.append((r, p))
    f = qopt.ours_file(path, Path(ours_dir) / f"{nm}.json")
    if f is not None:
        d = json.load(open(f)); best = [(g, tuple(q)) for g, q in d["best"]]
        r2, p2 = det_gates(best, hl, n, tag="ours", **kw); r2["ours_T"], r2["ours_CX"] = d["t"], d["cx"]; cands.append((r2, p2))
    cands = [c for c in cands if c[1] is not None and c[0]["verified"] is not False]
    if not cands: return dict(name=nm, feasible=False)
    cands.sort(key=lambda c: (c[0]["T"], c[0]["CX"]))
    res, prog = cands[0]; res["name"] = nm
    res["all"] = [(c[0]["start"], c[0]["T"], c[0]["CX"]) for c in cands]
    res["qasm"] = to_qasm3(prog, n, res["h"], qasm_out)
    return res


if __name__ == "__main__":
    for p in sys.argv[1:]:
        r = improve(p, log=lambda s: print(s, flush=True)); r.pop("qasm", None); print(json.dumps(r), flush=True)
