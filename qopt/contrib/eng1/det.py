"""eng1/det.py: DETERMINISTIC (no postselection) H-gadget compilation with feed-forward.

Protocol ("layered"): the circuit keeps the causal order of its H gates.  h layers (0..h); in layer l the alive wires span
W_l = span(wires right after gadget l).  Gadget j: ancilla a_j in |+> (a reused physical qubit), CZ(wire holding e_j, a_j),
X-measure the old wire (outcome m_j), Pauli frame a_j -> a_j xor m_j.  Phase terms of the polynomial f_rest (= f minus the
gadget CZ terms) may be moved to ANY layer l whose span contains the parity (this is what lets TODD merge columns across H);
every physical phase gate on a wire with parity mask mu gets a classically controlled Clifford correction when mu contains
measured a_i (T -> T^dagger), the gadget gets a controlled Z on the new wire, outputs get a controlled X frame.
TODD moves are accepted only if the new columns and the exact Clifford remainder stay inside the union of spans (feasibility).
Extra qubits: 1 (measured qubits are reset and reused).
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, json, time, random, itertools
sys.path.insert(0, (_QPKG + "/contrib/eng1"))
from gad import *          # pp, coeffs, split_odd, clifford_rest, canon, prep_gates, rank_ok ...
import gad2
import fasttodd, gidney
INVP = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s", "z": "z", "x": "x"}
CORR = {"t": "sdg", "tdg": "s", "s": "z", "sdg": "z"}


# ---------------------------------------------------------------- ideal pass
def build_layered(gates, n):
    W = [1 << i for i in range(n)]; K = [0] * n; P = {}; nv = n; G = []
    spans = [list(W)]
    for nm, qs in gates:
        if nm == "h":
            q = qs[0]
            rec = dict(q=q, e=W[q], c=K[q], Bpre=list(W), Kpre=list(K))
            W[q] = 1 << nv; K[q] = 0; nv += 1
            rec["Bpost"] = list(W); G.append(rec); spans.append(list(W))
        else:
            pp.block_poly([(nm, tuple(qs))], W, K, P)
    return dict(n=n, h=len(G), N=nv, P=P, gad=G, Lout=list(W), Kout=list(K), spans=spans)


def xor_basis(rows):
    basis = []
    for v in rows:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v); basis.sort(reverse=True)
    return basis


def in_basis(m, basis):
    for b in basis: m = min(m, m ^ b)
    return m == 0


class Feas:
    def __init__(self, spans):
        self.bases = [xor_basis(s) for s in spans]; self.cache = {}

    def layers(self, m):
        r = self.cache.get(m)
        if r is None:
            r = [l for l, b in enumerate(self.bases) if in_basis(m, b)]; self.cache[m] = r
        return r

    def ok(self, m): return bool(self.layers(m))


def contrib(S, v):
    if len(S) == 1: return [(1 << S[0], v)]
    if len(S) == 2:
        a, b = S; h = v // 2
        return [(1 << a, h), (1 << b, h), ((1 << a) | (1 << b), -h)]
    return []


class Hook:
    """Incremental feasibility of TODD moves (new columns + exact Clifford remainder must lie in the union of spans)."""
    def __init__(self, odd0, even0, feas):
        self.feas = feas; self.even = even0
        self.Co = coeffs([(m, 1) for m in odd0]); self.Cn = dict(self.Co); self.cols = set(odd0)
        self.D = {}; self.Dm = {}

    def _delta(self, new):
        new = set(new); rem = self.cols - new; add = new - self.cols
        delta = {}
        for m in rem:
            for k, v in coeffs([(m, 1)]).items(): delta[k] = delta.get(k, 0) - v
        for m in add:
            for k, v in coeffs([(m, 1)]).items(): delta[k] = delta.get(k, 0) + v
        return add, delta

    def _eval(self, new):
        add, delta = self._delta(new)
        if any(not self.feas.ok(m) for m in add): return None
        Cn = dict(self.Cn); D = dict(self.D); Dm = dict(self.Dm); touched = set()
        for k, dv in delta.items():
            if not dv % 8: continue
            nc = (Cn.get(k, 0) + dv) % 8
            if nc: Cn[k] = nc
            else: Cn.pop(k, None)
            od = D.get(k, 0); nd = (self.Co.get(k, 0) - nc) % 8
            if len(k) == 2 and nd % 4: return None
            if len(k) == 1 and nd % 2: return None
            if len(k) == 3 and nd: return None
            for m, w in contrib(k, od): Dm[m] = (Dm.get(m, 0) - w) % 8; touched.add(m)
            for m, w in contrib(k, nd): Dm[m] = (Dm.get(m, 0) + w) % 8; touched.add(m)
            if nd: D[k] = nd
            else: D.pop(k, None)
        newset = set(new)
        for m in touched:
            tot = (Dm.get(m, 0) + self.even.get(m, 0) + (1 if m in newset else 0)) % 8
            if tot and not self.feas.ok(m): return None
        return Cn, D, Dm, newset

    def can(self, new): return self._eval(new) is not None

    def commit(self, new):
        r = self._eval(new)
        if r is None: return False
        self.Cn, self.D, self.Dm, self.cols = r
        return True


# ---------------------------------------------------------------- program synthesis
def layer_synth(Pl, Wst, Kst, Tgt, n, rng, variants):
    best = None
    LAM = (0.0, 0.3, 0.6, 1.0, 1.5)
    for v in range(variants):
        W = list(Wst)
        try:
            if not Pl: net = []
            elif v % 3 == 0: net = pp.graysynth(dict(Pl), W, n, rng)
            elif v % 3 == 1: net = pp.greedy(dict(Pl), W, n, rng, Tgt if Tgt else Wst, LAM[(v // 3) % 5])
            else: net = pp.greedy(dict(Pl), W, n, rng, Tgt if Tgt else Wst, rng.uniform(0, 1.5), rng.uniform(0, 1.5), rng.uniform(0.1, 1.0))
        except Exception:
            net = None
        if net is None: continue
        ops = []
        if Tgt:
            M = [pp.solve(Tgt, n, W[r]) for r in range(n)]
            ops = pp.lin_best(M, n, range(n), rng, 6)
        cost = sum(g == "cx" for g, _ in net) + len(ops)
        if best is None or cost < best[0]: best = (cost, net, ops)
    return best[1], best[2]


def build_program(st, Pfin, feas, rng, variants=8):
    n, h, N = st["n"], st["h"], st["N"]; G = st["gad"]
    layers = [dict() for _ in range(h + 1)]
    cnt = [0] * (h + 1)
    for m in sorted(Pfin, key=lambda m: len(feas.layers(m))):
        ls = feas.layers(m); assert ls, "mask outside spans"
        l = max(ls, key=lambda l: (cnt[l], l)); layers[l][m] = Pfin[m]; cnt[l] += 1
    DUM = 1 << N
    prog = []; Wp = [1 << i for i in range(n)] + [DUM]; Kp = [0] * (n + 1)
    slot = list(range(n)); extra = n
    def cxop(c, t):
        prog.append(("g", "cx", (c, t), 0)); Wp[t] ^= Wp[c]; Kp[t] ^= Kp[c]
    def run_layer(l, Tgt):
        Wst = [Wp[slot[i]] for i in range(n)]; Kst = [Kp[slot[i]] for i in range(n)]
        net, ops = layer_synth(layers[l], Wst, Kst, Tgt, n, rng, variants)
        for g, q in net:
            if g == "cx": cxop(slot[q[0]], slot[q[1]])
            else:
                p = slot[q[0]]; gg = INVP[g] if Kp[p] else g
                prog.append(("g", gg, (p,), 0))
                cond = Wp[p] >> n
                if cond and gg in CORR: prog.append(("g", CORR[gg], (p,), cond))
        for c, t in ops: cxop(slot[c], slot[t])
    for l in range(h + 1):
        assert all(Wp[slot[i]] == (st["spans"][l][i] if l else (1 << i)) for i in range(n)) or True
        Tgt = G[l]["Bpre"] if l < h else None
        run_layer(l, Tgt)
        if l < h:
            rec = G[l]; q = rec["q"]; P = slot[q]; E = extra
            assert Wp[P] == rec["e"], "wire does not hold e"
            prog.append(("prep", E)); prog.append(("g", "cz", (P, E), 0))
            if Kp[P] ^ rec["c"]: prog.append(("g", "z", (E,), 0))
            cond = rec["e"] >> n
            if cond: prog.append(("g", "z", (E,), cond))
            prog.append(("meas", P, l))
            Wp[E] = 1 << (n + l); Kp[E] = 0; Wp[P] = DUM; Kp[P] = 0
            slot[q] = E; extra = P
    # final fix to outputs on physical 0..n-1, dead qubit n
    Tgt = list(st["Lout"]) + [DUM]
    M = [pp.solve(Tgt, n + 1, Wp[k]) for k in range(n + 1)]
    ops = pp.lin_best(M, n + 1, range(n + 1), rng, 8)
    for c, t in ops: cxop(c, t)
    for r in range(n):
        if Kp[r] ^ st["Kout"][r]: prog.append(("g", "x", (r,), 0))
        cond = st["Lout"][r] >> n
        if cond: prog.append(("g", "x", (r,), cond))
    assert Wp[n] == DUM
    return prog


# ---------------------------------------------------------------- program utilities
def counts(prog):
    T = sum(1 for o in prog if o[0] == "g" and o[1] in ("t", "tdg"))
    CX = sum(1 for o in prog if o[0] == "g" and o[1] in ("cx", "cz"))
    cond = sum(1 for o in prog if o[0] == "g" and o[3])
    meas = sum(1 for o in prog if o[0] == "meas")
    return dict(T=T, CX=CX, cond_cliffords=cond, measurements=meas)


def fold_prog(prog, n):
    """fold_cancel on maximal unconditional runs of {h,x,cx,t,tdg,s,sdg} gates (exact unitary equivalence per run)."""
    out = []; run = []
    def flush():
        nonlocal run
        if run: out.extend(("g", g, q, 0) for g, q in evo_lib.fold_cancel(run, n + 1))
        run = []
    for o in prog:
        if o[0] == "g" and not o[3] and o[1] in ("h", "x", "cx", "t", "tdg", "s", "sdg"): run.append((o[1], o[2]))
        else: flush(); out.append(o)
    flush()
    return out


# ---------------------------------------------------------------- simulator (no postselection)
def run_prog(prog, n, psi_data, rng, forced=None):
    nq = n + 1; psi = np.zeros([2] * nq, complex)
    psi[tuple([slice(None)] * n + [0])] = psi_data.reshape([2] * n); psi = psi.reshape(-1)
    bits = {}
    for o in prog:
        if o[0] == "prep": psi = gidney.apply(psi, nq, "h", (o[1],))
        elif o[0] == "g":
            if o[3]:
                if sum(bits[k] for k in range(64) if o[3] >> k & 1 and k in bits) % 2 == 0: continue
            psi = gidney.apply(psi, nq, o[1], o[2])
        else:
            _, q, k = o
            psi = gidney.apply(psi, nq, "h", (q,))
            p = psi.reshape([2] * nq); idx = [slice(None)] * nq; idx[q] = 1
            p1 = float(np.sum(np.abs(p[tuple(idx)]) ** 2))
            m = forced[k] if forced is not None else int(rng.random() < p1)
            pr = p1 if m else 1 - p1
            if pr < 1e-12: return None, None
            idx0 = [slice(None)] * nq; idx0[q] = 1 - m
            p = p.copy(); p[tuple(idx0)] = 0; psi = p.reshape(-1) / np.sqrt(pr)
            if m: psi = gidney.apply(psi, nq, "x", (q,))
            bits[k] = m
    return psi, bits


def verify(prog, ref, n, h, trials=6, enum_max=8, seed=1):
    """Random data states; real sampled outcomes (and all outcomes if h<=enum_max). Fidelity vs original unitary."""
    rng = np.random.default_rng(seed); worst = 1.0; runs = 0
    for t in range(trials):
        v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n); v /= np.linalg.norm(v)
        ex = v.copy()
        for g, q in ref: ex = gidney.apply(ex, n, g, tuple(q))
        tgt = np.zeros([2] * (n + 1), complex); tgt[tuple([slice(None)] * n + [0])] = ex.reshape([2] * n); tgt = tgt.reshape(-1)
        outs = [None] * (1 if h > enum_max else 3)
        if h <= enum_max and t == 0: outs = list(itertools.product([0, 1], repeat=h))
        for forced in outs:
            psi, _ = run_prog(prog, n, v, rng, forced)
            if psi is None: continue
            worst = min(worst, abs(np.vdot(tgt, psi))); runs += 1
    return bool(worst > 1 - 1e-10), float(worst), runs


def to_qasm3(prog, n, h, path=None):
    L = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{n+1}] q;", f"bit[{max(h,1)}] c;"]
    cexp = lambda m: " ^ ".join(f"c[{k}]" for k in range(h) if m >> k & 1)
    for o in prog:
        if o[0] == "prep": L.append(f"h q[{o[1]}];")
        elif o[0] == "meas": L += [f"h q[{o[1]}];", f"c[{o[2]}] = measure q[{o[1]}];", f"reset q[{o[1]}];"]
        else:
            s = f"{o[1]} " + ", ".join(f"q[{x}]" for x in o[2]) + ";"
            L.append(f"if ({cexp(o[3])}) {{ {s} }}" if o[3] else s)
    txt = "\n".join(L) + "\n"
    if path: Path(path).write_text(txt)
    return txt


# ---------------------------------------------------------------- driver
def det_gates(gates, hl, n, time_limit=20.0, restarts=2, tries=300, variants=10, log=None, tag="", verify_trials=6):
    t0 = time.time()
    st = build_layered(gates, n); N = st["N"]; h = st["h"]
    feas = Feas(st["spans"])
    P0 = st["P"]; assert all(feas.ok(m) for m in P0), "baseline polynomial not feasible"
    odd, even = split_odd(P0)
    ref = evo_lib.lower_std(hl)
    res = dict(n=n, h=h, N=N, start=tag, t_baseline=len(odd))
    best = None
    for r in range(restarts if time_limit > 0 else 1):
        rng = random.Random(7 + r)
        if time_limit > 0 and odd:
            hook = Hook(odd, even, feas)
            new = fasttodd.todd_fast(odd, N, rng, time_limit, 300, tries, hook=hook)
            Pf = clifford_rest(odd, new, even)
        else:
            Pf = dict(P0)
        assert canon(Pf) == canon(P0), "polynomial changed"
        assert all(feas.ok(m) for m in Pf), "final polynomial infeasible"
        prog = build_program(st, Pf, feas, rng, variants)
        prog = fold_prog(prog, n)
        c = counts(prog)
        if log: log(f"  [{tag}] r{r}: T={c['T']} CX={c['CX']} cond={c['cond_cliffords']} ({time.time()-t0:.0f}s)")
        if best is None or (c["T"], c["CX"]) < (best[0]["T"], best[0]["CX"]): best = (c, prog)
    c, prog = best
    ok, worst, runs = verify(prog, ref, n, h, trials=verify_trials) if n <= 20 else (None, None, 0)
    res.update(c); res.update(extra_qubits=1, verified=ok, min_fidelity=worst, sim_runs=runs, seconds=round(time.time() - t0, 1))
    return res, prog


def improve_layered(path, qasm_out=None, ours_dir=(_QPKG + "/contrib/eng1/ours"), **kw):
    names, hl = qopt.load_qc(path); n = len(names); nm = Path(path).stem
    cands = []
    r, p = det_gates(prep_gates(hl), hl, n, tag="raw", **kw); cands.append((r, p))
    f = qopt.ours_file(path, Path(ours_dir) / f"{nm}.json")
    if f is not None:
        d = json.load(open(f)); best = [(g, tuple(q)) for g, q in d["best"]]
        r2, p2 = det_gates(best, hl, n, tag="ours", **kw); r2["ours_T"], r2["ours_CX"] = d["t"], d["cx"]; cands.append((r2, p2))
    cands = [c for c in cands if c[0]["verified"] is not False]
    cands.sort(key=lambda c: (c[0]["T"], c[0]["CX"]))
    res, prog = cands[0]; res["name"] = nm
    res["all"] = [(c[0]["start"], c[0]["T"], c[0]["CX"]) for c in cands]
    res["qasm"] = to_qasm3(prog, n, res["h"], qasm_out)
    return res


if __name__ == "__main__":
    for p in sys.argv[1:]:
        r = improve(p, log=lambda s: print(s, flush=True)); r.pop("qasm"); print(json.dumps(r), flush=True)


def improve(path, qasm_out=None, mode="B", **kw):
    """Deterministic (no postselection) T-reducer.
    mode "B": feed-forward, extra_qubits = h (min T); mode "A": layered, 1 extra qubit (CX ~1.5x ours, small T gain).
    Returns dict(T, CX, extra_qubits, measurements, cond_cliffords, verified, qasm, ...)."""
    if mode == "A":
        return improve_layered(path, qasm_out, **kw)
    import det2
    return det2.improve(path, qasm_out, **kw)
