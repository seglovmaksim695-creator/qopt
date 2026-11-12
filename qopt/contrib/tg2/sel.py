"""tg2/sel.py: SELECTIVE Hadamard gadgetization + windowed TODD + qubit reuse.

Idea: keep a subset K of the internal H gates as real H gates.  They cut the circuit into windows (segments); inside a window all
H gates are gadgetized, the window's phase polynomial is TODD-reduced and the window is compiled deterministically with
eng1's protocol B (det2.corrections / topo_order / build_prog).  Windows are independent, so
    T = sum T_w,  CX = sum CX_w,  extra qubits = max_w h_w   (measured qubits are reset by a conditional X and reused).
Window costs are cached for all (a,b); a Pareto DP over the choice of kept H's (and of the candidate per window) gives the
(T, CX, extra) front.  Everything is exported as OpenQASM 3 in the subset accepted by verify_qasm.py (AND-only conditions;
XOR conditions are expanded into sequential ifs).
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, os, json, time, random, pickle, subprocess, itertools
from pathlib import Path
sys.path.insert(0, (_QPKG + "/contrib/eng1"))
import det2, det, gad, fasttodd
from gad import build, rank_ok, split_odd, clifford_rest, check_alg, synth, prep_gates, qopt, evo_lib, pp
ROOT = Path(_QPKG)
HERE = Path(__file__).resolve().parent
OURS = Path((_QPKG + "/contrib/eng1/ours"))


def prune(items, key):
    """keep nondominated items (key -> tuple, all minimised)."""
    items = sorted(items, key=key); out = []
    for it in items:
        k = key(it)
        if not any(all(a <= b for a, b in zip(key(o), k)) for o in out): out.append(it)
    return out


# ------------------------------------------------------------------ window candidates
def seg_cands(gates, n, tl=2.0, restarts=2, seeds=6, tries=200, seed0=0):
    """Trailing X/Z Paulis of the window are peeled off and appended unconditionally (det2's protocol B mishandles output
    X constants) ."""
    gates = list(gates); trail = []
    while gates and gates[-1][0] in ("x", "z"): trail.append(gates.pop())
    res = _seg_cands(gates, n, tl, restarts, seeds, tries, seed0)
    if res and trail:
        res = [(T, CX, e, prog + [("g", g, tuple(q), 0) for g, q in reversed(trail)], pr, nb) for (T, CX, e, prog, pr, nb) in res]
    return res


def _seg_cands(gates, n, tl=2.0, restarts=2, seeds=6, tries=200, seed0=0):
    """-> list of (T, CX, h, prog) nondominated, or None if protocol B is infeasible for this window."""
    N, P, W, K = build(gates, n); h = N - n
    if not rank_ok(W, N): return None
    Cs = det2.corrections(P, W, n, N) if h else []
    order = det2.topo_order(Cs, n, h) if h else []
    if order is None: return None
    odd, even = split_odd(P); out = []
    for r in range(restarts + 1):
        rng = random.Random(1000 + seed0 + r)
        new = list(odd) if r == 0 else (fasttodd.todd_fast(odd, N, rng, tl, 300, tries) if odd else [])
        P2 = clifford_rest(odd, new, even)
        variants = [synth(P2, W, K, N, n, seeds=seeds, base_seed=r * 100)]
        for place in ("start", "end"):
            try: variants.append(det2.synth_split(P2, new, W, K, N, n, seeds, r * 100, place))
            except Exception: pass
        for body in variants:
            if body is None or not check_alg(body, N, n, P, W, K): continue
            prog = det2.build_prog(body, Cs, order, W, n, N, P, random.Random(5))
            c = det2.counts(prog)
            out.append((c["T"], c["CX"], h, prog, "B", h))
    if h > 0:
        out += seg_cands_A(gates, n, tl, restarts, tries, seed0)
    return prune(out, lambda x: (x[0], x[1], x[2])) if out else None


def seg_cands_A(gates, n, tl, restarts, tries, seed0):
    """layered protocol (eng1 det.py): 1 extra qubit, windows with feasible TODD."""
    out = []
    try:
        st = det.build_layered(gates, n); N = st["N"]; feas = det.Feas(st["spans"])
        P0 = st["P"]; odd, even = split_odd(P0)
        for r in range(restarts + 1):
            rng = random.Random(2000 + seed0 + r)
            if r == 0: Pf = dict(P0)
            else:
                hook = det.Hook(odd, even, feas)
                new = fasttodd.todd_fast(odd, N, rng, tl, 300, tries, hook=hook) if odd else []
                Pf = clifford_rest(odd, new, even)
            prog = det.build_program(st, Pf, feas, rng, 8)
            prog = det.fold_prog(prog, n)
            c = det.counts(prog)
            out.append((c["T"], c["CX"], 1, prog, "A", st["h"]))
    except Exception as e:
        print("A-fail", repr(e)[:80], flush=True)
    return out


def windows(gates):
    hpos = [i for i, (g, q) in enumerate(gates) if g == "h"]
    return hpos


def push_paulis(gs):
    """Remove X gates by propagating Pauli X/Z to the end of the window (det2's protocol B is unsound when a gadget wire carries an
    X constant).  Exact up to global phase.  Windows with ccz are left untouched."""
    if not any(g == "x" for g, q in gs) or any(g == "ccz" for g, q in gs): return gs
    xm = {}; zm = {}; out = []
    g1 = lambda d, q: d.get(q, 0)
    flip = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s"}
    for g, q in gs:
        if g == "x": xm[q[0]] = g1(xm, q[0]) ^ 1; continue
        if g == "cx":
            c, t = q; xm[t] = g1(xm, t) ^ g1(xm, c); zm[c] = g1(zm, c) ^ g1(zm, t); out.append((g, q))
        elif g == "cz":
            a, b = q; za, zb = g1(zm, a) ^ g1(xm, b), g1(zm, b) ^ g1(xm, a); zm[a] = za; zm[b] = zb; out.append((g, q))
        elif g == "h":
            x_, z_ = g1(xm, q[0]), g1(zm, q[0]); xm[q[0]] = z_; zm[q[0]] = x_; out.append((g, q))
        elif g in flip: out.append((flip[g] if g1(xm, q[0]) else g, q))
        else: out.append((g, q))
    for q in sorted(set(xm) | set(zm)):
        if g1(xm, q): out.append(("x", (q,)))
        if g1(zm, q): out.append(("z", (q,)))
    return out


def seg_gates(gates, hpos, a, b):
    s = hpos[a] + 1 if a >= 0 else 0
    e = hpos[b] if b < len(hpos) else len(gates)
    return push_paulis([(g, tuple(q)) for g, q in gates[s:e]])


def build_cache(gates, n, Hmax, cache_path, tl=2.0, restarts=2, log=print, deadline=None):
    hpos = windows(gates); h = len(hpos)
    cache = pickle.load(open(cache_path, "rb")) if Path(cache_path).exists() else {}
    if not cache.get("__v3__"):                  # purge windows that contain X gates (old, unsound B candidates)
        cache.pop("__v2__", None)
        for k in [k for k in cache if isinstance(k, tuple)]:
            a, b = k; s_ = hpos[a] + 1 if a >= 0 else 0; e_ = hpos[b] if b < h else len(gates)
            if any(g == "x" for g, q in gates[s_:e_]): del cache[k]
        cache["__v3__"] = True
    t0 = time.time()
    for L in range(0, Hmax + 1):                 # window gadget count (increasing, so deadline keeps small windows)
        for a in range(-1, h):
            b = a + L + 1
            if b > h or (a, b) in cache: continue
            if deadline and time.time() > deadline: log("deadline"); break
            c = seg_cands(seg_gates(gates, hpos, a, b), n, tl, restarts)
            cache[(a, b)] = c
            pickle.dump(cache, open(cache_path, "wb"))
        log(f"  L={L} done {time.time()-t0:.0f}s; feasible windows {sum(1 for k, v in cache.items() if isinstance(k, tuple) and v)}/{sum(1 for k in cache if isinstance(k, tuple))}")
    return cache


# ------------------------------------------------------------------ DP
def pareto_dp(h, cache, Hmax, maxfront=150):
    front = {-1: [(0, 0, 0, ())]}
    for b in range(0, h + 1):
        cur = []
        for a in range(max(-1, b - 1 - Hmax), b):
            cs = cache.get((a, b))
            if not cs or a not in front: continue
            for ci, (T, CX, hh, _, _p, _nb) in enumerate(cs):
                for (T0, C0, h0, path) in front[a]:
                    cur.append((T0 + T, C0 + CX, max(h0, hh), path + ((a, b, ci),)))
        cur = prune(cur, lambda x: (x[0], x[1], x[2]))
        if len(cur) > maxfront:
            cur = sorted(cur, key=lambda x: (x[0], x[1]))[:: max(1, len(cur) // maxfront)]
        front[b] = cur
    return front.get(h, [])


# ------------------------------------------------------------------ assembly
def assemble(gates, n, cache, path):
    hpos = windows(gates); h = len(hpos); ops = []; off = 0; hmax = 0
    for (a, b, ci) in path:
        T, CX, hh, prog, proto, nb = cache[(a, b)][ci]
        for o in prog:
            if o[0] == "g" and o[3] and not isinstance(o[3], tuple): o = ("g", o[1], o[2], (o[3], 0))
            if o[0] == "prep": ops.append(o)
            elif o[0] == "meas":
                ops.append(("meas", o[1], o[2])); ops.append(("g", "x", (o[1],), (1 << o[2], 0)))
            else:
                c = o[3]; c = (c[0], c[1]) if c else 0
                ops.append(("g", o[1], o[2], c))
        off = max(off, nb); hmax = max(hmax, hh)
        if b < h: ops.append(("g", "h", (gates[hpos[b]][1][0],), 0))
    return ops, off, hmax


def fold_ops(ops, nq):
    out = []; run = []
    def flush():
        nonlocal run
        if run: out.extend(("g", g, q, 0) for g, q in evo_lib.fold_cancel(run, nq))
        run = []
    for o in ops:
        if o[0] == "g" and not o[3] and o[1] in ("h", "x", "cx", "t", "tdg", "s", "sdg"): run.append((o[1], o[2]))
        else: flush(); out.append(o)
    flush()
    return out


def counts(ops):
    return dict(T=sum(1 for o in ops if o[0] == "g" and o[1] in ("t", "tdg")),
                CX=sum(1 for o in ops if o[0] == "g" and o[1] in ("cx", "cz")),
                measurements=sum(1 for o in ops if o[0] == "meas"),
                cond_gates=sum(1 for o in ops if o[0] == "g" and o[3]))


def emit_qasm(ops, nq, nbits, path=None):
    L = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"bit[{max(nbits,1)}] c;", f"qubit[{nq}] q;"]
    bits = lambda m: [k for k in range(m.bit_length()) if m >> k & 1]
    def ifg(bs, g, qs):
        s = f"{g} " + ", ".join(f"q[{x}]" for x in qs) + ";"
        if not bs: return s
        return "if (" + " && ".join(f"(c[{k}])" for k in sorted(bs)) + ") { " + s + " }"
    for o in ops:
        if o[0] == "prep": L.append(f"h q[{o[1]}];")
        elif o[0] == "meas": L += [f"h q[{o[1]}];", f"c[{o[2]}] = measure q[{o[1]}];"]
        else:
            g, qs, cond = o[1], o[2], o[3]
            if not cond: L.append(ifg([], g, qs)); continue
            b1 = bits(cond[0]); b2 = bits(cond[1]) if cond[1] else None
            terms = [{a, b} for a in b1 for b in b2] if b2 else [{a} for a in b1]
            if g in ("x", "z", "cz", "cx", "h"):
                for t in terms: L.append(ifg(t, g, qs))
            elif g in ("s", "sdg"):
                for t in terms: L.append(ifg(t, g, qs))
                for t1, t2 in itertools.combinations(terms, 2): L.append(ifg(t1 | t2, "z", qs))
            else: raise ValueError("conditional " + g)
    txt = "\n".join(L) + "\n"
    if path: Path(path).write_text(txt)
    return txt


def verify(qasm_path, qc_path, trials=30):
    r = subprocess.run([sys.executable, str(ROOT / "verify_qasm.py"), str(qasm_path), str(qc_path)], cwd=ROOT,
                       capture_output=True, text=True, timeout=3000)
    try: return float(r.stdout.strip().split()[-1])
    except Exception: return None


# ------------------------------------------------------------------ driver
def get_start(name, start):
    names, hl = qopt.load_qc(str(ROOT / "bench_qc" / f"{name}.qc")); n = len(names)
    if start == "ours":
        d = json.load(open(OURS / f"{name}.json")); gates = [(g, tuple(q)) for g, q in d["best"]]
        ref = (d["t"], d["cx"])
    else:
        gates = prep_gates(hl); ref = None
    return n, gates, ref


def pareto_points(name, start, Hmax_list=(0, 1, 2, 3, 4, 6, 8, 10), cache=None):
    n, gates, ref = get_start(name, start)
    h = len(windows(gates)); pts = []
    for Hm in Hmax_list:
        for (T, CX, hm, path) in pareto_dp(h, cache, Hm):
            pts.append(dict(T=T, CX=CX, extra=hm, path=path, start=start, kept=h - sum(max(0, b - a - 1) for a, b, _ in path)))
    return prune(pts, lambda p: (p["T"], p["CX"], p["extra"])), ref


def realize(name, start, cache, pt, qasm_out, do_verify=True):
    n, gates, ref = get_start(name, start)
    ops, nb, hmax = assemble(gates, n, cache, pt["path"]); nq = n + hmax
    ops = fold_ops(ops, nq)
    c = counts(ops); emit_qasm(ops, nq, nb, qasm_out)
    fid = verify(qasm_out, ROOT / "bench_qc" / f"{name}.qc") if do_verify and nq <= 22 else None
    return dict(c, extra_qubits=hmax, fidelity=fid, qasm=qasm_out, bits=nb)


def improve(path, qasm_out=None, cx_cap=None, starts=("ours", "raw"), Hmax=8, tl=2.0, restarts=2, log=print, deadline=None, t_slack=2):
    """-> dict(T, CX, extra_qubits, measurements, fidelity, pareto[...]).  Selection: cx_cap given -> min T s.t. CX<=cap;
    else min CX among points with T <= T_min + t_slack."""
    name = Path(path).stem; (_QCACHE / "tg2").mkdir(parents=True, exist_ok=True)
    allpts = []; caches = {}
    for st in starts:
        if st == "ours" and qopt.ours_file(path, OURS / f"{name}.json") is None: continue
        n, gates, ref = get_start(name, st); cp = _QCACHE / "tg2" / f"{name}_{st}.pkl"
        caches[st] = build_cache(gates, n, Hmax, cp, tl, restarts, log, deadline)
        pts, ref = pareto_points(name, st, tuple(range(0, Hmax + 1)), caches[st]); allpts += pts
    allpts = prune(allpts, lambda p: (p["T"], p["CX"], p["extra"]))
    tmin = min(p["T"] for p in allpts)
    ok = [p for p in allpts if (p["CX"] <= cx_cap if cx_cap is not None else p["T"] <= tmin + t_slack)]
    if not ok: ok = allpts
    pt = min(ok, key=(lambda p: (p["T"], p["CX"])) if cx_cap is not None else (lambda p: (p["CX"], p["T"])))
    out = realize(name, pt["start"], caches[pt["start"]], pt, qasm_out or str(HERE / f"{name}_sel.qasm"))
    out["pareto"] = [(p["T"], p["CX"], p["extra"], p["start"]) for p in sorted(allpts, key=lambda p: (p["T"], p["CX"]))]
    return out


if __name__ == "__main__":
    nm = sys.argv[1]; Hm = int(sys.argv[2]) if len(sys.argv) > 2 else 8; tl = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
    starts = tuple(sys.argv[4].split(",")) if len(sys.argv) > 4 else ("ours", "raw")
    r = improve(str(ROOT / "bench_qc" / f"{nm}.qc"), Hmax=Hm, tl=tl, starts=starts, log=lambda s: print(f"[{nm}] {s}", flush=True))
    json.dump(r, open(HERE / f"res_{nm}.json", "w"))
    print(json.dumps({k: v for k, v in r.items() if k != "pareto"}), flush=True)
    print("PARETO", r["pareto"], flush=True)
