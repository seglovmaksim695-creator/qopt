"""eng1: drivers that start from (a) raw Toffoli-level circuit, (b) our current best circuit (cached by cache_ours.py)."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, json, time, os, random
sys.path.insert(0, (_QPKG + "/contrib/eng1"))
import gad
from fasttodd import todd_fast
from gad import *


def improve_low(gates, ref, n, time_limit=20.0, restarts=3, seeds=6, log=None, tag="", fast=True, equiv_time=0.0,
                todd_tl=3.0, tries=None):
    """gates: any list over {h,x,cx,t,tdg,s,sdg,z,cz,ccz}; ref: lowered reference (n qubits) for the numerical check."""
    import equiv
    t0 = time.time()
    N, P, W, K = build(gates, n)
    assert rank_ok(W, N)
    odd, even = split_odd(P)
    res = dict(n=n, h=N - n, N=N, t_poly=len(odd), start=tag)
    best = None
    for r in range(restarts):
        rng = random.Random(100 + r)
        new = todd_fast(odd, N, rng, time_limit, 300, tries) if fast else todd(odd, N, rng, time_limit)
        P2 = clifford_rest(odd, new, even); steps = []
        if equiv_time > 0:
            P2, steps = equiv.search(P2, W, N, n, rng, equiv_time, todd_tl, log=log)
        body = synth(P2, W, K, N, n, seeds=seeds, base_seed=r * 100)
        if body is None: continue
        k = (evo_lib.tcount(body), evo_lib.cx(body))
        if log: log(f"  [{tag}] restart {r}: cols {len(odd)}->{k[0]} T={k[0]} CX={k[1]} steps={len(steps)} ({time.time()-t0:.0f}s)")
        if best is None or k < best[0]: best = (k, body, steps)
    k, body, steps = best
    res.update(T=k[0], CX=k[1], extra_qubits=N - n, measurements=N - n, equiv_steps=len(steps))
    c1, W1, K1 = extract(body, N)
    ok = (c1 == equiv.certify_chain(P, W, n, N, steps)) and W1 == list(W) and K1[:n] == list(K[:n])
    res.update(alg_proof=bool(ok), sim=check_sim(body, ref, N, n), seconds=round(time.time() - t0, 1))
    return res, body


def improve_path(path, ours_dir=(_QPKG + "/contrib/eng1/ours"), **kw):
    names, hl = qoptct.load_qc(path); n = len(names); nm = Path(path).stem
    ref = evo_lib.lower_std(hl)
    cands = []
    r, b = improve_low(prep_gates(hl), ref, n, tag="raw", **kw); cands.append((r, b))
    f = qoptct.ours_file(path, Path(ours_dir) / f"{nm}.json")
    if f is not None:
        d = json.load(open(f)); best = [(g, tuple(q)) for g, q in d["best"]]
        r2, b2 = improve_low(best, ref, n, tag="ours", **kw)
        r2["ours_T"], r2["ours_CX"] = d["t"], d["cx"]; cands.append((r2, b2))
    cands.sort(key=lambda x: (x[0]["T"], x[0]["CX"]))
    res, body = cands[0]
    res["name"] = nm; res["all"] = [(c[0]["start"], c[0]["T"], c[0]["CX"]) for c in cands]
    return res, body


if __name__ == "__main__":
    for p in sys.argv[1:]:
        r, b = improve_path(p, log=lambda s: print(s, flush=True)); print(json.dumps(r), flush=True)


def improve(path, qasm_out=None, **kw):
    """Entry point for qopt: improve(path) -> dict(T, CX, extra_qubits, measurements, alg_proof, sim, ...).
    Circuit semantics: n data wires + h ancillas (|+>, i.e. H on |0>), the returned body, then H and Z-measurement of the
    h ancillas; the unitary is obtained on outcome 0^h (postselected / Clifford fix-up, see REPORT.md)."""
    r, body = improve_path(path, **kw)
    names, _ = qoptct.load_qc(path); n = len(names); N = n + r["extra_qubits"]
    r["body_len"] = len(body)
    if qasm_out: to_qasm3(qasm_out, body, n, N)
    return r


def to_qasm3(path, body, n, N):
    L = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{N}] q;", f"bit[{N-n}] c;"]
    L += [f"h q[{a}];" for a in range(n, N)]
    Q = {"h": "h", "x": "x", "z": "z", "s": "s", "sdg": "sdg", "t": "t", "tdg": "tdg", "cx": "cx", "cz": "cz"}
    for g, q in body: L.append(f"{Q[g]} " + ", ".join(f"q[{x}]" for x in q) + ";")
    for a in range(n, N): L += [f"h q[{a}];", f"c[{a-n}] = measure q[{a}];"]
    L.append("// unitary = branch with all c == 0")
    Path(path).write_text("\n".join(L) + "\n")
