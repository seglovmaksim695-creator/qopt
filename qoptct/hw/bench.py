"""Benchmark driver of ``qoptct hw``.

    python -m qoptct.hw.bench run   --backend fake_torino --out DIR [--suite std] [--shard 0/3] [--seeds 10]
    python -m qoptct.hw.bench sim   --backend fake_torino --out DIR [--shots 4000]     # noisy simulation (Aer)
    python -m qoptct.hw.bench table --out DIR [--md FILE]

``run`` writes one JSON per benchmark (metrics of every strategy, verification status, compile time) and the compiled
circuits (QPY).  Strategies:  qiskit_o3 (seed 0), qiskit_o3_best (best of --seeds seeds by expected fidelity),
qopt+qiskit, hw, qopt+hw  (see compiler.py); ``tket`` rows are added by scripts/hw_tket.py when present.
"""
import argparse
import json
import math
import os
import sys
import time

LABELS = ["qiskit_o3", "qiskit_o3_best", "tket", "qopt+qiskit", "hw", "qopt+hw"]


def _dump(c, d, label, dev, logical, verify_max):
    from qiskit import qpy
    from .check import verify_compiled
    from .compiler import valid_for
    ver = c.verified
    if ver is None:
        ver = verify_compiled(logical, c.circuit, c.init_map, c.final_map, max_qubits=verify_max) \
            if valid_for(c.circuit, dev) else False
    fn = os.path.join(d, label.replace("+", "_") + ".qpy")
    with open(fn, "wb") as f:
        qpy.dump(c.circuit, f)
    m = dict(c.metrics)
    return {"metrics": m, "verified": ver, "seconds": round(c.seconds, 2), "init": list(map(int, c.init_map)),
            "final": list(map(int, c.final_map)), "note": c.note, "file": os.path.basename(fn)}


def run(args):
    from .device import Device
    from . import benchlib
    from .compiler import compile_hw, _key
    from .metrics import logical_stats
    dev = Device(args.backend)
    benches = benchlib.suite(args.suite)
    if args.only:
        want = set(args.only.split(","))
        benches = [b for b in benches if b.name in want]
    si, sn = (int(x) for x in args.shard.split("/"))
    benches = [b for k, b in enumerate(sorted(benches, key=lambda b: b.circuit.size())) if k % sn == si]
    os.makedirs(args.out, exist_ok=True)
    for b in benches:
        d = os.path.join(args.out, b.name)
        jf = os.path.join(d, "result.json")
        if os.path.exists(jf) and not args.force:
            continue
        os.makedirs(d, exist_ok=True)
        t0 = time.time()
        try:
            r = compile_hw(b.circuit, dev, seeds=tuple(range(args.seeds)), seconds=args.seconds, verify=False,
                           qopt_seconds=args.qopt_seconds, alternatives=_alts(b), log=lambda s: print("   ", s))
        except Exception as e:
            print(f"{b.name}: FAILED {type(e).__name__}: {e}")
            continue
        rows = {}
        key = _key("fid")
        q = [c for c in r.candidates if c.name == "qiskit"]
        s0 = next(c for c in q if c.note == "seed=0")
        rows["qiskit_o3"] = _dump(s0, d, "qiskit_o3", dev, b.circuit, args.verify_max)
        qb = min(q, key=key)
        qb.seconds = sum(c.seconds for c in q)
        rows["qiskit_o3_best"] = _dump(qb, d, "qiskit_o3_best", dev, b.circuit, args.verify_max)
        for nm in ("qopt+qiskit", "hw", "qopt+hw"):
            cs = [c for c in r.candidates if c.name == nm]
            if cs:
                cb = min(cs, key=key)
                if nm == "qopt+qiskit":
                    cb.seconds = sum(c.seconds for c in cs)
                rows[nm] = _dump(cb, d, nm, dev, b.circuit, args.verify_max)
        from qiskit import qpy
        with open(os.path.join(d, "logical.qpy"), "wb") as f:
            qpy.dump(b.circuit, f)
        if b.prep is not None:
            with open(os.path.join(d, "prep.qpy"), "wb") as f:
                qpy.dump(b.prep, f)
        res = {"name": b.name, "family": b.family, "backend": dev.name, "logical": logical_stats(b.circuit),
               "rows": rows, "seeds": args.seeds, "seconds": args.seconds, "total_seconds": round(time.time() - t0, 1)}
        with open(jf, "w") as f:
            json.dump(res, f, indent=1)
        line = " | ".join(f"{k}: {v['metrics']['twoq']}/{v['metrics']['fid']:.3f}/{'ok' if v['verified'] else v['verified']}"
                          for k, v in rows.items())
        print(f"{b.name} n={b.circuit.num_qubits} [{time.time() - t0:.0f}s] {line}")
        sys.stdout.flush()


def _alts(b):
    """Published offline qoptct result (results/best/unitary/NAME.qasm) of a .qc benchmark, if present."""
    if not b.name.startswith("qc_"):
        return []
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    f = os.path.join(root, "results", "best", "unitary", b.name[3:] + ".qasm")
    if not os.path.exists(f):
        return []
    try:
        from qiskit import qasm3, QuantumCircuit
        c = qasm3.load(f)
        out = QuantumCircuit(c.num_qubits, global_phase=c.global_phase)
        for ins in c.data:
            out.append(ins.operation, [c.find_bit(q).index for q in ins.qubits])
        return [out]
    except Exception as e:
        print(f"    offline result {f}: {type(e).__name__}: {e}")
        return []


def _load(path):
    from qiskit import qpy
    with open(path, "rb") as f:
        return qpy.load(f)[0]


def sim(args):
    """Noisy simulation with the noise model of the backend (gate errors, readout errors, T1/T2 on scheduled
    idle times).  Metric: Hellinger fidelity between the noisy and the ideal output distribution."""
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector, hellinger_fidelity
    from qiskit.transpiler import PassManager
    from qiskit.transpiler.passes import ALAPScheduleAnalysis, PadDelay
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel
    from .device import Device
    dev = Device(args.backend)
    nm = NoiseModel.from_backend(dev.backend)
    simr = AerSimulator(noise_model=nm, method="statevector", max_parallel_threads=1, seed_simulator=args.seed)
    sched = PassManager([ALAPScheduleAnalysis(target=dev.target), PadDelay(target=dev.target)])
    for name in sorted(os.listdir(args.out)):
        d = os.path.join(args.out, name)
        jf = os.path.join(d, "result.json")
        if not os.path.exists(jf):
            continue
        res = json.load(open(jf))
        if args.only and name not in args.only.split(","):
            continue
        n = res["logical"]["qubits"]
        if n > args.max_qubits:
            continue
        logical = _load(os.path.join(d, "logical.qpy"))
        pf = os.path.join(d, "prep.qpy")
        prep = _load(pf) if os.path.exists(pf) else QuantumCircuit(n)
        ideal = Statevector(prep.compose(logical)).probabilities_dict()
        ideal = {k: v for k, v in ideal.items() if v > 1e-12}
        changed = False
        for label, row in res["rows"].items():
            if "sim" in row and row["sim"].get("shots") == args.shots and not args.force:
                continue
            if row["metrics"]["used_qubits"] > args.max_used or row["verified"] is not True:
                continue
            comp = _load(os.path.join(d, row["file"]))
            c = QuantumCircuit(dev.n, n)
            for ins in prep.data:
                c.append(ins.operation, [row["init"][prep.find_bit(q).index] for q in ins.qubits])
            c.compose(comp, inplace=True)
            c.measure(row["final"], list(range(n)))
            try:
                c = sched.run(c)
            except Exception as e:
                print(f"{name}/{label}: scheduling failed ({e}); simulated without idle noise")
            t0 = time.time()
            counts = simr.run(c, shots=args.shots).result().get_counts()
            hf = hellinger_fidelity(counts, ideal)
            top = max(ideal, key=ideal.get)
            row["sim"] = {"hellinger": hf, "p_top": counts.get(top, 0) / args.shots, "ideal_top": ideal[top],
                          "shots": args.shots, "seconds": round(time.time() - t0, 1)}
            changed = True
            print(f"{name:22s} {label:15s} hellinger={hf:.4f} p_top={row['sim']['p_top']:.4f} (ideal {ideal[top]:.3f})")
            sys.stdout.flush()
        if changed:
            with open(jf, "w") as f:
                json.dump(res, f, indent=1)


def load_results(out):
    rs = []
    for name in sorted(os.listdir(out)):
        jf = os.path.join(out, name, "result.json")
        if os.path.exists(jf):
            rs.append(json.load(open(jf)))
    return rs


def _gmean(xs):
    xs = [x for x in xs if x and x > 0]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float("nan")


def table(args):
    rs = load_results(args.out)
    rs.sort(key=lambda r: (r["family"], r["logical"]["qubits"], r["name"]))
    labels = [l for l in LABELS if any(l in r["rows"] for r in rs)]
    out = []
    w = out.append
    w(f"Backend: {rs[0]['backend'] if rs else '?'}; {len(rs)} circuits.  Cell: two-qubit gates / two-qubit depth / "
      f"expected fidelity (calibration).  `!` = not verified (too many qubits), `X` = verification FAILED.\n")
    w("| circuit | n | " + " | ".join(labels) + " | best new vs O3-best: 2q, fid |")
    w("|---|---|" + "---|" * len(labels) + "---|")
    fam = {}
    for r in rs:
        cells = []
        for l in labels:
            row = r["rows"].get(l)
            if not row:
                cells.append("-")
                continue
            m = row["metrics"]
            flag = "" if row["verified"] is True else (" !" if row["verified"] is None else " X")
            cells.append(f"{m['twoq']} / {m['depth2q']} / {m['fid']:.3f}{flag}")
        ref = r["rows"]["qiskit_o3_best"]["metrics"]
        ref0 = r["rows"]["qiskit_o3"]["metrics"]
        new = [r["rows"][l] for l in ("qopt+qiskit", "hw", "qopt+hw") if l in r["rows"] and r["rows"][l]["verified"] is not False]
        hwonly = [r["rows"][l] for l in ("hw",) if l in r["rows"]]
        if new:
            bn = max(new, key=lambda x: x["metrics"]["fid"])["metrics"]
            r2 = bn["twoq"] / max(ref["twoq"], 1)
            rf = bn["fid"] / max(ref["fid"], 1e-300)
            cells.append(f"{(r2 - 1) * 100:+.0f}%, x{rf:.2f}")
            f = fam.setdefault(r["family"], {"r2": [], "rf": [], "r2_0": [], "rf_0": [], "h2": [], "hf": [], "d2": []})
            f["r2"].append(r2); f["rf"].append(rf)
            f["r2_0"].append(bn["twoq"] / max(ref0["twoq"], 1)); f["rf_0"].append(bn["fid"] / max(ref0["fid"], 1e-300))
            f["d2"].append(bn["depth2q"] / max(ref["depth2q"], 1))
            if hwonly:
                f["h2"].append(hwonly[0]["metrics"]["twoq"] / max(ref["twoq"], 1))
                f["hf"].append(hwonly[0]["metrics"]["fid"] / max(ref["fid"], 1e-300))
        else:
            cells.append("-")
        w(f"| {r['name']} | {r['logical']['qubits']} | " + " | ".join(cells) + " |")
    w("")
    w("Geometric means per family (best verified qoptct strategy relative to Qiskit O3):\n")
    w("| family | circuits | 2q vs O3 best-of-seeds | 2q-depth vs O3 best | fidelity vs O3 best | 2q vs O3 seed 0 | fidelity vs O3 seed 0 | router alone (hw): 2q, fidelity vs O3 best |")
    w("|---|---|---|---|---|---|---|---|")
    for k, f in sorted(fam.items()):
        w(f"| {k} | {len(f['r2'])} | {(_gmean(f['r2']) - 1) * 100:+.1f}% | {(_gmean(f['d2']) - 1) * 100:+.1f}% | x{_gmean(f['rf']):.2f} | "
          f"{(_gmean(f['r2_0']) - 1) * 100:+.1f}% | x{_gmean(f['rf_0']):.2f} | {(_gmean(f['h2']) - 1) * 100:+.1f}%, x{_gmean(f['hf']):.2f} |")
    sims = [r for r in rs if any("sim" in v for v in r["rows"].values())]
    if sims:
        w("")
        w("Noisy simulation (Aer, noise model of the backend, idle T1/T2 included; Hellinger fidelity to the ideal "
          "distribution / probability of the most likely ideal outcome):\n")
        w("| circuit | n | " + " | ".join(labels) + " |")
        w("|---|---|" + "---|" * len(labels))
        for r in sims:
            cells = []
            for l in labels:
                s = r["rows"].get(l, {}).get("sim")
                cells.append(f"{s['hellinger']:.3f} / {s['p_top']:.3f}" if s else "-")
            w(f"| {r['name']} | {r['logical']['qubits']} | " + " | ".join(cells) + " |")
    text = "\n".join(out)
    if args.md:
        with open(args.md, "w") as f:
            f.write(text + "\n")
    print(text)


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m qoptct.hw.bench", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("--backend", default="fake_torino"); r.add_argument("--out", required=True)
    r.add_argument("--suite", default="std"); r.add_argument("--shard", default="0/1")
    r.add_argument("--seeds", type=int, default=10); r.add_argument("--seconds", type=float, default=20)
    r.add_argument("--only"); r.add_argument("--force", action="store_true")
    r.add_argument("--verify-max", type=int, default=26); r.add_argument("--qopt-seconds", type=float, default=60)
    s = sp.add_parser("sim")
    s.add_argument("--backend", default="fake_torino"); s.add_argument("--out", required=True)
    s.add_argument("--shots", type=int, default=4000); s.add_argument("--max-qubits", type=int, default=12)
    s.add_argument("--max-used", type=int, default=20); s.add_argument("--only"); s.add_argument("--seed", type=int, default=1)
    s.add_argument("--force", action="store_true")
    t = sp.add_parser("table")
    t.add_argument("--out", required=True); t.add_argument("--md")
    a = p.parse_args(argv)
    {"run": run, "sim": sim, "table": table}[a.cmd](a)


if __name__ == "__main__":
    main()
