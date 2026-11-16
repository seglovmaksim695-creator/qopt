"""One short run on a real IBM Quantum backend: Qiskit level 3 against ``qoptct hw`` on a few small circuits.

NOT RUN BY ANY TEST OR BENCHMARK.  Nothing is submitted unless you pass ``--submit`` AND set QOPT_ALLOW_REAL=1.

    # 1. plan (reads the current calibration of the backend; no job, no QPU time)
    python scripts/hw_real_run.py plan   --backend ibm_torino --out runs/real
    # 2. submit ONE Sampler job (default 8 circuits x 2000 shots, about 10-20 s of QPU time)
    QOPT_ALLOW_REAL=1 python scripts/hw_real_run.py submit --out runs/real --submit
    # 3. fetch the result and print the table
    python scripts/hw_real_run.py result --out runs/real

``plan --backend fake_torino`` works offline (fake backend) and lets you look at the circuits first; ``submit``
refuses a plan made for a fake backend.

For every benchmark two circuits are run with the same input state: A = best of ``--seeds`` Qiskit level-3
transpilations (by expected fidelity), B = the circuit chosen by ``qoptct.hw.compile_hw`` among its own strategies
(Qiskit candidates excluded).  Metric: Hellinger fidelity between the measured and the ideal distribution, and the
probability of the most likely ideal outcome.  Account: ~/.qiskit/quantum-start.json (name "quantum-start").
"""
import argparse
import json
import os
import sys
from pathlib import Path

BENCHES = ["qft_5", "add_cdkm_2", "qaoa_full_5", "heis_ring_6"]


def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    return QiskitRuntimeService(filename=str(Path.home() / ".qiskit/quantum-start.json"), name="quantum-start")


def get_backend(name):
    if name.lower().startswith("fake") or ":" in name:
        from qoptct.hw.device import get_backend as gb
        return gb(name), True
    return service().backend(name), False


def build(prep, comp, init, fin, n, n_phys):
    from qiskit import QuantumCircuit
    c = QuantumCircuit(n_phys, n)
    if prep is not None:
        for ins in prep.data:
            c.append(ins.operation, [init[prep.find_bit(q).index] for q in ins.qubits])
    c.compose(comp, inplace=True)
    c.measure(list(fin), list(range(n)))
    return c


def plan(a):
    from qiskit import qpy, QuantumCircuit, transpile
    from qiskit.quantum_info import Statevector
    from qoptct.hw import benchlib
    from qoptct.hw.compiler import compile_hw
    from qoptct.hw.device import Device
    bk, fake = get_backend(a.backend)
    dev = Device(bk)
    os.makedirs(a.out, exist_ok=True)
    want = a.benches.split(",") if a.benches else BENCHES
    benches = {b.name: b for b in benchlib.suite("tiny") + benchlib.suite("std")}
    items, circs = [], []
    for name in want:
        b = benches[name]
        n = b.circuit.num_qubits
        r = compile_hw(b.circuit, dev, seeds=tuple(range(a.seeds)), seconds=a.seconds)
        ref = r.best["qiskit"]
        new = [c for k, c in r.best.items() if k != "qiskit" and c.verified is True]
        if ref.verified is not True or not new:
            print(f"{name}: skipped (verification: qiskit={ref.verified}, new candidates={len(new)})")
            continue
        new = max(new, key=lambda c: c.metrics["fid"])
        prep = b.prep if b.prep is not None else QuantumCircuit(n)
        ideal = {k: float(v) for k, v in Statevector(prep.compose(b.circuit)).probabilities_dict().items() if v > 1e-12}
        for label, c in (("qiskit_o3", ref), (new.name, new)):
            circ = build(prep, c.circuit, c.init_map, c.final_map, n, dev.n)
            circ = transpile(circ, bk, optimization_level=0, initial_layout=list(range(dev.n)))  # prep gates -> native
            circs.append(circ)
            items.append({"bench": name, "label": label, "n": n, "twoq": c.metrics["twoq"], "depth2q": c.metrics["depth2q"],
                          "exp_fid": c.metrics["fid"], "ideal": ideal})
            print(f"{name:14s} {label:12s} 2q={c.metrics['twoq']:4d} 2q-depth={c.metrics['depth2q']:4d} expected fidelity={c.metrics['fid']:.3f}")
    with open(os.path.join(a.out, "circuits.qpy"), "wb") as f:
        qpy.dump(circs, f)
    json.dump({"backend": dev.name, "fake": fake, "shots": a.shots, "items": items}, open(os.path.join(a.out, "plan.json"), "w"), indent=1)
    print(f"plan: {len(circs)} circuits x {a.shots} shots for {dev.name} written to {a.out} (nothing submitted)")


def submit(a):
    from qiskit import qpy
    P = json.load(open(os.path.join(a.out, "plan.json")))
    if P["fake"]:
        sys.exit("the plan was made for a fake backend; make it for a real one (plan --backend ibm_...)")
    if not a.submit or os.environ.get("QOPT_ALLOW_REAL") != "1":
        sys.exit("refusing to submit: pass --submit and set QOPT_ALLOW_REAL=1 (this spends QPU time of the account)")
    jf = os.path.join(a.out, "job.json")
    if os.path.exists(jf):
        sys.exit(f"{jf} exists: a job was already submitted for this plan")
    with open(os.path.join(a.out, "circuits.qpy"), "rb") as f:
        circs = qpy.load(f)
    from qiskit_ibm_runtime import SamplerV2
    bk = service().backend(P["backend"])
    job = SamplerV2(mode=bk, options={"max_execution_time": a.max_seconds}).run(circs, shots=P["shots"])
    json.dump({"job_id": job.job_id(), "backend": bk.name}, open(jf, "w"), indent=1)
    print("submitted", job.job_id())


def result(a):
    from qiskit.quantum_info import hellinger_fidelity
    P = json.load(open(os.path.join(a.out, "plan.json")))
    J = json.load(open(os.path.join(a.out, "job.json")))
    job = service().job(J["job_id"])
    print("status", job.status())
    res = job.result()
    for it, pub in zip(P["items"], res):
        counts = next(iter(pub.data.values())).get_counts()
        shots = sum(counts.values())
        top = max(it["ideal"], key=it["ideal"].get)
        it["hellinger"] = hellinger_fidelity(counts, it["ideal"])
        it["p_top"] = counts.get(top, 0) / shots
        it["counts"] = counts
        print(f"{it['bench']:14s} {it['label']:12s} 2q={it['twoq']:4d} hellinger={it['hellinger']:.3f} p_top={it['p_top']:.3f} "
              f"(ideal {it['ideal'][top]:.3f}; expected fidelity {it['exp_fid']:.3f})")
    json.dump(P, open(os.path.join(a.out, "results.json"), "w"), indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    q = sp.add_parser("plan")
    q.add_argument("--backend", required=True); q.add_argument("--out", required=True)
    q.add_argument("--benches"); q.add_argument("--seeds", type=int, default=5); q.add_argument("--seconds", type=float, default=15)
    q.add_argument("--shots", type=int, default=2000)
    q = sp.add_parser("submit")
    q.add_argument("--out", required=True); q.add_argument("--submit", action="store_true")
    q.add_argument("--max-seconds", type=int, default=120)
    q = sp.add_parser("result")
    q.add_argument("--out", required=True)
    a = p.parse_args()
    {"plan": plan, "submit": submit, "result": result}[a.cmd](a)


if __name__ == "__main__":
    main()
