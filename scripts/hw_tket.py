"""Add a ``tket`` row to the results of ``python -m qoptct.hw.bench run`` (run in an environment with pytket and
pytket-qiskit; qoptct must be importable for the metrics and the equivalence check).

    python scripts/hw_tket.py --backend fake_torino --out DIR [--levels 2,3] [--only a,b]

The pass is pytket-qiskit's own default compilation pass for IBM backends (``IBMQBackend.pass_from_info``: noise-aware
placement, routing, synthesis; optimisation level 2 = FullPeepholeOptimise, level 3 = GreedyPauliSimp), built
offline from the calibration snapshot of the fake backend.  The row holds the better of the requested levels
(by expected fidelity); nothing is sent to IBM.
"""
import argparse
import json
import math
import os
import sys
import time

from qiskit import QuantumCircuit, qpy


def backend_info(dev):
    from pytket.extensions.qiskit import IBMQBackend
    b = dev.backend
    return IBMQBackend._get_backend_info(b.configuration(), b.properties())


def tk_compile(logical, info, level, n_phys, timeout):
    from pytket.extensions.qiskit import IBMQBackend, qiskit_to_tk
    from pytket.circuit import OpType
    from pytket.predicates import CompilationUnit
    tk = qiskit_to_tk(logical)
    cu = CompilationUnit(tk)
    IBMQBackend.pass_from_info(info, level, timeout).apply(cu)
    c = cu.circuit
    out = QuantumCircuit(n_phys)
    for cmd in c.get_commands():
        t = cmd.op.type
        if t in (OpType.noop, OpType.Barrier, OpType.Phase):
            continue
        qs = [q.index[0] for q in cmd.qubits]
        if t == OpType.Rz:
            out.rz(float(cmd.op.params[0]) * math.pi, qs[0])
        elif t == OpType.SX:
            out.sx(qs[0])
        elif t == OpType.X:
            out.x(qs[0])
        elif t == OpType.CZ:
            out.cz(*qs)
        elif t == OpType.ECR:
            out.ecr(*qs)
        elif t == OpType.CX:
            out.cx(*qs)
        elif t in (OpType.noop, OpType.Barrier, OpType.Phase):
            continue
        else:
            raise ValueError(f"unexpected tket op {t}")
    lq = {q: q.index[0] for q in tk.qubits}  # single register: tket qubit -> logical index
    init = [None] * logical.num_qubits
    fin = [None] * logical.num_qubits
    for q, node in cu.initial_map.items():
        if q in lq:
            init[lq[q]] = node.index[0]
    for q, node in cu.final_map.items():
        if q in lq:
            fin[lq[q]] = node.index[0]
    # a logical qubit without gates may be left unplaced by tket: park it on a free physical qubit
    used = {q.index[0] for q in c.qubits if q.reg_name == "node"} | {x for x in init + fin if x is not None}
    if any(q.reg_name != "node" for cmd in c.get_commands() for q in cmd.qubits):
        raise ValueError("tket left gates on unplaced qubits")
    free = (p for p in range(n_phys) if p not in used)
    for i in range(logical.num_qubits):
        if init[i] is None or fin[i] is None:
            init[i] = fin[i] = next(free)
    return out, init, fin


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--backend", default="fake_torino"); p.add_argument("--out", required=True)
    p.add_argument("--levels", default="2,3"); p.add_argument("--only"); p.add_argument("--timeout", type=int, default=120)
    p.add_argument("--verify-max", type=int, default=26); p.add_argument("--force", action="store_true")
    a = p.parse_args()
    from qoptct.hw.device import Device
    from qoptct.hw.metrics import analyze
    from qoptct.hw.check import verify_compiled
    from qoptct.hw.compiler import valid_for
    dev = Device(a.backend)
    info = backend_info(dev)
    for name in sorted(os.listdir(a.out)):
        d = os.path.join(a.out, name)
        jf = os.path.join(d, "result.json")
        if not os.path.exists(jf) or (a.only and name not in a.only.split(",")):
            continue
        res = json.load(open(jf))
        if "tket" in res["rows"] and not a.force:
            continue
        with open(os.path.join(d, "logical.qpy"), "rb") as f:
            logical = qpy.load(f)[0]
        best = None
        t0 = time.time()
        for lv in (int(x) for x in a.levels.split(",")):
            try:
                c, init, fin = tk_compile(logical, info, lv, dev.n, a.timeout)
            except Exception as e:
                print(f"{name}: tket level {lv} failed: {type(e).__name__}: {str(e)[:200]}")
                continue
            m = analyze(c, dev, fin)
            if best is None or m["fid"] > best[3]["fid"]:
                best = (c, init, fin, m, lv)
        if best is None:
            continue
        c, init, fin, m, lv = best
        ver = verify_compiled(logical, c, init, fin, max_qubits=a.verify_max) if valid_for(c, dev) else False
        with open(os.path.join(d, "tket.qpy"), "wb") as f:
            qpy.dump(c, f)
        res = json.load(open(jf))
        res["rows"]["tket"] = {"metrics": m, "verified": ver, "seconds": round(time.time() - t0, 2), "init": init,
                               "final": fin, "note": f"pytket level {lv} (best of {a.levels})", "file": "tket.qpy"}
        with open(jf, "w") as f:
            json.dump(res, f, indent=1)
        print(f"{name}: tket L{lv} 2q={m['twoq']} d2q={m['depth2q']} fid={m['fid']:.4f} verified={ver} ({time.time() - t0:.0f}s)")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
