"""``qoptct hw``: compile a circuit for a noisy IBM-like device (heavy-hex, CZ/ECR + RZ, SX, X).

    qoptct hw --backend fake_torino in.qasm -o out.qasm [--report out.json] [--measure]
    qoptct hw --list-backends

Input: OpenQASM 2 / 3 or a Feynman .qc file (unitary circuit; final measurements are removed).  The candidates
(Qiskit level 3 for every seed, qoptct logical optimization + Qiskit, the qoptct router with and without the logical
optimization) are compared by the expected fidelity computed from the calibration snapshot of the backend, the
best one is verified against the input (state-vector check on the used qubits) and written.  Fake backends only
replay stored calibration data: this command never contacts IBM Quantum and never submits a job.

Exit code 0: verified result written; 2: result written but not verified (too many used qubits); 3: input error.
"""
import argparse
import json
import sys


def load_circuit(path):
    from qiskit import QuantumCircuit
    if path.endswith(".qc"):
        from .benchlib import from_qc
        return from_qc(path)
    text = open(path).read()
    if "OPENQASM 3" in text[:400]:
        from qiskit import qasm3
        c = qasm3.loads(text)
    else:
        from qiskit import qasm2
        c = qasm2.loads(text, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    c = c.remove_final_measurements(inplace=False)
    out = QuantumCircuit(c.num_qubits, global_phase=c.global_phase)
    for ins in c.data:
        if ins.operation.name == "barrier":
            continue
        out.append(ins.operation, [c.find_bit(q).index for q in ins.qubits])
    return out


def with_measurements(circ, final_map):
    from qiskit import QuantumCircuit
    out = QuantumCircuit(circ.num_qubits, len(final_map), global_phase=circ.global_phase)
    out.compose(circ, inplace=True)
    out.measure(list(final_map), list(range(len(final_map))))
    return out


def main_hw(argv=None):
    p = argparse.ArgumentParser(prog="qoptct hw", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="?", help="circuit: OpenQASM 2/3 (.qasm) or Feynman .qc")
    p.add_argument("--backend", default="fake_torino", help="fake_torino (Heron, CZ), fake_fez, fake_marrakesh, "
                   "fake_sherbrooke (Eagle, ECR), ...; or line:N / ring:N / grid:RxC / heavyhex:D")
    p.add_argument("-o", "--output", help="compiled circuit on the physical qubits: OpenQASM 3 (.qasm) or .qpy")
    p.add_argument("--report", help="JSON report: metrics of every candidate, qubit maps, verification")
    p.add_argument("--seeds", type=int, default=5, help="number of Qiskit level-3 reference runs (default 5)")
    p.add_argument("--seconds", type=float, default=20, help="search time of the qoptct router per logical variant (default 20)")
    p.add_argument("--qopt-seconds", type=float, default=60, help="time for the logical CX/T optimization (default 60)")
    p.add_argument("--objective", choices=["fid", "twoq", "depth"], default="fid",
                   help="fid: expected fidelity from calibration (default); twoq: two-qubit gates; depth: two-qubit depth")
    p.add_argument("--strategies", default="qiskit,qopt+qiskit,hw,qopt+hw", help="comma-separated subset")
    p.add_argument("--no-qopt", action="store_true", help="skip the logical CX/T optimization")
    p.add_argument("--no-verify", action="store_true")
    p.add_argument("--max-verify-qubits", type=int, default=26)
    p.add_argument("--measure", action="store_true", help="append measurements of the logical qubits (bit i = qubit i)")
    p.add_argument("--list-backends", action="store_true")
    a = p.parse_args(argv)
    if a.list_backends:
        try:
            from qiskit_ibm_runtime import fake_provider as fp
            names = sorted(n for n in dir(fp) if n.startswith("Fake") and "Provider" not in n and not n.endswith("V1"))
        except ImportError:
            names = []
        print("\n".join(names + ["line:N", "ring:N", "grid:RxC", "heavyhex:D", "full:N"]))
        return 0
    if not a.input:
        p.error("input circuit required")
    try:
        qc = load_circuit(a.input)
    except Exception as e:
        print(f"error: cannot read {a.input}: {type(e).__name__}: {e}", file=sys.stderr)
        return 3
    from .compiler import compile_hw
    from .device import Device
    try:
        dev = Device(a.backend)
    except (ValueError, ImportError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 3
    try:
        r = compile_hw(qc, dev, seeds=tuple(range(a.seeds)), seconds=a.seconds, objective=a.objective,
                       strategies=a.strategies.split(","), verify=not a.no_verify, qopt=not a.no_qopt,
                       qopt_seconds=a.qopt_seconds, max_verify_qubits=a.max_verify_qubits,
                       log=lambda s: print(s, file=sys.stderr))
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 3
    print(f"backend {dev.name}: {dev.n} qubits, native two-qubit gate {dev.twoq}; circuit: {qc.num_qubits} qubits, {qc.size()} gates")
    print(f"{'strategy':14s} {'2q':>6s} {'2q-depth':>8s} {'depth':>6s} {'exp.fidelity':>12s} {'duration,us':>11s}  verified")
    for nm, c in sorted(r.best.items(), key=lambda kv: -kv[1].metrics["fid"]):
        m = c.metrics
        mark = "*" if c.circuit is r.circuit else " "
        ver = {True: "yes", False: "FAILED", None: "-"}[c.verified]
        print(f"{mark}{nm:13s} {m['twoq']:6d} {m['depth2q']:8d} {m['depth']:6d} {m['fid']:12.4f} {m['duration_us']:11.1f}  {ver}")
    print(f"selected: {r.strategy}; logical qubit i starts on physical init_map[i], is measured on final_map[i]")
    print(f"init_map  = {r.init_map}\nfinal_map = {r.final_map}")
    out = with_measurements(r.circuit, r.final_map) if a.measure else r.circuit
    if a.output:
        if a.output.endswith(".qpy"):
            from qiskit import qpy
            with open(a.output, "wb") as f:
                qpy.dump(out, f)
        else:
            from qiskit import qasm3
            with open(a.output, "w") as f:
                f.write(qasm3.dumps(out))
    if a.report:
        rep = {"backend": dev.name, "input": a.input, "selected": r.strategy, "verified": r.verified,
               "init_map": r.init_map, "final_map": r.final_map, "metrics": r.metrics, "candidates": r.table}
        with open(a.report, "w") as f:
            json.dump(rep, f, indent=1)
    return 0 if (r.verified is True or a.no_verify) else 2


if __name__ == "__main__":
    sys.exit(main_hw())
