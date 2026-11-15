"""bandc: scalable "logical AND + measurement-based uncomputation" compiler (lab agent b1), up to ~1500 qubits.

API::

    from qoptct.contrib import bandc
    cands = bandc.candidates("circuit.qc", clean=True, budget_sec=120)
    # -> [{"qasm": <OpenQASM 3 text>, "note": str}, ...]   (Pareto points in (T, CX, extra qubits), best T first)

The compiler (``compiler.py``) works on classical reversible circuits (x / cx / ccx / z / cz / ccz with H-frames).
Every returned candidate has been verified inside the worker:

1. ``lcheck`` -- EXACT check of the logical operations for all inputs (GF(2) polynomials), plus two op-level canaries
   (drop a CZ of a measurement correction, drop a CX) that must be rejected;
2. a check of the emitted QASM itself: the statevector checker ``qoptct.verify_qasm`` for <= 16 qubits, otherwise
   ``vpt`` (sampled path sum at the classical output point, c10 engine), plus QASM canaries
   (t -> tdg must be rejected; at least one of four dropped-CX mutants must be rejected).

The search and the checks run in a child process with a hard time limit (``budget_sec``) and a memory limit
(env ``QOPT_BANDC_MEM_MB``, default 3000); on timeout the process is killed and the candidates verified so far are
returned. Unsupported, non-classical or too large circuits (> 1600 wires / 40000 gates) give ``[]``.
Models: strict (``clean=False``): all .qc wires carry arbitrary data; clean: wires not listed in ``.i`` start in |0>.
"""
import json
import os
import subprocess
import sys
import tempfile
import time

__all__ = ["candidates"]


def _pareto(rows):
    keys = [(r["T"], r["CX"], r["extra"]) for r in rows]
    out = []
    for r, k in zip(rows, keys):
        if any(o != k and o[0] <= k[0] and o[1] <= k[1] and o[2] <= k[2] for o in keys): continue
        if any((x["T"], x["CX"], x["extra"]) == k for x in out): continue
        out.append(r)
    return sorted(out, key=lambda r: (r["T"], r["CX"], r["extra"]))


def candidates(qc_path: str, clean: bool, budget_sec: float = 60.0) -> list:
    """Verified candidates for the .qc circuit (see module docstring). Never raises; [] if nothing found."""
    try:
        if not os.path.isfile(qc_path) or budget_sec < 3: return []
        budget = float(budget_sec)
        mem = int(os.environ.get("QOPT_BANDC_MEM_MB", "3000"))
        fd, out = tempfile.mkstemp(suffix=".jsonl", prefix="bandc_"); os.close(fd)
        try:
            cmd = [sys.executable, "-m", "qoptct.contrib.bandc.worker", os.path.abspath(qc_path), out,
                   "--budget", "%.1f" % max(1.0, budget - 2.0), "--mem-mb", str(mem)]
            if clean: cmd.append("--clean")
            t0 = time.time()
            p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL,
                                 start_new_session=True)
            try:
                p.wait(timeout=max(1.0, budget - (time.time() - t0)))
            except subprocess.TimeoutExpired:
                p.kill(); p.wait()
            rows = []
            with open(out) as f:
                for line in f:
                    try:
                        r = json.loads(line)
                        if isinstance(r, dict) and isinstance(r.get("qasm"), str) and r["qasm"].startswith("OPENQASM"):
                            rows.append(r)
                    except ValueError:
                        continue          # a line cut by the kill
        finally:
            try: os.unlink(out)
            except OSError: pass
        return [{"qasm": r["qasm"], "note": r.get("note", "bandc")} for r in _pareto(rows)]
    except Exception:  # noqa: BLE001  (contract: never raise)
        return []
