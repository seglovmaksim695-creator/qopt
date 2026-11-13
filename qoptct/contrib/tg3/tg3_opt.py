"""tg3/improve.py: CX post-optimiser for deterministic (feed-forward) circuits made by eng1/det.py (mode B).  T is never changed.

improve(path, qasm_out) -> dict(T, CX, extra_qubits, measurements, fidelity, ...)
  1. det circuit (eng1 det.improve mode B, or an existing qasm via det_qasm=...);
  2. its body is re-read as (phase polynomial f, wire map W); T-parities are untouched;
  3. gauge search (core.py/search*.py):  F1 free initial ancilla wires (|+> is invariant under CX networks),
     F2 free measured-garbage rows (CX onto an X-measured wire = Pauli kickback), corrections re-derived (det2);
  4. eng2 region resynthesis on the body (accepted only if the polynomial is identical) + duplicate cancellation;
  5. export OpenQASM 3, verify with verify_qasm.py (N<=22) + canary (a removed correction must fail).
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, json, time, re
from pathlib import Path
sys.path.insert(0, (_QPKG + "/contrib/tg3"))
from search2 import *
from search_fin import finalize2

WORK = _QCACHE / "tg3_work"; WORK.mkdir(parents=True, exist_ok=True)
QC = (_QPKG + "/../bench_qc")


def canary(qasm, ctx, tries=6):
    """remove one conditional correction at a time; at least one removal must be detected (fidelity < 1-1e-9)."""
    import verify_qasm
    lines = open(qasm).read().split("\n")
    idx = [i for i, L in enumerate(lines) if L.startswith("if (") and (" cz " in L or " s " in L or " z " in L or " x " in L)]
    detected = 0; tested = 0; worst = 1.0
    for i in idx[:tries]:
        tmp = WORK / "canary.qasm"; tmp.write_text("\n".join(lines[:i] + lines[i + 1:]))
        f = verify_qasm.check(str(tmp), ctx["path"], trials=12)
        tested += 1; worst = min(worst, f)
        if f < 1 - 1e-9: detected += 1
    return dict(canary_tested=tested, canary_detected=detected, canary_worst_fid=worst)


def improve(path, qasm_out=None, seconds=300.0, seed=1, det_qasm=None, verify=True, log=print):
    name = Path(path).stem; t0 = time.time()
    if det_qasm is None:
        sys.path.insert(0, (_QPKG + "/contrib/eng1"))
        import det
        det_qasm = str(WORK / f"{name}_det.qasm")
        r0 = det.improve(path, qasm_out=det_qasm, mode="B")
        log(f"det: T={r0.get('T')} CX={r0.get('CX')}")
    ctx = setup_qasm(det_qasm, path)
    T0, CX0, M0 = cx_count_qasm(det_qasm)
    qasm_out = qasm_out or str(WORK / f"{name}_tg3.qasm")
    c, Wi, We = run_staged(ctx, seconds, seed, log=log)
    body, prog = finalize2(ctx, Wi, We, log=log)
    cnt = post.counts(prog)
    if cnt["CX"] >= CX0 or cnt["T"] > T0:       # never return something worse
        Path(qasm_out).write_text(Path(det_qasm).read_text()); cnt = dict(T=T0, CX=CX0, measurements=M0)
    else:
        post.to_qasm3(prog, ctx["n"], ctx["N"], qasm_out)
    T1, CX1, M1 = cx_count_qasm(qasm_out)
    res = dict(name=name, T=T1, CX=CX1, T_before=T0, CX_before=CX0, extra_qubits=ctx["h"], measurements=M1,
               fidelity=None, seconds=round(time.time() - t0), qasm=qasm_out)
    if verify and ctx["N"] <= 22:
        import verify_qasm
        res["fidelity"] = verify_qasm.check(qasm_out, path)
        res.update(canary(qasm_out, ctx))
    return res


if __name__ == "__main__":
    nm, secs = sys.argv[1], float(sys.argv[2])
    dq = sys.argv[3] if len(sys.argv) > 3 else None
    print(json.dumps(improve(f"{QC}/{nm}.qc", seconds=secs, det_qasm=dq)))
