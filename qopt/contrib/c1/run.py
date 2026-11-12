"""c1/run.py: improve(path, qasm_out) -> dict(T, CX, extra_qubits, fidelity).

Default (mode="c2"): takes agent c2's measured circuit for this benchmark (agents/c2/out, the point with minimal T, or `src=`)
and re-synthesises its unconditional phase/linear blocks with the c1 synthesiser (resyn.py, resyn2.py, resyn4.py:
free Clifford absorption, |+>-start freedom, X-measured-end freedom -> conditional Z, per-wire barriers).  T never grows;
the variant with the smallest CX that passes the strict verify_qasm (N<=22) is written.
mode="gadget": full H-gadget protocol-B circuit with the c1 body synthesiser (pipe.py); stored results in best_gad/ are used
when present (use_best=True)."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, re, glob, shutil
from pathlib import Path
sys.path.insert(0, (_QPKG + "/contrib/c1"))
sys.path.insert(0, _QPKG)
HERE = Path((_QPKG + "/contrib/c1"))


def count_qasm(path):
    T = CX = 0; nq = 0
    for L in open(path):
        m = re.search(r"qubit\[(\d+)\]", L)
        if m: nq = int(m.group(1))
        CX += len(re.findall(r"\b(cx|cz)\s+q", L)); T += len(re.findall(r"\b(t|tdg)\s+q", L))
    return T, CX, nq


def improve(path, qasm_out, mode="c2", src=None, use_best=True, seconds=300.0):
    import verify_qasm, qopt
    name = Path(path).stem; n = len(qopt.load_qc(path)[0])
    if mode == "clean":
        # clean-ancilla model (wires not in .i start in |0>): functional constructions from best_clean/
        src = HERE / "best_clean" / f"{name}.qasm"
        if not src.exists(): raise FileNotFoundError("no clean construction for " + name)
        shutil.copy(src, qasm_out)
        T, CX, nq = count_qasm(qasm_out)
        fid = verify_qasm.check(qasm_out, path, clean=True) if nq <= 22 else None
        return dict(T=T, CX=CX, extra_qubits=nq - n, fidelity=fid, model="clean")
    if mode == "gadget":
        best = HERE / "best_gad" / f"{name}.qasm"
        if not (use_best and best.exists()):
            import pipe
            best = Path(pipe.main(name, seconds, 1, 1, log=lambda s: None)["qasm"])
        shutil.copy(best, qasm_out)
    else:
        if src is None:
            fs = glob.glob(f"{_QPKG}/contrib/c2/out/{name}_T*.qasm")
            if not fs: raise FileNotFoundError("no c2 circuit for " + name)
            src = min(fs, key=lambda f: count_qasm(f)[:2])
        import resyn, resyn2, resyn4
        cands = [src]
        w = _QCACHE / "c1_work"; w.mkdir(parents=True, exist_ok=True)
        for tag, fn in (("a", lambda o: resyn.run(src, o, 60, log=lambda s: None)),
                        ("b", lambda o: resyn2.run2(src, path, o, 16, 24, log=lambda s: None)),
                        ("c", lambda o: resyn4.run4(src, path, o, 16, 24, log=lambda s: None))):
            o = str(w / f"run_{name}_{tag}.qasm")
            try:
                fn(o); cands.append(o)
            except Exception:
                pass
        cands.sort(key=lambda f: count_qasm(f)[:2])
        for c in cands:
            if count_qasm(c)[2] > 22 or verify_qasm.check(c, path) > 1 - 1e-9:
                shutil.copy(c, qasm_out); break
    T, CX, nq = count_qasm(qasm_out)
    fid = verify_qasm.check(qasm_out, path) if nq <= 22 else None
    return dict(T=T, CX=CX, extra_qubits=nq - n, fidelity=fid)


if __name__ == "__main__":
    print(improve(sys.argv[1], sys.argv[2], *(sys.argv[3:4] or [])))
