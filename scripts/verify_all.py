#!/usr/bin/env python3
"""Re-verify every circuit listed in results/table.csv and check its T/CX/qubit counts against the table.

  python scripts/verify_all.py [--jobs 2] [--only CIRCUIT] [--model unitary|strict|clean] [--kind best|pareto]

Checkers (independent of the optimizer's data structures):
  * <= 22 qubits: qopt.verify_qasm.check -- state-vector simulation with measurements; strict model by default,
    clean=True for the clean-ancilla model; unitary circuits are checked the same way (no measurements).
    30 random input states (10 for circuits with 20-22 qubits, which take minutes each).
  * > 22 qubits, measurement-based (strict): qopt.contrib.c10.verify_big.check -- sum-over-paths, randomized over
    basis inputs and measurement branches.
  * > 22 qubits, unitary: PyZX verify_equality (600 s limit).
A row passes if the checker accepts (F > 1 - 1e-9) and T, CX and the number of qubits counted from the file
equal the table. Exit code 0 only if every row passes."""
import argparse, csv, re, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REL = Path(__file__).resolve().parents[1]
THRESH = 1 - 1e-9


def counts(text):
    T = len(re.findall(r"\b(t|tdg)\s+q", text))
    CX = len(re.findall(r"\b(cx|cz)\s+q", text))
    N = int(re.search(r"qubit\[(\d+)\]", text).group(1))
    return T, CX, N


def check_row(row):
    import qopt
    from qopt import verify_qasm
    f = REL / row["file"]; qc = REL / "bench_qc" / f"{row['circuit']}.qc"
    t0 = time.time()
    text = f.read_text()
    T, CX, N = counts(text)
    out = dict(row, T_found=T, CX_found=CX, N_found=N)
    out["counts_ok"] = (T == int(row["T"]) and CX == int(row["CX"]) and N == int(row["total_qubits"]))
    try:
        if N <= 22:
            F = float(verify_qasm.check(str(f), str(qc), trials=30 if N < 20 else 10, clean=(row["model"] == "clean")))
            ok, how = F > THRESH, "verify_qasm" + (" clean" if row["model"] == "clean" else "")
        elif row["model"] == "strict":
            from qopt.contrib.c10 import verify_big
            r = verify_big.check(str(f), str(qc), verbose=False)
            F, ok, how = float(r.get("F", 0.0)), bool(r.get("ok")), "verify_big"
        elif row["model"] == "unitary":
            import evo_lib
            names, hl = qopt.load_qc(str(qc))
            gates = [(l.split()[0], tuple(int(x) for x in re.findall(r"q\[(\d+)\]", l)))
                     for l in text.splitlines() if l.strip() and not l.startswith(("OPENQASM", "include", "qubit", "bit"))]
            import bench_real as br  # PyZX only: a numerical check at 24-26 qubits needs several GB of memory
            n = len(names); ref = evo_lib.lower_std(hl)
            ok = bool(qopt.limited(600, lambda: br.ours_to_zx(gates, n).verify_equality(br.ours_to_zx(ref, n))))
            F, how = (1.0 if ok else 0.0), "pyzx verify_equality"
        else:
            ok, F, how = False, 0.0, "no checker for this model above 22 qubits"
    except Exception as e:  # noqa
        ok, F, how = False, 0.0, "error: " + repr(e)[:160]
    out.update(verified=ok, F_found=F, checker_used=how, sec=round(time.time() - t0, 1))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--only"); ap.add_argument("--model"); ap.add_argument("--kind")
    ap.add_argument("--table", default=str(REL / "results" / "table.csv"))
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.table)))
    # each file once (the best file is also listed under its Pareto name)
    seen, uniq = set(), []
    for r in rows:
        if (r["model"], r["file"]) not in seen:
            seen.add((r["model"], r["file"])); uniq.append(r)
    rows = [r for r in uniq if (not a.only or r["circuit"] == a.only) and (not a.model or r["model"] == a.model)
            and (not a.kind or r["kind"] == a.kind)]
    rows.sort(key=lambda r: -int(r["total_qubits"]))  # big circuits first so the pool stays busy
    bad = 0
    with ProcessPoolExecutor(a.jobs) as ex:
        for r in ex.map(check_row, rows):
            good = r["verified"] and r["counts_ok"]
            bad += not good
            print(f"{'PASS' if good else 'FAIL'} {r['model']:8s} {r['kind']:6s} {r['circuit']:15s} T={r['T_found']:>4} "
                  f"CX={r['CX_found']:>5} N={r['N_found']:>3} F={r['F_found']:.15f} {r['checker_used']} {r['sec']}s"
                  + ("" if r["counts_ok"] else f"  COUNTS DIFFER FROM TABLE (T={r['T']}, CX={r['CX']}, N={r['total_qubits']})"),
                  flush=True)
    print(f"{len(rows) - bad}/{len(rows)} passed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
