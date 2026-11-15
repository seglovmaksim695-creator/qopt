# SPDX-License-Identifier: Apache-2.0
"""Benchmark driver:  python -m qoptct.contrib.cxpass.bench QASM_DIR QC_DIR [--clean] [--noref] [--budget S] [--prove S] [-o OUT]

For every NAME.qasm in QASM_DIR with a reference QC_DIR/NAME.qc: improve(ref=NAME.qc), then an independent check of the
result (qoptct.verify.prove when available, PROVED/REFUTED/UNKNOWN, plus the dense simulator for <= 22 qubits).
Prints one JSON line per circuit and a markdown table at the end.
"""
import argparse
import json
import pathlib
import time

from . import improve_report, ir


def _independent(text, qc, clean, prove_sec):
    res = {}
    try:
        from qoptct.verify import prove, check_dense
    except ImportError:
        return res
    N = ir.parse(text)[0]
    if N <= 22:
        try:
            res["dense"] = bool(check_dense(text, qc, clean=clean, trials=30)["ok"])
        except Exception as e:  # pragma: no cover
            res["dense"] = repr(e)[:80]
    if prove_sec > 0:
        try:
            r = prove(text, qc, clean=clean, budget_sec=prove_sec)
            res["prove"] = r.get("status"); res["prove_method"] = r.get("method")
        except Exception as e:  # pragma: no cover
            res["prove"] = repr(e)[:80]
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("qasm_dir"); ap.add_argument("qc_dir")
    ap.add_argument("--clean", action="store_true"); ap.add_argument("--budget", type=float, default=300)
    ap.add_argument("--prove", type=float, default=300); ap.add_argument("-o", "--out")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--noref", action="store_true", help="call improve() without the reference (as the core does)")
    a = ap.parse_args(argv)
    rows = []
    for q in sorted(pathlib.Path(a.qasm_dir).glob("*.qasm")):
        nm = q.stem
        if a.only and nm not in a.only:
            continue
        qc = pathlib.Path(a.qc_dir) / f"{nm}.qc"
        if not qc.exists():
            continue
        txt = q.read_text()
        try:
            N = ir.parse(txt)[0]
        except ir.ParseError as e:
            print(json.dumps({"name": nm, "error": str(e)[:120]}), flush=True); continue
        t0 = time.time()
        out, rep = improve_report(txt, a.clean, a.budget, ref=None if a.noref else str(qc))
        row = {"name": nm, "N": N, "in": rep["in"], "out": rep["out"], "verified": rep["verified"],
               "sec": round(time.time() - t0, 1), "passes": [s["pass"] for s in rep["steps"] if s.get("verified")],
               "error": rep.get("error")}
        if out != txt:
            row.update(_independent(out, str(qc), a.clean, a.prove))
            if a.out:
                pathlib.Path(a.out).mkdir(parents=True, exist_ok=True)
                (pathlib.Path(a.out) / f"{nm}.qasm").write_text(out)
        rows.append(row)
        print(json.dumps(row), flush=True)
    print("| circuit | qubits | before T/CX | after T/CX | passes | check |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        chk = "-" if r["out"] == r["in"] else f"own+{'dense ' + str(r.get('dense')) if 'dense' in r else ''} prove {r.get('prove')}"
        print(f"| {r['name']} | {r['N']} | {r['in'][0]}/{r['in'][1]} | {r['out'][0]}/{r['out'][1]} | "
              f"{','.join(r['passes']) or '-'} | {chk} |")


if __name__ == "__main__":
    main()
