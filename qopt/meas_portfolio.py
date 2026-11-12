"""Этап 5 (2): режим с измерениями + портфель qopt на кусках уровня Тоффоли между вставками AND/MZ.
Каждый кусок оптимизируется qopt.optimize_hl (доказывается сам), затем свёртка фаз по сегменту между измерениями,
итог — симуляцией с измерениями (gidney.verify, до 24 кубитов). python meas_portfolio.py схема.qc ..."""
import json, os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "portfolio"))
os.environ.setdefault("QOPT_PROG_SEC", "60")
import qopt, gidney, evo_lib


def build_low(names, items, n):
    """Куски уровня Тоффоли между AND/MZ -> qopt.optimize_hl (каждый доказывается), сегменты -> свёртка фаз."""
    segs, cur, piece = [], [], []

    def flush():
        if piece:
            if any(g in ("ccx", "ccz") for g, _ in piece):
                _, ref, best, how = qopt.optimize_hl(names, list(piece))
                cur.extend(best if how else ref)
            else:
                cur.extend(evo_lib.lower_std(piece))
            piece.clear()
    for g, q in items:
        if g == "AND":
            flush(); cur.extend(gidney.and_gate(*q))
        elif g == "MZ":
            flush(); cur.append(("h", (q[2],))); segs.append((cur, q)); cur = []
        else:
            piece.append((g, q))
    flush(); segs.append((cur, None))
    return [(qopt.limited(qopt.FOLD_SEC, evo_lib.fold_cancel, s, n), mz) for s, mz in segs]


def run(path):
    names, hl0 = qopt.load_qc(path); n = len(names); ins = gidney.inputs_of(path)
    anc = {i for i, x in enumerate(names) if x not in ins}
    items, r1, r2 = gidney.transform(evo_lib.fold_h(hl0), n, anc)
    low = build_low(names, items, n)
    ok = gidney.verify(hl0, low, n, anc) if n <= 24 else None
    return dict(T=sum(evo_lib.tcount(s) for s, _ in low), CX=sum(evo_lib.cx(s) for s, _ in low), verified=ok)


if __name__ == "__main__":
    g = json.load(open(ROOT / "results" / "gidney.json")); out = {}
    for p in sys.argv[1:]:
        nm = Path(p).stem
        try:
            r = run(p)
        except Exception as e:
            r = {"err": repr(e)[:100]}
        r.update(T_plain=g.get(nm, {}).get("T_meas"), CX_plain=g.get(nm, {}).get("CX_meas"))
        out[nm] = r; print(nm, json.dumps(r, ensure_ascii=False), flush=True)
        (ROOT / "results" / "meas_portfolio.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
