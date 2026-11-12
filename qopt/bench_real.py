"""Этап 12, пункт 2+3 плана автора: реальные алгоритмы (набор Feynman/Amy–Maslov–Mosca, Clifford+T).
Сравнение CX и времени: исходная, Qiskit opt3, PyZX full_optimize, наш блочный (с сетью / без сети),
и наш блочный поверх результата PyZX. Проверка равенства — PyZX verify_equality (+ Operator для n<=10)."""
import glob, json, os, sys, time, warnings
from multiprocessing import Pool
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent
sys.setrecursionlimit(100000)
RESTARTS = int(os.environ.get("RESTARTS", "0")); SKIP_P = float(os.environ.get("SKIP_P", "0.15"))
VARIANTS = os.environ.get("VARIANTS", "ours,nonet,pz")
OUT = os.environ.get("OUT", "bench_real")
KS = (5, 4, 3) if os.environ.get("MODEL5") else (4, 3)


def zx_to_ours(c):
    out = []
    for g in c.to_basic_gates().gates:
        nm = type(g).__name__
        if nm == "HAD": out.append(("h", (g.target,)))
        elif nm == "NOT": out.append(("x", (g.target,)))
        elif nm == "CNOT": out.append(("cx", (g.control, g.target)))
        elif nm == "CZ": out += [("h", (g.target,)), ("cx", (g.control, g.target)), ("h", (g.target,))]
        elif nm in ("T", "S", "Z", "ZPhase"):
            if abs(float(g.phase) * 4 - round(float(g.phase) * 4)) > 1e-9:
                raise ValueError("фаза не кратна pi/4")
            k = int(round(float(g.phase) * 4)) % 8  # фаза в единицах pi/4
            seq = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}[k]
            out += [(s, (g.target,)) for s in seq]
        elif nm == "XPhase":
            k = int(round(float(g.phase) * 4)) % 8
            seq = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}[k]
            out += [("h", (g.target,))] + [(s, (g.target,)) for s in seq] + [("h", (g.target,))]
        else:
            raise ValueError("команда " + nm)
    return out


def ours_to_zx(c, n):
    import pyzx as zx
    z = zx.Circuit(n)
    for nm, qb in c:
        if nm == "h": z.add_gate("HAD", qb[0])
        elif nm == "x": z.add_gate("NOT", qb[0])
        elif nm == "cx": z.add_gate("CNOT", qb[0], qb[1])
        elif nm in ("t", "tdg"): z.add_gate("T", qb[0], adjoint=(nm == "tdg"))
        elif nm in ("s", "sdg"): z.add_gate("S", qb[0], adjoint=(nm == "sdg"))
    return z


def check(c1, c0, n):
    """Сквозная проверка: n<=10 — матрицы целиком; n<=22 — 3 случайных состояния (с общей фазой);
    больше — 'blocks' (каждый заменённый блок уже проверен точно при замене)."""
    from qiskit.quantum_info import Operator, Statevector, random_statevector
    import block_synth as bs
    if n <= 10:
        return bool(Operator(bs.to_qc(c1, n)).equiv(Operator(bs.to_qc(c0, n))))
    if n > 22:  # PyZX: True — доказано; False/неудача — НЕ доказано (None), а не «верно»
        try:
            r = ours_to_zx(c1, n).verify_equality(ours_to_zx(c0, n))
        except Exception:
            r = False
        if r:
            return True
        if n <= 26:  # PyZX не доказал — численно (complex128, 2 случайных состояния, общая фаза), как у агента
            sys.path.insert(0, str(ROOT / "proto2")); from verify import numeric_equal
            return True if numeric_equal(c1, c0, n) else False
        return None
    q1, q0 = bs.to_qc(c1, n), bs.to_qc(c0, n)
    ph = None
    for seed in range(3):
        v = random_statevector(2 ** n, seed=seed)
        o = np.vdot(v.evolve(q0).data, v.evolve(q1).data)
        if abs(abs(o) - 1) > 1e-11 or (ph is not None and abs(o - ph) > 1e-11):  # Г1: 1e-6 пропускал ошибки при n=20–22
            return False
        ph = o
    return True


def cx(c):
    return sum(g[0] == "cx" for g in c)


def _init(model4):
    import block_synth as bs
    bs._init(64, 256, model4, False, 0.5, False)


def job(args):
    """PORTFOLIO=1: несколько нарезок (без/с перестановками, разные зерна), берётся лучшая — нарезка неустойчива."""
    name, n, circ, no_net, cycles = args
    if not os.environ.get("PORTFOLIO"):
        return job1(args, 0)
    t0 = time.time(); best = None
    import block_synth as _bs
    _bs.CACHE.clear(); os.environ["KEEP_CACHE"] = "1"
    variants = ((False, 0), (True, 0), (False, 1)) if RESTARTS else ((False, 0), (True, 0))
    for commute, seed in variants:
        if commute: os.environ["COMMUTE"] = "1"
        else: os.environ.pop("COMMUTE", None)
        r = job1(args, seed)
        if best is None or r[3] and __import__("block_synth").key(r[3]) < __import__("block_synth").key(best[3]):
            best = r
    os.environ.pop("COMMUTE", None); os.environ.pop("KEEP_CACHE", None)
    return best[0], best[1], best[2], best[3], time.time() - t0


def job1(args, seed):
    name, n, circ, no_net, cycles = args
    import block_synth as bs
    bs.W["no_net"] = no_net
    if not os.environ.get("KEEP_CACHE"): bs.CACHE.clear()  # память только внутри схемы — результат не зависит от порядка схем по процессам (шум ±25)
    t0 = time.time(); best = circ
    GC = bool(os.environ.get("GCANCEL"))
    if GC:
        from global_cancel import global_cancel
        if os.environ.get("RTOF"):  # свёртка фаз по всей схеме (агент «нестандартные идеи»)
            sys.path.insert(0, str(ROOT / "proto2")); from phase_fold import fold_cancel
            _gc = global_cancel
            global_cancel = lambda c: _gc(fold_cancel(c, n))
        best = global_cancel(best)
    for _ in range(cycles):
        improved = False
        if GC:
            c2 = global_cancel(best)
            if bs.key(c2) < bs.key(best):
                best, improved = c2, True
        for k in KS:
            for flip in (0, 1):
                c = bs.inv_circ(best) if flip else best
                c2 = bs.one_pass(c, min(k, n))
                c2 = bs.inv_circ(c2) if flip else c2
                if bs.key(c2) < bs.key(best):
                    best, improved = c2, True
        if not improved:
            break
    rng = np.random.default_rng(seed); fails = 0
    while fails < RESTARTS:  # случайные перезапуски нарезки, пока RESTARTS подряд без выигрыша
        k = KS[int(rng.integers(len(KS)))]; flip = int(rng.integers(2))
        c = bs.inv_circ(best) if flip else best
        c2 = bs.one_pass(c, min(k, n), rng, SKIP_P)
        c2 = bs.inv_circ(c2) if flip else c2
        if bs.key(c2) < bs.key(best):
            best, fails = c2, 0
        else:
            fails += 1
    return name, no_net, len(circ), best, time.time() - t0


def main():
    import pyzx as zx
    from qiskit import QuantumCircuit, transpile
    model4 = sys.argv[1] if len(sys.argv) > 1 else ""
    rows, jobs, base = {}, [], {}
    for f in sorted(glob.glob(str(ROOT / "bench_qc" / "*.qc"))):
        name = os.path.basename(f)[:-3]
        if os.environ.get("ONLY") and name not in os.environ["ONLY"].split(","): continue
        z = zx.Circuit.load(f).to_basic_gates(); n = z.qubits
        orig = zx_to_ours(z)
        t = time.time(); zo = zx.optimize.full_optimize(z.copy()); tz = time.time() - t
        pz = zx_to_ours(zo)
        qc = QuantumCircuit(n)
        for nm, qb in orig: getattr(qc, nm)(*qb)
        t = time.time(); qcx = transpile(qc, basis_gates=["u", "cx"], optimization_level=3, seed_transpiler=0).count_ops().get("cx", 0); tq = time.time() - t
        base[name] = {"n": n, "orig": cx(orig), "qiskit": qcx, "t_qiskit": tq, "pyzx": cx(pz), "t_pyzx": tz, "z": z, "orig_c": orig}
        V = VARIANTS.split(",")
        jobs += ([(name, n, orig, False, 4)] if "ours" in V else []) + ([(name, n, orig, True, 4)] if "nonet" in V else []) \
            + ([(name + "+pz", n, pz, False, 4)] if "pz" in V else [])
        ep = ROOT / "results" / ("evo_starts_local.json" if os.environ.get("EVO_STARTS") == "local" else "evo_starts.json")
        if os.environ.get("EVO") and ep.exists():
            _evo = json.loads(ep.read_text())
            if name in _evo:
                ev = [(g, tuple(q)) for g, q in _evo[name]]
                base[name]["evo_start"] = cx(ev)
                jobs.append((name + "+evo", n, ev, False, 4))
        rp = ROOT / "proto2" / f"combo2_{name}.pkl"
        if os.environ.get("RTOF") and rp.exists():
            import pickle
            _n, _raw, rt = pickle.load(open(rp, "rb"))
            base[name]["rt_start"] = cx(rt)
            jobs.append((name + "+rt", n, rt, False, 4))
        print(f"{name}: n={n} CX {cx(orig)} | Qiskit {qcx} ({tq:.1f} с) | PyZX {cx(pz)} ({tz:.1f} с)", flush=True)
    jobs.sort(key=lambda j: -len(j[2]))
    with Pool(8, initializer=_init, initargs=(model4,)) as pool:
        for name, no_net, ln, best, sec in pool.imap_unordered(job, jobs):
            b = name.replace("+pz", "").replace("+rt", "").replace("+evo", "")
            ok = check(best, base[b]["orig_c"], base[b]["n"])
            key = "pz_ours" if name.endswith("+pz") else ("rt_ours" if name.endswith("+rt") else ("evo_ours" if name.endswith("+evo") else ("nonet" if no_net else "ours")))
            base[b][key], base[b]["t_" + key], base[b]["ok_" + key] = cx(best), sec, ok
            base[b]["circ_" + key] = [[g, list(q)] for g, q in best]  # сохранить схему для перепроверки
            print(f"  {name} {'без сети' if no_net else 'с сетью'}: CX {cx(best)} за {sec:.0f} с, проверка {ok}", flush=True)
    out = {k: {kk: vv for kk, vv in v.items() if kk not in ("z", "orig_c")} for k, v in base.items()}
    (ROOT / "results" / f"{OUT}.json").write_text(json.dumps(out, indent=1))
    print("ИТОГ таблица:", flush=True)
    print("схема n | исх | Qiskit | PyZX | наш(сеть) | наш(без сети) | PyZX+наш | время наш/без/Qiskit/PyZX | проверка", flush=True)
    for k, v in out.items():
        print(f"{k} {v['n']} | {v['orig']} | {v['qiskit']} | {v['pyzx']} | {v.get('ours')} | {v.get('nonet')} | {v.get('pz_ours')} | "
              f"{v.get('t_ours',0):.0f}/{v.get('t_nonet',0):.0f}/{v['t_qiskit']:.1f}/{v['t_pyzx']:.1f} | "
              f"{v.get('ok_ours')}/{v.get('ok_nonet')}/{v.get('ok_pz_ours')}", flush=True)
    tot = lambda key: sum(v.get(key, 0) or 0 for v in out.values())
    for v in out.values():  # лучшее из проверенных стартов
        cands = [v[k] for k in ("ours", "rt_ours", "evo_ours") if v.get(k) is not None and v.get("ok_" + k) is True]
        v["best_ours"] = min(cands) if cands else None
    print("СУММА лучшего из стартов (только проверенные):", tot("best_ours"),
          "| rt-старт с сетью:", tot("rt_ours"), "| без проверки:", [k for k, v in out.items() if v["best_ours"] is None], flush=True)
    (ROOT / "results" / f"{OUT}.json").write_text(json.dumps(out, indent=1))
    print(f"СУММА CX: наш {tot('ours')}, Qiskit {tot('qiskit')}, исходные {tot('orig')}; время наш {tot('t_ours'):.0f} с", flush=True)


if __name__ == "__main__":
    main()
