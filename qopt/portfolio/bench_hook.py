"""Встраивание принятых кандидатов в bench_real.py (ярус 2: с сетью, на машине Максима; на 24 тестовых —
только ИТОГОВЫЙ замер раз в K кругов, решения по нему не принимаются).
В bench_real.job1, внутри цикла `for _ in range(cycles)`, после global_cancel:
    import sys; sys.path.insert(0, EVOLVE_DIR); import bench_hook
    for f in bench_hook.passes():
        c2 = f(best, n)
        if bs.key(c2) < bs.key(best): best, improved = c2, True
и в main, где добавляется +rt старт, для каждого g in bench_hook.starts():
jobs.append((name+'+'+g.NAME, n, bench_hook.run_start(g, hl, n), False, 4)), где hl = rtof.load_qc(f)[1] (форма CCZ+H).
run_start, как песочница, подаёт fold_h(hl) (цели CCZ выбраны), а кандидатам с HL_RAW=True — исходную форму. Проверка равенства в bench_real (check) остаётся прежней — она и есть ворота яруса 2.
Код кандидатов уже прошёл AST-проверку; здесь он исполняется как обычный модуль (доверие после ворот)."""
import importlib.util, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def _elite(kind):
    import population as P
    for r in P.load():
        if r.get("elite") and r.get("kind") == kind and r.get("name") != "identity":
            yield r


def _mod(path):
    spec = importlib.util.spec_from_file_location(Path(path).stem, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def passes():
    import evo_lib, json
    out = [_mod(r["program"]).run for r in _elite("pass")]
    for r in _elite("rules"):
        rules = json.loads(Path(r["program"]).read_text())
        out.append(lambda c, n, R=rules: evo_lib.apply_rules(c, R))
    return out


def starts():
    return [_mod(r["program"]) for r in _elite("lower")]


def run_start(m, hl, n):
    import evo_lib
    hl = [(a, tuple(b)) for a, b in hl]
    return m.run(hl if getattr(m, "HL_RAW", False) else evo_lib.fold_h(hl), n)
