"""Мост к внешнему оптимизатору VarTODD (ext/VarTodd, C++; собран с разрешения Максима 12.11.2026).

install() подменяет agents/eng1/fasttodd.todd_fast: запускается и наш FastTODD, и VarTODD;
результат VarTODD принимается только если сигнатурный тензор совпадает (тот же фазовый полином
с точностью до клиффордовой части) и, если есть hook (ограничения протокола), он допускает решение.
Берётся более короткий список T-чётностей.
"""
import os, subprocess, tempfile, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
BIN = Path(os.environ.get("QOPT_VARTODD_BIN", "/nonexistent/vartodd_demo"))  # user-built VarTODD binary (optional)
STATS = {"calls": 0, "var_better": 0, "var_fail": 0, "saved_t": 0}


def available():
    return BIN.is_file() and os.access(BIN, os.X_OK)


def run_vartodd(cols, N, seconds=60, seeds=(4, 11, 23)):
    """cols: список int-масок T-чётностей на N проводах. Возвращает лучший (самый короткий) список или None."""
    import fasttodd
    if not cols:
        return []
    A = np.ascontiguousarray(fasttodd.cols_to_A(cols, N).T.astype(np.int8))  # строки = члены, как в VarTODD
    best = None
    with tempfile.TemporaryDirectory(prefix="vt_") as d:
        src = Path(d) / "in.npy"; np.save(src, A)
        t_end = time.time() + seconds
        for s in seeds:
            left = t_end - time.time()
            if left <= 1:
                break
            out = Path(d) / f"out{s}.npy"
            try:
                subprocess.run([str(BIN), "-f", str(src), "-o", str(out), "-s", str(s)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=left, check=True)  # падение = их внутренняя проверка нашла ошибку; такой результат не берём
                B = np.load(out).astype(np.uint8)
            except Exception:
                continue
            new = fasttodd.A_to_cols(B.T) if B.size else []
            if best is None or len(new) < len(best):
                best = new
    return best


def install(seconds=60):
    import fasttodd
    if getattr(fasttodd, "_vartodd_installed", False) or not available():
        return False
    orig = fasttodd.todd_fast

    def todd_fast(cols, N, rng, time_limit=30.0, plateau_max=300, tries=None, log=None, hook=None):
        ours = orig(cols, N, rng, time_limit, plateau_max, tries, log=log, hook=hook)
        if hook is not None:  # слоёный протокол A хранит состояние в hook — там VarTODD не подставляем
            return ours
        STATS["calls"] += 1
        try:
            v = run_vartodd(fasttodd.norm_cols(cols), N, seconds=min(seconds, max(5.0, time_limit * 2)))
        except Exception:
            v = None
        if v is None or len(v) >= len(ours):
            return ours
        if not np.array_equal(fasttodd.tens(v, N), fasttodd.tens(fasttodd.norm_cols(cols), N)):
            STATS["var_fail"] += 1
            return ours
        STATS["var_better"] += 1; STATS["saved_t"] += len(ours) - len(v)
        return v

    fasttodd.todd_fast = todd_fast
    fasttodd._vartodd_installed = True
    return True
