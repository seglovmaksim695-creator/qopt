"""Засев v2: combo2 (proto2/run_combo2.py, перенесён в evo_lib.combo2) — сильнейший известный метод.
Получает hl в исходной форме CCZ+H (HL_RAW=True) и сам выбирает цели CCZ (fold_h + локальный поиск),
пары CCX с относительной фазой (R/R4), варианты непарных CCX, свёртку фаз. На реальных схемах:
barenco_tof_10 110 CX, qcla_mod_7 264 CX. max_evals — потолок оценок (тайм-аут песочницы 30 с на схему)."""
NAME = "combo2"; KIND = "lower"; HL_RAW = True
from evo_lib import combo2


def run(hl, n):
    return combo2(hl, n, rounds=3, max_evals=1500)
