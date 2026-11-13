"""finalize with eng2 region resynthesis post-pass"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, time
sys.path.insert(0, (_QPKG + "/contrib/tg3"))
from core import *
import improve as e2


def finalize2(ctx, Winit, Wend, seeds=40, base_seed=1000, e2_rounds=6, e2_budget=2.0, log=None):
    Cs, order = corr_info(ctx, Wend)
    best = synth3(ctx, Winit, Wend, seeds=seeds, base_seed=base_seed)
    body = best[1]
    assert check_body(ctx, body, Winit, Wend)
    N = ctx["N"]; T0 = evo_lib.tcount(body)
    stall = 0
    for r in range(e2_rounds):
        imp = False
        for part in (1, 2):
            try:
                b2, _ = e2.improve_gates(body, N, seed=r * 7 + part, budget=e2_budget, part=part)
            except Exception as ex:
                if log: log(f"e2 err {ex}")
                continue
            b2 = evo_lib.fold_cancel(b2, N)
            if evo_lib.tcount(b2) <= T0 and cost_body(b2) < cost_body(body) and check_body(ctx, b2, Winit, Wend):
                if log: log(f"e2 round {r} part {part}: {cost_body(body)} -> {cost_body(b2)}")
                body = b2; imp = True
        if not imp: break
    prog = final_prog(ctx, body, Wend, Cs, order)
    return body, prog
