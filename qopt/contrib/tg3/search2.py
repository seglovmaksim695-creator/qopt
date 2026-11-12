"""tg3/search2.py: staged search (random P2-based initial wires, init-annealing, garbage-annealing, joint)."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, time, random, math
sys.path.insert(0, (_QPKG + "/contrib/tg3"))
from search import *


def random_init(ctx, rng):
    n, N = ctx["n"], ctx["N"]
    cands = [m for m in ctx["P2"] if m >> n]
    rng.shuffle(cands)
    rows = [1 << i for i in range(n)]; extra = []
    for m in cands:
        if len(extra) == N - n: break
        if rank_ok(rows + extra + [m] + [0] * (N - n - len(extra) - 1), len(rows) + len(extra) + 1) if False else _indep(rows + extra, m):
            extra.append(m)
    # fill the rest with unit ancilla rows
    for j in range(n, N):
        if len(extra) == N - n: break
        if _indep(rows + extra, 1 << j): extra.append(1 << j)
    W = rows + extra
    return W if rank_ok(W, N) else identity_init(ctx)


def _indep(rows, m):
    basis = []
    for v in rows:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    for b in basis: m = min(m, m ^ b)
    return m != 0


def anneal(S, Winit, Wend, seconds, mode, T0=0.8):
    S.mi, S.mg = ("i" in mode), ("g" in mode)
    t_start = time.process_time(); t_end = t_start + seconds
    cur = S.cost(Winit, Wend); best = (cur, Winit, Wend)
    while time.process_time() < t_end:
        nb = S.neighbors(Winit, Wend)
        if nb is None: continue
        c = S.cost(*nb)
        if c is None: continue
        frac = (time.process_time() - t_start) / seconds; T = T0 * (1 - frac) + 0.03
        if c <= cur or S.rng.random() < math.exp(-(c - cur) / T):
            Winit, Wend = nb; cur = c
            if c < best[0]: best = (c, Winit, Wend)
    return best


def run_staged(ctx, seconds, seed=0, log=None, body_seeds=None):
    if body_seeds is None: body_seeds = 4 if ctx['N'] <= 20 else 2
    S = Search(ctx, seed=seed, body_seeds=body_seeds); rng = random.Random(seed + 99)
    t0 = time.process_time(); We0 = list(ctx["W"])
    pool = [(S.cost(identity_init(ctx), We0), identity_init(ctx), We0)]
    for _ in range(40):
        Wi = random_init(ctx, rng); c = S.cost(Wi, We0)
        if c is not None: pool.append((c, Wi, We0))
    pool.sort(key=lambda x: x[0]); pool = pool[:4]
    if log: log(f"pool costs {[p[0] for p in pool]}")
    best = None
    for c, Wi, We in pool:
        r = anneal(S, Wi, We, seconds * 0.5 / len(pool), "i")
        if best is None or r[0] < best[0]: best = r
    if log: log(f"after init-anneal {best[0]} ({time.process_time()-t0:.0f}s)")
    r = anneal(S, best[1], best[2], seconds * 0.2, "g")
    if r[0] <= best[0]: best = r
    if log: log(f"after garb-anneal {best[0]} ({time.process_time()-t0:.0f}s)")
    r = anneal(S, best[1], best[2], seconds * 0.3, "ig", T0=0.4)
    if r[0] <= best[0]: best = r
    if log: log(f"after joint {best[0]} ({time.process_time()-t0:.0f}s)")
    return best
