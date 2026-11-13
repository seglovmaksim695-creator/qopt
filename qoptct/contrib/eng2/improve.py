"""eng2: пост-проход ресинтеза CX-фазовых областей поверх лучшей схемы. Вход/выход: плоский список [(имя,(кубиты))]."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, json, time, random
sys.path.insert(0, _QPKG); sys.path.insert(0, (_QPKG + '/portfolio')); sys.path.insert(0, (_QPKG + '/contrib/eng2'))
import evo_lib
from regions import split_regions
import synth as S

def tcnt(g): return sum(1 for a, _ in g if a in ('t', 'tdg'))
def cxc(g): return sum(1 for a, _ in g if a == 'cx')

import os
BM = int(os.environ.get('BM', '1'))
INV = {"h": "h", "x": "x", "s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t", "cx": "cx"}
def _inv(c): return [(INV[a], q) for a, q in reversed(c)]

def resynth_region(gates, qlist, rng, budget, fast=False):
    base_cx, base_t = cxc(gates), tcnt(gates)
    m, terms, wires = S.poly_of(gates, qlist)
    gi = _inv(gates); mi, termsi, wiresi = S.poly_of(gi, qlist)
    best = None; bk = (base_cx, base_t)
    t0 = time.time(); k = 0
    configs = [(4, 1.0, 0.0), (8, 0.7, 0.3), (12, 1.0, 0.5), (6, 0.5, 0.6), (16, 1.0, 1.0), (24, 0.8, 1.5)]
    while time.time() - t0 < budget:
        beam, alpha, noise = configs[(k // 2) % len(configs)]; inv = k % 2; k += 1
        beam = beam * BM
        if m > 9: beam = min(beam, 2)
        if inv: res = S.synth(mi, termsi, wiresi, rng, beam=beam, alpha=alpha, noise=noise, limit=base_cx)
        else: res = S.synth(m, terms, wires, rng, beam=beam, alpha=alpha, noise=noise, limit=base_cx)
        if res is None:
            if k >= 400: break
            continue
        if inv: res = _inv(res)
        t2, w2 = S.apply_circ(res, m)
        if t2 != terms or [tuple(w) for w in w2] != [tuple(w) for w in wires]:
            raise RuntimeError('region mismatch')
        key = (cxc(res), tcnt(res))
        if key < bk: bk = key; best = res
        if m <= 3 and k >= 8: break
    if best is None: return None
    return [(nm, tuple(qlist[i] for i in qs)) for nm, qs in best]

def improve_gates(gates, n, seed=0, budget=2.0, maxm=12, log=None, part=1):
    rng = random.Random(seed)
    from regions import split_regions2
    items = split_regions(gates, n) if part == 1 else split_regions2(gates, n); out = []; gained = 0
    for it in items:
        if it[0] == 'G': out.append(it[1]); continue
        _, rid, ql, gs = it
        cand = None
        if cxc(gs) >= 3 and len(ql) <= maxm:
            cand = resynth_region(gs, ql, rng, budget)
        if cand is not None:
            gained += cxc(gs) - cxc(cand)
            out += cand
        else: out += gs
    return out, gained

def improve(names, best_gates, ref_gates, seed=0, budget=3.0, rounds=3, maxm=13):
    """names: список имён кубитов (или число n). best_gates: лучшая схема (плоский список). ref_gates: evo_lib.lower_std(hl) исходной.
    Возвращает схему, ДОКАЗАННО равную ref (qoptct.prove) и не хуже по (CX, T); иначе возвращает best_gates без изменений."""
    import qoptct
    n = len(names) if not isinstance(names, int) else names
    best = [(a, tuple(q)) for a, q in best_gates]; cur = best
    for r in range(rounds):
        imp = False
        for part in (1, 2):
            new, _ = improve_gates(cur, n, seed, budget, maxm, part=part); seed += 1
            new = evo_lib.fold_cancel(new, n)
            if evo_lib.key(new) < evo_lib.key(cur): cur = new; imp = True
        if not imp: break
    if evo_lib.key(cur) >= evo_lib.key(best) or evo_lib.tcount(cur) > evo_lib.tcount(best): return best
    try:
        ok = qoptct.prove(cur, ref_gates, n)
    except Exception:
        ok = None
    return cur if ok is True else best
