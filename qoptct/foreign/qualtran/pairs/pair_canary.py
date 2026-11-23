"""pair_canary -- канарейки для зеркальной проверки пары Калиски: портим ОДИН элемент ОБРАТНОЙ программы (кусок-итерация, макро-узел
или верхний уровень) и требуем, чтобы проверка НЕ дала PASS. Отдельно «равносильные» перестановки (должны проходить).
python -m pair_canary n p [seed] [out.json]"""
import sys, json, random, copy, time
from . import pair_compose as PC
from . import mirror as M
from .mirror import K_X, K_CX, K_AND, K_ERASE


def clone_piece(P, ops):
    Q = copy.copy(P)
    Q.ops = ops
    return Q


def piece_canaries(ops, rng, nwires):
    """-> список (имя, ожидание, новая программа куска). ожидание: 'FAIL' (должно быть не PASS) или 'EQUIV' (равносильно)."""
    out = []
    N = len(ops)
    wires = sorted({w for o in ops for w in M.wires(o)})

    def pick(kind):
        idx = [i for i, o in enumerate(ops) if o[0] == kind]
        if not idx:  # в абстрактном куске Монтгомери нет вентилей этого вида: берём cx (название строки может говорить «x»)
            idx = [i for i, o in enumerate(ops) if o[0] == K_CX]
        return rng.choice(idx)
    # 1 drop gate (cx)
    i = pick(K_CX); o = ops[:]; del o[i]; out.append(('piece: drop cx #%d' % i, 'FAIL', o))
    # 2 drop x
    i = pick(K_X); o = ops[:]; del o[i]; out.append(('piece: drop x #%d' % i, 'FAIL', o))
    # 3 drop AND
    i = pick(K_AND); o = ops[:]; del o[i]; out.append(('piece: drop AND #%d' % i, 'FAIL', o))
    # 4 переставить два несоседних некоммутирующих вентиля
    for tries in range(2):
        for _ in range(2000):
            i = rng.randrange(N - 3)
            j = rng.randrange(i + 2, min(N, i + 60))
            if (ops[i] != ops[j] and not M.commute(ops[i], ops[j]) and ops[i][0] in (K_X, K_CX) and ops[j][0] in (K_X, K_CX)
                    and nwires not in M.wires(ops[i]) and nwires not in M.wires(ops[j])):
                o = ops[:]; o[i], o[j] = o[j], o[i]
                out.append(('piece: swap non-adjacent #%d,#%d' % (i, j), 'FAIL', o))
                break
    # 5 change control wire of cx
    for tries in range(2):
        i = pick(K_CX)
        k, a, b, t, l = ops[i]
        na = rng.choice([w for w in wires if w not in (a, b)])
        o = ops[:]; o[i] = (k, na, b, t, l)
        out.append(('piece: change control wire of cx #%d (%d->%d)' % (i, a, na), 'FAIL', o))
    # 6 испортить erase: сменить control a
    for tries in range(2):
        i = pick(K_ERASE)
        k, a, b, t, l = ops[i]
        na = rng.choice([w for w in wires if w not in (a, b, t)])
        o = ops[:]; o[i] = (k, na, b, t, l)
        out.append(('piece: erase #%d control a %d->%d (correction)' % (i, a, na), 'FAIL', o))
    # 7 AND -> ERASE той же формы (подмена вида)
    i = pick(K_AND); k, a, b, t, l = ops[i]; o = ops[:]; o[i] = (K_ERASE, a, b, t, -1)
    out.append(('piece: AND #%d replace by erase' % i, 'FAIL', o))
    # 8 target AND на другой провод
    i = pick(K_AND); k, a, b, t, l = ops[i]
    nt = rng.choice([w for w in wires if w not in (a, b, t)])
    o = ops[:]; o[i] = (k, a, b, nt, l)
    out.append(('piece: target AND #%d %d->%d' % (i, t, nt), 'FAIL', o))
    # 9 reverse cx
    i = pick(K_CX); k, a, b, t, l = ops[i]; o = ops[:]; o[i] = (k, b, a, t, l)
    out.append(('piece: reverse cx #%d' % i, 'FAIL', o))
    # 10 extra gate (дубль x)
    i = pick(K_X); o = ops[:]; o.insert(i, ops[i]); out.append(('piece: extra x #%d' % i, 'FAIL', o))
    # равносильные: соседние коммутирующие
    cnt = 0
    for _ in range(5000):
        i = rng.randrange(N - 1)
        if ops[i] != ops[i + 1] and M.commute(ops[i], ops[i + 1]) and ops[i][0] in (K_CX, K_X) and ops[i + 1][0] in (K_CX, K_X):
            o = ops[:]; o[i], o[i + 1] = o[i + 1], o[i]
            out.append(('piece: EQUIVALENT: adjacent commuting #%d,#%d' % (i, i + 1), 'EQUIV', o))
            cnt += 1
            if cnt == 2:
                break
    # AND: поменять местами контроли a,b
    i = pick(K_AND); k, a, b, t, l = ops[i]; o = ops[:]; o[i] = (k, b, a, t, l)
    out.append(('piece: EQUIVALENT: controls a,b of AND #%d' % i, 'EQUIV', o))
    return out


def macro_canaries(C, rng):
    n = C.n
    out = []
    bidx = [i for i, m in enumerate(C.macros) if m['kind'] == 'b']
    Pb = C.Pb
    inw = {r: w for r, w in Pb.ins.items()}

    def mutated(fn):
        ms = list(C.macros)
        fn(ms)
        return ms

    # 11 сдвиг провода на стыке итераций: у обратного макро-узла j подменить глобальный провод входа of register на провод соседа
    for tries in range(2):
        j = rng.choice(bidx)
        reg = rng.choice([r for r in inw if len(inw[r]) >= 2])
        w1, w2 = rng.sample(inw[reg], 2)

        def f(ms, j=j, w1=w1, w2=w2):
            m = dict(ms[j]); g = dict(m['g']); g[w1], g[w2] = g[w2], g[w1]; m['g'] = g; ms[j] = m
        out.append(('macro: joint, backward node #%d: swap wires %d,%d of register %s' % (j, w1, w2, reg), 'FAIL', mutated(f)))
    # 12 перестановка бит junk: у двух соседних обратных узлов поменять провод of register m
    for tries in range(2):
        k = rng.randrange(len(bidx) - 1)
        j1, j2 = bidx[k], bidx[k + 1]
        reg = 'm' if 'm' in inw else rng.choice(list(inw))
        w = inw[reg][0]

        def f(ms, j1=j1, j2=j2, w=w):
            a, b = dict(ms[j1]), dict(ms[j2]); ga, gb = dict(a['g']), dict(b['g'])
            ga[w], gb[w] = gb[w], ga[w]; a['g'], b['g'] = ga, gb; ms[j1], ms[j2] = a, b
        out.append(('macro: junk, nodes #%d,#%d: swap bit %s' % (j1, j2, reg), 'FAIL', mutated(f)))
    # 13 порядок макро-узлов: поменять два соседних обратных узла местами
    k = rng.randrange(len(bidx) - 1)
    j1, j2 = bidx[k], bidx[k + 1]

    def f(ms, j1=j1, j2=j2):
        ms[j1], ms[j2] = ms[j2], ms[j1]
    out.append(('macro: swap order of backward nodes #%d,#%d' % (j1, j2), 'FAIL', mutated(f)))
    return out


def top_canaries(C, rng):
    out = []
    T = C.top_ops
    cb = {C.top_cut[i] for i, m in enumerate(C.macros) if m['kind'] == 'b'}.pop()
    Fe = C.F_end
    regions = [('B_pre', Fe, cb), ('B_post', cb, len(T))]
    for nm, lo, hi in regions:
        if hi - lo < 4:
            continue
        i = rng.randrange(lo, hi)
        t = T[:]; del t[i]
        out.append(('top %s: drop gate #%d' % (nm, i), 'FAIL', t))
        i = rng.randrange(lo, hi)
        k, a, b, tt, l = T[i]
        wires = sorted({w for o in T[lo:hi] for w in M.wires(o)})
        t = T[:]
        if k == K_X:
            t[i] = (k, rng.choice([w for w in wires if w != a]), 0, 0, -1)
        elif k == K_CX:
            t[i] = (k, rng.choice([w for w in wires if w not in (a, b)]), b, 0, -1)
        else:
            t[i] = (k, rng.choice([w for w in wires if w not in (a, b, tt)]), b, tt, l)
        out.append(('top %s: change wire of gate #%d' % (nm, i), 'FAIL', t))
    for _ in range(2000):
        lo, hi = regions[0][1], regions[0][2]
        if hi - lo < 4:
            lo, hi = regions[1][1], regions[1][2]
        i = rng.randrange(lo, hi - 1)
        if T[i] != T[i + 1] and not M.commute(T[i], T[i + 1]):
            t = T[:]; t[i], t[i + 1] = t[i + 1], t[i]
            out.append(('top: swap dependent adjacent #%d,#%d' % (i, i + 1), 'FAIL', t))
            break
    return out


def main(n, p, seed=1, outf=None, only=None):
    """Run the canaries of the Kaliski pair.  ``only``: set of job indexes (default all)."""
    rng = random.Random(seed)
    t0 = time.time()
    C = PC.PairComposed(n, p)
    base = PC.check_composed(C)
    res = dict(n=n, p=p, seed=seed, base_mirror_ok=base['mirror_ok'], base_stats=base.get('piece_stats'), rows=[])
    print('base', base['mirror_ok'], round(time.time() - t0, 1), flush=True)
    piece_res = PC.piece_pi(C.Pf, C.Pb)[0]
    jobs = []
    for nm, ex, newops in piece_canaries(C.Pb.ops, rng, C.Pb.nsrc):
        jobs.append((nm, ex, dict(Pb=clone_piece(C.Pb, newops))))
    for nm, ex, ms in macro_canaries(C, rng):
        jobs.append((nm, ex, dict(macros=ms, piece_res=piece_res)))
    for nm, ex, t in top_canaries(C, rng):
        jobs.append((nm, ex, dict(T=t, piece_res=piece_res)))
    import gc
    res['jobs_total'] = len(jobs)
    for idx, (nm, ex, kw) in enumerate(jobs):
        if only is not None and idx not in only:
            continue
        gc.collect()
        t1 = time.time()
        try:
            r = PC.check_composed(C, **kw)
            ok = bool(r['mirror_ok'])
            why = r.get('error') or r.get('piece_msg') if not ok else 'ok'
            if not ok and not r.get('error') and r.get('piece_mirror_ok') and not r.get('top_ok', True):
                why = [m for m in r.get('top_match', []) if not m[1]]
        except Exception as e:  # noqa: BLE001
            ok, why = False, 'exception %s: %s' % (type(e).__name__, str(e)[:200])
        verdict = ('PASS' if ok else 'NOT PASS')
        good = (not ok) if ex == 'FAIL' else ok
        row = dict(index=idx, name=nm, expect=ex, verdict=verdict, good=good, why=str(why)[:300], sec=round(time.time() - t1, 1))
        res['rows'].append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)
    nf = [r for r in res['rows'] if r['expect'] == 'FAIL']
    res['n_fail_canaries'] = len(nf)
    res['n_caught'] = sum(1 for r in nf if r['verdict'] != 'PASS')
    res['false_pass'] = [r['name'] for r in nf if r['verdict'] == 'PASS']
    res['equiv'] = [(r['name'], r['verdict']) for r in res['rows'] if r['expect'] == 'EQUIV']
    res['sec'] = round(time.time() - t0, 1)
    print(json.dumps({k: v for k, v in res.items() if k != 'rows'}, ensure_ascii=False), flush=True)
    if outf:
        json.dump(res, open(outf, 'w'), ensure_ascii=False, indent=1)
    return res
