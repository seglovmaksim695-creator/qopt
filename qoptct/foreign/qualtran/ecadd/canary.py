"""Канарейки: порча НАСТОЯЩЕЙ переведённой программы шага; каждая обязана быть не PROVED.
Для каждого REFUTED — свидетель: конкретное исполнение испорченной программы (csim) на контрпримере.
Если порча вдруг PROVED — проверка csim на случайных входах: нарушение => ЛОЖНОЕ PROVED (HOLES)."""
import argparse, copy, json, math, random, re, time
from qoptct.verify.bdd import ZERO
from .core import build_program, lower_step, Run, Fail, lsb
from .csim import simulate
from .regions import regions
from .specs import STEPS
from .ispecs import ISTEPS
from . import ec


def _wrap_and(ops, i):
    """and(a,b,t) и её erase -> обе с инверсией a (согласованно): меняется только предикат лестницы."""
    _, a, b, t = ops[i]
    j = next(j for j in range(i + 1, len(ops)) if ops[j][0] == 'erase' and ops[j][1] == t)
    out = ops[:i] + [('x', a), ops[i], ('x', a)] + ops[i + 1:j] + [('x', a), ops[j], ('x', a)] + ops[j + 1:]
    return out


def make_mutations(k, tr, low, rnd, count):
    ops = tr['ops']
    reg = regions(tr)
    live_in = [q for ws in low['inputs'].values() for q in ws]
    muts = []

    def add(name, targeted, fn):
        muts.append((name, targeted, fn))

    def in_region(pattern):
        idx = []
        for key, (lo, hi, _) in reg.items():
            if re.search(pattern, key):
                idx.extend(range(lo, hi))
        return idx

    # --- целевые: стыки с ящиками
    for li, leaf in enumerate(low['leaves']):
        if leaf['name'] != 'opaque':
            continue
        kind = leaf['kind']
        if kind.endswith('_adj'):
            r = 'target' if kind == 'MM_adj' else 'junk'
            def fn(tr, low, li=li, r=r):
                L = list(low['leaves']); L[li] = copy.deepcopy(L[li]); w = L[li]['regs_in'][r]; w[0], w[-1] = w[-1], w[0]
                return tr['ops'], L
            add(f'joint:{kind}#{li}:swap-bits-{r}', True, fn)
            if kind == 'MM_adj':
                def fn2(tr, low, li=li):
                    L = list(low['leaves']); L[li] = copy.deepcopy(L[li]); w = L[li]['regs_in']['x']; w[0], w[1] = w[1], w[0]
                    return tr['ops'], L
                add(f'joint:{kind}#{li}:swap-bits-x', True, fn2)
        else:
            def fn(tr, low, li=li):
                L = list(low['leaves']); L[li] = copy.deepcopy(L[li]); w = L[li]['regs_in']['x']
                L[li]['regs_out']['x'] = w  # THRU: тот же список
                w[-1], w[-2] = w[-2], w[-1]
                return tr['ops'], L
            add(f'joint:{kind}#{li}:swap-low-bits-x', True, fn)
    # --- целевые: поправки стираний в важных областях
    erases = [i for i, op in enumerate(ops) if op[0] == 'erase']
    hot = {1: r'^(20|21|22|23):', 2: r'ControlledViaAnd|ModSub|CModSub', 3: r'CModSub|ModDbl|CModAdd',
           4: r'ModSub', 5: r'ControlledViaAnd|CModNeg|ModAdd|CModSub', 6: r'^(0|1):|CModSub|CModAdd|Equals|MultiAnd'}[k]
    hot_idx = set(in_region(hot))
    pool = [i for i in erases if i in hot_idx] or erases
    for i in rnd.sample(pool, min(7, len(pool))):
        def fn(tr, low, i=i):
            o = list(tr['ops']); t, ((cz, a, b),) = o[i][1], o[i][2]
            c = next(q for q in live_in if q not in (a, b, t))
            o[i] = ('erase', t, (('cz', a, c),)); return o, low['leaves']
        add(f'erase-correction:op{i}', True, fn)
    # --- целевые: редкие входы (проверка на ноль / равенство: инверсия одного управления лестницы)
    zero_regions = in_region(r'MultiControlX|Equals|ModNeg|ControlledViaAnd|Controlled$|MultiAnd')
    ands = [i for i in zero_regions if ops[i][0] == 'and']
    for i in rnd.sample(ands, min(6, len(ands))):
        add(f'rare-input:flip-ladder-control-op{i}', True, lambda tr, low, i=i: (_wrap_and(tr['ops'], i), low['leaves']))
    # --- целевые: перенос через все биты (последние cx/and сумматоров)
    add_regions = [key for key in reg if re.search(r'ModSub|ModAdd|CModSub|CModAdd|ModDbl|ModNeg|CModNeg|Adjoint', key)]
    for key in add_regions[:5]:
        lo, hi, _ = reg[key]
        cand = [i for i in range(lo, hi) if ops[i][0] == 'cx']
        if cand:
            i = cand[len(cand) // 2]
            add(f'carry:drop-cx-mid-of-{key}', True, lambda tr, low, i=i: (tr['ops'][:i] + tr['ops'][i + 1:], low['leaves']))
    # --- случайные
    while len(muts) < count + 8:
        i = rnd.randrange(len(ops))
        kind = rnd.choice(['drop', 'xins', 'retarget'])
        if ops[i][0] == 'opaque':
            continue
        if kind == 'drop':
            add(f'random:drop-op{i}', False, lambda tr, low, i=i: (tr['ops'][:i] + tr['ops'][i + 1:], low['leaves']))
        elif kind == 'xins':
            q = rnd.choice(live_in)
            add(f'random:x-on-{q}-after-op{i}', False,
                lambda tr, low, i=i, q=q: (tr['ops'][:i + 1] + [('x', q)] + tr['ops'][i + 1:], low['leaves']))
        elif ops[i][0] == 'cx':
            q = rnd.choice([w for w in live_in if w not in ops[i][1:]])
            add(f'random:cx-control-op{i}->{q}', False,
                lambda tr, low, i=i, q=q: (tr['ops'][:i] + [('cx', q, tr['ops'][i][2])] + tr['ops'][i + 1:], low['leaves']))
    return muts


def random_witness(prog, k, n, p, patched, rnd, trials=60):
    """Конкретное исполнение испорченной программы на случайных входах из области (ящики — случайные значения <p)."""
    fwd = [it[1]['kind'] for it in prog['items'] if it[0] == 'opaque' and it[1]['kind'] in ('MM', 'KINV')]
    for _ in range(trials):
        ins = {nm: (rnd.randrange(2) if len(ws) == 1 else rnd.randrange(p)) for nm, ws in prog['inputs'].items()}
        ov = [dict(target=rnd.randrange(p), qrom_indices=rnd.randrange(1 << n), reduced=rnd.randrange(2)) if kd == 'MM'
              else dict(x=rnd.randrange(p), junk=rnd.randrange(1 << (4 * n))) for kd in fwd]
        post = _post_int(k, n, p, patched, ins, ov)
        _, viol, _ = simulate(prog, ins, n, p, opaque_values=ov, post=post)
        if viol and (post is not None or viol[0][0] in ('and', 'erase', 'MM_adj:pairing', 'KINV_adj:pairing',
                                                         'MM:domain', 'KINV:domain')):
            return viol
    return None


def neutral_or_hole(prog, k, n, p, patched, rnd, trials=300):
    """Конкретная проверка PROVED-порчи на входах ИЗ области предусловия: случайные (фильтр ispecs.pre) и входы
    шага из цепочки на точках кривых (ящики — смысловые MM/KINV). Нарушение => ложное PROVED."""
    from .lemmas_check import semantic
    from .indep import chain
    MM, KINV = semantic(n, p)
    orc = lambda kind, args: (dict(x=args['x'], y=args['y'], target=MM(args['x'], args['y']), qrom_indices=0, reduced=0)
                              if kind == 'MM' else dict(x=KINV(args['x']), junk=0))
    spec = ISTEPS[k](n, p, patched)
    tried, viol_count, samples = 0, 0, []
    for kind, A, Bc, P, Q in ec.scenarios(p, rnd, 40):
        lr = ec.lam_r(P, A, n, p)
        ins, _ = chain(n, p, patched, P, Q, lr, orc)
        samples.append(ins[k])
    for _ in range(trials):
        samples.append({nm: (rnd.randrange(2) if len(ws) == 1 else rnd.randrange(p)) for nm, ws in prog['inputs'].items()})
    for v in samples:
        v = {nm: v[nm] for nm in prog['inputs']}
        o = []
        cl_all = None
        for _ in range(4):
            cl = spec.calls(v, o + [dict(x=0, target=0)] * 4)
            if len(o) == len(cl):
                break
            kind, args = cl[len(o)]
            o.append(orc(kind, args))
        if not spec.pre(v, o):
            continue
        tried += 1
        _, viol, _ = simulate(prog, v, n, p, opaque_values=o, post=spec.post(v, o))
        if viol:
            viol_count += 1
    return dict(tried=tried, violations=viol_count)


def opaque_values_from(cex, prog):
    vals = []
    c = 0
    for it in prog['items']:
        if it[0] == 'opaque' and it[1]['kind'] in ('MM', 'KINV'):
            regs = ('target', 'qrom_indices', 'reduced') if it[1]['kind'] == 'MM' else ('x', 'junk')
            vals.append({r: int(cex.get(f'o{c}.{r}', 0)) for r in regs})
            c += 1
    return vals


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--step', type=int, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True)
    ap.add_argument('--count', type=int, default=24)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--out', default=None)
    ap.add_argument('--only', default=None, help='только порчи с этими номерами (через запятую): по одной на процесс при n=256')
    a = ap.parse_args(argv)
    n, p = a.n, eval(a.p, {}, {})
    only = {int(t) for t in a.only.split(',')} if a.only else None
    from qoptct.foreign.qualtran.ecadd.run_step import step_bloq, version_info
    from qoptct.foreign.qualtran.ecadd.indep import oracle
    vi = version_info()
    patched = vi['version'] == 'patched'
    rnd = random.Random(a.seed)
    low, tr = lower_step(step_bloq(a.step, n, p))
    muts = make_mutations(a.step, tr, low, rnd, a.count)
    orc = oracle(n, p) if n <= 64 else None
    rows, t0 = [], time.time()
    for idx, (name, targeted, fn) in enumerate(muts):
        if only is not None and idx not in only:
            continue
        row = dict(name=name, targeted=targeted, index=idx)
        try:
            prog = build_program(None, mutate=fn, lowered=(low, tr))
            R = Run(prog, STEPS[a.step](n, p, patched), n, p)
            try:
                row['status'] = R.run()
                # PROVED у порчи: это либо порча, не меняющая смысла на области предусловия, либо ДЫРА доказателя.
                v = neutral_or_hole(prog, a.step, n, p, patched, rnd)
                row['neutral_check'] = v
                if v['violations']:
                    row['status'] = 'HOLE'
                    import os
                    os.makedirs('HOLES', exist_ok=True)
                    json.dump(dict(row, n=n, p=str(p)), open(f'HOLES/hole_step{a.step}_n{n}_{len(rows)}.json', 'w'),
                              ensure_ascii=False, indent=1)
            except Fail as e:
                row['status'] = 'REFUTED'
                row['why'] = str(e)
                cex = e.detail.get('counterexample')
                if (getattr(R, 'viol', None) not in (None, ZERO) and getattr(R, 'A_d', None) is not None
                        and n <= 64):                       # при n=256 подсчёт долей слишком дорог по памяти
                    B = R.B
                    bad = B.AND(R.A_d, R.viol)
                    lv = range(B.nvars)
                    frac = B.satcount(bad, lv) / max(1, B.satcount(R.A_d, lv))
                    row['bad_fraction_log2'] = round(math.log2(frac), 1) if frac > 0 else None
                if cex:
                    ins = {k2: int(cex[k2]) for k2 in prog['inputs']}
                    _, viol, _ = simulate(prog, ins, n, p, opaque_values=opaque_values_from(cex, prog),
                                          post=_post_int(a.step, n, p, patched, ins, opaque_values_from(cex, prog)))
                    row['witness'] = bool(viol)
                    row['witness_kind'] = viol[0][0] if viol else None
                    row['cex'] = {k2: str(v) for k2, v in ins.items()}
                if not row.get('witness'):          # нет контрпримера (напр. отказ пары ящиков) — ищем случайно
                    w = random_witness(prog, a.step, n, p, patched, rnd)
                    row['witness'], row['witness_kind'] = bool(w), (w[0][0] if w else None)
                    row['witness_by'] = 'random search (csim)'
        except Exception as e:
            row['status'] = 'ERROR'
            row['why'] = f'{type(e).__name__}: {str(e)[:200]}'
        rows.append(row)
        print(json.dumps({k2: v for k2, v in row.items() if k2 != 'cex'}, ensure_ascii=False), flush=True)
        R = prog = None                                     # освободить BDD порчи до следующей
        import gc as _gc; _gc.collect()
        if row['status'] == 'HOLE':
            break
        good = [r for r in rows if r['status'] == 'REFUTED' and r.get('witness')]
        if len(good) >= a.count and sum(r['targeted'] for r in good) >= (a.count + 1) // 2:
            break
    proved = [r for r in rows if r['status'] == 'PROVED']
    summ = dict(step=a.step, n=n, p=str(p), version=vi['version'], total=len(rows),
                refuted_with_witness=sum(1 for r in rows if r['status'] == 'REFUTED' and r.get('witness')),
                targeted_refuted=sum(1 for r in rows if r['status'] == 'REFUTED' and r.get('witness') and r['targeted']),
                refuted_without_witness=sum(1 for r in rows if r['status'] == 'REFUTED' and not r.get('witness')),
                proved_neutral=len(proved), holes=sum(1 for r in rows if r['status'] == 'HOLE'),
                errors=sum(1 for r in rows if r['status'] == 'ERROR'),
                seconds=round(time.time() - t0, 1))
    print(json.dumps(summ, ensure_ascii=False))
    if a.out:
        json.dump(dict(summary=summ, rows=rows), open(a.out, 'w'), ensure_ascii=False, indent=1)


def _post_int(k, n, p, patched, ins, ovals):
    spec = ISTEPS[k](n, p, patched)
    try:
        return spec.post(ins, ovals) if spec.pre(ins, ovals) else None
    except Exception:
        return None


if __name__ == '__main__':
    main()
