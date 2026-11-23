"""Независимая симуляция: классический симулятор САМОГО Qualtran на раскрытом шаге (ящики MM/KINV — их собственным
on_classical_vals) против целочисленной спецификации ispecs. Плюс прогон цепочки шести шагов на точках кривых."""
import argparse, json, random, time, sys
from . import ec
from .ispecs import ISTEPS
from .opaque import opaque_kind

NAMES = {1: ('a', 'b', 'x', 'y'), 2: ('f1', 'ctrl', 'a', 'b', 'x', 'y', 'lam_r'), 3: ('ctrl', 'a', 'b', 'x', 'y', 'lam'),
         4: ('x', 'y', 'lam'), 5: ('ctrl', 'a', 'b', 'x', 'y', 'lam_r', 'lam'),
         6: ('f1', 'f2', 'f3', 'f4', 'ctrl', 'a', 'b', 'x', 'y')}


def oracle(n, p):
    from qualtran.bloqs.mod_arithmetic import DirtyOutOfPlaceMontgomeryModMul as MMc, KaliskiModInverse as KI
    mm, ki = MMc(bitsize=n, window_size=1, mod=p), KI(bitsize=n, mod=p)

    def call(kind, args):
        if kind == 'MM':
            return dict(zip(('x', 'y', 'target', 'qrom_indices', 'reduced'), mm.call_classically(x=args['x'], y=args['y'])))
        return dict(zip(('x', 'junk'), ki.call_classically(x=args['x'])))
    return call


def expected(spec, v, orc):
    o = []
    for _ in range(4):
        cl = spec.calls(v, o + [dict(x=0, target=0)] * 4)
        if len(o) == len(cl):
            break
        kind, args = cl[len(o)]
        o.append(orc(kind, args))
    return o, spec.post(v, o)


def chain(n, p, patched, P, Q, lr, orc, keep=None):
    """Шесть int-спецификаций подряд. Возвращает (входы шагов, итог)."""
    a, b = ec.enc(P, n, p)
    x, y = ec.enc(Q, n, p)
    S = {k: ISTEPS[k](n, p, patched) for k in ISTEPS}
    ins = {}
    v1 = dict(a=a, b=b, x=x, y=y); ins[1] = v1
    o1, r1 = expected(S[1], v1, orc)
    v2 = dict(f1=r1['f1'], ctrl=r1['ctrl'], a=a, b=b, x=x, y=y, lam_r=lr); ins[2] = v2
    o2, r2 = expected(S[2], v2, orc)
    v3 = dict(ctrl=r2['ctrl'], a=a, b=b, x=r2['x'], y=r2['y'], lam=r2['lam']); ins[3] = v3
    o3, r3 = expected(S[3], v3, orc)
    v4 = dict(x=r3['x'], y=r3['y'], lam=r3['lam']); ins[4] = v4
    o4, r4 = expected(S[4], v4, orc)
    v5 = dict(ctrl=r3['ctrl'], a=a, b=b, x=r4['x'], y=r4['y'], lam_r=lr, lam=r4['lam']); ins[5] = v5
    o5, r5 = expected(S[5], v5, orc)
    pre5 = S[5].pre(v5, o5)
    v6 = dict(f1=r2['f1'], f2=r1['f2'], f3=r1['f3'], f4=r1['f4'], ctrl=r5['ctrl'], a=a, b=b, x=r5['x'], y=r5['y'])
    ins[6] = v6
    pre6 = S[6].pre(v6, [])
    r6 = S[6].post(v6, [])
    return ins, dict(pre5=pre5, pre6=pre6, out=(r6['x'], r6['y']))


def flat_step(k, n, p):
    from qoptct.foreign.qualtran.ecadd.run_step import step_bloq
    st = step_bloq(k, n, p)
    return st.decompose_bloq().flatten(lambda bi: opaque_kind(bi.bloq) is None)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True)
    ap.add_argument('--steps', default='1,2,3,4,5,6')
    ap.add_argument('--scen', type=int, default=40)
    ap.add_argument('--rand', type=int, default=40)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default=None)
    ap.add_argument('--semantic-oracle', action='store_true',
                    help='значения ящиков — по формулам (MM=u*v*R^-1, KINV=u^-1*R^2), а не классикой Qualtran '
                         '(нужно для исходной 0.7.0 при n=256: её классика умножителя падает, находка 04)')
    a = ap.parse_args(argv)
    n, p = a.n, eval(a.p, {}, {})
    from qoptct.foreign.qualtran.ecadd.run_step import version_info
    vi = version_info()
    patched = vi['version'] == 'patched'
    rnd = random.Random(a.seed)
    orc = oracle(n, p)
    if a.semantic_oracle:
        from qoptct.foreign.qualtran.ecadd.lemmas_check import semantic
        MMs, KIs = semantic(n, p)
        orc = lambda kind, args: (dict(x=args['x'], y=args['y'], target=MMs(args['x'], args['y']), qrom_indices=0, reduced=0)
                                  if kind == 'MM' else dict(x=KIs(args['x']), junk=0))
    t0 = time.time()
    res = dict(n=n, p=str(p), **vi, method='Qualtran classical simulation of flattened step vs ispecs', steps={},
               semantic_oracle=a.semantic_oracle)
    # 1) цепочка на сценариях: входы шагов + проверка итогового сложения
    samples = {k: [] for k in range(1, 7)}
    chain_rows = []
    for kind, A, B, P, Q in ec.scenarios(p, rnd, a.scen):
        lr = ec.lam_r(P, A, n, p)
        if P is None and rnd.random() < 0.5 and kind != 'G3like':
            lr = rnd.randrange(p)
        ins, fin = chain(n, p, patched, P, Q, lr, orc)
        want = ec.enc(ec.add(P, Q, A, p), n, p)
        chain_rows.append(dict(kind=kind, ok=(fin['out'] == want and fin['pre5'] and fin['pre6']), **{k: fin[k] for k in ('pre5', 'pre6')}))
        for k in samples:
            samples[k].append(('chain:' + kind, ins[k]))
    res['chain'] = dict(total=len(chain_rows), ok=sum(r['ok'] for r in chain_rows),
                        bad=sorted({r['kind'] for r in chain_rows if not r['ok']}))
    # 2) случайные входы (с редкими значениями) — для шагов 1-5
    special = [0, 1, p - 1, p - 2, (p - 1) // 2, (p + 1) // 2]
    for k in (1, 2, 3, 4, 5):
        spec = ISTEPS[k](n, p, patched)
        got = 0
        while got < a.rand:
            v = {nm: rnd.randrange(2) if nm in ('f1', 'f2', 'f3', 'f4', 'ctrl') else
                 (rnd.choice(special) if rnd.random() < 0.3 else rnd.randrange(p)) for nm in NAMES[k]}
            if k == 4:
                v['y'] = rnd.randrange(1 << n)
            if k == 5:
                o, _ = expected(ISTEPS[5](n, p, patched), dict(v, lam=0), orc)
                K, Z = o[0]['x'], o[1]['target']
                v['lam'] = v['ctrl'] * (Z ^ (v['lam_r'] if K == 0 else 0))
            samples[k].append(('rand', v)); got += 1
    csim_prog = {}
    for k in [int(s) for s in a.steps.split(',')]:
        spec = ISTEPS[k](n, p, patched)
        cb = flat_step(k, n, p)
        rights = [r.name for r in cb.signature.rights()]
        stat = dict(samples=0, skipped_pre=0, mismatch=0, errors=0, examples=[])
        for tag, v in samples[k]:
            o, want = expected(spec, v, orc)
            if not spec.pre(v, o):
                stat['skipped_pre'] += 1
                continue
            stat['samples'] += 1
            try:
                out = cb.call_classically(**{nm: v[nm] for nm in NAMES[k]})
                got = {nm: int(val) for nm, val in zip(rights, out)}
                if got != want:
                    stat['mismatch'] += 1
                    if len(stat['examples']) < 3:
                        stat['examples'].append(dict(tag=tag, v={kk: str(vv) for kk, vv in v.items()},
                                                     diff={kk: (str(got.get(kk)), str(want.get(kk))) for kk in want if got.get(kk) != want[kk]}))
            except Exception as e:
                stat['errors'] += 1
                if len(stat['examples']) < 3:
                    stat['examples'].append(dict(tag=tag, error=f'{type(e).__name__}: {str(e)[:200]}'))
                # запасной независимый путь: конкретный исполнитель переведённой программы (ew/csim) с ящиками Qualtran
                try:
                    if k not in csim_prog:
                        from qoptct.foreign.qualtran.ecadd.core import build_program
                        from qoptct.foreign.qualtran.ecadd.run_step import step_bloq
                        csim_prog[k] = build_program(step_bloq(k, n, p))
                    from qoptct.foreign.qualtran.ecadd.csim import simulate
                    outs, viol, _ = simulate(csim_prog[k], {nm: v[nm] for nm in NAMES[k]}, n, p, orc=orc, post=want)
                    stat['csim_fallback_ok' if not viol else 'csim_fallback_bad'] = \
                        stat.get('csim_fallback_ok' if not viol else 'csim_fallback_bad', 0) + 1
                except Exception as e2:
                    stat['csim_fallback_error'] = f'{type(e2).__name__}: {str(e2)[:120]}'
        stat['seconds'] = round(time.time() - t0, 1)
        res['steps'][k] = stat
        print(k, stat, flush=True)
    res['wall_seconds'] = round(time.time() - t0, 1)
    print(json.dumps(dict(chain=res['chain']), ensure_ascii=False))
    if a.out:
        json.dump(res, open(a.out, 'w'), ensure_ascii=False, indent=1)


if __name__ == '__main__':
    main()
