"""E2/E3 glue in Z3 (SMT-LIB2, integers): the six step contracts (zspecs -- the same record as specs.py/ispecs.py)
+ explicit lemmas.  Every obligation is checked separately: UNSAT of its negation.  The wiring between the steps is
read from the real ECAdd (see wiring.py) and must equal the reference wiring."""
import argparse, json, os, subprocess, time
from .zspecs import SmtOps, IntOps, ZSTEPS

WIRING = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'wiring_ecadd.json')))
# wiring canaries: the glue must refute them (shows that the wiring is really read)
WMUT = {'swap-xy-step4': lambda w: w['steps']['4'].update(x=w['steps']['4']['y'], y=w['steps']['4']['x']),
        'step6-x-from-4': lambda w: w['steps']['6'].update(x=[4, 'x']),
        'out-y-from-4': lambda w: w['out'].update(y=[4, 'y']),
        'step5-x-from-3': lambda w: w['steps']['5'].update(x=[3, 'x'])}

from . import _z3
LEMMAS = ['MM0', 'K0', 'C1', 'C2', 'EC_DBL', 'EC_TAN', 'EC_RNOTO', 'EC_RNOTP']
LEMMA_TEXT = {
    'MM0': 'MM(u,v)=0 <=> u=0 or v=0 (u,v<p, p prime; MM(u,v)=u*v*R^-1 mod p, the meaning of contract C)',
    'K0': 'KINV(u)=0 <=> u=0 (KINV(u)=u^-1*R^2 mod p for u!=0, KINV(0)=0, the meaning of contract D)',
    'C1': 'MM(u, MM(KINV(u), v)) = v for u!=0',
    'C2': 'MM(KINV(u), MM(u, w)) = w for u!=0',
    'EC_DBL': 'P,Q on the curve and finite, a=x, b!=-y => y=b and b!=0',
    'EC_TAN': 'ordinary branch and binding: x_R=a => Lam=lam_r (the line through P,Q is tangent to the curve at P)',
    'EC_RNOTO': 'ordinary branch and binding: (x_R,y_R) != (0,0) (R is on the curve, (0,0) is not)',
    'EC_RNOTP': 'ordinary branch and binding: (x_R,y_R) != (a,b) (otherwise Q=O)',
}


def build(n, p, versions, lemmas, extra=None, wiring=None):
    Wg = wiring or WIRING
    o = SmtOps(n, p)
    L = o.lines
    a, b, x, y, lr = (o.word(k) for k in ('a', 'b', 'x', 'y', 'lam_r'))
    V = {k: o.word(k) for k in ('K2', 'Z2', 'Z1', 'Z3', 'Z3p', 'K5', 'Z5')}
    obl, calls = [], []

    def step(k, W, O):
        pre, cl, post = ZSTEPS[k](o, W, O, versions[k])
        for j, t in enumerate(pre):
            obl.append((f'pre:step{k}#{j}', t))
        for j, (kind, args) in enumerate(cl):
            for r, w in args.items():
                obl.append((f'arg<p:step{k}.{kind}{j}.{r}', o.lt(w, p)))
            calls.append((f'step{k}.{kind}{j}', kind, args, O[j]))
        return post

    # The wiring between the steps comes from the real composite block of Qualtran (ecadd/wiring.py -> data/wiring_ecadd.json,
    # checked for both versions at n=16 and 256), not hand-written: the input of a step = the output of the step it is connected to.
    inp = dict(a=a, b=b, x=x, y=y, lam_r=lr)
    Obox = {1: [], 2: [V['K2'], V['Z2']], 3: [V['Z1']], 4: [V['Z3'], V['Z3p']], 5: [V['K5'], V['Z5']], 6: []}
    posts = {}
    src = lambda s: inp[s[1]] if s[0] == 'in' else posts[s[0]][s[1]]
    for k in Wg['order']:
        posts[k] = step(k, {r: src(s) for r, s in Wg['steps'][str(k)].items()}, Obox[k])
    OUT = {r: src(s) for r, s in Wg['out'].items()}
    e1, e2 = o.eq(a, x), o.eq(b, o.modneg(y))
    e3, e4 = o.and_(o.eqz(a), o.eqz(b)), o.and_(o.eqz(x), o.eqz(y))
    gen = o.and_(o.not_(e3), o.not_(e4), o.not_(o.and_(e1, e2)))
    Lam = o.ite(e1, lr, V['Z2'])
    XR = o.modsub(o.modsub(V['Z3'], a), x)
    YR = o.modsub(V['Z3p'], b)
    cm = {c[0]: c for c in calls}
    k2, m2, m4a, m4b = cm['step2.KINV0'], cm['step2.MM1'], cm['step4.MM0'], cm['step4.MM1']
    bind = o.and_(o.eq(k2[2]['x'], o.modsub(x, a)), o.eq(m2[2]['x'], V['K2']), o.eq(m2[2]['y'], o.modsub(y, b)),
                  o.eq(m4a[2]['x'], Lam), o.eq(m4a[2]['y'], Lam), o.eq(m4b[2]['x'], o.modsub(a, XR)),
                  o.eq(m4b[2]['y'], Lam))
    gb = o.and_(gen, bind)
    obl.append(('bind', f'(=> {gen} {bind})'))
    if 'MM0' in lemmas:
        for tag, kind, args, out in calls:
            if kind == 'MM':
                L.append(f"(assert (= (or (= {args['x']} 0) (= {args['y']} 0)) (= {out} 0)))  ; MM0 {tag}")
    if 'K0' in lemmas:
        for tag, kind, args, out in calls:
            if kind == 'KINV':
                L.append(f"(assert (= (= {args['x']} 0) (= {out} 0)))  ; K0 {tag}")
    if 'C1' in lemmas:   # MM(u, MM(KINV(u), v)) = v  (step 3 after step 2)
        tA, _, aA, oA = cm['step3.MM0']
        u = k2[2]['x']
        L.append(f"(assert (=> (and (= {aA['x']} {u}) (not (= {u} 0)) (= {aA['y']} {m2[3]}) (= {m2[2]['x']} {k2[3]})) "
                 f"(= {oA} {m2[2]['y']})))  ; C1")
    if 'C2' in lemmas:   # MM(KINV(u), MM(u, w)) = w  (step 5 after step 4)
        k5, m5 = cm['step5.KINV0'], cm['step5.MM1']
        u = k5[2]['x']
        L.append(f"(assert (=> (and (= {m5[2]['x']} {k5[3]}) (= {m4b[2]['x']} {u}) (not (= {u} 0)) "
                 f"(= {m5[2]['y']} {m4b[3]})) (= {m5[3]} {m4b[2]['y']})))  ; C2")
    if 'EC_DBL' in lemmas:
        L.append(f"(assert (=> (and {e1} (not {e2}) (not {e3}) (not {e4})) (and (= {y} {b}) (not (= {b} 0)))))  ; EC_DBL")
    if 'EC_TAN' in lemmas:
        L.append(f"(assert (=> (and {gb} (= {XR} {a})) (= {Lam} {lr})))  ; EC_TAN")
    if 'EC_RNOTO' in lemmas:
        L.append(f"(assert (=> {gb} (not (and (= {XR} 0) (= {YR} 0)))))  ; EC_RNOTO")
    if 'EC_RNOTP' in lemmas:
        L.append(f"(assert (=> {gb} (not (and (= {XR} {a}) (= {YR} {b})))))  ; EC_RNOTP")
    fx = o.ite(e3, x, o.ite(e4, a, o.ite(o.and_(e1, e2), 0, XR)))
    fy = o.ite(e3, y, o.ite(e4, b, o.ite(o.and_(e1, e2), 0, YR)))
    for k2_, t in (('x', fx), ('y', fy), ('a', a), ('b', b), ('lam_r', lr)):
        obl.append((f'out:{k2_}', f'(= {OUT[k2_]} {t})'))
    F = dict(e1=e1, e2=e2, e3=e3, e4=e4, a=a, b=b, x=x, y=y, lam_r=lr)
    if extra is not None:
        for t in extra(F):
            L.append(f'(assert {t})  ; extra-pre')
    return o, obl


def groups(F):
    e1, e2, e3, e4 = F['e1'], F['e2'], F['e3'], F['e4']
    b, x, y, lr = F['b'], F['x'], F['y'], F['lam_r']
    g1 = f'(and (not {e3}) (not {e4}) (not {e1}) {e2})'
    g2 = f'(and (not {e3}) {e1} (= {b} {y}) (not (= {b} 0)) (= {lr} 0))'
    g3 = f'(and {e3} (not (= {x} 0)) (not (= {y} 0)) (= {lr} 0))'
    return g1, g2, g3


def solve(lines, checks, timeout_s=600):
    return _z3.solve(lines, checks, timeout_s)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True)
    ap.add_argument('--mode', default='E2', choices=['E2', 'E3', 'E3fail'])
    ap.add_argument('--orig', default='')
    ap.add_argument('--drop', default='')
    ap.add_argument('--group', type=int, default=0)
    ap.add_argument('--out', default=None)
    ap.add_argument('--smt', default=None)
    ap.add_argument('--wiring', default=None, help='wiring JSON read from the real ECAdd (ecadd.wiring --out); '
                    'default: the reference wiring shipped in data/')
    ap.add_argument('--wmut', default=None, help='wiring canary: ' + ', '.join(WMUT))
    ap.add_argument('--mutant', default=None, help='canary: corruption of one contract (see zmutants.py)')
    a = ap.parse_args(argv)
    if a.mutant:
        from .zmutants import apply
        apply(a.mutant)
    n, p = a.n, eval(a.p, {}, {})
    orig = {int(s) for s in a.orig.split(',') if s}
    versions = {k: k not in orig for k in range(1, 7)}
    lem = [l for l in LEMMAS if l not in a.drop.split(',')]
    extra = None
    if a.mode == 'E3':
        extra = lambda F: [f'(not {g})' for g in groups(F)]
    elif a.mode == 'E3fail':
        extra = lambda F: [groups(F)[a.group - 1]]
    wiring, wiring_src = None, 'reference'
    if a.wiring:
        live = json.load(open(a.wiring))
        wiring, wiring_src = {k: live[k] for k in ('order', 'steps', 'out')}, 'live ECAdd.decompose_bloq()'
    if a.wmut:
        import copy
        wiring = copy.deepcopy(wiring or WIRING)
        WMUT[a.wmut](wiring)
    key = lambda w: (w['order'], w['steps'], w['out'])
    same_as_ref = wiring is None or key(wiring) == key(WIRING)
    o, obl = build(n, p, versions, lem, extra, wiring)
    if a.mode == 'E3fail':
        flagtags = {'pre:step5#5', 'pre:step6#4', 'pre:step6#5', 'pre:step6#6', 'pre:step6#7', 'pre:step6#8'}
        bad = '(or ' + ' '.join(f'(not {t})' for tag, t in obl if tag in flagtags) + ')'
        checks = [('always_fails', bad), ('nonvacuous', 'false')]
    else:
        checks = obl
    ans, sec, err, txt = solve(o.lines, checks)
    if a.smt:
        open(a.smt, 'w').write(txt)
    rows = [dict(tag=t, z3=r) for (t, _), r in zip(checks, ans)]
    if a.mode == 'E3fail':
        ok = len(ans) == 2 and ans[0] == 'unsat' and ans[1] == 'sat'
        status = 'PROVED_ALWAYS_FAILS' if ok else 'NOT_PROVED'
    else:
        ok = len(ans) == len(checks) and all(r == 'unsat' for r in ans)
        status = 'PROVED' if ok else ('REFUTED' if 'sat' in ans else 'UNKNOWN')
    if not same_as_ref and status == 'PROVED':
        status = 'NOT_PROVED'   # the wiring read from Qualtran differs from the reference: refuse
    res = dict(wiring_source=wiring_src, wiring_equals_reference=same_as_ref, n=n, p=str(p), mode=a.mode, versions=versions, lemmas=lem, group=a.group, status=status, mutant=a.mutant, wmut=a.wmut, wiring_checked_in=[c['file'] for c in WIRING['checked_in']],
               z3_version=_z3.version(), checks=len(checks), answers={r: ans.count(r) for r in set(ans)}, z3_seconds=sec,
               failed=[r['tag'] for r in rows if r['z3'] != 'unsat'][:10], stderr=err)
    print(json.dumps(res, ensure_ascii=False))
    if a.out:
        json.dump(dict(res, rows=rows), open(a.out, 'w'), ensure_ascii=False, indent=1)


if __name__ == '__main__':
    main()
