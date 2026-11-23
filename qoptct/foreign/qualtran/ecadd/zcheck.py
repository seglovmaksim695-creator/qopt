"""Cross-check of the three records of the step specifications at sample points: specs.py (ROBDD, used in E1), zspecs.py
(IntOps; the same record gives the Z3 glue formulas), ispecs.py.  Points: random, rare values, chains on curves."""
import argparse, json, random, time
from qoptct.verify.bdd import BDD, ONE, ZERO
from .specs import STEPS
from .ispecs import ISTEPS
from .zspecs import ZSTEPS, IntOps
from . import ec
from .indep import NAMES

BOX = {1: [], 2: [('KINV', 'x'), ('MM', 'target')], 3: [('MM', 'target')], 4: [('MM', 'target'), ('MM', 'target')],
       5: [('KINV', 'x'), ('MM', 'target')], 6: []}


def bdd_eval(k, n, p, patched, v, ovals):
    B = BDD()
    W, asg = {}, {}
    def word(val, w):
        # constant edges: the same code of specs.py but without BDD growth (constant folding): fast at n=256
        return [ONE if (val >> i) & 1 else ZERO for i in range(w)]
    for nm in NAMES[k]:
        W[nm] = word(int(v[nm]), 1 if nm in ('f1', 'f2', 'f3', 'f4', 'ctrl') else n)
    O = [{key: word(int(ov), n)} for (kind, key), ov in zip(BOX[k], ovals)]
    S = STEPS[k](n, p, patched=patched)
    if hasattr(S, 'definitions'):
        S.N = S.definitions(B, W)          # named words of step 3, by their definitions
    def val(w):
        assert all(e in (ONE, ZERO) for e in w)
        return sum(int(e == ONE) << i for i, e in enumerate(w))
    pre = S.pre(B, W, O)
    pre = [bool(e == ONE) for e in (pre if isinstance(pre, list) else [pre])]
    calls = [(kind, {r: val(w) for r, w in args.items()}) for kind, args in S.calls(B, W, O)]
    post = {nm: val(w) for nm, w in S.post(B, W, O).items()}
    return pre, calls, post


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True)
    ap.add_argument('--count', type=int, default=200)
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    n, p = a.n, eval(a.p, {}, {})
    rnd = random.Random(a.seed)
    io = IntOps(n, p)
    t0 = time.time()
    res = {}
    special = [0, 1, p - 1, p - 2, (p - 1) // 2]
    for patched in (True, False):
        for k in range(1, 7):
            bad, cnt = [], 0
            for t in range(a.count):
                v = {nm: (rnd.randrange(2) if nm in ('f1', 'f2', 'f3', 'f4', 'ctrl') else
                          (rnd.choice(special) if rnd.random() < 0.3 else rnd.randrange(p))) for nm in NAMES[k]}
                if k == 6 and rnd.random() < 0.3:   # some inputs of step 6 with x=y=0 / x=a etc.
                    v['x'], v['y'] = rnd.choice([(0, 0), (v['a'], v['b']), (v['a'], (p - v['b']) % p)])
                if k == 5 and rnd.random() < 0.5:   # some with P5 satisfied
                    pass
                ovals = [rnd.choice([0, rnd.randrange(p)]) for _ in BOX[k]]
                if k == 5 and rnd.random() < 0.5:
                    K, Z = ovals
                    v['lam'] = v['ctrl'] * (Z ^ (v['lam_r'] if K == 0 else 0))
                zp, zc, zpost = ZSTEPS[k](io, v, ovals, patched)
                ip = ISTEPS[k](n, p, patched)
                o_i = [{key: ov} for (kind_, key), ov in zip(BOX[k], ovals)]
                ipre, ipost = ip.pre(v, o_i), ip.post(v, o_i)
                bp, bc, bpost = bdd_eval(k, n, p, patched, v, ovals)
                zpost_i = {kk: int(vv) for kk, vv in zpost.items()}
                zc_i = [(kind_, {r: int(w) for r, w in args.items()}) for kind_, args in zc]
                ok = ([bool(x) for x in zp] == bp and zc_i == bc and zpost_i == bpost
                      and bool(ipre) == all(bp) and (not ipre or ipost == {kk: zpost_i[kk] for kk in ipost}))
                cnt += 1
                if not ok and len(bad) < 3:
                    bad.append(dict(v={kk: str(vv) for kk, vv in v.items()}, ovals=[str(x) for x in ovals]))
            res[f'{"patched" if patched else "orig"}:step{k}'] = dict(points=cnt, mismatches=len(bad), examples=bad)
            print(k, patched, cnt, len(bad), round(time.time() - t0, 1), flush=True)
    print(json.dumps({kk: (vv['points'], vv['mismatches']) for kk, vv in res.items()}))
    if a.out:
        json.dump(dict(n=n, p=str(p), seconds=round(time.time() - t0, 1), steps=res), open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
