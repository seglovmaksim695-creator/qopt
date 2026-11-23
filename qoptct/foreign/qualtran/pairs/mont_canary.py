"""Canaries for the mirror check of the Montgomery pair (w=1): corrupt ONE element of the backward program and
require that the check does NOT pass.  python -m qoptct.foreign.qualtran.pairs canary-montgomery --n N --w 1 --p P"""
import sys, json, random, time, gc
from . import pair_compose as PC
from . import pair_canary as PK
from . import mont_compose as MC
from . import mirror as M
from .mirror import K_X, K_CX, K_AND, K_ERASE, K_LD, K_CLR


def token_canaries(ops, rng):
    out = []
    def pick(kind):
        return rng.choice([i for i, o in enumerate(ops) if o[0] == kind])
    for _ in range(2):
        i = pick(K_LD); k, a, b, c, l = ops[i]; o = ops[:]; o[i] = (k, a, b, c ^ 1 if c else 1, l)
        out.append(('token: code LD #%d changed %d->%d' % (i, c, o[i][3]), 'FAIL', o))
    for _ in range(2):
        i = pick(K_CLR); k, a, b, c, l = ops[i]; o = ops[:]; o[i] = (k, a, b, c ^ 1 if c else 1, l)
        out.append(('token: code CLR #%d changed %d->%d' % (i, c, o[i][3]), 'FAIL', o))
    i = pick(K_LD); o = ops[:]; del o[i]; out.append(('token: drop LD #%d' % i, 'FAIL', o))
    i = pick(K_CLR); o = ops[:]; del o[i]; out.append(('token: drop CLR #%d' % i, 'FAIL', o))
    i = pick(K_LD); k, a, b, c, l = ops[i]; o = ops[:]; o[i] = (K_CLR, a, b, c, l)
    out.append(('token: LD #%d replace by CLR (break the mirror)' % i, 'FAIL', o))
    return out


def macro_canaries2(C, rng):
    """like pair_canary.macro_canaries, but the junk canary takes a wire that REALLY differs between two neighbouring nodes (otherwise the mutation is empty)"""
    out = [c for c in PK.macro_canaries(C, rng) if 'junk' not in c[0]]
    bidx = [i for i, m in enumerate(C.macros) if m['kind'] == 'b']
    Pb = C.Pb
    cnt = 0
    for _ in range(200):
        k = rng.randrange(len(bidx) - 1); j1, j2 = bidx[k], bidx[k + 1]
        reg = rng.choice(list(Pb.ins)); w = rng.choice(Pb.ins[reg])
        if C.macros[j1]['g'][w] == C.macros[j2]['g'][w]:
            continue
        def f(ms, j1=j1, j2=j2, w=w):
            a, b = dict(ms[j1]), dict(ms[j2]); ga, gb = dict(a['g']), dict(b['g'])
            ga[w], gb[w] = gb[w], ga[w]; a['g'], b['g'] = ga, gb; ms[j1], ms[j2] = a, b
        ms = list(C.macros); f(ms)
        out.append(('macro: junk, nodes #%d,#%d: swap local wire %d (reg. %s)' % (j1, j2, w, reg), 'FAIL', ms))
        cnt += 1
        if cnt == 3:
            break
    return out


def main(n, w, p, seed=1, outf=None, only=None):
    """Run the canaries of the Montgomery pair.  ``only``: set of job indexes (default all)."""
    rng = random.Random(seed); t0 = time.time()
    C = MC.MontPair(n, w, p)
    base = PC.check_composed(C)
    res = dict(n=n, p=p, seed=seed, base_mirror_ok=base['mirror_ok'], rows=[])
    print('base', base['mirror_ok'], flush=True)
    piece_res = PC.piece_pi(C.Pf, C.Pb)[0]
    jobs = []
    ops = C.Pb.ops
    for nm, ex, newops in PK.piece_canaries(ops, rng, C.Pb.nsrc) + token_canaries(ops, rng):
        jobs.append((nm, ex, dict(Pb=PK.clone_piece(C.Pb, newops))))
    for nm, ex, ms in macro_canaries2(C, rng):
        jobs.append((nm, ex, dict(macros=ms, piece_res=piece_res)))
    for nm, ex, t in PK.top_canaries(C, rng):
        jobs.append((nm, ex, dict(T=t, piece_res=piece_res)))
    res['jobs_total'] = len(jobs)
    for idx, (nm, ex, kw) in enumerate(jobs):
        if only is not None and idx not in only:
            continue
        gc.collect(); t1 = time.time()
        try:
            r = PC.check_composed(C, **kw); ok = bool(r['mirror_ok'])
            why = (r.get('error') or r.get('piece_msg')) if not ok else 'ok'
            if not ok and not r.get('error') and r.get('piece_mirror_ok') and not r.get('top_ok', True):
                why = [m for m in r.get('top_match', []) if not m[1]]
        except Exception as e:
            ok, why = False, 'exception %s: %s' % (type(e).__name__, str(e)[:200])
        row = dict(index=idx, name=nm, expect=ex, verdict='PASS' if ok else 'NOT PASS', why=str(why)[:300], sec=round(time.time() - t1, 1))
        res['rows'].append(row); print(json.dumps(row, ensure_ascii=False), flush=True)
    nf = [r for r in res['rows'] if r['expect'] == 'FAIL']
    res['n_fail_canaries'] = len(nf); res['n_caught'] = sum(1 for r in nf if r['verdict'] != 'PASS')
    res['false_pass'] = [r['name'] for r in nf if r['verdict'] == 'PASS']
    res['equiv'] = [(r['name'], r['verdict']) for r in res['rows'] if r['expect'] == 'EQUIV']
    res['sec'] = round(time.time() - t0, 1)
    print(json.dumps({k: v for k, v in res.items() if k != 'rows'}, ensure_ascii=False), flush=True)
    if outf:
        json.dump(res, open(outf, 'w'), ensure_ascii=False, indent=1)
    return res
