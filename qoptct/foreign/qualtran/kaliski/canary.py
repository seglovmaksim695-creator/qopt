"""canary -- damaged variants of the REAL composed circuit (one operation replaced).

A canary counts only if (1) the prover (the same QChain client, the same Walker) does NOT give CHAIN_PROVED and
(2) a real witness is found: an input x on which the damaged circuit, in the independent simulation (simulate_flat), violates the
condition of a gadget, leaves a dirty wire, or gives an output/junk different from the specification.
    python -m qoptct.foreign.qualtran.kaliski.canary --n N --p P --out FILE [--only NAME] [--random K]
"""
import argparse, copy, hashlib, json, random, sys, time
from array import array

from qoptct.verify.compose.lift import K_X, K_CX, K_AND, K_ERASE
from qoptct.verify.compose.walker import simulate_flat

from .lower import Composed, OpsView
from .chain import QChain, ProofFailure, build_marks
from . import pq as qkal_pq
from .run import samples, expected, qualtran_classical, word


class PatchedOps(OpsView):
    def __init__(self, C, patches):
        OpsView.__init__(self, C)
        self.patches = dict(patches)

    def part_array(self, j):
        arr = OpsView.part_array(self, j)
        g0, nops = self.C.parts[j][0], self.C.parts[j][1]
        for i, op in self.patches.items():
            if g0 <= i < g0 + nops:
                arr[5 * (i - g0):5 * (i - g0) + 5] = array('i', op)
        return arr


def mutant(C, patches):
    M = copy.copy(C)
    M.ops = PatchedOps(C, patches)
    return M


def op_at(C, i):
    return tuple(C.ops[5 * i + f] for f in range(5))


def seg_ranges(marks):
    out, lo, it = {}, 0, -1
    for pos, lab, info in marks:
        if lab == 'S1':
            it += 1
        out[(it if lab.startswith('S') else None, lab)] = (lo, pos)
        lo = pos
    return out


def find(C, lo, hi, pred):
    return [i for i in range(lo, hi) if pred(op_at(C, i))]


def witness(M, n, p, xs):
    """simulate the mutant; -> (x, reason) or None."""
    X = M.inputs['x']
    init = {w: sum(1 << j for j, x in enumerate(xs) if (x >> i) & 1) for i, w in enumerate(X)}
    vals, err = simulate_flat(M.lifted, init, len(xs))
    if err is not None:
        mk = err[2]
        j = (mk & -mk).bit_length() - 1
        return xs[j], 'op %d: %s' % (err[0], err[1])
    Xo, J = M.outputs['x'], M.outputs['junk']
    keep = set(Xo) | set(J)
    ref = qualtran_classical(n, p, xs)
    for j, x in enumerate(xs):
        o, jk = word(vals, Xo, j), word(vals, J, j)
        if o != expected(n, p, x):
            return x, 'wrong output x^-1*2^(2N) mod p'
        if ref[x] != (o, jk):
            return x, 'junk (history bits) differs from Qualtran'
    dirty = 0
    for w in range(M.nq):
        if w not in keep and vals[w]:
            dirty |= vals[w]
    if dirty:
        j = (dirty & -dirty).bit_length() - 1
        return xs[j], 'a wire outside x/junk is not zero at the end'
    return None


def detection_iter(n, p, x):
    out, junk, states, bad = qkal_pq.run(n, p, x)
    ts = [(junk >> (2 * n - 1 - k)) & 1 for k in range(2 * n)]
    return ts.index(1)


def targeted(C, marks, n, p):
    """the targeted canaries: (name, what, patches, extra inputs)."""
    R = seg_ranges(marks)
    scr = C.scratch
    info = {(k if lab.startswith('S') else None, lab): inf for k, (pos, lab, inf) in
            zip(_iters(marks), marks)}
    cands = []
    nop = lambda t: (K_CX, scr, t, 0, -1)              # CX from the clean scratch wire: deletes a gate
    # T1: first iteration, S6: CX(s0 -> a) retargeted to b
    lo, hi = R[(0, 'S6')]
    i = hi - 2
    op = op_at(C, i)
    a, b = info[(0, 'S6')]['ab']
    cands.append(('T1 joint: first iteration, clearing of a', 'iteration 0, S6, op %d %r: CX target a -> b' % (i, op),
                  {i: (K_CX, op[1], b if op[2] == a else a, 0, -1)}, []))
    # T2: last iteration, S6: inside ModDbl -- the CX that clears the sign bit (X, CX, X on lower_bit) retargeted
    k = 2 * n - 1
    lo, hi = R[(k, 'S6')]
    xs_ = [j for j in range(lo + 3, hi - 3) if op_at(C, j)[0] == K_X and op_at(C, j + 2)[0] == K_X and
           op_at(C, j + 1)[0] == K_CX and op_at(C, j)[1] == op_at(C, j + 1)[1] == op_at(C, j + 2)[1]]
    j = xs_[-1] + 1
    op = op_at(C, j)
    cands.append(('T2 joint: last iteration, clearing of the sign in ModDbl', 'iteration %d, S6, op %d %r: target -> u0' % (k, j, op),
                  {j: (K_CX, op[1], info[(k, 'S6')]['ins']['u'][0], 0, -1)}, []))
    # T3, T4: the detection iteration of x0 = p-1, S1: CX(m -> f) and CX(m -> t) deleted
    x0 = p - 1
    k0 = detection_iter(n, p, x0)
    lo, hi = R[(k0, 'S1')]
    m, f, t = (info[(k0, 'S1')]['ins'][r][0] for r in ('m', 'f', 'is_terminal'))
    jf = find(C, lo, hi, lambda o: o[0] == K_CX and o[1] == m and o[2] == f)[0]
    jt = find(C, lo, hi, lambda o: o[0] == K_CX and o[1] == m and o[2] == t)[0]
    cands.append(('T3 iteration that detects v=0 (x=p-1): flag f', 'iteration %d, S1, op %d: CX(m->f) deleted' % (k0, jf),
                  {jf: nop(f)}, [x0]))
    cands.append(('T4 iteration that detects v=0 (x=p-1): history bit t', 'iteration %d, S1, op %d: CX(m->t) deleted' % (k0, jt),
                  {jt: nop(t)}, [x0]))
    # T5: iteration 3, S2: CX(m -> b) deleted
    lo, hi = R[(3, 'S2')]
    m, b = info[(3, 'S2')]['ins']['m'][0], info[(3, 'S2')]['ins']['b'][0]
    j = find(C, lo, hi, lambda o: o[0] == K_CX and o[1] == m and o[2] == b)[0]
    cands.append(('T5 flag m -> b (iteration 3)', 'iteration 3, S2, op %d: CX(m->b) deleted' % j, {j: nop(b)}, []))
    # T6: iteration 5, S3: CX(ctrl -> a) deleted
    lo, hi = R[(5, 'S3')]
    a = info[(5, 'S3')]['ins']['a'][0]
    j = find(C, lo, hi, lambda o: o[0] == K_CX and o[2] == a)[0]
    cands.append(('T6 flag a in step 3 (iteration 5)', 'iteration 5, S3, op %d: CX(c->a) deleted' % j, {j: nop(a)}, []))
    # T7, T8: measurement-based erasures: the correction wire changed
    for nm, (kk, lab) in (('T7 erase correction: MultiAnd of step 1 (iteration 1)', (1, 'S1')),
                          ('T8 erase correction: CAdd of step 5 (iteration 10)', (10, 'S5'))):
        lo, hi = R[(kk, lab)]
        js = find(C, lo, hi, lambda o: o[0] == K_ERASE)
        j = js[len(js) // 2]
        op = op_at(C, j)
        other = info[(kk, lab)]['ins']['f'][0]
        if other in op[1:4]:
            other = info[(kk, lab)]['ins']['v'][1]
        cands.append((nm, 'iteration %d, %s, op %d %r: correction wire %d -> %d' % (kk, lab, j, op, op[2], other),
                      {j: (K_ERASE, op[1], other, op[3], op[4])}, []))
    # T8b: (T8 turned out EQUIVALENT at N=256: the walker proved it and no witness exists) -- the FIRST erasure of S5 in
    # iteration 10 with the other correction wire replaced by a bit of u
    lo, hi = R[(10, 'S5')]
    j = find(C, lo, hi, lambda o: o[0] == K_ERASE)[0]
    op = op_at(C, j)
    other = [w for w in info[(10, 'S5')]['ins']['u'] if w not in op[1:4]][1]
    cands.append(('T8b erase correction: first erase of step 5 (iteration 10)', 'iteration 10, S5, op %d %r: correction wire %d -> %d'
                  % (j, op, op[1], other), {j: (K_ERASE, other, op[2], op[3], op[4])}, []))
    # T9: the last reduction AddK(p+1): one X of the constant loader deleted
    flo, fhi = R[(None, 'fin')]
    js = find(C, flo, fhi, lambda o: o[0] == K_X)
    rset = set(info[(None, 'fin')]['r'])
    js = [j for j in js if op_at(C, j)[1] not in rset and op_at(C, j)[1] not in info[(None, 'fin')]['u']
          and op_at(C, j)[1] not in info[(None, 'fin')]['s']]
    j = js[0]
    cands.append(('T9 last reduction AddK(p+1): constant loader', 'fin, op %d %r: X deleted' % (j, op_at(C, j)),
                  {j: nop(op_at(C, j)[1])}, []))
    # T10: the x = 0 branch: one CX(t0 -> r_i) of C[XorK(p)] deleted
    t0 = info[(None, 'fin')]['t0']
    js = find(C, flo, fhi, lambda o: o[0] == K_CX and o[1] == t0 and o[2] in rset)
    j = js[0]
    cands.append(('T10 branch x=0: C[XorK(p)] on t0', 'fin, op %d %r: CX(t0->r) deleted' % (j, op_at(C, j)),
                  {j: nop(op_at(C, j)[2])}, [0]))
    # T11: last iteration, S6 ModDbl: the controlled +p correction -- one CX deleted
    lo, hi = R[(2 * n - 1, 'S6')]
    js = find(C, lo, hi, lambda o: o[0] == K_CX)
    j = js[len(js) // 2]
    cands.append(('T11 last iteration: middle of step 6 (ModDbl)', 'iteration %d, S6, op %d %r: CX deleted'
                  % (2 * n - 1, j, op_at(C, j)), {j: nop(op_at(C, j)[2])}, []))
    # T13: carry through all bits: the LAST AND of the final adder AddK(p+1) (the highest carry) gets another control
    js = find(C, flo, fhi, lambda o: o[0] == K_AND)
    j = js[-1]
    op = op_at(C, j)
    c = [w for w in info[(None, 'fin')]['u'] if w not in op[1:4]][0]
    cands.append(('T13 carry through all bits: top carry of the adder AddK(p+1)', 'fin, op %d %r: AND control %d -> %d'
                  % (j, op, op[1], c), {j: (K_AND, c, op[2], op[3], -1)}, []))
    # T12: the loader u = p: one X deleted
    ilo, ihi = R[(None, 'init')]
    useti = set(info[(None, 'init')]['u'])
    j = find(C, ilo, ihi, lambda o: o[0] == K_X and o[1] in useti)[0]
    cands.append(('T12 joint: loader u = p', 'init, op %d %r: X deleted' % (j, op_at(C, j)), {j: nop(op_at(C, j)[1])}, []))
    return cands


def _iters(marks):
    it = -1
    for pos, lab, info in marks:
        if lab == 'S1':
            it += 1
        yield it


def random_canaries(C, marks, n, p, K, seed=20261124):
    rnd = random.Random(seed)
    R = seg_ranges(marks)
    keys = [key for key in R if key[0] is not None]
    out = []
    while len(out) < K:
        key = rnd.choice(keys)
        lo, hi = R[key]
        j = rnd.randrange(lo, hi)
        op = op_at(C, j)
        kind = op[0]
        info = [inf for (pos, lab, inf), it in zip(marks, _iters(marks)) if (it, lab) == key][0]
        regw = [w for ws in info['ins'].values() for w in ws]
        if kind == K_X:
            new = (K_CX, C.scratch, op[1], 0, -1)
            what = 'X deleted'
        elif kind == K_CX:
            if rnd.random() < 0.5:
                new = (K_CX, op[2], op[1], 0, -1)
                what = 'CX control/target swapped'
            else:
                c = rnd.choice([w for w in regw if w not in (op[1], op[2])])
                new = (K_CX, op[1], c, 0, -1)
                what = 'CX target -> %d' % c
        elif kind == K_AND:
            c = rnd.choice([w for w in regw if w not in op[1:4]])
            new = (K_AND, c, op[2], op[3], -1)
            what = 'AND control -> %d' % c
        else:
            c = rnd.choice([w for w in regw if w not in op[1:4]])
            new = (K_ERASE, op[1], c, op[3], op[4])
            what = 'ERASE correction -> %d' % c
        out.append(('R%02d random: iteration %d, %s' % (len(out) + 1, key[0], key[1]),
                    'op %d %r: %s' % (j, op, what), {j: new}, []))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--only', default='')
    ap.add_argument('--random', type=int, default=10)
    ap.add_argument('--no-prove', action='store_true')
    a = ap.parse_args(argv)
    n, p = a.n, a.p
    t0 = time.time()
    C = Composed(n, p)
    marks = build_marks(C)
    xs = samples(p, n, 127)
    cands = targeted(C, marks, n, p) + random_canaries(C, marks, n, p, a.random)
    from .run import manifest
    res = dict(n=n, p=p, started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), canaries=[], manifest=manifest())
    for name, what, patches, extra in cands:
        if a.only and a.only not in name:
            continue
        t = time.time()
        M = mutant(C, patches)
        wit = witness(M, n, p, sorted(set(xs + extra)))
        row = dict(name=name, what=what, patches={str(k): v for k, v in patches.items()},
                   witness=None if wit is None else dict(x=wit[0], reason=wit[1]), sim_sec=round(time.time() - t, 1))
        if not a.no_prove:
            t = time.time()
            ch = QChain(M, lifted=M.lifted, marks=marks)
            try:
                ch.run()
                row['prove'] = 'CHAIN_PROVED'
            except ProofFailure as e:
                row['prove'] = 'NOT_PROVED'
                row['failure'] = str(e)[:300]
            row['prove_sec'] = round(time.time() - t, 1)
        row['counted'] = row.get('prove') == 'NOT_PROVED' and row['witness'] is not None
        res['canaries'].append(row)
        print(json.dumps(row, ensure_ascii=False)[:700], flush=True)
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1, ensure_ascii=False)
    res['counted'] = sum(r['counted'] for r in res['canaries'])
    res['targeted_counted'] = sum(r['counted'] for r in res['canaries'] if r['name'].startswith('T'))
    res['sec'] = round(time.time() - t0, 1)
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1, ensure_ascii=False)
    print('counted %d (targeted %d) of %d; %.0fs' % (res['counted'], res['targeted_counted'], len(res['canaries']), res['sec']))


if __name__ == '__main__':
    main()
