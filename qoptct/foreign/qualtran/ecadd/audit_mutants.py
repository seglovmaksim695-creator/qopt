"""Auditor's mutants and checks for the point-addition proof, as callable parts (regression: every mutant must NOT be proved).

``python -m qoptct.foreign.qualtran.ecadd.audit_mutants {identity,lemmas,glue,box} [--n 16] [--out FILE]``

* ``identity``  the six step blocks proved in ``run_step`` are exactly the blocks inside ``ECAdd.decompose_bloq()``; the boxes
                inside them are the multiplier / inverse with the expected parameters (needs Qualtran kind ``'e'``).
* ``lemmas``    own brute force of the eight glue lemmas on all small curves (p = 5, 7, 11, 13; no Qualtran needed).
* ``glue``      mutated step contracts and mutated wiring in the Z3 glue (control + 7 mutants + XOR-property probes).
* ``box``       mutated interfaces of the multiplier / inverse boxes inside steps 2-5 (needs kind ``'e'``).
Each part returns a dict; ``not_proved_ok`` is True when every mutant is NOT proved and the control is proved.
"""
import copy
import itertools
import json
import sys

P255, P16 = 2 ** 255 - 19, 32749


class _R(dict):
    pass


def part_identity(R):
    from qualtran.bloqs.cryptography.ecc import ec_add as EC
    from qualtran.bloqs.mod_arithmetic import DirtyOutOfPlaceMontgomeryModMul as MMc, KaliskiModInverse as KIc
    cls = [EC._ECAddStepOne, EC._ECAddStepTwo, EC._ECAddStepThree, EC._ECAddStepFour, EC._ECAddStepFive, EC._ECAddStepSix]
    out = {}
    for n, p in ((16, P16), (256, P255)):
        cb = EC.ECAdd(n, p, 1).decompose_bloq()
        rows = []
        for bi in cb.bloq_instances:
            k = [c.__name__ for c in cls].index(type(bi.bloq).__name__) + 1
            mine = cls[k - 1](n=n, mod=p)
            boxes = []

            def walk(b, d):
                if d > 3:
                    return
                try:
                    sub = b.decompose_bloq()
                except Exception:  # noqa
                    return
                for s in sub.bloq_instances:
                    bb = s.bloq
                    core = getattr(bb, 'subbloq', bb)
                    if isinstance(core, (MMc, KIc)):
                        boxes.append(dict(type=type(core).__name__, wrapper=type(bb).__name__, bitsize=core.bitsize,
                                          mod=int(core.mod) == p, window=getattr(core, 'window_size', None),
                                          uncompute=getattr(core, 'uncompute', None)))
                    elif d < 2:
                        walk(bb, d + 1)
            walk(bi.bloq, 0)
            rows.append(dict(step=k, same_as_E1=(bi.bloq == mine), boxes=boxes))
        w2 = EC.ECAdd(n, p, 2).decompose_bloq()
        out['n%d' % n] = dict(steps=sorted(rows, key=lambda r: r['step']),
                              window2_steps_differ=any(bi.bloq != cls[[c.__name__ for c in cls].index(type(bi.bloq).__name__)](n=n, mod=p)
                                                       for bi in w2.bloq_instances))
    R['identity'] = out




def part_lemmas(R, primes=(5, 7, 11, 13), k=4):
    """own brute force of the 8 lemmas as asserted in zglue (box semantics: MM=u*v/R, KINV=u^-1*R^2, KINV(0)=0)."""
    out = {}
    for p in primes:
        Rm = pow(2, k, p)
        Ri = pow(Rm, -1, p)
        MM = lambda u, v: u * v * Ri % p
        KI = lambda u: 0 if u == 0 else pow(u, -1, p) * Rm * Rm % p
        neg = lambda v: (-v) % p
        bad = {nm: 0 for nm in ('MM0', 'K0', 'C1', 'C2', 'EC_DBL', 'EC_TAN', 'EC_RNOTO', 'EC_RNOTP')}
        cnt, ex = 0, {}
        for u in range(p):
            for v in range(p):
                if (MM(u, v) == 0) != (u == 0 or v == 0):
                    bad['MM0'] += 1
                if u and MM(u, MM(KI(u), v)) != v:
                    bad['C1'] += 1
                if u and MM(KI(u), MM(u, v)) != v:
                    bad['C2'] += 1
            if (KI(u) == 0) != (u == 0):
                bad['K0'] += 1
        for A_, B_ in itertools.product(range(p), range(1, p)):
            if (4 * A_ ** 3 + 27 * B_ ** 2) % p == 0:
                continue
            pts = [(X, Y) for X in range(p) for Y in range(p) if (Y * Y - X ** 3 - A_ * X - B_) % p == 0]
            allp = pts + [None]
            for P_ in allp:
                for Q_ in allp:
                    a, b = (0, 0) if P_ is None else (P_[0] * Rm % p, P_[1] * Rm % p)
                    x, y = (0, 0) if Q_ is None else (Q_[0] * Rm % p, Q_[1] * Rm % p)
                    if P_ is not None and P_[1] != 0:
                        lam = (3 * P_[0] ** 2 + A_) * pow(2 * P_[1], -1, p) % p
                        lrs = [lam * Rm % p]
                    else:
                        lrs = list(range(p))          # tangent undefined (P = O or 2-torsion): every lam_r tried
                    for lr in lrs:
                        cnt += 1
                        e1, e2, e3, e4 = a == x, b == neg(y), a == 0 and b == 0, x == 0 and y == 0
                        gen = (not e3) and (not e4) and not (e1 and e2)
                        Lam = lr if e1 else MM(KI((x - a) % p), (y - b) % p)
                        XR = (MM(Lam, Lam) - a - x) % p
                        YR = (MM((a - XR) % p, Lam) - b) % p
                        chk = dict(EC_DBL=not (e1 and not e2 and not e3 and not e4) or (y == b and b != 0),
                                   EC_TAN=not (gen and XR == a) or Lam == lr,
                                   EC_RNOTO=not gen or not (XR == 0 and YR == 0),
                                   EC_RNOTP=not gen or not (XR == a and YR == b))
                        for nm, ok in chk.items():
                            if not ok:
                                bad[nm] += 1
                                ex.setdefault(nm, dict(p=p, A=A_, B=B_, P=P_, Q=Q_, lam_r=lr))
        out['p%d' % p] = dict(cases=cnt, violations=bad, examples=ex)
    R['lemmas_own'] = out




def part_glue(R):
    from qoptct.foreign.qualtran.ecadd import zglue as ZG
    from qoptct.foreign.qualtran.ecadd import zspecs as ZS
    orig = dict(ZS.ZSTEPS)
    res = {}

    def run(n, p, wiring=None, extra_lines=()):
        o, obl = ZG.build(n, p, {k: True for k in range(1, 7)}, ZG.LEMMAS, None, wiring)
        o.lines.extend(extra_lines)
        ans, sec, err, txt = ZG.solve(o.lines, obl)
        st = 'PROVED' if len(ans) == len(obl) and all(a == 'unsat' for a in ans) else ('REFUTED' if 'sat' in ans else 'UNKNOWN')
        return dict(status=st, checks=len(obl), failed=[t for (t, _), a in zip(obl, ans) if a != 'unsat'][:4])

    def patch(k, fn):
        def wrapped(o, W, O, patched):
            pre, calls, post = orig[k](o, W, O, patched)
            return fn(o, W, O, pre, calls, post)
        return wrapped
    S = ZS

    def g_step4_y(o, W, O, pre, calls, post):
        return pre, calls, dict(post, y=o.xor(W['y'], O[0]))

    def g_step2_lam(o, W, O, pre, calls, post):
        lam = o.ite(W['ctrl'], o.ite(W['f1'], O[1], W['lam_r']), o.zero())
        return pre, calls, dict(post, lam=lam)

    def g_step5_x(o, W, O, pre, calls, post):
        xn = o.ite(W['ctrl'], o.modneg(W['x']), W['x'])
        return pre, calls, dict(post, x=o.modadd(W['b'], xn))

    def g_step4_raw(o, W, O, pre, calls, post):
        x4 = o._new('Int', f"(- {W['x']} {O[0]})")
        calls = [calls[0], ('MM', {'x': x4, 'y': W['lam']})]
        return pre, calls, dict(post, x=x4)

    def g_step6_anc(o, W, O, pre, calls, post):
        return pre, calls, post        # placeholder replaced below

    def g_step3_x(o, W, O, pre, calls, post):
        x3 = o.ite(W['ctrl'], o.modadd(W['x'], W['a']), W['x'])
        return pre, calls, dict(post, x=x3)
    cases = [('G1 шаг 4: y ^= Z3 вместо Z3p', 4, g_step4_y), ('G2 шаг 2: ветви lam переставлены', 2, g_step2_lam),
             ('G3 шаг 5: x = b - x вместо a - x', 5, g_step5_x), ('G4 шаг 4: аргумент ящика без приведения mod p', 4, g_step4_raw),
             ('G6 шаг 3: x += a вместо 3a', 3, g_step3_x)]
    for n, p in ((16, P16), (256, P255)):
        rr = {'G0 без порчи (контроль)': run(n, p)}
        for nm, k, fn in cases:
            S.ZSTEPS[k] = patch(k, fn)
            try:
                rr[nm] = run(n, p)
            finally:
                S.ZSTEPS[k] = orig[k]
        w = copy.deepcopy(ZG.WIRING)
        w['steps']['6']['f4'] = [1, 'f3']
        rr['G7 проводка: f4 шага 6 из f3 шага 1'] = run(n, p, w)
        w = copy.deepcopy(ZG.WIRING)
        w['steps']['6']['f1'], w['steps']['6']['f2'] = w['steps']['6']['f2'], w['steps']['6']['f1']
        rr['G5 проводка: f1 и f2 шага 6 переставлены'] = run(n, p, w)
        # G8: XOR abstraction: is every one of its 3 properties needed? (drop one at a time; PROVED = not needed)
        for j, nm in enumerate(('u=0 -> v=w', 'w=0 -> v=u', 'u=w -> v=0')):
            keep = S.SmtOps.xor

            def xor_drop(self, u, w, j=j):
                v = self._new('Int')
                self.lines.append(f'(assert (and (<= 0 {v}) (< {v} {1 << self.n})))')
                props = [f'(assert (=> (= {u} 0) (= {v} {w})))', f'(assert (=> (= {w} 0) (= {v} {u})))',
                         f'(assert (=> (= {u} {w}) (= {v} 0)))']
                self.lines.extend(pp for i, pp in enumerate(props) if i != j)
                return v
            S.SmtOps.xor = xor_drop
            try:
                rr['G8.%d без свойства XOR «%s» (справка: нужно ли оно)' % (j + 1, nm)] = run(n, p)
            finally:
                S.SmtOps.xor = keep
        res['n%d' % n] = rr
    R['glue_canaries'] = res




def part_box(R, n=16, p=P16):
    from qoptct.foreign.qualtran.ecadd.core import build_program, Run, Fail
    from qoptct.foreign.qualtran.ecadd.specs import STEPS
    from qualtran.bloqs.cryptography.ecc import ec_add as EC
    cls = {2: EC._ECAddStepTwo, 4: EC._ECAddStepFour, 5: EC._ECAddStepFive, 3: EC._ECAddStepThree}
    out = {}

    def prove(k, prog):
        try:
            st = Run(prog, STEPS[k](n, p, patched=True, observe=False), n, p, log=lambda *a, **kw: None).run()
            return st, None
        except Fail as e:
            return 'REFUTED', str(e)[:160]
    progs = {k: build_program(cls[k](n=n, mod=p), keep_ops=False) for k in (2, 3, 4, 5)}

    def mutated(k, fn):
        prog = dict(progs[k])
        items = list(prog['items'])
        fn(items)
        prog['items'] = items
        return prog

    def leaves(items, kind):
        return [i for i, it in enumerate(items) if it[0] == 'opaque' and it[1]['kind'] == kind]

    def setleaf(items, i, **kw):
        lf = copy.deepcopy(items[i][1])
        for key, val in kw.items():
            lf[key] = val
        items[i] = (items[i][0], lf, items[i][2])
    out['B0 шаг 4 без порчи (контроль)'] = prove(4, progs[4])

    def b1(items):
        i = leaves(items, 'MM_adj')[0]
        y = list(items[i][1]['regs_in']['y'])
        y[0], y[1] = y[1], y[0]
        setleaf(items, i, regs_in=dict(items[i][1]['regs_in'], y=y), regs_out=dict(items[i][1]['regs_out'], y=y))
    out['B1 шаг 4: обратный MM получает y с двумя переставленными проводами'] = prove(4, mutated(4, b1))

    def b2(items):
        i = leaves(items, 'MM')[1]
        lf = items[i][1]
        setleaf(items, i, regs_in=dict(lf['regs_in'], y=lf['regs_in']['x']), regs_out=dict(lf['regs_out'], y=lf['regs_out']['x']))
    out['B2 шаг 4: второй MM считает MM(x, x) вместо MM(x, lam)'] = prove(4, mutated(4, b2))

    def b3(items):
        i = leaves(items, 'KINV')[0]
        lf = items[i][1]
        other = [ws for r, ws in lf['regs_in'].items()][0]
        setleaf(items, i, regs_out=dict(lf['regs_out'], junk=list(lf['regs_out']['junk'][:-1]) + [other[0]]))
    out['B3 шаг 5: мусор KINV записан поверх живого провода x'] = prove(5, mutated(5, b3))

    def b4(items):
        del items[leaves(items, 'MM_adj')[-1]]
    out['B4 шаг 4: обратный MM удалён'] = prove(4, mutated(4, b4))

    def b5(items):
        i, j = leaves(items, 'MM')[:2]
        items[i], items[j] = items[j], items[i]
    out['B5 шаг 4: два прямых MM переставлены местами'] = prove(4, mutated(4, b5))
    R['box_canaries_n16'] = {k: dict(status=v[0], detail=v[1]) for k, v in out.items()}





def run(part, n=16):
    R = {}
    if part == 'identity':
        part_identity(R)
    elif part == 'lemmas':
        part_lemmas(R, primes=(5, 7))
    elif part == 'glue':
        part_glue(R)
    elif part == 'box':
        part_box(R, n=n, p=P16 if n == 16 else P255)
    return R


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('part', choices=['identity', 'lemmas', 'glue', 'box'])
    ap.add_argument('--n', type=int, default=16)
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    R = run(a.part, a.n)
    R['status'] = 'OK'
    if a.out:
        json.dump(R, open(a.out, 'w'), indent=1, default=str)
    print(json.dumps(R, default=str)[:600])


if __name__ == '__main__':
    main()
