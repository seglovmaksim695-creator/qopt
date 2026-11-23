"""Legacy BDD glue of ECAdd from the six step contracts (specs.py) + explicit lemmas (each a separate mathematical
statement, see LEMMAS), on ROBDD at full width n.  It did not scale (nested arithmetic) and is superseded by the Z3 glue
(``zglue``); it is kept as an independent small-n cross-check and is not used by the API.

The state is words as edges.  The output of every box (MM target, KINV x) is a fresh variable; its link to the arguments is
only through the listed lemmas.  Checked: the precondition of every step (ranges, P5, P6) and the arguments of every box
(<p) hold; the result equals FORMULA (point addition in Montgomery form by the textbook formulas, branches O/-P)."""
import argparse, json, time
from qoptct.verify.bdd import BDD, ONE, ZERO
from . import bv
from .bv import eq, is_zero, ite, modsub, lt_const, all_of
from .core import disjoint
from .specs import STEPS

LEMMAS = {
    'MM0': 'MM(u,v)=0 <=> u=0 или v=0 (u,v<p; p простое; MM(u,v)=u*v*R^-1 mod p — контракт C)',
    'K0': 'KINV(u)=0 <=> u=0 (u<p; KINV(u)=u^-1*R^2 mod p при u!=0, KINV(0)=0 — контракт D)',
    'CONG': 'одинаковые аргументы => одинаковый выход (MM, KINV — функции)',
    'C1': 'MM(u, MM(KINV(u), v)) = v при u!=0',
    'C2': 'MM(KINV(u), MM(u, w)) = w при u!=0',
    'EC_DBL': 'P,Q на кривой, конечны, a=x, b!=-y => y=b и b!=0',
    'EC_TAN': 'ветвь ctrl: x_R=a => наклон Lam = lam_r (касательная в P)',
    'EC_RNOTO': 'ветвь ctrl: (x_R,y_R) != (0,0) (R конечна и лежит на кривой, (0,0) не на кривой)',
    'EC_RNOTP': 'ветвь ctrl: (x_R,y_R) != (a,b) (иначе Q=O)',
}


class GlueFail(Exception):
    def __init__(self, msg, **kw):
        super().__init__(msg)
        self.detail = kw


class Glue:
    def __init__(self, n, p, versions=None, lemmas=None, extra=None, log=print, full_lemmas=False):
        """versions: {шаг: patched?}; extra: функция (B, W, F) -> список конъюнктов (доп. предусловие, E3)."""
        self.n, self.p = n, p
        self.versions = versions or {k: True for k in range(1, 7)}
        self.lemmas = set(LEMMAS) if lemmas is None else set(lemmas)
        self.extra, self.log, self.full_lemmas = extra, log, full_lemmas
        self.stats = dict(checks=0)

    def holds(self, d):
        """d == ZERO при допущении A (конъюнкция частей). Части добавляются «слоями» по общим переменным с d:
        если d несовместно уже с частью условий — несовместно и со всем A (достаточное условие, корректно);
        если не хватило всей связной компоненты — это настоящее нарушение (остальные части от d не зависят)."""
        if d == ZERO:
            return True
        self.stats['checks_nontrivial' if 'checks_nontrivial' in self.stats else 'checks'] += 1
        B = self.B
        used, frontier = set(), B.support(d)
        while True:
            new = [i for i, (e, s_) in enumerate(self.conj) if i not in used and not s_.isdisjoint(frontier)]
            if not new:
                return False
            used.update(new)
            key = frozenset(used)
            Ak = self.cache.get(key)
            if Ak is None:
                Ak = ONE
                for i in sorted(used):
                    Ak = B.AND(Ak, self.conj[i][0])
                self.cache[key] = Ak
            if disjoint(B, Ak, d):
                return True
            frontier = set()
            for i in new:
                frontier |= self.conj[i][1]

    def add(self, e, tag):
        if e != ONE:
            self.conj.append((e, frozenset(self.B.support(e))))
            self.tags.append(tag)
            self.cache = {}

    def pop_to(self, k):
        del self.conj[k:]
        del self.tags[k:]
        self.cache = {}

    def const(self, f):
        """f — константа при текущем допущении? -> ONE/ZERO или None."""
        if f in (ONE, ZERO):
            return f
        if self.holds(f):
            return ZERO
        if self.holds(f ^ 1):
            return ONE
        return None

    def run(self, final_check='correct'):
        """Разбор случаев по флагам: в каждом случае все однобитные входы шагов — константы (доказанные при
        допущении случая), поэтому формулы шагов упрощаются. Случаи покрывают всё (каждое разбиение — f=0 и f=1).
        final_check: 'correct' — корректность; 'incorrect' — на ВСЕХ входах допущения нарушено освобождение."""
        self.t0 = time.time()
        n, p = self.n, self.p
        B = self.B = BDD()
        self.conj, self.tags, self.cache = [], [], {}
        names = ['a', 'b', 'x', 'y', 'lam_r', 'K2', 'Z2', 'Z1', 'Z3', 'Z3p', 'K5', 'Z5']
        cols = {k: [None] * n for k in names}
        for i in reversed(range(n)):
            for k in names:
                cols[k][i] = B.new_var()[1]
        self.cols, self.names = cols, names
        self.W0 = {k: cols[k] for k in names[:5]}
        self.V = {k: cols[k] for k in names[5:]}
        for k in names:
            self.add(lt_const(B, cols[k], p), f'range:{k}')
        self.S = {k: STEPS[k](n, p, patched=self.versions[k]) for k in range(1, 7)}
        a, b, x, y, lr = (self.W0[k] for k in ('a', 'b', 'x', 'y', 'lam_r'))
        self.E = dict(e1=eq(B, a, x), e2=eq(B, b, bv.modneg(B, y, p)),
                      e3=B.AND(is_zero(B, a), is_zero(B, b)), e4=B.AND(is_zero(B, x), is_zero(B, y)))
        if self.extra is not None:
            for e in self.extra(B, self.W0, self.E):
                self.add(e, 'extra-pre')
        self.final_check = final_check
        self.result = dict(n=n, p=str(p), versions=self.versions, lemmas=sorted(self.lemmas), cases=[],
                           status=None)
        self.case_rows = []
        try:
            self._cases({}, ['e1', 'e2', 'e3', 'e4'])
        except _Stop:
            pass
        r = self.result
        if r['status'] is None:
            r['status'] = 'PROVED' if final_check == 'correct' else 'PROVED_INCORRECT'
        r.update(cases=self.case_rows, seconds=round(time.time() - self.t0, 1), bdd_nodes=B.nodes(),
                 checks=self.stats['checks'])
        return r

    def _cases(self, fixed, todo):
        B = self.B
        if todo:
            k = todo[0]
            for c in (ZERO, ONE):
                mark = len(self.conj)
                self.add(B.XOR(self.E[k], c) ^ 1, f'case:{k}={int(c == ONE)}')
                self._cases(dict(fixed, **{k: c}), todo[1:])
                self.pop_to(mark)
            return
        # непустота случая (только диапазоны и предикаты случая)
        A = ONE
        for e, _ in self.conj:
            A = B.AND(A, e)
            if A == ZERO:
                break
        label = ','.join(f'{k}={int(v == ONE)}' for k, v in fixed.items())
        if A == ZERO:
            self.case_rows.append(dict(case=label, status='vacuous'))
            return
        self._chain(fixed, label)

    def _chain(self, Ec, label):
        """Шаги 1..6 при константных e1..e4 (Ec); ветвление по прочим флагам — через _flag."""
        B, n, p, S, V = self.B, self.n, self.p, self.S, self.V
        a, b, x, y, lr = (self.W0[k] for k in ('a', 'b', 'x', 'y', 'lam_r'))
        mark = len(self.conj)
        obligations, calls = [], []

        def flags(W):
            for k2, w in W.items():
                if len(w) == 1 and w[0] not in (ONE, ZERO):
                    c = self.const(w[0])
                    if c is None:
                        raise _Split(w[0])
                    W[k2] = [c]
            return W

        def step(k, W, O):
            W = flags(W)
            pre = S[k].pre(B, W, O)
            for j, e in enumerate(pre if isinstance(pre, list) else [pre]):
                obligations.append((f'pre:step{k}#{j}', e ^ 1))
            for j, (kind, args) in enumerate(S[k].calls(B, W, O)):
                for r, w in args.items():
                    obligations.append((f'arg<p:step{k}.{kind}{j}.{r}', lt_const(B, w, p) ^ 1))
                out = O[j]['x'] if kind == 'KINV' else O[j]['target']
                calls.append((f'step{k}.{kind}{j}', kind, args, out))
            return S[k].post(B, W, O)

        try:
            r1 = step(1, dict(a=a, b=b, x=x, y=y), [])
            r2 = step(2, dict(f1=r1['f1'], ctrl=r1['ctrl'], a=a, b=b, x=x, y=y, lam_r=lr),
                      [dict(x=V['K2']), dict(target=V['Z2'])])
            r3 = step(3, dict(ctrl=r2['ctrl'], a=a, b=b, x=r2['x'], y=r2['y'], lam=r2['lam']), [dict(target=V['Z1'])])
            r4 = step(4, dict(x=r3['x'], y=r3['y'], lam=r3['lam']), [dict(target=V['Z3']), dict(target=V['Z3p'])])
            r5 = step(5, dict(ctrl=r3['ctrl'], a=a, b=b, x=r4['x'], y=r4['y'], lam_r=lr, lam=r4['lam']),
                      [dict(x=V['K5']), dict(target=V['Z5'])])
            W6 = flags(dict(f1=r2['f1'], f2=r1['f2'], f3=r1['f3'], f4=r1['f4'], ctrl=r5['ctrl'], a=a, b=b,
                            x=r5['x'], y=r5['y']))
        except _Split as sp:
            for c in (ZERO, ONE):
                m2 = len(self.conj)
                self.add(B.XOR(sp.f, c) ^ 1, f'split')
                self._chain(Ec, label + f'|split={int(c == ONE)}')
                self.pop_to(m2)
            return
        r6 = step(6, W6, [])
        e1, e2, e3, e4 = (Ec[k] for k in ('e1', 'e2', 'e3', 'e4'))
        gen = all_of(B, e3 ^ 1, e4 ^ 1, B.AND(e1, e2) ^ 1)
        Lam = ite(B, e1, lr, V['Z2'])
        XR = modsub(B, modsub(B, V['Z3'], a, p), x, p)
        YR = modsub(B, V['Z3p'], b, p)
        cmap = {c[0]: c for c in calls}
        k2, m2, m4a, m4b = cmap['step2.KINV0'], cmap['step2.MM1'], cmap['step4.MM0'], cmap['step4.MM1']
        bind = all_of(B, eq(B, k2[2]['x'], modsub(B, x, a, p)), eq(B, m2[2]['x'], V['K2']),
                      eq(B, m2[2]['y'], modsub(B, y, b, p)), eq(B, m4a[2]['x'], Lam), eq(B, m4a[2]['y'], Lam),
                      eq(B, m4b[2]['x'], modsub(B, a, XR, p)), eq(B, m4b[2]['y'], Lam))
        gb = B.AND(gen, bind)
        obligations.append(('bind', B.AND(gen, bind ^ 1)))
        L = self.lemmas
        mms = [c for c in calls if c[1] == 'MM']
        kis = [c for c in calls if c[1] == 'KINV']
        if 'MM0' in L:
            for tag, _, args, out in mms:
                self.add(B.XOR(B.OR(is_zero(B, args['x']), is_zero(B, args['y'])), is_zero(B, out)) ^ 1, f'MM0:{tag}')
        if 'K0' in L:
            for tag, _, args, out in kis:
                self.add(B.XOR(is_zero(B, args['x']), is_zero(B, out)) ^ 1, f'K0:{tag}')

        def kinv_of(word):
            for tag, _, args, out in kis:
                if word == out:
                    return args['x']
            return None
        if 'C1' in L:   # MM(u, MM(KINV(u), v)) = v ; пара: шаг 3 после шага 2
            tA, _, aA, oA = cmap['step3.MM0']
            tB, _, aB, oB = cmap['step2.MM1']
            u_k = kinv_of(aB['x'])
            ante = all_of(B, eq(B, aA['x'], u_k), is_zero(B, u_k) ^ 1, eq(B, aA['y'], oB))
            self.add(B.OR(ante ^ 1, eq(B, oA, aB['y'])), f'C1:{tA}({tB})')
        if 'C2' in L:   # MM(KINV(u), MM(u, w)) = w ; пара: шаг 5 после шага 4
            tA, _, aA, oA = cmap['step5.MM1']
            tB, _, aB, oB = cmap['step4.MM1']
            u_k = kinv_of(aA['x'])
            ante = all_of(B, eq(B, aB['x'], u_k), is_zero(B, u_k) ^ 1, eq(B, aA['y'], oB))
            self.add(B.OR(ante ^ 1, eq(B, oA, aB['y'])), f'C2:{tA}({tB})')
        if 'EC_DBL' in L:
            ante = all_of(B, e1, e2 ^ 1, e3 ^ 1, e4 ^ 1)
            self.add(B.OR(ante ^ 1, B.AND(eq(B, y, b), is_zero(B, b) ^ 1)), 'EC_DBL')
        if 'EC_TAN' in L:
            self.add(B.OR(B.AND(gb, eq(B, XR, a)) ^ 1, eq(B, Lam, lr)), 'EC_TAN')
        if 'EC_RNOTO' in L:
            self.add(B.OR(gb ^ 1, B.AND(is_zero(B, XR), is_zero(B, YR)) ^ 1), 'EC_RNOTO')
        if 'EC_RNOTP' in L:
            self.add(B.OR(gb ^ 1, B.AND(eq(B, XR, a), eq(B, YR, b)) ^ 1), 'EC_RNOTP')
        zero = bv.const(0, n)
        fx = ite(B, e3, x, ite(B, e4, a, ite(B, B.AND(e1, e2), zero, XR)))
        fy = ite(B, e3, y, ite(B, e4, b, ite(B, B.AND(e1, e2), zero, YR)))
        obligations.append(('out:x', bv.eq(B, r6['x'], fx) ^ 1))
        obligations.append(('out:y', bv.eq(B, r6['y'], fy) ^ 1))
        obligations.append(('out:a', bv.eq(B, r6['a'], a) ^ 1))
        obligations.append(('out:b', bv.eq(B, r6['b'], b) ^ 1))
        row = dict(case=label)
        if self.final_check == 'correct':
            for tag, d in obligations:
                if not self.holds(d):
                    A = ONE
                    for e, _ in self.conj:
                        A = B.AND(A, e)
                    asg = B.pick_sat(B.AND(A, d)) or {}
                    cex = {k: sum(1 << i for i, e in enumerate(self.cols[k]) if asg.get(B.var[e >> 1], 0))
                           for k in self.names}
                    row.update(status='REFUTED', failed=tag, counterexample={k: str(v) for k, v in cex.items()})
                    self.case_rows.append(row)
                    self.result.update(status='REFUTED', failed_case=label, failed=tag,
                                       counterexample=row['counterexample'])
                    self.pop_to(mark)
                    raise _Stop()
            row['status'] = 'PROVED'
        else:
            bad = ZERO
            for tag, d in obligations:
                if tag in ('pre:step5#5', 'pre:step6#4', 'pre:step6#5', 'pre:step6#6', 'pre:step6#7', 'pre:step6#8'):
                    bad = B.OR(bad, d)
            ok = self.holds(bad ^ 1)
            row['status'] = 'ALWAYS_FAILS' if ok else 'NOT_ALWAYS_FAILS'
            if not ok:
                self.result['status'] = 'NOT_ALWAYS_INCORRECT'
        row['seconds'] = round(time.time() - self.t0, 1)
        self.case_rows.append(row)
        self.log(f"[glue {row['seconds']:7.1f}s] {label}: {row['status']}", flush=True)
        self.pop_to(mark)


class _Split(Exception):
    def __init__(self, f):
        self.f = f


class _Stop(Exception):
    pass


def groups(B, W, F, p):
    """Группы находки 05 (README, раздел «Итог проверки 22.11») как предикаты входов."""
    a, b, x, y, lr = (W[k] for k in ('a', 'b', 'x', 'y', 'lam_r'))
    finite = B.AND(F['e3'] ^ 1, F['e4'] ^ 1)
    g1 = all_of(B, finite, F['e1'] ^ 1, F['e2'])
    g2 = all_of(B, F['e3'] ^ 1, F['e1'], eq(B, b, y), is_zero(B, b) ^ 1, is_zero(B, lr))
    g3 = all_of(B, F['e3'], is_zero(B, x) ^ 1, is_zero(B, y) ^ 1, is_zero(B, lr))
    return g1, g2, g3


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True)
    ap.add_argument('--mode', default='E2', choices=['E2', 'E3', 'E3fail', 'canary'])
    ap.add_argument('--orig', default='', help='шаги в исходной версии, напр. 1,2,6')
    ap.add_argument('--drop', default='', help='убрать леммы (канарейки склейки)')
    ap.add_argument('--group', type=int, default=0)
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    n, p = a.n, eval(a.p, {}, {})
    orig = {int(s) for s in a.orig.split(',') if s}
    versions = {k: k not in orig for k in range(1, 7)}
    lem = set(LEMMAS) - {s for s in a.drop.split(',') if s}
    extra, final = None, 'correct'
    if a.mode == 'E3':
        extra = lambda B, W, F: [g ^ 1 for g in groups(B, W, F, p)]
    elif a.mode == 'E3fail':
        extra = lambda B, W, F: [groups(B, W, F, p)[a.group - 1]]
        final = 'incorrect'
    G = Glue(n, p, versions=versions, lemmas=lem, extra=extra)
    try:
        res = G.run(final_check=final)
    except GlueFail as e:
        res = dict(status='ERROR', detail=str(e))
    res.update(mode=a.mode, orig=sorted(orig), dropped=sorted(set(LEMMAS) - lem), group=a.group)
    print(json.dumps(res, ensure_ascii=False)[:2000])
    if a.out:
        json.dump(res, open(a.out, 'w'), ensure_ascii=False, indent=1)


if __name__ == '__main__':
    main()
