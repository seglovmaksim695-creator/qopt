"""E4: ECWindowAddR (исправленная) = выборка QROAMClean (K=1) + ECAdd (ящик с контрактом E2) + снятие выборки
QROAMCleanAdjoint (правило C1 из shared/rules, без изменений). Доказатель ew (ROBDD qoptct) для ВСЕХ ctrl (w бит) и
x, y < p проверяет: на входе ящика ECAdd стоят ровно (T_a[c], T_b[c], x, y, T_lam[c]) — таблица посчитана ЗДЕСЬ
независимо (c·R обычной арифметикой точек, затем форма Монтгомери t·2^n mod p), каждая and/erase корректна,
после снятия выборки все провода, кроме ctrl, x, y, равны 0, ctrl не меняется, x, y = выходы ящика.
Вместе с E2 (ECAdd: (a,b,x,y,lam_r) -> (a,b,P+Q,lam_r), остальное чисто) это даёт ECWindowAddR: (c, Q) -> (c, Q + c·R).
Примеры:  python -m ew.e4 prove --n 16 --p 32749 --w 4 --out results/e4_n16_w4.json
          python -m ew.e4 canary --n 256 --p '2**255-19' --w 4 --count 24 --out results/e4_canary_n256.json
          python -m ew.e4 indep --n 256 --p '2**255-19' --w 4 --out results/e4_indep_n256.json"""
import argparse, copy, json, random, resource, sys, time
from qoptct.verify.bdd import ONE, ZERO
from . import bv
from .core import Run, Fail, lsb, build_program, MSB_ROOT

CURVE_A, CURVE_B = 2, 3            # y^2 = x^3 + 2x + 3 (B != 0: (0,0) не на кривой), выбор произвольный


# ------------------------------------------------------------------ независимая арифметика точек (None = O)
def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q, s = q // 2, s + 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % p, i + 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def on_curve(P, p):
    return P is None or (P[1] * P[1] - (P[0] ** 3 + CURVE_A * P[0] + CURVE_B)) % p == 0


def pt_add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + CURVE_A) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def pt_mul(c, P, p):
    R = None
    for _ in range(c):
        R = pt_add(R, P, p)
    return R


def base_point(p, w, start=5):
    assert (4 * CURVE_A ** 3 + 27 * CURVE_B ** 2) % p != 0
    x = start
    while True:
        y = sqrt_mod(x ** 3 + CURVE_A * x + CURVE_B, p)
        if y:
            R = (x, min(y, p - y))
            mults = [pt_mul(c, R, p) for c in range(1, 2 ** w)]
            if all(M is not None and M[1] != 0 for M in mults):   # c·R != O, нет 2-кручения в таблице
                return R
        x += 1


def table(n, p, w, R):
    mont = lambda t: t * (1 << n) % p
    Ta, Tb, Tl = [0], [0], [0]
    P = None
    for c in range(1, 2 ** w):
        P = pt_add(P, R, p)
        x, y = P
        assert on_curve(P, p) and y != 0
        Ta.append(mont(x)); Tb.append(mont(y)); Tl.append(mont((3 * x * x + CURVE_A) * pow(2 * y, -1, p) % p))
    return Ta, Tb, Tl


# ------------------------------------------------------------------ Qualtran-блок и разложение
def win_bloq(n, p, w, R):
    from qualtran.bloqs.cryptography.ecc.ec_add_r import ECWindowAddR
    from qualtran.bloqs.cryptography.ecc.ec_point import ECPoint
    return ECWindowAddR(n=n, R=ECPoint(R[0], R[1], mod=p, curve_a=CURVE_A), add_window_size=w, mul_window_size=1)


def make_rule():
    from qoptct.foreign.qualtran.bridge.qroam_clean_adjoint import chained_atom_rule
    from qoptct.foreign.qualtran.bridge.qualtran_layout import AtomExpansion
    from qualtran import Side
    from qualtran.bloqs.cryptography.ecc.ec_add import ECAdd
    from .opaque import state_atom
    crule = chained_atom_rule()

    def rule(bloq, inputs, allocate, path):
        st = state_atom(bloq, inputs, allocate, path)
        if st is not None:
            return st
        if type(bloq) is ECAdd:
            if any(r.side != Side.THRU for r in bloq.signature) or bloq.window_size != 1:
                raise ValueError('ECAdd: ожидались THRU-регистры и window_size=1')
            regs = {name: [int(q) for q in arr.flat] for name, arr in inputs.items()}
            qubits = sorted({q for w in regs.values() for q in w})
            leaf = dict(name='opaque', kind='ECADD', regs_in=regs, regs_out={k: list(v) for k, v in regs.items()},
                        qubits=qubits, bloq_repr=repr(bloq), path=path)
            return AtomExpansion(dict(inputs), [leaf], ())
        return crule(bloq, inputs, allocate, path)
    return rule


def lower_e4(bloq):
    from qoptct.foreign.qualtran.bridge.qualtran_layout import lower_bloqgraph
    from .translate_ew import translate_ew
    q = bloq.qrom
    if any(int(k) != 0 for k in q.log_block_sizes):
        raise Fail('QROAMClean с K>1 (SwapWithZero мостом не переводится)', log_block_sizes=str(q.log_block_sizes))
    low = lower_bloqgraph(bloq, extra_atom_rule=make_rule(), max_leaves=5_000_000, max_nodes=20_000_000)
    tr = translate_ew(low)
    for key in ('wire_names', 'graph_events', 'source_wire_names'):
        low.pop(key, None)
        tr['layout'].pop(key, None)
    return low, tr


# ------------------------------------------------------------------ спецификация и прогон
class E4Spec:
    name = 'e4'
    range_words = ('x', 'y')

    def __init__(self, n, p, T, w):
        self.n, self.p, self.T, self.w = n, p, T, w
        self.cnames = [f'ctrl[{i}]' for i in range(w)]       # мост раскладывает массив ctrl на ctrl[0..w-1]

    def sel(self, W):
        # bb.join: ctrl[0] — старший бит; слово LSB-first: [ctrl[w-1], ..., ctrl[0]]
        return [W[f'ctrl[{self.w - 1 - j}]'][0] for j in range(self.w)]

    def lookup(self, B, sel, data):
        word = bv.const(data[0], self.n)
        for c in range(1, len(data)):
            word = bv.ite(B, bv.eq(B, sel, bv.const(c, len(sel))), bv.const(data[c], self.n), word)
        return word

    def pre(self, B, W, O):
        return [bv.lt_const(B, W['x'], self.p), bv.lt_const(B, W['y'], self.p)]

    def calls(self, B, W, O):
        Ta, Tb, Tl = self.T
        s = self.sel(W)
        return [('ECADD', dict(a=self.lookup(B, s, Ta), b=self.lookup(B, s, Tb),
                               lam_r=self.lookup(B, s, Tl), x=W['x'], y=W['y']))]

    def post(self, B, W, O):
        return dict({k: W[k] for k in self.cnames}, x=O[0]['x'], y=O[0]['y'])


class RunE4(Run):
    def _alloc(self):
        B, n = self.B, self.n
        ins = self.prog['inputs']
        fwd = [it[1] for it in self.prog['items'] if it[0] == 'opaque']
        if [l['kind'] for l in fwd] != ['ECADD']:
            raise Fail('ожидался ровно один ящик ECADD', got=[l['kind'] for l in fwd])
        cn = self.spec.cnames
        if set(ins) != set(cn) | {'x', 'y'} or len(ins['x']) != n or len(ins['y']) != n or any(len(ins[k]) != 1 for k in cn):
            raise Fail('сигнатура ECWindowAddR', got={k: len(v) for k, v in ins.items()})
        self.fwd_leaves = fwd
        ctrl = {k: [B.new_var()[1]] for k in cn}             # самый низ (при прогоне — разбор по значениям ctrl)
        cols = {k: [None] * n for k in ('x', 'y', 'X', 'Y')}
        for i in (range(n) if MSB_ROOT else reversed(range(n))):
            for k in ('x', 'y', 'X', 'Y'):
                cols[k][i] = B.new_var()[1]
        self.W = dict(ctrl, x=cols['x'], y=cols['y'])
        self.O = [{'x': cols['X'], 'y': cols['Y']}]
        self.N = {}
        self.spec.N = self.N
        self.nvars = B.nvars

    def _opaque(self, leaf, fwd_i):
        B, p, val = self.B, self.p, self.val
        if leaf['kind'] != 'ECADD' or fwd_i != 0 or self.fwd_leaves[0] is not leaf:
            raise Fail('неожиданный ящик', kind=leaf['kind'])
        rin = {r: lsb(ws) for r, ws in leaf['regs_in'].items()}
        rout = {r: lsb(ws) for r, ws in leaf['regs_out'].items()}
        regs = ('a', 'b', 'x', 'y', 'lam_r')
        if set(rin) != set(regs) or rin != rout:
            raise Fail('ECADD: сигнатура (все регистры THRU)')
        cur = {r: [val[q] for q in ws] for r, ws in rin.items()}
        ekind, exp = self.expected[0]
        for r in regs:
            if not self.word_eq_holds(cur[r], exp[r]):
                self.fail(f'ящик ECADD: вход {r} не равен заявленному', bv.eq(B, cur[r], exp[r]) ^ 1)
            lt = bv.lt_const(B, cur[r], p)
            if not self.holds(lt ^ 1):
                self.fail(f'ящик ECADD: вход {r} вне области <p', lt ^ 1)
        for r in ('x', 'y'):                       # контракт E2: a, b, lam_r не меняются; x, y — результат ящика
            for q, e in zip(rout[r], self.O[0][r]):
                val[q] = e
        return 1


def setup(n, p, w):
    R = base_point(p, w)
    T = table(n, p, w, R)
    bloq = win_bloq(n, p, w, R)
    return R, T, bloq


def version_info():
    import hashlib, os
    from qoptct.foreign.qualtran._variant import detect_variant
    from .run_step import code_hashes
    v = detect_variant()
    hs = code_hashes()
    here = os.path.dirname(os.path.abspath(__file__))
    rule = os.path.join(os.path.dirname(here), 'bridge', 'qroam_clean_adjoint.py')
    hs['ecadd/e4.py'] = hashlib.sha256(open(os.path.join(here, 'e4.py'), 'rb').read()).hexdigest()[:16]
    hs['bridge/qroam_clean_adjoint.py'] = hashlib.sha256(open(rule, 'rb').read()).hexdigest()[:16]
    return dict(qualtran_file=v['qualtran_file'], kind=v['kind'], files_sha256_16=v['files_sha256_16'],
                version='patched(e4)' if v['kind'] == 'e4' else v['kind'], code_sha256_16=hs)


def prove(n, p, w, lowered=None, mutate=None, log=print, setup_=None):
    R, T, bloq = setup_ or setup(n, p, w)
    low, tr = lowered or lower_e4(bloq)
    prog = build_program(None, lowered=(low, tr), mutate=mutate, keep_ops=False)
    run = RunE4(prog, E4Spec(n, p, T, w), n, p, log=log)
    return run.run(), run, prog


# ------------------------------------------------------------------ конкретное исполнение (свидетель для канареек)
def exec_concrete(ops, leaves, N, low, T, n, p, w, ctrl, x, y, X, Y):
    v = bytearray(N)
    put = lambda ws, val: [v.__setitem__(q, (val >> j) & 1) for j, q in enumerate(ws)]
    get = lambda ws: sum(v[q] << j for j, q in enumerate(ws))
    for i in range(w):
        v[low['inputs'][f'ctrl[{i}]'][0]] = (ctrl >> (w - 1 - i)) & 1
    put(low['inputs']['x'], x); put(low['inputs']['y'], y)
    for op in ops:
        k = op[0]
        if k == 'x':
            v[op[1]] ^= 1
        elif k == 'cx':
            v[op[2]] ^= v[op[1]]
        elif k == 'and':
            if v[op[3]]:
                return 'and'
            v[op[3]] = v[op[1]] & v[op[2]]
        elif k == 'erase':
            (cz, a, b), = op[2]
            if v[op[1]] != (v[a] & v[b]):
                return 'erase'
            v[op[1]] = 0
        else:
            leaf = leaves[op[1]]
            rin = {r: lsb(ws) for r, ws in leaf['regs_in'].items()}
            want = dict(a=T[0][ctrl], b=T[1][ctrl], lam_r=T[2][ctrl], x=x, y=y)
            if any(get(rin[r]) != want[r] for r in want):
                return 'box-input'
            put(rin['x'], X); put(rin['y'], Y)
    outs = low['outputs']
    if any(v[outs[f'ctrl[{i}]'][0]] != (ctrl >> (w - 1 - i)) & 1 for i in range(w)) or get(outs['x']) != X or get(outs['y']) != Y:
        return 'output'
    ow = {q for ws in outs.values() for q in ws}
    if any(v[q] for q in range(N) if q not in ow):
        return 'nonzero'
    return None


def witness(ops, leaves, N, low, T, n, p, w, rnd, trials=40):
    for t in range(trials):
        for c in range(2 ** w):
            x, y, X, Y = (rnd.randrange(p) for _ in range(4))
            if t == 0:
                x = y = 0
            r = exec_concrete(ops, leaves, N, low, T, n, p, w, c, x, y, X, Y)
            if r:
                return r, c
    return None, None


def mutations(tr, low, rnd, count):
    from .regions import regions
    ops, reg = tr['ops'], regions(tr)
    live_in = [q for ws in low['inputs'].values() for q in ws]
    muts = []
    add = lambda name, targeted, fn: muts.append((name, targeted, fn))
    li = next(i for i, l in enumerate(low['leaves']) if l['name'] == 'opaque')
    for r in ('a', 'b', 'lam_r', 'x', 'y'):                       # стыки с ящиком
        def fn(tr, low, r=r):
            L = list(low['leaves']); L[li] = copy.deepcopy(L[li]); ws = L[li]['regs_in'][r]
            ws[-1], ws[-2] = ws[-2], ws[-1]; L[li]['regs_out'][r] = ws
            return tr['ops'], L
        add(f'joint:ECADD:swap-low-bits-{r}', True, fn)
    keys = list(reg)
    adj_idx = [i for k in keys if 'Adjoint' in k for i in range(reg[k][0], reg[k][1])]
    fwd_idx = [i for k in keys if 'Adjoint' not in k and 'QROAM' in k for i in range(reg[k][0], reg[k][1])]
    er = [i for i in adj_idx if ops[i][0] == 'erase']
    for i in rnd.sample(er, min(7, len(er))):                       # поправки стираний снятия выборки
        def fn(tr, low, i=i):
            o = list(tr['ops']); t, ((cz, a, b),) = o[i][1], o[i][2]
            c = next(q for q in live_in if q not in (a, b, t))
            o[i] = ('erase', t, (('cz', a, c),)); return o, low['leaves']
        add(f'erase-correction:op{i}', True, fn)
    ands = [i for i in fwd_idx if ops[i][0] == 'and']
    for i in rnd.sample(ands, min(4, len(ands))):                   # редкие входы: одна ветка лестницы выборки
        def fn(tr, low, i=i):
            o = tr['ops']; _, a, b, t = o[i]
            j = next(j for j in range(i + 1, len(o)) if o[j][0] == 'erase' and o[j][1] == t)
            return o[:i] + [('x', a), o[i], ('x', a)] + o[i + 1:j] + [('x', a), o[j], ('x', a)] + o[j + 1:], low['leaves']
        add(f'rare-input:flip-ladder-control-op{i}', True, fn)
    cxs = [i for i in fwd_idx if ops[i][0] == 'cx']
    for i in rnd.sample(cxs, min(3, len(cxs))):                     # один бит таблицы
        add(f'table-bit:drop-cx-op{i}', True, lambda tr, low, i=i: (tr['ops'][:i] + tr['ops'][i + 1:], low['leaves']))
    while len(muts) < count:
        i = rnd.randrange(len(ops))
        if ops[i][0] == 'opaque':
            continue
        kind = rnd.choice(['drop', 'xins', 'retarget'])
        if kind == 'drop':
            add(f'random:drop-op{i}', False, lambda tr, low, i=i: (tr['ops'][:i] + tr['ops'][i + 1:], low['leaves']))
        elif kind == 'xins':
            q = rnd.choice(live_in)
            add(f'random:x-on-{q}-after-op{i}', False,
                lambda tr, low, i=i, q=q: (tr['ops'][:i + 1] + [('x', q)] + tr['ops'][i + 1:], low['leaves']))
        elif ops[i][0] == 'cx':
            q = rnd.choice([u for u in live_in if u not in ops[i][1:]])
            add(f'random:cx-control-op{i}->{q}', False,
                lambda tr, low, i=i, q=q: (tr['ops'][:i] + [('cx', q, tr['ops'][i][2])] + tr['ops'][i + 1:], low['leaves']))
    return muts[:max(count, len(muts))]


# ------------------------------------------------------------------ независимая проверка средствами Qualtran
def indep(n, p, w, seed=11, rand=12):
    import numpy as np
    from qualtran import QMontgomeryUInt
    R, T, bloq = setup(n, p, w)
    rows = []
    q = bloq.qrom
    data = [list(map(int, np.asarray(d).ravel())) for d in q.data]   # таблица Qualtran
    rows.append(dict(check='table == независимая (c·R, Монтгомери)', ok=data == [list(t) for t in T]))
    dt = QMontgomeryUInt(n, p)
    mont = lambda t: t * (1 << n) % p
    rows.append(dict(check='uint_to_montgomery == t·2^n mod p', ok=all(int(dt.uint_to_montgomery(t)) == mont(t)
                                                                      for t in (0, 1, 2, p - 1, R[0], R[1]))))
    rnd = random.Random(seed)
    pts = [None, R, pt_mul(2, R, p)]
    for c in range(1, 2 ** w):
        P = pt_mul(c, R, p)
        pts += [P, (P[0], (-P[1]) % p)]               # Q = c·R (удвоение) и Q = -c·R (результат O)
    while len(pts) < 2 ** w * 2 + 3 + rand:
        xx = rnd.randrange(p)
        yy = sqrt_mod(xx ** 3 + CURVE_A * xx + CURVE_B, p)
        if yy is not None:
            pts.append((xx, yy))
    cb = bloq.decompose_bloq()
    bad, n_ok, errs = [], 0, []
    for i, Q in enumerate(pts):
        c = (i * 7 + 3) % (2 ** w) if i >= 3 else i % (2 ** w)
        if Q is not None and i >= 3 and i < 3 + 2 * (2 ** w - 1):
            c = (i - 3) // 2 + 1                                   # точно c, дающее удвоение / O
        want = pt_add(Q, pt_mul(c, R, p), p)
        xm, ym = (0, 0) if Q is None else (mont(Q[0]), mont(Q[1]))
        wx, wy = (0, 0) if want is None else (mont(want[0]), mont(want[1]))
        bits = np.array([(c >> (w - 1 - j)) & 1 for j in range(w)])
        try:
            out = cb.call_classically(ctrl=bits, x=xm, y=ym)
            got = (int(out[1]), int(out[2]))
            ok = got == (wx, wy) and list(map(int, out[0])) == list(bits)
        except Exception as e:
            ok, got = None, None
            errs.append(f'{type(e).__name__}: {str(e)[:120]}')
        if ok:
            n_ok += 1
        elif ok is False:
            bad.append(dict(c=c, Q=str(Q), got=str(got)))
    rows.append(dict(check='разложение ECWindowAddR (Qualtran call_classically) == Q + c·R', total=len(pts), ok_count=n_ok,
                     bad=bad[:5], errors=sorted(set(errs))[:3]))
    return dict(R=[str(R[0]), str(R[1])], curve=dict(a=CURVE_A, b=CURVE_B), rows=rows)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prove', 'canary', 'indep'])
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--p', required=True)
    ap.add_argument('--w', type=int, default=4)
    ap.add_argument('--count', type=int, default=24)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    n, p, w = a.n, eval(a.p, {}, {}), a.w
    t0 = time.time()
    res = dict(cmd=a.cmd, n=n, p=str(p), w=w, mul_window_size=1, curve=f'y^2=x^3+{CURVE_A}x+{CURVE_B}', **version_info(),
               assumptions=['E2: ECAdd(n,p,1) исправленная = сложение точек, a,b,lam_r не меняются, остальное чисто',
                            'C1 правило QROAMCleanAdjoint (A1-A4 цели C)', 'QROAMClean K=1 раскрыт мостом (как у C)'])
    try:
        if a.cmd == 'prove':
            from .regions import regions
            S = setup(n, p, w)
            low, tr = lower_e4(S[2])
            res['regions'] = {k: v[2] for k, v in regions(tr).items()}
            res['R'] = [str(S[0][0]), str(S[0][1])]
            st, run, prog = prove(n, p, w, lowered=(low, tr), setup_=S)
            res.update(status=st, stats=run.stats, logical_ops=prog['nops'])
        elif a.cmd == 'indep':
            res.update(indep(n, p, w))
            r3 = res['rows'][2]
            res['status'] = 'OK' if all(r.get('ok', True) for r in res['rows'][:2]) and r3['ok_count'] == r3['total'] else 'CHECK'
        else:
            rnd = random.Random(a.seed)
            S = setup(n, p, w)
            low, tr = lower_e4(S[2])
            rows = []
            for name, targeted, fn in mutations(tr, low, rnd, a.count):
                t1 = time.time()
                try:
                    st, run, prog = prove(n, p, w, lowered=(low, tr), mutate=fn, log=lambda *x, **k: None, setup_=S)
                    why = ''
                except Fail as e:
                    st, why = 'REFUTED', str(e)[:120]
                except Exception as e:
                    st, why = 'ERROR', f'{type(e).__name__}: {str(e)[:120]}'
                ops, leaves = fn(tr, low)
                wk, wc = witness(ops, leaves, tr['nq'], low, S[1], n, p, w, rnd)
                row = dict(name=name, targeted=targeted, status=st, why=why, witness=wk is not None, witness_kind=wk,
                           witness_ctrl=wc, seconds=round(time.time() - t1, 1))
                if st == 'PROVED' and wk is not None:
                    row['HOLE'] = True
                rows.append(row)
                print(json.dumps(row, ensure_ascii=False), flush=True)
            res.update(rows=rows, summary=dict(total=len(rows), refuted=sum(r['status'] == 'REFUTED' for r in rows),
                       refuted_with_witness=sum(r['status'] == 'REFUTED' and r['witness'] for r in rows),
                       targeted=sum(r['targeted'] for r in rows),
                       proved_neutral=sum(r['status'] == 'PROVED' and not r['witness'] for r in rows),
                       holes=sum(bool(r.get('HOLE')) for r in rows), errors=sum(r['status'] == 'ERROR' for r in rows)))
            res['status'] = 'OK' if res['summary']['holes'] == 0 and res['summary']['errors'] == 0 else 'CHECK'
    except Fail as e:
        res.update(status='REFUTED', detail=str(e), info={k: str(v)[:300] for k, v in e.detail.items()})
    except Exception as e:
        import traceback
        res.update(status='ERROR', detail=f'{type(e).__name__}: {e}', tb=traceback.format_exc()[-2500:])
    res['wall_seconds'] = round(time.time() - t0, 2)
    res['peak_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('rows',)}, ensure_ascii=False)[:2500])
    if a.out:
        json.dump(res, open(a.out, 'w'), ensure_ascii=False, indent=1)


if __name__ == '__main__':
    main()
