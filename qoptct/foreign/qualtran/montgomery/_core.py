"""Core of the Montgomery multiplication proof (ported verbatim from the round-3 scripts c2_lib / c2_prove).

The multiplier ``DirtyOutOfPlaceMontgomeryModMul`` is lowered gate by gate with the strict bridge and cut at the
boundaries of the Qualtran nodes: ``X(reduced)`` + window 0, windows 1..N-1, ``LessThanConstant`` and the
controlled ``AddK``.  Every segment is proved with the unchanged qoptct prover ``kprove`` (BDD path sum) plus the exact
partial-domain guard, and cross-checked by the independent simulator.  Wire glue, live-wire bookkeeping, the
interval invariants and the telescoping identity are checked by code.
"""
import os
import sys
import time
import random
from hashlib import sha256

from ..bridge.qroam_clean_adjoint import chained_atom_rule
from ..bridge.qualtran_layout import lower_bloqgraph
from ..bridge.rules import translate
from ..bridge import sim
from ..bridge.emit import emit_qasm
from ..bridge.strict_proof import prove_program, digest_json
from qoptct.verify.ksym import SymInt
from qoptct.verify.bdd import ONE, ZERO


def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


def prime_for(n):
    """n=256: 2^255-19; иначе наибольшее простое p <= 2^(n-1) (тогда 2p-1 < 2^n)"""
    if n == 256: return 2 ** 255 - 19
    p = 2 ** (n - 1)
    while not is_prime(p): p -= 1
    return p


def build(n, w, p):
    from qualtran.bloqs.mod_arithmetic.mod_multiplication import DirtyOutOfPlaceMontgomeryModMul
    return DirtyOutOfPlaceMontgomeryModMul(bitsize=n, window_size=w, mod=p)


def lower_translate(bloq):
    low = lower_bloqgraph(bloq, extra_atom_rule=chained_atom_rule(), max_leaves=5_000_000, max_nodes=5_000_000)
    tr = translate(low, low['inputs'], low['outputs'], light=True)
    return low, tr


def reference(n, w, p, x, y):
    """точная рекуррентность Монтгомери (правильная версия: вычитание p). -> dict outputs"""
    N = 1 << w
    inv = pow(p, 2 ** (w - 1) - 1, N)
    T = 0; hist = 0
    for i in range(0, n, w):
        xd = (x >> i) & (N - 1)
        T += xd * y
        m = T & (N - 1)
        T = (T + ((-m * inv) % N) * p) >> w
        hist = (hist << w) | m
    red = int(T >= p)
    return dict(x=x, y=y, target=T - p * red, qrom_indices=hist, reduced=red, T=T)


def reference_orig(n, w, p, x, y):
    """неисправленная версия: при T>=p прибавляется +p (по модулю 2^n)"""
    r = reference(n, w, p, x, y)
    t = r['T'] + p if r['reduced'] else r['T']
    r['target'] = t % (1 << n)
    return r


def run_basis(low, tr, n, x, y):
    ins = low['inputs']; outs = low['outputs']
    init = 0
    for k in range(n):
        init |= ((x >> k) & 1) << ins['x'][k]
        init |= ((y >> k) & 1) << ins['y'][k]
    val, ph = sim.run(tr['nq'], tr['ops'], init)
    res = {}
    for name, wires in outs.items():
        res[name] = sum(((val >> w_) & 1) << k for k, w_ in enumerate(wires))
    rest = val
    for wires in outs.values():
        for w_ in wires: rest &= ~(1 << w_)
    return res, ph, rest


def compact(ops, words_in, words_out):
    order = []; seen = set()

    def add(q):
        if q not in seen:
            seen.add(q); order.append(q)
    for d in (words_in, words_out):          # КАНОНИЧЕСКАЯ нумерация: сначала слова контракта, потом по первому появлению в операциях
        for ws in d.values():
            for q in ws: add(q)
    for o in ops:
        if o[0] == 'x': add(o[1])
        elif o[0] == 'cx': add(o[1]); add(o[2])
        elif o[0] == 'and': add(o[1]); add(o[2]); add(o[3])
        elif o[0] == 'erase':
            add(o[1])
            for c in o[2]:
                for q in c[1:]: add(q)
    mp = {q: i for i, q in enumerate(order)}

    def mo(o):
        if o[0] == 'x': return ('x', mp[o[1]])
        if o[0] == 'cx': return ('cx', mp[o[1]], mp[o[2]])
        if o[0] == 'and': return ('and', mp[o[1]], mp[o[2]], mp[o[3]])
        return ('erase', mp[o[1]], tuple((c[0],) + tuple(mp[q] for q in c[1:]) for c in o[2]))
    cops = [mo(o) for o in ops]
    wi = {k: [mp[q] for q in v] for k, v in words_in.items()}
    wo = {k: [mp[q] for q in v] for k, v in words_out.items()}
    return len(order), cops, wi, wo, mp


def make_layout(nq, wi, wo, ops):
    free = {q for v in wi.values() for q in v}
    return dict(nq=nq, inputs=wi, outputs=wo, initial_zero=[q for q in range(nq) if q not in free], garbage=[], free_wires=[],
                partial_gate_domains_required=any(o[0] in ('and', 'erase') for o in ops), domain_policy='exact_logical_AND_ERASE_v1')


class Params:
    def __init__(s, n, w, p):
        s.n, s.w, s.p = n, w, p
        s.N = 1 << w
        s.inv = pow(p, 2 ** (w - 1) - 1, s.N)
        s.c = [(-m * s.inv) % s.N for m in range(s.N)]
        for m in range(s.N):   # целочисленный факт: p*c(m)+m = 0 mod 2^w
            assert (p * s.c[m] + m) % s.N == 0
        s.table = [s.c[m] * p for m in range(s.N)]


def win_eval(P, xd, y, T):
    S = T + xd * y
    m = S & (P.N - 1)
    Z = S + P.table[m]
    assert Z & (P.N - 1) == 0
    return m, Z >> P.w


def make_win(P, first, wi, wo):
    n, w, p = P.n, P.w, P.p

    def pre(S):
        y = S.reg(wi['y'])
        c = y.lt(p)
        if not first:
            T = S.reg(wi['T'])
            c = S.B.AND(c, T.lt(y + p))
        return c

    def spec(S):
        xs = S.reg(wi['x']); y = S.reg(wi['y'])
        T = S.reg(wi['T']) if not first else S.const(0, n + w)
        Ssum = T
        for k in range(w):
            term = SymInt(S.B, [S.B.AND(xs.bits[k], b) for b in y.bits]).shl(k)
            Ssum = Ssum + term
        m = SymInt(S.B, Ssum.bits[:w])
        cp = S.lookup(m, P.table, n + w + 1)
        Z = Ssum + cp
        Tn = Z.shr(w)
        out = {}
        out.update(S.out(wo['x'], xs)); out.update(S.out(wo['y'], y)); out.update(S.out(wo['T'], Tn)); out.update(S.out(wo['m'], m))
        if first: out[wo['r'][0]] = ONE
        return out
    return pre, spec


def make_lt(P, wi, wo):
    p = P.p

    def spec(S):
        t = S.reg(wi['t']); r = S.reg(wi['r'])
        out = S.out(wo['t'], t)
        out[wo['r'][0]] = S.B.XOR(r.bits[0], t.lt(p))
        return out
    return None, spec


def make_add(P, wi, wo):
    n, p = P.n, P.p

    def spec(S):
        t = S.reg(wi['t']); r = S.reg(wi['r'])
        add = (t + S.const((1 << n) - p, n)).resize(n)
        res = S.mux(r.bits[0], add, t)
        out = S.out(wo['t'], res)
        out[wo['r'][0]] = r.bits[0]
        return out
    return None, spec


def bits_set(val, wires):
    return sum(((val >> k) & 1) << q for k, q in enumerate(wires))


def sim_check(kind, P, nq, ops, wi, wo, first, rng, count):
    """независимая симуляция bridge.sim; ожидаемое считает отдельный целочисленный код"""
    n, w, p = P.n, P.w, P.p
    cases = []
    if kind == 'win':
        ys = sorted({0, 1, p - 1, p // 2} | {rng.randrange(p) for _ in range(count)})
        for y in ys:
            Ts = [0] if first else sorted({0, 1, y + p - 1, p - 1, p, max(0, y + p - 2), rng.randrange(y + p)})
            for T in Ts:
                for xd in sorted({0, P.N - 1, rng.randrange(P.N)}):
                    cases.append((xd, y, T))
    else:
        for _ in range(count):
            t = rng.choice([0, 1, p - 1, p, p + 1, 2 * p - 1, (1 << n) - 1, rng.randrange(1 << n), rng.randrange(2 * p)])
            cases.append((t, rng.randrange(2)))
    bad = 0
    first_bad = None
    for cs in cases:
        init = 0
        if kind == 'win':
            xd, y, T = cs
            init |= bits_set(xd, wi['x']) | bits_set(y, wi['y'])
            if not first: init |= bits_set(T, wi['T'])
            m, Tn = win_eval(P, xd, y, T)
            exp = bits_set(xd, wo['x']) | bits_set(y, wo['y']) | bits_set(Tn, wo['T']) | bits_set(m, wo['m'])
            if first: exp |= 1 << wo['r'][0]
        else:
            t, r = cs
            init = bits_set(t, wi['t']) | bits_set(r, wi['r'])
            if kind == 'lt':
                exp = bits_set(t, wo['t']) | ((r ^ int(t < p)) << wo['r'][0])
            else:
                exp = bits_set((t + r * ((1 << n) - p)) % (1 << n), wo['t']) | (r << wo['r'][0])
        val, ph = sim.run(nq, ops, init)
        if val != exp or any(ph):
            bad += 1
            if first_bad is None: first_bad = cs
    return len(cases), bad, first_bad


def run_segment(label, kind, P, ops_g, wi_g, wo_g, first, rng, budget, max_mb, tag='', cache=None):
    nq, cops, wi, wo, mp = compact(ops_g, wi_g, wo_g)
    layout = make_layout(nq, wi, wo, cops)
    qasm = emit_qasm(nq, cops)
    if kind == 'win':
        pre, spec = make_win(P, first, wi, wo)
    elif kind == 'lt':
        pre, spec = make_lt(P, wi, wo)
    else:
        pre, spec = make_add(P, wi, wo)
    contract = dict(source='C2 segment', kind=kind, n=P.n, w=P.w, p=P.p, first=first, qualtran=tag,
                    formula={'win': "S=T+xd*y; m=S mod 2^w; T'=(S+c(m)*p)>>w; outs x,y,T',m; pre y<p and T<y+p (first: T=0)",
                             'lt': "r'=r xor [t<p], t unchanged", 'add': "t'=(t+r*(2^n-p)) mod 2^n, r unchanged"}[kind])
    cnt, bad, first_bad = sim_check(kind, P, nq, cops, wi, wo, first, rng, 12)
    key = (sha256(qasm.encode()).hexdigest(), digest_json(layout), digest_json(contract))
    if cache is not None and key in cache and not bad:
        r = dict(cache[key]); r.update(label=label, dedup_of=cache[key]['label'], sim_cases=cnt, sim_bad=bad)
        return r
    t0 = time.time()
    receipt = dict(status='PASS', qasm_sha256=sha256(qasm.encode()).hexdigest(), contract_sha256=digest_json(contract),
                   layout_sha256=digest_json(layout), independent_inputs=cnt, all_wires_checked=True, phase_checked=True)
    if bad:
        return dict(label=label, kind=kind, status='SIM_MISMATCH', sim_cases=cnt, sim_bad=bad, nq=nq, ops=len(cops),
                    qasm_sha256=key[0], contract_sha256=key[2], witness=first_bad)
    v = prove_program(qasm, layout, spec, pre=pre, validation=receipt, contract=contract, budget=budget, max_mb=max_mb)
    r = dict(label=label, kind=kind, status=v['status'], detail=(v.get('detail') or '')[:160], nq=nq, ops=len(cops), sim_cases=cnt, sim_bad=bad,
                wall=round(time.time() - t0, 2), peak_mb=v.get('peak_rss_mb'), qasm_sha256=receipt['qasm_sha256'], contract_sha256=receipt['contract_sha256'],
                domain=(v.get('partial_domains') or {}).get('status'))
    if cache is not None: cache[key] = r
    return r


def segmentation(tr, low, n, w):
    top = {}; order = []
    for t in tr['trace']:
        key = t['leaf']['path'].split('/')[1]
        if key not in top:
            top[key] = [t['logical_start'], t['logical_end']]; order.append(key)
        else:
            top[key][1] = t['logical_end']
    wins = [k for k in order if 'SingleWindowModMul' in k]
    assert len(wins) == n // w, (len(wins), n // w)
    lt = [k for k in order if 'LessThanConstant' in k]; ad = [k for k in order if 'Controlled' in k]
    xg = [k for k in order if 'XGate' in k]
    assert len(lt) == 1 and len(ad) == 1 and len(xg) == 1 and len(order) == len(wins) + 3 and order[0] == xg[0], order
    pos = 0
    for k in order:
        assert top[k][0] == pos, (k, top[k], pos)
        pos = top[k][1]
    assert pos == len(tr['ops'])
    return top, order, wins, lt[0], ad[0], xg[0]


def support(ops, wi, wo):
    sup = set()
    for o in ops:
        if o[0] == 'x': sup.add(o[1])
        elif o[0] == 'cx': sup.update(o[1:3])
        elif o[0] == 'and': sup.update(o[1:4])
        else:
            sup.add(o[1])
            for c in o[2]: sup.update(c[1:])
    for d in (wi, wo):
        for ws in d.values(): sup.update(ws)
    return sup


def bookkeeping(segs, low, n, w):
    """G5: учёт живых проводов. Вход сегмента обязан быть живым; не-входные провода сегмента обязаны быть свежими (никто раньше не трогал)."""
    res = dict(ok=True, notes=[])
    prog_in = {q for k in ('x', 'y') for q in low['inputs'][k]}
    live = set(prog_in); touched = set(prog_in) & set()
    seen_prog_in = set(prog_in)
    for label, kind, sops, wi, wo, first in segs:
        ins = {q for v in wi.values() for q in v}; outs = {q for v in wo.values() for q in v}
        sup = support(sops, wi, wo)
        if not ins <= live: res['ok'] = False; res['notes'].append(f'{label}: вход не живой {sorted(ins - live)[:5]}')
        bad_fresh = (sup - ins) & (touched | seen_prog_in)
        if bad_fresh: res['ok'] = False; res['notes'].append(f'{label}: не-входные провода не свежие {sorted(bad_fresh)[:5]}')
        live = (live - ins) | outs
        touched |= sup
    outw = {q for k in ('x', 'y', 'target', 'qrom_indices', 'reduced') for q in low['outputs'][k]}
    lastT = segs[len(segs) - 3][4]['T']
    extra = live - outw
    if extra != set(lastT[n:]): res['ok'] = False; res['notes'].append(f'живые провода в конце != выходы + старшие w провода T: {sorted(extra)[:6]}')
    if not outw <= live: res['ok'] = False; res['notes'].append('не все выходные провода живы в конце')
    res['final_high_T_wires_must_be_zero'] = sorted(extra)
    return res


def invariant_checks(P, cases=(True, False)):
    """BDD: (1) pre => T'<y+p (инвариант окна), T'<2^n; (2) первое окно (T=0) => то же; (3) Z mod 2^w = 0 (целочисленный факт таблицы)"""
    from qoptct.verify.bdd import BDD, ZERO, ONE
    from qoptct.verify.ksym import SymInt, lookup
    n, w, p = P.n, P.w, P.p
    out = {}
    for first in cases:
        if ((not first) and n > 192) or w >= 4:
            # BDD-проверка generic-окна при n=256 не помещается в 3 ГБ; даём точное целочисленное доказательство (интервальная арифметика):
            # T<y+p, y<p, x<2^w, c(m)<2^w  =>  Z=T+x*y+c*p <= 2^w*(y+p)-1  =>  (Z делится на 2^w)  T'=Z/2^w < y+p < 2p <= 2^n.
            N_ = 1 << w
            tab_ok = all(0 <= P.c[m] < N_ for m in range(N_))
            low_ok = all((m + P.c[m] * p) % N_ == 0 for m in range(N_))   # (T+x*y) mod 2^w = m  =>  Z mod 2^w = 0
            bound_ok = (2 * p <= (1 << n))
            out['first' if first else 'generic'] = dict(inv_T_lt_y_plus_p=bool(tab_ok and low_ok), Tn_lt_2n=bool(bound_ok), low_bits_zero=bool(low_ok), method='analytic_interval_arithmetic_not_BDD')
            continue
        B = BDD()
        mk = lambda k: SymInt(B, [B.new_var()[1] for _ in range(k)])
        xs = mk(w)
        # порядок как у спецификации: биты слов чередуются; здесь свои переменные (порядок для проверки импликации)
        y = SymInt(B, []); T = SymInt(B, [])
        yb, Tb = [], []
        for i in range(n + w):
            if i < n: yb.append(B.new_var()[1])
            Tb.append(B.new_var()[1])
        y = SymInt(B, yb); T = SymInt(B, Tb) if not first else SymInt.const(B, 0, n + w)
        pre = y.lt(p)
        if not first: pre = B.AND(pre, T.lt(y + p))
        Ssum = T
        for k in range(w):
            Ssum = Ssum + SymInt(B, [B.AND(xs.bits[k], b) for b in y.bits]).shl(k)
        _syms = [xs, y, T, Ssum]
        _rm = B.gc([pre] + [e for s_ in _syms for e in s_.bits]); pre = _rm(pre)
        for s_ in _syms: s_.bits = [_rm(e) for e in s_.bits]
        if os.environ.get('INV_DEBUG'): print('STAGE', 'ssum-gc', B.nodes(), __import__('resource').getrusage(0).ru_maxrss//10**6, flush=True)
        m = SymInt(B, Ssum.bits[:w])
        Z = Ssum + lookup(B, m, P.table, n + w + 1)
        Tn = Z.shr(w)
        lowz = ZERO
        for b in Z.bits[:w]: lowz = B.OR(lowz, b) if hasattr(B, 'OR') else (lowz ^ b ^ B.AND(lowz, b))
        B.c_and.clear(); B.c_xor.clear()
        if os.environ.get('INV_DEBUG'): print('STAGE', 'tn', B.nodes(), __import__('resource').getrusage(0).ru_maxrss//10**6, flush=True)
        r1 = B.AND(pre, Tn.lt(y + p) ^ 1)
        B.c_and.clear(); B.c_xor.clear()
        if os.environ.get('INV_DEBUG'): print('STAGE', 'r1', B.nodes(), __import__('resource').getrusage(0).ru_maxrss//10**6, flush=True)
        r2 = B.AND(pre, (Tn.shr(n).eq(0)) ^ 1) if len(Tn.shr(n).bits) else ZERO
        B.c_and.clear(); B.c_xor.clear()
        out['first' if first else 'generic'] = dict(inv_T_lt_y_plus_p=(r1 == ZERO), Tn_lt_2n=(r2 == ZERO), low_bits_zero=(B.AND(pre, lowz) == ZERO))
    return out


def telescoping(n, w):
    """точная символьная проверка: 2^n*T_N = x*y + p*sum c_i 2^(w i), T_{i+1}=(T_i+x_i*y+c_i*p)/2^w, T_0=0 (дроби, мономы x_i*y и c_i*p)"""
    from fractions import Fraction
    N = n // w
    T = {}
    for i in range(N):
        new = {k: v / (1 << w) for k, v in T.items()}
        new[('xy', i)] = new.get(('xy', i), 0) + Fraction(1, 1 << w)
        new[('cp', i)] = new.get(('cp', i), 0) + Fraction(1, 1 << w)
        T = new
    lhs = {k: v * (1 << n) for k, v in T.items()}
    rhs = {}
    for i in range(N):
        rhs[('xy', i)] = Fraction(1 << (w * i)); rhs[('cp', i)] = Fraction(1 << (w * i))
    return lhs == rhs


def build_segs(n, w, p, build_fn=None):
    t0 = time.time()
    low, tr = lower_translate((build_fn or build)(n, w, p))
    print(f'n={n} w={w} p={p} nq={tr["nq"]} ops={len(tr["ops"])} lower {time.time()-t0:.1f}s', flush=True)
    top, order, wins, ltk, adk, xgk = segmentation(tr, low, n, w)
    ops = tr['ops']
    root = low['graph_events'][-1]['path'].split('/')[0]
    evd = {e['path']: e for e in low['graph_events']}
    rwire = ops[top[xgk][0]][1]
    segs = []
    for i, k in enumerate(wins):
        e = evd[root + '/' + k]
        ins, outs = e['inputs_msb'], e['outputs_msb']
        wi = {'x': ins['x'], 'y': ins['y'][::-1], 'T': ins['target'][::-1]}
        wo = {'x': ins['x'], 'y': ins['y'][::-1], 'T': outs['target'][::-1], 'm': outs['qrom_index'][::-1]}
        a0, b0 = top[k]
        if i == 0:
            a0 = top[xgk][0]; wo['r'] = [rwire]
            wi = {kk: v for kk, v in wi.items() if kk != 'T'}
        segs.append(('win%d' % i, 'win', ops[a0:b0], wi, wo, i == 0))
    e = evd[root + '/' + ltk]
    tw = e['inputs_msb']['x'][::-1]
    assert len(tw) == n
    segs.append(('lt', 'lt', ops[top[ltk][0]:top[ltk][1]], {'t': tw, 'r': [rwire]}, {'t': tw, 'r': [rwire]}, False))
    segs.append(('add', 'add', ops[top[adk][0]:top[adk][1]], {'t': tw, 'r': [rwire]}, {'t': tw, 'r': [rwire]}, False))
    return dict(low=low, tr=tr, top=top, order=order, wins=wins, segs=segs, rwire=rwire, tw=tw, ops=ops)
