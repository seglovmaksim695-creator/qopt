"""Proof of the generic Montgomery window by its Qualtran sub-nodes (ported from c2_split of CFIX).

A window = w x ``CAdd`` + ``Xor`` + [``QROAMClean`` + ``Add`` + ``QROAMCleanAdjoint``].  Each sub-node is proved
separately and the composition is checked against the contract of the window (BDD, or interval arithmetic plus exact
integer runs for w >= 4).  The CFIX checks are part of this module and are NOT optional:

* coverage: every operation of the window belongs to exactly one sub-node, the children of the window tile its span
  with no gaps and no overlaps (``tile_info``, ``window_checks``);
* freshness: the wires that a sub-node assumes to be zero are neither inputs of the window nor touched before.
"""
import os
import time
import random
from hashlib import sha256

from ..bridge import sim
from ..bridge.emit import emit_qasm
from ..bridge.strict_proof import prove_program, digest_json
from qoptct.verify.bdd import BDD, ZERO, ONE
from qoptct.verify.ksym import SymInt, SpecContext, lookup
from ._core import (compact, make_layout, make_win, support, bits_set, win_eval)

def children_info(low, tr, top, wins, order):
    root = low['graph_events'][-1]['path'].split('/')[0]
    evd = {e['path']: e for e in low['graph_events']}
    byw = {}
    for t in tr['trace']:
        ps = t['leaf']['path'].split('/')
        if len(ps) < 3 or ps[1] not in wins: continue
        d = byw.setdefault(ps[1], {})
        k = ps[2]
        if k not in d: d[k] = [t['logical_start'], t['logical_end']]
        else: d[k][1] = t['logical_end']
    res = []
    for wk in wins:
        ch = []
        for k, (a, b) in sorted(byw[wk].items(), key=lambda kv: kv[1][0]):
            nm = k.split(':', 1)[1]
            e = evd[root + '/' + wk + '/' + k]
            ch.append(dict(name=nm, a=a, b=b, ins=e['inputs_msb'], outs=e['outputs_msb']))
        e = evd[root + '/' + wk]
        res.append(dict(key=wk, children=ch, wins_in=e['inputs_msb'], wins_out=e['outputs_msb'], span=top[wk]))
    return res


def build_sub(P, ch_list, ops):
    """ch_list: дети окна. -> список сегментов: (label, kind, ops, wi, wo, descr)"""
    segs = []
    c = ch_list
    names = [x['name'] for x in c if x['name'] not in ('Join', 'Split')]
    w_ = P.w
    assert names == ['CAdd'] * w_ + ['Xor', 'QROAMClean', 'Add', 'QROAMCleanAdjointWrapper'], names
    real = [x for x in c if x['name'] not in ('Join', 'Split')]
    for i in range(w_):
        x = real[i]
        wi = {'c': x['ins']['ctrl'], 'a': x['ins']['a'][::-1], 'b': x['ins']['b'][::-1]}
        wo = {'c': x['outs']['ctrl'], 'a': x['outs']['a'][::-1], 'b': x['outs']['b'][::-1]}
        segs.append(('cadd%d' % i, 'cadd', ops[x['a']:x['b']], wi, wo))
    x = real[w_]
    segs.append(('xor', 'xor', ops[x['a']:x['b']], {'x': x['ins']['x'][::-1], 'y': x['ins']['y'][::-1]}, {'x': x['outs']['x'][::-1], 'y': x['outs']['y'][::-1]}))
    q0, ad, q1 = real[w_ + 1], real[w_ + 2], real[w_ + 3]
    a0, b0 = q0['a'], q1['b']
    assert q0['b'] <= ad['a'] and ad['b'] <= q1['a']
    wi = {'m': q0['ins']['selection'][::-1], 'b': ad['ins']['b'][::-1]}
    wo = {'m': q1['outs']['selection'][::-1], 'b': ad['outs']['b'][::-1]}
    segs.append(('qaa', 'qaa', ops[a0:b0], wi, wo))
    return segs


def sub_ranges(win):
    """CFIX: диапазоны операций подузлов окна (глобальные индексы): CAdd*w, Xor, [QROAMClean..Adjoint]"""
    real = [x for x in win['children'] if x['name'] not in ('Join', 'Split')]
    w_ = sum(1 for x in real if x['name'] == 'CAdd')
    rs = [('cadd%d' % i, real[i]['a'], real[i]['b']) for i in range(w_)]
    rs.append(('xor', real[w_]['a'], real[w_]['b']))
    rs.append(('qaa', real[w_ + 1]['a'], real[w_ + 3]['b']))
    return rs


def tile_info(win):
    """CFIX: дети окна (включая нулевые Join/Split) покрывают отрезок окна ровно: без дыр, пересечений и выходов за границы"""
    notes = []
    a0, b0 = win['span']
    for x in win['children']:
        if x['name'] in ('Join', 'Split') and x['b'] != x['a']: notes.append('%s has ops' % x['name'])
        if x['b'] < x['a']: notes.append('child %s negative length' % x['name'])
    pos = a0
    for x in sorted(win['children'], key=lambda c: (c['a'], c['b'])):
        if x['a'] != pos:
            notes.append('children: gap/overlap before %s at %d (expected %d)' % (x['name'], x['a'], pos))
        pos = max(pos, x['b']) if x['a'] == pos else x['b']
    if pos != b0: notes.append('children end at %d, window span ends at %d' % (pos, b0))
    return notes


def window_checks(P, win, wops, off, segs):
    """CFIX (закрытие дыры AUD-C-01), вызывается после build_sub:
    (1) покрытие: длина операций окна == его отрезок; дети тайлят отрезок; отрезки подузлов тайлят отрезок;
        каждая операция окна принадлежит РОВНО ОДНОМУ подузлу; операции каждого подузла (после среза) — ровно его диапазон;
    (2) свежесть: идём по подузлам в порядке операций; НЕвходные провода подузла (он их считает нулевыми) не входы окна
        (ни один ключ wins_in) и не тронуты ранее (ни окном, ни прежними подузлами)."""
    notes = []
    a0, b0 = win['span']
    if len(wops) != b0 - a0 or off != a0: notes.append('window ops %d (off %d) != span %d..%d' % (len(wops), off, a0, b0))
    notes += tile_info(win)
    cnt = [0] * (b0 - a0)
    rng_ = sub_ranges(win)
    for lab, a, b in rng_:
        if a < a0 or b > b0 or b < a:
            notes.append('%s range %d..%d outside window %d..%d' % (lab, a, b, a0, b0)); continue
        for i in range(a - a0, b - a0): cnt[i] += 1
    miss = [i + a0 for i, c in enumerate(cnt) if c == 0]
    dup = [i + a0 for i, c in enumerate(cnt) if c > 1]
    if miss: notes.append('%d window ops belong to no sub-node (first: %s)' % (len(miss), miss[:4]))
    if dup: notes.append('%d window ops belong to >1 sub-nodes (first: %s)' % (len(dup), dup[:4]))
    segd = {lab: (a, b) for lab, a, b in rng_}
    for (lab, kind, o, swi, swo) in segs:
        a, b = segd[lab]
        if len(o) != b - a: notes.append('%s: %d ops cut, range has %d' % (lab, len(o), b - a))
    wi = win['wins_in']
    touched = set(q for v in wi.values() for q in v)
    for (lab, kind, o, swi, swo) in segs:
        ins = {q for v in swi.values() for q in v}
        sup = support(o, swi, swo)
        bad = (sup - ins) & touched
        if bad:
            notes.append('%s: non-input (assumed zero) wires of sub-node are window inputs or already touched: %s' % (lab, sorted(bad)[:4]))
        touched |= sup
    return dict(ok=not notes, notes=notes[:8], n_subs=len(segs), n_ops=b0 - a0)


def apply_sub(P, kind, S_B, rd, wi, wo):
    """возвращает {wire: edge} выходов по входным регистрам (rd(wires)->SymInt); для spec сегмента и композиции"""
    B = S_B
    out = {}
    if kind == 'cadd':
        c = rd(wi['c']).bits[0]; a = rd(wi['a']); b = rd(wi['b'])
        am = SymInt(B, [B.AND(c, t) for t in a.bits]).resize(len(b.bits))
        r = (b + am).resize(len(b.bits))
        for q, e in zip(wo['b'], r.bits): out[q] = e
        for q, e in zip(wo['a'], a.bits): out[q] = e
        out[wo['c'][0]] = c
    elif kind == 'xor':
        x = rd(wi['x']); y = rd(wi['y'])
        for q, e1, e2 in zip(wo['y'], x.bits, y.bits): out[q] = B.XOR(e1, e2)
        for q, e in zip(wo['x'], x.bits): out[q] = e
    else:
        m = rd(wi['m']); b = rd(wi['b'])
        r = (b + lookup(B, m, P.table, len(b.bits))).resize(len(b.bits))
        for q, e in zip(wo['b'], r.bits): out[q] = e
        for q, e in zip(wo['m'], m.bits): out[q] = e
    return out


def sim_sub(P, kind, nq, ops, wi, wo, rng, count):
    cases = []
    bad = 0
    for _ in range(count):
        if kind == 'cadd':
            nb = len(wi['b']); a = rng.randrange(1 << len(wi['a'])); b = rng.choice([0, (1 << nb) - 1, rng.randrange(1 << nb)]); c = rng.randrange(2)
            init = bits_set(a, wi['a']) | bits_set(b, wi['b']) | bits_set(c, wi['c'])
            exp = bits_set(a, wo['a']) | bits_set((b + c * a) % (1 << nb), wo['b']) | bits_set(c, wo['c'])
        elif kind == 'xor':
            x = rng.randrange(1 << len(wi['x'])); y = rng.randrange(1 << len(wi['y']))
            init = bits_set(x, wi['x']) | bits_set(y, wi['y']); exp = bits_set(x, wo['x']) | bits_set(x ^ y, wo['y'])
        else:
            nb = len(wi['b']); m = rng.randrange(P_N(P)); b = rng.choice([0, (1 << nb) - 1, rng.randrange(1 << nb), rng.randrange(2 * P.p)])
            init = bits_set(m, wi['m']) | bits_set(b, wi['b']); exp = bits_set(m, wo['m']) | bits_set((b + P.table[m]) % (1 << nb), wo['b'])
        val, ph = sim.run(nq, ops, init)
        if val != exp or any(ph): bad += 1
    return count, bad


def P_N(P): return 1 << P.w


def run_sub(label, kind, P, ops_g, wi_g, wo_g, rng, budget, max_mb, tag, cache):
    nq, cops, wi, wo, mp = compact(ops_g, wi_g, wo_g)
    layout = make_layout(nq, wi, wo, cops)
    qasm = emit_qasm(nq, cops)
    formula = {'cadd': "b'=(b+c*a) mod 2^len(b), a,c unchanged", 'xor': "y'=y xor x, x unchanged",
               'qaa': "m unchanged, b'=(b+c(m)*p) mod 2^len(b); QROAM target wires are fresh and consumed (erased) inside"}[kind]
    contract = dict(source='C2 split sub-segment', kind=kind, n=P.n, w=P.w, p=P.p, qualtran=tag, formula=formula)
    cnt, bad = sim_sub(P, kind, nq, cops, wi, wo, rng, 12)
    key = (sha256(qasm.encode()).hexdigest(), digest_json(layout), digest_json(contract))
    if key in cache and not bad:
        r = dict(cache[key]); r.update(label=label, dedup_of=cache[key]['label'], sim_cases=cnt, sim_bad=bad); return r
    if bad: return dict(label=label, kind=kind, status='SIM_MISMATCH', sim_cases=cnt, sim_bad=bad, nq=nq, ops=len(cops))
    t0 = time.time()
    receipt = dict(status='PASS', qasm_sha256=key[0], contract_sha256=digest_json(contract), layout_sha256=digest_json(layout), independent_inputs=cnt, all_wires_checked=True, phase_checked=True)

    def spec(S):
        return apply_sub(P, kind, S.B, S.reg, wi, wo)
    v = prove_program(qasm, layout, spec, pre=None, validation=receipt, contract=contract, budget=budget, max_mb=max_mb)
    r = dict(label=label, kind=kind, status=v['status'], detail=(v.get('detail') or '')[:160], nq=nq, ops=len(cops), sim_cases=cnt, sim_bad=bad,
             wall=round(time.time() - t0, 2), peak_mb=v.get('peak_rss_mb'), qasm_sha256=key[0], contract_sha256=receipt['contract_sha256'],
             domain=(v.get('partial_domains') or {}).get('status'))
    cache[key] = r
    return r


def compose_check(P, win, subs_w):
    """BDD: при pre окна композиция подконтрактов == контракт окна (по всем выходным проводам x,y,T,m)"""
    n, w, p = P.n, P.w, P.p
    wi = {'x': win['wins_in']['x'], 'y': win['wins_in']['y'][::-1], 'T': win['wins_in']['target'][::-1]}
    wo = {'x': wi['x'], 'y': wi['y'], 'T': win['wins_out']['target'][::-1], 'm': win['wins_out']['qrom_index'][::-1]}
    B = BDD()
    xs = [B.new_var()[1] for _ in range(w)]
    ye, Te = [], []
    for i in range(n + w):
        if i < n: ye.append(B.new_var()[1])
        Te.append(B.new_var()[1])
    edge = {}
    for q, e in zip(wi['x'], xs): edge[q] = e
    for q, e in zip(wi['y'], ye): edge[q] = e
    for q, e in zip(wi['T'], Te): edge[q] = e
    S = SpecContext(B, lambda q: edge[q], set())
    pre, spec = make_win(P, False, wi, wo)
    import resource
    dbg = lambda t: print('CSTAGE', t, B.nodes(), resource.getrusage(0).ru_maxrss // 10**6, flush=True) if os.environ.get('INV_DEBUG') else None
    prec = pre(S); dbg('pre'); want = spec(S); dbg('want')
    state = dict(edge)

    def squeeze():   # сборка мусора BDD (экономия памяти), корни: state, want, prec
        nonlocal prec, want, state
        rm = B.gc([prec] + list(want.values()) + list(state.values()))
        prec = rm(prec); want = {q: rm(e) for q, e in want.items()}; state = {q: rm(e) for q, e in state.items()}
    squeeze()
    rd = lambda ws: SymInt(B, [state.get(q, ZERO) for q in ws])
    for (label, kind, o, swi, swo) in subs_w:
        res = apply_sub(P, kind, B, rd, swi, swo)
        state.update(res); dbg('sub ' + label)
        squeeze(); dbg('gc')
    bad = []
    for q in set(q for k in ('x', 'y', 'T', 'm') for q in wo[k]):
        got = state.get(q, ZERO)
        d = B.AND(prec, B.XOR(got, want[q]))
        if d != ZERO: bad.append(q)
    return dict(ok=not bad, bad_wires=bad[:8], bdd_nodes=B.nodes())


def compose_analytic(P, win, subs_w, rng, first=False):
    """w>=4: BDD-композиция не помещается в 3 ГБ. Вместо неё: (1) проверка кодом положения проводов (срезы T, y, x, m), (2) точная целочисленная
    проверка отсутствия переполнения срезов (интервальная арифметика на pre: y<p, T<y+p), (3) прогон композиции подконтрактов на целых числах против контракта окна."""
    n, w, p = P.n, P.w, P.p
    notes = []
    Tw = win['wins_in']['target'][::-1]; yw = win['wins_in']['y'][::-1]; xw = win['wins_in']['x']
    wo_T = win['wins_out']['target'][::-1]; wo_m = win['wins_out']['qrom_index'][::-1]
    segd = {lab: (kind, swi, swo) for (lab, kind, o, swi, swo) in subs_w}
    for k in range(w):
        kind, swi, swo = segd['cadd%d' % k]
        if swi['c'] != [xw[k]]: notes.append(f'cadd{k}: ctrl wire != x[{k}]')
        if swi['a'] != yw: notes.append(f'cadd{k}: a != y')
        if swi['b'] != Tw[k:k + n + 1]: notes.append(f'cadd{k}: b != T[{k}:{k + n + 1}]')
        if swo['b'] != swi['b'] or swo['a'] != swi['a']: notes.append(f'cadd{k}: outputs wires differ from inputs')
    kind, swi, swo = segd['xor']
    if swi['x'] != Tw[:w]: notes.append('xor: x != T[:w]')
    if swo['y'] != wo_m or swi['y'] != wo_m: notes.append('xor: y != qrom_index wires')
    kind, swi, swo = segd['qaa']
    if swi['m'] != wo_m or swo['m'] != wo_m: notes.append('qaa: selection != qrom_index wires (copy of m)')
    if swi['b'] != Tw[:n + w] or swo['b'] != Tw[:n + w]: notes.append('qaa: b != T[:n+w]')
    if wo_T[:n] != Tw[w:n + w]: notes.append('window out T[:n] != T[w:]')
    if wo_T[n:] != Tw[:w]: notes.append('window out T high wires != old low wires T[:w] (cyclic shift)')
    # (2) границы без переполнения срезов
    N_ = 1 << w
    Tmax = 0 if first else 2 * p - 2   # T <= y+p-1 <= 2p-2 (первое окно: T=0)
    for k in range(w):
        Tk = Tmax + (2 ** (k + 1) - 1) * (p - 1)
        if not Tk < 2 ** (n + k + 1): notes.append(f'overflow possible in cadd{k}')
    if not Tmax + (N_ - 1) * (p - 1) + max(P.table) < 2 ** (n + w): notes.append('overflow possible in qaa')
    # (3) целочисленная композиция против контракта окна
    cases = 0
    for _ in range(3000):
        y = rng.choice([0, 1, p - 1, rng.randrange(p)]); T = 0 if first else rng.choice([0, y + p - 1, p, rng.randrange(y + p)]); xd = rng.choice([0, N_ - 1, rng.randrange(N_)])
        b = T
        for k in range(w):
            sl = (b >> k) & ((1 << (n + 1)) - 1); sl = (sl + ((xd >> k) & 1) * y) % (1 << (n + 1))
            b = (b & ~(((1 << (n + 1)) - 1) << k)) | (sl << k)
        m = b & (N_ - 1)
        b = (b + P.table[m]) % (1 << (n + w))
        mm, Tn = win_eval(P, xd, y, T)
        if m != mm or (b >> w) != Tn or (b & (N_ - 1)) != 0: notes.append(f'composition != contract at x={xd} y={y} T={T}'); break
        cases += 1
    return dict(ok=not notes, notes=notes[:6], method='analytic_wires+interval+int_sim', int_cases=cases)
