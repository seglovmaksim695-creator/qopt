"""c12: распознавание MCX/ESOP-схем в .qc и синтез кандидатов.
  * лестница логических И в свежие анциллы (и чистые провода) с разборкой в глубину, снятие X-измерением + CZ-поправка;
  * варианты «последний уровень — настоящий Тоффоли на цель» (меньше анцилл);
  * конструкции с малым числом анцилл (Хаттар-Гидни 2024, Маслов 2015, Ито 2015) из Qiskit как кандидаты, с СВОЕЙ
    переписью внешней пары RCCX на чистой анцилле в И + измерение (−2 CX на анциллу) и своей CX-оптимизацией (свёртка фаз).
recognize(path, clean=False) -> dict(n, names, ins, anc, spec=[(цель, [моном=кортеж проводов, ...]), ...]) или None
candidates(path, clean=False) -> {(T, CX, extra): qasm_text}
Модель: clean=False — все провода .qc данные (любое состояние), |0> только у НОВЫХ кубитов;
        clean=True  — провода .qc не из строки .i стартуют в |0> и обязаны вернуться в |0>."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, random, itertools
sys.path[:0] = [_QPKG, (_QPKG + '/portfolio'), (_QPKG + '/contrib/c2')]
import qoptct, evo_lib, andc

MAXQ = 22
MAXMONO = 150  # больше мономов — это уже не «несколько MCX», лестница на них заведомо хуже общего компилятора c2
TAGS = {}


class Fail(Exception): pass


# ---------------------------------------------------------------- распознавание
def _mul(p, q):
    r = set()
    for a in p:
        for b in q:
            r ^= {a | b}
    return frozenset(r)


def anf(path, clean=False):
    """-> (names, n, ins_idx, anc_set, V): V[i] — АНФ провода i на выходе (frozenset масок мономов; маска 0 = константа 1).
    None, если схема не классическая перестановка (есть фазовые вентили)."""
    names, hl = qoptct.load_qc(path); n = len(names)
    try:
        lg = andc.logical(hl, n)
    except Exception:
        return None
    ins = andc.inputs_idx(path, names)
    anc = set(range(n)) - set(ins) if (clean and ins is not None) else set()
    V = [frozenset() if i in anc else frozenset({1 << i}) for i in range(n)]
    for g, q in lg:
        if g == 'x': V[q[0]] = V[q[0]] ^ frozenset({0})
        elif g == 'cx': V[q[1]] = V[q[1]] ^ V[q[0]]
        elif g == 'ccx': V[q[2]] = V[q[2]] ^ _mul(V[q[0]], V[q[1]])
        else: return None
    return names, n, ins, anc, V


def recognize(path, clean=False):
    """Функция = набор MCX/ESOP-членов: выход_i = вход_i ⊕ P_i(входы), остальные провода тождественны.
    Цели ставятся в порядке: j раньше i, если P_j зависит от x_i (читает ещё не изменённый вход). Цикл -> None."""
    r = anf(path, clean)
    if r is None: return None
    names, n, ins, anc, V = r
    D = {}
    for i in range(n):
        if i in anc:
            if V[i]: return None
            continue
        d = V[i] ^ frozenset({1 << i})
        if d: D[i] = d
    if not D: return None
    dep = {j: {i for m in D[j] for i in D if (m >> i) & 1} for j in D}
    if any(j in dep[j] for j in D): return None
    order, left = [], set(D)
    while left:
        cand = [j for j in left if not any(j in dep[k] for k in left if k != j)]
        if not cand: return None
        j = min(cand); order.append(j); left.discard(j)
    spec = []
    for j in order:
        monos = sorted((tuple(i for i in range(n) if (m >> i) & 1) for m in D[j]), key=lambda t: (len(t), t))
        spec.append((j, monos))
    return dict(n=n, names=names, ins=ins, anc=sorted(anc), spec=spec)


# ---------------------------------------------------------------- лестница И
def ccx_ops(a, b, t):
    return list(evo_lib.lower_gate('ccx', (a, b, t)))


def ladder(rec, max_extra, direct2=False, deep='and', shuffle=None, free=None):
    """Лестница И: узел = И(родитель, переменная) в свежую анциллу/чистый провод (4 T, 3 CX), копия в цель CX,
    разборка в глубину: H+измерение, при 1 — CZ(родитель, переменная), сброс (1 CX, 0 T). Узлы с общим префиксом
    (разделяемые мономы) считаются один раз.
    deep='ccx': последний уровень — настоящий Тоффоли на цель (7 T, 6 CX, на анциллу меньше);
    direct2: двухконтрольные мономы — Тоффоли напрямую."""
    n = rec['n']; ops = []; pool = list(free or []); state = {'nxt': n}
    lim = n + max_extra

    def alloc():
        if pool: return pool.pop(0)
        if state['nxt'] >= lim: raise Fail('анцилл не хватает')
        state['nxt'] += 1; return state['nxt'] - 1

    for t, monos in rec['spec']:
        freq = {}
        for m in monos:
            for v in m: freq[v] = freq.get(v, 0) + 1
        if shuffle is None: key = lambda v: (-freq[v], v)
        else:
            rk = {v: shuffle.random() for v in freq}; key = lambda v: (-freq[v] + 1.5 * rk[v], v)
        trie = {}
        for m in monos:
            node = trie
            for v in sorted(m, key=key): node = node.setdefault(v, {})
            node[None] = True

        def visit(node, wire, depth):
            if node.get(None):
                ops.append(('x', (t,)) if depth == 0 else ('cx', (wire, t)))
            for v, c in [(v, c) for v, c in node.items() if v is not None]:
                if depth == 0: visit(c, v, 1); continue
                leaf = list(c.keys()) == [None]
                if leaf and ((depth == 1 and direct2) or (depth >= 2 and deep == 'ccx')):
                    ops.extend(ccx_ops(wire, v, t)); continue
                w = alloc()
                ops.append(('AND', (wire, v, w)))
                visit(c, w, depth + 1)
                ops.append(('MZ', (w, (('cz', (wire, v)), ('x', (w,))))))
                pool.append(w)
        visit(trie, None, 0)
    return ops, state['nxt']


# ---------------------------------------------------------------- Qiskit-кандидаты (Хаттар-Гидни и др.)
CT = ('h', 's', 'sdg', 't', 'tdg', 'x', 'z', 'cx', 'cz')


def _qiskit_ops(sub, qmap):
    out = []
    for ins in sub.data:
        nm = ins.operation.name
        if nm not in CT: raise Fail('вентиль ' + nm)
        out.append((nm, tuple(qmap[sub.find_bit(b).index] for b in ins.qubits)))
    return out


def RCCX(a, b, c):
    return [('h', (a,)), ('t', (a,)), ('cx', (c, a)), ('tdg', (a,)), ('cx', (b, a)), ('t', (a,)),
            ('cx', (c, a)), ('tdg', (a,)), ('h', (a,))]


def erase_rewrite(ops, cleanq):
    """Для чистых анцилл a: первый блок RCCX(b,c->a) -> И (+sdg), последний -> H+измерение+CZ(b,c)+сброс. −2 CX на анциллу.
    (Корректно: после внешнего RCCX середина действует точно, а a = b·c, поэтому снятие измерением точно.)"""
    ops = list(ops)
    for a in cleanq:
        idx = [i for i, (g, q) in enumerate(ops) if g not in ('MZ', 'AND') and a in q]
        if len(idx) < 18: continue
        f, l = idx[:9], idx[-9:]
        if f != list(range(f[0], f[0] + 9)) or l != list(range(l[0], l[0] + 9)): continue
        blk = ops[f[0]:f[0] + 9]
        cxs = [q for g, q in blk if g == 'cx']
        if len(cxs) != 3: continue
        c, b = cxs[0][0], cxs[1][0]
        if blk != RCCX(a, b, c) or ops[l[0]:l[0] + 9] != blk: continue
        ops[l[0]:l[0] + 9] = [('MZ', (a, (('cz', (b, c)), ('x', (a,)))))]
        ops[f[0]:f[0] + 9] = [('AND', (b, c, a))]
    return ops


def qiskit_variants(k):
    import qiskit.synthesis as S
    fns = {'n_clean_m15': (S.synth_mcx_n_clean_m15, True), '1_clean_kg24': (S.synth_mcx_1_clean_kg24, True),
           '2_clean_kg24': (S.synth_mcx_2_clean_kg24, True), 'n_dirty_i15': (S.synth_mcx_n_dirty_i15, False),
           '1_dirty_kg24': (S.synth_mcx_1_dirty_kg24, False), '2_dirty_kg24': (S.synth_mcx_2_dirty_kg24, False)}
    res = {}
    for nm, (f, cl) in fns.items():
        try:
            sub = f(k)
        except Exception:
            continue
        if any(i.operation.name not in CT for i in sub.data): continue  # вращения не Clifford+T — не берём
        res[nm] = (sub, cl)
    return res


_VC = {}


def qiskit_mcx(rec, vname, max_extra, free=None, rewrite=True, use_clean_wires=False, pre=0):
    """Каждый моном — один MCX готовым синтезом. Анциллы: «грязные» — любые провода вне моном∪цель (данные в любом состоянии);
    «чистые» — свежие кубиты (или чистые провода .qc при use_clean_wires)."""
    n = rec['n']; ops = []; extra_used = 0
    for t, monos in rec['spec']:
        for m in monos:
            k = len(m)
            if k == 0: ops.append(('x', (t,))); continue
            if k == 1: ops.append(('cx', (m[0], t))); continue
            if k == 2: ops.extend(ccx_ops(m[0], m[1], t)); continue
            chain = []; head = []; tail = []
            m = list(m)
            if pre >= 2 and k - pre + 1 >= 3:
                # префикс: И первых pre контролей лестницей в свежие анциллы (4(pre-1) T), далее MCX на (узел + остаток)
                chain = [n + i for i in range(pre - 1)]
                head.append(('AND', (m[0], m[1], chain[0])))
                for j in range(2, pre): head.append(('AND', (chain[j - 2], m[j], chain[j - 1])))
                for j in range(pre - 1, 0, -1):
                    par = m[0] if j == 1 else chain[j - 2]
                    tail.append(('MZ', (chain[j - 1], (('cz', (par, m[j])), ('x', (chain[j - 1],))))))
                mm = [chain[-1]] + m[pre:]
            elif pre >= 2: raise Fail('префикс не нужен')
            else: mm = m
            kk = len(mm)
            if kk not in _VC: _VC[kk] = qiskit_variants(kk)
            if vname not in _VC[kk]: raise Fail('нет варианта')
            sub, cl = _VC[kk][vname]; need = sub.num_qubits - kk - 1
            base = len(chain)
            if cl:
                pool = [w for w in (free or []) if w not in m and w != t] if use_clean_wires else []
                ancs = pool[:need]; fresh = need - len(ancs)
                if base + fresh > max_extra: raise Fail('мало анцилл')
                extra_used = max(extra_used, base + fresh)
                ancs = ancs + [n + base + i for i in range(fresh)]; cleanq = list(ancs)
            else:
                others = [w for w in range(n) if w not in m and w != t]
                if need > len(others): raise Fail('мало грязных проводов')
                if base > max_extra: raise Fail('мало анцилл')
                extra_used = max(extra_used, base)
                ancs = others[:need]; cleanq = []
            part = _qiskit_ops(sub, list(mm) + [t] + ancs)
            if rewrite and cleanq: part = erase_rewrite(part, cleanq)
            ops.extend(head + part + tail)
    return ops, n + extra_used


# ---------------------------------------------------------------- кандидаты
def emit(ops, n, out, tag):
    ops, N2 = andc.compact(ops, n)
    for fold in (False, True):
        try:
            txt = andc.to_qasm(ops, N2, fold=fold)
        except Exception:
            continue
        T, CX = andc.count(txt); k = (T, CX, N2 - n)
        if k not in out:
            out[k] = txt; TAGS[k] = tag + ('+fold' if fold else '')


def candidates(path, clean=False, seeds=4, use_qiskit=True):
    rec = recognize(path, clean)
    out = {}
    if rec is None: return out
    n = rec['n']; lim = MAXQ - n
    if sum(len(ms) for _, ms in rec['spec']) > MAXMONO or n > MAXQ: return out
    free = rec['anc'] if clean else None
    for ma in sorted({lim, 6, 5, 4, 3, 2, 1, 0}):
        if ma > lim: continue
        for d2, deep in itertools.product((False, True), ('and', 'ccx')):
            for s in range(seeds):
                try:
                    ops, N = ladder(rec, ma, d2, deep, shuffle=None if s == 0 else random.Random(s), free=free)
                except Fail:
                    continue
                emit(ops, n, out, f'ladder(ma={ma},d2={d2},{deep},s={s})')
    if use_qiskit:
        for vn in ('n_clean_m15', '1_clean_kg24', '2_clean_kg24', 'n_dirty_i15', '1_dirty_kg24', '2_dirty_kg24'):
            for rw in (True, False):
                for ucw in ((False, True) if clean else (False,)):
                    for pre in (0, 2, 3, 4, 5, 6):
                        try:
                            ops, N = qiskit_mcx(rec, vn, lim, free=free, rewrite=rw, use_clean_wires=ucw, pre=pre)
                        except Fail:
                            continue
                        emit(ops, n, out, f'qiskit {vn} rewrite={rw} cleanwires={ucw} pre={pre}')
    return out
