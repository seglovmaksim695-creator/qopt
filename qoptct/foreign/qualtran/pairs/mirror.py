"""Mirror matching of the backward program with the forward one (commuting neighbours allowed).

A program = a list of operations (kind, a, b, t, l), kind: 0 x, 1 cx(a->b), 2 and(a,b->t), 3 erase(t; a,b, l).
The mirror of a forward program P: reverse order; x->x; cx->cx; and(a,b->t)->erase(a,b;t); erase(a,b;t,l)->and(a,b->t)[+cx(l->t)].
Matching: every operation b of the backward program is looked up among the not yet matched operations m of the mirror so that
  (1) m and b have the same kind and their wires agree with the renaming pi (injective both ways, pi = id on the joint wires);
  (2) m commutes with all not yet matched mirror operations that stand EARLIER than it.
Commutation (exact, for the classical action of the gates): the shared wires of two operations are, in both, only "read"
(controls of cx, a and b of and/erase, the linear correction l), or, in both, an "xor target" (target of cx or x).  Everything
else does not commute.
Success = the backward program is obtained from the mirror by swapping pairs of commuting neighbours and renaming by pi.
This is enough for the identity of the pair: classical values on the wires of every operation are preserved.
"""
from collections import deque

K_X, K_CX, K_AND, K_ERASE = 0, 1, 2, 3
K_LD, K_CLR = 10, 11     # абстрактные: загрузка таблицы QROM (t ^= b0^s*(b0^b1)) и её снятие; код (b0,b1) в поле t


def roles(op):
    k, a, b, t, l = op
    if k == K_X:
        return ((a, 'x'),)
    if k == K_CX:
        return ((a, 'r'), (b, 'x'))
    if k == K_LD:
        return ((a, 'r'), (b, 'x'))
    if k == K_CLR:
        return ((a, 'r'), (b, 'w'))
    r = ((a, 'r'), (b, 'r'), (t, 'w'))
    return r + ((l, 'r'),) if l >= 0 else r


def commute(o1, o2):
    r2 = dict(roles(o2))
    for w, ro in roles(o1):
        if w in r2:
            ro2 = r2[w]
            if not (ro == ro2 and ro in ('r', 'x')):
                return False
    return True


def mirror(op):
    k, a, b, t, l = op
    if k == K_X:
        return [(K_X, a, 0, 0, -1)]
    if k == K_CX:
        return [(K_CX, a, b, 0, -1)]
    if k == K_LD:
        return [(K_CLR, a, b, t, l)]
    if k == K_CLR:
        return [(K_LD, a, b, t, l)]
    if k == K_AND:
        return [(K_ERASE, a, b, t, -1)]
    r = [(K_AND, a, b, t, -1)]
    if l >= 0:
        r.append((K_CX, l, t, 0, -1))
    return r


def wires(op):
    k, a, b, t, l = op
    if k == K_X:
        return (a,)
    if k in (K_CX, K_LD, K_CLR):
        return (a, b)
    return (a, b, t) if l < 0 else (a, b, t, l)


def mirror_prog(ops):
    out = []
    for op in reversed(ops):
        out.extend(mirror(op))
    return out


def fingerprints(prog, rounds=80, anchors=None):
    anchors = anchors or {}
    """провод -> цвет: инвариант перенумерации и коммутации (уточнение Вейсфейлера-Лемана по мультимножествам);
    только подсказка для выбора пары, на корректность проверки не влияет (проверяется каждая операция)."""
    inc = {}
    for op in prog:
        rl = roles(op)
        for w, ro in rl:
            inc.setdefault(w, []).append((op[0], ro, rl))
    col = {w: hash(tuple(sorted((k, ro) for k, ro, _ in v))) for w, v in inc.items()}
    for q, lab in anchors.items():
        col[q] = hash(('anchor', lab))
    ncls = len(set(col.values()))
    for _ in range(rounds):
        new = {}
        for w, v in inc.items():
            items = []
            for k, ro, rl in v:
                items.append((k, ro, tuple(sorted((r2, col.get(w2)) for w2, r2 in rl if w2 != w))))
            new[w] = hash(tuple(sorted(items, key=repr)))
        for q, lab in anchors.items():
            new[q] = hash(('anchor', lab))
        col = new
        n2 = len(set(col.values()))
        if n2 == ncls:
            break
        ncls = n2
    return col


def build_items(prog, scratch=()):
    """список элементов: одиночная операция или БЛОК на проводе-scratch моста (AND ...->t ... ERASE t, все операции между
    ними). Блок возвращает scratch в 0 и действует только на прочие провода; роли блока: 'r' если провод только читается,
    'x' если только xor-мишень, иначе 'm' (не коммутирует). Элемент = (ops(tuple), roles(dict))."""
    scr = set(scratch)
    items, i, N = [], 0, len(prog)
    while i < N:
        op = prog[i]
        if op[0] == K_AND and op[3] in scr:
            t = op[3]
            j = i + 1
            while j < N and not (prog[j][0] == K_ERASE and prog[j][3] == t):
                if prog[j][0] in (K_AND, K_ERASE) and prog[j][3] == t:
                    raise ValueError('вложенный AND на scratch в позиции %d' % j)
                j += 1
            if j == N:
                raise ValueError('AND на scratch без erase, позиция %d' % i)
            blk = tuple(prog[i:j + 1])
            rl = {}
            for o in blk:
                for w, ro in roles(o):
                    if w == t:
                        continue
                    rl.setdefault(w, set()).add(ro)
            net = {w: (next(iter(v)) if len(v) == 1 and next(iter(v)) in ('r', 'x') else 'm') for w, v in rl.items()}
            items.append((blk, net))
            i = j + 1
        else:
            items.append(((op,), {w: ('w' if ro == 'w' else ro) for w, ro in roles(op)}))
            i += 1
    return items


def commute_items(r1, r2):
    for w, ro in r1.items():
        if w in r2:
            ro2 = r2[w]
            if not (ro == ro2 and ro in ('r', 'x')):
                return False
    return True


def match(fwd, bwd, join, window=4096, pi0=None, scratch=(), scratch_b=None):
    """join: провода с pi = id (стык). pi0: дополнительные заданные пары (прямой провод -> обратный).
    scratch / scratch_b: провод(а) scratch моста в прямой / обратной программе (pi между ними задаётся автоматически)."""
    sb = list(scratch if scratch_b is None else scratch_b)
    return match_items(build_items(mirror_prog(fwd), scratch), build_items(list(bwd), sb), join, window, pi0, scratch, sb)


def match_items(mir_items, bwd_items, join, window=4096, pi0=None, scratch=(), scratch_b=()):
    """-> dict(ok, pi, msg, stats). join: провода стыка (pi = id на них). Блоки сопоставляются поэлементно по порядку внутри."""
    st = dict(nb=len(bwd_items), seeded=0, reordered=0, ab_swaps=0, max_window=0, blocks=sum(1 for it in bwd_items if len(it[0]) > 1))
    if len(mir_items) != len(bwd_items):
        return dict(ok=False, pi=None, msg='число элементов: зеркало %d, обратная %d' % (len(mir_items), len(bwd_items)), stats=st)
    mir_ops = [o for it in mir_items for o in it[0]]
    bwd_ops = [o for it in bwd_items for o in it[0]]
    pi = {q: q for q in join}
    if pi0:
        pi.update(pi0)
    for a, b in zip(scratch, scratch_b):
        pi[a] = b
    inv = {v: u for u, v in pi.items()}
    assert len(inv) == len(pi), 'заданная перенумерация не инъективна'
    am = {u: ('s', u) for u in pi}
    ab = {v: ('s', u) for u, v in pi.items()}
    fm, fb = fingerprints(mir_ops, anchors=am), fingerprints(bwd_ops, anchors=ab)
    # провода с единственным цветом в обеих программах: пара определена однозначно (подсказка; итог всё равно проверяется поэлементно)
    cm, cb = {}, {}
    for w, c in fm.items():
        cm.setdefault(c, []).append(w)
    for w, c in fb.items():
        cb.setdefault(c, []).append(w)
    nseed = 0
    for c, lu in cm.items():
        lv = cb.get(c)
        if lv and len(lu) == 1 and len(lv) == 1 and lu[0] not in pi and lv[0] not in inv:
            pi[lu[0]] = lv[0]
            inv[lv[0]] = lu[0]
            nseed += 1

    def fits_ops(mo, bo, newmap, use_fp=True):
        """операции mo (зеркало) и bo (обратная) согласуются с pi + newmap; допускается a<->b у AND/ERASE."""
        if mo[0] != bo[0]:
            return None
        if mo[0] >= K_LD and mo[3] != bo[3]:
            return None
        wm, wb = wires(mo), wires(bo)
        if len(wm) != len(wb):
            return None
        for sw in (0, 1):
            if sw and mo[0] not in (K_AND, K_ERASE):
                break
            w2 = ((wm[1], wm[0]) + wm[2:]) if sw else wm
            loc, linv, ok = {}, {}, True
            for u, v in zip(w2, wb):
                cu = newmap.get(u, loc.get(u, pi.get(u)))
                cv = newmap_inv.get(v, linv.get(v, inv.get(v)))
                if (cu is not None and cu != v) or (cv is not None and cv != u):
                    ok = False
                    break
                if use_fp and cu is None and fm.get(u) != fb.get(v):
                    ok = False
                    break
                loc[u] = v
                linv[v] = u
            if ok:
                return loc, linv, sw
        return None

    pend = deque()
    nxt = 0
    N = len(mir_items)
    st['seeded'] = nseed
    for i, bit in enumerate(bwd_items):
        while nxt < N and len(pend) < window:
            pend.append(mir_items[nxt]); nxt += 1
        bops, brl = bit
        found = None
        for use_fp in (True, False):
          if found is not None:
            break
          for j, mit in enumerate(pend):
            mops, mrl = mit
            if len(mops) == len(bops):
                newmap, newmap_inv = {}, {}
                swc, okb = 0, True
                for mo, bo in zip(mops, bops):
                    r = fits_ops(mo, bo, newmap, use_fp)
                    if r is None:
                        okb = False
                        break
                    loc, linv, sw = r
                    newmap.update(loc); newmap_inv.update(linv)
                    swc += sw
                if okb and all(commute_items(pend[q][1], mrl) for q in range(j)):
                    found = (j, newmap, newmap_inv, swc)
                    used = use_fp
                    break
            if j >= window - 1:
                break
        if found is not None and not used:
            st['fp_overrides'] = st.get('fp_overrides', 0) + 1
        if found is None:
            return dict(ok=False, pi=None, msg='нет пары для элемента %d: %r' % (i, bops[:2]), stats=st, dbg=(pi, inv, list(pend), bit, fm, fb))
        j, newmap, newmap_inv, swc = found
        if j:
            st['reordered'] += 1
            st['max_window'] = max(st['max_window'], j)
        st['ab_swaps'] += swc
        pi.update(newmap)
        inv.update(newmap_inv)
        del pend[j]
    return dict(ok=True, pi=pi, msg='ok', stats=st)
