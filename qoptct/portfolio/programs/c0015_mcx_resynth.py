"""mcx_resynth: распознать в уровне Тоффоли многоуправляемый X с грязными вспомогательными и синтезировать заново.
1) Цели CCZ — по «кадру H» (fold_h + нечётность H на проводе); после expand схема классическая (только ccx/cx/x),
   её действие вычисляется точно в АНФ (многочлены над GF(2)). Если провод t: t ^= c1·…·ck (k ≥ 3), а все прочие
   провода тождественны, то это C^kX, и остальные n-k-1 ≥ k-2 провода — грязные вспомогательные.
2) Синтез Барэнко как A·(M~ A M~^{-1}): A = CCX(c_k, d_{k-3} -> t) — точный (6 CX), M~ — V-цепочка d_0 ^= c1·…·c_{k-1}
   (2k-5 CCX) целиком из Тоффоли с относительной фазой (rccx, 3 CX): M~ = D·M, D диагональна и не задевает t,
   A диагональна вне t => M~ A M~^{-1} = M A M точно; мусор на d_1.. откатывается M~^{-1}.
3) Ориентация каждого rccx (порядок управляющих = «управляющий-вспомогательный первым») и вид двух A подобраны так,
   что соседние rccx цепочки делят CX после fold_cancel: итог 8k-8 CX (k=6: 40, k=7: 48, k=10: 72 против 12k-18).
   Поверх — короткий координатный спуск по вариантам. Иначе (не C^kX) — стандартное разложение (портфель берёт min)."""
NAME = "mcx_resynth"; KIND = "lower"; HL_RAW = True
from evo_lib import (fold_cancel, key, inv_circ, rccx, ccx_std, lower_std, fold_h_choice, expand)


def frame_choice(g):
    par = {}; ch = {}
    for i, (nm, qs) in enumerate(g):
        if nm == "h":
            par[qs[0]] = par.get(qs[0], 0) ^ 1
        elif nm in ("ccz", "cz"):
            odd = [q for q in qs if par.get(q, 0)]
            if len(odd) == 1:
                ch[i] = odd[0]
    return ch


CAP = 4096


def pmul(P, Q):
    r = set()
    for m in P:
        for k in Q:
            r ^= {m | k}
    return r


def anf(g, n):
    """Многочлены проводов (множества одночленов-битовых масок) или None, если схема не классическая/слишком велика."""
    poly = [{1 << q} for q in range(n)]
    for nm, qs in g:
        if nm == "x":
            poly[qs[0]] = poly[qs[0]] ^ {0}
        elif nm == "cx":
            poly[qs[1]] = poly[qs[1]] ^ poly[qs[0]]
        elif nm == "ccx":
            poly[qs[2]] = poly[qs[2]] ^ pmul(poly[qs[0]], poly[qs[1]])
            if len(poly[qs[2]]) > CAP:
                return None
        else:
            return None
    return poly


def detect(hl, n):
    g0 = [(x[0], tuple(x[1])) for x in hl]
    ch = fold_h_choice(g0); ch.update(frame_choice(g0))
    for choice in (ch, fold_h_choice(g0)):
        g = expand(g0, choice)
        poly = anf(g, n)
        if poly is None:
            continue
        moved = [q for q in range(n) if poly[q] != {1 << q}]
        if len(moved) != 1:
            continue
        t = moved[0]; rest = poly[t] ^ {1 << t}
        if len(rest) != 1:
            continue
        m = next(iter(rest))
        if m & (1 << t):
            continue
        C = [q for q in range(n) if m >> q & 1]
        if len(C) < 3:
            continue
        D = [q for q in range(n) if q != t and q not in C]
        if len(D) < len(C) - 2:
            continue
        return C, t, D
    return None


def rv(qs, v):
    a, b, t = qs
    if v & 1: a, b = b, a
    c = rccx(a, b, t)
    return inv_circ(c) if v & 2 else c


def sv(qs, v):
    a, b, t = qs
    if v & 1: a, b = b, a
    c = ccx_std(a, b, t)
    return inv_circ(c) if v & 2 else c


def chain(C, t, D):
    k = len(C); A = (C[k - 1], D[k - 3], t)
    down = [(C[i + 1], D[i - 1], D[i]) for i in range(k - 3, 0, -1)]
    return A, down + [(C[0], C[1], D[0])] + down[::-1]


def build(C, t, D, va, vm, order, n):
    A, M = chain(C, t, D)
    m = []
    for j, q in enumerate(M):
        m += rv(q, vm[j])
    a1 = sv(A, va[0]); a2 = sv(A, va[1])
    c = a1 + m + a2 + inv_circ(m) if order == 0 else m + a2 + inv_circ(m) + a1
    return fold_cancel(c, n)


def synth(C, t, D, n, max_evals=400):
    A, M = chain(C, t, D); L = len(M)
    starts = [(1, [2, 0], [1] * L), (0, [0, 2], [1] * L), (1, [2, 0], [3] * L), (0, [0, 0], [0] * L)]
    best = None
    for o, va, vm in starts:
        s = key(build(C, t, D, va, vm, o, n))
        if best is None or s < best[0]:
            best = (s, o, list(va), list(vm))
    b, o, va, vm = best; cnt = 0
    for _ in range(2):
        imp = False
        for arr in (vm, va):
            for j in range(len(arr)):
                cur = arr[j]
                for v in range(4):
                    if v == cur or cnt >= max_evals: continue
                    arr[j] = v; cnt += 1; s = key(build(C, t, D, va, vm, o, n))
                    if s < b: b, cur, imp = s, v, True
                arr[j] = cur
        if not imp or cnt >= max_evals: break
    return build(C, t, D, va, vm, o, n)


def run(hl, n):
    base = fold_cancel(lower_std([(x[0], tuple(x[1])) for x in hl]), n)
    det = None
    try:
        det = detect(hl, n)
    except Exception:
        det = None
    if det is None:
        return base
    C, t, D = det
    best = None
    # две раскладки управляющих/вспомогательных (по номерам и обратная) — порядок влияет на fold_cancel слабо
    for Cs, Ds in ((sorted(C), sorted(D)), (sorted(C)[::-1], sorted(D)[::-1])):
        c = synth(Cs, t, Ds[:len(Cs) - 2], n, max_evals=150 if len(C) > 9 else 300)
        if best is None or key(c) < key(best):
            best = c
    return best if key(best) < key(base) else base
