"""Прототип 2: свёртка фаз (phase folding, Amy/Nam rotation merging) по всей схеме в формализме суммы по путям.
Каждый провод несёт аффинную чётность путевых переменных (битовая маска + константа); H даёт новую переменную.
Фазы (T,S,Z, …) на одинаковой чётности складываются в первое вхождение, остальные удаляются. Точно (до общей фазы).
Затем global_cancel проекта сокращает CX, ставшие соседними."""
from rtof import global_cancel, cx

PH = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}
SEQ = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["s", "s"], 5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}


def phase_fold(c, n):
    wire = [(1 << q, 0) for q in range(n)]  # (маска, константа)
    nv = n
    first = {}   # маска -> индекс выхода
    acc = {}     # маска -> суммарная фаза (по модулю 8), в «положительной» полярности
    out = []
    for nm, qs in c:
        if nm in PH:
            q = qs[0]; m, k = wire[q]
            v = PH[nm] * (-1 if k else 1) % 8  # фаза на ¬p = −фаза на p (+ общая фаза)
            if m in first:
                acc[m] = (acc[m] + v) % 8
            else:
                first[m] = len(out); acc[m] = v
                out.append(("PH", m))
            continue
        if nm == "h":
            wire[qs[0]] = (1 << nv, 0); nv += 1
        elif nm == "x":
            m, k = wire[qs[0]]; wire[qs[0]] = (m, k ^ 1)
        elif nm == "cx":
            a, b = qs
            wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        out.append((nm, qs))
    # второй проход: развернуть PH в команды на проводе, где эта чётность лежит в момент first[m]
    wire = [(1 << q, 0) for q in range(n)]; nv = n
    res = []
    for nm, qs in out:
        if nm == "PH":
            m = qs
            q = next(i for i in range(n) if wire[i][0] == m)
            v = acc[m] * (-1 if wire[q][1] else 1) % 8
            res += [(s, (q,)) for s in SEQ[v]]
            continue
        if nm == "h":
            wire[qs[0]] = (1 << nv, 0); nv += 1
        elif nm == "x":
            m, k = wire[qs[0]]; wire[qs[0]] = (m, k ^ 1)
        elif nm == "cx":
            a, b = qs
            wire[b] = (wire[b][0] ^ wire[a][0], wire[b][1] ^ wire[a][1])
        res.append((nm, qs))
    return res


def tcount(c):
    return sum(nm in ("t", "tdg") for nm, _ in c)


def fold_cancel(c, n, iters=4):
    best = global_cancel(c)
    for _ in range(iters):
        c2 = global_cancel(phase_fold(best, n))
        if (cx(c2), tcount(c2), len(c2)) >= (cx(best), tcount(best), len(best)):
            break
        best = c2
    return best
