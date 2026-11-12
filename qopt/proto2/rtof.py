"""Прототип: оптимизация на уровне Тоффоли до разложения в Clifford+T.
Идея: пары Tof(a,b;t) ... Tof(a,b;t), между которыми a,b,t встречаются только «диагонально»
(управляющие CX/Tof, Z/CZ/CCZ/фазы), заменяются на Tof с относительной фазой: R=D·Tof (3 CX) и R^† (3 CX).
Итог: D^† M D = M, точное равенство всей операции. Затем стандартное разложение + global_cancel проекта.
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[1])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, types, time
from pathlib import Path

# ---- разбор .qc (Feynman) -------------------------------------------------
def load_qc(path):
    names, gates = None, []
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == ".v":
            names = p[1:]
            idx = {nm: i for i, nm in enumerate(names)}
            continue
        if p[0].startswith(".") or p[0] in ("BEGIN", "END"):
            continue
        op, qs = p[0], [idx[q] for q in p[1:]]
        if op == "H": gates.append(("h", (qs[0],)))
        elif op == "X": gates.append(("x", (qs[0],)))
        elif op in ("Z", "Zd"):
            gates.append({1: "z", 2: "cz", 3: "ccz"}[len(qs)] and ({1: "z", 2: "cz", 3: "ccz"}[len(qs)], tuple(qs)))
        elif op == "tof":
            gates.append({1: "x", 2: "cx", 3: "ccx"}[len(qs)] and ({1: "x", 2: "cx", 3: "ccx"}[len(qs)], tuple(qs)))
        else:
            raise ValueError(op)
    return len(names), gates


# ---- H·CCZ·H -> CCX, H·CZ·H -> CX ---------------------------------------------
def _wire_nb(g, i, t, d):
    k = i + d
    while 0 <= k < len(g):
        if g[k] and t in g[k][1]:
            return k
        k += d
    return None


def cancel_h(g):
    """H_t H_t, соседние на проводе t (между ними нет команд на t) -> убрать. Точно."""
    g = list(g)
    last = {}
    for i, x in enumerate(g):
        if x is None:
            continue
        nm, qs = x
        if nm == "h":
            t = qs[0]
            if t in last and g[last[t]] == ("h", (t,)):
                g[last[t]] = g[i] = None
                last.pop(t)
                continue
        for q in qs:
            last[q] = i
    return [x for x in g if x]


def fold_h(g):
    """CCZ(a,b,t) = H_t CCX(a,b;t) H_t. Цель t выбираем там, где по соседству на проводе стоят H
    (в т.ч. «виртуально»: H·CCZ·M·CCZ·H у Баренко), затем сокращаем соседние H на проводах."""
    g = list(g); out = []
    for i in range(len(g)):
        nm, qs = g[i]
        if nm not in ("ccz", "cz"):
            out.append(g[i]); continue
        best, bs = None, (0, 0)
        for t in qs:
            sc = 0; dist = 0
            for d in (-1, 1):
                k = _wire_nb(g, i, t, d)
                while k is not None and g[k][0] in ("ccz", "cz") and t in g[k][1] and set(g[k][1]) == set(qs):
                    k = _wire_nb(g, k, t, d)  # пропустить такой же CCZ (он тоже станет H·CCX·H)
                if k is not None and g[k] == ("h", (t,)):
                    sc += 1; dist -= abs(k - i)
            sc = (sc, dist)
            if sc > bs:
                best, bs = t, sc
        if best is not None and bs[0] == 2:
            ctrl = tuple(q for q in qs if q != best)
            out += [("h", (best,)), ("ccx" if nm == "ccz" else "cx", ctrl + (best,)), ("h", (best,))]
            CHOICE[i] = best
        else:
            out.append(g[i])
    return cancel_h(out)


CHOICE = {}


def diag_on(gate, q):
    """Команда действует на кубит q диагонально (в вычислительном базисе)?"""
    nm, qs = gate
    if q not in qs:
        return True
    if nm in ("z", "cz", "ccz", "t", "tdg", "s", "sdg"):
        return True
    if nm in ("cx", "ccx"):
        return q != qs[-1]
    return False  # h, x, цель cx/ccx


def pair_toffolis(g, use4=True):
    """Жадно: для каждого CCX ищем следующий CCX с теми же управляющими, такой что между ними
    команды на опоре S фазы D диагональны. Варианты:
      'R'  — Маргулус (3 CX), D на {a,b,t}: нужна та же цель и диагональность на {a,b,t};
      'R4' — RTOF4 (4 CX), D = CS^†(a,b) только на управляющих: цель любая, диагональность только на {a,b}.
    Возвращает {индекс: (вид, сопряжённый?)}."""
    mark = {}
    for i, (nm, qs) in enumerate(g):
        if nm != "ccx" or i in mark:
            continue
        C, S = set(qs[:2]), set(qs)
        ok3 = True
        for j in range(i + 1, len(g)):
            h = g[j]
            if h[0] == "ccx" and set(h[1][:2]) == C and j not in mark:
                if ok3 and h[1][-1] == qs[-1]:
                    mark[i], mark[j] = ("R", 0, qs), ("R", 1, qs); break
                if use4:
                    mark[i], mark[j] = ("R4", 0, qs), ("R4", 1, qs[:2] + h[1][2:]); break
            ok3 = ok3 and all(diag_on(h, q) for q in S)
            if not all(diag_on(h, q) for q in C) or (not use4 and not ok3):
                break
    return mark


# ---- разложение в Clifford+T ----------------------------------------------------
def ccx_std(a, b, t):
    return [("h", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)), ("cx", (b, t)),
            ("tdg", (t,)), ("cx", (a, t)), ("t", (b,)), ("t", (t,)), ("h", (t,)), ("cx", (a, b)),
            ("t", (a,)), ("tdg", (b,)), ("cx", (a, b))]


def rccx(a, b, t):  # = D·Tof, D диагональна (проверено)
    return [("h", (t,)), ("t", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)),
            ("cx", (b, t)), ("tdg", (t,)), ("h", (t,))]


def rccx4(a, b, t):  # = CS^†(a,b)·Tof (проверено), D не зависит от цели
    return [("h", (t,)), ("t", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)),
            ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("h", (t,))]


INV = {"h": "h", "x": "x", "cx": "cx", "t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s", "z": "z"}


def inv(c):
    return [(INV[n], q) for n, q in reversed(c)]


def lower(g, mark=None):
    mark = mark or {}
    out = []
    for i, (nm, qs) in enumerate(g):
        if nm == "ccx":
            a, b, t = qs
            kind, dg, qq = mark.get(i, (None, 0, qs))
            a, b, t = qq
            if kind is None:
                out += ccx_std(a, b, t)
            else:
                c = rccx(a, b, t) if kind == "R" else rccx4(a, b, t)
                out += inv(c) if dg else c
        elif nm == "ccz":
            a, b, t = qs
            out += [("h", (t,))] + ccx_std(a, b, t) + [("h", (t,))]
        elif nm == "cz":
            out += [("h", (qs[1],)), ("cx", qs), ("h", (qs[1],))]
        elif nm == "z":
            out += [("s", qs), ("s", qs)]
        else:
            out.append((nm, qs))
    return out


# ---- global_cancel проекта (без torch) ---------------------------------------
def _load_gc():
    # если настоящий block_synth уже загружен (bench_real) — не подменять его заглушкой (ошибка 8 нояб.)
    if "block_synth" in sys.modules and hasattr(sys.modules["block_synth"], "W"):
        sys.path.insert(0, _QPKG)
        import global_cancel
        return global_cancel.global_cancel
    src = Path((_QPKG + "/block_synth.py")).read_text()
    ns = {}
    start = src.index("DIAG = "); end = src.index("def take_block_commute")
    exec(src[start:end], ns)
    m = types.ModuleType("block_synth"); m.commutes, m.DIAG = ns["commutes"], ns["DIAG"]
    sys.modules["block_synth"] = m
    sys.path.insert(0, _QPKG)
    import global_cancel
    return global_cancel.global_cancel


global_cancel = _load_gc()


def cx(c):
    return sum(n == "cx" for n, _ in c)


# ---- выбор целей CCZ->H·CCX·H локальным поиском по итоговому числу CX -----------
def expand(g, choice):
    out = []
    for i, x in enumerate(g):
        t = choice.get(i)
        if t is None:
            out.append(x)
        else:
            nm, qs = x
            ctrl = tuple(q for q in qs if q != t)
            out += [("h", (t,)), ("ccx" if nm == "ccz" else "cx", ctrl + (t,)), ("h", (t,))]
    return cancel_h(out)


def score(g2, use4=True):
    c = global_cancel(lower(g2, pair_toffolis(g2, use4)))
    return (cx(c), len(c))


def optimize_targets(g, rounds=3, use4=True):
    idx = [i for i, x in enumerate(g) if x[0] in ("ccz", "cz")]
    CHOICE.clear(); fold_h(g); choice = dict(CHOICE)
    best = score(expand(g, choice), use4)
    for _ in range(rounds):
        improved = False
        for i in idx:
            cur = choice.get(i)
            for t in [None] + list(g[i][1]):
                if t == cur:
                    continue
                choice[i] = t
                s = score(expand(g, choice), use4)
                if s < best:
                    best, cur, improved = s, t, True
            choice[i] = cur
        if not improved:
            break
    return expand(g, choice), best
