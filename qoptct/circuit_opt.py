"""Оптимизатор небольших схем: чистит найденные решения перед тем, как они попадут в учебник.

Решение хранится как список ходов, которые «разматывают» цель обратно в |0...0>.
Каждое изменение проверяется ТОЧНО: прогоняем цель через новый список ходов на быстром
моторе qfast и требуем, чтобы вероятность |0...0> осталась >= 0.99. Цель улучшения:
меньше CX, при равенстве — короче.

Проходы (повторяются, пока что-то улучшается; 2б — удаление пары «CX + любая команда»):
  1. удаление одной команды (убирает «мусор»: отменяющие друг друга пары после сдвига,
     фазовые команды на кубитах, где они ничего не меняют, и т.п.);
  2. удаление двух соседних команд;
  3. замена окна из 2–3 команд одной командой (например, T,T -> S).
Проверка: .venv/bin/python circuit_opt.py --qubits 6   (на учебнике, без записи)
"""

import argparse
import json
import time

import numpy as np

SUCCESS = 0.99


def _ok(eng, state, actions):
    s = state
    for a in actions:
        s = eng.apply1(s, a)
    from qfast import amp0
    import numpy as _np
    s = _np.asarray(s, complex)
    return abs(amp0(eng, s)) ** 2 / _np.vdot(s, s).real >= 1 - 1e-6  # точно (Г1 Н3), а не >= 0,99


def _key(eng, actions):
    return (int(eng.is_cx[actions].sum()) if actions else 0, len(actions))


def optimize(eng, target, actions, max_passes=4, time_limit=2.0):
    """Возвращает улучшенный список ходов (или исходный). Всегда проверено точно."""
    target = np.asarray(target, complex)
    best = list(actions)
    if not _ok(eng, target, best):
        return best
    t0 = time.time()
    n_act = eng.n_actions
    for _ in range(max_passes):
        improved = False
        # 1) и 2): удаление одной или двух соседних команд
        for width in (1, 2):
            i = 0
            while i + width <= len(best):
                cand = best[:i] + best[i + width:]
                if _ok(eng, target, cand):
                    best = cand; improved = True
                else:
                    i += 1
                if time.time() - t0 > time_limit:
                    return best
        # 2б) удаление любой пары команд, где первая — CX (например, две CX, разделённые
        #     командами, которые с ними коммутируют) — главный источник экономии CX
        cx_pos = [i for i, a in enumerate(best) if eng.is_cx[a]]
        found = False
        for i in cx_pos:
            for j in range(len(best)):
                if j == i:
                    continue
                cand = [a for k, a in enumerate(best) if k not in (i, j)]
                if _ok(eng, target, cand):
                    best = cand; improved = found = True
                    break
            if found or time.time() - t0 > time_limit:
                break
        # 3) замена окна из 2–3 команд одной командой (с проверкой, что стало лучше)
        for width in (2, 3):
            i = 0
            while i + width <= len(best):
                window_key = _key(eng, best[i:i + width])
                pre = target
                for a in best[:i]:
                    pre = eng.apply1(pre, a)
                done = False
                for r in range(n_act):
                    if (int(eng.is_cx[r]), 1) >= window_key:
                        continue
                    cand = best[:i] + [r] + best[i + width:]
                    if _ok(eng, pre, [r] + best[i + width:]):
                        best = cand; improved = done = True
                        break
                if not done:
                    i += 1
                if time.time() - t0 > time_limit:
                    return best
        if not improved:
            break
    return best


def main():
    from qfast import Engine
    p = argparse.ArgumentParser()
    p.add_argument("--qubits", type=int, default=6)
    p.add_argument("--solutions", default=None)
    p.add_argument("--limit", type=int, default=300)
    a = p.parse_args()
    n = a.qubits
    eng = Engine(n)
    T = np.load(f"data/train_{n}q.npz")["targets"]
    sol = json.loads(open(a.solutions or f"data/ei_solutions_{n}q.json").read())
    keys = list(sol)[:a.limit]
    cx0 = l0 = cx1 = l1 = better = 0
    t = time.time()
    for k in keys:
        acts = sol[k][1]
        new = optimize(eng, T[int(k)], acts)
        assert _ok(eng, np.asarray(T[int(k)], complex), new)
        c0, c1 = _key(eng, acts), _key(eng, new)
        cx0 += c0[0]; l0 += c0[1]; cx1 += c1[0]; l1 += c1[1]
        better += c1 < c0
    m = len(keys)
    print(f"{n} кубитов, {m} решений учебника: CX {cx0 / m:.2f} -> {cx1 / m:.2f}, "
          f"длина {l0 / m:.1f} -> {l1 / m:.1f}; улучшено {better} ({better / m:.0%}); "
          f"{(time.time() - t) / m * 1000:.0f} мс на решение")


if __name__ == "__main__":
    main()
