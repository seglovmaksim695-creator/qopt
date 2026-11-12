"""Экспертная итерация v2 — на быстром моторе qfast (годится для 6–10+ кубитов).

Отличия от expert_iter.py:
  - нет плотных матриц: память и время растут как 2^n, а не 4^n;
  - на каждом шаге поиска дёшево считается только амплитуда |0..0> для ВСЕХ ходов,
    полное новое состояние — только для лучших кандидатов (3 x ширина поиска);
  - опционально плотная подсказка gamma * (средняя вероятность кубитов быть в |0>),
    полезна для многих кубитов, где точная близость к цели почти ноль до самого конца;
  - ранняя остановка, когда экзамен достиг цели (--stop-at);
  - ограничение числа учебных примеров по памяти.
Файлы моделей/учебника те же, что у expert_iter.py (можно продолжать старые запуски).

Запуск:  caffeinate -i .venv/bin/python -u ei2.py --qubits 6 --rounds 40 --stop-at 0.9
"""

import argparse
import json
import random
import os
import time
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from qfast import ChoiEngine, Engine, state_sig


def make_engine(n, task=None):
    """task="unitary" — синтез операций (этап 10): операция n кубитов как состояние 2n кубитов."""
    return ChoiEngine(n) if task == "unitary" else Engine(n)
from train_agent import QNet

ROOT = Path(__file__).resolve().parent
SUCCESS = 0.99
_net = _eng = None
_cfg = {}
_val = None


def atomic_text(path, text):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def atomic_checkpoint(value, path):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        torch.save(value, stream)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def default_max_steps(n):
    return {2: 30, 3: 40, 4: 60, 5: 60}.get(n, 12 * n)


def obs_batch(states):
    k = np.argmax(np.abs(states), axis=1)
    ph = states[np.arange(len(states)), k]
    v = states * (np.abs(ph) / ph)[:, None]
    return np.concatenate([v.real, v.imag], axis=1).astype(np.float32)


def _init(state_dict, n, hidden, cfg):
    global _net, _eng, _cfg
    torch.set_num_threads(1)
    _eng = make_engine(n, cfg.get("task"))
    if cfg.get("arch") in ("gnn", "gnn2", "gnn3", "gnnu", "gnnh", "gnnhp"):  # сеть, не зависящая от числа кубитов (qgnn.py)
        import qgnn
        _net = {"gnn": qgnn.QGNN, "gnn2": qgnn.QGNN2, "gnn3": qgnn.QGNN3, "gnnu": qgnn.QGNNU, "gnnh": qgnn.QGNNH, "gnnhp": qgnn.QGNNHP}[cfg["arch"]](**cfg.get("gnn_kw", {}))
    else:
        _net = QNet(2 * _eng.dim, _eng.n_actions, hidden=hidden)
    _net.load_state_dict(state_dict); _net.eval()
    global _val
    _val = None
    if cfg.get("value"):  # оценщик перспективности (этап 10г)
        import qgnn
        vck = torch.load(cfg["value"], map_location="cpu")
        _val = qgnn.QValueU(**vck["kw"]); _val.load_state_dict(vck["state_dict"]); _val.eval()
    _cfg = cfg


def beam(target, width, bound=None):
    """bound — известное решение. Число: его длина (дальше искать незачем). Пара (cx, длина) — режим «по CX»:
    ветки, где CX уже больше известного, отбрасываются, а длина ограничена с запасом (длина × CXSLACK),
    чтобы не пропустить более длинную схему с меньшим числом CX."""
    cx_cap = None
    if isinstance(bound, (tuple, list)):
        cx_cap, bound = bound[0], int(np.ceil(bound[1] * _cfg.get("len_slack", 1.5)))
    eng, cfg = _eng, _cfg
    beta, gamma, max_steps = cfg["beta"], cfg["gamma"], cfg["max_steps"]
    inv = np.array(eng.inverse_action)
    states = np.asarray(target, complex)[None, :]
    acts, lps = [[]], np.zeros(1)
    best, stop = None, (min(max_steps, bound) if bound else max_steps)
    for depth in range(1, max_steps + 1):
        with torch.no_grad():
            logp = torch.log_softmax(_net(torch.from_numpy(obs_batch(states))), 1).numpy()
        fid = np.abs(eng.amp0_all(states)) ** 2
        score = lps[:, None] + logp
        last = np.array([a[-1] if a else -1 for a in acts])
        has = last >= 0
        score[np.flatnonzero(has), inv[last[has]]] = -np.inf
        rank = score + beta * fid
        ok = np.isfinite(rank)
        for i, a in np.argwhere(ok & (fid >= SUCCESS)):
            na = acts[i] + [int(a)]
            key = (int(eng.is_cx[na].sum()), len(na))
            if best is None or key < best[0]:
                best = (key, na); stop = min(stop, depth + cfg["extra"])
        rank = np.where(ok & (fid < SUCCESS), rank, -np.inf)
        if cx_cap is not None:  # ветки, где CX уже больше, чем в известном решении, — бесполезны
            cxn = np.array([int(eng.is_cx[a].sum()) if a else 0 for a in acts])[:, None] + eng.is_cx[None, :]
            rank = np.where(cxn > cx_cap, -np.inf, rank)
        k = min(3 * width, int(np.isfinite(rank).sum()))
        if depth >= stop or k == 0:
            break
        flat = np.argpartition(-rank, k - 1, axis=None)[:k]
        bi, ai = np.divmod(flat, eng.n_actions)
        ns = eng.apply(states[bi], ai)
        r2 = rank[bi, ai] + (gamma * eng.p0_per_qubit(ns) if gamma else 0.0)
        if _val is not None:  # оценщик: потрачено CX + ожидается ещё CX — чем меньше, тем перспективнее
            with torch.no_grad():
                v = _val(torch.from_numpy(obs_batch(ns))).numpy()
            spent = np.array([int(eng.is_cx[acts[b]].sum()) if acts[b] else 0 for b in bi]) + eng.is_cx[ai]
            r2 = r2 - cfg.get("lam", 1.0) * (spent + v)
        new_s, new_a, new_l, seen = [], [], [], set()
        for j in np.argsort(-r2):
            s = ns[j]
            sig = state_sig(s, 4)
            if sig in seen:
                continue
            seen.add(sig)
            new_s.append(s); new_a.append(acts[bi[j]] + [int(ai[j])]); new_l.append(score[bi[j], ai[j]])
            if len(new_s) >= width:
                break
        states, acts, lps = np.array(new_s), new_a, np.array(new_l)
    return best


def _solve(args):
    idx, target, width = args[:3]
    bound = args[3] if len(args) > 3 else None
    res = beam(target, width, bound)
    if res is not None and _cfg.get("optimize"):
        from circuit_opt import optimize
        acts = optimize(_eng, target, res[1])
        res = ((int(_eng.is_cx[acts].sum()), len(acts)), acts)
    return idx, res


def run_pool(state_dict, n, hidden, cfg, jobs, workers):
    with Pool(workers, initializer=_init, initargs=(state_dict, n, hidden, cfg)) as pool:
        return pool.map(_solve, jobs, chunksize=4)


def exact_ok(eng, target, acts):
    s = np.asarray(target, complex)
    for a in acts:
        s = eng.apply1(s, a)
    return abs(amp0_any(eng, s)) ** 2 / np.vdot(s, s).real >= 1 - 1e-6


def amp0_any(eng, s):
    return eng.amp0(s) if hasattr(eng, "amp0") else s[0]


def perm_maps(eng, perm):
    """Перестановка кубитов q -> perm[q]: индексы амплитуд и номера действий."""
    n = eng.n
    idx = np.arange(eng.dim)
    jmap = np.zeros(eng.dim, dtype=np.intp)
    full = list(perm) + ([n + int(p) for p in perm] if hasattr(eng, "n_total") else [])  # Г3: двойники вместе с рабочими
    for q in range(len(full)):
        jmap |= ((idx >> q) & 1) << full[q]
    amap = np.array([eng.gates.index((nm, tuple(int(perm[x]) for x in qb))) for nm, qb in eng.gates])
    return jmap, amap


def train_epochs(net, opt, X, Y, epochs, eng, perm_aug, np_rng):
    """Обучение подражанием; perm_aug>0 — к каждому пакету добавляются копии с переставленными кубитами."""
    Xt, Yt = torch.from_numpy(X), torch.from_numpy(Y)
    Mt = torch.from_numpy(FRONT_M) if (FRONT_M is not None and len(FRONT_M) == len(X)) else None
    D = eng.dim
    maps = [perm_maps(eng, np_rng.permutation(eng.n)) for _ in range(32)] if perm_aug else []
    loss = torch.tensor(0.0)
    net.train()
    for _ in range(epochs):
        order = torch.randperm(len(Xt))
        for b in range(0, len(Xt), 512):
            j = order[b:b + 512]
            xb, yb = Xt[j], Yt[j]
            mb = Mt[j] if Mt is not None else None
            if perm_aug:
                xs, ys, ms = [xb], [yb], [mb]
                for _ in range(perm_aug):
                    jmap, amap = maps[int(np_rng.integers(len(maps)))]
                    xp = torch.empty_like(xb)
                    jm = torch.from_numpy(jmap)
                    xp[:, jm] = xb[:, :D]
                    xp[:, D + jm] = xb[:, D:]
                    # заново приводим глобальную фазу, как obs_batch (самая большая амплитуда — вещественная)
                    z = torch.complex(xp[:, :D], xp[:, D:])
                    k = torch.argmax(z.abs(), dim=1)
                    ph = z[torch.arange(len(z)), k]
                    z = z * (ph.abs() / ph)[:, None]
                    xp = torch.cat([z.real, z.imag], 1)
                    xs.append(xp); ys.append(torch.from_numpy(amap)[yb])
                    if mb is not None:
                        mp = torch.zeros_like(mb); mp[:, torch.from_numpy(amap)] = mb; ms.append(mp)
                xb, yb = torch.cat(xs), torch.cat(ys)
                if mb is not None:
                    mb = torch.cat(ms)
            lg = net(xb)
            if mb is not None:  # метки-множества: −log Σ_{a∈фронт} p(a)
                loss = (torch.logsumexp(lg, 1) - torch.logsumexp(lg.masked_fill(~mb, -1e9), 1)).mean()
            else:
                loss = nn.functional.cross_entropy(lg, yb)
            opt.zero_grad(); loss.backward(); opt.step()
    net.eval()
    return loss


FRONT_M = None  # маски фронта для последнего вызова trajectories (FRONT=1)


def front_masks(eng, acts):
    from block_synth import commutes
    g = [(eng.gates[a][0], tuple(eng.gates[a][1])) for a in acts]
    out = []
    for i in range(len(acts)):
        m = np.zeros(eng.n_actions, bool)
        for j in range(i, len(acts)):
            if all(commutes(g[j], g[t]) for t in range(i, j)):
                m[acts[j]] = True
        out.append(m)
    return out


def trajectories(eng, targets, solutions, cap, rng):
    global FRONT_M
    keys = list(solutions)
    rng.shuffle(keys)
    X, Y, M, total = [], [], [], 0
    use_front = bool(os.environ.get("FRONT"))
    for k in keys:
        s = np.asarray(targets[int(k)], complex)
        if use_front:
            M += front_masks(eng, solutions[k][1])
        for a in solutions[k][1]:
            X.append(obs_batch(s[None])[0]); Y.append(a)
            s = eng.apply1(s, a)
        total += len(solutions[k][1])
        if total >= cap:
            break
    FRONT_M = np.array(M, bool) if use_front else None
    return np.array(X, np.float32), np.array(Y, np.int64)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--qubits", type=int, required=True)
    p.add_argument("--rounds", type=int, default=40)
    p.add_argument("--per-round", type=int, default=4000)
    p.add_argument("--width", type=int, default=128)
    p.add_argument("--eval-n", type=int, default=400)
    p.add_argument("--epochs", type=int, default=6)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--hidden", type=int, default=1024)
    p.add_argument("--beta", type=float, default=3.0)
    p.add_argument("--gamma", type=float, default=0.0)
    p.add_argument("--max-steps", type=int, default=0)
    # --- для экспериментов (по умолчанию поведение прежнее) ---
    p.add_argument("--run-name", default="", help="отдельные файлы в exp/ (эксперимент не трогает основные)")
    p.add_argument("--init-from", default="", help="стартовая модель для нового эксперимента")
    p.add_argument("--init-solutions", default="", help="стартовый учебник для нового эксперимента")
    p.add_argument("--optimize", action="store_true", help="чистить решения оптимизатором circuit_opt")
    p.add_argument("--curriculum", default="easy", choices=["easy", "zone"],
                   help="easy: лёгкие нерешённые; zone: цели с ожидаемым успехом ~30-70%%")
    p.add_argument("--time-budget", type=float, default=0, help="минут работы (0 = без ограничения)")
    p.add_argument("--final-eval", action="store_true", help="в конце — экзамен последней модели на всём тесте")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--pretrain-epochs", type=int, default=0, help="сначала обучить сеть на стартовом учебнике")
    p.add_argument("--perm-aug", type=int, default=0, help="копий с переставленными кубитами на пакет")
    p.add_argument("--warmup", type=int, default=0, help="первые N раундов — только самые лёгкие цели")
    p.add_argument("--final-width", type=int, default=128)
    p.add_argument("--widen", type=int, default=0,
                   help="расширение поиска: неудачи повторить с этой шириной (обучение и экзамен)")
    p.add_argument("--value", default="", help="оценщик перспективности (models/value_u<n>q.pt)")
    p.add_argument("--lam", type=float, default=1.0, help="вес оценщика")
    p.add_argument("--data-tag", default="", help="другие данные: data/<u>train<tag>_<n>q.npz и test (этап 11: _blocks)")
    p.add_argument("--bound", nargs="?", const="len", default=None, choices=["len", "cx"],
                   help="при обучении не искать решения длиннее уже известного (ускорение поиска)")
    p.add_argument("--task", default="state", choices=["state", "unitary"],
                   help="unitary — синтез операций (данные data/utrain_*, data/utest_*)")
    p.add_argument("--replay", type=int, nargs="*", default=[],
                   help="смешанное обучение (только для GNN): подмешивать учебники этих размеров, напр. 4 5 6 7")
    p.add_argument("--replay-n", type=int, default=20000, help="примеров каждого старого размера за раунд")
    p.add_argument("--gen-fill", default="", help="учебник из генераторов: брать решение ТОЛЬКО для целей, "
                                                 "которые поиск не решил")
    p.add_argument("--stop-at", type=float, default=1.01, help="остановиться, когда экзамен >= этого 2 раунда подряд")
    a = p.parse_args()
    n = a.qubits
    torch.set_num_threads(4)
    eng = make_engine(n, a.task)
    obs_dim = 2 * eng.dim
    PFX = "u" if a.task == "unitary" else ""
    cap = int(min(600_000, 1.5e9 / (obs_dim * 4)))  # не больше ~1,5 ГБ учебных примеров в памяти

    if a.run_name:
        out = ROOT / "exp"; out.mkdir(exist_ok=True)
        latest = out / f"{a.run_name}_{n}q.pt"
        sol_file = out / f"{a.run_name}_solutions_{n}q.json"
        log_file = out / f"{a.run_name}_log_{n}q.jsonl"
        best_file = out / f"{a.run_name}_best_{n}q.pt"
        journal_file = out / f"{a.run_name}.log"
    else:
        latest = ROOT / "models" / f"ei_{n}q.pt"
        sol_file = ROOT / "data" / f"ei_solutions_{n}q.json"
        log_file = ROOT / "results" / f"ei_log_{n}q.jsonl"
        best_file = ROOT / "models" / f"ei_best_{n}q.pt"
        journal_file = ROOT / "JOURNAL.md"
    stats_file = sol_file.with_name(sol_file.name.replace("solutions", "stats"))
    if a.init_from and not latest.exists():
        ck0 = torch.load(a.init_from)
        torch.save({**ck0, "round": 0, "best_test": 0.0}, latest)
    # стартовый учебник — НЕЗАВИСИМО от стартовой модели (ошибка до 4 нояб.: без --init-from он не грузился)
    if a.init_solutions and not sol_file.exists():
        import shutil
        shutil.copy(a.init_solutions, sol_file)
    if latest.exists():
        ck = torch.load(latest)
        hidden, max_steps = ck.get("hidden", 512), a.max_steps or ck.get("max_steps", default_max_steps(n))
        start, best_test = ck["round"] + 1, ck["best_test"]
        best_cx = ck.get("best_cx", float("inf"))
    else:
        ck, hidden, max_steps = None, a.hidden, a.max_steps or default_max_steps(n)
        start, best_test, best_cx = 1, 0.0, float("inf")
    if a.run_name:
        torch.manual_seed(a.seed)
    ARCH = {"arch": ck["arch"], "gnn_kw": ck["gnn_kw"]} if (ck and ck.get("arch") in ("gnn", "gnn2", "gnn3", "gnnu", "gnnh", "gnnhp")) else {}
    if ARCH:  # сеть, не зависящая от числа кубитов (этап 6)
        import qgnn
        net = {"gnn": qgnn.QGNN, "gnn2": qgnn.QGNN2, "gnn3": qgnn.QGNN3, "gnnu": qgnn.QGNNU, "gnnh": qgnn.QGNNH, "gnnhp": qgnn.QGNNHP}[ARCH["arch"]](**ARCH["gnn_kw"])
    else:
        net = QNet(obs_dim, eng.n_actions, hidden=hidden)
    if ck:
        net.load_state_dict(ck["state_dict"])
    cfg = {"beta": a.beta, "gamma": a.gamma, "max_steps": max_steps, "extra": 3, "optimize": a.optimize, **ARCH,
           "task": a.task, "data_tag": a.data_tag, **({"value": a.value, "lam": a.lam} if a.value else {})}

    solutions = {}
    if sol_file.exists():
        solutions = {k: (tuple(v[0]), v[1]) for k, v in json.loads(sol_file.read_text()).items()}
    tr = np.load(ROOT / "data" / f"{PFX}train{a.data_tag}_{n}q.npz"); T = tr["targets"]; DIFF = tr["qiskit_cx"]
    te = np.load(ROOT / "data" / f"{PFX}test{a.data_tag}_{n}q.npz"); TT = te["targets"]
    sp = te["statprep_cx"] if "statprep_cx" in te.files else None
    eval_idx = np.random.default_rng(7).choice(len(TT), size=min(a.eval_n, len(TT)), replace=False)
    opt = torch.optim.Adam(net.parameters(), lr=3e-4)
    if ck and "opt" in ck:
        opt.load_state_dict(ck["opt"])
    rng = random.Random(a.seed)
    np_rng = np.random.default_rng(a.seed)
    stats = json.loads(stats_file.read_text()) if stats_file.exists() else {}  # цель -> [попыток, успехов]
    t_start = time.time()
    loss = torch.tensor(0.0)
    hits = 0
    print(f"старт: в учебнике {len(solutions)} решений; предобучение {a.pretrain_epochs} эпох; "
          f"перестановки {a.perm_aug}; разогрев {a.warmup}; оптимизатор {a.optimize}", flush=True)
    if a.init_solutions and not solutions:
        raise SystemExit("ОШИБКА: --init-solutions задан, но учебник пуст — останавливаюсь")
    if a.pretrain_epochs and not solutions:
        raise SystemExit("ОШИБКА: --pretrain-epochs задан, но учиться не на чем — останавливаюсь")
    if a.pretrain_epochs and solutions:
        X, Y = trajectories(eng, T, solutions, cap, rng)
        loss = train_epochs(net, opt, X, Y, a.pretrain_epochs, eng, a.perm_aug, np_rng)
        print(f"предобучение на учебнике: {len(solutions)} решений, {len(X)} примеров, "
              f"{(time.time() - t_start) / 60:.1f} мин", flush=True)
        del X, Y
    easy = np.argsort(DIFF, kind="stable")[:max(1, len(T) // 10)]  # 10% самых лёгких целей
    REPLAY = {}
    if a.replay:
        assert ARCH, "--replay имеет смысл только для сети, не зависящей от числа кубитов"
        for m in a.replay:
            em = Engine(m)
            Tm = np.load(ROOT / "data" / f"train_{m}q.npz")["targets"]
            sm = {k: (tuple(v[0]), v[1]) for k, v in json.loads((ROOT / "data" / f"ei_solutions_{m}q.json").read_text()).items()}
            Xm, Ym = trajectories(em, Tm, sm, 5 * a.replay_n, rng)
            REPLAY[m] = (em, Xm, Ym)
        print("смешанное обучение: " + ", ".join(f"{m}q {len(v[1])} прим." for m, v in REPLAY.items()), flush=True)
    GEN = {}
    if a.gen_fill:
        GEN = {k: (tuple(v[0]), v[1]) for k, v in json.loads(Path(a.gen_fill).read_text()).items()}
        print(f"учебник для нерешённых: {len(GEN)} решений в запасе", flush=True)

    for rnd in range(start, start + a.rounds):
        t0 = time.time()
        # учебный план: из нерешённых берём случайно среди самых лёгких (по CX у Qiskit)
        if a.curriculum == "zone":
            pick = zone_pick(stats, DIFF, solutions, a.per_round // 2, np_rng)
        else:
            # учебный план: из нерешённых берём случайно среди самых лёгких (по CX у Qiskit)
            unsolved = sorted((i for i in range(len(T)) if str(i) not in solutions), key=lambda i: DIFF[i])
            k = min(len(unsolved), a.per_round // 2)
            pick = rng.sample(unsolved[:2 * k], k)
        if a.warmup and rnd - start < a.warmup:
            pick = [int(x) for x in np_rng.choice(easy, size=min(a.per_round, len(easy)), replace=False)]
        else:
            pick += rng.sample(range(len(T)), a.per_round - len(pick))
        sd = {k: v.cpu() for k, v in net.state_dict().items()}
        found = improved = 0
        jobs = [(i, T[i], a.width, (((solutions[str(i)][0][0], len(solutions[str(i)][1])) if a.bound == "cx"
                                      else len(solutions[str(i)][1])) if (a.bound and str(i) in solutions) else None))
                for i in pick]
        results = run_pool(sd, n, hidden, cfg, jobs, a.workers)
        if a.widen:  # расширение: неудачи — ещё раз, шире
            failed = {r[0] for r in results if r[1] is None}
            if failed:
                again = dict(run_pool(sd, n, hidden, cfg, [(j[0], j[1], a.widen, *j[3:]) for j in jobs
                                                           if j[0] in failed], a.workers))
                results = [(i, r if r is not None else again.get(i)) for i, r in results]
        for idx, res in results:
            st = stats.setdefault(str(idx), [0, 0])
            st[0] += 1
            if res is None:
                if str(idx) in GEN and str(idx) not in solutions:
                    solutions[str(idx)] = GEN[str(idx)]  # поиск не справился — берём решение генератора
                continue
            st[1] += 1
            found += 1
            old = solutions.get(str(idx))
            if not exact_ok(eng, T[idx], res[1]):  # Г1 Н3: принимать в учебник только точные решения
                continue
            if old is None or tuple(res[0]) < tuple(old[0]):
                solutions[str(idx)] = (tuple(res[0]), res[1]); improved += 1
        atomic_text(sol_file, json.dumps({k: [list(v[0]), v[1]] for k, v in solutions.items()}))
        atomic_text(stats_file, json.dumps(stats))
        t_search = time.time() - t0

        if solutions:
            X, Y = trajectories(eng, T, solutions, cap, rng)
            loss = train_epochs(net, opt, X, Y, a.epochs, eng, a.perm_aug, np_rng)
            del X, Y
            for m, (em, Xm, Ym) in REPLAY.items():  # смешанное обучение: повторяем старые размеры
                sel = np_rng.choice(len(Xm), size=min(a.replay_n, len(Xm)), replace=False)
                train_epochs(net, opt, Xm[sel], Ym[sel], 1, em, 0, np_rng)

        sd = {k: v.cpu() for k, v in net.state_dict().items()}
        res = run_pool(sd, n, hidden, cfg, [(int(i), TT[i], 64) for i in eval_idx], a.workers)
        if a.widen:
            fail = [int(i) for i, r in res if r is None]
            if fail:
                again = dict(run_pool(sd, n, hidden, cfg, [(i, TT[i], a.widen) for i in fail], a.workers))
                res = [(i, r if r is not None else again.get(int(i))) for i, r in res]
        ok = [(i, r) for i, r in res if r is not None]
        solved = len(ok) / len(eval_idx)
        cx = float(np.mean([r[0][0] for _, r in ok])) if ok else 0.0
        sp_cx = float(np.mean([sp[i] for i, _ in ok if sp[i] >= 0] or [np.nan])) if (ok and sp is not None) else float("nan")
        rec = {"round": rnd, "time": datetime.now().isoformat(timespec="minutes"), "engine": "qfast",
               "train_solved_total": len(solutions), "train_total": len(T), "found": found,
               "improved": improved, "loss": float(loss.detach()), "test_solved": solved, "test_cx": cx,
               "statprep_cx_same": sp_cx, "minutes": (time.time() - t0) / 60,
               "search_minutes": t_search / 60, "gamma": a.gamma, "width": a.width}
        with open(log_file, "a") as f:
            f.write(json.dumps(rec) + "\n")
        # лучшая модель: больше решено, а при равном проценте — меньше CX (до 7 нояб. учитывался только процент,
        # и при 100% best переставал обновляться — замечание автора)
        better = solved > best_test or (solved == best_test and np.isfinite(cx) and cx < best_cx - 1e-9)
        if better:
            best_cx = cx
        atomic_checkpoint({"state_dict": net.state_dict(), "round": rnd, "best_test": max(best_test, solved),
                    "best_cx": best_cx, "hidden": hidden, "max_steps": max_steps, **ARCH,
                    "opt": opt.state_dict()}, latest)  # состояние Adam — продолжение без сброса (П5)
        mark = ""
        if better:
            best_test = solved
            atomic_checkpoint({"state_dict": net.state_dict(), "round": rnd, "test_solved": solved,
                        "hidden": hidden, "max_steps": max_steps, **ARCH}, best_file)
            mark = " ★ новый рекорд"
        line = (f"- [авто] {rec['time']} {n}q раунд {rnd}: экзамен {solved:.1%}, CX {cx:.2f} vs {'оптимизатор' if a.data_tag else 'Qiskit'} {sp_cx:.2f}; "
                f"в учебнике {len(solutions)}/{len(T)}; {rec['minutes']:.1f} мин (поиск {t_search / 60:.1f}){mark}")
        print(line, flush=True)
        with open(journal_file, "a") as f:
            f.write(line + "\n")
        if a.time_budget and (time.time() - t_start) / 60 >= a.time_budget:
            print(f"Время ({a.time_budget:.0f} мин) вышло — останавливаюсь.", flush=True)
            break
        hits = hits + 1 if solved >= a.stop_at else 0
        if hits >= 2:
            print(f"Цель экзамена {a.stop_at:.0%} достигнута — останавливаюсь.", flush=True)
            break
    if a.final_eval:
        sd = {k: v.cpu() for k, v in net.state_dict().items()}
        fe = final_eval(n, hidden, cfg, sd, a.workers, a.final_width)
        fe.update({"run": a.run_name, "minutes": (time.time() - t_start) / 60,
                   "train_solved_total": len(solutions), "optimize": a.optimize, "curriculum": a.curriculum,
                   "seed": a.seed,
                   "textbook_cx": float(np.mean([v[0][0] for v in solutions.values()])) if solutions else None,
                   "textbook_len": float(np.mean([v[0][1] for v in solutions.values()])) if solutions else None})
        (ROOT / "exp").mkdir(exist_ok=True)
        with open(ROOT / "exp" / "final_results.jsonl", "a") as f:
            f.write(json.dumps(fe) + "\n")
        print("ИТОГ:", json.dumps(fe, ensure_ascii=False), flush=True)


def zone_pick(stats, diff, solutions, k, np_rng):
    """Подборщик: цели с ожидаемым успехом около 50% (зона 30-70%).
    Для целей, которые уже пробовали, — своя оценка (успехи+1)/(попытки+2);
    для новых — средний успех по их уровню трудности (CX у Qiskit)."""
    lv_att, lv_ok = {}, {}
    for key, (att, ok) in stats.items():
        lv = int(diff[int(key)])
        lv_att[lv] = lv_att.get(lv, 0) + att
        lv_ok[lv] = lv_ok.get(lv, 0) + ok
    p = np.empty(len(diff))
    for i in range(len(diff)):
        st = stats.get(str(i))
        if st:
            p[i] = (st[1] + 1) / (st[0] + 2)
        else:
            lv = int(diff[i])
            p[i] = (lv_ok.get(lv, 0) + 1) / (lv_att.get(lv, 0) + 2)
    w = np.exp(-((p - 0.5) / 0.2) ** 2)
    sure = [int(x) for x in solutions if p[int(x)] > 0.9]
    w[sure] *= 0.1  # уже уверенно решённые — реже
    w /= w.sum()
    return [int(x) for x in np_rng.choice(len(diff), size=min(k, len(diff)), replace=False, p=w)]


def final_eval(n, hidden, cfg, state_dict, workers, width=128):
    te = np.load(ROOT / "data" / f"{'u' if cfg.get('task') == 'unitary' else ''}test{cfg.get('data_tag', '')}_{n}q.npz")
    TT = te["targets"]
    res = dict(run_pool(state_dict, n, hidden, {**cfg, "optimize": False},
                        [(i, TT[i], width) for i in range(len(TT))], workers))
    ok = [i for i in range(len(TT)) if res[i] is not None]
    cx = [res[i][0][0] for i in ok]
    out = {"tests": len(TT), "solved": len(ok) / len(TT), "agent_cx": float(np.mean(cx)) if cx else None}
    if "statprep_cx" in te.files and ok:
        out["statprep_cx_same"] = float(np.mean([c for c in te["statprep_cx"][ok] if c >= 0] or [np.nan]))
    return out


if __name__ == "__main__":
    main()
