"""Экспертная итерация: поиск находит решения -> агент учится на лучших -> поиск становится сильнее.

Подробное описание метода: METHOD.md. Журнал: JOURNAL.md.

Запуск (продолжает с последнего сохранения, если оно есть):
    caffeinate -i .venv/bin/python -u expert_iter.py --rounds 30
"""

import argparse
import json
import random
import time
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from qenv import INVERSE, SUCCESS_FIDELITY, StatePrepEnv
from train_agent import QNet

ROOT = Path(__file__).resolve().parent
HIDDEN = 512
_net = _env = _mats = _iscx = None
_fast_gates = None


def _init(state_dict, n, max_steps, hidden=512, fast_gates=False):
    global _net, _env, _mats, _iscx, _fast_gates
    torch.set_num_threads(1)
    _env = StatePrepEnv(n_qubits=n, max_steps=max_steps)
    _net = QNet(_env.obs_dim, _env.n_actions, hidden=hidden)
    _net.load_state_dict(state_dict); _net.eval()
    _mats = np.stack([g[2] for g in _env.gates])
    _iscx = np.array([g[0] == "cx" for g in _env.gates])
    _fast_gates = None
    if fast_gates:
        from fast_gates import FastGateBatch
        _fast_gates = FastGateBatch(_mats)


def _obs(states):
    k = np.argmax(np.abs(states), axis=1)
    ph = states[np.arange(len(states)), k]
    v = states * (np.abs(ph) / ph)[:, None]
    return np.concatenate([v.real, v.imag], axis=1).astype(np.float32)


def beam(target, width, beta=3.0, extra=3):
    """Поиск лучом. Ранг варианта = сумма log-вероятностей ходов сети + beta * близость к цели."""
    env = _env
    inv = env.inverse_action
    states = np.asarray(target, complex)[None, :]
    acts = [[]]; lps = np.zeros(1)
    best, stop = None, env.max_steps
    for depth in range(1, env.max_steps + 1):
        with torch.no_grad():
            logp = torch.log_softmax(_net(torch.from_numpy(_obs(states))), 1).numpy()
        ns = (_fast_gates(states) if _fast_gates is not None
              else np.einsum("aij,bj->bai", _mats, states))  # все ходы из всех вариантов
        fid = np.abs(ns[:, :, 0]) ** 2
        score = lps[:, None] + logp
        for i, a in enumerate(acts):
            if a:
                score[i, inv[a[-1]]] = -np.inf
        rank = score + beta * fid
        order = np.argsort(-rank, axis=None)
        new_s, new_a, new_l, seen = [], [], [], set()
        for flat in order:
            i, a = divmod(int(flat), env.n_actions)
            if not np.isfinite(rank[i, a]):
                break
            na = acts[i] + [a]
            if fid[i, a] >= SUCCESS_FIDELITY:
                key = (int(_iscx[na].sum()), len(na))
                if best is None or key < best[0]:
                    best = (key, na); stop = min(stop, depth + extra)
                continue
            s = ns[i, a]
            k = int(np.argmax(np.abs(s)))
            sig = tuple(np.round(s * (abs(s[k]) / s[k]), 5).view(float))
            if sig in seen:
                continue
            seen.add(sig)
            new_s.append(s); new_a.append(na); new_l.append(score[i, a])
            if len(new_s) >= width:
                break
        if depth >= stop or not new_s:
            break
        states, acts, lps = np.array(new_s), new_a, np.array(new_l)
    return best


def _solve(args):
    idx, target, width = args
    return idx, beam(target, width)


def run_pool(state_dict, n, max_steps, jobs, workers, hidden=512, fast_gates=False):
    with Pool(workers, initializer=_init, initargs=(state_dict, n, max_steps, hidden, fast_gates)) as pool:
        return pool.map(_solve, jobs, chunksize=8)


def trajectories(env, targets, solutions):
    """Из найденных решений делаем учебные примеры (состояние -> правильный ход)."""
    X, Y = [], []
    for idx, (_, actions) in solutions.items():
        env.reset(target=targets[int(idx)])
        for a in actions:
            X.append(env.observation()); Y.append(a)
            env.step(a)
    return np.array(X, np.float32), np.array(Y, np.int64)


def journal(text):
    with open(ROOT / "JOURNAL.md", "a") as f:
        f.write(text + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--qubits", type=int, default=3)
    p.add_argument("--rounds", type=int, default=30)
    p.add_argument("--per-round", type=int, default=4000, help="сколько учебных целей решать за раунд")
    p.add_argument("--width", type=int, default=96, help="ширина поиска при сборе решений")
    p.add_argument("--eval-n", type=int, default=400)
    p.add_argument("--epochs", type=int, default=4)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--hidden", type=int, default=512, help="ширина сети (если нет сохранения)")
    p.add_argument("--fast-gates", action="store_true", help="ускоренный расчёт разреженных гейтов")
    a = p.parse_args()
    n = a.qubits
    max_steps = {2: 30, 3: 40}.get(n, 60)
    torch.set_num_threads(4)

    latest = ROOT / "models" / f"ei_{n}q.pt"
    global HIDDEN
    HIDDEN = torch.load(latest).get("hidden", a.hidden) if latest.exists() else a.hidden
    env = StatePrepEnv(n_qubits=n, max_steps=max_steps)
    net = QNet(env.obs_dim, env.n_actions, hidden=HIDDEN)
    sol_file = ROOT / "data" / f"ei_solutions_{n}q.json"
    log_file = ROOT / "results" / f"ei_log_{n}q.jsonl"
    if latest.exists():
        ck = torch.load(latest); net.load_state_dict(ck["state_dict"])
        start, best_test = ck["round"] + 1, ck["best_test"]
    else:
        init = ROOT / "models" / f"hard_{n}q.pt"
        if init.exists():
            net.load_state_dict(torch.load(init)["state_dict"])
        # иначе — с нуля: поначалу поиск держится на близости к цели, сеть учится по найденному
        start, best_test = 1, 0.0
    solutions = {}
    if sol_file.exists():
        solutions = {k: (tuple(v[0]), v[1]) for k, v in json.loads(sol_file.read_text()).items()}

    tr = np.load(ROOT / "data" / f"train_{n}q.npz"); T, L = tr["targets"], tr["qiskit_cx"]
    te = np.load(ROOT / "data" / f"test_{n}q.npz")
    rng_eval = np.random.default_rng(7)
    eval_idx = rng_eval.choice(len(te["targets"]), size=min(a.eval_n, len(te["targets"])), replace=False)
    opt = torch.optim.Adam(net.parameters(), lr=3e-4)

    loss = torch.tensor(0.0)
    for rnd in range(start, start + a.rounds):
        t0 = time.time()
        # 1) собираем решения: половина — ещё не решённые цели, половина — случайные (ищем короче)
        unsolved = [i for i in range(len(T)) if str(i) not in solutions]
        pick = random.sample(unsolved, min(len(unsolved), a.per_round // 2))
        pick += random.sample(range(len(T)), a.per_round - len(pick))
        sd = {k: v.cpu() for k, v in net.state_dict().items()}
        found = improved = 0
        for idx, res in run_pool(sd, n, max_steps, [(i, T[i], a.width) for i in pick], a.workers, HIDDEN, a.fast_gates):
            if res is None:
                continue
            found += 1
            old = solutions.get(str(idx))
            if old is None or tuple(res[0]) < tuple(old[0]):
                solutions[str(idx)] = (tuple(res[0]), res[1]); improved += 1
        sol_file.write_text(json.dumps({k: [list(v[0]), v[1]] for k, v in solutions.items()}))

        # 2) учимся на лучших найденных решениях (подражание "эксперту" = поиску)
        X, Y = trajectories(env, T, solutions)
        Xt, Yt = torch.from_numpy(X), torch.from_numpy(Y)
        net.train()
        for _ in range(a.epochs):
            perm = torch.randperm(len(Xt))
            for b in range(0, len(Xt), 512):
                j = perm[b:b + 512]
                loss = nn.functional.cross_entropy(net(Xt[j]), Yt[j])
                opt.zero_grad(); loss.backward(); opt.step()
        net.eval()

        # 3) экзамен на целях, которых агент никогда не видел
        sd = {k: v.cpu() for k, v in net.state_dict().items()}
        res = run_pool(sd, n, max_steps, [(int(i), te["targets"][i], 64) for i in eval_idx], a.workers, HIDDEN, a.fast_gates)
        ok = [(i, r) for i, r in res if r is not None]
        solved = len(ok) / len(eval_idx)
        cx = float(np.mean([r[0][0] for _, r in ok])) if ok else 0.0
        at_min = sum(r[0][0] == n - 1 for _, r in ok)
        sp_cx = float(np.mean([te["statprep_cx"][i] for i, _ in ok])) if ok else 0.0
        hard = [i for i in eval_idx if te["qiskit_cx"][i] >= 9]
        ok_set = {i for i, _ in ok}
        hard_rate = sum(i in ok_set for i in hard) / max(1, len(hard))

        rec = {"round": rnd, "time": datetime.now().isoformat(timespec="minutes"),
               "train_solved_total": len(solutions), "train_total": len(T), "found": found,
               "improved": improved, "examples": len(X), "loss": float(loss.detach()),
               "test_solved": solved, "test_hard_solved": hard_rate, "test_cx": cx,
               "test_at_min": at_min, "statprep_cx_same": sp_cx, "minutes": (time.time() - t0) / 60}
        with open(log_file, "a") as f:
            f.write(json.dumps(rec) + "\n")
        (ROOT / "models").mkdir(exist_ok=True)
        torch.save({"state_dict": net.state_dict(), "round": rnd, "best_test": max(best_test, solved),
                    "hidden": HIDDEN, "max_steps": max_steps}, latest)
        mark = ""
        if solved > best_test:
            best_test = solved
            torch.save({"state_dict": net.state_dict(), "round": rnd, "test_solved": solved,
                        "hidden": HIDDEN, "max_steps": max_steps}, ROOT / "models" / f"ei_best_{n}q.pt")
            mark = " ★ новый рекорд"
        line = (f"- [авто] {rec['time']} раунд {rnd}: экзамен {solved:.1%} (трудные {hard_rate:.0%}), "
                f"CX {cx:.2f} vs Qiskit {sp_cx:.2f}; в учебнике {len(solutions)}/{len(T)} целей; "
                f"{rec['minutes']:.1f} мин{mark}")
        print(line, flush=True)
        journal(line)


if __name__ == "__main__":
    main()
