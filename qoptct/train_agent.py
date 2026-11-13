"""Обучение ИИ-агента собирать квантовые программы (алгоритм Double DQN).

Запуск:
    .venv/bin/python train_agent.py                 # 2 кубита, по умолчанию
    .venv/bin/python train_agent.py --qubits 3 --steps 600000

Обучение идёт ступенями (учебный план): сначала простые цели (1–2 команды),
когда агент решает >90% — цели усложняются.
Результат: models/agent_<N>q.pt и график results/training_<N>q.png
"""

import argparse
import json
import random
import time
from collections import deque
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from qenv import StatePrepEnv

ROOT = Path(__file__).resolve().parent


class QNet(nn.Module):
    def __init__(self, obs_dim, n_actions, hidden=256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(obs_dim, hidden), nn.ReLU(),
            nn.Linear(hidden, hidden), nn.ReLU(),
            nn.Linear(hidden, hidden), nn.ReLU(),
            nn.Linear(hidden, n_actions))

    def forward(self, x):
        return self.net(x)


def choose_action(qnet, obs, env, epsilon):
    banned = env.inverse_action[env.last_action()] if env.last_action() is not None else None
    if random.random() < epsilon:
        choices = [a for a in range(env.n_actions) if a != banned]
        return random.choice(choices)
    with torch.no_grad():
        q = qnet(torch.from_numpy(obs).unsqueeze(0))[0]
    if banned is not None:
        q[banned] = -1e9  # бессмысленно сразу отменять предыдущую команду
    return int(q.argmax())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--qubits", type=int, default=2)
    p.add_argument("--steps", type=int, default=300_000)
    p.add_argument("--max-depth", type=int, default=None, help="предел сложности целей")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    n = args.qubits
    depth_limit = args.max_depth or (8 if n == 2 else 10)
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    torch.set_num_threads(4)

    env = StatePrepEnv(n_qubits=n, max_steps=4 * depth_limit, max_depth=2, seed=args.seed)
    qnet = QNet(env.obs_dim, env.n_actions)
    target_net = QNet(env.obs_dim, env.n_actions)
    target_net.load_state_dict(qnet.state_dict())
    opt = torch.optim.Adam(qnet.parameters(), lr=5e-4)
    gamma, batch_size = 0.97, 256

    buf_size = 200_000
    B_obs = np.zeros((buf_size, env.obs_dim), np.float32)
    B_next = np.zeros((buf_size, env.obs_dim), np.float32)
    B_act = np.zeros(buf_size, np.int64)
    B_rew = np.zeros(buf_size, np.float32)
    B_done = np.zeros(buf_size, np.float32)
    ptr, size = 0, 0

    recent = deque(maxlen=500)
    log = []
    obs = env.reset()
    t0 = time.time()
    episodes = 0

    for step in range(1, args.steps + 1):
        epsilon = max(0.05, 1.0 - step / (0.4 * args.steps))
        a = choose_action(qnet, obs, env, epsilon)
        nxt, r, done, info = env.step(a)
        B_obs[ptr], B_act[ptr], B_rew[ptr], B_next[ptr] = obs, a, r, nxt
        B_done[ptr] = float(info["success"])  # обрыв по лимиту шагов — не конец в смысле ценности
        ptr = (ptr + 1) % buf_size
        size = min(size + 1, buf_size)
        obs = nxt

        if done:
            episodes += 1
            recent.append(info["success"])
            obs = env.reset()
            # учебный план: усложняем цели, когда агент справляется
            if len(recent) == recent.maxlen and np.mean(recent) > 0.9 and env.max_depth < depth_limit:
                env.max_depth += 1
                recent.clear()
                print(f"  ↑ шаг {step}: сложность целей повышена до {env.max_depth} команд", flush=True)

        if size >= 5_000 and step % 2 == 0:
            idx = np.random.randint(0, size, batch_size)
            o = torch.from_numpy(B_obs[idx]); nx = torch.from_numpy(B_next[idx])
            act = torch.from_numpy(B_act[idx]); rew = torch.from_numpy(B_rew[idx])
            dn = torch.from_numpy(B_done[idx])
            with torch.no_grad():
                best = qnet(nx).argmax(1, keepdim=True)
                target_q = rew + gamma * (1 - dn) * target_net(nx).gather(1, best).squeeze(1)
            q = qnet(o).gather(1, act.unsqueeze(1)).squeeze(1)
            loss = nn.functional.smooth_l1_loss(q, target_q)
            opt.zero_grad(); loss.backward()
            nn.utils.clip_grad_norm_(qnet.parameters(), 10)
            opt.step()

        if step % 2_000 == 0:
            target_net.load_state_dict(qnet.state_dict())

        if step % 10_000 == 0:
            sr = float(np.mean(recent)) if recent else 0.0
            log.append({"step": step, "success_rate": sr, "depth": env.max_depth, "epsilon": epsilon})
            print(f"шаг {step:>7} | успех {sr:5.1%} | сложность {env.max_depth} | "
                  f"эпизодов {episodes} | {time.time() - t0:5.0f} c", flush=True)

    (ROOT / "models").mkdir(exist_ok=True)
    (ROOT / "results").mkdir(exist_ok=True)
    model_path = ROOT / "models" / f"agent_{n}q.pt"
    torch.save({"state_dict": qnet.state_dict(), "qubits": n, "obs_dim": env.obs_dim,
                "n_actions": env.n_actions, "depth_reached": env.max_depth,
                "depth_limit": depth_limit}, model_path)
    (ROOT / "results" / f"training_{n}q.json").write_text(json.dumps(log, indent=1))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots(figsize=(9, 4.5))
    ax1.plot([e["step"] for e in log], [e["success_rate"] for e in log], label="доля решённых задач")
    ax1.set_xlabel("шаг обучения"); ax1.set_ylabel("успех"); ax1.set_ylim(0, 1.05)
    ax2 = ax1.twinx()
    ax2.step([e["step"] for e in log], [e["depth"] for e in log], color="tab:orange", label="сложность целей")
    ax2.set_ylabel("сложность (команд)")
    fig.legend(loc="lower right"); ax1.set_title(f"Обучение агента: {n} кубита")
    fig.tight_layout(); fig.savefig(ROOT / "results" / f"training_{n}q.png", dpi=120)
    print(f"Готово. Модель: {model_path}")


if __name__ == "__main__":
    main()
