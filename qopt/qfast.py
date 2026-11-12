"""Быстрый «мотор» для любого числа кубитов (без плотных матриц).

Каждая наша команда (H, X, S, S†, T, T†, CX) меняет каждую амплитуду, глядя
не более чем на 2 старые амплитуды. Поэтому новое состояние считается как
    new[i] = A[a, i] * s[P1[a, i]] + B[a, i] * s[P2[a, i]]
Это O(2^n) на команду вместо O(4^n), и не нужны матрицы 2^n x 2^n.

Порядок действий совпадает с qenv.build_gate_set (проверяется в selftest()).
"""

import numpy as np

SQ2 = 1 / np.sqrt(2)
PHASE = {"s": 1j, "sdg": -1j, "t": np.exp(1j * np.pi / 4), "tdg": np.exp(-1j * np.pi / 4)}
NAMES_1Q = ["h", "x", "s", "sdg", "t", "tdg"]
INVERSE = {"h": "h", "x": "x", "s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t", "cx": "cx"}


class Engine:
    def __init__(self, n):
        self.n = n
        self.dim = D = 2 ** n
        idx = np.arange(D)
        gates, P1, P2, A, B = [], [], [], [], []
        for q in range(n):
            m = 1 << q
            bit = (idx >> q) & 1
            for name in NAMES_1Q:
                gates.append((name, (q,)))
                if name == "h":
                    P1.append(idx & ~m); P2.append(idx | m)
                    A.append(np.full(D, SQ2, complex)); B.append(np.where(bit, -SQ2, SQ2).astype(complex))
                elif name == "x":
                    P1.append(idx ^ m); P2.append(idx)
                    A.append(np.ones(D, complex)); B.append(np.zeros(D, complex))
                else:
                    P1.append(idx); P2.append(idx)
                    A.append(np.where(bit, PHASE[name], 1).astype(complex)); B.append(np.zeros(D, complex))
        for c in range(n):
            for t in range(n):
                if c == t:
                    continue
                gates.append(("cx", (c, t)))
                P1.append(np.where((idx >> c) & 1, idx ^ (1 << t), idx)); P2.append(idx)
                A.append(np.ones(D, complex)); B.append(np.zeros(D, complex))
        self.gates = gates
        self.n_actions = len(gates)
        self.P1, self.P2 = np.array(P1), np.array(P2)
        self.A, self.B = np.array(A), np.array(B)
        self.is_cx = np.array([g[0] == "cx" for g in gates])
        self.inverse_action = [gates.index((INVERSE[nm], qb)) for nm, qb in gates]
        # для быстрого расчёта только амплитуды |0...0> после каждого хода
        self.p1_0, self.p2_0 = self.P1[:, 0], self.P2[:, 0]
        self.a_0, self.b_0 = self.A[:, 0], self.B[:, 0]

    def amp0_all(self, states):
        """(B, D) -> (B, actions): амплитуда |0..0> после каждого возможного хода. Очень дёшево."""
        return states[:, self.p1_0] * self.a_0 + states[:, self.p2_0] * self.b_0

    def apply(self, states, actions):
        """states (K, D), actions (K,) -> новые состояния (K, D)."""
        r = np.arange(len(actions))[:, None]
        a = np.asarray(actions)
        return (self.A[a] * states[r, self.P1[a]] + self.B[a] * states[r, self.P2[a]])

    def apply1(self, state, action):
        return self.A[action] * state[self.P1[action]] + self.B[action] * state[self.P2[action]]

    def p0_per_qubit(self, states):
        """Для каждого состояния: средняя по кубитам вероятность, что кубит в |0>. Плотный сигнал близости."""
        prob = np.abs(states) ** 2
        idx = np.arange(self.dim)
        out = np.zeros(len(states))
        for q in range(self.n):
            out += prob[:, ((idx >> q) & 1) == 0].sum(1)
        return out / self.n


def selftest(n=3):
    from qenv import build_gate_set
    eng = Engine(n)
    dense = build_gate_set(n)
    assert [(g[0], g[1]) for g in dense] == eng.gates, "порядок действий не совпадает"
    rng = np.random.default_rng(0)
    s = rng.normal(size=(4, 2 ** n)) + 1j * rng.normal(size=(4, 2 ** n))
    err = 0.0
    for a, (_, _, m) in enumerate(dense):
        err = max(err, np.abs(eng.apply(s, [a] * 4) - s @ m.T).max())
    err = max(err, np.abs(eng.amp0_all(s) - np.stack([s @ m[0] for _, _, m in dense], 1)).max())
    return err


if __name__ == "__main__":
    for n in (2, 3, 4):
        print(n, "кубита: максимальная ошибка", selftest(n))


def amp0(eng, s):
    """Амплитуда «цели» у состояния s (для обычной задачи — амплитуда |0...0>)."""
    return eng.amp0(s) if hasattr(eng, "amp0") else s[0]



def state_sig(t, dec=6):
    """Устойчивая подпись состояния до глобальной фазы: опорная амплитуда — первая из максимальных (после округления),
    целые числа вместо float (нет -0.0). Отчёт критика: прежняя подпись промахивалась в 264/300 одинаковых операций."""
    t = np.asarray(t, complex)
    a = np.round(np.abs(t), 5)
    m0 = int(np.argmax(a))
    v = t * (abs(t[m0]) / t[m0])
    return np.rint(np.stack([v.real, v.imag]) * 10 ** dec).astype(np.int64).tobytes()


class ChoiEngine:
    """Синтез ОПЕРАЦИЙ (этап 10). Операция U на n кубитах записывается как состояние на 2n кубитах:
    |U> = (U ⊗ I)|Φ>, где |Φ> — n пар «рабочий кубит q — двойник n+q» в состоянии Белла (|00>+|11>)/√2.
    Команды разрешены ТОЛЬКО на рабочих кубитах 0..n-1, порядок и номера действий — как у Engine(n).
    Цель «размотки» — снова |Φ> (операция стала тождественной). Совпадение = |<Φ|ψ>|^2 = |tr(U)/2^n|^2.
    Интерфейс как у Engine: apply, apply1, amp0_all, is_cx, inverse_action, gates, n (=n рабочих), dim (=4^n)."""

    def __init__(self, n):
        full = Engine(2 * n)
        small = Engine(n)
        self.n, self.n_total = n, 2 * n
        self.dim = full.dim
        self.gates = small.gates
        self.n_actions = small.n_actions
        self.is_cx, self.inverse_action = small.is_cx, small.inverse_action
        sel = [full.gates.index(g) for g in small.gates]  # те же команды, но в пространстве 2n кубитов
        self.P1, self.P2, self.A, self.B = full.P1[sel], full.P2[sel], full.A[sel], full.B[sel]
        g = np.ones(1, complex)
        for q in range(n):  # пара (q, n+q) в состоянии Белла
            pass
        idx = np.arange(self.dim)
        ok = np.ones(self.dim, bool)
        for q in range(n):
            ok &= ((idx >> q) & 1) == ((idx >> (n + q)) & 1)
        self.goal = ok.astype(complex) / np.sqrt(ok.sum())
        # W[a] = g_a^† |Φ>  ->  <Φ| g_a ψ> = <W[a] | ψ>
        self.Wc = np.stack([self.apply1(self.goal, self.inverse_action[a]) for a in range(self.n_actions)]).conj()

    def apply(self, states, actions):
        r = np.arange(len(actions))[:, None]
        a = np.asarray(actions)
        return self.A[a] * states[r, self.P1[a]] + self.B[a] * states[r, self.P2[a]]

    def apply1(self, state, action):
        return self.A[action] * state[self.P1[action]] + self.B[action] * state[self.P2[action]]

    def amp0_all(self, states):
        return states @ self.Wc.T

    def amp0(self, s):
        return np.vdot(self.goal, s)

    def choi(self, program):
        """Состояние |U> для программы (список действий), U = g_last ... g_first."""
        s = self.goal.copy()
        for a in program:
            s = self.apply1(s, a)
        return s


def choi_selftest(n=2):
    """Проверка: совпадение |tr(U)/2^n|^2 и размотка случайной программы обратно в цель."""
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    eng = ChoiEngine(n)
    rng = np.random.default_rng(0)
    prog = [int(x) for x in rng.integers(0, eng.n_actions, 12)]
    s = eng.choi(prog)
    qc = QuantumCircuit(n)
    for a in prog:
        nm, qb = eng.gates[a]
        getattr(qc, nm)(*qb)
    U = Operator(qc).data
    f1 = abs(eng.amp0(s)) ** 2
    f2 = abs(np.trace(U) / 2 ** n) ** 2
    back = s
    for a in reversed(prog):
        back = eng.apply1(back, eng.inverse_action[a])
    return f1, f2, abs(eng.amp0(back)) ** 2
