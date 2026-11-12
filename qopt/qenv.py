"""Тренажёр для ИИ: подготовка заданного состояния кубитов.

Как устроено:
  - Есть цель: квантовое состояние (например, запутанное состояние 2 кубитов).
  - Агент по одной выбирает команды (гейты) из набора H, X, S, S†, T, T†, CX.
  - Хитрость для удобства обучения: агент "разматывает" целевое состояние
    обратно в |00...0>. Найденную цепочку потом переворачиваем — получается
    программа, которая ИЗ |00...0> строит цель.
  - Награда: за приближение к цели, штраф за каждую лишнюю команду,
    бонус за успех. Так агент учится решать точно и коротко.

Порядок кубитов — как в Qiskit (кубит 0 — младший бит).
"""

import numpy as np

SQ2 = 1 / np.sqrt(2)
I2 = np.eye(2, dtype=complex)
ONE_QUBIT = {
    "h": np.array([[1, 1], [1, -1]], dtype=complex) * SQ2,
    "x": np.array([[0, 1], [1, 0]], dtype=complex),
    "s": np.diag([1, 1j]).astype(complex),
    "sdg": np.diag([1, -1j]).astype(complex),
    "t": np.diag([1, np.exp(1j * np.pi / 4)]).astype(complex),
    "tdg": np.diag([1, np.exp(-1j * np.pi / 4)]).astype(complex),
}
INVERSE = {"h": "h", "x": "x", "s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t", "cx": "cx"}

SUCCESS_FIDELITY = 0.99


def one_qubit_matrix(u, q, n):
    m = np.array([[1]], dtype=complex)
    for i in range(n - 1, -1, -1):
        m = np.kron(m, u if i == q else I2)
    return m


def cx_matrix(c, t, n):
    dim = 2 ** n
    m = np.zeros((dim, dim), dtype=complex)
    for b in range(dim):
        out = b ^ (1 << t) if (b >> c) & 1 else b
        m[out, b] = 1
    return m


def build_gate_set(n):
    """Список действий: (имя, кубиты, матрица)."""
    gates = []
    for q in range(n):
        for name, u in ONE_QUBIT.items():
            gates.append((name, (q,), one_qubit_matrix(u, q, n)))
    for c in range(n):
        for t in range(n):
            if c != t:
                gates.append(("cx", (c, t), cx_matrix(c, t, n)))
    return gates


class StatePrepEnv:
    def __init__(self, n_qubits=2, max_steps=20, max_depth=4, seed=None,
                 gate_cost=0.03, cx_cost=0.03):
        self.gate_cost = gate_cost  # штраф за обычную команду
        self.cx_cost = cx_cost      # штраф за CX (на железе она самая "шумная")
        self.n = n_qubits
        self.dim = 2 ** n_qubits
        self.max_steps = max_steps
        self.max_depth = max_depth
        self.rng = np.random.default_rng(seed)
        self.gates = build_gate_set(n_qubits)
        self.n_actions = len(self.gates)
        self.obs_dim = 2 * self.dim
        # индекс обратного действия для каждого действия
        self.inverse_action = []
        for name, qubits, _ in self.gates:
            inv = INVERSE[name]
            self.inverse_action.append(
                next(i for i, g in enumerate(self.gates) if g[0] == inv and g[1] == qubits))

    # ---------- цели ----------
    def random_target(self, depth):
        """Случайная цель = результат случайной программы длины depth."""
        while True:
            state = np.zeros(self.dim, dtype=complex)
            state[0] = 1
            program = []
            for _ in range(depth):
                a = int(self.rng.integers(self.n_actions))
                state = self.gates[a][2] @ state
                program.append(a)
            if abs(state[0]) ** 2 < SUCCESS_FIDELITY:  # цель не должна совпадать с |0>
                return state, program

    # ---------- эпизод ----------
    def reset(self, target=None, depth=None):
        if target is None:
            depth = depth or int(self.rng.integers(1, self.max_depth + 1))
            target, self.reference_program = self.random_target(depth)
        else:
            self.reference_program = None
        self.target = np.asarray(target, dtype=complex)
        self.state = self.target.copy()
        self.applied = []
        self.steps = 0
        return self.observation()

    def observation(self):
        v = self.state
        k = int(np.argmax(np.abs(v)))
        v = v * (abs(v[k]) / v[k])  # убираем глобальную фазу — она физически не важна
        return np.concatenate([v.real, v.imag]).astype(np.float32)

    def fidelity(self):
        return float(abs(self.state[0]) ** 2)

    def step(self, action):
        before = self.fidelity()
        self.state = self.gates[action][2] @ self.state
        self.applied.append(action)
        self.steps += 1
        after = self.fidelity()
        success = after >= SUCCESS_FIDELITY
        done = success or self.steps >= self.max_steps
        cost = self.cx_cost if self.gates[action][0] == "cx" else self.gate_cost
        reward = (after - before) - cost + (1.0 if success else 0.0)
        return self.observation(), reward, done, {"success": success, "fidelity": after}

    def last_action(self):
        return self.applied[-1] if self.applied else None

    # ---------- результат ----------
    def solution(self):
        """Программа, которая из |0...0> строит цель: [(имя, кубиты), ...]."""
        out = []
        for a in reversed(self.applied):
            name, qubits, _ = self.gates[a]
            out.append((INVERSE[name], qubits))
        return out


def to_qiskit(program, n):
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(n)
    for name, qubits in program:
        getattr(qc, name)(*qubits)
    return qc
