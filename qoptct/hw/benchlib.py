"""Benchmark circuits for ``qoptct hw`` (algorithmic circuits on 5-30 qubits, no measurements).

``suite(kind)`` returns a list of :class:`Bench` (name, family, circuit, prep).  ``prep`` is an optional layer of
one-qubit gates that prepares the input state used in the noisy simulation (a random basis state for reversible
arithmetic, the Fourier state of an integer for the QFT, nothing for circuits that prepare their own state).
"""
import math
import os
from dataclasses import dataclass

import numpy as np
from qiskit import QuantumCircuit, transpile

BASIS = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "sx", "sxdg", "rz", "rx", "ry", "p", "u", "cx", "cz", "cp",
         "rzz", "rxx", "ryy", "ccx", "swap"]
QC_FILES = ["tof_3", "tof_4", "tof_5", "tof_10", "barenco_tof_3", "barenco_tof_4", "barenco_tof_5", "barenco_tof_10",
            "mod5_4", "vbe_adder_3", "rc_adder_6", "adder_8", "gf2^4_mult", "gf2^5_mult", "mod_mult_55", "mod_red_21",
            "csla_mux_3", "qcla_com_7", "qcla_mod_7", "csum_mux_9", "grover_5", "qft_4", "ham15-low"]


@dataclass
class Bench:
    name: str
    family: str
    circuit: QuantumCircuit
    prep: QuantumCircuit = None


def flat(qc):
    """Same circuit written with standard one-/two-qubit gates and ccx only (no boxes, one register)."""
    out = transpile(qc, basis_gates=BASIS, optimization_level=0)
    res = QuantumCircuit(out.num_qubits, global_phase=out.global_phase)
    for ins in out.data:
        if ins.operation.name == "barrier":
            continue
        res.append(ins.operation, [out.find_bit(q).index for q in ins.qubits])
    return res


def from_qc(path):
    """QuantumCircuit of a Feynman .qc file (via qoptct.load_qc)."""
    import qoptct
    names, hl = qoptct.load_qc(path)
    qc = QuantumCircuit(len(names))
    for g, q in hl:
        getattr(qc, g)(*q)
    return qc


def x_prep(n, seed):
    rng = np.random.default_rng(seed)
    p = QuantumCircuit(n)
    for i in range(n):
        if rng.random() < 0.5:
            p.x(i)
    return p


def product_prep(state, n):
    """One-qubit gates preparing a product state given as a state vector (little endian)."""
    from qiskit.quantum_info import Statevector, partial_trace
    sv = Statevector(state)
    p = QuantumCircuit(n)
    for i in range(n):
        rho = partial_trace(sv, [j for j in range(n) if j != i]).data
        w, v = np.linalg.eigh(rho)
        a, b = v[:, -1]
        th = 2 * math.atan2(abs(b), abs(a))
        ph = (np.angle(b) - np.angle(a)) if abs(a) > 1e-9 and abs(b) > 1e-9 else 0.0
        p.u(th, ph, 0.0, i)
    return p


# ---- families ------------------------------------------------------------------------------------------------
def qft(n):
    from qiskit.synthesis import synth_qft_full
    return flat(synth_qft_full(n, do_swaps=False))


def qft_prep(n, k=None):
    from qiskit.quantum_info import Statevector
    if n > 14:
        return None
    k = (0b1011011101101 % (1 << n)) | 1 if k is None else k
    sv = Statevector.from_int(k, 1 << n).evolve(qft(n).inverse())
    return product_prep(sv.data, n)


def grover(n, iters=2):
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for _ in range(iters):
        qc.h(n - 1); qc.mcx(list(range(n - 1)), n - 1); qc.h(n - 1)          # oracle: phase flip of |1...1>
        qc.h(range(n)); qc.x(range(n))
        qc.h(n - 1); qc.mcx(list(range(n - 1)), n - 1); qc.h(n - 1)          # diffusion
        qc.x(range(n)); qc.h(range(n))
    return flat(qc)


def _graph(kind, n, seed=7):
    import rustworkx as rx
    if kind == "ring":
        return [(i, (i + 1) % n) for i in range(n)]
    if kind == "line":
        return [(i, i + 1) for i in range(n - 1)]
    if kind == "full":
        return [(i, j) for i in range(n) for j in range(i + 1, n)]
    if kind == "reg3":
        import random
        rnd = random.Random(seed)
        while True:  # random 3-regular graph by the pairing model
            pts = [i for i in range(n) for _ in range(3)]
            rnd.shuffle(pts)
            es = {(min(a, b), max(a, b)) for a, b in zip(pts[::2], pts[1::2])}
            if len(es) == 3 * n // 2 and all(a != b for a, b in es):
                g = rx.PyGraph(); g.add_nodes_from(range(n)); g.add_edges_from_no_data(sorted(es))
                if rx.is_connected(g):
                    return sorted(es)
    if kind.startswith("grid"):
        r, c = n
        es = []
        for i in range(r):
            for j in range(c):
                if j + 1 < c: es.append((i * c + j, i * c + j + 1))
                if i + 1 < r: es.append((i * c + j, (i + 1) * c + j))
        return es
    raise ValueError(kind)


def tfim(edges, n, steps=3, dt=0.2, h=1.0):
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for _ in range(steps):
        for a, b in edges:
            qc.rzz(2 * dt, a, b)
        for q in range(n):
            qc.rx(2 * dt * h, q)
    return qc


def heisenberg(edges, n, steps=2, dt=0.15):
    qc = QuantumCircuit(n)
    for q in range(0, n, 2):
        qc.x(q)  # Neel state
    for _ in range(steps):
        for a, b in edges:
            qc.rxx(2 * dt, a, b); qc.ryy(2 * dt, a, b); qc.rzz(2 * dt, a, b)
    return qc


def qaoa(edges, n, p=2, seed=3):
    rng = np.random.default_rng(seed)
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for _ in range(p):
        g, b = rng.uniform(0.2, 1.2, 2)
        for a, c in edges:
            qc.rzz(g, a, c)
        for q in range(n):
            qc.rx(2 * b, q)
    return qc


def vqe_full(n, reps=2, seed=5):
    rng = np.random.default_rng(seed)
    qc = QuantumCircuit(n)
    for r in range(reps + 1):
        for q in range(n):
            qc.ry(rng.uniform(0, 2 * math.pi), q); qc.rz(rng.uniform(0, 2 * math.pi), q)
        if r < reps:
            for a in range(n):
                for b in range(a + 1, n):
                    qc.cx(a, b)
    return qc


def qvol(n, depth=None, seed=9):
    from qiskit.quantum_info import random_unitary
    rng = np.random.default_rng(seed)
    qc = QuantumCircuit(n)
    for d in range(depth or n):
        perm = rng.permutation(n)
        for i in range(0, n - 1, 2):
            qc.unitary(random_unitary(4, seed=int(rng.integers(1 << 30))), [int(perm[i]), int(perm[i + 1])])
    return qc


def bv(n, secret=None):
    secret = secret or (0b101101110110101101101 % (1 << (n - 1))) | 1
    qc = QuantumCircuit(n)
    qc.x(n - 1); qc.h(range(n))
    for i in range(n - 1):
        if secret >> i & 1:
            qc.cx(i, n - 1)
    qc.h(range(n - 1))
    return qc


def qpe(m, theta=0.3141):
    from qiskit.synthesis import synth_qft_full
    qc = QuantumCircuit(m + 1)
    qc.x(m); qc.h(range(m))
    for j in range(m):
        qc.cp(2 * math.pi * theta * (1 << j), j, m)
    qc.compose(synth_qft_full(m, do_swaps=False, inverse=True), range(m), inplace=True)
    return flat(qc)


def _arith(kind, n):
    import qiskit.synthesis as syn
    if kind == "cdkm":
        return flat(syn.adder_ripple_c04(n, kind="full"))
    if kind == "vbe":
        return flat(syn.adder_ripple_v95(n, kind="full"))
    if kind == "draper":
        return flat(syn.adder_qft_d00(n, kind="half"))
    if kind == "mult_h18":
        return flat(syn.multiplier_cumulative_h18(n))
    if kind == "mult_qft":
        return flat(syn.multiplier_qft_r17(n))
    if kind == "cmp":
        from qiskit.circuit.library import IntegerComparatorGate
        qc = QuantumCircuit(2 * n)
        try:
            from qiskit.synthesis import synth_integer_comparator_2s
            return flat(synth_integer_comparator_2s(n, (1 << n) // 3 + 1))
        except Exception:
            qc = QuantumCircuit(n + 1)
            qc.append(IntegerComparatorGate(n, (1 << n) // 3 + 1), range(n + 1))
            return flat(qc)
    raise ValueError(kind)


def suite(kind="std", root=None):
    """Benchmark list.  kind: 'tiny' (tests), 'std' (full table), 'sim' (small circuits for the noisy simulation)."""
    root = root or os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "bench_qc")
    out = []

    def add(name, fam, fn, prep=None):
        try:
            c = fn()
        except Exception as e:  # a generator that is missing in this Qiskit version is skipped, loudly
            print(f"benchlib: skip {name}: {type(e).__name__}: {e}")
            return
        p = prep(c) if callable(prep) else prep
        out.append(Bench(name, fam, c, p))

    xp = lambda c: x_prep(c.num_qubits, c.num_qubits)  # noqa: E731
    if kind == "tiny":
        add("qft_5", "qft", lambda: qft(5), lambda c: qft_prep(5))
        add("add_cdkm_2", "arith", lambda: _arith("cdkm", 2), xp)
        add("qaoa_full_5", "qaoa", lambda: qaoa(_graph("full", 5), 5, 1))
        add("heis_ring_6", "ham", lambda: heisenberg(_graph("ring", 6), 6, 1))
        return out
    small = kind == "sim"
    for f in QC_FILES:
        p = os.path.join(root, f + ".qc")
        if os.path.exists(p):
            add("qc_" + f, "arith-qc", lambda p=p: from_qc(p), xp)
    for n in (2, 4, 6, 8):
        add(f"add_cdkm_{n}", "arith", lambda n=n: _arith("cdkm", n), xp)
    for n in (2, 4):
        add(f"add_vbe_{n}", "arith", lambda n=n: _arith("vbe", n), xp)
    for n in (3, 5, 8):
        add(f"add_draper_{n}", "arith", lambda n=n: _arith("draper", n), xp)
    for n in (4, 8):
        add(f"cmp_{n}", "arith", lambda n=n: _arith("cmp", n), xp)
    add("mult_h18_2", "arith", lambda: _arith("mult_h18", 2), xp)
    add("mult_qft_2", "arith", lambda: _arith("mult_qft", 2), xp)
    add("mult_qft_3", "arith", lambda: _arith("mult_qft", 3), xp)
    for n in (5, 8, 10, 14, 18, 22):
        add(f"qft_{n}", "qft", lambda n=n: qft(n), lambda c, n=n: qft_prep(n))
    for m in (6, 10):
        add(f"qpe_{m}", "qft", lambda m=m: qpe(m))
    for n in (4, 5, 6):
        add(f"grover_{n}", "grover", lambda n=n: grover(n, 1 if n < 5 else 2))
    for n in (8, 12, 16, 24):
        add(f"tfim_ring_{n}", "ham", lambda n=n: tfim(_graph("ring", n), n))
    for r, c in ((3, 3), (3, 4), (4, 5)):
        add(f"tfim_grid_{r}x{c}", "ham", lambda r=r, c=c: tfim(_graph("grid", (r, c)), r * c))
    for n in (8, 12, 16):
        add(f"heis_ring_{n}", "ham", lambda n=n: heisenberg(_graph("ring", n), n))
    for r, c in ((3, 3), (3, 4), (4, 4)):
        add(f"heis_grid_{r}x{c}", "ham", lambda r=r, c=c: heisenberg(_graph("grid", (r, c)), r * c))
    for n in (8, 12, 16, 22):
        add(f"qaoa_reg3_{n}", "qaoa", lambda n=n: qaoa(_graph("reg3", n), n, 2))
    for n in (5, 6, 8, 10, 12):
        add(f"qaoa_full_{n}", "qaoa", lambda n=n: qaoa(_graph("full", n), n, 1))
    for n in (6, 10):
        add(f"vqe_full_{n}", "misc", lambda n=n: vqe_full(n))
    for n in (6, 10):
        add(f"qvol_{n}", "misc", lambda n=n: qvol(n))
    for n in (10, 20):
        add(f"bv_{n}", "misc", lambda n=n: bv(n))
    if small:
        out = [b for b in out if b.circuit.num_qubits <= 10]
    return out
