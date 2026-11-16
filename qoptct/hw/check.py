"""Equivalence check of a compiled (laid out + routed) circuit against the logical circuit.

Claim that is checked:  for every state |psi> of the logical qubits, placed on the physical qubits ``init_map``
(all other physical qubits in |0>), the compiled circuit produces U|psi> on the physical qubits ``final_map`` and
leaves all other qubits in |0>  (up to a global phase).

Method (:func:`verify_compiled`): for a few random product states |psi> the circuit
``prepare(psi) -> compiled -> inverse(logical) on final_map -> inverse(prepare) on final_map`` is simulated as a
state vector on the used physical qubits; the amplitude of |0...0> must have modulus 1.  If the two circuits differ,
a random product state passes with probability 0 (the states on which two different unitaries agree up to a phase
form a set of measure zero), so the test is a numerical proof up to floating-point tolerance; it is exponential in
the number of used qubits (limit ``max_qubits``, default 26; above that the result is None = not checked).
"""
import math

import numpy as np
from qiskit import QuantumCircuit


def _compact(comp, init_map, final_map):
    used = set(init_map) | set(final_map)
    for ins in comp.data:
        if ins.operation.name in ("barrier", "delay"):
            continue
        used.update(comp.find_bit(q).index for q in ins.qubits)
    order = sorted(used)
    return order, {p: i for i, p in enumerate(order)}


def _zero_amp(circ, k):
    try:
        from qiskit_aer import AerSimulator
    except ImportError:
        AerSimulator = None
    if AerSimulator is None or k <= 10:
        from qiskit.quantum_info import Statevector
        from qiskit import transpile
        if k > 20:
            raise RuntimeError("install qiskit-aer to verify circuits on more than 20 used qubits")
        c = transpile(circ, basis_gates=["u", "cx", "cz", "ecr", "rz", "sx", "x", "h", "cp", "rzz", "swap"],
                      optimization_level=0)
        return complex(Statevector.from_label("0" * k).evolve(c).data[0])
    from qiskit import transpile
    sim = AerSimulator(method="statevector", precision="single" if k > 24 else "double", max_parallel_threads=1,
                       fusion_enable=True)
    c = circ.copy()
    c.save_amplitudes([0])
    c = transpile(c, sim, optimization_level=0)
    res = sim.run(c, shots=1).result()
    return complex(res.data(0)["amplitudes"][0])


def verify_compiled(logical, comp, init_map, final_map, trials=2, max_qubits=26, seed=1234):
    """True / False / None (too many used qubits).  ``logical``: circuit on n qubits without measurements;
    ``comp``: physical circuit; ``init_map[i]`` / ``final_map[i]``: physical qubit of logical qubit i at the
    start / at the end."""
    n = logical.num_qubits
    init_map, final_map = list(init_map), list(final_map)
    if len(init_map) != n or len(final_map) != n or len(set(init_map)) != n or len(set(final_map)) != n:
        return False
    order, idx = _compact(comp, init_map, final_map)
    k = len(order)
    if k > max_qubits:
        return None
    body = QuantumCircuit(k)
    for ins in comp.data:
        if ins.operation.name in ("barrier", "delay"):
            continue
        if ins.operation.name in ("measure", "reset") or ins.clbits:
            raise ValueError("verify_compiled needs circuits without measurements")
        body.append(ins.operation, [idx[comp.find_bit(q).index] for q in ins.qubits])
    inv = QuantumCircuit(k)
    inv.compose(logical.inverse(), qubits=[idx[p] for p in final_map], inplace=True)
    rng = np.random.default_rng(seed)
    tol = 2e-3 if k > 24 else 1e-6
    for _ in range(trials):
        ang = rng.uniform(0, 2 * math.pi, size=(n, 3))
        c = QuantumCircuit(k)
        for i in range(n):
            c.u(ang[i, 0], ang[i, 1], ang[i, 2], idx[init_map[i]])
        c.compose(body, inplace=True)
        c.compose(inv, inplace=True)
        for i in range(n):
            c.u(-ang[i, 0], -ang[i, 2], -ang[i, 1], idx[final_map[i]])
        a = _zero_amp(c, k)
        if abs(abs(a) - 1.0) > tol:
            return False
    return True


def maps_from_transpiled(tc, n):
    """(init_map, final_map) of the first n logical qubits of a circuit returned by ``qiskit.transpile``."""
    lay = tc.layout
    if lay is None:
        return list(range(n)), list(range(n))
    init = lay.initial_index_layout(filter_ancillas=True)[:n]
    fin = lay.final_index_layout(filter_ancillas=True)[:n]
    return list(init), list(fin)
