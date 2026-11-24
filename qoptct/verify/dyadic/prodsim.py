# SPDX-License-Identifier: Apache-2.0
"""Product-state simulator: an independent numeric check at large n (not a proof; shares no code with the path sum).

On a basis input, the QFT, Draper adders, constant adders and QFT multipliers keep a product state: every
controlled phase has all of its qubits but at most one in a basis state, every CX / MCX has basis controls.  The
state is kept as one complex 2-vector per qubit plus a global factor; a gate that would entangle stops the run
(-> None).  Gate matrices are written from their definitions (Qiskit conventions), angles theta = pi * r.
"""
import cmath
import math

import numpy as np

__all__ = ['run_product', 'check_arith', 'check_qft']

_S2 = 1 / math.sqrt(2)
TOL = 1e-9


def _mat1(name, th):
    if name == 'h':
        return np.array([[_S2, _S2], [_S2, -_S2]], dtype=complex)
    if name == 'x':
        return np.array([[0, 1], [1, 0]], dtype=complex)
    if name in ('p', 'u1', 'phase'):
        return np.diag([1, cmath.exp(1j * th[0])])
    if name == 'rz':
        return np.diag([cmath.exp(-0.5j * th[0]), cmath.exp(0.5j * th[0])])
    if name in ('z', 's', 'sdg', 't', 'tdg'):
        a = {'z': math.pi, 's': math.pi / 2, 'sdg': -math.pi / 2, 't': math.pi / 4, 'tdg': -math.pi / 4}[name]
        return np.diag([1, cmath.exp(1j * a)])
    if name == 'ry':
        c, s = math.cos(th[0] / 2), math.sin(th[0] / 2)
        return np.array([[c, -s], [s, c]], dtype=complex)
    if name == 'rx':
        c, s = math.cos(th[0] / 2), math.sin(th[0] / 2)
        return np.array([[c, -1j * s], [-1j * s, c]])
    return None


def _basis(v):
    """0 or 1 if the qubit is |0> or |1> up to a phase (within TOL), else None."""
    if abs(v[1]) < TOL:
        return 0
    if abs(v[0]) < TOL:
        return 1
    return None


def run_product(c, bits):
    """Run a DCircuit on the basis state `bits` (list of 0/1 per qubit) -> (list of 2-vectors, global factor) or
    None when a gate would entangle (or is not handled)."""
    st = [np.array([1, 0], dtype=complex) if not b else np.array([0, 1], dtype=complex) for b in bits]
    g = cmath.exp(1j * math.pi * float(c.gphase))
    for nm, qs, rs in c.ops:
        th = [math.pi * float(r) for r in rs]
        if nm == 'gphase':
            g *= cmath.exp(1j * th[0])
            continue
        if nm == 'swap':
            st[qs[0]], st[qs[1]] = st[qs[1]], st[qs[0]]
            continue
        U = _mat1(nm, th)
        if U is not None:
            st[qs[0]] = U @ st[qs[0]]
            continue
        if nm in ('cp', 'cu1', 'cphase', 'mcp', 'mcphase', 'mcu1', 'cz', 'ccz', 'mcz', 'cs', 'csdg'):
            ang = th[0] if th else {'cs': math.pi / 2, 'csdg': -math.pi / 2}.get(nm, math.pi)
            free = [q for q in qs if _basis(st[q]) is None]
            if len(free) > 1:
                return None
            if any(_basis(st[q]) == 0 for q in qs if q not in free):
                continue                                   # some control is 0: identity
            if free:
                st[free[0]] = np.diag([1, cmath.exp(1j * ang)]) @ st[free[0]]
            else:
                g *= cmath.exp(1j * ang)
            continue
        if nm in ('cx', 'cnot', 'ccx', 'toffoli', 'mcx', 'mct', 'c3x', 'c4x'):
            cs, t = qs[:-1], qs[-1]
            vals = [_basis(st[q]) for q in cs]
            if None in vals:
                return None
            if all(vals):
                st[t] = np.array([st[t][1], st[t][0]])
            continue
        return None
    return st, g


def check_arith(c, spec, trials, rng):
    """ArithSpec / ClassicalSpec-like check on random basis inputs: the output must be the basis state of the
    classical map with amplitude of modulus 1.  -> (inputs checked, failures, entangling runs)."""
    n = c.n
    bad = ent = 0
    for _ in range(trials):
        bits = [rng.randint(0, 1) for _ in range(n)]
        zero = getattr(spec, 'zero_wires', ())
        for q in zero:
            bits[q] = 0
        r = run_product(c, bits)
        if r is None:
            ent += 1
            continue
        st, g = r
        b = sum(v << q for q, v in enumerate(bits))
        if spec.kind == 'arith':
            t = spec.target
            s = sum(bits[q] << j for j, q in enumerate(t))
            v = (s + spec.G(lambda q: bits[q])) % (1 << len(t))
            want = list(bits)
            for j, q in enumerate(t):
                want[q] = (v >> j) & 1
        else:
            o = spec.table[spec.index(b)]
            want = [(o >> q) & 1 for q in range(n)]
        amp = g
        for q in range(n):
            amp *= st[q][want[q]]
        if abs(abs(amp) - 1) > 1e-7:
            bad += 1
    return trials, bad, ent


def check_qft(c, spec, trials, rng):
    """QFTSpec check on random basis inputs: QFT|x> is the product state with qubit y_k = (|0> + e^{2 pi i X 2^k /
    2^m}|1>)/sqrt 2 on its output wire; compared qubit by qubit and in the total phase."""
    n = c.n
    bad = ent = 0
    w = spec.wires
    m = len(w)
    outp = spec.out_positions()
    for _ in range(trials):
        bits = [rng.randint(0, 1) for _ in range(n)]
        r = run_product(c, bits)
        if r is None:
            ent += 1
            continue
        st, g = r
        amp = g
        if spec.inverse:
            # inverse QFT of a basis state: output is a product state with phases e^{-2 pi i Y 2^j / 2^m}
            Y = sum(bits[outp[k]] << k for k in range(m))
            exp = {w[j]: np.array([1, cmath.exp(-2j * math.pi * ((Y << j) % (1 << m)) / (1 << m))]) * _S2
                   for j in range(m)}
        else:
            X = sum(bits[w[j]] << j for j in range(m))
            exp = {outp[k]: np.array([1, cmath.exp(2j * math.pi * ((X << k) % (1 << m)) / (1 << m))]) * _S2
                   for k in range(m)}
        for q in range(n):
            if q in exp:
                amp *= np.vdot(exp[q], st[q])
            else:
                amp *= st[q][bits[q]]
        if abs(abs(amp) - 1) > 1e-7:
            bad += 1
    return trials, bad, ent
