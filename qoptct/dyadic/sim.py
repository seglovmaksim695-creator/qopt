# SPDX-License-Identifier: Apache-2.0
"""Independent dense simulator for qoptct.dyadic (cross-checks and counterexample confirmation; not a proof).

It shares nothing with the path-sum engine except the gate list of a DCircuit: every gate is applied by its matrix,
written here from the gate's definition (Qiskit conventions, qubit 0 = least significant bit), with complex128
numbers and angle theta = pi * r (r the exact dyadic angle of the circuit).  Specifications are built from their
mathematical definition (QFT matrix, permutation of an ArithSpec / PermutationSpec), not from their path sums.
"""
import cmath
import math

import numpy as np

__all__ = ['unitary', 'spec_unitary', 'compare', 'apply_op', 'MAX_DENSE']

MAX_DENSE = 12          # qubits for full unitaries
_S2 = 1 / math.sqrt(2)


def _u(th, ph, la):
    c, s = math.cos(th / 2), math.sin(th / 2)
    return np.array([[c, -cmath.exp(1j * la) * s], [cmath.exp(1j * ph) * s, cmath.exp(1j * (ph + la)) * c]])


def _one(name, th):
    if name in ('id', 'i'):
        return np.eye(2, dtype=complex)
    if name == 'x':
        return np.array([[0, 1], [1, 0]], dtype=complex)
    if name == 'y':
        return np.array([[0, -1j], [1j, 0]])
    if name == 'z':
        return np.diag([1, -1]).astype(complex)
    if name == 'h':
        return np.array([[_S2, _S2], [_S2, -_S2]], dtype=complex)
    if name == 's':
        return np.diag([1, 1j])
    if name == 'sdg':
        return np.diag([1, -1j])
    if name == 't':
        return np.diag([1, cmath.exp(1j * math.pi / 4)])
    if name == 'tdg':
        return np.diag([1, cmath.exp(-1j * math.pi / 4)])
    if name == 'sx':
        return 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
    if name == 'sxdg':
        return 0.5 * np.array([[1 - 1j, 1 + 1j], [1 + 1j, 1 - 1j]])
    if name in ('p', 'u1', 'phase'):
        return np.diag([1, cmath.exp(1j * th[0])])
    if name == 'rz':
        return np.diag([cmath.exp(-0.5j * th[0]), cmath.exp(0.5j * th[0])])
    if name == 'rx':
        c, s = math.cos(th[0] / 2), math.sin(th[0] / 2)
        return np.array([[c, -1j * s], [-1j * s, c]])
    if name == 'ry':
        c, s = math.cos(th[0] / 2), math.sin(th[0] / 2)
        return np.array([[c, -s], [s, c]], dtype=complex)
    if name in ('u', 'u3'):
        return _u(th[0], th[1], th[2])
    if name == 'u2':
        return _u(math.pi / 2, th[0], th[1])
    raise ValueError('no matrix for %s' % name)


# controlled families: name -> (base single-qubit gate, number of controls or None = all but the last qubit)
_CTRL = {'cx': ('x', 1), 'cnot': ('x', 1), 'cy': ('y', 1), 'cz': ('z', 1), 'cs': ('s', 1), 'csdg': ('sdg', 1),
         'csx': ('sx', 1), 'cp': ('p', 1), 'cu1': ('p', 1), 'cphase': ('p', 1), 'crz': ('rz', 1),
         'ccx': ('x', 2), 'toffoli': ('x', 2), 'ccz': ('z', 2), 'c3x': ('x', 3), 'c4x': ('x', 4),
         'mcx': ('x', None), 'mct': ('x', None), 'mcz': ('z', None), 'mcp': ('p', None), 'mcphase': ('p', None),
         'mcu1': ('p', None)}


def _axis(n, q):
    return n - 1 - q


def _apply1(st, n, U, q, ctrl=()):
    idx = [slice(None)] * st.ndim
    for c in ctrl:
        idx[_axis(n, c)] = 1
    idx = tuple(idx)
    sub = st[idx]
    ax = _axis(n, q) - sum(1 for c in ctrl if _axis(n, c) < _axis(n, q))
    st[idx] = np.moveaxis(np.tensordot(U, sub, axes=([1], [ax])), 0, ax)


def _swap(st, n, a, b, ctrl=()):
    idx = [slice(None)] * st.ndim
    for c in ctrl:
        idx[_axis(n, c)] = 1
    idx = tuple(idx)
    sub = st[idx]
    pa = _axis(n, a) - sum(1 for c in ctrl if _axis(n, c) < _axis(n, a))
    pb = _axis(n, b) - sum(1 for c in ctrl if _axis(n, c) < _axis(n, b))
    st[idx] = np.swapaxes(sub, pa, pb).copy()


def apply_op(st, n, name, qs, rs):
    """Apply one gate to the array st of shape (2,)*n (+ extra axes), axis n-1-q = qubit q."""
    th = [math.pi * float(r) for r in rs]
    if name == 'gphase':
        st *= cmath.exp(1j * th[0])
        return
    if name == 'swap':
        _swap(st, n, qs[0], qs[1])
        return
    if name in ('cswap', 'fredkin'):
        _swap(st, n, qs[1], qs[2], ctrl=(qs[0],))
        return
    if name == 'dcx':
        _apply1(st, n, _one('x', ()), qs[1], ctrl=(qs[0],))
        _apply1(st, n, _one('x', ()), qs[0], ctrl=(qs[1],))
        return
    if name == 'rzz':
        a, b = _axis(n, qs[0]), _axis(n, qs[1])
        for va in (0, 1):
            for vb in (0, 1):
                idx = [slice(None)] * st.ndim
                idx[a], idx[b] = va, vb
                st[tuple(idx)] *= cmath.exp((0.5j if va != vb else -0.5j) * th[0])
        return
    if name in _CTRL:
        base, k = _CTRL[name]
        ctrl = tuple(qs[:-1]) if k is None else tuple(qs[:k])
        _apply1(st, n, _one(base, th), qs[-1], ctrl=ctrl)
        return
    _apply1(st, n, _one(name, th), qs[0])


def unitary(c):
    """Dense unitary of a DCircuit (n <= MAX_DENSE); column b = U|b>."""
    n = c.n
    if n > MAX_DENSE:
        raise ValueError('too many qubits for the dense simulator')
    d = 1 << n
    st = np.eye(d, dtype=complex).reshape((2,) * n + (d,))
    for nm, qs, rs in c.ops:
        apply_op(st, n, nm, qs, rs)
    U = st.reshape(d, d)
    if c.gphase:
        U = U * cmath.exp(1j * math.pi * float(c.gphase))
    return U


def _bit(b, q):
    return (b >> q) & 1


def spec_unitary(spec, n):
    """Dense unitary of a specification from its definition."""
    d = 1 << n
    U = np.zeros((d, d), dtype=complex)
    if spec.kind == 'qft':
        w = spec.wires
        m = len(w)
        outp = spec.out_positions()
        M = 1 << m
        for b in range(d):
            rest = b
            for q in w:
                rest &= ~(1 << q)
            if not spec.inverse:
                X = sum(_bit(b, w[j]) << j for j in range(m))
                for Y in range(M):
                    o = rest
                    for k in range(m):
                        if (Y >> k) & 1:
                            o |= 1 << outp[k]
                    U[o, b] += cmath.exp(2j * math.pi * X * Y / M) / math.sqrt(M)
            else:
                Y = sum(_bit(b, outp[k]) << k for k in range(m))
                for X in range(M):
                    o = rest
                    for j in range(m):
                        if (X >> j) & 1:
                            o |= 1 << w[j]
                    U[o, b] += cmath.exp(-2j * math.pi * X * Y / M) / math.sqrt(M)
        return U
    if spec.kind == 'arith':
        t = spec.target
        m = len(t)
        for b in range(d):
            s = sum(_bit(b, t[j]) << j for j in range(m))
            v = (s + spec.G(lambda q: _bit(b, q))) % (1 << m)
            o = b
            for j in range(m):
                o = (o & ~(1 << t[j])) | (((v >> j) & 1) << t[j])
            U[o, b] = 1
        return U
    if spec.kind == 'classical':
        for b in range(d):
            if any((b >> q) & 1 for q in spec.zero_wires):
                continue                      # columns outside the domain stay 0 (compared only on the domain)
            U[spec.table[spec.index(b)], b] = 1
        return U
    raise ValueError('unknown specification kind')


def compare(U, V, tol=1e-8):
    """-> (equal up to a global phase, global phase angle in radians or None, max deviation)."""
    M = V.conj().T @ U
    d = M.shape[0]
    k = int(np.argmax(np.abs(np.diag(M))))
    lam = M[k, k]
    if abs(lam) < 0.5:
        return False, None, float(np.max(np.abs(M - np.eye(d))))
    lam = lam / abs(lam)
    dev = float(np.max(np.abs(M - lam * np.eye(d))))
    return dev < tol, cmath.phase(lam), dev
