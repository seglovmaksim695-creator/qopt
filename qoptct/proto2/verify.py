"""Проверка равенства: (1) PyZX verify_equality (True = точное доказательство переписыванием ZX);
(2) численно: случайные состояния, сравнение с точностью до общей фазы (numpy, до ~26 кубитов)."""
import numpy as np

S2 = 1 / np.sqrt(2)
PH = {"t": np.exp(1j * np.pi / 4), "tdg": np.exp(-1j * np.pi / 4), "s": 1j, "sdg": -1j, "z": -1}


def run(c, n, psi):
    psi = psi.reshape([2] * n)
    for nm, qs in c:
        ax = [q for q in qs]
        if nm == "h":
            q = ax[0]
            a = np.take(psi, 0, axis=q); b = np.take(psi, 1, axis=q)
            psi = np.stack([(a + b) * S2, (a - b) * S2], axis=q)
        elif nm == "x":
            psi = np.flip(psi, axis=ax[0])
        elif nm in PH:
            sl = [slice(None)] * n; sl[ax[0]] = 1
            psi[tuple(sl)] *= PH[nm]
        elif nm == "cx":
            c_, t = ax
            sl = [slice(None)] * n; sl[c_] = 1
            sub = psi[tuple(sl)]
            tt = t if t < c_ else t - 1
            psi[tuple(sl)] = np.flip(sub, axis=tt).copy()
        else:
            raise ValueError(nm)
    return psi.reshape(-1)


def numeric_equal(c1, c0, n, seeds=2, dtype=np.complex128):
    ph = None
    for s in range(seeds):
        r = np.random.default_rng(s)
        v = (r.standard_normal(2 ** n) + 1j * r.standard_normal(2 ** n)).astype(dtype)
        v /= np.linalg.norm(v)
        a = run(c0, n, v.copy()); b = run(c1, n, v.copy())
        o = np.vdot(a, b)
        tol = 1e-11 if dtype == np.complex128 else 1e-3  # Г1: ошибка на одном базисном состоянии ~2^-n
        if abs(abs(o) - 1) > tol or (ph is not None and abs(o - ph) > tol):
            return False
        ph = o
    return True


def to_zx(c, n):
    import pyzx as zx
    z = zx.Circuit(n)
    for nm, qb in c:
        if nm == "h": z.add_gate("HAD", qb[0])
        elif nm == "x": z.add_gate("NOT", qb[0])
        elif nm == "cx": z.add_gate("CNOT", qb[0], qb[1])
        elif nm in ("t", "tdg"): z.add_gate("T", qb[0], adjoint=(nm == "tdg"))
        elif nm in ("s", "sdg"): z.add_gate("S", qb[0], adjoint=(nm == "sdg"))
        elif nm == "z": z.add_gate("Z", qb[0])
    return z


def zx_equal(c1, c0, n):
    return to_zx(c1, n).verify_equality(to_zx(c0, n))
