# SPDX-License-Identifier: Apache-2.0
"""Verification of cxpass results.

* ``check_vs_ref``  -- state-vector check against the .qc reference in the strict or clean model, with the same figure of merit
  as qoptct.verify_qasm.check (random data states, sampled measurement outcomes, ancillas must end in a basis state), but with
  support for XOR conditions and ``reset``.  Up to MAX_SV qubits.
* ``check_vs_prog`` -- equivalence with the input program (used when no reference is given): both programs are run on random
  states over the wires on which the input program is deterministic; the full output states must agree up to a phase.
* ``check_big``     -- strict model, more than MAX_SV qubits: path-sum verifier of lab agent c10 (qoptct.contrib.c10.verify_big).
* ``verify``        -- dispatcher; prefers ``qoptct.verify.prove`` when that module is installed.
"""
import os
import tempfile
import time

import numpy as np

from . import ir

MAX_SV = 22
_S2 = 1 / np.sqrt(2)
_PH = {"z": -1.0, "s": 1j, "sdg": -1j, "t": np.exp(1j * np.pi / 4), "tdg": np.exp(-1j * np.pi / 4)}


def _ax(psi, N, q):
    return psi.reshape(1 << q, 2, 1 << (N - q - 1))


def apply(psi, N, g, qs):
    """In-place style application on a flat state vector (qubit 0 = most significant axis)."""
    if g == "h":
        v = _ax(psi, N, qs[0]); a = v[:, 0, :].copy(); b = v[:, 1, :]
        v[:, 0, :] = (a + b) * _S2; v[:, 1, :] = (a - b) * _S2
    elif g == "x":
        v = _ax(psi, N, qs[0]); a = v[:, 0, :].copy(); v[:, 0, :] = v[:, 1, :]; v[:, 1, :] = a
    elif g in _PH:
        _ax(psi, N, qs[0])[:, 1, :] *= _PH[g]
    elif g == "cx":
        c, t = qs; v = psi.reshape([2] * N); idx = [slice(None)] * N; idx[c] = 1
        sub = v[tuple(idx)]; ax = t - (1 if c < t else 0)
        v[tuple(idx)] = np.flip(sub, axis=ax).copy()
    elif g == "cz":
        v = psi.reshape([2] * N); idx = [slice(None)] * N; idx[qs[0]] = 1; idx[qs[1]] = 1
        v[tuple(idx)] *= -1
    else:
        raise ValueError(g)
    return psi


def _measure(psi, N, q, rng):
    v = _ax(psi, N, q)
    p1 = float(np.sum(np.abs(v[:, 1, :]) ** 2)); tot = float(np.sum(np.abs(psi) ** 2))
    m = int(rng.random() * tot < p1)
    v[:, 1 - m, :] = 0
    nrm = np.sqrt(p1 if m else tot - p1)
    if nrm > 0:
        psi /= nrm
    return m


def run(N, items, psi, rng):
    bits = {}
    for it in items:
        if it[0] == "G":
            apply(psi, N, it[1], it[2])
        elif it[0] == "M":
            bits[it[2]] = _measure(psi, N, it[1], rng)
        elif it[0] == "R":
            if _measure(psi, N, it[1], rng):
                apply(psi, N, "x", (it[1],))
        else:
            ok = all(sum(bits[b] for b in g) & 1 for g in ir.pcond(it[1]))
            if ok:
                for g, qs in it[2]:
                    apply(psi, N, g, qs)
    return psi


def _ref(ref_qc):
    import qoptct
    import qoptct.verify_qasm  # noqa: F401  (puts the flat internal modules evo_lib, gidney on sys.path)
    import evo_lib
    import gidney
    names, hl = qoptct.load_qc(ref_qc)
    ins = gidney.inputs_of(ref_qc)
    clean_anc = [i for i, nm in enumerate(names) if ins is not None and nm not in ins]
    return len(names), evo_lib.lower_std(hl), clean_anc


def check_vs_ref(text, ref_qc, clean=False, trials=12, seed=7, deadline=None):
    """Worst fidelity (1.0 = equal) of the program against the reference, as in verify_qasm.check(strict=True)."""
    N, items = ir.parse(text)
    n, ref, clean_anc = _ref(ref_qc)
    if N < n or N > MAX_SV:
        return None
    rng = np.random.default_rng(seed); worst = 1.0
    for k in range(trials):
        if deadline and time.time() > deadline and k >= 3:
            break
        v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
        if clean:
            vv = v.reshape([2] * n)
            for a in clean_anc:
                sl = [slice(None)] * n; sl[a] = 1; vv[tuple(sl)] = 0
            v = vv.reshape(-1)
        v /= np.linalg.norm(v)
        want = v.copy()
        for g, q in ref:
            want = apply(want, n, g, tuple(q)) if g in ("h", "x", "cx", "cz") or g in _PH else _gidney_apply(want, n, g, q)
        psi = np.zeros(2 ** N, complex); psi.reshape(2 ** n, 2 ** (N - n))[:, 0] = v
        out = run(N, items, psi, rng)
        o = out.reshape(2 ** n, 2 ** (N - n))
        col = int(np.argmax(np.sum(np.abs(o) ** 2, axis=0)))
        got = o[:, col]
        nrm = np.linalg.norm(out)
        if nrm == 0:
            return 0.0
        worst = min(worst, abs(np.vdot(want, got)) / nrm)
    return worst


def _gidney_apply(psi, n, g, q):
    import qoptct.verify_qasm  # noqa: F401
    import gidney
    return gidney.apply(psi, n, g, tuple(q))


def _random_state(N, D, rng):
    psi = np.zeros(2 ** N, complex); v = psi.reshape([2] * N)
    sub = rng.normal(size=[2] * len(D)) + 1j * rng.normal(size=[2] * len(D))
    idx = [0] * N
    for q in D:
        idx[q] = slice(None)
    v[tuple(idx)] = sub
    psi /= np.linalg.norm(psi)
    return psi


def _overlap(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return abs(np.vdot(a, b)) / (na * nb)


def _is_det(N, items, T, seed, tol, inputs=1, runs=4):
    """Same output state (up to a phase) for `runs` independent samplings of the measurement outcomes."""
    for s in range(inputs):
        rng = np.random.default_rng(seed + 17 * s + len(T))
        psi = _random_state(N, T, rng)
        a = run(N, items, psi.copy(), np.random.default_rng(seed + 100 + s))
        for r in range(1, runs):
            b = run(N, items, psi.copy(), np.random.default_rng(seed + 100 * (r + 1) + s))
            if _overlap(a, b) < 1 - tol:
                return False
    return True


def deterministic_wires(N, items, candidates, seed=3, tol=1e-9, greedy=True, deadline=None):
    """Set D of wires on which the program is deterministic for random input states (other wires start in |0>).
    The longest deterministic prefix of `candidates` is found by bisection (qoptct puts data wires first); with
    greedy=True the remaining candidates are then tried one by one (clean model: clean wires inside the data range)."""
    cand = list(candidates)
    lo, hi = 0, len(cand)
    if _is_det(N, items, cand, seed, tol):
        lo = hi
    else:
        hi -= 1
        while lo < hi:  # largest k with cand[:k] deterministic (assumes monotonicity; checked again below)
            mid = (lo + hi + 1) // 2
            if _is_det(N, items, cand[:mid], seed, tol):
                lo = mid
            else:
                hi = mid - 1
    D = cand[:lo]
    if greedy:
        for q in cand[lo + 1:]:
            if deadline and time.time() > deadline:
                break
            if _is_det(N, items, D + [q], seed, tol):
                D = D + [q]
    while D and not _is_det(N, items, D, seed + 1, tol, inputs=2):
        D = D[:-1]
    return D


def check_vs_prog(text, ref_text, D, trials=10, seed=11, deadline=None):
    N, items = ir.parse(text); N0, items0 = ir.parse(ref_text)
    if N != N0 or N > MAX_SV:
        return None
    worst = 1.0
    for k in range(trials):
        if deadline and time.time() > deadline and k >= 3:
            break
        rng = np.random.default_rng(seed + k)
        psi = _random_state(N, D, rng)
        a = run(N, items0, psi.copy(), np.random.default_rng(seed + 1000 + k))
        b = run(N, items, psi.copy(), np.random.default_rng(seed + 2000 + k))
        worst = min(worst, _overlap(a, b))
    return worst


def check_big(text, ref_qc, m_trials=12, x_trials=3):
    from qoptct.contrib.c10 import verify_big
    with tempfile.NamedTemporaryFile("w", suffix=".qasm", delete=False) as f:
        f.write(text); p = f.name
    try:
        r = verify_big.check(p, ref_qc, x_trials=x_trials, m_trials=m_trials, verbose=False)
    finally:
        os.unlink(p)
    return 1.0 if r.get("ok") else 0.0


def _official(text, ref_qc, clean, budget):
    """qoptct.verify.prove (verify team), if available: True / False / None (unknown or missing)."""
    try:
        from qoptct.verify import prove
    except Exception:
        return None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".qasm", delete=False) as f:
            f.write(text); p = f.name
        try:
            r = prove(p, ref_qc, clean=clean, budget_sec=budget)
        finally:
            os.unlink(p)
        st = (r or {}).get("status")
        return True if st == "PROVED" else (False if st == "REFUTED" else None)
    except Exception:
        return None


TOL = 1e-9


def _dense(text, ref_qc, clean, trials, seed=20261114):
    """qoptct.verify.check_dense (independent state-vector simulator, reads XOR conditions and reset):
    True / False, or None when that module is missing or the circuit is too large for it."""
    try:
        from qoptct.verify import check_dense
    except Exception:
        return None
    try:
        r = check_dense(text, ref_qc, clean=clean, trials=trials, seed=seed)
    except Exception:
        return None
    if r.get("ok"):
        return True
    return None if r.get("F") is None else False


def verify(text, clean, ref_qc=None, ref_text=None, D=None, budget=30.0, trials=12, official=True, strong=False):
    """True only if the program was checked equal to the reference (ref_qc) or to the input program (ref_text, over wires D)."""
    try:
        return _verify(text, clean, ref_qc, ref_text, D, budget, trials, official, strong)
    except Exception:
        return False


def _verify(text, clean, ref_qc, ref_text, D, budget, trials, official, strong=False):
    deadline = time.time() + max(1.0, budget)
    N, _ = ir.parse(text)
    if ref_qc:
        if official:
            r = _official(text, ref_qc, clean, budget)
            if r is not None:
                return r
        if N <= MAX_SV:
            f = check_vs_ref(text, ref_qc, clean=clean, trials=trials, deadline=deadline)
            if not (f is not None and f > 1 - TOL):
                return False
            # final check: a second, independent simulator with many more sampled measurement branches
            return _dense(text, ref_qc, clean, 60) is not False if strong else True
        if not clean:
            return check_big(text, ref_qc, m_trials=40 if strong else 12, x_trials=4 if strong else 3) > 1 - TOL
        return _dense(text, ref_qc, clean, 60 if strong else 20) is True   # clean, 23..26 qubits (else unknown -> no)
    if ref_text is not None and D is not None and N <= MAX_SV:
        f = check_vs_prog(text, ref_text, D, trials=trials, deadline=deadline)
        return f is not None and f > 1 - TOL
    return False
