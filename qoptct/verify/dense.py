# SPDX-License-Identifier: Apache-2.0
"""Dense sampled check (state vector, up to ~24 qubits) using the fast simulator qoptct.verify.fastsim.

Same property as qoptct.verify_qasm.check (random input states, Born-rule measurement outcomes, extra qubits must
end in one basis state, no renormalisation in the strict comparison), plus reset and XOR/NOT conditions.
"""
import numpy as np

from . import fastsim
from ._parse import QasmFile, parse_qasm_text, parse_qc, eval_cond, clean_wires

MAX_QUBITS = 26


def _measure(psi, N, q, rng, forced=None):
    p = psi.reshape(1 << q, 2, 1 << (N - q - 1))
    p1 = float(np.sum(np.abs(p[:, 1, :]) ** 2))
    tot = float(np.sum(np.abs(psi) ** 2))
    r = forced if forced is not None else int(rng.random() < p1 / tot)
    p[:, 1 - r, :] = 0
    pr = p1 if r else tot - p1
    if pr <= 0:
        return psi, r, 0.0
    psi /= np.sqrt(pr / tot)
    return psi, r, pr


def run(ops, N, psi, rng):
    c = {}
    pend = []
    for op in ops:
        if op[0] == 'G' and op[1] != 'reset':
            pend.append((op[1], op[2]))
            continue
        if pend:
            psi = fastsim.apply_gates(psi, N, pend)
            pend = []
        if op[0] == 'G':                       # reset = measure, then flip to |0>
            q = op[2][0]
            psi, r, _ = _measure(psi, N, q, rng)
            if r:
                psi = fastsim.apply_gates(psi, N, [('x', (q,))])
        elif op[0] == 'M':
            psi, r, _ = _measure(psi, N, op[1], rng)
            c[op[2]] = r
        elif eval_cond(op[1], c):
            psi = fastsim.apply_gates(psi, N, op[2])
    if pend:
        psi = fastsim.apply_gates(psi, N, pend)
    return psi, c


def check_dense(qasm, qc, clean=False, trials=30, seed=1, strict=True):
    """-> dict(ok, F) with F the worst fidelity over `trials` random input states."""
    with QasmFile(qasm) as path:
        txt = open(path).read()
    try:
        N, nb, ops = parse_qasm_text(txt)
        names, ins, ref = parse_qc(qc)
    except (ValueError, KeyError) as e:
        return dict(ok=False, F=0.0, error=str(e))
    n = len(names)
    if N < n:
        return dict(ok=False, F=0.0, error='fewer qubits than the reference')
    if N > MAX_QUBITS:
        return dict(ok=False, F=None, error='too many qubits for the dense check (%d > %d)' % (N, MAX_QUBITS))
    anc = clean_wires(names, ins, clean)
    rng = np.random.default_rng(seed)
    worst = 1.0
    for _ in range(trials):
        v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
        vv = v.reshape([2] * n)
        for a in anc:
            sl = [slice(None)] * n
            sl[a] = 1
            vv[tuple(sl)] = 0
        v = vv.reshape(-1)
        v /= np.linalg.norm(v)
        want = fastsim.apply_gates(v.copy(), n, ref)
        psi = np.zeros([2] * N, complex)
        psi[(slice(None),) * n + (0,) * (N - n)] = v.reshape([2] * n)
        out, _ = run(ops, N, psi.reshape(-1), rng)
        o = out.reshape([2] * N)
        idx = [slice(None)] * N
        for q in range(n, N):
            idx[q] = int(np.argmax([np.sum(np.abs(np.take(o, b, axis=q)) ** 2) for b in (0, 1)]))
        got = o[tuple(idx)].reshape(-1)
        if strict:
            f = abs(np.vdot(want, got)) / np.linalg.norm(out)
        else:
            f = abs(np.vdot(want, got / np.linalg.norm(got)))
        worst = min(worst, float(f))
    return dict(ok=bool(worst > 1 - 1e-9), F=worst, trials=trials)
