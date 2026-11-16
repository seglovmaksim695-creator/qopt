# SPDX-License-Identifier: Apache-2.0
"""dyn.lower: AND / MZ form -> Qiskit circuits.

    to_qiskit(C, erase)   the circuit to run. erase = 'meas' | 'unitary' | function(site_index, op) -> one of them.
        AND(a, b, t)            : 3 CX + 4 T (t must be |0>; exact, no relative phase on that subspace)
        MZ(t, corr), 'meas'     : h t; measure t; if (m) {x t; corr}   (dynamic circuit; two-qubit gates only if m = 1)
        MZ(t, corr), 'unitary'  : corr = [cz a b] -> the inverse AND gadget (3 CX); otherwise one (multi-)controlled X
                                  per monomial of corr.
        via_target=True: the correction cz(a, b) of a 'meas' erasure is routed through the just-measured wire t
        (x t; cx a t; cz t b; cx a t): t sits between a and b after the AND, a and b are usually not adjacent.
    reference(C)          the ideal permutation (AND -> ccx, MZ -> (m)cx per monomial), used for checks and as the
                          noise-free inverse in the superposition test.
"""
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister


def _and(qc, a, b, t):
    qc.h(t); qc.t(t); qc.cx(b, t); qc.tdg(t); qc.cx(a, t); qc.t(t); qc.cx(b, t); qc.tdg(t); qc.h(t); qc.sdg(t)


def _and_inv(qc, a, b, t):
    qc.s(t); qc.h(t); qc.t(t); qc.cx(b, t); qc.tdg(t); qc.cx(a, t); qc.t(t); qc.cx(b, t); qc.tdg(t); qc.h(t)


def _mono_x(qc, term, t):
    ws = list(term[1:])
    if len(ws) == 1:
        qc.cx(ws[0], t)
    elif len(ws) == 2:
        qc.ccx(ws[0], ws[1], t)
    else:
        qc.mcx(ws, t)


def _mono_z(qc, term):
    ws = list(term[1:])
    if len(ws) == 1:
        qc.z(ws[0])
    elif len(ws) == 2:
        qc.cz(ws[0], ws[1])
    elif len(ws) == 3:
        qc.ccz(*ws)
    else:
        qc.h(ws[-1]); qc.mcx(ws[:-1], ws[-1]); qc.h(ws[-1])


def n_sites(C):
    return sum(1 for o in C.ops if o[0] == 'MZ')


def to_qiskit(C, erase='meas', via_target=True):
    q = QuantumRegister(C.N, 'q'); qc = QuantumCircuit(q)
    pick = erase if callable(erase) else (lambda i, op: erase)
    modes = []; k = 0
    for op in C.ops:
        if op[0] == 'MZ':
            modes.append(pick(k, op)); k += 1
    nmeas = sum(1 for m in modes if m == 'meas')
    creg = None
    if nmeas:
        creg = ClassicalRegister(nmeas, 'm'); qc.add_register(creg)
    k = 0; km = 0
    for op in C.ops:
        g = op[0]
        if g == 'x':
            qc.x(op[1])
        elif g == 'z':
            qc.z(op[1])
        elif g == 'cx':
            qc.cx(op[1], op[2])
        elif g == 'cz':
            qc.cz(op[1], op[2])
        elif g == 'ccx':
            qc.ccx(op[1], op[2], op[3])
        elif g == 'AND':
            _and(qc, op[1], op[2], op[3])
        elif g == 'MZ':
            t, corr = op[1], op[2]; mode = modes[k]; k += 1
            if mode == 'unitary':
                if len(corr) == 1 and corr[0][0] == 'cz':
                    _and_inv(qc, corr[0][1], corr[0][2], t)
                else:
                    for c in corr:
                        _mono_x(qc, c, t)
            elif mode == 'meas':
                qc.h(t); qc.measure(t, creg[km])
                with qc.if_test((creg[km], 1)):
                    qc.x(t)
                    for c in corr:
                        if via_target and c[0] == 'cz':
                            qc.cx(c[1], t); qc.cz(t, c[2]); qc.cx(c[1], t)
                        else:
                            _mono_z(qc, c)
                km += 1
            else:
                raise ValueError('erase mode %r' % (mode,))
        else:
            raise ValueError('op %r' % (op,))
    qc.metadata = dict(erase_modes=modes)
    return qc


def reference(C):
    qc = QuantumCircuit(C.N)
    for op in C.ops:
        g = op[0]
        if g == 'x':
            qc.x(op[1])
        elif g == 'cx':
            qc.cx(op[1], op[2])
        elif g in ('ccx', 'AND'):
            qc.ccx(op[1], op[2], op[3])
        elif g == 'MZ':
            for c in op[2]:
                _mono_x(qc, c, op[1])
        elif g in ('z', 'cz'):
            getattr(qc, g)(*op[1:])
        else:
            raise ValueError('op %r' % (op,))
    return qc


def check_quantum(block, erase='meas', via_target=True, trials=3, shots=8, seed=0, body=None):
    """exact gate-level check of a lowered variant on random product input states (noise-free Aer, random
    measurement outcomes): prepare random single-qubit states on the superposable data wires, run the variant, run
    the inverse of the ideal reference, undo the preparation; the all-zero state must come back with probability 1
    in every shot. -> worst probability seen (1.0 = equivalent on the tested states, phases included)."""
    import numpy as np
    from qiskit_aer import AerSimulator
    C = block.C; rng = np.random.default_rng(seed); worst = 1.0
    sim = AerSimulator(method='statevector', max_parallel_threads=1)
    body = body if body is not None else to_qiskit(C, erase, via_target); ref_inv = reference(C).inverse()
    for _ in range(trials):
        qc = QuantumCircuit(*body.qregs, *body.cregs)
        ang = {w: rng.uniform(0, 2 * np.pi, 3) for w in block.xfree}
        fixed = {}
        d = block.inputs(1, int(rng.integers(1 << 30)))[0]
        for k, ws in block.regs.items():
            for j, w in enumerate(ws):
                if w not in ang:
                    fixed[w] = (d[k] >> j) & 1
        for w, a in ang.items():
            qc.u(a[0], a[1], a[2], w)
        for w, b in fixed.items():
            if b:
                qc.x(w)
        qc.compose(body, inplace=True)
        qc.compose(ref_inv, qubits=range(C.N), inplace=True)
        for w, b in fixed.items():
            if b:
                qc.x(w)
        for w, a in ang.items():
            qc.u(-a[0], -a[2], -a[1], w)
        qc.save_amplitudes_squared([0], label='p', pershot=True)
        r = sim.run(qc, shots=shots, seed_simulator=int(rng.integers(1 << 30))).result()
        ps = [float(x[0]) for x in r.data(0)['p']]
        worst = min(worst, min(ps))
    return worst
