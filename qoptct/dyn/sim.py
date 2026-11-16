# SPDX-License-Identifier: Apache-2.0
"""dyn.sim: noisy trajectory simulation (qiskit-aer) of a scheduled circuit (hw.Sched).

The noise comes from the same HwModel numbers as the estimate (a local Aer noise model on the active qubits):
    after every sx / x and two-qubit gate a depolarizing channel with the calibrated average error;
    every idle period (including the waiting during a mid-circuit measurement, the feed-forward latency and the
    untaken branch of a conditional) a thermal relaxation channel (T1, T2);
    every mid-circuit readout is flipped with the calibrated readout error (so the fix-up is applied wrongly).
Final readout is ideal: the success probability is read from the state (save_amplitudes_squared), averaged over
trajectories.

Two tests of a compiled block:
    Z  basis inputs: probability that all data wires hold the right answer and all other qubits are |0>;
    X  superposition input (|+> on the superposable data wires): after the noisy circuit the noise-free inverse of
       the ideal reference is applied and the preparation undone; the probability of all zeros equals the fidelity of
       the noisy output with the ideal output state (sensitive to phase errors, e.g. a wrong fix-up).
"""
import numpy as np

from .lower import reference


def _compact(S):
    qs = sorted(S.qubits)
    return {q: i for i, q in enumerate(qs)}


def idle_error(t1, t2, dt):
    """Pauli-twirled thermal relaxation for an idle period dt: pX = pY = (1 - exp(-dt/T1))/4,
    pZ = (1 - exp(-dt/T2))/2 - pX. Same average fidelity as the exact channel (1/2 + exp(-dt/T1)/6 + exp(-dt/T2)/3)
    and the exact decay of coherences; the T1 decay is symmetrized (|0> -> |1> as likely as |1> -> |0>).
    (The exact channel needs Kraus operators for T2 > T1, and qiskit-aer 0.17 / 0.18 fails on trajectory
    simulations whose noise model mixes such channels with other errors: "QuantumError: Kraus is empty".)"""
    import math
    from qiskit_aer.noise import pauli_error
    t2 = min(t2, 2 * t1)
    px = (1 - math.exp(-dt / t1)) / 4
    pz = max((1 - math.exp(-dt / t2)) / 2 - px, 0.0)
    return pauli_error([('I', 1 - 2 * px - pz), ('X', px), ('Y', px), ('Z', pz)])


def _emit(qc, items, cm, model, noise, clbits):
    from qiskit.circuit.library import IGate
    for it in items:
        k = it[0]
        if k == 'g':
            _, name, qs, params, err = it
            getattr(qc, name)(*params, *[cm[q] for q in qs])
            if err > 0 and name != 'rz':
                noise[('g', name, tuple(cm[q] for q in qs))] = err
        elif k == 'idle':
            _, q, dt = it
            if q not in cm:
                continue
            lab = 'idle_%d_%d' % (q, round(dt * 1e10))
            noise[('i', lab, cm[q])] = (q, dt)
            qc.append(IGate(label=lab), [cm[q]])
        elif k == 'meas':
            qc.measure(cm[it[1]], clbits[it[2]])
        elif k == 'if':
            _, c, body, els = it
            with qc.if_test((clbits[c], 1)) as else_:
                _emit(qc, body, cm, model, noise, clbits)
            with else_:
                _emit(qc, els, cm, model, noise, clbits)
                if not els:
                    qc.rz(0.0, 0)


def noisy_body(S, model, n_clbits):
    """-> (circuit on compact qubits; idle periods are labelled identity gates, compact map {physical: index},
    Aer NoiseModel: depolarizing errors on the gates, thermal relaxation on the idle labels, readout errors)."""
    from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error
    cm = _compact(S)
    qc = QuantumCircuit(QuantumRegister(len(cm), 'q'))
    clbits = []
    if n_clbits:
        cr = ClassicalRegister(n_clbits, 'm'); qc.add_register(cr); clbits = list(cr)
    noise = {}
    _emit(qc, S.items, cm, model, noise, clbits)
    nm = NoiseModel(basis_gates=['rz', 'sx', 'x', 'id', model.twoq, 'cx', 'ccx', 'h'])
    for key, val in noise.items():
        if key[0] == 'g':
            _, name, qs = key; d = 2 ** len(qs)
            e = depolarizing_error(min(val * d / (d - 1), 1.0), len(qs))
            nm.add_quantum_error(e, name, list(qs))
            if len(qs) == 2 and ('g', name, qs[::-1]) not in noise:
                nm.add_quantum_error(e, name, list(qs[::-1]))
        else:
            _, lab, i = key; q, dt = val
            nm.add_quantum_error(idle_error(model.t1[q], model.t2[q], dt), lab, [i])
    for q, i in cm.items():
        e = model.ro[q]
        if e > 0:
            nm.add_readout_error(ReadoutError([[1 - e, e], [e, 1 - e]]), [i])
    return qc, cm, nm


def _ideal(src):
    """copy of a circuit whose gates carry the label 'ideal' (the noise model does not touch them)."""
    from qiskit import QuantumCircuit
    out = QuantumCircuit(*src.qregs)
    for inst in src.data:
        op = inst.operation.to_mutable(); op.label = 'ideal'
        out.append(op, inst.qubits)
    return out


def _ix(qc, q):
    from qiskit.circuit.library import XGate
    qc.append(XGate(label='ideal'), [q])


def _ih(qc, q):
    from qiskit.circuit.library import HGate
    qc.append(HGate(label='ideal'), [q])


def run_tests(block, comp, model, shots=400, n_inputs=8, seed=0, method=None, tests=('Z', 'X'), x_repeat=None):
    """comp: result of hw.compile_circuit for a lowering of `block`. -> dict(pZ, pZ_err, pX, pX_err, n_active,
    method, per_input)."""
    from qiskit import QuantumCircuit
    from qiskit_aer import AerSimulator
    S = comp['sched']; tqc = comp['tqc']
    body, cm, nm = noisy_body(S, model, tqc.num_clbits)
    n = len(cm); init = comp['layout_init']; fin = comp['layout_final']
    if method is None:
        method = 'statevector' if n <= 16 else 'matrix_product_state'
    for w in block.data:                                    # a data wire the compiler never touched
        for L in (init, fin):
            if L[w] not in cm:
                raise ValueError('data wire %d is not active in the compiled circuit' % w)
    circs = []; kinds = []
    if 'Z' in tests:
        for d in block.inputs(n_inputs, seed):
            qc = QuantumCircuit(*body.qregs, *body.cregs)
            for k, ws in block.regs.items():
                for j, w in enumerate(ws):
                    if (d[k] >> j) & 1:
                        _ix(qc, cm[init[w]])
            qc.compose(body, inplace=True)
            e = block.fn(d)
            for k, ws in block.regs.items():                # map the expected answer to |0...0> (noise-free X)
                for j, w in enumerate(ws):
                    if (e[k] >> j) & 1:
                        _ix(qc, cm[fin[w]])
            qc.save_amplitudes_squared([0], label='p')
            circs.append(qc); kinds.append('Z')
    if 'X' in tests:
        ref_inv = _ideal(reference(block.C).inverse())
        pos = []
        for w in range(block.C.N):
            p = fin[w]
            if p not in cm:
                pos = None; break
            pos.append(cm[p])
        if pos is None:                                     # an untouched ancilla: give the reference a spare slot
            raise ValueError('an ancilla wire is not active in the compiled circuit')
        xf = set(block.xfree)
        for d in block.inputs(max(1, min(2, n_inputs)), seed + 1)[:1 if xf == set(block.data) else 2]:
            qc = QuantumCircuit(*body.qregs, *body.cregs)
            fixed = [(w, (d[k] >> j) & 1) for k, ws in block.regs.items() for j, w in enumerate(ws) if w not in xf]
            for w in xf:
                _ih(qc, cm[init[w]])
            for w, b in fixed:
                if b:
                    _ix(qc, cm[init[w]])
            qc.compose(body, inplace=True)
            qc.compose(ref_inv, qubits=pos, inplace=True)
            for w, b in fixed:
                if b:
                    _ix(qc, pos[w])
            for w in xf:
                _ih(qc, pos[w])
            qc.save_amplitudes_squared([0], label='p')
            for _ in range(x_repeat or max(1, n_inputs // 4)):   # more trajectories for the single X input
                circs.append(qc.copy()); kinds.append('X')
    sim = AerSimulator(method=method, noise_model=nm, max_parallel_threads=1,
                       max_parallel_experiments=1, max_parallel_shots=1)
    res = sim.run(circs, shots=shots, seed_simulator=seed).result()
    out = dict(n_active=n, method=method, shots=shots)
    for kind in ('Z', 'X'):
        ps = [float(res.data(i)['p'][0]) for i, k in enumerate(kinds) if k == kind]
        if ps:
            m = float(np.mean(ps))
            # trajectory noise: p per trajectory is in [0, 1]; a conservative standard error
            out['p' + kind] = m; out['p%s_err' % kind] = float(np.sqrt(max(m * (1 - m), 1e-9) / (shots * len(ps))))
            out['n' + kind] = len(ps)
    return out
