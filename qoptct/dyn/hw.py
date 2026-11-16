# SPDX-License-Identifier: Apache-2.0
"""dyn.hw: hardware cost model for static and dynamic circuits.

HwModel.from_backend(backend) reads a calibration (BackendV2 target): per qubit T1, T2, readout error, measurement
length; per gate error and length. Two numbers are not in the calibration and are parameters:
    ff_latency   classical feed-forward latency between the end of a measurement and the first conditional gate
                 (IBM: about 600 ns);
    p_branch     probability that the conditional body runs (1/2 for an X-basis erasure).

schedule(tqc, model, live) puts a transpiled circuit (if_else allowed) on a time line: every instruction starts when
its qubits (and, for a conditional, its classical bit + ff_latency) are ready; the gaps are explicit idle periods.
Only the qubits of a conditional wait for it (IBM runs independent conditionals in parallel). Both branches of a
conditional take the same time. The result (Sched) carries the time line, the counts and the estimated fidelity

    F_est = prod_gates (1 - e_g) * prod_mid-circuit measurements (1 - e_ro) * prod_idle F_idle(T1, T2, dt),
    F_idle = 1/2 + exp(-dt/T1)/6 + exp(-dt/T2)/3      (average fidelity of thermal relaxation),

where the factors inside a conditional body enter with weight p_branch (and the idle of the other branch with
1 - p_branch). A wrong mid-circuit readout is counted as a failure (the fix-up is applied wrongly). Final readout is
not included (it is the same for every variant).
"""
import copy
import math


class HwModel:
    def __init__(s, name, nq, t1, t2, ro, meas_dur, g1, g2, ff_latency=600e-9, p_branch=0.5, twoq='cz',
                 backend=None):
        s.name = name; s.nq = nq; s.t1 = t1; s.t2 = t2; s.ro = ro; s.meas_dur = meas_dur
        s.g1 = g1                # {(name, q): (err, dur)}
        s.g2 = g2                # {(a, b): (err, dur)}  (directed, both directions present if both calibrated)
        s.ff_latency = ff_latency; s.p_branch = p_branch; s.twoq = twoq; s.backend = backend

    @classmethod
    def from_backend(cls, backend, ff_latency=600e-9, p_branch=0.5):
        t = backend.target; nq = backend.num_qubits
        med = lambda xs, d: (sorted(xs)[len(xs) // 2] if xs else d)
        qp = t.qubit_properties or [None] * nq
        t1r = [getattr(p, 't1', None) for p in qp]; t2r = [getattr(p, 't2', None) for p in qp]
        m1 = med([x for x in t1r if x], 100e-6); m2 = med([x for x in t2r if x], 100e-6)
        t1 = [x or m1 for x in t1r]; t2 = [min(x or m2, 2 * a) for x, a in zip(t2r, t1)]
        ro = [0.01] * nq; md = [1.5e-6] * nq
        for qa, p in t['measure'].items():
            if p is not None:
                if p.error is not None:
                    ro[qa[0]] = min(float(p.error), 0.5)
                if p.duration is not None:
                    md[qa[0]] = float(p.duration)
        g1 = {}
        for nm in ('sx', 'x', 'rz', 'id'):
            if nm in t.operation_names:
                for qa, p in t[nm].items():
                    e = 0.0 if nm == 'rz' else min(float(p.error), 0.75) if (p and p.error is not None) else 3e-4
                    d = 0.0 if nm == 'rz' else float(p.duration) if (p and p.duration is not None) else 3e-8
                    g1[(nm, qa[0])] = (e, d)
        twoq = next(nm for nm in ('cz', 'ecr', 'cx', 'rzz') if nm in t.operation_names)
        g2 = {}
        for qa, p in t[twoq].items():
            e = min(float(p.error), 0.9375) if (p and p.error is not None) else 5e-3
            d = float(p.duration) if (p and p.duration is not None) else 7e-8
            g2[tuple(qa)] = (e, d)
        return cls(getattr(backend, 'name', 'backend'), nq, t1, t2, ro, md, g1, g2, ff_latency, p_branch, twoq,
                   backend)

    def scaled(s, ro=1.0, meas_dur=1.0, coh=1.0, g2=1.0, g1=1.0, ff_latency=None, name=None):
        """a hypothetical device: readout error x ro, measurement length x meas_dur, T1 and T2 x coh, two-qubit
        error x g2, one-qubit error x g1, feed-forward latency = ff_latency seconds (None: unchanged)."""
        m = copy.copy(s)
        m.ro = [min(x * ro, 0.5) for x in s.ro]; m.meas_dur = [x * meas_dur for x in s.meas_dur]
        m.t1 = [x * coh for x in s.t1]; m.t2 = [x * coh for x in s.t2]
        m.g2 = {k: (min(e * g2, 0.9375), d) for k, (e, d) in s.g2.items()}
        m.g1 = {k: (min(e * g1, 0.75), d) for k, (e, d) in s.g1.items()}
        if ff_latency is not None:
            m.ff_latency = ff_latency
        m.name = name or '%s[ro x%g, meas x%g, T x%g, 2q x%g, ff %g ns]' % (s.name, ro, meas_dur, coh, g2,
                                                                        m.ff_latency * 1e9)
        return m

    def idle_fid(s, q, dt):
        if dt <= 0:
            return 1.0
        return 0.5 + math.exp(-dt / s.t1[q]) / 6 + math.exp(-dt / s.t2[q]) / 3

    def summary(s, qubits=None):
        qs = list(qubits) if qubits is not None else list(range(s.nq))
        med = lambda xs: sorted(xs)[len(xs) // 2]
        e2 = [e for (a, b), (e, d) in s.g2.items() if a in qs and b in qs] or [e for e, d in s.g2.values()]
        d2 = [d for (a, b), (e, d) in s.g2.items() if a in qs and b in qs] or [d for e, d in s.g2.values()]
        return dict(name=s.name, t1_us=med([s.t1[q] for q in qs]) * 1e6, t2_us=med([s.t2[q] for q in qs]) * 1e6,
                    ro=med([s.ro[q] for q in qs]), meas_ns=med([s.meas_dur[q] for q in qs]) * 1e9,
                    e2q=med(e2), t2q_ns=med(d2) * 1e9, ff_ns=s.ff_latency * 1e9)


class Sched:
    """time line of a compiled circuit. items: ('g', name, qubits, params, err) | ('idle', q, dt) |
    ('meas', q, clbit) | ('if', clbit, body_items, else_items) | ('reset', q)."""

    def __init__(s):
        s.items = []; s.duration = 0.0; s.logF = 0.0; s.n2q = 0; s.n2q_cond = 0; s.n1q = 0; s.n_meas = 0
        s.n_if = 0; s.idle_time = 0.0; s.qubits = set(); s.parts = dict(g2=0.0, g1=0.0, ro=0.0, idle=0.0)

    @property
    def fid(s):
        return math.exp(s.logF)

    @property
    def n2q_exp(s):
        return s.n2q + s.n2q_cond * 0.5

    def stats(s):
        return dict(n2q=s.n2q, n2q_cond=s.n2q_cond, n1q=s.n1q, n_meas=s.n_meas, duration_us=s.duration * 1e6,
                    idle_us=s.idle_time * 1e6, qubits=len(s.qubits), f_est=s.fid,
                    loss_2q=1 - math.exp(s.parts['g2']), loss_1q=1 - math.exp(s.parts['g1']),
                    loss_ro=1 - math.exp(s.parts['ro']), loss_idle=1 - math.exp(s.parts['idle']))


def _log1m(e):
    return math.log(max(1.0 - e, 1e-12))


def schedule(tqc, model, live=(), final=()):
    """tqc: circuit on physical qubits in the backend basis (if_else allowed, no other control flow).
    live: physical qubits holding data from t = 0 (they idle until first used); other qubits are |0> until their
    first gate. final: physical qubits holding data at the end (padded to the end of the circuit)."""
    S = Sched()
    qi = {q: i for i, q in enumerate(tqc.qubits)}; ci = {c: i for i, c in enumerate(tqc.clbits)}
    tq = {q: 0.0 for q in live}; tc = {}

    def add_idle(items, q, dt, w):
        if dt > 1e-12:
            items.append(('idle', q, dt)); S.idle_time += dt * w
            lf = w * math.log(model.idle_fid(q, dt)); S.logF += lf; S.parts['idle'] += lf

    def gate(items, name, qs, params, times, w, cond):
        if name == 'barrier':
            return
        if name == 'delay':
            raise ValueError('explicit delays are not supported; schedule() inserts the idle periods itself')
        if len(qs) == 1:
            e, d = model.g1.get((name, qs[0]), (0.0, 0.0) if name == 'rz' else (3e-4, 3e-8))
        else:
            if tuple(qs) in model.g2:
                e, d = model.g2[tuple(qs)]
            elif tuple(qs[::-1]) in model.g2:
                e, d = model.g2[tuple(qs[::-1])]
            else:
                raise ValueError('two-qubit gate on uncoupled qubits %r' % (qs,))
        present = [times[q] for q in qs if q in times]
        start = max(present) if present else 0.0
        if not present and not cond:
            start = 0.0
        for q in qs:
            if q in times:
                add_idle(items, q, start - times[q], w)
            times[q] = start + d
            S.qubits.add(q)
        items.append(('g', name, tuple(qs), tuple(float(p) for p in params), e))
        lf = w * _log1m(e); S.logF += lf
        if len(qs) == 2:
            S.parts['g2'] += lf
            if cond:
                S.n2q_cond += 1
            else:
                S.n2q += 1
        elif name != 'rz':
            S.parts['g1'] += lf; S.n1q += 1

    for inst in tqc.data:
        op = inst.operation; name = op.name; qs = [qi[q] for q in inst.qubits]
        if name == 'measure':
            q = qs[0]; c = ci[inst.clbits[0]]
            start = tq.get(q, 0.0); tq[q] = start + model.meas_dur[q]; tc[c] = tq[q]; S.qubits.add(q)
            S.items.append(('meas', q, c)); S.n_meas += 1
            lf = _log1m(model.ro[q]); S.logF += lf; S.parts['ro'] += lf
        elif name == 'reset':
            raise ValueError('reset is not supported')
        elif name == 'if_else':
            cond = op.condition
            if not (isinstance(cond, tuple) and cond[0] in ci and cond[1] == 1) or op.params[1] is not None:
                raise ValueError('only `if (bit == 1) {...}` without else is supported')
            c = ci[cond[0]]; body = op.params[0]; S.n_if += 1
            bq = {b: qs[i] for i, b in enumerate(body.qubits)}
            present = [tq[q] for q in qs if q in tq]
            start = max(present + [tc.get(c, 0.0) + model.ff_latency])
            for q in qs:
                if q in tq:
                    add_idle(S.items, q, start - tq[q], 1.0)
                tq[q] = start; S.qubits.add(q)
            loc = {q: 0.0 for q in qs}; items = []; w = model.p_branch
            for bi in body.data:
                if bi.operation.name in ('if_else', 'measure', 'reset'):
                    raise ValueError('nested control flow / measurement inside a conditional is not supported')
                gate(items, bi.operation.name, [bq[b] for b in bi.qubits], bi.operation.params, loc, w, True)
            D = max(loc.values()) if loc else 0.0
            els = []
            for q in qs:
                add_idle(items, q, D - loc[q], w)
                add_idle(els, q, D, 1.0 - w)
                tq[q] = start + D
            S.items.append(('if', c, items, els))
        else:
            gate(S.items, name, qs, op.params, tq, 1.0, False)
    T = max(tq.values()) if tq else 0.0
    for q in final:
        if q in tq:
            add_idle(S.items, q, T - tq[q], 1.0)
    S.duration = T
    return S


def compile_circuit(qc, model, seeds=(0, 1, 2, 3), level=3, live_virtual=(), initial_layouts=()):
    """transpile qc for model.backend with several seeds (and optional initial layouts), schedule each, return the
    best by estimated fidelity: dict(tqc, sched, layout_init, layout_final, seed)."""
    from qiskit import transpile
    best = None
    trials = [(s, None) for s in seeds] + [(seeds[0], L) for L in initial_layouts if L is not None]
    for seed, L in trials:
        try:
            tqc = transpile(qc, model.backend, optimization_level=level, seed_transpiler=seed, initial_layout=L)
        except Exception as e:                              # a layout that does not fit, etc.
            err = e; continue
        if tqc.layout is not None:
            init = tqc.layout.initial_index_layout(filter_ancillas=True)
            fin = tqc.layout.final_index_layout(filter_ancillas=True)
        else:
            init = fin = list(range(qc.num_qubits))
        try:
            S = schedule(tqc, model, live=[init[w] for w in live_virtual], final=[fin[w] for w in live_virtual])
        except ValueError as e:
            err = e; continue
        if best is None or S.logF > best['sched'].logF:
            best = dict(tqc=tqc, sched=S, layout_init=init, layout_final=fin, seed=seed, from_layout=L is not None)
    if best is None:
        raise RuntimeError('compilation failed: %s' % err)
    return best
