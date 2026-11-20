# SPDX-License-Identifier: Apache-2.0
"""dyn.choose: pick between measurement-based erasure and reverse computation with the hardware model.

    variants(family, size)                 the candidate lowerings of a block:
        U   every erasure by reverse computation (static circuit)
        M   every erasure by measurement + feed-forward (dynamic circuit); the cz fix-up through the measured wire
        Md  the same with the direct cz(a, b) fix-up
        H   hybrid: single ANDs by reverse computation, multi-term erasures (lookup outputs) by measurement
        T   Toffoli baseline without logical ANDs (Cuccaro adder)
    choose(family, size, model)            compile all, estimate, (optionally simulate), return the best + the table
    site_rule(model, n_live)               the analytic per-site rule (no compilation)
    DynErasePass / choose_circuit          the same decision for a Qiskit circuit that already contains
                                           `h t; measure t -> c; if (c) {cz a, b; x t}` erasures

Dynamic variants (M, Md, H) are selected only with allow_dynamic=True (CLI: --allow-dynamic); see DYN_WARNING.
"""
import math
import warnings

from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.transpiler.basepasses import TransformationPass

from . import gen, hw, lower


DYN_WARNING = ('qoptct.dyn: allow_dynamic is on. On a real chip (ibm_kingston, 20.11.2026, one job per experiment) the '
               'model strongly underestimated the cost of a mid-circuit measurement (adder2: predicted 0.918, '
               'measured 0.525) and the hybrid lost to ordinary uncomputation (0.682 against 0.772). Dynamic '
               'variants (M, H) are not selected by default; the model estimate for them is optimistic.')


def warn_dynamic():
    """the warning printed whenever dynamic (measurement-based) variants may be selected."""
    warnings.warn(DYN_WARNING, UserWarning, stacklevel=3)


def resolve_backend(name):
    """'FakeBoston' (qiskit_ibm_runtime.fake_provider) or a backend object."""
    if not isinstance(name, str):
        return name
    import qiskit_ibm_runtime.fake_provider as fp
    if not hasattr(fp, name):
        raise ValueError('unknown fake backend %r' % name)
    return getattr(fp, name)()


def _is_and_site(op):
    return len(op) == 3 and len(op[2]) == 1 and op[2][0][0] == 'cz'


def variants(family, size, **kw):
    """-> list of (label, block, erase, via_target)."""
    out = []
    for B in gen.build(family, size, **kw):
        ns = lower.n_sites(B.C)
        if B.variant == 'toffoli':
            out.append(('T', B, 'unitary', True)); continue
        if ns == 0:
            out.append(('U', B, 'unitary', True)); continue
        if B.variant == 'out-meas':
            out.append(('M', B, 'meas', True)); out.append(('Md', B, 'meas', False))
            out.append(('H', B, lambda i, op: 'unitary' if _is_and_site(op) else 'meas', True))
        elif B.variant == 'out-unitary':
            out.append(('U', B, 'unitary', True))
        else:
            out.append(('U', B, 'unitary', True)); out.append(('M', B, 'meas', True))
            out.append(('Md', B, 'meas', False))
    return out


def evaluate(block, erase, via_target, model, seeds=(0, 1, 2, 3), layouts=(), router=None):
    qc = lower.to_qiskit(block.C, erase, via_target)
    comp = hw.compile_circuit(qc, model, seeds=seeds, live_virtual=block.data, initial_layouts=layouts,
                              router=router)
    comp['logical'] = qc
    return comp


def choose(family, size, model, seeds=(0, 1, 2, 3), simulate=False, shots=300, n_inputs=8, sim_max_qubits=24,
           check=True, router=None, allow_dynamic=False, **kw):
    """-> dict(best=label, rows=[...], comps={label: comp}, blocks={label: block}). Each row: label, qubits, 2q
    counts, duration, estimate (+ simulated pZ / pX). router: see hw.compile_circuit (e.g. hw.hw_router()).
    Default: 'best' is the best STATIC variant (U or T); the dynamic ones (M, Md, H) are still compiled and listed in
    the table ('best_model' is the best of all by the model estimate), but are 'best' only with allow_dynamic=True
    (which warns: on ibm_kingston the model was too optimistic for them)."""
    from . import sim
    if allow_dynamic:
        warn_dynamic()
    rows = []; comps = {}; blocks = {}; layout_u = None
    for label, B, erase, vt in variants(family, size, **kw):
        if check:
            ok, why, nin, exh = B.check(4096)
            if not ok:
                raise RuntimeError('%s %s: classical check failed: %s' % (B.name, label, why))
        comp = evaluate(B, erase, vt, model, seeds, [layout_u] if layout_u is not None else [], router)
        if label == 'U':
            layout_u = comp['layout_init']
        S = comp['sched']; st = S.stats()
        row = dict(circuit=B.name, variant=label, logical_qubits=B.n, sites=lower.n_sites(B.C),
                   depth=comp['tqc'].depth(), routed=comp['routed'], **st)
        if simulate and len(S.qubits) <= sim_max_qubits:
            r = sim.run_tests(B, comp, model, shots=shots, n_inputs=n_inputs)
            row.update({k: r[k] for k in ('pZ', 'pZ_err', 'pX', 'pX_err', 'methodZ', 'methodX') if k in r})
        rows.append(row); comps[label] = comp; blocks[label] = B
    def top(rs):
        return max(rs, key=lambda r: (round(r['f_est'], 9), r['variant'] in ('U', 'T')))['variant']   # ties: static
    best_model = top(rows)
    best = best_model if allow_dynamic else top([r for r in rows if r['variant'] in ('U', 'T')])
    return dict(best=best, best_model=best_model, allow_dynamic=allow_dynamic, rows=rows, comps=comps, blocks=blocks, model=model.summary())


def rescore(result, model):
    """re-estimate the compiled variants of a choose() result under another model (same compiled circuits).
    -> {label: f_est}."""
    out = {}
    for label, comp in result['comps'].items():
        B = result['blocks'][label]
        S = hw.schedule(comp['tqc'], model, live=[comp['layout_init'][w] for w in B.data],
                        final=[comp['layout_final'][w] for w in B.data])
        out[label] = S.fid
    return out


def _best_split(result, M2):
    sc = rescore(result, M2)
    st = max(v for k, v in sc.items() if k in ('U', 'T'))
    dy = max([v for k, v in sc.items() if k not in ('U', 'T')] or [0.0])
    return st, dy


def breakeven(result, model, coh=1.0, g2=1.0, ro=1.0, hi_ns=None):
    """measurement chain (measurement length + feed-forward latency, ns) at which the best dynamic variant of a
    choose() result gets as good as the best static one, with T1/T2 x coh, two-qubit error x g2, readout error x ro.
    The chain is scaled as a whole (length and latency by the same factor). Bisection (the dynamic estimate grows
    monotonically when the chain gets shorter). -> ns, or 0.0 if even a zero-length chain does not win, or None if
    the dynamic variant already wins at the current chain length."""
    s = model.summary(); chain = s['meas_ns'] + s['ff_ns']

    def gap(t):
        f = t / chain
        M2 = model.scaled(ro=ro, coh=coh, g2=g2, meas_dur=f, ff_latency=model.ff_latency * f)
        st, dy = _best_split(result, M2)
        return dy - st
    if gap(chain) >= 0:
        return None
    if gap(1e-3) < 0:
        return 0.0
    lo, hi = 1e-3, chain
    for _ in range(40):
        mid = (lo + hi) / 2
        if gap(mid) >= 0:
            lo = mid
        else:
            hi = mid
    return lo


def site_rule(model, n_live, corr_2q=3, unit_2q=3, qubits=None):
    """analytic cost (infidelity) of erasing one AND, without compilation:
        measurement: e_ro + n_live * rate * (t_meas + t_ff) + p_branch * (corr_2q * e_2q + n_live * rate * corr_2q * t_2q)
        reverse    : unit_2q * e_2q + n_live * rate * unit_2q * t_2q
    rate = 1/(6 T1) + 1/(3 T2) per idle qubit; n_live = qubits that wait for the erasure (in a ripple circuit: all
    live qubits). -> dict(meas, unitary, better)."""
    s = model.summary(qubits)
    rate = 1 / (6 * s['t1_us'] * 1e-6) + 1 / (3 * s['t2_us'] * 1e-6)
    tm = s['meas_ns'] * 1e-9 + model.ff_latency; t2 = s['t2q_ns'] * 1e-9; e2 = s['e2q']
    cm = s['ro'] + n_live * rate * tm + model.p_branch * (corr_2q * e2 + n_live * rate * corr_2q * t2)
    cu = unit_2q * e2 + n_live * rate * unit_2q * t2
    return dict(meas=cm, unitary=cu, better='meas' if cm < cu else 'unitary')


# ------------------------------------------------------------------ Qiskit circuits
def _match_sites(qc):
    """find `h t; measure t -> c; if (c == 1) {z / cz / ccz fix-ups not touching t; x t}`.
    -> list of (index_h, index_measure, index_if, t, terms)."""
    data = qc.data; last = {}; sites = []; pending = {}
    for i, inst in enumerate(data):
        name = inst.operation.name
        if name == 'measure':
            t = inst.qubits[0]; j = last.get(t)
            if j is not None and data[j].operation.name == 'h':
                pending[inst.clbits[0]] = (j, i, t)
        elif name == 'if_else':
            cond = inst.operation.condition
            if isinstance(cond, tuple) and cond[0] in pending and cond[1] == 1 and inst.operation.params[1] is None:
                j, m, t = pending.pop(cond[0]); body = inst.operation.params[0]
                bq = {b: inst.qubits[k] for k, b in enumerate(body.qubits)}
                terms = []; okb = last.get(t) == m; nx = 0
                for bi in body.data:
                    qs = [bq[b] for b in bi.qubits]; bn = bi.operation.name
                    if bn == 'x' and qs == [t]:
                        nx += 1
                    elif bn in ('z', 'cz', 'ccz') and t not in qs:
                        terms.append(tuple(qs))
                    else:
                        okb = False
                if okb and nx == 1:
                    sites.append((j, m, i, t, terms))
        for q in inst.qubits:
            last[q] = i
        if name != 'if_else':
            for c in inst.clbits:
                if name != 'measure':
                    pending.pop(c, None)
    return sites


def to_unitary_erasure(qc, hybrid=False):
    """replace every matched measurement-based erasure by reverse computation. Precondition (the same one that
    makes the measurement-based erasure itself correct): the erased qubit holds exactly the XOR of the fix-up
    monomials (t = a AND b for a single cz(a, b)). -> (new circuit, number of sites replaced)."""
    sites = _match_sites(qc)
    if hybrid:                                              # only the single-AND sites; multi-term ones stay
        sites = [x for x in sites if len(x[4]) == 1 and len(x[4][0]) == 2]
    if not sites:
        return qc, 0
    drop = {}; at = {}
    for j, m, i, t, terms in sites:
        drop[j] = drop[m] = True; at[i] = (t, terms)
    out = qc.copy_empty_like()
    for i, inst in enumerate(qc.data):
        if i in at:
            t, terms = at[i]
            if len(terms) == 1 and len(terms[0]) == 2:
                lower._and_inv(out, terms[0][0], terms[0][1], t)
            else:
                for ws in terms:
                    lower._mono_x(out, ('m',) + tuple(ws), t)
        elif i not in drop:
            out.append(inst)
    return out, len(sites)


def choose_circuit(qc, model, seeds=(0, 1, 2, 3), mode='auto', allow_dynamic=False):
    """qc: a logical circuit with measurement-based erasures. mode: 'auto' (compile both forms, keep the one with
    the higher estimated fidelity), 'rule' (analytic site rule, no compilation), 'unitary', 'keep'.
    Without allow_dynamic=True the result is always the static (reverse-computation) form, whatever the mode says
    ('keep' is the exception: the caller asks explicitly to leave the circuit as it is). With allow_dynamic=True the
    measurement-based and hybrid forms may be chosen by the model, with a warning (see DYN_WARNING).
    -> (circuit, report)."""
    alt, n = to_unitary_erasure(qc)
    rep = dict(sites=n, mode=mode, allow_dynamic=allow_dynamic)
    if n == 0 or mode == 'keep':
        rep['chosen'] = 'keep'; return qc, rep
    if mode == 'unitary' or not allow_dynamic:
        rep['chosen'] = 'unitary'; return alt, rep
    warn_dynamic()
    if mode == 'rule':
        r = site_rule(model, qc.num_qubits); rep.update(rule=r)
        rep['chosen'] = 'keep' if r['better'] == 'meas' else 'unitary'
        return (qc if r['better'] == 'meas' else alt), rep
    live = list(range(qc.num_qubits))
    a = hw.compile_circuit(qc, model, seeds=seeds, live_virtual=live)
    b = hw.compile_circuit(alt, model, seeds=seeds, live_virtual=live, initial_layouts=[])
    rep.update(f_est_meas=a['sched'].fid, f_est_unitary=b['sched'].fid, stats_meas=a['sched'].stats(),
               stats_unitary=b['sched'].stats())
    cands = {'keep': (a['sched'].fid, qc), 'unitary': (b['sched'].fid, alt)}
    hyb, nh = to_unitary_erasure(qc, hybrid=True)           # hybrid: single ANDs reversed, multi-term erasures measured
    if 0 < nh < n:
        c = hw.compile_circuit(hyb, model, seeds=seeds, live_virtual=live, initial_layouts=[])
        rep['f_est_hybrid'] = c['sched'].fid; cands['hybrid'] = (c['sched'].fid, hyb)
    # the ordinary (static) uncomputation is the default: another form is taken only when its estimate is strictly better
    order = ('unitary', 'hybrid', 'keep')
    best = max((k for k in order if k in cands), key=lambda k: (round(cands[k][0], 9), -order.index(k)))
    rep['chosen'] = best
    return cands[best][1], rep


class DynErasePass(TransformationPass):
    """Hardware-aware choice between measurement-based erasure and reverse computation.

    Run it on the logical circuit (before layout). ``backend`` (or ``model``: a dyn.hw.HwModel) supplies the
    calibration. By default (allow_dynamic=False) the result is always the ordinary (static) uncomputation: on a real
    chip (ibm_kingston, 20.11.2026) the model was far too optimistic for measurement-based erasure and the hybrid
    lost to the static form. With allow_dynamic=True (a warning is issued) mode='auto' compiles the forms (static,
    hybrid, measurement-based) and takes another form than the static one only when its estimate is better; 'rule'
    uses the analytic site rule; 'unitary' always converts. ``self.last`` holds the report."""

    def __init__(self, backend=None, model=None, mode='auto', seeds=(0, 1, 2, 3), ff_latency=600e-9,
                 allow_dynamic=False):
        super().__init__()
        if model is None and backend is not None:
            model = hw.HwModel.from_backend(backend, ff_latency=ff_latency)
        if model is None and mode == 'auto':
            mode = 'unitary'                                 # no calibration: the ordinary (static) uncomputation
        if model is None and mode == 'rule':
            raise ValueError('mode=%r needs a backend or a model' % mode)
        self.model = model; self.mode = mode; self.seeds = tuple(seeds); self.last = None
        self.allow_dynamic = allow_dynamic

    def run(self, dag):
        qc = dag_to_circuit(dag)
        out, self.last = choose_circuit(qc, self.model, self.seeds, self.mode, self.allow_dynamic)
        return dag if out is qc else circuit_to_dag(out)
