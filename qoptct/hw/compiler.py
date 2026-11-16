"""``compile_hw``: compile a logical circuit for a noisy IBM-like device and pick the best verified candidate.

Candidates (a portfolio, as everywhere in qoptct; every candidate that may be returned is verified):

* ``qiskit``       -- Qiskit ``transpile(optimization_level=3)`` (the reference; one run per seed),
* ``qopt+qiskit``  -- qoptct's proven CX/T optimization of the Clifford+T parts first, then Qiskit O3,
* ``hw``           -- the qoptct router (commutation-aware, SWAP absorption, calibration-weighted; router.py)
                      with layout search, then block consolidation / native synthesis / Qiskit's O3 optimization
                      stage and a calibration-aware relocation of the used region,
* ``qopt+hw``      -- the same router on the qoptct-optimized logical circuit.

The winner is the candidate with the highest expected fidelity from the calibration data (metrics.py); with
``objective="twoq"`` or ``"depth"`` the smallest two-qubit count / depth instead.
"""
import time
import random
from dataclasses import dataclass, field

from qiskit import QuantumCircuit, transpile

from .device import Device
from .metrics import analyze
from .check import verify_compiled, maps_from_transpiled
from . import router as R

_PM = {}


def _stages(dev):
    if id(dev) not in _PM:
        from qiskit.transpiler import generate_preset_pass_manager, PassManager
        from qiskit.transpiler.passes import (ConsolidateBlocks, UnitarySynthesis, Optimize1qGatesDecomposition,
                                              CommutativeCancellation, RemoveIdentityEquivalent, GateDirection)
        t = dev.target
        pm = generate_preset_pass_manager(3, dev.backend, seed_transpiler=0)
        pre = PassManager([ConsolidateBlocks(target=t, force_consolidate=True), UnitarySynthesis(target=t)])
        loop = PassManager([ConsolidateBlocks(target=t), UnitarySynthesis(target=t), RemoveIdentityEquivalent(target=t),
                            Optimize1qGatesDecomposition(target=t), CommutativeCancellation(target=t),
                            GateDirection(None, target=t), Optimize1qGatesDecomposition(target=t)])
        _PM[id(dev)] = (pre, pm.translation, loop)
    return _PM[id(dev)]


def _opt_loop(circ, loop):
    """Qiskit's level-3 optimization loop (same passes, run to a fixed point of (size, depth)); no layout change."""
    prev = None
    for _ in range(10):
        circ = loop.run(circ)
        cur = (circ.size(), circ.depth())
        if cur == prev:
            break
        prev = cur
    return circ


def flatten(qc):
    """(circuit of one-/two-qubit gates without SWAPs, wire): higher gates are expanded through their definitions and
    SWAP gates are removed by relabelling; the state of input wire i ends on wire ``wire[i]`` of the result."""
    n = qc.num_qubits
    out = QuantumCircuit(n, global_phase=qc.global_phase)
    wire = list(range(n))

    # SWAP elision: a gate addressed to original wire w acts on the wire that currently holds w's state
    def rec_top(circ):
        for ins in circ.data:
            op = ins.operation
            if op.name in ("barrier", "delay"):
                continue
            if ins.clbits or op.name in ("measure", "reset") or getattr(op, "condition", None) is not None:
                raise ValueError("compile_hw needs a unitary circuit (no measurements / classical control); "
                                 "remove final measurements first (QuantumCircuit.remove_final_measurements)")
            qs = [wire[circ.find_bit(q).index] for q in ins.qubits]
            if op.name == "swap":
                a, b = (circ.find_bit(q).index for q in ins.qubits)
                wire[a], wire[b] = wire[b], wire[a]
                continue
            sub = QuantumCircuit(op.num_qubits)
            sub.append(op, range(op.num_qubits))
            _expand(sub, qs, out)

    rec_top(qc)
    return out, wire


def _expand(circ, qmap, out):
    for ins in circ.data:
        op = ins.operation
        if op.name in ("barrier", "delay"):
            continue
        qs = [qmap[circ.find_bit(q).index] for q in ins.qubits]
        if op.num_qubits <= 2:
            try:
                op.to_matrix()
                out.append(op, qs)
                continue
            except Exception:
                pass
        d = op.definition
        if d is None:
            raise ValueError(f"cannot expand gate {op.name}")
        out.global_phase += d.global_phase
        _expand(d, qs, out)


def cleanup(qc):
    """Light logical clean-up that keeps the gate structure visible to the router (inverse pairs, commuting
    cancellations)."""
    from qiskit.transpiler import PassManager
    from qiskit.transpiler.passes import InverseCancellation, CommutativeCancellation
    from qiskit.circuit.library import CXGate, CZGate, HGate, XGate, ZGate, SGate, SdgGate, TGate, TdgGate
    pm = PassManager([InverseCancellation([CXGate(), CZGate(), HGate(), XGate(), ZGate(), (SGate(), SdgGate()),
                                           (TGate(), TdgGate())]), CommutativeCancellation()])
    try:
        return pm.run(pm.run(qc))
    except Exception:
        return qc


def qopt_logical(qc, seconds=60.0, cxpass=False):
    """qoptct's proven CX/T optimization of the Clifford+T parts of a logical circuit (QoptPass), limited to about
    ``seconds`` (the searches stop there and return the best proven result so far; a hard stop follows at twice the
    time).  Returns (circuit, info); the circuit is the input itself when nothing was improved."""
    from qiskit.transpiler import PassManager
    import qoptct
    from qoptct import qopt_core
    from qoptct.qiskit_qopt import QoptPass
    p = QoptPass(cxpass=cxpass)
    old = qopt_core._END
    qopt_core._END = time.time() + seconds
    try:
        out = qoptct.limited(2 * seconds, lambda: PassManager([p]).run(qc))
    except qopt_core._TO:
        return qc, {"skipped": f"time limit ({seconds:g} s)"}
    finally:
        qopt_core._END = old
    return out, (p.last or {})


def finish(circ, dev):
    """Routed circuit with SWAPs (physical qubits) -> native gates: block consolidation, exact two-qubit synthesis,
    translation and Qiskit's level-3 optimization loop.  The qubit mapping is not changed."""
    pre, tr, opt = _stages(dev)
    return _opt_loop(tr.run(pre.run(circ)), opt)


def valid_for(circ, dev):
    t = dev.target
    for ins in circ.data:
        name = ins.operation.name
        if name in ("barrier", "delay"):
            continue
        if not t.instruction_supported(name, tuple(circ.find_bit(q).index for q in ins.qubits)):
            return False
    return True


def relocate(circ, init_map, final_map, dev, max_maps=3000, seconds=5.0):
    """Move the compiled circuit to the isomorphic region of the device with the best expected fidelity
    (same idea as Qiskit's VF2PostLayout, scored with metrics.analyze terms).  Returns (circ, init, final)."""
    import math
    import rustworkx as rx
    use2, use1 = {}, {}
    for ins in circ.data:
        qs = [circ.find_bit(q).index for q in ins.qubits]
        if len(qs) == 2:
            k = (min(qs), max(qs))
            use2[k] = use2.get(k, 0) + 1
        elif ins.operation.name not in ("rz", "barrier", "delay"):
            use1[qs[0]] = use1.get(qs[0], 0) + 1
    nodes = sorted({q for e in use2 for q in e} | set(use1) | set(init_map) | set(final_map))
    if not use2:
        return circ, init_map, final_map
    idx = {p: i for i, p in enumerate(nodes)}
    g = rx.PyGraph()
    g.add_nodes_from(nodes)
    g.add_edges_from_no_data([(idx[a], idx[b]) for a, b in use2])
    fin = set(final_map)
    dir_needed = dev.twoq == "ecr"

    def score(m):  # m: used physical -> new physical
        s = 0.0
        for (a, b), c in use2.items():
            k = (min(m[a], m[b]), max(m[a], m[b]))
            if k not in dev.edges:
                return math.inf
            s += c * dev.elog[k]
        for a, c in use1.items():
            s += c * -math.log(max(1e-9, 1 - dev.err1[m[a]]))
        for a in fin:
            s += -math.log(max(1e-9, 1 - dev.ro_err[m[a]]))
        return s

    ident = {p: p for p in nodes}
    best, bs = ident, score(ident)
    t0 = time.time()
    good = rx.PyGraph()
    good.add_nodes_from(range(dev.n))
    good.add_edges_from_no_data(sorted(dev.edges))
    cnt = 0
    for mp in rx.vf2_mapping(good, g, subgraph=True, induced=False, id_order=False):
        m = {nodes[j]: i for i, j in mp.items()}
        s = score(m)
        if s < bs - 1e-12:
            best, bs = m, s
        cnt += 1
        if cnt >= max_maps or time.time() - t0 > seconds:
            break
    if best is ident:
        return circ, init_map, final_map
    out = QuantumCircuit(dev.n, global_phase=circ.global_phase)
    for ins in circ.data:
        out.append(ins.operation, [best[circ.find_bit(q).index] for q in ins.qubits])
    if dir_needed and not valid_for(out, dev):
        _, tr, opt = _stages(dev)
        out = _opt_loop(tr.run(out), opt)
    return out, [best[p] for p in init_map], [best[p] for p in final_map]


@dataclass
class Cand:
    name: str
    circuit: QuantumCircuit
    init_map: list
    final_map: list
    metrics: dict
    seconds: float = 0.0
    verified: object = None
    note: str = ""


@dataclass
class HwResult:
    circuit: QuantumCircuit
    init_map: list
    final_map: list
    metrics: dict
    strategy: str
    verified: object
    candidates: list = field(default_factory=list)
    best: dict = field(default_factory=dict)   # best candidate per strategy
    logical: QuantumCircuit = None

    @property
    def table(self):
        rows = []
        for c in self.candidates:
            m = c.metrics
            rows.append({"strategy": c.name, "twoq": m["twoq"], "depth": m["depth"], "depth2q": m["depth2q"],
                         "fid": m["fid"], "esp": m["esp"], "duration_us": m["duration_us"],
                         "verified": c.verified, "seconds": round(c.seconds, 2), "note": c.note})
        return rows


def qiskit_o3(qc, dev, seed=0, name="qiskit"):
    t0 = time.time()
    tc = transpile(qc, dev.backend, optimization_level=3, seed_transpiler=seed)
    init, fin = maps_from_transpiled(tc, qc.num_qubits)
    return Cand(name, tc, init, fin, analyze(tc, dev, fin), time.time() - t0, note=f"seed={seed}")


def _key(objective):
    if objective == "twoq":
        return lambda c: (c.metrics["twoq"], c.metrics["depth2q"], -c.metrics["fid"])
    if objective == "depth":
        return lambda c: (c.metrics["depth2q"], c.metrics["twoq"], -c.metrics["fid"])
    return lambda c: (-c.metrics["fid"], c.metrics["twoq"])


VARIANTS = {
    # name: (commutation-aware dependencies, Params keywords)
    "absorb": (True, dict(mu=1.0, floating=False, tau=None)),
    "float": (True, dict(mu=1.0, floating=True, tau=None)),
    "mirror": (True, dict(mu=1.0, floating=True, tau=0.0)),
    "mirror+": (True, dict(mu=1.0, floating=True, tau=1.5)),
    "mirror/2": (True, dict(mu=0.5, floating=True, tau=0.0)),
    "strict": (False, dict(mu=1.0, floating=False, tau=None)),
    "sabre": (False, dict(mu=0.0, floating=False, tau=None)),
}
LAYOUTS = (("greedy", 2), ("greedy", 1), ("line", 0), ("chain", 0), ("chain", 1))


def route_hw(qc, dev, layouts=(), seconds=20.0, trials=None, seed=0, top=4, name="hw", variants=None, noise=True,
             do_relocate=True, objective="fid", log=None):
    """Candidates of the qoptct router for the logical circuit ``qc`` (any gates; flattened here).

    The search is a small bandit over (router variant, initial-layout kind): every arm is tried once, then most
    trials go to the arms with the best results so far.  Layout kinds: ``greedy`` (interaction-graph placement on a
    good region, refined by forward/backward routing passes), ``line`` (logical qubits in index order on a low-error
    path of the device -- the natural layout of QFT-like and all-to-all circuits), ``chain`` (order from the
    interaction graph on such a path), ``given`` (the layouts passed in, e.g. those found by Qiskit).
    The ``top`` routed circuits (router's own cost estimate) are finished to native gates and scored."""
    t_start = time.time()
    flat, wire = flatten(qc)
    flat = cleanup(flat)
    nl = flat.num_qubits
    sets = {}

    def nodeset(commute):
        if commute not in sets:
            nodes = R.to_nodes(flat, commute=commute)
            rn = R.reverse_nodes(nodes)
            sets[commute] = (nodes, R.build_dag(nodes, nl), rn, R.build_dag(rn, nl))
        return sets[commute]

    nodes = nodeset(True)[0]
    if not any(nd.two for nd in nodes):
        return []
    adj = R.interaction(nodes, nl)
    rng = random.Random(seed)
    paths = [pth for pth in (R.find_path(dev, nl, rng, tries=30) for _ in range(4)) if pth]
    lay_list = [list(l) for l in layouts]
    vnames = list(variants or VARIANTS)
    arms = []
    for v in vnames:
        for kind, ref in LAYOUTS:
            if kind in ("line", "chain") and not paths:
                continue
            arms.append([v, kind, ref, None, 0])     # variant, layout kind, refinement passes, best cost, pulls
        for ref in (0, 1):
            if lay_list:
                arms.append([v, "given", ref, None, 0])
    routed = []
    n_tr = 0
    given_i = 0
    while True:
        if n_tr < len(arms):
            arm = arms[n_tr]
        else:
            ranked = sorted(arms, key=lambda x: (x[3] is None, x[3]))
            r = rng.random()
            arm = ranked[0] if r < 0.45 else rng.choice(ranked[:5]) if r < 0.85 else rng.choice(arms)
        v, kind, ref = arm[0], arm[1], arm[2]
        commute, kw = VARIANTS[v]
        nodes, dag, rnodes, rdag = nodeset(commute)
        par = R.Params(noise=noise, seed=seed * 100003 + n_tr, w_ext=rng.choice((0.25, 0.5, 0.5, 1.0)),
                       ext=rng.choice((20, 20, 40)), near=rng.choice((0.0, 0.0, 1.0)) if kw["tau"] is None else 0.0, **kw)
        if kind == "greedy":
            l0 = R.greedy_layout(adj, dev, rng)
        elif kind == "line":
            l0 = R.line_layout(list(range(nl)), rng.choice(paths))
            if rng.random() < 0.3:
                l0 = R.line_layout(list(range(nl))[::-1], rng.choice(paths))
        elif kind == "chain":
            l0 = R.line_layout(R.linear_order(adj, rng), rng.choice(paths))
        else:
            l0 = lay_list[given_i % len(lay_list)]
            given_i += 1
        arm[4] += 1
        try:
            l = list(l0)
            for _ in range(ref):
                _, lf, _, _ = R.route(nodes, dag, dev, l, par, want_ops=False)
                _, l, _, _ = R.route(rnodes, rdag, dev, lf, par, want_ops=False)
            ops, lf, cost, nsw = R.route(nodes, dag, dev, l, par)
            routed.append((cost, n_tr, ops, l, lf, par, nsw, commute, f"{v}/{kind}{ref}"))
            if arm[3] is None or cost < arm[3]:
                arm[3] = cost
        except RuntimeError as e:
            if log:
                log(f"router trial failed: {e}")
        n_tr += 1
        if trials is not None and n_tr >= trials:
            break
        if trials is None and time.time() - t_start > seconds and n_tr >= 2:
            break
        if n_tr >= 20000:
            break
    routed.sort(key=lambda r: (r[0], r[1]))
    t_route = time.time() - t_start
    cands = []
    seen = set()
    for cost, _, ops, l, lf, par, nsw, commute, arm_name in routed:
        sig = (round(cost, 6), nsw, tuple(l))
        if sig in seen:
            continue
        seen.add(sig)
        circ = finish(R.to_circuit(nodeset(commute)[0], ops, dev.n), dev)
        init, fin = list(l), [lf[wire[i]] for i in range(nl)]
        if do_relocate:
            circ, init, fin = relocate(circ, init, fin, dev)
        cands.append(Cand(name, circ, init, fin, analyze(circ, dev, fin), 0.0,
                          note=f"{arm_name}; swaps={nsw}; est={cost:.1f}; trials={n_tr}; route_s={t_route:.0f}"))
        if len(cands) >= top:
            break
    total = time.time() - t_start
    for c in cands:
        c.seconds = total
    cands.sort(key=_key(objective))
    return cands


def compile_hw(qc, backend, seeds=(0,), seconds=20.0, trials=None, objective="fid", strategies=None, verify=True,
               max_verify_qubits=26, qopt=True, qopt_seconds=60.0, alternatives=(), log=None):
    """Compile the unitary circuit ``qc`` for ``backend`` (name, backend object or :class:`Device`).

    seeds: seeds of the Qiskit O3 reference runs (each is a candidate).  seconds / trials: budget of the router
    search per logical variant.  strategies: subset of ("qiskit", "qopt+qiskit", "hw", "qopt+hw"); default all
    (the qopt variants only when the circuit has Clifford+T parts that qoptct improves).
    alternatives: equivalent logical circuits (same qubit order) to compile as well, e.g. offline qoptct results;
    they enter as "qopt+qiskit" / "qopt+hw" candidates and are verified against ``qc`` like everything else.
    Returns :class:`HwResult`; ``result.circuit`` is on the physical qubits of the device, the logical qubit i
    starts on ``init_map[i]`` and must be measured on ``final_map[i]``.
    """
    dev = backend if isinstance(backend, Device) else Device(backend)
    if qc.num_clbits and any(i.operation.name == "measure" for i in qc.data):
        qc = qc.remove_final_measurements(inplace=False)
    strategies = set(strategies or ("qiskit", "qopt+qiskit", "hw", "qopt+hw"))
    cands = []
    say = log or (lambda s: None)
    alts = [a.remove_final_measurements(inplace=False) if a.num_clbits else a for a in alternatives]
    if qopt and strategies & {"qopt+qiskit", "qopt+hw"}:
        t0 = time.time()
        try:
            l2, info = qopt_logical(qc, qopt_seconds)
            if info.get("replaced"):
                alts.append(l2)
                say(f"qopt logical: CX {info.get('cx_before')} -> {info.get('cx_after')}, T {info.get('t_before')} -> "
                    f"{info.get('t_after')} ({time.time() - t0:.1f} s)")
            else:
                say(f"qopt logical: {info.get('skipped', 'no change')}")
        except Exception as e:  # never break a compilation because of the optional logical step
            say(f"qopt logical failed: {type(e).__name__}: {e}")
    s0 = seeds[0] if seeds else 0
    lay = []
    if "qiskit" in strategies:
        for sd in seeds:
            c = qiskit_o3(qc, dev, sd)
            cands.append(c); lay.append(c.init_map)
    if "hw" in strategies:
        cands += route_hw(qc, dev, lay[:2], seconds, trials, seed=s0, objective=objective, log=log)
    for k, alt in enumerate(alts):
        lay = []
        if alt.num_qubits != qc.num_qubits:
            say(f"alternative {k}: different number of qubits, skipped")
            continue
        if "qopt+qiskit" in strategies:
            for sd in seeds:
                c = qiskit_o3(alt, dev, sd, "qopt+qiskit")
                c.note += f" alt={k}"
                cands.append(c); lay.append(c.init_map)
        if "qopt+hw" in strategies:
            for c in route_hw(alt, dev, lay[:2], seconds, trials, seed=s0, name="qopt+hw", objective=objective, log=log):
                c.note += f" alt={k}"
                cands.append(c)
    if not cands:
        raise ValueError("no candidate (empty strategy list?)")
    key = _key(objective)
    best = {}
    for c in sorted(cands, key=key):
        best.setdefault(c.name, c)
    if verify:
        for c in best.values():
            c.verified = verify_compiled(qc, c.circuit, c.init_map, c.final_map, max_qubits=max_verify_qubits) \
                if valid_for(c.circuit, dev) else False
    order = sorted(best.values(), key=key)
    win = next((c for c in order if c.verified is True), None) if verify else order[0]
    if win is None:
        # nothing could be verified (too many qubits): fall back to the reference, flagged as unverified
        win = next((c for c in order if c.verified is None), order[0])
    return HwResult(win.circuit, win.init_map, win.final_map, win.metrics, win.name, win.verified, cands, best, qc)
