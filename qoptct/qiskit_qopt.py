"""Qiskit integration of qoptct (project "qopt"): verified CX/T optimization of Clifford+T circuits.

Two ways to use it::

    # 1. as a transpiler pass
    from qiskit.transpiler import PassManager
    from qoptct.qiskit_qopt import QoptPass
    out = PassManager([QoptPass()]).run(qc)

    # 2. as transpiler stage plugins (entry points "qiskit.transpiler.optimization" and
    #    "qiskit.transpiler.init", both under the name "qoptct")
    from qiskit import transpile
    out = transpile(qc, basis_gates=["h", "x", "s", "sdg", "t", "tdg", "cx"], optimization_method="qoptct")
    out = transpile(qc, basis_gates=[...], init_method="qoptct", optimization_method="qoptct")  # Toffoli level too

How the pass works: the circuit is split into maximal runs of consecutive Clifford+T instructions
(h x z s sdg t tdg cx cz ccx ccz, no classical bits). Every run that contains two-qubit (or larger) gates is
optimized by ``qoptct.optimize_hl`` and replaced only if

* the new run is PROVEN equal to the old one (full unitary up to 10 qubits, otherwise PyZX, then a numerical
  check up to 26 qubits; see ``qoptct.prove``),
* it is strictly cheaper by (CX count, T count, size) than the old run written in Clifford+T,
* every new instruction is allowed by the target / basis gates / coupling map that the pass was given.

Anything else (measurements, barriers, control flow, parametrized or unsupported gates) is left untouched and
splits the runs, so the pass can be applied to any circuit; when nothing applies the circuit is returned
unchanged. For runs on at most 10 qubits the global phase is corrected exactly as well.

``QoptMeasPass`` is the measurement-based variant (Gidney logical AND with measurement-based uncomputation)
for qubits that are known to start in |0>.
"""
import cmath

import numpy as np
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.transpiler.preset_passmanagers.plugin import PassManagerStagePlugin

import qoptct

SUP = {"h", "x", "z", "s", "sdg", "t", "tdg", "cx", "cz", "ccx", "ccz"}
MULTI = {"cx", "cz", "ccx", "ccz"}
PHASE_MAX_QUBITS = 10


def to_hl(qc):
    """Toffoli-level gate list [(name, (qubit indices))] of `qc`, or None if it contains anything else."""
    hl = []
    for ins in qc.data:
        g = ins.operation.name
        if g not in SUP or ins.clbits or getattr(ins.operation, "params", None):
            return None
        hl.append((g, tuple(qc.find_bit(q).index for q in ins.qubits)))
    return hl


def from_low(n, c):
    qc = QuantumCircuit(n)
    for g, q in c:
        getattr(qc, g)(*q)
    return qc


def _supported(ins):
    op = ins.operation
    return (op.name in SUP and not ins.clbits and not getattr(op, "params", None)
            and getattr(op, "condition", None) is None)


def _cost(c):
    cx = sum(g == "cx" for g, _ in c)
    t = sum(g in ("t", "tdg") for g, _ in c)
    return cx, t, len(c)


class QoptPass(TransformationPass):
    """Verified CX/T optimization of the Clifford+T parts of a circuit.

    Args:
        net: also use neural block resynthesis (needs PyTorch: ``pip install "qoptct[net]"``).
        target: optional :class:`~qiskit.transpiler.Target`; new instructions must be supported by it.
        basis_gates: optional list of allowed gate names (used when no target is given).
        coupling_map: optional :class:`~qiskit.transpiler.CouplingMap` (used when no target is given).
        min_multi: runs with fewer multi-qubit gates than this are not touched (default 1).

    After ``run`` the attribute ``last`` holds a summary: ``runs`` (candidate runs), ``replaced``,
    ``verified`` (True if at least one run was replaced; every replacement is proven), ``method``,
    ``cx_before``/``cx_after`` and ``t_before``/``t_after`` of the replaced runs, or ``skipped`` with the reason.
    """

    def __init__(self, net=False, target=None, basis_gates=None, coupling_map=None, min_multi=1):
        super().__init__()
        self.net = net
        self.target = target
        self.basis = set(basis_gates) if basis_gates else None
        self.cmap = coupling_map
        self.min_multi = min_multi
        self.last = None
        if net:
            qoptct.require_torch()

    # -- target checks ---------------------------------------------------------------------------------------
    def _allowed(self, name, qargs):
        if self.target is not None:
            if name not in self.target.operation_names:
                return False
            return self.target.instruction_supported(name, tuple(qargs))
        if self.basis is not None and name not in self.basis:
            return False
        if self.cmap is not None and len(qargs) == 2:
            return tuple(qargs) in set(self.cmap.get_edges())
        return True

    # -- main ------------------------------------------------------------------------------------------------
    def run(self, dag):
        qc = dag_to_circuit(dag)
        runs, cur = [], []
        for ins in qc.data:
            if _supported(ins):
                cur.append(ins)
            else:
                if cur:
                    runs.append(cur)
                runs.append(ins)
                cur = []
        if cur:
            runs.append(cur)
        cand = [r for r in runs if isinstance(r, list) and sum(i.operation.name in MULTI for i in r) >= self.min_multi]
        if not cand:
            self.last = {"skipped": "no Clifford+T part with multi-qubit gates", "runs": 0, "replaced": 0,
                         "verified": False}
            return dag
        cand_ids = {id(r) for r in cand}
        out = qc.copy_empty_like()
        stats = {"runs": len(cand), "replaced": 0, "method": [], "cx_before": 0, "cx_after": 0,
                 "t_before": 0, "t_after": 0}
        for r in runs:
            if not isinstance(r, list):
                out.append(r.operation, r.qubits, r.clbits)
                continue
            new = self._optimize_run(qc, r, stats) if id(r) in cand_ids else None
            if new is None:
                for ins in r:
                    out.append(ins.operation, ins.qubits, ins.clbits)
            else:
                ops, phase = new
                for g, qs in ops:
                    getattr(out, g)(*qs)
                out.global_phase += phase
        stats["verified"] = stats["replaced"] > 0
        if not stats["replaced"]:
            stats["skipped"] = "no proven improvement"
        stats["method"] = ", ".join(m for m in stats["method"] if m) or None
        self.last = stats
        if not stats["replaced"]:
            return dag
        return circuit_to_dag(out)

    def _optimize_run(self, qc, run, stats):
        """Optimized replacement [(name, (Qubit, ...))] and global-phase correction for one run, or None."""
        qubits = sorted({q for ins in run for q in ins.qubits}, key=lambda q: qc.find_bit(q).index)
        loc = {q: i for i, q in enumerate(qubits)}
        hl = [(ins.operation.name, tuple(loc[q] for q in ins.qubits)) for ins in run]
        n = len(qubits)
        names = [str(i) for i in range(n)]
        try:
            _, ref, best, how = qoptct.optimize_hl(names, hl, self.net)
        except ImportError:
            raise
        except Exception:  # never break a transpilation: keep the run as it is
            return None
        if how is None:
            return None
        old_low = [(g, tuple(q)) for g, q in ref]
        if _cost(best) >= _cost(old_low):
            return None
        ops = [(g, tuple(qubits[i] for i in q)) for g, q in best]
        phys = [(g, tuple(qc.find_bit(q).index for q in qs)) for g, qs in ops]
        if not all(self._allowed(g, qs) for g, qs in phys):
            return None
        phase = 0.0
        if n <= PHASE_MAX_QUBITS:
            from qiskit.quantum_info import Operator
            sub = QuantumCircuit(n)
            for g, q in hl:
                getattr(sub, g)(*q)
            u_old = Operator(sub).data
            u_new = Operator(from_low(n, best)).data
            ov = np.vdot(u_new, u_old)  # = tr(U_new^dagger U_old) = e^{i phi} * dim if U_old = e^{i phi} U_new
            if abs(abs(ov) - u_old.shape[0]) > 1e-6:  # must be equal up to a global phase (already proven)
                return None
            phase = cmath.phase(ov)
        cb, tb, _ = _cost(old_low); ca, ta, _ = _cost(best)
        stats["replaced"] += 1; stats["method"].append(how)
        stats["cx_before"] += cb; stats["cx_after"] += ca; stats["t_before"] += tb; stats["t_after"] += ta
        return ops, phase


def _default_optimization(cfg, level):
    """Qiskit's own optimization stage for this configuration: the dedicated Clifford+T stage when the basis is
    Clifford+T and the installed Qiskit has it (Qiskit >= 2.1), otherwise the general one."""
    from qiskit.transpiler.preset_passmanagers import builtin_plugins as bp
    try:
        from qiskit.transpiler.passmanager_config import PassManagerCliffordTConfig
        from qiskit.transpiler.preset_passmanagers.common import is_clifford_t_basis
        if is_clifford_t_basis(cfg.basis_gates, cfg.target):
            return bp.OptimizeCliffordTPassManager().pass_manager(
                PassManagerCliffordTConfig._from_legacy_config(cfg), optimization_level=level)
    except (ImportError, AttributeError):
        pass
    return bp.OptimizationPassManager().pass_manager(cfg, level)


class QoptOptimizationPlugin(PassManagerStagePlugin):
    """Optimization-stage plugin ``optimization_method="qoptct"``.

    Runs Qiskit's own optimization stage of the requested level (for a Clifford+T basis the dedicated Clifford+T
    stage of Qiskit >= 2.1; nothing for level 0; level 2 if no level is given) and then
    :class:`QoptPass` restricted to the target / basis gates / coupling map of the transpilation, so the output
    stays valid for the backend. Parts of the circuit that are not in Clifford+T are left to Qiskit.
    """

    def pass_manager(self, pass_manager_config, optimization_level=None):
        from qiskit.transpiler import PassManager
        level = 2 if optimization_level is None else optimization_level
        pm = _default_optimization(pass_manager_config, level) or PassManager()
        pm.append(QoptPass(target=pass_manager_config.target,
                           basis_gates=pass_manager_config.basis_gates,
                           coupling_map=pass_manager_config.coupling_map))
        return pm


def _default_init(cfg, level):
    """Qiskit's own init stage for this configuration (dedicated Clifford+T variant when available)."""
    from qiskit.transpiler.preset_passmanagers import builtin_plugins as bp
    try:
        from qiskit.transpiler.passmanager_config import PassManagerCliffordTConfig
        from qiskit.transpiler.preset_passmanagers.common import is_clifford_t_basis
        if is_clifford_t_basis(cfg.basis_gates, cfg.target):
            return bp.CliffordTInitPassManager().pass_manager(
                PassManagerCliffordTConfig._from_legacy_config(cfg), optimization_level=level)
    except (ImportError, AttributeError):
        pass
    return bp.DefaultInitPassManager().pass_manager(cfg, level)


class QoptInitPlugin(PassManagerStagePlugin):
    """Init-stage plugin ``init_method="qoptct"`` (optional companion of the optimization plugin).

    Runs :class:`QoptPass` on the circuit as the user wrote it, i.e. still at the Toffoli level (ccx, ccz, cx, ...),
    BEFORE Qiskit decomposes the Toffolis; this is where qopt's rewriting programs find most savings. No
    layout exists yet, so the pass is unconstrained here; Qiskit's own init stage, layout, routing and
    translation follow as usual. Use it together with the optimization plugin::

        transpile(qc, basis_gates=..., init_method="qoptct", optimization_method="qoptct")
    """

    def pass_manager(self, pass_manager_config, optimization_level=None):
        from qiskit.transpiler import PassManager
        level = 2 if optimization_level is None else optimization_level
        pm = PassManager([QoptPass()])
        default = _default_init(pass_manager_config, level)
        if default is not None:
            pm.append(default.to_flow_controller())
        return pm


class QoptMeasPass(TransformationPass):
    """Measurement-based mode: Gidney logical AND with measurement-based uncomputation.

    ``ancillas`` are indices of qubits that are guaranteed to be |0> at the input. A Toffoli into a clean qubit
    costs 4 T; erasing a known product a*b becomes measure + if_test(cz, x). Proof: GF(2) tracking or a
    simulation with measurements (gidney.py). The whole circuit must consist of the gates listed in ``SUP``.
    If the result is not proven or the T count does not drop, the circuit is returned unchanged.
    """

    def __init__(self, ancillas):
        super().__init__(); self.anc = set(ancillas); self.last = None

    def run(self, dag):
        import gidney, evo_lib
        from qiskit import ClassicalRegister
        qc = dag_to_circuit(dag); hl = to_hl(qc)
        if not hl:
            self.last = {"skipped": "unsupported instructions"}; return dag
        n = qc.num_qubits; hl = evo_lib.fold_h(hl)
        items, r1, r2 = gidney.transform(hl, n, self.anc)
        if not (r1 or r2):
            self.last = {"skipped": "no logical-AND sites"}; return dag
        import meas_portfolio
        names = [str(i) for i in range(n)]
        low = meas_portfolio.build_low(names, items, n)
        ok = gidney.anf_check(items, n, self.anc)
        if ok is None:
            ok = gidney.verify(to_hl(qc), low, n, self.anc)
        T = sum(evo_lib.tcount(s) for s, _ in low); T0 = sum(g in ("t", "tdg") for g, _ in evo_lib.lower_std(to_hl(qc)))
        self.last = {"t_before": T0, "t_after": T, "verified": ok, "measurements": r2}
        if ok is not True or T >= T0:
            return dag
        out = QuantumCircuit(n); creg = ClassicalRegister(max(r2, 1), "m"); out.add_register(creg); k = 0
        for seg, mz in low:
            for g, q in seg: getattr(out, g)(*q)
            if mz:
                a, b, t = mz; out.measure(t, creg[k])
                with out.if_test((creg[k], 1)):
                    out.cz(a, b); out.x(t)
                k += 1
        return circuit_to_dag(out)
