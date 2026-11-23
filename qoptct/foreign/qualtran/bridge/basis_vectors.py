"""Exact four Qualtran computational-basis state/effect atoms.

IntState/IntEffect are NOT recognized here: their own source BloqGraph must
decompose to these four atoms. Effects are partial maps with an OPEN input
promise. Retiring a wire is never proof of its promised value. The full
physical final-zero contract, together with checked no-later-use, discharges
that promise. domain_guard alone covers AND/ERASE, not these effects.
"""

from copy import deepcopy
from numbers import Integral

from .qualtran_layout import AtomExpansion, LayoutError, controlled_and_atom_rule
from .sim_parse import parse_qasm
from .strict_proof import digest_json


class BasisVectorError(LayoutError):
    pass


def _wire(array, path):
    import numpy as np
    array = np.asarray(array, dtype=object)
    if array.shape != (1,):
        raise BasisVectorError(f"{path}: basis atom needs one physical QBit")
    q = array[0]
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 0:
        raise BasisVectorError(f"{path}: invalid physical wire {q!r}")
    return int(q)


def _op_wires(op):
    if op[0] in ("x", "cx", "and"):
        return set(op[1:])
    if op[0] == "erase":
        return {op[1]} | {q for correction in op[2] for q in correction[1:]}
    raise BasisVectorError(f"Unknown compiled atom while checking retired-wire use: {op}")


class BasisVectorRules:
    """One fresh callback instance per lower_bloqgraph call.

    Use extra_atom_rule=book.chain_with_a1. After ordinary translate(), call
    book.bind(lowered, compiled), persist its obligations in the contract, and
    require the normal full physical output proof. bind returns structural
    evidence only; semantic_status remains UNPROVED.
    """

    def __init__(self):
        self.obligations = []
        self._paths = set()
        self._bound = False

    def __call__(self, bloq, inputs, allocate, path):
        from qualtran import QBit, Side
        from qualtran.bloqs.basic_gates import ZeroState, OneState, ZeroEffect, OneEffect
        table = {ZeroState: (True, 0), OneState: (True, 1),
                 ZeroEffect: (False, 0), OneEffect: (False, 1)}
        if type(bloq) not in table:
            return None
        if self._bound or path in self._paths:
            raise BasisVectorError("A BasisVectorRules instance cannot be reused across source graphs")
        is_state, bit = table[type(bloq)]
        if type(bloq.state) is not bool or bloq.state != is_state or type(bloq.bit) is not bool or int(bloq.bit) != bit or type(bloq.n) is not int or bloq.n != 1:
            raise BasisVectorError(f"{path}: source basis atom fields differ from exact class semantics")
        regs = tuple(bloq.signature)
        if len(regs) != 1 or regs[0].name != "q" or type(regs[0].dtype) is not QBit or regs[0].shape or regs[0].side != (Side.RIGHT if is_state else Side.LEFT):
            raise BasisVectorError(f"{path}: unexpected source basis-vector signature")
        if set(inputs) != (set() if is_state else {"q"}):
            raise BasisVectorError(f"{path}: source state/effect input signature mismatch")
        if is_state:
            outputs = allocate(regs, path)
            if set(outputs) != {"q"}:
                raise BasisVectorError(f"{path}: allocator changed the source signature")
            q = _wire(outputs["q"], path)
            retired = ()
        else:
            q = _wire(inputs["q"], path)
            outputs, retired = {}, (q,)
        leaf = dict(name="x", qubits=[q], path=path,
                    basis_vector_rule=type(bloq).__name__, basis_value=bit)
        self._paths.add(path)
        self.obligations.append(dict(path=path, source_type=type(bloq).__name__, wire=q,
                                     kind="state" if is_state else "effect", source_bit=bit,
                                     initial_zero_required=is_state, required_input_bit=None if is_state else bit,
                                     emitted_x=bool(bit), final_zero_required=not is_state,
                                     no_later_use_required=not is_state, semantic_status="UNPROVED"))
        return AtomExpansion(outputs, [leaf] if bit else [], retired)

    def chain_with_a1(self, bloq, inputs, allocate, path):
        expansion = self(bloq, inputs, allocate, path)
        return expansion if expansion is not None else controlled_and_atom_rule(bloq, inputs, allocate, path)

    def bind(self, lowered, compiled):
        """Bind the OPEN promises to actual graph, trace and parsed QASM.

        A wrong source effect is not rejected here merely by fiat. It leaves
        the corresponding physical wire nonzero, which the complete output
        specification must refute. No wire can be recycled or touched later.
        """
        if self._bound:
            raise BasisVectorError("bind may be called once per source graph")
        if lowered.get("source_backend") != "qualtran_bloqgraph" or lowered.get("allocation_policy") != "fresh physical wire for each allocation; never recycle":
            raise BasisVectorError("Basis effects require the explicit fresh-wire BloqGraph backend")
        nq, parsed = parse_qasm(compiled["qasm"])
        if nq != compiled["nq"] or parsed != compiled["ops"]:
            raise BasisVectorError("Actual QASM does not match the compiled operations")
        layout = compiled["layout"]
        if layout["inputs"] != lowered["inputs"] or layout["outputs"] != lowered["outputs"]:
            raise BasisVectorError("Compiled layout changed the source signature")
        if layout.get("domain_policy") != "exact_logical_AND_ERASE_v1" or layout.get("garbage") or layout.get("free_wires"):
            raise BasisVectorError("Effects require a full exact logical contract without unspecified wires")
        output_wires = {q for wires in layout["outputs"].values() for q in wires}
        initial_zero = set(layout["initial_zero"])
        retired = set(lowered["retired_wires"])
        events = lowered["graph_events"]
        indices = {event["path"]: i for i, event in enumerate(events)}
        if len(indices) != len(events):
            raise BasisVectorError("Ambiguous repeated source event path")
        basis_types = {"ZeroState", "OneState", "ZeroEffect", "OneEffect"}
        if {event["path"] for event in events if event["bloq_type"] in basis_types} != self._paths:
            raise BasisVectorError("Not all source basis-vector atoms have recorded obligations")
        leaves, trace = lowered["leaves"], compiled["trace"]
        if len(leaves) != len(trace):
            raise BasisVectorError("Source/compiled trace size differs")
        previous_end, previous_event = 0, -1
        for i, (leaf, row) in enumerate(zip(leaves, trace)):
            event_index = indices.get(leaf.get("path"), -1)
            if (row.get("source_index") != i or row.get("leaf") != leaf
                    or row["logical_start"] != previous_end or event_index < 0 or event_index < previous_event
                    or not previous_end <= row["logical_end"] <= len(parsed)):
                raise BasisVectorError("Source traversal order or exact trace was changed")
            previous_end, previous_event = row["logical_end"], event_index
        if previous_end != len(parsed):
            raise BasisVectorError("Compiled trace does not cover the entire QASM")
        seen_states = set()
        bound = []
        for obligation in self.obligations:
            item = deepcopy(obligation)
            q, path = item["wire"], item["path"]
            index = indices[path]
            event = events[index]
            expected_inputs = {} if item["kind"] == "state" else {"q": [q]}
            expected_outputs = {"q": [q]} if item["kind"] == "state" else {}
            if event["inputs_msb"] != expected_inputs or event["outputs_msb"] != expected_outputs:
                raise BasisVectorError(f"{path}: recorded source wiring differs")
            own = [row for row in trace if row["leaf"]["path"] == path]
            if item["emitted_x"]:
                if len(own) != 1 or parsed[own[0]["logical_start"]:own[0]["logical_end"]] != [("x", q)]:
                    raise BasisVectorError(f"{path}: expected exactly one X")
            elif own:
                raise BasisVectorError(f"{path}: zero-vector atom unexpectedly emitted operations")
            if item["kind"] == "state":
                if q not in initial_zero or q in seen_states:
                    raise BasisVectorError(f"{path}: state allocation is not fresh and initially zero")
                seen_states.add(q)
                if any(q in leaf["qubits"] for leaf in leaves if indices[leaf["path"]] < index):
                    raise BasisVectorError(f"{path}: allocated state wire was used before creation")
            else:
                if q not in retired or q in output_wires or q not in lowered["final_zero"]:
                    raise BasisVectorError(f"{path}: effect wire missing from full final-zero obligation")
                cut = max((row["logical_end"] for row in trace if indices[row["leaf"]["path"]] <= index), default=0)
                if any(q in _op_wires(op) for op in parsed[cut:]):
                    raise BasisVectorError(f"{path}: retired effect wire is used later")
                item.update(logical_boundary_after_effect=cut,
                            discharge="full physical nonoutput q[%d]=0 after this boundary; no later use verified" % q)
            item["structural_status"] = "PASS"
            bound.append(item)
        self._bound = True
        return dict(status="PASS", scope="source basis-vector structure only", semantic_status="UNPROVED",
                    obligations=bound, obligations_sha256=digest_json(bound),
                    layout_sha256=digest_json(layout), source_backend="qualtran_bloqgraph",
                    requires_all_physical_nonoutputs_zero=True,
                    domain_guard_alone_sufficient=False)
