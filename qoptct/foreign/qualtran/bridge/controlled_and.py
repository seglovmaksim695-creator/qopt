"""A1: explicit rule for singly controlled Qualtran And / And-dagger.

This translates a partial isometry; it does not prove the calling context.
And: target=0. And-dagger: target=[c=cv]*[a=cv1]*[b=cv2]. Both need a clean scratch wire.
Measurement operations get no quantum control; their phase correction
accounts for the external control through an additional computed AND.
"""

from numbers import Integral


class ControlledAndError(ValueError):
    pass


def _cv(value, path):
    if not isinstance(value, Integral) or isinstance(value, bool) or value not in (0, 1):
        raise ControlledAndError(f"{path}: the control value must be an explicit 0 or 1: {value!r}")
    return int(value)


def match_gate(gate, qubits, path="root"):
    """Return the leaf for the exact source class, None for any other type.

    qubits: full signature order: external control, a, b, target.
    Works for the physical wires of a Bloq graph and for Cirq Qids before indexing.
    A known controlled-And of an unsupported form is rejected explicitly.
    """
    import cirq
    from qualtran import CtrlSpec, QBit
    from qualtran._infra.controlled import Controlled
    from qualtran.bloqs.mcmt import And
    from qualtran.bloqs.mcmt.controlled_via_and import ControlledViaAnd

    if type(gate) in (Controlled, ControlledViaAnd):
        if type(gate.subbloq) is not And:
            return None
        sub = gate.subbloq
        spec = gate.ctrl_spec
        if type(spec) is not CtrlSpec or spec.is_symbolic():
            raise ControlledAndError(f"{path}: symbolic or unknown CtrlSpec")
        if spec.num_bits != 1 or len(spec.qdtypes) != 1 or type(spec.qdtypes[0]) is not QBit:
            raise ControlledAndError(f"{path}: A1 supports exactly one external QBit; {spec!r}")
        cv = _cv(spec.get_single_ctrl_val(), path)
    elif type(gate) is cirq.ControlledGate:
        if type(gate.sub_gate) is not And:
            return None
        sub = gate.sub_gate
        if gate.num_controls() != 1 or tuple(gate.control_qid_shape) != (2,):
            raise ControlledAndError(f"{path}: A1 needs one binary external control")
        if type(gate.control_values) is not cirq.ProductOfSums:
            raise ControlledAndError(f"{path}: only ProductOfSums is allowed, not an arbitrary control formula")
        terms = tuple(gate.control_values.expand())
        if len(terms) != 1 or len(terms[0]) != 1:
            raise ControlledAndError(f"{path}: exactly one activation value 0/1 is required")
        cv = _cv(terms[0][0], path)
    else:
        return None

    if type(sub.uncompute) is not bool:
        raise ControlledAndError(f"{path}: the uncompute parameter must be a bool")
    cvs = [cv, _cv(sub.cv1, path), _cv(sub.cv2, path)]
    qubits = list(qubits)
    if len(qubits) != 4 or len(set(qubits)) != 4:
        raise ControlledAndError(f"{path}: four distinct wires of controlled-And are needed")
    if tuple(cirq.qid_shape(gate)) != (2, 2, 2, 2):
        raise ControlledAndError(f"{path}: non-binary form of the source gate")
    return dict(name="controlled_erase" if sub.uncompute else "controlled_and",
                qubits=qubits, controls=cvs, path=path,
                source_class=f"{type(gate).__module__}.{type(gate).__name__}",
                context_domain_required=True)


def match(op, path="root"):
    """Cirq wrapper; a separate classical control is not accepted."""
    import cirq
    if isinstance(op, (cirq.ClassicallyControlledOperation, cirq.TaggedOperation)):
        return None
    qubits = getattr(op, "qubits", ())
    if any(getattr(q, "dimension", None) != 2 for q in qubits):
        raise ControlledAndError(f"{path}: non-binary Qids are not supported")
    return match_gate(getattr(op, "gate", None), qubits, path)


def _params(leaf):
    if not isinstance(leaf, dict) or leaf.get("name") not in ("controlled_and", "controlled_erase"):
        raise ControlledAndError("A1: unknown leaf")
    path = leaf.get("path", "root")
    wires = leaf.get("qubits", ())
    values = leaf.get("controls", ())
    if len(wires) != 4 or len(set(wires)) != 4 or len(values) != 3:
        raise ControlledAndError(f"{path}: wrong form of an A1 leaf")
    return tuple(wires), tuple(_cv(v, path) for v in values), path


def require_domain(leaf, values, algebra, *, ancilla=None):
    """Exact domain check in the canonical Boolean algebra of the calling code.

    values[wire] is a BDD edge or a bit mask over all inputs; algebra has
    F, T, and_, not_. Equality of expressions must be exact and canonical.
    With plain 0/1 values only the concrete input is checked, not all inputs of the circuit.
    The values object must describe the current state, not an old boundary.
    """
    (c, a, b, t), cvs, path = _params(leaf)
    try:
        product = algebra.T
        for w, cv in zip((c, a, b), cvs):
            literal = values[w] if cv else algebra.not_(values[w])
            product = algebra.and_(product, literal)
        expected = algebra.F if leaf["name"] == "controlled_and" else product
        if values[t] != expected:
            raise ControlledAndError(f"{path}: domain not proved: target={'0' if leaf['name'] == 'controlled_and' else 'controlled predicate'}")
        if ancilla is not None and values[ancilla] != algebra.F:
            raise ControlledAndError(f"{path}: scratch is not proved to be zero")
    except KeyError as exc:
        raise ControlledAndError(f"{path}: no expression for wire {exc}") from exc
    return True


def expand(leaf, *, ancilla, nq, known_zero):
    """Return ops and the open obligations; gives no arithmetic verdict.

    nq is the number of existing wires. ancilla==nq creates one new clean
    wire; ancilla<nq is allowed only if it is in known_zero.
    For compute the target must be in known_zero before the translation.
    For erase the full t=predicate must be proved by the exact checker of the circuit.
    """
    (c, a, b, t), cvs, path = _params(leaf)
    if type(nq) is not int or nq < 4 or any(type(w) is not int or not 0 <= w < nq for w in (c, a, b, t)):
        raise ControlledAndError(f"{path}: invalid source wire indices")
    if type(ancilla) is not int or not 0 <= ancilla <= nq or ancilla in (c, a, b, t):
        raise ControlledAndError(f"{path}: scratch must be a separate existing wire or one new wire")
    zero = set(known_zero)
    total_nq = max(nq, ancilla + 1)
    if any(type(w) is not int or not 0 <= w < total_nq for w in zero):
        raise ControlledAndError(f"{path}: wrong set of known zeros")
    if ancilla < nq and ancilla not in zero:
        raise ControlledAndError(f"{path}: existing scratch is not proved clean")
    compute = leaf["name"] == "controlled_and"
    if compute and t not in zero:
        raise ControlledAndError(f"{path}: compute target is not proved clean")
    flips = [("x", w) for w, cv in zip((c, a, b), cvs) if cv == 0]
    operations = [("and", c, a, ancilla)]
    if compute:
        operations.append(("and", ancilla, b, t))
    else:
        operations.append(("erase", t, (("cz", ancilla, b),)))
    operations.append(("erase", ancilla, (("cz", c, a),)))
    obligations = [dict(kind="zero", wire=ancilla,
                        evidence="fresh wire" if ancilla == nq else "known_zero")]
    if compute:
        obligations.append(dict(kind="zero", wire=t, evidence="known_zero"))
    else:
        obligations.append(dict(kind="equals_controlled_predicate", wire=t,
                                literals=[dict(wire=w, value=cv) for w, cv in zip((c, a, b), cvs)],
                                evidence=None,
                                discharge="exact domain_guard or segment proof: t XOR (s*b) = 0 under pre; kprove alone insufficient"))
    return dict(nq=total_nq, ops=flips + operations + list(reversed(flips)),
                scratch=ancilla, newly_allocated=[ancilla] if ancilla == nq else [],
                clean_on_exit=[ancilla] + ([] if compute else [t]),
                domain_obligations=obligations, context_proved=False,
                measurement_count=1 if compute else 2,
                source_path=path, source_rule="A1.single_controlled_and.v1")
