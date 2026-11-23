""""Black box" rule for lower_bloqgraph: the Montgomery multiplier and the Kaliski inverse are not expanded.

Leaf {'name': 'opaque', 'kind': MM|MM_adj|KINV|KINV_adj, 'regs_in', 'regs_out', 'bloq'}.
The wires of the registers are lists of physical numbers in MSB-first order (as in lower_bloqgraph).
"""
from qoptct.foreign.qualtran.bridge.qualtran_layout import AtomExpansion


def opaque_kind(bloq):
    from qualtran.bloqs.mod_arithmetic import DirtyOutOfPlaceMontgomeryModMul, KaliskiModInverse
    t = type(bloq)
    if t is DirtyOutOfPlaceMontgomeryModMul:
        return 'MM_adj' if bloq.uncompute else 'MM'
    if t is KaliskiModInverse:
        return 'KINV_adj' if bloq.uncompute else 'KINV'
    return None


def state_atom(bloq, inputs, allocate, path):
    """ZeroState / IntState: fresh zero wires (+ x on the one bits of val). MSB-first."""
    from qualtran.bloqs.basic_gates import ZeroState, IntState
    t = type(bloq)
    if t not in (ZeroState, IntState):
        return None
    rights = list(bloq.signature.rights())
    if len(rights) != 1 or list(bloq.signature.lefts()):
        raise ValueError(f'{path}: unexpected state signature')
    out = allocate(rights, path)
    (name, arr), = out.items()
    wires = [int(q) for q in arr.flat]
    val = 0 if t is ZeroState else int(bloq.val)
    if not 0 <= val < (1 << len(wires)):
        raise ValueError(f'{path}: IntState value out of range')
    leaves = [{'name': 'x', 'qubits': [w], 'path': path}
              for i, w in enumerate(wires) if (val >> (len(wires) - 1 - i)) & 1]
    return AtomExpansion(out, leaves)


def opaque_atom_rule(bloq, inputs, allocate, path):
    from qualtran import Side
    st = state_atom(bloq, inputs, allocate, path)
    if st is not None:
        return st
    kind = opaque_kind(bloq)
    if kind is None:
        return None
    rights_only = [r for r in bloq.signature if r.side == Side.RIGHT]
    allocated = allocate(rights_only, path)
    regs_in = {name: [int(q) for q in arr.flat] for name, arr in inputs.items()}
    outputs = {}
    for r in bloq.signature.rights():
        outputs[r.name] = allocated[r.name] if r.side == Side.RIGHT else inputs[r.name]
    regs_out = {name: [int(q) for q in arr.flat] for name, arr in outputs.items()}
    retired = tuple(q for r in bloq.signature if r.side == Side.LEFT for q in regs_in[r.name])
    qubits = sorted({q for w in list(regs_in.values()) + list(regs_out.values()) for q in w})
    leaf = dict(name='opaque', kind=kind, regs_in=regs_in, regs_out=regs_out, qubits=qubits,
                bloq_repr=repr(bloq), path=path)
    return AtomExpansion(outputs, [leaf], retired)
