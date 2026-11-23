"""Two explicit sources: the original Cirq export and the exact BloqGraph.

Both use only the own decompositions of Qualtran. The old lower is not
replaced by the graph path: it needs the graph for the mandatory THRU guard.
"""
from .rules import TranslationError


def lower_graph(bloq, max_leaves=200000, *, max_nodes=1000000, max_depth=128,
                extra_atom_rule=None):
    """Return a self-contained qualtran_bloqgraph with exact physical soquets.

    Shaped arrays and the full LEFT/RIGHT/THRU interface are kept. A1 is plugged
    in explicitly; the semantics of its atoms and the call domain are checked in the later
    translation/proof. No result is attributed to the Cirq export.
    """
    from .qualtran_layout import lower_bloqgraph, controlled_and_atom_rule
    return lower_bloqgraph(bloq, extra_atom_rule=extra_atom_rule or controlled_and_atom_rule,
                          max_leaves=max_leaves, max_nodes=max_nodes, max_depth=max_depth)


def lower(bloq, max_leaves=200000):
    """The old Cirq path; it must check THRU BEFORE the export itself.

    An unknown BloqGraph, a non-topological order or a non-identity THRU
    permutation give a strict refusal, even if Cirq can export something.
    There is no fallback that changes the source behind the back of the client.
    """
    from .qualtran_layout import require_cirq_thru_identity
    layout_guard = require_cirq_thru_identity(lower_graph(bloq, max_leaves=max_leaves))
    import cirq
    import numpy as np
    from qualtran._infra.gate_with_registers import get_named_qubits
    from qualtran.bloqs.mcmt import And
    from qualtran.bloqs.basic_gates import XGate, CNOT, Toffoli
    from qualtran.cirq_interop import BloqAsCirqGate
    for reg in bloq.signature:
        if reg.shape:
            raise TranslationError(f'Register {reg.name}: shape {reg.shape} is not supported by the layout yet')
    inregs=get_named_qubits(bloq.signature.lefts())
    class CleanManager(cirq.SimpleQubitManager):
        def qborrow(self,n,dim=2):
            raise TranslationError('A borrowed dirty ancilla (qborrow) is not allowed by the clean model')
    manager=CleanManager(prefix='bridge_source')
    context=cirq.DecompositionContext(manager)
    try:
        circuit,outregs=bloq.as_composite_bloq().to_cirq_circuit_and_quregs(qubit_manager=manager,**inregs)
    except Exception as e:
        raise TranslationError(f'Export of {bloq}: {type(e).__name__}: {e}') from e
    leaves=[]; native=[]
    def walk(op,path,depth):
        if depth>128:raise TranslationError(f'{path}: decomposition depth >128')
        if isinstance(op,cirq.ClassicallyControlledOperation):
            raise TranslationError(f'{path}: separate classical control {op}')
        gate=op.gate
        name=None; extra={}
        from .controlled_and import match as match_controlled_and
        controlled_leaf=match_controlled_and(op,path)
        if controlled_leaf is not None:
            name=controlled_leaf['name']
            extra={k:v for k,v in controlled_leaf.items() if k not in ('name','qubits','path')}
        elif type(gate) is BloqAsCirqGate and type(gate.bloq) in (XGate,CNOT,Toffoli):
            name={XGate:'x',CNOT:'cx',Toffoli:'ccx'}[type(gate.bloq)]
        elif type(gate) is And:
            name='erase' if gate.uncompute else 'and'
            extra['controls']=[int(gate.cv1),int(gate.cv2)]
        else:
            for cls, atom, nm in ((cirq.XPowGate,cirq.X,'x'),(cirq.CXPowGate,cirq.CNOT,'cx'),(cirq.CCXPowGate,cirq.CCX,'ccx'),(cirq.SwapPowGate,cirq.SWAP,'swap')):
                if type(gate) in (cls,type(atom)) and gate.exponent in (-1,1) and gate.global_shift==0:
                    name=nm;break
        if name is not None:
            leaves.append(dict(name=name,qubits=list(op.qubits),path=path,**extra));native.append(op)
            if len(leaves)>max_leaves:raise TranslationError(f'{path}: leaves >{max_leaves}')
            return
        # No cirq.decompose fallback, unitary synthesis or optimizers.
        if gate is not None and isinstance(gate,(cirq.EigenGate,cirq.MeasurementGate,cirq.ResetChannel)):
            raise TranslationError(f'{path}: unsupported source gate {gate!r}')
        try:
            children=cirq.decompose_once(op,default=None,context=context)
        except Exception as e:
            raise TranslationError(f'{path}: decomposition of {op}: {type(e).__name__}: {e}') from e
        if children is None:
            raise TranslationError(f'{path}: no own decomposition {op}')
        for j,child in enumerate(cirq.flatten_to_ops(children)):
            walk(child,f'{path}/{j}:{child.gate}',depth+1)
    for j,op in enumerate(circuit.all_operations()):walk(op,f'{j}:{op.gate}',0)
    wires=[]
    for regs in (inregs,outregs):
        for arr in regs.values():
            for q in np.asarray(arr).flat:
                if q not in wires:wires.append(q)
    for op in native:
        for q in op.qubits:
            if q not in wires:wires.append(q)
    indexes={q:i for i,q in enumerate(wires)}
    for leaf in leaves:leaf['qubits']=[indexes[q] for q in leaf['qubits']]
    def mapped(regs):return {k:[indexes[q] for q in np.asarray(v).ravel()[::-1]] for k,v in regs.items()}
    return dict(nq=len(wires),leaves=leaves,inputs=mapped(inregs),outputs=mapped(outregs),
                wire_names=list(map(str,wires)),source=circuit,native_ops=native,wire_objects=wires,bloq=bloq,
                source_backend='qualtran_cirq_export', layout_guard=layout_guard)
