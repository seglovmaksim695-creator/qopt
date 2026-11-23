"""Opt-in exact TwoBitCSwap atom; no core translator registration.

The source has native quantum tensors and a native Clifford+T circuit.
This callback emits CX(y,x), CCX(ctrl,x,y), CX(y,x). The unchanged
translator implements CCX using its clean AND/CX/ERASE scratch rule.
Structural binding is NOT a semantic or whole-program proof.
"""
from copy import deepcopy
from hashlib import sha256
from numbers import Integral

from .basis_vectors import BasisVectorRules
from .qualtran_layout import AtomExpansion, LayoutError
from .sim_parse import parse_qasm
from .strict_proof import check_layout, digest_json

RULE_ID = 'qualtran.TwoBitCSwap.exact-fredkin.v1'


class FredkinError(LayoutError):
    pass


def match(bloq, inputs, path):
    """Return physical ctrl/x/y for the exact native type, otherwise None."""
    import numpy as np
    from qualtran import QBit, Side
    from qualtran.bloqs.basic_gates import TwoBitCSwap

    if type(bloq) is not TwoBitCSwap:
        return None
    regs = tuple(bloq.signature)
    if (tuple(r.name for r in regs) != ('ctrl', 'x', 'y')
            or any(type(r.dtype) is not QBit or r.shape or r.side != Side.THRU for r in regs)
            or set(inputs) != {'ctrl', 'x', 'y'}):
        raise FredkinError(f'{path}: TwoBitCSwap signature differs from three THRU QBits')
    wires = []
    for name in ('ctrl', 'x', 'y'):
        arr = np.asarray(inputs[name], dtype=object)
        if arr.shape != (1,):
            raise FredkinError(f'{path}/{name}: require exactly one physical QBit')
        q = arr[0]
        if isinstance(q, bool) or not isinstance(q, Integral) or q < 0:
            raise FredkinError(f'{path}/{name}: invalid physical wire {q!r}')
        wires.append(int(q))
    if len(set(wires)) != 3:
        raise FredkinError(f'{path}: aliased Fredkin wires')
    return dict(zip(('ctrl', 'x', 'y'), wires))


def _leaves(wires, path):
    c, x, y = (wires[name] for name in ('ctrl', 'x', 'y'))
    return [dict(name=name, qubits=qs, path=path, fredkin_rule=RULE_ID, fredkin_step=step)
            for step, (name, qs) in enumerate((('cx', [y, x]),
                                               ('ccx', [c, x, y]),
                                               ('cx', [y, x])))]


class FredkinRules:
    """One fresh book per source graph; chain existing basis/A1 explicitly.

    lower_graph(source, extra_atom_rule=book.chain_with_basis_and_a1)
    then translate(...), then book.bind(low, prog). The caller must still
    require the independently executed rule receipt, exact domain_guard,
    and complete physical input/output proof for the actual context.
    """
    def __init__(self):
        self.basis = BasisVectorRules()
        self.instances = []
        self._paths = set()
        self._bound = False

    def __call__(self, bloq, inputs, allocate, path):
        wires = match(bloq, inputs, path)
        if wires is None:
            return None
        if self._bound or path in self._paths:
            raise FredkinError('FredkinRules may not be reused across source graphs')
        if not isinstance(path, str) or not path:
            raise FredkinError('A nonempty actual source path is required')
        self._paths.add(path)
        self.instances.append(dict(path=path, wires=wires, source_type='TwoBitCSwap',
                                   source_signature='ctrl:QBit, x:QBit, y:QBit; all THRU'))
        # No allocation, retirement, hidden remapping, or source modification.
        return AtomExpansion({name: array.copy() for name, array in inputs.items()},
                             _leaves(wires, path))

    def chain_with_basis_and_a1(self, bloq, inputs, allocate, path):
        result = self(bloq, inputs, allocate, path)
        return result if result is not None else self.basis.chain_with_a1(bloq, inputs, allocate, path)

    def bind(self, low, prog):
        if self._bound:
            raise FredkinError('FredkinRules.bind may be called only once')
        if low.get('source_backend') != 'qualtran_bloqgraph':
            raise FredkinError('Fredkin rule requires the explicit source BloqGraph backend')
        nq, ops = parse_qasm(prog['qasm'])
        if nq != prog['nq'] or ops != prog['ops']:
            raise FredkinError('Literal QASM differs from logical operations')
        layout = prog['layout']
        check_layout(layout, nq)
        if (layout.get('source_nq') != low['nq'] or nq < low['nq']
                or layout['inputs'] != low['inputs'] or layout['outputs'] != low['outputs']
                or layout.get('domain_policy') != 'exact_logical_AND_ERASE_v1'):
            raise FredkinError('Source layout or exact domain policy changed')
        events = {event['path']: event for event in low['graph_events']}
        if len(events) != len(low['graph_events']):
            raise FredkinError('Repeated source event path')
        if {path for path, event in events.items() if event['bloq_type'] == 'TwoBitCSwap'} != self._paths:
            raise FredkinError('Not every source Fredkin atom belongs to this rule book')
        if len(prog['trace']) != len(low['leaves']):
            raise FredkinError('Source/compiled trace length differs')
        cursor = 0
        for index, (leaf, row) in enumerate(zip(low['leaves'], prog['trace'])):
            if (row['source_index'] != index or row['leaf'] != leaf
                    or row['logical_start'] != cursor
                    or not cursor <= row['logical_end'] <= len(ops)):
                raise FredkinError('Noncontiguous or changed actual source trace')
            cursor = row['logical_end']
        if cursor != len(ops):
            raise FredkinError('Trace does not cover all actual operations')
        bound = []
        for item in self.instances:
            path, wires = item['path'], item['wires']
            c, x, y = (wires[name] for name in ('ctrl', 'x', 'y'))
            physical = {name: [q] for name, q in wires.items()}
            if events[path]['inputs_msb'] != physical or events[path]['outputs_msb'] != physical:
                raise FredkinError(f'{path}: actual source atom wiring differs')
            rows = [row for row in prog['trace'] if row['leaf']['path'] == path]
            if [row['leaf'] for row in rows] != _leaves(wires, path):
                raise FredkinError(f'{path}: exact three-leaf expansion changed')
            if ([row['logical_end'] - row['logical_start'] for row in rows] != [1, 3, 1]
                    or rows[0]['logical_end'] != rows[1]['logical_start']
                    or rows[1]['logical_end'] != rows[2]['logical_start']):
                raise FredkinError(f'{path}: noncontiguous Fredkin primitive expansion')
            first = ops[rows[1]['logical_start']]
            if len(first) != 4 or first[:3] != ('and', c, x):
                raise FredkinError(f'{path}: CCX expansion does not start with the expected AND')
            scratch = first[3]
            if (scratch in wires.values() or scratch < low['nq']
                    or scratch not in layout['bridge_scratch'] or scratch not in layout['initial_zero']
                    or any(scratch in word for word in layout['outputs'].values())):
                raise FredkinError(f'{path}: CCX scratch is not separate, initially clean and final-zero')
            start, end = rows[0]['logical_start'], rows[-1]['logical_end']
            expected = [('cx', y, x), ('and', c, x, scratch), ('cx', scratch, y),
                        ('erase', scratch, (('cz', c, x),)), ('cx', y, x)]
            if ops[start:end] != expected:
                raise FredkinError(f'{path}: actual AND/CX/ERASE implementation differs')
            bound.append(dict(item, logical_start=start, logical_end=end, scratch=scratch,
                              required_scratch_initial_and_final=0,
                              exact_partial_domains_required=True, semantic_status='UNPROVED'))
        if bound and not layout.get('partial_gate_domains_required'):
            raise FredkinError('Missing exact partial-domain requirement')
        basis_binding = self.basis.bind(low, prog)
        self._bound = True
        return dict(status='PASS', semantic_status='UNPROVED', rule=RULE_ID,
                    scope='Structural binding only; not quantum or arithmetic proof',
                    instances=deepcopy(bound), basis_vectors=basis_binding,
                    layout_sha256=digest_json(layout), qasm_sha256=sha256(prog['qasm'].encode()).hexdigest(),
                    requires_current_rule_test_receipt=True,
                    requires_exact_partial_domains=True,
                    requires_all_physical_nonoutputs_zero=True)
