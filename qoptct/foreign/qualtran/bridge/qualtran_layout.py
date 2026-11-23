"""A3: exact BloqGraph wiring and total, declarative classical contracts.

This is a separate source backend, NOT a repair of Qualtran's Cirq export.
No circuit simulation, arithmetic proof, synthesis or gate reordering occurs here.
Physical wires are never recycled. A Free creates a final zero obligation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from numbers import Integral
from typing import Any, Callable, Mapping

from .rules import TranslationError


class LayoutError(TranslationError):
    pass


def _integer(value, label, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise LayoutError(f'{label}: expected concrete integer >= {minimum}: {value!r}')
    return int(value)


def _shape(reg):
    shape = tuple(_integer(x, f'{reg.name}.shape', 1) for x in reg.shape)
    bits = _integer(reg.bitsize, f'{reg.name}.bitsize', 1)
    if reg.dtype.num_cbits != 0 or reg.dtype.num_qubits != bits:
        raise LayoutError(f'{reg.name}: only pure quantum registers are supported')
    return shape, bits


def element_name(name, idx):
    # Brackets are reserved. Thus scalar 'a[0]' cannot alias shaped a[0].
    if not isinstance(name, str) or not name or '[' in name or ']' in name:
        raise LayoutError(f'ambiguous register name: {name!r}')
    return name if not idx else name + '[' + ','.join(map(str, idx)) + ']'


def _validated_regs(regs, values, label):
    import numpy as np
    regs = tuple(regs)
    if len({r.name for r in regs}) != len(regs) or set(values) != {r.name for r in regs}:
        raise LayoutError(f'{label}: register names do not match signature')
    result = {}
    for reg in regs:
        shape, bits = _shape(reg)
        element_name(reg.name, ())
        array = np.asarray(values[reg.name], dtype=object)
        if array.shape != shape + (bits,):
            raise LayoutError(f'{label}/{reg.name}: {array.shape} != {shape + (bits,)}')
        for q in array.flat:
            _integer(q, f'{label}/{reg.name}: physical wire')
        result[reg.name] = array.copy()
    flat = [int(q) for arr in result.values() for q in arr.flat]
    if len(flat) != len(set(flat)):
        raise LayoutError(f'{label}: aliased physical wires')
    return result


def flatten_registers(regs, physical):
    """Return scalar element names -> LSB wires, plus lossless shape metadata.

    Qualtran arrays have shape register.shape + (bitsize,) and MSB-first bits.
    We reverse ONLY the last axis, never the order of array elements.
    """
    regs = tuple(regs)
    physical = _validated_regs(regs, physical, 'flatten')
    words, descriptors = {}, []
    for reg in regs:
        shape, bits = _shape(reg)
        elements = []
        for idx in reg.all_idxs():
            name = element_name(reg.name, idx)
            if name in words:
                raise LayoutError(f'duplicate element name: {name}')
            msb = [int(q) for q in physical[reg.name][idx]]
            words[name] = list(reversed(msb))
            elements.append({'name': name, 'index': list(idx), 'wires_msb': msb})
        descriptors.append({'name': reg.name, 'shape': list(shape), 'bitsize': bits,
                            'dtype': str(reg.dtype), 'side': reg.side.name,
                            'elements': elements})
    return words, descriptors


@dataclass
class AtomExpansion:
    """Trusted atom-rule result. This interface does not establish rule correctness."""
    outputs: Mapping[str, Any]
    leaves: list[dict] = field(default_factory=list)
    retired: tuple[int, ...] = ()


def controlled_and_atom_rule(bloq, inputs, allocate, path):
    """Adapter to the separately audited A1 module; opt in via extra_atom_rule.

    Import errors are explicit failures. This never silently falls back to an
    unverified coherent-control decomposition of measurement uncomputation.
    """
    from qualtran import Side
    from .controlled_and import match_gate

    # First probe the exact matcher using placeholder distinct physical IDs.
    regs = tuple(bloq.signature)
    size = 0
    for reg in regs:
        shape, bits = _shape(reg)
        count = bits
        for dim in shape:
            count *= dim
        size += count
    if match_gate(bloq, list(range(size)), path) is None:
        return None
    only_rights = [r for r in regs if r.side == Side.RIGHT]
    allocated = allocate(only_rights, path)
    all_regs = dict(inputs, **allocated)
    qubits = [int(q) for reg in regs for q in all_regs[reg.name].flat]
    leaf = match_gate(bloq, qubits, path)
    if leaf is None:
        raise LayoutError(f'{path}: inconsistent controlled-And matcher')
    outputs = {r.name: all_regs[r.name] for r in bloq.signature.rights()}
    consumed = tuple(int(q) for r in regs if r.side == Side.LEFT for q in inputs[r.name].flat)
    return AtomExpansion(outputs, [leaf], consumed)


def basic_atom_rule(bloq, inputs, allocate, path):
    """Only exact registered types. Unknown subclasses must decompose or fail."""
    import numpy as np
    from qualtran.bloqs.basic_gates import XGate, CNOT, Toffoli
    from qualtran.bloqs.bookkeeping import Allocate, Free, Split, Join, Cast, Partition
    from qualtran.bloqs.mcmt import And

    typ = type(bloq)
    rights = tuple(bloq.signature.rights())
    if typ is Allocate:
        if bloq.dirty:
            raise LayoutError(f'{path}: dirty Allocate is forbidden')
        return AtomExpansion(allocate(rights, path))
    if typ is Free:
        if bloq.dirty:
            raise LayoutError(f'{path}: dirty Free is forbidden')
        return AtomExpansion({}, retired=tuple(int(q) for q in inputs['reg'].flat))
    if typ in (Split, Join, Cast):
        reg, = rights
        shape, bits = _shape(reg)
        return AtomExpansion({reg.name: inputs['reg'].reshape(shape + (bits,))})
    if typ is Partition:
        # Signature order is the source's tensor order; do not sort register names.
        bits = [q for reg in bloq.signature.lefts() for q in inputs[reg.name].flat]
        out, start = {}, 0
        for reg in rights:
            shape, width = _shape(reg)
            size = width
            for dim in shape:
                size *= dim
            out[reg.name] = np.asarray(bits[start:start + size], dtype=object).reshape(shape + (width,))
            start += size
        if start != len(bits):
            raise LayoutError(f'{path}: Partition changed bit count')
        return AtomExpansion(out)
    if typ in (XGate, CNOT, Toffoli):
        # The declared signature order is also the official atom's qubit order.
        qubits = [int(q) for reg in bloq.signature.lefts() for q in inputs[reg.name].flat]
        name = {XGate: 'x', CNOT: 'cx', Toffoli: 'ccx'}[typ]
        expected = {'x': 1, 'cx': 2, 'ccx': 3}[name]
        if len(qubits) != expected:
            raise LayoutError(f'{path}: atom signature changed')
        return AtomExpansion(inputs, [{'name': name, 'qubits': qubits, 'path': path}])
    if typ is And:
        cvs = [_integer(v, f'{path}: control value') for v in (bloq.cv1, bloq.cv2)]
        if any(v not in (0, 1) for v in cvs):
            raise LayoutError(f'{path}: nonbinary And controls')
        ctrl = inputs['ctrl'].copy()
        if ctrl.shape != (2, 1):
            raise LayoutError(f'{path}: And signature changed')
        if bloq.uncompute:
            target = inputs['target']
            output = {'ctrl': ctrl}
        else:
            target = allocate([r for r in rights if r.name == 'target'], path)['target']
            output = {'ctrl': ctrl, 'target': target}
        qubits = [int(q) for q in ctrl.flat] + [int(q) for q in target.flat]
        if len(qubits) != 3:
            raise LayoutError(f'{path}: And target is not one qubit')
        return AtomExpansion(output, [{'name': 'erase' if bloq.uncompute else 'and',
                                     'qubits': qubits, 'controls': cvs, 'path': path}],
                             (qubits[-1],) if bloq.uncompute else ())
    return None


def lower_bloqgraph(bloq, *, extra_atom_rule: Callable | None = None,
                    max_leaves=200000, max_nodes=1000000, max_depth=128, record_events=True):
    """Lower actual soquet connections, preserving each unflattened builder order.

    A CompositeBloq whose instance indices are not a topological order is rejected.
    No optimizer/greedy scheduling is used. This intentionally refuses flattened
    graphs whose renumbering has destroyed the original construction order.
    All allocations use fresh IDs; final-zero obligations include every dead wire.
    """
    import numpy as np
    from qualtran import CompositeBloq, LeftDangle, RightDangle, Soquet, Side

    max_leaves = _integer(max_leaves, 'max_leaves', 1)
    max_nodes = _integer(max_nodes, 'max_nodes', 1)
    max_depth = _integer(max_depth, 'max_depth', 1)
    wire_names, leaves, events, permutations, retired = [], [], [], [], set()
    node_count = 0

    def allocate(regs, path):
        out = {}
        for reg in regs:
            shape, bits = _shape(reg)
            array = np.empty(shape + (bits,), dtype=object)
            for idx in reg.all_idxs():
                for bit in range(bits):
                    array[idx + (bit,)] = len(wire_names)
                    wire_names.append(f'{path}/{element_name(reg.name, idx)}:msb{bit}')
            out[reg.name] = array
        return out

    def expand(current, inputs, path, depth):
        nonlocal node_count
        node_count += 1
        if node_count > max_nodes or depth > max_depth:
            raise LayoutError(f'{path}: graph size/depth limit exceeded')
        inputs = _validated_regs(current.signature.lefts(), inputs, path + ':inputs')
        in_set = {int(q) for arr in inputs.values() for q in arr.flat}
        if in_set & retired:
            raise LayoutError(f'{path}: re-use of retired physical wire')
        first_new = len(wire_names)
        expansion = basic_atom_rule(current, inputs, allocate, path)
        if expansion is None and extra_atom_rule is not None:
            expansion = extra_atom_rule(current, inputs, allocate, path)
            if expansion is not None and not isinstance(expansion, AtomExpansion):
                raise LayoutError(f'{path}: atom callback must return AtomExpansion or None')
        if expansion is None:
            try:
                cb = current if isinstance(current, CompositeBloq) else current.decompose_bloq()
            except Exception as exc:
                raise LayoutError(f'{path}: no supported atom or source BloqGraph: {exc}') from exc
            if not isinstance(cb, CompositeBloq) or cb.signature != current.signature:
                raise LayoutError(f'{path}: decomposition changed signature')
            incoming = {}
            for cxn in cb.connections:
                if cxn.right in incoming:
                    raise LayoutError(f'{path}: multiple connections into one soquet')
                incoming[cxn.right] = cxn.left
            values = {}
            for reg in cb.signature.lefts():
                for idx in reg.all_idxs():
                    values[Soquet(LeftDangle, reg, idx)] = tuple(inputs[reg.name][idx])

            def receive(owner, regs):
                out = {}
                for reg in regs:
                    shape, bits = _shape(reg)
                    array = np.empty(shape + (bits,), dtype=object)
                    for idx in reg.all_idxs():
                        destination = Soquet(owner, reg, idx)
                        if destination not in incoming:
                            raise LayoutError(f'{path}: missing input connection {destination}')
                        predecessor = incoming.pop(destination)
                        if predecessor not in values:
                            raise LayoutError(f'{path}: reused soquet or non-topological instance order: {predecessor}')
                        wires = values.pop(predecessor)
                        if len(wires) != bits:
                            raise LayoutError(f'{path}: connection changes element width')
                        array[idx] = wires
                    out[reg.name] = array
                return out

            instances = sorted(cb.bloq_instances, key=lambda b: b.i)
            ids = [_integer(b.i, f'{path}: instance index') for b in instances]
            if len(ids) != len(set(ids)):
                raise LayoutError(f'{path}: duplicate instance index')
            for binst in instances:
                subinputs = receive(binst, binst.bloq.signature.lefts())
                suboutputs = expand(binst.bloq, subinputs, f'{path}/{binst.i}:{type(binst.bloq).__name__}', depth + 1)
                for reg in binst.bloq.signature.rights():
                    for idx in reg.all_idxs():
                        values[Soquet(binst, reg, idx)] = tuple(suboutputs[reg.name][idx])
            outputs = receive(RightDangle, cb.signature.rights())
            if values or incoming:
                raise LayoutError(f'{path}: unused or dangling soquets')
        else:
            outputs = expansion.outputs
            legal = in_set | set(range(first_new, len(wire_names)))
            for leaf in expansion.leaves:
                qs = leaf.get('qubits', [])
                if (not qs or len(set(qs)) != len(qs)
                        or any(type(q) is not int or q not in legal for q in qs)):
                    raise LayoutError(f'{path}: atom emitted illegal physical wires')
                leaves.append(dict(leaf, path=leaf.get('path', path)))
                if len(leaves) > max_leaves:
                    raise LayoutError(f'{path}: too many leaves')
            new_retired = {_integer(q, f'{path}: retired wire') for q in expansion.retired}
            if not new_retired <= legal or new_retired & retired:
                raise LayoutError(f'{path}: invalid retired wires')
            retired.update(new_retired)
        outputs = _validated_regs(current.signature.rights(), outputs, path + ':outputs')
        out_set = {int(q) for arr in outputs.values() for q in arr.flat}
        legal = in_set | set(range(first_new, len(wire_names)))
        if not out_set <= legal or out_set & retired:
            raise LayoutError(f'{path}: output uses nonexistent or retired physical wire')
        missing = legal - out_set - retired
        if missing:
            raise LayoutError(f'{path}: atom/decomposition silently drops physical wires {sorted(missing)}')
        for reg in current.signature:
            if reg.side == Side.THRU:
                before = [int(q) for q in inputs[reg.name].flat]
                after = [int(q) for q in outputs[reg.name].flat]
                if before != after:
                    positions = {q: i for i, q in enumerate(before)}
                    permutations.append({'path': path, 'bloq_type': type(current).__name__,
                                         'register': reg.name, 'shape': list(reg.shape),
                                         'before_msb': before, 'after_msb': after,
                                         'permutation': [positions.get(q) for q in after]})
        if record_events is True or (record_events == 'top' and path.count('/') == 1):  # 'top': record only the top level
          events.append({'path': path, 'bloq_type': type(current).__name__,
                       'inputs_msb': {k: [int(q) for q in v.flat] for k, v in inputs.items()},
                       'outputs_msb': {k: [int(q) for q in v.flat] for k, v in outputs.items()}})
        return outputs

    input_physical = allocate(bloq.signature.lefts(), 'input')
    output_physical = expand(bloq, input_physical, type(bloq).__name__, 0)
    inputs, input_desc = flatten_registers(bloq.signature.lefts(), input_physical)
    outputs, output_desc = flatten_registers(bloq.signature.rights(), output_physical)
    input_wires = {q for word in inputs.values() for q in word}
    output_wires = {q for word in outputs.values() for q in word}
    return {'nq': len(wire_names), 'leaves': leaves, 'inputs': inputs, 'outputs': outputs,
            'wire_names': wire_names, 'bloq': bloq, 'source_backend': 'qualtran_bloqgraph',
            'source_semantics': 'actual decompose_bloq connections; not Cirq export',
            'registers': {'inputs': input_desc, 'outputs': output_desc},
            'initial_zero': sorted(set(range(len(wire_names))) - input_wires),
            'final_zero': sorted(set(range(len(wire_names))) - output_wires),
            'retired_wires': sorted(retired), 'garbage': [], 'free_wires': [],
            'thru_permutations': permutations, 'graph_events': events,
            'operation_order': 'recursive source BloqInstance.i; reject non-topological indices',
            'allocation_policy': 'fresh physical wire for each allocation; never recycle'}


def require_cirq_thru_identity(lowered):
    """Structural guard for the OLD fixed-wire Cirq decomposition path.

    A nonidentity THRU mapping needs an explicit exporter correspondence proof.
    The legacy bridge has none, so it must fail before an arithmetic prover.
    This conservative gate does not certify phase or other export semantics.
    """
    if lowered.get('source_backend') != 'qualtran_bloqgraph':
        raise LayoutError('THRU audit requires a real BloqGraph result')
    failures = lowered.get('thru_permutations')
    if not isinstance(failures, list):
        raise LayoutError('missing structural THRU audit')
    if failures:
        f = failures[0]
        raise LayoutError(f"THRU permutation lost by unaudited fixed-wire Cirq path: {f['path']}/"
                          f"{f['register']}: {f['before_msb']} -> {f['after_msb']}")
    return {'status': 'PASS', 'check': 'THRU physical identity only', 'arithmetic_prover_used': False}


# Expressions below are an independent, serializable integer specification.
# They never refer to source gates, source call_classically, or simulated outputs.
def ref(name):
    return ('ref', name)


def expr(op, *args):
    return (op, *args)


_ARITY = {'ref': 1, 'add': 2, 'sub': 2, 'mul': 2, 'mod': 2, 'div': 2,
          'xor': 2, 'and': 2, 'or': 2, 'shl': 2, 'shr': 2,
          'lt': 2, 'eq': 2, 'select': 3}


def _check_expr(node, names):
    if type(node) is int:
        return
    if not isinstance(node, (tuple, list)) or not node or node[0] not in _ARITY:
        raise LayoutError(f'not an explicit supported expression: {node!r}')
    op, *args = node
    if len(args) != _ARITY[op]:
        raise LayoutError(f'{op}: wrong expression arity')
    if op == 'ref':
        if not isinstance(args[0], str) or args[0] not in names:
            raise LayoutError(f'unknown/forward expression reference: {args[0]!r}')
    else:
        for arg in args:
            _check_expr(arg, names)


def _eval_expr(node, values):
    if type(node) is int:
        return node
    op, *args = node
    if op == 'ref':
        return values[args[0]]
    if op == 'select':
        condition = _eval_expr(args[0], values)
        if condition not in (0, 1):
            raise LayoutError('select condition must be Boolean')
        return _eval_expr(args[1] if condition else args[2], values)
    a, b = [_eval_expr(arg, values) for arg in args]
    if op in ('div', 'mod') and b <= 0:
        raise LayoutError('division/modulus requires positive divisor')
    if op in ('shl', 'shr') and b < 0:
        raise LayoutError('negative shift')
    if op == 'add': return a + b
    if op == 'sub': return a - b
    if op == 'mul': return a * b
    if op == 'div': return a // b
    if op == 'mod': return a % b
    if op == 'xor': return a ^ b
    if op == 'and': return a & b
    if op == 'or': return a | b
    if op == 'shl': return a << b
    if op == 'shr': return a >> b
    if op == 'lt': return int(a < b)
    if op == 'eq': return int(a == b)
    raise LayoutError(f'unknown expression operation {op}')


@dataclass
class ClassicalContract:
    """Every output has an expression. No existential/ignored junk is possible."""
    input_widths: dict[str, int]
    output_widths: dict[str, int]
    outputs: dict[str, Any]
    statements: list[tuple[str, Any]] = field(default_factory=list)
    precondition: Any = 1
    description: str = ''

    def validate(self):
        for side in (self.input_widths, self.output_widths):
            for name, width in side.items():
                if not isinstance(name, str) or not name:
                    raise LayoutError('contract variable names must be nonempty strings')
                _integer(width, f'contract width {name}', 1)
        if set(self.outputs) != set(self.output_widths):
            raise LayoutError('missing/extra output function; arbitrary garbage is forbidden')
        names = set(self.input_widths)
        _check_expr(self.precondition, names)
        for name, value in self.statements:
            if not isinstance(name, str) or not name or name in names:
                raise LayoutError(f'duplicate/invalid temporary name: {name!r}')
            _check_expr(value, names)
            names.add(name)
        for value in self.outputs.values():
            _check_expr(value, names)
        return self

    def evaluate(self, inputs):
        self.validate()
        if set(inputs) != set(self.input_widths):
            raise LayoutError('input values must cover the entire input signature')
        for name, value in inputs.items():
            if type(value) is not int or not 0 <= value < 1 << self.input_widths[name]:
                raise LayoutError(f'input {name}: outside unsigned bit-vector domain')
        values = dict(inputs)
        pre = _eval_expr(self.precondition, values)
        if pre not in (0, 1):
            raise LayoutError('precondition must evaluate to a Boolean')
        if not pre:
            raise LayoutError('input violates the explicit contract precondition')
        for name, value in self.statements:
            values[name] = _eval_expr(value, values)
        result = {name: _eval_expr(value, values) for name, value in self.outputs.items()}
        for name, value in result.items():
            if not 0 <= value < 1 << self.output_widths[name]:
                raise LayoutError(f'output function {name} overflows; truncation must be explicit')
        return result

    def bind(self, lowered):
        self.validate()
        _integer(lowered['nq'], 'nq')
        if lowered.get('garbage') or lowered.get('free_wires'):
            raise LayoutError('arbitrary garbage/free_wires are forbidden')
        for field_name, widths in (('inputs', self.input_widths), ('outputs', self.output_widths)):
            words = lowered[field_name]
            if set(widths) != set(words) or any(widths[k] != len(words[k]) for k in words):
                raise LayoutError(f'{field_name}: contract does not cover the exact shaped signature')
            flat = [q for word in words.values() for q in word]
            if (len(flat) != len(set(flat))
                    or any(type(q) is not int or not 0 <= q < lowered['nq'] for q in flat)):
                raise LayoutError(f'{field_name}: invalid or aliased physical layout')
        inputs = {q for word in lowered['inputs'].values() for q in word}
        outputs = {q for word in lowered['outputs'].values() for q in word}
        all_wires = set(range(lowered['nq']))
        for field_name in ('initial_zero', 'final_zero'):
            zeros = lowered.get(field_name)
            if (not isinstance(zeros, (tuple, list)) or len(zeros) != len(set(zeros))
                    or any(type(q) is not int or q not in all_wires for q in zeros)):
                raise LayoutError(f'{field_name}: invalid full clean obligation')
        if set(lowered.get('initial_zero', ())) != all_wires - inputs:
            raise LayoutError('missing full clean-input obligation')
        if set(lowered.get('final_zero', ())) != all_wires - outputs:
            raise LayoutError('missing full clean-output obligation')
        return self

    def expected_physical_bits(self, lowered, inputs):
        self.bind(lowered)
        values = self.evaluate(inputs)
        bits = [0] * lowered['nq']
        for name, word in lowered['outputs'].items():
            for i, wire in enumerate(word):
                bits[wire] = (values[name] >> i) & 1
        return bits

    def to_dict(self):
        self.validate()
        return {'format': 'a3_classical_contract_v1', 'description': self.description,
                'input_widths': self.input_widths, 'output_widths': self.output_widths,
                'precondition': self.precondition, 'statements': self.statements,
                'outputs': self.outputs, 'junk_policy': 'every output explicit; all other wires zero'}


def montgomery_junk_contract(n, window, modulus):
    """Independent arithmetic definition for DirtyOutOfPlaceMontgomeryModMul.

    x,y are canonical Montgomery encodings in [0,p). Each window stores the low
    w bits BEFORE reduction. Its correction is -p^{-1} times those bits mod 2^w.
    The first window is the most significant word of qrom_indices.
    This contract does not call or inspect the candidate's classical action.
    """
    n = _integer(n, 'n', 1)
    w = _integer(window, 'window', 1)
    p = _integer(modulus, 'modulus', 3)
    if p % 2 != 1 or w > n or n % w or 2 * p - 1 >= 1 << n:
        raise LayoutError('require odd p, w|n, and 2p-1 < 2^n')
    base = 1 << w
    # Python modular inverse is specification arithmetic, not a source gate rule.
    minus_inverse = (-pow(p, -1, base)) % base
    statements = [('t0', 0), ('h0', 0)]
    for i in range(n // w):
        t, h = f't{i}', f'h{i}'
        u, m, correction = f'u{i}', f'm{i}', f'k{i}'
        digit = expr('mod', expr('shr', ref('x'), i * w), base)
        statements.extend([
            (u, expr('add', ref(t), expr('mul', digit, ref('y')))),
            (m, expr('mod', ref(u), base)),
            (correction, expr('mod', expr('mul', ref(m), minus_inverse), base)),
            (f't{i + 1}', expr('div', expr('add', ref(u), expr('mul', ref(correction), p)), base)),
            (f'h{i + 1}', expr('add', expr('mul', ref(h), base), ref(m))),
        ])
    t, h = ref(f't{n // w}'), ref(f'h{n // w}')
    reduced = expr('sub', 1, expr('lt', t, p))
    return ClassicalContract(
        {'x': n, 'y': n},
        {'x': n, 'y': n, 'target': n, 'qrom_indices': n, 'reduced': 1},
        {'x': ref('x'), 'y': ref('y'), 'target': expr('sub', t, expr('mul', reduced, p)),
         'qrom_indices': h, 'reduced': reduced}, statements,
        expr('and', expr('lt', ref('x'), p), expr('lt', ref('y'), p)),
        'Canonical x,y<p. Montgomery product; explicit window-history h and final reduction flag c.',
    ).validate()
