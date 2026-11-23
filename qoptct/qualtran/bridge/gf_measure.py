"""A2: the exact flow of the measurement-based uncomputation of Qualtran GF2MulMBUC.

No synthesis from the arithmetic meaning: lower_mbuc reads the actual
CompositeBloq. Every source operation has a trace record. A classical Parity
is substituted as an unchanged XOR expression into the later conditions; quantum
operations are not moved. This module changes none of the common bridge modules.
"""
from __future__ import annotations

from functools import reduce
import math
import operator
import re


class GFMeasureError(ValueError):
    pass


def lower_mbuc(bloq):
    """Return nq, ops, trace, layout and the actual CompositeBloq in _source.

    ops: measure_x(q,c), parity((c,...),new_c), classical_cz(c,a,b,cv),
    discard_classical(c). Quantum numbers are the LSB of the integer state; the register
    lists of the layout are LSB-first. A Soquet inside Qualtran stores MSB-first.
    """
    import numpy as np
    from qualtran import CBit, LeftDangle, Soquet
    from qualtran.bloqs.basic_gates import CZ, Discard, MeasureX
    from qualtran.bloqs.bookkeeping import Join, Split
    from qualtran.bloqs.gf_arithmetic.gf2_multiplication import GF2MulMBUC, Parity
    from qualtran.bloqs.mcmt.classically_controlled import ClassicallyControlled

    if type(bloq) is not GF2MulMBUC:
        raise GFMeasureError('The original GF2MulMBUC is required')
    cb = bloq.decompose_bloq()
    n = int(bloq.bitsize)
    mapping, inputs, offset = {}, {}, 0
    for reg in cb.signature.lefts():
        if reg.shape or int(reg.bitsize) != n:
            raise GFMeasureError('Unexpected input signature of GF2MulMBUC')
        word = tuple(range(offset, offset + n))
        inputs[reg.name] = list(word)
        mapping[Soquet(LeftDangle, reg)] = ('q', tuple(reversed(word)))
        offset += n
    if set(inputs) != {'x', 'y', 'result'}:
        raise GFMeasureError('The signature of the original GF2MulMBUC changed')
    ops, trace, next_c = [], [], 0

    def mapped(soq):
        if isinstance(soq, np.ndarray):
            return [mapped(s) for s in soq.flat]
        try:
            return mapping[soq]
        except KeyError as exc:
            raise GFMeasureError('The source used a soquet that is not defined yet') from exc

    def qbit(value):
        if not isinstance(value, tuple) or value[0] != 'q' or len(value[1]) != 1:
            raise GFMeasureError('One quantum bit expected')
        return value[1][0]

    def cbit(value):
        if not isinstance(value, tuple) or value[0] != 'c':
            raise GFMeasureError('One classical bit expected')
        return value[1]

    def bind(soqs, value):
        if isinstance(soqs, np.ndarray):
            if len(soqs.flat) != len(value):
                raise GFMeasureError('The size of the output array changed')
            for s, v in zip(soqs.flat, value):
                mapping[s] = v
        else:
            mapping[soqs] = value

    # iter_bloqsoqs gives a valid topological sort, which may reorder
    # independent operations. i keeps the bb.add order of the source.
    rows = sorted(cb.iter_bloqsoqs(), key=lambda row: row[0].i)
    for binst, in_soqs, out_soqs in rows:
        b = binst.bloq
        values = {key: mapped(value) for key, value in in_soqs.items()}
        start = len(ops)
        if type(b) is Split:
            ref = values['reg']
            if ref[0] != 'q' or len(ref[1]) != int(b.dtype.num_qubits):
                raise GFMeasureError('Wrong quantum layout of Split')
            produced = {'reg': [('q', (q,)) for q in ref[1]]}
            rule = 'split_layout_only'
        elif type(b) is Join:
            bits = tuple(qbit(ref) for ref in values['reg'])
            if len(bits) != int(b.dtype.num_qubits):
                raise GFMeasureError('Wrong quantum layout of Join')
            produced = {'reg': ('q', bits)}
            rule = 'join_layout_only'
        elif type(b) is MeasureX:
            cid = next_c
            next_c += 1
            ops.append(('measure_x', qbit(values['q']), cid))
            produced = {'c': ('c', cid)}
            rule = 'measure_x_reset_consumed_wire'
        elif type(b) is Parity:
            ids = tuple(cbit(ref) for ref in values['x'])
            if len(ids) != int(b.n):
                raise GFMeasureError('Wrong length of the source Parity')
            cid = next_c
            next_c += 1
            ops.append(('parity', ids, cid))
            produced = {'x': values['x'], 'parity': ('c', cid)}
            rule = 'classical_xor_expression'
        elif type(b) is ClassicallyControlled and type(b.subbloq) is CZ:
            spec = b.ctrl_spec
            if (len(spec.qdtypes) != 1 or type(spec.qdtypes[0]) is not CBit
                    or tuple(spec.concrete_shapes) != ((),)):
                raise GFMeasureError('Only one scalar classical control is allowed')
            cv = int(np.asarray(spec.cvs[0]))
            if cv not in (0, 1) or len(b.ctrl_reg_names) != 1:
                raise GFMeasureError('Unsupported classical control')
            ctrl = cbit(values[b.ctrl_reg_names[0]])
            ops.append(('classical_cz', ctrl, qbit(values['q1']), qbit(values['q2']), cv))
            produced = values
            rule = 'classical_conditional_cz'
        elif type(b) is Discard:
            ops.append(('discard_classical', cbit(values['c'])))
            produced = {}
            rule = 'discard_classical_trace'
        else:
            raise GFMeasureError(f'binst={binst.i}: no exact rule {type(b).__name__}: {b}')
        registers = tuple(b.signature.rights())
        if len(registers) != len(out_soqs) or set(produced) != {r.name for r in registers}:
            raise GFMeasureError('The rule does not cover all outputs of the source operation')
        for reg, soq in zip(registers, out_soqs):
            bind(soq, produced[reg.name])
        trace.append({'binst_i': binst.i, 'source_type': type(b).__name__,
                      'source_repr': str(b), 'rule': rule,
                      'op_start': start, 'op_stop': len(ops)})
    outputs = {}
    for key, soq in cb.final_soqs().items():
        ref = mapped(soq)
        if not isinstance(ref, tuple) or ref[0] != 'q':
            raise GFMeasureError('Unexpected classical output of GF2MulMBUC')
        outputs[key] = list(reversed(ref[1]))
    if outputs != {key: inputs[key] for key in ('x', 'y')}:
        raise GFMeasureError('The output layout of the source did not keep x,y')
    program = {'nq': offset, 'nc_ssa': next_c, 'ops': ops, 'trace': trace,
               'layout': {'nq': offset, 'inputs': inputs, 'outputs': outputs,
                          'garbage': [],
                          'domain_policy': 'total_measurement_instruments_v1',
                          'partial_gate_domains_required': False,
                          'source_backend': 'qualtran_bloqgraph_measurement',
                          'consumed_quantum': inputs['result'],
                          'initial_zero': [], 'final_zero': inputs['result'],
                          'bit_order': 'LSB-first'},
               'n': n, 'field_poly': int(bloq.qgf.gf_type.irreducible_poly),
               '_source': cb, '_bloq': bloq}
    validate_ir(program)
    return program


def validate_ir(program):
    nq = program['nq']
    if type(nq) is not int or nq <= 0:
        raise GFMeasureError('Invalid number of qubits')
    live, defined, consumed = set(), set(), set()
    for index, op in enumerate(program['ops']):
        kind = op[0]

        def qubits(*qs):
            if len(set(qs)) != len(qs) or any(type(q) is not int or not 0 <= q < nq or q in consumed for q in qs):
                raise GFMeasureError(f'Operation {index}: wrong or already consumed qubit')

        def new(cid):
            if type(cid) is not int or cid < 0 or cid in defined:
                raise GFMeasureError('Classical values must have unique SSA numbers')
            defined.add(cid)
            live.add(cid)

        def reads(*ids):
            if any(type(cid) is not int or cid not in live for cid in ids):
                raise GFMeasureError('Read of an undefined or already discarded CBit')

        if kind == 'measure_x' and len(op) == 3:
            _, q, cid = op
            qubits(q)
            new(cid)
            consumed.add(q)
        elif kind == 'parity' and len(op) == 3:
            _, ids, cid = op
            reads(*ids)
            new(cid)
        elif kind == 'classical_cz' and len(op) == 5:
            _, cid, a, b, cv = op
            reads(cid)
            qubits(a, b)
            if type(cv) is not int or cv not in (0, 1):
                raise GFMeasureError('The control value must be 0 or 1')
        elif kind == 'discard_classical' and len(op) == 2:
            reads(op[1])
            live.remove(op[1])
        else:
            raise GFMeasureError(f'Unknown A2 operation {op}')
    if live:
        raise GFMeasureError('Undeclared live classical outputs')
    return True


def emit_qasm(program):
    """Exact QASM for the unchanged common kprove, not for the old erase parser.

    Parity and Discard are operations of the classical graph. Their record stays
    in a comment and in the trace; no new quantum gates are created for them.
    """
    validate_ir(program)
    measurements = sum(op[0] == 'measure_x' for op in program['ops'])
    lines = ['OPENQASM 3;', 'include "stdgates.inc";',
             f'qubit[{program["nq"]}] q;', f'bit[{max(1, measurements)}] c;']
    expressions, next_measure = {}, 0
    for index, op in enumerate(program['ops']):
        kind = op[0]
        lines.append(f'// A2 op {index}: {kind}')
        if kind == 'measure_x':
            _, q, cid = op
            expressions[cid] = f'c[{next_measure}]'
            lines.extend((f'h q[{q}];', f'c[{next_measure}] = measure q[{q}];',
                          f'if (c[{next_measure}]) {{x q[{q}];}}'))
            next_measure += 1
        elif kind == 'parity':
            _, ids, cid = op
            # No sorting, no cancelling of equal terms and no recognition of a
            # polynomial: exactly the XOR of the inputs of the source Parity in the source order.
            expressions[cid] = '(' + ' ^ '.join(expressions[x] for x in ids) + ')' if ids else '0'
            lines.append(f'// SSA {cid} := {expressions[cid]}')
        elif kind == 'classical_cz':
            _, cid, a, b, cv = op
            condition = expressions[cid] if cv == 1 else f'!({expressions[cid]})'
            lines.append(f'if ({condition}) {{cz q[{a}], q[{b}];}}')
        elif kind == 'discard_classical':
            del expressions[op[1]]
            lines.append(f'// discard SSA {op[1]}: no quantum operation')
    return '\n'.join(lines) + '\n'


def basis_branch(program, initial, outcomes):
    """Independent exact basis branch: (bits, sign, number_measurements)."""
    validate_ir(program)
    if type(initial) is not int or not 0 <= initial < 1 << program['nq']:
        raise GFMeasureError('Initial state out of range')
    value, parity, klass, cursor = initial, 0, {}, 0
    for op in program['ops']:
        if op[0] == 'measure_x':
            _, q, cid = op
            if cursor >= len(outcomes) or outcomes[cursor] not in (0, 1):
                raise GFMeasureError('Every binary measurement outcome is required')
            m = outcomes[cursor]
            cursor += 1
            parity ^= m & ((value >> q) & 1)
            value &= ~(1 << q)
            klass[cid] = m
        elif op[0] == 'parity':
            klass[op[2]] = reduce(operator.xor, (klass[c] for c in op[1]), 0)
        elif op[0] == 'classical_cz':
            _, cid, a, b, cv = op
            if klass[cid] == cv:
                parity ^= ((value >> a) & 1) & ((value >> b) & 1)
        elif op[0] == 'discard_classical':
            del klass[op[1]]
    if cursor != len(outcomes):
        raise GFMeasureError('Surplus measurement outcomes')
    return value, (-1 if parity else 1), cursor


def _condition(text):
    """Own restricted XOR/NOT parser; qoptct is not imported."""
    tokens = re.findall(r'c\[\d+\]|[()!^01]', text)
    if ''.join(tokens) != re.sub(r'\s+', '', text):
        raise GFMeasureError('Condition outside the A2 language: ' + text)
    position = 0

    def atom():
        nonlocal position
        if position >= len(tokens):
            raise GFMeasureError('Unfinished classical condition')
        token = tokens[position]
        position += 1
        if token == '!':
            return ('not', atom())
        if token == '(':
            result = expr()
            if position >= len(tokens) or tokens[position] != ')':
                raise GFMeasureError('Unclosed bracket of the condition')
            position += 1
            return result
        if token in ('0', '1'):
            return ('const', int(token))
        if token.startswith('c['):
            return ('bit', int(token[2:-1]))
        raise GFMeasureError('Wrong atom of the condition')

    def expr():
        nonlocal position
        result = atom()
        while position < len(tokens) and tokens[position] == '^':
            position += 1
            result = ('xor', result, atom())
        return result

    tree = expr()
    if position != len(tokens):
        raise GFMeasureError('Surplus text of the condition')
    return tree


def _condition_bits(tree):
    if tree[0] == 'bit':
        return {tree[1]}
    if tree[0] == 'const':
        return set()
    return set().union(*(_condition_bits(t) for t in tree[1:]))


def _eval_condition(tree, values):
    kind = tree[0]
    if kind == 'bit':
        return values[tree[1]]
    if kind == 'const':
        return tree[1]
    if kind == 'not':
        return 1 ^ _eval_condition(tree[1], values)
    return _eval_condition(tree[1], values) ^ _eval_condition(tree[2], values)


def parse_qasm(text, *, allow_canary_gates=False):
    """Independent strict parser of exactly the A2 QASM actually emitted."""
    lines = [line.split('//', 1)[0].strip() for line in text.splitlines()]
    lines = [line for line in lines if line]
    if lines[:2] != ['OPENQASM 3;', 'include "stdgates.inc";'] or len(lines) < 4:
        raise GFMeasureError('Unknown header of the A2 QASM')
    qm = re.fullmatch(r'qubit\[(\d+)\] q;', lines[2])
    cm = re.fullmatch(r'bit\[(\d+)\] c;', lines[3])
    if qm is None or cm is None:
        raise GFMeasureError('No explicit registers q and c')
    nq, nc = int(qm[1]), int(cm[1])
    if nq <= 0 or nc <= 0:
        raise GFMeasureError('Empty physical register')
    result, measured = [], set()

    def gate(statement):
        names = r'h|x|cz|z|ccx' if allow_canary_gates else r'h|x|cz'
        match = re.fullmatch(r'(' + names + r') (q\[\d+\](?:, q\[\d+\])*);', statement)
        if match is None:
            raise GFMeasureError('Gate outside the A2 language: ' + statement)
        qs = tuple(map(int, re.findall(r'q\[(\d+)\]', match[2])))
        arity = {'h': 1, 'x': 1, 'z': 1, 'cz': 2, 'ccx': 3}
        if len(qs) != arity[match[1]] or len(set(qs)) != len(qs) or any(q >= nq for q in qs):
            raise GFMeasureError('Wrong gate arguments')
        return (match[1], *qs)

    for line in lines[4:]:
        measurement = re.fullmatch(r'c\[(\d+)\] = measure q\[(\d+)\];', line)
        condition = re.fullmatch(r'if \((.*)\) \{(.*)\}', line)
        if measurement is not None:
            c, q = map(int, measurement.groups())
            if c >= nc or q >= nq or c in measured:
                raise GFMeasureError('Measurement out of range or repeated write of a CBit')
            measured.add(c)
            result.append(('measure', q, c))
        elif condition is not None:
            tree = _condition(condition[1])
            if not _condition_bits(tree) <= measured:
                raise GFMeasureError('The condition reads a bit that is not measured yet')
            body = gate(condition[2])
            if body[0] not in ('x', 'cz'):
                raise GFMeasureError('Unsupported conditional gate')
            result.append(('if', tree, body))
        else:
            result.append(gate(line))
    return nq, nc, result


def vector_branch(qasm, state, outcomes, *, allow_canary_gates=False):
    """Non-normalized branch of the actual QASM. Norm squared is the probability of the branch."""
    import numpy as np
    nq, _, operations = parse_qasm(qasm, allow_canary_gates=allow_canary_gates)
    vector = np.array(state, dtype=np.complex128, copy=True)
    if vector.shape != (1 << nq,):
        raise GFMeasureError('Wrong length of the source vector')
    indices = np.arange(len(vector))
    classical, cursor = {}, 0

    def apply_gate(op):
        nonlocal vector
        if op[0] == 'x':
            vector = vector[indices ^ (1 << op[1])]
        elif op[0] == 'h':
            zero = indices[(indices & (1 << op[1])) == 0]
            one = zero | (1 << op[1])
            a, b = vector[zero].copy(), vector[one].copy()
            vector[zero], vector[one] = (a + b) / math.sqrt(2), (a - b) / math.sqrt(2)
        elif op[0] == 'cz':
            active = ((indices >> op[1]) & 1) & ((indices >> op[2]) & 1)
            vector *= 1 - 2 * active
        elif op[0] == 'z':
            vector *= 1 - 2 * ((indices >> op[1]) & 1)
        elif op[0] == 'ccx':
            active = ((indices >> op[1]) & 1) & ((indices >> op[2]) & 1)
            vector = vector[indices ^ (active << op[3])]
        else:
            raise GFMeasureError('Unknown vector operation')

    for op in operations:
        if op[0] == 'measure':
            if cursor >= len(outcomes) or outcomes[cursor] not in (0, 1):
                raise GFMeasureError('Every binary measurement outcome is required')
            bit = outcomes[cursor]
            cursor += 1
            classical[op[2]] = bit
            vector[((indices >> op[1]) & 1) != bit] = 0
        elif op[0] == 'if':
            if _eval_condition(op[1], classical):
                apply_gate(op[2])
        else:
            apply_gate(op)
    if cursor != len(outcomes):
        raise GFMeasureError('Surplus measurement outcomes')
    return vector


def multiply_poly(x, y, poly, n):
    """Independent textbook reference; not used when translating a circuit."""
    out = 0
    while y:
        if y & 1:
            out ^= x
        y >>= 1
        x <<= 1
        if x & (1 << n):
            x ^= poly
    return out
