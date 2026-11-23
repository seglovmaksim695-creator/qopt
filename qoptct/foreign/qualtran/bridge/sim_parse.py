"""Independent strict parse of the actual bridge text for simulation.

There are no imports of emit, rules or qoptct. The AND template is written independently.
Only the canonical dialect is supported, not general OpenQASM.
"""
import re


class SimulationParseError(ValueError):
    pass


_NUMBER = r'(?:0|[1-9][0-9]*)'
_WIRE = r'q\[' + _NUMBER + r'\]'
_GATE = re.compile(r'([a-z]+) (' + _WIRE + r'(?:, ' + _WIRE + r')*);')
_INDEX = re.compile(r'q\[(' + _NUMBER + r')\]')
_QREG = re.compile(r'qubit\[(' + _NUMBER + r')\] q;')
_CREG = re.compile(r'bit\[(' + _NUMBER + r')\] c;')
_MEASURE = re.compile(r'c\[(' + _NUMBER + r')\] = measure q\[(' + _NUMBER + r')\];')
_IF = re.compile(r'if \(c\[(' + _NUMBER + r')\]\) \{(.*)\}')
_ARITY = {'x': 1, 'cx': 2, 'h': 1, 't': 1, 'tdg': 1, 'sdg': 1, 'z': 1, 'cz': 2}


def parse_qasm(text):
    """Return (nq, ops), checking the whole text and the indices without gaps.

    An erase without CZ/Z is syntactically allowed: its wrong phase is found by run.
    Declared but unused classical bits are allowed. Real measurements
    must go to c[0], c[1], ...; reuse is forbidden.
    """
    if not isinstance(text, str):
        raise SimulationParseError('QASM text expected, not a path or an object')
    lines = text.splitlines()
    if len(lines) < 4 or lines[:2] != ['OPENQASM 3;', 'include "stdgates.inc";']:
        raise SimulationParseError('The standard four-line bridge header was expected')
    qmatch = _QREG.fullmatch(lines[2])
    cmatch = _CREG.fullmatch(lines[3])
    if not qmatch or not cmatch:
        raise SimulationParseError('Non-canonical qubit/bit declaration')
    nq, nb = int(qmatch.group(1)), int(cmatch.group(1))
    if nq < 1 or nb < 1:
        raise SimulationParseError('The qubit and bit registers must have a positive size')

    def fail(index, message):
        raise SimulationParseError('Line %d: %s' % (index + 1, message))

    def check_wires(index, qs):
        if any(q < 0 or q >= nq for q in qs):
            fail(index, 'wire out of range 0..%d: %r' % (nq - 1, qs))
        if len(set(qs)) != len(qs):
            fail(index, 'repeated wire: %r' % (qs,))

    def gate(index, line=None):
        if line is None:
            if index >= len(lines):
                fail(index, 'truncated gadget')
            line = lines[index]
        match = _GATE.fullmatch(line)
        if not match:
            fail(index, 'unsupported or non-canonical gate: %r' % line)
        name = match.group(1)
        qs = tuple(int(q) for q in _INDEX.findall(match.group(2)))
        if name not in _ARITY or len(qs) != _ARITY[name]:
            fail(index, 'unsupported name or number of operands: %r' % line)
        check_wires(index, qs)
        return name, qs

    ops = []
    i = 4
    next_measurement = 0
    while i < len(lines):
        name, qs = gate(i)
        if name == 'x':
            ops.append(('x', qs[0]))
            i += 1
            continue
        if name == 'cx':
            ops.append(('cx', *qs))
            i += 1
            continue
        if name != 'h':
            fail(i, 'gate outside an allowed gadget')
        target = qs[0]
        if i + 1 >= len(lines):
            fail(i, 'truncated gadget after H')
        if lines[i + 1] == 't q[%d];' % target:
            if i + 10 > len(lines):
                fail(i, 'truncated ten-line AND gadget')
            name_b, operands_b = gate(i + 2)
            name_a, operands_a = gate(i + 4)
            if name_b != 'cx' or name_a != 'cx':
                fail(i, 'the AND gadget has no controlling CX')
            b, target_b = operands_b
            a, target_a = operands_a
            if target_b != target or target_a != target:
                fail(i, 'the CX targets do not match the target of the AND gadget')
            check_wires(i, (a, b, target))
            expected = [
                'h q[%d];' % target,
                't q[%d];' % target,
                'cx q[%d], q[%d];' % (b, target),
                'tdg q[%d];' % target,
                'cx q[%d], q[%d];' % (a, target),
                't q[%d];' % target,
                'cx q[%d], q[%d];' % (b, target),
                'tdg q[%d];' % target,
                'h q[%d];' % target,
                'sdg q[%d];' % target,
            ]
            for j, expected_line in enumerate(expected):
                if lines[i + j] != expected_line:
                    fail(i + j, 'the text does not match the AND gadget')
            ops.append(('and', a, b, target))
            i += 10
            continue

        measurement = _MEASURE.fullmatch(lines[i + 1])
        if measurement is None:
            fail(i + 1, 'after H an AND or a measurement erase was expected')
        k, measured_target = map(int, measurement.groups())
        check_wires(i + 1, (measured_target,))
        if measured_target != target:
            fail(i + 1, 'another wire is measured')
        if k != next_measurement or not 0 <= k < nb:
            fail(i + 1, 'measurements must use the declared c[0], c[1], ... in order')
        if i + 2 >= len(lines):
            fail(i + 2, 'no conditional block after the measurement')
        correction = _IF.fullmatch(lines[i + 2])
        if correction is None or int(correction.group(1)) != k:
            fail(i + 2, 'unsupported erase condition')
        body = correction.group(2)
        pieces = body.split(';')
        if pieces[-1] != '' or not pieces[:-1]:
            fail(i + 2, 'empty or truncated conditional block')
        statements = [piece.strip() + ';' for piece in pieces[:-1]]
        if ' '.join(statements) != body:
            fail(i + 2, 'non-canonical text of the conditional block')
        if statements[-1] != 'x q[%d];' % target:
            fail(i + 2, 'the erase does not end with X on the measured target')
        corrections = []
        for statement in statements[:-1]:
            correction_name, correction_wires = gate(i + 2, statement)
            if correction_name not in ('z', 'cz') or target in correction_wires:
                fail(i + 2, 'unsupported erase phase correction')
            corrections.append((correction_name, *correction_wires))
        ops.append(('erase', target, tuple(corrections)))
        next_measurement += 1
        i += 3
    return nq, ops
