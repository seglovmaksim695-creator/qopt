"""Independent execution of the logical operations and of the phases of every measurement branch.

The least significant bit of an integer is wire 0. The phase of branch m is
(-1)**sum(m[j]*phase[j]); the common amplitude modulus is 2**(-len(phase)/2).
The code imports neither qoptct, Qualtran, Cirq nor Qiskit.
"""
from .sim_parse import SimulationParseError, parse_qasm


class SimulationError(ValueError):
    pass


def run(nq, ops, initial=0):
    """Return (final bits, phase coefficients per measurement outcome).

    With a nonzero erase coefficient the circuit is still executed: the coefficient
    is needed to detect a spoiled phase correction. A non-clean AND target
    is rejected because the classical model of this gadget does not apply then.
    """
    if not isinstance(nq, int) or isinstance(nq, bool) or nq < 0:
        raise SimulationError('The number of wires must be a non-negative integer')
    if not isinstance(initial, int) or isinstance(initial, bool) or not 0 <= initial < (1 << nq):
        raise SimulationError('Initial state outside the range of the wires')
    value = initial
    phases = []

    def wires(*qs):
        if any(not isinstance(q, int) or isinstance(q, bool) or not 0 <= q < nq for q in qs):
            raise SimulationError('Wire index out of range: %r' % (qs,))
        if len(set(qs)) != len(qs):
            raise SimulationError('Repeated wire in an operation: %r' % (qs,))

    def bit(q):
        return (value >> q) & 1

    for index, op in enumerate(ops):
        try:
            if not isinstance(op, (tuple, list)) or not op:
                raise SimulationError('Invalid operation: %r' % (op,))
            if op[0] == 'x' and len(op) == 2:
                t = op[1]
                wires(t)
                value ^= 1 << t
            elif op[0] == 'cx' and len(op) == 3:
                a, t = op[1:]
                wires(a, t)
                value ^= bit(a) << t
            elif op[0] == 'and' and len(op) == 4:
                a, b, t = op[1:]
                wires(a, b, t)
                if bit(t):
                    raise SimulationError('The AND target q[%d] is not zero' % t)
                value ^= (bit(a) & bit(b)) << t
            elif op[0] == 'erase' and len(op) == 3:
                t, corrections = op[1:]
                wires(t)
                coefficient = bit(t)
                for correction in corrections:
                    if correction[0] == 'z' and len(correction) == 2:
                        a = correction[1]
                        wires(t, a)
                        coefficient ^= bit(a)
                    elif correction[0] == 'cz' and len(correction) == 3:
                        a, b = correction[1:]
                        wires(t, a, b)
                        coefficient ^= bit(a) & bit(b)
                    else:
                        raise SimulationError('Unknown phase correction: %r' % (correction,))
                phases.append(coefficient)
                value &= ~(1 << t)
            else:
                raise SimulationError('Unknown operation: %r' % (op,))
        except (IndexError, TypeError) as exc:
            raise SimulationError('Operation %d: wrong format %r' % (index, op)) from exc
        except SimulationError as exc:
            raise SimulationError('Operation %d: %s' % (index, exc)) from exc
    return value, phases


def run_batch(nq, ops, inputs):
    """A simple sequential batch, without threads or child processes."""
    return [run(nq, ops, initial) for initial in inputs]
