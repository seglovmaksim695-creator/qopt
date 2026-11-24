"""Record a ProjectQ program command by command.

The program is run on a ProjectQ ``MainEngine`` whose last engine is a plain recorder.  In front of the recorder the
engine chain is [AutoReplacer, InstructionFilter]: ProjectQ's own decomposition rules (``projectq.setups.decompositions``
and ``projectq.libs.math``) bring every command into the accepted gate set, so the translation never interprets a
high-level gate itself.  Two gate sets exist:

* ``mode='nct'``        X with at most ``max_controls`` controls (default 2: X, CNOT, Toffoli); allocate / deallocate;
* ``mode='clifford_t'`` X with at most one control, H, S, S^dagger, T, T^dagger, Z, CZ.

A single-qubit gate is accepted only if its 2x2 matrix is EXACTLY (no global phase) the matrix of one of the atoms.
Everything else (measurement, dirty ancillas, rotations, unknown gates, a decomposition that does not exist) is an
``Unsupported`` error, which the callers turn into the verdict UNKNOWN: nothing is guessed.

A program is a function ``program(eng)`` that returns a dictionary::

    {'inputs':  {'x': qureg, 'y': qureg},      # registers whose initial value is free (all 2^n basis states)
     'outputs': {'x': qureg, 'result': qureg}, # registers that carry the specified result
     'garbage': {'g': qureg}}                   # optional: registers whose final value is declared unspecified

Wires of the recording are numbered in the order of their allocation.
"""
from __future__ import annotations

import hashlib
import json

__all__ = ['Unsupported', 'ProjectQUnavailable', 'Recording', 'record', 'need_projectq', 'atom_of']


class Unsupported(ValueError):
    """A command the strict translation does not understand (verdict UNKNOWN, never PROVED)."""


class ProjectQUnavailable(ImportError):
    pass


def need_projectq():
    try:
        import projectq  # noqa: F401
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise ProjectQUnavailable("ProjectQ is not installed: pip install 'qoptct[projectq]'") from e


_ATOMS = None


def _atoms():
    global _ATOMS
    if _ATOMS is None:
        import numpy as np
        r = 2 ** -0.5
        e8 = complex(r, r)
        _ATOMS = {
            'x': np.array([[0, 1], [1, 0]], dtype=complex),
            'z': np.array([[1, 0], [0, -1]], dtype=complex),
            'h': np.array([[r, r], [r, -r]], dtype=complex),
            's': np.array([[1, 0], [0, 1j]], dtype=complex),
            'sdg': np.array([[1, 0], [0, -1j]], dtype=complex),
            't': np.array([[1, 0], [0, e8]], dtype=complex),
            'tdg': np.array([[1, 0], [0, e8.conjugate()]], dtype=complex),
        }
    return _ATOMS


def atom_of(gate):
    """Name of the atom whose matrix is exactly the matrix of ``gate`` (no global phase allowed), or None."""
    import numpy as np
    try:
        m = np.array(gate.matrix, dtype=complex)
    except Exception:  # noqa: BLE001 - gates without a matrix
        return None
    if m.shape != (2, 2):
        return None
    for name, a in _atoms().items():
        if np.allclose(m, a, rtol=0, atol=1e-12):
            return name
    return None


class Recording:
    """Neutral list of operations of one recorded program.

    ``ops``: ``('alloc', w)``, ``('free', w)``, ``('x', controls, target)`` with ``controls`` a tuple of wires,
    ``('h'|'s'|'sdg'|'t'|'tdg'|'z', w)``, ``('cz', a, b)``.  ``inputs`` / ``outputs`` / ``garbage``: register name ->
    list of wires (LSB first, as the qureg was given).  ``ids``: ProjectQ qubit id of every wire.
    """

    def __init__(self, ops, n_wires, inputs, outputs, garbage, ids, mode, meta=None):
        self.ops = ops
        self.n_wires = n_wires
        self.inputs = inputs
        self.outputs = outputs
        self.garbage = garbage
        self.ids = ids
        self.mode = mode
        self.meta = meta or {}

    def copy(self, ops=None):
        return Recording(list(self.ops if ops is None else ops), self.n_wires, {k: list(v) for k, v in self.inputs.items()},
                         {k: list(v) for k, v in self.outputs.items()}, {k: list(v) for k, v in self.garbage.items()},
                         list(self.ids), self.mode, dict(self.meta))

    def sha256(self):
        blob = json.dumps([self.mode, self.n_wires, self.ops, sorted(self.inputs.items()), sorted(self.outputs.items()),
                           sorted(self.garbage.items())], separators=(',', ':'), default=list).encode()
        return hashlib.sha256(blob).hexdigest()

    def free_input_wires(self):
        """Input wires in a fixed order: registers sorted by name, bits LSB first."""
        return [w for name in sorted(self.inputs) for w in self.inputs[name]]

    def counts(self):
        c = {}
        for op in self.ops:
            k = op[0] if op[0] != 'x' else 'x%d' % len(op[1])
            c[k] = c.get(k, 0) + 1
        return c


def _make_filter(mode, max_controls):
    from projectq.ops import AllocateQubitGate, DeallocateQubitGate, FlushGate

    def accept(eng, cmd):
        g = cmd.gate
        if isinstance(g, (AllocateQubitGate, DeallocateQubitGate, FlushGate)):
            return True
        if len(cmd.qubits) != 1 or len(cmd.qubits[0]) != 1:
            return False
        a = atom_of(g)
        nc = len(cmd.control_qubits)
        if mode == 'nct':
            return a == 'x' and (max_controls is None or nc <= max_controls)
        if a in ('x', 'z'):
            return nc <= 1
        return a in ('h', 's', 'sdg', 't', 'tdg') and nc == 0
    return accept


def _ops_of(cmds):
    """Commands -> (neutral ops, number of wires, id of every wire)."""
    from projectq.ops import AllocateQubitGate, DeallocateQubitGate, FlushGate
    wire = {}
    ids = []
    ops = []
    freed = set()
    for k, c in enumerate(cmds):
        g = c.gate
        qs = [q.id for reg in c.qubits for q in reg]
        if isinstance(g, FlushGate):
            continue
        if isinstance(g, AllocateQubitGate):
            if c.control_qubits or len(qs) != 1 or qs[0] in wire:
                raise Unsupported(f'command {k}: bad or repeated allocation')
            wire[qs[0]] = len(ids)
            ids.append(qs[0])
            ops.append(('alloc', wire[qs[0]]))
            continue
        for q in qs + [q.id for q in c.control_qubits]:
            if q not in wire:
                raise Unsupported(f'command {k}: qubit {q} is used before it is allocated')
            if wire[q] in freed:
                raise Unsupported(f'command {k}: qubit {q} is used after it was freed')
        if isinstance(g, DeallocateQubitGate):
            if c.control_qubits or len(qs) != 1:
                raise Unsupported(f'command {k}: bad deallocation')
            freed.add(wire[qs[0]])
            ops.append(('free', wire[qs[0]]))
            continue
        if len(qs) != 1:
            raise Unsupported(f'command {k}: {type(g).__name__} on {len(qs)} qubits is not an atom')
        a = atom_of(g)
        if a is None:
            raise Unsupported(f'command {k}: {type(g).__name__} is not one of the exact atoms X, Z, H, S, Sdg, T, Tdg')
        ctl = [q.id for q in c.control_qubits]
        state = str(c.control_state) if ctl else ''
        if len(state) != len(ctl) or any(ch not in '01' for ch in state):
            raise Unsupported(f'command {k}: unreadable control state {c.control_state!r}')
        tw = wire[qs[0]]
        cw = [wire[q] for q in ctl]
        by_wire = {}
        for w, ch in zip(cw, state):                    # the same control twice is the same condition (AND is idempotent);
            if by_wire.setdefault(w, ch) != ch:         # opposite control states on one wire never fire: not an atom
                raise Unsupported(f'command {k}: one wire is a control with both states')
        if len(by_wire) != len(cw):
            cw, state = list(by_wire), ''.join(by_wire.values())
        if tw in by_wire:
            raise Unsupported(f'command {k}: the target is also a control')
        flips = [w for w, ch in zip(cw, state) if ch == '0']
        for w in flips:
            ops.append(('x', (), w))
        if a == 'x':
            ops.append(('x', tuple(cw), tw))
        elif a == 'z' and len(cw) == 1:
            ops.append(('cz', cw[0], tw))
        elif a == 'z' and not cw:
            ops.append(('z', tw))
        elif not cw and a in ('h', 's', 'sdg', 't', 'tdg'):
            ops.append((a, tw))
        else:
            raise Unsupported(f'command {k}: controlled {a} is not supported')
        for w in reversed(flips):
            ops.append(('x', (), w))
    return ops, len(ids), ids


def _recorder():
    from projectq.cengines import BasicEngine

    class Rec(BasicEngine):
        def __init__(self):
            super().__init__()
            self.cmds = []

        def is_available(self, cmd):
            return True

        def receive(self, command_list):
            self.cmds.extend(command_list)
    return Rec()


def _registers(layout, key, required):
    regs = layout.get(key) if isinstance(layout, dict) else None
    if regs is None:
        if required:
            raise Unsupported(f"the program must return a dictionary with the key {key!r}")
        return {}
    if not isinstance(regs, dict) or any(not isinstance(k, str) for k in regs):
        raise Unsupported(f'layout[{key!r}] must map names to registers')
    return {k: [q.id for q in v] for k, v in regs.items()}


def record(program, mode='nct', max_controls=2, extra_modules=()):
    """Run ``program(eng)`` and return a :class:`Recording`.  Raises :class:`Unsupported` for anything that is not an atom."""
    need_projectq()
    import projectq
    import projectq.libs.math
    import projectq.setups.decompositions
    from projectq import MainEngine
    from projectq.cengines import AutoReplacer, DecompositionRuleSet, InstructionFilter
    if mode not in ('nct', 'clifford_t'):
        raise ValueError("mode must be 'nct' or 'clifford_t'")
    rules = DecompositionRuleSet(modules=[projectq.libs.math, projectq.setups.decompositions] + list(extra_modules))
    rec = _recorder()
    eng = MainEngine(rec, engine_list=[AutoReplacer(rules), InstructionFilter(_make_filter(mode, max_controls))])
    try:
        layout = program(eng)
        eng.flush()
    except Unsupported:
        raise
    except Exception as e:  # noqa: BLE001 - the engine raises NoGateDecompositionError and friends
        raise Unsupported(f'ProjectQ could not bring the program into the {mode} gate set: {type(e).__name__}: {e}') from e
    cmds = list(rec.cmds)
    ids_in = _registers(layout, 'inputs', True)
    ids_out = _registers(layout, 'outputs', True)
    ids_gar = _registers(layout, 'garbage', False)
    ops, n, ids = _ops_of(cmds)
    idw = {q: w for w, q in enumerate(ids)}
    try:
        inputs = {k: [idw[q] for q in v] for k, v in ids_in.items()}
        outputs = {k: [idw[q] for q in v] for k, v in ids_out.items()}
        garbage = {k: [idw[q] for q in v] for k, v in ids_gar.items()}
    except KeyError as e:
        raise Unsupported(f'a register of the layout has the qubit {e.args[0]}, which was never allocated') from e
    for name, regs in (('inputs', inputs), ('outputs', outputs), ('garbage', garbage)):
        flat = [w for v in regs.values() for w in v]
        if len(flat) != len(set(flat)):
            raise Unsupported(f'registers of {name} overlap')
    iw = {w for v in inputs.values() for w in v}
    if any(w in iw for v in garbage.values() for w in v):
        raise Unsupported('a garbage register overlaps an input register')
    layout = None  # release the registers; the recording is already taken
    try:
        eng.flush(deallocate_qubits=True)
    except Exception:  # noqa: BLE001
        pass
    return Recording(ops, n, inputs, outputs, garbage, ids, mode, dict(projectq=_pq_version()))


def _pq_version():
    try:
        from importlib.metadata import version
        return version('projectq')
    except Exception:  # noqa: BLE001
        return None
