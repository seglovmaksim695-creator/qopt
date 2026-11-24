"""Independent check: run the genuine program on ProjectQ's own ``ClassicalSimulator``.

Nothing of the recorder, the translation or the provers is used.  The program runs on
``MainEngine(ClassicalSimulator())`` (no decomposition: math gates such as AddQuantum are executed by their own
``math_fun``), the input registers are written right after their allocation, and the output registers are read after
the final flush.  This tests the whole chain (program -> decomposition -> translation) against ProjectQ's own meaning of
the program on the sampled inputs.  It cannot see the deallocated wires (the simulator drops them).
"""
from __future__ import annotations

import random

from ._record import Unsupported, need_projectq

__all__ = ['run_native', 'native_check']


def _injector(values):
    from projectq.cengines import BasicEngine
    from projectq.ops import AllocateQubitGate
    from projectq.types import WeakQubitRef

    class Injector(BasicEngine):
        def __init__(self):
            super().__init__()
            self.seen = 0

        def is_available(self, cmd):
            return self.next_engine.is_available(cmd)

        def receive(self, command_list):
            for cmd in command_list:
                self.send([cmd])
                if isinstance(cmd.gate, AllocateQubitGate):
                    q = cmd.qubits[0][0]
                    if q.id in values:
                        self.main_engine.backend.write_bit(WeakQubitRef(self.main_engine, q.id), values[q.id])
                        self.seen += 1
    return Injector()


def run_native(program, values=None, ids=None, decompose=None):
    """Run ``program`` on the classical simulator.  ``values``: input register name -> int (written after allocation);
    ``ids``: {name: [qubit ids]} of the input registers (from a first run).  -> (output values, input ids, output ids)."""
    need_projectq()
    from projectq import MainEngine
    from projectq.backends import ClassicalSimulator
    inj = {}
    if values is not None:
        for name, v in values.items():
            for i, q in enumerate(ids[name]):
                inj[q] = (v >> i) & 1
    chain = [_injector(inj)] if inj else []
    if decompose is not None:             # ProjectQ's own decomposition into X/CNOT/Toffoli, then the classical simulator
        import projectq.libs.math
        import projectq.setups.decompositions
        from projectq.cengines import AutoReplacer, DecompositionRuleSet, InstructionFilter
        from ._record import _make_filter
        rules = DecompositionRuleSet(modules=[projectq.libs.math, projectq.setups.decompositions])
        chain += [AutoReplacer(rules), InstructionFilter(_make_filter('nct', decompose))]
    eng = MainEngine(ClassicalSimulator(), engine_list=chain)
    layout = program(eng)
    eng.flush()
    outs = {name: eng.backend.read_register(reg) for name, reg in layout['outputs'].items()}
    in_ids = {name: [q.id for q in reg] for name, reg in layout['inputs'].items()}
    return outs, in_ids


def confirm_decomposition(program, values, max_controls=2):
    """Run ``program`` on the classical simulator AFTER ProjectQ's decomposition on one input.  -> output values, or None."""
    try:
        _, ids = run_native(program)
        got, _ = run_native(program, values, ids, decompose=max_controls)
        return got
    except Exception:  # noqa: BLE001
        return None


def _vectors(widths, count, seed):
    rng = random.Random(seed)
    names = sorted(widths)
    out = [{n: 0 for n in names}, {n: (1 << widths[n]) - 1 for n in names}]
    for n in names:
        out.append({m: (1 if m == n else 0) for m in names})
    while len(out) < count:
        out.append({n: rng.getrandbits(widths[n]) for n in names})
    return out[:count]


def native_check(program, spec, widths, count=24, seed=7, exhaustive_bits=10):
    """Compare the native run with ``spec.fn`` on ``count`` inputs (all inputs when there are <= ``exhaustive_bits`` input
    bits).  -> dict(ok, tested, detail)."""
    dec = None                            # None: the simulator runs the gates itself; 2: after ProjectQ's decomposition
    try:
        try:
            _, in_ids = run_native(program)
        except ValueError as e:
            if 'Only support' not in str(e):
                raise
            dec = 2                       # e.g. Swap: the simulator has no rule, ProjectQ's own decomposition then
            _, in_ids = run_native(program, decompose=dec)
    except Unsupported:
        raise
    except Exception as e:  # noqa: BLE001
        return dict(ok=False, tested=0, detail=f'native run failed: {type(e).__name__}: {e}')
    names = sorted(widths)
    total = sum(widths.values())
    if total <= exhaustive_bits:
        vecs = []
        for lane in range(1 << total):
            o, v = 0, {}
            for n in names:
                v[n] = (lane >> o) & ((1 << widths[n]) - 1)
                o += widths[n]
            vecs.append(v)
    else:
        vecs = _vectors(widths, count, seed)
    tested = 0
    for v in vecs:
        if spec.promise is not None and not spec.promise(v):
            continue
        try:
            got, ids2 = run_native(program, v, in_ids, decompose=dec)
        except Exception as e:  # noqa: BLE001
            return dict(ok=False, tested=tested, detail=f'native run failed on {v}: {type(e).__name__}: {e}')
        if ids2 != in_ids:
            return dict(ok=False, tested=tested, detail='the program allocates different qubit ids on different runs')
        want = spec.fn(v)
        for name, val in want.items():
            if got.get(name) != val:
                return dict(ok=False, tested=tested, detail=f'on input {v} register {name} is {got.get(name)}, specification says {val}',
                            witness=v)
        tested += 1
    return dict(ok=True, tested=tested, via='decomposition' if dec else 'gates', detail=f'{tested} inputs agree with ProjectQ ClassicalSimulator'
                                               + (' (all inputs)' if total <= exhaustive_bits else ' (sampled)'))
