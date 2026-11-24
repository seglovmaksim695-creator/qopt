"""Whole-domain check of a classical (x, cx, and, erase) program against a Python function of the input registers.

The program is simulated on ALL 2^n inputs at once (bit-parallel, see ``_sim``).  The claim that is proved (status PROVED):

* for every input in the domain (``promise`` true), every output register holds ``fn(inputs)[name]``;
* every other wire (not an output, not a pinned garbage wire) ends in |0> (this includes freed ancillas and the
  scratch wires of the AND chains);
* every pinned garbage wire ends in exactly the value of the pin, on every input of the whole space.

The AND gadget needs a clean target and its erasure needs ``target = p AND q`` on every lane; a violation is UNKNOWN.
"""
from __future__ import annotations

from ._logical import Logical
from ._record import Unsupported
from ._gf2 import GarbagePin
from ._sim import PreconditionError, masks_full, simulate

__all__ = ['ClassicalSpec', 'pin_classical', 'check_classical', 'expected_lanes']


class ClassicalSpec:
    """``fn(values) -> {output name: int}`` where ``values`` maps every input register name to its integer value (bit i =
    i-th qubit of the register).  ``promise(values) -> bool`` restricts the claim to a domain (default: everything)."""

    def __init__(self, fn, promise=None, name='classical'):
        self.fn, self.promise, self.name = fn, promise, name


def _offsets(lg):
    off, k = {}, 0
    for name in sorted(lg.inputs):
        off[name] = (k, len(lg.inputs[name]))
        k += len(lg.inputs[name])
    return off, k


def expected_lanes(lg, spec, n):
    """-> (expected lane integer per output wire, promise mask integer) over the 2^n lanes."""
    import numpy as np
    off, _ = _offsets(lg)
    L = 1 << n
    wires = {}
    for name, ws in lg.outputs.items():
        for i, w in enumerate(ws):
            if w in wires:
                raise ValueError('two output registers share a wire')
            wires[w] = np.zeros(L, dtype=np.uint8)
    prom = np.ones(L, dtype=np.uint8)
    for lane in range(L):
        vals = {name: (lane >> o) & ((1 << ln) - 1) for name, (o, ln) in off.items()}
        if spec.promise is not None and not spec.promise(vals):
            prom[lane] = 0
            continue
        res = spec.fn(vals)
        for name, ws in lg.outputs.items():
            v = res[name]
            if v < 0 or v >> len(ws):
                raise ValueError(f'the specification returned {v} for {name}, which does not fit in {len(ws)} bits')
            for i, w in enumerate(ws):
                wires[w][lane] = (v >> i) & 1

    def pack(a):
        return int.from_bytes(np.packbits(a, bitorder='little').tobytes(), 'little')
    return {w: pack(a) for w, a in wires.items()}, pack(prom)


def _init(lg, n):
    masks, full, _ = masks_full(n)
    return dict(zip(lg.free_wires, masks)), full


def pin_classical(lg, max_input_bits=20):
    """Pin of the golden program: every wire that is not an output and ends non-zero, over all 2^n inputs."""
    n = len(lg.free_wires)
    if n > max_input_bits:
        raise ValueError('too many input bits for a full-lane pin')
    init, full = _init(lg, n)
    st = simulate(lg.ops, lg.nq, init, full)
    keep = {w for ws in lg.outputs.values() for w in ws}
    return GarbagePin('full', n, {w: st[w] for w in range(lg.nq) if w not in keep and st[w]})


def check_classical(lg: Logical, spec: ClassicalSpec, pin=None, max_input_bits=20):
    """-> dict(status PROVED | REFUTED | UNKNOWN, detail, witness?).  Does not run the independent native check."""
    n = len(lg.free_wires)
    out = dict(method='exhaustive-lanes', input_bits=n, domain='all inputs' if spec.promise is None else 'inputs satisfying the promise')
    if lg.garbage:
        return dict(out, status='UNKNOWN', detail='declared garbage registers are not supported by the classical check (use a pin)')
    if not lg.is_classical():
        return dict(out, status='UNKNOWN', detail='the program has Clifford+T atoms; it is not a classical program')
    if n > max_input_bits:
        return dict(out, status='UNKNOWN', detail=f'{n} input bits exceed max_input_bits={max_input_bits}')
    if pin is not None and (pin.kind != 'full' or pin.nbits != n):
        return dict(out, status='UNKNOWN', detail='the pin does not belong to this program')
    from ...verify.exact2 import gadget_check
    if gadget_check() > 1e-9:
        return dict(out, status='UNKNOWN', detail='the AND gadget check failed')
    try:
        exp, prom = expected_lanes(lg, spec, n)
    except ValueError as e:
        return dict(out, status='UNKNOWN', detail=str(e))
    init, full = _init(lg, n)
    try:
        st = simulate(lg.parsed_ops(), lg.nq, init, full)
    except Unsupported as e:
        return dict(out, status='UNKNOWN', detail=str(e))
    except PreconditionError as e:
        return dict(out, status='UNKNOWN', detail=f'AND/erasure precondition does not hold on every input: {e}')
    off, _ = _offsets(lg)

    def witness(lane):
        return {name: (lane >> o) & ((1 << ln) - 1) for name, (o, ln) in off.items()}
    for w, e in exp.items():                      # the outputs
        d = (st[w] ^ e) & prom
        if d:
            lane = (d & -d).bit_length() - 1
            return dict(out, status='REFUTED', detail=f'output wire {w} differs on input {witness(lane)}', witness=witness(lane))
    pinned = pin.wires if pin is not None else {}
    for w in range(lg.nq):
        if w in exp:
            continue
        if w in pinned:
            if st[w] != pinned[w]:
                d = st[w] ^ pinned[w]
                lane = (d & -d).bit_length() - 1
                return dict(out, status='REFUTED', detail=f'pinned garbage wire {w} differs from the pin on input {witness(lane)}',
                            witness=witness(lane))
        elif st[w] & prom:
            d = st[w] & prom
            lane = (d & -d).bit_length() - 1
            return dict(out, status='REFUTED', detail=f'wire {w} must end in 0 but does not, on input {witness(lane)}', witness=witness(lane))
    for w in pinned:
        if w in exp:
            return dict(out, status='UNKNOWN', detail='a pinned wire is an output wire')
    return dict(out, status='PROVED', detail=f'all {1 << n} inputs: outputs equal the specification, clean wires clean, '
                                             f'{len(pinned)} garbage wires equal the pin', garbage_wires=len(pinned),
                pin=None if pin is None else pin.digest)
