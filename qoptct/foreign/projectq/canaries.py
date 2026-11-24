"""Canaries (mutants): deliberately wrong versions of a program.  None of them may be PROVED.

Mutants are made on the recording (wires keep their numbers) or on the logical program / its OpenQASM text:

group ``junction``  an operation at the first or the last use of a wire (the seam between two blocks of the circuit);
group ``erasure``   the measurement-based erasure of a logical AND: wrong correction pair, missing correction, AND -> CNOT;
group ``phase``     a T/Tdg/S/Sdg of the AND gadget replaced (a wrong phase); not a classical circuit any more;
group ``rare``      a change that is visible on a small part of the inputs (a 3-control gate added to an output wire,
                    a Toffoli control replaced by another wire);
group ``other``     dropped / retargeted / reversed / duplicated gates, flipped output or garbage bits.
"""
from __future__ import annotations

import random

from ._logical import to_logical

__all__ = ['make_canaries', 'run_canaries']


def _uses(ops):
    first, last = {}, {}
    for i, op in enumerate(ops):
        if op[0] in ('alloc', 'free'):
            continue
        ws = list(op[1]) + [op[2]] if op[0] == 'x' else list(op[1:])
        for w in ws:
            first.setdefault(w, i)
            last[w] = i
    return first, last


def _rec_mutant(rec, ops, desc, group):
    m = rec.copy(ops)
    return dict(desc=desc, group=group, obj=m)


def make_canaries(rec, count=24, seed=1, pin=None):
    """-> list of dicts(desc, group, obj); ``obj`` is a Recording or a Logical (what ``check_projectq`` takes)."""
    rng = random.Random(seed)
    ops = rec.ops
    first, last = _uses(ops)
    cx = [i for i, o in enumerate(ops) if o[0] == 'x' and len(o[1]) == 1]
    cc = [i for i, o in enumerate(ops) if o[0] == 'x' and len(o[1]) == 2]
    seam_first = sorted({i for i in first.values() if i in set(cx)})
    seam_last = sorted({i for i in last.values() if i in set(cx)})
    free_in = rec.free_input_wires()
    res_wires = [w for name in sorted(rec.outputs) for w in rec.outputs[name]]
    gar = sorted(pin.wires) if pin is not None else []
    lg = to_logical(rec)
    erases = [i for i, o in enumerate(lg.ops) if o[0] == 'erase']
    ands = [i for i, o in enumerate(lg.ops) if o[0] == 'and']
    out = []

    def drop(i):
        return ops[:i] + ops[i + 1:]

    def sub(i, op):
        return ops[:i] + [op] + ops[i + 1:]

    def pick(lst, k):
        return lst[k % len(lst)] if lst else None

    makers = []
    # junction
    makers += [lambda k: ('junction', f'drop cx #{pick(seam_first, k * 7 + 1)} (first use of a wire)', drop(pick(seam_first, k * 7 + 1)))] if seam_first else []
    makers += [lambda k: ('junction', f'drop cx #{pick(seam_last, k * 5 + 2)} (last use of a wire)', drop(pick(seam_last, k * 5 + 2)))] if seam_last else []

    def retarget_seam(k):
        i = pick(seam_last, k * 3)
        c, t = ops[i][1][0], ops[i][2]
        t2 = (t + 1) % rec.n_wires
        t2 = t2 if t2 != c else (t2 + 1) % rec.n_wires
        return ('junction', f'retarget cx #{i} (last use of a wire)', sub(i, ('x', (c,), t2)))
    if seam_last:
        makers.append(retarget_seam)
    # rare
    if len(free_in) >= 3 and res_wires:
        def rare3(k):
            a, b, c = rng.sample(free_in, 3)
            t = res_wires[k % len(res_wires)]
            return ('rare', f'append a 3-control X ({a},{b},{c}) -> {t}', ops + [('x', (a, b, c), t)])
        makers.append(rare3)
    if cc:
        def ctrl_sub(k):
            i = pick(cc, k * 11 + 3)
            c0, c1 = ops[i][1]
            others = [x for x in range(rec.n_wires) if x not in (c0, c1, ops[i][2])]
            if not others:
                return ('rare', f'Toffoli #{i} -> CNOT (one control dropped)', sub(i, ('x', (c1,), ops[i][2])))
            w = rng.choice([x for x in free_in if x in others] or others)
            return ('rare', f'replace a control of the Toffoli #{i} by wire {w}', sub(i, ('x', (c0, w), ops[i][2])))
        makers.append(ctrl_sub)

        def cc_drop_ctl(k):
            i = pick(cc, k * 13 + 5)
            return ('rare', f'Toffoli #{i} -> CNOT (one control dropped)', sub(i, ('x', (ops[i][1][1],), ops[i][2])))
        makers.append(cc_drop_ctl)
    # other
    if cx:
        makers.append(lambda k: ('other', f'drop cx #{cx[rng.randrange(len(cx))]}', None))
        makers.append(lambda k: ('other', 'reverse a cx', None))
        makers.append(lambda k: ('other', 'duplicate a cx', None))
    if cc:
        makers.append(lambda k: ('other', 'drop a Toffoli', None))
    if res_wires:
        makers.append(lambda k: ('other', f'flip result wire {res_wires[k % len(res_wires)]} at the end', ops + [('x', (), res_wires[k % len(res_wires)])]))
    if free_in and res_wires:
        makers.append(lambda k: ('other', 'append cx input -> result', ops + [('x', (rng.choice(free_in),), res_wires[k % len(res_wires)])]))
    if gar:
        makers.append(lambda k: ('other', f'flip garbage wire {gar[k % len(gar)]} at the end', ops + [('x', (), gar[k % len(gar)])]))
        makers.append(lambda k: ('other', f'cx from an input into garbage wire {gar[(k * 3) % len(gar)]}', ops + [('x', (rng.choice(free_in),), gar[(k * 3) % len(gar)])]))
    k = 0
    tries = 0
    while len(out) < count and tries < count * 20:
        tries += 1
        mk = makers[k % len(makers)] if makers else None
        k += 1
        if mk is None:
            break
        group, desc, new = mk(k)
        if new is None:
            if desc.startswith('drop cx'):
                i = int(desc.split('#')[1])
                new = drop(i)
            elif desc == 'reverse a cx':
                i = cx[rng.randrange(len(cx))]
                c, t = ops[i][1][0], ops[i][2]
                desc, new = f'reverse cx #{i}', sub(i, ('x', (t,), c))
            elif desc == 'duplicate a cx':
                i = cx[rng.randrange(len(cx))]
                desc, new = f'duplicate cx #{i}', ops[:i + 1] + ops[i:]
            elif desc == 'drop a Toffoli':
                i = cc[rng.randrange(len(cc))]
                desc, new = f'drop Toffoli #{i}', drop(i)
        out.append(_rec_mutant(rec, new, desc, group))
    # logical / text level
    extra = []
    for j in range(max(4, count // 4)):
        if erases:
            i = erases[(j * 17 + 3) % len(erases)]
            o = lg.ops[i]
            a, b = o[2][0][1:]
            c = free_in[j % len(free_in)]
            c = c if c not in (a, b, o[1]) else (c + 1)
            m = lg.copy(lg.ops[:i] + [('erase', o[1], (('cz', a, c),))] + lg.ops[i + 1:])
            extra.append(dict(desc=f'erasure #{i}: wrong correction pair', group='erasure', obj=m))
        if erases and j % 2 == 0:
            i = erases[(j * 29 + 1) % len(erases)]
            o = lg.ops[i]
            m = lg.copy(lg.ops[:i] + [('erase', o[1], ())] + lg.ops[i + 1:])
            extra.append(dict(desc=f'erasure #{i}: correction removed', group='erasure', obj=m))
        if ands:
            i = ands[(j * 23 + 2) % len(ands)]
            o = lg.ops[i]
            m = lg.copy(lg.ops[:i] + [('cx', o[1], o[3])] + lg.ops[i + 1:])
            extra.append(dict(desc=f'AND #{i} replaced by a CNOT', group='erasure', obj=m))
            m = lg.copy()
            m.qasm_edits = [('sdg q[', 's q[', (j * 31 + 7) % max(1, len(ands)))]
            extra.append(dict(desc=f'phase: Sdg of AND gadget {(j * 31 + 7) % len(ands)} replaced by S', group='phase', obj=m))
            m = lg.copy()
            m.qasm_edits = [('tdg q[', 't q[', (j * 37 + 11) % max(1, 2 * len(ands)))]
            extra.append(dict(desc=f'phase: a Tdg of the AND gadgets replaced by T (occurrence {(j * 37 + 11) % (2 * len(ands))})', group='phase', obj=m))
    merged = []
    for i in range(max(len(out), len(extra))):          # alternate, so that a short run still has every group
        if i < len(extra):
            merged.append(extra[i])
        if i < len(out):
            merged.append(out[i])
    return merged


def _same(a, b, mask):
    if mask is None:
        return a == b
    return all(((x ^ y) & mask) == 0 for x, y in zip(a, b))


def _final_state(lg, kind):
    from ._sim import masks_full, masks_w2, simulate
    free = lg.free_wires if kind == 'full' else list(lg.inputs['x']) + list(lg.inputs['y'])
    masks, full, _ = masks_full(len(free)) if kind == 'full' else masks_w2(len(free))
    return simulate(lg.parsed_ops(), lg.nq, dict(zip(free, masks)), full)


def run_canaries(rec, spec, pin=None, count=24, seed=1, check=None, extra=12):
    """Run the canaries through ``check`` (default: ``check_projectq`` with the independent native run switched off).

    A mutant that is PROVED is tested for EQUIVALENCE with the golden program (identical final value on every wire, on all
    inputs or on all weight <= 2 inputs, which decides degree <= 2 circuits): an equivalent mutant is not a mutant (e.g. a
    Toffoli control that is always 1) and is reported with status ``EQUIVALENT``; a non-equivalent PROVED mutant is a
    false PROVED.  Stops after ``count`` non-equivalent canaries.
    -> list of dict(desc, group, status, detail)."""
    from .check import check_projectq
    chk = check or (lambda obj: check_projectq(obj, spec, pin=pin, native=False))
    golden = to_logical(rec)
    gold_state = None
    dom_mask = None
    if pin is None or pin.kind == 'full':           # a promise limits the claim: equivalence is judged on the promised inputs
        try:
            from ._classical import expected_lanes
            if getattr(spec, 'promise', None) is not None:
                dom_mask = expected_lanes(golden, spec, len(golden.free_wires))[1]
        except Exception:  # noqa: BLE001
            dom_mask = None
    rows, decided = [], 0
    for c in make_canaries(rec, count + extra, seed, pin):
        if decided >= count:
            break
        try:
            r = chk(c['obj'])
            row = dict(desc=c['desc'], group=c['group'], status=r['status'], detail=str(r.get('detail'))[:160])
            if r['status'] == 'PROVED':
                kind = pin.kind if pin is not None else 'full'
                obj = c['obj']
                lg = obj if hasattr(obj, 'parsed_ops') else to_logical(obj)
                if gold_state is None:
                    gold_state = _final_state(golden, kind)
                if lg.nq == golden.nq and _same(gold_state, _final_state(lg, kind), dom_mask):
                    row['status'] = 'EQUIVALENT'
                    row['detail'] = 'identical final value on every wire: not a real mutant'
        except Exception as e:  # noqa: BLE001 - a crash is not a proof either
            row = dict(desc=c['desc'], group=c['group'], status='ERROR', detail=f'{type(e).__name__}: {e}'[:160])
        rows.append(row)
        if row['status'] != 'EQUIVALENT':
            decided += 1
    return rows
