# SPDX-License-Identifier: Apache-2.0
"""Canaries: damaged copies of a real lowered program, and damaged proof data, that must NOT be proved.

The program is lowered once; every canary works on a copy (``clone_program``) and is proved by a fresh ``Engine``.
Two families:

* **joints** (the audit lessons of round 3): swapped wires of a child, wires of two children exchanged, an operation
  inserted between two children, an operation moved across a child boundary, a child with a wrong contract, a
  precondition of the root that is too weak (a child is then called outside its domain), a child contract with a
  narrower domain, a wrong loop invariant, a wrong loop step;
* **gates** (targeted single-operation mutants of the real gates): an erase correction removed or moved (phase), a
  CX control changed in an adder (rare inputs: the carry), an X dropped from a constant, two operations exchanged,
  the final operation removed.  A gate canary counts only if an independent simulation finds an input on which the
  damaged program differs from the specification (``witness``); otherwise it is reported as ``equivalent``.

``run_canaries(target, n)`` returns a dict with every canary, its status and witness; ``ok`` is True only when no
canary is PROVED and at least ``min_count`` canaries were made.
"""
from __future__ import annotations

import copy
import json
import random
import time
from dataclasses import replace
from pathlib import Path

from . import expr as X
from .compose import Engine
from .tree import Program

__all__ = ['run_canaries', 'canaries_child', 'clone_program']


def clone_program(P):
    ad = P.adapter
    ad2 = copy.copy(ad)
    ad2.events = dict(ad.events)
    ad2.spans = dict(ad.spans)
    ad2._nodes = {}
    ad2._bloq = ad._bloq
    ad2._decomp = ad._decomp
    ad2.registry = ad.registry
    P2 = Program(list(P.ops), P.nq, None, ad2, meta=dict(P.meta))
    P2.root = ad2.node(ad2.root_path)
    return P2


def _mut_event(P, path):
    e = P.adapter.events[path]
    e2 = dict(e, inputs_msb={k: list(v) for k, v in e['inputs_msb'].items()},
              outputs_msb={k: list(v) for k, v in e['outputs_msb'].items()})
    P.adapter.events[path] = e2
    P.adapter._nodes.pop(path, None)
    return e2


def _refresh(P):
    P.adapter._nodes = {}
    P._support = {}
    P.root = P.adapter.node(P.adapter.root_path)


def insert_op(P, pos, op, owner):
    """insert ``op`` at ``pos``, attributed to ``owner`` (and its ancestors) and to no other node"""
    P.ops.insert(pos, op)
    anc = {owner}
    p = owner
    while '/' in p:
        p = p.rsplit('/', 1)[0]
        anc.add(p)
    new = {}
    for path, ((a, b), ok) in P.adapter.spans.items():
        if path in anc:
            new[path] = ((min(a, pos), max(b, pos) + 1), ok and a <= pos <= b)
        elif a >= pos:
            new[path] = ((a + 1, b + 1), ok)
        elif b <= pos:
            new[path] = ((a, b), ok)
        else:                       # pos falls strictly inside a node that does not own the new operation
            new[path] = ((a, b + 1), False)
    for path in anc:
        if path not in new:
            new[path] = ((pos, pos + 1), True)
    P.adapter.spans = new
    _refresh(P)


def _proper_children(P, node):
    return [k for k in P.children(node) if k.kind == 'op']


def _leaves(P, node, out, depth=0, limit=40):
    """registry leaves (primitives) below node, in program order"""
    from . import registry as R
    for k in _proper_children(P, node):
        c = R.contract_for(k.obj)
        if c is not None and c.kind == 'primitive' and k.nops:
            out.append(k)
        elif k.nops and depth < 6:
            _leaves(P, k, out, depth + 1, limit)
        if len(out) >= limit:
            return out
    return out


# ------------------------------------------------------------------------------------------------ joint canaries
def joint_canaries(P0, spec, rng):
    """list of (name, family, function(P) -> description, spec override or None)"""
    out = []
    root = P0.root
    kids = _proper_children(P0, root)
    big = [k for k in kids if k.nops]
    multi = [k for k in big if any(len(ws) >= 2 for ws in k.ins.values())]

    def swap_in(path, reg_pick=0):
        def f(P):
            e = _mut_event(P, path)
            names = [k for k, v in e['inputs_msb'].items() if len(v) >= 2]
            nm = names[reg_pick % len(names)]
            v = e['inputs_msb'][nm]
            v[0], v[-1] = v[-1], v[0]
            if nm in e['outputs_msb'] and len(e['outputs_msb'][nm]) == len(v):
                w = e['outputs_msb'][nm]
                w[0], w[-1] = w[-1], w[0]
            _refresh(P)
            return f'{path}: register {nm}: its first and last wire exchanged (input and output)'
        return f
    # only children for which the exchange changes the function (by their contract): exchanging two wires of a
    # bitwise NOT, say, gives the same program, and such a canary would rightly be proved
    eff = [(k, j) for k in multi for j in range(2) if _swap_effective(k, j, rng)]
    if eff:
        out.append(('J1_swap_wires_of_a_child', 'joint', swap_in(eff[0][0].path, eff[0][1]), None))
        later = [e for e in eff[len(eff) // 2:] if e[0].path != eff[0][0].path]
        if later:
            out.append(('J1b_swap_wires_of_a_later_child', 'joint', swap_in(later[0][0].path, later[0][1]), None))

    def swap_out(path):
        def f(P):
            e = _mut_event(P, path)
            names = [k for k, v in e['outputs_msb'].items() if len(v) >= 2]
            nm = names[-1]
            w = e['outputs_msb'][nm]
            w[0], w[1] = w[1], w[0]
            _refresh(P)
            return f'{path}: output register {nm}: two wires exchanged (only on the output side)'
        return f
    if multi:
        out.append(('J2_swap_output_wires_of_a_child', 'joint', swap_out(multi[-1].path), None))

    pair = None
    for i in range(len(big)):
        for j in range(i + 1, len(big)):
            wa, wb = set(big[i].in_wires()), set(big[j].in_wires())
            if wa - wb and wb - wa:
                pair = (big[i], big[j])
                break
        if pair:
            break
    if pair is not None:
        a, b = pair

        def exch(P):
            ea, eb = _mut_event(P, a.path), _mut_event(P, b.path)
            wa = {q for v in ea['inputs_msb'].values() for q in v}
            wb = {q for v in eb['inputs_msb'].values() for q in v}
            ra = max(ea['inputs_msb'], key=lambda k: len([q for q in ea['inputs_msb'][k] if q not in wb]))
            rb = max(eb['inputs_msb'], key=lambda k: len([q for q in eb['inputs_msb'][k] if q not in wa]))
            qa = [q for q in ea['inputs_msb'][ra] if q not in wb][-1]
            qb = [q for q in eb['inputs_msb'][rb] if q not in wa][-1]
            for e, r, old, new in ((ea, ra, qa, qb), (eb, rb, qb, qa)):
                for side in ('inputs_msb', 'outputs_msb'):
                    if r in e[side]:
                        e[side][r] = [new if q == old else q for q in e[side][r]]
            _refresh(P)
            return f'{a.path}.{ra} and {b.path}.{rb}: wires {qa} and {qb} exchanged between the two children'
        out.append(('J3_wires_exchanged_between_children', 'joint', exch, None))

    if len(big) >= 2:
        k1 = big[0]

        def extra_between(P):
            pos = k1.span[1]
            q = P.root.in_wires()[0]
            insert_op(P, pos, ('x', q), P.adapter.root_path)
            return f'X on input wire {q} inserted at {pos}, between {k1.path} and the next child (owned by the root)'
        out.append(('J4_extra_operation_between_children', 'joint', extra_between, None))

        def extra_inside(P):
            pos = k1.span[0] + max(1, k1.nops // 2)
            q = k1.in_wires()[0]
            insert_op(P, pos, ('x', q), k1.path)
            return f'X on wire {q} inserted inside {k1.path} at {pos}'
        out.append(('J5_extra_operation_inside_a_child', 'joint', extra_inside, None))

    if big:
        first, last = big[0], big[-1]

        def extra_before(P):
            pos = first.span[0]
            q = P.root.in_wires()[-1]
            insert_op(P, pos, ('x', q), P.adapter.root_path)
            return f'X on input wire {q} inserted at {pos}, before the first child {first.path} (owned by the root)'
        out.append(('J4b_extra_operation_before_the_first_child', 'joint', extra_before, None))

        def extra_after(P):
            pos = last.span[1]
            q = P.root.out_wires()[0]
            insert_op(P, pos, ('x', q), P.adapter.root_path)
            return f'X on output wire {q} inserted at {pos}, after the last child {last.path} (owned by the root)'
        out.append(('J4c_extra_operation_after_the_last_child', 'joint', extra_after, None))

    # a composite child: an operation between two of its children; for the moved boundary (J7) a pair of
    # neighbouring grandchildren where the first one, without its last operation, no longer meets its contract
    from . import registry as R

    def composed(k):            # a child proved from its own children (a primitive is proved on its gates as a whole)
        c = R.contract_for(k.obj)
        return c is None or c.kind != 'primitive'
    comp = None
    for k in sorted(big, key=lambda k: not composed(k)):
        sub = [c for c in _proper_children(P0, k) if c.nops]
        if len(sub) >= 2:
            comp = (k, sub)
            break
    j7 = None
    for k in [k for k in big if composed(k)]:
        sub_k = [c for c in _proper_children(P0, k) if c.nops]
        for j in range(len(sub_k) - 1):
            s0, s1 = sub_k[j], sub_k[j + 1]
            if s0.span[1] == s1.span[0] and s0.nops >= 2 and \
                    _leaf_differs(P0, s0, lambda o: o.__delitem__(len(o) - 1), rng):
                j7 = (k, [s0, s1])
                break
        if j7:
            break
    if comp is not None:
        k, sub = comp

        def extra_between_sub(P, sub=sub, k=k):
            pos = sub[0].span[1]
            q = sub[0].out_wires()[0]
            insert_op(P, pos, ('x', q), k.path)
            return f'X on wire {q} inserted between {sub[0].path} and {sub[1].path} (owned by {k.path})'
        out.append(('J6_extra_operation_between_grandchildren', 'joint', extra_between_sub, None))

    if j7 is not None:
        k, sub = j7

        def move_boundary(P, sub=sub):
            (a0, b0), ok0 = P.adapter.spans[sub[0].path]
            (a1, b1), ok1 = P.adapter.spans[sub[1].path]
            if b0 != a1 or b0 - a0 < 2:
                return 'skipped'
            # the last operation of the first grandchild is attributed to the second one
            for path in list(P.adapter.spans):
                (a, b), ok = P.adapter.spans[path]
                if path == sub[0].path or path.startswith(sub[0].path + '/'):
                    if b == b0:
                        P.adapter.spans[path] = ((a, b - 1), ok)
                if path == sub[1].path:
                    P.adapter.spans[path] = ((a1 - 1, b1), ok)
            _refresh(P)
            return f'the last operation of {sub[0].path} attributed to {sub[1].path}'
        out.append(('J7_operation_moved_across_a_boundary', 'joint', move_boundary, None))

    # wrong contracts of children (registry entries replaced)
    leaves = _leaves(P0, root, [], limit=60)
    wrongs = []
    for lf in leaves:
        from . import registry as R
        c = R.contract_for(lf.obj)
        bad = _wrong_contract(c)
        if bad is not None and lf.type not in [w[0] for w in wrongs]:
            wrongs.append((lf.type, lf.key, bad))
        if len(wrongs) >= 3:
            break
    for i, (tname, key, bad) in enumerate(wrongs):
        out.append((f'J8_wrong_contract_of_child_{tname}', 'joint', _override_contract(key, bad), None))
    # a child that claims to leave a register unchanged (or a new register zero)
    idle = []
    for lf in leaves:
        from . import registry as R
        bad = _idle_contract(R.contract_for(lf.obj))
        if bad is not None and lf.type not in [w[0] for w in idle]:
            idle.append((lf.type, lf.key, bad))
        if len(idle) >= 2:
            break
    for tname, key, bad in idle:
        out.append((f'J8b_child_claims_no_effect_{tname}', 'joint', _override_contract(key, bad), None))

    # domain: the root precondition weakened (children called outside their domains)
    if spec.pre != 1:
        weak = _weaken(spec.pre)
        for j, wp in enumerate(weak[:2]):
            sp = replace(spec, pre=wp)
            out.append((f'J9_root_precondition_weakened_{j}', 'joint', lambda P: 'root precondition weakened', sp))
    # domain: a child contract with a narrower domain
    for lf in leaves:
        from . import registry as R
        c = R.contract_for(lf.obj)
        nar = _narrow(c)
        if nar is not None:
            out.append((f'J10_child_domain_narrowed_{lf.type}', 'joint', _override_contract(lf.key, nar), None))
            break
    # loops: a wrong invariant, a wrong step
    if spec.folds:
        f = spec.folds[0]
        if f.invariant != 1:
            inv2 = _weaken_inv(f.invariant)
            if inv2 is not None:
                f2 = replace(f, invariant=inv2)
                out.append(('J11_wrong_loop_invariant', 'joint', lambda P: 'loop invariant replaced',
                            replace(spec, folds=[f2] + list(spec.folds[1:]))))
        bc = f.body_contract
        bad = _wrong_body(bc)
        if bad is not None:
            f3 = replace(f, body_contract=bad)
            out.append(('J12_wrong_loop_step', 'joint', lambda P: 'loop step (body contract) replaced',
                        replace(spec, folds=[f3] + list(spec.folds[1:]))))
        # the chain of the loop broken: an iteration reads its carried register from other wires
        iters = [k for k in kids if k.type == f.body]
        if len(iters) >= 3 and f.state:
            it = iters[2]
            st = f.state[0]

            def break_chain(P, it=it, st=st):
                e = _mut_event(P, it.path)
                v = e['inputs_msb'][st]
                v[0], v[1] = v[1], v[0]
                if st in e['outputs_msb']:
                    w = e['outputs_msb'][st]
                    w[0], w[1] = w[1], w[0]
                _refresh(P)
                return f'{it.path}: carried register {st} read from exchanged wires'
            out.append(('J13_loop_chain_broken', 'joint', break_chain, None))
    return out


def _swap_effective(k, reg_pick, rng, cases=24):
    """does exchanging the most and the least significant wire of an input register of child k (input and output
    side, as swap_in does) change its function?  Judged on its registry contract; False without a contract"""
    from . import registry as R
    from .leaf import sample_inputs
    c = R.contract_for(k.obj)
    if c is None:
        return False
    names = [nm for nm, v in k.ins.items() if len(v) >= 2]
    if not names:
        return False
    nm = names[reg_pick % len(names)]
    w = c.inputs.get(nm)
    if not w or w < 2:
        return False

    def sw(v):
        hi, lo = (v >> (w - 1)) & 1, v & 1
        return (v & ~((1 << (w - 1)) | 1)) | (lo << (w - 1)) | hi
    try:
        samples = sample_inputs(c, random.Random(rng.randrange(1 << 30)), cases)
    except Exception:  # noqa: BLE001
        return False
    for v in samples:
        v2 = dict(v, **{nm: sw(v[nm])})
        if not c.in_domain(v2):
            return True
        a, b = c.evaluate(v), c.evaluate(v2)
        if nm in b and c.outputs.get(nm) == w:
            b = dict(b, **{nm: sw(b[nm])})
        if a != b:
            return True
    return False


def _leaf_differs(P0, lf, edit, rng, cases=200):
    """does the edit (``edit(ops)`` changes the list of the leaf's operations in place) change leaf lf, judged by
    simulation (every wire, phases) against its registry contract on sample inputs of its domain?"""
    from . import registry as R
    from .leaf import sample_inputs, simulate
    from ..qualtran.montgomery._core import compact
    c = R.contract_for(lf.obj)
    if c is None:
        return False
    a, b = lf.span
    ops = list(P0.ops[a:b])
    edit(ops)
    try:
        nq, cops, wi, wo, _ = compact(ops, dict(lf.ins), dict(lf.outs))
        cnt, bad = simulate(nq, cops, wi, wo, c, sample_inputs(c, random.Random(rng.randrange(1 << 30)), cases))
    except Exception:  # noqa: BLE001
        return True          # the damaged leaf does not even run: a change
    return bool(bad)


def _override_contract(key, bad):
    def f(P):
        from . import registry as R
        ad = P.adapter
        base = ad.registry if ad.registry is not None else R

        class Reg:
            @staticmethod
            def contract_for(b):
                c = base.contract_for(b)
                k = type(b).__name__ + ':' + __import__('hashlib').sha256(repr(b).encode()).hexdigest()[:16]
                return bad if k == key else c
        ad.registry = Reg
        _refresh(P)
        return f'contract of {key} replaced by a wrong one: {bad.name}'
    return f


def _first_out(c):
    for k, e in c.out.items():
        if e != ('var', k):
            return k
    return None


def _wrong_contract(c):
    if c is None:
        return None
    k = _first_out(c)
    if k is None:
        return None
    w = c.outputs[k]
    e = X.mod(X.add(c.out[k], 1), 1 << w)
    return replace(c, out=dict(c.out, **{k: e}), name=c.name + ' [WRONG: +1 on ' + k + ']')


def _idle_contract(c):
    """the contract with its first computed output replaced by 'unchanged' (or 0 for a new register)"""
    if c is None:
        return None
    k = _first_out(c)
    if k is None:
        return None
    e = X.var(k) if k in c.inputs and c.inputs[k] == c.outputs[k] else 0
    return replace(c, out=dict(c.out, **{k: e}), name=c.name + ' [WRONG: ' + k + ' unchanged]')


def _narrow(c):
    """the domain of the widest input register cut down to {0, 1}: the parent then calls the child outside it"""
    if c is None or not c.inputs:
        return None
    k = max(c.inputs, key=lambda r: c.inputs[r])
    w = c.inputs[k]
    if w < 3:
        return None
    pre = X.and_(c.pre, X.lt(X.var(k), 2)) if c.pre != 1 else X.lt(X.var(k), 2)
    return replace(c, pre=pre, name=c.name + f' [NARROW: {k} < 2]')


def _weaken(pre):
    """conjuncts dropped one at a time; for a single condition: always true"""
    if type(pre) is tuple and pre[0] == 'and':
        return [pre[2], pre[1]]
    return [1]


def _weaken_inv(inv):
    """replace the first comparison bound a < b + c by a < b (a wrong invariant)"""
    def go(e):
        if type(e) is tuple and e[0] == 'lt' and type(e[2]) is tuple and e[2][0] == 'add':
            return ('lt', e[1], e[2][1])
        if type(e) is tuple and e[0] == 'and':
            return ('and', go(e[1]), e[2])
        return None
    return go(inv)


def _wrong_body(bc):
    for nm, e in bc.defs:
        if type(e) is tuple and e[0] == 'div' and type(e[1]) is tuple and e[1][0] == 'add':
            inner = e[1]
            if type(inner[2]) is tuple and inner[2][0] == 'lut':
                tab = list(inner[2][2])
                if len(tab) > 1:
                    tab[1] += 2
                    e2 = ('div', ('add', inner[1], ('lut', inner[2][1], tuple(tab))), e[2])
                    defs = [(a, e2 if a == nm else b) for a, b in bc.defs]
                    return replace(bc, defs=defs, name=bc.name + ' [WRONG table entry]')
    return None


# ------------------------------------------------------------------------------------------------ gate canaries
def gate_canaries(P0, rng, count=14):
    """targeted single-operation mutants of the real program (positions chosen inside registry leaves)"""
    out = []
    leaves = _leaves(P0, P0.root, [], limit=80)
    by_type = {}
    for lf in leaves:
        by_type.setdefault(lf.type, []).append(lf)
    ops = P0.ops

    def positions(lf, kind):
        a, b = lf.span
        return [i for i in range(a, b) if ops[i][0] == kind]

    def differs(lf, edit):
        return _leaf_differs(P0, lf, edit, rng)

    def order(idx):
        """candidate positions, the preferred one first"""
        return idx[len(idx) // 2:] + idx[:len(idx) // 2][::-1]

    # phase: erase corrections (only where removing / moving the correction changes the leaf)
    g1 = g2 = None
    for lf in leaves:
        a0 = lf.span[0]
        for i in order(positions(lf, 'erase'))[:8]:
            op = ops[i]
            if not op[2]:
                continue
            if g1 is None and differs(lf, lambda o, i=i, op=op: o.__setitem__(i - a0, ('erase', op[1], ()))):
                g1 = (lf, i)
            others = [q for q in lf.in_wires() if q not in (op[2][0][1], op[2][0][2], op[1])]
            if g2 is None and others and len(op[2]) == 1 and differs(
                    lf, lambda o, i=i, op=op, q=others[0]: o.__setitem__(i - a0, ('erase', op[1],
                                                                                 (('cz', op[2][0][1], q),)))):
                g2 = (lf, i, others[0])
            if g1 and g2:
                break
        if g1 and g2:
            break
    if g1:
        lf, i = g1

        def drop_corr(P, i=i):
            op = P.ops[i]
            P.ops[i] = ('erase', op[1], ())
            P._support = {}
            return f'op {i}: erase correction removed (phase)'
        out.append((f'G1_phase_correction_removed_{lf.type}', 'gate', drop_corr, None))
    if g2:
        lf, i, q = g2

        def move_corr(P, i=i, q=q):
            op = P.ops[i]
            (cz, a, b), = op[2]
            P.ops[i] = ('erase', op[1], (('cz', a, q),))
            P._support = {}
            return f'op {i}: erase correction moved to wire {q} (phase)'
        out.append((f'G2_phase_correction_moved_{lf.type}', 'gate', move_corr, None))
    # rare inputs: a carry CX near the top of an adder, the AND of a high carry
    for tname in ('CAdd', 'Add', 'Controlled', 'AddK', 'LessThanConstant'):
        for lf in by_type.get(tname, [])[:1]:
            a0 = lf.span[0]
            cx = positions(lf, 'cx')
            ands = positions(lf, 'and')
            for i in [cx[-k] for k in range(3, min(len(cx), 12) + 1)] if len(cx) > 4 else []:
                op = ops[i]
                if op[1] != op[2] and differs(lf, lambda o, i=i, op=op: o.__setitem__(i - a0, ('cx', op[2], op[1]))):
                    def flip_cx(P, i=i):
                        op = P.ops[i]
                        P.ops[i] = ('cx', op[2], op[1])
                        P._support = {}
                        return f'op {i}: control and target of a CX near the top of the carry chain exchanged'
                    out.append((f'G3_carry_cx_reversed_{tname}', 'gate', flip_cx, None))
                    break
            for i in [ands[-k] for k in range(2, min(len(ands), 10) + 1)] if len(ands) > 3 else []:
                op = ops[i]
                if differs(lf, lambda o, i=i, op=op: o.__setitem__(i - a0, ('cx', op[1], op[3]))):
                    def drop_and(P, i=i):
                        op = P.ops[i]
                        P.ops[i] = ('cx', op[1], op[3])
                        P._support = {}
                        return (f'op {i}: the AND of a high carry replaced by a CX of one control '
                                f'(wrong only on rare inputs)')
                    out.append((f'G4_high_carry_and_weakened_{tname}', 'gate', drop_and, None))
                    break
    # an X dropped (constants are loaded with X)
    g5 = None
    for lf in leaves:
        a0 = lf.span[0]
        xs = positions(lf, 'x')
        if len(xs) >= 2:
            for i in order(xs)[:8]:
                if differs(lf, lambda o, i=i: o.__delitem__(i - a0)):
                    g5 = (lf, i)
                    break
        if g5:
            break
    if g5:
        lf, i = g5

        def drop_x(P, i=i):
            del P.ops[i]
            _shift_after_delete(P, i)
            return f'op {i}: an X removed'
        out.append((f'G5_x_removed_{lf.type}', 'gate', drop_x, None))
    # two neighbouring operations of a leaf exchanged (they must not commute)
    g6 = None
    for lf in leaves:
        a, b = lf.span
        for i in range(a, b - 1):
            o1, o2 = ops[i], ops[i + 1]
            if o1[0] == 'cx' and o2[0] == 'cx' and (o1[2] == o2[1] or o2[2] == o1[1]):
                def ed(o, i=i, a=a):
                    o[i - a], o[i + 1 - a] = o[i + 1 - a], o[i - a]
                if differs(lf, ed):
                    g6 = (lf, i)
                    break
        if g6:
            break
    if g6:
        lf, i = g6

        def swap2(P, i=i):
            P.ops[i], P.ops[i + 1] = P.ops[i + 1], P.ops[i]
            P._support = {}
            return f'ops {i}, {i + 1}: two non-commuting CX exchanged'
        out.append((f'G6_noncommuting_swap_{lf.type}', 'gate', swap2, None))
    # the last operation of the program removed (the final reduction)
    def drop_last(P):
        i = len(P.ops) - 1
        del P.ops[i]
        _shift_after_delete(P, i)
        return f'op {i} (the last one) removed'
    out.append(('G7_last_operation_removed', 'gate', drop_last, None))
    # random single-operation mutants inside random leaves; a mutant is kept only if the damaged leaf differs from
    # its contract on some input (simulation), so that an equivalent rewrite is not taken for a canary
    pool = [lf for lf in leaves if lf.nops >= 3]
    for j in range(max(0, count - len(out))):
        chosen = None
        for _try in range(30):
            if not pool:
                break
            lf = rng.choice(pool)
            i = rng.randrange(lf.span[0], lf.span[1])
            op = ops[i]
            wires = sorted(set(lf.in_wires()) | set(lf.out_wires()))
            new = op
            for _ in range(20):
                if op[0] == 'cx':
                    t = rng.choice(wires)
                    new = ('cx', op[1], t) if t != op[1] else op
                elif op[0] == 'x':
                    new = ('x', rng.choice(wires))
                elif op[0] == 'and':
                    a = rng.choice(wires)
                    new = ('and', a, op[2], op[3]) if a not in (op[2], op[3]) else op
                else:
                    new = ('erase', op[1], ())
                if new != op:
                    break
            if new != op and _leaf_differs(P0, lf, lambda o, i=i, a0=lf.span[0], new=new: o.__setitem__(i - a0, new),
                                           rng):
                chosen = (lf, i, op, new)
                break
        if chosen is None:
            continue
        lf, i, op, new = chosen

        def rnd(P, i=i, lf=lf, op=op, new=new):
            if P.ops[i] != op:
                return 'skipped'
            P.ops[i] = new
            P._support = {}
            return f'op {i} in {lf.path}: {op} -> {new}'
        out.append((f'G8_random_{j}_{lf.type}', 'gate', rnd, None))
    return out


def _shift_after_delete(P, i):
    new = {}
    for path, ((a, b), ok) in P.adapter.spans.items():
        if b <= i:
            new[path] = ((a, b), ok)
        elif a > i:
            new[path] = ((a - 1, b - 1), ok)
        else:
            new[path] = ((a, b - 1), ok) if b - 1 > a else None
    P.adapter.spans = {k: v for k, v in new.items() if v is not None}
    _refresh(P)


# ------------------------------------------------------------------------------------------------ runner
def _witness(P, spec, rng, E, cases=40):
    """an input on which the damaged program differs (independent simulation): end to end on the domain of the
    specification, else on the domain of the first leaf whose proof failed (the original leaf is proved to satisfy
    that contract on its whole domain, so a difference there is a real change of the program)"""
    from .api import end_to_end
    from .leaf import sample_inputs, simulate
    from ..qualtran.montgomery._core import compact
    r = end_to_end(P, spec, cases=cases, seed=rng.randrange(1 << 30))
    if r.get('mismatches'):
        r['level'] = 'end-to-end'
        return r
    for node, contract, lr in E.failed_leaves:
        if lr.get('witness') and lr.get('method') == 'independent_simulation':
            return dict(mismatches=1, cases=lr.get('sim_cases'), first=lr['witness'], level=f'leaf {node.path}')
        try:
            nq, cops, wi, wo, _ = compact(list(P.node_ops(node)), dict(node.ins), dict(node.outs))
            samples = sample_inputs(contract, random.Random(rng.randrange(1 << 30)), 400)
            cnt, bad = simulate(nq, cops, wi, wo, contract, samples)
        except Exception as e:  # noqa: BLE001
            continue
        if bad:
            return dict(mismatches=1, cases=cnt, first=bad, level=f'leaf {node.path}')
    r['level'] = 'none found'
    return r


def run_canaries(target, n, p=None, *, w=1, only=None, out=None, seed=5, log=None, min_count=20):
    from .api import build_target, _versions, _variant, source_digest
    from .qualtran_adapter import lower
    from . import registry as R
    say = log or (lambda m: print(m, flush=True))
    t0 = time.time()
    from .api import target_spec
    bloq, kind, params = build_target(target, n, p, w=w)
    spec = target_spec(target, bloq, params)
    rng = random.Random(seed)
    P0 = lower(bloq)
    say(f'lowered {target} {params}: {len(P0.ops)} ops ({time.time() - t0:.0f} s)')
    cands = joint_canaries(P0, spec, rng) + gate_canaries(P0, rng)
    rows = []
    for name, family, fn, spec_over in cands:
        if only and not any(name.startswith(o) for o in only):
            continue
        P = clone_program(P0)
        t1 = time.time()
        try:
            desc = fn(P)
        except Exception as e:  # noqa: BLE001
            rows.append(dict(name=name, family=family, status='SKIPPED', detail=f'mutation failed: {e}'))
            continue
        if desc == 'skipped':
            rows.append(dict(name=name, family=family, status='SKIPPED', detail='not applicable'))
            continue
        sp = spec_over if spec_over is not None else spec
        E = Engine(P, seed=seed)
        try:
            r = E.prove(P.root, sp)
            st, det = r['status'], (r.get('detail') or '')[:300]
        except Exception as e:  # noqa: BLE001
            st, det = 'ERROR', f'{type(e).__name__}: {e}'[:300]
        row = dict(name=name, family=family, mutation=desc, status=st, detail=det, seconds=round(time.time() - t1, 1))
        if family == 'gate':
            wt = _witness(P, spec, rng, E)
            row['witness'] = dict(mismatches=wt.get('mismatches'), cases=wt.get('cases'), first=wt.get('first'),
                                  level=wt.get('level'))
            row['counts'] = bool(wt.get('mismatches')) and st != 'PROVED'
            if not wt.get('mismatches') and st != 'PROVED':
                row['note'] = 'no witness found by the simulation (possibly equivalent or rare); not counted'
        else:
            row['counts'] = st != 'PROVED'
        rows.append(row)
        say(f"{name}: {st} {('witness' if row.get('witness', {}).get('mismatches') else '')} {det[:140]}")
    proved = [r['name'] for r in rows if r['status'] == 'PROVED']
    counted = [r for r in rows if r.get('counts')]
    res = dict(schema='qoptct.auto.canaries/1', target=target, params=params, versions=_versions(),
               qualtran_variant=_variant(), source_sha256=source_digest(), program_ops_sha256=P0.meta['ops_sha256'],
               total=len(rows),
               counted=len(counted), joints=sum(1 for r in counted if r['family'] == 'joint'),
               gates=sum(1 for r in counted if r['family'] == 'gate'), proved=proved, rows=rows,
               ok=(not proved and len(counted) >= min_count), seconds=round(time.time() - t0, 1))
    res['not_proved'] = sum(1 for r in rows if r['status'] not in ('PROVED', 'SKIPPED'))
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(res, indent=1, default=str))
    return res


def canaries_child(target, n, p=None, *, w=1, only=None, out=None, timeout=None):
    from ..qualtran._variant import run_json
    from .cli import variant_of
    args = ['canaries', target, '--n', str(n), '--w', str(w), '--inprocess'] + (['--p', str(p)] if p else []) + \
        (['--only', only] if only else [])
    res = run_json(variant_of(target), 'qoptct.foreign.auto', args, timeout=timeout)
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(res, indent=1, default=str))
    return res
