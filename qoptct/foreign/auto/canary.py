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
    if multi:
        out.append(('J1_swap_wires_of_a_child', 'joint', swap_in(multi[0].path), None))
        if len(multi) > 2:
            out.append(('J1b_swap_wires_of_a_later_child', 'joint', swap_in(multi[len(multi) // 2].path, 1), None))

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

    # a composite child: an operation between two of its children
    comp = None
    for k in big:
        sub = [c for c in _proper_children(P0, k) if c.nops]
        if len(sub) >= 2:
            comp = (k, sub)
            break
    if comp is not None:
        k, sub = comp

        def extra_between_sub(P):
            pos = sub[0].span[1]
            q = sub[0].out_wires()[0]
            insert_op(P, pos, ('x', q), k.path)
            return f'X on wire {q} inserted between {sub[0].path} and {sub[1].path} (owned by {k.path})'
        out.append(('J6_extra_operation_between_grandchildren', 'joint', extra_between_sub, None))

        def move_boundary(P):
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

    # phase: erase corrections
    for lf in leaves:
        er = positions(lf, 'erase')
        if er:
            i = er[len(er) // 2]

            def drop_corr(P, i=i):
                op = P.ops[i]
                P.ops[i] = ('erase', op[1], ())
                P._support = {}
                return f'op {i}: erase correction removed (phase)'
            out.append((f'G1_phase_correction_removed_{lf.type}', 'gate', drop_corr, None))

            def move_corr(P, i=i, lf=lf):
                op = P.ops[i]
                (cz, a, b), = op[2]
                others = [q for q in lf.in_wires() if q not in (a, b, op[1])]
                if not others:
                    return 'skipped'
                P.ops[i] = ('erase', op[1], (('cz', a, others[0]),))
                P._support = {}
                return f'op {i}: erase correction moved to wire {others[0]} (phase)'
            out.append((f'G2_phase_correction_moved_{lf.type}', 'gate', move_corr, None))
            break
    # rare inputs: a carry CX near the top of an adder
    for tname in ('CAdd', 'Add', 'Controlled', 'AddK', 'LessThanConstant'):
        for lf in by_type.get(tname, [])[:1]:
            cx = positions(lf, 'cx')
            ands = positions(lf, 'and')
            if len(cx) > 4:
                i = cx[-3]

                def flip_cx(P, i=i):
                    op = P.ops[i]
                    P.ops[i] = ('cx', op[2], op[1]) if op[1] != op[2] else op
                    P._support = {}
                    return f'op {i}: control and target of a CX near the top of the carry chain exchanged'
                out.append((f'G3_carry_cx_reversed_{tname}', 'gate', flip_cx, None))
            if len(ands) > 3:
                i = ands[-2]

                def drop_and(P, i=i):
                    op = P.ops[i]
                    P.ops[i] = ('cx', op[1], op[3])
                    P._support = {}
                    return f'op {i}: the AND of a high carry replaced by a CX of one control (wrong only on rare inputs)'
                out.append((f'G4_high_carry_and_weakened_{tname}', 'gate', drop_and, None))
    # an X dropped (constants are loaded with X)
    for lf in leaves:
        xs = positions(lf, 'x')
        if len(xs) >= 2:
            i = xs[len(xs) // 2]

            def drop_x(P, i=i):
                del P.ops[i]
                _shift_after_delete(P, i)
                return f'op {i}: an X removed'
            out.append((f'G5_x_removed_{lf.type}', 'gate', drop_x, None))
            break
    # two neighbouring operations of a leaf exchanged (they must not commute)
    for lf in leaves:
        a, b = lf.span
        for i in range(a, b - 1):
            o1, o2 = ops[i], ops[i + 1]
            if o1[0] == 'cx' and o2[0] == 'cx' and (o1[2] == o2[1] or o2[2] == o1[1]):
                def swap2(P, i=i):
                    P.ops[i], P.ops[i + 1] = P.ops[i + 1], P.ops[i]
                    P._support = {}
                    return f'ops {i}, {i + 1}: two non-commuting CX exchanged'
                out.append((f'G6_noncommuting_swap_{lf.type}', 'gate', swap2, None))
                break
        else:
            continue
        break
    # the last operation of the program removed (the final reduction)
    def drop_last(P):
        i = len(P.ops) - 1
        del P.ops[i]
        _shift_after_delete(P, i)
        return f'op {i} (the last one) removed'
    out.append(('G7_last_operation_removed', 'gate', drop_last, None))
    # random single-operation mutants inside random leaves
    pool = [lf for lf in leaves if lf.nops >= 3]
    for j in range(max(0, count - len(out))):
        if not pool:
            break
        lf = rng.choice(pool)
        i = rng.randrange(lf.span[0], lf.span[1])

        def rnd(P, i=i, lf=lf, seed=rng.random()):
            r = random.Random(seed)
            op = P.ops[i]
            wires = sorted(set(lf.in_wires()) | set(lf.out_wires()))
            new = op
            for _ in range(20):
                if op[0] == 'cx':
                    t = r.choice(wires)
                    new = ('cx', op[1], t) if t != op[1] else op
                elif op[0] == 'x':
                    new = ('x', r.choice(wires))
                elif op[0] == 'and':
                    a = r.choice(wires)
                    new = ('and', a, op[2], op[3]) if a not in (op[2], op[3]) else op
                else:
                    new = ('erase', op[1], ())
                if new != op:
                    break
            if new == op:
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
