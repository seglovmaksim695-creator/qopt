# SPDX-License-Identifier: Apache-2.0
"""Direct proof of one node (a leaf of the proof tree) against its contract.

The operations of the node are cut out of the program, renumbered canonically (``montgomery._core.compact``, the
function the 0.7 proofs used, so the hashes of equal pieces are equal) and proved with the UNCHANGED qoptct prover:
``bridge.strict_proof.prove_program`` = kprove (ROBDD path sum, every input and every measurement outcome, exact
phase) + the exact partial-domain guard of the AND / ERASE gates.  The specification is built from the contract
expressions (``expr.bdd_word``): every output wire gets the bit of the exact value, every other wire must be 0, and
an output value that does not fit its register on the domain is a failed proof.

Before the prover, an independent simulation (``bridge.sim``, phases of every outcome) runs the same program on
sample inputs from the domain (edge values, constants of the contract and their neighbours, random values; models of
the precondition found by Z3 when random sampling misses it) and compares ALL wires with ``Contract.evaluate``.
Nodes without inputs or without outputs (allocations, basis states, effects) are proved by exhaustive simulation
instead (at most ``max_exhaustive_bits`` input bits); so is a node whose domain has at most ``small_domain`` points
(the domain is enumerated exactly by Z3: the precondition with every found point excluded is unsat), e.g. the
measurement-based uncomputation of a table lookup, whose precondition ``target = data[selection]`` makes the ROBDD of
the phase grow exponentially for kprove.  Exhaustive simulation of every point of the domain with the phase
coefficient of every measurement is a complete check of the classical action and of the phase of every outcome.
"""
from __future__ import annotations

import random
import time
from hashlib import sha256

from . import expr as X
from .contract import ContractError

__all__ = ['prove_leaf', 'sample_inputs', 'simulate', 'enumerate_domain', 'LeafError']


class LeafError(ValueError):
    pass


def _constants(c):
    out = set()
    stack = [c.pre] + [e for _, e in c.defs] + list(c.out.values())
    seen = set()
    while stack:
        e = stack.pop()
        if type(e) is int:
            out.add(e)
            continue
        if id(e) in seen:
            continue
        seen.add(id(e))
        if e[0] == 'lut':
            out.update(e[2])
            stack.append(e[1])
        elif e[0] in ('var', 'ref', 'val', 'sym', 'fold'):
            continue
        else:
            stack.extend(a for a in e[1:] if not isinstance(a, str))
    return out


def sample_inputs(contract, rng, count, *, use_z3=True, tries=None):
    """Up to ``count`` distinct input dictionaries inside the domain of the contract (edge values first)."""
    widths = contract.inputs
    if not widths:
        return [{}] if contract.in_domain({}) else []
    consts = sorted({v for c in _constants(contract) for v in (c - 1, c, c + 1) if v >= 0})
    pool = {}
    for k, w in widths.items():
        top = (1 << w) - 1
        vals = {0, 1 & top, top, (1 << (w - 1)) if w > 1 else 1, top - 1 if top else 0}
        vals |= {v for v in consts if v <= top}
        vals |= {(top - v) for v in consts if v <= top}
        pool[k] = sorted(vals)
    out, seen = [], set()

    def push(vals):
        key = tuple(sorted(vals.items()))
        if key not in seen and contract.in_domain(vals):
            seen.add(key)
            out.append(vals)

    names = list(widths)
    push({k: 0 for k in names})
    push({k: (1 << widths[k]) - 1 for k in names})
    for _ in range(tries or 40 * count):
        if len(out) >= count:
            break
        push({k: (rng.choice(pool[k]) if rng.random() < 0.5 else rng.getrandbits(widths[k])) for k in names})
    if len(out) < count and use_z3 and contract.pre != 1:
        for vals in _z3_models(contract, rng, count - len(out), seen):
            if contract.in_domain(vals):
                out.append(vals)
    return out


def _z3_models(contract, rng, count, seen):
    try:
        import z3
    except Exception:
        return []
    widths = contract.inputs
    vs = {k: z3.Int('in_' + k) for k in widths}
    refs = {}

    def leaf(lf):
        if lf[0] == 'var':
            return vs[lf[1]]
        return T.term(refs[lf[1]])

    def lb(lf):
        if lf[0] == 'var':
            return (0, (1 << widths[lf[1]]) - 1)
        return T.bounds(refs[lf[1]])
    T = X.Z3Terms(z3, leaf, lb)
    try:
        for nm, e in contract.defs:
            refs[nm] = e
        pre = T.cond(contract.pre)
    except X.ExprError:
        return []
    out = []
    names = list(widths)
    for _ in range(4 * count):
        if len(out) >= count:
            break
        s = z3.Solver()
        s.set('timeout', 5000)
        s.add(pre)
        for k, w in widths.items():
            s.add(vs[k] >= 0, vs[k] < (1 << w))
        s.push()
        for k in names:
            if rng.random() < 0.5:
                s.add(vs[k] == rng.getrandbits(widths[k]))
        if s.check() != z3.sat:
            s.pop()
            for k in names:
                if rng.random() < 0.3:
                    s.add(vs[k] >= rng.getrandbits(max(1, widths[k] - 1)))
            if s.check() != z3.sat:
                continue
        m = s.model()
        vals = {k: m.eval(vs[k], model_completion=True).as_long() for k in names}
        key = tuple(sorted(vals.items()))
        if key not in seen:
            seen.add(key)
            out.append(vals)
    return out


def enumerate_domain(contract, limit=256, timeout_ms=20000):
    """Every input of the domain of ``contract`` when there are at most ``limit`` of them: a list of dicts, exact
    (Z3: the precondition with all found points blocked is unsat), or None (more points, or Z3 undecided)."""
    if contract.pre == 1:
        total = 1
        for w in contract.inputs.values():
            total <<= w
        if total > limit:
            return None
    try:
        import z3
    except Exception:
        return None
    widths = contract.inputs
    vs = {k: z3.Int('in_' + k) for k in widths}
    refs = dict(contract.defs)

    def leaf(lf):
        return vs[lf[1]] if lf[0] == 'var' else T.term(refs[lf[1]])

    def lb(lf):
        return (0, (1 << widths[lf[1]]) - 1) if lf[0] == 'var' else T.bounds(refs[lf[1]])
    T = X.Z3Terms(z3, leaf, lb)
    try:
        pre = T.cond(contract.pre)
    except X.ExprError:
        return None
    s = z3.Solver()
    s.set('timeout', timeout_ms)
    s.add(pre)
    for k, w in widths.items():
        s.add(vs[k] >= 0, vs[k] < (1 << w))
    out = []
    while True:
        r = s.check()
        if r == z3.unsat:
            return out
        if r != z3.sat or len(out) >= limit:
            return None
        m = s.model()
        vals = {k: m.eval(vs[k], model_completion=True).as_long() for k in widths}
        out.append(vals)
        s.add(z3.Or([vs[k] != vals[k] for k in widths]) if widths else z3.BoolVal(False))


def _bits(value, wires):
    v = 0
    for j, q in enumerate(wires):
        v |= ((value >> j) & 1) << q
    return v


def _read(state, wires):
    return sum(((state >> q) & 1) << j for j, q in enumerate(wires))


def simulate(nq, cops, wi, wo, contract, samples):
    """-> (cases, first failure or None).  Every wire is compared; phases of all outcomes must be 0."""
    from ..qualtran.bridge import sim
    out_wires = {q for ws in wo.values() for q in ws}
    for vals in samples:
        init = 0
        for k, ws in wi.items():
            init |= _bits(vals[k], ws)
        want = contract.evaluate(vals)
        try:
            got, ph = sim.run(nq, cops, init)
        except sim.SimulationError as e:
            return len(samples), dict(input=vals, why=f'not executable: {e}')
        exp = 0
        for k, ws in wo.items():
            exp |= _bits(want[k], ws)
        if got != exp or any(ph):
            rest = got
            for q in out_wires:
                rest &= ~(1 << q)
            gotv = {k: _read(got, ws) for k, ws in wo.items()}
            return len(samples), dict(input=vals, expected=want, got=gotv, phase=any(ph),
                                      other_wires_nonzero=bool(rest))
    return len(samples), None


def _spec_functions(contract, wi, wo):
    from qoptct.verify.bdd import ZERO
    widths = contract.inputs

    def words(S):
        w = {}
        for k, ws in wi.items():
            r = S.reg(ws)
            w[k] = X.BddWord(list(r.bits) + [ZERO], 0, (1 << widths[k]) - 1)
        return w

    def evaluator(S):
        ws = words(S)
        defs = {}

        def leaf(x):
            if x[0] == 'var':
                return ws[x[1]]
            if x[0] == 'ref':
                return defs[x[1]]
            raise LeafError(f'leaf {x[0]} cannot appear in a leaf contract')
        for nm, e in contract.defs:
            defs[nm] = X.bdd_word(e, S.B, leaf)
        return leaf

    def pre_fn(S):
        leaf = evaluator(S)
        r = X.bdd_word(contract.pre, S.B, leaf)
        if r.lo < 0 or r.hi > 1:
            raise LeafError('the precondition is not 0/1')
        return r.bits[0]

    def spec(S):
        B = S.B
        leaf = evaluator(S)
        prew = None
        if contract.pre != 1:
            pw = X.bdd_word(contract.pre, B, leaf)
            if pw.lo < 0 or pw.hi > 1:
                raise LeafError('the precondition is not 0/1')
            prew = pw.bits[0]
        out = {}
        for k, e in contract.out.items():
            wv = X.bdd_word(e, B, leaf)
            width = contract.outputs[k]
            bits = wv.ext(max(len(wv.bits), width + 1))
            for b in bits[width:]:
                if b != ZERO and (prew is None or B.AND(prew, b) != ZERO):
                    raise LeafError(f'output {k}: its value does not fit {width} bits on the domain '
                                    f'(the contract must truncate explicitly)')
            for j, q in enumerate(wo[k]):
                out[q] = bits[j]
        return out
    return spec, (pre_fn if contract.pre != 1 else None)


def prove_leaf(ops, ins, outs, contract, *, budget=600, max_mb=1500, seed=11, cases=24, label='',
               max_exhaustive_bits=16, small_domain=256):
    """Prove ``ops`` (the operations of one node) against ``contract``.

    ins / outs: {register: LSB-first physical wires} of the node, in the order of the source signature.
    Returns dict(status = PROVED | REFUTED | UNKNOWN | ERROR, method, qasm_sha256, layout_sha256, contract_sha256,
    sim_cases, witness, detail, ...).  Nothing but 'PROVED' is a proof."""
    from ..qualtran.montgomery._core import compact, make_layout
    from ..qualtran.bridge.emit import emit_qasm
    from ..qualtran.bridge.strict_proof import prove_program, digest_json
    t0 = time.time()
    res = dict(label=label, status='ERROR', method=None, contract=contract.name)
    try:
        contract.validate()
        res['contract_sha256'] = contract.digest()
        for side, regs, widths in (('inputs', ins, contract.inputs), ('outputs', outs, contract.outputs)):
            if {k: len(v) for k, v in regs.items()} != widths:
                raise LeafError(f'{side}: node registers {{{", ".join(f"{k}:{len(v)}" for k, v in regs.items())}}} '
                                f'do not match the contract {widths}')
        nq, cops, wi, wo, _mp = compact(list(ops), dict(ins), dict(outs))
        layout = make_layout(nq, wi, wo, cops)
        qasm = emit_qasm(nq, cops)
        qh = sha256(qasm.encode()).hexdigest()
        res.update(nq=nq, ops=len(cops), qasm_sha256=qh, layout_sha256=digest_json(layout))
        rng = random.Random(seed)
        in_bits = sum(contract.inputs.values())
        if not contract.inputs or not contract.outputs:
            if in_bits > max_exhaustive_bits:
                res.update(status='UNKNOWN', detail=f'no inputs or no outputs and {in_bits} input bits: '
                                                    f'too many to enumerate')
                return res
            samples = []
            names = list(contract.inputs)
            for code in range(1 << in_bits):
                vals, sh = {}, 0
                for k in names:
                    vals[k] = (code >> sh) & ((1 << contract.inputs[k]) - 1)
                    sh += contract.inputs[k]
                if contract.in_domain(vals):
                    samples.append(vals)
            if not samples:
                res.update(status='UNKNOWN', detail='the domain is empty')
                return res
            cnt, bad = simulate(nq, cops, wi, wo, contract, samples)
            res.update(sim_cases=cnt, method='exhaustive_simulation_all_outcomes_with_phase')
            if bad:
                res.update(status='REFUTED', witness=bad, detail='exhaustive simulation differs from the contract')
            else:
                res.update(status='PROVED', detail=f'all {cnt} inputs of the domain, every outcome, no phase')
            res['wall'] = round(time.time() - t0, 3)
            return res
        if contract.pre != 1 and small_domain:
            dom = enumerate_domain(contract, limit=small_domain)
            if dom is not None:
                dom = [d for d in dom if contract.in_domain(d)] if dom else dom
                if not dom:
                    res.update(status='UNKNOWN', detail='the domain is empty')
                    return res
                cnt, bad = simulate(nq, cops, wi, wo, contract, dom)
                res.update(sim_cases=cnt, method='exhaustive_domain_simulation_all_outcomes_with_phase',
                           domain_points=cnt, domain_exact='z3: precondition with every found point blocked is unsat')
                if bad:
                    res.update(status='REFUTED', witness=bad, detail='simulation of a domain point differs')
                else:
                    res.update(status='PROVED', detail=f'every one of the {cnt} inputs of the domain (enumerated '
                                                       f'exactly by Z3), every outcome, no phase')
                res['wall'] = round(time.time() - t0, 3)
                return res
        samples = sample_inputs(contract, rng, cases)
        if not samples:
            res.update(status='UNKNOWN', detail='no sample input inside the domain was found')
            return res
        cnt, bad = simulate(nq, cops, wi, wo, contract, samples)
        res['sim_cases'] = cnt
        if bad:
            res.update(status='REFUTED', witness=bad, method='independent_simulation',
                       detail='the independent simulation differs from the contract', wall=round(time.time() - t0, 3))
            return res
        cdict = dict(source='qoptct.foreign.auto leaf contract', label=label, contract=contract.to_dict())
        receipt = dict(status='PASS', qasm_sha256=qh, contract_sha256=digest_json(cdict),
                       layout_sha256=digest_json(layout), independent_inputs=cnt, all_wires_checked=True,
                       phase_checked=True)
        spec, pre = _spec_functions(contract, wi, wo)
        v = prove_program(qasm, layout, spec, pre=pre, validation=receipt, contract=cdict, budget=budget,
                          max_mb=max_mb)
        st = v.get('status')
        res.update(status=st if st in ('PROVED', 'REFUTED', 'UNKNOWN') else 'ERROR', method='kprove_bdd_exact_domains',
                   detail=(v.get('detail') or '')[:300], prover_sha256=v.get('prover_sha256'),
                   domain=(v.get('partial_domains') or {}).get('status'), witness=v.get('witness'),
                   peak_mb=v.get('peak_rss_mb'))
    except (LeafError, ContractError, X.ExprError) as e:
        res.update(status='UNKNOWN', detail=f'{type(e).__name__}: {e}')
    res['wall'] = round(time.time() - t0, 3)
    return res
