# SPDX-License-Identifier: Apache-2.0
"""``prove_bloq``: the whole pipeline for one Qualtran bloq, and the named targets of the command line."""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

from . import expr as X
from .compose import Engine
from .contract import ContractError
from .leaf import sample_inputs

__all__ = ['prove_bloq', 'TARGETS', 'build_target', 'end_to_end', 'spec_selfcheck', 'SCHEMA', 'ASSUMPTIONS']

SCHEMA = 'qoptct.auto/1'

ASSUMPTIONS = [
    'the strict bridge of 0.7 (lower_bloqgraph + translate) and its rules A1 (controlled And) and C1 '
    '(QROAMCleanAdjoint as a primitive with a contract, assumptions A1-A4 of bridge/qroam_clean_adjoint.py); '
    'TwoBitCSwap, where it occurs, by the opt-in rule bridge/fredkin_rule.py (CX.CCX.CX)',
    'the qoptct prover kprove (ROBDD path sum) with the exact AND/ERASE domain guard, and the gadget lemma for AND/erase',
    'Z3 (linear integer arithmetic) for the composition of the contracts of the children',
    'the specification of the root is the statement (formulas in qoptct/foreign/auto/specs.py); it is compared with '
    'call_classically of the same Qualtran on a small instance (spec_check of the verdict), not proved equal to it',
]


def _versions():
    from importlib.metadata import version, PackageNotFoundError
    out = {'python': sys.version.split()[0]}
    for nm in ('qualtran', 'z3-solver', 'numpy'):
        try:
            out[nm] = version(nm)
        except PackageNotFoundError:
            out[nm] = None
    try:
        import z3
        out['z3'] = z3.get_version_string()
    except Exception:
        out['z3'] = None
    return out


# files of this package that do not take part in a verdict (the canary generator, the registry of results and the
# cross-check with 0.7): a change there does not change what a verdict means; canary files record canary.py apart
_NOT_IN_VERDICT = ('canary.py', 'results.py', 'crosscheck.py')


def source_digest():
    """SHA-256 of the sources a verdict depends on: this package (without the canary generator and the results
    registry), the bridge and the qoptct provers"""
    from hashlib import sha256
    import os
    here = Path(os.path.abspath(__file__)).parent
    pkg = here.parent
    files = sorted(f for f in here.glob('*.py') if f.name not in _NOT_IN_VERDICT) + \
        sorted((pkg / 'qualtran' / 'bridge').glob('*.py')) + \
        sorted((pkg.parent / 'verify').rglob('*.py')) + [pkg / 'qualtran' / 'montgomery' / '_core.py']
    h = sha256()
    for f in files:
        h.update(str(f.relative_to(pkg.parent)).encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def _variant():
    try:
        from ..qualtran._variant import detect_variant
        return detect_variant()
    except Exception as e:  # noqa: BLE001
        return dict(kind='unknown', error=f'{type(e).__name__}: {e}')


def end_to_end(P, spec, *, cases=16, seed=7, max_ops=3_000_000):
    """Independent check: run the WHOLE lowered program (bridge.sim, phases of every outcome) on sample inputs of the
    domain and compare every wire with the specification (Python integers)."""
    from ..qualtran.bridge import sim
    if len(P.ops) > max_ops:
        return dict(skipped=f'{len(P.ops)} operations > {max_ops}')
    root = P.root
    rng = random.Random(seed)
    samples = sample_inputs(spec, rng, cases)
    bad, first = 0, None
    t0 = time.time()
    out_wires = {q for ws in root.outs.values() for q in ws}
    for vals in samples:
        init = 0
        for k, ws in root.ins.items():
            for j, q in enumerate(ws):
                init |= ((vals[k] >> j) & 1) << q
        want = spec.evaluate(vals)
        try:
            got, ph = sim.run(P.nq, P.ops, init)
        except sim.SimulationError as e:
            bad += 1
            first = first or dict(input=vals, why=str(e)[:200])
            continue
        exp = 0
        for k, ws in root.outs.items():
            for j, q in enumerate(ws):
                exp |= ((want[k] >> j) & 1) << q
        if got != exp or any(ph):
            bad += 1
            if first is None:
                rest = got
                for q in out_wires:
                    rest &= ~(1 << q)
                first = dict(input={k: str(v) for k, v in vals.items()}, phase=bool(any(ph)), other_wires=bool(rest),
                             outputs_differ=[k for k, ws in root.outs.items()
                                             if sum(((got >> q) & 1) << j for j, q in enumerate(ws)) != want[k]])
    return dict(cases=len(samples), mismatches=bad, first=first, seconds=round(time.time() - t0, 1))


def _summary(E):
    """distinct pieces that were proved: type, how, program hash, operations, how many times used"""
    rows = {}
    for r in E.records:
        k = (r.get('type'), r.get('mode'), r.get('qasm_sha256'), r.get('contract_name'))
        d = rows.setdefault(k, dict(type=k[0], mode=k[1], qasm_sha256=k[2], contract=k[3], ops=r.get('ops'),
                                    status=r.get('status'), method=r.get('method'), count=0))
        d['count'] += 1
        if r.get('status') != 'PROVED':
            d['status'] = r.get('status')
            d['detail'] = (r.get('detail') or '')[:300]
    return sorted(rows.values(), key=lambda d: (d['mode'] != 'compose', -(d['ops'] or 0)))


def prove_bloq(bloq, spec=None, *, target=None, params=None, expect_variant=None, e2e_cases=16, seed=7,
               leaf_budget=600, leaf_max_mb=1500, z3_timeout_ms=120_000, log=None, mutate=None, out=None,
               e2e_max_ops=3_000_000, spec_check=None):
    """Prove ``bloq`` against ``spec`` (default: the registry specification of its type) with the generic
    compositional mechanism.  Returns the verdict dict; ``verdict['verdict']`` is ``'PROVED'`` only if the root
    composition is proved, the independent end-to-end simulation agrees and Z3 is present.  ``mutate(P)`` damages
    the lowered program (canaries only)."""
    from . import registry as R
    from .qualtran_adapter import lower
    t0 = time.time()
    say = log or (lambda m: None)
    v = dict(schema=SCHEMA, target=target, params=params or {}, verdict='NOT_PROVED', reasons=[],
             versions=_versions(), qualtran_variant=_variant(), source_sha256=source_digest(),
             assumptions=list(ASSUMPTIONS))
    if spec_check is not None:
        v['spec_check'] = spec_check
        if not spec_check.get('ok'):
            v['reasons'].append(f"the specification differs from call_classically on a small instance: "
                                f"{str(spec_check.get('first_mismatch') or spec_check.get('why'))[:300]}")
    if expect_variant is not None and v['qualtran_variant'].get('kind') != expect_variant:
        v['reasons'].append(f"Qualtran variant {v['qualtran_variant'].get('kind')!r}, the proof needs {expect_variant!r}")
    try:
        spec = spec if spec is not None else R.contract_for(bloq)
    except Exception as e:  # noqa: BLE001
        spec = None
        v['reasons'].append(f'no specification: {type(e).__name__}: {e}')
    if spec is None:
        v['reasons'].append('no specification for this bloq (the root needs one)')
        return _finish(v, t0, out)
    try:
        spec.validate()
    except (ContractError, X.ExprError) as e:
        v['reasons'].append(f'invalid specification: {e}')
        return _finish(v, t0, out)
    v['spec'] = dict(name=spec.name, sha256=spec.digest(), note=spec.note)
    if v['versions'].get('z3') is None:
        v['reasons'].append('z3 is not available (the verdict is never PROVED without it)')
        return _finish(v, t0, out)
    say('lowering ...')
    t1 = time.time()
    try:
        P = lower(bloq, log=say)
    except Exception as e:  # noqa: BLE001
        v['reasons'].append(f'the bridge refused the bloq: {type(e).__name__}: {str(e)[:400]}')
        return _finish(v, t0, out)
    v['program'] = dict(P.meta)
    v['lower_seconds'] = round(time.time() - t1, 1)
    if mutate is not None:
        v['mutation'] = mutate(P)
        v['program']['mutated_ops_sha256'] = __import__('qoptct.foreign.auto.qualtran_adapter',
                                                        fromlist=['ops_digest']).ops_digest(P.ops)
    E = Engine(P, z3_timeout_ms=z3_timeout_ms, leaf_budget=leaf_budget, leaf_max_mb=leaf_max_mb, log=say, seed=seed)
    t2 = time.time()
    r = E.prove(P.root, spec)
    v['prove_seconds'] = round(time.time() - t2, 1)
    v['root'] = {k: r.get(k) for k in ('status', 'mode', 'detail', 'failure', 'obligations', 'z3_calls', 'deferred',
                                       'children', 'folds', 'qasm_sha256', 'layout_sha256', 'seconds')}
    v['pieces'] = _summary(E)
    v['stats'] = dict(E.stats, z3_seconds=round(E.z3_seconds, 2))
    if r['status'] != 'PROVED':
        v['reasons'].append(f"composition: {r['status']}: {(r.get('detail') or '')[:600]}")
    say('end-to-end simulation ...')
    e2e = end_to_end(P, spec, cases=e2e_cases, seed=seed, max_ops=e2e_max_ops)
    v['end_to_end'] = e2e
    if e2e.get('mismatches'):
        v['reasons'].append(f"end-to-end simulation differs from the specification on {e2e['mismatches']} of "
                            f"{e2e['cases']} inputs")
        v['witness'] = e2e.get('first')
    if not v['reasons']:
        v['verdict'] = 'PROVED'
    return _finish(v, t0, out)


def _finish(v, t0, out):
    v['seconds'] = round(time.time() - t0, 1)
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(v, indent=1, default=str))
    return v


# ------------------------------------------------------------------------------------------------ named targets
def _p_default(n):
    if n == 256:
        return 2 ** 255 - 19
    p = (1 << (n - 1)) - 1
    from .specs import is_probable_odd_prime
    while not is_probable_odd_prime(p):
        p -= 2
    return p


def _k_default(p):
    """the constant of the *K targets when none is given: 3^101 mod p (a fixed, unremarkable residue)"""
    return pow(3, 101, p)


def build_target(name, n, p=None, w=1, cv=1, k=None):
    """(bloq, expected Qualtran variant, params) for a named target"""
    from qualtran import QUInt
    from qualtran.bloqs import mod_arithmetic as M
    from qualtran.bloqs import arithmetic as A
    if name in ('subtract', 'subtractfrom', 'negate'):
        cls = {'subtract': A.Subtract, 'subtractfrom': A.SubtractFrom, 'negate': A.Negate}[name]
        return cls(QUInt(n)), None, dict(n=n)
    p = p or _p_default(n)
    params = dict(n=n, p=p)
    if name in ('cmodaddk', 'ctrlscalemodadd', 'cmodmulk'):
        k = _k_default(p) if k is None else k
        params['k'] = k
        if name == 'cmodaddk':
            return M.CModAddK(k=k, mod=p, bitsize=n), None, params
        if name == 'ctrlscalemodadd':
            return M.CtrlScaleModAdd(k=k, mod=p, bitsize=n), None, params
        return M.CModMulK(QUInt(n), k=k, mod=p), None, params
    if name == 'montgomery':
        from qualtran.bloqs.mod_arithmetic.mod_multiplication import DirtyOutOfPlaceMontgomeryModMul
        params['w'] = w
        return DirtyOutOfPlaceMontgomeryModMul(bitsize=n, window_size=w, mod=p), 'c', params
    if name == 'modadd':
        return M.ModAdd(n, p), None, params
    if name == 'cmodadd':
        params['cv'] = cv
        return M.CModAdd(QUInt(n), p, cv=cv), None, params
    if name == 'modsub':
        return M.ModSub(QUInt(n), p), None, params
    if name == 'cmodsub':
        params['cv'] = cv
        return M.CModSub(QUInt(n), p, cv=cv), None, params
    if name == 'modneg':
        return M.ModNeg(QUInt(n), p), None, params
    if name == 'cmodneg':
        params['cv'] = cv
        return M.CModNeg(QUInt(n), p, cv=cv), None, params
    if name == 'moddbl':
        return M.ModDbl(QUInt(n), p), None, params
    if name.startswith('ecadd_step'):
        from qualtran.bloqs.cryptography.ecc.ec_add import ECAdd
        k = int(name[len('ecadd_step'):])
        steps = [bi.bloq for bi in sorted(ECAdd(n, p, w).decompose_bloq().bloq_instances, key=lambda b: b.i)
                 if type(bi.bloq).__name__.startswith('_ECAddStep')]
        return steps[k - 1], 'e', params
    raise ValueError(f'unknown target {name}')


TARGETS = ('montgomery', 'modadd', 'cmodadd', 'modsub', 'cmodsub', 'modneg', 'cmodneg', 'moddbl', 'ecadd_step1',
           'subtract', 'subtractfrom', 'negate', 'cmodaddk', 'ctrlscalemodadd', 'cmodmulk')

# the classical action of Qualtran's DirtyOutOfPlaceMontgomeryModMul is defined (and returns zeros otherwise) only
# for 0 < x, y < p; the specification is compared with it there (it also asserts the closed form x*y*2^-n mod p)
_CLASSICAL_DOMAIN = {'montgomery': lambda v, p: 0 < v['x'] < p and 0 < v['y'] < p}


def spec_selfcheck(target, n=None, *, w=1, cv=1, max_inputs=1 << 16):
    """Compare the specification of a named target with ``call_classically`` of the Qualtran in use, on every input
    of its domain at a small width (n = 6 for Montgomery, 4 otherwise).  Must run with the target's Qualtran variant."""
    import itertools
    from . import registry as R
    n = n or ((6 + w - 1) // w * w if target == 'montgomery' else 4)
    bloq, _kind, params = build_target(target, n, None, w=w, cv=cv)
    spec = target_spec(target, bloq, params)
    res = dict(target=target, params={k: str(v) for k, v in params.items()}, checked=0, ok=False, first_mismatch=None)
    lefts, rights = list(bloq.signature.lefts()), list(bloq.signature.rights())
    names = [nm for r in lefts for nm in R.contract_names(r)]
    if set(names) != set(spec.inputs):
        res['why'] = f'specification inputs {sorted(spec.inputs)} != registers {names}'
        return res
    total = 1
    for k in names:
        total *= 1 << spec.inputs[k]
    if total > max_inputs:
        res['why'] = f'{total} inputs > {max_inputs}'
        return res
    dom = _CLASSICAL_DOMAIN.get(target)
    for vals in itertools.product(*[range(1 << spec.inputs[k]) for k in names]):
        v = dict(zip(names, vals))
        if not spec.in_domain(v) or (dom is not None and not dom(v, params['p'])):
            continue
        want = spec.evaluate(v)
        try:
            out = bloq.call_classically(**{r.name: R.to_classical(r, v) for r in lefts})
        except Exception as e:  # noqa: BLE001
            res['first_mismatch'] = dict(input=v, why=f'call_classically raised {type(e).__name__}: {e}'[:300])
            return res
        out = out if isinstance(out, tuple) else (out,)
        got = {}
        for r, val in zip(rights, out):
            got.update(R.from_classical(r, val))
        res['checked'] += 1
        if got != want:
            res['first_mismatch'] = dict(input=v, specification=want, qualtran=got)
            return res
    res['ok'] = res['checked'] > 0
    return res


def target_spec(name, bloq, params):
    """the specification of a named target (the registry one, or a specification by name)"""
    from . import registry as R
    from . import specs as S
    if name == 'ecadd_step1':
        return S.ecadd_step1_spec(params['n'], params['p'])
    return R.contract_for(bloq)
