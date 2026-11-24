"""Evidence runs: prove one block, run its canaries, and write a result dictionary (the files of results/projectq/).

Every result is reproduced by the command in its ``reproduce`` field.  A result is registered (``registry``) only if the
proof is PROVED, every canary that is not equivalent to the golden program is NOT proved, and at least ``MIN_CANARIES``
canaries were decided.
"""
from __future__ import annotations

import collections
import time

__all__ = ['MIN_CANARIES', 'prove_multiplier', 'prove_point_addition', 'prove_library']

MIN_CANARIES = 20


def _canary_block(rec, spec, pin, count, seed):
    from .canaries import run_canaries
    rows = run_canaries(rec, spec, pin, count, seed)
    decided = [r for r in rows if r['status'] != 'EQUIVALENT']
    groups = collections.Counter(r['group'] for r in decided)
    targeted = sum(groups[g] for g in ('junction', 'erasure', 'phase', 'rare'))
    false_proved = [r for r in decided if r['status'] == 'PROVED']
    return dict(requested=count, decided=len(decided), equivalent=len(rows) - len(decided), not_proved=len(decided) - len(false_proved),
                false_proved=false_proved, by_group=dict(groups), targeted=targeted,
                by_status=dict(collections.Counter(r['status'] for r in decided)),
                ok=(not false_proved and len(decided) >= MIN_CANARIES and 2 * targeted >= len(decided)), rows=rows)


def _base(rid, name, block, claim, r, rec, lg, t0, reproduce, canaries):
    c = r.get('circuit') or {}
    return dict(schema='qoptct.projectq.result/1', id=rid, name=name, block=block, claim=claim, status=r['status'], detail=r.get('detail'),
                method=r.get('method'), domain=r.get('domain'), projectq=c.get('projectq'),
                program=dict(name='logical_sha256', value=lg.sha256()), recording_sha256=rec.sha256(),
                circuit=dict(mode=rec.mode, wires=rec.n_wires, logical_wires=lg.nq, logical_ops=len(lg.ops), counts=rec.counts()),
                independent=r.get('independent'), native=r.get('native'), exact2=r.get('exact2'), canaries=canaries,
                assumptions=r['assumptions'], reproduce=reproduce, wall_seconds=round(time.time() - t0, 1))


def prove_multiplier(m, source=None, canaries=24, seed=1):
    from . import binary_ecc as B
    from .check import check_projectq, pin_garbage
    from ._logical import to_logical
    from ._record import record
    t0 = time.time()
    prog, spec = B.multiplier_program(m, source), B.multiplier_spec(m)
    rec = record(prog)
    lg = to_logical(rec)
    pin = pin_garbage(lg, spec)
    r = check_projectq(prog, spec, pin=pin)
    can = _canary_block(rec, spec, pin, canaries, seed) if r['status'] == 'PROVED' else None
    out = _base(f'binary_ecc_multiplier/m{m}', 'binary_ecc_multiplier', f'Binary_ECC Multiplication/Mul_x{m}.py',
                f"result = x*y mod f in GF(2^{m}), f = {hex(spec.poly)}, x and y unchanged, every wire that is not x, y, result or garbage ends in |0>, "
                f"the {len(pin)} garbage wires end in the values of the pin", r, rec, lg, t0,
                f'python -m qoptct.foreign.projectq prove multiplier --m {m} --source <Binary_ECC folder>', can)
    out.update(m=m, poly=hex(spec.poly), source=dict(file=f'Multiplication/Mul_x{m}.py', sha256=B.file_sha256(B.multiplier_path(m, source))),
               garbage=dict(wires=len(pin), pin_sha256=pin.digest, kind=pin.kind,
                            note='garbage wires are compared with the golden circuit, not proved meaningful'))
    return out


def prove_point_addition(m, source=None, canaries=24, seed=1, variant='FLT_in'):
    from . import binary_ecc as B
    from .check import check_projectq, pin_garbage
    from ._logical import to_logical
    from ._record import record
    t0 = time.time()
    prog = B.point_addition_program(m, variant, source)
    spec = B.point_addition_spec(m, prog)
    rec = record(prog)
    lg = to_logical(rec)
    pin = pin_garbage(lg, spec)
    r = check_projectq(prog, spec, pin=pin, native_samples=48)
    can = _canary_block(rec, spec, pin, canaries, seed) if r['status'] == 'PROVED' else None
    cst = prog.constants
    out = _base(f'binary_ecc_point_addition/{variant}_m{m}', 'binary_ecc_point_addition',
                f'Binary_ECC Point_additions/{variant}/x^{m}/Point_addition_x{m}.py',
                f'q = 0: (x1, y1) unchanged; q = 1 and x1 != x2: (x1, y1) <- (x1, y1) + (x2, y2) on y^2 + xy = x^3 + a x^2 + b over GF(2^{m}); '
                f'constants a={cst["a"]}, x2={hex(cst["x2"])}, y2={hex(cst["y2"])}; all {1 << (2 * m + 1)} inputs; clean wires clean; '
                f'{len(pin)} garbage wires equal the pin', r, rec, lg, t0,
                f'python -m qoptct.foreign.projectq prove point-addition --m {m} --source <Binary_ECC folder>', can)
    out.update(m=m, poly=hex(sum(1 << e for e in B.POLYNOMIALS[m])), constants=cst,
               source=dict(file=f'Point_additions/{variant}/x^{m}/Point_addition_x{m}.py',
                           sha256=B.file_sha256(B.point_addition_path(m, variant, source))),
               garbage=dict(wires=len(pin), pin_sha256=pin.digest, kind=pin.kind,
                            note='garbage wires are compared with the golden circuit, not proved meaningful'))
    return out


def prove_library(name, n, canaries=24, seed=1, mode='nct'):
    from .library import block
    from .check import check_projectq
    from ._logical import to_logical
    from ._record import record
    t0 = time.time()
    prog, spec = block(name, n)
    rec = record(prog)
    lg = to_logical(rec)
    r = check_projectq(prog, spec, mode=mode)
    can = _canary_block(rec, spec, None, canaries, seed) if r['status'] == 'PROVED' and mode == 'nct' else None
    out = _base(f'projectq_library/{name}_n{n}' + ('' if mode == 'nct' else '_clifford_t'), 'projectq_library', f'projectq.libs.math: {name}',
                f'{spec.name}: all outputs equal the arithmetic specification' + ('' if spec.promise is None else ' on the promised domain')
                + ', every other wire ends in |0>', r, rec, lg, t0, f'python -m qoptct.foreign.projectq prove library --block {name} --n {n}'
                + ('' if mode == 'nct' else ' --clifford-t'), can)
    out.update(n=n, block_name=name, source=None, garbage=dict(wires=0, pin_sha256=None, kind=None, note='no garbage'))
    return out
