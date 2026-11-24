"""check_projectq: prove a ProjectQ program against a specification.

    from qoptct.foreign.projectq import check_projectq, ClassicalSpec, GF2MulSpec

    def program(eng):                       # any ProjectQ code, here: a quantum adder
        from projectq.ops import All, Measure
        from projectq.libs.math import AddQuantum
        a = eng.allocate_qureg(3); b = eng.allocate_qureg(3)
        AddQuantum | (a, b)
        return dict(inputs=dict(a=a, b=b), outputs=dict(a=a, b=b))

    r = check_projectq(program, ClassicalSpec(lambda v: dict(a=v['a'], b=(v['a'] + v['b']) % 8)))
    r['status']        # 'PROVED' | 'REFUTED' | 'UNKNOWN'

What PROVED means
-----------------
The program is run on a ProjectQ engine whose decomposition rules bring it into NCT (``mode='nct'``) or Clifford+T
(``mode='clifford_t'``); every command is translated by exactly one rule (``_logical``).  For a classical program
(``ClassicalSpec``) all 2^n inputs are simulated bit-parallel; for ``GF2MulSpec`` the exact prover ``exact2`` proves the
product for all inputs at any size.  Every wire that is not a specified output (and not a pinned garbage wire) must
end in |0>.  A Clifford+T program is proved equal to its NCT version by the exact path-sum prover.  PROVED additionally
requires the independent run of the genuine program on ProjectQ's ClassicalSimulator to agree (when the program is given
as a function).  REFUTED comes with a witness input.  UNKNOWN: anything the translation or the provers do not cover.

Assumptions recorded in every result (``assumptions``): the exact atom matrices (no global phase), the rules of
``_logical`` and the AND gadget lemma of ``qoptct.verify.exact2``.
"""
from __future__ import annotations

from ._record import Recording, Unsupported, ProjectQUnavailable, record, need_projectq
from ._logical import Logical, to_logical
from ._classical import ClassicalSpec, check_classical, pin_classical
from ._gf2 import GarbagePin, check_gf2_mul, pin_gf2, mulmod

__all__ = ['check_projectq', 'ClassicalSpec', 'GF2MulSpec', 'GarbagePin', 'pin_garbage', 'summarize']

ASSUMPTIONS = [
    'the matrices of the atoms X, Z, H, S, Sdg, T, Tdg are compared exactly, without a global phase',
    'one explicit rule per command (qoptct.foreign.projectq._logical); ProjectQ decomposition rules are not audited, '
    'their output is what is proved',
    'the AND gadget and its measurement-based erasure are the lemma checked by qoptct.verify.exact2.gadget_check',
    'a freed wire is required to end in |0> like every other clean wire',
]


class GF2MulSpec:
    """Specification of a GF(2^m) multiplier: inputs x, y (m qubits each), outputs x, y unchanged and result = x*y mod f."""

    def __init__(self, m, poly):
        self.m, self.poly = m, poly

    def as_classical(self):
        p = self.poly
        return ClassicalSpec(lambda v: {'x': v['x'], 'y': v['y'], 'result': mulmod(v['x'], v['y'], p)}, name=f'gf2^{self.m}')


def _circuit_info(rec, lg):
    return dict(mode=rec.mode, wires=rec.n_wires, logical_wires=lg.nq, logical_ops=len(lg.ops), counts=rec.counts(),
                recording_sha256=rec.sha256(), logical_sha256=lg.sha256(), projectq=rec.meta.get('projectq'))


def _as_logical(obj):
    if isinstance(obj, Logical):
        return None, obj
    if isinstance(obj, Recording):
        return obj, to_logical(obj)
    raise TypeError('expected a program, a Recording or a Logical')


def pin_garbage(program, spec, *, mode='nct', max_controls=2, max_input_bits=20):
    """Pin of the garbage of the golden ``program`` (a function, Recording or Logical) for ``spec``."""
    rec, lg = (None, None)
    if callable(program):
        rec = record(program, mode, max_controls)
        lg = to_logical(rec)
    else:
        rec, lg = _as_logical(program)
    if isinstance(spec, GF2MulSpec):
        return pin_gf2(lg, spec.m)
    return pin_classical(lg, max_input_bits)


def _widths(lg):
    return {k: len(v) for k, v in lg.inputs.items()}


def check_projectq(program, spec, *, mode='nct', max_controls=2, pin=None, max_input_bits=20, native=True, native_samples=24,
                   budget_sec=300.0, seed=7):
    """Prove ``program`` (a function ``program(eng)``, a :class:`Recording` or a :class:`Logical`) against ``spec``.

    ``pin``: a :class:`GarbagePin` of the golden program when the circuit leaves unspecified values on some wires
    (see ``pin_garbage``); without a pin every wire that is not a specified output must end in |0>.
    ``native``: also compare the genuine program with ProjectQ's ClassicalSimulator (needs ``program`` to be a function).
    -> dict(status, method, detail, circuit, independent, assumptions, witness?)."""
    res = dict(status='UNKNOWN', detail='', assumptions=list(ASSUMPTIONS))
    try:
        if callable(program):
            need_projectq()
            rec = record(program, 'nct' if mode == 'clifford_t' else mode, None if mode == 'clifford_t' else max_controls)
            lg = to_logical(rec)
        else:
            rec, lg = _as_logical(program)
            mode = lg.mode
    except ProjectQUnavailable:
        raise
    except Unsupported as e:
        return dict(res, detail=f'unsupported: {e}')
    except Exception as e:  # noqa: BLE001 - never an exception to the caller; UNKNOWN is never a proof
        return dict(res, detail=f'{type(e).__name__}: {e}')
    try:
        if mode == 'clifford_t':
            return _check_ct(program, spec, rec, lg, res, max_input_bits, native, native_samples, budget_sec, seed)
        r = _check_nct(program, spec, rec, lg, pin, max_input_bits, native, native_samples, seed)
    except Unsupported as e:
        return dict(res, detail=f'unsupported: {e}')
    except Exception as e:  # noqa: BLE001
        return dict(res, detail=f'{type(e).__name__}: {e}')
    res.update(r)
    if rec is not None:
        res['circuit'] = _circuit_info(rec, lg)
    return res


def _native(program, spec, lg, native, samples, seed):
    if not native or not callable(program):
        return dict(ok=None, tested=0, detail='not run (the program is not a function)' if not callable(program) else 'switched off')
    from ._native import native_check
    return native_check(program, spec.as_classical() if isinstance(spec, GF2MulSpec) else spec, _widths(lg),
                        count=samples, seed=seed)


def _check_nct(program, spec, rec, lg, pin, max_input_bits, native, samples, seed):
    if isinstance(spec, GF2MulSpec):
        r = check_gf2_mul(lg, spec.m, spec.poly, pin)
        domain = 'all inputs'
        cspec = spec
    elif isinstance(spec, ClassicalSpec):
        r = check_classical(lg, spec, pin, max_input_bits)
        domain = r.get('domain')
        cspec = spec
    else:
        raise TypeError('spec must be a ClassicalSpec or a GF2MulSpec')
    r = dict(r)
    r['domain'] = domain
    if r['status'] == 'PROVED' or (r['status'] == 'UNKNOWN' and r.get('exact2', {}).get('status') == 'PROVED'):
        nat = _native(program, cspec, lg, native, samples, seed)
        r['native'] = nat
        if nat['ok'] is False:
            r['status'] = 'UNKNOWN'
            r['detail'] = 'proved by the provers, but ProjectQ ClassicalSimulator disagrees (not trusted): ' + nat['detail']
    elif r['status'] == 'REFUTED' and callable(program) and native and isinstance(r.get('witness'), dict) and \
            set(r['witness']) == set(lg.inputs) and all(isinstance(v, int) for v in r['witness'].values()):
        from ._native import confirm_decomposition
        got = confirm_decomposition(program, r['witness'])
        want = (cspec.as_classical() if isinstance(cspec, GF2MulSpec) else cspec).fn(r['witness'])
        r['native'] = dict(ok=None, tested=1, decomposed_run=got, specification=want,
                           detail='ProjectQ decomposition run on the witness: '
                                  + ('the same wrong value (the decomposition itself differs from the specification)'
                                     if got is not None and any(got.get(k) != v for k, v in want.items()) else 'not confirmed'))
    return r


def _check_ct(program, spec, rec, lg_nct, res, max_input_bits, native, samples, budget_sec, seed):
    """Clifford+T: (1) the NCT version (multi-controlled X) is proved against ``spec``, (2) the Clifford+T version is proved
    equal to it.  Both must be PROVED."""
    if not callable(program):
        return dict(res, detail='clifford_t mode needs the program as a function')
    from ._ct import check_clifford_t
    ref = record(program, 'nct', None)
    lg_ref = to_logical(ref)
    first = _check_nct(program, spec, ref, lg_ref, None, max_input_bits, native, samples, seed)
    rct = record(program, 'clifford_t')
    lct = to_logical(rct)
    second = check_clifford_t(lct, ref, budget_sec=budget_sec)
    st = 'PROVED' if first['status'] == 'PROVED' and second['status'] == 'PROVED' else (
        'REFUTED' if 'REFUTED' in (first['status'], second['status']) else 'UNKNOWN')
    return dict(res, status=st, method='nct-spec + pathsum(clifford_t = nct)', detail=f"nct level: {first['status']} ({first['detail']}); "
                f"clifford_t = nct: {second['status']} ({second['detail']})", nct_level=first, clifford_t_level=second,
                circuit=_circuit_info(rct, lct), witness=first.get('witness') or second.get('witness'),
                domain=first.get('domain'), independent=first.get('native'))


def summarize(r):
    """One line for a result of check_projectq."""
    return f"{r['status']}: {r.get('detail', '')}"[:400]
