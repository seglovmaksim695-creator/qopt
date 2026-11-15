# SPDX-License-Identifier: Apache-2.0
"""qoptct.verify -- exact and fast verification of qoptct circuits against their .qc reference.

    prove(qasm_or_path, ref_qc_path, clean=False, budget_sec=120.0)
        -> {"status": "PROVED" | "REFUTED" | "UNKNOWN", "method": str, "detail": str, ...}

PROVED means: for EVERY basis input and EVERY combination of measurement outcomes, each branch of the circuit equals
c_m * (U_ref (x) |a(m)>) (strict model; clean model: on inputs whose clean wires are 0). It is established by the
exact path-sum prover (qoptct.verify.pathsum: symbolic outcomes, normal form, trace miter in Z[w]) or, for unitary
circuits, by MQT QCEC when it is installed (optional extra "qcec"). REFUTED comes with a counterexample when one is
found. UNKNOWN means "not proved" (time budget, width limit, unsupported gate) and must never be treated as PROVED.

Fast sampled checks (not proofs): check_dense (state vector, <= 26 qubits), check_sparse (independent sparse
simulator), check_big (c10 path-sum sampler for 30-70+ qubits, strict model), check_qcec, and qoptct.verify.arith for
arithmetic circuits at 64-2048 bits. quick_check picks one of them by size.
"""
import time

from ._parse import QasmFile, parse_qasm_text, parse_qc, is_unitary_program
from ._limits import run_limited
from .branchsim import check_sparse, confirm_branch
from .dense import check_dense
from .qcec import check_qcec, available as qcec_available

__all__ = ['prove', 'quick_check', 'check_dense', 'check_sparse', 'check_big', 'check_qcec', 'qcec_available',
           'confirm_branch', 'main_prove']

STATUSES = ('PROVED', 'REFUTED', 'UNKNOWN')


def _prove_child(path, qc, clean, budget, max_width):
    from .pathsum import prove_exact
    r = prove_exact(path, qc, clean=clean, budget_sec=budget, max_width=max_width)
    return {k: v for k, v in r.items() if k in (
        'verdict', 'method', 'trace', 'trace_meta', 'reason', 'error', 'counterexample', 'confirmed', 'confirm_dev',
        'confirm_error', 'time', 'n', 'N', 'meas', 'ys_after_rules', 'normalization')}


def _detail_of(r):
    parts = []
    if r.get('method') == 'normal-form' and r['verdict'] == 'PROVED':
        parts.append('path sum reduced to the normal form for all inputs and outcomes')
    if r.get('trace_meta') and isinstance(r['trace_meta'], dict):
        tm = r['trace_meta']
        if 'S' in tm:
            parts.append('trace miter S=%s, r^2=%s, width %s' % (tm.get('S'), tm.get('r2'), tm.get('width')))
        elif 'why' in tm:
            parts.append('trace miter: %s' % tm['why'])
    if r.get('reason'):
        parts.append('violations: %s' % (r['reason'],))
    if r.get('counterexample'):
        parts.append('counterexample %s%s' % (r['counterexample'],
                                              ' (confirmed by the sparse simulator)' if r.get('confirmed') else ''))
    if r.get('error'):
        parts.append(str(r['error']))
    return '; '.join(parts)


def prove(qasm_or_path, ref_qc_path, clean=False, budget_sec=120.0, use_qcec=True, max_width=22):
    """Exact equivalence check of a qoptct OpenQASM 3 circuit (text or path) against a .qc reference.
    budget_sec: wall-clock limit for the whole call (None = no limit, runs in-process). The prover runs in a separate
    Python process that is killed at the limit. max_width bounds the trace-miter tables (2^w * 32 bytes each;
    22 -> 128 MB). Never raises: bad input, time-outs and internal errors give UNKNOWN, never PROVED."""
    t0 = time.time()
    out = dict(status='UNKNOWN', method='none', detail='')
    try:
        with QasmFile(qasm_or_path) as path:
            try:
                txt = open(path).read()
                N, nb, ops = parse_qasm_text(txt)
                names, ins, ref = parse_qc(ref_qc_path)
            except (ValueError, KeyError, OSError) as e:
                out.update(method='parse', detail='cannot parse: %s' % e)
                return out
            n = len(names)
            unitary = is_unitary_program(ops) and N == n
            qc_ok = use_qcec and unitary and qcec_available()
            # the whole call stays within budget_sec: the prover gets 85% (55% if QCEC may follow), a hard kill
            # slightly later, QCEC the remainder
            pb = None if budget_sec is None else (budget_sec * (0.55 if qc_ok else 0.85))
            hard = None if pb is None else pb * 1.1 + 1.0
            r = run_limited(_prove_child, (path, ref_qc_path, clean, pb, max_width), budget_sec=hard)
            if r[0] == 'ok':
                pr = r[1]
                v = pr.get('verdict')
                out['proof'] = pr
                if v in ('PROVED', 'REFUTED'):
                    meth = {'normal-form': 'pathsum-normal-form', 'trace': 'pathsum-trace-miter',
                            'counterexample': 'pathsum-counterexample'}.get(pr.get('method'), 'pathsum')
                    out.update(status=v, method=meth, detail=_detail_of(pr))
                    if pr.get('counterexample'):
                        out['counterexample'] = pr['counterexample']
                else:
                    why = 'internal inconsistency (reported, not trusted)' if v == 'INTERNAL' else _detail_of(pr)
                    out.update(method='pathsum', detail='path-sum prover: %s %s' % (v, why))
            elif r[0] == 'timeout':
                out.update(method='pathsum', detail='path-sum prover: time budget exceeded')
            else:
                out.update(method='pathsum', detail='path-sum prover error: %s' % r[1])
            if out['status'] == 'UNKNOWN' and qc_ok:
                left = None if budget_sec is None else max(3.0, budget_sec - (time.time() - t0) - 3.0)
                q = check_qcec(txt, ref_qc_path, budget_sec=left if left is not None else 600.0)
                if q['status'] in ('PROVED', 'REFUTED'):
                    out.update(status=q['status'], method='qcec', detail=q['detail'])
                else:
                    out['detail'] = (out['detail'] + '; ' if out['detail'] else '') + q['detail']
    except FileNotFoundError as e:
        out.update(method='parse', detail='file not found: %s' % e)
    except Exception as e:      # never raise to the optimizer: an unexpected error is UNKNOWN, not PROVED
        out.update(method='error', detail='%s: %s' % (type(e).__name__, e))
    out['time'] = round(time.time() - t0, 2)
    return out


def check_big(qasm_or_path, ref_qc_path, x_trials=4, m_trials=4, seed=1, cap=2048):
    """Sampled path-sum check for large strict-model circuits (qoptct.contrib.c10.verify_big, v0.4 fixed version).
    -> dict(ok, F, sampled, ...). Not a proof."""
    from ..contrib.c10 import verify_big as VB
    with QasmFile(qasm_or_path) as path:
        try:
            return VB.check(path, ref_qc_path, x_trials=x_trials, m_trials=m_trials, seed=seed, verbose=False, cap=cap)
        except (AssertionError, ValueError, KeyError) as e:
            return dict(ok=False, F=None, error='%s: %s' % (type(e).__name__, e))


def quick_check(qasm_or_path, ref_qc_path, clean=False, seed=1):
    """Fast sampled check chosen by size: dense (<= 22 qubits), else sparse, else check_big (strict model only).
    -> dict(ok, F, engine, ...). Not a proof; use prove()."""
    with QasmFile(qasm_or_path) as path:
        txt = open(path).read()
        try:
            N, nb, ops = parse_qasm_text(txt)
        except ValueError as e:
            return dict(ok=False, F=None, engine='parse', error=str(e))
        if N <= 22:
            r = check_dense(path, ref_qc_path, clean=clean, seed=seed)
            r['engine'] = 'dense'
            return r
        r = check_sparse(path, ref_qc_path, clean=clean, seed=seed, cap=1 << 18)
        r['engine'] = 'sparse'
        if r.get('F') is not None or clean:
            return r
        r = check_big(path, ref_qc_path, seed=seed)
        r['engine'] = 'big'
        return r


def main_prove(argv=None):
    from .cli import main_prove as _m
    return _m(argv)
