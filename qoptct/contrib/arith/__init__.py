# SPDX-License-Identifier: Apache-2.0
"""arith: recognition of classical reversible arithmetic in .qc circuits and re-synthesis on logical ANDs
(Gidney 2018: 4 T per AND, uncomputation by X-basis measurement + Clifford correction, 0 T).

    candidates(qc_path, clean, budget_sec=60.0) -> list of {"qasm": str, "note": str}

Pipeline: the .qc circuit is lowered to x/cx/ccx (H-frames), its full truth table is computed bit-sliced over all
inputs (up to MAX_VARS input bits), the function is recognised by probing (one- and two-bit inputs) and confirmed
on ALL inputs against the definition of the family (adder, out-of-place adder, comparator, addition mod 2^n-1,
hidden weighted bit), a circuit is generated from the definition, and every candidate is checked
  1. exactly on all inputs while it is built (values of all wires, every measurement-based erase),
  2. by an independent sparse simulation of the emitted OpenQASM text (forced measurement outcomes),
  3. by qoptct.verify.prove (if available) or the dense checker qoptct.verify_qasm (small circuits).
Only candidates passing all checks are returned; nothing is raised to the caller.

Models: clean=True — wires not on the .i line start in |0> and must return to |0> (the reference may leave
results on them); clean=False (strict) — all wires carry arbitrary data, only new qubits are clean.
"""
import os
import tempfile
import time

from . import tt
from .ir import BuildError, lower_ops, peephole, to_qasm, count, check_text
from . import fam, hwb

MAX_VARS = 24
DENSE_MAX_QUBITS = 18


def _families():
    return fam.FAMILIES + [('hwb', hwb.rec_hwb, hwb.variants_hwb)]


def _external(qasm, qc_path, clean, budget):
    """-> (ok, how). A REFUTED / failed external check rejects the candidate."""
    try:
        from qoptct import verify as _v      # qoptct >= 0.4
    except Exception:
        _v = None
    if _v is not None:
        try:
            r = _v.prove(qasm, qc_path, clean=clean, budget_sec=max(5.0, budget))
            st = r.get('status')
            if st == 'REFUTED':
                return False, 'verify.prove REFUTED: %s' % r.get('detail', '')
            if st == 'PROVED':
                return True, 'verify.prove PROVED (%s)' % r.get('method', '')
        except Exception:
            pass
        try:
            r = _v.check_sparse(qasm, qc_path, clean=clean, trials=4)
            if r.get('ok'):
                return True, 'verify.prove UNKNOWN; verify.check_sparse F=%.12f' % r.get('F', 0.0)
            if r.get('F') is not None:
                return False, 'verify.check_sparse failed F=%s' % r.get('F')
            # simulator limits (no verdict): fall through to the remaining checks
        except Exception:
            pass
    N = int(qasm.split('qubit[', 1)[1].split(']', 1)[0])
    if N <= DENSE_MAX_QUBITS:
        try:
            from qoptct import verify_qasm
            fd, p = tempfile.mkstemp(suffix='.qasm')
            with os.fdopen(fd, 'w') as f:
                f.write(qasm)
            try:
                F = verify_qasm.check(p, qc_path, trials=3, clean=clean)
            finally:
                os.unlink(p)
            if F < 1 - 1e-9:
                return False, 'verify_qasm F=%.3g' % F
            return True, 'verify_qasm F=1-%.1e' % (1 - F)
        except Exception as e:
            return False, 'verify_qasm error %r' % (e,)
    return True, 'internal checks only'


def build_all(qc_path, clean, budget_sec=60.0, external=True):
    """-> list of dicts with qasm, note, T, CX, N, extra (all checks passed)."""
    t0 = time.time()
    T = tt.load(qc_path, clean, MAX_VARS)
    if T is None:
        return []
    res = []
    seen = set()
    for name, rec, variants in _families():
        if time.time() - t0 > budget_sec:
            break
        try:
            hyps = rec(T)
        except Exception:
            hyps = []
        for h in hyps:
            for vname, make in variants(T, h):
                if time.time() - t0 > budget_sec:
                    break
                try:
                    bd = make()
                    bd.finish()
                except Exception:
                    continue
                body, nb = lower_ops(bd.ops)
                body = peephole(body)
                qasm = to_qasm(body, bd.N, nb, T.plus)
                Tc, CX = count(qasm)
                key = (Tc, CX, bd.N)
                if key in seen:
                    continue
                seen.add(key)
                res.append(dict(qasm=qasm, body=body, nb=nb, T=Tc, CX=CX, N=bd.N, extra=bd.N - T.n, family=name,
                                nand=bd.nand, variant=vname))
    res = [r for r in res if not any(_dominates(s, r) for s in res)]
    good = []
    for r in res:
        ok, msg = check_text(r.pop('body'), r.pop('nb'), T)
        if not ok:
            continue
        r['note'] = 'arith/%s/%s: %d AND, T=%d CX=%d, +%d qubits (%s model); exact check on all %d inputs; ' \
                    'text simulation %s' % (r['family'], r['variant'], r['nand'], r['T'], r['CX'], r['extra'],
                                            'clean' if clean else 'strict', 1 << T.k, msg)
        good.append(r)
    res = good
    if external:
        good = []
        for i, r in enumerate(res):
            left = budget_sec - (time.time() - t0)
            ok, how = _external(r['qasm'], qc_path, clean, max(5.0, left / (len(res) - i)))
            if ok:
                r['note'] += '; ' + how
                good.append(r)
        res = good
    return res


def _dominates(s, r):
    a = (s['T'], s['CX'], s['extra'])
    b = (r['T'], r['CX'], r['extra'])
    return a != b and all(x <= y for x, y in zip(a, b))


def candidates(qc_path, clean, budget_sec=60.0):
    try:
        return [{'qasm': r['qasm'], 'note': r['note']} for r in build_all(qc_path, clean, budget_sec)]
    except Exception:
        return []
