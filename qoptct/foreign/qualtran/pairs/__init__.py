"""Mirror check of "forward + backward = identity" pairs of Qualtran blocks.

``check_pair_kaliski(n, p)``
    ``KaliskiModInverse(n, p)`` followed by ``KaliskiModInverse(n, p, uncompute=True)``: the pair is the identity on
    all wires (exact phase, every measurement outcome), for ``0 <= x < p``.
``check_pair_montgomery(n, w, p)``
    ``DirtyOutOfPlaceMontgomeryModMul(n, w, p)`` followed by its ``uncompute=True`` twin (fixed Qualtran, kind ``'c'``).

Method (gate level, no semantics of the blocks is trusted): both halves are lowered gate by gate by the strict bridge to
``x / cx / and / erase``; the backward half must be the MIRROR of the forward half (reverse order, ``and`` <-> ``erase``)
up to an injective renaming of wires and swaps of commuting gates.  Because every gate is reversible on the classical
values it acts on, a mirror program followed by its original returns every wire to its start value.  For large n the
check runs per macro node (window / iteration) and the assembly is compared element by element with the full lowering at
small n (``compare_full``).  The forward block itself is NOT verified here: its contract is a separate result
(``qoptct.foreign.qualtran.ecadd.contracts``).

A verdict is a dict: ``status`` ``PROVED`` only if the mirror holds (and, when requested, the assembly equals the full
lowering), ``NOT_PROVED`` otherwise.  Work runs in a child process on the right Qualtran variant (Kaliski: installed
0.7.0; Montgomery: patched kind ``'c'``).
"""
from __future__ import annotations

__all__ = ['check_pair_kaliski', 'check_pair_montgomery']

_PKG = 'qoptct.foreign.qualtran.pairs'


def _verdict(res, kind, extra_required=()):
    from .._variant import expected_kind_name
    want = expected_kind_name(kind)
    ok = res.get('mirror_ok') is True and res.get('qualtran_kind') == want
    for k in extra_required:
        ok = ok and res.get(k) is True
    out = dict(res)
    out['status'] = 'PROVED' if ok else 'NOT_PROVED'
    if res.get('qualtran_kind') != want:
        out['detail'] = f"child ran on Qualtran variant {res.get('qualtran_kind')!r}, expected {want!r}"
    return out


def check_pair_kaliski(n, p, *, compare_full=None, timeout=None):
    """Mirror check of the Kaliski pair.  ``compare_full``: also compare with the full lowering (default: n <= 16)."""
    from .._variant import run_json
    from ..domain import check_domain
    dom = check_domain('pair-kaliski', n, p)
    if not dom['ok']:
        return dict(status='NOT_PROVED', api='check_pair_kaliski', n=n, p=str(p), domain=dom,
                    detail='outside the domain of the proof: ' + '; '.join(dom['problems']))
    full = (n <= 16) if compare_full is None else bool(compare_full)
    args = ['run-kaliski', '--n', n, '--p', str(p)] + (['--full'] if full else [])
    res = run_json(None, _PKG, args, timeout=timeout)
    if res.get('status') in ('ERROR', 'UNKNOWN'):
        return dict(res, api='check_pair_kaliski', n=n, p=str(p))
    out = _verdict(res, None, ('equals_full_lowering',) if full else ())
    out.update(api='check_pair_kaliski', n=n, p=str(p), compare_full=full)
    return out


def check_pair_montgomery(n, w, p, *, compare_full=None, timeout=None):
    """Mirror check of the Montgomery multiplier pair (window ``w``).  ``compare_full``: default n <= 8."""
    from .._variant import run_json
    from ..domain import check_domain
    dom = check_domain('pair-montgomery', n, p)
    if not dom['ok']:
        return dict(status='NOT_PROVED', api='check_pair_montgomery', n=n, w=w, p=str(p), domain=dom,
                    detail='outside the domain of the proof: ' + '; '.join(dom['problems']))
    full = (n <= 8) if compare_full is None else bool(compare_full)
    args = ['run-montgomery', '--n', n, '--w', w, '--p', str(p)] + (['--full'] if full else [])
    res = run_json('c', _PKG, args, timeout=timeout)
    if res.get('status') in ('ERROR', 'UNKNOWN'):
        return dict(res, api='check_pair_montgomery', n=n, w=w, p=str(p))
    out = _verdict(res, 'c', ('equals_full_lowering',) if full else ())
    out.update(api='check_pair_montgomery', n=n, w=w, p=str(p), compare_full=full)
    return out
