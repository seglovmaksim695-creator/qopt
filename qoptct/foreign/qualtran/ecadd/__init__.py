"""Proof of point addition for Qualtran's ``ECAdd`` and ``ECWindowAddR`` (compositional, exact, for ALL inputs).

Public API
----------
``prove_ecadd(n, p, *, patched=True)``
    Six steps of ``ECAdd`` (each proved by the ROBDD prover of this package with the Montgomery multiplier and the
    Kaliski inverse as boxes with contracts) + gluing of the six step contracts in Z3, with the wiring between
    the steps read from the real ``ECAdd.decompose_bloq()``.  ``patched=True`` is the fixed scheme (Qualtran 0.7.0 +
    flag fix E05 + PR 1997 + exact-integer table; applied to a COPY of Qualtran, never to the installed one) and
    must come out as the point sum.  ``patched=False`` is the original 0.7.0 scheme: it is proved correct exactly
    outside the groups G1-G3 of finding 05 (status ``PROVED_OUTSIDE_G1_G3``).
``prove_ecwindowaddr(n, p, window=1)``
    ``ECWindowAddR`` (fixed scheme): ``(c, Q) -> (c, Q + c*R)``; the table lookup and its removal are proved gate by
    gate, the ``ECAdd`` box is taken with the contract of ``prove_ecadd``.
``contracts_status(n, p)``
    For which (n, p) the contracts of the multiplier and the inverse are proved, for which only assumed.

Every verdict is a dict with ``status`` (``PROVED`` / ``NOT_PROVED``), the per-part results, the Qualtran variant, the
Z3 version, and ``contracts_status``.  A verdict is only ``PROVED`` when every part is ``PROVED``.
Heavy parts run in child processes against the patched copy of Qualtran (``qoptct.foreign.qualtran.patches``).
"""
from __future__ import annotations

__all__ = ['prove_ecadd', 'prove_ecwindowaddr', 'contracts_status', 'CONTRACTS']


def __getattr__(name):
    if name in ('prove_ecadd', 'prove_ecwindowaddr'):
        from . import api
        return getattr(api, name)
    if name in ('contracts_status', 'CONTRACTS'):
        from . import contracts
        return getattr(contracts, name)
    raise AttributeError(name)
