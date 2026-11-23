"""Proof client for the Kaliski modular inverse of Qualtran (``KaliskiModInverse``).

``prove_kaliski(n, p)``  ->  verdict dictionary, ``status`` is ``'PROVED'`` or ``'NOT PROVED (<failed checks>)'``.

The circuit is lowered gate by gate with the strict bridge (``qoptct.foreign.qualtran.bridge``), composed from the real pieces
(``lower``), and proved with the unchanged qoptct Walker (``chain``, with the wire-flow check) under assumptions that Z3 closes
for all n, p, x (``lemmas_z3``); ``verdict`` demands every item of the criterion.  Needs Qualtran 0.7.0 and Z3.
Command line: ``python -m qoptct.foreign.qualtran.kaliski --help``.
"""

__all__ = ['prove_kaliski', 'judge']


def __getattr__(name):
    if name in __all__:
        from . import verdict
        return getattr(verdict, name)
    raise AttributeError(name)
