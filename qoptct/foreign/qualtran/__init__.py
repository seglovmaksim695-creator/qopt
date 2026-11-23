"""Qualtran -> qoptct: a strict gate-by-gate bridge and proof clients for Qualtran blocks.

Qualtran is an OPTIONAL dependency (``pip install "qoptct[qualtran]"``).  Importing this
package never imports Qualtran; the modules that need it import it lazily.

Main entry points
-----------------
``qoptct.foreign.qualtran.check_bloq(bloq, ...)``
    Lower a Qualtran bloq with the strict bridge and prove it against Qualtran's own
    classical action (``call_classically``) on the whole input domain, phase and all
    measurement outcomes included.  Verdict: PROVED / REFUTED (with a witness) / UNKNOWN.
``qoptct.foreign.qualtran.patches.patched_qualtran(kind)``
    A copy of the installed Qualtran with the patches the proofs were done with.
``qoptct.foreign.qualtran.bridge``
    The bridge itself (layout lowering, rules, strict proof, domain guard, rules A1 and C1).
"""
from ._cache import cache_dir

__all__ = ['cache_dir', 'check_bloq']


def __getattr__(name):
    if name == 'check_bloq':
        from .check import check_bloq
        return check_bloq
    raise AttributeError(name)
