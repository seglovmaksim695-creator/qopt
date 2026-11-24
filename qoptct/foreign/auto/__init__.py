# SPDX-License-Identifier: Apache-2.0
"""Generic hierarchical prover for large foreign circuits (``qoptct.foreign.auto``, new in 0.8).

The mechanism descends the decomposition of a block and proves it from its pieces, with no code written for the
block itself:

* a **leaf** (a primitive of the registry: addition, subtraction, comparison, copy, XOR, swap, table lookup, basis
  states ...) is proved on its own gates by the unchanged qoptct prover (kprove ROBDD + exact AND/ERASE domains,
  ``leaf``), or by exhaustive simulation when its domain is tiny;
* a **composite** node is proved by composing the contracts of its children in Z3 (linear integer arithmetic,
  ``compose``), with the checks the audits of round 3 asked for: coverage of every operation, freshness of every
  wire a child takes as zero, wire flow (every child reads its real wires), the domain of every child, the end state
  of every wire;
* a **loop** of identical children needs an invariant and a step given as data (``contract.Fold``); without them
  the loop is unrolled and, if Z3 cannot follow it, the result is UNKNOWN, never PROVED.

The core (``expr``, ``contract``, ``tree``, ``compose``, ``leaf``) does not import Qualtran; ``qualtran_adapter`` turns a
Qualtran bloq into the library-independent ``tree.Program`` (the 0.7 bridge does the lowering), ``registry`` holds the
contracts of Qualtran primitives, ``specs`` the specifications of composite blocks.  See ``docs/AUTO.md``.

    from qoptct.foreign.auto import prove_bloq
    v = prove_bloq(CModAdd(QUInt(256), 2**255 - 19))      # v['verdict'] == 'PROVED' or 'NOT_PROVED'

    python -m qoptct.foreign.auto prove cmodadd --n 256 --out results/auto/cmodadd_n256.json
"""
from __future__ import annotations

__all__ = ['prove_bloq', 'prove_target', 'Contract', 'Fold', 'Engine', 'TARGETS']


def __getattr__(name):
    if name in ('prove_bloq', 'TARGETS', 'build_target'):
        from . import api
        return getattr(api, name)
    if name == 'prove_target':
        from .cli import prove_target
        return prove_target
    if name in ('Contract', 'Fold'):
        from . import contract
        return getattr(contract, name)
    if name == 'Engine':
        from .compose import Engine
        return Engine
    raise AttributeError(name)
