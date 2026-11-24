# SPDX-License-Identifier: Apache-2.0
"""qoptct.loop -- the closed circle "propose -> prove -> accumulate" (see docs/LOOP.md).

Proposers suggest improved circuits for a target (a reference .qc); the unchanged provers of qoptct.verify check each
candidate against the reference; only proved improvements enter the journal (journal.jsonl, hash chain) and the
best-known table is built from the journal alone. `qoptct loop verify` re-checks everything.
"""
def main(argv=None):
    from .cli import main as _m
    return _m(argv)


__all__ = ['main']
