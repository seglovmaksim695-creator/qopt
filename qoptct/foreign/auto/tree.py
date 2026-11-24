# SPDX-License-Identifier: Apache-2.0
"""The library-independent view of a lowered program: a flat list of logical operations plus a tree of nodes.

A library adapter (``qualtran_adapter`` today; a Qiskit or ProjectQ adapter would build the same objects from the
instruction hierarchy) produces a ``Program``:

* ``ops``  -- the logical operations of the qoptct language (``('x', t)``, ``('cx', a, t)``, ``('and', a, b, t)``,
  ``('erase', t, (('cz', a, b),))``), in program order;
* ``root`` -- a ``Node`` whose descendants own contiguous, nested ranges of ``ops``.

A ``Node`` records what the source library says about one call: its type, its parameters (``key``), the physical
wires of every input and output register (LSB first), and the range of operations ``span`` that the lowering
attributed to it (``None``: no operations).  Children are produced lazily by the adapter.  ``kind`` is
``'op'`` (an ordinary call), ``'wires'`` (pure renaming: split, join, cast), ``'alloc'`` (fresh zero wires) or
``'free'`` (wires that must be zero and are dropped).  Nothing here is trusted beyond the wiring of the source: every
claim a node makes is checked (spans by the coverage check, interfaces by the proofs of the children).
"""
from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ['Node', 'Program', 'op_wires', 'support']


def op_wires(op):
    k = op[0]
    if k == 'x':
        return (op[1],)
    if k == 'cx':
        return (op[1], op[2])
    if k == 'and':
        return (op[1], op[2], op[3])
    if k == 'erase':
        ws = [op[1]]
        for c in op[2]:
            ws.extend(c[1:])
        return tuple(ws)
    raise ValueError(f'unknown operation {op!r}')


def support(ops):
    s = set()
    for op in ops:
        s.update(op_wires(op))
    return s


@dataclass
class Node:
    path: str
    type: str
    kind: str                      # op | wires | alloc | free
    ins: dict                      # contract input name -> LSB-first physical wires
    outs: dict                     # contract output name -> LSB-first physical wires
    span: tuple | None             # [a, b) of Program.ops, or None
    key: str = ''                  # type + parameters (equal keys: equal calls)
    obj: object = None             # the library object (adapter only)
    depth: int = 0
    info: dict = field(default_factory=dict)

    @property
    def nops(self):
        return 0 if self.span is None else self.span[1] - self.span[0]

    def in_wires(self):
        return [q for ws in self.ins.values() for q in ws]

    def out_wires(self):
        return [q for ws in self.outs.values() for q in ws]

    def widths_in(self):
        return {k: len(v) for k, v in self.ins.items()}

    def widths_out(self):
        return {k: len(v) for k, v in self.outs.items()}


class Program:
    """ops + root + an adapter that knows how to list children.  ``adapter.children(node)`` -> list of Node in the
    order of the source (instance order); ``adapter.contract(node)`` -> a registry Contract or None."""

    def __init__(self, ops, nq, root, adapter, *, meta=None):
        self.ops, self.nq, self.root, self.adapter = ops, nq, root, adapter
        self.meta = dict(meta or {})
        self._support = {}

    def node_ops(self, node):
        return [] if node.span is None else self.ops[node.span[0]:node.span[1]]

    def node_support(self, node):
        k = node.span
        if k is None:
            return set()
        s = self._support.get(k)
        if s is None:
            s = support(self.ops[k[0]:k[1]])
            if k[1] - k[0] < 2_000_000:
                self._support[k] = s
        return s

    def children(self, node):
        return self.adapter.children(node)
