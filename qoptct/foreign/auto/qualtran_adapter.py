# SPDX-License-Identifier: Apache-2.0
"""Qualtran adapter: the lowered program and the tree of calls of a Qualtran bloq, for ``qoptct.foreign.auto``.

The program is lowered exactly as in the 0.7 proofs: ``bridge.qualtran_layout.lower_bloqgraph`` (the real
``decompose_bloq`` connections, fresh wires, rules A1 controlled-And and C1 ``QROAMCleanAdjoint``; the basis-state
atoms ``ZeroState``/``OneState``/``ZeroEffect``/``OneEffect`` of ``bridge.basis_vectors``) and
``bridge.rules.translate``.  The hash of the program is ``montgomery._stages.ops_digest`` of 0.7.

The tree comes from the graph events of the lowering (path ``Root/i:Type/j:Type``, the physical wires of every register
of every call) and from the trace of the translation (which operations each call produced).  A node's operations must
form one contiguous range (checked here and refused by the prover otherwise).  Registers are turned into the
integers of the contracts like this: a scalar register is an unsigned integer (LSB-first wires); a one-dimensional bit
array is an integer whose element 0 is the most significant bit (as ``bb.join`` / ``QDType.from_bits`` of Qualtran);
any other shaped register is split into elements ``name[i,j]``.
"""
from __future__ import annotations

import gc
from hashlib import sha256

from .tree import Node, Program

__all__ = ['lower', 'QualtranAdapter', 'AdapterError', 'register_words']

WIRES = ('Split', 'Join', 'Cast', 'Partition')


class AdapterError(ValueError):
    pass


def register_words(regs, flat):
    """{contract name: LSB-first wires} for the registers ``regs`` (signature order) from the event's flat lists"""
    out = {}
    for reg in regs:
        if reg.name not in flat:
            raise AdapterError(f'register {reg.name} missing in the graph event')
        ws = [int(q) for q in flat[reg.name]]
        shape = tuple(int(s) for s in reg.shape)
        b = int(reg.bitsize)
        size = b
        for s in shape:
            size *= s
        if len(ws) != size:
            raise AdapterError(f'register {reg.name}: {len(ws)} wires, expected {size}')
        if not shape or (len(shape) == 1 and b == 1):
            out[reg.name] = list(reversed(ws))
        else:
            for k, idx in enumerate(reg.all_idxs()):
                out[reg.name + '[' + ','.join(map(str, idx)) + ']'] = list(reversed(ws[k * b:(k + 1) * b]))
    if len({q for v in out.values() for q in v}) != sum(len(v) for v in out.values()):
        raise AdapterError('aliased wires in one side of a call')
    return out


def _spans(trace):
    """path -> ((start, end), contiguous) for every prefix of every leaf path that produced operations"""
    spans = {}
    stack = []

    def close(entry):
        pf, st, cnt, end = entry
        if pf in spans:
            a, b, c, _ = spans[pf]
            spans[pf] = (min(a, st), max(b, end), c + cnt, False)
        else:
            spans[pf] = (st, end, cnt, True)
    prev = 0
    for t in trace:
        a, b = t['logical_start'], t['logical_end']
        if a != prev:
            raise AdapterError('the trace of the translation is not contiguous')
        prev = b
        if b == a:
            continue
        parts = t['leaf']['path'].split('/')
        prefixes, pf = [], ''
        for i, c in enumerate(parts):
            pf = c if i == 0 else pf + '/' + c
            prefixes.append(pf)
        k = 0
        while k < len(stack) and k < len(prefixes) and stack[k][0] == prefixes[k]:
            k += 1
        for entry in reversed(stack[k:]):
            close(entry)
        del stack[k:]
        for d in range(k, len(prefixes)):
            stack.append([prefixes[d], a, 0, a])
        for entry in stack:
            entry[2] += b - a
            entry[3] = b
    for entry in reversed(stack):
        close(entry)
    return {pf: ((a, b), ok and b - a == c) for pf, (a, b, c, ok) in spans.items()}


def ops_digest(ops):
    """the program hash of 0.7 (``montgomery._stages.ops_digest``)"""
    h = sha256()
    for o in ops:
        h.update(repr(o).encode())
    return h.hexdigest()


class QualtranAdapter:
    def __init__(self, bloq, low, tr, registry=None):
        self.root_bloq = bloq
        self.registry = registry
        self.events = {e['path']: e for e in low['graph_events']}
        if len(self.events) != len(low['graph_events']):
            raise AdapterError('repeated event paths')
        self.root_path = type(bloq).__name__
        if self.root_path not in self.events:
            raise AdapterError('no event for the root')
        self.kids = {}
        for p in self.events:
            if '/' in p:
                parent = p.rsplit('/', 1)[0]
                self.kids.setdefault(parent, []).append(p)
        for v in self.kids.values():
            v.sort(key=lambda q: int(q.rsplit('/', 1)[1].split(':', 1)[0]))
        self.spans = _spans(tr['trace'])
        self._bloq = {self.root_path: bloq}
        self._decomp = {}
        self._nodes = {}

    # ------------------------------------------------------------------ bloq objects of the paths
    def bloq_of(self, path):
        b = self._bloq.get(path)
        if b is not None:
            return b
        parent, comp = path.rsplit('/', 1)
        idx, tname = comp.split(':', 1)
        insts = self._instances(parent)
        b = insts.get(int(idx))
        if b is None or type(b).__name__ != tname:
            raise AdapterError(f'{path}: the decomposition of the parent does not give this call')
        self._bloq[path] = b
        return b

    def _instances(self, path):
        d = self._decomp.get(path)
        if d is None:
            from qualtran import CompositeBloq
            pb = self.bloq_of(path)
            cb = pb if isinstance(pb, CompositeBloq) else pb.decompose_bloq()
            d = {int(bi.i): bi.bloq for bi in cb.bloq_instances}
            self._decomp[path] = d
        return d

    # ------------------------------------------------------------------ nodes
    def node(self, path):
        nd = self._nodes.get(path)
        if nd is not None:
            return nd
        e = self.events.get(path)
        if e is None:
            raise AdapterError(f'no graph event {path}')
        b = self.bloq_of(path)
        sig = b.signature
        if set(e['inputs_msb']) != {r.name for r in sig.lefts()} or set(e['outputs_msb']) != {r.name for r in sig.rights()}:
            raise AdapterError(f'{path}: the event registers differ from the signature')
        ins = register_words(sig.lefts(), e['inputs_msb'])
        outs = register_words(sig.rights(), e['outputs_msb'])
        tname = type(b).__name__
        kind = 'wires' if tname in WIRES else 'alloc' if tname == 'Allocate' else 'free' if tname == 'Free' else 'op'
        sp = self.spans.get(path)
        span, contiguous = (sp if sp is not None else (None, True))
        try:
            text = repr(b)
        except Exception:
            text = tname
        key = tname + ':' + sha256(text.encode()).hexdigest()[:16]
        nd = Node(path=path, type=tname, kind=kind, ins=ins, outs=outs, span=span, key=key, obj=b,
                  depth=path.count('/'), info=dict(contiguous=contiguous))
        self._nodes[path] = nd
        return nd

    def children(self, node):
        return [self.node(p) for p in self.kids.get(node.path, [])]

    def contract(self, node):
        if self.registry is None:
            from . import registry as R
            self.registry = R
        return self.registry.contract_for(node.obj)

    def release(self, node):
        """forget the cached children of a proved subtree (memory)"""
        pre = node.path + '/'
        for p in [p for p in self._nodes if p.startswith(pre)]:
            del self._nodes[p]
        for p in [p for p in self._decomp if p == node.path or p.startswith(pre)]:
            del self._decomp[p]


def lower(bloq, *, registry=None, max_leaves=8_000_000, max_nodes=30_000_000, basis_states=True, log=None):
    """Lower ``bloq`` with the strict bridge (as in 0.7) and return a ``Program`` with its Qualtran adapter."""
    from ..qualtran.bridge.qualtran_layout import lower_bloqgraph
    from ..qualtran.bridge.qroam_clean_adjoint import chained_atom_rule
    from ..qualtran.bridge.rules import translate
    chain = chained_atom_rule()
    if basis_states:
        from ..qualtran.bridge.basis_vectors import BasisVectorRules
        bv = BasisVectorRules()

        def rule(b, inputs, allocate, path):
            r = bv(b, inputs, allocate, path)
            return r if r is not None else chain(b, inputs, allocate, path)
    else:
        rule = chain
    low = lower_bloqgraph(bloq, extra_atom_rule=rule, max_leaves=max_leaves, max_nodes=max_nodes)
    tr = translate(low, low['inputs'], low['outputs'], light=True)
    if log:
        log(f'lowered: nq={tr["nq"]} ops={len(tr["ops"])} events={len(low["graph_events"])}')
    ad = QualtranAdapter(bloq, low, tr, registry=registry)
    meta = dict(nq=tr['nq'], ops=len(tr['ops']), ops_sha256=ops_digest(tr['ops']),
                bridge_scratch=tr['layout'].get('bridge_scratch'), source_nq=low['nq'],
                program_inputs={k: list(v) for k, v in low['inputs'].items()},
                program_outputs={k: list(v) for k, v in low['outputs'].items()},
                retired_wires=len(low.get('retired_wires', ())), events=len(low['graph_events']))
    ops = tr['ops']
    del low['graph_events'][:]
    low.clear()
    tr['trace'] = None
    gc.collect()
    root = ad.node(ad.root_path)
    P = Program(ops, tr['nq'], root, ad, meta=meta)
    # the root must own every operation, its interface must be the program interface
    if root.span is not None and root.span != (0, len(ops)):
        raise AdapterError('the root does not own all operations')
    return P
