# SPDX-License-Identifier: Apache-2.0
"""Qualtran adapter: the lowered program and the tree of calls of a Qualtran bloq, for ``qoptct.foreign.auto``.

The program is lowered exactly as in the 0.7 proofs: ``bridge.qualtran_layout.lower_bloqgraph`` (the real
``decompose_bloq`` connections, fresh wires, rules A1 controlled-And and C1 ``QROAMCleanAdjoint``; the basis-state
atoms ``ZeroState``/``OneState``/``ZeroEffect``/``OneEffect`` of ``bridge.basis_vectors``) and
``bridge.rules.translate``.  The hash of the program is ``montgomery._stages.ops_digest`` of 0.7.

The tree comes from the graph events of the lowering (path ``Root/i:Type/j:Type``, the physical wires of every register
of every call) and from the trace of the translation (which operations each call produced).  Only the calls of the
root are recorded globally; the calls inside a call come from lowering that call alone and aligning its operations
with the program (checked operation by operation, an injective renaming of wires).  A node's operations must form one
contiguous range (checked here and refused by the prover otherwise).  Registers are turned into the
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


def _spans(trace, keep=None):
    """path -> ((start, end), contiguous) for every prefix of every leaf path that produced operations
    (only the paths in ``keep`` when given)"""
    spans = {}
    stack = []

    def close(entry):
        pf, st, cnt, end = entry
        if keep is not None and pf not in keep:
            return
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
    """Nodes of a lowered program.  Only the calls of the root are recorded by the global lowering; the calls inside
    a call are obtained by lowering that call's bloq alone and aligning its operations with the global ones (the same
    operations up to an injective renaming of wires, checked operation by operation).  This keeps the memory of very
    large programs small and every derived wire is checked against the real program."""

    def __init__(self, bloq, root_event, top_events, top_spans, rule_factory, registry=None):
        self.root_bloq = bloq
        self.registry = registry
        self.rule_factory = rule_factory
        self.root_path = type(bloq).__name__
        self.events = {self.root_path: root_event}
        self.events.update(top_events)
        self.spans = dict(top_spans)
        self.kids = {}
        if top_events:
            self.kids[self.root_path] = sorted(top_events, key=lambda q: int(q.rsplit('/', 1)[1].split(':', 1)[0]))
        else:
            self.kids[self.root_path] = []
        self.expanded = {self.root_path}
        self._bloq = {self.root_path: bloq}
        self._decomp = {}
        self._nodes = {}
        self._prim = {}
        self._sub = {}
        self.ops = None
        self.nq = None
        self._synth = 0
        self.alignments = 0

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
        nd = Node(path=path, type=tname, kind=kind, ins=ins, outs=outs, span=span, key=bloq_key(b), obj=b,
                  depth=path.count('/'), info=dict(contiguous=contiguous))
        self._nodes[path] = nd
        return nd

    def children(self, node):
        if node.path not in self.expanded:
            self._expand(node)
        return [self.node(p) for p in self.kids.get(node.path, [])]

    def contract(self, node):
        if self.registry is None:
            from . import registry as R
            self.registry = R
        return self.registry.contract_for(node.obj)

    # ------------------------------------------------------------------ calls inside a call
    def _expand(self, node):
        if node.kind != 'op':
            self.expanded.add(node.path)
            self.kids[node.path] = []
            return
        from qualtran import CompositeBloq
        b = node.obj
        try:
            b.decompose_bloq() if not isinstance(b, CompositeBloq) else None
        except Exception:
            self.expanded.add(node.path)            # an atom: no children
            self.kids[node.path] = []
            return
        key = node.key
        sub = self._sub.get(key)
        if sub is None:
            sub = _lower_alone(b, self.rule_factory)
            if len(self._sub) > 64:
                self._sub.clear()
            self._sub[key] = sub
        a, z = node.span if node.span is not None else (0, 0)
        if len(sub['ops']) != z - a:
            raise AdapterError(f'{node.path}: lowered alone it has {len(sub["ops"])} operations, in the program {z - a}')
        m, inv = {}, {}

        def bind(l, g):
            if m.get(l, g) != g or inv.get(g, l) != l:
                raise AdapterError(f'{node.path}: the operations differ from the call lowered alone (wire {l} -> {g})')
            m[l] = g
            inv[g] = l
        for side, words in (('ins', node.ins), ('outs', node.outs)):
            loc = sub[side]
            if set(loc) != set(words):
                raise AdapterError(f'{node.path}: interface of the call lowered alone differs')
            for name, ws in loc.items():
                if len(ws) != len(words[name]):
                    raise AdapterError(f'{node.path}: register {name} has another width when lowered alone')
                for l, g in zip(ws, words[name]):
                    bind(l, g)
        gops = self.ops
        for i, op in enumerate(sub['ops']):
            g = gops[a + i]
            if op[0] != g[0] or len(op) != len(g):
                raise AdapterError(f'{node.path}: operation {a + i} differs from the call lowered alone')
            if op[0] == 'erase':
                bind(op[1], g[1])
                if len(op[2]) != len(g[2]):
                    raise AdapterError(f'{node.path}: erase correction differs')
                for c1, c2 in zip(op[2], g[2]):
                    if c1[0] != c2[0] or len(c1) != len(c2):
                        raise AdapterError(f'{node.path}: erase correction differs')
                    for l, gg in zip(c1[1:], c2[1:]):
                        bind(l, gg)
            else:
                for l, gg in zip(op[1:], g[1:]):
                    bind(l, gg)

        def gw(l):
            if l not in m:
                # a wire of the call that no operation touches and that is not on its interface: it carries 0
                # wherever it goes; it gets a fresh name outside the program
                self._synth += 1
                m[l] = self.nq + self._synth
            return m[l]
        rootp = sub['root']
        kids = []
        for lp, e in sub['events'].items():
            gp = node.path + lp[len(rootp):]
            self.events[gp] = dict(path=gp, bloq_type=e['bloq_type'],
                                   inputs_msb={k: [gw(q) for q in v] for k, v in e['inputs_msb'].items()},
                                   outputs_msb={k: [gw(q) for q in v] for k, v in e['outputs_msb'].items()})
            if lp in sub['spans']:
                (la, lb), ok = sub['spans'][lp]
                self.spans[gp] = ((la + a, lb + a), ok)
            kids.append(gp)
        self.kids[node.path] = sorted(kids, key=lambda q: int(q.rsplit('/', 1)[1].split(':', 1)[0]))
        self.expanded.add(node.path)
        self.alignments += 1

    def release(self, node):
        """forget the cached children of a proved subtree (memory)"""
        pre = node.path + '/'
        for d in (self._nodes, self.events, self.spans):
            for p in [p for p in d if p.startswith(pre)]:
                del d[p]
        for p in [p for p in self.kids if p == node.path or p.startswith(pre)]:
            del self.kids[p]
        self.expanded = {p for p in self.expanded if not (p == node.path or p.startswith(pre))}


def bloq_key(b):
    try:
        text = repr(b)
    except Exception:
        text = type(b).__name__
    return type(b).__name__ + ':' + sha256(text.encode()).hexdigest()[:16]


def _rule_factory(basis_states=True):
    """the atom rules of the strict bridge, a fresh book per lowering: TwoBitCSwap as CX.CCX.CX (the opt-in rule
    ``bridge.fredkin_rule``, used by the 0.7 Kaliski client), basis states, controlled And (A1), QROAMCleanAdjoint (C1)"""
    def make():
        from ..qualtran.bridge.qroam_clean_adjoint import chained_atom_rule
        from ..qualtran.bridge.fredkin_rule import FredkinRules
        chain = chained_atom_rule()
        fredkin = FredkinRules()
        bv = None
        if basis_states:
            from ..qualtran.bridge.basis_vectors import BasisVectorRules
            bv = BasisVectorRules()

        def rule(b, inputs, allocate, path):
            r = fredkin(b, inputs, allocate, path)
            if r is None and bv is not None:
                r = bv(b, inputs, allocate, path)
            return r if r is not None else chain(b, inputs, allocate, path)
        return rule
    return make


def _lower_raw(bloq, rule_factory, max_leaves=8_000_000, max_nodes=30_000_000):
    from ..qualtran.bridge.qualtran_layout import lower_bloqgraph
    from ..qualtran.bridge.rules import translate
    import os
    low = lower_bloqgraph(bloq, extra_atom_rule=rule_factory(), max_leaves=max_leaves, max_nodes=max_nodes,
                          record_events='top')
    old_skip = os.environ.get('BRIDGE_SKIP_QASM')
    os.environ['BRIDGE_SKIP_QASM'] = '1'       # no QASM text of the whole program (memory), as the 0.7 clients
    try:
        tr = translate(low, low['inputs'], low['outputs'], light=True)
    finally:
        if old_skip is None:
            os.environ.pop('BRIDGE_SKIP_QASM', None)
        else:
            os.environ['BRIDGE_SKIP_QASM'] = old_skip
    return low, tr


def _root_flat(regs, words):
    """event-style MSB-first flat lists of the root registers from the LSB-first element words of the lowering"""
    from ..qualtran.bridge.qualtran_layout import element_name
    out = {}
    for reg in regs:
        flat = []
        for idx in reg.all_idxs():
            flat.extend(reversed(words[element_name(reg.name, idx)]))
        out[reg.name] = flat
    return out


def _split_lowering(bloq, low, tr):
    root = type(bloq).__name__
    events = {e['path']: e for e in low['graph_events'] if e['path'].count('/') == 1}
    if len(events) != sum(1 for e in low['graph_events'] if e['path'].count('/') == 1):
        raise AdapterError('repeated event paths')
    spans = _spans(tr['trace'], keep=set(events) | {root})
    sig = bloq.signature
    root_event = dict(path=root, bloq_type=root, inputs_msb=_root_flat(sig.lefts(), low['inputs']),
                      outputs_msb=_root_flat(sig.rights(), low['outputs']))
    return root, root_event, events, spans


def _lower_alone(bloq, rule_factory):
    """a call lowered alone (local wire numbering): ops, interface words, its top-level events and spans"""
    low, tr = _lower_raw(bloq, rule_factory)
    root, root_event, events, spans = _split_lowering(bloq, low, tr)
    sig = bloq.signature
    out = dict(root=root, ops=tr['ops'], ins=register_words(sig.lefts(), root_event['inputs_msb']),
               outs=register_words(sig.rights(), root_event['outputs_msb']), events=events, spans=spans)
    low.clear()
    tr['trace'] = None
    return out


def lower(bloq, *, registry=None, max_leaves=8_000_000, max_nodes=30_000_000, basis_states=True, log=None):
    """Lower ``bloq`` with the strict bridge (as in 0.7) and return a ``Program`` with its Qualtran adapter."""
    rf = _rule_factory(basis_states)
    low, tr = _lower_raw(bloq, rf, max_leaves=max_leaves, max_nodes=max_nodes)
    if log:
        log(f'lowered: nq={tr["nq"]} ops={len(tr["ops"])}')
    meta = dict(nq=tr['nq'], ops=len(tr['ops']), ops_sha256=ops_digest(tr['ops']),
                bridge_scratch=tr['layout'].get('bridge_scratch'), source_nq=low['nq'],
                program_inputs={k: list(v) for k, v in low['inputs'].items()},
                program_outputs={k: list(v) for k, v in low['outputs'].items()},
                retired_wires=len(low.get('retired_wires', ())))
    root, root_event, events, spans = _split_lowering(bloq, low, tr)
    ops = tr['ops']
    low.clear()
    tr['trace'] = None
    gc.collect()
    ad = QualtranAdapter(bloq, root_event, events, spans, rf, registry=registry)
    ad.ops = ops
    ad.nq = tr['nq']
    if ops:
        ad.spans[root] = ((0, len(ops)), True)
    rnode = ad.node(root)
    P = Program(ops, tr['nq'], rnode, ad, meta=meta)
    return P
