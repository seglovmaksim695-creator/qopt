"""lower -- the lifted program of the Qualtran circuit KaliskiModInverse, composed from the real pieces.

The full expansion by the bridge at N=256 is too large (2N iterations).  Therefore:
  * the top graph of the public bloq is expanded by the bridge as usual (lower_bloqgraph, rules A1/Fredkin/basis states),
    but the node _KaliskiIteration becomes a macro leaf: it allocates as many new wires, in the same order, as its full expansion
    would, and hands the outputs out by the real wire map of the piece;
  * one iteration is expanded by the bridge separately (the same Qualtran object as in all 2N nodes), translated by
    translate(), printed to QASM and lifted by the unchanged qoptct lift;
  * program = the lift of the top level, with the lifted piece inserted in place of every macro leaf under an injective renaming
    of the wires (inputs -> wires of the macro leaf, own -> the allocated block, scratch -> common).
check_against_full(n, p) checks for small N that the composed program is BYTE-EQUAL to the lift of the full expansion.
It proves nothing: it is a translation of the circuit (trusted part: bridge + composition; the composition is checked by comparison).
"""
from array import array
from bisect import bisect_right
import hashlib

import numpy as np

from ..bridge.from_qualtran import lower_graph
from ..bridge.fredkin_rule import FredkinRules
from ..bridge.qualtran_layout import AtomExpansion
from ..bridge.rules import translate
from qoptct.verify.compose.lift import lift, K_X, K_CX, K_AND, K_ERASE

BIG = 10 ** 9


class LowerError(Exception):
    pass


def need(c, msg):
    if not c:
        raise LowerError(msg)


def _top_child(path):
    """'_KaliskiIteration/2:_KaliskiIterationStep1/1:MultiAnd/...' -> '2:_KaliskiIterationStep1'"""
    parts = path.split('/')
    return parts[1] if len(parts) > 1 else None


class Piece:
    """one iteration: expansion, translation, lift; the map of the steps in local numbering."""

    def __init__(self, n, p):
        from qualtran.bloqs.mod_arithmetic.mod_division import _KaliskiIteration
        self.n, self.p = n, p
        self.bloq = _KaliskiIteration(n, p)
        book = FredkinRules()
        self.low = low = lower_graph(self.bloq, max_leaves=BIG, max_nodes=BIG, extra_atom_rule=book.chain_with_basis_and_a1)
        self.prog = prog = translate(low, low['inputs'], low['outputs'])
        self.binding = book.bind(low, prog)
        self.qasm_sha256 = hashlib.sha256(prog['qasm'].encode()).hexdigest()
        self.L = L = lift(prog['qasm'])
        need(L.n == len(prog['ops']) and L.nq == prog['nq'], 'lift differs from translate')
        self.nsrc = low['nq']                       # source wires; scratch (if any) = nsrc
        self.scratch = prog['layout']['bridge_scratch']
        need(self.scratch in ([], [self.nsrc]), 'unexpected scratch layout')
        self.nin = sum(len(e['wires_msb']) for d in low['registers']['inputs'] for e in d['elements'])
        need(sorted(q for d in low['registers']['inputs'] for e in d['elements'] for q in e['wires_msb'])
             == list(range(self.nin)), 'inputs are not wires 0..nin-1')
        # op ranges of the six steps (contiguous, in order)
        rng = []
        for row in prog['trace']:
            ch = _top_child(row['leaf']['path'])
            if rng and rng[-1][0] == ch:
                need(rng[-1][2] == row['logical_start'], 'non-contiguous step')
                rng[-1][2] = row['logical_end']
            else:
                need(all(r[0] != ch for r in rng), 'step appears twice: %s' % ch)
                rng.append([ch, row['logical_start'], row['logical_end']])
        names = [r[0].split(':')[1] for r in rng]
        need(names == ['_KaliskiIterationStep%d' % j for j in range(1, 7)], 'unexpected step order %r' % names)
        need(rng[0][1] == 0 and rng[-1][2] == L.n, 'steps do not cover the piece')
        ev = {e['path']: e for e in low['graph_events']}
        self.steps = []
        for ch, lo, hi in rng:
            e = ev['_KaliskiIteration/' + ch]
            self.steps.append(dict(name=ch.split(':')[1], lo=lo, hi=hi, ins=e['inputs_msb'], outs=e['outputs_msb']))
        top = ev['_KaliskiIteration']
        self.ins, self.outs = top['inputs_msb'], top['outputs_msb']
        self.allocs = [e['outputs_msb']['reg'][0] for e in low['graph_events']
                       if e['bloq_type'] == 'Allocate' and e['path'].count('/') == 1]
        self.frees = [e['inputs_msb']['reg'][0] for e in low['graph_events']
                      if e['bloq_type'] == 'Free' and e['path'].count('/') == 1]
        need(len(self.allocs) == 2 and self.frees == self.allocs, 'a, b: unexpected allocation')


class Composed:
    """the composed program: nq, ops (lazy view), n; plus the map of the segments in global numbering."""

    def __init__(self, n, p, piece=None, keep_full_leaves=False):
        from qualtran import Register, QAny
        from qualtran.bloqs.mod_arithmetic.mod_division import KaliskiModInverse, _KaliskiIteration
        self.n, self.p = n, p
        self.piece = P = piece or Piece(n, p)
        nanc = P.nsrc - P.nin
        macros = []
        book = FredkinRules()
        reg_in = P.low['registers']['inputs']

        def rule(bloq, inputs, allocate, path):
            if type(bloq) is _KaliskiIteration:
                need(bloq == P.bloq, 'iteration bloq differs from the lowered piece')
                anc = allocate([Register('anc', QAny(nanc))], path + '/macro')['anc'].ravel() if nanc else []
                g = {}
                for d in reg_in:
                    loc = d['elements'][0]['wires_msb']
                    glob = [int(q) for q in inputs[d['name']].flat]
                    need(len(loc) == len(glob), 'macro: register width')
                    g.update(zip(loc, glob))
                for j in range(nanc):
                    g[P.nin + j] = int(anc[j])
                need(len(g) == P.nsrc and len(set(g.values())) == P.nsrc, 'macro: map not injective')
                outs = {}
                for d in P.low['registers']['outputs']:
                    loc = d['elements'][0]['wires_msb']
                    outs[d['name']] = np.asarray([g[w] for w in loc], dtype=object).reshape(inputs[d['name']].shape)
                macros.append(dict(path=path, g=g))
                leaf = dict(name='macro', qubits=sorted(g.values()), path=path, macro=len(macros) - 1)
                return AtomExpansion(outs, [leaf], tuple(g[w] for w in P.low['retired_wires']))
            return book.chain_with_basis_and_a1(bloq, inputs, allocate, path)
        self.bloq = KaliskiModInverse(n, p)
        self.top = top = lower_graph(self.bloq, max_leaves=BIG, max_nodes=BIG, extra_atom_rule=rule)
        need(len(macros) == 2 * n, 'expected 2n iteration nodes')
        # translate the top without the macro leaves; remember where they go
        leaves, mpos = [], []
        for lf in top['leaves']:
            if lf['name'] == 'macro':
                mpos.append((len(leaves), lf['macro']))
            else:
                leaves.append(lf)
        need([m for _, m in mpos] == list(range(2 * n)), 'macro order')
        # wires that may be non-zero for the static checks of translate(): the input and everything a macro touches
        freew = sorted({q for w in top['inputs'].values() for q in w} |
                       {q for lf in top['leaves'] if lf['name'] == 'macro' for q in lf['qubits']})
        tprog = translate(dict(nq=top['nq'], leaves=leaves), {'all': freew}, {'all': freew})
        TL = lift(tprog['qasm'])
        need(TL.n == len(tprog['ops']), 'top lift')
        starts = [row['logical_start'] for row in tprog['trace']] + [len(tprog['ops'])]
        self.nsrc = top['nq']
        uses_scr = bool(tprog['layout']['bridge_scratch']) or bool(P.scratch)
        need(tprog['layout']['bridge_scratch'] in ([], [self.nsrc]), 'top scratch')
        self.nq = self.nsrc + (1 if uses_scr else 0)
        self.scratch = self.nsrc if uses_scr else None
        # wire maps of the pieces (arrays local -> global)
        self.gmaps = []
        for m in macros:
            g = np.zeros(P.L.nq, dtype=np.int64)
            for a, b in m['g'].items():
                g[a] = b
            if P.scratch:
                g[P.nsrc] = self.scratch
            self.gmaps.append(g)
        # parts: (global start op, nops, source array, source start op, map or None)
        parts, cur, prev = [], 0, 0
        top_cut = [starts[i] for i, _ in mpos]
        self.iter_start = []
        for k, cut in enumerate(top_cut):
            if cut > prev:
                parts.append((cur, cut - prev, TL.ops, prev, None))
                cur += cut - prev
            prev = cut
            self.iter_start.append(cur)
            parts.append((cur, P.L.n, P.L.ops, 0, k))
            cur += P.L.n
        if TL.n > prev:
            parts.append((cur, TL.n - prev, TL.ops, prev, None))
            cur += TL.n - prev
        self.parts, self.nops = parts, cur
        self.ops = OpsView(self)
        self.n_ops = cur
        self.TL, self.tprog, self.top_starts, self.mpos = TL, tprog, starts, mpos
        self.inputs, self.outputs = top['inputs'], top['outputs']
        self.top_events = top['graph_events']
        self.macros = macros
        self.nmeas = TL.nmeas + 2 * n * P.L.nmeas

    # a duck-typed Lifted for the walker / simulate_flat
    @property
    def lifted(self):
        return _LiftedView(self)

    def materialize(self):
        out = array('i')
        for k in range(len(self.parts)):
            out.extend(self.ops.part_array(k))
        return out

    def top_range(self, top_index):
        """op range (global, [lo, hi)) of the leaves of one top-level node 'KaliskiModInverse/<i>:...'; None if no ops."""
        lo = hi = None
        for j, row in enumerate(self.tprog['trace']):
            pth = row['leaf']['path'].split('/')
            if len(pth) > 1 and pth[1].split(':')[0] == str(top_index) and row['logical_end'] > row['logical_start']:
                a = self.glob_top_op(row['logical_start'])
                b = self.glob_top_op(row['logical_end'] - 1) + 1
                need(b - a == row['logical_end'] - row['logical_start'], 'a macro inside a top leaf')
                lo = a if lo is None else lo
                need(hi is None or hi == a, 'top node not contiguous')
                hi = b
        return None if lo is None else (lo, hi)

    def glob_top_op(self, i):
        """global index of op i of the top program (a macro inserted at position c precedes top op c)."""
        nb = sum(1 for (pos, _) in self.mpos if self.top_starts[pos] <= i)
        return i + nb * self.piece.L.n


class _LiftedView:
    def __init__(self, C):
        self.nq, self.ops, self.n, self.nmeas = C.nq, C.ops, C.nops, C.nmeas
        self.style = 1


class OpsView:
    """ops[i] of the composed program (5 ints per op), materialising one renamed part at a time."""

    def __init__(self, C):
        self.C = C
        self.starts = [5 * pt[0] for pt in C.parts]
        self.cache = (0, 0, None)

    def __len__(self):
        return 5 * self.C.nops

    def part_array(self, j):
        g0, nops, src, s0, k = self.C.parts[j]
        if k is None:
            return src[5 * s0:5 * (s0 + nops)]
        a = np.frombuffer(src, dtype=np.int32)[5 * s0:5 * (s0 + nops)].reshape(nops, 5).astype(np.int64)
        g = self.C.gmaps[k]
        kind = a[:, 0]
        b = a.copy()
        b[:, 1] = g[a[:, 1]]
        m2 = kind != K_X
        b[m2, 2] = g[a[m2, 2]]
        m3 = kind >= K_AND
        b[m3, 3] = g[a[m3, 3]]
        m4 = m3 & (a[:, 4] >= 0)
        b[m4, 4] = g[a[m4, 4]]
        return array('i', b.astype(np.int32).ravel().tobytes())

    def __getitem__(self, i):
        lo, hi, arr = self.cache
        if lo <= i < hi:
            return arr[i - lo]
        if not 0 <= i < 5 * self.C.nops:
            raise IndexError(i)
        j = bisect_right(self.starts, i) - 1
        arr = self.part_array(j)
        lo = self.starts[j]
        self.cache = (lo, lo + len(arr), arr)
        return arr[i - lo]


def full_lift(n, p):
    """full expansion by the bridge (small N only): (lifted, lowered)."""
    from qualtran.bloqs.mod_arithmetic.mod_division import KaliskiModInverse
    book = FredkinRules()
    low = lower_graph(KaliskiModInverse(n, p), max_leaves=BIG, max_nodes=BIG, extra_atom_rule=book.chain_with_basis_and_a1)
    prog = translate(low, low['inputs'], low['outputs'])
    book.bind(low, prog)
    return lift(prog['qasm']), low


def check_against_full(n, p):
    C = Composed(n, p)
    FL, flow = full_lift(n, p)
    mine = C.materialize()
    ok = (FL.nq == C.nq and FL.n == C.nops and FL.ops == mine and flow['inputs'] == C.inputs
          and flow['outputs'] == C.outputs)
    return dict(n=n, p=p, identical=ok, nq=(FL.nq, C.nq), nops=(FL.n, C.nops),
                sha_full=hashlib.sha256(FL.ops.tobytes()).hexdigest(), sha_composed=hashlib.sha256(mine.tobytes()).hexdigest())


if __name__ == '__main__':
    import sys, json, time
    for a in sys.argv[1:]:
        n, p = map(int, a.split(':'))
        t = time.time()
        r = check_against_full(n, p)
        r['sec'] = round(time.time() - t, 2)
        print(json.dumps(r), flush=True)
