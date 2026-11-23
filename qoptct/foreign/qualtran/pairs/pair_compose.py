"""Kaliski pair (forward + backward) in macro form: 2n forward iterations and 2n backward (Adjoint) iterations become macro
leaves; one forward and one backward iteration are expanded by the bridge separately.
Mirror check by parts (see check_pair):
  (1) piece: the backward iteration = the pi_p-mirror of the forward iteration (pi_p is given on the interface: forward output -> backward input);
  (2) assembly of the global pi from the pieces over the 2n macro nodes, check of compatibility and injectivity;
  (3) top level: tail of the forward ~ head of the backward, head of the forward ~ tail of the backward, with the global pi as the seed;
  (4) pi = id on the joint wires, pi(input x) = output x.
For small n the assembly is compared ELEMENT BY ELEMENT with the full expansion (materialize == lift_pair)."""
import sys, json, time, hashlib
import numpy as np
from qoptct.foreign.qualtran.bridge.from_qualtran import lower_graph
from qoptct.foreign.qualtran.bridge.fredkin_rule import FredkinRules
from qoptct.foreign.qualtran.bridge.rules import translate
from qoptct.verify.compose.lift import lift
from qoptct.foreign.qualtran.kaliski.lower import (BIG, LowerError as CErr, need, lower_full, piece_interface,
                                                   macro_leaf)
from . import mirror as M
from .pair_kal import pair_bloq, ops_list, join_wires, lift_pair


class PieceG:
    def __init__(self, n, p, adjoint):
        from qualtran.bloqs.mod_arithmetic.mod_division import _KaliskiIteration
        self.n, self.p, self.adjoint = n, p, adjoint
        b = _KaliskiIteration(n, p)
        self.bloq = b.adjoint() if adjoint else b
        self.L, self.low, self.prog, _binding = lower_full(self.bloq)
        self.ops = ops_list(self.L)
        self.nsrc, self.scratch, self.nin = piece_interface(self.low, self.prog)
        low, prog = self.low, self.prog
        self.ins = {d['name']: [q for e in d['elements'] for q in e['wires_msb']] for d in low['registers']['inputs']}
        self.outs = {d['name']: [q for e in d['elements'] for q in e['wires_msb']] for d in low['registers']['outputs']}
        self.sha = hashlib.sha256(prog['qasm'].encode()).hexdigest()


class PairComposed:
    def __init__(self, n, p, pieces=None):
        from qualtran import Adjoint
        from qualtran.bloqs.mod_arithmetic.mod_division import _KaliskiIteration
        self.n, self.p = n, p
        self.Pf, self.Pb = pieces or (PieceG(n, p, False), PieceG(n, p, True))
        itb = _KaliskiIteration(n, p)
        macros, book = [], FredkinRules()
        Pf, Pb = self.Pf, self.Pb

        def rule(bloq, inputs, allocate, path):
            P, kind = None, None
            if type(bloq) is _KaliskiIteration:
                P, kind = Pf, 'f'
                need(bloq == Pf.bloq, 'forward iteration differs')
            elif type(bloq) is Adjoint and bloq.subbloq == itb:
                P, kind = Pb, 'b'
                need(bloq == Pb.bloq, 'adjoint iteration differs')
            if P is None:
                return book.chain_with_basis_and_a1(bloq, inputs, allocate, path)
            return macro_leaf(P, inputs, allocate, path, macros, kind=kind)
        self.bloq = pair_bloq(n, p)
        self.top = top = lower_graph(self.bloq, max_leaves=BIG, max_nodes=BIG, extra_atom_rule=rule)
        need(len(macros) == 4 * n, 'expected 4n iteration nodes, got %d' % len(macros))
        leaves, mpos = [], []
        for lf in top['leaves']:
            if lf['name'] == 'macro':
                mpos.append((len(leaves), lf['macro']))
            else:
                leaves.append(lf)
        need([m for _, m in mpos] == list(range(4 * n)), 'macro order')
        freew = sorted({q for w in top['inputs'].values() for q in w} |
                       {q for lf in top['leaves'] if lf['name'] == 'macro' for q in lf['qubits']})
        tprog = translate(dict(nq=top['nq'], leaves=leaves), {'all': freew}, {'all': freew})
        TL = lift(tprog['qasm'])
        need(TL.n == len(tprog['ops']), 'top lift')
        starts = [row['logical_start'] for row in tprog['trace']] + [len(tprog['ops'])]
        self.nsrc = top['nq']
        uses_scr = bool(tprog['layout']['bridge_scratch']) or bool(Pf.scratch) or bool(Pb.scratch)
        need(tprog['layout']['bridge_scratch'] in ([], [self.nsrc]), 'top scratch')
        self.nq = self.nsrc + (1 if uses_scr else 0)
        self.scratch = [self.nsrc] if uses_scr else []
        self.top_ops = ops_list(TL)
        self.top_cut = [starts[i] for i, _ in mpos]
        self.macros = macros
        # which half: the trace of the top level
        halves = [row['leaf']['path'].split('/')[1].split(':')[0] for row in tprog['trace']]
        need(set(halves) <= {'0', '1'}, 'halves')
        first_b = next((i for i, h in enumerate(halves) if h == '1'), len(halves))
        need(all(h == '0' for h in halves[:first_b]) and all(h == '1' for h in halves[first_b:]), 'halves not contiguous')
        self.F_end = starts[first_b]
        self.top_low, self.tprog = top, tprog
        for m in macros:
            m['half'] = 'f' if m['path'].split('/')[1].split(':')[0] == '0' else 'b'
            need((m['half'] == 'f') == (m['kind'] == 'f'), 'half/kind mismatch')
            # global scratch wire of the piece
        self.nmac = len(macros)

    def gfull(self, m):
        P = self.Pf if m['kind'] == 'f' else self.Pb
        g = np.zeros(P.nsrc + (1 if P.scratch else 0), dtype=np.int64)
        for a, b in m['g'].items():
            g[a] = b
        if P.scratch:
            g[P.nsrc] = self.scratch[0]
        return g

    def materialize(self):
        out, pos = [], 0
        for k, cut in enumerate(self.top_cut):
            out.extend(self.top_ops[pos:cut])
            pos = cut
            m = self.macros[k]
            P = self.Pf if m['kind'] == 'f' else self.Pb
            g = self.gfull(m)
            for (kk, a, b, t, l) in P.ops:
                if kk == 0:
                    out.append((0, int(g[a]), 0, 0, -1))
                elif kk == 1:
                    out.append((1, int(g[a]), int(g[b]), 0, -1))
                else:
                    out.append((kk, int(g[a]), int(g[b]), int(g[t]), int(g[l]) if l >= 0 else -1))
        out.extend(self.top_ops[pos:])
        return out


def piece_pi(Pf, Pb):
    """mirror of the piece: pi_p (forward local wire -> backward local wire), the interface is given."""
    seed = {}
    for r, wf in Pf.outs.items():
        wb = Pb.ins[r]
        need(len(wf) == len(wb), 'width of register %s' % r)
        for a, b in zip(wf, wb):
            seed[a] = b
    for r, wf in Pf.ins.items():
        wb = Pb.outs[r]
        need(len(wf) == len(wb), 'width of register %s' % r)
        for a, b in zip(wf, wb):
            need(seed.get(a, b) == b, 'interface: wire %d is matched to two' % a)
            seed[a] = b
    need(len(set(seed.values())) == len(seed), 'the interface is not injective')
    res = M.match(Pf.ops, Pb.ops, [], pi0=seed, scratch=list(Pf.scratch), scratch_b=list(Pb.scratch))
    return res, seed


def check_composed(C, Pb=None, macros=None, T=None, piece_res=None, quick=False):
    """Pb: (possibly corrupted) backward iteration; macros: (possibly corrupted) list of macro nodes; T: top operations.
    piece_res: ready result of matching the piece (for canaries that do not change the piece)."""
    n = C.n
    Pb = Pb or C.Pb
    macros = macros if macros is not None else C.macros
    T = T if T is not None else C.top_ops
    r = dict(n=n, p=C.p)
    res = piece_res or piece_pi(C.Pf, Pb)[0]
    res.pop('dbg', None)
    r['piece_mirror_ok'] = res['ok']
    r['piece_msg'] = res['msg']
    r['piece_stats'] = res['stats']
    if not res['ok']:
        r['mirror_ok'] = False
        return r
    pip = res['pi']
    fm = [m for m in macros if m['kind'] == 'f']
    bm = [m for m in macros if m['kind'] == 'b']
    pi = np.full(C.nq + 1, -1, dtype=np.int64)
    inv = np.full(C.nq + 1, -1, dtype=np.int64)
    npi = [0]

    def put(u, v, where):
        a, b = int(pi[u]), int(inv[v])
        if (a != -1 and a != v) or (b != -1 and b != u):
            raise CErr('renumbering conflict (%s): %d -> %d, was %r / inverse %r' % (where, u, v, a, b))
        if a == -1:
            npi[0] += 1
        pi[u] = v
        inv[v] = u
    try:
        nw = getattr(C, 'nwin', 2 * n)
        need(len(fm) == len(bm) == nw, 'number of macro nodes')
        for k in range(nw):
            mf, mb = fm[k], bm[nw - 1 - k]
            for w, v in pip.items():
                if w in mf['g'] and v in mb['g']:
                    put(mf['g'][w], mb['g'][v], 'macro %d' % k)
        for q in C.scratch:
            put(q, q, 'scratch')
        r['pi_macro_size'] = npi[0]
        cf = {C.top_cut[i] for i, m in enumerate(macros) if m['kind'] == 'f'}
        cb = {C.top_cut[i] for i, m in enumerate(macros) if m['kind'] == 'b'}
        need(len(cf) == 1 and len(cb) == 1, 'macro nodes are not consecutive')
        cf, cb = cf.pop(), cb.pop()
        Fe = C.F_end
        F_pre, F_post, B_pre, B_post = T[:cf], T[cf:Fe], T[Fe:cb], T[cb:]
        # operations before the forward macro nodes that do not touch their wires commute with them: move them to the start of F_post
        if all('tw' in m for m in fm):
            twf = set().union(*[m['tw'] for m in fm])
            k = len(F_pre)
            while k > 0 and not (set(M.wires(F_pre[k - 1])) & twf):
                k -= 1
            r['moved_across_macros'] = len(F_pre) - k
            F_pre, F_post = F_pre[:k], F_pre[k:] + F_post
        r['top_parts'] = dict(F_pre=len(F_pre), F_post=len(F_post), B_pre=len(B_pre), B_post=len(B_post))
        J = join_wires(C.top_low)
        for q in J:
            put(q, q, 'joint')
        scr = list(C.scratch)
        touched = {w for part in (F_pre, F_post, B_pre, B_post) for o in part for w in M.wires(o)}
        seed_top = {u: int(pi[u]) for u in touched if pi[u] != -1}
        ok_all, msgs = True, []
        for nm, fw, bw in (('F_post~B_pre', F_post, B_pre), ('F_pre~B_post', F_pre, B_post)):
            res2 = M.match(fw, bw, [], pi0=seed_top, scratch=scr, scratch_b=scr)
            res2.pop('dbg', None)
            ok_all &= res2['ok']
            msgs.append((nm, res2['ok'], res2['msg'], res2['stats']))
            if res2['ok']:
                for u, v in res2['pi'].items():
                    put(u, v, nm)
        r['top_match'] = msgs
        r['top_ok'] = ok_all
        if ok_all:
            r['pi_id_on_join'] = all(int(pi[q]) == q for q in J)
            r['pi_xin_is_xout'] = all([int(pi[q]) for q in C.top_low['inputs'][k]] == list(C.top_low['outputs'][k]) for k in C.top_low['inputs'])
            r['pi_injective'] = True  # injectivity is checked at every put (array inv)
        r['mirror_ok'] = bool(res['ok'] and ok_all and r.get('pi_id_on_join') and r.get('pi_xin_is_xout') and r.get('pi_injective'))
    except CErr as e:
        r['mirror_ok'] = False
        r['error'] = str(e)
    return r


def check_pair(n, p, compare_full=False):
    t0 = time.time()
    C = PairComposed(n, p)
    r0 = dict(piece_fwd_ops=len(C.Pf.ops), piece_bwd_ops=len(C.Pb.ops), top_ops=len(C.top_ops), macros=C.nmac,
              qasm_sha_fwd_piece=C.Pf.sha, qasm_sha_bwd_piece=C.Pb.sha)
    if compare_full:
        L, low, prog = lift_pair(n, p)
        r0['equals_full_lowering'] = (ops_list(L) == C.materialize()) and L.nq == C.nq
    r = check_composed(C)
    r.update(r0)
    r['sec'] = round(time.time() - t0, 1)
    return r


if __name__ == '__main__':
    n, p = int(sys.argv[1]), int(sys.argv[2])
    cmp_full = len(sys.argv) > 3 and sys.argv[3] == 'full'
    try:
        r = check_pair(n, p, compare_full=cmp_full)
    except CErr as e:
        r = dict(n=n, p=p, mirror_ok=False, error=str(e))
    print(json.dumps(r, ensure_ascii=False, default=str))
