"""Pair KaliskiModInverse(uncompute=False) + (uncompute=True) as one composite block; full expansion by the bridge.
Mirror check: the second half of the program == the pi-renamed mirror of the first.
  mirror: reverse order; x->x; cx->cx; AND(a,b->t) -> ERASE(a,b;t); ERASE(a,b;t,l) -> AND(a,b->t)[; CX(l->t)].
  pi is built during matching (injective both ways); required: pi = id on the wires that carry data at the joint
  (output x and junk of the forward block), and pi(input wires of x) = output wires of x of the pair."""
import sys, json, time
from qoptct.verify.compose.lift import lift, K_X, K_CX, K_AND, K_ERASE


def pair_bloq(n, p, bwd_uncompute=True):
    from qualtran import BloqBuilder, Register, QMontgomeryUInt
    from qualtran.bloqs.mod_arithmetic.mod_division import KaliskiModInverse
    bb = BloqBuilder()
    x = bb.add_register(Register('x', QMontgomeryUInt(n, p)))
    x, junk = bb.add(KaliskiModInverse(n, p), x=x)
    x = bb.add(KaliskiModInverse(n, p, uncompute=bwd_uncompute), x=x, junk=junk)
    return bb.finalize(x=x)


def lift_pair(n, p, bloq=None):
    """Full expansion of the pair by the bridge (the one lowering routine of the Kaliski client): ``(lifted, low, prog)``."""
    from ..kaliski.lower import lower_full
    L, low, prog, _binding = lower_full(bloq or pair_bloq(n, p))
    return L, low, prog


def ops_list(L):
    a = L.ops
    return [tuple(a[5 * i:5 * i + 5]) for i in range(L.n)]


def wires(op):
    k, a, b, t, l = op
    if k == K_X:
        return (a,)
    if k == K_CX:
        return (a, b)
    return (a, b, t) if l < 0 else (a, b, t, l)


def mirror(op):
    k, a, b, t, l = op
    if k == K_X:
        return [(K_X, a, 0, 0, -1)]
    if k == K_CX:
        return [(K_CX, a, b, 0, -1)]
    if k == K_AND:
        return [(K_ERASE, a, b, t, -1)]
    r = [(K_AND, a, b, t, -1)]
    if l >= 0:
        r.append((K_CX, l, t, 0, -1))
    return r


def mirror_prog(ops):
    out = []
    for op in reversed(ops):
        out.extend(mirror(op))
    return out


def match(fwd, bwd):
    """match bwd with the mirror of fwd -> (ok, pi, msg, nswap). The order of the wires must agree (the controls a,b of AND/ERASE
    are symmetric by the meaning of the gate; a<->b substitution is allowed and counted in nswap)."""
    mir = mirror_prog(fwd)
    if len(mir) != len(bwd):
        return False, None, 'lengths: mirror %d, backward %d' % (len(mir), len(bwd)), 0
    pi, inv = {}, {}
    nswap = 0

    def fits(wm, wb):
        return all(pi.get(u, v) == v and inv.get(v, u) == u for u, v in zip(wm, wb))
    for i, (m, b) in enumerate(zip(mir, bwd)):
        if m[0] != b[0]:
            return False, None, 'kind of gate, position %d: %r versus %r' % (i, m, b), nswap
        wm, wb = wires(m), wires(b)
        if m[0] in (K_AND, K_ERASE) and not fits(wm, wb):
            wm2 = (wm[1], wm[0]) + wm[2:]
            if fits(wm2, wb):
                wm = wm2
                nswap += 1
        if not fits(wm, wb):
            return False, None, 'renumbering conflict, position %d: %r versus %r' % (i, m, b), nswap
        for u, v in zip(wm, wb):
            pi[u] = v
            inv[v] = u
    return True, pi, 'ok', nswap


def split_pair(prog):
    tr = {}
    for row in prog['trace']:
        k = row['leaf']['path'].split('/')[1]
        tr.setdefault(k, [row['logical_start'], row['logical_end']])[1] = row['logical_end']
    keys = list(tr)
    assert len(keys) == 2 and tr[keys[0]][1] == tr[keys[1]][0], keys
    return tr[keys[0]][1]


def join_wires(low):
    ev = {e['path']: e for e in low['graph_events']}
    e0 = [e for pth, e in ev.items() if pth.count('/') == 1 and pth.split('/')[1].split(':')[0] == '0'][0]
    out = set()
    for v in e0['outputs_msb'].values():
        out |= set(int(q) for q in (v if isinstance(v, (list, tuple)) else [v]))
    return sorted(out)


def check(n, p, bwd_uncompute=True):
    L, low, prog = lift_pair(n, p, pair_bloq(n, p, bwd_uncompute))
    ops = ops_list(L)
    F = split_pair(prog)
    fwd, bwd = ops[:F], ops[F:]
    C = join_wires(low) + list(prog['layout']['bridge_scratch'])
    from . import mirror as M
    res = M.match(fwd, bwd, C, scratch=list(prog['layout']['bridge_scratch']))
    ok, pi, msg = res['ok'], res['pi'], res['msg']
    r = dict(n=n, p=p, nops=len(ops), fwd=len(fwd), bwd=len(bwd), nq=L.nq, mirror_ok=ok, msg=msg, **res['stats'])
    if ok:
        r['join_wires'] = len(C)
        r['pi_id_on_join'] = all(pi.get(q, q) == q for q in C)
        r['pi_xin_is_xout'] = [pi.get(q) for q in low['inputs']['x']] == list(low['outputs']['x'])
        r['pi_injective'] = len(set(pi.values())) == len(pi)
    return r


if __name__ == '__main__':
    n, p = int(sys.argv[1]), int(sys.argv[2])
    t = time.time()
    r = check(n, p)
    r['sec'] = round(time.time() - t, 1)
    print(json.dumps(r))
