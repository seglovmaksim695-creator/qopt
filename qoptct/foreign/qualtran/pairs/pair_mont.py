"""Pair DirtyOutOfPlaceMontgomeryModMul (uncompute=False) + (uncompute=True), window w: assembly and expansion by the bridge
(rule QROAMCleanAdjoint of the bridge; Qualtran variant: the FIXED one, kind 'c')."""
import sys, json, time
from qoptct.foreign.qualtran.bridge.qroam_clean_adjoint import chained_atom_rule
from qoptct.foreign.qualtran.bridge.qualtran_layout import lower_bloqgraph
from qoptct.foreign.qualtran.bridge.rules import translate
from qoptct.verify.compose.lift import lift
from .pair_kal import ops_list


def pair_bloq(n, w, p, bwd_uncompute=True):
    from qualtran import BloqBuilder, Register, QMontgomeryUInt, QAny, QBit
    from qualtran.bloqs.mod_arithmetic.mod_multiplication import DirtyOutOfPlaceMontgomeryModMul as D
    bb = BloqBuilder()
    x = bb.add_register(Register('x', QMontgomeryUInt(n, p)))
    y = bb.add_register(Register('y', QMontgomeryUInt(n, p)))
    f = D(bitsize=n, window_size=w, mod=p)
    x, y, t, q, r = bb.add(f, x=x, y=y)
    x, y = bb.add(D(bitsize=n, window_size=w, mod=p, uncompute=bwd_uncompute), x=x, y=y, target=t, qrom_indices=q, reduced=r)
    return bb.finalize(x=x, y=y)


def lower_pair(n, w, p, bloq=None):
    bloq = bloq or pair_bloq(n, w, p)
    low = lower_bloqgraph(bloq, extra_atom_rule=chained_atom_rule(), max_leaves=50_000_000, max_nodes=50_000_000)
    tr = translate(low, low['inputs'], low['outputs'])
    return low, tr


if __name__ == '__main__':
    n, w, p = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    t = time.time()
    low, tr = lower_pair(n, w, p)
    print('nq', tr['nq'], 'ops', len(tr['ops']), 'inputs', {k: v[:3] for k, v in low['inputs'].items()}, 'sec', round(time.time() - t, 1))
    try:
        L = lift(tr['qasm'])
        print('lift ok', L.n, L.nq, L.nmeas)
    except Exception as e:  # noqa: BLE001
        print('lift FAIL', repr(e)[:300])
    import collections
    seg = collections.OrderedDict()
    for row in tr['trace']:
        k = '/'.join(row['leaf']['path'].split('/')[:3])
        seg.setdefault(k, [row['logical_start'], row['logical_end']])[1] = row['logical_end']
    for k, v in list(seg.items())[:40]:
        print(k, v[1] - v[0])
    print(len(seg))
