# SPDX-License-Identifier: Apache-2.0
"""Block resynthesis post-pass on QASM text (b10).
 1) segmentation into CX/CZ/X/phase blocks (segment2);
 2) per block: (V0) exact equality (original path with recomputed Clifford weights, GraySynth restarts, windows),
    (V1) + absorption of an ancilla permutation, (V2) + freedom of fresh |+> wires;
 3) the cheapest CX+CZ variant is kept; phase-polynomial equality is checked by coefficients (+ brute force for n <= 14).
The whole circuit is verified by the caller."""
import re, sys, time, random
from . import blocks as B, big, gray, window, perm, plus, cliff, synth2


def relabel_raw(line, rl):
    return re.sub(r'q\[(\d+)\]', lambda m: 'q[%d]' % rl.get(int(m.group(1)), int(m.group(1))), line)


def synth_paths(n, O, terms, A, rng, tries=12, orig_ops=None, win_tries=0, deadline=None):
    """Best (ops, sol, total) among the original path (if A matches) and GraySynth restarts; then windows."""
    best = None
    cands = []
    if orig_ops is not None: cands.append(list(orig_ops))
    for mode in ('max', 'bal'):
        for r in range(tries):
            if deadline and time.time() > deadline: break
            o = gray.gray_path(n, O, A, rng, mode=mode, tailtries=2)
            cands.append(o)
    for ops in cands:
        vis, fin = big.path_rows(n, ops)
        if fin != list(A): continue
        sol = big.solve_clifford(n, terms, ops, iters=1500)
        if sol is None: continue
        tot = len(ops) + sol['ncz']
        if best is None or tot < best[2]: best = (ops, sol, tot)
    if best is None: return None
    if win_tries:
        ops, tot, sol = window.window_pass(n, best[0], O, terms, rng, tries=win_tries, deadline=deadline)
        best = (ops, sol, tot)
    return best


def block_variants(ws, gates, nd, fresh, rng, tries=12, win_tries=0, small_exact=True):
    """Best variant: dict(cost, gates, pi, label) in local wire indices of the block."""
    p = B.block_poly(ws, gates)
    n = p['n']; terms = p['terms']; A = list(p['rows']); off = p['off']
    cost0 = p['ncx']
    idl = {w: i for i, w in enumerate(ws)}
    orig_ops = [(idl[a[0]], idl[a[1]]) for nm, a in gates if nm == 'cx']
    O = [m for m, w in terms.items() if w & 1]
    best = dict(cost=cost0, gates=None, pi=None, label='orig')
    r0 = synth_paths(n, O, terms, A, rng, tries=tries, orig_ops=orig_ops, win_tries=win_tries)
    if r0 is not None and r0[2] < best['cost']:
        g = big.build_big(n, off, r0[0], r0[1])
        q = B.block_poly(list(range(n)), g)
        if B.poly_equiv(dict(p, rows=A), q) and B.brute_equal(dict(p), q) in (True, None):
            best = dict(cost=r0[2], gates=g, pi=None, label='V0 exact')
    # V1, V2
    anc = [i for i in range(n) if ws[i] >= nd]
    if anc:
        pi = perm.best_perm(A, 0 if False else 0) if False else None
        first_anc = min(anc)
        assert anc == list(range(first_anc, n))
        pi = perm.best_perm(A, first_anc, rng)
        Ap = perm.apply_perm_rows(A, pi, first_anc)
        offp = list(off)
        for w in range(n): offp[pi.get(w, w)] = off[w]
        for variant in ('V1', 'V2'):
            if variant == 'V1':
                Av, Ov, tv = Ap, O, terms
            else:
                P = [i for i in range(n) if ws[i] in fresh]
                if not P: continue
                S0 = plus.make_S0(n, P, Ap)
                Ov, Av = plus.transform(S0, O, Ap)
                tv = {}
                for m, w in terms.items():
                    mm = gray.coords(S0, [m])[0]
                    tv[mm] = (tv.get(mm, 0) + w) % 8
            r = synth_paths(n, Ov, tv, Av, rng, tries=tries, orig_ops=None, win_tries=win_tries)
            if r is None: continue
            ops, sol, tot = r
            if tot < best['cost']:
                g = big.build_big(n, offp, ops, sol)
                q = B.block_poly(list(range(n)), g)
                ok = (q['rows'] == list(Av) and cliff.coeffs(q['terms']) == cliff.coeffs(tv))
                if ok:
                    best = dict(cost=tot, gates=g, pi=pi, label=variant + (' n=%d' % n))
    return best


def process(text, nd, seed=1, tries=12, win_tries=0, minblock=1, verbose=False, small_exact_n=9, deadline=None):
    head, nq, items = B.parse_qasm(text)
    seg = B.segment2(items)
    rng = random.Random(seed)
    touched = {}
    rl = {}
    out = []
    n_before = B.count(seg)
    hist = {}
    for it in seg:
        if it[0] == 'raw':
            ws_ = [rl.get(w, w) if isinstance(w, int) else w for w in it[2]]
            line = relabel_raw(it[1], rl)
            out.append(('raw', line, ws_))
            for w in it[2]:
                if isinstance(w, int): hist.setdefault(w, []).append(it[1].split()[0])
            continue
        ws, gates = it[1], it[2]
        # |+> freedom only for ancillas (wires >= nd start in |0>) whose first gate is this H
        fresh = {w for w in ws if w >= nd and hist.get(w, None) == ['h']}
        for w in ws: hist.setdefault(w, []).append('blk')
        p = B.block_poly(ws, gates)
        if p['ncx'] < minblock or (deadline is not None and time.time() > deadline):
            out.append(('block', [rl.get(w, w) for w in ws], [(nm, [rl.get(w, w) for w in a]) for nm, a in gates]))
            continue
        best = None
        if p['n'] <= small_exact_n and not any(w >= nd for w in ws) and not any(nm == 'cz' for nm, _ in gates):
            res = synth2.exact2(p['n'], p['terms'], p['ncx'], node_limit=150000, timeout=20, A=p['rows'])
            if isinstance(res, tuple):
                ops, V, u = res
                g = synth2.build_gates2(p['n'], p['terms'], p['off'], ops, V, u)
                q = B.block_poly(list(range(p['n'])), g)
                if B.poly_equiv(dict(p), q):
                    best = dict(cost=len(ops), gates=g, pi=None, label='exact')
        if best is None:
            big_ = p['n'] > small_exact_n or any(w >= nd for w in ws)
            best = block_variants(ws, gates, nd, fresh, rng, tries=tries if big_ else 3, win_tries=win_tries if big_ else 0)
        if best['gates'] is None:
            out.append(('block', [rl.get(w, w) for w in ws], [(nm, [rl.get(w, w) for w in a]) for nm, a in gates]))
            continue
        if verbose and best['cost'] < p['ncx']:
            print('  block n=%d: %d -> %d  (%s) cx=%d cz=%d' % (p['n'], p['ncx'], best['cost'], best['label'], sum(1 for g in best['gates'] if g[0]=='cx'), sum(1 for g in best['gates'] if g[0]=='cz')), flush=True)
        oldphys = [rl.get(ws[i], ws[i]) for i in range(p['n'])]
        gl = [(nm, [oldphys[i] for i in a]) for nm, a in best['gates']]
        out.append(('block', list(oldphys), gl))
        if best['pi']:
            pi = best['pi']
            for i in range(p['n']):
                rl[ws[i]] = oldphys[pi.get(i, i)]
    return head, nq, out, n_before
