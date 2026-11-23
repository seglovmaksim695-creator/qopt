"""The stages of the Montgomery proof; each stage runs in its own process (see ``__init__``) so that the memory of the
lowering is released before the proofs start.

lower  - lower + translate the multiplier with the strict bridge, cut it into segments and sub-nodes, run every
         structural check (wire glue, live-wire bookkeeping, coverage and freshness of every window, window hashes,
         end-to-end simulation at small n), write the compact cache of the segments that are proved.
inv    - interval invariants of the window (BDD, or interval arithmetic where BDD does not fit in memory).
prove  - kprove (BDD) + exact partial-domain guard for the segments; the generic window by its sub-nodes.
"""
import json
import os
import pickle
import random
import sys
import time
from hashlib import sha256

from . import _core as C
from . import _split as SP


def _qualtran_info(kind):
    import qualtran
    from importlib.metadata import version
    from pathlib import Path
    info = dict(version=version('qualtran'), file=None, patch_kind=kind, patch_files_sha256=None)
    root = Path(qualtran.__file__).resolve().parents[1]
    info['file'] = 'qualtran/__init__.py in ' + ('patched copy' if kind else 'installed Qualtran')
    if kind:
        from ..patches import verify_copy, PatchError
        try:
            info['patch_files_sha256'] = verify_copy(kind, root)
        except PatchError as e:
            raise RuntimeError('Qualtran in this process is not the patched copy %r: %s' % (kind, e))
    else:
        from ..patches import _MANIFEST
        # the original is only acceptable if it really is the original 0.7.0 file (the sign of the final reduction)
        from hashlib import sha256 as _s
        tgt = root / 'qualtran' / 'bloqs' / 'mod_arithmetic' / 'mod_multiplication.py'
        info['mod_multiplication_sha256'] = _s(tgt.read_bytes()).hexdigest()
    return info


def window_ops(segs, info, base0, i):
    """operations of window i (window 0: without the ``X(reduced)`` prefix that stands before its span)"""
    if i == 0:
        return segs[0][2][info[0]['span'][0] - base0:]
    return segs[i][2]


class Sl:
    """slice by GLOBAL operation indices of a window whose operations start at global index ``off``"""

    def __init__(s, ops, off):
        s.ops, s.off = ops, off

    def __getitem__(s, sl):
        return s.ops[sl.start - s.off:sl.stop - s.off]


def glue_check(segs, low, n, nw, tw, rwire):
    glue = dict(ok=True, notes=[])
    for i in range(nw - 1):
        if segs[i][4]['T'] != segs[i + 1][3]['T']: glue['ok'] = False; glue['notes'].append(f'T wires mismatch {i}->{i+1}')
        if segs[i][3]['y'] != segs[i + 1][3]['y']: glue['ok'] = False; glue['notes'].append('y wires differ')
    lastT = segs[nw - 1][4]['T']
    if lastT[:n] != tw: glue['ok'] = False; glue['notes'].append('final target wires != last window T[:n]')
    flat_m = [q for i in reversed(range(nw)) for q in segs[i][4]['m']]
    if flat_m != low['outputs']['qrom_indices']: glue['ok'] = False; glue['notes'].append('qrom_indices wires != concat(m_i), window 0 most significant')
    if low['outputs']['target'] != tw: glue['ok'] = False; glue['notes'].append('output target wires differ')
    if low['outputs']['reduced'] != [rwire]: glue['ok'] = False; glue['notes'].append('reduced wire differs')
    if low['outputs']['x'] != low['inputs']['x'] or low['outputs']['y'] != low['inputs']['y']: glue['ok'] = False; glue['notes'].append('x/y passthrough')
    xs = [q for i in range(nw) for q in segs[i][3]['x']]
    if xs != low['inputs']['x']: glue['ok'] = False; glue['notes'].append('x windows != x register (LSB-first)')
    if segs[0][3]['y'] != low['inputs']['y']: glue['ok'] = False; glue['notes'].append('y register')
    return glue


def end_to_end(low, tr, n, w, p, rng, count=24):
    """simulate the whole lowered program on input pairs and compare with the exact recurrence (fixed contract)"""
    cases = [(0, 0), (0, p - 1), ((1 << n) - 1, p - 1), ((1 << n) - 1, 0), (1, 1)]
    red = [0, 0]
    tries = 0
    while min(red) < count // 2 and tries < 4000:
        tries += 1
        x, y = rng.randrange(1 << n), rng.randrange(p)
        r = C.reference(n, w, p, x, y)['reduced']
        if red[r] < count // 2:
            red[r] += 1; cases.append((x, y))
    bad = 0; first = None; n_red = 0
    for x, y in cases:
        res, ph, rest = C.run_basis(low, tr, n, x, y)
        e = C.reference(n, w, p, x, y)
        n_red += e['reduced']
        if any(ph) or rest or any(res[k] != e[k] for k in ('x', 'y', 'target', 'qrom_indices', 'reduced')):
            bad += 1
            if first is None:
                first = dict(x=x, y=y, got={k: res[k] for k in ('target', 'reduced')}, expected={k: e[k] for k in ('target', 'reduced')},
                             T_N_ge_p=bool(e['T'] >= p), T_N=e['T'], bad_phase=bool(any(ph)), bad_rest=bool(rest))
    return dict(cases=len(cases), cases_with_reduction=n_red, mismatches=bad, first_witness=first)


def ops_digest(ops):
    h = sha256()
    for o in ops:
        h.update(repr(o).encode())
    return h.hexdigest()


def stage_lower(n, w, p, kind, out_path, mutant=None, e2e_max_n=64, seed=1):
    t0 = time.time()
    qinfo = _qualtran_info(kind)
    P = C.Params(n, w, p)
    rng = random.Random(seed)
    B = C.build_segs(n, w, p)
    low, tr, top, order, wins, rwire, tw = (B[k] for k in ('low', 'tr', 'top', 'order', 'wins', 'rwire', 'tw'))
    segs = [tuple(s) for s in B['segs']]
    nw = len(wins)
    info = SP.children_info(low, tr, top, wins, order)
    base0 = top[order[0]][0]
    low_io = dict(inputs=low['inputs'], outputs=low['outputs'])
    meta = dict(nq=tr['nq'], ops=len(tr['ops']), ops_sha256=ops_digest(tr['ops']), lower_seconds=round(time.time() - t0, 1))
    e2e = None
    if mutant is None and n <= e2e_max_n:
        e2e = end_to_end(low, tr, n, w, p, rng)
    elif mutant is None:
        e2e = dict(skipped='n > %d (the end-to-end simulation of the whole program is only run at small n)' % e2e_max_n)
    del B, tr, low
    if mutant is not None:
        from . import mutants as MUT
        segs, info = MUT.apply(mutant, P, segs, info, base0, nw)
    checks = structural_checks(P, segs, info, base0, nw, low_io, tw, rwire)
    sm = [(s[0], s[1], s[3], s[4], s[5]) for s in segs]
    slim = dict(schema=1, n=n, w=w, p=p, kind=kind, mutant=mutant, qualtran=qinfo, meta=meta, low=low_io, rwire=rwire, tw=tw, nw=nw,
                seg_meta=sm, seg_ops={'win0': segs[0][2], 'lt': segs[nw][2], 'add': segs[nw + 1][2]},
                win1_ops=(segs[1][2] if nw > 1 else None), win0_ops=window_ops(segs, info, base0, 0), info=info, base0=base0,
                checks=checks, e2e=e2e)
    tmp = out_path + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(slim, f)
    os.replace(tmp, out_path)
    digest = sha256(open(out_path, 'rb').read()).hexdigest()
    open(out_path + '.sha256', 'w').write(digest)
    print('LOWER_DONE', dict(nq=meta['nq'], ops=meta['ops'], glue=checks['glue']['ok'], bookkeeping=checks['bookkeeping']['ok'],
                              windows_ok=checks['windows_ok'], prefix_ok=checks['prefix']['ok'], same_hash=checks['same_hash'],
                              same_structure=checks['same_structure'], e2e=e2e and {k: e2e.get(k) for k in ('cases', 'mismatches')}), flush=True)


def structural_checks(P, segs, info, base0, nw, low_io, tw, rwire):
    n, w = P.n, P.w
    low = low_io
    glue = glue_check(segs, low, n, nw, tw, rwire)
    bk = C.bookkeeping(segs, low, n, w)
    wins = []
    tiles = {}
    for i, wnf in enumerate(info):
        wo_ = window_ops(segs, info, base0, i)
        off = wnf['span'][0]
        t = SP.tile_info(wnf)
        if t: tiles[wnf['key']] = t
        try:
            subs = SP.build_sub(P, wnf['children'], Sl(wo_, off))
            chk = SP.window_checks(P, wnf, wo_, off, subs)
        except Exception as e:   # a changed structure of the window is a failed check, never a pass
            chk = dict(ok=False, notes=['structure of the window cannot be cut into sub-nodes: %r' % (e,)])
        wins.append(dict(key=wnf['key'], ok=bool(chk['ok']), notes=chk['notes']))
    # the operation X(reduced) before window 0 (CFIX F-CFIX-1): exactly [X(reduced)], one wire, untouched by window 0
    pn = []
    a_ = info[0]['span'][0] - base0
    rw = low['outputs']['reduced']
    s0 = segs[0][2]
    if s0[:a_] != [('x', rwire)] or len(rw) != 1 or rw[0] != rwire:
        pn.append('prefix ops before window 0 are not exactly [X(reduced)]: %s' % (s0[:a_][:3],))
    w0ops = s0[a_:]
    if rwire in C.support(w0ops, {}, {}) or any(rwire in v for v in info[0]['wins_in'].values()):
        pn.append('reduced wire touched by window 0')
    whash = [None] * nw
    for i in range(1, nw):
        label, kind, ops_, wi_, wo_, first = segs[i]
        nq_, c_, a2, b2, _m = C.compact(ops_, wi_, wo_)
        whash[i] = (sha256(C.emit_qasm(nq_, c_).encode()).hexdigest(), C.digest_json(C.make_layout(nq_, a2, b2, c_)))
    same_hash = all(h == whash[1] for h in whash[1:]) if nw > 1 else True
    sig = lambda wv: [(c['name'], c['b'] - c['a']) for c in wv['children']]
    same_structure = all(sig(wv) == sig(info[1]) for wv in info[1:]) if nw > 1 else True
    return dict(glue=glue, bookkeeping=bk, windows=wins, windows_ok=all(x['ok'] for x in wins), tiles_bad=tiles,
                prefix=dict(ok=not pn, notes=pn), same_hash=bool(same_hash), same_structure=bool(same_structure),
                generic_window_hash=(whash[1] if nw > 1 else None), n_windows=nw)


def stage_inv(n, w, p, part, out_path):
    P = C.Params(n, w, p)
    r = C.invariant_checks(P, (part == 'first',))
    json.dump(r, open(out_path, 'w'))
    print('INV_DONE', part, r, flush=True)


def _brief(r):
    keys = ('label', 'kind', 'status', 'detail', 'ops', 'nq', 'wall', 'peak_mb', 'qasm_sha256', 'contract_sha256', 'domain', 'sim_cases', 'sim_bad', 'dedup_of', 'witness')
    return {k: r.get(k) for k in keys if k in r}


def stage_prove(in_path, out_path, budget, max_mb, seed):
    D = pickle.load(open(in_path, 'rb'))
    n, w, p, nw = D['n'], D['w'], D['p'], D['nw']
    P = C.Params(n, w, p)
    rng = random.Random(seed)
    info = D['info']
    tag = 'Qualtran %s (%s)' % (D['qualtran']['version'], 'patch ' + D['qualtran']['patch_kind'] if D['qualtran']['patch_kind'] else 'original')
    cache = {}
    out = dict(whole=[], win0_split=None, generic=None, reasons=[])
    meta = {s[0]: s for s in D['seg_meta']}   # label -> (label, kind, wi, wo, first)

    def whole(label):
        _, kind, wi, wo, first = meta[label]
        t = time.time()
        r = C.run_segment(label, kind, P, D['seg_ops'][label], wi, wo, first, rng, budget, max_mb, tag, cache)
        print(label, {k: r.get(k) for k in ('status', 'ops', 'nq', 'wall', 'domain', 'sim_bad')}, flush=True)
        out['whole'].append(_brief(r))
        return r

    # ---- window 0 (with X(reduced)): the whole segment where it fits (w <= 2), else (or when UNKNOWN) by sub-nodes
    r0 = whole('win0') if w <= 2 else None
    if r0 is None or r0['status'] == 'UNKNOWN':
        w0 = info[0]
        o0 = D['win0_ops']; off0 = w0['span'][0]
        segs0 = SP.build_sub(P, w0['children'], Sl(o0, off0))
        rr = []
        for (label, kind, o, swi, swo) in segs0:
            r = SP.run_sub('win0_' + label, kind, P, o, swi, swo, rng, budget, max_mb, tag, cache)
            print('win0_' + label, {k: r.get(k) for k in ('status', 'ops', 'nq', 'wall', 'domain', 'sim_bad', 'dedup_of')}, flush=True)
            rr.append(_brief(r))
        c0 = SP.compose_analytic(P, w0, segs0, rng, first=True)
        ok0 = all(r['status'] == 'PROVED' for r in rr) and c0['ok']
        out['win0_split'] = dict(results=rr, compose=c0, ok=ok0)
        print('win0_split', ok0, c0, flush=True)
    out['win0_ok'] = (r0 is not None and r0['status'] == 'PROVED') or bool(out['win0_split'] and out['win0_split']['ok'])
    if (r0 is not None and r0['status'] != 'PROVED' and r0['status'] != 'UNKNOWN'):
        out['reasons'].append('window 0: whole-segment proof status %s' % r0['status'])
    elif not out['win0_ok']:
        out['reasons'].append('window 0 not proved')
    # ---- LessThanConstant and the controlled AddK
    rl = whole('lt'); ra = whole('add')
    for r, nm in ((rl, 'LessThanConstant'), (ra, 'controlled AddK')):
        if r['status'] != 'PROVED': out['reasons'].append('%s: %s' % (nm, r['status']))
    # ---- generic window (window 1; windows 2.. have the same hash, checked in the lower stage) by its sub-nodes
    if nw > 1:
        win = info[1]
        o1 = D['win1_ops']; off = win['span'][0]
        segs1 = SP.build_sub(P, win['children'], Sl(o1, off))
        res = []
        for (label, kind, o, swi, swo) in segs1:
            r = SP.run_sub(label, kind, P, o, swi, swo, rng, budget, max_mb, tag, cache)
            print(label, {k: r.get(k) for k in ('status', 'detail', 'ops', 'nq', 'wall', 'domain', 'sim_cases', 'sim_bad', 'dedup_of')}, flush=True)
            res.append(_brief(r))
        comp = SP.compose_analytic(P, win, segs1, rng) if w >= 4 else SP.compose_check(P, win, segs1)
        okg = all(r['status'] == 'PROVED' for r in res) and comp['ok']
        out['generic'] = dict(window='window 1 (representative of windows 1..%d)' % (nw - 1), subnodes=res, compose=comp, ok=okg)
        print('compose', comp, flush=True)
        if not okg: out['reasons'].append('generic window: sub-nodes %s, composition ok=%s' % ([r['status'] for r in res], comp['ok']))
    json.dump(out, open(out_path, 'w'), indent=1, default=str)
    print('PROVE_DONE', dict(reasons=out['reasons']), flush=True)
