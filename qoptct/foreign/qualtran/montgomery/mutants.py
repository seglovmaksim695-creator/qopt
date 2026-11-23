"""Structural mutants of the lowered Montgomery multiplier, used by the regression tests.

Every mutant must make ``prove_montgomery`` return NOT_PROVED.  They are the holes found by the audit of round 3
(H1, H2) and the mutants written while closing them (M1..M5); the checks of coverage and freshness of the sub-node
wires (``_split.window_checks``) are what catches them.

====  =====================================================================================================
H1    the first internal wire of QROAMClean+Add+Adjoint of EVERY generic window is replaced by the window wire x_0
H2    an extra X on T[0] between the first CAdd and Xor (an operation that belongs to no sub-node), every window
M1    an extra X on T[0] at the very beginning of the window, before all sub-nodes
M1b   an extra X on T[0] at the very end of the window, after all sub-nodes
M2    the internal wire of QROAM+Add+Adjoint is the input wire y[0] of the window
M3    the internal wire of QROAM+Add+Adjoint is the internal wire of the first CAdd (touched earlier)
M4    the operations are not changed, but the last operation of the first CAdd is also counted in Xor
M5    the operation X(reduced) before window 0 acts on another wire
====  =====================================================================================================
"""
import copy

from ._core import compact, make_layout, emit_qasm, digest_json, sha256
from ._split import sub_ranges

TAGS = ('H1', 'H2', 'M1', 'M1b', 'M2', 'M3', 'M4', 'M5')


def ren(ops, mp):
    g = lambda q: mp.get(q, q)
    out = []
    for o in ops:
        if o[0] == 'x':
            out.append(('x', g(o[1])))
        elif o[0] == 'cx':
            out.append(('cx', g(o[1]), g(o[2])))
        elif o[0] == 'and':
            out.append(('and', g(o[1]), g(o[2]), g(o[3])))
        else:
            out.append(('erase', g(o[1]), tuple((c[0],) + tuple(g(q) for q in c[1:]) for c in o[2])))
    return out


def _wires_of(s):
    return {q for v in s[3].values() for q in v} | {q for v in s[4].values() for q in v}


def _first_anc(s, a, b):
    ws = _wires_of(s)
    return next(q for o in s[2][a:b] for q in ([o[1]] if o[0] in ('x', 'erase') else list(o[1:3]) if o[0] == 'cx' else list(o[1:4])) if q not in ws)


def apply(tag, P, segs, info, base0, nw):
    """return (segs, info) of the mutant ``tag``; ``segs`` entries are (label, kind, ops, wi, wo, first)"""
    if tag not in TAGS:
        raise ValueError('unknown mutant %r; known: %s' % (tag, TAGS))
    segs = list(segs)
    if tag == 'M5':
        s = segs[0]
        a_ = info[0]['span'][0] - base0
        assert a_ >= 1 and s[2][0][0] == 'x'
        other = next(q for q in s[3]['x'])   # a wire that is not the reduced wire
        segs[0] = (s[0], s[1], [('x', other)] + list(s[2][1:]), s[3], s[4], s[5])
        return segs, info
    off1 = info[1]['span'][0]
    R_ = {r[0]: r for r in sub_ranges(info[1])}
    rel = lambda lab: (R_[lab][1] - off1, R_[lab][2] - off1)
    qa, qb = rel('qaa')
    for i in range(1, nw):
        s = segs[i]; o = s[2]
        if tag == 'H1':
            o2 = o[:qa] + ren(o[qa:qb], {_first_anc(s, qa, qb): s[3]['x'][0]}) + o[qb:]
        elif tag == 'M2':
            o2 = o[:qa] + ren(o[qa:qb], {_first_anc(s, qa, qb): s[3]['y'][0]}) + o[qb:]
        elif tag == 'M3':
            ca, cb = rel('cadd0')
            o2 = o[:qa] + ren(o[qa:qb], {_first_anc(s, qa, qb): _first_anc(s, ca, cb)}) + o[qb:]
        elif tag == 'H2':
            cut = rel('cadd0')[1]; o2 = o[:cut] + [('x', s[3]['T'][0])] + o[cut:]
        elif tag == 'M1':
            o2 = [('x', s[3]['T'][0])] + o
        elif tag == 'M1b':
            o2 = o + [('x', s[3]['T'][0])]
        else:   # M4: operations unchanged
            o2 = o
        segs[i] = (s[0], s[1], o2, s[3], s[4], s[5])
    inf = info
    if tag in ('H2', 'M1', 'M4'):
        inf = copy.deepcopy(info)
        for wv in inf[1:]:
            if tag == 'M4':
                for c in wv['children']:
                    if c['name'] == 'Xor': c['a'] -= 1
                continue
            cb = [c for c in wv['children'] if c['name'] == 'CAdd'][0]['b'] if tag == 'H2' else wv['span'][0]
            for c in wv['children']:
                if c['a'] >= cb: c['a'] += 1; c['b'] += 1
    return segs, inf
