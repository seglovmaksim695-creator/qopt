# SPDX-License-Identifier: Apache-2.0
"""Merge of conditional CZ corrections acting on the same pair of wire forms (lab agent a7):
cz(u,v)^(C1) cz(u,v)^(C2) = cz(u,v)^(C1 xor C2), placed at the last occurrence.

Wire values are tracked as affine forms over GF(2) of path variables (H, measurement and reset start a fresh variable).
A CZ is diagonal, so the phase (-1)^(C u v) may be applied at any moment where two wires hold u and v and the bits of C are
known. Wires touched by a conditional non-diagonal gate are "poisoned" (their value is no longer an affine form).
Requires normalised bits (one measurement = one bit, see ir.normalize).
"""
from .ir import pcond, fcond, DIAG


def merge(ops, n):
    form = [1 << i for i in range(n)]; const = [0] * n; poison = [False] * n; nl = [n]

    def fresh(q):
        form[q] = 1 << nl[0]; nl[0] += 1; const[q] = 0; poison[q] = False

    terms = {}   # key -> list of (op index, gate index, condition bits)
    for i, o in enumerate(ops):
        if o[0] == 'G':
            nm, qs = o[1], o[2]
            if nm == 'cx':
                form[qs[1]] ^= form[qs[0]]; const[qs[1]] ^= const[qs[0]]; poison[qs[1]] |= poison[qs[0]]
            elif nm == 'x':
                const[qs[0]] ^= 1
            elif nm == 'h':
                fresh(qs[0])
            elif nm not in DIAG:
                for q in qs:
                    fresh(q); poison[q] = True
        elif o[0] in ('M', 'R'):
            fresh(o[1])
        else:
            C = pcond(o[1])
            for gi, (nm, qs) in enumerate(o[2]):
                if nm == 'cz' and len(C) == 1:
                    a, b = qs
                    if const[a] != 0 or const[b] != 0 or poison[a] or poison[b]:
                        continue
                    key = frozenset((form[a], form[b]))
                    if len(key) < 2:
                        continue
                    terms.setdefault(key, []).append((i, gi, frozenset(C[0])))
                elif nm not in DIAG:
                    # later gates of the same body see the changed wire
                    for q in (qs[1:] if nm == 'cx' else qs):
                        poison[q] = True
    ops = [list(o) if o[0] == 'IF' else o for o in ops]
    drop = set(); newcond = {}
    for key, lst in terms.items():
        if len(lst) < 2:
            continue
        cb = frozenset()
        for _, _, c in lst:
            cb = cb ^ c
        li, lg, _ = lst[-1]
        for (i, g, c) in lst[:-1]:
            drop.add((i, g))
        if not cb:
            drop.add((li, lg))
        else:
            newcond[(li, lg)] = fcond([tuple(sorted(cb))])
    out = []
    for i, o in enumerate(ops):
        if o[0] in ('G', 'M', 'R'):
            out.append(o); continue
        # split the body into runs with the same condition, keeping gate order
        cur_c, cur_b = None, []
        for gi, g in enumerate(o[2]):
            if (i, gi) in drop:
                continue
            c = newcond.get((i, gi), o[1])
            if c != cur_c and cur_b:
                out.append(('IF', cur_c, cur_b)); cur_b = []
            cur_c = c; cur_b.append(g)
        if cur_b:
            out.append(('IF', cur_c, cur_b))
    return out
