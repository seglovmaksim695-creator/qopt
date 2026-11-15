# SPDX-License-Identifier: Apache-2.0
"""CPO: optimiser of conditional phase corrections (batched uncomputation + XOR conditions), lab agent a10.

Corrections of measurement-based uncomputation are conditional diagonal phases (-1)^{m f}, f a quadratic form of wire values.
The pass carries them forward as polynomials of path variables (wire values are tracked affinely over GF(2): CX changes a value,
H / measurement introduce a fresh variable), collects pending corrections of different measurements and emits them at a moment
tau (between creation and the first H after which the form is no longer expressible), where
  * the form is rewritten in terms of the CURRENT wire values (the number of CZ depends on the moment);
  * corrections of different measurements on the same wire pair are merged into one CZ with an XOR condition."""
import re, sys

HDR = ("OPENQASM", "include", "qubit", "bit", "//")
DIAG = {'z', 's', 'sdg', 't', 'tdg', 'cz', 'ccz', 'rz', 'u1', 'p', 'id'}


def gate(s):
    s = s.strip().rstrip(";").strip(); name, args = s.split(None, 1)
    return name, tuple(int(x) for x in re.findall(r"q\[(\d+)\]", args))


def parse_cond(s):
    gs = []
    for g in s.split("&&"):
        g = g.replace("(", "").replace(")", "")
        gs.append([int(x) for x in re.findall(r"c\[(\d+)\]", g)])
    return gs


def parse(txt):
    """Bits are renamed: every measurement gets a fresh bit (reuse of c[k] in the input is allowed)."""
    head, items = [], []; cur = {}; nb = [0]
    for L in txt.splitlines():
        Ls = L.strip()
        if not Ls: continue
        if Ls.startswith(HDR):
            if re.match(r"bit\[\d+\] c;", Ls): continue
            head.append(Ls); continue
        m = re.match(r"c\[(\d+)\] = measure q\[(\d+)\];", Ls)
        if m:
            cur[int(m.group(1))] = nb[0]; items.append(('M', int(m.group(2)), nb[0])); nb[0] += 1; continue
        m = re.match(r"if \((.*)\) \{(.*)\}$", Ls)
        if m:
            body = [gate(b) for b in m.group(2).split(";") if b.strip()]
            items.append(('IF', [[cur[k] for k in g] for g in parse_cond(m.group(1))], body)); continue
        items.append(('G',) + gate(Ls))
    head.append(f"bit[{max(nb[0], 1)}] c;")
    return head, items


def tog(d, k):
    if k in d: d.discard(k)
    else: d.add(k)


def bits(m):
    r = []
    while m:
        lo = m & -m; r.append(lo.bit_length() - 1); m ^= lo
    return r


def prod_poly(a, b):
    """(affine a) * (affine b) -> set of monomials (frozensets of variables), a = (mask, c)."""
    f = set(); (ma, ca), (mb, cb) = a, b
    for i in bits(ma):
        for j in bits(mb): tog(f, frozenset((i, j)))
    if ca:
        for j in bits(mb): tog(f, frozenset((j,)))
    if cb:
        for i in bits(ma): tog(f, frozenset((i,)))
    if ca and cb: tog(f, frozenset())
    return f


def lin_poly(a):
    f = set(); ma, ca = a
    for i in bits(ma): tog(f, frozenset((i,)))
    if ca: tog(f, frozenset())
    return f


class Basis:
    """Wire values of a snapshot vals (list of (mask, c)) + lazy completion to a basis (symbols < 0).
    A variable is an XOR of symbols + const; a form is expressible if no completion symbol appears."""
    def __init__(self, vals, excl=(), order=None):
        self.piv = {}; self.nz = 0; self.E = {}
        for q in (order if order is not None else range(len(vals))):
            m, c = vals[q]
            if q in excl or not m: continue
            self._ins(m, {q}, c)

    def _ins(self, v, s, cc):
        while v:
            p = v.bit_length() - 1
            if p in self.piv:
                bv, bs, bc = self.piv[p]; v ^= bv; s = s ^ bs; cc ^= bc
            else:
                self.piv[p] = (v, s, cc); return

    def expr(self, i):
        if i in self.E: return self.E[i]
        v = 1 << i; s = set(); cc = 0
        while v:
            p = v.bit_length() - 1
            if p not in self.piv:
                self.nz += 1; self.piv[p] = (v, {-self.nz}, 0)
            bv, bs, bc = self.piv[p]; v ^= bv; s ^= bs; cc ^= bc
        self.E[i] = (s, cc)
        return self.E[i]

    def ok(self, f):
        pairs, sing = to_wires(f, self)
        return not any(a < 0 or b < 0 for a, b in pairs) and not any(q < 0 for q in sing)


def to_wires(f, B):
    """Polynomial over variables -> (pairs of symbols, single symbols) in basis B; the constant is dropped."""
    pairs, sing = set(), set()
    for mon in f:
        if len(mon) == 0: continue
        if len(mon) == 1:
            (i,) = mon; R, c = B.expr(i)
            for q in R: tog(sing, q)
            continue
        i, j = tuple(mon); (Ri, ci), (Rj, cj) = B.expr(i), B.expr(j)
        for a in Ri:
            for b in Rj:
                if a == b: tog(sing, a)
                else: tog(pairs, (min(a, b), max(a, b)))
        if ci:
            for b in Rj: tog(sing, b)
        if cj:
            for a in Ri: tog(sing, a)
    return pairs, sing


def merged(group, B):
    P, S = {}, {}
    for cond, f in group:
        pairs, sing = to_wires(f, B)
        for p in pairs: P[p] = P.get(p, frozenset()) ^ cond
        for q in sing: S[q] = S.get(q, frozenset()) ^ cond
    P = {p: c for p, c in P.items() if c}; S = {q: c for q, c in S.items() if c}
    return P, S


NORD = 3

import numpy as _np
_H = _np.array([[1, 1], [1, -1]]) / _np.sqrt(2)
_P1 = {'t': _np.exp(1j * _np.pi / 4), 'tdg': _np.exp(-1j * _np.pi / 4), 's': 1j, 'sdg': -1j, 'z': -1}
SDIAG = {'t', 'tdg', 's', 'sdg', 'z'}


def _sim(gates, Q, init):
    """State on qubits Q (list) for a basis input init (dict q->bit); gates: (name, qs)."""
    k = len(Q); pos = {q: i for i, q in enumerate(Q)}
    psi = _np.zeros(2 ** k, complex); idx = sum(init[q] << (k - 1 - pos[q]) for q in Q); psi[idx] = 1
    for g, qs in gates:
        psi = psi.reshape([2] * k)
        if g == 'h':
            psi = _np.moveaxis(_np.tensordot(_H, psi, axes=([1], [pos[qs[0]]])), 0, pos[qs[0]])
        elif g in _P1:
            sl = [slice(None)] * k; sl[pos[qs[0]]] = 1; psi = psi.copy(); psi[tuple(sl)] *= _P1[g]
        elif g == 'x':
            psi = _np.flip(psi, pos[qs[0]])
        elif g == 'cx':
            c, t = pos[qs[0]], pos[qs[1]]; psi = psi.copy()
            sl = [slice(None)] * k; sl[c] = 1
            sub = psi[tuple(sl)]; tt = t if t < c else t - 1
            psi[tuple(sl)] = _np.flip(sub, tt)
        elif g == 'cz':
            a, b = pos[qs[0]], pos[qs[1]]; psi = psi.copy()
            sl = [slice(None)] * k; sl[a] = 1; sl[b] = 1; psi[tuple(sl)] *= -1
        else:
            return None
        psi = psi.reshape(-1)
    return psi


def detect_and(items, t0, tq, maxc=3, maxlen=40):
    """Item t0 = 'h tq' on a wire in |0>. Detects a logical-AND gadget up to the second h tq (+ trailing phases on tq).
    -> (included indices, window end, C (list), ANF monomials over C) or None"""
    inc = [t0]; gates = [('h', (tq,))]; C = []; end = None
    for j in range(t0 + 1, min(len(items), t0 + maxlen)):
        it = items[j]
        if it[0] != 'G':
            qs = {it[1]} if it[0] == 'M' else {q for _, b in [(0, x) for x in it[2]] for q in b[1]}
            if tq in qs or qs & set(C): return None
            continue
        g, qs = it[1], it[2]
        if tq in qs:
            for q in qs:
                if q != tq and q not in C: C.append(q)
            inc.append(j); gates.append((g, qs))
            if g == 'h': end = j; break
        elif set(qs) & set(C):
            for q in qs:
                if q not in C: C.append(q)
            inc.append(j); gates.append((g, qs))
        if len(C) > maxc: return None
    if end is None: return None
    span = end
    for j in range(end + 1, min(len(items), end + 12)):
        it = items[j]
        if it[0] != 'G':
            if it[0] == 'M' and it[1] == tq: break
            if it[0] == 'IF' and any(tq in qs for _, qs in it[2]): break
            continue
        g, qs = it[1], it[2]
        if tq in qs:
            if g in SDIAG and len(qs) == 1: inc.append(j); gates.append((g, qs)); span = j; continue
            break
        if set(qs) & set(C): break
    Q = [tq] + C; k = len(C); ph = None; g_tab = {}
    for cb in range(2 ** k):
        init = {tq: 0}; init.update({C[i]: (cb >> i) & 1 for i in range(k)})
        psi = _sim(gates, Q, init)
        if psi is None: return None
        nz = _np.nonzero(_np.abs(psi) > 1e-9)[0]
        if len(nz) != 1: return None
        idx = nz[0]; bits = {Q[i]: (idx >> (len(Q) - 1 - i)) & 1 for i in range(len(Q))}
        if any(bits[C[i]] != init[C[i]] for i in range(k)): return None
        a = psi[idx]
        if ph is None: ph = a
        elif abs(a - ph) > 1e-9: return None
        g_tab[cb] = bits[tq]
    f = [g_tab[cb] for cb in range(2 ** k)]
    for i in range(k):
        for cb in range(2 ** k):
            if cb >> i & 1: f[cb] ^= f[cb ^ (1 << i)]
    anf = []
    for cb in range(2 ** k):
        if f[cb]:
            mon = tuple(i for i in range(k) if cb >> i & 1)
            if len(mon) > 2: return None
            anf.append(mon)
    return inc, span, C, anf


def optimize(txt, maxwin=120, n=None, use_and=True):
    head, items = parse(txt)
    N = int(re.search(r"qubit\[(\d+)\]", txt).group(1))
    nv = [0]
    def fresh():
        nv[0] += 1; return (1 << (nv[0] - 1), 0)
    val = [fresh() if (n is None or q < n) else (0, 0) for q in range(N)]
    rel = {}
    blocked = set()
    skip = set()
    snaps = []
    ins = {}
    pend = []                      # [cond(frozenset), f(set), t0]
    measq = {}
    stats = dict(before=0, after=0)

    Bc = {}
    def basis_at(tau, o=0):
        if (tau, o) not in Bc:
            vals = snaps[tau] if tau < len(snaps) else list(val)
            order = None if o == 0 else (list(range(len(vals)))[::-1] if o == 1 else
                                          sorted(range(len(vals)), key=lambda q: (bin(vals[q][0]).count('1'), q)))
            Bc[(tau, o)] = Basis(vals, order=order)
        return Bc[(tau, o)]

    def best_sub(B, f):
        """Best substitution z -> P (z = product held by a wire): -> (cost, f') or None"""
        if not rel: return (0, f) if B.ok(f) else None
        vs = set().union(*f) if f else set()
        R = [z for z in rel if z in vs]
        for z in rel:
            if z not in R and B.expr(z)[0] and all(x >= 0 for x in B.expr(z)[0]): R.append(z)
        R = R[:7]
        best = None
        for msk in range(1 << len(R)):
            g = set(f)
            for i, z in enumerate(R):
                if msk >> i & 1:
                    tog(g, frozenset((z,)))
                    for mon in rel[z]: tog(g, mon)
            pairs, sing = to_wires(g, B)
            if any(a < 0 or b < 0 for a, b in pairs) or any(q < 0 for q in sing): continue
            if best is None or len(pairs) < best[0]: best = (len(pairs), g)
        return best

    def cost_at(tau, fs):
        best = None
        for o in range(NORD):
            B = basis_at(tau, o)
            fs2 = []
            for c, f in fs:
                r = best_sub(B, f)
                if r is None: break
                fs2.append((c, r[1]))
            else:
                P, S = merged(fs2, B)
                if best is None or len(P) < best[0]: best = (len(P), P, S)
        return best if best is not None else (10 ** 6, {}, {})

    def flush(group_idx, tf):
        group = [pend[i] for i in group_idx]
        if len(Bc) > 4000: Bc.clear()
        wins = []
        for c, f, t0 in group:
            cands = list(range(t0, tf + 1))
            if len(cands) > maxwin:
                st_ = len(cands) / maxwin
                cands = sorted({cands[int(i * st_)] for i in range(maxwin)} | set(cands[:8]) | set(cands[-8:]))
            w = [tau for tau in cands if tau not in blocked and best_sub(basis_at(tau), f) is not None]
            assert w, "correction not expressible at any moment of its window"
            wins.append(w)
        place = []
        for (c, f, t0), w in zip(group, wins):
            place.append(min(w, key=lambda tau: (cost_at(tau, [(c, f)])[0], -tau)))
        def total(pl):
            byt = {}
            for k, tau in enumerate(pl): byt.setdefault(tau, []).append((group[k][0], group[k][1]))
            return sum(cost_at(tau, fs)[0] for tau, fs in byt.items())
        cur = total(place); imp = True
        while imp:
            imp = False
            taus = sorted(set(place))
            for k in range(len(group)):
                ws = set(wins[k])
                for tau in taus:
                    if tau == place[k] or tau not in ws: continue
                    old = place[k]; place[k] = tau; c2 = total(place)
                    if c2 < cur: cur = c2; imp = True; break
                    place[k] = old
            common = set(wins[0]).intersection(*map(set, wins[1:])) if len(wins) > 1 else set()
            for tau in sorted(common):
                c2 = total([tau] * len(group))
                if c2 < cur: cur = c2; place = [tau] * len(group); imp = True
        byt = {}
        for k, tau in enumerate(place): byt.setdefault(tau, []).append((group[k][0], group[k][1]))
        for tau, fs in byt.items():
            _, P, S = cost_at(tau, fs)
            bycond = {}
            for p, c in P.items(): bycond.setdefault(c, []).append(('cz', p))
            for q, c in S.items(): bycond.setdefault(c, []).append(('z', (q,)))
            lines = []
            for c, gl in sorted(bycond.items(), key=lambda kv: sorted(kv[0])):
                cond = " ^ ".join(f"c[{k}]" for k in sorted(c))
                body = " ".join(f"{g} " + ", ".join(f"q[{x}]" for x in qs) + ";" for g, qs in gl)
                lines.append(f"if ({cond}) {{{body}}}")
            ins.setdefault(tau, []).extend(lines)
            stats['after'] += len(P)
        for i in sorted(group_idx, reverse=True): del pend[i]

    def flush_broken(t, excl):
        B = Basis(val, excl)
        bad = [i for i, (c, f, t0) in enumerate(pend) if best_sub(B, f) is None]
        if bad: flush(bad, t)

    def flush_all(t):
        if pend: flush(list(range(len(pend))), t)

    out_items = []
    for t, it in enumerate(items):
        snaps.append(list(val))
        keep = it
        if t in skip:
            out_items.append(keep); continue
        if use_and and it[0] == 'G' and it[1] == 'h' and val[it[2][0]] == (0, 0):
            r = detect_and(items, t, it[2][0])
            if r is not None:
                inc, span, C, anf = r
                z = fresh(); zi = nv[0] - 1
                P = set()
                for mon in anf:
                    if len(mon) == 0: tog(P, frozenset())
                    elif len(mon) == 1: [tog(P, m) for m in lin_poly(val[C[mon[0]]])]
                    else: [tog(P, m) for m in prod_poly(val[C[mon[0]]], val[C[mon[1]]])]
                rel[zi] = P
                val[it[2][0]] = z
                skip.update(inc); blocked.update(range(t + 1, span + 1))
                out_items.append(keep); continue
        if it[0] == 'G':
            g, qs = it[1], it[2]
            if g == 'cx': val[qs[1]] = (val[qs[1]][0] ^ val[qs[0]][0], val[qs[1]][1] ^ val[qs[0]][1])
            elif g == 'x': val[qs[0]] = (val[qs[0]][0], val[qs[0]][1] ^ 1)
            elif g in DIAG: pass
            elif g == 'swap': val[qs[0]], val[qs[1]] = val[qs[1]], val[qs[0]]
            elif g == 'h':
                flush_broken(t, {qs[0]}); val[qs[0]] = fresh()
            else:
                flush_all(t)
                for q in qs: val[q] = fresh()
        elif it[0] == 'M':
            q, k = it[1], it[2]
            flush_broken(t, {q}); val[q] = fresh(); measq[k] = q
        else:
            groups, body = it[1], it[2]
            simple = len(groups) == 1
            cond = frozenset(groups[0]) if simple else None
            if simple and len(groups[0]) != len(set(groups[0])): simple = False
            phase = [(g, qs) for g, qs in body if g in ('cz', 'z')]
            rest = [(g, qs) for g, qs in body if g not in ('cz', 'z')]
            ok_rest = all(g == 'x' and len(cond or ()) == 1 and qs[0] == measq.get(next(iter(cond))) for g, qs in rest) if simple else False
            if simple and ok_rest and phase:
                f = set()
                for g, qs in phase:
                    stats['before'] += (g == 'cz')
                    pf = prod_poly(val[qs[0]], val[qs[1]]) if g == 'cz' else lin_poly(val[qs[0]])
                    for mon in pf: tog(f, mon)
                pend.append([cond, f, t + 1])
                keep = ('IF', groups, rest) if rest else None
                for g, qs in rest:
                    flush_broken(t, {qs[0]}); val[qs[0]] = (0, 0)
            else:
                stats['before'] += sum(1 for g, _ in body if g == 'cz'); stats['after'] += sum(1 for g, _ in body if g == 'cz')
                if any(g not in DIAG for g, _ in body):
                    flush_all(t)
                    for g, qs in body:
                        if g not in DIAG:
                            for q in qs: val[q] = fresh()
        out_items.append(keep)
    snaps.append(list(val))
    for b in blocked:
        pass
    flush_all(len(items))
    lines = list(head)
    def fmt(it):
        if it[0] == 'M': return f"c[{it[2]}] = measure q[{it[1]}];"
        if it[0] == 'G': return f"{it[1]} " + ", ".join(f"q[{x}]" for x in it[2]) + ";"
        cond = " && ".join(" ^ ".join(f"c[{k}]" for k in g) for g in it[1])
        body = " ".join(f"{g} " + ", ".join(f"q[{x}]" for x in qs) + ";" for g, qs in it[2])
        return f"if ({cond}) {{{body}}}"
    for t, it in enumerate(out_items):
        for L in ins.get(t, []): lines.append(L)
        if it is not None: lines.append(fmt(it))
    for L in ins.get(len(items), []): lines.append(L)
    return "\n".join(lines) + "\n", stats


def count(txt):
    return len(re.findall(r"\b(t|tdg) q", txt)), len(re.findall(r"\b(cx|cz) q", txt))
