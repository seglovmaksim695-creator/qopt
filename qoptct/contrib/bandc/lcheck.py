"""EXACT check of the compiler's logical operations for ALL inputs at once, independent of the compiler code.

Wires and the global phase are polynomials over GF(2) (ANF as sets of monomial masks; ANF is canonical, so equality
of functions == equality of ANFs). Semantics of the compiler ops:
  x / cx / z / cz / CCX / CCZ -- as usual (phase: exponent P of the sign (-1)^P);
  AND(a, b, t) -- requires W[t] == 0 identically; afterwards W[t] = W[a] W[b], no phase
                  (rccx + sdg gadget on |0>);
  MZ(a, corr)  -- X-measurement of a: factor (-1)^(m W[a]); if m = 1, corr is applied. Correct iff the phase of the
                  CZ/Z gates in corr (on the current values, without wire a) equals W[a] (up to a constant, which is a
                  global phase of the branch) and corr contains x a (reset to 0);
  CUT(w, _, k, t) -- the value of wire w is renamed to a fresh variable; its definition is taken from the
                  reference (sound: compared syntactically).
The reference is built directly from the .qc file (H-frames: h ccz h = ccx etc.); its phase is a polynomial too.
Result: data wires == reference, ancillas == 0, phase == reference phase (mod 2, exactly).
The lower level (to_qasm: rccx, peephole, phase folding) is NOT covered here -- see vpt.py / verify_qasm.
"""
import time


class TooBig(Exception):
    """the check exceeded its size/time budget (inconclusive, never a pass)"""


MAX_PRODUCT = 4_000_000     # |A|*|B| term pairs per multiplication
MAX_TERMS = 400_000         # monomials per polynomial


class P:
    """polynomial over GF(2): frozenset of monomial masks"""
    __slots__ = ('s',)

    def __init__(self, s=frozenset()): self.s = frozenset(s)

    def __xor__(self, o): return P(self.s ^ o.s)

    def __mul__(self, o):
        if len(self.s) * len(o.s) > MAX_PRODUCT: raise TooBig('polynomial product too large')
        r = set()
        for a in self.s:
            for b in o.s:
                m = a | b
                if m in r: r.discard(m)
                else: r.add(m)
        if len(r) > MAX_TERMS: raise TooBig('polynomial too large')
        return P(r)

    def __eq__(self, o): return self.s == o.s

    def __hash__(self): return hash(self.s)

    def iszero(self): return not self.s


ONE = P({0})


def var(i): return P({1 << i})


def load_qc(path):
    """independent .qc parser: (names, gates, inputs)"""
    names = None; ins = None; gates = []; body = False
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line: continue
        if line.startswith('.v'): names = line.split()[1:]; continue
        if line.startswith('.i'): ins = line.split()[1:]; continue
        if line.startswith('.o'): continue
        if line == 'BEGIN': body = True; continue
        if line == 'END': break
        if not body: continue
        t = line.split(); q = [names.index(x) for x in t[1:]]
        if t[0] == 'tof' or (t[0] == 'X' and len(q) == 1):
            gates.append(('tof', q))
        elif t[0] in ('Z', 'Zd'): gates.append(('z', list(dict.fromkeys(q))))
        elif t[0] == 'H' and len(q) == 1: gates.append(('h', q))
        else: raise ValueError('non-classical/unknown .qc gate: ' + line)
    return names, gates, ins


def reference(path, clean, cuts=None):
    """cuts: {k: t} -- after the k-th (non-H) gate the value of wire t is replaced by a fresh variable
    (id = n + index); also returns defs {k: (definition, id)}."""
    names, gates, ins = load_qc(path); n = len(names)
    W = [(var(i) if (not clean or ins is None or names[i] in ins) else P()) for i in range(n)]
    ph = P(); fr = [0] * n; defs = {}; kk = -1; cuts = cuts or {}
    for g, q in gates:
        if g == 'h': fr[q[0]] ^= 1; continue
        kk += 1
        if kk > 0 and (kk - 1) in cuts: _cut(W, cuts[kk - 1], kk - 1, defs, n)
        if g == 'tof':
            c, t = q[:-1], q[-1]
            if any(fr[x] for x in c): raise ValueError('control inside an H-frame')
            prod = ONE
            for x in c: prod = prod * W[x]
            if fr[t]: ph = ph ^ (prod * W[t])     # inside an H-frame of t: H X H = Z -> controlled Z
            else: W[t] = W[t] ^ prod
        elif g == 'z':
            F = [x for x in q if fr[x]]
            if not F:
                prod = ONE
                for x in q: prod = prod * W[x]
                ph = ph ^ prod
            elif len(F) == 1:
                t = F[0]; prod = ONE
                for x in q:
                    if x != t: prod = prod * W[x]
                W[t] = W[t] ^ prod
            else: raise ValueError('two H-frames in one gate')
        else: raise ValueError('gate ' + g)
    if kk in cuts: _cut(W, cuts[kk], kk, defs, n)
    if any(fr): raise ValueError('open H-frame')
    return W, ph, n, defs


def _cut(W, t, k, defs, n):
    y = n + len(defs); defs[k] = (W[t], y); W[t] = var(y)


def check_ops(ops, N, qc_path, clean=False, deadline=None):
    """-> {'ok': bool, 'why'?: str}. Raises TooBig when the budget is exceeded (treat as 'not verified')."""
    cuts = {q[2]: q[3] for g, q in ops if g == 'CUT'}
    Wr, phr, n, defs = reference(qc_path, clean, cuts)
    names, gates, ins = load_qc(qc_path)
    W = [(var(i) if (not clean or ins is None or names[i] in ins) else P()) for i in range(n)] + [P() for _ in range(N - n)]
    ph = P(); nand = nmz = 0
    for step, (g, q) in enumerate(ops):
        if deadline is not None and step % 64 == 0 and time.time() > deadline: raise TooBig('lcheck timeout')
        if g == 'x': W[q[0]] = W[q[0]] ^ ONE
        elif g == 'cx': W[q[1]] = W[q[1]] ^ W[q[0]]
        elif g == 'z': ph = ph ^ W[q[0]]
        elif g == 'cz': ph = ph ^ (W[q[0]] * W[q[1]])
        elif g == 'CCX': W[q[2]] = W[q[2]] ^ (W[q[0]] * W[q[1]])
        elif g == 'CCZ': ph = ph ^ (W[q[0]] * W[q[1]] * W[q[2]])
        elif g == 'AND':
            a, b, t = q
            if not W[t].iszero(): return dict(ok=False, why=f'AND target {t} not 0')
            if t in (a, b): return dict(ok=False, why='AND target = control')
            W[t] = W[a] * W[b]; nand += 1
        elif g == 'CUT':
            w, _, k, t = q
            if k not in defs or not (W[w] == defs[k][0]): return dict(ok=False, why=f'CUT {k}: wire {w} != definition')
            W[w] = var(defs[k][1])
        elif g == 'MZ':
            a, corr = q; cp = P(); rx = False
            for gg, qq in corr:
                if gg == 'x' and tuple(qq) == (a,): rx = True; continue
                if a in qq: return dict(ok=False, why='correction touches the measured wire')
                if gg == 'cz': cp = cp ^ (W[qq[0]] * W[qq[1]])
                elif gg == 'z': cp = cp ^ W[qq[0]]
                else: return dict(ok=False, why='correction gate ' + gg)
            if not rx: return dict(ok=False, why='no reset in correction')
            # a difference by the constant 1 is the phase (-1)^m, depending only on the outcome (global per branch)
            if not (cp == W[a] or (cp ^ ONE) == W[a]): return dict(ok=False, why=f'MZ {a}: correction != value')
            W[a] = P(); nmz += 1
        else: return dict(ok=False, why='op ' + g)
    for i in range(n):
        if not (W[i] == Wr[i]): return dict(ok=False, why=f'wire {i} mismatch')
    for i in range(n, N):
        if not W[i].iszero(): return dict(ok=False, why=f'ancilla {i} dirty')
    if not (ph == phr): return dict(ok=False, why='phase mismatch')
    return dict(ok=True, nand=nand, nmz=nmz)
