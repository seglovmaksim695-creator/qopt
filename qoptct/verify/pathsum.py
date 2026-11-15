# SPDX-License-Identifier: Apache-2.0
"""Exact verification of Clifford+T circuits with measurements and classical feed-forward.

Model. A circuit (qoptct OpenQASM 3 subset) acts on N qubits; wires 0..n-1 are the wires of the .qc reference, the
others are ancillas starting in |0>. A branch m (vector of all measurement outcomes) is an operator K_m. The circuit is
CORRECT (strict model) if for EVERY branch m:  K_m = c_m * (U_ref (x) |a(m)>), c_m a number (possibly 0) and a(m) a
basis state of the ancillas. Clean model: the same on inputs whose clean wires (.qc wires not on the .i line) are 0.
The unitary model is the special case without measurements.

Method: sum-over-paths (Amy, arXiv:1805.06908) with SYMBOLIC measurement outcomes:
  * state = 2^{e/2} sum_y w^{F(x,m,y)} |W(x,m,y)>,  w = exp(i pi/4); x inputs, m outcomes (free), y path variables;
  * wires W are F2 polynomials in ANF, the phase F is a Z8 polynomial over Boolean monomials (both canonical);
  * Z-measurement of wire v with outcome m is the factor delta(v + m) = 1/2 sum_z (-1)^{z (v + m)};
    reset = measurement with a symbolic outcome followed by a conditional X; `if (cond)` = gate controlled by
    the F2 polynomial of the condition;
  * U_ref^dagger is appended and the path sum is reduced with the exact rules Elim, HH and omega.
If no path variables remain, the normal form (data wires = x, ancillas = functions of m only, phase depends on m
only) is a PROOF for all x and m; any deviation is canonical and yields a concrete counterexample (x, m).
If the rules get stuck, the trace miter (qoptct.verify.miter) decides exactly. No step can produce a false PROVED.
"""
import random
import re
import time

from ._parse import (ONE, ZERO, bits, pmul, peval, cond_poly, parse_qasm, parse_qc, fold_h, fold_ops, nbits_of)

__all__ = ['PathSum', 'build', 'prove_exact']


class PathSum:
    def __init__(self):
        self.nv = 0
        self.kind = []          # 'x' input, 'm' outcome, 'y' path variable
        self.label = []
        self.ph = {}            # monomial -> coefficient mod 8 (phase w^F)
        self.w = []             # wires: frozenset of monomials
        self.guards = []        # constraints P = 0 (on free variables only)
        self.e2 = 0             # amplitude factor sqrt(2)^e2
        self.zero = False
        self.ys = set()
        self.stats = dict(h=0, meas=0, elim=0, hh=0, omega=0, guard=0)
        self.zpick = 'occ'
        self.frozen = set()     # cut variables, temporarily excluded from the rules
        self.desc = False
        self.cut_K = None       # ANF size threshold for cutting a wire (None = off)
        self.cuttab = {}        # polynomial -> cut variable
        self.deadline = None    # wall-clock limit (time.time() value) checked inside long loops

    def var(self, kind, label=None):
        v = self.nv
        self.nv += 1
        self.kind.append(kind)
        self.label.append(label)
        if kind == 'y':
            self.ys.add(v)
        return v

    # --- phase
    def addc(self, mo, c):
        v = (self.ph.get(mo, 0) + c) & 7
        if v:
            self.ph[mo] = v
        else:
            self.ph.pop(mo, None)

    def add_lift(self, P, c, extra=0):
        """F += c * [P] with [P] in {0,1} as an integer: XOR_i b_i = sum_S (-2)^{|S|-1} prod_S b_i (|S| <= 3 is enough mod 8)."""
        c &= 7
        if not c:
            return
        it = list(P)
        L = len(it)
        for i in range(L):
            self.addc(it[i] | extra, c)
        if c & 3:
            c2 = (-2 * c) & 7
            for i in range(L):
                if self.deadline is not None and L > 256 and time.time() > self.deadline:
                    raise TimeoutError('path-sum reduction exceeded the time budget')
                for j in range(i + 1, L):
                    self.addc(it[i] | it[j] | extra, c2)
            if c & 1:
                c3 = (4 * c) & 7
                for i in range(L):
                    if self.deadline is not None and L > 64 and time.time() > self.deadline:
                        raise TimeoutError('path-sum reduction exceeded the time budget')
                    for j in range(i + 1, L):
                        ij = it[i] | it[j] | extra
                        for k in range(j + 1, L):
                            self.addc(ij | it[k], c3)

    # --- gates (cond: F2 control polynomial, ONE = unconditional)
    def gate(self, name, qs, cond=ONE):
        w = self.w
        if name == 'x':
            w[qs[0]] = w[qs[0]] ^ cond
        elif name in ('cx', 'cnot'):
            c, t = qs
            w[t] = w[t] ^ pmul(cond, w[c])
            self._cut(t)
        elif name in ('ccx', 'toffoli', 'mcx'):
            *cs, t = qs
            if self.cut_K:
                est = len(cond)
                for c in cs:
                    est *= max(1, len(w[c]))
                if est > self.cut_K:
                    for c in cs:
                        if len(w[c]) > 1:
                            self._cut(c, force=True)
            P = cond
            for c in cs:
                P = pmul(P, w[c])
            w[t] = w[t] ^ P
            self._cut(t)
        elif name in ('z', 's', 'sdg', 't', 'tdg'):
            k = {'z': 4, 's': 2, 'sdg': 6, 't': 1, 'tdg': 7}[name]
            if self.cut_K and len(w[qs[0]]) * len(cond) > (8 if k & 1 else self.cut_K):
                self._cut(qs[0], force=True)
            self.add_lift(pmul(cond, w[qs[0]]), k)
        elif name in ('cz', 'ccz', 'mcz'):
            if self.cut_K:
                est = len(cond)
                for c in qs:
                    est *= max(1, len(w[c]))
                if est > self.cut_K:
                    for c in qs:
                        if len(w[c]) > 1:
                            self._cut(c, force=True)
            P = cond
            for q in qs:
                P = pmul(P, w[q])
            self.add_lift(P, 4)
        elif name == 'swap':
            a, b = qs
            if cond == ONE:
                w[a], w[b] = w[b], w[a]
            else:
                for (c, t) in ((a, b), (b, a), (a, b)):
                    self.gate('cx', (c, t), cond)
        elif name == 'h':
            if cond != ONE:
                raise ValueError('conditional H is not supported')
            self.hadamard(qs[0])
        elif name == 'y':   # Y = i X Z
            self.add_lift(cond, 2)
            self.gate('z', qs, cond)
            self.gate('x', qs, cond)
        elif name == 'id':
            pass
        else:
            raise ValueError('unsupported gate: ' + name)

    def _cut(self, q, force=False):
        """If the ANF of a wire exceeds cut_K: wire := c, with c = P enforced by delta(c + P) = 1/2 sum_z (-1)^{z(c+P)}
        (c and z are summed and frozen until the end of the construction). Equal polynomials share one variable."""
        K = self.cut_K
        if not K:
            return
        P = self.w[q]
        if len(P) <= (1 if force else K):
            return
        c = self.cuttab.get(P)
        if c is None:
            c = self.var('y', 'cut')
            z = self.var('y', 'cutz')
            zb = 1 << z
            for mo in P ^ frozenset([1 << c]):
                self.addc(mo | zb, 4)
            self.e2 -= 2
            self.frozen |= {c, z}
            self.cuttab[P] = c
            self.stats['cut'] = self.stats.get('cut', 0) + 1
        self.w[q] = frozenset([1 << c])

    def hadamard(self, q):
        y = self.var('y', 'h')
        yb = 1 << y
        for mo in self.w[q]:
            self.addc(mo | yb, 4)
        self.w[q] = frozenset([yb])
        self.e2 -= 1
        self.stats['h'] += 1

    def measure(self, q, label=None):
        """Z-measurement; returns the outcome variable."""
        mv = self.var('m', label)
        mb = 1 << mv
        z = self.var('y', 'meas')
        zb = 1 << z
        for mo in self.w[q] ^ frozenset([mb]):
            self.addc(mo | zb, 4)
        self.e2 -= 2
        self.w[q] = frozenset([mb])
        self.stats['meas'] += 1
        return mv

    # --- substitution of a path variable z := R
    def subst(self, z, R):
        zb = 1 << z
        hit = [(mo, c) for mo, c in self.ph.items() if mo & zb]
        for mo, c in hit:
            del self.ph[mo]
        dl = self.deadline
        for mo, c in hit:
            if dl is not None and time.time() > dl:
                raise TimeoutError('path-sum reduction exceeded the time budget')
            self.add_lift(R, c, mo ^ zb)

        def sp(P):
            if not any(mo & zb for mo in P):
                return P
            r = set()
            for mo in P:
                if mo & zb:
                    for k in pmul(frozenset([mo ^ zb]), R):
                        if k in r:
                            r.remove(k)
                        else:
                            r.add(k)
                else:
                    if mo in r:
                        r.remove(mo)
                    else:
                        r.add(mo)
            return frozenset(r)
        self.w = [sp(P) for P in self.w]
        self.guards = [sp(P) for P in self.guards]
        if self.cuttab:
            self.cuttab = {P: c for P, c in self.cuttab.items() if not any(mo & zb for mo in P)}

    # --- reduction (rules Elim, HH, omega)
    def reduce(self, max_pass=100000, deadline=None):
        if deadline is None:
            deadline = self.deadline
        for _ in range(max_pass):
            if self.zero:
                return
            if deadline is not None and time.time() > deadline:
                raise TimeoutError('path-sum reduction exceeded the time budget')
            changed = False
            busy = 0
            for P in self.w:
                for mo in P:
                    busy |= mo
            for P in self.guards:
                for mo in P:
                    busy |= mo
            ymask = sum(1 << y for y in self.ys)
            occ = {}
            for mo, c in self.ph.items():
                mm = mo & ymask
                while mm:
                    lo = mm & -mm
                    mm ^= lo
                    occ.setdefault(lo.bit_length() - 1, []).append(mo)
            dirty = 0
            for y in sorted(self.ys, reverse=self.desc):
                if deadline is not None and time.time() > deadline:
                    raise TimeoutError('path-sum reduction exceeded the time budget')
                yb = 1 << y
                if busy & yb or dirty & yb or y not in self.ys or y in self.frozen:
                    continue
                ms = occ.get(y, [])
                ms = [mo for mo in ms if mo in self.ph]
                if any(not (mo & yb) for mo in ms):
                    continue
                if not ms:
                    self.ys.discard(y)
                    self.e2 += 2
                    self.stats['elim'] += 1
                    changed = True
                    continue
                cs = [self.ph[mo] for mo in ms]
                if all(c == 4 for c in cs):
                    Q = frozenset(mo ^ yb for mo in ms)
                    if Q == ONE:
                        self.zero = True
                        return
                    cnt = {}
                    for mo in Q:
                        for v in bits(mo & ymask):
                            cnt[v] = cnt.get(v, 0) + 1
                    cand = [v for v, k in cnt.items() if k == 1 and (1 << v) in Q and v != y and v not in self.frozen]
                    if cand:
                        if self.zpick == 'new':
                            z = max(cand)
                        elif self.zpick == 'old':
                            z = min(cand)
                        else:
                            z = min(cand, key=lambda v: (len(occ.get(v, ())), v))
                        for mo in ms:
                            del self.ph[mo]
                        R = Q ^ frozenset([1 << z])
                        touched = 0
                        for mo in ms:
                            touched |= mo
                        for mo in occ.get(z, ()):
                            touched |= mo
                        for mo in R:
                            touched |= mo
                        self.subst(z, R)
                        self.ys.discard(y)
                        self.ys.discard(z)
                        self.e2 += 2
                        self.stats['hh'] += 1
                        changed = True
                        dirty |= touched
                        busy = 0
                        for P in self.w:
                            for mo in P:
                                busy |= mo
                        continue
                    if not any(mo & ymask for mo in Q):
                        # constraint on free variables only
                        for mo in ms:
                            del self.ph[mo]
                        self.guards.append(Q)
                        self.ys.discard(y)
                        self.e2 += 2
                        self.stats['guard'] += 1
                        changed = True
                        for mo in Q:
                            busy |= mo
                        continue
                    continue
                if yb in self.ph and self.ph[yb] in (2, 6) and all(self.ph[mo] == 4 for mo in ms if mo != yb):
                    c = self.ph[yb]
                    Q = frozenset(mo ^ yb for mo in ms if mo != yb)
                    for mo in ms:
                        del self.ph[mo]
                    self.addc(0, 1 if c == 2 else 7)
                    self.add_lift(Q, 6 if c == 2 else 2)
                    self.ys.discard(y)
                    self.e2 += 1
                    self.stats['omega'] += 1
                    changed = True
                    for mo in ms:
                        dirty |= mo
                    for mo in Q:
                        dirty |= mo
                    continue
            if not changed:
                return

    def size(self):
        return dict(ys=len(self.ys), ph=len(self.ph), wmax=max((len(P) for P in self.w), default=0),
                    guards=len(self.guards))


INV = {'t': 'tdg', 'tdg': 't', 's': 'sdg', 'sdg': 's'}


MAX_ANF = 200000   # monomials per wire; beyond this the construction is abandoned (UNKNOWN)


def build(qasm, qc, clean=False, reduce_every=1, zpick='occ', cut_K=None, fold=True, keep_frozen=False, deadline=None,
          max_anf=MAX_ANF):
    """Path sum of (U_ref^dagger (x) I) * circuit, reduced. -> (PathSum, info)."""
    names, ins, ref = parse_qc(qc)
    n = len(names)
    N, ops = parse_qasm(qasm)
    if fold:
        ref = fold_h(ref)
        ops = fold_ops(ops)
    if N < n:
        raise ValueError('the circuit has fewer qubits than the reference')
    ps = PathSum()
    ps.deadline = deadline
    ps.zpick = zpick
    ps.cut_K = cut_K
    xv = {}
    anc_clean = set()
    if clean and ins is not None:
        anc_clean = {i for i, nm in enumerate(names) if nm not in ins}
    for q in range(N):
        if q < n and q not in anc_clean:
            v = ps.var('x', names[q])
            xv[q] = v
            ps.w.append(frozenset([1 << v]))
        else:
            ps.w.append(ZERO)
    cvar = {}

    def cv(k):
        if k not in cvar:
            raise ValueError('condition on the unmeasured bit c[%d]' % k)
        return 1 << cvar[k]
    since = 0

    def tick():
        if deadline is not None and time.time() > deadline:
            raise TimeoutError('path-sum construction exceeded the time budget')
        if max_anf and max((len(P) for P in ps.w), default=0) > max_anf:
            raise TimeoutError('wire polynomial exceeds %d monomials' % max_anf)
    for op in ops:
        tick()
        if op[0] == 'G':
            name, qs = op[1], op[2]
            if max(qs, default=-1) >= N:
                raise ValueError('qubit index out of range')
            if name == 'reset':
                q = qs[0]
                mv = ps.measure(q, 'reset')
                ps.gate('x', (q,), frozenset([1 << mv]))
            else:
                ps.gate(name, qs)
        elif op[0] == 'M':
            _, q, k = op
            cvar[k] = ps.measure(q, 'c[%d]' % k)
            since += 1
            if reduce_every and since >= reduce_every:
                ps.reduce(deadline=deadline)
                since = 0
        else:
            cond = cond_poly(op[1], cv)
            for name, qs in op[2]:
                ps.gate(name, qs, cond)
    if reduce_every:
        ps.reduce(deadline=deadline)
    mid = ps.size()
    for name, qs in reversed(ref):
        tick()
        ps.gate(INV.get(name, name), qs)
    ps.reduce(deadline=deadline)
    if ps.frozen and not keep_frozen:
        ps.frozen = set()
        ps.cuttab = {}
        ps.desc = True
        ps.reduce(deadline=deadline)
        ps.desc = False
    return ps, dict(n=n, N=N, xv=xv, anc_clean=anc_clean, cvar=cvar, mid=mid)


def violations(ps, info):
    """Violations of the normal form (when no path variables are left) -> list of (kind, wire, polynomial/data)."""
    n, N, xv = info['n'], info['N'], info['xv']
    xmask = sum(1 << v for v in xv.values())
    out = []
    for q in range(N):
        P = ps.w[q]
        if q < n:
            want = frozenset([1 << xv[q]]) if q in xv else ZERO
            if P != want:
                out.append(('data', q, P ^ want))
        else:
            if any(mo & xmask for mo in P):
                out.append(('anc', q, P))
    xph = {mo: c for mo, c in ps.ph.items() if mo & xmask}
    if xph:
        out.append(('phase', None, xph))
    for G in ps.guards:
        if any(mo & xmask for mo in G):
            out.append(('guard', None, G))
    return out


def count_ok(ps):
    """Self-consistency: sum_m |c_m|^2 = 1 (when all constraints are affine). None if not applicable."""
    mvars = [v for v in range(ps.nv) if ps.kind[v] == 'm']
    rows = []
    for G in ps.guards:
        if any(len(bits(mo)) > 1 for mo in G):
            return None
        rows.append(sum(mo for mo in G if mo))
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            basis.append(r)
    nm = len(mvars) - len(basis)
    return ps.e2 + nm == 0, ps.e2 + nm


def _assign_search(ps, info, viol, rng, tries=4000, skip=0):
    """Find a concrete assignment (x, x', m) violating the normal form and satisfying the constraints."""
    xv = info['xv']
    xs = list(xv.values())
    ms = [v for v in range(ps.nv) if ps.kind[v] == 'm']
    typ, q, P = viol

    def rand_val(bias=None):
        val = 0
        for v in xs + ms:
            if rng.random() < 0.5:
                val |= 1 << v
        if bias is not None:
            val |= bias
        return val

    def guards_ok(val):
        return all(peval(G, val) == 0 for G in ps.guards)

    def phase(val):
        s = 0
        for mo, c in ps.ph.items():
            if mo & val == mo:
                s += c
        return s & 7
    xmask = sum(1 << v for v in xs)
    seeds = [None]
    if typ in ('data', 'anc', 'guard', 'phase'):
        seeds += sorted(P, key=lambda mo: bin(mo).count('1'))[:20]
    for t in range(skip, skip + tries):
        sd = seeds[t % len(seeds)]
        if sd is not None and t < len(seeds) * 3:
            val = sd if t < len(seeds) else rand_val(sd)
        else:
            val = rand_val()
        if not guards_ok(val):
            continue
        if typ == 'data':
            if peval(P, val):
                return val, val & ~xmask
        elif typ in ('anc', 'phase'):
            val2 = val & ~xmask
            if not guards_ok(val2):
                continue
            if typ == 'phase':
                if phase(val) != phase(val2):
                    return val, val2
            elif peval(P, val) != peval(P, val2):
                return val, val2
        elif typ == 'guard':
            val2 = val & ~xmask
            if guards_ok(val) != guards_ok(val2):
                return val, val2
            if peval(P, val) != peval(P, val2):
                return val, val2
    return None


def to_concrete(ps, info, val, nbits):
    xv = info['xv']
    x = {q: (val >> v) & 1 for q, v in xv.items()}
    m = [0] * nbits
    for k, v in info['cvar'].items():
        m[k] = (val >> v) & 1
    return x, m


def prove_exact(qasm, qc, clean=False, do_confirm=True, seed=1, use_trace=True, max_width=22, cex_tries=30,
                cut_K=None, keep_frozen=False, budget_sec=None):
    """Exact check for all inputs and all measurement outcomes.
    -> dict(verdict=PROVED|REFUTED|UNKNOWN|ERROR|INTERNAL, method=normal-form|trace|counterexample, ...).
    A REFUTED verdict from the normal form or the trace miter is exact; a counterexample (x, m) is attached and, when
    possible, confirmed by the independent branch simulator (qoptct.verify.branchsim)."""
    from .branchsim import confirm_branch
    t0 = time.time()
    deadline = None if not budget_sec else t0 + budget_sec
    res = {}
    try:
        names, ins, _ = parse_qc(qc)
        nb = nbits_of(qasm)
    except (ValueError, KeyError, IndexError, OSError) as e:
        return dict(verdict='ERROR', error=str(e), time=round(time.time() - t0, 2))
    n = len(names)
    cleanw = {i for i, nm in enumerate(names) if nm not in ins} if (clean and ins is not None) else set()
    xwires = [q for q in range(n) if q not in cleanw]

    def confirm(xs, m):
        try:
            return confirm_branch(qasm, qc, xs, m, n)
        except Exception as e:     # confirmation is auxiliary; never fatal
            res['confirm_error'] = repr(e)[:200]
            return None

    def random_cex(tries, until=None, hard=None):
        rng = random.Random(seed)
        for _ in range(tries):
            if ((hard or deadline) and time.time() > (hard or deadline)) or (until and time.time() > until):
                return
            xa = {q: rng.randint(0, 1) for q in xwires}
            xb = {q: rng.randint(0, 1) for q in xwires}
            m = [rng.randint(0, 1) for _ in range(nb)]
            dev = confirm([xa, xb], m)
            if dev is None:
                return
            if dev > 1e-6:
                res.update(counterexample=dict(x=[xa.get(q, 0) for q in range(n)], m=m),
                           confirm_dev=dev, confirmed=True)
                return

    try:
        # the construction gets 80% of the budget; the rest is kept for a counterexample search if it gives up
        ps, info = build(qasm, qc, clean=clean, cut_K=cut_K, keep_frozen=keep_frozen,
                         deadline=None if not budget_sec else t0 + 0.8 * budget_sec)
    except TimeoutError as e:
        res.update(verdict='UNKNOWN', error=str(e))
        if do_confirm:
            random_cex(cex_tries, hard=deadline or time.time() + 15.0)
            if res.get('confirmed'):
                res.update(verdict='REFUTED', method='counterexample')
        res['time'] = round(time.time() - t0, 2)
        return res
    except (ValueError, KeyError, IndexError) as e:
        return dict(verdict='ERROR', error=str(e), time=round(time.time() - t0, 2))
    res.update(n=info['n'], N=info['N'], meas=ps.stats['meas'], stats=ps.stats, size=ps.size(), e2=ps.e2)

    early = False
    if ps.ys and not ps.zero and do_confirm:
        # rules got stuck: a quick independent counterexample search first (a confirmed one settles the case,
        # and the trace miter can be expensive)
        random_cex(12, time.time() + (min(10.0, 0.15 * budget_sec) if budget_sec else 10.0))
        early = bool(res.get('confirmed'))
    tv = None
    if not ps.zero and use_trace and not early:
        from .miter import trace_check
        try:
            tv, tmeta = trace_check(ps, info, max_width=max_width, deadline=deadline)
        except TimeoutError:
            tv, tmeta = 'UNKNOWN', dict(why='time budget')
        res['trace'] = tv
        res['trace_meta'] = tmeta

    if ps.zero:
        res.update(verdict='REFUTED', method='normal-form', reason='the path sum is identically 0 (no branch is realised)')
    elif ps.ys:
        res['method'] = 'trace'
        res['ys_after_rules'] = len(ps.ys)
        res['verdict'] = tv if tv in ('PROVED', 'REFUTED', 'INTERNAL') else 'UNKNOWN'

        if early:
            res.update(verdict='REFUTED', method='counterexample')
        elif tv not in ('PROVED', 'REFUTED') and do_confirm:
            # the miter did not finish: look for an explicit counterexample on wires free of path variables
            ym = sum(1 << y for y in ps.ys)
            cand = [v for v in violations(ps, info) if v[0] in ('data', 'anc') and not any(mo & ym for mo in v[2])]
            for attempt in range(8):
                rng = random.Random(seed + attempt)
                cex = None
                for v in cand:
                    cex = _assign_search(ps, info, v, rng, skip=attempt * 37)
                    if cex:
                        break
                if not cex:
                    break
                xa, m = to_concrete(ps, info, cex[0], nb)
                xb, _ = to_concrete(ps, info, cex[1], nb)
                dev = confirm([xa, xb], m)
                if dev is None:
                    break
                if dev > 1e-6:
                    res.update(counterexample=dict(x=[xa.get(q, 0) for q in range(info['n'])], m=m),
                               confirm_dev=dev, confirmed=True)
                    break
            if not res.get('confirmed'):
                random_cex(cex_tries)
            if res.get('confirmed'):
                # only an independently confirmed counterexample turns an unfinished miter into REFUTED
                res.update(verdict='REFUTED', method='counterexample')
        elif tv == 'REFUTED' and do_confirm:
            random_cex(cex_tries)
    else:
        res['method'] = 'normal-form'
        viol = violations(ps, info)
        if not viol:
            co = count_ok(ps)
            res.update(verdict='PROVED', normalization=co)
            if co is not None and not co[0]:
                res.update(verdict='INTERNAL', reason='normalisation sum|c_m|^2 = 1 violated: %s' % (co,))
        else:
            res.update(verdict='REFUTED', reason=[(v[0], v[1]) for v in viol[:5]])
            res['counterexample'] = None
            for attempt in range(12 if do_confirm else 1):
                rng = random.Random(seed + attempt)
                cex = None
                for v in viol:
                    cex = _assign_search(ps, info, v, rng, skip=attempt * 37)
                    if cex:
                        break
                if not cex:
                    continue
                xa, m = to_concrete(ps, info, cex[0], nb)
                xb, _ = to_concrete(ps, info, cex[1], nb)
                res['counterexample'] = dict(x=[xa.get(q, 0) for q in range(info['n'])],
                                             x2=[xb.get(q, 0) for q in range(info['n'])], m=m)
                if not do_confirm:
                    break
                dev = confirm([xa, xb], m)
                if dev is None:
                    break
                res['confirm_dev'] = dev
                res['confirmed'] = dev > 1e-6
                if res['confirmed']:
                    break
        if tv in ('PROVED', 'REFUTED') and res['verdict'] in ('PROVED', 'REFUTED') and tv != res['verdict']:
            res['verdict'] = 'INTERNAL'
            res['reason'] = 'normal form and trace miter disagree'
    res['time'] = round(time.time() - t0, 2)
    return res
