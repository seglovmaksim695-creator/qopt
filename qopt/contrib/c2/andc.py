"""c2: компилятор «логическое И + стирание измерением» для обратимых схем (уровень Тоффоли).
Идея: каждое нелинейное произведение L_a*L_b считается логическим И в свежую анциллу (4 T, 3 CX),
стирается X-измерением с клиффордовой поправкой (CZ/Z по ТЕКУЩИМ значениям проводов, 0 T).
Значения всех проводов отслеживаются как многочлены над GF(2) от входов (точная алгебра), поэтому
переиспользование произведений, «перенос носителя» логического провода в анциллу, линейные выражения
через живые провода и квадратичные поправки находятся автоматически.
Данные (все провода .qc) могут быть в любом состоянии; доп. анциллы стартуют в |0> (как в verify_qasm)."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, random, re
sys.path[:0] = [_QPKG, (_QPKG + '/portfolio')]
import qopt, evo_lib


class Fail(Exception): pass


def logical(hl, n):
    """H-рамки -> логическая классическая схема из x,z,cx,cz,ccx,ccz."""
    fr = [0] * n; out = []
    for g, q in hl:
        q = tuple(q)
        if g == 'h': fr[q[0]] ^= 1; continue
        if g in ('ccz', 'cz', 'z'):
            F = [x for x in q if fr[x]]
            if not F: out.append((g, q))
            elif len(F) == 1:
                t = F[0]; c = tuple(x for x in q if x != t)
                out.append(({'ccz': 'ccx', 'cz': 'cx', 'z': 'x'}[g], c + (t,)))
            else: raise Fail('2 frames')
        elif g == 'x':
            out.append(('z' if fr[q[0]] else 'x', q))
        elif g in ('cx', 'ccx'):
            if any(fr[x] for x in q[:-1]): raise Fail('ctrl in frame')
            out.append((({'cx': 'cz', 'ccx': 'ccz'}[g]) if fr[q[-1]] else g, q))
        else: raise Fail(g)
    if any(fr): raise Fail('open frame')
    return out


# ---------- многочлены: int-битсет по мономам (моном = маска переменных) ----------
ONE = 1
def var(i): return 1 << (1 << i)
_mc = {}
def mul(p, q):
    if p > q: p, q = q, p
    k = (p, q)
    r = _mc.get(k)
    if r is not None: return r
    A = []; x = p
    while x:
        lo = x & -x; A.append(lo.bit_length() - 1); x ^= lo
    B = []; x = q
    while x:
        lo = x & -x; B.append(lo.bit_length() - 1); x ^= lo
    r = 0
    for a in A:
        for b in B: r ^= 1 << (a | b)
    if len(_mc) > 200000: _mc.clear()
    _mc[k] = r
    return r


def solve(target, cols):
    """target как XOR подмножества cols (список (ключ, вектор)); порядок cols = приоритет. -> список ключей или None."""
    basis = {}  # pivot -> (vec, set_of_keys)
    for key, v in cols:
        s = {key}
        while v:
            p = v.bit_length() - 1
            if p in basis:
                bv, bs = basis[p]; v ^= bv; s = s ^ bs
            else:
                basis[p] = (v, s); break
    v = target; s = set()
    while v:
        p = v.bit_length() - 1
        if p not in basis: return None
        bv, bs = basis[p]; v ^= bv; s ^= bs
    return list(s)


def solve_min(target, cols, wt):
    """как solve, но затем уменьшает вес wt(ключ) решения перебором векторов ядра (жадно)."""
    basis = {}; kernel = []
    for key, v in cols:
        s = {key}
        while v:
            p = v.bit_length() - 1
            if p in basis:
                bv, bs = basis[p]; v ^= bv; s = s ^ bs
            else:
                basis[p] = (v, s); break
        if not v: kernel.append(frozenset(s))
    v = target; s = set()
    while v:
        p = v.bit_length() - 1
        if p not in basis: return None
        bv, bs = basis[p]; v ^= bv; s ^= bs
    W = lambda z: sum(wt(k) for k in z)
    cur = W(s); imp = True
    while imp:
        imp = False
        for k in kernel:
            z = s ^ k; w = W(z)
            if w < cur: s, cur, imp = z, w, True
    return list(s)


class Compiler:
    def __init__(self, n, lg, max_anc, rng, eps=0.0, keep_prob=1.0, lazy=True, linmax=3, inputs=None):
        self.n, self.lg, self.rng, self.eps, self.keep_prob = n, lg, rng, eps, keep_prob
        self.lazy, self.linmax, self.pend, self.guard, self.k = lazy, linmax, {}, False, 0
        self.N = n + max_anc
        # clean-модель: провода не из inputs стартуют в |0> (значение 0), переменные — только входы
        if inputs is None: inputs = list(range(n))
        vi = {q: k for k, q in enumerate(inputs)}
        init = [var(vi[i]) if i in vi else 0 for i in range(n)]
        self.V = init + [0] * max_anc
        self.L = list(init)
        self.car = list(range(n))
        self.live = set(range(n))   # физические провода с данными (данные всегда живые)
        self.free = list(range(n, self.N))
        self.ops = []; self.nand = 0; self.used = set()
        # будущие нужды
        L = list(self.L); seq = []
        for g, q in lg:
            need = set()
            if g in ('ccx', 'ccz'):
                a, b, c = q
                need |= {mul(L[a], L[b]), mul(L[a], L[c]), mul(L[b], L[c])}
                if g == 'ccx': L[c] ^= mul(L[a], L[b])
            elif g == 'cx': L[q[1]] ^= L[q[0]]
            elif g == 'x': L[q[0]] ^= ONE
            need |= set(L)
            seq.append(need)
        self.final = list(L)
        fut = [set(L)]
        for s in reversed(seq): fut.append(fut[-1] | s)
        self.fut = list(reversed(fut))  # fut[k] = нужды начиная с шага k (включая итог)

    def carriers(self): return set(self.car)

    def emit(self, g, q): self.ops.append((g, q))

    def lin_rep(self, target, excl=(), maxw=None):
        ws = sorted(w for w in self.live if w not in excl)
        if target == 0: return [], False
        if target == ONE: return [], True
        for t2, cst in ((target, False), (target ^ ONE, True)):
            for w in ws:
                if self.V[w] == t2: return [w], cst
        best = None
        for t2, cst in ((target, False), (target ^ ONE, True)):
            d = {}
            for w in ws: d.setdefault(self.V[w], w)
            for w in ws:
                u = d.get(t2 ^ self.V[w])
                if u is not None and u != w: return [w, u], cst
        if maxw is not None and maxw <= 2: return None
        for t2, cst in ((target, False), (target ^ ONE, True)):
            r = solve(t2, [(w, self.V[w]) for w in ws])
            if r is not None and (best is None or len(r) < len(best[0])): best = (r, cst)
        if best and maxw is not None and len(best[0]) > maxw: return None
        return best

    def quad_rep(self, target, excl):
        """target = const + sum lin + sum quad по живым проводам кроме excl -> (lin, pairs) или None."""
        r = self.lin_rep(target, excl)
        if r is not None: return r[0], []
        ws = sorted(w for w in self.live if w not in excl)
        cols = [(('l', w), self.V[w]) for w in ws] + [(('c',), ONE)]
        for i, a in enumerate(ws):
            for b in ws[i + 1:]:
                cols.append((('q', a, b), mul(self.V[a], self.V[b])))
        s = solve_min(target, cols, lambda k: 1 if k[0] == 'q' else 0)
        if s is None: return None
        return [k[1] for k in s if k[0] == 'l'], [(k[1], k[2]) for k in s if k[0] == 'q']

    def mz(self, a):
        r = self.quad_rep(self.V[a], {a})
        if r is None: return False
        lin, pairs = r
        corr = [('cz', p) for p in pairs] + [('z', (w,)) for w in lin] + [('x', (a,))]
        self.emit('MZ', (a, tuple(corr)))
        self.V[a] = 0; self.live.discard(a); self.free.append(a)
        return True

    def garbage(self):
        cs = self.carriers()
        return [w for w in self.live if w >= self.n and w not in cs]

    def ok_state(self, skip=None):
        """весь мусор (кроме skip) выразим квадратично, все отложенные дельты представимы."""
        for a in self.garbage():
            if a == skip: continue
            if self.quad_rep(self.V[a], {a}) is None: return False
        for t, D in self.pend.items():
            if D and self.car[t] != skip and self.quad_rep(D, {self.car[t]}) is None: return False
        return True

    def safe_mz(self, a):
        if self.quad_rep(self.V[a], {a}) is None: return False
        self.live.discard(a); ok = self.ok_state(); self.live.add(a)
        if not ok: return False
        return self.mz(a)

    def mz_force(self, a):
        """стереть a; сначала реализовать отложенные дельты и стереть мусор, которым без a не выразиться."""
        self.live.discard(a)
        dep_p = [t for t, D in self.pend.items() if D and self.quad_rep(D, {self.car[t]}) is None]
        dep_g = [b for b in self.garbage() if self.quad_rep(self.V[b], {b}) is None]
        self.live.add(a)
        for t in dep_p: self.realize(t, self.k)
        for b in dep_g:
            if b in self.live and b not in self.carriers(): self.mz(b)
        if a in self.live: self.mz(a)

    def before_modify(self, w, newval):
        if self.guard: return
        old = self.V[w]
        bad_g, bad_p = [], []
        pre_g = {a: self.quad_rep(self.V[a], {a}) is not None for a in self.garbage() if a != w}
        pre_p = {t: self.quad_rep(D, {self.car[t]}) is not None for t, D in self.pend.items() if D and self.car[t] != w}
        self.V[w] = newval
        for a, ok in pre_g.items():
            if ok and self.quad_rep(self.V[a], {a}) is None: bad_g.append(a)
        for t, ok in pre_p.items():
            if ok and self.quad_rep(self.pend[t], {self.car[t]}) is None: bad_p.append(t)
        self.V[w] = old
        self.guard = True
        try:
            for t in bad_p: self.realize(t, self.k)
            for a in bad_g:
                if a in self.live: self.mz_force(a)
        finally:
            self.guard = False

    def cx(self, c, t):
        self.before_modify(t, self.V[t] ^ self.V[c])
        self.emit('cx', (c, t)); self.V[t] ^= self.V[c]

    def rawcx(self, c, t):
        self.emit('cx', (c, t)); self.V[t] ^= self.V[c]

    def helper(self, target, excl):
        """добавить одно И(W_i, W_j) в свежую анциллу, чтобы target стал квадратично выразим."""
        ws = sorted(w for w in self.live if w not in excl)
        if not self.free: return False
        a = self.free[-1]
        for i, x in enumerate(ws):
            for y in ws[i + 1:]:
                self.V[a] = mul(self.V[x], self.V[y]); self.live.add(a)
                ok = self.quad_rep(target, set(excl)) is not None
                self.live.discard(a); self.V[a] = 0
                if ok:
                    b = self.fresh(); assert b == a
                    self.emit('AND', (x, y, a)); self.nand += 1; self.V[a] = mul(self.V[x], self.V[y])
                    return True
        return False

    def fresh(self):
        if not self.free:
            for a in sorted(self.garbage()):
                if self.safe_mz(a): break
        if not self.free: raise Fail('no ancilla')
        a = self.free.pop(); self.live.add(a); return a

    def try_AND(self, a, b):
        if not self.free:
            for x in sorted(self.garbage()):
                if self.safe_mz(x): break
        if not self.free: return None
        return self.AND(a, b)

    def AND(self, a, b):
        t = self.fresh(); self.emit('AND', (a, b, t)); self.nand += 1
        self.V[t] = mul(self.V[a], self.V[b]); return t

    def cleanup(self, k):
        fut = self.fut[k]
        for a in sorted(self.garbage()):
            if self.V[a] in fut and self.rng.random() < self.keep_prob: continue
            self.safe_mz(a)

    @staticmethod
    def decompose(pairs):
        Q = set(pairs); prods = []; lin = set()
        while Q:
            cnt = {}
            for i, j in Q: cnt[i] = cnt.get(i, 0) + 1; cnt[j] = cnt.get(j, 0) + 1
            i, j = min(Q, key=lambda p: (cnt[p[0]] + cnt[p[1]], p))
            A = {k for p in Q for k in p if i in p and j not in p and k != i}
            B = {k for p in Q for k in p if j in p and i not in p and k != j}
            for p in list(Q):
                if i in p or j in p: Q.discard(p)
            for a in A:
                for b in B:
                    if a == b: lin ^= {a}
                    else: Q ^= {(min(a, b), max(a, b))}
            prods.append(({i} | B, {j} | A))
        return prods, lin

    def realize(self, t, k):
        D = self.pend.get(t, 0)
        if not D: return
        newL = self.L[t]; ct = self.car[t]; old = self.V[ct]
        cs = self.carriers()
        for w in sorted(self.live):
            if w in cs and w != ct: continue
            if w < self.n and w != t: continue
            if self.V[w] in (newL, newL ^ ONE):
                if w != ct and ct >= self.n:
                    self.car[t] = w
                    ok = self.quad_rep(self.V[ct], {ct}) is not None
                    self.car[t] = ct
                    if not ok: continue
                if self.V[w] != newL: self.emit('x', (w,)); self.V[w] ^= ONE
                self.car[t] = w; self.pend[t] = 0; return
        lr = self.lin_rep(D, {ct}, maxw=self.linmax)
        if newL in (0, ONE) and old not in (0, ONE):
            # стирание измерением: новое значение — константа, старое квадратично выразимо через остальные
            r = self.quad_rep(old, {ct})
            if r is not None and (lr is None or len(lr[0]) > len(r[1])):
                self.before_modify(ct, 0)
                corr = [('cz', p) for p in r[1]] + [('z', (w,)) for w in r[0]] + [('x', (ct,))]
                self.emit('MZ', (ct, tuple(corr))); self.V[ct] = 0
                if newL == ONE: self.emit('x', (ct,)); self.V[ct] = ONE
                self.pend[t] = 0; return
        if lr is not None:
            for w in lr[0]: self.cx(w, ct)
            if lr[1]: self.emit('x', (ct,)); self.V[ct] ^= ONE
            self.pend[t] = 0; return
        qr = self.quad_rep(D, {ct})
        if qr is None and self.helper(D, {ct}): qr = self.quad_rep(D, {ct})
        if qr is None: raise Fail('unrep')
        lin, pairs = qr
        prods, lx = self.decompose([tuple(sorted(p)) for p in pairs])
        lin = set(lin) ^ lx
        move = old in self.fut[min(k + 1, len(self.fut) - 1)]
        if self.rng.random() < self.eps: move = not move
        if old in (0, ONE): move = False
        tgt = ct
        if not move: self.before_modify(ct, newL)
        for idx, (Ls, Ms) in enumerate(prods):
            u = min(Ls, key=lambda w: (w in Ms, w)); v = min(Ms, key=lambda w: (w in Ls, w))
            assert u not in Ms and v not in Ls
            mods = [(w, u) for w in Ls if w != u] + [(w, v) for w in Ms if w != v]
            for c_, t_ in mods: self.rawcx(c_, t_)
            if self.V[tgt] in (0, ONE) and not (idx == 0 and move):
                # цель сейчас чистая: И прямо в неё (без анциллы и без снятия)
                if self.V[tgt] == ONE: self.emit('x', (tgt,)); self.V[tgt] = 0
                self.emit('AND', (u, v, tgt)); self.nand += 1; self.V[tgt] = mul(self.V[u], self.V[v])
                for c_, t_ in reversed(mods): self.rawcx(c_, t_)
                continue
            g = self.try_AND(u, v)
            if g is None:
                self.emit('CCX', (u, v, tgt)); self.V[tgt] ^= mul(self.V[u], self.V[v])
            elif idx == 0 and move:
                self.rawcx(tgt, g); self.car[t] = g; tgt = g
            else:
                self.rawcx(g, tgt)
                if mods: self.mz(g)
            for c_, t_ in reversed(mods): self.rawcx(c_, t_)
        for w in sorted(lin): self.rawcx(w, tgt)
        if self.V[tgt] != newL:
            assert self.V[tgt] ^ ONE == newL
            self.emit('x', (tgt,)); self.V[tgt] ^= ONE
        self.pend[t] = 0

    def add(self, t, delta, k):
        lz = self.lazy if isinstance(self.lazy, bool) else (self.rng.random() < self.lazy)
        newp = self.pend.get(t, 0) ^ delta
        if not lz or (newp and self.quad_rep(newp, {self.car[t]}) is None):
            self.realize(t, k); newp = delta
        self.pend[t] = newp; self.L[t] ^= delta
        if not lz: self.realize(t, k)

    def run(self):
        for k, (g, q) in enumerate(self.lg):
            self.k = k
            if g in ('z', 'cz', 'ccz'): rd = q
            elif g in ('cx', 'ccx'): rd = q[:-1]
            else: rd = ()
            for w in rd: self.realize(w, k)
            c = self.car
            if g == 'x':
                self.emit('x', (c[q[0]],)); self.V[c[q[0]]] ^= ONE; self.L[q[0]] ^= ONE
            elif g == 'z': self.emit('z', (c[q[0]],))
            elif g == 'cz': self.emit('cz', (c[q[0]], c[q[1]]))
            elif g == 'cx': self.add(q[1], self.V[c[q[0]]], k)
            elif g == 'ccx': self.add(q[2], mul(self.V[c[q[0]]], self.V[c[q[1]]]), k)
            elif g == 'ccz':
                opts = []
                for (x, y, z) in ((q[0], q[1], q[2]), (q[0], q[2], q[1]), (q[1], q[2], q[0])):
                    r = self.lin_rep(mul(self.L[x], self.L[y]), {c[z]}, maxw=3)
                    opts.append((len(r[0]) if r else 99, self.rng.random(), x, y, z, r))
                opts.sort(key=lambda o: o[:2])
                w, _, x, y, z, r = opts[0]
                if r is not None and w <= 3:
                    for u in r[0]: self.emit('cz', (u, c[z]))
                    if r[1]: self.emit('z', (c[z],))
                else:
                    a = self.try_AND(c[x], c[y])
                    if a is None: self.emit('CCZ', (c[x], c[y], c[z]))
                    else: self.emit('cz', (a, c[z]))
            else: raise Fail(g)
            self.cleanup(k + 1)
        self.k = len(self.lg)
        for t in range(self.n): self.realize(t, self.k)
        for w in range(self.n):
            if self.car[w] != w:
                d = self.V[w] ^ self.L[w]
                r = self.lin_rep(d, {w})
                if r is not None:
                    for u in r[0]: self.cx(u, w)
                    if r[1]: self.emit('x', (w,)); self.V[w] ^= ONE
                else:
                    a = self.car[w]
                    for (x, y) in ((a, w), (w, a), (a, w)): self.rawcx(x, y)
                self.car[w] = w
        for _ in range(6):
            prog = True
            while prog and self.garbage():
                prog = False
                for a in sorted(self.garbage()):
                    if self.safe_mz(a): prog = True
            gs = [a for a in sorted(self.garbage()) if self.quad_rep(self.V[a], {a}) is None]
            if not gs: continue
            if not self.helper(self.V[gs[0]], {gs[0]}): break
        if self.garbage(): raise Fail('stuck garbage')
        for w in range(self.n):
            if self.V[w] != self.final[w]: raise Fail('final mismatch')
        return self.ops



def to_items(ops):
    it = []
    for g, q in ops:
        if g == 'AND':
            a, b, t = q
            for gg, qq in evo_lib.rccx(a, b, t) + [('sdg', (t,))]: it.append(('G', gg, tuple(qq)))
        elif g == 'MZ':
            it.append(('G', 'h', (q[0],))); it.append(('M', q[0], q[1]))
        elif g in ('CCX', 'CCZ'):
            for gg, qq in evo_lib.lower_gate(g.lower(), q): it.append(('G', gg, tuple(qq)))
        else: it.append(('G', g, tuple(q)))
    return it


FOLDABLE = {'h', 'x', 'cx', 't', 'tdg', 's', 'sdg', 'z'}


def fold_items(items, N):
    """свёртка фаз (evo_lib.fold_cancel) на унитарных кусках между измерениями/прочими элементами."""
    out, seg = [], []
    def flush():
        if seg:
            out.extend(('G', g, tuple(q)) for g, q in evo_lib.fold_cancel([(g, q) for _, g, q in seg], N))
            seg.clear()
    for it in items:
        if it[0] == 'G' and it[1] in FOLDABLE: seg.append(it)
        else: flush(); out.append(it)
    flush()
    return out


def to_qasm(ops, N, path=None, peephole=True, fold=False):
    import peep
    items = to_items(ops)
    if fold: items = fold_items(items, N)
    if peephole: items = peep.run(items)
    lines = []; nb = 0
    for it in items:
        if it[0] == 'M':
            a, corr = it[1], it[2]
            lines.append(f"c[{nb}] = measure q[{a}];")
            body = " ".join(f"{gg} " + ", ".join(f"q[{x}]" for x in qq) + ";" for gg, qq in corr)
            lines.append(f"if (c[{nb}]) {{{body}}}"); nb += 1
        else:
            lines.append(f"{it[1]} " + ", ".join(f"q[{x}]" for x in it[2]) + ";")
    txt = "OPENQASM 3;\ninclude \"stdgates.inc\";\n" + f"qubit[{N}] q;\nbit[{max(nb,1)}] c;\n" + "\n".join(lines) + "\n"
    if path: open(path, 'w').write(txt)
    return txt


def count(txt):
    T = len(re.findall(r"\b(t|tdg) q", txt)); CX = len(re.findall(r"\b(cx|cz) q", txt))
    return T, CX


def compact(ops, n):
    """перенумеровать использованные анциллы подряд -> (ops, N)."""
    used = sorted({x for g, q in ops for x in (q[:3] if g == 'AND' else ((q[0],) + tuple(y for _, qq in q[1] for y in qq)) if g == 'MZ' else q) if x >= n})
    mp = {u: n + i for i, u in enumerate(used)}
    f = lambda x: mp.get(x, x)
    out = []
    for g, q in ops:
        if g == 'MZ': out.append((g, (f(q[0]), tuple((gg, tuple(map(f, qq))) for gg, qq in q[1]))))
        else: out.append((g, tuple(map(f, q))))
    return out, n + len(used)


def inputs_idx(path, names):
    ins = None
    for line in open(path):
        if line.startswith('.i'): ins = line.split()[1:]
    return None if ins is None else [i for i, x in enumerate(names) if x in ins]


def compile_path(path, seed=0, max_anc=None, eps=0.0, keep_prob=1.0, lazy=True, linmax=3, clean=False):
    names, hl = qopt.load_qc(path); n = len(names)
    lg = logical(hl, n)
    if max_anc is None: max_anc = 22 - n
    C = Compiler(n, lg, max_anc, random.Random(seed), eps=eps, keep_prob=keep_prob, lazy=lazy, linmax=linmax,
                 inputs=inputs_idx(path, names) if clean else None)
    ops = C.run()
    ops, N = compact(ops, n)
    return ops, N, n
