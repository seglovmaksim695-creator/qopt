# SPDX-License-Identifier: Apache-2.0
"""engine: reversible-circuit builder of the exact point addition (AND + measurement erase model of Gidney 2018),
ported from the lab agents a9/k4 (arith.py + circk.py). See qoptct.crypto.ecpoint."""
from .._errors import CheckError

ZERO, ONE = ('c', 0), ('c', 1)
def W(w): return ('w', w)


class Circ:
    def __init__(self):
        self.ops = []; self.nw = 0; self.names = []; self.free = []; self.data = []
    def reg(self, name, n):
        ws = [self._new(f"{name}{i}") for i in range(n)]; self.data += ws; return ws
    def _new(self, nm):
        w = self.nw; self.nw += 1; self.names.append(nm); return w
    def alloc(self):
        return self.free.pop() if self.free else self._new(f"t{self.nw}")
    def release(self, w):
        assert w not in self.free; self.free.append(w)
    # --- вентили
    def x(self, q): self.ops.append(('x', q))
    def cx(self, c, t):
        assert c != t; self.ops.append(('cx', c, t))
    def ccx(self, a, b, t): self.ops.append(('ccx', a, b, t))
    def AND(self, a, b, t=None):
        if t is None: t = self.alloc()
        assert len({a, b, t}) == 3, (a, b, t)
        self.ops.append(('and', a, b, t)); return t
    def erase(self, t, a, b, lin=(), release=True):
        assert t not in (a, b) and t not in lin and a != b
        self.ops.append(('erase', t, a, b, tuple(lin)))
        if release: self.release(t)
    # --- подсчёт
    def counts(self):
        k = {}
        for o in self.ops: k[o[0]] = k.get(o[0], 0) + 1
        return dict(AND=k.get('and', 0), CCX=k.get('ccx', 0), ERASE=k.get('erase', 0), qubits=self.nw)


def inverse_ops(ops):
    out = []
    for o in reversed(ops):
        k = o[0]
        if k in ('x', 'cx', 'ccx'): out.append(o)
        elif k == 'and': out.append(('erase', o[3], o[1], o[2], ()))
        elif k == 'erase':
            t, a, b, lin = o[1:]
            out.append(('and', a, b, t))
            for w in lin: out.append(('cx', w, t))
        else: raise ValueError(k)
    return out


def emit_inverse(C, start):
    """заменить C.ops[start:] обратной последовательностью (функция обращается; стоимость = число стираний прямой)."""
    seg = C.ops[start:]; del C.ops[start:]; C.ops += inverse_ops(seg)


# ---------------- классическая симуляция ----------------

def simulate(C, init, S):
    M = (1 << S) - 1; v = [0] * C.nw
    for w, m in init.items(): v[w] = m
    for i, o in enumerate(C.ops):
        k = o[0]
        if k == 'x': v[o[1]] ^= M
        elif k == 'cx': v[o[2]] ^= v[o[1]]
        elif k == 'ccx': v[o[3]] ^= v[o[1]] & v[o[2]]
        elif k == 'and':
            if v[o[3]]: raise CheckError(f"op {i}: AND target not zero")
            v[o[3]] = v[o[1]] & v[o[2]]
        elif k == 'erase':
            t, a, b, lin = o[1:]; val = v[a] & v[b]
            for w in lin: val ^= v[w]
            if v[t] != val: raise CheckError(f"op {i}: erase condition violated")
            v[t] = 0
        else: raise ValueError(k)
    return v

def pack(vals, wires):
    m = [0] * len(wires)
    for s, x in enumerate(vals):
        for i in range(len(wires)):
            if (x >> i) & 1: m[i] |= 1 << s
    return dict(zip(wires, m))

def unpack(v, wires, S):
    out = [0] * S
    for i, w in enumerate(wires):
        m = v[w]
        for s in range(S):
            if (m >> s) & 1: out[s] |= 1 << i
    return out


# ---------------- перенос (MAJ) с распространением констант ----------------
class Rec:
    """рецепт бита переноса MAJ(x,y,c): kind in const/alias/and/or/gid."""
    def __init__(self, kind, *args): self.kind = kind; self.args = args; self.t = None; self.fold = False

def plan(x, y, c):
    ops = [x, y, c]
    for i in range(3):
        for j in range(i + 1, 3):
            if ops[i] == ops[j]: return Rec('op', ops[i])
    cs = [o for o in ops if o[0] == 'c']
    if len(cs) >= 2:   # 0 и 1 -> третий
        return Rec('op', [o for o in ops if o[0] == 'w'][0] if len(cs) == 2 else ZERO)
    if len(cs) == 1:
        u, v = [o[1] for o in ops if o[0] == 'w']
        return Rec('and' if cs[0][1] == 0 else 'or', u, v)
    return Rec('gid', x[1], y[1], c[1])   # x — провод-цель (можно «сложить» с c)

def compute(C, r, fold_ok=False, t=None):
    """вычислить перенос; вернуть операнд. t — заданный нулевой провод (иначе новый)."""
    k = r.kind
    if k == 'op':
        o = r.args[0]
        if t is None: return o
        if o == ONE: C.x(t)
        elif o[0] == 'w': C.cx(o[1], t)
        r.t = t; return W(t)
    if k == 'and':
        u, v = r.args; r.t = C.AND(u, v, t); return W(r.t)
    if k == 'or':
        u, v = r.args; C.x(u); C.x(v); r.t = C.AND(u, v, t); C.x(u); C.x(v); C.x(r.t); return W(r.t)
    x, y, c = r.args
    C.cx(c, x); C.cx(c, y); r.t = C.AND(x, y, t); C.cx(c, y)
    if fold_ok: r.fold = True
    else: C.cx(c, x)
    C.cx(c, r.t); return W(r.t)

def uncompute(C, r, t=None, release=True):
    """стереть значение переноса с провода t (по умолчанию r.t). Для 'op' — CX/X (детерминированно)."""
    k = r.kind; t = r.t if t is None else t
    if k == 'op':
        o = r.args[0]
        if o == ONE: C.x(t)
        elif o[0] == 'w': C.cx(o[1], t)
        if release: C.release(t)
        return
    if k == 'and':
        u, v = r.args; C.erase(t, u, v, release=release); return
    if k == 'or':
        u, v = r.args; C.x(t); C.x(u); C.x(v); C.erase(t, u, v, release=release); C.x(u); C.x(v); return
    x, y, c = r.args
    if not r.fold: C.cx(c, x)
    C.cx(c, y); C.erase(t, x, y, lin=(c,), release=release); C.cx(c, y)
    if not r.fold: C.cx(c, x)


# ---------------- блоки ----------------
def add(C, a, b, cout=None):
    """b += a (mod 2^n); если cout (нулевой провод) задан — перенос наружу в cout. И: n-1 (без cout) / n (с cout). Гидни 2018."""
    n = len(a); c = [None] * (n + 1)
    last = n if cout is not None else n - 1
    for i in range(last):
        t = cout if i == n - 1 else C.alloc()
        if c[i] is None: C.AND(a[i], b[i], t)
        else:
            C.cx(c[i], a[i]); C.cx(c[i], b[i]); C.AND(a[i], b[i], t); C.cx(c[i], t)
        c[i + 1] = t
    for i in reversed(range(n)):
        if i < last:
            if i + 1 < n:
                C.erase(c[i + 1], a[i], b[i], lin=(() if c[i] is None else (c[i],)))
            if c[i] is not None: C.cx(c[i], a[i])
            C.cx(a[i], b[i])
        else:
            C.cx(a[i], b[i])
            if c[i] is not None: C.cx(c[i], b[i])

def sub(C, a, b, bout=None):
    """b -= a (mod 2^n); bout ^= [b < a] (заём)."""
    for w in b: C.x(w)
    add(C, a, b, bout)
    for w in b: C.x(w)

def add_k(C, k, b, ctrl=None, known=None):
    """b += k (или ctrl·k), k — классическая константа, по модулю 2^len(b). known: {провод: константа} на входе."""
    m = len(b); known = known or {}
    X = [('c', known[w]) if w in known else W(w) for w in b]
    Y = [(W(ctrl) if ctrl is not None else ONE) if (k >> i) & 1 else ZERO for i in range(m)]
    c = [ZERO]; R = [None]
    for i in range(m - 1):
        r = plan(X[i], Y[i], c[i]); c.append(compute(C, r, fold_ok=True)); R.append(r)
    for i in reversed(range(m)):
        if i + 1 < m and R[i + 1].kind != 'op': uncompute(C, R[i + 1])
        if Y[i] == ONE: C.x(b[i])
        elif Y[i][0] == 'w': C.cx(Y[i][1], b[i])
        if i > 0 and R[i].kind == 'gid' and R[i].fold:
            pass  # unreachable: fold относится к биту i (цель x_i) в R[i+1]
        folded = (i + 1 < m) and R[i + 1].kind == 'gid' and R[i + 1].fold
        if not folded:
            if c[i] == ONE: C.x(b[i])
            elif c[i][0] == 'w': C.cx(c[i][1], b[i])

def carry_chain(C, X, Y, cin=ZERO):
    """переносы c_1..c_{m-1} для X+Y+cin (без сумм); вернуть (c, R) — c[m-1] последний вычисленный, R для стирания."""
    m = len(X); c = [cin]; R = [None]
    for i in range(m - 1):
        r = plan(X[i], Y[i], c[i]); c.append(compute(C, r, fold_ok=False)); R.append(r)
    return c, R

def chain_undo(C, R):
    for r in reversed(R[1:]):
        if r.kind != 'op': uncompute(C, r)

def clear_cout(C, X, Y, f, cin=ZERO, release=True):
    """f обязан равняться переносу наружу X+Y+cin; стереть f (измерением; И: m-1 минус свободные)."""
    m = len(X); c, R = carry_chain(C, X, Y, cin)
    rf = plan(X[m - 1], Y[m - 1], c[m - 1]); uncompute(C, rf, t=f, release=release)
    chain_undo(C, R)

def compute_cout(C, X, Y, f, cin=ZERO):
    """f (нулевой) := перенос наружу X+Y+cin. И: m (минус свободные)."""
    m = len(X); c, R = carry_chain(C, X, Y, cin)
    rf = plan(X[m - 1], Y[m - 1], c[m - 1]); compute(C, rf, t=f)
    chain_undo(C, R)


# ---------------- модульная арифметика (p нечётно, p < 2^n) ----------------
def mod_add(C, a, b, p):
    """b <- a + b mod p (a, b < p). И = 4n-3."""
    n = len(a); z = C.alloc()
    add(C, a, b, cout=z)                                # n
    add_k(C, (1 << (n + 1)) - p, b + [z])               # n-1 : b:z -= p, z = [a+b < p]
    add_k(C, p, b, ctrl=z)                              # n-1 : если z, вернуть p
    C.x(z)                                              # z = [b' < a] = перенос наружу (~b') + a
    for w in b: C.x(w)
    clear_cout(C, [W(w) for w in b], [W(w) for w in a], z)   # n-1
    for w in b: C.x(w)

def mod_sub(C, a, b, p):
    """b <- b - a mod p = обращение mod_add. И = 4n-3."""
    s = len(C.ops); mod_add(C, a, b, p); emit_inverse(C, s)

def ctrl_mod_add(C, q, a, b, p, sub=False):
    """b <- b ± q·a mod p. И = n + (4n-3)."""
    t = [C.AND(q, w) for w in a]
    (mod_sub if sub else mod_add)(C, t, b, p)
    for w, tw in zip(reversed(a), reversed(t)): C.erase(tw, q, w)

def mod_double(C, b, p):
    """b <- 2b mod p, на месте с переименованием: возвращает новый список проводов результата
    (старший провод b[n-1] освобождается как нулевая анцилла)."""
    n = len(b); lo = C.alloc()
    y = [lo] + list(b)
    add_k(C, (1 << (n + 1)) - p, y, known={lo: 0})      # 2b - p, знак в b[n-1]
    z = b[n - 1]; r = [lo] + list(b[:n - 1])
    add_k(C, p, r, ctrl=z, known={lo: 1})
    C.cx(lo, z); C.x(z)                                 # z = 1 - lsb(r) -> 0
    C.release(z)
    return r


# ---------------- compact op storage (k4/circk) ----------------
from array import array

KIND = {'x': 0, 'cx': 1, 'ccx': 2, 'and': 3, 'erase': 4}
NAME = ['x', 'cx', 'ccx', 'and', 'erase']


def enc(o):
    k = KIND[o[0]]
    if k == 4:
        lin = o[4]; assert len(lin) <= 1
        return (4, o[1], o[2], o[3], lin[0] if lin else -1)
    r = list(o[1:]) + [-1] * (4 - (len(o) - 1))
    return (k, *r)


def dec(t):
    k = t[0]
    if k == 0: return ('x', t[1])
    if k == 1: return ('cx', t[1], t[2])
    if k in (2, 3): return (NAME[k], t[1], t[2], t[3])
    return ('erase', t[1], t[2], t[3], () if t[4] < 0 else (t[4],))


class OpStore:
    def __init__(self): self.a = array('i')
    def __len__(self): return len(self.a) // 5
    def append(self, o): self.a.extend(enc(o))
    def __iadd__(self, ops):
        for o in ops: self.append(o)
        return self
    def __iter__(self):
        a = self.a
        for i in range(0, len(a), 5): yield dec(a[i:i + 5])
    def __getitem__(self, s):
        if isinstance(s, slice):
            st, sp, _ = s.indices(len(self))
            return [dec(self.a[5 * i:5 * i + 5]) for i in range(st, sp)]
        return dec(self.a[5 * s:5 * s + 5])
    def __delitem__(self, s):
        assert isinstance(s, slice) and s.stop is None
        del self.a[5 * s.start:]
    def append_inverse(self, s, e):
        """дописать обращение сегмента [s,e) потоково."""
        a = self.a
        for i in range(e - 1, s - 1, -1):
            k = a[5 * i]
            if k in (0, 1, 2): a.extend(a[5 * i:5 * i + 5])
            elif k == 3: a.extend((4, a[5 * i + 3], a[5 * i + 1], a[5 * i + 2], -1))
            else:
                t, x, y, l = a[5 * i + 1], a[5 * i + 2], a[5 * i + 3], a[5 * i + 4]
                a.extend((3, x, y, t, -1))
                if l >= 0: a.extend((1, l, t, -1, -1))
    def counts(self):
        a = self.a; c = [0] * 5
        for i in range(0, len(a), 5): c[a[i]] += 1
        return c


class CircK(Circ):
    def __init__(self):
        super().__init__(); self.ops = OpStore()
    def counts(self):
        c = self.ops.counts()
        return dict(AND=c[3], CCX=c[2], ERASE=c[4], CXops=c[1], qubits=self.nw,
                    T=4 * c[3] + 7 * c[2], CX=c[1] + 3 * c[3] + c[4] + 6 * c[2])


def compute_apply_uncompute(C, fn, apply):
    s = len(C.ops); res, live = fn(); e = len(C.ops)
    apply(res)
    if isinstance(C.ops, OpStore): C.ops.append_inverse(s, e)
    else: C.ops += inverse_ops(C.ops[s:e])
    for w in reversed(live): C.release(w)
