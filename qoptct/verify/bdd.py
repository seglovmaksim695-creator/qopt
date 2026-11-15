# SPDX-License-Identifier: Apache-2.0
"""Small ROBDD package with complement edges (pure Python), for kprove.

Edge = int: (node_index << 1) | complement. ONE = 0 (terminal, regular), ZERO = 1.
Variable level = creation index; a LARGER level is closer to the root (newest variables on top).
Canonical form: the hi (then) edge of every node is regular; equal functions <=> equal edges.
"""
import sys
from array import array

sys.setrecursionlimit(1 << 20)

ONE = 0
ZERO = 1


class BDD:
    def __init__(self):
        # compact storage: typed arrays (8 bytes per field) and single-int keys (memory ~1/3 of tuples/lists)
        self.var = array('q', [-1])
        self.hi = array('q', [0])
        self.lo = array('q', [0])
        self.unique = {}
        self.nvars = 0
        self.c_and = {}
        self.c_xor = {}
        self.hook = None

    # ------------------------------------------------------------------ basics
    def nodes(self):
        return len(self.var)

    def level(self, e):
        return self.var[e >> 1]

    def new_var(self):
        v = self.nvars
        self.nvars += 1
        return v, self.mk(v, ONE, ZERO)

    def lit(self, v):
        return self.mk(v, ONE, ZERO)

    def mk(self, v, hi, lo):
        if hi == lo:
            return hi
        c = hi & 1
        if c:
            hi ^= 1
            lo ^= 1
        key = (v << 80) | (hi << 40) | lo
        n = self.unique.get(key)
        if n is None:
            n = len(self.var)
            if (n & 0xFFFF) == 0 and self.hook is not None:
                self.hook()           # memory / time guard of the caller (may raise)
            self.var.append(v)
            self.hi.append(hi)
            self.lo.append(lo)
            self.unique[key] = n
        return (n << 1) | c

    def cof(self, e, v):
        n = e >> 1
        if self.var[n] != v:
            return e, e
        c = e & 1
        return self.hi[n] ^ c, self.lo[n] ^ c

    # ------------------------------------------------------------------ boolean ops
    def AND(self, f, g):
        if f == g:
            return f
        if f ^ g == 1:
            return ZERO
        if f == ONE:
            return g
        if g == ONE:
            return f
        if f == ZERO or g == ZERO:
            return ZERO
        if f > g:
            f, g = g, f
        key = (f << 40) | g
        r = self.c_and.get(key)
        if r is not None:
            return r
        var = self.var
        vf, vg = var[f >> 1], var[g >> 1]
        v = vf if vf > vg else vg
        f1, f0 = self.cof(f, v)
        g1, g0 = self.cof(g, v)
        r = self.mk(v, self.AND(f1, g1), self.AND(f0, g0))
        self.c_and[key] = r
        return r

    def OR(self, f, g):
        return self.AND(f ^ 1, g ^ 1) ^ 1

    def XOR(self, f, g):
        c = (f ^ g) & 1
        f &= ~1
        g &= ~1
        if f == g:
            return ZERO ^ c
        if f == ONE:
            return g ^ 1 ^ c
        if g == ONE:
            return f ^ 1 ^ c
        if f > g:
            f, g = g, f
        key = (f << 40) | g
        r = self.c_xor.get(key)
        if r is None:
            var = self.var
            vf, vg = var[f >> 1], var[g >> 1]
            v = vf if vf > vg else vg
            f1, f0 = self.cof(f, v)
            g1, g0 = self.cof(g, v)
            r = self.mk(v, self.XOR(f1, g1), self.XOR(f0, g0))
            self.c_xor[key] = r
        return r ^ c

    def ITE(self, c, t, e):
        return self.OR(self.AND(c, t), self.AND(c ^ 1, e))

    def MAJ(self, a, b, c):
        return self.OR(self.AND(a, b), self.AND(c, self.XOR(a, b)))

    def restrict(self, f, v, val, memo=None):
        if memo is None:
            memo = {}
        return self._restrict(f, v, val, memo)

    def _restrict(self, f, v, val, memo):
        n = f >> 1
        lv = self.var[n]
        if lv < v:
            return f
        c = f & 1
        if lv == v:
            return (self.hi[n] if val else self.lo[n]) ^ c
        r = memo.get(n)
        if r is None:
            r = self.mk(lv, self._restrict(self.hi[n], v, val, memo), self._restrict(self.lo[n], v, val, memo))
            memo[n] = r
        return r ^ c

    def exists(self, f, levels, memo=None):
        """exists over the variables in `levels` (a set)."""
        if memo is None:
            memo = {}
        n = f >> 1
        if n == 0:
            return f
        if f in memo:
            return memo[f]
        c = f & 1
        lv = self.var[n]
        h = self.exists(self.hi[n] ^ c, levels, memo)
        l = self.exists(self.lo[n] ^ c, levels, memo)
        r = self.OR(h, l) if lv in levels else self.mk(lv, h, l)
        memo[f] = r
        return r

    def compose(self, f, v, g):
        """f with variable v replaced by function g."""
        if self.var[f >> 1] < v:
            return f
        f1 = self.restrict(f, v, 1)
        f0 = self.restrict(f, v, 0)
        if f1 == f0:
            return f1
        return self.ITE(g, f1, f0)

    def support(self, f, minlevel=-1):
        seen = set()
        out = set()
        stack = [f >> 1]
        var, hi, lo = self.var, self.hi, self.lo
        while stack:
            n = stack.pop()
            if n in seen or var[n] < minlevel or n == 0:
                continue
            seen.add(n)
            out.add(var[n])
            stack.append(hi[n] >> 1)
            stack.append(lo[n] >> 1)
        return out

    def size(self, f):
        seen = set()
        stack = [f >> 1]
        while stack:
            n = stack.pop()
            if n in seen or n == 0:
                continue
            seen.add(n)
            stack.append(self.hi[n] >> 1)
            stack.append(self.lo[n] >> 1)
        return len(seen)

    def pick_sat(self, f, prefer=0):
        """One satisfying assignment {level: bit} of f (only variables on the chosen path), or None."""
        if f == ZERO:
            return None
        asg = {}
        while f != ONE:
            n = f >> 1
            c = f & 1
            h, l = self.hi[n] ^ c, self.lo[n] ^ c
            if (l != ZERO and prefer == 0) or h == ZERO:
                asg[self.var[n]] = 0
                f = l
            else:
                asg[self.var[n]] = 1
                f = h
        return asg

    def eval(self, f, asg):
        while (f >> 1) != 0:
            n = f >> 1
            c = f & 1
            f = (self.hi[n] if asg.get(self.var[n], 0) else self.lo[n]) ^ c
        return f == ONE

    def satcount(self, f, levels):
        """Number of assignments of the variables `levels` (f must depend only on them) with f = 1."""
        lv = sorted(levels)
        pos = {v: i for i, v in enumerate(lv)}
        memo = {}

        def rec(e):
            # returns count over variables with position < pos(top)... use: count over vars below level of e
            n = e >> 1
            if n == 0:
                return (1 if e == ONE else 0), -1
            if e in memo:
                return memo[e], pos[self.var[n]]
            c = e & 1
            p = pos[self.var[n]]
            tot = 0
            for ch in (self.hi[n] ^ c, self.lo[n] ^ c):
                k, pc = rec(ch)
                tot += k << (p - 1 - pc)
            memo[e] = tot
            return tot, p
        k, p = rec(f)
        return k << (len(lv) - 1 - p)

    # ------------------------------------------------------------------ garbage collection
    def gc(self, roots):
        """Keep only nodes reachable from roots (list of edges). Returns a function mapping old edges to new edges."""
        live = bytearray(len(self.var))
        live[0] = 1
        stack = [r >> 1 for r in roots]
        hi, lo = self.hi, self.lo
        while stack:
            n = stack.pop()
            if live[n]:
                continue
            live[n] = 1
            stack.append(hi[n] >> 1)
            stack.append(lo[n] >> 1)
        remap = array('q', bytes(8 * len(self.var)))
        nv, nh, nl = array('q', [-1]), array('q', [0]), array('q', [0])
        uniq = {}
        for n in range(1, len(self.var)):
            if live[n]:
                h = (remap[hi[n] >> 1] << 1) | (hi[n] & 1)
                l = (remap[lo[n] >> 1] << 1) | (lo[n] & 1)
                remap[n] = len(nv)
                nv.append(self.var[n])
                nh.append(h)
                nl.append(l)
                uniq[(self.var[n] << 80) | (h << 40) | l] = remap[n]
        self.var, self.hi, self.lo, self.unique = nv, nh, nl, uniq
        self.c_and = {}
        self.c_xor = {}
        return lambda e: (remap[e >> 1] << 1) | (e & 1)


class Z8:
    """Z8-valued function as three bit-BDDs (b0, b1, b2)."""

    @staticmethod
    def zero():
        return (ZERO, ZERO, ZERO)

    @staticmethod
    def add_term(B, acc, c, f):
        c &= 7
        if not c or f == ZERO:
            return acc
        t = (f if c & 1 else ZERO, f if c & 2 else ZERO, f if c & 4 else ZERO)
        a0, a1, a2 = acc
        s0 = B.XOR(a0, t[0])
        k0 = B.AND(a0, t[0])
        s1 = B.XOR(B.XOR(a1, t[1]), k0)
        k1 = B.MAJ(a1, t[1], k0)
        s2 = B.XOR(B.XOR(a2, t[2]), k1)
        return (s0, s1, s2)
