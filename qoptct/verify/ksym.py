# SPDX-License-Identifier: Apache-2.0
"""ksym — word-level specifications for kprove: unsigned integers whose bits are BDD functions of the circuit inputs.

    def spec(S):                         # S: context given by kprove
        a = S.reg(range(0, 256))         # input wires 0..255, LSB first
        b = S.reg(range(256, 512))
        s = (a + b).mod_reduce(N)        # or: t = a + b; s = S.mux(t >= N, t - N, t)
        return {**S.out(range(256, 512), s)}   # expected output bits of wires 256..511; other data wires: unchanged
                                                # (clean wires: 0) unless given
    kprove(qasm, spec=spec, n=513, clean_wires=[512], pre=...)
Every operation is a plain Boolean construction (ripple adders / borrow comparators / muxes) on canonical BDDs;
the library is self-tested exhaustively at small widths (python ksym.py).
"""
from .bdd import ONE, ZERO


class SymInt:
    def __init__(self, B, bits):
        self.B = B
        self.bits = list(bits)

    # ------------------------------------------------------------ construction
    @staticmethod
    def const(B, value, width):
        return SymInt(B, [ONE if (value >> i) & 1 else ZERO for i in range(width)])

    def width(self):
        return len(self.bits)

    def resize(self, w):
        b = self.bits[:w] + [ZERO] * max(0, w - len(self.bits))
        return SymInt(self.B, b)

    def _lift(self, o, w=None):
        if isinstance(o, int):
            return SymInt.const(self.B, o, w or max(self.width(), o.bit_length()))
        return o

    # ------------------------------------------------------------ arithmetic (results keep the carry: width+1)
    def add(self, o, cin=ZERO):
        o = self._lift(o)
        B = self.B
        w = max(self.width(), o.width())
        a, b = self.resize(w).bits, o.resize(w).bits
        c = cin
        out = []
        for i in range(w):
            out.append(B.XOR(B.XOR(a[i], b[i]), c))
            c = B.MAJ(a[i], b[i], c)
        return SymInt(B, out + [c])

    def __add__(self, o):
        return self.add(o)

    def sub(self, o):
        """(self - o) mod 2^w with w = max width; plus borrow flag as top bit (1 iff self < o)."""
        o = self._lift(o)
        w = max(self.width(), o.width())
        nb = SymInt(self.B, [x ^ 1 for x in o.resize(w).bits])
        r = self.resize(w).add(nb, cin=ONE)
        return SymInt(self.B, r.bits[:w] + [r.bits[w] ^ 1])

    def __sub__(self, o):
        return self.sub(o)

    def lt(self, o):
        """BDD of self < o (unsigned)."""
        return self.sub(o).bits[-1]

    def __lt__(self, o):
        return self.lt(o)

    def __ge__(self, o):
        return self.lt(o) ^ 1

    def eq(self, o):
        o = self._lift(o)
        B = self.B
        w = max(self.width(), o.width())
        r = ONE
        for x, y in zip(self.resize(w).bits, o.resize(w).bits):
            r = B.AND(r, B.XOR(x, y) ^ 1)
        return r

    def mod_reduce(self, N, width=None):
        """self mod N assuming self < 2N (one conditional subtraction)."""
        w = width or max(self.width(), N.bit_length())
        t = self.resize(w + 1)
        d = t.sub(N)
        ge = d.bits[-1] ^ 1
        return mux(self.B, ge, d.resize(w), t.resize(w))

    def shl(self, k, width=None):
        b = [ZERO] * k + self.bits
        return SymInt(self.B, b[:width] if width else b)

    def shr(self, k):
        return SymInt(self.B, self.bits[k:])

    def __xor__(self, o):
        o = self._lift(o)
        w = max(self.width(), o.width())
        return SymInt(self.B, [self.B.XOR(x, y) for x, y in zip(self.resize(w).bits, o.resize(w).bits)])

    def __and__(self, o):
        o = self._lift(o)
        w = max(self.width(), o.width())
        return SymInt(self.B, [self.B.AND(x, y) for x, y in zip(self.resize(w).bits, o.resize(w).bits)])

    def __invert__(self):
        return SymInt(self.B, [x ^ 1 for x in self.bits])

    def neg_mod(self, N):
        """(N - self) mod N, for self < N."""
        w = max(self.width(), N.bit_length())
        z = self.eq(0)
        d = SymInt.const(self.B, N, w).sub(self).resize(w)
        return mux(self.B, z, SymInt.const(self.B, 0, w), d)

    def mul_const(self, c, width=None):
        """self * c (shift-add); BDDs of products can grow exponentially - use for small widths only."""
        w = width or self.width() + c.bit_length()
        acc = SymInt.const(self.B, 0, w)
        i = 0
        while c >> i:
            if (c >> i) & 1:
                acc = acc.add(self.shl(i, w)).resize(w)
            i += 1
        return acc


def mux(B, c, a, b):
    """c ? a : b, bitwise."""
    w = max(a.width(), b.width())
    return SymInt(B, [B.ITE(c, x, y) for x, y in zip(a.resize(w).bits, b.resize(w).bits)])


def lookup(B, addr, table, width):
    """table[addr] (addr: SymInt) as a mux tree."""
    vals = [SymInt.const(B, v, width) for v in table]
    for bit in addr.bits:
        if len(vals) == 1:
            break
        vals = [mux(B, bit, vals[i + 1] if i + 1 < len(vals) else SymInt.const(B, 0, width), vals[i])
                for i in range(0, len(vals), 2)]
    return vals[0]


class SpecContext:
    """Passed to spec(S). S.x(q): input literal of wire q (ZERO for clean wires). S.reg(wires) -> SymInt (LSB first)."""

    def __init__(self, B, xedge, clean_wires):
        self.B = B
        self._x = xedge
        self.clean = set(clean_wires)

    def x(self, q):
        return ZERO if q in self.clean else self._x(q)

    def reg(self, wires):
        return SymInt(self.B, [self.x(q) for q in wires])

    def const(self, v, w):
        return SymInt.const(self.B, v, w)

    def mux(self, c, a, b):
        return mux(self.B, c, a, b)

    def lookup(self, addr, table, width):
        return lookup(self.B, addr, table, width)

    def out(self, wires, val):
        """{wire: bit} for an output register (val truncated / zero-extended to len(wires))."""
        wires = list(wires)
        v = val.resize(len(wires))
        return dict(zip(wires, v.bits))


# ------------------------------------------------------------------ exhaustive self-test
if __name__ == '__main__':
    import itertools
    from .bdd import BDD
    for w in (1, 2, 3, 4):
        B = BDD()
        va = [B.new_var()[1] for _ in range(w)]
        vb = [B.new_var()[1] for _ in range(w)]
        lv = list(range(2 * w))
        A, Bv = SymInt(B, va), SymInt(B, vb)
        for N in range(1, 1 << w):
            exprs = {'add': (A + Bv, lambda a, b: a + b), 'sub': (A - Bv, lambda a, b: ((a - b) % (1 << w)) | ((a < b) << w)),
                     'modadd': ((A + Bv).mod_reduce(N, w), lambda a, b: (a + b) % N if a < N and b < N else None),
                     'neg': (A.neg_mod(N), lambda a, b: (-a) % N if a < N else None),
                     'mulc': (A.mul_const(N, 2 * w), lambda a, b: a * N),
                     'lut': (lookup(B, A, [(7 * i + N) % 16 for i in range(1 << w)], 4), lambda a, b: (7 * a + N) % 16)}
            for nm, (s, f) in exprs.items():
                for a, b in itertools.product(range(1 << w), repeat=2):
                    want = f(a, b)
                    if want is None:
                        continue
                    asg = {i: (a >> i) & 1 for i in range(w)}
                    asg.update({w + i: (b >> i) & 1 for i in range(w)})
                    got = sum(int(B.eval(bit, asg)) << i for i, bit in enumerate(s.bits))
                    assert got == want, (nm, w, N, a, b, got, want)
            for a, b in itertools.product(range(1 << w), repeat=2):
                asg = {i: (a >> i) & 1 for i in range(w)}
                asg.update({w + i: (b >> i) & 1 for i in range(w)})
                assert B.eval(A.lt(Bv), asg) == (a < b) and B.eval(A.eq(Bv), asg) == (a == b)
    print('ksym self-test ok')
