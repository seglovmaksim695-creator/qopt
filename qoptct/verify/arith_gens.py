# SPDX-License-Identifier: Apache-2.0
"""Generators of large arithmetic circuits (OpenQASM 3, qoptct dialect) and their classical references, used to
exercise qoptct.verify.arith at 64-2048 bits (ported from the lab, agent a8)."""
import numpy as np


class QB:
    def __init__(s, N):
        s.N = N; s.L = []; s.nb = 0

    def g(s, name, *qs):
        s.L.append("%s %s;" % (name, ", ".join("q[%d]" % q for q in qs)))

    def and_compute(s, a, b, t):
        """logical AND (Gidney 2018, fig. 3): t (|0>) := a AND b"""
        s.g("h", t); s.g("t", t); s.g("cx", a, t); s.g("cx", b, t); s.g("cx", t, a); s.g("cx", t, b)
        s.g("tdg", a); s.g("tdg", b); s.g("t", t); s.g("cx", t, a); s.g("cx", t, b); s.g("h", t); s.g("s", t)

    def and_uncompute(s, a, b, t):
        k = s.nb; s.nb += 1
        s.g("h", t); s.L.append("c[%d] = measure q[%d];" % (k, t))
        s.L.append("if (c[%d]) {cz q[%d], q[%d]; x q[%d];}" % (k, a, b, t))

    def text(s):
        return "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[%d] q;\nbit[%d] c;\n" % (s.N, max(s.nb, 1)) + "\n".join(s.L) + "\n"


def gidney_adder(n):
    """a: 0..n-1, b: n..2n-1, z: 2n; ancillas t_1..t_n = 2n+1..3n. b := a+b mod 2^n, z ^= carry."""
    A = lambda i: i; Bq = lambda i: n + i; Z = 2 * n; T = lambda i: 2 * n + i   # T(i), i=1..n: carry c_i
    q = QB(3 * n + 1)
    for i in range(n):
        if i >= 1:
            q.g("cx", T(i), A(i)); q.g("cx", T(i), Bq(i))
        q.and_compute(A(i), Bq(i), T(i + 1))
        if i >= 1: q.g("cx", T(i), T(i + 1))
    q.g("cx", T(n), Z)
    for i in range(n - 1, -1, -1):
        if i >= 1: q.g("cx", T(i), T(i + 1))
        q.and_uncompute(A(i), Bq(i), T(i + 1))
        q.g("cx", A(i), Bq(i))
        if i >= 1:
            q.g("cx", T(i), Bq(i)); q.g("cx", T(i), A(i))
    return q.text()


def gidney_adder_ref(n):
    def f(X):
        D = X.copy()
        for r in range(len(X)):
            a = int("".join(str(int(v)) for v in X[r, :n][::-1]), 2)
            b = int("".join(str(int(v)) for v in X[r, n:2 * n][::-1]), 2)
            s = a + b
            for i in range(n): D[r, n + i] = (s >> i) & 1
            D[r, 2 * n] = X[r, 2 * n] ^ ((s >> n) & 1)
        return D
    return (2 * n + 1, f, list(range(2 * n + 1)))


def comparator(n):
    """z ^= carry c_n(a, b) (comparison a > ~b); a and b are restored. a:0..n-1, b:n..2n-1, z:2n"""
    A = lambda i: i; Bq = lambda i: n + i; Z = 2 * n; T = lambda i: 2 * n + i
    q = QB(3 * n + 1)
    for i in range(n):
        if i >= 1:
            q.g("cx", T(i), A(i)); q.g("cx", T(i), Bq(i))
        q.and_compute(A(i), Bq(i), T(i + 1))
        if i >= 1: q.g("cx", T(i), T(i + 1))
    q.g("cx", T(n), Z)
    for i in range(n - 1, -1, -1):
        if i >= 1: q.g("cx", T(i), T(i + 1))
        q.and_uncompute(A(i), Bq(i), T(i + 1))
        if i >= 1:
            q.g("cx", T(i), Bq(i)); q.g("cx", T(i), A(i))
    return q.text()


def comparator_ref(n):
    def f(X):
        D = X.copy()
        for r in range(len(X)):
            a = int("".join(str(int(v)) for v in X[r, :n][::-1]), 2)
            b = int("".join(str(int(v)) for v in X[r, n:2 * n][::-1]), 2)
            D[r, 2 * n] = X[r, 2 * n] ^ ((a + b) >> n & 1)
        return D
    return (2 * n + 1, f, list(range(2 * n + 1)))


def const_adder(n, cst):
    """b := b + cst mod 2^n (cst a classical constant), z ^= carry. b:0..n-1, z:n. Carries via AND (bit cst=1: c' = b|c)."""
    Bq = lambda i: i; Z = n; T = lambda i: n + i       # T(i): carry c_i, i=1..n
    q = QB(2 * n + 1)
    cb = [(cst >> i) & 1 for i in range(n)]
    # c_{i+1} = maj(cb_i, b_i, c_i): cb=0: b_i & c_i ; cb=1: b_i | c_i = b_i ^ c_i ^ b_i c_i
    for i in range(n):
        if i == 0:
            if cb[0]: q.g("cx", Bq(0), T(1))              # c_1 = b_0 (cb=1, c_0=0)
            # cb=0: c_1 = 0 (T(1) stays |0>)
        else:
            if cb[i] == 0:
                q.and_compute(Bq(i), T(i), T(i + 1))
            else:
                # c' = b ^ c ^ (b c): AND into the carry wire, then the two xors
                q.and_compute(Bq(i), T(i), T(i + 1)); q.g("cx", Bq(i), T(i + 1)); q.g("cx", T(i), T(i + 1))
    q.g("cx", T(n), Z)
    for i in range(n - 1, -1, -1):
        if i == 0:
            if cb[0]: q.g("cx", Bq(0), T(1))
        else:
            if cb[i] == 1: q.g("cx", T(i), T(i + 1)); q.g("cx", Bq(i), T(i + 1))
            q.and_uncompute(Bq(i), T(i), T(i + 1))
        # sum: b_i ^= cb_i ^ c_i, applied after carry i+1 is uncomputed (c_i is still present)
        if cb[i]: q.g("x", Bq(i))
        if i >= 1: q.g("cx", T(i), Bq(i))
    return q.text()


def const_adder_ref(n, cst):
    def f(X):
        D = X.copy()
        for r in range(len(X)):
            b = int("".join(str(int(v)) for v in X[r, :n][::-1]), 2)
            s = b + cst
            for i in range(n): D[r, i] = (s >> i) & 1
            D[r, n] = X[r, n] ^ ((s >> n) & 1)
        return D
    return (n + 1, f, list(range(n + 1)))


def schoolbook_mult(n):
    """not implemented (see multiplier())."""
    raise NotImplementedError


def build(name, n, *extra):
    if name == "gidney_adder": return gidney_adder(n), gidney_adder_ref(n)
    if name == "comparator": return comparator(n), comparator_ref(n)
    if name == "const_adder": return const_adder(n, extra[0]), const_adder_ref(n, extra[0])
    raise KeyError(name)


def adder_qc(n, sum_bits=True):
    """.qc reference (Toffoli) for gidney_adder / comparator: wires a0.. b0.. z c1..cn (c are clean ancillas, not on .i)"""
    nm = ["a%d" % i for i in range(n)] + ["b%d" % i for i in range(n)] + ["z"] + ["c%d" % i for i in range(1, n + 1)]
    L = [".v " + " ".join(nm), ".i " + " ".join(nm[:2 * n + 1]), "", "BEGIN"]
    def carry(i):  # c_{i+1} ^= maj(a_i, b_i, c_i)
        a, b, c1 = "a%d" % i, "b%d" % i, "c%d" % (i + 1)
        L.append("tof %s %s %s" % (a, b, c1))
        if i >= 1:
            ci = "c%d" % i
            L.append("tof %s %s %s" % (a, ci, c1)); L.append("tof %s %s %s" % (b, ci, c1))
    for i in range(n): carry(i)
    L.append("tof c%d z" % n)
    for i in range(n - 1, -1, -1):
        carry(i)
        if sum_bits:
            L.append("tof a%d b%d" % (i, i))
            if i >= 1: L.append("tof c%d b%d" % (i, i))
    L.append("END")
    return "\n".join(L) + "\n"


def lookup(k, tab, w):
    """Table lookup (unary iteration, Babbush et al. 2018 + Gidney logical AND): t ^= tab[addr].
    addr: 0..k-1 (addr[0] most significant), t: k..k+w-1, ancillas k+w..k+w+k-2 (level l>=1 -> k+w+l-1)."""
    N = k + w + max(k - 1, 0)
    q = QB(N); A = lambda l: l; Tt = lambda j: k + j; AN = lambda l: k + w + l - 1

    def rec(l, ctrl, prefix):
        if l == k:
            v = tab[prefix]
            for j in range(w):
                if (v >> j) & 1: q.g("cx", ctrl, Tt(j))
            return
        a, new = A(l), AN(l)
        q.and_compute(ctrl, a, new); rec(l + 1, new, prefix * 2 + 1)
        q.g("cx", ctrl, new); rec(l + 1, new, prefix * 2)
        q.g("cx", ctrl, new); q.and_uncompute(ctrl, a, new)
    if k == 1:
        q.g("x", 0);
        for j in range(w):
            if (tab[0] >> j) & 1: q.g("cx", 0, Tt(j))
        q.g("x", 0)
        for j in range(w):
            if (tab[1] >> j) & 1: q.g("cx", 0, Tt(j))
    else:
        q.g("x", 0); rec(1, 0, 0); q.g("x", 0); rec(1, 0, 1)
    return q.text()


def lookup_ref(k, tab, w):
    """reference as a table (independent of the circuit): (n, f, inw)"""
    def f(X):
        D = X.copy()
        for r in range(len(X)):
            addr = 0
            for l in range(k): addr = addr * 2 + int(X[r, l])
            for j in range(w): D[r, k + j] = X[r, k + j] ^ ((tab[addr] >> j) & 1)
        return D
    return (k + w, f, list(range(k + w)))


def multiplier(n):
    """p := a*b (p starts at 0), schoolbook multiplication: n times (temp = a AND b_i; p[i..i+n] += temp with an AND ripple adder).
    a: 0..n-1, b: n..2n-1, p: 2n..4n-1, temp: 4n..5n-1, carries c_1..c_n: 5n..6n-1 (c_j -> 5n+j-1)."""
    A = lambda j: j; Bq = lambda j: n + j; P = lambda j: 2 * n + j; TMP = lambda j: 4 * n + j; C = lambda j: 5 * n + j - 1
    q = QB(6 * n)
    for i in range(n):
        for j in range(n): q.and_compute(A(j), Bq(i), TMP(j))
        # ripple-add: p_i..p_{i+n-1} += tmp, carry-out ^= p_{i+n}
        pp = lambda j: P(i + j)
        for j in range(n):
            if j >= 1:
                q.g("cx", C(j), TMP(j)); q.g("cx", C(j), pp(j))
            q.and_compute(TMP(j), pp(j), C(j + 1))
            if j >= 1: q.g("cx", C(j), C(j + 1))
        q.g("cx", C(n), P(i + n))
        for j in range(n - 1, -1, -1):
            if j >= 1: q.g("cx", C(j), C(j + 1))
            q.and_uncompute(TMP(j), pp(j), C(j + 1))
            q.g("cx", TMP(j), pp(j))
            if j >= 1:
                q.g("cx", C(j), pp(j)); q.g("cx", C(j), TMP(j))
        for j in range(n): q.and_uncompute(A(j), Bq(i), TMP(j))
    return q.text()


def multiplier_ref(n):
    def f(X):
        D = X.copy()
        for r in range(len(X)):
            a = int("".join(str(int(v)) for v in X[r, :n][::-1]), 2)
            b = int("".join(str(int(v)) for v in X[r, n:2 * n][::-1]), 2)
            p = a * b
            for j in range(2 * n): D[r, 2 * n + j] = (p >> j) & 1
        return D
    return (4 * n, f, list(range(2 * n)))
