# SPDX-License-Identifier: Apache-2.0
"""exact2 -- exact equivalence proof for measurement-based circuits whose wires stay of GF(2) degree <= 2
(GF(2^m) multipliers and other bilinear maps), against a classical .qc reference, at any size.

    from qoptct.verify.exact2 import prove_degree2
    r = prove_degree2("mult.qasm", "school.qc", clean=False)    # -> dict(status, method='exact2', detail, ...)

Accepted circuit language (anything else: status ERROR, no verdict):
    cx q[c], q[t];      x q[t];
    the logical-AND gadget (Gidney 2018; target must be |0>), exactly the 10 lines of AND_BODY;
    the AND erasure    h q[a]; c[k] = measure q[a]; if (c[k]) {cz q[p], q[q]; x q[a];}   (requires a = p AND q).
Reference: .qc with tof / cx / x (any number of controls) and no H.

Model (as kprove): strict (clean=False): every reference wire is a free input; clean=True: wires not on the .i line
start at 0.  Circuit wires beyond the reference start at 0 and must end at 0.

Proof (status PROVED is a statement for ALL basis inputs and ALL measurement outcomes):
 1. A degree bound is propagated per wire, in the circuit (CX: max, X: same, AND: sum, erasure: 0) and in the
    reference (Toffoli: max(target, sum of the controls)).  If any wire of either exceeds degree 2, exact2 does not
    apply (status UNKNOWN).
 2. Every quantity that must vanish -- AND target before the AND, (erased wire) XOR (p AND q) before an erasure,
    final circuit wire XOR final reference wire, final ancilla -- is then a polynomial of degree <= 2 over GF(2) in
    the input bits.  A polynomial of degree <= 2 that vanishes on all inputs of Hamming weight <= 2 is zero
    (its coefficient of a monomial S, |S| <= 2, is the XOR of its values on the inputs supported in S -- Moebius
    inversion), so the circuit and the reference are run bitsliced on exactly the 1 + n + n(n-1)/2 inputs of
    weight <= 2 (n = number of free input bits), with every precondition enforced.
 3. Phases.  On a basis state whose AND target is 0, the AND gadget maps |p, q, 0> to |p, q, pq> with amplitude 1
    (checked once per call on the 8-dimensional state space, `gadget_check`).  The erasure branch with outcome s
    has amplitude (-1)^(s*a) / sqrt2 before and (-1)^(s*(a + pq)) / sqrt2 = 1/sqrt2 after the correction, since
    a = pq.  Hence every branch maps each basis input x to 2^(-k/2) |f(x), 0> with an amplitude independent of x;
    by step 2, f = the reference map.  So branch K_s = c_s * U_ref on all inputs: the same normal form that
    kprove's BDD engine proves.
A final mismatch with all preconditions held is a genuine counterexample (status REFUTED, the input of weight <= 2
is returned); a violated precondition gives UNKNOWN (the circuit may still be correct in a weaker sense).
Cost: (1 + n + n(n-1)/2)-bit integers per wire; GF(2^163) strict multiplier (n = 489): ~1.2e5 lanes, seconds.
"""
import re
import time

from ._parse import parse_qc

AND_BODY = ["h q[{a}];", "t q[{a}];", "cx q[{q}], q[{a}];", "tdg q[{a}];", "cx q[{p}], q[{a}];", "t q[{a}];",
            "cx q[{q}], q[{a}];", "tdg q[{a}];", "h q[{a}];", "sdg q[{a}];"]

_CX = re.compile(r"^cx q\[(\d+)\], ?q\[(\d+)\];$")
_X = re.compile(r"^x q\[(\d+)\];$")
_H = re.compile(r"^h q\[(\d+)\];$")
_MEAS = re.compile(r"^c\[(\d+)\] = measure q\[(\d+)\];$")
_IF = re.compile(r"^if \(c\[(\d+)\]\) \{cz q\[(\d+)\], ?q\[(\d+)\]; x q\[(\d+)\];\}$")


def parse(text):
    """-> (N, ops) with ops ('c', c, t) | ('x', t) | ('a', p, q, a) | ('e', p, q, a); ValueError otherwise."""
    lines = []
    for s in text.split("\n"):
        s = s.split("//")[0].strip()
        if s:
            lines.append(s)
    if lines[:2] != ["OPENQASM 3;", 'include "stdgates.inc";'] and lines[:2] != ["OPENQASM 3.0;", 'include "stdgates.inc";']:
        raise ValueError("header: OPENQASM 3 + stdgates.inc expected")
    mt = re.fullmatch(r"qubit\[(\d+)\] q;", lines[2]) if len(lines) > 2 else None
    if not mt:
        raise ValueError("qubit[N] q; expected")
    N = int(mt.group(1))
    i = 3
    nbits = 0
    if i < len(lines):
        mb = re.fullmatch(r"bit\[(\d+)\] c;", lines[i])
        if mb:
            nbits = int(mb.group(1))
            i += 1
    used = set()
    ops = []

    def q_ok(*qs):
        for x in qs:
            if not 0 <= x < N:
                raise ValueError("qubit index out of range: %d" % x)
        if len(set(qs)) != len(qs):
            raise ValueError("repeated qubit in one gate")
    while i < len(lines):
        s = lines[i]
        mt = _CX.match(s)
        if mt:
            c, t = int(mt.group(1)), int(mt.group(2))
            q_ok(c, t)
            ops.append(('c', c, t))
            i += 1
            continue
        mt = _X.match(s)
        if mt:
            t = int(mt.group(1))
            q_ok(t)
            ops.append(('x', t))
            i += 1
            continue
        mt = _H.match(s)
        if not mt:
            raise ValueError("unsupported line %d: %s" % (i, s))
        a = int(mt.group(1))
        if i + 1 < len(lines) and _MEAS.match(lines[i + 1]):
            mm = _MEAS.match(lines[i + 1])
            k, a1 = int(mm.group(1)), int(mm.group(2))
            mi = _IF.match(lines[i + 2]) if i + 2 < len(lines) else None
            if not mi:
                raise ValueError("unsupported measurement construct at line %d" % (i + 1))
            k2, p, q, a2 = map(int, mi.groups())
            if a1 != a or a2 != a or k2 != k or k in used or k >= max(nbits, 1):
                raise ValueError("malformed erasure at line %d" % i)
            used.add(k)
            q_ok(p, q, a)
            ops.append(('e', p, q, a))
            i += 3
            continue
        if i + 9 >= len(lines):
            raise ValueError("truncated AND gadget at line %d" % i)
        m2, m4 = _CX.match(lines[i + 2]), _CX.match(lines[i + 4])
        if not (m2 and m4):
            raise ValueError("unsupported construct at line %d" % i)
        q, p = int(m2.group(1)), int(m4.group(1))
        if lines[i:i + 10] != [t.format(a=a, p=p, q=q) for t in AND_BODY]:
            raise ValueError("unsupported construct (not the AND gadget) at line %d" % i)
        q_ok(p, q, a)
        ops.append(('a', p, q, a))
        i += 10
    return N, ops


def gadget_check():
    """the AND gadget on |p, q, 0> -> |p, q, pq> with amplitude 1, for all p, q (dense 8-dim state).  -> max error."""
    import cmath
    import math
    s2 = 1 / math.sqrt(2)
    w = cmath.exp(1j * math.pi / 4)
    one = {"t": w, "tdg": w.conjugate(), "sdg": -1j}
    P, Q, A = 0, 1, 2

    def apply(st, g, qs):
        out = [0j] * 8
        for k, amp in enumerate(st):
            if amp == 0:
                continue
            if g == "cx":
                c, t = qs
                k2 = k ^ (1 << t) if k >> c & 1 else k
                out[k2] += amp
            elif g == "h":
                t = qs[0]
                b = k >> t & 1
                out[k & ~(1 << t)] += amp * s2
                out[k | (1 << t)] += amp * s2 * (-1 if b else 1)
            else:
                t = qs[0]
                out[k] += amp * (one[g] if k >> t & 1 else 1)
        return out
    body = [("h", (A,)), ("t", (A,)), ("cx", (Q, A)), ("tdg", (A,)), ("cx", (P, A)), ("t", (A,)), ("cx", (Q, A)),
            ("tdg", (A,)), ("h", (A,)), ("sdg", (A,))]
    err = 0.0
    for p in (0, 1):
        for q in (0, 1):
            st = [0j] * 8
            st[p << P | q << Q] = 1
            for g, qs in body:
                st = apply(st, g, qs)
            want = p << P | q << Q | (p & q) << A
            for k in range(8):
                err = max(err, abs(st[k] - (1 if k == want else 0)))
    return err


def prove_degree2(qasm, qc, clean=False, max_lanes=5_000_000):
    """exact proof for degree-<=2 circuits (see the module docstring).  qasm, qc: paths (or texts with a newline).
    -> dict(status PROVED | REFUTED | UNKNOWN | ERROR, method='exact2', detail, lanes, time[, counterexample])."""
    t0 = time.time()

    def out(status, detail, **kw):
        d = dict(status=status, method='exact2', detail=detail, time=round(time.time() - t0, 3))
        d.update(kw)
        return d
    try:
        text = qasm if "\n" in qasm else open(qasm).read()
        N, ops = parse(text)
        if "\n" in qc:
            import os
            import tempfile
            fd, path = tempfile.mkstemp(suffix=".qc")
            with os.fdopen(fd, "w") as fh:
                fh.write(qc)
            try:
                names, ins, ref = parse_qc(path)
            finally:
                os.unlink(path)
        else:
            names, ins, ref = parse_qc(qc)
    except (ValueError, KeyError, OSError) as e:
        return out('ERROR', 'exact2 does not apply: %s' % e)
    n = len(names)
    if N < n:
        return out('ERROR', 'circuit has fewer qubits than the reference')
    for g, qs in ref:
        if g not in ('x', 'cx', 'mcx'):
            return out('UNKNOWN', 'exact2 does not apply: reference gate %s is not a Toffoli' % g)
    gerr = gadget_check()
    if gerr > 1e-9:
        return out('ERROR', 'AND gadget check failed (%.3g)' % gerr)
    if clean and ins is not None:
        free = [i for i, nm in enumerate(names) if nm in ins]
    else:
        free = list(range(n))
    nf = len(free)
    L = 1 + nf + nf * (nf - 1) // 2
    if L > max_lanes:
        return out('UNKNOWN', 'exact2: %d lanes > max_lanes' % L)
    # lane 0: zero input; 1..nf: unit vectors; then all pairs
    var = [1 << (1 + v) for v in range(nf)]
    k = 1 + nf
    for i in range(nf):
        mi = 0
        for j in range(i + 1, nf):
            mi |= 1 << (k + j - i - 1)
            var[j] |= 1 << (k + j - i - 1)
        var[i] |= mi
        k += nf - i - 1
    full = (1 << L) - 1

    def init():
        st = [0] * N
        deg = [0] * N
        for v, w in enumerate(free):
            st[w] = var[v]
            deg[w] = 1
        return st, deg
    # reference
    rs, rd = init()
    for g, qs in ref:
        t = qs[-1]
        if g == 'x':
            rs[t] ^= full
            continue
        v = full
        d = 0
        for c in qs[:-1]:
            v &= rs[c]
            d += rd[c]
        rs[t] ^= v
        rd[t] = max(rd[t], d)
        if rd[t] > 2:
            return out('UNKNOWN', 'exact2 does not apply: reference wire %s has degree > 2' % names[t])
    # circuit
    st, deg = init()
    for i, o in enumerate(ops):
        kk = o[0]
        if kk == 'c':
            st[o[2]] ^= st[o[1]]
            deg[o[2]] = max(deg[o[2]], deg[o[1]])
        elif kk == 'x':
            st[o[1]] ^= full
        elif kk == 'a':
            p, q, a = o[1:]
            if deg[p] + deg[q] > 2:
                return out('UNKNOWN', 'exact2 does not apply: op %d (AND) reaches degree > 2' % i)
            if st[a]:
                return out('UNKNOWN', 'op %d: AND target not |0> (precondition violated)' % i)
            st[a] = st[p] & st[q]
            deg[a] = deg[p] + deg[q]
        else:
            p, q, a = o[1:]
            if deg[p] + deg[q] > 2:
                return out('UNKNOWN', 'exact2 does not apply: op %d (erasure) reaches degree > 2' % i)
            if st[a] != st[p] & st[q]:
                return out('UNKNOWN', 'op %d: erased wire != p AND q (precondition violated)' % i)
            st[a] = 0
            deg[a] = 0
    bad = 0
    where = None
    for w in range(N):
        want = rs[w] if w < n else 0
        if st[w] != want:
            d = st[w] ^ want
            bad |= d
            if where is None:
                where = names[w] if w < n else 'ancilla q[%d]' % w
    if bad:
        lane = (bad & -bad).bit_length() - 1
        x = [(var[v] >> lane) & 1 for v in range(nf)]
        cex = {names[free[v]]: 1 for v in range(nf) if x[v]}
        return out('REFUTED', 'wire %s differs from the reference on an input of weight <= 2 (all preconditions held)'
                   % where, lanes=L, counterexample=cex)
    return out('PROVED', 'degree <= 2 (circuit and reference); all %d inputs of weight <= 2 agree with every '
               'AND / erasure precondition (%d ops) => all inputs, all measurement outcomes' % (L, len(ops)),
               lanes=L)


def main(argv=None):
    import argparse
    import json
    import sys
    ap = argparse.ArgumentParser(prog='python -m qoptct.verify.exact2', description=__doc__.split('\n')[0])
    ap.add_argument('qasm')
    ap.add_argument('qc')
    ap.add_argument('--clean', action='store_true')
    a = ap.parse_args(argv)
    r = prove_degree2(a.qasm, a.qc, clean=a.clean)
    print(json.dumps(r))
    return {'PROVED': 0, 'REFUTED': 1, 'ERROR': 3}.get(r['status'], 2)


if __name__ == '__main__':
    import sys
    sys.exit(main())
