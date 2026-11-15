# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""qcirc: the measurement-based circuit object returned by the GF(2^m) generators of qoptct v0.5.

Operations (ops list):
    ('c', c, t)       CX  c -> t
    ('x', t)          X
    ('a', p, q, a)    logical AND (Gidney): q[a] <- q[p] AND q[q]; q[a] must be |0>     4 T, 3 CX (as emitted)
    ('e', p, q, a)    measurement erasure of q[a] = q[p] AND q[q]: H, measure, if 1: CZ(p, q), X(a)      0 T
    ('r', a, poly)    measurement erasure of a wire that equals a quadratic polynomial Q of other wires:
                      H, measure, if 1: (-1)^Q as Z / CZ gates, X(a)                                 0 T
                      poly = tuple of monomials (w,) or (u, v); a constant monomial () is a global phase of
                      that branch and is dropped (the branch is still U up to a phase).
The OpenQASM 3 text uses exactly the constructs of qoptct.contrib.gf2m.check (and, for 'r', a longer
``if (c[k]) {...}`` body), which qoptct.verify.prove reads.

On basis states every construct acts classically and without phase provided its precondition holds (AND target 0,
erased wire equal to the product / polynomial).  `simulate` runs the circuit bitsliced (one Python int per wire, one
bit per input lane) and enforces every precondition, so a passing simulation is an exact statement about the
simulated inputs (all measurement outcomes included).  `degree2_exact` turns it into a proof for all inputs when
every wire stays of degree <= 2 (see its docstring).
"""
import re


class PreconditionError(AssertionError):
    pass


AND_BODY = ["h q[{a}];", "t q[{a}];", "cx q[{q}], q[{a}];", "tdg q[{a}];", "cx q[{p}], q[{a}];", "t q[{a}];",
            "cx q[{q}], q[{a}];", "tdg q[{a}];", "h q[{a}];", "sdg q[{a}];"]


class GF2Circuit:
    """ops + qubit count + named registers.  regs: {name: [wire, ...]} (bit i of the register value on wire[i]).
    meta: free-form description (field, modulus, algorithm, ...)."""

    def __init__(self, nq, ops=None, regs=None, meta=None):
        self.nq = nq
        self.ops = [] if ops is None else ops
        self.regs = {} if regs is None else regs
        self.meta = {} if meta is None else meta

    # -- building
    def cx(self, c, t):
        self.ops.append(('c', c, t))

    def x(self, t):
        self.ops.append(('x', t))

    def AND(self, p, q, a):
        self.ops.append(('a', p, q, a))

    def erase(self, p, q, a):
        self.ops.append(('e', p, q, a))

    def rerase(self, a, poly):
        const = sum(1 for mm in poly if not len(mm)) & 1
        self.ops.append(('r', a, tuple(tuple(mm) for mm in poly if len(mm)), const))

    # -- costs
    def counts(self):
        """AND (= Toffoli-equivalents), T, CX_lin (plain cx), CZ (in erasures), CX = all cx + cz in the text
        (3 cx per AND gadget), measurements, qubits."""
        nA = ncx = ne = ncz = nz = 0
        for o in self.ops:
            k = o[0]
            if k == 'c':
                ncx += 1
            elif k == 'a':
                nA += 1
            elif k == 'e':
                ne += 1
                ncz += 1
            elif k == 'r':
                ne += 1
                for mm in o[2]:
                    if len(mm) == 2:
                        ncz += 1
                    else:
                        nz += 1
        return dict(AND=nA, T=4 * nA, CX=ncx + 3 * nA + ncz, CX_lin=ncx, CZ=ncz, Z=nz, measurements=ne,
                    qubits=self.nq)

    # -- output
    def qasm(self):
        nb = sum(1 for o in self.ops if o[0] in 'er')
        L = ["OPENQASM 3;", 'include "stdgates.inc";', "qubit[%d] q;" % self.nq, "bit[%d] c;" % max(1, nb)]
        k = 0
        for o in self.ops:
            t = o[0]
            if t == 'c':
                L.append("cx q[%d], q[%d];" % (o[1], o[2]))
            elif t == 'x':
                L.append("x q[%d];" % o[1])
            elif t == 'a':
                p, q, a = o[1:]
                L += [s.format(a=a, p=p, q=q) for s in AND_BODY]
            elif t == 'e':
                p, q, a = o[1:]
                L += ["h q[%d];" % a, "c[%d] = measure q[%d];" % (k, a),
                      "if (c[%d]) {cz q[%d], q[%d]; x q[%d];}" % (k, p, q, a)]
                k += 1
            else:
                a, poly = o[1], o[2]
                body = []
                for mm in poly:
                    body.append("z q[%d];" % mm[0] if len(mm) == 1 else "cz q[%d], q[%d];" % mm)
                body.append("x q[%d];" % a)
                L += ["h q[%d];" % a, "c[%d] = measure q[%d];" % (k, a), "if (c[%d]) {%s}" % (k, " ".join(body))]
                k += 1
        return "\n".join(L) + "\n"

    def reference_qc(self, inputs=None):
        """Toffoli-level reference (.qc text) of the same classical map: AND -> tof, erasures -> uncompute.
        inputs: wire list for the .i line (default: all wires of all registers)."""
        nm = ["q%d" % i for i in range(self.nq)]
        if inputs is None:
            inputs = sorted(w for ws in self.regs.values() for w in ws)
        L = [".v " + " ".join(nm), ".i " + " ".join(nm[i] for i in inputs), "", "BEGIN"]
        for o in self.ops:
            k = o[0]
            if k == 'x':
                L.append("tof %s" % nm[o[1]])
            elif k == 'c':
                L.append("tof %s %s" % (nm[o[1]], nm[o[2]]))
            elif k in 'ae':
                L.append("tof %s %s %s" % (nm[o[1]], nm[o[2]], nm[o[3]]))
            else:
                a, poly, const = o[1:]
                for mm in poly:
                    L.append("tof " + " ".join(nm[w] for w in mm) + " " + nm[a])
                if const:
                    L.append("tof %s" % nm[a])
        L.append("END")
        return "\n".join(L) + "\n"

    # -- simulation
    def run_lanes(self, state, full):
        """bitsliced run on `state` (list of nq ints, modified in place); full = all-lanes mask."""
        st = state
        for i, o in enumerate(self.ops):
            k = o[0]
            if k == 'c':
                st[o[2]] ^= st[o[1]]
            elif k == 'x':
                st[o[1]] ^= full
            elif k == 'a':
                p, q, a = o[1:]
                if st[a]:
                    raise PreconditionError("op %d: AND target not zero" % i)
                st[a] = st[p] & st[q]
            elif k == 'e':
                p, q, a = o[1:]
                if st[a] != st[p] & st[q]:
                    raise PreconditionError("op %d: erasure precondition violated" % i)
                st[a] = 0
            else:
                a, poly, const = o[1:]
                v = full if const else 0
                for mm in poly:
                    v ^= st[mm[0]] if len(mm) == 1 else st[mm[0]] & st[mm[1]]
                if st[a] != v:
                    raise PreconditionError("op %d: register-erasure precondition violated" % i)
                st[a] = 0
        return st

    def simulate(self, values, lanes):
        """values: {register name: list of `lanes` ints}; other wires start at 0.  Returns {name: list of ints} at
        the end; raises PreconditionError on a violated precondition or a dirty non-register wire at the end."""
        full = (1 << lanes) - 1
        st = [0] * self.nq
        for name, vals in values.items():
            ws = self.regs[name]
            for i, w in enumerate(ws):
                x = 0
                for ln, v in enumerate(vals):
                    if v >> i & 1:
                        x |= 1 << ln
                st[w] = x
        self.run_lanes(st, full)
        regw = set(w for ws in self.regs.values() for w in ws)
        for w in range(self.nq):
            if w not in regw and st[w]:
                raise PreconditionError("ancilla q[%d] not returned to 0" % w)
        out = {}
        for name, ws in self.regs.items():
            vals = [0] * lanes
            for i, w in enumerate(ws):
                x = st[w]
                while x:
                    ln = (x & -x).bit_length() - 1
                    x &= x - 1
                    vals[ln] |= 1 << i
            out[name] = vals
        return out


# ---------------------------------------------------------------------------------------------- parsing

_CX = re.compile(r"^cx q\[(\d+)\], q\[(\d+)\];$")
_X = re.compile(r"^x q\[(\d+)\];$")
_H = re.compile(r"^h q\[(\d+)\];$")
_MEAS = re.compile(r"^c\[(\d+)\] = measure q\[(\d+)\];$")
_IF = re.compile(r"^if \(c\[(\d+)\]\) \{cz q\[(\d+)\], q\[(\d+)\]; x q\[(\d+)\];\}$")


def parse_qasm(text):
    """OpenQASM 3 text of the four constructs (cx, x, AND gadget, erasure) -> (nq, ops).  Raises ValueError on
    anything else (used for the circuits of qoptct.contrib.gf2m.compile)."""
    lines = [s.strip() for s in text.split("\n") if s.strip()]
    if lines[:2] != ["OPENQASM 3;", 'include "stdgates.inc";']:
        raise ValueError("header")
    N = int(re.fullmatch(r"qubit\[(\d+)\] q;", lines[2]).group(1))
    i = 3
    if i < len(lines) and re.fullmatch(r"bit\[(\d+)\] c;", lines[i]):
        i += 1
    ops = []
    while i < len(lines):
        s = lines[i]
        mt = _CX.match(s)
        if mt:
            ops.append(('c', int(mt.group(1)), int(mt.group(2))))
            i += 1
            continue
        mt = _X.match(s)
        if mt:
            ops.append(('x', int(mt.group(1))))
            i += 1
            continue
        mt = _H.match(s)
        if not mt:
            raise ValueError("unexpected line: " + s)
        a = int(mt.group(1))
        if i + 2 < len(lines) and _MEAS.match(lines[i + 1]):
            mi = _IF.match(lines[i + 2])
            if not mi:
                raise ValueError("unexpected erasure body: " + lines[i + 2])
            _, p, q, a2 = map(int, mi.groups())
            if a2 != a:
                raise ValueError("erasure wire mismatch")
            ops.append(('e', p, q, a))
            i += 3
            continue
        m2 = _CX.match(lines[i + 2])
        m4 = _CX.match(lines[i + 4])
        q = int(m2.group(1))
        p = int(m4.group(1))
        if lines[i:i + 10] != [t.format(a=a, p=p, q=q) for t in AND_BODY]:
            raise ValueError("unexpected AND gadget at line %d" % i)
        ops.append(('a', p, q, a))
        i += 10
    return N, ops


# ---------------------------------------------------------------------------------------------- exact check

def degree2_exact(circ, in_wires, ref):
    """EXACT check (all inputs, all measurement outcomes) of a circuit whose wires all stay of degree <= 2.

    in_wires: input variables (wires holding a free bit at the start); all other wires start at 0.
    ref(lane_vars) -> {wire: expected final int} for the listed output wires, computed by the caller on the same
    lanes (lane_vars[w] = int of the variable on wire w); every wire not in `ref` must end equal to its start value
    (input wires restored, ancillas back to 0).

    Argument (as lab k5/exact2): a degree bound is propagated per wire (CX: max, AND: sum, erasures: 0); if every
    wire keeps degree <= 2, every compared quantity -- AND targets, erased wire xor its product/polynomial, final
    wire xor reference -- is a polynomial of degree <= 2 over GF(2) in the inputs (the reference must be of degree
    <= 2 as well, which holds for a bilinear map), and a polynomial of degree <= 2 that vanishes on all inputs of
    Hamming weight <= 2 vanishes identically (Moebius inversion).  So all inputs of weight <= 2 are simulated
    (bitsliced).  -> (ok, reason)."""
    n = len(in_wires)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    L = 1 + n + len(pairs)
    var = [0] * n
    for v in range(n):
        var[v] |= 1 << (1 + v)
    for k, (i, j) in enumerate(pairs):
        var[i] |= 1 << (1 + n + k)
        var[j] |= 1 << (1 + n + k)
    full = (1 << L) - 1
    st = [0] * circ.nq
    deg = [0] * circ.nq
    lane_vars = {}
    for v, w in enumerate(in_wires):
        st[w] = var[v]
        deg[w] = 1
        lane_vars[w] = var[v]
    for i, o in enumerate(circ.ops):
        k = o[0]
        if k == 'c':
            st[o[2]] ^= st[o[1]]
            deg[o[2]] = max(deg[o[2]], deg[o[1]])
        elif k == 'x':
            st[o[1]] ^= full
        elif k == 'a':
            p, q, a = o[1:]
            if st[a]:
                return False, "op %d: AND target not zero" % i
            st[a] = st[p] & st[q]
            deg[a] = deg[p] + deg[q]
            if deg[a] > 2:
                return False, "op %d: degree > 2" % i
        elif k == 'e':
            p, q, a = o[1:]
            if deg[p] + deg[q] > 2:
                return False, "op %d: degree > 2 in erasure" % i
            if st[a] != st[p] & st[q]:
                return False, "op %d: erasure precondition" % i
            st[a] = 0
            deg[a] = 0
        else:
            a, poly, const = o[1:]
            v = full if const else 0
            for mm in poly:
                if len(mm) == 1:
                    v ^= st[mm[0]]
                else:
                    if deg[mm[0]] + deg[mm[1]] > 2:
                        return False, "op %d: degree > 2 in register erasure" % i
                    v ^= st[mm[0]] & st[mm[1]]
            if st[a] != v:
                return False, "op %d: register-erasure precondition" % i
            st[a] = 0
            deg[a] = 0
    want = ref(lane_vars)
    start = {w: lane_vars.get(w, 0) for w in range(circ.nq)}
    for w in range(circ.nq):
        exp = want[w] if w in want else start[w]
        if st[w] != exp:
            return False, "wire q[%d] wrong at the end" % w
    return True, "exact (degree <= 2, all %d inputs of weight <= 2)" % L
