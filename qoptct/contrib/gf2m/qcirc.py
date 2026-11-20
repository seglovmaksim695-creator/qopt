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
    ('k', wt, wa, wb, ands, lins)   assertion, no gate: the wires wt hold the bilinear form
                      sum (alpha.wa)(beta.wb) gamma + sum (u.wa) g  (checked by `simulate`; uncompute in
                      `reference_qc`)
    ('m', a, id)      H, measure q[a] into a new classical bit, if 1: X(a)                           0 T
    ('p', id, ws)     if the outcome of measurement id is 1: Z(ws[0]) or CZ(ws[0], ws[1])
                      ('k' + 'm' * |wt| + 'p' ... = erasure of a product register, qoptct.contrib.gf2m.itoh)
The OpenQASM 3 text uses exactly the constructs of qoptct.contrib.gf2m.check (and, for 'r', a longer
``if (c[k]) {...}`` body), which qoptct.verify.prove reads.

On basis states every construct acts classically and without phase provided its precondition holds (AND target 0,
erased wire equal to the product / polynomial).  `simulate` runs the circuit bitsliced (one Python int per wire, one
bit per input lane) and enforces every precondition, so a passing simulation is an exact statement about the
simulated inputs (all measurement outcomes included).  `degree2_exact` turns it into a proof for all inputs when
every wire stays of degree <= 2 (see its docstring).  For the 'k'/'m'/'p' register erasure the simulation checks
the classical precondition ('k'); that the conditional phases equal (-1)^(s.t) is checked exactly and for any m by
`ErasurePhaseChecker` / `erasure_phase_check` (symbolic phase polynomial of the emitted gates; its |0> assumptions on
work wires are checked by `simulate(..., zero_at=checker.assumed_zero)`), and whole circuits are proved by
qoptct.verify.kprove at small m (tests).  The argument for every m is in the docstring of itoh.
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

    # register erasure by bilinear terms (qoptct.contrib.gf2m.itoh)
    def check_product(self, wt, wa, wb, ands, lins):
        """assertion (emits no gate): wires wt hold sum_k (alpha_k.wa)(beta_k.wb) gamma_k + sum (u.wa) g."""
        self.ops.append(('k', tuple(wt), tuple(wa), tuple(wb), tuple(ands), tuple(lins)))

    def measure_reset(self, a):
        """H, measure into a fresh classical bit, conditional X (wire back to |0>); -> measurement id."""
        self._mid = getattr(self, "_mid", 0) + 1
        self.ops.append(('m', a, self._mid))
        return self._mid

    def cphase(self, mid, wires):
        """if outcome `mid` is 1: Z (one wire) or CZ (two wires)."""
        self.ops.append(('p', mid, tuple(wires)))

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
            elif k == 'm':
                ne += 1
            elif k == 'p':
                if len(o[2]) == 2:
                    ncz += 1
                else:
                    nz += 1
        return dict(AND=nA, T=4 * nA, CX=ncx + 3 * nA + ncz, CX_lin=ncx, CZ=ncz, Z=nz, measurements=ne,
                    qubits=self.nq)

    # -- output
    def qasm(self):
        nb = sum(1 for o in self.ops if o[0] in 'erm')
        L = ["OPENQASM 3;", 'include "stdgates.inc";', "qubit[%d] q;" % self.nq, "bit[%d] c;" % max(1, nb)]
        k = 0
        bit_of = {}
        for o in self.ops:
            t = o[0]
            if t == 'k':
                continue
            if t == 'm':
                a = o[1]
                bit_of[o[2]] = k
                L += ["h q[%d];" % a, "c[%d] = measure q[%d];" % (k, a), "if (c[%d]) {x q[%d];}" % (k, a)]
                k += 1
                continue
            if t == 'p':
                ws = o[2]
                g = "z q[%d];" % ws[0] if len(ws) == 1 else "cz q[%d], q[%d];" % ws
                L.append("if (c[%d]) {%s}" % (bit_of[o[1]], g))
                continue
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
            elif k == 'k':                       # uncompute the product register by its monomials
                wt, wa, wb, ands, lins = o[1:]
                for u, g in lins:
                    for j in range(len(wa)):
                        if u >> j & 1:
                            for i in range(len(wt)):
                                if g >> i & 1:
                                    L.append("tof %s %s" % (nm[wa[j]], nm[wt[i]]))
                for al, be, ga in ands:
                    for j in range(len(wa)):
                        if not al >> j & 1:
                            continue
                        for jj in range(len(wb)):
                            if not be >> jj & 1:
                                continue
                            ctl = (nm[wa[j]],) if wa[j] == wb[jj] else (nm[wa[j]], nm[wb[jj]])
                            for i in range(len(wt)):
                                if ga >> i & 1:
                                    L.append("tof " + " ".join(ctl) + " " + nm[wt[i]])
            elif k in 'mp':
                continue
            else:
                a, poly, const = o[1:]
                for mm in poly:
                    L.append("tof " + " ".join(nm[w] for w in mm) + " " + nm[a])
                if const:
                    L.append("tof %s" % nm[a])
        L.append("END")
        return "\n".join(L) + "\n"

    # -- simulation
    def run_lanes(self, state, full, zero_at=None):
        """bitsliced run on `state` (list of nq ints, modified in place); full = all-lanes mask.
        zero_at: {op index: wires} -- wires that must be 0 (all lanes) just before that op (the assumptions of
        ErasurePhaseChecker)."""
        st = state
        for i, o in enumerate(self.ops):
            k = o[0]
            if zero_at and i in zero_at:
                for w in zero_at[i]:
                    if st[w]:
                        raise PreconditionError("op %d: wire q[%d] assumed |0> by the erasure phase check is not 0"
                                                % (i, w))
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
            elif k == 'k':
                wt, wa, wb, ands, lins = o[1:]
                want = [0] * len(wt)

                def par(ws, mask):
                    v = 0
                    j = 0
                    while mask:
                        if mask & 1:
                            v ^= st[ws[j]]
                        mask >>= 1
                        j += 1
                    return v
                for u, g in lins:
                    v = par(wa, u)
                    for t in range(len(wt)):
                        if g >> t & 1:
                            want[t] ^= v
                for al, be, ga in ands:
                    v = par(wa, al) & par(wb, be)
                    for t in range(len(wt)):
                        if ga >> t & 1:
                            want[t] ^= v
                if [st[w] for w in wt] != want:
                    raise PreconditionError("op %d: register-erasure precondition violated" % i)
            elif k == 'm':
                st[o[1]] = 0                     # measured and reset (the outcome only steers phases)
            elif k == 'p':
                pass                             # conditional phase: no action on basis states
            else:
                a, poly, const = o[1:]
                v = full if const else 0
                for mm in poly:
                    v ^= st[mm[0]] if len(mm) == 1 else st[mm[0]] & st[mm[1]]
                if st[a] != v:
                    raise PreconditionError("op %d: register-erasure precondition violated" % i)
                st[a] = 0
        return st

    def simulate(self, values, lanes, zero_at=None):
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
        self.run_lanes(st, full, zero_at)
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


# ---------------------------------------------------------------------------------------------- erasure phases

class ErasurePhaseChecker:
    """EXACT symbolic check of the phase corrections of every product-register erasure, for any m.

    A register erasure is the op group  ('k', wt, wa, wb, ands, lins), then the X measurements ('m') of all wires of
    wt, and CX / conditional Z, CZ ('c', 'x', 'p') until the next op of another kind.  The assertion 'k' says that
    before the group the wires wt hold  K(A, B) = sum_k gamma_k (alpha_k.A)(beta_k.B) + sum g (u.A)  (A, B = the
    values on wa, wb; this classical precondition is checked by `simulate` / the stream check, not here).
    The X measurement of t = K(A, B) leaves the factor (-1)^(s.t) in branch s; the group is a correct erasure iff the
    conditional phases it applies equal exactly (-1)^(s.K(A, B)) for all s, A, B (then every branch is the identity
    on A, B times a data-independent amplitude).

    The check: inside the group only CX / X act (enforced), so every wire holds an affine GF(2) form of the
    variables A, B at the start of the group (wires outside wa, wb start at 0: fresh / released wires, which the
    classical checks guarantee); a 'p' op with outcome s_i on wires holding L1 (L2) contributes s_i*L1 (s_i*L1*L2).
    The collected polynomial  sum_i s_i Q_i(A, B)  is compared with  sum_i s_i K_i(A, B)  coefficient by coefficient
    in the canonical multilinear form (x^2 = x; constants = data-independent branch phases, ignored).  This is a proof
    for all inputs and all outcomes of the phase part; the wires measured must not be read inside the group before
    their measurement.  -> feed(op) per op, result() -> (ok, reason)."""

    def __init__(self):
        self.seg = None
        self.n_groups = 0
        self.n_phase = 0
        self.err = None
        self.i = -1
        self.assumed_zero = {}       # op index of the 'k' op -> wires assumed |0> at the start of that erasure

    def _fail(self, why):
        if self.err is None:
            self.err = "op %d: %s" % (self.i, why)

    def feed(self, o):
        self.i += 1
        if self.err is not None:
            return
        k = o[0]
        if k == 'k':
            self._close()
            self._open(o)
            return
        if self.seg is None:
            if k == 'p':
                self._fail("conditional phase outside a register erasure")
            return
        S = self.seg
        if k == 'm':
            a = o[1]
            if a not in S['opaque'] or a in S['measured']:
                self._fail("measurement of a wire that is not an unmeasured bit of the erased register")
                return
            S['mids'][o[2]] = S['wt_index'][a]
            S['last_rel'] = self.i
            S['measured'].add(a)
            S['opaque'].discard(a)
            S['form'][a] = (0, 0)                # H, measure, conditional X: back to |0>
        elif k == 'c':
            c, t = o[1], o[2]
            if c in S['opaque'] or t in S['opaque']:
                self._fail("CX on an unmeasured wire of the erased register")
                return
            fc = self._get(S, c)
            ft = self._get(S, t)
            S['form'][t] = (ft[0] ^ fc[0], ft[1] ^ fc[1])
        elif k == 'x':
            t = o[1]
            if t in S['opaque']:
                self._fail("X on an unmeasured wire of the erased register")
                return
            ft = self._get(S, t)
            S['form'][t] = (ft[0], ft[1] ^ 1)
        elif k == 'p':
            mid, ws = o[1], o[2]
            if mid not in S['mids']:
                self._fail("phase conditioned on a measurement outside this erasure")
                return
            if any(w in S['opaque'] for w in ws):
                self._fail("phase on an unmeasured wire of the erased register")
                return
            i = S['mids'][mid]
            S['last_rel'] = self.i
            self.n_phase += 1
            Q = S['Q'][i]
            if len(ws) == 1:
                m1, c1 = self._get(S, ws[0])
                Q['lin'] ^= m1
            else:
                m1, c1 = self._get(S, ws[0])
                m2, c2 = self._get(S, ws[1])
                if c1:
                    Q['lin'] ^= m2
                if c2:
                    Q['lin'] ^= m1
                rows = Q['rows']
                x = m1
                while x:
                    u = (x & -x).bit_length() - 1
                    x &= x - 1
                    rows[u] = rows.get(u, 0) ^ m2
        else:
            self._close()

    def _get(self, S, w):
        f = S['form'].get(w)
        if f is None:                            # first touch of a wire outside wa, wb, wt: assumed |0>
            S['zero'].add(w)
            S['ztouch'][w] = self.i
            f = S['form'][w] = (0, 0)
        return f

    def _open(self, o):
        wt, wa, wb, ands, lins = o[1:]
        var = {}
        for w in list(wa) + list(wb):
            if w not in var:
                var[w] = len(var)
        if set(wt) & set(var):
            self._fail("erased register overlaps its operands")
            return
        form = {w: (1 << v, 0) for w, v in var.items()}
        m = len(wt)
        tgt = [dict(rows={}, lin=0) for _ in range(m)]

        def vmask(ws, mask):
            r = 0
            j = 0
            while mask:
                if mask & 1:
                    r |= 1 << var[ws[j]]
                mask >>= 1
                j += 1
            return r
        for u, g in lins:
            um = vmask(wa, u)
            for i in range(m):
                if g >> i & 1:
                    tgt[i]['lin'] ^= um
        for al, be, ga in ands:
            am, bm = vmask(wa, al), vmask(wb, be)
            for i in range(m):
                if ga >> i & 1:
                    rows = tgt[i]['rows']
                    x = am
                    while x:
                        u = (x & -x).bit_length() - 1
                        x &= x - 1
                        rows[u] = rows.get(u, 0) ^ bm
        self.seg = dict(form=form, opaque=set(wt), measured=set(), mids={}, wt=list(wt),
                        wt_index={w: i for i, w in enumerate(wt)}, Q=[dict(rows={}, lin=0) for _ in range(m)],
                        tgt=tgt, start=self.i, zero=set(), ztouch={}, last_rel=self.i)
        self.assumed_zero[self.i] = self.seg['zero']

    @staticmethod
    def _canon_zero(rows, lin):
        """is the multilinear form sum_u x_u * (rows[u] . x) + lin . x zero (x^2 = x)?"""
        upper = {}
        for u, r in rows.items():
            if r >> u & 1:
                lin ^= 1 << u
                r ^= 1 << u
            hi = r >> (u + 1) << (u + 1)
            if hi:
                upper[u] = upper.get(u, 0) ^ hi
            lo = r ^ hi
            while lo:
                v = (lo & -lo).bit_length() - 1
                lo &= lo - 1
                upper[v] = upper.get(v, 0) ^ (1 << u)
        return lin == 0 and not any(upper.values())

    def _close(self):
        S = self.seg
        if S is None:
            return
        self.seg = None
        # CX / X after the last measurement or phase of the group belong to the following computation (the group
        # is only closed by an op of another kind): wires first touched there are no assumption of this check
        S['zero'].difference_update([w for w, t in S['ztouch'].items() if t > S['last_rel']])
        if len(S['measured']) != len(S['wt']):
            self._fail("register erasure from op %d: not every bit measured" % S['start'])
            return
        for i in range(len(S['wt'])):
            Q, T = S['Q'][i], S['tgt'][i]
            rows = dict(Q['rows'])
            for u, r in T['rows'].items():
                rows[u] = rows.get(u, 0) ^ r
            if not self._canon_zero(rows, Q['lin'] ^ T['lin']):
                self._fail("register erasure from op %d: phase polynomial of bit %d != asserted product"
                           % (S['start'], i))
                return
        self.n_groups += 1

    def result(self):
        self._close()
        if self.err is not None:
            return False, self.err
        return True, "exact: %d register erasures, %d conditional phases (all inputs, all outcomes)" % (
            self.n_groups, self.n_phase)


def erasure_phase_check(circ_or_ops):
    """ErasurePhaseChecker over a GF2Circuit (or an ops list) -> (ok, reason)."""
    ops = circ_or_ops.ops if hasattr(circ_or_ops, 'ops') else circ_or_ops
    ch = ErasurePhaseChecker()
    for o in ops:
        ch.feed(o)
    return ch.result()


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
        elif k in 'kmp':
            return False, "op %d: register-erasure ops are not covered by the degree-2 argument" % i
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
