# SPDX-License-Identifier: Apache-2.0
"""Independent sparse simulator (state = dict basis-int -> amplitude; bit q of the key is qubit q).

It shares no code with the path-sum engine except the text parser, and serves two purposes:
  * check_sparse: a fast sampled check of a circuit against its .qc reference on random superpositions of inputs
    (measurement outcomes sampled with the Born rule); the extra qubits must end in one basis state;
  * confirm_branch: independent confirmation of a counterexample (x, m) found by the exact prover.
Gates: x y z s sdg t tdg h cx cz ccx ccz mcx mcz swap id reset; `c[k] = measure q[i]`; `if (...)` with && & ^ ! ( ).
"""
import cmath
import math
import random

from ._parse import parse_qasm_text, parse_qc, eval_cond, clean_wires

W8 = [cmath.exp(1j * math.pi * k / 4) for k in range(8)]
PH = {'z': 4, 's': 2, 'sdg': 6, 't': 1, 'tdg': 7}
R2 = 1 / math.sqrt(2)


class TooManyTerms(RuntimeError):
    pass


def apply(st, name, q, cap=1 << 20):
    """Apply one unitary gate to a sparse state; returns a new dict."""
    if name == 'h':
        b = 1 << q[0]
        out = {}
        for k, a in st.items():
            k0 = k & ~b
            k1 = k | b
            out[k0] = out.get(k0, 0) + a * R2
            out[k1] = out.get(k1, 0) + (-a if k & b else a) * R2
        out = {k: v for k, v in out.items() if abs(v) > 1e-14}
        if len(out) > cap:
            raise TooManyTerms('sparse state exceeds %d terms' % cap)
        return out
    out = {}
    if name in PH:
        f = W8[PH[name]]
        b = 1 << q[0]
        for k, a in st.items():
            out[k] = a * f if k & b else a
        return out
    if name in ('x', 'cx', 'ccx', 'mcx'):
        cm = 0
        for c in q[:-1]:
            cm |= 1 << c
        tb = 1 << q[-1]
        for k, a in st.items():
            out[k ^ tb if k & cm == cm else k] = a
        return out
    if name in ('cz', 'ccz', 'mcz'):
        cm = 0
        for c in q:
            cm |= 1 << c
        for k, a in st.items():
            out[k] = -a if k & cm == cm else a
        return out
    if name == 'y':     # Y = i X Z
        b = 1 << q[0]
        for k, a in st.items():
            out[k ^ b] = a * (-1j if k & b else 1j)
        return out
    if name == 'swap':
        i, j = q
        for k, a in st.items():
            bi, bj = (k >> i) & 1, (k >> j) & 1
            if bi != bj:
                k ^= (1 << i) | (1 << j)
            out[k] = a
        return out
    if name == 'id':
        return dict(st)
    raise ValueError('unsupported gate: ' + name)


def _project(st, q, b):
    bb = 1 << q
    return {k: a for k, a in st.items() if bool(k & bb) == bool(b)}


def _norm2(st):
    return sum(abs(a) ** 2 for a in st.values())


def run_sampled(ops, st, rng, cap=1 << 20):
    """Run with Born-rule outcomes (state renormalised after each measurement / reset)."""
    cbits = {}
    for op in ops:
        if op[0] == 'G':
            if op[1] == 'reset':
                q = op[2][0]
                p1 = sum(abs(a) ** 2 for k, a in st.items() if k >> q & 1) / max(_norm2(st), 1e-300)
                r = 1 if rng.random() < p1 else 0
                st = _project(st, q, r)
                nrm = math.sqrt(_norm2(st))
                st = {(k & ~(1 << q)): a / nrm for k, a in st.items()}
            else:
                st = apply(st, op[1], op[2], cap)
        elif op[0] == 'M':
            q, c = op[1], op[2]
            p1 = sum(abs(a) ** 2 for k, a in st.items() if k >> q & 1) / max(_norm2(st), 1e-300)
            r = 1 if rng.random() < p1 else 0
            st = _project(st, q, r)
            nrm = math.sqrt(_norm2(st))
            st = {k: a / nrm for k, a in st.items()}
            cbits[c] = r
        else:
            if eval_cond(op[1], cbits):
                for name, q in op[2]:
                    st = apply(st, name, q, cap)
    return st, cbits


def run_branch_scaled(ops, st, outcomes, cap=1 << 20):
    """Run the branch with the given outcomes (list indexed by classical bit).

    The (unnormalised) branch vector is returned as (state, lg) with branch = 2**lg * state: after every measurement the
    state is rescaled to max |amplitude| = 1 and the scale is kept as a base-2 logarithm, so chains of thousands of
    measurements neither underflow nor lose small terms (the 1e-14 cut-off of `apply` is relative to that maximum).
    A branch of amplitude zero gives ({}, None). A reset is allowed only on a wire with a definite value (ValueError)."""
    cbits = {}
    lg = 0.0
    for op in ops:
        if op[0] == 'G':
            if op[1] == 'reset':
                q = op[2][0]
                vals = {(k >> q) & 1 for k in st}
                if len(vals) > 1:
                    raise ValueError('reset of a superposed wire q[%d] in a fixed branch' % q)
                st = {(k & ~(1 << q)): a for k, a in st.items()}
            else:
                st = apply(st, op[1], op[2], cap)
        elif op[0] == 'M':
            q, c = op[1], op[2]
            r = outcomes[c]
            st = _project(st, q, r)
            cbits[c] = r
            mx = max((abs(a) for a in st.values()), default=0.0)
            if mx < 1e-13:
                return {}, None
            st = {k: a / mx for k, a in st.items() if abs(a) > 1e-14 * mx}
            lg += math.log2(mx)
        else:
            if eval_cond(op[1], cbits):
                for name, q in op[2]:
                    st = apply(st, name, q, cap)
    return st, lg


def run_branch(ops, st, outcomes, cap=1 << 20):
    """Run the branch with the given outcomes (list indexed by classical bit); returns the unnormalised branch vector.
    Internally renormalised after every measurement (see run_branch_scaled, which also returns the scale as log2 and
    should be used for very long chains, where 2**lg itself underflows)."""
    out, lg = run_branch_scaled(ops, st, outcomes, cap)
    if lg is None:
        return {}
    f = 2.0 ** lg
    return {k: a * f for k, a in out.items()}


def apply_ref(st, gates, cap=1 << 20):
    for name, q in gates:
        st = apply(st, name, q, cap)
    return st


def _key(x, n):
    return sum((x.get(q, 0) & 1) << q for q in range(n))


def confirm_branch(qasm, qc, xs, m, n, cap=1 << 18):
    """Independent confirmation of a counterexample: for branch m and basis inputs xs (dicts wire -> bit), the branch
    vector must be t * (U_ref|x> (x) |a>) with one basis ancilla state a per input and one number t for all inputs.
    Returns the largest deviation (> 1e-6 means the counterexample is confirmed)."""
    txt = open(qasm).read()
    N, nb, ops = parse_qasm_text(txt)
    names, ins, ref = parse_qc(qc)
    m = list(m) + [0] * max(0, nb - len(m))
    worst = 0.0
    ts = []
    mask = (1 << n) - 1
    for x in xs:
        k0 = _key(x, n)
        out, lg = run_branch_scaled(ops, {k0: 1.0 + 0j}, m, cap)   # renormalised: long measurement chains are fine
        want = apply_ref({k0: 1.0 + 0j}, ref, cap)
        tot = _norm2(out)
        if lg is None or tot < 1e-12:
            ts.append(None)
            continue
        wa = {}
        for k, a in out.items():
            wa[k >> n] = wa.get(k >> n, 0) + abs(a) ** 2
        anc = max(wa, key=wa.get)
        worst = max(worst, (tot - wa[anc]) / tot)
        col = {k & mask: a for k, a in out.items() if k >> n == anc}
        i = max(want, key=lambda k: abs(want[k]))
        t = col.get(i, 0j) / want[i]
        if abs(t) < 1e-9:
            worst = max(worst, 1.0)
            ts.append(None)
            continue
        keys = set(col) | set(want)
        dev = max(abs(col.get(k, 0j) - t * want.get(k, 0j)) for k in keys) / abs(t)
        worst = max(worst, dev)
        ts.append((math.log2(abs(t)) + lg, t / abs(t)))     # t = 2**(log2|t|) * phase, kept in log form
    nz = [t for t in ts if t is not None]
    if nz and len(nz) < len(ts):
        worst = max(worst, 1.0)
    for lgt, ph in nz:
        d = lgt - nz[0][0]
        worst = max(worst, abs(2.0 ** d - 1) if d < 60 else 1e6, abs(ph - nz[0][1]))
    return float(min(worst, 1e6))


def check_sparse(qasm, qc, clean=False, trials=8, terms=64, seed=1, cap=1 << 20, tol=1e-9):
    """Sampled check: random superpositions of `terms` basis inputs (all inputs if there are at most `terms`),
    Born-rule outcomes; the output must equal (U_ref|psi>) (x) |a> for one basis state a of the extra qubits.
    -> dict(ok, F (worst overlap), trials, error?). Not a proof: use qoptct.verify.prove for that."""
    from ._parse import QasmFile
    with QasmFile(qasm) as path:
        txt = open(path).read()
    try:
        N, nb, ops = parse_qasm_text(txt)
        names, ins, ref = parse_qc(qc)
    except (ValueError, KeyError) as e:
        return dict(ok=False, F=0.0, error=str(e))
    n = len(names)
    if N < n:
        return dict(ok=False, F=0.0, error='fewer qubits than the reference')
    free = [q for q in range(n) if q not in clean_wires(names, ins, clean)]
    rng = random.Random(seed)
    worst = 1.0
    try:
        for _ in range(trials):
            if len(free) <= 20 and (1 << len(free)) <= terms:
                xs = list(range(1 << len(free)))
            else:
                xs = list({rng.getrandbits(len(free)) for _ in range(terms)})
            amp = {}
            for xi in xs:
                k = 0
                for j, q in enumerate(free):
                    if xi >> j & 1:
                        k |= 1 << q
                amp[k] = complex(rng.gauss(0, 1), rng.gauss(0, 1))
            nrm = math.sqrt(_norm2(amp))
            amp = {k: a / nrm for k, a in amp.items()}
            out, _ = run_sampled(ops, dict(amp), rng, cap)
            want = apply_ref(dict(amp), ref, cap)
            wa = {}
            for k, a in out.items():
                wa[k >> n] = wa.get(k >> n, 0) + abs(a) ** 2
            anc = max(wa, key=wa.get) if wa else 0
            ov = sum(want.get(k & ((1 << n) - 1), 0j).conjugate() * a for k, a in out.items() if k >> n == anc)
            worst = min(worst, abs(ov))
    except TooManyTerms as e:
        return dict(ok=False, F=None, error=str(e))
    except ValueError as e:
        return dict(ok=False, F=0.0, error=str(e))
    return dict(ok=bool(worst > 1 - tol), F=worst, trials=trials)
