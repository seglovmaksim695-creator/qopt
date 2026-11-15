# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""Exact verification of the measurement-based circuits this module emits.

The emitted OpenQASM 3 uses only four constructs:
  * ``cx q[i], q[j];`` and ``x q[j];``                       (classical linear gates)
  * the logical-AND gadget (4 T, target must be |0>):        q[a] <- q[p] AND q[q]
        h a; t a; cx q,a; tdg a; cx p,a; t a; cx q,a; tdg a; h a; sdg a;
  * the measurement-based erasure (0 T, requires q[a] = q[p] AND q[q]):
        h a; c[k] = measure a; if (c[k]) {cz p, q; x a;}

On computational basis states the AND gadget maps |p, q, 0> -> |p, q, pq> without any phase, and the erasure maps
|p, q, pq> -> |p, q, 0> (outcome 1 picks up (-1)^{pq}, undone by the conditional CZ).  Both facts are exact
identities (checked numerically in tests/test_gf2m.py).  Hence a circuit built from these constructs acts on every
basis input as a classical reversible map without phases, provided every AND finds its target equal to 0 and every
erasure finds its target equal to the AND of its operands.  This checker tracks every wire as a GF(2) polynomial in
the input variables (symbolically, i.e. for ALL inputs at once), enforces both preconditions exactly and compares
the final polynomials with the ANF of the reference .qc (whose phase polynomial must be zero).  Extra qubits must
start and end in |0>.  Any other instruction makes the check fail.
"""
import re

from .anf import pmul, initial

AND_BODY = ["h q[{a}];", "t q[{a}];", "cx q[{q}], q[{a}];", "tdg q[{a}];", "cx q[{p}], q[{a}];", "t q[{a}];",
            "cx q[{q}], q[{a}];", "tdg q[{a}];", "h q[{a}];", "sdg q[{a}];"]

_CX = re.compile(r"^cx q\[(\d+)\], q\[(\d+)\];$")
_X = re.compile(r"^x q\[(\d+)\];$")
_H = re.compile(r"^h q\[(\d+)\];$")
_MEAS = re.compile(r"^c\[(\d+)\] = measure q\[(\d+)\];$")
_IF = re.compile(r"^if \(c\[(\d+)\]\) \{cz q\[(\d+)\], q\[(\d+)\]; x q\[(\d+)\];\}$")


def structural_check(qasm, P_ref, n, inp, clean):
    """-> (ok: bool, reason: str).  P_ref: final polynomials of the n wires of the .qc (same model)."""
    lines = [L.strip() for L in qasm.split("\n") if L.strip()]
    if lines[:2] != ["OPENQASM 3;", 'include "stdgates.inc";']:
        return False, "header"
    m = re.fullmatch(r"qubit\[(\d+)\] q;", lines[2])
    if not m:
        return False, "qubit declaration"
    N = int(m.group(1))
    if N < n:
        return False, "fewer qubits than the .qc"
    i = 3
    nbits = 0
    if i < len(lines):
        mb = re.fullmatch(r"bit\[(\d+)\] c;", lines[i])
        if mb:
            nbits = int(mb.group(1))
            i += 1
    vals = initial(n, inp, clean) + [set() for _ in range(N - n)]
    while i < len(lines):
        L = lines[i]
        mc = _CX.match(L)
        if mc:
            c, t = int(mc.group(1)), int(mc.group(2))
            if c == t or max(c, t) >= N:
                return False, "bad cx: " + L
            vals[t] = vals[t] ^ vals[c]
            i += 1
            continue
        mx = _X.match(L)
        if mx:
            t = int(mx.group(1))
            if t >= N:
                return False, "bad x"
            vals[t] = vals[t] ^ {0}
            i += 1
            continue
        mh = _H.match(L)
        if not mh:
            return False, "unexpected instruction: " + L
        a = int(mh.group(1))
        if a >= N:
            return False, "bad h"
        # erasure?
        if i + 2 < len(lines) and _MEAS.match(lines[i + 1]):
            k, qa = map(int, _MEAS.match(lines[i + 1]).groups())
            mi = _IF.match(lines[i + 2])
            if not mi or qa != a:
                return False, "malformed erasure at line %d" % i
            k2, p, q, a2 = map(int, mi.groups())
            if k2 != k or a2 != a or k >= nbits or len({p, q, a}) < 3 or max(p, q) >= N:
                return False, "malformed erasure at line %d" % i
            if vals[a] != pmul(vals[p], vals[q]):
                return False, "erasure precondition violated at line %d" % i
            vals[a] = set()
            i += 3
            continue
        # AND gadget
        if i + 9 >= len(lines):
            return False, "truncated AND gadget"
        m2 = _CX.match(lines[i + 2])
        m4 = _CX.match(lines[i + 4])
        if not m2 or not m4:
            return False, "malformed AND gadget at line %d" % i
        q = int(m2.group(1))
        p = int(m4.group(1))
        want = [s.format(a=a, p=p, q=q) for s in AND_BODY]
        if lines[i:i + 10] != want or len({p, q, a}) < 3 or max(p, q) >= N:
            return False, "malformed AND gadget at line %d" % i
        if vals[a]:
            return False, "AND target not zero at line %d" % i
        vals[a] = pmul(vals[p], vals[q])
        i += 10
    for w in range(n):
        if vals[w] != P_ref[w]:
            return False, "wire %d differs from the reference" % w
    for w in range(n, N):
        if vals[w]:
            return False, "extra qubit %d not returned to |0>" % w
    return True, "exact symbolic check (all inputs)"


def counts(qasm):
    T = len(re.findall(r"^\s*(t|tdg) q", qasm, re.M))
    CX = len(re.findall(r"\b(cx|cz) q", qasm))
    return T, CX


def numeric_check(qasm, qc_path, clean, trials=3):
    """independent state-vector check with qoptct.verify_qasm (small circuits only). -> fidelity or None."""
    import os
    import tempfile
    from qoptct import verify_qasm
    fd, path = tempfile.mkstemp(suffix=".qasm")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(qasm)
        return float(verify_qasm.check(path, qc_path, trials=trials, clean=clean))
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass
