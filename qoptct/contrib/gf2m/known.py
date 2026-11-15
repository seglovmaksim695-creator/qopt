# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""Table of verified constructions for specific classical functions, recognised by their ANF (not by file name).

Entry: {"fp": fingerprint, "n": wires, "clean": bool, "src": str, "ops": "a1,2,17;e1,2,17,0;c3,4;x5;..."}
ops: a p,q,t = logical AND t <- p q ; e p,q,t = erase t (= p q) ; c c,t = CNOT ; x t = NOT.
The circuit is re-emitted with the gadgets of compile.Emit and verified before use.
"""
import hashlib
import json
import re

from .compile import Emit


def fingerprint(P, n, inp, clean):
    canon = dict(n=n, inp=list(inp), clean=bool(clean), P=[sorted(p) for p in P])
    return hashlib.sha256(json.dumps(canon, separators=(",", ":")).encode()).hexdigest()


def encode(qasm):
    """our gadget-form qasm -> compact op string (raises ValueError on anything else)."""
    from .check import AND_BODY
    lines = [L.strip() for L in qasm.split("\n") if L.strip()]
    nq = int(re.fullmatch(r"qubit\[(\d+)\] q;", lines[2]).group(1))
    i = 4 if lines[3].startswith("bit[") else 3
    ops = []
    while i < len(lines):
        L = lines[i]
        m = re.fullmatch(r"cx q\[(\d+)\], q\[(\d+)\];", L)
        if m:
            ops.append("c%s,%s" % m.groups())
            i += 1
            continue
        m = re.fullmatch(r"x q\[(\d+)\];", L)
        if m:
            ops.append("x%s" % m.group(1))
            i += 1
            continue
        a = int(re.fullmatch(r"h q\[(\d+)\];", L).group(1))
        if "measure" in lines[i + 1]:
            m = re.fullmatch(r"if \(c\[\d+\]\) \{cz q\[(\d+)\], q\[(\d+)\]; x q\[(\d+)\];\}", lines[i + 2])
            ops.append("e%s,%s,%s" % m.groups())
            i += 3
            continue
        q = int(re.match(r"cx q\[(\d+)\]", lines[i + 2]).group(1))
        p = int(re.match(r"cx q\[(\d+)\]", lines[i + 4]).group(1))
        if lines[i:i + 10] != [s.format(a=a, p=p, q=q) for s in AND_BODY]:
            raise ValueError("not a gadget circuit")
        ops.append("a%d,%d,%d" % (p, q, a))
        i += 10
    return nq, ";".join(ops)


def decode(nq, ops):
    E = Emit(nq)
    for op in ops.split(";"):
        k, args = op[0], [int(x) for x in op[1:].split(",")]
        if k == "c":
            E.cx(*args)
        elif k == "x":
            E.x(*args)
        elif k == "a":
            E.AND(*args)
        elif k == "e":
            E.erase(*args)
        else:
            raise ValueError(op)
    return E.text()
