# SPDX-License-Identifier: Apache-2.0
"""Small intermediate representation for the OpenQASM 3 subset written by qoptct.

Items:
    ('G', name, qubits)            unconditional gate (h x z s sdg t tdg cx cz)
    ('M', qubit, bit)              c[bit] = measure q[qubit];
    ('R', qubit)                   reset q[qubit];
    ('IF', cond, body)             if (cond) { body }   cond = string in canonical form, body = [(name, qubits)]

Conditions are conjunctions (&&) of XOR groups:  (c[0] ^ c[3]) && (c[5]).  ``normalize`` renames the classical bits so that
every measurement writes a fresh bit; passes that reason about conditions rely on that.
"""
import re

GATES1 = {"h", "x", "z", "s", "sdg", "t", "tdg"}
GATES2 = {"cx", "cz"}
DIAG = {"z", "s", "sdg", "t", "tdg", "cz"}


class ParseError(ValueError):
    pass


def _gate(s):
    s = s.strip().rstrip(";").strip()
    m = re.fullmatch(r"([a-z]+)\s+(q\[\d+\](?:\s*,\s*q\[\d+\])*)", s)
    if not m:
        raise ParseError("bad gate: " + s)
    name = m.group(1); qs = tuple(int(x) for x in re.findall(r"q\[(\d+)\]", m.group(2)))
    if not ((name in GATES1 and len(qs) == 1) or (name in GATES2 and len(qs) == 2 and qs[0] != qs[1])):
        raise ParseError("unsupported gate: " + s)
    return name, qs


def parse_cond(s):
    """'(c[0] ^ c[2]) && c[3]' -> [[0, 2], [3]] (AND of XOR groups). Anything else raises ParseError."""
    s = s.strip()
    if not re.fullmatch(r"[\s()c\[\]\d&^]+", s):
        raise ParseError("unsupported condition: " + s)
    groups = []
    for g in s.split("&&"):
        g = g.replace("(", " ").replace(")", " ").strip()
        if "&" in g:
            # a single '&' is a bitwise AND of single bits: treat as separate groups
            parts = [p.strip() for p in g.split("&")]
        else:
            parts = [g]
        for p in parts:
            toks = [t.strip() for t in p.split("^")]
            bits = []
            for t in toks:
                m = re.fullmatch(r"c\[(\d+)\]", t)
                if not m:
                    raise ParseError("unsupported condition: " + s)
                bits.append(int(m.group(1)))
            groups.append(bits)
    return groups


def fmt_cond(groups):
    return " && ".join("(" + " ^ ".join(f"c[{b}]" for b in g) + ")" for g in groups)


def parse(text):
    """-> (N, items). Bits are renamed to fresh indices per measurement (normalisation)."""
    N = None; items = []; cur = {}; nb = 0
    for raw in text.splitlines():
        L = raw.strip()
        if not L or L.startswith("//"):
            continue
        if L.startswith(("OPENQASM", "include")):
            continue
        m = re.fullmatch(r"qubit\[(\d+)\]\s+q;", L)
        if m:
            if N is not None:
                raise ParseError("two qubit registers")
            N = int(m.group(1)); continue
        if re.fullmatch(r"bit\[\d+\]\s+c;", L):
            continue
        m = re.fullmatch(r"c\[(\d+)\]\s*=\s*measure\s+q\[(\d+)\];", L)
        if m:
            cur[int(m.group(1))] = nb; items.append(("M", int(m.group(2)), nb)); nb += 1; continue
        m = re.fullmatch(r"reset\s+q\[(\d+)\];", L)
        if m:
            items.append(("R", int(m.group(1)))); continue
        m = re.fullmatch(r"if\s*\((.*)\)\s*\{(.*)\}", L)
        if m:
            groups = parse_cond(m.group(1))
            try:
                groups = [[cur[b] for b in g] for g in groups]
            except KeyError:
                raise ParseError("condition on an unmeasured bit: " + L)
            # drop duplicated bits inside a XOR group (b ^ b = 0)
            gg = []
            for g in groups:
                s = set()
                for b in g:
                    s ^= {b}
                gg.append(sorted(s))
            body = [_gate(b) for b in m.group(2).split(";") if b.strip()]
            if not body:
                continue
            if any(not g for g in gg):
                continue  # condition is identically 0
            items.append(("IF", fmt_cond(gg), body)); continue
        items.append(("G",) + _gate(L))
    if N is None:
        raise ParseError("no qubit register")
    for it in items:
        qs = (it[1],) if it[0] in ("M", "R") else (it[2] if it[0] == "G" else tuple(q for _, b in it[2] for q in b))
        if any(q >= N for q in qs):
            raise ParseError("qubit index out of range")
    return N, items


def nbits(items):
    return 1 + max([it[2] for it in items if it[0] == "M"], default=0)


def emit(N, items):
    nb = sum(1 for it in items if it[0] == "M")
    L = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{N}] q;"]
    if nb:
        L.append(f"bit[{nbits(items)}] c;")
    for it in items:
        if it[0] == "G":
            L.append(f"{it[1]} " + ", ".join(f"q[{q}]" for q in it[2]) + ";")
        elif it[0] == "M":
            L.append(f"c[{it[2]}] = measure q[{it[1]}];")
        elif it[0] == "R":
            L.append(f"reset q[{it[1]}];")
        else:
            body = " ".join(f"{g} " + ", ".join(f"q[{q}]" for q in qs) + ";" for g, qs in it[2])
            L.append(f"if ({it[1]}) {{{body}}}")
    return "\n".join(L) + "\n"


def normalize(text):
    N, items = parse(text)
    return emit(N, items)


def counts_text(text):
    """(T, CX) of a QASM text: T = t/tdg gates, CX = cx/cz gates (conditional ones included)."""
    return (len(re.findall(r"\b(?:t|tdg)\s+q\[", text)), len(re.findall(r"\b(?:cx|cz)\s+q\[", text)))


def counts(ops):
    t = cx = 0
    for o in ops:
        gs = [(o[1], o[2])] if o[0] == "G" else (o[2] if o[0] == "IF" else [])
        for n, q in gs:
            if n in ("t", "tdg"):
                t += 1
            elif n in ("cx", "cz"):
                cx += 1
    return t, cx


# condition helpers used by the a7 passes (merge, E-transform)
def pcond(s):
    return [tuple(sorted(set(g))) for g in parse_cond(s)]


def fcond(groups):
    return " && ".join("(" + " ^ ".join(f"c[{b}]" for b in g) + ")" for g in groups)
