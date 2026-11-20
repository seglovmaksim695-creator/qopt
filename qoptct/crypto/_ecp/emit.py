# SPDX-License-Identifier: Apache-2.0
"""a9/emit.py — вывод схем Circ в OpenQASM 3 (формат verify_qasm / verify_big) и в Feynman .qc (эталон: and/erase -> tof),
плюс независимый разреженный симулятор квантового состояния (sparse_check) для больших схем: суперпозиция K базисных входов
со случайными комплексными амплитудами, измерения с выборкой исходов, сравнение с классической функцией.
"""
import re, cmath, math, random

AND_G = lambda a, b, t: [("h", (t,)), ("t", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)),
                         ("cx", (b, t)), ("tdg", (t,)), ("h", (t,)), ("sdg", (t,))]
CCX_G = lambda a, b, t: [("h", (t,)), ("cx", (b, t)), ("tdg", (t,)), ("cx", (a, t)), ("t", (t,)), ("cx", (b, t)),
                         ("tdg", (t,)), ("cx", (a, t)), ("t", (b,)), ("t", (t,)), ("h", (t,)), ("cx", (a, b)),
                         ("t", (a,)), ("tdg", (b,)), ("cx", (a, b))]


def g(name, qs): return f"{name} " + ", ".join(f"q[{x}]" for x in qs) + ";"


def to_qasm(C, path=None):
    L = []; nb = 0
    for o in C.ops:
        k = o[0]
        if k == 'x': L.append(g('x', (o[1],)))
        elif k == 'cx': L.append(g('cx', o[1:]))
        elif k == 'ccx': L += [g(a, b) for a, b in CCX_G(*o[1:])]
        elif k == 'and': L += [g(a, b) for a, b in AND_G(*o[1:])]
        elif k == 'erase':
            t, a, b, lin = o[1:]
            L.append(g('h', (t,))); L.append(f"c[{nb}] = measure q[{t}];")
            body = [g('cz', (a, b))] + [g('z', (w,)) for w in lin] + [g('x', (t,))]
            L.append(f"if (c[{nb}]) {{" + " ".join(body) + "}"); nb += 1
    txt = 'OPENQASM 3;\ninclude "stdgates.inc";\n' + f"qubit[{C.nw}] q;\nbit[{max(nb, 1)}] c;\n" + "\n".join(L) + "\n"
    if path: open(path, 'w').write(txt)
    return txt


def count(txt):
    T = len(re.findall(r"\b(t|tdg) q", txt)); CX = len(re.findall(r"\b(cx|cz) q", txt))
    return T, CX


def to_qc(C, path, inputs=None):
    """эталон уровня Тоффоли: and -> tof, erase -> tof (+ cnot от lin). inputs: провода-данные (.i)."""
    nm = [f"q{i}" for i in range(C.nw)]
    ins = C.data if inputs is None else inputs
    L = [".v " + " ".join(nm), ".i " + " ".join(nm[i] for i in ins), "", "BEGIN"]
    for o in C.ops:
        k = o[0]
        if k == 'x': L.append(f"tof {nm[o[1]]}")
        elif k == 'cx': L.append(f"tof {nm[o[1]]} {nm[o[2]]}")
        elif k in ('ccx', 'and'): L.append(f"tof {nm[o[1]]} {nm[o[2]]} {nm[o[3]]}")
        elif k == 'erase':
            t, a, b, lin = o[1:]
            for w in lin: L.append(f"tof {nm[w]} {nm[t]}")
            L.append(f"tof {nm[a]} {nm[b]} {nm[t]}")
    L.append("END")
    open(path, 'w').write("\n".join(L) + "\n")


# ---------------- разреженный симулятор (независим от Circ: читает только текст qasm) ----------------
W8 = [cmath.exp(1j * math.pi * k / 4) for k in range(8)]
PH = {'z': 4, 's': 2, 'sdg': 6, 't': 1, 'tdg': 7}


def parse(txt):
    prog = []
    for L in txt.splitlines():
        L = L.strip()
        if not L or L.startswith(("OPENQASM", "include", "qubit", "bit")): continue
        m = re.match(r"c\[(\d+)\] = measure q\[(\d+)\];", L)
        if m: prog.append(("M", int(m.group(2)), int(m.group(1)))); continue
        m = re.match(r"if \(c\[(\d+)\]\) \{(.*)\}$", L)
        if m:
            body = [b.strip() for b in m.group(2).split(";") if b.strip()]
            prog.append(("IF", int(m.group(1)), [_gate(b) for b in body])); continue
        prog.append(("G",) + _gate(L))
    return prog


def _gate(s):
    s = s.rstrip(";"); name, args = s.split(None, 1)
    return (name, tuple(int(x) for x in re.findall(r"q\[(\d+)\]", args)))


def _apply(st, name, q):
    if name == 'x':
        b = 1 << q[0]; return {k ^ b: v for k, v in st.items()}
    if name == 'cx':
        c, t = 1 << q[0], 1 << q[1]; return {(k ^ t if k & c else k): v for k, v in st.items()}
    if name == 'cz':
        m = (1 << q[0]) | (1 << q[1]); return {k: (-v if k & m == m else v) for k, v in st.items()}
    if name in PH:
        b = 1 << q[0]; f = W8[PH[name]]; return {k: (v * f if k & b else v) for k, v in st.items()}
    if name == 'h':
        b = 1 << q[0]; r = 1 / math.sqrt(2); out = {}
        for k, v in st.items():
            k0 = k & ~b; k1 = k | b
            out[k0] = out.get(k0, 0) + v * r
            out[k1] = out.get(k1, 0) + (-v * r if k & b else v * r)
        return {k: v for k, v in out.items() if abs(v) > 1e-13}
    raise ValueError(name)


def run(prog, st, rng):
    c = {}
    for op in prog:
        if op[0] == "G": st = _apply(st, op[1], op[2])
        elif op[0] == "IF":
            if c.get(op[1], 0):
                for nm, q in op[2]: st = _apply(st, nm, q)
        else:
            _, q, k = op; b = 1 << q
            p1 = sum(abs(v) ** 2 for kk, v in st.items() if kk & b)
            m = int(rng.random() < p1); pm = p1 if m else 1 - p1
            st = {kk: v / math.sqrt(pm) for kk, v in st.items() if bool(kk & b) == bool(m)}; c[k] = m
    return st


def sparse_check(txt, encode, func, inputs, K=8, trials=3, seed=1):
    """encode(inp)->базисный индекс (int), func(inp)->ожидаемый выходной индекс (анциллы в 0).
    inputs: список допустимых входов. Возвращает худшую |<want|got>| (без перенормировки по анциллам)."""
    prog = parse(txt); rng = random.Random(seed); worst = 1.0
    for _ in range(trials):
        xs = rng.sample(inputs, min(K, len(inputs)))
        amps = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in xs]
        nr = math.sqrt(sum(abs(a) ** 2 for a in amps)); amps = [a / nr for a in amps]
        st = {}
        for x, a in zip(xs, amps): st[encode(x)] = st.get(encode(x), 0) + a
        out = run(prog, st, rng)
        want = {}
        for x, a in zip(xs, amps): want[func(x)] = want.get(func(x), 0) + a
        ov = abs(sum(v.conjugate() * out.get(k, 0) for k, v in want.items()))
        worst = min(worst, ov)
    return worst
