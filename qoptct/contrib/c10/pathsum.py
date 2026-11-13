"""c10: пути-интегральное (sum-over-paths) представление детерминированных схем Clifford+T с измерениями.
Схема: |0>/базисный вход, H на свежих проводах (|+>), CX/CZ/S/T/Z/X, пара "h q; measure q" (X-измерение), условные
поправки if (c[a] ^ c[b] && c[d]) {...}. Исходы измерений ЗАДАНЫ (ветвь), поэтому условия вычисляются сразу.
Состояние ветви: 2^{-nH/2} * sum_y w^{F(y)} |W y + e>, где w = exp(i pi/4), F — многочлен над Z8:
  F = gph + sum_mask lin[mask]*P(mask) + 4*sum_{(a,b)} P(a)P(b),  P(mask) = XOR переменных y в маске.
"""
import re, random

W8 = 8


def parse_cond(s):
    s = s.strip()
    assert re.fullmatch(r"[\s()c\[\]\d&^!]+", s), "неизвестное условие: " + s
    groups = []
    for g in s.split("&&"):
        g = g.strip()
        while g.startswith("(") and g.endswith(")") and g.count("(") == 1:
            g = g[1:-1].strip()
        g = g.replace("(", "").replace(")", "")
        assert "!" not in g
        groups.append([int(x) for x in re.findall(r"c\[(\d+)\]", g)])
    return groups


def parse(path):
    """-> (N, nbits, ops). ops: ('g',name,qs) | ('h',q) | ('m',q,k) | ('if',groups,[(name,qs)...])"""
    txt = open(path).read()
    N = int(re.search(r"qubit\[(\d+)\]", txt).group(1))
    mb = re.search(r"(?<!qu)bit\[(\d+)\]", txt)
    nb = int(mb.group(1)) if mb else 0
    ops = []
    for L in txt.splitlines():
        L = L.strip()
        if not L or L.startswith(("OPENQASM", "include", "qubit", "bit", "//")):
            continue
        m = re.match(r"c\[(\d+)\] = measure q\[(\d+)\];", L)
        if m:
            ops.append(("m", int(m.group(2)), int(m.group(1)))); continue
        m = re.match(r"if \((.*)\) \{(.*)\}$", L)
        if m:
            body = [b.strip() for b in m.group(2).split(";") if b.strip()]
            ops.append(("if", parse_cond(m.group(1)), [_gate(b) for b in body])); continue
        name, qs = _gate(L.rstrip(";"))
        if name == "h":
            ops.append(("h", qs[0]))
        else:
            ops.append(("g", name, qs))
    return N, nb, ops


def _gate(s):
    s = s.rstrip(";").strip(); name, args = s.split(None, 1)
    return (name, tuple(int(x) for x in re.findall(r"q\[(\d+)\]", args)))


class PS:
    """Построитель пути-интеграла. wires[q] = (mask, const)."""

    def __init__(self, N, init_bits=None):
        self.N = N
        init_bits = init_bits or {}
        self.w = [(0, init_bits.get(q, 0)) for q in range(N)]
        self.k = 0
        self.lin = {}      # mask -> coef mod 8
        self.quad = {}     # (ma, mb) -> 1  (xor-накопление), коэффициент 4
        self.gph = 0
        self.nH = 0
        self.nmeas = 0
        self.alive = [True] * N   # False после измерения
        self.consistent = True
        self.cons = []   # (mask, bit): P(mask) = bit  (Z-измерения проводов с переменными)

    # ---- элементарные
    def add_lin(self, mask, c, const=0):
        """w^{c * (const XOR P(mask))}"""
        c %= 8
        if mask == 0:
            self.gph = (self.gph + c * const) % 8; return
        if const:
            self.gph = (self.gph + c) % 8; c = (-c) % 8
        v = (self.lin.get(mask, 0) + c) % 8
        if v: self.lin[mask] = v
        else: self.lin.pop(mask, None)

    def add_quad(self, wa, wb):
        """(-1)^{va*vb}"""
        (ma, ca), (mb, cb) = wa, wb
        # (ca^Pa)(cb^Pb) = ca cb ^ ca Pb ^ cb Pa ^ Pa Pb  (над Z2, т.к. знак (-1)^.)
        if ca and cb: self.gph = (self.gph + 4) % 8
        if ca: self.add_lin(mb, 4)
        if cb: self.add_lin(ma, 4)
        if ma == 0 or mb == 0:
            return
        if ma == mb:
            self.add_lin(ma, 4); return   # P*P = P
        key = (ma, mb) if ma < mb else (mb, ma)
        if key in self.quad: del self.quad[key]
        else: self.quad[key] = 1

    # ---- вентили
    def gate(self, name, qs):
        w = self.w
        if name == "cx":
            c, t = qs; w[t] = (w[t][0] ^ w[c][0], w[t][1] ^ w[c][1])
        elif name == "x":
            q = qs[0]; w[q] = (w[q][0], w[q][1] ^ 1)
        elif name in ("t", "tdg", "s", "sdg", "z"):
            c = {"t": 1, "tdg": 7, "s": 2, "sdg": 6, "z": 4}[name]; m, cn = w[qs[0]]; self.add_lin(m, c, cn)
        elif name == "cz":
            self.add_quad(w[qs[0]], w[qs[1]])
        elif name == "reset":
            q = qs[0]
            if w[q][0] != 0: raise ValueError("reset провода с переменными (не измеренного) не поддержан")
            w[q] = (0, 0); self.alive[q] = True
        elif name == "swap":
            a, b = qs; w[a], w[b] = w[b], w[a]
        else:
            raise ValueError("вентиль не поддержан: " + name)

    def hadamard(self, q):
        """H на проводе: новая переменная."""
        m, c = self.w[q]
        v = 1 << self.k; self.k += 1; self.nH += 1
        if m == 0 and c == 0:
            self.w[q] = (v, 0); return
        # sum_v (-1)^{w v}: фаза между P(w)^c и v
        self.add_quad((m, c), (v, 0))
        self.w[q] = (v, 0)

    def measure_after_h(self, q, outcome):
        """h q; measure q с исходом outcome: множитель (-1)^{outcome * w}, провод -> константа outcome."""
        m, c = self.w[q]
        self.nH += 1; self.nmeas += 1
        if outcome: self.add_lin(m, 4, c)
        self.w[q] = (0, outcome); self.alive[q] = False

    def measure_z(self, q, outcome):
        """Измерение в Z (без H): ограничение w = outcome. Только константные провода."""
        m, c = self.w[q]
        if m != 0:
            self.cons.append((m, outcome ^ c))
        elif c != outcome: self.consistent = False
        self.w[q] = (0, outcome); self.alive[q] = False


def build(prog, init_bits, outcomes=None, rng=None, ref_check=False):
    """prog = (N, nb, ops). outcomes: список из nb битов (если None — случайно). Возвращает (PS, c)."""
    N, nb, ops = prog
    rng = rng or random.Random(0)
    c = [0] * nb
    ps = PS(N, init_bits)
    if outcomes is None: outcomes = [rng.randint(0, 1) for _ in range(nb)]
    i = 0
    while i < len(ops):
        op = ops[i]
        if op[0] == "h":
            q = op[1]
            nxt = ops[i + 1] if i + 1 < len(ops) else None
            if nxt and nxt[0] == "m" and nxt[1] == q:
                k = nxt[2]; c[k] = outcomes[k]; ps.measure_after_h(q, outcomes[k]); i += 2; continue
            ps.hadamard(q)
        elif op[0] == "g":
            ps.gate(op[1], op[2])
        elif op[0] == "if":
            val = all((sum(c[b] for b in g) & 1) for g in op[1])
            if val:
                for name, qs in op[2]: ps.gate(name, qs)
        elif op[0] == "m":
            k = op[2]; c[k] = outcomes[k]; ps.measure_z(op[1], outcomes[k])
        i += 1
    return ps, c


def stats(ps):
    import collections
    odd = [m for m, c in ps.lin.items() if c & 1]
    wt = collections.Counter(bin(m).count("1") for m in ps.lin)
    return dict(k=ps.k, nH=ps.nH, lin=len(ps.lin), odd=len(odd), quad=len(ps.quad),
                wt=dict(sorted(wt.items())))


if __name__ == "__main__":
    import sys
    pr = parse(sys.argv[1])
    ps, c = build(pr, {})
    print(stats(ps))
    print("alive wires", [q for q in range(ps.N) if ps.alive[q]])
    for q in range(ps.N):
        if ps.alive[q]:
            print(q, bin(ps.w[q][0]).count("1"), ps.w[q][1])
