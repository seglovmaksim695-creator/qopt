"""c8/hyb.py: hybrid compiler. Classical reversible blocks (Toffoli/CNOT/X, possibly with H-sandwiches) -> c2 (AND + X-measurement);
non-classical gates and the H gates that cannot be absorbed into frames stay REAL and cut the circuit into blocks.
Blocks are compiled independently (ancillas start/end in |0>), Pareto-combined (T, CX, extra) and spliced as OpenQASM 3 text."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, random, itertools, re, json
sys.path[:0] = [_QPKG, (_QPKG + '/portfolio'), (_QPKG + '/contrib/c2')]
import qopt, evo_lib, andc, sweep2

NONCL = ('t', 'tdg', 's', 'sdg', 'h')


class Need(Exception):
    def __init__(self, qs): self.qs = set(qs)


def build(hl, n, real):
    """-> list of ('L'|'R', g, q).  Raises Need(set of H indices to make real)."""
    fr = [0] * n; src = [None] * n; last = [-1] * n; out = []
    def need(qs): raise Need({src[q] for q in qs if src[q] is not None})
    for i, (g, q) in enumerate(hl):
        q = tuple(q)
        if g == 'h':
            if i in real: out.append(('R', 'h', q))
            else: fr[q[0]] ^= 1; src[q[0]] = i; last[q[0]] = i
            continue
        if g in ('ccz', 'cz', 'z'):
            F = [x for x in q if fr[x]]
            if len(F) <= 1:
                if not F: out.append(('L', g, q))
                else:
                    t = F[0]; c = tuple(x for x in q if x != t)
                    out.append(('L', {'ccz': 'ccx', 'cz': 'cx', 'z': 'x'}[g], c + (t,)))
            else:
                keep = max(F, key=lambda x: last[x]); need([x for x in F if x != keep])
        elif g == 'x':
            out.append(('L', 'z' if fr[q[0]] else 'x', q))
        elif g in ('cx', 'ccx'):
            if any(fr[x] for x in q[:-1]): need([x for x in q[:-1] if fr[x]])
            out.append(('L', ({'cx': 'cz', 'ccx': 'ccz'}[g]) if fr[q[-1]] else g, q))
        elif g in NONCL:
            F = [x for x in q if fr[x]]
            if F: need(F)
            out.append(('R', g, q))
        else: raise ValueError(g)
    if any(fr): need([x for x in range(n) if fr[x]])
    return out


def segment(hl, n):
    real = set()
    while True:
        try: return build(hl, n, real), real
        except Need as e:
            if not e.qs - real: raise RuntimeError("no progress")
            real |= e.qs


def blocks_of(items):
    seq = []; cur = []
    for k, g, q in items:
        if k == 'L': cur.append((g, q))
        else:
            if cur: seq.append(('B', cur)); cur = []
            seq.append(('R', (g, q)))
    if cur: seq.append(('B', cur))
    return seq


def block_cands(lg0, n, lim, variants=12, seed=1):
    rng = random.Random(seed); out = {}
    for v in range(variants):
        lg = lg0 if (v == 0 or len(lg0) < 3) else sweep2.shuffle(lg0, rng)
        for lazy, lm, ma in itertools.product((True, False, 0.2, 0.4, 0.6), (3, 1), sorted({lim, 4, 3, 2, 1, 0})):
            if ma > lim: continue
            try:
                C = andc.Compiler(n, lg, ma, random.Random(v), lazy=lazy, linmax=lm)
                ops, N = andc.compact(C.run(), n)
            except Exception: continue
            for fold in (False, True):
                txt = andc.to_qasm(ops, N, fold=fold); T, CX = andc.count(txt)
                k = (T, CX, N - n)
                if k not in out: out[k] = txt
    return out


def prune(pts):
    pts = sorted(pts, key=lambda p: p[:3]); res = []
    for p in pts:
        if not any(all(a <= b for a, b in zip(o[:3], p[:3])) for o in res): res.append(p)
    return res


def body_lines(txt):
    L = txt.split("\n"); i = max(k for k, l in enumerate(L) if l.startswith(("qubit", "bit")))
    return [l for l in L[i + 1:] if l.strip()]


def shift_bits(lines, off):
    return [re.sub(r"c\[(\d+)\]", lambda m: f"c[{int(m.group(1)) + off}]", l) for l in lines]


def compile_hybrid(path, variants=12, max_total=22, log=print):
    names, hl = qopt.load_qc(path); n = len(names); lim = max_total - n
    items, real = segment(hl, n)
    seq = blocks_of(items)
    log(f"{len(real)} real H, {len(seq)} segments, blocks {[len(b) for k, b in seq if k=='B']}")
    realT = sum(1 for k, g in seq if k == 'R' and g[0] in ('t', 'tdg'))
    parts = []
    for k, b in seq:
        if k == 'R': parts.append(('R', b)); continue
        parts.append(('B', block_cands(b, n, lim, variants)))
    # Pareto DP
    front = [(realT, 0, 0, [])]
    for k, b in parts:
        if k == 'R': continue
        new = []
        for (T, CX, e, ch) in front:
            for (t, c, ex), txt in b.items(): new.append((T + t, CX + c, max(e, ex), ch + [(t, c, ex)]))
        front = prune(new)[:400]
    return n, parts, front


def emit(n, parts, choice, path=None, fold=True):
    ex = max(c[2] for c in choice) if choice else 0; N = n + ex
    lines = []; nb = 0; bi = 0
    for k, b in parts:
        if k == 'R':
            g, q = b; lines.append(f"{g} " + ", ".join(f"q[{x}]" for x in q) + ";")
        else:
            txt = b[choice[bi]]; bi += 1
            bl = body_lines(txt); cnt = len(re.findall(r"= measure", txt))
            lines += shift_bits(bl, nb); nb += cnt
    if fold: lines = fold_lines(lines, N)
    txt = "OPENQASM 3;\ninclude \"stdgates.inc\";\n" + f"qubit[{N}] q;\nbit[{max(nb,1)}] c;\n" + "\n".join(lines) + "\n"
    if path: open(path, "w").write(txt)
    return txt, N


PLAIN = re.compile(r"^(h|x|cx|t|tdg|s|sdg|z) (q\[\d+\](?:, q\[\d+\])?);$")


def fold_lines(lines, N):
    out = []; run = []
    def flush():
        if run:
            for g, q in evo_lib.fold_cancel(list(run), N): out.append(f"{g} " + ", ".join(f"q[{x}]" for x in q) + ";")
        run.clear()
    for l in lines:
        m = PLAIN.match(l)
        if m: run.append((m.group(1), tuple(int(x) for x in re.findall(r"\d+", m.group(2)))))
        else: flush(); out.append(l)
    flush()
    return out


def count_txt(txt):
    return andc.count(txt)
