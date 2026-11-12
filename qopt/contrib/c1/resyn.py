"""c1/resyn.py: CX re-synthesis of the unconditional phase/linear pieces of a measured circuit (e.g. agents/c2/out/*.qasm).

The qasm is cut at every non-unitary line (measure, if-blocks, ...).  Inside each unitary segment the gates {cx,x,cz,t,tdg,s,
sdg,z} are grouped into maximal H-free blocks with commutation (a gate is pulled across an H on a different wire; layers()
from the phase-polynomial lowering program).  Each block = phase polynomial f over its input wires + affine map (W, K).  It is
re-synthesised by the c1 synthesiser: T-parity network (graysynth / end-aware greedy, many seeds) + linear fix + FREE Clifford
absorption (S on any parity present on a wire at any time; CZ only for the remainder) + X for constants.  A block is replaced
only if cx+cz strictly drops, T does not grow, and the block re-extracts to exactly the same (canonical f, W, K).
CLI: python resyn.py in.qasm qc_path out.qasm [seeds]"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, re, random
sys.path.insert(0, (_QPKG + "/contrib/c1"))
from csyn import *
layers = pp.layers

GATES = {"h", "x", "cx", "cz", "t", "tdg", "s", "sdg", "z"}


def parse(path):
    head, items = [], []
    nq = None
    for L in open(path):
        s = L.strip()
        if not s: continue
        m = re.match(r"qubit\[(\d+)\]", s)
        if s.startswith(("OPENQASM", "include")) or re.match(r"(qubit|bit)\[", s):
            head.append(s)
            if m: nq = int(m.group(1))
            continue
        m = re.match(r"^(\w+)\s+(q\[\d+\](?:\s*,\s*q\[\d+\])*)\s*;$", s)
        if m and m.group(1) in GATES:
            items.append(("g", m.group(1), tuple(int(x) for x in re.findall(r"q\[(\d+)\]", m.group(2)))))
        else:
            items.append(("raw", s))
    return head, items, nq


def block_cost(gs):
    return sum(1 for g, _ in gs if g in ("cx", "cz"))


def resynth_block(gs, nq, seeds=60, rng_base=0):
    """gs: list of (g, qs) H-free.  returns better list or None."""
    if block_cost(gs) == 0: return None
    wires = sorted({q for _, qs in gs for q in qs})
    loc = {w: i for i, w in enumerate(wires)}; N = len(wires)
    lg = [(g, tuple(loc[q] for q in qs)) for g, qs in gs]
    W = [1 << i for i in range(N)]; K = [0] * N; P = {}
    pp.block_poly(lg, W, K, P)
    P = {m: w % 8 for m, w in P.items() if w % 8}
    canP = canon(P)
    PT = {m: 1 for m in t_parities(P)}
    T0 = sum(1 for g, _ in gs if g in ("t", "tdg"))
    if len(PT) > T0: return None
    best = None; Wi = [1 << i for i in range(N)]
    for s in range(seeds):
        rng = random.Random(rng_base + s)
        net, Wn = net_variants(PT, Wi, W, N, rng, s) if PT else ([], list(Wi))
        if net is None: continue
        M = [pp.solve(W, N, Wn[r]) for r in range(N)]
        ops = pp.lin_best(M, N, range(N), rng, 6)
        body = clifford_absorb(list(net) + [("cx", o) for o in ops], Wi, P, N, canP=canP)
        if body is None: continue
        body = body + [("x", (q,)) for q in range(N) if K[q]]
        c = cost(body)
        if best is None or c < best[0]: best = (c, body)
        if not PT and s >= 8: break
    if best is None or best[0] >= block_cost(gs): return None
    body = best[1]
    # exact algebraic check
    W2 = [1 << i for i in range(N)]; K2 = [0] * N; P2 = {}
    pp.block_poly(body, W2, K2, P2)
    if canon(P2) != canP or W2 != W or K2 != K: return None
    if sum(1 for g, _ in body if g in ("t", "tdg")) > T0: return None
    return [(g, tuple(wires[q] for q in qs)) for g, qs in body]


def run(path_in, path_out, seeds=60, log=print):
    head, items, nq = parse(path_in)
    out = []; seg = []; stats = dict(blocks=0, improved=0, before=0, after=0)
    def flush():
        nonlocal seg
        if not seg: return
        for kind, x in layers(seg):
            if kind == "H":
                out.extend(("g", "h", (q,)) for q in x)
            else:
                stats["blocks"] += 1; stats["before"] += block_cost(x)
                nb = resynth_block(x, nq, seeds, rng_base=stats["blocks"] * 1000)
                if nb is not None:
                    stats["improved"] += 1; x = nb
                stats["after"] += block_cost(x)
                out.extend(("g", g, qs) for g, qs in x)
        seg = []
    for it in items:
        if it[0] == "g": seg.append((it[1], it[2]))
        else: flush(); out.append(it)
    flush()
    lines = list(head)
    for it in out:
        if it[0] == "raw": lines.append(it[1])
        else: lines.append(f"{it[1]} " + ", ".join(f"q[{q}]" for q in it[2]) + ";")
    open(path_out, "w").write("\n".join(lines) + "\n")
    log(str(stats))
    return stats


def counts(path):
    T = CX = 0
    for L in open(path):
        CX += len(re.findall(r"\b(cx|cz)\s+q", L)); T += len(re.findall(r"\b(t|tdg)\s+q", L))
    return T, CX


if __name__ == "__main__":
    import verify_qasm
    a = sys.argv
    run(a[1], a[3], int(a[4]) if len(a) > 4 else 60)
    print("before", counts(a[1]), "after", counts(a[3]))
    print("F", verify_qasm.check(a[3], a[2]))
