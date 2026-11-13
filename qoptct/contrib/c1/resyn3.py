"""c1/resyn3.py: LAZY LINEAR FRAME re-synthesis of a measured circuit (c2 outputs).

Logical circuit = original qasm.  Physical wires differ from logical ones by an invertible linear frame F (+ constants) that is
carried forward instead of being undone.  Barriers: H gates and every non-unitary line (measure, if-blocks); they need F to
be trivial ONLY on the wires they touch (row q exact, q not mixed into other rows, equal constant) -> partial linear fix.
Between barriers (maximal blocks with commutation), the block's phase polynomial (over global path variables) is synthesised
from the CURRENT physical rows: T-parity network (graysynth / end-aware greedy) + partial fix + FREE Clifford absorption
(S on any parity present on a wire during the block, CZ only for the remainder) + constant-aware emission.  At the end the
full frame is removed.  Every block is checked algebraically (phases + physical rows); the whole circuit by verify_qasm.
CLI: python resyn3.py in.qasm qc out.qasm [variants]"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, re, random
sys.path.insert(0, (_QPKG + "/contrib/c1"))
from resyn import parse, counts
from csyn import *

DIAG = {"t", "tdg", "s", "sdg", "z", "cz"}
NEG = {"t": "tdg", "tdg": "t", "s": "sdg", "sdg": "s", "z": "z"}


def raw_wires(s):
    return tuple(sorted({int(q) for q in re.findall(r"q\[(\d+)\]", s)}))


def comm2(x, y):
    if x[0] == "raw" or y[0] == "raw":
        if x[0] == "raw" and y[0] == "raw": return False
        return not (set(x[2]) & set(y[2]))
    return pp.comm((x[1], x[2]), (y[1], y[2]))


def layers2(ops):
    R = list(ops); out = []
    while R:
        for want_bar in (False, True):
            rest, took = [], []
            for x in R:
                isbar = x[0] == "raw" or x[1] == "h"
                if isbar == want_bar and all(comm2(x, y) for y in rest): took.append(x)
                else: rest.append(x)
            R = rest
            if took: out.append(("X" if want_bar else "B", took))
    return out


def absorb_general(body0, Wstart, Ptarget, NW, NV, canP):
    """clifford_absorb with separate wire count NW and variable count NV (indexing)."""
    P0 = {}; W = list(Wstart); K = [0] * NW
    pp.block_poly(body0, W, K, P0)
    R = sub_canon(canP, canon(P0))
    B = target_B(R, NV)
    if B is None: return None
    seen = {}; W = list(Wstart)
    for q in range(NW):
        if W[q] and W[q] not in seen: seen[W[q]] = (0, q)
    snaps = [(0, list(W))]
    for idx, (g, qs) in enumerate(body0):
        if g == "cx":
            W[qs[1]] ^= W[qs[0]]; t = qs[1]
            if W[t] not in seen: seen[W[t]] = (idx + 1, t)
            snaps.append((idx + 1, list(W)))
    gens = list(seen.items()); E = Elim()
    for gi, (l, _) in enumerate(gens): E.add(llT(l, NV), 1 << gi)
    r, comb = E.full_reduce(B); czs = []
    if r:
        red = {}
        for pos, Ws in snaps:
            for a in range(NW):
                for b in range(a + 1, NW):
                    v, _ = E.full_reduce(abT(Ws[a], Ws[b], NV))
                    if v and v not in red: red[v] = (pos, a, b, Ws[a], Ws[b])
        if r in red: czs = [red[r]]
        else:
            for v in list(red):
                if (r ^ v) in red: czs = [red[v], red[r ^ v]]; break
            if not czs:
                E2 = Elim(); lst = list(red.items())
                for j, (v, _) in enumerate(lst): E2.add(v, 1 << j)
                rr, cc = E2.reduce(r)
                if rr: return None
                czs = [lst[j][1] for j in bitsof(cc)]
                czs = shrink_cz(r, czs, red, E, NV)
        B2 = B
        for c in czs: B2 ^= abT(c[3], c[4], NV)
        r2, comb = E.reduce(B2)
        if r2: return None
    ins = {}
    for gi in bitsof(comb):
        l, (pos, w) = gens[gi]; ins.setdefault(pos, []).append(("s", (w,)))
    for (pos, a, b, _, _) in czs: ins.setdefault(pos, []).append(("cz", (a, b)))
    body = []
    for idx in range(len(body0) + 1):
        body += ins.get(idx, [])
        if idx < len(body0): body.append(body0[idx])
    P1 = {}; W = list(Wstart); K = [0] * NW
    pp.block_poly(body, W, K, P1)
    D = sub_canon(canP, canon(P1)); z = 0
    for k, v in D.items():
        if len(k) != 1 or v != 4: return None
        z |= 1 << k[0]
    pre = []
    if z:
        c = pp.solve(Wstart, NW, z)
        if c is None: return None
        pre = [("z", (q,)) for q in bitsof(c)]
    return merge_phases(pre + body)


def const_emit(body, KW):
    """physical constants KW (modified in place): phase gates on a wire with constant 1 are negated; cz gets Z fixes."""
    out = []
    for g, qs in body:
        if g == "cx":
            KW[qs[1]] ^= KW[qs[0]]; out.append((g, qs))
        elif g == "x":
            KW[qs[0]] ^= 1; out.append((g, qs))
        elif g == "cz":
            a, b = qs; out.append((g, qs))
            if KW[a]: out.append(("z", (b,)))
            if KW[b]: out.append(("z", (a,)))
        elif g in NEG:
            out.append((NEG[g] if KW[qs[0]] else g, qs))
        else:
            out.append((g, qs))
    return out


def fix_ops(W, L, NW, Q, rng, tries=6):
    M = [pp.solve(L, NW, W[r]) for r in range(NW)]
    return pp.lin_best(M, NW, Q, rng, tries)


LAM = [0.0]


def run3(path_in, qc, path_out, variants=24, seed=0, log=print):
    head, items, nq = parse(path_in)
    NW = nq
    ops = []
    for it in items:
        ops.append(("g", it[1], it[2]) if it[0] == "g" else ("raw", it[1], raw_wires(it[1])))
    lay = layers2(ops)
    rng = random.Random(seed)
    L = [1 << q for q in range(NW)]; KL = [0] * NW; W = list(L); KW = [0] * NW; nv = NW
    out = []; stats = dict(blocks=0, fallback=0)
    def emit(gs):
        for g, qs in gs:
            if g == "cx": W[qs[1]] ^= W[qs[0]]
            out.append(("g", g, qs))
    def fix(Q):
        ops_ = fix_ops(W, L, NW, Q, rng)
        emit(const_emit([("cx", o) for o in ops_], KW))
        for q in Q:
            if KW[q] != KL[q]: emit([("x", (q,))]); KW[q] ^= 1
    for li, (kind, x) in enumerate(lay):
        if kind == "B":
            gs = [(o[1], o[2]) for o in x]
            P = {}; Lb = list(L); KLb = list(KL)
            pp.block_poly(gs, L, KL, P)
            P = {m: w % 8 for m, w in P.items() if w % 8}
            nxt = lay[li + 1][1] if li + 1 < len(lay) else None
            Q = sorted({q for o in nxt for q in o[2]}) if nxt else list(range(NW))
            canP = canon(P); PT = {m: 1 for m in t_parities(P)}
            T0 = sum(1 for g, _ in gs if g in ("t", "tdg"))
            best = None
            for v in range(variants):
                r2 = random.Random(seed * 100003 + li * 1009 + v)
                W2 = list(W)
                if PT:
                    try:
                        if v % 3 == 0: net = pp.graysynth(dict(PT), W2, NW, r2)
                        else: net = pp.greedy(dict(PT), W2, NW, r2, L, r2.uniform(0, 1.5), r2.uniform(0, 1.0), r2.uniform(0.1, 1.0))
                    except Exception:
                        net = None
                    if net is None: continue
                else:
                    net = []
                M = [pp.solve(L, NW, W2[r]) for r in range(NW)]
                fo = pp.lin_best(M, NW, Q, r2, 4)
                body0 = list(net) + [("cx", o) for o in fo]
                if P:
                    body = absorb_general(body0, W, P, NW, nv, canP)
                    if body is None: continue
                else:
                    body = body0
                if sum(1 for g, _ in body if g in ("t", "tdg")) > T0: continue
                c = cost(body)
                if LAM[0]:
                    Wd = list(W)
                    for g, qs in body:
                        if g == "cx": Wd[qs[1]] ^= Wd[qs[0]]
                    c += LAM[0] * len(fix_ops(Wd, L, NW, list(range(NW)), r2, 2))
                if best is None or c < best[0]: best = (c, body)
                if not PT and v >= 4: break
            ok = False
            if best is not None and LAM[0]:
                # alternative: remove the frame now and emit the original block (debt 0)
                alt = len(fix_ops(W, Lb, NW, list(range(NW)), rng, 4)) + sum(1 for g, _ in gs if g in ("cx", "cz"))
                if alt < best[0]: best = None
            if best is not None:
                body = best[1]
                # algebraic check from the physical start rows
                P1 = {}; Wc = list(W); Kc = [0] * NW
                pp.block_poly(body, Wc, Kc, P1)
                if canon(P1) == canP:
                    ok = True
            if ok:
                stats["blocks"] += 1
                emit(const_emit(body, KW))
            else:
                stats["fallback"] += 1
                L2, KL2 = list(L), list(KL)
                L[:] = Lb; KL[:] = KLb
                fix(list(range(NW)))
                L[:] = L2; KL[:] = KL2
                for g, qs in gs:
                    if g == "cx": KW[qs[1]] ^= KW[qs[0]]
                    elif g == "x": KW[qs[0]] ^= 1
                emit(gs)
            # the constants of the logical side changed by x gates inside the block: physical KW follows emitted gates
            continue
        # barrier group
        Q = sorted({q for o in x for q in o[2]})
        fix(Q)
        for o in x:
            if o[0] == "raw":
                out.append(("raw", o[1]))
                for q in o[2]:
                    if re.search(rf"\bx q\[{q}\]", o[1]) or re.search(rf"measure q\[{q}\]", o[1]):
                        L[q] = W[q] = 1 << nv; nv += 1; KL[q] = KW[q] = 0
            else:
                q = o[2][0]
                out.append(("g", "h", (q,)))
                L[q] = W[q] = 1 << nv; nv += 1; KL[q] = KW[q] = 0
    fix(list(range(NW)))
    assert W == L and KW == KL
    lines = list(head)
    for it in out:
        lines.append(it[1] if it[0] == "raw" else f"{it[1]} " + ", ".join(f"q[{q}]" for q in it[2]) + ";")
    open(path_out, "w").write("\n".join(lines) + "\n")
    log(str(stats))
    return stats


if __name__ == "__main__":
    import verify_qasm
    a = sys.argv
    best = None
    for lam in (0.0, 0.5, 1.0):
        for sd in range(int(a[5]) if len(a) > 5 else 2):
            LAM[0] = lam
            tmp = a[3] + ".tmp"
            run3(a[1], a[2], tmp, int(a[4]) if len(a) > 4 else 24, seed=sd, log=lambda s: None)
            c = counts(tmp)
            print("lam", lam, "seed", sd, c, flush=True)
            if best is None or c < best[0]:
                best = (c, open(tmp).read())
    open(a[3], "w").write(best[1])
    print("before", counts(a[1]), "after", counts(a[3]))
    print("F", verify_qasm.check(a[3], a[2]))
