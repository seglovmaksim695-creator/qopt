"""c1/resyn2.py: resyn.py + two measurement-related freedoms of a block (taken from the tg3 gauge analysis):
 F1 (start): a wire that is in |+> at the block start (|0> + H; |0> = fresh ancilla or measured+reset wire) is invariant under
     any CX that targets it, so the block may start from rows Winit[a] = e_a + (other rows)  (the CX are never executed);
 F2 (end): a wire a that is X-measured right after the block (only its H in between) may end with W[a] + sum W[c]; the missing
     CX(c->a) before an X-measurement equal a conditional Z on c after the measurement (added as `if (c[k]) {z q[c];}`).
     c must not get an H in the layer between the block and the measurement.
Each block: many random (Winit, Wend) gauges x synthesis seeds; replaced only if cx+cz drops; exact algebraic check of
(f, W, K) with the gauges; the whole circuit is then checked by verify_qasm (strict)."""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, re, random
sys.path.insert(0, (_QPKG + "/contrib/c1"))
from resyn import *


def synth_gauge(gs, free_start, free_end, seeds, variants, rng_base=0):
    wires = sorted({q for _, qs in gs for q in qs} | set(free_end))
    loc = {w: i for i, w in enumerate(wires)}; N = len(wires)
    lg = [(g, tuple(loc[q] for q in qs)) for g, qs in gs]
    W = [1 << i for i in range(N)]; K = [0] * N; P = {}
    pp.block_poly(lg, W, K, P)
    P = {m: w % 8 for m, w in P.items() if w % 8}
    canP = canon(P); PT = {m: 1 for m in t_parities(P)}
    T0 = sum(1 for g, _ in gs if g in ("t", "tdg"))
    if len(PT) > T0: return None
    fs = [loc[a] for a in free_start if a in loc]
    fe = {loc[a]: [loc[c] for c in cs if c in loc and c != a] for a, cs in free_end.items() if a in loc}
    base = block_cost(gs); best = None
    rng0 = random.Random(rng_base)
    for v in range(variants):
        Wi = [1 << i for i in range(N)]; We = list(W); extra = {}
        if v > 0:
            for a in fs:
                if rng0.random() < 0.7:
                    for c in rng0.sample(range(N), min(N, rng0.randint(1, 2))):
                        if c != a: Wi[a] ^= 1 << c
            for a, cs in fe.items():
                if cs and rng0.random() < 0.7:
                    S = set(rng0.sample(cs, min(len(cs), rng0.randint(1, 2))))
                    for c in S: We[a] ^= W[c]
                    if S: extra[a] = S
            if not rank_ok(Wi, N) or not rank_ok(We, N): continue
        for s in range(seeds):
            rng = random.Random(rng_base + 7919 * v + s)
            net, Wn = net_variants(PT, Wi, We, N, rng, s) if PT else ([], list(Wi))
            if net is None: continue
            M = [pp.solve(We, N, Wn[r]) for r in range(N)]
            ops = pp.lin_best(M, N, range(N), rng, 6)
            body = clifford_absorb(list(net) + [("cx", o) for o in ops], Wi, P, N, canP=canP)
            if body is None: continue
            body = body + [("x", (q,)) for q in range(N) if K[q]]
            c = cost(body)
            if best is None or c < best[0]: best = (c, body, Wi, We, dict(extra))
            if not PT and s >= 6: break
    if best is None or best[0] >= base: return None
    c, body, Wi, We, extra = best
    pre = mat_cx(Wi, N)
    c1, W2, K2 = extract(pre + [g for g in body], N)
    if c1 != canP or W2 != We or K2 != K: return None
    if sum(1 for g, _ in body if g in ("t", "tdg")) > T0: return None
    return [(g, tuple(wires[q] for q in qs)) for g, qs in body], {wires[a]: [wires[c] for c in S] for a, S in extra.items()}


def run2(path_in, qc, path_out, seeds=20, variants=30, log=print):
    import qopt
    nd = len(qopt.load_qc(qc)[0])
    head, items, nq = parse(path_in)
    state = {q: ("zero" if q >= nd else "other") for q in range(nq)}
    out = []; seg = []; stats = dict(blocks=0, improved=0, before=0, after=0, condz=0)
    reset_ok = [None]
    pending_z = {}     # wire a -> list of c  (to add after the measurement of a)
    def flush(next_raw):
        nonlocal seg
        if not seg: return
        lay = layers(seg)
        meas_a = None
        if next_raw is not None:
            m = re.match(r"c\[(\d+)\]\s*=\s*measure\s+q\[(\d+)\]", next_raw)
            if m: meas_a = int(m.group(2))
        lastB = max((i for i, (k, _) in enumerate(lay) if k == "B"), default=-1)
        for i, (kind, x) in enumerate(lay):
            if kind == "H":
                for q in x:
                    state[q] = {"zero": "plus", "plus": "zero"}.get(state[q], "other")
                out.extend(("g", "h", (q,)) for q in x)
                continue
            fstart = [q for q in sorted({q for _, qs in x for q in qs}) if state[q] == "plus"]
            fend = {}
            if i == lastB and meas_a is not None and i + 1 < len(lay) and lay[i + 1][0] == "H" and meas_a in lay[i + 1][1] \
                    and i + 2 == len(lay):
                hset = set(lay[i + 1][1])
                fend = {meas_a: [q for q in sorted({q for _, qs in x for q in qs}) if q not in hset and q != meas_a]}
            stats["blocks"] += 1; stats["before"] += block_cost(x)
            r = synth_gauge(x, fstart, fend, seeds, variants, rng_base=stats["blocks"] * 1000)
            if r is not None:
                nb, extra = r; stats["improved"] += 1; x = nb
                for a, cs in extra.items():
                    pending_z.setdefault(a, []).extend(cs); stats["condz"] += len(cs)
            stats["after"] += block_cost(x)
            for _, qs in x:
                for q in qs: state[q] = "other"
            out.extend(("g", g, qs) for g, qs in x)
        seg = []
    for idx, it in enumerate(items):
        if it[0] == "g": seg.append((it[1], it[2])); continue
        flush(it[1]); out.append(it)
        m = re.match(r"c\[(\d+)\]\s*=\s*measure\s+q\[(\d+)\]", it[1])
        if m:
            k, a = int(m.group(1)), int(m.group(2))
            for c in pending_z.pop(a, []):
                out.append(("raw", f"if (c[{k}]) {{z q[{c}];}}"))
            nxt = items[idx + 1] if idx + 1 < len(items) else None
            state[a] = "zero" if (nxt and nxt[0] == "raw" and re.search(rf"c\[{k}\]\)", nxt[1]) and f"x q[{a}];" in nxt[1]) else "other"
            reset_ok[0] = a if state[a] == "zero" else None
        else:
            # any other raw line: wires mentioned become 'other' (conservative), except the reset of the just-measured wire
            for q in re.findall(r"q\[(\d+)\]", it[1]):
                if int(q) != reset_ok[0]: state[int(q)] = "other"
            reset_ok[0] = None
    flush(None)
    assert not pending_z, pending_z
    lines = list(head)
    for it in out:
        lines.append(it[1] if it[0] == "raw" else f"{it[1]} " + ", ".join(f"q[{q}]" for q in it[2]) + ";")
    open(path_out, "w").write("\n".join(lines) + "\n")
    log(str(stats))
    return stats


if __name__ == "__main__":
    import verify_qasm
    a = sys.argv
    run2(a[1], a[2], a[3], int(a[4]) if len(a) > 4 else 20, int(a[5]) if len(a) > 5 else 30)
    print("before", counts(a[1]), "after", counts(a[3]))
    print("F", verify_qasm.check(a[3], a[2]))
