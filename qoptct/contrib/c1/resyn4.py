"""c1/resyn4.py: resyn2 (exact blocks + F1 |+>-start freedom + F2 X-measured-end freedom) with PER-WIRE barriers:
non-unitary lines (measure / if-blocks) block only the wires they touch (layers2 from resyn3), so unitary blocks extend across
measurements of unrelated wires.  CLI: python resyn4.py in.qasm qc out.qasm [seeds variants rounds]"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys, re, random
sys.path.insert(0, (_QPKG + "/contrib/c1"))
from resyn2 import synth_gauge, parse, counts, block_cost
from resyn3 import layers2, raw_wires


def run4(path_in, qc, path_out, seeds=16, variants=24, base=0, log=print):
    import qoptct
    nd = len(qoptct.load_qc(qc)[0])
    head, items, nq = parse(path_in)
    ops = [("g", it[1], it[2]) if it[0] == "g" else ("raw", it[1], raw_wires(it[1])) for it in items]
    lay = layers2(ops)
    state = {q: ("zero" if q >= nd else "other") for q in range(nq)}
    out = []; stats = dict(blocks=0, improved=0, before=0, after=0, condz=0)
    pending = {}
    for li, (kind, x) in enumerate(lay):
        if kind == "B":
            gs = [(o[1], o[2]) for o in x]
            bw = sorted({q for _, qs in gs for q in qs})
            fstart = [q for q in bw if state[q] == "plus"]
            fend = {}
            if li + 1 < len(lay):
                grp = lay[li + 1][1]; touched = set()
                for j, o in enumerate(grp):
                    if o[0] == "raw":
                        m = re.match(r"c\[(\d+)\]\s*=\s*measure\s+q\[(\d+)\]", o[1])
                        if m:
                            a = int(m.group(2))
                            # wire a must have exactly one H before its measurement inside this group and nothing else
                            prev = [p for p in grp[:j] if a in p[2]]
                            if len(prev) == 1 and prev[0][0] == "g" and prev[0][1] == "h" and a in bw:
                                fend[a] = [c for c in bw if c != a and c not in touched]
                    touched |= set(o[2])
            stats["blocks"] += 1; stats["before"] += block_cost(gs)
            r = synth_gauge(gs, fstart, fend, seeds, variants, rng_base=base + stats["blocks"] * 1000)
            if r is not None:
                gs, extra = r; stats["improved"] += 1
                for a, cs in extra.items():
                    pending.setdefault(a, []).extend(cs); stats["condz"] += len(cs)
            stats["after"] += block_cost(gs)
            for _, qs in gs:
                for q in qs: state[q] = "other"
            out.extend(("g", g, qs) for g, qs in gs)
            continue
        reset_ok = None
        for o in x:
            if o[0] == "g":
                q = o[2][0]; state[q] = {"zero": "plus", "plus": "zero"}.get(state[q], "other")
                out.append(o); continue
            out.append(("raw", o[1]))
            m = re.match(r"c\[(\d+)\]\s*=\s*measure\s+q\[(\d+)\]", o[1])
            if m:
                k, a = int(m.group(1)), int(m.group(2))
                for c in pending.pop(a, []): out.append(("raw", f"if (c[{k}]) {{z q[{c}];}}"))
                state[a] = "measured"; reset_ok = (k, a)
            else:
                for q in o[2]:
                    if reset_ok and q == reset_ok[1] and re.search(rf"c\[{reset_ok[0]}\]\)", o[1]) and f"x q[{q}];" in o[1] \
                            and state[q] == "measured":
                        state[q] = "zero"
                    else:
                        state[q] = "other"
                reset_ok = None
        for q in range(nq):
            if state[q] == "measured": state[q] = "other"
    assert not pending, pending
    lines = list(head)
    for it in out:
        lines.append(it[1] if it[0] == "raw" else f"{it[1]} " + ", ".join(f"q[{q}]" for q in it[2]) + ";")
    open(path_out, "w").write("\n".join(lines) + "\n")
    log(str(stats))
    return stats


if __name__ == "__main__":
    import verify_qasm, shutil
    a = sys.argv
    seeds = int(a[4]) if len(a) > 4 else 16; var = int(a[5]) if len(a) > 5 else 24; rounds = int(a[6]) if len(a) > 6 else 2
    cur = a[1]; best = counts(a[1]); shutil.copy(a[1], a[3])
    for r in range(rounds):
        tmp = a[3] + f".r{r}"
        run4(cur, a[2], tmp, seeds, var, base=r * 100000)
        c = counts(tmp)
        f = verify_qasm.check(tmp, a[2])
        print("round", r, c, "F", f, flush=True)
        if f > 1 - 1e-9 and c < best:
            best = c; shutil.copy(tmp, a[3]); cur = tmp
        else:
            break
    print("before", counts(a[1]), "after", counts(a[3]))
    print("F", verify_qasm.check(a[3], a[2]))
