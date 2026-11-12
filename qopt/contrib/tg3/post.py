"""tg3/post.py: program-level post-pass for deterministic circuits.

Rule R-CZ (phase kickback onto an X-measured wire):  a (conditional) CZ(g, d), where g is a wire that is later measured in the X
basis and only diagonal gates touch g (and d) in between, is equal to: measure g (bit m_g), then Z on d conditioned on cond AND m_g.
CZ costs 1 CX, the conditional Z costs 0.  Program ops: ("prep",q) | ("g",name,qubits,cond) | ("meas",q,k).
cond: 0 or a tuple of masks (each mask = XOR of classical bits, masks are AND-ed).
"""
import os as _qos, pathlib as _qpl
_QPKG = str(_qpl.Path(__file__).resolve().parents[2])
_QCACHE = _qpl.Path(_qos.environ.get('QOPT_CACHE', _qpl.Path.home() / '.cache' / 'qopt'))
import sys
sys.path.insert(0, (_QPKG + "/contrib/tg3"))
DIAG = {"cz", "s", "sdg", "z", "t", "tdg"}


def norm_cond(c):
    if not c: return ()
    return tuple(m for m in c if m)


def free_cz_pass(prog, n):
    return cancel_pairs(list(prog))   # R-CZ is INVALID for X-measurement (CZ kickback needs Z-basis), disabled


def _free_cz_pass_invalid(prog, n):
    prog = list(prog)
    measpos = {o[1]: i for i, o in enumerate(prog) if o[0] == "meas"}
    measbit = {o[1]: o[2] for o in prog if o[0] == "meas"}
    drop = set(); ins = {}      # ins[meas index] -> list of ops to insert after it
    for i, o in enumerate(prog):
        if o[0] != "g" or o[1] != "cz" or i in drop: continue
        a, b = o[2]
        cands = [w for w in (a, b) if w in measpos and measpos[w] > i]
        if not cands: continue
        g = min(cands, key=lambda w: measpos[w]); d = b if g == a else a
        mp = measpos[g]; ok = True
        for j in range(i + 1, mp):
            p = prog[j]
            if p[0] == "meas":
                if p[1] == d: ok = False; break
                continue
            if p[0] == "prep": continue
            if (g in p[2] or d in p[2]) and p[1] not in DIAG: ok = False; break
        if not ok: continue
        # d may not be measured before g's measurement (checked above); d must be alive after mp
        if d in measpos and measpos[d] < mp: continue
        drop.add(i)
        ins.setdefault(mp, []).append(("g", "z", (d,), norm_cond(o[3]) + (1 << measbit[g],)))
    out = []
    for i, o in enumerate(prog):
        if i in drop: continue
        out.append(o)
        for x in ins.get(i, []): out.append(x)
    return cancel_pairs(out)


def cancel_pairs(prog):
    """cancel identical conditional diagonal gates that are adjacent in the per-wire order (z z, cz cz) and merge s s -> z."""
    out = []
    for o in prog:
        if o[0] == "g" and o[1] in ("z", "cz") and o[3]:
            key = (o[1], tuple(sorted(o[2])), tuple(sorted(norm_cond(o[3]))))
            # find previous identical gate with no non-diagonal op on its wires in between
            k = len(out) - 1; found = None
            while k >= 0:
                p = out[k]
                if p[0] == "g" and p[1] == o[1] and tuple(sorted(p[2])) == key[1] and tuple(sorted(norm_cond(p[3]))) == key[2]:
                    found = k; break
                if p[0] == "meas" and p[1] in o[2]: break
                if p[0] == "g" and set(p[2]) & set(o[2]) and p[1] not in DIAG: break
                k -= 1
            if found is not None: del out[found]; continue
        out.append(o)
    return out


def counts(prog):
    T = sum(1 for o in prog if o[0] == "g" and o[1] in ("t", "tdg"))
    CX = sum(1 for o in prog if o[0] == "g" and o[1] in ("cx", "cz"))
    CXu = sum(1 for o in prog if o[0] == "g" and o[1] in ("cx", "cz") and not o[3])
    return dict(T=T, CX=CX, CX_uncond=CXu, measurements=sum(1 for o in prog if o[0] == "meas"))


def bits(m): return [k for k in range(m.bit_length()) if m >> k & 1]


def to_qasm3(prog, n, N, path=None):
    h = N - n
    L = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{N}] q;", f"bit[{max(h,1)}] c;"]
    ce = lambda m: " ^ ".join(f"c[{k}]" for k in bits(m))
    for o in prog:
        if o[0] == "prep": L.append(f"h q[{o[1]}];")
        elif o[0] == "meas": L += [f"h q[{o[1]}];", f"c[{o[2]}] = measure q[{o[1]}];"]
        else:
            s = f"{o[1]} " + ", ".join(f"q[{x}]" for x in o[2]) + ";"
            cond = o[3]
            if cond:
                s = "if (" + " && ".join(f"({ce(m)})" for m in norm_cond(cond)) + ") { " + s + " }"
            L.append(s)
    txt = "\n".join(L) + "\n"
    if path: open(path, "w").write(txt)
    return txt
