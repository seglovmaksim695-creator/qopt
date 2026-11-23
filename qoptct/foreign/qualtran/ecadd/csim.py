"""Concrete (not symbolic) executor of the same translated program: bits are Python integers.
Boxes take their values from a dictionary (a counterexample of the prover) or from a Qualtran oracle.
Returns the list of violations."""
from .core import lsb


def word(vals, ws):
    return sum(vals[q] << i for i, q in enumerate(ws))


def setw(vals, ws, v):
    for i, q in enumerate(ws):
        vals[q] = (v >> i) & 1


def simulate(prog, inputs, n, p, opaque_values=None, orc=None, post=None):
    """inputs: name -> int. opaque_values: list of dicts for the forward boxes (in order) or None (then orc)."""
    N = prog['N']
    vals = [0] * N
    for name, ws in prog['inputs'].items():
        setw(vals, ws, inputs[name])
    viol, recs, fwd = [], [], 0
    for idx, op in enumerate(prog['ops']):
        k = op[0]
        if k == 'x':
            vals[op[1]] ^= 1
        elif k == 'cx':
            vals[op[2]] ^= vals[op[1]]
        elif k == 'and':
            a, b, t = op[1:]
            if vals[t]:
                viol.append(('and', idx))
            vals[t] = vals[a] & vals[b]
        elif k == 'erase':
            t, ((cz, a, b),) = op[1], op[2]
            if vals[t] != (vals[a] & vals[b]):
                viol.append(('erase', idx))
            vals[t] = 0
        elif k == 'opaque':
            leaf = prog['leaves'][op[1]]
            kind = leaf['kind']
            rin = {r: lsb(ws) for r, ws in leaf['regs_in'].items()}
            rout = {r: lsb(ws) for r, ws in leaf['regs_out'].items()}
            cur = {r: word(vals, ws) for r, ws in rin.items()}
            if kind in ('MM', 'KINV'):
                for r in (('x', 'y') if kind == 'MM' else ('x',)):
                    if cur[r] >= p:
                        viol.append((kind + ':domain', idx))
                if opaque_values is not None:
                    out = opaque_values[fwd]
                else:
                    out = orc(kind, cur)
                regs = ('target', 'qrom_indices', 'reduced') if kind == 'MM' else ('x', 'junk')
                for r in regs:
                    setw(vals, rout[r], out[r])
                recs.append(dict(kind=kind, inp=cur, out={r: out[r] for r in regs}, used=False))
                fwd += 1
            else:
                fk = kind[:-4]
                regs = ('target', 'qrom_indices', 'reduced') if fk == 'MM' else ('x', 'junk')
                for rec in reversed(recs):
                    if rec['kind'] != fk or rec['used']:
                        continue
                    same = all(cur[r] == rec['out'][r] for r in regs)
                    if fk == 'MM':
                        same = same and cur['x'] == rec['inp']['x'] and cur['y'] == rec['inp']['y']
                    if same:
                        rec['used'] = True
                        break
                else:
                    viol.append((kind + ':pairing', idx))
                    rec = None
                if fk == 'MM':
                    for r in regs:
                        setw(vals, rin[r], 0)
                else:
                    setw(vals, rin['junk'], 0)
                    setw(vals, rin['x'], rec['inp']['x'] if rec else 0)
    outs = {name: word(vals, ws) for name, ws in prog['outputs'].items()}
    outw = {q for ws in prog['outputs'].values() for q in ws}
    if any(vals[q] for q in range(N) if q not in outw):
        viol.append(('nonzero_wire', -1))
    if post is not None:
        for k2, v in post.items():
            if outs.get(k2) != v:
                viol.append(('output:' + k2, -1))
    return outs, viol, recs
