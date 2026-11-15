# SPDX-License-Identifier: Apache-2.0
"""Block layer of the b10 resynthesis pass: QASM -> blocks of CX/CZ/X/phase gates, their phase polynomial and linear map.
Block gates: cx, cz, x, z, s, sdg, t, tdg; everything else (h, measurements, if ...) is a raw line that closes blocks on its wires.
Phase weights are in units of pi/4: t=1, s=2, z=4, sdg=6, tdg=7.
"""
import re

from . import cliff
W = {'t': 1, 's': 2, 'z': 4, 'sdg': 6, 'tdg': 7}
_q = re.compile(r'q\[(\d+)\]')
_c = re.compile(r'c\[(\d+)\]')


def parse_qasm(text):
    """-> (header_lines, nq, items) ; item = ('g', name, wires) | ('raw', line, wires)"""
    head, items, nq = [], [], 0
    for ln in text.splitlines():
        s = ln.strip()
        if not s: continue
        if s.startswith('OPENQASM') or s.startswith('include'):
            head.append(ln); continue
        m = re.match(r'qubit\[(\d+)\]', s)
        if m:
            nq = int(m.group(1)); head.append(ln); continue
        m = re.match(r'(cx|cz|x|z|s|sdg|t|tdg)\s+(.*);$', s)
        if m:
            ws = [int(a) for a in _q.findall(m.group(2))]
            items.append(('g', m.group(1), ws)); continue
        ws = [int(a) for a in _q.findall(s)] + ['c%s' % a for a in _c.findall(s)]
        items.append(('raw', s, ws))
    return head, nq, items


def segment(items):
    """Split into blocks: list of ('block', sorted wires, gates) | ('raw', line, wires)."""
    out = []
    wire_blk = {}  # wire -> block dict
    def flush(b):
        for w in b['wires']: wire_blk.pop(w, None)
        out.append(('block', sorted(b['wires']), b['gates']))
    for it in items:
        if it[0] == 'g':
            _, name, ws = it
            bl = []
            for w in ws:
                b = wire_blk.get(w)
                if b is None:
                    b = {'wires': {w}, 'gates': []}; wire_blk[w] = b
                if b not in bl: bl.append(b)
            main = bl[0]
            for b in bl[1:]:
                main['wires'] |= b['wires']; main['gates'] += b['gates']
                for w in b['wires']: wire_blk[w] = main
            main['gates'].append((name, ws))
        else:
            seen = []
            for w in it[2]:
                b = wire_blk.get(w)
                if b is not None and all(b is not x for x in seen): seen.append(b)
            for b in seen: flush(b)
            out.append(('raw', it[1], it[2]))
    seen = []
    for b in list(wire_blk.values()):
        if all(b is not x for x in seen): seen.append(b)
    for b in seen: flush(b)
    return out


def block_poly(wires, gates):
    """-> dict(terms mask->weight mod8, rows (masks over local init vars), off (bits), ncx)"""
    idx = {w: i for i, w in enumerate(wires)}
    n = len(wires)
    rows = [1 << i for i in range(n)]
    off = [0] * n
    terms = {}
    ncx = 0
    for name, ws in gates:
        if name == 'cx':
            c, t = idx[ws[0]], idx[ws[1]]
            rows[t] ^= rows[c]; off[t] ^= off[c]; ncx += 1
        elif name == 'cz':
            a, b = idx[ws[0]], idx[ws[1]]
            ra, rb = rows[a], rows[b]
            for m, w in ((ra, 2), (rb, 2), (ra ^ rb, 6)):
                terms[m] = (terms.get(m, 0) + w) % 8
            if off[a]: terms[rb] = (terms.get(rb, 0) + 4) % 8
            if off[b]: terms[ra] = (terms.get(ra, 0) + 4) % 8
            ncx += 1
        elif name == 'x':
            off[idx[ws[0]]] ^= 1
        else:
            i = idx[ws[0]]
            w = W[name]
            if off[i]: w = -w
            m = rows[i]
            terms[m] = (terms.get(m, 0) + w) % 8
    terms = {m: w for m, w in terms.items() if w}
    return dict(n=n, terms=terms, rows=rows, off=off, ncx=ncx)


def poly_equal(p, q):
    return p['n'] == q['n'] and p['terms'] == q['terms'] and p['rows'] == q['rows'] and p['off'] == q['off']


_PG = {0: [], 1: ['t'], 2: ['s'], 3: ['s', 't'], 4: ['z'], 5: ['z', 't'], 6: ['sdg'], 7: ['tdg']}


def phase_gates(w, wire):
    """Weight mod 8 -> list of single-qubit phase gates (name, [wire])."""
    return [(g, [wire]) for g in _PG[w % 8]]


def emit_qasm(head, nq, items):
    L = list(head)
    for it in items:
        if it[0] == 'raw': L.append(it[1])
        else:
            for name, ws in it[2]:
                L.append(f"{name} " + ", ".join(f"q[{w}]" for w in ws) + ";")
    return "\n".join(L) + "\n"


def count(items_or_gates):
    T = CX = 0
    for it in items_or_gates:
        if it[0] == 'raw':
            s = it[1].split()[0]
            if s in ('cx', 'cz'): CX += 1
        else:
            for name, ws in it[2]:
                if name in ('t', 'tdg'): T += 1
                if name in ('cx', 'cz'): CX += 1
    return T, CX


def poly_equiv(p, q):
    """Functional equality of phase polynomials mod 8 (coefficients of monomials of degree <= 3) and of the affine map."""
    return p['n'] == q['n'] and p['rows'] == q['rows'] and p['off'] == q['off'] and cliff.coeffs(p['terms']) == cliff.coeffs(q['terms'])


def brute_equal(p, q):
    """Brute-force check over all x (n <= 14): phase values mod 8 agree on every input. None if n is too large."""
    n = p['n']
    if n > 14: return None
    tp = list(p['terms'].items()); tq = list(q['terms'].items())
    for x in range(1 << n):
        fa = sum(w for m, w in tp if bin(m & x).count('1') & 1)
        fb = sum(w for m, w in tq if bin(m & x).count('1') & 1)
        if (fa - fb) % 8: return False
    return True


def segment2(items):
    """Extended segmentation: a block is a connected component of wire segments (between H / raw items); components are only
    merged if the quotient graph stays acyclic, otherwise the CX/CZ becomes a raw barrier. Same output format as segment, topological order."""
    import heapq
    parent = []; succ = []; gates = []; first = []; kind = []; rawitem = []
    def new_node(k, it=None, idx=0):
        parent.append(len(parent)); succ.append(set()); gates.append([]); first.append(idx); kind.append(k); rawitem.append(it)
        return len(parent) - 1
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def reach(a, b):
        a = find(a); b = find(b)
        stack = [a]; seen = {a}; first_step = True
        while stack:
            x = stack.pop()
            for y in succ[x]:
                y = find(y)
                if y == b: return True
                if y not in seen: seen.add(y); stack.append(y)
        return False
    cur = {}; last = {}
    def edge(a, b):
        succ[find(a)].add(b)
    def union(a, b):
        a = find(a); b = find(b)
        if a == b: return a
        parent[b] = a
        succ[a] |= succ[b]; gates[a] += gates[b]; first[a] = min(first[a], first[b])
        return a
    def barrier(it, ws, idx):
        r = new_node('raw', it, idx)
        for w in ws:
            if cur.get(w) is not None: edge(cur[w], r); cur[w] = None
            elif last.get(w) is not None: edge(last[w], r)
        for w in ws: last[w] = r
    def seg_of(w, idx):
        if cur.get(w) is None:
            s = new_node('blk', None, idx)
            if last.get(w) is not None: edge(last[w], s)
            cur[w] = s
        return cur[w]
    for idx, it in enumerate(items):
        if it[0] == 'g':
            _, name, ws = it
            if name in ('cx', 'cz'):
                sa, sb = seg_of(ws[0], idx), seg_of(ws[1], idx)
                ca, cb = find(sa), find(sb)
                if ca != cb and (reach(ca, cb) or reach(cb, ca)):
                    barrier(('raw', f"{name} q[{ws[0]}], q[{ws[1]}];", ws), ws, idx)
                    continue
                c = union(sa, sb)
                for w in ws: cur[w] = c
                gates[c].append((idx, name, ws))
            else:
                s = seg_of(ws[0], idx)
                gates[find(s)].append((idx, name, ws))
        else:
            barrier(it, it[2], idx)
    roots = sorted({find(x) for x in range(len(parent))})
    roots = [r for r in roots if kind[r] == 'raw' or gates[r]]
    indeg = {r: 0 for r in roots}; adj = {r: set() for r in roots}
    for r in roots:
        for y in succ[r]:
            y = find(y)
            if y != r and y in adj and y not in adj[r]:
                adj[r].add(y); indeg[y] += 1
    hp = [(first[r], r) for r in roots if indeg[r] == 0]; heapq.heapify(hp)
    out = []; done = 0
    while hp:
        _, r = heapq.heappop(hp); done += 1
        if kind[r] == 'raw': out.append(rawitem[r])
        else:
            gl = sorted(gates[r]); ws = sorted({w for _, _, a in gl for w in a})
            out.append(('block', ws, [(nm, a) for _, nm, a in gl]))
        for y in adj[r]:
            indeg[y] -= 1
            if indeg[y] == 0: heapq.heappush(hp, (first[y], y))
    assert done == len(roots), 'cycle in quotient graph'
    return out
