"""Routing for ``qoptct hw``: a SABRE-style SWAP router with three additions that Qiskit's default router lacks.

1. Commutation-aware dependencies.  Every gate gets, per wire, an axis label: ``z`` if it is block-diagonal in the
   Z basis of that wire (rz, cz, cp, rzz, the control of cx, ...), ``x`` / ``y`` likewise, otherwise generic.  Two
   gates whose labels agree on every shared wire commute, so a gate only depends on the gates of the previous
   *different-axis* group on each wire.  The front layer therefore contains every gate that may legally be executed
   next (all ZZ terms of a QAOA layer, all controlled phases onto one QFT qubit, ...), and the router may execute
   them in whatever order suits the device.
2. SWAP absorption.  A SWAP placed next to a gate on the same pair of qubits is later merged with it into one
   two-qubit block that costs at most 3 native gates.  The router tracks the open block (4x4 unitary) on every
   coupler and charges a SWAP its exact marginal cost (number of CZ of ``SWAP * block`` minus that of ``block``,
   from the KAK decomposition): 0 after a generic block (e.g. XX+YY+ZZ), 1 after a CX / CP / RZZ, 3 on a fresh pair.
3. Calibration-aware distances.  Distances are shortest paths weighted by -log(1 - error) of each coupler (in units
   of an average SWAP), and unusable couplers / qubits are removed.

The output is a list of operations on physical qubits; ``to_circuit`` writes them as a Qiskit circuit with explicit
SWAP gates (block consolidation and translation to native gates happen afterwards, see compiler.py).
"""
import math
import random

import numpy as np

_SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)
_Z = np.diag([1.0, -1.0]).astype(complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_I2 = np.eye(2, dtype=complex)
_DEC = None


def kak_cost(u):
    """Number of CZ (= CX = ECR) gates needed for the 4x4 unitary u (0..3)."""
    global _DEC
    if _DEC is None:
        from qiskit.circuit.library import CZGate
        from qiskit.synthesis import TwoQubitBasisDecomposer
        _DEC = TwoQubitBasisDecomposer(CZGate())
    return int(_DEC.num_basis_gates(np.ascontiguousarray(u)))


_H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
_V = np.array([[1, 0], [0, 1j]], dtype=complex) @ _H          # V Z V^dagger = Y
_BAS = {"x": _H, "y": _V}
_ODD = [np.array([[bool((i ^ j) & (1 << w)) for j in range(4)] for i in range(4)]) for w in (0, 1)]
_CH = {(lab, w): (np.kron(_I2, m) if w == 0 else np.kron(m, _I2)) for lab, m in _BAS.items() for w in (0, 1)}


def _axis1(m):
    if abs(m[0, 1]) < 1e-9 and abs(m[1, 0]) < 1e-9:
        return "z"
    for lab in ("x", "y"):
        t = _BAS[lab].conj().T @ m @ _BAS[lab]
        if abs(t[0, 1]) < 1e-9 and abs(t[1, 0]) < 1e-9:
            return lab
    return None


def _axes2(m):
    """Per wire (0 = first qubit = least significant bit): the Pauli axis in whose eigenbasis m is block diagonal."""
    out = []
    for w in (0, 1):
        lab = None
        if np.abs(m[_ODD[w]]).max() < 1e-9:
            lab = "z"
        else:
            for name in ("x", "y"):
                c = _CH[(name, w)]
                if np.abs((c.conj().T @ m @ c)[_ODD[w]]).max() < 1e-9:
                    lab = name
                    break
        out.append(lab)
    return out


class Node:
    __slots__ = ("i", "q", "op", "mat", "axes", "two")

    def __init__(self, i, q, op, mat, axes):
        self.i, self.q, self.op, self.mat, self.axes, self.two = i, q, op, mat, axes, len(q) == 2


def to_nodes(qc, commute=True):
    """Nodes of a circuit made of one- and two-qubit gates (matrix and commutation axes per wire)."""
    nodes = []
    cache = {}
    for ins in qc.data:
        op = ins.operation
        if op.name in ("barrier", "delay"):
            continue
        q = tuple(qc.find_bit(b).index for b in ins.qubits)
        if len(q) > 2 or ins.clbits:
            raise ValueError(f"router needs one- and two-qubit unitary gates, got {op.name}")
        key = (op.name, tuple(float(p) for p in op.params)) if op.name != "unitary" else None
        if key is not None and key in cache:
            mat, axes = cache[key]
        else:
            mat = np.asarray(op.to_matrix(), dtype=complex)
            axes = [_axis1(mat)] if len(q) == 1 else _axes2(mat)
            if key is not None:
                cache[key] = (mat, axes)
        nodes.append(Node(len(nodes), q, op, mat, list(axes) if commute else [None] * len(q)))
    return nodes


def build_dag(nodes, nq):
    """(npred, succs): gate j depends on gate i iff they share a wire on which they are in different consecutive
    commutation groups (see module docstring).  Transitive dependencies follow through the group members."""
    cur_axis = [None] * nq
    cur = [[] for _ in range(nq)]
    prev = [[] for _ in range(nq)]
    preds = [set() for _ in nodes]
    for nd in nodes:
        for w, a in zip(nd.q, nd.axes):
            if a is not None and a == cur_axis[w]:
                preds[nd.i].update(prev[w])
                cur[w].append(nd.i)
            else:
                preds[nd.i].update(cur[w])
                prev[w] = cur[w]
                cur[w] = [nd.i]
                cur_axis[w] = a
    succs = [[] for _ in nodes]
    for j, ps in enumerate(preds):
        for i in ps:
            succs[i].append(j)
    return [len(p) for p in preds], succs


def reverse_nodes(nodes):
    out = []
    for k, nd in enumerate(reversed(nodes)):
        out.append(Node(k, nd.q, None, nd.mat.conj().T, nd.axes))
    return out


class _Block:
    """Live two-qubit block on a coupler: product ``u`` of its member gates, native cost, commutation axes of ``u``
    on its two wires, member op records.  Every op emitted after the block that touches one of its wires and is
    not a member commutes with ``u`` -- so the block can still be moved forward and merged with a later gate or
    SWAP on the same coupler ("floating" block)."""
    __slots__ = ("lo", "hi", "u", "cost", "swc", "ax", "items", "clean", "stamp")

    def __init__(self, lo, hi):
        self.lo, self.hi, self.u, self.cost, self.swc = lo, hi, np.eye(4, dtype=complex), 0, None
        self.ax, self.items, self.clean, self.stamp = (None, None), [], True, 0


class Params:
    """Router parameters.  mu: weight of the SWAP marginal cost (0 = plain SABRE); w_ext: look-ahead weight;
    ext: look-ahead size; near: exponent of the 1/distance weighting of front gates; noise: use calibration-weighted
    distances; floating: let blocks float over commuting gates (commutation-aware SWAP absorption)."""

    def __init__(self, mu=1.0, w_ext=0.5, ext=20, near=0.0, noise=True, decay=0.001, seed=0, floating=True,
                 tau=None):
        self.mu, self.w_ext, self.ext, self.near, self.noise, self.decay, self.seed = mu, w_ext, ext, near, noise, decay, seed
        self.floating = floating
        self.tau = tau  # mirror rule: right after a gate, SWAP its pair if the SWAP score is below tau (None = off)

    def __repr__(self):
        return (f"mu={self.mu},w={self.w_ext},ext={self.ext},near={self.near},noise={int(self.noise)},"
                f"fl={int(self.floating)},tau={self.tau},s={self.seed}")


def route(nodes, dag, dev, l2p, par, want_ops=True):
    """Route ``nodes`` from the initial mapping ``l2p`` (logical -> physical).  Returns (ops, final l2p, ncz, swaps,
    wlog) where ops = [("g", node index, physical qubits) | ("s", p, q)] in an order in which the members of every
    merged block are adjacent, ncz = number of native two-qubit gates of all blocks (KAK counts) and wlog = the same
    sum weighted by -log(1 - error) of each coupler (minus the log-fidelity of the two-qubit gates)."""
    npred0, succs = dag
    npred = list(npred0)
    rng = random.Random(par.seed)
    D = dev.wdist if par.noise else dev.hop
    hop = dev.hop
    l2p = list(l2p)
    p2l = {p: l for l, p in enumerate(l2p)}
    nl = len(l2p)
    ops = []
    live = {}
    wcost = 0.0
    ncz = 0
    wtrue = 0.0
    med = dev.elog_med
    nbr = dev.nbr
    floating = par.floating
    stamp = 0

    def ew(e):
        return min(dev.elog[e] / med, 8.0) if par.noise else 1.0

    def emit(ws, mat, axes, rec):
        """Append one op (gate or SWAP) on physical wires ``ws``; maintain the live blocks."""
        nonlocal wcost, stamp, ncz, wtrue
        stamp += 1
        two = len(ws) == 2
        join = None
        if two:
            e = (ws[0], ws[1]) if ws[0] < ws[1] else (ws[1], ws[0])
            join = live.get(e)
        clash = []
        for w, a in zip(ws, axes):
            for v in nbr[w]:
                e2 = (w, v) if w < v else (v, w)
                b = live.get(e2)
                if b is None or b is join:
                    continue
                if floating and a is not None and a == (b.ax[0] if w == b.lo else b.ax[1]):
                    b.clean = False
                    continue
                clash.append(b)
        if clash:
            if not two:
                join = max(clash, key=lambda b: b.stamp)
            for b in clash:
                if b is not join:
                    live.pop((b.lo, b.hi), None)
        if join is None:
            if not two:
                if want_ops:
                    ops.append(rec)
                return
            join = _Block(e[0], e[1])
            live[e] = join
        elif want_ops and not join.clean:
            moved = []
            for r in join.items:
                r[0] = False
                nr = [True, r[1], r[2], r[3]]
                moved.append(nr)
                ops.append(nr)
            join.items = moved
        if want_ops:
            ops.append(rec)
            join.items.append(rec)
        join.clean = True
        join.stamp = stamp
        if two:
            m = mat if ws[0] < ws[1] else _SWAP @ mat @ _SWAP
            join.u = m @ join.u
            c = kak_cost(join.u)
            wcost += (c - join.cost) * ew((join.lo, join.hi))
            ncz += c - join.cost
            wtrue += (c - join.cost) * dev.elog[(join.lo, join.hi)]
            join.cost = c
        else:
            join.u = (np.kron(_I2, mat) if ws[0] == join.lo else np.kron(mat, _I2)) @ join.u
        join.swc = None
        join.ax = _axes2(join.u) if floating else (None, None)

    def swap_cost(p, q):
        e = (p, q) if p < q else (q, p)
        b = live.get(e)
        if b is None:
            return 3.0 * ew(e)
        if b.swc is None:
            b.swc = kak_cost(_SWAP @ b.u) - b.cost
        return b.swc * ew(e)

    front = set()                      # ready two-qubit nodes (adjacent or not)
    qfront = [set() for _ in range(nl)]
    exe = []                           # ready two-qubit nodes that were adjacent when last looked at
    ready = [i for i, c in enumerate(npred) if c == 0]
    remaining = len(nodes)
    ext_dirty = True
    qext = [[] for _ in range(nl)]
    n_ext = 0

    def drain():
        """Execute all ready one-qubit gates; move ready two-qubit gates to the front."""
        nonlocal remaining, ext_dirty
        while ready:
            i = ready.pop()
            nd = nodes[i]
            if nd.two:
                front.add(i); qfront[nd.q[0]].add(i); qfront[nd.q[1]].add(i)
                if dev.has_edge(l2p[nd.q[0]], l2p[nd.q[1]]):
                    exe.append(i)
                ext_dirty = True
                continue
            ws = (l2p[nd.q[0]],)
            emit(ws, nd.mat, nd.axes, [True, "g", i, ws])
            remaining -= 1
            for j in succs[i]:
                npred[j] -= 1
                if npred[j] == 0:
                    ready.append(j)

    def execute2(i):
        nonlocal remaining, ext_dirty
        nd = nodes[i]
        ws = (l2p[nd.q[0]], l2p[nd.q[1]])
        emit(ws, nd.mat, nd.axes, [True, "g", i, ws])
        front.discard(i); qfront[nd.q[0]].discard(i); qfront[nd.q[1]].discard(i)
        remaining -= 1
        ext_dirty = True
        for j in succs[i]:
            npred[j] -= 1
            if npred[j] == 0:
                ready.append(j)
        return ws

    def recheck(ls):
        for l in ls:
            for i in qfront[l]:
                nd = nodes[i]
                if dev.has_edge(l2p[nd.q[0]], l2p[nd.q[1]]):
                    exe.append(i)

    def build_ext():
        nonlocal n_ext, ext_dirty
        for l in range(nl):
            qext[l] = []
        n_ext = 0
        ext_dirty = False
        if par.ext <= 0:
            return
        seen = set()
        layer = list(front)
        dec = {}
        while layer and n_ext < par.ext:
            nxt = []
            for i in layer:
                for j in succs[i]:
                    dec[j] = dec.get(j, npred[j]) - 1
                    if dec[j] == 0 and j not in seen:
                        seen.add(j)
                        nxt.append(j)
                        nd = nodes[j]
                        if nd.two:
                            qext[nd.q[0]].append(j); qext[nd.q[1]].append(j)
                            n_ext += 1
                            if n_ext >= par.ext:
                                break
                if n_ext >= par.ext:
                    break
            layer = nxt

    def do_swap(p, q):
        emit((p, q), _SWAP, (None, None), [True, "s", p, q])
        a, b = p2l.get(p), p2l.get(q)
        if a is not None:
            l2p[a] = q
        if b is not None:
            l2p[b] = p
        p2l.pop(p, None); p2l.pop(q, None)
        if a is not None:
            p2l[q] = a
        if b is not None:
            p2l[p] = b
        recheck([x for x in (a, b) if x is not None])

    decay = {}
    stall = 0
    nswap = 0
    limit = 12 + int(2 * math.sqrt(nl))
    near = par.near
    wext = 0.0

    def score(p, q):
        a, b = p2l.get(p), p2l.get(q)
        dF = 0.0
        dE = 0.0
        for l, pn in ((a, q), (b, p)):
            if l is None:
                continue
            po = l2p[l]
            for i in qfront[l]:
                nd = nodes[i]
                o = nd.q[1] if nd.q[0] == l else nd.q[0]
                if o == a or o == b:
                    continue  # both ends move: distance unchanged
                pp = l2p[o]
                d0 = D[po, pp]
                dF += (D[pn, pp] - d0) / (d0 ** near if near else 1.0)
            if wext:
                for i in qext[l]:
                    nd = nodes[i]
                    o = nd.q[1] if nd.q[0] == l else nd.q[0]
                    if o == a or o == b:
                        continue
                    pp = l2p[o]
                    dE += D[pn, pp] - D[po, pp]
        sc = 3.0 * (1.0 + max(decay.get(p, 0.0), decay.get(q, 0.0))) * (dF + wext * dE)
        if par.mu:
            sc += par.mu * swap_cost(p, q)
        return sc

    def prep_ext():
        nonlocal wext
        if ext_dirty:
            build_ext()
        wext = par.w_ext * min(len(front), n_ext) / n_ext if n_ext else 0.0

    while remaining:
        drain()
        did = False
        while exe:
            i = exe.pop()
            if i not in front:
                continue
            nd = nodes[i]
            if not dev.has_edge(l2p[nd.q[0]], l2p[nd.q[1]]):
                continue
            p, q = execute2(i)
            did = True
            drain()
            if par.tau is not None and front:
                prep_ext()
                if score(p, q) < par.tau:
                    do_swap(p, q)
                    nswap += 1
        if did:
            stall = 0
            decay.clear()
            continue
        if not remaining:
            break
        if not front:
            raise RuntimeError("router: no executable gate and empty front (cyclic dependencies?)")
        if stall >= limit:  # release valve: bring the closest blocked gate together along a shortest path
            import rustworkx as rx
            i = min(front, key=lambda k: hop[l2p[nodes[k].q[0]], l2p[nodes[k].q[1]]])
            pa, pb = l2p[nodes[i].q[0]], l2p[nodes[i].q[1]]
            path = list(rx.dijkstra_shortest_paths(dev.graph, pa, target=pb, weight_fn=float)[pb])
            for k in range(len(path) - 2):
                do_swap(path[k], path[k + 1]); nswap += 1
            stall = 0
            decay.clear()
            continue
        prep_ext()
        cands = set()
        for i in front:
            for l in nodes[i].q:
                p = l2p[l]
                for q in nbr[p]:
                    cands.add((p, q) if p < q else (q, p))
        best, best_s = [], None
        for p, q in cands:
            sc = score(p, q)
            if best_s is None or sc < best_s - 1e-9:
                best, best_s = [(p, q)], sc
            elif sc < best_s + 1e-9:
                best.append((p, q))
        p, q = best[0] if len(best) == 1 else rng.choice(sorted(best))
        do_swap(p, q)
        nswap += 1
        stall += 1
        decay[p] = decay.get(p, 0.0) + par.decay
        decay[q] = decay.get(q, 0.0) + par.decay
    return [tuple(r[1:]) for r in ops if r[0]], l2p, ncz, nswap, wtrue


def to_circuit(nodes, ops, n_phys):
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(n_phys)
    for o in ops:
        if o[0] == "s":
            qc.swap(o[1], o[2])
        else:
            qc.append(nodes[o[1]].op, list(o[2]))
    return qc


# ---- initial layouts -------------------------------------------------------------------------------------------
def interaction(nodes, nl):
    w = {}
    for nd in nodes:
        if nd.two:
            a, b = nd.q
            k = (a, b) if a < b else (b, a)
            w[k] = w.get(k, 0) + 1
    adj = [dict() for _ in range(nl)]
    for (a, b), c in w.items():
        adj[a][b] = c; adj[b][a] = c
    return adj


def qubit_quality(dev):
    """Lower is better: -log fidelity of readout + mean of the couplers of each qubit."""
    out = []
    for q in range(dev.n):
        es = [dev.elog[(min(q, r), max(q, r))] for r in dev.nbr[q]]
        if not es or q in dev.bad_qubits:
            out.append(math.inf)
            continue
        out.append(-math.log(max(1e-9, 1 - dev.ro_err[q])) + 3 * sum(es) / len(es))
    return out


def greedy_layout(adj, dev, rng, jitter=0.3):
    """Place logical qubits one by one (most connected to the placed ones first) on free physical qubits next to
    the occupied region, minimizing the weighted distance to the placed partners plus a qubit-quality term."""
    nl = len(adj)
    qual = qubit_quality(dev)
    good = sorted((q for q in range(dev.n) if math.isfinite(qual[q])), key=lambda q: qual[q])
    start = rng.choice(good[: max(1, len(good) // 3)])
    order0 = max(range(nl), key=lambda l: (sum(adj[l].values()), rng.random()))
    l2p = {order0: start}
    used = {start}
    med = dev.elog_med
    while len(l2p) < nl:
        l = max((x for x in range(nl) if x not in l2p),
                key=lambda x: (sum(c for o, c in adj[x].items() if o in l2p), sum(adj[x].values()), rng.random()))
        frontier = {q for p in used for q in dev.nbr[p] if q not in used and math.isfinite(qual[q])}
        if not frontier:
            frontier = {q for q in good if q not in used}
        best, bs = None, None
        for q in frontier:
            s = sum(c * dev.wdist[q, l2p[o]] for o, c in adj[l].items() if o in l2p)
            s += 0.2 * qual[q] / med + jitter * rng.random()
            if not adj[l]:
                s += 0.0
            if bs is None or s < bs:
                best, bs = q, s
        l2p[l] = best
        used.add(best)
    return [l2p[l] for l in range(nl)]


def find_path(dev, n, rng, tries=60):
    """A simple path of n usable physical qubits with a low total error (randomized depth-first search)."""
    qual = qubit_quality(dev)
    good = [q for q in range(dev.n) if math.isfinite(qual[q])]
    good.sort(key=lambda q: qual[q])
    best, bs = None, None
    for t in range(tries):
        start = rng.choice(good[: max(1, len(good) // 2)])
        path, used = [start], {start}
        stack = [iter(sorted((q for q in dev.nbr[start] if math.isfinite(qual[q])), key=lambda q: qual[q] * (0.5 + rng.random())))]
        steps = 0
        while stack and len(path) < n and steps < 4000:
            steps += 1
            nxt = next(stack[-1], None)
            if nxt is None:
                stack.pop(); used.discard(path.pop())
                continue
            if nxt in used:
                continue
            path.append(nxt); used.add(nxt)
            stack.append(iter(sorted((q for q in dev.nbr[nxt] if q not in used and math.isfinite(qual[q])),
                                     key=lambda q: qual[q] * (0.5 + rng.random()))))
        if len(path) == n:
            sc = sum(dev.elog[(min(a, b), max(a, b))] for a, b in zip(path, path[1:])) * 3 + \
                sum(-math.log(max(1e-9, 1 - dev.ro_err[q])) for q in path)
            if bs is None or sc < bs:
                best, bs = list(path), sc
    return best


def linear_order(adj, rng=None):
    """Order of the logical qubits for a line: grow a chain at both ends by the strongest interaction."""
    nl = len(adj)
    rest = set(range(nl))
    a = max(rest, key=lambda l: (sum(adj[l].values()), -l))
    chain = [a]
    rest.discard(a)
    while rest:
        best = None
        for end, side in ((chain[0], 0), (chain[-1], 1)):
            for l in rest:
                w = adj[end].get(l, 0) * 4 + sum(c for o, c in adj[l].items() if o not in rest) + (rng.random() * 0.5 if rng else 0)
                if best is None or w > best[0]:
                    best = (w, l, side)
        _, l, side = best
        chain.insert(0, l) if side == 0 else chain.append(l)
        rest.discard(l)
    return chain


def line_layout(order, path):
    l2p = [None] * len(order)
    for k, l in enumerate(order):
        l2p[l] = path[k]
    return l2p
