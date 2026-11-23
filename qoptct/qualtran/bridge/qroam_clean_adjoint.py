"""Rule C1: QROAMCleanAdjoint as a PRIMITIVE WITH A CONTRACT.

Contract (from the Qualtran documentation and Berry et al. 2019, appendix C):
  if the target block equals data[sel] (with the same data as the paired QROAMClean; with num_controls>0
  ctrl AND data[sel]), the primitive returns ALL target wires to 0, does not change sel and control and
  leaves no phase.

How the rule works (ONLY existing bridge leaves: x, cx, and, erase; the qoptct provers are unchanged):
  1. a "unitary" indicator u_s = [ctrl=1 AND sel=s] is built for every needed s (a tree of And);
  2. for EVERY bit of the target block: p := XOR of the indicators s with data[s] bit = 1 (a chain of CX/X);
     then ERASE(target_bit | p AND one): the X measurement of the bit and the correction Z(p) on outcome 1;
     then p is returned to 0 by the same CX/X;
  3. the indicators are erased in reverse order (ERASE with a CZ correction).
Physically this is exactly the reference measurement-based uncomputation: outcome m, phase (-1)^(m*data[sel]).
The precondition "target == data[sel]" is the correctness condition of every ERASE; it is proved by
domain_guard / kprove / the Walker (if not proved the result is NOT PROVED).

ASSUMPTIONS (to be listed in every result):
  A1. Trust in Qualtran: the real QROAMCleanAdjoint (which has no decomposition) implements this
      contract; we prove properties of OUR implementation of the contract, not of Qualtran gates.
  A2. Bit order: target element e, value v, bit k (least significant = 0) sits on the wire
      inputs[reg][e + (bitsize-1-k,)] (most significant first, like all Qualtran registers in the bridge).
  A3. The data come from bloq.data_or_shape (for the pair with QROAMClean this is batched_data_permuted).
      For a selection outside data_shape the rule takes data=0 (the QROM loads nothing there);
      if this is reachable with a nonzero target, the ERASE condition is not proved (UNKNOWN).
  A4. The data values are non-negative integers < 2^bitsize (checked; otherwise refusal).
"""
from __future__ import annotations

RULE_VERSION = 'qroam_clean_adjoint_rule_v1_20261123'
ASSUMPTIONS = (
    'A1 trust in Qualtran: QROAMCleanAdjoint implements the contract (it has no decomposition)',
    'A2 target bit order: most significant first on inputs[reg][e+(bitsize-1-k,)]',
    'A3 data = bloq.data_or_shape (batched_data_permuted of the paired QROAMClean); sel out of range -> data=0',
    'A4 data are non-negative integers < 2^bitsize',
)
MAX_INDICATORS = 4096


def _int(v, what):
    try:
        iv = int(v)
    except Exception as e:  # noqa
        raise ValueError(f'{what}: not an integer: {v!r}') from e
    if iv != v:
        raise ValueError(f'{what}: not an integer: {v!r}')
    return iv


def qroam_clean_adjoint_atom_rule(bloq, inputs, allocate, path, *, max_indicators=MAX_INDICATORS,
                                  _canary=None):
    """extra_atom_rule for bridge.qualtran_layout.lower_bloqgraph. Returns an AtomExpansion or None.

    _canary -- ONLY for negative controls (None in normal use):
      'wrong_table' (the table is inverted in one bit), 'wrong_sel' (the two lowest sel wires
      are swapped), 'extra_phase' (an extra Z correction: an ERASE with a shifted constant)."""
    import numpy as np
    from qualtran import Register, QAny, Side
    from qualtran.bloqs.data_loading.qroam_clean import QROAMCleanAdjoint
    from .qualtran_layout import AtomExpansion, LayoutError
    if type(bloq) is not QROAMCleanAdjoint:
        return None
    if not bloq.has_data():
        raise LayoutError(f'{path}: QROAMCleanAdjoint has no data (only a shape): the rule is impossible')
    ctrl_regs = list(bloq.control_registers)
    sel_regs = list(bloq.selection_registers)
    tgt_regs = list(bloq.target_registers)
    datas = bloq.data_or_shape
    if len(datas) != len(tgt_regs):
        raise LayoutError(f'{path}: number of tables != number of target registers')
    # literals: (wire, ...); ctrl is always cv=1
    lits = []          # wires in the order: ctrl..., sel_0 (LSB..MSB), sel_1 ...
    for r in ctrl_regs:
        for q in inputs[r.name].flat:
            lits.append(int(q))
    n_ctrl = len(lits)
    sel_slices = []
    for r in sel_regs:
        arr = inputs[r.name]
        bits = [int(q) for q in arr.reshape(-1)][::-1]     # LSB first
        if _canary == 'wrong_sel' and len(bits) >= 2 and r is sel_regs[0]:
            bits[0], bits[1] = bits[1], bits[0]
        sel_slices.append((len(lits), len(bits)))
        lits.extend(bits)
    data_shape = tuple(int(x) for x in bloq.data_shape)
    nsel = len(sel_regs)
    if len(data_shape) != nsel:
        raise LayoutError(f'{path}: data_shape does not match selection_registers')

    def assignment(s):
        a = [1] * n_ctrl
        for (start, width), val in zip(sel_slices, s):
            a.extend((val >> k) & 1 for k in range(width))
        return tuple(a)

    # bit tables: for every table i, element e, bit k -> the set of s with bit 1
    work = []   # (reg, e, k, [assignments])
    needed = set()
    first_bit_done = False
    for reg, d in zip(tgt_regs, datas):
        d = np.asarray(d, dtype=object)
        b = int(reg.dtype.num_qubits)
        tshape = tuple(reg.shape)
        if d.shape != data_shape + tshape:
            raise LayoutError(f'{path}: data shape {d.shape} != {data_shape + tshape}')
        for e in np.ndindex(*tshape) if tshape else [()]:
            for k in range(b):
                S = []
                for s in np.ndindex(*data_shape):
                    v = _int(d[s + e], f'{path}: data')
                    if v < 0 or v >> b:
                        raise LayoutError(f'{path}: data value {v} outside 0..2^{b}-1')
                    bit = (v >> k) & 1
                    if _canary == 'wrong_table' and not first_bit_done and s == tuple(0 for _ in data_shape):
                        bit ^= 1
                    if bit:
                        S.append(assignment(s))
                first_bit_done = True
                work.append((reg, e, k, b, S))
                needed.update(S)
    if len(needed) > max_indicators:
        raise LayoutError(f'{path}: too many indicators {len(needed)} > {max_indicators}')

    ops_and = []        # (nodeA, nodeB, wire)   node=(wire,cv)
    # --- first count how many new wires are needed, then allocate
    # tree: recursively over the range of literals
    plan_counter = [0]
    created = []        # (a_node, b_node, key) in creation order; the wire is assigned later

    def build(lo, hi, assigns):
        """-> {partial assignment: node}, node=('lit',wire,cv) or ('new',index)"""
        if hi - lo == 1:
            return {a: ('lit', lits[lo], a[0]) for a in assigns}
        mid = (lo + hi) // 2
        left = build(lo, mid, {a[:mid - lo] for a in assigns})
        right = build(mid, hi, {a[mid - lo:] for a in assigns})
        out = {}
        for a in sorted(assigns):
            idx = len(created)
            created.append((left[a[:mid - lo]], right[a[mid - lo:]]))
            out[a] = ('new', idx)
        return out

    top = build(0, len(lits), set(needed)) if needed and lits else {}
    if not lits and needed:
        raise LayoutError(f'{path}: neither selection nor control wires')
    n_new = len(created)
    scr = [Register('one', QAny(1), side=Side.RIGHT), Register('p', QAny(1), side=Side.RIGHT)]
    if n_new:
        scr.append(Register('ind', QAny(n_new), side=Side.RIGHT))
    got = allocate(scr, path + ':qroam_adj_scratch')
    one = int(got['one'].reshape(-1)[0])
    p = int(got['p'].reshape(-1)[0])
    ind = [int(q) for q in got['ind'].reshape(-1)] if n_new else []

    def node_ref(nd):
        return (ind[nd[1]], 1) if nd[0] == 'new' else (nd[1], nd[2])

    leaves = []

    def L(name, qs, **kw):
        leaves.append(dict(name=name, qubits=list(qs), path=path, **kw))

    L('x', [one])
    for idx, (na, nb) in enumerate(created):          # nodes are created so that children come before parents
        (wa, ca), (wb, cb) = node_ref(na), node_ref(nb)
        L('and', [wa, wb, ind[idx]], controls=[ca, cb])
    for (reg, e, k, b, S) in work:
        t = int(inputs[reg.name][tuple(e) + (b - 1 - k,)])
        contrib = []
        for a in S:
            w, cv = node_ref(top[a])
            contrib.append(('cx', [w, p]))
            if cv == 0:
                contrib.append(('x', [p]))
        for nm, qs in contrib:
            L(nm, qs)
        if _canary == 'extra_phase' and (reg is tgt_regs[0]) and e == ((0,) * len(reg.shape)) and k == 0:
            L('x', [p])           # an extra constant -> an extra phase / a violated condition
        L('erase', [p, one, t], controls=[1, 1])
        if _canary == 'extra_phase' and (reg is tgt_regs[0]) and e == ((0,) * len(reg.shape)) and k == 0:
            L('x', [p])
        for nm, qs in reversed(contrib):
            L(nm, qs)
    for idx in reversed(range(n_new)):
        na, nb = created[idx]
        (wa, ca), (wb, cb) = node_ref(na), node_ref(nb)
        L('erase', [wa, wb, ind[idx]], controls=[ca, cb])
    L('x', [one])
    consumed = [int(q) for r in tgt_regs for q in inputs[r.name].flat]
    retired = tuple(consumed + [one, p] + ind)
    outputs = {r.name: inputs[r.name] for r in bloq.signature.rights() if r.side == Side.THRU}
    return AtomExpansion(outputs, leaves, retired)


def chained_atom_rule(**kw):
    """Ready-made extra_atom_rule: A1 (controlled And) + QROAMCleanAdjoint."""
    from .qualtran_layout import controlled_and_atom_rule

    def rule(bloq, inputs, allocate, path):
        r = controlled_and_atom_rule(bloq, inputs, allocate, path)
        if r is not None:
            return r
        return qroam_clean_adjoint_atom_rule(bloq, inputs, allocate, path, **kw)
    return rule
