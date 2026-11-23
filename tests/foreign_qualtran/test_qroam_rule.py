"""Rule C1 (QROAMCleanAdjoint):
1) basis simulation with bridge.sim (all selections x all target values), phases of every outcome;
2) dense simulation (numpy) on a superposition entangled with an external register, all measurement branches;
3) canaries: wrong_table / wrong_sel / extra_phase must NOT pass.
"""
import itertools
import json
import sys
import time

import numpy as np
import pytest

pytest.importorskip("qualtran")
pytest.importorskip("cirq")

from qoptct.foreign.qualtran.bridge.qroam_clean_adjoint import qroam_clean_adjoint_atom_rule, chained_atom_rule
from qoptct.foreign.qualtran.bridge.qualtran_layout import lower_bloqgraph
from qoptct.foreign.qualtran.bridge.rules import translate
from qoptct.foreign.qualtran.bridge import sim
from qualtran.bloqs.data_loading.qroam_clean import QROAMClean


def build(data, log_block, num_controls=0, tb=None):
    q = QROAMClean.build_from_data(np.array(data, dtype=object) if not isinstance(data, np.ndarray) else data,
                                   log_block_sizes=log_block, num_controls=num_controls)
    return q, q.adjoint().qroam_clean_adjoint_bloq


def lower(adj, canary=None):
    rule = lambda b, i, a, p: qroam_clean_adjoint_atom_rule(b, i, a, p, _canary=canary)
    low = lower_bloqgraph(adj, extra_atom_rule=rule)
    tr = translate(low, low['inputs'], low['outputs'])
    return low, tr


def words(low):
    """name -> LSB-first wire list for every input element"""
    return low['inputs']


def set_bits(val, wires):
    return sum(((val >> k) & 1) << w for k, w in enumerate(wires))


def case_inputs(q, adj, low):
    """iterator (sel_tuple, ctrl, target_values dict elem-> int)"""
    pass


def basis_test(adj, q, canary=None, name=''):
    low, tr = lower(adj, canary)
    nq, ops = tr['nq'], tr['ops']
    ins = low['inputs']
    sel_names = [r.name for r in adj.selection_registers]
    ctl_names = [r.name for r in adj.control_registers]
    tgt = [(r, e) for r in adj.target_registers for e in (np.ndindex(*r.shape) if r.shape else [()])]
    from qoptct.foreign.qualtran.bridge.qualtran_layout import element_name
    datas = adj.data_or_shape
    data_shape = tuple(adj.data_shape)
    b = {r.name: r.dtype.num_qubits for r in adj.target_registers}
    n_ok = n_viol = 0
    bad = []
    sel_ranges = [range(1 << r.dtype.num_qubits) for r in adj.selection_registers]
    tgt_bits = sum(b[r.name] for r, e in tgt)
    ctrl_vals = list(itertools.product((0, 1), repeat=len(ctl_names)))
    total = 0
    import random
    rnd = random.Random(1)
    for sel in itertools.product(*sel_ranges):
        for cv in ctrl_vals:
            active = all(cv)
            inrange = all(s < d for s, d in zip(sel, data_shape))
            def expect(r, e):
                di = list(adj.target_registers).index(r)
                return int(datas[di][tuple(sel) + e]) if (active and inrange) else 0
            good = {(r.name, e): expect(r, e) for r, e in tgt}
            # входы: правильный, и набор неправильных (все, если мало)
            cands = [good]
            if tgt_bits <= 7:
                for vals in itertools.product(*[range(1 << b[r.name]) for r, e in tgt]):
                    cands.append({(r.name, e): v for (r, e), v in zip(tgt, vals)})
            else:
                for _ in range(24):
                    d = dict(good); k = rnd.choice(list(d)); d[k] ^= 1 << rnd.randrange(b[k[0]]); cands.append(d)
            for c in cands:
                init = 0
                for n_, v in zip(sel_names, sel):
                    init |= set_bits(v, ins[n_])
                for n_, v in zip(ctl_names, cv):
                    init |= set_bits(v, ins[n_])
                for (rn, e), v in c.items():
                    init |= set_bits(v, ins[element_name(rn, e)])
                val, ph = sim.run(nq, ops, init)
                is_good = (c == good)
                keep = init
                for (rn, e) in c:
                    for w in ins[element_name(rn, e)]:
                        keep &= ~(1 << w)
                total += 1
                if is_good:
                    if val == keep and not any(ph):
                        n_ok += 1
                    else:
                        bad.append(('good_not_clean', sel, cv))
                else:
                    if any(ph):
                        n_viol += 1
                    else:
                        bad.append(('violation_not_detected', sel, cv, str(c)))
    return dict(name=name, canary=canary, total=total, good_ok=n_ok, violations_detected=n_viol, bad=bad[:5], n_bad=len(bad),
                nq=nq, ops=len(ops))


def dense_test(adj, q, canary=None, name=''):
    """dense simulation: sum_s amp_s |s>|s>_ref|target=data[s]>, all measurement outcomes."""
    low, tr = lower(adj, canary)
    nq, ops = tr['nq'], tr['ops']
    ins = low['inputs']
    from qoptct.foreign.qualtran.bridge.qualtran_layout import element_name
    sel_names = [r.name for r in adj.selection_registers]
    ctl_names = [r.name for r in adj.control_registers]
    assert len(sel_names) == 1
    sb = adj.selection_registers[0].dtype.num_qubits
    datas = adj.data_or_shape
    data_len = adj.data_shape[0]
    tgt = [(r, e) for r in adj.target_registers for e in (np.ndindex(*r.shape) if r.shape else [()])]
    nref = sb + len(ctl_names)
    N = nq + nref
    rng = np.random.default_rng(7)
    init = np.zeros(1 << N, dtype=complex)
    ideal = np.zeros(1 << N, dtype=complex)
    for s in range(1 << sb):
        for cv in itertools.product((0, 1), repeat=len(ctl_names)):
            amp = rng.normal() + 1j * rng.normal()
            idx = set_bits(s, ins[sel_names[0]])
            for n_, v in zip(ctl_names, cv):
                idx |= set_bits(v, ins[n_])
            ref = s | (sum(v << k for k, v in enumerate(cv)) << sb)
            full_ref = ref << nq
            idx2 = idx
            for r, e in tgt:
                di = list(adj.target_registers).index(r)
                v = int(datas[di][(s,) + e]) if (all(cv) and s < data_len) else 0
                idx2 |= set_bits(v, ins[element_name(r.name, e)])
            init[idx2 | full_ref] += amp
            ideal[idx | full_ref] += amp
    # проверим порядок бит: индекс = sum bit<<wire
    E = sum(1 for o in ops if o[0] == 'erase')
    branches = [(init, 1.0)]
    for op in ops:
        new = []
        for psi, w in branches:
            if op[0] == 'x':
                t = op[1]
                idx = np.arange(1 << N)
                new.append((psi[idx ^ (1 << t)], w))
            elif op[0] == 'cx':
                a, t = op[1:]
                idx = np.arange(1 << N)
                m = (idx >> a) & 1
                new.append((psi[idx ^ (m << t)], w))
            elif op[0] == 'and':
                a, b_, t = op[1:]
                idx = np.arange(1 << N)
                assert not np.any(psi[((idx >> t) & 1) == 1]), 'and target nonzero'
                m = ((idx >> a) & 1) & ((idx >> b_) & 1)
                new.append((psi[idx ^ (m << t)], w))
            elif op[0] == 'erase':
                t, corr = op[1:]
                idx = np.arange(1 << N)
                # H на t
                lo = idx & ~(1 << t)
                hi = lo | (1 << t)
                bit = (idx >> t) & 1
                Hpsi = (psi[lo] + np.where(bit == 0, 1, -1) * psi[hi]) / np.sqrt(2)
                for m in (0, 1):
                    ph = Hpsi.copy()
                    ph[bit != m] = 0
                    if m:
                        sign = np.ones(1 << N)
                        for c in corr:
                            if c[0] == 'z':
                                sign *= 1 - 2 * ((idx >> c[1]) & 1)
                            else:
                                sign *= 1 - 2 * (((idx >> c[1]) & 1) & ((idx >> c[2]) & 1))
                        ph = ph * sign
                        ph = ph[idx ^ (1 << t)]   # X на t
                    new.append((ph, w))
            else:
                raise ValueError(op)
        branches = new
    worst = 0.0
    scale = 2 ** (-E / 2)
    for psi, w in branches:
        worst = max(worst, float(np.max(np.abs(psi - scale * ideal))))
    return dict(name=name, canary=canary, branches=len(branches), max_dev=worst, nq=nq, E=E, ok=worst < 1e-9)


CASES = [
    ('d8_K1', [3, 5, 6, 1, 7, 2, 0, 4], (0,), 0),
    ('d8_K2', [3, 5, 6, 1, 7, 2, 0, 4], (1,), 0),
    ('d5_K1', [3, 1, 2, 3, 1], (0,), 0),
    ('d8_K4', [1, 2, 3, 0, 3, 2, 1, 1], (2,), 0),
    ('ctrl_d4_K1', [2, 3, 1, 0], (0,), 1),
    ('ctrl_d8_K2', [3, 5, 6, 1, 7, 2, 0, 4], (1,), 1),
]
DENSE = [('d4_K1', [2, 3, 1, 0], (0,), 0), ('d4_K2', [2, 3, 1, 0], (1,), 0), ('ctrl_d4_K1', [2, 3, 1, 0], (0,), 1)]


CANARY_CASES = (CASES[0], CASES[1], CASES[4])


@pytest.mark.parametrize('name,data,lb,nc', CASES, ids=[c[0] for c in CASES])
def test_rule_basis(name, data, lb, nc):
    q, adj = build(data, lb, nc)
    r = basis_test(adj, q, None, name)
    assert r['n_bad'] == 0 and r['good_ok'] > 0, r


@pytest.mark.parametrize('name,data,lb,nc', DENSE, ids=[c[0] for c in DENSE])
def test_rule_dense(name, data, lb, nc):
    q, adj = build(data, lb, nc)
    r = dense_test(adj, q, None, name)
    assert r['ok'], r


@pytest.mark.parametrize('canary', ['wrong_table', 'wrong_sel', 'extra_phase'])
@pytest.mark.parametrize('name,data,lb,nc', CANARY_CASES, ids=[c[0] for c in CANARY_CASES])
def test_spoiled_rule_is_caught_basis(canary, name, data, lb, nc):
    q, adj = build(data, lb, nc)
    r = basis_test(adj, q, canary, name)
    assert r['n_bad'] > 0, ('spoiled rule NOT caught', canary, name, r)


@pytest.mark.parametrize('canary', ['wrong_table', 'wrong_sel', 'extra_phase'])
def test_spoiled_rule_is_caught_dense(canary):
    name, data, lb, nc = DENSE[0]
    q, adj = build(data, lb, nc)
    r = dense_test(adj, q, canary, name)
    assert not r['ok'], ('spoiled rule NOT caught', canary, r)
