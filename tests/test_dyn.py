# SPDX-License-Identifier: Apache-2.0
"""qoptct.dyn: measurement-based erasure vs reverse computation on noisy hardware."""
import math

import pytest

from qoptct.dyn import gen, lower, hw, choose

aer = pytest.importorskip('qiskit_aer')


def _toy_model(nq=6, ro=0.01, meas=2e-6, ff=6e-7):
    g1 = {}
    for q in range(nq):
        g1[('sx', q)] = (1e-4, 3e-8); g1[('x', q)] = (1e-4, 3e-8); g1[('rz', q)] = (0.0, 0.0)
    g2 = {(a, a + 1): (2e-3, 7e-8) for a in range(nq - 1)}
    return hw.HwModel('toy', nq, [2e-4] * nq, [1.5e-4] * nq, [ro] * nq, [meas] * nq, g1, g2, ff_latency=ff)


@pytest.mark.parametrize('fam,size', [('adder', 2), ('adder', 5), ('cmp', 1), ('cmp', 4), ('modadd', 2),
                                      ('modadd', 3), ('lookup', 3), ('lookadd', 2), ('lookadd', 3), ('mult', 2)])
def test_blocks_classical_all_inputs(fam, size):
    for B in gen.build(fam, size):
        ok, why, n, exhaustive = B.check(1 << 16)
        assert ok, (B.name, B.variant, why)
        assert exhaustive


def test_classical_check_catches_wrong_fixup():
    B = gen.comparator(2)
    i = next(k for k, o in enumerate(B.C.ops) if o[0] == 'MZ')
    t, corr = B.C.ops[i][1], B.C.ops[i][2]
    B.C.ops[i] = ('MZ', t, (('z', corr[0][1]),))
    ok, why, _, _ = B.check()
    assert not ok and 'correction' in why


@pytest.mark.parametrize('fam,size', [('adder', 3), ('cmp', 2), ('modadd', 2), ('lookup', 3), ('lookadd', 2)])
def test_lowerings_are_exact(fam, size):
    """every lowering equals the ideal permutation on random product states, for random measurement outcomes."""
    for label, B, erase, vt in choose.variants(fam, size):
        assert lower.check_quantum(B, erase, vt, trials=2, shots=6) > 1 - 1e-9, (B.name, label)


def test_quantum_check_catches_missing_fixup():
    from qiskit import QuantumCircuit
    B = gen.adder(2)
    good = lower.to_qiskit(B.C, 'meas', via_target=False)
    bad = QuantumCircuit(*good.qregs, *good.cregs)
    for inst in good.data:
        if inst.operation.name == 'if_else':
            with bad.if_test((inst.operation.condition[0], 1)):
                bad.x(inst.qubits[[q for q in inst.operation.params[0].qubits].index(
                    next(bi.qubits[0] for bi in inst.operation.params[0].data if bi.operation.name == 'x'))])
        else:
            bad.append(inst)
    assert lower.check_quantum(B, body=bad, trials=3, shots=16) < 0.99


def test_schedule_counts_measurement_and_feedforward_idle():
    from qiskit import QuantumCircuit
    M = _toy_model()
    qc = QuantumCircuit(3, 1)
    qc.cz(0, 1); qc.cz(1, 2); qc.measure(1, 0)
    with qc.if_test((qc.clbits[0], 1)):
        qc.cz(0, 1)
    qc.sx(2)
    S = hw.schedule(qc, M, live=[0, 2], final=[0, 2])
    t_cz = 7e-8
    assert S.n2q == 2 and S.n2q_cond == 1 and S.n_meas == 1
    # cz, cz, measure, feed-forward, conditional cz
    assert math.isclose(S.duration, 2 * t_cz + 2e-6 + 6e-7 + t_cz, rel_tol=1e-9)
    idle0 = sum(it[2] for it in S.items if it[0] == 'idle' and it[1] == 0)
    assert math.isclose(idle0, t_cz + 2e-6 + 6e-7, rel_tol=1e-9)        # qubit 0 waits for the outcome
    # expected fidelity: gates, readout, idles; the conditional gate enters with weight 1/2
    f_cond = 0.5 * math.log(1 - 2e-3)
    assert S.logF < 2 * math.log(1 - 2e-3) + math.log(1 - 0.01) + f_cond
    # a static circuit of the same gates is shorter and has no readout term
    qs = QuantumCircuit(3); qs.cz(0, 1); qs.cz(1, 2); qs.cz(0, 1); qs.sx(2)
    S2 = hw.schedule(qs, M, live=[0, 2], final=[0, 2])
    assert S2.duration < 1e-6 and S2.fid > S.fid and S2.parts['ro'] == 0


def test_site_rule_flips_with_fast_accurate_measurement():
    M = _toy_model()
    assert choose.site_rule(M, 10)['better'] == 'unitary'
    M2 = _toy_model(ro=1e-4, meas=2e-8, ff=1e-8)
    assert choose.site_rule(M2, 10)['better'] == 'meas'


def test_pass_converts_measurement_erasure_and_keeps_semantics():
    from qiskit.transpiler import PassManager
    B = gen.comparator(2)
    qc = lower.to_qiskit(B.C, 'meas', via_target=False)
    assert qc.count_ops().get('if_else', 0) == 2
    p = choose.DynErasePass(mode='unitary')
    out = PassManager([p]).run(qc)
    assert p.last['sites'] == 2 and p.last['chosen'] == 'unitary'
    ops = out.count_ops()
    assert 'if_else' not in ops and 'measure' not in ops
    assert lower.check_quantum(B, body=out, trials=2, shots=4) > 1 - 1e-9
    keep = PassManager([choose.DynErasePass(mode='keep')]).run(qc)
    assert keep.count_ops().get('if_else', 0) == 2


def test_pass_reads_qoptct_meas_dialect():
    """the form written by qoptct.crypto.andcirc (h; measure; if (c) {cz; x}) is recognised."""
    from qiskit import qasm3
    pytest.importorskip('qiskit_qasm3_import')
    from qoptct.crypto.andcirc import Circ
    C = Circ(); a = C.new(2); b = C.new(2); C.add(b, a)
    qc = qasm3.loads(C.qasm(first=a + b))
    out, n = choose.to_unitary_erasure(qc)
    assert n == 1 and 'if_else' not in out.count_ops()


def test_choose_on_fake_backend_and_ideal_measurement():
    pytest.importorskip('qiskit_ibm_runtime')
    M = hw.HwModel.from_backend(choose.resolve_backend('FakeBoston'))
    r = choose.choose('adder', 3, M, seeds=(0, 1), simulate=True, shots=60, n_inputs=4)
    labels = {row['variant'] for row in r['rows']}
    assert labels == {'U', 'M', 'Md', 'T'}
    assert r['best'] in ('U', 'T')                       # today's calibration: the static circuit wins
    for row in r['rows']:
        assert 0 < row['f_est'] < 1 and abs(row['pZ'] - row['f_est']) < 0.15
    # a hypothetical device with an instant, perfect measurement chain: the dynamic variant is at least as good
    M2 = M.scaled(ro=0.0, meas_dur=1e-3, ff_latency=1e-10)
    sc = choose.rescore(r, M2)
    assert max(sc['M'], sc['Md']) > sc['U']


def test_cli_choose_smoke(capsys):
    pytest.importorskip('qiskit_ibm_runtime')
    from qoptct.dyn.__main__ import main
    main(['choose', 'cmp', '1', '--backend', 'FakeBoston', '--seeds', '1'])
    out = capsys.readouterr().out
    assert 'best by estimate:' in out and 'cmp1' in out
