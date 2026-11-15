"""qoptct.verify.branchsim: renormalisation after measurements (0.5.0).

Before 0.5.0 a fixed branch was simulated without renormalisation; every X-basis measurement halves the norm, terms
below 1e-14 were dropped, and after ~50 measurements the branch "vanished" (confirm_branch then silently reported no
deviation, so a real counterexample on a long circuit could not be confirmed)."""
from qoptct.verify._parse import parse_qasm_text
from qoptct.verify.branchsim import run_branch, run_branch_scaled, confirm_branch

K = 120          # measurements in the chain


def _chain(drop=None):
    """|a> (x) |0>: K times {cx a->anc; h anc; measure anc -> c[k]; if (c[k]) z a; reset anc via x if (c[k])}.
    Each step is the identity on a up to the global factor 2^-1/2 (phase (-1)^(a m) fixed by the z correction)."""
    L = ['OPENQASM 3;', 'include "stdgates.inc";', 'qubit[2] q;', 'bit[%d] c;' % K]
    for k in range(K):
        L += ['cx q[0], q[1];', 'h q[1];', 'c[%d] = measure q[1];' % k]
        if k != drop:
            L.append('if (c[%d]) {z q[0];}' % k)
        L.append('if (c[%d]) {x q[1];}' % k)
    return '\n'.join(L) + '\n'


QC = '.v a\n.i a\n\nBEGIN\nEND\n'


def test_long_chain_keeps_the_branch():
    N, nb, ops = parse_qasm_text(_chain())
    assert nb == K
    for a in (0, 1):
        for m in ([0] * K, [1] * K, [k & 1 for k in range(K)]):
            st, lg = run_branch_scaled(ops, {a: 1.0 + 0j}, m)
            assert lg is not None and abs(lg + K / 2) < 1e-9          # |amplitude| = 2^(-K/2) exactly
            assert set(st) == {a} and abs(abs(st[a]) - 1) < 1e-9
            assert abs(st[a] - 1) < 1e-9                              # phase corrected
    out = run_branch(ops, {1: 1.0 + 0j}, [1] * K)                     # 2^-60: still representable, not dropped
    assert set(out) == {1} and abs(abs(out[1]) - 2.0 ** (-K / 2)) < 1e-25


def test_confirm_branch_long_chain(tmp_path):
    qc = tmp_path / 'id.qc'
    qc.write_text(QC)
    good = tmp_path / 'good.qasm'
    good.write_text(_chain())
    bad = tmp_path / 'bad.qasm'
    bad.write_text(_chain(drop=110))                                  # missing phase correction deep in the chain
    xs = [{0: 0}, {0: 1}]
    m1 = [1] * K
    assert confirm_branch(str(good), str(qc), xs, m1, 1) < 1e-9
    assert confirm_branch(str(bad), str(qc), xs, m1, 1) > 1e-6        # relative phase -1 between a=0 and a=1
    assert confirm_branch(str(bad), str(qc), xs, [0] * K, 1) < 1e-9   # outcome 0 at step 110: no correction needed
