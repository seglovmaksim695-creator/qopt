# SPDX-License-Identifier: Apache-2.0
"""Optional equivalence check of UNITARY circuits with MQT QCEC (pip install "qoptct[qcec]").

Both circuits are built directly as mqt.core QuantumComputation objects (no text round trip). The ZX checker runs
first, then the decision-diagram alternating checker; only an `equivalent*` answer is accepted as a proof, and
`not_equivalent` is reported as a refutation. QCEC's own timeout does not always interrupt a heavy DD run, so the
check runs in a child process with a hard limit (qoptct.verify._limits).
Validation in the lab (agent b3): 56 mutants on 8 circuits, 0 false `equivalent`; all `not_equivalent` were correct.
"""
import time

from ._parse import parse_qasm_text, parse_qc, is_unitary_program
from ._limits import run_limited

_1Q = {'h', 'x', 'y', 'z', 's', 'sdg', 't', 'tdg'}


def available():
    try:
        import mqt.qcec  # noqa: F401
        from mqt.core.ir import QuantumComputation  # noqa: F401
        return True
    except Exception:
        return False


def _build(n, gates):
    from mqt.core.ir import QuantumComputation
    qc = QuantumComputation(n)
    for g, q in gates:
        if g in _1Q:
            getattr(qc, g)(q[0])
        elif g in ('cx', 'ccx', 'mcx'):
            if len(q) == 2:
                qc.cx(q[0], q[1])
            else:
                qc.mcx(set(q[:-1]), q[-1])
        elif g in ('cz', 'ccz', 'mcz'):
            if len(q) == 2:
                qc.cz(q[0], q[1])
            else:
                qc.mcz(set(q[:-1]), q[-1])
        elif g == 'swap':
            qc.swap(q[0], q[1])
        elif g == 'id':
            pass
        else:
            raise ValueError('unsupported gate for QCEC: ' + g)
    return qc


def _verify(qasm_text, qc_path, timeout):
    from mqt import qcec
    N, nb, ops = parse_qasm_text(qasm_text)
    names, ins, ref = parse_qc(qc_path)
    n = len(names)
    if not is_unitary_program(ops) or N != n:
        return dict(status='UNKNOWN', detail='QCEC is used only for unitary circuits on the reference wires')
    a = _build(n, ref)
    b = _build(n, [(op[1], op[2]) for op in ops])
    t0 = time.time()
    last = None
    for name, kw in (('zx', dict(run_zx_checker=True, run_alternating_checker=False)),
                     ('dd_alternating', dict(run_zx_checker=False, run_alternating_checker=True))):
        left = max(1.0, timeout - (time.time() - t0))
        r = qcec.verify(a, b, timeout=left / 2 if name == 'zx' else left, run_simulation_checker=False,
                        run_construction_checker=False, parallel=False, **kw)
        e = str(r.equivalence).split('.')[-1]
        last = e
        if e.startswith('equivalent'):
            return dict(status='PROVED', detail='QCEC %s checker: %s' % (name, e))
        if e == 'not_equivalent':
            return dict(status='REFUTED', detail='QCEC %s checker: not_equivalent' % name)
    return dict(status='UNKNOWN', detail='QCEC: %s' % last)


def check_qcec(qasm_text, qc_path, budget_sec=60.0):
    """-> dict(status=PROVED|REFUTED|UNKNOWN, detail)."""
    if not available():
        return dict(status='UNKNOWN', detail='mqt.qcec is not installed (pip install "qoptct[qcec]")')
    r = run_limited(_verify, (qasm_text, qc_path, max(1.0, budget_sec - 1.0)), budget_sec=budget_sec + 2)
    if r[0] == 'ok':
        return r[1]
    if r[0] == 'timeout':
        return dict(status='UNKNOWN', detail='QCEC timed out after %.0f s' % budget_sec)
    return dict(status='UNKNOWN', detail='QCEC error: ' + r[1])
