"""Run an SMT-LIB script with Z3: the ``z3`` binary if there is one, else the ``z3-solver`` Python package.

Every check is "assert the negation of the obligation, expect ``unsat``".  The answers are returned as strings
(``sat`` / ``unsat`` / ``unknown``); a missing Z3 gives ``RuntimeError`` (a proof never silently passes without Z3).
"""
import os
import shutil
import subprocess
import sys
import time


def find_binary():
    """Path of the z3 binary: $QOPT_Z3, then PATH, then next to the interpreter; ``None`` if absent."""
    cand = os.environ.get('QOPT_Z3')
    if cand and os.path.isfile(cand):
        return cand
    found = shutil.which('z3')
    if found:
        return found
    near = os.path.join(os.path.dirname(sys.executable), 'z3')
    return near if os.path.isfile(near) else None


def version():
    """String like ``'Z3 version 5.3.3 (binary /opt/homebrew/bin/z3)'`` or ``'z3-solver 5.3.3'``; ``None`` if no Z3."""
    b = find_binary()
    if b:
        try:
            r = subprocess.run([b, '--version'], capture_output=True, text=True, timeout=30)
            return f'{r.stdout.strip()} (binary)'
        except Exception:  # noqa: BLE001
            return None
    try:
        import z3
        return f'z3-solver {z3.get_version_string()} (python)'
    except ImportError:
        return None


def solve(lines, checks, timeout_s=600):
    """``lines``: SMT-LIB declarations and assertions; ``checks``: [(tag, term)].  -> (answers, seconds, stderr, script)."""
    script = ['(set-option :produce-models true)', '(set-logic QF_LIA)'] + list(lines)
    for _tag, t in checks:
        script += ['(push 1)', f'(assert (not {t}))', '(check-sat)', '(pop 1)']
    txt = '\n'.join(script) + '\n'
    t0 = time.time()
    b = find_binary()
    if b:
        r = subprocess.run([b, '-in', f'-T:{timeout_s}'], input=txt, capture_output=True, text=True)
        ans = [x.strip() for x in r.stdout.splitlines() if x.strip() in ('sat', 'unsat', 'unknown')]
        return ans, round(time.time() - t0, 1), r.stderr[-500:], txt
    try:
        import z3
    except ImportError as e:
        raise RuntimeError('Z3 is required (z3 binary on PATH / $QOPT_Z3, or pip install z3-solver)') from e
    ans = []
    base = '\n'.join(['(set-logic QF_LIA)'] + list(lines))
    for _tag, t in checks:
        s = z3.Solver()
        s.set('timeout', int(timeout_s * 1000))
        s.from_string(base + f'\n(assert (not {t}))\n')
        ans.append({z3.sat: 'sat', z3.unsat: 'unsat'}.get(s.check(), 'unknown'))
    return ans, round(time.time() - t0, 1), '', txt
