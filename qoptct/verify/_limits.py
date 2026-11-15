# SPDX-License-Identifier: Apache-2.0
"""Run a function in a separate Python process with a hard wall-clock limit (the process is killed on timeout).

A fresh interpreter is started with `python -m qoptct.verify._worker` (not multiprocessing), so the caller's
__main__ module is never re-imported: prove() is safe to call from scripts without an `if __name__ == "__main__"`
guard, from daemonic workers and from notebooks."""
import os
import pickle
import subprocess
import sys
import traceback

_PKG_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _inline(fn, args, kwargs):
    try:
        return ('ok', fn(*args, **kwargs))
    except Exception as e:
        return ('err', '%s: %s' % (type(e).__name__, e), traceback.format_exc(limit=5))


def run_limited(fn, args=(), kwargs=None, budget_sec=None):
    """-> ('ok', value) | ('err', message, traceback) | ('timeout', None).
    fn must be a module-level function of an importable module; args/kwargs/result must be picklable.
    budget_sec None or <= 0: run in-process without a hard limit."""
    kwargs = kwargs or {}
    if budget_sec is None or budget_sec <= 0 or os.environ.get('QOPT_VERIFY_INPROCESS'):
        return _inline(fn, args, kwargs)
    payload = pickle.dumps((fn.__module__, fn.__qualname__, tuple(args), dict(kwargs)))
    env = dict(os.environ)
    env['PYTHONPATH'] = _PKG_PARENT + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    try:
        p = subprocess.Popen([sys.executable, '-m', 'qoptct.verify._worker'], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    except OSError:
        return _inline(fn, args, kwargs)
    try:
        out, err = p.communicate(payload, timeout=budget_sec)
    except subprocess.TimeoutExpired:
        p.kill()
        p.communicate()
        return ('timeout', None)
    try:
        return pickle.loads(out)
    except Exception:
        tail = err.decode('utf-8', 'replace')[-400:] if err else ''
        return ('err', 'worker failed (exit code %s): %s' % (p.returncode, tail), '')
