"""Which Qualtran is imported, and how to run a module against a patched copy of it.

``detect_variant()`` classifies the Qualtran that is *currently imported* by the SHA-256 of the four files the
patches touch: ``'original'`` (pristine 0.7.0), one of the patch kinds ``'c'``, ``'e'``, ``'e4'`` (see
``qoptct.foreign.qualtran.patches``) or ``'other'`` (anything else; every proof client then refuses to proceed).

``run_module(kind, module, args)`` starts ``python -m module args`` in a fresh interpreter whose ``PYTHONPATH``
starts with the patched copy for ``kind`` (``None`` = the installed, unpatched Qualtran).  The installed Qualtran
is never modified.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

__all__ = ['detect_variant', 'expected_kind_name', 'repo_root', 'child_env', 'run_module', 'run_json']

_THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMBA_NUM_THREADS')


def repo_root():
    """Directory that contains the ``qoptct`` package (must be importable in a child process)."""
    import qoptct
    return Path(qoptct.__file__).resolve().parent.parent


def _manifest():
    from . import patches
    return patches._MANIFEST, patches._sha


def expected_kind_name(kind):
    """Name used in results for a requested kind (``None`` -> ``'original'``)."""
    return 'original' if kind is None else kind


def detect_variant():
    """Classify the imported Qualtran.  Returns a dict with ``kind`` and the SHA-256 (16 hex) of the touched files."""
    import qualtran
    man, sha = _manifest()
    root = Path(qualtran.__file__).resolve().parent
    targets = {}
    for name, info in man['files'].items():
        targets.setdefault(info['target'], []).append((name, info))
    actual = {}
    for target in targets:
        path = root / target
        actual[target] = sha(path) if path.is_file() else None
    files = {t: (h or 'missing')[:16] for t, h in actual.items()}
    kind = 'other'
    if all(actual[t] == targets[t][0][1]['sha256_original'] for t in targets):
        kind = 'original'
    else:
        for k, spec in man['kinds'].items():
            want = {t: targets[t][0][1]['sha256_original'] for t in targets}
            for name in spec['files']:
                info = man['files'][name]
                want[info['target']] = info['sha256_patched']
            if all(actual[t] == want[t] for t in targets):
                kind = k
                break
    return dict(kind=kind, qualtran_file=str(Path(qualtran.__file__).resolve()),
                qualtran_version=getattr(qualtran, '__version__', None), files_sha256_16=files)


def child_env(kind, extra=None):
    """Environment of a child process: patched copy first in PYTHONPATH, single-threaded numerics."""
    from .patches import patched_qualtran
    env = dict(os.environ)
    parts = []
    if kind is not None:
        parts.append(str(patched_qualtran(kind)))
    parts.append(str(repo_root()))
    if env.get('PYTHONPATH'):
        parts.append(env['PYTHONPATH'])
    env['PYTHONPATH'] = os.pathsep.join(parts)
    env.setdefault('PYTHONDONTWRITEBYTECODE', '1')
    for k in _THREADS:
        env.setdefault(k, '1')
    if extra:
        env.update(extra)
    return env


def run_module(kind, module, args, *, timeout=None, env_extra=None):
    """Run ``python -m module *args`` against Qualtran ``kind``; returns the CompletedProcess (text output)."""
    argv = [sys.executable, '-B', '-s', '-m', module] + [str(a) for a in args]
    return subprocess.run(argv, env=child_env(kind, env_extra), capture_output=True, text=True, timeout=timeout)


def run_json(kind, module, args, *, timeout=None, env_extra=None):
    """Like ``run_module`` but the module must accept ``--out FILE``; returns the parsed JSON written there.

    A crash, a missing or unreadable file gives ``{'status': 'ERROR', ...}`` (never a verdict of success).
    """
    with tempfile.TemporaryDirectory(prefix='qoptct_') as td:
        out = os.path.join(td, 'out.json')
        try:
            r = run_module(kind, module, list(args) + ['--out', out], timeout=timeout, env_extra=env_extra)
        except subprocess.TimeoutExpired:
            return dict(status='UNKNOWN', detail=f'timeout after {timeout} s', module=module)
        try:
            with open(out) as f:
                return json.load(f)
        except Exception as e:  # noqa: BLE001
            return dict(status='ERROR', detail=f'no result file ({type(e).__name__}); exit code {r.returncode}',
                        stderr=(r.stderr or '')[-1500:], stdout=(r.stdout or '')[-500:], module=module)
