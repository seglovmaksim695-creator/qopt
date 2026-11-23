"""Patched copies of Qualtran 0.7.0, exactly as the proofs were made.

The proofs of Montgomery multiplication and of point addition were done with Qualtran
0.7.0 plus small fixes.  The installed Qualtran is never modified.  Instead
``patched_qualtran(kind)`` builds (once) a COPY of the installed package in ``$QOPT_CACHE``
and overlays the patched files on it.  Put the returned directory first in ``PYTHONPATH`` of
a fresh interpreter (or use ``patched_env`` / ``run_with_patched_qualtran`` /
``use_patched_qualtran``).

Kinds
-----
``'c'``   PR 1997 (the final reduction of the Montgomery multiplier subtracts p instead of
          adding it) and QROM/Montgomery tables as exact Python integers (no int64 overflow
          for widths above 63 bits).  Used for the Montgomery multiplication proofs.
          Files: data_loading/qrom_base.py, data_loading/select_swap_qrom.py,
          mod_arithmetic/mod_multiplication.py.
``'e'``   the ECAdd flag fix (finding E05) plus PR 1997 and the exact-integer Montgomery
          table.  Used for the point addition proofs at n=16.
          Files: cryptography/ecc/ec_add.py, mod_arithmetic/mod_multiplication.py.
``'e4'``  kind ``'e'`` plus the QROM files of kind ``'c'``.  Used for point addition at n=256.

Every patched file has a unified diff in ``patches/diff`` (for reading) and the full file in
``patches/src`` (what is applied).  ``patches/manifest.json`` stores the SHA-256 of each
original 0.7.0 file and of each patched file; a mismatch of an original file (another
Qualtran version or a modified installation) is refused.
"""
from __future__ import annotations

import contextlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from hashlib import sha256
from importlib.metadata import version as _version
from pathlib import Path

from .._cache import cache_dir

__all__ = ['KINDS', 'patched_qualtran', 'patched_env', 'run_with_patched_qualtran',
           'use_patched_qualtran', 'verify_copy', 'PatchError']

_HERE = Path(__file__).resolve().parent
_MANIFEST = json.loads((_HERE / 'manifest.json').read_text())
KINDS = tuple(_MANIFEST['kinds'])


class PatchError(RuntimeError):
    pass


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _installed_root():
    spec = importlib.util.find_spec('qualtran')  # does not import qualtran
    if spec is None or not spec.submodule_search_locations:
        raise PatchError('Qualtran is not installed (pip install "qoptct[qualtran]")')
    return Path(list(spec.submodule_search_locations)[0])


def _plan(kind):
    if kind not in _MANIFEST['kinds']:
        raise PatchError(f'unknown patch kind {kind!r}; known: {KINDS}')
    return [(name, _MANIFEST['files'][name]) for name in _MANIFEST['kinds'][kind]['files']]


def verify_copy(kind, root):
    """Check that ``root`` (a directory holding ``qualtran/``) has exactly the patched files of ``kind``.

    Returns a dict ``{relative target: sha256}`` of the patched files; raises ``PatchError`` on any mismatch.
    """
    root = Path(root)
    out = {}
    for name, info in _plan(kind):
        path = root / 'qualtran' / info['target']
        if not path.is_file() or _sha(path) != info['sha256_patched']:
            raise PatchError(f'{info["target"]}: the copy does not carry patch {name}')
        out[info['target']] = info['sha256_patched']
    return out


def patched_qualtran(kind):
    """Return the directory (containing ``qualtran/``) of the patched copy for ``kind``; build it if needed.

    The copy lives in ``$QOPT_CACHE/qualtran_patched/<kind>-<version>/``.  Files that are patched are
    checked against the manifest; the installed Qualtran must be the version the patches were made for
    and its original files must hash as expected.
    """
    plan = _plan(kind)
    installed = _installed_root()
    ver = _version('qualtran')
    if ver != _MANIFEST['qualtran_version']:
        raise PatchError(f'patches are made for Qualtran {_MANIFEST["qualtran_version"]}, installed {ver}')
    for name, info in plan:
        if _sha(installed / info['target']) != info['sha256_original']:
            raise PatchError(f'{info["target"]}: the installed file is not the original Qualtran {ver} file')
    base = cache_dir('qualtran_patched') / f'{kind}-{ver}'
    try:
        verify_copy(kind, base)
        return base
    except PatchError:
        pass
    if base.exists():
        shutil.rmtree(base)
    tmp = base.with_name(base.name + f'.tmp{os.getpid()}')
    if tmp.exists():
        shutil.rmtree(tmp)
    shutil.copytree(installed, tmp / 'qualtran', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name, info in plan:
        shutil.copyfile(_HERE / 'src' / name, tmp / 'qualtran' / info['target'])
    verify_copy(kind, tmp)
    tmp.rename(base)
    return base


def patched_env(kind, env=None):
    """A copy of ``env`` (default ``os.environ``) whose ``PYTHONPATH`` starts with the patched copy."""
    env = dict(os.environ if env is None else env)
    path = str(patched_qualtran(kind))
    env['PYTHONPATH'] = path + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    return env


def run_with_patched_qualtran(kind, argv, **kwargs):
    """``subprocess.run(argv, env=patched_env(kind))``; ``argv`` is a list (e.g. ``[sys.executable, 'x.py']``)."""
    env = patched_env(kind, kwargs.pop('env', None))
    return subprocess.run(argv, env=env, **kwargs)


@contextlib.contextmanager
def use_patched_qualtran(kind):
    """In-process variant: only valid BEFORE ``qualtran`` is imported in this interpreter.

    Prefer a fresh process (``run_with_patched_qualtran``).  Raises if Qualtran is already imported from
    elsewhere, because modules loaded from the installed copy cannot be swapped safely.
    """
    path = str(patched_qualtran(kind))
    loaded = sys.modules.get('qualtran')
    if loaded is not None and not str(Path(loaded.__file__).resolve()).startswith(path):
        raise PatchError('qualtran is already imported from the installed copy; start a fresh process')
    sys.path.insert(0, path)
    try:
        yield Path(path)
    finally:
        with contextlib.suppress(ValueError):
            sys.path.remove(path)
