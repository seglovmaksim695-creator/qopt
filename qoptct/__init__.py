"""qoptct -- a verified Clifford+T circuit optimizer (T count and CX count). Project name: qopt.

The distribution and the import package are called ``qoptct`` (the name ``qopt`` on PyPI belongs to an
unrelated project):  ``pip install qoptct``  /  ``import qoptct``.

Command line:  qoptct optimize INPUT.qc [-o OUT] [--mode cx|balanced|t] [--clean] [--meas] [--net] [--post] ...
Python:        import qoptct; names, ref, best, how = qoptct.optimize("circuit.qc")
Qiskit:        from qoptct.qiskit_qopt import QoptPass
               transpile(qc, optimization_method="qoptct")   # transpiler stage plugin

PyTorch is optional and only needed for the neural mode (``--net``):  pip install "qoptct[net]".

Internal modules use flat imports (``import evo_lib``) and ``import qoptct`` for the core functions; the package
directory is therefore put on sys.path when qoptct is imported, and the public names of ``qoptct.qopt_core``
(load_qc, save_qc, prove, optimize, ... and the time limits PROG_SEC, PROVE_SEC, FOLD_SEC, PYZX_SEC) are
re-exported here. Caches go to $QOPT_CACHE (default ~/.cache/qopt).
"""
import os as _os
import sys as _sys

__version__ = "0.6.0"
_PKG = _os.path.dirname(_os.path.abspath(__file__))
if _PKG not in _sys.path:
    _sys.path.insert(0, _PKG)



class _TopLevelNameGuard:
    """Keep the optional dependency ``qualtran`` importable.

    The package directory is on sys.path (see above), so the sub-package ``qoptct/qualtran`` would otherwise
    answer ``import qualtran``.  For the top-level name ``qualtran`` only, search sys.path without it.
    """

    @staticmethod
    def find_spec(name, path=None, target=None):
        if name != "qualtran" or path is not None:
            return None
        import importlib.machinery as _m
        return _m.PathFinder.find_spec(name, [p for p in _sys.path if _os.path.abspath(p or ".") != _PKG])


if _TopLevelNameGuard not in _sys.meta_path:
    _sys.meta_path.insert(0, _TopLevelNameGuard)

from . import qopt_core as _core  # noqa: E402

globals().update({k: v for k, v in vars(_core).items() if not k.startswith("__") and k not in ("_os", "_sys")})
