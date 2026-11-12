"""qopt — a verified Clifford+T circuit optimizer (T count and CX count).

Command line:  qopt optimize INPUT.qc [-o OUT] [--mode cx|balanced|t] [--clean] [--meas] [--net] [--post] ...
Python:        import qopt; names, ref, best, how = qopt.optimize("circuit.qc")
Qiskit:        from qopt.qiskit_qopt import QoptPass

Internal modules use flat imports (``import evo_lib``) and ``import qopt`` for the core functions; the package
directory is therefore put on sys.path when qopt is imported, and the public names of ``qopt.qopt_core``
(load_qc, save_qc, prove, optimize, ... and the time limits PROG_SEC, PROVE_SEC, FOLD_SEC, PYZX_SEC) are
re-exported here. Caches go to $QOPT_CACHE (default ~/.cache/qopt).
"""
import os as _os
import sys as _sys

__version__ = "0.2.0"
_PKG = _os.path.dirname(_os.path.abspath(__file__))
if _PKG not in _sys.path:
    _sys.path.insert(0, _PKG)

from . import qopt_core as _core  # noqa: E402

globals().update({k: v for k, v in vars(_core).items() if not k.startswith("__") and k not in ("_os", "_sys")})
