"""ProjectQ -> qoptct: a strict gate-by-gate bridge and proof clients for ProjectQ programs.

ProjectQ is an OPTIONAL dependency (``pip install "qoptct[projectq]"``).  Importing this package never imports ProjectQ;
the modules that need it import it lazily.

Main entry points
-----------------
``check_projectq(program, spec, ...)``
    Run a ProjectQ program through ProjectQ's own decomposition engine into NCT (X, CNOT, Toffoli) or Clifford+T,
    translate it command by command, and prove it against ``spec`` (``ClassicalSpec`` for a whole-domain classical
    function, ``GF2MulSpec`` for GF(2^m) multipliers at any size).  Verdict: PROVED / REFUTED (with a witness) / UNKNOWN.
``record(program, mode)`` / ``to_logical(recording)``
    The recording and the logical qoptct program (x, cx, and, erase, Clifford+T atoms) with its OpenQASM text.
``registry``
    Proved results about the Binary_ECC circuits and ProjectQ library blocks (``qoptct projectq results``).
"""

__all__ = ['check_projectq', 'ClassicalSpec', 'GF2MulSpec', 'GarbagePin', 'pin_garbage', 'record', 'to_logical', 'Recording',
           'Logical', 'Unsupported', 'ProjectQUnavailable']

_LAZY = {
    'check_projectq': 'check', 'ClassicalSpec': 'check', 'GF2MulSpec': 'check', 'GarbagePin': 'check', 'pin_garbage': 'check',
    'record': '_record', 'Recording': '_record', 'Unsupported': '_record', 'ProjectQUnavailable': '_record',
    'to_logical': '_logical', 'Logical': '_logical',
}


def __getattr__(name):
    if name in _LAZY:
        import importlib
        return getattr(importlib.import_module('.' + _LAZY[name], __name__), name)
    raise AttributeError(name)
