"""qoptct.hw -- compilation for noisy IBM-like hardware (heavy-hex coupling, native CZ/ECR + RZ, SX, X).

    from qoptct.hw import compile_hw
    res = compile_hw(qc, "fake_torino")          # res.circuit, res.metrics, res.table, res.verified

Command line:  qoptct hw --backend fake_torino in.qasm -o out.qasm

Needs qiskit-ibm-runtime for the fake backends (calibration snapshots) and, optionally, qiskit-aer for the noisy
simulation:  pip install "qoptct[hw]".  Nothing here ever submits a job to real hardware.
"""
from .device import Device, get_backend  # noqa: F401


def __getattr__(name):
    if name in ("compile_hw", "HwResult", "qiskit_o3"):
        from . import compiler
        return getattr(compiler, name)
    if name in ("analyze",):
        from . import metrics
        return getattr(metrics, name)
    if name in ("verify_compiled",):
        from . import check
        return getattr(check, name)
    raise AttributeError(name)
