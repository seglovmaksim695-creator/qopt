"""F4: pytket ZXGraphlikeOptimisation leaves an implicit wire permutation.
Exporting the circuit without replace_implicit_wire_swaps() gives a different unitary.
Requires: pip install pytket==2.18.5 (separate environment), numpy."""
import numpy as np
from pytket import Circuit, OpType
from pytket.passes import AutoRebase, ZXGraphlikeOptimisation

# tof_3-like Toffoli network on 5 qubits, written in Clifford+T
c = Circuit(5)
for a, b, t in ((0, 1, 4), (2, 4, 3), (0, 1, 4)):
    c.CCX(a, b, t)
ref = c.get_unitary()

opt = c.copy()
AutoRebase({OpType.CX, OpType.Rz, OpType.H}).apply(opt)
ZXGraphlikeOptimisation().apply(opt)
print("implicit permutation after ZX pass:", opt.implicit_qubit_permutation())

naive = Circuit(opt.n_qubits)          # what a gate-by-gate QASM export contains
for cmd in opt.get_commands():
    naive.add_gate(cmd.op, cmd.args)
fixed = opt.copy(); fixed.replace_implicit_wire_swaps()


def same(u, v):
    k = np.argmax(np.abs(u.ravel()))
    ph = v.ravel()[k] / u.ravel()[k]
    return abs(abs(ph) - 1) < 1e-9 and np.allclose(u * ph, v, atol=1e-8)


print("gate-by-gate export equals the input:", same(ref, naive.get_unitary()))
print("after replace_implicit_wire_swaps():  ", same(ref, fixed.get_unitary()))
