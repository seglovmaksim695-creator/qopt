"""Minimal reproduction, Qiskit only: VBERippleCarryAdder(kind='fixed') is not (a+b) mod 2^n."""
import numpy as np, qiskit
from qiskit.circuit.library import VBERippleCarryAdder
from qiskit.quantum_info import Statevector
print('qiskit', qiskit.__version__)
for n in (2, 3, 4):
    qc = VBERippleCarryAdder(n, kind='fixed')      # registers a(n), b(n), helper(n-1)
    bad = []
    for a in range(1 << n):
        for b in range(1 << n):
            sv = Statevector.from_int(a | (b << n), 2 ** qc.num_qubits).evolve(qc)
            k = int(np.argmax(np.abs(sv.data)))
            if (k >> n) & ((1 << n) - 1) != (a + b) % (1 << n):
                bad.append((a, b, (k >> n) & ((1 << n) - 1)))
    print('n=%d: %d of %d inputs give b != (a+b) mod 2^n; first: (a,b,got) = %s' % (n, len(bad), 4 ** n, bad[:3]))
