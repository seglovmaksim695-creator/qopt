"""F2: Qrisp 0.9.10, mcx(method="jones") compiled with compile_mcm=True.
|+,+,0> -> mcx(jones) -> ordinary Toffoli (undoes it) -> H on both controls must give |000> with certainty.
Uses Qrisp's own simulator. Requires: pip install qrisp==0.9.10 (separate environment)."""
import warnings
warnings.filterwarnings("ignore")
from qrisp import QuantumVariable, mcx, h
from qrisp.simulator import run as qrun

for meth in ("jones", "gray"):
    for mcm in (False, True):
        qv = QuantumVariable(3)
        h(qv[0]); h(qv[1])
        mcx([qv[0], qv[1]], qv[2], method=meth)
        mcx([qv[0], qv[1]], qv[2], method="gray")
        h(qv[0]); h(qv[1])
        qc = qv.qs.compile(compile_mcm=mcm)
        for i in range(3):
            qc.measure(qc.qubits[i])
        print(f"method={meth:5s} compile_mcm={mcm!s:5s} counts:", qrun(qc.transpile(), 2000, None))
print("expected: only the all-zero outcome in every line")
