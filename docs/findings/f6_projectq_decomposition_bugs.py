"""Minimal reproduction, ProjectQ 0.8.0 only (no qoptct code).  The gate is run twice on ProjectQ's ClassicalSimulator:
(a) as the gate itself (its math function), (b) after ProjectQ's own decomposition rules into X/CNOT/Toffoli.
They must agree; for ComparatorQuantum (n = 1 and n >= 3) and MultiplyQuantum (b >= 2) they do not."""
import projectq.libs.math
import projectq.setups.decompositions
from projectq import MainEngine
from projectq.backends import ClassicalSimulator
from projectq.cengines import AutoReplacer, DecompositionRuleSet, InstructionFilter
from projectq.libs.math import ComparatorQuantum, MultiplyQuantum
from projectq.ops import XGate


def run(gate, widths, values, decompose):
    chain = []
    if decompose:
        rules = DecompositionRuleSet(modules=[projectq.libs.math, projectq.setups.decompositions])
        chain = [AutoReplacer(rules), InstructionFilter(lambda eng, cmd: isinstance(cmd.gate, XGate) and len(cmd.control_qubits) <= 2
                 or cmd.gate.__class__.__name__ in ('AllocateQubitGate', 'DeallocateQubitGate', 'FlushGate'))]
    eng = MainEngine(ClassicalSimulator(), engine_list=chain)
    regs = [eng.allocate_qureg(w) for w in widths]
    eng.flush()
    for r, v in zip(regs, values):
        eng.backend.write_register(r, v)
    gate | tuple(regs)
    eng.flush()
    return [eng.backend.read_register(r) for r in regs]


# ComparatorQuantum | (a, b, c): c ^= (a > b), a and b unchanged.  n = 3, a = 1, b = 0:
print('comparator n=3 gate     :', run(ComparatorQuantum, [3, 3, 1], [1, 0, 0], False), 'expected [1, 0, 1]')
print('comparator n=3 decomposed:', run(ComparatorQuantum, [3, 3, 1], [1, 0, 0], True))
# MultiplyQuantum | (a, b, c) with c of 2n+1 qubits, c = 0: c = a*b.  n = 2, a = 2, b = 2:
print('multiply n=2 gate       :', run(MultiplyQuantum, [2, 2, 5], [2, 2, 0], False), 'expected [2, 2, 4]')
print('multiply n=2 decomposed :', run(MultiplyQuantum, [2, 2, 5], [2, 2, 0], True))
