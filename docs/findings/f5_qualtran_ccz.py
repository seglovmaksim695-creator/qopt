"""F5: Qualtran 0.7.0 -- the T cost of an imported Cirq gate depends on which Toffoli-type gate is used.
cirq.CCX / cirq.TOFFOLI become Qualtran's Toffoli (4 T in Qualtran's cost model); cirq.CCZ is decomposed by Cirq (7 T).
Requires: pip install qualtran==0.7.0 cirq (separate environment)."""
import cirq
from qualtran.cirq_interop import cirq_optree_to_cbloq
from qualtran.resource_counting import QECGatesCost, get_cost_value

q = cirq.LineQubit.range(3)
for name, op in (("CCX", cirq.CCX(*q)), ("TOFFOLI", cirq.TOFFOLI(*q)), ("CCZ", cirq.CCZ(*q))):
    gc = get_cost_value(cirq_optree_to_cbloq(cirq.Circuit(op)), QECGatesCost())
    print(f"cirq.{name:8s} -> {gc}  total_t_count = {gc.total_t_count()}")
