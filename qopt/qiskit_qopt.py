"""Плагин Qiskit: проход транспилятора QoptPass — проверяемая оптимизация CX/T схем Clifford+T.
    from qiskit.transpiler import PassManager
    from qiskit_qopt import QoptPass
    out = PassManager([QoptPass()]).run(qc)
Поддерживаемые команды: h x z s sdg t tdg cx cz ccx ccz (иначе схема возвращается без изменений).
Результат доказан (как в qopt.py); если доказать не удалось — возвращается исходная схема."""
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.transpiler.basepasses import TransformationPass
import qopt

SUP = {"h", "x", "z", "s", "sdg", "t", "tdg", "cx", "cz", "ccx", "ccz"}


def to_hl(qc):
    hl = []
    for ins in qc.data:
        g = ins.operation.name
        if g not in SUP or ins.clbits:
            return None
        hl.append((g, tuple(qc.find_bit(q).index for q in ins.qubits)))
    return hl


def from_low(n, c):
    qc = QuantumCircuit(n)
    for g, q in c:
        getattr(qc, g)(*q)
    return qc


class QoptPass(TransformationPass):
    def __init__(self, net=False):
        super().__init__(); self.net = net; self.last = None

    def run(self, dag):
        qc = dag_to_circuit(dag); hl = to_hl(qc)
        if hl is None or not hl:
            self.last = {"skipped": "неподдерживаемые команды"}; return dag
        names = [str(i) for i in range(qc.num_qubits)]
        _, ref, best, how = qopt.optimize_hl(names, hl, self.net)
        self.last = {"method": how, "verified": how is not None}
        if how is None:
            return dag
        out = from_low(qc.num_qubits, best)
        new = circuit_to_dag(out)
        # сохранить исходные регистры: копируем команды на исходные кубиты
        res = dag.copy_empty_like()
        res.compose(new, qubits=list(dag.qubits))
        return res


class QoptMeasPass(TransformationPass):
    """Режим с измерениями: логическое И Гидни. ancillas — индексы кубитов, которые на входе гарантированно |0>.
    Тоффоли в чистый кубит -> 4 T; стирание известного a·b -> measure + if_test(cz, x). Доказательство: GF(2) или
    симуляция с измерениями (gidney.py). Если не доказано или T не меньше — схема не меняется."""

    def __init__(self, ancillas):
        super().__init__(); self.anc = set(ancillas); self.last = None

    def run(self, dag):
        import gidney, evo_lib
        from qiskit import ClassicalRegister
        qc = dag_to_circuit(dag); hl = to_hl(qc)
        if not hl:
            self.last = {"skipped": "неподдерживаемые команды"}; return dag
        n = qc.num_qubits; hl = evo_lib.fold_h(hl)
        items, r1, r2 = gidney.transform(hl, n, self.anc)
        if not (r1 or r2):
            self.last = {"skipped": "нет мест"}; return dag
        import meas_portfolio
        names = [str(i) for i in range(n)]
        low = meas_portfolio.build_low(names, items, n)
        ok = gidney.anf_check(items, n, self.anc)
        if ok is None:
            ok = gidney.verify(to_hl(qc), low, n, self.anc)
        T = sum(evo_lib.tcount(s) for s, _ in low); T0 = sum(g in ("t", "tdg") for g, _ in evo_lib.lower_std(to_hl(qc)))
        self.last = {"t_before": T0, "t_after": T, "verified": ok, "measurements": r2}
        if ok is not True or T >= T0:
            return dag
        out = QuantumCircuit(n); creg = ClassicalRegister(max(r2, 1), "m"); out.add_register(creg); k = 0
        for seg, mz in low:
            for g, q in seg: getattr(out, g)(*q)
            if mz:
                a, b, t = mz; out.measure(t, creg[k])
                with out.if_test((creg[k], 1)):
                    out.cz(a, b); out.x(t)
                k += 1
        return circuit_to_dag(out)
