"""Metrics of a circuit compiled for a device: two-qubit gate count, depth, duration and the expected fidelity
computed from the calibration snapshot.

Expected fidelity (``fid``) = product of (1 - error) over all gates  x  product of (1 - readout error) over the
final positions of the logical qubits  x  T1/T2 damping of every used qubit over its idle time (as-soon-as-possible
schedule with the calibrated gate durations; a qubit starts to count at its first gate and all qubits are read out
together at the end).  ``esp`` is the same without the idle term.  This is the usual first-order estimate: it
ignores crosstalk and coherent errors, so it ranks circuits; it does not predict the absolute device result.
Check against 8 real IBM jobs (ibm_fez, ibm_marrakesh): the error was understated by a factor of about 1.6; the
ranking of variants was right for 93.7 % of the distinguishable pairs (Fez 100 %, Marrakesh 76.4 %, only 6 tasks).
Use it for a preliminary choice between variants, not for absolute fidelity numbers. No correction factor is applied.
"""
import math

SKIP = {"barrier", "delay", "measure", "id", "reset"}


def analyze(circ, device, final_map=None):
    """Metrics dict for a physical circuit ``circ`` (qubit index = physical qubit of ``device``).

    ``final_map``: physical qubit of every logical qubit at the end (where it is measured); default: every used
    qubit is measured."""
    t_q = {}
    first = {}
    busy = {}
    lg_gate = 0.0
    n2 = n1 = 0
    bad = 0
    twoq_names = {}
    for ins in circ.data:
        op = ins.operation
        name = op.name
        if name in SKIP:
            continue
        qs = [circ.find_bit(q).index for q in ins.qubits]
        if len(qs) == 1:
            q = qs[0]
            if name == "rz":
                dur, err = 0.0, 0.0
            else:
                dur, err = device.dur1[q], device.err1[q]
                n1 += 1
            start = t_q.get(q, 0.0)
        elif len(qs) == 2:
            a, b = qs
            k = (a, b) if a < b else (b, a)
            if k not in device.err2:
                bad += 1
                err, dur = 1.0, 0.0
            else:
                err, dur = device.err2[k], device.dur2[k]
            n2 += 1
            twoq_names[name] = twoq_names.get(name, 0) + 1
            start = max(t_q.get(a, 0.0), t_q.get(b, 0.0))
        else:
            raise ValueError(f"gate {name} on {len(qs)} qubits in a compiled circuit")
        lg_gate += math.log(max(1e-12, 1.0 - err))
        for q in qs:
            first.setdefault(q, start)
            busy[q] = busy.get(q, 0.0) + dur
            t_q[q] = start + dur
    t_end = max(t_q.values(), default=0.0)
    lg_idle = 0.0
    for q, t0 in first.items():
        lg_idle += device.idle_logfid(q, (t_end - t0) - busy[q])
    meas = list(final_map) if final_map is not None else sorted(first)
    lg_ro = sum(math.log(max(1e-12, 1.0 - device.ro_err[q])) for q in meas)
    return {
        "twoq": n2, "oneq": n1, "depth": circ.depth(lambda i: i.operation.name not in SKIP),
        "depth2q": circ.depth(lambda i: i.operation.num_qubits == 2 and i.operation.name not in SKIP),
        "duration_us": round(t_end * 1e6, 3), "used_qubits": len(first), "offmap_2q": bad,
        "log_gate": lg_gate, "log_readout": lg_ro, "log_idle": lg_idle,
        "esp": math.exp(lg_gate + lg_ro), "fid": math.exp(lg_gate + lg_ro + lg_idle),
        "twoq_names": twoq_names,
    }


def logical_stats(circ):
    """Gate statistics of a logical circuit (before compilation)."""
    ops = dict(circ.count_ops())
    n2 = sum(1 for i in circ.data if i.operation.num_qubits == 2 and i.operation.name not in SKIP)
    n3 = sum(1 for i in circ.data if i.operation.num_qubits >= 3 and i.operation.name not in SKIP)
    return {"qubits": circ.num_qubits, "size": circ.size(), "twoq": n2, "multiq": n3, "depth": circ.depth(), "ops": ops}
