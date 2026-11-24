# SPDX-License-Identifier: Apache-2.0
"""Cost vector of an OpenQASM 3 circuit written by qoptct, and the Pareto order used by the loop."""
import re

KEYS = ('AND', 'T', 'CX', 'qubits')


def qasm_cost(txt):
    """AND (= T/4, as in results/crypto/index.json), T (t + tdg), CX (cx + cz), qubits (highest wire used + 1)."""
    t = len(re.findall(r'^(?:t|tdg) q', txt, re.M))
    cx = len(re.findall(r'\bcx q', txt))
    cz = len(re.findall(r'\bcz q', txt))
    used = [int(x) for x in re.findall(r'q\[(\d+)\]', txt)]
    return {'AND': t // 4, 'T': t, 'CX': cx + cz, 'qubits': (max(used) + 1) if used else 0}


def dominates(a, b):
    """a is no worse than b in every key and better in at least one."""
    return all(a[k] <= b[k] for k in KEYS) and any(a[k] < b[k] for k in KEYS)


def sort_key(c):
    return (c['AND'], c['CX'], c['qubits'], c['T'])


def front_insert(front, cost):
    """Pareto front update. Returns (new_front, accepted). A cost equal to a member or dominated is not accepted."""
    for f in front:
        if f == cost or dominates(f, cost):
            return list(front), False
    return [f for f in front if not dominates(cost, f)] + [dict(cost)], True


def best_of(front):
    return min(front, key=sort_key) if front else None


def pad_qubits(txt, n):
    """Raise the declared register size to at least n (the reference may have more wires than the circuit uses)."""
    m = re.search(r'qubit\[(\d+)\] q;', txt)
    if m and int(m.group(1)) < n:
        txt = txt[:m.start()] + 'qubit[%d] q;' % n + txt[m.end():]
    return txt
