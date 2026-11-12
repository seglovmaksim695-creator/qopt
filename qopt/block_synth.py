"""Этап 11: блочный синтез. Схему на n кубитах режем на блоки по <=k (4/3) кубитов,
каждый блок пересобирают обученные агенты синтеза 2q/3q/4q (exp/u_k_best), берём, если CX меньше.
Проходы повторяются (k=4, k=3, и на обращённой схеме — другие границы блоков), пока есть выигрыш.
Итог проверяется точным сравнением матриц (Qiskit Operator.equiv)."""
import argparse, json, os, time, warnings
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import torch
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent
INV = {"h": "h", "x": "x", "s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t", "cx": "cx"}
CTX = {}
CACHE = {}
W = {}


def _init(width, widen, model4='', no_net=False, opt_time=0.5, uniform=False):
    import qgnn
    from qfast import ChoiEngine
    torch.set_num_threads(1)
    W["w"], W["widen"], W["no_net"], W["opt"], W["uniform"] = width, widen, no_net, opt_time, uniform
    for k in (2, 3, 4) + ((5,) if os.environ.get("MODEL5") else ()):
        mp = {4: model4, 3: os.environ.get("MODEL3", ""), 5: os.environ.get("MODEL5", "")}.get(k, "")
        ck = torch.load((ROOT / mp) if mp else (ROOT / "models" / f"u_{k}_best_{k}q.pt"), map_location="cpu", weights_only=False)
        net = {"gnnu": qgnn.QGNNU, "gnnh": qgnn.QGNNH, "gnnhp": qgnn.QGNNHP}[ck.get("arch", "gnnu")](**ck["gnn_kw"])
        net.load_state_dict(ck["state_dict"]); net.eval()
        eng = ChoiEngine(k)
        cfg = {"beta": 3.0, "gamma": 0.0, "max_steps": ck["max_steps"], "extra": 3, "task": "unitary"}
        val = None
        vp = ROOT / "models" / f"value_u{k}q.pt"
        if vp.exists():
            vck = torch.load(vp, map_location="cpu", weights_only=False)
            val = qgnn.QValueU(**vck["kw"]); val.load_state_dict(vck["state_dict"]); val.eval(); cfg["lam"] = 0.5
        gi = {(nm, tuple(qb)): i for i, (nm, qb) in enumerate(eng.gates)}
        CTX[k] = (net, eng, cfg, val, gi)


def key(c):
    return (sum(g[0] == "cx" for g in c), len(c))


def inv_circ(c):
    return [(INV[nm], qb) for nm, qb in reversed(c)]


DIAG = {"s", "sdg", "t", "tdg"}


def commutes(g1, g2):
    """Точные правила перестановки команд Clifford+T (И4): диагональные — через управляющий CX, X — через цель,
    CX с общими только управляющими или только целевыми кубитами."""
    (n1, q1), (n2, q2) = g1, g2
    if not set(q1) & set(q2):
        return True
    if n1 != "cx" and n2 != "cx":
        return (n1 in DIAG and n2 in DIAG) or (n1 == n2 == "x") or (n1 == n2 == "h")
    if n1 == "cx" and n2 == "cx":
        return q1[0] != q2[1] and q2[0] != q1[1]
    one, cx = (g1, g2) if n2 == "cx" else (g2, g1)
    return (one[0] in DIAG and one[1][0] == cx[1][0]) or (one[0] == "x" and one[1][0] == cx[1][1])


def take_block_commute(rem, k, rng=None, p=0.0):
    """Как take_block, но команда может «перепрыгнуть» пропущенные, если коммутирует со всеми пропущенными на своих кубитах."""
    Q, take, skip = set(), [], []
    for g in rem:
        qs = set(g[1])
        if rng is not None and Q and rng.random() < p:
            skip.append(g); continue
        if len(Q | qs) <= k and all(commutes(g, s) for s in skip if set(s[1]) & qs):
            take.append(g); Q |= qs
        else:
            skip.append(g)
    return sorted(Q), take, skip


def take_block(rem, k, rng=None, p=0.0):
    if os.environ.get("COMMUTE"):
        return take_block_commute(rem, k, rng, p)
    Q, B, take, skip = set(), set(), [], []
    for g in rem:
        qs = set(g[1])
        if rng is not None and Q and rng.random() < p:  # случайный пропуск — другие границы блоков
            skip.append(g); B |= qs; continue
        if not (qs & B) and len(Q | qs) <= k:
            take.append(g); Q |= qs
        else:
            skip.append(g); B |= qs
    return sorted(Q), take, skip


def resynth(qs, gates):
    import ei2
    from circuit_opt import optimize
    k = len(qs)
    if k < 2 or not any(g[0] == "cx" for g in gates):
        return gates
    net, eng, cfg, val, gi = CTX[k]
    m = {q: i for i, q in enumerate(qs)}
    prog = [gi[(nm, tuple(m[q] for q in qb))] for nm, qb in gates]
    target = eng.choi(prog)
    from qfast import state_sig
    sig = (k, state_sig(target))
    if W.get("use_cache", True) and sig in CACHE:  # память: этот блок уже решали
        W["hits"] = W.get("hits", 0) + 1
        best = CACHE[sig]
        return [(eng.gates[eng.inverse_action[a]][0], tuple(qs[i] for i in eng.gates[eng.inverse_action[a]][1]))
                for a in reversed(best)]
    raw = [eng.inverse_action[a] for a in reversed(prog)]  # известная размотка (точная)
    best = optimize(eng, target, raw, time_limit=W["opt"])
    if not exact(eng, target, best):
        best = raw
    ei2._net, ei2._eng, ei2._cfg, ei2._val = net, eng, cfg, val
    if W.get("uniform"):  # перебор без сети: равные вероятности всех команд, без оценщика
        ei2._net, ei2._val = Uniform(eng.n_actions), None
    for w in (() if W["no_net"] else (W["w"], W["widen"])):
        if not w:
            continue
        res = ei2.beam(target, w, bound=len(best))
        if res is not None:
            acts = optimize(eng, target, res[1], time_limit=1.0)
            if exact(eng, target, acts) and key([eng.gates[a] for a in acts]) < key([eng.gates[a] for a in best]):
                best = acts
            break
    CACHE[sig] = list(best)
    return [(eng.gates[eng.inverse_action[a]][0], tuple(qs[i] for i in eng.gates[eng.inverse_action[a]][1]))
            for a in reversed(best)]


class Uniform(torch.nn.Module):
    def __init__(self, na):
        super().__init__(); self.na = na

    def forward(self, x):
        return torch.zeros(x.shape[0], self.na)


def exact(eng, target, acts):
    """Точная проверка блока: размотка должна дать тождество (совпадение 1 с точностью 1e-9), а не >= 0,99."""
    s = target
    for a in acts:
        s = eng.apply1(s, a)
    return abs(eng.amp0(s)) ** 2 >= 1 - 1e-9


def one_pass(circ, k, rng=None, p=0.0):
    rem, out = circ, []
    while rem:
        qs, take, rem = take_block(rem, k, rng, p)
        out += resynth(qs, take)
    return out


def to_qc(c, n):
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(n)
    for nm, qb in c:
        getattr(qc, nm)(*qb)
    return qc


def solve(args):
    idx, n, circ, max_cycles = args
    from qiskit.quantum_info import Operator
    t0 = time.time()
    best = circ
    for _ in range(max_cycles):
        improved = False
        for k in (4, 3):
            for flip in (0, 1):
                c = inv_circ(best) if flip else best
                c2 = one_pass(c, min(k, n))
                c2 = inv_circ(c2) if flip else c2
                if key(c2) < key(best):
                    best, improved = c2, True
        if not improved:
            break
    ok = Operator(to_qc(best, n)).equiv(Operator(to_qc(circ, n)))
    return idx, {"cx": key(best)[0], "len": len(best), "equiv": bool(ok), "sec": time.time() - t0, "circ": best}


def gen(n, count, seed0=9_000_000):
    from qiskit import transpile
    from qfast import Engine, ChoiEngine
    eng = Engine(n)
    ce = ChoiEngine(n) if n <= 5 else None
    one = [i for i, (nm, _) in enumerate(eng.gates) if nm != "cx"]
    cxs = [i for i, (nm, _) in enumerate(eng.gates) if nm == "cx"]
    out, seed = [], seed0
    while len(out) < count:
        rng = np.random.default_rng(seed); seed += 1
        prog = []
        for _ in range(int(rng.integers(n, 4 * n + 1))):
            for _ in range(int(rng.integers(1, n + 1))):
                prog.append(int(rng.choice(one)))
            prog.append(int(rng.choice(cxs)))
        if ce is not None and abs(ce.amp0(ce.choi(prog))) ** 2 >= 0.99:
            continue
        c = [(eng.gates[a][0], tuple(eng.gates[a][1])) for a in prog]
        qcx = transpile(to_qc(c, n), basis_gates=["u", "cx"], optimization_level=3,
                        seed_transpiler=0).count_ops().get("cx", 0)
        if qcx == 0:
            continue
        out.append((c, qcx))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--qubits", type=int, default=5)
    p.add_argument("--count", type=int, default=300)
    p.add_argument("--width", type=int, default=64)
    p.add_argument("--widen", type=int, default=0)
    p.add_argument("--cycles", type=int, default=3)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--model4", default="")
    p.add_argument("--no-net", action="store_true")
    p.add_argument("--uniform", action="store_true", help="тот же поиск, но без сети (равные вероятности)")
    p.add_argument("--opt-time", type=float, default=0.5)
    p.add_argument("--tag", default="")
    a = p.parse_args()
    n = a.qubits
    t0 = time.time()
    data = gen(n, a.count)
    print(f"{n}q: {len(data)} схем, CX в исходной {np.mean([key(c)[0] for c, _ in data]):.1f}, "
          f"Qiskit opt3 {np.mean([q for _, q in data]):.1f} ({time.time()-t0:.0f} с)", flush=True)
    with Pool(a.workers, initializer=_init, initargs=(a.width, a.widen, a.model4, a.no_net, a.opt_time, a.uniform)) as pool:
        res = dict(pool.imap_unordered(solve, [(i, n, c, a.cycles) for i, (c, _) in enumerate(data)]))
    orig = np.array([key(c)[0] for c, _ in data]); qk = np.array([q for _, q in data])
    ours = np.array([res[i]["cx"] for i in range(len(data))])
    eq = np.mean([res[i]["equiv"] for i in range(len(data))])
    out = {"uniform": a.uniform, "tag": a.tag, "no_net": a.no_net, "opt_time": a.opt_time, "model4": a.model4, "qubits": n, "count": len(data), "width": a.width, "widen": a.widen, "orig_cx": orig.mean(),
           "qiskit_opt3_cx": qk.mean(), "block_cx": ours.mean(), "equiv": eq,
           "better_than_qiskit": float(np.mean(ours < qk)), "equal": float(np.mean(ours == qk)),
           "worse": float(np.mean(ours > qk)), "minutes": (time.time() - t0) / 60,
           "sec_per_circuit": float(np.mean([res[i]["sec"] for i in range(len(data))]))}
    print("ИТОГ:", json.dumps(out), flush=True)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / f"block{a.tag}_{n}q.json").write_text(json.dumps(
        {"summary": out, "circuits": [{"orig": data[i][0], "qiskit_cx": int(data[i][1]), **res[i]} for i in range(len(data))]}))


if __name__ == "__main__":
    main()
