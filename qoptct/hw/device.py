"""Device model for ``qoptct hw``: coupling graph and calibration data of an IBM-like backend.

A :class:`Device` is built from a Qiskit ``Target`` (fake backends of ``qiskit_ibm_runtime.fake_provider`` carry a
snapshot of real calibration data; a real backend object works the same way, but nothing in this package ever
submits a job).  Synthetic devices (``line:N``, ``ring:N``, ``grid:RxC``, ``heavyhex:D``) are available for tests.
"""
import math

import numpy as np

BAD = 0.3  # a gate / readout error above this marks the qubit or coupler as unusable ("faulty" in the calibration)

_ALIASES = {"torino": "FakeTorino", "fez": "FakeFez", "marrakesh": "FakeMarrakesh", "sherbrooke": "FakeSherbrooke",
            "kyiv": "FakeKyiv", "brisbane": "FakeBrisbane", "kingston": "FakeKingston", "aachen": "FakeAachen",
            "osaka": "FakeOsaka", "kawasaki": "FakeKawasaki", "quebec": "FakeQuebec", "cusco": "FakeCusco"}


def get_backend(name):
    """Backend object for a name such as ``fake_torino`` / ``torino`` / ``FakeTorino`` or a synthetic
    ``line:8``, ``ring:8``, ``grid:3x4``, ``heavyhex:3`` (uniform errors).  A backend object is returned as is."""
    if not isinstance(name, str):
        return name
    low = name.lower()
    if ":" in low:
        return _synthetic(low)
    key = low.replace("fake_", "").replace("fake", "").replace("ibm_", "").replace("_", "")
    try:
        from qiskit_ibm_runtime import fake_provider as fp
    except ImportError as e:  # pragma: no cover
        raise ImportError("fake IBM backends need qiskit-ibm-runtime:  pip install 'qoptct[hw]'") from e
    cls = getattr(fp, _ALIASES.get(key, ""), None)
    if cls is None:
        for attr in dir(fp):
            if attr.lower() == "fake" + key or attr.lower() == "fake" + key + "v2":
                cls = getattr(fp, attr)
                break
    if cls is None:
        raise ValueError(f"unknown backend {name!r} (examples: fake_torino, fake_fez, fake_sherbrooke, line:8)")
    return cls()


def _synthetic(spec):
    from qiskit.providers.fake_provider import GenericBackendV2
    from qiskit.transpiler import CouplingMap
    kind, arg = spec.split(":", 1)
    if kind == "line":
        cm = CouplingMap.from_line(int(arg))
    elif kind == "ring":
        cm = CouplingMap.from_ring(int(arg))
    elif kind == "grid":
        r, c = arg.split("x")
        cm = CouplingMap.from_grid(int(r), int(c))
    elif kind == "heavyhex":
        cm = CouplingMap.from_heavy_hex(int(arg))
    elif kind == "full":
        cm = CouplingMap.from_full(int(arg))
    else:
        raise ValueError(f"unknown synthetic device {spec!r}")
    return GenericBackendV2(cm.size(), basis_gates=["cz", "rz", "sx", "x", "id"], coupling_map=cm, seed=11,
                            noise_info=True)


class Device:
    """Coupling graph + calibration snapshot.

    Attributes: ``n`` qubits, ``twoq`` (name of the native two-qubit gate), ``edges`` (set of unordered pairs
    (a, b), a < b, usable couplers only), ``err2[(a, b)]``, ``dur2[(a, b)]``, ``err1[q]``, ``dur1[q]`` (sx),
    ``ro_err[q]``, ``ro_dur[q]``, ``t1[q]``, ``t2[q]`` (seconds; inf when unknown), ``nbr[q]`` (neighbours),
    ``hop`` (hop distance matrix), ``wdist`` (noise-weighted distance in units of "average SWAPs").
    """

    def __init__(self, backend):
        self.backend = get_backend(backend)
        self.name = getattr(self.backend, "name", None) or type(self.backend).__name__
        t = self.target = self.backend.target
        n = self.n = t.num_qubits
        names = set(t.operation_names)
        self.twoq = next((g for g in ("cz", "ecr", "cx") if g in names), None)
        if self.twoq is None:
            raise ValueError("backend has none of the two-qubit gates cz / ecr / cx")

        def prop(gate, qargs, attr, default):
            try:
                p = t[gate].get(tuple(qargs))
            except KeyError:
                return default
            v = getattr(p, attr, None) if p is not None else None
            return default if v is None else float(v)

        self.err1 = [prop("sx", (q,), "error", 0.0) for q in range(n)]
        self.dur1 = [prop("sx", (q,), "duration", 0.0) for q in range(n)]
        self.ro_err = [prop("measure", (q,), "error", 0.0) for q in range(n)]
        self.ro_dur = [prop("measure", (q,), "duration", 0.0) for q in range(n)]
        qp = t.qubit_properties or [None] * n
        self.t1 = [float(getattr(p, "t1", None) or math.inf) if p is not None else math.inf for p in qp]
        self.t2 = [float(getattr(p, "t2", None) or math.inf) if p is not None else math.inf for p in qp]
        self.err2, self.dur2, self.direction = {}, {}, {}
        for qargs in t[self.twoq]:
            a, b = qargs
            e = prop(self.twoq, (a, b), "error", 0.0)
            d = prop(self.twoq, (a, b), "duration", 0.0)
            k = (min(a, b), max(a, b))
            if k not in self.err2 or e < self.err2[k]:
                self.err2[k], self.dur2[k], self.direction[k] = e, d, (a, b)
        self.bad_qubits = {q for q in range(n) if self.ro_err[q] >= BAD or self.err1[q] >= BAD}
        self.edges = {k for k, e in self.err2.items() if e < BAD and not (set(k) & self.bad_qubits)}
        self.nbr = [[] for _ in range(n)]
        for a, b in sorted(self.edges):
            self.nbr[a].append(b); self.nbr[b].append(a)
        self._dist()

    def _dist(self):
        import rustworkx as rx
        n = self.n
        g = rx.PyGraph()
        g.add_nodes_from(range(n))
        es = sorted(self.edges)
        logs = [-math.log(max(1e-9, 1.0 - self.err2[e])) for e in es]
        med = float(np.median(logs)) if logs else 1.0
        med = med if med > 0 else 1.0
        self.elog = {e: l for e, l in zip(es, logs)}  # -log fidelity of one native 2q gate on the coupler
        self.elog_med = med
        for e, l in zip(es, logs):
            # weight of one hop in units of an average SWAP; never below 0.5 so that hop count still matters
            g.add_edge(e[0], e[1], max(0.5, min(l / med, 8.0)))
        self.graph = g
        self.hop = np.asarray(rx.distance_matrix(g, null_value=1e6), dtype=float)
        w = rx.all_pairs_dijkstra_path_lengths(g, lambda x: float(x))
        wd = np.full((n, n), 1e6)
        for a in range(n):
            wd[a, a] = 0.0
            for b, d in w[a].items():
                wd[a, b] = d
        self.wdist = wd

    def has_edge(self, a, b):
        return ((a, b) if a < b else (b, a)) in self.edges

    def idle_logfid(self, q, t):
        """log of the average fidelity of an idle period t (seconds) on qubit q (T1/T2 damping)."""
        if t <= 0:
            return 0.0
        f = 0.5 + math.exp(-t / self.t1[q]) / 6.0 + math.exp(-t / self.t2[q]) / 3.0
        return math.log(max(f, 1e-12))
