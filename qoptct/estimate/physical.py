# SPDX-License-Identifier: Apache-2.0
"""Simple, explicit physical cost models for fault-tolerant execution on a surface-code machine.

Input: logical counts of an algorithm -- Toffoli (= CCZ = logical AND) count T, logical qubit count Q and, optionally,
the reaction depth D (number of measure -> classical decision -> next operation steps on the critical path).
Output: physical qubits, run time per run and expected run time (including retries), space-time volume.

Every model is a short closed formula taken from one paper; the parameters and the place in the paper are listed
below and in the function docstrings. The models are calibrated to reproduce the published numbers of their paper
when they are given the paper's own logical counts (see PAPER and tests/test_estimate_physical.py). They are NOT a
layout or a compiler: applied to other circuits they assume that the circuit has the same structure (same mix of
arithmetic, same parallelism, same reaction/Toffoli ratio) as the circuit of the paper.

Common physical assumptions (all three papers): planar nearest-neighbour grid, uniform circuit noise p = 1e-3,
surface code cycle 1 us, control-system reaction time 10 us.

Models
------
ge2019        Gidney, Ekera, "How to factor 2048 bit RSA integers in 8 hours using 20 million noisy qubits",
              arXiv:1905.09749 (Quantum 5, 433, 2021).
              * logical error per d x d x 1 spacetime cell: 0.1 (100 p)^((d+1)/2)   (ancillary estimate_costs.py,
                topological_error_per_unit_cell; Sec. 2.13 states ~10^-(ceil(d/2)+1) at p = 1e-3, the same at d=27);
              * physical qubits per logical patch 2 (d+1)^2                         (Fig. 8; Sec. 2.14);
              * CCZ factory of Gidney-Fowler arXiv:1812.01238: 15-to-1 T factories at level-1 distance d1 feeding
                a CCZ factory at d2; footprint/depth formulas of estimate_costs.py (logical_factory_dimensions);
                at d1=15, d2=27 the factory is 13 x 7 patches and makes one CCZ every 5 d2 cycles (135 us);
              * distillation error 28 (35 e0^3 + topo)^2 + topo per CCZ        (compute_distillation_error);
              * run time: the adders are reaction limited (Sec. 2.11, Eq. 5): run time = D * reaction time, with the
                measurement depth of Eq. (3) 300 n_e n + 0.5 n_e n lg n; enough factories are placed to keep up;
              * choice of d1, d2: minimise (physical qubits)^1.2 * time / (1 - P_fail) with P_fail < 0.8
                (Sec. 3.1, rank_estimate);
              * data/routing area: GE's layout uses 9728 patches for 3n + 0.002 n lg n = 6189 logical qubits at
                n = 2048 (Sec. 2.14, Fig. 9: 3 register rows + adder + routing per piece) -> factor 1.572.
              ge2019_rsa_reference() re-implements GE's exact layout and parameter search for RSA (their published
              formulas) and is the reference against which the generic ge2019() is checked.
              Published (Table 3/Table 1, n = 2048, p = 1e-3): 20 Mqubits, 5.1 h per run, 31 % retry risk,
              expected 7.4 h (0.31 days), 2.7e9 Toffoli (Table 3 'B tofs'), 4.1 / 5.9 megaqubit-days.

gidney2025    Gidney, "How to factor 2048 bit RSA integers with less than a million noisy qubits",
              arXiv:2505.15917, Sec. 3.2, Fig. 6, Fig. 7, Table 5, App. A.
              * cold (idle) storage in yoked surface codes: 430 physical qubits per logical qubit at 1e-15 logical
                error per qubit-round (Fig. 6, 2D yoking);
              * hot storage: d = 25 patches, 2 (d+1)^2 = 1352 physical qubits each; 131 hot qubits (Sec. 3.2);
              * compute region 7 x 18 hot patches = 6 factories of 3 x 4 patches + 54 patches of workspace (Fig. 7);
              * magic state cultivation (T state at 1e-7, 30000 qubit-rounds) + 8T -> CCZ: 114.7 rounds per CCZ
                per factory, rounded up to 150; six factories -> one CCZ per 25 us = one lattice surgery layer
                at d = 25 (Sec. 3.2);
              * schedule: the paper books 2 ms for an addition that needs >= 32 CCZ x 25 us = 0.8 ms and 2 ms for a
                lookup that needs >= 1.575 ms (App. A.2, A.4) -> slack factor 2.5 over the CCZ supply rate;
                their own numbers give 12.07 h / (6.5e9 / 9.2) = 61.5 us per Toffoli, our 2.5 x 25 us = 62.5 us;
              * T in Table 5 is expected Toffolis per factoring (all E(shots) = 9.2 shots, caption of Table 5);
              * logical failure per shot: 1600 qubits x 12 h x 1e-15 -> 93.3 % success; expected time
                4.63 days / 0.933 = 4.96 days ("less than a week").
              Published: 897 864 physical qubits ("less than a million"), 4.96 days, 6.5e9 Toffoli, 1399 qubits.
              Simplifications: logical error per round at other distances scaled by Lambda = 10 per +2 of d
              (assumption, consistent with GE's formula at p = 1e-3); cold storage scaled as ((d+1)/26)^2 (assumption);
              the circuit is assumed to have the structure of Gidney's (few hot qubits, the rest idle in cold storage),
              which is NOT true for an arbitrary circuit.

litinski2023  Litinski, "How to compute a 256-bit elliptic curve private key with only 50 million Toffoli gates",
              arXiv:2306.08585, active-volume architecture of Litinski-Nickerson arXiv:2211.15465.
              * active volume ~ 60 blocks per Toffoli for arithmetic circuits (Sec. 3); per key (2.8+3.8/k) e9
                blocks, 109e6 Toffoli at k = 1 (Fig. 2) -> 6.6e9 blocks;
              * Q memory modules + Q workspace modules (2 Q in all), 2 d^2 physical qubits each (1152 at d = 24);
              * Q blocks are executed per logical cycle of d code cycles; time x 10/9 for the error budget (Sec. 2);
              * logical error 10^(-d/2) per block at p = 10 % of threshold (Sec. 2.1), i.e. (p/p_th)^(d/2);
              * reaction limit: depth <= 2 T, 36 min per key at 10 us (Sec. 2.3) -- applied only if D is given.
              Published: 6000 modules x 1152 qubits = 6.9e6 qubits, one key per 9.7 min at a 10 us code cycle
              (58 s at 1 us); 2D-local baseline (architecture='baseline'): 4 logical cycles per Toffoli, d = 28,
              2 Q patches of 2 d^2 -> 9.4e6 qubits, 3.8 h at 1 us (Sec. 2.1).
              Note: the 10 % of threshold assumption is optimistic for p = 1e-3 hardware (threshold ~1 %): it matches.

AND gates (Gidney logical AND with measurement-based uncomputation) cost one CCZ/Toffoli state each; the
uncomputation costs none. So 'toffoli' = number of Toffoli + number of logical AND computations.
"""
import math
from dataclasses import dataclass, field, asdict

US = 1e-6

# ---------------------------------------------------------------------------------------------------------------
# Published reference numbers (for the reproduction tests and the CLI --check)
PAPER = {
    "ge2019": dict(ref="arXiv:1905.09749, Table 3 / Table 1 (n=2048, p=1e-3)", toffoli=2.7e9, qubits=6189,
                   physical_qubits=20e6, runtime_h=5.1, expected_runtime_h=7.4, retry_risk=0.31,
                   megaqubit_days=4.1),
    "gidney2025": dict(ref="arXiv:2505.15917, Sec. 3.2 / Table 5 (n=2048)", toffoli=6.5e9, qubits=1399, shots=9.2,
                       physical_qubits=897864, expected_runtime_days=4.96, shot_hours=12.07, success=0.933),
    "litinski2023": dict(ref="arXiv:2306.08585, abstract / Sec. 2.2 (k=1)", toffoli=109e6, qubits=3000,
                         physical_qubits=6000 * 1152, runtime_min_10us=9.7, runtime_s_1us=58.0, code_distance=24),
    "litinski2023_baseline": dict(ref="arXiv:2306.08585, Sec. 2.1", toffoli=109e6, qubits=3000,
                                  physical_qubits=9.4e6, runtime_h_1us=3.8, code_distance=28),
}

# GE 2019, Eq. (3) and Eq. (2) with n_e = 1.5 n, n = 2048: measurement (reaction) depth per Toffoli of their circuit
_GE_N = 2048
_GE_DEPTH = 300 * 1.5 * _GE_N ** 2 + 0.5 * 1.5 * _GE_N ** 2 * math.log2(_GE_N)
_GE_TOF = 0.3 * _GE_N ** 3 + 0.0005 * _GE_N ** 3 * math.log2(_GE_N)
GE_DEPTH_PER_TOFFOLI = _GE_DEPTH / _GE_TOF          # ~0.733
GE_ROUTING_OVERHEAD = 9728 / 6189                     # data + routing patches per logical qubit, GE layout n=2048


@dataclass
class PhysicalEstimate:
    model: str
    toffoli: float
    logical_qubits: int
    reaction_depth: float            # the depth used (given or assumed)
    physical_qubits: float
    runtime_s: float                 # one run (one shot), or the whole factoring for gidney2025 without retries
    expected_runtime_s: float        # including retries after logical failure (and all shots)
    failure_prob: float              # probability that one run fails (logical + distillation + ...)
    code_distance: int
    reaction_limited: bool
    factories: int = 0
    params: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    @property
    def megaqubit_days(self):
        """Expected space-time volume in megaqubit-days."""
        return self.physical_qubits / 1e6 * self.expected_runtime_s / 86400

    def to_dict(self):
        d = asdict(self)
        d["megaqubit_days"] = self.megaqubit_days
        return d

    def summary(self):
        t = self.expected_runtime_s
        tt = "%.2f days" % (t / 86400) if t >= 2 * 86400 else ("%.2f h" % (t / 3600) if t >= 3600 else "%.1f s" % t)
        return ("%s: %.3g physical qubits, expected time %s (run %.3g s, P_fail %.1f %%), d=%d, %.3g megaqubit-days%s"
                % (self.model, self.physical_qubits, tt, self.runtime_s, 100 * self.failure_prob, self.code_distance,
                   self.megaqubit_days, ", reaction-limited" if self.reaction_limited else ""))


def _check_inputs(toffoli, qubits, reaction_depth, p, cycle_us, reaction_us):
    for name, v in (("toffoli", toffoli), ("qubits", qubits), ("p", p), ("cycle_us", cycle_us),
                    ("reaction_us", reaction_us)):
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0:
            raise ValueError("%s must be a positive finite number, got %r" % (name, v))
    if qubits < 1 or int(qubits) != qubits:
        raise ValueError("qubits must be a positive integer, got %r" % (qubits,))
    if reaction_depth is not None:
        if isinstance(reaction_depth, bool) or not isinstance(reaction_depth, (int, float)) or \
                not math.isfinite(reaction_depth) or reaction_depth <= 0:
            raise ValueError("reaction_depth must be a positive finite number or None, got %r" % (reaction_depth,))
    if p >= 0.01:
        raise ValueError("p = %g is at or above the surface code threshold (~1e-2); the models do not apply" % p)


def _union(*ps):
    t = 1.0
    for q in ps:
        if q >= 1:
            return 1.0
        t *= 1 - q
    return 1 - t


# ---------------------------------------------------------------------------------------------------------------
# GE 2019
def _ge_cell_error(d, p):
    """Logical error per d x d x 1 spacetime cell (GE19 estimate_costs.py)."""
    return 0.1 * (100 * p) ** ((d + 1) / 2)


def _ge_factory(d1, d2):
    """(width, height, depth) of the Gidney-Fowler CCZ factory in d2 patches / d2 cycles (GE19 estimate_costs.py)."""
    t1_height = 4 * d1 / d2
    t1_width = 8 * d1 / d2
    t1_depth = 5.75 * d1 / d2
    ccz_depth, ccz_height, ccz_width = 5, 6, 3
    storage_width = 2 * d1 / d2
    t1_factories = int(math.ceil((8 / ccz_depth) / (1 / t1_depth)))
    col = t1_height * math.ceil(t1_factories / 2)
    width = int(math.ceil(t1_width * 2 + ccz_width + storage_width))
    height = int(math.ceil(max(ccz_height, col)))
    return width, height, max(ccz_depth, t1_depth)


def _ge_ccz_error(d1, d2, p):
    """Error per distilled CCZ state (GE19 compute_distillation_error, after Gidney-Fowler 1812.01238)."""
    e0 = p + 100 * _ge_cell_error(d1 // 2, p)
    e1 = 35 * e0 ** 3 + 1100 * _ge_cell_error(d1, p)
    return 28 * e1 ** 2 + 1000 * _ge_cell_error(d2, p)


def ge2019(toffoli, qubits, reaction_depth=None, *, p=1e-3, cycle_us=1.0, reaction_us=10.0,
           routing_overhead=GE_ROUTING_OVERHEAD, max_error=0.8, post_process_error=0.01, algorithm_error=0.0,
           d1_options=(9, 11, 13, 15, 17, 19, 21, 23), d2_options=tuple(range(9, 81, 2)), code_distance=None):
    """Gidney-Ekera 2019 style estimate (see module docstring for the formulas and their places in the paper).

    reaction_depth: measurement (reaction) depth of one run in GE's sense, Eq. (3): the run time in units of the
    reaction time (GE's lookups run at code speed, not reaction limited, and are counted as their time equivalent);
    None -> GE_DEPTH_PER_TOFFOLI * toffoli (the ratio of GE's own circuit, Eq. 3 / Eq. 2). Run time =
    reaction_depth * reaction time; the number of CCZ factories is the smallest that delivers `toffoli` CCZ states
    in that time. Distances d1, d2 minimise GE's skewed expected volume (Sec. 3.1); like GE we accept a large retry
    risk (max_error 0.8) when that minimises the expected volume.
    algorithm_error: failure probability of the algorithm itself (e.g. GE's approximate-arithmetic 'deviation'
    error, 0.004 at n = 2048), combined with the topological, distillation and post-processing errors.
    GE's own optimum at n = 2048 (d1 = 15, d2 = 27) has 31 % risk = topological 0.18 + distillation 0.15 + 1 %
    post-processing; this generic search usually prefers a larger d1 (fewer distillation failures).
    code_distance: fix d2 instead of searching it."""
    _check_inputs(toffoli, qubits, reaction_depth, p, cycle_us, reaction_us)
    if not (0 <= algorithm_error < 1) or not (0 <= post_process_error < 1):
        raise ValueError("algorithm_error and post_process_error must be in [0, 1)")
    notes = []
    if reaction_depth is None:
        reaction_depth = GE_DEPTH_PER_TOFFOLI * toffoli
        notes.append("reaction depth not given: assumed %.3f x Toffoli count (ratio of GE19's own circuit)"
                     % GE_DEPTH_PER_TOFFOLI)
    tc, tr = cycle_us * US, reaction_us * US
    d2s = (code_distance,) if code_distance else d2_options
    best = None
    for d1 in d1_options:
        for d2 in d2s:
            if d1 >= d2:
                continue
            fw, fh, fdepth = _ge_factory(d1, d2)
            ccz_time = fdepth * d2 * tc
            t_react = reaction_depth * tr
            n_fact = max(1, int(math.ceil(toffoli * ccz_time / t_react)))
            runtime = max(t_react, toffoli * ccz_time / n_fact)
            data = int(math.ceil(qubits * routing_overhead))
            patches = data + n_fact * fw * fh
            phys = patches * 2 * (d2 + 1) ** 2
            topo = data * (runtime / tc) * _ge_cell_error(d2, p)
            dist = toffoli * _ge_ccz_error(d1, d2, p)
            err = _union(topo, dist, post_process_error, algorithm_error)
            if err >= max_error:
                continue
            score = phys ** 1.2 * runtime / (1 - err)
            if best is None or score < best[0]:
                best = (score, d1, d2, n_fact, phys, runtime, err, topo, dist, fw * fh)
    if best is None:
        raise ValueError("no code distance in the search range meets max_error=%g" % max_error)
    _, d1, d2, n_fact, phys, runtime, err, topo, dist, farea = best
    return PhysicalEstimate(
        model="ge2019", toffoli=toffoli, logical_qubits=int(qubits), reaction_depth=reaction_depth,
        physical_qubits=phys, runtime_s=runtime, expected_runtime_s=runtime / (1 - err), failure_prob=err,
        code_distance=d2, reaction_limited=True, factories=n_fact,
        params=dict(p=p, cycle_us=cycle_us, reaction_us=reaction_us, l1_distance=d1, routing_overhead=routing_overhead,
                    factory_patches=farea, topological_error=topo, distillation_error=dist,
                    post_process_error=post_process_error, algorithm_error=algorithm_error, max_error=max_error),
        notes=notes)


def ge2019_rsa_reference(n=2048, p=1e-3, n_e=None):
    """Re-implementation of GE19's exact RSA layout and parameter search (formulas of the paper's ancillary file
    estimate_costs.py: eh_rsa, estimate_algorithm_cost, board_logical_dimensions). Returns a dict with the chosen
    parameters, Toffoli count, logical patches, physical qubits, hours and error. Used to check ge2019()."""
    if n_e is None:  # eh_rsa: Ekera-Hastad with s = 1, delta = 20
        m = math.ceil(n / 2) - 1
        n_e = m + 2 * (m - 20)
    tc, tr = 1 * US, 10 * US
    best = None
    for d1 in (5, 7, 9, 11, 13, 15, 17, 19, 21, 23):
        for d2 in range(9, 51, 2):
            for g_exp in (4, 5, 6):
                for g_mul in (4, 5, 6):
                    if g_mul > g_exp:
                        continue
                    for sep in (512, 768, 1024, 1536, 2048):
                        if n % sep:
                            continue
                        for dev_off in range(2, 10):
                            for tt in ((False, True) if (d1 == 15 and d2 >= 31) else (False,)):
                                r = _ge_rsa_cost(n, n_e, p, d1, d2, g_exp, g_mul, sep, dev_off, tt, tc, tr)
                                if r is None:
                                    continue
                                score = r["megaqubits"] ** 1.2 * r["hours"] / (1 - r["error"])
                                if best is None or score < best[0]:
                                    best = (score, r)
    return best[1]


def _ge_rsa_cost(n, n_e, p, d1, d2, g_exp, g_mul, sep, dev_off, tt, tc, tr):
    pad = int(math.ceil(math.log2(n * n * n_e))) + dev_off
    pieces = int(math.ceil(n / sep))
    piece_len = sep + pad
    mul_in = n + pad + 2
    loops = int(math.ceil(n_e * 2 * mul_in / (g_exp * g_mul)))
    dev_err = 4 * math.sqrt(loops * pieces / 2 ** pad)
    adder_depth = piece_len * 2 - 1
    lookup_depth = 2 ** (g_exp + g_mul) - 1
    unlookup_depth = 2 * math.sqrt(lookup_depth)
    if tt:
        fw, fh, fdepth = 24, 16, 6
    else:
        fw, fh, fdepth = _ge_factory(d1, d2)
    ccz_time = fdepth * tc * d2
    pairs = int(math.ceil(ccz_time / tr / 2))
    width = (fw + 1) * pairs + 1
    reg_h = int(math.ceil(piece_len / (width - 2)))
    height = fh * 2 + 3 * 2 + 3 + 6 + reg_h * 3
    piece_dist = fh * fw * pairs * 2
    logical = width * height * pieces
    dist_area = piece_dist * pieces
    tof = (adder_depth * pieces + lookup_depth + unlookup_depth) * loops
    loop_time = adder_depth * tr + lookup_depth * d2 * tc / 2 + unlookup_depth * d2 * tc / 2
    total = loop_time * loops
    topo = (logical - dist_area) * (total / tc) * _ge_cell_error(d2, p)
    dist = tof * (4 * 9 * 1e-17 if tt else _ge_ccz_error(d1, d2, p))
    err = _union(topo, dist, dev_err, 1e-2)
    if err >= 0.8:
        return None
    phys = logical * 2 * (d2 + 1) ** 2
    return dict(n=n, n_e=n_e, d1=d1, d2=d2, g_exp=g_exp, g_mul=g_mul, runway_sep=sep, dev_off=dev_off,
                tt_distillation=tt, toffoli=tof, logical_patches=logical, distillation_patches=dist_area,
                factories=pairs * 2 * pieces, physical_qubits=phys, megaqubits=phys / 1e6, hours=total / 3600,
                error=err, expected_hours=total / 3600 / (1 - err), adder_reaction_depth=adder_depth * loops,
                deviation_error=dev_err, topological_error=topo, distillation_error=dist,
                # run time in units of the reaction time (adders reaction limited + lookups at code speed): this is
                # the 'reaction depth' to give ge2019(); it is close to GE's measurement depth, Eq. (3)
                reaction_depth=total / tr)


# ---------------------------------------------------------------------------------------------------------------
# Gidney 2025
def gidney2025(toffoli, qubits, reaction_depth=None, *, shots=1.0, p=1e-3, cycle_us=1.0, reaction_us=10.0,
               hot_qubits=131, workspace_patches=54, factories=6, factory_patches=12, rounds_per_ccz=150,
               cold_qubits_per_logical=430, schedule_slack=2.5, error_per_round_d25=1e-15, lam=10.0,
               max_error=0.1, code_distance=None):
    """Gidney 2025 style estimate (yoked cold storage + hot d-patches + cultivation CCZ factories); see the module
    docstring for the formulas and their places in the paper.

    toffoli: Toffolis of the whole computation (all shots, as in Table 5); shots: expected number of shots it is
    split into (9.2 for RSA-2048 in Table 5; the logical failure probability is per shot). reaction_depth (same
    units as toffoli): only a lower bound on the time, depth * reaction time; the paper's time is set by CCZ supply.
    code_distance: hot patch distance; None -> the smallest odd d >= 15 with per-shot logical failure <= max_error."""
    _check_inputs(toffoli, qubits, reaction_depth, p, cycle_us, reaction_us)
    if not (isinstance(shots, (int, float)) and math.isfinite(shots) and shots >= 1):
        raise ValueError("shots must be a number >= 1, got %r" % (shots,))
    notes = []
    if abs(p - 1e-3) > 1e-12:
        notes.append("Gidney 2025 numbers (430 q/logical, cultivation 150 rounds) are for p = 1e-3; for p = %g only "
                     "the d-scaling via Lambda = %g is changed, which is a rough extrapolation" % (p, lam))
    lam_eff = lam * (1e-3 / p) if abs(p - 1e-3) > 1e-12 else lam
    tc, tr = cycle_us * US, reaction_us * US
    hot = min(hot_qubits, qubits)
    cold = qubits - hot
    compute = workspace_patches + factories * factory_patches
    ds = (code_distance,) if code_distance else tuple(range(15, 61, 2))
    chosen = None
    for d in ds:
        eps = error_per_round_d25 * lam_eff ** (-(d - 25) / 2)
        cold_pp = cold_qubits_per_logical * ((d + 1) / 26) ** 2
        hot_pp = 2 * (d + 1) ** 2
        phys = cold * cold_pp + (hot + compute) * hot_pp
        t_ccz = max(rounds_per_ccz * tc / factories, d * tc)
        t_tof = schedule_slack * t_ccz
        t_shot = toffoli / shots * t_tof
        rl = False
        if reaction_depth is not None and reaction_depth / shots * tr > t_shot:
            t_shot, rl = reaction_depth / shots * tr, True
        fail = 1 - math.exp(-(cold + hot + compute) * (t_shot / tc) * eps)
        chosen = (d, phys, t_shot, fail, rl, t_tof)
        if fail <= max_error:
            break
    d, phys, t_shot, fail, rl, t_tof = chosen
    if fail > max_error:
        notes.append("no distance up to %d meets max_error=%g" % (d, max_error))
    if reaction_depth is None:
        notes.append("reaction depth not given: time set by CCZ supply only (as in the paper)")
    return PhysicalEstimate(
        model="gidney2025", toffoli=toffoli, logical_qubits=int(qubits),
        reaction_depth=reaction_depth if reaction_depth is not None else float("nan"),
        physical_qubits=phys, runtime_s=t_shot, expected_runtime_s=t_shot * shots / (1 - fail), failure_prob=fail,
        code_distance=d, reaction_limited=rl, factories=factories,
        params=dict(p=p, cycle_us=cycle_us, reaction_us=reaction_us, shots=shots, hot_qubits=hot, cold_qubits=cold,
                    compute_patches=compute, time_per_toffoli_us=t_tof / US, schedule_slack=schedule_slack,
                    cold_qubits_per_logical=cold_qubits_per_logical * ((d + 1) / 26) ** 2, max_error=max_error),
        notes=notes + ["runtime_s is one shot; expected_runtime_s = all shots / (1 - P_fail)",
                       "assumes Gidney's structure: %d hot qubits, the other %d idle in yoked cold storage" % (hot, cold)])


# ---------------------------------------------------------------------------------------------------------------
# Litinski 2023
def litinski2023(toffoli, qubits, reaction_depth=None, *, architecture="active", p=1e-3, p_threshold=1e-2,
                 cycle_us=1.0, reaction_us=10.0, blocks_per_toffoli=60.0, cycles_per_toffoli=4.0,
                 error_margin=10 / 9, max_error=0.1, code_distance=None):
    """Litinski 2023 estimate; see the module docstring.

    architecture='active': active-volume architecture (2 Q modules, Q blocks per logical cycle, 60 blocks/Toffoli);
    architecture='baseline': 2D-local surface code (2 Q patches, 4 logical cycles per Toffoli), Sec. 2.1.
    qubits = Q memory qubits (Litinski: 3000 for one key; the device has 2 Q modules/patches).
    reaction_depth: if given, time >= depth * reaction time (Sec. 2.3: depth <= 2 T)."""
    _check_inputs(toffoli, qubits, reaction_depth, p, cycle_us, reaction_us)
    if architecture not in ("active", "baseline"):
        raise ValueError("architecture must be 'active' or 'baseline'")
    tc, tr = cycle_us * US, reaction_us * US
    ratio = p / p_threshold
    if architecture == "active":
        cycles = toffoli * blocks_per_toffoli / qubits
    else:
        cycles = toffoli * cycles_per_toffoli
    blocks_at_risk = 2 * qubits * cycles
    ds = (code_distance,) if code_distance else tuple(range(8, 101, 2))
    for d in ds:
        fail = min(1.0, blocks_at_risk * ratio ** (d / 2))
        if fail <= max_error:
            break
    notes = []
    if fail > max_error:
        notes.append("no distance up to %d meets max_error=%g" % (d, max_error))
    phys = 2 * qubits * 2 * d ** 2
    runtime = cycles * d * tc * error_margin
    rl = False
    if reaction_depth is not None and reaction_depth * tr > runtime:
        runtime, rl = reaction_depth * tr, True
    if reaction_depth is None:
        notes.append("reaction depth not given: not applied; Litinski's bound depth <= 2 T would give >= %.3g s"
                     % (2 * toffoli * tr))
    notes.append("expected_runtime_s includes Litinski's 10/9 margin only (no separate retry factor)")
    return PhysicalEstimate(
        model="litinski2023" + ("" if architecture == "active" else "_baseline"), toffoli=toffoli,
        logical_qubits=int(qubits), reaction_depth=reaction_depth if reaction_depth is not None else float("nan"),
        physical_qubits=phys, runtime_s=runtime, expected_runtime_s=runtime, failure_prob=fail, code_distance=d,
        reaction_limited=rl, factories=0,
        params=dict(p=p, p_threshold=p_threshold, cycle_us=cycle_us, reaction_us=reaction_us,
                    architecture=architecture, logical_cycles=cycles,
                    active_volume_blocks=toffoli * blocks_per_toffoli if architecture == "active" else None,
                    max_error=max_error),
        notes=notes)


MODELS = {"ge2019": ge2019, "gidney2025": gidney2025, "litinski2023": litinski2023}


def estimate(model, toffoli, qubits, reaction_depth=None, **kw):
    """estimate('ge2019' | 'gidney2025' | 'litinski2023', toffoli, qubits, reaction_depth=None, **model_params)."""
    if model not in MODELS:
        raise ValueError("unknown model %r (known: %s)" % (model, ", ".join(MODELS)))
    return MODELS[model](toffoli, qubits, reaction_depth, **kw)
