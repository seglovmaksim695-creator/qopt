# SPDX-License-Identifier: Apache-2.0
"""Command line for `qoptct estimate` (also `python -m qoptct.estimate`).

    qoptct estimate --toffoli 4.91e9 --qubits 1400 [--depth D] [--model gidney2025|ge2019|litinski2023|all]
                    [--shots 9.2] [--p 1e-3] [--cycle-us 1] [--reaction-us 10] [--json]
    qoptct estimate --check        # reproduce the published numbers of the three papers
    qoptct estimate --attack rsa2048-g25|ecdlp-p256|ecdlp-k1|ecdlp-b163|...|aes128 [--variant L3|L2|L0]
                    [--source both|published|qoptct] [--json]     # logical model + physical estimate
    qoptct estimate --list-attacks | --check-attacks              # table published / qoptct blocks; reproduction

Exit code: 0 = estimate printed, 3 = usage or input error (one `error: ...` line on stderr)."""
import argparse
import json
import math
import sys

from . import logical as lg
from . import physical as ph


def build_parser(prog="qoptct estimate"):
    ap = argparse.ArgumentParser(prog=prog, description="Physical cost (surface code) of a fault-tolerant computation "
                                 "from its logical Toffoli/AND count and logical qubit count. Models: ge2019 "
                                 "(arXiv:1905.09749), gidney2025 (arXiv:2505.15917), litinski2023 (arXiv:2306.08585).")
    ap.add_argument("--toffoli", type=float, help="Toffoli + logical AND count (CCZ states)")
    ap.add_argument("--qubits", type=float, help="logical qubits")
    ap.add_argument("--depth", type=float, default=None, help="reaction depth (optional; see the model docstrings)")
    ap.add_argument("--model", default="all", choices=["all"] + sorted(ph.MODELS), help="model (default: all)")
    ap.add_argument("--shots", type=float, default=1.0, help="gidney2025: number of shots the Toffoli count is "
                    "spread over (Table 5 counts Toffolis per factoring = 9.2 shots for RSA-2048); default 1")
    ap.add_argument("--p", type=float, default=1e-3, help="physical error rate (default 1e-3)")
    ap.add_argument("--cycle-us", type=float, default=1.0, help="surface code cycle in us (default 1)")
    ap.add_argument("--reaction-us", type=float, default=10.0, help="reaction time in us (default 10)")
    ap.add_argument("--json", action="store_true", help="print JSON")
    ap.add_argument("--check", action="store_true", help="reproduce the papers' published numbers and exit")
    ap.add_argument("--attack", default=None, help="end-to-end logical model + physical estimate of a published "
                    "algorithm (qoptct.estimate.logical): %s; --toffoli/--qubits are then not needed" %
                    " | ".join(sorted(lg.ATTACKS)))
    ap.add_argument("--variant", default=None, help="--attack ecdlp-p256|ecdlp-k1: Karatsuba levels of the lab point "
                    "addition, L3 (default, most qubits), L2, L0 (fewest qubits)")
    ap.add_argument("--source", default="both", choices=["both", "published", "qoptct"],
                    help="--attack: published numbers, numbers with qoptct blocks, or both (default)")
    ap.add_argument("--list-attacks", action="store_true", help="table 'published / with qoptct blocks' for all attacks")
    ap.add_argument("--check-attacks", action="store_true",
                    help="reproduce the published tables of the logical models (Table 5 of Gidney 2025, ...) and exit")
    return ap


def _fmt_time(s):
    if s >= 2 * 86400:
        return "%.2f days" % (s / 86400)
    if s >= 3600:
        return "%.2f h" % (s / 3600)
    if s >= 120:
        return "%.1f min" % (s / 60)
    return "%.1f s" % s


def check_rows():
    """[(model, quantity, ours, published)] for the reproduction of the papers."""
    rows = []
    ref = ph.ge2019_rsa_reference()
    rows += [("ge2019 (exact layout)", "physical qubits", ref["physical_qubits"], 20e6),
             ("ge2019 (exact layout)", "hours per run", ref["hours"], 5.1),
             ("ge2019 (exact layout)", "expected hours", ref["expected_hours"], 7.4),
             ("ge2019 (exact layout)", "Toffoli", ref["toffoli"], 2.7e9)]
    g = ph.ge2019(ref["toffoli"], 6189, ref["reaction_depth"])
    rows += [("ge2019 (generic)", "physical qubits", g.physical_qubits, 20e6),
             ("ge2019 (generic)", "hours per run", g.runtime_s / 3600, 5.1),
             ("ge2019 (generic)", "expected hours", g.expected_runtime_s / 3600, 7.4)]
    g = ph.ge2019(ref["toffoli"], 6189, ref["reaction_depth"], code_distance=27, d1_options=(15,),
                  algorithm_error=ref["deviation_error"])
    rows += [("ge2019 (generic, GE's d)", "physical qubits", g.physical_qubits, 20e6),
             ("ge2019 (generic, GE's d)", "retry risk", g.failure_prob, 0.31),
             ("ge2019 (generic, GE's d)", "expected hours", g.expected_runtime_s / 3600, 7.4)]
    x = ph.gidney2025(6.5e9, 1399, shots=9.2)
    rows += [("gidney2025", "physical qubits", x.physical_qubits, 897864),
             ("gidney2025", "hours per shot", x.runtime_s / 3600, 12.07),
             ("gidney2025", "expected days", x.expected_runtime_s / 86400, 4.96)]
    l10 = ph.litinski2023(109e6, 3000, cycle_us=10)
    l1 = ph.litinski2023(109e6, 3000, cycle_us=1)
    b = ph.litinski2023(109e6, 3000, architecture="baseline")
    rows += [("litinski2023", "physical qubits", l10.physical_qubits, 6000 * 1152),
             ("litinski2023", "minutes per key (10 us cycle)", l10.runtime_s / 60, 9.7),
             ("litinski2023", "seconds per key (1 us cycle)", l1.runtime_s, 58.0),
             ("litinski2023 baseline", "physical qubits", b.physical_qubits, 9.4e6),
             ("litinski2023 baseline", "hours per key (1 us)", b.runtime_s / 3600, 3.8)]
    return rows


def main_estimate(argv=None):
    ap = build_parser()
    try:
        a = ap.parse_args(argv)
    except SystemExit as e:
        return 3 if e.code else 0
    if a.check:
        rows = check_rows()
        if a.json:
            print(json.dumps([dict(model=m, quantity=q, ours=o, published=p) for m, q, o, p in rows]))
        else:
            for m, q, o, p in rows:
                print("%-24s %-30s ours %-11.4g paper %-11.4g (%+.1f %%)" % (m, q, o, p, 100 * (o / p - 1)))
        return 0
    if a.check_attacks:
        rows = lg.reproduce_published()
        if a.json:
            print(json.dumps([dict(name=n, ours=g, published=w, tolerance=t, ok=ok, source=s)
                              for n, g, w, t, ok, s in rows]))
        else:
            for n, g, w, t, ok, s in rows:
                print("%-4s %-36s ours %-12.5g published %-12.5g tol %.1f %%  (%s)" % ("ok" if ok else "FAIL", n, g, w,
                                                                                    100 * t, s))
        return 0 if all(r[4] for r in rows) else 1
    if a.list_attacks:
        print(lg.attack_table())
        return 0
    if a.attack is not None:
        try:
            rep = lg.attack_report(a.attack, a.variant, a.source)
        except ValueError as e:
            print("error: %s" % e, file=sys.stderr)
            return 3
        print(json.dumps(rep, default=str) if a.json else lg.format_report(rep))
        return 0
    if a.toffoli is None or a.qubits is None:
        print("error: --toffoli and --qubits are required (or --check, --attack, --list-attacks)", file=sys.stderr)
        return 3
    if not math.isfinite(a.qubits) or a.qubits != int(a.qubits):
        print("error: --qubits must be an integer", file=sys.stderr)
        return 3
    models = sorted(ph.MODELS) if a.model == "all" else [a.model]
    out = []
    for m in models:
        kw = dict(p=a.p, cycle_us=a.cycle_us, reaction_us=a.reaction_us)
        if m == "gidney2025":
            kw["shots"] = a.shots
        try:
            out.append(ph.estimate(m, a.toffoli, int(a.qubits), a.depth, **kw))
        except ValueError as e:
            print("error: %s: %s" % (m, e), file=sys.stderr)
            return 3
    if a.json:
        print(json.dumps([r.to_dict() for r in out], default=str))
        return 0
    for r in out:
        print("%-13s %8.3g physical qubits   d=%-3d expected time %-11s (one run %s, P_fail %.1f %%)   "
              "%.3g megaqubit-days%s" % (r.model, r.physical_qubits, r.code_distance, _fmt_time(r.expected_runtime_s),
                                         _fmt_time(r.runtime_s), 100 * r.failure_prob, r.megaqubit_days,
                                         "   [reaction-limited]" if r.reaction_limited else ""))
        for n in r.notes:
            print("              note: " + n)
    return 0


if __name__ == "__main__":
    sys.exit(main_estimate())
