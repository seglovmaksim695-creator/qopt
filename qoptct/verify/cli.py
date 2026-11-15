# SPDX-License-Identifier: Apache-2.0
"""Command line for `qoptct prove` (the integrator registers main_prove as a sub-command).

    qoptct prove CIRCUIT.qasm REFERENCE.qc [--clean] [--budget SEC] [--no-qcec] [--quick] [--json]

Exit code: 0 = PROVED, 1 = REFUTED, 2 = UNKNOWN (not proved), 3 = usage or input error (missing file, malformed
circuit; one `error: ...` line on stderr).
"""
import argparse
import json
import sys

EXIT = {'PROVED': 0, 'REFUTED': 1, 'UNKNOWN': 2}


def build_parser(prog='qoptct prove'):
    ap = argparse.ArgumentParser(prog=prog, description='Exact verification of a qoptct OpenQASM 3 circuit against its '
                                 '.qc reference, for all inputs and all measurement outcomes.')
    ap.add_argument('qasm', help='circuit (OpenQASM 3 subset written by qoptct)')
    ap.add_argument('qc', help='reference circuit (.qc)')
    ap.add_argument('--clean', action='store_true',
                    help='clean-ancilla model: .qc wires not on the .i line start in |0> (and must end as the reference)')
    ap.add_argument('--budget', type=float, default=120.0, help='time budget in seconds (default 120; 0 = no limit)')
    ap.add_argument('--no-qcec', action='store_true', help='do not fall back to MQT QCEC for unitary circuits')
    ap.add_argument('--quick', action='store_true', help='also run the fast sampled check (not a proof)')
    ap.add_argument('--json', action='store_true', help='print the full result as JSON')
    return ap


def main_prove(argv=None):
    from . import prove, quick_check
    ap = build_parser()
    try:
        a = ap.parse_args(argv)
    except SystemExit as e:
        return 3 if e.code else 0
    r = prove(a.qasm, a.qc, clean=a.clean, budget_sec=a.budget or None, use_qcec=not a.no_qcec)
    det = str(r.get('detail', ''))
    if r['status'] == 'UNKNOWN' and (r.get('method') == 'parse' or 'unsupported gate' in det):
        # input error (missing file, malformed circuit): exit code 3, as for a usage error -- never 1 (REFUTED)
        if a.json:
            print(json.dumps(r, default=str))
        print('error: ' + det, file=sys.stderr)
        return 3
    if a.quick:
        q = quick_check(a.qasm, a.qc, clean=a.clean)
        r['quick'] = {k: q.get(k) for k in ('ok', 'F', 'engine', 'error', 'sampled') if k in q}
    if a.json:
        print(json.dumps(r, default=str))
    else:
        print('%s  (%s, %.2f s)' % (r['status'], r['method'], r.get('time', 0.0)))
        if r.get('detail'):
            print('  ' + r['detail'])
        if 'quick' in r:
            print('  quick check: %s' % r['quick'])
    return EXIT.get(r['status'], 2)


if __name__ == '__main__':
    sys.exit(main_prove())
