# SPDX-License-Identifier: Apache-2.0
"""qoptct loop run | verify | table | status

    qoptct loop run --state DIR --target SPEC [--target SPEC ...] [--hours H] [--proposers lookup,optimize]
                    [--canaries 20] [--budget SEC] [--max-steps N] [--prove-all]
    qoptct loop verify --state DIR [--no-reprove] [--budget SEC]      exit 0 = journal holds, 1 = problems found
    qoptct loop table --state DIR                                    the best-known table (from the journal)
    qoptct loop status --state DIR                                   counters of the proposals so far

Target specs: crypto:<name>, lookup:c6w16s1, lookupadd:c5w12s1, qc:<path>[:clean]. Long runs: start through
lab/tools/one (memory and priority limits); stop with Ctrl-C or kill, the next `run` continues from the journal.
"""
import argparse
import sys

from .engine import run, stats
from .journal import State, table_from, write_table
from .verify import verify


def main(argv=None):
    ap = argparse.ArgumentParser(prog='qoptct loop', description='propose -> prove -> accumulate')
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run')
    r.add_argument('--state', required=True)
    r.add_argument('--target', action='append', required=True, help='target spec (repeatable)')
    r.add_argument('--hours', type=float, default=1.0)
    r.add_argument('--proposers', default='lookup,optimize')
    r.add_argument('--canaries', type=int, default=20)
    r.add_argument('--budget', type=float, default=120.0, help='seconds per proof')
    r.add_argument('--canary-sec', type=float, default=600.0)
    r.add_argument('--job-timeout', type=float, default=1800.0)
    r.add_argument('--max-steps', type=int, default=None)
    r.add_argument('--prove-all', action='store_true', help='prove candidates that do not improve the front too (statistics)')
    v = sp.add_parser('verify')
    v.add_argument('--state', required=True)
    v.add_argument('--no-reprove', action='store_true')
    v.add_argument('--budget', type=float, default=120.0)
    for n in ('table', 'status'):
        t = sp.add_parser(n)
        t.add_argument('--state', required=True)
    try:
        a = ap.parse_args(argv)
    except SystemExit as e:
        return 3 if e.code else 0
    if a.cmd == 'run':
        try:
            s = run(a.state, a.target, names=tuple(x for x in a.proposers.split(',') if x), hours=a.hours,
                    canaries=a.canaries, budget=a.budget, canary_sec=a.canary_sec, job_timeout=a.job_timeout,
                    prove_all=a.prove_all, max_steps=a.max_steps, log=lambda m: print(m, flush=True))
        except (RuntimeError, ValueError) as e:
            print('error: %s' % e, file=sys.stderr)
            return 3
        print('attempts %(attempts)d, outcomes %(outcomes)s, journal records %(journal_records)d, '
              'improvements %(improvements)d' % s)
        return 0
    st = State(a.state)
    if a.cmd == 'verify':
        ok, probs, n = verify(a.state, reprove=not a.no_reprove, budget=a.budget,
                              log=(lambda m: print(m, flush=True)) if not a.no_reprove else None)
        if ok:
            print('OK: %d records verified' % n)
            return 0
        print('PROBLEMS (%d):' % len(probs))
        for p in probs:
            print('  - ' + p)
        return 1
    if a.cmd == 'table':
        tab = write_table(st) if st.read()[1] == [] else table_from(st.read()[0])
        for t in sorted(tab):
            b = tab[t]['best']
            print('%-40s AND %4d  T %5d  CX %6d  qubits %4d   (baseline AND %d CX %d, %d records)'
                  % (t, b['AND'], b['T'], b['CX'], b['qubits'], tab[t]['baseline']['AND'], tab[t]['baseline']['CX'],
                     tab[t]['records']))
        return 0
    s = stats(st)
    print('attempts %(attempts)d, outcomes %(outcomes)s, journal records %(journal_records)d, improvements %(improvements)d' % s)
    return 0


if __name__ == '__main__':
    sys.exit(main())
