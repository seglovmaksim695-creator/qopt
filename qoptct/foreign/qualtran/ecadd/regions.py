"""Map: top-level blocks of a step -> ranges of translated operations (for the targeted canaries)."""
import collections, re, sys


def regions(tr, depth=1):
    out = collections.OrderedDict()
    for t in tr['trace']:
        parts = t['leaf']['path'].split('/')
        key = '/'.join(parts[1:1 + depth]) if len(parts) > 1 else parts[0]
        lo, hi = t['logical_start'], t['logical_end']
        if key in out:
            out[key][1] = hi
            out[key][2] += hi - lo
        else:
            out[key] = [lo, hi, hi - lo]
    return out


if __name__ == '__main__':
    from qoptct.foreign.qualtran.ecadd.core import lower_step
    from qoptct.foreign.qualtran.ecadd.run_step import step_bloq
    k, n, p = int(sys.argv[1]), int(sys.argv[2]), eval(sys.argv[3])
    low, tr = lower_step(step_bloq(k, n, p))
    for key, (lo, hi, cnt) in regions(tr, int(sys.argv[4]) if len(sys.argv) > 4 else 1).items():
        print(f'{key:45s} ops[{lo}:{hi}] {cnt}')
