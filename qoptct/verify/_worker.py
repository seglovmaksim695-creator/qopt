# SPDX-License-Identifier: Apache-2.0
"""Worker process for qoptct.verify._limits.run_limited: reads (module, function, args, kwargs) as a pickle from
stdin, runs it and writes ('ok', result) or ('err', message, traceback) as a pickle to stdout."""
import importlib
import os
import pickle
import sys
import traceback


def main():
    out = os.fdopen(os.dup(sys.stdout.fileno()), 'wb')
    sys.stdout = sys.stderr            # anything printed by the function must not corrupt the result stream
    try:
        mod, name, args, kwargs = pickle.loads(sys.stdin.buffer.read())
        fn = importlib.import_module(mod)
        for part in name.split('.'):
            fn = getattr(fn, part)
        res = ('ok', fn(*args, **kwargs))
    except BaseException as e:
        res = ('err', '%s: %s' % (type(e).__name__, e), traceback.format_exc(limit=5))
    try:
        data = pickle.dumps(res)
    except Exception as e:
        data = pickle.dumps(('err', 'result not picklable: %s' % e, ''))
    out.write(data)
    out.flush()


if __name__ == '__main__':
    main()
