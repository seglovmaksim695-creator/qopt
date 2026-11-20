import sys, random, time
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__file__), '..', 'tests'))
import test_verify_algebra as T
from qoptct.verify import algebra as A
cases = [('clean', 3, 0), ('clean', 4, 1), ('mul', 3, 0), ('mul', 4, 0), ('mul', 5, 1), ('mul', 6, 2), ('mul', 7, 3), ('sq', 4, 0), ('sq', 5, 0), ('sq', 6, 1), ('sq', 8, 2), ('sq', 9, 1), ('fold', 5, 3), ('fold', 6, 5), ('fold', 7, 1), ('modmul', 5, 3), ('modmul', 6, 5)]
count = int(sys.argv[1])
tot = dict(PROVED=0, REFUTED=0, UNKNOWN=0, ERROR=0, wrong=0, proved_correct=0, FALSE_PROVED=0, FALSE_REFUTED=0)
for case in cases:
    if case[0] == 'clean':
        N, ops, spec = T._mul(case[1], case[2], 'mul', clean=True)
    elif case[0] in ('mul', 'sq'):
        N, ops, spec = T._mul(case[1], case[2], case[0])
    elif case[0] == 'fold':
        N, ops, spec = T._fold(case[1], case[2])
    else:
        N, ops, spec = T._modmul(case[1], case[2])
    rng = random.Random(12345 + len(ops))
    st = dict(PROVED=0, REFUTED=0, UNKNOWN=0, ERROR=0, wrong=0, proved_correct=0, FALSE_PROVED=0, FALSE_REFUTED=0)
    t0 = time.time()
    for new in T.mutants(N, ops, rng, count):
        truth = T.exhaustive_ok(N, new, spec)
        r = A.prove_ops(N, new, spec, budget_sec=30)
        st[r['status']] += 1
        if not truth:
            st['wrong'] += 1
            if r['status'] == 'PROVED':
                st['FALSE_PROVED'] += 1
                print('FALSE PROVED', case, r['detail'], flush=True)
        elif r['status'] == 'PROVED':
            st['proved_correct'] += 1
        elif r['status'] == 'REFUTED':
            st['FALSE_REFUTED'] += 1
            print('FALSE REFUTED', case, r['detail'], flush=True)
        if r['status'] == 'ERROR':
            print('ERROR', case, r['detail'], flush=True)
    for k in tot:
        tot[k] += st[k]
    print(case, 'ops', len(ops), st, '%.0fs' % (time.time() - t0), flush=True)
print('TOTAL', tot)
