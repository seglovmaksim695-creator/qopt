"""Негативные тесты вердикта: пустой файл лемм, чужой файл лемм, старый результат без проверок потока -> не PROVED."""
import json, os, shutil, subprocess, sys
R = '/Users/max/quantum-start/verify_bridge3_20261123'
D, DF = R + '/D/results', R + '/DFIX'
out = []
def trial(name, lemma_obj=None, result_src=None):
    vt = DF + '/tmp/vt'
    shutil.rmtree(vt, ignore_errors=True)
    os.makedirs(vt + '/tools'); os.makedirs(vt + '/results')
    shutil.copy(DF + '/tools/qkal_verdict.py', vt + '/tools/')
    os.symlink(DF + '/src', vt + '/src')
    for f in ('lemmas_final.json', 'primes.json', 'canaries_n256_p25519.json', 'canaries_n256_p25519_T8b.json', 'canaries_n256_p25519_extra.json'):
        shutil.copy(DF + '/results/' + f, vt + '/results/')
    if lemma_obj is not None:
        json.dump(lemma_obj, open(vt + '/results/lemmas_final.json', 'w'))
    shutil.copy(result_src or DF + '/results/dfix_n64_p2305843009213693951.json', vt + '/results/dfix_n64_p2305843009213693951.json')
    r = subprocess.run([sys.executable, '-B', vt + '/tools/qkal_verdict.py', 'dfix'], capture_output=True, text=True)
    try:
        v = json.load(open(vt + '/results/verdict_dfix.json'))['rows'][0]['status']
    except Exception as e:
        v = 'вердикт упал: %s %s' % (e, r.stderr[-200:])
    out.append((name, v)); print(name, '->', v)
good = json.load(open(DF + '/results/lemmas_final.json'))
trial('контроль: настоящие файлы')
trial('пустой файл лемм', dict(ok=True, implications=0, proved=0, controls=0, controls_ok=0, rows=[], z3='x', source_sha256='0'))
trial('пустые строки, но счётчики как у настоящего', dict(good, rows=[]))
trial('чужой файл лемм (старая попытка D lemmas_try1)', json.load(open(D + '/lemmas_try1.json')))
g2 = json.loads(json.dumps(good)); g2['source_sha256'] = '0' * 64
trial('хеш исходника лемм другой', g2)
g3 = json.loads(json.dumps(good)); g3['rows'][5]['sha256'] = '0' * 64
trial('хеш одного запроса подменён', g3)
trial('старый результат D (без проверок потока), n=64', None, D + '/final2_n64_p2305843009213693951.json')
json.dump(out, open(DF + '/results/verdict_negtests.json', 'w'), ensure_ascii=False, indent=1)
