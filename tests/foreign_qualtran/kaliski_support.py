"""Shared helpers of the Kaliski tests (not a test module)."""
import json
import pathlib

REPO = pathlib.Path(__file__).resolve().parents[2]
RESULTS = REPO / 'results' / 'qualtran'
REFERENCE = pathlib.Path(__file__).resolve().parent / 'data' / 'kaliski_dfix_reference.json'

P25519 = 2 ** 255 - 19
P256 = 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1

# (n, p): odd primes p < 2^n - 1, the small pairs of the fast tests
SMALL_PAIRS = [(4, 13), (4, 7), (5, 29), (6, 61), (7, 113), (8, 251), (10, 1021), (12, 4093), (16, 32749)]


def reference():
    return json.loads(REFERENCE.read_text())['pairs']


def lemma_result():
    """The Z3 lemma result for the CURRENT lemma source: the committed file if it is valid, else computed and cached."""
    from qoptct.foreign.qualtran.kaliski import lemmas_z3, verdict
    f = RESULTS / 'kaliski_lemmas.json'
    if f.exists():
        lem = json.loads(f.read_text())
        if verdict.lemmas_ok(lem):
            return lem
    return lemmas_z3.ensure_lemmas()


def canary_result():
    f = RESULTS / 'kaliski_canaries_n256_p25519.json'
    return json.loads(f.read_text()) if f.exists() else None
