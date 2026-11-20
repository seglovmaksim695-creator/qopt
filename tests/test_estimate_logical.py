"""qoptct.estimate.logical: published tables are reproduced; the qoptct variants are pinned (regression); the CLI.

Tolerances: Gidney 2025 Table 5 is printed with 2 significant digits and the model must hit them EXACTLY in all 7 rows
(Toffoli), P_deviant, E(shots) and qubits exactly; GE19 Table 1 2 %; Litinski Fig. 9c 0.6 % ('109M' rounding); Haner
2020 Table 1 2 %; Roetteler 2017 their own regression 2.5 %, qubits exactly; Banegas 2020 Table 5 and Taguchi-Takayasu
Table 6 exactly; Jaques 2020 Table 8 exactly, Table 10 0.5 %.

The block counts (BLOCKS) are compared with live qoptct.crypto counters: Kaliski (current count 12n^2 - 8n - 3) with
the special-prime reduction, GF(2^m) multiplier and squarer, narrow lookup."""
import json
import random

import pytest

from qoptct.estimate import logical as lg
from qoptct.estimate.cli import main_estimate


# ---------------------------------------------------------------------------------------------- published numbers
def test_reproduce_published_all_rows():
    rows = lg.reproduce_published()
    assert len(rows) >= 75
    bad = [r for r in rows if not r[4]]
    assert not bad, bad


@pytest.mark.parametrize("n", sorted(lg.G25_ROWS))
def test_gidney2025_table5_row(n):
    row = lg.G25_ROWS[n]
    r = lg.gidney2025_rsa(n)
    assert float('%.2g' % r.toffoli) == row[9]                  # Table 5 prints 2 digits
    assert r.qubits == row[10]
    assert round(r.extra['shots'], 1) == row[8]
    assert round(r.extra['p_deviant'] * 100, 2) == round(row[7] * 100, 2)
    assert abs(lg.gidney2025_rsa_table3(n).toffoli / row[9] - 1) < 0.025


def test_gidney2025_2048_value():
    r = lg.gidney2025_rsa(2048)
    assert abs(r.toffoli / 6.497e9 - 1) < 1e-3 and r.qubits == 1399 and r.status == 'published'


def test_ge2019_tables():
    assert abs(lg.ge2019_rsa(2048).toffoli / 2.698e9 - 1) < 1e-3
    assert lg.ge2019_rsa(2048).qubits == 6189
    assert abs(lg.ge2019_rsa(1024, dev_off=5).toffoli / 0.41e9 - 1) < 0.02
    assert abs(lg.ge2019_rsa(3072, cmul=4, dev_off=6).toffoli / 9.9e9 - 1) < 0.02


def test_litinski_fig9c():
    assert lg.litinski2023_ecdlp(256).extra['pointadd'] == 126.5 * 256 ** 2 + 189 * 256
    for k in (1, 2, 4):
        assert abs(lg.litinski2023_ecdlp(256, k=k).toffoli / ((44 * k + 65) * 1e6 / k) - 1) < 0.006


def test_binary_and_aes_published():
    assert lg.taguchi_takayasu2023_ecdlp(163).toffoli == 13175432
    assert lg.banegas2020_ecdlp(163).toffoli == 293095880
    assert lg.jaques2020_aes(128, 1).extra['oracle'] == 13727
    assert abs(lg.jaques2020_aes(128, 1).toffoli * 4 / (1.32 * 2 ** 79) - 1) < 0.005


# ---------------------------------------------------------------------------------------------- qoptct variants
def test_rsa_g25_qoptct_2048():
    r = lg.rsa_g25_qoptct(2048)
    assert r.status == 'qoptct' and r.qubits == 1400
    assert abs(r.toffoli / 4.920e9 - 1) < 2e-4                   # 4.91e9 in the lab model (see the loop3 start lookup)
    assert 0.2425 < r.extra['gain'] < 0.2430                     # -24.3 % against the code accounting
    assert r.extra['narrow'] == (41, 33) and r.extra['acc_bits'] == 29
    assert 'no end-to-end circuit' in r.note
    # lab K_CLAIMS gains for the other sizes (-26.3, -22.4, -21.2 %), to 0.3 percentage points
    for n, g in ((1024, 0.263), (4096, 0.224), (8192, 0.212)):
        assert abs(lg.rsa_g25_qoptct(n).extra['gain'] - g) < 0.003, n
    assert 'extrapolation' in lg.rsa_g25_qoptct(1024).note
    # sensitivity: the typical (not best) narrow lookup for the 21-bit words costs 34 ANDs
    assert lg.rsa_g25_qoptct(2048, narrow_loop3=34).toffoli > r.toffoli


def test_ecdlp_prime_qoptct_numbers_and_change_against_lab():
    r = lg.ecdlp_prime_qoptct('p256', 'L3')
    assert r.extra['pointadd'] == 3041421 - 4 * 2307 == 3032193
    assert r.toffoli == 13 * (3032193 + 2 ** 16 + 768) == 40280461 and r.qubits == 6742
    old = lg.ecdlp_prime_qoptct('p256', 'L3', kaliski='legacy')
    assert old.toffoli == 40400425 and old.extra['pointadd'] == 3041421       # 4.04e7 of K_SUMMARY
    assert old.toffoli - r.toffoli == 13 * 4 * 2307
    l0 = lg.ecdlp_prime_qoptct('p256', 'L0')
    assert l0.toffoli == 13 * (3284291 + 65536 + 768) and l0.qubits == 4368     # 4.37e7 -> 4.36e7
    k1 = lg.ecdlp_prime_qoptct('secp256k1', 'L3')
    assert k1.extra['pointadd'] == 3033954 - 4 * 2307 and k1.qubits == 6934
    assert abs(r.extra['gain'] - 0.6314) < 5e-4
    assert lg.pointadd_ands('p256', 'L3', 'legacy') == 3041421
    assert lg.pointadd_ands_with_mulsq('p256', 'L3') < r.extra['pointadd']       # optional estimate, small


def test_ecdlp_gf2m_qoptct():
    want = {163: 6063080, 233: 14045148, 283: 23072728, 571: 107461640}
    for m, t in want.items():
        r = lg.ecdlp_gf2m_qoptct(m)
        assert r.toffoli == t and r.status == 'qoptct'
        assert 0.52 < r.extra['gain'] < 0.55
    assert lg.ecdlp_gf2m_qoptct(163).qubits == 1959
    assert lg.ecdlp_gf2m_qoptct(409).toffoli == 820 * (2 * 11 * 2501 + 817)
    w = lg.ecdlp_gf2m_windowed_estimate(163)
    assert w.status == 'estimate' and w.toffoli == 720766 and 'ESTIMATE' in w.note
    assert lg.ecdlp_gf2m_windowed_estimate(163, conservative=True).toffoli == 1061300


def test_aes_qoptct():
    r = lg.aes_oracle_qoptct(128, 1)
    assert r.extra['oracle'] == 10046 and r.qubits == 1696 and 'NOT independently verified' in r.note
    assert abs(r.toffoli / 1.4555e23 - 1) < 1e-3
    assert abs(r.extra['vs_jaques_with_meuli_sbox'] - 0.2229) < 5e-4
    assert lg.aes_oracle_qoptct(192, 2).status == 'estimate'


def test_withdrawn_variants_are_absent():
    for name in ('K1_CHAIN', 'K1_PAR', 'ge2019_rsa_qoptct', 'GE19_AND'):
        assert not hasattr(lg, name)
    import inspect
    assert 'blk' not in inspect.signature(lg.ge2019_rsa).parameters       # the GE-frame substitution is not offered
    assert all('ge19' not in k for k in lg.ATTACKS)


# ---------------------------------------------------------------------------------------------- blocks vs the library
def test_blocks_kaliski_current_count():
    from qoptct.crypto import kaliski_and_count
    n = 256
    assert kaliski_and_count(n) == 12 * n * n - 8 * n - 3
    assert kaliski_and_count(n, legacy=True) - kaliski_and_count(n) == 9 * n + 3 == 2307
    lib = lg.library_blocks(gf2m_ms=())
    assert lib['kaliski_special'] == lg.BLOCKS['kaliski_special'] == {'p256': 664966, 'secp256k1': 657503}
    assert lib['kaliski_special_legacy'] == lg.BLOCKS['kaliski_special_legacy']
    assert lg.BLOCKS['kaliski_special_legacy']['p256'] - lg.BLOCKS['kaliski_special']['p256'] == 2307


def test_blocks_gf2m_library():
    lib = lg.library_blocks(gf2m_ms=(163, 571))
    for m in (163, 571):
        assert lib['gf2m_mult'][m] == lg.BLOCKS['gf2m_mult'][m]
        assert lib['gf2m_square_relabel_cx'][m] == lg.BLOCKS['gf2m_square_relabel_cx'][m]
    # the point-addition result does not change when the library numbers are passed in
    assert lg.ecdlp_gf2m_qoptct(163, blocks=lib).toffoli == lg.ecdlp_gf2m_qoptct(163).toffoli


@pytest.mark.slow
def test_blocks_gf2m_library_all_fields():
    lib = lg.library_blocks()
    assert lib['gf2m_mult'] == lg.BLOCKS['gf2m_mult']
    assert lib['gf2m_square_relabel_cx'] == lg.BLOCKS['gf2m_square_relabel_cx'] == {163: 408, 233: 153, 283: 713,
                                                                                     409: 204, 571: 1538}


def test_blocks_narrow_lookup_constants():
    from qoptct.crypto import lookup
    rng = random.Random(10000)
    t32 = [rng.getrandbits(32) for _ in range(64)]
    c32 = lookup.narrow_lookup(t32, 6, 32)
    assert c32.n_and <= lg.BLOCKS['narrow_lookup']['loop1_c6_w32']          # 40 (41 in 2 of 40 tables)
    t21 = [rng.getrandbits(21) for _ in range(64)]
    c21 = lookup.narrow_lookup(t21, 6, 21)
    assert c21.n_and <= 35                                                  # 33 typical (best case used in the model)


def test_blocks_fold_mulsq_counts_used_in_notes():
    from qoptct.crypto import multiply, square
    from qoptct.crypto.fold import CURVES, fold_and_count
    assert multiply(256, 3).counts()['AND'] == 37826 and multiply(256, 3, legacy=True).counts()['AND'] == 37880
    assert square(256, 3).counts()['AND'] == 22145 and square(256, 3, legacy=True).counts()['AND'] == 22332
    assert fold_and_count(CURVES['p256']) == 2744 and fold_and_count(CURVES['secp256k1']) == 2175


# ---------------------------------------------------------------------------------------------- report and physical link
def test_attack_report_physical_link():
    rep = lg.attack_report('rsa2048-g25')
    pub, q = rep['rows']
    assert pub['physical']['model'] == 'gidney2025' and 4.9 < pub['physical']['expected_runtime_s'] / 86400 < 5.2
    assert 3.6 < q['physical']['expected_runtime_s'] / 86400 < 3.9            # K_SUMMARY: 5.02 -> 3.74 days
    assert q['physical']['physical_qubits'] < 0.95e6
    ec = lg.attack_report('ecdlp-p256', 'L0')
    assert ec['rows'][1]['qubits'] == 4368 and ec['rows'][1]['physical']['model'] == 'litinski2023'
    aes = lg.attack_report('aes128', physical=True)
    assert aes['rows'][0]['physical'] is None
    assert len(lg.attack_report('ecdlp-b163', source='qoptct')['rows']) == 1


def test_attack_report_errors():
    with pytest.raises(ValueError):
        lg.attack_report('nope')
    with pytest.raises(ValueError):
        lg.attack_report('rsa2048-g25', 'L3')
    with pytest.raises(ValueError):
        lg.attack_report('ecdlp-p256', 'L9')
    with pytest.raises(ValueError):
        lg.attack_report('aes128', source='x')


def test_attack_table_markdown():
    t = lg.attack_table()
    assert t.count('\n') == 2 + len(lg.ATTACKS) - 1 + 0 and 'rsa2048-g25' in t and 'aes128' in t


# ---------------------------------------------------------------------------------------------- CLI
def test_cli_attack_text(capsys):
    assert main_estimate(['--attack', 'rsa2048-g25']) == 0
    out = capsys.readouterr().out
    assert 'published' in out and 'qoptct' in out and 'physical (gidney2025)' in out and '-24.3 %' in out
    assert main_estimate(['--attack', 'ecdlp-p256', '--variant', 'L0', '--source', 'qoptct']) == 0
    out = capsys.readouterr().out
    assert '4368 logical qubits' in out and '\n  published' not in out


def test_cli_attack_json_and_list(capsys):
    assert main_estimate(['--attack', 'ecdlp-k1', '--json']) == 0
    d = json.loads(capsys.readouterr().out)
    assert d['name'] == 'ecdlp-k1' and len(d['rows']) == 2 and d['rows'][1]['qubits'] == 6934
    assert main_estimate(['--list-attacks']) == 0
    assert 'ecdlp-b571' in capsys.readouterr().out


def test_cli_check_attacks(capsys):
    assert main_estimate(['--check-attacks']) == 0
    out = capsys.readouterr().out
    assert 'FAIL' not in out and 'G25 RSA-3072 Toffoli' in out
    assert main_estimate(['--check-attacks', '--json']) == 0
    rows = json.loads(capsys.readouterr().out)
    assert all(r['ok'] for r in rows)


def test_cli_attack_errors(capsys):
    assert main_estimate(['--attack', 'nope']) == 3
    assert main_estimate(['--attack', 'aes128', '--variant', 'L3']) == 3
    assert main_estimate(['--attack', 'rsa2048-g25', '--source', 'bogus']) == 3
    assert main_estimate(['--qubits', '10']) == 3                              # old behaviour unchanged
    capsys.readouterr()
    assert main_estimate(['--toffoli', '1e9', '--qubits', '2000', '--model', 'gidney2025']) == 0
