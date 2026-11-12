# Benchmark circuits: provenance and licensing

The 23 `.qc` files in this directory are the standard reversible/Clifford+T benchmark suite distributed with the
**Feynman** toolkit by Matthew Amy: <https://github.com/meamy/feynman>, directory `benchmarks/qc`
(repository license: BSD-3-Clause; the Feynman README refers to M. Amy's PhD thesis *Formal Methods in Quantum
Circuit Design*, University of Waterloo, 2019, for the benchmark suite). The suite was popularised by
Amy, Maslov, Mosca, *Polynomial-time T-depth optimization of Clifford+T circuits via matroid partitioning*,
IEEE TCAD (2014), doi:10.1109/TCAD.2014.2341953, arXiv:1303.2042, and is used in most later T-count /
CNOT-count papers (e.g. Nam et al., npj Quantum Information (2018), doi:10.1038/s41534-018-0072-4,
arXiv:1710.07345). Several circuits originate from earlier reversible-logic benchmark collections — the original
authors hold the rights to them.

**Licensing caveat.** These files are third-party data and are *not* covered by the Apache-2.0 license of qopt.
They are redistributed here as they are in the Feynman repository. Before a public release the maintainer should
confirm that redistribution is acceptable (see RELEASE_NOTES_RU.md checklist).

**Not included:** `hwb6.qc` (and `hwb8`, `hwb10`, `gf2^16_mult` used in development) carry the notice
"Distributing this file without an express written permission or erasing this note is prohibited"
(Saeedi, Saheb Zamani, Sedighi, Sasanian, arXiv:1004.4320; Maslov). They are therefore not redistributed.
`scripts/fetch_benchmarks.sh` downloads `hwb6.qc` from the Feynman repository for local use only. Optimized
hwb6 circuits are likewise not included in `results/` in this version.

## Relation to the Feynman files

The files here were used, unchanged, to produce every number in `results/`. Compared with the current Feynman
`master` they differ only textually: inverse multi-qubit Z gates are written `Z` instead of `Zd` (identical
operators, since multi-controlled Z is self-inverse), and whitespace on the `BEGIN` line. This was checked by
comparing normalized gate lists (`Zd`→`Z`, blank lines and trailing spaces removed) for all 24 files on
2026-11-12: identical.

| file | qubits | SHA-256 |
|---|---|---|
| adder_8.qc | 24 | a1853896f8f51ed4829fc8d0d529fc727c4988dd19740644af40c8ba0763cb5a |
| barenco_tof_10.qc | 19 | e90adfdbe05ab723508d40fc7d988536fb20c6d5037f92eef0eaeb5cf1ce12e4 |
| barenco_tof_3.qc | 5 | a53b2eee4fe21333a436822471d6f038244cf4a1d6fdd242f13531e70153160e |
| barenco_tof_4.qc | 7 | ba665da53fa3d7456f546fbb7945457b0b7129204c1a91e7d2a8e0133619aa96 |
| barenco_tof_5.qc | 9 | 5bfd37ccb42361b70a41d62533b73bf83c7d5f57d40f8c9632f21b97546017ad |
| csla_mux_3.qc | 15 | 9d52e8a719c86b757b13c4bd0001c6e69be4afb16acfb773ce8b4428957d2d5d |
| csum_mux_9.qc | 30 | cce08cd1773fb8eee05c9e6fea124038d1d0605c120646fd58e94ce613112a7a |
| gf2^4_mult.qc | 12 | 1f0b3f6637d1d73d0c6e69030e77198ed89b9111aa0a5334458084f97e5fce33 |
| gf2^5_mult.qc | 15 | 6ae44de8a94f4998dfb44bf33330cd0106b9f508e63f205ab7b49e8905a21e0a |
| grover_5.qc | 9 | 757738a2052dd9c40d8667bec405323b4e321054e0dd4d8758be5d152fcf519a |
| ham15-low.qc | 17 | 4d1528b0b75f6db2ee9e1c83b943152f6a8d7fc8f2f5ba3d5bd2ba57a8fb59f1 |
| mod5_4.qc | 5 | 9d32b6df9264366676c3e259211d43c8686b7b3cbcb215e0bfa9d81c3919b3a6 |
| mod_mult_55.qc | 9 | 6740ed207932ae523291aaca0f0d38e5a0fbe903abd10871b2b1923fbbb4b92c |
| mod_red_21.qc | 11 | de0573cbb25fbca8b9546e485627901a9fd4135cc2ddfb59abc0e288f29f6656 |
| qcla_com_7.qc | 24 | bef451acf9cc9c08139c382845cdb583d0796a1e58dbb181efa1302ff0547940 |
| qcla_mod_7.qc | 26 | f001d6a6fbf472d4c3f2f168e45ea1c06f3d38623b0b8d5b397f4d47a3fd83c1 |
| qft_4.qc | 5 | 4637ef09983f9b54c3d4d51baaf3603f77493a24d9b236b5a9f3590a43f6b175 |
| rc_adder_6.qc | 14 | 337176687452af345b8c078bbe63fe91a8306451c603fe9e627dc4a9395df70b |
| tof_10.qc | 19 | 8958838674a8025346f3eed2219497c61039f012beef7ec22c4afcfd72980ae0 |
| tof_3.qc | 5 | 10caf87ef6ab65b8f798fa002fa4e2ebb5b44d1e3f2c95a21540121a9613ba7d |
| tof_4.qc | 7 | 5778913a5c34ef226c9496f0fb8bb772b119b7c3455b31a83492e16913fd85f4 |
| tof_5.qc | 9 | fde7a805c93176d6db2d2eb93d9aa1f81d58b6537c62ec5f4c37e27a4d25977f |
| vbe_adder_3.qc | 10 | d2a3157bacc9a13fbf5cc07fd491545b7bf13105454ef952abad5b7538621b8d |

The `.i` line of each file lists the input wires; the other wires are the intended ancillas (used by the
clean-ancilla model, see README).
