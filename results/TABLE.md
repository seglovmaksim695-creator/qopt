# Results

All circuits listed here are files in `results/` and were re-verified for this release (`scripts/verify_all.py`).
T = number of t/tdg gates, CX = number of cx+cz gates (conditional ones included), (+A) = extra qubits beyond the
wires of the input `.qc`. Input counts use the standard 7-T, 6-CX decomposition of each Toffoli/CCZ.

* **unitary** — `qopt --mode cx` (best of the development runs): a plain Clifford+T circuit, best = minimum CX.
* **strict** — measurement-based; every input wire is data in an arbitrary state; fresh ancillas start in |0>
  and must end in a basis state. Best = minimum T, then CX, then extra qubits.
* **clean** — measurement-based; wires not in the `.i` line start in |0> and must be returned to |0>.
  Best = minimum T, then CX, then extra qubits. (Every strict circuit is also valid here.)

| circuit | qubits | Toffolis | input T/CX | unitary T/CX | strict T/CX (+A) | clean T/CX (+A) |
|---|---|---|---|---|---|---|
| adder_8 | 24 | 57 | 399/409 | 193/223 | — | — |
| barenco_tof_10 | 19 | 32 | 224/192 | 72/72 | 59/48 (+3) | 59/48 (+3) |
| barenco_tof_3 | 5 | 4 | 28/24 | 16/14 | 8/9 (+2) | 8/9 (+2) |
| barenco_tof_4 | 7 | 8 | 56/48 | 24/24 | 12/13 (+3) | 12/13 (+3) |
| barenco_tof_5 | 9 | 12 | 84/72 | 32/32 | 16/17 (+4) | 16/17 (+4) |
| csla_mux_3 | 15 | 10 | 70/80 | 64/64 | 40/70 (+1) | 28/84 (+3) |
| csum_mux_9 | 30 | 28 | 196/168 | 84/112 | — | — |
| gf2^4_mult | 12 | 16 | 112/99 | 68/86 | 58/115 (+0) | 36/79 (+1) |
| gf2^5_mult | 15 | 25 | 175/154 | 115/132 | 90/183 (+0) | 56/124 (+1) |
| grover_5 | 9 | 48 | 336/288 | 166/160 | 124/155 (+4) | 124/155 (+4) |
| ham15-low | 17 | 23 | 161/236 | 97/179 | — | 68/177 (+1) |
| mod5_4 | 5 | 4 | 28/28 | 7/10 | 4/11 (+1) | 4/11 (+0) |
| mod_mult_55 | 9 | 7 | 49/48 | 35/39 | 20/35 (+2) | 12/24 (+1) |
| mod_red_21 | 11 | 17 | 119/105 | 73/61 | 44/58 (+5) | 28/48 (+0) |
| qcla_com_7 | 24 | 29 | 203/186 | 95/111 | 59/565 (+37) | 59/565 (+37) |
| qcla_mod_7 | 26 | 59 | 413/382 | 237/256 | — | — |
| qft_4 | 5 | 2 | 69/46 | 67/37 | 53/265 (+42) | 53/265 (+42) |
| rc_adder_6 | 14 | 11 | 77/93 | 47/65 | 24/54 (+6) | 24/53 (+6) |
| tof_10 | 19 | 17 | 119/102 | 71/54 | 59/71 (+3) | 36/37 (+1) |
| tof_3 | 5 | 3 | 21/18 | 15/12 | 8/10 (+2) | 8/9 (+1) |
| tof_4 | 7 | 5 | 35/30 | 23/18 | 12/15 (+3) | 12/13 (+1) |
| tof_5 | 9 | 7 | 49/42 | 31/24 | 16/20 (+4) | 16/17 (+1) |
| vbe_adder_3 | 10 | 10 | 70/70 | 24/28 | 12/33 (+3) | 12/33 (+3) |

## Pareto fronts (T, CX, extra qubits) of the measurement-based models

| model | circuit | T | CX | +A | qubits | checker | source | method |
|---|---|---|---|---|---|---|---|---|
| strict | barenco_tof_10 | 59 | 48 | 3 | 22 | verify_qasm strict | `agents/c12/out/barenco_tof_10_T59_CX48_A3.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_10 | 63 | 50 | 2 | 21 | verify_qasm strict | `agents/c12/out/barenco_tof_10_T63_CX50_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_10 | 67 | 52 | 1 | 20 | verify_qasm strict | `agents/c12/out/barenco_tof_10_T67_CX52_A1.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_10 | 72 | 74 | 0 | 19 | verify_qasm strict | `agents/c12/out/barenco_tof_10_T72_CX74_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_3 | 8 | 9 | 2 | 7 | verify_qasm strict | `agents/c12/out/barenco_tof_3_T8_CX9_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_3 | 11 | 9 | 1 | 6 | verify_qasm strict | `agents/c1/best_c2/barenco_tof_3_T11_CX9_A1.qasm` | AND compiler + CX resynthesis |
| strict | barenco_tof_3 | 16 | 14 | 0 | 5 | verify_qasm strict | `agents/c8/qasm/barenco_tof_3_T16_CX14_e0.qasm` | partial H gadgets + TODD + CX polish |
| strict | barenco_tof_4 | 12 | 13 | 3 | 10 | verify_qasm strict | `agents/c12/out/barenco_tof_4_T12_CX13_A3.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_4 | 15 | 14 | 2 | 9 | verify_qasm strict | `agents/c12/out/barenco_tof_4_T15_CX14_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_4 | 19 | 16 | 1 | 8 | verify_qasm strict | `agents/c12/out/barenco_tof_4_T19_CX16_A1.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_4 | 24 | 24 | 0 | 7 | verify_qasm strict | `agents/c8/qasm/barenco_tof_4_T24_CX24_e0.qasm` | partial H gadgets + TODD + CX polish |
| strict | barenco_tof_5 | 16 | 17 | 4 | 13 | verify_qasm strict | `agents/c12/out/barenco_tof_5_T16_CX17_A4.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_5 | 19 | 18 | 3 | 12 | verify_qasm strict | `agents/c12/out/barenco_tof_5_T19_CX18_A3.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_5 | 23 | 20 | 2 | 11 | verify_qasm strict | `agents/c12/out/barenco_tof_5_T23_CX20_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_5 | 27 | 22 | 1 | 10 | verify_qasm strict | `agents/c12/out/barenco_tof_5_T27_CX22_A1.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | barenco_tof_5 | 32 | 34 | 0 | 9 | verify_qasm strict | `agents/c12/out/barenco_tof_5_T32_CX34_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | csla_mux_3 | 40 | 70 | 1 | 16 | verify_qasm strict | `agents/c2/out/csla_mux_3_T40_CX70_A1.qasm` | AND compiler |
| strict | csla_mux_3 | 56 | 80 | 0 | 15 | verify_qasm strict | `agents/c2/out/csla_mux_3_T56_CX80_A0.qasm` | AND compiler |
| strict | csla_mux_3 | 60 | 75 | 0 | 15 | verify_qasm strict | `agents/c8/hq/csla_mux_3_T60_CX75_e0.qasm` | hybrid AND compiler (genuine H kept) |
| strict | csla_mux_3 | 64 | 70 | 0 | 15 | verify_qasm strict | `agents/c2/out/csla_mux_3_T64_CX70_A0.qasm` | AND compiler |
| strict | gf2^4_mult | 58 | 115 | 0 | 12 | verify_qasm strict | `agents/c12/out/gf2^4_mult_T58_CX115_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | gf2^4_mult | 64 | 83 | 1 | 13 | verify_qasm strict | `agents/c2/out/gf2^4_mult_T64_CX83_A1.qasm` | AND compiler |
| strict | gf2^4_mult | 68 | 99 | 0 | 12 | verify_qasm strict | `agents/c2/out/gf2^4_mult_T68_CX99_A0.qasm` | AND compiler |
| strict | gf2^5_mult | 90 | 183 | 0 | 15 | verify_qasm strict | `agents/c12/out/gf2^5_mult_T90_CX183_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | gf2^5_mult | 100 | 129 | 1 | 16 | verify_qasm strict | `agents/c2/out/gf2^5_mult_T100_CX129_A1.qasm` | AND compiler |
| strict | gf2^5_mult | 115 | 154 | 0 | 15 | verify_qasm strict | `agents/c2/out/gf2^5_mult_T115_CX154_A0.qasm` | AND compiler |
| strict | grover_5 | 124 | 155 | 4 | 13 | verify_qasm strict | `agents/c8/hq/grover_5_T124_CX155_e4.qasm` | hybrid AND compiler (genuine H kept) |
| strict | grover_5 | 133 | 158 | 3 | 12 | verify_qasm strict | `agents/c8/hq/grover_5_T133_CX158_e3.qasm` | hybrid AND compiler (genuine H kept) |
| strict | grover_5 | 157 | 177 | 2 | 11 | verify_qasm strict | `agents/c8/hq/grover_5_T157_CX177_e2.qasm` | hybrid AND compiler (genuine H kept) |
| strict | grover_5 | 184 | 244 | 0 | 9 | verify_qasm strict | `agents/c8/hq/grover_5_T184_CX244_e0.qasm` | hybrid AND compiler (genuine H kept) |
| strict | grover_5 | 185 | 208 | 1 | 10 | verify_qasm strict | `agents/c8/hq/grover_5_T185_CX208_e1.qasm` | hybrid AND compiler (genuine H kept) |
| strict | mod5_4 | 4 | 11 | 1 | 6 | verify_qasm strict | `agents/c1/best_c2/mod5_4_T4_CX11_A1.qasm` | AND compiler + CX resynthesis |
| strict | mod5_4 | 8 | 12 | 0 | 5 | verify_qasm strict | `agents/c8/qasm/mod5_4_T8_CX12_e0.qasm` | partial H gadgets + TODD + CX polish |
| strict | mod_mult_55 | 20 | 35 | 2 | 11 | verify_qasm strict | `agents/c2/out/mod_mult_55_T20_CX35_A2.qasm` | AND compiler |
| strict | mod_mult_55 | 24 | 39 | 1 | 10 | verify_qasm strict | `agents/c2/out/mod_mult_55_T24_CX39_A1.qasm` | AND compiler |
| strict | mod_mult_55 | 28 | 38 | 1 | 10 | verify_qasm strict | `agents/c2/out/mod_mult_55_T28_CX38_A1.qasm` | AND compiler |
| strict | mod_mult_55 | 33 | 55 | 0 | 9 | verify_qasm strict | `agents/c2/out/mod_mult_55_T33_CX55_A0.qasm` | AND compiler |
| strict | mod_mult_55 | 35 | 40 | 0 | 9 | verify_qasm strict | `agents/c2/out/mod_mult_55_T35_CX40_A0.qasm` | AND compiler |
| strict | mod_red_21 | 44 | 58 | 5 | 16 | verify_qasm strict | `agents/c2/out/mod_red_21_T44_CX58_A5.qasm` | AND compiler |
| strict | mod_red_21 | 47 | 58 | 4 | 15 | verify_qasm strict | `agents/c1/best_c2/mod_red_21_T47_CX58_A4.qasm` | AND compiler + CX resynthesis |
| strict | mod_red_21 | 51 | 63 | 3 | 14 | verify_qasm strict | `agents/c1/best_c2/mod_red_21_T51_CX63_A3.qasm` | AND compiler + CX resynthesis |
| strict | mod_red_21 | 58 | 68 | 2 | 13 | verify_qasm strict | `agents/c1/best_c2/mod_red_21_T58_CX68_A2.qasm` | AND compiler + CX resynthesis |
| strict | mod_red_21 | 69 | 79 | 1 | 12 | verify_qasm strict | `agents/c1/best_c2/mod_red_21_T69_CX79_A1.qasm` | AND compiler + CX resynthesis |
| strict | qcla_com_7 | 59 | 565 | 37 | 61 | verify_big (sum-over-paths) | `agents/c10/qcla_com_7_final.qasm` | H gadgets + TODD + tg3 CX polish |
| strict | qft_4 | 53 | 265 | 42 | 47 | verify_big (sum-over-paths) | `agents/c10/qft_4_final.qasm` | H gadgets + TODD + tg3 CX polish |
| strict | qft_4 | 63 | 44 | 1 | 6 | verify_qasm strict | `agents/c8/hq/qft_4_T63_CX44_e1.qasm` | hybrid AND compiler (genuine H kept) |
| strict | qft_4 | 67 | 46 | 0 | 5 | verify_qasm strict | `agents/c8/hq/qft_4_T67_CX46_e0.qasm` | hybrid AND compiler (genuine H kept) |
| strict | rc_adder_6 | 24 | 54 | 6 | 20 | verify_qasm strict | `agents/c8/hq/rc_adder_6_T24_CX54_e6.qasm` | hybrid AND compiler (genuine H kept) |
| strict | rc_adder_6 | 31 | 60 | 4 | 18 | verify_qasm strict | `agents/c8/hq/rc_adder_6_T31_CX60_e4.qasm` | hybrid AND compiler (genuine H kept) |
| strict | rc_adder_6 | 35 | 64 | 3 | 17 | verify_qasm strict | `agents/c2/out/rc_adder_6_T35_CX64_A3.qasm` | AND compiler |
| strict | rc_adder_6 | 39 | 67 | 2 | 16 | verify_qasm strict | `agents/c2/out/rc_adder_6_T39_CX67_A2.qasm` | AND compiler |
| strict | rc_adder_6 | 43 | 70 | 1 | 15 | verify_qasm strict | `agents/c2/out/rc_adder_6_T43_CX70_A1.qasm` | AND compiler |
| strict | rc_adder_6 | 47 | 73 | 0 | 14 | verify_qasm strict | `agents/c2/out/rc_adder_6_T47_CX73_A0.qasm` | AND compiler |
| strict | tof_10 | 59 | 71 | 3 | 22 | verify_qasm strict | `agents/c2/out/tof_10_T59_CX71_A3.qasm` | AND compiler |
| strict | tof_10 | 63 | 76 | 2 | 21 | verify_qasm strict | `agents/c2/out/tof_10_T63_CX76_A2.qasm` | AND compiler |
| strict | tof_10 | 67 | 81 | 1 | 20 | verify_qasm strict | `agents/c2/out/tof_10_T67_CX81_A1.qasm` | AND compiler |
| strict | tof_10 | 71 | 86 | 0 | 19 | verify_qasm strict | `agents/c2/out/tof_10_T71_CX86_A0.qasm` | AND compiler |
| strict | tof_3 | 8 | 10 | 2 | 7 | verify_qasm strict | `agents/c2/out/tof_3_T8_CX10_A2.qasm` | AND compiler |
| strict | tof_3 | 11 | 10 | 1 | 6 | verify_qasm strict | `agents/c1/best_c2/tof_3_T11_CX10_A1.qasm` | AND compiler + CX resynthesis |
| strict | tof_3 | 15 | 12 | 0 | 5 | verify_qasm strict | `agents/c8/qasm/tof_3_T15_CX12_e0.qasm` | partial H gadgets + TODD + CX polish |
| strict | tof_4 | 12 | 15 | 3 | 10 | verify_qasm strict | `agents/c2/out/tof_4_T12_CX15_A3.qasm` | AND compiler |
| strict | tof_4 | 15 | 15 | 2 | 9 | verify_qasm strict | `agents/c1/best_c2/tof_4_T15_CX15_A2.qasm` | AND compiler + CX resynthesis |
| strict | tof_4 | 19 | 20 | 1 | 8 | verify_qasm strict | `agents/c1/best_c2/tof_4_T19_CX20_A1.qasm` | AND compiler + CX resynthesis |
| strict | tof_4 | 59 | 54 | 0 | 7 | verify_qasm strict | `agents/c12/out/tof_4_T59_CX54_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| strict | tof_5 | 16 | 20 | 4 | 13 | verify_qasm strict | `agents/c2/out/tof_5_T16_CX20_A4.qasm` | AND compiler |
| strict | tof_5 | 19 | 20 | 3 | 12 | verify_qasm strict | `agents/c1/best_c2/tof_5_T19_CX20_A3.qasm` | AND compiler + CX resynthesis |
| strict | tof_5 | 23 | 25 | 2 | 11 | verify_qasm strict | `agents/c1/best_c2/tof_5_T23_CX25_A2.qasm` | AND compiler + CX resynthesis |
| strict | tof_5 | 27 | 28 | 1 | 10 | verify_qasm strict | `agents/c1/best_c2/tof_5_T27_CX28_A1.qasm` | AND compiler + CX resynthesis |
| strict | tof_5 | 31 | 24 | 0 | 9 | verify_qasm strict | `agents/c8/qasm/tof_5_T31_CX24_e0.qasm` | partial H gadgets + TODD + CX polish |
| strict | vbe_adder_3 | 12 | 33 | 3 | 13 | verify_qasm strict | `agents/c1/best_c2/vbe_adder_3_T12_CX33_A3.qasm` | AND compiler + CX resynthesis |
| strict | vbe_adder_3 | 19 | 35 | 2 | 12 | verify_qasm strict | `agents/c1/best_c2/vbe_adder_3_T19_CX35_A2.qasm` | AND compiler + CX resynthesis |
| strict | vbe_adder_3 | 23 | 36 | 0 | 10 | verify_qasm strict | `agents/c8/qasm/vbe_adder_3_T23_CX36_e0.qasm` | partial H gadgets + TODD + CX polish |
| strict | vbe_adder_3 | 24 | 32 | 0 | 10 | verify_qasm strict | `agents/c8/qasm/vbe_adder_3_T24_CX32_e0.qasm` | partial H gadgets + TODD + CX polish |
| clean | barenco_tof_10 | 59 | 48 | 3 | 22 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_10_T59_CX48_A3.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_10 | 63 | 50 | 2 | 21 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_10_T63_CX50_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_10 | 67 | 52 | 1 | 20 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_10_T67_CX52_A1.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_10 | 72 | 74 | 0 | 19 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_10_T72_CX74_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_3 | 8 | 9 | 2 | 7 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_3_T8_CX9_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_3 | 11 | 9 | 1 | 6 | strict verdict (strict implies clean) | `agents/c1/best_c2/barenco_tof_3_T11_CX9_A1.qasm` | AND compiler + CX resynthesis |
| clean | barenco_tof_3 | 16 | 14 | 0 | 5 | strict verdict (strict implies clean) | `agents/c8/qasm/barenco_tof_3_T16_CX14_e0.qasm` | partial H gadgets + TODD + CX polish |
| clean | barenco_tof_4 | 12 | 13 | 3 | 10 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_4_T12_CX13_A3.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_4 | 15 | 14 | 2 | 9 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_4_T15_CX14_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_4 | 19 | 16 | 1 | 8 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_4_T19_CX16_A1.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_4 | 24 | 24 | 0 | 7 | strict verdict (strict implies clean) | `agents/c8/qasm/barenco_tof_4_T24_CX24_e0.qasm` | partial H gadgets + TODD + CX polish |
| clean | barenco_tof_5 | 16 | 17 | 4 | 13 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_5_T16_CX17_A4.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_5 | 19 | 18 | 3 | 12 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_5_T19_CX18_A3.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_5 | 23 | 20 | 2 | 11 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_5_T23_CX20_A2.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_5 | 27 | 22 | 1 | 10 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_5_T27_CX22_A1.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | barenco_tof_5 | 32 | 34 | 0 | 9 | strict verdict (strict implies clean) | `agents/c12/out/barenco_tof_5_T32_CX34_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | csla_mux_3 | 28 | 84 | 3 | 18 | verify_qasm clean | `agents/c1/best_clean/csla_mux_3.qasm` | hand-derived construction (clean ancillas) |
| clean | csla_mux_3 | 40 | 60 | 1 | 16 | verify_qasm clean | `agents/c2/out_clean/csla_mux_3_T40_CX60_A1.qasm` | AND compiler (clean ancillas) |
| clean | csla_mux_3 | 52 | 66 | 0 | 15 | verify_qasm clean | `agents/c2/out_clean/csla_mux_3_T52_CX66_A0.qasm` | AND compiler (clean ancillas) |
| clean | csla_mux_3 | 53 | 58 | 0 | 15 | verify_qasm clean | `agents/c2/out_clean/csla_mux_3_T53_CX58_A0.qasm` | AND compiler (clean ancillas) |
| clean | gf2^4_mult | 36 | 79 | 1 | 13 | verify_qasm clean | `agents/c1/best_clean/gf2^4_mult.qasm` | hand-derived construction (clean ancillas) |
| clean | gf2^4_mult | 58 | 115 | 0 | 12 | strict verdict (strict implies clean) | `agents/c12/out/gf2^4_mult_T58_CX115_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | gf2^4_mult | 64 | 75 | 1 | 13 | verify_qasm clean | `agents/c2/out_clean/gf2^4_mult_T64_CX75_A1.qasm` | AND compiler (clean ancillas) |
| clean | gf2^4_mult | 68 | 99 | 0 | 12 | strict verdict (strict implies clean) | `agents/c2/out/gf2^4_mult_T68_CX99_A0.qasm` | AND compiler |
| clean | gf2^4_mult | 69 | 90 | 0 | 12 | verify_qasm clean | `agents/c2/out_clean/gf2^4_mult_T69_CX90_A0.qasm` | AND compiler (clean ancillas) |
| clean | gf2^4_mult | 70 | 87 | 0 | 12 | verify_qasm clean | `agents/c2/out_clean/gf2^4_mult_T70_CX87_A0.qasm` | AND compiler (clean ancillas) |
| clean | gf2^5_mult | 56 | 124 | 1 | 16 | verify_qasm clean | `agents/c1/best_clean/gf2^5_mult.qasm` | hand-derived construction (clean ancillas) |
| clean | gf2^5_mult | 90 | 183 | 0 | 15 | strict verdict (strict implies clean) | `agents/c12/out/gf2^5_mult_T90_CX183_A0.qasm` | MCX recognition (AND ladder / KG24 rewrite) |
| clean | gf2^5_mult | 100 | 119 | 1 | 16 | verify_qasm clean | `agents/c2/out_clean/gf2^5_mult_T100_CX119_A1.qasm` | AND compiler (clean ancillas) |
| clean | gf2^5_mult | 108 | 139 | 0 | 15 | verify_qasm clean | `agents/c2/out_clean/gf2^5_mult_T108_CX139_A0.qasm` | AND compiler (clean ancillas) |
| clean | grover_5 | 124 | 155 | 4 | 13 | strict verdict (strict implies clean) | `agents/c8/hq/grover_5_T124_CX155_e4.qasm` | hybrid AND compiler (genuine H kept) |
| clean | grover_5 | 133 | 158 | 3 | 12 | strict verdict (strict implies clean) | `agents/c8/hq/grover_5_T133_CX158_e3.qasm` | hybrid AND compiler (genuine H kept) |
| clean | grover_5 | 157 | 177 | 2 | 11 | strict verdict (strict implies clean) | `agents/c8/hq/grover_5_T157_CX177_e2.qasm` | hybrid AND compiler (genuine H kept) |
| clean | grover_5 | 184 | 244 | 0 | 9 | strict verdict (strict implies clean) | `agents/c8/hq/grover_5_T184_CX244_e0.qasm` | hybrid AND compiler (genuine H kept) |
| clean | grover_5 | 185 | 208 | 1 | 10 | strict verdict (strict implies clean) | `agents/c8/hq/grover_5_T185_CX208_e1.qasm` | hybrid AND compiler (genuine H kept) |
| clean | ham15-low | 68 | 177 | 1 | 18 | verify_qasm clean | `agents/c2/out_clean/ham15-low_T68_CX177_A1.qasm` | AND compiler (clean ancillas) |
| clean | ham15-low | 89 | 186 | 0 | 17 | verify_qasm clean | `agents/c2/out_clean/ham15-low_T89_CX186_A0.qasm` | AND compiler (clean ancillas) |
| clean | mod5_4 | 4 | 11 | 0 | 5 | verify_qasm clean | `agents/c2/out_clean/mod5_4_T4_CX11_A0.qasm` | AND compiler (clean ancillas) |
| clean | mod_mult_55 | 12 | 24 | 1 | 10 | verify_qasm clean | `agents/c2/out_clean/mod_mult_55_T12_CX24_A1.qasm` | AND compiler (clean ancillas) |
| clean | mod_mult_55 | 15 | 24 | 0 | 9 | verify_qasm clean | `agents/c2/out_clean/mod_mult_55_T15_CX24_A0.qasm` | AND compiler (clean ancillas) |
| clean | mod_red_21 | 28 | 48 | 0 | 11 | verify_qasm clean | `agents/c1/best_clean/mod_red_21.qasm` | hand-derived construction (clean ancillas) |
| clean | mod_red_21 | 40 | 46 | 1 | 12 | verify_qasm clean | `agents/c2/out_clean/mod_red_21_T40_CX46_A1.qasm` | AND compiler (clean ancillas) |
| clean | qcla_com_7 | 59 | 565 | 37 | 61 | strict verdict (strict implies clean) | `agents/c10/qcla_com_7_final.qasm` | H gadgets + TODD + tg3 CX polish |
| clean | qft_4 | 53 | 265 | 42 | 47 | strict verdict (strict implies clean) | `agents/c10/qft_4_final.qasm` | H gadgets + TODD + tg3 CX polish |
| clean | qft_4 | 63 | 44 | 1 | 6 | strict verdict (strict implies clean) | `agents/c8/hq/qft_4_T63_CX44_e1.qasm` | hybrid AND compiler (genuine H kept) |
| clean | qft_4 | 67 | 46 | 0 | 5 | strict verdict (strict implies clean) | `agents/c8/hq/qft_4_T67_CX46_e0.qasm` | hybrid AND compiler (genuine H kept) |
| clean | rc_adder_6 | 24 | 53 | 6 | 20 | verify_qasm clean | `agents/c2/out_clean/rc_adder_6_T24_CX53_A6.qasm` | AND compiler (clean ancillas) |
| clean | rc_adder_6 | 24 | 69 | 5 | 19 | verify_qasm clean | `agents/c1/best_clean/rc_adder_6.qasm` | hand-derived construction (clean ancillas) |
| clean | rc_adder_6 | 28 | 62 | 5 | 19 | verify_qasm clean | `agents/c2/out_clean/rc_adder_6_T28_CX62_A5.qasm` | AND compiler (clean ancillas) |
| clean | rc_adder_6 | 31 | 59 | 3 | 17 | verify_qasm clean | `agents/c2/out_clean/rc_adder_6_T31_CX59_A3.qasm` | AND compiler (clean ancillas) |
| clean | rc_adder_6 | 35 | 65 | 2 | 16 | verify_qasm clean | `agents/c2/out_clean/rc_adder_6_T35_CX65_A2.qasm` | AND compiler (clean ancillas) |
| clean | rc_adder_6 | 39 | 72 | 1 | 15 | verify_qasm clean | `agents/c2/out_clean/rc_adder_6_T39_CX72_A1.qasm` | AND compiler (clean ancillas) |
| clean | rc_adder_6 | 43 | 68 | 0 | 14 | verify_qasm clean | `agents/c2/out_clean/rc_adder_6_T43_CX68_A0.qasm` | AND compiler (clean ancillas) |
| clean | tof_10 | 36 | 37 | 1 | 20 | verify_qasm clean | `agents/c2/out_clean/tof_10_T36_CX37_A1.qasm` | AND compiler (clean ancillas) |
| clean | tof_10 | 39 | 38 | 0 | 19 | verify_qasm clean | `agents/c2/out_clean/tof_10_T39_CX38_A0.qasm` | AND compiler (clean ancillas) |
| clean | tof_3 | 8 | 9 | 1 | 6 | verify_qasm clean | `agents/c2/out_clean/tof_3_T8_CX9_A1.qasm` | AND compiler (clean ancillas) |
| clean | tof_3 | 11 | 10 | 0 | 5 | verify_qasm clean | `agents/c2/out_clean/tof_3_T11_CX10_A0.qasm` | AND compiler (clean ancillas) |
| clean | tof_4 | 12 | 13 | 1 | 8 | verify_qasm clean | `agents/c2/out_clean/tof_4_T12_CX13_A1.qasm` | AND compiler (clean ancillas) |
| clean | tof_4 | 15 | 14 | 0 | 7 | verify_qasm clean | `agents/c2/out_clean/tof_4_T15_CX14_A0.qasm` | AND compiler (clean ancillas) |
| clean | tof_5 | 16 | 17 | 1 | 10 | verify_qasm clean | `agents/c2/out_clean/tof_5_T16_CX17_A1.qasm` | AND compiler (clean ancillas) |
| clean | tof_5 | 19 | 18 | 0 | 9 | verify_qasm clean | `agents/c2/out_clean/tof_5_T19_CX18_A0.qasm` | AND compiler (clean ancillas) |
| clean | vbe_adder_3 | 12 | 33 | 3 | 13 | strict verdict (strict implies clean) | `agents/c1/best_c2/vbe_adder_3_T12_CX33_A3.qasm` | AND compiler + CX resynthesis |
| clean | vbe_adder_3 | 19 | 35 | 2 | 12 | strict verdict (strict implies clean) | `agents/c1/best_c2/vbe_adder_3_T19_CX35_A2.qasm` | AND compiler + CX resynthesis |
| clean | vbe_adder_3 | 23 | 36 | 0 | 10 | strict verdict (strict implies clean) | `agents/c8/qasm/vbe_adder_3_T23_CX36_e0.qasm` | partial H gadgets + TODD + CX polish |
| clean | vbe_adder_3 | 24 | 32 | 0 | 10 | strict verdict (strict implies clean) | `agents/c8/qasm/vbe_adder_3_T24_CX32_e0.qasm` | partial H gadgets + TODD + CX polish |

## Unitary results: provenance

| circuit | T | CX | source | method | checker |
|---|---|---|---|---|---|
| adder_8 | 193 | 223 | `results/qopt_bench.json:adder_8:circ_evo_ours` | portfolio (evolved programs) | numerical (2 random states, complex128, tol 1e-11) |
| barenco_tof_10 | 72 | 72 | `agents/eng1/ours/barenco_tof_10.json` | qopt --mode cx (cached) | verify_qasm strict |
| barenco_tof_3 | 16 | 14 | `results/qopt_bench.json:barenco_tof_3:circ_ours` | portfolio + neural block resynthesis (--net) | verify_qasm strict |
| barenco_tof_4 | 24 | 24 | `results/qopt_bench.json:barenco_tof_4:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| barenco_tof_5 | 32 | 32 | `results/qopt_bench.json:barenco_tof_5:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| csla_mux_3 | 64 | 64 | `results/qopt_bench.json:csla_mux_3:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| csum_mux_9 | 84 | 112 | `results/qopt_bench.json:csum_mux_9:circ_ours` | portfolio + neural block resynthesis (--net) | PyZX verify_equality |
| gf2^4_mult | 68 | 86 | `results/qopt_bench.json:gf2^4_mult:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| gf2^5_mult | 115 | 132 | `results/qopt_bench.json:gf2^5_mult:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| grover_5 | 166 | 160 | `results/qopt_bench.json:grover_5:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| ham15-low | 97 | 179 | `results/qopt_bench.json:ham15-low:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| mod5_4 | 7 | 10 | `results/qopt_bench.json:mod5_4:circ_ours` | portfolio + neural block resynthesis (--net) | verify_qasm strict |
| mod_mult_55 | 35 | 39 | `results/qopt_bench.json:mod_mult_55:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| mod_red_21 | 73 | 61 | `results/qopt_bench.json:mod_red_21:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| qcla_com_7 | 95 | 111 | `results/qopt_bench.json:qcla_com_7:circ_evo_ours` | portfolio (evolved programs) | PyZX verify_equality |
| qcla_mod_7 | 237 | 256 | `results/qopt_bench.json:qcla_mod_7:circ_rt_ours` | portfolio + rt blocks | PyZX verify_equality |
| qft_4 | 67 | 37 | `results/qopt_bench.json:qft_4:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| rc_adder_6 | 47 | 65 | `results/qopt_bench.json:rc_adder_6:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| tof_10 | 71 | 54 | `results/qopt_bench.json:tof_10:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| tof_3 | 15 | 12 | `results/qopt_bench.json:tof_3:circ_ours` | portfolio + neural block resynthesis (--net) | verify_qasm strict |
| tof_4 | 23 | 18 | `results/qopt_bench.json:tof_4:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| tof_5 | 31 | 24 | `results/qopt_bench.json:tof_5:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |
| vbe_adder_3 | 24 | 28 | `results/qopt_bench.json:vbe_adder_3:circ_evo_ours` | portfolio (evolved programs) | verify_qasm strict |

Not included: hwb6 (benchmark file forbids redistribution, see `bench_qc/SOURCES.md`).
Sources are paths inside the (unpublished) research project; they identify which development stage produced
each circuit. The files in `results/` are byte-identical copies (SHA-256 in `results/candidates.jsonl`).
