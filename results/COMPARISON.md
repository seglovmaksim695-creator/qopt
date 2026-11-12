# Comparison with other tools

Third-party numbers come from runs of each tool on the same `.qc` files (scripts in the research project,
agent c11, 2026-11-12). Every circuit a tool produced was checked with `verify_qasm` in the strict and the
clean model; the column *model* is the strongest model in which it passed. *Estimate* rows are resource
estimates without a circuit (Q# Resource Estimator, Qualtran, cirq_ft): they have no CX count and were not verified.

Important caveats:

* "as-is" = the tool receives the gate list of the `.qc` file. "mcx_*", "and_ladder" etc. = we rewrote the
  circuit by hand as a multi-controlled X / ladder of logical ANDs using the tool's own primitives — this needs
  knowledge of what the circuit computes, which the tools do not infer.
* Unitary tools (Qiskit, pytket, Cirq, PyZX, pyLIQTR, PennyLane as-is) solve a harder problem (no measurements,
  no extra qubits); compare them with our *unitary* column, not with the measurement-based ones.
* Published T-count records (TODD/FastTODD, AlphaTensor-Quantum) use H-gadget models without classical
  feed-forward on the output; they are a different model and are listed in the paper, not here.
* Qrisp `mcx(method='jones')` points are invalid (docs/FINDINGS.md F2) and were removed.

| circuit | tool | variant | T | CX | extra qubits | model | source |
|---|---|---|---|---|---|---|---|
| adder_8 | **qopt 0.2.0** | unitary best | **193** | **223** | +0 | unitary | this release, re-verified |
| barenco_tof_10 | **qopt 0.2.0** | unitary best | **72** | **72** | +0 | unitary | this release, re-verified |
| barenco_tof_10 | **qopt 0.2.0** | strict best | **59** | **48** | +3 | strict | this release, re-verified |
| barenco_tof_10 | **qopt 0.2.0** | clean best | **59** | **48** | +3 | clean | this release, re-verified |
| barenco_tof_3 | **qopt 0.2.0** | unitary best | **16** | **14** | +0 | unitary | this release, re-verified |
| barenco_tof_3 | **qopt 0.2.0** | strict best | **8** | **9** | +2 | strict | this release, re-verified |
| barenco_tof_3 | **qopt 0.2.0** | clean best | **8** | **9** | +2 | clean | this release, re-verified |
| barenco_tof_3 | Cirq 1.7 | decompose | 28 | 32 | +0 | strict | run, circuit verified |
| barenco_tof_3 | PennyLane 0.45.1 | and_ladder | 8 | 9 | +2 | strict | run, circuit verified |
| barenco_tof_3 | PyZX 0.10.8 | full_reduce | 16 | 33 | +0 | strict | run, circuit verified |
| barenco_tof_3 | Q# RE (qsharp 1.31) | as-is | 16 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_3 | Q# RE (qsharp 1.31) | mcx_controlledX | 12 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_3 | Q# RE (qsharp 1.31) | mcx_AND_ladder | 8 | — | 1 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_3 | Qiskit 2.6.1 | fresh_1_clean_kg24 | 15 | 12 | +1 | strict | run, circuit verified |
| barenco_tof_3 | Qrisp 0.9.10 | as-is | 28 | 24 | +0 | strict | run, circuit verified |
| barenco_tof_3 | Qualtran 0.7.0 | mcx_and | 8 | 15 | +2 | strict | run, circuit verified |
| barenco_tof_3 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 28 | — | 0 | estimate | run, estimate only |
| barenco_tof_3 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 16 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| barenco_tof_3 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 28 | 24 | +0 | strict | run, circuit verified |
| barenco_tof_3 | pyLIQTR/cirq_ft | MultiControlPauli | 8 | — | 2 | estimate (t_complexity) | run, estimate only |
| barenco_tof_3 | pytket 2.18.5 | zxperm_fpo | 16 | 37 | +0 | strict | run, circuit verified |
| barenco_tof_4 | **qopt 0.2.0** | unitary best | **24** | **24** | +0 | unitary | this release, re-verified |
| barenco_tof_4 | **qopt 0.2.0** | strict best | **12** | **13** | +3 | strict | this release, re-verified |
| barenco_tof_4 | **qopt 0.2.0** | clean best | **12** | **13** | +3 | clean | this release, re-verified |
| barenco_tof_4 | Cirq 1.7 | decompose | 56 | 64 | +0 | strict | run, circuit verified |
| barenco_tof_4 | PennyLane 0.45.1 | and_ladder | 12 | 13 | +3 | strict | run, circuit verified |
| barenco_tof_4 | PyZX 0.10.8 | full_reduce | 28 | 52 | +0 | strict | run, circuit verified |
| barenco_tof_4 | Q# RE (qsharp 1.31) | as-is | 32 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_4 | Q# RE (qsharp 1.31) | mcx_controlledX | 20 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_4 | Q# RE (qsharp 1.31) | mcx_AND_ladder | 12 | — | 1 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_4 | Qiskit 2.6.1 | fresh_1_clean_kg24 | 23 | 18 | +1 | strict | run, circuit verified |
| barenco_tof_4 | Qrisp 0.9.10 | as-is | 56 | 48 | +0 | strict | run, circuit verified |
| barenco_tof_4 | Qualtran 0.7.0 | mcx_and | 12 | 22 | +3 | strict | run, circuit verified |
| barenco_tof_4 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 56 | — | 0 | estimate | run, estimate only |
| barenco_tof_4 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 32 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| barenco_tof_4 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 56 | 48 | +0 | strict | run, circuit verified |
| barenco_tof_4 | pyLIQTR/cirq_ft | MultiControlPauli | 12 | — | 3 | estimate (t_complexity) | run, estimate only |
| barenco_tof_4 | pytket 2.18.5 | zxperm_fpo | 28 | 66 | +0 | strict | run, circuit verified |
| barenco_tof_5 | **qopt 0.2.0** | unitary best | **32** | **32** | +0 | unitary | this release, re-verified |
| barenco_tof_5 | **qopt 0.2.0** | strict best | **16** | **17** | +4 | strict | this release, re-verified |
| barenco_tof_5 | **qopt 0.2.0** | clean best | **16** | **17** | +4 | clean | this release, re-verified |
| barenco_tof_5 | Cirq 1.7 | decompose | 84 | 96 | +0 | strict | run, circuit verified |
| barenco_tof_5 | PennyLane 0.45.1 | and_ladder | 16 | 17 | +4 | strict | run, circuit verified |
| barenco_tof_5 | PyZX 0.10.8 | full_reduce | 40 | 68 | +0 | strict | run, circuit verified |
| barenco_tof_5 | Q# RE (qsharp 1.31) | as-is | 48 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_5 | Q# RE (qsharp 1.31) | mcx_controlledX | 28 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_5 | Q# RE (qsharp 1.31) | mcx_AND_ladder | 16 | — | 1 | estimate (CCZ = 4 T states) | run, estimate only |
| barenco_tof_5 | Qiskit 2.6.1 | fresh_1_clean_kg24 | 31 | 24 | +1 | strict | run, circuit verified |
| barenco_tof_5 | Qrisp 0.9.10 | as-is | 84 | 72 | +0 | strict | run, circuit verified |
| barenco_tof_5 | Qualtran 0.7.0 | mcx_and | 16 | 29 | +4 | strict | run, circuit verified |
| barenco_tof_5 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 84 | — | 0 | estimate | run, estimate only |
| barenco_tof_5 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 48 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| barenco_tof_5 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 84 | 72 | +0 | strict | run, circuit verified |
| barenco_tof_5 | pyLIQTR/cirq_ft | MultiControlPauli | 16 | — | 4 | estimate (t_complexity) | run, estimate only |
| barenco_tof_5 | pytket 2.18.5 | zxperm_fpo | 40 | 106 | +0 | strict | run, circuit verified |
| csla_mux_3 | **qopt 0.2.0** | unitary best | **64** | **64** | +0 | unitary | this release, re-verified |
| csla_mux_3 | **qopt 0.2.0** | strict best | **40** | **70** | +1 | strict | this release, re-verified |
| csla_mux_3 | **qopt 0.2.0** | clean best | **28** | **84** | +3 | clean | this release, re-verified |
| csla_mux_3 | Cirq 1.7 | decompose | 70 | 100 | +0 | strict | run, circuit verified |
| csla_mux_3 | PennyLane 0.45.1 | as-is | 70 | 80 | +0 | strict | run, circuit verified |
| csla_mux_3 | PyZX 0.10.8 | full_reduce | 62 | 161 | +0 | strict | run, circuit verified |
| csla_mux_3 | Q# RE (qsharp 1.31) | as-is | 40 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| csla_mux_3 | Qiskit 2.6.1 | as-is | 64 | 72 | +0 | strict | run, circuit verified |
| csla_mux_3 | Qrisp 0.9.10 | as-is | 70 | 80 | +0 | strict | run, circuit verified |
| csla_mux_3 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 70 | — | 0 | estimate | run, estimate only |
| csla_mux_3 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 40 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| csla_mux_3 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 70 | 80 | +0 | strict | run, circuit verified |
| csla_mux_3 | pytket 2.18.5 | zxperm_fpo | 62 | 178 | +0 | strict | run, circuit verified |
| csum_mux_9 | **qopt 0.2.0** | unitary best | **84** | **112** | +0 | unitary | this release, re-verified |
| gf2^4_mult | **qopt 0.2.0** | unitary best | **68** | **86** | +0 | unitary | this release, re-verified |
| gf2^4_mult | **qopt 0.2.0** | strict best | **58** | **115** | +0 | strict | this release, re-verified |
| gf2^4_mult | **qopt 0.2.0** | clean best | **36** | **79** | +1 | clean | this release, re-verified |
| gf2^4_mult | Cirq 1.7 | decompose | 112 | 131 | +0 | strict | run, circuit verified |
| gf2^4_mult | PennyLane 0.45.1 | as-is | 112 | 99 | +0 | strict | run, circuit verified |
| gf2^4_mult | PyZX 0.10.8 | full_reduce | 68 | 274 | +0 | strict | run, circuit verified |
| gf2^4_mult | Q# RE (qsharp 1.31) | as-is | 64 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| gf2^4_mult | Qiskit 2.6.1 | as-is | 96 | 99 | +0 | strict | run, circuit verified |
| gf2^4_mult | Qrisp 0.9.10 | as-is | 112 | 99 | +0 | strict | run, circuit verified |
| gf2^4_mult | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 112 | — | 0 | estimate | run, estimate only |
| gf2^4_mult | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 64 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| gf2^4_mult | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 112 | 99 | +0 | strict | run, circuit verified |
| gf2^4_mult | pytket 2.18.5 | zxperm_fpo | 68 | 302 | +0 | strict | run, circuit verified |
| gf2^5_mult | **qopt 0.2.0** | unitary best | **115** | **132** | +0 | unitary | this release, re-verified |
| gf2^5_mult | **qopt 0.2.0** | strict best | **90** | **183** | +0 | strict | this release, re-verified |
| gf2^5_mult | **qopt 0.2.0** | clean best | **56** | **124** | +1 | clean | this release, re-verified |
| gf2^5_mult | Cirq 1.7 | decompose | 175 | 204 | +0 | strict | run, circuit verified |
| gf2^5_mult | PennyLane 0.45.1 | as-is | 175 | 154 | +0 | strict | run, circuit verified |
| gf2^5_mult | PyZX 0.10.8 | full_reduce | 115 | 518 | +0 | strict | run, circuit verified |
| gf2^5_mult | Q# RE (qsharp 1.31) | as-is | 100 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| gf2^5_mult | Qiskit 2.6.1 | as-is | 155 | 154 | +0 | strict | run, circuit verified |
| gf2^5_mult | Qrisp 0.9.10 | as-is | 175 | 154 | +0 | strict | run, circuit verified |
| gf2^5_mult | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 175 | — | 0 | estimate | run, estimate only |
| gf2^5_mult | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 100 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| gf2^5_mult | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 175 | 154 | +0 | strict | run, circuit verified |
| gf2^5_mult | pytket 2.18.5 | zxperm_fpo | 115 | 509 | +0 | strict | run, circuit verified |
| grover_5 | **qopt 0.2.0** | unitary best | **166** | **160** | +0 | unitary | this release, re-verified |
| grover_5 | **qopt 0.2.0** | strict best | **124** | **155** | +4 | strict | this release, re-verified |
| grover_5 | **qopt 0.2.0** | clean best | **124** | **155** | +4 | clean | this release, re-verified |
| grover_5 | Cirq 1.7 | decompose | 336 | 384 | +0 | strict | run, circuit verified |
| grover_5 | PennyLane 0.45.1 | as-is | 336 | 288 | +0 | strict | run, circuit verified |
| grover_5 | PyZX 0.10.8 | full_reduce | 166 | 301 | +0 | strict | run, circuit verified |
| grover_5 | Q# RE (qsharp 1.31) | as-is | 192 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| grover_5 | Qiskit 2.6.1 | as-is | 296 | 288 | +0 | strict | run, circuit verified |
| grover_5 | Qrisp 0.9.10 | as-is | 336 | 288 | +0 | strict | run, circuit verified |
| grover_5 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 336 | — | 0 | estimate | run, estimate only |
| grover_5 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 192 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| grover_5 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 336 | 288 | +0 | strict | run, circuit verified |
| grover_5 | pytket 2.18.5 | zxperm | 166 | 391 | +0 | strict | run, circuit verified |
| ham15-low | **qopt 0.2.0** | unitary best | **97** | **179** | +0 | unitary | this release, re-verified |
| ham15-low | **qopt 0.2.0** | clean best | **68** | **177** | +1 | clean | this release, re-verified |
| ham15-low | Cirq 1.7 | decompose | 161 | 282 | +0 | strict | run, circuit verified |
| ham15-low | PennyLane 0.45.1 | as-is | 161 | 236 | +0 | strict | run, circuit verified |
| ham15-low | PyZX 0.10.8 | full_reduce | 97 | 322 | +0 | strict | run, circuit verified |
| ham15-low | Q# RE (qsharp 1.31) | as-is | 92 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| ham15-low | Qiskit 2.6.1 | as-is | 147 | 236 | +0 | strict | run, circuit verified |
| ham15-low | Qrisp 0.9.10 | as-is | 161 | 236 | +0 | strict | run, circuit verified |
| ham15-low | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 161 | — | 0 | estimate | run, estimate only |
| ham15-low | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 92 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| ham15-low | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 161 | 236 | +0 | strict | run, circuit verified |
| ham15-low | pytket 2.18.5 | zxperm | 97 | 636 | +0 | strict | run, circuit verified |
| mod5_4 | **qopt 0.2.0** | unitary best | **7** | **10** | +0 | unitary | this release, re-verified |
| mod5_4 | **qopt 0.2.0** | strict best | **4** | **11** | +1 | strict | this release, re-verified |
| mod5_4 | **qopt 0.2.0** | clean best | **4** | **11** | +0 | clean | this release, re-verified |
| mod5_4 | Cirq 1.7 | decompose | 28 | 36 | +0 | strict | run, circuit verified |
| mod5_4 | PennyLane 0.45.1 | as-is | 28 | 28 | +0 | strict | run, circuit verified |
| mod5_4 | PyZX 0.10.8 | full_reduce | 8 | 20 | +0 | strict | run, circuit verified |
| mod5_4 | Q# RE (qsharp 1.31) | as-is | 16 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| mod5_4 | Qiskit 2.6.1 | as-is | 22 | 28 | +0 | strict | run, circuit verified |
| mod5_4 | Qrisp 0.9.10 | as-is | 28 | 28 | +0 | strict | run, circuit verified |
| mod5_4 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 28 | — | 0 | estimate | run, estimate only |
| mod5_4 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 16 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| mod5_4 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 28 | 28 | +0 | strict | run, circuit verified |
| mod5_4 | pytket 2.18.5 | zxperm_fpo | 8 | 27 | +0 | strict | run, circuit verified |
| mod_mult_55 | **qopt 0.2.0** | unitary best | **35** | **39** | +0 | unitary | this release, re-verified |
| mod_mult_55 | **qopt 0.2.0** | strict best | **20** | **35** | +2 | strict | this release, re-verified |
| mod_mult_55 | **qopt 0.2.0** | clean best | **12** | **24** | +1 | clean | this release, re-verified |
| mod_mult_55 | Cirq 1.7 | decompose | 49 | 62 | +0 | strict | run, circuit verified |
| mod_mult_55 | PennyLane 0.45.1 | as-is | 49 | 48 | +0 | strict | run, circuit verified |
| mod_mult_55 | PyZX 0.10.8 | full_reduce | 35 | 89 | +0 | strict | run, circuit verified |
| mod_mult_55 | Q# RE (qsharp 1.31) | as-is | 28 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| mod_mult_55 | Qiskit 2.6.1 | as-is | 45 | 48 | +0 | strict | run, circuit verified |
| mod_mult_55 | Qrisp 0.9.10 | as-is | 49 | 48 | +0 | strict | run, circuit verified |
| mod_mult_55 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 49 | — | 0 | estimate | run, estimate only |
| mod_mult_55 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 28 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| mod_mult_55 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 49 | 48 | +0 | strict | run, circuit verified |
| mod_mult_55 | pytket 2.18.5 | zxperm_fpo | 39 | 104 | +0 | strict | run, circuit verified |
| mod_red_21 | **qopt 0.2.0** | unitary best | **73** | **61** | +0 | unitary | this release, re-verified |
| mod_red_21 | **qopt 0.2.0** | strict best | **44** | **58** | +5 | strict | this release, re-verified |
| mod_red_21 | **qopt 0.2.0** | clean best | **28** | **48** | +0 | clean | this release, re-verified |
| mod_red_21 | Cirq 1.7 | decompose | 119 | 139 | +0 | strict | run, circuit verified |
| mod_red_21 | PennyLane 0.45.1 | as-is | 119 | 105 | +0 | strict | run, circuit verified |
| mod_red_21 | PyZX 0.10.8 | full_reduce | 73 | 149 | +0 | strict | run, circuit verified |
| mod_red_21 | Q# RE (qsharp 1.31) | as-is | 68 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| mod_red_21 | Qiskit 2.6.1 | as-is | 107 | 105 | +0 | strict | run, circuit verified |
| mod_red_21 | Qrisp 0.9.10 | as-is | 119 | 105 | +0 | strict | run, circuit verified |
| mod_red_21 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 119 | — | 0 | estimate | run, estimate only |
| mod_red_21 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 68 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| mod_red_21 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 119 | 105 | +0 | strict | run, circuit verified |
| mod_red_21 | pytket 2.18.5 | zxperm | 77 | 177 | +0 | strict | run, circuit verified |
| qcla_com_7 | **qopt 0.2.0** | unitary best | **95** | **111** | +0 | unitary | this release, re-verified |
| qcla_com_7 | **qopt 0.2.0** | strict best | **59** | **565** | +37 | strict | this release, re-verified |
| qcla_com_7 | **qopt 0.2.0** | clean best | **59** | **565** | +37 | clean | this release, re-verified |
| qcla_mod_7 | **qopt 0.2.0** | unitary best | **237** | **256** | +0 | unitary | this release, re-verified |
| qft_4 | **qopt 0.2.0** | unitary best | **67** | **37** | +0 | unitary | this release, re-verified |
| qft_4 | **qopt 0.2.0** | strict best | **53** | **265** | +42 | strict | this release, re-verified |
| qft_4 | **qopt 0.2.0** | clean best | **53** | **265** | +42 | clean | this release, re-verified |
| qft_4 | Cirq 1.7 | decompose | 69 | 50 | +0 | strict | run, circuit verified |
| qft_4 | PennyLane 0.45.1 | as-is | 69 | 46 | +0 | strict | run, circuit verified |
| qft_4 | PyZX 0.10.8 | full_reduce | 67 | 55 | +0 | strict | run, circuit verified |
| qft_4 | Q# RE (qsharp 1.31) | as-is | 63 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| qft_4 | Qiskit 2.6.1 | as-is | 409 | 45 | +0 | strict | run, circuit verified |
| qft_4 | Qrisp 0.9.10 | as-is | 69 | 46 | +0 | strict | run, circuit verified |
| qft_4 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 69 | — | 0 | estimate | run, estimate only |
| qft_4 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 63 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| qft_4 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 69 | 46 | +0 | strict | run, circuit verified |
| qft_4 | pytket 2.18.5 | zxperm | 67 | 72 | +0 | strict | run, circuit verified |
| rc_adder_6 | **qopt 0.2.0** | unitary best | **47** | **65** | +0 | unitary | this release, re-verified |
| rc_adder_6 | **qopt 0.2.0** | strict best | **24** | **54** | +6 | strict | this release, re-verified |
| rc_adder_6 | **qopt 0.2.0** | clean best | **24** | **53** | +6 | clean | this release, re-verified |
| rc_adder_6 | Cirq 1.7 | decompose | 77 | 115 | +0 | strict | run, circuit verified |
| rc_adder_6 | PennyLane 0.45.1 | as-is | 77 | 93 | +0 | strict | run, circuit verified |
| rc_adder_6 | PyZX 0.10.8 | full_reduce | 47 | 103 | +0 | strict | run, circuit verified |
| rc_adder_6 | Q# RE (qsharp 1.31) | as-is | 44 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| rc_adder_6 | Qiskit 2.6.1 | as-is | 67 | 83 | +0 | strict | run, circuit verified |
| rc_adder_6 | Qrisp 0.9.10 | as-is | 77 | 93 | +0 | strict | run, circuit verified |
| rc_adder_6 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 77 | — | 0 | estimate | run, estimate only |
| rc_adder_6 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 44 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| rc_adder_6 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 77 | 93 | +0 | strict | run, circuit verified |
| rc_adder_6 | pytket 2.18.5 | zxperm_fpo | 55 | 153 | +0 | strict | run, circuit verified |
| tof_10 | **qopt 0.2.0** | unitary best | **71** | **54** | +0 | unitary | this release, re-verified |
| tof_10 | **qopt 0.2.0** | strict best | **59** | **71** | +3 | strict | this release, re-verified |
| tof_10 | **qopt 0.2.0** | clean best | **36** | **37** | +1 | clean | this release, re-verified |
| tof_10 | Cirq 1.7 | decompose | 119 | 136 | +0 | strict | run, circuit verified |
| tof_10 | PennyLane 0.45.1 | as-is | 119 | 102 | +0 | strict | run, circuit verified |
| tof_10 | PennyLane 0.45.1 | and_ladder | 36 | 37 | +1 | clean | run, circuit verified |
| tof_10 | PyZX 0.10.8 | full_reduce | 71 | 161 | +0 | strict | run, circuit verified |
| tof_10 | Q# RE (qsharp 1.31) | as-is | 68 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_10 | Q# RE (qsharp 1.31) | mcx_controlledX | 68 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_10 | Q# RE (qsharp 1.31) | mcx_AND_ladder | 36 | — | 1 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_10 | Qiskit 2.6.1 | as-is | 103 | 102 | +0 | strict | run, circuit verified |
| tof_10 | Qiskit 2.6.1 | mcx_n_clean_m15 | 71 | 54 | +0 | clean | run, circuit verified |
| tof_10 | Qrisp 0.9.10 | as-is | 119 | 102 | +0 | strict | run, circuit verified |
| tof_10 | Qrisp 0.9.10 | mcx_gidney_ladder | 36 | 64 | +1 | clean | run, circuit verified |
| tof_10 | Qualtran 0.7.0 | mcx_and | 36 | 64 | +1 | clean | run, circuit verified |
| tof_10 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 119 | — | 0 | estimate | run, estimate only |
| tof_10 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 68 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| tof_10 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 119 | 102 | +0 | strict | run, circuit verified |
| tof_10 | pyLIQTR/cirq_ft | MultiControlPauli | 36 | — | 9 | estimate (t_complexity) | run, estimate only |
| tof_10 | pytket 2.18.5 | zxperm | 71 | 168 | +0 | strict | run, circuit verified |
| tof_3 | **qopt 0.2.0** | unitary best | **15** | **12** | +0 | unitary | this release, re-verified |
| tof_3 | **qopt 0.2.0** | strict best | **8** | **10** | +2 | strict | this release, re-verified |
| tof_3 | **qopt 0.2.0** | clean best | **8** | **9** | +1 | clean | this release, re-verified |
| tof_3 | Cirq 1.7 | decompose | 21 | 24 | +0 | strict | run, circuit verified |
| tof_3 | PennyLane 0.45.1 | as-is | 21 | 18 | +0 | strict | run, circuit verified |
| tof_3 | PennyLane 0.45.1 | and_ladder | 8 | 9 | +1 | clean | run, circuit verified |
| tof_3 | PyZX 0.10.8 | full_reduce | 15 | 26 | +0 | strict | run, circuit verified |
| tof_3 | Q# RE (qsharp 1.31) | as-is | 12 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_3 | Q# RE (qsharp 1.31) | mcx_controlledX | 12 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_3 | Q# RE (qsharp 1.31) | mcx_AND_ladder | 8 | — | 1 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_3 | Qiskit 2.6.1 | mcx_n_clean_m15 | 15 | 12 | +0 | strict | run, circuit verified |
| tof_3 | Qrisp 0.9.10 | mcx_maslov | 15 | 12 | +0 | strict | run, circuit verified |
| tof_3 | Qrisp 0.9.10 | mcx_gidney_ladder | 8 | 15 | +1 | clean | run, circuit verified |
| tof_3 | Qualtran 0.7.0 | mcx_and | 8 | 15 | +1 | clean | run, circuit verified |
| tof_3 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 21 | — | 0 | estimate | run, estimate only |
| tof_3 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 12 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| tof_3 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 21 | 18 | +0 | strict | run, circuit verified |
| tof_3 | pyLIQTR/cirq_ft | MultiControlPauli | 8 | — | 2 | estimate (t_complexity) | run, estimate only |
| tof_3 | pytket 2.18.5 | zxperm | 15 | 37 | +0 | strict | run, circuit verified |
| tof_4 | **qopt 0.2.0** | unitary best | **23** | **18** | +0 | unitary | this release, re-verified |
| tof_4 | **qopt 0.2.0** | strict best | **12** | **15** | +3 | strict | this release, re-verified |
| tof_4 | **qopt 0.2.0** | clean best | **12** | **13** | +1 | clean | this release, re-verified |
| tof_4 | Cirq 1.7 | decompose | 35 | 40 | +0 | strict | run, circuit verified |
| tof_4 | PennyLane 0.45.1 | as-is | 35 | 30 | +0 | strict | run, circuit verified |
| tof_4 | PennyLane 0.45.1 | and_ladder | 12 | 13 | +1 | clean | run, circuit verified |
| tof_4 | PyZX 0.10.8 | full_reduce | 23 | 42 | +0 | strict | run, circuit verified |
| tof_4 | Q# RE (qsharp 1.31) | as-is | 20 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_4 | Q# RE (qsharp 1.31) | mcx_controlledX | 20 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_4 | Q# RE (qsharp 1.31) | mcx_AND_ladder | 12 | — | 1 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_4 | Qiskit 2.6.1 | as-is | 31 | 30 | +0 | strict | run, circuit verified |
| tof_4 | Qiskit 2.6.1 | mcx_n_clean_m15 | 23 | 18 | +0 | clean | run, circuit verified |
| tof_4 | Qrisp 0.9.10 | as-is | 35 | 30 | +0 | strict | run, circuit verified |
| tof_4 | Qrisp 0.9.10 | mcx_gidney_ladder | 12 | 22 | +1 | clean | run, circuit verified |
| tof_4 | Qualtran 0.7.0 | mcx_and | 12 | 22 | +1 | clean | run, circuit verified |
| tof_4 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 35 | — | 0 | estimate | run, estimate only |
| tof_4 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 20 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| tof_4 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 35 | 30 | +0 | strict | run, circuit verified |
| tof_4 | pyLIQTR/cirq_ft | MultiControlPauli | 12 | — | 3 | estimate (t_complexity) | run, estimate only |
| tof_4 | pytket 2.18.5 | zxperm | 23 | 54 | +0 | strict | run, circuit verified |
| tof_5 | **qopt 0.2.0** | unitary best | **31** | **24** | +0 | unitary | this release, re-verified |
| tof_5 | **qopt 0.2.0** | strict best | **16** | **20** | +4 | strict | this release, re-verified |
| tof_5 | **qopt 0.2.0** | clean best | **16** | **17** | +1 | clean | this release, re-verified |
| tof_5 | Cirq 1.7 | decompose | 49 | 56 | +0 | strict | run, circuit verified |
| tof_5 | PennyLane 0.45.1 | as-is | 49 | 42 | +0 | strict | run, circuit verified |
| tof_5 | PennyLane 0.45.1 | and_ladder | 16 | 17 | +1 | clean | run, circuit verified |
| tof_5 | PyZX 0.10.8 | full_reduce | 31 | 59 | +0 | strict | run, circuit verified |
| tof_5 | Q# RE (qsharp 1.31) | as-is | 28 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_5 | Q# RE (qsharp 1.31) | mcx_controlledX | 28 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_5 | Q# RE (qsharp 1.31) | mcx_AND_ladder | 16 | — | 1 | estimate (CCZ = 4 T states) | run, estimate only |
| tof_5 | Qiskit 2.6.1 | as-is | 43 | 42 | +0 | strict | run, circuit verified |
| tof_5 | Qiskit 2.6.1 | mcx_n_clean_m15 | 31 | 24 | +0 | clean | run, circuit verified |
| tof_5 | Qrisp 0.9.10 | as-is | 49 | 42 | +0 | strict | run, circuit verified |
| tof_5 | Qrisp 0.9.10 | mcx_gidney_ladder | 16 | 29 | +1 | clean | run, circuit verified |
| tof_5 | Qualtran 0.7.0 | mcx_and | 16 | 29 | +1 | clean | run, circuit verified |
| tof_5 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 49 | — | 0 | estimate | run, estimate only |
| tof_5 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 28 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| tof_5 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 49 | 42 | +0 | strict | run, circuit verified |
| tof_5 | pyLIQTR/cirq_ft | MultiControlPauli | 16 | — | 4 | estimate (t_complexity) | run, estimate only |
| tof_5 | pytket 2.18.5 | zxperm | 31 | 73 | +0 | strict | run, circuit verified |
| vbe_adder_3 | **qopt 0.2.0** | unitary best | **24** | **28** | +0 | unitary | this release, re-verified |
| vbe_adder_3 | **qopt 0.2.0** | strict best | **12** | **33** | +3 | strict | this release, re-verified |
| vbe_adder_3 | **qopt 0.2.0** | clean best | **12** | **33** | +3 | clean | this release, re-verified |
| vbe_adder_3 | Cirq 1.7 | decompose | 70 | 90 | +0 | strict | run, circuit verified |
| vbe_adder_3 | PennyLane 0.45.1 | as-is | 70 | 70 | +0 | strict | run, circuit verified |
| vbe_adder_3 | PyZX 0.10.8 | full_reduce | 24 | 71 | +0 | strict | run, circuit verified |
| vbe_adder_3 | Q# RE (qsharp 1.31) | as-is | 40 | — | 0 | estimate (CCZ = 4 T states) | run, estimate only |
| vbe_adder_3 | Qiskit 2.6.1 | as-is | 56 | 58 | +0 | strict | run, circuit verified |
| vbe_adder_3 | Qrisp 0.9.10 | as-is | 70 | 70 | +0 | strict | run, circuit verified |
| vbe_adder_3 | Qualtran 0.7.0 | as-is (CCZ via Cirq, 7 T) | 70 | — | 0 | estimate | run, estimate only |
| vbe_adder_3 | Qualtran 0.7.0 | as-is, Toffoli bloq (4 T) | 40 | — | — | estimate (Toffoli = 4 T; needs ancilla + measurement) | computed from documented cost |
| vbe_adder_3 | pyLIQTR 1.0.1 | clifford_plus_t_direct_transform | 70 | 70 | +0 | strict | run, circuit verified |
| vbe_adder_3 | pytket 2.18.5 | zxperm_fpo | 24 | 81 | +0 | strict | run, circuit verified |
