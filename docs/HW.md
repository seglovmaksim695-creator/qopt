# `qoptct hw` — compiling for noisy IBM-like hardware

Status: experimental (branch `f-route`). Everything below was measured on *fake* backends, i.e. stored calibration
snapshots of real IBM devices (`qiskit_ibm_runtime.fake_provider`), and with a noisy simulation. The only real-chip
measurement is the small check in the section "Check on a real chip" below.

## What it does

Input: a unitary circuit (OpenQASM 2/3, a Feynman `.qc` file or a `QuantumCircuit`). Output: an equivalent circuit
on the physical qubits of the device in its native gates (CZ or ECR, RZ, SX, X), the qubit maps
(`init_map[i]`: where logical qubit *i* starts, `final_map[i]`: where it must be measured), metrics and the
verification status.

```
qoptct hw --backend fake_torino in.qasm -o out.qasm --report out.json [--seeds 5] [--seconds 20]
          [--objective fid|twoq|depth] [--strategies qiskit,qopt+qiskit,hw,qopt+hw] [--measure]
qoptct hw --list-backends
```
```python
from qoptct.hw import compile_hw
r = compile_hw(qc, "fake_torino", seeds=range(5), seconds=20)
```

Exit code of the command: 0 = verified result written, 2 = written but not verified (more than 26 used qubits),
3 = input error.

## Method

**Metric** (`qoptct/hw/metrics.py`). Two-qubit gate count, two-qubit depth, duration, and the expected fidelity
`fid` = product of (1 − error) over all gates × product of (1 − readout error) over the final positions of the
logical qubits × T1/T2 damping of every used qubit during its idle time (as-soon-as-possible schedule with the
calibrated durations; average fidelity of an idle period *t*: 1/2 + e^(−t/T1)/6 + e^(−t/T2)/3). It is the usual
first-order estimate: no crosstalk, no coherent errors. It is used to rank circuits, not to predict device results.

**Candidates** (`qoptct/hw/compiler.py`). The best candidate by `fid` that passes verification is returned.

| name | what |
|---|---|
| `qiskit` | `transpile(optimization_level=3)`, one run per seed. Always included. |
| `qopt+qiskit` | qoptct's proven CX/T optimization of the Clifford+T parts of the logical circuit (`QoptPass`), then Qiskit level 3. |
| `hw` | the qoptct router (below), then block consolidation, exact two-qubit synthesis, Qiskit's level-3 optimization loop, relocation. |
| `qopt+hw` | the same router on the qoptct-optimized logical circuit. |

**Router** (`qoptct/hw/router.py`). A SABRE-style SWAP router with these additions:

1. *Commutation-aware dependencies.* Each gate gets, per wire, an axis label (`z` if it is block-diagonal in the Z
   basis of that wire: RZ, CZ, CP, RZZ, the control of CX; `x`, `y` likewise). Gates whose labels agree on every
   shared wire commute, so a gate depends only on the previous different-axis group on each wire. All ZZ terms of a
   QAOA layer, or all controlled phases onto one QFT qubit, are in the front layer together and are executed in the
   order that suits the device.
2. *Exact SWAP cost.* A SWAP next to a gate on the same pair is later merged with it into one block of at most 3
   native gates. The router keeps the 4×4 unitary of the live block on each coupler and charges a SWAP its marginal
   cost from the KAK decomposition: +0 after a generic block (XX+YY+ZZ), +1 after CX / CP / RZZ, +3 on a fresh pair.
3. *Floating blocks.* A block stays live while every later gate on its wires commutes with it, so a SWAP can merge
   with a gate that was emitted long before; the emitted order is rearranged so that the block is contiguous.
4. *Mirror rule.* Right after a gate, its pair is swapped if the SWAP score is below a threshold. On a line this
   reproduces (and slightly improves) the known linear-nearest-neighbour QFT (3 native gates per qubit pair).
5. *Layouts.* Interaction-graph placement on a low-error region refined by forward/backward passes; logical qubits
   in index order (or in an order derived from the interaction graph) on a low-error path of the device; the
   layouts found by Qiskit.
6. *Search.* A small bandit over (router variant, layout kind) within the time budget; routing one 20-qubit circuit
   takes 0.01–0.2 s in pure Python, so hundreds to thousands of trials fit in 20 s. The best trials by gate count
   and by calibration-weighted cost are finished to native gates and scored with the metric above.
7. *Relocation.* The finished circuit is moved to the isomorphic region of the device with the best calibration
   (same idea as Qiskit's VF2PostLayout).

**Verification** (`qoptct/hw/check.py`). For random product states |ψ⟩ of the logical qubits placed on `init_map`
(all other used qubits in |0⟩) the circuit *prepare → compiled → inverse(logical) on `final_map` → inverse(prepare)*
is simulated as a state vector on the used physical qubits; the amplitude of |0…0⟩ must have modulus 1. Two
different unitaries agree on a random product state with probability 0, so this is a numerical proof up to
floating-point tolerance. Limit: 26 used qubits (state vector); above that the status is "not verified".

## Results

fake_torino, 76 circuits (arithmetic, QFT, Grover, Hamiltonian simulation, QAOA, ...), 5-30 qubits; geometric means of the
best qoptct strategy against Qiskit O3 best of 10 seeds (model figures from the chip model; on a real chip the measured
gain in fidelity was 0-5 % on small circuits, see below): arithmetic -34% two-qubit gates / x1.47 expected fidelity,
Toffoli-level arithmetic -41% / x1.85, Hamiltonian simulation -19% / x1.48, QAOA -24% / x1.38, QFT -16% / x1.38,
Grover -14% / x1.25, misc -29% / x1.09. Higher expected fidelity than Qiskit on 70 of 75 verified circuits, than tket
on 72 of 74. Aer noisy simulation (29 small circuits): Hellinger fidelity x1.14, agrees with the model's ranking on 27.
fake_sherbrooke (ECR), 31 smallest circuits: same picture. Full tables: `REPORT.md` and `runs/*_table.md` (not in git).

## Check on a real chip (ibm_kingston, 20.11.2026, one job per experiment)

Two experiments, one Sampler job each, success criteria written down before submission
(`docs/real_runs/2026-11-20_ibm_kingston/`). Statistical error: about 0.01 (2000 shots) and 0.025 (400 shots).

* **`qoptct hw` against Qiskit O3 (best of 5 seeds), 4 small circuits.** qoptct is better by +5 % (QFT-5) and +4 %
  (adder) in Hellinger fidelity, with fewer two-qubit gates (25 against 30; 27 against 43); on QAOA-5 and the Heisenberg
  chain the difference is within noise. Mean ratio 1.02. By the pre-set criterion (success needed at least 3 of 4
  circuits and a mean ratio of at least 1.05) the result is a **draw**, not a success.
* **Measurement-based erasure on the chip loses to ordinary uncomputation on every block checked** (adder2, adder3,
  cmp2); the hybrid lost as well (lookadd2x2, X test: 0.682 against 0.772). The model badly underestimates the cost of a
  measurement in the middle of a circuit (adder2: predicted 0.918, measured 0.525). The hybrid won only in the model.
  Therefore dynamic variants are no longer chosen by the model: only by the explicit option `--allow-dynamic`
  (`allow_dynamic=True`), which prints a warning.
* **The calibration estimate stays optimistic** (adder after `qoptct hw`: estimate 0.877, measured 0.781). Use it to rank
  variants, not as an absolute number.

Protocol, results and the plan of the real run: `docs/real_runs/2026-11-20_ibm_kingston/`.

## Limitations

* Check against real runs: on 8 real IBM jobs (ibm_fez, ibm_marrakesh; 286 result rows, one calibration snapshot per chip) the calibration estimate understated the error by a factor of about 1.6 (1.6-1.7 for the models compared). As a tool for a preliminary choice between variants it is useful: it ordered 93.7 % of the distinguishable pairs correctly (100 % on Fez, 76.4 % on Marrakesh, where only 6 tasks were available). The absolute fidelity numbers must not be trusted. A fitted correction factor exists but did not pass the accuracy criterion and is not built into the code.

* Unitary circuits only: no mid-circuit measurement, no classical control (a `ValueError` is raised). Final
  measurements are removed and must be put back on `final_map` (`--measure` does it).
* The expected fidelity is a model. The noisy simulation agrees with it in ranking, but both use the same
  calibration snapshot and neither contains crosstalk, leakage or calibration drift. A first hardware run
  (`scripts/hw_real_run.py`, ibm_kingston, 20.11.2026, 4 small circuits) gave a draw by the pre-set criterion: +4-5 % on
  two circuits, within noise on two; larger circuits were not run on the chip.
* Verification stops at 26 used qubits; larger results are returned with status "not verified" (exit code 2).
* The router is pure Python. Its search time is a fixed budget (default 20 s per logical variant), while Qiskit
  level 3 needs 0.05–1 s per seed for these circuits; with equal wall time Qiskit could run more seeds than the 10
  used here (the best-of-10 column shows how little more seeds help it).
* The gains are specific to circuit classes (see the table). For generic two-qubit unitaries without structure
  (quantum volume) and for the dense Heisenberg grids the router is on par with Qiskit or slightly worse; the
  portfolio then returns the Qiskit candidate.
* The logical CX/T optimization helps only Clifford+T (Toffoli-level) circuits and is time-limited
  (`--qopt-seconds`, default 60); for two GF(2^m) multipliers the CX-minimal logical circuit routes *worse* than the
  original (its CNOT structure is less local), and the portfolio discards it.
* Synthetic devices (`line:N`, `heavyhex:D`, ...) have random, not realistic, calibration values; they exist for tests.
