# Compositional proof of the Kaliski inversion circuit: what the machine checks and what is proved on paper

Circuit: `qoptct.crypto.kaliski.kaliski_inverse(p, reduction, clean)` (default construction, `legacy=False`),
p an odd prime, n = bitlen(p). Input: x in [1, p) on wires 0..n-1, every other wire |0>.

**Theorem.** If `verify_kaliski(qasm, p, reduction, clean)` returns PROVED, then for every x in [1, p) and every
string of measurement outcomes the circuit `qasm` maps |x>|0..0> to

* clean=False: a basis state whose wires `out` hold K = -x^{-1} 2^{2n} mod p, `u` hold 1, `v` hold 0 (the other
  wires hold garbage that is a function of x only);
* clean=True:  |x> |K> |0..0>  (x on wires 0..n-1, K on wires n..2n-1),

with amplitude exactly 2^{-m/2} for each of the 2^m outcome strings (m = number of measurements). Hence on every
superposition of inputs x in [1, p) the circuit is the isometry |x>|0> -> |x, K, ...> for every outcome: there is no
residual phase and no dependence of the final state on the outcomes.

The proof has a machine part (sections 1-3: performed by the checker on the concrete text, at full width, for all
inputs at once) and a paper part (sections 4-6: properties of Kaliski's algorithm as integer arithmetic). Section 7
lists exactly what is trusted.

## 1. Lifting (machine) and the gadget lemmas (machine, exact arithmetic)

`compose.lift` accepts the text only if it is a concatenation of the line groups `x`, `cx`, the 10-line AND gadget and
the 3-line ERASE gadget (measurement in the X basis, classically controlled `cz a b; [z l;] x t`), with distinct wires
inside each gadget. Result: a program of logical operations x(a), cx(a,b), and(a,b -> t), erase(t | a,b,l).

`check_gadget_lemmas()` verifies in exact arithmetic over Q(zeta_8), on the very templates the lifter matches:

* AND: |a,b,0> -> |a,b,ab> with amplitude 1, for all a, b.
* ERASE: for every basis state with t = ab xor l and each outcome r in {0,1}: the result is (1/sqrt2)|a,b,l,0>.

**Lemma 1 (classical reduction; paper, by linearity).** Let a basis input be such that in the classical run of the
lifted program every `and` finds its target equal to 0 and every `erase` finds t = ab xor l. Then, for every outcome
string, the circuit maps this basis state to 2^{-m/2}|classical result>. *Proof.* Induction over the operations: x and
cx are permutations of basis states; by the AND lemma the gadget acts on the current basis state as t := ab; by the
ERASE lemma each outcome maps it to (1/sqrt2) times the basis state with t := 0. []

So it remains to prove a purely classical statement about the lifted program P for all x in [1, p):
(O1) every and-target is 0, (O2) every erase condition holds, (O3) the final values are as specified.

## 2. Segments and contracts (machine: ROBDD proofs on the real ops at the real width)

The certificate (`KaliskiCircuit.hints()`: positions of marks in the op list and names of freshly allocated wires) cuts
the forward part of P into consecutive segments: `init`; for each Kaliski step k: `cmp`, `swap`, `sub`, `add`; for the
generic reduction n+1 segments `div`; for the special reduction `sp_const`, `sp_addL`, `sp_term`..., `sp_q`....
The certificate is not trusted: the walker checks that the segments are consecutive, cover all forward ops, and that
every wire the certificate calls fresh is a valid wire, owned by no register, and known to be zero.

For every segment the walker (`compose.kal.Chain`) forms its **contract** `spec(al, val)`: the values of the wires at
the exit as bit-vector expressions of the values at the entry (`compose.bv`: ripple definitions of +, -, <=, mux),
plus a precondition that consists of

* wires known to be zero at the entry (set Z, maintained by the walker: initially all wires except x; after each
  segment a wire is in Z iff its exit value was proved to be 0 under the precondition),
* *facts* established by earlier segments as proved postconditions, about wires not touched since
  (`lt`: value < constant; `nand`; `imp`; `geodd`: A_0 -> B <= A; `wdef`: w = wspec(A, B, w_prev)),
* *mathematical assumptions* M1..M6 on the entry state (section 4).

`_seg_prove` then creates one ROBDD variable for every entry wire not in Z (all register bits: up to 4n+2 variables),
executes the real ops of the segment symbolically and proves, as ROBDD identities, for ALL entry values satisfying
the precondition: (i) each and-target is 0, (ii) each erase condition holds, (iii) every touched wire ends with the
contract value (wires outside the contract: unchanged, ancillas: 0), (iv) the postfacts, and that the precondition is
satisfiable (no vacuous proof). In clean mode also (v): the mirror slice of the inverse half (located by counting
ops; erase-with-l mirrors to and + cx), started from the exit state, satisfies (i), (ii) and returns every wire to
its entry value.

The contracts (frame registers A, B: n wires; R: n wires; S: n+k wires, at most 2n; w_prev; all values integers):

| segment | contract (exact integer semantics)                                   | precondition                     |
|---------|-----------------------------------------------------------------------|----------------------------------|
| init    | A = p, S = 1 (B = x, R = 0)                                           | -                                |
| cmp k   | w = wspec(A,B,w_prev) = [u odd] and ([v even] or [u <= v]), (u,v) = (B,A) if w_prev else (A,B) | M1: A, B < 2^wa   |
| swap k  | d = w_prev xor w; if d: A<->B, R<->S. Postfacts: A_0 -> B_0; A_0 -> B <= A; not(R_0 and S_0) | wdef (or M5 at k = 2n-2); M2: d -> A_0 = B_0 (k >= 1); M4: d -> S < 2^n |
| sub k   | e = A_0; A = A - e*B (no borrow); then A_0 = 0 is proved and A := A/2 | facts imp, geodd; M3: e -> B < 2^wb |
| add k   | R = R + e*S (fits n bits); then S := 2S (fresh zero wire below)       | fact nand; M4: e -> S < 2^n, R + e*S < 2^n; M6 when the top wire of S is dropped |
| last (k = 2n-1) | nothing but A := A/2, S := 2S                                 | M5: A_0 = 0, wspec = 1           |
| div j   | W = window (n+1 bits): rem = W mod p (= W - p[W >= p]), f = [W < p], top wire 0. Postfact rem < p | fact: upper n bits of W < p (so W < 2p) |
| sp_const, sp_addL, sp_term | literal additions modulo 2^len (section 5)         | -                                |
| sp_q    | f = [A >= K_q], A = A - f*K_q. Postfact A < max(K_q, B - K_q)         | fact A < B                       |

Here wa = min(n, max(1, 2n-1-k)), wb = min(wa, max(1, floor((2n-k)/2))). The circuit uses the truncated widths and
omits the gates of P1-P5 of the generator docstring; the ROBDD proof shows that under the precondition the truncated
circuit equals the full-width exact contract.

## 3. Composition (machine bookkeeping + a two-line induction)

**Lemma 2.** Let x be an input for which all mathematical assumptions hold along the run of P_std (the composition of
the contracts, section 4). Then O1 and O2 hold for the forward part of P and after segment i the wires hold exactly
the state of P_std after step i (registers on the wires tracked by the walker, all other wires 0 or garbage wires of
P_std: w_k, e_k, f_j). *Proof.* Induction over segments. The entry state of segment i satisfies its precondition:
Z-wires are zero and facts hold by the induction hypothesis (they were proved at the exit of an earlier segment and
their wires were not touched since: the walker deletes a fact as soon as a segment touches one of its wires), the
mathematical assumptions hold by hypothesis. By the segment proof, (i)-(iv) hold for this entry state. []

**Lemma 3 (clean mode).** After the forward part and the n copy gates `cx out_i -> q[n+i]` (checked literally; the
target wires are in Z and are touched by no segment or mirror segment), the mirror segments in reverse order return
the state to the entry state of segment 0 with K on wires n..2n-1. *Proof.* Downward induction with (v): the state
before the mirror of segment i is the exit state of segment i (on all wires but the copy), whose entry state satisfied
the precondition. The walker checks that the mirror slices tile the inverse half exactly. []

## 4. P_std and the mathematical assumptions (paper; tested by execution)

P_std is Kaliski's algorithm with 2n steps on integers (u, v, r, s), started at (p, x, 0, 1). Step k:

    if u even:            u = u/2,      s = 2s          (w_k = 0)
    elif v even:          v = v/2,      r = 2r          (w_k = 1)
    elif u > v:           u = (u-v)/2,  r = r+s, s = 2s (w_k = 0)
    else (v >= u):        v = (v-u)/2,  s = s+r, r = 2r (w_k = 1)

i.e. w_k = [u odd] and ([v even] or [u <= v]) says "the step operates on v". The circuit keeps the registers in the
frame of w: (A, B, R, S) = (u, v, r, s) if w = 0 and (v, u, s, r) if w = 1; one swap on d = w_{k-1} xor w_k changes
the frame, after which every step reads: e = A_0; A = (A - eB)/2; R = R + eS; S = 2S. This is literally the chain of
contracts of section 2; `pstd_test` executes this chain.

Invariants of the true registers (each is checked against the four cases above):

* (I0) u >= 1, v >= 0, gcd(u, v) = gcd(p, x) = 1; if v = 0 then u = 1.
* (I1) u*s + v*r = p.
* (I2) x*r = -u*2^k and x*s = v*2^k (mod p) after k steps.
* (I3) r < 2^k and s <= 2^k after k steps (r_0 = 0 < 1; r+s < 2^k + 2^k; 2r < 2^{k+1}).
* (I4) with h_u, h_v the numbers of steps that operated on u, v (h_u + h_v = k): u < 2^{n-h_u}, v < 2^{n-h_v}
  (every step at least halves its register).
* (I5) while v >= 1 (before the step): h_u, h_v <= n-1 (from I4 and u, v >= 1), and r, s <= p < 2^n (from I1 with
  u, v >= 1). v becomes 0 only by the step (1,1) -> (1,0), which operates on v; afterwards (u,v) = (1,0), every step
  operates on v (w = 1, e = 0) and only doubles r.

The assumptions used by the contracts:

* **M1 (width)**: before step k, u, v < 2^{2n-1-k} (and < 2^n). While v >= 1: u < 2^{n-h_u} = 2^{n-k+h_v} <=
  2^{2n-1-k} by I5, same for v. If v = 0: (u,v) = (1,0) fits one bit.
* **M2 (P3)**: for k >= 1, d_k = 1 implies A_0 = B_0 (= 1). The register B not operated in step k-1 is odd (if one
  register is even it is the operated one; both even contradicts gcd = 1; v = 0 forces u = 1 odd with v operated).
  If A is even it is operated again (d = 0). So d = 1 implies A odd.
* **M3 (short v)**: if e = 1 after the swap of step k, then B < 2^{floor((2n-k)/2)}. Then both are odd and A >= B;
  b = bitlen(B) satisfies b <= n - h_B and b <= bitlen(A) <= n - h_A, so 2b <= 2n - k.
* **M4 (r, s fit n bits)**: (rs) d_k = 1 implies S < 2^n: if v >= 1 by I5; if v = 0 then k >= 1 and
  w_{k-1} = w_k = 1, so d = 0. (sfit) e = 1 implies both odd, v >= 1, S <= p. (rfit) R + eS < 2^n: for e = 0, R <= p
  (if v = 0 the operated register is v, so R = s = p by I1); for e = 1 the sum is the new r (or s) of a state with
  u, v' >= 1 resp. u >= 1, so it is <= p by I1.
* **M5 (the last two steps)**: v = 0 after at most 2n-1 steps (h_u, h_v <= n-1 before the step (1,1) -> (1,0), so
  that step has index <= 2n-2). Hence at step 2n-2 the state is (1,0) or (1,1): w_{2n-2} = 1; at step 2n-1 the state
  is (1,0): A = v = 0 is even, w = 1, d = 0, e = 0.
* **M6 (top bit of r)**: when the doubled register reaches 2n+1 wires its top wire is dropped; for k < 2n-1 that wire
  was never touched (the walker knows it is zero); at the very end it is bit 2n of r, zero by I3 (r < 2^{2n}).

**Result.** After 2n steps (u, v) = (1, 0) (M5) and by I2 x*r = -2^{2n} (mod p), i.e. r = K (mod p), 0 <= r < 2^{2n}.

`pstd_test` executes P_std bit-sliced and checks every assumption, every fact and the result:
all odd primes p < 2^11 and ALL x in [1, p) (lab run), sampled x at n = 256. This is a test of the paper part, not a proof.

## 5. Reductions

**Generic (long division; machine + two lines).** R2 = r with one fresh zero wire on top (2n+1 wires). For
j = n..0 the window W = R2[j..j+n] satisfies "upper n bits < p" (first window: tautology from p > 2^{n-1} and the zero
top wire, checked by ROBDD; later: the postfact of the previous `div`), so W < 2p and the contract gives
rem = W mod p in the low n wires, zero on the top wire. By induction (Horner) after the step j the wires R2[j..j+n-1]
hold floor(r / 2^j) mod p; at j = 0: out = r mod p = K.

**Special (folding; paper lemma with machine-checked numeric side conditions + ROBDD segments).**
*Fold lemma.* Let V < Bd be the value of the current register, m = n, L = V mod 2^m, H = floor(V / 2^m) <= Hmax =
floor((Bd-1)/2^m). Let A (len bits) be computed by the literal modular operations
A = K; A += L; A[e_k..] +-= H[0..len-e_k-1] (mod 2^{len-e_k}) for the terms (s_k, e_k). If
(a) sum s_k 2^{e_k} = 2^m (mod p), (b) K = 0 (mod p), (c) K >= sum_{s_k<0} Hmax 2^{e_k},
(d) Bn := K + 2^m + sum_{s_k>0} Hmax 2^{e_k} <= 2^len, then A = X := K + L + sum s_k H 2^{e_k} exactly,
0 <= A < Bn and A = V (mod p). *Proof.* Each literal operation is addition modulo 2^len of the corresponding term
(truncating H to len-e_k bits does not change (H << e_k) mod 2^len), so A = X mod 2^len; 0 <= X < Bn <= 2^len by
(c), (d); X = L + H 2^m = V (mod p) by (a), (b). [] The walker checks (a)-(d) with exact integers for every round
and that the register of H has bitlen(Hmax) wires; the literal operations are ROBDD-proved segments.
The quotient search `sp_q` (K_q = p 2^j, checked: multiple of p) is ROBDD-proved with bounds as carried facts; the
walker requires the final bound to be <= p, so out = A < p and A = r (mod p): out = K.

## 6. Primality

M2 and I0 need gcd(p, x) = 1 for all x in [1, p), i.e. p prime. The checker does not test primality;
secp256k1 and P-256 are standard primes. For a composite p the theorem holds for the x coprime to p on which the
assumptions hold (not claimed).

## 7. Trusted base

* the gadget lemma checker (60 lines of exact arithmetic) and Lemma 1;
* the strict lifter (`compose.lift`);
* the ROBDD package `qoptct.verify.bdd` (canonicity: equal functions have equal edges) and `compose.bv`;
* the walker `compose.kal` (contracts, Z and fact bookkeeping, coverage checks) and Lemmas 2, 3;
* the paper proofs of sections 4, 5 (M1..M6, I0..I5, fold lemma, Horner);
* NOT trusted: the generator `qoptct.crypto.kaliski`, its docstring, the certificate, the OpenQASM writer.
