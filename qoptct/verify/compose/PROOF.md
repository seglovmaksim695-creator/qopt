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

The proof has three layers. (1) Sections 1-3: machine proofs performed by the checker on the concrete text, at full
width, for all inputs at once (ROBDD). (2) Section 4: the properties M1..M6 of Kaliski's algorithm as integer
arithmetic, used as preconditions in layer 1; they are proved on paper AND by the SMT solver Z3 for every n, p, x
(`compose.lemmas_z3`, 297 implications), and tested by execution. (3) Short paper arguments that glue the layers:
Lemmas 1-3, the Horner argument and the fold lemma of section 5. Section 7 lists exactly what is trusted.

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
contract value (wires outside the contract: unchanged, ancillas: 0), (iv) the postfacts and the `proved` properties
of the contract (e.g. "no borrow", which turns a subtraction modulo 2^w into the exact integer subtraction), and
that the precondition is satisfiable (no vacuous proof).

Some preconditions are used as *parametrisations* instead of side conditions (argument `subst` of `Chain.seg`): e.g.
for the fact "A_0 -> B_0" the entry value of the wire B_0 is taken as `B_0 or A_0`. The prover checks, by ROBDD,
that the substituted values satisfy the precondition identically and that every entry state satisfying the
preconditions is a fixed point of the substitution; therefore "for all values of the variables" over the substituted
entry state is exactly "for all entry states satisfying these preconditions". (This only makes the circuit and the
contract equal as ROBDDs, so that no conjunction with the precondition has to be built per wire.) In clean mode also (v): the mirror slice of the inverse half (located by counting
ops; erase-with-l mirrors to and + cx), started from the exit state, satisfies (i), (ii) and returns every wire to
its entry value.

The contracts (frame registers A, B: n wires; R: n wires; S: n+k wires, at most 2n; w_prev; all values integers):

| segment | contract (exact integer semantics)                                   | precondition                     |
|---------|-----------------------------------------------------------------------|----------------------------------|
| init    | A = p, S = 1 (B = x, R = 0)                                           | -                                |
| cmp k   | w = wspec(A,B,w_prev) = [u odd] and ([v even] or [u <= v]), (u,v) = (B,A) if w_prev else (A,B) | M1: A, B < 2^wa   |
| swap k  | d = w_prev xor w; if d: A<->B, R<->S. Postfacts: A_0 -> B_0; A_0 -> B <= A; not(R_0 and S_0) | wdef (or M5 at k = 2n-2); M2: d -> A_0 = B_0 (k >= 1); M4: d -> S < 2^n |
| sub k   | e = A_0; A = A - e*B on wa bits, "no borrow" proved (so it is the integer subtraction); then A_0 = 0 is proved and A := A/2 | facts imp, geodd; M3: e -> B < 2^wb |
| add k   | R = R + e*S (fits n bits); then S := 2S (fresh zero wire below)       | fact nand; M4: e -> S < 2^n, R + e*S < 2^n; M6 when the top wire of S is dropped |
| last (k = 2n-1) | nothing but A := A/2, S := 2S                                 | M5: A_0 = 0, wspec = 1           |
| div j   | W = window (n+1 bits): f = bit n of (W - p) mod 2^(n+1), rem = (f ? W : W - p) on n bits, top wire 0; proved: f = [W < p] (so rem = W mod p). Postfact rem < p | fact: upper n bits of W < p (so W < 2p) |
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

### 4.1 Machine check of M1..M6 by Z3 (`python -m qoptct.verify.compose.lemmas_z3`)

The paper proofs above are replaced by an induction over k whose every implication is discharged by Z3, for all
integers n >= 2, all odd p with 2^(n-1) < p < 2^n and all 1 <= x < p at once. State before step k, in the frame
of wp = w_{k-1}: integers A, B, R, S, with ghost values a, b (exponents of the bounds) and cR, cS (quotients).
With (u, v, r, s) = (B, A, S, R) if wp else (A, B, R, S) and sg = -1 if wp else 1:

    Inv:  0 <= k <= 2n;  A, B, R, S >= 0;  A*S + B*R = p;  u >= 1;  r < 2^k;  s <= 2^k;
          x*R + sg*A*2^k = p*cR;   x*S - sg*B*2^k = p*cS;   k >= 1 -> B odd;   k = 0 -> not wp;
          either  v >= 1, v < p, a + b = 2n - k, 0 <= a, b <= n, A < 2^a, B < 2^b
          or      v = 0, u = 1, wp, k >= 1.

Proved (each as an unsatisfiable query; 2 frame bits x 4 branches of Kaliski's step x {v >= 1, v = 0} cases):

* INIT: Inv holds for (A, B, R, S) = (p, x, 0, 1), k = 0 (a = b = n, cR = 1, cS = 0).
* STEP: Inv(k), k <= 2n-1 imply, for the frame step exactly as in the contracts (w = wspec, d = wp xor w, swap,
  e = A_0, A = (A - eB)/2, R = R + eS, S = 2S): M1 (A, B < 2^wa), M2 (d -> A, B odd), M3 (e -> B < 2^wb),
  M4 (d -> R, S < 2^n; e -> S < 2^n; R + eS < 2^n), M5 (k = 2n-2 -> w; k = 2n-1 -> w, A even, d = 0, e = 0),
  M6 (2S < 2^(2n)), exactness of the step (e -> B odd and B <= A; A - eB even and >= 0), and Inv(k+1).
* FINAL: Inv(2n) implies u = 1, v = 0, 0 <= r < 2^(2n), x*r + 2^(2n) = p*c, i.e. r = K (mod p).

Assumed in these queries (all true): pow2(0) = 1, pow2(i+1) = 2 pow2(i), pow2 monotone (pow2(i) is 2^i; "a register
fits w bits" is "value < pow2(w)", with wa, wb computed by the same formulas as in the walker); p odd,
pow2(n-1) < p < pow2(n); one instance of the primality of p (if B(R+S) = p and 1 <= B < p then B = 1; it is needed
exactly in the step (u, u) -> (u, 0), where it gives u = 1). Three one-line lemmas of integer arithmetic
(y = z -> y*t = z*t) are proved by Z3 separately and used as instances, so that the main queries are linear over
the monomials. Guards against a vacuous proof: Z3 evaluates the text of Inv on all states of all real runs for
p = 5..31, 61 (pow2 as a table): it holds; 10 negative controls (hypotheses alone, slightly stronger claims, the tie
step without primality) are NOT unsat.

The reductions of section 5 are covered by four more queries (plain integer arithmetic): DIV (one step of the long
division keeps H = p*q + rem, 0 <= rem < p: the Horner induction), FOLD (the fold lemma, with the positive and
negative terms summed into Cpos, Cneg), QUOT (one step of the quotient search) and an auxiliary product bound.

What Z3 does not see: that the SMT text of the step equals the contracts of `compose.kal` (it is a transcription:
compare `lemmas_z3._case` with `Chain.run`; `pstd_test` executes the contracts themselves and checks the same
assertions), and that the literal additions modulo 2^len of a folding round compose to X mod 2^len.

### 4.2 Tests by execution

`pstd_test` executes P_std (the contracts themselves) bit-sliced and checks every assumption, every fact and the
result: all 307 odd primes p < 2^11 with ALL x in [1, p) (577 728 runs), and 4096 structured and random x for each
n = 256 prime. With `lifted=` it also runs the real circuit on the same inputs and compares every wire with P_std.

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

M2 and I0 need gcd(p, x) = 1 for all x in [1, p), i.e. p prime. Since the v0.6 audit fixes the checker TESTS this first
(`compose.primality`, result key `p_primality`):
* a composite p -> status ERROR (never PROVED; found by the v0.6 audit: p = 1021*1031 used to get PROVED);
* the named primes secp256k1, P-256, Curve25519, P-384, P-521 -> 'known prime';
* p < 2^64 -> deterministic Miller-Rabin (bases 2..37 are exact below 3.18e23), certain;
* any larger p (e.g. rand256) -> Miller-Rabin with 12 fixed + 40 pseudo-random bases, a PROBABILISTIC test; PROVED
  then carries the note 'PROVED under the condition that p is prime (only a probabilistic primality test was run)'.
The theorem is for x in [1, p) and prime p; x = 0 and x >= p are not covered.

## 7. Trusted base

* the gadget lemma checker (60 lines of exact arithmetic) and Lemma 1;
* the strict lifter (`compose.lift`);
* the ROBDD package `qoptct.verify.bdd` (canonicity: equal functions have equal edges) and `compose.bv`;
* the walker `compose.kal` (contracts, Z and fact bookkeeping, coverage checks) and Lemmas 2, 3;
* for M1..M6 and the reductions: Z3 (the 297 queries of `compose.lemmas_z3`), the three axioms of pow2, the primality of p (tested, see section 6), and the
  transcription of the frame step into SMT; the paper proofs of section 4 are then redundant;
* the glue arguments: Lemmas 1-3, and in section 5 the composition of the Z3-proved steps (Horner; fold rounds,
  whose numeric side conditions are checked by the walker with exact integers);
* NOT trusted: the generator `qoptct.crypto.kaliski`, its docstring, the certificate, the OpenQASM writer.

---

# Part II. The exact elliptic-curve point addition P -> P + Q (`compose.ecadd`, `compose.ecmath`)

Circuit: `qoptct.crypto.ecpoint.ec_point_add(p, a, b, Q)` (default construction, `legacy=False`): wires 0..n-1 hold
x, n..2n-1 hold y, the point at infinity O is encoded (0, 0); Q is a classical point. Checker:
`verify_ecadd(qasm, p, a, b, Q)` / `qoptct compose ecverify`.

**Theorem (status PROVED).** Let p >= 5 be prime, E: y^2 = x^3 + a x + b non-singular over F_p with (H1) no point
of order 2 (x^3 + a x + b has no root in F_p), Q = (xq, yq) in E with (H2) x(2Q) != xq (i.e. 2Q != +-Q; with H1:
the order of Q is not 1, 2, 3, 4), and let the certificate name a point D in E with (H3) x_D != xq, D != -2Q,
x(D + Q) != xq. Then for EVERY P in E(F_p) and for P = O, and for every string of measurement outcomes, the
circuit maps |P>|0...0> to 2^(-m/2) |ec_add(P, Q)>|0...0>: every AND target was |0>, every measurement-based
erasure is exact, no residual phase, every ancilla returns to |0>. `ec_add` is the textbook affine addition
(`compose.ecadd.ec_add`); it is the specification. Inputs that are not points of the curve are outside the domain.
H1, H2, H3, the primality of p (section 6) and the non-singularity are CHECKED by the verifier with exact integers
(`check_parameters`; H1 by gcd(x^p - x, x^3 + a x + b) = 1); a violation gives status ERROR.

## 8. Segments of the point addition

The lifted circuit (section 1) is walked by `ecadd.ECChain`. Besides the ROBDD segments of section 2 there are:

**Algebra segments** (products, squares, fold reductions). The slice of the real ops between two marks is handed to
`qoptct.verify.algebra.prove_ops` with a word-level spec built by the walker: input words (free), output words,
garbage wires, claims (`out == x*y`; `r + p*(...) == value`, `r <= p-1`). The walker first checks that every wire
mentioned by the slice or named as output / garbage, other than the input words, is in Z (known |0>).

*Lemma 4 (interface).* If these wires are |0> and the prover returns PROVED, then on every entry state the slice
satisfies every gadget precondition, the claims hold between the entry values of the input words and the exit
values of the output words, the kept input wires are unchanged, garbage wires hold some classical function of the
inputs, and every other wire of the slice ends at |0>. *This is exactly the statement of the algebra prover*
(its docstring; it does not look at wires outside the slice, and wires outside the slice are not touched). The
walker then removes the output and garbage wires from Z and forgets every fact about them. In the ROBDD world the
outputs get NO value: every later ROBDD segment is proved for all values of these wires; only the proved bound
`r < p` is carried (fact `lt`). A claim `r + p*(integer expression) == V`, `r <= p-1` says r = V mod p whatever the
expression is, so the spec needs no trust in the generator's formulas.

**Mirror segments.** Every block "compute F / apply A / uncompute G" ends with a mirror segment.

*Lemma 5 (mirror).* Let F = ops[s:e), A = ops[e:e'), G = ops[e':e'') where G is the op-by-op inverse of F in reverse
order (x, cx: themselves; and(a,b -> t): erase of t against a*b; erase of t against a*b xor l: and(a,b -> t)
followed by cx l -> t). Let sigma be a classical state on which F runs with every gadget precondition true, ending
in sigma'; let A run validly from sigma' to sigma'' with sigma''(w) = sigma'(w) for every wire w that F mentions
in ANY role (control, target, operand, correction). Then G runs from sigma'' with every gadget precondition true
and ends in tau with tau(w) = sigma(w) for the wires mentioned by F and tau(w) = sigma''(w) elsewhere.
*Proof.* Induction over the ops of F from the last to the first; the ops of F and G read and write only wires
mentioned by F. If F = F0 ; o and rho is the state before o in the run of F, the state before the image of o in G
agrees with the state after o on these wires. o = x or cx: the same gate undoes it. o = and(a, b -> t): after o,
t = a*b, so the erase condition of its image holds and t returns to 0 = rho(t). o = erase(t; a, b, l): rho(t) =
rho(a)rho(b) xor rho(l) and t = 0 after o, so the target of the image AND is 0 and and + cx restore rho(t). []
The walker checks (i) G is this inverse, op by op (integer comparison of the lifted ops); (ii) no contract of an
apply segment has a post on a wire mentioned by F -- and each apply segment is a ROBDD segment, in which every
touched wire without a post is PROVED to end with its entry value. This is the condition whose absence was hole 1
of the v0.6 audit. Then Z and the facts are restored for the wires of F: the wires of F that were |0> before F are
|0> again; facts of the state before F that do not speak about wires written by A hold again.
With Lemma 1 (gadgets) the mirrored erasures have no phase.

**Zero-length segments** carry a named assumption (precondition) on the current state.

**The chain W_ec** (contracts; R + k etc. are mod p, registers < p are carried as `lt` facts):

| step (marks) | contract | proved by |
|---|---|---|
| `w_pre` | flags f_O, f_Q, f_-Q, f_-2Q := [P = O, Q, -Q, -2Q]; (X, Y) := D if a flag is set | ROBDD, from fact `ecdom` |
| `addc` x2 | X := X - xq, Y := Y - yq | ROBDD |
| division 1 (compute) | `cm_*`: A := K + sum s 2^e X (shift-and-add, K = 0 mod p); reduction -> xs = 2^(2n) X mod p; **G1** xs != 0; Kaliski core -> r_ext (section 2-4); reduction -> K; product `mul` + reduction -> T = Y*K mod p | ROBDD; algebra; ROBDD + M1..M6; algebra |
| apply | Lr := T (Lr was 0) | ROBDD |
| mirror | Lemma 5 | syntactic |
| product (Lr, X); apply Y := Y + T; mirror | | algebra; ROBDD; syntactic |
| square Lr; apply X := -X, X := X + T; mirror | | algebra; ROBDD; syntactic |
| `addc` | X := X - 3xq | ROBDD |
| product (Lr, X); apply Y := Y + T; mirror | | |
| `neg` | X := -X | ROBDD |
| division 2; apply: **G2** T = Lr, Lr := 0; mirror | | |
| `neg`, `addc` x2 | X := -X + xq, Y := Y - yq | ROBDD |
| `w_post` | **G3** a flag is set -> (X, Y) = D + Q; **G4** no flag -> (X, Y) is in E, not Q, -Q, 2Q, O. Images of the flagged inputs (Q, 2Q, O, -Q), all four flags erased to 0 | ROBDD |
| end | every wire except X, Y is in Z | bookkeeping |

For a prime without a one-fold reduction (tests at small n) the reductions are the signed-term rounds of section 5
(ROBDD); the bound of a product, (p-1)^2 + 1, is then admitted from `out == x*y` and x, y <= p-1.
*Shift-and-add lemma* (as the fold lemma of section 5): the literal additions of X[: L-e] into A[e:] modulo
2^(L-e) are A := A + s 2^e X modulo 2^L; with K >= sum_{s<0} xmax 2^e and K + sum_{s>0} xmax 2^e < 2^L (checked) the
integer K + sum s 2^e X lies in [0, 2^L), so A holds it exactly; it is = 2^(2n) X mod p because K = 0 and
sum s 2^e = 2^(2n) (mod p) (checked).

## 9. What the walk proves, and the assumptions

Machine-proved by the walk, at the real width, for all entry values of every segment: **the circuit refines W_ec**:
on every input state on which the assumptions below hold along the run of W_ec, every gadget precondition holds and
the wires follow W_ec; at the end only X, Y are non-zero. The assumptions (reported under `assumptions`):

* `D-input` -- the domain: (x, y) is O = (0, 0) or a point of E. Used through its consequences `ecdom`: x, y < p;
  y = 0 -> x = 0 (H1); x = xq and y != 0 -> y = +-yq; x = x(2Q) and y != 0 -> y = +-y(2Q) (y^2 is determined by x).
* `G1` (twice) the Kaliski input is not 0; `G2` the second division returns the content of Lr; `G3`, `G4` as in
  the table. They are discharged by Theorem T below.
* `M1..M6` -- Kaliski (sections 4, 4.1; Z3 for all n, p, x in [1, p)). Their hypothesis x in [1, p) is G1 plus the
  proved bound xs < p.

## 10. Theorem T: W_ec computes ec_add(P, Q) and G1..G4 hold

Notation: f(x) = x^3 + a x + b; for P = (x, y) in E with x != xq: lambda = (y - yq)/(x - xq) and
chord(P, Q) = (x', y'), x' = lambda^2 - x - xq, y' = lambda (xq - x') - yq (= ec_add(P, Q) by definition, the two
formulas for y' agree because lambda (x - xq) = y - yq). mu = (3 xq^2 + a)/(2 yq) (yq != 0 by H1), tangent line
l(x) = mu (x - xq) + yq, 2Q = (xT, yT) with xT = mu^2 - 2 xq, yT = mu (xq - xT) - yq = -l(xT).
Polynomial identities, checked by exact integer expansion (`ecmath.check_identities`; true over Z, hence in F_p):
I1 (chord(P, Q) satisfies the curve equation, modulo y^2 = f(x)), I2 (4 yq^2 (f(x) - l(x)^2) =
4 yq^2 (x - xq)^2 (x - xT)), I3 (modulo the curve: x' = xq is equivalent to y = l(x), after multiplying by
(x - xq)^2 and 2 yq).

(a) A point of E with abscissa xq is (xq, yq) or (xq, p - yq), and these have different parities of y (p odd,
yq != 0); the same for 2Q. No point has y = 0 (H1). [-> `ecdom`]
(b) chord(P, Q) is in E (I1).
(c) If P in E, x != xq and x' = xq, then P = -2Q. *Proof.* By I3, y = l(x). Then 0 = f(x) - y^2 = f(x) - l(x)^2 =
(x - xq)^2 (x - xT) by I2, so x = xT and y = l(xT) = -yT. []
(d) If P in E, x != xq and chord(P, Q) = 2Q, then P = -2Q. *Proof.* y' = yT gives lambda (xq - xT) = mu (xq - xT);
xT != xq (H2), so lambda = mu and y = yq + lambda (x - xq) = l(x); conclude as in (c). []
(e) chord(-2Q, Q) = -Q: -2Q = (xT, l(xT)) lies on the tangent, so lambda = mu, x' = mu^2 - xT - xq = xq,
y' = -yq.

*Run of W_ec on P in E or O.*
1. `w_pre` (machine, from `ecdom`): the flags are the indicators of P = O, Q, -Q, -2Q -- four different points
   (yq != 0, H2) -- and P' = D if P is one of them, else P. So P' in E, x_P' != xq, P' != -2Q (H3, resp. (a)).
2. dx = x_P' - xq != 0, so xs = 2^(2n) dx mod p is in [1, p): **G1**. By the Kaliski contract (sections 2-5)
   r_ext = -xs^{-1} 2^(2n) = -dx^{-1} (mod p), hence Lr = -dy/dx = -lambda. Then Y = dy + Lr dx = 0;
   X = Lr^2 - dx - 3 xq = x' - xq; Y = Lr (x' - xq) = lambda (xq - x') = y' + yq; X = xq - x'.
3. xq - x' != 0: otherwise P' = -2Q by (c). **G1** (second). The second quotient is
   -(y' + yq)/(xq - x') = -lambda (xq - x')/(xq - x') = -lambda = Lr: **G2**. Finally (X, Y) = (x', y') =
   chord(P', Q).
4. If a flag is set, P' = D and chord(D, Q) = D + Q (H3: x_D != xq): **G3**. Otherwise R = chord(P, Q) is in E (b),
   R != O-code since y' != 0 (H1), x' != xq by (c) so R != +-Q, and R != 2Q by (d): **G4**.
5. `w_post` (machine, from G3, G4): the result is Q, 2Q, O, -Q for P = O, Q, -Q, -2Q -- these are ec_add(P, Q) by
   the definition of ec_add (for -2Q by (e)) -- and chord(P, Q) = ec_add(P, Q) otherwise; all flags are erased. []

Evidence beyond the written proof (finite, not a proof): `ecmath.exhaustive_T` executes an independent field-level
model of W_ec on every curve / Q / safe D / point for all primes p <= 31 (992 curves for p = 17..31 with D = 3Q,
all safe D for p <= 13: 36 thousand cases, 0.9 million points, no violation); `wec_test` executes the contract
chain itself (the code that the ROBDD proofs use) on all points of small curves and on 64 points at n = 256
(exceptional points, D, D + Q, random) and compares the real circuit on the same points.

## 11. Machine / paper / trusted, for the point addition

* **Machine, on the given text, full width, all values:** lifting and gadget lemmas; every ROBDD segment (flags,
  substitution, modular additions / negations, shift-and-add rounds, Kaliski steps, special reductions, images and
  flag erasure); every algebra sub-goal (5 products / squares and 9 fold reductions at n = 256, forward halves);
  the 5 mirror comparisons and their frame conditions; coverage of all ops; all ancillas in Z at the end; the
  numeric side conditions (H1-H3, primality, fold / shift-and-add conditions).
* **Machine, in general:** M1..M6 (Z3); I1..I3 (exact polynomial expansion).
* **Paper:** Lemmas 1-3 (Part I), Lemma 4 (reading of the algebra prover's statement), Lemma 5, the shift-and-add
  lemma, the case analysis of Theorem T (uses only I1..I3 and the field axioms; no appeal to the group law), and
  that the contracts of `ecadd.py` are the program W_ec of section 10 (tested by `wec_test`).
* **Trusted in addition to section 7:** `qoptct.verify.algebra` (word-level prover; audited in v0.6, not checked by
  an independent tool), the client `compose.ecadd`, the 40-line polynomial arithmetic of `compose.ecmath`.
* **Not trusted:** the generators (`qoptct.crypto.ecpoint`, `_ecp`, `kaliski`, `mulsq`, `fold`), the certificate
  (marks, D), the OpenQASM writer.
* **Not covered:** the windowed addition (quantum table lookup of Q), the legacy construction (`legacy=True`),
  inputs outside the curve, curves with a point of order 2, Q of order <= 4.
