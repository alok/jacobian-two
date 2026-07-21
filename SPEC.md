# Jacobian Two: research specification

## Mission

Build a public, Lean-first research repository that does two things without
conflating them:

1. independently certify the polynomial map announced by Levent Alpöge on
   2026-07-20 UTC (2026-07-19 Pacific) as a counterexample to the complex
   Jacobian conjecture in dimension three; and
2. use the structure of that example to pursue a concrete obstruction theorem
   for the still-separate two-dimensional case.

The repository is tracked by Linear issue
[ALOK-792](https://linear.app/aloksingh/issue/ALOK-792/formalize-the-3d-counterexample-and-establish-the-jc2-research).
The Galois/frontier continuation is tracked by
[ALOK-794](https://linear.app/aloksingh/issue/ALOK-794/audit-the-galois-frontier-and-determine-the-3d-maps-nonproper-value).
The arbitrary-degree affine-coordinate continuation is tracked by
[ALOK-795](https://linear.app/aloksingh/issue/ALOK-795/prove-the-arbitrary-degree-affine-coordinate-plane-keller-theorem).
The constant-leading quadratic-coordinate continuation is tracked by
[ALOK-798](https://linear.app/aloksingh/issue/ALOK-798/prove-the-constant-leading-quadratic-coordinate-keller-theorem).
The variable-leading quadratic direct formalization is tracked by
[ALOK-799](https://linear.app/aloksingh/issue/ALOK-799).
The generic-degree-six monodromy attack is tracked by
[ALOK-801](https://linear.app/aloksingh/issue/ALOK-801/classify-and-attack-generic-degree-six-plane-keller-monodromy).

## Exact input claim

Let `F = (P, Q, R)` with

```text
P = (1 + x*y)^3*z + y^2*(1 + x*y)*(4 + 3*x*y)
Q = y + 3*x*(1 + x*y)^2*z + 3*x*y^2*(4 + 3*x*y)
R = 2*x - 3*x^2*y - x^3*z.
```

The screenshot claims:

- the formal Jacobian determinant of `F` is the constant polynomial `-2`; and
- the three distinct rational points
  `(0, 0, -1/4)`, `(1, -3/2, 13/2)`, and
  `(-1, 3/2, 13/2)` all map to `(-1/4, 0, 0)`.

These are finite algebraic identities. They must be proved exactly, not inferred
from floating-point samples or from the authority of the announcement.

## Claim-status ledger

| Claim | Required status | Evidence |
|---|---|---|
| The screenshot is transcribed correctly | `SOURCE-CHECKED` | screenshot, original post, and two independent transcriptions agree |
| `det JF = -2` identically | `LEAN-CERTIFIED` | formal partial derivatives of multivariate polynomials and a determinant identity |
| The displayed fiber has three distinct points | `LEAN-CERTIFIED` | exact rational evaluation and pairwise inequality |
| The universal-in-dimension formulation is false | `DERIVED` | the certified dimension-three counterexample suffices |
| Counterexamples exist in every dimension `n >= 3` | `DERIVED` | append identity coordinates |
| Whether the two-dimensional conjecture is false | `OPEN` | no implication from the dimension-three example is known or assumed here |
| The Galois-case theorem settles every plane Keller map | `FALSE` | constant Jacobian gives finite separable, not normal/Galois |
| Generic degree six is the first unresolved plane sheet degree | `SOURCE-CHECKED` | Campbell, Orevkov, Domrina, Żołądek, and Borisov literature audit |
| A degree-six plane Keller counterexample must have monodromy `A6` or `S6` in its transitive degree-six action | `DERIVED / EXACTLY ENUMERATED` | affine fixed-sheet inertia leaves seven groups; the refined identity `sum(e*d+delta)=5`, exact class costs, and the irreducible-branch obstruction eliminate the other five |
| A one-dicritical degree-six counterexample can only have `A6` with `(e,d)=(3,1)` or `S6` with `(2,1)`, and its branch normalization is noninjective | `DERIVED / EXACTLY ENUMERATED` | finite-normalization deck rigidity, local multiplicity, purity, Lin--Zaidenberg topology, and exact subgroup witnesses |
| In the one-dicritical `A6` passport, the two jump units occur at one unique local-degree-five point, and its branch has Puiseux pair `(2,5)` | `DERIVED / HOSTILE-AUDITED` | one-chain source topology, Hirzebruch--Jung cyclic lift, quadratic-jet characters, exact local intersection parity, and the surviving `A5 -> A6` permutation fixture |
| The forced `A6` degree-five cusp block is locally analytically realizable and has a unique Orevkov case-(b) normal form up to the classified equivalence | `SOURCE-CLASSIFIED / EXACTLY CHECKED` | Orevkov 2026 Theorem 2 and Section 5; `a6_cusp_germ.py` checks `(d1,d2,N,n)=(2,5,5,3)`, the square Jacobian, full three-component cusp pullback, and componentwise degree-one normalizations |
| In the one-dicritical `S6` passport, the finite fiber has one of nine exact rows, every jump block has local group `S_(2+kappa)`, and every jump is `T(2,3)`, `T(3,4)`, or `T(4,5)` | `DERIVED / EXACTLY ENUMERATED` | Orevkov's three-unit budget, constant multiplicity, transposition normal generation, simple ramification, and hostile analytic germs |
| In both surviving one-dicritical passports, Orevkov's minimal `L_C` is empty and the affine finite normalization `W` is smooth finite flat of rank six, with `W-D=A2`, `D=A1`, `Pic(W)=Z[D]`, and `K_W=(e-1)D` | `DERIVED / INDEPENDENTLY HOSTILE-AUDITED` | Orevkov's relative-minimality statement, the index-two invariant-order contradiction, the index-three cyclic-cover parity theorem at local degrees three and five, miracle flatness, divisor localization, and finite-map Riemann--Hurwitz |
| On a jointly minimal one-dicritical source resolution, the dicritical is a type-3 `(-1)`-leaf of label `e`, its unique neighbor is type 2 of label `e-1`, and the path to the original infinity line contains `1--0--(-1)` | `DERIVED / INDEPENDENTLY HOSTILE-AUDITED` | Orevkov's finite-chain lemma, simultaneous Castelnuovo contraction for the source blowdown and resolved map, logarithmic Keller adjunction, and Borisov's canonical-label tree theorem |
| On that minimal leaf edge, determinant labels satisfy `d_E=d_A-1<0` with forced parity, but the canonical-plus-determinant constraints still admit an explicit infinite family | `DERIVED / EXACT HOSTILE FAMILY` | Borisov's determinant/content recurrences, Hodge index, and exact intersection-matrix replay in `canonical_leaf_graph.py` |
| If the one-dicritical branch has exactly two characteristic pairs at infinity and Orevkov's condition `(*)` holds, neither degree-six `A6` nor `S6` passport is possible | `DERIVED / INDEPENDENTLY HOSTILE-AUDITED / CONDITIONAL` | Orevkov's equations (10)--(15): the product identity gives `N>=8` for `A6`; the sole `S6` near-miss forces the impossible integral edge equation `3*Q_tilde=14`; `two_pair_infinity.py` replays the arithmetic and Belyi fixture |
| Projection formula and Hodge index do not by themselves force Orevkov's condition `(*)` in degree six | `DERIVED MATRIX THEOREM / EXACT HOSTILE BLOWUP TREES` | for split preimages the scalar relation becomes `Qm=qn`, with `Q-(q/N)nn^t` negative semidefinite when `q>0`; unimodular trees realize the equations with `A6` `3+3` and `S6` `2+2+2` inertia, and a stronger labeled `A6` fixture also retains local canonical pullback, valency, monomial, and primitive-monodromy data |
| Degree six, finite flatness over `A2`, geometric monodromy `S6`, and exact global boundary transport do not imply condition `(*)` without the Keller hypothesis | `EXACT FINITE HOSTILE COVER` | `(x,z)->(x^2+z,x^3-x^4z)` is free rank six, extends finitely from `P1 x P1` to `P2`, has primitive polynomial monodromy with a transposition, and pulls the sole boundary line back to two noncontracted components; its nonconstant Jacobian keeps the Keller frontier open |
| The one-boundary surface packages `K_W=D`, `K_W=2D` and the saturated two-boundary package `K_W=D2+2D3` are internally consistent and therefore cannot alone eliminate the surviving passports | `DERIVED / EXPLICIT HOSTILE MODELS` | two Hirzebruch-surface complements and one triple-blowup complement, with exact Picard, canonical, intersection, and Nakai--Moishezon calculations |
| In all three smooth finite-flat packages, the inverse different is a nontrivial source line, the algebra is not globally Frobenius or monogenic, and the trace matrix factors as `T=Phi^* H Phi` with exact determinants `b`, `b^2`, or `b2*b3^2` and normalization-module cokernels | `DERIVED / FINITE-DUALITY AUDITED` | finite Grothendieck duality, different/discriminant norm, divisor-step lattices, Quillen--Suslin, and the Tate determinant map |
| The quadratic trace-lattice, corank, trace-unit, symmetric-determinant, and theta conditions alone do not eliminate the one-dicritical `A6` or `S6` survivors | `EXACT HOSTILE MATRICES` | explicit Bezout matrices realize the `A6` cusp/collision table and a full three-cusp `S6` passport; `trace_hostile_matrices.py` checks determinants, kernel vectors, ranks, and norm-six vectors; no multiplication law is asserted |
| The explicit three-cusp hostile `S6` trace curve cannot carry connected six-sheet transposition monodromy | `DERIVED TOPOLOGICAL EXCLUSION / EXHAUSTIVELY CHECKED` | its quartic projection gives four meridian generators, whereas transitive transpositions on six sheets need at least five; the exact Zariski--van Kamp replay checks all `15^4` assignments, finds `735` representations and zero transitive images; this excludes only that curve, not the `S6` passport |
| A curve in the degree-minimal singular-one-pair width-five `S6` slot cannot carry connected six-sheet monodromy | `EXACT GEOMETRY / COMPUTER-ASSISTED TOPOLOGICAL EXCLUSION` | `(t^5+t^4,t^7+t^5)` has one `T(4,5)` cusp, six nodes, a `(2,7)` infinity pair, and delta `6+6+3=15`; Sage gives complement group `Z`, and the dependency-free `15^5` transposition replay finds only 15 order-two images and none transitive; this excludes that curve, not every width-five curve |
| Under the additional exactly-one-genuine-pair hypothesis, an `S6` branch cannot have normalization target degree at most eleven; the first infinity-group survivor is affine degree pair `(5,12)` | `DERIVED TOPOLOGICAL EXCLUSION / EXHAUSTIVELY CHECKED / CONDITIONAL` | the affine link group `pi1(T(m,n))` surjects onto the global complement; width gives `m>=5`, and exact centerless torus-group censuses find no transposition-meridian `S6` quotient for all ten coprime singular pairs through `n=11`; `(5,12)` has exactly `720` such quotients and an explicit star-transposition witness, but no branch curve or Keller cover is constructed |
| Under the four standing polynomial-normalization, one-genuine-pair, cusp-only, and smooth-collision hypotheses, a one-dicritical `A6` branch must have normalization-collision delta at least ten; at the coarse link stage equality leaves only affine degrees `(4,9)` and projective infinity pair `(5,9)` | `EXACT ALGEBRA / COMPUTER-ASSISTED TOPOLOGY / CONDITIONAL` | genus and exact link censuses reduce the low-delta candidates; complete family audits exclude every candidate through `Delta=7`, and exact torus-link censuses exclude `Delta=8,9`; at `Delta=10`, only `T(4,9)` has a single-three-cycle generating meridian, with exactly `720` qualifying pairs; that coarse scan alone neither derives the four hypotheses nor constructs a curve, cover, or Keller map, and it does not prove `JC(2)` |
| The clean locus of the conditional `(4,9)`, delta-ten `A6` family is excluded | `EXACT FAMILY ALGEBRA / COMPUTER-ASSISTED TOPOLOGY / CONDITIONAL` | every normalization has the five-parameter form `(t^2+k*t^3+t^4, a*t^5+b*t^6+c*t^7+d*t^8+t^9)` up to the residual involution; an exact member has one `T(2,5)` cusp, ten nodes, `T(5,9)` at infinity, cyclic complement, and only 40 diagonal `C3` assignments in the exhaustive `40^4` replay; the clean parameter locus is nonempty Zariski open and connected, and family-wide propagation uses proper projective Whitney--Thom triviality |
| The two dominant conditional delta-ten degeneration divisors, contact-two plus eight nodes and ordinary-triple plus seven nodes, are excluded | `EXACT INCIDENCE ALGEBRA / COMPUTER-ASSISTED TOPOLOGY / CONDITIONAL` | localized minor ideals prove that each valid incidence is one irreducible four-dimensional dominant component; exact rational representatives have cyclic complements and only 40 diagonal `C3` assignments in their full `40^4` raw-presentation replays; checked Sage regenerates both curves, singular schemes, and presentations, while family-wide propagation uses proper projective Whitney--Thom triviality |
| All six displayed generic or dominant components of the expected codimension-two conditional delta-ten profiles are excluded | `EXACT INCIDENCE ALGEBRA / COMPUTER-ASSISTED TOPOLOGY / CONDITIONAL / BOUNDARIES OPEN` | exact rank and saturation certificates treat `C3+7N`, `C2^2+6N`, `T112+6N`, `C2+T111+5N`, `Q0+4N`, and `T111^2+4N`; six rational representatives have cyclic complements and only 40 diagonal `C3` assignments in exhaustive `40^4` replays; topology transport for `T112` and the mixed chart is proved over their smooth labeled incidences by finite-etale labeling, relative embedded resolution, and proper Whitney--Thom isotopy; the three exposed residual factors cannot support a threefold, and the immersed `P`-critical `T112` and mixed boundaries have dimension at most two; a generated 22-row split allocation ledger has exact clean witnesses and full `k=2` to `k=-2` transport, but most split rank-drop subloci, denominator/overlap charts, component intersections, and lower-dimensional residual taxonomy remain open, so this neither completes the delta-ten audit nor proves `JC(2)` |
| The dense clean nonsplit `C4+6N` surface, the first expected codimension-three delta-ten profile audited, is excluded | `EXACT INCIDENCE ALGEBRA / COMPUTER-ASSISTED TOPOLOGY / CONDITIONAL / BOUNDARIES OPEN` | `H=H'=H''=H'''=0` has determinant `24*(k+2s)^4*(s^2+ks+1)^5*F`; Sage saturation makes the valid residual compatibility scheme zero-dimensional of length ten with rank-three affine-line fibers, so the Cramer graph is the sole nonsplit surface; an exact member has a `T(2,5)` cusp, one contact four, six nodes, a cyclic complement, and zero `A6` assignments in the exhaustive `40^4` replay; on the smooth irreducible clean Cramer open, finite-etale labeling, four relative contact blowups, and proper Whitney--Thom isotopy propagate the exclusion; split, pair-denominator, cusp-pair, diagonal, residual, and deeper intersections remain open |
| The conditional delta-ten collision ledger has 35 non-node local atoms and 145 global candidate profiles | `EXACT COMBINATORICS / ELIMINATION TARGET` | integral smooth-branch contact trees give 9 two-branch, 15 triple, and 11 quadruple atoms; the expected-codimension histogram is `(1,2,6,14,29,38,32,16,6,1)`, with six expected codimension-two profiles, 38 expected endpoints, and 55 overdetermined profiles that remain candidates until saturated elimination proves their status |
| The local cusp, collision, orientation, Nielsen-class, and spin constraints do not eliminate the surviving conditional `T(4,9)` passport | `EXACT HOSTILE CENSUS` | all `720` qualifying torus pairs split into two inner `A6` orbits fused by odd relabeling; every pair has one compatible forced cusp/collision decomposition and an exact `2.A6` lift, so these data are consistency checks rather than an obstruction |
| The conditional `(3,8)`, delta-five `A6` family is exhausted, including its codimension-two residual and every valid exceptional point | `EXACT FAMILY ALGEBRA / COMPUTER-ASSISTED TOPOLOGY / CONDITIONAL` | exact primary decomposition gives four valid rational residual curves `A,B,N,E`; exact factor identities leave only `P4`, `PNT`, and the two conjugate `P32` points; the four generic and four exceptional Sage presentations each have an exhaustive `40^3` replay with exactly 40 cyclic images and no `A6`; Zariski--van Kamp extraction and Whitney--Thom propagation remain explicit computer-assisted dependencies |
| The full conditional `(3,10)`, delta-seven `A6` family is excluded, including the triple-image wall, every positive-dimensional repeated-collision stratum, and every finite endpoint | `EXACT FAMILY ALGEBRA / COMPUTER-ASSISTED TOPOLOGY / CONDITIONAL` | the four-parameter normal form, collision-septic coefficient slice, root-partition ledger, resultants, saturations, Groebner bases, incidence-component checks, Hurwitz enumerations, and `40^3` permutation replays are exact; stored Zariski--van Kamp presentations, proper Whitney--Thom propagation on connected equisingular strata, and finite-etale/Riemann-existence transport with tame inertia across arithmetic endpoint embeddings are explicit dependencies, and not every original presentation extraction is regenerated by a checked script |
| The exact specialized `A6` trace and middle lattices admit the required associative `5+1` cusp algebra and `3+3` collision algebra | `EXACT LOCAL HOSTILE ALGEBRAS` | `a6_local_multiplication.py` checks both Gram-basis changes, perfect middle forms, divisor-section factorizations, and kernel alignment; the global 56-entry cubic-tensor system is written explicitly but neither solved nor contradicted |
| The forced finite `A6` cusp and collision relations alone do not produce a `2.A6` spin obstruction | `EXACT HOSTILE LIFTS` | rational Clifford arithmetic checks the cusp five-braid and collision commutator upstairs, enumerates the generated group of order `720`, and exhibits transitive all-3-cycle product-one completions with spin signs `+1` and `-1`; no actual infinity word is asserted |
| In the saturated two-curve `S6` transposition-plus-3-cycle passport, both branch normalizations self-collide; the exact collision rows are `3+3`, `2+2+1+1` or `2+2+2`, and cross-intersections are `2+3+1` | `DERIVED / EXACT HOSTILE FIXTURE` | Nguyen's singular-component theorem, zero excess, Orevkov's local embedding, constant multiplicity six, and four exact meridians generating `S6` |
| In that saturated `S6` passport, both minimal constant chains are empty and the finite normalization is smooth finite flat, with boundary `D2 disjoint-union D3`, `Pic(W)=Z[D2]+Z[D3]`, and `K_W=D2+2D3` | `DERIVED / INDEPENDENTLY HOSTILE-AUDITED` | componentwise zero excess, the index-two and index-three cyclic endpoint obstructions, complement purity, miracle flatness, and divisor localization |
| The first residual transformed `(72,108)` Newton configuration needs at least three nonzero strict-interior coefficients, and the second needs at least four | `DERIVED / EXHAUSTIVELY CHECKED` | all `7504` first-case supports of size at most two and all `3683` second-case supports of size at most three have replayed exact certificates; the five zero-product exceptions have unit-ideal certificates |
| The universal cubic discriminant coefficient expression is `-4*Q` | `LEAN-CERTIFIED` | explicit standard formula; semantic SymPy discriminant checked independently |
| The announced map has fiber counts `3/1/0` on the three stated strata | `DERIVED` | simple projective-root/source bijection plus exact symbolic checks |
| Its image is `C^3 \ Gamma` | `DERIVED` | complete fiber stratification |
| Its nonproper-value set is exactly `V(Q)` | `DERIVED` | explicit escaping family and projective compactness argument |
| A Keller map `(P,e(x)y+f(x))` with arbitrary `P` is an automorphism | `LEAN-CERTIFIED` | actual bivariate Jacobian, degree descent, both explicit inverse charts |
| A Keller map `(P,eps*y^2+g(x)y+f(x))` with `eps != 0` and arbitrary `P` is an automorphism | `LEAN-CERTIFIED` | affine discriminant, completed-square coordinates, normal form, and both explicit inverse laws |
| A Keller map with one coordinate of `y`-degree at most two is an automorphism | `KNOWN / LEAN-CERTIFIED` | Moskowicz, Theorem 2.7; independently kernel-checked here, including the affine branch |
| The direct reduction from variable-leading quadratic `Q` to the scalar-leading theorem | `DERIVED / LEAN-CERTIFIED` | Lean certifies odd-degree descent, the UFD square shape, fraction-field centering, the exact `k/h` transport, recurrence descent in `K[F]`, specialization, gcd-normalized noncancellation, and the final unit conclusion |
| Historical novelty of the full fiber/nonproper theorem | `UNKNOWN` | derived here; same-day sources compared, priority not established |
| Historical priority and the discovery account | `ANNOUNCED` | public posts are new and not a peer-reviewed historical record |

## Milestone 1: kernel-checked screenshot certificate

Implement the announced coordinates as integer-coefficient
`MvPolynomial (Fin 3) ℤ`. Define the formal Jacobian matrix with
`MvPolynomial.pderiv`, prove its determinant equals the constant `-2`, and
base-change the pointwise certificate through `Int.castRingHom`.

Separately evaluate the map at the three rational points, prove their images
are equal, prove the inputs are pairwise distinct, and derive that the induced
map is not injective. Where useful, transport the rational certificate to
`Complex` explicitly.

Acceptance gate:

- `lake build` passes;
- no `sorry`, `admit`, or custom `axiom` occurs in the certificate;
- `#print axioms` reports only standard foundations, if any;
- an independent typed SymPy program reproduces both identities with exact
  rational arithmetic; and
- adversarial tests fail after perturbing a coefficient or a collision point.

## Milestone 2: a two-dimensional obstruction theorem

The three-dimensional map is affine in `z`. The first tractable question is
whether the same low-complexity strategy can work in two dimensions.

For a characteristic-zero field `K`, consider

```text
G(x, y) = (A(x)*y + B(x), C(x)*y + D(x))
```

with `A, B, C, D : K[X]`. Its Jacobian splits as

```text
det JG = (A'*C - A*C')*y + (B'*C - A*D').
```

The first formal target is:

> If `det JG` is a nonzero constant, then `A' = 0` and `C' = 0`.

Proof plan: the constant term gives an explicit Bézout identity for `A` and
`C`, hence they are coprime. The vanishing `y` coefficient is their vanishing
Wronskian. Mathlib's coprime-Wronskian theorem then forces both derivatives to
vanish. This rules out a direct two-dimensional analogue in which both
coordinates are affine in one source variable.

The completed extension proves that the resulting plane map is injective: a
linear combination of the outputs recovers an affine polynomial in `x` with
nonzero slope, after which either nonzero constant `y`-slope recovers `y`.

The explicit reconstruction is now completed for this ansatz.  The broader
quadratic-in-one-variable class is treated separately in Milestone 4.

## Milestone 3: Galois audit and the exact nonproper-value set

The classical Galois-case theorem has a one-way hypothesis:

> A Keller map whose induced function-field extension is Galois is a
> polynomial automorphism.

It does **not** say that every two-dimensional Keller extension is Galois.
Consequently it does not settle `JC(2)`.  The literature audit must distinguish
the field-extension degree from the ordinary total degree of the coordinate
polynomials, and it must record the accepted result that a hypothetical plane
counterexample has topological (generic-fiber) degree at least six.

For the announced three-dimensional map, write a target as `(a,b,c)` and put

```text
p(T) = c*T^3 - 2*T^2 + b*T - 2*a
r(T) = p'(T) = 3*c*T^2 - 4*T + b
Q(a,b,c) = 27*a^2*c^2 - 18*a*b*c + 16*a + b^3*c - b^2.
```

On the chart `x != 0`, the fiber coordinate `T = y + 1/x` satisfies
`p(T)=0`, `r(T)=2/x`, and the source is reconstructed by

```text
x = 2/r,
y = T - r/2,
z = 5*r^2/4 - 3*T*r/2 - c*r^3/8.
```

The research target is the following complete determination of the Jelonek
set (the nonproper-value set) of the announced map:

> A target `(a,b,c)` is a nonproper value if and only if `Q(a,b,c)=0`.

The two directions must be checked separately.

1. If a bounded target sequence has an unbounded source sequence and `T`
   stays bounded, reconstruction forces `r -> 0`; the limiting polynomial has
   a repeated root, hence `Q=0`.  If `T -> infinity`, eliminate `c` using
   `p(T)=0`; the exact identities

   ```text
   x = T / (T^2 - b*T + 3*a),
   y = b - 3*a/T
   ```

   and an inverse-power expansion for `z` show that the source instead tends
   to the finite point `(0,b,a-4*b^2)`, a contradiction.
2. If `Q=0`, the cubic (or its nonzero quadratic specialization at `c=0`)
   has a repeated finite root `T0`.  Perturb `T0`, define `a` from `p(T)=0`,
   and use reconstruction.  Then `r(T) -> 0`, so `x=2/r` escapes to infinity
   while the targets converge to `(a,b,c)`.

Acceptance gate:

- Lean certifies the cubic, derivative, discriminant, reconstruction, and
  large-`T` cancellation identities with no `sorry`, `admit`, or custom axiom;
- an independent typed SymPy checker verifies those identities and hostile
  perturbations exactly;
- the readable proof handles the `x=0`, bounded-`T`, and unbounded-`T` cases
  explicitly; and
- no historical-priority claim is made without a literature review.

The complete set-theoretic consequence to be recorded is stronger than the
initial target.  Let

```text
Gamma = { (a,b,c) : 3*b*c = 4 and b^2 = 12*a }.
```

Simple projective roots of the homogenized cubic correspond bijectively to
source points.  Therefore the target strata have fiber cardinalities `3`,
`1`, and `0` according as `Q != 0`, `Q = 0` off `Gamma`, or the target lies on
`Gamma`.  In particular the exact image is `C^3 \ Gamma`.  The algebraic
identities are machine checked; the projective root count and compactness
argument remain explicitly identified as derived mathematical prose.

## Milestone 4: a quadratic-in-one-variable plane fragment

For a characteristic-zero field, consider

```text
(x,y) |-> (a(x)*y^2 + b(x)*y + c(x), e(x)*y + f(x)).
```

Equating the coefficients of a nonzero constant Jacobian gives

```text
a'*e - 2*a*e' = 0,
b'*e - b*e' = 2*a*f',
c'*e - b*f' = k.
```

On the chart `e != 0`, the target normal form is

```text
r = eps*y + f(x),
F(x,y) = (lam*r^2 + mu*r + alpha*x + beta, r),
alpha*eps = k.
```

This has a displayed polynomial inverse.  The `e=0` chart must be handled
separately: the equations force `a=0`, `b` to be a nonzero constant, and `f`
to be affine with nonzero slope, which is again triangular and invertible.

Acceptance gate:

- Lean derives the normal form from the coefficient identities;
- Lean proves both inverse identities and the full bijectivity theorem,
  including the `e=0` branch;
- no theorem is described as the full plane Jacobian conjecture; and
- prior mathematical coverage of this bounded-degree class is cited, so the
  result is presented as formalization rather than novelty.

## Milestone 5: arbitrary degree with one affine coordinate

Let `P : (K[X])[X]` be an arbitrary polynomial in `y` with coefficients in
`K[x]`, and let

```text
Q(x,y) = e(x)*y + f(x).
```

Over a characteristic-zero field, the target theorem is:

> If the formal Jacobian of `(P,Q)` is a nonzero constant `k`, then `(P,Q)`
> is a polynomial automorphism.

The theorem must cover both charts and expose the inverse-producing normal
forms.

1. If `e != 0`, induct on the `y`-degree of `P`.  If `a_m(x)` is the leading
   coefficient, the top Jacobian equation is

   ```text
   a_m' * e = m * a_m * e'.
   ```

   Characteristic zero forces `a_m = lambda * e^m`.  Subtracting
   `lambda * Q^m` lowers the `y`-degree and leaves the Jacobian unchanged.
   Iteration gives

   ```text
   P(x,y) = G(Q(x,y)) + alpha*x + beta,
   e = eps,
   alpha*eps = k,
   ```

   with `alpha` and `eps` nonzero.  This has the displayed inverse

   ```text
   x = (u - G(v) - beta) / alpha,
   y = (v - f(x)) / eps.
   ```

2. If `e = 0`, the Jacobian equation is `-P_y*f' = k`.  Since a product of
   polynomials is a nonzero scalar, `P_y` and `f'` are nonzero constants.
   Thus

   ```text
   P(x,y) = eta*y + c(x),
   f(x) = delta*x + gamma,
   -eta*delta = k,
   ```

   which is triangular with an explicit inverse.

Acceptance gate:

- Lean defines the two formal partial derivations on `(K[X])[X]` and proves
  the stated normal forms from an equality of the actual formal Jacobian;
- Lean proves both inverse laws and full bijectivity in both charts;
- all `y`-degrees are allowed, including the constant and zero-polynomial
  edge cases;
- exact executable tests check generated normal forms and reject perturbed
  nonconstant-slope fixtures;
- `lake build`, the exact Python suite, strict type checking, and an axiom/
  unfinished-proof audit all pass; and
- the literature note distinguishes mathematical coverage from any claim of
  historical novelty and states explicitly that this does not solve `JC(2)`.

## Milestone 6: arbitrary degree with a constant-leading quadratic coordinate

Take the first genuinely nonlinear second coordinate after Milestone 5:

```text
P in K[x,y] arbitrary,
Q(x,y) = eps*y^2 + g(x)*y + f(x),
eps in K^*.
```

Put

```text
s = 2*eps*y + g(x),
Delta = g(x)^2 - 4*eps*f(x).
```

The identities

```text
s^2 - 4*eps*Q = Delta,
J(s,Q) = Delta'/2
```

lead to the now-certified theorem:

> If `J(P,Q)=k` for `k in K^*`, then there are nonzero scalars `A` and
> `lambda`, a scalar `B`, and `G in K[T]` such that
>
> ```text
> Delta = A*x + B,
> lambda*A/2 = k,
> P = G(Q) + lambda*s.
> ```

Indeed, evaluate the Jacobian on the polynomial critical section
`y=-g/(2*eps)`, where `Q_y=s=0`.  This gives

```text
P_y(x,-g/(2*eps)) * Delta'(x) = 4*eps*k.
```

Both factors are therefore units, so `Delta'=A` is a nonzero scalar.  The
pair `(s,Q)` is then a polynomial coordinate system with inverse

```text
x = (s^2 - 4*eps*Q - B) / A,
y = (s - g(x)) / (2*eps).
```

In these coordinates, the Jacobian equation says that the derivative of `P`
with respect to `s` is the scalar `lambda=2*k/A`, yielding the normal form.
For target coordinates `(u,v)`, the displayed inverse is

```text
s = (u - G(v)) / lambda,
x = (s^2 - 4*eps*v - B) / A,
y = (s - g(x)) / (2*eps).
```

Acceptance gate:

- the result is hostile-audited for signs, characteristic, zero-slope, and
  degree-drop edge cases before implementation;
- Lean starts from the actual formal bivariate Jacobian and proves the affine
  discriminant, coordinate-pair, normal-form, and two-sided-inverse theorems;
- the proof reuses the bundled derivations and arbitrary-degree descent from
  `AffineCoordinate.lean` rather than restating coefficient hypotheses;
- an exact typed symbolic checker validates high-degree examples and a
  perturbed non-Keller fixture;
- the literature boundary is sourced before any novelty language is used;
  absent such evidence, the result is presented only as a formalized
  structural fragment; and
- the repository continues to state that the unrestricted degree-six sheet
  stratum and `JC(2)` remain open.

Status: complete.  The formal statement starts from the actual bivariate
Jacobian, proves the original-coordinate normal form, constructs the displayed
inverse, and proves both inverse laws.  The independent exact checker exercises
a first coordinate of total degree 84 and rejects a one-coefficient
perturbation whose determinant becomes `15-36*x`.

## Milestone 7: direct certificate for a variable-leading quadratic coordinate

Let

```text
P in K[x,y] arbitrary,
Q(x,y) = a(x)*y^2 + g(x)*y + f(x),
a != 0,
J(P,Q) = k in K^*.
```

The automorphy conclusion is a known theorem, not a new open case.
Moskowicz's Theorem 2.7 applies directly because the invariant attached to the
leading `y`-coefficient is

```text
gcd(2, deg_x(a)) in {1,2}.
```

The repository's research target is a separate, elementary, machine-checked
route that exposes why polynomiality rules out a variable leading coefficient:

1. subtract scalar powers of `Q` until the remaining mate has odd `y`-degree
   `n`;
2. derive `p_n^2=c*a^n`, then use unique factorization to write
   `a=eps*h^2` and `p_n=lambda*h^n`;
3. over `K(x)`, put `U=h*y+g/(2*eps*h)` and split uniquely as
   `P=H(Q)+U*B(Q)`;
4. transport the Jacobian identity to derive the coefficient recurrence for
   `B`;
5. specialize at `y=0` and use a gcd-normalized unique-survivor denominator
   argument to prove `h | g`; and
6. use the constant term of the recurrence to prove `h | k`, hence `h` is a
   unit and `a` is a nonzero scalar.

Current Lean status: steps 1 and 2 are certified from the actual bivariate
Jacobian.  Step 3 is instantiated over `K(x)`, including the quotient-rule
derivation, its constant field, affine transport, parity decomposition, and
degree/leading-coefficient control.  In step 4, Lean certifies both the
centered Jacobian identity with right-hand side `k/h`, the resulting
coefficient recurrence, its exact solution in `K[F]`, and the required degree
and nonzero-leading-coefficient descent.  Steps 5 and 6 are now certified as
well: specialization supplies a base-ring membership, gcd reduction exposes
the unique surviving numerator modulo the reduced denominator, and the
constant recurrence coefficient forces the remaining denominator to be a
unit.  The final theorem proves bijectivity, including the affine `a=0`
branch.  A typed exact checker separately covers the rational-mate trap, the
top equation, the central-binomial recurrence, and the forced infinite-tail
obstruction.

Simon--Weimann's coordinate criterion supplies a useful known corollary: after
automorphy, the leading coefficient `a` is constant and
`deg_x(g^2-4*a*f)=1`.  This locates the completed direct Lean theorem in its
known structural context.

## Milestone 8: exact monodromy obstruction at generic degree six

Let

```text
K = C(P,Q),
L = C(x,y),
[L:K] = 6,
G = Gal(M/K)
```

for the Galois closure `M` of a hypothetical plane Keller counterexample.
The target is a complete necessary-condition theorem, not a claimed solution
of `JC(2)`:

1. prove that every generic branch meridian fixes an affine sheet and that
   such meridians normally generate `G`;
2. enumerate all sixteen transitive subgroups of `S6` exactly and retain only
   those normally generated by fixed-point elements;
3. apply Orevkov's exact defect budget and the finite-normalization deck group
   to the remaining imprimitive groups;
4. use purity and local monodromy to handle singular branch germs without
   silently assuming they are smooth; and
5. refine the defect budget by tangential Riemann--Hurwitz, eliminate every
   sole ramified branch moving more than three sheets, and record the exact
   unrestricted and one-dicritical passports which remain.

Acceptance gate:

- every finite-group assertion is rebuilt from checked GAP generators by a
  dependency-free typed implementation and independently compared with
  GAP 4.14 through Sage;
- same-cycle-type conjugacy classes remain distinct throughout the audit;
- Orevkov's local multiplicity is never confused with global sheet degree;
- any use of an irreducible local branch explicitly proves that it is the
  entire local discriminant, with a hostile extra-branch countermodel guarding
  the distinction; and
- surviving groups are presented only as necessary abstract passports, never
  as realized Keller maps.

Status: complete.  Fixed-sheet inertia reduces sixteen groups to seven.
Orevkov's budget, boundary deck symmetry, purity, and the local
double-transposition lemma first eliminate `6T7`, `6T8`, and `6T11`.
Tangential Riemann--Hurwitz sharpens the exact budget to

```text
sum_E (e_E*d_E + delta_E) = 5,  delta_E >= 0.
```

The resulting moved-sheet budget and irreducible-branch topology obstruction
also eliminate `6T12=A5`, `6T14=S5`, the double-transposition `A6` profile,
and the four-cycle and same-branch `(3)(2)` `S6` profiles.  Thus every
hypothetical six-sheet counterexample has monodromy

```text
6T15 = A6  or  6T16 = S6.
```

In the unrestricted `A6` case there is one ramified `(e,d)=(3,1)` branch with
noninjective normalization.  In the unrestricted `S6` case the ramified
profiles are one transposition branch, two distinct transposition branches,
or distinct transposition and 3-cycle branches.

With one dicritical component, deck rigidity first leaves the same four
groups.  Local multiplicity sharpens the `(e,d)=(2,2)` and `(4,1)` passports,
and Lin--Zaidenberg topology then eliminates both: an injectively normalized
branch would have identical local and global complement groups, contradicting
the intransitive local sheet action.  Only `6T15=A6` with `(3,1)` and
`6T16=S6` with `(2,1)` remain, and their branch normalization must identify
points.  The surviving `A6` and `S6` profiles and `JC(2)` remain open.

For the one-dicritical `A6` survivor, the raw jump partition `1+1` is now
excluded.  Orevkov's one-chain structure initially gives at most one singular
source point on the boundary.  At a smooth point the Jacobian divisor is
twice the smooth ramification prime; rank-one parity and a corank-two quadratic-jet
argument both rule out local degree four.  Hence the two jump units form one
unique local degree five.  A smooth such point has rank one and a forced
`(2,5)` cusp.  Pulling the exceptional singular-source alternative to its
Hirzebruch--Jung universal cyclic cover forces an `A1` quotient, then an exact
intersection count forces an even contact order where deck invariance forces
an odd one.  The same argument applies to a contracted generic
local-degree-three point: index three forces every contracted endpoint to have
even local degree, whereas the exact boundary degrees are three and five.
Thus `L_C` is empty and the whole affine finite normalization is smooth.  An
exact hostile fixture
verifies that its knot group can still map onto local `A5` and extend with a
collision meridian to global `A6`, so this remains a reduction rather than an
elimination.

For the one-dicritical `S6` survivor, writing `kappa(t)=mu_t-2` gives
`sum kappa=3` and the exact fiber formula
`#F^-1(y)=6-2*#nu^-1(y)-kappa_y`.  The resulting nine rows are executable.
The connected block at `t` is `S_(2+kappa(t))`, because it is transitive and
normally generated by a transposition.  A nonempty contracted chain would
lift to a small cyclic quotient cover whose invariant coordinate functions
start in degree two, while its pulled-back Jacobian has order one.  This is
impossible, so `L_C` is empty.  Every jump therefore has branch valuations
`(m-1,m)` for `m=2+kappa`, yielding `T(2,3)`, `T(3,4)`, or `T(4,5)`.
Explicit germs `(a,z) -> (a,z^m+a*z)` realize all three types, so the result is
an exact local classification rather than a contradiction.

In both one-dicritical survivors, the normalization `W -> A2` is now a smooth
finite-flat map of rank six with `W-D=A2`, `D=A1`,
`Pic(W)=Z[D]`, and `K_W=(e-1)D`.  Quillen--Suslin makes the underlying
rank-six algebra module free, but does not split its multiplication.

On Orevkov's original pure source-blowup resolution, emptiness of `L_C` also
makes the unique dicritical `E` a leaf.  Its logarithmic Keller coefficient is
its normal ramification index `e`; adjunction therefore gives its unique
neighbor the label

```text
e*(-E^2) - 1.
```

That label is positive and congruent to `-1 mod e`, so the neighbor is type 2.
The path to the original line at infinity must cross `1--0--(-1)`.  This is a
new finite-graph constraint, not an elimination.  A further joint minimization
contracts every boundary `(-1)`-curve collapsed by both the source blowdown
and the resolved map.  In that model no corner blowup can remain above the
leaf edge, so `E^2=-1` and the neighbor label is exactly `e-1`: `2--1` for
`S6` and `3--2` for `A6`.  Later target blowups must still transform this leaf
configuration rather than reuse it unchanged.

The smooth pair invariants do not close the argument.  Explicit complements
of ample sections in Hirzebruch surfaces realize both one-boundary canonical
classes with `W-D=A2`, `D=A1`, and `Pic(W)=Z[D]`.  A separate three-blowup
construction realizes the saturated disjoint two-boundary package and proves
affineness by an explicit ample divisor.  None carries a finite rank-six
Keller morphism.  Consequently the next gate is genuinely the multiplication,
trace, discriminant, and monodromy of the finite algebra, not surface topology
or the displayed Picard/canonical classes in isolation.

Finite duality does sharpen the remaining algebra.  The relative dualizing
module `Hom_R(A,R)` is the inverse-different line and represents `D`, `2D`, or
`D2+2D3` in the free source Picard group.  It is therefore not free as an
`A`-line even though it is free over `R=C[P,Q]`.  Consequently `A/R` is not a
globally Frobenius algebra, is not monogenic, and admits no single global
square complete-intersection presentation; the finite morphism remains
syntomic and locally complete intersection.

Writing the ramification divisor as `2E+O` gives a free half-different lattice
and the exact symmetric factorization

```text
T = Phi^* H Phi,
det(T) = det(H) * det(Phi)^2.
```

The three trace-determinant rows are `b`, `b^2`, and `b2*b3^2`.  Divisor sequences
identify the cokernels with the branch-normalization modules and fix every
matrix corank at the known cusps and collisions.  All those ranks remain
possible, so this is a concrete algebraic endpoint rather than an exclusion.

## Milestone 9: sparse support in the residual `(72,108)` polygons

Keep this coordinate-degree problem separate from generic sheet degree.
Starting from the two transformed Newton-polygon pairs in Guccione--Guccione--
Horruitiner--Valqui, Proposition 4.3, exhaust every support with at most two
nonzero strict-interior coefficients in the first case and at most three in
the second, while allowing all boundary lattice coefficients.

Acceptance gate:

- rebuild every bracket coefficient from
  `(i*l-j*k)*p_ij*q_kl` using exact integers;
- require nonzero coefficients only at exact polygon vertices and at selected
  interior support points;
- replay every forced-zero certificate from rebuilt equations before counting
  it;
- include hand-checkable boundary-only identities, exact algebraic replays for
  any propagation exceptions, and hostile full-polygon and next-support-size
  fixtures; and
- do not infer satisfiability from the point where zero-product propagation
  stops.

Status: complete at the stated sparse boundary.  The first configuration has
`122` possible interior coefficients and all `7504` supports of size at most
two are impossible.  The second has `28` possible interior coefficients and
all `3683` supports of size at most three are impossible.  Zero-product
propagation proves `3678` of the second-case supports; exact Groebner
unit-ideal replays prove the five exceptional triples.  Hence Case 1 needs at
least three nonzero strict-interior coefficients and Case 2 needs at least
four.  The full coefficient systems remain open.

## Repository shape

```text
JacobianTwo/
  Counterexample.lean       # formal polynomial and collision certificate
  ConstantLeadingQuadratic.lean # arbitrary P, scalar-leading quadratic Q
  VariableLeadingQuadratic.lean # degree descent, UFD shape, and centering
  RatFuncDerivative.lean # quotient-rule derivative and constant field of K(x)
  FractionRingDerivative.lean # derivative on the fraction-ring representation
  QuadraticParityJacobian.lean # centered parity Jacobian identity
  QuadraticParityExtraction.lean # odd/even separation and recurrence extraction
  QuadraticTransportJacobian.lean # exact affine Jacobian factor k/h
  QuadraticReduction.lean # actual Keller pair to centered equation
  QuadraticRecurrence.lean # coefficient equations for the odd part
  QuadraticRecurrencePrimitive.lean # exact descent in K[F] with degree control
  QuadraticDenominator.lean # noncancellation and unit endpoints
  QuadraticClearedEval.lean # exact common numerators for C(N/q^2)
  QuadraticDenominatorDescent.lean # valuation-free unique-survivor proof
  QuadraticSpecialization.lean # specialization at the quadratic center
  QuadraticDirectPackage.lean # end-to-end odd-degree certificate
  QuadraticGCDDescent.lean # reduced denominator and h | g
  QuadraticUnitConclusion.lean # terminal recurrence forces h to be a unit
  VariableLeadingQuadraticConclusion.lean # full at-most-quadratic bijectivity
  AffineInOneVariable.lean  # JC(2) obstruction theorem
  QuadraticInOneVariable.lean # quadratic-in-y normal form and inverse
  CubicFiber.lean           # cubic, discriminant, reconstruction, infinity algebra
  Basic.lean                # shared definitions
scripts/
  verify.py                 # independent exact symbolic checker
  nonproper.py              # exact cubic-fiber/nonproper algebra
  constant_leading_quadratic.py # high-degree quadratic-coordinate checker
  variable_leading_quadratic.py # denominator and infinite-tail fixtures
  six_sheet_monodromy.py # exact transitive-group and local-orbit certificate
  newton_72_108.py       # exhaustive sparse Newton-support certificates
tests/
  test_verify.py            # positive and adversarial fixtures
  test_nonproper.py         # discriminant, strata, and hostile fixtures
  test_constant_leading_quadratic.py # normal form and perturbation checks
  test_variable_leading_quadratic.py # hostile arbitrary-leading fixtures
  test_six_sheet_monodromy.py # group catalogue, budgets, and local passports
  test_newton_72_108.py   # sparse-support exhaustion and hostile boundaries
tools/
  check_six_sheet_gap.sage # independent Sage/GAP group cross-check
docs/
  audit.md                  # readable derivation, provenance, claim boundary
  galois-frontier.md        # Galois theorem and first open sheet degree
  nonproper-set.md          # complete fiber/image/nonproper proof
  constant-leading-quadratic.md # proof and exact Lean declaration map
  variable-leading-quadratic.md # known theorem and direct proof target
  six-sheet-monodromy.md # primitive degree-six monodromy obstruction
  newton-72-108-sparse.md # residual coordinate-degree sparsity theorem
  research-log.md           # dated attempts and open obligations
```

## Reproducibility policy

- Pin Lean and mathlib to the same released version.
- Use `uv` for every Python command and dependency.
- Treat Lean as the authoritative certificate and SymPy as an independent
  implementation cross-check.
- Record negative results and failed reductions; do not silently turn them into
  conjectural prose.
- Do not claim a result about `JC(2)` merely because `JC(3)` is false.

## Source boundary

The originating public record is
[Alpöge's X post](https://x.com/__alpoge__/status/2079028340955197566).
The repository may cite expository follow-ups for context, but its decisive
algebraic claims must remain independently checkable from the formula above.
