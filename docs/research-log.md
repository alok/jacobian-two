# Research log

## 2026-07-20: source audit

Status: complete.

- Transcribed the polynomial map and all four displayed points from the
  screenshot.
- Matched the transcription against Levent Alpöge's public announcement and
  independent same-day expositions.
- Classified provenance claims as `ANNOUNCED`, because the public record is
  only hours old.
- Classified the determinant and collision as finite identities suitable for
  exact independent verification.

## 2026-07-20: exact counterexample certificate

Status: Lean-certified.

- Defined the three coordinate polynomials over `ℤ`.
- Constructed the formal Jacobian using `MvPolynomial.pderiv`.
- Proved its determinant is the constant polynomial `-2`.
- Transported the pointwise determinant statement to arbitrary commutative
  rings.
- Proved the three complex points are pairwise distinct and have the same
  image.
- Derived noninjectivity.
- Audited theorem dependencies with `#print axioms`; only Lean's standard
  foundations occur.

## 2026-07-20: independent exact checker

Status: complete.

- Added a typed SymPy implementation using exact rational arithmetic.
- Added a positive determinant-and-fiber test.
- Added a coefficient perturbation that destroys the constant determinant.
- Added a point perturbation that destroys the common fiber.
- Added a duplicate-input fixture that rejects a vacuous three-point claim.
- Checked the Python under strict mypy settings.

## 2026-07-20: structural factorization

Status: algebraic core Lean-certified.

- Introduced `u=1+xy`, `H=u^2 z+y^2(4+3xy)`, and `s=x/u`.
- Reduced fibers on the `u != 0` chart to
  `2*a*s^3-b*s^2+2*s-c=0`.
- Derived the rational reconstruction formulas and the identity
  `Phi'(s)=2*(1-b*s+3*a*s^2)`.
- Factored the displayed fiber as `-s*(s-2)*(s+2)/2`.
- Exhibited the escaping curve `(s,-1/s,5/s^2)`, whose image is `(0,2/s,0)`.
- Added the reciprocal `T=1/s` cubic, derivative, reconstruction, and
  large-`T` cancellation identities to `JacobianTwo/CubicFiber.lean`.

The chart restriction matters. These formulas illuminate the example but do
not replace the global polynomial determinant proof.

## 2026-07-20: Galois misconception and first open sheet degree

Status: source-audited.

- Checked the classical Galois-case theorem against Campbell, Razar, and
  Bass–Connell–Wright.
- Isolated the missing hypothesis: a Keller extension is finite separable but
  need not be normal.
- Kept generic fiber degree separate from ordinary coordinate degree.
- Located the accepted low-sheet boundary: Orevkov excludes degree three and
  Żołądek excludes through degree five, so degree six is the first unresolved
  generic fiber degree.
- Rejected the withdrawn/corrected Bartenwerfer plane proof as a dependency.
- Treated a 2024 prime-degree preprint as unconfirmed rather than silently
  promoting it to an accepted theorem.

## 2026-07-20: complete fibers, image, and nonproper set

Status: exact algebra machine-checked; projective/topological proof derived.

- Proved that the standard universal cubic discriminant coefficient
  expression is `-4Q` and added an explicit integral-coefficient Bézout
  certificate.
- Proved finite simple-root reconstruction and all cancellations needed at
  the projective root at infinity.
- Identified the singular/triple-root curve
  `Gamma={3*b*c=4, b^2=12*a}`.
- Derived fiber counts `3/1/0` on the three target strata and the exact image
  `C^3 \ Gamma`.
- Constructed an exact escaping family through every point of `V(Q)`.
- Proved conversely, by projective-root compactness and reconstruction, that
  targets outside `V(Q)` are not asymptotic values.
- Concluded `S_F=V(Q)`.

The full fiber/nonproper theorem is labeled “derived here; historical priority
not established.” Same-day comparison sources already contain the cubic,
reconstruction, discriminant, and generic `S_3` computation.

## 2026-07-20: first plane obstruction

Status: Lean-certified.

Tested the closest plane analogue

```text
(x,y) |-> (A(x)*y+B(x), C(x)*y+D(x)).
```

A nonzero constant determinant makes `A` and `C` coprime and gives them a
vanishing Wronskian. Lean proves `A'=C'=0`; in characteristic zero both slopes
are constant. A second Lean theorem combines the outputs to recover `x` and
then `y`, proving that the whole map is injective. This is a meaningful
noncollision result toward `JC(2)`, but it is far from the full conjecture.
The later quadratic-in-one-variable module contains this affine class and
packages it inside an explicit two-sided-inverse theorem.

## 2026-07-20: quadratic-in-one-variable plane fragment

Status: Lean-certified formalization of a known positive class.

For maps

```text
(x,y) |-> (a(x)*y^2+b(x)*y+c(x), e(x)*y+f(x)),
```

the three constant-Jacobian coefficient equations force a triangular normal
form.  On `e != 0`, writing `r=eps*y+f(x)` gives

```text
F(x,y) = (lam*r^2+mu*r+alpha*x+beta, r),
alpha*eps = k.
```

The Lean module supplies the displayed inverse and proves both inverse laws.
The complementary `e=0` case is triangular after the equations force `a=0`,
`b` constant nonzero, and `f` affine with nonzero slope.  The final theorem
combines the two branches.  This narrows a known class; it does not reach the
generic-degree-six frontier.

## 2026-07-21: arbitrary-degree affine-coordinate theorem

Status: Lean-certified formalization of a known positive class.

For an arbitrary first coordinate `P in K[x,y]` and

```text
Q(x,y) = e(x)*y + f(x),
```

Lean now defines both bivariate partial derivatives as bundled derivations and
starts from the actual formal identity `J(P,Q)=k`.  Strong induction removes
the top `y`-power of `P` by subtracting a scalar multiple of a power of `Q`.
The final theorem covers both charts:

```text
e != 0:  e=eps, P=G(Q)+alpha*x+beta, alpha*eps=k;
e = 0:   P=eta*y+c(x), Q=delta*x+gamma, -eta*delta=k.
```

Both normal forms have explicit kernel-checked two-sided inverses.  An exact
SymPy checker exercises a degree-seven outer polynomial and rejects a
nonconstant-slope perturbation.  The type-`(m,1)` elimination is known from
Sabatini's 2024 theorem; the repository contribution is the
characteristic-zero field-uniform formalization and explicit packaging, not a
novelty claim.

## 2026-07-21: constant-leading quadratic-coordinate theorem

Status: Lean-certified structural positive class; no novelty claim.

For arbitrary `P in K[x,y]` and

```text
Q = eps*y^2 + g(x)*y + f(x),  eps != 0,
s = 2*eps*y + g(x),
Delta = g(x)^2 - 4*eps*f(x),
```

the actual identity `J(P,Q)=k != 0` forces `Delta=A*x+B` with `A != 0`.
The completed-square pair `(s,Q)` is then a polynomial coordinate system.
After transporting the Jacobian equation through this coordinate change, Lean
proves

```text
P = G(Q) + lambda*s,  lambda*A/2 = k,
```

and verifies both laws of the resulting explicit inverse.  The independent
exact checker uses a degree-six `g`, a degree-seven `G`, and a first coordinate
of total degree 84.  Changing one coefficient of `f` makes the determinant
`15-36*x`, so the hostile fixture confirms that the affine-discriminant step is
doing real work.  This reaches arbitrary degree in the other coordinate but
still assumes that the coefficient of `y^2` is a nonzero scalar.

## 2026-07-21: variable-leading quadratic checkpoint

Status: known theorem; independently derived direct proof; partially
Lean-certified.

For

```text
Q = a(x)*y^2 + g(x)*y + f(x),  a != 0,
J(P,Q) = k != 0,
```

Moskowicz's Theorem 2.7 already proves that `(P,Q)` is an automorphism.  Its
invariant is `gcd(2,deg_x(a))`, which is either `1` or the prime `2`, so the
theorem matches an arbitrary characteristic-zero field and places no degree
bound on `P`.  Simon--Weimann's coordinate criterion then gives the stronger
corollaries that `a` is a nonzero scalar and `deg_x(g^2-4*a*f)=1`.  The exact
theorem statement is therefore not novel.

The repository nevertheless has a useful new certificate route.  From the
actual bivariate Jacobian, Lean now proves the top-coefficient differential
equation, all even target-shear steps needed to reach odd degree, the identity
`p_n^2=c*a^n`, and the UFD consequence

```text
a=epsilon*h^2,  p_n=lambda*h^n.
```

Over any field, Lean also proves existence and uniqueness of

```text
R(U) = H(epsilon*U^2+F) + U*B(epsilon*U^2+F),
```

together with the exact degree and leading coefficient of `B` when `R` has
odd degree.  The remaining formal boundary is to instantiate this over
`K(x)`, transport the `x`-derivation at fixed `Q`, extract the recurrence, and
formalize the denominator/valuation proof that first forces `h | g` and then
forces `h` to be a unit.  Exact hostile fixtures distinguish rational from
polynomial mates and exhibit the nonterminating odd coefficient tail for
`Q=x+x^2*y^2`.

## 2026-07-21: fraction-field quadratic bridge

Status: known theorem; recurrence descent and its direct-proof infrastructure
are Lean-certified; specialization and denominator integration remain.

The `K(x)` boundary above is now crossed.  The repository constructs the
quotient-rule derivation on `RatFunc K`, proves that its kernel is exactly the
constant field in characteristic zero, and transports it to the generic
`FractionRing K[X]` used by the affine chart.  The substitution

```text
y = (U-rho)/h
```

is certified coefficientwise.  A general composition theorem proves that the
centered Jacobian is multiplied by `1/h`; consequently an original Keller
equation `J(P,Q)=k` becomes exactly `k/h`.  The parity identity and uniqueness
then separate the odd part and produce both `D(H)=0` and the complete
coefficient recurrence for `B`.  Independent denominator lemmas certify the
unique-survivor obstruction and the implication “a unit numerator divided by
`h` lies in the base ring, therefore `h` is a unit.”

The recurrence gap is now closed: Lean constructs an explicit primitive for
each downward step, proves every coefficient is an evaluation of a polynomial
in `K[F]`, and certifies exact degree growth with a nonzero leading
coefficient.  The remaining direct-proof gap is narrower: use specialization
of the original polynomial mate at `y=0` to supply the unique-survivor
hypotheses forcing `h | g`, and then prove that the terminal `k/h` equation
lies in `K[x]`.  No automorphy claim is promoted until those links are
kernel-checked.

## 2026-07-21: complete variable-leading quadratic certificate

Status: known theorem; independently derived direct proof now fully
Lean-certified.

The two remaining links above are now closed.  Specialization of the original
mate at the quadratic center proves that `rho * B(f)` belongs to the embedded
polynomial ring.  If `h` did not divide `g`, Lean reduces `g/h` by its gcd,
writes the centered residual as `N/q^2`, and clears the common denominator
`q^(2*m+1)`.  Every numerator summand is divisible by `q` except the exact
leading term of the zero-index coefficient, whose nondivisibility follows
from coprimality and the recurrence's degree/leading-coefficient certificate.
This contradiction forces `h | g` without a valuation API.

Once `h | g`, the center, residual, and zero-index recurrence coefficient are
all embedded polynomials.  The terminal recurrence equation therefore puts
`k/h` in the base ring; its nonzero scalar numerator forces `h` to be a unit.
The general target-shear and UFD steps then prove that the original quadratic
coefficient `a` is a unit.  The theorem
`variableLeadingQuadratic_bijective_full` reduces the nonzero branch to the
constant-leading result and the `a=0` branch to the affine result.

This completes a machine-checked direct proof of a theorem already implied by
Moskowicz's Theorem 2.7.  It does not resolve unrestricted `JC(2)` or the
generic six-sheet frontier.

## 2026-07-21: primitive monodromy forced at generic sheet degree six

Status: derived here and exactly enumerated; mathematical priority not
established; `JC(2)` remains open.

The Galois audit was pushed into the first accepted open sheet degree.  For a
hypothetical degree-six plane Keller map, the omitted target set is finite.
Consequently, every generic nonproperness-curve meridian has an affine local
inverse branch and fixes at least one of the six sheets.  Since component
meridians normally generate the complement group, the monodromy group is
normally generated by fixed-point elements.  Exact enumeration of all sixteen
transitive degree-six groups leaves

```text
6T7, 6T8, 6T11, 6T12, 6T14, 6T15, 6T16.
```

The finite normalization then supplies more geometry.  A nontrivial field
deck group forces an unramified dicritical prime; otherwise every deck
automorphism preserves the entire ramified boundary, restricts to a finite
order polynomial automorphism of `A2`, and is killed by local invertibility at
a fixed point.  Combining this with Orevkov's exact degree-six defect budget
eliminates `6T11` and forces the saturated profiles used for `6T8` and `6T7`.
The `6T8` profile has a smooth embedded branch line, whose complement has
cyclic fundamental group.  The two possible `6T7` realizations are handled
separately: the zero-jump realization again forces a smooth line, while the
quadratic tangential realization has one local degree-four component with
double-transposition inertia.  Purity makes its unibranch image the entire
local discriminant; an irreducible branch meridian cannot normally generate a
transitive four-sheet group from `(12)(34)`.  This removes the last
imprimitive action.

The universal surviving list at this intermediate milestone was

```text
6T12 = A5, 6T14 = S5, 6T15 = A6, 6T16 = S6.
```

Under the additional one-dicritical hypothesis, finite-normalization deck
rigidity first gives the same list.  The local analysis forces the
`(e,d)=(2,2)` passport to spend all three jump units at one multiplicity-five
critical point, and the `(4,1)` passport to spend its one jump unit at one
multiplicity-five point.  In either case the branch normalization is
injective.  A finite bijective normalization makes the branch topologically
contractible; Lin--Zaidenberg puts it in the monomial form `X^m=Y^n`, whose
weighted radial action identifies the global complement group with its local
knot-complement group.  The local finite cover has more than one sheet orbit,
whereas the global cover is connected.  This eliminates both passports.

Consequently a one-dicritical degree-six counterexample can only have
`6T15=A6` with `(e,d)=(3,1)` or `6T16=S6` with `(2,1)`.  In either case the
normalization `A1 -> B` must identify points.  For `A6`, each collision is
exactly two smooth multiplicity-three branches, exhausts all six sheets, and
lies outside the image of the original affine Keller map.

The pure-Python certificate reconstructs all groups, conjugacy classes,
normal closures, normalizers, deck quotients, block systems, budget profiles,
and relevant local subgroup orbits.  Sage/GAP 4.14 independently reproduces
the catalogue and every conjugacy-class invariant.  The next milestone below
refines this four-group boundary further.

## 2026-07-21: tangential defect budget leaves only `A6` and `S6`

Status: derived here, independently hostile-audited, and exactly enumerated;
mathematical priority not established; generic degree six remains open.

For each dicritical `E`, Orevkov's affine-line structure factors its boundary
map as a polynomial `A1 -> A1` of tangential degree `d_E`.  If `a_c` is its
local tangential degree, local sheet counting gives
`mu_(pi(c)) >= e_E*a_c`, while Riemann--Hurwitz gives
`sum_c(a_c-1)=d_E-1`.  Orevkov's exact degree-six identity therefore sharpens
to

```text
sum_E (e_E*d_E + delta_E) = 5,  delta_E >= 0.
```

Thus a generic branch meridian costs at least the number of sheets it moves.
The exact class enumeration and total cost five leave only a short list of
normally generating primitive profiles.  If the sole ramified branch moves
four or five sheets, two normalization points cannot collide within degree
six, so the normalization is injective.  Lin--Zaidenberg then makes the branch
contractible and identifies its local and global complement groups; the local
finite cover is disconnected because every boundary local degree is at most
five, while the global cover is connected.  This contradiction eliminates
`A5`, `S5`, the `A6` double-transposition profile, and the `S6` four-cycle and
same-branch `(3)(2)` profiles.

The unconditional frontier is now

```text
6T15 = A6 or 6T16 = S6.
```

The `A6` case has one ramified `(e,d)=(3,1)` branch with noninjective
normalization.  Every normalization collision is exactly two smooth
multiplicity-three boundary points, exhausts all six sheets, and is omitted by
the affine map.  The `S6` ramified profiles are exactly one transposition
branch, two distinct transposition branches, or distinct transposition and
3-cycle branches.  The executable certificate now enumerates these complete
cost-compatible profiles and preserves the forbidden index-five classes as
hostile fixtures before applying Orevkov's final remark.

## 2026-07-21: the one-dicritical `A6` jump is unique

Status: derived here and hostile-audited; conditional local theorem, not an
elimination of `A6`.

Orevkov's one-dicritical boundary is one linear `L_C-E` chain, so contraction
creates at most one singular finite source point on the dicritical image.  At
every smooth source point the local Jacobian divisor is twice the unique
smooth ramification prime.  In rank one this forces every local degree to be
odd.  In corank two, local degree four would give two coprime quadratic jets
whose Wronskian is a doubled line, contradicting Riemann--Hurwitz for the
induced degree-two map of `P1`.  Hence a smooth local-degree-four point is
impossible.

The `A6,(e,d)=(3,1)` defect partition `1+1` would require two such points, at
least one smooth, and is therefore excluded.  There is a unique jump point of
local degree five.  If its source is smooth, a second initial-form/blow-up
argument rules out corank two; integration of the doubled Jacobian along the
ramification curve then gives target valuations `(2,5)`.  The intrinsic
singularity is a `(2,5)` cusp.

This does not close the passport.  The exact permutations `(345)` and `(123)`
satisfy the `T(2,5)` five-braid relation and generate local `A5`; adjoining
the collision meridian `(456)` generates global `A6`.  The regression test
locks this hostile stopping fixture.  This milestone initially left a
dichotomy between the smooth `(2,5)` cusp source and the unique
normal-surface singularity obtained by contracting `L_C`; the later
Hirzebruch--Jung parity milestone below eliminates the latter.

## 2026-07-21: exact one-dicritical `S6` local passports

Status: derived here and hostile-audited; exhaustive local classification,
not an elimination of `S6`.

For the surviving one-dicritical `(e,d)=(2,1)` passport, set
`kappa(t)=mu_t-2`.  Orevkov gives `sum kappa=3`, while constant multiplicity
gives the exact target-fiber formula

```text
#F^-1(y) = 6 - 2*#nu^-1(y) - kappa_y.
```

This produces nine and only nine block rows, from the generic `2+4 affine`
row through the omitted `2+2+2` triple collision.  The executable certificate
locks every row and its omission flag.  The connected local block at a branch
with jump `kappa` is transitive and normally generated by a transposition, so
its group is exactly `S_(2+kappa)`.  Positive jump is equivalent to intrinsic
branch singularity.

At that stage, the one boundary chain left at most one exceptional source
point.  At any other jump, reduced simple ramification forces critical-image valuations
`(m-1,m)` for `m=2+kappa`.  Thus every nonexceptional jump branch is exactly
`T(2,3)`, `T(3,4)`, or `T(4,5)`.  A double-branched-cover argument gives an
`A_(2+kappa)` quotient even at the exceptional point.  The three jump
partitions are therefore `3`, `2+1`, and `1+1+1`, with explicit knot and
collision constraints.

The finite analytic germs `(a,z) -> (a,z^m+a*z)` realize all three smooth
torus rows.  A separate two-double-cover model realizes a zero-jump double
collision with two affine sheets and arbitrary branch contact.  These hostile
models show that local multiplicity and knot parity alone cannot finish the
case; the next input must couple the rows through global braid or
compactification-graph data.

## 2026-07-21: universal collisions in the two-curve `S6` passports

Status: derived here, source-checked, and exactly fixture-tested; a global
collision theorem, not an elimination of `S6`.

Nguyen Van Chau's Theorem 4.4(E3) says that every irreducible component of a
Keller map's deficient-fiber/nonproper set is singular.  When a branch
dicritical has `d=1` and zero excess, Orevkov makes its normalization map an
immersion everywhere.  Its required singularity must therefore identify two
distinct normalization points.

In the two-transposition profile the residual defect is one, so at least one
of the two ramified curves has zero excess and must self-collide.  In the
transposition-plus-3-cycle profile the costs `2+3=5` saturate the budget:
there are no index-one dicriticals, both excesses vanish, and both branch
curves self-collide.  Constant multiplicity gives the exact rows `3+3`,
`2+2+1+1` or `2+2+2`, and `2+3+1` at cross-intersections.

The disjoint local meridians `(123),(456),(14),(25)` nevertheless generate
all of `S6`.  This hostile fixture proves that the new collision census stops
at the global braid/splice/canonical realization problem rather than at a
permutation contradiction.

## 2026-07-21: the saturated `S6` finite normalization is smooth

Status: derived here and independently hostile-audited; sharpens the
transposition-plus-3-cycle passport without eliminating it.

The saturated `(e,d)=(2,1)+(3,1)` profile spends all five refined defect
units, so both componentwise excesses vanish and there are no other
dicriticals.  Tangential degree one then makes every finite local degree on
the two boundary primes exactly two and three, respectively.  Orevkov's two
separate minimal `L_C-E` chains satisfy the sole-local-ramification hypotheses
of the cyclic endpoint theorem.  Index two is impossible outright; index
three would force its local degree three to be even.  Hence both minimal
constant chains are empty.

Writing `D2,D3` for the disjoint source boundary primes, the finite
normalization now satisfies

```text
W - (D2 union D3) = A2,
D2 = D3 = A1,
Pic(W) = Z[D2] + Z[D3],
K_W = D2 + 2*D3.
```

It is smooth and finite flat of rank six, and its pushforward algebra is free
as a target module.  Target self-collisions and cross-intersections remain:
they identify images of distinct source points, not the disjoint source
divisors.  The remaining obstruction lies in the rank-six algebra
multiplication, trace form, and compactification, not in module freeness.

## 2026-07-21: the exceptional one-dicritical `A6` source is impossible

Status: derived here and independently hostile-audited; eliminates one local
alternative, not the smooth `(2,5)` passport.

If the unique local-degree-five point came from contracting the linear
constant chain, its normal source germ would be a Hirzebruch--Jung cyclic
quotient.  Pulling the finite germ back to the universal quasi-etale cover
gives a smooth-source map `h` of degree `5n` with ramification divisor twice
one smooth lifted endpoint curvette.  Its quadratic jets have Wronskian a
doubled line.  Riemann--Hurwitz and deck characters force the jets to be
`x^2,xy`, the cyclic cover to have order two, and its involution to be `-I`.

The lifted map is therefore even and has degree ten.  The curve `v=0` has two
smooth transverse branches.  The Jacobian chain rule makes their exact
intersection contributions `2` and `2r`, so local degree gives `r=4`.
Invariance of the unique branch tangent to `x=0` makes its graphing series
odd, forcing `r` odd: a contradiction.  Thus the degree-five source is
necessarily smooth and the earlier analysis makes its target branch exactly
the `(2,5)` cusp.

## 2026-07-21: the entire one-dicritical finite normalization is smooth

Status: derived here and independently hostile-audited; eliminates every
contracted finite source point in both surviving one-dicritical passports,
but does not eliminate either passport.

The Hirzebruch--Jung argument above does not depend on local degree five.  For
generic ramification index three it forces the universal quotient cover to
have order two and gives the exact contact equation `r=mu-1`; deck invariance
requires `r` odd, so every contracted endpoint must have even local degree.
The `A6` boundary census has only local degrees three and five.  Therefore its
constant chain `L_C` is empty.

For generic ramification index two, the pulled-back Jacobian has order one.
The small cyclic deck representation has no invariant linear form, so both
coordinate functions start in degree at least two and their Jacobian has
order at least two.  This contradiction is independent of local degree and
makes `L_C` empty in the `S6` passport as well.  Every `S6` jump is therefore
exactly `T(2,3)`, `T(3,4)`, or `T(4,5)`, without an exceptional-source caveat.

In either passport the affine finite normalization is now a smooth affine
surface `W`, with `W-D=A2` and `D=A1`.  Miracle flatness makes
`W -> A2` finite flat of rank six; divisor localization gives
`Pic(W)=Z[D]`; and finite-map Riemann--Hurwitz gives `K_W=(e-1)D`.  The
underlying pushforward module is free by Quillen--Suslin, but this supplies no
splitting of the rank-six algebra and no contradiction by itself.

## 2026-07-21: a forced canonical-label corridor

Status: derived here and independently hostile-audited; exact restriction on
Orevkov's original one-dicritical source-blowup graph, not an elimination.

Once `L_C` is empty, Orevkov's finite-chain lemma makes the dicritical `E` a
leaf.  Borisov's logarithmic Keller formula assigns it augmented-canonical
label `e`, not `e-1`.  Writing `m=-E^2`, boundary adjunction gives

```text
label(neighbor(E)) = e*m - 1.
```

The neighbor lies over infinity.  Its positive label rules out type 1, whose
label would be a positive normal degree times the target label `-1`; hence it
is type 2.  Borisov's canonical-label tree theorem then forces the path from
this positive vertex to the original infinity line of label `-2` to contain
consecutive labels `1--0--(-1)`.  Thus the `S6` neighbor label is positive odd,
while the `A6` neighbor label is positive and congruent to two modulo three.

The self-intersection `m` is not bounded by this argument.  The leaf equation
also belongs to the original pure-blowup resolution: later target blowups and
re-resolution can change the strict transform's self-intersection and
adjacency, although its divisorial canonical label remains invariant.

## 2026-07-21: joint minimality fixes the dicritical self-intersection

Status: derived here and independently hostile-audited; a without-loss-of-
generality sharpening of the canonical corridor, not an elimination.

Simultaneously contract any boundary `(-1)`-curve which is exceptional for the
source blowdown and contracted by the resolved Keller map.  Such a contraction
preserves the smooth SNC source boundary and both maps descend; iteration
terminates by Picard rank.  In the resulting jointly minimal model, the
dicritical leaf cannot have been born at a boundary crossing, and a later
touch of the leaf would leave a positive-label terminal `(-1)`-curve of type 2,
contradicting joint minimality.  Therefore

```text
E^2 = -1,             label(neighbor(E)) = e-1.
```

The minimal leaf edge is consequently `2--1` in the `S6` passport and `3--2`
in the `A6` passport.  A corner blowup changes `m` to `m+1` and the adjacent
label from `e*m-1` to `e*(m+1)-1`, proving that the earlier unbounded
progressions were nonminimal blowup orbits rather than genuine graph moduli.
This fixes the first edge but still gives no canonical-label contradiction.

Determinant labels sharpen the endpoint to `d_E=d_A-1<0`, with `d_E` odd for
`S6` and even for `A6`.  An exact hostile blowup family nevertheless realizes

```text
S6: d_E = -(n+1)(n+2)-1,
A6: d_E = -(n+1)(n+2)-2
```

for every `n>=0`, with no positive nondicritical `(-1)`-curve.  The checker
`scripts/canonical_leaf_graph.py` rebuilds the intersection matrices and
determinants.  Therefore the next graph attack needs degree-six pullback,
type-1 covering, or Belyi data; canonical and determinant labels alone still
leave an infinite hostile family.

## 2026-07-21: hostile surface models isolate the finite-algebra gap

Status: explicit consistency models; they are not finite covers or Jacobian
counterexamples.

Complements of smooth ample sections `S=C0+aF` on Hirzebruch surfaces give
smooth affine pairs with `W-D=A2`, `D=A1`, `Pic(W)=Z[D]`, and

```text
K_W = (2*a-d-2)D.
```

The choices `(d,a)=(1,2)` and `(0,2)` realize `K_W=D` and `K_W=2D`, exactly
the one-dicritical `S6` and `A6` packages.  A separate construction on the
triple blowup of `P1 x P1` realizes two disjoint affine-line boundaries with
`Pic(W)=Z[D2]+Z[D3]` and `K_W=D2+2D3`.  Its boundary supports an explicit
ample divisor of square `94`, so the complement is affine.

These examples carry no finite flat degree-six map, Keller coordinates,
branch curves, or monodromy.  They close off a tempting but invalid strategy:
smoothness, boundary topology, the Picard group, and the canonical class alone
cannot eliminate the remaining passports.  Any successful next step must use
the finite rank-six algebra and its trace, different, discriminant, or
compatible monodromy.

## 2026-07-21: the inverse different is nontrivial and nonmonogenic

Status: exact finite-duality reduction; excludes simple global algebra
presentations but does not eliminate a passport.

For `R=C[P,Q]` and `A=O(W)`, finite duality identifies

```text
Hom_R(A,R) = O_W(R_ram).
```

The three ramification divisors are `D`, `2D`, and `D2+2D3`.  Their classes
are nonzero in the already computed free Picard groups.  Thus the rank-six
finite algebra is Gorenstein but not globally Frobenius, not monogenic, and
not one global square complete intersection over `R`.  This does not conflict
with the morphism being syntomic: the complete-intersection presentations
exist locally, while the nontrivial inverse-different line obstructs their
gluing to one standard global presentation.

Writing `R_ram=2E+O` and setting `M=rho_*O_W(E)` gives free rank-six lattices
and an exact trace factorization

```text
T = Phi^* H Phi.
```

Here `H` is symmetric, `det(Phi)` is the equation of the doubled branch part,
and `det(H)` is the equation of the odd branch part.  The exact discriminants
are `b`, `b^2`, and `b2*b3^2`.  The cokernels are pushforwards of boundary
normal lines, fixing the matrix coranks at every known cusp and collision.
The `A6` square-root matrix has corank two at both its `T(2,5)` cusp and its
separate `3+3` collision; the one-boundary `S6` trace matrix has coranks two,
three, and four at its three torus-knot jumps.  These ranks remain compatible.

Ordinary Euler characteristic also closes to an identity: it gives `2` in
both one-boundary cases and `3` in the saturated case, exactly matching the
surface decomposition.  The surviving target is therefore the multiplication
law plus compactified symmetric matrix data, not a scalar Euler deficit.

## 2026-07-21: Orevkov's classified germ realizes the forced `A6` cusp

Status: source-classified and exactly checked hostile local model; not a
Keller map or a global rank-six cover.

Orevkov's 2026 classification applies to the surviving degree-five local
block.  The case-(b) parameters `(k1,k2;l1,l2)=(2,1;2,0)` give exactly
`(d1,d2,N,n)=(2,5,5,3)`.  His derivative normal form integrates to

```text
u = (3*x^5 - 10*y*x^3 + 15*y^2*x)/8,
v = y,
J = (15/8)*(x^2-y)^2.
```

The critical parabola maps as `x -> (x^5,x^2)`.  The complete pullback is

```text
64*(u^2-v^5) = (x^2-y)^3*(9*x^4-33*x^2*y+64*y^2).
```

The residual factor splits into two smooth parabolas, and exact reduction in
their quadratic parameter proves that each maps with normalization degree
one to the cusp.  Together with the critical component, this verifies the
componentwise-bijective hypothesis rather than silently assuming it.  The
quintic fiber discriminant is exactly
`1,036,800,000*(u^2-v^5)^2`, so the required `A6` square discriminant also
survives.  The new typed checker and hostile perturbation test replay every
identity.

This is a sharp negative research result: the cusp forced by the global
degree-six argument is an actual classified finite-map germ.  Any elimination
must couple it to the separate `3+3` collision or to global multiplication,
trace, compactification, or splice data.

## 2026-07-21: two-pair infinity stratum eliminated

Status: derived and independently hostile-audited conditional theorem; not a
full elimination of either passport.

Orevkov's Section 2.4 splice equations apply when the one-dicritical target
branch has exactly two characteristic pairs at infinity and condition `(*)`
holds: some target boundary component has exactly one noncontracted
irreducible preimage.  In his notation the surviving passports have
`m=1`, `n=3` or `2`, and `N=6`.  Equations (11), (14), and (15) give

```text
N = m^2*x*d1*d2^2*R1*S1,
x = a+n,
R1 >= 2.
```

For `A6`, this forces `N>=8`.  For `S6`, equality first forces
`x=3, R1=2, d1=d2=S1=1`.  Equations (12)--(13) then determine the two boundary
stars, and the only surviving arm count reaches equation (10) as

```text
3*Q_tilde = 14,
```

impossible for an integral edge determinant.  The preceding Belyi passport
`(2,2,2),(5,1),(4,1,1)` is itself realizable, so the splice determinant is the
decisive input.  A typed exact checker replays the bound, both possible second
arm counts, the fractional determinant, and a transitive product-one
permutation triple.

The conclusion cannot be promoted yet: neither exactly two infinity
characteristic pairs nor condition `(*)` follows from the current smoothness,
Picard, monodromy, trace, or minimal-leaf results.

## 2026-07-21: exact split transport does not force condition `(*)`

Status: derived matrix theorem plus exact hostile blowup lattices; no
compactified map.

For a target boundary component with several noncontracted source preimages,
let `n` be their normal multiplicities, `m` their tangential degrees, `q` the
target dual square, and `Q` the selected block of the inverse source
intersection matrix.  Projection formula replaces Orevkov's scalar
determinant ratio by

```text
Q*m=q*n,   n.m=N,   m^t*Q*m=N*q.
```

For `q>0`, Hodge index gives the sharp inequality

```text
Q-(q/N)*n*n^t <= 0.
```

Orevkov's negative-intersection determinant convention reverses this sign:
his positive first-star determinant `R1*D1` means raw dual square
`q=-R1*D1<0`.  Thus the Hodge inequality is unavailable at the star used in
the conditional proof; the projection-formula matrix identity is the robust
part.

The scalar formula is recovered when there is one noncontracted preimage, but
the matrix version has genuine split solutions at degree six.  Three explicit
boundary blowups give an integral dual block `[[1,2],[2,1]]` with
`m=(1,1), n=(3,3)`, hence `A6`-compatible `3+3` inertia.  Four blowups give
`[[-2,2],[2,0]]` with `m=(1,2), n=(2,2)`, hence `S6`-compatible `2+2+2`
inertia.  Both full boundary intersection matrices are unimodular, and both
Hodge remainders are negative semidefinite with the required kernel.

A larger `A6` stopping fixture survives more structure.  Its selected target
component has valency two, augmented-canonical label `-1`, and dual square
`q=2`.  Two disjoint valency-two source components have labels `(-5,-1)`,
normal degrees `(5,1)`, tangential degrees `(1,1)`, and dual block
`[[6,4],[4,-2]]`.  Thus the local type-one canonical pullback rule holds
exactly, the monomial germs `(x^5,y)` and `(x,y)` supply local models, and a
5-cycle together with `(456)` generates all `A6`.

These models do not include the dicritical, canonical labels, or local
holomorphic maps.  They prove that determinant transport, Hodge index,
integrality, degree, and inertia parity alone cannot force condition `(*)`.
The next argument must couple the split matrix to canonical-discrepancy or
local-map data across the entire boundary, not just at one selected component.

## 2026-07-21: a genuine finite `S6` cover globally violates `(*)`

Status: exact finite flat hostile cover; not Keller and not the surviving
branch passport.

The affine map

```text
(x,z) -> (P,Q)=(x^2+z, x^3-x^4*z)
```

has coordinate algebra

```text
C[x,z] = C[P,Q][x]/(x^6-P*x^4+x^3-Q),
```

so it is finite flat of rank six.  Three basepoint-free sections of
`O_P1xP1(3,1)` extend it to a finite degree-six morphism
`P1 x P1 -> P2`.  The target boundary has one line `T`, and

```text
phi^*T=2*C1+C2,
```

with both source components noncontracted of tangential degrees `(1,4)`.
Thus condition `(*)` fails for the complete boundary of an actual finite
cover, not merely for an abstract matrix.

The six-sheet polynomial is indecomposable over `C(P)`, so its polynomial
monodromy is primitive.  Its exact discriminant is `Q^2*H(P,Q)`, with `H`
squarefree in `Q`; generic inertia along `H=0` is a transposition.  Hence the
geometric monodromy is all of `S6`.  After one corner and one free source
blowup, the boundary even has numerical labels `2--1--0--(-1)` and determinant
labels `-3,-2` on the first edge, while retaining exact pullback transport.

The decisive missing condition is visible in the Jacobian

```text
x^2*(4*x*z-2*x^3-3),
```

which is not constant.  The cover has affine branching and its numerical
label-two leaf is contracted to a boundary point, not a type-three
dicritical.  Therefore the model does not threaten the conjecture; it proves
that any derivation of `(*)` must use affine unramifiedness, the genuine
Keller type assignment, or the forced branch passport.

## 2026-07-21: exact matrices close the trace-only route

Status: exact hostile consistency models; no finite algebra or Keller map.

For an `A6` rational branch with parametrization
`(P,Q)=(t^2+t^3,t^5)`, an explicit symmetric `3 x 3` Bezout matrix `K` has
determinant the branch equation and cokernel the normalization module.  Setting

```text
Phi = diag(K,I3),
H   = [[0,I3],[I3,0]],
T   = Phi^t H Phi
```

gives `det(Phi)=b`, `det(H)=-1`, `det(T)=-b^2`, exact generic and special
coranks, and a primitive vector of norm six.  The curve has the forced
`T(2,5)` cusp and smooth two-preimage collisions.

A separate symmetric `5 x 5` Bezout matrix for

```text
P=t^4-6*t^2,
Q=t^5-5*t^3
```

realizes the full `S6` `1+1+1` jump partition: three `T(2,3)` cusps, three
normalization nodes, exact normalization-module coranks, and the unit line.
Its homogenization is compatible with the rational normalization's theta
characteristic `O_P1(-1)`, as predicted by Beauville's symmetric-resolution
theory.

The new exact checker replays all determinants, parametrized kernel vectors,
ranks, and unit norms.  The models omit the essential datum: a unital
commutative associative rank-six multiplication whose regular trace is the
displayed matrix.  The next algebraic attack is therefore the integral
symmetric cubic trace tensor, not another quadratic determinant invariant.

## 2026-07-21: global topology excludes the hostile `S6` trace curve

Status: exact exclusion of one explicit curve; the `S6` passport remains.

The trace curve

```text
P=t^4-6*t^2,
Q=t^5-5*t^3
```

has a degree-four projection to the `P`-axis.  Zariski--van Kamp therefore
generates its affine complement with four fiber meridians.  If all four map
to transpositions of six sheets, their edge graph has only four edges and
cannot be connected.  The generated group is consequently intransitive.

An independent exact replay stores Sage's six-relator, four-generator
presentation and checks all `15^4=50625` transposition assignments.  Exactly
`735` satisfy the relations: `15` generate order two and `720` generate order
`120`; none is transitive.  Thus the quadratic trace matrix remains a valid
consistency model for its advertised invariants, but its branch curve cannot
underlie a connected `S6` cover.

This does not eliminate the one-dicritical `S6` frontier.  It adds a necessary
global constraint missing from the trace lattice: a one-component
transposition branch must have affine projection width at least five.  The
next target is the smallest width-five curve compatible with the forced cusp
and normalization-collision budget.

## 2026-07-21: the smallest width-five `S6` near-miss is cyclic

Status: exact geometric fixture and computer-assisted global exclusion; the
`S6` passport remains.

The curve

```text
(P,Q)=(t^5+t^4,t^7+t^5)
```

attains the transposition-generator lower bound `width=5`.  Exact
normalization arithmetic finds the complete singularity package: one finite
`T(4,5)` jump-three cusp, six transverse nodes, and one `(2,7)` pair at
infinity.  Their deltas are `6+6+3=15`, equal to the arithmetic genus of the
degree-seven closure.  Thus this lies in the degree-minimal
singular-one-pair slot, not a fixture that failed a local budget.

It fails globally.  Sage's five-meridian Zariski--van Kamp presentation has an
explicit Tietze isomorphism to `Z`, sending all five geometric meridians to one
generator.  The dependency-free checker independently exhausts all
`15^5=759375` transposition assignments: exactly 15 satisfy the ten relators,
all have order-two image, and none is transitive.  A connected six-sheet cover
cannot use this branch curve.

This does not classify all width-five curves.  The remaining `S6` target is a
different width-five curve with noncyclic complement, or a higher-width or
multi-pair infinity configuration, coupled to the Keller compactification.

## 2026-07-21: one-pair `S6` infinity starts at `(5,12)`

Status: exact conditional stratum elimination and sharp group witness; not a
solution of `JC(2)`.

Assume the polynomial normalization degrees are `m<n` and there is exactly
one genuine singular pair at infinity.  The projective pair is `(n-m,n)`,
but the link on a large affine sphere is `T(m,n)`.  Leidy--Maxim's
at-infinity epimorphism sends this torus-knot group onto the global affine
complement.  In an `S6` quotient the common central power in

```text
<a,b | a^m=b^n>
```

dies.  A Bezout meridian `a^u*b^v`, with `n*u+m*v=1`, must map to a
transposition.

Width gives `m>=5`; a genuine singular pair gives `n-m>=2` and
`gcd(m,n)=1`.  There are ten candidate degree pairs with `n<=11`.  Exact
enumeration of all ordered `S6` images satisfying the two power constraints
and the meridian constraint finds no generating pair.  The sole substantial
near-miss is `T(5,8)`: its `735` compatible pairs give `15` groups of order
two and `720` groups of order `120`, never `S6`.  Hand analysis explains the
same result.  A transitive image would make `a` a 5-cycle fixing one sheet.
If the meridian avoids that sheet the image fixes it; if it contains that
sheet, the forced image of `b` is a 6-cycle and violates `b^8=1`.

The obstruction is sharp.  For `(m,n)=(5,12)`, the census gives

```text
compatible pairs: 2175
generated orders: 2:15, 120:1440, 720:720
```

One witness is

```text
a=(2 3 4 5 6), b=(1 2 5 3 6 4), a^-2*b^5=(1 2).
```

Conjugating `(1 2)` by the powers of `a` produces all five star
transpositions, proving that the image is `S6`.  Hence `(5,12)` is the first
one-pair infinity-group target.  This does not construct its branch curve or
cover, and it does not touch multi-pair infinity.

## 2026-07-21: the coarse one-pair `A6` audit first reaches delta five

Status: conditional, computer-assisted exclusion; not a solution of `JC(2)`.

Assume polynomial normalization degrees `a<d`, one genuine projective pair
`(d-a,d)` at infinity, the forced finite `T(2,5)` cusp, and otherwise only
smooth normalization collisions of total delta `Delta`.  Rational genus gives

```text
(d-1)*(a-1)=2*(Delta+2).
```

For `Delta=1`, the sole singular-infinity candidate is `(a,d)=(2,7)`.
The large affine link is `T(2,7)`, not the small projective `T(5,7)` knot;
an exact `A6` quotient census finds no generating pair with a single-3-cycle
meridian.  For `Delta=2`, `(2,9)` fails the same census, while `(3,5)` has
`720` coarse link quotients and requires a global complement calculation.

Every remaining `(3,5)` parametrization is affine-equivalent to

```text
(P,Q)=(t^2+t^3,c*t^4+t^5).
```

Exact collision elimination finds the only exceptional parameters
`1/2,5/6,1,5`.  The first three violate the forced singularity/fiber
hypotheses; `c=5` is the sole valid contact-two fiber.  Over the complement of
those values, the properly compactified family is Whitney equisingular, so
Thom isotopy transports the complement group from `c=0`, where Sage's exact
van Kamp presentation simplifies to `Z`.  At `c=5`, all `40^3` assignments of
the three geometric meridians to 3-cycles give only `40` cyclic and `720`
`A5` images, never `A6`.

The same genus and large-link census then excludes `Delta=3`: its only degree
pair is `(2,11)`.  At `Delta=4`, both `(2,13)` and `(3,7)` fail.  At
`Delta=5`, `(2,15)` has `2880` `A6`-generating torus pairs but none with the
required single-3-cycle meridian; `(3,8)` has exactly `720` correct pairs and
is the first coarse survivor.  At this stage, therefore, `Delta>=5`, with
`(3,8)` forced at equality under the stated hypotheses.  The next log entry
records the residual audit that eliminates this coarse survivor.

This also kills
the explicit hostile `A6` trace curve globally while preserving its original
purpose as a quadratic trace-lattice consistency model.  The remaining
unrestricted problems are to force the one-pair/finite-singularity hypotheses
from the Keller compactification and to handle multi-pair infinity.

## 2026-07-21: the conditional `A6` delta-five family is exhausted

Status: conditional exact family algebra and finite permutation replay, plus
computer-assisted Zariski--van Kamp extraction and equisingular topology.

After a polynomial target shear, every conditional `(3,8)` curve has the
three-parameter normal form

```text
P=t^2+t^3,
Q=alpha*t^5+beta*t^6+gamma*t^7+t^8.
```

The collision quintic `H(s)` gives exact equations for every degeneration.
The generic open has the forced cusp and five nodes.  Its two admissible
codimension-one walls are an irreducible ordinary-triple divisor `T=0` and an
irreducible contact-two divisor `G=0`; the other factors either destroy the
`T(2,5)` cusp, create extra intrinsic ramification, or add a smooth preimage to
the cusp image.

A clean generic model `(alpha,beta,gamma)=(1,1,0)` has exactly one `T(2,5)`
cusp, five transverse nodes, a `(5,8)` infinity pair, and genus accounting
`2+5+14=21`.  Sage nevertheless gives affine complement `Z`.  Representatives
of the ordinary-triple and contact-two walls also have complement `Z`.  For
each of the three exact presentations, all `40^3` single-3-cycle assignments
reduce to 40 cyclic images and never `A6`.

Plane-curve equisingularity and proper Thom isotopy propagate the three
representative computations across the generic open and both generic walls.
Every valid survivor is therefore confined to
`(Sing(G) union (G intersect T))` away from the explicitly invalid factors.
Exact primary decomposition splits this codimension-two locus into four valid
rational curves:

- `A`, with a generic contact-three collision and two nodes;
- `B`, with two generic contact-two collisions and one node;
- `N`, with a generic nonordinary triple and one node; and
- `E`, with a generic ordinary triple and a separate contact-two collision.

The old description of `(alpha,beta,gamma)=(20/3,1,2/3)` as an ordinary
triple plus a separate tangency was wrong.  It is the `u=-3` representative
of `N`, hence a nonordinary triple with one tangent branch pair plus a node.

Exact factor identities leave only four valid exceptional parameter points:
one contact-four point `P4`, one higher nonordinary triple `PNT`, and the two
conjugate contact-three/contact-two points `P32+/-` over `QQ(sqrt(6))`.  Sage
produces presentations for the four generic residual types and all four
exceptional points.  Each dependency-free `40^3` replay has exactly 40
satisfying assignments, all with cyclic image `C3`, and none with image `A6`.
The checker also verifies that `A,B` lie on `Sing(G)`, that `N,E` lie on
`(G,T)`, and that every exceptional factor maps either to an invalid divisor
or one of the four audited points.

Thus collision delta five is impossible under the one-pair and
finite-singularity hypotheses.  Delta six has only `(2,17)`, whose link
census fails.  Delta seven has `(2,19)`, `(3,10)`, and `(4,7)`; the exact
single-3-cycle `A6` counts are `0`, `720`, and `0`.  The sharpened conditional
frontier is therefore `Delta>=7`, with equality forcing affine degrees
`(3,10)` and projective infinity pair `(7,10)` at the coarse link stage.

This does not prove the one-pair or finite-singularity hypotheses, construct
the surviving `(3,10)` branch, eliminate the unrestricted `A6` passport, or
prove the plane Jacobian conjecture.  The algebraic identities and finite
group enumerations are exact; deriving the stored presentations and
transporting representatives along strata remain the computer-assisted
Zariski--van Kamp and Whitney--Thom layers.

## 2026-07-21: the generic conditional `A6` delta-seven family is cyclic

Status: exact four-parameter family algebra and finite permutation replay,
plus computer-assisted Zariski--van Kamp extraction and equisingular topology;
the collision walls remain open.

At the next conditional equality, target shears give the complete normal form

```text
P=t^2+t^3,
Q=alpha*t^5+beta*t^7+gamma*t^8+delta*t^9+t^10,
alpha != 0.
```

The collision polynomial is septic.  Its discriminant is `C*G`, the endpoint
resultant is `-alpha*L`, the tangent resultant is `-alpha*L*G`, and the
triple-image resultant is `C^2*T`.  Thus the fully nondegenerate open is
`alpha*L*C*T*G != 0`.

The exact point `(alpha,beta,gamma,delta)=(2,0,0,0)` gives
`Q=2*t^5+t^10`.  It has one `T(2,5)` cusp, seven distinct transverse nodes,
and a `(7,10)` branch at infinity, with complete genus accounting
`36=2+7+27`.  Sage simplifies its affine complement to `Z`.  The independent
`40^3` replay finds exactly 40 satisfying single-3-cycle assignments, all
with image `C3`, and no `A6` image.

The nondegenerate open is connected.  After a finite base change labels the
seven node sections, proper projective Whitney--Thom transport propagates the
cyclic complement throughout it.  Consequently any conditional delta-seven
survivor lies on `G=0` or `T=0`.  This does not classify those walls, derive
the one-pair and finite-singularity hypotheses, construct an `A6` cover or
Keller map, or prove the plane Jacobian conjecture.

## 2026-07-21: the forced `A6` relations survive the spin cover

Status: exact hostile lift; no compactification or Keller map.

The natural six-sheet action lifts through the spin double cover `2.A6`.
Using exact rational Clifford arithmetic, the canonical order-three lifts of

```text
r=(345), s=(123), b=(456)
```

satisfy both forced finite relations: the `T(2,5)` five-braid relation lifts
with no central error, and the disjoint `3+3` collision meridians commute.
The lifts generate all `720` elements of `2.A6`; in particular this is not a
splitting of `A6`.

The cusp preferred longitude does retain a central minus sign.  That sign is
not yet obstructive because an even-order infinity lift can absorb it.  A
stronger hostile certificate removes even that parity loophole: two explicit
transitive product-one six-tuples made only of 3-cycles extend the same forced
prefix, one with Fried--Serre spin `+1` and the other with spin `-1`.

Therefore the finite cusp/collision data do not determine a spin component.
Any renewed spin attack must first derive the actual infinity word and its
peripheral framing from the Keller compactification.  The checker enumerates
the upstairs and downstairs groups and replays every Clifford identity; a
sign-flipped order-six lift is an adversarial failure fixture.

## 2026-07-21: local multiplication survives both forced `A6` events

Status: exact local hostile algebras and a finite global equation system; no
global algebra, cover, or Keller map.

For the explicit `A6` trace matrix `T=[[0,K],[K,0]]`, a completely symmetric
cubic tensor has `56` independent polynomial entries.  Multiplication is
recovered using `adj(K)/b`.  Polynomial divisibility, the unit equations,
denominator-cleared WDVV equations, ordinary-regular-trace equations, and the
additional `Phi/H` middle-lattice equations now give an exact finite target
for a global enhancement.

The normalization parametrization `(P,Q)=(t^2+t^3,t^5)` makes the divisibility
conditions concrete.  Exact calculation gives

```text
adj(K)(t) = -t^4*Phi5(t)*(1,t,t^2)^t*(1,t,t^2),
```

and membership in `im(K)` becomes one polynomial identity on the
normalization.  This retains the value-and-jet information that disappears if
one merely specializes the vanishing adjugate at a singular point.

Both tempting pointwise contradictions fail.  An explicit basis identifies
the cusp trace form with the ordinary trace of
`C[z]/(z^5) x C`; another identifies the collision trace form with that of
`C[z]/(z^3) x C[w]/(w^3)`.  Each algebra also has a perfect Frobenius middle
form and a divisor section whose squared factorization is the ordinary trace,
with multiplication-section kernel equal to `ker(Phi)`.  A separate monic
rank-six family over the normalization realizes the transitions
`3+1+1+1 -> 5+1` and `3+1+1+1 -> 3+3`.

These are residue-fiber and normalization-level stopping fixtures.  They do
not solve the global cubic system or descend it across the singular branch.
The remaining obstruction must couple normalization jets, the middle module,
connected monodromy, and infinity globally.

## 2026-07-21: sparse-interior obstruction for `(72,108)`

Status: derived here and exhaustively checked; finite sparsity theorem only.

The separate coordinate-degree frontier below maximum degree `125` reduces,
by Guccione--Guccione--Horruitiner--Valqui, to `(72,108)` and its transpose.
Their Proposition 4.3 gives two transformed Newton-polygon pairs with
`[P,Q]=x^2`.  Every boundary lattice coefficient was allowed.  Case 1
enumerates every exact support with at most two nonzero strict-interior
coefficients; replayed zero-product certificates contradict all `7504`.
Case 2 enumerates every support with at most three: zero-product propagation
contradicts `3678`, and exact Groebner unit-ideal replays contradict the five
exceptional triples, for `3683/3683` total.  Therefore Case 1 needs at least
three nonzero interior coefficients and Case 2 needs at least four.

The full polygons deliberately stop the propagation procedure, as does a
named four-interior Case 2 fixture.  Those hostile checks prevent the sparse
lower bounds from being mislabeled as an elimination of `(72,108)`.

## 2026-07-21: the full conditional `A6` delta-seven family is excluded

Status: exact family and component algebra, exhaustive finite-group replay,
and computer-assisted topology under four standing one-pair hypotheses; not a
proof of the plane Jacobian conjecture.

The degree-`(3,10)` normal form has collision septic `H` in an exact
four-dimensional affine coefficient slice.  The complete root-partition
ledger now exhausts the triple-image wall `T=0`, every
positive-dimensional repeated-collision stratum on `G=0`, every excess-four
endpoint, and all coarser partitions.  The complete trigonal fiber on `T=0`
gives three global generating meridians.  For every contact profile allowed
by the delta budget, all `40^3` single-three-cycle assignments fail to
generate `A6`.

The positive-dimensional component audit is no longer based only on sample
curves.  Constant-rank incidence bundles cover the generic double and triple
strata.  Exact saturation handles the unordered two-double chart and shows
that its denominator branches lie on `C=0` or `L=T=0`.  The contact-four
chart loses no point at its apparent denominator.  The mixed
`(3,2,1,1)` incidence has one geometrically irreducible genus-two component;
its dense rational chart misses two valid conjugate points, but primary
decomposition places both on the same component.  The unrestricted
three-double incidence has four prime supports, exactly one of which survives
away from `C=0` and the already audited `T=0` wall.

The nine finite endpoint embeddings lie in three irreducible arithmetic
orbits of degrees `(2,4,3)`.  A checked presentation for one embedding of each
orbit has only cyclic `C3` images.  Riemann existence converts a hypothetical
finite topological quotient into a finite etale cover; field automorphisms
transport the labeled `A6` cover and divisorial tame inertia.  This excludes
the required quotient at every arithmetic conjugate, including the uncomputed
real `P33` embedding.  It does not assert homeomorphism or equality of the
discrete complement groups.

Consequently, under polynomial normalization, exactly one genuine infinity
pair, only the forced `T(2,5)` intrinsic cusp, and smooth-branch collisions for
all other finite singularities, collision delta seven is impossible.  The
exact link census then excludes deltas eight and nine.  At delta ten, only
affine degrees `(4,9)` survive the coarse single-three-cycle `A6` test, with
projective infinity pair `(5,9)` and exactly 720 qualifying pairs.  Thus the
combined conditional frontier is

```text
Delta >= 10,
Delta = 10  =>  (a,d) = (4,9),  (d-a,d) = (5,9).
```

The algebraic identities, primary decompositions, endpoint orbit arithmetic,
and permutation censuses are exactly reproducible.  Dependency-free tests
replay every stored relation word.  Checked Sage scripts regenerate the
generic-open and three excess-three presentations and certify all incidence
components; the triple-wall, excess-one/two, and finite-endpoint van Kamp
extractions remain trusted computer-assisted inputs.  Family-wide transport
also depends on proper Whitney--Thom triviality, and the endpoint argument on
Riemann existence and tame-inertia compatibility.

No cover or Keller map is constructed at `(4,9)`.  The next two checkpoints
exclude the clean branch-curve family and the generic opens of its two dominant
degeneration walls; lower-dimensional strata remain.  The four standing
hypotheses are not derived for an arbitrary Keller branch, and the unrestricted
`A6` passport, the `S6` passport, and therefore `JC(2)` remain open.

## 2026-07-21: the clean conditional `A6` delta-ten family is excluded

Status: exact family algebra, exhaustive finite-group replay, and
computer-assisted topology under the same four standing hypotheses; not a
proof of the plane Jacobian conjecture.

Every degree-`(4,9)` normalization is polynomially equivalent over `C` to

```text
P = t^2 + k*t^3 + t^4,
Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9,  a != 0,
```

with residual involution `(k,a,b,c,d) -> (-k,a,-b,c,-d)`.  Exact unordered-pair
elimination gives a monic collision decic.  The apparent denominator fibers
`k=0,+2,-2` are legitimate reducible incidence fibers with degree splits
`2+8`, `4+6`, and `4+6`; they are not silently discarded.

For `(k,a,b,c,d)=(1,1,0,0,0)`, the curve has one `T(2,5)` cusp, ten reduced
nodes, one `T(5,9)` branch at infinity, and genus balance `28=2+10+16`.
Sage regenerates its four-generator affine van Kamp presentation and
simplifies it to `Z`.  The dependency-free replay checks all
`40^4=2,560,000` single-three-cycle assignments: exactly 40 survive, all
diagonal with image `C3`, and none generates `A6`.  The clean parameter locus
is a nonempty Zariski-open subset of irreducible affine five-space, hence
connected; proper projective Whitney--Thom triviality propagates the complement
result through that locus.  Thus every remaining conditional delta-ten
candidate lies on a degeneration wall.

An independent hostile census prevents a false shortcut.  There are exactly
720 qualifying `T(4,9) -> A6` pairs, forming two inner orbits of 360 fused by
odd sheet relabeling.  Every pair has a unique compatible forced
`T(2,5)`-cusp plus collision decomposition, and all 720 admit exact compatible
spin lifts in `2.A6`.  Therefore the local, peripheral, orientation, and spin
constraints do not eliminate the frontier; global complement topology is
essential.

## 2026-07-21: both dominant conditional delta-ten walls are excluded

Status: exact incidence and component algebra, exhaustive finite-group replay,
and computer-assisted topology under the same four standing hypotheses; not
a complete wall audit and not a proof of the plane Jacobian conjecture.

The valid contact incidence is `H=T=0` over pair-sum coordinates `(k,s)`.
The gcd of the six coefficient minors in `(a,b,c,d)` is exactly

```text
s*(k+2*s)^3*(2*s^2+3*k*s+4).
```

After removing the denominator and diagonal factors, a Rabinowitsch
localization of the residual minors is the unit ideal.  Thus the incidence is
an affine two-plane bundle over an irreducible two-dimensional base and has
one irreducible four-dimensional dominant image.  Every raw rank-drop
component is separately identified as `a=0`, an extra critical point, a cusp
image collision, or a `k=0,+2,-2` split-chart artifact.

For ordinary triples, the elementary symmetric parameters of three source
points lie on the irreducible surface

```text
sigma2^2 - sigma1*sigma3 - sigma2 = 0.
```

Both the four-column coefficient matrix and its five-column augmentation have
localized minor ideal equal to one on the valid distinct-source open.  The
ordinary-triple incidence is therefore another irreducible four-dimensional
dominant component.

Exact representatives on both components were certified.  The contact member
has one `A3` contact and eight nodes; the triple member has one ordinary triple
and seven nodes.  Both retain the forced `T(2,5)` cusp and `T(5,9)` infinity
branch, and both close the genus budget at 28.  Checked Sage scripts
independently regenerate their primitive implicit curves, primary
decompositions, and four-generator raw van Kamp presentations.  Both
complements simplify to `Z`; both exhaustive `40^4` replays leave exactly 40
diagonal `C3` images and zero `A6` images.  Proper projective Whitney--Thom
transport excludes the connected generic equisingular open of each divisor.

The complete finite collision ledger has 35 non-node local contact-tree atoms
and 145 global profiles.  Exactly two have expected codimension one—the two
divisors just excluded.  The next layer has six expected codimension-two
profiles; there are 38 expected endpoints and 55 overdetermined candidates
that still require exact elimination.  Expected dimension is explicitly not
used as an emptiness proof.

## 2026-07-21: six codimension-two open-chart components are identified

Status: exact component algebra, exhaustive finite-group replay, manual Sage
singular schemes and presentations, and computer-assisted topology under the
same four standing hypotheses.  All six displayed generic or dominant
components are excluded at the stated level.  The `T112` and mixed
contact-plus-triple topology steps are proved over their smooth labeled
incidence spaces by finite-etale labeling, relative embedded resolution, and
proper Whitney--Thom isotopy.
Residual-rank, split-chart, removed pair or `P`-projection/critical-fiber, and
deeper-intersection loci remain open; a removed `P`-critical fiber is not
automatically a singular curve branch when `Q'` is nonzero.  This is not a
proof of the plane Jacobian conjecture.

The six expected codimension-two profiles were attacked by their actual
target-fiber incidences rather than only by roots of the collision decic:

- `C3+7N`: `H=H'=H''=0` has rank three on the valid pair chart.  Its only
  localized rank-drop base is inconsistent by a nonzero augmented minor, so
  the valid incidence has one irreducible three-dimensional image.
- `C2^2+6N`: the ordered equations at two distinct double roots give a
  rational Cramer threefold dominating the full `(k,u,v)` ordered base and
  invariant under the free root swap.  The determinant/augmentation gcd
  removes the residual rank factor as another dominant divisor.  Compatible
  residual-rank subvarieties and the split `k=0,+2,-2` intersections could
  support other components and remain open.
- `T112+6N`: the labeled triple-root base
  `sigma2^2-sigma1*sigma3-sigma2=0` is absolutely irreducible and smooth on
  the displayed `P`-unramified localization.  The removed
  triple-source/fourth-root loci are `P`-critical-fiber boundaries, not
  automatically singular branches of `(P,Q)`; the full clean profile there
  is not covered.  Sage saturation makes the maximal-minor ideal the unit
  ideal on the displayed chart.  The sample projection fiber is reduced of
  length two, proving generic finiteness.  Degree two additionally uses the
  geometric fact that the clean target has a unique `T112` fiber and hence
  only the two orientations of its tangent pair.  That optional degree-two
  statement is not needed for topology propagation: the labeled incidence is
  smooth and irreducible, and its nonempty clean open has a simultaneous
  embedded resolution after finite-etale labeling.  The generic chart is
  excluded.
- `C2+T111+5N`: the combined contact and triple equations define a rational
  Cramer threefold on the displayed dense open.  The sole nonboundary rank
  factor is irreducible over `QQ` and coprime to a specialized augmented
  minor, so it has no additional component dominating that divisor.  The
  cyclic sample is exact.  Its smooth irreducible clean Cramer open has the
  same relative-resolution and proper-isotopy transport, so the generic chart
  is excluded.  Compatible lower-dimensional rank, split,
  overlap, and removed triple/fourth-root `P`-critical loci remain open; the
  last are not automatically singular branches when `Q'` is nonzero.
- `Q0+4N`: constancy of `Q mod (P-h)` gives three affine-linear equations of
  rank three for `h!=0`; the first equation recovers `h` uniquely.  The valid
  incidence therefore projects to one irreducible threefold.
- `T111^2+4N`: two omitted roots give a dense rational Cramer graph.  The
  same-fiber determinant factor is invalid for separate target triples, while
  the residual factor `k*u*v-u-v` is generically inconsistent.  Compatible
  pieces supported on that factor and split boundaries remain open.

For every profile, an exact rational representative has the forced
`T(2,5)` cusp, collision delta ten, the fixed `T(5,9)` infinity branch, and
genus balance 28.  Manual Sage 10.8 checkers regenerate each primitive curve,
singular scheme, raw four-generator van Kamp presentation, and a
simplification to `Z`.  Every dependency-free `40^4` replay leaves exactly
40 diagonal `C3` images and no `A6` image.  The repository CI runs the Python
certificates and tests only, not Sage.  Neither CI nor the Sage checkers
certify proper projective Whitney--Thom transport.

This maps the displayed dimension-three open-chart layer, not its full
closure.  The next exact obligations are the three residual coefficient-rank loci,
all six true split-incidence charts, removed pair and
`P`-projection/critical-fiber intersections, and the fourteen expected
codimension-three profiles.  The detailed statement and commands are in
`docs/a6-delta-ten-codim-two.md`.

## Negative results and guarded boundaries

- Freezing `z` and selecting two outputs does not inherit a constant plane
  Jacobian. There is no direct dimensional descent.
- A low-degree brute-force search would not address the known plane frontier.
- The three-dimensional example does not justify saying that `JC(2)` is false.
- Finite separability does not justify calling the plane function-field
  extension Galois; normality is a separate condition.
- The degree-six frontier has been reduced to `A6` and `S6` monodromy, not
  solved.
- Same-day social posts are not enough to settle historical attribution, even
  though the displayed finite certificate is fully checkable.

## Next research obligations

1. Attack the surviving `A6` and `S6` generic-degree-six profiles with
   compactification, splice, adjunction, and Belyi constraints.  In particular,
   try to force Orevkov's condition `(*)` or replace it with a determinant
   inequality stable under several noncontracted preimages; do not relabel the
   conditional two-pair theorem as a solution of the unrestricted conjecture.
2. Continue the conditional delta-ten elimination with the residual-rank and
   split-chart closures from `docs/a6-delta-ten-codim-two.md`, then audit the
   fourteen expected codimension-three profiles without using expected
   dimension as emptiness.
3. Formalize more of the projective simple-root/fiber correspondence if a
   useful reusable algebraic-geometry interface is available in mathlib.
4. Extend the exact support analysis beyond two interior coefficients in the
   residual `(72,108)` and `(108,72)` cases, using elimination only after the
   zero-product core is exhausted.
5. Keep every claim tagged as certified, derived, announced, experimental, or
   open.
