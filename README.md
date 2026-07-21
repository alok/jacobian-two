# Jacobian Two

An exact, Lean-first audit of the new three-dimensional Jacobian
counterexample, together with a deliberately honest attack on the surviving
plane conjecture.

## Short answer: the plane case is still open

Galois theory proves a conditional theorem: a Keller map whose induced
function-field extension is **already Galois** is a polynomial automorphism.
For a plane map `F=(P,Q)`, a nonzero constant Jacobian makes

\[
  k(P,Q)\subset k(x,y)
\]

finite and separable, but it does not make the extension normal.  That missing
normality hypothesis is the gap in the tempting “Galois theory solves the
plane” argument.

Here “degree” means generic fiber size, or equivalently the degree of this
function-field extension.  Degrees through five are excluded for a
hypothetical plane counterexample; degree six is the first unresolved
small-sheet case in the accepted literature.  The detailed audit and primary
sources are in [`docs/galois-frontier.md`](docs/galois-frontier.md).

The resulting dimensional boundary is now:

\[
  \mathrm{JC}(1)\text{ is true},\qquad
  \mathrm{JC}(2)\text{ remains open},\qquad
  \mathrm{JC}(n)\text{ is false for }n\ge 3.
\]

This repository does **not** claim to solve `JC(2)` or its degree-six frontier.
It records exactly what was proved, what was derived, and what remains open.

## New progress at the first open sheet degree

The degree-six frontier can now be narrowed substantially.  Let `G <= S6` be
the monodromy group of a hypothetical six-sheet plane Keller counterexample.
Affine local inverse branches force every branch meridian to fix a sheet, and
those meridians normally generate `G`.  Exact enumeration of all sixteen
transitive degree-six groups initially leaves seven possibilities.  Orevkov's
exact defect budget, the finite-normalization boundary, deck symmetry, purity
of the branch locus, and a local double-transposition obstruction eliminate
the three imprimitive survivors.  A new Riemann--Hurwitz refinement then
charges a dicritical of normal index `e` and tangential degree `d` at least
`e*d` units in Orevkov's budget:

\[
  \sum_E(e_Ed_E+\delta_E)=5,\qquad \delta_E\ge0.
\]

Together with an irreducible-branch normalization obstruction, this also
eliminates `A5` and `S5`.  Therefore

\[
  \boxed{G\in\{6T15=A_6,\;6T16=S_6\}}.
\]

Equivalently, the function-field extension of any hypothetical six-sheet
counterexample must be primitive.  This applies without a one-dicritical
hypothesis.  More precisely, the unrestricted ramified profiles are now:

- `A6`: one `(e,d)=(3,1)` branch; its normalization is noninjective, and every
  normalization collision is two smooth multiplicity-three branches which
  use all six sheets and give an omitted affine value;
- `S6`: one transposition branch, two distinct transposition branches, or
  distinct transposition and 3-cycle branches.  In the two-transposition
  profile at least one zero-jump branch curve must self-collide.  In the
  saturated transposition-plus-3-cycle profile both curves self-collide, with
  exact fiber rows `3+3`, `2+2+1+1` or `2+2+2`, and `2+3+1` at a
  cross-intersection.  Its two minimal constant chains are also impossible:
  the finite normalization is smooth and finite flat, with disjoint boundary
  `D2,D3=A1`, `Pic(W)=Z[D2]+Z[D3]`, and `K_W=D2+2D3`.

Under the additional assumption of exactly one dicritical component, global
branch-curve topology gives the same two groups and narrows the local types:

- only `6T15=A6` with `(e,d)=(3,1)` or `6T16=S6` with `(2,1)` remains;
- the branch normalization `A1 -> B` must be noninjective, so the nonproperness
  curve has a finite multibranch singularity; and
- in the `A6` passport, every normalization collision consists of two smooth
  multiplicity-three branches, uses all six sheets, and is omitted by the
  original affine map; moreover, the two jump units concentrate at a unique
  local-degree-five point, and its branch is a `(2,5)` cusp.  Orevkov's 2026
  local classification contains exactly this block: an explicit finite
  degree-five polynomial germ has Jacobian
  `(15/8)(x^2-y)^2`, maps its critical parabola onto the cusp, and satisfies
  the required componentwise-bijective pullback condition.  Thus local
  analytic classification confirms the survivor rather than eliminating it;
- in the `S6` passport, the complete nine-row fiber census is known, each jump
  block has local group `S_(2+kappa)`, and every jump is exactly one of
  `T(2,3)`, `T(3,4)`, or `T(4,5)`; and
- in **both** passports, a cyclic endpoint cover contradicts any nonempty
  contracted constant chain in Orevkov's minimal model.  Thus its `L_C` is
  empty and the affine finite
  normalization `W` is smooth, finite flat of rank six over `A2`, with
  `W-D=A2`, `D=A1`, `Pic(W)=Z[D]`, and `K_W=(e-1)D`; and
- on Orevkov's original source-blowup resolution, the dicritical is then a
  type-3 leaf of augmented-canonical label `e`.  If its self-intersection is
  `-m`, its unique neighbor is type 2 with label `e*m-1`, and the path toward
  the original line at infinity necessarily crosses consecutive labels
  `1--0--(-1)`.  Jointly minimizing the source and extended map forces
  `m=1`, so the minimal `S6` and `A6` leaf edges are exactly `2--1` and
  `3--2`.  Their determinant labels satisfy `d_E=d_A-1<0` with forced parity.
  An exact hostile blowup family shows that these canonical-plus-determinant
  constraints alone still admit infinitely many abstract trees.

One significant compactification stratum is now completely excluded.  If the
irreducible dicritical image has exactly two characteristic pairs at infinity
and satisfies Orevkov's condition `(*)`—some target boundary component has
only one noncontracted irreducible preimage—then his splice product identity
forces degree at least eight in the `A6` case.  Equality is arithmetically
possible at the first `S6` star, but the final edge equation becomes
`3*Q_tilde=14`.  Thus neither passport survives in that conditional stratum.
The extra infinity hypotheses are not currently consequences of the earlier
one-dicritical reductions, so this is not a full elimination.  See the
[two-pair infinity note](docs/two-pair-infinity-elimination.md).

The natural attempt to remove condition `(*)` is now exact and known to be
insufficient on its own.  With several noncontracted preimages, Orevkov's
scalar determinant ratio becomes `Q*m=q*n`; Hodge index adds the sharp matrix
inequality `Q-(q/6)*n*n^t <= 0` when `q>0`.  Explicit unimodular trees produced
by boundary blowups satisfy all of these identities with split `A6` `3+3`
and `S6` `2+2+2` inertia.  Thus determinant transport, Hodge index, degree,
and parity do not force a unique preimage.  A stronger `A6` fixture also
retains valency two, the local augmented-canonical pullback labels `(-5,-1)`,
monomial normal forms, and primitive `A6` generation.  What remains must
couple all boundary components through global pullback, effectivity, or splice
data.  See the
[split-boundary note](docs/split-boundary-transport.md).
Orevkov's relevant first-star determinant has the opposite raw-intersection
sign, so the Hodge inequality is unavailable there rather than stronger.

There is now a global algebraic stopping model as well.  The finite flat map
`(x,z) -> (x^2+z, x^3-x^4*z)` has degree six, extends to a finite morphism
`P1 x P1 -> P2`, and has geometric monodromy `S6`.  The sole target boundary
line pulls back as `2*C1+C2`, where both components are noncontracted of
degrees one and four; condition `(*)` therefore genuinely fails for this
cover.  Two source blowups even reproduce the numerical `2--1--0--(-1)`
corridor and determinant parity.  Its Jacobian is nonconstant and it has
affine branching, so it is not a Keller map.  This isolates affine
unramifiedness and the genuine dicritical type assignment as indispensable
remaining inputs.  See the
[finite split-cover note](docs/s6-split-boundary-cover.md).

The eliminated one-dicritical types `(2,2)` and `(4,1)` would have injective
normalization.  Lin--Zaidenberg then makes the branch a monomial contractible
curve whose weighted-orbit product identifies the homotopy types of its local
and global complements; its intransitive local six-sheet action cannot equal
the transitive global monodromy.

These are necessary conditions, not constructions or an exclusion of `A6`
and `S6`.  The refined identity and universal elimination are proved in the
[refined budget note](docs/refined-six-sheet-budget.md), with the broader
geometric setup in the
[six-sheet monodromy note](docs/six-sheet-monodromy.md).  The universal
two-curve `S6` collision theorem is in the
[two-curve collision note](docs/s6-two-curve-collisions.md), which also proves
the saturated smooth-normalization package.  The
contracted-source obstruction and finite-flat consequences are in the
  [smooth-normalization note](docs/one-dicritical-source-smoothness.md).  The
  [A6 local note](docs/a6-one-dicritical-local.md) now includes the exact
  Orevkov germ and its complete three-component cusp pullback.  The
  [canonical-label note](docs/one-dicritical-leaf-labels.md) derives the forced
leaf-to-infinity corridor and states its compactification-model boundary.  The
dependency-free
[Python certificate](scripts/six_sheet_monodromy.py) rebuilds the exact groups,
classes, normal closures, normalizers, deck groups, blocks, and local subgroup
orbits; an optional [Sage/GAP checker](tools/check_six_sheet_gap.sage) verifies
the catalogue independently.

The smooth surface package is not itself contradictory.  Explicit
Hirzebruch-surface complements realize the one-boundary patterns
`K_W=D` and `K_W=2D`, and an explicit triple blowup realizes the disjoint
two-boundary pattern `K_W=D2+2D3`, including the required Picard groups and
`A2` interiors.  These are hostile consistency models, not finite covers.
They prove that the next obstruction must use the multiplication, trace,
discriminant, or monodromy of the finite rank-six algebra.  See the
[hostile-model note](docs/smooth-normalization-hostile-models.md).

The finite algebra does add a new exact obstruction.  Its relative dualizing
module is the nontrivial inverse-different line
`O_W(D)`, `O_W(2D)`, or `O_W(D2+2D3)`.  Hence the rank-six algebra is locally
Gorenstein but not globally Frobenius, not monogenic, and not one global
square complete intersection over `C[P,Q]`.  Splitting the ramification
divisor as `2E+O` factors the trace Gram matrix exactly as

```text
T = Phi^* H Phi,
```

with `H` symmetric.  The trace determinants are respectively `b`, `b^2`, and
`b2*b3^2`.  The cokernels are the
branch-normalization modules, so the exact
matrix coranks at every cusp and collision are fixed.  Those coranks remain
consistent; this is a concrete matrix target, not yet a contradiction.  See
the [trace-lattice note](docs/finite-flat-trace-lattices.md).

That consistency is now witnessed by exact hostile matrices.  A symmetric
`A6` model has `det(Phi)=b`, unimodular middle form, `det(T)=-b^2`, the required
generic/cusp/collision coranks, and a primitive norm-six vector.  A symmetric
`S6` model realizes the entire three-`T(2,3)` jump partition, three allowed
normalization collisions, and the same unit-line condition.  Its projective
homogenization carries the expected theta characteristic.  These models have
no commutative rank-six multiplication law; they prove that the next
obstruction must use the cubic trace tensor, integrality, and associativity,
not just the quadratic trace form.  See the
[hostile-matrix note](docs/trace-hostile-matrices.md).

Global topology now kills that particular `S6` matrix curve.  Projection to
one affine coordinate has four strands, so its complement is generated by
four meridians.  Four transpositions cannot act transitively on six sheets:
their four-edge graph on the sheets cannot be connected.  An exact
Zariski--van Kamp replay checks all `15^4` transposition assignments; `735`
satisfy the relators and none is transitive.  This is a real exclusion of the
explicit trace curve, but not of the `S6` passport: any replacement branch
curve must have projection width at least five.  See the
[trace-curve topology note](docs/s6-trace-curve-topology.md).

The first curve at equality is now excluded too.  The parametrization
`(t^5+t^4,t^7+t^5)` has exact width five, one allowed `T(4,5)` jump cusp, six
transverse normalization nodes, one `(2,7)` pair at infinity, and complete
genus accounting.  Nevertheless, exact Zariski--van Kamp simplification gives
complement group `Z`: all five geometric meridians coincide.  A separate
`15^5` transposition census finds only 15 cyclic images and no transitive one.
This kills a near-miss in the degree-minimal singular-one-pair slot, not all
width-five `S6` curves.  See the
[width-five near-miss note](docs/s6-width-five-near-miss.md).

The full group at infinity now removes much more of that stratum.  Assume a
polynomial normalization has degrees `m<n` and exactly one genuine singular
pair at infinity.  Its affine link is `T(m,n)`, whose group surjects onto the
global complement.  Requiring a transposition meridian and an exact `S6`
quotient eliminates every width-compatible pair with `n<=11`; an exhaustive
permutation census verifies the finite exceptional cases.  The bound is
sharp: `(m,n)=(5,12)` has `720` valid `S6` images, including an explicit
5-cycle/6-cycle witness.  Thus `(5,12)` is the first one-pair
infinity-topology target, not a constructed cover or Keller map.  Multi-pair
infinity remains open.  See the
[S6 one-pair infinity note](docs/s6-one-pair-infinity.md).

The explicit `A6` trace curve is globally dead as well: its affine complement
group is `Z`, so single-3-cycle meridians cannot generate `A6`.  More broadly,
the conditional one-pair audit now works under four standing hypotheses: the
branch normalization is `A1` and is represented polynomially; its projective
closure has exactly one genuine Puiseux pair at infinity; its only intrinsic
finite singularity is the forced `T(2,5)` cusp, with no additional
normalization preimage over the cusp image; and every other finite
singularity is a collision of smooth normalization points.  Genus and
link-at-infinity arithmetic reduce collision delta `Delta<=2` to three degree
pairs.  Two fail exact torus-group censuses.  The third has the exhaustive
family `(t^2+t^3,c*t^4+t^5)`; proper equisingular isotopy reduces its generic
part to the cyclic curve, while the only valid exceptional fiber has only
`C3` and `A5` three-cycle images.  Continuing the exact genus and large-link
census, together with the family audits below, excludes every `Delta<=7`.
The large-link census then excludes `Delta=8,9`.  Consequently

\[
  \boxed{\Delta\ge10}
\]

under the four hypotheses.  At the coarse link stage equality leaves only
affine degrees `(a,d)=(4,9)`, with projective infinity pair `(5,9)`.

That equality family is now explicit.  Every normalization is polynomially
equivalent to

\[
  P=t^2+kt^3+t^4,
  \qquad
  Q=at^5+bt^6+ct^7+dt^8+t^9,
  \qquad a\ne0,
\]

up to one residual involution.  The exact member
`(t^2+t^3+t^4,t^5+t^9)` has one `T(2,5)` cusp, ten nodes, a `T(5,9)` branch
at infinity, and cyclic affine complement.  Its stored four-generator van
Kampen presentation has exactly 40 single-three-cycle assignments, all with
image `C3`, and no `A6` image.  Proper projective Whitney--Thom transport
therefore excludes the entire connected clean stratum.  A hostile independent
census also shows why local link data alone cannot finish the boundary: all
720 admissible `T(4,9) -> A6` pairs have compatible cusp, collision,
orientation, and `2.A6` spin data.

The two dominant degeneration divisors are now excluded as well.  Exact
localized Gröbner calculations prove that the contact-two and ordinary-triple
incidences each have one irreducible four-dimensional dominant component.
Exact representatives on both components have cyclic complement; both raw
`40^4` replays again have only 40 diagonal `C3` images.  Whitney--Thom
transport excludes their connected generic equisingular opens.  Consequently
any remaining conditional delta-ten survivor lies in a lower-dimensional
degeneration stratum.

The six expected codimension-two profiles now have exact component-level
certificates as well.  The valid generic components of contact-three and
ordinary-quadruple incidence are irreducible; the displayed `P`-unramified
`T112` incidence is an irreducible affine-line bundle, and the displayed
Cramer calculations identify rational threefolds for two contacts, a separate
contact and ordinary triple, and two ordinary triples.  An exact
representative of each has cyclic complement, and every raw `40^4` replay
again leaves only 40 diagonal `C3` images and no `A6` image.  All six
displayed generic or dominant components are excluded at that level.  For
`T112` and the mixed contact-plus-triple chart, the proof works over the
smooth irreducible labeled incidence space: a finite-etale cover labels the
clean collision sections, relative blowups give a simultaneous embedded
resolution, and proper Whitney--Thom isotopy transports the cyclic sample.
This avoids any assumption that the coefficient-image threefold is smooth or
normal and does not require the `T112` image map to have degree exactly two.

The three residual coefficient-rank factors are now closed at the threefold
level: the mixed factor has rank exactly three on its valid divisor, the
two-triple residual is empty after localization, and all two-contact residual
incidences have dimension at most two.  The immersed `P`-critical `T112` and
mixed boundaries likewise have dimension at most two.  On the true split
charts, a generated 22-row allocation ledger has exact clean witnesses at
`k=0,+2` and full replay under transport to `k=-2`.  That ledger is not a
global split exclusion: most coefficient-rank-drop subloci, denominator and
overlap charts, component intersections, and the lower-dimensional residual
taxonomy remain open.

The first of the fourteen expected codimension-three profiles is now attacked
as well.  For `C4+6N`, the determinant-nonzero nonsplit incidence is one
rational surface.  Exact Sage saturation shows that its residual determinant
curve has a finite length-ten compatibility scheme with rank-three
affine-line fibers, so it hides no second surface.  An exact member has one
contact-four point, six nodes, cyclic complement, and no `A6` assignment.
Four relative contact blowups and proper Whitney--Thom transport exclude the
dense clean Cramer surface.  The split, denominator, cusp-pair, diagonal, and
residual boundary pieces remain open.  The complete combinatorial ledger
still contains 145 candidate collision
profiles, including 55 overdetermined profiles deliberately retained until
exact saturation proves them empty, invalid, or contained elsewhere.  No
`A6` cover or Keller map has been constructed, the four hypotheses remain
unproved for arbitrary Keller branches, and this is not a proof of `JC(2)`.
See the [A6 one-pair note](docs/a6-one-pair-infinity.md), the
[generic delta-ten audit](docs/a6-delta-ten-generic.md), the
[dominant delta-ten wall audit](docs/a6-delta-ten-walls.md), and the
[codimension-two checkpoint](docs/a6-delta-ten-codim-two.md), and the
[codimension-three checkpoint](docs/a6-delta-ten-codim-three.md).

The delta-five equality family is now fully exhausted.  Every conditional
`(3,8)` curve is polynomially equivalent to
`(t^2+t^3,alpha*t^5+beta*t^6+gamma*t^7+t^8)`.  Exact collision resultants split
its valid parameter space into a generic five-node stratum, an ordinary-triple
divisor, and a contact-two divisor.  Sage's presentation simplifier reports
complement group `Z` for clean representatives of all three; each `40^3`
single-3-cycle census has only 40 cyclic images.  Exact primary decomposition
splits the residual into four valid
rational curves.  Their four generic presentations and the four presentations
at every valid exceptional point again have exactly 40 cyclic images and no
`A6` image.  Exact algebra and finite permutation replay are separated from
the computer-assisted Zariski--van Kamp extraction and Whitney--Thom transport.
Consequently `Delta=5` is conditionally impossible; this does not eliminate
the unrestricted `A6` passport.  See the
[A6 delta-five family note](docs/a6-delta-five-family.md).

At the historical `Delta=7` coarse-link stage, every conditional survivor had
affine degrees `(3,10)` and the four-parameter normal form
`(t^2+t^3, alpha*t^5+beta*t^7+gamma*t^8+delta*t^9+t^10)`.  The exact member
`(t^2+t^3,2*t^5+t^10)` has the forced cusp, seven transverse nodes, a `(7,10)`
infinity pair, and complete genus accounting, but its affine complement is
`Z`; the exhaustive `40^3` replay has no `A6` image.  Proper Whitney--Thom
transport excludes the connected nondegenerate four-parameter open.  The
subsequent wall audit now excludes the full triple-image wall, every
positive-dimensional repeated-collision stratum, and every finite endpoint.

The normal form, coefficient-slice algebra, partition ledger, resultants,
saturations, Groebner bases, and finite permutation censuses are exact and
replayable.  The topology layer separately depends on stored
Zariski--van Kamp presentations, proper Whitney--Thom propagation over
connected equisingular strata, and finite-etale/Riemann-existence transport
with tame inertia across arithmetic endpoint embeddings.  The stored
presentations are replayed exactly, but not every original presentation
extraction is regenerated by a checked script.  Subject to those dependencies
and the four standing hypotheses, `Delta=7` is therefore fully excluded.  The
multi-pair case and unrestricted `A6` passport remain open.  See the
[generic delta-seven note](docs/a6-delta-seven-generic.md) and the
[complete wall audit](docs/a6-delta-seven-walls.md).

The first cubic-cover obstruction also survives an exact lift test.  In the
natural spin double cover `2.A6`, the canonical order-three lifts of the
forced `T(2,5)` cusp meridians satisfy the five-braid relation exactly, and
the two disjoint degree-three collision meridians still commute.  These lifts
generate all of `2.A6`.  More decisively, the forced prefix has transitive
product-one completions made entirely of 3-cycles with both Fried--Serre spin
signs.  Thus finite local monodromy does not determine a spin obstruction;
one must derive an infinity word and framing from the Keller compactification.
See the [spin-lift note](docs/a6-spin-lift.md).

The multiplication enhancement is now an exact finite target as well.  A
completely symmetric cubic tensor with `56` polynomial entries must satisfy
explicit divisibility, unit, WDVV associativity, ordinary-regular-trace, and
middle-lattice equations.  Normalization turns every apparent division by
the branch equation into one polynomial identity in `t`.  Both forced special
fibers nevertheless survive: the exact trace data admit
`C[z]/(z^5) x C` at the cusp and
`C[z]/(z^3) x C[w]/(w^3)` at a collision, including compatible perfect
middle forms and divisor sections.  Thus any multiplication obstruction is
global, not pointwise.  See the
[multiplication-tensor note](docs/a6-multiplication-tensor.md).

## A separate sparse obstruction at coordinate degree `(72,108)`

Generic sheet degree and coordinate degree are different invariants.  On the
coordinate-degree side, Guccione--Guccione--Horruitiner--Valqui reduce the
remaining sub-`125` problem to `(72,108)` and its transpose, with two explicit
transformed Newton-polygon configurations satisfying `[P,Q]=x^2`.

An exhaustive exact support calculation now proves that the first
configuration needs at least three nonzero coefficients strictly inside its
two Newton polygons, while the second needs at least four.  All boundary
lattice coefficients remain arbitrary; only the exact polygon vertices are
assumed nonzero.  The checker certifies all `7504` first-case supports with at
most two interior terms and all `3683` second-case supports with at most three.
Of the latter, `3678` have replayed zero-product certificates and the five
remaining triples have exact unit-ideal certificates.  Hostile fixtures
confirm that the method stops on the full polygons and on a named four-term
second-case support.  This is a sparse-support lower bound, not an elimination
of `(72,108)`.  See the [Newton-polygon note](docs/newton-72-108-sparse.md) and
its [exact certificate](scripts/newton_72_108.py).

## The screenshot is exactly correct

For

\[
\begin{aligned}
A&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
B&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
C&=2x-3x^2y-x^3z,
\end{aligned}
\]

the map `F=(A,B,C)` satisfies the polynomial identity

\[
  \det JF=-2.
\]

The three distinct rational points

\[
\left(0,0,-\tfrac14\right),\quad
\left(1,-\tfrac32,\tfrac{13}{2}\right),\quad
\left(-1,\tfrac32,\tfrac{13}{2}\right)
\]

all map exactly to `(-1/4,0,0)`.  This is a complete finite certificate of a
counterexample in dimension three.  Appending identity coordinates gives the
same conclusion in every higher dimension.

[`JacobianTwo/Counterexample.lean`](JacobianTwo/Counterexample.lean) constructs
the formal Jacobian from actual `MvPolynomial.pderiv` entries, proves its
determinant, proves the three-point collision and pairwise distinctness, and
derives noninjectivity over `ℂ`.  The independent typed SymPy checker
[`scripts/verify.py`](scripts/verify.py) recomputes both identities using exact
rational arithmetic and includes hostile transcription fixtures.

## A solved follow-on problem: every fiber and every escaping value

For a target `(a,b,c)`, introduce the reciprocal fiber coordinate

\[
  T=y+\frac1x
\]

on `x != 0`.  It satisfies

\[
  p(T)=cT^3-2T^2+bT-2a=0,
  \qquad
  p'(T)=\frac2x.
\]

Define

\[
  Q(a,b,c)=27a^2c^2-18abc+16a+b^3c-b^2
\]

and

\[
  \Gamma=\{3bc=4,\ b^2=12a\}.
\]

The exact fiber calculation gives

| Target stratum | Number of source points |
|---|---:|
| `Q != 0` | 3 |
| `Q = 0` and target not in `Gamma` | 1 |
| target in `Gamma` | 0 |

Consequently,

\[
  F(\mathbb C^3)=\mathbb C^3\setminus\Gamma.
\]

Moreover, the complete nonproper-value set—the targets approached by images
of source sequences escaping to infinity—is

\[
  S_F=V(Q).
\]

The proof includes both directions.  Every point of `V(Q)` has an explicit
escaping family, while projective-root compactness plus exact reconstruction
shows that no point outside `V(Q)` can be an asymptotic value.  See
[`docs/nonproper-set.md`](docs/nonproper-set.md) for the complete argument.

[`JacobianTwo/CubicFiber.lean`](JacobianTwo/CubicFiber.lean) certifies the
fiber cubic, its derivative, the standard universal cubic discriminant
coefficient expression `-4Q`, an explicit Bézout common-root certificate,
finite-root reconstruction, and all large-`T` cancellation identities.
[`scripts/nonproper.py`](scripts/nonproper.py)
independently checks the remaining exact algebra: the infinity chart,
repeated- and triple-root parameterizations, singular-locus elimination, and
the escaping family.  The compactness argument is written explicitly in the
mathematical note rather than being mislabeled as kernel-checked topology.

These consequences are labeled **derived here; historical priority not
established**.  Same-day sources already contained the cubic, reconstruction,
discriminant, and generic `S_3` calculation; this repository makes no
literature-priority claim for the fuller stratification.

## Certified positive fragments in the plane

The strongest plane result in this repository has no degree bound on its
first coordinate.  Let

\[
  F=(P,Q),\qquad P\in K[x,y],\qquad Q=e(x)y+f(x),
\]

over a characteristic-zero field.  If the actual formal Jacobian is a
nonzero scalar, Lean proves that `F` is a polynomial automorphism.  In the
`e != 0` chart it derives

\[
  e=\varepsilon\in K^\times,\qquad
  P=G(Q)+\alpha x+\beta,\qquad
  \alpha\varepsilon=k,
\]

and in the `e=0` chart it derives the complementary triangular form.  Both
charts have kernel-checked explicit inverses.  See
[`JacobianTwo/AffineCoordinate.lean`](JacobianTwo/AffineCoordinate.lean) and
the [proof and literature note](docs/affine-coordinate.md).

This is a characteristic-zero algebraic formalization of the known
type-`(m,1)` reduction, not a new mathematical class.  Sabatini's published
real theorem uses the same leading-power elimination.  The repository's
field-uniform statement deliberately assumes a genuinely constant Jacobian.

The next certified class allows a genuinely quadratic coordinate.  Let

\[
  Q=\varepsilon y^2+g(x)y+f(x),\qquad \varepsilon\in K^\times,
  \qquad s=2\varepsilon y+g(x).
\]

For arbitrary `P`, a nonzero constant identity `J(P,Q)=k` forces the
discriminant `Delta=g^2-4*epsilon*f` to be affine, say `Delta=A*x+B` with
`A != 0`, and Lean proves the original-coordinate normal form

\[
  P=G(Q)+\lambda s,\qquad \lambda A/2=k.
\]

It also proves both laws for the explicit polynomial inverse

\[
  \sigma=(u-G(v))/\lambda,\quad
  x=(\sigma^2-4\varepsilon v-B)/A,\quad
  y=(\sigma-g(x))/(2\varepsilon).
\]

See
[`JacobianTwo/ConstantLeadingQuadratic.lean`](JacobianTwo/ConstantLeadingQuadratic.lean)
and the [proof, certificate map, and literature boundary](docs/constant-leading-quadratic.md).
This theorem has no degree bound on `P`; it does require the quadratic leading
coefficient of `Q` to be a nonzero scalar.

The scalar-leading hypothesis is not needed in the known mathematical theorem.
If

\[
  Q=a(x)y^2+g(x)y+f(x),\qquad a\ne0,
\]

and `J(P,Q)` is a nonzero scalar, Moskowicz's Theorem 2.7 already implies that
`(P,Q)` is an automorphism: its invariant
`gcd(2, deg_x(a))` is either `1` or the prime `2`.  Simon--Weimann's coordinate
criterion then implies, after scalar extension if necessary, that `a` is a
nonzero scalar and that `deg_x(g^2-4*a*f)=1`.

[`JacobianTwo/VariableLeadingQuadratic.lean`](JacobianTwo/VariableLeadingQuadratic.lean)
develops an independent direct certificate for that known theorem.  Lean now
certifies the top-coefficient equation, target-shear descent to an odd
`y`-degree, the identity `p_n^2=c*a^n`, the UFD shape
`a=epsilon*h^2, p_n=lambda*h^n`, and the unique even--odd decomposition over a
field.  The fraction-field layer also certifies the quotient-rule derivation
on `K(x)`, its constant field, the affine substitution `y=(U-rho)/h`, the
exact Jacobian factor `k/h`, parity extraction, and the full coefficient
recurrence.  Lean also constructs its primitive explicitly, proves every
coefficient lies in `K[F]`, and tracks the exact degree and nonzero leading
coefficient through the downward descent.  Specialization at `y=0`, a
valuation-free gcd normalization, and a unique-survivor denominator theorem
then force `h | g`; the terminal recurrence forces `h` to be a unit.  The
resulting theorem `variableLeadingQuadratic_bijective_full` is a complete
kernel-checked proof for every coordinate of the displayed at-most-quadratic
form, including its affine branch.  The complete direct proof, hostile
fixtures, and exact literature boundary are in
[`docs/variable-leading-quadratic.md`](docs/variable-leading-quadratic.md).

Two earlier bounded modules expose useful intermediate mechanisms.  First,

\[
  (x,y)\longmapsto(A(x)y+B(x),\ C(x)y+D(x))
\]

cannot be a noninjective Keller map.  A nonzero constant determinant forces
the variable slopes `A` and `C` to be constant, after which a linear
combination of the outputs recovers `x` and then `y`.
[`JacobianTwo/AffineInOneVariable.lean`](JacobianTwo/AffineInOneVariable.lean)
contains the proof.

Second,

\[
  (x,y)\longmapsto
  (a(x)y^2+b(x)y+c(x),\ e(x)y+f(x))
\]

is reduced by its constant-Jacobian coefficient equations to an explicit
triangular normal form in the nonzero-`e` chart, with a displayed polynomial
inverse.  The complementary `e=0` chart is handled separately in the final
theorem.  See
[`JacobianTwo/QuadraticInOneVariable.lean`](JacobianTwo/QuadraticInOneVariable.lean).
These bounded results are now subsumed by the arbitrary-degree theorem, but
their shorter coefficient proofs remain useful.  None of these statements is
a proof of general `JC(2)` or of the generic-degree-six frontier.

## Reproduce the certificates

The Lean toolchain and mathlib revision are pinned.  Python dependencies are
locked by `uv.lock`.

```bash
lake build
uv run --frozen python -m scripts.verify
uv run --frozen python -m scripts.nonproper
uv run --frozen python -m scripts.affine_coordinate
uv run --frozen python -m scripts.constant_leading_quadratic
uv run --frozen python -m scripts.variable_leading_quadratic --depth 9
uv run --frozen python -m scripts.six_sheet_monodromy
uv run --frozen python -m scripts.canonical_leaf_graph
uv run --frozen python -m scripts.newton_72_108
uv run --frozen pytest
uv run --frozen mypy
# Optional independent finite-group cross-check:
sage tools/check_six_sheet_gap.sage
```

The Lean source contains no `sorry`, `admit`, or custom axiom.  CI runs the
Lean build, every displayed `uv` command, and the unfinished-proof check.  The
optional Sage/GAP replay is an additional independent local cross-check.

## Reading map

- [`SPEC.md`](SPEC.md) is the research specification and claim-status ledger.
- [`docs/galois-frontier.md`](docs/galois-frontier.md) explains the Galois
  misconception and identifies generic degree six as the first open frontier.
- [`docs/six-sheet-monodromy.md`](docs/six-sheet-monodromy.md) proves the
  degree-six monodromy filters and the conditional one-dicritical passports.
- [`docs/refined-six-sheet-budget.md`](docs/refined-six-sheet-budget.md)
  proves `sum(e*d+delta)=5`, eliminates `A5` and `S5` universally, and lists
  the exact surviving `A6`/`S6` ramified profiles.
- [`docs/a6-one-dicritical-local.md`](docs/a6-one-dicritical-local.md)
  rules out the `A6` jump partition `1+1`, proves the unique
  local-degree-five point, and classifies its smooth-source branch as a
  `(2,5)` cusp.
- [`docs/a6-one-pair-infinity.md`](docs/a6-one-pair-infinity.md) derives the
  conditional genus/link framework and records the combined one-pair
  frontier.
- [`docs/a6-delta-five-family.md`](docs/a6-delta-five-family.md) exhausts the
  earlier conditional `(3,8)` equality family, including its residual curves
  and exceptional points.
- [`docs/a6-delta-seven-generic.md`](docs/a6-delta-seven-generic.md) proves
  cyclicity on the connected nondegenerate open of the conditional `(3,10)`
  family.
- [`docs/a6-delta-seven-walls.md`](docs/a6-delta-seven-walls.md) exhausts the
  conditional `(3,10)` collision walls and finite endpoints, derives
  `Delta>=10`, and separates exact algebra and finite-group replay from the
  stored-presentation, Whitney--Thom, and finite-etale dependencies.
- [`docs/a6-delta-ten-generic.md`](docs/a6-delta-ten-generic.md) derives the
  complete normalized `(4,9)` family, treats the exceptional pair-incidence
  charts honestly, and excludes its connected clean stratum by an exact
  cyclic-complement representative and exhaustive `40^4` replay.
- [`docs/a6-delta-ten-walls.md`](docs/a6-delta-ten-walls.md) proves that the
  contact-two and ordinary-triple incidences are the two irreducible dominant
  wall components, excludes their generic equisingular opens, and gives the
  exact 145-profile ledger for the deeper audit.
- [`docs/a6-delta-ten-codim-two.md`](docs/a6-delta-ten-codim-two.md) excludes
  all six displayed generic or dominant components, while retaining the
  residual-rank, split-chart, removed
  `P`-projection/critical-fiber, and deeper-intersection loci as explicit open
  obligations.
- [`docs/a6-delta-ten-propagation.md`](docs/a6-delta-ten-propagation.md)
  verifies the connected-clean-open, finite-etale labeling, simultaneous
  embedded-resolution, and proper-isotopy steps for the `T112` and mixed
  contact-plus-triple dominant charts.
- [`docs/a6-delta-ten-codim-three.md`](docs/a6-delta-ten-codim-three.md)
  excludes the dense clean nonsplit `C4+6N` Cramer surface, including an exact
  residual-rank saturation, contact-four singular scheme, cyclic complement,
  and four-blowup topology transport, while retaining every omitted boundary.
- [`docs/one-dicritical-source-smoothness.md`](docs/one-dicritical-source-smoothness.md)
  eliminates the complete contracted constant chain in both surviving
  one-dicritical passports and derives the smooth finite-flat normalization.
- [`docs/one-dicritical-leaf-labels.md`](docs/one-dicritical-leaf-labels.md)
  fixes the adjacent canonical label `e*(-E^2)-1`, jointly minimizes the leaf
  to `E^2=-1`, proves the neighbor is type 2, and forces a `1--0--(-1)`
  transition toward the negative core.
- [`docs/smooth-normalization-hostile-models.md`](docs/smooth-normalization-hostile-models.md)
  constructs explicit affine surface pairs realizing every smooth/Picard/
  canonical package above, thereby isolating the missing finite-cover data.
- [`docs/finite-flat-trace-lattices.md`](docs/finite-flat-trace-lattices.md)
  derives the inverse different, exact discriminants, the symmetric
  half-different factorization, normalization-module cokernels, and the
  nonmonogenic obstruction.
- [`docs/s6-one-dicritical-local.md`](docs/s6-one-dicritical-local.md)
  gives the complete `S6` fiber census, symmetric local blocks, jump
  trichotomy, torus-knot types, and hostile analytic models.
- [`docs/s6-two-curve-collisions.md`](docs/s6-two-curve-collisions.md)
  proves the unrestricted collision rows and the saturated two-boundary
  smooth finite-flat normalization.
- [`docs/newton-72-108-sparse.md`](docs/newton-72-108-sparse.md) gives the
  exact sparse-support obstruction in the separate residual coordinate-degree
  configurations.
- [`docs/nonproper-set.md`](docs/nonproper-set.md) proves the complete fiber,
  image, and nonproper-set theorem.
- [`docs/affine-coordinate.md`](docs/affine-coordinate.md) proves the
  arbitrary-degree affine-coordinate normal form and marks its exact
  literature boundary.
- [`docs/constant-leading-quadratic.md`](docs/constant-leading-quadratic.md)
  proves the constant-leading quadratic-coordinate normal form and displays
  its polynomial inverse.
- [`docs/variable-leading-quadratic.md`](docs/variable-leading-quadratic.md)
  gives the known arbitrary-leading quadratic theorem, an independently
  derived direct proof, and its complete Lean certificate.
- [`docs/audit.md`](docs/audit.md) gives a hand-checkable structural derivation
  of the original screenshot.
- [`docs/research-log.md`](docs/research-log.md) records completed work,
  negative results, and remaining obligations.

## Sources and provenance

- Levent Alpöge's [original public announcement on X][announcement]
- L. Andrew Campbell's [Galois-case theorem][campbell]
- S. Yu. Orevkov's [three-sheet theorem][orevkov]
- A. V. Domrina's [four-sheet theorem][domrina]
- Henryk Żołądek's [result through generic degree five][zoladek]
- Alexander Borisov's [Keller-map compactification framework][borisov]
- Guccione--Guccione--Horruitiner--Valqui's [(72,108) reduction][guccione]
- Vered Moskowicz's [quadratic-coordinate antecedent][moskowicz]
- Denis Simon and Martin Weimann's [coordinate/discriminant criterion][simon-weimann]
- Marco Sabatini's [type-`(m,1)` triangular reduction][sabatini]
- Zihan Zhang's [direct-consequences note][consequences]

The announcement and expository note establish provenance.  The finite
claims made here are supported by the repository's reproducible Lean and exact
symbolic certificates.

[announcement]: https://x.com/__alpoge__/status/2079028340955197566
[campbell]: https://doi.org/10.1007/BF01349234
[orevkov]: https://doi.org/10.1070/IM1987v029n03ABEH000984
[domrina]: https://doi.org/10.1070/im2000v064n01ABEH000273
[zoladek]: https://doi.org/10.1016/j.top.2008.04.001
[borisov]: https://www.combinatorics.org/ojs/index.php/eljc/article/view/v27i3p54
[guccione]: https://arxiv.org/abs/2204.14178
[moskowicz]: https://arxiv.org/abs/1810.08202
[simon-weimann]: https://doi.org/10.1216/JCA-2018-10-4-559
[sabatini]: https://doi.org/10.4064/cm9195-1-2024
[consequences]: https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html

## License

Apache-2.0. See [`LICENSE`](LICENSE).
