# Conditional `A6` delta-ten codimension-two checkpoint

## Claim boundary

This note continues the [generic delta-ten audit](a6-delta-ten-generic.md)
and the [codimension-one wall audit](a6-delta-ten-walls.md).  Every statement
below is conditional on the same four additional hypotheses for the
one-dicritical `A6` branch:

1. its normalization is `A1` and is represented by a polynomial map;
2. its projective closure has exactly one genuine Puiseux pair at infinity;
3. its only intrinsic finite singularity is the forced `T(2,5)` cusp, and no
   other normalization point maps to the cusp image; and
4. every other finite singularity is a collision of smooth normalization
   points.

Under those hypotheses, collision delta ten forces the normalized family

\[
  P=t^2+kt^3+t^4,
  \qquad
  Q=at^5+bt^6+ct^7+dt^8+t^9,
  \qquad a\ne0. \tag{0.1}
\]

The clean locus and both codimension-one divisors are already excluded.  This
checkpoint treats the six profiles of expected codimension two.  Four exact
component exclusions are now in hand.  The `T112` and mixed
contact-plus-triple calculations identify the expected threefolds and cyclic
members, but their generic exclusions retain explicit geometric propagation
obligations described below.  None of the calculations classifies every
component supported on a denominator, coefficient-rank, split `k=0,+2,-2`,
overlap, removed `P`-projection/critical-fiber chart, or deeper-degeneration
boundary.  A removed `P`-critical fiber is not automatically a singular
branch of `(P,Q)` when `Q'` is nonzero.  Those residual loci remain in the
audit.

The topology conclusion also has two explicit computer-assisted inputs:

- manual Sage 10.8 checkers regenerate an affine Zariski--van Kamp
  presentation for one exact rational point of each identified component and
  simplify its complement to `Z`; and
- component-wide conclusions separately use proper projective Whitney--Thom
  triviality on a verified connected clean equisingular open.

The dependency-free replay then checks all `40^4=2,560,000` assignments of
the four raw meridians to single three-cycles.  For every representative,
exactly forty assignments survive, all with diagonal image `C3`; none has
image `A6`.

This is a conditional branch-curve elimination, not a construction or
classification of degree-six Keller maps.  It does not derive the four
hypotheses, exclude the unrestricted `A6` or `S6` passports, or prove
`JC(2)`.

## 1. Status of the six profiles

Here `N` denotes an ordinary node, `C_m` a two-branch contact of intersection
multiplicity `m`, `T111` an ordinary triple point, `T112` a three-branch
point with contact multiplicities `(1,1,2)`, and `Q0` an ordinary quadruple
point.

| profile | component proved to exist | topology conclusion on that component | boundary not covered |
|---|---|---|---|
| `C3 + 7N` | one irreducible threefold on the valid pair chart | generic clean component excluded | removed projection-chart and split pair-chart loci, including `s=0` and pair-denominator behavior, plus deeper intersections |
| `C2^2 + 6N` | Cramer component dominating the full `(k,u,v)` ordered base, modulo the free swap | that dominant clean component excluded | residual-rank closures and split-fiber intersections may support other components |
| `T112 + 6N` | irreducible labeled affine-line bundle on the displayed `P`-unramified chart and generically finite threefold image | conditional generic-chart exclusion supported, not CAS-certified | geometric two-orientation labeling and connected clean Whitney--Thom propagation; removed `P`-critical-fiber and deeper loci, which are not automatically singular curve branches |
| `C2 + T111 + 5N` | dense rational Cramer component on the displayed open | cyclic sample proved; generic-chart exclusion still conditional | separately verify a connected equisingular clean open and proper Whitney--Thom propagation; residual-rank, split, overlap, removed triple/fourth-root `P`-critical, and deeper lower loci |
| `Q0 + 4N` | one irreducible threefold on the valid four-source fiber chart | generic clean component excluded | removed critical-fiber, non-clean, and deeper intersections |
| `T111^2 + 4N` | rational threefold dominating the dense valid Cramer open | dense Cramer component excluded | compatible residual-rank-supported and split boundaries |

“Excluded” means that the component's clean equisingular open cannot carry
the required connected six-sheet `A6` quotient with single-three-cycle branch
meridians.  It does not mean that the corresponding profile has been removed
scheme-theoretically from every boundary chart.  “Conditional generic
exclusion supported” means that the exact incidence and cyclic member are in
hand, while the stated labeling or connected-propagation argument remains a
separate mathematical dependency; neither CI nor Sage certifies that step.

## 2. Contact three: `C3 + 7N`

Let `H(s)` be the collision decic on the valid unordered-pair chart.  A
contact-three collision is a triple root, so its incidence equations are

\[
  H(s)=H'(s)=H''(s)=0. \tag{2.1}
\]

They are affine-linear in `(a,b,c,d)`.  The gcd of the maximal coefficient
minors is

\[
  2(k+2s)^3(s^2+ks+1)^3. \tag{2.2}
\]

After localizing away from the pair denominator, the cusp-image factor, the
diagonal factor, and `s=0`, the only possible rank-drop base is

\[
  k=-6s,\qquad 5s^2=3. \tag{2.3}
\]

An augmented minor reduces there to the nonzero rational number
`148635648/625`.  Thus those fibers are inconsistent rather than extra
components.  The valid incidence is an affine-line bundle over an
irreducible surface.  Monicity of `H` in `s` makes the ambient hypersurface
`H=0` finite over coefficient space; restricting to the localized valid open
is therefore quasi-finite, and its image closure is one irreducible
threefold.  This chart does not classify the removed `s=0` or
pair-denominator loci, including their split pair-incidence behavior.

The exact member

\[
 (k,a,b,c,d)=\left(1,3,0,\frac12,-\frac{21}{16}\right) \tag{2.4}
\]

has the forced cusp, one contact of order three, seven nodes, and the fixed
`T(5,9)` branch at infinity.  Its checked complement is cyclic, so the
generic clean component is excluded.

## 3. Two contacts: `C2^2 + 6N`

For two ordered prospective double roots `u,v` of `H`, impose

\[
  H(u)=H'(u)=H(v)=H'(v)=0. \tag{3.1}
\]

The four equations are affine-linear in `(a,b,c,d)`.  The visible part of
their determinant is

\[
 (k+2u)^2(k+2v)^2(u-v)^4
 (u^2+ku+1)(v^2+kv+1). \tag{3.2}
\]

On the nonvanishing open, Cramer's rule gives a rational ordered threefold
dominating the full `(k,u,v)` Cramer base.
The swap `u<->v` is free and the coefficient solution is invariant, so it
descends to unordered coordinates.  The gcd of the determinant with all four
augmented Cramer numerators contains only

\[
  (u-v)^4(u^2+ku+1)(v^2+kv+1), \tag{3.3}
\]

not the residual determinant factor.  Hence that residual factor does not
support another component dominating its divisor.  This gcd calculation does
not rule out compatible subvarieties of higher codimension on the residual
rank locus.

The exact member

\[
  (k,a,b,c,d)=\left(1,3,\frac32,3,0\right) \tag{3.4}
\]

has contact roots `u=-1`, `v=1`, distinct targets `(0,-1/2)` and
`(-2,7/2)`, and six nodes.  Sage verifies the singular scheme and a cyclic
complement.  At each split value `k=0,+2,-2`, the two genuine pair-component
discriminants are nonzero, squarefree, and coprime, of degree pairs `(2,13)`,
`(6,9)`, and `(6,9)`.  That rules out an automatic double-contact divisor on
a whole split fiber; it does not classify their codimension-two intersections.

Therefore that full-base dominant clean component is excluded.  Residual-rank
closures and split-fiber intersections remain open and could support other
components of the same collision profile.

## 4. One tangent triple: `T112 + 6N`

Write the three source parameters as `u,v,w` and put

\[
 \sigma_1=u+v+w,
 \quad \sigma_2=uv+uw+vw,
 \quad \sigma_3=uvw. \tag{4.1}
\]

On `\sigma_2\ne0`, they lie in one `P`-fiber precisely on the surface

\[
 E=\sigma_2^2-\sigma_1\sigma_3-\sigma_2=0, \tag{4.2}
\]

with fourth source `-\sigma_3/\sigma_2`.  The polynomial `E` is absolutely
irreducible.  Its singular ideal becomes the unit ideal after localizing away
from zero and repeated sources, `\sigma_2=0`, and coincidence with the fourth
source, so the displayed `P`-unramified root base is smooth.  Coincidence with
the fourth root is a `P`-projection/critical-fiber boundary.  It is not
automatically an invalid or singular branch of the parametrized curve when
`Q'` remains nonzero, and the full clean `T112` profile on that locus is not
covered by this chart calculation.

For a labeled tangent pair, say `(u,v)`, two target-equality equations and
one slope-equality equation form a `3 x 4` affine-linear system in
`(a,b,c,d)`.  Sage saturates the ideal of its four maximal minors by the exact
`P`-unramified localization product.  The raw rank-drop ideal has dimension
one, but the saturation is the unit ideal (with recorded exponent three),
proving rank three everywhere on that displayed chart.  The labeled incidence
is therefore an irreducible affine-line bundle of dimension three there.

The three tangent-pair choices form one orbit under relabeling.  On the exact
sample, the saturated projection fiber is reduced of length two, with basis

```text
v^2 + v + 6/25,
u + v + 1,
w - 1/5.
```

The reduced length-two computation proves generic finiteness.  Identifying
its degree as two additionally uses the geometric fact that a clean target in
this stratum has a unique `T112` fiber and that the only labels above it are
the two orientations of its unique tangent pair.  With that labeling input,
the clean unlabeled quotient and its threefold image are irreducible and
connected.

The rational member has

\[
 (k,a,b,c,d)=
 \left(2,\frac{2089}{11875},\frac{3542}{2375},
          \frac{1788}{475},\frac{338}{95}\right), \tag{4.3}
\]

with triple sources `-3/5,-2/5,1/5`; the first two have contact exactly two.
The remaining collision packets contain two and four nodes.  Sage finds
singular-scheme lengths `4+6+2+4` and a cyclic complement.  These facts
support exclusion of the generic clean `T112` component conditional on the
unique-fiber/two-orientation labeling argument and proper projective
Whitney--Thom propagation on the connected clean open.  The manual Sage
checker certifies the saturation, singular scheme, and presentation, not the
propagation theorem.  The removed `P`-critical-fiber loci and deeper
intersections are outside the claim; they are not declared geometrically
invalid.

## 5. A separate contact and triple: `C2 + T111 + 5N`

Parameterize an ordinary triple `P`-fiber by its second and third elementary
symmetric functions `(q,r)`, with

\[
 p=\frac{q(q-1)}r,
 \qquad
 k=\frac{r^2+q^2-q^3}{qr}. \tag{5.1}
\]

Two equations make `Q` constant on the triple cubic.  A separate prospective
contact with pair sum `w` supplies `H(w)=T(w)=0`.  The resulting four
equations are affine-linear in `(a,b,c,d)`.

After inverting the displayed triple, contact, split, overlap, discriminant,
fourth-root, and residual-rank factors, the `4 x 4` matrix is invertible.  The
triple-source/fourth-root factor is a `P`-critical-fiber boundary, not an
automatic singularity of `(P,Q)` if `Q'` is nonzero; that portion of the clean
profile is outside the displayed Cramer chart.
Cramer's rule gives one rational threefold over an irreducible open in
`(q,r,w)`.  The one nonboundary residual-rank factor is irreducible and has
gcd one with a selected augmented minor after the exact specialization
`(q,r)=(2,-1)`.  The factor is proved irreducible over `QQ`; no absolute
irreducibility claim is used.  It supplies no second component dominating
that rank divisor.  Compatible lower-dimensional intersections on the
divisor are not eliminated by this test.

The exact member, after harmless target scaling, is

\[
 X=2t^4+3t^3+2t^2,
 \qquad
 Y=t^9-2t^7+t^6-2t^5. \tag{5.2}
\]

It has the forced cusp, a contact-two fiber, a separate ordinary triple, and
five nodes.  The independent Sage checker verifies singular-scheme lengths
`4+3+4+5`, regenerates the raw presentation, and simplifies it to `Z`.  This
proves a cyclic member on the dense rational Cramer component.  A generic
clean exclusion is supported only conditional on separately verifying that
the relevant equisingular clean open is connected and then applying proper
projective Whitney--Thom propagation.  Neither the Python CI nor the manual
Sage checker proves that transport.  Residual-rank, split, overlap, and deeper
lower-dimensional loci remain open.

## 6. One ordinary quadruple: `Q0 + 4N`

Four distinct source points have one `P`-value exactly when they are the
roots of `P-h`.  Requiring `Q` to be constant on those four roots says that
the three nonconstant coefficients of

\[
  Q\bmod(P-h) \tag{6.1}
\]

vanish.  These equations are affine-linear in `(a,b,c,d)`.  Their four
maximal minors are

```text
h^2*(h+1),
-h^2*k*(2*h+1),
h^2*(h*k^2+2*h+1),
-2*h^3*k.
```

They have rank three for every `h!=0`; the localized maximal-minor ideal is
the unit ideal.  Moreover, the first remainder equation is `h` times a monic
linear equation in `h`, so `h` is recovered uniquely from the coefficient
point on the valid open.  The incidence is an affine-line bundle over an
irreducible open in `(k,h)`, and its projection has one irreducible
three-dimensional image.

For

\[
 (k,a,b,c,d)=\left(1,-2,-2,-\frac12,0\right),\qquad h=1, \tag{6.2}
\]

the fiber quartic has discriminant `-279`, the `Q` remainder is `-1/2`, and
the four branch slopes are distinct.  The residual collision factor gives
four nodes.  Sage verifies one length-nine ordinary quadruple point, the
forced length-four cusp, four reduced nodes, and a cyclic complement.  The
generic clean ordinary-quadruple component is excluded.

## 7. Two ordinary triples: `T111^2 + 4N`

For an ordinary triple fiber, choose the fourth, omitted root `e`.  The other
three roots form

\[
  \frac{P(t)-P(e)}{t-e}. \tag{7.1}
\]

They have one common `Q`-value when the linear and quadratic coefficients of
the remainder of `Q` modulo that cubic vanish.  Two omitted roots `u,v` give
four affine-linear equations in `(a,b,c,d)`.  Their determinant factors as

\[
 -(u-v)^2(u^2+ku+1)^2(v^2+kv+1)^2
 (kuv-u-v)\,B, \tag{7.2}
\]

where

\[
 B=\frac{P(u)-P(v)}{u-v}. \tag{7.3}
\]

The factor `B=0` puts the two omitted roots in the same `P`-fiber and is
invalid for two separate triple targets.  The residual factor
`A=kuv-u-v` is generically inconsistent: at `(k,u,v)=(3/2,1,2)`, the
coefficient matrix has rank three while the augmented matrix has rank four.
That one hostile point proves that `A` is not a compatible dominant divisor;
it does not classify compatible lower-dimensional pieces contained in `A=0`.

The complement of all determinant factors is one rational three-dimensional
Cramer graph, symmetric under `u<->v`.  The exact member

\[
 (k,a,b,c,d)=\left(3,5,\frac43,6,5\right) \tag{7.4}
\]

has omitted roots `-1,1`, two distinct triple target values `1/3,-125/3`,
and four residual nodes.  Sage verifies the cusp, two ordinary triple points,
four nodes, and a cyclic complement.  This excludes the component dominating
the dense valid Cramer open.  Residual `A`-supported components and split
boundaries remain open.

## 8. What the checkpoint changes

The conditional delta-ten filtration is now:

| layer | profile count | status |
|---|---:|---|
| clean, dimension five | 1 | excluded |
| codimension-one, dimension four | 2 | both generic components excluded |
| expected codimension two, dimension three | 6 | four identified generic/dominant components excluded; `T112` and `C2+T111` generic exclusions conditionally supported with the propagation qualifications above |
| expected codimension at least three | 136 | not audited component by component |

The 145-profile ledger remains an elimination ledger, not a list of proven
components.  In particular, its 55 overdetermined profiles are not empty by
dimension count alone.  The next exact work is to saturate and decompose:

1. the residual-rank loci in `C2^2`, `C2+T111`, and `T111^2`;
2. the split `k=0,+2,-2` charts for all six profiles;
3. removed pair and `P`-projection/critical-fiber loci that can re-enter
   through a different chart and are not automatically invalid;
4. the remaining `T112` labeling/propagation and mixed-component connected
   equisingular propagation arguments;
5. intersections and specializations of the six identified components; and
6. the fourteen expected codimension-three profiles, followed by their
   endpoints and overdetermined containments.

Even a complete audit of that list would establish only the stated
conditional one-pair, single-three-cycle `A6` result.  It would not settle
the unrestricted generic-degree-six problem.

## 9. Reproduction

Run the dependency-free algebra, exact geometry, stored-relation replay, and
finite permutation censuses with:

```bash
uv run python -m scripts.a6_delta_ten_contact_three
uv run python -m scripts.a6_delta_ten_double_contact
uv run python -m scripts.a6_delta_ten_t112
uv run python -m scripts.a6_delta_ten_contact_triple
uv run python -m scripts.a6_delta_ten_quadruple
uv run python -m scripts.a6_delta_ten_double_triple

uv run pytest -q \
  tests/test_a6_delta_ten_contact_three.py \
  tests/test_a6_delta_ten_double_contact.py \
  tests/test_a6_delta_ten_t112.py \
  tests/test_a6_delta_ten_contact_triple.py \
  tests/test_a6_delta_ten_quadruple.py \
  tests/test_a6_delta_ten_double_triple.py

uv run mypy --no-incremental
```

Regenerate the six primitive curves, singular schemes, raw presentations,
cyclic simplifications, and the `T112` saturated incidence calculation with
Sage 10.8:

```bash
sage tools/check_a6_delta_ten_contact_three.sage
sage tools/check_a6_delta_ten_double_contact.sage
sage tools/check_a6_delta_ten_t112.sage
sage tools/check_a6_delta_ten_contact_triple.sage
sage tools/check_a6_delta_ten_quadruple.sage
sage tools/check_a6_delta_ten_double_triple.sage
```

These Sage commands are manual checkers; GitHub CI runs the Python
certificates and tests, not Sage.  Sage's Zariski--van Kamp implementation
remains a computer-assisted dependency.  The presentation replay itself,
including the complete `40^4` single-three-cycle census, is dependency-free
once the relation words have been generated.  Neither the Python jobs nor the
Sage scripts certify proper Whitney--Thom propagation.
