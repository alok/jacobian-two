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
checkpoint treats the six profiles of expected codimension two.  All six
displayed generic or dominant incidence components are now excluded.  For
`T112` and the mixed contact-plus-triple chart, the separate
[propagation note](a6-delta-ten-propagation.md) works over the smooth labeled
incidence spaces and supplies the previously missing connected-clean-open,
finite-etale labeling, simultaneous-resolution, and proper-isotopy argument.
The three previously exposed residual coefficient-rank factors can no longer
support a threefold: the mixed factor has coefficient rank three on its valid
residual divisor, the two-triple residual is incompatible after exact
localization, and every two-contact residual incidence component has
dimension at most two.  The removed `P`-critical triple charts for `T112` and
the mixed profile likewise have dimension at most two when the critical
branch is immersed.  A generated split-chart ledger now covers 22 generic
allocation types—eleven at `k=0` and eleven `k=+/-2` involution-orbit
types—representing 33 actual rows over all three split values.  It supplies
exact clean witnesses at `k=0,+2` plus full transport to `k=-2`.  Focused audits of the
[contact systems](a6-delta-ten-split-contact-rank.md), the
[triple/mixed systems](a6-delta-ten-split-t112-mixed-rank.md), and the
[global-fiber systems](a6-delta-ten-split-rank-strata.md) now classify every
coefficient-rank stratum in every representative row.  None supports a threefold.  Exact
total-pair and total-two-pair jet comparisons place the eight clean
maximal-rank contact charts in their global excluded components; analogous
labeled-triple, total mixed-base, and global-fiber arguments treat the other
eleven clean split rows.  Exact dominant arcs additionally put the three
prescribed overlap incidences in the algebraic closures of the ordinary
components, but do not yet transport their complement topology.  Five finite
rank-drop orbit representatives—eight actual schemes after `k=-2`
transport—also still need topology.  Denominator
charts, component intersections, lower-dimensional residual taxonomies, and
deeper degenerations also remain in the audit.

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
| `C3 + 7N` | one irreducible threefold on the valid pair chart | generic component and all three clean split rank-open charts excluded | two fixed `k=0` overlap incidences are algebraically contained but their topology and deeper intersections remain open |
| `C2^2 + 6N` | Cramer component dominating the full ordered base, modulo the free swap | dominant component and all five clean split rank-open charts excluded | `k=0` overlap-plus-contact is algebraically contained; its topology and finite residual affine-line fibers of representative ordered lengths `4,6,6` remain open |
| `T112 + 6N` | irreducible labeled affine-line bundle and generically finite threefold image | dominant component and all four clean true-split tangent loci excluded | immersed `P`-critical and deeper intersections remain |
| `C2 + T111 + 5N` | dense rational Cramer component on the displayed open | dominant component and all three clean split rank-open loci excluded | immersed `P`-critical boundary, two representative finite length-four affine-line schemes (three after transport), overlap intersections, and lower-dimensional taxonomy remain |
| `Q0 + 4N` | one irreducible threefold on the valid four-source fiber chart | generic component and clean split fibers excluded | removed critical-fiber, non-clean, and deeper intersections remain |
| `T111^2 + 4N` | rational threefold dominating the dense valid Cramer open | dense component and clean split fibers excluded | same-target, cusp-fiber, non-clean, and deeper intersections remain |

“Excluded” means that the component's clean equisingular open cannot carry
the required connected six-sheet `A6` quotient with single-three-cycle branch
meridians.  It does not mean that the corresponding profile has been removed
scheme-theoretically from every boundary chart.  The exact family-wide
argument for the last two rows is mathematical rather than CAS-certified:
the [propagation note](a6-delta-ten-propagation.md) verifies the hypotheses for
simultaneous embedded resolution and Thom's first isotopy lemma.

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
support another component dominating its divisor.  The follow-up residual
checker closes the higher-codimension threat as well.  In unordered root
coordinates `sigma=u+v`, `pi=uv`, rank-three compatibility is a proper closed
subset of the residual surface and therefore has base dimension at most one.
Every rank-two point is singular on that surface.  Projection resultants put
every positive-dimensional singular locus on `pi=0`, the split curves
`pi=sigma-1`, `pi=-sigma-1`, or the diagonal `4*pi=sigma^2`; exact
restriction makes each invalid or empty on the current chart.  A vertical
projection Gröbner basis is the unit ideal.  Finally, normalized
coefficient-function Wronskians force the two jet rows at a valid contact
root to have rank at least two.  Thus rank-three bases have dimension at most
one, rank-two bases are finite, and rank one is impossible.  Every compatible
residual incidence component has dimension at most two.

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

Therefore that full-base dominant clean component is excluded, and its
residual determinant cannot hide a second threefold.  The compatible
residual curves and points have not been assigned to adjacent profiles.
Split-fiber rank strata and their intersections remain open.

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
the two orientations of its unique tangent pair.  That degree-two
identification is not needed for topology transport: the clean labeled
incidence itself is smooth and irreducible, and its image closure is
irreducible.

The rational member has

\[
 (k,a,b,c,d)=
 \left(2,\frac{2089}{11875},\frac{3542}{2375},
          \frac{1788}{475},\frac{338}{95}\right), \tag{4.3}
\]

with triple sources `-3/5,-2/5,1/5`; the first two have contact exactly two.
The remaining collision packets contain two and four nodes.  Sage finds
singular-scheme lengths `4+6+2+4` and a cyclic complement.  The
[propagation proof](a6-delta-ten-propagation.md) verifies that the exact clean
locus in the smooth labeled incidence is a nonempty irreducible open, labels
its residual nodes after a finite-etale cover, constructs a relative embedded
resolution, and applies proper Whitney--Thom isotopy.  This excludes the
generic clean `T112` chart without using the degree-two labeling claim.  The
manual Sage checker certifies the saturation, singular scheme, and
presentation, not the propagation theorem.  The removed `P`-critical-fiber loci and deeper
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
irreducibility claim is used.  A manual Sage saturation of the coefficient
rank-at-most-two minors by the complete valid-chart localizer is the unit
ideal, with recorded exponent nine.  Hence the coefficient rank is exactly
three everywhere on the valid residual divisor.  Compatibility is a proper
closed subset of that divisor, so its base dimension is at most one and its
affine-line fibers give total incidence dimension at most two.  The remaining
compatible curves and points are not decomposed, but no same-profile
threefold is supported there.

The exact member, after harmless target scaling, is

\[
 X=2t^4+3t^3+2t^2,
 \qquad
 Y=t^9-2t^7+t^6-2t^5. \tag{5.2}
\]

It has the forced cusp, a contact-two fiber, a separate ordinary triple, and
five nodes.  The independent Sage checker verifies singular-scheme lengths
`4+3+4+5`, regenerates the raw presentation, and simplifies it to `Z`.  The
[propagation proof](a6-delta-ten-propagation.md) identifies the displayed
Cramer incidence with a smooth irreducible open of affine three-space and
then verifies finite-etale labeling, relative embedded resolution, and
proper isotopy on its nonempty clean open.  Thus the generic dominant Cramer
chart is excluded.  Neither the Python CI nor the manual Sage checker proves
that theorem step.  The residual factor cannot support a threefold; split
rank strata, overlap intersections, and the lower-dimensional residual
taxonomy remain open.

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
invalid for two separate triple targets.  On the residual factor
`A=kuv-u-v`, substitute `k=(u+v)/(uv)`.  After dividing only the nonzero
root-difference, cusp-fiber, and same-fiber factors and applying a
Rabinowitsch localization, the four normalized augmented determinants
generate the unit ideal.  Thus there is no compatible point anywhere on the
valid residual chart.  The hostile family `k=0`, `v=-u` really is compatible,
but the exact equations put all of it on `B=0`; it is an ordinary-quadruple
boundary rather than two separate triple targets.

The complement of all determinant factors is one rational three-dimensional
Cramer graph, symmetric under `u<->v`.  The exact member

\[
 (k,a,b,c,d)=\left(3,5,\frac43,6,5\right) \tag{7.4}
\]

has omitted roots `-1,1`, two distinct triple target values `1/3,-125/3`,
and four residual nodes.  Sage verifies the cusp, two ordinary triple points,
four nodes, and a cyclic complement.  This excludes the component dominating
the dense valid Cramer open.  The valid residual `A` chart is empty; split
rank strata and component intersections remain open.

## 8. The immersed `P`-critical triple boundary

The fourth-source boundary removed from the ordinary triple-root chart is not
automatically invalid.  If `e` is the critical source, then

\[
 k=-\frac{4e^2+2}{3e},\qquad
 P(t)-P(e)=(t-e)^2
 \left(t^2-\frac{2(1-e^2)}{3e}t+\frac{e^2-1}{3}\right). \tag{8.1}
\]

On

\[
 e(e^2-1)(2e^2-1)(2e^2+1)\ne0,                    \tag{8.2}
\]

the three distinct sources are the critical source and the two roots of the
quadratic.  The critical plane branch is immersed exactly when `Q'(e)!=0`.
An exact saturation of the labeled fourth-root boundary gives three prime
curves, one for each choice of critical label.

For `T112`, the critical branch cannot belong to the tangent pair: the
cross-multiplied tangent equation then reduces to `Q'(e)=0`.  The two
noncritical branches may be tangent.  Their two target-equality equations and
one tangency equation have coefficient rank three throughout the valid
critical curve, so the incidence is an affine-line bundle of dimension two.
The exceptional algebraic value `e^2+2=0` still forces `Q'(e)=100`, rather
than a nonimmersed branch.

For `C2+T111`, the combined system is generically rank four over the
two-dimensional `(e,w)` base.  The common curve of its rank-at-most-two
minors is only the pair-denominator curve; off that curve the rank-at-most-two
base is finite, and the two triple-equality rows always have rank two.  The
denominator curve carries a genuine pair only at `e=+/-1/2`, the split
values.  On the true vertical split chart, `rho=-1/4` changes the profile to
`T112`, while the other determinant roots are augmented-rank inconsistent.
Thus the same-profile critical mixed incidence also has dimension at most
two.  Neither critical boundary supports a hidden codimension-two threefold.

The exact hostile members at `k=-3`, critical source `e=2`, have `Q'(2)=588`
for `T112+6N` and `Q'(2)=504/5` for `C2+T111+5N`.  Their manual Sage singular
schemes have length decompositions `(4,1)+(6,1)+(6,6)` and
`(4,1)+(3,1)+(4,1)+(5,5)`, respectively.  These fixtures prove nonemptiness;
no complement presentation is claimed for the resulting dimension-two
critical pieces.

## 9. The true split-chart generic ledger

At `k=0,+2,-2`, the unordered-pair incidence is reducible.  The genuine
collision degrees are

\[
 (\deg V,\deg W)=(2,8)\quad(k=0),\qquad
 (4,6)\quad(k=\pm2).                                  \tag{9.1}
\]

A finite generator now constructs every generic allocation from these
component budgets, the edge vectors for contacts/triples/quadruples, the
distinct-contact rules, and the unique clean `k=0` overlap rule.  It produces
exactly 22 orbit-type rows: eleven at `k=0` and eleven allocation types at
`k=+/-2`.  Instantiating both signs gives 33 actual rows.
Three tempting `k=0` rows are impossible:

- a triple root on the monic quadratic `V` component;
- two distinct double roots on that same component; and
- an ordinary-triple `V` edge plus a separate contact `V` edge.

Every allowed row has an exact clean witness at `k=0` or `k=2`.  The checker
verifies coefficient and augmented ranks, exact prescribed root
multiplicities, squarefree residual roots and target values, cross-component
target separation, separation from the special targets, and avoidance of
the cusp, extra-critical, and unintended-overlap walls.  Exact full-family
transport

\[
 t\mapsto-t,\qquad Y\mapsto-Y,\qquad
 (a,b,c,d)\mapsto(a,-b,c,-d)                           \tag{9.2}
\]

replays all eleven `k=-2` witnesses and every clean predicate.  The
legitimate `k=0` component overlap is handled explicitly: its generic length
is `1+1`, while the two contact-three allocations have lengths `2+1` and
`1+2`.

Three principal rank-drop systems are globally saturated on their valid
opens: `k=0` graph `C3`, `k=0` vertical/graph double contact, and `k=2`
vertical/vertical double contact.  Their saturation exponents are `3,3,4`.

The focused true-component audits then treat all 22 orbit types, with exact
transport supplying the eleven additional `k=-2` rows.  Every
valid `T112` matrix has rank exactly three.  The mixed residual compatibility
schemes have reduced lengths four and four.  The three `C2^2` residual
determinants have reduced compatible bases of ordered lengths `4,6,6`, with
rank exactly three.  The ordinary-quadruple matrix has rank three on every
valid split fiber, while the two-triple residual is either a same-target
boundary or incompatible.  All other apparent rank drops lie on exact cusp,
diagonal, repeated-root, overlap, or singular-fiber boundaries.  Hence every
valid stratum satisfies

\[
 \dim(\text{base rank locus})+4-\operatorname{rank}\le 2.            \tag{9.3}
\]

The all-allocation aggregate compares exact keys with the generated ledger,
so the conclusion is not inferred from the number 22.  Total-pair jet
derivations and a flat integral ordered two-pair base identify the eight
maximal-rank contact systems with restrictions of the global incidences.
Analogous exact row transformations handle all seven `T112`/mixed systems,
and the global-fiber argument handles the remaining four rows.  Proper
isotopy therefore excludes every clean split rank-open locus.  Exact rational
arcs dominate all three overlap allocation incidences from the ordinary
components, proving algebraic containment but not boundary topology.  No
topology is yet asserted for those overlaps or for the five finite compatible
orbit representatives, which become eight actual schemes over all split
values.

## 10. What the checkpoint changes

The conditional delta-ten filtration is now:

| layer | profile count | status |
|---|---:|---|
| clean, dimension five | 1 | excluded |
| codimension-one, dimension four | 2 | both generic components excluded |
| expected codimension two, dimension three | 6 | all six displayed generic/dominant components excluded; all residual and immersed `P`-critical threefold threats closed; all 22 split orbit types (33 actual rows) classified; every clean rank-open split locus excluded; overlap incidences algebraically contained but their topology, finite rank-drop topology, and deeper boundaries remain |
| expected codimension three, dimension two | 14 | dense clean nonsplit surfaces for `C4+6N`, `C2+C3+5N`, `C2+Q0+2N`, and `C3+T111+4N` are audited; the `C4` residual is bounded below surface dimension, both contact compatibility rank-drop curves are exactly inconsistent, and the `C2+Q0` component has full valid rank four; omitted pair, singular-fiber, same-target, and non-clean loci remain; see the codimension-three note |
| expected codimension at least four | 122 | not audited component by component |

The 145-profile ledger remains an elimination ledger, not a list of proven
components.  In particular, its 55 overdetermined profiles are not empty by
dimension count alone.  The next exact work is to:

1. determine complement topology at the two `C3` overlap fibers and the
   `C2^2` overlap-plus-contact surface, whose algebraic component closure is
   now exact, and treat the five finite compatible orbit representatives
   (eight actual schemes);
2. classify removed pair-denominator and non-clean overlap charts, plus
   intersections of the bounded residual and `P`-critical pieces with other
   walls;
3. assign every compatible lower-dimensional residual curve and point to an
   adjacent collision profile or prove its clean topology directly;
4. audit intersections and specializations of the six identified components;
5. treat the remaining ten expected codimension-three profiles and the
   omitted boundaries of `C4+6N`, `C2+C3+5N`, `C2+Q0+2N`, and
   `C3+T111+4N`; and
6. finish the endpoint and overdetermined containment audit.

Even a complete audit of that list would establish only the stated
conditional one-pair, single-three-cycle `A6` result.  It would not settle
the unrestricted generic-degree-six problem.

## 11. Reproduction

Run the dependency-free algebra, exact geometry, stored-relation replay, and
finite permutation censuses with:

```bash
uv run python -m scripts.a6_delta_ten_contact_three
uv run python -m scripts.a6_delta_ten_double_contact
uv run python -m scripts.a6_delta_ten_t112
uv run python -m scripts.a6_delta_ten_contact_triple
uv run python -m scripts.a6_delta_ten_quadruple
uv run python -m scripts.a6_delta_ten_double_triple
uv run python -m scripts.a6_delta_ten_pcritical_triples
uv run python -m scripts.a6_delta_ten_split_codim_two
uv run python -m scripts.a6_delta_ten_residual_rank

uv run pytest -q \
  tests/test_a6_delta_ten_contact_three.py \
  tests/test_a6_delta_ten_double_contact.py \
  tests/test_a6_delta_ten_t112.py \
  tests/test_a6_delta_ten_contact_triple.py \
  tests/test_a6_delta_ten_quadruple.py \
  tests/test_a6_delta_ten_double_triple.py \
  tests/test_a6_delta_ten_pcritical_triples.py \
  tests/test_a6_delta_ten_split_codim_two.py \
  tests/test_a6_delta_ten_residual_rank.py

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
sage tools/check_a6_delta_ten_pcritical_triples.sage
sage tools/check_a6_delta_ten_contact_triple_residual.sage
```

These Sage commands are manual checkers; GitHub CI runs the Python
certificates and tests, not Sage.  Sage's Zariski--van Kamp implementation
remains a computer-assisted dependency.  The presentation replay itself,
including the complete `40^4` single-three-cycle census, is dependency-free
once the relation words have been generated.  Neither the Python jobs nor the
Sage scripts certify proper Whitney--Thom propagation; the exact hypotheses
and theorem argument are documented separately in
[the propagation note](a6-delta-ten-propagation.md).
