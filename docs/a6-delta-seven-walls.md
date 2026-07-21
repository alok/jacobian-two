# The conditional `A6` delta-seven collision walls are excluded

## Claim boundary

This note continues the [generic delta-seven audit](a6-delta-seven-generic.md).
It works under four additional hypotheses on the unique one-dicritical branch:

1. its normalization is `A1` and is represented by a polynomial map;
2. its projective closure has exactly one genuine Puiseux pair at infinity;
3. its only intrinsic finite singularity is the forced `T(2,5)` cusp, with no
   additional normalization preimage over the cusp image; and
4. every other finite singularity is a collision of smooth normalization
   points.

These hypotheses are not known for an arbitrary degree-six Keller map.  The
unrestricted `A6` passport and the `S6` case remain open.  Everything below is
therefore a conditional, computer-assisted family exclusion, not a proof of
the plane Jacobian conjecture.

Subject to those hypotheses and the topology dependencies stated below, the
audit excludes collision delta seven completely.  The large-link obstruction
then eliminates deltas eight and nine, so the sharpened conditional bound is

\[
  \boxed{\Delta\ge10},
\]

with equality at the coarse link stage forcing affine degrees `(4,9)` and
projective infinity pair `(5,9)`.

After source and target automorphisms, the only degree pair surviving the
large-link test at collision delta seven is represented by

\[
  P(t)=t^2+t^3,
  \qquad
  Q(t)=\alpha t^5+\beta t^7+\gamma t^8+\delta t^9+t^{10},
  \qquad \alpha\ne0. \tag{1}
\]

The generic note derives the monic collision septic

\[
\begin{aligned}
H(s)={}&s^7+6s^6+(9-\gamma+3\delta)s^5
 +(-\beta-2\gamma+9\delta)s^4\\
&+(-4\beta+2\gamma+6\delta-5)s^3
 +(\alpha-3\beta+4\gamma-\delta)s^2\\
&+(\alpha+\beta)s-\alpha,
\end{aligned} \tag{2}
\]

whose roots encode unordered pairs of distinct normalization points with the
same image.  Its discriminant factors as

\[
  \operatorname{Disc}_s(H)=C\,G. \tag{3}
\]

Here `C=0` is invalid because it adds the smooth point `t=-1` over the cusp
image.  The other two relevant divisors are:

- `T=0`, where three normalization points have one common target; and
- `G=0`, where the collision septic acquires a repeated root.

The nondegenerate open `alpha*L*C*T*G != 0` was already excluded by a cyclic
complement computation.  This note records the audit of both remaining
walls.

## 1. Why partitions of seven roots are the right finite ledger

Write

\[
  H=s^7+6s^6+h_5s^5+h_4s^4+h_3s^3+h_2s^2+h_1s+h_0.
\]

The four-parameter family (2) is exactly the affine coefficient slice

\[
\begin{aligned}
h_3-4h_4+10h_5-85&=0,\\
h_0-2h_1+3h_2-11h_4+34h_5-306&=0. \tag{4}
\end{aligned}
\]

Conversely, the parameters are recovered linearly:

\[
\begin{aligned}
\alpha&=-h_0, & \beta&=h_1+h_0,\\
\delta&=(h_4+h_1+h_0+18-2h_5)/3,
&\gamma&=h_4+h_1+h_0+27-3h_5. \tag{5}
\end{aligned}
\]

Thus a root-multiplicity ansatz for `H`, followed by the two exact equations
(4), is an exhaustive algebraic description of each collision stratum.  It
does not depend on a guessed parametrization of the original four
coefficients.

On the valid locus the tangency polynomial `D` satisfies

\[
  (s+1)D-s(3s+4)H'=(9s+10)H. \tag{6}
\]

Consequently, the multiplicity of a valid root of `H` is the contact order of
the corresponding two smooth branches.  If

\[
  e=\sum_i(m_i-1)=7-\#\{\text{distinct roots of }H\},
\]

then the complete partition ledger is:

| Excess | Root partitions |
|---:|---|
| 0 | `(1,1,1,1,1,1,1)` |
| 1 | `(2,1,1,1,1,1)` |
| 2 | `(3,1,1,1,1)`, `(2,2,1,1,1)` |
| 3 | `(4,1,1,1)`, `(3,2,1,1)`, `(2,2,2,1)` |
| 4 | `(5,1,1)`, `(4,2,1)`, `(3,3,1)`, `(3,2,2)` |
| 5 | `(6,1)`, `(5,2)`, `(4,3)` |
| 6 | `(7)` |

The executable partition enumeration in
`scripts/a6_delta_seven_discriminant_wall.py` prevents a stratum from being
silently omitted.

## 2. The triple-image wall `T=0`

An exact three-parameter chart writes the collision polynomial as

\[
  H(s)=\bigl(s^3+2s^2+s+x\bigr)R_4(s). \tag{7}
\]

The three roots of the cubic factor are precisely the three points in one
fiber of the trigonal map `P(t)=t^2+t^3`.  Its discriminant is

\[
  -x(27x-4).
\]

On the valid locus these roots are distinct: `x=0` forces `alpha=0`, while
`27x-4=0` forces `L=0`.  By the standard Zariski--van Kamp construction for
this complete degree-three vertical fiber, the three corresponding geometric
meridians generate the global complement monodromy image.  That generation
statement, not merely the local singularity type, is the decisive global
input.

The pullback of the repeated-collision divisor factors exactly as

\[
  G|_{T=0}=(27x-4)\,D_4\,\rho^2, \tag{8}
\]

where `D4` is the discriminant of the residual quartic and `rho` is its
resultant with the cubic factor.  This also resolves the intersections of the
two walls, including the locus with two triple fibers.

The local braid on the three complete-fiber meridians depends on a contact
profile `(q,q,k)`, with `1 <= q <= k`.  In standard `B3` generators it is

\[
  \Delta_{B_3}^{2q}\sigma_1^{2(k-q)},
  \qquad \Delta_{B_3}=\sigma_1\sigma_2\sigma_1.
\]

The delta-seven budget permits exactly

\[
  (q,k)=(1,1),(1,2),(1,3),(1,4),(1,5),(2,2),(2,3). \tag{9}
\]

For every profile, the checker enumerates all `40^3` assignments of single
three-cycles in `A6`, imposes the exact Hurwitz fixed-point equations, and
finds no generating `A6` triple.  The two stored global presentations—one
ordinary-triple sample and one double-triple sample—each simplify
algebraically to `Z`.  The full-wall exclusion comes from the complete-fiber
generation and Hurwitz argument, not from extrapolating those two samples.

The complete-fiber condition cannot be dropped.  The hostile fixture in
`scripts/a6_delta_seven_triple_wall.py` gives three meridians satisfying the
ordinary-triple braid and generating a group of order nine, then adds a fourth
single-three-cycle meridian and generates all `A6`.  The local relation alone
is therefore not an obstruction.

## 3. Positive-dimensional strata of `G=0`

The generic discriminant divisor and its two generic singular components
give the excess-one and excess-two partitions.  The three excess-three
partitions are one-dimensional.  Exact representatives are:

| Partition | Geometry | Exact sample `(alpha,beta,gamma,delta)` | Constrained `A6` replay |
|---|---|---|---:|
| `(2,1,1,1,1,1)` | one contact two, five nodes | `(-27,0,0,0)` | 40 images, all `C3` |
| `(3,1,1,1,1)` | one contact three, four nodes | `(-4,2,9/4,2)` | 40 images, all `C3` |
| `(2,2,1,1,1)` | two contact twos, three nodes | `(-2,0,0,0)` | 40 images, all `C3` |
| `(4,1,1,1)` | one contact four, three nodes | `(-21,-145/4,-35/2,-45/4)` | 40 images, all `C3` |
| `(3,2,1,1)` | contacts three and two, two nodes | `(-968/19683,6820/6561,18745/6561,56410/19683)` | 40 images, all `C3` |
| `(2,2,2,1)` | three contact twos, one node | `(-3/8,5/2,145/32,105/32)` | 40 images, all `C3` |

Every row has complete genus accounting

\[
  36=27+2+\sum(\text{finite contact orders})+
  \#(\text{remaining nodes}). \tag{10}
\]

These representatives cover components, not merely root partitions.  For
`(2,1,1,1,1,1)`, the incidence equations `H(r)=H'(r)=0` form an affine
rank-two bundle over `r != -1`; the gcd of the maximal minors is `r+1`.  For
`(3,1,1,1,1)`, adjoining `H''(r)=0` gives an affine rank-one bundle over the
same base; the maximal-minor gcd is `2(r+1)^3`.  At the sole rank drop,
`r=-1`, the collision polynomial has `H(-1)=-C`, so the locus is invalid.
On the open where the repeated root is unique, each incidence map is
birational to its image.  The two images are therefore irreducible.

For `(2,2,1,1,1)`, use the elementary symmetric coordinates `a=r+q` and
`b=rq` of the two repeated roots.  The coefficient system for the residual
cubic has determinant

\[
  (a+b+1)D_2(a,b).
\]

The first factor means that one repeated root is `-1`, hence `C=0`.  On
`D2=0`, the compatibility eliminant is

\[
  (a+2)^4(3a+5)^2.
\]

The `a=-2` branch forces `b=1`, so the two nominal double roots merge at
`-1` and again give `C=0`.  The other branch has

\[
  a=-5/3,\qquad b=4/9,
\]

and is a genuine two-double line, but `L=T=0` identically.  Thus the rational
chart loses no valid component.  Away from these boundaries it is a dense
open of an irreducible surface over the `(a,b)`-plane.

For excess three, the contact-four incidence has one valid rational
component.  At its apparent denominator `v=-13/3`, the unrestricted slice
equations reduce to the incompatible equations

\[
  c_0+4c_1=34,
  \qquad 29c_0+116c_1=4736,
\]

so no point is lost.  The unrestricted three-double incidence has four
components: the stored rational component is the only one with a separable
repeated cubic away from `C=T=0`; the other components lie on `C=0`, on the
already excluded complete-fiber wall `T=0`, or on a deeper repeated-root
locus.

The mixed `(3,2,1,1)` incidence ideal has one main one-dimensional prime and
three components supported on `h=-1` or `k=-1`; the latter all lie on `C=0`.
The main component is the graph closure over an absolutely irreducible
genus-two quintic.  Its rational formula for `q` has two valid conjugate
boundary points, characterized exactly by

\[
  261h^2+402h+125=0,
  \qquad 12h+9k+17=0,
  \qquad 264h+387q-40=0.
\]

At both points `alpha`, `C`, `L`, and `T` are nonzero; the triple root, double
root, and two residual roots remain distinct.  These points belong to the
same prime graph closure, not to extra components.  The Sage component
checker verifies the minimal-prime decompositions, absolute irreducibility,
geometric genus two, and every denominator-boundary assertion above.

After invalid and deeper loci are removed, each relevant complex irreducible
curve or bundle base remains connected.  A finite etale base change used to
label collision sections may have several connected components, but every
component surjects onto the connected base.  Proper Whitney--Thom triviality
on those components therefore propagates one representative calculation
throughout each stratum.

For all six representatives, the dependency-free Python modules replay the
stored signed relation words against all `40^3=64,000`
single-three-cycle assignments.  Exactly 40 satisfy each presentation, and
each image is `C3`, so none is the required connected `A6` image.  The checked
Sage script regenerates the three excess-three presentations.  The
excess-one and excess-two presentation extractions remain trusted
computer-assisted input.

These global presentations are essential.  A lone contact-two relation and a
lone contact-three relation each admit 5,760 `A6`-generating triples; two
contact-two relations still admit 1,440.  Local contact type is not enough.

The propagation theorem and the presentation-extraction provenance are
explicit dependencies of the family-wide conclusion.

## 4. Finite excess-four endpoints

There are four possible three-root partitions.  Unrestricted ansätze in the
coefficient slice (4), not just the preceding rational charts, give the
following complete result.

### 4.1 Partition `(5,1,1)`

Put

\[
  H=(s-v)^5\bigl(s^2+(5v+6)s+q\bigr). \tag{11}
\]

Eliminating `q` gives

\[
  -5(v+1)^7(9v^2+51v+34). \tag{12}
\]

The entire discarded fiber `v=-1` has `C=0`.  Introduce the saturation
inverse `u_5` by `u_5(v+1)=1`.  Saturating by `v+1` gives the exact basis

\[
  8u_5-9v-42,
  \qquad 2q+25v+17,
  \qquad P_5(v)=9v^2+51v+34. \tag{13}
\]

The discriminant of `P5` is 1,377, so there are two real endpoints.  Both
are checked directly.

### 4.2 Partition `(4,2,1)`

Put

\[
  H=(s-v)^4(s-w)^2(s+6+4v+2w). \tag{14}
\]

Eliminating `w` gives

\[
-40(v+1)^{17}(3v+1)^2(15v+11)
 P_{42}(v), \tag{15}
\]

where

\[
  P_{42}(v)=3v^4+28v^3+80v^2+68v+17. \tag{16}
\]

The discarded branches are exhaustive:

- `v=-1`: every `w` has `C=0`;
- `v=-1/3`: the branch gcd forces `w=-4/3`, and `T=0`;
- `v=-11/15`: the branch gcd forces `w=-1`, and `C=0`.

Introduce the saturation inverse `u_42` by

\[
  u_{42}(v+1)(3v+1)(15v+11)=1.
\]

After saturating those factors, the basis is

\[
\begin{aligned}
120u_{42}+3931v^3+33809v^2+80057v+30467,\\
2w-3v^3-16v^2-13v,
\qquad P_{42}(v). \tag{17}
\end{aligned}
\]

The discriminant is `-108800`, so `P42` has two real roots and one nonreal
conjugate pair.  The two real embeddings and one nonreal embedding are
checked directly; complex conjugation covers the fourth point.

### 4.3 Partition `(3,3,1)`

Put

\[
  H=(s^2-rs+q)^3(s+3r+6). \tag{18}
\]

Eliminating `q` gives

\[
 -(r+2)^8(3r+5)^3P_{33}(r), \tag{19}
\]

where

\[
  P_{33}(r)=9r^3-93r^2-413r-391. \tag{20}
\]

At the two discarded factors, exact branch gcds force `(r,q)=(-2,1)` or
`(-5/3,2/3)`; both have `C=0`.  After eliminating the auxiliary saturation
inverse, the basis is

\[
  27q+6r^2+38r+34,
  \qquad P_{33}(r). \tag{21}
\]

The cubic discriminant is `-13996800`: there is one real embedding and a
nonreal conjugate pair.  One nonreal embedding is checked directly.  The
finite-etale transport argument in Section 6 carries the nonexistence of the
required finite `A6` quotient to the other two embeddings.  It does not claim
that the three discrete complement groups are isomorphic.

### 4.4 Partition `(3,2,2)`

Let

\[
\sigma=(-6-3r)/2,
\qquad
H=(s-r)^3(s^2-\sigma s+q)^2. \tag{22}
\]

The exact resultant is

\[
 80(r+1)^8(r+2)^3(3r+2)^2(5r+2). \tag{23}
\]

Saturating all four factors gives the unit ideal.  The specialized slice gcds
leave exactly four points:

| `(r,q)` | Result |
|---|---|
| `(-1,1/2)` | degenerates to `(5,2)` and has `C=0` |
| `(-2,1)` | genuine `(3,2,2)`, but lies on `T=0` |
| `(-2/3,1)` | degenerates to `(4,3)` and has `C=0` |
| `(-2/5,7/5)` | genuine `(3,2,2)`, but has `C=0` |

The last point is a useful hostile fixture: a boundary-only parametrization
misses it, although it remains invalid.  The unrestricted ansatz is what
makes the absence claim exhaustive.

## 5. Coarser partitions are absent or invalid

For the two-root partitions `(6,1)`, `(5,2)`, and `(4,3)`, the gcds of the two
slice equations are respectively

\[
  (v+1)^3,
  \qquad 5(v+1)^2,
  \qquad v+1. \tag{24}
\]

Their sole common root is `v=-1`, where `C=0` in every case.  For `(7)`, the
coefficient of `s^6` fixes the root at `-6/7`; the two slice equations evaluate
to

\[
  5/343,
  \qquad 47970/823543, \tag{25}
\]

so no such point lies in the family.

The unrestricted finite-endpoint ansätze in Sections 4--5 use only constant
denominators.  To land the `P5` and `P42` points in the rational contact-four
chart, one divides by `3v+13`; its resultants with `P5` and `P42` are `-162`
and `324`, so none of those endpoints lies on that chart boundary.  The
contact-four, two-double, and mixed-incidence denominator boundaries are
audited separately in Section 3.

## 6. Endpoint topology and arithmetic embeddings

The endpoint polynomials have nine geometric points:

| Component | Geometric points | Stored presentation embeddings | Covered directly or by complex conjugation | Finite-etale arithmetic orbits |
|---|---:|---:|---:|---:|
| `P5` | 2 | 2 | 2 | 1 |
| `P42` | 4 | 3 | 4 | 1 |
| `P33` | 3 | 1 | 2 | 1 |

Each representative presentation uses the three geometric meridians of a
generic vertical fiber.  The exact `40^3` replay has 40 satisfying
assignments, all with generated group of order three, and no `A6` assignment.

The three endpoint polynomials are irreducible over `QQ`; the executable
certificate checks degrees `(2,4,3)` and irreducibility.  Let `K` be any one
of these number fields and let

\[
  U=\mathbb A^2_K\setminus B
\]

be the corresponding branch complement.  A hypothetical surjection from the
topological fundamental group of one complex embedding of `U` to `A6`
algebraizes, by Riemann existence, to a connected finite etale cover.  A
`QQbar` automorphism carrying one embedding of `K` to another transports that
finite etale cover and its labeled `A6` action.  Divisorial tame inertia along
`B` transports as well.  A chosen inertia generator may be replaced by a
unit power in `Zhat`, but an order-three image can only change from `g` to
`g^2`; both are single three-cycles.

Here the endpoint parametrization and `B` are defined over `K`; the image
closure is geometrically irreducible because it is dominated by the
geometrically irreducible normalization `A1`.  Retaining the labeled `A6`
action during transport is essential: it avoids ambiguity from the
exceptional outer automorphisms of `A6`, which need not preserve the two
order-three cycle types.

Applying the inverse field automorphism gives the converse.  Therefore the
existence or nonexistence of the required finite `A6` quotient with
single-three-cycle inertia is constant across an irreducible arithmetic
orbit.  One directly checked embedding of `P5`, `P42`, and `P33` consequently
covers all nine endpoint points.

This is deliberately narrower than a topological-conjugacy claim.  Conjugate
complex varieties need not be homeomorphic and their discrete fundamental
groups need not be isomorphic; only their finite-etale cover data needed here
is transported.

## 7. Consequence and the next coarse frontier

The generic open, the full `T` wall, every positive-dimensional `G` stratum,
and every finite endpoint are now covered.  Under the four standing
hypotheses this eliminates collision delta seven completely.

The exact genus-and-link scan then skips the next two integers:

| Collision delta | Affine degree candidates | Single-three-cycle `A6` survivors |
|---:|---|---|
| 8 | `(2,21)`, `(3,11)` | none |
| 9 | `(2,23)` | none |
| 10 | `(2,25)`, `(3,13)`, `(4,9)`, `(5,7)` | `(4,9)` only |

At `T(2,25)` there are 2,880 generating `A6` pairs, but their meridians have
cycle type `(4,2)` or `(5,1)`, never a single three-cycle.  At `T(4,9)` there
are exactly 720 qualifying pairs.  Thus the next coarse conditional survivor
has

\[
  \Delta=10,
  \qquad (a,d)=(4,9),
  \qquad (d-a,d)=(5,9) \text{ at infinity}. \tag{26}
\]

This section is a lower-bound and degree-filter statement.  The later
[generic delta-ten audit](a6-delta-ten-generic.md) constructs the full
normalized `(4,9)` family and excludes its clean locus by a cyclic-complement
representative.  The later [delta-ten wall audit](a6-delta-ten-walls.md)
also excludes the connected generic opens of its two dominant degeneration
divisors.  It does not construct a cover or Keller map; lower-dimensional
delta-ten strata remain.

## 8. What is exact and what remains theorem-dependent

The following layers are exact and replayable:

- the normal form and collision polynomial identities;
- coefficient-slice inversion, resultants, saturations, Gröbner bases, gcds,
  discriminants, and separation factors;
- the finite Hurwitz and permutation-group enumerations;
- every `40^3` replay against a stored presentation; and
- the degree-candidate and torus-link censuses through delta ten.

The following are explicit computer-assisted or mathematical dependencies:

- Sage's Zariski--van Kamp extraction of the stored complement presentations;
- proper projective Whitney--Thom propagation over connected equisingular
  strata after finite labeling base change;
- the standard local braid identification
  `Delta_B3^(2q) sigma_1^(2(k-q))` and the Zariski--van Kamp fact that the
  three meridians in a complete trigonal fiber generate the global monodromy
  image;
- finite-etale base change, Riemann existence, and tame-inertia compatibility
  for transport across arithmetic endpoint embeddings; and
- the four standing one-pair and finite-singularity hypotheses.

Hostile fixtures prevent the exact local relations, Galois conjugacy, or
singularity counts from being promoted past those dependencies.

For the arithmetic transport layer, use the Riemann existence comparison in
[SGA 1, Expose XII, Theorem 5.1](https://arxiv.org/abs/math/0206203); a modern
finite-cover comparison is also recorded in the
[Stacks Project, Section 58.14](https://stacks.math.columbia.edu/tag/0BTU).
The use here concerns only finite covers and divisorial profinite inertia, not
the full discrete topology.

## 9. Reproduction

Run the exact Python certificates and their finite-group replays with:

```bash
uv run python -m scripts.a6_delta_seven_generic
uv run python -m scripts.a6_delta_seven_triple_wall
uv run python -m scripts.a6_delta_seven_discriminant_wall
uv run python -m scripts.a6_delta_seven_deeper_wall
uv run python -m scripts.a6_delta_seven_finite_wall
uv run python -m scripts.a6_post_delta_seven_frontier

uv run pytest -q \
  tests/test_a6_delta_seven_generic.py \
  tests/test_a6_delta_seven_triple_wall.py \
  tests/test_a6_delta_seven_discriminant_wall.py \
  tests/test_a6_delta_seven_deeper_wall.py \
  tests/test_a6_delta_seven_finite_wall.py \
  tests/test_a6_post_delta_seven_frontier.py

uv run mypy --no-incremental
```

The checked Sage scripts currently include:

```bash
sage tools/check_a6_delta_seven_generic.sage
sage tools/check_a6_delta_seven_components.sage
sage tools/check_a6_delta_seven_deeper_wall.sage
```

The Python tests are dependency-free once the exact signed relation words have
been stored.  The checked Sage scripts regenerate the generic-open and three
excess-three presentations and independently check the positive-dimensional
incidence components.  The stored triple-wall, excess-one/two, and
finite-endpoint presentations are replayed exactly, but their original
Zariski--van Kamp extraction is not currently regenerated by a checked script
in this repository.
