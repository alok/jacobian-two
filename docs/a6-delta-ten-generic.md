# The generic conditional `A6` delta-ten family is cyclic

## Claim boundary

This note continues the
[delta-seven wall audit](a6-delta-seven-walls.md) under the same four
additional hypotheses on the unique one-dicritical branch:

1. its normalization is `A1` and is represented by a polynomial map;
2. its projective closure has exactly one genuine Puiseux pair at infinity;
3. its only intrinsic finite singularity is the forced `T(2,5)` cusp, whose
   image has no additional normalization preimage; and
4. every other finite singularity is a collision of smooth normalization
   points.

The completed delta-seven audit and the exact large-link census leave, at
collision delta ten, only affine degrees `(4,9)` and projective infinity pair
`(5,9)`.  This note derives the complete five-parameter polynomial-normalization
family for that degree pair and proves that its clean stratum has cyclic affine
complement.  It follows that no clean member supports connected `A6` monodromy
when every geometric meridian must be a single three-cycle.

The conclusion is conditional and computer-assisted.  Proper projective
Whitney--Thom triviality is used to propagate one exact complement computation
over the clean parameter locus.  This note does **not** classify the
delta-ten degeneration walls, derive the four hypotheses for arbitrary Keller
branches, eliminate an unrestricted `A6` passport or the `S6` case, construct a
Keller map, or prove the plane Jacobian conjecture.

## 1. Complete degree-`(4,9)` normal form

Move the forced cusp preimage and image to zero.  Before the final source and
target scalings, write the lower-pole coordinate as

\[
  P_0=p_2t^2+p_3t^3+p_4t^4,
  \qquad p_2p_4\ne0.
\]

It has order two at the cusp.  Indeed, a coordinate defining the tangent-line
direction of a `T(2,5)` cusp has order at least five, which is impossible for a
nonzero quartic.

Write the other coordinate as

\[
  Q_0=q_2t^2+q_3t^3+\cdots+q_9t^9,
  \qquad q_9\ne0.
\]

First subtract the tangent term:

\[
  R=Q_0-\frac{q_2}{p_2}P_0.
\]

The `T(2,5)` condition says that the coefficient of `t^3` in `R` vanishes.
If `r_4=[t^4]R`, the polynomial target shear

\[
  \widetilde Q=R-\frac{r_4}{p_2^2}P_0^2
\]

then has the form

\[
  \widetilde Q=A_5t^5+A_6t^6+A_7t^7+A_8t^8+q_9t^9,
  \qquad A_5\ne0. \tag{1.1}
\]

Choose `rho` with `rho^2=p_2/p_4`, rescale the source by `t -> rho*t`,
and scale the two target coordinates.  Every candidate is thereby
polynomially equivalent over `C` to

\[
  \boxed{
  P=t^2+kt^3+t^4,
  \qquad
  Q=at^5+bt^6+ct^7+dt^8+t^9,
  \qquad a\ne0.} \tag{1.2}
\]

Here

\[
  k=\frac{p_3}{\sqrt{p_2p_4}},
\]

and, for the chosen square root,

\[
  a=\frac{A_5}{q_9\rho^4},\qquad
  b=\frac{A_6}{q_9\rho^3},\qquad
  c=\frac{A_7}{q_9\rho^2},\qquad
  d=\frac{A_8}{q_9\rho}. \tag{1.3}
\]

Changing the square-root choice gives the residual involution

\[
  \boxed{
  (k,a,b,c,d)\longmapsto(-k,a,-b,c,-d).} \tag{1.4}
\]

It is induced by `t -> -t` together with `Q -> -Q`.  In particular, `k` is
a genuine normal-form modulus.  Fixing `P=t^2+t^4` or
`P=t^2+t^3+t^4` would omit most of the family.

The pole-degree filtration also explains completeness.  A polynomial target
automorphism preserving the ordered pole degrees `(4,9)` reduces, by the
Jung--van der Kulk degree argument, to translations and a triangular map

\[
  X\longmapsto \lambda X,
  \qquad
  Y\longmapsto \mu Y+h(X),
  \qquad \deg h\le2,
\]

after the two coordinates have been ordered by pole degree.  These are exactly
the translations, tangent shear, quadratic shear, and scalings used above.
Over a nonclosed field the displayed chart may require adjoining
`sqrt(p_2/p_4)`; over `C` the only remaining ambiguity is (1.4).

## 2. The generic collision decic

For distinct normalization parameters `t,u`, put

\[
  s=t+u,\qquad r=tu.
\]

The first divided difference is

\[
  \frac{P(t)-P(u)}{t-u}
  =(2s+k)r-s(s^2+ks+1),
\]

up to an overall sign.  Set

\[
  D=2s+k,\qquad q=s^2+ks+1. \tag{2.1}
\]

Thus the unordered-pair incidence is

\[
  Dr=sq. \tag{2.2}
\]

Let `h_0=1`, `h_1=s`, and

\[
  h_n=sh_{n-1}-rh_{n-2}.
\]

Then `h_n=(t^(n+1)-u^(n+1))/(t-u)`, and the second divided difference is

\[
  G(s,r)=ah_4+bh_5+ch_6+dh_7+h_8. \tag{2.3}
\]

Exact elimination gives

\[
  \operatorname{Res}_r(Dr-sq,G)=s^2H(s), \tag{2.4}
\]

or, on the chart `D != 0`,

\[
  G\left(s,\frac{sq}{D}\right)=\frac{s^2H(s)}{D^4}. \tag{2.5}
\]

The factor `s^2` is the forced cusp diagonal, not an additional collision.
The monic collision decic is

\[
\begin{aligned}
H(s)={}&s^{10}+6ks^9+(4dk-2c+9k^2+8)s^8\\
&+(-4b+10dk^2+8d-3ck+3k^3+18k)s^7\\
&+(-4a-12bk+6dk^3+28dk+3ck^2+2c+6)s^6\\
&+(-16ak-9bk^2-8b+dk^4+16dk^2+16d\\
&\qquad+4ck^3+21ck-3k^3-18k)s^5\\
&+(-17ak^2-16a-2bk^3+2dk^3+4dk\\
&\qquad+ck^4+18ck^2+18c-9k^2-16)s^4\\
&+(-7ak^3-20ak+6bk^2+12b-2dk^2-8d\\
&\qquad+4ck^3+15ck-6k)s^3\\
&+(-ak^4-8ak^2+4a+2bk^3+12bk-4dk\\
&\qquad+3ck^2-2c+1)s^2\\
&+(-ak^3+4ak+3bk^2-ck)s+ak^2.
\end{aligned} \tag{2.6}
\]

For `k notin {0,+2,-2}`, every root of `H` has `D != 0`, because

\[
  \boxed{
  H(-k/2)=\frac{k^2(k^2-4)^4}{1024}.} \tag{2.7}
\]

A root `s_0` then determines one unordered pair: use

\[
  r_0=\frac{s_0(s_0^2+ks_0+1)}{2s_0+k}
\]

and take the roots of `z^2-s_0z+r_0`.  Multiplicities in `H` record the
scheme-theoretic collision budget.  Ten simple roots give ten distinct pairs,
but do not alone prove ten distinct target points: an ordinary triple target,
for example, supplies three different unordered pairs.

## 3. Exact invalidity and tangency identities

Two resultants isolate failures of the standing hypotheses.  Define

\[
  C=\operatorname{Res}_z
  \left(
    z^2+kz+1,
    z^4+dz^3+cz^2+bz+a
  \right), \tag{3.1}
\]

and

\[
  L=\operatorname{Res}_z
  \left(
    4z^2+3kz+2,
    9z^4+8dz^3+7cz^2+6bz+5a
  \right). \tag{3.2}
\]

The factor `C` is the product of the second-coordinate values at the two
nonzero roots of `P(t)=0`; hence `C=0` adds a smooth normalization point over
the cusp image.  The factor `L` is the resultant of `P'(t)/t` and `Q'(t)/t^4`;
hence `L=0` creates an intrinsic finite critical point away from the forced
cusp.  Both loci are invalid under the standing hypotheses.

The exact collision resultants are

\[
  \boxed{
  \operatorname{Res}_s(H,q)=(k^2-4)^4C} \tag{3.3}
\]

and

\[
  \boxed{
  \operatorname{Res}_s\bigl(H,-s(2s^2+3ks+4)\bigr)
  =ak^2(k^2-4)^4L.} \tag{3.4}
\]

Indeed, on the generic chart the pair discriminant is

\[
  (t-u)^2=s^2-4r
  =-\frac{s(2s^2+3ks+4)}{D}. \tag{3.5}
\]

Thus, away from the three split fibers, `a*C*L != 0` ensures that roots of
`H` neither hit the cusp image nor become diagonal pairs.

For tangent directions, let

\[
  K(t,u)=\frac{P'(t)Q'(u)-P'(u)Q'(t)}{t-u}.
\]

On (2.2), exact reduction has the form

\[
  K=\frac{s^2qT(s)}{D^5}, \tag{3.6}
\]

where `T` is characterized without printing its large expansion by the
checked syzygy

\[
\begin{aligned}
qT
&-s(2s^2+3ks+4)D H'\\
&=(9k^2s-2ks^2+10k-4s^3-16s)H.
\end{aligned} \tag{3.7}
\]

Taking resultants in (3.7) gives the useful derived identity

\[
  C\operatorname{Res}_s(H,T)
  =-ak^4(k^2-4)^4L\operatorname{Disc}_s(H). \tag{3.8}
\]

Consequently, on the generic valid chart, a nonzero decic discriminant makes
all ten pair collisions transverse.

The common first-coordinate value of a collision pair is

\[
  X(s)=-\frac{s(s+k)(s^2+ks+1)^2}{(2s+k)^2}. \tag{3.9}
\]

Therefore distinct collision targets can be certified by the discriminant in
`X` of

\[
  \operatorname{Res}_s
  \left(
    H(s),
    (2s+k)^2X+s(s+k)(s^2+ks+1)^2
  \right). \tag{3.10}
\]

This resultant has degree ten in `X` and leading coefficient
`k^4*(k^2-4)^8` on the generic chart.

## 4. The three split incidence fibers are legitimate

The values `k=0,+2,-2` are not invalid.  At them, the first-coordinate pair
incidence becomes reducible, so division by `D` loses a component.  The exact
component formulas recover the full collision budget.

### 4.1 The fiber `k=0`: degrees `2+8`

Equation (2.2) factors as

\[
  s(2r-s^2-1)=0. \tag{4.1}
\]

On the component `s=0`, the second divided difference is

\[
  G(0,r)=r^2V_0(r),
  \qquad
  V_0(r)=r^2-cr+a. \tag{4.2}
\]

The factor `r^2` is the forced diagonal at the cusp.  The legitimate vertical
component therefore contributes a degree-two collision polynomial.

On the graph component `r=(s^2+1)/2`, one has

\[
  16G=W_0(s), \tag{4.3}
\]

where

\[
\begin{aligned}
W_0(s)={}&s^8+(8-2c)s^6+(-4b+8d)s^5\\
&+(-4a+2c+6)s^4+(-8b+16d)s^3\\
&+(-16a+18c-16)s^2+(12b-8d)s\\
&+4a-2c+1.
\end{aligned} \tag{4.4}
\]

Thus the split collision budget is `2+8=10`.  The generic decic specializes
to

\[
  H|_{k=0}=s^2W_0(s). \tag{4.5}
\]

The displayed `s^2` in (4.5) is a chart-cancellation artifact, not a substitute
for `V_0`.  The two incidence components meet at `(s,r)=(0,1/2)`, and

\[
  W_0(0)=4V_0(1/2).
\]

If that common point is also a collision, the same pair lies on both
components and belongs to a degeneration stratum rather than the clean open.

### 4.2 The fibers `k=2*epsilon`: degrees `4+6`

Let `epsilon` be `+1` or `-1`.  Equation (2.2) factors as

\[
  (s+\varepsilon)\bigl(2r-s(s+\varepsilon)\bigr)=0. \tag{4.6}
\]

On the involution component `s=-epsilon`, the degree-four collision polynomial
is

\[
\begin{aligned}
V_\varepsilon(r)={}&r^4+(4\varepsilon d-c-10)r^3\\
&+(a-3\varepsilon b-10\varepsilon d+6c+15)r^2\\
&+(-3a+4\varepsilon b+6\varepsilon d-5c-7)r\\
&+a-\varepsilon b-\varepsilon d+c+1.
\end{aligned} \tag{4.7}
\]

On the graph component `r=s(s+epsilon)/2`, the exact identity is

\[
  16G=s^2W_\varepsilon(s), \tag{4.8}
\]

where

\[
\begin{aligned}
W_\varepsilon(s)={}&s^6+8\varepsilon s^5
 +(6-2c+8\varepsilon d)s^4\\
&+(-4b+16d+2\varepsilon c-16\varepsilon)s^3\\
&+(-4a-8\varepsilon b-8\varepsilon d+18c+1)s^2\\
&+(-16\varepsilon a+12b-2\varepsilon c)s+4a.
\end{aligned} \tag{4.9}
\]

After removing the forced `s^2` cusp diagonal, the graph component has degree
six.  Hence the split collision budget is `4+6=10`.  The generic decic obeys

\[
  H|_{k=2\varepsilon}=(s+\varepsilon)^4W_\varepsilon(s). \tag{4.10}
\]

Again, the fourth power in (4.10) records cancellation of the missing
incidence component; it is not four new collision pairs.  The components meet
at `(s,r)=(-epsilon,0)`, the possible extra preimage of the cusp image, and

\[
  W_\varepsilon(-\varepsilon)=16V_\varepsilon(0).
\]

Points on these three fibers belong to the clean family whenever the split
polynomials have the required distinctness, separation, and transversality.
The equations `k=0,+2,-2` themselves must never be listed as invalidity or
degeneration factors.

## 5. Exact genus accounting

The pole degrees `4` and `9` are coprime.  If the parametrization factored
through a nontrivial polynomial covering of degree `e`, then `e` would divide
both pole degrees.  Hence `e=1`, and `P1` is the normalization.

A degree-nine plane curve has arithmetic genus

\[
  p_a=\frac{(9-1)(9-2)}2=28. \tag{5.1}
\]

At infinity, with `w=1/t`,

\[
  \left(\frac P Q,\frac1Q\right)
  =\left(w^5\cdot\text{unit},w^9\cdot\text{unit}\right).
\]

The unique infinity branch therefore has pair `(5,9)` and contributes

\[
  \delta_\infty=\frac{(5-1)(9-1)}2=16. \tag{5.2}
\]

The forced `T(2,5)` cusp contributes two.  The remaining finite collision
delta is consequently

\[
  \boxed{28-16-2=10.} \tag{5.3}
\]

This independently explains the generic decic and the exceptional `8+2` and
`6+4` splittings.

## 6. A canonical clean member

Take

\[
  \boxed{
  P=t^2+t^3+t^4,
  \qquad
  Q=t^5+t^9.} \tag{6.1}
\]

Thus `(k,a,b,c,d)=(1,1,0,0,0)`.  Its implicit equation, derived exactly as
`Res_t(P-X,Q-Y)`, is

\[
\begin{aligned}
F(X,Y)={}&-X^9-4X^8-8X^7+21X^6Y-4X^6\\
&+41X^5Y-X^5+27X^4Y^2+20X^4Y\\
&+22X^3Y^2+5X^3Y+9X^2Y^3+16X^2Y^2\\
&+5XY^3+4XY^2+Y^4-Y^3+Y^2.
\end{aligned} \tag{6.2}
\]

The collision decic is

\[
\begin{aligned}
H={}&s^{10}+6s^9+17s^8+21s^7+2s^6-37s^5\\
&-58s^4-33s^3-4s^2+3s+1.
\end{aligned} \tag{6.3}
\]

The exact symbolic certificate gives

\[
\begin{array}{c|r}
\text{quantity}&\text{value}\\ \hline
\operatorname{Disc}(H)&-407351195013757923\\
\operatorname{Res}(H,2s+1)&81\\
\operatorname{Res}(H,-s(2s^2+3s+4))&335421\\
C&1\\
L&4141\\
\operatorname{Res}(H,T)&136634145182709696290583.
\end{array} \tag{6.4}
\]

In particular, the ten unordered pairs are distinct, off the denominator and
diagonal loci, away from the cusp image, and transverse.

Their distinct `X`-coordinates are the roots of

\[
\begin{aligned}
N(X)={}&X^{10}-10X^9+3X^8+191X^7-226X^6-712X^5\\
&-668X^4-331X^3-95X^2-15X-1,
\end{aligned} \tag{6.5}
\]

because the collision-coordinate resultant is `6561*N(X)`.  Moreover,

\[
\begin{aligned}
\operatorname{Disc}(N)
&=-766610929107006671875\\
&=-5^6\,11^2\,23^3\,89\,293\,1277993\ne0.
\end{aligned} \tag{6.6}
\]

Thus no two collision pairs share a target.  The curve has exactly one
`T(2,5)` cusp, ten reduced transverse nodes, and the `(5,9)` branch at
infinity.  Its complete genus accounting is

\[
  28=\underbrace{2}_{T(2,5)}
    +\underbrace{10}_{\text{nodes}}
    +\underbrace{16}_{(5,9)\text{ at infinity}}. \tag{6.7}
\]

The checked Sage singular-scheme computation independently finds that
`(F,F_X,F_Y)` has length fourteen and radical length eleven.  Its primary
decomposition consists of a length-four component at the origin and one
reduced length-ten node component.  The node component contains `N(X)`, while
the radical of the whole affine singular scheme contains `X*N(X)`.

## 7. The canonical affine complement is cyclic

Sage 10.8's unsimplified affine Zariski--van Kamp calculation for (6.2) has
four geometric vertical-fiber meridians and thirteen signed Tietze relators.
The exact relation words are stored as `SAMPLE_RELATIONS` in
`scripts/a6_delta_ten_generic.py`.

Sage's explicit simplification isomorphism has codomain

\[
  \langle z\mid\ \rangle\cong\mathbb Z, \tag{7.1}
\]

sends each of the four geometric meridians to `z`, and has a checked section
sending `z` back to the first raw generator.  Therefore

\[
  \boxed{\pi_1(\mathbb A^2\setminus V(F))\cong\mathbb Z.} \tag{7.2}
\]

The dependency-free finite-group checker nevertheless replays the raw
presentation directly.  It exhausts all

\[
  40^4=2{,}560{,}000
\]

assignments of the four generators to the forty single three-cycles in `A6`.
Exactly forty assignments satisfy all thirteen relators, and every satisfying
image has order three.  No satisfying assignment generates `A6`.

This census excludes the prescribed single-three-cycle `A6` passport.  It
does not exclude arbitrary conjugacy classes in `A6`, and it says nothing by
itself about unrestricted `S6` monodromy.

## 8. The clean parameter locus is connected

Let `U` be the locus in `C^5` with coordinates `(k,a,b,c,d)` on which:

1. `a*C*L != 0`;
2. the off-diagonal unordered-pair scheme is reduced of length ten;
3. its ten points have distinct target images; and
4. every corresponding pair of normalization branches is transverse.

The fixed leading terms make the cusp and infinity types constant once
`a != 0`.  The failures in the remaining conditions are algebraic: on the
generic incidence chart they are detected by the discriminants and resultants
in Sections 2 and 3, and on `k=0,+2,-2` they are detected by the finitely many
component discriminants, pairwise resultants, and tangent resultants built
from `V_0,W_0,V_epsilon,W_epsilon`.  Equivalently, these are the discriminant
and noninjectivity loci of the relative finite double-point scheme.  Hence
`U` is Zariski open.

The canonical member (6.1) lies in `U`, so `U` is nonempty.  Since `C^5` is
irreducible, every nonempty Zariski-open subset is irreducible and connected;
its complex analytic space is path-connected.  Notice that this definition
does not remove `k=0,+2,-2`: any clean points on those split incidence fibers
belong to the same `U`.

Over `U`, the projective curve has constant embedded singularity data: one
`T(2,5)` cusp, ten nodes, and one `(5,9)` branch at infinity.  After a finite
base change labels the ten node sections, apply simultaneous resolution and
proper projective Whitney--Thom isotopy to

\[
  (\mathbb P^2,\overline B_\lambda\cup L_\infty),
  \qquad \lambda\in U. \tag{8.1}
\]

All affine complements over `U` are therefore homeomorphic to the canonical
complement and have fundamental group `Z`.  It follows that the entire clean
stratum is excluded from the required connected `A6` passport.

This propagation is a theorem dependency, not a consequence of the
singularity count alone.  A reference for the equisingularity and
simultaneous-resolution input is Joseph Lipman's
[“Equisingularity and simultaneous resolution of singularities”](https://arxiv.org/abs/math/9802010).

## 9. What is exact, what is regenerated, and what remains

The following layers are exact and dependency-free once the stored relation
words are present:

- the complete normal form and its residual involution;
- the pair-incidence resultant, collision decic, denominator identity, and
  exceptional `8+2` and `6+4` component identities;
- the cusp-image, extra-critical, diagonal, tangency, and target-separation
  checks;
- the implicit equation and all exact geometry of the canonical member; and
- the exhaustive `40^4` replay of the stored four-meridian presentation.

The Python replay does **not** derive the relation words.  In this case,
however, the checked Sage source does more than trust an old transcript:
`tools/check_a6_delta_ten_generic.sage` reruns the affine Zariski--van Kamp
calculation, asserts that the newly produced raw relations equal the stored
tuple exactly, checks the simplification isomorphism and its section, and
independently checks the affine singular-scheme primary decomposition.  The
presentation is therefore regenerated by checked source, while the
mathematical correctness of Sage's Zariski--van Kamp implementation remains
an explicit computer-assisted dependency.

The other explicit dependency is the proper Whitney--Thom propagation in
Section 8.  The four standing branch hypotheses remain assumptions.

The follow-up [delta-ten wall audit](a6-delta-ten-walls.md) proves that the
contact-two and ordinary-triple loci are the two irreducible dominant
degeneration divisors and excludes their connected generic equisingular opens
by exact cyclic-complement representatives.  What remains at collision delta
ten is therefore lower-dimensional: deeper repeated roots, nonordinary
multiple fibers, intersections of the two divisors, split-fiber degenerations,
and finite endpoints.  The split equations `k=0,+2,-2` are not themselves
walls; only degenerations detected within their two-component incidence
formulas belong to this remaining audit.

Even a complete audit of those walls would still be conditional on the four
standing hypotheses and on the prescribed single-three-cycle `A6` passport.
The unrestricted `A6` and `S6` cases, and therefore `JC(2)`, remain open.

## 10. Reproduction

Run the exact symbolic certificates and finite-group replay with:

```bash
uv run python -m scripts.a6_delta_ten_generic
uv run pytest -q tests/test_a6_delta_ten_generic.py
uv run mypy --no-incremental \
  scripts/a6_delta_ten_generic.py \
  tests/test_a6_delta_ten_generic.py
```

Regenerate the raw presentation, cyclic simplification, and singular-scheme
decomposition with:

```bash
sage tools/check_a6_delta_ten_generic.sage
```

The ordinary Python CI runs the symbolic certificate, presentation replay,
tests, and type checker.  Sage is a separate checked regeneration step.
