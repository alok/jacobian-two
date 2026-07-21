# The generic conditional `A6` delta-seven family is cyclic

## Claim boundary

Assume the conditional one-pair setting used in the preceding `A6` audit:

1. the branch normalization is `A1` and is given by a polynomial map;
2. the projective closure has exactly one genuine Puiseux pair at infinity;
3. the only intrinsic finite singularity is the forced `T(2,5)` cusp, whose
   image has no additional normalization preimage; and
4. every other finite singularity is a collision of smooth normalization
   points.

The genus and large-link calculation leaves one candidate at collision delta
seven: affine degrees `(3,10)`, projective infinity pair `(7,10)`.  This note
proves that the **nondegenerate open stratum** of its complete four-parameter
normal-form family has cyclic affine complement.  Hence that open stratum
cannot support the required connected `A6` monodromy with single-3-cycle
meridians.

The propagation from one exact presentation to the entire open stratum uses
proper projective Whitney--Thom equisingular triviality after a finite base
change labels the seven node sections.  That standard topology theorem is an
explicit dependency, not something certified by the symbolic checker.

This is a conditional, computer-assisted generic-stratum exclusion.  It does
not classify the collision-degeneration walls, derive the four hypotheses for
arbitrary Keller branches, eliminate the unrestricted `A6` passport, or prove
the plane Jacobian conjecture.

## 1. Complete degree-`(3,10)` normal form

Move the forced cusp preimage and image to zero and use the same source and
target scaling as in the lower-delta audit.  The first coordinate is

\[
  P=t^2+t^3.
\]

After removing the tangent term, write the monic degree-ten coordinate as

\[
  Q_0=c_4t^4+c_5t^5+c_6t^6+c_7t^7+c_8t^8+c_9t^9+t^{10}.
\]

The polynomial target shear

\[
  Q=Q_0-c_4P^2-(c_6-c_4)P^3
\]

removes both the fourth- and sixth-order terms.  Thus every member is
polynomially equivalent to

\[
  \boxed{
  P=t^2+t^3,
  \qquad
  Q=\alpha t^5+\beta t^7+\gamma t^8+\delta t^9+t^{10},
  \qquad \alpha\ne0.} \tag{1.1}
\]

The exact coefficient change is

\[
\begin{aligned}
\alpha&=c_5-2c_4,\\
\beta&=c_7-3(c_6-c_4),\\
\gamma&=c_8-3(c_6-c_4),\\
\delta&=c_9-(c_6-c_4).
\end{aligned}
\]

The condition `alpha != 0` is exactly what keeps the local characteristic
pair at the origin equal to `(2,5)`.

## 2. Collision septic and the nondegenerate open

For distinct parameters `t,u`, put `s=t+u` and `r=tu`.  Equality of the first
coordinates gives

\[
  r=s^2+s,
  \qquad
  (t-u)^2=-s(3s+4). \tag{2.1}
\]

Modulo that relation, the second divided difference is exactly

\[
  \frac{Q(t)-Q(u)}{t-u}=-s^2H(s),
\]

where

\[
\begin{aligned}
H(s)={}&s^7+6s^6+(9-\gamma+3\delta)s^5
 +(-\beta-2\gamma+9\delta)s^4\\
&+(-4\beta+2\gamma+6\delta-5)s^3
 +(\alpha-3\beta+4\gamma-\delta)s^2\\
&+(\alpha+\beta)s-\alpha. \tag{2.2}
\end{aligned}
\]

The apparent `s^2` factor is the forced cusp diagonal, not an extra
collision.  Three short factors control the invalid and degenerate loci:

\[
\begin{aligned}
C&=\alpha+\beta-\gamma+\delta-1,\\
L&=1215\alpha+756\beta-576\gamma+432\delta-320,\\
T&=\alpha^2-\alpha\beta\gamma+3\alpha\beta\delta-\alpha\beta
   +\alpha\gamma^2-3\alpha\gamma\delta+7\alpha\delta-9\alpha\\
&\quad-2\beta^2\gamma+6\beta^2\delta-6\beta^2
   +5\beta\gamma^2-20\beta\gamma\delta+22\beta\gamma
   +15\beta\delta^2-32\beta\delta+17\beta\\
&\quad-3\gamma^3+15\gamma^2\delta-17\gamma^2
   -21\gamma\delta^2+46\gamma\delta-25\gamma
   +9\delta^3-29\delta^2+31\delta-11.
\end{aligned} \tag{2.3}
\]

The checker proves

\[
\begin{aligned}
H(-1)&=-C,\\
\operatorname{Res}_s(H,-s(3s+4))&=-\alpha L,\\
Q'(-2/3)&=\frac{16L}{19683}. \tag{2.4}
\end{aligned}
\]

Thus `C=0` adds the smooth preimage `t=-1` to the cusp image, while `L=0`
simultaneously creates an extra intrinsic critical point and a diagonal
collision.  Both are invalid under the hypotheses.

Reducing `Q` modulo `t^3+t^2-X`, then taking the resultant of its quadratic
and linear coefficients, gives

\[
  C^2T. \tag{2.5}
\]

Away from the invalid cusp factor, `T=0` is exactly the triple-image wall.

The septic discriminant factors over `QQ` as

\[
  \operatorname{Disc}(H)=C\,G, \tag{2.6}
\]

where the two factors have total degrees one and ten, each with exponent one.
The degree-ten polynomial is intentionally defined by this quotient rather
than printed as a 16,000-character expansion.

For tangent directions, divide

\[
  P'(t)Q'(u)-P'(u)Q'(t)
\]

by `t-u` and reduce by the first-coordinate collision relation.  The exact
answer is

\[
  -s^2(s+1)D(s), \tag{2.7}
\]

not merely `D(s)`.  The factors `s=0,-1` are already excluded on the valid
open.  The displayed checker stores the full degree-seven `D` and proves

\[
  \operatorname{Res}_s(H,D)=-\alpha L G. \tag{2.8}
\]

Consequently the fully nondegenerate parameter set is

\[
  \boxed{
  U=\{(\alpha,\beta,\gamma,\delta)\in\mathbb C^4:
  \alpha LCTG\ne0\}.} \tag{2.9}
\]

Every curve over `U` has the forced cusp, seven distinct transverse nodes,
no triple target, no extra intrinsic critical point, and the fixed infinity
pair `(7,10)`.

## 3. A small exact member has seven nodes

Take

\[
  \boxed{P=t^2+t^3,\qquad Q=2t^5+t^{10}.} \tag{3.1}
\]

At this point

\[
  (\alpha,L,C,T,G)=(2,2110,1,-25,464000000),
\]

so the member lies in `U`.  Its collision and reduced tangency polynomials are

\[
\begin{aligned}
H={}&s^7+6s^6+9s^5-5s^3+2s^2+2s-2,\\
D={}&30s^7+170s^6+250s^5+20s^4-110s^3+30s^2+30s-20.
\end{aligned}
\]

Exact arithmetic gives

\[
\begin{array}{c|r}
\text{quantity}&\text{value}\\ \hline
\operatorname{Disc}(H)&464000000=2^{10}5^6\cdot29\\
\operatorname{Res}(H,-s(3s+4))&-4220\\
\operatorname{Res}(H,D)&-1958080000000\\
Q'(-2/3)&33760/19683\\
Q(-1)&-1.
\end{array} \tag{3.2}
\]

The seven node `X`-coordinates are the roots of

\[
  N(X)=X^7-24X^6-31X^5-30X^4-65X^3+28X^2+18X+2. \tag{3.3}
\]

Here

\[
  \operatorname{Disc}(N)=113281250000000000=2^{10}5^{18}\cdot29,
  \qquad N(0)=2. \tag{3.4}
\]

Thus the seven collision targets are distinct and avoid the cusp image, while
the tangent resultant proves that every collision is transverse.

The genus accounting closes exactly.  A degree-ten plane curve has arithmetic
genus `36`; the `(7,10)` infinity branch contributes `27`, the `T(2,5)` cusp
contributes `2`, and the seven nodes contribute `7`:

\[
  36=27+2+7. \tag{3.5}
\]

The implicit equation is

\[
\begin{aligned}
f(X,Y)={}&-X^{10}-10X^8+15X^6Y-20X^6-4X^5Y-4X^5\\
&+50X^4Y+10X^3Y^2+10X^3Y-25X^2Y^2+Y^3+Y^2.
\end{aligned} \tag{3.6}
\]

The symbolic certificate derives this equation as the parameter resultant
and substitutes `(P,Q)` back into it exactly.

## 4. Its affine complement is cyclic

Sage 10.8's unsimplified affine Zariski--van Kamp calculation uses three
geometric fiber meridians and returns the following signed Tietze relators:

```text
((2,1,-2,-1),
 (2,1,-2,-1),
 (-3,1,3,-1),
 (-3,-2,3,1,-3,2,3,1,-3,2,3,1,-3,-2,3,-1,-3,-2,3,-1),
 (3,2,-3,-2),
 (-2,-1,-2,-1,3,1,2,1,2,-1,-2,-1,-3,1,2,1),
 (-2,-1,3,1),
 (2,1,-2,-1),
 (3,2,-3,-2))
```

Sage's exact simplification isomorphism has codomain

\[
  \langle z\mid\ \rangle\cong\mathbb Z
\]

and sends all three geometric meridians to `z`.  This can also be read from
the raw presentation: three short relators make the generators commute, and
the fourth and seventh independent relators identify them.

The dependency-free checker exhausts all

\[
  40^3=64000
\]

assignments of the three meridians to single 3-cycles in `A6`.  Exactly `40`
satisfy the relators, and every satisfying image has order `3`.  Therefore
there is no `A6` image.

## 5. Propagation and the remaining frontier

The set `U` is the complement of a complex algebraic hypersurface in
`C^4`, hence is path-connected.  Over it, the local embedded types are fixed:
one `T(2,5)` cusp, seven nodes, and one `(7,10)` branch at infinity.  After a
finite base change labels the node sections, apply the proper projective
Whitney--Thom equisingular triviality theorem to

\[
  (\mathbb P^2,\overline{B}_{\lambda}\cup L_\infty),
  \qquad \lambda\in U.
\]

All affine complements over `U` are therefore homeomorphic to the clean
member's complement and have fundamental group `Z`.  Under the stated
hypotheses, a conditional delta-seven survivor must consequently lie on

\[
  G=0\quad\text{or}\quad T=0
\]

inside the valid locus `alpha*L*C != 0`: a repeated-collision or triple-image
wall, including their deeper intersections.  Those walls are not classified
here.

## 6. Reproduction

The dependency-free exact algebra and finite-group replay run with:

```bash
uv run python -m scripts.a6_delta_seven_generic
uv run pytest -q tests/test_a6_delta_seven_generic.py
uv run mypy --no-incremental scripts/a6_delta_seven_generic.py tests/test_a6_delta_seven_generic.py
```

The Sage source calculation, including the raw presentation and explicit
simplification isomorphism, runs with:

```bash
sage tools/check_a6_delta_seven_generic.sage
```

The propagation step uses the plane-curve equisingularity/simultaneous-
resolution theorem followed by proper Whitney--Thom isotopy; see Joseph
Lipman's
[“Equisingularity and simultaneous resolution of singularities”](https://arxiv.org/abs/math/9802010).
