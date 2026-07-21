# Eliminating the exceptional source in the one-dicritical `A6` passport

## Claim boundary

Assume the one-dicritical degree-six `A6` passport from
[`a6-one-dicritical-local.md`](a6-one-dicritical-local.md).  Its unique
positive-jump point has local degree five.  The earlier local analysis left
two alternatives:

1. the source point is smooth, in which case the target branch is a `(2,5)`
   cusp; or
2. the source point is the unique normal-surface singularity obtained by
   contracting Orevkov's linear constant boundary chain.

This note eliminates alternative 2.  It does **not** eliminate the smooth
`(2,5)` cusp, the one-dicritical `A6` passport, generic degree six, or `JC(2)`.

The local statement proved is:

> Let `rho:(X,x)->(C2,0)` be the exceptional local-degree-five finite germ.
> Suppose `(X,x)` is obtained by contracting the connected rational linear
> chain adjacent to the unique dicritical, and the ramification curve `D` is
> the endpoint curvette with generic index three.  Then such a germ cannot
> exist.

## 1. The Hirzebruch--Jung cover

Discarding unnecessary `(-1)` curves turns the contracted exceptional divisor
into a negative-definite linear chain of rational curves.  The
Hirzebruch--Jung classification identifies `(X,x)` as a cyclic quotient
singularity.  Let

\[
q:(\mathbb C^2,0)\longrightarrow(X,x)
\tag{1.1}
\]

be its universal quasi-etale cyclic cover, of order `n>1`.

The strict transform of `D` meets an endpoint of the chain.  In the toric
Hirzebruch--Jung model an endpoint curvette lifts to one coordinate axis.
Equivalently, the endpoint-arrow plumbing has a solid-torus complement; its
lift to the universal `S3` link is again a solid torus, so the lifted
plane-curve link is one unknot.  Thus

\[
\widetilde D=q^{-1}(D)
\]

is one smooth irreducible branch.  Straighten it equivariantly to

\[
\widetilde D=\{x=0\}.
\tag{1.2}
\]

Put

\[
h=\rho\circ q=(u,v):(\mathbb C^2,0)\longrightarrow(\mathbb C^2,0).
\]

The local degree multiplies:

\[
\deg_0 h=5n.
\tag{1.3}
\]

The cover `q` is etale in codimension one.  The ramification divisor of
`rho` has coefficient `e-1=2` along `D` and no other local prime.  Pullback
therefore gives

\[
J_h=\varepsilon x^2
\tag{1.4}
\]

for a unit `epsilon`.

## 2. Quadratic jets force the `A1` quotient

The cyclic deck representation is small: if it fixed a linear form, it would
contain a pseudoreflection and the quotient would be smooth.  Both components
of `h` are deck invariant, so

\[
dh_0=0.
\]

Equation (1.4) has order two.  Consequently both `u` and `v` have nonzero
quadratic initial forms `U,V`, and

\[
J(U,V)=c x^2\ne0.
\tag{2.1}
\]

If `U,V` were coprime, they would define a degree-two map

\[
[U:V]:\mathbb P^1\longrightarrow\mathbb P^1.
\]

Its Wronskian (2.1) would put the entire degree-two ramification divisor over
one point.  That is impossible: a degree-two map has two simple ramification
points by Riemann--Hurwitz.  The quadratics cannot be proportional because
their Wronskian is nonzero.  Hence their gcd is exactly one linear form.
After linear source and target changes,

\[
U=x^2,\qquad V=xy.
\tag{2.2}
\]

This normal form also determines the deck group.  For a deck transformation
`g`, invariance of `U,V` makes their gcd `x` a semi-invariant:

\[
g^*x=\chi_gx.
\]

Write `U=xA`, `V=xB`.  The independent linear forms `A,B` transform by
`chi_g^-1`, so `g` acts on the whole cotangent space as the scalar
`chi_g^-1`.  Since `x` is itself in that cotangent space, it simultaneously
has character `chi_g`.  Therefore

\[
\chi_g^2=1
\]

for every deck transformation.  Faithfulness and `n>1` give

\[
\boxed{n=2,\qquad g=-I.}
\tag{2.3}
\]

Thus the only candidate contracted singularity is `A1`, and both `u` and `v`
are even power series.

## 3. The `A1` parity contradiction

Equations (1.3), (1.4), and (2.2) now read

\[
\deg_0h=10,\qquad
J_h=\varepsilon x^2,\qquad
u_2=x^2,\qquad v_2=xy,
\tag{3.1}
\]

with

\[
h(-x,-y)=h(x,y).
\tag{3.2}
\]

Because the tangent cone of `v=0` is the reduced crossing `xy=0`, Hensel
factorization gives exactly two smooth transverse branches.

### 3.1 The branch tangent to `y=0`

Parametrize it as

\[
\gamma_y(t)=(t,\alpha(t)).
\]

Along this branch `v_y=t+O(t^2)`.  Differentiating `v o gamma_y=0` and using
the Jacobian determinant gives

\[
\frac d{dt}(u\circ\gamma_y)=\frac{J_h\circ\gamma_y}{v_y\circ\gamma_y}.
\]

The numerator has order two and the denominator order one.  Hence this branch
contributes exactly two to `I_0(u,v)`.

### 3.2 The branch tangent to `x=0`

Write

\[
\gamma_x(t)=(\beta(t),t),
\qquad r=\operatorname{ord}\beta.
\]

The series `beta` is not identically zero.  Otherwise this branch would be
the ramification curve `x=0`; the same chain-rule identity below would make
`u` constant on it, so the finite map `h` would contract a curve.

Now `v_x o gamma_x=t+O(t^2)`.  From `v o gamma_x=0`,

\[
v_x\beta'+v_y=0,
\]

and therefore

\[
\frac d{dt}(u\circ\gamma_x)
=-\frac{J_h\circ\gamma_x}{v_x\circ\gamma_x}.
\tag{3.3}
\]

The numerator has exact order `2r` and the denominator exact order one.
Thus this branch contributes `2r`.  Local intersection multiplicity is the
sum over the two reduced branches, so (3.1) gives

\[
10=I_0(u,v)=2+2r,
\qquad r=4.
\tag{3.4}
\]

On the other hand, the involution `-I` preserves the unique branch tangent to
`x=0`.  Its graph obeys

\[
\beta(-t)=-\beta(t).
\tag{3.5}
\]

Therefore `beta` is odd and its finite order must be odd.  This contradicts
`r=4` in (3.4).

The exceptional contracted-source germ cannot exist.

## 4. Exact surviving local frontier

The earlier smooth-source analysis and the contradiction above combine to
give one alternative, not a dichotomy:

\[
\boxed{
\text{the unique positive-jump source point is smooth and rank one, and its
target branch has Puiseux pair }(2,5).}
\tag{4.1}
\]

The exact permutations `(345)` and `(123)` still realize the local `A5`
five-braid relation, and adjoining `(456)` generates global `A6`.  Therefore
(4.1) survives abstract monodromy.  The next contradiction must couple this
smooth cusp to the separate omitted `3+3` normalization collision through
braid, splice, or compactification-graph data.

## Sources

- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587--596.  Lemma 2.1 gives
  the linear boundary chain, Lemma 3.1 gives the generic index-three local
  form, and Lemma 5.2 supplies the endpoint plumbing/local-embedding input.
- Egbert Brieskorn,
  [“Rationale Singularitaten komplexer Flachen”](https://doi.org/10.1007/BF01425318),
  *Inventiones Mathematicae* 4 (1968), 336--358.  The Hirzebruch--Jung
  classification identifies the rational linear-chain singularity with a
  cyclic quotient; the endpoint lift is the standard toric curvette model.

The quadratic-jet character reduction and the final intersection/parity
contradiction are derived in this repository.  No claim of historical
priority is made.
