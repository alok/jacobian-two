# The one-dicritical `A6` local frontier

## Claim boundary

Assume that a hypothetical topological-degree-six plane Keller map has exactly
one dicritical component and has the surviving passport

\[
G=6T15=A_6,
\qquad
(e,d)=(3,1).
\]

Orevkov's exact defect budget gives two positive-jump units.  This note proves
that they cannot split as `1+1`:

\[
\boxed{
\text{there is a unique boundary point }x_*
\text{ with }\mu_{x_*}=5,
\quad
\mu_x=3\text{ for every }x\ne x_*.
}
\tag{0.1}
\]

The cyclic endpoint theorem in
[`one-dicritical-source-smoothness.md`](one-dicritical-source-smoothness.md)
eliminates the normal-surface singularity created by contracting the constant
boundary chain at local degree three or five.  Consequently `L_C` is empty,
the entire affine finite normalization is smooth, and the corresponding
unibranch target singularity at `x_*` has Puiseux pair `(2,5)`.

This is a strict reduction, not an elimination of the `A6` passport.  An exact
permutation fixture shows that the `(2,5)` knot group has the required local
`A5` monodromy and that a collision meridian extends it to global `A6`.

## 1. One possible singular source point

Orevkov's Lemma 2.1 says that every connected component of `L_FC` is a linear
chain

\[
L_C-E,
\]

with exactly one dicritical `E` at its endpoint and a connected, possibly
empty, chain `L_C` of constant components before it.  With only one
dicritical, there is only one such component.  Contracting `L_C` therefore
creates at most one nonsmooth finite point on the image `D` of `E` in the
finite normalization.  Away from that point, both the source surface and `D`
are smooth.  Orevkov chooses the contraction minimally: every point produced
from a nonempty `L_C` chain is a singular source point, since otherwise the
corresponding blowups were unnecessary.  Thus a smooth source point lies off
the contracted chain, where the contraction is an isomorphism and `D` is the
image of the smooth curve `E`.

The one-dicritical defect identity is

\[
\sum_{x\in D}(\mu_x-3)=2. \tag{1.1}
\]

Its only positive partitions are `2` and `1+1`.  The latter would require two
distinct local-degree-four points.  At least one would be a smooth source
point.  The next section proves that such a point cannot exist.

## 2. A smooth point cannot have local degree four

Let

\[
\rho:(X,x)\longrightarrow(\mathbb C^2,0)
\]

be the finite map germ at a smooth point of `D`.  The only ramification prime
is the smooth curve `D`, and its generic index is three.  Orevkov's generic
local form `(s,t) -> (s,t^3)` gives coefficient two along `D`; characteristic
zero tameness and purity rule out any other local ramification divisor.  Thus

\[
\det d\rho=\varepsilon d^2, \tag{2.1}
\]

where `d=0` is a reduced local equation of `D` and `epsilon` is a unit.  Put

\[
m=\mu_x\rho.
\]

### 2.1 Rank one

If `rank(d rho_x)=1`, choose source and target coordinates so that

\[
\rho(X,z)=(X,g(X,z)).
\]

Then `g_z=epsilon*d^2`.  Finiteness gives

\[
m=I_0(X,g)=\operatorname{ord}_z g(0,z).
\]

Therefore

\[
\begin{aligned}
m-1
&=\operatorname{ord}_z g_z(0,z)\\
&=2\operatorname{ord}_z d(0,z)\\
&=2I_0(D,\{X=0\}).
\end{aligned}
\tag{2.2}
\]

Thus `m` is odd, and in particular `m` is not four.

### 2.2 Corank two

Suppose instead that `d rho_x=0` and `m=4`.  Write `rho=(u,v)`.  Both
coordinate functions have order at least two, and local Bezout gives

\[
I_0(u,v)\ge\operatorname{ord}(u)\operatorname{ord}(v)\ge4.
\]

Equality forces both orders to be two and their quadratic initial forms `U,V`
to be coprime.  Hence

\[
[U:V]:\mathbb P^1\longrightarrow\mathbb P^1
\]

is a degree-two morphism.  The quadratic initial form of the Jacobian is its
Wronskian

\[
J(U,V)=U_XV_z-U_zV_X.
\]

Equation (2.1) says that this nonzero quadratic is one doubled line,
`c L^2`.  It would make the entire ramification divisor of a degree-two map of
`P1` a single point of multiplicity two.  This is impossible: a degree-two
map has local ramification contribution at most one at each point, while
Riemann--Hurwitz requires total contribution two.

Thus a smooth source point cannot have local degree four in either rank case.

## 3. The jump is uniquely `mu=5`

If (1.1) had partition `1+1`, its two local-degree-four points could not both
be the at-most-one exceptional source point.  Section 2 excludes the other
one.  Hence only partition `2` remains, proving (0.1).

Several consequences are immediate.

- The point `x_*` cannot share its target value with another normalization
  point: the local degrees would total at least `5+3>6`.
- Its target germ is unibranch.
- That branch germ cannot be smooth.  A smooth branch complement has cyclic
  fundamental group, but a 3-cycle cannot generate a transitive action on the
  connected five-sheet local component.
- Every other noncollision normalization point has zero jump, so Orevkov's
  Lemma 5.2 makes its target branch smooth.
- The normalization collision required by the global topology occurs at a
  different target value.  It consists of two smooth local-degree-three
  branches and is omitted by the affine map.

Thus the finite branch curve has exactly one intrinsic unibranch singularity,
in addition to its multibranch normalization-collision singularities.

## 4. A smooth degree-five source point has rank one

Suppose the source surface is smooth at `x_*` and, toward a contradiction,
`d rho_(x_*)=0`.  Local degree five again forces both coordinate functions to
have order two.  Their quadratic initial forms cannot be coprime, since that
would give intersection four.  They cannot be proportional either: a target
linear combination would then have order at least three, and Bezout would
give intersection at least six.  After linear source and target changes, the
quadratic jets have the form

\[
u=x^2+u_{\ge3},
\qquad
v=xz+v_{\ge3}. \tag{4.1}
\]

Their quadratic Jacobian is `2x^2`, so the tangent line of `D` is `x=0`.
Because the cubic term of the square `epsilon*d^2` is divisible by `x`, the
pure `z^3` coefficient in the cubic Jacobian must vanish.  Direct expansion
of

\[
u_xv_z-u_zv_x
\]

shows that this coefficient is minus three times the pure `z^3` coefficient
of the cubic term of `u`.  Hence that coefficient of `u` is zero.

Blow up the common tangent using `x=zw`.  After removing the common factor
`z^2`, the strict transforms are

\[
\widetilde u=w^2+z,u_3(w,1)+O(z^2),
\qquad
\widetilde v=w+z,v_3(w,1)+O(z^2). \tag{4.2}
\]

The vanishing just proved says that `tilde u` has no linear term at the common
infinitely-near point.  The germ `tilde v=0` is smooth there, so the strict
transform intersection multiplicity is at least two.  The blow-up
intersection formula then gives

\[
I_0(u,v)\ge2\cdot2+2=6,
\]

contradicting local degree five.  Therefore a smooth `x_*` has rank one.

## 5. The smooth subcase forces a `(2,5)` cusp

Use rank-one coordinates

\[
\rho(X,z)=(X,g(X,z)).
\]

Equation (2.2) with `m=5` gives

\[
I_0(D,\{X=0\})=2.
\]

Write the smooth ramification divisor as

\[
D:\quad X=\phi(z),
\qquad
\phi(z)=az^2+O(z^3),
\qquad a\ne0,
\]

and factor the Jacobian as

\[
g_z(X,z)=\varepsilon(X,z)(X-\phi(z))^2.
\]

After replacing the second target coordinate by `g(X,z)-g(X,0)`, the image
of `D` is parametrized by

\[
X=\phi(z),
\]

\[
g(\phi(z),z)=
\int_0^z
\varepsilon(\phi(z),t)
\bigl(\phi(z)-\phi(t)\bigr)^2,dt.
\tag{5.1}
\]

The leading term of (5.1) is

\[
\varepsilon(0,0)a^2z^5
\int_0^1(1-s^2)^2,ds
=\frac8{15}\varepsilon(0,0)a^2z^5.
\]

It is nonzero.  The two target coordinates therefore have orders `(2,5)`.
Since these are coprime, the branch has Puiseux pair `(2,5)` and link the
torus knot `T(2,5)`.

The smoothness of `D` is essential here, not cosmetic.  If it were omitted,
the finite rank-one germ

\[
(X,z)\longmapsto
\left(
X,
\frac{z^5}{5}-\frac23X^3z^3+X^6z
\right)
\]

would be a counterexample to the `(2,5)` conclusion: its Jacobian is

\[
(z^2-X^3)^2,
\]

but its singular critical curve `z^2=X^3`, parametrized by
`(X,z)=(t^2,t^3)`, maps with orders `(2,15)`.  Section 1 is exactly what rules
out this hostile germ in Orevkov's minimal one-chain model.

At this stage of the local argument there are two formal alternatives: a
smooth source with a `(2,5)` cusp, or the unique point produced by contracting
`L_C`.  The Hirzebruch--Jung lift in
[`a6-exceptional-source.md`](a6-exceptional-source.md) eliminates the latter at
local degree five.  Its generalized form in
[`one-dicritical-source-smoothness.md`](one-dicritical-source-smoothness.md)
also excludes a contracted generic degree-three point and proves that `L_C` is
empty.  Thus the smooth `(2,5)` case is mandatory.

## 6. The forced cusp germ actually exists

Orevkov's 2026 classification of smooth finite map germs branched over a
one-Puiseux-pair curve applies to this local block.  There is one hypothesis
worth making explicit.  The local five-sheet monodromy is `A5`, and a
peripheral longitude commutes with its 3-cycle meridian.  In the natural
five-point action,

\[
C_{A_5}((123))=\langle(123)\rangle.
\]

Consequently the longitude fixes each of the meridian's two fixed sheets,
not merely their set.  The two unramified components over the cusp therefore
have tangential degree one.  The ramification component also has tangential
degree one by the one-dicritical passport.  Thus the restriction of the
finite germ to every irreducible component of the cusp pullback is bijective,
which is Orevkov's condition `(A2)`; smoothness of the finite normalization
gives `(A1)`.

Case (b) of Orevkov's Theorem 2 has parameters

\[
m_1=k_2\ell_2+1,
\quad m_2=k_1\ell_1+1,
\quad d_1=k_1m_1,
\quad d_2=k_2m_2,
\quad N=m_1m_2,
\quad n=\ell_1+\ell_2+1.
\]

The specialization

\[
(k_1,k_2;\ell_1,\ell_2)=(2,1;2,0)
\]

gives

\[
(m_1,m_2;d_1,d_2;N,n)=(1,5;2,5;5,3).
\]

This is precisely the degree-five local block, target cusp
`u^2=v^5`, and ramification order three required here.  Orevkov's explicit
normal form in this row has `v=y` and

\[
u_x=\lambda(x^2-y)^2,
\qquad
u(0,y)=0,
\qquad
u(1,1)=1.
\]

Integration and the last normalization give `lambda=15/8` and

\[
F(x,y)=
\left(
\frac{3x^5-10yx^3+15y^2x}{8},
y
\right).
\tag{6.1}
\]

This polynomial map is finite of degree five, and

\[
J(F)=\frac{15}{8}(x^2-y)^2.
\tag{6.2}
\]

Its critical curve `R={y=x^2}` maps by

\[
x\longmapsto (u,v)=(x^5,x^2),
\]

so its image is exactly the `(2,5)` cusp.  More strongly, the complete
pullback factors as

\[
64(u^2-v^5)
=(x^2-y)^3(9x^4-33x^2y+64y^2).
\tag{6.3}
\]

Over `C`, the residual quadratic gives two distinct smooth components
`y=c_+x^2` and `y=c_-x^2`, where

\[
c_\pm=\frac{33\pm9i\sqrt{15}}{128}.
\]

They are nonzero because the parameter polynomial
`64c^2-33c+9` has nonzero constant term, and they are distinct because its
discriminant is `-1215`.  Put

\[
A(c)=\frac{3-10c+15c^2}{8}.
\]

Exact reduction modulo `64c^2-33c+9` gives `A(c)^2=c^5`.  Hence
`beta=A(c)/c^2` satisfies `beta^2=c` and `beta^5=A(c)`.  The parameter
`t=beta*x` identifies each residual restriction with the cusp normalization

\[
t\longmapsto(t^5,t^2).
\]

The critical restriction has the same form with `t=x`.  Thus (6.3) checks
condition `(A2)` directly: the three reduced pullback components all map
bijectively to the cusp, while the critical component occurs with
multiplicity three.

The typed exact checker
[`a6_cusp_germ.py`](../scripts/a6_cusp_germ.py) verifies the theorem
parameters, degree, Jacobian, critical image, full factorization, residual
discriminant, and the congruence `A(c)^2=c^5`.  It also treats the source
coordinate as a quintic over the target:

\[
p(X)=3X^5-10vX^3+15v^2X-8u.
\]

Its exact polynomial discriminant is

\[
\operatorname{disc}_X(p)
=1{,}036{,}800{,}000\,(u^2-v^5)^2.
\tag{6.4}
\]

Thus the square discriminant required by the `A6` passport is realized too;
the local trace discriminant cannot exclude the germ.  A hostile coefficient
perturbation breaks the certificate.

This closes a route, rather than the passport: local analytic classification
cannot eliminate the forced cusp because (6.1) realizes it exactly.  The germ
is not a Keller map—its Jacobian vanishes on `R`—and it supplies only the
five-sheet ramified block.  A hypothetical rank-six finite normalization has
one additional unramified local sheet.  The ramification lies on the boundary
omitted from the original affine `A2`, so there is no conflict with the
constant-Jacobian condition on that open set.

## 7. Abstract monodromy does not eliminate the cusp

The `T(2,5)` knot group has the meridional Artin presentation

\[
\langle r,s\mid rsrsr=srsrs\rangle.
\]

The permutations

\[
r=(345),
\qquad
s=(123)
\]

satisfy the five-braid relation and generate `A5` on the first five sheets.
Both are 3-cycles, exactly as required by the local inertia.  Adjoining the
collision meridian

\[
b=(456)
\]

generates `A6` on all six sheets.  The regression test
`test_a6_torus_2_5_local_monodromy_is_a_hostile_survivor` checks the relation,
the group orders `60` and `360`, the local orbit profile `(5,1)`, and equality
with the certificate's `6T15` group.

This fixture is deliberately hostile: the forced cusp type and the global
group are mutually compatible at the level of abstract local monodromy.  A
valid elimination must use more global Keller geometry.

## 8. Collision contact still carries data

Let two smooth zero-jump normalization branches collide with intersection
multiplicity `q`.  At any smooth source point over one branch, the local
degree-three component is analytically

\[
(u,z)\longmapsto(u,z^3).
\]

Writing the other branch as `v=u^q` times a unit and absorbing a cube root of
that unit, its pullback is

\[
z^3=u^q. \tag{8.1}
\]

It has `gcd(3,q)` puncture branches.  Thus the residual degree-three pullback
has a 3-cycle around the collision when `3` does not divide `q`, but trivial
local monodromy with three puncture branches when `3` divides `q`.  A
normalization collision does not automatically totally ramify the residual
cover.

The polynomial parametrization

\[
t\longmapsto\bigl(t^2,t^5(t^2-1)\bigr)
\]

is a useful geometric countermodel: its image

\[
y^2=x^5(x-1)^2
\]

has a `(2,5)` cusp at `t=0`, a normalization collision `t=1~-1`, and one
place at infinity.  It is not asserted to be a Keller branch curve.  It shows
that the newly forced curve singularities are not contradictory by
themselves.

## 9. The next exact target

The one-dicritical `A6` problem has been reduced to a smooth finite-flat
normalization, one degree-five local `A5` point with a `(2,5)` target cusp,
and at least one separate `3+3` collision.  Every exceptional cyclic-plumbing
source is impossible by the Hirzebruch--Jung lift and parity argument in
[`one-dicritical-source-smoothness.md`](one-dicritical-source-smoothness.md).
Orevkov's 2026 theorem and (6.1) now show that the surviving cusp block is
analytically realizable and unique up to the equivalence recorded in his
classification.  The next attack must couple the cusp and the collision
Kummer laws (8.1) through the
Domrina--Orevkov splice equations or Borisov's canonical and determinant
labels on the compactification graph.

None of the arguments in this note excludes the remaining `S6` passport.

## Primary source

- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587–596.  Lemma 2.1 gives
  the linear `L_C-E` chain, Lemma 4.2 gives the exact defect budget, and Lemma
  5.2 gives smooth local embeddings at zero-jump points.
- S. Yu. Orevkov,
  [“On germs of mappings `C^2 -> C^2`”](https://www.math.univ-toulouse.fr/~orevkov/k-en.pdf)
  (2026).  Conditions `(A1)` and `(A2)` are stated in Section 1; Theorem 2
  classifies the numerical cases; and Section 5 gives the derivative normal
  form used to derive (6.1).

The local-degree-four exclusion, the unique-jump theorem, and the smooth
`(2,5)` classification are derived in this note.  The separate
Hirzebruch--Jung argument proves that the smooth case is mandatory.  No claim
of historical priority is made.
