# The conditional `A6` delta-five family is exhausted

## Claim boundary

The conditional one-pair audit originally found collision delta five as the
first remaining `A6` budget and forced affine pole degrees `(3,8)` at
equality.  This note exhausts that family.  Using exact collision algebra,
equisingular-family topology, Sage Zariski--van Kamp presentations, and a
dependency-free finite-group replay, it proves that:

- the stored generic three-parameter presentation simplifies in Sage to a
  cyclic affine complement;
- the stored ordinary-triple presentation simplifies in Sage to a cyclic
  complement;
- the stored contact-two presentation simplifies in Sage to a cyclic
  complement;
- the codimension-two residual has four rational irreducible curves and only
  four valid exceptional parameter points; and
- none of the eight residual presentations admits an `A6` image when every
  geometric meridian is a single 3-cycle.

Consequently collision delta five is impossible under the hypotheses.  The
delta-six large link is impossible as well; at delta seven only affine degrees
`(3,10)` pass the link test.  Thus the sharpened conditional conclusion is

\[
  \boxed{\Delta\ge7,\qquad
  \Delta=7\Longrightarrow(a,d)=(3,10).}
\]

This remains a conditional, computer-assisted theorem.  The hypotheses that
the branch has polynomial normalization, exactly one genuine infinity pair,
only the forced intrinsic `T(2,5)` cusp, and otherwise only smooth-branch
normalization collisions have not been derived for every Keller branch.  The
result does not eliminate the unrestricted `A6` passport or prove the plane
Jacobian conjecture.

The evidence has three deliberately separate layers.  The collision
factorizations, primary-component memberships, exceptional-factor audit, and
finite permutation censuses are exact algebraic computations.  Sage's
Zariski--van Kamp implementation supplies the eight stored presentations.
Transport from one presentation to every point of an open stratum is the
computer-assisted topological layer: it additionally uses simultaneous
resolution and proper Whitney--Thom isotopy.  The dependency-free replay
checks the stored relators exactly, but does not independently derive those
relators or the equisingular transport theorem.

## 1. Three-parameter normal form

Move the forced finite `T(2,5)` cusp preimage and image to zero.  The same
normalization used in the lower-delta audit makes the degree-three coordinate

\[
  P=t^2+t^3.
\]

A degree-eight second coordinate with that cusp has the form

\[
  Q=c_4t^4+c_5t^5+c_6t^6+c_7t^7+t^8,
  \qquad c_5-2c_4\ne0.
\]

The triangular target automorphism `Y -> Y-c4*X^2` removes the fourth-order
term without changing the affine complement, the pole degree, or the
infinity pair.  Thus every curve in the conditional family is equivalent to

\[
  \boxed{
    P=t^2+t^3,
    \qquad
    Q=\alpha t^5+\beta t^6+\gamma t^7+t^8,
    \qquad \alpha\ne0.
  } \tag{1.1}
\]

## 2. Exact collision stratification

For distinct parameters `t,u`, put `s=t+u` and `r=tu`.  Equality of the first
coordinates gives `r=s^2+s`.  Modulo that relation, equality of the second
coordinates is

\[
  s^2H_{\alpha,\beta,\gamma}(s)=0,
\]

where the non-diagonal collision polynomial is

\[
\begin{aligned}
H(s)={}&s^5+(\gamma+2)s^4+(4\gamma-2)s^3\\
&+(-\alpha+2\beta+3\gamma-4)s^2\\
&+(-\alpha+3\beta-\gamma)s+\alpha.
\end{aligned} \tag{2.1}
\]

Set

\[
  L=135\alpha-108\beta+84\gamma-64.
\]

The exact pair and derivative identities are

\[
  \operatorname{Res}\bigl(H,-s(3s+4)\bigr)=-\alpha L,
  \qquad
  Q'(-2/3)=\frac{16L}{2187}. \tag{2.2}
\]

Thus `alpha=0` destroys the `T(2,5)` cusp, while `L=0` simultaneously creates
an extra intrinsic ramification and a diagonal collision.  Both are invalid
under the hypotheses.

The reduced tangent determinant `D(s)` satisfies

\[
  \operatorname{Res}(H,D)
  =\alpha^3L\operatorname{Disc}(H). \tag{2.3}
\]

The discriminant factors exactly as

\[
  \operatorname{Disc}(H)
  =(\alpha-\beta+\gamma-1)G(\alpha,\beta,\gamma), \tag{2.4}
\]

where exact factorization gives `G` total degree seven and irreducible over
the algebraic closure.  The linear factor is invalid rather than a legitimate
tangency wall:

\[
  H(-1)=\alpha-\beta+\gamma-1.
\]

Its vanishing means `nu(0)=nu(-1)`, so the forced cusp image gains the smooth
normalization preimage `t=-1`.

The geometric irreducibility assertion has an exact rational-incidence
certificate; it is not inferred merely from irreducibility over `QQ`.  Away
from `s=-3`, the equations `H(s)=H'(s)=0` solve as

\[
\alpha=
\frac{s^2(4\gamma s^2+13\gamma s+11\gamma+6s^3+14s^2-12)}{s+3},
\]

\[
\beta=
\frac{2\gamma s^4+5\gamma s^3-\gamma s^2-7\gamma s+\gamma
      +3s^5+5s^4-6s^3-6s^2+8s}{s+3}.
\]

Both substituted equations vanish identically.  Eliminating `s` from their
cleared numerators gives `24*(gamma-6)*G`; the extra `gamma=6` line is the
`s=-3` limit of the same incidence.  Thus an irreducible `(s,gamma)` parameter
space dominates all of `G`, proving geometric irreducibility.  The ordinary
CI checker verifies these rational identities, and the optional Sage checker
also factors `G` directly over `QQbar`.

A separate collision phenomenon occurs when three normalization parameters
have the same image.  Reducing `Q` modulo `t^3+t^2-X` gives the exact divisor

\[
\begin{aligned}
T={}&\alpha\gamma-\alpha-2\beta\gamma+3\beta\\
&+2\gamma^2-5\gamma+3.
\end{aligned} \tag{2.5}
\]

Indeed, the resultant of the nonconstant coefficients of that remainder is

\[
  -(\alpha-\beta+\gamma-1)^2T.
\]

The expression

\[
  T=(\gamma-1)\alpha-(2\gamma-3)\beta
    +(\gamma-1)(2\gamma-3)
\]

is primitive and linear in `(alpha,beta)`, so `T` is irreducible.

On the connected open set

\[
  \alpha L\operatorname{Disc}(H)T\ne0, \tag{2.6}
\]

the finite curve has exactly the forced cusp and five distinct transverse
nodes.  Its infinity chart is uniformly

\[
  \left(\frac P Q,\frac1Q\right)
  =\left(
    w^5\frac{1+w}{1+\gamma w+\beta w^2+\alpha w^3},
    w^8\frac1{1+\gamma w+\beta w^2+\alpha w^3}
  \right),
\]

so the `(5,8)` embedded infinity type is constant.

## 3. A clean curve is locally perfect; Sage simplifies its complement to `Z`

Take

\[
  (\alpha,\beta,\gamma)=(1,1,0),
  \qquad
  Q=t^5+t^6+t^8. \tag{3.1}
\]

Its collision polynomial is

\[
  H=s^5+2s^4-2s^3-3s^2+2s+1.
\]

The exact checks

\[
  \operatorname{Disc}(H)=-4903,
  \qquad
  \operatorname{Res}(H,-s(3s+4))=37,
\]

and

\[
  \operatorname{Res}(H,D)=181411=37\cdot4903
\]

give five distinct transverse pairs.  Their five `X`-coordinates are the
distinct roots of

\[
  X^5+3X^4-3X^3-35X^2+12X-1,
\]

whose discriminant is `-5^6*4903`; its value at zero is `-1`, so no node
shares the cusp image.  Also

\[
  Q'(-2/3)=-\frac{592}{2187}\ne0,
\]

so the cusp is the only intrinsic finite singularity.

The implicit equation is

\[
\begin{aligned}
b={}&-X^8-5X^7-6X^6+8X^5Y+X^5+19X^4Y-3X^3Y\\
&-15X^2Y^2+9XY^2+Y^3-Y^2.
\end{aligned}
\]

At infinity the pair is `(5,8)` and contributes delta 14.  The complete genus
accounting is

\[
  21=p_a
  =\underbrace{2}_{T(2,5)}
  +\underbrace{5}_{\text{nodes}}
  +\underbrace{14}_{(5,8)\text{ at infinity}}.
\]

Thus this is a clean realization of every conditional geometric requirement.
Nevertheless Sage 10.8 gives

\[
  \pi_1(\mathbb A^2-B)\cong\mathbb Z,
\]

with all three geometric fiber meridians mapping to the same generator.  An
independent replay of all `40^3=64000` single-3-cycle assignments finds only
40 homomorphisms, all with image `C3`, and none with image `A6`.

The large link `T(3,8)` admitted 720 correct `A6` pairs.  This calculation
shows that none factors through the actual global complement of (3.1).

## 4. The generic family is excluded

Over the open set (2.6), the finite cusp and five node sections, their
intersection data, and the infinity branch all have constant embedded type.
After a finite base change labels the collision sections, the standard
plane-curve equisingularity/simultaneous-resolution theorem supplies a
Whitney stratification.  Proper Thom isotopy then transports the projective
pair, and hence the affine complement group, throughout the connected open
set.  The clean cyclic computation therefore excludes every generic member.

This topological-triviality theorem is an explicit dependency, just as in the
lower-delta family audit; connectedness plus a singularity count alone would
not be sufficient.

## 5. Both admissible codimension-one walls are excluded

After removing `alpha=0`, `L=0`, and the invalid cusp-collision factor, there
are two admissible codimension-one walls.

On the ordinary-triple divisor `T=0`, take

\[
  (\alpha,\beta,\gamma)=(3,0,0),
  \qquad
  H=(s^2-3)(s^3+2s^2+s-1).
\]

The cubic factor gives the three pair sums of one ordinary triple point; the
remaining two pairs are nodes.  Sage's presentation simplifier again reports
complement group `Z`, and the decisive exact `40^3` replay gives only 40
cyclic images.

On the tangency divisor `G=0`, take

\[
  (\alpha,\beta,\gamma)=(-4,1,3),
  \qquad
  H=(s+2)^2(s^3+s^2+2s-1).
\]

This has one contact-two pair and three nodes.  Its complement and
representation census are again `Z`, 40 images of order three, and no `A6`.

Exact factorization verifies that `G` is geometrically irreducible, while `T`
is irreducible by (2.5).  Their valid generic opens are connected.  Constant
embedded collision and infinity data plus the same Whitney--Thom argument
therefore propagate the two representative exclusions across both
codimension-one strata.

## 6. Exact primary decomposition of the residual

Write

\[
  C=\alpha-\beta+\gamma-1,
  \qquad
  V=\{\alpha LC\ne0\}.
\]

Exact primary decomposition over `QQ` gives three components in each half of
the residual:

\[
\sqrt{I(\operatorname{Sing}G)}=I_A\cap I_B\cap
  \langle\alpha-\gamma+2,\;\beta-2\gamma+3\rangle, \tag{6.1}
\]

\[
\sqrt{\langle G,T\rangle}=I_N\cap I_E\cap
  \langle27\alpha-8\gamma+12,\;27\beta-31\gamma+31\rangle. \tag{6.2}
\]

The last line in (6.1) lies identically in `C=0`; the last line in (6.2)
lies identically in `L=0`.  They are invalid under the hypotheses.  The four
valid primes are rational curves.  The rational parametrizations below land
on the stated primary components and have dense image: the checker substitutes
`A,B` into `G` and all three first partials, substitutes `N,E` into `(G,T)`,
and verifies a nonconstant map to each one-dimensional component.  They make
the geometry transparent; no birational-normalization claim is needed.

### 6.1 Curve `A`: contact three plus two nodes

With normalization parameter `h`, set

\[
\begin{aligned}
\alpha_A&=-\frac{2h^3(3h^2+12h+13)}{3h+11},\\
\beta_A&=-\frac{3h^5+9h^4-4h^3-24h^2+6h-4}{3h+11},\\
\gamma_A&=-\frac{6(h^2+3h-2)}{3h+11}.
\end{aligned} \tag{6.3}
\]

Then

\[
H=(s-h)^3\frac{R_h(s)}{3h+11}, \tag{6.4}
\]

where

\[
R_h=(3h+11)s^2+(3h^2+21h+34)s+6h^2+24h+26.
\]

Generically this is one contact-three collision and two nodes.  The exact
validity factors are

\[
L=-\frac{2(3h+1)^2(3h+4)^3}{3h+11},\quad
C=-\frac{3(h+1)^5}{3h+11},\quad
T=\frac{(h+3)^3(3h+1)^3}{(3h+11)^2}. \tag{6.5}
\]

### 6.2 Curve `B`: two contact-two collisions plus one node

With parameter `m`, put

\[
\begin{aligned}
p&=-\frac{2(m^2-4m+5)}{m-3},\\
q&=-\frac{2m^3-6m^2+3m+3}{m-3},\\
r&=-\frac{2(m-2)}{m-3}.
\end{aligned}
\]

There are explicit rational functions
`(alpha_B,beta_B,gamma_B)` stored in the checker for which

\[
  H=(s^2+ps+q)^2(s-r). \tag{6.6}
\]

A small valid representative is

\[
  (\alpha,\beta,\gamma)=\left(\frac16,\frac{127}{108},2\right),
  \qquad
  H=\frac{(3s+2)(6s^2+10s+3)^2}{108}. \tag{6.7}
\]

### 6.3 Curve `N`: a nonordinary triple plus one node

With parameter `u`, set

\[
\begin{aligned}
\alpha_N&=\frac{(u+1)^2(u+2)^2(2u+1)}u,\\
\beta_N&=\frac{(u^2+3u+1)(u^3+4u^2+5u+3)}u,\\
\gamma_N&=\frac{u^2+4u+1}{u}.
\end{aligned} \tag{6.8}
\]

Then

\[
H=\frac{(s+u+2)^2(su+2u+1)
  (s^2-su+u^2+2u+1)}u. \tag{6.9}
\]

The double root is one of the three pair sums of a triple collision.  Thus
the generic topology is a nonordinary triple with one tangent branch pair,
plus one separate node.  In particular, the earlier description of
`(20/3,1,2/3)` as an ordinary triple plus a separate tangency was incorrect;
that point is the `u=-3` representative of (6.9).

### 6.4 Curve `E`: an ordinary triple plus a separate tangency

With `v=beta`, set

\[
  \gamma_E=6,\qquad \alpha_E=\frac95v-9.
\]

Then

\[
H=(s+3)^2
  \frac{v+5s^3+10s^2+5s-5}{5}. \tag{6.10}
\]

The cubic factor contains the three pair sums of an ordinary triple point;
the fixed double root `s=-3` is a separate contact-two collision.  Both
`(-9,0,6)` and `(2,55/9,6)` lie in the connected generic open of this line.

## 7. The four generic residual curves are excluded

For exact rational representatives of `A`, `B`, `N`, and `E`, Sage 10.8
produces presentations on three geometric fiber meridians.  The dependency-
free checker exhausts all `40^3=64000` assignments of those meridians to
single 3-cycles.  In every case exactly 40 assignments satisfy the relators,
and every satisfying image has order three:

| residual curve | generic finite collisions | satisfying images | `A6` images |
|---|---|---:|---:|
| `A` | contact 3 + node + node | 40 copies of `C3` | 0 |
| `B` | contact 2 + contact 2 + node | 40 copies of `C3` | 0 |
| `N` | nonordinary triple + node | 40 copies of `C3` | 0 |
| `E` | ordinary triple + separate contact 2 | 40 copies of `C3` | 0 |

Each rational curve minus its finitely many invalid and exceptional values is
a connected complex curve.  On that open, the collision equivalence
relation, contact orders, forced finite cusp, and infinity pair `(5,8)` are
constant.  After a finite base change labels the collision sections,
simultaneous resolution and proper Whitney--Thom isotopy transport the exact
representative complement throughout the stratum.  This theorem dependency
is essential: no conclusion is inferred from a singularity count alone.

## 8. Every exceptional point

The factorization identities leave only three exceptional geometric types,
comprising four parameter points.

### 8.1 Contact four plus one node

At `h=-13/3` (equivalently `m=8/3`),

\[
(\alpha,\beta,\gamma)=
\left(-\frac{114244}{81},-\frac{63407}{81},\frac{34}{3}\right),
\qquad
H=\frac{(s-4)(3s+13)^4}{81}. \tag{8.1}
\]

Sage's presentation simplifier reports a noncyclic affine complement.
Regardless of that whole-group simplification, the decisive complete `40^3`
replay again finds exactly 40 order-three images and no `A6` image.

### 8.2 Contact three plus contact two

The two conjugate points are parameterized by

\[
h_\pm=-2\pm\frac{2\sqrt6}{3},\qquad
3h_\pm^2+12h_\pm+4=0, \tag{8.2}
\]

with

\[
\begin{aligned}
\alpha_\pm&=4272\mp1744\sqrt6,\\
\beta_\pm&=2384\mp\frac{8756}{9}\sqrt6,\\
\gamma_\pm&=-8\pm4\sqrt6.
\end{aligned}
\]

If `k=-3(3h+2)/2`, then

\[
H=(s-h)^3(s-k)^2. \tag{8.3}
\]

Sage was run separately on both embeddings of `QuadraticField(6)`.  Sirocco
emitted abort/retry warnings before Sage's exact Puiseux fallback completed;
the final raw presentations differ, and Sage's simplifier reports both as
noncyclic.  Independently of those whole-group observations, each exact
permutation replay has only 40 cyclic images and no `A6` image.

### 8.3 A higher nonordinary triple

At `h=-3`, `u=1`, and `v=65`, all three residual curves meet:

\[
  (\alpha,\beta,\gamma)=(108,65,6),
  \qquad H=(s+3)^3(s^2-s+4). \tag{8.4}
\]

This is a three-branch point with pairwise contacts `(3,1,1)`.  Sage's
presentation simplifier reports a cyclic complement, and its decisive
representation census again has 40 `C3` images and no `A6` image.

### 8.4 Why this list is complete

On curve `A`, the two possible changes are

\[
R_h(h)=\frac{2(h+1)^2(3h+13)}{3h+11},\qquad
\operatorname{Disc}R_h=3(h+1)^2(3h^2+12h+4).
\]

Together with (6.5), these leave exactly (8.1), (8.2), and (8.4) on `V`.
On curve `B`,

\[
\operatorname{Disc}(s^2+ps+q)=
\frac{4(m-2)(m-1)^2(3m-8)}{(m-3)^2},
\]

\[
r^2+pr+q=-\frac{(m-1)^2(2m^2-12m+15)}{(m-3)^2},
\]

which give precisely (8.1) and (8.2) after invalid factors are removed.  The
checker performs that removal as polynomial arithmetic, not by inspecting a
list of roots.  If `I_*` is the product of invalid factors and `K_*` the
product of all topology-change factors, it verifies

\[
\begin{aligned}
\operatorname{factor}(K_A/\gcd(K_A,I_A))
  &=(h+3)(3h+13)(3h^2+12h+4),\\
\operatorname{factor}(K_B/\gcd(K_B,I_B))
  &=(3m-8)(2m^2-12m+15),\\
\operatorname{factor}(K_N/\gcd(K_N,I_N))&=u-1,\\
\operatorname{factor}(K_E/\gcd(K_E,I_E))&=v-65.
\end{aligned} \tag{8.5}
\]

The factors `I_*` and `K_*` are stored explicitly in the Python certificate.
For
(6.9), the double and simple collision roots differ by

\[
  \frac{1-u^2}{u}.
\]

The root `u=-1` lies in `C=0`, while `u=1` is exactly (8.4).  For (6.10),
the cubic discriminant factors as `-25(v-5)(27v-155)`: the first root lies
in `C=0`, the second in `L=0`, and the cubic meets `s=-3` only at `v=65`,
again (8.4).  The checker also identifies (8.1) and both embeddings of (8.2)
from the `A` and `B` parameters.  Hence no valid residual parameter is
unclassified or represented only by extrapolation.

Combining Sections 3--8 proves the conditional exclusion

\[
  \boxed{\Delta=5\text{ is impossible}.} \tag{8.6}
\]

## 9. The next one-pair frontier

The genus equation at delta six leaves only `(2,17)`.  Its large affine link
has no `A6`-generating quotient with a single-3-cycle meridian.  At delta
seven the candidates are

\[
  (2,19),\qquad(3,10),\qquad(4,7).
\]

Exact `A6` censuses give respectively `0`, `720`, and `0` suitable generating
pairs.  Therefore the full result of this conditional audit is

\[
  \boxed{
  \Delta\ge7,\qquad
  \Delta=7\Longrightarrow(a,d)=(3,10),
  \quad(d-a,d)=(7,10).
  } \tag{9.1}
\]

Passing the large-link test does not construct the degree-`(3,10)` branch,
an `A6` cover, or a Keller map.  It merely names the next conditional family.

## Reproduction

Run the original generic and codimension-one checker with:

```bash
uv run python -m scripts.a6_delta_five_family
```

Run the residual algebra, all eight dependency-free presentation censuses,
and the delta-six/seven link corollary with:

```bash
uv run python -m scripts.a6_delta_five_residual
```

Regenerate the original three Sage presentations with:

```bash
sage tools/check_a6_delta_five_family.sage
```

Regenerate the residual primary decompositions, six rational presentations,
and both quadratic-number-field presentations with:

```bash
sage tools/check_a6_delta_five_residual.sage
```

For the exact primary decompositions, component memberships, rational
implicit specializations, and quadratic parameter identities without the
long Zariski--van Kamp run, use:

```bash
sage tools/check_a6_delta_five_residual.sage --algebra-only
```

The exact relators are stored and exhaustively replayed in the
dependency-free checker.  Their Zariski--van Kamp extraction remains a Sage
computation.  The family-topology step uses the standard equivalence between
plane-curve equisingularity and simultaneous resolution, followed by proper
Whitney--Thom isotopy; see Joseph Lipman's
[“Equisingularity and simultaneous resolution of singularities”](https://arxiv.org/abs/math/9802010).
