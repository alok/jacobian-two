# The `A6` delta-five family reduces to codimension two

## Claim boundary

The conditional one-pair audit proves that collision delta five is the first
remaining `A6` budget and that affine pole degrees `(3,8)` are forced at
equality.  This note attacks that family.  It proves, using exact collision
algebra, equisingular-family topology, and Sage Zariski--van Kamp
presentations, that:

- the generic three-parameter stratum has cyclic affine complement;
- the generic ordinary-triple divisor has cyclic complement; and
- the generic contact-two divisor has cyclic complement.

Consequently every admissible survivor lies in an explicit
codimension-at-least-two collision-degeneration locus.  Those deeper strata
are not fully classified here.  This is therefore a conditional,
computer-assisted family reduction—not an elimination of delta five, `A6`, or
the plane Jacobian conjecture.

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

## 3. A clean curve is locally perfect and globally cyclic

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
remaining two pairs are nodes.  Sage again gives complement group `Z`, and
the exact `40^3` replay gives only 40 cyclic images.

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

## 6. The remaining locus

Every admissible survivor is now confined to

\[
  \boxed{
  \left(\operatorname{Sing}G\ \cup\ (G\cap T)\right)
  \cap
  \{\alpha L(\alpha-\beta+\gamma-1)\ne0\}.
  }
\]

This codimension-at-least-two locus contains higher-contact collisions, two
simultaneous multiple collisions, and nonordinary combinations of a triple
collision with tangency.  The list is exhaustive under the stated
hypotheses: a triple root or two double roots of `H` lies in the singular
locus of the universal discriminant, hence in `Sing(G)` on the valid set; a
nonordinary triple collision has a tangent branch pair and lies in `G intersect
T`; and four distinct preimages cannot collide because `P` has degree three.

The residual locus is genuinely nonempty.  For example,

\[
  (\alpha,\beta,\gamma)=\left(\frac{16}{5},\frac{32}{5},\frac{24}{5}\right)
\]

has `alpha*L*C=768/5`, `T=-5`, and

\[
  H=\frac{(s+2)^3(5s^2+4s+2)}5.
\]

It is a valid point of `Sing(G) minus T`, representing a contact-three
collision and two nodes.  Likewise,

\[
  (\alpha,\beta,\gamma)=\left(\frac{20}{3},1,\frac23\right)
\]

has `alpha*L*C=250880/9`, lies in `G intersect T`, and satisfies

\[
  H=\frac{(s-1)^2(3s+5)(s^2+3s+4)}3.
\]

It represents an ordinary triple collision together with a separate
contact-two pair.  The components containing these examples have not all
been classified.  A full delta-five exclusion must finish that finite
stratification rather than extrapolate from generic points.

## Reproduction

Run the dependency-free algebra and permutation replay with:

```bash
uv run python -m scripts.a6_delta_five_family
```

Regenerate all three affine complement presentations and their explicit
isomorphisms to `Z` with:

```bash
sage tools/check_a6_delta_five_family.sage
```

The exact relators are stored in both checkers.  The family-topology step uses
the standard equivalence between plane-curve equisingularity and simultaneous
resolution; see Joseph Lipman's
[“Equisingularity and simultaneous resolution of singularities”](https://arxiv.org/abs/math/9802010).
