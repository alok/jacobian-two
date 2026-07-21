# A conditional one-pair obstruction for the `A6` survivor

## Claim boundary

This note proves a new conditional narrowing, not the plane Jacobian
conjecture.  Assume that the branch curve of the one-dicritical degree-six
`A6` survivor has all of the following properties:

1. its normalization is `A1`, with polynomial parametrization
   `nu(t)=(P(t),Q(t))`;
2. after an affine target change, `deg P=a<d=deg Q`, `gcd(a,d)=1`, and the
   projective closure has one polynomial end with one genuine singular
   Puiseux pair `(d-a,d)` at infinity;
3. the only intrinsic finite singularity is the forced `T(2,5)` cusp, whose
   image has no additional normalization preimage; and
4. every other finite singularity is a collision of distinct smooth
   normalization points, with total delta/contact sum `Delta>=1`.

Then the exact calculation below proves

\[
  \boxed{\Delta\ge7,}
\]

and equality can pass the large-link test only for affine degrees `(3,10)`.

The one-pair hypothesis and the restriction on finite singularities have not
been derived for every Keller compactification.  The complement presentations
are produced by Sage's Zariski--van Kamp implementation and replayed by a
dependency-free permutation checker.  The family exclusions through delta
five are therefore conditional and computer-assisted.  The symbolic
collision algebra and finite permutation replays are exact; presentation
extraction and Whitney--Thom propagation across equisingular strata are
separate computer-assisted topology dependencies.  The `Delta=3,4,6`
exclusions and the delta-five and delta-seven degree reductions use finite
exact group censuses once the same geometric hypotheses are assumed.

## 1. Genus leaves finitely many degree pairs at each delta

The degree-`d` projective closure has arithmetic genus

\[
  p_a=\frac{(d-1)(d-2)}2.
\]

Its infinity branch has local pair `(d-a,d)`, hence

\[
  \delta_\infty=\frac{(d-a-1)(d-1)}2.
\]

The finite delta is `2+Delta`: the forced `T(2,5)` cusp contributes two.
Since the normalization has genus zero,

\[
\frac{(d-1)(d-2)}2
=2+\Delta+\frac{(d-a-1)(d-1)}2,
\]

or equivalently

\[
  \boxed{(d-1)(a-1)=2(\Delta+2).} \tag{1.1}
\]

Exact divisor enumeration gives:

| collision delta | affine degrees `(a,d)` | projective pair `(d-a,d)` |
|---:|---|---|
| 1 | `(2,7)` | `(5,7)` |
| 2 | `(2,9)` | `(7,9)` |
| 2 | `(3,5)` | `(2,5)` |
| 3 | `(2,11)` | `(9,11)` |
| 4 | `(2,13)` | `(11,13)` |
| 4 | `(3,7)` | `(4,7)` |
| 5 | `(2,15)` | `(13,15)` |
| 5 | `(3,8)` | `(5,8)` |
| 6 | `(2,17)` | `(15,17)` |
| 7 | `(2,19)` | `(17,19)` |
| 7 | `(3,10)` | `(7,10)` |
| 7 | `(4,7)` | `(3,7)` |

If smooth infinity were allowed, `Delta=1` would also contain `(3,4)`.
Excluding it is exactly where the genuine singular-pair hypothesis enters.

## 2. The large affine link is not the projective local knot

Two superficially similar torus knots occur here, and interchanging them
would give a false proof.  In the `Q!=0` chart near the projective point at
infinity,

\[
  \left(\frac P Q,\frac1Q\right)
  \sim (w^{d-a},w^d),
\]

so the small projective-germ link is `T(d-a,d)`.  On a sufficiently large
affine sphere, however, the phases of the leading terms trace

\[
  \theta\longmapsto(e^{ia\theta},e^{id\theta}),
\]

so the link at infinity is `T(a,d)`.  It is the latter link complement that
surjects onto `pi_1(A2-B)` by the Zariski--Lefschetz theorem at infinity.

Write its group and meridian as

\[
  G_{a,d}=\langle x,y\mid x^a=y^d\rangle,
  \qquad m=x^uy^v,
  \qquad du+av=1. \tag{2.1}
\]

In a quotient onto centerless `A6`, the common central power dies, so
`x^a=y^d=1`.  The one-dicritical `A6` branch meridian must be a single
3-cycle.  Exhausting all relevant pairs in `A6` gives:

| affine link | meridian `(u,v)` | `A6`-generating pairs with single-3 meridian |
|---|---|---:|
| `T(2,7)` | `(1,-3)` | 0 |
| `T(2,9)` | `(1,-4)` | 0 |
| `T(3,5)` | `(-1,2)` | 720 |
| `T(2,11)` | `(1,-5)` | 0 |
| `T(2,13)` | `(1,-6)` | 0 |
| `T(3,7)` | `(1,-2)` | 0 |
| `T(2,15)` | `(1,-7)` | 0 |
| `T(3,8)` | `(-1,3)` | 720 |
| `T(2,17)` | `(1,-8)` | 0 |
| `T(2,19)` | `(1,-9)` | 0 |
| `T(3,10)` | `(1,-3)` | 720 |
| `T(4,7)` | `(-1,2)` | 0 |

Thus `Delta=1` is excluded, and one of the two `Delta=2` cases is excluded.
The coarse link test genuinely does not kill `(3,5)`.  For example,

\[
  x=(1\,6\,3)(2\,5\,4),\qquad
  y=(1\,5\,4\,3\,2)
\]

satisfy `x^3=y^5=1`, generate `A6`, and have
`x^{-1}y^2=(1 5 6)`.

The small projective `T(2,5)` group itself has no such global `A6` quotient.
That fact is not an obstruction: the forced finite cusp has local image `A5`
on five sheets, and the separate collision meridian is what enlarges global
monodromy to `A6`.

## 3. Every remaining `Delta=2` curve lies in one family

Move the forced cusp preimage and image to zero.  A degree-three first
coordinate with cusp order two can be scaled to

\[
  P=t^2+t^3.
\]

After removing the tangent term from the degree-five second coordinate, the
`T(2,5)` semigroup forbids a `t^3` term.  Scaling its leading coefficient
therefore gives the complete affine normal form

\[
  \nu_c(t)=(t^2+t^3,\;ct^4+t^5). \tag{3.1}
\]

For a collision `t!=u`, put `s=t+u` and `r=tu`.  Equality of the first
coordinates gives `r=s^2+s`; equality of the second then reduces exactly to

\[
  H_c(s)=s^2+(c+1)s+(2c-1)=0. \tag{3.2}
\]

The pair discriminant and the discriminant of `H_c` are

\[
  (t-u)^2=-s(3s+4),
  \qquad
  \operatorname{Disc}_s(H_c)=(c-1)(c-5). \tag{3.3}
\]

The complete exceptional-parameter audit is:

| `c` | event | status under the hypotheses |
|---:|---|---|
| `1/2` | the origin changes from `T(2,5)` to `T(2,7)` | invalid |
| `5/6` | an extra `T(2,3)` ramification occurs at `t=-2/3` | invalid |
| `1` | the cusp preimage `0` collides with the smooth preimage `-1` | invalid |
| `5` | one collision of two smooth branches has contact exactly two | valid exceptional fiber |

For every other `c`, (3.1) has exactly two transverse nodes.  One exact
nontransversality eliminant is

\[
 -(c-5)(c-1)(2c-1)^3(6c-5), \tag{3.4}
\]

so the table contains every degeneration relevant here.  At `c=5`, the two
parameters satisfy `t^2+3t+6=0`.  Their second derivatives `d^2Q/dP^2` reduce
modulo this quadratic to

\[
  -\frac{9t+61}{128},

\]

whose two values differ.  The contact is therefore exactly two, not higher.

## 4. Complement topology kills the whole family

For the generic representative `c=0`, the implicit equation is

\[
  -X^5+5X^3Y-5XY^2+Y^3+Y^2=0. \tag{4.1}
\]

Sage 10.8's unsimplified affine Zariski--van Kamp presentation has three
geometric fiber meridians `a,b,z`.  Its first three relators imply in order

\[
  [z,b]=1,\qquad z=b,\qquad a=b.

\]

The fourth is then redundant, so

\[
  \pi_1(\mathbb A^2-B_0)\cong\mathbb Z. \tag{4.2}

\]

Let

\[
  U=\mathbb C\setminus\{1/2,5/6,1,5\}.

\]

Over `U`, the exact exceptional-parameter calculation gives constant embedded
Puiseux and intersection data at the finite cusp and two nodes.  At infinity,
with `w=1/t`, the family is

\[
  \left(\frac P Q,\frac1Q\right)
  =\left(
    w^2\frac{1+w}{1+cw},
    w^5\frac1{1+cw}
  \right),
\]

so its embedded infinity data are constant as well.  The standard
equisingularity theorem for reduced plane-curve families therefore supplies a
Whitney stratification of
`(P2, closure(B_c) union L_infinity)` after a finite base change labels the two
node sections.  The projective family is proper, and Thom's first isotopy
lemma then makes the embedded pairs topologically locally trivial over `U`.
Since `U` is path-connected, every generic fiber has the complement group
(4.2).  This equisingularity-to-Whitney implication is a theorem dependency of
the conditional result; it is not asserted merely from singularity counts.

The remaining valid fiber `c=5` has equation

\[
 -X^5-100X^4-10X^3Y+40X^2Y+15XY^2+Y^3-4Y^2=0. \tag{4.3}
\]

Its three-meridian van Kamp presentation is stored exactly in
[`a6_one_pair_infinity.py`](../scripts/a6_one_pair_infinity.py).  The checker
exhausts all `40^3=64000` assignments of its geometric meridians to single
3-cycles.  Exactly `760` satisfy the relators: `40` have image `C3`, and `720`
have image `A5`.  None has image `A6`.

For comparison, the same replay at `c=0` finds only the `40` cyclic images.
Consequently neither the generic nor the contact-two part of the exhaustive
family supports the required global monodromy.  This excludes `Delta=2` and
completes the only family-level calculation needed below.

## 5. The coarse delta-five survivor is eliminated

Equation (1.1) gives one candidate at `Delta=3`, namely `(2,11)`, and two at
`Delta=4`, namely `(2,13)` and `(3,7)`.  Their large affine links have no
`A6`-generating quotient at all: in the first two cases the second torus
generator must be trivial because `A6` has no element of order `11` or `13`;
the exact census also kills `T(3,7)`.

At `Delta=5`, equation (1.1) gives `(2,15)` and `(3,8)`.  The first link has
`2880` generating pairs in `A6`, but their meridians are split evenly between
cycle types `(4,2)` and `(5,1)`; none is the required single 3-cycle.  The
second link has `4320` generating pairs, of which exactly `720` have a
single-3-cycle meridian.  One witness is

\[
  x=(1\,2\,3)(4\,5\,6),\qquad
  y=(1\,2)(3\,4\,6\,5),\qquad
  x^{-1}y^3=(2\,3\,4).
\]

The large-link calculation alone therefore gives the intermediate statement

\[
  \boxed{
    \Delta\ge5,
    \quad
    \Delta=5\Longrightarrow(a,d)=(3,8)
    \text{ after the large-link test}.
  }
\]

The [delta-five family audit](a6-delta-five-family.md) then treats every
degree-`(3,8)` curve in the normalized three-parameter family.  It excludes
the generic open, both codimension-one walls, all four generic residual
curves, and every valid exceptional point.  The last eight presentation
censuses each exhaust `40^3` single-3-cycle assignments and find exactly 40
cyclic images, never `A6`.  Thus `Delta=5` is impossible under the same
hypotheses.

At `Delta=6`, equation (1.1) leaves only `(2,17)`, and its large affine link
has zero suitable `A6` images.  At `Delta=7`, the three candidates are
`(2,19)`, `(3,10)`, and `(4,7)`; their exact suitable-image counts are `0`,
`720`, and `0`.  The strongest conclusion of the combined audit is therefore

\[
  \boxed{
    \Delta\ge7,
    \qquad
    \Delta=7\Longrightarrow(a,d)=(3,10),
    \quad(d-a,d)=(7,10).
  }
\]

The [generic delta-seven audit](a6-delta-seven-generic.md) gives the complete
four-parameter normal form and excludes its nondegenerate open by an exact
cyclic-complement representative plus proper Whitney--Thom transport.  A
conditional equality survivor must therefore lie on the repeated-collision
or triple-image walls `G=0` or `T=0`.  No `A6` cover or Keller map is
constructed.  Those walls, multi-pair infinity, and branches violating the
finite-singularity hypotheses remain outside this theorem.

## Reproduction and source boundary

Run the dependency-free replay with:

```bash
uv run python -m scripts.a6_one_pair_infinity
uv run python -m scripts.a6_delta_five_family
uv run python -m scripts.a6_delta_five_residual
uv run python -m scripts.a6_delta_seven_generic
```

Together these verify the genus candidates, symbolic family and residual
identities, all relevant torus-quotient censuses through delta seven, and all
stored presentation censuses.  Extracting the presentations themselves uses
Sage 10.8/Sirocco; reproduce that extraction with
`sage tools/check_a6_one_pair_infinity.sage` and
`sage tools/check_a6_delta_five_residual.sage`, and reproduce the generic
delta-seven presentation with
`sage tools/check_a6_delta_seven_generic.sage`.  The large-link epimorphism is
the standard Zariski--Lefschetz theorem at infinity.  The family steps use
proper Whitney--Thom isotopy, not mere constancy of singularity counts.  These
topological dependencies are distinct from the exact symbolic and finite
permutation replays.

Useful source:

- Eva Elduque and Laurentiu Maxim,
  [“Higher-order degrees of affine plane curve complements”](https://people.math.wisc.edu/~lmaxim/hoa2.pdf),
  Section 2, for the link-at-infinity epimorphism onto the affine curve
  complement group.
- Joseph Lipman,
  [“Equisingularity and simultaneous resolution of singularities”](https://arxiv.org/abs/math/9802010),
  for the equivalence between plane-curve equisingularity and simultaneous
  resolution used in the generic-family step.
