# Conditional `A6` delta-ten codimension-three checkpoint

## Claim boundary

This note starts the fourteen expected codimension-three profiles in the
conditional one-dicritical `A6` delta-ten audit.  It treats the valid
**nonsplit** unordered-pair/fiber charts of four selected profiles,

\[
 C_4+6N
 \qquad
 C_2+C_3+5N
 \qquad
 C_2+Q_0+2N
 \qquad\text{and}\qquad
 C_3+T_{111}+4N.                                                   \tag{0.1}
\]

Here `Cm` is a two-branch contact of intersection multiplicity `m`.

For `C4+6N`, the exact incidence calculation proves that the chart has one
rational Cramer surface and that its residual determinant curve hides only a
one-dimensional incidence, not a second surface.  For `C2+C3+5N`, the five
jet equations have a 409-term compatibility residual in `(k,u,v)`.  Sage
derives and factors it over `QQ`; one smooth rational point proves that the
irreducible residual is geometrically irreducible and has a nonempty
rank-four clean open.  Exact saturation finds a degree-thirty coefficient
rank-drop curve but makes the augmented rank-at-most-three ideal the unit
ideal, so that entire curve is inconsistent and hides no incidence.  The
full-localizer singular-Jacobian saturation is also the unit ideal, proving
that the compatibility surface is smooth on the valid chart.

For `C2+Q0+2N`, the augmented determinant separates a same-target boundary
from one genuine irreducible compatibility quartic.  Its projective curve has
arithmetic genus three and exactly two ordinary nodes, both on the removed
split fibers, so the valid base is a smooth genus-one open times the fiber
parameter.  Full-localizer saturation makes every coefficient and augmented
rank-at-most-three ideal the unit ideal.  A Cramer chart and two independent
coefficient-image tangent vectors prove that its image is a codimension-three
surface rather than a collapsed curve.

For `C3+T111+4N`, three contact jets and two ordinary-triple remainder
equations give a five-by-five augmented determinant.  After removing the
explicit singular-fiber and contact/triple overlap factors, one irreducible
62-term degree-nine compatibility polynomial remains.  Its full-valid
singular saturation is the unit ideal.  The coefficient-rank-at-most-three
ideal is a degree-fourteen curve, but the augmented-rank-at-most-three ideal
is the unit ideal, so every point of that curve is inconsistent.  Two exact
coefficient-image tangent vectors prove that the rank-four incidence maps to
a genuine codimension-three surface.

An exact rational member of each of the four dominant surfaces has the prescribed
contacts and residual nodes, cyclic affine complement, and no required `A6`
three-cycle quotient.  On precisely defined nonempty clean opens,
finite-etale labeling, relative contact blowups, and proper Whitney--Thom
isotopy propagate that complement topology.  Thus all four dominant clean
nonsplit surfaces are excluded.  Manual Sage 10.8 regenerates the incidence
factorizations, implicit curves, singular schemes, and van Kamp
presentations; it does not certify the propagation theorem.

The true split fibers `k=0,+2,-2`, the pair-denominator, cusp-pair, and
diagonal charts, and their intersections are not classified here.  This note
is conditional on the same four branch hypotheses as the preceding
[codimension-two checkpoint](a6-delta-ten-codim-two.md), and it does not prove
`JC(2)`.

## 1. The fourteen profiles

The exact contact-tree ledger gives the following expected codimension-three
rows:

1. `C2 + C3 + 5N`;
2. `C2 + Q0 + 2N`;
3. `C2 + T111^2 + 2N`;
4. `C2 + T112 + 4N`;
5. `C2^2 + T111 + 3N`;
6. `C2^3 + 4N`;
7. `C3 + T111 + 4N`;
8. `C4 + 6N`;
9. `Q2 + 3N`;
10. `T111 + Q0 + N`;
11. `T111 + T112 + 3N`;
12. `T111^3 + N`;
13. `T113 + 5N`; and
14. `T222 + 4N`.

Expected codimension is only a work-ordering device.  It is not used to
declare any incidence empty or to bound a coefficient fiber after its matrix
loses rank.

## 2. The contact-four incidence

In the normalized family

\[
 P=t^2+kt^3+t^4,
 \qquad
 Q=at^5+bt^6+ct^7+dt^8+t^9,                                      \tag{2.1}
\]

let `H(s)` be the collision decic on the nonsplit unordered-pair chart.  An
exact contact four is a quadruple root, so impose

\[
 H(s)=H'(s)=H''(s)=H'''(s)=0.                                     \tag{2.2}
\]

These four equations are affine-linear in `(a,b,c,d)`.  Their coefficient
determinant is

\[
 24(k+2s)^4(s^2+ks+1)^5F(k,s),                                   \tag{2.3}
\]

where

\[
\begin{aligned}
F={}&3k^5s+35k^4s^2+13k^4+88k^3s^3+164k^3s\\
 &+88k^2s^4+412k^2s^2+52k^2+32ks^5\\
 &+416ks^3+64ks+160s^4.
\end{aligned}                                                     \tag{2.4}
\]

Exact factorization proves that `F` is irreducible over `QQ`.  Work on

\[
 D\!\left(k(k^2-4)s(k+2s)(s^2+ks+1)(2s^2+3ks+4)\right).           \tag{2.5}
\]

This localization explicitly removes the three split fibers, the
pair-denominator, cusp-pair, and diagonal loci.  On `F!=0`, Cramer's rule
identifies the incidence with an open of the irreducible `(k,s)` plane, so it
is one rational surface.  Because `H` is monic in `s`, projection to
coefficient space is quasi-finite on the valid incidence and generically
finite onto a codimension-three image.

## 3. The residual determinant cannot hide another surface

The factor `F=0` needs an actual rank audit.  A dimension expectation alone
would miss a large coefficient fiber.  The independent Sage checker forms
the complete coefficient and augmented determinantal ideals over `QQ[k,s]`
and saturates by (2.5).  It proves:

| exact valid-chart ideal | result |
|---|---|
| coefficient rank at most two | unit ideal, saturation exponent `2` |
| compatibility on `F=0` | dimension `0`, scheme length `10`, exponent `4` |
| augmented rank at most two | unit ideal, exponent `2` |
| augmented rank at most one | unit ideal, exponent `1` |

Thus every compatible valid point of `F=0` has coefficient and augmented
rank exactly three.  There are only finitely many such base points, and each
has an affine-line coefficient fiber.  Their total incidence dimension is

\[
 0+(4-3)=1.                                                        \tag{3.1}
\]

Consequently the residual curve contains neither another contact-four
surface nor any larger incidence component.  This conclusion is stronger
than showing one augmented minor is merely generically nonzero on `F`.

## 4. An exact clean member

Take

\[
 (k,a,b,c,d)=
 \left(1,\frac97,\frac{27}{14},3,\frac{12}{7}\right).              \tag{4.1}
\]

The incidence determinant is `-168`, the valid localizer is `-9`, and

\[
 H(s)=\frac{(s+1)^4}{14}
 \left(14s^6+28s^5+54s^4+52s^3+56s^2+21s+18\right).               \tag{4.2}
\]

The residual sextic has nonzero discriminant

\[
 -\frac{2627098335909}{2582630848},                                \tag{4.3}
\]

is disjoint from the contact root, and supplies six residual nodes.  The
contact sources are the roots of `t^2+t+1`.  After scaling `Q` by fourteen,
their common image is `(0,-1)`.  The resultant of the contact quadratic with
`P'` is `3`, so both branches are immersed and `X=P(t)` is a local parameter.
Their first three `d/dX` jets agree, while the fourth-jet difference is

\[
 984(2t+1),                                                        \tag{4.4}
\]

whose resultant with `t^2+t+1` is `3`.  Hence their local intersection
multiplicity is exactly four.

The primitive implicit curve has affine Jacobian components

```text
((4,1), (7,1), (6,6)).
```

These are the length-four forced `T(2,5)` cusp, the length-seven
contact-four singularity, and six reduced nodes.  Together with the fixed
`T(5,9)` infinity branch, the projective delta balance is

\[
 2+4+6+16=28=p_a(\text{degree }9).                                \tag{4.5}
\]

## 5. Complement calculation

Sage regenerates a raw affine Zariski--van Kamp presentation with four
geometric meridians and ten relations.  Its exact simplification is
infinite cyclic.  The dependency-free replay checks every one of the
`40^4=2,560,000` assignments of those meridians to three-cycles in `A6`:

```text
satisfying assignments: 40
generated image orders: {3: 40}
A6 assignments: 0
```

Thus the exact member cannot carry the required connected six-sheet `A6`
quotient.  The extraction of the presentation remains a computer-assisted
Sage dependency; the finite permutation replay is exact once the stored
relation words are present.

## 6. Propagation across the clean Cramer surface

On the determinant-nonzero chart, Cramer's rule identifies the incidence
with

\[
 D\!\left(k(k^2-4)s(k+2s)(s^2+ks+1)
          (2s^2+3ks+4)F(k,s)\right)\subset\mathbb A^2_{k,s}.        \tag{6.1}
\]

This base is smooth and irreducible.  Shrink it by the exact nonvanishing
conditions that define the clean profile:

- `a!=0` for the forced cusp, together with nonzero cusp-image and
  extra-critical resultants;
- nonzero contact-pair discriminant and `Res(pair,P')`;
- nonzero fourth-jet difference and `H''''(s)`;
- nonzero residual-sextic leading coefficient and discriminant;
- nonzero residual boundary and tangency resultants;
- squarefree residual node-target eliminant; and
- separation of node, contact, and cusp target sections.

These are principal open conditions.  The sample (4.1) lies in their
intersection: among the recorded values are `a=9/7`, cusp-image factor
`1/196`, extra-critical factor `11421/49`, pair discriminant `-3`,
`Res(pair,P')=3`, fourth-jet norm `3`, `H''''(-1)=492/7`, nonzero residual
discriminant (4.3), and nonzero node-target and special-target separations.
Hence the clean incidence is a nonempty irreducible complex surface and is
analytically path-connected.

After a finite-etale ordering cover, the two contact branches and all six
residual nodes are global labeled sections.  Resolve each node once.  The
contact-four section requires **four** relative blowups: the first three
successively lower the contact order to one, but after the third the two
branches and the newest exceptional divisor form a triple crossing; the
fourth separates those three directions and gives relative simple normal
crossings.  The cusp and infinity sections retain their fixed `(2,5)` and
`(5,9)` embedded types and have fixed resolution sequences.

The resulting projective family is smooth and proper over the labeling
cover, and the reduced total transform of the curve together with
`L_infinity` is a relative SNC divisor.  The same simultaneous-resolution
and proper-isotopy argument proved in
[the propagation note](a6-delta-ten-propagation.md) applies.  Thom's first
isotopy lemma makes the resolved complements locally trivial; blowing up is
an isomorphism off the total transform, so the original affine complements
are homeomorphic.  The labeling cover is surjective and its relevant
component is path-connected.  Therefore

\[
 \pi_1(\mathbb C^2\setminus B_b)\cong\mathbb Z                     \tag{6.2}
\]

throughout this dense clean Cramer open.  The dominant nonsplit `C4+6N`
surface cannot realize the required `A6` passport.

This transport does not include the finite residual-compatibility scheme on
`F=0` or any boundary omitted by (6.1).

## 7. Contact two plus contact three

Choose distinct unordered-pair sums `u` and `v`, with the contact two at
`u` and the contact three at `v`.  The incidence equations are

\[
 H(u)=H'(u)=H(v)=H'(v)=H''(v)=0.                                 \tag{7.1}
\]

Their `5 x 5` augmented determinant factors exactly as

\[
 -2(u-v)^6(u^2+ku+1)(v^2+kv+1)^3R(k,u,v).                       \tag{7.2}
\]

The checked Sage source derives `R` from (7.1); it is not copied into the
repository as an uncheckable 409-term literal.  Its exact shape is

```text
terms:        409
total degree: 18
degrees:      deg_k=13, deg_u=9, deg_v=10
QQ factors:   one, with exponent one
```

The valid open also removes `k(k^2-4)uv(u-v)`, both pair denominators,
cusp-pair factors and diagonal factors, the resultant that makes the two
source-pair quadratics intersect, and `k+u+v`, which makes their target
`P`-values coincide.  No one of these factors is being declared impossible;
each names a separate boundary chart.

At

\[
 (k,u,v)=\left(\frac23,1,-1\right),                              \tag{7.3}
\]

the residual vanishes and, with the normalization in (7.2), its gradient is

\[
 \left(
 \frac{1073741824}{59049},
 \frac{536870912}{59049},
 \frac{1073741824}{177147}
 \right).                                                        \tag{7.4}
\]

The coefficient and augmented matrices both have rank four there, and the
valid localizer is nonzero.  Because `R` is irreducible over `QQ`, a smooth
`QQ`-point also proves geometric irreducibility: a rational point on several
distinct Galois-conjugate geometric components would lie in their
intersection and would be singular.  Therefore the compatible graph over a
nonempty open of `R=0` is one irreducible surface.  Monicity of `H` makes the
map to coefficient space quasi-finite after retaining the labeled contact
roots, so its image closure is a codimension-three component.

Expected dimension is not used to dispose of the locus where all `4 x 4`
coefficient minors vanish.  After stripping only declared localizer units,
that coefficient-rank-at-most-three locus is a genuine curve: its projective
Hilbert polynomial is `30*t-51`.  Saturation exponent `4` records the exact
coefficient calculation; saturation by the product of every remaining valid
factor makes no further change and returns exponent `0`.  But the normalized augmented-rank-at-most-three
ideal saturates to the unit ideal, with exponent `3`.  Hence every point of
the degree-thirty rank-drop curve is inconsistent.  Redundant coefficient
and augmented rank-at-most-two saturations are also unit ideals, both with
exponent `3`.  Every valid compatible base therefore has coefficient and
augmented rank exactly four and a unique coefficient solution; there is no
second surface or lower-rank affine fiber on this chart.

The raw singular locus of `R=0` has dimension one, so it must not be ignored.
After the invertible change `x=2u+k`, `y=2v+k`, exact modular-reconstruction
saturation of `(R,R_k,R_u,R_v)` by the complete valid localizer gives the
unit ideal with exponent `2`.  Thus every singular point lies on one of the
explicitly removed split, overlap, denominator, cusp-pair, diagonal, or
same-target boundaries; the valid compatibility surface itself is smooth.

## 8. An exact `C2+C3+5N` member

The smooth point (7.3) has the unique coefficient solution

\[
 (a,b,c,d)=(3,1,3,0).                                            \tag{8.1}
\]

Thus

\[
 P=t^2+\frac23t^3+t^4,
 \qquad
 Q=3t^5+t^6+3t^7+t^9.                                           \tag{8.2}
\]

The collision and tangency polynomials factor as

\[
 H(s)=\frac{(s-1)^2(s+1)^3}{9}
 \left(9s^5+27s^4+45s^3+53s^2+46s+12\right),                   \tag{8.3}
\]

and

\[
 T(s)=\frac{4(s-1)(s+1)^2}{9}
 \left(
 81s^8+270s^7+477s^6+516s^5+313s^4
 +24s^3-59s^2-66s-20
 \right).                                                        \tag{8.4}
\]

The residual quintic has discriminant `256979755008` and values `192,-8`
at `s=1,-1`.  Its resultants with the pair denominator, cusp-pair factor,
diagonal factor, and tangency polynomial are respectively

```text
-1024/27, 2048/9, 28672, -42392524029926834176.
```

The two source-pair quadratics and their images after replacing `X` by `3P`
are

| contact | pair quadratic | image `(X,Y)` | `Res(pair,(3P)')` |
|---|---|---:|---:|
| `C2` at `s=1` | `t^2-t+1` | `(-5,3)` | `252` |
| `C3` at `s=-1` | `t^2+t+1` | `(-1,-1)` | `36` |

Both pair discriminants are `-3`; their mutual source resultant and cleared
target-separation numerator are `4096/81` and `-4096/243`.  Exact `d/dX`
jets agree through order one at the `C2` point and through order two at the
`C3` point.  The first differing jets have nonzero norms `27/2401` and
`3/16`, proving that the contact orders are exactly two and three.

With `X=3P` and `Y=Q`, Sage regenerates the primitive irreducible nonic and
its affine Jacobian primary pieces:

```text
((5,5), (5,1), (3,1), (4,1)).
```

These are five reduced nodes, the `C3` point, the `C2` point, and the forced
`T(2,5)` cusp.  The Jacobian algebra has length `17` and radical length `8`.
The projective delta balance is

\[
 2+2+3+5+16=28=p_a(\text{degree }9).                            \tag{8.5}
\]

Sage also regenerates the exact four-generator, ten-relation affine van Kamp
presentation and simplifies it to `Z`.  The independent finite replay gives

```text
assignments:          2,560,000
satisfying:           40
generated orders:     {3: 40}
A6 assignments:       0
```

## 9. Propagation on the dominant compatibility surface

Intersect the smooth valid compatibility surface `R=0` with the principal clean
conditions certified by the sample: nonzero cusp and extra-critical factors,
distinct immersed source pairs and targets, exact contact jets, squarefree
residual collision and node-target polynomials, and separation from the cusp
and both contact images.  This is a nonempty open of a geometrically
irreducible complex surface, hence it is analytically path-connected.

On a finite-etale labeling cover, the four contact branches and five nodes
become disjoint global sections.  Resolve each node once, the contact-two
section with two relative blowups, and the contact-three section with three.
Together with the fixed cusp and infinity resolution sequences, this gives a
relative SNC divisor in a smooth proper projective family.  The proper
Whitney--Thom argument from the propagation note then transports the affine
complement.  Consequently the dense clean rank-four `C2+C3+5N` surface has
cyclic complement and cannot realize the required `A6` passport.

This transport makes no assertion about singular/equisingular or other
removed boundary strata.

## 10. Contact two plus an ordinary quadruple

Let `h` be the common `P`-value of four distinct source points and let
`u` be the sum of a separate unordered source pair.  Write `r1,r2,r3`
for the nonconstant coefficients of the remainder of `Q` modulo `P-h`.
The five equations

\[
 r_1=r_2=r_3=H(u)=H'(u)=0                                      \tag{10.1}
\]

are affine-linear in `(a,b,c,d)`.  Their augmented determinant is

\[
 h^2(u^2+ku+1)A(k,u)G(k,h,u)^2,                                \tag{10.2}
\]

where

\[
\begin{aligned}
A={}&k^3u+5k^2u^2+3k^2+4ku^3+16ku+12u^2+4,\\
G={}&h(k+2u)^2+u(u+k)(u^2+ku+1)^2.
\end{aligned}                                                   \tag{10.3}
\]

The factor `G=0` says exactly that the contact pair has `P`-value `h`;
it is the profile-changing boundary where the pair lands on the quadruple
target.  The genuine compatibility component is `A=0`.

The full valid localizer used in the checker is

\[
\begin{aligned}
L={}&k(k-2)(k+2)hu(k+2u)(u^2+ku+1)(2u^2+3ku+4)\\
 &\quad\cdot
 (256h^2+27hk^4-144hk^2+128h-4k^2+16)G.                        \tag{10.4}
\end{aligned}
\]

The penultimate factor is the reduced discriminant of `P-h`; the complete
quartic discriminant is `-h` times it.  Thus (10.4) removes the true split
fibers, zero or singular quadruple fibers, invalid or diagonal contact pairs,
and the same-target overlap.  These removed loci are named boundaries, not
declared empty.

Exact factorization makes `A` irreducible over `QQ`.  Its projective
quartic has arithmetic genus three and exactly two ordinary nodes,
`(-2:1:1)` and `(2:-1:1)`; both lie over the removed split fibers.
There are no singularities at infinity.  Consequently the valid
compatibility base is a smooth open of a geometrically irreducible genus-one
curve times the open `h`-line.  The rational point

\[
 (k,h,u)=(-4,1,1)                                                \tag{10.5}
\]

has `(A_k,A_h,A_u)=(4,0,8)` and `L=-45563904`.  It both proves that the
valid open is nonempty and supplies the smooth rational point needed to pass
from `QQ`-irreducibility to geometric irreducibility.

Expected rank is not used to dismiss hidden coefficient fibers.  Exact
modular-reconstruction saturation by the complete `L` gives:

| valid-chart ideal on `A=0` | result | exponent |
|---|---:|---:|
| singular ideal `(A,A_k,A_u)` | unit | `1` |
| coefficient rank at most three | unit | `1` |
| augmented rank at most three | unit | `1` |
| coefficient rank at most two | unit before saturation | `0` |
| augmented rank at most two | unit before saturation | `0` |

Every valid compatible base therefore has coefficient and augmented rank
exactly four and a unique coefficient vector.

On the Cramer chart `(k+2u)(ku+2)!=0`, four rows solve explicitly for
`(a,b,c,d)`; the fifth residual is

\[
 \frac{(u^2+ku+1)A(k,u)G(k,h,u)}{(k+2u)(ku+2)}.                 \tag{10.6}
\]

The chart-only factor `ku+2` is not inverted in the global saturation.
At (10.5), the `h` tangent and the quartic tangent
`(dk,du)=(2,-1)` map to

\[
 (0,7,-18,8,-1),\qquad(2,-4,16,-6,2)                           \tag{10.7}
\]

in `(k,a,b,c,d)`.  They are independent, so the coefficient image has
dimension two and codimension three.  The compatibility incidence has not
collapsed to a curve.

## 11. An exact `C2+Q0+2N` member

The unique coefficients over (10.5) are

\[
 (a,b,c,d)=(7,-19,13,-6),                                      \tag{11.1}
\]

so

\[
 P=t^4-4t^3+t^2,qquad
 Q=t^9-6t^8+13t^7-19t^6+7t^5.                                 \tag{11.2}
\]

The quadruple fiber is `P=1`.  Its discriminant is `-4944`, and

\[
 Q+1=(P-1)(t^5-2t^4+4t^3-t^2-1).                              \tag{11.3}
\]

At a root of `P-1`, the last factor is the branch slope `dQ/dP`.  Its
four values have eliminant

\[
 m^4-538m^3-3098m^2-5898m-3727
\]

with discriminant `-626582784`.  Thus the four immersed branches have
distinct tangents and meet in one ordinary quadruple point at `(1,-1)`.

The collision polynomial factors as

\[
\begin{aligned}
H(s)={}&(s-1)^2(s^2-10s+7)\\
 &\cdot(s^6-12s^5+50s^4-80s^3+37s^2-20s+16).                 \tag{11.4}
\end{aligned}
\]

The sextic is exactly the six off-diagonal pair sums in the quadruple fiber.
It is squarefree, has discriminant `56316985344`, and every one of its
roots maps to `X=1`.  The contact pair is `t^2-t+1`, is disjoint from the
quartic fiber with resultant `4`, and maps to `(3,-7)`.
Its two branches have equal first graph derivative `dQ/dP=-4`; their
second-derivative difference is `(2t-1)/98`, with norm `3/9604`.
Hence the contact is exactly `C2`.

The residual quadratic in (11.4) has discriminant `72` and supplies two
nodes.  Their `X`-coordinates satisfy

\[
 X^2+564X-476=0,                                                 \tag{11.5}
\]

whose discriminant is `320000` and whose values at the cusp, quadruple,
and contact targets `X=0,1,3` are `-476,89,1225`.  Thus all special
targets are distinct.  The forced-cusp and extra-critical resultants are
`-2` and `-311472`.

Sage regenerates the primitive irreducible nonic and the affine Jacobian
primary lengths

```text
((2,2), (9,1), (3,1), (4,1)).
```

They are respectively the two reduced nodes, the ordinary quadruple point,
the `C2` point, and the forced `T(2,5)` cusp.  The Jacobian algebra has
length `18` and radical length `5`.  The projective delta balance is

\[
 2+2+6+2+16=28=p_a(\text{degree }9).                            \tag{11.6}
\]

The exact affine van Kamp presentation has four generators and nine
relations and simplifies to `Z`.  The independent finite replay gives

```text
assignments:          2,560,000
satisfying:           40
generated orders:     {3: 40}
A6 assignments:       0
```

## 12. Propagation on the compatibility surface

Intersect the smooth geometrically irreducible valid surface from Section 10
with the principal clean conditions checked by the sample: nonzero cusp and
extra-critical factors, a reduced four-source fiber with four distinct
slopes, an immersed exact contact pair at a separate target, two squarefree
residual nodes, and separation of every labeled source and target section.
The result is a nonempty irreducible complex surface and therefore
analytically path-connected.

After a finite-etale ordering cover, the four quadruple branches, two contact
branches, and two nodes become labeled sections.  Resolve each node once.
One relative blowup resolves the ordinary quadruple section because its four
tangent directions are distinct; two relative blowups resolve the contact
two section.  Together with the fixed cusp and infinity resolution
sequences, the reduced total transform is a relative SNC divisor in a smooth
proper projective family.  Proper Whitney--Thom isotopy transports the
cyclic affine complement from (11.2) throughout this dense clean open.
Consequently the dominant clean nonsplit `C2+Q0+2N` surface cannot carry
the required connected six-sheet `A6` quotient.

This transport does not include the split fibers, singular quadruple fibers,
invalid contact charts, same-target overlap `G=0`, non-clean
specializations, or deeper intersections.

## 13. Contact three plus an ordinary triple

Choose the fourth, omitted root `e` of a four-point `P`-fiber.  The cubic
containing the other three roots is

\[
 C_e(t)=\frac{P(t)-P(e)}{t-e}
 =t^3+(e+k)t^2+(e^2+ke+1)t+e(e^2+ke+1).             \tag{13.1}
\]

The linear and quadratic coefficients of `Q mod C_e` vanish exactly when
those three roots have one common `Q`-value.  If `w` is the unordered-pair
sum of a separate contact-three event, the full incidence equations are

\[
 H(w)=H'(w)=H''(w)=0,
 \qquad
 [t](Q\bmod C_e)=[t^2](Q\bmod C_e)=0.               \tag{13.2}
\]

They are affine-linear in `(a,b,c,d)`.  With the right-hand-side convention
used by the certificate, the augmented determinant is

\[
 2(e^2+ke+1)^2(w^2+kw+1)^3O(k,e,w)^3R(k,e,w),        \tag{13.3}
\]

where

\[
\begin{aligned}
O={}&e^2k+2e^2w+ek^2+3ekw+2ew^2\\
   &+k^2w+2kw^2+k+w^3+w.                             \tag{13.4}
\end{aligned}
\]

The residual `R` has 62 terms, total degree nine, and variable degrees
`(6,3,7)` in `(k,e,w)`.  Sage factors it as one exponent-one factor over
`QQ`.  The separate factor

\[
 B=e^2k+2e^2w-ekw-2ew^2+kw^2+w^3+w                 \tag{13.5}
\]

is also a required removed boundary.  Indeed, the exact source and target
identities are

\[
 \operatorname{Res}_t(C_e,S_w)=B O^2,
 \qquad
 P(e)D_w-N_w=B O,                                   \tag{13.6}
\]

where `S_w` is the cleared contact-pair quadratic and `N_w/D_w` its target
`X`-coordinate.  Thus `B=0` or `O=0` makes the prescribed contact and triple
events overlap in source or target data; neither is silently counted as a
clean separate-event profile.

The complete base localizer used in the rank and smoothness audit is

\[
\begin{aligned}
L={}&k(k-2)(k+2)e w(k+2w)(w^2+kw+1)(2w^2+3kw+4)\\
 &\cdot(e^2+ke+1)(4e^2+3ke+2)\,T_D\,B\,O,             \tag{13.7}\\
T_D={}&16e^4+8e^3k-5e^2k^2+16e^2+3ek^3-12ek-k^2+4.
\end{aligned}
\]

Here

\[
 \operatorname{Disc}(C_e)=-(e^2+ke+1)T_D,
 \qquad P'(e)=e(4e^2+3ke+2),                          \tag{13.8}
\]

so every factor in (13.7) has an explicit fiber, pair, or overlap meaning.

Exact full-localizer saturation gives:

1. `(R,R_k,R_e,R_w):L^infinity=(1)`, with saturation exponent one;
2. the coefficient-rank-at-most-three ideal is a curve with Hilbert
   polynomial `14t-21`, also at exponent one;
3. the augmented-rank-at-most-three saturation is `(1)`, at exponent one;
4. both coefficient and augmented rank-at-most-two saturations are `(1)`,
   again at exponent one.

Therefore the valid compatibility surface is smooth, and every compatible
base has coefficient and augmented rank exactly four.  The degree-fourteen
coefficient-rank-drop curve is entirely inconsistent rather than a hidden
positive-dimensional coefficient fiber.

At

\[
 (k,e,w)=(-4,-1/2,1),                                 \tag{13.9}
\]

the compatibility gradient is `(-84,-112,0)`, the localizer is `7223580`,
and the unique coefficients are

\[
 (a,b,c,d)=\left(\frac{39}{2},-\frac{409}{8},
                  \frac{109}{4},-\frac{31}{4}\right). \tag{13.10}
\]

The two compatibility-surface tangent directions in base coordinates

\[
 (4,-3,0),\qquad(0,0,1)                               \tag{13.11}
\]

lift to coefficient derivatives

\[
\begin{aligned}
&\left(\frac{4147}{14},-\frac{20431}{28},337,-\frac{571}{14}\right),\\
&\left(\frac{7293}{28},-\frac{4686}{7},297,-\frac{1023}{28}\right).
\end{aligned}                                          \tag{13.12}
\]

Their images in `(k,a,b,c,d)` are independent.  Hence the coefficient image
has tangent dimension two at this point and is genuinely a codimension-three
surface, not a collapsed curve.

## 14. An exact `C3+T111+4N` member

Use target coordinates `X=P` and `Y=8Q`.  The sample is

\[
\begin{aligned}
X(t)&=t^4-4t^3+t^2,\\
Y(t)&=8t^9-62t^8+218t^7-409t^6+156t^5.                \tag{14.1}
\end{aligned}
\]

The ordinary-triple cubic is

\[
 C_e=t^3-\frac92t^2+\frac{13}{4}t-\frac{13}{8},
 \qquad \operatorname{Disc}(C_e)=-\frac{637}{4}.       \tag{14.2}
\]

Its three sources map to

\[
 (X,Y)=\left(\frac{13}{16},-\frac{2197}{128}\right).   \tag{14.3}
\]

The omitted fourth source has `Q`-value differing by `63/128`, so it does
not join the triple target.  For the slope calculation we use the unscaled
second coordinate `Q`; multiplication by eight does not affect whether the
slopes are distinct.  Writing

\[
 Q-Q_0=C_eW,
\]

gives

\[
 W=t^6-\frac{13}{4}t^5+\frac{75}{8}t^4+\frac{13}{4}t^3
   -\frac{13}{8}t^2-\frac{169}{64}t-\frac{169}{128}.   \tag{14.4}
\]

The tangent slope at a root `r` of `C_e` is **not** merely `W(r)`; because
`X-X_0=C_e(t)(t-e)`, it is `W(r)/(r-e)`.  Eliminating `r` from the corrected
slope equation gives the primitive cubic

\[
 1179648m^3-643069952m^2-6362836480m-15599304175,      \tag{14.5}
\]

whose discriminant is

\[
 -89062908728555597524369408000000\ne0.                \tag{14.6}
\]

Thus the three branches are smooth with distinct tangent directions: this
is an ordinary triple point.

The contact pair is `t^2-t+1`.  Both sources map to `(X,Y)=(3,-199)`, its
`P'` resultant is `84`, and the first two branch-graph jet differences
vanish.  The third is

\[
 \frac{165}{5488}t-\frac{165}{10976},
\]

with norm `81675/120472576`, so the contact is exactly `C3`.

The collision and tangency polynomials factor as

\[
\begin{aligned}
H={}&\frac12(s-1)^3
 (2s^3-18s^2+47s-26)
 (s^4-12s^3+24s^2+7s+24),\\
T={}&2(s-1)^2
 (18s^9-432s^8+4033s^7-18938s^6+47963s^5\\
 &\qquad-64986s^4+47578s^3-31966s^2+21980s-6240).
                                                               \tag{14.7}
\end{aligned}
\]

Their gcd is `(s-1)^2`.  The cubic pair-sum factor is the three pairs inside
the ordinary triple; its `X`-eliminant is `36(16X-13)^3`.  The residual
quartic gives four nodes, with primitive `X`-polynomial

\[
 16X^4+9324X^3-60343X^2+117047X-60306,                \tag{14.8}
\]

whose discriminant is `-4321892177659568`.  It is nonzero at the cusp,
contact, and triple targets.

The primitive implicit equation has 18 terms and is irreducible over `QQ`.
Sage regenerates its affine Jacobian algebra with length `17` and radical
length `7`.  Its primary-component `(length, radical length)` pairs are

\[
 (4,4),\quad(4,1),\quad(5,1),\quad(4,1),               \tag{14.9}
\]

corresponding respectively to four nodes, the ordinary triple, the `C3`
point, and the forced cusp.  The projective delta balance is

\[
 4+3+3+2+16=28=p_a(9),                                \tag{14.10}
\]

including the fixed `T(5,9)` branch at infinity.

Sage's raw affine van Kamp presentation has four generators and ten
relations.  Its checked simplification has one generator and no relations,
and every original meridian maps to that generator.  Thus

\[
 \pi_1(\mathbb C^2\setminus C)\cong\mathbb Z.          \tag{14.11}
\]

The dependency-free replay of the stored raw words reports

```text
assignments:          2560000
satisfying:           40
generated orders:     {3: 40}
A6 assignments:       0
```

## 15. Propagation on the new compatibility surface

Intersect the smooth geometrically irreducible valid surface from Section 13
with the open conditions checked by the sample: squarefree triple cubic,
three distinct triple slopes, an exact immersed contact three, four reduced
nodes, nonzero cusp and extra-critical factors, and separation of every
labeled source and target section.  This is a nonempty irreducible complex
surface, hence analytically path-connected.

After a finite-etale labeling cover, resolve each node once and the ordinary
triple section once.  Three relative blowups resolve the contact-three
section: after two the strict branches are transverse, while the third
removes the resulting triple crossing with the newest exceptional divisor.
Together with the fixed cusp and infinity sequences, the reduced total
transform is a relative SNC divisor in a smooth proper projective family.
Proper Whitney--Thom isotopy therefore transports the cyclic affine
complement from (14.1) through this dense clean open.  It cannot carry the
required connected six-sheet `A6` quotient.

This transport excludes neither the split fibers nor the removed
singular-fiber, overlap/same-target, non-clean, and deeper boundary charts.

## 16. Remaining boundary

The result currently removes the dominant clean nonsplit surfaces for
`C4+6N`, `C2+C3+5N`, `C2+Q0+2N`, and `C3+T111+4N`.  It also bounds the `C4`
residual determinant incidence below surface dimension.  It does not
classify:

1. `C4+6N` on the true split pair charts `k=0,+2,-2`;
2. the pair-denominator, cusp-pair, or diagonal charts;
3. intersections of those loci with `F=0` or other collision walls;
4. specializations where residual nodes merge or meet a contact or cusp
   target;
5. non-clean equisingular specializations on the `C2+C3` or `C2+Q0`
   compatibility surfaces;
6. the singular-fiber and same-target `C2+Q0` boundaries;
7. the singular-fiber, overlap/same-target, and non-clean
   `C3+T111` boundaries; or
8. any of the other ten codimension-three profiles.

## 17. Reproduction

Run the exact Python certificate, tests, types, and lint with:

```bash
uv run python -m scripts.a6_delta_ten_contact_four
uv run python -m scripts.a6_delta_ten_contact_two_three
uv run python -m scripts.a6_delta_ten_contact_quadruple
uv run python -m scripts.a6_delta_ten_contact_three_triple
uv run pytest -q \
  tests/test_a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_two_three.py \
  tests/test_a6_delta_ten_contact_quadruple.py \
  tests/test_a6_delta_ten_contact_three_triple.py
uv run mypy --no-incremental \
  scripts/a6_delta_ten_contact_four.py \
  scripts/a6_delta_ten_contact_two_three.py \
  scripts/a6_delta_ten_contact_quadruple.py \
  scripts/a6_delta_ten_contact_three_triple.py \
  tests/test_a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_two_three.py \
  tests/test_a6_delta_ten_contact_quadruple.py \
  tests/test_a6_delta_ten_contact_three_triple.py
uv run ruff check \
  scripts/a6_delta_ten_contact_four.py \
  scripts/a6_delta_ten_contact_two_three.py \
  scripts/a6_delta_ten_contact_quadruple.py \
  scripts/a6_delta_ten_contact_three_triple.py \
  tests/test_a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_two_three.py \
  tests/test_a6_delta_ten_contact_quadruple.py \
  tests/test_a6_delta_ten_contact_three_triple.py
```

Replay the determinantal saturations, primitive implicit curve, singular
scheme, raw presentation, and cyclic simplification with Sage 10.8:

```bash
sage tools/check_a6_delta_ten_contact_four.sage
sage tools/check_a6_delta_ten_contact_two_three.sage
sage tools/check_a6_delta_ten_contact_quadruple.sage
sage tools/check_a6_delta_ten_contact_three_triple.sage </dev/null
```

The Sage checkers are manual; GitHub CI runs the Python certificates and test
suite.
