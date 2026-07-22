# Conditional `A6` delta-ten codimension-three checkpoint

## Claim boundary

This note starts the fourteen expected codimension-three profiles in the
conditional one-dicritical `A6` delta-ten audit.  It treats the valid
**nonsplit** unordered-pair charts of two selected profiles,

\[
 C_4+6N
 \qquad\text{and}\qquad
 C_2+C_3+5N.                                                       \tag{0.1}
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

An exact rational member of each dominant surface has the prescribed
contacts and residual nodes, cyclic affine complement, and no required `A6`
three-cycle quotient.  On precisely defined nonempty clean opens,
finite-etale labeling, relative contact blowups, and proper Whitney--Thom
isotopy propagate that complement topology.  Thus both dominant clean
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

## 10. Remaining boundary

The result currently removes the dominant clean nonsplit surfaces for
`C4+6N` and `C2+C3+5N`.  It also bounds the `C4` residual determinant
incidence below surface dimension.  It does not classify:

1. `C4+6N` on the true split pair charts `k=0,+2,-2`;
2. the pair-denominator, cusp-pair, or diagonal charts;
3. intersections of those loci with `F=0` or other collision walls;
4. specializations where residual nodes merge or meet a contact or cusp
   target;
5. non-clean equisingular specializations on the `C2+C3` compatibility
   surface; or
6. any of the other twelve codimension-three profiles.

## 11. Reproduction

Run the exact Python certificate, tests, types, and lint with:

```bash
uv run python -m scripts.a6_delta_ten_contact_four
uv run python -m scripts.a6_delta_ten_contact_two_three
uv run pytest -q \
  tests/test_a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_two_three.py
uv run mypy --no-incremental \
  scripts/a6_delta_ten_contact_four.py \
  scripts/a6_delta_ten_contact_two_three.py \
  tests/test_a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_two_three.py
uv run ruff check \
  scripts/a6_delta_ten_contact_four.py \
  scripts/a6_delta_ten_contact_two_three.py \
  tests/test_a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_two_three.py
```

Replay the determinantal saturations, primitive implicit curve, singular
scheme, raw presentation, and cyclic simplification with Sage 10.8:

```bash
sage tools/check_a6_delta_ten_contact_four.sage
sage tools/check_a6_delta_ten_contact_two_three.sage
```

The Sage checker is manual; GitHub CI runs the Python certificate and test
suite.
