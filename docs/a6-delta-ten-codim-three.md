# Conditional `A6` delta-ten codimension-three checkpoint

## Claim boundary

This note starts the fourteen expected codimension-three profiles in the
conditional one-dicritical `A6` delta-ten audit.  It treats the valid
**nonsplit** unordered-pair chart of the first selected profile,

\[
 C_4+6N,                                                            \tag{0.1}
\]

where `C4` is a two-branch contact of intersection multiplicity four.

The exact incidence calculation proves that this chart has one rational
surface and that its residual determinant curve hides only a
one-dimensional incidence, not a second surface.  An exact rational member
has the prescribed contact, six nodes, cyclic affine complement, and no
required `A6` three-cycle quotient.  On a precisely defined nonempty clean
open of the Cramer surface, finite-etale labeling, four relative contact
blowups, and proper Whitney--Thom isotopy propagate that complement topology.
Thus the dominant clean nonsplit surface is excluded.  Manual Sage 10.8
regenerates the determinantal saturations, implicit curve, singular scheme,
and van Kamp presentation; it does not certify the propagation theorem.

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

## 7. Remaining boundary

The result currently removes the dominant nonsplit contact-four surface and
its residual determinant curve from the algebraic incidence audit.  It does
not classify:

1. `C4+6N` on the true split pair charts `k=0,+2,-2`;
2. the pair-denominator, cusp-pair, or diagonal charts;
3. intersections of those loci with `F=0` or other collision walls;
4. specializations where the six residual nodes merge or meet the contact or
   cusp targets; or
5. any of the other thirteen codimension-three profiles.

## 8. Reproduction

Run the exact Python certificate, tests, types, and lint with:

```bash
uv run python -m scripts.a6_delta_ten_contact_four
uv run pytest -q tests/test_a6_delta_ten_contact_four.py
uv run mypy --no-incremental \
  scripts/a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_four.py
uv run ruff check \
  scripts/a6_delta_ten_contact_four.py \
  tests/test_a6_delta_ten_contact_four.py
```

Replay the determinantal saturations, primitive implicit curve, singular
scheme, raw presentation, and cyclic simplification with Sage 10.8:

```bash
sage tools/check_a6_delta_ten_contact_four.sage
```

The Sage checker is manual; GitHub CI runs the Python certificate and test
suite.
