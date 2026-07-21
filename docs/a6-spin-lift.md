# The forced `A6` local relations lift to `2.A6`

## Status and claim boundary

This note gives an exact **hostile consistency certificate**, not a Keller map
and not an elimination of the one-dicritical `A6` passport.

The forced finite data consist of a `T(2,5)` cusp and a separate collision of
two local degree-three branches.  In the natural six-sheet action, choose

```text
r = (345),   s = (123),   b = (456).
```

The exact result is:

> The canonical odd-order lifts of `r`, `s`, and `b` to the spin double cover
> `2.A6` satisfy the cusp five-braid relation and the collision commutator
> relation.  They generate all of `2.A6`.  Moreover, the same forced prefix
> has product-one completions by 3-cycles with both spin signs.

Consequently, the presently forced finite monodromy does not yield a spin
obstruction.  A successful spin argument would have to derive the actual
infinity word and its peripheral framing from a Keller compactification.

## The double cover and exact Clifford model

The natural permutation representation sends `A6` into `SO(6)`.  Pulling back
`Spin(6) -> SO(6)` gives the nonsplit central extension

```text
1 -> {+1,-1} -> 2.A6 -> A6 -> 1.
```

Here `2.A6` is isomorphic to `SL(2,9)` and has order `720`.  The full Schur
multiplier of `A6` has order six, so the universal Schur cover is `6.A6`, not
the double cover used here.  The distinction matters.

Work in the Euclidean Clifford algebra with `e_i^2=1`, and set

```text
alpha_ij = (e_i-e_j)/sqrt(2).
```

For a 3-cycle `(ijk)=(ij)(jk)`, its canonical order-three lift is

```text
widehat(ijk) = alpha_ij alpha_jk.
```

The checker avoids irrational coefficients: it multiplies the two root
differences first and divides by two, so every coefficient is rational and
all identities are exact.

Put

```text
R = alpha_34 alpha_45,
S = alpha_12 alpha_23,
B = alpha_45 alpha_56.
```

All three have order three.  The other preimage of any one of them has order
six.

## Cusp and collision relations

The `T(2,5)` meridians obey

```text
rsrsr = srsrs.
```

Exact Clifford multiplication gives the stronger equality

```text
RSRSR = SRSRS = (-e12-e24+e15-e45)/2,
(RSRSR)^2 = -1.
```

Thus there is no hidden central sign in the five-braid relation itself.

At the separate `3+3` collision, `s=(123)` and `b=(456)` commute.  Their
canonical lifts commute as well:

```text
SB = BS.
```

The lifts do not define a splitting of `A6`.  Indeed,

```text
RB = alpha_34 alpha_56,
(RB)^2 = -1.
```

Since `r,s,b` generate `A6`, the lifts generate a group projecting onto all
`360` elements of `A6` and containing the central kernel.  Its exact
enumerated order is therefore `720`.

## The cusp longitude retains a framing sign

The central element of the `T(2,5)` knot group is `(rs)^5`, and the canonical
lifts satisfy

```text
(RS)^5 = -1.
```

For the standard preferred longitude

```text
lambda = (rs)^5 r^(-10),
```

the lifted value is `-R^(-1)`, while the permutation value is `r^(-1)`.
This is a genuine local framing sign.  It is not by itself an obstruction:
the knot-group kernel may map to the central element.  One would need a
global infinity theorem forcing the opposite lift of the peripheral word.

## An even-order infinity completion

For the simplest ordered prefix `p=rsb`,

```text
p = (1243)(56).
```

If `Q=RSB`, exact multiplication gives `Q^4=-1`, so `Q` has order eight.  The
element `I=Q^(-1)` lifts `p^(-1)` and satisfies

```text
RSBI = 1.
```

This is a direct product-one lift.  It is not a proposed Keller infinity
word.  Moreover, an even-order permutation has no distinguished odd-order
lift, so changing `I` to `-I` changes the total central sign.  Finite local
data alone cannot choose between them.

## Both Fried--Serre spin signs occur with 3-cycles only

The ambiguity persists even if every branch cycle has odd order.  Let

```text
c = (356) = r b^(-1) r^(-1).
```

The following two tuples both multiply to one and generate `A6`:

```text
g_plus  = (r,s,b,s^(-1),r^(-1),c),
g_minus = (r,s,b,s^(-1),r,b).
```

Using the unique order-three lift of each entry gives

```text
spin(g_plus)  = +1,
spin(g_minus) = -1.
```

For the first tuple, conjugation makes the lifted product telescope to one.
For the second it reduces to `(RB)^2=-1`.  Thus the forced prefix is
compatible with either component of the relevant 3-cycle Nielsen space.

## Reproducible certificate

Run

```bash
uv run python -m scripts.a6_spin_lift
uv run pytest -q tests/test_a6_spin_lift.py
```

[`a6_spin_lift.py`](../scripts/a6_spin_lift.py) implements the rational
Clifford algebra, enumerates the generated groups, checks the cusp and
collision relations, verifies the longitude and infinity identities, and
checks both product-one tuples.  The adversarial test flips the collision
lift to its order-six preimage and verifies that the canonical certificate
then fails.

## What this closes and what remains

This calculation closes the simplest spin-only attack:

- the five-braid relation lifts without a central error;
- the disjoint collision meridians commute upstairs;
- the cusp longitude sign can be absorbed at infinity; and
- even all-odd-order completions occur with both spin signs.

It does **not** show that either completion is realized by a polynomial
Keller map.  The remaining useful target is a compactification theorem that
determines the infinity peripheral word and its canonical lift, or couples
that framing to the trace lattice and boundary intersection form.

## Sources

- Michael D. Fried, [*Alternating Groups and Moduli Space Lifting
  Invariants*](https://www.math.uci.edu/~mfried/paplist-cov/hf-can0611591.pdf),
  especially Section 1.1.4 for the spin pullback, unique odd-order lifts, and
  the lifting invariant.
- [ATLAS of Finite Group Representations: `A6`](https://brauer.maths.qmul.ac.uk/Atlas/v3/alt/A6/),
  for the group order, multiplier, and `2.A6 = SL(2,9)` identification.
- Jean-Pierre Serre,
  [*L'invariant de Witt de la forme `Tr(x^2)`*](https://eudml.org/doc/139995),
  for the relation between the trace form and the spin-cover class.
