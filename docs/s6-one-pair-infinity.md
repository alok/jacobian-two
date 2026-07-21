# The one-pair `S6` infinity frontier starts at `(5,12)`

## Claim boundary

Assume that the irreducible transposition branch of the one-dicritical `S6`
passport has polynomial normalization

\[
  \nu(t)=(P(t),Q(t)),\qquad \deg P=m<\deg Q=n,
\]

and exactly one genuine characteristic pair at infinity.  The width bound and
the group of its affine link at infinity exclude every possible degree pair
with `n<=11`.  The first pair that passes this infinity-group test is

\[
  \boxed{(m,n)=(5,12)}.
\]

This is a conditional stratum theorem.  It does **not** prove that every
Keller branch has one characteristic pair, treat multi-pair or multi-component
infinity, construct the `(5,12)` curve or cover, or solve `JC(2)`.

## 1. The hypotheses and the two different knot pairs

For one polynomial end with `m<n`, the projective chart at infinity has
orders

\[
  (n-m,n).
\]

Thus a genuine single singular pair requires

\[
  n-m\ge2,\qquad \gcd(m,n)=1. \tag{1.1}
\]

The second condition is equivalent to `gcd(n-m,n)=1`.  It is important not
to confuse this projective pair with the link on a large sphere in affine
space.  The leading parametrization there is `(t^m,t^n)`, so the affine link
is the torus knot

\[
  T(m,n). \tag{1.2}
\]

For a connected six-sheet cover with generic transposition inertia, a generic
projection of the branch curve supplies `m` meridian generators.  Their
transposition supports form a graph on six sheets.  Transitivity requires that
graph to be connected, hence

\[
  m\ge5. \tag{1.3}
\]

Equations (1.1)--(1.3) give a finite list for each upper bound on `n`.

## 2. The group at infinity surjects onto the global complement

Let `B` be any affine plane curve and let `L_infty` be its link on a
sufficiently large sphere.  The inclusion of the link exterior induces an
epimorphism

\[
  \pi_1(S^3-L_\infty)\twoheadrightarrow
  \pi_1(\mathbb C^2-B). \tag{2.1}
\]

Leidy and Maxim give a direct Zariski--Lefschetz proof of (2.1): a generic
line can be chosen in a tubular neighborhood of the line at infinity, and
its complement already surjects onto the global affine complement.  See
Theorem 4.7 of their paper
[“Higher-order Alexander invariants of plane algebraic curves”](https://arxiv.org/abs/math/0509462).

For the one-pair end, the source of (2.1) is

\[
  G_{m,n}=\langle a,b\mid a^m=b^n\rangle. \tag{2.2}
\]

The common power `z=a^m=b^n` is central.  Suppose the global monodromy is all
of `S6`.  The image of `z` lies in the center of `S6`, so it is trivial:

\[
  A^m=B^n=1. \tag{2.3}
\]

Choose Bezout exponents `u,v` with

\[
  nu+mv=1.
\]

A meridian of `T(m,n)` is `a^u b^v`, up to conjugacy and orientation.
Consequently

\[
  M=A^uB^v \tag{2.4}

\]

must be a transposition.  Changing the Bezout solution does not change (2.4)
under (2.3).  Conjugacy or inversion also does not change its cycle type.

The obstruction is therefore finite and exact: enumerate all ordered pairs
`(A,B)` in `S6` satisfying (2.3), retain those for which (2.4) is a
transposition, and ask whether they generate **all** of `S6`.

## 3. Complete census through degree eleven

The constraints `m>=5`, `n-m>=2`, and `gcd(m,n)=1` give exactly the following
ten pairs with `n<=11`.  `Compatible` counts ordered pairs `(A,B)` satisfying
the power and meridian conditions.  The last column records the exact orders
of the generated subgroups.

| `(m,n)` | `(u,v)` | `#A` | `#B` | compatible | generated orders |
|---|---:|---:|---:|---:|---|
| `(5,7)` | `(-2,3)` | 145 | 1 | 0 | none |
| `(5,8)` | `(2,-3)` | 145 | 256 | 735 | `2:15, 120:720` |
| `(5,9)` | `(-1,2)` | 145 | 81 | 0 | none |
| `(7,9)` | `(-3,4)` | 1 | 81 | 0 | none |
| `(7,10)` | `(-2,3)` | 1 | 220 | 15 | `2:15` |
| `(5,11)` | `(1,-2)` | 145 | 1 | 0 | none |
| `(6,11)` | `(-1,2)` | 396 | 1 | 15 | `2:15` |
| `(7,11)` | `(2,-3)` | 1 | 1 | 0 | none |
| `(8,11)` | `(3,-4)` | 256 | 1 | 15 | `2:15` |
| `(9,11)` | `(-4,5)` | 81 | 1 | 0 | none |

Here `#A` and `#B` are the numbers of elements satisfying the two power
conditions.  The census is over actual ordered permutations, not conjugacy
classes.  No row has a subgroup of order `720`, so none gives an `S6`
quotient.

Most rows have short hand proofs.  If `n` is 7 or 11, then `B=1` because
`S6` has no element of order 7 or 11.  If `m=7`, then `A=1`.  The image is
therefore cyclic.  In `(5,9)`, both permitted factors are even: `A` has order
one or five and `B` has order one or three.  Their image lies in `A6`, and in
particular their meridian cannot be a transposition.

The only substantial near-miss is `(5,8)`.  Here `A` must be a 5-cycle for a
transitive image, and it fixes one sheet `f`.  Since `B^8=1`, the order of `B`
divides four.  With the displayed Bezout choice,

\[
  M=A^2B^{-3}=A^2B,
  \qquad B=A^{-2}M. \tag{3.1}
\]

If the transposition `M` avoids `f`, then both generators fix `f`; these are
the 720 order-120 rows.  If `M` contains `f`, the product of the 5-cycle
`A^{-2}` and that connecting transposition is a 6-cycle, contradicting
`B^8=1`.  This proves the exclusion without relying on the enumeration.

## 4. Sharpness at `(5,12)`

At `n=12`, the two possible pairs are `(5,12)` and `(7,12)`.  The second is
again cyclic because `A^7=1` forces `A=1`.  The first genuinely survives:

| `(m,n)` | compatible | generated orders | `S6` pairs |
|---|---:|---|---:|
| `(5,12)` | 2175 | `2:15, 120:1440, 720:720` | 720 |
| `(7,12)` | 15 | `2:15` | 0 |

One exact witness is

\[
  A=(2\ 3\ 4\ 5\ 6),
  \qquad B=(1\ 2\ 5\ 3\ 6\ 4). \tag{4.1}
\]

Then `ord(A)=5`, `ord(B)=6`, and the Bezout solution
`12*(-2)+5*5=1` gives

\[
  A^{-2}B^5=(1\ 2). \tag{4.2}
\]

This is visibly enough to generate `S6`: `A` fixes sheet one and cycles the
other five, so conjugating `(1 2)` by the powers of `A` produces every star
transposition `(1 i)`.  Those five transpositions generate `S6`.

Thus the infinity argument is sharp at degree 12.  Passing it is only a
necessary condition.  A future `(5,12)` target must still realize the forced
finite cusp/collision passport, its global affine complement quotient, and
all Keller compactification equations.

## 5. Reproduction

The checker uses only the Python standard library and the repository's small
permutation engine:

```text
uv run python -m scripts.s6_one_pair_infinity
uv run pytest -q tests/test_s6_one_pair_infinity.py
```

It enumerates all 720 elements of `S6`, filters the two power conditions,
checks the meridian cycle type, closes every compatible generated subgroup,
and independently verifies the witness and its five star transpositions.
