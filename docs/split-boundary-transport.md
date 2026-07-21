# Split boundary determinant transport

## Status and claim boundary

Orevkov's determinant identity for one selected target boundary component is
scalar when that component has exactly one noncontracted irreducible
preimage.  That hypothesis is his condition `(*)` in the two-characteristic-
pair argument.

This note derives the exact matrix replacement when there are several
noncontracted preimages.  For positive dual square, Hodge index supplies a
sharp negative-semidefinite inequality.  The main finding is adverse but
decisive for proof strategy:

> The matrix identity and the Hodge inequality do not force condition `(*)`
> in degree six.  Explicit unimodular boundary trees obtained by three and
> four blowups support split data with generic inertia `3+3` and `2+2+2`,
> respectively.

These are exact intersection-lattice stopping models.  They are not
holomorphic maps, full compactifications, or Keller counterexamples.  They
show that a proof must add canonical-discrepancy, local-map, connectivity, or
global splice information; determinant transport by itself cannot replace
`(*)`.

## 1. The matrix transport identity

Let


```text
F : X_tilde -> X
```

be a proper generically finite map of smooth projective surfaces of degree
`N`.  Let `D` be a compact SNC divisor on `X` with nondegenerate intersection
matrix, put `D_tilde=F^{-1}(D)_red`, and assume the intersection matrix of
`D_tilde` is nondegenerate.  Let the chosen component of `D` be `T` and
suppose the noncontracted components of `D_tilde` over `T` are

```text
C_1,...,C_r.
```

Write


```text
F^*T = n_1 C_1 + ... + n_r C_r + contracted components,
F_*C_i = m_i T.
```

Thus `n_i` is the normal multiplicity and `m_i` is the tangential degree.
Counting a generic point of `T` gives


```text
n.m = sum_i n_i*m_i = N.                         (1.1)
```

Let `T^vee` be the divisor class supported on `D` that is
intersection-dual to `T`.  For every `C_i`, let `C_i^vee` be its corresponding
dual class in `D_tilde`.  Set


```text
q      = (T^vee)^2,
Q_ij   = C_i^vee . C_j^vee.
```

Equivalently, `q` is one entry of the inverse target intersection matrix and
`Q` is the selected block of the inverse source intersection matrix.
Orevkov defines determinants using the **negative** intersection matrix.  If
`delta_T=det(L-T)/det(L)` in that convention, then

```text
delta_T = -q.
```

Projection formula gives two useful class identities.  Testing against every
component of `D_tilde` gives


```text
F^*T^vee = sum_i m_i C_i^vee.                    (1.2)
```

Testing the pushforward of a source dual class against every component of `D`
component gives


```text
F_*C_i^vee = n_i T^vee.                          (1.3)
```

Intersecting (1.2) with `C_i^vee` and using (1.3) yields


```text
Q m = q n.                                       (1.4)
```

Multiplying by `m^t`, or simply squaring (1.2), gives


```text
m^t Q m = N q.                                   (1.5)
```

When `r=1`, equation (1.4) becomes


```text
Q_11 = (n_1/m_1) q,
```

which is exactly the determinant-ratio factor in Orevkov's Proposition
1.2(b).  With `r>1`, it is a coupled linear system.  There is no valid scalar
cancellation component by component.

The use of `D_tilde`, rather than every component of a larger compactification
boundary, is essential.  A dicritical mapping onto an affine branch curve is
not a component of `F^{-1}(D)`; inserting it into the inverse intersection
matrix would add dual-basis terms and invalidate (1.2) unless its pushforward
intersections were included separately.

## 2. The sharp Hodge inequality when `q>0`

Assume `q>0`.  The class in (1.2) has square `Nq>0`.  By Hodge index, its
orthogonal complement is negative definite on numerical divisor classes.
For any column vector `x`, subtract its projection onto `m` with respect to
`Q`.  Equations (1.1) and (1.4) give


```text
x^t Q x <= (q/N)*(n.x)^2.                        (2.1)
```

Equivalently,


```text
Q - (q/N) n n^t  is negative semidefinite.       (2.2)
```

The matrix in (2.2) kills `m`, so the estimate is sharp.  This is the natural
split-preimage replacement for the scalar determinant relation.  The sign
assumption is essential: for `q<0`, Hodge index supplies no such upper bound.
In particular, Orevkov's first two-pair star has
`delta_T=R_1*D_1>0`, hence `q=-R_1*D_1<0`; inequality (2.2) does **not** apply
at that star.  There the matrix projection identity (1.4) remains valid, but
one cannot attach a Hodge upper bound by silently reversing determinant signs.

The question is whether integrality and degree six make (1.1), (1.4), and
(2.2) impossible when `r>1`.  The following blowup trees answer no.

## 3. A three-blowup `A6` stopping model

Start from a boundary line of square `+1` and perform:


1. a blowup at a smooth point of component `0`;
2. a blowup at the intersection of components `0` and `1`;
3. a blowup at a smooth point of the new component `2`.

The resulting tree has


```text
weights = (-1,-2,-2,-1),
edges   = (0--2, 1--2, 2--3).
```

Its raw intersection matrix and inverse are


```text
A = [[-1, 0, 1, 0],
     [ 0,-2, 1, 0],
     [ 1, 1,-2, 1],
     [ 0, 0, 1,-1]],

A^-1 = [[1,1,2,2],
        [1,0,1,1],
        [2,1,2,2],
        [2,1,2,1]].
```

In particular, `det(A)=-1`: this is an actual unimodular boundary lattice
arising from blowups, not an arbitrary rational matrix.

Select source components `0` and `3`.  Their dual block is


```text
Q = [[1,2],
     [2,1]].
```

For a target line with `q=1`, take


```text
m = (1,1),
n = (3,3).
```

Then, exactly,


```text
n.m     = 6,
Qm      = n,
m^tQm   = 6,
Q-nn^t/6 = [[-1/2, 1/2],
            [ 1/2,-1/2]].
```

The last matrix has eigenvalues `0,-1`.  Its generic normal partition is
`(3,3)`, an even permutation type compatible with `A6`.

## 4. A four-blowup `S6` stopping model

Again start from the `+1` line, and perform:


1. two smooth blowups on component `0`;
2. a smooth blowup on component `1`;
3. a blowup at the intersection of components `1` and `3`.

The result is


```text
weights = (-1,-3,-1,-2,-1),
edges   = (0--1, 0--2, 1--4, 3--4),
```

with inverse intersection matrix


```text
A^-1 = [[1,1,1, 1, 2],
        [1,0,1, 0, 0],
        [1,1,0, 1, 2],
        [1,0,1,-1,-1],
        [2,0,2,-1,-2]].
```

Here `det(A)=1`.  Selecting components `4` and `2` gives


```text
Q = [[-2,2],
     [ 2,0]].
```

Take


```text
m = (1,2),
n = (2,2).
```

Then


```text
n.m     = 6,
Qm      = n,
m^tQm   = 6,
Q-nn^t/6 = [[-8/3, 4/3],
            [ 4/3,-2/3]].
```

This remainder has determinant zero and negative trace.  The normal indices,
repeated according to tangential degree, give `(2,2,2)`, the three-
transposition type available in `S6`.

## 5. A stronger labeled `A6` local fixture

The preceding examples intentionally retain only intersection and inertia
data.  A larger exact model shows that adding the local augmented-canonical
pullback rule still does not force uniqueness at one chosen component.

On the target, perform the blowups

```text
free(0), corner(0,1), corner(0,2),
free(3), free(4), free(5), free(6), free(7).
```

The final weights and augmented-canonical labels are

```text
weights = (-2,-2,-2,-2,-2,-2,-2,-2,-1),
labels  = (-2,-1,-3,-5,-4,-3,-2,-1, 0).
```

Choose target component `T=7`.  It has valency two, label `-1`, and raw dual
square `q=2`.

On the source, perform

```text
free(0), corner(0,1), free(2),
corner(2,3), free(1), corner(1,5).
```

This gives

```text
weights = (-1,-4,-3,-2,-1,-2,-1),
labels  = (-2,-1,-3,-2,-5, 0,-1).
```

Choose the disjoint valency-two components `(4,6)`.  Their inverse block and
degree vectors are

```text
Q = [[6, 4],
     [4,-2]],
m = (1,1),
n = (5,1).
```

Then

```text
Qm=2n,   n.m=6,   m^tQm=12=6q,
Q-(q/6)nn^t = (-7/3)*[[ 1,-1],
                      [-1, 1]].
```

The selected source labels are exactly

```text
(-5,-1) = (5,1)*label(T),
```

which is the type-one logarithmic canonical pullback rule.  At a valency-two
boundary chart, the monomial germs `(u,v)=(x^5,y)` and `(u,v)=(x,y)` realize
the two normal indices and tangential degree one.  This is a local
compatibility statement, not a gluing of all charts into a surface map.

Finally, the generic inertia can be `(12345)` with sheet six fixed.  Together
with `(456)` it generates a transitive group of order `360`, namely `A6`.
Thus Hodge transport, blowup integrality, local canonical labels, valency-two
monomial behavior, and primitive `A6` monodromy are simultaneously
consistent at the selected component.

The model still assigns no images or pullbacks to the remaining target and
source components.  It therefore neither violates nor proves condition `(*)`
for a complete compactification; it isolates the genuinely global gluing
that is missing.  In particular, the source tree here is only the vertical
divisor `D_tilde`.  No label-three dicritical leaf is inserted into its inverse
matrix: such a leaf lies outside `F^{-1}(D)` and must be coupled through the
full compactification by separate equations.

## 6. Exact certificate

Run


```bash
uv run python -m scripts.split_boundary_transport
uv run pytest -q tests/test_split_boundary_transport.py
```

[`split_boundary_transport.py`](../scripts/split_boundary_transport.py)
replays the blowups rather than merely storing the final graphs.  It computes
both full intersection matrices and inverses, extracts the selected dual
blocks, checks `Qm=qn`, degree, pullback square, Hodge semidefiniteness, and
generic inertia parity.  An adversarial test perturbs a normal degree and
checks that transport and degree fail.  It also reconstructs the labeled
fixture, verifies both canonical-label pullbacks and valencies, and enumerates
the generated permutation group.

## 7. What the models do and do not show

They prove a precise limitation:

- projection formula plus Hodge index does not force a unique
  noncontracted preimage;
- unimodularity of a boundary obtained by blowups does not fix the failure;
- degree six and the relevant `A6`/`S6` inertia parity do not fix it.
- even the local type-one canonical-label rule and primitive `A6` generation
  do not force uniqueness at a selected valency-two component.

They do **not** prove that either decorated tree is realized by the resolved
map of a Keller pair.  In particular, they do not yet include:

- the dicritical and its forced minimal edge `3--2` or `2--1`;
- the canonical and Jacobian discrepancy labels on every component;
- local holomorphic maps at the contracted connectors;
- the full target branch splice diagram; or
- the affine finite algebra and its global monodromy.

Those are now the only plausible sources of a split-preimage obstruction.
The next useful lemma must couple at least one of them to (1.4), rather than
reusing the determinant ratio in scalar form.

## Primary source

- S. Yu. Orevkov,
  [*Counter-examples to the Jacobian Conjecture at
  Infinity*](https://www.math.univ-toulouse.fr/~orevkov/jci-e.pdf),
  Proposition 1.2 for the unique-preimage determinant transform and Section
  2.4 for condition `(*)` and the two-characteristic-pair splice equations.

The matrix generalization, Hodge consequence, and the two explicit degree-six
blowup fixtures are derived and checked here.
