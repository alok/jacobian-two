# Split `C3` and `C2^2` coefficient-rank strata

Status: exact symbolic rank classification on the true `k=0,+2` split
charts.  The `k=-2` statement follows from the already checked full-family
involution.  This note does not compute complement topology or prove `JC(2)`.

## 1. Exact clean root factors

Write `u,v` for selected roots of the true split collision components.  The
root coordinates and the factors removed from the clean profile are:

| split | component | root coordinate | clean factor | reason |
|---|---|---:|---:|---|
| `k=0` | `V` | `r` | `r(2r-1)` | diagonal; component overlap |
| `k=0` | `W` | `s` | `s(s^2+1)(s^2+2)` | overlap; pair containing `t=0`; diagonal |
| `k=2` | `V` | `r` | `r(4r-1)` | overlap; diagonal |
| `k=2` | `W` | `s` | `s(s+1)(s+2)` | diagonal; overlap; diagonal |

For two roots on the same component, the clean localizer also contains
`u-v`.  These factors come from the actual unordered source pair, not from a
generic decic specialized after cancellation.

At `k=0`, the graph relation is

```text
r = (s^2+1)/2.
```

Thus `s^2+1=0` gives the pair `{0,s}`.  Any collision there contains the
forced cusp source `t=0`.  At `k=2`, the unique component intersection is
`(r,s)=(0,-1)` and is exactly the cusp-image wall.

## 2. The five allowed `C3 + 7N` systems

The allowed allocations are

```text
k=0: W, overlap-V, overlap-W
k=2: V, W.
```

For an ordinary component root, the three equations are the value and first
two root derivatives.  Their coefficient-maximal-minor gcds are

```text
k=0, W: 64*(s^2+1)^3
k=2, V: 2
k=2, W: 64.
```

The two fixed `k=0` overlap systems have ranks two and maximal-minor gcds
`4` and `8`.  Consequently every system has its expected rank on its exact
clean open.

The only raw rank loss is `k=0, W` at `s^2+1=0`.  It is not an artifact:
the coefficient and augmented ranks are both two, producing an affine
two-plane coefficient fiber over each of the two algebraic root values.
However, the pair is `{0,+i}` or `{0,-i}` and its image is the cusp image.
Moreover the graph collision polynomial has order at least four there, not
the clean order three.  The exact member

```text
(a,b,c,d) = (1,0,2,0),  s=i
```

is retained as a hostile fixture.  It explains why removing `s^2+1` is a
mathematical profile condition rather than cosmetic saturation.

## 3. Visible-only `C2^2 + 6N` rank factors

Three allowed systems have no residual clean rank divisor:

```text
k=0, VW:        det = 256*v^3*(v^2+1)
k=0, overlap+W: maximal-minor gcd = 32*v^2*(v^2+1)
k=2, VV:        det = (u-v)^4.
```

Every displayed factor is contained in the exact overlap/cusp/diagonal
localizer.  Hence all three systems have maximal coefficient rank on the
clean open.

This localization is necessary.  Exact compatible hostile members include:

* `k=0, (u,v)=(1/2,0)`, `(a,b,c,d)=(1/4,0,1,0)`: the chosen `V` and `W`
  roots are the same component-intersection pair;
* `k=0, u=v=1`, `(a,b,c,d)=(1,1,1,0)`: one ordinary graph contact has been
  labeled twice, giving coefficient and augmented rank two; and
* `k=2, u=v=1`, `(a,b,c,d)=(1,0,1,0)`: the analogous repeated vertical
  contact also has rank two.

Without `u-v`, the repeated-contact base is one-dimensional with an affine
two-plane coefficient fiber.  Its apparent dimension three is exactly the
spurious component that a valid two-distinct-contact chart must remove.

## 4. The three residual determinants

After dividing only the visible determinant factors, three irreducible
polynomials over `QQ` remain.

### 4.1 `k=0, WW`

```text
F0 = 3*u^4*v^2 + u^4 + 4*u^3*v^3 + 4*u^3*v
   + 3*u^2*v^4 + 20*u^2*v^2 + 5*u^2
   + 4*u*v^3 + 20*u*v + v^4 + 5*v^2 + 10.
```

The full determinant is

```text
512*(u-v)^4*(u^2+1)*(v^2+1)*F0.
```

### 4.2 `k=2, VW`

The determinant is `8*F2VW`, where `F2VW` is the exact irreducible
degree-eight, 29-term polynomial exported by
`scripts/a6_delta_ten_split_contact_rank.py`.  It has no visible clean factor.

### 4.3 `k=2, WW`

```text
F2WW = u^2*v + u*v^2 + 2*u^2 + 8*u*v + 2*v^2
      + 11*u + 11*v + 13,
```

and the determinant is `1024*(u-v)^4*F2WW`.

## 5. Exact localized compatibility schemes

On a residual determinant, consistency is equivalent to vanishing of the
four augmented maximal minors using the constant column.  Divide their
common factor, adjoin the residual polynomial, and call the resulting ideal
`I`.  If `h` is the clean localizer, the following reduced lexicographic
bases are certified by the ideal sandwich

```text
I subset J,       h^e J subset I,       J + (h) = (1).
```

Therefore `I : h^infinity = J`.

For `k=0, WW`, `e=2` and

```text
38*u + 4*v^3 + 71*v,
4*v^4 + 71*v^2 + 361.
```

The base is reduced of length four.  It is invariant under `u <-> v`, is
disjoint from the diagonal, and therefore gives two unordered base orbits.

For `k=2, VW`, `e=6` and

```text
414*u + 6*v^5 + 146*v^4 + 877*v^3
      + 1934*v^2 + 1527*v + 234,
6*v^6 + 74*v^5 + 367*v^4 + 932*v^3
      + 1296*v^2 + 954*v + 297.
```

The base is reduced of length six.

For `k=2, WW`, `e=1` and

```text
3*u - 46*v^5 - 336*v^4 - 1017*v^3
    - 1380*v^2 - 742*v - 117,
2*v^6 + 18*v^5 + 69*v^4 + 135*v^3
    + 134*v^2 + 60*v + 9.
```

The base is reduced of length six.  It is swap-invariant and gives three
unordered base orbits.

Squarefreeness follows exactly from gcd one between each displayed
univariate polynomial and its derivative.

## 6. Fiber dimensions

The ideals of all coefficient `3 x 3` minors become the unit ideal after
clean localization.  The exact saturation exponents for

```text
k=0 WW, k=2 VW, k=2 WW
```

are respectively

```text
2, 2, 1.
```

The determinant is zero and rank at most two is impossible, so every point
of all three compatible bases has coefficient rank exactly three.  Four
coefficients minus rank three gives an affine-line fiber.  Since every base
is zero-dimensional, every residual incidence component has dimension one:

```text
dim(base) + dim(fiber) = 0 + 1 = 1.
```

Thus no residual determinant hides a second dimension-two split incidence
surface.

The boundary fixture

```text
k=2, (u,v)=(0,-2), (a,b,c,d)=(-3,-7,-4,1)
```

is compatible of rank three but simultaneously selects the component
overlap and a diagonal pair.  It shows why raw compatibility before clean
localization contains genuine points belonging to different profiles.

## 7. Rank-open component closure

The rank calculation alone does not identify a split locus with the global
component carrying the cyclic sample.  That comparison is supplied by the
total unordered-pair surface

\[
 F=(2s+k)p-s(s^2+ks+1)=0.                           \tag{7.1}
\]

Let `C` be the coefficient of `t` in the remainder of `Q` modulo
`t^2-s*t+p`, and define the fiber-tangent derivation

\[
 D=F_p\partial_s-F_s\partial_p.                     \tag{7.2}
\]

Then `C=D(C)=0` is the global contact-two system and
`C=D(C)=D^2(C)=0` is the global contact-three system when the source roots
are distinct, the fixed-`k` pair fiber is smooth, both branches are immersed,
and the common first coordinate is unramified.  These are explicit clean-open
conditions, not consequences of the row equations alone.  On the principal
pair chart, the total rows are exactly a triangular transform of
`(H,H',H'')`, with determinant `s^6/(2s+k)^9`; two ordered contact blocks have
determinant `u^4*v^4/((2u+k)^7*(2v+k)^7)`.  Exact row transformations also
identify them with the true component jets.  Their determinants on the three maximal-rank `C3` charts
are

```text
k=0, W: s^3/512
k=2, V: -8*r^3
k=2, W: s^6*(s+1)^3/512.
```

On the five maximal-rank two-contact charts the corresponding four-row
determinants are

```text
k=0, VW:  -u^4*v*(2u-1)/128
k=0, WW:   u*v/16384
k=2, VV:   4*u*v
k=2, VW:  -u*v^4*(v+1)/64
k=2, WW:   u^4*v^4*(u+1)*(v+1)/16384.
```

Every factor is a unit on the exact clean localizer.  Power-membership
certificates prove this on the **entire** declared open; the stored clean
witness is an additional nonemptiness check.  Thus the row comparisons hold
universally there, including the pair-denominator vertical component at
`k=2`.

The total pair surface is geometrically integral and flat over the `k` line.
Its ordered self-fiber-product is integral too.  Indeed its algebra is flat
over `QQ[k]`, while its generic fiber is the tensor product of two independent
primitive-linear pair algebras and is a domain.  A zero divisor would vanish
after generic localization and hence create `QQ[k]`-torsion, contradicting
flatness.  Its only singular points are the three double-overlap points where
both ordered pairs equal the component intersection; the clean rank-open
localizers remove them.  The clean rank-open incidence is therefore the same
smooth irreducible global incidence that contains the exact cyclic sample.
The separate finite-etale labeling, simultaneous-resolution, and proper
Whitney--Thom argument in the propagation note excludes all eight rank-open
split contact charts.

Exact rational Cramer arcs now dominate both fixed `C3` overlap planes and
the `C2^2` overlap-plus-contact incidence from the ordinary global
components.  This proves algebraic containment of all three overlap
allocations.  It does not prove that complement topology is unchanged at
those boundary points.  The overlap topology and the rank-three affine-line
fibers over the three representative residual ordered bases of total length
`4+6+6=16` remain open.  Transport to `k=-2` turns the three representative
schemes into five actual contact schemes of total ordered length `28`.

## 8. Reproduction and boundary

Run the dependency-free exact certificate and its tests with

```bash
uv run python -m scripts.a6_delta_ten_split_contact_rank
uv run pytest -q tests/test_a6_delta_ten_split_contact_rank.py
uv run python -m scripts.a6_delta_ten_split_contact_closure
uv run python -m scripts.a6_delta_ten_split_c22_overlap_closure
uv run pytest -q tests/test_a6_delta_ten_split_contact_closure.py
uv run pytest -q tests/test_a6_delta_ten_split_c22_overlap_closure.py
```

Independently replay the saturations in Sage with

```bash
sage tools/check_a6_delta_ten_split_contact_rank.sage
sage tools/check_a6_delta_ten_split_c22_overlap_closure.sage
```

The exact conclusion is deliberately bounded.  Every coefficient-rank-drop
locus of the allowed split `C3` and `C2^2` allocations is either an explicit
overlap/cusp/diagonal boundary or one of the finite rank-three residual bases
above, and every clean maximal-rank locus lies in an already-excluded global
component.  All three overlap incidences are also algebraically contained,
but no complement calculation or proper-isotopy extension is yet supplied
at those boundary points or on the affine-line fibers.  No conclusion about
unrestricted `A6` or `JC(2)` is made.
