# The refined six-sheet defect budget

## Claim boundary

Let

\[
F=(P,Q):\mathbb A^2_{\mathbb C}\longrightarrow\mathbb A^2_{\mathbb C},
\qquad \det JF\in\mathbb C^\times,
\]

have topological degree six.  This note proves the following necessary
condition for a hypothetical counterexample:

\[
\boxed{G=6T15=A_6\quad\text{or}\quad G=6T16=S_6.}
\tag{0.1}
\]

Here `G` is the monodromy group in its transitive action on the six sheets.
The conclusion is unconditional: it does **not** assume that there is only one
dicritical component.  It strengthens the earlier primitive-group frontier
`A5`, `S5`, `A6`, or `S6`.

The argument has two new ingredients.

1. Riemann--Hurwitz refines Orevkov's exact degree-six identity to

   \[
   \boxed{\sum_E(e_Ed_E+\delta_E)=5,\qquad \delta_E\ge0.}
   \tag{0.2}
   \]

   Thus tangential degree is not free.
2. If the finite normalization has one irreducible branch curve and that
   curve's normalization is injective, Lin--Zaidenberg topology contradicts
   connected global monodromy.

The finite-group enumeration is exact and executable.  The geometric bridge
is proved below.  No surviving abstract profile is asserted to come from a
Keller map, and (0.1) is not a proof of the plane Jacobian conjecture.

## 1. Orevkov's terms and tangential degree

Resolve the rational extension of `F` at infinity as in Orevkov.  For a
dicritical component `E`, write

\[
e_E=\mu_E f^*
\]

for its generic normal local multiplicity.  Orevkov's Lemma 2.1 gives

\[
E^\circ=E\setminus L_\infty\cong\mathbb A^1.
\]

Let `B_E` be the affine image curve and factor the nonconstant restriction
through its normalization:

\[
E^\circ\mathop{\longrightarrow}^{r_E}B_E^\nu
\mathop{\longrightarrow}^{\nu_E}B_E.
\]

The map `r_E` is finite.  Its extension between projective completions shows
that the completion of `B_E^nu` is rational.  Every point missing from
`B_E^nu` must lift to the unique point missing from `E^circ`; since the target
normalization is affine, exactly one point is missing.  Hence

\[
B_E^\nu\cong\mathbb A^1.
\]

After choosing coordinates, `r_E` is therefore a polynomial.  Denote its
degree by `d_E`; this is the tangential degree of the dicritical.

For `c in E^circ`, let `a_c` be the local degree of `r_E` at `c`, and let
`x=pi(c)` be the corresponding point of Orevkov's contracted source.  Local
degree multiplicativity gives

\[
\mu_x f^*\ge e_Ea_c. \tag{1.1}
\]

One can see (1.1) by taking a nearby target point away from the branch curve.
The `a_c` nearby tangential points which coalesce at `c` each contribute
`e_E` normal sheets.  This argument remains valid at the point produced by
contracting a constant boundary chain because Orevkov's `mu_x` is precisely
the local topological multiplicity of the contracted map.

## 2. The exact Riemann--Hurwitz refinement

For a characteristic-zero polynomial of degree `d_E`, finite
Riemann--Hurwitz says

\[
\sum_{c\in E^\circ}(a_c-1)=d_E-1. \tag{2.1}
\]

Only finitely many summands are nonzero.  Define

\[
\delta_E=
\sum_{c\in E^\circ}
\bigl(\mu_{\pi(c)}f^*-e_Ea_c\bigr). \tag{2.2}
\]

Equation (1.1) makes `delta_E` a nonnegative integer.  The bracket belonging
to `E` in Orevkov's Lemma 4.2 becomes

\[
\begin{aligned}
C_E
&=e_E+
  \sum_{c\in E^\circ}(\mu_{\pi(c)}f^*-e_E)\\
&=e_E+
  e_E\sum_c(a_c-1)+
  \sum_c(\mu_{\pi(c)}f^*-e_Ea_c)\\
&=e_Ed_E+\delta_E.
\end{aligned}
\tag{2.3}
\]

Orevkov's degree-six identity is `sum_E C_E=6-1`.  Substitution proves

\[
\boxed{
\sum_E(e_Ed_E+\delta_E)=5,
\qquad
\sum_Ee_Ed_E\le5.
}
\tag{2.4}
\]

This does not change the meaning of Orevkov's original generic term: that
term is still `e_E`, not `e_E d_E`.  The extra
`e_E(d_E-1)` is forced inside his finite jump sum.  For example, one
`(e,d)=(2,2)` dicritical starts with outer term two but necessarily spends two
additional jump units, so its minimum total cost is four.

An `e=1` dicritical is generically unramified and contributes no nontrivial
cycle to a generic branch meridian.  It still costs `d_E+delta_E` in (2.4).

## 3. Moved sheets are defect cost

Let `B_i` be an irreducible component of the reduced branch locus of the
finite normalization, and let `sigma_i` be a generic meridian permutation.
Every ramified dicritical `E` dominating `B_i` supplies `d_E` disjoint cycles
of length `e_E`.  Consequently

\[
c_i:=\#\operatorname{supp}(\sigma_i)
=\sum_{\substack{E\to B_i\\e_E>1}}e_Ed_E. \tag{3.1}
\]

Summing (3.1) over the branch components and using (2.4) yields the support
budget

\[
\boxed{\sum_i c_i\le5.} \tag{3.2}
\]

Every `sigma_i` fixes an affine inverse sheet, so `c_i<=5`, and the conjugacy
classes of the `sigma_i` normally generate the transitive monodromy group.
A five-cycle is nevertheless unavailable: it requires a dicritical with
`e=5=N-1`, which Orevkov excludes in the final remark of his paper.

The exact fixed-sheet class data for the four primitive groups are:

| group | cycle type | class size | normal-closure order | cost `c_i` |
|---|---|---:|---:|---:|
| `6T12=A5` | `(2,2,1,1)` | 15 | 60 | 4 |
| `6T12=A5` | `(5,1)`, two classes | 12 each | 60 | forbidden |
| `6T14=S5` | `(2,2,1,1)` | 15 | 60 | 4 |
| `6T14=S5` | `(4,1,1)` | 30 | 120 | 4 |
| `6T14=S5` | `(5,1)` | 24 | 60 | forbidden |
| `6T15=A6` | `(2,2,1,1)` | 45 | 360 | 4 |
| `6T15=A6` | `(3,1,1,1)` | 40 | 360 | 3 |
| `6T15=A6` | `(5,1)`, two classes | 72 each | 360 | forbidden |
| `6T16=S6` | `(2,1,1,1,1)` | 15 | 720 | 2 |
| `6T16=S6` | `(2,2,1,1)` | 45 | 360 | 4 |
| `6T16=S6` | `(3,1,1,1)` | 40 | 360 | 3 |
| `6T16=S6` | `(3,2,1)` | 120 | 720 | 5 |
| `6T16=S6` | `(4,1,1)` | 90 | 720 | 4 |
| `6T16=S6` | `(5,1)` | 144 | 360 | forbidden |

The executable certificate retains conjugacy classes separately, including
the two five-cycle classes where applicable, before applying the index-five
exclusion.

## 4. The irreducible-branch obstruction

Let

\[
\rho:W_0\longrightarrow\mathbb A^2
\]

be the normalization of the affine target in the source function field.
Assume that its reduced branch locus is one irreducible curve `B` and that

\[
\nu:B^\nu\longrightarrow B
\]

is injective.  Then no degree-six Keller map can produce this configuration.

First, a ramified dicritical supplies a finite map `A1 -> B^nu`, so
`B^nu isomorphic to A1` by the argument in Section 1.  A finite injective
normalization of complex curves is a homeomorphism.  Hence `B` is
topologically contractible.

The Lin--Zaidenberg theorem, in Gurjar and Miyanishi's formulation, gives a
polynomial change of target coordinates after which

\[
B=\{X^m=Y^n\},\qquad \gcd(m,n)=1. \tag{4.1}
\]

The positive weighted dilations

\[
(X,Y)\longmapsto(r^nX,r^mY)
\]

identify the global complement with the product of a small link complement
and a radial interval.  They give the same product for the complement in a
small ball about the weighted vertex.  Thus the inclusion of the local
complement into the global complement induces an isomorphism on fundamental
groups.

Choose a point `x` on a ramified boundary prime over that vertex.  Its
Orevkov bracket contains

\[
e_E+(\mu_x-e_E)=\mu_x.
\]

All other terms are nonnegative and the total budget is five, so
`mu_x<=5`.  Constant multiplicity six then supplies another point of the
finite fiber.  Since `rho` is finite and proper, a sufficiently small target
ball has inverse image separated into disjoint neighborhoods of these fiber
points.  The induced six-sheet cover of the local branch complement is
disconnected, and its monodromy is intransitive.

On the other hand, `W_0` is integral.  The global cover over
`A2 - B` is connected, so its monodromy is transitive.  The local-to-global
fundamental-group isomorphism identifies the two monodromy images, a
contradiction.  Therefore

\[
\boxed{
\text{an irreducible finite-normalization branch curve must have}
\text{ noninjective normalization.}
}
\tag{4.2}
\]

Now define the generic ramified sheet contribution of an irreducible branch
curve by

\[
r_B=\sum_{\substack{E\to B\\e_E>1}}e_Ed_E.
\]

Every normalization point over a target value consumes at least `r_B` local
sheets.  If two normalization points collided, they would consume at least
`2r_B`.  Consequently

\[
2r_B>6\quad\Longrightarrow\quad \nu:B^\nu\to B
\text{ is injective}. \tag{4.3}

Combining (4.2) and (4.3), a sole ramified branch component cannot move four
or five sheets.

## 5. Elimination of `A5` and `S5`

For `A5`, the only permitted normally generating fixed-sheet class is the
double transposition.  It is realized either by one `(e,d)=(2,2)` dicritical
or by two `(2,1)` dicriticals dominating the same target curve.  In both cases
the refined cost and generic ramified sheet contribution are four.  Budget
five leaves no room for another ramified branch curve.  Thus the branch locus
is irreducible, while `2r_B=8>6` makes its normalization injective.  This
contradicts Section 4.

For `S5`, double transpositions normally generate only `A5`, and five-cycles
are forbidden.  Full `S5` monodromy therefore requires a `(4,1,1)` meridian,
realized by one `(e,d)=(4,1)` dicritical.  Its refined cost is four, so again
there is no second ramified branch curve.  Its normalization is injective
because it moves four sheets, giving the same contradiction.

Hence

\[
\boxed{6T12=A_5\text{ and }6T14=S_5\text{ are impossible}.} \tag{5.1}
\]

The same argument removes the double-transposition profile from `A6`, and it
removes both the lone four-cycle profile and the same-branch `(3)(2)` profile
from `S6`.

## 6. Exact surviving ramified profiles

The normally generating branch-class multisets of refined cost at most five
are exhaustively enumerated by
[`ramified_branch_profiles`](../scripts/six_sheet_monodromy.py).  After the
irreducible-branch obstruction, the survivors are:

| group | distinct ramified target curves | minimum cost | residual budget |
|---|---|---:|---:|
| `A6` | one `(3,1,1,1)` curve | 3 | 2 |
| `S6` | one transposition curve | 2 | 3 |
| `S6` | two distinct transposition curves | 4 | 1 |
| `S6` | distinct transposition and 3-cycle curves | 5 | 0 |

For `A6`, the unique ramified dicritical has `(e,d)=(3,1)`.  The reduced
branch locus is one irreducible curve `B`, so (4.2) forces
`B^nu -> B` to be noninjective.  Every normalization collision consists of
exactly two boundary points, each of local multiplicity three.  Orevkov's
zero-jump lemma makes both branches smooth; they exhaust all six sheets, so
the collision value is omitted by the original affine map.

For `S6`, every surviving profile contains an actual transposition branch.
In the two-curve profiles the target curves are distinct.  If an index-two and
an index-three dicritical dominated the same curve, the generic inertia would
be the eliminated `(3,2,1)` class.

The residual budget in the table is available for `e=1` dicriticals and for
the excesses `delta_E`; it is not evidence that any profile is geometrically
realizable.

Under the additional one-dicritical hypothesis, the `A6` residual budget two
cannot split into two excess-one points.  The proof and the forced smooth
`(2,5)` cusp are in
[`a6-one-dicritical-local.md`](a6-one-dicritical-local.md).  The
Hirzebruch--Jung lift and parity argument in
[`a6-exceptional-source.md`](a6-exceptional-source.md) eliminates the
contracted singular-source alternative.

For the one-dicritical `S6` survivor, the complete local fiber census and
smooth-source torus-knot classification are in
[`s6-one-dicritical-local.md`](s6-one-dicritical-local.md).

For the two-curve `S6` survivors, Nguyen Van Chau's singular-component
theorem forces at least one zero-jump transposition curve to self-collide.  In
the saturated transposition-plus-3-cycle profile both branch curves
self-collide, and constant multiplicity gives the exact rows `3+3`,
`2+2+1+1` or `2+2+2`, and `2+3+1` at cross-intersections.  See
[`s6-two-curve-collisions.md`](s6-two-curve-collisions.md).

## 7. Why this still does not finish degree six

The surviving actions are primitive.  Their point stabilizers

\[
A_5<A_6,\qquad S_5<S_6
\]

are maximal, so the degree-six function-field extension has no proper
intermediate field.  There is no degree-two or degree-three factor to which a
known low-sheet theorem can be applied.

Ordinary permutation Riemann--Hurwitz also does not finish the problem: both
`A6` and `S6` possess transitive genus-zero branch-cycle triples.  Such
abstract rational covers do not satisfy the Keller compactification, but they
show that product-one, transitivity, primitivity, and genus alone are too
coarse.  The missing input is the full boundary geometry: canonical labels,
determinant labels, adjacency, and compatibility among the rational maps over
the boundary graph, as in the Domrina--Orevkov and Borisov frameworks.

## 8. Reproduction

Run the finite certificate with:

```bash
uv run python scripts/six_sheet_monodromy.py
uv run pytest tests/test_six_sheet_monodromy.py
```

The tests lock:

- the exact fixed-sheet conjugacy-class table;
- the equality between moved sheets and minimum refined cost;
- all normally generating class multisets of total cost at most five;
- the irreducible-branch filter on those profiles;
- the final unrestricted group list `6T15`, `6T16`.

Sage/GAP remains an optional independent cross-check; the runtime certificate
has no GAP dependency.

## Primary sources

- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587–596.  Lemma 2.1 supplies
  the affine-line dicriticals, Lemma 3.1 supplies the generic local form, and
  Lemma 4.2 supplies the exact defect identity.
- R. V. Gurjar and M. Miyanishi,
  [“On contractible curves in the complex affine plane”](https://doi.org/10.2748/tmj/1178225344),
  *Tohoku Mathematical Journal* 48 (1996), 459–469.  Theorem 2 is the
  Lin--Zaidenberg normal form used in Section 4.
- Alexander Borisov,
  [“Frameworks for two-dimensional Keller maps”](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v27i3p54/pdf/),
  *Electronic Journal of Combinatorics* 27 (2020), P3.54, for the remaining
  compactification-label and boundary-graph framework.

Identity (0.2), the support-budget application, and the resulting degree-six
group elimination are derived in this repository.  No claim of historical
priority is made.
