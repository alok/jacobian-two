# Why the hostile `S6` trace curve is globally impossible

## Claim boundary

The symmetric trace matrix in
[`trace-hostile-matrices.md`](trace-hostile-matrices.md) is an exact quadratic
consistency model.  Its particular rational branch curve, however, cannot be
the branch curve of a connected six-sheet cover with transposition inertia.
This eliminates that explicit curve and matrix as a global `S6` candidate; it
does **not** eliminate the one-dicritical `S6` passport or the degree-six
Jacobian frontier.

## 1. The curve has a four-strand projection

The normalization is

\[
  P=t^4-6t^2,\qquad Q=t^5-5t^3,
\]

with implicit equation

\[
\begin{aligned}
b(P,Q)={}&P^5+10P^4+25P^3+10P^2Q^2\\
&+90PQ^2-Q^4+216Q^2.
\end{aligned}
\]

Projection to the `P`-axis has degree four: `P(t)` is a monic quartic, and
equivalently `b` has degree four in `Q` with constant leading coefficient.
Zariski--van Kamp therefore presents the affine complement using four fiber
meridians.

There is also a presentation-independent obstruction.  If meridians
\(g_1,\ldots,g_w\) map to transpositions in `S_n`, regard each transposition
as an edge on the `n` sheets.  The generated group preserves every connected
component of this edge graph.  A transitive action therefore requires the
graph to be connected, and hence

\[
  w\ge n-1.
\]

Here `w=4` and `n=6`, so transitivity is impossible.

This argument is specific to the `S6` survivor because its generic branch
inertia is a single transposition.  It is not an obstruction to `A6`, whose
generic inertia consists of a 3-cycle.

## 2. Exact Zariski--van Kamp replay

As an independent audit, Sage's Zariski--van Kamp implementation gives four
generators `x0,x1,x2,x3` and the six relators encoded in
[`s6_trace_curve_topology.py`](../scripts/s6_trace_curve_topology.py).  The
source computation is reproducible with:

```sage
from sage.schemes.curves.zariski_vankampen import fundamental_group

R.<P,Q> = QQ[]
b = (P^5 + 10*P^4 + 25*P^3 + 10*P^2*Q^2
     + 90*P*Q^2 - Q^4 + 216*Q^2)
G = fundamental_group(b, simplified=True, projective=False, puiseux=True)
print(G.relations())
```

The dependency-free checker then exhausts all

\[
  15^4=50{,}625
\]

assignments of the four generators to transpositions of `S6`.  Exactly `735`
satisfy all six relators.  Their generated-group orders are:

| group order | assignments |
|---:|---:|
| 2 | 15 |
| 120 | 720 |

None is transitive.  The order-120 images are conjugate sheet-stabilizing
copies of `S5`, exactly as the four-edge theorem predicts.

## 3. Consequence for the trace attack

The matrix still proves what it was built to prove: determinant,
normalization-cokernel, corank, trace-unit, symmetry, and theta data do not by
themselves encode connected six-sheet monodromy.  The new topological check
shows that a useful hostile trace model must additionally satisfy a global
complement-group constraint.  For a one-component `S6` branch, any polynomial
normalization used in a future model must have affine projection width at
least five.

The next question was therefore not whether this matrix secretly defines a
Keller cover—it does not—but what happens at width five.  The first exact
degree-minimal curve at equality is also cyclic and hence impossible, as shown
in the [width-five near-miss note](s6-width-five-near-miss.md).  That result
still does not classify every rational branch curve compatible with the forced
`S6` cusp/collision budget and width at least five.
