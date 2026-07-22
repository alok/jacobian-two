# True-split rank audit for `T112` and `C2+T111`

This note closes the coefficient-rank question for seven allocation charts in
the conditional `A6`, delta-ten filtration:

- `T112+6N`, tangent edge on `V` or `W`, at `k=0` and `k=2`; and
- `C2+T111+5N`, separate contact on `W` at `k=0`, or on `V` or `W` at
  `k=2`.

The result is deliberately local to these true-split incidence charts.  It
does not compute their affine complement groups, classify every intersection
with adjacent profiles, treat the other fifteen split allocation rows, or
establish an unrestricted Jacobian-conjecture consequence.

## 1. Why the component equations must be rebuilt

For

\[
 P(t)=t^2+k t^3+t^4,
 \qquad
 Q(t)=a t^5+b t^6+c t^7+d t^8+t^9,
\]

the unordered-pair incidence is reducible at `k=0,+2,-2`.  Specializing the
cancelled generic collision decic is not a safe way to recover its
multiplicities or coefficient-rank strata.  The audit instead uses the true
collision factors

\[
 V_k(r),\qquad W_k(s),
\]

of degrees `(2,8)` at `k=0` and `(4,6)` at `k=2`.

At `k=0`, label a four-source `P`-fiber by

\[
 (z,-z,w,-w),\qquad z^2+w^2+1=0.                    \tag{1.1}
\]

The vertical triple edge has `r=-z^2`; the two graph edges have
`s=z+w` and `s=w-z`.  At `k=2`, use

\[
 (z,-1-z,w,-1-w),\qquad z(z+1)+w(w+1)=0.           \tag{1.2}
\]

The vertical root is `r=-z(z+1)` and the two graph roots are `s=z+w` and
`s=w-z-1`.

The executable certificate verifies the following identities exactly.  At
`k=0`,

\[
 Q(z)-Q(-z)=2z^5V_0(-z^2),
\]

while at `k=2`,

\[
 Q(z)-Q(-1-z)=(2z+1)V_2(-z(z+1)).
\]

On either triple base,

\[
 16\bigl(Q(z)-Q(w)\bigr)
  =(z-w)m_k(z+w)W_k(z+w),                            \tag{1.3}
\]

where `m_0(s)=1` and `m_2(s)=s^2`; the analogous identity holds for the
second graph edge.  Thus `V=0` plus one `W=0` is exactly the pair of
independent triple-equality equations on the displayed open.

Differentiating `V` or `W` along its pair component gives the tangency
equation.  This is legitimate only when the pair is distinct, the component
target coordinate is unramified, and the collision factor used in (1.3) is
nonzero.  Those factors are included explicitly in the localizers below;
the calculation never silently divides by them.

## 2. Exact incidence systems

Write

\[
 E_V=V_k(r_T),\qquad E_W=W_k(s_T),
\]

for the vertical and selected graph edge of the triple.

The four `T112` systems are

\[
 (E_V,E_W,\partial_rV_k(r_T))=0
\]

for a tangent vertical edge, and

\[
 (E_V,E_W,\partial_sW_k(s_T))=0
\]

for a tangent graph edge.  Each is a three-row affine-linear system in
`(a,b,c,d)` over the one-dimensional triple base.

For a separate contact parameter `h`, the three mixed systems are

\[
 (E_V,E_W,V_k(h),\partial_rV_k(h))=0                 \tag{2.1}
\]

or

\[
 (E_V,E_W,W_k(h),\partial_sW_k(h))=0.                \tag{2.2}
\]

They are four-row affine-linear systems over the two-dimensional base given
by the triple conic and `h`.

The base localizers remove only genuine invalid boundaries.  For example,
the `k=0` triple localizer is

\[
 z w(z-w)(z+w)(2z^2+1),                              \tag{2.3}
\]

and the `k=2` one is

\[
 z(z+1)w(w+1)(2z+1)(2w+1)(z-w)(z+w+1).              \tag{2.4}
\]

For a separate graph contact at `k=0`, the additional factor is

\[
 h(h^2+1)(h^2+2)
 (h^2-(z+w)^2)(h^2-(w-z)^2).                         \tag{2.5}
\]

The first three factors remove the component overlap, ramified target, and
pair diagonal; the last two remove equality with the triple target.  The
certificate proves the target-difference factorization rather than treating
it as a heuristic.

At `k=2`, the vertical-contact factor is

\[
 h(1-4h)\bigl(h^2-z^2(z+1)^2\bigr),                 \tag{2.6}
\]

and the graph-contact factor, with `q1=z+w+1` and `q2=w-z`, is

\[
 h(h+1)(h+2)\bigl(2(h+1)^2-1\bigr)
 \bigl((h+1)^2-q_1^2\bigr)
 \bigl((h+1)^2-q_2^2\bigr).                         \tag{2.7}
\]

No coefficient determinant, augmented determinant, or residual
compatibility factor occurs in these localizers.

## 3. Saturated rank results

The independent Sage checker constructs every coefficient and augmented
matrix from the component formulas, cross-checks it entry by entry against
the typed Python module, and performs the following exact ideal saturations.

| chart | generic rank | valid coefficient-rank drop | valid compatibility | largest residual incidence |
|---|---:|---|---|---:|
| `T112`, `k=0`, tangent `V` | 3 | empty, exponent 1 | empty, exponent 1 | none |
| `T112`, `k=0`, tangent `W` | 3 | empty, exponent 1 | empty, exponent 1 | none |
| `T112`, `k=2`, tangent `V` | 3 | empty, exponent 2 | empty, exponent 2 | none |
| `T112`, `k=2`, tangent `W` | 3 | empty, exponent 2 | empty, exponent 2 | none |
| mixed, `k=0`, contact `W` | 4 | one rank-three curve | reduced length 4 | dimension 1 |
| mixed, `k=2`, contact `V` | 4 | one rank-three curve | reduced length 4 | dimension 1 |
| mixed, `k=2`, contact `W` | 4 | two rank-three curves | empty, exponent 2 | none |

For every mixed chart, all coefficient `3 x 3` minors and all augmented
`3 x 3` minors saturate to the unit ideal with exponent one.  Hence no valid
base has coefficient rank at most two.  The surviving compatibility points
therefore have coefficient rank exactly three and affine-line solution
fibers.

The two nonempty reduced compatibility schemes have exact Gröbner bases

\[
 \left(z^2+\frac{11}{12},\ h^2+\frac13,\ w-\frac h2\right)             \tag{3.1}
\]

for the `k=0/W` chart, and

\[
 \left(
 z^2+z-h+\frac14,\quad h^2-h+\frac12,\quad w-h+\frac32
 \right)                                                               \tag{3.2}
\]

for the `k=2/V` chart.  Each has length four and is radical.  Adding its
one-dimensional coefficient kernel gives total residual incidence dimension
one.  In particular, neither can be a hidden threefold.

On the maximal-rank opens, a `T112` system has a one-dimensional coefficient
fiber over a curve, while a mixed system has a unique coefficient point over
a surface.  Both generic incidences therefore have dimension two.  Combining
the maximal-rank and residual bounds proves:

> None of the seven audited true-split allocation charts supports a
> three-dimensional incidence component.

## 4. Hostile boundaries and `k=-2`

Every system has a compatible rank-two point before localization.  Examples
include

- `z=-i/sqrt(2), w=i/sqrt(2)` for the `k=0` vertical-tangent system;
- `z=0, w=i` for the `k=0` graph-tangent system;
- `z=0, w=-1` for both `k=2` tangent systems; and
- `(z,w,h)=(0,i,i)`, `(0,-1,0)`, and `(0,-1,-1)` for the three mixed
  systems.

In every case the coefficient and augmented ranks are both two, but the
valid localizer is zero: a source repeats, a critical source is used, a
component overlap occurs, or the contact is not separate.  These hostile
fixtures demonstrate that the raw rank ideals cannot be promoted without
localization.

The `k=-2` conclusion is not inferred merely from equality of collision
degrees.  The checker applies

\[
 t\longmapsto-t,
 \qquad (a,b,c,d)\longmapsto(a,-b,c,-d),              \tag{4.1}
\]

to all four `k=2` systems.  It verifies the triple-base equation, both full
polynomial-family identities, every `V/W` value row, and every component
derivative row.  Graph parameters change by `s -> -s`; vertical parameters
are fixed.  Thus the rank and compatibility results transport to `k=-2` for
the full incidence equations.

## 5. Executable artifacts and remaining boundary

Run the typed certificate and focused tests with

```bash
uv run python -m scripts.a6_delta_ten_split_t112_mixed_rank
uv run pytest -q tests/test_a6_delta_ten_split_t112_mixed_rank.py
```

Run the independent ideal calculation with Sage 10.8:

```bash
sage tools/check_a6_delta_ten_split_t112_mixed_rank.sage
```

The focused module deliberately proves only these seven rank statements.
The companion contact and global-fiber certificates now treat the other
fifteen rows, while `a6_delta_ten_split_rank_all` checks exact key-for-key
coverage of the generated 22-orbit-type ledger, representing 33 actual rows.

There is also a separate component-closure certificate.  For the mixed
profile, retain both the sum `s` and product `p` of the contact pair.  With
triple coordinates `q,r`, put `n=r^2+q^2-q^3`.  The total source base is

\[
 (2qrs+n)p-s\bigl(qr(s^2+1)+ns\bigr)=0.             \tag{5.1}
\]

It is primitive and linear in `p`: its two coefficients have gcd one.  Hence
the total base is geometrically irreducible and includes the `k=2` vertical
pair-denominator chart instead of deleting it.  Exact invertible row
transformations compare the total target/tangency rows with all three split
mixed systems.  Analogous transformations identify all four split `T112`
systems with restrictions of the existing irreducible labeled incidence.
Independent Rabinowitsch saturations prove that all seven transformation
determinants are units on their full valid split localizations; the stored
clean witnesses separately prove nonemptiness.  The total mixed base is
smooth after localizing by `q*r`, and an exact four-row diagonal bridge
recovers the historical mixed Cramer system on its ordinary chart.
Consequently every clean maximal-rank locus in these seven allocations lies
in the same connected equisingular incidence as the cyclic sample, and the
existing proper-isotopy exclusion applies.

What remains here is smaller but real: the two reduced length-four compatible
mixed orbit representatives have rank-three affine-line coefficient fibers
whose topology has not been connected to the global graph.  The `k=2`
representative transports to a third actual scheme at `k=-2`, so the actual
mixed total length is twelve rather than the representative length eight.
Deeper intersections and
assignment of lower-dimensional compatible pieces to adjacent profiles also
remain open.  Neither this note nor the aggregate rank theorem proves
`JC(2)`.
