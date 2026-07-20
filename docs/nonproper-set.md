# Cubic fibers, the omitted curve, and the exact nonproper-value set

## Claim boundary

This note derives exact consequences of the announced three-dimensional map.
The polynomial identities are independently checked by Lean and SymPy; the
fiber count and sequential nonproperness arguments are written out here.

The same-day public record already contains versions of the cubic,
reconstruction, discriminant, and generic `S_3` calculation.  The complete
fiber stratification, omitted curve, and equality of the nonproper-value set
with the discriminant hypersurface were not found in an authoritative earlier
source.  They are therefore labeled **derived here; historical priority not
established**, not advertised as literature-first results.

Write the announced map as

\[
 F=(A,B,C):\mathbb C^3_{x,y,z}\longrightarrow\mathbb C^3_{a,b,c}.
\]

## The reciprocal fiber coordinate

On `x != 0`, put

\[
 t=y+\frac1x.
\]

Exact simplification of the displayed coordinate formulas gives

\[
 b=4t+\frac2x-3ct^2,
 \qquad
 2a=ct^3-2t^2+bt.
\]

Consequently `t` is a root of

\[
 p_{a,b,c}(T)=cT^3-2T^2+bT-2a,
\]

and

\[
 r=p'(t)=3ct^2-4t+b=\frac2x.
\]

Conversely, a finite root `t` with `r != 0` reconstructs one source point:

\[
 x=\frac2r,
 \qquad
 y=t-\frac r2,
 \qquad
 z=\frac54r^2-\frac32tr-\frac c8r^3.
\]

Direct substitution sends this point to `(a,b,c)`.  It is the unique source
associated to that root.

## The projective root at infinity

Homogenize the cubic:

\[
 \bar p(T,U)=cT^3-2T^2U+bTU^2-2aU^3.
\]

The projective point `[1:0]` is a root exactly when `c=0`.  In the local
coordinate `w=U/T`, the polynomial is

\[
 c-2w+bw^2-2aw^3,
\]

whose derivative at `w=0` is `-2`; the root at infinity is always simple.

On the corresponding source chart,

\[
 F(0,y,z)=(z+4y^2,y,0).
\]

Thus the infinity root over `(a,b,0)` corresponds to the finite source point

\[
 (0,b,a-4b^2).
\]

It follows that the fiber of `F` is in exact bijection with the simple
projective roots of `\bar p`.

## Discriminant, generic Galois group, and singular curve

The universal degree-three coefficient discriminant—equivalently, the
discriminant of the homogenized binary cubic—is

\[
 \operatorname{Disc}_T(p)=-4Q(a,b,c),
\]

where

\[
 Q=27a^2c^2-18abc+16a+b^3c-b^2.
\]

When `c=0`, the dehomogenized polynomial has degree two.  Its ordinary
quadratic discriminant is `b^2-16a`, whereas the specialized binary-cubic
discriminant above is `4(b^2-16a)`.  The factor comes from the additional
simple projective root at infinity; the vanishing locus is the same.  Over the
generic function field `c` is nonzero, so the displayed formula is also the
ordinary cubic discriminant there.

### Generic field extension and Galois closure

Let `K=C(a,b,c)` and `L=C(x,y,z)`.  The nonzero Jacobian makes `a,b,c`
algebraically independent.  The reconstruction formulas show that

\[
  L=K(T).
\]

The cubic `p` is irreducible over `K`.  A short algebraic proof is to put
`k=C(b,c)` and work in `k[a,T]`.  The polynomial is primitive as a polynomial
in `T` because its `T^2` coefficient is the unit `-2`, and

\[
  k[a,T]/(p)\cong k[T],
  \qquad
  a\longmapsto\frac{cT^3-2T^2+bT}{2}.
\]

Thus `(p)` is prime, so Gauss's lemma gives irreducibility in `K[T]` and
`[L:K]=3`.

The polynomial `Q`, viewed as a quadratic in `a`, has discriminant

\[
  \operatorname{Disc}_a(Q)=-4(3bc-4)^3.
\]

The odd exponent makes this nonsquare in `C(b,c)`, and Gauss's lemma shows
that `Q` is an irreducible polynomial.  Hence `-4Q` is not a square in `K`.
An irreducible cubic with nonsquare discriminant has Galois group `S_3`, so
the splitting field has degree six while the cubic source extension `L/K` is
nonnormal and has no nontrivial `K`-automorphism.

This is the precise Galois-theoretic behavior of the example.  It does not
contradict the classical Galois-case theorem: the source extension is not
Galois.  “Galois group `S_3`” here means the function-field splitting field,
or generic geometric monodromy—not the Galois group of a numerical cubic over
the algebraically closed field `C`.

Let `Sigma=V(Q)`.  Its singular locus is

\[
 \Gamma=\{3bc=4,\ b^2=12a\}
 =\left\{\left(\frac{t^2}{3},2t,\frac{2}{3t}\right):t\in\mathbb C^*\right\}.
\]

For completeness, the partial derivatives are

\[
 Q_a=54ac^2-18bc+16,
\]

\[
 Q_b=-18ac+3b^2c-2b,
 \qquad
 Q_c=54a^2c-18ab+b^3.
\]

If all three vanish, then `c != 0`, since `c=0` would make `Q_a=16`.
Writing `lambda=bc` and `mu=ac^2`, the first two equations give

\[
 54\mu-18\lambda+16=0,
 \qquad
 -18\mu+3\lambda^2-2\lambda=0.
\]

Eliminating `mu` yields `(3 lambda-4)^2=0`, hence `3bc=4`; substituting back
gives `b^2=12a`.  Conversely, substitution of these two equations makes
`Q,Q_a,Q_b,Q_c` all zero.

As an independent exact check, let
`J=(Q,Q_a,Q_b,Q_c)` in `Q[a,b,c]` and use lexicographic order `a>b>c`.
The test suite computes the following Gröbner basis, up to unit scaling:

\[
 \left(16a-b^3c,\ (3bc-4)^2\right)
\]

This basis is deliberately nonradical.  Over `C`,

\[
  \sqrt J=(12a-b^2,\ 3bc-4),
\]

which is the ideal of `Gamma`.

## Complete fiber counts and image

If `Q != 0`, the projective cubic has three distinct roots, all simple.  If
`Q=0` but the cubic is not a cube, its multiplicity pattern is `2+1`, so it
has exactly one simple root.  A triple root cannot occur at infinity because
the infinity root is always simple.  Matching

\[
 \bar p(T,U)=c(T-tU)^3
\]

coefficient by coefficient shows that the triple-root locus is precisely
`Gamma`.

Therefore

\[
 \#F^{-1}(a,b,c)=
 \begin{cases}
 3,&Q(a,b,c)\ne0,\\
 1,&Q(a,b,c)=0\text{ and }(a,b,c)\notin\Gamma,\\
 0,&(a,b,c)\in\Gamma.
 \end{cases}
\]

In particular,

\[
 F(\mathbb C^3)=\mathbb C^3\setminus\Gamma.
\]

The screenshot target `(-1/4,0,0)` has `Q=-4`, so its three displayed
preimages are the generic three-root case, not a point on the asymptotic
hypersurface.

## Every point of `Q=0` is an asymptotic value

Every point of `Sigma` has a repeated finite root `t`.  Solving `p(t)=p'(t)=0`
gives a polynomial parametrization

\[
 q_0(t,c)=(t^2-ct^3,\ 4t-3ct^2,\ c).
\]

For nonzero `epsilon`, define

\[
 X_\varepsilon=
 \left(
 \frac2\varepsilon,
 t-\frac\varepsilon2,
 \frac54\varepsilon^2-\frac32t\varepsilon-\frac c8\varepsilon^3
 \right).
\]

Exact substitution gives

\[
 F(X_\varepsilon)=
 \left(
 t^2-ct^3+\frac{\varepsilon t}{2},
 4t-3ct^2+\varepsilon,
 c
 \right).
\]

As `epsilon -> 0`, the source escapes because its first coordinate is
`2/epsilon`, while its image tends to `q_0(t,c)`.  Hence

\[
 V(Q)\subseteq S_F,
\]

where `S_F` denotes the set of limits of images of source sequences escaping
to infinity.

## Nothing outside `Q=0` is an asymptotic value

Suppose, toward a contradiction, that a source sequence escapes,

\[
 \lVert(x_n,y_n,z_n)\rVert\longrightarrow\infty,
\]

while

\[
 F(x_n,y_n,z_n)=(a_n,b_n,c_n)\longrightarrow(a,b,c)
\]

with `Q(a,b,c) != 0`.

If `x_n=0` along an infinite subsequence, then

\[
 y_n=b_n,
 \qquad
 z_n=a_n-4b_n^2,
\]

so that source subsequence is bounded, contradicting escape to infinity.

Otherwise pass to a subsequence with `x_n != 0`, and put

\[
 t_n=y_n+\frac1{x_n}.
\]

The projective roots `[t_n:1]` have a convergent subsequence in
`\mathbb P^1(\mathbb C)`.  Their limit is a projective root of the limiting
cubic.  Since `Q != 0`, it is simple.

If the limiting root is finite, then `r_n=p'_n(t_n)` tends to a nonzero
number.  The reconstruction formulas make that source subsequence bounded,
again a contradiction.

If the limiting root is infinity, then `t_n != 0` eventually.  Moreover,
`1+x_ny_n=x_nt_n != 0`, so the reciprocal chart applies.  Set

\[
  s_n=\frac1{t_n}=\frac{x_n}{1+x_ny_n}\longrightarrow0.
\]

With

\[
 D_n=1-b_ns_n+3a_ns_n^2,
\]

the reciprocal reconstruction identities are

\[
 x_n=\frac{s_n}{D_n},
 \qquad
 y_n=b_n-3a_ns_n,
\]

\[
 z_n=a_nD_n^3-y_n^2D_n(D_n+3).
\]

Here `D_n -> 1`, so in fact

\[
 (x_n,y_n,z_n)\longrightarrow(0,b,a-4b^2).
\]

This source subsequence converges, again contradicting escape.  Thus no source
sequence can escape over a target outside `V(Q)`, and

\[
 \boxed{S_F=V(Q).}
\]

The hypersurface `V(Q)` is not a critical-value locus for `F`: the polynomial
Jacobian of `F` is the everywhere nonzero constant `-2`.  It is the universal
cubic discriminant locus and, for this nonfinite polynomial map, exactly where
sheets escape to infinity.  This note does not separately prove a
scheme-theoretic branch-divisor statement for the normalization or Galois
closure.

## Reproduction

The exact independent checker verifies the discriminant, both reconstruction
charts, the repeated- and triple-root parameterizations, and the escaping
family.  The pytest suite separately checks the singular-ideal Gröbner basis
and hostile perturbations:

```bash
uv run --frozen python -m scripts.nonproper
uv run --frozen pytest tests/test_nonproper.py
```

The Lean module `JacobianTwo/CubicFiber.lean` certifies the cubic and formal
derivative, the standard universal discriminant coefficient expression, a
Bézout common-root certificate,
finite-root reconstruction, the infinity source, repeated- and triple-root
parameterizations, the exact escaping-family image, and the large-`T`
cancellations.  It does **not** formalize the forward source-to-projective-root
map, projective root counts, singular-locus radical, or sequential compactness.
Those steps are exposed in this note rather than mislabeled as kernel-checked
algebraic geometry or topology.

## Same-day comparison sources

- [MathOverflow question on the cubic and `S_3` structure](https://mathoverflow.net/questions/513387/galois-structure-of-the-new-counterexample-to-the-jacobian-conjecture-an-explic)
- [Interactive counterexample explainer](https://jacobianfun.org/jacobian-explained)
