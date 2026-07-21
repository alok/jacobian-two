# The arbitrary-degree affine-coordinate theorem

Status: **Lean-certified.** The mathematics below is a known positive class,
not a solution of the plane Jacobian conjecture and not a novelty claim.

## Statement

Let `K` be a field of characteristic zero. Write a polynomial in two
variables as a polynomial in `y` with coefficients in `K[x]`, and consider

\[
  F=(P,Q),\qquad P\in K[x,y],\qquad Q=e(x)y+f(x).
\]

Use the Jacobian convention

\[
  J(P,Q)=P_xQ_y-P_yQ_x.
\]

If `J(P,Q)=k` for a nonzero `k in K`, then `F` is a polynomial
automorphism. More precisely, exactly one of the following charts applies.

### Nonzero `y`-slope

If `e` is not the zero polynomial, the Jacobian equation forces it to be a
nonzero scalar. There are

\[
  \varepsilon,\alpha\in K^\times,\qquad \beta\in K,
  \qquad G\in K[T]
\]

such that

\[
  e=\varepsilon,\qquad
  P=G(Q)+\alpha x+\beta,\qquad
  \alpha\varepsilon=k.
\]

Writing the target coordinates as `(u,v)`, an inverse is

\[
  x=\frac{u-G(v)-\beta}{\alpha},\qquad
  y=\frac{v-f(x)}{\varepsilon}.
\]

The harmless ambiguity between the constant term of `G` and `beta` can be
removed by requiring `G(0)=0`.

### Zero `y`-slope

If `e=0`, there are

\[
  \eta,\delta\in K^\times,\qquad \gamma\in K,
  \qquad c\in K[x]
\]

such that

\[
  P=\eta y+c(x),\qquad
  Q=\delta x+\gamma,\qquad
  -\eta\delta=k.
\]

This chart has inverse

\[
  x=\frac{v-\gamma}{\delta},\qquad
  y=\frac{u-c(x)}{\eta}.
\]

## Degree descent

Write

\[
  P(x,y)=\sum_{j=0}^m p_j(x)y^j,
  \qquad p_m\ne0.
\]

Since

\[
  Q_x=e'(x)y+f'(x),\qquad Q_y=e(x),
\]

the coefficient of `y^j` in the Jacobian is

\[
  ep_j'-j e'p_j-(j+1)f'p_{j+1}.
\]

For positive `m`, the top coefficient must vanish:

\[
  e p_m'=m e'p_m. \tag{1}
\]

The one-variable differential lemma used here says that, over a
characteristic-zero field, (1) and `e != 0` imply

\[
  p_m=\lambda e^m
\]

for some scalar `lambda`. This can be proved by comparing degrees and leading
coefficients; it needs no algebraic closure, root factorization, or
squarefreeness assumption.

The leading `y`-coefficient of `Q^m` is `e^m`, while

\[
  J(Q^m,Q)=0.
\]

Therefore replacing `P` by `P-lambda*Q^m` lowers its `y`-degree without
changing its Jacobian. Strong induction ends with a remainder `h(x)` and
gives

\[
  P=G(Q)+h(x),\qquad e(x)h'(x)=k.
\]

The last equality is the decisive Keller rigidity: a product of two
polynomials equal to a nonzero scalar has two scalar factors. Hence `e` is a
nonzero scalar and `h=alpha*x+beta`, which proves the first normal form.

When `e=0`, the full Jacobian equation is already

\[
  -P_y f'=k.
\]

Both factors must be nonzero scalars. The characteristic-zero derivative
kernel then gives the second normal form directly.

## What Lean certifies

[`JacobianTwo/AffineCoordinate.lean`](../JacobianTwo/AffineCoordinate.lean)
models `K[x,y]` as `(K[X])[Y]` and defines both formal partial derivatives as
bundled derivations. Its premise is an equality of the actual formal
Jacobian, not a manually supplied list of coefficient equations.

The certificate covers:

- the leading-coefficient identity extracted from the formal Jacobian;
- the one-variable differential lemma;
- strong degree descent with exact leading-term cancellation;
- both scalar normal forms, including `e=0`;
- explicit left and right inverse identities; and
- bijectivity of the original map over every characteristic-zero field.

The independent exact checker
[`scripts/affine_coordinate.py`](../scripts/affine_coordinate.py) exercises a
degree-seven outer polynomial, both inverse compositions, the complementary
triangular chart, and a nonconstant-slope perturbation whose determinant is
`5*x+15` rather than a scalar.

## Literature boundary

Marco Sabatini's Theorem 4 treats real analytic maps of type `(m,1)` whose
Jacobian is independent of `y`. Its proof performs the same induction: the
top relation gives `p_m=c_m*q_1^m`, and the triangular target shear

\[
  (u,v)\longmapsto(u-c_m v^m,v)
\]

removes the top power. See the
[published journal record](https://www.impan.pl/en/publishing-house/journals-and-series/colloquium-mathematicum/all/175/1/115518/global-injectivity-of-planar-non-singular-maps-that-are-polynomial-in-one-variable),
[DOI](https://doi.org/10.4064/cm9195-1-2024), and
[author preprint](https://arxiv.org/abs/2302.05394).

Accordingly, the underlying positive class is known. The artifact supplied
here is a field-uniform formalization of its constant-Jacobian specialization,
packaged as one normal-form theorem with explicit inverses. We do not import
stronger real-variable conclusions from the paper: the repository theorem
uses the genuinely constant identity `J(P,Q)=k`, not merely a determinant
that is nowhere zero on the real plane.

## Search consequence and exact limit

Every map in this class has generic function-field degree one. Thus a
hypothetical generic-degree-six plane counterexample cannot have a target
coordinate of the form `e(x)y+f(x)` in the chosen source coordinates.

The coordinate-invariant sieve is slightly stronger. If an invertible source
coordinate change and an invertible linear target change make any target
linear combination affine in one source variable, this theorem applies and
the map is an automorphism.

This does **not** solve the generic-degree-six stratum. No argument here says
that an arbitrary degree-six Keller map admits such a coordinate. The next
research problem is to find a broader, invariant condition at infinity that
either forces this affine-coordinate reduction or excludes its first
genuinely nonlinear alternative.

## Sharp hypotheses

- Characteristic zero is essential. In characteristic `p`,
  `(x-x^p,y)` has Jacobian one but is not an automorphism.
- `k != 0` is essential. At Jacobian zero, nonconstant `e` occurs freely.
- In the nonzero-slope chart, `P` may initially be independent of `y`; the
  theorem then forces `P=(k/eps)*x+beta` up to a constant `G`.
- In the zero-slope chart, `P` cannot be independent of `y`, because that
  would make the Jacobian zero.
- Repeated roots of `e` are not a hidden case: the proof never makes a
  rootwise argument, and the terminal unit equation excludes every
  nonconstant `e` at once.
