# Exact audit and structural analysis

## 1. Scope and verdict

This audit separates a finite calculation from a much larger historical
claim.

- The displayed formula, determinant, and three-point collision are exact
  algebraic statements. Lean and SymPy both verify them.
- Those statements refute the complex Jacobian conjecture in dimension three.
- Padding with identity coordinates refutes it in every dimension `n >= 3`.
- They do not imply a counterexample in dimension two. The plane case remains
  the research target of this repository.
- The attribution and discovery story come from a same-day public
  announcement, not yet a mature peer-reviewed record.

Multiplying the first output coordinate by `-1/2` changes the determinant from
`-2` to `1` and preserves noninjectivity. Thus the example also has the usual
determinant-one normalization.

## 2. Literal transcription

Write `F = (P, Q, R)`, where

```text
P = (1 + x*y)^3*z + y^2*(1 + x*y)*(4 + 3*x*y)
Q = y + 3*x*(1 + x*y)^2*z + 3*x*y^2*(4 + 3*x*y)
R = 2*x - 3*x^2*y - x^3*z.
```

The formal determinant is obtained by differentiating all three coordinate
polynomials with respect to `x`, `y`, and `z`, then expanding the determinant.
The result is exactly `-2`; no evaluation or random sampling is involved.

The displayed points and target are

```text
p0 = ( 0,    0,   -1/4)
p1 = ( 1, -3/2,   13/2)
p2 = (-1,  3/2,   13/2)
t  = (-1/4, 0,       0).
```

Direct rational substitution gives `F(p0) = F(p1) = F(p2) = t`.

## 3. Why the determinant collapses

There is useful structure behind the large expanded identity. Set

\[
u=1+xy,\qquad H=u^2z+y^2(4+3xy),\qquad s=\frac{x}{u}.
\]

On the rational chart `u != 0`, the first two outputs become

\[
P=uH,\qquad Q=y+3sP.
\]

Writing the output values as `(a,b,c)=(P,Q,R)`, the third coordinate can be
rewritten as

\[
c=2s-bs^2+2as^3.
\]

Equivalently, every point in this chart lying over `(a,b,c)` produces a root
of the cubic

\[
\Phi_{a,b,c}(s)=2as^3-bs^2+2s-c.
\]

Define

\[
D=1-bs+3as^2.
\]

Then

\[
\Phi'_{a,b,c}(s)=2D.
\]

For a root with `D != 0`, the input can be reconstructed rationally:

\[
x=\frac{s}{D},\qquad y=b-3as,
\]

and, using `u=1/D`,

\[
z=aD^3-y^2D(D+3).
\]

This factorization explains both phenomena at once. A generic fiber is
controlled by a cubic, while the derivative of that cubic is tied to the
denominator of the reconstruction. The full polynomial determinant identity
is nevertheless certified directly in Lean and does not rely on excluding the
chart boundary.

## 4. The displayed fiber by hand

For the target `(a,b,c)=(-1/4,0,0)`, the fiber equation is

\[
\Phi(s)=-\frac12s^3+2s
       =-\frac12s(s-2)(s+2).
\]

Its roots are `0`, `2`, and `-2`. Substituting them into the reconstruction
formulas produces exactly the three points in the screenshot. For these roots,
`D` is respectively `1`, `-2`, and `-2`, so all three lie in the valid chart.

This supplies a short human explanation of the collision; the Lean theorem
checks the original coordinates directly.

## 5. Escape at infinity

The curve

\[
(x,y,z)=\left(s,-\frac1s,\frac5{s^2}\right)
\]

satisfies

\[
F\left(s,-\frac1s,\frac5{s^2}\right)=\left(0,\frac2s,0\right).
\]

As `|s|` tends to infinity, the source escapes to infinity while the image
tends to zero. Thus the map is not proper. This is compatible with its
everywhere nonzero local Jacobian: local invertibility does not force global
injectivity without a global condition such as properness.

## 6. Consequences and non-consequences

The example proves `not JC(3)`. If `n > 3`, the map

\[
(v,w)\longmapsto(F(v),w),\qquad
v\in\mathbb C^3,\quad w\in\mathbb C^{n-3},
\]

has the same nonzero constant determinant and the same collision, so
`not JC(n)`.

The implication cannot be run downward. Restricting `z` to a constant and
dropping one output coordinate does not preserve the constant-Jacobian
hypothesis. A three-variable counterexample therefore gives no automatic
plane counterexample.

## 7. The affine-in-one-variable plane obstruction

To test the closest low-complexity analogue, consider

\[
G(x,y)=(A(x)y+B(x),C(x)y+D(x)).
\]

Its determinant is

\[
\det JG=(A'C-AC')y+(B'C-AD').
\]

Suppose this equals a nonzero constant `k`. Comparing coefficients gives

\[
A'C-AC'=0,\qquad B'C-AD'=k.
\]

The second identity is an explicit Bezout identity for `A` and `C` after
scaling by `k^{-1}`, so `A` and `C` are coprime. The first identity says their
Wronskian vanishes. Mathlib's coprime-Wronskian theorem forces

\[
A'=C'=0.
\]

In characteristic zero, `A` and `C` are therefore constant polynomials. This
already rules out the direct plane imitation with genuinely varying slopes.

The formal proof goes one step further. Put `alpha=A`, `gamma=C`, now regarded
as field constants, and define

\[
E(x)=\gamma B(x)-\alpha D(x).
\]

The determinant identity gives `E'(x)=k`, so `E(x)-kx` is constant. The same
linear combination of two equal output pairs shows `E(x1)=E(x2)`; since
`k != 0`, this forces `x1=x2`. At least one of `alpha,gamma` is nonzero, or the
determinant would vanish, and the corresponding output then forces `y1=y2`.
Lean certifies the resulting injectivity theorem.

The formal result is more general than the complex case: the derivative
statement holds over any field, while characteristic zero is used to turn zero
derivative into constant slopes and complete the injectivity argument.

## 8. Current plane frontier

The result above treats one rigid ansatz, not arbitrary plane polynomial maps.
A useful external baseline is the 2022 analysis by Guccione, Guccione,
Horruitiner, and Valqui: among hypothetical degree pairs with maximum degree
below `125`, all are excluded except `(72,108)` and `(108,72)`. Equivalently, a
plane counterexample must either have one of those exceptional degree pairs or
have maximum degree at least `125`.

That makes the next honest targets structural rather than brute-force:

1. finish an explicit inverse theorem for the affine-in-one-variable class;
2. determine which part of the three-dimensional cubic-fiber mechanism
   fundamentally requires a third variable; and
3. study the residual `(72,108)` Newton-polygon regime without claiming that a
   low-degree search addresses the full plane problem.

## 9. Independent and adversarial verification

The repository uses two implementations:

- Lean constructs the formal Jacobian from `MvPolynomial.pderiv` and proves
  the determinant and collision inside the kernel.
- SymPy independently expands the determinant and evaluates the points using
  exact rationals.

The SymPy tests also check three hostile fixtures:

- replacing the coefficient `4` by `5` makes the determinant nonconstant; and
- replacing the last point's `13/2` by `15/2` breaks the common fiber; and
- duplicating one valid point is rejected as a three-distinct-preimage
  certificate even though all listed images still equal the target.

These tests are not extra mathematical evidence once the Lean theorem is
trusted. Their role is engineering assurance that the transcription and the
independent checker are not vacuous.

## 10. Sources

- [Original announcement by Levent Alpöge](https://x.com/__alpoge__/status/2079028340955197566)
- [Direct consequences of the three-dimensional counterexample](https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html)
- [Increasing the degree of a possible counterexample from 100 to 108](https://arxiv.org/abs/2204.14178)
