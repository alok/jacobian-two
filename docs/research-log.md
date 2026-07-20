# Research log

## 2026-07-20: source audit

Status: complete.

- Transcribed the polynomial map and all four displayed points from the
  screenshot.
- Matched the transcription against Levent Alpöge's public announcement and
  independent same-day expositions.
- Classified provenance claims as `ANNOUNCED`, because the public record is
  only hours old.
- Classified the determinant and collision as finite identities suitable for
  exact independent verification.

## 2026-07-20: exact counterexample certificate

Status: Lean-certified.

- Defined the three coordinate polynomials over `ℤ`.
- Constructed the formal Jacobian using `MvPolynomial.pderiv`.
- Proved its determinant is the constant polynomial `-2`.
- Transported the pointwise determinant statement to arbitrary commutative
  rings.
- Proved the three complex points are pairwise distinct and have the same
  image.
- Derived noninjectivity.
- Audited theorem dependencies with `#print axioms`; only Lean's standard
  foundations occur.

## 2026-07-20: independent exact checker

Status: complete.

- Added a typed SymPy implementation using exact rational arithmetic.
- Added a positive determinant-and-fiber test.
- Added a coefficient perturbation that destroys the constant determinant.
- Added a point perturbation that destroys the common fiber.
- Checked the Python under strict mypy settings.

## 2026-07-20: structural factorization

Status: exact derivation recorded; not yet formalized in Lean.

- Introduced `u=1+xy`, `H=u^2 z+y^2(4+3xy)`, and `s=x/u`.
- Reduced fibers on the `u != 0` chart to
  `2*a*s^3-b*s^2+2*s-c=0`.
- Derived the rational reconstruction formulas and the identity
  `Phi'(s)=2*(1-b*s+3*a*s^2)`.
- Factored the displayed fiber as `-s*(s-2)*(s+2)/2`.
- Exhibited the escaping curve `(s,-1/s,5/s^2)`, whose image is `(0,2/s,0)`.

The chart restriction matters. These formulas illuminate the example but do
not replace the global polynomial determinant proof.

## 2026-07-20: first plane obstruction

Status: Lean-certified.

Tested the closest plane analogue

```text
(x,y) |-> (A(x)*y+B(x), C(x)*y+D(x)).
```

A nonzero constant determinant makes `A` and `C` coprime and gives them a
vanishing Wronskian. Lean proves `A'=C'=0`; in characteristic zero both slopes
are constant. A second Lean theorem combines the outputs to recover `x` and
then `y`, proving that the whole map is injective. This is a meaningful
noncollision result toward `JC(2)`, but it is far from the full conjecture and
does not yet package an explicit polynomial inverse.

## Negative results and guarded boundaries

- Freezing `z` and selecting two outputs does not inherit a constant plane
  Jacobian. There is no direct dimensional descent.
- A low-degree brute-force search would not address the known plane frontier.
- The three-dimensional example does not justify saying that `JC(2)` is false.
- Same-day social posts are not enough to settle historical attribution, even
  though the displayed finite certificate is fully checkable.

## Next research obligations

1. Prove an explicit inverse formula for the whole affine-in-one-variable
   plane class under a nonzero constant determinant.
2. Formalize the cubic-fiber reconstruction and its denominator/derivative
   identity.
3. Connect any broader plane ansatz to the known degree-pair restrictions,
   especially the residual `(72,108)` and `(108,72)` cases.
4. Keep every claim tagged as certified, derived, announced, experimental, or
   open.
