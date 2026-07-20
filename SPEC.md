# Jacobian Two: research specification

## Mission

Build a public, Lean-first research repository that does two things without
conflating them:

1. independently certify the polynomial map announced by Levent Alpöge on
   2026-07-20 as a counterexample to the complex Jacobian conjecture in
   dimension three; and
2. use the structure of that example to pursue a concrete obstruction theorem
   for the still-separate two-dimensional case.

The repository is tracked by Linear issue
[ALOK-792](https://linear.app/aloksingh/issue/ALOK-792/formalize-the-3d-counterexample-and-establish-the-jc2-research).

## Exact input claim

Let `F = (P, Q, R)` with

```text
P = (1 + x*y)^3*z + y^2*(1 + x*y)*(4 + 3*x*y)
Q = y + 3*x*(1 + x*y)^2*z + 3*x*y^2*(4 + 3*x*y)
R = 2*x - 3*x^2*y - x^3*z.
```

The screenshot claims:

- the formal Jacobian determinant of `F` is the constant polynomial `-2`; and
- the three distinct rational points
  `(0, 0, -1/4)`, `(1, -3/2, 13/2)`, and
  `(-1, 3/2, 13/2)` all map to `(-1/4, 0, 0)`.

These are finite algebraic identities. They must be proved exactly, not inferred
from floating-point samples or from the authority of the announcement.

## Claim-status ledger

| Claim | Required status | Evidence |
|---|---|---|
| The screenshot is transcribed correctly | `SOURCE-CHECKED` | screenshot, original post, and two independent transcriptions agree |
| `det JF = -2` identically | `LEAN-CERTIFIED` | formal partial derivatives of multivariate polynomials and a determinant identity |
| The displayed fiber has three distinct points | `LEAN-CERTIFIED` | exact rational evaluation and pairwise inequality |
| The all-dimensional Jacobian conjecture is false | `DERIVED` | the certified dimension-three counterexample suffices |
| Counterexamples exist in every dimension `n >= 3` | `DERIVED` | append identity coordinates |
| The two-dimensional conjecture is false | `OPEN` | no implication from the dimension-three example is known or assumed here |
| Historical priority and the discovery account | `ANNOUNCED` | public posts are new and not a peer-reviewed historical record |

## Milestone 1: kernel-checked screenshot certificate

Implement the announced coordinates as `MvPolynomial (Fin 3) R` over a
commutative ring. Define the formal Jacobian matrix with
`MvPolynomial.pderiv`, and prove its determinant equals the constant `-2`.

Separately evaluate the map at the three rational points, prove their images
are equal, prove the inputs are pairwise distinct, and derive that the induced
map is not injective. Where useful, transport the rational certificate to
`Complex` explicitly.

Acceptance gate:

- `lake build` passes;
- no `sorry`, `admit`, or custom `axiom` occurs in the certificate;
- `#print axioms` reports only standard foundations, if any;
- an independent typed SymPy program reproduces both identities with exact
  rational arithmetic; and
- adversarial tests fail after perturbing a coefficient or a collision point.

## Milestone 2: a two-dimensional obstruction theorem

The three-dimensional map is affine in `z`. The first tractable question is
whether the same low-complexity strategy can work in two dimensions.

For a characteristic-zero field `K`, consider

```text
G(x, y) = (A(x)*y + B(x), C(x)*y + D(x))
```

with `A, B, C, D : K[X]`. Its Jacobian splits as

```text
det JG = (A'*C - A*C')*y + (B'*C - A*D').
```

The first formal target is:

> If `det JG` is a nonzero constant, then `A' = 0` and `C' = 0`.

Proof plan: the constant term gives an explicit Bézout identity for `A` and
`C`, hence they are coprime. The vanishing `y` coefficient is their vanishing
Wronskian. Mathlib's coprime-Wronskian theorem then forces both derivatives to
vanish. This rules out a direct two-dimensional analogue in which both
coordinates are affine in one source variable.

The next target, after the slope theorem, is to construct an explicit
polynomial inverse for this entire ansatz and thereby package the obstruction
as a genuine positive fragment of `JC(2)`.

## Repository shape

```text
JacobianTwo/
  Counterexample.lean       # formal polynomial and collision certificate
  AffineInOneVariable.lean  # JC(2) obstruction theorem
  Basic.lean                # shared definitions
scripts/
  verify.py                 # independent exact symbolic checker
tests/
  test_verify.py            # positive and adversarial fixtures
docs/
  audit.md                  # readable derivation, provenance, claim boundary
  research-log.md           # dated attempts and open obligations
```

## Reproducibility policy

- Pin Lean and mathlib to the same released version.
- Use `uv` for every Python command and dependency.
- Treat Lean as the authoritative certificate and SymPy as an independent
  implementation cross-check.
- Record negative results and failed reductions; do not silently turn them into
  conjectural prose.
- Do not claim a result about `JC(2)` merely because `JC(3)` is false.

## Source boundary

The originating public record is
[Alpöge's X post](https://x.com/__alpoge__/status/2079028340955197566).
The repository may cite expository follow-ups for context, but its decisive
algebraic claims must remain independently checkable from the formula above.
