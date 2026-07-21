# Jacobian Two: research specification

## Mission

Build a public, Lean-first research repository that does two things without
conflating them:

1. independently certify the polynomial map announced by Levent Alpöge on
   2026-07-20 UTC (2026-07-19 Pacific) as a counterexample to the complex
   Jacobian conjecture in dimension three; and
2. use the structure of that example to pursue a concrete obstruction theorem
   for the still-separate two-dimensional case.

The repository is tracked by Linear issue
[ALOK-792](https://linear.app/aloksingh/issue/ALOK-792/formalize-the-3d-counterexample-and-establish-the-jc2-research).
The Galois/frontier continuation is tracked by
[ALOK-794](https://linear.app/aloksingh/issue/ALOK-794/audit-the-galois-frontier-and-determine-the-3d-maps-nonproper-value).
The arbitrary-degree affine-coordinate continuation is tracked by
[ALOK-795](https://linear.app/aloksingh/issue/ALOK-795/prove-the-arbitrary-degree-affine-coordinate-plane-keller-theorem).
The constant-leading quadratic-coordinate continuation is tracked by
[ALOK-798](https://linear.app/aloksingh/issue/ALOK-798/prove-the-constant-leading-quadratic-coordinate-keller-theorem).

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
| The universal-in-dimension formulation is false | `DERIVED` | the certified dimension-three counterexample suffices |
| Counterexamples exist in every dimension `n >= 3` | `DERIVED` | append identity coordinates |
| Whether the two-dimensional conjecture is false | `OPEN` | no implication from the dimension-three example is known or assumed here |
| The Galois-case theorem settles every plane Keller map | `FALSE` | constant Jacobian gives finite separable, not normal/Galois |
| Generic degree six is the first unresolved plane sheet degree | `SOURCE-CHECKED` | Campbell, Orevkov, Żołądek, and Borisov literature audit |
| The universal cubic discriminant coefficient expression is `-4*Q` | `LEAN-CERTIFIED` | explicit standard formula; semantic SymPy discriminant checked independently |
| The announced map has fiber counts `3/1/0` on the three stated strata | `DERIVED` | simple projective-root/source bijection plus exact symbolic checks |
| Its image is `C^3 \ Gamma` | `DERIVED` | complete fiber stratification |
| Its nonproper-value set is exactly `V(Q)` | `DERIVED` | explicit escaping family and projective compactness argument |
| A Keller map `(P,e(x)y+f(x))` with arbitrary `P` is an automorphism | `LEAN-CERTIFIED` | actual bivariate Jacobian, degree descent, both explicit inverse charts |
| Historical novelty of the full fiber/nonproper theorem | `UNKNOWN` | derived here; same-day sources compared, priority not established |
| Historical priority and the discovery account | `ANNOUNCED` | public posts are new and not a peer-reviewed historical record |

## Milestone 1: kernel-checked screenshot certificate

Implement the announced coordinates as integer-coefficient
`MvPolynomial (Fin 3) ℤ`. Define the formal Jacobian matrix with
`MvPolynomial.pderiv`, prove its determinant equals the constant `-2`, and
base-change the pointwise certificate through `Int.castRingHom`.

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

The completed extension proves that the resulting plane map is injective: a
linear combination of the outputs recovers an affine polynomial in `x` with
nonzero slope, after which either nonzero constant `y`-slope recovers `y`.

The explicit reconstruction is now completed for this ansatz.  The broader
quadratic-in-one-variable class is treated separately in Milestone 4.

## Milestone 3: Galois audit and the exact nonproper-value set

The classical Galois-case theorem has a one-way hypothesis:

> A Keller map whose induced function-field extension is Galois is a
> polynomial automorphism.

It does **not** say that every two-dimensional Keller extension is Galois.
Consequently it does not settle `JC(2)`.  The literature audit must distinguish
the field-extension degree from the ordinary total degree of the coordinate
polynomials, and it must record the accepted result that a hypothetical plane
counterexample has topological (generic-fiber) degree at least six.

For the announced three-dimensional map, write a target as `(a,b,c)` and put

```text
p(T) = c*T^3 - 2*T^2 + b*T - 2*a
r(T) = p'(T) = 3*c*T^2 - 4*T + b
Q(a,b,c) = 27*a^2*c^2 - 18*a*b*c + 16*a + b^3*c - b^2.
```

On the chart `x != 0`, the fiber coordinate `T = y + 1/x` satisfies
`p(T)=0`, `r(T)=2/x`, and the source is reconstructed by

```text
x = 2/r,
y = T - r/2,
z = 5*r^2/4 - 3*T*r/2 - c*r^3/8.
```

The research target is the following complete determination of the Jelonek
set (the nonproper-value set) of the announced map:

> A target `(a,b,c)` is a nonproper value if and only if `Q(a,b,c)=0`.

The two directions must be checked separately.

1. If a bounded target sequence has an unbounded source sequence and `T`
   stays bounded, reconstruction forces `r -> 0`; the limiting polynomial has
   a repeated root, hence `Q=0`.  If `T -> infinity`, eliminate `c` using
   `p(T)=0`; the exact identities

   ```text
   x = T / (T^2 - b*T + 3*a),
   y = b - 3*a/T
   ```

   and an inverse-power expansion for `z` show that the source instead tends
   to the finite point `(0,b,a-4*b^2)`, a contradiction.
2. If `Q=0`, the cubic (or its nonzero quadratic specialization at `c=0`)
   has a repeated finite root `T0`.  Perturb `T0`, define `a` from `p(T)=0`,
   and use reconstruction.  Then `r(T) -> 0`, so `x=2/r` escapes to infinity
   while the targets converge to `(a,b,c)`.

Acceptance gate:

- Lean certifies the cubic, derivative, discriminant, reconstruction, and
  large-`T` cancellation identities with no `sorry`, `admit`, or custom axiom;
- an independent typed SymPy checker verifies those identities and hostile
  perturbations exactly;
- the readable proof handles the `x=0`, bounded-`T`, and unbounded-`T` cases
  explicitly; and
- no historical-priority claim is made without a literature review.

The complete set-theoretic consequence to be recorded is stronger than the
initial target.  Let

```text
Gamma = { (a,b,c) : 3*b*c = 4 and b^2 = 12*a }.
```

Simple projective roots of the homogenized cubic correspond bijectively to
source points.  Therefore the target strata have fiber cardinalities `3`,
`1`, and `0` according as `Q != 0`, `Q = 0` off `Gamma`, or the target lies on
`Gamma`.  In particular the exact image is `C^3 \ Gamma`.  The algebraic
identities are machine checked; the projective root count and compactness
argument remain explicitly identified as derived mathematical prose.

## Milestone 4: a quadratic-in-one-variable plane fragment

For a characteristic-zero field, consider

```text
(x,y) |-> (a(x)*y^2 + b(x)*y + c(x), e(x)*y + f(x)).
```

Equating the coefficients of a nonzero constant Jacobian gives

```text
a'*e - 2*a*e' = 0,
b'*e - b*e' = 2*a*f',
c'*e - b*f' = k.
```

On the chart `e != 0`, the target normal form is

```text
r = eps*y + f(x),
F(x,y) = (lam*r^2 + mu*r + alpha*x + beta, r),
alpha*eps = k.
```

This has a displayed polynomial inverse.  The `e=0` chart must be handled
separately: the equations force `a=0`, `b` to be a nonzero constant, and `f`
to be affine with nonzero slope, which is again triangular and invertible.

Acceptance gate:

- Lean derives the normal form from the coefficient identities;
- Lean proves both inverse identities and the full bijectivity theorem,
  including the `e=0` branch;
- no theorem is described as the full plane Jacobian conjecture; and
- prior mathematical coverage of this bounded-degree class is cited, so the
  result is presented as formalization rather than novelty.

## Milestone 5: arbitrary degree with one affine coordinate

Let `P : (K[X])[X]` be an arbitrary polynomial in `y` with coefficients in
`K[x]`, and let

```text
Q(x,y) = e(x)*y + f(x).
```

Over a characteristic-zero field, the target theorem is:

> If the formal Jacobian of `(P,Q)` is a nonzero constant `k`, then `(P,Q)`
> is a polynomial automorphism.

The theorem must cover both charts and expose the inverse-producing normal
forms.

1. If `e != 0`, induct on the `y`-degree of `P`.  If `a_m(x)` is the leading
   coefficient, the top Jacobian equation is

   ```text
   a_m' * e = m * a_m * e'.
   ```

   Characteristic zero forces `a_m = lambda * e^m`.  Subtracting
   `lambda * Q^m` lowers the `y`-degree and leaves the Jacobian unchanged.
   Iteration gives

   ```text
   P(x,y) = G(Q(x,y)) + alpha*x + beta,
   e = eps,
   alpha*eps = k,
   ```

   with `alpha` and `eps` nonzero.  This has the displayed inverse

   ```text
   x = (u - G(v) - beta) / alpha,
   y = (v - f(x)) / eps.
   ```

2. If `e = 0`, the Jacobian equation is `-P_y*f' = k`.  Since a product of
   polynomials is a nonzero scalar, `P_y` and `f'` are nonzero constants.
   Thus

   ```text
   P(x,y) = eta*y + c(x),
   f(x) = delta*x + gamma,
   -eta*delta = k,
   ```

   which is triangular with an explicit inverse.

Acceptance gate:

- Lean defines the two formal partial derivations on `(K[X])[X]` and proves
  the stated normal forms from an equality of the actual formal Jacobian;
- Lean proves both inverse laws and full bijectivity in both charts;
- all `y`-degrees are allowed, including the constant and zero-polynomial
  edge cases;
- exact executable tests check generated normal forms and reject perturbed
  nonconstant-slope fixtures;
- `lake build`, the exact Python suite, strict type checking, and an axiom/
  unfinished-proof audit all pass; and
- the literature note distinguishes mathematical coverage from any claim of
  historical novelty and states explicitly that this does not solve `JC(2)`.

## Milestone 6: arbitrary degree with a constant-leading quadratic coordinate

Take the first genuinely nonlinear second coordinate after Milestone 5:

```text
P in K[x,y] arbitrary,
Q(x,y) = eps*y^2 + g(x)*y + f(x),
eps in K^*.
```

Put

```text
s = 2*eps*y + g(x),
Delta = g(x)^2 - 4*eps*f(x).
```

The identities

```text
s^2 - 4*eps*Q = Delta,
J(s,Q) = Delta'/2
```

suggest the exact target theorem:

> If `J(P,Q)=k` for `k in K^*`, then there are nonzero scalars `A` and
> `lambda`, a scalar `B`, and `G in K[T]` such that
>
> ```text
> Delta = A*x + B,
> lambda*A/2 = k,
> P = G(Q) + lambda*s.
> ```

Indeed, evaluate the Jacobian on the polynomial critical section
`y=-g/(2*eps)`, where `Q_y=s=0`.  This gives

```text
P_y(x,-g/(2*eps)) * Delta'(x) = 4*eps*k.
```

Both factors are therefore units, so `Delta'=A` is a nonzero scalar.  The
pair `(s,Q)` is then a polynomial coordinate system with inverse

```text
x = (s^2 - 4*eps*Q - B) / A,
y = (s - g(x)) / (2*eps).
```

In these coordinates, the Jacobian equation says that the derivative of `P`
with respect to `s` is the scalar `lambda=2*k/A`, yielding the normal form.
For target coordinates `(u,v)`, the displayed inverse is

```text
s = (u - G(v)) / lambda,
x = (s^2 - 4*eps*v - B) / A,
y = (s - g(x)) / (2*eps).
```

Acceptance gate:

- the result is hostile-audited for signs, characteristic, zero-slope, and
  degree-drop edge cases before implementation;
- Lean starts from the actual formal bivariate Jacobian and proves the affine
  discriminant, coordinate-pair, normal-form, and two-sided-inverse theorems;
- the proof reuses the bundled derivations and arbitrary-degree descent from
  `AffineCoordinate.lean` rather than restating coefficient hypotheses;
- an exact typed symbolic checker validates high-degree examples and a
  perturbed non-Keller fixture;
- the literature boundary is sourced before any novelty language is used;
  absent such evidence, the result is presented only as a formalized
  structural fragment; and
- the repository continues to state that the unrestricted degree-six sheet
  stratum and `JC(2)` remain open.

## Repository shape

```text
JacobianTwo/
  Counterexample.lean       # formal polynomial and collision certificate
  AffineInOneVariable.lean  # JC(2) obstruction theorem
  QuadraticInOneVariable.lean # quadratic-in-y normal form and inverse
  CubicFiber.lean           # cubic, discriminant, reconstruction, infinity algebra
  Basic.lean                # shared definitions
scripts/
  verify.py                 # independent exact symbolic checker
  nonproper.py              # exact cubic-fiber/nonproper algebra
tests/
  test_verify.py            # positive and adversarial fixtures
  test_nonproper.py         # discriminant, strata, and hostile fixtures
docs/
  audit.md                  # readable derivation, provenance, claim boundary
  galois-frontier.md        # Galois theorem and first open sheet degree
  nonproper-set.md          # complete fiber/image/nonproper proof
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
