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
- Added a duplicate-input fixture that rejects a vacuous three-point claim.
- Checked the Python under strict mypy settings.

## 2026-07-20: structural factorization

Status: algebraic core Lean-certified.

- Introduced `u=1+xy`, `H=u^2 z+y^2(4+3xy)`, and `s=x/u`.
- Reduced fibers on the `u != 0` chart to
  `2*a*s^3-b*s^2+2*s-c=0`.
- Derived the rational reconstruction formulas and the identity
  `Phi'(s)=2*(1-b*s+3*a*s^2)`.
- Factored the displayed fiber as `-s*(s-2)*(s+2)/2`.
- Exhibited the escaping curve `(s,-1/s,5/s^2)`, whose image is `(0,2/s,0)`.
- Added the reciprocal `T=1/s` cubic, derivative, reconstruction, and
  large-`T` cancellation identities to `JacobianTwo/CubicFiber.lean`.

The chart restriction matters. These formulas illuminate the example but do
not replace the global polynomial determinant proof.

## 2026-07-20: Galois misconception and first open sheet degree

Status: source-audited.

- Checked the classical Galois-case theorem against Campbell, Razar, and
  Bass–Connell–Wright.
- Isolated the missing hypothesis: a Keller extension is finite separable but
  need not be normal.
- Kept generic fiber degree separate from ordinary coordinate degree.
- Located the accepted low-sheet boundary: Orevkov excludes degree three and
  Żołądek excludes through degree five, so degree six is the first unresolved
  generic fiber degree.
- Rejected the withdrawn/corrected Bartenwerfer plane proof as a dependency.
- Treated a 2024 prime-degree preprint as unconfirmed rather than silently
  promoting it to an accepted theorem.

## 2026-07-20: complete fibers, image, and nonproper set

Status: exact algebra machine-checked; projective/topological proof derived.

- Proved that the standard universal cubic discriminant coefficient
  expression is `-4Q` and added an explicit integral-coefficient Bézout
  certificate.
- Proved finite simple-root reconstruction and all cancellations needed at
  the projective root at infinity.
- Identified the singular/triple-root curve
  `Gamma={3*b*c=4, b^2=12*a}`.
- Derived fiber counts `3/1/0` on the three target strata and the exact image
  `C^3 \ Gamma`.
- Constructed an exact escaping family through every point of `V(Q)`.
- Proved conversely, by projective-root compactness and reconstruction, that
  targets outside `V(Q)` are not asymptotic values.
- Concluded `S_F=V(Q)`.

The full fiber/nonproper theorem is labeled “derived here; historical priority
not established.” Same-day comparison sources already contain the cubic,
reconstruction, discriminant, and generic `S_3` computation.

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
noncollision result toward `JC(2)`, but it is far from the full conjecture.
The later quadratic-in-one-variable module contains this affine class and
packages it inside an explicit two-sided-inverse theorem.

## 2026-07-20: quadratic-in-one-variable plane fragment

Status: Lean-certified formalization of a known positive class.

For maps

```text
(x,y) |-> (a(x)*y^2+b(x)*y+c(x), e(x)*y+f(x)),
```

the three constant-Jacobian coefficient equations force a triangular normal
form.  On `e != 0`, writing `r=eps*y+f(x)` gives

```text
F(x,y) = (lam*r^2+mu*r+alpha*x+beta, r),
alpha*eps = k.
```

The Lean module supplies the displayed inverse and proves both inverse laws.
The complementary `e=0` case is triangular after the equations force `a=0`,
`b` constant nonzero, and `f` affine with nonzero slope.  The final theorem
combines the two branches.  This narrows a known class; it does not reach the
generic-degree-six frontier.

## 2026-07-21: arbitrary-degree affine-coordinate theorem

Status: Lean-certified formalization of a known positive class.

For an arbitrary first coordinate `P in K[x,y]` and

```text
Q(x,y) = e(x)*y + f(x),
```

Lean now defines both bivariate partial derivatives as bundled derivations and
starts from the actual formal identity `J(P,Q)=k`.  Strong induction removes
the top `y`-power of `P` by subtracting a scalar multiple of a power of `Q`.
The final theorem covers both charts:

```text
e != 0:  e=eps, P=G(Q)+alpha*x+beta, alpha*eps=k;
e = 0:   P=eta*y+c(x), Q=delta*x+gamma, -eta*delta=k.
```

Both normal forms have explicit kernel-checked two-sided inverses.  An exact
SymPy checker exercises a degree-seven outer polynomial and rejects a
nonconstant-slope perturbation.  The type-`(m,1)` elimination is known from
Sabatini's 2024 theorem; the repository contribution is the
characteristic-zero field-uniform formalization and explicit packaging, not a
novelty claim.

## 2026-07-21: constant-leading quadratic-coordinate theorem

Status: Lean-certified structural positive class; no novelty claim.

For arbitrary `P in K[x,y]` and

```text
Q = eps*y^2 + g(x)*y + f(x),  eps != 0,
s = 2*eps*y + g(x),
Delta = g(x)^2 - 4*eps*f(x),
```

the actual identity `J(P,Q)=k != 0` forces `Delta=A*x+B` with `A != 0`.
The completed-square pair `(s,Q)` is then a polynomial coordinate system.
After transporting the Jacobian equation through this coordinate change, Lean
proves

```text
P = G(Q) + lambda*s,  lambda*A/2 = k,
```

and verifies both laws of the resulting explicit inverse.  The independent
exact checker uses a degree-six `g`, a degree-seven `G`, and a first coordinate
of total degree 84.  Changing one coefficient of `f` makes the determinant
`15-36*x`, so the hostile fixture confirms that the affine-discriminant step is
doing real work.  This reaches arbitrary degree in the other coordinate but
still assumes that the coefficient of `y^2` is a nonzero scalar.

## 2026-07-21: variable-leading quadratic checkpoint

Status: known theorem; independently derived direct proof; partially
Lean-certified.

For

```text
Q = a(x)*y^2 + g(x)*y + f(x),  a != 0,
J(P,Q) = k != 0,
```

Moskowicz's Theorem 2.7 already proves that `(P,Q)` is an automorphism.  Its
invariant is `gcd(2,deg_x(a))`, which is either `1` or the prime `2`, so the
theorem matches an arbitrary characteristic-zero field and places no degree
bound on `P`.  Simon--Weimann's coordinate criterion then gives the stronger
corollaries that `a` is a nonzero scalar and `deg_x(g^2-4*a*f)=1`.  The exact
theorem statement is therefore not novel.

The repository nevertheless has a useful new certificate route.  From the
actual bivariate Jacobian, Lean now proves the top-coefficient differential
equation, all even target-shear steps needed to reach odd degree, the identity
`p_n^2=c*a^n`, and the UFD consequence

```text
a=epsilon*h^2,  p_n=lambda*h^n.
```

Over any field, Lean also proves existence and uniqueness of

```text
R(U) = H(epsilon*U^2+F) + U*B(epsilon*U^2+F),
```

together with the exact degree and leading coefficient of `B` when `R` has
odd degree.  The remaining formal boundary is to instantiate this over
`K(x)`, transport the `x`-derivation at fixed `Q`, extract the recurrence, and
formalize the denominator/valuation proof that first forces `h | g` and then
forces `h` to be a unit.  Exact hostile fixtures distinguish rational from
polynomial mates and exhibit the nonterminating odd coefficient tail for
`Q=x+x^2*y^2`.

## 2026-07-21: fraction-field quadratic bridge

Status: known theorem; new direct-proof infrastructure is Lean-certified;
denominator integration remains.

The `K(x)` boundary above is now crossed.  The repository constructs the
quotient-rule derivation on `RatFunc K`, proves that its kernel is exactly the
constant field in characteristic zero, and transports it to the generic
`FractionRing K[X]` used by the affine chart.  The substitution

```text
y = (U-rho)/h
```

is certified coefficientwise.  A general composition theorem proves that the
centered Jacobian is multiplied by `1/h`; consequently an original Keller
equation `J(P,Q)=k` becomes exactly `k/h`.  The parity identity and uniqueness
then separate the odd part and produce both `D(H)=0` and the complete
coefficient recurrence for `B`.  Independent denominator lemmas certify the
unique-survivor obstruction and the implication “a unit numerator divided by
`h` lies in the base ring, therefore `h` is a unit.”

The remaining direct-proof gap is narrower: formalize the exact solution of
the recurrence in `K[F]`, use the specialization of the original polynomial
mate at `y=0` to supply the unique-survivor hypotheses forcing `h | g`, and
then prove that the terminal `k/h` equation lies in `K[x]`.  No automorphy
claim is promoted until those links are kernel-checked.

## Negative results and guarded boundaries

- Freezing `z` and selecting two outputs does not inherit a constant plane
  Jacobian. There is no direct dimensional descent.
- A low-degree brute-force search would not address the known plane frontier.
- The three-dimensional example does not justify saying that `JC(2)` is false.
- Finite separability does not justify calling the plane function-field
  extension Galois; normality is a separate condition.
- The degree-six frontier has been identified, not solved.
- Same-day social posts are not enough to settle historical attribution, even
  though the displayed finite certificate is fully checkable.

## Next research obligations

1. Complete the direct variable-leading certificate.  The coefficient descent,
   odd leading-square identity, UFD square shape, fraction-field transport,
   `k/h` Jacobian equation, and abstract recurrence are Lean-certified.  The
   remaining work is the recurrence solution and the two polynomiality/
   denominator links forcing `h | g` and then `h` to be a unit.
2. Formalize more of the projective simple-root/fiber correspondence if a
   useful reusable algebraic-geometry interface is available in mathlib.
3. Connect broader plane ansätze to the known degree-pair restrictions,
   especially the residual `(72,108)` and `(108,72)` cases.
4. Keep every claim tagged as certified, derived, announced, experimental, or
   open.
