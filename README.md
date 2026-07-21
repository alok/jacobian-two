# Jacobian Two

An exact, Lean-first audit of the new three-dimensional Jacobian
counterexample, together with a deliberately honest attack on the surviving
plane conjecture.

## Short answer: the plane case is still open

Galois theory proves a conditional theorem: a Keller map whose induced
function-field extension is **already Galois** is a polynomial automorphism.
For a plane map `F=(P,Q)`, a nonzero constant Jacobian makes

\[
  k(P,Q)\subset k(x,y)
\]

finite and separable, but it does not make the extension normal.  That missing
normality hypothesis is the gap in the tempting “Galois theory solves the
plane” argument.

Here “degree” means generic fiber size, or equivalently the degree of this
function-field extension.  Degrees through five are excluded for a
hypothetical plane counterexample; degree six is the first unresolved
small-sheet case in the accepted literature.  The detailed audit and primary
sources are in [`docs/galois-frontier.md`](docs/galois-frontier.md).

The resulting dimensional boundary is now:

\[
  \mathrm{JC}(1)\text{ is true},\qquad
  \mathrm{JC}(2)\text{ remains open},\qquad
  \mathrm{JC}(n)\text{ is false for }n\ge 3.
\]

This repository does **not** claim to solve `JC(2)` or its degree-six frontier.
It records exactly what was proved, what was derived, and what remains open.

## New progress at the first open sheet degree

The degree-six frontier can now be narrowed substantially.  Let `G <= S6` be
the monodromy group of a hypothetical six-sheet plane Keller counterexample.
Affine local inverse branches force every branch meridian to fix a sheet, and
those meridians normally generate `G`.  Exact enumeration of all sixteen
transitive degree-six groups initially leaves seven possibilities.  Orevkov's
exact defect budget, the finite-normalization boundary, deck symmetry, purity
of the branch locus, and a local double-transposition obstruction eliminate
the three imprimitive survivors.  A new Riemann--Hurwitz refinement then
charges a dicritical of normal index `e` and tangential degree `d` at least
`e*d` units in Orevkov's budget:

\[
  \sum_E(e_Ed_E+\delta_E)=5,\qquad \delta_E\ge0.
\]

Together with an irreducible-branch normalization obstruction, this also
eliminates `A5` and `S5`.  Therefore

\[
  \boxed{G\in\{6T15=A_6,\;6T16=S_6\}}.
\]

Equivalently, the function-field extension of any hypothetical six-sheet
counterexample must be primitive.  This applies without a one-dicritical
hypothesis.  More precisely, the unrestricted ramified profiles are now:

- `A6`: one `(e,d)=(3,1)` branch; its normalization is noninjective, and every
  normalization collision is two smooth multiplicity-three branches which
  use all six sheets and give an omitted affine value;
- `S6`: one transposition branch, two distinct transposition branches, or
  distinct transposition and 3-cycle branches.

Under the additional assumption of exactly one dicritical component, global
branch-curve topology gives the same two groups and narrows the local types:

- only `6T15=A6` with `(e,d)=(3,1)` or `6T16=S6` with `(2,1)` remains;
- the branch normalization `A1 -> B` must be noninjective, so the nonproperness
  curve has a finite multibranch singularity; and
- in the `A6` passport, every normalization collision consists of two smooth
  multiplicity-three branches, uses all six sheets, and is omitted by the
  original affine map; moreover, the two jump units concentrate at a unique
  local-degree-five point, whose branch is a `(2,5)` cusp when its source is
  smooth.

The eliminated one-dicritical types `(2,2)` and `(4,1)` would have injective
normalization.  Lin--Zaidenberg then makes the branch a monomial contractible
curve whose weighted-orbit product identifies the homotopy types of its local
and global complements; its intransitive local six-sheet action cannot equal
the transitive global monodromy.

These are necessary conditions, not constructions or an exclusion of `A6`
and `S6`.  The refined identity and universal elimination are proved in the
[refined budget note](docs/refined-six-sheet-budget.md), with the broader
geometric setup in the
[six-sheet monodromy note](docs/six-sheet-monodromy.md).  The dependency-free
[Python certificate](scripts/six_sheet_monodromy.py) rebuilds the exact groups,
classes, normal closures, normalizers, deck groups, blocks, and local subgroup
orbits; an optional [Sage/GAP checker](tools/check_six_sheet_gap.sage) verifies
the catalogue independently.

## A separate sparse obstruction at coordinate degree `(72,108)`

Generic sheet degree and coordinate degree are different invariants.  On the
coordinate-degree side, Guccione--Guccione--Horruitiner--Valqui reduce the
remaining sub-`125` problem to `(72,108)` and its transpose, with two explicit
transformed Newton-polygon configurations satisfying `[P,Q]=x^2`.

An exhaustive exact support calculation now proves that the first
configuration needs at least three nonzero coefficients strictly inside its
two Newton polygons, while the second needs at least four.  All boundary
lattice coefficients remain arbitrary; only the exact polygon vertices are
assumed nonzero.  The checker certifies all `7504` first-case supports with at
most two interior terms and all `3683` second-case supports with at most three.
Of the latter, `3678` have replayed zero-product certificates and the five
remaining triples have exact unit-ideal certificates.  Hostile fixtures
confirm that the method stops on the full polygons and on a named four-term
second-case support.  This is a sparse-support lower bound, not an elimination
of `(72,108)`.  See the [Newton-polygon note](docs/newton-72-108-sparse.md) and
its [exact certificate](scripts/newton_72_108.py).

## The screenshot is exactly correct

For

\[
\begin{aligned}
A&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
B&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
C&=2x-3x^2y-x^3z,
\end{aligned}
\]

the map `F=(A,B,C)` satisfies the polynomial identity

\[
  \det JF=-2.
\]

The three distinct rational points

\[
\left(0,0,-\tfrac14\right),\quad
\left(1,-\tfrac32,\tfrac{13}{2}\right),\quad
\left(-1,\tfrac32,\tfrac{13}{2}\right)
\]

all map exactly to `(-1/4,0,0)`.  This is a complete finite certificate of a
counterexample in dimension three.  Appending identity coordinates gives the
same conclusion in every higher dimension.

[`JacobianTwo/Counterexample.lean`](JacobianTwo/Counterexample.lean) constructs
the formal Jacobian from actual `MvPolynomial.pderiv` entries, proves its
determinant, proves the three-point collision and pairwise distinctness, and
derives noninjectivity over `ℂ`.  The independent typed SymPy checker
[`scripts/verify.py`](scripts/verify.py) recomputes both identities using exact
rational arithmetic and includes hostile transcription fixtures.

## A solved follow-on problem: every fiber and every escaping value

For a target `(a,b,c)`, introduce the reciprocal fiber coordinate

\[
  T=y+\frac1x
\]

on `x != 0`.  It satisfies

\[
  p(T)=cT^3-2T^2+bT-2a=0,
  \qquad
  p'(T)=\frac2x.
\]

Define

\[
  Q(a,b,c)=27a^2c^2-18abc+16a+b^3c-b^2
\]

and

\[
  \Gamma=\{3bc=4,\ b^2=12a\}.
\]

The exact fiber calculation gives

| Target stratum | Number of source points |
|---|---:|
| `Q != 0` | 3 |
| `Q = 0` and target not in `Gamma` | 1 |
| target in `Gamma` | 0 |

Consequently,

\[
  F(\mathbb C^3)=\mathbb C^3\setminus\Gamma.
\]

Moreover, the complete nonproper-value set—the targets approached by images
of source sequences escaping to infinity—is

\[
  S_F=V(Q).
\]

The proof includes both directions.  Every point of `V(Q)` has an explicit
escaping family, while projective-root compactness plus exact reconstruction
shows that no point outside `V(Q)` can be an asymptotic value.  See
[`docs/nonproper-set.md`](docs/nonproper-set.md) for the complete argument.

[`JacobianTwo/CubicFiber.lean`](JacobianTwo/CubicFiber.lean) certifies the
fiber cubic, its derivative, the standard universal cubic discriminant
coefficient expression `-4Q`, an explicit Bézout common-root certificate,
finite-root reconstruction, and all large-`T` cancellation identities.
[`scripts/nonproper.py`](scripts/nonproper.py)
independently checks the remaining exact algebra: the infinity chart,
repeated- and triple-root parameterizations, singular-locus elimination, and
the escaping family.  The compactness argument is written explicitly in the
mathematical note rather than being mislabeled as kernel-checked topology.

These consequences are labeled **derived here; historical priority not
established**.  Same-day sources already contained the cubic, reconstruction,
discriminant, and generic `S_3` calculation; this repository makes no
literature-priority claim for the fuller stratification.

## Certified positive fragments in the plane

The strongest plane result in this repository has no degree bound on its
first coordinate.  Let

\[
  F=(P,Q),\qquad P\in K[x,y],\qquad Q=e(x)y+f(x),
\]

over a characteristic-zero field.  If the actual formal Jacobian is a
nonzero scalar, Lean proves that `F` is a polynomial automorphism.  In the
`e != 0` chart it derives

\[
  e=\varepsilon\in K^\times,\qquad
  P=G(Q)+\alpha x+\beta,\qquad
  \alpha\varepsilon=k,
\]

and in the `e=0` chart it derives the complementary triangular form.  Both
charts have kernel-checked explicit inverses.  See
[`JacobianTwo/AffineCoordinate.lean`](JacobianTwo/AffineCoordinate.lean) and
the [proof and literature note](docs/affine-coordinate.md).

This is a characteristic-zero algebraic formalization of the known
type-`(m,1)` reduction, not a new mathematical class.  Sabatini's published
real theorem uses the same leading-power elimination.  The repository's
field-uniform statement deliberately assumes a genuinely constant Jacobian.

The next certified class allows a genuinely quadratic coordinate.  Let

\[
  Q=\varepsilon y^2+g(x)y+f(x),\qquad \varepsilon\in K^\times,
  \qquad s=2\varepsilon y+g(x).
\]

For arbitrary `P`, a nonzero constant identity `J(P,Q)=k` forces the
discriminant `Delta=g^2-4*epsilon*f` to be affine, say `Delta=A*x+B` with
`A != 0`, and Lean proves the original-coordinate normal form

\[
  P=G(Q)+\lambda s,\qquad \lambda A/2=k.
\]

It also proves both laws for the explicit polynomial inverse

\[
  \sigma=(u-G(v))/\lambda,\quad
  x=(\sigma^2-4\varepsilon v-B)/A,\quad
  y=(\sigma-g(x))/(2\varepsilon).
\]

See
[`JacobianTwo/ConstantLeadingQuadratic.lean`](JacobianTwo/ConstantLeadingQuadratic.lean)
and the [proof, certificate map, and literature boundary](docs/constant-leading-quadratic.md).
This theorem has no degree bound on `P`; it does require the quadratic leading
coefficient of `Q` to be a nonzero scalar.

The scalar-leading hypothesis is not needed in the known mathematical theorem.
If

\[
  Q=a(x)y^2+g(x)y+f(x),\qquad a\ne0,
\]

and `J(P,Q)` is a nonzero scalar, Moskowicz's Theorem 2.7 already implies that
`(P,Q)` is an automorphism: its invariant
`gcd(2, deg_x(a))` is either `1` or the prime `2`.  Simon--Weimann's coordinate
criterion then implies, after scalar extension if necessary, that `a` is a
nonzero scalar and that `deg_x(g^2-4*a*f)=1`.

[`JacobianTwo/VariableLeadingQuadratic.lean`](JacobianTwo/VariableLeadingQuadratic.lean)
develops an independent direct certificate for that known theorem.  Lean now
certifies the top-coefficient equation, target-shear descent to an odd
`y`-degree, the identity `p_n^2=c*a^n`, the UFD shape
`a=epsilon*h^2, p_n=lambda*h^n`, and the unique even--odd decomposition over a
field.  The fraction-field layer also certifies the quotient-rule derivation
on `K(x)`, its constant field, the affine substitution `y=(U-rho)/h`, the
exact Jacobian factor `k/h`, parity extraction, and the full coefficient
recurrence.  Lean also constructs its primitive explicitly, proves every
coefficient lies in `K[F]`, and tracks the exact degree and nonzero leading
coefficient through the downward descent.  Specialization at `y=0`, a
valuation-free gcd normalization, and a unique-survivor denominator theorem
then force `h | g`; the terminal recurrence forces `h` to be a unit.  The
resulting theorem `variableLeadingQuadratic_bijective_full` is a complete
kernel-checked proof for every coordinate of the displayed at-most-quadratic
form, including its affine branch.  The complete direct proof, hostile
fixtures, and exact literature boundary are in
[`docs/variable-leading-quadratic.md`](docs/variable-leading-quadratic.md).

Two earlier bounded modules expose useful intermediate mechanisms.  First,

\[
  (x,y)\longmapsto(A(x)y+B(x),\ C(x)y+D(x))
\]

cannot be a noninjective Keller map.  A nonzero constant determinant forces
the variable slopes `A` and `C` to be constant, after which a linear
combination of the outputs recovers `x` and then `y`.
[`JacobianTwo/AffineInOneVariable.lean`](JacobianTwo/AffineInOneVariable.lean)
contains the proof.

Second,

\[
  (x,y)\longmapsto
  (a(x)y^2+b(x)y+c(x),\ e(x)y+f(x))
\]

is reduced by its constant-Jacobian coefficient equations to an explicit
triangular normal form in the nonzero-`e` chart, with a displayed polynomial
inverse.  The complementary `e=0` chart is handled separately in the final
theorem.  See
[`JacobianTwo/QuadraticInOneVariable.lean`](JacobianTwo/QuadraticInOneVariable.lean).
These bounded results are now subsumed by the arbitrary-degree theorem, but
their shorter coefficient proofs remain useful.  None of these statements is
a proof of general `JC(2)` or of the generic-degree-six frontier.

## Reproduce the certificates

The Lean toolchain and mathlib revision are pinned.  Python dependencies are
locked by `uv.lock`.

```bash
lake build
uv run --frozen python -m scripts.verify
uv run --frozen python -m scripts.nonproper
uv run --frozen python -m scripts.affine_coordinate
uv run --frozen python -m scripts.constant_leading_quadratic
uv run --frozen python -m scripts.variable_leading_quadratic --depth 9
uv run --frozen python -m scripts.six_sheet_monodromy
uv run --frozen python -m scripts.newton_72_108
uv run --frozen pytest
uv run --frozen mypy
# Optional independent finite-group cross-check:
sage tools/check_six_sheet_gap.sage
```

The Lean source contains no `sorry`, `admit`, or custom axiom.  CI runs the
Lean build, every displayed `uv` command, and the unfinished-proof check.  The
optional Sage/GAP replay is an additional independent local cross-check.

## Reading map

- [`SPEC.md`](SPEC.md) is the research specification and claim-status ledger.
- [`docs/galois-frontier.md`](docs/galois-frontier.md) explains the Galois
  misconception and identifies generic degree six as the first open frontier.
- [`docs/six-sheet-monodromy.md`](docs/six-sheet-monodromy.md) proves the
  degree-six monodromy filters and the conditional one-dicritical passports.
- [`docs/refined-six-sheet-budget.md`](docs/refined-six-sheet-budget.md)
  proves `sum(e*d+delta)=5`, eliminates `A5` and `S5` universally, and lists
  the exact surviving `A6`/`S6` ramified profiles.
- [`docs/a6-one-dicritical-local.md`](docs/a6-one-dicritical-local.md)
  rules out the `A6` jump partition `1+1`, proves the unique
  local-degree-five point, and classifies its smooth-source branch as a
  `(2,5)` cusp.
- [`docs/newton-72-108-sparse.md`](docs/newton-72-108-sparse.md) gives the
  exact sparse-support obstruction in the separate residual coordinate-degree
  configurations.
- [`docs/nonproper-set.md`](docs/nonproper-set.md) proves the complete fiber,
  image, and nonproper-set theorem.
- [`docs/affine-coordinate.md`](docs/affine-coordinate.md) proves the
  arbitrary-degree affine-coordinate normal form and marks its exact
  literature boundary.
- [`docs/constant-leading-quadratic.md`](docs/constant-leading-quadratic.md)
  proves the constant-leading quadratic-coordinate normal form and displays
  its polynomial inverse.
- [`docs/variable-leading-quadratic.md`](docs/variable-leading-quadratic.md)
  gives the known arbitrary-leading quadratic theorem, an independently
  derived direct proof, and its complete Lean certificate.
- [`docs/audit.md`](docs/audit.md) gives a hand-checkable structural derivation
  of the original screenshot.
- [`docs/research-log.md`](docs/research-log.md) records completed work,
  negative results, and remaining obligations.

## Sources and provenance

- Levent Alpöge's [original public announcement on X][announcement]
- L. Andrew Campbell's [Galois-case theorem][campbell]
- S. Yu. Orevkov's [three-sheet theorem][orevkov]
- A. V. Domrina's [four-sheet theorem][domrina]
- Henryk Żołądek's [result through generic degree five][zoladek]
- Alexander Borisov's [Keller-map compactification framework][borisov]
- Guccione--Guccione--Horruitiner--Valqui's [(72,108) reduction][guccione]
- Vered Moskowicz's [quadratic-coordinate antecedent][moskowicz]
- Denis Simon and Martin Weimann's [coordinate/discriminant criterion][simon-weimann]
- Marco Sabatini's [type-`(m,1)` triangular reduction][sabatini]
- Zihan Zhang's [direct-consequences note][consequences]

The announcement and expository note establish provenance.  The finite
claims made here are supported by the repository's reproducible Lean and exact
symbolic certificates.

[announcement]: https://x.com/__alpoge__/status/2079028340955197566
[campbell]: https://doi.org/10.1007/BF01349234
[orevkov]: https://doi.org/10.1070/IM1987v029n03ABEH000984
[domrina]: https://doi.org/10.1070/im2000v064n01ABEH000273
[zoladek]: https://doi.org/10.1016/j.top.2008.04.001
[borisov]: https://www.combinatorics.org/ojs/index.php/eljc/article/view/v27i3p54
[guccione]: https://arxiv.org/abs/2204.14178
[moskowicz]: https://arxiv.org/abs/1810.08202
[simon-weimann]: https://doi.org/10.1216/JCA-2018-10-4-559
[sabatini]: https://doi.org/10.4064/cm9195-1-2024
[consequences]: https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html

## License

Apache-2.0. See [`LICENSE`](LICENSE).
