# Jacobian Two

An exact, Lean-first audit of the newly announced three-dimensional
counterexample to the Jacobian conjecture, followed by a deliberately separate
attack on the still-open two-dimensional case.

## Verdict on the screenshot

**The algebra in the screenshot is true.** For

\[
\begin{aligned}
P&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
Q&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
R&=2x-3x^2y-x^3z,
\end{aligned}
\]

the map \(F=(P,Q,R):\mathbb C^3\to\mathbb C^3\) satisfies

\[
\det JF=-2
\]

as a polynomial identity. The three distinct points

\[
\left(0,0,-\tfrac14\right),\quad
\left(1,-\tfrac32,\tfrac{13}{2}\right),\quad
\left(-1,\tfrac32,\tfrac{13}{2}\right)
\]

all map exactly to \((-\tfrac14,0,0)\).

This is a counterexample in **dimension three**. Appending identity coordinates
gives counterexamples in every dimension at least three. It does **not** settle
the plane conjecture: `JC(2)` remains open.

The public announcement is extremely recent. Historical priority and the
discovery narrative should therefore be treated as announced rather than as a
settled scholarly record. The finite identities above do not depend on that
narrative: they are independently certified in this repository.

## What is machine checked

[`JacobianTwo/Counterexample.lean`](JacobianTwo/Counterexample.lean) represents
the coordinates as integer-coefficient multivariate polynomials and proves:

- `formalJacobian_det`: the determinant formed from actual
  `MvPolynomial.pderiv` entries is the constant polynomial `-2`;
- `jacobianAt_det`: after evaluation in any commutative ring, the determinant
  is still `-2` at every point;
- `three_distinct_preimages`: the three exact complex points are pairwise
  distinct and have the displayed common image; and
- `complex_not_injective`: the induced map on `ℂ³` is not injective.

There are no `sorry`, `admit`, or custom axioms in these proofs. A separately
implemented typed SymPy checker in [`scripts/verify.py`](scripts/verify.py)
recomputes the determinant and fiber using exact rational arithmetic. Its test
suite deliberately perturbs both a coefficient and a collision point to make
sure transcription errors are detected.

## A genuine result toward the two-dimensional case

The three-dimensional example is affine in `z`, suggesting the plane ansatz

\[
G(x,y)=\bigl(A(x)y+B(x),\ C(x)y+D(x)\bigr).
\]

Its Jacobian determinant is

\[
(A'C-AC')y+(B'C-AD').
\]

[`JacobianTwo/AffineInOneVariable.lean`](JacobianTwo/AffineInOneVariable.lean)
proves that, over a characteristic-zero field, a nonzero constant determinant
forces both \(A\) and \(C\) to be constant. The proof extracts a Bezout identity
from the constant term, obtains coprimality of \(A,C\), and combines it with
the vanishing Wronskian. The same module then proves that the resulting plane
map is injective. Thus the most direct affine-in-one-variable shadow of the
three-dimensional construction cannot reproduce its noninjective mechanism.

This is a narrow positive noncollision result toward `JC(2)`, not a proof of
`JC(2)`. Packaging the reconstruction as an explicit polynomial inverse is a
separate next step.

## Reproduce everything

The Lean toolchain and mathlib revision are pinned. Python dependencies are
locked by `uv.lock`.

```bash
lake build
uv run --frozen python -m scripts.verify
uv run --frozen pytest
uv run --frozen mypy
```

The exact checker prints:

```text
det JF = -2
F(0, 0, -1/4) = (-1/4, 0, 0)
F(1, -3/2, 13/2) = (-1/4, 0, 0)
F(-1, 3/2, 13/2) = (-1/4, 0, 0)
determinant verified: True
collision verified: True
distinct inputs verified: True
```

## Read the argument

- [`docs/audit.md`](docs/audit.md) gives a hand-checkable structural
  derivation, the cubic fiber equation, and the exact claim boundary.
- [`docs/research-log.md`](docs/research-log.md) records completed work,
  negative results, and the next plane target.
- [`SPEC.md`](SPEC.md) is the research specification and claim-status ledger.

## Sources and provenance

- Levent Alpöge's [original public announcement on X][announcement]
- Zihan Zhang's [direct-consequences note][consequences]
- Guccione, Guccione, Horruitiner, and Valqui's
  [degree-bound paper for the plane problem][degree-bound]

The announcement and expository note are provenance sources. The Lean proof
and exact SymPy checker are the evidence for the two algebraic assertions made
by the screenshot.

[announcement]: https://x.com/__alpoge__/status/2079028340955197566
[consequences]: https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html
[degree-bound]: https://arxiv.org/abs/2204.14178

## License

Apache-2.0. See [`LICENSE`](LICENSE).
