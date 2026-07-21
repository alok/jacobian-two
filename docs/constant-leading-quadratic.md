# A constant-leading quadratic coordinate

Status: **Lean-certified.** The module
[`JacobianTwo/ConstantLeadingQuadratic.lean`](../JacobianTwo/ConstantLeadingQuadratic.lean)
starts from the actual formal bivariate Jacobian, proves the original-coordinate
normal form below, and kernel-checks both laws of the displayed inverse.  This
is a structural positive class for the plane Jacobian problem, not a solution
of `JC(2)` and not a claim of mathematical novelty.

## Statement and conventions

Let `K` be a field of characteristic zero. For `P in K[x,y]`, take

\[
  Q(x,y)=\varepsilon y^2+g(x)y+f(x),
  \qquad \varepsilon\in K^\times,
  \qquad f,g\in K[x],
\]

and use the sign convention

\[
  J(P,Q)=P_xQ_y-P_yQ_x.
\]

Define the completed-square coordinate and the `y`-discriminant by

\[
  s=2\varepsilon y+g(x),
  \qquad
  \Delta=g(x)^2-4\varepsilon f(x).
\]

Milestone 6 proves the following exact normal-form theorem.

> If `J(P,Q)=k` for some `k in K^times`, then there are uniquely determined
> `A in K^times`, `B in K`, `lambda in K^times`, and `G in K[T]` such that
>
> \[
>   \Delta=Ax+B,
>   \qquad
>   \frac{\lambda A}{2}=k,
>   \qquad
>   P=G(Q)+\lambda s.
> \]
>
> Consequently `(P,Q)` is a polynomial automorphism.

The uniqueness assertion is harmless but useful for hostile checking: `A`
and `B` are the coefficients of `Delta`, `lambda=2k/A`, and `G` is unique
because `(s,Q)` is a polynomial coordinate system.

For target coordinates `(u,v)`, put

\[
  \sigma=\frac{u-G(v)}{\lambda},
  \qquad
  \xi=\frac{\sigma^2-4\varepsilon v-B}{A},
  \qquad
  \eta=\frac{\sigma-g(\xi)}{2\varepsilon}.
\]

Then the polynomial inverse is

\[
  (u,v)\longmapsto(\xi,\eta).
\]

All displayed divisions are by nonzero elements of the ground field, so this
really is a polynomial map over `K`.

## 1. The two discriminant identities

Direct expansion gives

\[
\begin{aligned}
  s^2-4\varepsilon Q
    &=(2\varepsilon y+g)^2
      -4\varepsilon(\varepsilon y^2+gy+f) \\
    &=g^2-4\varepsilon f
     =\Delta. \tag{1}
\end{aligned}
\]

Also

\[
  s_x=g',\quad s_y=2\varepsilon,
  \quad Q_x=g'y+f',\quad Q_y=2\varepsilon y+g=s.
\]

Therefore

\[
\begin{aligned}
  J(s,Q)
    &=g's-2\varepsilon(g'y+f') \\
    &=gg'-2\varepsilon f'
     =\frac{\Delta'}{2}. \tag{2}
\end{aligned}
\]

The factors and signs in (1) and (2) are essential. In particular, with the
opposite Jacobian convention the last expression would be `-Delta'/2`.

## 2. The critical section forces an affine discriminant

Because `2*epsilon` is a unit, the critical section of `Q` with respect to
`y` is the polynomial section

\[
  y_0(x)=-\frac{g(x)}{2\varepsilon}.
\]

On this section, `Q_y=s=0`, while

\[
  Q_x(x,y_0)
    =f'-\frac{gg'}{2\varepsilon}
    =-\frac{\Delta'}{4\varepsilon}.
\]

Evaluating the Keller identity there gives

\[
\begin{aligned}
  k
    &=J(P,Q)(x,y_0) \\
    &=-P_y(x,y_0)Q_x(x,y_0)
     =\frac{P_y(x,y_0)\Delta'(x)}{4\varepsilon}.
\end{aligned}
\]

Equivalently,

\[
  P_y(x,y_0)\Delta'(x)=4\varepsilon k. \tag{3}
\]

Both factors on the left of (3) lie in `K[x]`, and their product is a
nonzero scalar. Since the only units of `K[x]` are the nonzero scalars,
`Delta'` is a nonzero scalar; call it `A`. Characteristic zero then gives

\[
  \Delta=Ax+B
\]

for a scalar `B`, and necessarily `A != 0`. This step is stronger than a
generic-degree argument: it excludes a constant or nonlinear discriminant
before any factorization or algebraic closure is used.

## 3. `(s,Q)` is an explicit polynomial coordinate pair

Equation (1) and `Delta=Ax+B` reconstruct the source variables from `(s,Q)`:

\[
  x=\frac{s^2-4\varepsilon Q-B}{A},
  \qquad
  y=\frac{s-g(x)}{2\varepsilon}. \tag{4}
\]

Thus

\[
  C:(x,y)\longmapsto(s,Q)
\]

is a polynomial automorphism, with (4) as its inverse. Its Jacobian is
`A/2` by (2), in agreement with the inverse formula.

This is worth checking in both directions rather than invoking an inverse
function heuristic. Starting from `(x,y)`, equation (1) makes the reconstructed
`x` equal `(Delta-B)/A=x`, after which the reconstructed `y` is immediate.
Starting from independent variables `(sigma,v)`, define `xi` and `eta` by
(4). Then `s(xi,eta)=sigma`, and

\[
  A\xi+B=\sigma^2-4\varepsilon v.
\]

Applying (1) at `(xi,eta)` gives `Q(xi,eta)=v`.

## 4. Constant Jacobian in the new coordinates

Write

\[
  H(s,t)=P\bigl(C^{-1}(s,t)\bigr),
  \qquad t=Q.
\]

The polynomial chain rule and (2) give

\[
  k=J_{x,y}(P,Q)
    =J_{s,t}(H,t)J_{x,y}(s,Q)
    =H_s\frac{A}{2}.
\]

Hence

\[
  H_s=\lambda:=\frac{2k}{A}. \tag{5}
\]

View `H` as an element of `K[t][s]`. In characteristic zero, the kernel of
formal differentiation with respect to `s` is exactly `K[t]`. Integrating
(5) algebraically therefore yields a unique `G in K[t]` with

\[
  H(s,t)=\lambda s+G(t).
\]

Returning to `(x,y)` proves

\[
  P=G(Q)+\lambda s,
  \qquad
  \frac{\lambda A}{2}=k.
\]

Since `A` and `k` are nonzero, so is `lambda`.

## 5. Direct verification of the displayed inverse

Let `(u,v)=(P(x,y),Q(x,y))`. The normal form first recovers

\[
  \sigma=\frac{u-G(v)}{\lambda}=s.
\]

Then (1) and `Delta=Ax+B` give

\[
  \xi
    =\frac{s^2-4\varepsilon Q-B}{A}
    =\frac{\Delta-B}{A}
    =x,
\]

and finally `eta=(s-g(x))/(2*epsilon)=y`. This proves the left-inverse
identity.

Conversely, begin with arbitrary `(u,v)` and define `(sigma,xi,eta)` as in
the statement. Section 3 gives

\[
  s(\xi,\eta)=\sigma,
  \qquad
  Q(\xi,\eta)=v.
\]

The normal form then gives

\[
  P(\xi,\eta)=G(v)+\lambda\sigma=u.
\]

This is the right-inverse identity. No appeal to injectivity, surjectivity,
or point counting is hidden in the argument.

## Hostile edge-case audit

- **Characteristic zero is used twice.** It turns `Delta'=A` into
  `Delta=Ax+B`, and it identifies the derivative kernel in `K[t][s]` with
  `K[t]`. The theorem is false in positive characteristic. For example, over
  a field of odd characteristic `p`, take

  \[
    Q=x+y^2,\qquad s=2y,\qquad P=s-s^p.
  \]

  Then `J(P,Q)=-2`, but in the polynomial coordinates `(s,Q)` the map is
  `(s,Q) -> (s-s^p,Q)`, which has collisions after base change to an
  algebraic closure and is not an automorphism. The missing `s^p` term is
  invisible to formal differentiation.

- **Characteristic two is not a removable notational nuisance.** The
  critical section and the coordinate formulas divide by `2*epsilon` and
  `4*epsilon`. It belongs to a different problem.

- **`epsilon != 0` is essential to this chart.** If the quadratic leading
  coefficient vanishes, `Q` is affine in `y` and Milestone 5 applies instead.
  A nonconstant leading coefficient `a(x)` in `a(x)y^2+g(x)y+f(x)` is not
  covered here.

- **`k != 0` is essential.** If `k=0`, equation (3) does not force either
  factor to be a unit; for example `P=Q` has zero Jacobian with essentially
  arbitrary discriminant.

- **There is no `A=0` Keller branch.** Under the stated hypotheses, (3)
  makes `Delta'=0` incompatible with `epsilon*k != 0`.

- **A `y`-independent `P` is automatically excluded.** In that case the
  left factor in (3) is zero, contradicting `k != 0`. No degree-zero base
  case may be silently admitted in a formal proof.

- **`g=0` is allowed.** Then `s=2*epsilon*y` and the same proof works; the
  discriminant condition simply says `-4*epsilon*f=Ax+B`.

- **The degree of `P` in `y` is genuinely unrestricted, but constrained by
  the conclusion.** If `G` is nonconstant, `deg_y(P)=2*deg(G)`; if `G` is
  constant, `deg_y(P)=1`. A hostile test should reject a purported Keller
  example of any other positive `y`-degree.

- **No algebraic closure is needed for the proof.** All unit, derivative,
  and inverse arguments take place over the original characteristic-zero
  field. This is broader than the hypotheses of the discriminant theorem
  cited below.

## Lean certificate and reproduction

The Lean module starts from an equality of the actual formal bivariate
Jacobian, not from equation (3) as an extra premise.  Its certificate consists
of:

1. definitions of `s`, `Q`, and `Delta` in `(K[X])[Y]`;
2. formal proofs of (1) and (2);
3. a substitution homomorphism for the critical section and a proof of (3);
4. the unit argument and `Delta=AX+B`, with `A != 0`;
5. both compositions of the coordinate inverse (4);
6. transport of `P` through that coordinate automorphism, followed by the
   characteristic-zero derivative-kernel argument;
7. both compositions of the final inverse `(u,v) -> (xi,eta)`; and
8. an axiom and unfinished-proof audit after `lake build`.

The principal declarations are:

- `affine_discriminant`, which obtains `Delta=A*X+B` with `A != 0` from the
  actual Jacobian equality;
- `coordinatePairInverse_left`, `coordinatePairInverse_right`, and
  `coordinatePair_bijective`, which certify the coordinate pair `(s,Q)`;
- `constantLeadingQuadratic_original_normal_form`, which proves
  `P=G(Q)+lambda*s` and `lambda*A/2=k` in the original coordinates;
- `constantLeadingQuadraticInverse_left` and
  `constantLeadingQuadraticInverse_right`, which certify the displayed final
  inverse; and
- `constantLeadingQuadratic_bijective`, the main automorphism theorem.

Reproduce the kernel and exact-symbolic checks with:

```bash
lake env lean -DwarningAsError=true JacobianTwo/ConstantLeadingQuadratic.lean
lake build +JacobianTwo.ConstantLeadingQuadratic
uv run --frozen python -m scripts.constant_leading_quadratic
uv run --frozen pytest tests/test_constant_leading_quadratic.py
uv run --frozen mypy scripts/constant_leading_quadratic.py \
  tests/test_constant_leading_quadratic.py
```

The axiom audit of the normal form, bijectivity theorem, and two inverse laws
reports only Lean's standard `propext`, `Classical.choice`, and `Quot.sound`;
there are no custom axioms or unfinished proofs.

## Literature boundary

There are two directly relevant prior results, and they overlap this argument
in different ways.

Denis Simon and Martin Weimann prove that, over an algebraically closed field
of characteristic zero, an irreducible polynomial that is monic in `y` and
has minimal `y`-discriminant is a coordinate polynomial. Their convention
allows any nonzero constant leading coefficient under the word “monic.” In
the present setting, after scaling by `epsilon`, `Q` is monic of `y`-degree
two. Once equation (3) gives `Delta=Ax+B` with `A != 0`, its discriminant has
`x`-degree one, the minimal value `deg_y(Q)-1`; the nonsquare affine
discriminant also makes the quadratic irreducible. Thus their theorem already
places this `Q` in the known minimal-discriminant coordinate class. See
Simon--Weimann, *Plane Curves With Minimal Discriminant*, Journal of
Commutative Algebra **10** (2018), 559--598,
[DOI 10.1216/JCA-2018-10-4-559](https://doi.org/10.1216/JCA-2018-10-4-559),
and the [author preprint, arXiv:1507.01091](https://arxiv.org/abs/1507.01091),
especially Theorems 1.2 and 3.3.

Marco Sabatini studies real maps whose two components are polynomial in `y`.
His published result says that a real polynomial `y`-quadratic map with
constant nonzero Jacobian is a composition of triangular maps (after the
normalizing translation used in the paper), hence has a polynomial inverse.
That covers the subcase here in which `P` is also `y`-quadratic. It does not,
by itself, state the arbitrary-`y`-degree first-coordinate normal form proved
above, and its ambient setting is real analytic/polynomial rather than an
arbitrary characteristic-zero field. See Sabatini, *Global injectivity of
planar non-singular maps that are polynomial in one variable*, Colloquium
Mathematicum **175** (2024), 137--151,
[DOI 10.4064/cm9195-1-2024](https://doi.org/10.4064/cm9195-1-2024), the
[journal record](https://www.impan.pl/en/publishing-house/journals-and-series/colloquium-mathematicum/all/175/1/115518/global-injectivity-of-planar-non-singular-maps-that-are-polynomial-in-one-variable),
and [arXiv:2302.05394](https://arxiv.org/abs/2302.05394).

These sources establish substantial prior mathematical coverage. This
repository therefore presents Milestone 6 as an explicit, field-uniform
normal-form reproduction and formalization, with no novelty or priority
claim. Most importantly, the certified theorem assumes that the quadratic
leading coefficient is already a nonzero scalar.  The separate
nonconstant-leading descent remains outside this module until its rational
function and denominator arguments are formalized.  Excluding the certified
ansatz does not settle the unrestricted generic-degree-six frontier and does
not solve the two-dimensional Jacobian conjecture.
