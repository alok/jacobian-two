# A variable-leading quadratic coordinate

Status: **DERIVED / FORMALIZATION TARGET.** The argument below has been
algebraically and adversarially checked, but it is **not Lean-certified**.
It is a proposed next formalization milestone, not a solution of the full
plane Jacobian conjecture and not a claim of mathematical novelty. In
particular, the conclusion belongs to the known low-degree positive territory;
the purpose of this note is to make the reduction precise enough to formalize
and to expose every place where a rational denominator could invalidate a
naive proof.

## Statement and conventions

Let `K` be a field of characteristic zero, let `R=K[x]`, and let
`L=K(x)` be its fraction field. Use the Jacobian convention

\[
  J(P,Q)=P_xQ_y-P_yQ_x.
\]

Consider

\[
  P\in K[x,y],\qquad
  Q=a(x)y^2+g(x)y+f(x),qquad
  0\ne a\in K[x].
\]

The derived theorem is:

> If `J(P,Q)=k` for some `k in K^times`, then `a in K^times`.
> Consequently the constant-leading quadratic theorem applies, so `(P,Q)` is
> a polynomial automorphism.

Together with the already certified affine-coordinate theorem, this would
show that a plane Keller map is invertible whenever one component has degree
at most two in one of the variables, while the other component has arbitrary
degree. Only the affine and constant-leading portions of that combined claim
are currently Lean-certified; the reduction in this note is not.

## Denominator ledger

Several operations below take place in `L[U]`, not in `K[x,U]`. These
distinctions are proof obligations, not cosmetic details.

- Characteristic zero makes `2`, every positive integer, and every binomial
  coefficient used below nonzero in `K`; because `K` is a field, the relevant
  integers are invertible.
- The constants of the derivation `d/dx` on `K(x)` are exactly `K`. This fails
  in positive characteristic and is used repeatedly.
- Division by `a`, `h`, and powers of them is initially allowed only in
  `K(x)`. It does not produce a polynomial coordinate change.
- The variable `U` introduced below is a polynomial in `(x,y)` only after one
  proves `h | g`; even then, recovering `y` from `U` still divides by `h`.
  Thus `h | g` alone does **not** make `(x,U)` a polynomial coordinate system.
- The final coefficient equation, not the rational change of variables,
  forces `h` to be a unit. Skipping that step leaves a real gap.
- Every valuation is the discrete exponent valuation associated to an
  irreducible polynomial of `K[x]`; it is trivial on `K^times`. Strictly
  different valuations, rather than an informal leading-term comparison,
  prevent the cancellations used in the contradiction.

## 1. Descend to an odd `y`-degree

Write

\[
  P=\sum_{i=0}^{n}p_i(x)y^i,
  \qquad p_n\ne0.
\]

Since

\[
  Q_y=2ay+g,
  \qquad
  Q_x=a'y^2+g'y+f',
\]

the coefficient of `y^(n+1)` in `J(P,Q)` is

\[
  2a p_n'-n a'p_n.
\]

The right side `k` has `y`-degree zero, so

\[
  2a p_n'-n a'p_n=0. \tag{1}
\]

The case `n=0` is impossible: (1) gives `p_0'=0`, hence `P` is a scalar and
its Jacobian with `Q` is zero. Therefore `n>0`.

If `n=2r` is even, (1) says

\[
  \left(\frac{p_n}{a^r}\right)'=0
\]

in `K(x)`. Hence `p_n=c a^r` for some `c in K^times`. The leading
`y`-coefficient of `Q^r` is `a^r`, and

\[
  J(Q^r,Q)=0.
\]

Replacing `P` by `P-cQ^r` therefore preserves `J(P,Q)=k` and strictly lowers
the `y`-degree. Repeating this operation must terminate. It cannot terminate
at degree zero by the preceding paragraph, so after a target shear
`P -> P-H_0(Q)` with `H_0 in K[T]`, the surviving degree is odd:

\[
  n=2m+1 \quad (m\ge0). \tag{2}
\]

This replacement is harmless for the desired conclusion about `a`; it is
also reversible by a polynomial target shear once invertibility is known.
For the rest of the proof, `P` denotes the odd-degree remainder.

## 2. The leading coefficient makes `a` a square up to a unit

For odd `n`, equation (1) gives

\[
  \left(\frac{p_n^2}{a^n}\right)'=0.
\]

Thus

\[
  p_n^2=c a^n \tag{3}
\]

for some `c in K^times`. Factor in the UFD `K[x]`. For every irreducible
`pi`, (3) gives

\[
  2v_\pi(p_n)=n v_\pi(a).
\]

Because `n` is odd, every `v_pi(a)` is even. Consequently there are

\[
  \varepsilon,\lambda\in K^\times,
  \qquad 0\ne h\in K[x]
\]

such that

\[
  a=\varepsilon h^2,
  \qquad
  p_n=\lambda h^n. \tag{4}
\]

No algebraic closure or root extraction in `K` is being assumed. The scalar
`epsilon` is simply the unit left after halving all irreducible exponents of
`a`; the second identity follows by comparing the exponents of `p_n` and
`h^n`, leaving another scalar unit.

Equation (4) does not yet imply that `h` is constant. Indeed, the top-degree
equation alone admits nonconstant examples; one is recorded in the hostile
fixtures below.

## 3. Complete the square over `K(x)` and split parity

Work temporarily over `L=K(x)`. Define

\[
  u_0=\frac{g}{2\varepsilon h},
  \qquad
  U=hy+u_0,
  \qquad
  F=f-\varepsilon u_0^2.
\]

Then

\[
  Q=\varepsilon U^2+F. \tag{5}
\]

At this stage `u_0` and `F` may be rational functions. The change
`y=(U-u_0)/h` is valid in `L[U]`, but not necessarily in `K[x,U]`.

The determinant of `(x,y) -> (x,U)` is `h`. A direct chain-rule expansion,
with `x`-derivatives in the new coordinates taken at fixed `U`, gives

\[
  J_{x,y}(P,Q)=h J_{x,U}(P,Q).
\]

Therefore

\[
  J_{x,U}(P,Q)=\frac{k}{h}. \tag{6}
\]

Because `Q=epsilon U^2+F` and `epsilon` is a unit, `L[U]` decomposes uniquely
as

\[
  L[U]=L[Q]\oplus U L[Q].
\]

Write the unique parity decomposition

\[
  P=H(x,Q)+U B(x,Q),
  \qquad H,B\in L[Q]. \tag{7}
\]

Here subscripts `x` on `H` and `B` mean coefficientwise differentiation at
fixed `Q`. From (4), the leading `U`-coefficient of `P` is `lambda`. Since
`n=2m+1`, if

\[
  B=\sum_{j=0}^{m}b_j(x)Q^j,
\]

then

\[
  b_m=\lambda\varepsilon^{-m}\in K^\times. \tag{8}
\]

Differentiate (7), use `Q_U=2epsilon U` and `Q_x=F'` at fixed `U`, and cancel
the terms involving `H_Q` and `B_Q`. With the stated Jacobian convention,
the result is

\[
\begin{aligned}
  J_{x,U}(P,Q)
    &=2\varepsilon U H_x+2\varepsilon U^2 B_x-F'B \\
    &=2\varepsilon U H_x+2(Q-F)B_x-F'B. \tag{9}
\end{aligned}
\]

The first summand in (9) is odd in `U`; all remaining terms, including the
right side `k/h`, are even. Uniqueness of the parity decomposition yields

\[
  H_x=0,
  \qquad
  2(Q-F)B_x-F'B=\frac{k}{h}. \tag{10}
\]

The derivative constants of `K(x)` are `K`, so the first identity in (10)
means

\[
  H\in K[Q]. \tag{11}
\]

This fact is essential later: it makes `H(f)` polynomial when we specialize
to `y=0`.

## 4. Solve the coefficient recurrence exactly

Insert

\[
  B=\sum_{j=0}^{m}b_jQ^j
\]

into the even equation in (10). Comparing powers of `Q` gives

\[
\begin{aligned}
  b_m'&=0, \\[-2pt]
  b_{j-1}'&=F b_j'+\frac{F'}2 b_j
    &&(1\le j\le m), \\[-2pt]
  -2F b_0'-F'b_0&=\frac{k}{h}. \tag{12}
\end{aligned}
\]

Starting with the nonzero constant `b_m` from (8), downward induction in
(12) proves

\[
  b_j=C_j(F),
  \qquad C_j\in K[T],
  \qquad \deg C_j=m-j. \tag{13}
\]

Here is the integration step with its constants exposed. If
`b_j=C_j(F)`, then

\[
  b_{j-1}'
    =\left(F C_j'(F)+\frac12 C_j(F)\right)F'.
\]

Choose `C_{j-1} in K[T]` with

\[
  C_{j-1}'(T)=T C_j'(T)+\frac12 C_j(T). \tag{14}
\]

Such an antiderivative exists because `K` has characteristic zero. The
difference `b_{j-1}-C_{j-1}(F)` has zero derivative in `K(x)`, so it lies in
`K` and can be absorbed into the constant term of `C_(j-1)`.

If the leading coefficient of `C_j` is `c` and its degree is `d`, equation
(14) makes the leading coefficient of `C_(j-1)` equal to

\[
  c\,\frac{2d+1}{2d+2},
\]

which is nonzero. Thus the degrees in (13) are exact, not merely upper
bounds. In particular,

\[
  \operatorname{lc}(C_0)
    =b_m\prod_{r=0}^{m-1}\frac{2r+1}{2r+2}
    =b_m\frac{\binom{2m}{m}}{4^m}
    \ne0. \tag{15}
\]

For `m=0`, the empty product is one: `B=b_0 in K^times`, and the last line of
(12) remains the decisive equation.

## 5. A valuation argument forces `h | g`

Suppose, for contradiction, that `h` does not divide `g`. The case `g=0`
cannot occur under this supposition because every polynomial divides zero.
Choose an irreducible `pi in K[x]` for which

\[
  d:=v_\pi(h)-v_\pi(g)>0.
\]

Then

\[
  v_\pi(u_0)=-d.
\]

Since `f` is a polynomial, `v_pi(f)>=0`, whereas
`v_pi(epsilon u_0^2)=-2d`. The valuations are different, so there can be no
cancellation in their difference:

\[
  v_\pi(F)=v_\pi(f-\varepsilon u_0^2)=-2d. \tag{16}
\]

Now specialize (7) at `y=0`. Then `U=u_0` and `Q=f`, so (11) gives

\[
  u_0B(x,f)=P(x,0)-H(f)\in K[x]. \tag{17}
\]

On the other hand, (13) writes the left side as

\[
  u_0B(x,f)=\sum_{j=0}^{m}u_0 C_j(F)f^j. \tag{18}
\]

Because `v_pi(F)<0` and `C_j` has exact degree `m-j` with nonzero scalar
leading coefficient,

\[
  v_\pi(C_j(F))=-2d(m-j).
\]

The `j=0` summand of (18) therefore has valuation

\[
  -(2m+1)d. \tag{19}
\]

For every `j>=1` whose summand is nonzero, its valuation exceeds (19) by

\[
  2dj+jv_\pi(f)>0.
\]

If `f=0`, those higher summands simply vanish. Thus the `j=0` term is the
unique term of least valuation; it cannot cancel. Equation (18) has negative
valuation `-(2m+1)d`, contradicting the polynomial membership in (17).
Hence

\[
  h\mid g. \tag{20}
\]

This specialization is the denominator-clearing step. Merely writing the
completed square over `K(x)` does not provide it.

## 6. The constant coefficient forces `h` to be a unit

From (20),

\[
  u_0=\frac{g}{2\varepsilon h}\in K[x],
  \qquad
  F=f-\varepsilon u_0^2\in K[x].
\]

Equations (13) now show that every `b_j=C_j(F)` lies in `K[x]`. In particular,
the left side of the final coefficient equation in (12) is a polynomial:

\[
  -2F b_0'-F'b_0\in K[x].
\]

But that equation says the polynomial equals `k/h` in `K(x)`. Therefore
`k/h in K[x]`, so `h` divides the nonzero scalar `k` in `K[x]`. The only such
divisors are units:

\[
  h\in K^\times.
\]

Finally, (4) gives

\[
  a=\varepsilon h^2\in K^\times,
\]

as claimed. Notice the order of the last two steps: after `h | g`, `U` is a
polynomial expression, but its inverse still contains `1/h`. The coefficient
recurrence first makes the left side of (12) polynomial; only then does the
Keller constant prove that `1/h` is a scalar and close the gap.

## Edge cases and scope checks

- **`a=0` is deliberately excluded.** Then `Q` is affine in `y` and belongs
  to the already Lean-certified affine-coordinate theorem.
- **`k=0` is deliberately excluded.** The descent still yields some top
  relations, but the final divisibility argument has no force; for instance,
  `P=Q` has zero Jacobian for arbitrary `a`, `g`, and `f`.
- **A `y`-independent mate is impossible.** This is the `n=0` branch excluded
  in Step 1, not an omitted base case.
- **The smallest odd branch is included.** When `n=1`, one has `m=0`; the
  empty central-binomial product is one, the valuation argument still works,
  and the last line of (12) still forces `h` to be a unit.
- **`g=0` is included.** In that case `h | g` is automatic, `u_0=0`, and
  Step 6 applies directly. One must not try to choose a valuation of the zero
  polynomial in Step 5.
- **`f=0` is included.** In Step 5 the terms with `j>=1` vanish after
  specialization, while the `j=0` term retains its strictly negative
  valuation.
- **No algebraic closure is assumed.** UFD factor exponents, rational-function
  differentiation, and irreducible-polynomial valuations all live over the
  original characteristic-zero field `K`.
- **Positive characteristic is outside the claim.** Derivative constants and
  the nonvanishing integer factors in the descent and recurrence can fail
  there.

## Hostile fixtures

These examples are regression tests for tempting but invalid shortcuts.

### A nonsingular quadratic with no polynomial Keller mate

Let

\[
  Q=x+x^2y^2.
\]

Its gradient never vanishes simultaneously: `Q_y=2x^2y`, while
`Q_x=1+2xy^2`; if `x=0` or `y=0`, the latter is one. Thus nonsingularity of a
single component is much weaker than having a polynomial Keller mate. The
derived theorem rules out such a mate because its leading coefficient
`a=x^2` is not constant.

There is also a direct infinite-tail obstruction. For

\[
  P=\sum_{j\ge0}p_j(x)y^j
\]

with finite support, the constant coefficient of `J(P,Q)=k` gives

\[
  p_1=-k.
\]

For `ell>=1`, the coefficient of `y^ell` is

\[
  2x^2p_{\ell-1}'-2x(\ell-1)p_{\ell-1}
    -(\ell+1)p_{\ell+1}=0. \tag{21}
\]

Taking `ell=2r` and inducting from `p_1=-k` gives

\[
  p_{2r+1}=-\frac{2r}{2r+1}x p_{2r-1}
  \qquad(r\ge1). \tag{22}
\]

Every factor in (22) is nonzero in characteristic zero, so all odd
coefficients are nonzero. This contradicts finite support and independently
confirms that no polynomial `P` can have constant nonzero Jacobian with this
`Q`.

### The top equation alone does not make `a` constant

For any positive integer `n`, take

\[
  a=x^2,
  \qquad p_n=x^n.
\]

Then

\[
  2a p_n'-n a'p_n
    =2x^2(nx^{n-1})-n(2x)x^n
    =0.
\]

Thus the leading-coefficient identity (1), even for odd `n`, permits a
nonconstant square. The lower coefficient equations and their denominator
constraints are indispensable.

### A rational Jacobian-one mate is not a polynomial mate

Let

\[
  Q=x+x^3y^2,
  \qquad
  R=-\frac{xy}{Q}.
\]

In the rational function field `K(x,y)`, quotient-rule cancellation gives

\[
\begin{aligned}
  J(R,Q)
    &=\frac{xQ_x-yQ_y}{Q} \\
    &=\frac{x(1+3x^2y^2)-y(2x^3y)}{x+x^3y^2}
      =1.
\end{aligned}
\]

The denominator is real: `R` is not in `K[x,y]`. This fixture shows why a
rational normal form or a rational inverse cannot be promoted to a polynomial
statement without the valuation and divisibility arguments in Steps 5 and 6.

### The forced coefficient tail in completed-square form

For `Q=x+x^2y^2`, one has `epsilon=1`, `h=x`, `u_0=0`, and `F=x`. The
recurrence (12) builds polynomial coefficients `b_j=C_j(x)`, but its constant
line would require a polynomial to equal `k/x`. Equivalently, the direct
recurrence (22) produces an infinite nonzero odd tail. These are two views of
the same obstruction: the rational calculation is formally consistent until
the nonunit denominator `h=x` meets the nonzero Keller constant.

## Formalization boundary

A Lean proof should not encode only the final recurrence as a premise. The
intended certificate starts from the formal bivariate Jacobian and verifies:

1. even-degree target-shear descent to an odd `y`-degree;
2. the UFD valuation consequence `a=epsilon*h^2`;
3. the parity decomposition in `K(x)[U]` and the chain-rule factor `h`;
4. the exact coefficient recurrence and central-binomial leading coefficient;
5. the specialization-at-`y=0` valuation contradiction forcing `h | g`; and
6. the final divisibility argument forcing `h` to be a unit.

Until those steps are kernel-checked, this document remains a derived proof
target. It may guide experiments and formalization, but it must not be cited
as a Lean theorem or as a new mathematical result.

## Literature and tracking boundary

Formalization of the variable-leading reduction is tracked in
[ALOK-799](https://linear.app/aloksingh/issue/ALOK-799). That issue identifier
is project tracking, not evidence for the theorem.

Two primary sources already ground the adjacent positive territory in this
repository:

- Denis Simon and Martin Weimann, *Plane Curves With Minimal Discriminant*,
  Journal of Commutative Algebra **10** (2018), 559--598,
  [DOI 10.1216/JCA-2018-10-4-559](https://doi.org/10.1216/JCA-2018-10-4-559)
  and [arXiv:1507.01091](https://arxiv.org/abs/1507.01091). Their monic
  minimal-discriminant theorem covers the constant-leading coordinate after
  the reduction, not the rational-denominator argument of Steps 1--6.
- Marco Sabatini, *Global injectivity of planar non-singular maps that are
  polynomial in one variable*, Colloquium Mathematicum **175** (2024),
  137--151,
  [DOI 10.4064/cm9195-1-2024](https://doi.org/10.4064/cm9195-1-2024) and
  [arXiv:2302.05394](https://arxiv.org/abs/2302.05394). His real bounded-degree
  results overlap important subcases, but they are not being quoted here as
  this field-uniform arbitrary-mate statement.

Those overlaps are why this repository makes no novelty or priority claim.
The proof in this note is independently derived, remains unformalized, and is
presented only as the next certificate target.
