# A conditional degree-six elimination at infinity

## Claim boundary

Assume all of the following for a hypothetical one-dicritical degree-six
plane Keller counterexample.

1. The irreducible image of the dicritical has **exactly two characteristic
   pairs at infinity**.
2. Orevkov's condition `(*)` holds: there is a component of the target
   compactification boundary whose full preimage contains exactly one
   irreducible component not contracted by the resolved map to a point.
   Additional contracted preimage components are allowed.

Then neither surviving one-dicritical passport is possible:

\[
\boxed{
(A_6,e=3,d=1)\text{ and }(S_6,e=2,d=1)
\text{ are both excluded in this stratum.}
}
\tag{0.1}
\]

This is a genuine conditional elimination, not a solution of the
one-dicritical problem.  Neither hypothesis is currently forced by the
monodromy, smooth finite-normalization, Picard, trace, or minimal-leaf
theorems in this repository.  In particular, “two characteristic pairs”
describes the unique branch at infinity; it does not refer to the finite
`T(2,5)` cusp and does not mean “at most two.”

The proof below is a degree-six specialization of the splice equations in
Section 2.4 of Orevkov's
[“Counter-examples to the Jacobian Conjecture at Infinity”](https://www.math.univ-toulouse.fr/~orevkov/jci-e.pdf).
It is a short hand derivation from the displayed equations, not an appeal to
the source's reported computer search.

## 1. Dictionary and product identity

Orevkov writes `m` for the degree of the resolved map on the dicritical and
`n` for its normal branching order.  In the notation used elsewhere here,

\[
m=d=1,
\qquad
n=e\in\{2,3\},
\qquad
N=6.
\tag{1.1}
\]

The distinction between `n` and the Jacobian vanishing order matters:
Orevkov's Jacobian order is `n-1`.

Write `S_i` for his tilded determinant
`\widetilde S_i^\infty`.  Equations (11) and (14), followed by cancellation
of positive factors, give

\[
n_1=md_2S_1,
\qquad
m_1=mxd_1d_2R_1.
\tag{1.2}
\]

Proposition 1.6 identifies the topological degree as `N=m_1n_1`, so

\[
N=m^2x d_1d_2^2R_1S_1.
\tag{1.3}
\]

Every factor in (1.3) is a positive integer.  The source explicitly records
`R_1>1`, while equation (15) is

\[
x=a+n,
\qquad a\ge1.
\tag{1.4}
\]

Therefore

\[
N\ge 2m^2(n+1).
\tag{1.5}
\]

Condition `(*)` is not decorative here.  It licenses Orevkov's Figures 4--6
and the determinant transport in (14), from which (1.2) and hence the finite
product bound follow.

## 2. The `A6` passport dies immediately

For `A6`, equations (1.1) and (1.5) give

\[
N\ge2(3+1)=8,
\]

contrary to `N=6`.  No further splice data are needed.

## 3. The `S6` near-miss

For `S6`, `n=2`, so (1.5) only says `N>=6`.  Equality in (1.3) forces every
factor to take its minimum:

\[
x=3,quad a=1,quad R_1=2,quad
d_1=d_2=S_1=m=1.
\tag{3.1}
\]

Equation (1.2) then gives

\[
m_1=6,
\qquad
n_1=1.
\tag{3.2}
\]

At the first boundary star, equations (12)--(13) read

\[
k_0R_1=m_1,
\qquad
d_1+k_1D_1=m_1,
\qquad
m'_1=k_1+k_0,
\qquad
m'_1+k'_1=m_1.
\tag{3.3}
\]

The arm counts are nonnegative.  Substituting (3.1)--(3.2) gives

\[
k_0=3,
\qquad
k_1D_1=5,
\qquad
m'_1=k_1+3\le6.
\]

Thus

\[
k_1=1,quad D_1=5,quad m'_1=4,quad k'_1=2.
\tag{3.4}
\]

At the second star, equations (11) and (13) give

\[
S_2=\frac{n_2}{m}=\frac{n}{m'_2},
\qquad
m'_2=k_2+1.
\tag{3.5}
\]

Hence `2/(k_2+1)` is a positive integer.  Even allowing the weaker convention
`k_2>=0`, only `k_2=0` or `k_2=1` remains.

If `k_2=0`, equation (12) gives `m_2=d_2=1`, while (3.5) gives `n_2=2`.
The ratio in equation (11),

\[
\frac{n_1}{m_2}=\frac{n_2}{m'_1},
\]

would say `1=2/4`, a contradiction.

Therefore `k_2=1`, so

\[
m'_2=2,quad S_2=n_2=1.
\tag{3.6}
\]

The same ratio gives `m_2=4`; equation (12) then gives

\[
D_2=3.
\tag{3.7}
\]

Since `d_1=d_2=1`, equation (11) also says

\[
\widetilde D_1=5,
\qquad
\widetilde D_2=3.
\tag{3.8}
\]

The decisive step is the exact edge-determinant equation (10):

\[
-x\widetilde Q_2
=S_1S_2-\widetilde D_1\widetilde D_2.
\tag{3.9}
\]

Substitution yields

\[
-3\widetilde Q_2=1-15=-14,
\qquad
3\widetilde Q_2=14.
\]

This is impossible because `\widetilde Q_2` is an integral edge
determinant.  The `S6` stratum is therefore excluded as well.

## 4. Belyi theory alone would not have killed the near-miss

Before (3.9), the first star has the forced degree-six Belyi passport

\[
(2,2,2),qquad(5,1),qquad(4,1,1).
\tag{4.1}
\]

It is realizable.  For example, on six letters take

\[
\sigma_0=(12)(34)(56),
\qquad
\sigma_1=(23456),
\qquad
\sigma_\infty=(\sigma_0\sigma_1)^{-1}.
\]

Their product is one, their cycle types are exactly (4.1), and they generate
a transitive group of order `120`.  Their Riemann--Hurwitz defects are
`3+4+3=10=2*6-2`.  Thus the final divisibility `3 | 14`, not passport
existence or genus, is the actual obstruction.

The dependency-free arithmetic certificate
[`two_pair_infinity.py`](../scripts/two_pair_infinity.py) checks the `A6`
lower bound, both nonnegative `k_2` cases, the forced fractional determinant
`\widetilde Q_2=14/3`, and the exact permutation near-miss.

## 5. Why this does not yet finish a passport

One dicritical, tangential degree one, smoothness of the affine finite
normalization, and primitive `A6` or `S6` monodromy do not currently imply
condition `(*)`.  A target boundary component can have several
noncontracted source components above it.  Nor has the infinity branch been
proved to have exactly two characteristic pairs: it may have fewer or more.

Without `(*)`, the determinant-product estimate used in Proposition 2.5
disappears.  Local Belyi stars of degree at most six can still be enumerated,
but the number and determinant weights of global Puiseux chains are not
bounded by the present argument.  The next useful global lemma would force a
single noncontracted preimage over some target boundary component, or replace
that condition with a determinant inequality that survives arbitrary
splitting.

## Primary source

- S. Yu. Orevkov,
  [“Counter-examples to the Jacobian Conjecture at Infinity”](https://www.math.univ-toulouse.fr/~orevkov/jci-e.pdf),
  Section 2.4.  Condition `(*)` introduces the splice diagrams; equations
  (10)--(15) and Proposition 2.5 supply the relations specialized above.

The specialization to the two degree-six passports and the explicit
two-case hand proof are derived here and were independently hostile-audited.
