# A sparse-interior obstruction in the residual `(72,108)` configurations

Status: **derived and exhaustively checked in this repository**. This is a
finite sparse-support lemma, not a proof of the two-dimensional Jacobian
conjecture and not a claim of peer-reviewed novelty.

## Source boundary

Guccione, Guccione, Horruitiner, and Valqui showed that, among the coordinate
degree pairs with maximum degree below `125`, their reductions leave only
`(72,108)` and its transpose. Proposition 4.3 of their paper transforms the
remaining `(8,28)` configuration into Laurent coordinates in which

\[
  [P,Q] = P_xQ_y-P_yQ_x=x^2
\]

and the Newton polygons are one of the following two pairs:

\[
\begin{aligned}
\text{Case 1: }N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16),(0,8)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24),(0,12)\};\\[3pt]
\text{Case 2: }N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24)\}.
\end{aligned}
\]

Primary source: [J. A. Guccione, J. J. Guccione, C. Horruitiner, and C.
Valqui, *Increasing the degree of a possible counterexample to the Jacobian
Conjecture from 100 to 108* (2022), Proposition
4.3](https://arxiv.org/abs/2204.14178).

These exponents belong to the transformed Newton-polygon problem. The
ordinary coordinate degree pair `(72,108)` and the generic sheet degree
`[k(x,y):k(P,Q)]` are different invariants.

## Result

Let `supp(P)` and `supp(Q)` denote the exponent sets of nonzero coefficients.
Over a field of characteristic zero, Proposition 4.3 Case 1 with
`[P,Q]=x^2` must satisfy

\[
\begin{aligned}
&\#\bigl(\operatorname{supp}(P)\cap\operatorname{int}N(P)\bigr)\\
&\qquad+\#\bigl(\operatorname{supp}(Q)\cap\operatorname{int}N(Q)\bigr)
\ge 3.
\end{aligned}
\]

In Case 2 the stronger lower bound is

\[
\begin{aligned}
&\#\bigl(\operatorname{supp}(P)\cap\operatorname{int}N(P)\bigr)\\
&\qquad+\#\bigl(\operatorname{supp}(Q)\cap\operatorname{int}N(Q)\bigr)
\ge 4.
\end{aligned}
\]

Thus adding just one or two interior coefficients is insufficient in either
case, and adding three is still insufficient in Case 2. All boundary lattice
coefficients are allowed in the audit; coefficients at the listed polygon
vertices are nonzero because those points are vertices of the exact Newton
polygons.

The lattice-point census is:

| configuration | polynomial | all lattice points | boundary | strict interior |
|---|---:|---:|---:|---:|
| Case 1 | `P` | 61 | 26 | 35 |
| Case 1 | `Q` | 125 | 38 | 87 |
| Case 2 | `P` | 25 | 18 | 7 |
| Case 2 | `Q` | 47 | 26 | 21 |

Consequently the exact support exhaustion checks

\[
  1+122+\binom{122}{2}=7504
\]

patterns in Case 1. Case 2 checks

\[
  1+28+\binom{28}{2}+\binom{28}{3}=3683
\]

patterns, including all `3276` supports with exactly three interior terms.
Zero-product propagation proves `3271` of those triples; the five exceptions
have the exact algebraic certificates below.

## Certificate mechanism

Write `p_ij` and `q_kl` for the coefficients of `x^i y^j` in `P` and
`x^k y^l` in `Q`. Their monomial pair contributes

\[
  (i\ell-jk)p_{ij}q_{k\ell}
\]

to the coefficient of `x^(i+k-1)y^(j+ell-1)` in `[P,Q]`. Every coefficient
other than `(2,0)` must therefore vanish.

The checker repeatedly uses only the following implication:

> If every term but `cuv` in one zero coefficient equation has already
> vanished, `c` is a nonzero integer, and `u` is known nonzero, then `v=0`.

Eventually the sole surviving term has two known-nonzero factors, which is a
contradiction. Each implication sequence is replayed from rebuilt coefficient
equations and a fresh zero state before it is included in the audit digest.

The boundary-only cases also have short hand certificates.

In Case 1, the coefficient at exponent `(9,16)` is exactly

\[
  -24p_{8,16}q_{2,1},
\]

which cannot vanish because both exponents are polygon vertices.

In Case 2, three coefficient equations suffice. The coefficient at `(9,15)`
is

\[
  -22p_{8,15}q_{2,1}=0,
\]

so `p_(8,15)=0`. The coefficient at `(19,35)` is

\[
  8p_{8,14}q_{12,22}-12p_{8,15}q_{12,21}=0,
\]

so `q_(12,22)=0`. Finally the coefficient at `(19,37)` is

\[
  24p_{8,14}q_{12,24}
  +4p_{8,15}q_{12,23}
  -16p_{8,16}q_{12,22}=0,
\]

leaving the nonzero vertex product `24p_(8,14)q_(12,24)`, a contradiction.

### The five exceptional three-interior supports in Case 2

For every Case 2 support under discussion, the `(2,0)` coefficient is exactly
`p_(1,0)q_(2,1)=1`. Rescaling `P` and `Q` inversely therefore lets us normalize

\[
p_{1,0}=q_{2,1}=1
\]

without changing the bracket, support, or any nonvanishing condition. Write
`E_(r,s)` for the coefficient of `x^r y^s` in `[P,Q]`; all the equations below
have right-hand side zero.

The zero-product checker stops on exactly these five triples:

1. `(P[1,1],Q[2,2],Q[2,3])`;
2. `(P[1,1],Q[2,2],Q[3,5])`;
3. `(P[1,1],Q[2,2],Q[5,9])`;
4. `(P[2,3],Q[3,4],Q[4,7])`;
5. `(P[4,7],Q[5,8],Q[8,15])`.

The first three have short identities. For the first support,

\[
\begin{aligned}
E_{2,1}&=2q_{2,2}-p_{1,1},\\
E_{2,2}&=3q_{2,3}-3p_{1,2},\\
E_{2,3}&=4q_{2,4}+p_{1,1}q_{2,3}-2p_{1,2}q_{2,2},\\
E_{2,4}&=2p_{1,1}q_{2,4}-p_{1,2}q_{2,3}.
\end{aligned}
\]

The first three equations successively give `p_(1,1)=2q_(2,2)`,
`p_(1,2)=q_(2,3)`, and `q_(2,4)=0`; the last becomes
`-q_(2,3)^2=0`, contrary to the selected coefficient being nonzero.

For the second support, after the freshly replayed one-term zeros remove
the nuisance terms,

\[
E_{2,1}=2q_{2,2}-p_{1,1},\qquad
E_{3,4}=5q_{3,5}-6p_{2,4},\qquad
E_{3,5}=2p_{1,1}q_{3,5}-4p_{2,4}q_{2,2}.
\]

Substitution leaves a nonzero rational multiple of
`q_(2,2)q_(3,5)`. For the third support the analogous equations are

\[
E_{2,1}=2q_{2,2}-p_{1,1},\qquad
E_{5,8}=9q_{5,9}-12p_{4,8},\qquad
E_{5,9}=4p_{1,1}q_{5,9}-8p_{4,8}q_{2,2},
\]

which leave a nonzero rational multiple of `q_(2,2)q_(5,9)`.

For the fourth support, `E_(3,3)` gives
`p_(2,3)=q_(3,4)=a != 0`. Exact successive reductions then give

\[
\begin{array}{rcl}
-5p_{3,4}+5q_{4,5}&=&0,\\
-2a q_{4,5}&=&0,\\
-8p_{4,6}+7q_{5,7}&=&0,\\
-a(q_{5,7}+2p_{4,6})&=&0,\\
-4a p_{5,8}&=&0,\\
9q_{6,9}&=&0,\\
-14p_{6,10}+11q_{7,11}&=&0,\\
a(q_{7,11}-6p_{6,10})&=&0,\\
-17p_{7,12}+13q_{8,13}&=&0,\\
-5p_{7,12}q_{8,13}&=&0.
\end{array}
\]

Hence all displayed boundary variables vanish. Finally

\[
-20p_{8,14}+15q_{9,15}=0,
\qquad
-6p_{8,14}q_{9,15}=0
\]

forces `p_(8,14)=0`, contradicting a Case 2 polygon vertex. Notice that the
selected `q_(4,7)` is not needed, so it cannot rescue the support.

For the fifth support put

\[
\alpha=p_{4,7}\ne0,\quad \beta=q_{5,8},\quad
a=p_{2,2},\quad b=q_{3,3}.
\]

Three coefficients are

\[
E_{5,7}=8\beta-10\alpha,\qquad
E_{3,2}=3b-2a,\qquad
E_{6,9}=6a\beta-9\alpha b.
\]

The exact syzygy

\[
4E_{6,9}-3aE_{5,7}+12\alpha E_{3,2}=6\alpha a
\]

forces `a=b=0`. Meanwhile, for `k=1,...,6`, the upper-boundary equations

\[
E_{k+1,2k}=-3\sum_{i=1}^{k}i\,p_{i,2i}
q_{k-i+2,2(k-i)+1}=0
\]

are triangular because `q_(2,1)=1`; they force

\[
p_{1,2}=p_{2,4}=p_{3,6}=p_{4,8}=p_{5,10}=p_{6,12}=0.
\]

The full `(9,16)` equation then reduces to `-24p_(8,16)=0`, contradicting
the other upper vertex of `N(P)`. The selected `q_(8,15)` is not needed.

The checker encodes the relevant output exponents and reconstructs every
equation from the full supports. It first replays all one-term zeros, adjoins
formal inverses only for coefficients known nonzero from the exact support,
and verifies over `QQ` that each of these five small ideals has Gröbner basis
`[1]`. Thus the exceptional supports are proved, not merely placed on an
allowlist.

## Reproduction and adversarial checks

Run the audit directly:

```bash
uv run python scripts/newton_72_108.py
```

The current exact summaries are:

```text
Proposition 4.3, case 1
  verified supports with <= 2 interior terms: 7504/7504
  zero-product/algebraic certificates: 7504/0
  maximum forced-zero steps: 21
  certificate sha256: 7c7a3824ac64f0f94989125079be5b4ae7a54817c7eaa4c3b3f2bf79b5bb1519
Proposition 4.3, case 2
  verified supports with <= 3 interior terms: 3683/3683
  zero-product/algebraic certificates: 3678/5
  maximum forced-zero steps: 22
  certificate sha256: ee64dd2ba8035ee29c4efd86ef7f80fc049f8ecd5d32aec638f719cfd31f83fb
```

Run the regression suite and strict type checker:

```bash
uv run pytest -q tests/test_newton_72_108.py
uv run mypy scripts/newton_72_108.py tests/test_newton_72_108.py
```

The hostile tests deliberately relax the ansatz to every lattice point in the
two full polygons while requiring only the vertices to be nonzero. The
zero-product procedure then stops without a contradiction in both cases. A
specific Case 2 support with four interior terms,
`(P[2,3],P[3,5],Q[3,4],Q[4,7])`, also stops the procedure. Those checks
protect the claim boundary: the audit proves the Case 1 `>= 3` and Case 2
`>= 4` sparsity obstructions and nothing about the unresolved full coefficient
systems.

## What remains open

The result eliminates the sparsest possible realizations of the residual
Newton polygons. It does not eliminate either full Proposition 4.3 case. The
next exact target is the finite coefficient system with at least three
interior terms in Case 1 or at least four in Case 2, or a structural argument
that forces a larger interior support before Gröbner elimination.

For comparison, the 2025 paper [*The Groebner basis and solution set of a
polynomial system related to the Jacobian
conjecture*](https://arxiv.org/abs/2506.05697) studies a different Laurent
series specialization (`n=3`, `3` not dividing `m`, and `nu_i=0` for `i>0`).
It does not settle these two Proposition 4.3 Newton-polygon systems.
