# Hostile models for the smooth normalization packages

## Claim boundary

The smooth-normalization reductions in this repository produce affine surface
pairs with one of the following forms:

\[
\begin{aligned}
&W\setminus D\cong \mathbb A^2,
&&D\cong\mathbb A^1,
&&\operatorname{Pic}(W)=\mathbb Z[D],
&&\omega_W\cong\mathcal O_W(D)
  \text{ or }\mathcal O_W(2D); \\[2mm]
&W\setminus(D_2\cup D_3)\cong\mathbb A^2,
&&D_2,D_3\cong\mathbb A^1\text{ disjoint},
&&\operatorname{Pic}(W)=\mathbb Z[D_2]\oplus\mathbb Z[D_3],
&&\omega_W\cong\mathcal O_W(D_2+2D_3).
\end{aligned}
\tag{0.1}
\]

This note constructs explicit smooth affine pairs having exactly these
properties.  Consequently, the pair-theoretic data in (0.1) are consistent
and cannot by themselves eliminate a degree-six Keller candidate.

These are deliberately hostile consistency models, **not** Jacobian
counterexamples.  In particular, the constructions below do not provide a
finite flat rank-six morphism `W -> A2`, a rank-six algebra over `C[P,Q]`,
coordinate functions satisfying \(dP\wedge dQ=c\,dx\wedge dy\), the prescribed
branch curves or inertia, or transitive `A6`/`S6` monodromy.  Those missing
finite-algebra and map data remain the essential restrictions.

## 1. A one-boundary Hirzebruch family

Let `X=F_d` be the Hirzebruch surface.  Write `C_0` for its negative section
and `F` for the class of a fiber, so

\[
C_0^2=-d,\qquad C_0\cdot F=1,\qquad F^2=0,
\qquad
K_X=-2C_0-(d+2)F.
\tag{1.1}
\]

Fix an integer `a>d`.  The linear system `|C_0+aF|` is base-point-free, and a
general member

\[
S\sim C_0+aF
\tag{1.2}
\]

is a smooth irreducible section.  It is ample because a divisor
`u C_0+vF` on `F_d` is ample exactly when `u>0` and `v>ud`.  Set

\[
W:=X\setminus S.
\tag{1.3}
\]

The complement of an effective ample divisor on a projective variety is
affine, so `W` is smooth affine.  Choose a fiber `F_0` and put

\[
D:=F_0\cap W=F_0\setminus(F_0\cap S).
\tag{1.4}
\]

Because a section meets every fiber once, `D` is isomorphic to `A1`.
Moreover,

\[
W\setminus D=X\setminus(S\cup F_0)\cong\mathbb A^2.
\tag{1.5}
\]

One way to see (1.5) is to restrict the ruling `X -> P1` to
`P1\setminus\{[F_0]\}=A1`.  Removing the section leaves an `A1`-bundle over
`A1`.  It is a torsor under a line bundle; every line bundle on `A1` is
trivial and its first cohomology vanishes, so this torsor is `A1 x A1`.

Divisor localization for the complement of the irreducible divisor `S`
gives

\[
\operatorname{Pic}(W)
\cong
\frac{\mathbb Z[C_0]\oplus\mathbb Z[F]}
     {\mathbb Z[C_0+aF]}
\cong \mathbb Z[D].
\tag{1.6}
\]

The last group is torsion-free because `(1,a)` is primitive.  In this
quotient `[C_0]=-a[F]`; hence (1.1) yields

\[
\begin{aligned}
[K_W]
  &= -2[C_0]-(d+2)[F] \\
  &= (2a-d-2)[D] \\
  &= (S^2-2)[D],
\end{aligned}
\qquad S^2=-d+2a.
\tag{1.7}
\]

Thus

\[
\boxed{
\omega_W\cong\mathcal O_W\bigl((2a-d-2)D\bigr),\quad
\operatorname{Pic}(W)=\mathbb Z[D],\quad
W\setminus D\cong\mathbb A^2.
}
\tag{1.8}
\]

Two small members of the family realize the exact one-boundary packages:

| `d` | `a` | `S^2=-d+2a` | `omega_W` |
|---:|---:|---:|:---|
| 1 | 2 | 3 | `O_W(D)` |
| 0 | 2 | 4 | `O_W(2D)` |

The first is the canonical-class pattern forced by the one-dicritical `S6`
case, and the second is the pattern forced by the one-dicritical `A6` case.

## 2. A two-boundary triple-blowup model

The two-boundary package also occurs on an explicit affine surface.  Begin
with

\[
X_0=\mathbb P^1\times\mathbb P^1,
\qquad
\operatorname{Pic}(X_0)=\mathbb Z[A]\oplus\mathbb Z[B],
\tag{2.1}
\]

where

\[
A^2=B^2=0,\qquad A\cdot B=1.
\tag{2.2}
\]

Let `S` be the diagonal, of class `A+B`, and let `F` be a fiber of class `A`.
Choose distinct points

\[
p_1,p_2\in F\setminus S.
\tag{2.3}
\]

Blow up `p_1` and `p_2`, with exceptional curves `E_1` and `E_2`.  Then blow
up a point

\[
p_3\in E_2\setminus(E_2\cap F'),
\tag{2.4}
\]

where `F'` is the strict transform of `F`; call the last exceptional curve
`E_3`.  On the resulting smooth projective surface `X`, use `E_2` for the
total-transform class of the second exceptional curve.  The relevant strict
transforms have classes

\[
S\sim A+B,\qquad
F'\sim A-E_1-E_2,
\qquad
E_2'\sim E_2-E_3.
\tag{2.5}
\]

Set

\[
H:=S\cup F'\cup E_2',
\qquad
W:=X\setminus H,
\tag{2.6}
\]

and retain the two exceptional curves not removed with `H`:

\[
D_2:=E_1\cap W,
\qquad
D_3:=E_3\cap W.
\tag{2.7}
\]

The curve `E_1` meets `H` only in its one point on `F'`, while `E_3` meets
`H` only in its one point on `E_2'`.  Thus

\[
D_2\cong D_3\cong\mathbb A^1,
\qquad D_2\cap D_3=\varnothing.
\tag{2.8}
\]

After removing `D_2` and `D_3`, all components in the reduced inverse image
of `F` have been deleted.  The blowups are therefore invisible on the
remaining open set:

\[
W\setminus(D_2\cup D_3)
\cong X_0\setminus(S\cup F)
\cong\mathbb A^2.
\tag{2.9}
\]

For the final isomorphism, move `F` to the fiber at infinity of the first
projection.  Over its affine complement, removing the diagonal leaves an
`A1`-bundle over `A1`, hence a trivial one.  Explicitly, if `x` is the affine
coordinate on the first factor and `[Y_0:Y_1]` are coordinates on the second,
then the diagonal has equation `Y_1-xY_0=0`; on its complement

\[
z=\frac{Y_0}{Y_1-xY_0}
\tag{2.10}
\]

gives the second affine coordinate, with inverse
`[Y_0:Y_1]=[z:1+xz]`.

### 2.1 Picard group and canonical class

In the total-transform basis,

\[
\operatorname{Pic}(X)
=\mathbb Z[A]\oplus\mathbb Z[B]
 \oplus\mathbb Z[E_1]\oplus\mathbb Z[E_2]\oplus\mathbb Z[E_3].
\tag{2.11}
\]

Localizing along the three irreducible components of `H` imposes the
relations

\[
A+B=0,\qquad A-E_1-E_2=0,\qquad E_2-E_3=0.
\tag{2.12}
\]

Consequently

\[
B=-A,qquad E_2=E_3,qquad A=E_1+E_3,
\tag{2.13}
\]

and no torsion is introduced.  Since `E_1\cap W=D_2` and
`E_3\cap W=D_3`, this proves

\[
\boxed{
\operatorname{Pic}(W)
\cong\mathbb Z[D_2]\oplus\mathbb Z[D_3].
}
\tag{2.14}
\]

The canonical class after the three blowups is

\[
K_X=-2A-2B+E_1+E_2+E_3.
\tag{2.15}
\]

Reducing (2.15) with (2.12) gives

\[
[K_W]=[E_1]+2[E_3]=[D_2]+2[D_3],
\tag{2.16}
\]

so

\[
\boxed{
\omega_W\cong\mathcal O_W(D_2+2D_3).
}
\tag{2.17}
\]

### 2.2 Affineness

It remains to prove that the open surface in (2.6) is affine.  The support
`H` carries the effective divisor

\[
\begin{aligned}
L&:=6S+3F'+E_2' \\
 &=9A+6B-3E_1-2E_2-E_3.
\end{aligned}
\tag{2.18}
\]

In the total-transform basis the only nonzero pairings are
`A.B=1` and `E_i^2=-1`.  Therefore

\[
L^2=2(9)(6)-3^2-2^2-1^2=94>0.
\tag{2.19}
\]

Its intersections with the boundary and exceptional curves are

\[
L\cdot S=15,qquad
L\cdot F'=1,qquad
L\cdot E_1=3,qquad
L\cdot E_2'=1,qquad
L\cdot E_3=1.
\tag{2.20}
\]

Now let `C` be any other irreducible curve on `X`, and let its image in `X_0`
have class `aA+bB`.  Since the exceptional curves were already covered, its
class has the form

\[
C\sim aA+bB-m_1E_1-m_2E_2-m_3E_3,
\tag{2.21}
\]

where `a,b>=0` and they are not both zero.  If the image were `F`, then `C`
would be `F'`, already covered by (2.20).  Otherwise, intersection with `F`
at the two distinct first-stage centers gives

\[
m_1+m_2\le (aA+bB)\cdot A=b.
\tag{2.22}
\]

The third center lies on the second exceptional curve, and the first strict
transform has total intersection `m_2` with that curve, so

\[
m_3\le m_2.
\tag{2.23}
\]

Equations (2.18), (2.22), and (2.23) now give

\[
\begin{aligned}
L\cdot C
 &=6a+9b-3m_1-2m_2-m_3 \\
 &\ge 6a+9b-3m_1-3m_2 \\
 &\ge 6a+6b>0.
\end{aligned}
\tag{2.24}
\]

Thus `L` has positive intersection with every irreducible curve and positive
self-intersection.  Nakai--Moishezon makes `L` ample.  Because `L` is an
effective ample divisor with support exactly `H`, its complement

\[
W=X\setminus\operatorname{Supp}(L)
\tag{2.25}
\]

is affine.

Combining (2.8), (2.9), (2.14), (2.17), and (2.25) realizes every geometric
entry of the two-boundary row of (0.1).

## 3. What the examples do and do not prove

The examples prove that smoothness, affineness, an `A2` interior, affine-line
boundary components, the displayed Picard groups, and the displayed
canonical classes have no internal contradiction.  Any elimination argument
must therefore use structure absent from an arbitrary affine pair.

For a hypothetical Keller map that extra structure includes a finite flat
rank-six morphism

\[
\rho:W\longrightarrow\mathbb A^2,
\tag{3.1}
\]

the multiplication and trace pairing on the locally free algebra
\(\rho_*\mathcal O_W\), its discriminant and branch divisor, the restriction
of \(\rho\) to the distinguished \(\mathbb A^2\) interior, and the required local inertia and global
monodromy.  None of (3.1) or these compatible algebraic data is constructed
here.  The hostile models narrow the target: the next obstruction must couple
the surface pair to its finite rank-six Keller algebra, rather than use the
surface package in isolation.
