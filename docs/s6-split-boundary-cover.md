# A finite degree-six `S6` cover with split boundary

## Status and purpose

This note constructs an actual finite flat map


\[
F:\mathbb A^2\longrightarrow\mathbb A^2
\]

of degree six whose geometric monodromy is `S6` and whose compactified target
boundary has two noncontracted irreducible preimages.  It proves a precise
limitation on the current compactification attack:

> Orevkov's condition `(*)` does not follow from degree six, finiteness over
> `A2`, `S6` monodromy, or boundary intersection transport—even when all of
> those data occur in one genuine algebraic cover.

The map is **not Keller**: its Jacobian is nonconstant and it is branched in
the affine plane.  It has the wrong branch passport and is not a
counterexample to the Jacobian conjecture.  A Keller-specific input remains
capable of forcing `(*)`.

## 1. The affine finite flat map

Define


\[
P=x^2+z,
\qquad
Q=x^3-x^4z.
\tag{1.1}
\]

Since `z=P-x^2`, the coordinate `x` satisfies


\[
x^6-Px^4+x^3-Q=0.
\tag{1.2}
\]

Conversely, (1.2) together with `z=P-x^2` recovers (1.1).  Hence


\[
\mathbb C[x,z]
\simeq
\mathbb C[P,Q][x]/(x^6-Px^4+x^3-Q).
\tag{1.3}
\]

The defining polynomial is monic of degree six.  The target module is
therefore free with basis


```text
1,x,x^2,x^3,x^4,x^5.
```

Thus `F` is finite flat of degree six.

Its Jacobian is


\[
J_F=-2x^5+4x^3z-3x^2
=x^2(4xz-2x^3-3).
\tag{1.4}
\]

This visibly records the central claim boundary: the cover has affine
ramification and is not Keller.

## 2. A finite compactification `P1 x P1 -> P2`

Use homogeneous coordinates `[X0:X1]` and `[Y0:Y1]` on
`P1 x P1`, and set


\[
D=X_0Y_0+X_1Y_1.
\]

The three sections


\[
\begin{aligned}
s_0&=X_0^2D,\\
s_1&=X_0^3Y_1+X_1^2D,\\
s_2&=X_1^3Y_0
\end{aligned}
\tag{2.1}
\]

all have bidegree `(3,1)`.

They have no common zero.  On `X0=0`,


\[
(s_1,s_2)=X_1^3(Y_1,Y_0),
\]

which never vanishes simultaneously.  Away from `X0=0`, a zero of `s0`
has `D=0`; parametrizing that divisor by


\[
(Y_0,Y_1)=(-X_1,X_0)
\]

gives


\[
(s_1,s_2)=(X_0^4,-X_1^4),
\]

again without a common zero.  Therefore


\[
\phi=[s_0:s_1:s_2]:\mathbb P^1\times\mathbb P^1\to\mathbb P^2
\tag{2.2}
\]

is a morphism.

The line bundle `O(3,1)` is ample, so `phi` contracts no curve.  Its image is
a closed surface in `P2`, hence all of `P2`; consequently `phi` is finite.
If `A` and `B` are the two ruling classes, then


\[
(3A+B)^2=6,
\]

so `phi` has degree six.

The affine chart `X0 != 0, D != 0` is `A2` with


\[
x=X_1/X_0,
\qquad
z=Y_1/D.
\]

After scaling `X0=D=1`, one has `(Y0,Y1)=(1-xz,z)`, and (2.1) becomes


\[
[s_0:s_1:s_2]=[1:x^2+z:x^3-x^4z].
\]

Thus (2.2) really extends the affine map (1.1), rather than merely sharing
its degree.

## 3. The sole target boundary is globally split

Let the target boundary be the line


\[
T=\{U_0=0\}\subset\mathbb P^2.
\]

Equation `s0=X0^2D` gives the complete pullback


\[
\phi^*T=2C_1+C_2,
\tag{3.1}
\]

where


\[
C_1=\{X_0=0\}\sim A,
\qquad
C_2=\{D=0\}\sim A+B.
\]

Their intersection matrix is


\[
M=
\begin{pmatrix}
0&1\\
1&2
\end{pmatrix}.
\tag{3.2}
\]

Both components dominate `T`.  Their tangential degrees are


\[
m_i=(3A+B)\mathbin\cdot C_i,
\qquad
m=(1,4),
\]

while the normal vector from (3.1) is


\[
n=(2,1).
\]

The full raw-intersection transport is exact:


\[
Mn=m,
\qquad
n^tm=6,
\qquad
n^tMn=6=T^2\deg\phi.
\tag{3.3}
\]

Because the target boundary has only the one component `T`, and its full
preimage has the two noncontracted components `C1,C2`, condition `(*)` fails
globally for this compactified cover.  There is no second target component
that could serve as a unique-preimage witness.

## 4. Exact geometric monodromy

Over `C(P,Q)`, the six sheets are the roots of


\[
f(X)=X^6-PX^4+X^3-Q.
\tag{4.1}
\]

### 4.1 The polynomial is indecomposable

Only degree patterns `2 o 3` and `3 o 2` are possible.

For a normalized `2 o 3` decomposition, write


\[
v=X^3+aX^2+bX,
\qquad
u=v^2+cv.
\]

Coefficient comparison with `X^6-PX^4+X^3` gives successively


\[
a=0,
\qquad
b=-P/2,
\qquad
c=1.
\]

The remaining residual is


\[
\frac{P^2}{4}X^2-\frac P2X\ne0
\]

over `C(P)`.  Thus no such decomposition exists.

For a normalized `3 o 2` decomposition, the inner polynomial is
`X^2+aX`.  The `X^5` coefficient forces `a=0`, after which the composition
has only even powers, contradicting the `X^3` term.  Therefore (4.1) is
indecomposable, and its polynomial monodromy is primitive.

### 4.2 A simple branch factor supplies a transposition

Exact discriminant expansion gives


\[
\operatorname{Disc}_X(f)=Q^2H(P,Q),
\tag{4.2}
\]

where


\[
\begin{aligned}
H={}&1024P^6Q+13824P^3Q^2-8640P^3Q-108P^3\\
  &+46656Q^3+34992Q^2+8748Q+729.
\end{aligned}
\]

The checker proves


\[
\gcd_{\mathbb C(P)[Q]}(H,\partial_QH)=1.
\]

Thus a generic component of `H=0` occurs to discriminant order one and has
inertia a transposition.  A primitive permutation group containing a
transposition is the full symmetric group.  Hence


\[
\boxed{G_{\mathrm{geom}}=S_6.}
\tag{4.3}
\]

## 5. Even the numerical `2--1` corridor can coexist

Blow up first the crossing `C1 intersection C2` and then a free point of the strict
transform of `C1`.  In the order


```text
(E, C1, Z, C2),
```

the boundary intersection matrix and augmented-canonical labels are


\[
M'=
\begin{pmatrix}
-1&1&0&0\\
1&-2&1&0\\
0&1&-1&1\\
0&0&1&1
\end{pmatrix},
\qquad
(2,1,0,-1).
\tag{5.1}
\]

The total transform of `T` and the tangential/contracted-degree vector are


\[
n'=(2,2,3,1),
\qquad
m'=(0,1,0,4).
\]

Again,


\[
M'n'=m',
\qquad
(n')^tm'=6,
\qquad
(n')^tM'n'=6.
\tag{5.2}
\]

Moreover, `det(M')=-1`, and the negative-intersection determinants obtained
by deleting `E` and `C1` are


\[
d_E=-3,
\qquad
d_{C_1}=-2,
\qquad
d_E=d_{C_1}-1.
\]

So the exact numerical label corridor


```text
2 -- 1 -- 0 -- (-1)
```

and the determinant parity left by the one-dicritical `S6` reduction are
arithmetically consistent with the split cover.

This is deliberately **not** a type assignment.  The new curve `E` is
contracted to a point of the target boundary; it is not a dicritical mapping
onto the affine branch.  The calculation blocks a numerical graph argument,
not a Keller-specific geometric one.

## 6. Reproducible certificate and exact boundary

Run


```bash
uv run python -m scripts.s6_split_cover
uv run pytest -q tests/test_s6_split_cover.py
```

[`s6_split_cover.py`](../scripts/s6_split_cover.py) checks the monic algebra
presentation, affine inverse coordinate, Jacobian, all basepoint-free chart
identities, both intersection-transport systems, discriminant factorization,
squarefreeness, and the two decomposition obstructions.

The construction proves that the following package is insufficient to force
condition `(*)`:

- an actual finite flat degree-six cover of `A2`;
- an actual finite compactification over `P2`;
- geometric monodromy `S6`;
- two-component boundary pullback with exact intersection transport; and
- the surviving numerical canonical/determinant corridor.

What it omits is exactly what a successful argument must use: affine
unramifiedness, the Keller log-canonical formula with a genuine type-three
dicritical, and the forced finite branch passport.

## Source context

Orevkov's condition `(*)` and its role in the two-characteristic-pair splice
argument are in S. Yu. Orevkov,
[*Counter-examples to the Jacobian Conjecture at
Infinity*](https://www.math.univ-toulouse.fr/~orevkov/jci-e.pdf), Section 2.4.
The explicit cover and all calculations in this note are derived and checked
here.
