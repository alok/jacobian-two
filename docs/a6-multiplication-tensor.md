# The missing multiplication tensor in the hostile `A6` lattice

## Claim boundary

The one-dicritical `A6` trace-lattice package has an explicit hostile
quadratic model

\[
T=\Phi^tH\Phi,
\qquad
\det T=-b^2.
\]

The purpose of this note is to ask the next necessary question: can this
quadratic lattice be enhanced to a unital commutative associative rank-six
algebra whose **ordinary** regular trace form is `T`?

There is no global answer here.  Instead, this note provides three exact
outputs.

1. It writes a necessary-and-sufficient polynomial system for the missing
   cubic multiplication tensor.
2. It reduces every apparent division by `b` to one explicit identity on the
   normalization.
3. It constructs the required associative algebras at the forced `5+1` cusp
   and `3+3` collision fibers, compatibly with the local `Phi/H` factorization.

The third item closes a tempting route: neither forced singular fiber gives a
pointwise multiplication contradiction.  The constructions are not a global
rank-six algebra, do not have global `A6` monodromy, and are not Keller maps.
The exact checker is
[`a6_local_multiplication.py`](../scripts/a6_local_multiplication.py), with
regressions in
[`test_a6_local_multiplication.py`](../tests/test_a6_local_multiplication.py).

## 1. The hostile quadratic data

Put

\[
R=\mathbb C[P,Q]
\]

and

\[
b=-P^5+5P^3Q-5PQ^2+Q^3+Q^2.
\tag{1.1}
\]

The symmetric normalization presentation is

\[
K=
\begin{pmatrix}
-P^2&P^2-Q&P-Q\\
P^2-Q&-Q&-P\\
P-Q&-P&P-1
\end{pmatrix},
\qquad
\det K=b.
\tag{1.2}
\]

In the ordered basis

\[
(e_1,e_2,e_3,f_1,f_2,f_3),
\]

set

\[
\Phi=\operatorname{diag}(K,I_3),
\qquad
H=
\begin{pmatrix}
0&I_3\\
I_3&0
\end{pmatrix}.
\tag{1.3}
\]

Then

\[
T=\Phi^tH\Phi
=
\begin{pmatrix}
0&K\\
K&0
\end{pmatrix},
\qquad
\det T=-b^2.
\tag{1.4}
\]

The trace-unit candidate is

\[
u=e_3-3f_2-3f_3,
\qquad
u^tTu=6.
\tag{1.5}
\]

Write `A=adj(K)`.  Exact expansion gives

\[
A=
\begin{pmatrix}
-P^2-PQ+Q&-P^3+2PQ-Q&-P^3+2PQ-Q^2\\
-P^3+2PQ-Q&-P^3+2PQ-Q^2&-Q(P^2+P-Q)\\
-P^3+2PQ-Q^2&-Q(P^2+P-Q)&-P^4+3P^2Q-Q^2
\end{pmatrix}
\tag{1.6}
\]

and the checker verifies

\[
KA=AK=bI_3.
\tag{1.7}
\]

## 2. A cubic tensor with only 56 unknown entries

Let

\[
C=(C_{ijk})\in\operatorname{Sym}^3((R^6)^*)
\tag{2.1}
\]

be a completely symmetric cubic tensor.  There are `56` independent
polynomial entries.  For every pair `i,j`, split its final-index slice into

\[
v_{ij}=
\begin{pmatrix}
C_{ij1}\\C_{ij2}\\C_{ij3}
\end{pmatrix},
\qquad
w_{ij}=
\begin{pmatrix}
C_{ij\bar1}\\C_{ij\bar2}\\C_{ij\bar3}
\end{pmatrix}.
\tag{2.2}
\]

Over the fraction field, the putative product is forced by `T`:

\[
e_i e_j=
\begin{pmatrix}
x_{ij}\\y_{ij}
\end{pmatrix},
\qquad
x_{ij}=\frac{Aw_{ij}}b,
\qquad
y_{ij}=\frac{Av_{ij}}b.
\tag{2.3}
\]

Indeed, lowering the output index by `T` gives

\[
C_{ij,*}=T(e_i e_j,-)
=
\begin{pmatrix}
Ky_{ij}\\Kx_{ij}
\end{pmatrix}.
\]

The following equations are necessary and sufficient for (2.3) to be a
polynomial, unital, commutative, associative algebra law with ordinary trace
matrix exactly `T`.

### 2.1 Polynomial lift

For all `i <= j`,

\[
Av_{ij}\in bR^3,
\qquad
Aw_{ij}\in bR^3.
\tag{2.4}
\]

### 2.2 Unit

For every `j,k`,

\[
\sum_i u_iC_{ijk}=T_{jk}.
\tag{2.5}
\]

Because `T` is nondegenerate away from `b=0`, this says that multiplication by
`u` is the identity over the fraction field.  Polynomiality then makes it an
identity over `R`.

### 2.3 Associativity

Define

\[
E_{ij,kl}
=v_{ij}^tAw_{kl}+w_{ij}^tAv_{kl}.
\tag{2.6}
\]

The denominator-cleared WDVV equations are

\[
E_{ij,kl}=E_{jk,il}
\tag{2.7}
\]

for all four indices.  Pairing an associator with the fourth basis vector
shows that these equations are precisely associativity over `Frac(R)`.  The
generic nondegeneracy of `T` and polynomiality extend the equality over `R`.

### 2.4 The Frobenius functional must be the ordinary trace

Symmetry of `C` only produces an invariant Frobenius functional.  It need not
be the ordinary trace of the regular representation.  The latter condition is

\[
\sum_{a,r=1}^3
A_{ar}
\left(C_{i a\bar r}+C_{i\bar a r}\right)
=b\tau_i,
\tag{2.8}
\]

where

\[
\tau=Tu=
\begin{pmatrix}
-3(P^2+P-2Q)\\
3(P+Q)\\
3\\
P-Q\\
-P\\
P-1
\end{pmatrix}.
\tag{2.9}
\]

Equations (2.4), (2.5), (2.7), and (2.8) are therefore an exact finite
polynomial presentation of the trace-compatible algebra problem.  One small
degree conclusion is immediate: in this fixed basis, the multiplication
coefficients cannot all have degree at most one, because the first regular
trace in (2.9) contains `-3*P^2`.  There is no geometrically justified upper
degree bound, so failure of a bounded ansatz would not prove nonexistence.

## 3. Divisibility is a normalization identity

The normalization is

\[
\nu:\mathbb A^1_t\longrightarrow B,
\qquad
(P,Q)=(t^2+t^3,t^5).
\tag{3.1}
\]

Set

\[
q(t)=
\begin{pmatrix}
1\\t\\t^2
\end{pmatrix},
\qquad
\Phi_5(t)=t^4+t^3+t^2+t+1.
\tag{3.2}
\]

Exact substitution gives

\[
K(t^2+t^3,t^5)q(t)=0
\tag{3.3}
\]

and the stronger factorization

\[
\operatorname{adj}K(t^2+t^3,t^5)
=-t^4\Phi_5(t)q(t)q(t)^t.
\tag{3.4}
\]

The quotient map in the normalization presentation is

\[
\psi:R^3\longrightarrow\mathbb C[t],
\]

\[
\psi(v)=
v_1(t^2+t^3,t^5)
+t v_2(t^2+t^3,t^5)
+t^2v_3(t^2+t^3,t^5).
\tag{3.5}
\]

Consequently,

\[
Av\in bR^3
\iff v\in\operatorname{im}K
\iff\psi(v)=0.
\tag{3.6}
\]

The first equivalence also follows directly from (1.7): if `Av=bw`, then
`b(v-Kw)=0`, hence `v=Kw` because `R` is a domain.  Conversely, `v=Kw`
implies `Av=bw`.

Thus each three-component divisibility constraint in (2.4) becomes one exact
polynomial identity in `t`.  This is useful for a bounded symbolic search:
substitution and coefficient comparison turn all lift constraints into linear
equations on the coefficients of `C`.

### 3.1 What happens at the cusp

At `t=0`, equation (3.5) and its first derivative give

\[
v_1(0,0)=0,
\qquad
v_2(0,0)=0.
\tag{3.7}
\]

Therefore

\[
v(0,0)\in\mathbb C(0,0,1)^t
=\operatorname{im}K(0,0).
\]

The next normalization-jet equation is

\[
\partial_Pv_1(0,0)+v_3(0,0)=0.
\tag{3.8}
\]

This explains why merely substituting the cusp into `adj(K)*v` is too weak:
the factor `t^4` in (3.4) makes the adjugate vanish there, while the value and
jet equations retain the missing module information.

### 3.2 What happens at a collision

The two collision targets are

\[
(P,Q)=(\rho,1),
\qquad
\rho^2+\rho-1=0.
\tag{3.9}
\]

Each target has two distinct primitive-fifth-root preimages on the
normalization.  For those two values `t_1,t_2`, equation (3.5) gives

\[
q(t_1)^tv(\rho,1)=0,
\qquad
q(t_2)^tv(\rho,1)=0.
\tag{3.10}
\]

The equations are independent and leave precisely the one-dimensional image
of `K(rho,1)`.  Here the `Phi_5(t)` factor in (3.4) records the two-branch
normalization collision.

## 4. The cusp admits the exact `5+1` algebra

At `(P,Q)=(0,0)`,

\[
K_0=
\begin{pmatrix}
0&0&0\\
0&0&0\\
0&0&-1
\end{pmatrix}.
\tag{4.1}
\]

Choose `delta` with `delta^2=-5`, and put

\[
a=\frac{5+\delta}{6},
\qquad
c=\frac{-5+\delta}{2},
\]

\[
p=ae_3+cf_3,
\qquad
q=u-p.
\tag{4.2}
\]

The ordered basis matrix

\[
L_{\mathrm{cusp}}
=[p,f_1,f_2,e_1,e_2,q]
\tag{4.3}
\]

has determinant `-delta`, and the checker proves

\[
L_{\mathrm{cusp}}^tT_0L_{\mathrm{cusp}}
=\operatorname{diag}(5,0,0,0,0,1).
\tag{4.4}
\]

Identify this basis with

\[
(1_5,z,z^2,z^3,z^4,1_1)
\]

in

\[
A_{\mathrm{cusp}}
=\mathbb C[z]/(z^5)\times\mathbb C.
\tag{4.5}
\]

Its unit maps to `p+q=u`, and its ordinary trace form is exactly (4.4).  The
only nonzero entries of its trace cubic in this canonical basis are

\[
C(p,p,p)=5,
\qquad
C(q,q,q)=1.
\]

This is exactly the multiplication type of a `5+1` fiber.

### 4.1 Compatibility with the middle lattice

Take the divisor section

\[
s=(z^2,1).
\tag{4.6}
\]

Multiplication by `s` has kernel `span(z^3,z^4)`, which (4.3) maps to

\[
\operatorname{span}(e_1,e_2)=\ker\Phi_0.
\]

Define a perfect Frobenius functional by

\[
\lambda_{\mathrm{cusp}}(a,r)=5[z^4]a+r.
\tag{4.7}
\]

Then

\[
\lambda_{\mathrm{cusp}}(s^2ab)
=\operatorname{Tr}_{A_{\mathrm{cusp}}/\mathbb C}(ab).
\tag{4.8}
\]

The checker constructs the Gram matrix of (4.7), checks that it is
nondegenerate, checks (4.8) as a matrix identity, and checks the kernel
alignment with `Phi_0`.

## 5. A collision admits the exact `3+3` algebra

Work over the collision residue field with

\[
\rho^2+\rho-1=0,
\qquad Q=1.
\]

Define

\[
r_1=e_1-\rho e_2,
\qquad
r_2=-e_1+e_3,
\]

\[
s_1=f_1-\rho f_2,
\qquad
s_2=-f_1+f_3.
\tag{5.1}
\]

The first pair uses `ker K(rho,1)` in the first block, and the second pair
uses the corresponding copy of `ker K(rho,1)` in the second block.  Together
they span the radical of `T_rho`; only `span(r_1,r_2)` is the kernel of
`Phi_rho`.  Put

\[
x=e_3-\frac12u,
\]

\[
p_+=\frac12u+ix,
\qquad
p_-=\frac12u-ix.
\tag{5.2}
\]

The exact identities are

\[
T_\rho(p_+,p_+)=3,
\qquad
T_\rho(p_-,p_-)=3,
\qquad
T_\rho(p_+,p_-)=0.
\]

Use the ordered basis

\[
L_{\mathrm{node}}
=[p_+,s_1,r_1,p_-,s_2,r_2].
\tag{5.3}
\]

Modulo `rho^2+rho-1`, its determinant is `3*i`, and

\[
L_{\mathrm{node}}^tT_\rho L_{\mathrm{node}}
=\operatorname{diag}(3,0,0,3,0,0).
\tag{5.4}
\]

Identify (5.3) with

\[
(1_z,z,z^2,1_w,w,w^2)
\]

in

\[
A_{\mathrm{node}}
=\mathbb C[z]/(z^3)\times\mathbb C[w]/(w^3).
\tag{5.5}
\]

Its unit maps to `p_++p_-=u`, and its ordinary trace form is (5.4).  Thus the
exact specialized quadratic data admit the required `3+3` multiplication.

For the middle lattice, use

\[
s=(z,w)
\tag{5.6}
\]

and

\[
\lambda_{\mathrm{node}}(a,c)
=3[z^2]a+3[w^2]c.
\tag{5.7}
\]

The kernel of multiplication by `s` is `span(z^2,w^2)`, which maps to
`span(r_1,r_2)=ker(Phi_rho)`, and

\[
\lambda_{\mathrm{node}}(s^2ab)
=\operatorname{Tr}_{A_{\mathrm{node}}/\mathbb C}(ab).
\tag{5.8}
\]

Again, the checker verifies the perfect middle Gram matrix, its trace
factorization, and the kernel alignment exactly.

The induced isometry from the image of multiplication by `s` to the image of
`Phi` extends to the full nondegenerate middle spaces by Witt extension.
Transporting the regular module action through such an extension supplies the
local `H`-self-adjoint action matrices.  Therefore even the full local
`Phi/H` package, rather than only `T`, has no pointwise contradiction.

## 6. The additional global middle-lattice equations

An algebra with trace matrix `T` does not automatically preserve the specific
middle lattice selected by `Phi`.  Let `L_i` be multiplication by `e_i` in the
source basis, written in `3+3` blocks as

\[
L_i=
\begin{pmatrix}
X_i&Y_i\\
Z_i&W_i
\end{pmatrix}.
\]

Over `R[b^-1]`, its induced action on the middle lattice must be

\[
N_i=\Phi L_i\Phi^{-1}
=
\begin{pmatrix}
KX_iA/b&KY_i\\
Z_iA/b&W_i
\end{pmatrix}.
\tag{6.1}
\]

Polynomial preservation requires

\[
b\mid KX_iA,
\qquad
b\mid Z_iA.
\tag{6.2}
\]

The `H`-balanced condition `N_i^tH=HN_i` becomes

\[
Z_iA/b\text{ is symmetric},
\qquad
KY_i\text{ is symmetric},
\tag{6.3}
\]

\[
W_i=(KX_iA/b)^t.
\tag{6.4}
\]

These equations must be coupled to commutativity, associativity, the unit,
and the regular trace.  If

\[
w=\Phi u=(P-Q,-P,P-1,0,-3,-3)^t,
\]

then the cyclic-section equations are

\[
N_iw=\Phi e_i.
\tag{6.5}
\]

Equations (2.4)--(2.8) and (6.1)--(6.5) are the exact global polynomial target
left by this calculation.

## 7. A hostile family along the normalization

The required fiber partitions are also compatible in one free associative
family over the normalization parameter.  Put

\[
g(t)=\Phi_5(t),
\]

\[
a(t)=1-g(t),
\qquad
c(t)=1-(1+t)g(t),
\]

and define

\[
\mathcal B
=\mathbb C[t,z]\big/
\left(z^3(z-a(t))(z-c(t))(z-1)\right).
\tag{7.1}
\]

The defining polynomial is monic of degree six, so `B` is free of rank six
over `C[t]`.  At `t=0`,

\[
z^3(z-a)(z-c)(z-1)=z^5(z-1),
\tag{7.2}
\]

which is the `5+1` partition.  Modulo `Phi_5(t)`,

\[
z^3(z-a)(z-c)(z-1)=z^3(z-1)^3,
\tag{7.3}
\]

which is the `3+3` partition at every primitive fifth root.  Generically the
displayed roots give `3+1+1+1`.

This family may have additional accidental collision parameters.  More
importantly, it is only a family over `C[t]`; no descent through the singular
branch, equality with the global trace matrix, connected `A6` monodromy, or
Keller realization is claimed.

## 8. What remains open

The calculation rules out a local rank or multiplication contradiction at the
two forced events.  It does not decide whether the global systems in Sections
2 and 6 have a solution.  A successful next obstruction must retain how the
normalization jets, the middle `A`-module lattice, connected monodromy, and the
compactification at infinity couple globally.

In particular:

- the cusp algebra is a residue-fiber model, not a global degree-six cover;
- the collision algebra is a residue-fiber model, not a descent across both
  normalization preimages;
- the family (7.1) is a normalization-level stopping fixture;
- no polynomial solution of the complete cubic system was constructed;
- no contradiction for the complete cubic system was proved;
- the two-dimensional Jacobian conjecture and its degree-six frontier remain
  open.
