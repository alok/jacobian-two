# The forced canonical-label corridor beside a one-dicritical leaf

## Claim boundary

Work on Orevkov's original smooth resolution

\[
Z\longrightarrow\mathbb P^2
\]

obtained only by source blowups, with boundary

\[
L=Z\setminus\mathbb A^2.
\]

For a boundary curve \(C\), write \(\ell(C)\) for its coefficient in the
augmented canonical divisor

\[
\overline K_Z=K_Z+L.
\]

Assume there is one dicritical component \(E\), its normal ramification index
is \(e\ge 2\), and the relatively minimal constant chain has been eliminated
as in
[`one-dicritical-source-smoothness.md`](one-dicritical-source-smoothness.md).
This note proves:

\[
\boxed{
\begin{gathered}
E\text{ is a type-3 leaf},\qquad \ell(E)=e,\\
\ell(A)=e(-E^2)-1>0,\qquad
\ell(A)\equiv-1\pmod e,
\end{gathered}
}
\tag{0.1}
\]

where \(A\) is the unique boundary neighbor of \(E\).  Moreover, \(A\) is
type 2, and the path from \(A\) to the strict transform of the original line
at infinity contains consecutive labels

\[
\boxed{1\;--\;0\;--\;(-1).} \tag{0.2}
\]

For the two surviving one-dicritical passports, the neighbor labels therefore
lie in the exact progressions

\[
S_6:\quad 1,3,5,\ldots,
\qquad
A_6:\quad 2,5,8,\ldots. \tag{0.3}
\]

This is a compactification-graph restriction, not an elimination of either
passport.  Leafhood and the displayed self-intersection equation refer to
Orevkov's original pure-blowup resolution.  Later target blowups and
re-resolution can give the strict transform of \(E\) more neighbors and change
its self-intersection.

## 1. Why the dicritical is a leaf

Orevkov's Lemma 2.1 says that every connected component of \(L_{FC}\) is a
linear chain, meets \(L_\infty\) once at an endpoint, and has the form

\[
L_C-E
\]

with one dicritical endpoint \(E\).  In the relatively minimal model the cyclic
endpoint theorem gives \(L_C=\varnothing\).  The finite component is therefore
\(E\) alone.  It meets the rest of the full boundary at exactly one point, so
\(E\) is a leaf and its unique neighbor \(A\) lies in \(L_\infty\).

The original line at infinity belongs to neither of the finite branches in
Orevkov's proof.  Thus \(E\) is exceptional for the source blowup map, and

\[
m=-E^2\ge1. \tag{1.1}
\]

## 2. Keller adjunction fixes the adjacent label

Borisov's logarithmic Keller formula is

\[
\overline K_Z
=\rho^*\overline K_Y
+\sum_{\operatorname{type}(C)=3}e_C C. \tag{2.1}
\]

The coefficient of a type-3 curve is its normal ramification index \(e_C\), not
\(e_C-1\).  Therefore

\[
\ell(E)=e. \tag{2.2}
\]

Adjunction on the smooth rational leaf gives

\[
(K_Z+E)\mathbin\cdot E=-2.
\]

Because the reduced boundary \(L-E\) meets \(E\) once,

\[
\overline K_Z\mathbin\cdot E=-1. \tag{2.3}
\]

Only \(E\) and its unique neighbor \(A\) contribute to this intersection.  If
\(a=\ell(A)\), equations (2.2)--(2.3) give

\[
eE^2+a=-1,
\qquad
\boxed{a=em-1.} \tag{2.4}
\]

Now \(e\ge2\) and \(m\ge1\), so

\[
a\ge e-1>0,
\qquad
a\equiv-1\pmod e,
\qquad
\gcd(a,e)=1. \tag{2.5}
\]

For \(e=2\), (2.4) gives all positive odd labels.  For \(e=3\), it gives the
positive labels congruent to two modulo three, proving (0.3).

## 3. The neighbor is type 2

Since \(A\) lies in \(L_\infty\), it is type 1 or type 2 in Borisov's
terminology.  The target of Orevkov's original construction is
\(\mathbb P^1\times\mathbb P^1\), whose boundary is the union of two fibers.
Hence

\[
K_Y+(Y\setminus\mathbb A^2)=-F_1-F_2, \tag{3.1}
\]

so both target boundary curves have augmented-canonical label \(-1\).
Borisov's type-1 pullback rule says that a type-1 source label is its positive
normal degree times the target label.  Such a label is negative, while (2.5)
makes \(\ell(A)\) positive.  Therefore

\[
\boxed{A\text{ is type 2}.} \tag{3.2}
\]

The same sign argument shows that every positive-label vertex on the segment
from \(A\) toward the negative core is type 2, as long as the graph remains in
this original target model.

## 4. A forced `1--0--(-1)` transition

The strict transform of the original line at infinity has label \(-2\).
Borisov's canonical-label tree theorem says:

1. the negative-label induced subgraph is connected;
2. positive vertices are separated from it by zero-label vertices; and
3. a zero-label vertex is adjacent only to labels `1` and `-1`.

The full boundary graph is a tree.  Its unique path from the positive vertex
\(A\) to the negative original-infinity vertex must therefore cross consecutive
labels

\[
1\;--\;0\;--\;(-1),
\]

proving (0.2).  This is an existence statement along the path.  It does not
say that \(A\) itself has label one, and it does not assert that the transition
is unique.

## 5. Exact executable endpoint

The typed object
[`DicriticalLeafLabelEndpoint`](../scripts/six_sheet_monodromy.py) records the
arithmetic part of (2.4).  Its regression test checks the first six possible
self-intersections in both passports:

\[
e=2:\quad(1,3,5,7,9,11),
\qquad
e=3:\quad(2,5,8,11,14,17).
\]

As with the cyclic endpoint certificate, this does not formalize the
geometric input: Orevkov's leaf theorem, Borisov's log-Keller coefficient, or
the type/path assertions.

## Model scope and next target

Equation (2.4) is asserted before any later target blowups.  If a Borisov
framework is produced by blowing up the target at the image of \(E\) at
infinity and resolving again, \(\ell(E)=e\) remains a valuation invariant, but leafhood,
\(E^2\), and the number of adjacent curves must be transformed rather than
silently reused.

Within the original model, every surviving one-dicritical graph must begin
with the corridor

\[
E_{\mathrm{type\ 3},\,\ell=e}
\;--\;
A_{\mathrm{type\ 2},\,\ell=em-1}
\;--\;\cdots\;--\;
1\;--\;0\;--\;(-1)\;--\;\cdots\;--\;(-2).
\]

The remaining graph attack must combine this corridor with determinant
labels, the \(A_6\) cusp/collision data or the exact \(S_6\) jump partition,
and the degree-six pushforward/pullback equations.  The canonical labels alone are
not yet contradictory and \(m\) is not bounded by this argument.

## Primary sources

- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587--596.  Lemma 2.1 gives
  the finite linear chain and its unique meeting with the infinity boundary.
- Alexander Borisov,
  [“Frameworks for two-dimensional Keller maps”](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v27i3p54/pdf/),
  *Electronic Journal of Combinatorics* 27(3) (2020), P3.54.  The log-Keller
  formula gives the type-3 coefficient `e`; the boundary adjunction
  calculation and type-1 pullback rule give (2.4) and (3.2).
- Alexander Borisov,
  [“On two invariants of divisorial valuations at infinity”](https://people.math.binghamton.edu/borisov/documents/papers/divisorialvaluations.pdf).
  Proposition 2.2 supplies connectivity of the negative region and the
  `1--0--(-1)` transition rule.

The leaf-neighbor congruence and its application to the two surviving
passports are derived in this repository.  No claim of historical priority is
made.
