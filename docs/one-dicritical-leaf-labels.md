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
On every such pure source-blowup resolution, this note proves:

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

There is a safe further minimization.  Call the resolved triple jointly
minimal if no boundary curve is at once exceptional for the source blowdown,
contracted by the extended map, and a \((-1)\)-curve.  Every candidate admits
such a model, and on it

\[
\boxed{E^2=-1,\qquad \ell(A)=e-1.} \tag{0.3}
\]

Thus the two minimal corridors begin

\[
S_6:\quad 2\;--\;1,
\qquad
A_6:\quad 3\;--\;2. \tag{0.4}
\]

This is a compactification-graph restriction, not an elimination of either
passport.  The general labels \(em-1\) for \(m>1\) are precisely the orbit
created by unnecessary corner blowups along the leaf edge; they are not
distinct minimal possibilities.  Leafhood and the displayed
self-intersection equations refer to a pure source-blowup resolution.  Later
target blowups and re-resolution can give the strict transform of \(E\) more
neighbors and change its self-intersection.

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
positive labels congruent to two modulo three.  Section 3 explains why joint
minimality selects the first member of each progression.

## 3. Joint minimization forces \(E^2=-1\)

Write \(\pi:Z\to\mathbb P^2\) for the source blowdown and
\(F:Z\to\mathbb P^1\times\mathbb P^1\) for the resolved polynomial map.
Suppose a boundary curve \(C\) is simultaneously \(\pi\)-exceptional,
\(F\)-contracted, and satisfies \(C^2=-1\).  Its Castelnuovo contraction is a
smooth point and can be taken as the next step in a boundary-blowdown
factorization; the reduced boundary remains simple normal crossings.  Both
\(\pi\) and \(F\) are constant on \(C\), so they descend across the contraction.

Repeating lowers the Picard rank and terminates.  It changes neither affine
plane, the target, the generic degree, the branch passport, nor the number of
dicriticals.  It also cannot create a constant component over the affine
target.  Thus joint minimality is a without-loss-of-generality condition, not
an additional hypothesis.

There is a short intrinsic proof on the jointly minimal model.  Choose the
chronologically last-created boundary component \(P\) with positive
augmented-canonical label.  Such a component exists because \(\ell(E)=e>0\).
Borisov's label-tree theorem excludes an edge directly joining a positive
vertex to a negative vertex, so every neighbor of \(P\) has nonnegative label.
If a later blowup touched \(P\), a free blowup would create label
\(\ell(P)+1>0\), while a corner blowup would create the sum of two labels, again
positive.  Either contradicts the choice of \(P\).  Thus \(P\) is untouched
after its creation and \(P^2=-1\).

The curve \(P\) is exceptional for \(\pi\), because the original infinity line
has label \(-2\).  It cannot be type 1, whose label is negative on this target;
it cannot be type 2 by joint minimality; and it cannot be type 4 because
\(L_C=\varnothing\).  Hence it is type 3.  The dicritical is unique, so \(P=E\)
and \(E^2=-1\).  Leaf adjunction (2.4) now gives

\[
\boxed{E^2=-1,\qquad \ell(A)=e-1,} \tag{3.1}
\]

which proves (0.3)--(0.4).  Conversely, a corner blowup at \(E\cap A\)
changes

\[
m\longmapsto m+1,
\qquad
em-1\longmapsto e(m+1)-1. \tag{3.2}
\]

This recovers the nonminimal progressions and shows exactly why they carry no
extra graph information.  Equivalently, leafhood says that \(E\) was born by
a free blowup from a parent of label \(e-1\); joint minimality prevents any
later corner-blowup cluster from remaining above that edge.

## 4. The neighbor is type 2

Since \(A\) lies in \(L_\infty\), it is type 1 or type 2 in Borisov's
terminology.  The target of Orevkov's original construction is
\(\mathbb P^1\times\mathbb P^1\), whose boundary is the union of two fibers.
Hence

\[
K_Y+(Y\setminus\mathbb A^2)=-F_1-F_2, \tag{4.1}
\]

so both target boundary curves have augmented-canonical label \(-1\).
Borisov's type-1 pullback rule says that a type-1 source label is its positive
normal degree times the target label.  Such a label is negative, while (2.5)
makes \(\ell(A)\) positive.  Therefore

\[
\boxed{A\text{ is type 2}.} \tag{4.2}
\]

The same sign argument shows that every positive-label vertex on the segment
from \(A\) toward the negative core is type 2, as long as the graph remains in
this original target model.

## 5. A forced `1--0--(-1)` transition

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
assert that the transition is unique.  In the jointly minimal \(S_6\) model,
\(A\) does itself have label one; in the \(A_6\) model it has label two.

## 6. Determinant labels at the minimal leaf

For a boundary component \(C\), let \(d_C\) be the determinant of the negative
intersection matrix after deleting \(C\).  For an edge \(CD\), write
\(d_{CD}\) for the corresponding edge-deletion determinant.  The determinant
of the full negative boundary matrix is \(-1\): the boundary classes form a
basis of \(\operatorname{Pic}(Z)\), and the intersection form is congruent to
\(\operatorname{diag}(1,-1,\ldots,-1)\).

The minimal leaf \(E\) was created by a free blowup on \(A\), and its edge was
never touched.  Borisov's determinant blowup formula therefore gives

\[
\boxed{d_E=d_{EA}=d_A-1.} \tag{6.1}
\]

Moreover \(d_E<0\).  Indeed, for the ample target boundary
\(H=F_1+F_2\), the divisor \(F^*H\) is supported on \(L-E\) and has positive
square \(6H^2=12\).  The Hodge index theorem makes the determinant of the
negative intersection form on \(L-E\) negative.  Finally, Borisov's content
\(c_C=\ell(C)+d_C\) is always odd.  Consequently

\[
\begin{array}{c|c|c}
&d_E&d_A\\ \hline
S_6&-1,-3,-5,\ldots&0,-2,-4,\ldots\\
A_6&-2,-4,-6,\ldots&-1,-3,-5,\ldots
\end{array} \tag{6.2}
\]

These restrictions still do not make the graph list finite.  There is an
explicit hostile family for every \(n\ge0\): start from the original infinity
line of label \(-2\), use free blowups to create labels \(-1,0,1\), and then
blow up the edge between the zero vertex and the nearest label-one vertex
\(n\) times.  Attach the minimal \(S_6\) edge by one final free blowup, or the
minimal \(A_6\) edge by two.  Every positive nondicritical vertex has
self-intersection at most \(-2\), while

\[
\begin{aligned}
S_6:&\quad d_E=-(n+1)(n+2)-1,\\
A_6:&\quad d_E=-(n+1)(n+2)-2.
\end{aligned} \tag{6.3}
\]

Thus canonical labels, self-intersections, determinant labels, and joint
minimality alone admit infinitely many abstract trees.  These trees are not
realized degree-six maps: they do not provide pullback matrices, type-1
covering degrees, Belyi maps, or Keller coordinates.

## 7. Exact executable endpoint

The typed object
[`DicriticalLeafLabelEndpoint`](../scripts/six_sheet_monodromy.py) records the
arithmetic part of (2.4).  Its regression test starts with the two jointly
minimal endpoints and applies five corner blowups, recovering

\[
e=2:\quad(1,3,5,7,9,11),
\qquad
e=3:\quad(2,5,8,11,14,17).
\]

As with the cyclic endpoint certificate, this does not formalize the
geometric input: Orevkov's leaf theorem, the joint-contraction argument,
Borisov's log-Keller coefficient, or the type/path assertions.

The separate exact checker
[`canonical_leaf_graph.py`](../scripts/canonical_leaf_graph.py) performs the
blowups in the infinite hostile family, rebuilds every intersection matrix,
and verifies (6.1)--(6.3) for arbitrary tested \(n\).  It certifies only the
stated graph arithmetic and deliberately assigns, rather than realizes, the
geometric component types.

## Model scope and next target

Equation (2.4) is asserted before any later target blowups.  If a Borisov
framework is produced by blowing up the target at the image of \(E\) at
infinity and resolving again, \(\ell(E)=e\) remains a valuation invariant, but leafhood,
\(E^2\), and the number of adjacent curves must be transformed rather than
silently reused.

On a jointly minimal pure source model, every surviving one-dicritical graph
must begin with the corridor

\[
E_{\mathrm{type\ 3},\,\ell=e}
\;--\;
A_{\mathrm{type\ 2},\,\ell=e-1}
\;--\;\cdots\;--\;
1\;--\;0\;--\;(-1)\;--\;\cdots\;--\;(-2).
\]

The remaining graph attack must combine this corridor with determinant
labels, the \(A_6\) cusp/collision data or the exact \(S_6\) jump partition,
and the degree-six pushforward/pullback equations.  Joint minimality removes
the old unbounded \(m\), while the hostile family proves that the remaining
canonical-plus-determinant data still do not give a finite enumeration.

## Primary sources

- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587--596.  Lemma 2.1 gives
  the finite linear chain and its unique meeting with the infinity boundary.
- Alexander Borisov,
  [“Frameworks for two-dimensional Keller maps”](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v27i3p54/pdf/),
  *Electronic Journal of Combinatorics* 27(3) (2020), P3.54.  The log-Keller
  formula gives the type-3 coefficient `e`; the boundary adjunction
  calculation and type-1 pullback rule give (2.4) and (4.2).
- Alexander Borisov,
  [“On two invariants of divisorial valuations at infinity”](https://people.math.binghamton.edu/borisov/documents/papers/divisorialvaluations.pdf).
  Proposition 2.2 supplies connectivity of the negative region and the
  `1--0--(-1)` transition rule; Lemma 3.7 gives (6.1), and Definition 4.1 with
  Lemma 4.1 gives the parity in (6.2).

The leaf-neighbor congruence, the joint-minimality argument, and their
application to the two surviving passports are derived in this repository.
No claim of historical priority is made.
