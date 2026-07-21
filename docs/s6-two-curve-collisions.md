# Universal collision constraints for the two-curve `S6` passports

## Claim boundary

Let

\[
F=(P,Q):\mathbb A^2_{\mathbb C}\longrightarrow\mathbb A^2_{\mathbb C}
\]

be a hypothetical Keller counterexample of generic degree six.  This note
studies the two unrestricted `S6` profiles left by the refined defect budget:

1. two distinct transposition curves; and
2. one transposition curve and one distinct 3-cycle curve.

The conclusions here do **not** assume that there is only one dicritical.
They force normalization collisions and determine their exact six-sheet fiber
profiles.  They do not eliminate either passport.  An explicit permutation
fixture at the end shows that the strongest collision pattern is compatible
with abstract `S6` monodromy.

## 1. Every nonproper component is singular

Write

\[
E_F=\{y\in\mathbb A^2:\#F^{-1}(y)<6\}.
\]

For a generically finite Keller map, this is also the nonproper-value set.
Outside the nonproper set the map is a proper unramified six-sheet cover, so
every fiber has six points.  Conversely, a nonproper value with six affine
inverse points would have six persistent inverse-function branches and an
additional branch escaping from infinity, contradicting the generic degree.

Nguyen Van Chau's Theorem 4.4 decomposes `E_F` into irreducible polynomially
parametrized curves and proves in item `(E3)` that every one of those curves
has a singularity.  Thus

\[
\boxed{\text{every irreducible component of the Keller nonproper set is
singular}.}
\tag{1.1}
\]

Every reduced branch component of the finite normalization is a nonproper
component: affine ramification is absent, so its ramification prime lies on
the normalization boundary.  The converse can fail when an index-one
dicritical contributes a nonproper curve without ramifying the finite map.

Now let a dicritical `E` dominate a branch curve `B`, and assume

\[
d_E=1,\qquad \delta_E=0.
\tag{1.2}
\]

Then `E^circ` is the normalization `B^nu=A1`.  Every local excess is
nonnegative, so their vanishing and Orevkov's Lemma 5.2 make the normalization
map an immersion at every finite normalization point.  If that map were also
injective, `B` would be a nonsingular copy of `A1`, contrary to (1.1).
Therefore

\[
\boxed{
d_E=1,\ \delta_E=0
\quad\Longrightarrow\quad
B^\nu\longrightarrow B\text{ identifies two distinct finite points}.
}
\tag{1.3}
\]

The singularity forced in (1.3) is a normalization self-collision, not an
intrinsic cusp: every individual branch is smooth and immersed.

## 2. Two transposition curves

Suppose the two ramified target components `B_2` and `B'_2` both have generic
transposition inertia.  Their dicritical realizations have

\[
(e,d)=(2,1),\qquad (e',d')=(2,1),
\]

and consume four of the five defect units.  The final unit is either an
index-one dicritical or total excess one on one of the two ramified
dicriticals.  At least one of the two ramified dicriticals consequently has
zero excess.  Applying (1.3) proves

\[
\boxed{
\text{every two-transposition passport has a zero-jump normalization
self-collision on at least one ramified curve}.}
\tag{2.1}
\]

If the residual unit is spent on an index-one dicritical, both transposition
curves have zero excess and both self-collide.  If it is one excess unit on a
ramified dicritical, the other transposition curve self-collides.

Two zero-jump normalization branches at such a collision consume `2+2`
sheets.  The remaining two sheets can be affine, can form a third
transposition branch of the same curve, or can come from the other
transposition curve.  An index-one boundary contribution can also affect the
last two sheets in the subcase where the residual unit is spent that way.
The universal conclusion is therefore (2.1), not a stronger unique fiber
census.

## 3. Transposition plus 3-cycle

Now suppose the ramified target curves are distinct curves `B_2` and `B_3`
with

\[
(e_2,d_2)=(2,1),\qquad (e_3,d_3)=(3,1).
\tag{3.1}
\]

Their refined costs already exhaust the whole budget:

\[
2\cdot1+3\cdot1=5.
\]

Hence

\[
\delta_{E_2}=\delta_{E_3}=0,
\qquad
\text{there are no index-one dicriticals}.
\tag{3.2}
\]

In particular the nonproper set and the reduced branch locus have exactly
these two components.  Equation (1.3) applies to both:

\[
\boxed{
B_2^\nu\to B_2\text{ and }B_3^\nu\to B_3
\text{ both have finite normalization self-collisions}.}
\tag{3.3}
\]

All their finite normalization branches are smooth local embeddings.  The
constant degree-six multiplicity formula now determines every relevant
fiber exactly.

### 3.1 Self-collisions of `B_3`

Each normalization branch has local degree three.  Two distinct branches
over one value consume all six sheets:

\[
3+3=6.
\]

Thus every `B_3` self-collision has exactly two branches, no affine inverse
point, and no branch of `B_2`.  It is an omitted value of the original affine
map.  Its local monodromy has two disjoint 3-cycle blocks:

\[
\boxed{
B_3\text{ self-collision}:\quad 3+3,\quad
F^{-1}(y)=\varnothing,\quad C_3\times C_3.}
\tag{3.4}
\]

### 3.2 Self-collisions of `B_2`

Two normalization branches consume four sheets.  There are exactly two
possibilities:

\[
\boxed{
B_2\text{ self-collision}:\quad
2+2+1+1\quad\text{or}\quad2+2+2.}
\tag{3.5}
\]

The first row has two affine inverse points.  In the second row a third
normalization branch consumes the last two sheets and the value is omitted.
A branch of `B_3` cannot pass through either value, since `2+2+3>6`.

### 3.3 Intersections of `B_2` and `B_3`

At a cross-intersection one smooth `B_2` branch consumes two sheets and one
smooth `B_3` branch consumes three.  Any additional boundary branch would
exceed degree six.  Exactly one affine inverse point remains:

\[
\boxed{
B_2\cap B_3:\quad 2+3+1.}
\tag{3.6}
\]

The local monodromy is `C2 x C3` on disjoint two- and three-sheet supports,
with the affine sheet fixed.

Combining (3.3)--(3.6), the saturated two-curve passport requires two
polynomially parametrized immersed rational curves, each with normalization
`A1` and each with a finite self-collision.  The 3-cycle collisions are
omitted `3+3` values; the transposition collisions have the two rows (3.5);
and every cross-intersection has the row (3.6).

## 4. Exact hostile permutation fixture

The collision census does not contradict abstract monodromy.  In `S6` take

\[
a=(123),\qquad b=(456),\qquad c=(14),\qquad d=(25).
\tag{4.1}
\]

The pair `a,b` consists of disjoint 3-cycles, as required at a `B_3`
self-collision.  The pair `c,d` consists of disjoint transpositions, as
required at a `B_2` self-collision.  Exact group closure gives

\[
|\langle a,b,c,d\rangle|=720,
\]

so the four meridians generate all of `S6`.  The Python regression test locks
the disjoint supports, orbit blocks, and group order.  Therefore conjugacy
classes, local products, transitivity, primitivity, and normal generation
cannot eliminate (3.1).

## 5. The remaining global realization problem

The sharp remaining question for this passport is whether two immersed
polynomial `A1` curves can realize (3.4)--(3.6), admit the required complement
monodromy onto `S6`, and simultaneously satisfy the Keller compactification's
canonical and determinant equations.

The exact fixture (4.1) shows that a contradiction must use data absent from
the bare permutation representation.  Candidate inputs are:

- a braid-monodromy presentation coupling all finite singular values;
- the splice diagram joining the finite collisions to the common behavior at
  infinity; or
- Borisov's augmented-canonical and determinant labels on the complete
  boundary graph.

No existence or nonexistence result for that global realization problem is
claimed here.

## Primary sources

- Nguyen Van Chau,
  [“Non-zero constant Jacobian polynomial maps of `C^2`”](https://doi.org/10.4064/ap-71-3-287-310),
  *Annales Polonici Mathematici* 71.3 (1999), 287--310.  Theorem 4.4,
  item `(E3)`, proves that every parametrized component of the deficient-fiber
  set has a singularity.
- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587--596.  Lemma 4.2 gives
  the defect budget and Lemma 5.2 gives the zero-jump local embedding.

The collision consequences, exact fiber rows, and `S6` stopping fixture are
derived in this repository.  No claim of historical priority is made.
