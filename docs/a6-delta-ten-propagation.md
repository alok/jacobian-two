# Propagation on the two remaining delta-ten codimension-two charts

## Claim boundary

This note closes the topology-propagation obligations left in the
[codimension-two checkpoint](a6-delta-ten-codim-two.md) for the two displayed
dominant incidence charts:

- `T112 + 6N` on the labeled, `P`-unramified triple-source chart; and
- `C2 + T111 + 5N` on the displayed mixed Cramer chart.

It proves that a nonempty dense clean open of each incidence has constant
affine-complement homeomorphism type.  The exact cyclic representative on
each chart therefore excludes that chart from the required connected
six-sheet `A6` quotient with single-three-cycle meridians.

The proof is carried out over the **labeled incidence space**, not over its
possibly singular coefficient-space image.  It does not cover a component
supported on a removed `P`-critical/fourth-source, residual-rank, split,
overlap, or deeper-degeneration locus.  It is conditional on the same four
one-dicritical branch hypotheses as the rest of the delta-ten audit, and it
does not prove `JC(2)`.

The algebraic input and the topology input are deliberately separated.  The
Python and Sage certificates prove the exact incidence ranks, clean rational
members, singular schemes, and presentations.  The family-wide conclusion
below uses simultaneous embedded resolution and Thom's first isotopy lemma;
that theorem step is not a CAS calculation.

## 1. Smooth irreducible labeled incidences

For `T112`, let

\[
 S_T=V(E)\cap D(g_T)\subset\mathbb A^3_{u,v,w},
 \qquad
 E=\sigma _2^2-\sigma _1\sigma _3-\sigma _2,
 \tag{1.1}
\]

where `g_T` is the complete `P`-unramified product used by the exact
certificate: the three displayed roots are nonzero and distinct,
`sigma2 != 0`, and none equals the fourth source `-sigma3/sigma2`.  The
certificate proves that `E` is absolutely irreducible, that its localized
singular ideal is the unit ideal, and that the three affine-linear target and
tangent equations have coefficient rank three everywhere on `S_T`.
Consequently

\[
 I_T\longrightarrow S_T                                            \tag{1.2}
\]

is a torsor under the rank-one kernel bundle of a surjection
`O^4 -> O^3`.  Thus `I_T` is a smooth geometrically irreducible threefold.
Nothing here requires the coefficient-space image to be smooth, normal, or
even injectively parametrized.

For the mixed profile, let

\[
 S_M=D(g_M)\subset\mathbb A^3_{q,r,w},                              \tag{1.3}
\]

where `g_M` is the complete product of the triple, contact, split, overlap,
discriminant, fourth-root, and residual-rank factors displayed by the Cramer
certificate.  On `S_M` its `4 x 4` coefficient matrix is invertible, so
Cramer's rule identifies the mixed incidence `I_M` with `S_M`.  Therefore
`I_M` is also a smooth geometrically irreducible threefold.  This conclusion
does not need the removed residual polynomial to be absolutely irreducible:
the complement of any hypersurface in affine three-space is an irreducible
open.

## 2. The exact clean opens are nonempty and connected

Inside either incidence, impose the following nonvanishing conditions:

1. the distinguished singularity has exactly the prescribed contact orders;
2. every residual collision is a transverse node;
3. residual node targets are pairwise distinct and avoid the distinguished
   targets and the forced cusp target;
4. there is no additional diagonal critical point or fourth-source target
   collision;
5. the `T(2,5)` cusp coefficient remains nonzero; and
6. the branch at infinity retains the fixed `T(5,9)` type.

These conditions define Zariski opens.  For example, away from the diagonal
an ordered residual collision is cut out by

\[
 P(t)-P(s)=0,\qquad Q(t)-Q(s)=0.                                   \tag{2.1}
\]

At such a pair, the Jacobian determinant in the two source variables is, up
to sign,

\[
 \det\!\begin{pmatrix}
 P'(t)&P'(s)\\
 Q'(t)&Q'(s)
 \end{pmatrix}.                                                     \tag{2.2}
\]

Its nonvanishing says exactly that the two normalization branches are
transverse.  After removing the diagonal, distinguished sections, and
projective boundary, the residual double-point scheme is finite etale of
constant degree.  Higher contact, coincident target sections, and extra
critical points are closed conditions detected by the exact discriminants,
resultants, and jet determinants already evaluated by the certificates.

Call the resulting opens `I_T^cl` and `I_M^cl`.  They are nonempty: the exact
`T112` representative has six transverse residual nodes in two checked
packets, and the mixed representative has five transverse residual nodes;
all recorded separation and localizer values are nonzero.  The fact that the
`T112` sample has `k=2` is harmless.  At that coefficient the generic pair
formula splits into its true vertical and graph packets, whose exact
discriminants and mutual resultants are nonzero.  The divisor `k^2-4=0` is
not part of the `P`-unramified root-base localization, so the sample is not
trapped on an equisingularity boundary.

A nonempty open of an irreducible variety is irreducible.  Hence both clean
incidences are irreducible, and their complex analytifications are connected
and path-connected.  This is the missing connectedness argument; a boolean
field in a generated certificate would not by itself prove it.

## 3. Finite etale labeling

Over either clean incidence the residual node-pair scheme is finite etale,
of degree six for `T112` and degree five for the mixed profile.  The source
schemes of the distinguished contact and triple branches are finite etale as
well, because their discriminants are units on the clean open.

Take an ordering cover, or a finite etale Galois closure,

\[
 B'\longrightarrow I^{\mathrm{cl}},                                \tag{3.1}
\]

on which all residual nodes, their two source preimages, and every
distinguished branch are labeled by global sections.  Choose a connected
component containing a lift of the cyclic sample.  A connected component of
a finite etale cover has open-and-closed image; since the base is connected,
that image is the whole base.  Thus this component is still surjective and
is analytically path-connected.

The cover is only a bookkeeping device for resolution centers.  It does not
assert that the coefficient-image map has degree two.  In particular, the
unique-`T112`-fiber/two-orientation argument is unnecessary for topology
transport.

## 4. Relative embedded resolution

Over `B'`, all singular points and required infinitely-near centers are
smooth disjoint sections.  The resolution sequence is uniform:

- blow up each node section once;
- blow up an ordinary triple section once;
- for a contact-two section, blow up the target and then the common
  infinitely-near tangent point;
- for `T112`, the first blowup separates the transverse branch and the second
  resolves the remaining tangent pair; and
- use the fixed embedded resolution sequences for the cusp and the union of
  the infinity branch with `L_infinity`.

The two fixed branches really are uniform.  At `t=0`,

\[
 P(t)=t^2\cdot\text{unit},\qquad Q(t)=t^5\cdot\text{unit},           \tag{4.1}
\]

so the cusp type is `T(2,5)`.  With `z=1/t` in the `Q`-chart at infinity,

\[
 \left(\frac P Q,\frac1Q\right)
   =\left(z^5\cdot\text{unit},z^9\cdot\text{unit}\right),          \tag{4.2}
\]

so the infinity pair is fixed at `(5,9)`.  The infinitely-near centers are
algebraic sections because they are determined by labeled tangent and higher
jet data whose leading coefficients are units on the clean open.

Successive blowups along smooth relative sections give a smooth proper
family

\[
 \widetilde{\mathcal P}\longrightarrow B'                           \tag{4.3}
\]

whose reduced total transform of
`overline(B) union L_infinity` is a relative simple-normal-crossings divisor
with constant labeled intersection graph.  This is an explicit embedded
equiresolution in the sense discussed by Lipman in
[Equisingularity and simultaneous resolution of
singularities](https://arxiv.org/abs/math/9802010), Section 4.

## 5. Proper isotopy and the affine complements

Stratify the resolved total space by the complement of the total-transform
divisor, the smooth part of each labeled divisor component, and every
intersection stratum.  A relative SNC divisor has a Whitney stratification
of this form, and every stratum maps submersively to `B'`.  The ambient map is
proper because it is obtained from `P2 x B'` by projective blowups.

Thom's first isotopy lemma now gives local topological triviality.  In the
form used here, it says that a proper map on a closed Whitney-stratified set
which is submersive on every stratum is locally trivial; see Mather,
[Notes on Topological Stability](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/matherj.pdf),
Proposition 11.1.

The blowup map is an isomorphism off the total transform.  Therefore the
locally trivial resolved complement is exactly the original affine
complement, not merely a space with the same singularity counts.  Since `B'`
is path-connected, all fibers satisfy

\[
 \mathbb C^2\setminus B_{b_0}\ \cong\
 \mathbb C^2\setminus B_{b_1}.                                     \tag{5.1}
\]

The labeling cover is surjective, so the same conclusion holds on the
original clean incidence open.

The exact Sage representatives on both charts have affine complement group
`Z`.  Hence every member of these dense clean incidence opens has complement
group `Z`, and neither can realize the required `A6` monodromy quotient.

## 6. Exact consequence and remaining boundary

The strongest conclusion is:

> The displayed `P`-unramified `T112` incidence component and the displayed
> mixed `C2+T111` Cramer component each contain a nonempty dense clean open
> with constant affine-complement homeomorphism type.  Their exact cyclic
> representatives exclude both dominant charts from the conditional
> single-three-cycle `A6` passport.

This closes the two propagation qualifications in the previous checkpoint.
It does **not** promote any of the following to solved:

- `P`-critical/fourth-source boundaries;
- compatible residual-rank loci;
- split and overlap charts not already represented by the clean incidence;
- intersections of different incidence components;
- deeper collision profiles; or
- any component not dominated by the two displayed incidences.

Constant singular-scheme length alone would not justify the result.  The
finite-etale residual scheme, labeled relative blowups, relative SNC divisor,
properness, and stratumwise submersion are the hypotheses that make the
topological transport valid.
