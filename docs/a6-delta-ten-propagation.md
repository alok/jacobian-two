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

## 7. Extension over the split-containing total bases

The previous argument used solved source charts.  Passing from those charts
to the total unordered-pair bases is not automatic: an equation-level
restriction at one witness would prove neither universal containment nor a
smooth base for isotopy.  The split closure certificates now supply the
missing checks.

### 7.1 Exact row bridges and regular source bases

For one unordered source pair retain both its sum `s` and product `p` and use

\[
 F(k,s,p)=(2s+k)p-s(s^2+ks+1).                                  \tag{7.1}
\]

This surface is smooth and geometrically integral.  If `C` is the
coefficient of `t` in `Q mod (t^2-st+p)` and

\[
 D=F_p\partial_s-F_s\partial_p,                                  \tag{7.2}
\]

then, on the principal chart

\[
 p=\frac{s(s^2+ks+1)}{2s+k},
\]

the old collision polynomial `H` satisfies

\[
 C=\frac{s^2}{(2s+k)^4}H,
 \qquad D=(2s+k)\frac d{ds}.                                     \tag{7.3}
\]

Consequently `(C,D(C),D^2(C))` is an exact lower-triangular transform
of `(H,H',H'')` with determinant

\[
 \frac{s^6}{(2s+k)^9}.                                           \tag{7.4}
\]

For two ordered contact pairs, the corresponding block determinant is

\[
 \frac{u^4v^4}{(2u+k)^7(2v+k)^7}.                                \tag{7.5}
\]

The certificates verify the full affine rows, including constant terms, and
give explicit inverses in the relevant localizations.  They also prove by
power-membership identities that all eight component-row determinants for
the non-overlap split `C3` and `C2^2` allocations are units on their **full**
declared clean opens, rather than merely nonzero at one witness.

The ordered two-pair base is integral by geometric integrality of its generic
pair curve, product integrality over `Qbar(k)`, flatness over `QQ[k]`, and
torsion descent.  Its entire singular support consists of the three points

\[
\begin{array}{c|c|c}
k&(s_1,p_1)&(s_2,p_2)\\ \hline
0&(0,1/2)&(0,1/2)\\
2&(-1,0)&(-1,0)\\
-2&(1,0)&(1,0),
\end{array}                                                       \tag{7.6}
\]

where both pairs are the same vertical/graph overlap.  Every clean
maximal-rank two-contact allocation removes these points.

For the mixed profile, retain the separate contact product as well.  With
`n=r^2+q^2-q^3`, the total source base is

\[
 (2qrs+n)p-s\bigl(qr(s^2+1)+ns\bigr)=0.                          \tag{7.7}
\]

It is primitive linear in `p`, hence geometrically integral, and its
singular ideal becomes the unit ideal after localizing by `qr`.  On the old
ordinary chart, its four target rows are exactly the historical mixed rows
under the diagonal transform

\[
 \operatorname{diag}\!\left(
   \frac{(2s+k)^4}{s^2},
   -\frac{(2s+k)^5}{s^2(s^2+ks+1)},
   1,1
 \right),                                                        \tag{7.8}
\]

whose determinant is a unit there.  Independent Rabinowitsch saturations
prove that the four `T112` and three mixed component transforms remain units
on every stated split clean localization.  Thus these are restrictions of
the same global incidences carrying the cyclic samples, including the mixed
vertical pair-denominator chart.

The source-base primality, smoothness, singular-support, and determinant-unit
claims are replayed independently by
[`check_a6_delta_ten_split_component_closure.sage`](../tools/check_a6_delta_ten_split_component_closure.sage).

### 7.2 Proper isotopy on the enlarged clean opens

For each total incidence, now remove the closed loci where the coefficient
rank drops, two distinguished pairs coincide, a source pair ceases to be
etale or immersed, the exact contact-order jet vanishes, a residual collision
ceases to be a transverse node, or two distinguished targets meet.  The
exact split witnesses show that every required rank-open allocation has a
point in this clean locus.  Sections 1--2 and the regularity checks above show
that each resulting labeled incidence is a nonempty open of a smooth
irreducible total incidence, hence its complex analytification is
path-connected.

After the finite-etale ordering cover of Section 3, every branch, residual
node pair, and infinitely-near center is a global smooth section.  The local
resolution sequences are constant: one blowup per node, two for each
contact-two point, three for a contact-three point, the `T112` sequence from
Section 4, and the fixed cusp and infinity sequences.  The highest required
jet difference is a unit on the clean open, so the last exceptional
directions remain distinct.  The centers for different target sections are
disjoint.  These relative blowups therefore give a proper smooth family with
a relative SNC total transform and constant labeled intersection graph.

The Whitney stratification by divisor components and intersections maps
submersively to the clean base.  Thom's first isotopy lemma applies exactly
as in Section 5.  Therefore the cyclic complements on the old ordinary
charts propagate over every clean maximal-rank split `C3`, `C2^2`, `T112`,
and mixed locus.

This conclusion deliberately does **not** manufacture topology from the row
identities.  The row bridges prove algebraic containment; the preceding
smoothness, finite-etale labeling, unit-jet, relative-SNC, properness, and
submersion argument is the separate topology step.  Exact arcs also place
all three prescribed overlap allocation incidences in the algebraic closures
of the ordinary components, but those overlap points are outside the smooth
rank-open argument above.  Their complement topology, and the finite
rank-three affine-line fibers, remain open.

## 8. Reuse for the nonsplit codimension-three contact/fiber surfaces

The same argument applies to the first codimension-three profile on the
determinant-nonzero `C4+6N` Cramer surface.  Its incidence is a principal open
of the `(k,s)` plane and hence smooth and irreducible.  The exact cyclic sample
lies in the open where the contact branches are immersed, the fourth jet is
nonzero, the six residual nodes form a finite-etale scheme with separated
targets, and the cusp and infinity types remain fixed.

There is one resolution detail that must not be hidden in the phrase
“uniform contact resolution.”  A two-branch contact of order four requires
four relative blowups.  After three, the strict transforms are transverse to
one another but they and the newest exceptional component meet at a triple
crossing; a fourth blowup separates the three directions and produces a
relative SNC divisor.  Once the branches and nodes are labeled by a finite
etale cover, all four centers are algebraic sections with nonzero defining
jets.  Proper Thom isotopy therefore transports the cyclic complement over
the clean Cramer open exactly as in Sections 3--5.

The `C2+C3+5N` compatibility surface uses the same mechanism on its smooth
rank-four clean open.  Sage derives an irreducible compatibility
hypersurface over `QQ`; a smooth rational point proves geometric
irreducibility and shows that the clean open is nonempty.  Exact saturation
also proves that the hypersurface is smooth on the full declared valid chart.
After a
finite-etale cover, the two branches at each contact target and the five
nodes are disjoint labeled sections.  Two relative blowups resolve the
contact-two section and three resolve the contact-three section.  Along with
the fixed node, cusp, and infinity resolutions, these centers give a relative
SNC divisor in a smooth proper family.  Proper Thom isotopy thus transports
the cyclic sample throughout this dense clean surface.  Exact saturation
separately proves that the coefficient-rank-drop curve is wholly
inconsistent.  Removed pair-chart and non-clean equisingular boundaries lie
outside the topology argument and remain open.

The `C2+Q0+2N` component is a principal open of the product of a smooth
genus-one compatibility curve with the quadruple-fiber parameter.  Its exact
rational member proves that the clean open with four distinct quadruple
slopes, an immersed exact contact pair, and two separated nodes is nonempty.
On a finite-etale ordering cover, one relative blowup resolves the ordinary
quadruple section, two resolve the contact-two section, and one resolves each
node.  The same fixed cusp/infinity resolutions and proper Thom argument
transport the cyclic complement throughout the connected clean surface.
The split, singular-fiber, same-target, and non-clean boundaries are outside
this argument.

The `C3+T111+4N` compatibility surface is likewise geometrically
irreducible and smooth on its full declared valid chart.  A smooth rational
point and two independent coefficient-image tangents make its clean
rank-four open a genuine nonempty surface.  After finite-etale labeling, one
relative blowup resolves the ordinary triple section, three resolve the
contact-three section, and one resolves each of the four node sections.  The
fixed cusp and infinity sequences then complete a relative SNC divisor in a
smooth proper family.  Proper Thom isotopy transports the cyclic complement
throughout this connected clean open.  Split fibers, singular four-point
fibers, pair/triple overlap and same-target factors, non-clean slope or node
collisions, and deeper intersections remain outside the transport theorem.

The full calculations and their omitted split, denominator, cusp-pair,
diagonal, singular-fiber, same-target, non-clean, and residual boundaries are recorded in
[the codimension-three checkpoint](a6-delta-ten-codim-three.md).
