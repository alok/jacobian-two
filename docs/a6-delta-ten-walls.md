# The two dominant conditional `A6` delta-ten walls are cyclic

## Claim boundary

This note continues the [generic delta-ten audit](a6-delta-ten-generic.md)
under the same four additional hypotheses on the one-dicritical `A6` branch:

1. its normalization is `A1` and is represented by a polynomial map;
2. its projective closure has exactly one genuine Puiseux pair at infinity;
3. its only intrinsic finite singularity is the forced `T(2,5)` cusp, whose
   image has no additional normalization preimage; and
4. every other finite singularity is a collision of smooth normalization
   points.

Every collision-delta-ten normalization is polynomially equivalent over `C`
to

\[
  P=t^2+kt^3+t^4,
  \qquad
  Q=at^5+bt^6+ct^7+dt^8+t^9,
  \qquad a\ne0. \tag{0.1}
\]

The clean ten-node locus was already excluded by a cyclic-complement
representative.  This note treats the two dominant degeneration divisors:

- one contact-two collision plus eight nodes; and
- one ordinary triple target plus seven nodes.

Exact localized Gröbner calculations prove that each valid incidence has one
irreducible four-dimensional dominant component in the five-dimensional
parameter space.  An exact rational member on each component has affine
complement `Z`; each raw four-meridian presentation has an exhaustive
`40^4` single-three-cycle replay with only forty diagonal `C3` images and no
`A6` image.  Proper projective Whitney--Thom propagation therefore excludes
the connected generic equisingular open of both divisors.

The follow-up [codimension-two checkpoint](a6-delta-ten-codim-two.md) now
excludes four identified generic or dominant components and supplies cyclic
members plus exact incidence algebra for the `T112` and mixed
contact-plus-triple candidates.  Those last two generic exclusions retain
explicit labeling or connected-equisingular propagation obligations.  This
still does not complete the delta-ten wall audit: compatible
lower-dimensional residual-rank loci, the split `k=0,+2,-2` charts,
removed pair or `P`-projection/critical-fiber loci, deeper intersections, and
endpoints remain to be eliminated or checked.  A removed `P`-critical fiber
is not automatically a singular curve branch when `Q'` is nonzero.  The
result is conditional and computer-assisted.  It does not derive the four
hypotheses, exclude unrestricted `A6` or `S6`, construct a Keller map, or
prove the plane Jacobian conjecture.

## 1. Pair collisions and target fibers are different data

For distinct normalization parameters `t,u`, put

\[
  s=t+u,\qquad r=tu,
\]

and let `H(s)` be the collision decic from the generic audit.  On the chart

\[
  p:=k+2s\ne0,
\]

the first divided difference determines

\[
  r=\frac{s(s^2+ks+1)}{k+2s}. \tag{1.1}
\]

Write

\[
  q=s^2+ks+1,
  \qquad
  \Delta=2s^2+3ks+4. \tag{1.2}
\]

Then

\[
  (t-u)^2=-\frac{s\Delta}{p}. \tag{1.3}

Thus `s*Delta=0` is diagonal, `p=0` is the missing pair-incidence chart, and
`q=0` places one member of the pair over the forced cusp image.  All three
must be separated from a genuine smooth-branch collision.

On this valid nonsplit chart, a root of `H` of multiplicity `m` represents
one unordered pair with intersection multiplicity `m`.  It does **not**
record how several different pairs are grouped over one target point.  In
particular, three source points over an ordinary triple point give three
distinct pair sums and therefore three simple roots of `H`.  A complete wall
audit must track both:

1. the length partition of the finite double-point scheme; and
2. the map from that scheme to the full target `(X,Y)`.

This distinction is also essential on `k=0,+2,-2`.  There the specialized
generic decic has forced cancellation powers, while the genuine pair scheme
splits into the `V/W` components derived in the generic note.  The equations
`k=0,+2,-2` themselves are not degeneration walls.

## 2. The dominant contact incidence has rank two everywhere

Let `T(s)` be the exact tangency polynomial.  The syzygy

\[
  qT=s\Delta p\,H'
    +(9k^2s-2ks^2+10k-4s^3-16s)H \tag{2.1}

shows that on

\[
  s\,p\,q\,\Delta\ne0, \tag{2.2}

the equations `H=T=0` are equivalent to `H=H'=0`.  They therefore describe
a genuine tangential collision of two distinct immersed branches.

Both `H` and `T` are affine-linear in `(a,b,c,d)`.  Form their coefficient
matrix

\[
  M=
  \begin{pmatrix}
    H_a&H_b&H_c&H_d\\
    T_a&T_b&T_c&T_d
  \end{pmatrix}. \tag{2.3}

For its six `2 x 2` minors `D_ij`, exact polynomial gcd gives

\[
  \gcd_{i<j}D_{ij}
  =s(k+2s)^3(2s^2+3ks+4)
  =sp^3\Delta. \tag{2.4}

After division by this common factor, a lexicographic Gröbner basis is

\[
\begin{aligned}
 &p^4,\\
 &p^3(s^2-1),\\
 &p^2(s^2-1)^2,\\
 &(k+2s^3)(s^2-1)^3,\\
 &s^2(s^2-1)^4.
\end{aligned} \tag{2.5}

Its residual zero set lies on `p=0` at `s=0,+1,-1`.  The executable
certificate also checks the stronger localized unit-ideal statement

\[
  \left(
    \frac{D_{ij}}{sp^3\Delta},
    1-zsp\Delta
  \right)=(1). \tag{2.6}

Consequently `M` has rank exactly two throughout the valid base

\[
  B_{\mathrm c}=D(spq\Delta)\subset\mathbb A^2_{k,s}. \tag{2.7}

The incidence

\[
  I_{\mathrm c}=V(H,T)|_{B_{\mathrm c}}
  \subset B_{\mathrm c}\times\mathbb A^4_{a,b,c,d} \tag{2.8}

is an affine rank-two bundle over an irreducible base.  Hence it is
irreducible of dimension four.

Because `H` is monic of degree ten in `s`, the projection of its closure to
parameter space is finite.  The image

\[
  W_{\mathrm c}\subset\mathbb A^5_{k,a,b,c,d} \tag{2.9}

is therefore a closed irreducible hypersurface.  It is the unique dominant
contact wall on the nonsplit valid chart.

### 2.1 Why the omitted rank boundaries are not extra valid divisors

The exact substitutions identify every common rank factor.

At `s=0`,

\[
  H=ak^2,\qquad T=10ak^3. \tag{2.10}

The branch `a=0` destroys the forced `T(2,5)` cusp, while `k=0` is the split
pair-incidence fiber.

At `p=0`,

\[
  H=s^2(s^2-1)^4,
  \qquad
  T=-36s^3(s^2-1)^4. \tag{2.11}

The only points are `(k,s)=(0,0),(-2,1),(2,-1)`, exactly the three split
chart artifacts.

At `Delta=0`, substituting

\[
  k=-\frac{2s^2+4}{3s}

\]

reduces both equations to the common linear form

\[
  L_s=80a+48bs+28cs^2+16ds^3+9s^4. \tag{2.12}

With `t=s/2`, this is precisely the condition `P'(t)=Q'(t)=0`; its image is
the invalid extra-critical locus `L=0`.

Finally, on `q=0`, exact reduction gives the value and tangent equations for
a smooth normalization point over the forced cusp image.  This is removed by
the standing cusp-image hypothesis.  Genuine contacts internal to a split
`V/W` fiber require an additional discriminant or tangent equation and hence
belong to the deeper audit; no whole hyperplane `k=0,+2,-2` is removed.

## 3. An exact cyclic contact-two member

Take

\[
  (k,a,b,c,d)
  =\left(1,-\frac{12}{5},-\frac{16}{5},0,0\right). \tag{3.1}

Then

\[
  P=t^2+t^3+t^4,
  \qquad
  Q=t^9-\frac{16}{5}t^6-\frac{12}{5}t^5. \tag{3.2}

The collision polynomial factors as

\[
  H(s)=\frac{(s+2)^2}{5}K_8(s), \tag{3.3}

where

\[
\begin{aligned}
K_8={}&5s^8+10s^7+25s^6+29s^5+54s^4\\
&+27s^3-21s^2-18s-3.
\end{aligned}

The exact checks include

\[
  K_8(-2)=1269,
  \qquad
  \operatorname{Disc}(K_8)=-52320062156396485500, \tag{3.4}

and

\[
  C=\frac{133}{25},
  \qquad
  L=\frac{852228}{25}. \tag{3.5}

The unique double collision has `(s,r)=(-2,2)`, normalization parameters
`-1+i,-1-i`, and image `(-2,-128)` after replacing `Q` by `5Q`.  The two
branches have common slope `96`, while their curvature difference is

\[
  \frac{2256}{169}(t+1)\ne0. \tag{3.6}

Thus the contact order is exactly two.  The residual octic gives eight
ordinary nodes with distinct target images.  The complete singularity and
genus ledger is

\[
  28=2\quad(T(2,5))
     +2\quad(A_3\text{ contact})
     +8\quad(\text{nodes})
     +16\quad(T(5,9)\text{ at infinity}). \tag{3.7}

At the incidence point, one coefficient minor is `131220`,

\[
  \gcd(H,T)=s+2,
  \qquad
  H''(-2)=\frac{2538}{5},
  \qquad
  T'(-2)=\frac{30456}{5}. \tag{3.8}

Therefore the sample is a regular point of the dominant component, and the
generic member of `W_c` has one unique contact-two collision.

Sage independently regenerates the primitive implicit equation, the affine
Jacobian primary decomposition, and the raw van Kamp presentation.  The
Jacobian components have `(length, radical degree)`

\[
  (4,1),\qquad(3,1),\qquad(8,8), \tag{3.9}

for the cusp, contact, and nodes.  The four-generator, twelve-relator group
simplifies to `Z`, with all four geometric meridians equal to its generator.
The dependency-free replay nevertheless checks the raw presentation directly:

\[
  40^4=2{,}560{,}000

\]

formal assignments reduce to exactly forty satisfying assignments, all
diagonal with image `C3`; none generates `A6`.

## 4. The ordinary-triple incidence is irreducible

Let `t,u,v` be three distinct nonzero parameters over one target and put

\[
  \sigma_1=t+u+v,
  \qquad
  \sigma_2=tu+tv+uv,
  \qquad
  \sigma_3=tuv. \tag{4.1}

If `w` is the fourth root of `P(z)-P(t)`, Vieta's formulas give

\[
  w=-\frac{\sigma_3}{\sigma_2},
  \qquad
  k=\frac{\sigma_3}{\sigma_2}-\sigma_1, \tag{4.2}

and the base equation

\[
  \boxed{
  \sigma_2^2-\sigma_1\sigma_3-\sigma_2=0.} \tag{4.3}

As a monic quadratic in `sigma_2`, its discriminant is

\[
  1+4\sigma_1\sigma_3. \tag{4.4}

This polynomial is squarefree and has degree one in `sigma_1`; it is not a
square in `C(sigma_1,sigma_3)`.  Gauss's lemma therefore proves that (4.3)
defines a geometrically irreducible surface.

The two equations `Q(t)=Q(u)` and `Q(t)=Q(v)` have coefficient rows

\[
  (t^n-u^n)_{n=5}^8,
  \qquad
  (t^n-v^n)_{n=5}^8 \tag{4.5}

in `(a,b,c,d)`.  Put

\[
  \Omega=tuv(t-u)(t-v)(u-v)(tu+tv+uv). \tag{4.6}

For the ideal `J_4` of the six coefficient minors and the ideal `J_5` of the
ten minors after adjoining the fixed `n=9` column, exact Rabinowitsch
localization gives

\[
  (\text{(4.3)},J_4,1-z\Omega)=(1),
  \qquad
  (\text{(4.3)},J_5,1-z\Omega)=(1). \tag{4.7}

Hence the coefficient matrix has rank two everywhere on the valid distinct
triple base; the inhomogeneous equality system is always consistent and cuts
an affine two-plane in `(a,b,c,d)`.  The ordered triple incidence is an
affine rank-two bundle over an open of an irreducible surface, so it is
irreducible of dimension four.  Its parameter-space image has an irreducible
four-dimensional closure

\[
  W_{\mathrm t}\subset\mathbb A^5. \tag{4.8}

On the standing valid locus, the polynomial map is generically injective and
its finite double-point scheme has length ten, so it has only finitely many
unordered triple target fibers.  The sample below has exactly one.  Thus the
incidence is generically finite over its image, and `W_t` is the unique
dominant ordinary-triple hypersurface.

## 5. An exact cyclic ordinary-triple member

Take

\[
  (k,a,b,c,d)
  =\left(2,\frac{114}{625},-\frac15,-3,-3\right). \tag{5.1}

Then

\[
  P=t^2(t+1)^2,
  \qquad
  Q=\frac{114}{625}t^5-\frac15t^6-3t^7-3t^8+t^9. \tag{5.2}

The three parameters

\[
  -\frac35,\qquad-\frac25,\qquad\frac15

\]

map to the ordinary triple point

\[
  (X,Y)=\left(\frac{36}{625},0\right). \tag{5.3}

The fourth point in this `P`-fiber is `-6/5`, but its second coordinate is

\[
  Q\left(-\frac65\right)=-\frac{653184}{78125}\ne0, \tag{5.4}

so this is not a quadruple point.

The value `k=2` is treated with the honest split pair incidence.  Its two
collision factors are

\[
  V(r)=\frac{(25r-6)(25r^3-469r^2+582r-144)}{625}, \tag{5.5}

and

\[
  W(s)=\frac{(5s+1)(5s+2)
  (25s^4+185s^3-413s^2-1497s+228)}{625}. \tag{5.6}

Both are reduced.  The three linear factors encode the three pairs at the
triple point; the residual cubic and quartic encode seven ordinary nodes.
Their discriminants and mutual target-separation resultants are nonzero.

At the triple point, with

\[
  x=X-\frac{36}{625},\qquad y=Y,

\]

the tangent cone is

\[
  \frac{93312}{15625}
  (18x+21875y)
  (168x-3125y)
  (1782x-3125y). \tag{5.7}

The three distinct factors prove that the branches are transverse.  The
validity factors are

\[
  C=\frac{746496}{390625},
  \qquad
  L=-\frac{9906624}{15625}. \tag{5.8}

The genus ledger is

\[
  28=2\quad(T(2,5))
     +3\quad(\text{ordinary triple})
     +7\quad(\text{nodes})
     +16\quad(T(5,9)\text{ at infinity}). \tag{5.9}

For the ordered triple, the value of `Omega` in (4.6) is
`-72/390625`, and one coefficient minor is `8652/48828125`.  The sample lies
on the dominant rank-two component even though its convenient representative
uses the legitimate `k=2` chart.

The independent Sage checker regenerates the primitive irreducible degree-nine
implicit curve, its singular scheme, and the exact raw presentation.  The
Jacobian algebra has total length fifteen and radical length nine, with

- a length-four cusp component;
- a length-four ordinary-triple component;
- one reduced degree-four node component; and
- one reduced degree-three node component.

Again, the affine complement simplifies to `Z`, all four raw meridians map to
one generator, and the exhaustive raw-presentation replay has exactly forty
diagonal `C3` images and zero `A6` images.

## 6. Propagation over the two dominant walls

Let `W_c^clean` be the open subset of the irreducible contact hypersurface on
which the curve has exactly one contact-two point, eight nodes, the forced
cusp, and the fixed infinity branch.  Let `W_t^clean` be the analogous open
of the irreducible ordinary-triple hypersurface.  Sections 3 and 5 prove that
both opens are nonempty.

Each nonempty Zariski-open subset of an irreducible complex variety is
irreducible and analytically connected.  After a finite base change labels
the collision sections, apply simultaneous resolution and proper projective
Whitney--Thom isotopy to

\[
  (\mathbb P^2,\overline B_\lambda\cup L_\infty). \tag{6.1}

The affine complements are constant over each connected equisingular open.
Therefore the exact cyclic representatives exclude the required
single-three-cycle `A6` passport throughout both dominant wall strata.

This propagation is a theorem dependency.  It is not inferred from
connectedness plus a constant singularity count alone.  A reference for the
equisingularity and simultaneous-resolution input is Joseph Lipman's
[“Equisingularity and simultaneous resolution of singularities”](https://arxiv.org/abs/math/9802010).

Combining this note with the generic audit excludes all three open
five-/four-dimensional layers now known at conditional delta ten:

| stratum | expected dimension | exact representative complement | `A6` images |
|---|---:|---|---:|
| ten nodes | 5 | `Z` | 0 |
| contact two + eight nodes | 4 | `Z` | 0 |
| ordinary triple + seven nodes | 4 | `Z` | 0 |

## 7. The complete finite candidate ledger

The pair-length partition alone has the 42 ordinary partitions of ten, but
target-fiber grouping produces a finer finite ledger.  A target fiber has at
most four smooth normalization branches because `P(t)-X` is quartic.  For
pairwise branch contacts `q_ij`,

\[
  \delta=\sum_{i<j}q_{ij}. \tag{7.1}

The contact-tree enumeration through delta ten has:

- 9 two-branch contacts `C2,...,C10`;
- 15 three-branch profiles `(q,q,k)`; and
- 11 four-branch contact trees.

Thus there are 35 non-node local atoms.  Taking multisets of these atoms and
filling the unused collision delta with nodes gives exactly 145 global
candidate profiles.  Their multijet *expected*-codimension distribution is

| expected codimension | profile count |
|---:|---:|
| 0 | 1 |
| 1 | 2 |
| 2 | 6 |
| 3 | 14 |
| 4 | 29 |
| 5 | 38 |
| 6 | 32 |
| 7 | 16 |
| 8 | 6 |
| 9 | 1 |

The two codimension-one profiles are exactly the walls excluded above.  The
next six expected codimension-two profiles are:

1. contact three plus seven nodes;
2. two contact-two points plus six nodes;
3. a triple with contacts `(1,1,2)` plus six nodes;
4. a separate contact two and ordinary triple plus five nodes;
5. an ordinary quadruple plus four nodes; and
6. two ordinary triple points plus four nodes.

The ledger also contains 38 expected zero-dimensional profiles and 55
overdetermined profiles of expected codimension greater than five.  These 55
are retained as candidates.  Special coefficient relations can defeat a
naive dimension count, so exact saturation or containment is required before
any is declared empty.

The six codimension-two profiles have now been audited on their displayed
generic or Cramer opens.  Exact incidence calculations and cyclic
representatives exclude the identified components for contact three,
ordinary quadruple incidence, two contacts over the full ordered Cramer base,
and the dense double-triple Cramer open.  For `T112`, the length-two Sage
fiber proves generic finiteness while degree two additionally uses the unique
`T112`/two-orientation labeling argument.  For a separate contact and
ordinary triple, the calculation establishes a dense rational Cramer
component and cyclic sample; connected equisingular propagation remains a
separate obligation.  No statement eliminates compatible lower-dimensional
pieces supported on residual rank factors, split charts, removed pair or
`P`-projection/critical-fiber loci, or deeper intersections.  In particular,
the removed triple/fourth-root loci are not automatically singular branches
of `(P,Q)`.  See the
[codimension-two checkpoint](a6-delta-ten-codim-two.md) for the six distinct
claim boundaries and reproduction commands.

One useful coefficient-slice check makes the next elimination finite.  For
fixed `k`, the collision-decic coefficients are affine-linear in
`(a,b,c,d)`.  Two `4 x 4` coefficient minors are

\[
  256k,
  \qquad
  128(3k^2+4). \tag{7.2}

They never vanish together over `C`, so the coefficient map has rank four
for every `k`.  Also the `s^9` coefficient is `6k`, recovering `k`.  This
supports a delta-seven-style root-partition elimination, but target-fiber
incidence and the split `V/W` charts remain separate required calculations.

## 8. Exact layers, topology dependencies, and what remains

The following are exact and dependency-free once the stored relation words
are present:

- the contact and triple incidence equations;
- all gcd, Gröbner, Rabinowitsch-localization, factor, discriminant, and
  resultant identities;
- the exact geometry and genus accounting of both rational representatives;
- the 35-atom and 145-profile combinatorial ledger;
- the coefficient-slice rank calculation; and
- both exhaustive `40^4` finite-group replays.

The checked Sage sources independently regenerate:

- both primitive implicit equations;
- both affine singular-scheme primary decompositions;
- both four-generator raw van Kamp presentations; and
- both cyclic simplification maps and sections.

Sage's Zariski--van Kamp implementation remains a computer-assisted
dependency.  Family-wide propagation separately depends on proper projective
Whitney--Thom triviality.

The next exact work package is the unresolved boundary of the six
codimension-two profiles in Section 7, followed by the fourteen expected
codimension-three profiles.  A complete audit must:

1. saturate and decompose the residual coefficient-rank loci left by the
   two-contact, contact-plus-triple, and double-triple Cramer calculations;
2. audit all six profiles on the true split `k=0,+2,-2` pair charts;
3. determine the removed pair, overlap, and
   `P`-projection/critical-fiber loci on appropriate charts rather than
   declaring them invalid; prove whether each is contained in a known
   component or gives a genuine deeper component;
4. finish the `T112` geometric labeling/propagation argument and separately
   verify the connected equisingular clean open for the mixed
   contact-plus-triple component;
5. regenerate one presentation on every resulting connected equisingular
   component;
6. transport finite endpoint obstructions across arithmetic conjugacy only
   with the finite-etale/Riemann-existence hypotheses stated explicitly; and
7. prove emptiness, invalidity, or containment for every overdetermined
   profile instead of relying on expected dimension.

Even completing that list would still prove only the stated conditional
one-pair, single-three-cycle `A6` theorem.  The unrestricted degree-six
frontier and `JC(2)` remain open.

## 9. Reproduction

Run the dependency-free certificates with:

```bash
uv run python -m scripts.a6_delta_ten_contact_wall
uv run python -m scripts.a6_delta_ten_triple_wall
uv run python -m scripts.a6_delta_ten_wall_components
uv run python -m scripts.a6_delta_ten_partition_ledger
uv run python -m scripts.a6_delta_ten_contact_three
uv run python -m scripts.a6_delta_ten_double_contact
uv run python -m scripts.a6_delta_ten_t112
uv run python -m scripts.a6_delta_ten_contact_triple
uv run python -m scripts.a6_delta_ten_quadruple
uv run python -m scripts.a6_delta_ten_double_triple

uv run pytest -q \
  tests/test_a6_delta_ten_contact_wall.py \
  tests/test_a6_delta_ten_triple_wall.py \
  tests/test_a6_delta_ten_wall_components.py \
  tests/test_a6_delta_ten_partition_ledger.py \
  tests/test_a6_delta_ten_contact_three.py \
  tests/test_a6_delta_ten_double_contact.py \
  tests/test_a6_delta_ten_t112.py \
  tests/test_a6_delta_ten_contact_triple.py \
  tests/test_a6_delta_ten_quadruple.py \
  tests/test_a6_delta_ten_double_triple.py

uv run mypy --no-incremental
```

Regenerate both exact curves, singular schemes, and presentations with Sage
10.8:

```bash
sage tools/check_a6_delta_ten_contact_wall.sage
sage tools/check_a6_delta_ten_triple_wall.sage
sage tools/check_a6_delta_ten_contact_three.sage
sage tools/check_a6_delta_ten_double_contact.sage
sage tools/check_a6_delta_ten_t112.sage
sage tools/check_a6_delta_ten_contact_triple.sage
sage tools/check_a6_delta_ten_quadruple.sage
sage tools/check_a6_delta_ten_double_triple.sage
```

The Sage commands are manual checkers.  GitHub CI runs the Python
certificates and tests, not Sage, and neither route certifies the separate
proper Whitney--Thom propagation arguments.
