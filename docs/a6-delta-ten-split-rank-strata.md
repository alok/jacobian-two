# Exact target for the delta-ten split rank strata

## Scope

This checkpoint audits the coefficient-rank strata supported on the three
true reducible pair fibers

\[
  k=0,+2,-2
\]

of the conditional `(4,9)`, one-pair, single-three-cycle `A6` family.  It is
the determinantal continuation of the 22-orbit-type (33-actual-row) generic
clean-witness ledger in
[`a6-delta-ten-codim-two.md`](a6-delta-ten-codim-two.md).

For each allowed allocation, let `B_r` denote the valid base on which the
affine-linear incidence system in `(a,b,c,d)` has coefficient rank `r` and
the same augmented rank.  The required dimension bound is

\[
  \dim B_r + 4-r \le 2.                         \tag{1.1}
\]

An incompatible locus, where the augmented rank is larger than the
coefficient rank, contributes no incidence.  Every use of (1.1) must be
proved after localizing the intended allocation away from coincident roots,
cusp fibers, component overlaps not prescribed by the allocation, repeated
special targets, extra critical sources, and profile-changing higher
collisions.

## Exact work packages

1. `C3+7N` and `C2^2+6N`: classify the maximal-minor and augmented-minor
   strata for every `V`, `W`, `VW`, `WW`, and legitimate `k=0` overlap
   allocation.
2. `T112+6N` and `C2+T111+5N`: work on labeled four-source fibers, distinguish
   which pair edge acquires the extra contact order, and saturate the rank
   loci away from `P`-critical and fourth-source boundaries.
3. `Q0+4N` and `T111^2+4N`: specialize the already global fiber-remainder
   incidence systems.  For `Q0`, prove rank three on every valid nonzero
   fiber.  For two triples, prove that the `k=0` residual rank factor is
   contained in the same-`P`-fiber boundary and that the `k=+/-2` residual
   factor has no compatible point on its valid localization.
4. Replay the full family involution

   \[
     (k,t,Y,b,d)\longmapsto(-k,-t,-Y,-b,-d)
   \]

   on every `k=2` statement before promoting it to `k=-2`.

Each promoted statement needs exact symbolic identities, an independently
checked witness on every nonempty clean stratum, and a hostile excluded
fixture whenever dropping a localizer would create a false component.

## The global-fiber profiles close on the split values

The first exact closure is available without constructing a new split-only
incidence.

### Ordinary quadruple

For `Q0+4N`, reduce `Q` modulo the complete fiber quartic `P-h`.  The four
maximal minors of the resulting `3 x 4` coefficient matrix are

```text
h^2*(h+1),
-h^2*k*(2*h+1),
h^2*(h*k^2+2*h+1),
-2*h^3*k.
```

After division by the common invalid factor `h^2`, Rabinowitsch localization
at `h != 0` gives the unit ideal separately at `k=0,+2,-2`.  Hence the rank
is exactly three on every valid split ordinary-quadruple fiber.  At `h=0`
the rank really falls to two, but

```text
disc_t(P-h)|k=0  = -16*h*(4*h+1)^2,
disc_t(P-h)|k=±2 = -16*h^2*(16*h-1),
```

so the hostile rank drop is a singular source fiber rather than `Q0`.
The fixed-`k` incidence has a one-dimensional `h` base and an affine-line
coefficient fiber, hence dimension two.

Moreover this split surface is not a new topology component.  Before fixing
`k`, the rank-three incidence is an affine-line bundle over an irreducible
open of the `(k,h)` plane.  Its exact clean locus is therefore irreducible
and analytically path-connected and contains its clean split points.  The
existing simultaneous-resolution and proper-isotopy argument transports the
cyclic ordinary-quadruple sample to them.

### Two ordinary triples

For `T111^2+4N`, use the two omitted roots `u,v`.  At `k=0`, the remaining
rank factor and the same-fiber factor specialize to

```text
A = -u-v,
B = (u+v)*(u^2+v^2+1).
```

Thus `A=0` is contained in `B=0`, where the two proposed triple targets are
actually one four-source fiber.  It contributes no two-triple member.  At
`k=+2` the residual factor is `2*u*v-u-v`; at `k=-2` it is
`-2*u*v-u-v`.  The global residual certificate solves
`k=(u+v)/(u*v)` and proves that the four normalized augmented determinants,
localized away from coincident roots, cusp fibers, and `B=0`, generate the
unit ideal.  Therefore neither split residual factor carries a compatible
valid point.

Two hostile fixtures make the boundary use explicit.  The compatible
`k=0, v=-u` family lies on `B=0`.  At

\[
 (k,u,v)=\left(2,-1,\frac13\right)
\]

the coefficient and augmented ranks are both three and the exact member

\[
 (a,b,c,d)=
 \left(\frac{16}{27},-\frac{646}{297},-\frac{373}{99},0\right)
\]

solves all four two-triple equations, but its first cusp-fiber factor is
zero.  Dropping either localizer would therefore create a false split
component.

Every valid split two-triple point consequently lies in the determinant-open
Cramer graph over an irreducible open of `(k,u,v)`.  Its fixed-`k` base has
dimension two and its coefficient fiber is a point.  As for `Q0`, the clean
split locus lies in the same path-connected global clean incidence as the
exact cyclic sample, so the existing proper-isotopy conclusion applies.

The executable certificate is
[`a6_delta_ten_split_rank_q0_t111.py`](../scripts/a6_delta_ten_split_rank_q0_t111.py).
It closes both the rank and topology obligations for the clean split
`Q0+4N` and `T111^2+4N` rows; intersections with removed singular-fiber,
cusp-fiber, same-target, or deeper-degeneration boundaries remain separate.

## Contact-only profiles

The five `C3+7N` and six `C2^2+6N` allocation systems are rebuilt from the
true vertical and graph factors.  Every visible rank loss lies on an exact
cusp, diagonal, overlap, or repeated-contact factor except for three
irreducible two-contact residual determinants.  After valid localization,
their compatible bases are reduced of ordered lengths

```text
k=0, WW: 4
k=2, VW: 6
k=2, WW: 6.
```

The coefficient rank is exactly three on all sixteen representative residual base points,
so their affine-line fibers have incidence dimension one.  Rank at most two
saturates to the unit ideal.  The sole compatible raw `C3` rank loss uses a
pair containing the forced cusp source and has collision order at least four,
so it is not a clean `C3` point.  This closes (1.1) for all eleven contact
allocations.

The rank-open topology is controlled by the total pair equation

\[
 F=(2s+k)p-s(s^2+ks+1)=0
\]

and its fiber-tangent derivation `D=F_p*d/ds-F_s*d/dp`.  Exact invertible row
transformations compare `C,D(C),D^2(C)` with both the historical
`H,H',H''` rows and the component contact jets.  Localized power-membership
identities prove that every transformation determinant is a unit on the full
declared clean chart.  The
total pair surface is geometrically integral and flat over the `k` line; its
ordered self-fiber-product is integral by flatness plus generic-fiber
integrality, and its singular support is exactly the three removed
double-overlap points.  Hence the three maximal-rank `C3` and five
maximal-rank `C2^2` split loci lie in the already-excluded global incidences.
Exact rational arcs prove that the two `C3` overlap allocations and the
`C2^2` overlap-plus-contact allocation are algebraically contained in those
components as well.  Their boundary topology and the sixteen representative
ordered residual affine-line fibers remain open.

The executable algebra and closure certificates are
[`a6_delta_ten_split_contact_rank.py`](../scripts/a6_delta_ten_split_contact_rank.py)
and
[`a6_delta_ten_split_contact_closure.py`](../scripts/a6_delta_ten_split_contact_closure.py).

## Triple and mixed profiles

On the labeled split four-source curves

\[
 (z,-z,w,-w),\quad z^2+w^2+1=0
\]

at `k=0`, and

\[
 (z,-1-z,w,-1-w),\quad z(z+1)+w(w+1)=0
\]

at `k=2`, all four `T112` systems have rank exactly three with no valid rank
drop.  The three mixed systems have generic rank four.  Their compatible
rank-three bases have reduced lengths `4,4,0`; rank at most two is empty.
Thus all seven generic split incidences have dimension two and every
exceptional compatible mixed incidence has dimension one.

For topology, the `T112` rows restrict from the existing irreducible labeled
incidence by exact invertible transformations.  The mixed vertical boundary
is captured by retaining contact sum and product simultaneously.  If
`n=r^2+q^2-q^3`, its total base equation is

\[
 (2qrs+n)p-s\bigl(qr(s^2+1)+ns\bigr)=0.
\]

The two coefficients in `p` have gcd one, so this is geometrically
irreducible; its singular ideal is empty after localizing by `q*r`.  An exact
four-row diagonal bridge recovers the historical mixed Cramer system, and
Rabinowitsch saturations prove that every split transformation determinant is
a unit on its full valid localization.  All three maximal-rank mixed split
loci are therefore in the same global Cramer component as the cyclic sample.
Only the two representative length-four rank-three affine-line schemes—three
after `k=-2` transport—remain topologically open.

The executable certificates are
[`a6_delta_ten_split_t112_mixed_rank.py`](../scripts/a6_delta_ten_split_t112_mixed_rank.py)
and
[`a6_delta_ten_split_component_closure.py`](../scripts/a6_delta_ten_split_component_closure.py).

## Exact all-allocation coverage

The aggregate
[`a6_delta_ten_split_rank_all.py`](../scripts/a6_delta_ten_split_rank_all.py)
derives allocation keys from the focused systems and compares them with the
generated ledger.  There are exactly 22 distinct involution-orbit keys, with
profile counts

```text
C3+7N:          5
C2^2+6N:        6
T112+6N:        4
C2+T111+5N:     3
Q0+4N:          2
T111^2+4N:      2.
```

The eleven `k=+/-2` keys each represent two instantiated rows, so the audit
covers 33 actual rows over `k=0,+2,-2`.  The sets agree exactly.  Therefore every true-split rank stratum in the six
displayed codimension-two profiles satisfies (1.1); this is not inferred
from a count or expected dimension.

Singular independently checks that the total pair, ordered two-pair, and
total mixed base ideals are prime of dimensions `2,3,3`:

```bash
sage tools/check_a6_delta_ten_split_component_closure.sage
```

## Completion boundary

Closing (1.1) for all 22 orbit types, hence all 33 actual rows, proves that the true split fibers hide
no incidence component of dimension at least three in these six displayed
codimension-two profiles.  Every clean rank-open split locus is also connected
to an already-excluded global incidence.  All three prescribed overlap
incidences are algebraically contained too, but the result does **not** yet
exclude their boundary topology or the five finite compatible orbit
representatives (eight actual rank-three schemes); those require topology.
It
also does not classify removed non-clean denominator/intersection charts,
finish the remaining delta-ten profiles, or settle generic degree six or
`JC(2)`.
