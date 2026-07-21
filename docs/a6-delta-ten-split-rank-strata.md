# Exact target for the delta-ten split rank strata

## Scope

This checkpoint audits the coefficient-rank strata supported on the three
true reducible pair fibers

\[
  k=0,+2,-2
\]

of the conditional `(4,9)`, one-pair, single-three-cycle `A6` family.  It is
the determinantal continuation of the 22-row generic clean-witness ledger in
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

## Completion boundary

Closing (1.1) for all 22 allocations proves that the true split fibers hide
no incidence component of dimension at least three in these six displayed
codimension-two profiles.  It does **not** compute their affine complement
groups, prove Whitney--Thom transport, classify intersections with removed
denominator or overlap charts, finish the remaining delta-ten profiles, or
settle generic degree six or `JC(2)`.
