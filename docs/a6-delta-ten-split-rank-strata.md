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

## Completion boundary

Closing (1.1) for all 22 allocations proves that the true split fibers hide
no incidence component of dimension at least three in these six displayed
codimension-two profiles.  It does **not** compute their affine complement
groups, prove Whitney--Thom transport, classify intersections with removed
denominator or overlap charts, finish the remaining delta-ten profiles, or
settle generic degree six or `JC(2)`.
