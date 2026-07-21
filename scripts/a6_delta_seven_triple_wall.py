"""Exact local exclusion of the conditional delta-seven triple-fiber wall.

The degree-``(3, 10)`` normal form studied in
``a6_delta_seven_generic.py`` has a divisor ``T = 0`` on which three
normalization points acquire one common image.  This module normalizes that
divisor, checks its collision factorization exactly, and performs the finite
Hurwitz census needed for its local ``A6`` obstruction.

The decisive geometric input is deliberately exposed rather than hidden in
the computation: the three branches at a valid ``T``-point are all three
sheets of the trigonal projection ``P(t) = t^2 + t^3``.  Their three local
meridians must therefore generate the global monodromy image.  Without that
complete-fiber statement, the same local braid can be extended to a
surjective ``A6`` representation; ``HOSTILE_LOCAL_TRIPLE`` records an exact
counterfixture.

This remains part of a conditional, computer-assisted attack.  It does not
prove the unrestricted ``A6`` passport or the plane Jacobian conjecture.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache
from itertools import product
from typing import Final

from sympy import Expr, Symbol, discriminant, expand, resultant

from scripts.a6_delta_seven_generic import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    CUSP_COLLISION_FACTOR,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    GAMMA,
    S,
    TRIPLE_COLLISION_FACTOR,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    alternating_group_six,
    three_cycle_presentation_census,
)
from scripts.six_sheet_monodromy import (
    A6_COLLISION_MERIDIAN,
    A6_TORUS_2_5_LOCAL_GENERATORS,
    IDENTITY,
    Permutation,
    compose,
    cycle_type,
    generated_group,
    inverse,
)

TX: Final = Symbol("x")
TA: Final = Symbol("a")
TB: Final = Symbol("b")
TY: Final = Symbol("y")

T_WALL_PARAMETERS: Final = {
    ALPHA: TX * (2 * TA + TB + 2 * TX - 1),
    BETA: (
        -3 * TA * TX
        - 2 * TA
        - 3 * TB * TX
        - TB
        - 3 * TX**2
        + 3 * TX
        + 1
    ),
    GAMMA: (
        -3 * TA * TX
        - 3 * TA
        - 3 * TB * TX
        - 2 * TB
        - 3 * TX**2
        + 3 * TX
        + 3
    ),
    DELTA: -TA * TX - TA - TB * TX - TB - TX**2 + TX + 3,
}

TRIPLE_FIBER_FACTOR: Final = S**3 + 2 * S**2 + S + TX
RESIDUAL_COLLISION_FACTOR: Final = (
    S**4
    + 4 * S**3
    + (6 - TB) * S**2
    + (4 - TA - 2 * TB - TX) * S
    + 1
    - 2 * TA
    - TB
    - 2 * TX
)

RESIDUAL_DISCRIMINANT_FACTOR: Final = (
    27 * TA**3
    - 144 * TA**2 * TB
    + 81 * TA**2 * TX
    + 256 * TA**2
    - 4 * TA * TB**3
    + 128 * TA * TB**2
    - 288 * TA * TB * TX
    + 81 * TA * TX**2
    + 512 * TA * TX
    + 16 * TB**4
    - 4 * TB**3 * TX
    + 128 * TB**2 * TX
    - 144 * TB * TX**2
    + 27 * TX**3
    + 256 * TX**2
)

MIXED_RESULTANT_FACTOR: Final = (
    TA**3 * TX
    - 2 * TA**3
    - 2 * TA**2 * TB * TX
    - TA**2 * TB
    + 6 * TA**2 * TX**2
    - 10 * TA**2 * TX
    + TA**2
    - 2 * TA * TB**2 * TX
    - 8 * TA * TB * TX**2
    + 12 * TA * TX**3
    - 16 * TA * TX**2
    + 2 * TA * TX
    - TB**3 * TX**2
    - TB**2 * TX**2
    - 8 * TB * TX**3
    + TB * TX**2
    + 8 * TX**4
    - 8 * TX**3
    + TX**2
)

T_WALL_G_PULLBACK: Final = (
    (27 * TX - 4)
    * RESIDUAL_DISCRIMINANT_FACTOR
    * MIXED_RESULTANT_FACTOR**2
)

DOUBLE_TRIPLE_SUBSTITUTION: Final = {TA: -TX - TY, TB: 1}
SECOND_TRIPLE_FIBER_FACTOR: Final = S**3 + 2 * S**2 + S + TY
DOUBLE_TRIPLE_COLLISION_FACTOR: Final = (
    TRIPLE_FIBER_FACTOR * SECOND_TRIPLE_FIBER_FACTOR * (S + 2)
)
DOUBLE_TRIPLE_G_PULLBACK: Final = -(
    (27 * TX - 4)
    * (27 * TY - 4)
    * (TX - 2) ** 2
    * (TY - 2) ** 2
    * (TX - TY) ** 6
)

# Sage 10.8 affine presentations for one ordinary-triple member and one
# two-triple member.  The local proof below does not depend on these global
# sample calculations, but keeping them makes the wall audit independently
# replayable against the van Kampen implementation.
TRIPLE_WALL_GENERIC_RELATIONS: Final = (
    (3, 2, -3, -2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
    (2, 1, -2, -1),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, 3, 1),
    (3, 2, -3, -2),
    (2, 3, 1, -3, -2, -1),
    (2, -1, -3, -2, 3, 1),
)

DOUBLE_TRIPLE_RELATIONS: Final = (
    (2, 3, 1, -3, -2, -1),
    (2, -1, -3, -2, 3, 1),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, 3, 1),
    (2, 3, 1, -3, -2, -1),
    (2, -1, -3, -2, 3, 1),
    (3, 2, -3, -2),
)

ALLOWED_TRIPLE_CONTACT_PROFILES: Final = (
    (1, 1),
    (1, 2),
    (1, 3),
    (1, 4),
    (1, 5),
    (2, 2),
    (2, 3),
)

# The first three elements satisfy the local ordinary-triple braid, but only
# generate a subgroup of order nine.  The fourth element enlarges their image
# to all of A6.  This prevents silently using the local obstruction when the
# three branches are not the whole trigonal fiber.
HOSTILE_LOCAL_TRIPLE: Final = (
    A6_COLLISION_MERIDIAN,
    A6_COLLISION_MERIDIAN,
    A6_TORUS_2_5_LOCAL_GENERATORS[1],
)
HOSTILE_EXTENSION: Final = A6_TORUS_2_5_LOCAL_GENERATORS[0]


def _hurwitz_sigma(
    images: tuple[Permutation, Permutation, Permutation],
    index: int,
) -> tuple[Permutation, Permutation, Permutation]:
    """Apply the right Hurwitz generator ``sigma_(index + 1)``.

    The replacement ``(a, b) -> (b, b^-1 a b)`` preserves the ordered
    product.  Only ``B3`` is needed here, so ``index`` is zero or one.
    """

    if index not in (0, 1):
        msg = "a B3 Hurwitz generator index must be zero or one"
        raise ValueError(msg)
    result = list(images)
    left = images[index]
    right = images[index + 1]
    result[index] = right
    result[index + 1] = compose(compose(inverse(right), left), right)
    return result[0], result[1], result[2]


def hurwitz_action(
    images: tuple[Permutation, Permutation, Permutation],
    positive_word: Sequence[int],
) -> tuple[Permutation, Permutation, Permutation]:
    """Apply a positive ``B3`` word written with generator indices 1 and 2."""

    result = images
    for letter in positive_word:
        if letter not in (1, 2):
            msg = "the local braid word must use positive B3 generators"
            raise ValueError(msg)
        result = _hurwitz_sigma(result, letter - 1)
    return result


def triple_contact_braid(q: int, k: int) -> tuple[int, ...]:
    """Return ``Delta^(2q) sigma_1^(2(k-q))`` for profile ``(q,q,k)``."""

    if q < 1 or k < q:
        msg = "a three-branch contact profile requires 1 <= q <= k"
        raise ValueError(msg)
    half_twist = (1, 2, 1)
    return half_twist * (2 * q) + (1,) * (2 * (k - q))


@dataclass(frozen=True, slots=True)
class TripleContactCensus:
    """Fixed single-three-cycle tuples for one local contact profile."""

    profile: tuple[int, int, int]
    delta: int
    assignments: int
    fixed_assignments: int
    generated_order_histogram: tuple[tuple[int, int], ...]

    @property
    def a6_assignments(self) -> int:
        """Count fixed tuples that generate all of ``A6``."""

        return dict(self.generated_order_histogram).get(360, 0)


@cache
def triple_contact_census(q: int, k: int) -> TripleContactCensus:
    """Exhaust the 40^3 local meridian assignments for profile ``(q,q,k)``."""

    braid = triple_contact_braid(q, k)
    three_cycles = tuple(
        element
        for element in alternating_group_six()
        if cycle_type(element) == (3, 1, 1, 1)
    )
    orders: Counter[int] = Counter()
    fixed = 0
    for images in product(three_cycles, repeat=3):
        image_triple = (images[0], images[1], images[2])
        if hurwitz_action(image_triple, braid) != image_triple:
            continue
        fixed += 1
        orders[len(generated_group(image_triple))] += 1
    return TripleContactCensus(
        profile=(q, q, k),
        delta=2 * q + k,
        assignments=len(three_cycles) ** 3,
        fixed_assignments=fixed,
        generated_order_histogram=tuple(sorted(orders.items())),
    )


@dataclass(frozen=True, slots=True)
class TripleWallAlgebraCertificate:
    """Exact normalization and discriminant identities on ``T = 0``."""

    transformed_a_identity: Expr
    transformed_b_identity: Expr
    transformed_c_identity: Expr
    transformed_k_identity: Expr
    collision_factor_identity: Expr
    triple_divisor_identity: Expr
    alpha_factor_identity: Expr
    critical_factor_identity: Expr
    triple_factor_discriminant_identity: Expr
    residual_factor_discriminant_identity: Expr
    mixed_resultant_identity: Expr
    discriminant_wall_identity: Expr
    double_triple_factor_identity: Expr
    double_triple_c_identity: Expr
    double_triple_alpha_identity: Expr
    double_triple_critical_identity: Expr
    double_triple_discriminant_wall_identity: Expr

    @property
    def verified(self) -> bool:
        """Whether all polynomial identities vanish exactly."""

        return all(
            value == 0
            for value in (
                self.transformed_a_identity,
                self.transformed_b_identity,
                self.transformed_c_identity,
                self.transformed_k_identity,
                self.collision_factor_identity,
                self.triple_divisor_identity,
                self.alpha_factor_identity,
                self.critical_factor_identity,
                self.triple_factor_discriminant_identity,
                self.residual_factor_discriminant_identity,
                self.mixed_resultant_identity,
                self.discriminant_wall_identity,
                self.double_triple_factor_identity,
                self.double_triple_c_identity,
                self.double_triple_alpha_identity,
                self.double_triple_critical_identity,
                self.double_triple_discriminant_wall_identity,
            )
        )


@cache
def exact_triple_wall_algebra_certificate() -> TripleWallAlgebraCertificate:
    """Build the exact birational ``T``-wall certificate."""

    alpha = T_WALL_PARAMETERS[ALPHA]
    beta = T_WALL_PARAMETERS[BETA]
    gamma = T_WALL_PARAMETERS[GAMMA]
    delta = T_WALL_PARAMETERS[DELTA]
    cusp_factor = expand(CUSP_COLLISION_FACTOR.subs(T_WALL_PARAMETERS))
    wall_h = expand(COLLISION_POLYNOMIAL.subs(T_WALL_PARAMETERS))
    double_h = expand(wall_h.subs(DOUBLE_TRIPLE_SUBSTITUTION))
    double_c = expand(cusp_factor.subs(DOUBLE_TRIPLE_SUBSTITUTION))
    return TripleWallAlgebraCertificate(
        transformed_a_identity=expand(beta + 3 * delta - 2 * gamma - 4 - TA),
        transformed_b_identity=expand(gamma - 3 * delta + 6 - TB),
        transformed_c_identity=expand(cusp_factor - TX * (TA + TX)),
        transformed_k_identity=expand(
            alpha
            + 3 * beta
            + 5 * delta
            - 4 * gamma
            - 6
            - TA
            - (1 - TB) * TX
        ),
        collision_factor_identity=expand(
            wall_h - TRIPLE_FIBER_FACTOR * RESIDUAL_COLLISION_FACTOR
        ),
        triple_divisor_identity=expand(
            TRIPLE_COLLISION_FACTOR.subs(T_WALL_PARAMETERS)
        ),
        alpha_factor_identity=expand(alpha - TX * (2 * TA + TB + 2 * TX - 1)),
        critical_factor_identity=expand(
            EXTRA_CRITICAL_FACTOR.subs(T_WALL_PARAMETERS)
            - (27 * TX - 4) * (54 * TA + 9 * TB + 54 * TX - 1)
        ),
        triple_factor_discriminant_identity=expand(
            discriminant(TRIPLE_FIBER_FACTOR, S) + TX * (27 * TX - 4)
        ),
        residual_factor_discriminant_identity=expand(
            discriminant(RESIDUAL_COLLISION_FACTOR, S)
            + (TA + TX) * RESIDUAL_DISCRIMINANT_FACTOR
        ),
        mixed_resultant_identity=expand(
            resultant(TRIPLE_FIBER_FACTOR, RESIDUAL_COLLISION_FACTOR, S)
            - MIXED_RESULTANT_FACTOR
        ),
        discriminant_wall_identity=expand(
            discriminant(wall_h, S) - cusp_factor * T_WALL_G_PULLBACK
        ),
        double_triple_factor_identity=expand(
            double_h - DOUBLE_TRIPLE_COLLISION_FACTOR
        ),
        double_triple_c_identity=expand(double_c + TX * TY),
        double_triple_alpha_identity=expand(
            alpha.subs(DOUBLE_TRIPLE_SUBSTITUTION) + 2 * TX * TY
        ),
        double_triple_critical_identity=expand(
            EXTRA_CRITICAL_FACTOR.subs(T_WALL_PARAMETERS).subs(
                DOUBLE_TRIPLE_SUBSTITUTION
            )
            + 2 * (27 * TX - 4) * (27 * TY - 4)
        ),
        double_triple_discriminant_wall_identity=expand(
            discriminant(double_h, S)
            - double_c * DOUBLE_TRIPLE_G_PULLBACK
        ),
    )


@dataclass(frozen=True, slots=True)
class TripleWallLocalCertificate:
    """Finite braid obstruction for every valid delta-seven ``T`` point."""

    profiles: tuple[TripleContactCensus, ...]
    hostile_local_order: int
    hostile_extended_order: int
    hostile_braid_fixed: bool
    generic_complement_census: ThreeCyclePresentationCensus
    double_triple_complement_census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether the full profile list has no ``A6``-generating tuple."""

        return bool(
            tuple((c.profile[0], c.profile[2]) for c in self.profiles)
            == ALLOWED_TRIPLE_CONTACT_PROFILES
            and all(c.assignments == 40**3 for c in self.profiles)
            and all(c.delta <= 7 for c in self.profiles)
            and all(c.a6_assignments == 0 for c in self.profiles)
            and self.hostile_local_order == 9
            and self.hostile_extended_order == 360
            and self.hostile_braid_fixed
            and self.generic_complement_census.generated_order_histogram
            == ((3, 40),)
            and self.double_triple_complement_census.generated_order_histogram
            == ((3, 40),)
        )


@cache
def exact_triple_wall_local_certificate() -> TripleWallLocalCertificate:
    """Build the exhaustive local obstruction and hostile counterfixture."""

    return TripleWallLocalCertificate(
        profiles=tuple(
            triple_contact_census(q, k)
            for q, k in ALLOWED_TRIPLE_CONTACT_PROFILES
        ),
        hostile_local_order=len(generated_group(HOSTILE_LOCAL_TRIPLE)),
        hostile_extended_order=len(
            generated_group((*HOSTILE_LOCAL_TRIPLE, HOSTILE_EXTENSION))
        ),
        hostile_braid_fixed=(
            hurwitz_action(HOSTILE_LOCAL_TRIPLE, triple_contact_braid(1, 1))
            == HOSTILE_LOCAL_TRIPLE
        ),
        generic_complement_census=three_cycle_presentation_census(
            TRIPLE_WALL_GENERIC_RELATIONS
        ),
        double_triple_complement_census=three_cycle_presentation_census(
            DOUBLE_TRIPLE_RELATIONS
        ),
    )


def main() -> None:
    """Print the exact algebra and finite-braid checkpoint."""

    algebra = exact_triple_wall_algebra_certificate()
    local = exact_triple_wall_local_certificate()
    print("triple-wall algebra verified:", algebra.verified)
    print(
        "profile histograms:",
        {
            census.profile: dict(census.generated_order_histogram)
            for census in local.profiles
        },
    )
    print("hostile extension orders:", local.hostile_local_order, local.hostile_extended_order)
    print("entire valid triple wall excluded locally:", local.verified)


if __name__ == "__main__":
    main()
