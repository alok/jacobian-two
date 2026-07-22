"""Certify the conditional delta-ten ``C3 + T111 + 4N`` profile.

In the normalized degree-``(4, 9)`` family

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``,

choose an omitted root ``e`` of an ordinary four-point ``P``-fiber.  The
other three roots form ``(P(t)-P(e))/(t-e)``; requiring ``Q`` to be constant
on them gives two affine-linear equations in ``(a,b,c,d)``.  A separate
contact of order three at unordered-pair sum ``w`` gives the three equations
``H(w)=H'(w)=H''(w)=0``.

The five-by-five augmented determinant is

``2*C(e)^2*C(w)^3*O(k,e,w)^3*R(k,e,w)``,

where the first three factors are removed fiber/overlap boundaries and ``R``
is an irreducible 62-term degree-nine polynomial.  Exact Sage saturation by
the complete displayed valid-chart localizer proves that ``R`` is smooth on
the valid chart.  It finds a degree-fourteen coefficient-rank-drop curve but
proves that augmented rank at most three is empty there; redundant rank-two
saturations are units as well.  Hence every valid compatible base has a
unique coefficient vector.

The smooth rational base ``(k,e,w)=(-4,-1/2,1)`` gives

``(a,b,c,d)=(39/2,-409/8,109/4,-31/4)``.

The resulting curve has one exact contact three, one ordinary triple point,
four nodes, the forced ``T(2,5)`` cusp, and the fixed ``T(5,9)`` infinity
branch.  Its exact Sage van Kamp presentation is stored below once generated;
the independent finite replay tests every ``40^4`` single-three-cycle image.

This is a conditional, computer-assisted profile certificate.  Split,
singular-fiber, overlap/same-target, non-clean, and deeper boundary charts
remain open.  Nothing here proves the plane Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from itertools import combinations
from typing import Final

from sympy import (
    Expr,
    Matrix,
    Poly,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    gcd,
    invert,
    linear_eq_to_matrix,
    rem,
    resultant,
    together,
)

from scripts.a6_delta_ten_double_triple import (
    GENERIC_OMITTED_CUBIC,
    GENERIC_TRIPLE_EQUATIONS,
    OMITTED_ROOT,
)
from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    COLLISION_X_DENOMINATOR,
    COLLISION_X_NUMERATOR,
    CUSP_IMAGE_FACTOR,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    PAIR_DENOMINATOR,
    PAIR_DIAGONAL_FACTOR,
    PAIR_QUADRATIC,
    S,
    T,
    TANGENCY_POLYNOMIAL,
    X,
    Y,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    n_generator_three_cycle_presentation_census,
)

TRIPLE_OMITTED_ROOT: Final = OMITTED_ROOT
CONTACT_THREE_SUM: Final = Symbol("w")
TRIPLE_SLOPE: Final = Symbol("m_three_triple")
LINEAR_PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)


def _contact_jets(root: Expr) -> tuple[Expr, ...]:
    """Return the three collision jets imposing contact order three."""

    return tuple(
        diff(COLLISION_POLYNOMIAL, S, order).subs(S, root)
        for order in range(3)
    )


CONTACT_THREE_TRIPLE_EQUATIONS: Final = (
    *_contact_jets(CONTACT_THREE_SUM),
    *GENERIC_TRIPLE_EQUATIONS,
)
(
    CONTACT_THREE_TRIPLE_COEFFICIENT_MATRIX,
    CONTACT_THREE_TRIPLE_RIGHT_HAND_SIDE,
) = linear_eq_to_matrix(CONTACT_THREE_TRIPLE_EQUATIONS, LINEAR_PARAMETERS)
CONTACT_THREE_TRIPLE_AUGMENTED_MATRIX: Final = (
    CONTACT_THREE_TRIPLE_COEFFICIENT_MATRIX.row_join(
        CONTACT_THREE_TRIPLE_RIGHT_HAND_SIDE,
    )
)

TRIPLE_CUSP_FIBER_FACTOR: Final = (
    TRIPLE_OMITTED_ROOT**2 + KAPPA * TRIPLE_OMITTED_ROOT + 1
)
TRIPLE_REDUCED_DISCRIMINANT: Final = (
    16 * TRIPLE_OMITTED_ROOT**4
    + 8 * TRIPLE_OMITTED_ROOT**3 * KAPPA
    - 5 * TRIPLE_OMITTED_ROOT**2 * KAPPA**2
    + 16 * TRIPLE_OMITTED_ROOT**2
    + 3 * TRIPLE_OMITTED_ROOT * KAPPA**3
    - 12 * TRIPLE_OMITTED_ROOT * KAPPA
    - KAPPA**2
    + 4
)
TRIPLE_OMITTED_CRITICAL_FACTOR: Final = (
    4 * TRIPLE_OMITTED_ROOT**2
    + 3 * KAPPA * TRIPLE_OMITTED_ROOT
    + 2
)
CONTACT_PAIR_DENOMINATOR: Final = KAPPA + 2 * CONTACT_THREE_SUM
CONTACT_PAIR_QUADRATIC: Final = (
    CONTACT_THREE_SUM**2 + KAPPA * CONTACT_THREE_SUM + 1
)
CONTACT_PAIR_DIAGONAL: Final = (
    2 * CONTACT_THREE_SUM**2
    + 3 * KAPPA * CONTACT_THREE_SUM
    + 4
)

CONTACT_TRIPLE_OVERLAP_B: Final = (
    TRIPLE_OMITTED_ROOT**2 * KAPPA
    + 2 * TRIPLE_OMITTED_ROOT**2 * CONTACT_THREE_SUM
    - TRIPLE_OMITTED_ROOT * KAPPA * CONTACT_THREE_SUM
    - 2 * TRIPLE_OMITTED_ROOT * CONTACT_THREE_SUM**2
    + KAPPA * CONTACT_THREE_SUM**2
    + CONTACT_THREE_SUM**3
    + CONTACT_THREE_SUM
)
CONTACT_TRIPLE_OVERLAP_O: Final = (
    TRIPLE_OMITTED_ROOT**2 * KAPPA
    + 2 * TRIPLE_OMITTED_ROOT**2 * CONTACT_THREE_SUM
    + TRIPLE_OMITTED_ROOT * KAPPA**2
    + 3 * TRIPLE_OMITTED_ROOT * KAPPA * CONTACT_THREE_SUM
    + 2 * TRIPLE_OMITTED_ROOT * CONTACT_THREE_SUM**2
    + KAPPA**2 * CONTACT_THREE_SUM
    + 2 * KAPPA * CONTACT_THREE_SUM**2
    + KAPPA
    + CONTACT_THREE_SUM**3
    + CONTACT_THREE_SUM
)

e = TRIPLE_OMITTED_ROOT
k = KAPPA
w = CONTACT_THREE_SUM
CONTACT_THREE_TRIPLE_COMPATIBILITY: Final = (
    6 * e**3 * k**5 * w
    + 42 * e**3 * k**4 * w**2
    + 22 * e**3 * k**4
    + 108 * e**3 * k**3 * w**3
    + 156 * e**3 * k**3 * w
    + 120 * e**3 * k**2 * w**4
    + 408 * e**3 * k**2 * w**2
    + 48 * e**3 * k * w**5
    + 464 * e**3 * k * w**3
    + 192 * e**3 * w**4
    - 3 * e**2 * k**6 * w
    - 36 * e**2 * k**5 * w**2
    - 11 * e**2 * k**5
    - 149 * e**2 * k**4 * w**3
    - 141 * e**2 * k**4 * w
    - 276 * e**2 * k**3 * w**4
    - 600 * e**2 * k**3 * w**2
    - 28 * e**2 * k**3
    - 240 * e**2 * k**2 * w**5
    - 1120 * e**2 * k**2 * w**3
    - 144 * e**2 * k**2 * w
    - 80 * e**2 * k * w**6
    - 984 * e**2 * k * w**4
    - 192 * e**2 * k * w**2
    - 8 * e**2 * k
    - 336 * e**2 * w**5
    - 80 * e**2 * w**3
    + 3 * e * k**6 * w**2
    + 28 * e * k**5 * w**3
    + 15 * e * k**5 * w
    + 87 * e * k**4 * w**4
    + 147 * e * k**4 * w**2
    + 12 * e * k**4
    + 126 * e * k**3 * w**5
    + 468 * e * k**3 * w**3
    + 126 * e * k**3 * w
    + 88 * e * k**2 * w**6
    + 696 * e * k**2 * w**4
    + 336 * e * k**2 * w**2
    + 16 * e * k**2
    + 24 * e * k * w**7
    + 504 * e * k * w**5
    + 360 * e * k * w**3
    + 24 * e * k * w
    + 144 * e * w**6
    + 144 * e * w**4
    - 3 * k**5 * w**2
    - 28 * k**4 * w**3
    - 6 * k**4 * w
    - 87 * k**3 * w**4
    - 54 * k**3 * w**2
    - 3 * k**3
    - 126 * k**2 * w**5
    - 132 * k**2 * w**3
    - 30 * k**2 * w
    - 88 * k * w**6
    - 132 * k * w**4
    - 48 * k * w**2
    - 4 * k
    - 24 * w**7
    - 48 * w**5
    - 24 * w**3
)
EXPECTED_AUGMENTED_DETERMINANT: Final = expand(
    2
    * TRIPLE_CUSP_FIBER_FACTOR**2
    * CONTACT_PAIR_QUADRATIC**3
    * CONTACT_TRIPLE_OVERLAP_O**3
    * CONTACT_THREE_TRIPLE_COMPATIBILITY
)


def _cleared_contact_pair() -> Expr:
    """Return the contact pair's source quadratic with denominator cleared."""

    return expand(
        CONTACT_PAIR_DENOMINATOR * T**2
        - CONTACT_THREE_SUM * CONTACT_PAIR_DENOMINATOR * T
        + CONTACT_THREE_SUM * CONTACT_PAIR_QUADRATIC
    )


CONTACT_THREE_CLEARED_PAIR: Final = _cleared_contact_pair()
CONTACT_TRIPLE_SOURCE_RESULTANT_IDENTITY: Final = expand(
    resultant(GENERIC_OMITTED_CUBIC, CONTACT_THREE_CLEARED_PAIR, T)
    - CONTACT_TRIPLE_OVERLAP_B * CONTACT_TRIPLE_OVERLAP_O**2
)
CONTACT_TRIPLE_TARGET_SEPARATION_IDENTITY: Final = expand(
    FAMILY_P.subs(T, TRIPLE_OMITTED_ROOT)
    * COLLISION_X_DENOMINATOR.subs(S, CONTACT_THREE_SUM)
    - COLLISION_X_NUMERATOR.subs(S, CONTACT_THREE_SUM)
    - CONTACT_TRIPLE_OVERLAP_B * CONTACT_TRIPLE_OVERLAP_O
)

CONTACT_THREE_TRIPLE_VALID_LOCALIZER: Final = expand(
    KAPPA
    * (KAPPA - 2)
    * (KAPPA + 2)
    * TRIPLE_OMITTED_ROOT
    * CONTACT_THREE_SUM
    * CONTACT_PAIR_DENOMINATOR
    * CONTACT_PAIR_QUADRATIC
    * CONTACT_PAIR_DIAGONAL
    * TRIPLE_CUSP_FIBER_FACTOR
    * TRIPLE_OMITTED_CRITICAL_FACTOR
    * TRIPLE_REDUCED_DISCRIMINANT
    * CONTACT_TRIPLE_OVERLAP_B
    * CONTACT_TRIPLE_OVERLAP_O
)

# Reproducible Sage 10.8 metadata.  The independent checker derives the
# determinant residual and every localized saturation from the matrices.
SAGE_DETERMINANT_FACTOR_DATA: Final = (2, 2, 3, 3, 1)
SAGE_COMPATIBILITY_RESIDUAL_SHAPE: Final = (62, 9, 6, 3, 7)
SAGE_COMPATIBILITY_RATIONAL_FACTOR_COUNT: Final = 1
SAGE_COMPATIBILITY_SAMPLE_GRADIENT: Final = (-84, -112, 0)
SAGE_VALID_SINGULAR_SATURATION: Final = (True, 1)
SAGE_COEFFICIENT_RANK_THREE_CURVE: Final = (1, 14, -21, 1)
SAGE_AUGMENTED_RANK_THREE_UNIT: Final = (True, 1)
SAGE_RANK_TWO_UNIT_EXPONENTS: Final = (1, 1)
SAGE_NORMALIZED_MINOR_COUNTS: Final = (5, 25, 40, 100)

SAMPLE_BASE_POINT: Final = {
    KAPPA: -4,
    TRIPLE_OMITTED_ROOT: Rational(-1, 2),
    CONTACT_THREE_SUM: 1,
}
SAMPLE_PARAMETERS: Final = {
    KAPPA: -4,
    ALPHA: Rational(39, 2),
    BETA: Rational(-409, 8),
    GAMMA: Rational(109, 4),
    DELTA: Rational(-31, 4),
}
SAMPLE_P: Final = expand(FAMILY_P.subs(SAMPLE_PARAMETERS))
SAMPLE_Q: Final = expand(FAMILY_Q.subs(SAMPLE_PARAMETERS))
SAMPLE_Y: Final = expand(8 * SAMPLE_Q)
SAMPLE_TRIPLE_CUBIC: Final = expand(
    GENERIC_OMITTED_CUBIC.subs(SAMPLE_BASE_POINT),
)
SAMPLE_TRIPLE_VALUE_X: Final = Rational(13, 16)
SAMPLE_TRIPLE_VALUE_Q: Final = Rational(-2197, 1024)
SAMPLE_TRIPLE_QUOTIENT: Final = (
    T**6
    - Rational(13, 4) * T**5
    + Rational(75, 8) * T**4
    + Rational(13, 4) * T**3
    - Rational(13, 8) * T**2
    - Rational(169, 64) * T
    - Rational(169, 128)
)
SAMPLE_TRIPLE_SLOPE_POLYNOMIAL: Final = (
    1179648 * TRIPLE_SLOPE**3
    - 643069952 * TRIPLE_SLOPE**2
    - 6362836480 * TRIPLE_SLOPE
    - 15599304175
)

SAMPLE_CONTACT_BRANCH: Final = T**2 - T + 1
SAMPLE_TRIPLE_PAIR_SUMS: Final = 2 * S**3 - 18 * S**2 + 47 * S - 26
SAMPLE_RESIDUAL_NODE_SUMS: Final = S**4 - 12 * S**3 + 24 * S**2 + 7 * S + 24
SAMPLE_COLLISION: Final = expand(
    Rational(1, 2)
    * (S - 1) ** 3
    * SAMPLE_TRIPLE_PAIR_SUMS
    * SAMPLE_RESIDUAL_NODE_SUMS
)
SAMPLE_TANGENCY_COFACTOR: Final = (
    18 * S**9
    - 432 * S**8
    + 4033 * S**7
    - 18938 * S**6
    + 47963 * S**5
    - 64986 * S**4
    + 47578 * S**3
    - 31966 * S**2
    + 21980 * S
    - 6240
)
SAMPLE_TANGENCY: Final = expand(2 * (S - 1) ** 2 * SAMPLE_TANGENCY_COFACTOR)

CONTACT_THREE_TRIPLE_NODE_X_POLYNOMIAL: Final = (
    16 * X**4
    + 9324 * X**3
    - 60343 * X**2
    + 117047 * X
    - 60306
)
CONTACT_THREE_TRIPLE_IMPLICIT: Final = (
    -4096 * X**9
    + 479504 * X**8
    - 49932856 * X**7
    - 6816 * X**6 * Y
    + 31532917 * X**6
    - 6558016 * X**5 * Y
    + 7276464 * X**5
    - 2600 * X**4 * Y**2
    + 4524178 * X**4 * Y
    - 288810 * X**3 * Y**2
    + 688298 * X**3 * Y
    - 40 * X**2 * Y**3
    + 218553 * X**2 * Y**2
    - 4262 * X * Y**3
    + 16774 * X * Y**2
    + Y**4
    + 3554 * Y**3
    - 299 * Y**2
)

# Exact Sage 10.8 affine van Kamp presentation.  Generator indices 1..4 are
# geometric meridians of a generic vertical fiber.
CONTACT_THREE_TRIPLE_RELATIONS: Final[tuple[tuple[int, ...], ...]] = (
    (2, 1, -2, -1),
    (-2, -1, -4, -3, 1, 3, 4, 1),
    (
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        1,
        3,
        4,
        3,
        -4,
        -3,
        -1,
        3,
        4,
        -3,
        1,
        3,
        4,
        3,
        -4,
        -3,
        -1,
        3,
        4,
        -3,
        1,
        3,
        4,
        3,
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        1,
        3,
        4,
        -3,
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        1,
        3,
        4,
        -3,
    ),
    (
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        1,
        3,
        4,
        -3,
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        1,
        3,
        4,
        -3,
        1,
        3,
        4,
        3,
    ),
    (-3, -2, -1, 3, -4, -3, 1, 2, 3, 4),
    (-2, -1, 2, 3, 4, -3, 1, 3, -4, -3),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        1,
        3,
        4,
        3,
        -4,
        -3,
        -1,
        3,
        4,
        -3,
        1,
        3,
        4,
        -3,
    ),
    (
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        -1,
        3,
        4,
        -3,
        1,
        3,
        4,
        3,
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        1,
        3,
        4,
        -3,
        1,
        3,
        4,
        -3,
    ),
    (
        -4,
        -3,
        1,
        3,
        4,
        -3,
        1,
        3,
        4,
        -3,
        1,
        3,
        4,
        -3,
        -1,
        3,
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        -1,
        3,
    ),
)

SAMPLE_BASE_TANGENTS: Final = ((4, -3, 0), (0, 0, 1))
SAMPLE_COEFFICIENT_TANGENTS: Final = (
    (
        Rational(4147, 14),
        Rational(-20431, 28),
        337,
        Rational(-571, 14),
    ),
    (
        Rational(7293, 28),
        Rational(-4686, 7),
        297,
        Rational(-1023, 28),
    ),
)


def _modulo_branch(rational_expression: Expr, branch: Expr) -> Expr:
    """Reduce a rational function modulo a squarefree source branch."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    return expand(rem(numerator * invert(denominator, branch), branch, T))


def _contact_jet_differences(order: int) -> tuple[Expr, ...]:
    """Compare successive ``dQ/dP`` jets on the two contact branches."""

    current = SAMPLE_Q
    p_derivative = diff(SAMPLE_P, T)
    differences: list[Expr] = []
    for _ in range(order):
        current = cancel(diff(current, T) / p_derivative)
        first = _modulo_branch(current, SAMPLE_CONTACT_BRANCH)
        second = _modulo_branch(
            current.subs(T, 1 - T),
            SAMPLE_CONTACT_BRANCH,
        )
        differences.append(expand(rem(first - second, SAMPLE_CONTACT_BRANCH, T)))
    return tuple(differences)


@dataclass(frozen=True, slots=True)
class DeltaTenContactThreeTripleCertificate:
    """Exact incidence, sample geometry, and complement certificate."""

    source_overlap_identity: Expr
    target_separation_identity: Expr
    determinant_factor_data: tuple[int, int, int, int, int]
    compatibility_residual_shape: tuple[int, int, int, int, int]
    compatibility_rational_factor_count: int
    compatibility_sample_value: Expr
    compatibility_sample_gradient: tuple[Expr, Expr, Expr]
    sage_valid_singular_saturation: tuple[bool, int]
    sage_coefficient_rank_three_curve: tuple[int, int, int, int]
    sage_augmented_rank_three_unit: tuple[bool, int]
    sage_rank_two_unit_exponents: tuple[int, int]
    sage_normalized_minor_counts: tuple[int, int, int, int]
    sample_coefficient_rank: int
    sample_augmented_rank: int
    sample_incidence_residuals: tuple[Expr, ...]
    sample_maximal_minors: tuple[Expr, ...]
    sample_valid_localizer: Expr
    tangent_incidence_residuals: tuple[tuple[Expr, ...], ...]
    coefficient_image_tangent_rank: int
    triple_cubic_discriminant: Expr
    triple_p_derivative_resultant: Expr
    triple_omitted_root_separation: Expr
    triple_q_remainder: Expr
    triple_omitted_q_difference: Expr
    triple_quotient_identity: Expr
    triple_slope_eliminant_identity: Expr
    triple_slope_discriminant: Expr
    contact_pair_discriminant: Expr
    contact_image_remainders: tuple[Expr, Expr]
    contact_p_derivative_resultant: Expr
    contact_jet_differences: tuple[Expr, ...]
    contact_jet_separation: Expr
    collision_identity: Expr
    tangency_identity: Expr
    collision_tangency_gcd: Expr
    collision_jets_at_contact: tuple[Expr, ...]
    tangency_jets_at_contact: tuple[Expr, ...]
    collision_factor_discriminants: tuple[Expr, Expr]
    collision_factor_resultants: tuple[Expr, Expr, Expr]
    triple_boundary_resultants: tuple[Expr, Expr, Expr, Expr]
    node_boundary_resultants: tuple[Expr, Expr, Expr, Expr]
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    triple_x_eliminant_identity: Expr
    node_x_eliminant_identity: Expr
    node_x_discriminant: Expr
    node_target_separations: tuple[Expr, Expr, Expr]
    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    implicit_content: Expr
    sage_jacobian_length: int
    sage_jacobian_radical_length: int
    sage_jacobian_components: tuple[tuple[int, int], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    topology_propagation_dependencies: tuple[str, str]
    split_and_deeper_boundaries_classified_here: bool
    arithmetic_genus: int
    cusp_delta: int
    contact_three_delta: int
    ordinary_triple_delta: int
    node_count: int
    infinity_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def compatibility_sample_is_smooth(self) -> bool:
        """Whether the residual vanishes simply at the rational sample."""

        return self.compatibility_sample_value == 0 and any(
            value != 0 for value in self.compatibility_sample_gradient
        )

    @property
    def compatibility_surface_geometrically_irreducible(self) -> bool:
        """Apply irreducibility plus the smooth-rational-point criterion."""

        return (
            self.compatibility_rational_factor_count == 1
            and self.compatibility_sample_is_smooth
        )

    @property
    def valid_rank_drop_incidence_empty(self) -> bool:
        """Whether all valid compatible bases have coefficient rank four."""

        return bool(
            self.sage_coefficient_rank_three_curve == (1, 14, -21, 1)
            and self.sage_augmented_rank_three_unit == (True, 1)
            and self.sage_rank_two_unit_exponents == (1, 1)
        )

    @property
    def valid_compatibility_surface_smooth(self) -> bool:
        """Whether exact full-localizer singular saturation is the unit ideal."""

        return self.sage_valid_singular_saturation == (True, 1)

    @property
    def dominant_incidence_dimension(self) -> int:
        """Return the verified coefficient-image tangent dimension."""

        return self.coefficient_image_tangent_rank

    @property
    def total_delta(self) -> int:
        """Return the projective singularity-genus contribution."""

        return (
            self.cusp_delta
            + self.contact_three_delta
            + self.ordinary_triple_delta
            + self.node_count
            + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether every exact incidence and topology check agrees."""

        return bool(
            self.source_overlap_identity == 0
            and self.target_separation_identity == 0
            and self.determinant_factor_data == (2, 2, 3, 3, 1)
            and self.compatibility_residual_shape == (62, 9, 6, 3, 7)
            and self.compatibility_surface_geometrically_irreducible
            and self.valid_compatibility_surface_smooth
            and self.valid_rank_drop_incidence_empty
            and self.sage_normalized_minor_counts == (5, 25, 40, 100)
            and self.sample_coefficient_rank == 4
            and self.sample_augmented_rank == 4
            and self.sample_incidence_residuals == (0, 0, 0, 0, 0)
            and self.sample_maximal_minors
            == (-2504320, 3333120, 0, Rational(29575, 4), Rational(159705, 4))
            and self.sample_valid_localizer == 7223580
            and all(
                all(value == 0 for value in residuals)
                for residuals in self.tangent_incidence_residuals
            )
            and self.dominant_incidence_dimension == 2
            and self.triple_cubic_discriminant == Rational(-637, 4)
            and self.triple_p_derivative_resultant == Rational(5733, 8)
            and self.triple_omitted_root_separation == Rational(9, 2)
            and self.triple_q_remainder == Rational(-2197, 1024)
            and self.triple_omitted_q_difference == Rational(63, 128)
            and self.triple_quotient_identity == 0
            and self.triple_slope_eliminant_identity == 0
            and self.triple_slope_discriminant
            == -89062908728555597524369408000000
            and self.contact_pair_discriminant == -3
            and self.contact_image_remainders == (3, Rational(-199, 8))
            and self.contact_p_derivative_resultant == 84
            and self.contact_jet_differences
            == (0, 0, Rational(165, 5488) * T - Rational(165, 10976))
            and self.contact_jet_separation == Rational(81675, 120472576)
            and self.collision_identity == 0
            and self.tangency_identity == 0
            and self.collision_tangency_gcd == expand((S - 1) ** 2)
            and self.collision_jets_at_contact == (0, 0, 0, 660, 10356)
            and self.tangency_jets_at_contact == (0, 0, -3960, -96456)
            and self.collision_factor_discriminants == (-2548, -156754683)
            and self.collision_factor_resultants == (7292322, -5, 44)
            and all(value != 0 for value in self.triple_boundary_resultants)
            and all(value != 0 for value in self.node_boundary_resultants)
            and self.cusp_image_factor == Rational(-299, 64)
            and self.extra_critical_factor == -894348
            and self.triple_x_eliminant_identity == 0
            and self.node_x_eliminant_identity == 0
            and self.node_x_discriminant == -4321892177659568
            and self.node_target_separations
            == (-60306, 792, Rational(-135043, 4096))
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_length == 17
            and self.sage_jacobian_radical_length == 7
            and self.sage_jacobian_components
            == ((4, 4), (4, 1), (5, 1), (4, 1))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.topology_propagation_dependencies
            == (
                "connected clean open",
                "proper projective Whitney-Thom triviality",
            )
            and not self.split_and_deeper_boundaries_classified_here
            and self.total_delta == self.arithmetic_genus == 28
            and bool(CONTACT_THREE_TRIPLE_RELATIONS)
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_contact_three_triple_certificate(
) -> DeltaTenContactThreeTripleCertificate:
    """Build the exact ``C3 + T111 + 4N`` certificate."""

    sample_coefficient = CONTACT_THREE_TRIPLE_COEFFICIENT_MATRIX.subs(
        SAMPLE_BASE_POINT,
    )
    sample_rhs = CONTACT_THREE_TRIPLE_RIGHT_HAND_SIDE.subs(SAMPLE_BASE_POINT)
    sample_vector = Matrix(
        [SAMPLE_PARAMETERS[parameter] for parameter in LINEAR_PARAMETERS],
    )
    sample_minors = tuple(
        expand(sample_coefficient.extract(rows, range(4)).det())
        for rows in combinations(range(5), 4)
    )

    full_sample = {**SAMPLE_BASE_POINT, **SAMPLE_PARAMETERS}
    base_variables = (KAPPA, TRIPLE_OMITTED_ROOT, CONTACT_THREE_SUM)
    base_jacobian = Matrix(
        [
            [diff(equation, variable) for variable in base_variables]
            for equation in CONTACT_THREE_TRIPLE_EQUATIONS
        ],
    ).subs(full_sample)
    parameter_jacobian = CONTACT_THREE_TRIPLE_COEFFICIENT_MATRIX.subs(
        full_sample,
    )
    tangent_residuals = tuple(
        tuple(
            expand(value)
            for value in (
                base_jacobian * Matrix(base_tangent)
                + parameter_jacobian * Matrix(coefficient_tangent)
            )
        )
        for base_tangent, coefficient_tangent in zip(
            SAMPLE_BASE_TANGENTS,
            SAMPLE_COEFFICIENT_TANGENTS,
            strict=True,
        )
    )
    image_tangents = Matrix(
        [
            (base_tangent[0], *coefficient_tangent)
            for base_tangent, coefficient_tangent in zip(
                SAMPLE_BASE_TANGENTS,
                SAMPLE_COEFFICIENT_TANGENTS,
                strict=True,
            )
        ],
    )

    sample_collision = expand(COLLISION_POLYNOMIAL.subs(SAMPLE_PARAMETERS))
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(SAMPLE_PARAMETERS))
    contact_jet_differences = _contact_jet_differences(3)
    sample_x_numerator = COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS)
    sample_x_denominator = COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    implicit_polynomial = Poly(CONTACT_THREE_TRIPLE_IMPLICIT, X, Y)

    if CONTACT_THREE_TRIPLE_RELATIONS:
        complement_census = n_generator_three_cycle_presentation_census(
            CONTACT_THREE_TRIPLE_RELATIONS,
            4,
        )
        cyclic_simplification = (1, 0, True)
    else:
        complement_census = ThreeCyclePresentationCensus(0, 0, ())
        cyclic_simplification = (0, 0, False)

    return DeltaTenContactThreeTripleCertificate(
        source_overlap_identity=CONTACT_TRIPLE_SOURCE_RESULTANT_IDENTITY,
        target_separation_identity=CONTACT_TRIPLE_TARGET_SEPARATION_IDENTITY,
        determinant_factor_data=SAGE_DETERMINANT_FACTOR_DATA,
        compatibility_residual_shape=SAGE_COMPATIBILITY_RESIDUAL_SHAPE,
        compatibility_rational_factor_count=(
            SAGE_COMPATIBILITY_RATIONAL_FACTOR_COUNT
        ),
        compatibility_sample_value=CONTACT_THREE_TRIPLE_COMPATIBILITY.subs(
            SAMPLE_BASE_POINT,
        ),
        compatibility_sample_gradient=tuple(
            diff(CONTACT_THREE_TRIPLE_COMPATIBILITY, variable).subs(
                SAMPLE_BASE_POINT,
            )
            for variable in base_variables
        ),
        sage_valid_singular_saturation=SAGE_VALID_SINGULAR_SATURATION,
        sage_coefficient_rank_three_curve=SAGE_COEFFICIENT_RANK_THREE_CURVE,
        sage_augmented_rank_three_unit=SAGE_AUGMENTED_RANK_THREE_UNIT,
        sage_rank_two_unit_exponents=SAGE_RANK_TWO_UNIT_EXPONENTS,
        sage_normalized_minor_counts=SAGE_NORMALIZED_MINOR_COUNTS,
        sample_coefficient_rank=sample_coefficient.rank(),
        sample_augmented_rank=sample_coefficient.row_join(sample_rhs).rank(),
        sample_incidence_residuals=tuple(
            expand(value) for value in sample_coefficient * sample_vector - sample_rhs
        ),
        sample_maximal_minors=sample_minors,
        sample_valid_localizer=CONTACT_THREE_TRIPLE_VALID_LOCALIZER.subs(
            SAMPLE_BASE_POINT,
        ),
        tangent_incidence_residuals=tangent_residuals,
        coefficient_image_tangent_rank=image_tangents.rank(),
        triple_cubic_discriminant=discriminant(SAMPLE_TRIPLE_CUBIC, T),
        triple_p_derivative_resultant=resultant(
            SAMPLE_TRIPLE_CUBIC,
            diff(SAMPLE_P, T),
            T,
        ),
        triple_omitted_root_separation=resultant(
            SAMPLE_TRIPLE_CUBIC,
            T - Rational(-1, 2),
            T,
        ),
        triple_q_remainder=rem(SAMPLE_Q, SAMPLE_TRIPLE_CUBIC, T),
        triple_omitted_q_difference=expand(
            SAMPLE_Q.subs(T, Rational(-1, 2)) - SAMPLE_TRIPLE_VALUE_Q,
        ),
        triple_quotient_identity=expand(
            SAMPLE_Q
            - SAMPLE_TRIPLE_VALUE_Q
            - SAMPLE_TRIPLE_CUBIC * SAMPLE_TRIPLE_QUOTIENT,
        ),
        triple_slope_eliminant_identity=expand(
            262144
            * resultant(
                SAMPLE_TRIPLE_CUBIC,
                SAMPLE_TRIPLE_QUOTIENT
                - TRIPLE_SLOPE * (T + Rational(1, 2)),
                T,
            )
            + SAMPLE_TRIPLE_SLOPE_POLYNOMIAL
        ),
        triple_slope_discriminant=discriminant(
            SAMPLE_TRIPLE_SLOPE_POLYNOMIAL,
            TRIPLE_SLOPE,
        ),
        contact_pair_discriminant=discriminant(SAMPLE_CONTACT_BRANCH, T),
        contact_image_remainders=(
            rem(SAMPLE_P, SAMPLE_CONTACT_BRANCH, T),
            rem(SAMPLE_Q, SAMPLE_CONTACT_BRANCH, T),
        ),
        contact_p_derivative_resultant=resultant(
            SAMPLE_CONTACT_BRANCH,
            diff(SAMPLE_P, T),
            T,
        ),
        contact_jet_differences=contact_jet_differences,
        contact_jet_separation=resultant(
            SAMPLE_CONTACT_BRANCH,
            contact_jet_differences[-1],
            T,
        ),
        collision_identity=expand(sample_collision - SAMPLE_COLLISION),
        tangency_identity=expand(sample_tangency - SAMPLE_TANGENCY),
        collision_tangency_gcd=gcd(
            Poly(sample_collision, S),
            Poly(sample_tangency, S),
        ).as_expr(),
        collision_jets_at_contact=tuple(
            diff(sample_collision, S, order).subs(S, 1) for order in range(5)
        ),
        tangency_jets_at_contact=tuple(
            diff(sample_tangency, S, order).subs(S, 1) for order in range(4)
        ),
        collision_factor_discriminants=(
            discriminant(SAMPLE_TRIPLE_PAIR_SUMS, S),
            discriminant(SAMPLE_RESIDUAL_NODE_SUMS, S),
        ),
        collision_factor_resultants=(
            resultant(SAMPLE_TRIPLE_PAIR_SUMS, SAMPLE_RESIDUAL_NODE_SUMS, S),
            resultant(SAMPLE_TRIPLE_PAIR_SUMS, S - 1, S),
            resultant(SAMPLE_RESIDUAL_NODE_SUMS, S - 1, S),
        ),
        triple_boundary_resultants=tuple(
            resultant(SAMPLE_TRIPLE_PAIR_SUMS, expression.subs(SAMPLE_PARAMETERS), S)
            for expression in (
                PAIR_DENOMINATOR,
                PAIR_QUADRATIC,
                PAIR_DIAGONAL_FACTOR,
                TANGENCY_POLYNOMIAL,
            )
        ),
        node_boundary_resultants=tuple(
            resultant(SAMPLE_RESIDUAL_NODE_SUMS, expression.subs(SAMPLE_PARAMETERS), S)
            for expression in (
                PAIR_DENOMINATOR,
                PAIR_QUADRATIC,
                PAIR_DIAGONAL_FACTOR,
                TANGENCY_POLYNOMIAL,
            )
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(SAMPLE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(SAMPLE_PARAMETERS),
        triple_x_eliminant_identity=expand(
            resultant(
                SAMPLE_TRIPLE_PAIR_SUMS,
                X * sample_x_denominator - sample_x_numerator,
                S,
            )
            - 36 * (16 * X - 13) ** 3
        ),
        node_x_eliminant_identity=expand(
            resultant(
                SAMPLE_RESIDUAL_NODE_SUMS,
                X * sample_x_denominator - sample_x_numerator,
                S,
            )
            - 46656 * CONTACT_THREE_TRIPLE_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(
            CONTACT_THREE_TRIPLE_NODE_X_POLYNOMIAL,
            X,
        ),
        node_target_separations=tuple(
            CONTACT_THREE_TRIPLE_NODE_X_POLYNOMIAL.subs(X, target)
            for target in (0, 3, Rational(13, 16))
        ),
        implicit_resultant_identity=expand(
            resultant(X - SAMPLE_P, Y - SAMPLE_Y, T)
            + CONTACT_THREE_TRIPLE_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            CONTACT_THREE_TRIPLE_IMPLICIT.subs({X: SAMPLE_P, Y: SAMPLE_Y}),
        ),
        implicit_content=implicit_polynomial.content(),
        sage_jacobian_length=17,
        sage_jacobian_radical_length=7,
        sage_jacobian_components=((4, 4), (4, 1), (5, 1), (4, 1)),
        sage_cyclic_simplification=cyclic_simplification,
        topology_propagation_dependencies=(
            "connected clean open",
            "proper projective Whitney-Thom triviality",
        ),
        split_and_deeper_boundaries_classified_here=False,
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=2,
        contact_three_delta=3,
        ordinary_triple_delta=3,
        node_count=4,
        infinity_delta=16,
        complement_census=complement_census,
    )


def main() -> int:
    """Print the exact profile certificate and fail on any regression."""

    certificate = exact_delta_ten_contact_three_triple_certificate()
    print("delta-ten C3 + T111 + 4N certificate:", certificate.verified)
    print("compatibility residual shape:", certificate.compatibility_residual_shape)
    print(
        "geometrically irreducible valid surface:",
        certificate.compatibility_surface_geometrically_irreducible,
    )
    print("valid rank-drop incidence empty:", certificate.valid_rank_drop_incidence_empty)
    print("coefficient-image dimension:", certificate.dominant_incidence_dimension)
    print("sample delta balance:", certificate.total_delta, "= 2 + 3 + 3 + 4 + 16")
    print("sample complement census:", certificate.complement_census)
    print("split and deeper boundaries remain open")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
