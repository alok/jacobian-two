"""Certify the conditional delta-ten ``C2 + C3 + 5N`` profile.

For the normalized degree-``(4, 9)`` family

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``,

choose two distinct valid unordered-pair sums ``u`` and ``v``.  A contact of
order two at ``u`` and one of order three at ``v`` impose

``H(u) = H'(u) = H(v) = H'(v) = H''(v) = 0``.

These five equations are affine-linear in ``(a,b,c,d)``.  The augmented
determinant factors as the expected collision-boundary factors times one
409-term compatibility polynomial ``R(k,u,v)``.  The independent Sage
checker derives ``R`` rather than trusting a stored expansion and proves it
is irreducible over ``QQ``.  Its exact localized determinantal audit finds a
genuine coefficient-rank-three base curve of degree thirty, but proves that
the augmented-rank-three ideal is the unit ideal.  Thus the entire rank-drop
curve is inconsistent and cannot support a hidden coefficient fiber.  The
full-localizer singular-Jacobian ideal is also the unit ideal, with exact
saturation exponent two, so the compatibility hypersurface is smooth on this
valid chart.  The
rational point ``(k,u,v)=(2/3,1,-1)`` is smooth, has coefficient rank four,
and lies off every declared valid-chart boundary.  Irreducibility over
``QQ`` plus a
smooth rational point makes ``R`` geometrically irreducible: if its geometric
components were distinct Galois conjugates, that rational point would lie on
all of them and would be singular.

The corresponding coefficient solution ``(a,b,c,d)=(3,1,3,0)`` gives one
exact contact two, one exact contact three, five nodes, the forced cusp, and
the fixed infinity branch.  Sage 10.8 regenerates the implicit equation and
singular scheme.  Its stored exact affine van Kamp presentation simplifies
to ``Z``; an exhaustive dependency-free replay of all ``40^4`` three-cycle
assignments has no ``A6`` image.

This is a conditional, computer-assisted profile certificate.  The split,
denominator, cusp-pair, diagonal, pair-overlap, same-target,
singular-compatibility, and deeper boundary charts named by the certificate
are not silently removed from the larger audit.  Nothing here proves the
plane Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
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
    rem,
    resultant,
    together,
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

CONTACT_TWO_SUM: Final = Symbol("u")
CONTACT_THREE_SUM: Final = Symbol("v")
CONTACT_TWO_THREE_LINEAR_PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)
CONTACT_TWO_THREE_ZERO_PARAMETERS: Final = dict.fromkeys(
    CONTACT_TWO_THREE_LINEAR_PARAMETERS,
    0,
)


def _jet_equations(root: Expr, multiplicity: int) -> tuple[Expr, ...]:
    """Return the collision jets imposing one contact at ``root``."""

    return tuple(
        diff(COLLISION_POLYNOMIAL, S, order).subs(S, root)
        for order in range(multiplicity)
    )


CONTACT_TWO_THREE_EQUATIONS: Final = (
    *_jet_equations(CONTACT_TWO_SUM, 2),
    *_jet_equations(CONTACT_THREE_SUM, 3),
)
CONTACT_TWO_THREE_COEFFICIENT_MATRIX: Final = Matrix(
    [
        [diff(equation, parameter) for parameter in CONTACT_TWO_THREE_LINEAR_PARAMETERS]
        for equation in CONTACT_TWO_THREE_EQUATIONS
    ]
)
CONTACT_TWO_THREE_CONSTANT_COLUMN: Final = Matrix(
    [
        equation.subs(CONTACT_TWO_THREE_ZERO_PARAMETERS)
        for equation in CONTACT_TWO_THREE_EQUATIONS
    ]
)
CONTACT_TWO_THREE_AUGMENTED_MATRIX: Final = (
    CONTACT_TWO_THREE_COEFFICIENT_MATRIX.row_join(
        CONTACT_TWO_THREE_CONSTANT_COLUMN,
    )
)


def _at_sum(expression: Expr, root: Expr) -> Expr:
    """Substitute one unordered-pair sum into a generic pair expression."""

    return expand(expression.subs(S, root))


def _cleared_pair_quadratic(root: Expr) -> Expr:
    """Return the pair-source quadratic with its denominator cleared."""

    denominator = _at_sum(PAIR_DENOMINATOR, root)
    quadratic = _at_sum(PAIR_QUADRATIC, root)
    return expand(
        denominator * T**2
        - root * denominator * T
        + root * quadratic
    )


CONTACT_TWO_CLEARED_PAIR: Final = _cleared_pair_quadratic(CONTACT_TWO_SUM)
CONTACT_THREE_CLEARED_PAIR: Final = _cleared_pair_quadratic(CONTACT_THREE_SUM)
CONTACT_TWO_THREE_PAIR_DISJOINT_FACTOR: Final = cancel(
    resultant(
        CONTACT_TWO_CLEARED_PAIR,
        CONTACT_THREE_CLEARED_PAIR,
        T,
    )
    / (CONTACT_TWO_SUM - CONTACT_THREE_SUM) ** 2
)
CONTACT_TWO_THREE_TARGET_SEPARATION_FACTOR: Final = (
    KAPPA + CONTACT_TWO_SUM + CONTACT_THREE_SUM
)
CONTACT_TWO_THREE_TARGET_SEPARATION_IDENTITY: Final = expand(
    _at_sum(COLLISION_X_NUMERATOR, CONTACT_TWO_SUM)
    * _at_sum(COLLISION_X_DENOMINATOR, CONTACT_THREE_SUM)
    - _at_sum(COLLISION_X_NUMERATOR, CONTACT_THREE_SUM)
    * _at_sum(COLLISION_X_DENOMINATOR, CONTACT_TWO_SUM)
    + (CONTACT_TWO_SUM - CONTACT_THREE_SUM)
    * CONTACT_TWO_THREE_TARGET_SEPARATION_FACTOR
    * CONTACT_TWO_THREE_PAIR_DISJOINT_FACTOR
)

CONTACT_TWO_THREE_VALID_LOCALIZER: Final = expand(
    KAPPA
    * (KAPPA**2 - 4)
    * CONTACT_TWO_SUM
    * CONTACT_THREE_SUM
    * (CONTACT_TWO_SUM - CONTACT_THREE_SUM)
    * CONTACT_TWO_THREE_TARGET_SEPARATION_FACTOR
    * CONTACT_TWO_THREE_PAIR_DISJOINT_FACTOR
    * _at_sum(PAIR_DENOMINATOR, CONTACT_TWO_SUM)
    * _at_sum(PAIR_DENOMINATOR, CONTACT_THREE_SUM)
    * _at_sum(PAIR_QUADRATIC, CONTACT_TWO_SUM)
    * _at_sum(PAIR_QUADRATIC, CONTACT_THREE_SUM)
    * _at_sum(PAIR_DIAGONAL_FACTOR, CONTACT_TWO_SUM)
    * _at_sum(PAIR_DIAGONAL_FACTOR, CONTACT_THREE_SUM)
)

# Reproducible Sage 10.8 metadata from
# ``tools/check_a6_delta_ten_contact_two_three.sage``.  The 409-term residual
# is intentionally derived in that checker, not duplicated here as an opaque
# literal.
SAGE_COMPATIBILITY_FACTOR: Final = (
    -2,
    6,
    1,
    3,
)
SAGE_COMPATIBILITY_RESIDUAL_SHAPE: Final = (409, 18, 13, 9, 10)
SAGE_COMPATIBILITY_RATIONAL_FACTOR_COUNT: Final = 1
SAGE_COMPATIBILITY_SAMPLE_GRADIENT: Final = (
    Rational(1073741824, 59049),
    Rational(536870912, 59049),
    Rational(1073741824, 177147),
)
SAGE_COEFFICIENT_RANK_THREE_CURVE: Final = (1, 30, -51, 4, 0)
SAGE_AUGMENTED_RANK_THREE_UNIT: Final = (True, 3)
SAGE_RANK_TWO_UNIT_EXPONENTS: Final = (3, 3)
SAGE_VALID_SINGULAR_SATURATION: Final = (True, 2)

CONTACT_TWO_THREE_PARAMETERS: Final = {
    KAPPA: Rational(2, 3),
    ALPHA: 3,
    BETA: 1,
    GAMMA: 3,
    DELTA: 0,
}
CONTACT_TWO_THREE_BASE_POINT: Final = {
    KAPPA: Rational(2, 3),
    CONTACT_TWO_SUM: 1,
    CONTACT_THREE_SUM: -1,
}
CONTACT_TWO_THREE_P: Final = expand(FAMILY_P.subs(CONTACT_TWO_THREE_PARAMETERS))
CONTACT_TWO_THREE_Q: Final = expand(FAMILY_Q.subs(CONTACT_TWO_THREE_PARAMETERS))
CONTACT_TWO_THREE_SCALED_P: Final = expand(3 * CONTACT_TWO_THREE_P)

CONTACT_TWO_BRANCH: Final = T**2 - T + 1
CONTACT_THREE_BRANCH: Final = T**2 + T + 1
CONTACT_TWO_THREE_RESIDUAL: Final = (
    9 * S**5
    + 27 * S**4
    + 45 * S**3
    + 53 * S**2
    + 46 * S
    + 12
)
CONTACT_TWO_THREE_COLLISION: Final = expand(
    (S - 1) ** 2 * (S + 1) ** 3 * CONTACT_TWO_THREE_RESIDUAL / 9
)
CONTACT_TWO_THREE_TANGENCY_COFACTOR: Final = (
    81 * S**8
    + 270 * S**7
    + 477 * S**6
    + 516 * S**5
    + 313 * S**4
    + 24 * S**3
    - 59 * S**2
    - 66 * S
    - 20
)
CONTACT_TWO_THREE_TANGENCY: Final = expand(
    Rational(4, 9)
    * (S - 1)
    * (S + 1) ** 2
    * CONTACT_TWO_THREE_TANGENCY_COFACTOR
)

CONTACT_TWO_THREE_NODE_X_POLYNOMIAL: Final = (
    -729 * X**5
    + 15309 * X**4
    - 16011 * X**3
    + 93727 * X**2
    - 14512 * X
    + 50112
)
CONTACT_TWO_THREE_IMPLICIT: Final = (
    -X**9
    - 6 * X**8
    - 309 * X**7
    - 354 * X**6 * Y
    + 19 * X**6
    + 1986 * X**5 * Y
    - 216 * X**5
    + 1215 * X**4 * Y**2
    - 552 * X**4 * Y
    + 7071 * X**3 * Y**2
    + 1728 * X**3 * Y
    + 13122 * X**2 * Y**3
    + 10497 * X**2 * Y**2
    + 486 * X * Y**3
    + 1647 * X * Y**2
    + 19683 * Y**4
    + 21032 * Y**3
    + 5832 * Y**2
)

# Sage 10.8 exact raw affine presentation.  Generator indices 1..4 are the
# geometric meridians of a generic vertical fiber.
CONTACT_TWO_THREE_RELATIONS: Final = (
    (2, 1, 2, 1, -2, -1, -2, -1),
    (4, 3, -4, -3),
    (-2, -1, -3, 1, 2, 1, -2, -1, 3, 1, 2, -1),
    (-2, -1, 4, 1),
    (4, 3, 4, 3, 4, 3, -4, -3, -4, -3, -4, -3),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (
        -2,
        -1,
        -2,
        -1,
        3,
        4,
        3,
        4,
        -3,
        -4,
        -3,
        1,
        2,
        1,
        2,
        -1,
        -2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        -4,
        -3,
        1,
        2,
        1,
    ),
    (-2, -1, 3, 4, -3, -4, -3, 1, 2, 1),
    (4, 3, -4, -3),
    (-4, -3, 1, 3, 4, -3, -1, 3),
)


def _modulo_branch(rational_expression: Expr, branch: Expr) -> Expr:
    """Reduce one rational jet modulo a contact-pair quadratic."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    denominator_inverse = invert(denominator, branch)
    return expand(rem(numerator * denominator_inverse, branch, T))


def _branch_jet_differences(
    branch: Expr,
    pair_sum: Expr,
    order: int,
) -> tuple[Expr, ...]:
    """Return successive ``dQ/dX`` differences between paired branches."""

    current = CONTACT_TWO_THREE_Q
    p_derivative = diff(CONTACT_TWO_THREE_SCALED_P, T)
    differences: list[Expr] = []
    for _ in range(order):
        current = cancel(diff(current, T) / p_derivative)
        reduced = _modulo_branch(current, branch)
        partner = _modulo_branch(current.subs(T, pair_sum - T), branch)
        differences.append(expand(rem(reduced - partner, branch, T)))
    return tuple(differences)


@dataclass(frozen=True, slots=True)
class DeltaTenContactTwoThreeCertificate:
    """Exact incidence, sample geometry, and complement data."""

    target_separation_identity: Expr
    compatibility_factor_exponents: tuple[int, int, int, int]
    compatibility_residual_shape: tuple[int, int, int, int, int]
    compatibility_rational_factor_count: int
    compatibility_sample_value: Expr
    compatibility_sample_gradient: tuple[Expr, Expr, Expr]
    sage_coefficient_rank_three_curve: tuple[int, int, int, int, int]
    sage_augmented_rank_three_unit: tuple[bool, int]
    sage_rank_two_unit_exponents: tuple[int, int]
    sage_valid_singular_saturation: tuple[bool, int]
    sample_coefficient_rank: int
    sample_augmented_rank: int
    sample_incidence_residuals: tuple[Expr, ...]
    sample_maximal_minors: tuple[Expr, ...]
    sample_valid_localizer: Expr
    collision_identity: Expr
    tangency_identity: Expr
    collision_tangency_gcd: Expr
    double_root_support: Expr
    triple_root_support: Expr
    collision_jets_at_contact_two: tuple[Expr, ...]
    collision_jets_at_contact_three: tuple[Expr, ...]
    tangency_jets_at_contact_two: tuple[Expr, ...]
    tangency_jets_at_contact_three: tuple[Expr, ...]
    residual_discriminant: Expr
    residual_contact_separations: tuple[Expr, Expr]
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    residual_boundary_resultants: tuple[Expr, Expr, Expr, Expr]
    contact_chart_values: tuple[tuple[Expr, Expr, Expr], ...]
    pair_disjoint_resultant: Expr
    pair_target_separation_numerator: Expr
    contact_pair_discriminants: tuple[Expr, Expr]
    contact_image_remainders: tuple[tuple[Expr, Expr], ...]
    contact_p_derivative_resultants: tuple[Expr, Expr]
    contact_jet_differences: tuple[tuple[Expr, ...], ...]
    contact_jet_separations: tuple[Expr, Expr]
    node_x_resultant_identity: Expr
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
    contact_two_delta: int
    contact_three_delta: int
    node_count: int
    infinity_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def compatibility_sample_is_smooth(self) -> bool:
        """Whether the rational compatibility point has nonzero gradient."""

        return self.compatibility_sample_value == 0 and any(
            value != 0 for value in self.compatibility_sample_gradient
        )

    @property
    def compatibility_surface_geometrically_irreducible(self) -> bool:
        """Apply the smooth-rational-point geometric irreducibility test."""

        return (
            self.compatibility_rational_factor_count == 1
            and self.compatibility_sample_is_smooth
        )

    @property
    def dominant_incidence_dimension(self) -> int:
        """Return base dimension three minus one compatibility equation."""

        return 2

    @property
    def valid_rank_drop_incidence_empty(self) -> bool:
        """Whether every valid coefficient-rank-drop base is inconsistent."""

        return bool(
            self.sage_coefficient_rank_three_curve == (1, 30, -51, 4, 0)
            and self.sage_augmented_rank_three_unit == (True, 3)
            and self.sage_rank_two_unit_exponents == (3, 3)
        )

    @property
    def valid_compatibility_surface_smooth(self) -> bool:
        """Whether the residual hypersurface is smooth on the valid chart."""

        return self.sage_valid_singular_saturation == (True, 2)

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return (
            self.cusp_delta
            + self.contact_two_delta
            + self.contact_three_delta
            + self.node_count
            + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether every exact incidence and sample invariant agrees."""

        return bool(
            self.target_separation_identity == 0
            and self.compatibility_factor_exponents == (-2, 6, 1, 3)
            and self.compatibility_residual_shape == (409, 18, 13, 9, 10)
            and self.compatibility_rational_factor_count == 1
            and self.compatibility_sample_gradient
            == SAGE_COMPATIBILITY_SAMPLE_GRADIENT
            and self.compatibility_surface_geometrically_irreducible
            and self.sage_coefficient_rank_three_curve == (1, 30, -51, 4, 0)
            and self.sage_augmented_rank_three_unit == (True, 3)
            and self.sage_rank_two_unit_exponents == (3, 3)
            and self.valid_rank_drop_incidence_empty
            and self.sage_valid_singular_saturation == (True, 2)
            and self.valid_compatibility_surface_smooth
            and self.sample_coefficient_rank == 4
            and self.sample_augmented_rank == 4
            and self.sample_incidence_residuals == (0, 0, 0, 0, 0)
            and all(value != 0 for value in self.sample_maximal_minors)
            and self.sample_valid_localizer != 0
            and self.dominant_incidence_dimension == 2
            and self.collision_identity == 0
            and self.tangency_identity == 0
            and self.collision_tangency_gcd == expand((S - 1) * (S + 1) ** 2)
            and self.double_root_support == expand((S - 1) * (S + 1) ** 2)
            and self.triple_root_support == S + 1
            and self.collision_jets_at_contact_two
            == (0, 0, Rational(1024, 3), Rational(11648, 3))
            and self.collision_jets_at_contact_three
            == (0, 0, 0, Rational(-64, 3), Rational(640, 3))
            and self.tangency_jets_at_contact_two
            == (0, Rational(8192, 3), 37888)
            and self.tangency_jets_at_contact_three
            == (0, 0, Rational(-256, 3), 1024)
            and self.residual_discriminant == 256979755008
            and self.residual_contact_separations == (192, -8)
            and self.cusp_image_factor == Rational(8, 27)
            and self.extra_critical_factor == 9072
            and self.residual_boundary_resultants
            == (
                Rational(-1024, 27),
                Rational(2048, 9),
                28672,
                -42392524029926834176,
            )
            and self.contact_chart_values
            == (
                (Rational(8, 3), Rational(8, 3), 8),
                (Rational(-4, 3), Rational(4, 3), 4),
            )
            and self.pair_disjoint_resultant == Rational(4096, 81)
            and self.pair_target_separation_numerator == Rational(-4096, 243)
            and self.contact_pair_discriminants == (-3, -3)
            and self.contact_image_remainders == ((-5, 3), (-1, -1))
            and self.contact_p_derivative_resultants == (252, 36)
            and self.contact_jet_differences
            == (
                (0, Rational(3, 49) * (2 * T - 1)),
                (0, 0, -T / 2 - Rational(1, 4)),
            )
            and self.contact_jet_separations
            == (Rational(27, 2401), Rational(3, 16))
            and self.node_x_resultant_identity == 0
            and self.node_x_discriminant
            == 74803083224844271333355450007552000000
            and self.node_target_separations == (16313472, 190400, 50112)
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_length == 17
            and self.sage_jacobian_radical_length == 8
            and self.sage_jacobian_components
            == ((5, 5), (5, 1), (3, 1), (4, 1))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.topology_propagation_dependencies
            == (
                "connected clean open",
                "proper projective Whitney-Thom triviality",
            )
            and not self.split_and_deeper_boundaries_classified_here
            and self.total_delta == self.arithmetic_genus == 28
            and len(CONTACT_TWO_THREE_RELATIONS) == 10
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_contact_two_three_certificate(
) -> DeltaTenContactTwoThreeCertificate:
    """Build the exact dominant-incidence and rational-member certificate."""

    sample_matrix = CONTACT_TWO_THREE_COEFFICIENT_MATRIX.subs(
        CONTACT_TWO_THREE_BASE_POINT,
    )
    sample_constant = CONTACT_TWO_THREE_CONSTANT_COLUMN.subs(
        CONTACT_TWO_THREE_BASE_POINT,
    )
    sample_parameter_vector = Matrix(
        [
            CONTACT_TWO_THREE_PARAMETERS[parameter]
            for parameter in CONTACT_TWO_THREE_LINEAR_PARAMETERS
        ]
    )
    sample_maximal_minors = tuple(
        expand(
            sample_matrix.extract(
                [row for row in range(5) if row != omitted_row],
                range(4),
            ).det()
        )
        for omitted_row in range(5)
    )

    sample_collision = expand(
        COLLISION_POLYNOMIAL.subs(CONTACT_TWO_THREE_PARAMETERS),
    )
    sample_tangency = expand(
        TANGENCY_POLYNOMIAL.subs(CONTACT_TWO_THREE_PARAMETERS),
    )
    derivatives = tuple(
        Poly(diff(sample_collision, S, order), S) for order in range(3)
    )
    double_root_support = gcd(derivatives[0], derivatives[1]).as_expr()
    triple_root_support = gcd(
        gcd(derivatives[0], derivatives[1]),
        derivatives[2],
    ).as_expr()

    contact_branches = (CONTACT_TWO_BRANCH, CONTACT_THREE_BRANCH)
    contact_sums = (1, -1)
    contact_orders = (2, 3)
    jet_differences = tuple(
        _branch_jet_differences(branch, pair_sum, order)
        for branch, pair_sum, order in zip(
            contact_branches,
            contact_sums,
            contact_orders,
            strict=True,
        )
    )

    sample_x_numerator = COLLISION_X_NUMERATOR.subs(
        CONTACT_TWO_THREE_PARAMETERS,
    )
    sample_x_denominator = COLLISION_X_DENOMINATOR.subs(
        CONTACT_TWO_THREE_PARAMETERS,
    )
    implicit_polynomial = Poly(CONTACT_TWO_THREE_IMPLICIT, X, Y)

    return DeltaTenContactTwoThreeCertificate(
        target_separation_identity=CONTACT_TWO_THREE_TARGET_SEPARATION_IDENTITY,
        compatibility_factor_exponents=SAGE_COMPATIBILITY_FACTOR,
        compatibility_residual_shape=SAGE_COMPATIBILITY_RESIDUAL_SHAPE,
        compatibility_rational_factor_count=(
            SAGE_COMPATIBILITY_RATIONAL_FACTOR_COUNT
        ),
        compatibility_sample_value=0,
        compatibility_sample_gradient=SAGE_COMPATIBILITY_SAMPLE_GRADIENT,
        sage_coefficient_rank_three_curve=SAGE_COEFFICIENT_RANK_THREE_CURVE,
        sage_augmented_rank_three_unit=SAGE_AUGMENTED_RANK_THREE_UNIT,
        sage_rank_two_unit_exponents=SAGE_RANK_TWO_UNIT_EXPONENTS,
        sage_valid_singular_saturation=SAGE_VALID_SINGULAR_SATURATION,
        sample_coefficient_rank=sample_matrix.rank(),
        sample_augmented_rank=sample_matrix.row_join(sample_constant).rank(),
        sample_incidence_residuals=tuple(
            expand(value)
            for value in sample_matrix * sample_parameter_vector + sample_constant
        ),
        sample_maximal_minors=sample_maximal_minors,
        sample_valid_localizer=CONTACT_TWO_THREE_VALID_LOCALIZER.subs(
            CONTACT_TWO_THREE_BASE_POINT,
        ),
        collision_identity=expand(
            sample_collision - CONTACT_TWO_THREE_COLLISION,
        ),
        tangency_identity=expand(
            sample_tangency - CONTACT_TWO_THREE_TANGENCY,
        ),
        collision_tangency_gcd=gcd(
            Poly(sample_collision, S),
            Poly(sample_tangency, S),
        ).as_expr(),
        double_root_support=double_root_support,
        triple_root_support=triple_root_support,
        collision_jets_at_contact_two=tuple(
            diff(sample_collision, S, order).subs(S, 1) for order in range(4)
        ),
        collision_jets_at_contact_three=tuple(
            diff(sample_collision, S, order).subs(S, -1) for order in range(5)
        ),
        tangency_jets_at_contact_two=tuple(
            diff(sample_tangency, S, order).subs(S, 1) for order in range(3)
        ),
        tangency_jets_at_contact_three=tuple(
            diff(sample_tangency, S, order).subs(S, -1) for order in range(4)
        ),
        residual_discriminant=discriminant(CONTACT_TWO_THREE_RESIDUAL, S),
        residual_contact_separations=(
            CONTACT_TWO_THREE_RESIDUAL.subs(S, 1),
            CONTACT_TWO_THREE_RESIDUAL.subs(S, -1),
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(CONTACT_TWO_THREE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(
            CONTACT_TWO_THREE_PARAMETERS,
        ),
        residual_boundary_resultants=(
            resultant(
                CONTACT_TWO_THREE_RESIDUAL,
                PAIR_DENOMINATOR.subs(CONTACT_TWO_THREE_PARAMETERS),
                S,
            ),
            resultant(
                CONTACT_TWO_THREE_RESIDUAL,
                PAIR_QUADRATIC.subs(CONTACT_TWO_THREE_PARAMETERS),
                S,
            ),
            resultant(
                CONTACT_TWO_THREE_RESIDUAL,
                PAIR_DIAGONAL_FACTOR.subs(CONTACT_TWO_THREE_PARAMETERS),
                S,
            ),
            resultant(CONTACT_TWO_THREE_RESIDUAL, sample_tangency, S),
        ),
        contact_chart_values=tuple(
            tuple(
                expression.subs(CONTACT_TWO_THREE_PARAMETERS).subs(S, pair_sum)
                for expression in (
                    PAIR_DENOMINATOR,
                    PAIR_QUADRATIC,
                    PAIR_DIAGONAL_FACTOR,
                )
            )
            for pair_sum in contact_sums
        ),
        pair_disjoint_resultant=resultant(
            CONTACT_TWO_CLEARED_PAIR,
            CONTACT_THREE_CLEARED_PAIR,
            T,
        ).subs(CONTACT_TWO_THREE_BASE_POINT),
        pair_target_separation_numerator=(
            _at_sum(COLLISION_X_NUMERATOR, CONTACT_TWO_SUM)
            * _at_sum(COLLISION_X_DENOMINATOR, CONTACT_THREE_SUM)
            - _at_sum(COLLISION_X_NUMERATOR, CONTACT_THREE_SUM)
            * _at_sum(COLLISION_X_DENOMINATOR, CONTACT_TWO_SUM)
        ).subs(CONTACT_TWO_THREE_BASE_POINT),
        contact_pair_discriminants=tuple(
            discriminant(branch, T) for branch in contact_branches
        ),
        contact_image_remainders=tuple(
            (
                rem(CONTACT_TWO_THREE_SCALED_P, branch, T),
                rem(CONTACT_TWO_THREE_Q, branch, T),
            )
            for branch in contact_branches
        ),
        contact_p_derivative_resultants=tuple(
            resultant(
                branch,
                diff(CONTACT_TWO_THREE_SCALED_P, T),
                T,
            )
            for branch in contact_branches
        ),
        contact_jet_differences=jet_differences,
        contact_jet_separations=(
            resultant(CONTACT_TWO_BRANCH, jet_differences[0][-1], T),
            resultant(CONTACT_THREE_BRANCH, jet_differences[1][-1], T),
        ),
        node_x_resultant_identity=expand(
            81
            * resultant(
                CONTACT_TWO_THREE_RESIDUAL,
                3 * sample_x_numerator - X * sample_x_denominator,
                S,
            )
            - 1048576 * CONTACT_TWO_THREE_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(
            CONTACT_TWO_THREE_NODE_X_POLYNOMIAL,
            X,
        ),
        node_target_separations=tuple(
            CONTACT_TWO_THREE_NODE_X_POLYNOMIAL.subs(X, target)
            for target in (-5, -1, 0)
        ),
        implicit_resultant_identity=expand(
            resultant(
                CONTACT_TWO_THREE_SCALED_P - X,
                CONTACT_TWO_THREE_Q - Y,
                T,
            )
            - CONTACT_TWO_THREE_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            CONTACT_TWO_THREE_IMPLICIT.subs(
                {X: CONTACT_TWO_THREE_SCALED_P, Y: CONTACT_TWO_THREE_Q},
            )
        ),
        implicit_content=implicit_polynomial.content(),
        sage_jacobian_length=17,
        sage_jacobian_radical_length=8,
        sage_jacobian_components=((5, 5), (5, 1), (3, 1), (4, 1)),
        sage_cyclic_simplification=(1, 0, True),
        topology_propagation_dependencies=(
            "connected clean open",
            "proper projective Whitney-Thom triviality",
        ),
        split_and_deeper_boundaries_classified_here=False,
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=2,
        contact_two_delta=2,
        contact_three_delta=3,
        node_count=5,
        infinity_delta=16,
        complement_census=n_generator_three_cycle_presentation_census(
            CONTACT_TWO_THREE_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact ``C2 + C3 + 5N`` certificate and fail on regression."""

    certificate = exact_delta_ten_contact_two_three_certificate()
    print("delta-ten C2 + C3 + 5N certificate:", certificate.verified)
    print(
        "compatibility residual shape:",
        certificate.compatibility_residual_shape,
    )
    print(
        "geometrically irreducible compatibility surface:",
        certificate.compatibility_surface_geometrically_irreducible,
    )
    print(
        "valid rank-drop incidence empty:",
        certificate.valid_rank_drop_incidence_empty,
    )
    print(
        "valid compatibility surface smooth:",
        certificate.valid_compatibility_surface_smooth,
    )
    print(
        "sample singularity delta balance:",
        certificate.total_delta,
        "= 2 + 2 + 3 + 5 + 16",
    )
    print("sample complement census:", certificate.complement_census)
    print("split and deeper boundaries remain open")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
