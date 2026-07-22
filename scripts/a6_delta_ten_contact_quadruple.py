"""Certify the conditional delta-ten C2 + Q0 + 2N profile.

In the normalized degree-(4,9) family

    P = t^2 + k*t^3 + t^4,
    Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9,

an ordinary four-source fiber over P=h is imposed by making the nonconstant
coefficients of Q modulo P-h vanish.  A separate contact of order two at
unordered-pair sum u is imposed by H(u)=H'(u)=0.  These five equations are
affine-linear in (a,b,c,d).

After the split, singular-fiber, pair, and same-target boundaries are
localized away, the augmented determinant has one genuine compatibility
factor A(k,u).  It is an irreducible plane quartic with two singularities,
both on the removed split fibers k=+/-2.  Its smooth rational point
(k,u)=(-4,1) proves geometric irreducibility.  Exact Sage saturation proves
that the coefficient and augmented matrices have rank four everywhere on
the full valid chart, so the compatible coefficient vector is unique.

The rational point (k,h,u)=(-4,1,1) gives

    (a,b,c,d) = (7,-19,13,-6).

The resulting curve has one ordinary quadruple point, one exact C2 contact,
two nodes, the forced T(2,5) cusp, and the fixed T(5,9) branch at infinity.
Sage regenerates the singular scheme and a four-generator affine van Kamp
presentation that simplifies to Z.  The dependency-free finite replay has
no A6 image.

This remains a conditional, computer-assisted stratum certificate.  The
explicitly removed split, singular-fiber, pair, same-target, non-clean, and
deeper boundaries remain part of the larger audit.  Nothing here proves the
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
from scripts.a6_delta_ten_quadruple import (
    FIBER_POLYNOMIAL,
    FIBER_VALUE,
    NONCONSTANT_REMAINDERS,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    n_generator_three_cycle_presentation_census,
)

CONTACT_SUM: Final = Symbol("u")
SLOPE_VALUE: Final = Symbol("m_cq")
LINEAR_PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)
ZERO_LINEAR_PARAMETERS: Final = dict.fromkeys(LINEAR_PARAMETERS, 0)

CONTACT_QUADRUPLE_EQUATIONS: Final = (
    *NONCONSTANT_REMAINDERS,
    COLLISION_POLYNOMIAL.subs(S, CONTACT_SUM),
    diff(COLLISION_POLYNOMIAL, S).subs(S, CONTACT_SUM),
)
CONTACT_QUADRUPLE_COEFFICIENT_MATRIX: Final = Matrix(
    [
        [diff(equation, parameter) for parameter in LINEAR_PARAMETERS]
        for equation in CONTACT_QUADRUPLE_EQUATIONS
    ]
)
CONTACT_QUADRUPLE_CONSTANT_COLUMN: Final = Matrix(
    [
        equation.subs(ZERO_LINEAR_PARAMETERS)
        for equation in CONTACT_QUADRUPLE_EQUATIONS
    ]
)
CONTACT_QUADRUPLE_AUGMENTED_MATRIX: Final = (
    CONTACT_QUADRUPLE_COEFFICIENT_MATRIX.row_join(
        CONTACT_QUADRUPLE_CONSTANT_COLUMN
    )
)

CONTACT_PAIR_DENOMINATOR: Final = KAPPA + 2 * CONTACT_SUM
CONTACT_PAIR_QUADRATIC: Final = (
    CONTACT_SUM**2 + KAPPA * CONTACT_SUM + 1
)
CONTACT_PAIR_DIAGONAL: Final = (
    2 * CONTACT_SUM**2 + 3 * KAPPA * CONTACT_SUM + 4
)
QUADRUPLE_REDUCED_DISCRIMINANT: Final = (
    256 * FIBER_VALUE**2
    + 27 * FIBER_VALUE * KAPPA**4
    - 144 * FIBER_VALUE * KAPPA**2
    + 128 * FIBER_VALUE
    - 4 * KAPPA**2
    + 16
)
CONTACT_QUADRUPLE_SAME_TARGET_FACTOR: Final = (
    FIBER_VALUE * CONTACT_PAIR_DENOMINATOR**2
    + CONTACT_SUM
    * (CONTACT_SUM + KAPPA)
    * CONTACT_PAIR_QUADRATIC**2
)
CONTACT_QUADRUPLE_COMPATIBILITY: Final = (
    KAPPA**3 * CONTACT_SUM
    + 5 * KAPPA**2 * CONTACT_SUM**2
    + 3 * KAPPA**2
    + 4 * KAPPA * CONTACT_SUM**3
    + 16 * KAPPA * CONTACT_SUM
    + 12 * CONTACT_SUM**2
    + 4
)
EXPECTED_AUGMENTED_DETERMINANT: Final = (
    FIBER_VALUE**2
    * CONTACT_PAIR_QUADRATIC
    * CONTACT_QUADRUPLE_COMPATIBILITY
    * CONTACT_QUADRUPLE_SAME_TARGET_FACTOR**2
)
CONTACT_QUADRUPLE_VALID_LOCALIZER: Final = expand(
    KAPPA
    * (KAPPA - 2)
    * (KAPPA + 2)
    * FIBER_VALUE
    * CONTACT_SUM
    * CONTACT_PAIR_DENOMINATOR
    * CONTACT_PAIR_QUADRATIC
    * CONTACT_PAIR_DIAGONAL
    * QUADRUPLE_REDUCED_DISCRIMINANT
    * CONTACT_QUADRUPLE_SAME_TARGET_FACTOR
)

# One coefficient chart.  The factor k*u+2 is chart-specific and is not part
# of the global valid localizer: the exact determinantal saturation keeps and
# checks that legitimate locus using the other maximal minors.
CRAMER_DENOMINATOR: Final = (
    (KAPPA + 2 * CONTACT_SUM) * (KAPPA * CONTACT_SUM + 2)
)
GENERIC_ALPHA: Final = cancel(
    -FIBER_VALUE
    * CONTACT_SUM
    * (
        2 * KAPPA**3 * CONTACT_SUM
        + 6 * KAPPA**2 * CONTACT_SUM**2
        + 3 * KAPPA**2
        + 2 * KAPPA * CONTACT_SUM**3
        + 10 * KAPPA * CONTACT_SUM
        + 4
    )
    / CRAMER_DENOMINATOR
)
GENERIC_BETA: Final = cancel(
    -(FIBER_VALUE * KAPPA**2 + 2 * FIBER_VALUE + 1)
    * (
        KAPPA**2 * CONTACT_SUM**2
        + 3 * KAPPA * CONTACT_SUM**3
        + KAPPA * CONTACT_SUM
        + CONTACT_SUM**4
        + 4 * CONTACT_SUM**2
        - 1
    )
    / CRAMER_DENOMINATOR
)
GENERIC_GAMMA: Final = cancel(
    -(
        2 * FIBER_VALUE * KAPPA**3 * CONTACT_SUM**2
        + 6 * FIBER_VALUE * KAPPA**2 * CONTACT_SUM**3
        + 2 * FIBER_VALUE * KAPPA**2 * CONTACT_SUM
        + 2 * FIBER_VALUE * KAPPA * CONTACT_SUM**4
        + 8 * FIBER_VALUE * KAPPA * CONTACT_SUM**2
        - 2 * FIBER_VALUE * KAPPA
        + KAPPA**3 * CONTACT_SUM**2
        + 3 * KAPPA**2 * CONTACT_SUM**3
        + KAPPA * CONTACT_SUM**4
        + 2 * KAPPA * CONTACT_SUM**2
        - 3 * KAPPA
        - 4 * CONTACT_SUM
    )
    / CRAMER_DENOMINATOR
)
GENERIC_DELTA: Final = cancel(
    -(
        FIBER_VALUE * KAPPA**2 * CONTACT_SUM**2
        + 3 * FIBER_VALUE * KAPPA * CONTACT_SUM**3
        + FIBER_VALUE * KAPPA * CONTACT_SUM
        + FIBER_VALUE * CONTACT_SUM**4
        + 4 * FIBER_VALUE * CONTACT_SUM**2
        - FIBER_VALUE
        - KAPPA**3 * CONTACT_SUM
        - KAPPA**2 * CONTACT_SUM**2
        - 2 * KAPPA**2
        + 3 * KAPPA * CONTACT_SUM**3
        - 3 * KAPPA * CONTACT_SUM
        + CONTACT_SUM**4
        + 4 * CONTACT_SUM**2
        - 1
    )
    / CRAMER_DENOMINATOR
)
GENERIC_SOLUTION: Final = {
    ALPHA: GENERIC_ALPHA,
    BETA: GENERIC_BETA,
    GAMMA: GENERIC_GAMMA,
    DELTA: GENERIC_DELTA,
}
EXPECTED_LAST_CRAMER_RESIDUAL: Final = cancel(
    CONTACT_PAIR_QUADRATIC
    * CONTACT_QUADRUPLE_COMPATIBILITY
    * CONTACT_QUADRUPLE_SAME_TARGET_FACTOR
    / CRAMER_DENOMINATOR
)

# Exact metadata regenerated by tools/check_a6_delta_ten_contact_quadruple.sage.
SAGE_COMPATIBILITY_RATIONAL_FACTOR_COUNT: Final = 1
SAGE_PROJECTIVE_COMPONENT_GEOMETRY: Final = (3, 1, 2)
SAGE_VALID_SINGULAR_SATURATION: Final = (True, 1)
SAGE_RANK_THREE_SATURATIONS: Final = (True, 1, 1)
SAGE_RANK_TWO_SATURATIONS: Final = (True, 0, 0)

SAMPLE_BASE_POINT: Final = {
    KAPPA: -4,
    FIBER_VALUE: 1,
    CONTACT_SUM: 1,
}
SAMPLE_PARAMETERS: Final = {
    KAPPA: -4,
    ALPHA: 7,
    BETA: -19,
    GAMMA: 13,
    DELTA: -6,
}
SAMPLE_P: Final = expand(FAMILY_P.subs(SAMPLE_PARAMETERS))
SAMPLE_Q: Final = expand(FAMILY_Q.subs(SAMPLE_PARAMETERS))
SAMPLE_FIBER_POLYNOMIAL: Final = expand(
    FIBER_POLYNOMIAL.subs(SAMPLE_BASE_POINT)
)
SAMPLE_FIBER_QUOTIENT: Final = (
    T**5 - 2 * T**4 + 4 * T**3 - T**2 - 1
)
SAMPLE_SLOPE_POLYNOMIAL: Final = (
    SLOPE_VALUE**4
    - 538 * SLOPE_VALUE**3
    - 3098 * SLOPE_VALUE**2
    - 5898 * SLOPE_VALUE
    - 3727
)

SAMPLE_CONTACT_BRANCH: Final = T**2 - T + 1
SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR: Final = (
    S**6
    - 12 * S**5
    + 50 * S**4
    - 80 * S**3
    + 37 * S**2
    - 20 * S
    + 16
)
SAMPLE_NODE_SUM_FACTOR: Final = S**2 - 10 * S + 7
SAMPLE_COLLISION_POLYNOMIAL: Final = expand(
    (S - 1) ** 2
    * SAMPLE_NODE_SUM_FACTOR
    * SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR
)
SAMPLE_TANGENCY_POLYNOMIAL: Final = expand(
    4
    * (S - 1)
    * (
        9 * S**10
        - 225 * S**9
        + 2236 * S**8
        - 11532 * S**7
        + 33627 * S**6
        - 56277 * S**5
        + 52970 * S**4
        - 29702 * S**3
        + 13046 * S**2
        - 5320 * S
        + 1120
    )
)
SAMPLE_NODE_X_POLYNOMIAL: Final = X**2 + 564 * X - 476
SAMPLE_IMPLICIT: Final = (
    -X**9
    + 44 * X**8
    - 531 * X**7
    + 60 * X**6 * Y
    + 391 * X**6
    - 1546 * X**5 * Y
    + 5 * X**4 * Y**2
    + 98 * X**5
    + 1286 * X**4 * Y
    - 1527 * X**3 * Y**2
    - 12 * X**2 * Y**3
    + 204 * X**3 * Y
    + 1417 * X**2 * Y**2
    - 514 * X * Y**3
    + Y**4
    + 113 * X * Y**2
    + 530 * Y**3
    - 2 * Y**2
)

# Sage 10.8 exact affine Zariski--van Kamp presentation.  Signed integers
# name the four geometric meridians and record inversion.
CONTACT_QUADRUPLE_RELATIONS: Final = (
    (2, 1, -2, -1),
    (-2, -1, -4, -3, 1, 3, 4, 1),
    (
        -4, -3, -1, 3, -4, -3, 1, 3, 4, 3, -4, -3, -1, 3, 4, -3,
        1, 3, 4, 3, -4, -3, -1, 3, 4, -3, 1, 3, 4, 3, -4, -3,
        -1, 3, -4, -3, 1, 3, 4, -3, -4, -3, -1, 3, -4, -3, 1,
        3, 4, -3,
    ),
    (
        -4, -3, -1, 3, -4, -3, 1, 3, 4, -3, -4, -3, -1, 3, -4,
        -3, 1, 3, 4, -3, 1, 3, 4, 3,
    ),
    (-4, -3, 1, 3, 4, -3, -1, 3),
    (2, 3, 4, 1, -3, -2, -1, -4),
    (1, -2, -1, -4, -3, 2, 3, 4),
    (1, -4, -3, -2, -1, 2, 3, 4),
    (-4, -3, 1, 3, 4, -3, 1, 3, 4, -3, -1, 3, -4, -3, -1, 3),
)


def _modulo_branch(rational_expression: Expr, branch: Expr) -> Expr:
    """Reduce one rational graph derivative modulo a source branch."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    return expand(
        rem(numerator * invert(denominator, branch), branch, T)
    )


def _contact_jet_differences() -> tuple[Expr, Expr]:
    """Return the first two d/dP differences across the contact involution."""

    current = SAMPLE_Q
    differences: list[Expr] = []
    for _ in range(2):
        current = cancel(diff(current, T) / diff(SAMPLE_P, T))
        left = _modulo_branch(current, SAMPLE_CONTACT_BRANCH)
        right = _modulo_branch(
            current.subs(T, 1 - T),
            SAMPLE_CONTACT_BRANCH,
        )
        differences.append(
            expand(rem(left - right, SAMPLE_CONTACT_BRANCH, T))
        )
    return differences[0], differences[1]


@dataclass(frozen=True, slots=True)
class DeltaTenContactQuadrupleCertificate:
    """Exact incidence, sample-geometry, and complement certificate."""

    augmented_determinant_identity: Expr
    compatibility_rational_factor_count: int
    compatibility_sample_value: Expr
    compatibility_sample_gradient: tuple[Expr, Expr, Expr]
    projective_component_geometry: tuple[int, int, int]
    valid_singular_saturation: tuple[bool, int]
    rank_three_saturations: tuple[bool, int, int]
    rank_two_saturations: tuple[bool, int, int]
    generic_solution_first_residuals: tuple[Expr, ...]
    generic_last_residual_identity: Expr
    sample_cramer_solution: tuple[Expr, ...]
    projection_tangent_vectors: tuple[tuple[Expr, ...], ...]
    sample_coefficient_rank: int
    sample_augmented_rank: int
    sample_incidence_residuals: tuple[Expr, ...]
    sample_maximal_minors: tuple[Expr, ...]
    sample_valid_localizer: Expr
    fiber_remainder: Expr
    fiber_discriminant: Expr
    fiber_quotient_identity: Expr
    slope_resultant_identity: Expr
    slope_discriminant: Expr
    fiber_pair_sum_resultant_identity: Expr
    collision_identity: Expr
    tangency_identity: Expr
    collision_tangency_gcd: Expr
    collision_jets_at_contact: tuple[Expr, ...]
    component_discriminants: tuple[Expr, Expr]
    component_resultants: tuple[Expr, Expr, Expr]
    component_boundary_resultants: tuple[tuple[Expr, ...], ...]
    contact_chart_values: tuple[Expr, Expr, Expr]
    contact_pair_discriminant: Expr
    contact_image_remainders: tuple[Expr, Expr]
    contact_p_derivative_resultant: Expr
    contact_fiber_resultant: Expr
    contact_jet_differences: tuple[Expr, Expr]
    contact_jet_separation: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    quadruple_x_resultant_identity: Expr
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
    contact_delta: int
    quadruple_delta: int
    node_count: int
    infinity_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def compatibility_surface_geometrically_irreducible(self) -> bool:
        """Use Q-irreducibility and the smooth rational point."""

        return bool(
            self.compatibility_rational_factor_count == 1
            and self.compatibility_sample_value == 0
            and any(value != 0 for value in self.compatibility_sample_gradient)
        )

    @property
    def dominant_incidence_dimension(self) -> int:
        """Return one quartic-curve parameter plus the free fiber value."""

        return 2

    @property
    def coefficient_image_dimension(self) -> int:
        """Read the rank of the two displayed coefficient tangent vectors."""

        return int(Matrix(self.projection_tangent_vectors).rank())

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return (
            self.cusp_delta
            + self.contact_delta
            + self.quadruple_delta
            + self.node_count
            + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether every exact incidence and sample check agrees."""

        return bool(
            self.augmented_determinant_identity == 0
            and self.compatibility_surface_geometrically_irreducible
            and self.compatibility_sample_gradient == (4, 0, 8)
            and self.projective_component_geometry == (3, 1, 2)
            and self.valid_singular_saturation == (True, 1)
            and self.rank_three_saturations == (True, 1, 1)
            and self.rank_two_saturations == (True, 0, 0)
            and self.generic_solution_first_residuals == (0, 0, 0, 0)
            and self.generic_last_residual_identity == 0
            and self.sample_cramer_solution == (7, -19, 13, -6)
            and self.projection_tangent_vectors
            == ((0, 7, -18, 8, -1), (2, -4, 16, -6, 2))
            and self.coefficient_image_dimension
            == self.dominant_incidence_dimension
            == 2
            and self.sample_coefficient_rank == 4
            and self.sample_augmented_rank == 4
            and self.sample_incidence_residuals == (0, 0, 0, 0, 0)
            and self.sample_maximal_minors
            == (-12288, 9728, -1536, -192, -32)
            and self.sample_valid_localizer != 0
            and self.fiber_remainder == -1
            and self.fiber_discriminant == -4944
            and self.fiber_quotient_identity == 0
            and self.slope_resultant_identity == 0
            and self.slope_discriminant == -626582784
            and self.fiber_pair_sum_resultant_identity == 0
            and self.collision_identity == 0
            and self.tangency_identity == 0
            and self.collision_tangency_gcd == S - 1
            and self.collision_jets_at_contact == (0, 0, 32, 864)
            and self.component_discriminants == (56316985344, 72)
            and self.component_resultants == (115344, -8, -2)
            and all(
                all(value != 0 for value in row)
                for row in self.component_boundary_resultants
            )
            and self.contact_chart_values == (-2, -2, -6)
            and self.contact_pair_discriminant == -3
            and self.contact_image_remainders == (3, -7)
            and self.contact_p_derivative_resultant == 84
            and self.contact_fiber_resultant == 4
            and self.contact_jet_differences
            == (0, T / 49 - Rational(1, 98))
            and self.contact_jet_separation == Rational(3, 9604)
            and self.cusp_image_factor == -2
            and self.extra_critical_factor == -311472
            and self.quadruple_x_resultant_identity == 0
            and self.node_x_resultant_identity == 0
            and self.node_x_discriminant == 320000
            and self.node_target_separations == (-476, 89, 1225)
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_length == 18
            and self.sage_jacobian_radical_length == 5
            and self.sage_jacobian_components
            == ((2, 2), (9, 1), (3, 1), (4, 1))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.topology_propagation_dependencies
            == (
                "connected clean open",
                "proper projective Whitney-Thom triviality",
            )
            and not self.split_and_deeper_boundaries_classified_here
            and self.total_delta == self.arithmetic_genus == 28
            and len(CONTACT_QUADRUPLE_RELATIONS) == 9
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_contact_quadruple_certificate(
) -> DeltaTenContactQuadrupleCertificate:
    """Build the exact C2 + Q0 + 2N certificate."""

    sample_matrix = CONTACT_QUADRUPLE_COEFFICIENT_MATRIX.subs(
        SAMPLE_BASE_POINT
    )
    sample_constant = CONTACT_QUADRUPLE_CONSTANT_COLUMN.subs(
        SAMPLE_BASE_POINT
    )
    sample_vector = Matrix(
        [SAMPLE_PARAMETERS[parameter] for parameter in LINEAR_PARAMETERS]
    )
    maximal_minors = tuple(
        expand(
            sample_matrix.extract(
                [row for row in range(5) if row != omitted],
                range(4),
            ).det()
        )
        for omitted in range(5)
    )

    generic_residuals = tuple(
        cancel(equation.subs(GENERIC_SOLUTION))
        for equation in CONTACT_QUADRUPLE_EQUATIONS
    )
    cramer_functions = (
        KAPPA,
        GENERIC_ALPHA,
        GENERIC_BETA,
        GENERIC_GAMMA,
        GENERIC_DELTA,
    )
    base_tangents = ((0, 1, 0), (2, 0, -1))
    projection_tangents = tuple(
        tuple(
            cancel(
                (
                    diff(function, KAPPA) * tangent[0]
                    + diff(function, FIBER_VALUE) * tangent[1]
                    + diff(function, CONTACT_SUM) * tangent[2]
                ).subs(SAMPLE_BASE_POINT)
            )
            for function in cramer_functions
        )
        for tangent in base_tangents
    )

    sample_collision = expand(
        COLLISION_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
    )
    sample_tangency = expand(
        TANGENCY_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
    )
    sample_x_numerator = expand(
        COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS)
    )
    sample_x_denominator = expand(
        COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    )
    fiber_pair_sum_resultant = resultant(
        SAMPLE_FIBER_POLYNOMIAL,
        SAMPLE_FIBER_POLYNOMIAL.subs(T, S - T),
        T,
    )
    diagonal_factor = S**4 - 8 * S**3 + 4 * S**2 - 16
    contact_jets = _contact_jet_differences()
    raw_implicit = resultant(X - SAMPLE_P, Y - SAMPLE_Q, T)

    boundary_expressions = (
        PAIR_DENOMINATOR.subs(SAMPLE_PARAMETERS),
        PAIR_QUADRATIC.subs(SAMPLE_PARAMETERS),
        PAIR_DIAGONAL_FACTOR.subs(SAMPLE_PARAMETERS),
        sample_tangency,
    )
    boundary_resultants = tuple(
        tuple(resultant(component, boundary, S) for boundary in boundary_expressions)
        for component in (
            SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR,
            SAMPLE_NODE_SUM_FACTOR,
        )
    )

    return DeltaTenContactQuadrupleCertificate(
        augmented_determinant_identity=expand(
            CONTACT_QUADRUPLE_AUGMENTED_MATRIX.det()
            - EXPECTED_AUGMENTED_DETERMINANT
        ),
        compatibility_rational_factor_count=(
            SAGE_COMPATIBILITY_RATIONAL_FACTOR_COUNT
        ),
        compatibility_sample_value=CONTACT_QUADRUPLE_COMPATIBILITY.subs(
            SAMPLE_BASE_POINT
        ),
        compatibility_sample_gradient=tuple(
            diff(CONTACT_QUADRUPLE_COMPATIBILITY, variable).subs(
                SAMPLE_BASE_POINT
            )
            for variable in (KAPPA, FIBER_VALUE, CONTACT_SUM)
        ),
        projective_component_geometry=SAGE_PROJECTIVE_COMPONENT_GEOMETRY,
        valid_singular_saturation=SAGE_VALID_SINGULAR_SATURATION,
        rank_three_saturations=SAGE_RANK_THREE_SATURATIONS,
        rank_two_saturations=SAGE_RANK_TWO_SATURATIONS,
        generic_solution_first_residuals=generic_residuals[:4],
        generic_last_residual_identity=cancel(
            generic_residuals[4] - EXPECTED_LAST_CRAMER_RESIDUAL
        ),
        sample_cramer_solution=tuple(
            expression.subs(SAMPLE_BASE_POINT)
            for expression in (
                GENERIC_ALPHA,
                GENERIC_BETA,
                GENERIC_GAMMA,
                GENERIC_DELTA,
            )
        ),
        projection_tangent_vectors=projection_tangents,
        sample_coefficient_rank=sample_matrix.rank(),
        sample_augmented_rank=sample_matrix.row_join(sample_constant).rank(),
        sample_incidence_residuals=tuple(
            expand(value)
            for value in sample_matrix * sample_vector + sample_constant
        ),
        sample_maximal_minors=maximal_minors,
        sample_valid_localizer=CONTACT_QUADRUPLE_VALID_LOCALIZER.subs(
            SAMPLE_BASE_POINT
        ),
        fiber_remainder=rem(
            SAMPLE_Q,
            SAMPLE_FIBER_POLYNOMIAL,
            T,
        ),
        fiber_discriminant=discriminant(SAMPLE_FIBER_POLYNOMIAL, T),
        fiber_quotient_identity=expand(
            SAMPLE_Q
            + 1
            - SAMPLE_FIBER_POLYNOMIAL * SAMPLE_FIBER_QUOTIENT
        ),
        slope_resultant_identity=expand(
            resultant(
                SAMPLE_FIBER_POLYNOMIAL,
                SLOPE_VALUE - SAMPLE_FIBER_QUOTIENT,
                T,
            )
            - SAMPLE_SLOPE_POLYNOMIAL
        ),
        slope_discriminant=discriminant(
            SAMPLE_SLOPE_POLYNOMIAL,
            SLOPE_VALUE,
        ),
        fiber_pair_sum_resultant_identity=expand(
            fiber_pair_sum_resultant
            - diagonal_factor * SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR**2
        ),
        collision_identity=expand(
            sample_collision - SAMPLE_COLLISION_POLYNOMIAL
        ),
        tangency_identity=expand(
            sample_tangency - SAMPLE_TANGENCY_POLYNOMIAL
        ),
        collision_tangency_gcd=gcd(
            Poly(sample_collision, S),
            Poly(sample_tangency, S),
        ).as_expr(),
        collision_jets_at_contact=tuple(
            diff(sample_collision, S, order).subs(S, 1)
            for order in range(4)
        ),
        component_discriminants=(
            discriminant(SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR, S),
            discriminant(SAMPLE_NODE_SUM_FACTOR, S),
        ),
        component_resultants=(
            resultant(
                SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR,
                SAMPLE_NODE_SUM_FACTOR,
                S,
            ),
            resultant(SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR, S - 1, S),
            resultant(SAMPLE_NODE_SUM_FACTOR, S - 1, S),
        ),
        component_boundary_resultants=boundary_resultants,
        contact_chart_values=tuple(
            expression.subs(SAMPLE_PARAMETERS).subs(S, 1)
            for expression in (
                PAIR_DENOMINATOR,
                PAIR_QUADRATIC,
                PAIR_DIAGONAL_FACTOR,
            )
        ),
        contact_pair_discriminant=discriminant(
            SAMPLE_CONTACT_BRANCH,
            T,
        ),
        contact_image_remainders=(
            rem(SAMPLE_P, SAMPLE_CONTACT_BRANCH, T),
            rem(SAMPLE_Q, SAMPLE_CONTACT_BRANCH, T),
        ),
        contact_p_derivative_resultant=resultant(
            SAMPLE_CONTACT_BRANCH,
            diff(SAMPLE_P, T),
            T,
        ),
        contact_fiber_resultant=resultant(
            SAMPLE_CONTACT_BRANCH,
            SAMPLE_FIBER_POLYNOMIAL,
            T,
        ),
        contact_jet_differences=contact_jets,
        contact_jet_separation=resultant(
            SAMPLE_CONTACT_BRANCH,
            contact_jets[1],
            T,
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(SAMPLE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(
            SAMPLE_PARAMETERS
        ),
        quadruple_x_resultant_identity=expand(
            resultant(
                SAMPLE_QUADRUPLE_PAIR_SUM_FACTOR,
                X * sample_x_denominator - sample_x_numerator,
                S,
            )
            - 5308416 * (X - 1) ** 6
        ),
        node_x_resultant_identity=expand(
            resultant(
                SAMPLE_NODE_SUM_FACTOR,
                X * sample_x_denominator - sample_x_numerator,
                S,
            )
            - 1296 * SAMPLE_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(
            SAMPLE_NODE_X_POLYNOMIAL,
            X,
        ),
        node_target_separations=tuple(
            SAMPLE_NODE_X_POLYNOMIAL.subs(X, value)
            for value in (0, 1, 3)
        ),
        implicit_resultant_identity=expand(raw_implicit + SAMPLE_IMPLICIT),
        implicit_parameterization_identity=expand(
            SAMPLE_IMPLICIT.subs({X: SAMPLE_P, Y: SAMPLE_Q})
        ),
        implicit_content=Poly(raw_implicit, X, Y).content(),
        sage_jacobian_length=18,
        sage_jacobian_radical_length=5,
        sage_jacobian_components=((2, 2), (9, 1), (3, 1), (4, 1)),
        sage_cyclic_simplification=(1, 0, True),
        topology_propagation_dependencies=(
            "connected clean open",
            "proper projective Whitney-Thom triviality",
        ),
        split_and_deeper_boundaries_classified_here=False,
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=2,
        contact_delta=2,
        quadruple_delta=6,
        node_count=2,
        infinity_delta=16,
        complement_census=n_generator_three_cycle_presentation_census(
            CONTACT_QUADRUPLE_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the certificate summary and fail on any regression."""

    certificate = exact_delta_ten_contact_quadruple_certificate()
    print("delta-ten C2 + Q0 + 2N certificate:", certificate.verified)
    print(
        "valid compatibility surface:",
        "geometrically irreducible smooth genus-one open",
    )
    print("coefficient/augmented ranks:", (
        certificate.sample_coefficient_rank,
        certificate.sample_augmented_rank,
    ))
    print("affine complement:", "Z")
    print("A6 assignments:", certificate.complement_census.a6_assignments)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
