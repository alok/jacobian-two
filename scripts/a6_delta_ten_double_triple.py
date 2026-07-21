"""Audit the conditional delta-ten ``T111^2 + 4N`` profile.

For one ordinary triple fiber of the quartic

``P = t^2 + k*t^3 + t^4``, choose the fourth, omitted root ``e``.  The other
three roots form the cubic ``(P(t)-P(e))/(t-e)``.  They have one common
``Q``-value exactly when the linear and quadratic coefficients of the
remainder of ``Q`` modulo that cubic vanish.  Two omitted roots ``u,v``
therefore give four affine-linear equations in ``(a,b,c,d)``.

The coefficient determinant is factored exactly.  Its dense valid open gives
one rational three-dimensional graph, symmetric under ``u <-> v``.  The
same-fiber factor is invalid for two separate triple targets, while the one
remaining rank factor is generically inconsistent by an exact augmented
minor.  Compatible lower-dimensional pieces of that rank factor and split
boundaries remain a separate audit; the theorem here concerns the component
dominating the valid Cramer open.

The rational member with ``(k,a,b,c,d)=(3,5,4/3,6,5)`` has two ordinary
triple fibers and four residual nodes.  The exact algebra below is
dependency-free.  Complement topology and propagation are explicit separate
computer-assisted dependencies, and no claim about unrestricted ``A6``,
``S6``, or the plane Jacobian conjecture is made.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import (
    Expr,
    Poly,
    Rational,
    Symbol,
    cancel,
    discriminant,
    expand,
    linear_eq_to_matrix,
    rem,
    resultant,
)

from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    n_generator_three_cycle_presentation_census,
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
)

OMITTED_ROOT: Final = Symbol("e")
FIRST_OMITTED_ROOT: Final = Symbol("u")
SECOND_OMITTED_ROOT: Final = Symbol("v")
SLOPE_VALUE: Final = Symbol("m_double_triple")
TARGET_X: Final = Symbol("x_double_triple")
COEFFICIENTS: Final = (ALPHA, BETA, GAMMA, DELTA)


def omitted_fiber_cubic(root: Expr) -> Expr:
    """Return the cubic containing the other roots of the ``P(root)`` fiber."""

    return cancel((FAMILY_P - FAMILY_P.subs(T, root)) / (T - root))


GENERIC_OMITTED_CUBIC: Final = omitted_fiber_cubic(OMITTED_ROOT)
GENERIC_Q_REMAINDER: Final = expand(rem(FAMILY_Q, GENERIC_OMITTED_CUBIC, T))
GENERIC_TRIPLE_EQUATIONS: Final = tuple(
    Poly(GENERIC_Q_REMAINDER, T).coeff_monomial(T**degree) for degree in (1, 2)
)
DOUBLE_TRIPLE_EQUATIONS: Final = tuple(
    equation.subs(OMITTED_ROOT, root)
    for root in (FIRST_OMITTED_ROOT, SECOND_OMITTED_ROOT)
    for equation in GENERIC_TRIPLE_EQUATIONS
)
INCIDENCE_MATRIX, INCIDENCE_RIGHT_HAND_SIDE = linear_eq_to_matrix(
    DOUBLE_TRIPLE_EQUATIONS,
    COEFFICIENTS,
)
INCIDENCE_AUGMENTED_MATRIX: Final = INCIDENCE_MATRIX.row_join(INCIDENCE_RIGHT_HAND_SIDE)

FIRST_CUSP_FIBER_FACTOR: Final = FIRST_OMITTED_ROOT**2 + KAPPA * FIRST_OMITTED_ROOT + 1
SECOND_CUSP_FIBER_FACTOR: Final = (
    SECOND_OMITTED_ROOT**2 + KAPPA * SECOND_OMITTED_ROOT + 1
)
RESIDUAL_RANK_FACTOR: Final = (
    KAPPA * FIRST_OMITTED_ROOT * SECOND_OMITTED_ROOT
    - FIRST_OMITTED_ROOT
    - SECOND_OMITTED_ROOT
)
SAME_FIBER_FACTOR: Final = (
    KAPPA * FIRST_OMITTED_ROOT**2
    + KAPPA * FIRST_OMITTED_ROOT * SECOND_OMITTED_ROOT
    + KAPPA * SECOND_OMITTED_ROOT**2
    + FIRST_OMITTED_ROOT**3
    + FIRST_OMITTED_ROOT**2 * SECOND_OMITTED_ROOT
    + FIRST_OMITTED_ROOT * SECOND_OMITTED_ROOT**2
    + FIRST_OMITTED_ROOT
    + SECOND_OMITTED_ROOT**3
    + SECOND_OMITTED_ROOT
)
EXPECTED_INCIDENCE_DETERMINANT: Final = expand(
    -(
        (FIRST_OMITTED_ROOT - SECOND_OMITTED_ROOT) ** 2
        * FIRST_CUSP_FIBER_FACTOR**2
        * SECOND_CUSP_FIBER_FACTOR**2
        * RESIDUAL_RANK_FACTOR
        * SAME_FIBER_FACTOR
    )
)

RANKDROP_TEST_POINT: Final = {
    KAPPA: Rational(3, 2),
    FIRST_OMITTED_ROOT: 1,
    SECOND_OMITTED_ROOT: 2,
}

SAMPLE_PARAMETERS: Final = {
    KAPPA: 3,
    ALPHA: 5,
    BETA: Rational(4, 3),
    GAMMA: 6,
    DELTA: 5,
}
SAMPLE_P: Final = expand(FAMILY_P.subs(SAMPLE_PARAMETERS))
SAMPLE_Q: Final = expand(FAMILY_Q.subs(SAMPLE_PARAMETERS))
SAMPLE_OMITTED_ROOTS: Final = (Rational(-1), Rational(1))
SAMPLE_TRIPLE_CUBICS: Final = tuple(
    expand(omitted_fiber_cubic(root).subs(SAMPLE_PARAMETERS))
    for root in SAMPLE_OMITTED_ROOTS
)
SAMPLE_TRIPLE_VALUES: Final = tuple(
    expand(rem(SAMPLE_Q, cubic, T)) for cubic in SAMPLE_TRIPLE_CUBICS
)

FIRST_TRIPLE_PAIR_FACTOR: Final = S**3 + 4 * S**2 + 3 * S - 3
SECOND_TRIPLE_PAIR_FACTOR: Final = S**3 + 8 * S**2 + 21 * S + 15
RESIDUAL_NODE_SUM_FACTOR: Final = 3 * S**4 + 18 * S**3 + 27 * S**2 + 5 * S - 3
SAMPLE_COLLISION_POLYNOMIAL: Final = expand(
    COLLISION_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
)
RESIDUAL_NODE_X_POLYNOMIAL: Final = (
    81 * TARGET_X**4 - 81 * TARGET_X**3 - 459 * TARGET_X**2 - 59 * TARGET_X + 6
)
RESIDUAL_NODE_X_ELIMINANT: Final = resultant(
    RESIDUAL_NODE_SUM_FACTOR,
    TARGET_X * COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    - COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS),
    S,
)

DOUBLE_TRIPLE_RELATIONS: Final = (
    (4, -3),
    (4, 3, -4, -3),
    (1, 3, 4, -3, -4, -3, -2, -1, 2, 3, 4, 3, -4, -3),
    (4, -3, -4, -3, 1, 2, 3, 4, 3, -4, -3, -2, -1, 3),
    (2, 1, -2, -1),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, -2, -1, 3, -4, -3, 1, 2, 1, 2, 1),
    (-2, -1, 3, -4, -3, 1, 2, 1, -2, -1, 3, 4, -3, 1, 2, -1),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-4, 1, 3, 4, -3, -1),
    (-3, 4, 1, 3, -1, -4),
)


@dataclass(frozen=True, slots=True)
class DeltaTenDoubleTripleCertificate:
    """Exact dominant-chart and rational-member data for ``T111^2 + 4N``."""

    incidence_determinant_identity: Expr
    same_fiber_identity: Expr
    rankdrop_coefficient_rank: int
    rankdrop_augmented_rank: int
    rankdrop_augmented_minor: Expr
    sample_equation_residuals: tuple[Expr, ...]
    sample_cubic_discriminants: tuple[Expr, Expr]
    sample_triple_values: tuple[Expr, Expr]
    sample_omitted_value_differences: tuple[Expr, Expr]
    sample_slope_discriminants: tuple[Expr, Expr]
    pair_sum_identities: tuple[Expr, Expr]
    collision_factor_identity: Expr
    collision_factor_discriminants: tuple[Expr, Expr, Expr]
    collision_factor_resultants: tuple[Expr, Expr, Expr]
    validity_resultants: tuple[tuple[Expr, Expr, Expr, Expr], ...]
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    node_x_eliminant_identity: Expr
    node_x_discriminant: Expr
    node_target_separations: tuple[Expr, Expr, Expr]
    sage_jacobian_components: tuple[tuple[int, int], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    complement_census: ThreeCyclePresentationCensus
    arithmetic_genus: int
    total_delta: int

    @property
    def verified(self) -> bool:
        """Whether every exact double-triple check agrees."""

        return bool(
            self.incidence_determinant_identity == 0
            and self.same_fiber_identity == 0
            and self.rankdrop_coefficient_rank == 3
            and self.rankdrop_augmented_rank == 4
            and self.rankdrop_augmented_minor == -938448
            and self.sample_equation_residuals == (0, 0, 0, 0)
            and self.sample_cubic_discriminants == (-87, -255)
            and self.sample_triple_values == (Rational(1, 3), Rational(-125, 3))
            and self.sample_omitted_value_differences == (-6, 60)
            and self.sample_slope_discriminants
            == (Rational(-1560896, 3), Rational(-685443400000, 3))
            and self.pair_sum_identities == (0, 0)
            and self.collision_factor_identity == 0
            and self.collision_factor_discriminants == (-87, -255, 2108025)
            and self.collision_factor_resultants == (360, 2175, 33675)
            and all(
                all(value != 0 for value in row) for row in self.validity_resultants
            )
            and self.cusp_image_factor == Rational(-5, 9)
            and self.extra_critical_factor == 6409
            and self.node_x_eliminant_identity == 0
            and self.node_x_discriminant == 402849850982400
            and self.node_target_separations == (-232, 28736, 6)
            and self.sage_jacobian_components == ((4, 1), (4, 1), (4, 1), (4, 4))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
            and self.arithmetic_genus == 28
            and self.total_delta == 28
        )


@cache
def exact_delta_ten_double_triple_certificate() -> DeltaTenDoubleTripleCertificate:
    """Build the exact dominant-chart and sample certificate."""

    incidence_determinant = expand(INCIDENCE_MATRIX.det())
    rankdrop_matrix = INCIDENCE_MATRIX.subs(RANKDROP_TEST_POINT)
    rankdrop_augmented = INCIDENCE_AUGMENTED_MATRIX.subs(RANKDROP_TEST_POINT)

    sample_equations = tuple(
        expand(
            equation.subs(SAMPLE_PARAMETERS).subs(
                {
                    FIRST_OMITTED_ROOT: SAMPLE_OMITTED_ROOTS[0],
                    SECOND_OMITTED_ROOT: SAMPLE_OMITTED_ROOTS[1],
                }
            )
        )
        for equation in DOUBLE_TRIPLE_EQUATIONS
    )
    triple_data: list[tuple[Expr, Expr, Expr]] = []
    for root, cubic, triple_value in zip(
        SAMPLE_OMITTED_ROOTS,
        SAMPLE_TRIPLE_CUBICS,
        SAMPLE_TRIPLE_VALUES,
        strict=True,
    ):
        quotient = cancel((SAMPLE_Q - triple_value) / cubic)
        slope_resultant = resultant(
            cubic,
            SLOPE_VALUE * (T - root) - quotient,
            T,
        )
        triple_data.append(
            (
                discriminant(cubic, T),
                expand(SAMPLE_Q.subs(T, root) - triple_value),
                discriminant(slope_resultant, SLOPE_VALUE),
            )
        )

    pair_factors = (
        FIRST_TRIPLE_PAIR_FACTOR,
        SECOND_TRIPLE_PAIR_FACTOR,
        RESIDUAL_NODE_SUM_FACTOR,
    )
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(SAMPLE_PARAMETERS))
    sample_denominator = expand(PAIR_DENOMINATOR.subs(SAMPLE_PARAMETERS))
    sample_diagonal = expand(PAIR_DIAGONAL_FACTOR.subs(SAMPLE_PARAMETERS))
    sample_cusp_pair = expand(PAIR_QUADRATIC.subs(SAMPLE_PARAMETERS))

    return DeltaTenDoubleTripleCertificate(
        incidence_determinant_identity=expand(
            incidence_determinant - EXPECTED_INCIDENCE_DETERMINANT
        ),
        same_fiber_identity=expand(
            SAMPLE_P.subs(T, FIRST_OMITTED_ROOT)
            - SAMPLE_P.subs(T, SECOND_OMITTED_ROOT)
            - (FIRST_OMITTED_ROOT - SECOND_OMITTED_ROOT)
            * SAME_FIBER_FACTOR.subs(KAPPA, 3)
        ),
        rankdrop_coefficient_rank=rankdrop_matrix.rank(),
        rankdrop_augmented_rank=rankdrop_augmented.rank(),
        rankdrop_augmented_minor=rankdrop_augmented[:, [1, 2, 3, 4]].det(),
        sample_equation_residuals=sample_equations,
        sample_cubic_discriminants=tuple(row[0] for row in triple_data),
        sample_triple_values=SAMPLE_TRIPLE_VALUES,
        sample_omitted_value_differences=tuple(row[1] for row in triple_data),
        sample_slope_discriminants=tuple(row[2] for row in triple_data),
        pair_sum_identities=(
            expand(FIRST_TRIPLE_PAIR_FACTOR + SAMPLE_TRIPLE_CUBICS[0].subs(T, -2 - S)),
            expand(SECOND_TRIPLE_PAIR_FACTOR + SAMPLE_TRIPLE_CUBICS[1].subs(T, -4 - S)),
        ),
        collision_factor_identity=expand(
            3 * SAMPLE_COLLISION_POLYNOMIAL
            - FIRST_TRIPLE_PAIR_FACTOR
            * SECOND_TRIPLE_PAIR_FACTOR
            * RESIDUAL_NODE_SUM_FACTOR
        ),
        collision_factor_discriminants=tuple(
            discriminant(factor, S) for factor in pair_factors
        ),
        collision_factor_resultants=(
            resultant(FIRST_TRIPLE_PAIR_FACTOR, SECOND_TRIPLE_PAIR_FACTOR, S),
            resultant(FIRST_TRIPLE_PAIR_FACTOR, RESIDUAL_NODE_SUM_FACTOR, S),
            resultant(SECOND_TRIPLE_PAIR_FACTOR, RESIDUAL_NODE_SUM_FACTOR, S),
        ),
        validity_resultants=tuple(
            (
                resultant(factor, sample_denominator, S),
                resultant(factor, sample_diagonal, S),
                resultant(factor, sample_cusp_pair, S),
                resultant(factor, sample_tangency, S),
            )
            for factor in pair_factors
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(SAMPLE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(SAMPLE_PARAMETERS),
        node_x_eliminant_identity=expand(
            RESIDUAL_NODE_X_ELIMINANT - 5625 * RESIDUAL_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(
            RESIDUAL_NODE_X_POLYNOMIAL,
            TARGET_X,
        ),
        node_target_separations=(
            RESIDUAL_NODE_X_POLYNOMIAL.subs(TARGET_X, -1),
            RESIDUAL_NODE_X_POLYNOMIAL.subs(TARGET_X, 5),
            RESIDUAL_NODE_X_POLYNOMIAL.subs(TARGET_X, 0),
        ),
        sage_jacobian_components=((4, 1), (4, 1), (4, 1), (4, 4)),
        sage_cyclic_simplification=(1, 0, True),
        complement_census=n_generator_three_cycle_presentation_census(
            DOUBLE_TRIPLE_RELATIONS,
            4,
        ),
        arithmetic_genus=28,
        total_delta=2 + 3 + 3 + 4 + 16,
    )


def main() -> int:
    """Print the exact double-triple certificate summary."""

    certificate = exact_delta_ten_double_triple_certificate()
    print("delta-ten T111^2 + 4N certificate:", certificate.verified)
    print("incidence determinant identity:", certificate.incidence_determinant_identity)
    print(
        "rank-drop ranks:",
        certificate.rankdrop_coefficient_rank,
        certificate.rankdrop_augmented_rank,
    )
    print("sample triple values:", certificate.sample_triple_values)
    print("sample collision discriminants:", certificate.collision_factor_discriminants)
    print("sample complement census:", certificate.complement_census)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
