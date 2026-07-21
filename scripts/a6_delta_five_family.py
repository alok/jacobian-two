"""Conditional codimension-two reduction of the ``A6`` delta-five family.

Under the one-pair hypotheses isolated in ``a6_one_pair_infinity.py``, the
first large-link survivor has affine degrees ``(3,8)`` and collision delta
five.  Polynomial changes reduce it to a three-parameter family.  Exact
collision algebra and three Zariski--van Kamp presentations exclude the
generic stratum and both admissible codimension-one walls.  Only deeper
collision degenerations remain.

This is a conditional, computer-assisted family reduction.  It is not a full
elimination of delta five, the ``A6`` passport, or the plane conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from sympy import (
    Expr,
    Poly,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    factor_list,
    gcd,
    rem,
    resultant,
)

from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    three_cycle_presentation_census,
)

T: Final = Symbol("t")
U: Final = Symbol("u")
S: Final = Symbol("s")
X: Final = Symbol("X")
Y: Final = Symbol("Y")
ALPHA: Final = Symbol("alpha")
BETA: Final = Symbol("beta")
GAMMA: Final = Symbol("gamma")

FAMILY_P: Final = T**2 + T**3
FAMILY_Q: Final = (
    ALPHA * T**5 + BETA * T**6 + GAMMA * T**7 + T**8
)
COLLISION_POLYNOMIAL: Final = (
    S**5
    + (GAMMA + 2) * S**4
    + (4 * GAMMA - 2) * S**3
    + (-ALPHA + 2 * BETA + 3 * GAMMA - 4) * S**2
    + (-ALPHA + 3 * BETA - GAMMA) * S
    + ALPHA
)
EXTRA_CRITICAL_FACTOR: Final = (
    135 * ALPHA - 108 * BETA + 84 * GAMMA - 64
)
CUSP_COLLISION_FACTOR: Final = ALPHA - BETA + GAMMA - 1
TRIPLE_COLLISION_FACTOR: Final = (
    ALPHA * GAMMA
    - ALPHA
    - 2 * BETA * GAMMA
    + 3 * BETA
    + 2 * GAMMA**2
    - 5 * GAMMA
    + 3
)
INCIDENCE_ALPHA_NUMERATOR: Final = S**2 * (
    4 * GAMMA * S**2
    + 13 * GAMMA * S
    + 11 * GAMMA
    + 6 * S**3
    + 14 * S**2
    - 12
)
INCIDENCE_BETA_NUMERATOR: Final = (
    2 * GAMMA * S**4
    + 5 * GAMMA * S**3
    - GAMMA * S**2
    - 7 * GAMMA * S
    + GAMMA
    + 3 * S**5
    + 5 * S**4
    - 6 * S**3
    - 6 * S**2
    + 8 * S
)
TANGENCY_POLYNOMIAL: Final = S**2 * (S + 1) * (
    -15 * ALPHA * S**2
    - 15 * ALPHA * S
    + 10 * ALPHA
    + 30 * BETA * S**2
    + 42 * BETA * S
    + 21 * GAMMA * S**4
    + 77 * GAMMA * S**3
    + 56 * GAMMA * S**2
    - 14 * GAMMA * S
    + 24 * S**5
    + 48 * S**4
    - 32 * S**3
    - 72 * S**2
)

CLEAN_PARAMETERS: Final = {ALPHA: 1, BETA: 1, GAMMA: 0}
TRIPLE_PARAMETERS: Final = {ALPHA: 3, BETA: 0, GAMMA: 0}
TANGENCY_PARAMETERS: Final = {ALPHA: -4, BETA: 1, GAMMA: 3}
CONTACT_THREE_PARAMETERS: Final = {
    ALPHA: Rational(16, 5),
    BETA: Rational(32, 5),
    GAMMA: Rational(24, 5),
}
MIXED_PARAMETERS: Final = {
    ALPHA: Rational(20, 3),
    BETA: 1,
    GAMMA: Rational(2, 3),
}

CLEAN_Q: Final = FAMILY_Q.subs(CLEAN_PARAMETERS)
CLEAN_BRANCH_EQUATION: Final = (
    -X**8
    - 5 * X**7
    - 6 * X**6
    + 8 * X**5 * Y
    + X**5
    + 19 * X**4 * Y
    - 3 * X**3 * Y
    - 15 * X**2 * Y**2
    + 9 * X * Y**2
    + Y**3
    - Y**2
)
CLEAN_COLLISION_FACTOR: Final = (
    T**10
    + 2 * T**9
    + 4 * T**8
    + 6 * T**7
    + 10 * T**6
    + 10 * T**5
    + 11 * T**4
    + 13 * T**3
    + 6 * T**2
    - T
    - 1
)
CLEAN_NODE_COORDINATE_POLYNOMIAL: Final = (
    X**5 + 3 * X**4 - 3 * X**3 - 35 * X**2 + 12 * X - 1
)

# Sage 10.8 unsimplified affine Zariski--van Kamp presentations.  All three
# generators are geometric fiber meridians.
CLEAN_RELATIONS: Final = (
    (3, 2, -3, -2),
    (-3, -2, 3, 1),
    (-3, -2, -3, 2, 3, 1, -3, -2, 3, 2, 3, -1),
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (3, 2, -3, -2),
    (-3, -2, -3, -2, 1, 2, 3, 2, 3, -2, -3, -2, -1, 2, 3, 2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
)
TRIPLE_RELATIONS: Final = (
    (3, 2, -3, -2),
    (2, 3, 1, -3, -2, -1),
    (2, -1, -3, -2, 3, 1),
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -3, -2, 1, 2, 3, 2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
)
TANGENCY_RELATIONS: Final = (
    (2, 1, -2, -1),
    (3, 2, -3, -2),
    (
        -3,
        -2,
        -1,
        2,
        3,
        2,
        -3,
        -2,
        1,
        2,
        3,
        2,
        -3,
        -2,
        1,
        2,
        3,
        2,
        -3,
        -2,
        -1,
        2,
        3,
        -2,
        -3,
        -2,
        -1,
        2,
        3,
        -2,
    ),
    (-3, -2, -1, 2, 3, -2, -3, -2, -1, 2, 3, -2, 1, 2, 3, 2),
    (
        -3,
        -2,
        -1,
        2,
        -3,
        -2,
        -1,
        2,
        3,
        -2,
        1,
        2,
        3,
        2,
        -3,
        -2,
        -1,
        2,
        -3,
        -2,
        1,
        2,
        3,
        -2,
        1,
        2,
        3,
        -2,
    ),
    (-3, -2, 1, 2, 3, -2, 1, 2, 3, -2, -1, 2, -3, -2, -1, 2),
)


@dataclass(frozen=True, slots=True)
class A6DeltaFiveCleanCurveCertificate:
    """Exact geometry and complement-image census for the clean model."""

    implicit_resultant_identity: Expr
    parametrized_branch_identity: Expr
    off_diagonal_gcd: Expr
    collision_resultant_identity: Expr
    pairing_resultant_identity: Expr
    collision_discriminant: Expr
    pair_diagonal_resultant: Expr
    pair_product_resultant: Expr
    tangency_resultant: Expr
    node_coordinate_resultant_identity: Expr
    node_coordinate_discriminant: Expr
    cusp_image_separation: Expr
    residual_derivative_value: Expr
    finite_cusp_pair: tuple[int, int]
    node_count: int
    infinity_pair: tuple[int, int]
    arithmetic_genus: int
    total_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether all exact clean-curve checkpoints agree."""

        return bool(
            self.implicit_resultant_identity == 0
            and self.parametrized_branch_identity == 0
            and self.off_diagonal_gcd == 1
            and self.collision_resultant_identity == 0
            and self.pairing_resultant_identity == 0
            and self.collision_discriminant == -4903
            and self.pair_diagonal_resultant == 37
            and self.pair_product_resultant == -1
            and self.tangency_resultant == 181411
            and self.node_coordinate_resultant_identity == 0
            and self.node_coordinate_discriminant == -(5**6) * 4903
            and self.cusp_image_separation == -1
            and self.residual_derivative_value == Rational(-592, 2187)
            and self.finite_cusp_pair == (2, 5)
            and self.node_count == 5
            and self.infinity_pair == (5, 8)
            and self.arithmetic_genus == 21
            and self.total_delta == 21
            and self.complement_census.assignments == 40**3
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


def exact_a6_delta_five_clean_curve_certificate(
) -> A6DeltaFiveCleanCurveCertificate:
    """Build the exact clean ``(3,8)`` curve certificate."""

    clean_h = COLLISION_POLYNOMIAL.subs(CLEAN_PARAMETERS)
    parameter_p_u = FAMILY_P.subs(T, U)
    clean_q_u = CLEAN_Q.subs(T, U)
    first_collision = cancel((FAMILY_P - parameter_p_u) / (T - U))
    second_collision = cancel((CLEAN_Q - clean_q_u) / (T - U))
    return A6DeltaFiveCleanCurveCertificate(
        implicit_resultant_identity=expand(
            resultant(X - FAMILY_P, Y - CLEAN_Q, T)
            - CLEAN_BRANCH_EQUATION
        ),
        parametrized_branch_identity=expand(
            CLEAN_BRANCH_EQUATION.subs({X: FAMILY_P, Y: CLEAN_Q})
        ),
        off_diagonal_gcd=gcd(first_collision, second_collision),
        collision_resultant_identity=expand(
            resultant(first_collision, second_collision, U)
            - T**4 * CLEAN_COLLISION_FACTOR
        ),
        pairing_resultant_identity=expand(
            resultant(
                clean_h,
                T**2 - S * T + S**2 + S,
                S,
            )
            - CLEAN_COLLISION_FACTOR
        ),
        collision_discriminant=discriminant(clean_h, S),
        pair_diagonal_resultant=resultant(clean_h, -S * (3 * S + 4), S),
        pair_product_resultant=resultant(clean_h, S**2 + S, S),
        tangency_resultant=resultant(
            clean_h,
            TANGENCY_POLYNOMIAL.subs(CLEAN_PARAMETERS),
            S,
        ),
        node_coordinate_resultant_identity=expand(
            resultant(CLEAN_COLLISION_FACTOR, X - FAMILY_P, T)
            - CLEAN_NODE_COORDINATE_POLYNOMIAL**2
        ),
        node_coordinate_discriminant=discriminant(
            CLEAN_NODE_COORDINATE_POLYNOMIAL,
            X,
        ),
        cusp_image_separation=CLEAN_NODE_COORDINATE_POLYNOMIAL.subs(X, 0),
        residual_derivative_value=diff(CLEAN_Q, T).subs(T, Rational(-2, 3)),
        finite_cusp_pair=(2, 5),
        node_count=5,
        infinity_pair=(8 - 3, 8),
        arithmetic_genus=(8 - 1) * (8 - 2) // 2,
        total_delta=2 + 5 + (5 - 1) * (8 - 1) // 2,
        complement_census=three_cycle_presentation_census(CLEAN_RELATIONS),
    )


@dataclass(frozen=True, slots=True)
class A6DeltaFiveFamilyCertificate:
    """Exact algebra and representative topology for the family strata."""

    collision_reduction_remainder: Expr
    tangency_reduction_remainder: Expr
    pair_resultant_identity: Expr
    derivative_factor_identity: Expr
    tangency_resultant_identity: Expr
    cusp_collision_identity: Expr
    triple_resultant_identity: Expr
    discriminant_factor_degrees: tuple[tuple[int, int], ...]
    incidence_identities: tuple[Expr, Expr, Expr]
    clean_values: tuple[Expr, Expr, Expr, Expr]
    triple_factorization_identity: Expr
    tangency_factorization_identity: Expr
    residual_factorization_identities: tuple[Expr, Expr]
    residual_valid_products: tuple[Expr, Expr]
    residual_locus_identities: tuple[Expr, ...]
    clean: A6DeltaFiveCleanCurveCertificate
    triple_census: ThreeCyclePresentationCensus
    tangency_census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether every exact generic and codimension-one checkpoint agrees."""

        return bool(
            self.collision_reduction_remainder == 0
            and self.tangency_reduction_remainder == 0
            and self.pair_resultant_identity == 0
            and self.derivative_factor_identity == 0
            and self.tangency_resultant_identity == 0
            and self.cusp_collision_identity == 0
            and self.triple_resultant_identity == 0
            and self.discriminant_factor_degrees == ((1, 1), (7, 1))
            and self.incidence_identities == (0, 0, 0)
            and self.clean_values == (1, -37, -1, 5)
            and self.triple_factorization_identity == 0
            and self.tangency_factorization_identity == 0
            and self.residual_factorization_identities == (0, 0)
            and self.residual_valid_products == (
                Rational(768, 5),
                Rational(250880, 9),
            )
            and self.residual_locus_identities == (0, 0, 0, 0, 0, 0)
            and self.clean.verified
            and self.triple_census.assignments == 40**3
            and self.triple_census.satisfying_assignments == 40
            and self.triple_census.generated_order_histogram == ((3, 40),)
            and self.triple_census.a6_assignments == 0
            and self.tangency_census.assignments == 40**3
            and self.tangency_census.satisfying_assignments == 40
            and self.tangency_census.generated_order_histogram == ((3, 40),)
            and self.tangency_census.a6_assignments == 0
        )


def exact_a6_delta_five_family_certificate() -> A6DeltaFiveFamilyCertificate:
    """Build the exact family-stratification certificate."""

    parameter_p_u = FAMILY_P.subs(T, U)
    family_q_u = FAMILY_Q.subs(T, U)
    first_collision = cancel((FAMILY_P - parameter_p_u) / (T - U))
    second_collision = cancel((FAMILY_Q - family_q_u) / (T - U))
    collision_remainder = rem(
        second_collision - (T + U) ** 2 * COLLISION_POLYNOMIAL.subs(S, T + U),
        first_collision,
        U,
    )
    velocity_quotient = cancel(
        (
            diff(FAMILY_P, T) * diff(family_q_u, U)
            - diff(parameter_p_u, U) * diff(FAMILY_Q, T)
        )
        / (T - U)
    )
    tangency_remainder = rem(
        velocity_quotient - TANGENCY_POLYNOMIAL.subs(S, T + U),
        first_collision,
        U,
    )
    family_discriminant = discriminant(COLLISION_POLYNOMIAL, S)
    degree_seven_factor = cancel(
        family_discriminant / CUSP_COLLISION_FACTOR
    )
    factorization = factor_list(family_discriminant)
    factor_degrees = tuple(
        (int(Poly(factor, ALPHA, BETA, GAMMA).total_degree()), exponent)
        for factor, exponent in factorization[1]
    )

    Z = Symbol("z")
    reduced_q = rem(
        FAMILY_Q.subs(T, Z),
        Z**3 + Z**2 - X,
        Z,
    )
    reduced_poly = Poly(reduced_q, Z)
    quadratic_coefficient = reduced_poly.coeff_monomial(Z**2)
    linear_coefficient = reduced_poly.coeff_monomial(Z)

    triple_h = COLLISION_POLYNOMIAL.subs(TRIPLE_PARAMETERS)
    tangency_h = COLLISION_POLYNOMIAL.subs(TANGENCY_PARAMETERS)
    incidence_alpha = INCIDENCE_ALPHA_NUMERATOR / (S + 3)
    incidence_beta = INCIDENCE_BETA_NUMERATOR / (S + 3)
    incidence_substitution = {ALPHA: incidence_alpha, BETA: incidence_beta}
    incidence_first_equation = (
        (S + 3) * ALPHA - INCIDENCE_ALPHA_NUMERATOR
    )
    incidence_second_equation = (
        (S + 3) * BETA - INCIDENCE_BETA_NUMERATOR
    )
    contact_three_h = COLLISION_POLYNOMIAL.subs(CONTACT_THREE_PARAMETERS)
    mixed_h = COLLISION_POLYNOMIAL.subs(MIXED_PARAMETERS)
    contact_three_gradient = tuple(
        diff(degree_seven_factor, parameter).subs(CONTACT_THREE_PARAMETERS)
        for parameter in (ALPHA, BETA, GAMMA)
    )
    return A6DeltaFiveFamilyCertificate(
        collision_reduction_remainder=expand(collision_remainder),
        tangency_reduction_remainder=expand(tangency_remainder),
        pair_resultant_identity=expand(
            resultant(COLLISION_POLYNOMIAL, -S * (3 * S + 4), S)
            + ALPHA * EXTRA_CRITICAL_FACTOR
        ),
        derivative_factor_identity=expand(
            diff(FAMILY_Q, T).subs(T, Rational(-2, 3))
            - 16 * EXTRA_CRITICAL_FACTOR / 2187
        ),
        tangency_resultant_identity=expand(
            resultant(COLLISION_POLYNOMIAL, TANGENCY_POLYNOMIAL, S)
            - ALPHA**3 * EXTRA_CRITICAL_FACTOR * family_discriminant
        ),
        cusp_collision_identity=expand(
            COLLISION_POLYNOMIAL.subs(S, -1) - CUSP_COLLISION_FACTOR
        ),
        triple_resultant_identity=expand(
            resultant(quadratic_coefficient, linear_coefficient, X)
            + CUSP_COLLISION_FACTOR**2 * TRIPLE_COLLISION_FACTOR
        ),
        discriminant_factor_degrees=factor_degrees,
        incidence_identities=(
            cancel(COLLISION_POLYNOMIAL.subs(incidence_substitution)),
            cancel(
                diff(COLLISION_POLYNOMIAL, S).subs(incidence_substitution)
            ),
            expand(
                resultant(
                    incidence_first_equation,
                    incidence_second_equation,
                    S,
                )
                - 24 * (GAMMA - 6) * degree_seven_factor
            ),
        ),
        clean_values=(
            ALPHA.subs(CLEAN_PARAMETERS),
            EXTRA_CRITICAL_FACTOR.subs(CLEAN_PARAMETERS),
            CUSP_COLLISION_FACTOR.subs(CLEAN_PARAMETERS),
            TRIPLE_COLLISION_FACTOR.subs(CLEAN_PARAMETERS),
        ),
        triple_factorization_identity=expand(
            triple_h - (S**2 - 3) * (S**3 + 2 * S**2 + S - 1)
        ),
        tangency_factorization_identity=expand(
            tangency_h - (S + 2) ** 2 * (S**3 + S**2 + 2 * S - 1)
        ),
        residual_factorization_identities=(
            expand(
                contact_three_h
                - (S + 2) ** 3 * (5 * S**2 + 4 * S + 2) / 5
            ),
            expand(
                mixed_h
                - (S - 1) ** 2 * (3 * S + 5) * (S**2 + 3 * S + 4) / 3
            ),
        ),
        residual_valid_products=(
            expand(
                (
                    ALPHA
                    * EXTRA_CRITICAL_FACTOR
                    * CUSP_COLLISION_FACTOR
                ).subs(CONTACT_THREE_PARAMETERS)
            ),
            expand(
                (
                    ALPHA
                    * EXTRA_CRITICAL_FACTOR
                    * CUSP_COLLISION_FACTOR
                ).subs(MIXED_PARAMETERS)
            ),
        ),
        residual_locus_identities=(
            degree_seven_factor.subs(CONTACT_THREE_PARAMETERS),
            *contact_three_gradient,
            degree_seven_factor.subs(MIXED_PARAMETERS),
            TRIPLE_COLLISION_FACTOR.subs(MIXED_PARAMETERS),
        ),
        clean=exact_a6_delta_five_clean_curve_certificate(),
        triple_census=three_cycle_presentation_census(TRIPLE_RELATIONS),
        tangency_census=three_cycle_presentation_census(TANGENCY_RELATIONS),
    )


def main() -> int:
    """Print the exact family reduction and fail on any regression."""

    certificate = exact_a6_delta_five_family_certificate()
    print(
        "clean (3,8) curve:",
        {
            "cusp": certificate.clean.finite_cusp_pair,
            "nodes": certificate.clean.node_count,
            "infinity": certificate.clean.infinity_pair,
            "delta": certificate.clean.total_delta,
            "images": dict(
                certificate.clean.complement_census.generated_order_histogram
            ),
        },
    )
    print(
        "delta-five strata representation images:",
        {
            "generic": dict(
                certificate.clean.complement_census.generated_order_histogram
            ),
            "ordinary triple": dict(
                certificate.triple_census.generated_order_histogram
            ),
            "contact two": dict(
                certificate.tangency_census.generated_order_histogram
            ),
        },
    )
    print(f"conditional A6 delta-five family certificate: {certificate.verified}")
    print("excluded: generic open and both admissible codimension-one strata")
    print("remaining: codimension-at-least-two collision degenerations")
    print("claim boundary: conditional/computer-assisted; JC(2) remains open")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
