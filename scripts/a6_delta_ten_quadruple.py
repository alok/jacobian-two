"""Certify the ordinary-quadruple codimension-two delta-ten incidence.

The normalized conditional family is

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``.

Four distinct normalization points have one common image precisely when they
are the four roots of one fiber polynomial ``P-h`` and the remainder of ``Q``
modulo that quartic is constant.  The latter condition gives three equations
which are affine-linear in ``(a,b,c,d)``.  This module proves exactly that the
coefficient matrix has rank three on ``h != 0``.  Hence, after also removing
the critical-fiber discriminant, the incidence is an affine-line bundle over
an irreducible open of the ``(k,h)`` plane.

An exact rational member is then checked to have one ordinary quadruple fiber
and four residual nodes, in addition to the forced ``T(2,5)`` cusp and the
``T(5,9)`` branch at infinity.  These are incidence and local-geometry
certificates.  Excluding the component as an ``A6`` branch also requires the
separate computer-assisted complement presentation and Whitney--Thom
transport; no plane Jacobian-conjecture claim is made here.
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
    groebner,
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

FIBER_VALUE: Final = Symbol("h")
SLOPE_VALUE: Final = Symbol("m")
TARGET_X: Final = Symbol("x_quad")
LOCALIZATION_VARIABLE: Final = Symbol("z_quad")

COEFFICIENTS: Final = (ALPHA, BETA, GAMMA, DELTA)
FIBER_POLYNOMIAL: Final = expand(FAMILY_P - FIBER_VALUE)
Q_FIBER_REMAINDER: Final = expand(rem(FAMILY_Q, FIBER_POLYNOMIAL, T))
NONCONSTANT_REMAINDERS: Final = tuple(
    Poly(Q_FIBER_REMAINDER, T).coeff_monomial(T**degree) for degree in (1, 2, 3)
)
COEFFICIENT_MATRIX: Final = Matrix(
    [
        [diff(equation, coefficient) for coefficient in COEFFICIENTS]
        for equation in NONCONSTANT_REMAINDERS
    ]
)
COEFFICIENT_MINORS: Final = tuple(
    expand(COEFFICIENT_MATRIX[:, columns].det())
    for columns in combinations(range(4), 3)
)
EXPECTED_COEFFICIENT_MINORS: Final = tuple(
    expand(polynomial)
    for polynomial in (
        FIBER_VALUE**2 * (FIBER_VALUE + 1),
        -(FIBER_VALUE**2) * KAPPA * (2 * FIBER_VALUE + 1),
        FIBER_VALUE**2 * (FIBER_VALUE * KAPPA**2 + 2 * FIBER_VALUE + 1),
        -2 * FIBER_VALUE**3 * KAPPA,
    )
)

# This chart solves the first three coefficient columns.  It is enough to
# exhibit the sample and show that projection to parameter space has the full
# expected three-dimensional image.  Rank on the entire h != 0 open is checked
# independently by all four minors and a localized Groebner certificate.
GENERIC_ALPHA: Final = cancel(
    FIBER_VALUE
    * (2 * DELTA * KAPPA - FIBER_VALUE - 2 * KAPPA**2 - 1)
    / (FIBER_VALUE + 1)
)
GENERIC_BETA: Final = cancel(
    (DELTA - KAPPA) * (FIBER_VALUE * KAPPA**2 + 2 * FIBER_VALUE + 1) / (FIBER_VALUE + 1)
)
GENERIC_GAMMA: Final = cancel(
    (
        2 * DELTA * FIBER_VALUE * KAPPA
        + DELTA * KAPPA
        - 2 * FIBER_VALUE * KAPPA**2
        + FIBER_VALUE
        - KAPPA**2
        + 1
    )
    / (FIBER_VALUE + 1)
)
GENERIC_SOLUTION: Final = {
    ALPHA: GENERIC_ALPHA,
    BETA: GENERIC_BETA,
    GAMMA: GENERIC_GAMMA,
}

SAMPLE_FIBER_VALUE: Final = Rational(1)
SAMPLE_PARAMETERS: Final = {
    KAPPA: 1,
    ALPHA: -2,
    BETA: -2,
    GAMMA: Rational(-1, 2),
    DELTA: 0,
}
SAMPLE_P: Final = expand(FAMILY_P.subs(SAMPLE_PARAMETERS))
SAMPLE_Q: Final = expand(FAMILY_Q.subs(SAMPLE_PARAMETERS))
SAMPLE_FIBER_POLYNOMIAL: Final = expand(SAMPLE_P - SAMPLE_FIBER_VALUE)
SAMPLE_Q_REMAINDER: Final = expand(rem(SAMPLE_Q, SAMPLE_FIBER_POLYNOMIAL, T))
SAMPLE_FIBER_QUOTIENT: Final = cancel(
    (SAMPLE_Q - SAMPLE_Q_REMAINDER) / SAMPLE_FIBER_POLYNOMIAL
)

QUADRUPLE_PAIR_FACTOR_LEFT: Final = S**3 + S + 1
QUADRUPLE_PAIR_FACTOR_RIGHT: Final = S**3 + 3 * S**2 + 4 * S + 1
QUADRUPLE_DIAGONAL_SUM_FACTOR: Final = (S + 2) * (S**3 + 4 * S - 8)
RESIDUAL_NODE_SUM_FACTOR: Final = 2 * S**4 + 6 * S**3 + 8 * S**2 - 3 * S - 4
SAMPLE_COLLISION_POLYNOMIAL: Final = expand(
    COLLISION_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
)
SAMPLE_PAIR_SUM_RESULTANT: Final = expand(
    resultant(
        SAMPLE_FIBER_POLYNOMIAL,
        SAMPLE_FIBER_POLYNOMIAL.subs(T, S - T),
        T,
    )
)

SAMPLE_SLOPE_RESULTANT: Final = resultant(
    SAMPLE_FIBER_POLYNOMIAL,
    SLOPE_VALUE - SAMPLE_FIBER_QUOTIENT,
    T,
)
EXPECTED_SLOPE_RESULTANT: Final = (
    (2 * SLOPE_VALUE + 5)
    * (8 * SLOPE_VALUE**3 + 72 * SLOPE_VALUE**2 + 194 * SLOPE_VALUE + 127)
    / 16
)

RESIDUAL_NODE_X_POLYNOMIAL: Final = (
    64 * TARGET_X**4 - 336 * TARGET_X**3 + 484 * TARGET_X**2 + 189 * TARGET_X - 588
)
RESIDUAL_NODE_X_ELIMINANT: Final = resultant(
    RESIDUAL_NODE_SUM_FACTOR,
    TARGET_X * COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    - COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS),
    S,
)

# Sage 10.8 regenerates this unsimplified affine Zariski--van Kamp
# presentation from the integral target rescaling ``(P, 2*Q)``.  Signed
# integers name the four geometric meridians and record inversion.
QUADRUPLE_RAW_RELATIONS: Final = (
    (-4, -3, -2, -1, 2, 3, 4, -3, 2, 3),
    (2, 1, -2, -1),
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
        -4,
        -3,
        -1,
        3,
        -4,
        -3,
        -1,
        3,
    ),
    (2, 3, 4, 1, -3, -2, -1, -4),
    (1, -2, -1, -4, -3, 2, 3, 4),
    (1, -4, -3, -2, -1, 2, 3, 4),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-4, -3, 1, 3, 4, -3, -1, 3),
    (-4, -3, -1, -3, 1, 3, 4, -3, 1, 3),
    (-3, 1, 3, -1),
)


def _polynomial_gcd(polynomials: tuple[Expr, ...]) -> Expr:
    """Return the normalized iterated polynomial gcd."""

    if not polynomials:
        msg = "a polynomial gcd needs at least one input"
        raise ValueError(msg)
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = gcd(common, polynomial)
    return common


@dataclass(frozen=True, slots=True)
class DeltaTenQuadrupleCertificate:
    """Exact rank, component, and local-geometry checks for ``Q0 + 4N``."""

    coefficient_minors: tuple[Expr, ...]
    coefficient_minor_gcd: Expr
    localized_minor_basis: tuple[Expr, ...]
    generic_solution_residuals: tuple[Expr, ...]
    sample_projection_derivatives: tuple[Expr, Expr, Expr]
    sample_remainder: Expr
    sample_fiber_discriminant: Expr
    sample_quotient_identity: Expr
    sample_slope_resultant_identity: Expr
    sample_slope_discriminant: Expr
    pair_sum_resultant_identity: Expr
    complementary_pair_identity: Expr
    collision_factor_identity: Expr
    collision_factor_discriminants: tuple[Expr, Expr, Expr]
    collision_factor_resultants: tuple[Expr, Expr, Expr]
    validity_resultants: tuple[tuple[Expr, Expr, Expr, Expr], ...]
    node_x_eliminant_identity: Expr
    node_x_discriminant: Expr
    node_x_at_quadruple: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    affine_delta_total: int
    projective_delta_total: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether every exact ordinary-quadruple check agrees."""

        return bool(
            self.coefficient_minors == EXPECTED_COEFFICIENT_MINORS
            and self.coefficient_minor_gcd == FIBER_VALUE**2
            and self.localized_minor_basis == (1,)
            and self.generic_solution_residuals == (0, 0, 0)
            and self.sample_projection_derivatives
            == (Rational(-3, 2), Rational(-1, 2), Rational(-1, 4))
            and self.sample_remainder == Rational(-1, 2)
            and self.sample_fiber_discriminant == -279
            and self.sample_quotient_identity == 0
            and self.sample_slope_resultant_identity == 0
            and self.sample_slope_discriminant == Rational(-9756351, 4096)
            and self.pair_sum_resultant_identity == 0
            and self.complementary_pair_identity == 0
            and self.collision_factor_identity == 0
            and self.collision_factor_discriminants == (-31, -31, -458892)
            and self.collision_factor_resultants == (-27, 153, -99)
            and all(
                all(value != 0 for value in row) for row in self.validity_resultants
            )
            and self.node_x_eliminant_identity == 0
            and self.node_x_discriminant == -265262438203392
            and self.node_x_at_quadruple == -187
            and self.cusp_image_factor == Rational(7, 4)
            and self.extra_critical_factor == 13392
            and self.affine_delta_total == 12
            and self.projective_delta_total == 28
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_quadruple_certificate() -> DeltaTenQuadrupleCertificate:
    """Build the dependency-free exact certificate for the ``Q0`` profile."""

    normalized_minors = tuple(
        cancel(minor / FIBER_VALUE**2) for minor in COEFFICIENT_MINORS
    )
    localized_basis = groebner(
        (
            *normalized_minors,
            1 - LOCALIZATION_VARIABLE * FIBER_VALUE,
        ),
        LOCALIZATION_VARIABLE,
        KAPPA,
        FIBER_VALUE,
        order="lex",
    )

    sample_solution_point = {
        KAPPA: 1,
        FIBER_VALUE: 1,
        DELTA: 0,
    }
    collision_factors = (
        QUADRUPLE_PAIR_FACTOR_LEFT,
        QUADRUPLE_PAIR_FACTOR_RIGHT,
        RESIDUAL_NODE_SUM_FACTOR,
    )
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(SAMPLE_PARAMETERS))
    sample_denominator = expand(PAIR_DENOMINATOR.subs(SAMPLE_PARAMETERS))
    sample_diagonal = expand(PAIR_DIAGONAL_FACTOR.subs(SAMPLE_PARAMETERS))
    sample_cusp_pair = expand(PAIR_QUADRATIC.subs(SAMPLE_PARAMETERS))

    return DeltaTenQuadrupleCertificate(
        coefficient_minors=tuple(expand(minor) for minor in COEFFICIENT_MINORS),
        coefficient_minor_gcd=_polynomial_gcd(COEFFICIENT_MINORS),
        localized_minor_basis=tuple(
            polynomial.as_expr() for polynomial in localized_basis.polys
        ),
        generic_solution_residuals=tuple(
            cancel(equation.subs(GENERIC_SOLUTION))
            for equation in NONCONSTANT_REMAINDERS
        ),
        sample_projection_derivatives=tuple(
            cancel(diff(expression, FIBER_VALUE).subs(sample_solution_point))
            for expression in (GENERIC_ALPHA, GENERIC_BETA, GENERIC_GAMMA)
        ),
        sample_remainder=SAMPLE_Q_REMAINDER,
        sample_fiber_discriminant=discriminant(SAMPLE_FIBER_POLYNOMIAL, T),
        sample_quotient_identity=expand(
            SAMPLE_Q
            - SAMPLE_Q_REMAINDER
            - SAMPLE_FIBER_POLYNOMIAL * SAMPLE_FIBER_QUOTIENT
        ),
        sample_slope_resultant_identity=expand(
            SAMPLE_SLOPE_RESULTANT - EXPECTED_SLOPE_RESULTANT
        ),
        sample_slope_discriminant=discriminant(EXPECTED_SLOPE_RESULTANT, SLOPE_VALUE),
        pair_sum_resultant_identity=expand(
            SAMPLE_PAIR_SUM_RESULTANT
            - QUADRUPLE_DIAGONAL_SUM_FACTOR
            * QUADRUPLE_PAIR_FACTOR_LEFT**2
            * QUADRUPLE_PAIR_FACTOR_RIGHT**2
        ),
        complementary_pair_identity=expand(
            QUADRUPLE_PAIR_FACTOR_RIGHT + QUADRUPLE_PAIR_FACTOR_LEFT.subs(S, -1 - S)
        ),
        collision_factor_identity=expand(
            2 * SAMPLE_COLLISION_POLYNOMIAL
            - QUADRUPLE_PAIR_FACTOR_LEFT
            * QUADRUPLE_PAIR_FACTOR_RIGHT
            * RESIDUAL_NODE_SUM_FACTOR
        ),
        collision_factor_discriminants=tuple(
            discriminant(factor, S) for factor in collision_factors
        ),
        collision_factor_resultants=(
            resultant(
                QUADRUPLE_PAIR_FACTOR_LEFT,
                QUADRUPLE_PAIR_FACTOR_RIGHT,
                S,
            ),
            resultant(
                QUADRUPLE_PAIR_FACTOR_LEFT,
                RESIDUAL_NODE_SUM_FACTOR,
                S,
            ),
            resultant(
                QUADRUPLE_PAIR_FACTOR_RIGHT,
                RESIDUAL_NODE_SUM_FACTOR,
                S,
            ),
        ),
        validity_resultants=tuple(
            (
                resultant(factor, sample_denominator, S),
                resultant(factor, sample_diagonal, S),
                resultant(factor, sample_cusp_pair, S),
                resultant(factor, sample_tangency, S),
            )
            for factor in collision_factors
        ),
        node_x_eliminant_identity=expand(
            RESIDUAL_NODE_X_ELIMINANT - 81 * RESIDUAL_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(RESIDUAL_NODE_X_POLYNOMIAL, TARGET_X),
        node_x_at_quadruple=RESIDUAL_NODE_X_POLYNOMIAL.subs(
            TARGET_X, SAMPLE_FIBER_VALUE
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(SAMPLE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(SAMPLE_PARAMETERS),
        affine_delta_total=2 + 6 + 4,
        projective_delta_total=2 + 6 + 4 + 16,
        complement_census=n_generator_three_cycle_presentation_census(
            QUADRUPLE_RAW_RELATIONS,
            4,
        ),
    )


def main() -> None:
    """Print the exact ordinary-quadruple certificate summary."""

    certificate = exact_delta_ten_quadruple_certificate()
    print("delta-ten ordinary-quadruple certificate:", certificate.verified)
    print("coefficient minors:", certificate.coefficient_minors)
    print("sample Q remainder:", certificate.sample_remainder)
    print("sample fiber discriminant:", certificate.sample_fiber_discriminant)
    print(
        "sample collision factor discriminants:",
        certificate.collision_factor_discriminants,
    )
    print("residual node target discriminant:", certificate.node_x_discriminant)
    print("raw presentation census:", certificate.complement_census)


if __name__ == "__main__":
    main()
