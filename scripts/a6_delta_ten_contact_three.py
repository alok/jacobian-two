"""Certify the conditional delta-ten ``C3 + 7N`` stratum.

In the normalized degree-``(4, 9)`` family, a genuine contact-three collision
is a triple root of the collision decic ``H`` on the off-diagonal pair chart.
The three equations ``H = H' = H'' = 0`` are affine-linear in the four
parameters ``(a,b,c,d)``.  Their coefficient matrix has rank three on the
valid chart except at two conjugate base points; an exact augmented minor is
nonzero there, so those rank-drop fibers are inconsistent rather than hidden
components.  The valid incidence is therefore an affine-line bundle over an
irreducible open surface.  Since ``H`` is monic in the pair-sum coordinate,
its projection to coefficient space is finite; the image closure is therefore
one three-dimensional codimension-two component.

The module also checks the rational member

``P = t^2 + t^3 + t^4`` and
``Q = 3*t^5 + 1/2*t^7 - 21/16*t^8 + t^9``.

After scaling the second target coordinate by sixteen, its affine curve has
the forced ``T(2,5)`` cusp, one two-branch contact of intersection
multiplicity three, and seven nodes.  A stored exact Sage 10.8 van Kamp
presentation simplifies to ``Z``; the dependency-free finite replay has only
forty diagonal ``C3`` images and no ``A6`` image.

This is a conditional, computer-assisted wall certificate.  It neither
classifies deeper wall intersections nor proves the plane Jacobian
conjecture.
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

CONTACT_THREE_LINEAR_PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)
CONTACT_THREE_LOCALIZATION_VARIABLE: Final = Symbol("z_contact_three")
CONTACT_THREE_VALID_LOCALIZER: Final = (
    S * PAIR_DENOMINATOR * PAIR_QUADRATIC * PAIR_DIAGONAL_FACTOR
)
CONTACT_THREE_EXPECTED_MINOR_GCD: Final = 2 * PAIR_DENOMINATOR**3 * PAIR_QUADRATIC**3
CONTACT_THREE_EXPECTED_RANKDROP_BASIS: Final = (
    5 * S**2 - 3,
    672 * CONTACT_THREE_LOCALIZATION_VARIABLE + 25,
    KAPPA + 6 * S,
)
CONTACT_THREE_RANKDROP_IDEAL: Final = (KAPPA + 6 * S, 5 * S**2 - 3)

CONTACT_THREE_CUSP_VALUE_FORM: Final = (
    ALPHA + BETA * S + GAMMA * S**2 + DELTA * S**3 + S**4
)
CONTACT_THREE_CUSP_SECOND_FORM: Final = (
    18 * ALPHA * S**4
    + 7 * ALPHA * S**2
    - ALPHA
    + 18 * BETA * S**5
    + 8 * BETA * S**3
    - 2 * BETA * S
    + 18 * GAMMA * S**6
    + 9 * GAMMA * S**4
    - 3 * GAMMA * S**2
    + 18 * DELTA * S**7
    + 10 * DELTA * S**5
    - 4 * DELTA * S**3
    + 18 * S**8
    + 11 * S**6
    - 5 * S**4
)

CONTACT_THREE_PARAMETERS: Final = {
    KAPPA: 1,
    ALPHA: 3,
    BETA: 0,
    GAMMA: Rational(1, 2),
    DELTA: Rational(-21, 16),
}
CONTACT_THREE_P: Final = expand(FAMILY_P.subs(CONTACT_THREE_PARAMETERS))
CONTACT_THREE_Q: Final = expand(FAMILY_Q.subs(CONTACT_THREE_PARAMETERS))
CONTACT_THREE_SCALED_Q: Final = expand(16 * CONTACT_THREE_Q)

CONTACT_THREE_RESIDUAL: Final = (
    16 * S**7 - 20 * S**5 - 74 * S**4 - 86 * S**3 - 33 * S**2 + 8 * S + 6
)
CONTACT_THREE_COLLISION: Final = expand((S + 2) ** 3 * CONTACT_THREE_RESIDUAL / 16)
CONTACT_THREE_TANGENCY_COFACTOR: Final = (
    72 * S**9
    + 180 * S**8
    + 100 * S**7
    - 357 * S**6
    - 920 * S**5
    - 1003 * S**4
    - 548 * S**3
    - 132 * S**2
    + 19 * S
    + 15
)
CONTACT_THREE_TANGENCY: Final = expand(
    (S + 2) ** 2 * CONTACT_THREE_TANGENCY_COFACTOR / 2
)

CONTACT_THREE_BRANCH: Final = T**2 + 2 * T + 2
CONTACT_THREE_NODE_X_POLYNOMIAL: Final = (
    16777216 * X**7
    - 76283904 * X**6
    - 5288038400 * X**5
    - 6332232000 * X**4
    - 3246794416 * X**3
    - 892022247 * X**2
    - 131282424 * X
    - 8496306
)
CONTACT_THREE_IMPLICIT: Final = (
    -65536 * X**9
    - 718463 * X**8
    - 8446608 * X**7
    + 264164 * X**6 * Y
    - 3767040 * X**6
    + 917752 * X**5 * Y
    - 628992 * X**5
    + 15430 * X**4 * Y**2
    + 397314 * X**4 * Y
    + 11048 * X**3 * Y**2
    + 65520 * X**3 * Y
    + 228 * X**2 * Y**3
    + 5288 * X**2 * Y**2
    + 16 * X * Y**3
    + 1544 * X * Y**2
    + Y**4
    + 3 * Y**3
    + 273 * Y**2
)

# Sage 10.8 exact raw affine presentation.  Generator indices 1..4 are the
# geometric meridians of a generic vertical fiber.
CONTACT_THREE_RELATIONS: Final = (
    (2, 1, -2, -1),
    (4, 3, 4, 3, 4, 3, -4, -3, -4, -3, -4, -3),
    (-4, -3, -4, -3, -4, -3, 2, 3, 4, 3, 4, 3),
    (
        -4,
        -3,
        -4,
        -3,
        -1,
        3,
        4,
        3,
        4,
        3,
        -4,
        -3,
        -4,
        -3,
        1,
        3,
        4,
        3,
        4,
        -3,
    ),
    (
        -4,
        -3,
        -4,
        -3,
        -2,
        3,
        4,
        3,
        4,
        3,
        -4,
        -3,
        -4,
        -3,
        2,
        3,
        4,
        3,
        4,
        -3,
    ),
    (2, 1, -2, -1),
    (
        -4,
        -3,
        -4,
        -3,
        1,
        2,
        1,
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
        -1,
        -2,
        -1,
        3,
        4,
        3,
    ),
    (
        -4,
        -3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        3,
        4,
        -3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        3,
    ),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, 3, 4, -3, -4, -3, 1, 2, 1),
    (-4, -3, 1, 3, 4, -3, -1, 3),
)


def _three_by_three_minors(matrix: Matrix) -> tuple[Expr, ...]:
    """Return the maximal minors of a three-row matrix."""

    if matrix.rows != 3 or matrix.cols < 3:
        msg = "maximal contact-three minors need a 3 x n matrix"
        raise ValueError(msg)
    return tuple(
        expand(matrix[:, columns].det())
        for columns in combinations(range(matrix.cols), 3)
    )


def _polynomial_gcd(polynomials: tuple[Expr, ...]) -> Expr:
    """Return the iterated gcd of a nonempty polynomial tuple."""

    if not polynomials:
        msg = "a polynomial gcd needs at least one input"
        raise ValueError(msg)
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = gcd(common, polynomial)
    return common


def _modulo_contact_branch(rational_expression: Expr) -> Expr:
    """Reduce a rational expression modulo the contact-pair quadratic."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    denominator_inverse = invert(denominator, CONTACT_THREE_BRANCH)
    return expand(
        rem(
            numerator * denominator_inverse,
            CONTACT_THREE_BRANCH,
            T,
        )
    )


@dataclass(frozen=True, slots=True)
class DeltaTenContactThreeCertificate:
    """Exact incidence, geometry, and complement data for ``C3 + 7N``."""

    incidence_minor_gcd_identity: Expr
    incidence_rankdrop_localized_basis: tuple[Expr, ...]
    incidence_augmented_inconsistency: Expr
    incidence_sample_rank_minor: Expr
    boundary_identities: tuple[Expr, ...]
    collision_identity: Expr
    tangency_identity: Expr
    collision_gcd: Expr
    triple_root_gcd: Expr
    third_derivative: Expr
    tangency_second_derivative: Expr
    residual_discriminant: Expr
    residual_contact_separation: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    denominator_resultant: Expr
    diagonal_resultant: Expr
    residual_tangency_resultant: Expr
    contact_chart_values: tuple[Expr, Expr, Expr]
    contact_pair_discriminant: Expr
    contact_image_remainders: tuple[Expr, Expr]
    contact_p_derivative_resultant: Expr
    contact_jet_values: tuple[Expr, Expr]
    contact_jet_differences: tuple[Expr, Expr, Expr]
    contact_third_jet_separation: Expr
    node_x_resultant_identity: Expr
    node_x_discriminant: Expr
    contact_node_separation: Expr
    cusp_node_separation: Expr
    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    implicit_content: Expr
    sage_jacobian_components: tuple[tuple[int, int], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    arithmetic_genus: int
    cusp_delta: int
    contact_delta: int
    node_count: int
    infinity_delta: int
    relation_count: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return (
            self.cusp_delta + self.contact_delta + self.node_count + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether every exact incidence and sample invariant agrees."""

        expected_node_discriminant = 26718085239582199662953168061403185268320248008154806376081623813523997016425064262148616945139712
        return bool(
            self.incidence_minor_gcd_identity == 0
            and self.incidence_rankdrop_localized_basis
            == CONTACT_THREE_EXPECTED_RANKDROP_BASIS
            and self.incidence_augmented_inconsistency == Rational(148635648, 625)
            and self.incidence_sample_rank_minor == -2007666
            and all(identity == 0 for identity in self.boundary_identities)
            and self.collision_identity == 0
            and self.tangency_identity == 0
            and self.collision_gcd == expand((S + 2) ** 2)
            and self.triple_root_gcd == S + 2
            and self.third_derivative == Rational(-3069, 4)
            and self.tangency_second_derivative == -9207
            and self.residual_discriminant == 168219817614587461632
            and self.residual_contact_separation == -2046
            and self.cusp_image_factor == Rational(273, 256)
            and self.extra_critical_factor == 21463
            and self.denominator_resultant == 81
            and self.diagonal_resultant == 5215509
            and self.residual_tangency_resultant == 62222061329145233117755350912
            and self.contact_chart_values == (-3, 3, 6)
            and self.contact_pair_discriminant == -4
            and self.contact_image_remainders == (-2, -464)
            and self.contact_p_derivative_resultant == 52
            and self.contact_jet_values == (448, -208)
            and self.contact_jet_differences == (0, 0, -32736 * (T + 1) / 2197)
            and self.contact_third_jet_separation == 1
            and self.node_x_resultant_identity == 0
            and self.node_x_discriminant == expected_node_discriminant
            and self.contact_node_separation == 83532198178
            and self.cusp_node_separation == -8496306
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_components == ((4, 1), (5, 1), (7, 7))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.contact_delta == 3
            and self.node_count == 7
            and self.infinity_delta == 16
            and self.total_delta == self.arithmetic_genus
            and self.relation_count == 11
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_contact_three_certificate() -> DeltaTenContactThreeCertificate:
    """Build the exact dominant-incidence and rational-member certificate."""

    derivatives = (
        COLLISION_POLYNOMIAL,
        diff(COLLISION_POLYNOMIAL, S),
        diff(COLLISION_POLYNOMIAL, S, 2),
    )
    coefficient_matrix = Matrix(
        [
            [
                diff(polynomial, parameter)
                for parameter in CONTACT_THREE_LINEAR_PARAMETERS
            ]
            for polynomial in derivatives
        ]
    )
    incidence_minors = _three_by_three_minors(coefficient_matrix)
    incidence_minor_gcd = _polynomial_gcd(incidence_minors)
    normalized_minors = tuple(
        cancel(minor / incidence_minor_gcd) for minor in incidence_minors
    )
    localized_basis = groebner(
        (
            *normalized_minors,
            1 - CONTACT_THREE_LOCALIZATION_VARIABLE * CONTACT_THREE_VALID_LOCALIZER,
        ),
        CONTACT_THREE_LOCALIZATION_VARIABLE,
        KAPPA,
        S,
        order="grevlex",
    )

    zero_parameters = dict.fromkeys(CONTACT_THREE_LINEAR_PARAMETERS, 0)
    constant_column = Matrix(
        [polynomial.subs(zero_parameters) for polynomial in derivatives]
    )
    augmented_matrix = coefficient_matrix.row_join(constant_column)
    augmented_ab_constant_minor = expand(augmented_matrix[:, [0, 1, 4]].det())
    rankdrop_basis = groebner(
        CONTACT_THREE_RANKDROP_IDEAL,
        KAPPA,
        S,
        order="lex",
    )

    sample_collision = expand(COLLISION_POLYNOMIAL.subs(CONTACT_THREE_PARAMETERS))
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(CONTACT_THREE_PARAMETERS))
    sample_coefficient_matrix = coefficient_matrix.subs(CONTACT_THREE_PARAMETERS).subs(
        S, -2
    )
    sample_rank_minor = expand(sample_coefficient_matrix[:, [0, 1, 2]].det())

    denominator_kappa = -2 * S
    cusp_kappa = -(S**2 + 1) / S
    denominator_first_form = GAMMA + 4 * DELTA * S + 11 * S**2 - 1
    denominator_second_form = (
        4 * ALPHA
        + 12 * BETA * S
        + 26 * GAMMA * S**2
        - 2 * GAMMA
        + 56 * DELTA * S**3
        - 16 * DELTA * S
        + 117 * S**4
        - 58 * S**2
        + 1
    )
    boundary_identities = (
        expand(COLLISION_POLYNOMIAL.subs(S, 0) - ALPHA * KAPPA**2),
        expand(
            COLLISION_POLYNOMIAL.subs(KAPPA, denominator_kappa) - S**2 * (S**2 - 1) ** 4
        ),
        expand(
            diff(COLLISION_POLYNOMIAL, S).subs(KAPPA, denominator_kappa)
            - 2 * S * (S**2 - 1) ** 3 * denominator_first_form
        ),
        expand(
            diff(COLLISION_POLYNOMIAL, S, 2).subs(
                KAPPA,
                denominator_kappa,
            )
            - 2 * (S**2 - 1) ** 2 * denominator_second_form
        ),
        cancel(
            S**2 * COLLISION_POLYNOMIAL.subs(KAPPA, cusp_kappa)
            - (S**2 - 1) ** 4 * CONTACT_THREE_CUSP_VALUE_FORM
        ),
        cancel(
            S**3 * diff(COLLISION_POLYNOMIAL, S).subs(KAPPA, cusp_kappa)
            - (S**2 - 1) ** 3 * (7 * S**2 + 1) * CONTACT_THREE_CUSP_VALUE_FORM
        ),
        cancel(
            S**4 * diff(COLLISION_POLYNOMIAL, S, 2).subs(KAPPA, cusp_kappa)
            - 2 * (S**2 - 1) ** 2 * CONTACT_THREE_CUSP_SECOND_FORM
        ),
        cancel(
            S**2
            - 4 * S * PAIR_QUADRATIC / PAIR_DENOMINATOR
            + S * PAIR_DIAGONAL_FACTOR / PAIR_DENOMINATOR
        ),
    )

    p_derivative = diff(CONTACT_THREE_P, T)
    x_derivatives: list[Expr] = []
    current_derivative = CONTACT_THREE_SCALED_Q
    for _ in range(3):
        current_derivative = cancel(diff(current_derivative, T) / p_derivative)
        x_derivatives.append(current_derivative)
    reduced_jets = tuple(
        _modulo_contact_branch(derivative) for derivative in x_derivatives
    )
    partner_jets = tuple(
        _modulo_contact_branch(derivative.subs(T, -2 - T))
        for derivative in x_derivatives
    )
    jet_differences = tuple(
        expand(rem(left - right, CONTACT_THREE_BRANCH, T))
        for left, right in zip(reduced_jets, partner_jets, strict=True)
    )

    sample_x_numerator = COLLISION_X_NUMERATOR.subs(CONTACT_THREE_PARAMETERS)
    sample_x_denominator = COLLISION_X_DENOMINATOR.subs(CONTACT_THREE_PARAMETERS)
    implicit_polynomial = Poly(CONTACT_THREE_IMPLICIT, X, Y)

    return DeltaTenContactThreeCertificate(
        incidence_minor_gcd_identity=expand(
            incidence_minor_gcd - CONTACT_THREE_EXPECTED_MINOR_GCD
        ),
        incidence_rankdrop_localized_basis=tuple(
            polynomial.as_expr() for polynomial in localized_basis.polys
        ),
        incidence_augmented_inconsistency=rankdrop_basis.reduce(
            augmented_ab_constant_minor
        )[1],
        incidence_sample_rank_minor=sample_rank_minor,
        boundary_identities=boundary_identities,
        collision_identity=expand(sample_collision - CONTACT_THREE_COLLISION),
        tangency_identity=expand(sample_tangency - CONTACT_THREE_TANGENCY),
        collision_gcd=gcd(
            Poly(sample_collision, S),
            Poly(sample_tangency, S),
        ).as_expr(),
        triple_root_gcd=gcd(
            gcd(
                Poly(sample_collision, S),
                Poly(diff(sample_collision, S), S),
            ),
            Poly(diff(sample_collision, S, 2), S),
        ).as_expr(),
        third_derivative=diff(sample_collision, S, 3).subs(S, -2),
        tangency_second_derivative=diff(sample_tangency, S, 2).subs(S, -2),
        residual_discriminant=discriminant(CONTACT_THREE_RESIDUAL, S),
        residual_contact_separation=CONTACT_THREE_RESIDUAL.subs(S, -2),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(CONTACT_THREE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(CONTACT_THREE_PARAMETERS),
        denominator_resultant=resultant(
            sample_collision,
            PAIR_DENOMINATOR.subs(CONTACT_THREE_PARAMETERS),
            S,
        ),
        diagonal_resultant=resultant(
            sample_collision,
            (-S * PAIR_DIAGONAL_FACTOR).subs(CONTACT_THREE_PARAMETERS),
            S,
        ),
        residual_tangency_resultant=resultant(
            CONTACT_THREE_RESIDUAL,
            sample_tangency,
            S,
        ),
        contact_chart_values=tuple(
            expression.subs(CONTACT_THREE_PARAMETERS).subs(S, -2)
            for expression in (
                PAIR_DENOMINATOR,
                PAIR_QUADRATIC,
                PAIR_DIAGONAL_FACTOR,
            )
        ),
        contact_pair_discriminant=(-2) ** 2 - 4 * 2,
        contact_image_remainders=(
            rem(CONTACT_THREE_P, CONTACT_THREE_BRANCH, T),
            rem(CONTACT_THREE_SCALED_Q, CONTACT_THREE_BRANCH, T),
        ),
        contact_p_derivative_resultant=resultant(
            CONTACT_THREE_BRANCH,
            p_derivative,
            T,
        ),
        contact_jet_values=(reduced_jets[0], reduced_jets[1]),
        contact_jet_differences=jet_differences,
        contact_third_jet_separation=resultant(
            CONTACT_THREE_BRANCH,
            T + 1,
            T,
        ),
        node_x_resultant_identity=expand(
            resultant(
                CONTACT_THREE_RESIDUAL,
                sample_x_numerator - X * sample_x_denominator,
                S,
            )
            + 9 * CONTACT_THREE_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(CONTACT_THREE_NODE_X_POLYNOMIAL, X),
        contact_node_separation=CONTACT_THREE_NODE_X_POLYNOMIAL.subs(X, -2),
        cusp_node_separation=CONTACT_THREE_NODE_X_POLYNOMIAL.subs(X, 0),
        implicit_resultant_identity=expand(
            resultant(
                CONTACT_THREE_P - X,
                CONTACT_THREE_SCALED_Q - Y,
                T,
            )
            - CONTACT_THREE_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            CONTACT_THREE_IMPLICIT.subs({X: CONTACT_THREE_P, Y: CONTACT_THREE_SCALED_Q})
        ),
        implicit_content=implicit_polynomial.content(),
        sage_jacobian_components=((4, 1), (5, 1), (7, 7)),
        sage_cyclic_simplification=(1, 0, True),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=2,
        contact_delta=3,
        node_count=7,
        infinity_delta=16,
        relation_count=len(CONTACT_THREE_RELATIONS),
        complement_census=n_generator_three_cycle_presentation_census(
            CONTACT_THREE_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact ``C3 + 7N`` certificate and fail on regression."""

    certificate = exact_delta_ten_contact_three_certificate()
    print("delta-ten C3 + 7N certificate:", certificate.verified)
    print(
        "valid rank-drop base is inconsistent:",
        certificate.incidence_augmented_inconsistency,
    )
    print(
        "sample singularity delta balance:",
        certificate.total_delta,
        "= 2 + 3 + 7 + 16",
    )
    print("sample complement census:", certificate.complement_census)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
