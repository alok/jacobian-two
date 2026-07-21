"""Certify one exact generic contact-two member at conditional delta ten.

The clean ``(4, 9)`` family has ten simple normalization collisions.  This
module moves to the generic repeated-collision divisor and checks the rational
member

``P = t^2 + t^3 + t^4`` and
``Q = -12/5*t^5 - 16/5*t^6 + t^9``.

Its collision decic has one double root and eight simple residual roots.  The
double root gives two smooth branches with common tangent and distinct
curvatures, while the residual collisions are transverse and have distinct
target values.  Thus the affine curve has the forced ``T(2,5)`` cusp, one
contact-two point, and eight nodes.  Together with the ``T(5,9)`` branch at
infinity these singularities exhaust the projective genus budget.

The module also stores the exact twelve-relator Sage 10.8 affine van Kamp
presentation and independently replays all ``40^4`` assignments of its four
geometric meridians to single three-cycles.  Only the forty diagonal ``C3``
assignments survive, so this member has no required ``A6`` quotient.

This is a conditional, computer-assisted certificate for one exact member of
the generic contact-two locus.  Propagation across the whole wall still needs
the appropriate connected-stratum and proper Whitney--Thom argument.  Nothing
here proves or disproves the two-dimensional Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final, TypeAlias

from sympy import (
    Expr,
    Poly,
    Rational,
    Symbol,
    diff,
    discriminant,
    expand,
    rem,
    resultant,
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
    PAIR_DIAGONAL_FACTOR,
    PAIR_INCIDENCE,
    R,
    S,
    SECOND_DIVIDED_DIFFERENCE,
    T,
    TANGENCY_POLYNOMIAL,
    X,
    Y,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    alternating_group_six,
    evaluate_signed_word,
    n_generator_three_cycle_presentation_census,
)
from scripts.six_sheet_monodromy import (
    IDENTITY,
    Permutation,
    cycle_type,
)

Assignment: TypeAlias = tuple[Permutation, ...]

LOCAL_X: Final = Symbol("x")
LOCAL_Y: Final = Symbol("y")

CONTACT_PARAMETERS: Final = {
    KAPPA: 1,
    ALPHA: Rational(-12, 5),
    BETA: Rational(-16, 5),
    GAMMA: 0,
    DELTA: 0,
}
CONTACT_P: Final = expand(FAMILY_P.subs(CONTACT_PARAMETERS))
CONTACT_Q: Final = expand(FAMILY_Q.subs(CONTACT_PARAMETERS))
# Scaling the second target coordinate by five preserves the curve complement
# and makes the implicit equation primitive over the integers.
CONTACT_SCALED_Q: Final = expand(5 * CONTACT_Q)

CONTACT_RESIDUAL_COLLISION: Final = (
    5 * S**8
    + 10 * S**7
    + 25 * S**6
    + 29 * S**5
    + 54 * S**4
    + 27 * S**3
    - 21 * S**2
    - 18 * S
    - 3
)
CONTACT_COLLISION_POLYNOMIAL: Final = expand(
    (S + 2) ** 2 * CONTACT_RESIDUAL_COLLISION / 5
)

CONTACT_TANGENCY_COFACTOR: Final = (
    60 * S**10
    + 270 * S**9
    + 675 * S**8
    + 1157 * S**7
    + 1561 * S**6
    + 1684 * S**5
    + 1161 * S**4
    + 193 * S**3
    - 271 * S**2
    - 152 * S
    - 20
)

CONTACT_NODE_X_POLYNOMIAL: Final = (
    15625 * X**8
    - 87500 * X**7
    + 332500 * X**6
    - 747150 * X**5
    + 1299030 * X**4
    - 493779 * X**3
    - 58605 * X**2
    + 1202586 * X
    - 636804
)

CONTACT_IMPLICIT: Final = (
    -625 * X**9
    + 12000 * X**8
    - 42600 * X**7
    + 4225 * X**6 * Y
    + 66496 * X**6
    - 14715 * X**5 * Y
    - 19152 * X**5
    + 675 * X**4 * Y**2
    + 20496 * X**4 * Y
    - 1797 * X**3 * Y**2
    - 3724 * X**3 * Y
    + 45 * X**2 * Y**3
    + 1410 * X**2 * Y**2
    - 60 * X * Y**3
    - 48 * X * Y**2
    + Y**4
    + 10 * Y**3
    + 133 * Y**2
)

CONTACT_BRANCH_POLYNOMIAL: Final = T**2 + 2 * T + 2
CONTACT_WEIGHTED_INITIAL: Final = (
    10309 * LOCAL_Y**2
    + 379176 * LOCAL_X**2 * LOCAL_Y
    + 3601440 * LOCAL_X**4
)

# Sage 10.8 exact raw affine presentation.  Generator indices 1..4 are the
# geometric meridians in a generic vertical fiber.  The unsimplified list is
# retained so that the finite-image replay does not trust Tietze reduction.
CONTACT_RELATIONS: Final = (
    (3, 2, -3, -2),
    (4, 3, -4, -3),
    (
        -4,
        -3,
        -2,
        -1,
        2,
        3,
        4,
        3,
        -4,
        -3,
        -2,
        1,
        2,
        3,
        4,
        3,
        -4,
        -3,
        -2,
        1,
        2,
        3,
        4,
        -3,
        -4,
        -3,
        -2,
        -1,
        2,
        3,
        4,
        -3,
    ),
    (-4, -3, 2, 3, 4, -3, -2, 3),
    (2, 1, -2, -1),
    (-4, -3, -2, -1, 3, -4, -3, 1, 2, 3, 4, 3),
    (
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        -3,
        1,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
    ),
    (-2, -1, 3, -4, -3, 1, 3, 4, -3, 1),
    (-4, -3, -1, 3, 4, 3, -4, -3, 1, 3, 4, -3),
    (-4, -3, 1, 3, 4, -3, -1, 3),
    (-4, -3, -2, 3, 4, 3, -4, -3, 2, 3, 4, -3),
    (-4, -3, 2, 3, 4, -3, -2, 3),
)


def _weighted_truncation(polynomial: Expr, maximum_weight: int) -> Expr:
    """Return terms of weight at most ``maximum_weight`` for weights (1,2)."""

    result: Expr = 0
    for monomial, coefficient in Poly(polynomial, LOCAL_X, LOCAL_Y).terms():
        x_degree, y_degree = monomial
        if x_degree + 2 * y_degree <= maximum_weight:
            result += coefficient * LOCAL_X**x_degree * LOCAL_Y**y_degree
    return expand(result)


@cache
def _presentation_pruning() -> tuple[tuple[int, ...], tuple[Assignment, ...]]:
    """Replay incremental relation pruning and retain the final assignments."""

    three_cycles = tuple(
        element
        for element in alternating_group_six()
        if cycle_type(element) == (3, 1, 1, 1)
    )
    partial_assignments: tuple[Assignment, ...] = ((),)
    stage_counts: list[int] = []
    for assigned_count in range(1, 5):
        newly_decidable = tuple(
            relation
            for relation in CONTACT_RELATIONS
            if relation
            and max(abs(letter) for letter in relation) == assigned_count
        )
        partial_assignments = tuple(
            partial + (image,)
            for partial in partial_assignments
            for image in three_cycles
            if all(
                evaluate_signed_word(relation, partial + (image,)) == IDENTITY
                for relation in newly_decidable
            )
        )
        stage_counts.append(len(partial_assignments))
    return tuple(stage_counts), partial_assignments


def _curvature_difference_identity() -> Expr:
    """Clear denominators in the exact difference of branch curvatures."""

    x_prime = diff(CONTACT_P, T)
    y_prime = diff(CONTACT_SCALED_Q, T)
    numerator = diff(CONTACT_SCALED_Q, T, 2) * x_prime - y_prime * diff(
        CONTACT_P, T, 2
    )
    denominator = x_prime**3
    partner_numerator = expand(numerator.subs(T, -2 - T))
    partner_denominator = expand(denominator.subs(T, -2 - T))
    expected_difference = Rational(2256, 169) * (T + 1)
    return rem(
        expand(
            numerator * partner_denominator
            - partner_numerator * denominator
            - expected_difference * denominator * partner_denominator
        ),
        CONTACT_BRANCH_POLYNOMIAL,
        T,
    )


@dataclass(frozen=True, slots=True)
class DeltaTenContactWallCertificate:
    """Exact geometry and finite-image audit for the contact-two member."""

    collision_factor_identity: Expr
    residual_collision_discriminant: Expr
    contact_residual_separation: Expr
    denominator_value: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    pair_diagonal_resultant: Expr
    tangency_factor_identity: Expr
    contact_tangency_cofactor_value: Expr
    residual_tangency_resultant: Expr
    contact_incidence_identities: tuple[Expr, Expr]
    contact_pair_discriminant: Expr
    contact_image_remainders: tuple[Expr, Expr]
    contact_slope_remainder: Expr
    curvature_denominator_resultants: tuple[Expr, Expr]
    curvature_difference_identity: Expr
    weighted_initial_identity: Expr
    weighted_initial_discriminant: Expr
    node_x_resultant_identity: Expr
    node_x_discriminant: Expr
    contact_node_separation: Expr
    cusp_node_separation: Expr
    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    implicit_content: Expr
    sage_jacobian_component_lengths: tuple[tuple[int, int], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    arithmetic_genus: int
    cusp_delta: int
    contact_delta: int
    node_count: int
    infinity_delta: int
    relation_count: int
    pruning_stage_survivors: tuple[int, ...]
    diagonal_satisfying_assignments: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return (
            self.cusp_delta
            + self.contact_delta
            + self.node_count
            + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether every exact sample and presentation invariant agrees."""

        expected_node_discriminant = -(
            2**2
            * 3**8
            * 5**23
            * 7**5
            * 19**5
            * 29
            * 61**16
            * 79**2
            * 787**2
            * 37215261907
        )
        return bool(
            self.collision_factor_identity == 0
            and self.residual_collision_discriminant
            == -52320062156396485500
            and self.contact_residual_separation == 1269
            and self.denominator_value == Rational(81, 1024)
            and self.cusp_image_factor == Rational(133, 25)
            and self.extra_critical_factor == Rational(852228, 25)
            and self.pair_diagonal_resultant == Rational(-828365616, 125)
            and self.tangency_factor_identity == 0
            and self.contact_tangency_cofactor_value == 10152
            and self.residual_tangency_resultant
            == Rational(29153413543712129429066042292, 6103515625)
            and self.contact_incidence_identities == (0, 0)
            and self.contact_pair_discriminant == -4
            and self.contact_image_remainders == (0, 0)
            and self.contact_slope_remainder == 0
            and self.curvature_denominator_resultants == (140608, 140608)
            and self.curvature_difference_identity == 0
            and self.weighted_initial_identity == 0
            and self.weighted_initial_discriminant == -4734540864
            and self.node_x_resultant_identity == 0
            and self.node_x_discriminant == expected_node_discriminant
            and self.contact_node_separation == 81847116
            and self.cusp_node_separation == -636804
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_component_lengths
            == ((4, 1), (3, 1), (8, 8))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.contact_delta == 2
            and self.node_count == 8
            and self.infinity_delta == 16
            and self.total_delta == self.arithmetic_genus
            and self.relation_count == 12
            and self.pruning_stage_survivors == (40, 160, 640, 40)
            and self.diagonal_satisfying_assignments == 40
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_contact_wall_certificate() -> DeltaTenContactWallCertificate:
    """Build the exact repeated-collision member and presentation certificate."""

    sample_collision = expand(COLLISION_POLYNOMIAL.subs(CONTACT_PARAMETERS))
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(CONTACT_PARAMETERS))
    sample_x_numerator = expand(
        COLLISION_X_NUMERATOR.subs(CONTACT_PARAMETERS)
    )
    sample_x_denominator = expand(
        COLLISION_X_DENOMINATOR.subs(CONTACT_PARAMETERS)
    )
    x_prime_cubed = diff(CONTACT_P, T) ** 3
    partner_x_prime_cubed = expand(x_prime_cubed.subs(T, -2 - T))
    translated_implicit = expand(
        CONTACT_IMPLICIT.subs(
            {
                X: LOCAL_X - 2,
                Y: LOCAL_Y - 128 + 96 * LOCAL_X,
            }
        )
    )
    stage_survivors, assignments = _presentation_pruning()
    diagonal_assignments = sum(
        all(image == images[0] for image in images[1:])
        for images in assignments
    )
    implicit_polynomial = Poly(CONTACT_IMPLICIT, X, Y)

    return DeltaTenContactWallCertificate(
        collision_factor_identity=expand(
            sample_collision - CONTACT_COLLISION_POLYNOMIAL
        ),
        residual_collision_discriminant=discriminant(
            CONTACT_RESIDUAL_COLLISION,
            S,
        ),
        contact_residual_separation=CONTACT_RESIDUAL_COLLISION.subs(S, -2),
        denominator_value=sample_collision.subs(S, Rational(-1, 2)),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(CONTACT_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(CONTACT_PARAMETERS),
        pair_diagonal_resultant=resultant(
            sample_collision,
            (-S * PAIR_DIAGONAL_FACTOR).subs(CONTACT_PARAMETERS),
            S,
        ),
        tangency_factor_identity=expand(
            sample_tangency
            - Rational(3, 5) * (S + 2) * CONTACT_TANGENCY_COFACTOR
        ),
        contact_tangency_cofactor_value=CONTACT_TANGENCY_COFACTOR.subs(S, -2),
        residual_tangency_resultant=resultant(
            CONTACT_RESIDUAL_COLLISION / 5,
            sample_tangency,
            S,
        ),
        contact_incidence_identities=(
            PAIR_INCIDENCE.subs(CONTACT_PARAMETERS).subs({S: -2, R: 2}),
            SECOND_DIVIDED_DIFFERENCE.subs(CONTACT_PARAMETERS).subs(
                {S: -2, R: 2}
            ),
        ),
        contact_pair_discriminant=(-2) ** 2 - 4 * 2,
        contact_image_remainders=(
            rem(CONTACT_P + 2, CONTACT_BRANCH_POLYNOMIAL, T),
            rem(CONTACT_SCALED_Q + 128, CONTACT_BRANCH_POLYNOMIAL, T),
        ),
        contact_slope_remainder=rem(
            diff(CONTACT_SCALED_Q, T) - 96 * diff(CONTACT_P, T),
            CONTACT_BRANCH_POLYNOMIAL,
            T,
        ),
        curvature_denominator_resultants=(
            resultant(CONTACT_BRANCH_POLYNOMIAL, x_prime_cubed, T),
            resultant(
                CONTACT_BRANCH_POLYNOMIAL,
                partner_x_prime_cubed,
                T,
            ),
        ),
        curvature_difference_identity=_curvature_difference_identity(),
        weighted_initial_identity=expand(
            _weighted_truncation(translated_implicit, 4)
            - CONTACT_WEIGHTED_INITIAL
        ),
        weighted_initial_discriminant=(
            379176**2 - 4 * 10309 * 3601440
        ),
        node_x_resultant_identity=expand(
            resultant(
                sample_collision,
                sample_x_numerator - X * sample_x_denominator,
                S,
            )
            - Rational(6561, 15625)
            * (X + 2) ** 2
            * CONTACT_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(CONTACT_NODE_X_POLYNOMIAL, X),
        contact_node_separation=CONTACT_NODE_X_POLYNOMIAL.subs(X, -2),
        cusp_node_separation=CONTACT_NODE_X_POLYNOMIAL.subs(X, 0),
        implicit_resultant_identity=expand(
            resultant(CONTACT_P - X, CONTACT_SCALED_Q - Y, T)
            - CONTACT_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            CONTACT_IMPLICIT.subs({X: CONTACT_P, Y: CONTACT_SCALED_Q})
        ),
        implicit_content=implicit_polynomial.content(),
        # Exact Sage primary-component (length, radical-degree) output:
        # forced cusp, contact-two point, and eight reduced nodes.
        sage_jacobian_component_lengths=((4, 1), (3, 1), (8, 8)),
        # Simplified presentation: one generator, no relators; all four raw
        # geometric meridians map to that generator.
        sage_cyclic_simplification=(1, 0, True),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=(2 - 1) * (5 - 1) // 2,
        contact_delta=2,
        node_count=8,
        infinity_delta=(5 - 1) * (9 - 1) // 2,
        relation_count=len(CONTACT_RELATIONS),
        pruning_stage_survivors=stage_survivors,
        diagonal_satisfying_assignments=diagonal_assignments,
        complement_census=n_generator_three_cycle_presentation_census(
            CONTACT_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact contact-wall certificate and fail on regression."""

    certificate = exact_delta_ten_contact_wall_certificate()
    print(
        "delta-ten contact member:",
        {
            "cusp": (2, 5),
            "contact orders": (2,),
            "nodes": certificate.node_count,
            "infinity": (5, 9),
            "delta": certificate.total_delta,
        },
    )
    print("presentation pruning:", certificate.pruning_stage_survivors)
    print(
        "single-three-cycle images:",
        dict(certificate.complement_census.generated_order_histogram),
    )
    print("exact contact-two member excluded from A6:", certificate.verified)
    print(
        "claim boundary: exact conditional/computer-assisted member only; "
        "wall propagation and JC(2) remain open"
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
