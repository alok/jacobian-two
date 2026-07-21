"""Certify one exact ordinary-triple member at conditional delta ten.

The normalized ``(4, 9)`` family has a legitimate split collision incidence
when ``k = 2``.  This module works directly on that exceptional fiber and
checks the rational member

``P = t^2 + 2*t^3 + t^4`` and
``Q = 114/625*t^5 - 1/5*t^6 - 3*t^7 - 3*t^8 + t^9``.

Three distinct normalization points have one common image and distinct
tangent lines.  The residual split-incidence factors give seven ordinary
nodes.  Together with the forced ``T(2,5)`` cusp and the ``T(5,9)`` branch at
infinity, these singularities exhaust the projective genus budget.

The module also stores the exact twelve-relator Sage 10.8 affine van Kamp
presentation and independently replays all ``40^4`` assignments of its four
geometric meridians to single three-cycles.  Only the forty diagonal ``C3``
assignments survive, so this member has no required ``A6`` quotient.

This is a conditional, computer-assisted certificate for one exact member of
the ordinary-triple locus.  Propagation across the whole wall still requires
the appropriate connected-stratum and proper Whitney--Thom argument.  It does
not classify deeper wall intersections and does not prove or disprove the
two-dimensional Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final, TypeAlias

from sympy import Expr, Poly, Rational, Symbol, diff, discriminant, expand, resultant

from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    CUSP_IMAGE_FACTOR,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    R,
    S,
    T,
    X,
    Y,
    exceptional_graph,
    exceptional_vertical,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    alternating_group_six,
    evaluate_signed_word,
    n_generator_three_cycle_presentation_census,
)
from scripts.six_sheet_monodromy import IDENTITY, Permutation, cycle_type

Assignment: TypeAlias = tuple[Permutation, ...]

LOCAL_U: Final = Symbol("u")
LOCAL_V: Final = Symbol("v")

TRIPLE_WALL_PARAMETERS: Final = {
    KAPPA: 2,
    ALPHA: Rational(114, 625),
    BETA: Rational(-1, 5),
    GAMMA: -3,
    DELTA: -3,
}
TRIPLE_WALL_P: Final = expand(FAMILY_P.subs(TRIPLE_WALL_PARAMETERS))
TRIPLE_WALL_Q: Final = expand(FAMILY_Q.subs(TRIPLE_WALL_PARAMETERS))

TRIPLE_PARAMETERS: Final = (
    Rational(-3, 5),
    Rational(-2, 5),
    Rational(1, 5),
)
FOURTH_FIBER_PARAMETER: Final = Rational(-6, 5)
TRIPLE_IMAGE_X: Final = Rational(36, 625)
TRIPLE_IMAGE_Y: Final = Rational(0)

# For k=2 the unordered-pair incidence is reducible.  The vertical component
# has s=-1 and coordinate r=tu; the graph component has
# r=s(s+1)/2.  The three linear factors below are exactly the three pairs at
# the ordinary triple point, while the residual cubic and quartic encode the
# seven remaining normalization collisions.
TRIPLE_VERTICAL_RESIDUAL: Final = 25 * R**3 - 469 * R**2 + 582 * R - 144
TRIPLE_VERTICAL_COLLISIONS: Final = expand(
    (25 * R - 6) * TRIPLE_VERTICAL_RESIDUAL / 625
)
TRIPLE_GRAPH_RESIDUAL: Final = 25 * S**4 + 185 * S**3 - 413 * S**2 - 1497 * S + 228
TRIPLE_GRAPH_COLLISIONS: Final = expand(
    (5 * S + 1) * (5 * S + 2) * TRIPLE_GRAPH_RESIDUAL / 625
)

TRIPLE_NODE_X_CUBIC: Final = 625 * X**3 - 190861 * X**2 + 203652 * X - 20736
TRIPLE_NODE_X_QUARTIC: Final = (
    390625 * X**4 + 332908125 * X**3 + 22873309446 * X**2 + 7957013400 * X + 579156480
)

TRIPLE_WALL_IMPLICIT: Final = (
    -152587890625 * X**9
    + 50791992187500 * X**8
    - 8773818750000 * X**7
    + 99914550781250 * X**6 * Y
    + 505313640000 * X**6
    + 183958593750000 * X**5 * Y
    + 39825439453125 * X**4 * Y**2
    - 9701462016 * X**5
    - 22186507500000 * X**4 * Y
    - 284648681640625 * X**3 * Y**2
    + 4577636718750 * X**2 * Y**3
    + 648518400000 * X**3 * Y
    + 106483359375000 * X**2 * Y**2
    + 14560546875000 * X * Y**3
    + 152587890625 * Y**4
    - 10259156250000 * X * Y**2
    + 421875000000 * Y**3
    + 291600000000 * Y**2
)
TRIPLE_IMPLICIT_RESULTANT_MULTIPLIER: Final = 152587890625

TRIPLE_TANGENT_CONE: Final = expand(
    Rational(93312, 15625)
    * (18 * LOCAL_U + 21875 * LOCAL_V)
    * (168 * LOCAL_U - 3125 * LOCAL_V)
    * (1782 * LOCAL_U - 3125 * LOCAL_V)
)

# Sage 10.8 exact raw affine presentation.  Generator indices 1..4 are the
# geometric meridians in a generic vertical fiber.  The raw relators are kept
# so that the finite-image audit does not trust the cyclic Tietze reduction.
TRIPLE_WALL_RELATIONS: Final = (
    (4, 3, -4, -3),
    (2, 1, -2, -1),
    (
        -4,
        -3,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        3,
        4,
        3,
        -4,
        -3,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        3,
        4,
        -3,
    ),
    (2, 1, -2, -1),
    (
        -4,
        -3,
        -2,
        -1,
        3,
        -4,
        -3,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        3,
        4,
        -3,
        1,
        2,
        3,
        4,
        3,
    ),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-3, -2, -1, 3, -4, -3, 1, 2, 3, 4),
    (-2, -1, 2, 3, 4, -3, 1, 3, -4, -3),
    (-2, -1, -2, -1, -2, -1, 3, -4, -3, 1, 2, 1, 2, 1, 2, 1),
    (
        -2,
        -1,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        1,
        2,
        1,
        -2,
        -1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        1,
        2,
        -1,
    ),
    (-2, -1, 3, -4, -3, 1, 2, 1, -2, -1, 3, 4, -3, 1, 2, -1),
    (-4, -3, 1, 3, 4, -3, -1, 3),
)


def _homogeneous_part(polynomial: Expr, degree: int) -> Expr:
    """Return the total-degree ``degree`` homogeneous part in ``u,v``."""

    result: Expr = 0
    for monomial, coefficient in Poly(polynomial, LOCAL_U, LOCAL_V).terms():
        u_degree, v_degree = monomial
        if u_degree + v_degree == degree:
            result += coefficient * LOCAL_U**u_degree * LOCAL_V**v_degree
    return expand(result)


def _terms_below_degree(polynomial: Expr, degree: int) -> Expr:
    """Return all terms of total degree strictly below ``degree``."""

    result: Expr = 0
    for monomial, coefficient in Poly(polynomial, LOCAL_U, LOCAL_V).terms():
        u_degree, v_degree = monomial
        if u_degree + v_degree < degree:
            result += coefficient * LOCAL_U**u_degree * LOCAL_V**v_degree
    return expand(result)


@cache
def _presentation_pruning() -> tuple[tuple[int, ...], tuple[Assignment, ...]]:
    """Replay incremental relation pruning and retain final assignments."""

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
            for relation in TRIPLE_WALL_RELATIONS
            if relation and max(abs(letter) for letter in relation) == assigned_count
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


@dataclass(frozen=True, slots=True)
class DeltaTenTripleWallCertificate:
    """Exact geometry and finite-image audit for one triple-wall member."""

    split_vertical_identity: Expr
    split_graph_identity: Expr
    split_vertical_discriminant: Expr
    split_graph_discriminant: Expr
    residual_vertical_discriminant: Expr
    residual_graph_discriminant: Expr
    triple_p_values: tuple[Expr, ...]
    triple_q_values: tuple[Expr, ...]
    triple_p_derivatives: tuple[Expr, ...]
    triple_slopes: tuple[Expr, ...]
    fourth_fiber_p_value: Expr
    fourth_fiber_q_value: Expr
    cusp_leading_coefficients: tuple[Expr, Expr]
    other_cusp_fiber_q_value: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    triple_lower_order_terms: Expr
    triple_tangent_cone_identity: Expr
    vertical_node_resultant_identity: Expr
    graph_node_resultant_identity: Expr
    cubic_node_discriminant: Expr
    quartic_node_discriminant: Expr
    node_component_resultant: Expr
    triple_node_separations: tuple[Expr, Expr]
    cusp_node_separations: tuple[Expr, Expr]
    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    implicit_content: Expr
    implicit_total_degree: int
    sage_jacobian_components: tuple[tuple[int, int, bool], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    arithmetic_genus: int
    cusp_delta: int
    triple_delta: int
    node_count: int
    infinity_delta: int
    infinity_orders: tuple[int, int]
    relation_count: int
    pruning_stage_survivors: tuple[int, ...]
    diagonal_satisfying_assignments: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return (
            self.cusp_delta + self.triple_delta + self.node_count + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether every exact sample and presentation invariant agrees."""

        return bool(
            self.split_vertical_identity == 0
            and self.split_graph_identity == 0
            and self.split_vertical_discriminant
            == Rational(953474574156827904, 19073486328125)
            and self.split_graph_discriminant
            == Rational(
                73833539268797569216535158849536,
                7450580596923828125,
            )
            and self.residual_vertical_discriminant == 12708946980
            and self.residual_graph_discriminant == 162050781640608000
            and self.triple_p_values == (TRIPLE_IMAGE_X,) * 3
            and self.triple_q_values == (TRIPLE_IMAGE_Y,) * 3
            and self.triple_p_derivatives
            == (Rational(12, 125), Rational(-12, 125), Rational(84, 125))
            and self.triple_slopes
            == (
                Rational(1782, 3125),
                Rational(168, 3125),
                Rational(-18, 21875),
            )
            and len(set(self.triple_slopes)) == 3
            and self.fourth_fiber_p_value == TRIPLE_IMAGE_X
            and self.fourth_fiber_q_value == Rational(-653184, 78125)
            and self.cusp_leading_coefficients == (Rational(1), Rational(114, 625))
            and self.other_cusp_fiber_q_value == Rational(-864, 625)
            and self.cusp_image_factor == Rational(746496, 390625)
            and self.extra_critical_factor == Rational(-9906624, 15625)
            and self.triple_lower_order_terms == 0
            and self.triple_tangent_cone_identity == 0
            and self.vertical_node_resultant_identity == 0
            and self.graph_node_resultant_identity == 0
            and self.cubic_node_discriminant == 922081535273396664720
            and self.quartic_node_discriminant
            == 403096338022260849807369901193666641448747609062500000000
            and self.node_component_resultant
            == 118871648755782266797686957816152064000000000000
            and self.triple_node_separations
            == (Rational(-150605568, 15625), Rational(695895141888, 625))
            and self.cusp_node_separations == (-20736, 579156480)
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.implicit_total_degree == 9
            and self.sage_jacobian_components
            == (
                (4, 1, False),
                (4, 1, False),
                (4, 4, True),
                (3, 3, True),
            )
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.triple_delta == 3
            and self.node_count == 7
            and self.infinity_delta == 16
            and self.infinity_orders == (5, 9)
            and self.total_delta == self.arithmetic_genus
            and self.relation_count == 12
            and self.pruning_stage_survivors == (40, 40, 1600, 40)
            and self.diagonal_satisfying_assignments == 40
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_triple_wall_certificate() -> DeltaTenTripleWallCertificate:
    """Build the exact ordinary-triple member and presentation certificate."""

    sample_vertical = expand(exceptional_vertical(1).subs(TRIPLE_WALL_PARAMETERS))
    sample_graph = expand(exceptional_graph(1).subs(TRIPLE_WALL_PARAMETERS))
    p_derivative = diff(TRIPLE_WALL_P, T)
    q_derivative = diff(TRIPLE_WALL_Q, T)
    p_derivative_values = tuple(
        p_derivative.subs(T, parameter) for parameter in TRIPLE_PARAMETERS
    )
    slopes = tuple(
        q_derivative.subs(T, parameter) / p_derivative.subs(T, parameter)
        for parameter in TRIPLE_PARAMETERS
    )
    translated_implicit = expand(
        TRIPLE_WALL_IMPLICIT.subs({X: LOCAL_U + TRIPLE_IMAGE_X, Y: LOCAL_V})
    )
    stage_survivors, assignments = _presentation_pruning()
    diagonal_assignments = sum(
        all(image == images[0] for image in images[1:]) for images in assignments
    )
    implicit_polynomial = Poly(TRIPLE_WALL_IMPLICIT, X, Y)

    return DeltaTenTripleWallCertificate(
        split_vertical_identity=expand(sample_vertical - TRIPLE_VERTICAL_COLLISIONS),
        split_graph_identity=expand(sample_graph - TRIPLE_GRAPH_COLLISIONS),
        split_vertical_discriminant=discriminant(sample_vertical, R),
        split_graph_discriminant=discriminant(sample_graph, S),
        residual_vertical_discriminant=discriminant(
            TRIPLE_VERTICAL_RESIDUAL,
            R,
        ),
        residual_graph_discriminant=discriminant(TRIPLE_GRAPH_RESIDUAL, S),
        triple_p_values=tuple(
            TRIPLE_WALL_P.subs(T, parameter) for parameter in TRIPLE_PARAMETERS
        ),
        triple_q_values=tuple(
            TRIPLE_WALL_Q.subs(T, parameter) for parameter in TRIPLE_PARAMETERS
        ),
        triple_p_derivatives=p_derivative_values,
        triple_slopes=slopes,
        fourth_fiber_p_value=TRIPLE_WALL_P.subs(
            T,
            FOURTH_FIBER_PARAMETER,
        ),
        fourth_fiber_q_value=TRIPLE_WALL_Q.subs(
            T,
            FOURTH_FIBER_PARAMETER,
        ),
        cusp_leading_coefficients=(
            Poly(TRIPLE_WALL_P, T).nth(2),
            Poly(TRIPLE_WALL_Q, T).nth(5),
        ),
        other_cusp_fiber_q_value=TRIPLE_WALL_Q.subs(T, -1),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(TRIPLE_WALL_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(TRIPLE_WALL_PARAMETERS),
        triple_lower_order_terms=_terms_below_degree(translated_implicit, 3),
        triple_tangent_cone_identity=expand(
            _homogeneous_part(translated_implicit, 3) - TRIPLE_TANGENT_CONE
        ),
        vertical_node_resultant_identity=expand(
            resultant(TRIPLE_VERTICAL_RESIDUAL, R**2 - X, R) + TRIPLE_NODE_X_CUBIC
        ),
        graph_node_resultant_identity=expand(
            resultant(
                TRIPLE_GRAPH_RESIDUAL,
                -S * (S + 2) * (S + 1) ** 2 - 4 * X,
                S,
            )
            - 256 * TRIPLE_NODE_X_QUARTIC
        ),
        cubic_node_discriminant=discriminant(TRIPLE_NODE_X_CUBIC, X),
        quartic_node_discriminant=discriminant(TRIPLE_NODE_X_QUARTIC, X),
        node_component_resultant=resultant(
            TRIPLE_NODE_X_CUBIC,
            TRIPLE_NODE_X_QUARTIC,
            X,
        ),
        triple_node_separations=(
            TRIPLE_NODE_X_CUBIC.subs(X, TRIPLE_IMAGE_X),
            TRIPLE_NODE_X_QUARTIC.subs(X, TRIPLE_IMAGE_X),
        ),
        cusp_node_separations=(
            TRIPLE_NODE_X_CUBIC.subs(X, 0),
            TRIPLE_NODE_X_QUARTIC.subs(X, 0),
        ),
        implicit_resultant_identity=expand(
            TRIPLE_IMPLICIT_RESULTANT_MULTIPLIER
            * resultant(TRIPLE_WALL_P - X, TRIPLE_WALL_Q - Y, T)
            - TRIPLE_WALL_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            TRIPLE_WALL_IMPLICIT.subs({X: TRIPLE_WALL_P, Y: TRIPLE_WALL_Q})
        ),
        implicit_content=implicit_polynomial.content(),
        implicit_total_degree=implicit_polynomial.total_degree(),
        # Exact Sage primary-component (length, radical-degree, reduced)
        # output: cusp, triple point, four nodes, and three nodes.
        sage_jacobian_components=(
            (4, 1, False),
            (4, 1, False),
            (4, 4, True),
            (3, 3, True),
        ),
        # Simplified presentation: one generator, no relators; all four raw
        # geometric meridians map to that generator.
        sage_cyclic_simplification=(1, 0, True),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=(2 - 1) * (5 - 1) // 2,
        triple_delta=3,
        node_count=7,
        infinity_delta=(5 - 1) * (9 - 1) // 2,
        infinity_orders=(5, 9),
        relation_count=len(TRIPLE_WALL_RELATIONS),
        pruning_stage_survivors=stage_survivors,
        diagonal_satisfying_assignments=diagonal_assignments,
        complement_census=n_generator_three_cycle_presentation_census(
            TRIPLE_WALL_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact triple-wall certificate and fail on regression."""

    certificate = exact_delta_ten_triple_wall_certificate()
    print(
        "delta-ten ordinary-triple member:",
        {
            "cusp": (2, 5),
            "triple branches": 3,
            "nodes": certificate.node_count,
            "infinity": certificate.infinity_orders,
            "delta": certificate.total_delta,
        },
    )
    print(
        "k=2 split collision degrees:",
        {"vertical": (1, 3), "graph": (1, 1, 4)},
    )
    print("presentation pruning:", certificate.pruning_stage_survivors)
    print(
        "single-three-cycle images:",
        dict(certificate.complement_census.generated_order_histogram),
    )
    print("exact ordinary-triple member excluded from A6:", certificate.verified)
    print(
        "claim boundary: exact conditional/computer-assisted member only; "
        "wall propagation, deeper intersections, and JC(2) remain open"
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
