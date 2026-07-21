"""Tests for the exact conditional delta-ten ordinary-triple member."""

from sympy import Rational

from scripts.a6_delta_ten_triple_wall import (
    TRIPLE_GRAPH_COLLISIONS,
    TRIPLE_VERTICAL_COLLISIONS,
    TRIPLE_WALL_PARAMETERS,
    TRIPLE_WALL_RELATIONS,
    exact_delta_ten_triple_wall_certificate,
)
from scripts.a6_delta_ten_generic import ALPHA, KAPPA


def test_split_kappa_two_collision_fiber_is_exact_and_reduced() -> None:
    certificate = exact_delta_ten_triple_wall_certificate()

    assert TRIPLE_WALL_PARAMETERS[KAPPA] == 2
    assert TRIPLE_WALL_PARAMETERS[ALPHA] == Rational(114, 625)
    assert TRIPLE_VERTICAL_COLLISIONS.as_poly() is not None
    assert TRIPLE_GRAPH_COLLISIONS.as_poly() is not None
    assert certificate.split_vertical_identity == 0
    assert certificate.split_graph_identity == 0
    assert certificate.split_vertical_discriminant == Rational(
        953474574156827904,
        19073486328125,
    )
    assert certificate.split_graph_discriminant == Rational(
        73833539268797569216535158849536,
        7450580596923828125,
    )
    assert certificate.residual_vertical_discriminant == 12708946980
    assert certificate.residual_graph_discriminant == 162050781640608000


def test_three_parameters_give_one_ordinary_triple_and_not_a_quadruple() -> None:
    certificate = exact_delta_ten_triple_wall_certificate()

    assert certificate.triple_p_values == (Rational(36, 625),) * 3
    assert certificate.triple_q_values == (0, 0, 0)
    assert certificate.triple_p_derivatives == (
        Rational(12, 125),
        Rational(-12, 125),
        Rational(84, 125),
    )
    assert certificate.triple_slopes == (
        Rational(1782, 3125),
        Rational(168, 3125),
        Rational(-18, 21875),
    )
    assert len(set(certificate.triple_slopes)) == 3
    assert certificate.fourth_fiber_p_value == Rational(36, 625)
    assert certificate.fourth_fiber_q_value == Rational(-653184, 78125)
    assert certificate.triple_lower_order_terms == 0
    assert certificate.triple_tangent_cone_identity == 0


def test_cusp_nodes_infinity_and_genus_budget_are_exact() -> None:
    certificate = exact_delta_ten_triple_wall_certificate()

    assert certificate.cusp_leading_coefficients == (1, Rational(114, 625))
    assert certificate.other_cusp_fiber_q_value == Rational(-864, 625)
    assert certificate.cusp_image_factor == Rational(746496, 390625)
    assert certificate.extra_critical_factor == Rational(-9906624, 15625)
    assert certificate.vertical_node_resultant_identity == 0
    assert certificate.graph_node_resultant_identity == 0
    assert certificate.cubic_node_discriminant == 922081535273396664720
    assert certificate.quartic_node_discriminant == (
        403096338022260849807369901193666641448747609062500000000
    )
    assert certificate.node_component_resultant == (
        118871648755782266797686957816152064000000000000
    )
    assert certificate.triple_node_separations == (
        Rational(-150605568, 15625),
        Rational(695895141888, 625),
    )
    assert certificate.cusp_node_separations == (-20736, 579156480)
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.implicit_total_degree == 9
    assert certificate.sage_jacobian_components == (
        (4, 1, False),
        (4, 1, False),
        (4, 4, True),
        (3, 3, True),
    )
    assert certificate.total_delta == 28
    assert certificate.total_delta == certificate.arithmetic_genus
    assert certificate.infinity_orders == (5, 9)


def test_raw_triple_wall_presentation_has_only_diagonal_c3_images() -> None:
    certificate = exact_delta_ten_triple_wall_certificate()

    assert len(TRIPLE_WALL_RELATIONS) == 12
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.pruning_stage_survivors == (40, 40, 1600, 40)
    assert certificate.diagonal_satisfying_assignments == 40
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0


def test_full_delta_ten_triple_wall_certificate_is_verified() -> None:
    assert exact_delta_ten_triple_wall_certificate().verified
