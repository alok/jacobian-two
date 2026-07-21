"""Tests for the exact conditional delta-ten contact-two member."""

from sympy import Rational

from scripts.a6_delta_ten_contact_wall import (
    CONTACT_RELATIONS,
    exact_delta_ten_contact_wall_certificate,
)


def test_contact_wall_collision_and_validity_factors_are_exact() -> None:
    certificate = exact_delta_ten_contact_wall_certificate()

    assert certificate.collision_factor_identity == 0
    assert certificate.residual_collision_discriminant == -52320062156396485500
    assert certificate.contact_residual_separation == 1269
    assert certificate.denominator_value == Rational(81, 1024)
    assert certificate.cusp_image_factor == Rational(133, 25)
    assert certificate.extra_critical_factor == Rational(852228, 25)
    assert certificate.pair_diagonal_resultant == Rational(-828365616, 125)


def test_double_collision_is_exactly_a_contact_two_point() -> None:
    certificate = exact_delta_ten_contact_wall_certificate()

    assert certificate.tangency_factor_identity == 0
    assert certificate.contact_tangency_cofactor_value == 10152
    assert certificate.residual_tangency_resultant == Rational(
        29153413543712129429066042292,
        6103515625,
    )
    assert certificate.contact_incidence_identities == (0, 0)
    assert certificate.contact_pair_discriminant == -4
    assert certificate.contact_image_remainders == (0, 0)
    assert certificate.contact_slope_remainder == 0
    assert certificate.curvature_denominator_resultants == (140608, 140608)
    assert certificate.curvature_difference_identity == 0
    assert certificate.weighted_initial_identity == 0
    assert certificate.weighted_initial_discriminant == -4734540864


def test_residual_collisions_are_eight_distinct_nodes() -> None:
    certificate = exact_delta_ten_contact_wall_certificate()

    expected_discriminant = -(
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
    assert certificate.node_x_resultant_identity == 0
    assert certificate.node_x_discriminant == expected_discriminant
    assert certificate.contact_node_separation == 81847116
    assert certificate.cusp_node_separation == -636804
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_component_lengths == (
        (4, 1),
        (3, 1),
        (8, 8),
    )
    assert certificate.total_delta == 28
    assert certificate.total_delta == certificate.arithmetic_genus


def test_raw_contact_wall_presentation_has_only_diagonal_c3_images() -> None:
    certificate = exact_delta_ten_contact_wall_certificate()

    assert len(CONTACT_RELATIONS) == 12
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.pruning_stage_survivors == (40, 160, 640, 40)
    assert certificate.diagonal_satisfying_assignments == 40
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0


def test_full_delta_ten_contact_wall_certificate_is_verified() -> None:
    assert exact_delta_ten_contact_wall_certificate().verified
