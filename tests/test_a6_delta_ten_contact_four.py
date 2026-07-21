"""Tests for the conditional delta-ten ``C4 + 6N`` profile."""

from sympy import Rational, expand

from scripts.a6_delta_ten_contact_four import (
    CONTACT_FOUR_RELATIONS,
    exact_delta_ten_contact_four_certificate,
)
from scripts.a6_delta_ten_generic import S, T


def test_contact_four_incidence_has_one_dominant_surface() -> None:
    certificate = exact_delta_ten_contact_four_certificate()

    assert certificate.incidence_determinant_identity == 0
    assert certificate.residual_factor_irreducible
    assert certificate.coefficient_rank_two_minor_gcd_identity == 0
    assert certificate.specialized_augmented_gcd == 1
    assert certificate.sage_coefficient_rank_two_saturation == (True, 2)
    assert certificate.sage_residual_compatibility_saturation == (0, 10, 4)
    assert certificate.sage_augmented_rank_two_saturation == (True, 2)
    assert certificate.sage_augmented_rank_one_saturation == (True, 1)
    assert certificate.cramer_open_component_count == 1
    assert certificate.cramer_open_component_dimension == 2
    assert certificate.coefficient_image_codimension == 3
    assert certificate.projection_quasi_finite_on_cramer_open
    assert certificate.projection_generically_finite_onto_image
    assert certificate.residual_coefficient_rank == 3
    assert certificate.residual_incidence_dimension == 1
    assert certificate.residual_surface_component_excluded
    assert not certificate.split_boundaries_classified_here
    assert certificate.topology_propagation_dependencies == (
        "connected clean open",
        "proper projective Whitney-Thom triviality",
    )
    assert certificate.sample_incidence_residuals == (0, 0, 0, 0)
    assert certificate.sample_incidence_determinant == -168
    assert certificate.sample_valid_localizer == -9
    assert certificate.sample_cusp_coefficient == Rational(9, 7)


def test_rational_member_has_one_exact_contact_four_root() -> None:
    certificate = exact_delta_ten_contact_four_certificate()

    assert certificate.collision_identity == 0
    assert certificate.tangency_identity == 0
    assert certificate.collision_tangency_gcd == expand((S + 1) ** 3)
    assert certificate.quadruple_root_gcd == S + 1
    assert certificate.fourth_derivative == Rational(492, 7)
    assert certificate.tangency_third_derivative == Rational(1476, 7)
    assert certificate.residual_discriminant == Rational(
        -2627098335909,
        2582630848,
    )
    assert certificate.residual_contact_separation == Rational(41, 14)
    assert certificate.cusp_image_factor == Rational(1, 196)
    assert certificate.extra_critical_factor == Rational(11421, 49)
    assert certificate.denominator_resultant == 81
    assert certificate.diagonal_resultant == Rational(8325909, 343)
    assert certificate.residual_tangency_resultant == Rational(
        763059488663600977727961,
        173625106649344,
    )
    assert certificate.contact_chart_values == (-1, 1, 3)
    assert certificate.contact_pair_discriminant == -3


def test_local_jets_prove_intersection_multiplicity_four() -> None:
    certificate = exact_delta_ten_contact_four_certificate()

    assert certificate.contact_image_remainders == (0, -1)
    assert certificate.contact_p_derivative_resultant == 3
    assert certificate.contact_jet_values == (6, -30, 108)
    assert certificate.contact_jet_differences == (
        0,
        0,
        0,
        984 * (2 * T + 1),
    )
    assert certificate.contact_fourth_jet_separation == 3


def test_six_nodes_and_projective_genus_budget_are_exact() -> None:
    certificate = exact_delta_ten_contact_four_certificate()

    assert certificate.node_x_resultant_identity == 0
    assert certificate.node_x_discriminant != 0
    assert certificate.contact_node_separation == 738
    assert certificate.cusp_node_separation == 738
    assert certificate.contact_cusp_target_separation == -1
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_components == ((4, 1), (7, 1), (6, 6))
    assert certificate.total_delta == 28
    assert certificate.total_delta == certificate.arithmetic_genus


def test_cyclic_complement_has_no_required_a6_image() -> None:
    certificate = exact_delta_ten_contact_four_certificate()

    assert len(CONTACT_FOUR_RELATIONS) == 10
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0
    assert certificate.verified
