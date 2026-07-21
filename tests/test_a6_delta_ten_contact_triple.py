"""Tests for the conditional delta-ten ``C2 + T111 + 5N`` stratum."""

from sympy import Rational

from scripts.a6_delta_ten_contact_triple import (
    CONTACT_TRIPLE_RELATIONS,
    exact_delta_ten_contact_triple_certificate,
)
from scripts.a6_delta_ten_generic import S, T


def test_combined_incidence_has_one_dominant_cramer_component() -> None:
    certificate = exact_delta_ten_contact_triple_certificate()

    assert certificate.incidence_determinant_identity == 0
    assert certificate.rank_factor_irreducible
    assert certificate.specialized_rank_identity == 0
    assert certificate.specialized_augmented_identity == 0
    assert certificate.specialized_augmented_gcd == 1
    assert certificate.cramer_open_component_count == 1
    assert certificate.cramer_open_component_dimension == 3
    assert certificate.coefficient_image_codimension == 2
    assert certificate.projection_quasi_finite_on_cramer_open
    assert certificate.projection_generically_finite_onto_image
    assert certificate.sage_residual_rank_two_saturation == (1, 9, True)
    assert certificate.residual_rank_threefold_excluded
    assert not certificate.residual_rank_intersections_classified
    assert not certificate.split_boundary_intersections_classified
    assert certificate.topology_propagation_dependencies == (
        "connected clean open",
        "proper projective Whitney-Thom triviality",
    )
    assert certificate.sample_incidence_residuals == (0, 0, 0, 0)
    assert certificate.sample_incidence_determinant == 276710448
    assert certificate.sample_linear_solution == (-2, 1, -2, 0)
    assert all(value != 0 for value in certificate.sample_open_factor_values)
    assert all(identity == 0 for identity in certificate.boundary_identities)


def test_collision_decic_splits_as_contact_triple_and_five_nodes() -> None:
    certificate = exact_delta_ten_contact_triple_certificate()

    assert certificate.collision_identity == 0
    assert certificate.tangency_identity == 0
    assert certificate.collision_tangency_gcd == S - 1
    assert certificate.contact_collision_second_derivative == Rational(2401, 2)
    assert certificate.contact_tangency_derivative == Rational(50421, 4)
    assert certificate.residual_discriminant == -2138993541120
    assert certificate.residual_factor_resultants == (-343, 63315, -14)
    assert certificate.tangency_factor_resultants == (
        -39058273766241112220775,
        Rational(1759080645, 1024),
    )
    assert certificate.pair_boundary_resultants == (
        Rational(-9, 2),
        Rational(21609, 1024),
        Rational(103243, 512),
        Rational(-15993922959, 4096),
    )
    assert certificate.cusp_image_factor == Rational(43, 2)
    assert certificate.extra_critical_factor == Rational(740151, 8)


def test_local_jets_give_exact_contact_two_and_ordinary_triple() -> None:
    certificate = exact_delta_ten_contact_triple_certificate()

    assert certificate.contact_image_remainders == (-5, -2)
    assert certificate.triple_image_remainders == (1, 4)
    assert certificate.fourth_source_image == (1, Rational(-31, 512))
    assert certificate.source_derivative_resultants == (237, -63)
    assert certificate.contact_jet_differences == (
        0,
        Rational(392, 6241) * (2 * T - 1),
    )
    assert certificate.contact_second_jet_separation == 3
    assert certificate.triple_slope_resultant_identity == 0
    assert certificate.triple_slope_discriminant == -344755200


def test_five_nodes_and_projective_genus_budget_are_exact() -> None:
    certificate = exact_delta_ten_contact_triple_certificate()

    assert certificate.node_x_resultant_identity == 0
    assert certificate.node_x_discriminant == (
        -17004216373814107665135371132504044296002359787520
    )
    assert certificate.node_target_separations == (
        -2746044819,
        325050165,
        115495936,
    )
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_components == (
        (4, 1),
        (3, 1),
        (4, 1),
        (5, 5),
    )
    assert certificate.total_delta == 28
    assert certificate.total_delta == certificate.arithmetic_genus


def test_cyclic_complement_has_no_required_a6_image() -> None:
    certificate = exact_delta_ten_contact_triple_certificate()

    assert len(CONTACT_TRIPLE_RELATIONS) == 11
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0
    assert certificate.verified
