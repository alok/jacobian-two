"""Exact regressions for maximal-rank split contact component closure."""

from sympy import cancel

from scripts.a6_delta_ten_split_contact_closure import (
    ORDERED_DOUBLE_CONTACT_BRIDGE_LOCALIZER,
    OVERLAP_LAMBDA,
    exact_contact_three_jet_bridge_certificate,
    exact_contact_three_overlap_degeneration_certificates,
    exact_maximal_rank_contact_embedding_certificates,
    exact_ordered_double_contact_jet_bridge_certificate,
    exact_split_contact_closure_certificate,
    exact_total_contact_base_certificate,
)
from scripts.a6_delta_ten_generic import KAPPA, PAIR_DENOMINATOR, S


def test_total_one_and_two_pair_bases_are_geometrically_integral() -> None:
    certificate = exact_total_contact_base_certificate()

    assert certificate.pair_surface_geometrically_integral
    assert certificate.pair_surface_flat_over_kappa
    assert certificate.generic_pair_product_degree == 1
    assert certificate.generic_pair_linear_coefficient_degree == 1
    assert certificate.generic_pair_constant_coefficient_degree == 3
    assert certificate.generic_pair_coefficient_gcd == 1
    assert certificate.generic_pair_coefficient_resultant == KAPPA**3 - 4 * KAPPA
    assert certificate.generic_pair_curve_primitive_linear
    assert certificate.generic_pair_curve_geometrically_integral
    assert certificate.ordered_generic_product_domain_consequence
    assert certificate.ordered_two_pair_base_flat_over_kappa
    assert certificate.ordered_two_pair_generic_fiber_domain
    assert certificate.ordered_two_pair_torsion_descent
    assert certificate.ordered_two_pair_base_domain
    assert certificate.verified


def test_total_contact_three_rows_are_exactly_the_old_h_jet_system() -> None:
    certificate = exact_contact_three_jet_bridge_certificate()
    unit = certificate.determinant_unit

    assert certificate.transformation.shape == (3, 3)
    assert certificate.target_scaling_residual == 0
    assert certificate.fiber_chain_rule_residuals == (0, 0)
    assert all(value == 0 for value in certificate.row_residuals)
    assert cancel(unit.element - S**6 / PAIR_DENOMINATOR**9) == 0
    assert cancel(unit.inverse - PAIR_DENOMINATOR**9 / S**6) == 0
    assert unit.ambient_factorization_residual == 0
    assert unit.element_denominator_membership_residual == 0
    assert unit.inverse_denominator_membership_residual == 0
    assert unit.verified
    assert certificate.verified


def test_two_total_contact_blocks_are_the_ordered_h_hprime_system() -> None:
    certificate = exact_ordered_double_contact_jet_bridge_certificate()
    unit = certificate.determinant_unit

    assert certificate.transformation.shape == (4, 4)
    assert all(value == 0 for value in certificate.row_residuals)
    assert unit.localizer == ORDERED_DOUBLE_CONTACT_BRIDGE_LOCALIZER
    assert unit.product_residual == 0
    assert unit.element_denominator_membership_residual == 0
    assert unit.inverse_denominator_membership_residual == 0
    assert unit.verified
    assert certificate.verified


def test_principal_c3_arcs_dominate_both_split_overlap_planes() -> None:
    certificates = exact_contact_three_overlap_degeneration_certificates()

    assert tuple(item.allocation for item in certificates) == (
        "overlap-W",
        "overlap-V",
    )
    assert certificates[0].cramer_determinant_leading_coefficient == (
        -3 * (OVERLAP_LAMBDA - 3) / 4
    )
    assert certificates[1].cramer_determinant_leading_coefficient == (
        -OVERLAP_LAMBDA / 2
    )
    for certificate in certificates:
        assert certificate.pair_incidence_residual == 0
        assert certificate.cramer_determinant != 0
        assert all(value == 0 for value in certificate.cramer_residuals)
        assert all(value == 0 for value in certificate.coefficient_limit_residuals)
        assert all(value == 0 for value in certificate.overlap_equation_residuals)
        assert all(value == 0 for value in certificate.dominance_residuals)
        assert certificate.verified


def test_eight_maximal_rank_component_jet_systems_are_global_restrictions() -> None:
    certificates = exact_maximal_rank_contact_embedding_certificates()

    assert len(certificates) == 8
    assert sum(item.profile == "C3+7N" for item in certificates) == 3
    assert sum(item.profile == "C2^2+6N" for item in certificates) == 5
    for certificate in certificates:
        assert certificate.transformation_determinant != 0
        assert all(value == 0 for value in certificate.row_residuals)
        assert certificate.determinant_unit.localized_ideal_is_unit
        assert certificate.determinant_unit.verified
        assert certificate.witness_transformation_determinant != 0
        assert certificate.clean_witness_verified
        assert certificate.verified


def test_only_overlap_and_residual_affine_line_topology_remains() -> None:
    certificate = exact_split_contact_closure_certificate()

    assert certificate.contact_three_jet_bridge_verified
    assert certificate.ordered_double_contact_jet_bridge_verified
    assert certificate.global_contact_three_sample_cyclic
    assert certificate.global_double_contact_sample_cyclic
    assert certificate.global_double_contact_cramer_verified
    assert certificate.maximal_rank_algebraically_contained
    assert certificate.proper_isotopy_extension_recorded
    assert certificate.maximal_rank_topology_closed
    assert certificate.overlap_algebraic_closure_count == 3
    assert certificate.exceptional_allocation_count == 3
    assert certificate.residual_affine_line_base_length == 16
    assert certificate.exceptional_topology_open
    assert certificate.verified
