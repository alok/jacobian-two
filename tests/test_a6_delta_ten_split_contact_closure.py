"""Exact regressions for maximal-rank split contact component closure."""

from scripts.a6_delta_ten_split_contact_closure import (
    exact_maximal_rank_contact_embedding_certificates,
    exact_split_contact_closure_certificate,
    exact_total_contact_base_certificate,
)


def test_total_one_and_two_pair_bases_are_geometrically_integral() -> None:
    certificate = exact_total_contact_base_certificate()

    assert certificate.pair_surface_geometrically_integral
    assert certificate.pair_surface_flat_over_kappa
    assert certificate.ordered_two_pair_base_flat_over_kappa
    assert certificate.ordered_two_pair_generic_fiber_domain
    assert certificate.ordered_two_pair_base_domain
    assert certificate.verified


def test_eight_maximal_rank_component_jet_systems_are_global_restrictions() -> None:
    certificates = exact_maximal_rank_contact_embedding_certificates()

    assert len(certificates) == 8
    assert sum(item.profile == "C3+7N" for item in certificates) == 3
    assert sum(item.profile == "C2^2+6N" for item in certificates) == 5
    for certificate in certificates:
        assert certificate.transformation_determinant != 0
        assert all(value == 0 for value in certificate.row_residuals)
        assert certificate.witness_transformation_determinant != 0
        assert certificate.clean_witness_verified
        assert certificate.verified


def test_only_overlap_and_residual_affine_line_topology_remains() -> None:
    certificate = exact_split_contact_closure_certificate()

    assert certificate.global_contact_three_sample_cyclic
    assert certificate.global_double_contact_sample_cyclic
    assert certificate.global_double_contact_cramer_verified
    assert certificate.maximal_rank_topology_closed
    assert certificate.exceptional_allocation_count == 3
    assert certificate.residual_affine_line_base_length == 16
    assert certificate.exceptional_topology_open
    assert certificate.verified
