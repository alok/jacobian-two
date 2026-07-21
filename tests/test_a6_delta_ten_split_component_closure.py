"""Exact regressions for true-split component containment."""

from scripts.a6_delta_ten_split_component_closure import (
    exact_split_component_closure_certificate,
    exact_split_component_embedding_certificates,
    exact_total_mixed_base_certificate,
)


def test_total_mixed_source_base_is_primitive_linear_and_irreducible() -> None:
    certificate = exact_total_mixed_base_certificate()

    assert certificate.degree_in_contact_product == 1
    assert certificate.coefficient_gcd == 1
    assert certificate.nonsplit_sample_identity == 0
    assert certificate.irreducible
    assert certificate.absolutely_irreducible
    assert certificate.verified


def test_all_clean_split_row_systems_restrict_from_global_incidences() -> None:
    certificates = exact_split_component_embedding_certificates()

    assert len(certificates) == 7
    assert sum(item.profile == "T112+6N" for item in certificates) == 4
    assert sum(item.profile == "C2+T111+5N" for item in certificates) == 3
    for certificate in certificates:
        assert all(value == 0 for value in certificate.global_row_residuals)
        assert certificate.total_base_identity == 0
        assert certificate.kappa_identity == 0
        assert certificate.transformation_determinant != 0
        assert certificate.witness_transformation_determinant != 0
        assert any(value != 0 for value in certificate.witness_base_gradient)
        assert certificate.clean_witness_verified
        assert certificate.verified


def test_only_the_rank_three_affine_line_fibers_remain_topologically_open() -> None:
    certificate = exact_split_component_closure_certificate()

    assert certificate.total_mixed_base.verified
    assert certificate.t112_global_incidence_excluded
    assert certificate.maximal_rank_split_topology_closed
    assert certificate.representative_exceptional_rank_three_base_length == 8
    assert certificate.exceptional_rank_three_fiber_topology_open
    assert certificate.verified
