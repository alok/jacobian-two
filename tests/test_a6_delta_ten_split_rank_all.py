"""Exact regression for complete true-split rank coverage."""

from scripts.a6_delta_ten_split_rank_all import (
    EXPECTED_PROFILE_COUNTS,
    exact_all_split_rank_certificate,
)


def test_all_twenty_two_generated_allocation_keys_are_audited_once() -> None:
    certificate = exact_all_split_rank_certificate()

    assert len(certificate.expected_keys) == 22
    assert len(certificate.audited_keys) == 22
    assert len(set(certificate.audited_keys)) == 22
    assert set(certificate.audited_keys) == set(certificate.expected_keys)
    assert certificate.profile_counts == EXPECTED_PROFILE_COUNTS
    assert certificate.allocation_coverage_exact


def test_no_true_split_rank_stratum_can_support_a_threefold() -> None:
    certificate = exact_all_split_rank_certificate()

    assert certificate.contact_rank_verified
    assert certificate.triple_mixed_rank_verified
    assert certificate.global_fiber_rank_verified
    assert certificate.maximum_split_incidence_dimension == 2
    assert certificate.maximum_residual_rankdrop_incidence_dimension == 1
    assert certificate.all_rank_strata_classified


def test_aggregate_preserves_the_topology_and_jc_claim_boundary() -> None:
    certificate = exact_all_split_rank_certificate()

    assert certificate.triple_mixed_maximal_rank_closure_verified
    assert certificate.contact_rank_open_closure_verified
    assert certificate.contact_exceptional_topology_open
    assert certificate.deeper_boundaries_open
    assert not certificate.proves_plane_jacobian_conjecture
    assert certificate.verified
