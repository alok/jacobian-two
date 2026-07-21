"""Tests for the finite delta-seven collision endpoints."""

from scripts.a6_delta_seven_finite_wall import (
    P33_COMPLEX_SIMPLIFIED_RELATIONS,
    P42_COMPLEX_SIMPLIFIED_RELATIONS,
    P42_REAL_ONE_RELATIONS,
    P42_REAL_ZERO_RELATIONS,
    P5_NEGATIVE_RELATIONS,
    P5_POSITIVE_RELATIONS,
    checked_endpoint_presentations,
    exact_finite_wall_algebra_certificate,
)
from scripts.a6_one_pair_infinity import three_cycle_presentation_census


def test_finite_partition_classification_is_exact() -> None:
    certificate = exact_finite_wall_algebra_certificate()

    assert certificate.contact_five_factor_remainder == 0
    assert certificate.contact_five_elimination_identity == 0
    assert len(certificate.contact_five_saturated_basis) == 3
    assert certificate.contact_five_discarded_cusp_identity == 0
    assert certificate.contact_five_discriminant == 1377
    assert certificate.contact_five_residual_discriminant_resultant == -50
    assert certificate.contact_five_residual_value_resultant == 4896
    assert all(value != 0 for value in certificate.contact_five_invalid_resultants)
    assert certificate.contact_four_chart_denominator_resultants == (-162, 324)
    assert certificate.contact_four_two_discriminant == -108800
    assert certificate.contact_four_two_elimination_identity == 0
    assert len(certificate.contact_four_two_saturated_basis) == 3
    assert certificate.contact_four_two_discarded_branch_gcds
    assert certificate.contact_four_two_discarded_wall_identities == (0, 0, 0)
    assert certificate.contact_four_two_factor_remainder == 0
    assert certificate.contact_four_two_distinct_root_resultant == 720
    assert certificate.contact_four_two_contact_separation_resultant == 936360000
    assert all(value != 0 for value in certificate.contact_four_two_invalid_resultants)
    assert certificate.contact_three_three_first_slice_remainder == 0
    assert certificate.contact_three_three_second_slice_remainder == 0
    assert certificate.contact_three_three_elimination_identity == 0
    assert len(certificate.contact_three_three_saturated_basis) == 3
    assert certificate.contact_three_three_discarded_branch_gcds
    assert certificate.contact_three_three_degenerate_cusp_values == (0, 0)
    assert certificate.contact_three_three_discriminant == -13996800
    assert all(value != 0 for value in certificate.contact_three_three_invalid_resultants)
    assert certificate.contact_three_two_two_elimination_identity == 0
    assert certificate.contact_three_two_two_saturated_basis == (1,)
    assert len(certificate.contact_three_two_two_branch_gcds) == 4
    assert certificate.contact_three_two_two_branch_wall_values == (0, 0, 0, 0)
    assert certificate.contact_three_two_two_boundary_identity == 0
    assert certificate.contact_three_two_two_t_value == 0
    assert certificate.contact_three_two_two_c_value == 0
    assert certificate.two_root_cusp_values == (0, 0, 0)
    assert certificate.verified


def test_every_stored_number_field_embedding_has_no_a6_image() -> None:
    endpoints = checked_endpoint_presentations()

    assert len(endpoints) == 6
    assert tuple(endpoint.partition for endpoint in endpoints) == (
        (5, 1, 1),
        (5, 1, 1),
        (4, 2, 1),
        (4, 2, 1),
        (4, 2, 1),
        (3, 3, 1),
    )
    assert all(endpoint.verified for endpoint in endpoints)


def test_endpoint_relation_words_have_stable_lengths() -> None:
    assert tuple(map(len, P5_NEGATIVE_RELATIONS)) == (4, 8, 40, 18, 20)
    assert tuple(map(len, P5_POSITIVE_RELATIONS)) == (16, 30, 160, 20, 6)
    assert tuple(map(len, P42_REAL_ZERO_RELATIONS)) == (30, 32, 40, 32, 6)
    assert tuple(map(len, P42_REAL_ONE_RELATIONS)) == (10, 8, 16, 6, 8)
    assert tuple(map(len, P42_COMPLEX_SIMPLIFIED_RELATIONS)) == (8, 8, 18, 30, 88)
    assert tuple(map(len, P33_COMPLEX_SIMPLIFIED_RELATIONS)) == (10, 10, 10, 16, 20)


def test_complex_quartic_simplification_keeps_only_cyclic_images() -> None:
    census = three_cycle_presentation_census(P42_COMPLEX_SIMPLIFIED_RELATIONS)

    assert census.assignments == 40**3
    assert census.satisfying_assignments == 40
    assert census.generated_order_histogram == ((3, 40),)
    assert census.a6_assignments == 0
