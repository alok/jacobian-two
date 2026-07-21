"""Tests for the smallest width-five S6 branch-curve near-miss."""

from scripts.s6_width_five_curve import (
    VAN_KAMPEN_RELATIONS,
    exact_s6_width_five_curve_certificate,
    exact_width_five_curve_arithmetic_certificate,
    width_five_transposition_census,
)


def test_width_five_curve_has_the_exact_allowed_geometric_budget() -> None:
    certificate = exact_width_five_curve_arithmetic_certificate()

    assert certificate.finite_cusp_pair == (4, 5)
    assert certificate.infinity_pair == (2, 7)
    assert certificate.finite_node_count == 6
    assert certificate.arithmetic_genus == 15
    assert certificate.total_delta == 15
    assert certificate.width == 5
    assert certificate.verified


def test_five_meridian_presentation_has_only_cyclic_transposition_images() -> None:
    census = width_five_transposition_census()

    assert len(VAN_KAMPEN_RELATIONS) == 10
    assert census.total_assignments == 15**5
    assert census.satisfying_assignments == 15
    assert census.generated_order_histogram == ((2, 15),)
    assert census.transitive_assignments == 0
    assert census.verified


def test_complete_width_five_certificate() -> None:
    certificate = exact_s6_width_five_curve_certificate()

    assert certificate.minimum_s6_width == 5
    assert certificate.width_is_minimal
    assert certificate.verified
