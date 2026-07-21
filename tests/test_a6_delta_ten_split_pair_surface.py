"""Exact tests for the total unordered-pair surface and split fibers."""

from scripts.a6_delta_ten_split_pair_surface import (
    exact_split_pair_surface_certificate,
)


def test_total_pair_surface_is_smooth_irreducible_and_flat() -> None:
    certificate = exact_split_pair_surface_certificate()

    assert certificate.linear_content_gcd == 1
    assert certificate.singular_ideal_basis == (1,)
    assert certificate.total_surface_irreducible
    assert certificate.total_surface_smooth
    assert certificate.flat_over_kappa
    assert certificate.verified


def test_only_three_vertical_pair_fibers_and_their_intersections_are_exact() -> None:
    certificate = exact_split_pair_surface_certificate()

    assert certificate.exceptional_base_identity == 0
    assert certificate.split_fiber_identities == (0, 0, 0)
    assert certificate.component_intersection_residuals == (
        (0, 0),
        (0, 0),
        (0, 0),
    )
