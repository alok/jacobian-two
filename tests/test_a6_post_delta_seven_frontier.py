"""Tests for the conditional large-link frontier after delta seven."""

from scripts.a6_one_pair_infinity import one_pair_degree_candidates
from scripts.a6_post_delta_seven_frontier import (
    exact_post_delta_seven_frontier_certificate,
)


def test_degree_candidates_through_delta_ten_are_exact() -> None:
    certificate = exact_post_delta_seven_frontier_certificate()

    assert tuple(
        candidate.affine_degrees
        for candidate in certificate.delta_eight_candidates
    ) == ((2, 21), (3, 11))
    assert tuple(
        candidate.affine_degrees
        for candidate in certificate.delta_nine_candidates
    ) == ((2, 23),)
    assert tuple(
        candidate.affine_degrees
        for candidate in certificate.delta_ten_candidates
    ) == ((2, 25), (3, 13), (4, 9), (5, 7))

    unrestricted_delta_eight = one_pair_degree_candidates(
        8,
        require_singular_infinity=False,
    )
    assert tuple(c.affine_degrees for c in unrestricted_delta_eight) == (
        (2, 21),
        (3, 11),
        (5, 6),
    )


def test_exact_censuses_skip_eight_and_nine_then_select_t49() -> None:
    certificate = exact_post_delta_seven_frontier_certificate()

    assert tuple(
        census.single_three_cycle_a6_pairs for census in certificate.censuses
    ) == (0, 0, 0, 0, 0, 720, 0)
    assert certificate.censuses[3].a6_meridian_histogram == (
        ((4, 2), 1440),
        ((5, 1), 1440),
    )
    assert certificate.censuses[5].a6_meridian_histogram == (
        ((3, 1, 1, 1), 720),
        ((3, 3), 720),
        ((5, 1), 2880),
    )


def test_full_post_delta_seven_certificate_is_verified() -> None:
    assert exact_post_delta_seven_frontier_certificate().verified
