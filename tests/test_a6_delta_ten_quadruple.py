"""Tests for the conditional delta-ten ordinary-quadruple incidence."""

from sympy import Rational

from scripts.a6_delta_ten_quadruple import (
    EXPECTED_COEFFICIENT_MINORS,
    FIBER_VALUE,
    exact_delta_ten_quadruple_certificate,
)


def test_quadruple_incidence_has_rank_three_on_every_nonzero_fiber() -> None:
    certificate = exact_delta_ten_quadruple_certificate()

    assert certificate.coefficient_minors == EXPECTED_COEFFICIENT_MINORS
    assert certificate.coefficient_minor_gcd == FIBER_VALUE**2
    assert certificate.localized_minor_basis == (1,)
    assert certificate.generic_solution_residuals == (0, 0, 0)
    assert certificate.sample_projection_derivatives == (
        Rational(-3, 2),
        Rational(-1, 2),
        Rational(-1, 4),
    )


def test_exact_member_has_one_reduced_four_source_fiber() -> None:
    certificate = exact_delta_ten_quadruple_certificate()

    assert certificate.sample_remainder == Rational(-1, 2)
    assert certificate.sample_fiber_discriminant == -279
    assert certificate.sample_quotient_identity == 0
    assert certificate.sample_slope_resultant_identity == 0
    assert certificate.sample_slope_discriminant == Rational(-9756351, 4096)
    assert certificate.pair_sum_resultant_identity == 0
    assert certificate.complementary_pair_identity == 0


def test_collision_decic_splits_as_six_quadruple_pairs_plus_four_nodes() -> None:
    certificate = exact_delta_ten_quadruple_certificate()

    assert certificate.collision_factor_identity == 0
    assert certificate.collision_factor_discriminants == (-31, -31, -458892)
    assert certificate.collision_factor_resultants == (-27, 153, -99)
    assert all(
        all(value != 0 for value in row) for row in certificate.validity_resultants
    )
    assert certificate.node_x_eliminant_identity == 0
    assert certificate.node_x_discriminant == -265262438203392
    assert certificate.node_x_at_quadruple == -187


def test_sample_stays_inside_the_standing_validity_open_and_closes_genus() -> None:
    certificate = exact_delta_ten_quadruple_certificate()

    assert certificate.cusp_image_factor == Rational(7, 4)
    assert certificate.extra_critical_factor == 13392
    assert certificate.affine_delta_total == 12
    assert certificate.projective_delta_total == 28


def test_raw_quadruple_presentation_has_only_diagonal_three_cycle_images() -> None:
    census = exact_delta_ten_quadruple_certificate().complement_census

    assert census.assignments == 40**4
    assert census.satisfying_assignments == 40
    assert census.generated_order_histogram == ((3, 40),)
    assert census.a6_assignments == 0


def test_full_delta_ten_ordinary_quadruple_certificate_is_verified() -> None:
    assert exact_delta_ten_quadruple_certificate().verified
