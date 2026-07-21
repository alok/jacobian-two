"""Tests for the conditional delta-ten two-triple profile."""

from sympy import Rational

from scripts.a6_delta_ten_double_triple import (
    exact_delta_ten_double_triple_certificate,
)


def test_double_triple_incidence_has_one_valid_dense_cramer_chart() -> None:
    certificate = exact_delta_ten_double_triple_certificate()

    assert certificate.incidence_determinant_identity == 0
    assert certificate.same_fiber_identity == 0
    assert certificate.rankdrop_coefficient_rank == 3
    assert certificate.rankdrop_augmented_rank == 4
    assert certificate.rankdrop_augmented_minor == -938448
    assert certificate.sample_equation_residuals == (0, 0, 0, 0)


def test_sample_has_two_separate_ordinary_triple_fibers() -> None:
    certificate = exact_delta_ten_double_triple_certificate()

    assert certificate.sample_cubic_discriminants == (-87, -255)
    assert certificate.sample_triple_values == (
        Rational(1, 3),
        Rational(-125, 3),
    )
    assert certificate.sample_omitted_value_differences == (-6, 60)
    assert certificate.sample_slope_discriminants == (
        Rational(-1560896, 3),
        Rational(-685443400000, 3),
    )
    assert certificate.pair_sum_identities == (0, 0)


def test_collision_decic_is_six_triple_pairs_plus_four_nodes() -> None:
    certificate = exact_delta_ten_double_triple_certificate()

    assert certificate.collision_factor_identity == 0
    assert certificate.collision_factor_discriminants == (-87, -255, 2108025)
    assert certificate.collision_factor_resultants == (360, 2175, 33675)
    assert all(
        all(value != 0 for value in row) for row in certificate.validity_resultants
    )
    assert certificate.cusp_image_factor == Rational(-5, 9)
    assert certificate.extra_critical_factor == 6409


def test_four_nodes_are_separate_and_genus_is_exact() -> None:
    certificate = exact_delta_ten_double_triple_certificate()

    assert certificate.node_x_eliminant_identity == 0
    assert certificate.node_x_discriminant == 402849850982400
    assert certificate.node_target_separations == (-232, 28736, 6)
    assert certificate.total_delta == certificate.arithmetic_genus == 28


def test_raw_presentation_is_cyclic_and_has_no_a6_image() -> None:
    certificate = exact_delta_ten_double_triple_certificate()

    assert certificate.sage_jacobian_components == (
        (4, 1),
        (4, 1),
        (4, 1),
        (4, 4),
    )
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.complement_census.assignments == 40**4
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0


def test_full_double_triple_certificate_is_verified() -> None:
    assert exact_delta_ten_double_triple_certificate().verified
