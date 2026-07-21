"""Tests for the conditional generic ``A6`` delta-seven exclusion."""

from sympy import Rational, expand

from scripts.a6_delta_seven_generic import (
    ALPHA,
    BETA,
    CLEAN_H,
    COLLISION_POLYNOMIAL,
    DELTA,
    GENERIC_RELATIONS,
    GAMMA,
    S,
    exact_a6_delta_seven_clean_certificate,
    exact_a6_delta_seven_family_certificate,
)
from scripts.a6_one_pair_infinity import three_cycle_presentation_census


def test_four_parameter_normal_form_and_collision_algebra_are_exact() -> None:
    certificate = exact_a6_delta_seven_family_certificate()

    assert certificate.normal_form_identity == 0
    assert certificate.collision_reduction_remainder == 0
    assert certificate.tangency_reduction_remainder == 0
    assert certificate.pair_resultant_identity == 0
    assert certificate.derivative_factor_identity == 0
    assert certificate.cusp_collision_identity == 0
    assert certificate.triple_resultant_identity == 0
    assert certificate.discriminant_factor_degrees == ((1, 1), (10, 1))
    assert certificate.tangency_resultant_identity == 0
    assert certificate.clean_open_values == (2, 2110, 1, -25, 464_000_000)
    assert certificate.verified


def test_clean_member_has_exact_cusp_node_and_genus_census() -> None:
    certificate = exact_a6_delta_seven_clean_certificate()

    assert certificate.finite_cusp_pair == (2, 5)
    assert certificate.node_count == 7
    assert certificate.infinity_pair == (7, 10)
    assert certificate.arithmetic_genus == 36
    assert certificate.total_delta == 36
    assert certificate.collision_discriminant == 464_000_000
    assert certificate.pair_diagonal_resultant == -4220
    assert certificate.pair_product_resultant == 2
    assert certificate.tangency_resultant == -1_958_080_000_000
    assert certificate.node_coordinate_discriminant == 113_281_250_000_000_000
    assert certificate.cusp_coordinate_separation == 2
    assert certificate.cusp_image_separation == -1
    assert certificate.residual_derivative_value == Rational(33760, 19683)
    assert certificate.verified


def test_clean_complement_has_only_cyclic_three_cycle_images() -> None:
    census = three_cycle_presentation_census(GENERIC_RELATIONS)

    assert census.assignments == 40**3
    assert census.satisfying_assignments == 40
    assert census.generated_order_histogram == ((3, 40),)
    assert census.a6_assignments == 0


def test_collision_formula_is_coefficient_sensitive() -> None:
    corrupted = COLLISION_POLYNOMIAL + ALPHA * S**4
    clean_substitution = {ALPHA: 2, BETA: 0, GAMMA: 0, DELTA: 0}

    assert expand(corrupted.subs(clean_substitution) - CLEAN_H) == 2 * S**4
