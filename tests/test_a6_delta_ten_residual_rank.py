"""Tests for the exact delta-ten residual-rank closure."""

from scripts.a6_delta_ten_generic import KAPPA
from scripts.a6_delta_ten_residual_rank import (
    ROOT_SUM,
    exact_delta_ten_residual_rank_certificate,
    exact_double_contact_residual_certificate,
    exact_double_triple_residual_certificate,
    normalized_double_triple_compatibility,
)


def test_two_triple_residual_is_empty_after_valid_chart_saturation() -> None:
    certificate = exact_double_triple_residual_certificate()
    normalized_minors, localizer = normalized_double_triple_compatibility()

    assert len(normalized_minors) == 4
    assert localizer != 0
    assert certificate.residual_substitution_identity == 0
    assert certificate.localized_groebner_basis == (1,)
    assert certificate.generic_residual_value == 0
    assert certificate.generic_localizer_value != 0
    assert certificate.generic_coefficient_rank == 3
    assert certificate.generic_augmented_rank == 4
    assert certificate.verified


def test_same_fiber_hostile_family_explains_two_triple_localization() -> None:
    certificate = exact_double_triple_residual_certificate()

    assert certificate.excluded_family_equation_identities == (0, 0, 0, 0)
    assert certificate.excluded_family_same_fiber_identity == 0
    assert certificate.hostile_equation_residuals == (0, 0, 0, 0)
    assert certificate.hostile_residual_value == 0
    assert certificate.hostile_same_fiber_value == 0
    assert certificate.hostile_cusp_values == (2, 2)
    assert certificate.hostile_root_difference == 2
    assert certificate.hostile_coefficient_rank == 3
    assert certificate.hostile_augmented_rank == 3


def test_two_contact_residual_has_only_a_proper_compatible_locus() -> None:
    certificate = exact_double_contact_residual_certificate()

    assert certificate.determinant_factor_remainder == 0
    assert certificate.residual_total_degree == 11
    assert certificate.residual_term_count == 79
    assert certificate.augmented_minor_gcd == 1

    # Generic incompatibility does not mean the residual is empty: retain
    # exact hostile fixtures of both behaviors.
    assert certificate.incompatible_residual_value == 0
    assert certificate.incompatible_base_localizer_value != 0
    assert certificate.incompatible_coefficient_rank == 3
    assert certificate.incompatible_augmented_rank == 4
    assert certificate.compatible_residual_value == 0
    assert certificate.compatible_base_localizer_value != 0
    assert certificate.compatible_coefficient_rank == 3
    assert certificate.compatible_augmented_rank == 3
    assert certificate.compatible_equation_residuals == (0, 0, 0, 0)
    assert certificate.compatible_extra_critical_value == 0


def test_rank_two_curves_are_confined_to_invalid_boundaries() -> None:
    certificate = exact_double_contact_residual_certificate()

    assert certificate.projection_curve_factor_identity == 0
    assert certificate.vertical_projection_groebner_basis == (1,)
    assert certificate.restricted_singular_curve_gcds == (
        1,
        (KAPPA + 2) ** 3,
        (KAPPA - 2) ** 3,
        2 * KAPPA * ROOT_SUM + ROOT_SUM**2 + 4,
    )
    assert certificate.split_curve_boundary_identities == (0, 0, 0, 0)
    assert certificate.diagonal_cusp_identity == 0


def test_rank_one_is_impossible_on_the_valid_contact_chart() -> None:
    certificate = exact_double_contact_residual_certificate()

    assert certificate.jet_division_remainders == (0, 0, 0, 0, 0, 0)
    assert certificate.jet_split_power_remainder == 0


def test_both_residuals_have_no_threefold_incidence_component() -> None:
    certificate = exact_delta_ten_residual_rank_certificate()

    assert certificate.double_triple.verified
    assert certificate.double_contact.maximum_valid_incidence_dimension == 2
    assert certificate.double_contact.verified
    assert certificate.verified
