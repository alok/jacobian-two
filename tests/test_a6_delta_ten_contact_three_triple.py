"""Tests for the conditional delta-ten ``C3 + T111 + 4N`` profile."""

from sympy import Rational, expand

from scripts.a6_delta_ten_contact_three_triple import (
    CONTACT_THREE_TRIPLE_RELATIONS,
    exact_delta_ten_contact_three_triple_certificate,
)
from scripts.a6_delta_ten_generic import S, T


def test_compatibility_surface_is_smooth_irreducible_and_rank_four() -> None:
    certificate = exact_delta_ten_contact_three_triple_certificate()

    assert certificate.source_overlap_identity == 0
    assert certificate.target_separation_identity == 0
    assert certificate.determinant_factor_data == (2, 2, 3, 3, 1)
    assert certificate.compatibility_residual_shape == (62, 9, 6, 3, 7)
    assert certificate.compatibility_rational_factor_count == 1
    assert certificate.compatibility_sample_value == 0
    assert certificate.compatibility_sample_gradient == (-84, -112, 0)
    assert certificate.compatibility_sample_is_smooth
    assert certificate.compatibility_surface_geometrically_irreducible
    assert certificate.sage_valid_singular_saturation == (True, 1)
    assert certificate.valid_compatibility_surface_smooth
    assert certificate.sage_coefficient_rank_three_curve == (1, 14, -21, 1)
    assert certificate.sage_augmented_rank_three_unit == (True, 1)
    assert certificate.sage_rank_two_unit_exponents == (1, 1)
    assert certificate.sage_normalized_minor_counts == (5, 25, 40, 100)
    assert certificate.valid_rank_drop_incidence_empty
    assert certificate.sample_coefficient_rank == 4
    assert certificate.sample_augmented_rank == 4
    assert certificate.sample_incidence_residuals == (0, 0, 0, 0, 0)
    assert certificate.sample_maximal_minors == (
        -2504320,
        3333120,
        0,
        Rational(29575, 4),
        Rational(159705, 4),
    )
    assert certificate.sample_valid_localizer == 7223580


def test_two_exact_tangents_prove_a_surface_in_coefficient_space() -> None:
    certificate = exact_delta_ten_contact_three_triple_certificate()

    assert certificate.tangent_incidence_residuals == (
        (0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0),
    )
    assert certificate.coefficient_image_tangent_rank == 2
    assert certificate.dominant_incidence_dimension == 2


def test_sample_has_an_ordinary_triple_and_exact_contact_three() -> None:
    certificate = exact_delta_ten_contact_three_triple_certificate()

    assert certificate.triple_cubic_discriminant == Rational(-637, 4)
    assert certificate.triple_p_derivative_resultant == Rational(5733, 8)
    assert certificate.triple_omitted_root_separation == Rational(9, 2)
    assert certificate.triple_q_remainder == Rational(-2197, 1024)
    assert certificate.triple_omitted_q_difference == Rational(63, 128)
    assert certificate.triple_quotient_identity == 0
    assert certificate.triple_slope_eliminant_identity == 0
    assert certificate.triple_slope_discriminant == (
        -89062908728555597524369408000000
    )
    assert certificate.contact_pair_discriminant == -3
    assert certificate.contact_image_remainders == (3, Rational(-199, 8))
    assert certificate.contact_p_derivative_resultant == 84
    assert certificate.contact_jet_differences == (
        0,
        0,
        Rational(165, 5488) * T - Rational(165, 10976),
    )
    assert certificate.contact_jet_separation == Rational(81675, 120472576)
    assert certificate.collision_identity == 0
    assert certificate.tangency_identity == 0
    assert certificate.collision_tangency_gcd == expand((S - 1) ** 2)
    assert certificate.collision_jets_at_contact == (0, 0, 0, 660, 10356)
    assert certificate.tangency_jets_at_contact == (0, 0, -3960, -96456)


def test_four_nodes_and_projective_genus_budget_are_exact() -> None:
    certificate = exact_delta_ten_contact_three_triple_certificate()

    assert certificate.collision_factor_discriminants == (-2548, -156754683)
    assert certificate.collision_factor_resultants == (7292322, -5, 44)
    assert certificate.triple_boundary_resultants == (
        -96,
        -39,
        -2352,
        -87404254689024000,
    )
    assert certificate.node_boundary_resultants == (
        864,
        -1242,
        -146016,
        -14836390085591121857052672,
    )
    assert certificate.cusp_image_factor == Rational(-299, 64)
    assert certificate.extra_critical_factor == -894348
    assert certificate.triple_x_eliminant_identity == 0
    assert certificate.node_x_eliminant_identity == 0
    assert certificate.node_x_discriminant == -4321892177659568
    assert certificate.node_target_separations == (
        -60306,
        792,
        Rational(-135043, 4096),
    )
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_length == 17
    assert certificate.sage_jacobian_radical_length == 7
    assert certificate.sage_jacobian_components == (
        (4, 4),
        (4, 1),
        (5, 1),
        (4, 1),
    )
    assert certificate.total_delta == certificate.arithmetic_genus == 28


def test_cyclic_complement_has_no_required_a6_image() -> None:
    certificate = exact_delta_ten_contact_three_triple_certificate()

    assert CONTACT_THREE_TRIPLE_RELATIONS
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0
    assert certificate.topology_propagation_dependencies == (
        "connected clean open",
        "proper projective Whitney-Thom triviality",
    )
    assert not certificate.split_and_deeper_boundaries_classified_here
    assert certificate.verified
