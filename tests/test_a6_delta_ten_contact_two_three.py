"""Tests for the conditional delta-ten ``C2 + C3 + 5N`` profile."""

from sympy import Rational, expand

from scripts.a6_delta_ten_contact_two_three import (
    CONTACT_TWO_THREE_RELATIONS,
    exact_delta_ten_contact_two_three_certificate,
)
from scripts.a6_delta_ten_generic import S, T


def test_compatibility_surface_has_a_smooth_rational_rank_four_point() -> None:
    certificate = exact_delta_ten_contact_two_three_certificate()

    assert certificate.target_separation_identity == 0
    assert certificate.compatibility_factor_exponents == (-2, 6, 1, 3)
    assert certificate.compatibility_residual_shape == (409, 18, 13, 9, 10)
    assert certificate.compatibility_rational_factor_count == 1
    assert certificate.compatibility_sample_value == 0
    assert certificate.compatibility_sample_gradient == (
        Rational(1073741824, 59049),
        Rational(536870912, 59049),
        Rational(1073741824, 177147),
    )
    assert certificate.compatibility_sample_is_smooth
    assert certificate.compatibility_surface_geometrically_irreducible
    assert certificate.sage_coefficient_rank_three_curve == (1, 30, -51, 4, 0)
    assert certificate.sage_augmented_rank_three_unit == (True, 3)
    assert certificate.sage_rank_two_unit_exponents == (3, 3)
    assert certificate.sage_valid_singular_saturation == (True, 2)
    assert certificate.valid_rank_drop_incidence_empty
    assert certificate.valid_compatibility_surface_smooth
    assert certificate.sample_coefficient_rank == 4
    assert certificate.sample_augmented_rank == 4
    assert certificate.sample_incidence_residuals == (0, 0, 0, 0, 0)
    assert all(value != 0 for value in certificate.sample_maximal_minors)
    assert certificate.sample_valid_localizer != 0
    assert certificate.dominant_incidence_dimension == 2


def test_sample_collision_has_exact_contact_two_and_contact_three() -> None:
    certificate = exact_delta_ten_contact_two_three_certificate()

    assert certificate.collision_identity == 0
    assert certificate.tangency_identity == 0
    assert certificate.collision_tangency_gcd == expand((S - 1) * (S + 1) ** 2)
    assert certificate.double_root_support == expand((S - 1) * (S + 1) ** 2)
    assert certificate.triple_root_support == S + 1
    assert certificate.collision_jets_at_contact_two == (
        0,
        0,
        Rational(1024, 3),
        Rational(11648, 3),
    )
    assert certificate.collision_jets_at_contact_three == (
        0,
        0,
        0,
        Rational(-64, 3),
        Rational(640, 3),
    )
    assert certificate.tangency_jets_at_contact_two == (
        0,
        Rational(8192, 3),
        37888,
    )
    assert certificate.tangency_jets_at_contact_three == (
        0,
        0,
        Rational(-256, 3),
        1024,
    )
    assert certificate.residual_discriminant == 256979755008
    assert certificate.residual_contact_separations == (192, -8)


def test_sample_contact_pairs_are_clean_and_locally_exact() -> None:
    certificate = exact_delta_ten_contact_two_three_certificate()

    assert certificate.cusp_image_factor == Rational(8, 27)
    assert certificate.extra_critical_factor == 9072
    assert certificate.residual_boundary_resultants == (
        Rational(-1024, 27),
        Rational(2048, 9),
        28672,
        -42392524029926834176,
    )
    assert certificate.contact_chart_values == (
        (Rational(8, 3), Rational(8, 3), 8),
        (Rational(-4, 3), Rational(4, 3), 4),
    )
    assert certificate.pair_disjoint_resultant == Rational(4096, 81)
    assert certificate.pair_target_separation_numerator == Rational(-4096, 243)
    assert certificate.contact_pair_discriminants == (-3, -3)
    assert certificate.contact_image_remainders == ((-5, 3), (-1, -1))
    assert certificate.contact_p_derivative_resultants == (252, 36)
    assert certificate.contact_jet_differences == (
        (0, Rational(3, 49) * (2 * T - 1)),
        (0, 0, -T / 2 - Rational(1, 4)),
    )
    assert certificate.contact_jet_separations == (
        Rational(27, 2401),
        Rational(3, 16),
    )


def test_five_nodes_and_projective_genus_budget_are_exact() -> None:
    certificate = exact_delta_ten_contact_two_three_certificate()

    assert certificate.node_x_resultant_identity == 0
    assert certificate.node_x_discriminant == (
        74803083224844271333355450007552000000
    )
    assert certificate.node_target_separations == (16313472, 190400, 50112)
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_length == 17
    assert certificate.sage_jacobian_radical_length == 8
    assert certificate.sage_jacobian_components == (
        (5, 5),
        (5, 1),
        (3, 1),
        (4, 1),
    )
    assert certificate.total_delta == certificate.arithmetic_genus == 28


def test_cyclic_complement_has_no_required_a6_image() -> None:
    certificate = exact_delta_ten_contact_two_three_certificate()

    assert len(CONTACT_TWO_THREE_RELATIONS) == 10
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
