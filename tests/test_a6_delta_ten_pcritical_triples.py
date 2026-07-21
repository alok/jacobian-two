"""Tests for the delta-ten ``P``-critical triple-boundary audit."""

from sympy import Rational

from scripts.a6_delta_ten_generic import S, T
from scripts.a6_delta_ten_pcritical_triples import (
    CRITICAL_BASE_LOCALIZER,
    CRITICAL_ROOT,
    exact_delta_ten_pcritical_triple_certificate,
)


def test_fourth_root_boundary_is_three_prime_labeled_curves() -> None:
    certificate = exact_delta_ten_pcritical_triple_certificate()

    assert certificate.critical_p_factorization_identity == 0
    assert all(
        all(remainder == 0 for remainder in row)
        for row in certificate.labeled_component_remainders
    )
    assert certificate.labeled_boundary_saturation_exponent == 1
    assert certificate.labeled_boundary_component_count == 3
    assert certificate.labeled_boundary_dimensions == (1, 1, 1)
    assert certificate.labeled_boundary_prime == (True, True, True)
    assert certificate.other_root_discriminant_identity == 0
    assert certificate.other_p_derivative_resultant_identity == 0
    assert CRITICAL_BASE_LOCALIZER.subs(CRITICAL_ROOT, 2) != 0


def test_immersed_critical_t112_boundary_has_dimension_two() -> None:
    certificate = exact_delta_ten_pcritical_triple_certificate()

    assert certificate.t111_minor_gcd_identity == 0
    assert certificate.t111_normalized_minor_gcd == 1
    assert certificate.t112_minor_gcd_identity == 0
    assert certificate.t112_normalized_minor_gcd == 1
    assert certificate.critical_tangent_reduction_identity == 0
    assert certificate.t112_immersion_determinant_identity == 0
    # At e^2+2=0 the immersion row is coefficient-dependent, but the three
    # incidence equations force Q'(e)=100 rather than zero.
    assert certificate.exceptional_immersion_residuals == (0, 0, 0, 100)
    assert certificate.t112_component_dimension == 2
    assert not certificate.t112_same_profile_threefold_exists


def test_critical_c2_t111_boundary_and_true_split_charts_are_lower_dimensional() -> (
    None
):
    certificate = exact_delta_ten_pcritical_triple_certificate()

    assert certificate.c2_incidence_determinant_identity == 0
    assert certificate.c2_rank_two_minor_gcd_identity == 0
    assert certificate.c2_normalized_rank_two_minor_gcd == 1
    assert certificate.pair_denominator_incidence_identity == 0
    assert certificate.c2_component_dimension == 2
    assert not certificate.c2_same_profile_threefold_exists
    assert certificate.split_vertical_determinant_identities == (0, 0)
    # rho=-1/4 is compatible but is the already-tangent triple pair; the two
    # other determinant roots are inconsistent augmented systems.
    assert certificate.split_vertical_rank_data == (
        ((3, 3), (3, 4), (3, 4)),
        ((3, 3), (3, 4), (3, 4)),
    )


def test_hostile_pcritical_t112_member_has_six_nodes() -> None:
    certificate = exact_delta_ten_pcritical_triple_certificate()

    assert certificate.t112_sample_incidence_residuals == (0, 0, 0)
    assert certificate.t112_sample_critical_derivative == (0, 588)
    assert certificate.t112_sample_images == (-4, 8, -4, 8)
    assert certificate.t112_sample_contact_jets == (
        0,
        Rational(23, 343) * (2 * T + 1),
    )
    assert certificate.t112_sample_collision_identity == 0
    assert certificate.t112_sample_collision_gcd == S + 1
    assert certificate.t112_sample_residual_discriminant == 6340872382174352400
    assert certificate.t112_sample_residual_resultants == (
        55641600,
        124085451579933126415376475000000,
    )
    assert certificate.t112_sample_node_x_identity == 0
    assert certificate.t112_sample_node_x_discriminant != 0
    assert certificate.t112_sage_jacobian_components == (
        (4, 1),
        (6, 1),
        (6, 6),
    )


def test_hostile_pcritical_c2_t111_member_has_five_nodes() -> None:
    certificate = exact_delta_ten_pcritical_triple_certificate()

    assert certificate.c2_sample_incidence_residuals == (0, 0, 0, 0)
    assert certificate.c2_sample_critical_derivative == (0, Rational(504, 5))
    assert certificate.c2_sample_images == (
        -4,
        Rational(272, 5),
        2,
        Rational(89, 5),
    )
    assert certificate.c2_sample_triple_slope_difference == (
        Rational(-58, 35) * (2 * T + 1)
    )
    assert certificate.c2_sample_contact_jets == (
        0,
        Rational(-14, 845) * (2 * T - 1),
    )
    assert certificate.c2_sample_collision_identity == 0
    assert certificate.c2_sample_collision_gcd == S - 1
    assert certificate.c2_sample_residual_discriminant == -697306407293437500
    assert certificate.c2_sample_residual_resultants == (
        7,
        225504000,
        -8638634891348924179843125000000,
    )
    assert certificate.c2_sample_node_x_identity == 0
    assert certificate.c2_sample_node_x_discriminant != 0
    assert certificate.c2_sage_jacobian_components == (
        (4, 1),
        (3, 1),
        (4, 1),
        (5, 5),
    )


def test_topology_is_an_explicit_uncomputed_dependency() -> None:
    certificate = exact_delta_ten_pcritical_triple_certificate()

    assert not certificate.topology_computed
    assert certificate.topology_propagation_dependencies == (
        "connected clean open",
        "proper projective Whitney-Thom triviality",
    )
    assert certificate.verified
