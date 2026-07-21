"""Tests for the conditional delta-ten ``C3 + 7N`` stratum."""

from sympy import Rational, expand

from scripts.a6_delta_ten_contact_three import (
    CONTACT_THREE_EXPECTED_RANKDROP_BASIS,
    CONTACT_THREE_RELATIONS,
    exact_delta_ten_contact_three_certificate,
)
from scripts.a6_delta_ten_generic import S, T


def test_contact_three_incidence_has_one_valid_dominant_component() -> None:
    certificate = exact_delta_ten_contact_three_certificate()

    assert certificate.incidence_minor_gcd_identity == 0
    assert (
        certificate.incidence_rankdrop_localized_basis
        == CONTACT_THREE_EXPECTED_RANKDROP_BASIS
    )
    assert certificate.incidence_augmented_inconsistency == Rational(
        148635648,
        625,
    )
    assert certificate.incidence_sample_rank_minor == -2007666
    assert all(identity == 0 for identity in certificate.boundary_identities)


def test_rational_member_has_one_exact_contact_three_root() -> None:
    certificate = exact_delta_ten_contact_three_certificate()

    assert certificate.collision_identity == 0
    assert certificate.tangency_identity == 0
    assert certificate.collision_gcd == expand((S + 2) ** 2)
    assert certificate.triple_root_gcd == S + 2
    assert certificate.third_derivative == Rational(-3069, 4)
    assert certificate.tangency_second_derivative == -9207
    assert certificate.residual_discriminant == 168219817614587461632
    assert certificate.residual_contact_separation == -2046
    assert certificate.cusp_image_factor == Rational(273, 256)
    assert certificate.extra_critical_factor == 21463
    assert certificate.denominator_resultant == 81
    assert certificate.diagonal_resultant == 5215509
    assert certificate.residual_tangency_resultant == 62222061329145233117755350912
    assert certificate.contact_chart_values == (-3, 3, 6)
    assert certificate.contact_pair_discriminant == -4


def test_local_jets_prove_intersection_multiplicity_three() -> None:
    certificate = exact_delta_ten_contact_three_certificate()

    assert certificate.contact_image_remainders == (-2, -464)
    assert certificate.contact_p_derivative_resultant == 52
    assert certificate.contact_jet_values == (448, -208)
    assert certificate.contact_jet_differences == (
        0,
        0,
        -32736 * (T + 1) / 2197,
    )
    assert certificate.contact_third_jet_separation == 1


def test_seven_nodes_and_projective_genus_budget_are_exact() -> None:
    certificate = exact_delta_ten_contact_three_certificate()

    assert certificate.node_x_resultant_identity == 0
    assert certificate.node_x_discriminant != 0
    assert certificate.contact_node_separation == 83532198178
    assert certificate.cusp_node_separation == -8496306
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_components == ((4, 1), (5, 1), (7, 7))
    assert certificate.total_delta == 28
    assert certificate.total_delta == certificate.arithmetic_genus


def test_cyclic_complement_has_no_required_a6_image() -> None:
    certificate = exact_delta_ten_contact_three_certificate()

    assert len(CONTACT_THREE_RELATIONS) == 11
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0
    assert certificate.verified
