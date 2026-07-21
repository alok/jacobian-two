"""Tests for the exact conditional delta-ten two-contact algebra."""

from sympy import Rational

from scripts.a6_delta_ten_double_contact import (
    ORDERED_EXPECTED_COMPATIBILITY_GCD,
    SAMPLE_CONTACT_SOURCE_QUADRATICS,
    SAMPLE_CONTACT_TARGETS,
    SAMPLE_IMPLICIT,
    SAMPLE_NODE_X_POLYNOMIAL,
    SAMPLE_PARAMETERS,
    SAMPLE_RELATIONS,
    exact_coincident_root_hostile_certificate,
    exact_double_contact_sample_certificate,
    exact_ordered_incidence_certificate,
    exact_split_boundary_certificates,
    ordered_compatibility_gcd,
)


def test_ordered_incidence_has_one_exact_swap_invariant_cramer_chart() -> None:
    certificate = exact_ordered_incidence_certificate()

    assert certificate.determinant_factor_remainder == 0
    assert certificate.residual_rank_swap_identity == 0
    assert certificate.cramer_swap_identities == (0, 0, 0, 0)
    assert certificate.compatibility_gcd_identity == 0
    assert ordered_compatibility_gcd() == ORDERED_EXPECTED_COMPATIBILITY_GCD
    assert certificate.sample_determinant != 0
    assert certificate.sample_parameter_remainders == (0, 0, 0, 0)
    assert certificate.verified


def test_visible_cusp_and_split_denominator_boundaries_are_exact() -> None:
    certificate = exact_ordered_incidence_certificate()

    assert certificate.denominator_root_identity == 0
    assert certificate.cusp_boundary_resultant_identity == 0
    assert certificate.critical_boundary_resultant_identity == 0
    assert certificate.tangency_syzygy_identity == 0


def test_split_fibers_have_no_common_or_nonreduced_discriminant_divisor() -> None:
    certificates = exact_split_boundary_certificates()

    assert tuple(certificate.kappa for certificate in certificates) == (0, 2, -2)
    assert tuple(
        (
            certificate.vertical_discriminant_degree,
            certificate.graph_discriminant_degree,
        )
        for certificate in certificates
    ) == ((2, 13), (6, 9), (6, 9))
    assert all(
        certificate.vertical_discriminant_squarefree for certificate in certificates
    )
    assert all(
        certificate.graph_discriminant_squarefree for certificate in certificates
    )
    assert all(
        certificate.common_discriminant_degree == 0 for certificate in certificates
    )
    assert all(certificate.verified for certificate in certificates)


def test_rational_member_has_two_exact_double_collision_roots() -> None:
    certificate = exact_double_contact_sample_certificate()

    assert SAMPLE_PARAMETERS
    assert certificate.collision_identity == 0
    assert certificate.repeated_gcd_identity == 0
    assert certificate.residual_discriminant == Rational(350724573, 64)
    assert certificate.residual_contact_separation == Rational(-405, 4)
    assert certificate.contact_second_derivatives == (-12, 540)


def test_rational_member_has_two_contact_two_fibers() -> None:
    certificate = exact_double_contact_sample_certificate()

    assert SAMPLE_CONTACT_SOURCE_QUADRATICS
    assert SAMPLE_CONTACT_TARGETS == ((0, Rational(-1, 2)), (-2, Rational(7, 2)))
    assert certificate.tangency_factor_identity == 0
    assert certificate.tangency_gcd_identity == 0
    assert certificate.tangency_gcd_degree == 2
    assert certificate.contact_tangency_cofactor_values == (18, 2430)
    assert certificate.residual_tangency_resultant == Rational(
        -181729682354700675,
        256,
    )
    assert certificate.contact_pair_products == (1, 1)
    assert certificate.contact_pair_discriminants == (-3, -3)
    assert certificate.contact_source_quadratic_identities == (0, 0)
    assert certificate.contact_q_remainders == (0, 0)
    assert certificate.contact_x_values == (0, -2)


def test_rational_member_has_six_separated_transverse_nodes() -> None:
    certificate = exact_double_contact_sample_certificate()

    assert SAMPLE_NODE_X_POLYNOMIAL.as_poly().degree() == 6
    assert certificate.denominator_resultant == -27
    assert certificate.pair_diagonal_resultant == -85293
    assert certificate.cusp_image_factor == Rational(1, 4)
    assert certificate.extra_critical_factor == 9477
    assert certificate.support_x_resultant_identity == 0
    assert certificate.support_x_discriminant == Rational(
        3**99 * 5**24 * 13**2 * 1443311,
        2**52,
    )
    assert SAMPLE_IMPLICIT
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_components == (
        (4, 1),
        (3, 1),
        (3, 1),
        (6, 6),
    )
    assert certificate.node_count == 6
    assert certificate.total_delta == 28
    assert certificate.total_delta == certificate.arithmetic_genus
    assert certificate.verified


def test_raw_presentation_is_cyclic_and_has_no_a6_three_cycle_image() -> None:
    certificate = exact_double_contact_sample_certificate()

    assert len(SAMPLE_RELATIONS) == 11
    assert tuple(len(relation) for relation in SAMPLE_RELATIONS) == (
        4,
        8,
        10,
        26,
        32,
        40,
        60,
        44,
        28,
        12,
        8,
    )
    assert certificate.relation_count == 11
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.complement_census.assignments == 40**4
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0


def test_coincident_root_hostile_fixture_requires_diagonal_saturation() -> None:
    certificate = exact_coincident_root_hostile_certificate()

    assert certificate.duplicated_incidence_values == (0, 0, 0, 0)
    assert certificate.repeated_gcd_identity == 0
    assert certificate.repeated_gcd_degree == 1
    assert certificate.verified
