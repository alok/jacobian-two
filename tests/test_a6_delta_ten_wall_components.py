"""Tests for the dominant delta-ten contact and triple wall components."""

from sympy import Rational, expand

from scripts.a6_delta_ten_generic import (
    PAIR_DENOMINATOR,
    PAIR_DIAGONAL_FACTOR,
    S,
)
from scripts.a6_delta_ten_wall_components import (
    CONTACT_EXPECTED_RESIDUAL_BASIS,
    TRIPLE_BASE_DISCRIMINANT,
    TRIPLE_BASE_POLYNOMIAL,
    exact_delta_ten_wall_components_certificate,
)


def test_contact_incidence_has_only_the_exact_rank_drop_factors() -> None:
    certificate = exact_delta_ten_wall_components_certificate()

    assert certificate.contact_minor_gcd_identity == 0
    assert certificate.contact_residual_groebner_basis == tuple(
        expand(polynomial) for polynomial in CONTACT_EXPECTED_RESIDUAL_BASIS
    )
    assert certificate.contact_localized_groebner_basis == (1,)
    assert S * PAIR_DENOMINATOR * PAIR_DIAGONAL_FACTOR != 0


def test_contact_rank_drop_components_and_dominant_sample_are_exact() -> None:
    certificate = exact_delta_ten_wall_components_certificate()

    assert certificate.contact_s_zero_identities == (0, 0)
    assert certificate.contact_denominator_identities == (0, 0)
    assert certificate.contact_diagonal_identities == (0, 0)
    assert certificate.contact_cusp_boundary_identities == (0, 0)
    assert certificate.contact_sample_ab_minor == 131220
    assert certificate.contact_sample_gcd == S + 2
    assert certificate.contact_sample_second_derivative == Rational(2538, 5)
    assert certificate.contact_sample_tangency_derivative == Rational(30456, 5)


def test_ordinary_triple_base_is_irreducible_and_sample_has_rank_two() -> None:
    certificate = exact_delta_ten_wall_components_certificate()

    assert certificate.triple_base_factorization == TRIPLE_BASE_POLYNOMIAL
    assert certificate.triple_base_discriminant == TRIPLE_BASE_DISCRIMINANT
    assert (
        certificate.triple_base_discriminant_squarefree_part
        == TRIPLE_BASE_DISCRIMINANT
    )
    assert certificate.triple_base_discriminant_p_degree == 1
    assert certificate.triple_sample_base_value == 0
    assert certificate.triple_sample_localizer != 0
    assert certificate.triple_sample_kappa == 2
    assert certificate.triple_sample_ab_minor == Rational(8652, 48828125)


def test_triple_equality_rank_never_drops_on_the_valid_base() -> None:
    certificate = exact_delta_ten_wall_components_certificate()

    assert certificate.triple_coefficient_localized_groebner_basis == (1,)
    assert certificate.triple_augmented_localized_groebner_basis == (1,)
    assert certificate.verified
