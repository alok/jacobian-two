"""Tests for the conditional ``A6`` delta-five residual exhaustion."""

from sympy import Rational, expand

from scripts.a6_delta_five_residual import (
    A_SAMPLE,
    B_SAMPLE,
    E_ALPHA,
    E_BETA,
    E_GAMMA,
    E_PARAMETER,
    E_SAMPLE,
    N_ALPHA,
    N_BETA,
    N_GAMMA,
    N_PARAMETER,
    N_SAMPLE,
    P32_MINUS_H,
    P32_PLUS_H,
    P4_PARAMETERS,
    PNT_PARAMETERS,
    exact_a6_delta_five_residual_certificate,
    exact_delta_six_seven_certificate,
)


def test_four_residual_curves_are_exact_and_valid() -> None:
    certificate = exact_a6_delta_five_residual_certificate()

    assert len(certificate.stratum_identities) == 16
    assert all(identity == 0 for identity in certificate.stratum_identities)
    assert len(certificate.residual_ideal_identities) == 12
    assert all(
        identity == 0 for identity in certificate.residual_ideal_identities
    )
    assert certificate.generic_valid_products == (
        Rational(768, 5),
        Rational(1, 1296),
        Rational(250880, 9),
        Rational(-27900),
    )
    assert A_SAMPLE != B_SAMPLE
    assert N_SAMPLE != E_SAMPLE
    assert tuple(
        value.subs(N_PARAMETER, -3)
        for value in (N_ALPHA, N_BETA, N_GAMMA)
    ) == N_SAMPLE
    assert tuple(
        value.subs(E_PARAMETER, 0)
        for value in (E_ALPHA, E_BETA, E_GAMMA)
    ) == E_SAMPLE


def test_exceptional_factor_audit_is_exhaustive() -> None:
    certificate = exact_a6_delta_five_residual_certificate()

    assert len(certificate.exhaustion_identities) == 14
    assert all(identity == 0 for identity in certificate.exhaustion_identities)
    assert len(certificate.exceptional_identities) == 36
    assert all(identity == 0 for identity in certificate.exceptional_identities)
    assert all(value != 0 for value in certificate.exceptional_valid_products)
    assert P4_PARAMETERS != PNT_PARAMETERS
    assert expand(3 * P32_PLUS_H**2 + 12 * P32_PLUS_H + 4) == 0
    assert expand(3 * P32_MINUS_H**2 + 12 * P32_MINUS_H + 4) == 0


def test_every_residual_presentation_has_only_cyclic_three_cycle_images() -> None:
    certificate = exact_a6_delta_five_residual_certificate()

    assert len(certificate.presentations) == 8
    for presentation in certificate.presentations:
        assert presentation.census.assignments == 40**3
        assert presentation.census.satisfying_assignments == 40
        assert presentation.census.generated_order_histogram == ((3, 40),)
        assert presentation.census.a6_assignments == 0
        assert presentation.verified


def test_delta_six_dies_and_delta_seven_has_one_link_survivor() -> None:
    certificate = exact_delta_six_seven_certificate()

    assert tuple(
        candidate.affine_degrees for candidate in certificate.delta_six_candidates
    ) == ((2, 17),)
    assert tuple(
        candidate.affine_degrees for candidate in certificate.delta_seven_candidates
    ) == ((2, 19), (3, 10), (4, 7))
    assert certificate.t2_17.single_three_cycle_a6_pairs == 0
    assert certificate.t2_19.single_three_cycle_a6_pairs == 0
    assert certificate.t3_10.single_three_cycle_a6_pairs == 720
    assert certificate.t4_7.single_three_cycle_a6_pairs == 0
    assert certificate.verified


def test_full_residual_certificate_is_verified() -> None:
    assert exact_a6_delta_five_residual_certificate().verified
