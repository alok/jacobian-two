"""Positive and adversarial tests for Orevkov's degree-five cusp germ."""

from sympy import Rational

from scripts.a6_cusp_germ import (
    A6_CUSP_PARAMETERS,
    TARGET_U,
    TARGET_V,
    X,
    Y,
    exact_a6_cusp_germ_certificate,
)


def test_orevkov_case_b_parameters_match_the_a6_cusp_block() -> None:
    parameters = A6_CUSP_PARAMETERS

    assert parameters.satisfies_hypotheses
    assert (parameters.m1, parameters.m2) == (1, 5)
    assert (
        parameters.d1,
        parameters.d2,
        parameters.degree,
        parameters.ramification_order,
    ) == (2, 5, 5, 3)


def test_a6_cusp_germ_certificate_is_exact() -> None:
    certificate = exact_a6_cusp_germ_certificate()
    u, v = certificate.polynomial_map

    assert u == 3 * X**5 / 8 - 5 * X**3 * Y / 4 + 15 * X * Y**2 / 8
    assert v == Y
    assert certificate.jacobian == 15 * (X**2 - Y) ** 2 / 8
    assert certificate.critical_image_u == X**5
    assert certificate.critical_image_v == X**2
    assert certificate.cusp_pullback == (
        (X**2 - Y) ** 3
        * (9 * X**4 - 33 * X**2 * Y + 64 * Y**2)
        / 64
    )
    assert certificate.residual_image_relation_remainder == 0
    assert certificate.fiber_discriminant == (
        1_036_800_000 * (TARGET_U**2 - TARGET_V**5) ** 2
    )
    assert certificate.verified


def test_mixed_coefficient_perturbation_breaks_the_local_model() -> None:
    certificate = exact_a6_cusp_germ_certificate(
        mixed_coefficient=Rational(-1, 2)
    )

    assert certificate.finite_degree_verified
    assert not certificate.ramification_verified
    assert not certificate.critical_component_verified
    assert not certificate.full_pullback_verified
    assert not certificate.verified
