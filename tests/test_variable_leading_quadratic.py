"""Exact regression tests for the variable-leading quadratic proof target."""

import pytest
from sympy import Integer, Poly, Rational, expand

from scripts.variable_leading_quadratic import (
    KAPPA,
    X,
    Y,
    central_binomial_certificate,
    gradient_unit_ideal_certificate,
    leading_equation_certificate,
    main,
    rational_mate_certificate,
    truncated_recurrence_certificate,
)


def test_rational_jacobian_one_mate_keeps_a_denominator() -> None:
    certificate = rational_mate_certificate()

    assert certificate.determinant == 1
    assert certificate.denominator == 1 + X**2 * Y**2
    assert certificate.denominator != 1
    assert certificate.verified


def test_hostile_quadratic_gradient_generates_the_unit_ideal() -> None:
    certificate = gradient_unit_ideal_certificate()

    assert certificate.derivative_x == 1 + 2 * X * Y**2
    assert certificate.derivative_y == 2 * X**2 * Y
    assert certificate.multiplier_x == 1 - 2 * X * Y**2
    assert certificate.multiplier_y == 2 * Y**3
    assert certificate.bezout_combination == 1
    assert certificate.verified


def test_nonconstant_square_solves_only_the_top_equation() -> None:
    for degree in (1, 2, 3, 8, 17):
        certificate = leading_equation_certificate(degree)

        assert certificate.leading_quadratic_coefficient == X**2
        assert certificate.leading_mate_coefficient == X**degree
        assert certificate.residual == 0
        assert certificate.verified


def test_central_binomial_product_is_exact_at_multiple_depths() -> None:
    expected = (
        Integer(1),
        Rational(1, 2),
        Rational(3, 8),
        Rational(5, 16),
    )
    for depth, expected_value in enumerate(expected):
        certificate = central_binomial_certificate(depth)

        assert certificate.recurrence_product == expected_value
        assert certificate.central_binomial_form == expected_value
        assert certificate.residual == 0
        assert certificate.verified

    assert central_binomial_certificate(12).verified


def test_every_finite_forced_tail_has_a_nonzero_top_residual() -> None:
    for depth in (0, 1, 2, 5, 9):
        certificate = truncated_recurrence_certificate(depth)

        assert len(certificate.coefficients) == depth + 1
        assert len(certificate.recurrence_residuals) == depth
        assert all(residual == 0
                   for residual in certificate.recurrence_residuals)
        for index in range(1, depth + 1):
            expected_coefficient = expand(
                -Rational(2 * index, 2 * index + 1)
                * X
                * certificate.coefficients[index - 1]
            )
            assert certificate.coefficients[index] == expected_coefficient

        assert certificate.top_residual == certificate.expected_top_residual
        assert certificate.top_residual != 0
        assert expand(
            certificate.determinant
            - KAPPA
            - certificate.expected_top_residual
        ) == 0

        residual_polynomial = Poly(certificate.top_residual, X, Y)
        assert residual_polynomial.monoms() == [
            (depth + 1, 2 * depth + 2),
        ]
        assert certificate.verified


def test_invalid_depths_and_zero_keller_constant_are_rejected() -> None:
    with pytest.raises(ValueError, match="positive"):
        leading_equation_certificate(0)
    with pytest.raises(ValueError, match="nonnegative"):
        central_binomial_certificate(-1)
    with pytest.raises(ValueError, match="nonnegative"):
        truncated_recurrence_certificate(-1)
    with pytest.raises(ValueError, match="nonzero"):
        truncated_recurrence_certificate(2, kappa=Integer(0))


def test_command_line_depth_runs_all_hostile_checks(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--depth", "4"]) == 0
    output = capsys.readouterr().out
    assert "forced recurrence depth: 4" in output
    assert "all hostile identities verified: True" in output
