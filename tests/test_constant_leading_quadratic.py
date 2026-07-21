"""Exact positive and hostile checks for Milestone 6."""

from sympy import Poly, expand

from scripts.constant_leading_quadratic import (
    G_COEFFICIENT,
    OUTER_POLYNOMIAL,
    T,
    U,
    V,
    X,
    Y,
    exact_certificate,
    one_coefficient_perturbation_determinant,
)


def test_high_degree_normal_form_and_explicit_inverse() -> None:
    certificate = exact_certificate()

    assert Poly(G_COEFFICIENT, X).degree() == 6
    assert Poly(OUTER_POLYNOMIAL, T).degree() == 7
    assert Poly(certificate.polynomial_map[0], X, Y).total_degree() == 84
    assert certificate.determinant == 15
    assert certificate.completed_square_residual == 0
    assert certificate.inverse_after_map == (X, Y)
    assert certificate.map_after_inverse == (U, V)
    assert certificate.verified


def test_one_coefficient_perturbation_is_not_keller() -> None:
    determinant = one_coefficient_perturbation_determinant()

    assert expand(determinant - (15 - 36 * X)) == 0
    assert X in determinant.free_symbols
