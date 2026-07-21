"""Exact positive and hostile tests for the affine-coordinate normal forms."""

from scripts.affine_coordinate import (
    U,
    V,
    X,
    Y,
    nonconstant_slope_determinant,
    nonzero_slope_certificate,
    zero_slope_certificate,
)


def test_arbitrary_degree_nonzero_slope_normal_form() -> None:
    certificate = nonzero_slope_certificate()

    assert certificate.determinant == 15
    assert certificate.inverse_after_map == (X, Y)
    assert certificate.map_after_inverse == (U, V)
    assert certificate.verified


def test_zero_slope_triangular_normal_form() -> None:
    certificate = zero_slope_certificate()

    assert certificate.determinant == 28
    assert certificate.inverse_after_map == (X, Y)
    assert certificate.map_after_inverse == (U, V)
    assert certificate.verified


def test_nonconstant_slope_perturbation_is_not_keller() -> None:
    determinant = nonconstant_slope_determinant()

    assert determinant == 5 * X + 15
    assert X in determinant.free_symbols
