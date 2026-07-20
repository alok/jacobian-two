"""Positive and adversarial checks for the independent SymPy verifier."""

from sympy import Integer, Rational

from scripts.verify import DISPLAYED_POINTS, TARGET, exact_certificate


def test_exact_screenshot_certificate() -> None:
    certificate = exact_certificate()

    assert certificate.determinant == -2
    assert certificate.images == (TARGET, TARGET, TARGET)
    assert certificate.verified


def test_perturbed_coefficient_breaks_constant_determinant() -> None:
    certificate = exact_certificate(displayed_four=5)

    assert certificate.determinant != -2
    assert not certificate.determinant_verified
    assert not certificate.verified


def test_perturbed_collision_point_breaks_fiber_certificate() -> None:
    bad_points = (
        DISPLAYED_POINTS[0],
        DISPLAYED_POINTS[1],
        (Integer(-1), Rational(3, 2), Rational(15, 2)),
    )
    certificate = exact_certificate(points=bad_points)

    assert certificate.determinant_verified
    assert certificate.images[-1] != TARGET
    assert not certificate.collision_verified
    assert not certificate.verified

