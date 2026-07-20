"""Positive and adversarial tests for the nonproper-set algebra."""

from sympy import Integer, Rational, cancel, discriminant, expand

from scripts.nonproper import (
    A,
    B,
    C,
    FIBER_POLYNOMIAL,
    NONPROPER_POLYNOMIAL,
    T,
    evaluate_announced_symbolically,
    exact_nonproper_algebra_certificate,
    singular_locus_groebner_basis,
    source_from_root,
    target_with_root,
)


def test_exact_nonproper_algebra_certificate() -> None:
    certificate = exact_nonproper_algebra_certificate()

    assert certificate.verified
    assert certificate.discriminant_of_locus_residual == 0


def test_zero_target_has_an_exact_escaping_family() -> None:
    """At ``epsilon=1/n``, the first source coordinate is exactly ``-n/2``."""

    n = Integer(11)
    epsilon = 1 / n
    t = epsilon
    source = source_from_root(t, Integer(0), Integer(0))
    target = target_with_root(t, Integer(0), Integer(0))

    assert source[0] == -n / 2
    assert target == (-epsilon**2, Integer(0), Integer(0))
    assert evaluate_announced_symbolically(source) == target


def test_off_locus_target_has_no_repeated_fiber_root() -> None:
    off_locus = {A: Integer(1), B: Integer(0), C: Integer(0)}

    assert NONPROPER_POLYNOMIAL.subs(off_locus) == 16
    assert discriminant(FIBER_POLYNOMIAL.subs(off_locus), T) != 0


def test_sign_perturbation_breaks_discriminant_identity() -> None:
    wrong_locus = NONPROPER_POLYNOMIAL + 2 * B**2
    residual = cancel(discriminant(FIBER_POLYNOMIAL, T) + 4 * wrong_locus)

    assert residual != 0
    assert B in residual.free_symbols


def test_singular_locus_eliminates_to_the_triple_root_curve() -> None:
    basis = singular_locus_groebner_basis()

    assert basis == (16 * A - B**3 * C, expand((3 * B * C - 4) ** 2))


def test_announced_collision_target_is_off_the_nonproper_hypersurface() -> None:
    announced_target = {A: Rational(-1, 4), B: Integer(0), C: Integer(0)}

    assert NONPROPER_POLYNOMIAL.subs(announced_target) == -4
