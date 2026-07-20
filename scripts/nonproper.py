"""Exact algebra for the announced map's nonproper-value set.

The topological statement is proved in ``docs/nonproper-set.md``.  This module
independently checks every polynomial and rational identity used by that
proof.  All calculations are symbolic and exact.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, cast

from sympy import (
    Expr,
    Integer,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    groebner,
)

from scripts.verify import X, Y, Z, announced_map

A: Final = Symbol("a")
B: Final = Symbol("b")
C: Final = Symbol("c")
T: Final = Symbol("T")
EPSILON: Final = Symbol("epsilon")
W: Final = Symbol("w")
TAU: Final = Symbol("tau")

FIBER_POLYNOMIAL: Final[Expr] = C * T**3 - 2 * T**2 + B * T - 2 * A
FIBER_DERIVATIVE: Final[Expr] = diff(FIBER_POLYNOMIAL, T)
NONPROPER_POLYNOMIAL: Final[Expr] = (
    27 * A**2 * C**2 - 18 * A * B * C + 16 * A + B**3 * C - B**2
)
HOMOGENIZED_FIBER: Final[Expr] = (
    C * T**3 - 2 * T**2 * W + B * T * W**2 - 2 * A * W**3
)


def _zero_rational_expression(expression: Expr) -> bool:
    """Return whether an exact rational expression simplifies to zero."""

    return bool(cancel(expression) == 0)


def source_from_root(t: Expr, b: Expr, c: Expr) -> tuple[Expr, Expr, Expr]:
    """Reconstruct a source point from a simple root of the fiber equation."""

    r = 3 * c * t**2 - 4 * t + b
    return (
        2 / r,
        t - r / 2,
        Rational(5, 4) * r**2 - Rational(3, 2) * t * r - c * r**3 / 8,
    )


def target_with_root(t: Expr, b: Expr, c: Expr) -> tuple[Expr, Expr, Expr]:
    """Return the unique target with coordinates ``b,c`` for which ``t`` is a root."""

    a = (c * t**3 - 2 * t**2 + b * t) / 2
    return (a, b, c)


def repeated_root_target(t: Expr, c: Expr) -> tuple[Expr, Expr, Expr]:
    """Parametrize the discriminant hypersurface by a repeated finite root."""

    return (t**2 - c * t**3, 4 * t - 3 * c * t**2, c)


def source_with_fiber_derivative(
    t: Expr, r: Expr, c: Expr
) -> tuple[Expr, Expr, Expr]:
    """Use ``t`` and a prescribed fiber derivative as reconstruction data."""

    return (
        2 / r,
        t - r / 2,
        Rational(5, 4) * r**2 - Rational(3, 2) * t * r - c * r**3 / 8,
    )


def evaluate_announced_symbolically(point: tuple[Expr, Expr, Expr]) -> tuple[Expr, Expr, Expr]:
    """Evaluate the announced map at a symbolic rational point."""

    substitution = {X: point[0], Y: point[1], Z: point[2]}
    return tuple(cancel(coordinate.subs(substitution)) for coordinate in announced_map())


def large_root_formulas() -> tuple[Expr, Expr, Expr]:
    """Return reconstruction after eliminating ``c`` with the fiber equation.

    These expressions make the removable ``T = infinity`` chart explicit.
    In particular they tend to ``(0, b, a - 4*b**2)`` for bounded ``a,b``.
    """

    c_from_root = (2 * T**2 - B * T + 2 * A) / T**3
    reconstructed = source_from_root(T, B, c_from_root)
    return tuple(cancel(coordinate) for coordinate in reconstructed)


def singular_locus_groebner_basis() -> tuple[Expr, ...]:
    """Compute a lex basis for ``(Q, dQ/da, dQ/db, dQ/dc)`` exactly."""

    generators = (
        NONPROPER_POLYNOMIAL,
        diff(NONPROPER_POLYNOMIAL, A),
        diff(NONPROPER_POLYNOMIAL, B),
        diff(NONPROPER_POLYNOMIAL, C),
    )
    basis = groebner(generators, A, B, C, order="lex")
    return tuple(cast(Expr, polynomial.as_expr()) for polynomial in basis.polys)


LARGE_ROOT_Z_EXPANSION: Final[Expr] = (
    A
    - 4 * B**2
    + (21 * A * B + 5 * B**3) / T
    - (27 * A**2 + 42 * A * B**2 + B**4) / T**2
    + (117 * A**2 * B + 11 * A * B**3) / T**3
    - (108 * A**3 + 45 * A**2 * B**2) / T**4
    + 81 * A**3 * B / T**5
    - 54 * A**4 / T**6
)


@dataclass(frozen=True)
class NonproperAlgebraCertificate:
    """Exact residuals for the algebraic core of the nonproper-set proof."""

    discriminant_residual: Expr
    discriminant_of_locus_residual: Expr
    reconstruction_residuals: tuple[Expr, Expr, Expr]
    large_x_residual: Expr
    large_y_residual: Expr
    large_z_residual: Expr
    repeated_root_locus_residual: Expr
    perturbation_root_residual: Expr
    perturbation_derivative_residual: Expr
    infinity_image_residuals: tuple[Expr, Expr, Expr]
    infinity_root_derivative_residual: Expr
    triple_root_factor_residual: Expr
    triple_locus_residuals: tuple[Expr, Expr, Expr, Expr]
    escaping_family_residuals: tuple[Expr, Expr, Expr]

    @property
    def verified(self) -> bool:
        """Whether every exact identity has zero residual."""

        residuals = (
            self.discriminant_residual,
            self.discriminant_of_locus_residual,
            *self.reconstruction_residuals,
            self.large_x_residual,
            self.large_y_residual,
            self.large_z_residual,
            self.repeated_root_locus_residual,
            self.perturbation_root_residual,
            self.perturbation_derivative_residual,
            *self.infinity_image_residuals,
            self.infinity_root_derivative_residual,
            self.triple_root_factor_residual,
            *self.triple_locus_residuals,
            *self.escaping_family_residuals,
        )
        return all(_zero_rational_expression(residual) for residual in residuals)


def exact_nonproper_algebra_certificate() -> NonproperAlgebraCertificate:
    """Compute all exact identities used to identify the nonproper-value set."""

    discriminant_residual = expand(
        discriminant(FIBER_POLYNOMIAL, T) + 4 * NONPROPER_POLYNOMIAL
    )
    discriminant_of_locus_residual = expand(
        discriminant(NONPROPER_POLYNOMIAL, A) + 4 * (3 * B * C - 4) ** 3
    )

    target = target_with_root(T, B, C)
    source = source_from_root(T, B, C)
    image = evaluate_announced_symbolically(source)
    reconstruction_residuals = tuple(
        cancel(image_coordinate - target_coordinate)
        for image_coordinate, target_coordinate in zip(image, target, strict=True)
    )

    large_x, large_y, large_z = large_root_formulas()
    large_x_expected = T / (T**2 - B * T + 3 * A)
    large_y_expected = B - 3 * A / T

    repeated_b = 4 * T - 3 * C * T**2
    repeated_a = T**2 - C * T**3
    repeated_root_locus_residual = expand(
        NONPROPER_POLYNOMIAL.subs({A: repeated_a, B: repeated_b})
    )

    perturbed_t = T + EPSILON
    perturbed_target = target_with_root(perturbed_t, repeated_b, C)
    perturbation_root_residual = expand(
        FIBER_POLYNOMIAL.subs(
            {
                A: perturbed_target[0],
                B: repeated_b,
                T: perturbed_t,
            }
        )
    )
    perturbed_derivative = FIBER_DERIVATIVE.subs(
        {B: repeated_b, T: perturbed_t}
    )
    perturbation_derivative_expected = EPSILON * (
        6 * C * T - 4 + 3 * C * EPSILON
    )

    infinity_source = (Integer(0), B, A - 4 * B**2)
    infinity_image = evaluate_announced_symbolically(infinity_source)
    infinity_target = (A, B, Integer(0))
    infinity_image_residuals = tuple(
        expand(image_coordinate - target_coordinate)
        for image_coordinate, target_coordinate in zip(
            infinity_image, infinity_target, strict=True
        )
    )
    infinity_local_polynomial = HOMOGENIZED_FIBER.subs({T: Integer(1)})
    infinity_root_derivative_residual = expand(
        diff(infinity_local_polynomial, W).subs(
            {C: Integer(0), W: Integer(0)}
        )
        + 2
    )

    triple_target = (TAU**2 / 3, 2 * TAU, 2 / (3 * TAU))
    triple_root_factor_residual = cancel(
        FIBER_POLYNOMIAL.subs(
            {A: triple_target[0], B: triple_target[1], C: triple_target[2]}
        )
        - 2 * (T - TAU) ** 3 / (3 * TAU)
    )
    triple_locus_residuals = tuple(
        cancel(
            expression.subs(
                {A: triple_target[0], B: triple_target[1], C: triple_target[2]}
            )
        )
        for expression in (
            NONPROPER_POLYNOMIAL,
            diff(NONPROPER_POLYNOMIAL, A),
            diff(NONPROPER_POLYNOMIAL, B),
            diff(NONPROPER_POLYNOMIAL, C),
        )
    )

    repeated_target = repeated_root_target(T, C)
    escaping_source = source_with_fiber_derivative(T, EPSILON, C)
    escaping_target = (
        repeated_target[0] + EPSILON * T / 2,
        repeated_target[1] + EPSILON,
        C,
    )
    escaping_image = evaluate_announced_symbolically(escaping_source)
    escaping_family_residuals = tuple(
        cancel(image_coordinate - target_coordinate)
        for image_coordinate, target_coordinate in zip(
            escaping_image, escaping_target, strict=True
        )
    )

    return NonproperAlgebraCertificate(
        discriminant_residual=discriminant_residual,
        discriminant_of_locus_residual=discriminant_of_locus_residual,
        reconstruction_residuals=reconstruction_residuals,
        large_x_residual=cancel(large_x - large_x_expected),
        large_y_residual=cancel(large_y - large_y_expected),
        large_z_residual=cancel(large_z - LARGE_ROOT_Z_EXPANSION),
        repeated_root_locus_residual=repeated_root_locus_residual,
        perturbation_root_residual=perturbation_root_residual,
        perturbation_derivative_residual=expand(
            perturbed_derivative - perturbation_derivative_expected
        ),
        infinity_image_residuals=infinity_image_residuals,
        infinity_root_derivative_residual=infinity_root_derivative_residual,
        triple_root_factor_residual=triple_root_factor_residual,
        triple_locus_residuals=triple_locus_residuals,
        escaping_family_residuals=escaping_family_residuals,
    )


def main() -> int:
    """Print the certificate status and fail if an identity does not hold."""

    certificate = exact_nonproper_algebra_certificate()
    print(f"disc_T(p) = {-4 * NONPROPER_POLYNOMIAL}")
    print(f"disc_a(Q) = {-4 * (3 * B * C - 4) ** 3}")
    print(f"large-root source = {large_root_formulas()}")
    print(f"nonproper algebra verified: {certificate.verified}")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
