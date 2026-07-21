"""Exact fixtures for the constant-leading quadratic-coordinate theorem.

The Lean development is the authoritative proof.  This independent SymPy
module generates a high-degree instance of the announced normal form, checks
its formal Jacobian, and verifies the displayed inverse in both directions.
It also changes one coefficient of ``f(x)`` to demonstrate that the constant
Jacobian identity is genuinely using the affine-discriminant hypothesis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, TypeAlias

from sympy import Expr, Integer, Matrix, Symbol, cancel, expand

PlaneMap: TypeAlias = tuple[Expr, Expr]

X: Final = Symbol("x")
Y: Final = Symbol("y")
U: Final = Symbol("u")
V: Final = Symbol("v")
T: Final = Symbol("t")

EPSILON: Final[Expr] = Integer(3)
DISCRIMINANT_SLOPE: Final[Expr] = Integer(10)
DISCRIMINANT_INTERCEPT: Final[Expr] = Integer(7)
LAMBDA: Final[Expr] = Integer(3)
EXPECTED_DETERMINANT: Final[Expr] = Integer(15)

INNER_POLYNOMIAL: Final[Expr] = X**6 - 2 * X**4 + 3 * X**2 - X + 5
G_COEFFICIENT: Final[Expr] = 4 * INNER_POLYNOMIAL
OUTER_POLYNOMIAL: Final[Expr] = (
    T**7 - 3 * T**5 + 2 * T**3 + 7 * T**2 - 11 * T + 13
)
DISCRIMINANT: Final[Expr] = (
    DISCRIMINANT_SLOPE * X + DISCRIMINANT_INTERCEPT
)
F: Final[Expr] = cancel(
    (G_COEFFICIENT**2 - DISCRIMINANT) / (4 * EPSILON)
)
S: Final[Expr] = 2 * EPSILON * Y + G_COEFFICIENT


def jacobian(polynomial_map: PlaneMap) -> Expr:
    """Return the exact formal Jacobian determinant in ``(x,y)`` order."""

    return cancel(Matrix(polynomial_map).jacobian((X, Y)).det())


def generated_map(*, f_x2_perturbation: int = 0) -> PlaneMap:
    """Generate the normal form, optionally changing one coefficient of ``f``.

    With no perturbation, ``g^2 - 4*eps*f = A*x + B``.  A nonzero argument
    adds that scalar to the coefficient of ``x^2`` in ``f`` and rebuilds the
    same normal form ``P = G_outer(Q) + lambda*s``.
    """

    f = F + Integer(f_x2_perturbation) * X**2
    second = EPSILON * Y**2 + G_COEFFICIENT * Y + f
    first = OUTER_POLYNOMIAL.subs(T, second) + LAMBDA * S
    return first, second


def displayed_inverse() -> PlaneMap:
    """Return the theorem's polynomial inverse in target variables ``(u,v)``."""

    recovered_s = (U - OUTER_POLYNOMIAL.subs(T, V)) / LAMBDA
    recovered_x = (
        recovered_s**2
        - 4 * EPSILON * V
        - DISCRIMINANT_INTERCEPT
    ) / DISCRIMINANT_SLOPE
    recovered_y = (
        recovered_s - G_COEFFICIENT.subs(X, recovered_x)
    ) / (2 * EPSILON)
    return recovered_x, recovered_y


@dataclass(frozen=True)
class ConstantLeadingQuadraticCertificate:
    """Exact identities for one high-degree generated normal form."""

    polynomial_map: PlaneMap
    inverse: PlaneMap
    determinant: Expr
    expected_determinant: Expr
    completed_square_residual: Expr
    inverse_after_map: PlaneMap
    map_after_inverse: PlaneMap

    @property
    def verified(self) -> bool:
        """Whether the determinant, coordinate identity, and inverses agree."""

        return (
            expand(self.determinant - self.expected_determinant) == 0
            and self.completed_square_residual == 0
            and self.inverse_after_map == (X, Y)
            and self.map_after_inverse == (U, V)
        )


def exact_certificate() -> ConstantLeadingQuadraticCertificate:
    """Compute the determinant and both explicit inverse compositions.

    The compositions are simplified in the intermediate ``s`` coordinate.
    For the target-side composition, ``Q`` is evaluated through the separately
    checked identity ``s^2 - 4*eps*Q = A*x + B``.  This is algebraically the
    original second coordinate while avoiding an enormous blind expansion of
    a degree-84 polynomial.
    """

    polynomial_map = generated_map()
    first, second = polynomial_map
    inverse = displayed_inverse()
    inverse_x, inverse_y = inverse

    completed_square_residual = expand(
        S**2 - 4 * EPSILON * second - DISCRIMINANT
    )

    source_s = cancel(
        (first - OUTER_POLYNOMIAL.subs(T, second)) / LAMBDA
    )
    source_x = cancel(
        (
            source_s**2
            - 4 * EPSILON * second
            - DISCRIMINANT_INTERCEPT
        )
        / DISCRIMINANT_SLOPE
    )
    source_y = cancel(
        (source_s - G_COEFFICIENT.subs(X, source_x)) / (2 * EPSILON)
    )

    target_s = cancel(
        2 * EPSILON * inverse_y + G_COEFFICIENT.subs(X, inverse_x)
    )
    target_q = cancel(
        (
            target_s**2
            - (DISCRIMINANT_SLOPE * inverse_x + DISCRIMINANT_INTERCEPT)
        )
        / (4 * EPSILON)
    )
    target_p = cancel(OUTER_POLYNOMIAL.subs(T, target_q) + LAMBDA * target_s)

    return ConstantLeadingQuadraticCertificate(
        polynomial_map=polynomial_map,
        inverse=inverse,
        determinant=jacobian(polynomial_map),
        expected_determinant=EXPECTED_DETERMINANT,
        completed_square_residual=completed_square_residual,
        inverse_after_map=(source_x, source_y),
        map_after_inverse=(target_p, target_q),
    )


def one_coefficient_perturbation_determinant() -> Expr:
    """Return the nonconstant determinant after replacing ``f`` by ``f+x^2``."""

    return jacobian(generated_map(f_x2_perturbation=1))


def main() -> int:
    """Print the exact positive and hostile certificates."""

    certificate = exact_certificate()
    perturbed_determinant = one_coefficient_perturbation_determinant()
    perturbation_rejected = X in perturbed_determinant.free_symbols
    print(f"det J(P,Q) = {certificate.determinant}")
    print(f"completed-square identity verified: "
          f"{certificate.completed_square_residual == 0}")
    print(f"inverse after map: {certificate.inverse_after_map}")
    print(f"map after inverse: {certificate.map_after_inverse}")
    print(f"normal form verified: {certificate.verified}")
    print(f"one-coefficient perturbation determinant: {perturbed_determinant}")
    print(f"hostile perturbation rejected: {perturbation_rejected}")
    return 0 if certificate.verified and perturbation_rejected else 1


if __name__ == "__main__":
    raise SystemExit(main())
