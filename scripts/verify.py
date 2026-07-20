"""Independently verify the announced polynomial map with exact arithmetic.

Lean is the authoritative certificate in this repository.  This module is a
small, separately implemented SymPy checker intended to catch transcription
and formalization mistakes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, TypeAlias

from sympy import Expr, Integer, Matrix, Rational, Symbol, expand, factor

Point: TypeAlias = tuple[Expr, Expr, Expr]
PolynomialMap: TypeAlias = tuple[Expr, Expr, Expr]

X: Final = Symbol("x")
Y: Final = Symbol("y")
Z: Final = Symbol("z")
VARIABLES: Final = (X, Y, Z)

TARGET: Final[Point] = (Rational(-1, 4), Integer(0), Integer(0))
DISPLAYED_POINTS: Final[tuple[Point, Point, Point]] = (
    (Integer(0), Integer(0), Rational(-1, 4)),
    (Integer(1), Rational(-3, 2), Rational(13, 2)),
    (Integer(-1), Rational(3, 2), Rational(13, 2)),
)


@dataclass(frozen=True)
class ExactCertificate:
    """The exact outputs produced by a symbolic verification run."""

    determinant: Expr
    images: tuple[Point, ...]
    target: Point

    @property
    def determinant_verified(self) -> bool:
        """Whether the formal determinant is identically ``-2``."""

        return bool(expand(self.determinant + 2) == 0)

    @property
    def collision_verified(self) -> bool:
        """Whether every supplied point maps to the displayed target."""

        return all(image == self.target for image in self.images)

    @property
    def verified(self) -> bool:
        """Whether both parts of the screenshot certificate hold."""

        return self.determinant_verified and self.collision_verified


def announced_map(*, displayed_four: int = 4) -> PolynomialMap:
    """Return the displayed map, optionally perturbing its coefficient ``4``.

    The keyword is exposed solely so the test suite can demonstrate that a
    one-character transcription error is detected.
    """

    u = 1 + X * Y
    tail = displayed_four + 3 * X * Y
    return (
        u**3 * Z + Y**2 * u * tail,
        Y + 3 * X * u**2 * Z + 3 * X * Y**2 * tail,
        2 * X - 3 * X**2 * Y - X**3 * Z,
    )


def evaluate(polynomial_map: PolynomialMap, point: Point) -> Point:
    """Evaluate a polynomial map at an exact point."""

    substitution = dict(zip(VARIABLES, point, strict=True))
    return tuple(expand(coordinate.subs(substitution)) for coordinate in polynomial_map)


def exact_certificate(
    *,
    displayed_four: int = 4,
    points: tuple[Point, ...] = DISPLAYED_POINTS,
) -> ExactCertificate:
    """Compute the exact determinant and the images of the supplied points."""

    polynomial_map = announced_map(displayed_four=displayed_four)
    matrix = Matrix(polynomial_map)
    determinant = factor(matrix.jacobian(VARIABLES).det())
    images = tuple(evaluate(polynomial_map, point) for point in points)
    return ExactCertificate(determinant=determinant, images=images, target=TARGET)


def main() -> int:
    """Print the exact certificate and return failure if either check fails."""

    certificate = exact_certificate()
    print(f"det JF = {certificate.determinant}")
    for point, image in zip(DISPLAYED_POINTS, certificate.images, strict=True):
        print(f"F{point} = {image}")
    print(f"determinant verified: {certificate.determinant_verified}")
    print(f"collision verified: {certificate.collision_verified}")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
