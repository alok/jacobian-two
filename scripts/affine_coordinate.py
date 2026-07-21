"""Exact fixtures for the arbitrary-degree affine-coordinate theorem.

The Lean development proves the theorem.  This independent SymPy module
checks high-degree instances of both normal forms, their displayed inverses,
and a hostile nonconstant-slope perturbation.
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


def evaluate(polynomial_map: PlaneMap, point: PlaneMap) -> PlaneMap:
    """Evaluate a plane polynomial map by simultaneous substitution."""

    substitution = {X: point[0], Y: point[1]}
    return tuple(expand(coordinate.subs(substitution, simultaneous=True))
                 for coordinate in polynomial_map)


def jacobian(polynomial_map: PlaneMap) -> Expr:
    """Return the exact formal Jacobian determinant in ``(x,y)`` order."""

    return cancel(Matrix(polynomial_map).jacobian((X, Y)).det())


@dataclass(frozen=True)
class AffineCoordinateCertificate:
    """Exact identities for the nonzero-slope normal form."""

    polynomial_map: PlaneMap
    inverse: PlaneMap
    determinant: Expr
    expected_determinant: Expr
    inverse_after_map: PlaneMap
    map_after_inverse: PlaneMap

    @property
    def verified(self) -> bool:
        """Whether the determinant and both inverse compositions are exact."""

        return (
            expand(self.determinant - self.expected_determinant) == 0
            and self.inverse_after_map == (X, Y)
            and self.map_after_inverse == (U, V)
        )


def nonzero_slope_certificate() -> AffineCoordinateCertificate:
    """Check a degree-seven nonzero-slope normal form."""

    eps = Integer(3)
    alpha = Integer(5)
    beta = Integer(-2)
    f = X**5 - 2 * X**2 + 7
    generator = eps * Y + f
    outer = T**7 - 3 * T**4 + 2 * T**2 + 11 * T
    first = outer.subs(T, generator) + alpha * X + beta
    polynomial_map = (first, generator)

    inverse_x = expand((U - outer.subs(T, V) - beta) / alpha)
    inverse_y = expand((V - f.subs(X, inverse_x)) / eps.subs(X, inverse_x))
    inverse = (inverse_x, inverse_y)

    recovered_x = cancel((first - outer.subs(T, generator) - beta) / alpha)
    recovered_y = cancel((generator - f.subs(X, recovered_x)) / eps)
    inverse_after_map = (recovered_x, recovered_y)
    recovered_v = cancel(eps * inverse_y + f.subs(X, inverse_x))
    recovered_u = cancel(outer.subs(T, recovered_v) + alpha * inverse_x + beta)
    map_after_inverse = (recovered_u, recovered_v)
    return AffineCoordinateCertificate(
        polynomial_map=polynomial_map,
        inverse=inverse,
        determinant=jacobian(polynomial_map),
        expected_determinant=alpha * eps,
        inverse_after_map=inverse_after_map,
        map_after_inverse=map_after_inverse,
    )


def nonconstant_slope_determinant() -> Expr:
    """Return the Jacobian after the hostile replacement ``eps=3+x``.

    Keeping the same normal-form construction gives determinant ``5*(3+x)``.
    It is an exact witness that a nonconstant slope fails the Keller premise;
    no rational inverse simplification is needed for this negative fixture.
    """

    eps = X + 3
    alpha = Integer(5)
    f = X**5 - 2 * X**2 + 7
    generator = eps * Y + f
    outer = T**7 - 3 * T**4 + 2 * T**2 + 11 * T
    polynomial_map = (outer.subs(T, generator) + alpha * X - 2, generator)
    return jacobian(polynomial_map)


@dataclass(frozen=True)
class ZeroSlopeCertificate:
    """Exact identities for the complementary triangular normal form."""

    determinant: Expr
    expected_determinant: Expr
    inverse_after_map: PlaneMap
    map_after_inverse: PlaneMap

    @property
    def verified(self) -> bool:
        """Whether the determinant and both inverse compositions are exact."""

        return (
            expand(self.determinant - self.expected_determinant) == 0
            and self.inverse_after_map == (X, Y)
            and self.map_after_inverse == (U, V)
        )


def zero_slope_certificate() -> ZeroSlopeCertificate:
    """Check the arbitrary-``c(x)`` triangular chart exactly."""

    eta = Integer(7)
    delta = Integer(-4)
    gamma = Integer(9)
    c = X**8 - 6 * X**3 + X - 1
    polynomial_map = (expand(eta * Y + c), delta * X + gamma)

    inverse_x = expand((V - gamma) / delta)
    inverse_y = expand((U - c.subs(X, inverse_x)) / eta)
    inverse = (inverse_x, inverse_y)
    inverse_after_map = tuple(
        expand(coordinate.subs({U: polynomial_map[0], V: polynomial_map[1]},
                               simultaneous=True))
        for coordinate in inverse
    )
    map_after_inverse = tuple(
        expand(coordinate.subs({X: inverse_x, Y: inverse_y}, simultaneous=True))
        for coordinate in polynomial_map
    )
    return ZeroSlopeCertificate(
        determinant=jacobian(polynomial_map),
        expected_determinant=-(eta * delta),
        inverse_after_map=inverse_after_map,
        map_after_inverse=map_after_inverse,
    )


def main() -> int:
    """Print the two exact certificates and reject the hostile fixture."""

    nonzero = nonzero_slope_certificate()
    zero = zero_slope_certificate()
    perturbed_determinant = nonconstant_slope_determinant()
    perturbation_rejected = X in perturbed_determinant.free_symbols
    print(f"nonzero-slope determinant: {nonzero.determinant}")
    print(f"nonzero-slope normal form verified: {nonzero.verified}")
    print(f"zero-slope determinant: {zero.determinant}")
    print(f"zero-slope normal form verified: {zero.verified}")
    print(f"nonconstant-slope determinant: {perturbed_determinant}")
    print(f"hostile perturbation rejected: {perturbation_rejected}")
    return 0 if nonzero.verified and zero.verified and perturbation_rejected else 1


if __name__ == "__main__":
    raise SystemExit(main())
