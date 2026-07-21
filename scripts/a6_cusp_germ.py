"""Exact hostile local germ for the surviving one-dicritical ``A6`` cusp.

Orevkov's 2026 classification of smooth finite map germs branched over a
one-Puiseux-pair curve contains the parameters

``(k1, k2; l1, l2) = (2, 1; 2, 0)``.

They give ``(d1, d2, N, n) = (2, 5, 5, 3)`` and hence exactly the local
degree-five, ramification-index-three block forced by the one-dicritical
``A6`` analysis.  This module derives and checks an explicit representative.
It is a hostile local model, not a Keller map or a global rank-six cover.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Final, TypeAlias

from sympy import (
    Expr,
    Integer,
    Poly,
    Rational,
    Symbol,
    discriminant,
    expand,
    factor,
)

PolynomialMap: TypeAlias = tuple[Expr, Expr]

X: Final = Symbol("x")
Y: Final = Symbol("y")
LAMBDA: Final = Symbol("lambda")
TARGET_U: Final = Symbol("U")
TARGET_V: Final = Symbol("V")


@dataclass(frozen=True)
class OrevkovCaseBParameters:
    """The integer parameters in case (b) of Orevkov's Theorem 2."""

    k1: int
    k2: int
    l1: int
    l2: int

    @property
    def m1(self) -> int:
        """Return ``k2*l2 + 1``."""

        return self.k2 * self.l2 + 1

    @property
    def m2(self) -> int:
        """Return ``k1*l1 + 1``."""

        return self.k1 * self.l1 + 1

    @property
    def d1(self) -> int:
        """Return the first exponent of the target branch equation."""

        return self.k1 * self.m1

    @property
    def d2(self) -> int:
        """Return the second exponent of the target branch equation."""

        return self.k2 * self.m2

    @property
    def degree(self) -> int:
        """Return the finite covering degree ``N``."""

        return self.m1 * self.m2

    @property
    def ramification_order(self) -> int:
        """Return the generic ramification order ``n``."""

        return self.l1 + self.l2 + 1

    @property
    def satisfies_hypotheses(self) -> bool:
        """Check every numerical hypothesis in Theorem 2(b)."""

        return (
            self.k1 >= 1
            and self.k2 >= 1
            and self.l1 >= 0
            and self.l2 >= 0
            and self.k1 + self.l2 > 1
            and self.k2 + self.l1 > 1
            and gcd(self.k1, self.k2) == 1
            and gcd(self.m1, self.m2) == 1
        )


A6_CUSP_PARAMETERS: Final = OrevkovCaseBParameters(2, 1, 2, 0)


def a6_cusp_germ(
    *,
    mixed_coefficient: Expr = Rational(-2, 3),
) -> PolynomialMap:
    """Return Orevkov's normalized degree-five germ.

    ``mixed_coefficient`` exposes one coefficient solely for the hostile test
    that demonstrates that the square-Jacobian certificate is nonvacuous.
    """

    primitive = (
        X**5 / Integer(5)
        + mixed_coefficient * Y * X**3
        + Y**2 * X
    )
    return (expand(Rational(15, 8) * primitive), Y)


@dataclass(frozen=True)
class A6CuspGermCertificate:
    """Exact symbolic outputs used to validate the hostile local model."""

    parameters: OrevkovCaseBParameters
    polynomial_map: PolynomialMap
    jacobian: Expr
    critical_image_u: Expr
    critical_image_v: Expr
    cusp_pullback: Expr
    residual_factor: Expr
    residual_parameter_polynomial: Expr
    residual_image_relation_remainder: Expr
    fiber_polynomial: Expr
    fiber_discriminant: Expr

    @property
    def parameters_verified(self) -> bool:
        """Check the exact case-(b) specialization relevant to ``A6``."""

        return (
            self.parameters.satisfies_hypotheses
            and (
                self.parameters.d1,
                self.parameters.d2,
                self.parameters.degree,
                self.parameters.ramification_order,
            )
            == (2, 5, 5, 3)
        )

    @property
    def finite_degree_verified(self) -> bool:
        """Check that ``v=y`` and ``u`` has nonzero constant leading term."""

        u, v = self.polynomial_map
        u_as_x_polynomial = Poly(u, X)
        return bool(
            v == Y
            and u_as_x_polynomial.degree() == 5
            and u_as_x_polynomial.LC() == Rational(3, 8)
        )

    @property
    def ramification_verified(self) -> bool:
        """Check the smooth critical curve and its multiplicity two."""

        return bool(
            expand(self.jacobian - Rational(15, 8) * (X**2 - Y) ** 2)
            == 0
        )

    @property
    def critical_component_verified(self) -> bool:
        """Check the normalization parametrization of the critical component."""

        return bool(
            self.critical_image_u == X**5
            and self.critical_image_v == X**2
            and expand(self.critical_image_u**2 - self.critical_image_v**5)
            == 0
        )

    @property
    def full_pullback_verified(self) -> bool:
        """Check the critical multiplicity and both residual components."""

        expected = (X**2 - Y) ** 3 * self.residual_factor / Integer(64)
        return bool(expand(self.cusp_pullback - expected) == 0)

    @property
    def residual_components_verified(self) -> bool:
        """Check exact evidence for two further bijective cusp components.

        Write either residual component as ``y=lambda*x^2``.  The quadratic
        for ``lambda`` is square-free and has no zero root.  If ``c(lambda)``
        is the coefficient in ``u=c(lambda)*x^5``, the checked congruence
        ``c(lambda)^2=lambda^5`` lets one choose a square root ``r`` of
        ``lambda`` and parametrize the target cusp by ``s=+/-r*x``.  Thus the
        restriction on each component is linear in its normalization
        parameter and hence locally bijective.
        """

        parameter_polynomial = Poly(
            self.residual_parameter_polynomial,
            LAMBDA,
        )
        return bool(
            parameter_polynomial.degree() == 2
            and parameter_polynomial.discriminant() == -1215
            and parameter_polynomial.eval(Integer(0)) == 9
            and self.residual_image_relation_remainder == 0
        )

    @property
    def discriminant_verified(self) -> bool:
        """Check the square ``A6`` discriminant with its exact scalar."""

        expected = (
            Integer(1_036_800_000) * (TARGET_U**2 - TARGET_V**5) ** 2
        )
        return bool(expand(self.fiber_discriminant - expected) == 0)

    @property
    def verified(self) -> bool:
        """Whether all exact local-germ checks pass."""

        return (
            self.parameters_verified
            and self.finite_degree_verified
            and self.ramification_verified
            and self.critical_component_verified
            and self.full_pullback_verified
            and self.residual_components_verified
            and self.discriminant_verified
        )


def exact_a6_cusp_germ_certificate(
    *,
    mixed_coefficient: Expr = Rational(-2, 3),
) -> A6CuspGermCertificate:
    """Build the exact certificate, optionally with a hostile perturbation."""

    u, v = a6_cusp_germ(mixed_coefficient=mixed_coefficient)
    jacobian = factor(u.diff(X) * v.diff(Y) - u.diff(Y) * v.diff(X))
    cusp_pullback = factor(u**2 - v**5)
    residual_factor = 9 * X**4 - 33 * X**2 * Y + 64 * Y**2
    residual_parameter_polynomial = (
        64 * LAMBDA**2 - 33 * LAMBDA + 9
    )
    residual_image_coefficient = expand(u.subs({X: Integer(1), Y: LAMBDA}))
    residual_image_relation = Poly(
        residual_image_coefficient**2 - LAMBDA**5,
        LAMBDA,
    )
    residual_image_relation_remainder = expand(
        residual_image_relation.rem(
            Poly(residual_parameter_polynomial, LAMBDA)
        ).as_expr()
    )
    fiber_polynomial = (
        3 * X**5
        - 10 * TARGET_V * X**3
        + 15 * TARGET_V**2 * X
        - 8 * TARGET_U
    )
    fiber_discriminant = factor(discriminant(fiber_polynomial, X))
    return A6CuspGermCertificate(
        parameters=A6_CUSP_PARAMETERS,
        polynomial_map=(u, v),
        jacobian=jacobian,
        critical_image_u=expand(u.subs(Y, X**2)),
        critical_image_v=expand(v.subs(Y, X**2)),
        cusp_pullback=cusp_pullback,
        residual_factor=residual_factor,
        residual_parameter_polynomial=residual_parameter_polynomial,
        residual_image_relation_remainder=residual_image_relation_remainder,
        fiber_polynomial=fiber_polynomial,
        fiber_discriminant=fiber_discriminant,
    )


def main() -> int:
    """Print the local certificate and fail if any exact check breaks."""

    certificate = exact_a6_cusp_germ_certificate()
    u, v = certificate.polynomial_map
    print(f"(u, v) = ({u}, {v})")
    print(f"det dF = {certificate.jacobian}")
    print(f"u^2-v^5 = {certificate.cusp_pullback}")
    print(f"disc_X(fiber) = {certificate.fiber_discriminant}")
    print(
        "(d1, d2, N, n) =",
        (
            certificate.parameters.d1,
            certificate.parameters.d2,
            certificate.parameters.degree,
            certificate.parameters.ramification_order,
        ),
    )
    print(f"certificate verified: {certificate.verified}")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
