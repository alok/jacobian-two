"""Exact hostile fixtures for the variable-leading quadratic proof target.

The reduction documented in ``docs/variable-leading-quadratic.md`` is not yet
Lean-certified.  This module is an independent SymPy regression checker for
the algebraic identities most likely to conceal a sign, truncation, or
denominator error.  Exact symbolic checks are evidence for those identities;
they are not a proof of the general theorem.
"""

from __future__ import annotations

from argparse import ArgumentParser
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final, TypeAlias

from sympy import Expr, Integer, Rational, Symbol, binomial, cancel, expand

PlaneMap: TypeAlias = tuple[Expr, Expr]

X: Final = Symbol("x")
Y: Final = Symbol("y")
KAPPA: Final = Symbol("kappa", nonzero=True)

NONMATE_QUADRATIC: Final[Expr] = X + X**2 * Y**2
RATIONAL_MATE_QUADRATIC: Final[Expr] = X + X**3 * Y**2


def jacobian(polynomial_map: PlaneMap) -> Expr:
    """Return ``P_x*Q_y - P_y*Q_x`` in the repository's sign convention."""

    first, second = polynomial_map
    return cancel(
        first.diff(X) * second.diff(Y)
        - first.diff(Y) * second.diff(X)
    )


@dataclass(frozen=True)
class RationalMateCertificate:
    """The rational Jacobian-one mate and its unavoidable denominator."""

    second_coordinate: Expr
    rational_mate: Expr
    determinant: Expr
    numerator: Expr
    denominator: Expr

    @property
    def verified(self) -> bool:
        """Whether the determinant is one and the mate is visibly nonpolynomial."""

        return (
            self.determinant == Integer(1)
            and self.denominator != Integer(1)
            and bool(self.denominator.free_symbols)
        )


def rational_mate_certificate() -> RationalMateCertificate:
    """Check ``J(-x*y/Q,Q)=1`` for ``Q=x+x^3*y^2`` exactly."""

    second = RATIONAL_MATE_QUADRATIC
    mate = cancel(-X * Y / second)
    numerator, denominator = mate.as_numer_denom()
    return RationalMateCertificate(
        second_coordinate=second,
        rational_mate=mate,
        determinant=jacobian((mate, second)),
        numerator=numerator,
        denominator=denominator,
    )


@dataclass(frozen=True)
class GradientIdealCertificate:
    """A Bezout certificate that the hostile quadratic has no critical point."""

    polynomial: Expr
    derivative_x: Expr
    derivative_y: Expr
    multiplier_x: Expr
    multiplier_y: Expr
    bezout_combination: Expr

    @property
    def verified(self) -> bool:
        """Whether the two partial derivatives generate the unit ideal."""

        return bool(self.bezout_combination == Integer(1))


def gradient_unit_ideal_certificate() -> GradientIdealCertificate:
    """Give an explicit unit-ideal witness for ``Q=x+x^2*y^2``.

    The identity

    ``(1-2*x*y^2)*Q_x + 2*y^3*Q_y = 1``

    proves over every characteristic-zero field that the two partial
    derivatives have no common zero, without relying on numerical sampling.
    """

    polynomial = NONMATE_QUADRATIC
    derivative_x = polynomial.diff(X)
    derivative_y = polynomial.diff(Y)
    multiplier_x = 1 - 2 * X * Y**2
    multiplier_y = 2 * Y**3
    bezout_combination = expand(
        multiplier_x * derivative_x + multiplier_y * derivative_y
    )
    return GradientIdealCertificate(
        polynomial=polynomial,
        derivative_x=derivative_x,
        derivative_y=derivative_y,
        multiplier_x=multiplier_x,
        multiplier_y=multiplier_y,
        bezout_combination=bezout_combination,
    )


@dataclass(frozen=True)
class LeadingEquationCertificate:
    """A nonconstant-square solution of the top coefficient equation."""

    degree: int
    leading_quadratic_coefficient: Expr
    leading_mate_coefficient: Expr
    residual: Expr

    @property
    def verified(self) -> bool:
        """Whether the equation vanishes although ``a=x^2`` is nonconstant."""

        return (
            self.degree >= 1
            and self.residual == Integer(0)
            and self.leading_quadratic_coefficient.diff(X) != Integer(0)
        )


def leading_equation_certificate(degree: int) -> LeadingEquationCertificate:
    """Check ``2*a*p_n' - n*a'*p_n=0`` for ``a=x^2,p_n=x^n``."""

    if degree < 1:
        msg = "the hostile leading degree must be positive"
        raise ValueError(msg)
    a = X**2
    p_n = X**degree
    residual = expand(
        2 * a * p_n.diff(X)
        - Integer(degree) * a.diff(X) * p_n
    )
    return LeadingEquationCertificate(
        degree=degree,
        leading_quadratic_coefficient=a,
        leading_mate_coefficient=p_n,
        residual=residual,
    )


@dataclass(frozen=True)
class CentralBinomialCertificate:
    """The exact product accumulated by the coefficient descent."""

    depth: int
    recurrence_product: Expr
    central_binomial_form: Expr
    residual: Expr

    @property
    def verified(self) -> bool:
        """Whether both nonzero rational expressions agree exactly."""

        return (
            self.depth >= 0
            and self.residual == Integer(0)
            and self.recurrence_product != Integer(0)
        )


def central_binomial_certificate(depth: int) -> CentralBinomialCertificate:
    """Check ``prod (2r+1)/(2r+2) = binomial(2m,m)/4^m``."""

    if depth < 0:
        msg = "the central-binomial depth must be nonnegative"
        raise ValueError(msg)
    recurrence_product: Expr = Integer(1)
    for index in range(depth):
        recurrence_product *= Rational(2 * index + 1, 2 * index + 2)
    central_binomial_form = cancel(
        binomial(2 * depth, depth) / Integer(4) ** depth
    )
    return CentralBinomialCertificate(
        depth=depth,
        recurrence_product=recurrence_product,
        central_binomial_form=central_binomial_form,
        residual=cancel(recurrence_product - central_binomial_form),
    )


@dataclass(frozen=True)
class TruncatedRecurrenceCertificate:
    """A finite forced tail and the nonzero obstruction at its top degree."""

    depth: int
    kappa: Expr
    coefficients: tuple[Expr, ...]
    truncated_mate: Expr
    determinant: Expr
    recurrence_residuals: tuple[Expr, ...]
    top_residual: Expr
    expected_top_residual: Expr

    @property
    def verified(self) -> bool:
        """Whether all forced equations vanish but the terminal one does not."""

        return (
            self.depth >= 0
            and self.kappa != Integer(0)
            and all(residual == Integer(0)
                    for residual in self.recurrence_residuals)
            and all(coefficient != Integer(0)
                    for coefficient in self.coefficients)
            and self.top_residual == self.expected_top_residual
            and self.top_residual != Integer(0)
        )


def truncated_recurrence_certificate(
    depth: int,
    *,
    kappa: Expr = KAPPA,
) -> TruncatedRecurrenceCertificate:
    """Build the forced odd tail for ``Q=x+x^2*y^2`` through ``depth``.

    ``depth=0`` retains only ``p_1=-kappa``. For each positive ``r``, the
    coefficient of ``y^(2r)`` forces

    ``p_(2r+1) = -(2r)/(2r+1) * x * p_(2r-1)``.

    Every imposed coefficient equation vanishes, but the first equation above
    the truncation remains a nonzero monomial. This supplies an exact finite
    witness to the infinite-tail obstruction at any requested depth.
    """

    if depth < 0:
        msg = "the recurrence depth must be nonnegative"
        raise ValueError(msg)
    if kappa == Integer(0):
        msg = "the Keller constant must be nonzero"
        raise ValueError(msg)

    coefficient_list: list[Expr] = [-kappa]
    recurrence_residual_list: list[Expr] = []
    for index in range(1, depth + 1):
        previous = coefficient_list[index - 1]
        current = expand(
            -Rational(2 * index, 2 * index + 1) * X * previous
        )
        coefficient_list.append(current)
        ell = 2 * index
        recurrence_residual_list.append(
            expand(
                2 * X**2 * previous.diff(X)
                - 2 * X * (ell - 1) * previous
                - (ell + 1) * current
            )
        )

    truncated_mate: Expr = Integer(0)
    for index, coefficient in enumerate(coefficient_list):
        truncated_mate += coefficient * Y ** (2 * index + 1)
    truncated_mate = expand(truncated_mate)

    determinant = expand(jacobian((truncated_mate, NONMATE_QUADRATIC)))
    top_residual = expand(determinant - kappa)
    top_coefficient = coefficient_list[-1]
    expected_top_residual = expand(
        -2 * (depth + 1) * X * top_coefficient * Y ** (2 * depth + 2)
    )
    return TruncatedRecurrenceCertificate(
        depth=depth,
        kappa=kappa,
        coefficients=tuple(coefficient_list),
        truncated_mate=truncated_mate,
        determinant=determinant,
        recurrence_residuals=tuple(recurrence_residual_list),
        top_residual=top_residual,
        expected_top_residual=expected_top_residual,
    )


def _argument_parser() -> ArgumentParser:
    """Build the small command-line interface for configurable depth."""

    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--depth",
        type=int,
        default=8,
        help="number of forced recurrence steps to check (default: 8)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Print all exact hostile certificates at the requested depth."""

    arguments = _argument_parser().parse_args(argv)
    depth = int(arguments.depth)
    if depth < 0:
        _argument_parser().error("--depth must be nonnegative")

    rational = rational_mate_certificate()
    gradient = gradient_unit_ideal_certificate()
    leading = leading_equation_certificate(2 * depth + 1)
    central = central_binomial_certificate(depth)
    recurrence = truncated_recurrence_certificate(depth)

    print(f"rational mate determinant: {rational.determinant}")
    print(f"rational mate denominator: {rational.denominator}")
    print(f"gradient Bezout combination: {gradient.bezout_combination}")
    print(f"top-equation residual: {leading.residual}")
    print(f"central-binomial residual: {central.residual}")
    print(f"forced recurrence depth: {recurrence.depth}")
    print(f"forced recurrence top residual: {recurrence.top_residual}")

    verified = (
        rational.verified
        and gradient.verified
        and leading.verified
        and central.verified
        and recurrence.verified
    )
    print(f"all hostile identities verified: {verified}")
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
