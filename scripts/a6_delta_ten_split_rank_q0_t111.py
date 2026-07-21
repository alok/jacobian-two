"""Close the split rank strata for ``Q0+4N`` and ``T111^2+4N``.

The true pair charts at ``k=0,+/-2`` are reducible, but the ordinary
quadruple and two-ordinary-triple incidence systems are most efficiently
audited before splitting the pair polynomial.

For ``Q0+4N``, reduction of ``Q`` modulo the full fiber quartic ``P-h`` gives
three affine-linear equations in ``(a,b,c,d)``.  The four maximal minors have
common factor ``h^2`` and, after dividing by it, generate the unit ideal for
each split value.  Thus the coefficient rank is exactly three on every valid
nonzero fiber.  The rank-two fixture at ``h=0`` lies on the singular fiber
discriminant and is not an ordinary quadruple.

For ``T111^2+4N``, the global omitted-root determinant specializes cleanly.
At ``k=0`` its residual rank factor ``-(u+v)`` is contained in the
same-``P``-fiber factor, hence is invalid for two distinct triple targets.  At
``k=+/-2`` the only remaining rank factor is the residual divisor already
proved incompatible after Rabinowitsch localization by
``a6_delta_ten_residual_rank``.  A compatible rank-three point at ``k=2`` is
kept as a hostile fixture: it lies on the cusp-fiber boundary and shows why
that localization cannot be dropped.

Consequently both split profile incidences have dimension at most two.  This
is an incidence-rank result only; it does not compute complement topology or
settle any other split allocation, generic degree six, or ``JC(2)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import Expr, Rational, Symbol, cancel, discriminant, expand, groebner

from scripts.a6_delta_ten_double_triple import (
    DOUBLE_TRIPLE_EQUATIONS,
    EXPECTED_INCIDENCE_DETERMINANT,
    FIRST_CUSP_FIBER_FACTOR,
    FIRST_OMITTED_ROOT,
    INCIDENCE_AUGMENTED_MATRIX as TRIPLE_AUGMENTED_MATRIX,
    INCIDENCE_MATRIX as TRIPLE_MATRIX,
    RESIDUAL_RANK_FACTOR,
    SAME_FIBER_FACTOR,
    SECOND_CUSP_FIBER_FACTOR,
    SECOND_OMITTED_ROOT,
)
from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    DELTA,
    GAMMA,
    KAPPA,
    T,
)
from scripts.a6_delta_ten_quadruple import (
    COEFFICIENT_MATRIX as QUADRUPLE_MATRIX,
    COEFFICIENT_MINORS as QUADRUPLE_MINORS,
    EXPECTED_COEFFICIENT_MINORS,
    FIBER_POLYNOMIAL,
    FIBER_VALUE,
)
from scripts.a6_delta_ten_residual_rank import (
    exact_double_triple_residual_certificate,
)
from scripts.a6_delta_ten_split_codim_two import (
    exact_split_plus_minus_transport_certificate,
)

LOCALIZATION_VARIABLE: Final = Symbol("z_split_q0")
SPLIT_VALUES: Final = (0, 2, -2)


def _localized_quadruple_minor_basis(kappa: int) -> tuple[Expr, ...]:
    """Return the rank-at-most-two ideal localized at ``h != 0``.

    Correctness: the coefficient matrix has rank at most two exactly when all
    four maximal minors vanish.  Each minor is divisible by ``h^2``; dividing
    removes the invalid rank collapse forced by the singular zero fiber.
    Rabinowitsch localization then computes the ideal on ``D(h)`` exactly.
    """

    normalized = tuple(
        expand(cancel(minor / FIBER_VALUE**2).subs(KAPPA, kappa))
        for minor in QUADRUPLE_MINORS
    )
    basis = groebner(
        (
            *normalized,
            1 - LOCALIZATION_VARIABLE * FIBER_VALUE,
        ),
        LOCALIZATION_VARIABLE,
        FIBER_VALUE,
        order="lex",
    )
    return tuple(polynomial.as_expr() for polynomial in basis.polys)


@dataclass(frozen=True, slots=True)
class SplitQuadrupleRankCertificate:
    """Exact rank certificate for the three split ordinary-quadruple charts."""

    global_minor_identities: tuple[Expr, ...]
    localized_rankdrop_bases: tuple[tuple[Expr, ...], ...]
    zero_fiber_ranks: tuple[int, ...]
    fiber_discriminants: tuple[Expr, ...]
    zero_fiber_discriminant_values: tuple[Expr, ...]
    maximum_valid_incidence_dimension: int

    @property
    def verified(self) -> bool:
        """Whether rank three holds on every valid split quadruple fiber."""

        return bool(
            self.global_minor_identities == (0, 0, 0, 0)
            and self.localized_rankdrop_bases == ((1,), (1,), (1,))
            and self.zero_fiber_ranks == (2, 2, 2)
            and self.zero_fiber_discriminant_values == (0, 0, 0)
            and self.maximum_valid_incidence_dimension == 2
        )


@cache
def exact_split_quadruple_rank_certificate() -> SplitQuadrupleRankCertificate:
    """Prove the split ``Q0+4N`` coefficient rank is always three.

    At fixed ``k``, the valid fiber-value base has dimension one and a rank
    three affine system has a one-dimensional coefficient fiber.  Hence the
    incidence dimension is ``1 + (4 - 3) = 2``.
    """

    discriminants = tuple(
        expand(discriminant(FIBER_POLYNOMIAL.subs(KAPPA, kappa), T))
        for kappa in SPLIT_VALUES
    )
    return SplitQuadrupleRankCertificate(
        global_minor_identities=tuple(
            expand(actual - expected)
            for actual, expected in zip(
                QUADRUPLE_MINORS,
                EXPECTED_COEFFICIENT_MINORS,
                strict=True,
            )
        ),
        localized_rankdrop_bases=tuple(
            _localized_quadruple_minor_basis(kappa) for kappa in SPLIT_VALUES
        ),
        zero_fiber_ranks=tuple(
            QUADRUPLE_MATRIX.subs({KAPPA: kappa, FIBER_VALUE: 0}).rank()
            for kappa in SPLIT_VALUES
        ),
        fiber_discriminants=discriminants,
        zero_fiber_discriminant_values=tuple(
            expression.subs(FIBER_VALUE, 0) for expression in discriminants
        ),
        maximum_valid_incidence_dimension=2,
    )


K2_HOSTILE_BASE: Final = {
    KAPPA: 2,
    FIRST_OMITTED_ROOT: -1,
    SECOND_OMITTED_ROOT: Rational(1, 3),
}
K2_HOSTILE_MEMBER: Final = {
    **K2_HOSTILE_BASE,
    ALPHA: Rational(16, 27),
    BETA: -Rational(646, 297),
    GAMMA: -Rational(373, 99),
    DELTA: 0,
}


@dataclass(frozen=True, slots=True)
class SplitDoubleTripleRankCertificate:
    """Exact rank and compatibility certificate for split two-triple charts."""

    determinant_specializations: tuple[Expr, ...]
    residual_specializations: tuple[Expr, ...]
    same_fiber_specializations: tuple[Expr, ...]
    k0_residual_same_fiber_identity: Expr
    global_valid_residual_basis: tuple[Expr, ...]
    k0_hostile_ranks: tuple[int, int]
    k0_hostile_same_fiber_value: Expr
    k2_hostile_equation_residuals: tuple[Expr, ...]
    k2_hostile_residual_value: Expr
    k2_hostile_same_fiber_value: Expr
    k2_hostile_cusp_values: tuple[Expr, Expr]
    k2_hostile_ranks: tuple[int, int]
    plus_minus_transport_verified: bool
    maximum_valid_incidence_dimension: int

    @property
    def verified(self) -> bool:
        """Whether no compatible split rank-drop incidence survives."""

        u = FIRST_OMITTED_ROOT
        v = SECOND_OMITTED_ROOT
        return bool(
            self.residual_specializations
            == (-u - v, 2 * u * v - u - v, -2 * u * v - u - v)
            and all(
                expand(actual - expected) == 0
                for actual, expected in zip(
                    self.same_fiber_specializations,
                    (
                        (u + v) * (u**2 + v**2 + 1),
                        (u + v + 1) * (u**2 + u + v**2 + v),
                        (u + v - 1) * (u**2 - u + v**2 - v),
                    ),
                    strict=True,
                )
            )
            and self.k0_residual_same_fiber_identity == 0
            and self.global_valid_residual_basis == (1,)
            and self.k0_hostile_ranks == (3, 3)
            and self.k0_hostile_same_fiber_value == 0
            and self.k2_hostile_equation_residuals == (0, 0, 0, 0)
            and self.k2_hostile_residual_value == 0
            and self.k2_hostile_same_fiber_value != 0
            and self.k2_hostile_cusp_values == (0, Rational(16, 9))
            and self.k2_hostile_ranks == (3, 3)
            and self.plus_minus_transport_verified
            and self.maximum_valid_incidence_dimension == 2
        )


@cache
def exact_split_double_triple_rank_certificate() -> SplitDoubleTripleRankCertificate:
    """Prove the split ``T111^2+4N`` rank strata add no component.

    On the full-rank valid base the two omitted roots give a two-dimensional
    base and the coefficient fiber is a point.  The global residual
    saturation excludes every compatible rank-three point; lower rank can
    occur only on already removed coincident, cusp, or same-fiber factors.
    """

    u = FIRST_OMITTED_ROOT
    v = SECOND_OMITTED_ROOT
    residual_certificate = exact_double_triple_residual_certificate()
    transport = exact_split_plus_minus_transport_certificate()
    determinant_specializations = tuple(
        expand(EXPECTED_INCIDENCE_DETERMINANT.subs(KAPPA, kappa))
        for kappa in SPLIT_VALUES
    )
    residual_specializations = tuple(
        expand(RESIDUAL_RANK_FACTOR.subs(KAPPA, kappa)) for kappa in SPLIT_VALUES
    )
    same_fiber_specializations = tuple(
        expand(SAME_FIBER_FACTOR.subs(KAPPA, kappa)) for kappa in SPLIT_VALUES
    )
    k2_matrix = TRIPLE_MATRIX.subs(K2_HOSTILE_BASE)
    k2_augmented = TRIPLE_AUGMENTED_MATRIX.subs(K2_HOSTILE_BASE)
    return SplitDoubleTripleRankCertificate(
        determinant_specializations=determinant_specializations,
        residual_specializations=residual_specializations,
        same_fiber_specializations=same_fiber_specializations,
        k0_residual_same_fiber_identity=expand(
            same_fiber_specializations[0]
            + residual_specializations[0] * (u**2 + v**2 + 1)
        ),
        global_valid_residual_basis=residual_certificate.localized_groebner_basis,
        k0_hostile_ranks=(
            residual_certificate.hostile_coefficient_rank,
            residual_certificate.hostile_augmented_rank,
        ),
        k0_hostile_same_fiber_value=(residual_certificate.hostile_same_fiber_value),
        k2_hostile_equation_residuals=tuple(
            expand(equation.subs(K2_HOSTILE_MEMBER))
            for equation in DOUBLE_TRIPLE_EQUATIONS
        ),
        k2_hostile_residual_value=RESIDUAL_RANK_FACTOR.subs(K2_HOSTILE_BASE),
        k2_hostile_same_fiber_value=SAME_FIBER_FACTOR.subs(K2_HOSTILE_BASE),
        k2_hostile_cusp_values=(
            FIRST_CUSP_FIBER_FACTOR.subs(K2_HOSTILE_BASE),
            SECOND_CUSP_FIBER_FACTOR.subs(K2_HOSTILE_BASE),
        ),
        k2_hostile_ranks=(k2_matrix.rank(), k2_augmented.rank()),
        plus_minus_transport_verified=transport.verified,
        maximum_valid_incidence_dimension=2,
    )


@dataclass(frozen=True, slots=True)
class SplitQ0DoubleTripleRankCertificate:
    """Combined exact closure of the two global-fiber split rank systems."""

    quadruple: SplitQuadrupleRankCertificate
    double_triple: SplitDoubleTripleRankCertificate

    @property
    def verified(self) -> bool:
        """Whether both profile bundles satisfy the dimension-two bound."""

        return self.quadruple.verified and self.double_triple.verified


@cache
def exact_split_q0_double_triple_rank_certificate() -> (
    SplitQ0DoubleTripleRankCertificate
):
    """Build the combined split ``Q0`` and two-triple certificate."""

    return SplitQ0DoubleTripleRankCertificate(
        quadruple=exact_split_quadruple_rank_certificate(),
        double_triple=exact_split_double_triple_rank_certificate(),
    )


def main() -> int:
    """Print the exact split global-fiber rank summary."""

    certificate = exact_split_q0_double_triple_rank_certificate()
    print("split Q0/two-triple rank certificate:", certificate.verified)
    print(
        "Q0 localized rank-drop bases:",
        certificate.quadruple.localized_rankdrop_bases,
    )
    print(
        "two-triple valid residual basis:",
        certificate.double_triple.global_valid_residual_basis,
    )
    print(
        "maximum valid incidence dimensions:",
        certificate.quadruple.maximum_valid_incidence_dimension,
        certificate.double_triple.maximum_valid_incidence_dimension,
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
