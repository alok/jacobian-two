"""Certify the global geometry behind the three true split pair fibers.

For an unordered pair of distinct normalization sources with sum ``s`` and
product ``r``, equality of the first coordinate is

``(2*s+k)*r - s*(s^2+k*s+1) = 0``.

This total pair-incidence surface is a smooth irreducible hypersurface.  Its
only fibers in which the coefficient of ``r`` and the constant term vanish
simultaneously occur at ``(k,s)=(0,0),(2,-1),(-2,1)``.  Those fibers split as
the exact vertical and graph components used by the delta-ten ledger.

Thus the vertical split branches are not separate irreducible components of
the total source-pair space.  This statement does not rule out new vertical
components after target-equality/contact equations are imposed: coefficient
rank and augmented compatibility still have to be audited profile by
profile.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

from sympy import Expr, Rational, diff, expand, gcd, groebner

from scripts.a6_delta_ten_generic import (
    KAPPA,
    PAIR_DENOMINATOR,
    PAIR_INCIDENCE,
    PAIR_QUADRATIC,
    R,
    S,
)


@dataclass(frozen=True, slots=True)
class SplitPairSurfaceCertificate:
    """Exact irreducibility, smoothness, and split-fiber identities."""

    linear_content_gcd: Expr
    singular_ideal_basis: tuple[Expr, ...]
    exceptional_base_identity: Expr
    split_fiber_identities: tuple[Expr, ...]
    component_intersection_residuals: tuple[tuple[Expr, Expr], ...]
    total_surface_irreducible: bool
    total_surface_smooth: bool
    flat_over_kappa: bool

    @property
    def verified(self) -> bool:
        """Whether the total pair surface and all three fibers agree."""

        return bool(
            self.linear_content_gcd == 1
            and self.singular_ideal_basis == (1,)
            and self.exceptional_base_identity == 0
            and self.split_fiber_identities == (0, 0, 0)
            and self.component_intersection_residuals == ((0, 0), (0, 0), (0, 0))
            and self.total_surface_irreducible
            and self.total_surface_smooth
            and self.flat_over_kappa
        )


@cache
def exact_split_pair_surface_certificate() -> SplitPairSurfaceCertificate:
    """Build the exact global unordered-pair surface certificate.

    Irreducibility follows because the defining polynomial is primitive and
    linear in ``r``: its two coefficients have gcd one in ``QQ[k,s]``.
    Smoothness is certified by the unit Jacobian ideal.  The coordinate ring
    is then a domain, hence torsion-free over the PID ``QQ[k]`` and therefore
    flat over the ``k`` line.
    """

    constant_coefficient = expand(-S * PAIR_QUADRATIC)
    singular_basis = groebner(
        (
            PAIR_INCIDENCE,
            diff(PAIR_INCIDENCE, KAPPA),
            diff(PAIR_INCIDENCE, S),
            diff(PAIR_INCIDENCE, R),
        ),
        R,
        S,
        KAPPA,
        order="lex",
    )
    expected_exceptional = expand(S * (S**2 - 1))
    exceptional_after_denominator = expand(constant_coefficient.subs(KAPPA, -2 * S))

    split_factors = (
        S * (2 * R - S**2 - 1),
        (S + 1) * (2 * R - S * (S + 1)),
        (S - 1) * (2 * R - S * (S - 1)),
    )
    split_values = (0, 2, -2)
    intersections = (
        {KAPPA: 0, S: 0, R: Rational(1, 2)},
        {KAPPA: 2, S: -1, R: 0},
        {KAPPA: -2, S: 1, R: 0},
    )
    vertical_factors = (S, S + 1, S - 1)
    graph_factors = (
        2 * R - S**2 - 1,
        2 * R - S * (S + 1),
        2 * R - S * (S - 1),
    )
    return SplitPairSurfaceCertificate(
        linear_content_gcd=expand(gcd(PAIR_DENOMINATOR, constant_coefficient)),
        singular_ideal_basis=tuple(
            polynomial.as_expr() for polynomial in singular_basis.polys
        ),
        exceptional_base_identity=expand(
            exceptional_after_denominator - expected_exceptional
        ),
        split_fiber_identities=tuple(
            expand(PAIR_INCIDENCE.subs(KAPPA, kappa) - factor)
            for kappa, factor in zip(split_values, split_factors, strict=True)
        ),
        component_intersection_residuals=tuple(
            (
                expand(vertical.subs(point)),
                expand(graph.subs(point)),
            )
            for point, vertical, graph in zip(
                intersections,
                vertical_factors,
                graph_factors,
                strict=True,
            )
        ),
        total_surface_irreducible=True,
        total_surface_smooth=True,
        flat_over_kappa=True,
    )


def main() -> int:
    """Print the exact total pair-surface summary."""

    certificate = exact_split_pair_surface_certificate()
    print("split pair surface certificate:", certificate.verified)
    print("singular ideal basis:", certificate.singular_ideal_basis)
    print("split fiber identities:", certificate.split_fiber_identities)
    print(
        "global geometry:",
        "irreducible=",
        certificate.total_surface_irreducible,
        "smooth=",
        certificate.total_surface_smooth,
        "flat-over-k=",
        certificate.flat_over_kappa,
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
