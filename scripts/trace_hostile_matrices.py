"""Exact hostile matrices for the surviving one-dicritical trace lattices.

The finite-normalization argument determines symmetric matrices, determinants,
cokernels, and pointwise coranks.  This module constructs explicit polynomial
matrices satisfying all of those quadratic constraints for representative
``A6`` and ``S6`` branch curves.  They deliberately carry no rank-six
commutative multiplication law and therefore are not finite Keller algebras.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from sympy import Expr, Matrix, Symbol, diag, expand, factor, sqrt, zeros

P: Final = Symbol("P")
Q: Final = Symbol("Q")
PARAMETER: Final = Symbol("t")


@dataclass(frozen=True)
class A6TraceHostileCertificate:
    """Exact outputs for a square-discriminant rank-six hostile matrix."""

    branch_equation: Expr
    parametrized_branch_equation: Expr
    presentation_determinant: Expr
    presentation_kernel: tuple[Expr, ...]
    generic_presentation_rank: int
    cusp_presentation_rank: int
    node_presentation_ranks: tuple[int, int]
    phi_determinant: Expr
    middle_determinant: Expr
    trace_determinant: Expr
    generic_phi_rank: int
    cusp_phi_rank: int
    generic_trace_rank: int
    cusp_trace_rank: int
    trace_unit_norm: Expr

    @property
    def verified(self) -> bool:
        """Whether every exact ``A6`` lattice constraint is realized."""

        return bool(
            self.parametrized_branch_equation == 0
            and expand(self.presentation_determinant - self.branch_equation)
            == 0
            and all(entry == 0 for entry in self.presentation_kernel)
            and self.generic_presentation_rank == 2
            and self.cusp_presentation_rank == 1
            and self.node_presentation_ranks == (1, 1)
            and expand(self.phi_determinant - self.branch_equation) == 0
            and self.middle_determinant == -1
            and expand(self.trace_determinant + self.branch_equation**2) == 0
            and self.generic_phi_rank == 5
            and self.cusp_phi_rank == 4
            and self.generic_trace_rank == 4
            and self.cusp_trace_rank == 2
            and self.trace_unit_norm == 6
        )


def exact_a6_trace_hostile_certificate(
    *,
    bottom_right_shift: int = 0,
) -> A6TraceHostileCertificate:
    """Build the hostile ``A6`` trace matrix and all of its exact checks.

    ``bottom_right_shift`` is exposed solely for an adversarial regression
    test; a nonzero value perturbs the symmetric presentation.
    """

    branch_equation = (
        -P**5
        + 5 * P**3 * Q
        - 5 * P * Q**2
        + Q**3
        + Q**2
    )
    parameter_p = PARAMETER**2 + PARAMETER**3
    parameter_q = PARAMETER**5
    presentation = Matrix(
        [
            [-P**2, P**2 - Q, P - Q],
            [P**2 - Q, -Q, -P],
            [P - Q, -P, P - 1 + bottom_right_shift],
        ]
    )
    normalization_generators = Matrix([1, PARAMETER, PARAMETER**2])
    parametrized_presentation = presentation.subs(
        {P: parameter_p, Q: parameter_q}
    )

    identity = Matrix.eye(3)
    zero = zeros(3)
    phi = diag(presentation, identity)
    middle = zero.row_join(identity).col_join(identity.row_join(zero))
    trace = phi.T * middle * phi
    trace_unit = Matrix([0, 0, 1, 0, -3, -3])

    generic_point = {P: 12, Q: 32}  # parameter t=2
    cusp_point = {P: 0, Q: 0}
    node_points = (
        {P: (-1 + sqrt(5)) / 2, Q: 1},
        {P: (-1 - sqrt(5)) / 2, Q: 1},
    )

    return A6TraceHostileCertificate(
        branch_equation=branch_equation,
        parametrized_branch_equation=factor(
            branch_equation.subs({P: parameter_p, Q: parameter_q})
        ),
        presentation_determinant=factor(presentation.det()),
        presentation_kernel=tuple(
            factor(entry)
            for entry in parametrized_presentation * normalization_generators
        ),
        generic_presentation_rank=int(presentation.subs(generic_point).rank()),
        cusp_presentation_rank=int(presentation.subs(cusp_point).rank()),
        node_presentation_ranks=(
            int(presentation.subs(node_points[0]).rank()),
            int(presentation.subs(node_points[1]).rank()),
        ),
        phi_determinant=factor(phi.det()),
        middle_determinant=factor(middle.det()),
        trace_determinant=factor(trace.det()),
        generic_phi_rank=int(phi.subs(generic_point).rank()),
        cusp_phi_rank=int(phi.subs(cusp_point).rank()),
        generic_trace_rank=int(trace.subs(generic_point).rank()),
        cusp_trace_rank=int(trace.subs(cusp_point).rank()),
        trace_unit_norm=factor((trace_unit.T * trace * trace_unit)[0]),
    )


@dataclass(frozen=True)
class S6TraceHostileCertificate:
    """Exact outputs for a full three-cusp ``S6`` trace-form model."""

    branch_equation: Expr
    parametrized_branch_equation: Expr
    derivative_gcd: Expr
    presentation_determinant: Expr
    presentation_kernel: tuple[Expr, ...]
    generic_rank: int
    cusp_ranks: tuple[int, int, int]
    node_ranks: tuple[int, int, int]
    trace_determinant: Expr
    trace_unit_norm: Expr

    @property
    def verified(self) -> bool:
        """Whether every exact ``S6`` trace-form constraint is realized."""

        return bool(
            self.parametrized_branch_equation == 0
            and self.derivative_gcd == PARAMETER * (PARAMETER**2 - 3)
            and expand(self.presentation_determinant - self.branch_equation)
            == 0
            and all(entry == 0 for entry in self.presentation_kernel)
            and self.generic_rank == 5
            and self.cusp_ranks == (4, 4, 4)
            and self.node_ranks == (4, 4, 4)
            and expand(self.trace_determinant - self.branch_equation) == 0
            and self.trace_unit_norm == 6
        )


def exact_s6_trace_hostile_certificate() -> S6TraceHostileCertificate:
    """Build the three-``T(2,3)`` hostile symmetric presentation."""

    parameter_p = PARAMETER**4 - 6 * PARAMETER**2
    parameter_q = PARAMETER**5 - 5 * PARAMETER**3
    branch_equation = (
        P**5
        + 10 * P**4
        + 25 * P**3
        + 10 * P**2 * Q**2
        + 90 * P * Q**2
        - Q**4
        + 216 * Q**2
    )
    presentation = Matrix(
        [
            [0, 6 * Q, -5 * P, -Q, P],
            [6 * Q, -5 * P, -Q, P, 0],
            [-5 * P, -Q, P - 30, 0, 6],
            [-Q, P, 0, 1, 0],
            [P, 0, 6, 0, -1],
        ]
    )
    normalization_generators = Matrix(
        [1, PARAMETER, PARAMETER**2, PARAMETER**3, PARAMETER**4]
    )
    parametrized_presentation = presentation.subs(
        {P: parameter_p, Q: parameter_q}
    )
    trace = diag(presentation, 1)
    trace_unit = Matrix([0, 0, 0, 0, 0, sqrt(6)])

    cusp_points = (
        {P: 0, Q: 0},
        {P: -9, Q: -6 * sqrt(3)},
        {P: -9, Q: 6 * sqrt(3)},
    )
    node_points = (
        {P: -5, Q: 0},
        {P: -4, Q: 2 * sqrt(2)},
        {P: -4, Q: -2 * sqrt(2)},
    )
    derivative_gcd = factor(
        parameter_p.diff(PARAMETER).as_poly(PARAMETER).gcd(
            parameter_q.diff(PARAMETER).as_poly(PARAMETER)
        ).monic().as_expr()
    )

    return S6TraceHostileCertificate(
        branch_equation=branch_equation,
        parametrized_branch_equation=factor(
            branch_equation.subs({P: parameter_p, Q: parameter_q})
        ),
        derivative_gcd=derivative_gcd,
        presentation_determinant=factor(presentation.det()),
        presentation_kernel=tuple(
            factor(entry)
            for entry in parametrized_presentation * normalization_generators
        ),
        generic_rank=int(trace.subs({P: -5, Q: -4}).rank()),
        cusp_ranks=(
            int(trace.subs(cusp_points[0]).rank()),
            int(trace.subs(cusp_points[1]).rank()),
            int(trace.subs(cusp_points[2]).rank()),
        ),
        node_ranks=(
            int(trace.subs(node_points[0]).rank()),
            int(trace.subs(node_points[1]).rank()),
            int(trace.subs(node_points[2]).rank()),
        ),
        trace_determinant=factor(trace.det()),
        trace_unit_norm=factor((trace_unit.T * trace * trace_unit)[0]),
    )


def main() -> int:
    """Print both hostile certificates and fail if either breaks."""

    a6 = exact_a6_trace_hostile_certificate()
    s6 = exact_s6_trace_hostile_certificate()
    print(
        "A6 hostile trace:",
        {
            "detPhi": a6.phi_determinant,
            "detH": a6.middle_determinant,
            "detT": a6.trace_determinant,
            "generic coranks (Phi,T)": (
                6 - a6.generic_phi_rank,
                6 - a6.generic_trace_rank,
            ),
            "cusp coranks (Phi,T)": (
                6 - a6.cusp_phi_rank,
                6 - a6.cusp_trace_rank,
            ),
            "unit norm": a6.trace_unit_norm,
        },
    )
    print(
        "S6 hostile trace:",
        {
            "detT": s6.trace_determinant,
            "generic corank": 6 - s6.generic_rank,
            "cusp coranks": tuple(6 - rank for rank in s6.cusp_ranks),
            "node coranks": tuple(6 - rank for rank in s6.node_ranks),
            "unit norm": s6.trace_unit_norm,
        },
    )
    print(f"A6 certificate verified: {a6.verified}")
    print(f"S6 certificate verified: {s6.verified}")
    return 0 if a6.verified and s6.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
