"""Record the exact closure of the split ``C2^2`` overlap-plus-``W`` locus.

The independent Sage checker
``tools/check_a6_delta_ten_split_c22_overlap_closure.sage`` constructs a
rational arc in the ordinary ordered two-contact Cramer incidence.  Its first
unordered pair approaches the vertical/graph intersection at ``k=0`` and its
second pair approaches a clean ``W`` pair.  The exact coefficient limit is
recorded here and compared with the complete affine-line solution of the true
split rank-three system.

This proves algebraic component containment.  It does not, by itself, prove
that affine-complement topology is constant at the overlap boundary.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import Expr, Symbol, cancel, diff

from scripts.a6_delta_ten_generic import ALPHA, BETA, DELTA, GAMMA
from scripts.a6_delta_ten_split_contact_rank import V, simple_double_contact_specs


ARC_PARAMETER: Final = Symbol("lambda_c22_overlap")
FREE_DELTA: Final = Symbol("delta_c22_overlap")

LIMIT_DENOMINATOR: Final = V**5 + 5 * V**3 - 8 * ARC_PARAMETER + 10 * V
LIMIT_COEFFICIENTS: Final = (
    cancel(
        (3 * V**7 + 14 * V**5 + 21 * V**3 - 8 * ARC_PARAMETER + 14 * V)
        / (4 * LIMIT_DENOMINATOR)
    ),
    cancel(
        (
            -3 * ARC_PARAMETER * V**8
            + V**9
            - 13 * ARC_PARAMETER * V**6
            - V**7
            - 9 * ARC_PARAMETER * V**4
            - 21 * V**5
            + 19 * ARC_PARAMETER * V**2
            - 35 * V**3
            - 2 * ARC_PARAMETER
        )
        / (4 * V * LIMIT_DENOMINATOR)
    ),
    cancel(
        (3 * V**7 + 15 * V**5 + 26 * V**3 - 16 * ARC_PARAMETER + 24 * V)
        / (2 * LIMIT_DENOMINATOR)
    ),
    cancel(
        (
            -3 * ARC_PARAMETER * V**8
            + 3 * V**9
            - 13 * ARC_PARAMETER * V**6
            - 3 * V**7
            + 5 * ARC_PARAMETER * V**4
            - 63 * V**5
            + 65 * ARC_PARAMETER * V**2
            - 105 * V**3
            - 6 * ARC_PARAMETER
        )
        / (8 * V * LIMIT_DENOMINATOR)
    ),
)

SPLIT_SOLUTION_DENOMINATOR: Final = V**6 + 5 * V**4 - 5 * V**2 + 15
COMPLETE_SPLIT_SOLUTION: Final = (
    cancel(
        -16 * FREE_DELTA * V / SPLIT_SOLUTION_DENOMINATOR
        + (3 * V**8 + 14 * V**6 - 70 * V**2 + 21) / (4 * SPLIT_SOLUTION_DENOMINATOR)
    ),
    cancel(
        FREE_DELTA
        * (2 * V**6 + 10 * V**4 + 10 * V**2 + 10)
        / SPLIT_SOLUTION_DENOMINATOR
        + (-(V**9) + V**7 + 21 * V**5 + 35 * V**3) / (2 * SPLIT_SOLUTION_DENOMINATOR)
    ),
    cancel(
        -32 * FREE_DELTA * V / SPLIT_SOLUTION_DENOMINATOR
        + (3 * V**8 + 15 * V**6 + 5 * V**4 - 75 * V**2 + 36)
        / (2 * SPLIT_SOLUTION_DENOMINATOR)
    ),
    FREE_DELTA,
)

LAMBDA_FOR_FREE_DELTA: Final = cancel(
    V**2
    * (
        -8 * FREE_DELTA * V**4
        - 40 * FREE_DELTA * V**2
        - 80 * FREE_DELTA
        + 3 * V**7
        - 3 * V**5
        - 63 * V**3
        - 105 * V
    )
    / (-64 * FREE_DELTA * V + 3 * V**8 + 13 * V**6 - 5 * V**4 - 65 * V**2 + 6)
)

# The independent Sage replay checks four total-pair rows, a generically
# invertible 4 x 4 Cramer matrix, regular epsilon limits, and equality with
# the formulas above.  Keep this as explicit external-check metadata rather
# than pretending that the Python formula comparison reconstructed the lift.
SAGE_ARC_LIFT_METADATA: Final = (4, 4, True, True)


@dataclass(frozen=True, slots=True)
class DoubleContactOverlapClosureCertificate:
    """Exact limit/dominance comparison for the overlap-plus-``W`` incidence."""

    split_equation_residuals: tuple[Expr, ...]
    dominance_residuals: tuple[Expr, ...]
    free_parameter_residual: Expr
    arc_parameter_derivative: Expr
    sage_arc_lift_metadata: tuple[int, int, bool, bool]
    topology_computed: bool

    @property
    def algebraically_contained(self) -> bool:
        """Whether the ordinary component dominates the full split incidence."""

        return bool(
            all(residual == 0 for residual in self.split_equation_residuals)
            and all(residual == 0 for residual in self.dominance_residuals)
            and self.free_parameter_residual == 0
            and self.arc_parameter_derivative != 0
            and self.sage_arc_lift_metadata == SAGE_ARC_LIFT_METADATA
            and self.sage_arc_lift_metadata[2:]
            == (
                True,
                True,
            )
        )

    @property
    def verified(self) -> bool:
        """Whether containment and the topology claim boundary both agree."""

        return self.algebraically_contained and not self.topology_computed


@cache
def exact_double_contact_overlap_closure_certificate() -> (
    DoubleContactOverlapClosureCertificate
):
    """Build the exact split-limit and dominant-parametrization certificate."""

    overlap = next(
        spec for spec in simple_double_contact_specs() if spec.allocation == "overlap+W"
    )
    coefficient_substitution = dict(
        zip(
            (ALPHA, BETA, GAMMA, DELTA),
            LIMIT_COEFFICIENTS,
            strict=True,
        )
    )
    dominant_limit = tuple(
        cancel(expression.subs(ARC_PARAMETER, LAMBDA_FOR_FREE_DELTA))
        for expression in LIMIT_COEFFICIENTS
    )
    return DoubleContactOverlapClosureCertificate(
        split_equation_residuals=tuple(
            cancel(equation.subs(coefficient_substitution))
            for equation in overlap.equations
        ),
        dominance_residuals=tuple(
            cancel(actual - expected)
            for actual, expected in zip(
                dominant_limit,
                COMPLETE_SPLIT_SOLUTION,
                strict=True,
            )
        ),
        free_parameter_residual=cancel(dominant_limit[3] - FREE_DELTA),
        arc_parameter_derivative=cancel(diff(LIMIT_COEFFICIENTS[3], ARC_PARAMETER)),
        sage_arc_lift_metadata=SAGE_ARC_LIFT_METADATA,
        topology_computed=False,
    )


def main() -> int:
    """Print the exact algebraic-containment conclusion and topology boundary."""

    certificate = exact_double_contact_overlap_closure_certificate()
    print(
        "C2^2 overlap+W algebraically contained:", certificate.algebraically_contained
    )
    print("topology computed:", certificate.topology_computed)
    print("proves JC(2): False")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
