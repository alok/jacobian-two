"""Topological obstruction for the explicit hostile ``S6`` trace curve.

The symmetric trace matrix in ``trace_hostile_matrices.py`` uses the rational
curve ``(P,Q)=(t^4-6t^2,t^5-5t^3)``.  Its quadratic trace data are correct,
but its affine complement cannot support the required transposition
monodromy on six sheets.

A Sage Zariski--van Kamp computation gives the four-meridian presentation
stored below.  This dependency-free module replays every homomorphism that
sends the four meridians to transpositions.  More generally, four
transpositions cannot generate a transitive group on six letters: their edge
graph has too few edges to be connected.  Thus this curve is a stopping model
for trace lattices, not a possible global ``S6`` branch curve.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import Final

from sympy import Expr, Symbol, expand

from scripts.six_sheet_monodromy import (
    IDENTITY,
    Permutation,
    compose,
    generated_group,
    inverse,
    is_transitive,
)

P: Final = Symbol("P")
Q: Final = Symbol("Q")
PARAMETER: Final = Symbol("t")

PARAMETER_P: Final = PARAMETER**4 - 6 * PARAMETER**2
PARAMETER_Q: Final = PARAMETER**5 - 5 * PARAMETER**3
BRANCH_EQUATION: Final = (
    P**5
    + 10 * P**4
    + 25 * P**3
    + 10 * P**2 * Q**2
    + 90 * P * Q**2
    - Q**4
    + 216 * Q**2
)

# Sage generators 1..4 correspond to x0..x3.  Negative indices mean inverse.
VAN_KAMPEN_RELATIONS: Final = (
    (-4, 2, 4, -2),
    (-4, 1, 4, -1),
    (3, 4, 3, -4, -3, -4),
    (-3, 1, 3, 1, -3, -1),
    (3, 2, 3, -2, -3, -2),
    (-3, -2, 3, 1, -3, 2, 3, -1),
)


def all_transpositions() -> tuple[Permutation, ...]:
    """Return the fifteen transpositions in the natural six-point action."""

    result: list[Permutation] = []
    for first in range(6):
        for second in range(first + 1, 6):
            image = list(IDENTITY)
            image[first], image[second] = image[second], image[first]
            result.append(tuple(image))
    return tuple(result)


def evaluate_signed_word(
    word: tuple[int, ...],
    images: tuple[Permutation, Permutation, Permutation, Permutation],
) -> Permutation:
    """Evaluate one signed relation with ``gh = g after h``."""

    result = IDENTITY
    for letter in word:
        generator = images[abs(letter) - 1]
        if letter < 0:
            generator = inverse(generator)
        result = compose(result, generator)
    return result


def assignment_satisfies_relations(
    images: tuple[Permutation, Permutation, Permutation, Permutation],
) -> bool:
    """Check every Zariski--van Kamp relator exactly."""

    return all(
        evaluate_signed_word(relation, images) == IDENTITY
        for relation in VAN_KAMPEN_RELATIONS
    )


@dataclass(frozen=True, slots=True)
class TranspositionAssignmentCensus:
    """Exact census of transposition-valued representations."""

    total_assignments: int
    relation_satisfying_assignments: int
    group_order_histogram: tuple[tuple[int, int], ...]
    transitive_assignments: int

    @property
    def verified(self) -> bool:
        """Whether the exact expected census is recovered."""

        return (
            self.total_assignments == 15**4
            and self.relation_satisfying_assignments == 735
            and self.group_order_histogram == ((2, 15), (120, 720))
            and self.transitive_assignments == 0
        )


def transposition_assignment_census() -> TranspositionAssignmentCensus:
    """Enumerate all four-meridian maps into the transpositions of ``S6``."""

    transpositions = all_transpositions()
    satisfying = 0
    transitive = 0
    orders: Counter[int] = Counter()
    for raw_images in product(transpositions, repeat=4):
        images = (
            raw_images[0],
            raw_images[1],
            raw_images[2],
            raw_images[3],
        )
        if not assignment_satisfies_relations(images):
            continue
        satisfying += 1
        group = generated_group(images)
        orders[len(group)] += 1
        if is_transitive(group):
            transitive += 1
    return TranspositionAssignmentCensus(
        total_assignments=len(transpositions) ** 4,
        relation_satisfying_assignments=satisfying,
        group_order_histogram=tuple(sorted(orders.items())),
        transitive_assignments=transitive,
    )


def minimum_transposition_generators(sheet_count: int) -> int:
    """Return the connected-edge lower bound for transitive generation.

    Regard each transposition as one edge on the sheet set.  The generated
    group preserves every connected component of this graph.  Transitivity
    therefore requires a connected graph, which needs at least ``n-1`` edges.
    """

    if sheet_count < 1:
        msg = "the sheet count must be positive"
        raise ValueError(msg)
    return sheet_count - 1


@dataclass(frozen=True, slots=True)
class S6TraceCurveTopologyCertificate:
    """Exact implicit-curve and complement-monodromy obstruction."""

    parametrized_branch_equation: Expr
    projection_width: int
    minimum_required_width: int
    census: TranspositionAssignmentCensus

    @property
    def width_obstructs_s6(self) -> bool:
        """Whether meridian generation is too narrow for six-sheet transitivity."""

        return self.projection_width < self.minimum_required_width

    @property
    def verified(self) -> bool:
        """Whether the curve identity, width theorem, and exact census agree."""

        return bool(
            self.parametrized_branch_equation == 0
            and self.projection_width == 4
            and self.minimum_required_width == 5
            and self.width_obstructs_s6
            and self.census.verified
        )


def exact_s6_trace_curve_topology_certificate(
) -> S6TraceCurveTopologyCertificate:
    """Build the exact global-topology obstruction for the hostile curve."""

    return S6TraceCurveTopologyCertificate(
        parametrized_branch_equation=expand(
            BRANCH_EQUATION.subs({P: PARAMETER_P, Q: PARAMETER_Q})
        ),
        projection_width=4,
        minimum_required_width=minimum_transposition_generators(6),
        census=transposition_assignment_census(),
    )


def main() -> int:
    """Print the exact topology census and fail on any regression."""

    certificate = exact_s6_trace_curve_topology_certificate()
    print(
        "trace-curve width:",
        certificate.projection_width,
        "required for transitive transpositions:",
        certificate.minimum_required_width,
    )
    print(
        "transposition representations:",
        {
            "solutions": certificate.census.relation_satisfying_assignments,
            "orders": dict(certificate.census.group_order_histogram),
            "transitive": certificate.census.transitive_assignments,
        },
    )
    print(f"S6 trace-curve topology certificate verified: {certificate.verified}")
    print("claim boundary: this explicit trace curve only; the S6 passport remains")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
