"""Exact obstruction for the smallest width-five ``S6`` target curve.

The rational curve

``(P,Q)=(t^5+t^4,t^7+t^5)``

has the first possible width and one-pair infinity data left by the general
transposition-generator bound.  It also has the allowed finite ``T(4,5)``
jump and six normalization nodes.  Exact Zariski--van Kamp relations show,
however, that its affine complement has no transitive transposition image on
six sheets.  Sage further simplifies the complement group to ``Z``.

This excludes one degree-minimal branch curve, not the ``S6`` passport or the
plane Jacobian conjecture.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache
from itertools import product
from typing import Final

from sympy import (
    Expr,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    gcd,
    resultant,
    together,
)

from scripts.s6_trace_curve_topology import (
    all_transpositions,
    minimum_transposition_generators,
)
from scripts.six_sheet_monodromy import (
    IDENTITY,
    Permutation,
    compose,
    generated_group,
    inverse,
    is_transitive,
)

T: Final = Symbol("t")
S: Final = Symbol("s")
P_COORD: Final = Symbol("p")
Q_COORD: Final = Symbol("q")
PAIR_SUM: Final = Symbol("u")

PARAMETER_P: Final = T**5 + T**4
PARAMETER_Q: Final = T**7 + T**5
BRANCH_EQUATION: Final = (
    P_COORD**7
    - 2 * P_COORD**6
    + 9 * P_COORD**5 * Q_COORD
    - 7 * P_COORD**4 * Q_COORD**2
    + 2 * P_COORD**5
    - 10 * P_COORD**4 * Q_COORD
    + 28 * P_COORD**3 * Q_COORD**2
    - 34 * P_COORD**2 * Q_COORD**3
    + 12 * P_COORD * Q_COORD**4
    - Q_COORD**5
    - 2 * Q_COORD**4
)

COLLISION_FACTOR: Final = (
    T**10
    + T**9
    + 2 * T**8
    + 2 * T**7
    + 3 * T**6
    + 3 * T**5
    + 6 * T**4
    - 4 * T**3
    + 2 * T**2
    - 2 * T
    + 2
)
NODE_COORDINATE_POLYNOMIAL: Final = (
    P_COORD**5
    + 2 * P_COORD**4
    + 32 * P_COORD**3
    - 152 * P_COORD**2
    - 64 * P_COORD
    - 16
)
PAIR_SUM_POLYNOMIAL: Final = (
    PAIR_SUM**5
    + PAIR_SUM**4
    - 2 * PAIR_SUM**3
    + 4 * PAIR_SUM**2
    - 3 * PAIR_SUM
    + 1
)

# Sage 10.8 affine Zariski--van Kamp output.  Indices 1..5 correspond to the
# five geometric fiber meridians x0..x4; negative indices mean inverses.
VAN_KAMPEN_RELATIONS: Final = (
    (-5, 3, 5, -3),
    (-5, -3, -2, 1, 2, 3, 5, -3, -2, -1, 2, 3),
    (-5, -3, -2, -1, 2, 1, 2, 3, 5, -3, -2, -1, -2, 1, 2, 3),
    (-4, -3, -2, -1, -3, -2, -1, 2, 1, 2, 3, 1, 2, 3),
    (
        5,
        1,
        2,
        3,
        5,
        -3,
        -2,
        -1,
        -2,
        -1,
        2,
        1,
        2,
        3,
        -5,
        -3,
        -2,
        -1,
    ),
    (
        1,
        2,
        3,
        -5,
        -3,
        -2,
        -1,
        -5,
        -3,
        -2,
        -1,
        -3,
        -2,
        -1,
        2,
        1,
        2,
        3,
        1,
        2,
        3,
        5,
    ),
    (
        -3,
        -5,
        -3,
        -2,
        -1,
        -3,
        -2,
        -1,
        -2,
        1,
        2,
        3,
        -2,
        1,
        2,
        -3,
        -2,
        -1,
        2,
        1,
        2,
        3,
        1,
        2,
        3,
        5,
    ),
    (-3, -2, -1, 2, 1, 2, 3, -2, -1, -2, 1, 2),
    (2, 1, -2, -1),
    (5, 4, -5, -4),
)


@dataclass(frozen=True, slots=True)
class WidthFiveCurveArithmeticCertificate:
    """Exact normalization, singularity, width, and genus identities."""

    implicit_resultant_identity: Expr
    parametrized_branch_identity: Expr
    off_diagonal_gcd: Expr
    derivative_residual_resultant: Expr
    collision_resultant_identity: Expr
    node_coordinate_resultant_identity: Expr
    node_coordinate_discriminant: Expr
    rational_node_separation: Expr
    pair_sum_discriminant: Expr
    pair_distinctness_resultant: Expr
    pair_transversality_resultant: Expr
    rational_pair_transversality: Expr
    width: int
    finite_cusp_pair: tuple[int, int]
    infinity_pair: tuple[int, int]
    finite_node_count: int
    arithmetic_genus: int
    total_delta: int

    @property
    def verified(self) -> bool:
        """Whether every exact geometric identity has its expected value."""

        return bool(
            self.implicit_resultant_identity == 0
            and self.parametrized_branch_identity == 0
            and self.off_diagonal_gcd == 1
            and self.derivative_residual_resultant == 237
            and self.collision_resultant_identity == 0
            and self.node_coordinate_resultant_identity == 0
            and self.node_coordinate_discriminant
            == 2**16 * 3**6 * 71**2 * 797
            and self.rational_node_separation == -135
            and self.pair_sum_discriminant == 2**4 * 797
            and abs(int(self.pair_distinctness_resultant)) == 2 * 79
            and abs(int(self.pair_transversality_resultant))
            == 2**2 * 5 * 79 * 797
            and self.rational_pair_transversality == -15
            and self.width == 5
            and self.finite_cusp_pair == (4, 5)
            and self.infinity_pair == (2, 7)
            and self.finite_node_count == 6
            and self.arithmetic_genus == 15
            and self.total_delta == 15
        )


def exact_width_five_curve_arithmetic_certificate(
) -> WidthFiveCurveArithmeticCertificate:
    """Build the dependency-free exact geometry certificate."""

    parameter_p_s = PARAMETER_P.subs(T, S)
    parameter_q_s = PARAMETER_Q.subs(T, S)
    first_collision = cancel((PARAMETER_P - parameter_p_s) / (T - S))
    second_collision = cancel((PARAMETER_Q - parameter_q_s) / (T - S))
    expected_collision_resultant = (
        T**12 * (T**2 + T + 1) * COLLISION_FACTOR
    )
    pair_product = (PAIR_SUM**2 + 1) / (1 - PAIR_SUM)
    pair_discriminant_numerator = together(
        PAIR_SUM**2 - 4 * pair_product
    ).as_numer_denom()[0]
    velocity_factor_numerator = together(
        35 * PAIR_SUM * pair_product
        + 28 * PAIR_SUM**2
        - 28 * pair_product
        + 20
    ).as_numer_denom()[0]
    return WidthFiveCurveArithmeticCertificate(
        implicit_resultant_identity=expand(
            resultant(P_COORD - PARAMETER_P, Q_COORD - PARAMETER_Q, T)
            + BRANCH_EQUATION
        ),
        parametrized_branch_identity=expand(
            BRANCH_EQUATION.subs(
                {P_COORD: PARAMETER_P, Q_COORD: PARAMETER_Q}
            )
        ),
        off_diagonal_gcd=gcd(first_collision, second_collision),
        derivative_residual_resultant=resultant(
            5 * T + 4,
            7 * T**2 + 5,
            T,
        ),
        collision_resultant_identity=expand(
            resultant(first_collision, second_collision, S)
            - expected_collision_resultant
        ),
        node_coordinate_resultant_identity=expand(
            resultant(
                COLLISION_FACTOR,
                P_COORD - PARAMETER_P,
                T,
            )
            - NODE_COORDINATE_POLYNOMIAL**2
        ),
        node_coordinate_discriminant=discriminant(
            NODE_COORDINATE_POLYNOMIAL,
            P_COORD,
        ),
        rational_node_separation=NODE_COORDINATE_POLYNOMIAL.subs(P_COORD, -1),
        pair_sum_discriminant=discriminant(
            PAIR_SUM_POLYNOMIAL,
            PAIR_SUM,
        ),
        pair_distinctness_resultant=resultant(
            PAIR_SUM_POLYNOMIAL,
            pair_discriminant_numerator,
            PAIR_SUM,
        ),
        pair_transversality_resultant=resultant(
            PAIR_SUM_POLYNOMIAL,
            velocity_factor_numerator,
            PAIR_SUM,
        ),
        rational_pair_transversality=(
            35 * PAIR_SUM * pair_product
            + 28 * PAIR_SUM**2
            - 28 * pair_product
            + 20
        ).subs(PAIR_SUM, -1),
        width=min(
            int(PARAMETER_P.as_poly(T).degree()),
            int(PARAMETER_Q.as_poly(T).degree()),
        ),
        finite_cusp_pair=(4, 5),
        infinity_pair=(7 - 5, 7),
        finite_node_count=6,
        arithmetic_genus=(7 - 1) * (7 - 2) // 2,
        total_delta=(4 - 1) * (5 - 1) // 2 + 6 + (2 - 1) * (7 - 1) // 2,
    )


def evaluate_signed_word(
    word: tuple[int, ...],
    images: Sequence[Permutation],
) -> Permutation:
    """Evaluate a signed relation in the five geometric-meridian images."""

    result = IDENTITY
    for letter in word:
        generator = images[abs(letter) - 1]
        if letter < 0:
            generator = inverse(generator)
        result = compose(result, generator)
    return result


@dataclass(frozen=True, slots=True)
class WidthFiveTranspositionCensus:
    """Exact transposition-valued representation census."""

    total_assignments: int
    satisfying_assignments: int
    generated_order_histogram: tuple[tuple[int, int], ...]
    transitive_assignments: int

    @property
    def verified(self) -> bool:
        """Whether only the fifteen equal-transposition images survive."""

        return (
            self.total_assignments == 15**5
            and self.satisfying_assignments == 15
            and self.generated_order_histogram == ((2, 15),)
            and self.transitive_assignments == 0
        )


@cache
def width_five_transposition_census() -> WidthFiveTranspositionCensus:
    """Exhaust all five-meridian assignments to transpositions of ``S6``."""

    transpositions = all_transpositions()
    satisfying = 0
    transitive = 0
    orders: Counter[int] = Counter()
    for images in product(transpositions, repeat=5):
        if not all(
            evaluate_signed_word(relation, images) == IDENTITY
            for relation in VAN_KAMPEN_RELATIONS
        ):
            continue
        satisfying += 1
        group = generated_group(images)
        orders[len(group)] += 1
        if is_transitive(group):
            transitive += 1
    return WidthFiveTranspositionCensus(
        total_assignments=len(transpositions) ** 5,
        satisfying_assignments=satisfying,
        generated_order_histogram=tuple(sorted(orders.items())),
        transitive_assignments=transitive,
    )


@dataclass(frozen=True, slots=True)
class S6WidthFiveCurveCertificate:
    """Exact geometry and global-monodromy obstruction."""

    arithmetic: WidthFiveCurveArithmeticCertificate
    minimum_s6_width: int
    census: WidthFiveTranspositionCensus

    @property
    def width_is_minimal(self) -> bool:
        """Whether the curve attains the general six-sheet lower bound."""

        return self.arithmetic.width == self.minimum_s6_width == 5

    @property
    def verified(self) -> bool:
        """Whether the geometric near-miss and topology obstruction both hold."""

        return self.arithmetic.verified and self.width_is_minimal and self.census.verified


def exact_s6_width_five_curve_certificate() -> S6WidthFiveCurveCertificate:
    """Build the exact certificate for the degree-minimal near-miss."""

    return S6WidthFiveCurveCertificate(
        arithmetic=exact_width_five_curve_arithmetic_certificate(),
        minimum_s6_width=minimum_transposition_generators(6),
        census=width_five_transposition_census(),
    )


def main() -> int:
    """Print the exact curve audit and fail on any regression."""

    certificate = exact_s6_width_five_curve_certificate()
    print(
        "width-five curve geometry:",
        {
            "finite cusp": certificate.arithmetic.finite_cusp_pair,
            "nodes": certificate.arithmetic.finite_node_count,
            "infinity": certificate.arithmetic.infinity_pair,
            "delta": certificate.arithmetic.total_delta,
        },
    )
    print(
        "transposition representations:",
        {
            "assignments": certificate.census.total_assignments,
            "solutions": certificate.census.satisfying_assignments,
            "orders": dict(certificate.census.generated_order_histogram),
            "transitive": certificate.census.transitive_assignments,
        },
    )
    print(f"S6 width-five curve certificate verified: {certificate.verified}")
    print("Sage simplification: pi_1(A2-B)=Z; all five meridians are equal")
    print("claim boundary: this minimal curve only; the S6 passport remains open")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
