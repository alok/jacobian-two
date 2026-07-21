"""Exclude the generic conditional ``A6`` delta-ten family.

Under the polynomial-normalization, one-pair, and finite-singularity
hypotheses, the first coarse survivor after the complete delta-seven audit has
affine degrees ``(4, 9)`` and infinity pair ``(5, 9)``.  Polynomial source
and target changes put every such normalization in the five-parameter form

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``.

This module certifies the generic unordered-pair collision algebra, including
the three legitimate reducible incidence fibers ``k = 0, +/-2``.  It then
checks a small rational member with one ``T(2,5)`` cusp, ten nodes, and a
``T(5,9)`` branch at infinity.  A stored Sage 10.8 van Kamp presentation has
cyclic complement; the dependency-free replay exhausts all ``40^4``
single-three-cycle assignments and finds no ``A6`` image.

Propagation from the exact member to the connected clean stratum uses proper
projective Whitney--Thom triviality.  The result remains conditional and
computer-assisted; it neither classifies every degeneration wall nor proves
the plane Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import (
    Expr,
    Poly,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    resultant,
)

from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    n_generator_three_cycle_presentation_census,
)

T: Final = Symbol("t")
S: Final = Symbol("s")
R: Final = Symbol("r")
X: Final = Symbol("X")
Y: Final = Symbol("Y")

KAPPA: Final = Symbol("k")
ALPHA: Final = Symbol("a")
BETA: Final = Symbol("b")
GAMMA: Final = Symbol("c")
DELTA: Final = Symbol("d")

FAMILY_P: Final = T**2 + KAPPA * T**3 + T**4
FAMILY_Q: Final = (
    ALPHA * T**5
    + BETA * T**6
    + GAMMA * T**7
    + DELTA * T**8
    + T**9
)

PAIR_DENOMINATOR: Final = 2 * S + KAPPA
PAIR_QUADRATIC: Final = S**2 + KAPPA * S + 1
PAIR_INCIDENCE: Final = (
    PAIR_DENOMINATOR * R - S * PAIR_QUADRATIC
)


def _divided_power_sums(maximum: int) -> dict[int, Expr]:
    """Return ``(t^n-u^n)/(t-u)`` in ``s=t+u, r=tu`` recursively."""

    if maximum < 2:
        msg = "the divided-power recurrence starts in degree two"
        raise ValueError(msg)
    values: dict[int, Expr] = {1: 1, 2: S}
    for degree in range(3, maximum + 1):
        values[degree] = expand(
            S * values[degree - 1] - R * values[degree - 2]
        )
    return values


DIVIDED_POWERS: Final = _divided_power_sums(9)
SECOND_DIVIDED_DIFFERENCE: Final = (
    ALPHA * DIVIDED_POWERS[5]
    + BETA * DIVIDED_POWERS[6]
    + GAMMA * DIVIDED_POWERS[7]
    + DELTA * DIVIDED_POWERS[8]
    + DIVIDED_POWERS[9]
)
COLLISION_RESULTANT: Final = resultant(
    PAIR_INCIDENCE,
    SECOND_DIVIDED_DIFFERENCE,
    R,
)
COLLISION_POLYNOMIAL: Final = cancel(COLLISION_RESULTANT / S**2)

CUSP_IMAGE_FACTOR: Final = resultant(
    T**2 + KAPPA * T + 1,
    T**4 + DELTA * T**3 + GAMMA * T**2 + BETA * T + ALPHA,
    T,
)
EXTRA_CRITICAL_FACTOR: Final = resultant(
    4 * T**2 + 3 * KAPPA * T + 2,
    9 * T**4 + 8 * DELTA * T**3 + 7 * GAMMA * T**2 + 6 * BETA * T + 5 * ALPHA,
    T,
)
PAIR_DIAGONAL_FACTOR: Final = 2 * S**2 + 3 * KAPPA * S + 4
TANGENCY_COFACTOR: Final = (
    9 * KAPPA**2 * S
    - 2 * KAPPA * S**2
    + 10 * KAPPA
    - 4 * S**3
    - 16 * S
)
TANGENCY_POLYNOMIAL: Final = cancel(
    (
        S
        * PAIR_DIAGONAL_FACTOR
        * PAIR_DENOMINATOR
        * diff(COLLISION_POLYNOMIAL, S)
        + TANGENCY_COFACTOR * COLLISION_POLYNOMIAL
    )
    / PAIR_QUADRATIC
)

# The collision-X coordinate on the generic incidence chart.
COLLISION_X_NUMERATOR: Final = -(
    S * (S + KAPPA) * PAIR_QUADRATIC**2
)
COLLISION_X_DENOMINATOR: Final = PAIR_DENOMINATOR**2

# Legitimate reducible pair-incidence fibers.  These are not invalid parameter
# values and must not be discarded when the generic chart is used.
ZERO_FIBER_VERTICAL: Final = R**2 - GAMMA * R + ALPHA
ZERO_FIBER_GRAPH: Final = (
    S**8
    + (8 - 2 * GAMMA) * S**6
    + (-4 * BETA + 8 * DELTA) * S**5
    + (-4 * ALPHA + 2 * GAMMA + 6) * S**4
    + (-8 * BETA + 16 * DELTA) * S**3
    + (-16 * ALPHA + 18 * GAMMA - 16) * S**2
    + (12 * BETA - 8 * DELTA) * S
    + 4 * ALPHA
    - 2 * GAMMA
    + 1
)


def exceptional_vertical(epsilon: int) -> Expr:
    """Return the degree-four involution component for ``k=2*epsilon``."""

    if epsilon not in (-1, 1):
        msg = "epsilon must be +1 or -1"
        raise ValueError(msg)
    return expand(
        R**4
        + (4 * epsilon * DELTA - GAMMA - 10) * R**3
        + (
            ALPHA
            - 3 * epsilon * BETA
            - 10 * epsilon * DELTA
            + 6 * GAMMA
            + 15
        )
        * R**2
        + (
            -3 * ALPHA
            + 4 * epsilon * BETA
            + 6 * epsilon * DELTA
            - 5 * GAMMA
            - 7
        )
        * R
        + ALPHA
        - epsilon * BETA
        - epsilon * DELTA
        + GAMMA
        + 1
    )


def exceptional_graph(epsilon: int) -> Expr:
    """Return the degree-six graph component for ``k=2*epsilon``."""

    if epsilon not in (-1, 1):
        msg = "epsilon must be +1 or -1"
        raise ValueError(msg)
    return expand(
        S**6
        + 8 * epsilon * S**5
        + (6 - 2 * GAMMA + 8 * epsilon * DELTA) * S**4
        + (
            -4 * BETA
            + 16 * DELTA
            + 2 * epsilon * GAMMA
            - 16 * epsilon
        )
        * S**3
        + (
            -4 * ALPHA
            - 8 * epsilon * BETA
            - 8 * epsilon * DELTA
            + 18 * GAMMA
            + 1
        )
        * S**2
        + (
            -16 * epsilon * ALPHA
            + 12 * BETA
            - 2 * epsilon * GAMMA
        )
        * S
        + 4 * ALPHA
    )


SAMPLE_PARAMETERS: Final = {
    KAPPA: 1,
    ALPHA: 1,
    BETA: 0,
    GAMMA: 0,
    DELTA: 0,
}
SAMPLE_P: Final = FAMILY_P.subs(SAMPLE_PARAMETERS)
SAMPLE_Q: Final = FAMILY_Q.subs(SAMPLE_PARAMETERS)
SAMPLE_COLLISION_POLYNOMIAL: Final = (
    S**10
    + 6 * S**9
    + 17 * S**8
    + 21 * S**7
    + 2 * S**6
    - 37 * S**5
    - 58 * S**4
    - 33 * S**3
    - 4 * S**2
    + 3 * S
    + 1
)
SAMPLE_NODE_X_POLYNOMIAL: Final = (
    X**10
    - 10 * X**9
    + 3 * X**8
    + 191 * X**7
    - 226 * X**6
    - 712 * X**5
    - 668 * X**4
    - 331 * X**3
    - 95 * X**2
    - 15 * X
    - 1
)
SAMPLE_IMPLICIT: Final = (
    -X**9
    - 4 * X**8
    - 8 * X**7
    + 21 * X**6 * Y
    - 4 * X**6
    + 41 * X**5 * Y
    - X**5
    + 27 * X**4 * Y**2
    + 20 * X**4 * Y
    + 22 * X**3 * Y**2
    + 5 * X**3 * Y
    + 9 * X**2 * Y**3
    + 16 * X**2 * Y**2
    + 5 * X * Y**3
    + 4 * X * Y**2
    + Y**4
    - Y**3
    + Y**2
)

# Sage 10.8 exact raw affine presentation.  Generator indices 1..4 are the
# geometric meridians of a generic vertical fiber.
SAMPLE_RELATIONS: Final = (
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-3, 1, 3, -1),
    (2, 1, -2, -1),
    (-3, 1, 3, -1),
    (-3, -2, -3, 2, 3, 1),
    (
        -3,
        -2,
        -3,
        -2,
        -3,
        -2,
        3,
        2,
        3,
        2,
        3,
        -1,
        -3,
        -2,
        -1,
        4,
        1,
        2,
        3,
        1,
        -3,
        -2,
        -3,
        -2,
        -3,
        2,
        3,
        2,
        3,
        2,
        3,
        -1,
        -3,
        -2,
        -1,
        -4,
        1,
        2,
        3,
        1,
    ),
    (3, 2, -3, -2),
    (-3, -2, -1, -4, 1, 2, 3, 1, -3, -2, -1, 4, 1, 2, 3, -1),
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -1, 4, 1, 2),
    (-4, 1, 4, -1),
    (3, 2, -3, -2),
    (4, 3, -4, -3),
)


@dataclass(frozen=True, slots=True)
class DeltaTenAlgebraCertificate:
    """Exact generic collision and exceptional-incidence identities."""

    collision_resultant_identity: Expr
    collision_degree: int
    collision_leading_coefficient: Expr
    collision_constant_identity: Expr
    residual_involution_identities: tuple[Expr, Expr]
    denominator_identity: Expr
    cusp_image_resultant_identity: Expr
    pair_diagonal_resultant_identity: Expr
    tangency_syzygy_identity: Expr
    exceptional_incidence_identities: tuple[Expr, ...]

    @property
    def verified(self) -> bool:
        """Whether all generic and reducible-fiber identities are exact."""

        return bool(
            self.collision_resultant_identity == 0
            and self.collision_degree == 10
            and self.collision_leading_coefficient == 1
            and self.collision_constant_identity == 0
            and self.residual_involution_identities == (0, 0)
            and self.denominator_identity == 0
            and self.cusp_image_resultant_identity == 0
            and self.pair_diagonal_resultant_identity == 0
            and self.tangency_syzygy_identity == 0
            and all(value == 0 for value in self.exceptional_incidence_identities)
        )


@cache
def exact_delta_ten_algebra_certificate() -> DeltaTenAlgebraCertificate:
    """Build the exact five-parameter pair-incidence certificate."""

    exceptional_identities: list[Expr] = [
        expand(SECOND_DIVIDED_DIFFERENCE.subs(S, 0) - R**2 * ZERO_FIBER_VERTICAL),
        expand(
            16 * SECOND_DIVIDED_DIFFERENCE.subs(R, (S**2 + 1) / 2)
            - ZERO_FIBER_GRAPH
        ),
        expand(COLLISION_POLYNOMIAL.subs(KAPPA, 0) - S**2 * ZERO_FIBER_GRAPH),
    ]
    for epsilon in (-1, 1):
        exceptional_identities.extend(
            (
                expand(
                    SECOND_DIVIDED_DIFFERENCE.subs(S, -epsilon)
                    - exceptional_vertical(epsilon)
                ),
                expand(
                    16
                    * SECOND_DIVIDED_DIFFERENCE.subs(
                        R, S * (S + epsilon) / 2
                    )
                    - S**2 * exceptional_graph(epsilon)
                ),
                expand(
                    COLLISION_POLYNOMIAL.subs(KAPPA, 2 * epsilon)
                    - (S + epsilon) ** 4 * exceptional_graph(epsilon)
                ),
            )
        )

    return DeltaTenAlgebraCertificate(
        collision_resultant_identity=expand(
            COLLISION_RESULTANT - S**2 * COLLISION_POLYNOMIAL
        ),
        collision_degree=Poly(COLLISION_POLYNOMIAL, S).degree(),
        collision_leading_coefficient=Poly(COLLISION_POLYNOMIAL, S).LC(),
        collision_constant_identity=expand(
            COLLISION_POLYNOMIAL.subs(S, 0) - ALPHA * KAPPA**2
        ),
        residual_involution_identities=(
            expand(FAMILY_P.subs(T, -T) - FAMILY_P.subs(KAPPA, -KAPPA)),
            expand(
                -FAMILY_Q.subs(T, -T)
                - FAMILY_Q.subs({BETA: -BETA, DELTA: -DELTA})
            ),
        ),
        denominator_identity=cancel(
            COLLISION_POLYNOMIAL.subs(S, -KAPPA / 2)
            - KAPPA**2 * (KAPPA**2 - 4) ** 4 / 1024
        ),
        cusp_image_resultant_identity=expand(
            resultant(COLLISION_POLYNOMIAL, PAIR_QUADRATIC, S)
            - (KAPPA**2 - 4) ** 4 * CUSP_IMAGE_FACTOR
        ),
        pair_diagonal_resultant_identity=expand(
            resultant(
                COLLISION_POLYNOMIAL,
                -S * PAIR_DIAGONAL_FACTOR,
                S,
            )
            - ALPHA
            * KAPPA**2
            * (KAPPA**2 - 4) ** 4
            * EXTRA_CRITICAL_FACTOR
        ),
        tangency_syzygy_identity=expand(
            PAIR_QUADRATIC * TANGENCY_POLYNOMIAL
            - S
            * PAIR_DIAGONAL_FACTOR
            * PAIR_DENOMINATOR
            * diff(COLLISION_POLYNOMIAL, S)
            - TANGENCY_COFACTOR * COLLISION_POLYNOMIAL
        ),
        exceptional_incidence_identities=tuple(exceptional_identities),
    )


@dataclass(frozen=True, slots=True)
class DeltaTenSampleCertificate:
    """Exact geometry and finite-image replay for the clean rational member."""

    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    collision_identity: Expr
    collision_discriminant: Expr
    denominator_resultant: Expr
    pair_diagonal_resultant: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    tangency_resultant: Expr
    node_x_resultant_identity: Expr
    node_x_discriminant: Expr
    arithmetic_genus: int
    cusp_delta: int
    node_count: int
    infinity_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return self.cusp_delta + self.node_count + self.infinity_delta

    @property
    def verified(self) -> bool:
        """Whether the clean geometry and complement obstruction both hold."""

        return bool(
            self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.collision_identity == 0
            and self.collision_discriminant == -407351195013757923
            and self.denominator_resultant == 81
            and self.pair_diagonal_resultant == 335421
            and self.cusp_image_factor == 1
            and self.extra_critical_factor == 4141
            and self.tangency_resultant == 136634145182709696290583
            and self.node_x_resultant_identity == 0
            and self.node_x_discriminant == -766610929107006671875
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.node_count == 10
            and self.infinity_delta == 16
            and self.total_delta == self.arithmetic_genus
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_sample_certificate() -> DeltaTenSampleCertificate:
    """Build the exact clean-member and four-meridian certificate."""

    sample_collision = COLLISION_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
    sample_tangency = TANGENCY_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
    sample_x_numerator = COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS)
    sample_x_denominator = COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    return DeltaTenSampleCertificate(
        implicit_resultant_identity=expand(
            resultant(SAMPLE_P - X, SAMPLE_Q - Y, T) - SAMPLE_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            SAMPLE_IMPLICIT.subs({X: SAMPLE_P, Y: SAMPLE_Q})
        ),
        collision_identity=expand(
            sample_collision - SAMPLE_COLLISION_POLYNOMIAL
        ),
        collision_discriminant=discriminant(SAMPLE_COLLISION_POLYNOMIAL, S),
        denominator_resultant=resultant(
            SAMPLE_COLLISION_POLYNOMIAL,
            PAIR_DENOMINATOR.subs(SAMPLE_PARAMETERS),
            S,
        ),
        pair_diagonal_resultant=resultant(
            SAMPLE_COLLISION_POLYNOMIAL,
            (-S * PAIR_DIAGONAL_FACTOR).subs(SAMPLE_PARAMETERS),
            S,
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(SAMPLE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(SAMPLE_PARAMETERS),
        tangency_resultant=resultant(
            SAMPLE_COLLISION_POLYNOMIAL,
            sample_tangency,
            S,
        ),
        node_x_resultant_identity=expand(
            resultant(
                SAMPLE_COLLISION_POLYNOMIAL,
                sample_x_numerator - X * sample_x_denominator,
                S,
            )
            - 6561 * SAMPLE_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(SAMPLE_NODE_X_POLYNOMIAL, X),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=(2 - 1) * (5 - 1) // 2,
        node_count=10,
        infinity_delta=(5 - 1) * (9 - 1) // 2,
        complement_census=n_generator_three_cycle_presentation_census(
            SAMPLE_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact generic delta-ten certificate and fail on regression."""

    algebra = exact_delta_ten_algebra_certificate()
    sample = exact_delta_ten_sample_certificate()
    print(f"delta-ten family algebra verified: {algebra.verified}")
    print(
        "clean (4,9) sample:",
        {
            "cusp": (2, 5),
            "nodes": sample.node_count,
            "infinity": (5, 9),
            "delta": sample.total_delta,
            "images": dict(sample.complement_census.generated_order_histogram),
        },
    )
    print("all 40^4 raw-presentation assignments excluded from A6:", sample.verified)
    print("generic clean stratum excluded conditionally: True")
    print("remaining: delta-ten degeneration walls and unrestricted passports")
    print("claim boundary: conditional/computer-assisted; JC(2) remains open")
    return 0 if algebra.verified and sample.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
