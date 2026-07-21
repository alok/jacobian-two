"""Exclude the generic conditional ``A6`` delta-seven family.

Under the one-pair and finite-singularity hypotheses used by
``a6_one_pair_infinity``, delta seven has one large-link survivor, with
affine degrees ``(3,10)``.  Polynomial target shears reduce that survivor to
a four-parameter family.  This module checks the family collision algebra
and one nondegenerate rational member exactly.  Sage supplies the member's
affine-complement presentation; the dependency-free finite-group replay
shows that it has no ``A6`` image with single-3-cycle meridians.

Proper Whitney--Thom equisingular propagation after a finite base change is
a theorem dependency, not an executable certificate here.  Consequently the
result is conditional and computer-assisted: it excludes the nondegenerate
open stratum, not its collision-degeneration walls, the unrestricted ``A6``
passport, or the plane Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import (
    Expr,
    Poly,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    factor_list,
    gcd,
    rem,
    resultant,
)

from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    three_cycle_presentation_census,
)

T: Final = Symbol("t")
U: Final = Symbol("u")
S: Final = Symbol("s")
Z: Final = Symbol("z")
X: Final = Symbol("X")
Y: Final = Symbol("Y")

ALPHA: Final = Symbol("alpha")
BETA: Final = Symbol("beta")
GAMMA: Final = Symbol("gamma")
DELTA: Final = Symbol("delta")

RAW_C4: Final = Symbol("c4")
RAW_C5: Final = Symbol("c5")
RAW_C6: Final = Symbol("c6")
RAW_C7: Final = Symbol("c7")
RAW_C8: Final = Symbol("c8")
RAW_C9: Final = Symbol("c9")

FAMILY_P: Final = T**2 + T**3
FAMILY_Q: Final = (
    ALPHA * T**5
    + BETA * T**7
    + GAMMA * T**8
    + DELTA * T**9
    + T**10
)

RAW_Q: Final = (
    RAW_C4 * T**4
    + RAW_C5 * T**5
    + RAW_C6 * T**6
    + RAW_C7 * T**7
    + RAW_C8 * T**8
    + RAW_C9 * T**9
    + T**10
)
NORMAL_ALPHA: Final = RAW_C5 - 2 * RAW_C4
NORMAL_BETA: Final = RAW_C7 - 3 * (RAW_C6 - RAW_C4)
NORMAL_GAMMA: Final = RAW_C8 - 3 * (RAW_C6 - RAW_C4)
NORMAL_DELTA: Final = RAW_C9 - (RAW_C6 - RAW_C4)
SHEARED_RAW_Q: Final = (
    RAW_Q
    - RAW_C4 * FAMILY_P**2
    - (RAW_C6 - RAW_C4) * FAMILY_P**3
)

# For a distinct collision put s=t+u and r=tu.  Equality of P gives
# r=s^2+s, and the second divided difference reduces to -s^2*H(s).
COLLISION_POLYNOMIAL: Final = (
    S**7
    + 6 * S**6
    + (9 - GAMMA + 3 * DELTA) * S**5
    + (-BETA - 2 * GAMMA + 9 * DELTA) * S**4
    + (-4 * BETA + 2 * GAMMA + 6 * DELTA - 5) * S**3
    + (ALPHA - 3 * BETA + 4 * GAMMA - DELTA) * S**2
    + (ALPHA + BETA) * S
    - ALPHA
)
PAIR_DISCRIMINANT: Final = -S * (3 * S + 4)
CUSP_COLLISION_FACTOR: Final = ALPHA + BETA - GAMMA + DELTA - 1
EXTRA_CRITICAL_FACTOR: Final = (
    1215 * ALPHA
    + 756 * BETA
    - 576 * GAMMA
    + 432 * DELTA
    - 320
)
TRIPLE_COLLISION_FACTOR: Final = (
    ALPHA**2
    - ALPHA * BETA * GAMMA
    + 3 * ALPHA * BETA * DELTA
    - ALPHA * BETA
    + ALPHA * GAMMA**2
    - 3 * ALPHA * GAMMA * DELTA
    + 7 * ALPHA * DELTA
    - 9 * ALPHA
    - 2 * BETA**2 * GAMMA
    + 6 * BETA**2 * DELTA
    - 6 * BETA**2
    + 5 * BETA * GAMMA**2
    - 20 * BETA * GAMMA * DELTA
    + 22 * BETA * GAMMA
    + 15 * BETA * DELTA**2
    - 32 * BETA * DELTA
    + 17 * BETA
    - 3 * GAMMA**3
    + 15 * GAMMA**2 * DELTA
    - 17 * GAMMA**2
    - 21 * GAMMA * DELTA**2
    + 46 * GAMMA * DELTA
    - 25 * GAMMA
    + 9 * DELTA**3
    - 29 * DELTA**2
    + 31 * DELTA
    - 11
)

# After dividing the velocity determinant by t-u and reducing by the
# P-collision relation, the exact answer is -s^2*(s+1)*D(s).  The factors
# s=0 and s=-1 are invalid diagonal/cusp collisions on the open stratum.
TANGENCY_POLYNOMIAL: Final = (
    15 * ALPHA * S**2
    + 15 * ALPHA * S
    - 10 * ALPHA
    - 21 * BETA * S**4
    - 77 * BETA * S**3
    - 56 * BETA * S**2
    + 14 * BETA * S
    - 24 * GAMMA * S**5
    - 48 * GAMMA * S**4
    + 32 * GAMMA * S**3
    + 72 * GAMMA * S**2
    + 72 * DELTA * S**5
    + 207 * DELTA * S**4
    + 135 * DELTA * S**3
    - 18 * DELTA * S**2
    + 30 * S**7
    + 170 * S**6
    + 250 * S**5
    + 20 * S**4
    - 110 * S**3
)

CLEAN_PARAMETERS: Final = {
    ALPHA: 2,
    BETA: 0,
    GAMMA: 0,
    DELTA: 0,
}
CLEAN_Q: Final = FAMILY_Q.subs(CLEAN_PARAMETERS)
CLEAN_H: Final = COLLISION_POLYNOMIAL.subs(CLEAN_PARAMETERS)
CLEAN_TANGENCY: Final = TANGENCY_POLYNOMIAL.subs(CLEAN_PARAMETERS)
CLEAN_COLLISION_FACTOR: Final = (
    (T**4 + T**3 + T**2 + T + 1)
    * (
        T**10
        + 5 * T**9
        + 15 * T**8
        + 15 * T**7
        + 5 * T**6
        + 2 * T**5
        + 10 * T**3
        + 10 * T**2
        + 2
    )
)
CLEAN_NODE_COORDINATE_POLYNOMIAL: Final = (
    X**7
    - 24 * X**6
    - 31 * X**5
    - 30 * X**4
    - 65 * X**3
    + 28 * X**2
    + 18 * X
    + 2
)
CLEAN_BRANCH_EQUATION: Final = (
    -X**10
    - 10 * X**8
    + 15 * X**6 * Y
    - 20 * X**6
    - 4 * X**5 * Y
    - 4 * X**5
    + 50 * X**4 * Y
    + 10 * X**3 * Y**2
    + 10 * X**3 * Y
    - 25 * X**2 * Y**2
    + Y**3
    + Y**2
)

# Sage 10.8, exact affine Zariski--van Kamp presentation for CLEAN_BRANCH.
# The three generators are geometric fiber meridians.
GENERIC_RELATIONS: Final = (
    (2, 1, -2, -1),
    (2, 1, -2, -1),
    (-3, 1, 3, -1),
    (
        -3,
        -2,
        3,
        1,
        -3,
        2,
        3,
        1,
        -3,
        2,
        3,
        1,
        -3,
        -2,
        3,
        -1,
        -3,
        -2,
        3,
        -1,
    ),
    (3, 2, -3, -2),
    (
        -2,
        -1,
        -2,
        -1,
        3,
        1,
        2,
        1,
        2,
        -1,
        -2,
        -1,
        -3,
        1,
        2,
        1,
    ),
    (-2, -1, 3, 1),
    (2, 1, -2, -1),
    (3, 2, -3, -2),
)


@dataclass(frozen=True, slots=True)
class A6DeltaSevenCleanCertificate:
    """Exact geometry and complement-image census for the clean member."""

    implicit_resultant_identity: Expr
    parametrized_branch_identity: Expr
    off_diagonal_gcd: Expr
    collision_resultant_identity: Expr
    pairing_resultant_identity: Expr
    collision_discriminant: Expr
    pair_diagonal_resultant: Expr
    pair_product_resultant: Expr
    tangency_resultant: Expr
    node_coordinate_resultant_identity: Expr
    node_coordinate_discriminant: Expr
    cusp_coordinate_separation: Expr
    cusp_image_separation: Expr
    residual_derivative_value: Expr
    finite_cusp_pair: tuple[int, int]
    node_count: int
    infinity_pair: tuple[int, int]
    arithmetic_genus: int
    total_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether every clean-member checkpoint agrees exactly."""

        return bool(
            self.implicit_resultant_identity == 0
            and self.parametrized_branch_identity == 0
            and self.off_diagonal_gcd == 1
            and self.collision_resultant_identity == 0
            and self.pairing_resultant_identity == 0
            and self.collision_discriminant == 464_000_000
            and self.pair_diagonal_resultant == -4220
            and self.pair_product_resultant == 2
            and self.tangency_resultant == -1_958_080_000_000
            and self.node_coordinate_resultant_identity == 0
            and self.node_coordinate_discriminant
            == 113_281_250_000_000_000
            and self.cusp_coordinate_separation == 2
            and self.cusp_image_separation == -1
            and self.residual_derivative_value == Rational(33760, 19683)
            and self.finite_cusp_pair == (2, 5)
            and self.node_count == 7
            and self.infinity_pair == (7, 10)
            and self.arithmetic_genus == 36
            and self.total_delta == 36
            and self.complement_census.assignments == 40**3
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_a6_delta_seven_clean_certificate(
) -> A6DeltaSevenCleanCertificate:
    """Build the exact nondegenerate ``(3,10)`` representative certificate."""

    parameter_p_u = FAMILY_P.subs(T, U)
    clean_q_u = CLEAN_Q.subs(T, U)
    first_collision = cancel((FAMILY_P - parameter_p_u) / (T - U))
    second_collision = cancel((CLEAN_Q - clean_q_u) / (T - U))
    return A6DeltaSevenCleanCertificate(
        implicit_resultant_identity=expand(
            resultant(X - FAMILY_P, Y - CLEAN_Q, T)
            - CLEAN_BRANCH_EQUATION
        ),
        parametrized_branch_identity=expand(
            CLEAN_BRANCH_EQUATION.subs({X: FAMILY_P, Y: CLEAN_Q})
        ),
        off_diagonal_gcd=gcd(first_collision, second_collision),
        collision_resultant_identity=expand(
            resultant(first_collision, second_collision, U)
            - T**4 * CLEAN_COLLISION_FACTOR
        ),
        pairing_resultant_identity=expand(
            resultant(
                CLEAN_H,
                T**2 - S * T + S**2 + S,
                S,
            )
            - CLEAN_COLLISION_FACTOR
        ),
        collision_discriminant=discriminant(CLEAN_H, S),
        pair_diagonal_resultant=resultant(
            CLEAN_H,
            PAIR_DISCRIMINANT,
            S,
        ),
        pair_product_resultant=resultant(CLEAN_H, S**2 + S, S),
        tangency_resultant=resultant(CLEAN_H, CLEAN_TANGENCY, S),
        node_coordinate_resultant_identity=expand(
            resultant(
                CLEAN_COLLISION_FACTOR,
                X - FAMILY_P,
                T,
            )
            - CLEAN_NODE_COORDINATE_POLYNOMIAL**2
        ),
        node_coordinate_discriminant=discriminant(
            CLEAN_NODE_COORDINATE_POLYNOMIAL,
            X,
        ),
        cusp_coordinate_separation=CLEAN_NODE_COORDINATE_POLYNOMIAL.subs(X, 0),
        cusp_image_separation=CLEAN_Q.subs(T, -1),
        residual_derivative_value=diff(CLEAN_Q, T).subs(
            T,
            Rational(-2, 3),
        ),
        finite_cusp_pair=(2, 5),
        node_count=7,
        infinity_pair=(10 - 3, 10),
        arithmetic_genus=(10 - 1) * (10 - 2) // 2,
        total_delta=2 + 7 + (7 - 1) * (10 - 1) // 2,
        complement_census=three_cycle_presentation_census(GENERIC_RELATIONS),
    )


@dataclass(frozen=True, slots=True)
class A6DeltaSevenFamilyCertificate:
    """Exact four-parameter collision algebra plus the clean member."""

    normal_form_identity: Expr
    collision_reduction_remainder: Expr
    tangency_reduction_remainder: Expr
    pair_resultant_identity: Expr
    derivative_factor_identity: Expr
    cusp_collision_identity: Expr
    triple_resultant_identity: Expr
    discriminant_factor_degrees: tuple[tuple[int, int], ...]
    tangency_resultant_identity: Expr
    clean_open_values: tuple[Expr, Expr, Expr, Expr, Expr]
    clean: A6DeltaSevenCleanCertificate

    @property
    def verified(self) -> bool:
        """Whether the family and representative identities all agree."""

        return bool(
            self.normal_form_identity == 0
            and self.collision_reduction_remainder == 0
            and self.tangency_reduction_remainder == 0
            and self.pair_resultant_identity == 0
            and self.derivative_factor_identity == 0
            and self.cusp_collision_identity == 0
            and self.triple_resultant_identity == 0
            and self.discriminant_factor_degrees == ((1, 1), (10, 1))
            and self.tangency_resultant_identity == 0
            and self.clean_open_values
            == (2, 2110, 1, -25, 464_000_000)
            and all(value != 0 for value in self.clean_open_values)
            and self.clean.verified
        )


@cache
def exact_a6_delta_seven_family_certificate(
) -> A6DeltaSevenFamilyCertificate:
    """Build the exact generic-family certificate.

    If ``Disc(H)=C*G``, the nondegenerate algebraic open used for the
    topology propagation is ``alpha*L*C*T*G != 0``.  The executable part
    proves all displayed identities and that the clean member lies in it.
    """

    parameter_p_u = FAMILY_P.subs(T, U)
    family_q_u = FAMILY_Q.subs(T, U)
    first_collision = cancel((FAMILY_P - parameter_p_u) / (T - U))
    second_collision = cancel((FAMILY_Q - family_q_u) / (T - U))
    collision_remainder = rem(
        second_collision
        + (T + U) ** 2 * COLLISION_POLYNOMIAL.subs(S, T + U),
        first_collision,
        U,
    )

    velocity_quotient = cancel(
        (
            diff(FAMILY_P, T) * diff(family_q_u, U)
            - diff(parameter_p_u, U) * diff(FAMILY_Q, T)
        )
        / (T - U)
    )
    tangency_remainder = rem(
        velocity_quotient
        + (T + U) ** 2
        * (T + U + 1)
        * TANGENCY_POLYNOMIAL.subs(S, T + U),
        first_collision,
        U,
    )

    reduced_q = rem(
        FAMILY_Q.subs(T, Z),
        Z**3 + Z**2 - X,
        Z,
    )
    reduced_poly = Poly(reduced_q, Z)
    quadratic_coefficient = reduced_poly.coeff_monomial(Z**2)
    linear_coefficient = reduced_poly.coeff_monomial(Z)

    family_discriminant = discriminant(COLLISION_POLYNOMIAL, S)
    reduced_discriminant = cancel(
        family_discriminant / CUSP_COLLISION_FACTOR
    )
    factorization = factor_list(
        family_discriminant,
        ALPHA,
        BETA,
        GAMMA,
        DELTA,
    )
    factor_degrees = tuple(
        sorted(
            (
                int(Poly(factor, ALPHA, BETA, GAMMA, DELTA).total_degree()),
                exponent,
            )
            for factor, exponent in factorization[1]
        )
    )

    normal_substitution = {
        ALPHA: NORMAL_ALPHA,
        BETA: NORMAL_BETA,
        GAMMA: NORMAL_GAMMA,
        DELTA: NORMAL_DELTA,
    }
    clean_substitution = CLEAN_PARAMETERS
    return A6DeltaSevenFamilyCertificate(
        normal_form_identity=expand(
            SHEARED_RAW_Q - FAMILY_Q.subs(normal_substitution)
        ),
        collision_reduction_remainder=expand(collision_remainder),
        tangency_reduction_remainder=expand(tangency_remainder),
        pair_resultant_identity=expand(
            resultant(
                COLLISION_POLYNOMIAL,
                PAIR_DISCRIMINANT,
                S,
            )
            + ALPHA * EXTRA_CRITICAL_FACTOR
        ),
        derivative_factor_identity=expand(
            diff(FAMILY_Q, T).subs(T, Rational(-2, 3))
            - 16 * EXTRA_CRITICAL_FACTOR / 19683
        ),
        cusp_collision_identity=expand(
            COLLISION_POLYNOMIAL.subs(S, -1)
            + CUSP_COLLISION_FACTOR
        ),
        triple_resultant_identity=expand(
            resultant(quadratic_coefficient, linear_coefficient, X)
            - CUSP_COLLISION_FACTOR**2 * TRIPLE_COLLISION_FACTOR
        ),
        discriminant_factor_degrees=factor_degrees,
        tangency_resultant_identity=expand(
            resultant(
                COLLISION_POLYNOMIAL,
                TANGENCY_POLYNOMIAL,
                S,
            )
            + ALPHA * EXTRA_CRITICAL_FACTOR * reduced_discriminant
        ),
        clean_open_values=(
            ALPHA.subs(clean_substitution),
            EXTRA_CRITICAL_FACTOR.subs(clean_substitution),
            CUSP_COLLISION_FACTOR.subs(clean_substitution),
            TRIPLE_COLLISION_FACTOR.subs(clean_substitution),
            reduced_discriminant.subs(clean_substitution),
        ),
        clean=exact_a6_delta_seven_clean_certificate(),
    )


def main() -> int:
    """Print the generic delta-seven reduction and fail on regression."""

    certificate = exact_a6_delta_seven_family_certificate()
    print(
        "clean conditional (3,10) curve:",
        {
            "cusp": certificate.clean.finite_cusp_pair,
            "nodes": certificate.clean.node_count,
            "infinity": certificate.clean.infinity_pair,
            "delta": certificate.clean.total_delta,
            "images": dict(
                certificate.clean.complement_census.generated_order_histogram
            ),
        },
    )
    print("nondegenerate open values alpha,L,C,T,G:", certificate.clean_open_values)
    print(f"generic delta-seven certificate: {certificate.verified}")
    print("excluded conditionally: the nondegenerate (3,10) open stratum")
    print("remaining: collision-degeneration walls and unrestricted passports")
    print("claim boundary: conditional/computer-assisted; JC(2) remains open")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
