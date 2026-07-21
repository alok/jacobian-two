"""Exact exclusion of the three one-dimensional delta-seven wall strata.

After the generic discriminant divisor and its two generic singular
components, the collision septic has excess three.  There are exactly three
root partitions:

* ``(4,1,1,1)`` -- one contact-four point and three nodes;
* ``(3,2,1,1)`` -- one contact-three point, one contact-two point, two nodes;
* ``(2,2,2,1)`` -- three contact-two points and one node.

The first and third strata are rational curves.  The middle stratum is an
absolutely irreducible genus-two curve; treating it as rational would be an
error.  This module stores exact normalizations/incidence equations, valid
rational representatives, certified Sage 10.8 van Kampen presentations, and
dependency-free replays of all ``40^3`` single-three-cycle assignments into
``A6``.  Every presentation has only 40 cyclic ``C3`` images and no ``A6``
image.

The propagation from one representative to its connected equisingular
stratum uses the same proper projective Whitney--Thom theorem dependency as
the preceding wall audit.  All conclusions remain conditional on the
one-pair and finite-singularity hypotheses; this is not a proof of ``JC(2)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import (
    Expr,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
)

from scripts.a6_delta_seven_discriminant_wall import (
    COEFFICIENT_SLICE_FIRST,
    COEFFICIENT_SLICE_INVERSE,
    COEFFICIENT_SLICE_SECOND,
    H0,
    H1,
    H2,
    H3,
    H4,
    H5,
)
from scripts.a6_delta_seven_generic import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    CUSP_COLLISION_FACTOR,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    S,
    TRIPLE_COLLISION_FACTOR,
    X,
    Y,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    three_cycle_presentation_census,
)

V: Final = Symbol("v")
MH: Final = Symbol("h")
MK: Final = Symbol("k")
MQ: Final = Symbol("q")
E: Final = Symbol("e")


def _collision_coefficients(polynomial: Expr) -> dict[Symbol, Expr]:
    """Extract ``h0`` through ``h5`` from a monic collision septic."""

    expanded = expand(polynomial)
    return {
        H0: expanded.coeff(S, 0),
        H1: expanded.coeff(S, 1),
        H2: expanded.coeff(S, 2),
        H3: expanded.coeff(S, 3),
        H4: expanded.coeff(S, 4),
        H5: expanded.coeff(S, 5),
    }


def _parameters_from_collision(polynomial: Expr) -> dict[Symbol, Expr]:
    """Recover the four normal-form parameters from a slice septic."""

    coefficients = _collision_coefficients(polynomial)
    return {
        parameter: cancel(expression.subs(coefficients))
        for parameter, expression in COEFFICIENT_SLICE_INVERSE.items()
    }


# ---------------------------------------------------------------------------
# The rational contact-four stratum (4,1,1,1)
# ---------------------------------------------------------------------------

CONTACT_FOUR_PARAMETERS: Final = {
    ALPHA: -V**4 * (V + 1) * (15 * V**2 + 68 * V + 85) / (3 * V + 13),
    BETA: 5
    * V**3
    * (9 * V**4 + 39 * V**3 + 9 * V**2 - 153 * V - 136)
    / (2 * (3 * V + 13)),
    GAMMA: 5
    * (
        9 * V**7
        + 39 * V**6
        + 9 * V**5
        - 141 * V**4
        - 75 * V**3
        + 51 * V**2
        - 5 * V
        + 1
    )
    / (2 * (3 * V + 13)),
    DELTA: 5
    * (
        3 * V**7
        + 13 * V**6
        + 3 * V**5
        - 47 * V**4
        - 28 * V**3
        - 20 * V
        + 4
    )
    / (2 * (3 * V + 13)),
}

CONTACT_FOUR_RESIDUAL_FACTOR: Final = (
    S**3
    + 2 * (2 * V + 3) * S**2
    + (15 * V**3 + 149 * V**2 + 403 * V + 289)
    / (2 * (3 * V + 13))
    * S
    + (V + 1) * (15 * V**2 + 68 * V + 85) / (3 * V + 13)
)
CONTACT_FOUR_FACTORIZATION: Final = (S - V) ** 4 * CONTACT_FOUR_RESIDUAL_FACTOR

CONTACT_FOUR_C_FACTOR: Final = -(
    (V + 1) ** 6 * (15 * V + 11) / (2 * (3 * V + 13))
)
CONTACT_FOUR_L_FACTOR: Final = -(
    5 * (3 * V + 1) ** 3 * (3 * V + 4) ** 4 / (3 * V + 13)
)
CONTACT_FOUR_T_FACTOR: Final = -(
    (V + 1) ** 3
    * (3 * V + 1) ** 6
    * (3 * V + 11) ** 4
    / (8 * (3 * V + 13) ** 3)
)
CONTACT_FOUR_DERIVATIVE_FACTOR: Final = (
    60 * (V + 1) ** 2 * (9 * V**2 + 51 * V + 34) / (3 * V + 13)
)
CONTACT_FOUR_RESIDUAL_DISCRIMINANT: Final = (
    5
    * (V + 1) ** 2
    * (3 * V + 1) ** 2
    * (15 * V + 11)
    * (3 * V**4 + 28 * V**3 + 80 * V**2 + 68 * V + 17)
    / (2 * (3 * V + 13) ** 3)
)

CONTACT_FOUR_SAMPLE_PARAMETERS: Final = {
    ALPHA: -21,
    BETA: Rational(-145, 4),
    GAMMA: Rational(-35, 2),
    DELTA: Rational(-45, 4),
}
CONTACT_FOUR_SAMPLE_FACTORIZATION: Final = (
    (S - 1) ** 4 * (2 * S + 3) * (2 * S**2 + 17 * S + 28) / 4
)
CONTACT_FOUR_BRANCH_EQUATION: Final = (
    64 * X**10
    - 146225 * X**9
    + 868720 * X**8
    - 833440 * X**7
    - 36540 * X**6 * Y
    - 3091200 * X**6
    + 250320 * X**5 * Y
    - 1467648 * X**5
    - 257600 * X**4 * Y
    - 2800 * X**3 * Y**2
    - 349440 * X**3 * Y
    + 17360 * X**2 * Y**2
    - 21120 * X * Y**2
    - 64 * Y**3
    + 3328 * Y**2
)
CONTACT_FOUR_RELATIONS: Final = (
    (2, 1, 2, 1, 2, 1, 2, 1, -2, -1, -2, -1, -2, -1, -2, -1),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, -2, -1, -2, -1, 3, 1, 2, 1, 2, 1),
    (-2, -1, -2, -1, 3, 1, 2, 1, 2, -1, -2, -1, -3, 1, 2, 1),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (3, 2, -3, -2),
)


# ---------------------------------------------------------------------------
# The genus-two contact-three/contact-two stratum (3,2,1,1)
# ---------------------------------------------------------------------------

MIXED_Q_DENOMINATOR: Final = (
    3 * MH**2 + 6 * MH * MK + 12 * MH + MK**2 + 8 * MK + 10
)
MIXED_Q_NUMERATOR: Final = (
    3 * MH**4
    + 18 * MH**3 * MK
    + 38 * MH**3
    + 18 * MH**2 * MK**2
    + 108 * MH**2 * MK
    + 132 * MH**2
    + 6 * MH * MK**3
    + 66 * MH * MK**2
    + 204 * MH * MK
    + 180 * MH
    + 8 * MK**3
    + 54 * MK**2
    + 120 * MK
    + 85
)
MIXED_Q: Final = MIXED_Q_NUMERATOR / MIXED_Q_DENOMINATOR
MIXED_PARAMETER_CURVE: Final = (
    3 * MH**4 * MK
    + 9 * MH**4
    + 18 * MH**3 * MK**2
    + 83 * MH**3 * MK
    + 87 * MH**3
    + 18 * MH**2 * MK**3
    + 168 * MH**2 * MK**2
    + 393 * MH**2 * MK
    + 273 * MH**2
    + 6 * MH * MK**4
    + 108 * MH * MK**3
    + 444 * MH * MK**2
    + 681 * MH * MK
    + 357 * MH
    + 22 * MK**4
    + 146 * MK**3
    + 366 * MK**2
    + 408 * MK
    + 170
)
MIXED_COLLISION_ANSATZ: Final = (
    (S - MH) ** 3
    * (S - MK) ** 2
    * (S**2 + (6 + 3 * MH + 2 * MK) * S + MQ)
)
MIXED_COLLISION_FAMILY: Final = cancel(MIXED_COLLISION_ANSATZ.subs(MQ, MIXED_Q))
MIXED_PARAMETERS: Final = _parameters_from_collision(MIXED_COLLISION_FAMILY)

MIXED_SAMPLE_SUBSTITUTION: Final = {
    MH: Rational(-2, 3),
    MK: Rational(-11, 9),
    MQ: Rational(1, 9),
}
MIXED_SAMPLE_PARAMETERS: Final = {
    ALPHA: Rational(-968, 19683),
    BETA: Rational(6820, 6561),
    GAMMA: Rational(18745, 6561),
    DELTA: Rational(56410, 19683),
}
MIXED_SAMPLE_FACTORIZATION: Final = (
    (3 * S + 2) ** 3
    * (9 * S + 11) ** 2
    * (9 * S**2 + 14 * S + 1)
    / 19683
)
MIXED_Y_SCALE: Final = 19683
MIXED_BRANCH_EQUATION: Final = (
    -7625597484987 * X**10
    + 3094601589800 * X**9
    - 472284263760 * X**8
    + 33943539520 * X**7
    - 176547210 * X**6 * Y
    - 1158502400 * X**6
    + 43414680 * X**5 * Y
    + 14992384 * X**5
    - 3202400 * X**4 * Y
    + 27600 * X**3 * Y**2
    + 77440 * X**3 * Y
    - 8295 * X**2 * Y**2
    + 640 * X * Y**2
    + Y**3
    - 16 * Y**2
)
MIXED_RELATIONS: Final = (
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, 3, 1, 2, -1, 3, 1, 2, -1, -3, 1, -2, -1, -3, 1),
    (-2, -1, -3, 1, 3, 1, 2, -1, -3, -1, 3, 1),
    (-3, 1, 3, 1, 3, 1, 3, -1, -3, -1, -3, -1),
    (3, -2),
    (3, 2, -3, -2),
)


# ---------------------------------------------------------------------------
# The rational three-contact-two stratum (2,2,2,1)
# ---------------------------------------------------------------------------

THREE_DOUBLE_CUBIC: Final = (
    S**3
    - E * S**2
    + (4 * E**2 + 16 * E + 17) * S
    + 12 * E**2
    + 57 * E
    + 68
)
THREE_DOUBLE_FACTORIZATION: Final = THREE_DOUBLE_CUBIC**2 * (S + 6 + 2 * E)
THREE_DOUBLE_PARAMETERS: Final = {
    ALPHA: -2 * (E + 3) * (12 * E**2 + 57 * E + 68) ** 2,
    BETA: 5
    * (12 * E**2 + 57 * E + 68)
    * (8 * E**3 + 62 * E**2 + 159 * E + 136),
    GAMMA: 5
    * (
        96 * E**5
        + 1200 * E**4
        + 5988 * E**3
        + 14930 * E**2
        + 18620 * E
        + 9301
    ),
    DELTA: 5
    * (
        32 * E**5
        + 400 * E**4
        + 1996 * E**3
        + 4977 * E**2
        + 6208 * E
        + 3102
    ),
}
THREE_DOUBLE_C_FACTOR: Final = -4 * (2 * E + 5) ** 5
THREE_DOUBLE_L_FACTOR: Final = -50 * (
    (3 * E + 7) * (3 * E + 8) ** 2 * (12 * E + 29) ** 2
)
THREE_DOUBLE_T_FACTOR: Final = -4 * (E + 2) ** 8 * (12 * E + 29) ** 2
THREE_DOUBLE_CUBIC_DISCRIMINANT: Final = -20 * (2 * E + 5) ** 2 * (
    3 * E**4 + 32 * E**3 + 151 * E**2 + 340 * E + 289
)
THREE_DOUBLE_SEPARATION: Final = -5 * (E + 2) * (2 * E + 5) ** 2

THREE_DOUBLE_SAMPLE_PARAMETERS: Final = {
    ALPHA: Rational(-3, 8),
    BETA: Rational(5, 2),
    GAMMA: Rational(145, 32),
    DELTA: Rational(105, 32),
}
THREE_DOUBLE_SAMPLE_FACTORIZATION: Final = (
    (2 * S + 3) * (4 * S**3 + 9 * S**2 + 5 * S + 2) ** 2 / 32
)
THREE_DOUBLE_Y_SCALE: Final = 32
THREE_DOUBLE_BRANCH_EQUATION: Final = (
    -32768 * X**10
    + 114055 * X**9
    - 131320 * X**8
    + 55840 * X**7
    - 2445 * X**6 * Y
    - 9840 * X**6
    + 4970 * X**5 * Y
    + 576 * X**5
    - 1980 * X**4 * Y
    + 5 * X**3 * Y**2
    + 240 * X**3 * Y
    - 90 * X**2 * Y**2
    + 35 * X * Y**2
    + Y**3
    - 4 * Y**2
)
THREE_DOUBLE_RELATIONS: Final = (
    (-3, 1, 3, 1, 3, 1, -3, -1, -3, -1),
    (
        -3,
        -1,
        -3,
        -2,
        3,
        1,
        3,
        1,
        -3,
        -1,
        -3,
        2,
        3,
        1,
        3,
        1,
        -3,
        -1,
        -3,
        2,
        3,
        1,
        3,
        -1,
        -3,
        -1,
        -3,
        -2,
        3,
        1,
        3,
        -1,
    ),
    (-3, -1, -3, 2, 3, 1, 3, -1, -3, -2, 3, 1),
    (-3, -1, -3, 2, 3, 1),
    (-3, -2, 3, 1, -3, 2, 3, 1, -3, 2, 3, -1, -3, -2, 3, -1),
    (3, 2, 3, 2, -3, -2, -3, -2),
)


@dataclass(frozen=True, slots=True)
class ExcessThreeAlgebraCertificate:
    """Exact incidence identities for all three excess-three partitions."""

    contact_four_factor_identity: Expr
    contact_four_c_identity: Expr
    contact_four_l_identity: Expr
    contact_four_t_identity: Expr
    contact_four_derivative_identity: Expr
    contact_four_residual_discriminant_identity: Expr
    mixed_first_slice_identity: Expr
    mixed_second_slice_identity: Expr
    mixed_sample_parameter_identities: tuple[Expr, ...]
    three_double_factor_identity: Expr
    three_double_c_identity: Expr
    three_double_l_identity: Expr
    three_double_t_identity: Expr
    three_double_discriminant_identity: Expr
    three_double_separation_identity: Expr

    @property
    def verified(self) -> bool:
        """Whether every normalization identity vanishes exactly."""

        return all(
            value == 0
            for value in (
                self.contact_four_factor_identity,
                self.contact_four_c_identity,
                self.contact_four_l_identity,
                self.contact_four_t_identity,
                self.contact_four_derivative_identity,
                self.contact_four_residual_discriminant_identity,
                self.mixed_first_slice_identity,
                self.mixed_second_slice_identity,
                *self.mixed_sample_parameter_identities,
                self.three_double_factor_identity,
                self.three_double_c_identity,
                self.three_double_l_identity,
                self.three_double_t_identity,
                self.three_double_discriminant_identity,
                self.three_double_separation_identity,
            )
        )


@cache
def exact_excess_three_algebra_certificate() -> ExcessThreeAlgebraCertificate:
    """Build the three exact one-dimensional wall normalizations."""

    contact_four_h = COLLISION_POLYNOMIAL.subs(CONTACT_FOUR_PARAMETERS)
    mixed_coefficients = _collision_coefficients(MIXED_COLLISION_ANSATZ)
    mixed_sample_parameters = {
        parameter: cancel(expression.subs(MIXED_SAMPLE_SUBSTITUTION))
        for parameter, expression in MIXED_PARAMETERS.items()
    }
    three_double_h = COLLISION_POLYNOMIAL.subs(THREE_DOUBLE_PARAMETERS)
    simple_root = -6 - 2 * E
    return ExcessThreeAlgebraCertificate(
        contact_four_factor_identity=cancel(
            contact_four_h - CONTACT_FOUR_FACTORIZATION
        ),
        contact_four_c_identity=cancel(
            CUSP_COLLISION_FACTOR.subs(CONTACT_FOUR_PARAMETERS)
            - CONTACT_FOUR_C_FACTOR
        ),
        contact_four_l_identity=cancel(
            EXTRA_CRITICAL_FACTOR.subs(CONTACT_FOUR_PARAMETERS)
            - CONTACT_FOUR_L_FACTOR
        ),
        contact_four_t_identity=cancel(
            TRIPLE_COLLISION_FACTOR.subs(CONTACT_FOUR_PARAMETERS)
            - CONTACT_FOUR_T_FACTOR
        ),
        contact_four_derivative_identity=cancel(
            diff(contact_four_h, S, 4).subs(S, V)
            - CONTACT_FOUR_DERIVATIVE_FACTOR
        ),
        contact_four_residual_discriminant_identity=cancel(
            discriminant(CONTACT_FOUR_RESIDUAL_FACTOR, S)
            - CONTACT_FOUR_RESIDUAL_DISCRIMINANT
        ),
        mixed_first_slice_identity=cancel(
            COEFFICIENT_SLICE_FIRST.subs(mixed_coefficients).subs(MQ, MIXED_Q)
        ),
        mixed_second_slice_identity=cancel(
            COEFFICIENT_SLICE_SECOND.subs(mixed_coefficients).subs(MQ, MIXED_Q)
            + (MH + 1) ** 3
            * (MK + 1)
            * MIXED_PARAMETER_CURVE
            / MIXED_Q_DENOMINATOR
        ),
        mixed_sample_parameter_identities=tuple(
            cancel(mixed_sample_parameters[parameter] - expected)
            for parameter, expected in MIXED_SAMPLE_PARAMETERS.items()
        ),
        three_double_factor_identity=expand(
            three_double_h - THREE_DOUBLE_FACTORIZATION
        ),
        three_double_c_identity=expand(
            CUSP_COLLISION_FACTOR.subs(THREE_DOUBLE_PARAMETERS)
            - THREE_DOUBLE_C_FACTOR
        ),
        three_double_l_identity=expand(
            EXTRA_CRITICAL_FACTOR.subs(THREE_DOUBLE_PARAMETERS)
            - THREE_DOUBLE_L_FACTOR
        ),
        three_double_t_identity=expand(
            TRIPLE_COLLISION_FACTOR.subs(THREE_DOUBLE_PARAMETERS)
            - THREE_DOUBLE_T_FACTOR
        ),
        three_double_discriminant_identity=expand(
            discriminant(THREE_DOUBLE_CUBIC, S)
            - THREE_DOUBLE_CUBIC_DISCRIMINANT
        ),
        three_double_separation_identity=expand(
            THREE_DOUBLE_CUBIC.subs(S, simple_root) - THREE_DOUBLE_SEPARATION
        ),
    )


@dataclass(frozen=True, slots=True)
class ExcessThreeSampleCertificate:
    """Geometry and exact three-cycle replay for a stratum representative."""

    name: str
    partition: tuple[int, ...]
    contact_orders: tuple[int, ...]
    node_count: int
    collision_identity: Expr
    implicit_identity: Expr
    cusp_factor: Expr
    critical_factor: Expr
    triple_factor: Expr
    residual_separation: Expr
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the projective genus-budget contribution."""

        return 27 + 2 + sum(self.contact_orders) + self.node_count

    @property
    def verified(self) -> bool:
        """Whether exact geometry and finite-image replay agree."""

        return bool(
            self.collision_identity == 0
            and self.implicit_identity == 0
            and self.cusp_factor != 0
            and self.critical_factor != 0
            and self.triple_factor != 0
            and self.residual_separation != 0
            and self.total_delta == 36
            and self.complement_census.assignments == 40**3
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


def _sample_certificate(
    *,
    name: str,
    parameters: dict[Symbol, Expr | int],
    factorization: Expr,
    branch_equation: Expr,
    y_scale: int,
    relations: tuple[tuple[int, ...], ...],
    partition: tuple[int, ...],
    contact_orders: tuple[int, ...],
    node_count: int,
    residual_separation: Expr,
) -> ExcessThreeSampleCertificate:
    """Build one exact excess-three representative certificate."""

    sample_q = FAMILY_Q.subs(parameters)
    return ExcessThreeSampleCertificate(
        name=name,
        partition=partition,
        contact_orders=contact_orders,
        node_count=node_count,
        collision_identity=expand(
            COLLISION_POLYNOMIAL.subs(parameters) - factorization
        ),
        implicit_identity=expand(
            branch_equation.subs({X: FAMILY_P, Y: y_scale * sample_q})
        ),
        cusp_factor=expand(CUSP_COLLISION_FACTOR.subs(parameters)),
        critical_factor=expand(EXTRA_CRITICAL_FACTOR.subs(parameters)),
        triple_factor=expand(TRIPLE_COLLISION_FACTOR.subs(parameters)),
        residual_separation=residual_separation,
        complement_census=three_cycle_presentation_census(relations),
    )


@cache
def exact_contact_four_sample_certificate() -> ExcessThreeSampleCertificate:
    """Return the exact contact-four representative."""

    residual = (2 * S + 3) * (2 * S**2 + 17 * S + 28) / 4
    return _sample_certificate(
        name="one contact four plus three nodes",
        parameters=CONTACT_FOUR_SAMPLE_PARAMETERS,
        factorization=CONTACT_FOUR_SAMPLE_FACTORIZATION,
        branch_equation=CONTACT_FOUR_BRANCH_EQUATION,
        y_scale=1,
        relations=CONTACT_FOUR_RELATIONS,
        partition=(4, 1, 1, 1),
        contact_orders=(4,),
        node_count=3,
        residual_separation=discriminant(residual, S) * residual.subs(S, 1),
    )


@cache
def exact_mixed_sample_certificate() -> ExcessThreeSampleCertificate:
    """Return the exact contact-three/contact-two representative."""

    residual = 9 * S**2 + 14 * S + 1
    return _sample_certificate(
        name="contact three plus contact two plus two nodes",
        parameters=MIXED_SAMPLE_PARAMETERS,
        factorization=MIXED_SAMPLE_FACTORIZATION,
        branch_equation=MIXED_BRANCH_EQUATION,
        y_scale=MIXED_Y_SCALE,
        relations=MIXED_RELATIONS,
        partition=(3, 2, 1, 1),
        contact_orders=(3, 2),
        node_count=2,
        residual_separation=(
            discriminant(residual, S)
            * residual.subs(S, Rational(-2, 3))
            * residual.subs(S, Rational(-11, 9))
        ),
    )


@cache
def exact_three_double_sample_certificate() -> ExcessThreeSampleCertificate:
    """Return the exact three-contact-two representative."""

    cubic = 4 * S**3 + 9 * S**2 + 5 * S + 2
    return _sample_certificate(
        name="three contact-two points plus one node",
        parameters=THREE_DOUBLE_SAMPLE_PARAMETERS,
        factorization=THREE_DOUBLE_SAMPLE_FACTORIZATION,
        branch_equation=THREE_DOUBLE_BRANCH_EQUATION,
        y_scale=THREE_DOUBLE_Y_SCALE,
        relations=THREE_DOUBLE_RELATIONS,
        partition=(2, 2, 2, 1),
        contact_orders=(2, 2, 2),
        node_count=1,
        residual_separation=discriminant(cubic, S)
        * cubic.subs(S, Rational(-3, 2)),
    )


def main() -> None:
    """Print the three exact one-dimensional wall checkpoints."""

    algebra = exact_excess_three_algebra_certificate()
    samples = (
        exact_contact_four_sample_certificate(),
        exact_mixed_sample_certificate(),
        exact_three_double_sample_certificate(),
    )
    print("excess-three incidence algebra verified:", algebra.verified)
    print(
        "presentation image histograms:",
        {
            sample.partition: sample.complement_census.generated_order_histogram
            for sample in samples
        },
    )
    print("all one-dimensional strata excluded:", all(s.verified for s in samples))


if __name__ == "__main__":
    main()
