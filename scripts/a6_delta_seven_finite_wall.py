"""Exact finite collision endpoints of the conditional delta-seven family.

The three one-dimensional excess-three walls meet in finitely many valid
partition-changing points.  This module classifies those points and stores
the exact number-field van Kampen relations needed to replay their
single-three-cycle images into ``A6``.

Two principles are kept separate:

* algebraic Galois conjugacy is not, by itself, used to identify complement
  groups at distinct real embeddings;
* coefficientwise complex conjugation *is* an antiholomorphic complement
  homeomorphism.  It reverses meridian orientation, but the inverse of a
  single three-cycle is again a single three-cycle.

The algebraic classification covers ``(5,1,1)``, ``(4,2,1)``, ``(3,3,1)``,
and every absent or invalid coarser partition.  The stored topology covers
both ``(5,1,1)`` points, all four ``(4,2,1)`` points, and the conjugate pair
of nonreal ``(3,3,1)`` points.  The real ``(3,3,1)`` embedding is deliberately
left open until its separate van Kamp computation finishes.  As in the rest
of the repository, this is a conditional, computer-assisted statement, not
an unconditional proof of the plane Jacobian conjecture.
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
    discriminant,
    expand,
    gcd,
    groebner,
    rem,
    resultant,
    together,
)

from scripts.a6_delta_seven_deeper_wall import (
    CONTACT_FOUR_C_FACTOR,
    CONTACT_FOUR_FACTORIZATION,
    CONTACT_FOUR_L_FACTOR,
    CONTACT_FOUR_PARAMETERS,
    CONTACT_FOUR_RESIDUAL_FACTOR,
    CONTACT_FOUR_T_FACTOR,
    E,
    THREE_DOUBLE_C_FACTOR,
    THREE_DOUBLE_CUBIC,
    THREE_DOUBLE_FACTORIZATION,
    THREE_DOUBLE_T_FACTOR,
    V,
    _collision_coefficients,
    _parameters_from_collision,
)
from scripts.a6_delta_seven_discriminant_wall import (
    COEFFICIENT_SLICE_FIRST,
    COEFFICIENT_SLICE_SECOND,
)
from scripts.a6_delta_seven_generic import (
    ALPHA,
    CUSP_COLLISION_FACTOR,
    EXTRA_CRITICAL_FACTOR,
    S,
    TRIPLE_COLLISION_FACTOR,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    three_cycle_presentation_census,
)

CONTACT_FIVE_POLYNOMIAL: Final = 9 * V**2 + 51 * V + 34
Q5: Final = Symbol("Q5")
U5: Final = Symbol("U5")
CONTACT_FIVE_ANSATZ: Final = (
    (S - V) ** 5 * (S**2 + (5 * V + 6) * S + Q5)
)
CONTACT_FIVE_RESIDUAL_FACTOR: Final = (
    S**2 + (5 * V + 6) * S - (25 * V + 17) / 2
)

CONTACT_FOUR_TWO_POLYNOMIAL: Final = (
    3 * V**4 + 28 * V**3 + 80 * V**2 + 68 * V + 17
)
W42: Final = Symbol("W42")
U42: Final = Symbol("U42")
CONTACT_FOUR_TWO_ANSATZ: Final = (
    (S - V) ** 4
    * (S - W42) ** 2
    * (S + 6 + 4 * V + 2 * W42)
)
CONTACT_FOUR_TWO_DOUBLE_ROOT: Final = V * (3 * V**2 + 16 * V + 13) / 2
CONTACT_FOUR_TWO_SIMPLE_ROOT: Final = -(
    3 * V**3 + 16 * V**2 + 17 * V + 6
)
CONTACT_FOUR_TWO_RESIDUAL_FACTORIZATION: Final = (
    (S - CONTACT_FOUR_TWO_DOUBLE_ROOT) ** 2
    * (S - CONTACT_FOUR_TWO_SIMPLE_ROOT)
)

R33: Final = Symbol("R33")
Q33: Final = Symbol("Q33")
W33: Final = Symbol("W33")
CONTACT_THREE_THREE_POLYNOMIAL: Final = (
    9 * R33**3 - 93 * R33**2 - 413 * R33 - 391
)
CONTACT_THREE_THREE_ANSATZ: Final = (
    (S**2 - R33 * S + Q33) ** 3 * (S + 3 * R33 + 6)
)
CONTACT_THREE_THREE_E2: Final = -2 * (
    3 * R33**2 + 19 * R33 + 17
) / 27
CONTACT_THREE_THREE_QUADRATIC: Final = (
    S**2 - R33 * S + CONTACT_THREE_THREE_E2
)
CONTACT_THREE_THREE_FACTORIZATION: Final = (
    CONTACT_THREE_THREE_QUADRATIC**3 * (S + 3 * R33 + 6)
)
CONTACT_THREE_THREE_PARAMETERS: Final = _parameters_from_collision(
    CONTACT_THREE_THREE_FACTORIZATION
)

R322: Final = Symbol("R322")
Q322: Final = Symbol("Q322")
U322: Final = Symbol("U322")
CONTACT_THREE_TWO_TWO_SIGMA: Final = (-6 - 3 * R322) / 2
CONTACT_THREE_TWO_TWO_ANSATZ: Final = (
    (S - R322) ** 3
    * (S**2 - CONTACT_THREE_TWO_TWO_SIGMA * S + Q322) ** 2
)


# ---------------------------------------------------------------------------
# Exact endpoint presentations already reproduced over their number fields.
# ---------------------------------------------------------------------------

P5_NEGATIVE_RELATIONS: Final = (
    (2, 1, -2, -1),
    (-3, -2, 3, 1, -3, 2, 3, -1),
    (
        -3,
        -2,
        -3,
        -2,
        3,
        2,
        3,
        1,
        -3,
        -2,
        -3,
        2,
        3,
        2,
        3,
        1,
        -3,
        -2,
        -3,
        2,
        3,
        2,
        3,
        1,
        -3,
        -2,
        -3,
        -2,
        3,
        2,
        3,
        -1,
        -3,
        -2,
        -3,
        -2,
        3,
        2,
        3,
        -1,
    ),
    (-3, -2, -3, -2, 3, 2, 3, -1, -3, -2, -3, -2, 3, 2, 3, 2, 3, 1),
    (3, 2, 3, 2, 3, 2, 3, 2, 3, 2, -3, -2, -3, -2, -3, -2, -3, -2, -3, -2),
)

_P5_PLUS_A: Final = (
    2,
    3,
    -2,
    -3,
    -2,
    -1,
    2,
    3,
    -2,
    1,
    2,
    3,
    2,
    -3,
    -2,
    1,
)
_P5_PLUS_B: Final = (
    2,
    3,
    -2,
    -3,
    -2,
    -1,
    2,
    -3,
    -2,
    1,
    2,
    3,
    2,
    -3,
    -2,
    -1,
)
_P5_PLUS_C: Final = (2, 3, -2, -3, -2, -1, 2, -3, -2, 1, 2, 3, 2)
P5_POSITIVE_RELATIONS: Final = (
    (-3, -2, -3, -2, 1, 2, 3, 2, 3, -2, -3, -2, -1, 2, 3, 2),
    (
        -3,
        -2,
        -1,
        2,
        3,
        2,
        -3,
        -2,
        1,
        2,
        3,
        2,
        -3,
        -2,
        1,
        2,
        3,
        2,
        -3,
        -2,
        -1,
        2,
        3,
        -2,
        -3,
        -2,
        -1,
        2,
        3,
        -2,
    ),
    (-3, -2, -1) + _P5_PLUS_A * 5 + _P5_PLUS_B * 4 + _P5_PLUS_C,
    (-3, -2, -1, 2, -3, -2, 1, 2, 3, 2, -3, -2, -1, 2, 3, -2, 1, 2, 3, -2),
    (-3, -2, -1, 2, 1, 2),
)

P42_REAL_ZERO_RELATIONS: Final = (
    (-3, -2, -1, 2, 3, 2, -3, -2, 1, 2, 3, 2, -3, -2, 1, 2, 3, 2, -3, -2, -1, 2, 3, -2, -3, -2, -1, 2, 3, -2),
    (-3, -2, -1, 2, 3, -2, -3, -2, -1, 2, 3, -2, 1, 2, 3, 2, -3, -2, 1, 2, 3, -2, -3, -2, -1, 2, -3, -2, 1, 2, 3, 2),
    (-3, -2, -1, 2, -3, -2, 1, 2, 3, 2, -3, -2, -1, 2, 3, -2, 1, 2, 3, 2, -3, -2, -1, 2, 3, -2, 1, 2, 3, -2, -3, -2, -1, 2, -3, -2, 1, 2, 3, -2),
    (-3, -2, 1, 2, 3, -2, 1, 2, 3, -2, 1, 2, 3, -2, 1, 2, 3, -2, -1, 2, -3, -2, -1, 2, -3, -2, -1, 2, -3, -2, -1, 2),
    (-3, -2, -1, 2, 1, 2),
)

P42_REAL_ONE_RELATIONS: Final = (
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-3, 1, 3, 1, 3, 1, 3, 1, 3, -1, -3, -1, -3, -1, -3, -1),
    (-3, -1, -3, 2, 3, 1),
    (3, 2, 3, 2, -3, -2, -3, -2),
)

_P42_COMPLEX_A: Final = (2, 1, -2, -1, 3, 1, 2, -1, -2, -1, 3, 1)
_P42_COMPLEX_MIDDLE: Final = (
    2,
    1,
    -2,
    -1,
    3,
    1,
    2,
    -1,
    -1,
    -2,
    -1,
    -3,
    1,
    2,
    1,
    -2,
    -1,
    3,
    1,
    2,
    1,
)
_P42_COMPLEX_C: Final = (1, -2, -1, -3, 1, 2, -1, -2, -1, -3, 1, 2)
_P42_COMPLEX_TAIL: Final = (
    1,
    -2,
    -1,
    -3,
    1,
    2,
    -1,
    -2,
    3,
    1,
    -3,
    -3,
    -3,
    -2,
    3,
    3,
    3,
    -1,
    -3,
)
P42_COMPLEX_SIMPLIFIED_RELATIONS: Final = (
    (-3, 1, 3, 1, 3, -1, -3, -1),
    (-3, -2, 3, 1, -3, 2, 3, -1),
    (-3, -2, -1, -3, 2, -1, -2, -1, -3, 1, 2, -1, 3, 1, 2, 1, 3, 1),
    (-2, -1, -3) + (1, 2, 1, -2, -1, 3) * 2 + (1, 2, 1) + (-2, -1, -3, 1, 2, -1) * 2,
    _P42_COMPLEX_A * 2
    + _P42_COMPLEX_MIDDLE
    + _P42_COMPLEX_C * 2
    + _P42_COMPLEX_TAIL,
)

P33_COMPLEX_SIMPLIFIED_RELATIONS: Final = (
    (3, 2, 3, -1, -1, -2, 1, 1, -3, -2),
    (-2, 1, 1, 2, 3, -1, 2, -3, -2, -1),
    (-3, -2, -1, 2, -3, -2, 1, 2, 3, 2),
    (-3, -3, -2, -1, 2, 2, 3, 2, -3, -2, -2, 1, 2, 3, 3, -2),
    (
        -1,
        2,
        -1,
        2,
        -3,
        1,
        2,
        2,
        3,
        1,
        2,
        3,
        -2,
        1,
        -2,
        -1,
        -2,
        -3,
        -2,
        -1,
    ),
)


@dataclass(frozen=True, slots=True)
class FiniteWallAlgebraCertificate:
    """Exact finite-partition and absence identities."""

    contact_five_factor_remainder: Expr
    contact_five_elimination_identity: Expr
    contact_five_saturated_basis: tuple[Expr, ...]
    contact_five_discarded_cusp_identity: Expr
    contact_five_discriminant: Expr
    contact_five_residual_discriminant_resultant: Expr
    contact_five_residual_value_resultant: Expr
    contact_five_invalid_resultants: tuple[Expr, ...]
    contact_four_chart_denominator_resultants: tuple[Expr, Expr]
    contact_four_two_discriminant: Expr
    contact_four_two_elimination_identity: Expr
    contact_four_two_saturated_basis: tuple[Expr, ...]
    contact_four_two_discarded_branch_gcds: tuple[Expr, Expr]
    contact_four_two_discarded_wall_identities: tuple[Expr, Expr, Expr]
    contact_four_two_factor_remainder: Expr
    contact_four_two_distinct_root_resultant: Expr
    contact_four_two_contact_separation_resultant: Expr
    contact_four_two_invalid_resultants: tuple[Expr, ...]
    contact_three_three_first_slice_remainder: Expr
    contact_three_three_second_slice_remainder: Expr
    contact_three_three_elimination_identity: Expr
    contact_three_three_saturated_basis: tuple[Expr, ...]
    contact_three_three_discarded_branch_gcds: tuple[Expr, Expr]
    contact_three_three_degenerate_cusp_values: tuple[Expr, Expr]
    contact_three_three_discriminant: Expr
    contact_three_three_invalid_resultants: tuple[Expr, ...]
    contact_three_two_two_elimination_identity: Expr
    contact_three_two_two_saturated_basis: tuple[Expr, ...]
    contact_three_two_two_branch_gcds: tuple[Expr, Expr, Expr, Expr]
    contact_three_two_two_branch_wall_values: tuple[Expr, Expr, Expr, Expr]
    contact_three_two_two_boundary_identity: Expr
    contact_three_two_two_t_value: Expr
    contact_three_two_two_c_value: Expr
    six_one_gcd: Expr
    five_two_gcd: Expr
    four_three_gcd: Expr
    two_root_cusp_values: tuple[Expr, Expr, Expr]
    seven_slice_values: tuple[Expr, Expr]

    @property
    def verified(self) -> bool:
        """Whether the finite partition list and invalid boundaries agree."""

        return bool(
            self.contact_five_factor_remainder == 0
            and self.contact_five_elimination_identity == 0
            and self.contact_five_saturated_basis
            == (
                8 * U5 - 9 * V - 42,
                2 * Q5 + 25 * V + 17,
                CONTACT_FIVE_POLYNOMIAL,
            )
            and self.contact_five_discarded_cusp_identity == 0
            and self.contact_five_discriminant == 1377
            and self.contact_five_residual_discriminant_resultant != 0
            and self.contact_five_residual_value_resultant != 0
            and all(value != 0 for value in self.contact_five_invalid_resultants)
            and all(
                value != 0
                for value in self.contact_four_chart_denominator_resultants
            )
            and self.contact_four_two_discriminant == -108800
            and self.contact_four_two_elimination_identity == 0
            and self.contact_four_two_saturated_basis
            == (
                120 * U42
                + 3931 * V**3
                + 33809 * V**2
                + 80057 * V
                + 30467,
                -3 * V**3 - 16 * V**2 - 13 * V + 2 * W42,
                CONTACT_FOUR_TWO_POLYNOMIAL,
            )
            and self.contact_four_two_discarded_branch_gcds
            == (W42 + Rational(4, 3), W42 + 1)
            and self.contact_four_two_discarded_wall_identities == (0, 0, 0)
            and self.contact_four_two_factor_remainder == 0
            and self.contact_four_two_distinct_root_resultant != 0
            and self.contact_four_two_contact_separation_resultant != 0
            and all(value != 0 for value in self.contact_four_two_invalid_resultants)
            and self.contact_three_three_first_slice_remainder == 0
            and self.contact_three_three_second_slice_remainder == 0
            and self.contact_three_three_elimination_identity == 0
            and self.contact_three_three_saturated_basis
            == (
                -171 * R33**2 + 2028 * R33 + 72 * W33 + 4763,
                6 * R33**2 + 38 * R33 + 27 * Q33 + 34,
                CONTACT_THREE_THREE_POLYNOMIAL,
            )
            and self.contact_three_three_discarded_branch_gcds
            == (
                Q33**2 - 2 * Q33 + 1,
                Q33 - Rational(2, 3),
            )
            and self.contact_three_three_degenerate_cusp_values == (0, 0)
            and self.contact_three_three_discriminant == -13996800
            and all(value != 0 for value in self.contact_three_three_invalid_resultants)
            and self.contact_three_two_two_elimination_identity == 0
            and self.contact_three_two_two_saturated_basis == (1,)
            and self.contact_three_two_two_branch_gcds
            == (
                4 * Q322**2 - 4 * Q322 + 1,
                4 * (Q322 - 1),
                Q322 - 1,
                Q322 - Rational(7, 5),
            )
            and self.contact_three_two_two_branch_wall_values == (0, 0, 0, 0)
            and self.contact_three_two_two_boundary_identity == 0
            and self.contact_three_two_two_t_value == 0
            and self.contact_three_two_two_c_value == 0
            and expand(self.six_one_gcd - (V + 1) ** 3) == 0
            and expand(self.five_two_gcd - 5 * (V + 1) ** 2) == 0
            and expand(self.four_three_gcd - (V + 1)) == 0
            and all(value == 0 for value in self.two_root_cusp_values)
            and all(value != 0 for value in self.seven_slice_values)
        )


def _numerator(expression: Expr) -> Expr:
    """Return the expanded numerator of a rational expression."""

    return expand(together(cancel(expression)).as_numer_denom()[0])


def _slice_values(polynomial: Expr) -> tuple[Expr, Expr]:
    """Evaluate the two coefficient-slice equations on a septic."""

    coefficients = _collision_coefficients(polynomial)
    return (
        cancel(COEFFICIENT_SLICE_FIRST.subs(coefficients)),
        cancel(COEFFICIENT_SLICE_SECOND.subs(coefficients)),
    )


def _partition_slice_gcd(first_multiplicity: int, second_multiplicity: int) -> Expr:
    """Gcd of the two slice equations for a two-root partition."""

    polynomial = _two_root_partition_polynomial(
        first_multiplicity,
        second_multiplicity,
    )
    first, second = _slice_values(polynomial)
    return gcd(_numerator(first), _numerator(second))


def _two_root_partition_polynomial(
    first_multiplicity: int,
    second_multiplicity: int,
) -> Expr:
    """Monic slice septic with the requested two-root multiplicities."""

    second_root = (-6 - first_multiplicity * V) / second_multiplicity
    return (S - V) ** first_multiplicity * (
        S - second_root
    ) ** second_multiplicity


def _two_root_cusp_value(
    first_multiplicity: int,
    second_multiplicity: int,
) -> Expr:
    """Evaluate the cusp-validity factor at the slice gcd root ``v=-1``."""

    parameters = _parameters_from_collision(
        _two_root_partition_polynomial(
            first_multiplicity,
            second_multiplicity,
        )
    )
    return cancel(CUSP_COLLISION_FACTOR.subs(parameters).subs(V, -1))


@cache
def exact_finite_wall_algebra_certificate() -> FiniteWallAlgebraCertificate:
    """Build the exact finite endpoint and absence certificate."""

    p5_general_first, p5_general_second = _slice_values(CONTACT_FIVE_ANSATZ)
    p5_saturated_basis = groebner(
        (
            p5_general_first,
            p5_general_second,
            1 - U5 * (V + 1),
        ),
        U5,
        Q5,
        V,
        order="lex",
    )
    contact_five_difference = cancel(
        CONTACT_FOUR_FACTORIZATION
        - (S - V) ** 5 * CONTACT_FIVE_RESIDUAL_FACTOR
    )
    contact_five_numerator = _numerator(contact_five_difference)

    p5_invalid = tuple(
        resultant(CONTACT_FIVE_POLYNOMIAL, _numerator(expression), V)
        for expression in (
            CONTACT_FOUR_PARAMETERS[ALPHA],
            CONTACT_FOUR_C_FACTOR,
            CONTACT_FOUR_L_FACTOR,
            CONTACT_FOUR_T_FACTOR,
        )
    )
    p42_invalid = tuple(
        resultant(CONTACT_FOUR_TWO_POLYNOMIAL, _numerator(expression), V)
        for expression in (
            CONTACT_FOUR_PARAMETERS[ALPHA],
            CONTACT_FOUR_C_FACTOR,
            CONTACT_FOUR_L_FACTOR,
            CONTACT_FOUR_T_FACTOR,
        )
    )
    p42_general_first, p42_general_second = _slice_values(
        CONTACT_FOUR_TWO_ANSATZ
    )
    p42_saturated_basis = groebner(
        (
            p42_general_first,
            p42_general_second,
            1 - U42 * (V + 1) * (3 * V + 1) * (15 * V + 11),
        ),
        U42,
        W42,
        V,
        order="lex",
    )

    p33_first, p33_second = _slice_values(CONTACT_THREE_THREE_FACTORIZATION)
    p33_general_first, p33_general_second = _slice_values(
        CONTACT_THREE_THREE_ANSATZ
    )
    p33_saturated_basis = groebner(
        (
            p33_general_first,
            p33_general_second,
            1 - W33 * (R33 + 2) * (3 * R33 + 5),
        ),
        W33,
        Q33,
        R33,
        order="lex",
    )
    p33_invalid = tuple(
        resultant(
            CONTACT_THREE_THREE_POLYNOMIAL,
            _numerator(expression),
            R33,
        )
        for expression in (
            CONTACT_THREE_THREE_PARAMETERS[ALPHA],
            CUSP_COLLISION_FACTOR.subs(CONTACT_THREE_THREE_PARAMETERS),
            EXTRA_CRITICAL_FACTOR.subs(CONTACT_THREE_THREE_PARAMETERS),
            TRIPLE_COLLISION_FACTOR.subs(CONTACT_THREE_THREE_PARAMETERS),
            discriminant(CONTACT_THREE_THREE_QUADRATIC, S),
            CONTACT_THREE_THREE_QUADRATIC.subs(S, -3 * R33 - 6),
        )
    )

    p322_general_first, p322_general_second = _slice_values(
        CONTACT_THREE_TWO_TWO_ANSATZ
    )
    p322_first_numerator = _numerator(p322_general_first)
    p322_second_numerator = _numerator(p322_general_second)
    p322_discarded_factor = (
        (R322 + 1)
        * (R322 + 2)
        * (3 * R322 + 2)
        * (5 * R322 + 2)
    )
    p322_saturated_basis = groebner(
        (
            p322_first_numerator,
            p322_second_numerator,
            1 - U322 * p322_discarded_factor,
        ),
        U322,
        Q322,
        R322,
        order="lex",
    )

    seven_polynomial = (S + Rational(6, 7)) ** 7
    return FiniteWallAlgebraCertificate(
        contact_five_factor_remainder=rem(
            contact_five_numerator,
            CONTACT_FIVE_POLYNOMIAL,
            V,
        ),
        contact_five_elimination_identity=expand(
            resultant(p5_general_first, p5_general_second, Q5)
            + 5 * (V + 1) ** 7 * CONTACT_FIVE_POLYNOMIAL
        ),
        contact_five_saturated_basis=tuple(
            expand(polynomial.as_expr())
            for polynomial in p5_saturated_basis.polys
        ),
        contact_five_discarded_cusp_identity=cancel(
            CUSP_COLLISION_FACTOR.subs(
                _parameters_from_collision(CONTACT_FIVE_ANSATZ.subs(V, -1))
            )
        ),
        contact_five_discriminant=discriminant(CONTACT_FIVE_POLYNOMIAL, V),
        contact_five_residual_discriminant_resultant=resultant(
            CONTACT_FIVE_POLYNOMIAL,
            _numerator(discriminant(CONTACT_FIVE_RESIDUAL_FACTOR, S)),
            V,
        ),
        contact_five_residual_value_resultant=resultant(
            CONTACT_FIVE_POLYNOMIAL,
            _numerator(CONTACT_FIVE_RESIDUAL_FACTOR.subs(S, V)),
            V,
        ),
        contact_five_invalid_resultants=p5_invalid,
        contact_four_chart_denominator_resultants=(
            resultant(CONTACT_FIVE_POLYNOMIAL, 3 * V + 13, V),
            resultant(CONTACT_FOUR_TWO_POLYNOMIAL, 3 * V + 13, V),
        ),
        contact_four_two_discriminant=discriminant(
            CONTACT_FOUR_TWO_POLYNOMIAL,
            V,
        ),
        contact_four_two_elimination_identity=expand(
            resultant(p42_general_first, p42_general_second, W42)
            + 40
            * (V + 1) ** 17
            * (3 * V + 1) ** 2
            * (15 * V + 11)
            * CONTACT_FOUR_TWO_POLYNOMIAL
        ),
        contact_four_two_saturated_basis=tuple(
            expand(polynomial.as_expr())
            for polynomial in p42_saturated_basis.polys
        ),
        contact_four_two_discarded_branch_gcds=(
            gcd(
                p42_general_first.subs(V, Rational(-1, 3)),
                p42_general_second.subs(V, Rational(-1, 3)),
                W42,
            ),
            gcd(
                p42_general_first.subs(V, Rational(-11, 15)),
                p42_general_second.subs(V, Rational(-11, 15)),
                W42,
            ),
        ),
        contact_four_two_discarded_wall_identities=(
            cancel(
                CUSP_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_FOUR_TWO_ANSATZ.subs(V, -1)
                    )
                )
            ),
            cancel(
                TRIPLE_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_FOUR_TWO_ANSATZ.subs(
                            {V: Rational(-1, 3), W42: Rational(-4, 3)}
                        )
                    )
                )
            ),
            cancel(
                CUSP_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_FOUR_TWO_ANSATZ.subs(
                            {V: Rational(-11, 15), W42: -1}
                        )
                    )
                )
            ),
        ),
        contact_four_two_factor_remainder=rem(
            _numerator(
                CONTACT_FOUR_RESIDUAL_FACTOR
                - CONTACT_FOUR_TWO_RESIDUAL_FACTORIZATION
            ),
            CONTACT_FOUR_TWO_POLYNOMIAL,
            V,
        ),
        contact_four_two_distinct_root_resultant=resultant(
            CONTACT_FOUR_TWO_POLYNOMIAL,
            _numerator(
                CONTACT_FOUR_TWO_DOUBLE_ROOT
                - CONTACT_FOUR_TWO_SIMPLE_ROOT
            ),
            V,
        ),
        contact_four_two_contact_separation_resultant=resultant(
            CONTACT_FOUR_TWO_POLYNOMIAL,
            _numerator(CONTACT_FOUR_RESIDUAL_FACTOR.subs(S, V)),
            V,
        ),
        contact_four_two_invalid_resultants=p42_invalid,
        contact_three_three_first_slice_remainder=rem(
            _numerator(p33_first),
            CONTACT_THREE_THREE_POLYNOMIAL,
            R33,
        ),
        contact_three_three_second_slice_remainder=rem(
            _numerator(p33_second),
            CONTACT_THREE_THREE_POLYNOMIAL,
            R33,
        ),
        contact_three_three_elimination_identity=expand(
            resultant(
                p33_general_first,
                p33_general_second,
                Q33,
            )
            + (R33 + 2) ** 8
            * (3 * R33 + 5) ** 3
            * CONTACT_THREE_THREE_POLYNOMIAL
        ),
        contact_three_three_saturated_basis=tuple(
            expand(polynomial.as_expr())
            for polynomial in p33_saturated_basis.polys
        ),
        contact_three_three_discarded_branch_gcds=(
            gcd(
                p33_general_first.subs(R33, -2),
                p33_general_second.subs(R33, -2),
                Q33,
            ),
            gcd(
                p33_general_first.subs(R33, Rational(-5, 3)),
                p33_general_second.subs(R33, Rational(-5, 3)),
                Q33,
            ),
        ),
        contact_three_three_degenerate_cusp_values=(
            cancel(
                CUSP_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_THREE_THREE_ANSATZ.subs(
                            {R33: -2, Q33: 1}
                        )
                    )
                )
            ),
            cancel(
                CUSP_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_THREE_THREE_ANSATZ.subs(
                            {R33: Rational(-5, 3), Q33: Rational(2, 3)}
                        )
                    )
                )
            ),
        ),
        contact_three_three_discriminant=discriminant(
            CONTACT_THREE_THREE_POLYNOMIAL,
            R33,
        ),
        contact_three_three_invalid_resultants=p33_invalid,
        contact_three_two_two_elimination_identity=expand(
            resultant(
                p322_first_numerator,
                p322_second_numerator,
                Q322,
            )
            - 80
            * (R322 + 1) ** 8
            * (R322 + 2) ** 3
            * (3 * R322 + 2) ** 2
            * (5 * R322 + 2)
        ),
        contact_three_two_two_saturated_basis=tuple(
            expand(polynomial.as_expr())
            for polynomial in p322_saturated_basis.polys
        ),
        contact_three_two_two_branch_gcds=(
            gcd(
                p322_first_numerator.subs(R322, -1),
                p322_second_numerator.subs(R322, -1),
                Q322,
            ),
            gcd(
                p322_first_numerator.subs(R322, -2),
                p322_second_numerator.subs(R322, -2),
                Q322,
            ),
            gcd(
                p322_first_numerator.subs(R322, Rational(-2, 3)),
                p322_second_numerator.subs(R322, Rational(-2, 3)),
                Q322,
            ),
            gcd(
                p322_first_numerator.subs(R322, Rational(-2, 5)),
                p322_second_numerator.subs(R322, Rational(-2, 5)),
                Q322,
            ),
        ),
        contact_three_two_two_branch_wall_values=(
            cancel(
                CUSP_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_THREE_TWO_TWO_ANSATZ.subs(
                            {R322: -1, Q322: Rational(1, 2)}
                        )
                    )
                )
            ),
            cancel(
                TRIPLE_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_THREE_TWO_TWO_ANSATZ.subs(
                            {R322: -2, Q322: 1}
                        )
                    )
                )
            ),
            cancel(
                CUSP_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_THREE_TWO_TWO_ANSATZ.subs(
                            {R322: Rational(-2, 3), Q322: 1}
                        )
                    )
                )
            ),
            cancel(
                CUSP_COLLISION_FACTOR.subs(
                    _parameters_from_collision(
                        CONTACT_THREE_TWO_TWO_ANSATZ.subs(
                            {R322: Rational(-2, 5), Q322: Rational(7, 5)}
                        )
                    )
                )
            ),
        ),
        contact_three_two_two_boundary_identity=expand(
            THREE_DOUBLE_CUBIC.subs(S, -6 - 2 * E)
            + 5 * (E + 2) * (2 * E + 5) ** 2
        ),
        contact_three_two_two_t_value=expand(
            THREE_DOUBLE_T_FACTOR.subs(E, -2)
        ),
        contact_three_two_two_c_value=expand(
            THREE_DOUBLE_C_FACTOR.subs(E, Rational(-5, 2))
        ),
        six_one_gcd=_partition_slice_gcd(6, 1),
        five_two_gcd=_partition_slice_gcd(5, 2),
        four_three_gcd=_partition_slice_gcd(4, 3),
        two_root_cusp_values=(
            _two_root_cusp_value(6, 1),
            _two_root_cusp_value(5, 2),
            _two_root_cusp_value(4, 3),
        ),
        seven_slice_values=_slice_values(seven_polynomial),
    )


@dataclass(frozen=True, slots=True)
class EndpointPresentationCertificate:
    """One directly checked number-field embedding of a finite endpoint."""

    name: str
    partition: tuple[int, ...]
    relations: tuple[tuple[int, ...], ...]
    complement_census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether the embedding has only cyclic three-cycle images."""

        return bool(
            self.complement_census.assignments == 40**3
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def checked_endpoint_presentations() -> tuple[EndpointPresentationCertificate, ...]:
    """Return checked embeddings; the real P33 endpoint is not yet stored."""

    data = (
        ("P5 negative real embedding", (5, 1, 1), P5_NEGATIVE_RELATIONS),
        ("P5 positive real embedding", (5, 1, 1), P5_POSITIVE_RELATIONS),
        ("P42 first real embedding", (4, 2, 1), P42_REAL_ZERO_RELATIONS),
        ("P42 second real embedding", (4, 2, 1), P42_REAL_ONE_RELATIONS),
        (
            "P42 nonreal embedding (conjugate covers the fourth)",
            (4, 2, 1),
            P42_COMPLEX_SIMPLIFIED_RELATIONS,
        ),
        (
            "P33 nonreal embedding (conjugate covers the third)",
            (3, 3, 1),
            P33_COMPLEX_SIMPLIFIED_RELATIONS,
        ),
    )
    return tuple(
        EndpointPresentationCertificate(
            name=name,
            partition=partition,
            relations=relations,
            complement_census=three_cycle_presentation_census(relations),
        )
        for name, partition, relations in data
    )


def main() -> None:
    """Print the exact finite-wall checkpoint currently stored."""

    algebra = exact_finite_wall_algebra_certificate()
    presentations = checked_endpoint_presentations()
    print("finite partition algebra verified:", algebra.verified)
    print(
        "checked endpoint histograms:",
        {
            endpoint.name: endpoint.complement_census.generated_order_histogram
            for endpoint in presentations
        },
    )
    print("all stored endpoints excluded:", all(p.verified for p in presentations))


if __name__ == "__main__":
    main()
