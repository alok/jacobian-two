"""Exact audit of the first conditional delta-seven discriminant strata.

For the degree-``(3, 10)`` normal form, a collision root of multiplicity
``m`` is a contact-``m`` pair on the valid locus.  The generic repeated-root
divisor and the two generic components of its singular locus therefore have
collision partitions

``(2,1,1,1,1,1)``, ``(3,1,1,1,1)``, and ``(2,2,1,1,1)``.

This module checks exact incidence charts, rank-drop boundaries, representative
curve geometry, Sage van Kampen presentations, and all single-three-cycle
images into ``A6``.  Each representative complement is infinite cyclic.  A
separate hostile census proves that the local contact relations themselves do
*not* imply that conclusion, so propagation across each connected
equisingular stratum remains an explicit Whitney--Thom theorem dependency.

The result is conditional on the one-pair and finite-singularity hypotheses
of ``a6_one_pair_infinity.py``.  It is not a proof of the plane Jacobian
conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from itertools import combinations, product
from typing import Final

from sympy import (
    Expr,
    Matrix,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    factor,
    gcd,
    groebner,
    resultant,
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
    T,
    TANGENCY_POLYNOMIAL,
    TRIPLE_COLLISION_FACTOR,
    X,
    Y,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    alternating_group_six,
    permutation_power,
    three_cycle_presentation_census,
)
from scripts.six_sheet_monodromy import (
    A6_COLLISION_MERIDIAN,
    A6_TORUS_2_5_LOCAL_GENERATORS,
    Permutation,
    compose,
    cycle_type,
    generated_group,
)

R: Final = Symbol("r")
QROOT: Final = Symbol("q")
FREE_GAMMA: Final = Symbol("g")
FREE_DELTA: Final = Symbol("d")
ROOT_SUM: Final = Symbol("a")
ROOT_PRODUCT: Final = Symbol("b")

H0: Final = Symbol("h0")
H1: Final = Symbol("h1")
H2: Final = Symbol("h2")
H3: Final = Symbol("h3")
H4: Final = Symbol("h4")
H5: Final = Symbol("h5")

COEFFICIENT_SLICE_FIRST: Final = H3 - 4 * H4 + 10 * H5 - 85
COEFFICIENT_SLICE_SECOND: Final = (
    H0 - 2 * H1 + 3 * H2 - 11 * H4 + 34 * H5 - 306
)
COEFFICIENT_SLICE_INVERSE: Final = {
    ALPHA: -H0,
    BETA: H1 + H0,
    DELTA: (H4 + H1 + H0 + 18 - 2 * H5) / 3,
    GAMMA: H4 + H1 + H0 + 27 - 3 * H5,
}

DOUBLE_INCIDENCE_DENOMINATOR: Final = (
    (R - 1) * (2 * R**3 + 7 * R**2 + 6 * R - 1)
)
DOUBLE_INCIDENCE_PARAMETERS: Final = {
    ALPHA: R**2
    / DOUBLE_INCIDENCE_DENOMINATOR
    * (
        FREE_DELTA
        * (3 * R**5 + 21 * R**4 + 36 * R**3 + 8 * R**2 - 13 * R + 1)
        - FREE_GAMMA * (R**5 + 7 * R**4 + 12 * R**3 + 4 * R**2 + 4)
        + 3 * R**7
        + 25 * R**6
        + 71 * R**5
        + 67 * R**4
        - 11 * R**3
        - 25 * R**2
        + 10 * R
    ),
    BETA: R
    / DOUBLE_INCIDENCE_DENOMINATOR
    * (
        FREE_DELTA * (9 * R**4 + 21 * R**3 - 3 * R**2 - 21 * R + 2)
        + FREE_GAMMA * (-3 * R**4 - 5 * R**3 + 6 * R**2 + 6 * R - 8)
        + 5 * R**6
        + 25 * R**5
        + 25 * R**4
        - 25 * R**3
        - 25 * R**2
        + 15 * R
    ),
    GAMMA: FREE_GAMMA,
    DELTA: FREE_DELTA,
}

TRIPLE_INCIDENCE_DENOMINATOR: Final = (
    3 * R**5 + 9 * R**4 - 4 * R**3 - 24 * R**2 + 6 * R - 4
)
TRIPLE_INCIDENCE_PARAMETERS: Final = {
    ALPHA: -R**3
    / TRIPLE_INCIDENCE_DENOMINATOR
    * (
        FREE_DELTA * (6 * R**2 + 24 * R + 26)
        + 3 * R**7
        + 29 * R**6
        + 81 * R**5
        + 81 * R**4
        + 36 * R**3
        + 66 * R**2
        + 60 * R
        - 20
    ),
    BETA: -R**2
    / TRIPLE_INCIDENCE_DENOMINATOR
    * (
        FREE_DELTA * (-9 * R**3 - 27 * R**2 + 12 * R + 78)
        + 15 * R**6
        + 45 * R**5
        - 25 * R**4
        - 135 * R**3
        + 60 * R**2
        + 200 * R
        - 60
    ),
    GAMMA: 1
    / TRIPLE_INCIDENCE_DENOMINATOR
    * (
        FREE_DELTA
        * (9 * R**5 + 27 * R**4 - 12 * R**3 - 72 * R**2 + 21 * R - 1)
        + 15 * R**7
        + 70 * R**6
        + 60 * R**5
        - 135 * R**4
        - 155 * R**3
        + 90 * R**2
        - 15 * R
    ),
    DELTA: FREE_DELTA,
}

TWO_DOUBLE_DENOMINATOR: Final = (
    3 * ROOT_SUM**3
    + ROOT_SUM**2 * ROOT_PRODUCT
    + 21 * ROOT_SUM**2
    + 10 * ROOT_SUM * ROOT_PRODUCT
    + 42 * ROOT_SUM
    + 2 * ROOT_PRODUCT**2
    + 12 * ROOT_PRODUCT
    + 26
)
TWO_DOUBLE_N0: Final = (
    11 * ROOT_SUM**4
    + 13 * ROOT_SUM**3 * ROOT_PRODUCT
    + 73 * ROOT_SUM**3
    + 4 * ROOT_SUM**2 * ROOT_PRODUCT**2
    + 65 * ROOT_SUM**2 * ROOT_PRODUCT
    + 183 * ROOT_SUM**2
    + 12 * ROOT_SUM * ROOT_PRODUCT**2
    + 112 * ROOT_SUM * ROOT_PRODUCT
    + 204 * ROOT_SUM
    - ROOT_PRODUCT**3
    + 11 * ROOT_PRODUCT**2
    + 65 * ROOT_PRODUCT
    + 85
)
TWO_DOUBLE_N1: Final = (
    24 * ROOT_SUM**4
    + 20 * ROOT_SUM**3 * ROOT_PRODUCT
    + 182 * ROOT_SUM**3
    + 4 * ROOT_SUM**2 * ROOT_PRODUCT**2
    + 122 * ROOT_SUM**2 * ROOT_PRODUCT
    + 514 * ROOT_SUM**2
    + 17 * ROOT_SUM * ROOT_PRODUCT**2
    + 252 * ROOT_SUM * ROOT_PRODUCT
    + 635 * ROOT_SUM
    - ROOT_PRODUCT**3
    + 23 * ROOT_PRODUCT**2
    + 169 * ROOT_PRODUCT
    + 289
)
TWO_DOUBLE_C2: Final = 6 + 2 * ROOT_SUM
TWO_DOUBLE_C1: Final = TWO_DOUBLE_N1 / TWO_DOUBLE_DENOMINATOR
TWO_DOUBLE_C0: Final = 2 * TWO_DOUBLE_N0 / TWO_DOUBLE_DENOMINATOR
TWO_DOUBLE_E5: Final = (
    TWO_DOUBLE_C1
    - 2 * ROOT_SUM * TWO_DOUBLE_C2
    + ROOT_SUM**2
    + 2 * ROOT_PRODUCT
)
TWO_DOUBLE_E4: Final = (
    TWO_DOUBLE_C0
    - 2 * ROOT_SUM * TWO_DOUBLE_C1
    + (ROOT_SUM**2 + 2 * ROOT_PRODUCT) * TWO_DOUBLE_C2
    - 2 * ROOT_SUM * ROOT_PRODUCT
)
TWO_DOUBLE_PARAMETERS: Final = {
    ALPHA: -ROOT_PRODUCT**2 * TWO_DOUBLE_C0,
    BETA: (
        ROOT_PRODUCT**2 * TWO_DOUBLE_C1
        + ROOT_PRODUCT
        * (ROOT_PRODUCT - 2 * ROOT_SUM)
        * TWO_DOUBLE_C0
    ),
}
TWO_DOUBLE_PARAMETERS[DELTA] = (
    TWO_DOUBLE_E4
    + TWO_DOUBLE_PARAMETERS[BETA]
    - 2 * TWO_DOUBLE_E5
    + 18
) / 3
TWO_DOUBLE_PARAMETERS[GAMMA] = (
    TWO_DOUBLE_E4 + TWO_DOUBLE_PARAMETERS[BETA] - 3 * TWO_DOUBLE_E5 + 27
)
TWO_DOUBLE_REPEATED_FACTOR: Final = S**2 - ROOT_SUM * S + ROOT_PRODUCT
TWO_DOUBLE_RESIDUAL_FACTOR: Final = (
    S**3 + TWO_DOUBLE_C2 * S**2 + TWO_DOUBLE_C1 * S + TWO_DOUBLE_C0
)

# The rational two-double chart divides by ``TWO_DOUBLE_DENOMINATOR``.  To
# prove that this loses no valid component, retain the two residual cubic
# coefficients before solving the coefficient-slice equations.  Their
# compatibility ideal on the denominator boundary has only two geometric
# possibilities.  The first has a repeated root at ``s = -1`` and hence lies
# on the forbidden cusp boundary.  The second is the line below; its extra
# critical factor vanishes identically.
TWO_DOUBLE_FREE_C0: Final = Symbol("c0")
TWO_DOUBLE_FREE_C1: Final = Symbol("c1")
TWO_DOUBLE_GENERAL_RESIDUAL_FACTOR: Final = (
    S**3
    + (6 + 2 * ROOT_SUM) * S**2
    + TWO_DOUBLE_FREE_C1 * S
    + TWO_DOUBLE_FREE_C0
)
TWO_DOUBLE_GENERAL_POLYNOMIAL: Final = (
    TWO_DOUBLE_REPEATED_FACTOR**2 * TWO_DOUBLE_GENERAL_RESIDUAL_FACTOR
)
TWO_DOUBLE_EXCEPTIONAL_SUM: Final = Rational(-5, 3)
TWO_DOUBLE_EXCEPTIONAL_PRODUCT: Final = Rational(4, 9)
TWO_DOUBLE_EXCEPTIONAL_C1: Final = 2 * TWO_DOUBLE_FREE_C0 + 1
TWO_DOUBLE_EXCEPTIONAL_FACTORIZATION: Final = (
    (3 * S + 1) ** 2
    * (3 * S + 4) ** 2
    * (
        6 * TWO_DOUBLE_FREE_C0 * S
        + 3 * TWO_DOUBLE_FREE_C0
        + 3 * S**3
        + 8 * S**2
        + 3 * S
    )
    / 243
)

# Exact Sage 10.8 affine presentations.  Generator indices 1..3 are geometric
# meridians of a generic vertical fiber.
GENERIC_G_RELATIONS: Final = (
    (-3, 1, 3, -1),
    (3, 2, -3, -2),
    (2, 1, -2, -1),
    (-2, -1, -3, 1, 3, 1),
    (-3, 1, 3, 1, 3, 1, -3, -1, -3, -1),
    (-3, 1, 3, -1),
    (-3, -2, 3, 1, -3, 2, 3, -1),
    (3, 2, 3, 2, -3, -2, -3, -2),
)

CONTACT_THREE_RELATIONS: Final = (
    (2, 1, -2, -1),
    (-3, 1, 3, -1),
    (-3, -2, 3, 1, -3, 2, 3, -1),
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
    (-3, -2, 3, -1, -3, -2, 3, 2, 3, 1),
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
        -1,
    ),
    (3, 2, 3, 2, 3, 2, -3, -2, -3, -2, -3, -2),
)

TWO_DOUBLE_RELATIONS: Final = (
    (2, 1, 2, 1, -2, -1, -2, -1),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-3, 1, 3, 1, 3, -1, -3, -1),
    (-3, -2, 1, 2),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-3, 1, 3, -1),
    (3, 2, -3, -2),
)

GENERIC_G_PARAMETERS: Final = {ALPHA: -27, BETA: 0, GAMMA: 0, DELTA: 0}
GENERIC_G_COLLISION_FACTORIZATION: Final = (
    (S + 3) ** 2
    * (S**2 + S - 1)
    * (S**3 - S**2 + 2 * S - 3)
)
GENERIC_G_BRANCH_EQUATION: Final = (
    -X**10
    + 135 * X**8
    + 15 * X**6 * Y
    - 3645 * X**6
    + 83 * X**5 * Y
    + 20412 * X**5
    - 675 * X**4 * Y
    + 10 * X**3 * Y**2
    + 3780 * X**3 * Y
    - 25 * X**2 * Y**2
    + 145 * X * Y**2
    + Y**3
    - 28 * Y**2
)

CONTACT_THREE_PARAMETERS: Final = {
    ALPHA: -4,
    BETA: 2,
    GAMMA: Rational(9, 4),
    DELTA: 2,
}
CONTACT_THREE_COLLISION_FACTORIZATION: Final = (
    (S + 2) ** 3 * (4 * S**4 + 3 * S**2 - 4 * S + 2) / 4
)
CONTACT_THREE_BRANCH_EQUATION: Final = (
    64 * X**10
    + 64 * X**9
    + 501 * X**8
    - 1200 * X**7
    + 112 * X**6 * Y
    + 3072 * X**6
    - 1216 * X**5 * Y
    - 3328 * X**5
    + 1656 * X**4 * Y
    - 256 * X**3 * Y**2
    - 4160 * X**3 * Y
    + 128 * X**2 * Y**2
    - 1024 * X * Y**2
    - 64 * Y**3
    + 208 * Y**2
)

TWO_DOUBLE_SAMPLE_PARAMETERS: Final = {
    ALPHA: -2,
    BETA: 0,
    GAMMA: 0,
    DELTA: 0,
}
TWO_DOUBLE_COLLISION_FACTORIZATION: Final = (
    (S**2 + S - 1) ** 2 * (S**3 + 4 * S**2 + 2 * S + 2)
)
TWO_DOUBLE_BRANCH_EQUATION: Final = (
    -X**10
    + 10 * X**8
    + 15 * X**6 * Y
    - 20 * X**6
    + 8 * X**5 * Y
    + 12 * X**5
    - 50 * X**4 * Y
    + 10 * X**3 * Y**2
    + 30 * X**3 * Y
    - 25 * X**2 * Y**2
    + 20 * X * Y**2
    + Y**3
    - 3 * Y**2
)

HOSTILE_CONTACT_TWO_TRIPLE: Final = (
    A6_COLLISION_MERIDIAN,
    A6_TORUS_2_5_LOCAL_GENERATORS[0],
    A6_TORUS_2_5_LOCAL_GENERATORS[1],
)

# ``(456), (346), (123)`` is the analogous witness for contact order three.
HOSTILE_CONTACT_THREE_TRIPLE: Final = (
    A6_COLLISION_MERIDIAN,
    (0, 1, 3, 5, 4, 2),
    A6_TORUS_2_5_LOCAL_GENERATORS[1],
)


def _coefficient_substitution() -> dict[Symbol, Expr]:
    """Return the six lower coefficients of the family collision septic."""

    polynomial = expand(COLLISION_POLYNOMIAL)
    return {
        H0: polynomial.coeff(S, 0),
        H1: polynomial.coeff(S, 1),
        H2: polynomial.coeff(S, 2),
        H3: polynomial.coeff(S, 3),
        H4: polynomial.coeff(S, 4),
        H5: polynomial.coeff(S, 5),
    }


def _parameter_coefficient_matrix(derivative_orders: tuple[int, ...]) -> Matrix:
    """Coefficient matrix of ``H^(i)(r)`` in the four family parameters."""

    parameters = (ALPHA, BETA, GAMMA, DELTA)
    equations = tuple(
        diff(COLLISION_POLYNOMIAL, S, order).subs(S, R)
        for order in derivative_orders
    )
    return Matrix(
        [
            [expand(equation).coeff(parameter) for parameter in parameters]
            for equation in equations
        ]
    )


def _minor_gcd(matrix: Matrix) -> Expr:
    """Return the exact gcd of all maximal minors of a short matrix."""

    size = matrix.rows
    minors = tuple(
        expand(matrix[:, columns].det())
        for columns in combinations(range(matrix.cols), size)
    )
    result = minors[0]
    for minor in minors[1:]:
        result = gcd(result, minor)
    return factor(result)


@dataclass(frozen=True, slots=True)
class DiscriminantIncidenceCertificate:
    """Exact coefficient-slice and repeated-root incidence identities."""

    first_slice_identity: Expr
    second_slice_identity: Expr
    inverse_collision_identity: Expr
    contact_syzygy_identity: Expr
    double_incidence_remainders: tuple[Expr, Expr]
    double_rank_minor_gcd: Expr
    triple_incidence_remainders: tuple[Expr, Expr, Expr]
    triple_rank_minor_gcd: Expr
    two_double_factor_identity: Expr
    two_double_determinant_identity: Expr
    two_double_compatibility_eliminant: Expr
    two_double_cusp_boundary_gcd: Expr
    two_double_exceptional_boundary_gcd: Expr
    two_double_exceptional_slice_remainders: tuple[Expr, Expr]
    two_double_exceptional_factor_identity: Expr
    two_double_exceptional_critical_factor: Expr
    two_double_exceptional_triple_factor: Expr

    @property
    def verified(self) -> bool:
        """Whether every incidence and rank-boundary identity is exact."""

        return bool(
            self.first_slice_identity == 0
            and self.second_slice_identity == 0
            and self.inverse_collision_identity == 0
            and self.contact_syzygy_identity == 0
            and self.double_incidence_remainders == (0, 0)
            and self.double_rank_minor_gcd == R + 1
            and self.triple_incidence_remainders == (0, 0, 0)
            and self.triple_rank_minor_gcd == 2 * (R + 1) ** 3
            and self.two_double_factor_identity == 0
            and self.two_double_determinant_identity == 0
            and self.two_double_compatibility_eliminant
            == (ROOT_SUM + 2) ** 4 * (3 * ROOT_SUM + 5) ** 2
            and self.two_double_cusp_boundary_gcd
            == (ROOT_PRODUCT - 1) ** 2
            and expand(
                self.two_double_exceptional_boundary_gcd
                - (ROOT_PRODUCT - Rational(4, 9))
            )
            == 0
            and self.two_double_exceptional_slice_remainders == (0, 0)
            and self.two_double_exceptional_factor_identity == 0
            and self.two_double_exceptional_critical_factor == 0
            and self.two_double_exceptional_triple_factor == 0
        )


@cache
def exact_discriminant_incidence_certificate(
) -> DiscriminantIncidenceCertificate:
    """Build the exact coefficient and incidence normalization certificate."""

    coefficients = _coefficient_substitution()
    inverse_h = COLLISION_POLYNOMIAL.subs(COEFFICIENT_SLICE_INVERSE)
    abstract_h = (
        S**7
        + 6 * S**6
        + H5 * S**5
        + H4 * S**4
        + H3 * S**3
        + H2 * S**2
        + H1 * S
        + H0
    )
    inverse_difference = expand(inverse_h - abstract_h)
    inverse_reduced = expand(
        inverse_difference.subs(
            {
                H3: 4 * H4 - 10 * H5 + 85,
                H2: (-H0 + 2 * H1 + 11 * H4 - 34 * H5 + 306) / 3,
            }
        )
    )

    double_h = COLLISION_POLYNOMIAL.subs(DOUBLE_INCIDENCE_PARAMETERS)
    triple_h = COLLISION_POLYNOMIAL.subs(TRIPLE_INCIDENCE_PARAMETERS)
    repeated_factor = TWO_DOUBLE_REPEATED_FACTOR**2
    two_double_h = COLLISION_POLYNOMIAL.subs(TWO_DOUBLE_PARAMETERS)

    general_two_double_coefficients = {
        coefficient: expand(TWO_DOUBLE_GENERAL_POLYNOMIAL).coeff(S, degree)
        for coefficient, degree in (
            (H0, 0),
            (H1, 1),
            (H2, 2),
            (H3, 3),
            (H4, 4),
            (H5, 5),
        )
    }
    general_two_double_slice = (
        expand(COEFFICIENT_SLICE_FIRST.subs(general_two_double_coefficients)),
        expand(COEFFICIENT_SLICE_SECOND.subs(general_two_double_coefficients)),
    )
    compatibility_basis = groebner(
        [TWO_DOUBLE_DENOMINATOR, TWO_DOUBLE_N0, TWO_DOUBLE_N1],
        ROOT_PRODUCT,
        ROOT_SUM,
        order="lex",
    )
    exceptional_substitution = {
        ROOT_SUM: TWO_DOUBLE_EXCEPTIONAL_SUM,
        ROOT_PRODUCT: TWO_DOUBLE_EXCEPTIONAL_PRODUCT,
        TWO_DOUBLE_FREE_C1: TWO_DOUBLE_EXCEPTIONAL_C1,
    }
    exceptional_h = factor(
        TWO_DOUBLE_GENERAL_POLYNOMIAL.subs(exceptional_substitution)
    )
    exceptional_coefficients = {
        coefficient: expand(exceptional_h).coeff(S, degree)
        for coefficient, degree in (
            (H0, 0),
            (H1, 1),
            (H2, 2),
            (H3, 3),
            (H4, 4),
            (H5, 5),
        )
    }
    exceptional_parameters = {
        parameter: factor(expression.subs(exceptional_coefficients))
        for parameter, expression in COEFFICIENT_SLICE_INVERSE.items()
    }

    ordered_equations = tuple(
        diff(COLLISION_POLYNOMIAL, S, order).subs(S, root)
        for root in (R, QROOT)
        for order in (0, 1)
    )
    ordered_matrix = Matrix(
        [
            [expand(equation).coeff(parameter) for parameter in (ALPHA, BETA, GAMMA, DELTA)]
            for equation in ordered_equations
        ]
    )
    expected_ordered_determinant = -(
        (R + 1)
        * (QROOT + 1)
        * (QROOT - R) ** 4
        * TWO_DOUBLE_DENOMINATOR.subs(
            {ROOT_SUM: R + QROOT, ROOT_PRODUCT: R * QROOT}
        )
    )

    return DiscriminantIncidenceCertificate(
        first_slice_identity=expand(COEFFICIENT_SLICE_FIRST.subs(coefficients)),
        second_slice_identity=expand(COEFFICIENT_SLICE_SECOND.subs(coefficients)),
        inverse_collision_identity=inverse_reduced,
        contact_syzygy_identity=expand(
            (S + 1) * TANGENCY_POLYNOMIAL
            - S * (3 * S + 4) * diff(COLLISION_POLYNOMIAL, S)
            - (9 * S + 10) * COLLISION_POLYNOMIAL
        ),
        double_incidence_remainders=tuple(
            cancel(diff(double_h, S, order).subs(S, R))
            for order in range(2)
        ),
        double_rank_minor_gcd=_minor_gcd(_parameter_coefficient_matrix((0, 1))),
        triple_incidence_remainders=tuple(
            cancel(diff(triple_h, S, order).subs(S, R))
            for order in range(3)
        ),
        triple_rank_minor_gcd=_minor_gcd(
            _parameter_coefficient_matrix((0, 1, 2))
        ),
        two_double_factor_identity=cancel(
            two_double_h - repeated_factor * TWO_DOUBLE_RESIDUAL_FACTOR
        ),
        two_double_determinant_identity=expand(
            ordered_matrix.det() - expected_ordered_determinant
        ),
        two_double_compatibility_eliminant=factor(
            compatibility_basis.polys[-1].as_expr()
        ),
        two_double_cusp_boundary_gcd=factor(
            gcd(
                gcd(
                    TWO_DOUBLE_DENOMINATOR.subs(ROOT_SUM, -2),
                    TWO_DOUBLE_N0.subs(ROOT_SUM, -2),
                ),
                TWO_DOUBLE_N1.subs(ROOT_SUM, -2),
            )
        ),
        two_double_exceptional_boundary_gcd=factor(
            gcd(
                gcd(
                    TWO_DOUBLE_DENOMINATOR.subs(
                        ROOT_SUM, TWO_DOUBLE_EXCEPTIONAL_SUM
                    ),
                    TWO_DOUBLE_N0.subs(
                        ROOT_SUM, TWO_DOUBLE_EXCEPTIONAL_SUM
                    ),
                ),
                TWO_DOUBLE_N1.subs(ROOT_SUM, TWO_DOUBLE_EXCEPTIONAL_SUM),
            )
        ),
        two_double_exceptional_slice_remainders=tuple(
            expand(equation.subs(exceptional_substitution))
            for equation in general_two_double_slice
        ),
        two_double_exceptional_factor_identity=expand(
            exceptional_h - TWO_DOUBLE_EXCEPTIONAL_FACTORIZATION
        ),
        two_double_exceptional_critical_factor=expand(
            EXTRA_CRITICAL_FACTOR.subs(exceptional_parameters)
        ),
        two_double_exceptional_triple_factor=expand(
            TRIPLE_COLLISION_FACTOR.subs(exceptional_parameters)
        ),
    )


@dataclass(frozen=True, slots=True)
class DiscriminantWallSampleCertificate:
    """Exact geometry and finite-image census for one wall representative."""

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
    arithmetic_genus: int
    total_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether geometry, genus, and ``A6`` replay all agree."""

        return bool(
            self.collision_identity == 0
            and self.implicit_identity == 0
            and self.cusp_factor != 0
            and self.critical_factor != 0
            and self.triple_factor != 0
            and self.residual_separation != 0
            and self.arithmetic_genus == 36
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
    relations: tuple[tuple[int, ...], ...],
    partition: tuple[int, ...],
    contact_orders: tuple[int, ...],
    node_count: int,
    residual_separation: Expr,
) -> DiscriminantWallSampleCertificate:
    """Build one exact representative certificate."""

    sample_q = FAMILY_Q.subs(parameters)
    return DiscriminantWallSampleCertificate(
        name=name,
        partition=partition,
        contact_orders=contact_orders,
        node_count=node_count,
        collision_identity=expand(
            COLLISION_POLYNOMIAL.subs(parameters) - factorization
        ),
        implicit_identity=expand(
            branch_equation.subs({X: FAMILY_P, Y: sample_q})
        ),
        cusp_factor=expand(CUSP_COLLISION_FACTOR.subs(parameters)),
        critical_factor=expand(EXTRA_CRITICAL_FACTOR.subs(parameters)),
        triple_factor=expand(TRIPLE_COLLISION_FACTOR.subs(parameters)),
        residual_separation=residual_separation,
        arithmetic_genus=36,
        total_delta=27 + 2 + sum(contact_orders) + node_count,
        complement_census=three_cycle_presentation_census(relations),
    )


@cache
def exact_generic_g_wall_certificate() -> DiscriminantWallSampleCertificate:
    """Return the contact-two plus five-node representative."""

    squarefree_quotient = (
        (S**2 + S - 1) * (S**3 - S**2 + 2 * S - 3)
    )
    return _sample_certificate(
        name="one contact two plus five nodes",
        parameters=GENERIC_G_PARAMETERS,
        factorization=GENERIC_G_COLLISION_FACTORIZATION,
        branch_equation=GENERIC_G_BRANCH_EQUATION,
        relations=GENERIC_G_RELATIONS,
        partition=(2, 1, 1, 1, 1, 1),
        contact_orders=(2,),
        node_count=5,
        residual_separation=discriminant(squarefree_quotient, S)
        * squarefree_quotient.subs(S, -3),
    )


@cache
def exact_contact_three_wall_certificate() -> DiscriminantWallSampleCertificate:
    """Return the contact-three plus four-node representative."""

    quartic = (4 * S**4 + 3 * S**2 - 4 * S + 2) / 4
    return _sample_certificate(
        name="one contact three plus four nodes",
        parameters=CONTACT_THREE_PARAMETERS,
        factorization=CONTACT_THREE_COLLISION_FACTORIZATION,
        branch_equation=CONTACT_THREE_BRANCH_EQUATION,
        relations=CONTACT_THREE_RELATIONS,
        partition=(3, 1, 1, 1, 1),
        contact_orders=(3,),
        node_count=4,
        residual_separation=discriminant(quartic, S) * quartic.subs(S, -2),
    )


@cache
def exact_two_double_wall_certificate() -> DiscriminantWallSampleCertificate:
    """Return the two-contact-two plus three-node representative."""

    repeated = S**2 + S - 1
    cubic = S**3 + 4 * S**2 + 2 * S + 2
    return _sample_certificate(
        name="two contact-two points plus three nodes",
        parameters=TWO_DOUBLE_SAMPLE_PARAMETERS,
        factorization=TWO_DOUBLE_COLLISION_FACTORIZATION,
        branch_equation=TWO_DOUBLE_BRANCH_EQUATION,
        relations=TWO_DOUBLE_RELATIONS,
        partition=(2, 2, 1, 1, 1),
        contact_orders=(2, 2),
        node_count=3,
        residual_separation=(
            discriminant(repeated, S)
            * discriminant(cubic, S)
            * resultant(repeated, cubic, S)
        ),
    )


def contact_relation_holds(left: Permutation, right: Permutation, order: int) -> bool:
    """Whether ``(xy)^order = (yx)^order`` for a two-branch contact."""

    if order < 1:
        msg = "a contact order must be positive"
        raise ValueError(msg)
    return permutation_power(compose(left, right), order) == permutation_power(
        compose(right, left), order
    )


@dataclass(frozen=True, slots=True)
class ContactRelationCensus:
    """Hostile ``A6`` census for one isolated local contact relation."""

    order: int
    satisfying_pairs: int
    a6_generating_triples: int


@cache
def contact_relation_census(order: int) -> ContactRelationCensus:
    """Exhaust contact pairs and all choices of a third three-cycle."""

    three_cycles = tuple(
        element
        for element in alternating_group_six()
        if cycle_type(element) == (3, 1, 1, 1)
    )
    pairs = tuple(
        (left, right)
        for left, right in product(three_cycles, repeat=2)
        if contact_relation_holds(left, right, order)
    )
    a6_triples = sum(
        len(generated_group((left, right, third))) == 360
        for left, right in pairs
        for third in three_cycles
    )
    return ContactRelationCensus(
        order=order,
        satisfying_pairs=len(pairs),
        a6_generating_triples=a6_triples,
    )


@cache
def two_contact_a6_assignments(first_order: int, second_order: int) -> int:
    """Count generating triples satisfying contacts on pairs ``(x,y),(x,z)``."""

    three_cycles = tuple(
        element
        for element in alternating_group_six()
        if cycle_type(element) == (3, 1, 1, 1)
    )
    return sum(
        contact_relation_holds(first, second, first_order)
        and contact_relation_holds(first, third, second_order)
        and len(generated_group((first, second, third))) == 360
        for first, second, third in product(three_cycles, repeat=3)
    )


def collision_partitions(excess: int) -> tuple[tuple[int, ...], ...]:
    """Return partitions of seven roots with ``sum(m_i - 1) = excess``."""

    if excess < 0 or excess > 6:
        msg = "collision excess must lie between zero and six"
        raise ValueError(msg)

    def descend(total: int, maximum: int) -> tuple[tuple[int, ...], ...]:
        if total == 0:
            return ((),)
        result: list[tuple[int, ...]] = []
        for head in range(min(total, maximum), 0, -1):
            result.extend(
                (head, *tail)
                for tail in descend(total - head, head)
            )
        return tuple(result)

    length = 7 - excess
    return tuple(partition for partition in descend(7, 7) if len(partition) == length)


def main() -> None:
    """Print the exact first-stage discriminant-wall checkpoint."""

    incidence = exact_discriminant_incidence_certificate()
    samples = (
        exact_generic_g_wall_certificate(),
        exact_contact_three_wall_certificate(),
        exact_two_double_wall_certificate(),
    )
    print("incidence algebra verified:", incidence.verified)
    print("sample complements:", {sample.name: sample.complement_census for sample in samples})
    print("all three strata excluded:", all(sample.verified for sample in samples))
    print("hostile contact two:", contact_relation_census(2))
    print("hostile contact three:", contact_relation_census(3))
    print("two-contact A6 assignments:", two_contact_a6_assignments(2, 2))


if __name__ == "__main__":
    main()
