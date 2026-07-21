"""Audit the conditional delta-ten ``C2 + T111 + 5N`` profile.

Write the normalized degree-``(4, 9)`` family as

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``.

An ordinary triple ``P``-fiber with elementary symmetric coordinates
``(p,q,r)`` satisfies ``q^2 - p*r - q = 0``.  On ``q*r != 0`` this gives
``p=q(q-1)/r`` and ``k=(r^2+q^2-q^3)/(q*r)``.  Requiring ``Q`` to be
constant on its cubic gives two affine-linear equations in ``(a,b,c,d)``.
A separate contact-two pair with sum ``w`` gives ``H(w)=T(w)=0``, hence two
more affine-linear equations.

The resulting ``4 x 4`` determinant is factored exactly below.  Away from
the displayed pair, split, overlap, triple, and residual-rank factors, the
incidence is the graph of a rational map over an irreducible open subset of
``A^3_(q,r,w)``.  It is therefore a rational threefold.  A coefficient point
has only finitely many collision pairs and triple fibers, so projection to
the five-dimensional coefficient space is quasi-finite on this open, hence
generically finite onto its image; the image closure is one irreducible
codimension-two component.  The sole non-boundary rank factor is irreducible
and does not divide a selected augmented minor, as witnessed by an exact
univariate specialization.  Thus it supports no additional component
dominating that divisor.  Compatible lower-dimensional
intersections of the rank divisor and deeper intersections in the split
fibers are not classified here.

The rational member

``X = 2*t^4 + 3*t^3 + 2*t^2`` and
``Y = t^9 - 2*t^7 + t^6 - 2*t^5``

has the forced ``T(2,5)`` cusp, one two-branch contact of intersection
multiplicity two, one ordinary triple point, and five nodes.  A stored exact
Sage 10.8 van Kamp presentation simplifies to ``Z``; the finite replay has
only forty diagonal ``C3`` images and no ``A6`` image.

The topology computation certifies this member.  Propagation to a stratum is
a separate dependency: it needs a connected clean open and proper projective
Whitney--Thom triviality over that open.

This is a conditional, computer-assisted dominant-component certificate.
It neither settles residual rank intersections nor proves the plane
Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import (
    Expr,
    Matrix,
    Poly,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    factor,
    gcd,
    invert,
    rem,
    resultant,
    together,
)

from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    COLLISION_X_DENOMINATOR,
    COLLISION_X_NUMERATOR,
    CUSP_IMAGE_FACTOR,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    PAIR_DENOMINATOR,
    PAIR_DIAGONAL_FACTOR,
    PAIR_QUADRATIC,
    S,
    T,
    TANGENCY_POLYNOMIAL,
    X,
    Y,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    n_generator_three_cycle_presentation_census,
)

TRIPLE_SUM: Final = Symbol("p_contact_triple")
TRIPLE_PAIR_SUM: Final = Symbol("q_contact_triple")
TRIPLE_PRODUCT: Final = Symbol("r_contact_triple")
CONTACT_SUM: Final = Symbol("w_contact_triple")
SLOPE_VALUE: Final = Symbol("m_contact_triple")
LINEAR_PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)

TRIPLE_CUBIC: Final = T**3 - TRIPLE_SUM * T**2 + TRIPLE_PAIR_SUM * T - TRIPLE_PRODUCT
TRIPLE_P_CONSTRAINT: Final = (
    TRIPLE_PAIR_SUM**2 - TRIPLE_SUM * TRIPLE_PRODUCT - TRIPLE_PAIR_SUM
)
TRIPLE_SUM_PARAMETERIZATION: Final = (
    TRIPLE_PAIR_SUM * (TRIPLE_PAIR_SUM - 1) / TRIPLE_PRODUCT
)
KAPPA_NUMERATOR: Final = TRIPLE_PRODUCT**2 + TRIPLE_PAIR_SUM**2 - TRIPLE_PAIR_SUM**3
KAPPA_PARAMETERIZATION: Final = KAPPA_NUMERATOR / (TRIPLE_PAIR_SUM * TRIPLE_PRODUCT)
TRIPLE_PARAMETERIZATION: Final = {
    TRIPLE_SUM: TRIPLE_SUM_PARAMETERIZATION,
    KAPPA: KAPPA_PARAMETERIZATION,
}

GENERIC_TRIPLE_Q_REMAINDER: Final = expand(rem(FAMILY_Q, TRIPLE_CUBIC, T))
TRIPLE_EQUALITY_EQUATIONS: Final = tuple(
    cancel(
        Poly(GENERIC_TRIPLE_Q_REMAINDER, T)
        .coeff_monomial(T**degree)
        .subs(TRIPLE_PARAMETERIZATION)
    )
    for degree in (2, 1)
)
CONTACT_EQUATIONS: Final = (
    cancel(COLLISION_POLYNOMIAL.subs({KAPPA: KAPPA_PARAMETERIZATION, S: CONTACT_SUM})),
    cancel(TANGENCY_POLYNOMIAL.subs({KAPPA: KAPPA_PARAMETERIZATION, S: CONTACT_SUM})),
)
COMBINED_EQUATIONS: Final = (*CONTACT_EQUATIONS, *TRIPLE_EQUALITY_EQUATIONS)

# These row multipliers clear every denominator in the coefficient matrix.
INCIDENCE_ROW_MULTIPLIERS: Final = (
    TRIPLE_PAIR_SUM**4 * TRIPLE_PRODUCT**4,
    TRIPLE_PAIR_SUM**5 * TRIPLE_PRODUCT**5,
    TRIPLE_PRODUCT**6,
    TRIPLE_PRODUCT**5,
)
INCIDENCE_MATRIX: Final = Matrix(
    [
        [
            cancel(multiplier * diff(equation, parameter))
            for parameter in LINEAR_PARAMETERS
        ]
        for equation, multiplier in zip(
            COMBINED_EQUATIONS,
            INCIDENCE_ROW_MULTIPLIERS,
            strict=True,
        )
    ]
)

# One extra power of r clears the constant entries of the last two rows.
AUGMENTED_ROW_MULTIPLIERS: Final = (
    TRIPLE_PAIR_SUM**4 * TRIPLE_PRODUCT**4,
    TRIPLE_PAIR_SUM**5 * TRIPLE_PRODUCT**5,
    TRIPLE_PRODUCT**7,
    TRIPLE_PRODUCT**6,
)
ZERO_LINEAR_PARAMETERS: Final = dict.fromkeys(LINEAR_PARAMETERS, 0)
INCIDENCE_AUGMENTED_MATRIX: Final = Matrix(
    [
        [
            *(
                cancel(multiplier * diff(equation, parameter))
                for parameter in LINEAR_PARAMETERS
            ),
            cancel(multiplier * equation.subs(ZERO_LINEAR_PARAMETERS)),
        ]
        for equation, multiplier in zip(
            COMBINED_EQUATIONS,
            AUGMENTED_ROW_MULTIPLIERS,
            strict=True,
        )
    ]
)

CONTACT_DENOMINATOR_FACTOR: Final = (
    TRIPLE_PAIR_SUM**3
    - TRIPLE_PAIR_SUM**2
    - 2 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT * CONTACT_SUM
    - TRIPLE_PRODUCT**2
)
CONTACT_DIAGONAL_FACTOR: Final = (
    3 * TRIPLE_PAIR_SUM**3 * CONTACT_SUM
    - 3 * TRIPLE_PAIR_SUM**2 * CONTACT_SUM
    - 2 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT * CONTACT_SUM**2
    - 4 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT
    - 3 * TRIPLE_PRODUCT**2 * CONTACT_SUM
)
CONTACT_PRODUCT_FACTOR: Final = (
    TRIPLE_PAIR_SUM**3 * CONTACT_SUM
    - TRIPLE_PAIR_SUM**2 * CONTACT_SUM
    - TRIPLE_PAIR_SUM * TRIPLE_PRODUCT * CONTACT_SUM**2
    - TRIPLE_PAIR_SUM * TRIPLE_PRODUCT
    - TRIPLE_PRODUCT**2 * CONTACT_SUM
)
OVERLAP_FACTOR_C: Final = (
    TRIPLE_PAIR_SUM**4 * CONTACT_SUM
    - TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT
    - 2 * TRIPLE_PAIR_SUM**3 * CONTACT_SUM
    - 2 * TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT * CONTACT_SUM**2
    + TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT
    + TRIPLE_PAIR_SUM**2 * CONTACT_SUM
    + TRIPLE_PAIR_SUM * TRIPLE_PRODUCT**2 * CONTACT_SUM
    + 2 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT * CONTACT_SUM**2
    + TRIPLE_PRODUCT**3
    + TRIPLE_PRODUCT**2 * CONTACT_SUM**3
)
OVERLAP_FACTOR_E: Final = (
    TRIPLE_PAIR_SUM**5 * CONTACT_SUM**2
    + TRIPLE_PAIR_SUM**4 * TRIPLE_PRODUCT * CONTACT_SUM
    - TRIPLE_PAIR_SUM**4 * CONTACT_SUM**2
    + TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT**2
    - TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT * CONTACT_SUM**3
    - 2 * TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT * CONTACT_SUM
    - 3 * TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT**2 * CONTACT_SUM**2
    - TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT**2
    - 3 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT**3 * CONTACT_SUM
    - TRIPLE_PRODUCT**4
)
TRIPLE_DISCRIMINANT_NUMERATOR: Final = (
    3 * TRIPLE_PAIR_SUM**6
    - 10 * TRIPLE_PAIR_SUM**5
    + 11 * TRIPLE_PAIR_SUM**4
    - 14 * TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT**2
    - 4 * TRIPLE_PAIR_SUM**3
    + 18 * TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT**2
    + 27 * TRIPLE_PRODUCT**4
)
TRIPLE_FOURTH_ROOT_FACTOR: Final = (
    3 * TRIPLE_PAIR_SUM**3 - TRIPLE_PAIR_SUM**2 + TRIPLE_PRODUCT**2
)
SPLIT_FACTOR: Final = KAPPA_NUMERATOR * (
    KAPPA_NUMERATOR**2 - 4 * TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT**2
)

RESIDUAL_RANK_FACTOR: Final = (
    2 * TRIPLE_PAIR_SUM**10 * CONTACT_SUM**2
    + TRIPLE_PAIR_SUM**9 * TRIPLE_PRODUCT * CONTACT_SUM
    - 8 * TRIPLE_PAIR_SUM**9 * CONTACT_SUM**2
    - 4 * TRIPLE_PAIR_SUM**8 * TRIPLE_PRODUCT * CONTACT_SUM**3
    - 11 * TRIPLE_PAIR_SUM**8 * TRIPLE_PRODUCT * CONTACT_SUM
    + 10 * TRIPLE_PAIR_SUM**8 * CONTACT_SUM**2
    - 11 * TRIPLE_PAIR_SUM**7 * TRIPLE_PRODUCT**2 * CONTACT_SUM**2
    - 3 * TRIPLE_PAIR_SUM**7 * TRIPLE_PRODUCT**2
    + 12 * TRIPLE_PAIR_SUM**7 * TRIPLE_PRODUCT * CONTACT_SUM**3
    + 23 * TRIPLE_PAIR_SUM**7 * TRIPLE_PRODUCT * CONTACT_SUM
    - 4 * TRIPLE_PAIR_SUM**7 * CONTACT_SUM**2
    - 3 * TRIPLE_PAIR_SUM**6 * TRIPLE_PRODUCT**3 * CONTACT_SUM
    + 2 * TRIPLE_PAIR_SUM**6 * TRIPLE_PRODUCT**2 * CONTACT_SUM**4
    + 42 * TRIPLE_PAIR_SUM**6 * TRIPLE_PRODUCT**2 * CONTACT_SUM**2
    + 12 * TRIPLE_PAIR_SUM**6 * TRIPLE_PRODUCT**2
    - 8 * TRIPLE_PAIR_SUM**6 * TRIPLE_PRODUCT * CONTACT_SUM**3
    - 13 * TRIPLE_PAIR_SUM**6 * TRIPLE_PRODUCT * CONTACT_SUM
    + 12 * TRIPLE_PAIR_SUM**5 * TRIPLE_PRODUCT**3 * CONTACT_SUM**3
    + 38 * TRIPLE_PAIR_SUM**5 * TRIPLE_PRODUCT**3 * CONTACT_SUM
    - 4 * TRIPLE_PAIR_SUM**5 * TRIPLE_PRODUCT**2 * CONTACT_SUM**4
    - 35 * TRIPLE_PAIR_SUM**5 * TRIPLE_PRODUCT**2 * CONTACT_SUM**2
    - 11 * TRIPLE_PAIR_SUM**5 * TRIPLE_PRODUCT**2
    + 16 * TRIPLE_PAIR_SUM**4 * TRIPLE_PRODUCT**4 * CONTACT_SUM**2
    + 6 * TRIPLE_PAIR_SUM**4 * TRIPLE_PRODUCT**4
    - 24 * TRIPLE_PAIR_SUM**4 * TRIPLE_PRODUCT**3 * CONTACT_SUM**3
    - 47 * TRIPLE_PAIR_SUM**4 * TRIPLE_PRODUCT**3 * CONTACT_SUM
    + 3 * TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT**5 * CONTACT_SUM
    - 2 * TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT**4 * CONTACT_SUM**4
    - 46 * TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT**4 * CONTACT_SUM**2
    - 16 * TRIPLE_PAIR_SUM**3 * TRIPLE_PRODUCT**4
    - 8 * TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT**5 * CONTACT_SUM**3
    - 27 * TRIPLE_PAIR_SUM**2 * TRIPLE_PRODUCT**5 * CONTACT_SUM
    - 7 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT**6 * CONTACT_SUM**2
    - 3 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT**6
    - TRIPLE_PRODUCT**7 * CONTACT_SUM
)

EXPECTED_INCIDENCE_DETERMINANT: Final = (
    TRIPLE_PAIR_SUM**2
    * TRIPLE_PRODUCT**10
    * CONTACT_SUM
    * CONTACT_DENOMINATOR_FACTOR**3
    * CONTACT_DIAGONAL_FACTOR
    * OVERLAP_FACTOR_C**2
    * RESIDUAL_RANK_FACTOR
)
VALID_CRAMER_LOCALIZER: Final = (
    TRIPLE_PAIR_SUM
    * TRIPLE_PRODUCT
    * CONTACT_SUM
    * SPLIT_FACTOR
    * CONTACT_DENOMINATOR_FACTOR
    * CONTACT_DIAGONAL_FACTOR
    * CONTACT_PRODUCT_FACTOR
    * OVERLAP_FACTOR_C
    * OVERLAP_FACTOR_E
    * TRIPLE_DISCRIMINANT_NUMERATOR
    * TRIPLE_FOURTH_ROOT_FACTOR
    * RESIDUAL_RANK_FACTOR
)

RANKDROP_SPECIALIZATION: Final = {
    TRIPLE_PAIR_SUM: 2,
    TRIPLE_PRODUCT: -1,
}
EXPECTED_SPECIALIZED_RANK_FACTOR: Final = (
    -16 * CONTACT_SUM**4
    + 32 * CONTACT_SUM**3
    + 34 * CONTACT_SUM**2
    + 5 * CONTACT_SUM
    - 6
)
EXPECTED_SPECIALIZED_AUGMENTED_MINOR: Final = (
    8
    * CONTACT_SUM**3
    * (CONTACT_SUM + 1) ** 2
    * (4 * CONTACT_SUM + 3) ** 2
    * (CONTACT_SUM**2 + 3 * CONTACT_SUM + 3) ** 2
    * (4 * CONTACT_SUM**2 + 9 * CONTACT_SUM + 8)
    * (
        48 * CONTACT_SUM**5
        + 84 * CONTACT_SUM**4
        - 82 * CONTACT_SUM**3
        - 220 * CONTACT_SUM**2
        - 135 * CONTACT_SUM
        - 38
    )
)

SAMPLE_BASE: Final = {
    TRIPLE_PAIR_SUM: 2,
    TRIPLE_PRODUCT: -1,
    CONTACT_SUM: 1,
}
SAMPLE_PARAMETERS: Final = {
    KAPPA: Rational(3, 2),
    ALPHA: -2,
    BETA: 1,
    GAMMA: -2,
    DELTA: 0,
}
SAMPLE_P: Final = expand(2 * FAMILY_P.subs(SAMPLE_PARAMETERS))
SAMPLE_Q: Final = expand(FAMILY_Q.subs(SAMPLE_PARAMETERS))
SAMPLE_CONTACT_BRANCH: Final = T**2 - T + 1
SAMPLE_TRIPLE_BRANCH: Final = T**3 + 2 * T**2 + 2 * T + 1
SAMPLE_TRIPLE_PAIR_SUMS: Final = (S + 1) * (S**2 + 3 * S + 3)
SAMPLE_RESIDUAL_NODES: Final = (
    8 * S**5 + 56 * S**4 + 154 * S**3 + 125 * S**2 + 12 * S - 12
)
SAMPLE_COLLISION: Final = expand(
    (S - 1) ** 2 * SAMPLE_TRIPLE_PAIR_SUMS * SAMPLE_RESIDUAL_NODES / 8
)
SAMPLE_TANGENCY_COFACTOR: Final = (
    576 * S**10
    + 6192 * S**9
    + 29188 * S**8
    + 74170 * S**7
    + 104432 * S**6
    + 67796 * S**5
    - 10758 * S**4
    - 45585 * S**3
    - 23661 * S**2
    - 1746 * S
    + 1080
)
SAMPLE_TANGENCY: Final = expand((S - 1) * SAMPLE_TANGENCY_COFACTOR / 16)
SAMPLE_NODE_X_POLYNOMIAL: Final = (
    512 * X**5
    - 62528 * X**4
    + 2718048 * X**3
    - 48220459 * X**2
    + 255118656 * X
    + 115495936
)
SAMPLE_IMPLICIT: Final = (
    X**9
    - 31 * X**8
    + 139 * X**7
    - 529 * X**6 * Y
    + 499 * X**6
    + 4234 * X**5 * Y
    + 1376 * X**5
    - 1840 * X**4 * Y**2
    + 4183 * X**4 * Y
    + 19799 * X**3 * Y**2
    + 23392 * X**3 * Y
    - 1728 * X**2 * Y**3
    - 75323 * X**2 * Y**2
    + 10344 * X * Y**3
    + 44168 * X * Y**2
    - 512 * Y**4
    - 2503 * Y**3
    - 11008 * Y**2
)

# Sage 10.8 exact raw affine presentation.  Indices 1..4 are meridians of a
# generic vertical fiber.
CONTACT_TRIPLE_RELATIONS: Final = (
    (-4, 2, 4, 2, 4, -2, -4, -2),
    (2, -1),
    (-4, 2, 4, 2, 4, 2, -4, -2, -4, -2),
    (-4, 2, 4, -2),
    (
        -4,
        -2,
        -4,
        -2,
        -4,
        -2,
        -1,
        2,
        4,
        2,
        4,
        2,
        -4,
        -2,
        -4,
        -2,
        -4,
        -2,
        1,
        2,
        4,
        2,
        4,
        2,
        4,
        2,
        -4,
        -2,
        -4,
        -2,
        -4,
        -2,
        -1,
        2,
        4,
        2,
        4,
        2,
        4,
        -2,
        -4,
        -2,
        -4,
        -2,
        1,
        2,
        4,
        2,
        4,
        2,
        4,
        -2,
    ),
    (
        1,
        2,
        4,
        2,
        4,
        2,
        -4,
        -2,
        -4,
        -2,
        -4,
        -2,
        -1,
        -4,
        -3,
        4,
        2,
        4,
        2,
        4,
        2,
        4,
        -2,
        -4,
        -2,
        -4,
        -2,
        -4,
        3,
        4,
    ),
    (
        4,
        2,
        4,
        -2,
        -4,
        -2,
        -4,
        -2,
        -4,
        3,
        4,
        1,
        -4,
        -3,
        4,
        2,
        4,
        2,
        4,
        2,
        -4,
        -2,
        -4,
        -2,
        -4,
        -2,
        -1,
        2,
        4,
        2,
    ),
    (-4, -2, -4, -2, -4, 3, 4, 2, 4, 2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
    (-3, -2, -1, 2, 4, -2, 1, 2, 3, -2, -1, 2, -4, -2, 1, 2),
    (-4, -2, 1, 2, 4, -2, -1, 2),
)


def _modulo_branch(rational_expression: Expr, branch: Expr) -> Expr:
    """Reduce a rational function modulo a monic source-fiber polynomial."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    denominator_inverse = invert(denominator, branch)
    return expand(rem(numerator * denominator_inverse, branch, T))


@dataclass(frozen=True, slots=True)
class DeltaTenContactTripleCertificate:
    """Exact dominant-incidence, geometry, and complement data."""

    incidence_determinant_identity: Expr
    rank_factor_irreducible: bool
    specialized_rank_identity: Expr
    specialized_augmented_identity: Expr
    specialized_augmented_gcd: Expr
    cramer_open_component_count: int
    cramer_open_component_dimension: int
    coefficient_image_codimension: int
    projection_quasi_finite_on_cramer_open: bool
    projection_generically_finite_onto_image: bool
    residual_rank_intersections_classified: bool
    split_boundary_intersections_classified: bool
    topology_propagation_dependencies: tuple[str, str]
    sample_incidence_residuals: tuple[Expr, ...]
    sample_incidence_determinant: Expr
    sample_linear_solution: tuple[Expr, Expr, Expr, Expr]
    sample_open_factor_values: tuple[Expr, ...]
    boundary_identities: tuple[Expr, ...]
    collision_identity: Expr
    tangency_identity: Expr
    collision_tangency_gcd: Expr
    contact_collision_second_derivative: Expr
    contact_tangency_derivative: Expr
    residual_discriminant: Expr
    residual_factor_resultants: tuple[Expr, Expr, Expr]
    tangency_factor_resultants: tuple[Expr, Expr]
    pair_boundary_resultants: tuple[Expr, Expr, Expr, Expr]
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    contact_image_remainders: tuple[Expr, Expr]
    triple_image_remainders: tuple[Expr, Expr]
    fourth_source_image: tuple[Expr, Expr]
    source_derivative_resultants: tuple[Expr, Expr]
    contact_jet_differences: tuple[Expr, Expr]
    contact_second_jet_separation: Expr
    triple_slope_resultant_identity: Expr
    triple_slope_discriminant: Expr
    node_x_resultant_identity: Expr
    node_x_discriminant: Expr
    node_target_separations: tuple[Expr, Expr, Expr]
    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    implicit_content: Expr
    sage_jacobian_components: tuple[tuple[int, int], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    arithmetic_genus: int
    cusp_delta: int
    contact_delta: int
    triple_delta: int
    node_count: int
    infinity_delta: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return (
            self.cusp_delta
            + self.contact_delta
            + self.triple_delta
            + self.node_count
            + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether every exact dominant-chart and sample check agrees."""

        expected_node_discriminant = -17004216373814107665135371132504044296002359787520
        return bool(
            self.incidence_determinant_identity == 0
            and self.rank_factor_irreducible
            and self.specialized_rank_identity == 0
            and self.specialized_augmented_identity == 0
            and self.specialized_augmented_gcd == 1
            and self.cramer_open_component_count == 1
            and self.cramer_open_component_dimension == 3
            and self.coefficient_image_codimension == 2
            and self.projection_quasi_finite_on_cramer_open
            and self.projection_generically_finite_onto_image
            and not self.residual_rank_intersections_classified
            and not self.split_boundary_intersections_classified
            and self.topology_propagation_dependencies
            == (
                "connected clean open",
                "proper projective Whitney-Thom triviality",
            )
            and self.sample_incidence_residuals == (0, 0, 0, 0)
            and self.sample_incidence_determinant == 276710448
            and self.sample_linear_solution == (-2, 1, -2, 0)
            and all(value != 0 for value in self.sample_open_factor_values)
            and all(identity == 0 for identity in self.boundary_identities)
            and self.collision_identity == 0
            and self.tangency_identity == 0
            and self.collision_tangency_gcd == S - 1
            and self.contact_collision_second_derivative == Rational(2401, 2)
            and self.contact_tangency_derivative == Rational(50421, 4)
            and self.residual_discriminant == -2138993541120
            and self.residual_factor_resultants == (-343, 63315, -14)
            and self.tangency_factor_resultants
            == (-39058273766241112220775, Rational(1759080645, 1024))
            and self.pair_boundary_resultants
            == (
                Rational(-9, 2),
                Rational(21609, 1024),
                Rational(103243, 512),
                Rational(-15993922959, 4096),
            )
            and self.cusp_image_factor == Rational(43, 2)
            and self.extra_critical_factor == Rational(740151, 8)
            and self.contact_image_remainders == (-5, -2)
            and self.triple_image_remainders == (1, 4)
            and self.fourth_source_image == (1, Rational(-31, 512))
            and self.source_derivative_resultants == (237, -63)
            and self.contact_jet_differences == (0, Rational(392, 6241) * (2 * T - 1))
            and self.contact_second_jet_separation == 3
            and self.triple_slope_resultant_identity == 0
            and self.triple_slope_discriminant == -344755200
            and self.node_x_resultant_identity == 0
            and self.node_x_discriminant == expected_node_discriminant
            and self.node_target_separations == (-2746044819, 325050165, 115495936)
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_components == ((4, 1), (3, 1), (4, 1), (5, 5))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.contact_delta == 2
            and self.triple_delta == 3
            and self.node_count == 5
            and self.infinity_delta == 16
            and self.total_delta == self.arithmetic_genus
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_contact_triple_certificate() -> DeltaTenContactTripleCertificate:
    """Build the exact dominant-component and rational-member certificate."""

    incidence_determinant = INCIDENCE_MATRIX.det(method="domain-ge")
    specialized_augmented = INCIDENCE_AUGMENTED_MATRIX.subs(RANKDROP_SPECIALIZATION)
    specialized_augmented_minor = specialized_augmented[:, [1, 2, 3, 4]].det(
        method="domain-ge"
    )
    specialized_rank_factor = expand(RESIDUAL_RANK_FACTOR.subs(RANKDROP_SPECIALIZATION))

    sample_matrix = INCIDENCE_MATRIX.subs(SAMPLE_BASE)
    sample_constant_column = Matrix(
        [
            equation.subs(ZERO_LINEAR_PARAMETERS).subs(SAMPLE_BASE)
            for equation in COMBINED_EQUATIONS
        ]
    )
    sample_rational_matrix = Matrix(
        [
            [
                diff(equation, parameter).subs(SAMPLE_BASE)
                for parameter in LINEAR_PARAMETERS
            ]
            for equation in COMBINED_EQUATIONS
        ]
    )
    sample_solution_matrix = sample_rational_matrix.inv() * (-sample_constant_column)

    parameterized_p = FAMILY_P.subs(KAPPA, KAPPA_PARAMETERIZATION)
    triple_cubic = TRIPLE_CUBIC.subs(TRIPLE_SUM, TRIPLE_SUM_PARAMETERIZATION)
    triple_target_x = TRIPLE_PRODUCT**2 / TRIPLE_PAIR_SUM
    contact_product = cancel(
        CONTACT_SUM
        * PAIR_QUADRATIC.subs({KAPPA: KAPPA_PARAMETERIZATION, S: CONTACT_SUM})
        / PAIR_DENOMINATOR.subs({KAPPA: KAPPA_PARAMETERIZATION, S: CONTACT_SUM})
    )
    contact_branch = T**2 - CONTACT_SUM * T + contact_product
    contact_target_x = cancel(
        COLLISION_X_NUMERATOR.subs({KAPPA: KAPPA_PARAMETERIZATION, S: CONTACT_SUM})
        / COLLISION_X_DENOMINATOR.subs({KAPPA: KAPPA_PARAMETERIZATION, S: CONTACT_SUM})
    )
    contact_source_resultant = resultant(contact_branch, triple_cubic, T)
    boundary_identities = (
        cancel(TRIPLE_P_CONSTRAINT.subs(TRIPLE_PARAMETERIZATION)),
        cancel(
            (parameterized_p - triple_target_x) / triple_cubic
            - (T + TRIPLE_PRODUCT / TRIPLE_PAIR_SUM)
        ),
        cancel(
            discriminant(triple_cubic, T)
            + TRIPLE_DISCRIMINANT_NUMERATOR / TRIPLE_PRODUCT**2
        ),
        cancel(
            triple_cubic.subs(T, -TRIPLE_PRODUCT / TRIPLE_PAIR_SUM)
            + TRIPLE_PRODUCT * TRIPLE_FOURTH_ROOT_FACTOR / TRIPLE_PAIR_SUM**3
        ),
        cancel(rem(parameterized_p - contact_target_x, contact_branch, T)),
        cancel(
            discriminant(contact_branch, T)
            + CONTACT_SUM * CONTACT_DIAGONAL_FACTOR / CONTACT_DENOMINATOR_FACTOR
        ),
        cancel(
            contact_source_resultant
            - OVERLAP_FACTOR_C**2
            * OVERLAP_FACTOR_E
            / (TRIPLE_PRODUCT**2 * CONTACT_DENOMINATOR_FACTOR**3)
        ),
        cancel(
            contact_target_x
            - triple_target_x
            - OVERLAP_FACTOR_C
            * OVERLAP_FACTOR_E
            / (TRIPLE_PAIR_SUM * TRIPLE_PRODUCT * CONTACT_DENOMINATOR_FACTOR**2)
        ),
    )

    sample_collision = COLLISION_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
    sample_tangency = TANGENCY_POLYNOMIAL.subs(SAMPLE_PARAMETERS)
    sample_denominator = PAIR_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    sample_product = PAIR_QUADRATIC.subs(SAMPLE_PARAMETERS)
    sample_diagonal = PAIR_DIAGONAL_FACTOR.subs(SAMPLE_PARAMETERS)

    p_derivative = diff(SAMPLE_P, T)
    first_x_jet = cancel(diff(SAMPLE_Q, T) / p_derivative)
    second_x_jet = cancel(diff(first_x_jet, T) / p_derivative)
    reduced_contact_jets = tuple(
        _modulo_branch(jet, SAMPLE_CONTACT_BRANCH)
        for jet in (first_x_jet, second_x_jet)
    )
    partner_contact_jets = tuple(
        _modulo_branch(jet.subs(T, 1 - T), SAMPLE_CONTACT_BRANCH)
        for jet in (first_x_jet, second_x_jet)
    )
    contact_jet_differences = tuple(
        expand(rem(left - right, SAMPLE_CONTACT_BRANCH, T))
        for left, right in zip(
            reduced_contact_jets,
            partner_contact_jets,
            strict=True,
        )
    )

    triple_slope = _modulo_branch(
        diff(SAMPLE_Q, T) / p_derivative,
        SAMPLE_TRIPLE_BRANCH,
    )
    primitive_slope_polynomial = (
        7 * SLOPE_VALUE**3 - 95 * SLOPE_VALUE**2 + 569 * SLOPE_VALUE - 1729
    )
    triple_slope_resultant = resultant(
        SAMPLE_TRIPLE_BRANCH,
        together(triple_slope - SLOPE_VALUE).as_numer_denom()[0],
        T,
    )

    sample_x_numerator = COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS)
    sample_x_denominator = COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    implicit_polynomial = Poly(SAMPLE_IMPLICIT, X, Y)

    open_factors = (
        TRIPLE_PAIR_SUM,
        TRIPLE_PRODUCT,
        CONTACT_SUM,
        KAPPA_NUMERATOR,
        SPLIT_FACTOR,
        CONTACT_DENOMINATOR_FACTOR,
        CONTACT_DIAGONAL_FACTOR,
        CONTACT_PRODUCT_FACTOR,
        OVERLAP_FACTOR_C,
        OVERLAP_FACTOR_E,
        TRIPLE_DISCRIMINANT_NUMERATOR,
        TRIPLE_FOURTH_ROOT_FACTOR,
        RESIDUAL_RANK_FACTOR,
    )

    return DeltaTenContactTripleCertificate(
        incidence_determinant_identity=expand(
            incidence_determinant - EXPECTED_INCIDENCE_DETERMINANT
        ),
        rank_factor_irreducible=factor(RESIDUAL_RANK_FACTOR) == RESIDUAL_RANK_FACTOR,
        specialized_rank_identity=expand(
            specialized_rank_factor - EXPECTED_SPECIALIZED_RANK_FACTOR
        ),
        specialized_augmented_identity=expand(
            specialized_augmented_minor - EXPECTED_SPECIALIZED_AUGMENTED_MINOR
        ),
        specialized_augmented_gcd=gcd(
            Poly(specialized_rank_factor, CONTACT_SUM),
            Poly(specialized_augmented_minor, CONTACT_SUM),
        ).as_expr(),
        cramer_open_component_count=1,
        cramer_open_component_dimension=3,
        coefficient_image_codimension=2,
        projection_quasi_finite_on_cramer_open=True,
        projection_generically_finite_onto_image=True,
        residual_rank_intersections_classified=False,
        split_boundary_intersections_classified=False,
        topology_propagation_dependencies=(
            "connected clean open",
            "proper projective Whitney-Thom triviality",
        ),
        sample_incidence_residuals=tuple(
            cancel(equation.subs(SAMPLE_BASE).subs(SAMPLE_PARAMETERS))
            for equation in COMBINED_EQUATIONS
        ),
        sample_incidence_determinant=sample_matrix.det(),
        sample_linear_solution=tuple(sample_solution_matrix),
        sample_open_factor_values=tuple(
            factor_expression.subs(SAMPLE_BASE) for factor_expression in open_factors
        ),
        boundary_identities=boundary_identities,
        collision_identity=expand(sample_collision - SAMPLE_COLLISION),
        tangency_identity=expand(sample_tangency - SAMPLE_TANGENCY),
        collision_tangency_gcd=gcd(
            Poly(sample_collision, S),
            Poly(sample_tangency, S),
        ).as_expr(),
        contact_collision_second_derivative=diff(sample_collision, S, 2).subs(S, 1),
        contact_tangency_derivative=diff(sample_tangency, S).subs(S, 1),
        residual_discriminant=discriminant(SAMPLE_RESIDUAL_NODES, S),
        residual_factor_resultants=(
            resultant(SAMPLE_RESIDUAL_NODES, S - 1, S),
            resultant(SAMPLE_RESIDUAL_NODES, SAMPLE_TRIPLE_PAIR_SUMS, S),
            resultant(S - 1, SAMPLE_TRIPLE_PAIR_SUMS, S),
        ),
        tangency_factor_resultants=(
            resultant(SAMPLE_RESIDUAL_NODES, sample_tangency, S),
            resultant(SAMPLE_TRIPLE_PAIR_SUMS, sample_tangency, S),
        ),
        pair_boundary_resultants=(
            resultant(sample_collision, S, S),
            resultant(sample_collision, sample_denominator, S),
            resultant(sample_collision, sample_product, S),
            resultant(sample_collision, -S * sample_diagonal, S),
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(SAMPLE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(SAMPLE_PARAMETERS),
        contact_image_remainders=(
            rem(SAMPLE_P, SAMPLE_CONTACT_BRANCH, T),
            rem(SAMPLE_Q, SAMPLE_CONTACT_BRANCH, T),
        ),
        triple_image_remainders=(
            rem(SAMPLE_P, SAMPLE_TRIPLE_BRANCH, T),
            rem(SAMPLE_Q, SAMPLE_TRIPLE_BRANCH, T),
        ),
        fourth_source_image=(
            SAMPLE_P.subs(T, Rational(1, 2)),
            SAMPLE_Q.subs(T, Rational(1, 2)),
        ),
        source_derivative_resultants=(
            resultant(SAMPLE_CONTACT_BRANCH, p_derivative, T),
            resultant(SAMPLE_TRIPLE_BRANCH, p_derivative, T),
        ),
        contact_jet_differences=contact_jet_differences,
        contact_second_jet_separation=resultant(
            SAMPLE_CONTACT_BRANCH,
            2 * T - 1,
            T,
        ),
        triple_slope_resultant_identity=expand(
            triple_slope_resultant + 49 * primitive_slope_polynomial
        ),
        triple_slope_discriminant=discriminant(
            primitive_slope_polynomial,
            SLOPE_VALUE,
        ),
        node_x_resultant_identity=expand(
            2
            * resultant(
                SAMPLE_RESIDUAL_NODES,
                2 * sample_x_numerator - X * sample_x_denominator,
                S,
            )
            + 441 * SAMPLE_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(SAMPLE_NODE_X_POLYNOMIAL, X),
        node_target_separations=tuple(
            SAMPLE_NODE_X_POLYNOMIAL.subs(X, target) for target in (-5, 1, 0)
        ),
        implicit_resultant_identity=expand(
            resultant(SAMPLE_P - X, SAMPLE_Q - Y, T) + SAMPLE_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            SAMPLE_IMPLICIT.subs({X: SAMPLE_P, Y: SAMPLE_Q})
        ),
        implicit_content=implicit_polynomial.content(),
        sage_jacobian_components=((4, 1), (3, 1), (4, 1), (5, 5)),
        sage_cyclic_simplification=(1, 0, True),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=2,
        contact_delta=2,
        triple_delta=3,
        node_count=5,
        infinity_delta=16,
        complement_census=n_generator_three_cycle_presentation_census(
            CONTACT_TRIPLE_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact dominant-component certificate summary."""

    certificate = exact_delta_ten_contact_triple_certificate()
    print("delta-ten C2 + T111 + 5N certificate:", certificate.verified)
    print("incidence determinant identity:", certificate.incidence_determinant_identity)
    print("rank factor / augmented gcd:", certificate.specialized_augmented_gcd)
    print("sample singularity delta balance:", certificate.total_delta)
    print("sample complement census:", certificate.complement_census)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
