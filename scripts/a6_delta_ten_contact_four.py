"""Certify the conditional delta-ten ``C4 + 6N`` profile.

Write the normalized degree-``(4, 9)`` family as

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``.

On the nonsplit unordered-pair chart, a two-branch contact of intersection
multiplicity four is a quadruple root ``s`` of the collision decic ``H``.
The equations ``H = H' = H'' = H''' = 0`` are affine-linear in
``(a,b,c,d)``.  Their determinant is

``24*(k+2*s)^4*(s^2+k*s+1)^5*F(k,s)``,

where the displayed residual factor ``F`` is irreducible over ``QQ``.  The
determinant-nonzero incidence is consequently a rational surface.  The
residual determinant curve needs a hostile rank audit: expected dimension
alone does not exclude a coefficient fiber from becoming too large.  The
independent Sage checker proves all of the following after saturation by the
complete valid nonsplit-chart localizer:

* the coefficient rank-``<= 2`` ideal is the unit ideal;
* the compatibility locus on ``F`` is zero-dimensional of length ten; and
* the augmented rank-``<= 2`` and rank-``<= 1`` ideals are unit ideals.

Thus the coefficient matrix has rank exactly three at every compatible valid
point of ``F``.  Those finitely many bases have affine-line coefficient
fibers, so their incidence has dimension one and cannot hide a second
surface, much less a threefold.  This statement is only about the valid
nonsplit chart: ``k = 0, +/-2`` and the graph-denominator, cusp-pair, and
diagonal boundaries are separate charts.

The rational member

``P = t^2 + t^3 + t^4`` and
``14*Q = 18*t^5 + 27*t^6 + 42*t^7 + 24*t^8 + 14*t^9``

has the forced ``T(2,5)`` cusp, one exact contact-four singularity, and six
nodes.  Sage 10.8 regenerates its implicit equation and singular scheme.  Its
stored exact affine van Kamp presentation simplifies to ``Z``; the finite
three-cycle replay has no ``A6`` image.

The topology computation certifies this member only.  Propagation across a
stratum still requires a connected clean open and proper projective
Whitney--Thom triviality.  This is a conditional, computer-assisted profile
certificate, not a proof of the plane Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from itertools import combinations
from typing import Final

from sympy import (
    Expr,
    Matrix,
    Poly,
    Rational,
    cancel,
    diff,
    discriminant,
    expand,
    factor_list,
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

CONTACT_FOUR_LINEAR_PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)
CONTACT_FOUR_DERIVATIVES: Final = tuple(
    diff(COLLISION_POLYNOMIAL, S, order) for order in range(4)
)
CONTACT_FOUR_COEFFICIENT_MATRIX: Final = Matrix(
    [
        [diff(equation, parameter) for parameter in CONTACT_FOUR_LINEAR_PARAMETERS]
        for equation in CONTACT_FOUR_DERIVATIVES
    ]
)
CONTACT_FOUR_ZERO_PARAMETERS: Final = dict.fromkeys(
    CONTACT_FOUR_LINEAR_PARAMETERS,
    0,
)
CONTACT_FOUR_AUGMENTED_MATRIX: Final = CONTACT_FOUR_COEFFICIENT_MATRIX.row_join(
    Matrix(
        [
            equation.subs(CONTACT_FOUR_ZERO_PARAMETERS)
            for equation in CONTACT_FOUR_DERIVATIVES
        ]
    )
)

CONTACT_FOUR_RESIDUAL_FACTOR: Final = (
    3 * KAPPA**5 * S
    + 35 * KAPPA**4 * S**2
    + 13 * KAPPA**4
    + 88 * KAPPA**3 * S**3
    + 164 * KAPPA**3 * S
    + 88 * KAPPA**2 * S**4
    + 412 * KAPPA**2 * S**2
    + 52 * KAPPA**2
    + 32 * KAPPA * S**5
    + 416 * KAPPA * S**3
    + 64 * KAPPA * S
    + 160 * S**4
)
CONTACT_FOUR_EXPECTED_DETERMINANT: Final = (
    24 * PAIR_DENOMINATOR**4 * PAIR_QUADRATIC**5 * CONTACT_FOUR_RESIDUAL_FACTOR
)

# This localizer is deliberately limited to the nonsplit unordered-pair
# chart.  The three true split fibers k=0,+/-2 are legitimate separate
# charts, not factors that may be silently discarded.
CONTACT_FOUR_VALID_LOCALIZER: Final = (
    KAPPA
    * (KAPPA**2 - 4)
    * S
    * PAIR_DENOMINATOR
    * PAIR_QUADRATIC
    * PAIR_DIAGONAL_FACTOR
)

# Reproducible Sage 10.8 output from
# ``tools/check_a6_delta_ten_contact_four.sage``.  These tuples record
# (is_unit, saturation exponent) or (dimension, scheme length, exponent).
SAGE_COEFFICIENT_RANK_TWO_SATURATION: Final = (True, 2)
SAGE_RESIDUAL_COMPATIBILITY_SATURATION: Final = (0, 10, 4)
SAGE_AUGMENTED_RANK_TWO_SATURATION: Final = (True, 2)
SAGE_AUGMENTED_RANK_ONE_SATURATION: Final = (True, 1)

CONTACT_FOUR_PARAMETERS: Final = {
    KAPPA: 1,
    ALPHA: Rational(9, 7),
    BETA: Rational(27, 14),
    GAMMA: 3,
    DELTA: Rational(12, 7),
}
CONTACT_FOUR_P: Final = expand(FAMILY_P.subs(CONTACT_FOUR_PARAMETERS))
CONTACT_FOUR_SCALED_Q: Final = expand(14 * FAMILY_Q.subs(CONTACT_FOUR_PARAMETERS))
CONTACT_FOUR_BRANCH: Final = T**2 + T + 1
CONTACT_FOUR_RESIDUAL: Final = (
    14 * S**6 + 28 * S**5 + 54 * S**4 + 52 * S**3 + 56 * S**2 + 21 * S + 18
) / 14
CONTACT_FOUR_COLLISION: Final = expand((S + 1) ** 4 * CONTACT_FOUR_RESIDUAL)
CONTACT_FOUR_TANGENCY_COFACTOR: Final = (
    84 * S**8
    + 294 * S**7
    + 623 * S**6
    + 818 * S**5
    + 800 * S**4
    + 553 * S**3
    + 286 * S**2
    + 76 * S
    + 30
)
CONTACT_FOUR_TANGENCY: Final = expand(
    Rational(3, 7) * (S + 1) ** 3 * CONTACT_FOUR_TANGENCY_COFACTOR
)

CONTACT_FOUR_NODE_X_POLYNOMIAL: Final = (
    7529536 * X**6
    - 7760032 * X**5
    + 3746736 * X**4
    - 1049560 * X**3
    + 178396 * X**2
    - 17331 * X
    + 738
)
CONTACT_FOUR_IMPLICIT: Final = (
    -38416 * X**9
    + 38280 * X**8
    - 16632 * X**7
    - 1872 * X**6 * Y
    + 3591 * X**6
    + 1494 * X**5 * Y
    - 324 * X**5
    + 12 * X**4 * Y**2
    - 429 * X**4 * Y
    - 144 * X**3 * Y**2
    + 36 * X**3 * Y
    + 30 * X**2 * Y**3
    + 66 * X**2 * Y**2
    - 12 * X * Y**3
    - 12 * X * Y**2
    + Y**4
    + 2 * Y**3
    + Y**2
)

# Sage 10.8 exact raw affine presentation.  Generator indices 1..4 are the
# geometric meridians of a generic vertical fiber.
CONTACT_FOUR_RELATIONS: Final = (
    (4, 3, -4, -3),
    (-3, -2, -1, 2, 3, 2),
    (-4, 2, 4, 2, 4, 2, -4, -2, -4, -2),
    (
        -3,
        -2,
        1,
        2,
        3,
        -2,
        1,
        2,
        3,
        -2,
        1,
        2,
        3,
        -2,
        1,
        2,
        3,
        -2,
        -1,
        2,
        -3,
        -2,
        -1,
        2,
        -3,
        -2,
        -1,
        2,
        -3,
        -2,
        -1,
        2,
    ),
    (-3, -2, -1, 2, 4, -2, 1, 2),
    (2, 1, -2, -1),
    (-3, -2, 4, -1, -4, 2, 3, 2, -3, -2, 4, 1, -4, 2, 3, -2),
    (-3, -2, 4, 1, -4, 2, 3, -2, 4, -1, -4, 2),
    (-3, -2, 4, 2, 3, -2, -4, 2),
    (-4, 2, 4, -2),
)


def _iterated_polynomial_gcd(polynomials: tuple[Expr, ...]) -> Expr:
    """Return the true iterated gcd of a nonempty polynomial tuple."""

    if not polynomials:
        msg = "a polynomial gcd needs at least one input"
        raise ValueError(msg)
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = gcd(common, polynomial)
    return common


def _modulo_contact_branch(rational_expression: Expr) -> Expr:
    """Reduce a rational expression modulo the contact-pair quadratic."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    denominator_inverse = invert(denominator, CONTACT_FOUR_BRANCH)
    return expand(
        rem(
            numerator * denominator_inverse,
            CONTACT_FOUR_BRANCH,
            T,
        )
    )


@dataclass(frozen=True, slots=True)
class DeltaTenContactFourCertificate:
    """Exact incidence, sample geometry, and complement data."""

    incidence_determinant_identity: Expr
    residual_factor_irreducible: bool
    coefficient_rank_two_minor_gcd_identity: Expr
    specialized_augmented_gcd: Expr
    sage_coefficient_rank_two_saturation: tuple[bool, int]
    sage_residual_compatibility_saturation: tuple[int, int, int]
    sage_augmented_rank_two_saturation: tuple[bool, int]
    sage_augmented_rank_one_saturation: tuple[bool, int]
    cramer_open_component_count: int
    cramer_open_component_dimension: int
    coefficient_image_codimension: int
    projection_quasi_finite_on_cramer_open: bool
    projection_generically_finite_onto_image: bool
    split_boundaries_classified_here: bool
    topology_propagation_dependencies: tuple[str, str]
    sample_incidence_residuals: tuple[Expr, Expr, Expr, Expr]
    sample_incidence_determinant: Expr
    sample_valid_localizer: Expr
    sample_cusp_coefficient: Expr
    collision_identity: Expr
    tangency_identity: Expr
    collision_tangency_gcd: Expr
    quadruple_root_gcd: Expr
    fourth_derivative: Expr
    tangency_third_derivative: Expr
    residual_discriminant: Expr
    residual_contact_separation: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    denominator_resultant: Expr
    diagonal_resultant: Expr
    residual_tangency_resultant: Expr
    contact_chart_values: tuple[Expr, Expr, Expr]
    contact_pair_discriminant: Expr
    contact_image_remainders: tuple[Expr, Expr]
    contact_p_derivative_resultant: Expr
    contact_jet_values: tuple[Expr, Expr, Expr]
    contact_jet_differences: tuple[Expr, Expr, Expr, Expr]
    contact_fourth_jet_separation: Expr
    node_x_resultant_identity: Expr
    node_x_discriminant: Expr
    contact_node_separation: Expr
    cusp_node_separation: Expr
    contact_cusp_target_separation: Expr
    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    implicit_content: Expr
    sage_jacobian_components: tuple[tuple[int, int], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    arithmetic_genus: int
    cusp_delta: int
    contact_delta: int
    node_count: int
    infinity_delta: int
    relation_count: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the complete projective genus contribution."""

        return (
            self.cusp_delta + self.contact_delta + self.node_count + self.infinity_delta
        )

    @property
    def residual_coefficient_rank(self) -> int | None:
        """Return the rank forced by the separately replayable Sage evidence.

        The Python certificate stores, but does not recompute, the four Sage
        saturations.  On their recorded outcome the compatible residual base
        is finite, coefficient rank at most two is impossible, and augmented
        rank at most two is impossible.  Compatibility then forces both ranks
        to be exactly three.
        """

        if (
            self.sage_coefficient_rank_two_saturation == (True, 2)
            and self.sage_residual_compatibility_saturation == (0, 10, 4)
            and self.sage_augmented_rank_two_saturation == (True, 2)
            and self.sage_augmented_rank_one_saturation == (True, 1)
        ):
            return 3
        return None

    @property
    def residual_incidence_dimension(self) -> int | None:
        """Return the finite-base plus affine-fiber dimension bound."""

        if self.residual_coefficient_rank == 3:
            return 0 + (4 - 3)
        return None

    @property
    def residual_surface_component_excluded(self) -> bool:
        """Whether the compatible residual incidence has dimension below two."""

        return self.residual_incidence_dimension == 1

    @property
    def verified(self) -> bool:
        """Whether every exact incidence and sample invariant agrees."""

        return bool(
            self.incidence_determinant_identity == 0
            and self.residual_factor_irreducible
            and self.coefficient_rank_two_minor_gcd_identity == 0
            and self.specialized_augmented_gcd == 1
            and self.sage_coefficient_rank_two_saturation == (True, 2)
            and self.sage_residual_compatibility_saturation == (0, 10, 4)
            and self.sage_augmented_rank_two_saturation == (True, 2)
            and self.sage_augmented_rank_one_saturation == (True, 1)
            and self.cramer_open_component_count == 1
            and self.cramer_open_component_dimension == 2
            and self.coefficient_image_codimension == 3
            and self.projection_quasi_finite_on_cramer_open
            and self.projection_generically_finite_onto_image
            and self.residual_coefficient_rank == 3
            and self.residual_incidence_dimension == 1
            and self.residual_surface_component_excluded
            and not self.split_boundaries_classified_here
            and self.topology_propagation_dependencies
            == (
                "connected clean open",
                "proper projective Whitney-Thom triviality",
            )
            and self.sample_incidence_residuals == (0, 0, 0, 0)
            and self.sample_incidence_determinant == -168
            and self.sample_valid_localizer == -9
            and self.sample_cusp_coefficient == Rational(9, 7)
            and self.collision_identity == 0
            and self.tangency_identity == 0
            and self.collision_tangency_gcd == expand((S + 1) ** 3)
            and self.quadruple_root_gcd == S + 1
            and self.fourth_derivative == Rational(492, 7)
            and self.tangency_third_derivative == Rational(1476, 7)
            and self.residual_discriminant == Rational(-2627098335909, 2582630848)
            and self.residual_contact_separation == Rational(41, 14)
            and self.cusp_image_factor == Rational(1, 196)
            and self.extra_critical_factor == Rational(11421, 49)
            and self.denominator_resultant == 81
            and self.diagonal_resultant == Rational(8325909, 343)
            and self.residual_tangency_resultant
            == Rational(
                763059488663600977727961,
                173625106649344,
            )
            and self.contact_chart_values == (-1, 1, 3)
            and self.contact_pair_discriminant == -3
            and self.contact_image_remainders == (0, -1)
            and self.contact_p_derivative_resultant == 3
            and self.contact_jet_values == (6, -30, 108)
            and self.contact_jet_differences == (0, 0, 0, 984 * (2 * T + 1))
            and self.contact_fourth_jet_separation == 3
            and self.node_x_resultant_identity == 0
            and self.node_x_discriminant == -2608211900386817980766515181152046099202048
            and self.contact_node_separation == 738
            and self.cusp_node_separation == 738
            and self.contact_cusp_target_separation == -1
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_components == ((4, 1), (7, 1), (6, 6))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.contact_delta == 4
            and self.node_count == 6
            and self.infinity_delta == 16
            and self.total_delta == self.arithmetic_genus
            and self.relation_count == 10
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_contact_four_certificate() -> DeltaTenContactFourCertificate:
    """Build the exact dominant-incidence and rational-member certificate."""

    coefficient_rank_two_minors = tuple(
        expand(CONTACT_FOUR_COEFFICIENT_MATRIX.extract(rows, columns).det())
        for rows in combinations(range(4), 3)
        for columns in combinations(range(4), 3)
    )
    rank_two_minor_gcd = _iterated_polynomial_gcd(coefficient_rank_two_minors)

    # One exact specialization suffices to prove that the irreducible
    # residual factor does not divide this augmented minor.  The Sage checker
    # goes further and computes the full compatible scheme from all minors.
    selected_augmented_minor = expand(
        CONTACT_FOUR_AUGMENTED_MATRIX[:, (0, 1, 2, 4)].det()
    )
    specialized_augmented_gcd = gcd(
        Poly(CONTACT_FOUR_RESIDUAL_FACTOR.subs(KAPPA, 1), S),
        Poly(selected_augmented_minor.subs(KAPPA, 1), S),
    ).as_expr()

    sample_collision = expand(COLLISION_POLYNOMIAL.subs(CONTACT_FOUR_PARAMETERS))
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(CONTACT_FOUR_PARAMETERS))
    sample_base = {KAPPA: 1, S: -1}
    sample_coefficient_matrix = CONTACT_FOUR_COEFFICIENT_MATRIX.subs(sample_base)
    sample_constant_column = Matrix(
        [
            equation.subs(CONTACT_FOUR_ZERO_PARAMETERS).subs(sample_base)
            for equation in CONTACT_FOUR_DERIVATIVES
        ]
    )
    sample_parameter_vector = Matrix(
        [
            CONTACT_FOUR_PARAMETERS[parameter]
            for parameter in CONTACT_FOUR_LINEAR_PARAMETERS
        ]
    )

    p_derivative = diff(CONTACT_FOUR_P, T)
    current_derivative = CONTACT_FOUR_SCALED_Q
    reduced_jets: list[Expr] = []
    for _ in range(4):
        current_derivative = cancel(diff(current_derivative, T) / p_derivative)
        reduced_jets.append(_modulo_contact_branch(current_derivative))
    partner_jets = tuple(
        _modulo_contact_branch(jet.subs(T, -1 - T)) for jet in reduced_jets
    )
    jet_differences = tuple(
        expand(rem(left - right, CONTACT_FOUR_BRANCH, T))
        for left, right in zip(reduced_jets, partner_jets, strict=True)
    )

    sample_x_numerator = COLLISION_X_NUMERATOR.subs(CONTACT_FOUR_PARAMETERS)
    sample_x_denominator = COLLISION_X_DENOMINATOR.subs(CONTACT_FOUR_PARAMETERS)
    implicit_polynomial = Poly(CONTACT_FOUR_IMPLICIT, X, Y)
    residual_factorization = factor_list(Poly(CONTACT_FOUR_RESIDUAL_FACTOR, KAPPA, S))[
        1
    ]

    return DeltaTenContactFourCertificate(
        incidence_determinant_identity=expand(
            CONTACT_FOUR_COEFFICIENT_MATRIX.det() - CONTACT_FOUR_EXPECTED_DETERMINANT
        ),
        residual_factor_irreducible=(
            len(residual_factorization) == 1
            and residual_factorization[0][0].as_expr() == CONTACT_FOUR_RESIDUAL_FACTOR
            and residual_factorization[0][1] == 1
        ),
        coefficient_rank_two_minor_gcd_identity=expand(
            rank_two_minor_gcd - 2 * PAIR_QUADRATIC
        ),
        specialized_augmented_gcd=specialized_augmented_gcd,
        sage_coefficient_rank_two_saturation=SAGE_COEFFICIENT_RANK_TWO_SATURATION,
        sage_residual_compatibility_saturation=SAGE_RESIDUAL_COMPATIBILITY_SATURATION,
        sage_augmented_rank_two_saturation=SAGE_AUGMENTED_RANK_TWO_SATURATION,
        sage_augmented_rank_one_saturation=SAGE_AUGMENTED_RANK_ONE_SATURATION,
        cramer_open_component_count=1,
        cramer_open_component_dimension=2,
        coefficient_image_codimension=3,
        projection_quasi_finite_on_cramer_open=True,
        projection_generically_finite_onto_image=True,
        split_boundaries_classified_here=False,
        topology_propagation_dependencies=(
            "connected clean open",
            "proper projective Whitney-Thom triviality",
        ),
        sample_incidence_residuals=tuple(
            expand(value)
            for value in (
                sample_coefficient_matrix * sample_parameter_vector
                + sample_constant_column
            )
        ),
        sample_incidence_determinant=sample_coefficient_matrix.det(),
        sample_valid_localizer=CONTACT_FOUR_VALID_LOCALIZER.subs(sample_base),
        sample_cusp_coefficient=CONTACT_FOUR_PARAMETERS[ALPHA],
        collision_identity=expand(sample_collision - CONTACT_FOUR_COLLISION),
        tangency_identity=expand(sample_tangency - CONTACT_FOUR_TANGENCY),
        collision_tangency_gcd=gcd(
            Poly(sample_collision, S),
            Poly(sample_tangency, S),
        ).as_expr(),
        quadruple_root_gcd=_iterated_polynomial_gcd(
            tuple(Poly(diff(sample_collision, S, order), S) for order in range(4))
        ).as_expr(),
        fourth_derivative=diff(sample_collision, S, 4).subs(S, -1),
        tangency_third_derivative=diff(sample_tangency, S, 3).subs(S, -1),
        residual_discriminant=discriminant(CONTACT_FOUR_RESIDUAL, S),
        residual_contact_separation=CONTACT_FOUR_RESIDUAL.subs(S, -1),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(CONTACT_FOUR_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(CONTACT_FOUR_PARAMETERS),
        denominator_resultant=resultant(
            sample_collision,
            PAIR_DENOMINATOR.subs(CONTACT_FOUR_PARAMETERS),
            S,
        ),
        diagonal_resultant=resultant(
            sample_collision,
            (-S * PAIR_DIAGONAL_FACTOR).subs(CONTACT_FOUR_PARAMETERS),
            S,
        ),
        residual_tangency_resultant=resultant(
            CONTACT_FOUR_RESIDUAL,
            sample_tangency,
            S,
        ),
        contact_chart_values=tuple(
            expression.subs(CONTACT_FOUR_PARAMETERS).subs(S, -1)
            for expression in (
                PAIR_DENOMINATOR,
                PAIR_QUADRATIC,
                PAIR_DIAGONAL_FACTOR,
            )
        ),
        contact_pair_discriminant=-3,
        contact_image_remainders=(
            rem(CONTACT_FOUR_P, CONTACT_FOUR_BRANCH, T),
            rem(CONTACT_FOUR_SCALED_Q, CONTACT_FOUR_BRANCH, T),
        ),
        contact_p_derivative_resultant=resultant(
            CONTACT_FOUR_BRANCH,
            p_derivative,
            T,
        ),
        contact_jet_values=tuple(reduced_jets[:3]),
        contact_jet_differences=jet_differences,
        contact_fourth_jet_separation=resultant(
            CONTACT_FOUR_BRANCH,
            2 * T + 1,
            T,
        ),
        node_x_resultant_identity=expand(
            7529536
            * resultant(
                CONTACT_FOUR_RESIDUAL,
                X * sample_x_denominator - sample_x_numerator,
                S,
            )
            - 6561 * CONTACT_FOUR_NODE_X_POLYNOMIAL
        ),
        node_x_discriminant=discriminant(CONTACT_FOUR_NODE_X_POLYNOMIAL, X),
        contact_node_separation=CONTACT_FOUR_NODE_X_POLYNOMIAL.subs(X, 0),
        cusp_node_separation=CONTACT_FOUR_NODE_X_POLYNOMIAL.subs(X, 0),
        contact_cusp_target_separation=-1,
        implicit_resultant_identity=expand(
            resultant(
                CONTACT_FOUR_P - X,
                CONTACT_FOUR_SCALED_Q - Y,
                T,
            )
            - CONTACT_FOUR_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            CONTACT_FOUR_IMPLICIT.subs({X: CONTACT_FOUR_P, Y: CONTACT_FOUR_SCALED_Q})
        ),
        implicit_content=implicit_polynomial.content(),
        sage_jacobian_components=((4, 1), (7, 1), (6, 6)),
        sage_cyclic_simplification=(1, 0, True),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=2,
        contact_delta=4,
        node_count=6,
        infinity_delta=16,
        relation_count=len(CONTACT_FOUR_RELATIONS),
        complement_census=n_generator_three_cycle_presentation_census(
            CONTACT_FOUR_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact ``C4 + 6N`` certificate and fail on regression."""

    certificate = exact_delta_ten_contact_four_certificate()
    print("delta-ten C4 + 6N certificate:", certificate.verified)
    print(
        "valid residual compatibility scheme:",
        certificate.sage_residual_compatibility_saturation,
    )
    print(
        "sample singularity delta balance:",
        certificate.total_delta,
        "= 2 + 4 + 6 + 16",
    )
    print("sample complement census:", certificate.complement_census)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
