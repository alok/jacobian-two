"""Exact algebra for the conditional delta-ten two-contact stratum.

The expected-codimension-two profile ``C2^2 + 6N`` has two distinct
two-branch contact-two target fibers and six ordinary nodes.  On the generic
pair-incidence chart, introduce ordered prospective contact roots ``u`` and
``v`` of the collision decic ``H`` and impose

``H(u) = H'(u) = H(v) = H'(v) = 0``.

These four equations are affine-linear in ``(a,b,c,d)``.  Their determinant
defines one nonempty rational ordered threefold on its nonvanishing open.  The
swap ``u <-> v`` is free after ``u-v`` is inverted, and the unique Cramer
solution is swap-invariant, so it descends to the unordered coordinates
``u+v`` and ``u*v``.

The determinant-boundary audit is deliberately precise about its limits:

* ``u=v`` repeats one contact equation and is not ``C2^2``;
* ``u^2+k*u+1=0`` or its ``v`` analogue is the cusp-image boundary on the
  nonsplit chart and is invalid under the standing hypotheses;
* ``2*u+k=0`` or ``2*v+k=0`` forces ``k=0,+/-2`` when it is a collision and
  belongs to the split incidence charts; and
* ``s*(2*s^2+3*k*s+4)=0`` is the coincident-source/critical boundary.  Away
  from the split values its collision resultant is the product of ``a`` and
  the extra-critical resultant, so both factors must be inverted on the
  valid two-branch contact chart; and
* the residual rank factor is generically incompatible with the augmented
  linear system, so it does not supply another dominant Cramer component.

For each split value ``k=0,+/-2``, the two true component discriminants are
nonzero, squarefree, and coprime.  This rules out a common or nonreduced
discriminant divisor that would automatically provide two contacts.  It does
not classify all codimension-two intersections inside the split fibers or
all lower-dimensional closures of the residual-rank boundary; those remain
explicit elimination gaps.

The rational member

``(k,a,b,c,d) = (1,3,3/2,3,0)``

has repeated collision roots ``u=-1`` and ``v=1``.  Exact resultants certify
two contact-two fibers, six transverse nodes, distinct collision targets,
the forced cusp, and the fixed infinity branch.  Sage 10.8 regenerates its
primitive implicit equation and singular scheme, and its exact eleven-relator
affine van Kamp presentation simplifies to an infinite cyclic group with an
explicit section.  An independent replay exhausts all ``40^4`` assignments
of the raw meridians to single three-cycles; only forty order-three images
survive, so no assignment generates ``A6``.

Thus the generic clean open of the component dominating the full ordered
Cramer base is excluded, conditional on the standing reduction and the
standard proper Whitney--Thom propagation.  This does not rule out components
supported on lower-dimensional rank-drop or split loci.  Those closures remain
open, and nothing here proves the plane Jacobian conjecture.
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
    resultant,
)
from sympy.polys.matrices import DomainMatrix

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
    R,
    S,
    T,
    TANGENCY_COFACTOR,
    TANGENCY_POLYNOMIAL,
    X,
    Y,
    ZERO_FIBER_GRAPH,
    ZERO_FIBER_VERTICAL,
    exceptional_graph,
    exceptional_vertical,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    n_generator_three_cycle_presentation_census,
)

U: Final = Symbol("u")
V: Final = Symbol("v")
ROOT_SUM: Final = Symbol("sigma")
ROOT_PRODUCT: Final = Symbol("pi")

PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)
ZERO_PARAMETERS: Final = {parameter: 0 for parameter in PARAMETERS}


def _domain_determinant(matrix: Matrix) -> Expr:
    """Compute a polynomial determinant without heuristic simplification."""

    determinant: Expr = DomainMatrix.from_Matrix(matrix).det().as_expr()
    return determinant


@cache
def ordered_incidence_system() -> tuple[tuple[Expr, ...], Matrix, Matrix]:
    """Return the four ordered equations, coefficient matrix, and right side."""

    equations = (
        expand(COLLISION_POLYNOMIAL.subs(S, U)),
        expand(diff(COLLISION_POLYNOMIAL, S).subs(S, U)),
        expand(COLLISION_POLYNOMIAL.subs(S, V)),
        expand(diff(COLLISION_POLYNOMIAL, S).subs(S, V)),
    )
    coefficient_matrix = Matrix(
        [
            [equation.coeff(parameter) for parameter in PARAMETERS]
            for equation in equations
        ]
    )
    right_side = Matrix([-equation.subs(ZERO_PARAMETERS) for equation in equations])
    return equations, coefficient_matrix, right_side


@cache
def ordered_determinant_and_cramer_numerators() -> tuple[Expr, tuple[Expr, ...]]:
    """Return the ordered determinant and four exact Cramer numerators."""

    _, coefficient_matrix, right_side = ordered_incidence_system()
    determinant = _domain_determinant(coefficient_matrix)
    numerators: list[Expr] = []
    for column in range(len(PARAMETERS)):
        replaced = coefficient_matrix.copy()
        replaced[:, column] = right_side
        numerators.append(_domain_determinant(replaced))
    return determinant, tuple(numerators)


ORDERED_VISIBLE_DETERMINANT_FACTOR: Final = (
    (KAPPA + 2 * U) ** 2
    * (KAPPA + 2 * V) ** 2
    * (U - V) ** 4
    * (U**2 + KAPPA * U + 1)
    * (V**2 + KAPPA * V + 1)
)


@cache
def ordered_residual_rank_factor() -> tuple[Expr, Expr]:
    """Divide the determinant by its geometric visible factors exactly."""

    determinant, _ = ordered_determinant_and_cramer_numerators()
    quotient, remainder = Poly(
        determinant,
        KAPPA,
        U,
        V,
    ).div(
        Poly(
            ORDERED_VISIBLE_DETERMINANT_FACTOR,
            KAPPA,
            U,
            V,
        )
    )
    return quotient.as_expr(), remainder.as_expr()


@cache
def ordered_compatibility_gcd() -> Expr:
    """Common divisor of the determinant and all augmented Cramer minors."""

    determinant, numerators = ordered_determinant_and_cramer_numerators()
    common = Poly(determinant, KAPPA, U, V)
    for numerator in numerators:
        common = common.gcd(Poly(numerator, KAPPA, U, V))
    result: Expr = common.as_expr()
    return factor(result)


ORDERED_EXPECTED_COMPATIBILITY_GCD: Final = (
    (U - V) ** 4 * (U**2 + KAPPA * U + 1) * (V**2 + KAPPA * V + 1)
)


def _swap_roots(expression: Expr) -> Expr:
    """Exchange the two ordered repeated roots simultaneously."""

    return expression.xreplace({U: V, V: U})


@cache
def dominant_parameterization() -> dict[Symbol, Expr]:
    """Return the unique Cramer solution on the determinant nonvanishing open."""

    determinant, numerators = ordered_determinant_and_cramer_numerators()
    return {
        parameter: numerator / determinant
        for parameter, numerator in zip(PARAMETERS, numerators, strict=True)
    }


@dataclass(frozen=True, slots=True)
class OrderedIncidenceCertificate:
    """Exact dominant-chart, quotient, and determinant-boundary checks."""

    determinant_factor_remainder: Expr
    residual_rank_swap_identity: Expr
    cramer_swap_identities: tuple[Expr, ...]
    compatibility_gcd_identity: Expr
    denominator_root_identity: Expr
    cusp_boundary_resultant_identity: Expr
    critical_boundary_resultant_identity: Expr
    tangency_syzygy_identity: Expr
    sample_determinant: Expr
    sample_parameter_remainders: tuple[Expr, ...]

    @property
    def verified(self) -> bool:
        """Whether the ordered dominant chart and visible boundaries agree."""

        return bool(
            self.determinant_factor_remainder == 0
            and self.residual_rank_swap_identity == 0
            and all(identity == 0 for identity in self.cramer_swap_identities)
            and self.compatibility_gcd_identity == 0
            and self.denominator_root_identity == 0
            and self.cusp_boundary_resultant_identity == 0
            and self.critical_boundary_resultant_identity == 0
            and self.tangency_syzygy_identity == 0
            and self.sample_determinant != 0
            and self.sample_parameter_remainders == (0, 0, 0, 0)
        )


@cache
def exact_ordered_incidence_certificate() -> OrderedIncidenceCertificate:
    """Build the exact generic ordered-incidence certificate."""

    determinant, numerators = ordered_determinant_and_cramer_numerators()
    residual_rank_factor, factor_remainder = ordered_residual_rank_factor()
    parameterization = dominant_parameterization()
    sample_root_substitution = {KAPPA: 1, U: -1, V: 1}
    sample_parameters = (3, Rational(3, 2), 3, 0)
    return OrderedIncidenceCertificate(
        determinant_factor_remainder=factor_remainder,
        residual_rank_swap_identity=expand(
            _swap_roots(residual_rank_factor) - residual_rank_factor
        ),
        cramer_swap_identities=tuple(
            expand(_swap_roots(numerator) - numerator) for numerator in numerators
        ),
        compatibility_gcd_identity=expand(
            ordered_compatibility_gcd() - ORDERED_EXPECTED_COMPATIBILITY_GCD
        ),
        denominator_root_identity=expand(
            COLLISION_POLYNOMIAL.subs(S, -KAPPA / 2)
            - KAPPA**2 * (KAPPA**2 - 4) ** 4 / 1024
        ),
        cusp_boundary_resultant_identity=expand(
            resultant(COLLISION_POLYNOMIAL, PAIR_QUADRATIC, S)
            - (KAPPA**2 - 4) ** 4 * CUSP_IMAGE_FACTOR
        ),
        critical_boundary_resultant_identity=expand(
            resultant(
                COLLISION_POLYNOMIAL,
                -S * PAIR_DIAGONAL_FACTOR,
                S,
            )
            - ALPHA * KAPPA**2 * (KAPPA**2 - 4) ** 4 * EXTRA_CRITICAL_FACTOR
        ),
        tangency_syzygy_identity=expand(
            PAIR_QUADRATIC * TANGENCY_POLYNOMIAL
            - S
            * PAIR_DIAGONAL_FACTOR
            * PAIR_DENOMINATOR
            * diff(COLLISION_POLYNOMIAL, S)
            - TANGENCY_COFACTOR * COLLISION_POLYNOMIAL
        ),
        sample_determinant=determinant.subs(sample_root_substitution),
        sample_parameter_remainders=tuple(
            cancel(
                parameterization[parameter].subs(sample_root_substitution) - expected
            )
            for parameter, expected in zip(
                PARAMETERS,
                sample_parameters,
                strict=True,
            )
        ),
    )


@dataclass(frozen=True, slots=True)
class SplitBoundaryCertificate:
    """Discriminant checks for one reducible pair-incidence fiber."""

    kappa: int
    vertical_discriminant_degree: int
    graph_discriminant_degree: int
    vertical_discriminant_squarefree: bool
    graph_discriminant_squarefree: bool
    common_discriminant_degree: int

    @property
    def verified(self) -> bool:
        """Whether neither split component has an automatic double-contact divisor."""

        expected_degrees = (2, 13) if self.kappa == 0 else (6, 9)
        return bool(
            (
                self.vertical_discriminant_degree,
                self.graph_discriminant_degree,
            )
            == expected_degrees
            and self.vertical_discriminant_squarefree
            and self.graph_discriminant_squarefree
            and self.common_discriminant_degree == 0
        )


@cache
def exact_split_boundary_certificates() -> tuple[SplitBoundaryCertificate, ...]:
    """Check squarefreeness and coprimality at ``k=0,+/-2``."""

    split_data = (
        (0, ZERO_FIBER_VERTICAL, ZERO_FIBER_GRAPH),
        (2, exceptional_vertical(1), exceptional_graph(1)),
        (-2, exceptional_vertical(-1), exceptional_graph(-1)),
    )
    result: list[SplitBoundaryCertificate] = []
    for kappa, vertical, graph in split_data:
        vertical_discriminant = Poly(
            discriminant(vertical, R),
            *PARAMETERS,
        )
        graph_discriminant = Poly(
            discriminant(graph, S),
            *PARAMETERS,
        )
        result.append(
            SplitBoundaryCertificate(
                kappa=kappa,
                vertical_discriminant_degree=vertical_discriminant.total_degree(),
                graph_discriminant_degree=graph_discriminant.total_degree(),
                vertical_discriminant_squarefree=(
                    vertical_discriminant.sqf_part().total_degree()
                    == vertical_discriminant.total_degree()
                ),
                graph_discriminant_squarefree=(
                    graph_discriminant.sqf_part().total_degree()
                    == graph_discriminant.total_degree()
                ),
                common_discriminant_degree=(
                    vertical_discriminant.gcd(graph_discriminant).total_degree()
                ),
            )
        )
    return tuple(result)


SAMPLE_PARAMETERS: Final = {
    KAPPA: 1,
    ALPHA: 3,
    BETA: Rational(3, 2),
    GAMMA: 3,
    DELTA: 0,
}
SAMPLE_P: Final = expand(FAMILY_P.subs(SAMPLE_PARAMETERS))
SAMPLE_Q: Final = expand(FAMILY_Q.subs(SAMPLE_PARAMETERS))
SAMPLE_CONTACT_ROOTS: Final = (-1, 1)
SAMPLE_CONTACT_SOURCE_QUADRATICS: Final = (
    T**2 + T + 1,
    T**2 - T + 1,
)
SAMPLE_CONTACT_TARGETS: Final = (
    (0, Rational(-1, 2)),
    (-2, Rational(7, 2)),
)
SAMPLE_RESIDUAL_COLLISION: Final = (
    2 * S**6 + 12 * S**5 + 26 * S**4 + 36 * S**3 + 32 * S**2 + 21 * S + 6
) / 2
SAMPLE_COLLISION_POLYNOMIAL: Final = expand((S**2 - 1) ** 2 * SAMPLE_RESIDUAL_COLLISION)
SAMPLE_TANGENCY_COFACTOR: Final = (
    12 * S**9
    + 78 * S**8
    + 199 * S**7
    + 289 * S**6
    + 238 * S**5
    + 91 * S**4
    - 15 * S**3
    - 40 * S**2
    - 32 * S
    - 10
)
SAMPLE_NODE_X_POLYNOMIAL: Final = (
    64 * X**6 + 160 * X**5 - 22352 * X**4 + 14232 * X**3 - 3540 * X**2 + 405 * X - 18
)
SAMPLE_SUPPORT_X_RESULTANT: Final = (
    Rational(729, 64) * X * (X + 2) * SAMPLE_NODE_X_POLYNOMIAL
)
# Scaling the second target coordinate by two preserves the complement and
# gives the primitive integral implicit equation regenerated by Sage 10.8.
SAMPLE_SCALED_Q: Final = expand(2 * SAMPLE_Q)
SAMPLE_IMPLICIT: Final = (
    -16 * X**9
    - 72 * X**8
    - 1176 * X**7
    - 144 * X**6 * Y
    + 387 * X**6
    + 882 * X**5 * Y
    - 36 * X**5
    + 60 * X**4 * Y**2
    - 273 * X**4 * Y
    - 84 * X**3 * Y**2
    + 24 * X**3 * Y
    + 18 * X**2 * Y**3
    + 54 * X**2 * Y**2
    - 12 * X * Y**3
    - 12 * X * Y**2
    + Y**4
    + 2 * Y**3
    + Y**2
)

# Sage 10.8's exact raw affine Zariski--van Kamp presentation.  Generator
# indices 1..4 are the geometric meridians of a generic vertical fiber.
SAMPLE_RELATIONS: Final = (
    (4, 3, -4, -3),
    (2, 1, 2, 1, -2, -1, -2, -1),
    (-2, -1, 3, 4, -3, -4, -3, 1, 2, 1),
    (
        -2,
        -1,
        3,
        4,
        -3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        1,
    ),
    (
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        -1,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        -1,
    ),
    (
        -2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        -4,
        -3,
        1,
    ),
    (
        -2,
        -1,
        3,
        4,
        -3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        1,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        -3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        -1,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
    ),
    (
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        -3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        1,
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        -1,
    ),
    (
        -2,
        -1,
        3,
        -4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        3,
        -4,
        -3,
        1,
        -2,
        -1,
        3,
        4,
        -3,
        1,
        2,
        -1,
        3,
        4,
        -3,
        -4,
        -3,
        1,
    ),
    (-2, -1, 3, 4, -3, 1, 2, -1, 3, -4, -3, 1),
    (-4, -3, 1, 3, 4, -3, -1, 3),
)


def _source_pair_data(pair_sum: int) -> tuple[Expr, Expr, Expr]:
    """Return pair product, source quadratic, and common ``Q`` remainder."""

    pair_product = cancel(
        (S * PAIR_QUADRATIC / PAIR_DENOMINATOR)
        .subs(SAMPLE_PARAMETERS)
        .subs(S, pair_sum)
    )
    source_quadratic = T**2 - pair_sum * T + pair_product
    common_q_remainder = Poly(SAMPLE_Q, T).rem(Poly(source_quadratic, T)).as_expr()
    return pair_product, source_quadratic, common_q_remainder


@dataclass(frozen=True, slots=True)
class DoubleContactSampleCertificate:
    """Exact geometry of the rational ``C2^2 + 6N`` representative."""

    collision_identity: Expr
    repeated_gcd_identity: Expr
    residual_discriminant: Expr
    residual_contact_separation: Expr
    denominator_resultant: Expr
    pair_diagonal_resultant: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    contact_second_derivatives: tuple[Expr, Expr]
    tangency_factor_identity: Expr
    tangency_gcd_identity: Expr
    tangency_gcd_degree: int
    contact_tangency_cofactor_values: tuple[Expr, Expr]
    residual_tangency_resultant: Expr
    contact_pair_products: tuple[Expr, Expr]
    contact_pair_discriminants: tuple[Expr, Expr]
    contact_source_quadratic_identities: tuple[Expr, Expr]
    contact_q_remainders: tuple[Expr, Expr]
    contact_x_values: tuple[Expr, Expr]
    support_x_resultant_identity: Expr
    support_x_discriminant: Expr
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
        """Total projective delta accounting for the sample."""

        return (
            self.cusp_delta + self.contact_delta + self.node_count + self.infinity_delta
        )

    @property
    def verified(self) -> bool:
        """Whether the exact collision geometry and genus budget agree."""

        return bool(
            self.collision_identity == 0
            and self.repeated_gcd_identity == 0
            and self.residual_discriminant == Rational(350724573, 64)
            and self.residual_contact_separation == Rational(-405, 4)
            and self.denominator_resultant == -27
            and self.pair_diagonal_resultant == -85293
            and self.cusp_image_factor == Rational(1, 4)
            and self.extra_critical_factor == 9477
            and self.contact_second_derivatives == (-12, 540)
            and self.tangency_factor_identity == 0
            and self.tangency_gcd_identity == 0
            and self.tangency_gcd_degree == 2
            and self.contact_tangency_cofactor_values == (18, 2430)
            and self.residual_tangency_resultant == Rational(-181729682354700675, 256)
            and self.contact_pair_products == (1, 1)
            and self.contact_pair_discriminants == (-3, -3)
            and self.contact_source_quadratic_identities == (0, 0)
            and self.contact_q_remainders == (0, 0)
            and self.contact_x_values == (0, -2)
            and self.support_x_resultant_identity == 0
            and self.support_x_discriminant
            == Rational(3**99 * 5**24 * 13**2 * 1443311, 2**52)
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.sage_jacobian_components == ((4, 1), (3, 1), (3, 1), (6, 6))
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.contact_delta == 4
            and self.node_count == 6
            and self.infinity_delta == 16
            and self.total_delta == self.arithmetic_genus
            and self.relation_count == 11
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_double_contact_sample_certificate() -> DoubleContactSampleCertificate:
    """Build the exact rational two-contact sample certificate."""

    sample_collision = expand(COLLISION_POLYNOMIAL.subs(SAMPLE_PARAMETERS))
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(SAMPLE_PARAMETERS))
    support_polynomial = expand((S**2 - 1) * SAMPLE_RESIDUAL_COLLISION)
    contact_data = tuple(_source_pair_data(root) for root in SAMPLE_CONTACT_ROOTS)
    collision_x = cancel(
        COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS)
        / COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS)
    )
    support_x_resultant = resultant(
        support_polynomial,
        COLLISION_X_NUMERATOR.subs(SAMPLE_PARAMETERS)
        - X * COLLISION_X_DENOMINATOR.subs(SAMPLE_PARAMETERS),
        S,
    )
    repeated_gcd = (
        gcd(
            Poly(sample_collision, S),
            Poly(diff(sample_collision, S), S),
        )
        .monic()
        .as_expr()
    )
    tangency_gcd = gcd(
        Poly(sample_collision, S),
        Poly(sample_tangency, S),
    ).monic()
    implicit_polynomial = Poly(SAMPLE_IMPLICIT, X, Y)
    return DoubleContactSampleCertificate(
        collision_identity=expand(sample_collision - SAMPLE_COLLISION_POLYNOMIAL),
        repeated_gcd_identity=expand(repeated_gcd - (S**2 - 1)),
        residual_discriminant=discriminant(SAMPLE_RESIDUAL_COLLISION, S),
        residual_contact_separation=resultant(
            SAMPLE_RESIDUAL_COLLISION,
            S**2 - 1,
            S,
        ),
        denominator_resultant=resultant(
            support_polynomial,
            PAIR_DENOMINATOR.subs(SAMPLE_PARAMETERS),
            S,
        ),
        pair_diagonal_resultant=resultant(
            support_polynomial,
            (-S * PAIR_DIAGONAL_FACTOR).subs(SAMPLE_PARAMETERS),
            S,
        ),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(SAMPLE_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(SAMPLE_PARAMETERS),
        contact_second_derivatives=tuple(
            diff(sample_collision, S, 2).subs(S, root) for root in SAMPLE_CONTACT_ROOTS
        ),
        tangency_factor_identity=expand(
            sample_tangency - 3 * (S**2 - 1) * SAMPLE_TANGENCY_COFACTOR
        ),
        tangency_gcd_identity=expand(tangency_gcd.as_expr() - (S**2 - 1)),
        tangency_gcd_degree=tangency_gcd.degree(),
        contact_tangency_cofactor_values=tuple(
            (3 * SAMPLE_TANGENCY_COFACTOR).subs(S, root)
            for root in SAMPLE_CONTACT_ROOTS
        ),
        residual_tangency_resultant=resultant(
            SAMPLE_RESIDUAL_COLLISION,
            sample_tangency,
            S,
        ),
        contact_pair_products=tuple(data[0] for data in contact_data),
        contact_pair_discriminants=tuple(
            factor(root**2 - 4 * data[0])
            for root, data in zip(
                SAMPLE_CONTACT_ROOTS,
                contact_data,
                strict=True,
            )
        ),
        contact_source_quadratic_identities=tuple(
            expand(data[1] - expected)
            for data, expected in zip(
                contact_data,
                SAMPLE_CONTACT_SOURCE_QUADRATICS,
                strict=True,
            )
        ),
        contact_q_remainders=tuple(
            expand(data[2] - target[1])
            for data, target in zip(
                contact_data,
                SAMPLE_CONTACT_TARGETS,
                strict=True,
            )
        ),
        contact_x_values=tuple(
            collision_x.subs(S, root) for root in SAMPLE_CONTACT_ROOTS
        ),
        support_x_resultant_identity=expand(
            support_x_resultant - SAMPLE_SUPPORT_X_RESULTANT
        ),
        support_x_discriminant=discriminant(
            SAMPLE_SUPPORT_X_RESULTANT,
            X,
        ),
        implicit_resultant_identity=expand(
            resultant(SAMPLE_P - X, SAMPLE_SCALED_Q - Y, T) - SAMPLE_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            SAMPLE_IMPLICIT.subs({X: SAMPLE_P, Y: SAMPLE_SCALED_Q})
        ),
        implicit_content=implicit_polynomial.content(),
        # Exact Sage primary-component (length, radical-degree) output:
        # cusp, two contact-two points, and six reduced nodes.
        sage_jacobian_components=((4, 1), (3, 1), (3, 1), (6, 6)),
        # The simplification has one generator and no relators; every raw
        # geometric meridian maps to that generator, with a checked section.
        sage_cyclic_simplification=(1, 0, True),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=(2 - 1) * (5 - 1) // 2,
        contact_delta=2 + 2,
        node_count=6,
        infinity_delta=(5 - 1) * (9 - 1) // 2,
        relation_count=len(SAMPLE_RELATIONS),
        complement_census=n_generator_three_cycle_presentation_census(
            SAMPLE_RELATIONS,
            4,
        ),
    )


HOSTILE_SINGLE_CONTACT_PARAMETERS: Final = {
    KAPPA: 1,
    ALPHA: Rational(-12, 5),
    BETA: Rational(-16, 5),
    GAMMA: 0,
    DELTA: 0,
}
HOSTILE_REPEATED_ROOT: Final = -2


@dataclass(frozen=True, slots=True)
class CoincidentRootHostileCertificate:
    """Fixture showing why ``u-v`` must be saturated."""

    duplicated_incidence_values: tuple[Expr, ...]
    repeated_gcd_identity: Expr
    repeated_gcd_degree: int

    @property
    def verified(self) -> bool:
        """Whether four duplicated equations encode only one contact root."""

        return bool(
            self.duplicated_incidence_values == (0, 0, 0, 0)
            and self.repeated_gcd_identity == 0
            and self.repeated_gcd_degree == 1
        )


@cache
def exact_coincident_root_hostile_certificate() -> CoincidentRootHostileCertificate:
    """Build the one-contact fixture placed twice in the ordered incidence."""

    hostile_collision = expand(
        COLLISION_POLYNOMIAL.subs(HOSTILE_SINGLE_CONTACT_PARAMETERS)
    )
    hostile_gcd = gcd(
        Poly(hostile_collision, S),
        Poly(diff(hostile_collision, S), S),
    ).monic()
    return CoincidentRootHostileCertificate(
        duplicated_incidence_values=(
            hostile_collision.subs(S, HOSTILE_REPEATED_ROOT),
            diff(hostile_collision, S).subs(S, HOSTILE_REPEATED_ROOT),
            hostile_collision.subs(S, HOSTILE_REPEATED_ROOT),
            diff(hostile_collision, S).subs(S, HOSTILE_REPEATED_ROOT),
        ),
        repeated_gcd_identity=expand(
            hostile_gcd.as_expr() - (S - HOSTILE_REPEATED_ROOT)
        ),
        repeated_gcd_degree=hostile_gcd.degree(),
    )


def main() -> int:
    """Print the exact two-contact checkpoint and its remaining gaps."""

    ordered = exact_ordered_incidence_certificate()
    splits = exact_split_boundary_certificates()
    sample = exact_double_contact_sample_certificate()
    hostile = exact_coincident_root_hostile_certificate()
    verified = bool(
        ordered.verified
        and all(certificate.verified for certificate in splits)
        and sample.verified
        and hostile.verified
    )
    print("ordered dominant incidence verified:", ordered.verified)
    print(
        "split discriminant boundaries:",
        tuple(
            (
                certificate.kappa,
                certificate.vertical_discriminant_degree,
                certificate.graph_discriminant_degree,
                certificate.common_discriminant_degree,
            )
            for certificate in splits
        ),
    )
    print(
        "rational C2^2 + 6N sample:",
        {
            "parameters": SAMPLE_PARAMETERS,
            "contacts": SAMPLE_CONTACT_TARGETS,
            "nodes": sample.node_count,
            "delta": sample.total_delta,
        },
    )
    print("coincident-root hostile fixture verified:", hostile.verified)
    print("sample raw-presentation census:", sample.complement_census)
    print("delta-ten double-contact algebra verified:", verified)
    print(
        "remaining: split codimension-two intersections, residual-rank "
        "closures, and the other delta-ten partition strata"
    )
    print(
        "claim boundary: the component dominating the full ordered Cramer "
        "base is conditionally excluded on its generic clean open; "
        "lower-dimensional closures remain open"
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
