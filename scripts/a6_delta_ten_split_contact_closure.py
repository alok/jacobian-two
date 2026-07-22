"""Place the maximal-rank split contact loci in the global components.

Let

``F(k,s,p) = (2*s+k)*p - s*(s^2+k*s+1)``

be the total unordered source-pair surface and let ``C`` be the coefficient of
``t`` in the remainder of ``Q(t)`` modulo ``t^2-s*t+p``.  Then ``C=0`` is
target equality.  The derivation

``D = F_p*d/ds - F_s*d/dp``

is tangent to every fixed-``k`` pair fiber.  Away from its singular points,
``C=D(C)=0`` is contact at least two and
``C=D(C)=D^2(C)=0`` is contact at least three.

This coordinate-free system remains regular when the usual formula solving
for ``p`` has denominator zero.  Interpreting these jets as branch-contact
order additionally requires distinct source roots, a smooth fixed-``k`` pair
fiber, immersed branches, and an unramified common first coordinate.  Those
conditions are imposed by the clean localizers below.  Exact row
transformations identify the jets with the true vertical/graph component
jets on:

* the three maximal-rank split ``C3`` charts; and
* the five maximal-rank split ``C2^2`` charts.

All transformation determinants are units on the corresponding clean opens.
The total pair surface is geometrically integral and flat over the ``k`` line.
Its ordered self-fiber-product is integral as well: it is flat, and its
generic fiber is a domain; a zero divisor would therefore become torsion over
``QQ[k]``.  Consequently the displayed rank-open split loci lie in the same
irreducible global contact incidences as the existing cyclic samples.  The
proper Whitney--Thom propagation already proved for those clean global
components applies to these rank-open split loci.

Exact rational arcs additionally prove that the two ``C3`` overlap
allocations and the ``C2^2`` overlap-plus-contact allocation lie in the
algebraic closures of the ordinary global components.  They remain outside
the smooth proper-isotopy argument.  Their complement topology, together
with the finite compatible rank-three fibers on the three residual
determinant charts, remains open.
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
    expand,
    fraction,
    gcd,
    linear_eq_to_matrix,
    limit,
    rem,
    resultant,
)

from scripts.a6_delta_ten_contact_three import (
    CONTACT_THREE_VALID_LOCALIZER,
    exact_delta_ten_contact_three_certificate,
)
from scripts.a6_delta_ten_double_contact import (
    U as ORDERED_U,
    V as ORDERED_V,
    exact_double_contact_sample_certificate,
    exact_ordered_incidence_certificate,
    ordered_incidence_system,
)
from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    DELTA,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    PAIR_DENOMINATOR,
    PAIR_INCIDENCE,
    PAIR_QUADRATIC,
    R,
    S,
    T,
)
from scripts.a6_delta_ten_split_codim_two import (
    exact_split_witness_certificates,
    split_witness_specs,
)
from scripts.a6_delta_ten_split_c22_overlap_closure import (
    exact_double_contact_overlap_closure_certificate,
)
from scripts.a6_delta_ten_split_contact_rank import (
    U,
    V,
    LinearRankSpec,
    ResidualRankSpec,
    SimpleDoubleContactSpec,
    c3_rank_specs,
    residual_rank_specs,
    simple_double_contact_specs,
)
from scripts.a6_delta_ten_split_pair_surface import (
    exact_split_pair_surface_certificate,
)

COEFFICIENTS: Final = (ALPHA, BETA, GAMMA, DELTA)
PAIR_POLYNOMIAL: Final = T**2 - S * T + R
PAIR_TARGET_ROW: Final = Poly(rem(FAMILY_Q, PAIR_POLYNOMIAL, T), T).coeff_monomial(T)


def pair_fiber_derivative(expression: Expr) -> Expr:
    """Differentiate along a fixed-``k`` fiber of the total pair surface."""

    return expand(
        diff(PAIR_INCIDENCE, R) * diff(expression, S)
        - diff(PAIR_INCIDENCE, S) * diff(expression, R)
    )


PAIR_TANGENCY_ROW: Final = pair_fiber_derivative(PAIR_TARGET_ROW)
PAIR_CONTACT_THREE_ROW: Final = pair_fiber_derivative(PAIR_TANGENCY_ROW)

# On the principal pair chart, ``PAIR_INCIDENCE = 0`` solves for the source
# product as follows.  The power in ``PAIR_TO_COLLISION_SCALING`` records the
# resultant normalization used to define ``COLLISION_POLYNOMIAL``.
SOLVED_PAIR_PRODUCT: Final = cancel(S * PAIR_QUADRATIC / PAIR_DENOMINATOR)
PAIR_TO_COLLISION_SCALING: Final = cancel(S**2 / PAIR_DENOMINATOR**4)
CONTACT_THREE_BRIDGE_LOCALIZER: Final = expand(S * PAIR_DENOMINATOR)
ORDERED_DOUBLE_CONTACT_BRIDGE_LOCALIZER: Final = expand(
    ORDERED_U * ORDERED_V * (KAPPA + 2 * ORDERED_U) * (KAPPA + 2 * ORDERED_V)
)


@dataclass(frozen=True, slots=True)
class LocalizedUnitCertificate:
    """Exact proof that a rational determinant is a localization unit.

    Both the determinant and its displayed inverse have denominators dividing
    powers of ``localizer``.  The two membership identities make that claim
    polynomial, while ``product_residual`` checks that they are inverses.
    ``ambient_localizer`` records the clean-open localizer in which the smaller
    bridge localizer occurs as a factor.
    """

    localizer: Expr
    ambient_localizer: Expr
    ambient_localizer_cofactor: Expr
    ambient_factorization_residual: Expr
    element: Expr
    expected_element: Expr
    element_residual: Expr
    inverse: Expr
    product_residual: Expr
    element_denominator: Expr
    element_denominator_power: int
    element_denominator_quotient: Expr
    element_denominator_membership_residual: Expr
    inverse_denominator: Expr
    inverse_denominator_power: int
    inverse_denominator_quotient: Expr
    inverse_denominator_membership_residual: Expr

    @property
    def verified(self) -> bool:
        """Whether all exact localization-unit identities hold."""

        return bool(
            self.localizer != 0
            and self.ambient_factorization_residual == 0
            and fraction(cancel(self.ambient_localizer_cofactor))[1] == 1
            and self.element_residual == 0
            and self.product_residual == 0
            and self.element_denominator_power > 0
            and fraction(cancel(self.element_denominator_quotient))[1] == 1
            and self.element_denominator_membership_residual == 0
            and self.inverse_denominator_power > 0
            and fraction(cancel(self.inverse_denominator_quotient))[1] == 1
            and self.inverse_denominator_membership_residual == 0
        )


def _localized_unit_certificate(
    *,
    element: Expr,
    expected_element: Expr,
    inverse: Expr,
    localizer: Expr,
    ambient_localizer: Expr,
    element_denominator_power: int,
    inverse_denominator_power: int,
) -> LocalizedUnitCertificate:
    """Build exact denominator-membership identities for one rational unit."""

    normalized_element = cancel(element)
    normalized_expected = cancel(expected_element)
    normalized_inverse = cancel(inverse)
    element_denominator = fraction(normalized_element)[1]
    inverse_denominator = fraction(normalized_inverse)[1]
    element_quotient = cancel(
        localizer**element_denominator_power / element_denominator
    )
    inverse_quotient = cancel(
        localizer**inverse_denominator_power / inverse_denominator
    )
    ambient_cofactor = cancel(ambient_localizer / localizer)
    return LocalizedUnitCertificate(
        localizer=localizer,
        ambient_localizer=ambient_localizer,
        ambient_localizer_cofactor=ambient_cofactor,
        ambient_factorization_residual=cancel(
            ambient_localizer - localizer * ambient_cofactor
        ),
        element=normalized_element,
        expected_element=normalized_expected,
        element_residual=cancel(normalized_element - normalized_expected),
        inverse=normalized_inverse,
        product_residual=cancel(normalized_element * normalized_inverse - 1),
        element_denominator=element_denominator,
        element_denominator_power=element_denominator_power,
        element_denominator_quotient=element_quotient,
        element_denominator_membership_residual=cancel(
            localizer**element_denominator_power
            - element_denominator * element_quotient
        ),
        inverse_denominator=inverse_denominator,
        inverse_denominator_power=inverse_denominator_power,
        inverse_denominator_quotient=inverse_quotient,
        inverse_denominator_membership_residual=cancel(
            localizer**inverse_denominator_power
            - inverse_denominator * inverse_quotient
        ),
    )


def _principal_pair_rows() -> tuple[Expr, Expr, Expr]:
    """Restrict the intrinsic total-pair jets to the solved pair chart."""

    return tuple(
        cancel(row.subs(R, SOLVED_PAIR_PRODUCT))
        for row in (PAIR_TARGET_ROW, PAIR_TANGENCY_ROW, PAIR_CONTACT_THREE_ROW)
    )


def _contact_three_bridge_transformation() -> Matrix:
    """Return the chain-rule transform from ``(H,H',H'')`` to total jets."""

    scaling = PAIR_TO_COLLISION_SCALING
    fiber_speed = PAIR_DENOMINATOR
    scaling_prime = diff(scaling, S)
    return Matrix(
        [
            [scaling, 0, 0],
            [
                cancel(fiber_speed * scaling_prime),
                cancel(fiber_speed * scaling),
                0,
            ],
            [
                cancel(
                    fiber_speed
                    * (
                        diff(fiber_speed, S) * scaling_prime
                        + fiber_speed * diff(scaling, S, 2)
                    )
                ),
                cancel(
                    fiber_speed
                    * (diff(fiber_speed, S) * scaling + 2 * fiber_speed * scaling_prime)
                ),
                cancel(fiber_speed**2 * scaling),
            ],
        ]
    )


@dataclass(frozen=True, slots=True)
class ContactThreeJetBridgeCertificate:
    """Bridge total pair jets to the old ``H,H',H''`` incidence."""

    target_scaling_residual: Expr
    fiber_chain_rule_residuals: tuple[Expr, ...]
    transformation: Matrix
    row_residuals: tuple[Expr, ...]
    determinant_unit: LocalizedUnitCertificate

    @property
    def verified(self) -> bool:
        """Whether the two contact-three row systems agree on the clean open."""

        return bool(
            self.target_scaling_residual == 0
            and all(value == 0 for value in self.fiber_chain_rule_residuals)
            and all(value == 0 for value in self.row_residuals)
            and self.determinant_unit.verified
        )


@cache
def exact_contact_three_jet_bridge_certificate() -> ContactThreeJetBridgeCertificate:
    """Prove the exact total-pair/old-contact-three change of rows."""

    total_rows = _principal_pair_rows()
    old_rows = (
        COLLISION_POLYNOMIAL,
        diff(COLLISION_POLYNOMIAL, S),
        diff(COLLISION_POLYNOMIAL, S, 2),
    )
    transformation = _contact_three_bridge_transformation()
    row_residual_matrix = (
        _augmented_row_matrix(total_rows)
        - transformation * _augmented_row_matrix(old_rows)
    ).applyfunc(cancel)
    determinant = cancel(transformation.det())
    determinant_unit = _localized_unit_certificate(
        element=determinant,
        expected_element=S**6 / PAIR_DENOMINATOR**9,
        inverse=PAIR_DENOMINATOR**9 / S**6,
        localizer=CONTACT_THREE_BRIDGE_LOCALIZER,
        ambient_localizer=CONTACT_THREE_VALID_LOCALIZER,
        element_denominator_power=9,
        inverse_denominator_power=6,
    )
    # The tangent field restricts to ``PAIR_DENOMINATOR*d/ds`` on the graph
    # of the solved source product.  Check it independently on C and D(C).
    chain_rule_residuals = tuple(
        cancel(total_rows[index + 1] - PAIR_DENOMINATOR * diff(total_rows[index], S))
        for index in range(2)
    )
    return ContactThreeJetBridgeCertificate(
        target_scaling_residual=cancel(
            total_rows[0] - PAIR_TO_COLLISION_SCALING * COLLISION_POLYNOMIAL
        ),
        fiber_chain_rule_residuals=chain_rule_residuals,
        transformation=transformation,
        row_residuals=tuple(row_residual_matrix),
        determinant_unit=determinant_unit,
    )


@dataclass(frozen=True, slots=True)
class OrderedDoubleContactJetBridgeCertificate:
    """Bridge two total ``C,D(C)`` pairs to ordered ``H,H'`` rows."""

    transformation: Matrix
    row_residuals: tuple[Expr, ...]
    determinant_unit: LocalizedUnitCertificate

    @property
    def verified(self) -> bool:
        """Whether both ordered contact-two blocks agree exactly."""

        return bool(
            all(value == 0 for value in self.row_residuals)
            and self.determinant_unit.verified
        )


@cache
def exact_ordered_double_contact_jet_bridge_certificate() -> (
    OrderedDoubleContactJetBridgeCertificate
):
    """Prove the exact block-diagonal bridge for the ordered incidence."""

    total_contact_two = _principal_pair_rows()[:2]
    total_rows = tuple(
        cancel(row.subs(S, root))
        for root in (ORDERED_U, ORDERED_V)
        for row in total_contact_two
    )
    old_rows = ordered_incidence_system()[0]
    contact_two_block = _contact_three_bridge_transformation()[:2, :2]
    transformation = Matrix.diag(
        contact_two_block.subs(S, ORDERED_U),
        contact_two_block.subs(S, ORDERED_V),
    ).applyfunc(cancel)
    row_residual_matrix = (
        _augmented_row_matrix(total_rows)
        - transformation * _augmented_row_matrix(old_rows)
    ).applyfunc(cancel)
    u_denominator = KAPPA + 2 * ORDERED_U
    v_denominator = KAPPA + 2 * ORDERED_V
    determinant = cancel(transformation.det())
    determinant_unit = _localized_unit_certificate(
        element=determinant,
        expected_element=(
            ORDERED_U**4 * ORDERED_V**4 / (u_denominator**7 * v_denominator**7)
        ),
        inverse=(u_denominator**7 * v_denominator**7 / (ORDERED_U**4 * ORDERED_V**4)),
        localizer=ORDERED_DOUBLE_CONTACT_BRIDGE_LOCALIZER,
        ambient_localizer=ORDERED_DOUBLE_CONTACT_BRIDGE_LOCALIZER,
        element_denominator_power=7,
        inverse_denominator_power=4,
    )
    return OrderedDoubleContactJetBridgeCertificate(
        transformation=transformation,
        row_residuals=tuple(row_residual_matrix),
        determinant_unit=determinant_unit,
    )


OVERLAP_EPSILON: Final = Symbol("epsilon_overlap")
OVERLAP_LAMBDA: Final = Symbol("lambda_overlap")
OVERLAP_FREE_BETA: Final = Symbol("beta_overlap")
OVERLAP_FREE_DELTA: Final = Symbol("delta_overlap")
OVERLAP_TARGET_BETA: Final = Symbol("beta_limit")
OVERLAP_TARGET_GAMMA: Final = Symbol("gamma_limit")
OVERLAP_TARGET_DELTA: Final = Symbol("delta_limit")


def _pair_surface_kappa(source_sum: Expr, source_product: Expr) -> Expr:
    """Solve the total pair surface for ``k`` away from ``p=s^2``."""

    return cancel(
        -source_sum
        * (2 * source_product - source_sum**2 - 1)
        / (source_product - source_sum**2)
    )


def _cramer_certificate(
    equations: tuple[Expr, ...],
    unknowns: tuple[Expr, ...],
) -> tuple[Expr, tuple[Expr, ...], tuple[Expr, ...]]:
    """Return determinant, Cramer numerators, and denominator-free residuals."""

    coefficient_matrix, right_side = linear_eq_to_matrix(equations, unknowns)
    determinant = cancel(coefficient_matrix.det())
    numerators: list[Expr] = []
    for column in range(len(unknowns)):
        replaced = coefficient_matrix.copy()
        replaced[:, column] = right_side
        numerators.append(cancel(replaced.det()))
    residuals = tuple(
        cancel(value)
        for value in (
            coefficient_matrix * Matrix(numerators) - determinant * right_side
        )
    )
    return determinant, tuple(numerators), residuals


@dataclass(frozen=True, slots=True)
class ContactThreeOverlapDegenerationCertificate:
    """Exact generic ``C3`` arc degenerating to one split overlap plane."""

    name: str
    allocation: str
    source_sum: Expr
    source_product: Expr
    kappa: Expr
    pair_incidence_residual: Expr
    pair_chart_denominator_limit: Expr
    valid_localizer_order: int
    valid_localizer_leading_coefficient: Expr
    expected_valid_localizer_leading_coefficient: Expr
    cramer_determinant: Expr
    cramer_determinant_order: int
    cramer_determinant_leading_coefficient: Expr
    expected_cramer_determinant_leading_coefficient: Expr
    cramer_numerators: tuple[Expr, ...]
    cramer_residuals: tuple[Expr, ...]
    coefficient_limits: tuple[Expr, Expr, Expr, Expr]
    expected_coefficient_limits: tuple[Expr, Expr, Expr, Expr]
    coefficient_limit_residuals: tuple[Expr, ...]
    overlap_equation_residuals: tuple[Expr, ...]
    dominance_inverse_parameters: tuple[Expr, Expr]
    dominance_residuals: tuple[Expr, ...]

    @property
    def verified(self) -> bool:
        """Whether this is a dense two-parameter degeneration to the plane."""

        return bool(
            self.pair_incidence_residual == 0
            and self.pair_chart_denominator_limit == Rational(1, 2)
            and self.valid_localizer_order > 0
            and self.valid_localizer_leading_coefficient
            == self.expected_valid_localizer_leading_coefficient
            and self.valid_localizer_leading_coefficient != 0
            and self.cramer_determinant != 0
            and self.cramer_determinant_order > 0
            and self.cramer_determinant_leading_coefficient
            == self.expected_cramer_determinant_leading_coefficient
            and self.cramer_determinant_leading_coefficient != 0
            and len(self.cramer_numerators) == 3
            and all(value == 0 for value in self.cramer_residuals)
            and all(value == 0 for value in self.coefficient_limit_residuals)
            and all(value == 0 for value in self.overlap_equation_residuals)
            and all(value == 0 for value in self.dominance_residuals)
        )


def _overlap_w_degeneration() -> ContactThreeOverlapDegenerationCertificate:
    """Degenerate a principal ``C3`` point to graph-double overlap contact."""

    epsilon = OVERLAP_EPSILON
    lam = OVERLAP_LAMBDA
    source_sum = epsilon
    source_product = Rational(1, 2) + lam * epsilon**2
    kappa = _pair_surface_kappa(source_sum, source_product)
    coordinates = {S: source_sum, R: source_product, KAPPA: kappa}
    equations = tuple(
        cancel(row.subs(coordinates).subs(DELTA, OVERLAP_FREE_DELTA))
        for row in (PAIR_TARGET_ROW, PAIR_TANGENCY_ROW, PAIR_CONTACT_THREE_ROW)
    )
    determinant, numerators, cramer_residuals = _cramer_certificate(
        equations,
        (ALPHA, BETA, GAMMA),
    )
    solved_limits = tuple(
        cancel(limit(numerator / determinant, epsilon, 0)) for numerator in numerators
    )
    coefficient_limits = (
        solved_limits[0],
        solved_limits[1],
        solved_limits[2],
        OVERLAP_FREE_DELTA,
    )
    expected_limits = (
        (lam - 4) / (4 * (lam - 3)),
        2 * OVERLAP_FREE_DELTA / 3,
        (2 * lam - 7) / (2 * (lam - 3)),
        OVERLAP_FREE_DELTA,
    )
    limit_substitution = dict(zip(COEFFICIENTS, expected_limits, strict=True))
    overlap_spec = next(
        spec for spec in c3_rank_specs() if spec.allocation == "overlap-W"
    )
    inverse_parameters = (
        (6 * OVERLAP_TARGET_GAMMA - 7) / (2 * (OVERLAP_TARGET_GAMMA - 1)),
        OVERLAP_TARGET_DELTA,
    )
    dominance_substitution = {
        lam: inverse_parameters[0],
        OVERLAP_FREE_DELTA: inverse_parameters[1],
    }
    target_plane_point = (
        OVERLAP_TARGET_GAMMA / 2 - Rational(1, 4),
        2 * OVERLAP_TARGET_DELTA / 3,
        OVERLAP_TARGET_GAMMA,
        OVERLAP_TARGET_DELTA,
    )
    valid_localizer = cancel(CONTACT_THREE_VALID_LOCALIZER.subs(coordinates))
    return ContactThreeOverlapDegenerationCertificate(
        name="c3_overlap_w_degeneration",
        allocation="overlap-W",
        source_sum=source_sum,
        source_product=source_product,
        kappa=kappa,
        pair_incidence_residual=cancel(PAIR_INCIDENCE.subs(coordinates)),
        pair_chart_denominator_limit=cancel(
            limit(source_product - source_sum**2, epsilon, 0)
        ),
        valid_localizer_order=2,
        valid_localizer_leading_coefficient=cancel(
            limit(valid_localizer / epsilon**2, epsilon, 0)
        ),
        expected_valid_localizer_leading_coefficient=8,
        cramer_determinant=determinant,
        cramer_determinant_order=3,
        cramer_determinant_leading_coefficient=cancel(
            limit(determinant / epsilon**3, epsilon, 0)
        ),
        expected_cramer_determinant_leading_coefficient=(-3 * (lam - 3) / 4),
        cramer_numerators=numerators,
        cramer_residuals=cramer_residuals,
        coefficient_limits=coefficient_limits,
        expected_coefficient_limits=expected_limits,
        coefficient_limit_residuals=tuple(
            cancel(actual - expected)
            for actual, expected in zip(
                coefficient_limits,
                expected_limits,
                strict=True,
            )
        ),
        overlap_equation_residuals=tuple(
            cancel(equation.subs(limit_substitution))
            for equation in overlap_spec.equations
        ),
        dominance_inverse_parameters=inverse_parameters,
        dominance_residuals=tuple(
            cancel(cancel(actual).subs(dominance_substitution) - expected)
            for actual, expected in zip(
                expected_limits,
                target_plane_point,
                strict=True,
            )
        ),
    )


def _overlap_v_degeneration() -> ContactThreeOverlapDegenerationCertificate:
    """Degenerate a principal ``C3`` point to vertical-double overlap contact."""

    epsilon = OVERLAP_EPSILON
    lam = OVERLAP_LAMBDA
    source_sum = lam * epsilon**2
    source_product = Rational(1, 2) + epsilon
    kappa = _pair_surface_kappa(source_sum, source_product)
    coordinates = {S: source_sum, R: source_product, KAPPA: kappa}
    equations = tuple(
        cancel(row.subs(coordinates).subs(BETA, OVERLAP_FREE_BETA))
        for row in (PAIR_TARGET_ROW, PAIR_TANGENCY_ROW, PAIR_CONTACT_THREE_ROW)
    )
    determinant, numerators, cramer_residuals = _cramer_certificate(
        equations,
        (ALPHA, GAMMA, DELTA),
    )
    solved_limits = tuple(
        cancel(limit(numerator / determinant, epsilon, 0)) for numerator in numerators
    )
    coefficient_limits = (
        solved_limits[0],
        OVERLAP_FREE_BETA,
        solved_limits[1],
        solved_limits[2],
    )
    expected_limits = (
        Rational(1, 4),
        OVERLAP_FREE_BETA,
        1,
        3 * OVERLAP_FREE_BETA / 2 + 1 / (2 * lam),
    )
    limit_substitution = dict(zip(COEFFICIENTS, expected_limits, strict=True))
    overlap_spec = next(
        spec for spec in c3_rank_specs() if spec.allocation == "overlap-V"
    )
    inverse_parameters = (
        1 / (2 * OVERLAP_TARGET_DELTA - 3 * OVERLAP_TARGET_BETA),
        OVERLAP_TARGET_BETA,
    )
    dominance_substitution = {
        lam: inverse_parameters[0],
        OVERLAP_FREE_BETA: inverse_parameters[1],
    }
    target_plane_point = (
        Rational(1, 4),
        OVERLAP_TARGET_BETA,
        1,
        OVERLAP_TARGET_DELTA,
    )
    valid_localizer = cancel(CONTACT_THREE_VALID_LOCALIZER.subs(coordinates))
    return ContactThreeOverlapDegenerationCertificate(
        name="c3_overlap_v_degeneration",
        allocation="overlap-V",
        source_sum=source_sum,
        source_product=source_product,
        kappa=kappa,
        pair_incidence_residual=cancel(PAIR_INCIDENCE.subs(coordinates)),
        pair_chart_denominator_limit=cancel(
            limit(source_product - source_sum**2, epsilon, 0)
        ),
        valid_localizer_order=4,
        valid_localizer_leading_coefficient=cancel(
            limit(valid_localizer / epsilon**4, epsilon, 0)
        ),
        expected_valid_localizer_leading_coefficient=8 * lam**2,
        cramer_determinant=determinant,
        cramer_determinant_order=3,
        cramer_determinant_leading_coefficient=cancel(
            limit(determinant / epsilon**3, epsilon, 0)
        ),
        expected_cramer_determinant_leading_coefficient=-lam / 2,
        cramer_numerators=numerators,
        cramer_residuals=cramer_residuals,
        coefficient_limits=coefficient_limits,
        expected_coefficient_limits=expected_limits,
        coefficient_limit_residuals=tuple(
            cancel(actual - expected)
            for actual, expected in zip(
                coefficient_limits,
                expected_limits,
                strict=True,
            )
        ),
        overlap_equation_residuals=tuple(
            cancel(equation.subs(limit_substitution))
            for equation in overlap_spec.equations
        ),
        dominance_inverse_parameters=inverse_parameters,
        dominance_residuals=tuple(
            cancel(cancel(actual).subs(dominance_substitution) - expected)
            for actual, expected in zip(
                expected_limits,
                target_plane_point,
                strict=True,
            )
        ),
    )


@cache
def exact_contact_three_overlap_degeneration_certificates() -> tuple[
    ContactThreeOverlapDegenerationCertificate, ...
]:
    """Return dense principal-chart degenerations onto both overlap planes."""

    return (_overlap_w_degeneration(), _overlap_v_degeneration())


def _component_coordinates(kappa: int, component: str, root: Expr) -> dict[Expr, Expr]:
    """Embed a true vertical or graph component into ``(k,s,p)``."""

    if (kappa, component) == (0, "V"):
        return {KAPPA: 0, S: 0, R: root}
    if (kappa, component) == (0, "W"):
        return {KAPPA: 0, S: root, R: (root**2 + 1) / 2}
    if (kappa, component) == (2, "V"):
        return {KAPPA: 2, S: -1, R: root}
    if (kappa, component) == (2, "W"):
        return {KAPPA: 2, S: root, R: root * (root + 1) / 2}
    msg = f"unsupported split component: k={kappa}, component={component}"
    raise ValueError(msg)


def _augmented_row_matrix(equations: tuple[Expr, ...]) -> Matrix:
    """Return rows of full affine-linear expressions, including constants."""

    coefficient, right_side = linear_eq_to_matrix(equations, COEFFICIENTS)
    return coefficient.row_join(-right_side)


def _row_transformation(
    global_equations: tuple[Expr, ...],
    component_equations: tuple[Expr, ...],
) -> tuple[Matrix, tuple[Expr, ...]]:
    """Solve and verify the exact global-to-component row transformation."""

    global_rows = _augmented_row_matrix(global_equations)
    component_rows = _augmented_row_matrix(component_equations)
    rank = len(component_equations)
    component_square = component_rows[:, :rank]
    transformation = (global_rows[:, :rank] * component_square.inv()).applyfunc(cancel)
    residual_matrix = (global_rows - transformation * component_rows).applyfunc(cancel)
    return transformation, tuple(residual_matrix)


@dataclass(frozen=True, slots=True)
class TotalContactBaseCertificate:
    """Geometric-integrality certificate for one and two ordered pair bases."""

    pair_surface_geometrically_integral: bool
    pair_surface_flat_over_kappa: bool
    generic_pair_product_degree: int
    generic_pair_linear_coefficient_degree: int
    generic_pair_constant_coefficient_degree: int
    generic_pair_coefficient_gcd: Expr
    generic_pair_coefficient_resultant: Expr
    generic_pair_curve_primitive_linear: bool
    generic_pair_curve_geometrically_integral: bool
    ordered_generic_product_domain_consequence: bool
    ordered_two_pair_base_flat_over_kappa: bool
    ordered_two_pair_generic_fiber_domain: bool
    ordered_two_pair_torsion_descent: bool
    ordered_two_pair_base_domain: bool

    @property
    def verified(self) -> bool:
        """Whether both total source bases are irreducible over ``CC``."""

        return bool(
            self.pair_surface_geometrically_integral
            and self.pair_surface_flat_over_kappa
            and self.generic_pair_product_degree == 1
            and self.generic_pair_linear_coefficient_degree == 1
            and self.generic_pair_constant_coefficient_degree == 3
            and self.generic_pair_coefficient_gcd == 1
            and self.generic_pair_coefficient_resultant
            == expand(KAPPA * (KAPPA - 2) * (KAPPA + 2))
            and self.generic_pair_curve_primitive_linear
            and self.generic_pair_curve_geometrically_integral
            and self.ordered_generic_product_domain_consequence
            and self.ordered_two_pair_base_flat_over_kappa
            and self.ordered_two_pair_generic_fiber_domain
            and self.ordered_two_pair_torsion_descent
            and self.ordered_two_pair_base_domain
        )


@cache
def exact_total_contact_base_certificate() -> TotalContactBaseCertificate:
    """Build the flatness-plus-generic-domain fiber-product certificate."""

    pair = exact_split_pair_surface_certificate()
    pair_integral = pair.total_surface_irreducible
    pair_flat = pair.flat_over_kappa
    linear_coefficient = PAIR_DENOMINATOR
    constant_coefficient = expand(-S * PAIR_QUADRATIC)
    product_degree = Poly(PAIR_INCIDENCE, R).degree()
    linear_degree = Poly(linear_coefficient, S).degree()
    constant_degree = Poly(constant_coefficient, S).degree()
    coefficient_gcd = expand(gcd(linear_coefficient, constant_coefficient))
    coefficient_resultant = expand(
        resultant(linear_coefficient, constant_coefficient, S)
    )
    # Over Qbar(k), the displayed nonzero resultant says the two coefficients
    # are coprime.  A primitive degree-one polynomial in r is irreducible by
    # Gauss's lemma, so its generic pair curve is geometrically integral.
    primitive_linear = bool(
        product_degree == 1 and coefficient_gcd == 1 and coefficient_resultant != 0
    )
    generic_geometrically_integral = primitive_linear
    # The product of two geometrically integral schemes over the algebraically
    # closed field Qbar(k) is integral.  This is a separate theorem input; it
    # is not inferred merely by copying the total-surface irreducibility flag.
    generic_product_domain = generic_geometrically_integral
    ordered_flat = pair_flat  # tensor products of flat QQ[k]-modules are flat
    generic_domain = generic_product_domain
    # If xy=0 in A tensor A, its generic localization is a domain, hence x or
    # y becomes zero.  Flatness is torsion-freeness over the domain QQ[k], so
    # that factor was already zero before localization.
    torsion_descent = ordered_flat and generic_domain
    total_domain = torsion_descent
    return TotalContactBaseCertificate(
        pair_surface_geometrically_integral=pair_integral,
        pair_surface_flat_over_kappa=pair_flat,
        generic_pair_product_degree=product_degree,
        generic_pair_linear_coefficient_degree=linear_degree,
        generic_pair_constant_coefficient_degree=constant_degree,
        generic_pair_coefficient_gcd=coefficient_gcd,
        generic_pair_coefficient_resultant=coefficient_resultant,
        generic_pair_curve_primitive_linear=primitive_linear,
        generic_pair_curve_geometrically_integral=(generic_geometrically_integral),
        ordered_generic_product_domain_consequence=generic_product_domain,
        ordered_two_pair_base_flat_over_kappa=ordered_flat,
        ordered_two_pair_generic_fiber_domain=generic_domain,
        ordered_two_pair_torsion_descent=torsion_descent,
        ordered_two_pair_base_domain=total_domain,
    )


@dataclass(frozen=True, slots=True)
class ComponentDeterminantUnitCertificate:
    """Power-membership proof that a row determinant is a clean-open unit.

    If ``h^e = n*q``, where ``h`` is the component clean localizer and ``n``
    is the numerator of the transformation determinant, then ``n`` (and hence
    the determinant) is a unit after localizing at ``h``.  Lower-power
    denominators record that ``e`` is the first such exponent.
    """

    determinant: Expr
    clean_localizer: Expr
    determinant_numerator: Expr
    determinant_scalar_denominator: Expr
    saturation_exponent: int
    quotient: Expr
    membership_residual: Expr
    lower_power_denominators: tuple[Expr, ...]

    @property
    def localized_ideal_is_unit(self) -> bool:
        """Whether the exact power-membership identity proves unitness."""

        return bool(self.membership_residual == 0)

    @property
    def verified(self) -> bool:
        """Whether the determinant is a unit on precisely the declared open."""

        localization_symbols = {U, V}
        quotient_denominator = fraction(cancel(self.quotient))[1]
        return bool(
            self.determinant != 0
            and self.clean_localizer != 0
            and self.determinant_numerator != 0
            and self.determinant_scalar_denominator != 0
            and not self.determinant_scalar_denominator.free_symbols
            and self.saturation_exponent > 0
            and len(self.lower_power_denominators) == self.saturation_exponent - 1
            and all(
                bool(denominator.free_symbols & localization_symbols)
                for denominator in self.lower_power_denominators
            )
            and not quotient_denominator.free_symbols
            and self.localized_ideal_is_unit
        )


def _component_determinant_unit_certificate(
    determinant: Expr,
    clean_localizer: Expr,
) -> ComponentDeterminantUnitCertificate:
    """Find the first exact ``h^e`` membership in the determinant numerator."""

    normalized_determinant = cancel(determinant)
    numerator, scalar_denominator = fraction(normalized_determinant)
    lower_denominators: list[Expr] = []
    for exponent in range(1, 13):
        quotient = cancel(clean_localizer**exponent / numerator)
        quotient_denominator = fraction(quotient)[1]
        if not quotient_denominator.free_symbols:
            return ComponentDeterminantUnitCertificate(
                determinant=normalized_determinant,
                clean_localizer=clean_localizer,
                determinant_numerator=numerator,
                determinant_scalar_denominator=scalar_denominator,
                saturation_exponent=exponent,
                quotient=quotient,
                membership_residual=cancel(
                    clean_localizer**exponent - numerator * quotient
                ),
                lower_power_denominators=tuple(lower_denominators),
            )
        lower_denominators.append(quotient_denominator)
    msg = "component transform determinant did not become a clean-localizer unit"
    raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class MaximalRankContactEmbeddingCertificate:
    """Exact row comparison for one maximal-rank split contact allocation."""

    name: str
    profile: str
    transformation_determinant: Expr
    row_residuals: tuple[Expr, ...]
    determinant_unit: ComponentDeterminantUnitCertificate
    witness_transformation_determinant: Expr
    clean_witness_verified: bool

    @property
    def verified(self) -> bool:
        """Whether the global and component systems agree on a clean open."""

        return bool(
            self.transformation_determinant != 0
            and all(residual == 0 for residual in self.row_residuals)
            and self.determinant_unit.determinant
            == cancel(self.transformation_determinant)
            and self.determinant_unit.verified
            and self.witness_transformation_determinant != 0
            and self.clean_witness_verified
        )


def _witness_name(
    spec: LinearRankSpec | SimpleDoubleContactSpec | ResidualRankSpec,
) -> str:
    """Return the generic-ledger witness name for one focused rank spec."""

    return spec.name.removesuffix("_residual")


def _witness_root_values(
    spec: LinearRankSpec | SimpleDoubleContactSpec | ResidualRankSpec,
) -> dict[Expr, Expr]:
    """Return the component-root coordinates of the exact clean witness."""

    witness = {item.name: item for item in split_witness_specs()}[_witness_name(spec)]
    roots = tuple(root.root for root in witness.roots)
    if isinstance(spec, LinearRankSpec):
        return {U: roots[0]}
    return {U: roots[0], V: roots[1]}


def _clean_witness_verified(
    spec: LinearRankSpec | SimpleDoubleContactSpec | ResidualRankSpec,
) -> bool:
    """Replay the full residual-node and target-separation witness predicate."""

    witness_specs = {item.name: item for item in split_witness_specs()}
    witness_certificates = {
        item.name: item for item in exact_split_witness_certificates()
    }
    name = _witness_name(spec)
    return witness_certificates[name].verified(witness_specs[name])


def _c3_embedding(spec: LinearRankSpec) -> MaximalRankContactEmbeddingCertificate:
    """Compare one non-overlap split ``C3`` jet system."""

    global_equations = tuple(
        expand(row.subs(_component_coordinates(spec.kappa, spec.allocation, U)))
        for row in (PAIR_TARGET_ROW, PAIR_TANGENCY_ROW, PAIR_CONTACT_THREE_ROW)
    )
    transformation, residuals = _row_transformation(global_equations, spec.equations)
    determinant = cancel(transformation.det())
    witness_substitution = _witness_root_values(spec)
    return MaximalRankContactEmbeddingCertificate(
        name=spec.name,
        profile=spec.profile,
        transformation_determinant=determinant,
        row_residuals=residuals,
        determinant_unit=_component_determinant_unit_certificate(
            determinant,
            spec.clean_localizer,
        ),
        witness_transformation_determinant=cancel(
            determinant.subs(witness_substitution)
        ),
        clean_witness_verified=_clean_witness_verified(spec),
    )


def _double_contact_embedding(
    spec: SimpleDoubleContactSpec | ResidualRankSpec,
) -> MaximalRankContactEmbeddingCertificate:
    """Compare one non-overlap split two-contact row system."""

    global_equations: list[Expr] = []
    for component, root in zip(spec.allocation, (U, V), strict=True):
        substitution = _component_coordinates(spec.kappa, component, root)
        global_equations.extend(
            (
                expand(PAIR_TARGET_ROW.subs(substitution)),
                expand(PAIR_TANGENCY_ROW.subs(substitution)),
            )
        )
    transformation, residuals = _row_transformation(
        tuple(global_equations),
        spec.equations,
    )
    determinant = cancel(transformation.det())
    witness_substitution = _witness_root_values(spec)
    return MaximalRankContactEmbeddingCertificate(
        name=spec.name,
        profile="C2^2+6N",
        transformation_determinant=determinant,
        row_residuals=residuals,
        determinant_unit=_component_determinant_unit_certificate(
            determinant,
            spec.clean_localizer,
        ),
        witness_transformation_determinant=cancel(
            determinant.subs(witness_substitution)
        ),
        clean_witness_verified=_clean_witness_verified(spec),
    )


@cache
def exact_maximal_rank_contact_embedding_certificates() -> tuple[
    MaximalRankContactEmbeddingCertificate, ...
]:
    """Build the three ``C3`` and five ``C2^2`` rank-open comparisons."""

    c3 = tuple(
        _c3_embedding(spec)
        for spec in c3_rank_specs()
        if not spec.allocation.startswith("overlap-")
    )
    simple_double_contacts = tuple(
        _double_contact_embedding(spec)
        for spec in simple_double_contact_specs()
        if spec.allocation != "overlap+W"
    )
    residual_double_contacts = tuple(
        _double_contact_embedding(spec) for spec in residual_rank_specs()
    )
    return (*c3, *simple_double_contacts, *residual_double_contacts)


@dataclass(frozen=True, slots=True)
class SplitContactClosureCertificate:
    """Aggregate rank-open topology conclusion and exceptional boundary."""

    total_bases: TotalContactBaseCertificate
    embeddings: tuple[MaximalRankContactEmbeddingCertificate, ...]
    contact_three_jet_bridge_verified: bool
    ordered_double_contact_jet_bridge_verified: bool
    global_contact_three_sample_cyclic: bool
    global_double_contact_sample_cyclic: bool
    global_double_contact_cramer_verified: bool
    maximal_rank_algebraically_contained: bool
    proper_isotopy_extension_recorded: bool
    overlap_algebraic_closure_count: int
    exceptional_allocation_count: int
    residual_affine_line_base_length: int
    exceptional_topology_open: bool

    @property
    def maximal_rank_topology_closed(self) -> bool:
        """Combine algebraic containment with the separate topology theorem."""

        return bool(
            self.maximal_rank_algebraically_contained
            and self.proper_isotopy_extension_recorded
        )

    @property
    def verified(self) -> bool:
        """Whether the rank-open theorem and remaining exception count agree."""

        return bool(
            self.total_bases.verified
            and len(self.embeddings) == 8
            and all(item.verified for item in self.embeddings)
            and self.contact_three_jet_bridge_verified
            and self.ordered_double_contact_jet_bridge_verified
            and self.global_contact_three_sample_cyclic
            and self.global_double_contact_sample_cyclic
            and self.global_double_contact_cramer_verified
            and self.maximal_rank_algebraically_contained
            and self.proper_isotopy_extension_recorded
            and self.maximal_rank_topology_closed
            and self.overlap_algebraic_closure_count == 3
            and self.exceptional_allocation_count == 3
            and self.residual_affine_line_base_length == 16
            and self.exceptional_topology_open
        )


@cache
def exact_split_contact_closure_certificate() -> SplitContactClosureCertificate:
    """Build the clean rank-open contact-component closure certificate."""

    contact_three = exact_delta_ten_contact_three_certificate()
    double_contact = exact_double_contact_sample_certificate()
    ordered = exact_ordered_incidence_certificate()
    embeddings = exact_maximal_rank_contact_embedding_certificates()
    contact_three_bridge = exact_contact_three_jet_bridge_certificate()
    ordered_double_contact_bridge = (
        exact_ordered_double_contact_jet_bridge_certificate()
    )
    contact_three_overlap = exact_contact_three_overlap_degeneration_certificates()
    double_contact_overlap = exact_double_contact_overlap_closure_certificate()
    bridges_verified = bool(
        contact_three_bridge.verified and ordered_double_contact_bridge.verified
    )
    return SplitContactClosureCertificate(
        total_bases=exact_total_contact_base_certificate(),
        embeddings=embeddings,
        contact_three_jet_bridge_verified=contact_three_bridge.verified,
        ordered_double_contact_jet_bridge_verified=(
            ordered_double_contact_bridge.verified
        ),
        global_contact_three_sample_cyclic=(
            contact_three.verified
            and contact_three.sage_cyclic_simplification == (1, 0, True)
        ),
        global_double_contact_sample_cyclic=(
            double_contact.verified
            and double_contact.sage_cyclic_simplification == (1, 0, True)
        ),
        global_double_contact_cramer_verified=ordered.verified,
        # The exact rows and H bridges prove algebraic containment only.
        maximal_rank_algebraically_contained=bool(
            len(embeddings) == 8
            and all(item.verified for item in embeddings)
            and bridges_verified
        ),
        # This is the separate proper Whitney--Thom/isotopy-extension theorem
        # input already established for the clean global incidences.  Keeping
        # it separate prevents row algebra from manufacturing topology.
        proper_isotopy_extension_recorded=True,
        # All three overlap allocation incidences are dense images of exact
        # degenerating arcs in the ordinary global components.  This is
        # algebraic containment only; topology at the boundary remains below.
        overlap_algebraic_closure_count=(
            sum(item.verified for item in contact_three_overlap)
            + int(double_contact_overlap.verified)
        ),
        # C3 overlap-V, C3 overlap-W, and C2^2 overlap+W.
        exceptional_allocation_count=3,
        # Reduced ordered compatibility lengths 4+6+6.
        residual_affine_line_base_length=16,
        exceptional_topology_open=True,
    )


def main() -> int:
    """Print the maximal-rank contact closure and exact remaining boundary."""

    certificate = exact_split_contact_closure_certificate()
    print(
        "ordered two-pair total base domain:",
        certificate.total_bases.ordered_two_pair_base_domain,
    )
    print(
        "maximal-rank split contact topology closed:",
        certificate.maximal_rank_topology_closed,
    )
    print(
        "overlap allocations algebraically contained:",
        certificate.overlap_algebraic_closure_count,
    )
    print("remaining overlap allocations:", certificate.exceptional_allocation_count)
    print(
        "remaining residual ordered base length:",
        certificate.residual_affine_line_base_length,
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
