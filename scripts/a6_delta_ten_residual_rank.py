"""Close the two residual rank divisors in the delta-ten audit.

Two determinant factors were left open by the dominant-chart calculations:

* ``k*u*v-u-v`` for the two-triple profile ``T111^2 + 4N``; and
* the residual factor of the two-contact profile ``C2^2 + 6N``.

For the two-triple profile, solve the residual equation for ``k`` and
localize away from the two cusp fibers, coincident omitted roots, and the
same-``P``-fiber boundary.  The four normalized augmented determinants then
generate the unit ideal after Rabinowitsch localization.  Thus there is no
compatible point on this valid residual chart.  The explicit compatible
family ``k=0, v=-u`` is retained as a hostile fixture: it lies entirely on
the excluded same-fiber boundary and shows why the localization is needed.

For the two-contact profile, replace the ordered roots by
``sigma=u+v, pi=u*v``.  Divisibility of the collision decic by
``(s^2-sigma*s+pi)^2`` gives four affine-linear equations in ``(a,b,c,d)``.
On the residual determinant hypersurface:

* one augmented determinant is coprime to the residual factor, so the
  compatible rank-three base has dimension at most one;
* a rank-two point is singular on the residual hypersurface.  Exact
  projection resultants show that every positive-dimensional singular locus
  lies over ``pi=0``, ``pi=sigma-1``, ``pi=-sigma-1``, or
  ``sigma^2=4*pi``.  Restriction to these curves respectively gives no curve,
  ``k=-2``, ``k=2``, or the coincident-root/cusp boundary.  Hence the valid
  rank-two base is finite; and
* the coefficient-function Wronskians rule out rank one at every valid root.

The largest possible valid incidence dimension is therefore two, not three.
There are compatible rank-three points on the enlarged nonsplit, cusp-free
determinant chart.  The stored hostile member also has the extra-critical
factor zero, so it is not asserted to realize the clean ``C2^2 + 6N``
singularity profile.

This closes only these residual determinant factors.  It does not classify
the split fibers, replace the topology/propagation inputs, or prove the plane
Jacobian conjecture.  Every calculation below is exact over ``QQ`` and uses
only SymPy.
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
    Symbol,
    cancel,
    diff,
    expand,
    factor,
    fraction,
    gcd,
    groebner,
    rem,
    resultant,
)
from sympy.polys.matrices import DomainMatrix

from scripts.a6_delta_ten_double_triple import (
    COEFFICIENTS as TRIPLE_COEFFICIENTS,
    DOUBLE_TRIPLE_EQUATIONS,
    FIRST_CUSP_FIBER_FACTOR,
    FIRST_OMITTED_ROOT,
    INCIDENCE_AUGMENTED_MATRIX as TRIPLE_AUGMENTED_MATRIX,
    INCIDENCE_MATRIX as TRIPLE_MATRIX,
    INCIDENCE_RIGHT_HAND_SIDE as TRIPLE_RIGHT_HAND_SIDE,
    RESIDUAL_RANK_FACTOR as TRIPLE_RESIDUAL_FACTOR,
    SAME_FIBER_FACTOR,
    SECOND_CUSP_FIBER_FACTOR,
    SECOND_OMITTED_ROOT,
)
from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    GAMMA,
    KAPPA,
    PAIR_QUADRATIC,
    S,
)

LOCALIZATION_VARIABLE: Final = Symbol("z_residual_rank")
ROOT_SUM: Final = Symbol("sigma")
ROOT_PRODUCT: Final = Symbol("pi")
CONTACT_COEFFICIENTS: Final = (ALPHA, BETA, GAMMA, DELTA)
ZERO_CONTACT_COEFFICIENTS: Final = {
    coefficient: 0 for coefficient in CONTACT_COEFFICIENTS
}


def _domain_determinant(matrix: Matrix) -> Expr:
    """Return a polynomial determinant without heuristic simplification."""

    determinant: Expr = DomainMatrix.from_Matrix(matrix).det().as_expr()
    return determinant


def _primitive_part(expression: Expr, *generators: Symbol) -> Expr:
    """Remove the nonzero rational content of a polynomial expression."""

    _, primitive = Poly(expression, *generators).primitive()
    result: Expr = primitive.as_expr()
    return result


def _common_gcd(expressions: tuple[Expr, ...]) -> Expr:
    """Return the gcd of a nonempty tuple of exact expressions."""

    if not expressions:
        msg = "a common gcd needs at least one expression"
        raise ValueError(msg)
    common = expressions[0]
    for expression in expressions[1:]:
        common = gcd(common, expression)
    return common


# ---------------------------------------------------------------------------
# The residual divisor for T111^2 + 4N
# ---------------------------------------------------------------------------

TRIPLE_RESIDUAL_KAPPA: Final = cancel(
    (FIRST_OMITTED_ROOT + SECOND_OMITTED_ROOT)
    / (FIRST_OMITTED_ROOT * SECOND_OMITTED_ROOT)
)


@cache
def normalized_double_triple_compatibility() -> tuple[tuple[Expr, ...], Expr]:
    """Return normalized augmented minors and the valid-chart localizer."""

    u = FIRST_OMITTED_ROOT
    v = SECOND_OMITTED_ROOT
    substitution = {KAPPA: TRIPLE_RESIDUAL_KAPPA}
    first_cusp = cancel(FIRST_CUSP_FIBER_FACTOR.subs(substitution))
    second_cusp = cancel(SECOND_CUSP_FIBER_FACTOR.subs(substitution))
    same_fiber = cancel(SAME_FIBER_FACTOR.subs(substitution))
    invertible_common_factor = cancel(
        (u - v) ** 2 * first_cusp**2 * second_cusp**2 * same_fiber
    )

    normalized: list[Expr] = []
    for column in range(len(TRIPLE_COEFFICIENTS)):
        replaced = TRIPLE_MATRIX.copy()
        replaced[:, column] = TRIPLE_RIGHT_HAND_SIDE
        determinant = _domain_determinant(replaced)
        numerator, _ = fraction(
            cancel(determinant.subs(substitution) / invertible_common_factor)
        )
        normalized.append(_primitive_part(numerator, u, v))

    first_cusp_numerator, _ = fraction(first_cusp)
    second_cusp_numerator, _ = fraction(second_cusp)
    same_fiber_numerator, _ = fraction(same_fiber)
    localizer = expand(
        u
        * v
        * (u - v)
        * first_cusp_numerator
        * second_cusp_numerator
        * same_fiber_numerator
    )
    return tuple(normalized), localizer


@dataclass(frozen=True, slots=True)
class DoubleTripleResidualCertificate:
    """Exact saturation and hostile fixtures for ``T111^2 + 4N``."""

    residual_substitution_identity: Expr
    normalized_compatibility_count: int
    localized_groebner_basis: tuple[Expr, ...]
    generic_residual_value: Expr
    generic_localizer_value: Expr
    generic_coefficient_rank: int
    generic_augmented_rank: int
    excluded_family_equation_identities: tuple[Expr, ...]
    excluded_family_same_fiber_identity: Expr
    hostile_equation_residuals: tuple[Expr, ...]
    hostile_residual_value: Expr
    hostile_same_fiber_value: Expr
    hostile_cusp_values: tuple[Expr, Expr]
    hostile_root_difference: Expr
    hostile_coefficient_rank: int
    hostile_augmented_rank: int

    @property
    def verified(self) -> bool:
        """Whether the valid residual chart is empty and fixtures agree."""

        return bool(
            self.residual_substitution_identity == 0
            and self.normalized_compatibility_count == 4
            and self.localized_groebner_basis == (1,)
            and self.generic_residual_value == 0
            and self.generic_localizer_value != 0
            and self.generic_coefficient_rank == 3
            and self.generic_augmented_rank == 4
            and self.excluded_family_equation_identities == (0, 0, 0, 0)
            and self.excluded_family_same_fiber_identity == 0
            and self.hostile_equation_residuals == (0, 0, 0, 0)
            and self.hostile_residual_value == 0
            and self.hostile_same_fiber_value == 0
            and self.hostile_cusp_values == (2, 2)
            and self.hostile_root_difference == 2
            and self.hostile_coefficient_rank == 3
            and self.hostile_augmented_rank == 3
        )


@cache
def exact_double_triple_residual_certificate() -> DoubleTripleResidualCertificate:
    """Prove that the valid two-triple residual divisor is incompatible."""

    u = FIRST_OMITTED_ROOT
    v = SECOND_OMITTED_ROOT
    normalized, localizer = normalized_double_triple_compatibility()
    localized_basis = groebner(
        (
            *normalized,
            1 - LOCALIZATION_VARIABLE * localizer,
        ),
        LOCALIZATION_VARIABLE,
        u,
        v,
        order="grevlex",
    )

    generic_point = {
        KAPPA: Rational(3, 2),
        u: 1,
        v: 2,
    }

    # A whole compatible family before saturation.  Its two omitted roots
    # have the same P-value, so it is an ordinary-quadruple boundary rather
    # than two separate triple target fibers.
    free_delta = Symbol("d_same_fiber")
    excluded_family = {
        KAPPA: 0,
        v: -u,
        ALPHA: -(u**4) - u**2,
        BETA: cancel(free_delta * (2 * u**4 + 2 * u**2 + 1) / (u**4 + u**2 + 1)),
        GAMMA: 1,
        DELTA: free_delta,
    }

    hostile_point = {
        KAPPA: 0,
        u: 1,
        v: -1,
    }
    hostile_member = {
        **hostile_point,
        ALPHA: -2,
        BETA: 0,
        GAMMA: 1,
        DELTA: 0,
    }
    hostile_matrix = TRIPLE_MATRIX.subs(hostile_point)
    hostile_augmented = TRIPLE_AUGMENTED_MATRIX.subs(hostile_point)

    return DoubleTripleResidualCertificate(
        residual_substitution_identity=cancel(
            TRIPLE_RESIDUAL_FACTOR.subs(KAPPA, TRIPLE_RESIDUAL_KAPPA)
        ),
        normalized_compatibility_count=len(normalized),
        localized_groebner_basis=tuple(
            polynomial.as_expr() for polynomial in localized_basis.polys
        ),
        generic_residual_value=TRIPLE_RESIDUAL_FACTOR.subs(generic_point),
        generic_localizer_value=localizer.subs(
            {
                u: generic_point[u],
                v: generic_point[v],
            }
        ),
        generic_coefficient_rank=TRIPLE_MATRIX.subs(generic_point).rank(),
        generic_augmented_rank=TRIPLE_AUGMENTED_MATRIX.subs(generic_point).rank(),
        excluded_family_equation_identities=tuple(
            cancel(equation.subs(excluded_family))
            for equation in DOUBLE_TRIPLE_EQUATIONS
        ),
        excluded_family_same_fiber_identity=expand(
            SAME_FIBER_FACTOR.subs({KAPPA: 0, v: -u})
        ),
        hostile_equation_residuals=tuple(
            expand(equation.subs(hostile_member))
            for equation in DOUBLE_TRIPLE_EQUATIONS
        ),
        hostile_residual_value=TRIPLE_RESIDUAL_FACTOR.subs(hostile_point),
        hostile_same_fiber_value=SAME_FIBER_FACTOR.subs(hostile_point),
        hostile_cusp_values=(
            FIRST_CUSP_FIBER_FACTOR.subs(hostile_point),
            SECOND_CUSP_FIBER_FACTOR.subs(hostile_point),
        ),
        hostile_root_difference=(u - v).subs(hostile_point),
        hostile_coefficient_rank=hostile_matrix.rank(),
        hostile_augmented_rank=hostile_augmented.rank(),
    )


# ---------------------------------------------------------------------------
# The residual divisor for C2^2 + 6N
# ---------------------------------------------------------------------------

CONTACT_ROOT_QUARTIC: Final = (S**2 - ROOT_SUM * S + ROOT_PRODUCT) ** 2
CONTACT_ROOT_DISCRIMINANT: Final = ROOT_SUM**2 - 4 * ROOT_PRODUCT
CONTACT_SPLIT_JET_PRODUCT: Final = KAPPA**2 + 2 * KAPPA * ROOT_SUM + 4 * ROOT_PRODUCT
CONTACT_CUSP_PRODUCT: Final = (
    KAPPA**2 * ROOT_PRODUCT
    + KAPPA * ROOT_PRODUCT * ROOT_SUM
    + KAPPA * ROOT_SUM
    + ROOT_PRODUCT**2
    - 2 * ROOT_PRODUCT
    + ROOT_SUM**2
    + 1
)
CONTACT_GLOBAL_SPLIT_FACTOR: Final = KAPPA * (KAPPA**2 - 4)
CONTACT_BASE_LOCALIZER: Final = expand(
    CONTACT_ROOT_DISCRIMINANT
    * CONTACT_SPLIT_JET_PRODUCT
    * CONTACT_CUSP_PRODUCT
    * CONTACT_GLOBAL_SPLIT_FACTOR
)


@cache
def unordered_contact_incidence_system() -> tuple[tuple[Expr, ...], Matrix, Matrix]:
    """Return the four remainder equations in ``(k, sigma, pi)``."""

    remainder = expand(rem(COLLISION_POLYNOMIAL, CONTACT_ROOT_QUARTIC, S))
    equations = tuple(
        Poly(remainder, S).coeff_monomial(S**degree) for degree in range(4)
    )
    coefficient_matrix = Matrix(
        [
            [equation.coeff(parameter) for parameter in CONTACT_COEFFICIENTS]
            for equation in equations
        ]
    )
    right_side = Matrix(
        [-equation.subs(ZERO_CONTACT_COEFFICIENTS) for equation in equations]
    )
    return equations, coefficient_matrix, right_side


@cache
def contact_residual_factor() -> tuple[Expr, Expr]:
    """Return the residual determinant factor and exact division remainder."""

    _, coefficient_matrix, _ = unordered_contact_incidence_system()
    determinant = _domain_determinant(coefficient_matrix)
    visible = CONTACT_SPLIT_JET_PRODUCT**2 * CONTACT_CUSP_PRODUCT
    quotient, remainder_polynomial = Poly(
        determinant,
        KAPPA,
        ROOT_SUM,
        ROOT_PRODUCT,
    ).div(
        Poly(
            visible,
            KAPPA,
            ROOT_SUM,
            ROOT_PRODUCT,
        )
    )
    return quotient.as_expr(), remainder_polynomial.as_expr()


CONTACT_PROJECTION_CURVE_FACTOR: Final = (
    ROOT_PRODUCT
    * (4 * ROOT_PRODUCT - ROOT_SUM**2) ** 2
    * (ROOT_PRODUCT - ROOT_SUM + 1) ** 12
    * (ROOT_PRODUCT + ROOT_SUM + 1) ** 12
)
CONTACT_JET_VISIBLE_FACTOR: Final = (KAPPA + 2 * S) ** 2 * PAIR_QUADRATIC


@cache
def normalized_contact_jet_wronskians() -> tuple[tuple[Expr, ...], tuple[Expr, ...]]:
    """Divide the six coefficient-function Wronskians by visible factors."""

    coefficient_functions = tuple(
        COLLISION_POLYNOMIAL.coeff(parameter) for parameter in CONTACT_COEFFICIENTS
    )
    normalized: list[Expr] = []
    remainders: list[Expr] = []
    for first, second in combinations(range(len(coefficient_functions)), 2):
        wronskian = expand(
            coefficient_functions[first] * diff(coefficient_functions[second], S)
            - coefficient_functions[second] * diff(coefficient_functions[first], S)
        )
        quotient, remainder_polynomial = Poly(wronskian, KAPPA, S).div(
            Poly(CONTACT_JET_VISIBLE_FACTOR, KAPPA, S)
        )
        normalized.append(quotient.as_expr())
        remainders.append(remainder_polynomial.as_expr())
    return tuple(normalized), tuple(remainders)


@dataclass(frozen=True, slots=True)
class DoubleContactResidualCertificate:
    """Exact dimension bounds for the ``C2^2 + 6N`` residual divisor."""

    determinant_factor_remainder: Expr
    residual_total_degree: int
    residual_term_count: int
    augmented_minor_gcd: Expr
    projection_curve_factor_identity: Expr
    vertical_projection_groebner_basis: tuple[Expr, ...]
    restricted_singular_curve_gcds: tuple[Expr, Expr, Expr, Expr]
    split_curve_boundary_identities: tuple[Expr, Expr, Expr, Expr]
    diagonal_cusp_identity: Expr
    jet_division_remainders: tuple[Expr, ...]
    jet_split_power_remainder: Expr
    incompatible_residual_value: Expr
    incompatible_base_localizer_value: Expr
    incompatible_coefficient_rank: int
    incompatible_augmented_rank: int
    compatible_residual_value: Expr
    compatible_base_localizer_value: Expr
    compatible_coefficient_rank: int
    compatible_augmented_rank: int
    compatible_equation_residuals: tuple[Expr, ...]
    compatible_extra_critical_value: Expr
    maximum_valid_incidence_dimension: int

    @property
    def verified(self) -> bool:
        """Whether no valid residual component can have dimension three."""

        return bool(
            self.determinant_factor_remainder == 0
            and self.residual_total_degree == 11
            and self.residual_term_count == 79
            and self.augmented_minor_gcd == 1
            and self.projection_curve_factor_identity == 0
            and self.vertical_projection_groebner_basis == (1,)
            and self.restricted_singular_curve_gcds
            == (
                1,
                (KAPPA + 2) ** 3,
                (KAPPA - 2) ** 3,
                2 * KAPPA * ROOT_SUM + ROOT_SUM**2 + 4,
            )
            and self.split_curve_boundary_identities == (0, 0, 0, 0)
            and self.diagonal_cusp_identity == 0
            and self.jet_division_remainders == (0, 0, 0, 0, 0, 0)
            and self.jet_split_power_remainder == 0
            and self.incompatible_residual_value == 0
            and self.incompatible_base_localizer_value != 0
            and self.incompatible_coefficient_rank == 3
            and self.incompatible_augmented_rank == 4
            and self.compatible_residual_value == 0
            and self.compatible_base_localizer_value != 0
            and self.compatible_coefficient_rank == 3
            and self.compatible_augmented_rank == 3
            and self.compatible_equation_residuals == (0, 0, 0, 0)
            and self.compatible_extra_critical_value == 0
            and self.maximum_valid_incidence_dimension == 2
        )


@cache
def exact_double_contact_residual_certificate() -> DoubleContactResidualCertificate:
    """Bound every valid two-contact residual incidence component by two."""

    equations, coefficient_matrix, right_side = unordered_contact_incidence_system()
    residual, determinant_remainder = contact_residual_factor()

    replaced = coefficient_matrix.copy()
    replaced[:, 0] = right_side
    augmented_minor = _domain_determinant(replaced)
    augmented_gcd = gcd(residual, augmented_minor)

    residual_derivatives = (
        diff(residual, KAPPA),
        diff(residual, ROOT_SUM),
        diff(residual, ROOT_PRODUCT),
    )
    kappa_projection_resultants = tuple(
        resultant(residual, derivative, KAPPA)
        for derivative in residual_derivatives[:2]
    )
    projection_gcd = _primitive_part(
        gcd(*kappa_projection_resultants),
        ROOT_SUM,
        ROOT_PRODUCT,
    )

    kappa_coefficients = tuple(Poly(residual, KAPPA).all_coeffs())
    vertical_projection_basis = groebner(
        kappa_coefficients,
        ROOT_SUM,
        ROOT_PRODUCT,
        order="grevlex",
    )

    singular_equations = (residual, *residual_derivatives)
    curve_substitutions = (
        {ROOT_PRODUCT: 0},
        {ROOT_PRODUCT: ROOT_SUM - 1},
        {ROOT_PRODUCT: -ROOT_SUM - 1},
        {ROOT_PRODUCT: ROOT_SUM**2 / 4},
    )
    restricted_gcds = tuple(
        factor(
            _primitive_part(
                _common_gcd(
                    tuple(
                        expand(equation.subs(substitution))
                        for equation in singular_equations
                    )
                ),
                KAPPA,
                ROOT_SUM,
            )
        )
        for substitution in curve_substitutions
    )

    normalized_wronskians, jet_remainders = normalized_contact_jet_wronskians()
    jet_basis = groebner(
        normalized_wronskians,
        KAPPA,
        S,
        order="grevlex",
    )
    _, jet_split_power_remainder = jet_basis.reduce((KAPPA + 2 * S) ** 4)

    incompatible_point = {
        KAPPA: Rational(-10, 3),
        ROOT_SUM: Rational(10, 3),
        ROOT_PRODUCT: Rational(59, 27),
    }
    incompatible_matrix = coefficient_matrix.subs(incompatible_point)
    incompatible_augmented = incompatible_matrix.row_join(
        right_side.subs(incompatible_point)
    )

    compatible_point = {
        KAPPA: -3,
        ROOT_SUM: 3,
        ROOT_PRODUCT: 3,
    }
    compatible_member = {
        **compatible_point,
        ALPHA: -20,
        BETA: 43,
        GAMMA: -20,
        DELTA: 0,
    }
    compatible_matrix = coefficient_matrix.subs(compatible_point)
    compatible_augmented = compatible_matrix.row_join(right_side.subs(compatible_point))

    diagonal_singular_factor = 2 * KAPPA * ROOT_SUM + ROOT_SUM**2 + 4
    return DoubleContactResidualCertificate(
        determinant_factor_remainder=determinant_remainder,
        residual_total_degree=Poly(
            residual,
            KAPPA,
            ROOT_SUM,
            ROOT_PRODUCT,
        ).total_degree(),
        residual_term_count=len(Poly(residual, KAPPA, ROOT_SUM, ROOT_PRODUCT).terms()),
        augmented_minor_gcd=augmented_gcd,
        projection_curve_factor_identity=expand(
            projection_gcd - CONTACT_PROJECTION_CURVE_FACTOR
        ),
        vertical_projection_groebner_basis=tuple(
            polynomial.as_expr() for polynomial in vertical_projection_basis.polys
        ),
        restricted_singular_curve_gcds=restricted_gcds,
        split_curve_boundary_identities=(
            expand(
                CONTACT_SPLIT_JET_PRODUCT.subs({ROOT_PRODUCT: ROOT_SUM - 1, KAPPA: -2})
            ),
            expand(CONTACT_CUSP_PRODUCT.subs({ROOT_PRODUCT: ROOT_SUM - 1, KAPPA: -2})),
            expand(
                CONTACT_SPLIT_JET_PRODUCT.subs({ROOT_PRODUCT: -ROOT_SUM - 1, KAPPA: 2})
            ),
            expand(CONTACT_CUSP_PRODUCT.subs({ROOT_PRODUCT: -ROOT_SUM - 1, KAPPA: 2})),
        ),
        diagonal_cusp_identity=expand(
            CONTACT_CUSP_PRODUCT.subs(
                ROOT_PRODUCT,
                ROOT_SUM**2 / 4,
            )
            - diagonal_singular_factor**2 / 16
        ),
        jet_division_remainders=jet_remainders,
        jet_split_power_remainder=jet_split_power_remainder,
        incompatible_residual_value=residual.subs(incompatible_point),
        incompatible_base_localizer_value=CONTACT_BASE_LOCALIZER.subs(
            incompatible_point
        ),
        incompatible_coefficient_rank=incompatible_matrix.rank(),
        incompatible_augmented_rank=incompatible_augmented.rank(),
        compatible_residual_value=residual.subs(compatible_point),
        compatible_base_localizer_value=CONTACT_BASE_LOCALIZER.subs(compatible_point),
        compatible_coefficient_rank=compatible_matrix.rank(),
        compatible_augmented_rank=compatible_augmented.rank(),
        compatible_equation_residuals=tuple(
            expand(equation.subs(compatible_member)) for equation in equations
        ),
        compatible_extra_critical_value=EXTRA_CRITICAL_FACTOR.subs(compatible_member),
        # rank 3: base <= 1 plus fiber 1; rank 2: finite base plus fiber 2.
        # Rank 1 is excluded by the localized Wronskian calculation.
        maximum_valid_incidence_dimension=2,
    )


@dataclass(frozen=True, slots=True)
class DeltaTenResidualRankCertificate:
    """Combined exact closure of both residual determinant factors."""

    double_triple: DoubleTripleResidualCertificate
    double_contact: DoubleContactResidualCertificate

    @property
    def verified(self) -> bool:
        """Whether both residual rank-divisor arguments are exact."""

        return self.double_triple.verified and self.double_contact.verified


@cache
def exact_delta_ten_residual_rank_certificate() -> DeltaTenResidualRankCertificate:
    """Build the combined exact residual-rank certificate."""

    return DeltaTenResidualRankCertificate(
        double_triple=exact_double_triple_residual_certificate(),
        double_contact=exact_double_contact_residual_certificate(),
    )


def main() -> int:
    """Print a compact summary of the exact residual-rank closure."""

    certificate = exact_delta_ten_residual_rank_certificate()
    print("delta-ten residual-rank certificate:", certificate.verified)
    print(
        "two-triple localized basis:",
        certificate.double_triple.localized_groebner_basis,
    )
    print(
        "two-contact projected singular curves:",
        factor(CONTACT_PROJECTION_CURVE_FACTOR),
    )
    print(
        "two-contact maximum valid incidence dimension:",
        certificate.double_contact.maximum_valid_incidence_dimension,
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
