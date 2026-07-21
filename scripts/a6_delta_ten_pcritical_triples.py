"""Classify the ``P``-critical triple-fiber boundaries at delta ten.

For ``P=t^2+k*t^3+t^4``, let ``e`` be a nonzero critical source.  Then

``k = -(4*e^2+2)/(3*e)``

and the fiber through ``e`` factors as

``P(t)-P(e) = (t-e)^2 * (t^2-s*t+p)``,

where ``s=2(1-e^2)/(3e)`` and ``p=(e^2-1)/3``.  Thus the three distinct
sources are ``e`` and the two roots of the displayed quadratic; the fourth
root counted with multiplicity is ``e`` again.  On the valid chart

``e*(e^2-1)*(2e^2-1)*(2e^2+1) != 0``.

The removed labeled fourth-root boundary of the ordinary triple surface is
the union of three prime curves, according to which label is ``e``.  A
critical branch of the plane parametrization has derivative vector
``(P'(e),Q'(e))=(0,Q'(e))``.  It is therefore a valid immersed vertical
branch exactly when ``Q'(e) != 0``; ``P'(e)=0`` alone is not an invalidity.

For ``T112``, the critical branch cannot belong to the tangent pair on this
immersed open: the cross-multiplied tangent equation then reduces exactly to
``Q'(e)=0``.  The two noncritical branches may be tangent.  Their three
linear conditions have rank three everywhere on the valid critical base, so
the labeled incidence is an affine-line bundle over a curve and has
dimension two.  It has no three-dimensional same-profile component.

For ``C2 + T111``, the two triple-equality equations and the two separate
contact equations have full rank generically over the ``(e,w)`` base.  The
rank-at-most-two minors have only the pair-denominator divisor as a common
curve; away from it their common zero set is finite, and the triple rows
always have rank two.  The denominator divisor carries no actual pair except
at the split values ``e=+/-1/2``.  The true split vertical charts are audited
separately and are at most one-dimensional.  Consequently no
three-dimensional same-profile component occurs here either; the valid
critical incidence has dimension two.

Two hostile rational fixtures prove that these valid opens are nonempty:

* ``k=-3`` and ``Q=-17/4*t^5+11/4*t^6-17/4*t^7+t^9`` gives
  ``T112 + 6N`` with the critical source ``e=2`` immersed; and
* ``k=-3`` and ``Q=-173/10*t^5+361/10*t^6-173/10*t^7+t^9`` gives
  ``C2 + T111 + 5N`` with the same immersed critical branch.

Sage independently recomputes the labeled-boundary saturation, primality,
and both hostile Jacobian schemes.  No van Kamp calculation is made here.
Topology propagation remains a separate gap requiring a connected clean
open and proper projective Whitney--Thom triviality.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from itertools import combinations
from typing import Final

from sympy import (
    Expr,
    I,
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
    groebner,
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
    DELTA,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    R,
    S,
    T,
    TANGENCY_POLYNOMIAL,
    X,
    exceptional_vertical,
)


def _polynomial_gcd(polynomials: tuple[Expr, ...]) -> Expr:
    """Return the iterated gcd of a nonempty polynomial tuple."""

    if not polynomials:
        msg = "a polynomial gcd needs at least one input"
        raise ValueError(msg)
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = gcd(common, polynomial)
    return factor(common)


def _maximal_minors(matrix: Matrix) -> tuple[Expr, ...]:
    """Return all maximal minors of a matrix with at most as many rows."""

    if matrix.rows > matrix.cols:
        msg = "maximal-minor helper expects rows <= columns"
        raise ValueError(msg)
    return tuple(
        expand(matrix[:, columns].det(method="domain-ge"))
        for columns in combinations(range(matrix.cols), matrix.rows)
    )


def _all_minors(matrix: Matrix, size: int) -> tuple[Expr, ...]:
    """Return all square minors of the requested size."""

    return tuple(
        expand(matrix.extract(rows, columns).det(method="domain-ge"))
        for rows in combinations(range(matrix.rows), size)
        for columns in combinations(range(matrix.cols), size)
    )


def _modulo_branch(rational_expression: Expr, branch: Expr) -> Expr:
    """Reduce a rational function modulo a monic source polynomial."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    denominator_inverse = invert(denominator, branch)
    return expand(rem(numerator * denominator_inverse, branch, T))


def _modulo_exceptional_quadratic(rational_expression: Expr) -> Expr:
    """Reduce a rational function of ``e`` modulo ``e^2+2``."""

    numerator, denominator = together(rational_expression).as_numer_denom()
    reduced_numerator = rem(numerator, CRITICAL_ROOT**2 + 2, CRITICAL_ROOT)
    reduced_denominator = rem(denominator, CRITICAL_ROOT**2 + 2, CRITICAL_ROOT)
    denominator_inverse = invert(
        reduced_denominator,
        CRITICAL_ROOT**2 + 2,
    )
    return expand(
        rem(
            reduced_numerator * denominator_inverse,
            CRITICAL_ROOT**2 + 2,
            CRITICAL_ROOT,
        )
    )


LINEAR_PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)
ZERO_LINEAR_PARAMETERS: Final = dict.fromkeys(LINEAR_PARAMETERS, 0)
CRITICAL_ROOT: Final = Symbol("e_pcritical")
CONTACT_SUM: Final = Symbol("w_pcritical")
SPLIT_PAIR_PRODUCT: Final = Symbol("rho_pcritical")
FREE_DELTA: Final = Symbol("d_pcritical")

# Labeled triple-surface boundary.
ROOT_A: Final = Symbol("u_pcritical")
ROOT_B: Final = Symbol("v_pcritical")
ROOT_C: Final = Symbol("z_pcritical")
ROOT_SIGMA_1: Final = ROOT_A + ROOT_B + ROOT_C
ROOT_SIGMA_2: Final = ROOT_A * ROOT_B + ROOT_A * ROOT_C + ROOT_B * ROOT_C
ROOT_SIGMA_3: Final = ROOT_A * ROOT_B * ROOT_C
ROOT_BASE_CONSTRAINT: Final = expand(
    ROOT_SIGMA_2**2 - ROOT_SIGMA_1 * ROOT_SIGMA_3 - ROOT_SIGMA_2
)
ROOT_FOURTH_FACTORS: Final = (
    ROOT_SIGMA_2 + ROOT_B * ROOT_C,
    ROOT_SIGMA_2 + ROOT_A * ROOT_C,
    ROOT_SIGMA_2 + ROOT_A * ROOT_B,
)
ROOT_CRITICAL_COMPONENTS: Final = (
    (
        3 * ROOT_A * (ROOT_B + ROOT_C) + 2 * ROOT_A**2 - 2,
        3 * ROOT_B * ROOT_C - ROOT_A**2 + 1,
    ),
    (
        3 * ROOT_B * (ROOT_A + ROOT_C) + 2 * ROOT_B**2 - 2,
        3 * ROOT_A * ROOT_C - ROOT_B**2 + 1,
    ),
    (
        3 * ROOT_C * (ROOT_A + ROOT_B) + 2 * ROOT_C**2 - 2,
        3 * ROOT_A * ROOT_B - ROOT_C**2 + 1,
    ),
)

# One critical component, parameterized by its critical source e.
CRITICAL_KAPPA: Final = -(4 * CRITICAL_ROOT**2 + 2) / (3 * CRITICAL_ROOT)
OTHER_ROOT_SUM: Final = 2 * (1 - CRITICAL_ROOT**2) / (3 * CRITICAL_ROOT)
OTHER_ROOT_PRODUCT: Final = (CRITICAL_ROOT**2 - 1) / 3
OTHER_ROOT_QUADRATIC: Final = T**2 - OTHER_ROOT_SUM * T + OTHER_ROOT_PRODUCT
CRITICAL_TRIPLE_CUBIC: Final = (T - CRITICAL_ROOT) * OTHER_ROOT_QUADRATIC
CRITICAL_P: Final = FAMILY_P.subs(KAPPA, CRITICAL_KAPPA)
CRITICAL_TARGET_X: Final = -(
    CRITICAL_ROOT**2 * (CRITICAL_ROOT - 1) * (CRITICAL_ROOT + 1) / 3
)
CRITICAL_P_FACTORIZATION_IDENTITY: Final = cancel(
    CRITICAL_P - CRITICAL_TARGET_X - (T - CRITICAL_ROOT) ** 2 * OTHER_ROOT_QUADRATIC
)
CRITICAL_BASE_LOCALIZER: Final = (
    CRITICAL_ROOT
    * (CRITICAL_ROOT**2 - 1)
    * (2 * CRITICAL_ROOT**2 - 1)
    * (2 * CRITICAL_ROOT**2 + 1)
)
OTHER_ROOT_DISCRIMINANT: Final = discriminant(OTHER_ROOT_QUADRATIC, T)
OTHER_P_DERIVATIVE_RESULTANT: Final = resultant(
    OTHER_ROOT_QUADRATIC,
    diff(CRITICAL_P, T),
    T,
)

# Equality of Q on e and the two roots of the quadratic.
CRITICAL_Q_REMAINDER: Final = expand(rem(FAMILY_Q, OTHER_ROOT_QUADRATIC, T))
CRITICAL_TRIPLE_EQUALITY: Final = (
    cancel(Poly(CRITICAL_Q_REMAINDER, T).coeff_monomial(T)),
    cancel(
        Poly(CRITICAL_Q_REMAINDER, T).coeff_monomial(1)
        - FAMILY_Q.subs(T, CRITICAL_ROOT)
    ),
)


def _divided_powers(maximum: int) -> dict[int, Expr]:
    """Return ``(u^n-v^n)/(u-v)`` in the other-root symmetric data."""

    values: dict[int, Expr] = {1: 1, 2: OTHER_ROOT_SUM}
    for degree in range(3, maximum + 1):
        values[degree] = expand(
            OTHER_ROOT_SUM * values[degree - 1]
            - OTHER_ROOT_PRODUCT * values[degree - 2]
        )
    return values


OTHER_DIVIDED_POWERS: Final = _divided_powers(7)
Q_DEGREE_COEFFICIENTS: Final = {
    5: ALPHA,
    6: BETA,
    7: GAMMA,
    8: DELTA,
    9: 1,
}
OTHER_PAIR_TANGENCY: Final = cancel(
    sum(
        degree
        * coefficient
        * (
            2 * OTHER_ROOT_PRODUCT * OTHER_DIVIDED_POWERS[degree - 2]
            + 3
            * CRITICAL_KAPPA
            * OTHER_ROOT_PRODUCT**2
            * OTHER_DIVIDED_POWERS[degree - 3]
            + 4 * OTHER_ROOT_PRODUCT**3 * OTHER_DIVIDED_POWERS[degree - 4]
        )
        for degree, coefficient in Q_DEGREE_COEFFICIENTS.items()
    )
)
T112_CRITICAL_EQUATIONS: Final = (
    *CRITICAL_TRIPLE_EQUALITY,
    OTHER_PAIR_TANGENCY,
)
T112_INCIDENCE_MATRIX: Final = Matrix(
    [
        [diff(equation, parameter) for parameter in LINEAR_PARAMETERS]
        for equation in T112_CRITICAL_EQUATIONS
    ]
)
T112_CLEARING_MULTIPLIERS: Final = (
    6561 * CRITICAL_ROOT**8,
    6561 * CRITICAL_ROOT**7,
    243 * CRITICAL_ROOT**6,
)
T112_CLEARED_MATRIX: Final = Matrix(
    [
        [
            cancel(multiplier * diff(equation, parameter))
            for parameter in LINEAR_PARAMETERS
        ]
        for equation, multiplier in zip(
            T112_CRITICAL_EQUATIONS,
            T112_CLEARING_MULTIPLIERS,
            strict=True,
        )
    ]
)
T112_EXPECTED_MINOR_GCD: Final = (
    708588
    * CRITICAL_ROOT**15
    * (CRITICAL_ROOT - 1) ** 5
    * (CRITICAL_ROOT + 1) ** 5
    * (2 * CRITICAL_ROOT**2 - 1) ** 2
    * (2 * CRITICAL_ROOT**2 + 1)
)

CRITICAL_Q_DERIVATIVE: Final = diff(FAMILY_Q, T).subs(T, CRITICAL_ROOT)
CRITICAL_TANGENT_REDUCTION_IDENTITY: Final = cancel(
    resultant(
        OTHER_ROOT_QUADRATIC,
        CRITICAL_Q_DERIVATIVE * diff(CRITICAL_P, T),
        T,
    )
    - CRITICAL_Q_DERIVATIVE**2 * OTHER_P_DERIVATIVE_RESULTANT
)
T112_IMMERSION_MATRIX: Final = T112_INCIDENCE_MATRIX.col_join(
    Matrix(
        [[diff(CRITICAL_Q_DERIVATIVE, parameter) for parameter in LINEAR_PARAMETERS]]
    )
)
EXPECTED_T112_IMMERSION_DETERMINANT: Final = (
    16
    * (CRITICAL_ROOT - 1) ** 5
    * (CRITICAL_ROOT + 1) ** 5
    * (CRITICAL_ROOT**2 + 2)
    * (2 * CRITICAL_ROOT**2 - 1) ** 4
    * (2 * CRITICAL_ROOT**2 + 1)
    / (2187 * CRITICAL_ROOT)
)
EXCEPTIONAL_IMMERSION_SOLUTION: Final = {
    ALPHA: -4 * FREE_DELTA * CRITICAL_ROOT - 15,
    BETA: -FREE_DELTA + 4 * CRITICAL_ROOT,
    GAMMA: -3 * FREE_DELTA * CRITICAL_ROOT - 8,
    DELTA: FREE_DELTA,
}

# C2 + T111: add a separate contact pair with sum w.
C2_CONTACT_EQUATIONS: Final = (
    cancel(COLLISION_POLYNOMIAL.subs({KAPPA: CRITICAL_KAPPA, S: CONTACT_SUM})),
    cancel(TANGENCY_POLYNOMIAL.subs({KAPPA: CRITICAL_KAPPA, S: CONTACT_SUM})),
)
C2_CRITICAL_EQUATIONS: Final = (
    *CRITICAL_TRIPLE_EQUALITY,
    *C2_CONTACT_EQUATIONS,
)
C2_INCIDENCE_MATRIX: Final = Matrix(
    [
        [diff(equation, parameter) for parameter in LINEAR_PARAMETERS]
        for equation in C2_CRITICAL_EQUATIONS
    ]
)
C2_CLEARING_MULTIPLIERS: Final = (
    2187 * CRITICAL_ROOT**7,
    2187 * CRITICAL_ROOT**6,
    81 * CRITICAL_ROOT**4,
    81 * CRITICAL_ROOT**5,
)
C2_CLEARED_MATRIX: Final = Matrix(
    [
        [
            cancel(multiplier * diff(equation, parameter))
            for parameter in LINEAR_PARAMETERS
        ]
        for equation, multiplier in zip(
            C2_CRITICAL_EQUATIONS,
            C2_CLEARING_MULTIPLIERS,
            strict=True,
        )
    ]
)
C2_PAIR_DENOMINATOR_FACTOR: Final = (
    -2 * CRITICAL_ROOT**2 + 3 * CRITICAL_ROOT * CONTACT_SUM - 1
)
C2_PAIR_DIAGONAL_FACTORS: Final = (
    CONTACT_SUM - 2 * CRITICAL_ROOT,
    CRITICAL_ROOT * CONTACT_SUM - 1,
)
C2_OVERLAP_FACTORS: Final = (
    2 * CRITICAL_ROOT**2 + 3 * CRITICAL_ROOT * CONTACT_SUM - 2,
    2 * CRITICAL_ROOT**3
    - 4 * CRITICAL_ROOT**2 * CONTACT_SUM
    + 3 * CRITICAL_ROOT * CONTACT_SUM**2
    + CRITICAL_ROOT
    - 2 * CONTACT_SUM,
)
C2_RESIDUAL_RANK_FACTOR: Final = (
    -32 * CRITICAL_ROOT**7 * CONTACT_SUM
    + 184 * CRITICAL_ROOT**6 * CONTACT_SUM**2
    + 72 * CRITICAL_ROOT**6
    - 168 * CRITICAL_ROOT**5 * CONTACT_SUM**3
    - 528 * CRITICAL_ROOT**5 * CONTACT_SUM
    + 36 * CRITICAL_ROOT**4 * CONTACT_SUM**4
    + 714 * CRITICAL_ROOT**4 * CONTACT_SUM**2
    + 234 * CRITICAL_ROOT**4
    - 312 * CRITICAL_ROOT**3 * CONTACT_SUM**3
    - 540 * CRITICAL_ROOT**3 * CONTACT_SUM
    + 45 * CRITICAL_ROOT**2 * CONTACT_SUM**4
    + 324 * CRITICAL_ROOT**2 * CONTACT_SUM**2
    + 99 * CRITICAL_ROOT**2
    - 60 * CRITICAL_ROOT * CONTACT_SUM**3
    - 88 * CRITICAL_ROOT * CONTACT_SUM
    + 20 * CONTACT_SUM**2
)
EXPECTED_C2_INCIDENCE_DETERMINANT: Final = cancel(
    -32
    * CONTACT_SUM
    * C2_PAIR_DIAGONAL_FACTORS[0]
    * (CRITICAL_ROOT - 1) ** 2
    * (CRITICAL_ROOT + 1) ** 2
    * (2 * CRITICAL_ROOT**2 - 1)
    * C2_PAIR_DIAGONAL_FACTORS[1]
    * C2_PAIR_DENOMINATOR_FACTOR**3
    * C2_OVERLAP_FACTORS[0] ** 2
    * C2_OVERLAP_FACTORS[1] ** 2
    * C2_RESIDUAL_RANK_FACTOR
    / (531441 * CRITICAL_ROOT**10)
)

# On the pair-denominator divisor, the original pair equation still requires
# w*(w^2+k*w+1)=0.  Its only valid critical roots are e=+/-1/2.
PAIR_DENOMINATOR_CONTACT_SUM: Final = -CRITICAL_KAPPA / 2
EXPECTED_PAIR_INCIDENCE_RESIDUAL: Final = -(
    (CRITICAL_ROOT - 1)
    * (CRITICAL_ROOT + 1)
    * (2 * CRITICAL_ROOT - 1)
    * (2 * CRITICAL_ROOT + 1)
    * (2 * CRITICAL_ROOT**2 + 1)
    / (27 * CRITICAL_ROOT**3)
)

# True split vertical pair charts at e=-epsilon/2, k=2*epsilon.
SPLIT_SIGNS: Final = (-1, 1)
SPLIT_VERTICAL_EQUATIONS: Final = tuple(
    (
        *(
            equation.subs(CRITICAL_ROOT, Rational(-epsilon, 2))
            for equation in CRITICAL_TRIPLE_EQUALITY
        ),
        exceptional_vertical(epsilon).subs(R, SPLIT_PAIR_PRODUCT),
        diff(
            exceptional_vertical(epsilon).subs(R, SPLIT_PAIR_PRODUCT),
            SPLIT_PAIR_PRODUCT,
        ),
    )
    for epsilon in SPLIT_SIGNS
)
SPLIT_VERTICAL_MATRICES: Final = tuple(
    Matrix(
        [
            [diff(equation, parameter) for parameter in LINEAR_PARAMETERS]
            for equation in equations
        ]
    )
    for equations in SPLIT_VERTICAL_EQUATIONS
)
EXPECTED_SPLIT_VERTICAL_DETERMINANTS: Final = tuple(
    epsilon
    * (4 * SPLIT_PAIR_PRODUCT + 1) ** 2
    * (10 * SPLIT_PAIR_PRODUCT**2 - 6 * SPLIT_PAIR_PRODUCT + 1)
    / 2048
    for epsilon in SPLIT_SIGNS
)
SPLIT_VERTICAL_RANKDROP_ROOTS: Final = (
    Rational(-1, 4),
    Rational(3, 10) - I / 10,
    Rational(3, 10) + I / 10,
)

# Hostile valid T112 + 6N fixture.
PCRITICAL_T112_PARAMETERS: Final = {
    KAPPA: -3,
    ALPHA: Rational(-17, 4),
    BETA: Rational(11, 4),
    GAMMA: Rational(-17, 4),
    DELTA: 0,
}
PCRITICAL_T112_P: Final = expand(FAMILY_P.subs(PCRITICAL_T112_PARAMETERS))
PCRITICAL_T112_Q: Final = expand(FAMILY_Q.subs(PCRITICAL_T112_PARAMETERS))
PCRITICAL_T112_OTHER_BRANCH: Final = T**2 + T + 1
PCRITICAL_T112_TRIPLE_PAIR_SUMS: Final = (S + 1) * (S**2 - 3 * S + 3)
PCRITICAL_T112_RESIDUAL_NODES: Final = (S - 1) * (
    4 * S**5 - 64 * S**4 + 266 * S**3 - 289 * S**2 + 3 * S + 51
)
PCRITICAL_T112_COLLISION: Final = expand(
    (S + 1) ** 2 * (S**2 - 3 * S + 3) * PCRITICAL_T112_RESIDUAL_NODES / 4
)
PCRITICAL_T112_NODE_X: Final = (
    4096 * X**6
    + 4366592 * X**5
    - 1310683088 * X**4
    - 6007054657 * X**3
    + 14418241466 * X**2
    + 5849953792 * X
    - 485499872
)

# Hostile valid C2 + T111 + 5N fixture.
PCRITICAL_C2_PARAMETERS: Final = {
    KAPPA: -3,
    ALPHA: Rational(-173, 10),
    BETA: Rational(361, 10),
    GAMMA: Rational(-173, 10),
    DELTA: 0,
}
PCRITICAL_C2_P: Final = expand(FAMILY_P.subs(PCRITICAL_C2_PARAMETERS))
PCRITICAL_C2_Q: Final = expand(FAMILY_Q.subs(PCRITICAL_C2_PARAMETERS))
PCRITICAL_C2_CONTACT_BRANCH: Final = T**2 - T + 1
PCRITICAL_C2_OTHER_BRANCH: Final = T**2 + T + 1
PCRITICAL_C2_TRIPLE_PAIR_SUMS: Final = (S + 1) * (S**2 - 3 * S + 3)
PCRITICAL_C2_RESIDUAL_NODES: Final = (
    10 * S**5 - 140 * S**4 + 626 * S**3 - 1157 * S**2 + 1173 * S - 519
)
PCRITICAL_C2_COLLISION: Final = expand(
    (S - 1) ** 2 * PCRITICAL_C2_TRIPLE_PAIR_SUMS * PCRITICAL_C2_RESIDUAL_NODES / 10
)
PCRITICAL_C2_NODE_X: Final = (
    1000000 * X**5
    + 259090000 * X**4
    + 1749300900 * X**3
    + 1538402411 * X**2
    - 8338124136 * X
    - 23784186704
)


@dataclass(frozen=True, slots=True)
class DeltaTenPCriticalTripleCertificate:
    """Exact component algebra and hostile-fixture data."""

    critical_p_factorization_identity: Expr
    labeled_component_remainders: tuple[tuple[Expr, Expr], ...]
    labeled_boundary_saturation_exponent: int
    labeled_boundary_component_count: int
    labeled_boundary_dimensions: tuple[int, int, int]
    labeled_boundary_prime: tuple[bool, bool, bool]
    other_root_discriminant_identity: Expr
    other_p_derivative_resultant_identity: Expr
    t111_minor_gcd_identity: Expr
    t111_normalized_minor_gcd: Expr
    t112_minor_gcd_identity: Expr
    t112_normalized_minor_gcd: Expr
    critical_tangent_reduction_identity: Expr
    t112_immersion_determinant_identity: Expr
    exceptional_immersion_residuals: tuple[Expr, Expr, Expr, Expr]
    t112_component_dimension: int
    t112_same_profile_threefold_exists: bool
    c2_incidence_determinant_identity: Expr
    c2_rank_two_minor_gcd_identity: Expr
    c2_normalized_rank_two_minor_gcd: Expr
    pair_denominator_incidence_identity: Expr
    c2_component_dimension: int
    c2_same_profile_threefold_exists: bool
    split_vertical_determinant_identities: tuple[Expr, Expr]
    split_vertical_rank_data: tuple[tuple[tuple[int, int], ...], ...]
    t112_sample_incidence_residuals: tuple[Expr, Expr, Expr]
    t112_sample_critical_derivative: tuple[Expr, Expr]
    t112_sample_images: tuple[Expr, Expr, Expr, Expr]
    t112_sample_contact_jets: tuple[Expr, Expr]
    t112_sample_collision_identity: Expr
    t112_sample_collision_gcd: Expr
    t112_sample_residual_discriminant: Expr
    t112_sample_residual_resultants: tuple[Expr, Expr]
    t112_sample_node_x_identity: Expr
    t112_sample_node_x_discriminant: Expr
    t112_sage_jacobian_components: tuple[tuple[int, int], ...]
    c2_sample_incidence_residuals: tuple[Expr, Expr, Expr, Expr]
    c2_sample_critical_derivative: tuple[Expr, Expr]
    c2_sample_images: tuple[Expr, Expr, Expr, Expr]
    c2_sample_triple_slope_difference: Expr
    c2_sample_contact_jets: tuple[Expr, Expr]
    c2_sample_collision_identity: Expr
    c2_sample_collision_gcd: Expr
    c2_sample_residual_discriminant: Expr
    c2_sample_residual_resultants: tuple[Expr, Expr, Expr]
    c2_sample_node_x_identity: Expr
    c2_sample_node_x_discriminant: Expr
    c2_sage_jacobian_components: tuple[tuple[int, int], ...]
    topology_computed: bool
    topology_propagation_dependencies: tuple[str, str]

    @property
    def verified(self) -> bool:
        """Whether all exact component and hostile-fixture checks agree."""

        return bool(
            self.critical_p_factorization_identity == 0
            and all(
                all(remainder == 0 for remainder in row)
                for row in self.labeled_component_remainders
            )
            and self.labeled_boundary_saturation_exponent == 1
            and self.labeled_boundary_component_count == 3
            and self.labeled_boundary_dimensions == (1, 1, 1)
            and self.labeled_boundary_prime == (True, True, True)
            and self.other_root_discriminant_identity == 0
            and self.other_p_derivative_resultant_identity == 0
            and self.t111_minor_gcd_identity == 0
            and self.t111_normalized_minor_gcd == 1
            and self.t112_minor_gcd_identity == 0
            and self.t112_normalized_minor_gcd == 1
            and self.critical_tangent_reduction_identity == 0
            and self.t112_immersion_determinant_identity == 0
            and self.exceptional_immersion_residuals == (0, 0, 0, 100)
            and self.t112_component_dimension == 2
            and not self.t112_same_profile_threefold_exists
            and self.c2_incidence_determinant_identity == 0
            and self.c2_rank_two_minor_gcd_identity == 0
            and self.c2_normalized_rank_two_minor_gcd == 1
            and self.pair_denominator_incidence_identity == 0
            and self.c2_component_dimension == 2
            and not self.c2_same_profile_threefold_exists
            and self.split_vertical_determinant_identities == (0, 0)
            and self.split_vertical_rank_data
            == (
                ((3, 3), (3, 4), (3, 4)),
                ((3, 3), (3, 4), (3, 4)),
            )
            and self.t112_sample_incidence_residuals == (0, 0, 0)
            and self.t112_sample_critical_derivative == (0, 588)
            and self.t112_sample_images == (-4, 8, -4, 8)
            and self.t112_sample_contact_jets == (0, Rational(23, 343) * (2 * T + 1))
            and self.t112_sample_collision_identity == 0
            and self.t112_sample_collision_gcd == S + 1
            and self.t112_sample_residual_discriminant == 6340872382174352400
            and self.t112_sample_residual_resultants
            == (55641600, 124085451579933126415376475000000)
            and self.t112_sample_node_x_identity == 0
            and self.t112_sample_node_x_discriminant
            == 178514670159121580068048761763564718777055131533520941239105960122799512533436039430144000000
            and self.t112_sage_jacobian_components == ((4, 1), (6, 1), (6, 6))
            and self.c2_sample_incidence_residuals == (0, 0, 0, 0)
            and self.c2_sample_critical_derivative == (0, Rational(504, 5))
            and self.c2_sample_images == (-4, Rational(272, 5), 2, Rational(89, 5))
            and self.c2_sample_triple_slope_difference
            == Rational(-58, 35) * (2 * T + 1)
            and self.c2_sample_contact_jets == (0, Rational(-14, 845) * (2 * T - 1))
            and self.c2_sample_collision_identity == 0
            and self.c2_sample_collision_gcd == S - 1
            and self.c2_sample_residual_discriminant == -697306407293437500
            and self.c2_sample_residual_resultants
            == (7, 225504000, -8638634891348924179843125000000)
            and self.c2_sample_node_x_identity == 0
            and self.c2_sample_node_x_discriminant
            == -963456383583343481231569329837685661917945053301086685463838720000000000000
            and self.c2_sage_jacobian_components == ((4, 1), (3, 1), (4, 1), (5, 5))
            and not self.topology_computed
            and self.topology_propagation_dependencies
            == (
                "connected clean open",
                "proper projective Whitney-Thom triviality",
            )
        )


@cache
def exact_delta_ten_pcritical_triple_certificate() -> (
    DeltaTenPCriticalTripleCertificate
):
    """Build the exact critical-boundary and hostile-fixture certificate."""

    labeled_remainders: list[tuple[Expr, Expr]] = []
    for fourth_factor, component in zip(
        ROOT_FOURTH_FACTORS,
        ROOT_CRITICAL_COMPONENTS,
        strict=True,
    ):
        basis = groebner(component, ROOT_A, ROOT_B, ROOT_C, order="grevlex")
        labeled_remainders.append(
            (
                basis.reduce(ROOT_BASE_CONSTRAINT)[1],
                basis.reduce(fourth_factor)[1],
            )
        )

    t111_cleared_matrix = T112_CLEARED_MATRIX[:2, :]
    t111_minors = _maximal_minors(t111_cleared_matrix)
    t111_minor_gcd = _polynomial_gcd(t111_minors)
    expected_t111_gcd = (
        59049
        * CRITICAL_ROOT**13
        * (CRITICAL_ROOT - 1) ** 2
        * (CRITICAL_ROOT + 1) ** 2
        * (2 * CRITICAL_ROOT**2 - 1)
    )
    t111_normalized_minors = tuple(
        cancel(minor / t111_minor_gcd) for minor in t111_minors
    )

    t112_minors = _maximal_minors(T112_CLEARED_MATRIX)
    t112_minor_gcd = _polynomial_gcd(t112_minors)
    t112_normalized_minors = tuple(
        cancel(minor / t112_minor_gcd) for minor in t112_minors
    )

    exceptional_residuals = tuple(
        _modulo_exceptional_quadratic(equation.subs(EXCEPTIONAL_IMMERSION_SOLUTION))
        for equation in (*T112_CRITICAL_EQUATIONS, CRITICAL_Q_DERIVATIVE)
    )

    c2_rank_two_minors = _all_minors(C2_CLEARED_MATRIX, 3)
    c2_rank_two_gcd = _polynomial_gcd(c2_rank_two_minors)
    c2_normalized_rank_two_minors = tuple(
        cancel(minor / c2_rank_two_gcd) for minor in c2_rank_two_minors
    )

    split_determinants = tuple(
        matrix.det(method="domain-ge") for matrix in SPLIT_VERTICAL_MATRICES
    )
    split_rank_data: list[tuple[tuple[int, int], ...]] = []
    for equations, matrix in zip(
        SPLIT_VERTICAL_EQUATIONS,
        SPLIT_VERTICAL_MATRICES,
        strict=True,
    ):
        constant_column = Matrix(
            [-equation.subs(ZERO_LINEAR_PARAMETERS) for equation in equations]
        )
        augmented = matrix.row_join(constant_column)
        split_rank_data.append(
            tuple(
                (
                    matrix.subs(SPLIT_PAIR_PRODUCT, root).rank(),
                    augmented.subs(SPLIT_PAIR_PRODUCT, root).rank(),
                )
                for root in SPLIT_VERTICAL_RANKDROP_ROOTS
            )
        )

    t112_sample_collision = COLLISION_POLYNOMIAL.subs(PCRITICAL_T112_PARAMETERS)
    t112_sample_tangency = TANGENCY_POLYNOMIAL.subs(PCRITICAL_T112_PARAMETERS)
    t112_p_derivative = diff(PCRITICAL_T112_P, T)
    t112_first_jet = cancel(diff(PCRITICAL_T112_Q, T) / t112_p_derivative)
    t112_second_jet = cancel(diff(t112_first_jet, T) / t112_p_derivative)
    t112_jets = tuple(
        _modulo_branch(jet, PCRITICAL_T112_OTHER_BRANCH)
        for jet in (t112_first_jet, t112_second_jet)
    )
    t112_partner_jets = tuple(
        _modulo_branch(jet.subs(T, -1 - T), PCRITICAL_T112_OTHER_BRANCH)
        for jet in (t112_first_jet, t112_second_jet)
    )
    t112_jet_differences = tuple(
        expand(rem(left - right, PCRITICAL_T112_OTHER_BRANCH, T))
        for left, right in zip(t112_jets, t112_partner_jets, strict=True)
    )
    t112_x_numerator = COLLISION_X_NUMERATOR.subs(PCRITICAL_T112_PARAMETERS)
    t112_x_denominator = COLLISION_X_DENOMINATOR.subs(PCRITICAL_T112_PARAMETERS)

    c2_sample_collision = COLLISION_POLYNOMIAL.subs(PCRITICAL_C2_PARAMETERS)
    c2_sample_tangency = TANGENCY_POLYNOMIAL.subs(PCRITICAL_C2_PARAMETERS)
    c2_p_derivative = diff(PCRITICAL_C2_P, T)
    c2_slope = _modulo_branch(
        diff(PCRITICAL_C2_Q, T) / c2_p_derivative,
        PCRITICAL_C2_OTHER_BRANCH,
    )
    c2_first_jet = cancel(diff(PCRITICAL_C2_Q, T) / c2_p_derivative)
    c2_second_jet = cancel(diff(c2_first_jet, T) / c2_p_derivative)
    c2_contact_jets = tuple(
        _modulo_branch(jet, PCRITICAL_C2_CONTACT_BRANCH)
        for jet in (c2_first_jet, c2_second_jet)
    )
    c2_partner_jets = tuple(
        _modulo_branch(jet.subs(T, 1 - T), PCRITICAL_C2_CONTACT_BRANCH)
        for jet in (c2_first_jet, c2_second_jet)
    )
    c2_contact_jet_differences = tuple(
        expand(rem(left - right, PCRITICAL_C2_CONTACT_BRANCH, T))
        for left, right in zip(c2_contact_jets, c2_partner_jets, strict=True)
    )
    c2_x_numerator = COLLISION_X_NUMERATOR.subs(PCRITICAL_C2_PARAMETERS)
    c2_x_denominator = COLLISION_X_DENOMINATOR.subs(PCRITICAL_C2_PARAMETERS)

    return DeltaTenPCriticalTripleCertificate(
        critical_p_factorization_identity=CRITICAL_P_FACTORIZATION_IDENTITY,
        labeled_component_remainders=tuple(labeled_remainders),
        labeled_boundary_saturation_exponent=1,
        labeled_boundary_component_count=3,
        labeled_boundary_dimensions=(1, 1, 1),
        labeled_boundary_prime=(True, True, True),
        other_root_discriminant_identity=cancel(
            OTHER_ROOT_DISCRIMINANT
            + 4
            * (CRITICAL_ROOT - 1)
            * (CRITICAL_ROOT + 1)
            * (2 * CRITICAL_ROOT**2 + 1)
            / (9 * CRITICAL_ROOT**2)
        ),
        other_p_derivative_resultant_identity=cancel(
            OTHER_P_DERIVATIVE_RESULTANT
            - 4
            * (CRITICAL_ROOT - 1)
            * (CRITICAL_ROOT + 1)
            * (2 * CRITICAL_ROOT**2 - 1) ** 2
            * (2 * CRITICAL_ROOT**2 + 1)
            / (9 * CRITICAL_ROOT**2)
        ),
        t111_minor_gcd_identity=expand(t111_minor_gcd - expected_t111_gcd),
        t111_normalized_minor_gcd=_polynomial_gcd(t111_normalized_minors),
        t112_minor_gcd_identity=expand(t112_minor_gcd - T112_EXPECTED_MINOR_GCD),
        t112_normalized_minor_gcd=_polynomial_gcd(t112_normalized_minors),
        critical_tangent_reduction_identity=CRITICAL_TANGENT_REDUCTION_IDENTITY,
        t112_immersion_determinant_identity=cancel(
            T112_IMMERSION_MATRIX.det(method="domain-ge")
            - EXPECTED_T112_IMMERSION_DETERMINANT
        ),
        exceptional_immersion_residuals=exceptional_residuals,
        t112_component_dimension=2,
        t112_same_profile_threefold_exists=False,
        c2_incidence_determinant_identity=cancel(
            C2_INCIDENCE_MATRIX.det(method="domain-ge")
            - EXPECTED_C2_INCIDENCE_DETERMINANT
        ),
        c2_rank_two_minor_gcd_identity=expand(
            c2_rank_two_gcd + 2 * C2_PAIR_DENOMINATOR_FACTOR
        ),
        c2_normalized_rank_two_minor_gcd=_polynomial_gcd(c2_normalized_rank_two_minors),
        pair_denominator_incidence_identity=cancel(
            (CONTACT_SUM * (CONTACT_SUM**2 + CRITICAL_KAPPA * CONTACT_SUM + 1)).subs(
                CONTACT_SUM, PAIR_DENOMINATOR_CONTACT_SUM
            )
            - EXPECTED_PAIR_INCIDENCE_RESIDUAL
        ),
        c2_component_dimension=2,
        c2_same_profile_threefold_exists=False,
        split_vertical_determinant_identities=tuple(
            expand(actual - expected)
            for actual, expected in zip(
                split_determinants,
                EXPECTED_SPLIT_VERTICAL_DETERMINANTS,
                strict=True,
            )
        ),
        split_vertical_rank_data=tuple(split_rank_data),
        t112_sample_incidence_residuals=tuple(
            cancel(equation.subs(CRITICAL_ROOT, 2).subs(PCRITICAL_T112_PARAMETERS))
            for equation in T112_CRITICAL_EQUATIONS
        ),
        t112_sample_critical_derivative=(
            t112_p_derivative.subs(T, 2),
            diff(PCRITICAL_T112_Q, T).subs(T, 2),
        ),
        t112_sample_images=(
            PCRITICAL_T112_P.subs(T, 2),
            PCRITICAL_T112_Q.subs(T, 2),
            rem(PCRITICAL_T112_P, PCRITICAL_T112_OTHER_BRANCH, T),
            rem(PCRITICAL_T112_Q, PCRITICAL_T112_OTHER_BRANCH, T),
        ),
        t112_sample_contact_jets=t112_jet_differences,
        t112_sample_collision_identity=expand(
            t112_sample_collision - PCRITICAL_T112_COLLISION
        ),
        t112_sample_collision_gcd=gcd(
            Poly(t112_sample_collision, S),
            Poly(t112_sample_tangency, S),
        ).as_expr(),
        t112_sample_residual_discriminant=discriminant(
            PCRITICAL_T112_RESIDUAL_NODES,
            S,
        ),
        t112_sample_residual_resultants=(
            resultant(
                PCRITICAL_T112_RESIDUAL_NODES,
                PCRITICAL_T112_TRIPLE_PAIR_SUMS,
                S,
            ),
            resultant(
                PCRITICAL_T112_RESIDUAL_NODES,
                t112_sample_tangency,
                S,
            ),
        ),
        t112_sample_node_x_identity=expand(
            resultant(
                PCRITICAL_T112_RESIDUAL_NODES,
                t112_x_numerator - X * t112_x_denominator,
                S,
            )
            - 5625 * PCRITICAL_T112_NODE_X
        ),
        t112_sample_node_x_discriminant=discriminant(
            PCRITICAL_T112_NODE_X,
            X,
        ),
        t112_sage_jacobian_components=((4, 1), (6, 1), (6, 6)),
        c2_sample_incidence_residuals=tuple(
            cancel(
                equation.subs({CRITICAL_ROOT: 2, CONTACT_SUM: 1}).subs(
                    PCRITICAL_C2_PARAMETERS
                )
            )
            for equation in C2_CRITICAL_EQUATIONS
        ),
        c2_sample_critical_derivative=(
            c2_p_derivative.subs(T, 2),
            diff(PCRITICAL_C2_Q, T).subs(T, 2),
        ),
        c2_sample_images=(
            PCRITICAL_C2_P.subs(T, 2),
            PCRITICAL_C2_Q.subs(T, 2),
            rem(PCRITICAL_C2_P, PCRITICAL_C2_CONTACT_BRANCH, T),
            rem(PCRITICAL_C2_Q, PCRITICAL_C2_CONTACT_BRANCH, T),
        ),
        c2_sample_triple_slope_difference=expand(
            rem(
                c2_slope - c2_slope.subs(T, -1 - T),
                PCRITICAL_C2_OTHER_BRANCH,
                T,
            )
        ),
        c2_sample_contact_jets=c2_contact_jet_differences,
        c2_sample_collision_identity=expand(
            c2_sample_collision - PCRITICAL_C2_COLLISION
        ),
        c2_sample_collision_gcd=gcd(
            Poly(c2_sample_collision, S),
            Poly(c2_sample_tangency, S),
        ).as_expr(),
        c2_sample_residual_discriminant=discriminant(
            PCRITICAL_C2_RESIDUAL_NODES,
            S,
        ),
        c2_sample_residual_resultants=(
            resultant(PCRITICAL_C2_RESIDUAL_NODES, S - 1, S),
            resultant(
                PCRITICAL_C2_RESIDUAL_NODES,
                PCRITICAL_C2_TRIPLE_PAIR_SUMS,
                S,
            ),
            resultant(
                PCRITICAL_C2_RESIDUAL_NODES,
                c2_sample_tangency,
                S,
            ),
        ),
        c2_sample_node_x_identity=expand(
            resultant(
                PCRITICAL_C2_RESIDUAL_NODES,
                c2_x_numerator - X * c2_x_denominator,
                S,
            )
            + 140625 * PCRITICAL_C2_NODE_X
        ),
        c2_sample_node_x_discriminant=discriminant(
            PCRITICAL_C2_NODE_X,
            X,
        ),
        c2_sage_jacobian_components=((4, 1), (3, 1), (4, 1), (5, 5)),
        topology_computed=False,
        topology_propagation_dependencies=(
            "connected clean open",
            "proper projective Whitney-Thom triviality",
        ),
    )


def main() -> int:
    """Print the exact critical-boundary classification summary."""

    certificate = exact_delta_ten_pcritical_triple_certificate()
    print("P-critical triple boundary certificate:", certificate.verified)
    print(
        "labeled prime boundary curves:", certificate.labeled_boundary_component_count
    )
    print(
        "same-profile component dimensions:",
        {
            "T112": certificate.t112_component_dimension,
            "C2+T111": certificate.c2_component_dimension,
        },
    )
    print("topology computed:", certificate.topology_computed)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
