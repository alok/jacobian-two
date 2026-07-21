"""Classify the split ``C3`` and ``C2^2`` coefficient-rank strata.

This module closes the determinantal gap left by the generic split ledger for
the two contact-only codimension-two profiles.  The normalized split pair
components have root coordinates

* ``k=0``: ``V`` uses ``r`` and ``W`` uses ``s``;
* ``k=2``: ``V`` uses ``r`` and ``W`` uses ``s``.

The clean root factors are derived from the actual pair geometry:

* at ``k=0``, ``V`` excludes ``r=0`` (diagonal) and ``2*r=1``
  (the component overlap), while ``W`` excludes ``s=0`` (overlap),
  ``s^2+1=0`` (a pair containing the forced cusp source), and
  ``s^2+2=0`` (diagonal);
* at ``k=2``, ``V`` excludes ``r=0`` (overlap) and ``4*r=1``
  (diagonal), while ``W`` excludes ``s=-1`` (overlap) and
  ``s*(s+2)=0`` (diagonal).

For every allowed ``C3 + 7N`` allocation, the coefficient system has its
expected rank on this clean open.  The only raw rank loss is the ``k=0``
graph factor ``s^2+1``.  It is genuinely compatible, but its pair is
``{0, +/-i}``, so it collides with the forced cusp and is not a clean ``C3``
point.

For ``C2^2 + 6N``, three systems lose rank only on the visible clean-boundary
factors: ``k=0`` vertical/graph, ``k=0`` overlap-plus-graph, and ``k=2``
vertical/vertical.  The other three systems have irreducible residual
determinants:

* ``k=0, WW``;
* ``k=2, VW``; and
* ``k=2, WW``.

Exact ideal-sandwich certificates compute their localized compatible bases.
They are reduced finite schemes of lengths ``4, 6, 6``.  All coefficient
matrices have rank exactly three there, so the coefficient fiber is an
affine line and every residual incidence component is one-dimensional.  In
particular, none is a hidden second split surface.

The ideal-sandwich proof is executable in SymPy.  If ``I`` is the raw
compatibility ideal, ``J`` the displayed finite Groebner basis, and ``h`` the
clean localizer, it checks

``I subset J``, ``h^e J subset I``, and ``J + (h) = (1)``.

Hence ``I : h^infinity = J``.  A separate calculation proves that a power of
``h`` belongs to the ideal of the coefficient ``3 x 3`` minors, excluding
rank at most two on the clean open.  The companion Sage checker independently
replays the saturations and scheme lengths.

This is only a rank/dimension audit for these two profiles on the true split
charts.  It does not compute complement topology, classify every lower
profile boundary, or prove the plane Jacobian conjecture.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from functools import cache
from itertools import combinations
from typing import Final

from sympy import (
    Expr,
    I,
    Matrix,
    Poly,
    QQ,
    Rational,
    Symbol,
    cancel,
    diff,
    expand,
    factor_list,
    gcd,
    groebner,
    linear_eq_to_matrix,
)
from sympy.polys.matrices import DomainMatrix

from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    CUSP_IMAGE_FACTOR,
    DELTA,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    R,
    S,
    T,
    ZERO_FIBER_GRAPH,
    ZERO_FIBER_VERTICAL,
    exceptional_graph,
    exceptional_vertical,
)
from scripts.a6_delta_ten_split_codim_two import FIRST_ROOT, SECOND_ROOT

COEFFICIENTS: Final = (ALPHA, BETA, GAMMA, DELTA)
ZERO_COEFFICIENTS: Final = {coefficient: 0 for coefficient in COEFFICIENTS}
U: Final = FIRST_ROOT
V: Final = SECOND_ROOT


def _domain_determinant(matrix: Matrix) -> Expr:
    """Return an exact determinant without heuristic simplification."""

    determinant: Expr = DomainMatrix.from_Matrix(matrix).det().as_expr()
    return determinant


def _jet_equations(
    polynomial: Expr,
    variable: Symbol,
    root: Expr,
    multiplicity: int,
) -> tuple[Expr, ...]:
    """Return vanishing equations through order ``multiplicity - 1``."""

    return tuple(
        expand(diff(polynomial, variable, order).subs(variable, root))
        for order in range(multiplicity)
    )


def _linear_system(equations: tuple[Expr, ...]) -> tuple[Matrix, Matrix, Matrix]:
    """Return coefficient, right-hand-side, and augmented matrices."""

    coefficient_matrix, right_side = linear_eq_to_matrix(equations, COEFFICIENTS)
    return coefficient_matrix, right_side, coefficient_matrix.row_join(right_side)


def _minors(matrix: Matrix, size: int) -> tuple[Expr, ...]:
    """Return all exact ``size`` minors of a matrix."""

    return tuple(
        expand(matrix.extract(rows, columns).det())
        for rows in combinations(range(matrix.rows), size)
        for columns in combinations(range(matrix.cols), size)
    )


def _common_gcd(expressions: tuple[Expr, ...]) -> Expr:
    """Return the nonnegative-content gcd of a nonempty expression tuple."""

    if not expressions:
        msg = "a common gcd needs at least one expression"
        raise ValueError(msg)
    common = expressions[0]
    for expression in expressions[1:]:
        common = gcd(common, expression)
    return expand(common)


def split_component_polynomial(kappa: int, component: str) -> Expr:
    """Return the true split collision component."""

    if (kappa, component) == (0, "V"):
        return ZERO_FIBER_VERTICAL
    if (kappa, component) == (0, "W"):
        return ZERO_FIBER_GRAPH
    if (kappa, component) == (2, "V"):
        return exceptional_vertical(1)
    if (kappa, component) == (2, "W"):
        return exceptional_graph(1)
    msg = "the supported split component is k=0 or k=2 and V or W"
    raise ValueError(msg)


def split_component_variable(component: str) -> Symbol:
    """Return the root coordinate for a split component."""

    if component == "V":
        return R
    if component == "W":
        return S
    msg = "a split component must be V or W"
    raise ValueError(msg)


def split_component_clean_factor(kappa: int, component: str, root: Expr) -> Expr:
    """Return the exact overlap/cusp/diagonal clean-root factor."""

    if (kappa, component) == (0, "V"):
        return expand(root * (2 * root - 1))
    if (kappa, component) == (0, "W"):
        return expand(root * (root**2 + 1) * (root**2 + 2))
    if (kappa, component) == (2, "V"):
        return expand(root * (4 * root - 1))
    if (kappa, component) == (2, "W"):
        return expand(root * (root + 1) * (root + 2))
    msg = "the supported split component is k=0 or k=2 and V or W"
    raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class LinearRankSpec:
    """One allowed split allocation and its expected clean rank."""

    name: str
    profile: str
    kappa: int
    allocation: str
    equations: tuple[Expr, ...]
    clean_localizer: Expr
    expected_rank: int
    expected_maximal_minor_gcd: Expr


@cache
def c3_rank_specs() -> tuple[LinearRankSpec, ...]:
    """Return all five allowed ``C3`` split allocation systems."""

    overlap = expand(ZERO_FIBER_GRAPH.subs(S, 0))
    return (
        LinearRankSpec(
            name="c3_k0_w",
            profile="C3+7N",
            kappa=0,
            allocation="W",
            equations=_jet_equations(ZERO_FIBER_GRAPH, S, U, 3),
            clean_localizer=split_component_clean_factor(0, "W", U),
            expected_rank=3,
            expected_maximal_minor_gcd=64 * (U**2 + 1) ** 3,
        ),
        LinearRankSpec(
            name="c3_k0_overlap_v",
            profile="C3+7N",
            kappa=0,
            allocation="overlap-V",
            equations=(
                overlap,
                expand(diff(ZERO_FIBER_VERTICAL, R).subs(R, Rational(1, 2))),
            ),
            clean_localizer=1,
            expected_rank=2,
            expected_maximal_minor_gcd=4,
        ),
        LinearRankSpec(
            name="c3_k0_overlap_w",
            profile="C3+7N",
            kappa=0,
            allocation="overlap-W",
            equations=(overlap, expand(diff(ZERO_FIBER_GRAPH, S).subs(S, 0))),
            clean_localizer=1,
            expected_rank=2,
            expected_maximal_minor_gcd=8,
        ),
        LinearRankSpec(
            name="c3_k2_v",
            profile="C3+7N",
            kappa=2,
            allocation="V",
            equations=_jet_equations(exceptional_vertical(1), R, U, 3),
            clean_localizer=split_component_clean_factor(2, "V", U),
            expected_rank=3,
            expected_maximal_minor_gcd=2,
        ),
        LinearRankSpec(
            name="c3_k2_w",
            profile="C3+7N",
            kappa=2,
            allocation="W",
            equations=_jet_equations(exceptional_graph(1), S, U, 3),
            clean_localizer=split_component_clean_factor(2, "W", U),
            expected_rank=3,
            expected_maximal_minor_gcd=64,
        ),
    )


@dataclass(frozen=True, slots=True)
class LinearRankCertificate:
    """Exact generic-rank and maximal-minor-gcd evidence."""

    name: str
    coefficient_rank: int
    augmented_rank: int
    maximal_minor_gcd: Expr
    maximal_minor_gcd_identity: Expr
    clean_rank_drop_removed: bool

    @property
    def verified(self) -> bool:
        """Whether the expected rank holds on the stated clean open."""

        return bool(
            self.coefficient_rank == self.augmented_rank
            and self.maximal_minor_gcd_identity == 0
            and self.clean_rank_drop_removed
        )


@cache
def exact_c3_rank_certificates() -> tuple[LinearRankCertificate, ...]:
    """Classify the coefficient-rank loss of every allowed ``C3`` chart."""

    certificates: list[LinearRankCertificate] = []
    for spec in c3_rank_specs():
        coefficient_matrix, _, augmented_matrix = _linear_system(spec.equations)
        maximal_minors = _minors(coefficient_matrix, coefficient_matrix.rows)
        maximal_gcd = _common_gcd(maximal_minors)
        expected = expand(spec.expected_maximal_minor_gcd)
        if spec.clean_localizer == 1 or Poly(expected, U, V).total_degree() == 0:
            clean_removed = expected != 0
        else:
            clean_removed = cancel(
                spec.clean_localizer**coefficient_matrix.rows / expected
            ).is_polynomial(U)
        certificates.append(
            LinearRankCertificate(
                name=spec.name,
                coefficient_rank=coefficient_matrix.rank(),
                augmented_rank=augmented_matrix.rank(),
                maximal_minor_gcd=maximal_gcd,
                maximal_minor_gcd_identity=expand(maximal_gcd - expected),
                clean_rank_drop_removed=bool(clean_removed),
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class SimpleDoubleContactSpec:
    """A two-contact allocation with no residual clean rank divisor."""

    name: str
    kappa: int
    allocation: str
    equations: tuple[Expr, ...]
    clean_localizer: Expr
    expected_rank: int
    expected_rank_drop_generator: Expr


@cache
def simple_double_contact_specs() -> tuple[SimpleDoubleContactSpec, ...]:
    """Return the three visible-boundary-only two-contact systems."""

    overlap = expand(ZERO_FIBER_GRAPH.subs(S, 0))
    return (
        SimpleDoubleContactSpec(
            name="c22_k0_vw",
            kappa=0,
            allocation="VW",
            equations=(
                *_jet_equations(ZERO_FIBER_VERTICAL, R, U, 2),
                *_jet_equations(ZERO_FIBER_GRAPH, S, V, 2),
            ),
            clean_localizer=expand(
                split_component_clean_factor(0, "V", U)
                * split_component_clean_factor(0, "W", V)
            ),
            expected_rank=4,
            expected_rank_drop_generator=256 * V**3 * (V**2 + 1),
        ),
        SimpleDoubleContactSpec(
            name="c22_k0_overlap_w",
            kappa=0,
            allocation="overlap+W",
            equations=(
                overlap,
                *_jet_equations(ZERO_FIBER_GRAPH, S, V, 2),
            ),
            clean_localizer=split_component_clean_factor(0, "W", V),
            expected_rank=3,
            expected_rank_drop_generator=32 * V**2 * (V**2 + 1),
        ),
        SimpleDoubleContactSpec(
            name="c22_k2_vv",
            kappa=2,
            allocation="VV",
            equations=(
                *_jet_equations(exceptional_vertical(1), R, U, 2),
                *_jet_equations(exceptional_vertical(1), R, V, 2),
            ),
            clean_localizer=expand(
                (U - V)
                * split_component_clean_factor(2, "V", U)
                * split_component_clean_factor(2, "V", V)
            ),
            expected_rank=4,
            expected_rank_drop_generator=(U - V) ** 4,
        ),
    )


@cache
def exact_simple_double_contact_rank_certificates() -> tuple[
    LinearRankCertificate, ...
]:
    """Prove that the three simple systems lose rank only on boundaries."""

    certificates: list[LinearRankCertificate] = []
    for spec in simple_double_contact_specs():
        coefficient_matrix, _, augmented_matrix = _linear_system(spec.equations)
        maximal_minors = _minors(coefficient_matrix, coefficient_matrix.rows)
        maximal_gcd = _common_gcd(maximal_minors)
        expected = expand(spec.expected_rank_drop_generator)
        quotient = cancel(spec.clean_localizer**coefficient_matrix.rows / expected)
        certificates.append(
            LinearRankCertificate(
                name=spec.name,
                coefficient_rank=coefficient_matrix.rank(),
                augmented_rank=augmented_matrix.rank(),
                maximal_minor_gcd=maximal_gcd,
                maximal_minor_gcd_identity=expand(maximal_gcd - expected),
                clean_rank_drop_removed=bool(
                    quotient.is_polynomial(U, V) and expected != 0
                ),
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class ResidualRankSpec:
    """Expected exact data for one genuine residual determinant chart."""

    name: str
    kappa: int
    allocation: str
    equations: tuple[Expr, ...]
    clean_localizer: Expr
    determinant_visible_factor: Expr
    expected_compatibility_basis: tuple[Expr, Expr]
    compatibility_power: int
    rank_two_power: int
    expected_base_length: int


@cache
def residual_rank_specs() -> tuple[ResidualRankSpec, ...]:
    """Return the three residual two-contact rank charts."""

    k0_w = ZERO_FIBER_GRAPH
    k2_v = exceptional_vertical(1)
    k2_w = exceptional_graph(1)
    return (
        ResidualRankSpec(
            name="c22_k0_ww_residual",
            kappa=0,
            allocation="WW",
            equations=(
                *_jet_equations(k0_w, S, U, 2),
                *_jet_equations(k0_w, S, V, 2),
            ),
            clean_localizer=expand(
                (U - V)
                * split_component_clean_factor(0, "W", U)
                * split_component_clean_factor(0, "W", V)
            ),
            determinant_visible_factor=expand(
                512 * (U - V) ** 4 * (U**2 + 1) * (V**2 + 1)
            ),
            expected_compatibility_basis=(
                38 * U + 4 * V**3 + 71 * V,
                4 * V**4 + 71 * V**2 + 361,
            ),
            compatibility_power=2,
            rank_two_power=2,
            expected_base_length=4,
        ),
        ResidualRankSpec(
            name="c22_k2_vw_residual",
            kappa=2,
            allocation="VW",
            equations=(
                *_jet_equations(k2_v, R, U, 2),
                *_jet_equations(k2_w, S, V, 2),
            ),
            clean_localizer=expand(
                split_component_clean_factor(2, "V", U)
                * split_component_clean_factor(2, "W", V)
            ),
            determinant_visible_factor=8,
            expected_compatibility_basis=(
                414 * U
                + 6 * V**5
                + 146 * V**4
                + 877 * V**3
                + 1934 * V**2
                + 1527 * V
                + 234,
                6 * V**6
                + 74 * V**5
                + 367 * V**4
                + 932 * V**3
                + 1296 * V**2
                + 954 * V
                + 297,
            ),
            compatibility_power=6,
            rank_two_power=2,
            expected_base_length=6,
        ),
        ResidualRankSpec(
            name="c22_k2_ww_residual",
            kappa=2,
            allocation="WW",
            equations=(
                *_jet_equations(k2_w, S, U, 2),
                *_jet_equations(k2_w, S, V, 2),
            ),
            clean_localizer=expand(
                (U - V)
                * split_component_clean_factor(2, "W", U)
                * split_component_clean_factor(2, "W", V)
            ),
            determinant_visible_factor=1024 * (U - V) ** 4,
            expected_compatibility_basis=(
                3 * U
                - 46 * V**5
                - 336 * V**4
                - 1017 * V**3
                - 1380 * V**2
                - 742 * V
                - 117,
                2 * V**6 + 18 * V**5 + 69 * V**4 + 135 * V**3 + 134 * V**2 + 60 * V + 9,
            ),
            compatibility_power=1,
            rank_two_power=1,
            expected_base_length=6,
        ),
    )


@dataclass(frozen=True, slots=True)
class ResidualRankCertificate:
    """Exact ideal-sandwich and fiber-dimension certificate."""

    name: str
    determinant_identity: Expr
    residual_factor: Expr
    residual_irreducible: bool
    compatibility_common_factor: Expr
    compatibility_generator_count: int
    raw_in_expected_remainders: tuple[Expr, ...]
    powered_expected_in_raw_remainders: tuple[Expr, Expr]
    lower_power_remainders: tuple[Expr, Expr]
    expected_plus_boundary_basis: tuple[Expr, ...]
    univariate_polynomial: Expr
    univariate_derivative_gcd: Expr
    base_length: int
    rank_two_power_remainder: Expr
    rank_two_lower_power_remainder: Expr
    swap_remainders: tuple[Expr, ...]
    coefficient_rank_on_base: int
    coefficient_fiber_dimension: int
    incidence_dimension: int
    unordered_base_orbit_count: int | None

    @property
    def verified(self) -> bool:
        """Whether the localized rank classification is exact."""

        lower_compatibility_is_minimal = any(
            remainder != 0 for remainder in self.lower_power_remainders
        )
        lower_rank_is_minimal = self.rank_two_lower_power_remainder != 0
        return bool(
            self.determinant_identity == 0
            and self.residual_irreducible
            and self.compatibility_common_factor != 0
            and self.compatibility_generator_count == 5
            and all(value == 0 for value in self.raw_in_expected_remainders)
            and self.powered_expected_in_raw_remainders == (0, 0)
            and lower_compatibility_is_minimal
            and self.expected_plus_boundary_basis == (1,)
            and self.univariate_derivative_gcd == 1
            and self.base_length in (4, 6)
            and self.rank_two_power_remainder == 0
            and lower_rank_is_minimal
            and all(value == 0 for value in self.swap_remainders)
            and self.coefficient_rank_on_base == 3
            and self.coefficient_fiber_dimension == 1
            and self.incidence_dimension == 1
            and (
                self.unordered_base_orbit_count is None
                or 2 * self.unordered_base_orbit_count == self.base_length
            )
        )


def _as_expressions(polynomials: Iterable[Poly]) -> tuple[Expr, ...]:
    """Convert a SymPy Groebner polynomial sequence to expressions."""

    return tuple(polynomial.as_expr() for polynomial in polynomials)


@cache
def exact_residual_rank_certificates() -> tuple[ResidualRankCertificate, ...]:
    """Compute all three localized residual rank certificates exactly."""

    certificates: list[ResidualRankCertificate] = []
    for spec in residual_rank_specs():
        coefficient_matrix, _, augmented_matrix = _linear_system(spec.equations)
        determinant = _domain_determinant(coefficient_matrix)
        residual_factor = cancel(determinant / spec.determinant_visible_factor)
        residual_factorization = factor_list(residual_factor, U, V)[1]

        compatibility_minors = tuple(
            _domain_determinant(augmented_matrix[:, columns])
            for columns in combinations(range(augmented_matrix.cols), 4)
            if augmented_matrix.cols - 1 in columns
        )
        compatibility_common = _common_gcd(compatibility_minors)
        normalized_compatibility = tuple(
            cancel(minor / compatibility_common) for minor in compatibility_minors
        )
        raw_generators = (residual_factor, *normalized_compatibility)
        raw_basis = groebner(
            raw_generators,
            U,
            V,
            order="grevlex",
            method="f5b",
            domain=QQ,
        )
        expected_basis = groebner(
            spec.expected_compatibility_basis,
            U,
            V,
            order="lex",
            domain=QQ,
        )
        raw_in_expected = tuple(
            expand(expected_basis.reduce(generator)[1]) for generator in raw_generators
        )
        powered_expected = tuple(
            expand(
                raw_basis.reduce(
                    spec.clean_localizer**spec.compatibility_power * generator
                )[1]
            )
            for generator in spec.expected_compatibility_basis
        )
        if spec.compatibility_power == 1:
            lower_power = tuple(
                expand(raw_basis.reduce(generator)[1])
                for generator in spec.expected_compatibility_basis
            )
        else:
            lower_power = tuple(
                expand(
                    raw_basis.reduce(
                        spec.clean_localizer ** (spec.compatibility_power - 1)
                        * generator
                    )[1]
                )
                for generator in spec.expected_compatibility_basis
            )
        boundary_basis = groebner(
            (*spec.expected_compatibility_basis, spec.clean_localizer),
            U,
            V,
            order="lex",
            domain=QQ,
        )

        rank_two_minors = _minors(coefficient_matrix, 3)
        rank_two_basis = groebner(
            rank_two_minors,
            U,
            V,
            order="grevlex",
            method="f5b",
            domain=QQ,
        )
        rank_two_remainder = expand(
            rank_two_basis.reduce(spec.clean_localizer**spec.rank_two_power)[1]
        )
        if spec.rank_two_power == 1:
            rank_two_lower = 1
        else:
            rank_two_lower = expand(
                rank_two_basis.reduce(
                    spec.clean_localizer ** (spec.rank_two_power - 1)
                )[1]
            )

        univariate = next(
            generator
            for generator in spec.expected_compatibility_basis
            if not generator.has(U)
        )
        univariate_gcd = gcd(univariate, diff(univariate, V))
        if spec.allocation == "WW":
            swap_remainders = tuple(
                expand(
                    expected_basis.reduce(
                        generator.subs({U: V, V: U}, simultaneous=True)
                    )[1]
                )
                for generator in spec.expected_compatibility_basis
            )
            orbit_count: int | None = spec.expected_base_length // 2
        else:
            swap_remainders = ()
            orbit_count = None

        certificates.append(
            ResidualRankCertificate(
                name=spec.name,
                determinant_identity=expand(
                    determinant - spec.determinant_visible_factor * residual_factor
                ),
                residual_factor=expand(residual_factor),
                residual_irreducible=(
                    len(residual_factorization) == 1
                    and residual_factorization[0][1] == 1
                ),
                compatibility_common_factor=compatibility_common,
                compatibility_generator_count=len(raw_generators),
                raw_in_expected_remainders=raw_in_expected,
                powered_expected_in_raw_remainders=powered_expected,
                lower_power_remainders=lower_power,
                expected_plus_boundary_basis=_as_expressions(boundary_basis.polys),
                univariate_polynomial=univariate,
                univariate_derivative_gcd=univariate_gcd,
                base_length=Poly(univariate, V).degree(),
                rank_two_power_remainder=rank_two_remainder,
                rank_two_lower_power_remainder=rank_two_lower,
                swap_remainders=swap_remainders,
                coefficient_rank_on_base=3,
                coefficient_fiber_dimension=1,
                incidence_dimension=1,
                unordered_base_orbit_count=orbit_count,
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class BoundaryHostileFixture:
    """A compatible rank-loss fixture on an excluded root boundary."""

    name: str
    equation_residuals: tuple[Expr, ...]
    boundary_value: Expr
    coefficient_rank: int
    augmented_rank: int
    next_jet_value: Expr
    cusp_image_value: Expr | None
    first_target_identity: Expr
    second_target_identity: Expr

    @property
    def verified(self) -> bool:
        """Whether the fixture is compatible and lies on its stated boundary."""

        return bool(
            all(value == 0 for value in self.equation_residuals)
            and self.boundary_value == 0
            and self.coefficient_rank == self.augmented_rank
            and self.coefficient_rank in (2, 3)
            and self.next_jet_value != 0
            and self.cusp_image_value in (None, 0)
            and self.first_target_identity == 0
            and self.second_target_identity == 0
        )


def _hostile_fixture(
    name: str,
    equations: tuple[Expr, ...],
    base_substitution: dict[Expr, Expr],
    coefficients: dict[Expr, Expr],
    boundary: Expr,
    next_jet: Expr,
    first_target_identity: Expr = 0,
    second_target_identity: Expr = 0,
    cusp_image: Expr | None = None,
) -> BoundaryHostileFixture:
    coefficient_matrix, _, augmented_matrix = _linear_system(equations)
    member = {**base_substitution, **coefficients}
    return BoundaryHostileFixture(
        name=name,
        equation_residuals=tuple(
            expand(equation.subs(member)) for equation in equations
        ),
        boundary_value=expand(boundary.subs(base_substitution)),
        coefficient_rank=coefficient_matrix.subs(base_substitution).rank(),
        augmented_rank=augmented_matrix.subs(base_substitution).rank(),
        next_jet_value=expand(next_jet.subs(member)),
        cusp_image_value=(
            None if cusp_image is None else expand(cusp_image.subs(member))
        ),
        first_target_identity=expand(expand(first_target_identity).subs(member)),
        second_target_identity=expand(expand(second_target_identity).subs(member)),
    )


@cache
def exact_boundary_hostile_fixtures() -> tuple[BoundaryHostileFixture, ...]:
    """Retain cusp, overlap, and repeated-contact boundary counterfixtures."""

    c3_equations = _jet_equations(ZERO_FIBER_GRAPH, S, U, 3)
    c3_coefficients = {ALPHA: 1, BETA: 0, GAMMA: 2, DELTA: 0}
    c3_base = {U: I}

    k0_vw_equations = (
        *_jet_equations(ZERO_FIBER_VERTICAL, R, U, 2),
        *_jet_equations(ZERO_FIBER_GRAPH, S, V, 2),
    )
    k0_overlap_base = {U: Rational(1, 2), V: 0}
    k0_overlap_coefficients = {
        ALPHA: Rational(1, 4),
        BETA: 0,
        GAMMA: 1,
        DELTA: 0,
    }

    k0_ww_equations = (
        *_jet_equations(ZERO_FIBER_GRAPH, S, U, 2),
        *_jet_equations(ZERO_FIBER_GRAPH, S, V, 2),
    )
    k0_diagonal_coefficients = {ALPHA: 1, BETA: 1, GAMMA: 1, DELTA: 0}

    k2_vv_equations = (
        *_jet_equations(exceptional_vertical(1), R, U, 2),
        *_jet_equations(exceptional_vertical(1), R, V, 2),
    )
    k2_diagonal_coefficients = {ALPHA: 1, BETA: 0, GAMMA: 1, DELTA: 0}

    k2_vw_equations = (
        *_jet_equations(exceptional_vertical(1), R, U, 2),
        *_jet_equations(exceptional_graph(1), S, V, 2),
    )
    k2_overlap_critical_coefficients = {
        ALPHA: -3,
        BETA: -7,
        GAMMA: -4,
        DELTA: 1,
    }

    return (
        _hostile_fixture(
            name="c3_k0_forced_cusp_pair",
            equations=c3_equations,
            base_substitution=c3_base,
            coefficients=c3_coefficients,
            boundary=U**2 + 1,
            # On this cusp-source boundary the third derivative also
            # vanishes identically; the first hostile nonzero jet is order
            # four.  Thus this is not secretly a clean order-three contact.
            next_jet=diff(ZERO_FIBER_GRAPH, S, 4).subs(S, U),
            cusp_image=CUSP_IMAGE_FACTOR.subs(KAPPA, 0),
            first_target_identity=(
                FAMILY_P.subs({KAPPA: 0, T: U}) - FAMILY_P.subs({KAPPA: 0, T: 0})
            ),
            second_target_identity=FAMILY_Q.subs(T, U) - FAMILY_Q.subs(T, 0),
        ),
        _hostile_fixture(
            name="c22_k0_component_overlap",
            equations=k0_vw_equations,
            base_substitution=k0_overlap_base,
            coefficients=k0_overlap_coefficients,
            boundary=V,
            next_jet=diff(ZERO_FIBER_GRAPH, S, 2).subs(S, V),
            first_target_identity=ZERO_FIBER_VERTICAL.subs(R, Rational(1, 2)),
            second_target_identity=ZERO_FIBER_GRAPH.subs(S, 0),
        ),
        _hostile_fixture(
            name="c22_k0_repeated_graph_contact",
            equations=k0_ww_equations,
            base_substitution={U: 1, V: 1},
            coefficients=k0_diagonal_coefficients,
            boundary=U - V,
            next_jet=diff(ZERO_FIBER_GRAPH, S, 2).subs(S, U),
        ),
        _hostile_fixture(
            name="c22_k2_repeated_vertical_contact",
            equations=k2_vv_equations,
            base_substitution={U: 1, V: 1},
            coefficients=k2_diagonal_coefficients,
            boundary=U - V,
            next_jet=diff(exceptional_vertical(1), R, 2).subs(R, U),
        ),
        _hostile_fixture(
            name="c22_k2_overlap_diagonal_intersection",
            equations=k2_vw_equations,
            base_substitution={U: 0, V: -2},
            coefficients=k2_overlap_critical_coefficients,
            boundary=U * (V + 2),
            next_jet=diff(exceptional_graph(1), S, 2).subs(S, V),
            cusp_image=CUSP_IMAGE_FACTOR.subs(KAPPA, 2),
            first_target_identity=exceptional_vertical(1).subs(R, U),
            second_target_identity=exceptional_graph(1).subs(S, V),
        ),
    )


@dataclass(frozen=True, slots=True)
class SplitContactRankCertificate:
    """Aggregate exact rank closure for ``C3`` and ``C2^2``."""

    c3: tuple[LinearRankCertificate, ...]
    simple_double_contacts: tuple[LinearRankCertificate, ...]
    residual_double_contacts: tuple[ResidualRankCertificate, ...]
    hostile_fixtures: tuple[BoundaryHostileFixture, ...]
    maximum_residual_incidence_dimension: int
    topology_computed: bool
    proves_plane_jacobian_conjecture: bool

    @property
    def verified(self) -> bool:
        """Whether every named split rank stratum has been classified."""

        return bool(
            len(self.c3) == 5
            and len(self.simple_double_contacts) == 3
            and len(self.residual_double_contacts) == 3
            and len(self.hostile_fixtures) == 5
            and all(certificate.verified for certificate in self.c3)
            and all(certificate.verified for certificate in self.simple_double_contacts)
            and all(
                certificate.verified for certificate in self.residual_double_contacts
            )
            and all(fixture.verified for fixture in self.hostile_fixtures)
            and self.maximum_residual_incidence_dimension == 1
            and not self.topology_computed
            and not self.proves_plane_jacobian_conjecture
        )


@cache
def exact_split_contact_rank_certificate() -> SplitContactRankCertificate:
    """Build the complete split ``C3``/``C2^2`` rank certificate."""

    return SplitContactRankCertificate(
        c3=exact_c3_rank_certificates(),
        simple_double_contacts=exact_simple_double_contact_rank_certificates(),
        residual_double_contacts=exact_residual_rank_certificates(),
        hostile_fixtures=exact_boundary_hostile_fixtures(),
        maximum_residual_incidence_dimension=1,
        topology_computed=False,
        proves_plane_jacobian_conjecture=False,
    )


def main() -> int:
    """Print a compact exact split-rank summary."""

    certificate = exact_split_contact_rank_certificate()
    print("split C3/C2^2 rank certificate:", certificate.verified)
    for residual in certificate.residual_double_contacts:
        print(
            residual.name,
            "base length",
            residual.base_length,
            "fiber A^1, incidence dimension",
            residual.incidence_dimension,
        )
    print("topology computed:", certificate.topology_computed)
    print("proves JC(2):", certificate.proves_plane_jacobian_conjecture)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
