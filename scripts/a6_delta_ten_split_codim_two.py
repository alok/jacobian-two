"""Classify the split pair charts at conditional delta ten.

For the normalized family

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``,

the unordered-pair incidence is reducible at ``k = 0, +/-2``.  The genuine
collision polynomials have degree splits ``2+8`` and ``4+6`` rather than the
misleading cancellation powers in the specialized generic decic.

This module gives an exact finite ledger for the six expected
codimension-two profiles on those true split charts.  It records:

* the exact ``V/W`` component identities, intersections, diagonals, and
  cusp/critical boundaries;
* every component-length allocation compatible with the split degree budget;
* exact coefficient-matrix ranks and clean algebraic witnesses for every
  allowed generic allocation at ``k=0`` and ``k=2``;
* the involutive transport from ``k=2`` to ``k=-2``; and
* the special legitimate ``k=0`` overlap contact, where the same unordered
  pair lies on both components and has local length two.

The witnesses certify the incidence geometry and nonemptiness of the listed
clean algebraic opens.  They do not compute affine complement groups or prove
Whitney--Thom propagation.  Compatible residual coefficient-rank loci and
the ``P``-critical triple fibers at ``k=+/-2`` remain deeper boundary audits;
no plane Jacobian-conjecture claim is made here.
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
    gcd,
    linear_eq_to_matrix,
    resultant,
)

from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    CUSP_IMAGE_FACTOR,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    PAIR_INCIDENCE,
    R,
    S,
    SECOND_DIVIDED_DIFFERENCE,
    T,
    ZERO_FIBER_GRAPH,
    ZERO_FIBER_VERTICAL,
    exceptional_graph,
    exceptional_vertical,
)

PARAMETERS: Final = (ALPHA, BETA, GAMMA, DELTA)
FIRST_ROOT: Final = Symbol("u_split")
SECOND_ROOT: Final = Symbol("v_split")
NODE_TARGET_X: Final = Symbol("x_split_node")


def split_vertical(kappa: int) -> Expr:
    """Return the true vertical collision component at a split ``k``."""

    if kappa == 0:
        return ZERO_FIBER_VERTICAL
    if kappa in (-2, 2):
        return exceptional_vertical(kappa // 2)
    msg = "the split parameter must be 0, +2, or -2"
    raise ValueError(msg)


def split_graph(kappa: int) -> Expr:
    """Return the true graph collision component at a split ``k``."""

    if kappa == 0:
        return ZERO_FIBER_GRAPH
    if kappa in (-2, 2):
        return exceptional_graph(kappa // 2)
    msg = "the split parameter must be 0, +2, or -2"
    raise ValueError(msg)


def vertical_target_x(kappa: int) -> Expr:
    """Return the common first target coordinate on the vertical component."""

    if kappa == 0:
        return R * (R - 1)
    if kappa in (-2, 2):
        return R**2
    msg = "the split parameter must be 0, +2, or -2"
    raise ValueError(msg)


def graph_target_x(kappa: int) -> Expr:
    """Return the common first target coordinate on the graph component."""

    if kappa == 0:
        return -((S**2 + 1) ** 2) / 4
    if kappa == 2:
        return -(S * (S + 1) ** 2 * (S + 2)) / 4
    if kappa == -2:
        return -(S * (S - 1) ** 2 * (S - 2)) / 4
    msg = "the split parameter must be 0, +2, or -2"
    raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class SplitChartCertificate:
    """Exact identities and boundary data for one reducible pair chart."""

    kappa: int
    pair_factor_identity: Expr
    vertical_restriction_identity: Expr
    graph_restriction_identity: Expr
    vertical_degree: int
    graph_degree: int
    overlap_identity: Expr
    overlap_condition: Expr
    overlap_cusp_identity: Expr
    vertical_diagonal_critical_identity: Expr
    graph_diagonal_critical_identity: Expr
    discriminant_degrees: tuple[int, int]
    discriminants_squarefree: tuple[bool, bool]
    discriminant_gcd_degree: int

    @property
    def verified(self) -> bool:
        """Whether the split component and boundary audit agrees exactly."""

        expected_degrees = (2, 8) if self.kappa == 0 else (4, 6)
        expected_discriminants = (2, 13) if self.kappa == 0 else (6, 9)
        return bool(
            self.pair_factor_identity == 0
            and self.vertical_restriction_identity == 0
            and self.graph_restriction_identity == 0
            and (self.vertical_degree, self.graph_degree) == expected_degrees
            and self.overlap_identity == 0
            and self.overlap_cusp_identity == 0
            and self.vertical_diagonal_critical_identity == 0
            and self.graph_diagonal_critical_identity == 0
            and self.discriminant_degrees == expected_discriminants
            and all(self.discriminants_squarefree)
            and self.discriminant_gcd_degree == 0
        )


@cache
def exact_split_chart_certificates() -> tuple[SplitChartCertificate, ...]:
    """Build the exact ``V/W`` certificates for ``k=0,+2,-2``."""

    certificates: list[SplitChartCertificate] = []
    for kappa in (0, 2, -2):
        vertical = split_vertical(kappa)
        graph = split_graph(kappa)
        if kappa == 0:
            pair_factor = S * (2 * R - S**2 - 1)
            vertical_restriction = SECOND_DIVIDED_DIFFERENCE.subs({KAPPA: 0, S: 0})
            expected_vertical = R**2 * vertical
            graph_restriction = 16 * SECOND_DIVIDED_DIFFERENCE.subs(
                {KAPPA: 0, R: (S**2 + 1) / 2}
            )
            expected_graph = graph
            overlap_identity = graph.subs(S, 0) - 4 * vertical.subs(R, Rational(1, 2))
            overlap_condition = expand(graph.subs(S, 0))
            overlap_cusp_identity = 0
            vertical_diagonal_critical_identity = 0
            graph_diagonal_critical_identity = expand(
                resultant(graph, S**2 + 2, S)
                - EXTRA_CRITICAL_FACTOR.subs(KAPPA, 0) / 16
            )
        else:
            epsilon = kappa // 2
            pair_factor = (S + epsilon) * (2 * R - S * (S + epsilon))
            vertical_restriction = SECOND_DIVIDED_DIFFERENCE.subs(
                {KAPPA: kappa, S: -epsilon}
            )
            expected_vertical = vertical
            graph_restriction = 16 * SECOND_DIVIDED_DIFFERENCE.subs(
                {KAPPA: kappa, R: S * (S + epsilon) / 2}
            )
            expected_graph = S**2 * graph
            overlap_identity = graph.subs(S, -epsilon) - 16 * vertical.subs(R, 0)
            overlap_condition = expand(vertical.subs(R, 0))
            overlap_cusp_identity = expand(
                CUSP_IMAGE_FACTOR.subs(KAPPA, kappa) - overlap_condition**2
            )
            vertical_critical_form = expand(256 * vertical.subs(R, Rational(1, 4)))
            graph_critical_form = expand(graph.subs(S, -2 * epsilon) / 4)
            vertical_diagonal_critical_identity = expand(
                EXTRA_CRITICAL_FACTOR.subs(KAPPA, kappa)
                - 16 * vertical_critical_form * graph_critical_form
            )
            graph_diagonal_critical_identity = 0

        vertical_discriminant = Poly(
            discriminant(vertical, R),
            *PARAMETERS,
        )
        graph_discriminant = Poly(
            discriminant(graph, S),
            *PARAMETERS,
        )
        certificates.append(
            SplitChartCertificate(
                kappa=kappa,
                pair_factor_identity=expand(
                    PAIR_INCIDENCE.subs(KAPPA, kappa) - pair_factor
                ),
                vertical_restriction_identity=expand(
                    vertical_restriction - expected_vertical
                ),
                graph_restriction_identity=expand(graph_restriction - expected_graph),
                vertical_degree=Poly(vertical, R).degree(),
                graph_degree=Poly(graph, S).degree(),
                overlap_identity=expand(overlap_identity),
                overlap_condition=overlap_condition,
                overlap_cusp_identity=expand(overlap_cusp_identity),
                vertical_diagonal_critical_identity=expand(
                    vertical_diagonal_critical_identity
                ),
                graph_diagonal_critical_identity=expand(
                    graph_diagonal_critical_identity
                ),
                discriminant_degrees=(
                    vertical_discriminant.total_degree(),
                    graph_discriminant.total_degree(),
                ),
                discriminants_squarefree=(
                    vertical_discriminant.sqf_part().total_degree()
                    == vertical_discriminant.total_degree(),
                    graph_discriminant.sqf_part().total_degree()
                    == graph_discriminant.total_degree(),
                ),
                discriminant_gcd_degree=vertical_discriminant.gcd(
                    graph_discriminant
                ).total_degree(),
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class SplitAllocation:
    """One allowed component-length allocation for a collision profile."""

    profile: str
    split_class: str
    allocation: str
    special_vertical_length: int
    special_graph_length: int
    vertical_node_count: int
    graph_node_count: int
    expected_dimension: int

    @property
    def total_delta(self) -> int:
        """Return the split collision length, including residual nodes."""

        return (
            self.special_vertical_length
            + self.special_graph_length
            + self.vertical_node_count
            + self.graph_node_count
        )

    @property
    def node_count(self) -> int:
        """Return the residual ordinary-node count."""

        return self.vertical_node_count + self.graph_node_count

    @property
    def verified(self) -> bool:
        """Whether this row fills the correct split degree budget."""

        vertical_budget, graph_budget = (2, 8) if self.split_class == "k=0" else (4, 6)
        return bool(
            self.special_vertical_length + self.vertical_node_count == vertical_budget
            and self.special_graph_length + self.graph_node_count == graph_budget
            and self.total_delta == 10
            and self.expected_dimension == 2
        )


def _allocation(
    profile: str,
    split_class: str,
    allocation: str,
    special: tuple[int, int],
    nodes: tuple[int, int],
) -> SplitAllocation:
    return SplitAllocation(
        profile=profile,
        split_class=split_class,
        allocation=allocation,
        special_vertical_length=special[0],
        special_graph_length=special[1],
        vertical_node_count=nodes[0],
        graph_node_count=nodes[1],
        expected_dimension=2,
    )


SPLIT_ALLOCATIONS: Final = (
    # k=0.  A non-overlap C3 cannot lie on the monic quadratic V component.
    _allocation("C3+7N", "k=0", "W", (0, 3), (2, 5)),
    _allocation("C3+7N", "k=0", "overlap-V", (2, 1), (0, 7)),
    _allocation("C3+7N", "k=0", "overlap-W", (1, 2), (1, 6)),
    # Two distinct double roots cannot both fit in the degree-two V component.
    _allocation("C2^2+6N", "k=0", "VW", (2, 2), (0, 6)),
    _allocation("C2^2+6N", "k=0", "WW", (0, 4), (2, 4)),
    _allocation("C2^2+6N", "k=0", "overlap+W", (1, 3), (1, 5)),
    # A generic triple fiber has one V edge and two W edges.
    _allocation("T112+6N", "k=0", "tangent-V", (2, 2), (0, 6)),
    _allocation("T112+6N", "k=0", "tangent-W", (1, 3), (1, 5)),
    # The triple already consumes one V edge, so a separate V contact does not
    # fit in the degree-two V budget.
    _allocation("C2+T111+5N", "k=0", "contact-W", (1, 4), (1, 4)),
    _allocation("Q0+4N", "k=0", "2V+4W", (2, 4), (0, 4)),
    _allocation("T111^2+4N", "k=0", "2V+4W", (2, 4), (0, 4)),
    # k=+/-2.  Both component budgets support every generic edge type.
    _allocation("C3+7N", "k=+/-2", "V", (3, 0), (1, 6)),
    _allocation("C3+7N", "k=+/-2", "W", (0, 3), (4, 3)),
    _allocation("C2^2+6N", "k=+/-2", "VV", (4, 0), (0, 6)),
    _allocation("C2^2+6N", "k=+/-2", "VW", (2, 2), (2, 4)),
    _allocation("C2^2+6N", "k=+/-2", "WW", (0, 4), (4, 2)),
    _allocation("T112+6N", "k=+/-2", "tangent-V", (2, 2), (2, 4)),
    _allocation("T112+6N", "k=+/-2", "tangent-W", (1, 3), (3, 3)),
    _allocation("C2+T111+5N", "k=+/-2", "contact-V", (3, 2), (1, 4)),
    _allocation("C2+T111+5N", "k=+/-2", "contact-W", (1, 4), (3, 2)),
    _allocation("Q0+4N", "k=+/-2", "2V+4W", (2, 4), (2, 2)),
    _allocation("T111^2+4N", "k=+/-2", "2V+4W", (2, 4), (2, 2)),
)


FORBIDDEN_K0_ALLOCATIONS: Final = (
    ("C3+7N", "V", "the monic quadratic V cannot have a triple root"),
    (
        "C2^2+6N",
        "VV",
        "two distinct double roots require V-degree at least four",
    ),
    (
        "C2+T111+5N",
        "contact-V",
        "the triple V edge plus a separate double V edge has length three",
    ),
)


@dataclass(frozen=True, slots=True)
class RootMultiplicity:
    """A prescribed root and exact multiplicity on ``V`` or ``W``."""

    component: str
    root: Expr
    multiplicity: int


@dataclass(frozen=True, slots=True)
class SplitWitnessSpec:
    """Exact linear incidence data for one allowed allocation."""

    name: str
    profile: str
    kappa: int
    allocation: str
    equations: tuple[Expr, ...]
    coefficient_values: tuple[Expr, Expr, Expr, Expr]
    base_dimension: int
    expected_rank: int
    node_count: int
    roots: tuple[RootMultiplicity, ...]
    special_target_x_values: tuple[Expr, ...]


def _q_equal(left: Expr, right: Expr) -> Expr:
    return expand(FAMILY_Q.subs(T, left) - FAMILY_Q.subs(T, right))


def _tangent(kappa: int, left: Expr, right: Expr) -> Expr:
    p_derivative = diff(FAMILY_P.subs(KAPPA, kappa), T)
    q_derivative = diff(FAMILY_Q, T)
    return expand(
        q_derivative.subs(T, left) * p_derivative.subs(T, right)
        - q_derivative.subs(T, right) * p_derivative.subs(T, left)
    )


def _component_contact_equations(
    polynomial: Expr,
    variable: Expr,
    root: Expr,
    multiplicity: int,
) -> tuple[Expr, ...]:
    return tuple(
        diff(polynomial, variable, order).subs(variable, root)
        for order in range(multiplicity)
    )


K0_FIBER_ONE: Final = (-4 * I / 5, 4 * I / 5, 3 * I / 5, -3 * I / 5)
K0_FIBER_TWO: Final = (-12 * I / 13, 12 * I / 13, 5 * I / 13, -5 * I / 13)
K2_FIBER_ONE: Final = (
    Rational(1, 5),
    Rational(-6, 5),
    Rational(-2, 5),
    Rational(-3, 5),
)
K2_FIBER_TWO: Final = (
    Rational(2, 13),
    Rational(-15, 13),
    Rational(-3, 13),
    Rational(-10, 13),
)


def _triple_equalities(roots: tuple[Expr, Expr, Expr, Expr]) -> tuple[Expr, Expr]:
    return (_q_equal(roots[0], roots[1]), _q_equal(roots[0], roots[2]))


def _quadruple_equalities(
    roots: tuple[Expr, Expr, Expr, Expr],
) -> tuple[Expr, Expr, Expr]:
    return (
        _q_equal(roots[0], roots[1]),
        _q_equal(roots[0], roots[2]),
        _q_equal(roots[0], roots[3]),
    )


def _coefficients(
    alpha: Expr | int,
    beta: Expr | int,
    gamma: Expr | int,
    delta: Expr | int,
) -> tuple[Expr, Expr, Expr, Expr]:
    return (Rational(alpha), Rational(beta), Rational(gamma), Rational(delta))


@cache
def split_witness_specs() -> tuple[SplitWitnessSpec, ...]:
    """Return exact clean witnesses for every allowed allocation."""

    specs: list[SplitWitnessSpec] = []

    def add(
        name: str,
        profile: str,
        kappa: int,
        allocation: str,
        equations: tuple[Expr, ...],
        coefficients: tuple[Expr, Expr, Expr, Expr],
        base_dimension: int,
        expected_rank: int,
        node_count: int,
        roots: tuple[RootMultiplicity, ...],
        special_x: tuple[Expr, ...],
    ) -> None:
        specs.append(
            SplitWitnessSpec(
                name=name,
                profile=profile,
                kappa=kappa,
                allocation=allocation,
                equations=equations,
                coefficient_values=coefficients,
                base_dimension=base_dimension,
                expected_rank=expected_rank,
                node_count=node_count,
                roots=roots,
                special_target_x_values=special_x,
            )
        )

    v0 = ZERO_FIBER_VERTICAL
    w0 = ZERO_FIBER_GRAPH
    v2 = exceptional_vertical(1)
    w2 = exceptional_graph(1)

    # C3 witnesses, including the two branches of the legitimate k=0 overlap.
    add(
        "c3_k0_w",
        "C3+7N",
        0,
        "W",
        _component_contact_equations(w0, S, 1, 3),
        _coefficients(13, Rational(-7, 2), 12, 1),
        1,
        3,
        7,
        (RootMultiplicity("W", Rational(1), 3),),
        (graph_target_x(0).subs(S, 1),),
    )
    overlap = expand(w0.subs(S, 0))
    add(
        "c3_k0_overlap_v",
        "C3+7N",
        0,
        "overlap-V",
        (overlap, diff(v0, R).subs(R, Rational(1, 2))),
        _coefficients(Rational(1, 4), 1, 1, 0),
        0,
        2,
        7,
        (
            RootMultiplicity("V", Rational(1, 2), 2),
            RootMultiplicity("W", Rational(0), 1),
        ),
        (Rational(-1, 4),),
    )
    add(
        "c3_k0_overlap_w",
        "C3+7N",
        0,
        "overlap-W",
        (overlap, diff(w0, S).subs(S, 0)),
        _coefficients(Rational(-1, 4), 2, 0, 3),
        0,
        2,
        7,
        (
            RootMultiplicity("V", Rational(1, 2), 1),
            RootMultiplicity("W", Rational(0), 2),
        ),
        (Rational(-1, 4),),
    )
    for allocation, polynomial, variable, coefficients in (
        ("V", v2, R, _coefficients(-13, -22, -14, -1)),
        ("W", w2, S, _coefficients(27, 15, 27, 0)),
    ):
        add(
            f"c3_k2_{allocation.lower()}",
            "C3+7N",
            2,
            allocation,
            _component_contact_equations(polynomial, variable, 1, 3),
            coefficients,
            1,
            3,
            7,
            (RootMultiplicity(allocation, Rational(1), 3),),
            (
                (vertical_target_x(2).subs(R, 1))
                if allocation == "V"
                else (graph_target_x(2).subs(S, 1)),
            ),
        )

    # C2^2 witnesses.
    pair_specs = (
        (
            "c22_k0_vw",
            0,
            "VW",
            v0,
            R,
            Rational(1),
            w0,
            S,
            Rational(1),
            _coefficients(1, Rational(-1, 2), 2, -1),
        ),
        (
            "c22_k0_ww",
            0,
            "WW",
            w0,
            S,
            Rational(1),
            w0,
            S,
            Rational(2),
            _coefficients(
                Rational(3981, 304),
                Rational(-261, 152),
                Rational(87, 8),
                Rational(675, 304),
            ),
        ),
        (
            "c22_k2_vv",
            2,
            "VV",
            v2,
            R,
            Rational(1),
            v2,
            R,
            Rational(2),
            _coefficients(17, 31, 24, 7),
        ),
        (
            "c22_k2_vw",
            2,
            "VW",
            v2,
            R,
            Rational(1),
            w2,
            S,
            Rational(1),
            _coefficients(3, 3, 3, 0),
        ),
        (
            "c22_k2_ww",
            2,
            "WW",
            w2,
            S,
            Rational(1),
            w2,
            S,
            Rational(2),
            _coefficients(
                Rational(397, 13),
                Rational(256, 13),
                Rational(378, 13),
                Rational(19, 13),
            ),
        ),
    )
    for (
        name,
        kappa,
        allocation,
        first_polynomial,
        first_variable,
        first_value,
        second_polynomial,
        second_variable,
        second_value,
        coefficients,
    ) in pair_specs:
        first_component = allocation[0]
        second_component = allocation[1]
        first_x = (
            vertical_target_x(kappa).subs(R, first_value)
            if first_component == "V"
            else graph_target_x(kappa).subs(S, first_value)
        )
        second_x = (
            vertical_target_x(kappa).subs(R, second_value)
            if second_component == "V"
            else graph_target_x(kappa).subs(S, second_value)
        )
        add(
            name,
            "C2^2+6N",
            kappa,
            allocation,
            (
                *_component_contact_equations(
                    first_polynomial, first_variable, first_value, 2
                ),
                *_component_contact_equations(
                    second_polynomial, second_variable, second_value, 2
                ),
            ),
            coefficients,
            2,
            4,
            6,
            (
                RootMultiplicity(first_component, first_value, 2),
                RootMultiplicity(second_component, second_value, 2),
            ),
            (first_x, second_x),
        )
    add(
        "c22_k0_overlap_w",
        "C2^2+6N",
        0,
        "overlap+W",
        (
            overlap,
            *_component_contact_equations(w0, S, Rational(1), 2),
        ),
        _coefficients(Rational(-3, 2), Rational(15, 4), Rational(-5, 2), 1),
        1,
        3,
        6,
        (
            RootMultiplicity("V", Rational(1, 2), 1),
            RootMultiplicity("W", Rational(0), 1),
            RootMultiplicity("W", Rational(1), 2),
        ),
        (Rational(-1, 4), graph_target_x(0).subs(S, 1)),
    )

    # T112 witnesses.  The first two entries in each fiber form the V pair;
    # the first and third form the selected W pair.
    t112_data = (
        (
            "t112_k0_v",
            0,
            "tangent-V",
            K0_FIBER_ONE,
            (K0_FIBER_ONE[0], K0_FIBER_ONE[1]),
            (
                Rational(256, 625),
                -Rational(1701, 60125) * I,
                Rational(32, 25),
                Rational(0),
            ),
            (
                RootMultiplicity("V", Rational(16, 25), 2),
                RootMultiplicity("W", -I / 5, 1),
                RootMultiplicity("W", 7 * I / 5, 1),
            ),
        ),
        (
            "t112_k0_w",
            0,
            "tangent-W",
            K0_FIBER_ONE,
            (K0_FIBER_ONE[0], K0_FIBER_ONE[2]),
            (
                Rational(1554736, 4329625),
                -Rational(87966, 4329625) * I,
                Rational(1040047, 865925),
                Rational(0),
            ),
            (
                RootMultiplicity("V", Rational(16, 25), 1),
                RootMultiplicity("W", -I / 5, 2),
                RootMultiplicity("W", 7 * I / 5, 1),
            ),
        ),
        (
            "t112_k2_v",
            2,
            "tangent-V",
            K2_FIBER_ONE,
            (K2_FIBER_ONE[0], K2_FIBER_ONE[1]),
            (
                -Rational(151971228, 107591875),
                -Rational(124855138, 21518375),
                -Rational(22783989, 4303675),
                Rational(0),
            ),
            (
                RootMultiplicity("V", -Rational(6, 25), 2),
                RootMultiplicity("W", -Rational(1, 5), 1),
                RootMultiplicity("W", -Rational(8, 5), 1),
            ),
        ),
        (
            "t112_k2_w",
            2,
            "tangent-W",
            K2_FIBER_ONE,
            (K2_FIBER_ONE[0], K2_FIBER_ONE[2]),
            (
                -Rational(438148, 1705125),
                -Rational(2910662, 1705125),
                -Rational(915311, 341025),
                Rational(0),
            ),
            (
                RootMultiplicity("V", -Rational(6, 25), 1),
                RootMultiplicity("W", -Rational(1, 5), 2),
                RootMultiplicity("W", -Rational(8, 5), 1),
            ),
        ),
    )
    for (
        name,
        kappa,
        allocation,
        fiber,
        tangent_pair,
        coefficients,
        t112_roots,
    ) in t112_data:
        triple_x = expand(FAMILY_P.subs({KAPPA: kappa, T: fiber[0]}))
        add(
            name,
            "T112+6N",
            kappa,
            allocation,
            (*_triple_equalities(fiber), _tangent(kappa, *tangent_pair)),
            coefficients,
            1,
            3,
            6,
            t112_roots,
            (triple_x,),
        )

    # Mixed contact plus ordinary triple.
    mixed_data = (
        (
            "mixed_k0_w",
            0,
            "contact-W",
            K0_FIBER_ONE,
            (-12 * I / 13, 5 * I / 13),
            (
                Rational(19324528, 79493375),
                Rational(3533716867, 7233897125) * I,
                Rational(16214067, 15898675),
                Rational(389566783, 556453625) * I,
            ),
            (
                RootMultiplicity("V", Rational(16, 25), 1),
                RootMultiplicity("W", -I / 5, 1),
                RootMultiplicity("W", 7 * I / 5, 1),
                RootMultiplicity("W", -7 * I / 13, 2),
            ),
        ),
        (
            "mixed_k2_v",
            2,
            "contact-V",
            K2_FIBER_ONE,
            (Rational(2), Rational(-3)),
            (
                Rational(824972, 40625),
                Rational(690283, 8125),
                Rational(717603, 8125),
                Rational(223897, 8125),
            ),
            (
                RootMultiplicity("V", -Rational(6, 25), 1),
                RootMultiplicity("W", -Rational(1, 5), 1),
                RootMultiplicity("W", -Rational(8, 5), 1),
                RootMultiplicity("V", Rational(-6), 2),
            ),
        ),
        (
            "mixed_k2_w",
            2,
            "contact-W",
            K2_FIBER_ONE,
            (Rational(2, 13), Rational(-3, 13)),
            (
                Rational(114404, 31581875),
                Rational(1599289, 6316375),
                Rational(7929321, 6316375),
                Rational(13090099, 6316375),
            ),
            (
                RootMultiplicity("V", -Rational(6, 25), 1),
                RootMultiplicity("W", -Rational(1, 5), 1),
                RootMultiplicity("W", -Rational(8, 5), 1),
                RootMultiplicity("W", -Rational(1, 13), 2),
            ),
        ),
    )
    for (
        name,
        kappa,
        allocation,
        fiber,
        contact_pair,
        coefficients,
        mixed_roots,
    ) in mixed_data:
        triple_x = expand(FAMILY_P.subs({KAPPA: kappa, T: fiber[0]}))
        contact_x = expand(FAMILY_P.subs({KAPPA: kappa, T: contact_pair[0]}))
        add(
            name,
            "C2+T111+5N",
            kappa,
            allocation,
            (
                *_triple_equalities(fiber),
                _q_equal(*contact_pair),
                _tangent(kappa, *contact_pair),
            ),
            coefficients,
            2,
            4,
            5,
            mixed_roots,
            (triple_x, contact_x),
        )

    # Ordinary quadruple witnesses.
    q0_data = (
        (
            "q0_k0",
            0,
            K0_FIBER_ONE,
            _coefficients(Rational(144, 625), Rational(337, 481), 1, 1),
            (
                RootMultiplicity("V", Rational(16, 25), 1),
                RootMultiplicity("V", Rational(9, 25), 1),
                RootMultiplicity("W", -I / 5, 1),
                RootMultiplicity("W", -7 * I / 5, 1),
                RootMultiplicity("W", 7 * I / 5, 1),
                RootMultiplicity("W", I / 5, 1),
            ),
        ),
        (
            "q0_k2",
            2,
            K2_FIBER_ONE,
            _coefficients(
                -Rational(203796, 413125),
                -Rational(1682, 661),
                -Rational(2127, 661),
                0,
            ),
            (
                RootMultiplicity("V", -Rational(6, 25), 1),
                RootMultiplicity("V", Rational(6, 25), 1),
                RootMultiplicity("W", -Rational(1, 5), 1),
                RootMultiplicity("W", -Rational(2, 5), 1),
                RootMultiplicity("W", -Rational(8, 5), 1),
                RootMultiplicity("W", -Rational(9, 5), 1),
            ),
        ),
    )
    for name, kappa, fiber, coefficients, q0_roots in q0_data:
        target_x = expand(FAMILY_P.subs({KAPPA: kappa, T: fiber[0]}))
        add(
            name,
            "Q0+4N",
            kappa,
            "2V+4W",
            _quadruple_equalities(fiber),
            coefficients,
            1,
            3,
            4,
            q0_roots,
            (target_x,),
        )

    # Two ordinary triple fibers.
    double_triple_data = (
        (
            "double_triple_k0",
            0,
            K0_FIBER_ONE,
            K0_FIBER_TWO,
            (
                Rational(2304, 4225),
                -Rational(4443841, 17576000) * I,
                Rational(6304, 4225),
                -Rational(391957, 1352000) * I,
            ),
            (
                RootMultiplicity("V", Rational(16, 25), 1),
                RootMultiplicity("V", Rational(144, 169), 1),
                RootMultiplicity("W", -I / 5, 1),
                RootMultiplicity("W", 7 * I / 5, 1),
                RootMultiplicity("W", -7 * I / 13, 1),
                RootMultiplicity("W", 17 * I / 13, 1),
            ),
        ),
        (
            "double_triple_k2",
            2,
            K2_FIBER_ONE,
            K2_FIBER_TWO,
            (
                Rational(1817532, 40919125),
                Rational(262163, 327353),
                Rational(108449379, 40919125),
                Rational(117773061, 40919125),
            ),
            (
                RootMultiplicity("V", -Rational(6, 25), 1),
                RootMultiplicity("V", -Rational(30, 169), 1),
                RootMultiplicity("W", -Rational(1, 5), 1),
                RootMultiplicity("W", -Rational(8, 5), 1),
                RootMultiplicity("W", -Rational(1, 13), 1),
                RootMultiplicity("W", -Rational(18, 13), 1),
            ),
        ),
    )
    for (
        name,
        kappa,
        first_fiber,
        second_fiber,
        coefficients,
        double_triple_roots,
    ) in double_triple_data:
        first_x = expand(FAMILY_P.subs({KAPPA: kappa, T: first_fiber[0]}))
        second_x = expand(FAMILY_P.subs({KAPPA: kappa, T: second_fiber[0]}))
        add(
            name,
            "T111^2+4N",
            kappa,
            "2V+4W",
            (*_triple_equalities(first_fiber), *_triple_equalities(second_fiber)),
            coefficients,
            2,
            4,
            4,
            double_triple_roots,
            (first_x, second_x),
        )

    return tuple(specs)


def _remove_prescribed_roots(
    polynomial: Expr,
    variable: Expr,
    roots: tuple[RootMultiplicity, ...],
    component: str,
) -> tuple[Expr, tuple[Expr, ...]]:
    residual = polynomial
    order_checks: list[Expr] = []
    for root_data in roots:
        if root_data.component != component:
            continue
        for order in range(root_data.multiplicity):
            order_checks.append(
                expand(diff(polynomial, variable, order).subs(variable, root_data.root))
            )
        order_checks.append(
            expand(
                diff(polynomial, variable, root_data.multiplicity).subs(
                    variable, root_data.root
                )
            )
        )
        residual = cancel(
            residual / (variable - root_data.root) ** root_data.multiplicity
        )
    return expand(residual), tuple(order_checks)


def _target_eliminant(
    residual: Expr,
    variable: Expr,
    target_expression: Expr,
) -> Poly:
    residual_poly = Poly(residual, variable)
    if residual_poly.degree() == 0:
        return Poly(1, NODE_TARGET_X)
    return Poly(
        resultant(
            residual_poly.as_expr(),
            NODE_TARGET_X - target_expression,
            variable,
        ),
        NODE_TARGET_X,
    )


@dataclass(frozen=True, slots=True)
class SplitWitnessCertificate:
    """Exact rank, multiplicity, boundary, and node-separation evidence."""

    name: str
    profile: str
    kappa: int
    allocation: str
    coefficient_rank: int
    augmented_rank: int
    base_dimension: int
    free_coefficient_dimension: int
    incidence_dimension: int
    equation_residuals: tuple[Expr, ...]
    alpha_nonzero: bool
    vertical_order_checks: tuple[Expr, ...]
    graph_order_checks: tuple[Expr, ...]
    vertical_residual_degree: int
    graph_residual_degree: int
    residual_squarefree: tuple[bool, bool]
    node_target_degrees: tuple[int, int]
    node_targets_squarefree: tuple[bool, bool]
    cross_node_target_gcd_degree: int
    special_targets_distinct: bool
    special_node_separations_nonzero: bool
    cusp_image_nonzero: bool
    extra_critical_nonzero: bool
    overlap_clean: bool

    @staticmethod
    def _orders_exact(
        checks: tuple[Expr, ...], roots: tuple[RootMultiplicity, ...], component: str
    ) -> bool:
        position = 0
        for root in roots:
            if root.component != component:
                continue
            values = checks[position : position + root.multiplicity + 1]
            if any(value != 0 for value in values[:-1]) or values[-1] == 0:
                return False
            position += root.multiplicity + 1
        return position == len(checks)

    def verified(self, spec: SplitWitnessSpec) -> bool:
        """Whether every exact clean-incidence check agrees."""

        return bool(
            self.name == spec.name
            and self.profile == spec.profile
            and self.kappa == spec.kappa
            and self.allocation == spec.allocation
            and self.coefficient_rank == spec.expected_rank
            and self.augmented_rank == spec.expected_rank
            and self.free_coefficient_dimension == 4 - spec.expected_rank
            and self.incidence_dimension == 2
            and all(value == 0 for value in self.equation_residuals)
            and self.alpha_nonzero
            and self._orders_exact(self.vertical_order_checks, spec.roots, "V")
            and self._orders_exact(self.graph_order_checks, spec.roots, "W")
            and self.vertical_residual_degree + self.graph_residual_degree
            == spec.node_count
            and all(self.residual_squarefree)
            and self.node_target_degrees
            == (self.vertical_residual_degree, self.graph_residual_degree)
            and all(self.node_targets_squarefree)
            and self.cross_node_target_gcd_degree == 0
            and self.special_targets_distinct
            and self.special_node_separations_nonzero
            and self.cusp_image_nonzero
            and self.extra_critical_nonzero
            and self.overlap_clean
        )


@cache
def exact_split_witness_certificates() -> tuple[SplitWitnessCertificate, ...]:
    """Build exact clean witnesses for all allowed split allocations."""

    certificates: list[SplitWitnessCertificate] = []
    for spec in split_witness_specs():
        matrix, right_side = linear_eq_to_matrix(spec.equations, PARAMETERS)
        augmented = matrix.row_join(right_side)
        substitution = dict(zip(PARAMETERS, spec.coefficient_values, strict=True))
        vertical = expand(split_vertical(spec.kappa).subs(substitution))
        graph = expand(split_graph(spec.kappa).subs(substitution))
        vertical_residual, vertical_orders = _remove_prescribed_roots(
            vertical,
            R,
            spec.roots,
            "V",
        )
        graph_residual, graph_orders = _remove_prescribed_roots(
            graph,
            S,
            spec.roots,
            "W",
        )
        vertical_residual_poly = Poly(vertical_residual, R)
        graph_residual_poly = Poly(graph_residual, S)
        vertical_nodes = _target_eliminant(
            vertical_residual,
            R,
            vertical_target_x(spec.kappa),
        )
        graph_nodes = _target_eliminant(
            graph_residual,
            S,
            graph_target_x(spec.kappa),
        )
        special_values = tuple(expand(value) for value in spec.special_target_x_values)
        special_targets_distinct = len(set(special_values)) == len(special_values)
        special_node_separations_nonzero = all(
            vertical_nodes.as_expr().subs(NODE_TARGET_X, target) != 0
            and graph_nodes.as_expr().subs(NODE_TARGET_X, target) != 0
            for target in special_values
        )
        if spec.kappa == 0:
            overlap_is_prescribed = any(
                root.component == "V" and root.root == Rational(1, 2)
                for root in spec.roots
            ) and any(root.component == "W" and root.root == 0 for root in spec.roots)
            overlap_value = expand(vertical.subs(R, Rational(1, 2)))
            overlap_clean = overlap_is_prescribed or overlap_value != 0
        else:
            overlap_clean = (
                split_vertical(spec.kappa).subs(R, 0).subs(substitution) != 0
            )
        certificates.append(
            SplitWitnessCertificate(
                name=spec.name,
                profile=spec.profile,
                kappa=spec.kappa,
                allocation=spec.allocation,
                coefficient_rank=matrix.rank(),
                augmented_rank=augmented.rank(),
                base_dimension=spec.base_dimension,
                free_coefficient_dimension=4 - matrix.rank(),
                incidence_dimension=spec.base_dimension + 4 - matrix.rank(),
                equation_residuals=tuple(
                    expand(equation.subs(substitution)) for equation in spec.equations
                ),
                alpha_nonzero=substitution[ALPHA] != 0,
                vertical_order_checks=vertical_orders,
                graph_order_checks=graph_orders,
                vertical_residual_degree=vertical_residual_poly.degree(),
                graph_residual_degree=graph_residual_poly.degree(),
                residual_squarefree=(
                    vertical_residual_poly.sqf_part().degree()
                    == vertical_residual_poly.degree(),
                    graph_residual_poly.sqf_part().degree()
                    == graph_residual_poly.degree(),
                ),
                node_target_degrees=(vertical_nodes.degree(), graph_nodes.degree()),
                node_targets_squarefree=(
                    vertical_nodes.sqf_part().degree() == vertical_nodes.degree(),
                    graph_nodes.sqf_part().degree() == graph_nodes.degree(),
                ),
                cross_node_target_gcd_degree=vertical_nodes.gcd(graph_nodes).degree(),
                special_targets_distinct=special_targets_distinct,
                special_node_separations_nonzero=special_node_separations_nonzero,
                cusp_image_nonzero=CUSP_IMAGE_FACTOR.subs(KAPPA, spec.kappa).subs(
                    substitution
                )
                != 0,
                extra_critical_nonzero=EXTRA_CRITICAL_FACTOR.subs(
                    KAPPA, spec.kappa
                ).subs(substitution)
                != 0,
                overlap_clean=bool(overlap_clean),
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class SplitRankBoundaryCertificate:
    """Generic rank and forbidden-allocation evidence on component charts."""

    k0_vertical_c3_ranks: tuple[int, int]
    k0_graph_c3_minor_gcd: Expr
    k2_vertical_c3_minor_gcd: Expr
    k2_graph_c3_minor_gcd: Expr
    k0_vw_double_contact_determinant: Expr
    k0_vv_double_contact_determinant: Expr
    k2_vv_double_contact_determinant: Expr
    plus_minus_vertical_symmetry: Expr
    plus_minus_graph_symmetry: Expr
    k0_overlap_value_identity: Expr
    k0_overlap_vertical_derivative: Expr
    k0_overlap_graph_derivative: Expr

    @property
    def verified(self) -> bool:
        """Whether the generic-rank and overlap audit agrees."""

        return bool(
            self.k0_vertical_c3_ranks == (2, 3)
            and expand(self.k0_graph_c3_minor_gcd - 64 * (S**2 + 1) ** 3) == 0
            and self.k2_vertical_c3_minor_gcd == 2
            and self.k2_graph_c3_minor_gcd == 64
            and expand(
                self.k0_vw_double_contact_determinant
                - 256 * SECOND_ROOT**3 * (SECOND_ROOT**2 + 1)
            )
            == 0
            and self.k0_vv_double_contact_determinant == 0
            and expand(
                self.k2_vv_double_contact_determinant - (FIRST_ROOT - SECOND_ROOT) ** 4
            )
            == 0
            and self.plus_minus_vertical_symmetry == 0
            and self.plus_minus_graph_symmetry == 0
            and self.k0_overlap_value_identity == 0
            and self.k0_overlap_vertical_derivative == 1 - GAMMA
            and self.k0_overlap_graph_derivative == 4 * (3 * BETA - 2 * DELTA)
        )


def _maximal_minor_gcd(matrix: Matrix) -> Expr:
    minors = [
        expand(matrix[:, columns].det())
        for columns in combinations(range(matrix.cols), matrix.rows)
    ]
    current = minors[0]
    for minor in minors[1:]:
        current = gcd(current, minor)
    return expand(current)


def _contact_matrix(
    first: tuple[Expr, Expr, Expr],
    second: tuple[Expr, Expr, Expr],
) -> Matrix:
    first_polynomial, first_variable, first_root = first
    second_polynomial, second_variable, second_root = second
    equations = (
        first_polynomial.subs(first_variable, first_root),
        diff(first_polynomial, first_variable).subs(first_variable, first_root),
        second_polynomial.subs(second_variable, second_root),
        diff(second_polynomial, second_variable).subs(second_variable, second_root),
    )
    matrix, _ = linear_eq_to_matrix(equations, PARAMETERS)
    return matrix


@cache
def exact_split_rank_boundary_certificate() -> SplitRankBoundaryCertificate:
    """Build the exact generic-rank, symmetry, and overlap certificate."""

    k0_v_c3_equations = _component_contact_equations(
        ZERO_FIBER_VERTICAL,
        R,
        FIRST_ROOT,
        3,
    )
    k0_v_c3_matrix, k0_v_c3_rhs = linear_eq_to_matrix(
        k0_v_c3_equations,
        PARAMETERS,
    )

    def c3_matrix(polynomial: Expr, variable: Expr) -> Matrix:
        equations = _component_contact_equations(
            polynomial,
            variable,
            FIRST_ROOT,
            3,
        )
        matrix, _ = linear_eq_to_matrix(equations, PARAMETERS)
        return matrix

    k0_vw_matrix = _contact_matrix(
        (ZERO_FIBER_VERTICAL, R, FIRST_ROOT),
        (ZERO_FIBER_GRAPH, S, SECOND_ROOT),
    )
    k0_vv_matrix = _contact_matrix(
        (ZERO_FIBER_VERTICAL, R, FIRST_ROOT),
        (ZERO_FIBER_VERTICAL, R, SECOND_ROOT),
    )
    k2_vv_matrix = _contact_matrix(
        (exceptional_vertical(1), R, FIRST_ROOT),
        (exceptional_vertical(1), R, SECOND_ROOT),
    )
    sign_flip = {BETA: -BETA, DELTA: -DELTA}
    overlap_form = expand(ZERO_FIBER_GRAPH.subs(S, 0))
    return SplitRankBoundaryCertificate(
        k0_vertical_c3_ranks=(
            k0_v_c3_matrix.rank(),
            k0_v_c3_matrix.row_join(k0_v_c3_rhs).rank(),
        ),
        k0_graph_c3_minor_gcd=_maximal_minor_gcd(c3_matrix(ZERO_FIBER_GRAPH, S)).subs(
            FIRST_ROOT, S
        ),
        k2_vertical_c3_minor_gcd=_maximal_minor_gcd(
            c3_matrix(exceptional_vertical(1), R)
        ),
        k2_graph_c3_minor_gcd=_maximal_minor_gcd(c3_matrix(exceptional_graph(1), S)),
        k0_vw_double_contact_determinant=expand(k0_vw_matrix.det()),
        k0_vv_double_contact_determinant=expand(k0_vv_matrix.det()),
        k2_vv_double_contact_determinant=expand(k2_vv_matrix.det()),
        plus_minus_vertical_symmetry=expand(
            exceptional_vertical(1).subs(sign_flip, simultaneous=True)
            - exceptional_vertical(-1)
        ),
        plus_minus_graph_symmetry=expand(
            exceptional_graph(1).subs(sign_flip, simultaneous=True).subs(S, -S)
            - exceptional_graph(-1)
        ),
        k0_overlap_value_identity=expand(
            overlap_form - 4 * ZERO_FIBER_VERTICAL.subs(R, Rational(1, 2))
        ),
        k0_overlap_vertical_derivative=expand(
            diff(ZERO_FIBER_VERTICAL, R).subs(R, Rational(1, 2))
        ),
        k0_overlap_graph_derivative=expand(diff(ZERO_FIBER_GRAPH, S).subs(S, 0)),
    )


@dataclass(frozen=True, slots=True)
class PrincipalRankSaturationCertificate:
    """Exact localization certificate for a principal rank-drop ideal.

    If ``h^e = q*g`` with ``g`` the rank-drop generator and ``h`` the boundary
    excluded from the clean chart, then ``(g) : h^infinity`` is the unit ideal.
    The nonzero division remainders for lower powers certify that ``e`` is the
    first such exponent, not merely an upper bound.
    """

    name: str
    rank_drop_generator: Expr
    excluded_boundary: Expr
    saturation_exponent: int
    quotient: Expr
    membership_identity: Expr
    lower_power_remainders: tuple[Expr, ...]

    @property
    def saturated_ideal_is_unit(self) -> bool:
        """Whether localization by the boundary removes the rank-drop locus."""

        return bool(self.membership_identity == 0)

    @property
    def verified(self) -> bool:
        """Whether the exact minimal-exponent saturation calculation agrees."""

        return bool(
            self.saturation_exponent > 0
            and len(self.lower_power_remainders) == self.saturation_exponent - 1
            and all(remainder != 0 for remainder in self.lower_power_remainders)
            and self.saturated_ideal_is_unit
        )


def _principal_rank_saturation(
    name: str,
    generator: Expr,
    excluded_boundary: Expr,
    variables: tuple[Symbol, ...],
) -> PrincipalRankSaturationCertificate:
    """Certify ``(generator) : excluded_boundary^infinity = (1)`` exactly."""

    generator_poly = Poly(generator, *variables)
    lower_remainders: list[Expr] = []
    for exponent in range(1, 9):
        quotient, remainder = Poly(excluded_boundary**exponent, *variables).div(
            generator_poly
        )
        if remainder.is_zero:
            quotient_expression = expand(quotient.as_expr())
            return PrincipalRankSaturationCertificate(
                name=name,
                rank_drop_generator=expand(generator),
                excluded_boundary=expand(excluded_boundary),
                saturation_exponent=exponent,
                quotient=quotient_expression,
                membership_identity=expand(
                    excluded_boundary**exponent - quotient_expression * generator
                ),
                lower_power_remainders=tuple(lower_remainders),
            )
        lower_remainders.append(expand(remainder.as_expr()))
    msg = f"no principal saturation exponent at most eight for {name}"
    raise ValueError(msg)


@cache
def exact_split_rank_saturation_certificates() -> tuple[
    PrincipalRankSaturationCertificate, ...
]:
    """Localize the three principal generic-rank obstructions exactly.

    These are precisely the charts in the finite ledger whose maximal-rank
    condition is controlled by one already excluded boundary factor.  The
    remaining mixed/graph double-contact rank divisors are deliberately not
    promoted here; their compatible residual-rank subloci remain part of the
    stated boundary gap.
    """

    rank = exact_split_rank_boundary_certificate()
    return (
        _principal_rank_saturation(
            "k0_graph_c3",
            rank.k0_graph_c3_minor_gcd,
            S**2 + 1,
            (S,),
        ),
        _principal_rank_saturation(
            "k0_vw_double_contact",
            rank.k0_vw_double_contact_determinant,
            SECOND_ROOT * (SECOND_ROOT**2 + 1),
            (SECOND_ROOT,),
        ),
        _principal_rank_saturation(
            "k2_vv_double_contact",
            rank.k2_vv_double_contact_determinant,
            FIRST_ROOT - SECOND_ROOT,
            (FIRST_ROOT, SECOND_ROOT),
        ),
    )


@dataclass(frozen=True, slots=True)
class SplitCodimensionTwoCertificate:
    """Aggregate finite ledger and exact witness certificate."""

    charts: tuple[SplitChartCertificate, ...]
    allocations: tuple[SplitAllocation, ...]
    forbidden_k0_allocations: tuple[tuple[str, str, str], ...]
    witnesses: tuple[SplitWitnessCertificate, ...]
    rank_boundaries: SplitRankBoundaryCertificate
    rank_saturations: tuple[PrincipalRankSaturationCertificate, ...]
    witness_keys_cover_allocations: bool
    k0_profile_allocation_counts: tuple[tuple[str, int], ...]
    kpm_profile_allocation_counts: tuple[tuple[str, int], ...]
    critical_triple_boundary_open: bool
    topology_not_computed: bool

    @property
    def verified(self) -> bool:
        """Whether every exact split-chart and finite-ledger check agrees."""

        specs = split_witness_specs()
        return bool(
            all(chart.verified for chart in self.charts)
            and all(allocation.verified for allocation in self.allocations)
            and self.forbidden_k0_allocations == FORBIDDEN_K0_ALLOCATIONS
            and len(self.witnesses) == len(specs) == len(self.allocations)
            and all(
                certificate.verified(spec)
                for certificate, spec in zip(self.witnesses, specs, strict=True)
            )
            and self.rank_boundaries.verified
            and len(self.rank_saturations) == 3
            and all(saturation.verified for saturation in self.rank_saturations)
            and self.witness_keys_cover_allocations
            and self.k0_profile_allocation_counts
            == (
                ("C2+T111+5N", 1),
                ("C2^2+6N", 3),
                ("C3+7N", 3),
                ("Q0+4N", 1),
                ("T111^2+4N", 1),
                ("T112+6N", 2),
            )
            and self.kpm_profile_allocation_counts
            == (
                ("C2+T111+5N", 2),
                ("C2^2+6N", 3),
                ("C3+7N", 2),
                ("Q0+4N", 1),
                ("T111^2+4N", 1),
                ("T112+6N", 2),
            )
            and self.critical_triple_boundary_open
            and self.topology_not_computed
        )


def _profile_counts(split_class: str) -> tuple[tuple[str, int], ...]:
    counts: dict[str, int] = {}
    for allocation in SPLIT_ALLOCATIONS:
        if allocation.split_class != split_class:
            continue
        counts[allocation.profile] = counts.get(allocation.profile, 0) + 1
    return tuple(sorted(counts.items()))


@cache
def exact_split_codimension_two_certificate() -> SplitCodimensionTwoCertificate:
    """Build the complete exact split-chart ledger."""

    witnesses = exact_split_witness_certificates()
    witness_keys = {
        (
            witness.profile,
            "k=0" if witness.kappa == 0 else "k=+/-2",
            witness.allocation,
        )
        for witness in witnesses
    }
    allocation_keys = {
        (allocation.profile, allocation.split_class, allocation.allocation)
        for allocation in SPLIT_ALLOCATIONS
    }
    return SplitCodimensionTwoCertificate(
        charts=exact_split_chart_certificates(),
        allocations=SPLIT_ALLOCATIONS,
        forbidden_k0_allocations=FORBIDDEN_K0_ALLOCATIONS,
        witnesses=witnesses,
        rank_boundaries=exact_split_rank_boundary_certificate(),
        rank_saturations=exact_split_rank_saturation_certificates(),
        witness_keys_cover_allocations=witness_keys == allocation_keys,
        k0_profile_allocation_counts=_profile_counts("k=0"),
        kpm_profile_allocation_counts=_profile_counts("k=+/-2"),
        # At k=+/-2 a P-critical fiber can contain three distinct source
        # points while the curve branch remains immersed if Q' is nonzero.
        # Those lower-dimensional triple/fourth-root loci are not classified
        # by the generic unramified allocation ledger above.
        critical_triple_boundary_open=True,
        # This module proves algebraic incidence and clean witnesses only.
        topology_not_computed=True,
    )


def main() -> int:
    """Print the exact split-chart classification and honest residual gaps."""

    certificate = exact_split_codimension_two_certificate()
    print("delta-ten split codimension-two certificate:", certificate.verified)
    print(
        "split component degrees:",
        tuple(
            (chart.kappa, chart.vertical_degree, chart.graph_degree)
            for chart in certificate.charts
        ),
    )
    print("k=0 allocation counts:", certificate.k0_profile_allocation_counts)
    print("k=+/-2 allocation counts:", certificate.kpm_profile_allocation_counts)
    print("exact clean allocation witnesses:", len(certificate.witnesses))
    print(
        "principal rank saturations:",
        tuple(
            (saturation.name, saturation.saturation_exponent)
            for saturation in certificate.rank_saturations
        ),
    )
    print("forbidden k=0 allocations:", certificate.forbidden_k0_allocations)
    print(
        "remaining: P-critical triple/fourth-root loci, compatible residual-rank "
        "subloci, component intersections beyond the k=0 overlap, and topology"
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
