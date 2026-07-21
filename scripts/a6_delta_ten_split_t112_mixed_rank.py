"""Close the split rank strata for ``T112`` and ``C2 + T111``.

At ``k=0`` and ``k=2`` the unordered-pair incidence of

``P=t^2+k*t^3+t^4``

is reducible.  The cancelled generic collision decic is therefore the wrong
object for a rank audit.  This module derives seven affine-linear incidence
systems directly from the true vertical and graph components:

* ``T112`` with the tangent edge on ``V`` or ``W`` at ``k=0,2``;
* ``C2+T111`` with the separate contact on ``W`` at ``k=0``; and
* ``C2+T111`` with the separate contact on ``V`` or ``W`` at ``k=2``.

A labeled triple fiber is written

``(z,-z,w,-w)`` with ``z^2+w^2+1=0`` at ``k=0``, and
``(z,-1-z,w,-1-w)`` with ``z(z+1)+w(w+1)=0`` at ``k=2``.

Triple equality is imposed by one value equation on each true component.
Tangency/contact is imposed by differentiating that component along its
one-dimensional pair chart.  The recorded localizers remove precisely the
repeated-source, ramified-pair, component-overlap, and same-target boundaries
needed for this derivative criterion to be equivalent to equality of branch
slopes.

The independent Sage checker proves:

* every valid ``T112`` matrix has rank exactly three, so its incidence is an
  affine-line bundle over a curve and has dimension two;
* every mixed matrix has generic rank four; its valid rank-three divisor is
  compatible over a reduced scheme of length four for ``k=0/W`` and
  ``k=2/V``, and nowhere for ``k=2/W``;
* coefficient rank at most two, and augmented rank at most two, saturate to
  the unit ideal in all three mixed charts.  The finite compatible points
  therefore have rank exactly three and affine-line coefficient fibers, so
  their total residual incidence has dimension one.

Thus none of these seven true-split allocation charts hides a
three-dimensional incidence component.  Exact hostile fixtures show why the
localizers cannot be omitted.  The full component equations, including their
derivatives, are transported from ``k=2`` to ``k=-2`` by
``t |-> -t`` and ``(a,b,c,d) |-> (a,-b,c,-d)``.

This is a conditional chart-level algebraic result.  It does not compute
complement groups, classify intersections with every other collision
profile, treat the other split allocations, or prove the plane Jacobian
conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final, Literal

from sympy import (
    Expr,
    I,
    Matrix,
    Poly,
    Rational,
    Symbol,
    diff,
    expand,
    linear_eq_to_matrix,
    rem,
    sqrt,
)

from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    DELTA,
    FAMILY_P,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    R,
    S,
    T,
)
from scripts.a6_delta_ten_split_codim_two import (
    COEFFICIENT_SIGN_FLIP,
    K0_FIBER_ONE,
    K2_FIBER_ONE,
    graph_target_x,
    split_graph,
    split_vertical,
    split_witness_specs,
    vertical_target_x,
)

Profile = Literal["T112+6N", "C2+T111+5N"]
Allocation = Literal["tangent-V", "tangent-W", "contact-V", "contact-W"]

COEFFICIENTS: Final = (ALPHA, BETA, GAMMA, DELTA)
ZERO_COEFFICIENTS: Final = dict.fromkeys(COEFFICIENTS, 0)
TRIPLE_ROOT: Final = Symbol("z_split_t112_mixed")
THIRD_ROOT: Final = Symbol("w_split_t112_mixed")
CONTACT_ROOT: Final = Symbol("h_split_t112_mixed")


@dataclass(frozen=True, slots=True)
class SplitTripleGeometry:
    """One labeled four-source ``P``-fiber on a true split chart."""

    kappa: int
    base_constraint: Expr
    sources: tuple[Expr, Expr, Expr, Expr]
    vertical_root: Expr
    graph_roots: tuple[Expr, Expr]
    valid_localizer: Expr

    @property
    def verified(self) -> bool:
        """Whether all four sources have the same first coordinate."""

        first_value = FAMILY_P.subs({KAPPA: self.kappa, T: self.sources[0]})
        fiber_residuals = tuple(
            _reduce_on_base(
                FAMILY_P.subs({KAPPA: self.kappa, T: source}) - first_value,
                self.base_constraint,
            )
            for source in self.sources[1:]
        )
        vertical_target_residual = _reduce_on_base(
            vertical_target_x(self.kappa).subs(R, self.vertical_root) - first_value,
            self.base_constraint,
        )
        graph_target_residuals = tuple(
            _reduce_on_base(
                graph_target_x(self.kappa).subs(S, root) - first_value,
                self.base_constraint,
            )
            for root in self.graph_roots
        )
        return bool(
            fiber_residuals == (0, 0, 0)
            and vertical_target_residual == 0
            and graph_target_residuals == (0, 0)
        )


def _reduce_on_base(expression: Expr, constraint: Expr) -> Expr:
    """Reduce exactly modulo the monic quadratic triple-base equation."""

    result = rem(
        Poly(expand(expression), THIRD_ROOT),
        Poly(constraint, THIRD_ROOT),
    )
    return expand(result.as_expr())


@cache
def split_triple_geometry(kappa: int) -> SplitTripleGeometry:
    """Return the exact labeled triple parameterization at ``k=0`` or ``2``."""

    z = TRIPLE_ROOT
    w = THIRD_ROOT
    if kappa == 0:
        sources = (z, -z, w, -w)
        constraint = z**2 + w**2 + 1
        vertical_root = -(z**2)
        graph_roots = (z + w, w - z)
        localizer = z * w * (z - w) * (z + w) * (2 * z**2 + 1)
    elif kappa == 2:
        sources = (z, -1 - z, w, -1 - w)
        constraint = z * (z + 1) + w * (w + 1)
        vertical_root = -z * (z + 1)
        graph_roots = (z + w, w - z - 1)
        localizer = (
            z
            * (z + 1)
            * w
            * (w + 1)
            * (2 * z + 1)
            * (2 * w + 1)
            * (z - w)
            * (z + w + 1)
        )
    else:
        msg = "the audited split triple charts are k=0 and k=2"
        raise ValueError(msg)
    return SplitTripleGeometry(
        kappa=kappa,
        base_constraint=expand(constraint),
        sources=sources,
        vertical_root=expand(vertical_root),
        graph_roots=tuple(expand(root) for root in graph_roots),
        valid_localizer=expand(localizer),
    )


@dataclass(frozen=True, slots=True)
class SplitIncidenceSpec:
    """One exact true-component incidence system and expected Sage result."""

    name: str
    profile: Profile
    kappa: int
    allocation: Allocation
    expected_rank: int
    base_dimension: int
    witness_name: str
    witness_contact_root: Expr
    rank_drop_saturation_exponent: int
    rank_drop_valid_dimension: int | None
    rank_drop_prime_count: int
    compatibility_saturation_exponent: int
    compatibility_valid_dimension: int | None
    compatibility_length: int
    compatibility_reduced: bool
    coefficient_rank_two_saturation_exponent: int | None
    augmented_rank_two_saturation_exponent: int | None

    @property
    def generic_coefficient_fiber_dimension(self) -> int:
        """Dimension of the affine solution fiber at maximal rank."""

        return 4 - self.expected_rank

    @property
    def generic_incidence_dimension(self) -> int:
        """Dimension of the maximal-rank incidence."""

        return self.base_dimension + self.generic_coefficient_fiber_dimension

    @property
    def residual_incidence_dimension_bound(self) -> int:
        """Largest incidence dimension supported on a rank-drop locus."""

        if self.profile == "T112+6N":
            return -1
        if self.compatibility_valid_dimension is None:
            return -1
        # Rank <= 2 is empty, so every compatible residual point has rank 3.
        return self.compatibility_valid_dimension + 1

    @property
    def hidden_threefold_excluded(self) -> bool:
        """Whether every generic and residual fiber has dimension below three."""

        return bool(
            self.generic_incidence_dimension == 2
            and self.residual_incidence_dimension_bound <= 1
        )


SPLIT_INCIDENCE_SPECS: Final = (
    SplitIncidenceSpec(
        "t112_k0_v",
        "T112+6N",
        0,
        "tangent-V",
        3,
        1,
        "t112_k0_v",
        0,
        1,
        None,
        0,
        1,
        None,
        0,
        False,
        None,
        None,
    ),
    SplitIncidenceSpec(
        "t112_k0_w",
        "T112+6N",
        0,
        "tangent-W",
        3,
        1,
        "t112_k0_w",
        0,
        1,
        None,
        0,
        1,
        None,
        0,
        False,
        None,
        None,
    ),
    SplitIncidenceSpec(
        "t112_k2_v",
        "T112+6N",
        2,
        "tangent-V",
        3,
        1,
        "t112_k2_v",
        0,
        2,
        None,
        0,
        2,
        None,
        0,
        False,
        None,
        None,
    ),
    SplitIncidenceSpec(
        "t112_k2_w",
        "T112+6N",
        2,
        "tangent-W",
        3,
        1,
        "t112_k2_w",
        0,
        2,
        None,
        0,
        2,
        None,
        0,
        False,
        None,
        None,
    ),
    SplitIncidenceSpec(
        "mixed_k0_w",
        "C2+T111+5N",
        0,
        "contact-W",
        4,
        2,
        "mixed_k0_w",
        -7 * I / 13,
        2,
        1,
        1,
        2,
        0,
        4,
        True,
        1,
        1,
    ),
    SplitIncidenceSpec(
        "mixed_k2_v",
        "C2+T111+5N",
        2,
        "contact-V",
        4,
        2,
        "mixed_k2_v",
        -6,
        2,
        1,
        1,
        2,
        0,
        4,
        True,
        1,
        1,
    ),
    SplitIncidenceSpec(
        "mixed_k2_w",
        "C2+T111+5N",
        2,
        "contact-W",
        4,
        2,
        "mixed_k2_w",
        -Rational(1, 13),
        2,
        1,
        2,
        2,
        None,
        0,
        False,
        1,
        1,
    ),
)


def _contact_localizer(spec: SplitIncidenceSpec) -> Expr:
    """Return the exact contact-chart and target-separation localizer."""

    if spec.profile == "T112+6N":
        return 1
    geometry = split_triple_geometry(spec.kappa)
    h = CONTACT_ROOT
    if spec.kappa == 0:
        s1, s2 = geometry.graph_roots
        return expand(h * (h**2 + 1) * (h**2 + 2) * (h**2 - s1**2) * (h**2 - s2**2))
    if spec.allocation == "contact-V":
        r = geometry.vertical_root
        return expand(h * (1 - 4 * h) * (h**2 - r**2))
    z = TRIPLE_ROOT
    w = THIRD_ROOT
    q1 = z + w + 1
    q2 = w - z
    return expand(
        h
        * (h + 1)
        * (h + 2)
        * (2 * (h + 1) ** 2 - 1)
        * ((h + 1) ** 2 - q1**2)
        * ((h + 1) ** 2 - q2**2)
    )


def valid_base_localizer(spec: SplitIncidenceSpec) -> Expr:
    """Return the complete base-only localizer used by the Sage audit."""

    return expand(
        split_triple_geometry(spec.kappa).valid_localizer * _contact_localizer(spec)
    )


def split_incidence_equations(spec: SplitIncidenceSpec) -> tuple[Expr, ...]:
    """Derive the affine-linear equations from the true ``V/W`` factors."""

    geometry = split_triple_geometry(spec.kappa)
    vertical = split_vertical(spec.kappa)
    graph = split_graph(spec.kappa)
    vertical_value = vertical.subs(R, geometry.vertical_root)
    graph_value = graph.subs(S, geometry.graph_roots[0])
    if spec.profile == "T112+6N":
        if spec.allocation == "tangent-V":
            tangent = diff(vertical, R).subs(R, geometry.vertical_root)
        else:
            tangent = diff(graph, S).subs(S, geometry.graph_roots[0])
        return tuple(
            expand(equation) for equation in (vertical_value, graph_value, tangent)
        )
    if spec.allocation == "contact-V":
        contact_value = vertical.subs(R, CONTACT_ROOT)
        contact_derivative = diff(vertical, R).subs(R, CONTACT_ROOT)
    else:
        contact_value = graph.subs(S, CONTACT_ROOT)
        contact_derivative = diff(graph, S).subs(S, CONTACT_ROOT)
    return tuple(
        expand(equation)
        for equation in (
            vertical_value,
            graph_value,
            contact_value,
            contact_derivative,
        )
    )


def split_incidence_matrices(spec: SplitIncidenceSpec) -> tuple[Matrix, Matrix]:
    """Return coefficient and augmented matrices for one split incidence."""

    matrix, right_hand_side = linear_eq_to_matrix(
        split_incidence_equations(spec),
        COEFFICIENTS,
    )
    return matrix, matrix.row_join(right_hand_side)


@dataclass(frozen=True, slots=True)
class ComponentDerivationCertificate:
    """Exact identities deriving triple equality from the true components."""

    kappa: int
    p_fiber_residuals: tuple[Expr, Expr, Expr]
    vertical_q_identity: Expr
    graph_q_identities: tuple[Expr, Expr]
    vertical_target_identity: Expr
    graph_target_identities: tuple[Expr, Expr]
    vertical_target_derivative: Expr
    graph_target_derivative: Expr
    target_derivative_localizer_identity: Expr

    @property
    def verified(self) -> bool:
        """Whether pair values and the derivative-validity factor agree."""

        return bool(
            self.p_fiber_residuals == (0, 0, 0)
            and self.vertical_q_identity == 0
            and self.graph_q_identities == (0, 0)
            and self.vertical_target_identity == 0
            and self.graph_target_identities == (0, 0)
            and self.vertical_target_derivative != 0
            and self.graph_target_derivative != 0
            and self.target_derivative_localizer_identity == 0
        )


@cache
def exact_component_derivation_certificates() -> tuple[
    ComponentDerivationCertificate, ...
]:
    """Prove that the displayed component rows are the actual pair equations."""

    certificates: list[ComponentDerivationCertificate] = []
    for kappa in (0, 2):
        geometry = split_triple_geometry(kappa)
        z, z_star, w, w_star = geometry.sources
        vertical = split_vertical(kappa)
        graph = split_graph(kappa)
        vertical_value = vertical.subs(R, geometry.vertical_root)
        graph_scale: Expr = Rational(1) if kappa == 0 else S**2
        vertical_scale = 2 * z**5 if kappa == 0 else 2 * z + 1
        vertical_q_identity = expand(
            FAMILY_Q.subs(T, z)
            - FAMILY_Q.subs(T, z_star)
            - vertical_scale * vertical_value
        )
        graph_q_identities: list[Expr] = []
        for left, right, graph_root in (
            (z, w, geometry.graph_roots[0]),
            (z_star, w, geometry.graph_roots[1]),
        ):
            identity = 16 * (FAMILY_Q.subs(T, left) - FAMILY_Q.subs(T, right)) - (
                left - right
            ) * graph_scale.subs(S, graph_root) * graph.subs(S, graph_root)
            graph_q_identities.append(
                _reduce_on_base(identity, geometry.base_constraint)
            )

        first_value = FAMILY_P.subs({KAPPA: kappa, T: z})
        p_residuals = tuple(
            _reduce_on_base(
                FAMILY_P.subs({KAPPA: kappa, T: source}) - first_value,
                geometry.base_constraint,
            )
            for source in (z_star, w, w_star)
        )
        vertical_target_identity = _reduce_on_base(
            vertical_target_x(kappa).subs(R, geometry.vertical_root) - first_value,
            geometry.base_constraint,
        )
        graph_target_identities = tuple(
            _reduce_on_base(
                graph_target_x(kappa).subs(S, root) - first_value,
                geometry.base_constraint,
            )
            for root in geometry.graph_roots
        )
        vertical_derivative = expand(
            diff(vertical_target_x(kappa), R).subs(R, geometry.vertical_root)
        )
        graph_derivative = _reduce_on_base(
            diff(graph_target_x(kappa), S).subs(S, geometry.graph_roots[0]),
            geometry.base_constraint,
        )
        if kappa == 0:
            derivative_factor = -2 * z * w * (z + w)
        else:
            derivative_factor = -((z + w + 1) * (2 * z + 1) * (2 * w + 1)) / 2
        certificates.append(
            ComponentDerivationCertificate(
                kappa=kappa,
                p_fiber_residuals=p_residuals,
                vertical_q_identity=vertical_q_identity,
                graph_q_identities=(graph_q_identities[0], graph_q_identities[1]),
                vertical_target_identity=vertical_target_identity,
                graph_target_identities=graph_target_identities,
                vertical_target_derivative=vertical_derivative,
                graph_target_derivative=graph_derivative,
                target_derivative_localizer_identity=_reduce_on_base(
                    graph_derivative - derivative_factor,
                    geometry.base_constraint,
                ),
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class TargetSeparationCertificate:
    """Exact factorization of the mixed same-target boundary."""

    name: str
    identity: Expr
    contact_derivative_factor: Expr
    contact_diagonal_factor: Expr

    @property
    def verified(self) -> bool:
        """Whether target separation and unramifiedness are nontrivial."""

        return bool(
            self.identity == 0
            and self.contact_derivative_factor != 0
            and self.contact_diagonal_factor != 0
        )


@cache
def exact_target_separation_certificates() -> tuple[TargetSeparationCertificate, ...]:
    """Factor the three contact-target separation conditions exactly."""

    results: list[TargetSeparationCertificate] = []
    for spec in SPLIT_INCIDENCE_SPECS:
        if spec.profile != "C2+T111+5N":
            continue
        geometry = split_triple_geometry(spec.kappa)
        triple_target = vertical_target_x(spec.kappa).subs(R, geometry.vertical_root)
        h = CONTACT_ROOT
        if spec.kappa == 0:
            s1, s2 = geometry.graph_roots
            contact_target = graph_target_x(0).subs(S, h)
            separation = (h**2 - s1**2) * (h**2 - s2**2)
            identity = _reduce_on_base(
                4 * (contact_target - triple_target) + separation,
                geometry.base_constraint,
            )
            derivative_factor = -h * (h**2 + 1)
            diagonal_factor = -(h**2 + 2)
        elif spec.allocation == "contact-V":
            contact_target = vertical_target_x(2).subs(R, h)
            separation = h**2 - geometry.vertical_root**2
            identity = expand(contact_target - triple_target - separation)
            derivative_factor = 2 * h
            diagonal_factor = 1 - 4 * h
        else:
            z = TRIPLE_ROOT
            w = THIRD_ROOT
            q1 = z + w + 1
            q2 = w - z
            contact_target = graph_target_x(2).subs(S, h)
            separation = ((h + 1) ** 2 - q1**2) * ((h + 1) ** 2 - q2**2)
            identity = _reduce_on_base(
                4 * (contact_target - triple_target) + separation,
                geometry.base_constraint,
            )
            derivative_factor = -((h + 1) * (2 * (h + 1) ** 2 - 1)) / 2
            diagonal_factor = -h * (h + 2)
        results.append(
            TargetSeparationCertificate(
                name=spec.name,
                identity=identity,
                contact_derivative_factor=expand(derivative_factor),
                contact_diagonal_factor=expand(diagonal_factor),
            )
        )
    return tuple(results)


MIXED_COMPATIBILITY_BASES: Final = {
    "mixed_k0_w": (
        TRIPLE_ROOT**2 + Rational(11, 12),
        CONTACT_ROOT**2 + Rational(1, 3),
        THIRD_ROOT - CONTACT_ROOT / 2,
    ),
    "mixed_k2_v": (
        TRIPLE_ROOT**2 + TRIPLE_ROOT - CONTACT_ROOT + Rational(1, 4),
        CONTACT_ROOT**2 - CONTACT_ROOT + Rational(1, 2),
        THIRD_ROOT - CONTACT_ROOT + Rational(3, 2),
    ),
    "mixed_k2_w": (1,),
}


@dataclass(frozen=True, slots=True)
class SplitIncidenceCertificate:
    """Python-side derivation, witness, and Sage result for one chart."""

    spec: SplitIncidenceSpec
    equation_count: int
    coefficient_rank: int
    augmented_rank: int
    witness_base_residual: Expr
    witness_localizer_value: Expr
    witness_equation_residuals: tuple[Expr, ...]
    mixed_compatibility_basis: tuple[Expr, ...]

    @property
    def verified(self) -> bool:
        """Whether the symbolic system and exact rank metadata agree."""

        expected_equations = 3 if self.spec.profile == "T112+6N" else 4
        expected_basis = MIXED_COMPATIBILITY_BASES.get(self.spec.name, ())
        return bool(
            self.equation_count == expected_equations
            and self.coefficient_rank == self.spec.expected_rank
            and self.augmented_rank == self.spec.expected_rank
            and self.witness_base_residual == 0
            and self.witness_localizer_value != 0
            and all(residual == 0 for residual in self.witness_equation_residuals)
            and self.mixed_compatibility_basis == expected_basis
            and self.spec.generic_incidence_dimension == 2
            and self.spec.hidden_threefold_excluded
            and (
                self.spec.profile == "T112+6N"
                or (
                    self.spec.coefficient_rank_two_saturation_exponent == 1
                    and self.spec.augmented_rank_two_saturation_exponent == 1
                )
            )
        )


def _witness_sources(kappa: int) -> tuple[Expr, Expr]:
    fiber = K0_FIBER_ONE if kappa == 0 else K2_FIBER_ONE
    return fiber[0], fiber[2]


@cache
def exact_split_incidence_certificates() -> tuple[SplitIncidenceCertificate, ...]:
    """Build all seven exact split incidence and rank certificates."""

    old_witnesses = {spec.name: spec for spec in split_witness_specs()}
    certificates: list[SplitIncidenceCertificate] = []
    for spec in SPLIT_INCIDENCE_SPECS:
        equations = split_incidence_equations(spec)
        matrix, augmented = split_incidence_matrices(spec)
        old_witness = old_witnesses[spec.witness_name]
        z_value, w_value = _witness_sources(spec.kappa)
        base_substitution = {
            TRIPLE_ROOT: z_value,
            THIRD_ROOT: w_value,
            CONTACT_ROOT: spec.witness_contact_root,
        }
        coefficient_substitution = dict(
            zip(COEFFICIENTS, old_witness.coefficient_values, strict=True)
        )
        substitution = base_substitution | coefficient_substitution
        geometry = split_triple_geometry(spec.kappa)
        certificates.append(
            SplitIncidenceCertificate(
                spec=spec,
                equation_count=len(equations),
                coefficient_rank=matrix.subs(base_substitution).rank(),
                augmented_rank=augmented.subs(base_substitution).rank(),
                witness_base_residual=expand(
                    geometry.base_constraint.subs(base_substitution)
                ),
                witness_localizer_value=expand(
                    valid_base_localizer(spec).subs(base_substitution)
                ),
                witness_equation_residuals=tuple(
                    expand(equation.subs(substitution)) for equation in equations
                ),
                mixed_compatibility_basis=MIXED_COMPATIBILITY_BASES.get(spec.name, ()),
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class PlusMinusFullSystemTransportCertificate:
    """Transport of every full ``k=2`` component incidence equation."""

    name: str
    base_identity: Expr
    p_identity: Expr
    q_identity: Expr
    equation_identities: tuple[Expr, ...]

    @property
    def verified(self) -> bool:
        """Whether the base, full family, and every incidence row transport."""

        return bool(
            self.base_identity == 0
            and self.p_identity == 0
            and self.q_identity == 0
            and all(identity == 0 for identity in self.equation_identities)
        )


@cache
def exact_plus_minus_full_system_transport_certificates() -> tuple[
    PlusMinusFullSystemTransportCertificate, ...
]:
    """Verify ``k=-2`` only after transporting every incidence equation."""

    certificates: list[PlusMinusFullSystemTransportCertificate] = []
    geometry = split_triple_geometry(2)
    z = TRIPLE_ROOT
    w = THIRD_ROOT
    minus_base = (-z) * (-z - 1) + (-w) * (-w - 1)
    minus_vertical = split_vertical(-2).subs(
        COEFFICIENT_SIGN_FLIP,
        simultaneous=True,
    )
    minus_graph = split_graph(-2).subs(
        COEFFICIENT_SIGN_FLIP,
        simultaneous=True,
    )
    for spec in SPLIT_INCIDENCE_SPECS:
        if spec.kappa != 2:
            continue
        plus_equations = split_incidence_equations(spec)
        minus_rows: list[Expr] = [
            minus_vertical.subs(R, geometry.vertical_root),
            minus_graph.subs(S, -geometry.graph_roots[0]),
        ]
        signs: list[int] = [1, 1]
        if spec.profile == "T112+6N":
            if spec.allocation == "tangent-V":
                minus_rows.append(
                    diff(minus_vertical, R).subs(R, geometry.vertical_root)
                )
                signs.append(1)
            else:
                minus_rows.append(
                    diff(minus_graph, S).subs(S, -geometry.graph_roots[0])
                )
                signs.append(-1)
        elif spec.allocation == "contact-V":
            minus_rows.extend(
                (
                    minus_vertical.subs(R, CONTACT_ROOT),
                    diff(minus_vertical, R).subs(R, CONTACT_ROOT),
                )
            )
            signs.extend((1, 1))
        else:
            minus_rows.extend(
                (
                    minus_graph.subs(S, -CONTACT_ROOT),
                    diff(minus_graph, S).subs(S, -CONTACT_ROOT),
                )
            )
            signs.extend((1, -1))
        certificates.append(
            PlusMinusFullSystemTransportCertificate(
                name=spec.name,
                base_identity=expand(minus_base - geometry.base_constraint),
                p_identity=expand(
                    FAMILY_P.subs({KAPPA: -2, T: -T}) - FAMILY_P.subs(KAPPA, 2)
                ),
                q_identity=expand(
                    FAMILY_Q.subs(COEFFICIENT_SIGN_FLIP, simultaneous=True).subs(T, -T)
                    + FAMILY_Q
                ),
                equation_identities=tuple(
                    expand(minus - sign * plus)
                    for minus, sign, plus in zip(
                        minus_rows,
                        signs,
                        plus_equations,
                        strict=True,
                    )
                ),
            )
        )
    return tuple(certificates)


@dataclass(frozen=True, slots=True)
class HostileSplitRankFixture:
    """A compatible rank-two point removed by the valid localizer."""

    name: str
    spec_name: str
    base_values: tuple[Expr, Expr, Expr]
    base_residual: Expr
    localizer_value: Expr
    coefficient_rank: int
    augmented_rank: int

    @property
    def verified(self) -> bool:
        """Whether the point is genuinely compatible but invalid."""

        return bool(
            self.base_residual == 0
            and self.localizer_value == 0
            and self.coefficient_rank == self.augmented_rank == 2
        )


HOSTILE_BASES: Final = {
    "t112_k0_v": (-I / sqrt(2), I / sqrt(2), 0),
    "t112_k0_w": (0, I, 0),
    "t112_k2_v": (0, -1, 0),
    "t112_k2_w": (0, -1, 0),
    "mixed_k0_w": (0, I, I),
    "mixed_k2_v": (0, -1, 0),
    "mixed_k2_w": (0, -1, -1),
}


@cache
def exact_hostile_split_rank_fixtures() -> tuple[HostileSplitRankFixture, ...]:
    """Replay one invalid compatible rank drop for every audited system."""

    fixtures: list[HostileSplitRankFixture] = []
    for spec in SPLIT_INCIDENCE_SPECS:
        z_value, w_value, h_value = HOSTILE_BASES[spec.name]
        substitution = {
            TRIPLE_ROOT: z_value,
            THIRD_ROOT: w_value,
            CONTACT_ROOT: h_value,
        }
        matrix, augmented = split_incidence_matrices(spec)
        geometry = split_triple_geometry(spec.kappa)
        fixtures.append(
            HostileSplitRankFixture(
                name=f"hostile_{spec.name}",
                spec_name=spec.name,
                base_values=(z_value, w_value, h_value),
                base_residual=expand(geometry.base_constraint.subs(substitution)),
                localizer_value=expand(valid_base_localizer(spec).subs(substitution)),
                coefficient_rank=matrix.subs(substitution).rank(),
                augmented_rank=augmented.subs(substitution).rank(),
            )
        )
    return tuple(fixtures)


@dataclass(frozen=True, slots=True)
class SplitT112MixedRankCertificate:
    """Aggregate exact true-split rank and compatibility certificate."""

    geometries: tuple[SplitTripleGeometry, ...]
    derivations: tuple[ComponentDerivationCertificate, ...]
    target_separations: tuple[TargetSeparationCertificate, ...]
    incidences: tuple[SplitIncidenceCertificate, ...]
    transports: tuple[PlusMinusFullSystemTransportCertificate, ...]
    hostile_fixtures: tuple[HostileSplitRankFixture, ...]
    audited_allocations: tuple[tuple[str, int, str], ...]
    unaudited_split_allocations_remain: bool
    topology_not_computed: bool

    @property
    def split_rank_strata_closed(self) -> bool:
        """Whether all seven named allocation charts have no hidden threefold."""

        return bool(
            len(self.geometries) == 2
            and all(geometry.verified for geometry in self.geometries)
            and len(self.derivations) == 2
            and all(derivation.verified for derivation in self.derivations)
            and len(self.target_separations) == 3
            and all(item.verified for item in self.target_separations)
            and len(self.incidences) == 7
            and all(incidence.verified for incidence in self.incidences)
            and all(
                incidence.spec.hidden_threefold_excluded
                for incidence in self.incidences
            )
            and len(self.transports) == 4
            and all(transport.verified for transport in self.transports)
            and len(self.hostile_fixtures) == 7
            and all(fixture.verified for fixture in self.hostile_fixtures)
            and len(self.audited_allocations) == 7
        )

    @property
    def verified(self) -> bool:
        """Compatibility alias for the conditional chart-level theorem."""

        return self.split_rank_strata_closed


@cache
def exact_split_t112_mixed_rank_certificate() -> SplitT112MixedRankCertificate:
    """Build the complete seven-chart split rank certificate."""

    return SplitT112MixedRankCertificate(
        geometries=tuple(split_triple_geometry(kappa) for kappa in (0, 2)),
        derivations=exact_component_derivation_certificates(),
        target_separations=exact_target_separation_certificates(),
        incidences=exact_split_incidence_certificates(),
        transports=exact_plus_minus_full_system_transport_certificates(),
        hostile_fixtures=exact_hostile_split_rank_fixtures(),
        audited_allocations=tuple(
            (spec.profile, spec.kappa, spec.allocation)
            for spec in SPLIT_INCIDENCE_SPECS
        ),
        # The other four codimension-two profiles and deeper intersections are
        # outside this deliberately focused certificate.
        unaudited_split_allocations_remain=True,
        topology_not_computed=True,
    )


def main() -> int:
    """Print the exact dimensions while preserving the remaining gaps."""

    certificate = exact_split_t112_mixed_rank_certificate()
    print("split T112/mixed rank strata closed:", certificate.verified)
    for incidence in certificate.incidences:
        spec = incidence.spec
        print(
            spec.name,
            "rank",
            spec.expected_rank,
            "generic incidence dimension",
            spec.generic_incidence_dimension,
            "compatible residual length",
            spec.compatibility_length,
            "residual incidence dimension bound",
            spec.residual_incidence_dimension_bound,
        )
    print("k=2 to k=-2 full-system transports:", len(certificate.transports))
    print("remaining: other split allocations, deeper intersections, and topology")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
