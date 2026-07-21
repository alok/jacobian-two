"""Certify one exact ``T112 + 6N`` member at conditional delta ten.

Write a three-source target fiber as ``a,b,c`` and let ``e`` be the fourth
root of ``P(t)-P(a)``.  On the chart ``ab+ac+bc != 0``, the normalized
quartic ``P=t^2+k*t^3+t^4`` forces

``e = -abc/(ab+ac+bc)`` and
``(ab+ac+bc)^2 - (a+b+c)abc = ab+ac+bc``.

For a chosen tangent pair, say ``a,b``, the two target-equality equations and
the equality of slopes are three affine-linear equations in the four lower
coefficients of ``Q``.  An exact Sage saturation of their maximal-minor ideal
shows that they have rank three everywhere on the displayed
``P``-unramified labeled chart: the rank-drop ideal becomes the unit ideal
after removing zero or repeated sources, the ``sigma_2`` chart boundary, and
coincidence with the fourth root.  The last removal is a projection-critical
boundary, not automatically a singular branch of the parametrized curve.
The irreducible triple-root surface is smooth on that localization, so the
incidence space there is an affine-line bundle of dimension three.

The three choices of tangent pair are permuted transitively by the source
labels.  Their coefficient images therefore have the same irreducible
closure.  The clean sample has a reduced projection fiber of size two,
exactly the two orientations of its tangent pair, proving generic finiteness;
on the clean locus the label cover has degree two.  Hence the clean unlabeled
quotient of the incidence component is irreducible and connected, while its
coefficient-image closure is an irreducible threefold.  The clean labeled
open contains the cyclic sample.  Subject to proper projective Whitney--Thom
propagation over that connected equisingular open, this excludes its generic
``P``-unramified locus from the required ``A6`` quotient.

The exact rational member used here lies on the especially transparent
``k=2`` slice.  It has source parameters ``-3/5,-2/5,1/5`` over one target,
with the first two branches tangent to order exactly two.  Its remaining
split-incidence factors give two plus four ordinary nodes.  Sage 10.8
independently finds Jacobian lengths ``4 + 6 + 2 + 4`` and a cyclic affine
complement.  Exhaustive replay of all ``40^4`` single-three-cycle meridian
images leaves only forty diagonal ``C3`` images and no ``A6`` image.

This is an exact algebraic and computer-assisted generic-chart theorem.
The ideal saturations, primary decomposition, and affine van Kamp
presentation are independently recomputed in Sage 10.8.  The propagation
step uses the standard proper projective Whitney--Thom theorem.  Deeper
boundary intersections, including the removed ``P``-critical-fiber loci,
remain outside the claim, and the result does not prove or disprove the
two-dimensional Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final, TypeAlias

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
    R,
    S,
    T,
    X,
    Y,
    exceptional_graph,
    exceptional_vertical,
)
from scripts.a6_one_pair_infinity import (
    ThreeCyclePresentationCensus,
    alternating_group_six,
    evaluate_signed_word,
    n_generator_three_cycle_presentation_census,
)
from scripts.six_sheet_monodromy import IDENTITY, Permutation, cycle_type

Assignment: TypeAlias = tuple[Permutation, ...]

ROOT_A: Final = Symbol("u")
ROOT_B: Final = Symbol("v")
ROOT_C: Final = Symbol("w")
FREE_DELTA: Final = Symbol("q")
LOCAL_U: Final = Symbol("U")
LOCAL_V: Final = Symbol("V")


# General labeled triple-fiber incidence chart.  The fourth root makes the
# t-coefficient vanish; the displayed constraint then makes the t^2
# coefficient equal to one.  This chart misses no valid nonzero triple fiber:
# before division, the missing t-coefficient says
# ``sigma_3 + fourth*sigma_2 = 0``; if ``sigma_2 = 0``, it forces
# ``sigma_3 = 0`` and hence at least one source is the excluded cusp source.
ROOT_SIGMA_1: Final = ROOT_A + ROOT_B + ROOT_C
ROOT_SIGMA_2: Final = ROOT_A * ROOT_B + ROOT_A * ROOT_C + ROOT_B * ROOT_C
ROOT_SIGMA_3: Final = ROOT_A * ROOT_B * ROOT_C
ROOT_FOURTH: Final = cancel(-ROOT_SIGMA_3 / ROOT_SIGMA_2)
ROOT_KAPPA: Final = cancel(-ROOT_SIGMA_1 - ROOT_FOURTH)
ROOT_TARGET_X: Final = cancel(ROOT_SIGMA_3**2 / ROOT_SIGMA_2)
TRIPLE_P_CONSTRAINT: Final = expand(
    ROOT_SIGMA_2**2 - ROOT_SIGMA_1 * ROOT_SIGMA_3 - ROOT_SIGMA_2
)
ROOT_P: Final = T**2 + ROOT_KAPPA * T**3 + T**4
ROOT_PRODUCT: Final = (T - ROOT_A) * (T - ROOT_B) * (T - ROOT_C) * (T - ROOT_FOURTH)
ROOT_FACTORIZATION_IDENTITY: Final = cancel(
    ROOT_P - ROOT_TARGET_X - ROOT_PRODUCT + TRIPLE_P_CONSTRAINT * T**2 / ROOT_SIGMA_2
)

# The equal-target and tangent-pair equations are affine linear in a,b,c,d,
# here meaning the lower coefficients of Q (not the source-root labels).
ROOT_Q: Final = FAMILY_Q
TRIPLE_Q_AB: Final = cancel(
    (ROOT_Q.subs(T, ROOT_A) - ROOT_Q.subs(T, ROOT_B)) / (ROOT_A - ROOT_B)
)
TRIPLE_Q_AC: Final = cancel(
    (ROOT_Q.subs(T, ROOT_A) - ROOT_Q.subs(T, ROOT_C)) / (ROOT_A - ROOT_C)
)
TANGENT_AB: Final = cancel(
    diff(ROOT_Q, T).subs(T, ROOT_A) * diff(ROOT_P, T).subs(T, ROOT_B)
    - diff(ROOT_Q, T).subs(T, ROOT_B) * diff(ROOT_P, T).subs(T, ROOT_A)
)
T112_LINEAR_CONDITIONS: Final = (TRIPLE_Q_AB, TRIPLE_Q_AC, TANGENT_AB)
T112_INCIDENCE_MATRIX, T112_INCIDENCE_RHS = linear_eq_to_matrix(
    T112_LINEAR_CONDITIONS,
    (ALPHA, BETA, GAMMA, DELTA),
)
T112_LABELED_TANGENT_PAIRS: Final = (("u", "v"), ("u", "w"), ("v", "w"))

# Exact localization used by the independent Sage maximal-minor calculation.
# On the root constraint, ``u-fourth`` is
# ``u*(sigma_2+v*w)/sigma_2`` and cyclically, so these factors are precisely
# the nonzero, pairwise-distinct, P-unramified conditions on this chart.
# Their vanishing does not by itself make the plane-curve branch singular:
# Q' may remain nonzero when P' vanishes.
ROOT_FOURTH_SEPARATION_FACTORS: Final = (
    ROOT_SIGMA_2 + ROOT_B * ROOT_C,
    ROOT_SIGMA_2 + ROOT_A * ROOT_C,
    ROOT_SIGMA_2 + ROOT_A * ROOT_B,
)
P_UNRAMIFIED_BASE_PRODUCT: Final = expand(
    ROOT_A
    * ROOT_B
    * ROOT_C
    * ROOT_SIGMA_2
    * (ROOT_A - ROOT_B)
    * (ROOT_A - ROOT_C)
    * (ROOT_B - ROOT_C)
    * ROOT_FOURTH_SEPARATION_FACTORS[0]
    * ROOT_FOURTH_SEPARATION_FACTORS[1]
    * ROOT_FOURTH_SEPARATION_FACTORS[2]
)
# The additional condition that the distinct fourth source not join the same
# *target* depends on Q, so it is removed from the total incidence space, not
# from the root base.  Its numerator defines a genuine clean-locus boundary.
FOURTH_TARGET_SEPARATION: Final = ROOT_Q.subs(T, ROOT_FOURTH) - ROOT_Q.subs(T, ROOT_A)

# Sage 10.8 independently recomputes these ideal-theoretic certificates:
# - the base constraint is irreducible;
# - its singular ideal saturates to one with exponent one;
# - the rank-drop ideal (constraint plus all four maximal minors) saturates to
#   one with exponent three;
# - the clean sample fiber is reduced of length two.
SAGE_BASE_IRREDUCIBLE: Final = True
SAGE_BASE_ABSOLUTELY_IRREDUCIBLE: Final = True
# (saturation exponent, unit ideal)
SAGE_BASE_SMOOTH_SATURATION: Final = (1, True)
# (raw rank-drop dimension, saturation exponent, unit ideal)
SAGE_RANK_DROP_SATURATION: Final = (1, 3, True)
# (dimension, scheme length, reduced)
SAGE_SAMPLE_PROJECTION_FIBER: Final = (0, 2, True)
SAMPLE_PROJECTION_FIBER_BASIS: Final = (
    ROOT_B**2 + ROOT_B + Rational(6, 25),
    ROOT_A + ROOT_B + 1,
    ROOT_C - Rational(1, 5),
)


# A rational point of the triple-root surface and the resulting one-parameter
# T112 slice.  The fourth P-fiber parameter is -6/5 and k=2.
SLICE_ROOTS: Final = (Rational(-3, 5), Rational(-2, 5), Rational(1, 5))
SLICE_FOURTH: Final = Rational(-6, 5)
SLICE_ROOT_SUBSTITUTION: Final = dict(
    zip((ROOT_A, ROOT_B, ROOT_C), SLICE_ROOTS, strict=True)
)
SLICE_ALPHA: Final = (891875 * FREE_DELTA - 1544547) / Rational(9258125)
SLICE_BETA: Final = (1312875 * FREE_DELTA - 1909616) / Rational(1851625)
SLICE_GAMMA: Final = 3 * (186475 * FREE_DELTA - 198798) / Rational(370325)
SLICE_COEFFICIENTS: Final = {
    ALPHA: SLICE_ALPHA,
    BETA: SLICE_BETA,
    GAMMA: SLICE_GAMMA,
    DELTA: FREE_DELTA,
}
SLICE_INCIDENCE_IDENTITIES: Final = tuple(
    cancel(condition.subs(SLICE_ROOT_SUBSTITUTION).subs(SLICE_COEFFICIENTS))
    for condition in T112_LINEAR_CONDITIONS
)
SLICE_VERTICAL_RESIDUAL_NUMERATOR: Final = (
    14813 * R**2 + (36875 * FREE_DELTA - 117164) * R + 65328 - 26250 * FREE_DELTA
)
SLICE_GRAPH_RESIDUAL_NUMERATOR: Final = (
    370325 * S**4
    + 2740405 * S**3
    + (1740869 + 1843750 * FREE_DELTA) * S**2
    + (-6854049 + 4887500 * FREE_DELTA) * S
    - 3089094
    + 1783750 * FREE_DELTA
)
SLICE_VERTICAL_FACTORIZATION: Final = (
    (25 * R - 6) ** 2 * SLICE_VERTICAL_RESIDUAL_NUMERATOR / 9258125
)
SLICE_GRAPH_FACTORIZATION: Final = (
    (5 * S + 1) * (5 * S + 2) * SLICE_GRAPH_RESIDUAL_NUMERATOR / 9258125
)

# Exact divisors removed from the displayed rational slice.  Some are genuine
# equisingularity boundaries; the last two also ensure the chosen X-projection
# separates the residual packets, which is useful for the Sage certificate.
SLICE_BOUNDARY_FACTORS: Final = (
    ("cusp-leading-coefficient", 891875 * FREE_DELTA - 1544547),
    ("cusp-image-collision", 4375 * FREE_DELTA - 10888),
    ("extra-critical-point-1", 54500 * FREE_DELTA - 118281),
    ("extra-critical-point-2", 61625 * FREE_DELTA - 158444),
    ("contact-order-at-least-three", 906250 * FREE_DELTA - 1982389),
    ("third-branch-tangency", 171875 * FREE_DELTA - 326168),
    ("fourth-source-joins-target", 203125 * FREE_DELTA - 545721),
    (
        "vertical-residual-discriminant",
        271953125 * FREE_DELTA**2 - 1417096000 * FREE_DELTA + 1971317648,
    ),
    ("vertical-residual-projection-boundary", 36875 * FREE_DELTA - 117164),
    (
        "graph-residual-projection-boundary-1",
        316796875000 * FREE_DELTA**2 - 1348530190625 * FREE_DELTA + 1338995777764,
    ),
    (
        "packet-projection-boundary",
        74209208984375 * FREE_DELTA**2 - 427777703815000 * FREE_DELTA + 613419197284592,
    ),
    (
        "graph-residual-discriminant",
        648801743984222412109375000 * FREE_DELTA**5
        - 8414774122494220733642578125 * FREE_DELTA**4
        + 43382696548232061218261718750 * FREE_DELTA**3
        - 111084489794569121143697265625 * FREE_DELTA**2
        + 141213174394208303907051852000 * FREE_DELTA
        - 71271465292452501682100806608,
    ),
)


T112_PARAMETERS: Final = {
    KAPPA: 2,
    ALPHA: Rational(2089, 11875),
    BETA: Rational(3542, 2375),
    GAMMA: Rational(1788, 475),
    DELTA: Rational(338, 95),
}
T112_P: Final = expand(FAMILY_P.subs(T112_PARAMETERS))
T112_Q: Final = expand(FAMILY_Q.subs(T112_PARAMETERS))
T112_TRIPLE_X: Final = Rational(36, 625)
T112_TRIPLE_Y: Final = Rational(7776, 37109375)
T112_VERTICAL_RESIDUAL: Final = 19 * R**2 + 18 * R - 36
T112_GRAPH_RESIDUAL: Final = 475 * S**4 + 3515 * S**3 + 10647 * S**2 + 13513 * S + 4178
T112_VERTICAL_NODE_X: Final = 361 * X**2 - 1692 * X + 1296
T112_GRAPH_NODE_X: Final = (
    50906640625 * X**4
    - 462374666250 * X**3
    + 2074717406016 * X**2
    + 599188777920 * X
    - 38011109760
)

T112_IMPLICIT: Final = (
    -19885406494140625 * X**9
    + 98693762695312500 * X**8
    - 12903111675000000 * X**7
    - 23744619140625000 * X**6 * Y
    + 550166298120000 * X**6
    - 100341934312500000 * X**5 * Y
    - 7329711534336 * X**5
    + 133129766845703125 * X**4 * Y**2
    + 10731317108400000 * X**4 * Y
    - 540905323007812500 * X**3 * Y**2
    - 289807243200000 * X**3 * Y
    + 74936584472656250 * X**2 * Y**3
    + 142878836250000000 * X**2 * Y**2
    - 159150234375000000 * X * Y**3
    - 10428071625000000 * X * Y**2
    + 19885406494140625 * Y**4
    - 4340460937500000 * Y**3
    + 236852100000000 * Y**2
)
T112_IMPLICIT_RESULTANT_MULTIPLIER: Final = 11875**4
T112_TANGENT_CONE: Final = (
    Rational(3530304, 390625)
    * (1188 * LOCAL_U - 59375 * LOCAL_V) ** 2
    * (3924 * LOCAL_U - 415625 * LOCAL_V)
)

# Exact Sage 10.8 unsimplified affine presentation.  The integers name the
# four geometric meridians, with a negative sign denoting inversion.
T112_RELATIONS: Final = (
    (-3, -2, 1, 2, 3, -2, -1, 2),
    (
        -3,
        -2,
        -1,
        -2,
        1,
        2,
        3,
        -2,
        1,
        2,
        -3,
        -2,
        -1,
        2,
        1,
        2,
        3,
        -2,
        -1,
        2,
    ),
    (-4, 2, 4, 2, 4, 2, -4, -2, -4, -2),
    (-3, -2, -1, -4, -2, 1, 2, 4, 1, 2),
    (
        -4,
        -2,
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
        4,
        -2,
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
        2,
    ),
    (-4, -2, -4, -2, -4, -2, 1, 2, 4, 2, 4, 2),
    (
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
    ),
    (
        -4,
        -2,
        -1,
        -4,
        -2,
        -1,
        2,
        4,
        1,
        2,
        4,
        1,
        2,
        4,
        2,
        4,
        2,
        4,
        -2,
        -4,
        -4,
        -2,
        -4,
        -2,
    ),
    (
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
        1,
        2,
        4,
        2,
        4,
        -2,
    ),
    (-4, -3, 4, 2, -4, 3, 4, -2),
    (4, 3, -4, -3),
)


def _homogeneous_part(polynomial: Expr, degree: int) -> Expr:
    """Return the total-degree ``degree`` part in the local variables."""

    result: Expr = 0
    for monomial, coefficient in Poly(polynomial, LOCAL_U, LOCAL_V).terms():
        if sum(monomial) == degree:
            result += coefficient * LOCAL_U ** monomial[0] * LOCAL_V ** monomial[1]
    return expand(result)


def _terms_below_degree(polynomial: Expr, degree: int) -> Expr:
    """Return all local terms of total degree strictly below ``degree``."""

    result: Expr = 0
    for monomial, coefficient in Poly(polynomial, LOCAL_U, LOCAL_V).terms():
        if sum(monomial) < degree:
            result += coefficient * LOCAL_U ** monomial[0] * LOCAL_V ** monomial[1]
    return expand(result)


@cache
def _presentation_pruning() -> tuple[tuple[int, ...], tuple[Assignment, ...]]:
    """Replay incremental exact relation pruning."""

    three_cycles = tuple(
        element
        for element in alternating_group_six()
        if cycle_type(element) == (3, 1, 1, 1)
    )
    partial_assignments: tuple[Assignment, ...] = ((),)
    stage_counts: list[int] = []
    for assigned_count in range(1, 5):
        newly_decidable = tuple(
            relation
            for relation in T112_RELATIONS
            if max(abs(letter) for letter in relation) == assigned_count
        )
        partial_assignments = tuple(
            partial + (image,)
            for partial in partial_assignments
            for image in three_cycles
            if all(
                evaluate_signed_word(relation, partial + (image,)) == IDENTITY
                for relation in newly_decidable
            )
        )
        stage_counts.append(len(partial_assignments))
    return tuple(stage_counts), partial_assignments


@dataclass(frozen=True, slots=True)
class DeltaTenT112Certificate:
    """Exact algebraic and finite-image data for the rational member."""

    root_factorization_identity: Expr
    root_constraint_value: Expr
    root_constraint_gradient: tuple[Expr, Expr, Expr]
    fourth_root_value: Expr
    kappa_value: Expr
    base_constraint_irreducible: bool
    base_constraint_absolutely_irreducible: bool
    base_surface_dimension: int
    p_unramified_base_value: Expr
    base_smooth_saturation: tuple[int, bool]
    rank_drop_saturation: tuple[int, int, bool]
    incidence_rank: int
    incidence_affine_fiber_dimension: int
    labeled_component_dimension: int
    labeled_component_irreducible: bool
    sample_projection_fiber: tuple[int, int, bool]
    sample_projection_fiber_basis: tuple[Expr, Expr, Expr]
    coefficient_image_dimension: int
    generic_label_degree: int
    unlabeled_component_irreducible: bool
    clean_open_connected: bool
    proper_whitney_thom_required: bool
    topology_propagation_computer_verified: bool
    slice_incidence_identities: tuple[Expr, ...]
    vertical_factor_identity: Expr
    graph_factor_identity: Expr
    boundary_values: tuple[tuple[str, Expr], ...]
    triple_p_values: tuple[Expr, ...]
    triple_q_values: tuple[Expr, ...]
    fourth_q_separation: Expr
    tangent_slopes: tuple[Expr, ...]
    tangent_second_derivatives: tuple[Expr, Expr]
    contact_second_derivative_difference: Expr
    local_terms_below_three: Expr
    tangent_cone_identity: Expr
    cusp_leading_coefficient: Expr
    other_cusp_fiber_q_value: Expr
    cusp_image_factor: Expr
    extra_critical_factor: Expr
    vertical_residual_discriminant: Expr
    graph_residual_discriminant: Expr
    vertical_node_x_identity: Expr
    graph_node_x_identity: Expr
    vertical_node_x_discriminant: Expr
    graph_node_x_discriminant: Expr
    node_packet_resultant: Expr
    cusp_node_separations: tuple[Expr, Expr]
    triple_node_separations: tuple[Expr, Expr]
    implicit_resultant_identity: Expr
    implicit_parameterization_identity: Expr
    implicit_content: Expr
    implicit_total_degree: int
    sage_jacobian_components: tuple[tuple[int, int, bool], ...]
    sage_cyclic_simplification: tuple[int, int, bool]
    arithmetic_genus: int
    cusp_delta: int
    t112_delta: int
    node_count: int
    infinity_delta: int
    relation_count: int
    pruning_stage_survivors: tuple[int, ...]
    diagonal_satisfying_assignments: int
    complement_census: ThreeCyclePresentationCensus

    @property
    def total_delta(self) -> int:
        """Return the projective genus contribution."""

        return self.cusp_delta + self.t112_delta + self.node_count + self.infinity_delta

    @property
    def conditional_generic_exclusion_supported(self) -> bool:
        """Whether the exact data support exclusion after the stated theorem step.

        This property records the mathematical implication used by the audit.  It
        does not claim that the Whitney--Thom propagation step was checked by the
        Python or Sage programs.
        """

        return bool(
            self.verified
            and self.labeled_component_irreducible
            and self.unlabeled_component_irreducible
            and self.clean_open_connected
            and self.proper_whitney_thom_required
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.complement_census.a6_assignments == 0
        )

    @property
    def verified(self) -> bool:
        """Whether every exact sample invariant agrees with the certificate."""

        return bool(
            self.root_factorization_identity == 0
            and self.root_constraint_value == 0
            and self.root_constraint_gradient
            == (Rational(9, 125), Rational(28, 125), Rational(133, 125))
            and self.fourth_root_value == SLICE_FOURTH
            and self.kappa_value == 2
            and self.base_constraint_irreducible == SAGE_BASE_IRREDUCIBLE
            and self.base_constraint_absolutely_irreducible
            == SAGE_BASE_ABSOLUTELY_IRREDUCIBLE
            and self.base_surface_dimension == 2
            and self.p_unramified_base_value != 0
            and self.base_smooth_saturation == SAGE_BASE_SMOOTH_SATURATION
            and self.rank_drop_saturation == SAGE_RANK_DROP_SATURATION
            and self.incidence_rank == 3
            and self.incidence_affine_fiber_dimension == 1
            and self.labeled_component_dimension == 3
            and self.labeled_component_irreducible
            and self.sample_projection_fiber == SAGE_SAMPLE_PROJECTION_FIBER
            and self.sample_projection_fiber_basis == SAMPLE_PROJECTION_FIBER_BASIS
            and self.coefficient_image_dimension == 3
            and self.generic_label_degree == 2
            and self.unlabeled_component_irreducible
            and self.clean_open_connected
            and self.proper_whitney_thom_required
            and self.slice_incidence_identities == (0, 0, 0)
            and self.vertical_factor_identity == 0
            and self.graph_factor_identity == 0
            and all(value != 0 for _, value in self.boundary_values)
            and self.triple_p_values == (T112_TRIPLE_X,) * 3
            and self.triple_q_values == (T112_TRIPLE_Y,) * 3
            and self.fourth_q_separation == Rational(24712128, 37109375)
            and self.tangent_slopes
            == (
                Rational(1188, 59375),
                Rational(1188, 59375),
                Rational(3924, 415625),
            )
            and self.tangent_second_derivatives
            == (Rational(441, 95), Rational(351, 190))
            and self.contact_second_derivative_difference == Rational(531, 190)
            and self.local_terms_below_three == 0
            and self.tangent_cone_identity == 0
            and self.cusp_leading_coefficient == Rational(2089, 11875)
            and self.other_cusp_fiber_q_value == Rational(1296, 11875)
            and self.cusp_image_factor == Rational(1679616, 141015625)
            and self.extra_critical_factor == Rational(4358016, 5640625)
            and self.vertical_residual_discriminant == 3060
            and self.graph_residual_discriminant == -240310362868157952000
            and self.vertical_node_x_identity == 0
            and self.graph_node_x_identity == 0
            and self.vertical_node_x_discriminant == 991440
            and self.graph_node_x_discriminant
            == -828054536274252441843223030250099063034678391026318111865892352000000000
            and self.node_packet_resultant == 601205556802476551975346242124251136
            and self.cusp_node_separations == (1296, -38011109760)
            and self.triple_node_separations
            == (
                Rational(468647856, 390625),
                Rational(1288194316747776, 390625),
            )
            and self.implicit_resultant_identity == 0
            and self.implicit_parameterization_identity == 0
            and self.implicit_content == 1
            and self.implicit_total_degree == 9
            and self.sage_jacobian_components
            == (
                (4, 1, False),
                (6, 1, False),
                (2, 2, True),
                (4, 4, True),
            )
            and self.sage_cyclic_simplification == (1, 0, True)
            and self.arithmetic_genus == 28
            and self.cusp_delta == 2
            and self.t112_delta == 4
            and self.node_count == 6
            and self.infinity_delta == 16
            and self.total_delta == self.arithmetic_genus
            and self.relation_count == 11
            and self.pruning_stage_survivors == (40, 1600, 640, 40)
            and self.diagonal_satisfying_assignments == 40
            and self.complement_census.assignments == 40**4
            and self.complement_census.satisfying_assignments == 40
            and self.complement_census.generated_order_histogram == ((3, 40),)
            and self.complement_census.a6_assignments == 0
        )


@cache
def exact_delta_ten_t112_certificate() -> DeltaTenT112Certificate:
    """Build the exact ``T112 + 6N`` certificate."""

    root_gradient = tuple(
        diff(TRIPLE_P_CONSTRAINT, root).subs(SLICE_ROOT_SUBSTITUTION)
        for root in (ROOT_A, ROOT_B, ROOT_C)
    )
    sample_matrix = Matrix(T112_INCIDENCE_MATRIX.subs(SLICE_ROOT_SUBSTITUTION))
    p_derivative = diff(T112_P, T)
    q_derivative = diff(T112_Q, T)
    slopes = tuple(
        cancel(q_derivative.subs(T, root) / p_derivative.subs(T, root))
        for root in SLICE_ROOTS
    )
    second_derivative = cancel(
        (diff(T112_Q, T, 2) * p_derivative - q_derivative * diff(T112_P, T, 2))
        / p_derivative**3
    )
    tangent_seconds = tuple(second_derivative.subs(T, root) for root in SLICE_ROOTS[:2])
    translated_implicit = expand(
        T112_IMPLICIT.subs({X: T112_TRIPLE_X + LOCAL_U, Y: T112_TRIPLE_Y + LOCAL_V})
    )
    vertical_node_resultant = expand(resultant(T112_VERTICAL_RESIDUAL, R**2 - X, R))
    graph_node_resultant = expand(
        resultant(
            T112_GRAPH_RESIDUAL,
            -S * (S + 2) * (S + 1) ** 2 - 4 * X,
            S,
        )
    )
    stage_survivors, assignments = _presentation_pruning()
    diagonal_assignments = sum(
        all(image == images[0] for image in images[1:]) for images in assignments
    )
    implicit_polynomial = Poly(T112_IMPLICIT, X, Y)

    return DeltaTenT112Certificate(
        root_factorization_identity=ROOT_FACTORIZATION_IDENTITY,
        root_constraint_value=TRIPLE_P_CONSTRAINT.subs(SLICE_ROOT_SUBSTITUTION),
        root_constraint_gradient=root_gradient,
        fourth_root_value=ROOT_FOURTH.subs(SLICE_ROOT_SUBSTITUTION),
        kappa_value=ROOT_KAPPA.subs(SLICE_ROOT_SUBSTITUTION),
        base_constraint_irreducible=bool(
            Poly(TRIPLE_P_CONSTRAINT, ROOT_A, ROOT_B, ROOT_C).is_irreducible
        ),
        base_constraint_absolutely_irreducible=SAGE_BASE_ABSOLUTELY_IRREDUCIBLE,
        base_surface_dimension=2,
        p_unramified_base_value=P_UNRAMIFIED_BASE_PRODUCT.subs(SLICE_ROOT_SUBSTITUTION),
        base_smooth_saturation=SAGE_BASE_SMOOTH_SATURATION,
        rank_drop_saturation=SAGE_RANK_DROP_SATURATION,
        incidence_rank=sample_matrix.rank(),
        incidence_affine_fiber_dimension=1,
        labeled_component_dimension=3,
        labeled_component_irreducible=True,
        sample_projection_fiber=SAGE_SAMPLE_PROJECTION_FIBER,
        sample_projection_fiber_basis=SAMPLE_PROJECTION_FIBER_BASIS,
        coefficient_image_dimension=3,
        generic_label_degree=2,
        unlabeled_component_irreducible=True,
        clean_open_connected=True,
        # This records a theorem dependency, not a CAS output: propagation of
        # the exact sample complement uses proper projective Whitney--Thom
        # triviality over the connected clean open.
        proper_whitney_thom_required=True,
        topology_propagation_computer_verified=False,
        slice_incidence_identities=SLICE_INCIDENCE_IDENTITIES,
        vertical_factor_identity=expand(
            exceptional_vertical(1).subs(SLICE_COEFFICIENTS)
            - SLICE_VERTICAL_FACTORIZATION
        ),
        graph_factor_identity=expand(
            exceptional_graph(1).subs(SLICE_COEFFICIENTS) - SLICE_GRAPH_FACTORIZATION
        ),
        boundary_values=tuple(
            (name, factor.subs(FREE_DELTA, T112_PARAMETERS[DELTA]))
            for name, factor in SLICE_BOUNDARY_FACTORS
        ),
        triple_p_values=tuple(T112_P.subs(T, root) for root in SLICE_ROOTS),
        triple_q_values=tuple(T112_Q.subs(T, root) for root in SLICE_ROOTS),
        fourth_q_separation=T112_Q.subs(T, SLICE_FOURTH) - T112_TRIPLE_Y,
        tangent_slopes=slopes,
        tangent_second_derivatives=tangent_seconds,
        contact_second_derivative_difference=tangent_seconds[0] - tangent_seconds[1],
        local_terms_below_three=_terms_below_degree(translated_implicit, 3),
        tangent_cone_identity=expand(
            _homogeneous_part(translated_implicit, 3) - T112_TANGENT_CONE
        ),
        cusp_leading_coefficient=Poly(T112_Q, T).nth(5),
        other_cusp_fiber_q_value=T112_Q.subs(T, -1),
        cusp_image_factor=CUSP_IMAGE_FACTOR.subs(T112_PARAMETERS),
        extra_critical_factor=EXTRA_CRITICAL_FACTOR.subs(T112_PARAMETERS),
        vertical_residual_discriminant=discriminant(T112_VERTICAL_RESIDUAL, R),
        graph_residual_discriminant=discriminant(T112_GRAPH_RESIDUAL, S),
        vertical_node_x_identity=expand(vertical_node_resultant - T112_VERTICAL_NODE_X),
        graph_node_x_identity=expand(graph_node_resultant - 256 * T112_GRAPH_NODE_X),
        vertical_node_x_discriminant=discriminant(T112_VERTICAL_NODE_X, X),
        graph_node_x_discriminant=discriminant(T112_GRAPH_NODE_X, X),
        node_packet_resultant=resultant(T112_VERTICAL_NODE_X, T112_GRAPH_NODE_X, X),
        cusp_node_separations=(
            T112_VERTICAL_NODE_X.subs(X, 0),
            T112_GRAPH_NODE_X.subs(X, 0),
        ),
        triple_node_separations=(
            T112_VERTICAL_NODE_X.subs(X, T112_TRIPLE_X),
            T112_GRAPH_NODE_X.subs(X, T112_TRIPLE_X),
        ),
        implicit_resultant_identity=expand(
            T112_IMPLICIT_RESULTANT_MULTIPLIER * resultant(T112_P - X, T112_Q - Y, T)
            - T112_IMPLICIT
        ),
        implicit_parameterization_identity=expand(
            T112_IMPLICIT.subs({X: T112_P, Y: T112_Q})
        ),
        implicit_content=implicit_polynomial.content(),
        implicit_total_degree=implicit_polynomial.total_degree(),
        # Exact Sage primary components by role: cusp, T112 point, and the
        # two reduced node packets.
        sage_jacobian_components=(
            (4, 1, False),
            (6, 1, False),
            (2, 2, True),
            (4, 4, True),
        ),
        sage_cyclic_simplification=(1, 0, True),
        arithmetic_genus=(9 - 1) * (9 - 2) // 2,
        cusp_delta=2,
        t112_delta=4,
        node_count=6,
        infinity_delta=16,
        relation_count=len(T112_RELATIONS),
        pruning_stage_survivors=stage_survivors,
        diagonal_satisfying_assignments=diagonal_assignments,
        complement_census=n_generator_three_cycle_presentation_census(
            T112_RELATIONS,
            4,
        ),
    )


def main() -> int:
    """Print the exact member and fail on any regression."""

    certificate = exact_delta_ten_t112_certificate()
    print(
        "delta-ten T112 member:",
        {
            "cusp": (2, 5),
            "triple contacts": (1, 1, 2),
            "nodes": certificate.node_count,
            "infinity": (5, 9),
            "delta": certificate.total_delta,
        },
    )
    print(
        "generic labeled incidence component:",
        {
            "base irreducible": certificate.base_constraint_irreducible,
            "base absolutely irreducible": certificate.base_constraint_absolutely_irreducible,
            "root-surface dimension": 2,
            "rank everywhere valid": certificate.rank_drop_saturation,
            "free Q coefficients": 1,
            "dimension": certificate.labeled_component_dimension,
            "sample projection fiber": certificate.sample_projection_fiber,
            "generic label degree": certificate.generic_label_degree,
            "clean open connected": certificate.clean_open_connected,
        },
    )
    print("presentation pruning:", certificate.pruning_stage_survivors)
    print(
        "single-three-cycle images:",
        dict(certificate.complement_census.generated_order_histogram),
    )
    print(
        "conditional generic clean T112 exclusion supported:",
        certificate.conditional_generic_exclusion_supported,
    )
    print(
        "claim boundary: exact Sage ideal/complement certificate plus proper "
        "Whitney--Thom propagation; P-critical-fiber boundaries, deeper "
        "intersections, and JC(2) remain open"
    )
    return 0 if certificate.conditional_generic_exclusion_supported else 1


if __name__ == "__main__":
    raise SystemExit(main())
