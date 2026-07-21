"""Conditional one-Puiseux-pair obstruction for the ``A6`` survivor.

Assume that the one-dicritical branch has polynomial normalization, one
genuine singular place at infinity, the forced finite ``T(2,5)`` cusp, and
only smooth-branch normalization collisions otherwise.  Genus accounting and
the large affine link reduce small collision budgets to a handful of degree
pairs.  Exact complement presentations then exclude collision delta at most
two, and further link-at-infinity censuses exclude deltas three and four.
The first surviving coarse degree pair occurs at delta five.

This is a conditional, computer-assisted frontier theorem.  It does not prove
that every Keller branch satisfies the one-pair or finite-singularity
hypotheses, and it does not solve the plane Jacobian conjecture.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache
from itertools import product
from math import gcd
from typing import Final

from sympy import (
    Expr,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    factor,
    invert,
    rem,
    resultant,
    together,
)

from scripts.six_sheet_monodromy import (
    IDENTITY,
    TRANSITIVE_GROUPS,
    CycleType,
    Permutation,
    compose,
    cycle_type,
    generated_group,
    inverse,
)

T: Final = Symbol("t")
U: Final = Symbol("u")
C: Final = Symbol("c")
S: Final = Symbol("s")
X: Final = Symbol("X")
Y: Final = Symbol("Y")

FAMILY_P: Final = T**2 + T**3
FAMILY_Q: Final = C * T**4 + T**5
FAMILY_IMPLICIT: Final = (
    -X**5
    - C**3 * X**4
    + C**2 * X**4
    - 3 * C * X**3 * Y
    + 5 * X**3 * Y
    + 2 * C**2 * X**2 * Y
    - 2 * C * X**2 * Y
    + 4 * C * X * Y**2
    - 5 * X * Y**2
    + Y**3
    - C * Y**2
    + Y**2
)

COLLISION_POLYNOMIAL: Final = S**2 + (C + 1) * S + 2 * C - 1
PAIR_DISCRIMINANT: Final = -S * (3 * S + 4)
INVALID_FAMILY_PARAMETERS: Final = (
    (1, 2),
    (5, 6),
    (1, 1),
)
EXCEPTIONAL_TANGENCY_PARAMETER: Final = (5, 1)

# Sage 10.8, Curve(FAMILY_IMPLICIT.subs(C, 0)).fundamental_group(
# simplified=False).  Indices 1..3 correspond to geometric fiber meridians.
GENERIC_FAMILY_RELATIONS: Final = (
    (3, 2, -3, -2),
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -3, -2, 1, 2, 3, 2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
)

# The same unsimplified computation at the valid contact-two fiber C=5.
CONTACT_TWO_RELATIONS: Final = (
    (
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
        1,
        2,
        3,
        -2,
        -1,
        2,
        1,
        2,
        -3,
        -2,
        -1,
        -2,
        1,
        2,
        -3,
        -2,
        -1,
        -2,
        1,
        2,
    ),
    (-3, -2, -1, -2, -1, 2, 1, 2, 3, -2, -1, 2, 1, 2),
    (2, 1, 2, 1, -2, -1, -2, -1),
)


@dataclass(frozen=True, slots=True)
class OnePairDegreeCandidate:
    """Affine degree data and the corresponding projective infinity pair."""

    affine_degrees: tuple[int, int]
    projective_pair: tuple[int, int]

    @property
    def affine_link_pair(self) -> tuple[int, int]:
        """Return the large-sphere link pair, which uses the affine degrees."""

        return self.affine_degrees


def one_pair_degree_candidates(
    collision_delta: int,
    *,
    require_singular_infinity: bool = True,
) -> tuple[OnePairDegreeCandidate, ...]:
    """Solve ``(d-1)(a-1)=2(Delta+2)`` with ``a<d`` and coprime degrees."""

    if collision_delta < 1:
        msg = "the normalization-collision delta must be positive"
        raise ValueError(msg)
    target = 2 * (collision_delta + 2)
    candidates: list[OnePairDegreeCandidate] = []
    for a_minus_one in range(1, target + 1):
        if target % a_minus_one:
            continue
        d_minus_one = target // a_minus_one
        a = a_minus_one + 1
        d = d_minus_one + 1
        if not a < d or gcd(a, d) != 1:
            continue
        projective_pair = (d - a, d)
        if require_singular_infinity and projective_pair[0] <= 1:
            continue
        candidates.append(
            OnePairDegreeCandidate(
                affine_degrees=(a, d),
                projective_pair=projective_pair,
            )
        )
    return tuple(candidates)


def permutation_power(permutation: Permutation, exponent: int) -> Permutation:
    """Raise a six-point permutation to an arbitrary integer power."""

    if exponent < 0:
        return permutation_power(inverse(permutation), -exponent)
    result = IDENTITY
    base = permutation
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = compose(result, base)
        base = compose(base, base)
        remaining >>= 1
    return result


def permutation_from_cycles(*cycles: tuple[int, ...]) -> Permutation:
    """Build a permutation from disjoint cycles on the one-based sheet set."""

    image = list(IDENTITY)
    support: set[int] = set()
    for cycle in cycles:
        if len(cycle) < 2:
            msg = "cycles must have length at least two"
            raise ValueError(msg)
        zero_based = tuple(point - 1 for point in cycle)
        if any(point < 0 or point >= 6 for point in zero_based):
            msg = "cycle point outside 1..6"
            raise ValueError(msg)
        if support.intersection(zero_based):
            msg = "cycles must be disjoint"
            raise ValueError(msg)
        support.update(zero_based)
        for source, target in zip(zero_based, zero_based[1:] + zero_based[:1]):
            image[source] = target
    return tuple(image)


A6_T35_WITNESS: Final[tuple[Permutation, Permutation]] = (
    permutation_from_cycles((1, 6, 3), (2, 5, 4)),
    permutation_from_cycles((1, 5, 4, 3, 2)),
)
A6_T38_WITNESS: Final[tuple[Permutation, Permutation]] = (
    permutation_from_cycles((1, 2, 3), (4, 5, 6)),
    permutation_from_cycles((1, 2), (3, 4, 6, 5)),
)


def meridian_exponents(a: int, d: int) -> tuple[int, int]:
    """Return ``(u,v)`` with ``d*u+a*v=1`` for the ``T(a,d)`` meridian."""

    if a < 1 or d < 1 or gcd(a, d) != 1:
        msg = "the torus-link degrees must be positive and coprime"
        raise ValueError(msg)
    old_r, r = d, a
    old_u, u = 1, 0
    old_v, v = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_u, u = u, old_u - quotient * u
        old_v, v = v, old_v - quotient * v
    return old_u, old_v


@cache
def alternating_group_six() -> tuple[Permutation, ...]:
    """Enumerate ``A6`` from its checked transitive-group generators."""

    fixture = next(
        entry for entry in TRANSITIVE_GROUPS if entry.identifier == "6T15"
    )
    group = generated_group(fixture.generators)
    if len(group) != 360:
        msg = "the stored 6T15 generators did not recover A6"
        raise RuntimeError(msg)
    return tuple(sorted(group))


@dataclass(frozen=True, slots=True)
class TorusA6Census:
    """Exact centerless ``T(a,d)`` quotient census inside ``A6``."""

    pair: tuple[int, int]
    meridian_exponents: tuple[int, int]
    x_choices: int
    y_choices: int
    generated_order_histogram: tuple[tuple[int, int], ...]
    a6_meridian_histogram: tuple[tuple[CycleType, int], ...]

    @property
    def single_three_cycle_a6_pairs(self) -> int:
        """Count generating pairs with the required one-dicritical meridian."""

        return dict(self.a6_meridian_histogram).get((3, 1, 1, 1), 0)


@cache
def torus_a6_census(a: int, d: int) -> TorusA6Census:
    """Enumerate all pairs satisfying ``x^a=y^d=1`` in centerless ``A6``."""

    u, v = meridian_exponents(a, d)
    group = alternating_group_six()
    x_choices = tuple(g for g in group if permutation_power(g, a) == IDENTITY)
    y_choices = tuple(g for g in group if permutation_power(g, d) == IDENTITY)
    orders: Counter[int] = Counter()
    meridians: Counter[CycleType] = Counter()
    for x in x_choices:
        for y in y_choices:
            generated = generated_group((x, y))
            order = len(generated)
            orders[order] += 1
            if order != 360:
                continue
            meridian = compose(
                permutation_power(x, u),
                permutation_power(y, v),
            )
            meridians[cycle_type(meridian)] += 1
    return TorusA6Census(
        pair=(a, d),
        meridian_exponents=(u, v),
        x_choices=len(x_choices),
        y_choices=len(y_choices),
        generated_order_histogram=tuple(sorted(orders.items())),
        a6_meridian_histogram=tuple(sorted(meridians.items())),
    )


def evaluate_signed_word(
    word: tuple[int, ...],
    images: Sequence[Permutation],
) -> Permutation:
    """Evaluate a signed Tietze word in the given generator images."""

    result = IDENTITY
    for letter in word:
        generator = images[abs(letter) - 1]
        if letter < 0:
            generator = inverse(generator)
        result = compose(result, generator)
    return result


@dataclass(frozen=True, slots=True)
class ThreeCyclePresentationCensus:
    """Images of geometric meridians constrained to single 3-cycles."""

    assignments: int
    satisfying_assignments: int
    generated_order_histogram: tuple[tuple[int, int], ...]

    @property
    def a6_assignments(self) -> int:
        """Return the number of representations with image all of ``A6``."""

        return dict(self.generated_order_histogram).get(360, 0)


@cache
def n_generator_three_cycle_presentation_census(
    relations: tuple[tuple[int, ...], ...],
    generator_count: int,
) -> ThreeCyclePresentationCensus:
    """Exhaust single-three-cycle images with incremental relation pruning."""

    if generator_count < 1:
        msg = "the presentation must have at least one generator"
        raise ValueError(msg)
    for relation in relations:
        for letter in relation:
            if letter == 0 or abs(letter) > generator_count:
                msg = "a relation references a generator outside the presentation"
                raise ValueError(msg)

    three_cycles = tuple(
        element
        for element in alternating_group_six()
        if cycle_type(element) == (3, 1, 1, 1)
    )
    partial_assignments: tuple[tuple[Permutation, ...], ...] = ((),)
    for assigned_count in range(1, generator_count + 1):
        newly_decidable = tuple(
            relation
            for relation in relations
            if relation and max(abs(letter) for letter in relation) == assigned_count
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

    orders: Counter[int] = Counter()
    for images in partial_assignments:
        orders[len(generated_group(images))] += 1
    return ThreeCyclePresentationCensus(
        assignments=len(three_cycles) ** generator_count,
        satisfying_assignments=len(partial_assignments),
        generated_order_histogram=tuple(sorted(orders.items())),
    )


@cache
def three_cycle_presentation_census(
    relations: tuple[tuple[int, ...], ...],
) -> ThreeCyclePresentationCensus:
    """Exhaust all 40^3 three-cycle assignments to a three-meridian group."""

    return n_generator_three_cycle_presentation_census(relations, 3)


@dataclass(frozen=True, slots=True)
class A6OnePairAlgebraCertificate:
    """Exact identities controlling the complete ``Delta=2`` family."""

    implicit_parameterization: Expr
    cusp_remainder: Expr
    collision_polynomial_discriminant: Expr
    pair_discriminant: Expr
    tangency_resultant: Expr
    extra_critical_resultant: Expr
    generic_implicit: Expr
    contact_two_implicit: Expr
    contact_two_curvature_remainder: Expr

    @property
    def verified(self) -> bool:
        """Check the symbolic exceptional-parameter and contact identities."""

        return bool(
            self.implicit_parameterization == 0
            and expand(
                self.cusp_remainder
                - (-C * T**6 + (1 - 2 * C) * T**5)
            )
            == 0
            and expand(
                self.collision_polynomial_discriminant - (C - 1) * (C - 5)
            )
            == 0
            and expand(self.pair_discriminant + S * (3 * S + 4)) == 0
            and expand(
                self.tangency_resultant
                + (C - 5) * (C - 1) * (2 * C - 1) ** 3 * (6 * C - 5)
            )
            == 0
            and expand(self.extra_critical_resultant - 2 * (6 * C - 5)) == 0
            and expand(
                self.generic_implicit
                - (-X**5 + 5 * X**3 * Y - 5 * X * Y**2 + Y**3 + Y**2)
            )
            == 0
            and expand(
                self.contact_two_implicit
                - (
                    -X**5
                    - 100 * X**4
                    - 10 * X**3 * Y
                    + 40 * X**2 * Y
                    + 15 * X * Y**2
                    + Y**3
                    - 4 * Y**2
                )
            )
            == 0
            and self.contact_two_curvature_remainder == 0
        )


def exact_a6_one_pair_algebra_certificate() -> A6OnePairAlgebraCertificate:
    """Build the exact normal-form, collision, and contact-two certificate."""

    collision_tangency = -S**2 * (S + 1) * (
        12 * C * S + 20 * C + 15 * S**2 + 15 * S - 10
    )
    tangent_resultant = factor(
        resultant(COLLISION_POLYNOMIAL, collision_tangency, S)
    )
    p_at_five = FAMILY_P
    q_at_five = FAMILY_Q.subs(C, 5)
    curvature = factor(
        diff(diff(q_at_five, T) / diff(p_at_five, T), T)
        / diff(p_at_five, T)
    )
    expected_curvature = -(9 * T + 61) / 128
    curvature_difference = together(curvature - expected_curvature)
    curvature_remainder = rem(
        curvature_difference.as_numer_denom()[0],
        T**2 + 3 * T + 6,
        T,
    )
    # Make the denominator invertibility explicit at the collision pair.
    denominator = curvature.as_numer_denom()[1]
    invert(rem(denominator, T**2 + 3 * T + 6, T), T**2 + 3 * T + 6)
    return A6OnePairAlgebraCertificate(
        implicit_parameterization=expand(
            FAMILY_IMPLICIT.subs({X: FAMILY_P, Y: FAMILY_Q})
        ),
        cusp_remainder=expand(FAMILY_Q - C * FAMILY_P**2),
        collision_polynomial_discriminant=factor(
            discriminant(COLLISION_POLYNOMIAL, S)
        ),
        pair_discriminant=PAIR_DISCRIMINANT,
        tangency_resultant=tangent_resultant,
        extra_critical_resultant=expand(
            resultant(2 + 3 * T, 4 * C + 5 * T, T)
        ),
        generic_implicit=expand(FAMILY_IMPLICIT.subs(C, 0)),
        contact_two_implicit=expand(FAMILY_IMPLICIT.subs(C, 5)),
        contact_two_curvature_remainder=expand(curvature_remainder),
    )


@dataclass(frozen=True, slots=True)
class A6OnePairInfinityCertificate:
    """Conditional exact exclusion through collision delta four."""

    delta_one_candidates: tuple[OnePairDegreeCandidate, ...]
    delta_two_candidates: tuple[OnePairDegreeCandidate, ...]
    delta_three_candidates: tuple[OnePairDegreeCandidate, ...]
    delta_four_candidates: tuple[OnePairDegreeCandidate, ...]
    delta_five_candidates: tuple[OnePairDegreeCandidate, ...]
    t27_census: TorusA6Census
    t29_census: TorusA6Census
    t35_census: TorusA6Census
    t2_11_census: TorusA6Census
    t2_13_census: TorusA6Census
    t3_7_census: TorusA6Census
    t2_15_census: TorusA6Census
    t3_8_census: TorusA6Census
    generic_family_census: ThreeCyclePresentationCensus
    contact_two_census: ThreeCyclePresentationCensus
    algebra: A6OnePairAlgebraCertificate

    @property
    def verified(self) -> bool:
        """Whether every exact arithmetic and finite-group checkpoint agrees."""

        return (
            tuple(c.affine_degrees for c in self.delta_one_candidates) == ((2, 7),)
            and tuple(c.affine_degrees for c in self.delta_two_candidates)
            == ((2, 9), (3, 5))
            and tuple(c.affine_degrees for c in self.delta_three_candidates)
            == ((2, 11),)
            and tuple(c.affine_degrees for c in self.delta_four_candidates)
            == ((2, 13), (3, 7))
            and tuple(c.affine_degrees for c in self.delta_five_candidates)
            == ((2, 15), (3, 8))
            and self.t27_census.single_three_cycle_a6_pairs == 0
            and self.t29_census.single_three_cycle_a6_pairs == 0
            and self.t35_census.single_three_cycle_a6_pairs == 720
            and self.t2_11_census.single_three_cycle_a6_pairs == 0
            and self.t2_13_census.single_three_cycle_a6_pairs == 0
            and self.t3_7_census.single_three_cycle_a6_pairs == 0
            and self.t2_15_census.single_three_cycle_a6_pairs == 0
            and self.t3_8_census.single_three_cycle_a6_pairs == 720
            and self.generic_family_census.assignments == 40**3
            and self.generic_family_census.satisfying_assignments == 40
            and self.generic_family_census.generated_order_histogram == ((3, 40),)
            and self.contact_two_census.assignments == 40**3
            and self.contact_two_census.satisfying_assignments == 760
            and self.contact_two_census.generated_order_histogram
            == ((3, 40), (60, 720))
            and self.generic_family_census.a6_assignments == 0
            and self.contact_two_census.a6_assignments == 0
            and self.algebra.verified
        )


def exact_a6_one_pair_infinity_certificate() -> A6OnePairInfinityCertificate:
    """Construct the exclusion and the first surviving degree-pair certificate."""

    return A6OnePairInfinityCertificate(
        delta_one_candidates=one_pair_degree_candidates(1),
        delta_two_candidates=one_pair_degree_candidates(2),
        delta_three_candidates=one_pair_degree_candidates(3),
        delta_four_candidates=one_pair_degree_candidates(4),
        delta_five_candidates=one_pair_degree_candidates(5),
        t27_census=torus_a6_census(2, 7),
        t29_census=torus_a6_census(2, 9),
        t35_census=torus_a6_census(3, 5),
        t2_11_census=torus_a6_census(2, 11),
        t2_13_census=torus_a6_census(2, 13),
        t3_7_census=torus_a6_census(3, 7),
        t2_15_census=torus_a6_census(2, 15),
        t3_8_census=torus_a6_census(3, 8),
        generic_family_census=three_cycle_presentation_census(
            GENERIC_FAMILY_RELATIONS
        ),
        contact_two_census=three_cycle_presentation_census(
            CONTACT_TWO_RELATIONS
        ),
        algebra=exact_a6_one_pair_algebra_certificate(),
    )


def main() -> int:
    """Print the conditional frontier certificate and fail on regression."""

    certificate = exact_a6_one_pair_infinity_certificate()
    print(
        "one-pair degree candidates:",
        {
            "Delta=1": [
                candidate.affine_degrees
                for candidate in certificate.delta_one_candidates
            ],
            "Delta=2": [
                candidate.affine_degrees
                for candidate in certificate.delta_two_candidates
            ],
            "Delta=3": [
                candidate.affine_degrees
                for candidate in certificate.delta_three_candidates
            ],
            "Delta=4": [
                candidate.affine_degrees
                for candidate in certificate.delta_four_candidates
            ],
            "Delta=5": [
                candidate.affine_degrees
                for candidate in certificate.delta_five_candidates
            ],
        },
    )
    print(
        "large-link A6-generating single-3 meridians:",
        {
            "T(2,7)": certificate.t27_census.single_three_cycle_a6_pairs,
            "T(2,9)": certificate.t29_census.single_three_cycle_a6_pairs,
            "T(3,5)": certificate.t35_census.single_three_cycle_a6_pairs,
            "T(2,11)": certificate.t2_11_census.single_three_cycle_a6_pairs,
            "T(2,13)": certificate.t2_13_census.single_three_cycle_a6_pairs,
            "T(3,7)": certificate.t3_7_census.single_three_cycle_a6_pairs,
            "T(2,15)": certificate.t2_15_census.single_three_cycle_a6_pairs,
            "T(3,8)": certificate.t3_8_census.single_three_cycle_a6_pairs,
        },
    )
    print(
        "Delta=2 complement images:",
        {
            "generic": dict(
                certificate.generic_family_census.generated_order_histogram
            ),
            "contact-two": dict(
                certificate.contact_two_census.generated_order_histogram
            ),
        },
    )
    print(f"conditional A6 one-pair certificate verified: {certificate.verified}")
    print("conclusion under stated hypotheses: collision delta is at least five")
    print("first surviving degree pair at Delta=5: (a,d)=(3,8)")
    print("claim boundary: one-pair/computer-assisted; JC(2) remains open")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
