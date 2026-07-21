"""Exact infinity-group obstruction for the one-pair ``S6`` survivor.

Assume an irreducible branch curve has polynomial normalization
``t |-> (P(t), Q(t))`` with degrees ``m < n``, one genuine singular pair at
infinity, and full six-sheet monodromy ``S6`` with transposition inertia.
The affine link at infinity is the torus knot ``T(m,n)``.  Its group

``<a,b | a^m = b^n>``

surjects onto the affine complement group.  In an ``S6`` quotient the common
central power dies, and a Bezout meridian ``a^u b^v`` must be a
transposition.  This dependency-free checker exhausts those finite
conditions.  Combined with the width-five lower bound, it excludes every
singular one-pair degree slot with ``n <= 11``.  The first link-group
survivor is ``(m,n)=(5,12)``, for which an explicit generating witness is
checked exactly.

This is an infinity-topology obstruction under the stated one-pair
hypotheses.  It does not construct a branch curve or Keller map, address
multi-pair infinity, or solve the plane Jacobian conjecture.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import cache
from itertools import permutations
from math import gcd
from typing import Final

from scripts.six_sheet_monodromy import (
    IDENTITY,
    CycleType,
    Permutation,
    compose,
    cycle_type,
    generated_group,
    inverse,
)

SHEET_COUNT: Final = 6
S6_ORDER: Final = 720
MINIMUM_TRANSPOSITION_WIDTH: Final = SHEET_COUNT - 1
TRANSPOSITION_TYPE: Final[CycleType] = (2, 1, 1, 1, 1)

SMALL_ONE_PAIR_DEGREES: Final = (
    (5, 7),
    (5, 8),
    (5, 9),
    (7, 9),
    (7, 10),
    (5, 11),
    (6, 11),
    (7, 11),
    (8, 11),
    (9, 11),
)
FIRST_SURVIVING_TARGET_DEGREES: Final = ((5, 12), (7, 12))

# pair -> (x choices, y choices, meridian-compatible pairs, order histogram)
EXPECTED_CENSUS_DATA: Final[
    dict[
        tuple[int, int],
        tuple[int, int, int, tuple[tuple[int, int], ...]],
    ]
] = {
    (5, 7): (145, 1, 0, ()),
    (5, 8): (145, 256, 735, ((2, 15), (120, 720))),
    (5, 9): (145, 81, 0, ()),
    (7, 9): (1, 81, 0, ()),
    (7, 10): (1, 220, 15, ((2, 15),)),
    (5, 11): (145, 1, 0, ()),
    (6, 11): (396, 1, 15, ((2, 15),)),
    (7, 11): (1, 1, 0, ()),
    (8, 11): (256, 1, 15, ((2, 15),)),
    (9, 11): (81, 1, 0, ()),
    (5, 12): (145, 576, 2175, ((2, 15), (120, 1440), (720, 720))),
    (7, 12): (1, 576, 15, ((2, 15),)),
}


@dataclass(frozen=True, slots=True)
class OnePairDegreeCandidate:
    """Affine normalization degrees and their projective infinity pair."""

    affine_degrees: tuple[int, int]
    projective_pair: tuple[int, int]

    @property
    def affine_link_pair(self) -> tuple[int, int]:
        """Return the torus-knot pair seen on a large affine sphere."""

        return self.affine_degrees


def one_pair_degree_candidates(
    maximum_target_degree: int,
    *,
    minimum_width: int = MINIMUM_TRANSPOSITION_WIDTH,
) -> tuple[OnePairDegreeCandidate, ...]:
    """Enumerate ``m<n``, ``gcd(m,n)=1``, and ``n-m>=2``.

    The width theorem supplies ``m>=5`` in the six-sheet transposition case.
    The gap and coprimality encode a genuine single projective pair
    ``(n-m,n)`` rather than a smooth or multi-pair end.
    """

    if maximum_target_degree < 1:
        msg = "the maximum target degree must be positive"
        raise ValueError(msg)
    if minimum_width < 1:
        msg = "the minimum width must be positive"
        raise ValueError(msg)
    candidates: list[OnePairDegreeCandidate] = []
    for n in range(minimum_width + 2, maximum_target_degree + 1):
        for m in range(minimum_width, n - 1):
            if gcd(m, n) != 1:
                continue
            candidates.append(
                OnePairDegreeCandidate(
                    affine_degrees=(m, n),
                    projective_pair=(n - m, n),
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
    """Build a permutation from disjoint one-based cycles."""

    image = list(IDENTITY)
    support: set[int] = set()
    for cycle in cycles:
        if len(cycle) < 2:
            msg = "cycles must have length at least two"
            raise ValueError(msg)
        zero_based = tuple(point - 1 for point in cycle)
        if any(point < 0 or point >= SHEET_COUNT for point in zero_based):
            msg = "cycle point outside 1..6"
            raise ValueError(msg)
        if support.intersection(zero_based):
            msg = "cycles must be disjoint"
            raise ValueError(msg)
        support.update(zero_based)
        for source, target in zip(zero_based, zero_based[1:] + zero_based[:1]):
            image[source] = target
    return tuple(image)


def meridian_exponents(m: int, n: int) -> tuple[int, int]:
    """Return ``(u,v)`` with ``n*u+m*v=1`` for a ``T(m,n)`` meridian."""

    if m < 1 or n < 1 or gcd(m, n) != 1:
        msg = "the torus-knot degrees must be positive and coprime"
        raise ValueError(msg)
    old_r, r = n, m
    old_u, u = 1, 0
    old_v, v = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_u, u = u, old_u - quotient * u
        old_v, v = v, old_v - quotient * v
    return old_u, old_v


@cache
def symmetric_group_six() -> tuple[Permutation, ...]:
    """Return all 720 permutations in the natural six-point action."""

    group = tuple(permutations(range(SHEET_COUNT)))
    if len(group) != S6_ORDER:
        msg = "the six-point permutation enumeration did not recover S6"
        raise RuntimeError(msg)
    return group


@dataclass(frozen=True, slots=True)
class TorusS6Census:
    """Exact centerless ``T(m,n)`` quotient census inside ``S6``."""

    pair: tuple[int, int]
    projective_pair: tuple[int, int]
    meridian_exponents: tuple[int, int]
    x_choices: int
    y_choices: int
    meridian_compatible_pairs: int
    generated_order_histogram: tuple[tuple[int, int], ...]

    @property
    def s6_pairs(self) -> int:
        """Count compatible ordered pairs that generate all of ``S6``."""

        return dict(self.generated_order_histogram).get(S6_ORDER, 0)

    @property
    def verified(self) -> bool:
        """Whether this census matches the stored exact regression data."""

        expected = EXPECTED_CENSUS_DATA.get(self.pair)
        actual = (
            self.x_choices,
            self.y_choices,
            self.meridian_compatible_pairs,
            self.generated_order_histogram,
        )
        return expected == actual


@cache
def torus_s6_census(m: int, n: int) -> TorusS6Census:
    """Exhaust all centerless torus-knot images with transposition meridian."""

    if not m < n:
        msg = "the affine degrees must satisfy m < n"
        raise ValueError(msg)
    u, v = meridian_exponents(m, n)
    group = symmetric_group_six()
    x_choices = tuple(
        element
        for element in group
        if permutation_power(element, m) == IDENTITY
    )
    y_choices = tuple(
        element
        for element in group
        if permutation_power(element, n) == IDENTITY
    )
    compatible = 0
    orders: Counter[int] = Counter()
    for x in x_choices:
        x_meridian_factor = permutation_power(x, u)
        for y in y_choices:
            meridian = compose(x_meridian_factor, permutation_power(y, v))
            if cycle_type(meridian) != TRANSPOSITION_TYPE:
                continue
            compatible += 1
            orders[len(generated_group((x, y)))] += 1
    return TorusS6Census(
        pair=(m, n),
        projective_pair=(n - m, n),
        meridian_exponents=(u, v),
        x_choices=len(x_choices),
        y_choices=len(y_choices),
        meridian_compatible_pairs=compatible,
        generated_order_histogram=tuple(sorted(orders.items())),
    )


S6_T5_12_WITNESS: Final[tuple[Permutation, Permutation]] = (
    permutation_from_cycles((2, 3, 4, 5, 6)),
    permutation_from_cycles((1, 2, 5, 3, 6, 4)),
)


@dataclass(frozen=True, slots=True)
class S6T512WitnessCertificate:
    """Exact witness that the ``T(5,12)`` infinity test is sharp."""

    x_fifth_power: Permutation
    y_twelfth_power: Permutation
    meridian: Permutation
    meridian_cycle_type: CycleType
    generated_order: int
    conjugate_star: frozenset[Permutation]

    @property
    def verified(self) -> bool:
        """Whether the relation, meridian, and hand generation proof hold."""

        expected_star = frozenset(
            permutation_from_cycles((1, point)) for point in range(2, 7)
        )
        return (
            self.x_fifth_power == IDENTITY
            and self.y_twelfth_power == IDENTITY
            and self.meridian == permutation_from_cycles((1, 2))
            and self.meridian_cycle_type == TRANSPOSITION_TYPE
            and self.generated_order == S6_ORDER
            and self.conjugate_star == expected_star
        )


def exact_t5_12_witness_certificate() -> S6T512WitnessCertificate:
    """Check the first surviving torus-knot quotient and its star generators."""

    x, y = S6_T5_12_WITNESS
    u, v = meridian_exponents(5, 12)
    meridian = compose(permutation_power(x, u), permutation_power(y, v))
    conjugate_star = frozenset(
        compose(
            compose(permutation_power(x, exponent), meridian),
            permutation_power(x, -exponent),
        )
        for exponent in range(5)
    )
    return S6T512WitnessCertificate(
        x_fifth_power=permutation_power(x, 5),
        y_twelfth_power=permutation_power(y, 12),
        meridian=meridian,
        meridian_cycle_type=cycle_type(meridian),
        generated_order=len(generated_group((x, y))),
        conjugate_star=conjugate_star,
    )


@dataclass(frozen=True, slots=True)
class S6OnePairInfinityCertificate:
    """Complete small-degree elimination and first-survivor certificate."""

    small_degree_pairs: tuple[tuple[int, int], ...]
    small_censuses: tuple[TorusS6Census, ...]
    degree_twelve_pairs: tuple[tuple[int, int], ...]
    degree_twelve_censuses: tuple[TorusS6Census, ...]
    witness: S6T512WitnessCertificate

    @property
    def all_degrees_through_eleven_excluded(self) -> bool:
        """Whether no small candidate admits an exact ``S6`` quotient."""

        return all(census.s6_pairs == 0 for census in self.small_censuses)

    @property
    def first_surviving_pair(self) -> tuple[int, int] | None:
        """Return the first degree-twelve pair with a valid ``S6`` quotient."""

        for census in self.degree_twelve_censuses:
            if census.s6_pairs:
                return census.pair
        return None

    @property
    def verified(self) -> bool:
        """Whether enumeration, elimination, and the sharp witness all agree."""

        return (
            self.small_degree_pairs == SMALL_ONE_PAIR_DEGREES
            and self.degree_twelve_pairs == FIRST_SURVIVING_TARGET_DEGREES
            and all(census.verified for census in self.small_censuses)
            and all(census.verified for census in self.degree_twelve_censuses)
            and self.all_degrees_through_eleven_excluded
            and self.first_surviving_pair == (5, 12)
            and self.witness.verified
        )


@cache
def exact_s6_one_pair_infinity_certificate() -> S6OnePairInfinityCertificate:
    """Build the exact one-pair ``S6`` infinity certificate."""

    small_pairs = tuple(
        candidate.affine_degrees for candidate in one_pair_degree_candidates(11)
    )
    degree_twelve_pairs = tuple(
        candidate.affine_degrees
        for candidate in one_pair_degree_candidates(12)
        if candidate.affine_degrees[1] == 12
    )
    return S6OnePairInfinityCertificate(
        small_degree_pairs=small_pairs,
        small_censuses=tuple(torus_s6_census(*pair) for pair in small_pairs),
        degree_twelve_pairs=degree_twelve_pairs,
        degree_twelve_censuses=tuple(
            torus_s6_census(*pair) for pair in degree_twelve_pairs
        ),
        witness=exact_t5_12_witness_certificate(),
    )


def main() -> int:
    """Print the exact small-degree table and fail on any regression."""

    certificate = exact_s6_one_pair_infinity_certificate()
    for census in certificate.small_censuses + certificate.degree_twelve_censuses:
        print(
            census.pair,
            {
                "projective pair": census.projective_pair,
                "meridian exponents": census.meridian_exponents,
                "compatible": census.meridian_compatible_pairs,
                "orders": dict(census.generated_order_histogram),
                "S6": census.s6_pairs,
            },
        )
    print(
        "first one-pair infinity-group survivor:",
        certificate.first_surviving_pair,
    )
    print(f"S6 one-pair infinity certificate verified: {certificate.verified}")
    print("claim boundary: the one-pair stratum only; (5,12) is not a cover")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
