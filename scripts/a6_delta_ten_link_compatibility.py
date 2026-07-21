"""Hostile compatibility certificate at the conditional delta-ten frontier.

The previous one-Puiseux-pair scan leaves affine degrees ``(4, 9)`` at
collision delta ten.  This module enumerates every centerless ``A6`` torus
representation with the required single-three-cycle meridian and asks whether
the already-forced finite cusp and collision data can fit it.  They can:

* the 720 pairs form two free inner ``A6`` orbits, fused by odd sheet
  relabeling;
* every pair has a unique ordered decomposition into the forced ``T(2, 5)``
  cusp meridians and a complementary collision meridian; and
* the corresponding canonical spin lifts satisfy the ``T(4, 9)`` relation
  and recover the meridian in all 720 cases.

Thus these local, peripheral, Nielsen-class, and spin data do **not** obstruct
the surviving ``(4, 9)`` case.  This is deliberately a hostile consistency
test.  It constructs neither an algebraic curve or cover realizing the data
nor a Keller map, and it does not prove or disprove the Jacobian conjecture.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import cache
from typing import Final, TypeAlias

from scripts.a6_one_pair_infinity import (
    alternating_group_six,
    meridian_exponents,
    permutation_from_cycles,
    permutation_power,
)
from scripts.a6_post_delta_seven_frontier import A6_T49_WITNESS
from scripts.a6_spin_lift import (
    CLIFFORD_MINUS_ONE,
    CLIFFORD_ONE,
    CliffordElement,
    generated_clifford_group,
    three_cycle_lift,
)
from scripts.six_sheet_monodromy import (
    IDENTITY,
    CycleType,
    Permutation,
    PermutationGroup,
    compose,
    cycle_type,
    generated_group,
    inverse,
)

TorusPair: TypeAlias = tuple[Permutation, Permutation]
PairPassport: TypeAlias = tuple[CycleType, CycleType, CycleType]

A6_ORDER: Final = 360
DOUBLE_A6_ORDER: Final = 720
FOUR_TWO: Final[CycleType] = (4, 2)
DOUBLE_THREE: Final[CycleType] = (3, 3)
SINGLE_THREE: Final[CycleType] = (3, 1, 1, 1)
FIVE_CYCLE: Final[CycleType] = (5, 1)

# This fixture joins the forced finite T(2,5) cusp and 3+3 collision in one
# representation.  Its derived (x,y) pair is one of the 720 T(4,9) pairs.
A6_DELTA_TEN_CUSP_FIRST: Final[Permutation] = permutation_from_cycles(
    (3, 4, 5)
)
A6_DELTA_TEN_CUSP_SECOND: Final[Permutation] = permutation_from_cycles(
    (1, 2, 3)
)
A6_DELTA_TEN_COLLISION: Final[Permutation] = permutation_from_cycles(
    (4, 5, 6)
)


def permutation_product(elements: tuple[Permutation, ...]) -> Permutation:
    """Compose an ordered tuple with the convention ``gh = g after h``."""

    result = IDENTITY
    for element in elements:
        result = compose(result, element)
    return result


A6_DELTA_TEN_X: Final[Permutation] = permutation_product(
    (
        A6_DELTA_TEN_CUSP_FIRST,
        A6_DELTA_TEN_CUSP_SECOND,
        A6_DELTA_TEN_COLLISION,
    )
)
A6_DELTA_TEN_Y: Final[Permutation] = inverse(
    compose(A6_DELTA_TEN_CUSP_SECOND, A6_DELTA_TEN_COLLISION)
)


def _conjugate(
    conjugator: Permutation,
    element: Permutation,
) -> Permutation:
    """Return ``conjugator * element * conjugator^-1``."""

    return compose(compose(conjugator, element), inverse(conjugator))


def _conjugate_pair(
    conjugator: Permutation,
    pair: TorusPair,
) -> TorusPair:
    """Relabel both monodromies by the same sheet permutation."""

    x, y = pair
    return (_conjugate(conjugator, x), _conjugate(conjugator, y))


@cache
def _generated_order(*generators: Permutation) -> int:
    """Return the exact order of a generated permutation group."""

    return len(generated_group(generators))


def _orbit_sizes(group: PermutationGroup) -> tuple[int, ...]:
    """Return the sorted orbit sizes in the natural six-point action."""

    remaining = set(range(6))
    sizes: list[int] = []
    while remaining:
        point = min(remaining)
        orbit = {element[point] for element in group}
        sizes.append(len(orbit))
        remaining.difference_update(orbit)
    return tuple(sorted(sizes))


@cache
def _is_five_sheet_a5_cusp(first: Permutation, second: Permutation) -> bool:
    """Check that two cusp meridians generate ``A5`` on five sheets."""

    group = generated_group((first, second))
    return len(group) == 60 and _orbit_sizes(group) == (1, 5)


def _nontrivial_cycles(permutation: Permutation) -> tuple[tuple[int, ...], ...]:
    """Recover oriented nontrivial cycles in canonical one-based notation."""

    seen: set[int] = set()
    cycles: list[tuple[int, ...]] = []
    for start in range(6):
        if start in seen or permutation[start] == start:
            seen.add(start)
            continue
        cycle: list[int] = []
        point = start
        while point not in seen:
            seen.add(point)
            cycle.append(point + 1)
            point = permutation[point]
        cycles.append(tuple(cycle))
    return tuple(cycles)


def _double_three_components(
    permutation: Permutation,
) -> tuple[Permutation, Permutation] | None:
    """Split a ``(3,3)`` permutation into its two oriented factors."""

    cycles = _nontrivial_cycles(permutation)
    if len(cycles) != 2 or any(len(cycle) != 3 for cycle in cycles):
        return None
    first, second = cycles
    return (
        permutation_from_cycles(first),
        permutation_from_cycles(second),
    )


def _three_cycle_tuple(permutation: Permutation) -> tuple[int, int, int]:
    """Return an oriented tuple representing a single three-cycle."""

    cycles = _nontrivial_cycles(permutation)
    if len(cycles) != 1 or len(cycles[0]) != 3:
        msg = "expected a single three-cycle"
        raise ValueError(msg)
    first, second, third = cycles[0]
    return first, second, third


def five_braid_holds(first: Permutation, second: Permutation) -> bool:
    """Check the meridional five-braid relation for a ``T(2,5)`` cusp."""

    left = permutation_product((first, second, first, second, first))
    right = permutation_product((second, first, second, first, second))
    return left == right


@cache
def qualifying_t49_pairs() -> tuple[TorusPair, ...]:
    """Enumerate the 720 generating ``T(4,9) -> A6`` survivor pairs.

    The passport forces ``x`` to type ``(4,2)`` and ``y`` to type ``(3,3)``.
    Restricting to those classes makes the exhaustive scan much smaller than
    iterating over every element whose order merely divides four or nine.
    """

    group = alternating_group_six()
    x_choices = tuple(element for element in group if cycle_type(element) == FOUR_TWO)
    y_choices = tuple(
        element for element in group if cycle_type(element) == DOUBLE_THREE
    )
    u, v = meridian_exponents(4, 9)
    pairs: list[TorusPair] = []
    for x in x_choices:
        for y in y_choices:
            meridian = compose(
                permutation_power(x, u),
                permutation_power(y, v),
            )
            if cycle_type(meridian) != SINGLE_THREE:
                continue
            if _generated_order(x, y) == A6_ORDER:
                pairs.append((x, y))
    return tuple(pairs)


@dataclass(frozen=True, slots=True)
class ForcedCuspCollisionDecomposition:
    """The forced local factorization attached to one torus pair."""

    x: Permutation
    y: Permutation
    cusp_first: Permutation
    cusp_second: Permutation
    collision: Permutation

    @property
    def cusp_relation_holds(self) -> bool:
        """Whether the two cusp meridians satisfy the five-braid relation."""

        return five_braid_holds(self.cusp_first, self.cusp_second)

    @property
    def cusp_group_order(self) -> int:
        """Return the group order generated by the cusp meridians."""

        return _generated_order(self.cusp_first, self.cusp_second)

    @property
    def cusp_orbit_sizes(self) -> tuple[int, ...]:
        """Return the cusp subgroup's natural sheet-orbit profile."""

        return _orbit_sizes(
            generated_group((self.cusp_first, self.cusp_second))
        )

    @property
    def collision_commutes(self) -> bool:
        """Whether the two local branches at the collision commute."""

        return compose(self.cusp_second, self.collision) == compose(
            self.collision,
            self.cusp_second,
        )

    @property
    def total_group_order(self) -> int:
        """Return the group order generated by all finite meridians."""

        return _generated_order(
            self.cusp_first,
            self.cusp_second,
            self.collision,
        )

    @property
    def reconstructs_pair(self) -> bool:
        """Whether the three local meridians reconstruct ``(x,y)`` exactly."""

        return bool(
            permutation_product(
                (self.cusp_first, self.cusp_second, self.collision)
            )
            == self.x
            and inverse(compose(self.cusp_second, self.collision)) == self.y
        )

    @property
    def verified(self) -> bool:
        """Whether every required cusp, collision, and reconstruction check holds."""

        return bool(
            cycle_type(self.cusp_first) == SINGLE_THREE
            and cycle_type(self.cusp_second) == SINGLE_THREE
            and cycle_type(self.collision) == SINGLE_THREE
            and self.cusp_relation_holds
            and self.cusp_group_order == 60
            and self.cusp_orbit_sizes == (1, 5)
            and self.collision_commutes
            and cycle_type(compose(self.cusp_second, self.collision))
            == DOUBLE_THREE
            and self.total_group_order == A6_ORDER
            and self.reconstructs_pair
        )


@cache
def forced_cusp_collision_decompositions(
    x: Permutation,
    y: Permutation,
) -> tuple[ForcedCuspCollisionDecomposition, ...]:
    """Return all ordered forced local decompositions of ``(x,y)``.

    Since ``y^-1`` has type ``(3,3)``, there are only two possible ordered
    assignments of its disjoint three-cycle components to the cusp and
    collision meridians.  Checking both makes uniqueness exhaustive.
    """

    u, v = meridian_exponents(4, 9)
    cusp_first = compose(
        permutation_power(x, u),
        permutation_power(y, v),
    )
    components = _double_three_components(inverse(y))
    if cycle_type(cusp_first) != SINGLE_THREE or components is None:
        return ()

    candidates: list[ForcedCuspCollisionDecomposition] = []
    for cusp_second, collision in (components, components[::-1]):
        candidate = ForcedCuspCollisionDecomposition(
            x=x,
            y=y,
            cusp_first=cusp_first,
            cusp_second=cusp_second,
            collision=collision,
        )
        if candidate.verified:
            candidates.append(candidate)
    return tuple(candidates)


def _clifford_order(element: CliffordElement, *, maximum: int = 24) -> int:
    """Return the exact positive order within a conservative finite bound."""

    power = CLIFFORD_ONE
    for order in range(1, maximum + 1):
        power = power * element
        if power == CLIFFORD_ONE:
            return order
    msg = f"Clifford element order exceeds checked bound {maximum}"
    raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class SpinCompatibility:
    """Canonical spin lifts and exact identities for one decomposition."""

    cusp_first_lift: CliffordElement
    cusp_second_lift: CliffordElement
    collision_lift: CliffordElement
    x_lift: CliffordElement
    y_lift: CliffordElement

    @property
    def local_relations_hold(self) -> bool:
        """Whether the cusp braid and collision commutation lift exactly."""

        first = self.cusp_first_lift
        second = self.cusp_second_lift
        collision = self.collision_lift
        return bool(
            first * second * first * second * first
            == second * first * second * first * second
            and second * collision == collision * second
        )

    @property
    def torus_relation_holds(self) -> bool:
        """Whether ``X^4 = Y^9`` is the nontrivial central spin element."""

        return bool(
            self.x_lift**4 == CLIFFORD_MINUS_ONE
            and self.y_lift**9 == CLIFFORD_MINUS_ONE
        )

    @property
    def meridian_is_canonical(self) -> bool:
        """Whether ``X Y^-2`` recovers the canonical lifted meridian."""

        # This chosen preimage of y has order six, so Y^-2 = Y^4.
        return bool(
            self.y_lift**6 == CLIFFORD_ONE
            and self.x_lift * (self.y_lift**4) == self.cusp_first_lift
        )

    @property
    def verified(self) -> bool:
        """Whether all exact local, torus, and meridian identities hold."""

        return bool(
            self.cusp_first_lift**3 == CLIFFORD_ONE
            and self.cusp_second_lift**3 == CLIFFORD_ONE
            and self.collision_lift**3 == CLIFFORD_ONE
            and self.local_relations_hold
            and self.torus_relation_holds
            and self.meridian_is_canonical
        )


def spin_compatibility(
    decomposition: ForcedCuspCollisionDecomposition,
) -> SpinCompatibility:
    """Build canonical order-three lifts and the order-six lift of ``y``."""

    first_lift = three_cycle_lift(*_three_cycle_tuple(decomposition.cusp_first))
    second_lift = three_cycle_lift(
        *_three_cycle_tuple(decomposition.cusp_second)
    )
    collision_lift = three_cycle_lift(
        *_three_cycle_tuple(decomposition.collision)
    )
    x_lift = first_lift * second_lift * collision_lift
    canonical_y_lift = (second_lift**2) * (collision_lift**2)
    # The negative preimage has order six and makes X^4 = Y^9 = -1.
    y_lift = -canonical_y_lift
    return SpinCompatibility(
        cusp_first_lift=first_lift,
        cusp_second_lift=second_lift,
        collision_lift=collision_lift,
        x_lift=x_lift,
        y_lift=y_lift,
    )


@dataclass(frozen=True, slots=True)
class A6DeltaTenLinkCompatibilityCertificate:
    """Exact census, orbit, decomposition, inversion, and spin outputs."""

    qualifying_pair_count: int
    pair_passport_histogram: tuple[tuple[PairPassport, int], ...]
    inner_orbit_sizes: tuple[int, int]
    inner_orbits_cover_all_pairs: bool
    odd_sheet_orbit_size: int
    odd_sheet_fusion_pass_count: int
    split_five_cycle_counts: tuple[int, int]
    cusp_torus_split_class_agreement_count: int
    decomposition_multiplicity_histogram: tuple[tuple[int, int], ...]
    unique_decomposition_count: int
    cusp_pass_count: int
    collision_pass_count: int
    reconstruction_pass_count: int
    finite_a6_pass_count: int
    simultaneous_inversion_pass_count: int
    simultaneous_inversion_within_inner_orbit_count: int
    x_only_inversion_five_cycle_count: int
    y_only_inversion_five_cycle_count: int
    local_spin_pass_count: int
    torus_spin_pass_count: int
    meridian_spin_pass_count: int
    spin_pass_count: int
    witness: ForcedCuspCollisionDecomposition
    witness_spin_orders: tuple[int, int, int, int, int]
    witness_spin_group_order: int
    witness_lifted_meridian: CliffordElement
    witness_infinity_longitude: CliffordElement
    witness_cusp_longitude: CliffordElement

    @property
    def verified(self) -> bool:
        """Whether every hostile compatibility assertion was certified."""

        expected_passes = DOUBLE_A6_ORDER
        return bool(
            self.qualifying_pair_count == expected_passes
            and self.pair_passport_histogram
            == (((FOUR_TWO, DOUBLE_THREE, SINGLE_THREE), expected_passes),)
            and self.inner_orbit_sizes == (A6_ORDER, A6_ORDER)
            and self.inner_orbits_cover_all_pairs
            and self.odd_sheet_orbit_size == expected_passes
            and self.odd_sheet_fusion_pass_count == expected_passes
            and self.split_five_cycle_counts == (A6_ORDER, A6_ORDER)
            and self.cusp_torus_split_class_agreement_count == expected_passes
            and self.decomposition_multiplicity_histogram
            == ((1, expected_passes),)
            and self.unique_decomposition_count == expected_passes
            and self.cusp_pass_count == expected_passes
            and self.collision_pass_count == expected_passes
            and self.reconstruction_pass_count == expected_passes
            and self.finite_a6_pass_count == expected_passes
            and self.simultaneous_inversion_pass_count == expected_passes
            and self.simultaneous_inversion_within_inner_orbit_count
            == expected_passes
            and self.x_only_inversion_five_cycle_count == expected_passes
            and self.y_only_inversion_five_cycle_count == expected_passes
            and self.local_spin_pass_count == expected_passes
            and self.torus_spin_pass_count == expected_passes
            and self.meridian_spin_pass_count == expected_passes
            and self.spin_pass_count == expected_passes
            and self.witness.verified
            and self.witness_spin_orders == (3, 3, 3, 8, 6)
            and self.witness_spin_group_order == DOUBLE_A6_ORDER
            and self.witness_infinity_longitude == CLIFFORD_MINUS_ONE
            and self.witness_cusp_longitude
            == -(self.witness_lifted_meridian**2)
        )


def _meridian(x: Permutation, y: Permutation) -> Permutation:
    """Return the geometric ``T(4,9)`` meridian ``x y^-2``."""

    u, v = meridian_exponents(4, 9)
    return compose(permutation_power(x, u), permutation_power(y, v))


@cache
def exact_a6_delta_ten_link_compatibility_certificate(
) -> A6DeltaTenLinkCompatibilityCertificate:
    """Build the exhaustive hostile compatibility certificate."""

    pairs = qualifying_t49_pairs()
    pair_set = frozenset(pairs)
    passport_counts: Counter[PairPassport] = Counter(
        (cycle_type(x), cycle_type(y), cycle_type(_meridian(x, y)))
        for x, y in pairs
    )

    a6 = alternating_group_six()
    first_inner_orbit = frozenset(
        _conjugate_pair(conjugator, A6_T49_WITNESS) for conjugator in a6
    )
    odd_relabeling = permutation_from_cycles((1, 2))
    second_seed = _conjugate_pair(odd_relabeling, A6_T49_WITNESS)
    second_inner_orbit = frozenset(
        _conjugate_pair(conjugator, second_seed) for conjugator in a6
    )
    inner_union = first_inner_orbit | second_inner_orbit
    odd_fusion_pass_count = sum(
        _conjugate_pair(odd_relabeling, pair) in second_inner_orbit
        for pair in first_inner_orbit
    ) + sum(
        _conjugate_pair(odd_relabeling, pair) in first_inner_orbit
        for pair in second_inner_orbit
    )

    base_five_cycle = permutation_from_cycles((1, 2, 3, 4, 5))
    first_split_five_class = frozenset(
        _conjugate(conjugator, base_five_cycle) for conjugator in a6
    )
    split_five_counts = Counter(
        compose(x, inverse(y)) in first_split_five_class for x, y in pairs
    )

    decompositions = tuple(
        forced_cusp_collision_decompositions(x, y) for x, y in pairs
    )
    multiplicities: Counter[int] = Counter(map(len, decompositions))
    unique = tuple(entries[0] for entries in decompositions if len(entries) == 1)

    cusp_pass_count = sum(
        decomposition.cusp_relation_holds
        and decomposition.cusp_group_order == 60
        and decomposition.cusp_orbit_sizes == (1, 5)
        for decomposition in unique
    )
    collision_pass_count = sum(
        decomposition.collision_commutes
        and cycle_type(
            compose(decomposition.cusp_second, decomposition.collision)
        )
        == DOUBLE_THREE
        for decomposition in unique
    )
    reconstruction_pass_count = sum(
        decomposition.reconstructs_pair for decomposition in unique
    )
    finite_a6_pass_count = sum(
        decomposition.total_group_order == A6_ORDER for decomposition in unique
    )
    cusp_torus_split_class_agreement_count = sum(
        (
            compose(decomposition.x, inverse(decomposition.y))
            in first_split_five_class
        )
        == (
            compose(decomposition.cusp_first, decomposition.cusp_second)
            in first_split_five_class
        )
        for decomposition in unique
    )

    simultaneous_inversion_pass_count = sum(
        (inverse(x), inverse(y)) in pair_set for x, y in pairs
    )
    simultaneous_inversion_within_inner_orbit_count = sum(
        (inverse(x), inverse(y)) in first_inner_orbit
        for x, y in first_inner_orbit
    ) + sum(
        (inverse(x), inverse(y)) in second_inner_orbit
        for x, y in second_inner_orbit
    )
    x_only_inversion_five_cycle_count = sum(
        cycle_type(_meridian(inverse(x), y)) == FIVE_CYCLE for x, y in pairs
    )
    y_only_inversion_five_cycle_count = sum(
        cycle_type(_meridian(x, inverse(y))) == FIVE_CYCLE for x, y in pairs
    )

    spin_data = tuple(spin_compatibility(decomposition) for decomposition in unique)
    local_spin_pass_count = sum(spin.local_relations_hold for spin in spin_data)
    torus_spin_pass_count = sum(spin.torus_relation_holds for spin in spin_data)
    meridian_spin_pass_count = sum(
        spin.meridian_is_canonical for spin in spin_data
    )
    spin_pass_count = sum(spin.verified for spin in spin_data)

    witness_entries = forced_cusp_collision_decompositions(
        A6_DELTA_TEN_X,
        A6_DELTA_TEN_Y,
    )
    if len(witness_entries) != 1:
        msg = "the explicit hostile witness lost its unique decomposition"
        raise RuntimeError(msg)
    witness = witness_entries[0]
    witness_spin = spin_compatibility(witness)
    witness_spin_orders = (
        _clifford_order(witness_spin.cusp_first_lift),
        _clifford_order(witness_spin.cusp_second_lift),
        _clifford_order(witness_spin.collision_lift),
        _clifford_order(witness_spin.x_lift),
        _clifford_order(witness_spin.y_lift),
    )
    witness_lifted_meridian = witness_spin.x_lift * (witness_spin.y_lift**4)
    # For T(4,9), z=X^4=Y^9 and lambda_infinity=z*m^-36.
    witness_infinity_longitude = (witness_spin.x_lift**4) * (
        witness_lifted_meridian ** ((-36) % 3)
    )
    # For the T(2,5) cusp, lambda_cusp=(RS)^5 R^-10.
    witness_cusp_longitude = (
        (witness_spin.cusp_first_lift * witness_spin.cusp_second_lift) ** 5
    ) * (witness_spin.cusp_first_lift ** ((-10) % 3))

    return A6DeltaTenLinkCompatibilityCertificate(
        qualifying_pair_count=len(pairs),
        pair_passport_histogram=tuple(sorted(passport_counts.items())),
        inner_orbit_sizes=(len(first_inner_orbit), len(second_inner_orbit)),
        inner_orbits_cover_all_pairs=(
            first_inner_orbit.isdisjoint(second_inner_orbit)
            and inner_union == pair_set
        ),
        odd_sheet_orbit_size=len(inner_union),
        odd_sheet_fusion_pass_count=odd_fusion_pass_count,
        split_five_cycle_counts=(
            split_five_counts[True],
            split_five_counts[False],
        ),
        cusp_torus_split_class_agreement_count=(
            cusp_torus_split_class_agreement_count
        ),
        decomposition_multiplicity_histogram=tuple(sorted(multiplicities.items())),
        unique_decomposition_count=len(unique),
        cusp_pass_count=cusp_pass_count,
        collision_pass_count=collision_pass_count,
        reconstruction_pass_count=reconstruction_pass_count,
        finite_a6_pass_count=finite_a6_pass_count,
        simultaneous_inversion_pass_count=simultaneous_inversion_pass_count,
        simultaneous_inversion_within_inner_orbit_count=(
            simultaneous_inversion_within_inner_orbit_count
        ),
        x_only_inversion_five_cycle_count=x_only_inversion_five_cycle_count,
        y_only_inversion_five_cycle_count=y_only_inversion_five_cycle_count,
        local_spin_pass_count=local_spin_pass_count,
        torus_spin_pass_count=torus_spin_pass_count,
        meridian_spin_pass_count=meridian_spin_pass_count,
        spin_pass_count=spin_pass_count,
        witness=witness,
        witness_spin_orders=witness_spin_orders,
        witness_spin_group_order=len(
            generated_clifford_group((witness_spin.x_lift, witness_spin.y_lift))
        ),
        witness_lifted_meridian=witness_lifted_meridian,
        witness_infinity_longitude=witness_infinity_longitude,
        witness_cusp_longitude=witness_cusp_longitude,
    )


def main() -> int:
    """Print the exact hostile certificate and fail if any check breaks."""

    certificate = exact_a6_delta_ten_link_compatibility_certificate()
    print(
        "T(4,9) Nielsen census:",
        {
            "pairs": certificate.qualifying_pair_count,
            "inner orbits": certificate.inner_orbit_sizes,
            "odd-sheet orbit": certificate.odd_sheet_orbit_size,
        },
    )
    print(
        "forced finite decomposition:",
        {
            "multiplicities": certificate.decomposition_multiplicity_histogram,
            "cusp passes": certificate.cusp_pass_count,
            "collision passes": certificate.collision_pass_count,
        },
    )
    print(
        "spin compatibility:",
        {
            "passes": certificate.spin_pass_count,
            "explicit lift group order": certificate.witness_spin_group_order,
            "infinity longitude": certificate.witness_infinity_longitude,
            "cusp longitude": certificate.witness_cusp_longitude,
        },
    )
    print(
        "claim boundary: local/peripheral/spin data do not obstruct; "
        "no curve, cover, or Keller map is constructed"
    )
    print(f"delta-ten compatibility certificate verified: {certificate.verified}")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
