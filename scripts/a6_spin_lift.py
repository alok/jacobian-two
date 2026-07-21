"""Exact Clifford certificate for the surviving ``A6`` local relations.

The natural six-point action embeds ``A6`` in ``SO(6)``.  Pulling back the
spin cover gives the nonsplit central extension ``2.A6``.  This module checks
that the forced ``T(2,5)`` cusp relation and the separate ``3+3`` collision
relation lift simultaneously, and that all-3-cycle product-one completions
exist with both spin signs.  It is a hostile lifting model, not a Keller
compactification or a determination of the actual infinity monodromy.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Final

from scripts.six_sheet_monodromy import (
    IDENTITY,
    Permutation,
    compose,
    cycle_type,
    generated_group,
    inverse,
    is_transitive,
)

CLIFFORD_DIMENSION: Final = 6


def _normalized_terms(
    coefficients: dict[int, Fraction],
) -> tuple[tuple[int, Fraction], ...]:
    """Drop zero coefficients and return a canonical immutable ordering."""

    return tuple(
        sorted(
            (mask, coefficient)
            for mask, coefficient in coefficients.items()
            if coefficient != 0
        )
    )


def _blade_product(first: int, second: int) -> tuple[int, int]:
    """Multiply two Euclidean Clifford basis blades exactly.

    A bit mask stores an increasing product of basis vectors.  Since
    ``e_i^2=1``, the output mask is XOR; the sign is the parity of the
    inversions between the two ordered blades.
    """

    inversions = 0
    remaining = first
    while remaining:
        lowest = remaining & -remaining
        inversions += (second & (lowest - 1)).bit_count()
        remaining ^= lowest
    return (-1 if inversions % 2 else 1, first ^ second)


@dataclass(frozen=True, slots=True)
class CliffordElement:
    """An exact rational element of the Euclidean Clifford algebra ``Cl_6``."""

    terms: tuple[tuple[int, Fraction], ...]

    @classmethod
    def from_dict(cls, coefficients: dict[int, Fraction]) -> CliffordElement:
        """Build a normalized element from a sparse coefficient dictionary."""

        return cls(_normalized_terms(coefficients))

    @classmethod
    def scalar(cls, value: int | Fraction) -> CliffordElement:
        """Return a scalar Clifford element."""

        coefficient = Fraction(value)
        return cls(()) if coefficient == 0 else cls(((0, coefficient),))

    @classmethod
    def blade(cls, *indices: int) -> CliffordElement:
        """Return ``e_i1 ... e_ik`` for increasing one-based indices."""

        if tuple(sorted(indices)) != indices or len(set(indices)) != len(indices):
            msg = "blade indices must be strictly increasing"
            raise ValueError(msg)
        if any(index < 1 or index > CLIFFORD_DIMENSION for index in indices):
            msg = "blade index outside 1..6"
            raise ValueError(msg)
        mask = sum(1 << (index - 1) for index in indices)
        return cls(((mask, Fraction(1)),))

    def __add__(self, other: CliffordElement) -> CliffordElement:
        coefficients = dict(self.terms)
        for mask, coefficient in other.terms:
            coefficients[mask] = coefficients.get(mask, Fraction(0)) + coefficient
        return CliffordElement.from_dict(coefficients)

    def __neg__(self) -> CliffordElement:
        return CliffordElement(
            tuple((mask, -coefficient) for mask, coefficient in self.terms)
        )

    def __sub__(self, other: CliffordElement) -> CliffordElement:
        return self + (-other)

    def __mul__(self, other: CliffordElement) -> CliffordElement:
        coefficients: dict[int, Fraction] = {}
        for first_mask, first_coefficient in self.terms:
            for second_mask, second_coefficient in other.terms:
                sign, output_mask = _blade_product(first_mask, second_mask)
                contribution = sign * first_coefficient * second_coefficient
                coefficients[output_mask] = (
                    coefficients.get(output_mask, Fraction(0)) + contribution
                )
        return CliffordElement.from_dict(coefficients)

    def __pow__(self, exponent: int) -> CliffordElement:
        if exponent < 0:
            msg = "use a known positive order to express inverses"
            raise ValueError(msg)
        result = CLIFFORD_ONE
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = result * base
            base = base * base
            remaining >>= 1
        return result

    def scale(self, scalar: int | Fraction) -> CliffordElement:
        """Multiply all coefficients by a rational scalar."""

        value = Fraction(scalar)
        return CliffordElement.from_dict(
            {
                mask: value * coefficient
                for mask, coefficient in self.terms
            }
        )


CLIFFORD_ZERO: Final = CliffordElement.scalar(0)
CLIFFORD_ONE: Final = CliffordElement.scalar(1)
CLIFFORD_MINUS_ONE: Final = CliffordElement.scalar(-1)


def root_difference(first: int, second: int) -> CliffordElement:
    """Return ``e_first-e_second`` without the ``sqrt(2)`` denominator."""

    if first == second:
        msg = "a root difference requires distinct indices"
        raise ValueError(msg)
    return CliffordElement.blade(first) - CliffordElement.blade(second)


def three_cycle_lift(first: int, second: int, third: int) -> CliffordElement:
    """Return the canonical order-three lift of ``(first second third)``."""

    if len({first, second, third}) != 3:
        msg = "a three-cycle requires distinct indices"
        raise ValueError(msg)
    return (
        root_difference(first, second) * root_difference(second, third)
    ).scale(Fraction(1, 2))


def three_cycle_permutation(
    first: int,
    second: int,
    third: int,
) -> Permutation:
    """Return a three-cycle in the six-point action."""

    if len({first, second, third}) != 3:
        msg = "a three-cycle requires distinct indices"
        raise ValueError(msg)
    if any(index < 1 or index > 6 for index in (first, second, third)):
        msg = "cycle point outside 1..6"
        raise ValueError(msg)
    image = list(IDENTITY)
    image[first - 1] = second - 1
    image[second - 1] = third - 1
    image[third - 1] = first - 1
    return tuple(image)


def clifford_product(elements: tuple[CliffordElement, ...]) -> CliffordElement:
    """Multiply an ordered tuple of Clifford elements."""

    result = CLIFFORD_ONE
    for element in elements:
        result = result * element
    return result


def permutation_product(elements: tuple[Permutation, ...]) -> Permutation:
    """Compose an ordered tuple with the convention ``gh=g after h``."""

    result = IDENTITY
    for element in elements:
        result = compose(result, element)
    return result


def generated_clifford_group(
    generators: tuple[CliffordElement, ...],
) -> frozenset[CliffordElement]:
    """Enumerate the finite subgroup generated by exact lift elements."""

    seen = {CLIFFORD_ONE}
    frontier = [CLIFFORD_ONE]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            candidate = element * generator
            if candidate not in seen:
                seen.add(candidate)
                frontier.append(candidate)
    return frozenset(seen)


@dataclass(frozen=True)
class A6SpinLiftCertificate:
    """Exact outputs for the forced local relations and hostile completions."""

    canonical_orders: tuple[int, int, int]
    braid_common_value: CliffordElement
    braid_relation_holds: bool
    braid_common_square: CliffordElement
    collision_commutes: bool
    collision_pair_square: CliffordElement
    lift_group_order: int
    downstairs_group_order: int
    prefix: Permutation
    prefix_lift_fourth_power: CliffordElement
    infinity_product: CliffordElement
    positive_tuple_product: Permutation
    negative_tuple_product: Permutation
    positive_spin: CliffordElement
    negative_spin: CliffordElement
    positive_tuple_generates_a6: bool
    negative_tuple_generates_a6: bool
    lifted_cusp_central_element: CliffordElement
    lifted_preferred_longitude: CliffordElement
    expected_preferred_longitude: CliffordElement

    @property
    def verified(self) -> bool:
        """Whether every exact hostile spin-lifting assertion holds."""

        expected_braid = (
            -CliffordElement.blade(1, 2)
            - CliffordElement.blade(2, 4)
            + CliffordElement.blade(1, 5)
            - CliffordElement.blade(4, 5)
        ).scale(Fraction(1, 2))
        return (
            self.canonical_orders == (3, 3, 3)
            and self.braid_relation_holds
            and self.braid_common_value == expected_braid
            and self.braid_common_square == CLIFFORD_MINUS_ONE
            and self.collision_commutes
            and self.collision_pair_square == CLIFFORD_MINUS_ONE
            and self.lift_group_order == 720
            and self.downstairs_group_order == 360
            and cycle_type(self.prefix) == (4, 2)
            and self.prefix_lift_fourth_power == CLIFFORD_MINUS_ONE
            and self.infinity_product == CLIFFORD_ONE
            and self.positive_tuple_product == IDENTITY
            and self.negative_tuple_product == IDENTITY
            and self.positive_spin == CLIFFORD_ONE
            and self.negative_spin == CLIFFORD_MINUS_ONE
            and self.positive_tuple_generates_a6
            and self.negative_tuple_generates_a6
            and self.lifted_cusp_central_element == CLIFFORD_MINUS_ONE
            and self.lifted_preferred_longitude
            == self.expected_preferred_longitude
        )


def _element_order(element: CliffordElement, *, maximum: int = 24) -> int:
    """Return the first positive power equal to one within a safe bound."""

    power = CLIFFORD_ONE
    for order in range(1, maximum + 1):
        power = power * element
        if power == CLIFFORD_ONE:
            return order
    msg = f"element order exceeds checked bound {maximum}"
    raise ValueError(msg)


def exact_a6_spin_lift_certificate(
    *,
    collision_lift_sign: int = 1,
) -> A6SpinLiftCertificate:
    """Build the exact spin certificate, optionally perturbing one lift sign."""

    if collision_lift_sign not in (-1, 1):
        msg = "the lift sign must be +1 or -1"
        raise ValueError(msg)

    r = three_cycle_permutation(3, 4, 5)
    s = three_cycle_permutation(1, 2, 3)
    b = three_cycle_permutation(4, 5, 6)
    r_lift = three_cycle_lift(3, 4, 5)
    s_lift = three_cycle_lift(1, 2, 3)
    b_lift = three_cycle_lift(4, 5, 6).scale(collision_lift_sign)

    left_braid = r_lift * s_lift * r_lift * s_lift * r_lift
    right_braid = s_lift * r_lift * s_lift * r_lift * s_lift
    collision_pair = r_lift * b_lift

    prefix = permutation_product((r, s, b))
    prefix_lift = r_lift * s_lift * b_lift
    infinity_lift = prefix_lift**7

    r_inverse = inverse(r)
    s_inverse = inverse(s)
    b_inverse = inverse(b)
    r_lift_inverse = r_lift**2
    s_lift_inverse = s_lift**2
    b_lift_inverse = b_lift**2
    c = compose(compose(r, b_inverse), r_inverse)
    c_lift = r_lift * b_lift_inverse * r_lift_inverse

    positive_tuple = (r, s, b, s_inverse, r_inverse, c)
    negative_tuple = (r, s, b, s_inverse, r, b)
    positive_lifts = (
        r_lift,
        s_lift,
        b_lift,
        s_lift_inverse,
        r_lift_inverse,
        c_lift,
    )
    negative_lifts = (
        r_lift,
        s_lift,
        b_lift,
        s_lift_inverse,
        r_lift,
        b_lift,
    )

    downstairs_group = generated_group((r, s, b))
    positive_group = generated_group(positive_tuple)
    negative_group = generated_group(negative_tuple)
    cusp_central_element = (r_lift * s_lift) ** 5
    preferred_longitude = cusp_central_element * (r_lift**2)

    return A6SpinLiftCertificate(
        canonical_orders=(
            _element_order(r_lift),
            _element_order(s_lift),
            _element_order(b_lift),
        ),
        braid_common_value=left_braid,
        braid_relation_holds=left_braid == right_braid,
        braid_common_square=left_braid**2,
        collision_commutes=s_lift * b_lift == b_lift * s_lift,
        collision_pair_square=collision_pair**2,
        lift_group_order=len(
            generated_clifford_group((r_lift, s_lift, b_lift))
        ),
        downstairs_group_order=len(downstairs_group),
        prefix=prefix,
        prefix_lift_fourth_power=prefix_lift**4,
        infinity_product=prefix_lift * infinity_lift,
        positive_tuple_product=permutation_product(positive_tuple),
        negative_tuple_product=permutation_product(negative_tuple),
        positive_spin=clifford_product(positive_lifts),
        negative_spin=clifford_product(negative_lifts),
        positive_tuple_generates_a6=(
            is_transitive(positive_group) and len(positive_group) == 360
        ),
        negative_tuple_generates_a6=(
            is_transitive(negative_group) and len(negative_group) == 360
        ),
        lifted_cusp_central_element=cusp_central_element,
        lifted_preferred_longitude=preferred_longitude,
        expected_preferred_longitude=-r_lift_inverse,
    )


def main() -> int:
    """Print the exact lift certificate and fail if any assertion breaks."""

    certificate = exact_a6_spin_lift_certificate()
    print(
        "forced local lifts:",
        {
            "orders": certificate.canonical_orders,
            "five-braid": certificate.braid_relation_holds,
            "collision commutes": certificate.collision_commutes,
            "lift group order": certificate.lift_group_order,
        },
    )
    print(
        "infinity completion:",
        {
            "prefix type": cycle_type(certificate.prefix),
            "lift fourth power": certificate.prefix_lift_fourth_power,
            "product": certificate.infinity_product,
        },
    )
    print(
        "all-3-cycle spin signs:",
        certificate.positive_spin,
        certificate.negative_spin,
    )
    print(
        "cusp longitude:",
        {
            "central element": certificate.lifted_cusp_central_element,
            "preferred longitude": certificate.lifted_preferred_longitude,
        },
    )
    print(f"spin certificate verified: {certificate.verified}")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
