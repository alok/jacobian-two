"""Exact arithmetic replay of a conditional degree-six splice obstruction.

Section 2.4 of Orevkov's *Counter-examples to the Jacobian Conjecture at
Infinity* treats a one-dicritical compactification whose irreducible target
branch has exactly two characteristic pairs at infinity and which satisfies
his additional condition (*).  Under those extra hypotheses, equations
(10)--(15) exclude both one-dicritical degree-six passports.

This module checks the specialized integer argument and the realizable Belyi
passport that survives until the final edge-determinant equation.  It does not
assert condition (*) for an arbitrary Keller compactification.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Final

from scripts.six_sheet_monodromy import (
    IDENTITY,
    CycleType,
    Permutation,
    compose,
    cycle_type,
    generated_group,
    inverse,
    is_transitive,
)

GENERIC_DEGREE: Final = 6
A6_NORMAL_INDEX: Final = 3
S6_NORMAL_INDEX: Final = 2
TANGENTIAL_DEGREE: Final = 1


def two_pair_product_lower_bound(
    normal_index: int,
    *,
    tangential_degree: int = TANGENTIAL_DEGREE,
) -> int:
    """Return the minimum allowed by Orevkov's product identity.

    The identity is

    ``N = m^2*x*d1*d2^2*R1*S1``.

    Here ``x=a+n >= n+1``, ``R1>=2``, and the remaining displayed factors
    are positive integers.  The argument is conditional on the Section 2.4
    splice diagrams and hence on its two extra infinity hypotheses.
    """

    if normal_index < 1 or tangential_degree < 1:
        msg = "normal and tangential degrees must be positive"
        raise ValueError(msg)
    return tangential_degree**2 * (normal_index + 1) * 2


@dataclass(frozen=True)
class S6FirstStar:
    """Values forced before the second target-boundary star is analyzed."""

    x: int = 3
    a: int = 1
    r1: int = 2
    d1: int = 1
    d2: int = 1
    s1: int = 1
    m: int = 1
    n: int = 2
    k0: int = 3
    k1: int = 1
    d_capital_1: int = 5
    m1_prime: int = 4
    k1_prime: int = 2

    @property
    def m1(self) -> int:
        """Return ``m*x*d1*d2*R1`` from Proposition 2.5."""

        return self.m * self.x * self.d1 * self.d2 * self.r1

    @property
    def n1(self) -> int:
        """Return ``m*d2*S1`` from Proposition 2.5."""

        return self.m * self.d2 * self.s1

    @property
    def product_degree(self) -> int:
        """Return ``m1*n1``."""

        return self.m1 * self.n1

    @property
    def verified(self) -> bool:
        """Check the product, first-star, and Riemann--Hurwitz equations."""

        return (
            self.x == self.a + self.n
            and self.product_degree == GENERIC_DEGREE
            and self.k0 * self.r1 == self.m1
            and self.d1 + self.k1 * self.d_capital_1 == self.m1
            and self.m1_prime == self.k1 + self.k0
            and self.m1_prime + self.k1_prime == self.m1
        )


S6_FIRST_STAR: Final = S6FirstStar()


@dataclass(frozen=True)
class S6SecondStarCase:
    """One nonnegative arm count allowed by ``S2=2/(k2+1)``."""

    k2: int

    @property
    def m2_prime(self) -> int:
        """Return equation (13), ``m2'=k2+1``."""

        return self.k2 + 1

    @property
    def n2(self) -> Fraction:
        """Return ``n2=S2=2/m2'`` because ``m=1``."""

        return Fraction(S6_NORMAL_INDEX, self.m2_prime)

    @property
    def m2_from_ratio(self) -> Fraction:
        """Solve ``n1/m2=n2/m1'`` for ``m2``."""

        return Fraction(
            S6_FIRST_STAR.n1 * S6_FIRST_STAR.m1_prime,
            1,
        ) / self.n2

    @property
    def d_capital_2(self) -> Fraction | None:
        """Solve ``d2+k2*D2=m2`` when ``k2`` is nonzero."""

        if self.k2 == 0:
            return None
        return (
            self.m2_from_ratio - S6_FIRST_STAR.d2
        ) / self.k2

    @property
    def star_equations_hold(self) -> bool:
        """Whether equations (11)--(13) survive before the edge equation."""

        if self.n2.denominator != 1 or self.m2_from_ratio.denominator != 1:
            return False
        if self.k2 == 0:
            return self.m2_from_ratio == S6_FIRST_STAR.d2
        d_capital_2 = self.d_capital_2
        return bool(
            d_capital_2 is not None
            and d_capital_2.denominator == 1
            and d_capital_2 > 0
        )

    @property
    def forced_q_tilde(self) -> Fraction | None:
        """Return the edge determinant forced by equation (10), if reached."""

        if not self.star_equations_hold:
            return None
        d_capital_2 = self.d_capital_2
        if d_capital_2 is None:
            return None
        # With d1=d2=1, the capital and tilded determinants coincide.
        return Fraction(
            S6_FIRST_STAR.d_capital_1 * d_capital_2 - self.n2,
            S6_FIRST_STAR.x,
        )

    @property
    def edge_integrality_holds(self) -> bool:
        """Whether equation (10) gives an integral edge determinant."""

        q_tilde = self.forced_q_tilde
        return q_tilde is not None and q_tilde.denominator == 1

    @property
    def survives(self) -> bool:
        """Whether this arm count survives all specialized equations."""

        return self.star_equations_hold and self.edge_integrality_holds


def s6_second_star_cases() -> tuple[S6SecondStarCase, ...]:
    """Return every nonnegative ``k2`` for which ``2/(k2+1)`` is integral."""

    return tuple(
        S6SecondStarCase(k2)
        for k2 in range(S6_NORMAL_INDEX)
        if S6_NORMAL_INDEX % (k2 + 1) == 0
    )


@dataclass(frozen=True)
class BelyiNearMiss:
    """A transitive product-one triple for the forced ``S6`` passport."""

    first: Permutation
    second: Permutation
    third: Permutation

    @property
    def cycle_types(self) -> tuple[CycleType, CycleType, CycleType]:
        """Return the three branch partitions."""

        return (
            cycle_type(self.first),
            cycle_type(self.second),
            cycle_type(self.third),
        )

    @property
    def group_order(self) -> int:
        """Return the order of the generated transitive group."""

        return len(generated_group((self.first, self.second)))

    @property
    def verified(self) -> bool:
        """Check product one, passport, transitivity, and group order."""

        group = generated_group((self.first, self.second))
        return (
            compose(compose(self.first, self.second), self.third) == IDENTITY
            and self.cycle_types
            == ((2, 2, 2), (5, 1), (4, 1, 1))
            and is_transitive(group)
            and len(group) == 120
        )


def s6_belyi_near_miss() -> BelyiNearMiss:
    """Build the passport that reaches the final determinant contradiction."""

    first: Permutation = (1, 0, 3, 2, 5, 4)
    second: Permutation = (0, 2, 3, 4, 5, 1)
    third = inverse(compose(first, second))
    return BelyiNearMiss(first, second, third)


@dataclass(frozen=True)
class ConditionalTwoPairCertificate:
    """The complete specialized degree-six conditional exclusion."""

    @property
    def a6_excluded(self) -> bool:
        """Whether the product bound exceeds six for normal index three."""

        return two_pair_product_lower_bound(A6_NORMAL_INDEX) > GENERIC_DEGREE

    @property
    def s6_excluded(self) -> bool:
        """Whether every second-star arm count contradicts an exact equation."""

        return S6_FIRST_STAR.verified and all(
            not case.survives for case in s6_second_star_cases()
        )

    @property
    def verified(self) -> bool:
        """Whether both passports are conditionally excluded."""

        return (
            self.a6_excluded
            and self.s6_excluded
            and s6_belyi_near_miss().verified
        )


def main() -> int:
    """Print the exact conditional certificate and fail if it breaks."""

    certificate = ConditionalTwoPairCertificate()
    print(
        "A6 product lower bound:",
        two_pair_product_lower_bound(A6_NORMAL_INDEX),
    )
    print(
        "S6 forced first star:",
        {
            "x": S6_FIRST_STAR.x,
            "R1": S6_FIRST_STAR.r1,
            "m1": S6_FIRST_STAR.m1,
            "n1": S6_FIRST_STAR.n1,
            "D1": S6_FIRST_STAR.d_capital_1,
            "m1_prime": S6_FIRST_STAR.m1_prime,
        },
    )
    for case in s6_second_star_cases():
        print(
            "S6 second-star case:",
            {
                "k2": case.k2,
                "n2": case.n2,
                "m2": case.m2_from_ratio,
                "D2": case.d_capital_2,
                "Q_tilde": case.forced_q_tilde,
                "survives": case.survives,
            },
        )
    near_miss = s6_belyi_near_miss()
    print(
        "Belyi near-miss:",
        near_miss.cycle_types,
        "group order",
        near_miss.group_order,
    )
    print(f"conditional certificate verified: {certificate.verified}")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
