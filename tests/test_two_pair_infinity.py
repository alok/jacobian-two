"""Tests for the exact two-characteristic-pair conditional exclusion."""

from fractions import Fraction

import pytest

from scripts.two_pair_infinity import (
    A6_NORMAL_INDEX,
    GENERIC_DEGREE,
    S6_FIRST_STAR,
    S6_NORMAL_INDEX,
    ConditionalTwoPairCertificate,
    s6_belyi_near_miss,
    s6_second_star_cases,
    two_pair_product_lower_bound,
)


def test_product_bound_immediately_excludes_the_a6_stratum() -> None:
    assert two_pair_product_lower_bound(A6_NORMAL_INDEX) == 8
    assert two_pair_product_lower_bound(S6_NORMAL_INDEX) == GENERIC_DEGREE
    with pytest.raises(ValueError, match="must be positive"):
        two_pair_product_lower_bound(0)


def test_s6_first_star_is_forced_exactly() -> None:
    assert S6_FIRST_STAR.verified
    assert (S6_FIRST_STAR.x, S6_FIRST_STAR.a, S6_FIRST_STAR.r1) == (3, 1, 2)
    assert (S6_FIRST_STAR.m1, S6_FIRST_STAR.n1) == (6, 1)
    assert (
        S6_FIRST_STAR.k0,
        S6_FIRST_STAR.k1,
        S6_FIRST_STAR.d_capital_1,
        S6_FIRST_STAR.m1_prime,
        S6_FIRST_STAR.k1_prime,
    ) == (3, 1, 5, 4, 2)


def test_both_nonnegative_s6_second_star_counts_fail() -> None:
    zero_arms, one_arm = s6_second_star_cases()

    assert zero_arms.k2 == 0
    assert zero_arms.n2 == 2
    assert zero_arms.m2_from_ratio == 2
    assert not zero_arms.star_equations_hold
    assert not zero_arms.survives

    assert one_arm.k2 == 1
    assert one_arm.n2 == 1
    assert one_arm.m2_from_ratio == 4
    assert one_arm.d_capital_2 == 3
    assert one_arm.star_equations_hold
    assert one_arm.forced_q_tilde == Fraction(14, 3)
    assert not one_arm.edge_integrality_holds
    assert not one_arm.survives


def test_s6_belyi_near_miss_is_not_the_obstruction() -> None:
    fixture = s6_belyi_near_miss()

    assert fixture.cycle_types == ((2, 2, 2), (5, 1), (4, 1, 1))
    assert fixture.group_order == 120
    assert fixture.verified


def test_complete_conditional_certificate() -> None:
    certificate = ConditionalTwoPairCertificate()

    assert certificate.a6_excluded
    assert certificate.s6_excluded
    assert certificate.verified
