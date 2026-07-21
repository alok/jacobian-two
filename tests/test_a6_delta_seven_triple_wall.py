"""Tests for the conditional delta-seven triple-fiber wall."""

import pytest
from sympy import expand

from scripts.a6_delta_seven_triple_wall import (
    ALLOWED_TRIPLE_CONTACT_PROFILES,
    DOUBLE_TRIPLE_COLLISION_FACTOR,
    DOUBLE_TRIPLE_G_PULLBACK,
    DOUBLE_TRIPLE_RELATIONS,
    DOUBLE_TRIPLE_SUBSTITUTION,
    HOSTILE_EXTENSION,
    HOSTILE_LOCAL_TRIPLE,
    TRIPLE_WALL_GENERIC_RELATIONS,
    T_WALL_PARAMETERS,
    exact_triple_wall_algebra_certificate,
    exact_triple_wall_local_certificate,
    hurwitz_action,
    triple_contact_braid,
    triple_contact_census,
)
from scripts.a6_delta_seven_generic import (
    COLLISION_POLYNOMIAL,
    CUSP_COLLISION_FACTOR,
    TRIPLE_COLLISION_FACTOR,
)
from scripts.a6_one_pair_infinity import three_cycle_presentation_census
from scripts.six_sheet_monodromy import generated_group


def test_triple_wall_normalization_and_discriminants_are_exact() -> None:
    certificate = exact_triple_wall_algebra_certificate()

    assert certificate.verified
    assert expand(TRIPLE_COLLISION_FACTOR.subs(T_WALL_PARAMETERS)) == 0
    assert certificate.discriminant_wall_identity == 0


def test_double_triple_locus_has_the_stated_exact_factorization() -> None:
    wall_h = COLLISION_POLYNOMIAL.subs(T_WALL_PARAMETERS)
    wall_c = CUSP_COLLISION_FACTOR.subs(T_WALL_PARAMETERS)

    assert expand(
        wall_h.subs(DOUBLE_TRIPLE_SUBSTITUTION)
        - DOUBLE_TRIPLE_COLLISION_FACTOR
    ) == 0
    assert expand(
        exact_triple_wall_algebra_certificate().double_triple_discriminant_wall_identity
    ) == 0
    assert DOUBLE_TRIPLE_G_PULLBACK != 0
    assert wall_c.subs(DOUBLE_TRIPLE_SUBSTITUTION) != 0


def test_every_delta_seven_three_branch_contact_profile_has_no_a6_image() -> None:
    certificate = exact_triple_wall_local_certificate()

    assert tuple(
        (census.profile[0], census.profile[2])
        for census in certificate.profiles
    ) == ALLOWED_TRIPLE_CONTACT_PROFILES
    assert all(census.delta <= 7 for census in certificate.profiles)
    assert all(census.assignments == 40**3 for census in certificate.profiles)
    assert all(census.a6_assignments == 0 for census in certificate.profiles)
    assert certificate.verified


def test_ordinary_triple_and_deeper_profile_histograms_are_stable() -> None:
    ordinary = triple_contact_census(1, 1)
    deepest = triple_contact_census(2, 2)

    assert ordinary.fixed_assignments == 1000
    assert ordinary.generated_order_histogram == ((3, 160), (9, 480), (12, 360))
    assert deepest.fixed_assignments == 11_080
    assert deepest.generated_order_histogram == (
        (3, 160),
        (9, 480),
        (12, 1800),
        (60, 8640),
    )


def test_complete_fiber_hypothesis_is_not_optional() -> None:
    assert (
        hurwitz_action(HOSTILE_LOCAL_TRIPLE, triple_contact_braid(1, 1))
        == HOSTILE_LOCAL_TRIPLE
    )
    assert len(generated_group(HOSTILE_LOCAL_TRIPLE)) == 9
    assert len(generated_group((*HOSTILE_LOCAL_TRIPLE, HOSTILE_EXTENSION))) == 360


def test_two_global_sage_samples_agree_with_the_local_obstruction() -> None:
    generic = three_cycle_presentation_census(TRIPLE_WALL_GENERIC_RELATIONS)
    double = three_cycle_presentation_census(DOUBLE_TRIPLE_RELATIONS)

    assert generic.generated_order_histogram == ((3, 40),)
    assert double.generated_order_histogram == ((3, 40),)
    assert generic.a6_assignments == double.a6_assignments == 0


def test_invalid_hurwitz_profiles_are_rejected() -> None:
    with pytest.raises(ValueError, match="1 <= q <= k"):
        triple_contact_braid(0, 1)
    with pytest.raises(ValueError, match="1 <= q <= k"):
        triple_contact_braid(2, 1)
