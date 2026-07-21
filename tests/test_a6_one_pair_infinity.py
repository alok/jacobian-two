"""Tests for the conditional A6 one-pair infinity obstruction."""

import pytest

from scripts.a6_one_pair_infinity import (
    A6_T35_WITNESS,
    A6_T38_WITNESS,
    CONTACT_TWO_RELATIONS,
    GENERIC_FAMILY_RELATIONS,
    exact_a6_one_pair_algebra_certificate,
    exact_a6_one_pair_infinity_certificate,
    meridian_exponents,
    one_pair_degree_candidates,
    permutation_power,
    three_cycle_presentation_census,
    torus_a6_census,
)
from scripts.six_sheet_monodromy import compose, cycle_type, generated_group


def test_genus_identity_reduces_the_two_smallest_collision_budgets() -> None:
    delta_one = one_pair_degree_candidates(1)
    delta_two = one_pair_degree_candidates(2)

    assert tuple(candidate.affine_degrees for candidate in delta_one) == ((2, 7),)
    assert tuple(candidate.projective_pair for candidate in delta_one) == ((5, 7),)
    assert tuple(candidate.affine_degrees for candidate in delta_two) == (
        (2, 9),
        (3, 5),
    )
    assert tuple(candidate.projective_pair for candidate in delta_two) == (
        (7, 9),
        (2, 5),
    )
    assert tuple(
        candidate.affine_degrees
        for candidate in one_pair_degree_candidates(
            1, require_singular_infinity=False
        )
    ) == ((2, 7), (3, 4))


def test_affine_link_quotients_kill_t27_and_t29_but_not_t35() -> None:
    t27 = torus_a6_census(2, 7)
    t29 = torus_a6_census(2, 9)
    t35 = torus_a6_census(3, 5)

    assert meridian_exponents(2, 7) == (1, -3)
    assert meridian_exponents(2, 9) == (1, -4)
    assert meridian_exponents(3, 5) == (-1, 2)
    assert t27.single_three_cycle_a6_pairs == 0
    assert t29.single_three_cycle_a6_pairs == 0
    assert t35.single_three_cycle_a6_pairs == 720


def test_genus_and_link_census_push_the_frontier_to_delta_five() -> None:
    expected = {
        3: ((2, 11),),
        4: ((2, 13), (3, 7)),
        5: ((2, 15), (3, 8)),
    }
    for collision_delta, degree_pairs in expected.items():
        assert tuple(
            candidate.affine_degrees
            for candidate in one_pair_degree_candidates(collision_delta)
        ) == degree_pairs

    assert torus_a6_census(2, 11).single_three_cycle_a6_pairs == 0
    assert torus_a6_census(2, 13).single_three_cycle_a6_pairs == 0
    assert torus_a6_census(3, 7).single_three_cycle_a6_pairs == 0
    assert torus_a6_census(2, 15).single_three_cycle_a6_pairs == 0
    assert torus_a6_census(3, 8).single_three_cycle_a6_pairs == 720


def test_first_two_large_link_survivors_have_exact_witnesses() -> None:
    for pair, witness in (((3, 5), A6_T35_WITNESS), ((3, 8), A6_T38_WITNESS)):
        a, d = pair
        x, y = witness
        u, v = meridian_exponents(a, d)
        meridian = compose(permutation_power(x, u), permutation_power(y, v))

        assert permutation_power(x, a) == tuple(range(6))
        assert permutation_power(y, d) == tuple(range(6))
        assert len(generated_group(witness)) == 360
        assert cycle_type(meridian) == (3, 1, 1, 1)


def test_delta_two_family_algebra_detects_every_exception() -> None:
    certificate = exact_a6_one_pair_algebra_certificate()

    assert certificate.collision_polynomial_discriminant != 0
    assert certificate.tangency_resultant != 0
    assert certificate.extra_critical_resultant != 0
    assert certificate.verified


def test_generic_and_contact_two_complements_have_no_a6_image() -> None:
    generic = three_cycle_presentation_census(GENERIC_FAMILY_RELATIONS)
    contact_two = three_cycle_presentation_census(CONTACT_TWO_RELATIONS)

    assert generic.satisfying_assignments == 40
    assert generic.generated_order_histogram == ((3, 40),)
    assert generic.a6_assignments == 0
    assert contact_two.satisfying_assignments == 760
    assert contact_two.generated_order_histogram == ((3, 40), (60, 720))
    assert contact_two.a6_assignments == 0


def test_complete_conditional_certificate() -> None:
    assert exact_a6_one_pair_infinity_certificate().verified


def test_invalid_genus_and_meridian_inputs_fail_closed() -> None:
    with pytest.raises(ValueError, match="positive"):
        one_pair_degree_candidates(0)
    with pytest.raises(ValueError, match="coprime"):
        meridian_exponents(2, 4)
