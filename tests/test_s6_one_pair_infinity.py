"""Tests for the exact one-pair S6 infinity-group obstruction."""

import pytest

from scripts.s6_one_pair_infinity import (
    FIRST_SURVIVING_TARGET_DEGREES,
    S6_ORDER,
    S6_T5_12_WITNESS,
    SMALL_ONE_PAIR_DEGREES,
    TRANSPOSITION_TYPE,
    exact_s6_one_pair_infinity_certificate,
    exact_t5_12_witness_certificate,
    meridian_exponents,
    one_pair_degree_candidates,
    permutation_from_cycles,
    permutation_power,
    torus_s6_census,
)
from scripts.six_sheet_monodromy import IDENTITY, compose, cycle_type, generated_group


def test_width_and_one_pair_hypotheses_give_the_complete_small_table() -> None:
    candidates = one_pair_degree_candidates(11)

    assert tuple(candidate.affine_degrees for candidate in candidates) == (
        SMALL_ONE_PAIR_DEGREES
    )
    assert tuple(candidate.projective_pair for candidate in candidates) == (
        (2, 7),
        (3, 8),
        (4, 9),
        (2, 9),
        (3, 10),
        (6, 11),
        (5, 11),
        (4, 11),
        (3, 11),
        (2, 11),
    )
    assert tuple(
        candidate.affine_degrees
        for candidate in one_pair_degree_candidates(12)
        if candidate.affine_degrees[1] == 12
    ) == FIRST_SURVIVING_TARGET_DEGREES


def test_every_small_torus_group_census_excludes_an_s6_quotient() -> None:
    expected_compatible = {
        (5, 7): 0,
        (5, 8): 735,
        (5, 9): 0,
        (7, 9): 0,
        (7, 10): 15,
        (5, 11): 0,
        (6, 11): 15,
        (7, 11): 0,
        (8, 11): 15,
        (9, 11): 0,
    }
    for pair in SMALL_ONE_PAIR_DEGREES:
        census = torus_s6_census(*pair)
        u, v = census.meridian_exponents

        assert pair[1] * u + pair[0] * v == 1
        assert census.meridian_compatible_pairs == expected_compatible[pair]
        assert census.s6_pairs == 0
        assert census.verified


def test_the_only_nontrivial_small_near_miss_stays_on_five_sheets() -> None:
    census = torus_s6_census(5, 8)

    assert census.x_choices == 145
    assert census.y_choices == 256
    assert census.generated_order_histogram == ((2, 15), (120, 720))
    assert census.s6_pairs == 0


def test_degree_twelve_census_and_explicit_witness_are_exact() -> None:
    t5_12 = torus_s6_census(5, 12)
    t7_12 = torus_s6_census(7, 12)
    witness = exact_t5_12_witness_certificate()
    x, y = S6_T5_12_WITNESS
    u, v = meridian_exponents(5, 12)
    meridian = compose(permutation_power(x, u), permutation_power(y, v))

    assert t5_12.meridian_compatible_pairs == 2175
    assert t5_12.generated_order_histogram == (
        (2, 15),
        (120, 1440),
        (720, 720),
    )
    assert t5_12.s6_pairs == 720
    assert t7_12.s6_pairs == 0
    assert permutation_power(x, 5) == IDENTITY
    assert permutation_power(y, 12) == IDENTITY
    assert meridian == permutation_from_cycles((1, 2))
    assert cycle_type(meridian) == TRANSPOSITION_TYPE
    assert len(generated_group((x, y))) == S6_ORDER
    assert witness.verified


def test_complete_s6_one_pair_certificate() -> None:
    certificate = exact_s6_one_pair_infinity_certificate()

    assert certificate.all_degrees_through_eleven_excluded
    assert certificate.first_surviving_pair == (5, 12)
    assert certificate.verified


def test_invalid_degree_and_cycle_inputs_fail_closed() -> None:
    with pytest.raises(ValueError, match="positive"):
        one_pair_degree_candidates(0)
    with pytest.raises(ValueError, match="positive"):
        one_pair_degree_candidates(12, minimum_width=0)
    with pytest.raises(ValueError, match="coprime"):
        meridian_exponents(5, 10)
    with pytest.raises(ValueError, match="m < n"):
        torus_s6_census(8, 5)
    with pytest.raises(ValueError, match="outside"):
        permutation_from_cycles((1, 7))
