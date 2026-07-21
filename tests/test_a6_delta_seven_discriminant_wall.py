"""Tests for the first delta-seven discriminant strata."""

import pytest
from sympy import Rational

from scripts.a6_delta_seven_discriminant_wall import (
    COEFFICIENT_SLICE_FIRST,
    COEFFICIENT_SLICE_SECOND,
    CONTACT_THREE_PARAMETERS,
    DOUBLE_INCIDENCE_PARAMETERS,
    GENERIC_G_PARAMETERS,
    HOSTILE_CONTACT_THREE_TRIPLE,
    HOSTILE_CONTACT_TWO_TRIPLE,
    R,
    TWO_DOUBLE_SAMPLE_PARAMETERS,
    collision_partitions,
    contact_relation_census,
    contact_relation_holds,
    exact_contact_three_wall_certificate,
    exact_discriminant_incidence_certificate,
    exact_generic_g_wall_certificate,
    exact_two_double_wall_certificate,
    two_contact_a6_assignments,
)
from scripts.a6_delta_seven_generic import (
    ALPHA,
    BETA,
    CUSP_COLLISION_FACTOR,
    DELTA,
    EXTRA_CRITICAL_FACTOR,
    GAMMA,
    TRIPLE_COLLISION_FACTOR,
)
from scripts.six_sheet_monodromy import generated_group


def test_coefficient_slice_and_incidence_normalizations_are_exact() -> None:
    certificate = exact_discriminant_incidence_certificate()

    assert certificate.first_slice_identity == 0
    assert certificate.second_slice_identity == 0
    assert certificate.inverse_collision_identity == 0
    assert certificate.contact_syzygy_identity == 0
    assert certificate.double_incidence_remainders == (0, 0)
    assert certificate.double_rank_minor_gcd == R + 1
    assert certificate.triple_incidence_remainders == (0, 0, 0)
    assert certificate.triple_rank_minor_gcd == 2 * (R + 1) ** 3
    assert certificate.two_double_factor_identity == 0
    assert certificate.two_double_determinant_identity == 0
    assert certificate.verified


@pytest.mark.parametrize(
    ("builder", "parameters", "partition"),
    [
        (
            exact_generic_g_wall_certificate,
            GENERIC_G_PARAMETERS,
            (2, 1, 1, 1, 1, 1),
        ),
        (
            exact_contact_three_wall_certificate,
            CONTACT_THREE_PARAMETERS,
            (3, 1, 1, 1, 1),
        ),
        (
            exact_two_double_wall_certificate,
            TWO_DOUBLE_SAMPLE_PARAMETERS,
            (2, 2, 1, 1, 1),
        ),
    ],
)
def test_wall_sample_has_exact_geometry_and_only_cyclic_three_cycle_images(
    builder: object,
    parameters: dict[object, object],
    partition: tuple[int, ...],
) -> None:
    certificate = builder()  # type: ignore[operator]

    assert certificate.partition == partition
    assert certificate.cusp_factor == CUSP_COLLISION_FACTOR.subs(parameters)
    assert certificate.critical_factor == EXTRA_CRITICAL_FACTOR.subs(parameters)
    assert certificate.triple_factor == TRIPLE_COLLISION_FACTOR.subs(parameters)
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.verified


def test_expected_sample_validity_values_are_stable() -> None:
    generic = exact_generic_g_wall_certificate()
    contact_three = exact_contact_three_wall_certificate()
    two_double = exact_two_double_wall_certificate()

    assert (generic.cusp_factor, generic.critical_factor, generic.triple_factor) == (
        -28,
        -33125,
        961,
    )
    assert (
        contact_three.cusp_factor,
        contact_three.critical_factor,
        contact_three.triple_factor,
    ) == (Rational(-13, 4), -4100, Rational(625, 64))
    assert (
        two_double.cusp_factor,
        two_double.critical_factor,
        two_double.triple_factor,
    ) == (-3, -2750, 11)


def test_local_contact_relations_do_not_exclude_a6() -> None:
    contact_two = contact_relation_census(2)
    contact_three = contact_relation_census(3)

    assert (contact_two.satisfying_pairs, contact_two.a6_generating_triples) == (
        520,
        5760,
    )
    assert (contact_three.satisfying_pairs, contact_three.a6_generating_triples) == (
        520,
        5760,
    )
    assert len(generated_group(HOSTILE_CONTACT_TWO_TRIPLE)) == 360
    assert contact_relation_holds(*HOSTILE_CONTACT_TWO_TRIPLE[:2], 2)
    assert len(generated_group(HOSTILE_CONTACT_THREE_TRIPLE)) == 360
    assert contact_relation_holds(*HOSTILE_CONTACT_THREE_TRIPLE[:2], 3)
    assert two_contact_a6_assignments(2, 2) == 1440


def test_collision_partition_frontier_is_exhaustive() -> None:
    assert collision_partitions(2) == (
        (3, 1, 1, 1, 1),
        (2, 2, 1, 1, 1),
    )
    assert collision_partitions(3) == (
        (4, 1, 1, 1),
        (3, 2, 1, 1),
        (2, 2, 2, 1),
    )
    assert collision_partitions(4) == (
        (5, 1, 1),
        (4, 2, 1),
        (3, 3, 1),
        (3, 2, 2),
    )
    assert collision_partitions(5) == (
        (6, 1),
        (5, 2),
        (4, 3),
    )
    assert collision_partitions(6) == ((7,),)

    with pytest.raises(ValueError, match="between zero and six"):
        collision_partitions(7)


def test_a_corrupted_incidence_parameter_is_detected() -> None:
    corrupted = dict(DOUBLE_INCIDENCE_PARAMETERS)
    corrupted[ALPHA] = corrupted[ALPHA] + R

    assert corrupted != DOUBLE_INCIDENCE_PARAMETERS
    assert corrupted[ALPHA] - DOUBLE_INCIDENCE_PARAMETERS[ALPHA] == R


def test_the_coefficient_slice_is_not_vacuous() -> None:
    assert COEFFICIENT_SLICE_FIRST != 0
    assert COEFFICIENT_SLICE_SECOND != 0
    assert set(GENERIC_G_PARAMETERS) == {ALPHA, BETA, GAMMA, DELTA}
