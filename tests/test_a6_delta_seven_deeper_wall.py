"""Tests for the three one-dimensional delta-seven collision walls."""

import pytest
from sympy import Rational, cancel

from scripts.a6_delta_seven_deeper_wall import (
    CONTACT_FOUR_PARAMETERS,
    CONTACT_FOUR_SAMPLE_PARAMETERS,
    E,
    MH,
    MIXED_PARAMETER_CURVE,
    MIXED_SAMPLE_PARAMETERS,
    MIXED_SAMPLE_SUBSTITUTION,
    MK,
    THREE_DOUBLE_SAMPLE_PARAMETERS,
    THREE_DOUBLE_PARAMETERS,
    V,
    exact_contact_four_sample_certificate,
    exact_excess_three_algebra_certificate,
    exact_mixed_sample_certificate,
    exact_three_double_sample_certificate,
)
from scripts.a6_delta_seven_generic import (
    CUSP_COLLISION_FACTOR,
    EXTRA_CRITICAL_FACTOR,
    TRIPLE_COLLISION_FACTOR,
)


def test_all_three_excess_three_normalizations_are_exact() -> None:
    certificate = exact_excess_three_algebra_certificate()

    assert certificate.verified
    assert certificate.contact_four_factor_identity == 0
    assert certificate.mixed_first_slice_identity == 0
    assert certificate.mixed_second_slice_identity == 0
    assert certificate.three_double_factor_identity == 0


@pytest.mark.parametrize(
    ("builder", "partition", "parameters"),
    [
        (
            exact_contact_four_sample_certificate,
            (4, 1, 1, 1),
            CONTACT_FOUR_SAMPLE_PARAMETERS,
        ),
        (
            exact_mixed_sample_certificate,
            (3, 2, 1, 1),
            MIXED_SAMPLE_PARAMETERS,
        ),
        (
            exact_three_double_sample_certificate,
            (2, 2, 2, 1),
            THREE_DOUBLE_SAMPLE_PARAMETERS,
        ),
    ],
)
def test_representative_geometry_and_three_cycle_replay(
    builder: object,
    partition: tuple[int, ...],
    parameters: dict[object, object],
) -> None:
    certificate = builder()  # type: ignore[operator]

    assert certificate.partition == partition
    assert certificate.cusp_factor == CUSP_COLLISION_FACTOR.subs(parameters)
    assert certificate.critical_factor == EXTRA_CRITICAL_FACTOR.subs(parameters)
    assert certificate.triple_factor == TRIPLE_COLLISION_FACTOR.subs(parameters)
    assert certificate.total_delta == 36
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0
    assert certificate.verified


def test_validity_values_for_all_three_samples_are_stable() -> None:
    contact_four = exact_contact_four_sample_certificate()
    mixed = exact_mixed_sample_certificate()
    three_double = exact_three_double_sample_certificate()

    assert (
        contact_four.cusp_factor,
        contact_four.critical_factor,
        contact_four.triple_factor,
    ) == (-52, -48020, -38416)
    assert (mixed.cusp_factor, mixed.critical_factor, mixed.triple_factor) == (
        Rational(-16, 19683),
        Rational(-40, 27),
        Rational(64, 43046721),
    )
    assert (
        three_double.cusp_factor,
        three_double.critical_factor,
        three_double.triple_factor,
    ) == (Rational(-1, 8), Rational(-625, 8), Rational(-1, 4096))


def test_mixed_parameter_space_is_not_misreported_as_rational() -> None:
    assert MIXED_PARAMETER_CURVE.subs(
        {MH: MIXED_SAMPLE_SUBSTITUTION[MH], MK: MIXED_SAMPLE_SUBSTITUTION[MK]}
    ) == 0
    assert MIXED_PARAMETER_CURVE.as_poly(MH, MK).total_degree() == 5

    # Sage independently computes geometric genus two for the projective
    # closure.  This hostile fixture keeps the executable model as a genuine
    # bivariate quintic rather than silently replacing it by a line chart.
    assert MIXED_PARAMETER_CURVE.coeff(MH, 4).has(MK)
    assert MIXED_PARAMETER_CURVE.coeff(MK, 4).has(MH)


def test_normalization_parameters_recover_the_displayed_samples() -> None:
    assert all(
        cancel(CONTACT_FOUR_PARAMETERS[key].subs(V, 1) - expected) == 0
        for key, expected in CONTACT_FOUR_SAMPLE_PARAMETERS.items()
    )
    assert all(
        cancel(THREE_DOUBLE_PARAMETERS[key].subs(E, Rational(-9, 4)) - expected)
        == 0
        for key, expected in THREE_DOUBLE_SAMPLE_PARAMETERS.items()
    )
