"""Exact tests for the split ``Q0`` and two-triple rank closure."""

from sympy import Rational, expand, factor

from scripts.a6_delta_ten_double_triple import (
    FIRST_OMITTED_ROOT,
    SECOND_OMITTED_ROOT,
)
from scripts.a6_delta_ten_quadruple import FIBER_VALUE
from scripts.a6_delta_ten_split_rank_q0_t111 import (
    exact_split_double_triple_rank_certificate,
    exact_split_q0_double_triple_rank_certificate,
    exact_split_quadruple_rank_certificate,
)


def test_split_quadruple_rankdrop_is_only_on_invalid_zero_fiber() -> None:
    certificate = exact_split_quadruple_rank_certificate()

    assert certificate.global_minor_identities == (0, 0, 0, 0)
    assert certificate.localized_rankdrop_bases == ((1,), (1,), (1,))
    assert certificate.zero_fiber_ranks == (2, 2, 2)
    assert certificate.zero_fiber_discriminant_values == (0, 0, 0)
    assert tuple(factor(value) for value in certificate.fiber_discriminants) == (
        -16 * FIBER_VALUE * (4 * FIBER_VALUE + 1) ** 2,
        -16 * FIBER_VALUE**2 * (16 * FIBER_VALUE - 1),
        -16 * FIBER_VALUE**2 * (16 * FIBER_VALUE - 1),
    )
    assert certificate.maximum_valid_incidence_dimension == 2
    assert certificate.verified


def test_split_two_triple_residual_is_boundary_or_incompatible() -> None:
    certificate = exact_split_double_triple_rank_certificate()
    u = FIRST_OMITTED_ROOT
    v = SECOND_OMITTED_ROOT

    assert certificate.residual_specializations == (
        -u - v,
        2 * u * v - u - v,
        -2 * u * v - u - v,
    )
    assert tuple(factor(value) for value in certificate.same_fiber_specializations) == (
        (u + v) * (u**2 + v**2 + 1),
        (u + v + 1) * (u**2 + u + v**2 + v),
        (u + v - 1) * (u**2 - u + v**2 - v),
    )
    assert certificate.k0_residual_same_fiber_identity == 0
    assert certificate.global_valid_residual_basis == (1,)
    assert certificate.maximum_valid_incidence_dimension == 2
    assert certificate.plus_minus_transport_verified
    assert certificate.verified


def test_hostile_compatible_rankdrop_points_are_rejected_by_boundaries() -> None:
    certificate = exact_split_double_triple_rank_certificate()

    # At k=0 the compatible residual family is the same-P-fiber boundary.
    assert certificate.k0_hostile_ranks == (3, 3)
    assert certificate.k0_hostile_same_fiber_value == 0

    # At k=2 this exact compatible residual point is removed only by the
    # first cusp-fiber factor.  It prevents silently dropping that localizer.
    assert certificate.k2_hostile_equation_residuals == (0, 0, 0, 0)
    assert certificate.k2_hostile_residual_value == 0
    assert certificate.k2_hostile_same_fiber_value == Rational(4, 27)
    assert certificate.k2_hostile_cusp_values == (0, Rational(16, 9))
    assert certificate.k2_hostile_ranks == (3, 3)


def test_combined_split_global_fiber_rank_certificate_is_exact() -> None:
    certificate = exact_split_q0_double_triple_rank_certificate()

    assert certificate.quadruple.verified
    assert certificate.double_triple.verified
    assert certificate.verified
    assert expand(certificate.double_triple.k0_residual_same_fiber_identity) == 0
