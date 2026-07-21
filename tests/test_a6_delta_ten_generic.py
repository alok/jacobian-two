"""Tests for the generic conditional ``A6`` delta-ten family."""

import pytest

from scripts.a6_delta_ten_generic import (
    ALPHA,
    COLLISION_POLYNOMIAL,
    KAPPA,
    S,
    ZERO_FIBER_GRAPH,
    exceptional_graph,
    exceptional_vertical,
    exact_delta_ten_algebra_certificate,
    exact_delta_ten_sample_certificate,
)
from scripts.a6_one_pair_infinity import (
    n_generator_three_cycle_presentation_census,
)


def test_complete_delta_ten_collision_algebra_is_exact() -> None:
    certificate = exact_delta_ten_algebra_certificate()

    assert certificate.collision_resultant_identity == 0
    assert certificate.collision_degree == 10
    assert certificate.collision_leading_coefficient == 1
    assert certificate.collision_constant_identity == 0
    assert certificate.residual_involution_identities == (0, 0)
    assert certificate.denominator_identity == 0
    assert certificate.cusp_image_resultant_identity == 0
    assert certificate.pair_diagonal_resultant_identity == 0
    assert certificate.tangency_syzygy_identity == 0
    assert all(value == 0 for value in certificate.exceptional_incidence_identities)
    assert certificate.verified


def test_exceptional_kappa_fibers_are_retained() -> None:
    assert (
        COLLISION_POLYNOMIAL.subs(KAPPA, 0) - S**2 * ZERO_FIBER_GRAPH
    ).expand() == 0
    for epsilon in (-1, 1):
        assert (
            COLLISION_POLYNOMIAL.subs(KAPPA, 2 * epsilon)
            - (S + epsilon) ** 4 * exceptional_graph(epsilon)
        ).expand() == 0
    assert exceptional_vertical(1) != exceptional_vertical(-1)
    assert exceptional_graph(1) != exceptional_graph(-1)
    assert exceptional_vertical(1).has(ALPHA)
    assert exceptional_graph(-1).has(ALPHA)

    with pytest.raises(ValueError, match="epsilon"):
        exceptional_vertical(0)
    with pytest.raises(ValueError, match="epsilon"):
        exceptional_graph(2)


def test_clean_delta_ten_member_has_exact_geometry_and_only_cyclic_images() -> None:
    certificate = exact_delta_ten_sample_certificate()

    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.collision_identity == 0
    assert certificate.collision_discriminant == -407351195013757923
    assert certificate.denominator_resultant == 81
    assert certificate.pair_diagonal_resultant == 335421
    assert certificate.cusp_image_factor == 1
    assert certificate.extra_critical_factor == 4141
    assert certificate.tangency_resultant == 136634145182709696290583
    assert certificate.node_x_resultant_identity == 0
    assert certificate.node_x_discriminant == -766610929107006671875
    assert certificate.total_delta == 28
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0
    assert certificate.verified


def test_n_generator_census_rejects_out_of_range_relations() -> None:
    with pytest.raises(ValueError, match="outside"):
        n_generator_three_cycle_presentation_census(((3,),), 2)
