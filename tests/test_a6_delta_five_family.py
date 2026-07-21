"""Tests for the conditional A6 delta-five family reduction."""

from scripts.a6_delta_five_family import (
    CLEAN_RELATIONS,
    TANGENCY_RELATIONS,
    TRIPLE_RELATIONS,
    exact_a6_delta_five_clean_curve_certificate,
    exact_a6_delta_five_family_certificate,
)
from scripts.a6_one_pair_infinity import three_cycle_presentation_census


def test_clean_delta_five_curve_has_complete_genus_and_node_census() -> None:
    certificate = exact_a6_delta_five_clean_curve_certificate()

    assert certificate.finite_cusp_pair == (2, 5)
    assert certificate.node_count == 5
    assert certificate.infinity_pair == (5, 8)
    assert certificate.arithmetic_genus == 21
    assert certificate.total_delta == 21
    assert certificate.collision_discriminant == -4903
    assert certificate.pair_diagonal_resultant == 37
    assert certificate.tangency_resultant == 181411
    assert certificate.verified


def test_all_three_representative_complements_have_only_cyclic_images() -> None:
    for relations in (CLEAN_RELATIONS, TRIPLE_RELATIONS, TANGENCY_RELATIONS):
        census = three_cycle_presentation_census(relations)

        assert census.assignments == 40**3
        assert census.satisfying_assignments == 40
        assert census.generated_order_histogram == ((3, 40),)
        assert census.a6_assignments == 0


def test_family_algebra_and_codimension_one_reduction() -> None:
    certificate = exact_a6_delta_five_family_certificate()

    assert certificate.discriminant_factor_degrees == ((1, 1), (7, 1))
    assert certificate.incidence_identities == (0, 0, 0)
    assert certificate.clean_values == (1, -37, -1, 5)
    assert certificate.verified


def test_the_remaining_codimension_two_locus_is_nonempty() -> None:
    certificate = exact_a6_delta_five_family_certificate()

    assert certificate.residual_factorization_identities == (0, 0)
    assert certificate.residual_valid_products != (0, 0)
    assert certificate.residual_locus_identities == (0, 0, 0, 0, 0, 0)
