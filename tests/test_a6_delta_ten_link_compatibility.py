"""Tests for the hostile delta-ten ``T(4,9)`` compatibility certificate."""

from scripts.a6_delta_ten_link_compatibility import (
    A6_DELTA_TEN_COLLISION,
    A6_DELTA_TEN_CUSP_FIRST,
    A6_DELTA_TEN_CUSP_SECOND,
    A6_DELTA_TEN_X,
    A6_DELTA_TEN_Y,
    exact_a6_delta_ten_link_compatibility_certificate,
    forced_cusp_collision_decompositions,
    qualifying_t49_pairs,
)
from scripts.a6_one_pair_infinity import permutation_from_cycles
from scripts.a6_spin_lift import CLIFFORD_MINUS_ONE
from scripts.six_sheet_monodromy import cycle_type, generated_group


def test_t49_pairs_form_two_inner_orbits_and_one_odd_sheet_orbit() -> None:
    certificate = exact_a6_delta_ten_link_compatibility_certificate()

    assert len(qualifying_t49_pairs()) == 720
    assert certificate.pair_passport_histogram == (
        (((4, 2), (3, 3), (3, 1, 1, 1)), 720),
    )
    assert certificate.inner_orbit_sizes == (360, 360)
    assert certificate.inner_orbits_cover_all_pairs
    assert certificate.odd_sheet_orbit_size == 720
    assert certificate.odd_sheet_fusion_pass_count == 720
    assert certificate.split_five_cycle_counts == (360, 360)
    assert certificate.cusp_torus_split_class_agreement_count == 720


def test_all_pairs_have_one_forced_cusp_collision_decomposition() -> None:
    certificate = exact_a6_delta_ten_link_compatibility_certificate()

    assert certificate.decomposition_multiplicity_histogram == ((1, 720),)
    assert certificate.unique_decomposition_count == 720
    assert certificate.cusp_pass_count == 720
    assert certificate.collision_pass_count == 720
    assert certificate.reconstruction_pass_count == 720
    assert certificate.finite_a6_pass_count == 720

    decompositions = forced_cusp_collision_decompositions(
        A6_DELTA_TEN_X,
        A6_DELTA_TEN_Y,
    )
    assert len(decompositions) == 1
    witness = decompositions[0]
    assert witness.cusp_first == A6_DELTA_TEN_CUSP_FIRST
    assert witness.cusp_second == A6_DELTA_TEN_CUSP_SECOND
    assert witness.collision == A6_DELTA_TEN_COLLISION
    assert witness.x == permutation_from_cycles((1, 2, 4, 3), (5, 6))
    assert witness.y == permutation_from_cycles((1, 3, 2), (4, 6, 5))
    assert witness.cusp_group_order == 60
    assert witness.cusp_orbit_sizes == (1, 5)
    assert witness.total_group_order == 360
    assert len(generated_group((witness.x, witness.y))) == 360
    assert witness.verified


def test_all_forced_decompositions_have_exact_spin_lifts() -> None:
    certificate = exact_a6_delta_ten_link_compatibility_certificate()

    assert certificate.local_spin_pass_count == 720
    assert certificate.torus_spin_pass_count == 720
    assert certificate.meridian_spin_pass_count == 720
    assert certificate.spin_pass_count == 720
    assert certificate.witness_spin_orders == (3, 3, 3, 8, 6)
    assert certificate.witness_spin_group_order == 720
    assert certificate.witness_infinity_longitude == CLIFFORD_MINUS_ONE
    assert certificate.witness_cusp_longitude == -(
        certificate.witness_lifted_meridian**2
    )


def test_orientation_checks_do_not_add_an_obstruction() -> None:
    certificate = exact_a6_delta_ten_link_compatibility_certificate()

    assert certificate.simultaneous_inversion_pass_count == 720
    assert certificate.simultaneous_inversion_within_inner_orbit_count == 720
    assert certificate.x_only_inversion_five_cycle_count == 720
    assert certificate.y_only_inversion_five_cycle_count == 720
    assert cycle_type(A6_DELTA_TEN_X) == (4, 2)
    assert cycle_type(A6_DELTA_TEN_Y) == (3, 3)


def test_full_delta_ten_compatibility_certificate_is_verified() -> None:
    assert exact_a6_delta_ten_link_compatibility_certificate().verified
