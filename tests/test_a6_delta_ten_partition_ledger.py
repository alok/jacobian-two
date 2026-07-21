"""Tests for the conditional delta-ten collision-profile ledger."""

import pytest

from scripts.a6_delta_ten_generic import KAPPA
from scripts.a6_delta_ten_partition_ledger import (
    EXPECTED_ATOM_KIND_COUNTS,
    EXPECTED_ATOM_NAMES_BY_KIND,
    EXPECTED_CODIMENSION_HISTOGRAM,
    EXPECTED_CODIMENSION_ONE_NAMES,
    EXPECTED_CODIMENSION_TWO_NAMES,
    FiberKind,
    coefficient_slice_rank_certificate,
    contact_profile_is_ultrametric,
    exact_partition_ledger_certificate,
    expected_codimension_histogram,
    expected_contact_tree_codimension,
    global_profiles,
    local_atoms,
    profiles_with_expected_codimension,
)


def test_all_local_contact_tree_atoms_are_enumerated() -> None:
    atoms = local_atoms()

    assert len(atoms) == 35
    assert (
        tuple((kind, sum(atom.kind == kind for atom in atoms)) for kind in FiberKind)
        == EXPECTED_ATOM_KIND_COUNTS
    )
    assert tuple(
        atom.name for atom in atoms if atom.kind == FiberKind.TWO_BRANCH
    ) == tuple(f"C{contact}" for contact in range(2, 11))
    assert (
        tuple(
            (
                kind,
                tuple(sorted(atom.name for atom in atoms if atom.kind == kind)),
            )
            for kind in FiberKind
        )
        == EXPECTED_ATOM_NAMES_BY_KIND
    )
    assert all(atom.delta == sum(atom.contacts) for atom in atoms)
    assert all(
        atom.expected_codimension
        == expected_contact_tree_codimension(int(atom.kind), atom.contacts)
        for atom in atoms
    )


def test_all_eleven_quadruple_contact_trees_have_the_exact_invariants() -> None:
    quadruples = {
        atom.name: (
            atom.pair_length_partition,
            atom.delta,
            atom.expected_codimension,
        )
        for atom in local_atoms()
        if atom.kind == FiberKind.QUADRUPLE
    }

    assert quadruples == {
        "Q0": ((1, 1, 1, 1, 1, 1), 6, 2),
        "Q2": ((2, 1, 1, 1, 1, 1), 7, 3),
        "Q3": ((3, 1, 1, 1, 1, 1), 8, 4),
        "Q4": ((4, 1, 1, 1, 1, 1), 9, 5),
        "Q5": ((5, 1, 1, 1, 1, 1), 10, 6),
        "Q2|2": ((2, 2, 1, 1, 1, 1), 8, 4),
        "Q2|3": ((3, 2, 1, 1, 1, 1), 9, 5),
        "Q2|4": ((4, 2, 1, 1, 1, 1), 10, 6),
        "Q3|3": ((3, 3, 1, 1, 1, 1), 10, 6),
        "Q[222]": ((2, 2, 2, 1, 1, 1), 9, 4),
        "Q[223]": ((3, 2, 2, 1, 1, 1), 10, 5),
    }


def test_contact_tree_validation_rejects_an_impossible_triple() -> None:
    assert contact_profile_is_ultrametric(3, (2, 2, 3))
    assert not contact_profile_is_ultrametric(3, (1, 2, 3))
    with pytest.raises(ValueError, match="contact tree"):
        expected_contact_tree_codimension(3, (1, 2, 3))


def test_global_multisets_exhaust_the_delta_ten_candidate_ledger() -> None:
    profiles = global_profiles()

    assert len(profiles) == 145
    assert len(set(profiles)) == 145
    assert all(
        profile.non_node_delta + profile.node_count == 10 for profile in profiles
    )
    assert expected_codimension_histogram() == EXPECTED_CODIMENSION_HISTOGRAM
    assert (
        tuple(
            len(profiles_with_expected_codimension(codimension))
            for codimension in range(10)
        )
        == EXPECTED_CODIMENSION_HISTOGRAM
    )


def test_target_grouping_is_not_confused_with_the_pair_partition() -> None:
    profiles_by_name = {profile.name: profile for profile in global_profiles()}
    ten_nodes = profiles_by_name["10N"]
    ordinary_triple = profiles_by_name["T111 + 7N"]
    ordinary_quadruple = profiles_by_name["Q0 + 4N"]

    assert ten_nodes.pair_length_partition == (1,) * 10
    assert ordinary_triple.pair_length_partition == (1,) * 10
    assert ordinary_quadruple.pair_length_partition == (1,) * 10
    assert ten_nodes.expected_codimension == 0
    assert ordinary_triple.expected_codimension == 1
    assert ordinary_quadruple.expected_codimension == 2


def test_exact_expected_divisor_and_codimension_two_names() -> None:
    assert (
        tuple(profile.name for profile in profiles_with_expected_codimension(1))
        == EXPECTED_CODIMENSION_ONE_NAMES
    )
    assert (
        tuple(profile.name for profile in profiles_with_expected_codimension(2))
        == EXPECTED_CODIMENSION_TWO_NAMES
    )


EXPECTED_ENDPOINT_NAMES = tuple(
    sorted(
        (
            "C2 + C3 + T112 + N",
            "C2 + C3^2 + 2N",
            "C2 + C4 + T111 + N",
            "C2 + C5 + 3N",
            "C2 + Q2|2",
            "C2 + Q3",
            "C2 + T111 + T113",
            "C2 + T112^2",
            "C2 + T114 + 2N",
            "C2 + T223 + N",
            "C2^2 + C3 + T111",
            "C2^2 + C4 + 2N",
            "C2^2 + T113 + N",
            "C2^2 + T222",
            "C2^3 + C3 + N",
            "C2^3 + T112",
            "C2^5",
            "C3 + C4 + 3N",
            "C3 + Q2",
            "C3 + T111 + T112",
            "C3 + T113 + 2N",
            "C3 + T222 + N",
            "C3^2 + T111 + N",
            "C4 + Q0",
            "C4 + T111^2",
            "C4 + T112 + 2N",
            "C5 + T111 + 2N",
            "C6 + 4N",
            "Q2|3 + N",
            "Q4 + N",
            "Q[223]",
            "T111 + T114 + N",
            "T111 + T223",
            "T112 + T113 + N",
            "T112 + T222",
            "T115 + 3N",
            "T224 + 2N",
            "T333 + N",
        )
    )
)


def test_all_38_expected_endpoint_names_are_retained() -> None:
    endpoints = profiles_with_expected_codimension(5)

    assert len(endpoints) == 38
    assert tuple(profile.name for profile in endpoints) == EXPECTED_ENDPOINT_NAMES


def test_overdetermined_profiles_are_counted_but_not_declared_empty() -> None:
    profiles = global_profiles()

    assert sum(profile.expected_codimension > 5 for profile in profiles) == 55
    assert (
        sum(
            len(profiles_with_expected_codimension(codimension))
            for codimension in range(6, 10)
        )
        == 55
    )


def test_coefficient_slice_has_rank_four_for_every_kappa() -> None:
    first_minor, second_minor, common_divisor = coefficient_slice_rank_certificate()

    assert first_minor == 256 * KAPPA
    assert (second_minor - 128 * (3 * KAPPA**2 + 4)).expand() == 0
    assert common_divisor == 1


def test_full_partition_ledger_certificate_is_verified() -> None:
    certificate = exact_partition_ledger_certificate()

    assert certificate.atom_kind_counts == EXPECTED_ATOM_KIND_COUNTS
    assert certificate.atom_names_by_kind == EXPECTED_ATOM_NAMES_BY_KIND
    assert certificate.profile_count == 145
    assert certificate.codimension_histogram == EXPECTED_CODIMENSION_HISTOGRAM
    assert certificate.codimension_one_names == EXPECTED_CODIMENSION_ONE_NAMES
    assert certificate.codimension_two_names == EXPECTED_CODIMENSION_TWO_NAMES
    assert certificate.endpoint_count == 38
    assert certificate.overdetermined_count == 55
    assert certificate.verified
