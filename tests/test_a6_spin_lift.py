"""Exact positive and adversarial tests for the ``2.A6`` lift certificate."""

import pytest

from scripts.a6_spin_lift import (
    CLIFFORD_MINUS_ONE,
    CLIFFORD_ONE,
    CliffordElement,
    exact_a6_spin_lift_certificate,
    root_difference,
    three_cycle_lift,
    three_cycle_permutation,
)


def test_euclidean_clifford_basis_relations() -> None:
    e1 = CliffordElement.blade(1)
    e2 = CliffordElement.blade(2)

    assert e1 * e1 == CLIFFORD_ONE
    assert e2 * e2 == CLIFFORD_ONE
    assert e1 * e2 == -(e2 * e1)
    assert (e1 * e2) ** 2 == CLIFFORD_MINUS_ONE


def test_canonical_three_cycle_lift_has_order_three() -> None:
    lift = three_cycle_lift(1, 2, 3)

    assert lift**3 == CLIFFORD_ONE
    assert (-lift) ** 3 == CLIFFORD_MINUS_ONE
    with pytest.raises(ValueError, match="distinct"):
        three_cycle_lift(1, 1, 2)
    with pytest.raises(ValueError, match="distinct"):
        root_difference(2, 2)
    with pytest.raises(ValueError, match="outside"):
        three_cycle_permutation(1, 2, 7)


def test_forced_a6_relations_and_both_spin_signs_are_exact() -> None:
    certificate = exact_a6_spin_lift_certificate()

    assert certificate.canonical_orders == (3, 3, 3)
    assert certificate.braid_relation_holds
    assert certificate.braid_common_square == CLIFFORD_MINUS_ONE
    assert certificate.collision_commutes
    assert certificate.collision_pair_square == CLIFFORD_MINUS_ONE
    assert certificate.lift_group_order == 720
    assert certificate.downstairs_group_order == 360
    assert certificate.prefix_lift_fourth_power == CLIFFORD_MINUS_ONE
    assert certificate.infinity_product == CLIFFORD_ONE
    assert certificate.positive_spin == CLIFFORD_ONE
    assert certificate.negative_spin == CLIFFORD_MINUS_ONE
    assert certificate.positive_tuple_generates_a6
    assert certificate.negative_tuple_generates_a6
    assert certificate.lifted_cusp_central_element == CLIFFORD_MINUS_ONE
    assert certificate.verified


def test_noncanonical_collision_lift_breaks_the_certificate() -> None:
    certificate = exact_a6_spin_lift_certificate(collision_lift_sign=-1)

    assert certificate.canonical_orders[-1] == 6
    assert not certificate.verified


def test_invalid_lift_sign_is_rejected() -> None:
    with pytest.raises(ValueError, match=r"must be \+1 or -1"):
        exact_a6_spin_lift_certificate(collision_lift_sign=0)
