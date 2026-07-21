"""Regression tests for the hostile one-dicritical boundary graphs."""

import pytest

from scripts.canonical_leaf_graph import (
    BoundaryBlowupTree,
    build_joint_minimal_leaf_graph,
    integer_determinant,
)


def test_exact_integer_determinant_rejects_nonsquare_input() -> None:
    assert integer_determinant(((2, 1), (1, 2))) == 3
    assert integer_determinant(()) == 1
    with pytest.raises(ValueError, match="square matrix"):
        integer_determinant(((1, 2),))


def test_corner_blowup_requires_an_edge() -> None:
    tree = BoundaryBlowupTree()
    child = tree.free_blowup(0)
    with pytest.raises(ValueError, match="existing boundary edge"):
        tree.corner_blowup(0, child + 1)


@pytest.mark.parametrize("ramification_index", [2, 3])
def test_joint_minimal_hostile_graph_family_is_exact(
    ramification_index: int,
) -> None:
    for zero_edge_blowups in range(20):
        certificate = build_joint_minimal_leaf_graph(
            ramification_index,
            zero_edge_blowups,
        )
        assert certificate.is_exact
        assert certificate.expected_dicritical_determinant == (
            -(zero_edge_blowups + 1) * (zero_edge_blowups + 2)
            - (ramification_index - 1)
        )


def test_hostile_graph_builder_rejects_parameters_outside_its_scope() -> None:
    with pytest.raises(ValueError, match="index must be two or three"):
        build_joint_minimal_leaf_graph(4, 0)
    with pytest.raises(ValueError, match="must be nonnegative"):
        build_joint_minimal_leaf_graph(2, -1)
