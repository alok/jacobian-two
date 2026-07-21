"""Exact and adversarial tests for split boundary determinant transport."""

import pytest
from sympy import Matrix, Rational

from scripts.split_boundary_transport import (
    BoundaryTree,
    SplitBoundaryModel,
    a6_canonical_split_certificate,
    a6_split_boundary_model,
    s6_split_boundary_model,
)


def test_boundary_blowups_rebuild_the_exact_unimodular_trees() -> None:
    a6 = a6_split_boundary_model()
    s6 = s6_split_boundary_model()

    assert a6.source_boundary.weights == (-1, -2, -2, -1)
    assert a6.source_boundary.edges == ((0, 2), (1, 2), (2, 3))
    assert a6.source_boundary.intersection_matrix.det() == -1
    assert s6.source_boundary.weights == (-1, -3, -1, -2, -1)
    assert s6.source_boundary.edges == ((0, 1), (0, 2), (1, 4), (3, 4))
    assert s6.source_boundary.intersection_matrix.det() == 1


def test_a6_split_model_satisfies_transport_hodge_and_parity() -> None:
    model = a6_split_boundary_model()

    assert model.dual_block == Matrix([[1, 2], [2, 1]])
    assert model.transport_defect == Matrix.zeros(2, 1)
    assert model.total_degree == 6
    assert model.pullback_square == 6
    assert model.hodge_remainder == Matrix(
        [[Rational(-1, 2), Rational(1, 2)],
         [Rational(1, 2), Rational(-1, 2)]]
    )
    assert model.inertia_partition == (3, 3)
    assert model.inertia_is_even
    assert model.verified


def test_s6_split_model_satisfies_transport_hodge_and_parity() -> None:
    model = s6_split_boundary_model()

    assert model.dual_block == Matrix([[-2, 2], [2, 0]])
    assert model.transport_defect == Matrix.zeros(2, 1)
    assert model.total_degree == 6
    assert model.pullback_square == 6
    assert model.hodge_remainder == Matrix(
        [[Rational(-8, 3), Rational(4, 3)],
         [Rational(4, 3), Rational(-2, 3)]]
    )
    assert model.inertia_partition == (2, 2, 2)
    assert not model.inertia_is_even
    assert model.verified


def test_a6_labeled_split_retains_canonical_and_local_monodromy_data() -> None:
    certificate = a6_canonical_split_certificate()
    source = certificate.source_model

    assert certificate.target_dual_square == 2
    assert certificate.target_boundary.labels[7] == -1
    assert certificate.target_boundary.valency(7) == 2
    assert source.dual_block == Matrix([[6, 4], [4, -2]])
    assert source.transport_defect == Matrix.zeros(2, 1)
    assert source.normal_degrees == (5, 1)
    assert source.tangential_degrees == (1, 1)
    assert source.source_boundary.labels[4] == -5
    assert source.source_boundary.labels[6] == -1
    assert source.source_boundary.valency(4) == 2
    assert source.source_boundary.valency(6) == 2
    assert certificate.canonical_label_defect == (0, 0)
    assert certificate.generated_group_order == 360
    assert certificate.verified


def test_split_transport_rejects_wrong_degree_assignment() -> None:
    valid = a6_split_boundary_model()
    perturbed = SplitBoundaryModel(
        name="wrong normal index",
        source_boundary=valid.source_boundary,
        noncontracted_components=valid.noncontracted_components,
        tangential_degrees=valid.tangential_degrees,
        normal_degrees=(3, 2),
        expected_source_inverse=valid.expected_source_inverse,
        expected_dual_block=valid.expected_dual_block,
    )

    assert perturbed.transport_defect != Matrix.zeros(2, 1)
    assert perturbed.total_degree == 5
    assert not perturbed.verified


def test_invalid_boundary_operations_and_split_vectors_are_rejected() -> None:
    boundary = BoundaryTree.positive_line()

    with pytest.raises(ValueError, match="outside"):
        boundary.smooth_blowup(1)
    with pytest.raises(ValueError, match="adjacent"):
        boundary.edge_blowup(0, 1)
    with pytest.raises(ValueError, match="at least two"):
        SplitBoundaryModel(
            name="not split",
            source_boundary=boundary,
            noncontracted_components=(0,),
            tangential_degrees=(1,),
            normal_degrees=(6,),
            expected_source_inverse=((1,),),
            expected_dual_block=((1,),),
        )
