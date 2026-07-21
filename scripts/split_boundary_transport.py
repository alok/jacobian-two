"""Exact split-preimage stopping models for boundary determinant transport.

Orevkov's scalar determinant formula assumes that a target boundary
component has only one noncontracted irreducible preimage.  With several such
components, the scalar ratio becomes a matrix identity.  This module checks
that identity and the resulting Hodge inequality, then constructs degree-six
``A6``- and ``S6``-compatible vertical boundary lattices that satisfy both.

The models are obtained by explicit blowups of a ``(+1)`` boundary line and
are interpreted as the reduced vertical divisor over the target boundary.  They
are intersection-lattice and generic-inertia fixtures only: no holomorphic
map, full Keller compactification, or polynomial counterexample is asserted.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from sympy import Expr, Matrix, Rational, expand

from scripts.six_sheet_monodromy import (
    Permutation,
    cycle_type,
    generated_group,
    is_transitive,
)

GENERIC_DEGREE: Final = 6


def _normalized_edge(first: int, second: int) -> tuple[int, int]:
    """Return an undirected edge in canonical order."""

    if first == second:
        msg = "a boundary edge needs two distinct vertices"
        raise ValueError(msg)
    return (first, second) if first < second else (second, first)


@dataclass(frozen=True, slots=True)
class BoundaryTree:
    """A weighted SNC boundary tree produced by point blowups."""

    weights: tuple[int, ...]
    edges: tuple[tuple[int, int], ...]
    labels: tuple[int, ...]

    def __post_init__(self) -> None:
        vertex_count = len(self.weights)
        if len(self.labels) != vertex_count:
            msg = "canonical labels must match the boundary vertices"
            raise ValueError(msg)
        normalized = tuple(
            sorted(_normalized_edge(first, second) for first, second in self.edges)
        )
        if normalized != self.edges or len(set(self.edges)) != len(self.edges):
            msg = "boundary edges must be unique and canonically sorted"
            raise ValueError(msg)
        if any(
            first < 0
            or second >= vertex_count
            for first, second in self.edges
        ):
            msg = "boundary edge endpoint outside the vertex set"
            raise ValueError(msg)
        if vertex_count > 0 and len(self.edges) != vertex_count - 1:
            msg = "the boundary graph must be a tree"
            raise ValueError(msg)

    @classmethod
    def positive_line(cls) -> BoundaryTree:
        """Return the initial projective boundary line of square ``+1``."""

        return cls((1,), (), (-2,))

    def smooth_blowup(self, vertex: int) -> BoundaryTree:
        """Blow up a smooth point of one boundary component."""

        if vertex < 0 or vertex >= len(self.weights):
            msg = "smooth blowup vertex outside the boundary"
            raise ValueError(msg)
        weights = list(self.weights)
        weights[vertex] -= 1
        exceptional = len(weights)
        weights.append(-1)
        labels = (*self.labels, self.labels[vertex] + 1)
        edges = tuple(
            sorted((*self.edges, _normalized_edge(vertex, exceptional)))
        )
        return BoundaryTree(tuple(weights), edges, labels)

    def edge_blowup(self, first: int, second: int) -> BoundaryTree:
        """Blow up the intersection of two adjacent boundary components."""

        edge = _normalized_edge(first, second)
        if edge not in self.edges:
            msg = "edge blowup requires adjacent boundary components"
            raise ValueError(msg)
        weights = list(self.weights)
        weights[first] -= 1
        weights[second] -= 1
        exceptional = len(weights)
        weights.append(-1)
        labels = (*self.labels, self.labels[first] + self.labels[second])
        remaining_edges = [candidate for candidate in self.edges if candidate != edge]
        remaining_edges.extend(
            (
                _normalized_edge(first, exceptional),
                _normalized_edge(second, exceptional),
            )
        )
        return BoundaryTree(
            tuple(weights),
            tuple(sorted(remaining_edges)),
            labels,
        )

    def valency(self, vertex: int) -> int:
        """Return the number of boundary neighbors of one component."""

        if vertex < 0 or vertex >= len(self.weights):
            msg = "valency vertex outside the boundary"
            raise ValueError(msg)
        return sum(vertex in edge for edge in self.edges)

    @property
    def intersection_matrix(self) -> Matrix:
        """Return the raw boundary intersection matrix."""

        matrix = Matrix.zeros(len(self.weights))
        for index, weight in enumerate(self.weights):
            matrix[index, index] = weight
        for first, second in self.edges:
            matrix[first, second] = 1
            matrix[second, first] = 1
        return matrix

    @property
    def inverse_intersection_matrix(self) -> Matrix:
        """Return the Gram matrix of the dual boundary divisor basis."""

        return self.intersection_matrix.inv()


@dataclass(frozen=True, slots=True)
class SplitBoundaryModel:
    """Several noncontracted source components over one target component.

    ``normal_degrees`` are the coefficients in the pullback of the target
    component.  ``tangential_degrees`` are the degrees of the source
    components over the target component.  ``dual_block`` is the appropriate
    block of the inverse source intersection matrix.
    """

    name: str
    source_boundary: BoundaryTree
    noncontracted_components: tuple[int, ...]
    tangential_degrees: tuple[int, ...]
    normal_degrees: tuple[int, ...]
    expected_source_inverse: tuple[tuple[int, ...], ...]
    expected_dual_block: tuple[tuple[int, ...], ...]
    target_dual_square: int = 1

    def __post_init__(self) -> None:
        rank = len(self.noncontracted_components)
        if not (
            rank == len(self.tangential_degrees)
            == len(self.normal_degrees)
        ):
            msg = "split-boundary vectors must have the same rank"
            raise ValueError(msg)
        if rank < 2:
            msg = "a split model must have at least two noncontracted components"
            raise ValueError(msg)
        if any(value <= 0 for value in (*self.tangential_degrees, *self.normal_degrees)):
            msg = "normal and tangential degrees must be positive"
            raise ValueError(msg)

    @property
    def dual_block(self) -> Matrix:
        """Return ``Q=(C_i^vee.C_j^vee)`` for the selected components."""

        inverse = self.source_boundary.inverse_intersection_matrix
        indices = list(self.noncontracted_components)
        return inverse.extract(indices, indices)

    @property
    def tangential_vector(self) -> Matrix:
        """Return the column vector ``m`` of tangential degrees."""

        return Matrix(self.tangential_degrees)

    @property
    def normal_vector(self) -> Matrix:
        """Return the column vector ``n`` of normal degrees."""

        return Matrix(self.normal_degrees)

    @property
    def total_degree(self) -> int:
        """Return ``n.m``, the generic degree over the target boundary."""

        return sum(
            normal * tangential
            for normal, tangential in zip(
                self.normal_degrees,
                self.tangential_degrees,
                strict=True,
            )
        )

    @property
    def transport_defect(self) -> Matrix:
        """Return the defect in the split transport identity ``Qm=q*n``."""

        return (
            self.dual_block * self.tangential_vector
            - self.target_dual_square * self.normal_vector
        ).applyfunc(expand)

    @property
    def pullback_square(self) -> Expr:
        """Return ``m^t Q m`` as an exact scalar."""

        value = (self.tangential_vector.T * self.dual_block * self.tangential_vector)[0]
        return expand(value)

    @property
    def hodge_remainder(self) -> Matrix:
        """Return ``Q-(q/N)nn^t``, which Hodge makes negative semidefinite."""

        normal = self.normal_vector
        return (
            self.dual_block
            - Rational(self.target_dual_square, self.total_degree)
            * normal
            * normal.T
        ).applyfunc(expand)

    @property
    def hodge_remainder_is_negative_semidefinite(self) -> bool:
        """Check negative semidefiniteness for the exact rank-two fixtures."""

        remainder = self.hodge_remainder
        if remainder.shape != (2, 2) or remainder != remainder.T:
            return False
        return bool(
            remainder[0, 0] <= 0
            and remainder[1, 1] <= 0
            and expand(remainder.det()) >= 0
        )

    @property
    def inertia_partition(self) -> tuple[int, ...]:
        """Expand normal indices according to the tangential degrees."""

        return tuple(
            sorted(
                normal
                for normal, tangential in zip(
                    self.normal_degrees,
                    self.tangential_degrees,
                    strict=True,
                )
                for _ in range(tangential)
            )
        )

    @property
    def inertia_is_even(self) -> bool:
        """Return the sign parity of the generic boundary permutation."""

        defect = sum(value - 1 for value in self.inertia_partition)
        return defect % 2 == 0

    @property
    def verified(self) -> bool:
        """Whether the blowup lattice and every transport identity are exact."""

        inverse = self.source_boundary.inverse_intersection_matrix
        expected_inverse = Matrix(self.expected_source_inverse)
        expected_block = Matrix(self.expected_dual_block)
        return bool(
            self.source_boundary.intersection_matrix.det() in (-1, 1)
            and inverse == expected_inverse
            and self.dual_block == expected_block
            and self.total_degree == GENERIC_DEGREE
            and self.transport_defect == Matrix.zeros(len(self.normal_degrees), 1)
            and self.pullback_square
            == self.total_degree * self.target_dual_square
            and self.hodge_remainder * self.tangential_vector
            == Matrix.zeros(len(self.normal_degrees), 1)
            and self.hodge_remainder_is_negative_semidefinite
        )


def a6_split_boundary_model() -> SplitBoundaryModel:
    """Return a two-component ``3+3`` split model after three blowups."""

    boundary = BoundaryTree.positive_line()
    boundary = boundary.smooth_blowup(0)
    boundary = boundary.edge_blowup(0, 1)
    boundary = boundary.smooth_blowup(2)
    return SplitBoundaryModel(
        name="A6 3+3 split",
        source_boundary=boundary,
        noncontracted_components=(0, 3),
        tangential_degrees=(1, 1),
        normal_degrees=(3, 3),
        expected_source_inverse=(
            (1, 1, 2, 2),
            (1, 0, 1, 1),
            (2, 1, 2, 2),
            (2, 1, 2, 1),
        ),
        expected_dual_block=((1, 2), (2, 1)),
    )


def s6_split_boundary_model() -> SplitBoundaryModel:
    """Return a three-transposition split model after four blowups."""

    boundary = BoundaryTree.positive_line()
    boundary = boundary.smooth_blowup(0)
    boundary = boundary.smooth_blowup(0)
    boundary = boundary.smooth_blowup(1)
    boundary = boundary.edge_blowup(1, 3)
    return SplitBoundaryModel(
        name="S6 2+2+2 split",
        source_boundary=boundary,
        noncontracted_components=(4, 2),
        tangential_degrees=(1, 2),
        normal_degrees=(2, 2),
        expected_source_inverse=(
            (1, 1, 1, 1, 2),
            (1, 0, 1, 0, 0),
            (1, 1, 0, 1, 2),
            (1, 0, 1, -1, -1),
            (2, 0, 2, -1, -2),
        ),
        expected_dual_block=((-2, 2), (2, 0)),
    )


def _a6_labeled_target_boundary() -> BoundaryTree:
    """Build the target tree for the canonical-label split fixture."""

    boundary = BoundaryTree.positive_line()
    boundary = boundary.smooth_blowup(0)
    boundary = boundary.edge_blowup(0, 1)
    boundary = boundary.edge_blowup(0, 2)
    boundary = boundary.smooth_blowup(3)
    boundary = boundary.smooth_blowup(4)
    boundary = boundary.smooth_blowup(5)
    boundary = boundary.smooth_blowup(6)
    return boundary.smooth_blowup(7)


def _a6_labeled_source_model() -> SplitBoundaryModel:
    """Build the source tree with split normal indices five and one."""

    boundary = BoundaryTree.positive_line()
    boundary = boundary.smooth_blowup(0)
    boundary = boundary.edge_blowup(0, 1)
    boundary = boundary.smooth_blowup(2)
    boundary = boundary.edge_blowup(2, 3)
    boundary = boundary.smooth_blowup(1)
    boundary = boundary.edge_blowup(1, 5)
    return SplitBoundaryModel(
        name="A6 labeled 5+1 split",
        source_boundary=boundary,
        noncontracted_components=(4, 6),
        tangential_degrees=(1, 1),
        normal_degrees=(5, 1),
        expected_source_inverse=(
            (1, 1, 2, 2, 4, 1, 2),
            (1, 0, 1, 1, 2, 0, 0),
            (2, 1, 2, 2, 4, 1, 2),
            (2, 1, 2, 1, 3, 1, 2),
            (4, 2, 4, 3, 6, 2, 4),
            (1, 0, 1, 1, 2, -1, -1),
            (2, 0, 2, 2, 4, -1, -2),
        ),
        expected_dual_block=((6, 4), (4, -2)),
        target_dual_square=2,
    )


@dataclass(frozen=True, slots=True)
class A6CanonicalSplitCertificate:
    """A split fixture retaining labels, local normal forms, and monodromy."""

    target_boundary: BoundaryTree
    target_component: int
    source_model: SplitBoundaryModel
    infinity_inertia: Permutation
    finite_three_cycle: Permutation

    @property
    def target_dual_square(self) -> Expr:
        """Return the selected target dual divisor square."""

        inverse = self.target_boundary.inverse_intersection_matrix
        return expand(inverse[self.target_component, self.target_component])

    @property
    def canonical_label_defect(self) -> tuple[int, ...]:
        """Check ``label(C_i)=n_i*label(T)`` for type-one local maps."""

        target_label = self.target_boundary.labels[self.target_component]
        return tuple(
            self.source_model.source_boundary.labels[component]
            - normal_degree * target_label
            for component, normal_degree in zip(
                self.source_model.noncontracted_components,
                self.source_model.normal_degrees,
                strict=True,
            )
        )

    @property
    def generated_group_order(self) -> int:
        """Return the order generated by the 5-cycle and a 3-cycle."""

        return len(generated_group((self.infinity_inertia, self.finite_three_cycle)))

    @property
    def verified(self) -> bool:
        """Whether every labeled local stopping assertion is exact."""

        target_inverse = self.target_boundary.inverse_intersection_matrix
        group = generated_group((self.infinity_inertia, self.finite_three_cycle))
        source = self.source_model.source_boundary
        return bool(
            self.target_boundary.weights
            == (-2, -2, -2, -2, -2, -2, -2, -2, -1)
            and self.target_boundary.labels
            == (-2, -1, -3, -5, -4, -3, -2, -1, 0)
            and self.target_boundary.intersection_matrix.det() == 1
            and target_inverse[self.target_component, self.target_component] == 2
            and self.target_boundary.valency(self.target_component) == 2
            and source.weights == (-1, -4, -3, -2, -1, -2, -1)
            and source.labels == (-2, -1, -3, -2, -5, 0, -1)
            and source.valency(4) == 2
            and source.valency(6) == 2
            and self.source_model.verified
            and self.canonical_label_defect == (0, 0)
            and cycle_type(self.infinity_inertia) == (5, 1)
            and cycle_type(self.finite_three_cycle) == (3, 1, 1, 1)
            and is_transitive(group)
            and len(group) == 360
        )


def a6_canonical_split_certificate() -> A6CanonicalSplitCertificate:
    """Return the exact labeled ``5+1`` split boundary certificate."""

    five_cycle: Permutation = (1, 2, 3, 4, 0, 5)
    three_cycle: Permutation = (0, 1, 2, 4, 5, 3)
    return A6CanonicalSplitCertificate(
        target_boundary=_a6_labeled_target_boundary(),
        target_component=7,
        source_model=_a6_labeled_source_model(),
        infinity_inertia=five_cycle,
        finite_three_cycle=three_cycle,
    )


def main() -> int:
    """Print both exact split fixtures and fail on a regression."""

    models = (a6_split_boundary_model(), s6_split_boundary_model())
    for model in models:
        print(
            model.name,
            {
                "weights": model.source_boundary.weights,
                "edges": model.source_boundary.edges,
                "Q": model.dual_block.tolist(),
                "m": model.tangential_degrees,
                "n": model.normal_degrees,
                "inertia": model.inertia_partition,
                "Hodge remainder": model.hodge_remainder.tolist(),
                "verified": model.verified,
            },
        )
    labeled = a6_canonical_split_certificate()
    print(
        "A6 labeled split",
        {
            "target q": labeled.target_dual_square,
            "target label": labeled.target_boundary.labels[labeled.target_component],
            "source labels": tuple(
                labeled.source_model.source_boundary.labels[index]
                for index in labeled.source_model.noncontracted_components
            ),
            "Q": labeled.source_model.dual_block.tolist(),
            "m": labeled.source_model.tangential_degrees,
            "n": labeled.source_model.normal_degrees,
            "group order": labeled.generated_group_order,
            "verified": labeled.verified,
        },
    )
    verified = all(model.verified for model in models) and labeled.verified
    print(f"split boundary certificate verified: {verified}")
    print("claim boundary: intersection lattices only; no compactified map")
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
