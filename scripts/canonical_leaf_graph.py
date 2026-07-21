"""Exact hostile graph family for the one-dicritical canonical corridor.

The graphs in this module are obtained by honest boundary blowups of the line
at infinity in ``P2``.  They certify that the augmented-canonical labels,
self-intersections, and determinant labels left by the joint-minimal leaf
theorem admit infinitely many abstract tree models.  They do not construct a
surface morphism, a Keller map, or compatible pullback/pushforward data.
"""

from __future__ import annotations

from dataclasses import dataclass, field


def integer_determinant(matrix: tuple[tuple[int, ...], ...]) -> int:
    """Return the exact determinant using fraction-free Bareiss elimination."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        msg = "the determinant requires a square matrix"
        raise ValueError(msg)
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    sign = 1
    previous_pivot = 1
    for column in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(column, size)
                if work[row][column] != 0
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for inner_column in range(column + 1, size):
                numerator = (
                    work[row][inner_column] * pivot
                    - work[row][column] * work[column][inner_column]
                )
                if numerator % previous_pivot != 0:
                    msg = "Bareiss division was not exact"
                    raise ArithmeticError(msg)
                work[row][inner_column] = numerator // previous_pivot
            work[row][column] = 0
        previous_pivot = pivot
    return sign * work[-1][-1]


@dataclass
class BoundaryBlowupTree:
    """Intersection tree and canonical labels under boundary blowups."""

    labels: list[int] = field(default_factory=lambda: [-2])
    self_intersections: list[int] = field(default_factory=lambda: [1])
    edges: set[tuple[int, int]] = field(default_factory=set)

    @staticmethod
    def _edge(first: int, second: int) -> tuple[int, int]:
        return (first, second) if first < second else (second, first)

    def free_blowup(self, parent: int) -> int:
        """Blow up a smooth point of one boundary component."""

        self.self_intersections[parent] -= 1
        child = len(self.labels)
        self.labels.append(self.labels[parent] + 1)
        self.self_intersections.append(-1)
        self.edges.add(self._edge(parent, child))
        return child

    def corner_blowup(self, first: int, second: int) -> int:
        """Blow up the crossing of two adjacent boundary components."""

        edge = self._edge(first, second)
        if edge not in self.edges:
            msg = "a corner blowup requires an existing boundary edge"
            raise ValueError(msg)
        self.edges.remove(edge)
        self.self_intersections[first] -= 1
        self.self_intersections[second] -= 1
        child = len(self.labels)
        self.labels.append(self.labels[first] + self.labels[second])
        self.self_intersections.append(-1)
        self.edges.add(self._edge(first, child))
        self.edges.add(self._edge(second, child))
        return child

    @property
    def intersection_matrix(self) -> tuple[tuple[int, ...], ...]:
        """Return the symmetric boundary intersection matrix."""

        size = len(self.labels)
        matrix = [[0] * size for _ in range(size)]
        for vertex, square in enumerate(self.self_intersections):
            matrix[vertex][vertex] = square
        for first, second in self.edges:
            matrix[first][second] = 1
            matrix[second][first] = 1
        return tuple(tuple(row) for row in matrix)

    @property
    def minus_intersection_determinant(self) -> int:
        """Return the determinant of the negative intersection matrix."""

        return integer_determinant(
            tuple(
                tuple(-entry for entry in row)
                for row in self.intersection_matrix
            )
        )

    def determinant_label(self, omitted_vertex: int) -> int:
        """Delete one vertex and take the negative-intersection determinant."""

        matrix = self.intersection_matrix
        minor = tuple(
            tuple(
                -entry
                for column, entry in enumerate(row)
                if column != omitted_vertex
            )
            for row_index, row in enumerate(matrix)
            if row_index != omitted_vertex
        )
        return integer_determinant(minor)


@dataclass(frozen=True)
class JointMinimalLeafGraphCertificate:
    """One member of the infinite hostile canonical/determinant graph family."""

    ramification_index: int
    zero_edge_blowups: int
    tree: BoundaryBlowupTree
    dicritical: int
    dicritical_neighbor: int

    @property
    def expected_dicritical_determinant(self) -> int:
        """Return the closed formula for the determinant label."""

        base = -(self.zero_edge_blowups + 1) * (
            self.zero_edge_blowups + 2
        )
        return base - (self.ramification_index - 1)

    @property
    def is_exact(self) -> bool:
        """Check every graph-only condition claimed by the hostile family."""

        positive_non_dicritical = (
            vertex
            for vertex, label in enumerate(self.tree.labels)
            if label > 0 and vertex != self.dicritical
        )
        return (
            self.tree.minus_intersection_determinant == -1
            and self.tree.labels[self.dicritical] == self.ramification_index
            and self.tree.self_intersections[self.dicritical] == -1
            and self.tree.labels[self.dicritical_neighbor]
            == self.ramification_index - 1
            and self.tree.self_intersections[self.dicritical_neighbor] == -2
            and self.tree.determinant_label(self.dicritical)
            == self.expected_dicritical_determinant
            and self.tree.determinant_label(self.dicritical_neighbor)
            == self.expected_dicritical_determinant + 1
            and all(
                self.tree.self_intersections[vertex] <= -2
                for vertex in positive_non_dicritical
            )
        )


def build_joint_minimal_leaf_graph(
    ramification_index: int,
    zero_edge_blowups: int,
) -> JointMinimalLeafGraphCertificate:
    """Build the ``S6`` or ``A6`` hostile graph with the requested tail length."""

    if ramification_index not in (2, 3):
        msg = "the surviving leaf index must be two or three"
        raise ValueError(msg)
    if zero_edge_blowups < 0:
        msg = "the number of zero-edge blowups must be nonnegative"
        raise ValueError(msg)

    tree = BoundaryBlowupTree()
    minus_one = tree.free_blowup(0)
    zero = tree.free_blowup(minus_one)
    nearest_one = tree.free_blowup(zero)
    for _ in range(zero_edge_blowups):
        nearest_one = tree.corner_blowup(zero, nearest_one)

    if ramification_index == 2:
        neighbor = nearest_one
        dicritical = tree.free_blowup(neighbor)
    else:
        neighbor = tree.free_blowup(nearest_one)
        dicritical = tree.free_blowup(neighbor)

    return JointMinimalLeafGraphCertificate(
        ramification_index=ramification_index,
        zero_edge_blowups=zero_edge_blowups,
        tree=tree,
        dicritical=dicritical,
        dicritical_neighbor=neighbor,
    )


def main() -> None:
    """Print the first six exact hostile determinant-label pairs."""

    for ramification_index in (2, 3):
        for zero_edge_blowups in range(6):
            certificate = build_joint_minimal_leaf_graph(
                ramification_index,
                zero_edge_blowups,
            )
            assert certificate.is_exact
            print(
                "e=",
                ramification_index,
                " n=",
                zero_edge_blowups,
                " d_E=",
                certificate.expected_dicritical_determinant,
                sep="",
            )


if __name__ == "__main__":
    main()
