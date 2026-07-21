"""Tests for the global obstruction to the explicit hostile S6 trace curve."""

import pytest
from sympy import expand

from scripts.s6_trace_curve_topology import (
    BRANCH_EQUATION,
    PARAMETER,
    PARAMETER_P,
    PARAMETER_Q,
    P,
    Q,
    VAN_KAMPEN_RELATIONS,
    all_transpositions,
    assignment_satisfies_relations,
    evaluate_signed_word,
    exact_s6_trace_curve_topology_certificate,
    minimum_transposition_generators,
    transposition_assignment_census,
)
from scripts.six_sheet_monodromy import IDENTITY


def test_trace_curve_parametrization_and_width_obstruction() -> None:
    certificate = exact_s6_trace_curve_topology_certificate()

    assert expand(
        BRANCH_EQUATION.subs(
            {P: PARAMETER_P, Q: PARAMETER_Q}
        )
    ) == 0
    assert PARAMETER in (PARAMETER_P + PARAMETER_Q).free_symbols
    assert certificate.projection_width == 4
    assert certificate.minimum_required_width == 5
    assert certificate.width_obstructs_s6
    assert certificate.verified


def test_four_meridian_transposition_census_is_exact() -> None:
    census = transposition_assignment_census()

    assert len(VAN_KAMPEN_RELATIONS) == 6
    assert census.total_assignments == 15**4
    assert census.relation_satisfying_assignments == 735
    assert dict(census.group_order_histogram) == {2: 15, 120: 720}
    assert census.transitive_assignments == 0
    assert census.verified


def test_signed_word_evaluation_and_relation_checker() -> None:
    transpositions = all_transpositions()
    images = (
        transpositions[0],
        transpositions[1],
        transpositions[2],
        transpositions[3],
    )

    assert evaluate_signed_word((1, -1, 2, -2), images) == IDENTITY
    assert assignment_satisfies_relations((transpositions[0],) * 4)


def test_transposition_edge_bound_rejects_invalid_sheet_counts() -> None:
    assert minimum_transposition_generators(1) == 0
    assert minimum_transposition_generators(6) == 5
    with pytest.raises(ValueError, match="positive"):
        minimum_transposition_generators(0)
