"""Positive and adversarial tests for the hostile trace matrices."""

from sympy import expand

from scripts.trace_hostile_matrices import (
    P,
    Q,
    exact_a6_trace_hostile_certificate,
    exact_s6_trace_hostile_certificate,
)


def test_a6_hostile_matrix_realizes_every_quadratic_constraint() -> None:
    certificate = exact_a6_trace_hostile_certificate()

    assert certificate.presentation_determinant == certificate.branch_equation
    assert certificate.phi_determinant == certificate.branch_equation
    assert certificate.middle_determinant == -1
    assert expand(
        certificate.trace_determinant + certificate.branch_equation**2
    ) == 0
    assert certificate.generic_presentation_rank == 2
    assert certificate.cusp_presentation_rank == 1
    assert certificate.node_presentation_ranks == (1, 1)
    assert (certificate.generic_phi_rank, certificate.cusp_phi_rank) == (5, 4)
    assert (certificate.generic_trace_rank, certificate.cusp_trace_rank) == (4, 2)
    assert certificate.trace_unit_norm == 6
    assert certificate.verified


def test_a6_matrix_perturbation_breaks_the_certificate() -> None:
    certificate = exact_a6_trace_hostile_certificate(bottom_right_shift=1)

    assert certificate.presentation_determinant != certificate.branch_equation
    assert certificate.trace_unit_norm != 6
    assert not certificate.verified


def test_s6_hostile_matrix_realizes_cusps_nodes_and_unit_line() -> None:
    certificate = exact_s6_trace_hostile_certificate()

    assert certificate.presentation_determinant == certificate.branch_equation
    assert certificate.trace_determinant == certificate.branch_equation
    assert certificate.generic_rank == 5
    assert certificate.cusp_ranks == (4, 4, 4)
    assert certificate.node_ranks == (4, 4, 4)
    assert certificate.trace_unit_norm == 6
    assert certificate.verified


def test_hostile_branch_equations_are_nonconstant() -> None:
    a6 = exact_a6_trace_hostile_certificate()
    s6 = exact_s6_trace_hostile_certificate()

    assert {P, Q}.issubset(a6.branch_equation.free_symbols)
    assert {P, Q}.issubset(s6.branch_equation.free_symbols)
