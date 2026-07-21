"""Exact tests for the hostile local ``A6`` multiplication fixtures."""

from sympy import I, Matrix, expand, sqrt, zeros

from scripts.a6_local_multiplication import (
    BRANCH_EQUATION,
    COLLISION_P,
    COLLISION_RELATION,
    CYCLOTOMIC_FIVE,
    FIBER_VARIABLE,
    P,
    PARAMETER,
    Q,
    a6_presentation_matrix,
    a6_trace_data,
    a6_trace_unit,
    exact_a6_local_multiplication_certificate,
    hostile_normalization_family,
)


def test_presentation_adjugate_and_normalization_factor_are_exact() -> None:
    certificate = exact_a6_local_multiplication_certificate()
    presentation = a6_presentation_matrix()

    assert certificate.presentation_determinant == BRANCH_EQUATION
    assert (
        presentation * presentation.adjugate()
        - BRANCH_EQUATION * Matrix.eye(3)
    ).applyfunc(expand) == zeros(3)
    assert certificate.parametrized_branch_equation == 0
    assert certificate.normalization_kernel == (0, 0, 0)
    assert certificate.normalized_adjugate_factor_identity == (0,) * 9


def test_cusp_supports_the_exact_five_plus_one_trace_algebra() -> None:
    certificate = exact_a6_local_multiplication_certificate()

    assert certificate.trace_unit_norm == 6
    assert certificate.cusp_basis_determinant == -I * sqrt(5)
    assert certificate.cusp_gram_identity == (0,) * 36
    assert certificate.cusp_unit_identity == (0,) * 6
    assert certificate.cusp_middle_determinant != 0
    assert certificate.cusp_middle_trace_identity == (0,) * 36
    assert certificate.cusp_phi_kernel_columns == ((0,) * 6, (0,) * 6)


def test_collision_supports_the_exact_three_plus_three_trace_algebra() -> None:
    certificate = exact_a6_local_multiplication_certificate()

    assert certificate.collision_basis_determinant == 3 * I
    assert certificate.collision_gram_identity == (0,) * 36
    assert certificate.collision_unit_identity == (0,) * 6
    assert certificate.collision_middle_determinant != 0
    assert certificate.collision_middle_trace_identity == (0,) * 36
    assert certificate.collision_phi_kernel_columns == (
        (0,) * 6,
        (0,) * 6,
    )


def test_hostile_family_has_the_forced_special_fiber_partitions() -> None:
    certificate = exact_a6_local_multiplication_certificate()
    family = hostile_normalization_family()

    assert certificate.hostile_family_degree == 6
    assert certificate.hostile_family_leading_coefficient == 1
    assert expand(family.subs(PARAMETER, 0)) == expand(
        FIBER_VARIABLE**5 * (FIBER_VARIABLE - 1)
    )
    assert certificate.hostile_family_cusp_identity == 0
    assert certificate.hostile_family_collision_identity == 0
    assert certificate.verified


def test_collision_and_branch_symbols_are_not_silently_specialized() -> None:
    phi, middle, trace = a6_trace_data()
    unit = a6_trace_unit()

    assert {P, Q}.issubset(BRANCH_EQUATION.free_symbols)
    assert COLLISION_RELATION == COLLISION_P**2 + COLLISION_P - 1
    assert CYCLOTOMIC_FIVE.subs(PARAMETER, 1) == 5
    assert trace == phi.T * middle * phi
    assert expand((unit.T * trace * unit)[0] - 6) == 0
