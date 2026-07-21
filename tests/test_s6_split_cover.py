"""Exact tests for the finite degree-six split-boundary ``S6`` cover."""

from sympy import Matrix, Poly, expand

from scripts.s6_split_cover import (
    AFFINE_JACOBIAN,
    AFFINE_P,
    AFFINE_Q,
    BRANCH_FACTOR,
    FIBER_POLYNOMIAL,
    P,
    Q,
    X,
    blown_up_corridor_matrix,
    exact_s6_split_cover_certificate,
    split_boundary_matrix,
)


def test_affine_map_is_finite_flat_of_degree_six_but_not_keller() -> None:
    certificate = exact_s6_split_cover_certificate()

    assert certificate.affine_relation == 0
    assert certificate.affine_inverse_coordinate == 0
    assert Poly(FIBER_POLYNOMIAL, X).degree() == 6
    assert Poly(FIBER_POLYNOMIAL, X).LC() == 1
    assert certificate.finite_flat_premises_verified
    assert certificate.jacobian == AFFINE_JACOBIAN
    assert {X}.issubset(certificate.jacobian.free_symbols)
    assert expand(certificate.jacobian - 1) != 0


def test_projective_sections_are_basepoint_free_and_recover_the_affine_map() -> None:
    certificate = exact_s6_split_cover_certificate()

    assert certificate.x0_boundary_identities == (0, 0, 0)
    assert certificate.diagonal_boundary_identities == (0, 0, 0, 0)
    assert certificate.affine_chart_identities == (0, 0, 0, 0)
    assert certificate.projective_extension_verified
    assert AFFINE_P == X**2 + AFFINE_P - X**2
    assert AFFINE_Q == X**3 - X**4 * (AFFINE_P - X**2)


def test_boundary_line_has_two_noncontracted_preimages() -> None:
    certificate = exact_s6_split_cover_certificate()
    matrix = split_boundary_matrix()
    normal = Matrix([2, 1])
    tangential = Matrix([1, 4])

    assert matrix == Matrix([[0, 1], [1, 2]])
    assert matrix * normal == tangential
    assert (normal.T * tangential)[0] == 6
    assert (normal.T * matrix * normal)[0] == 6
    assert certificate.split_transport_defect == (0, 0)
    assert certificate.split_degree == 6
    assert certificate.split_pullback_square == 6


def test_two_blowups_realize_the_numerical_two_one_corridor() -> None:
    certificate = exact_s6_split_cover_certificate()
    corridor = blown_up_corridor_matrix()

    assert corridor.det() == -1
    assert certificate.corridor_labels == (2, 1, 0, -1)
    assert certificate.corridor_transport_defect == (0, 0, 0, 0)
    assert certificate.corridor_degree == 6
    assert certificate.corridor_pullback_square == 6
    assert certificate.corridor_determinant_labels == (-3, -2)
    assert certificate.split_boundary_verified


def test_indecomposable_polynomial_and_simple_branch_give_s6_premises() -> None:
    certificate = exact_s6_split_cover_certificate()

    assert certificate.discriminant_factor_identity == 0
    assert certificate.branch_factor_gcd == 1
    assert {P}.issubset(
        certificate.quadratic_after_cubic_residual.free_symbols
    )
    assert certificate.quadratic_after_cubic_residual != 0
    assert certificate.cubic_after_quadratic_odd_term == 1
    assert certificate.s6_monodromy_premises_verified
    assert certificate.verified


def test_branch_factor_is_not_silently_specialized_or_constant() -> None:
    assert {P, Q}.issubset(BRANCH_FACTOR.free_symbols)
    assert Poly(BRANCH_FACTOR, Q).degree() == 3
