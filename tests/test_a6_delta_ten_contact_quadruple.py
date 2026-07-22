"""Tests for the conditional delta-ten C2 + Q0 + 2N profile."""

from sympy import Rational, expand

from scripts.a6_delta_ten_contact_quadruple import (
    CONTACT_QUADRUPLE_RELATIONS,
    SAMPLE_NODE_X_POLYNOMIAL,
    exact_delta_ten_contact_quadruple_certificate,
)
from scripts.a6_delta_ten_generic import S, T, X


def test_compatibility_component_is_one_smooth_rank_four_surface() -> None:
    certificate = exact_delta_ten_contact_quadruple_certificate()

    assert certificate.augmented_determinant_identity == 0
    assert certificate.compatibility_rational_factor_count == 1
    assert certificate.compatibility_sample_value == 0
    assert certificate.compatibility_sample_gradient == (4, 0, 8)
    assert certificate.compatibility_surface_geometrically_irreducible
    assert certificate.projective_component_geometry == (3, 1, 2)
    assert certificate.valid_singular_saturation == (True, 1)
    assert certificate.rank_three_saturations == (True, 1, 1)
    assert certificate.rank_two_saturations == (True, 0, 0)
    assert certificate.sample_coefficient_rank == 4
    assert certificate.sample_augmented_rank == 4
    assert certificate.sample_maximal_minors == (
        -12288,
        9728,
        -1536,
        -192,
        -32,
    )
    assert certificate.sample_valid_localizer != 0


def test_cramer_chart_has_a_codimension_three_coefficient_image() -> None:
    certificate = exact_delta_ten_contact_quadruple_certificate()

    assert certificate.generic_solution_first_residuals == (0, 0, 0, 0)
    assert certificate.generic_last_residual_identity == 0
    assert certificate.sample_cramer_solution == (7, -19, 13, -6)
    assert certificate.sample_incidence_residuals == (0, 0, 0, 0, 0)
    assert certificate.projection_tangent_vectors == (
        (0, 7, -18, 8, -1),
        (2, -4, 16, -6, 2),
    )
    assert certificate.coefficient_image_dimension == 2
    assert certificate.dominant_incidence_dimension == 2


def test_sample_has_one_reduced_ordinary_quadruple_fiber() -> None:
    certificate = exact_delta_ten_contact_quadruple_certificate()

    assert certificate.fiber_remainder == -1
    assert certificate.fiber_discriminant == -4944
    assert certificate.fiber_quotient_identity == 0
    assert certificate.slope_resultant_identity == 0
    assert certificate.slope_discriminant == -626582784
    assert certificate.fiber_pair_sum_resultant_identity == 0
    assert certificate.quadruple_x_resultant_identity == 0


def test_sample_has_one_exact_contact_two_and_two_nodes() -> None:
    certificate = exact_delta_ten_contact_quadruple_certificate()

    assert certificate.collision_identity == 0
    assert certificate.tangency_identity == 0
    assert certificate.collision_tangency_gcd == S - 1
    assert certificate.collision_jets_at_contact == (0, 0, 32, 864)
    assert certificate.component_discriminants == (56316985344, 72)
    assert certificate.component_resultants == (115344, -8, -2)
    assert all(
        all(value != 0 for value in row)
        for row in certificate.component_boundary_resultants
    )
    assert certificate.contact_chart_values == (-2, -2, -6)
    assert certificate.contact_pair_discriminant == -3
    assert certificate.contact_image_remainders == (3, -7)
    assert certificate.contact_p_derivative_resultant == 84
    assert certificate.contact_fiber_resultant == 4
    assert certificate.contact_jet_differences == (
        0,
        T / 49 - Rational(1, 98),
    )
    assert certificate.contact_jet_separation == Rational(3, 9604)
    assert certificate.node_x_resultant_identity == 0
    assert certificate.node_x_discriminant == 320000
    assert certificate.node_target_separations == (-476, 89, 1225)
    assert expand(SAMPLE_NODE_X_POLYNOMIAL.subs(X, 1)) == 89


def test_singular_scheme_and_genus_budget_are_exact() -> None:
    certificate = exact_delta_ten_contact_quadruple_certificate()

    assert certificate.cusp_image_factor == -2
    assert certificate.extra_critical_factor == -311472
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.sage_jacobian_length == 18
    assert certificate.sage_jacobian_radical_length == 5
    assert certificate.sage_jacobian_components == (
        (2, 2),
        (9, 1),
        (3, 1),
        (4, 1),
    )
    assert certificate.total_delta == certificate.arithmetic_genus == 28


def test_cyclic_complement_excludes_the_required_a6_image() -> None:
    certificate = exact_delta_ten_contact_quadruple_certificate()

    assert len(CONTACT_QUADRUPLE_RELATIONS) == 9
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.complement_census.assignments == 2_560_000
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0
    assert certificate.topology_propagation_dependencies == (
        "connected clean open",
        "proper projective Whitney-Thom triviality",
    )
    assert not certificate.split_and_deeper_boundaries_classified_here
    assert certificate.verified
