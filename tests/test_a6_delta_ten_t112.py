"""Tests for the exact conditional delta-ten ``T112 + 6N`` member."""

from sympy import Rational

from scripts.a6_delta_ten_t112 import (
    FOURTH_TARGET_SEPARATION,
    ROOT_FACTORIZATION_IDENTITY,
    SAGE_BASE_ABSOLUTELY_IRREDUCIBLE,
    SAGE_BASE_SMOOTH_SATURATION,
    SAGE_RANK_DROP_SATURATION,
    SAGE_SAMPLE_PROJECTION_FIBER,
    SAMPLE_PROJECTION_FIBER_BASIS,
    SLICE_ROOT_SUBSTITUTION,
    SLICE_BOUNDARY_FACTORS,
    SLICE_INCIDENCE_IDENTITIES,
    T112_LABELED_TANGENT_PAIRS,
    T112_PARAMETERS,
    T112_RELATIONS,
    exact_delta_ten_t112_certificate,
)


def test_general_labeled_incidence_component_is_globally_irreducible() -> None:
    """The valid base is smooth and the Q-system has rank three everywhere."""

    certificate = exact_delta_ten_t112_certificate()

    assert ROOT_FACTORIZATION_IDENTITY == 0
    assert certificate.root_constraint_value == 0
    assert certificate.root_constraint_gradient == (
        Rational(9, 125),
        Rational(28, 125),
        Rational(133, 125),
    )
    assert certificate.fourth_root_value == Rational(-6, 5)
    assert certificate.kappa_value == 2
    assert certificate.base_constraint_irreducible
    assert (
        certificate.base_constraint_absolutely_irreducible
        == SAGE_BASE_ABSOLUTELY_IRREDUCIBLE
    )
    assert certificate.base_surface_dimension == 2
    assert certificate.valid_base_value != 0
    assert certificate.base_smooth_saturation == SAGE_BASE_SMOOTH_SATURATION
    assert certificate.rank_drop_saturation == SAGE_RANK_DROP_SATURATION
    assert certificate.incidence_rank == 3
    assert certificate.incidence_affine_fiber_dimension == 1
    assert certificate.labeled_component_dimension == 3
    assert certificate.labeled_component_irreducible
    assert T112_LABELED_TANGENT_PAIRS == (("u", "v"), ("u", "w"), ("v", "w"))


def test_projection_is_generically_the_two_orientation_label_cover() -> None:
    """The clean sample fiber proves generic finiteness and degree two."""

    certificate = exact_delta_ten_t112_certificate()

    assert certificate.sample_projection_fiber == SAGE_SAMPLE_PROJECTION_FIBER
    assert certificate.sample_projection_fiber_basis == SAMPLE_PROJECTION_FIBER_BASIS
    assert certificate.coefficient_image_dimension == 3
    assert certificate.generic_label_degree == 2
    assert certificate.unlabeled_component_irreducible
    assert certificate.clean_open_connected
    assert certificate.proper_whitney_thom_required


def test_rational_slice_solves_incidence_and_avoids_named_boundaries() -> None:
    """The one-parameter formulas solve all equations before specialization."""

    certificate = exact_delta_ten_t112_certificate()

    assert SLICE_INCIDENCE_IDENTITIES == (0, 0, 0)
    assert certificate.vertical_factor_identity == 0
    assert certificate.graph_factor_identity == 0
    assert len(SLICE_BOUNDARY_FACTORS) == 12
    assert tuple(name for name, _ in certificate.boundary_values) == tuple(
        name for name, _ in SLICE_BOUNDARY_FACTORS
    )
    assert all(value != 0 for _, value in certificate.boundary_values)


def test_three_sources_give_exact_t112_and_not_a_fourth_branch() -> None:
    """Two slopes agree to first order only and the third is transverse."""

    certificate = exact_delta_ten_t112_certificate()

    assert certificate.triple_p_values == (Rational(36, 625),) * 3
    assert certificate.triple_q_values == (Rational(7776, 37109375),) * 3
    assert certificate.fourth_q_separation == Rational(24712128, 37109375)
    assert (
        FOURTH_TARGET_SEPARATION.subs(SLICE_ROOT_SUBSTITUTION).subs(T112_PARAMETERS)
        == certificate.fourth_q_separation
    )
    assert certificate.tangent_slopes == (
        Rational(1188, 59375),
        Rational(1188, 59375),
        Rational(3924, 415625),
    )
    assert certificate.tangent_second_derivatives == (
        Rational(441, 95),
        Rational(351, 190),
    )
    assert certificate.contact_second_derivative_difference == Rational(531, 190)
    assert certificate.local_terms_below_three == 0
    assert certificate.tangent_cone_identity == 0


def test_affine_singularities_exhaust_the_delta_ten_genus_budget() -> None:
    """The exact eliminants and Sage metadata give T112 plus six nodes."""

    certificate = exact_delta_ten_t112_certificate()

    assert certificate.cusp_leading_coefficient == Rational(2089, 11875)
    assert certificate.other_cusp_fiber_q_value == Rational(1296, 11875)
    assert certificate.cusp_image_factor != 0
    assert certificate.extra_critical_factor != 0
    assert certificate.vertical_residual_discriminant == 3060
    assert certificate.graph_residual_discriminant != 0
    assert certificate.vertical_node_x_identity == 0
    assert certificate.graph_node_x_identity == 0
    assert certificate.vertical_node_x_discriminant != 0
    assert certificate.graph_node_x_discriminant != 0
    assert certificate.node_packet_resultant != 0
    assert all(value != 0 for value in certificate.cusp_node_separations)
    assert all(value != 0 for value in certificate.triple_node_separations)
    assert certificate.implicit_resultant_identity == 0
    assert certificate.implicit_parameterization_identity == 0
    assert certificate.implicit_content == 1
    assert certificate.implicit_total_degree == 9
    assert certificate.sage_jacobian_components == (
        (4, 1, False),
        (6, 1, False),
        (2, 2, True),
        (4, 4, True),
    )
    assert certificate.total_delta == certificate.arithmetic_genus == 28


def test_raw_t112_presentation_has_only_diagonal_c3_images() -> None:
    """Exhaust the raw presentation rather than trusting its cyclic reduction."""

    certificate = exact_delta_ten_t112_certificate()

    assert len(T112_RELATIONS) == 11
    assert certificate.sage_cyclic_simplification == (1, 0, True)
    assert certificate.pruning_stage_survivors == (40, 1600, 640, 40)
    assert certificate.diagonal_satisfying_assignments == 40
    assert certificate.complement_census.assignments == 40**4
    assert certificate.complement_census.satisfying_assignments == 40
    assert certificate.complement_census.generated_order_histogram == ((3, 40),)
    assert certificate.complement_census.a6_assignments == 0


def test_full_delta_ten_t112_certificate_is_verified() -> None:
    """Every exact algebraic and computer-assisted invariant agrees."""

    assert exact_delta_ten_t112_certificate().verified
    assert exact_delta_ten_t112_certificate().generic_component_excluded
