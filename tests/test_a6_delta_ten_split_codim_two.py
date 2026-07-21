"""Exact regression tests for the split delta-ten codimension-two ledger."""

from sympy import Rational, expand

from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    CUSP_IMAGE_FACTOR,
    DELTA,
    GAMMA,
    KAPPA,
)
from scripts.a6_delta_ten_split_codim_two import (
    FIRST_ROOT,
    FORBIDDEN_K0_ALLOCATIONS,
    PARAMETERS,
    SECOND_ROOT,
    SPLIT_ALLOCATIONS,
    exact_split_chart_certificates,
    exact_split_codimension_two_certificate,
    exact_split_rank_boundary_certificate,
    exact_split_rank_saturation_certificates,
    exact_split_witness_certificates,
    split_witness_specs,
)


def test_true_split_charts_have_exact_degrees_and_overlap_boundaries() -> None:
    charts = {chart.kappa: chart for chart in exact_split_chart_certificates()}

    assert set(charts) == {0, 2, -2}
    assert (charts[0].vertical_degree, charts[0].graph_degree) == (2, 8)
    assert (charts[2].vertical_degree, charts[2].graph_degree) == (4, 6)
    assert (charts[-2].vertical_degree, charts[-2].graph_degree) == (4, 6)
    assert charts[0].overlap_condition == 4 * ALPHA - 2 * GAMMA + 1
    assert charts[2].overlap_condition == ALPHA - BETA + GAMMA - DELTA + 1
    assert charts[-2].overlap_condition == ALPHA + BETA + GAMMA + DELTA + 1
    assert charts[0].discriminant_degrees == (2, 13)
    assert charts[2].discriminant_degrees == (6, 9)
    assert charts[-2].discriminant_degrees == (6, 9)
    assert all(
        chart.discriminants_squarefree == (True, True) for chart in charts.values()
    )
    assert all(chart.discriminant_gcd_degree == 0 for chart in charts.values())
    assert all(chart.verified for chart in charts.values())

    # At k=+/-2 the component intersection is exactly the extra cusp-image
    # boundary, so it cannot be silently counted as a clean extra node.
    for kappa in (-2, 2):
        assert charts[kappa].overlap_cusp_identity == 0
        assert (
            expand(
                CUSP_IMAGE_FACTOR.subs(KAPPA, kappa)
                - charts[kappa].overlap_condition ** 2
            )
            == 0
        )


def test_finite_allocation_ledger_fills_each_component_budget() -> None:
    assert len(SPLIT_ALLOCATIONS) == 22
    assert sum(row.split_class == "k=0" for row in SPLIT_ALLOCATIONS) == 11
    assert sum(row.split_class == "k=+/-2" for row in SPLIT_ALLOCATIONS) == 11
    assert all(row.expected_dimension == 2 for row in SPLIT_ALLOCATIONS)
    assert all(row.total_delta == 10 for row in SPLIT_ALLOCATIONS)
    assert all(row.verified for row in SPLIT_ALLOCATIONS)

    allowed_keys = {
        (row.profile, row.split_class, row.allocation) for row in SPLIT_ALLOCATIONS
    }
    for profile, allocation, _reason in FORBIDDEN_K0_ALLOCATIONS:
        assert (profile, "k=0", allocation) not in allowed_keys

    assert FORBIDDEN_K0_ALLOCATIONS == (
        ("C3+7N", "V", "the monic quadratic V cannot have a triple root"),
        (
            "C2^2+6N",
            "VV",
            "two distinct double roots require V-degree at least four",
        ),
        (
            "C2+T111+5N",
            "contact-V",
            "the triple V edge plus a separate double V edge has length three",
        ),
    )


def test_generic_rank_obstructions_and_principal_saturations_are_exact() -> None:
    rank = exact_split_rank_boundary_certificate()

    # This inconsistent augmented system is the algebraic obstruction to a
    # non-overlap C3 on the k=0 quadratic vertical component.
    assert rank.k0_vertical_c3_ranks == (2, 3)
    assert rank.k0_vv_double_contact_determinant == 0
    assert (
        expand(
            rank.k0_vw_double_contact_determinant
            - 256 * SECOND_ROOT**3 * (SECOND_ROOT**2 + 1)
        )
        == 0
    )
    assert (
        expand(rank.k2_vv_double_contact_determinant - (FIRST_ROOT - SECOND_ROOT) ** 4)
        == 0
    )
    assert rank.plus_minus_vertical_symmetry == 0
    assert rank.plus_minus_graph_symmetry == 0
    assert rank.verified

    saturations = {
        certificate.name: certificate
        for certificate in exact_split_rank_saturation_certificates()
    }
    assert set(saturations) == {
        "k0_graph_c3",
        "k0_vw_double_contact",
        "k2_vv_double_contact",
    }
    assert saturations["k0_graph_c3"].saturation_exponent == 3
    assert saturations["k0_graph_c3"].quotient == Rational(1, 64)
    assert saturations["k0_vw_double_contact"].saturation_exponent == 3
    assert (
        expand(
            saturations["k0_vw_double_contact"].quotient
            - (SECOND_ROOT**2 + 1) ** 2 / 256
        )
        == 0
    )
    assert saturations["k2_vv_double_contact"].saturation_exponent == 4
    assert saturations["k2_vv_double_contact"].quotient == 1
    assert all(
        certificate.saturated_ideal_is_unit for certificate in saturations.values()
    )
    assert all(certificate.verified for certificate in saturations.values())


def test_every_allowed_allocation_has_a_clean_exact_witness() -> None:
    specs = split_witness_specs()
    witnesses = exact_split_witness_certificates()

    assert len(specs) == len(witnesses) == len(SPLIT_ALLOCATIONS) == 22
    for spec, witness in zip(specs, witnesses, strict=True):
        assert witness.verified(spec), spec.name
        assert witness.coefficient_rank == witness.augmented_rank == spec.expected_rank
        assert witness.incidence_dimension == 2
        assert witness.vertical_residual_degree + witness.graph_residual_degree == (
            spec.node_count
        )
        assert witness.residual_squarefree == (True, True)
        assert witness.node_targets_squarefree == (True, True)
        assert witness.cross_node_target_gcd_degree == 0
        assert witness.special_targets_distinct
        assert witness.special_node_separations_nonzero
        assert witness.cusp_image_nonzero
        assert witness.extra_critical_nonzero
        assert witness.overlap_clean

    aggregate = exact_split_codimension_two_certificate()
    assert aggregate.witness_keys_cover_allocations
    assert aggregate.critical_triple_boundary_open
    assert aggregate.topology_not_computed
    assert aggregate.verified


def test_overlap_contacts_are_legitimate_but_boundary_candidates_are_rejected() -> None:
    rank = exact_split_rank_boundary_certificate()
    assert rank.k0_overlap_value_identity == 0
    assert rank.k0_overlap_vertical_derivative == 1 - GAMMA
    assert rank.k0_overlap_graph_derivative == 4 * (3 * BETA - 2 * DELTA)

    specs = {spec.name: spec for spec in split_witness_specs()}
    overlap_vertical = specs["c3_k0_overlap_v"]
    overlap_graph = specs["c3_k0_overlap_w"]
    assert tuple(
        (root.component, root.root, root.multiplicity)
        for root in overlap_vertical.roots
    ) == (("V", Rational(1, 2), 2), ("W", 0, 1))
    assert tuple(
        (root.component, root.root, root.multiplicity) for root in overlap_graph.roots
    ) == (("V", Rational(1, 2), 1), ("W", 0, 2))

    # Hostile fixture: this point solves the k=2 vertical C3 equations, but it
    # lies at V(0)=0, hence also on the cusp-image boundary.  Rank and contact
    # equations alone must not promote it to a clean witness.
    boundary_coefficients = {
        ALPHA: -1,
        BETA: -1,
        GAMMA: 1,
        DELTA: 2,
    }
    k2_vertical = specs["c3_k2_v"]
    assert all(
        expand(equation.subs(boundary_coefficients)) == 0
        for equation in k2_vertical.equations
    )
    assert CUSP_IMAGE_FACTOR.subs(KAPPA, 2).subs(boundary_coefficients) == 0

    # The checked witness solves the same incidence equations away from every
    # cusp, critical, overlap, residual-collision, and target-collision wall.
    clean_coefficients = dict(
        zip(PARAMETERS, k2_vertical.coefficient_values, strict=True)
    )
    assert all(
        expand(equation.subs(clean_coefficients)) == 0
        for equation in k2_vertical.equations
    )
    assert CUSP_IMAGE_FACTOR.subs(KAPPA, 2).subs(clean_coefficients) != 0
