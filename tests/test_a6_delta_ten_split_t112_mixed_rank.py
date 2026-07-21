"""Exact regressions for the focused true-split rank audit."""

from sympy import Rational, expand

from scripts.a6_delta_ten_split_t112_mixed_rank import (
    CONTACT_ROOT,
    MIXED_COMPATIBILITY_BASES,
    SPLIT_INCIDENCE_SPECS,
    THIRD_ROOT,
    TRIPLE_ROOT,
    exact_component_derivation_certificates,
    exact_hostile_split_rank_fixtures,
    exact_plus_minus_full_system_transport_certificates,
    exact_split_incidence_certificates,
    exact_split_t112_mixed_rank_certificate,
    exact_target_separation_certificates,
    split_incidence_equations,
    split_incidence_matrices,
    split_triple_geometry,
    valid_base_localizer,
)


def test_true_component_equations_derive_both_split_triple_fibers() -> None:
    geometries = tuple(split_triple_geometry(kappa) for kappa in (0, 2))
    derivations = exact_component_derivation_certificates()

    assert all(geometry.verified for geometry in geometries)
    assert all(derivation.verified for derivation in derivations)
    assert geometries[0].base_constraint == TRIPLE_ROOT**2 + THIRD_ROOT**2 + 1
    assert geometries[0].vertical_root == -(TRIPLE_ROOT**2)
    assert geometries[0].graph_roots == (
        TRIPLE_ROOT + THIRD_ROOT,
        THIRD_ROOT - TRIPLE_ROOT,
    )
    assert geometries[1].base_constraint == expand(
        TRIPLE_ROOT * (TRIPLE_ROOT + 1) + THIRD_ROOT * (THIRD_ROOT + 1)
    )
    assert geometries[1].vertical_root == expand(-TRIPLE_ROOT * (TRIPLE_ROOT + 1))
    assert geometries[1].graph_roots == (
        TRIPLE_ROOT + THIRD_ROOT,
        THIRD_ROOT - TRIPLE_ROOT - 1,
    )


def test_all_seven_component_systems_have_exact_expected_generic_rank() -> None:
    certificates = exact_split_incidence_certificates()

    assert len(SPLIT_INCIDENCE_SPECS) == len(certificates) == 7
    assert sum(spec.profile == "T112+6N" for spec in SPLIT_INCIDENCE_SPECS) == 4
    assert sum(spec.profile == "C2+T111+5N" for spec in SPLIT_INCIDENCE_SPECS) == 3
    for certificate in certificates:
        spec = certificate.spec
        equations = split_incidence_equations(spec)
        coefficient, augmented = split_incidence_matrices(spec)
        assert len(equations) == (3 if spec.profile == "T112+6N" else 4)
        assert coefficient.rows == augmented.rows == len(equations)
        assert coefficient.cols == 4
        assert augmented.cols == 5
        assert certificate.coefficient_rank == spec.expected_rank
        assert certificate.augmented_rank == spec.expected_rank
        assert certificate.witness_base_residual == 0
        assert certificate.witness_localizer_value != 0
        assert all(value == 0 for value in certificate.witness_equation_residuals)
        assert valid_base_localizer(spec) != 0
        assert spec.generic_incidence_dimension == 2
        assert spec.hidden_threefold_excluded
        assert certificate.verified


def test_sage_rank_and_compatibility_results_bound_every_residual_fiber() -> None:
    specs = {spec.name: spec for spec in SPLIT_INCIDENCE_SPECS}

    for name in ("t112_k0_v", "t112_k0_w"):
        assert specs[name].rank_drop_saturation_exponent == 1
        assert specs[name].compatibility_saturation_exponent == 1
        assert specs[name].rank_drop_valid_dimension is None
        assert specs[name].residual_incidence_dimension_bound == -1
    for name in ("t112_k2_v", "t112_k2_w"):
        assert specs[name].rank_drop_saturation_exponent == 2
        assert specs[name].compatibility_saturation_exponent == 2
        assert specs[name].rank_drop_valid_dimension is None
        assert specs[name].residual_incidence_dimension_bound == -1

    k0_graph = specs["mixed_k0_w"]
    assert k0_graph.rank_drop_valid_dimension == 1
    assert k0_graph.rank_drop_prime_count == 1
    assert k0_graph.compatibility_valid_dimension == 0
    assert k0_graph.compatibility_length == 4
    assert k0_graph.compatibility_reduced
    assert k0_graph.residual_incidence_dimension_bound == 1
    assert MIXED_COMPATIBILITY_BASES["mixed_k0_w"] == (
        TRIPLE_ROOT**2 + Rational(11, 12),
        CONTACT_ROOT**2 + Rational(1, 3),
        THIRD_ROOT - CONTACT_ROOT / 2,
    )

    k2_vertical = specs["mixed_k2_v"]
    assert k2_vertical.rank_drop_valid_dimension == 1
    assert k2_vertical.rank_drop_prime_count == 1
    assert k2_vertical.compatibility_valid_dimension == 0
    assert k2_vertical.compatibility_length == 4
    assert k2_vertical.compatibility_reduced
    assert k2_vertical.residual_incidence_dimension_bound == 1
    assert MIXED_COMPATIBILITY_BASES["mixed_k2_v"] == (
        TRIPLE_ROOT**2 + TRIPLE_ROOT - CONTACT_ROOT + Rational(1, 4),
        CONTACT_ROOT**2 - CONTACT_ROOT + Rational(1, 2),
        THIRD_ROOT - CONTACT_ROOT + Rational(3, 2),
    )

    k2_graph = specs["mixed_k2_w"]
    assert k2_graph.rank_drop_valid_dimension == 1
    assert k2_graph.rank_drop_prime_count == 2
    assert k2_graph.compatibility_valid_dimension is None
    assert k2_graph.compatibility_length == 0
    assert not k2_graph.compatibility_reduced
    assert k2_graph.residual_incidence_dimension_bound == -1
    assert MIXED_COMPATIBILITY_BASES["mixed_k2_w"] == (1,)

    for spec in (k0_graph, k2_vertical, k2_graph):
        assert spec.rank_drop_saturation_exponent == 2
        assert spec.compatibility_saturation_exponent == 2
        assert spec.coefficient_rank_two_saturation_exponent == 1
        assert spec.augmented_rank_two_saturation_exponent == 1


def test_contact_localizers_are_exact_same_target_and_unramified_boundaries() -> None:
    certificates = exact_target_separation_certificates()

    assert {certificate.name for certificate in certificates} == {
        "mixed_k0_w",
        "mixed_k2_v",
        "mixed_k2_w",
    }
    assert all(certificate.identity == 0 for certificate in certificates)
    assert all(
        certificate.contact_derivative_factor != 0 for certificate in certificates
    )
    assert all(certificate.contact_diagonal_factor != 0 for certificate in certificates)
    assert all(certificate.verified for certificate in certificates)


def test_kminus_two_is_used_only_after_full_incidence_transport() -> None:
    transports = exact_plus_minus_full_system_transport_certificates()

    assert len(transports) == 4
    assert {transport.name for transport in transports} == {
        "t112_k2_v",
        "t112_k2_w",
        "mixed_k2_v",
        "mixed_k2_w",
    }
    for transport in transports:
        assert transport.base_identity == 0
        assert transport.p_identity == 0
        assert transport.q_identity == 0
        assert all(identity == 0 for identity in transport.equation_identities)
        assert transport.verified


def test_hostile_rank_two_points_are_compatible_but_removed() -> None:
    fixtures = exact_hostile_split_rank_fixtures()

    assert len(fixtures) == len(SPLIT_INCIDENCE_SPECS) == 7
    assert {fixture.spec_name for fixture in fixtures} == {
        spec.name for spec in SPLIT_INCIDENCE_SPECS
    }
    for fixture in fixtures:
        assert fixture.base_residual == 0
        assert fixture.localizer_value == 0
        assert fixture.coefficient_rank == fixture.augmented_rank == 2
        assert fixture.verified


def test_aggregate_closes_only_the_named_split_rank_strata() -> None:
    certificate = exact_split_t112_mixed_rank_certificate()

    assert len(certificate.audited_allocations) == 7
    assert all(
        incidence.spec.hidden_threefold_excluded for incidence in certificate.incidences
    )
    assert certificate.split_rank_strata_closed
    assert certificate.verified
    assert certificate.unaudited_split_allocations_remain
    assert certificate.topology_not_computed
