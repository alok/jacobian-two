"""Exact tests for the split ``C3``/``C2^2`` rank-stratum closure."""

from sympy import Poly, expand, factor

from scripts.a6_delta_ten_split_contact_rank import (
    U,
    V,
    c3_rank_specs,
    exact_boundary_hostile_fixtures,
    exact_c3_rank_certificates,
    exact_residual_rank_certificates,
    exact_simple_double_contact_rank_certificates,
    exact_split_contact_rank_certificate,
    residual_rank_specs,
    simple_double_contact_specs,
    split_component_clean_factor,
)


def test_clean_root_factors_are_the_actual_overlap_cusp_diagonal_factors() -> None:
    assert expand(split_component_clean_factor(0, "V", U) - U * (2 * U - 1)) == 0
    assert (
        expand(split_component_clean_factor(0, "W", U) - U * (U**2 + 1) * (U**2 + 2))
        == 0
    )
    assert expand(split_component_clean_factor(2, "V", U) - U * (4 * U - 1)) == 0
    assert expand(split_component_clean_factor(2, "W", U) - U * (U + 1) * (U + 2)) == 0


def test_all_five_allowed_c3_systems_have_expected_clean_rank() -> None:
    specs = c3_rank_specs()
    certificates = exact_c3_rank_certificates()

    assert len(specs) == len(certificates) == 5
    assert [spec.allocation for spec in specs] == [
        "W",
        "overlap-V",
        "overlap-W",
        "V",
        "W",
    ]
    assert [certificate.coefficient_rank for certificate in certificates] == [
        3,
        2,
        2,
        3,
        3,
    ]
    assert all(
        certificate.coefficient_rank == certificate.augmented_rank
        for certificate in certificates
    )
    assert all(
        certificate.maximal_minor_gcd_identity == 0 for certificate in certificates
    )
    assert all(certificate.clean_rank_drop_removed for certificate in certificates)
    assert all(certificate.verified for certificate in certificates)


def test_visible_two_contact_rank_loss_is_confined_to_clean_boundaries() -> None:
    specs = simple_double_contact_specs()
    certificates = exact_simple_double_contact_rank_certificates()

    assert [spec.name for spec in specs] == [
        "c22_k0_vw",
        "c22_k0_overlap_w",
        "c22_k2_vv",
    ]
    assert [factor(spec.expected_rank_drop_generator) for spec in specs] == [
        256 * V**3 * (V**2 + 1),
        32 * V**2 * (V**2 + 1),
        (U - V) ** 4,
    ]
    assert [certificate.coefficient_rank for certificate in certificates] == [
        4,
        3,
        4,
    ]
    assert all(
        certificate.maximal_minor_gcd_identity == 0 for certificate in certificates
    )
    assert all(certificate.clean_rank_drop_removed for certificate in certificates)
    assert all(certificate.verified for certificate in certificates)


def test_three_residual_determinants_are_irreducible_and_have_finite_bases() -> None:
    specs = residual_rank_specs()
    certificates = exact_residual_rank_certificates()

    assert len(specs) == len(certificates) == 3
    assert [spec.allocation for spec in specs] == ["WW", "VW", "WW"]
    assert [
        Poly(certificate.residual_factor, U, V).total_degree()
        for certificate in certificates
    ] == [
        6,
        8,
        3,
    ]
    assert [
        len(Poly(certificate.residual_factor, U, V).terms())
        for certificate in certificates
    ] == [
        12,
        29,
        8,
    ]
    assert all(certificate.determinant_identity == 0 for certificate in certificates)
    assert all(certificate.residual_irreducible for certificate in certificates)
    assert [certificate.base_length for certificate in certificates] == [4, 6, 6]
    assert [certificate.unordered_base_orbit_count for certificate in certificates] == [
        2,
        None,
        3,
    ]


def test_ideal_sandwiches_and_rank_two_saturations_are_exact() -> None:
    certificates = exact_residual_rank_certificates()

    assert [
        certificate.compatibility_generator_count for certificate in certificates
    ] == [
        5,
        5,
        5,
    ]
    assert all(
        all(remainder == 0 for remainder in certificate.raw_in_expected_remainders)
        for certificate in certificates
    )
    assert all(
        certificate.powered_expected_in_raw_remainders == (0, 0)
        for certificate in certificates
    )
    assert all(
        any(remainder != 0 for remainder in certificate.lower_power_remainders)
        for certificate in certificates
    )
    assert all(
        certificate.expected_plus_boundary_basis == (1,) for certificate in certificates
    )
    assert all(
        certificate.univariate_derivative_gcd == 1 for certificate in certificates
    )

    assert all(
        certificate.rank_two_power_remainder == 0 for certificate in certificates
    )
    assert all(
        certificate.rank_two_lower_power_remainder != 0 for certificate in certificates
    )
    assert [certificate.coefficient_rank_on_base for certificate in certificates] == [
        3,
        3,
        3,
    ]
    assert [
        certificate.coefficient_fiber_dimension for certificate in certificates
    ] == [
        1,
        1,
        1,
    ]
    assert [certificate.incidence_dimension for certificate in certificates] == [
        1,
        1,
        1,
    ]
    assert all(certificate.verified for certificate in certificates)


def test_hostile_fixtures_make_the_clean_localizations_necessary() -> None:
    fixtures = exact_boundary_hostile_fixtures()

    assert [fixture.name for fixture in fixtures] == [
        "c3_k0_forced_cusp_pair",
        "c22_k0_component_overlap",
        "c22_k0_repeated_graph_contact",
        "c22_k2_repeated_vertical_contact",
        "c22_k2_overlap_diagonal_intersection",
    ]
    assert [fixture.coefficient_rank for fixture in fixtures] == [2, 3, 2, 2, 3]
    assert all(
        fixture.coefficient_rank == fixture.augmented_rank for fixture in fixtures
    )
    assert all(fixture.boundary_value == 0 for fixture in fixtures)
    assert all(
        all(value == 0 for value in fixture.equation_residuals) for fixture in fixtures
    )
    assert all(fixture.next_jet_value != 0 for fixture in fixtures)
    assert all(fixture.verified for fixture in fixtures)


def test_aggregate_closes_rank_not_topology_or_jc2() -> None:
    certificate = exact_split_contact_rank_certificate()

    assert certificate.maximum_residual_incidence_dimension == 1
    assert not certificate.topology_computed
    assert not certificate.proves_plane_jacobian_conjecture
    assert certificate.verified
