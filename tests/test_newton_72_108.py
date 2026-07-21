"""Regression tests for the residual ``(72, 108)`` sparse-support audit."""

from dataclasses import replace
from itertools import combinations
from math import comb

import sympy as sp

from scripts.newton_72_108 import (
    CASE_1,
    CASE_2,
    CASE_2_ALGEBRAIC_CERTIFICATES,
    CONFIGURATIONS,
    BracketTerm,
    Coefficient,
    ZeroProductCertificate,
    audit_sparse_interior_bound,
    boundary_supports,
    bracket_equations,
    find_relaxed_full_polygon_certificate,
    find_zero_product_certificate,
    interior_coefficients,
    polygon_boundary_lattice_points,
    polygon_interior_lattice_points,
    polygon_lattice_points,
    verify_algebraic_support_certificate,
    verify_zero_product_certificate,
)


def test_exact_polygon_lattice_counts() -> None:
    assert tuple(
        (
            len(polygon_lattice_points(configuration.p_vertices)),
            len(polygon_boundary_lattice_points(configuration.p_vertices)),
            len(polygon_interior_lattice_points(configuration.p_vertices)),
            len(polygon_lattice_points(configuration.q_vertices)),
            len(polygon_boundary_lattice_points(configuration.q_vertices)),
            len(polygon_interior_lattice_points(configuration.q_vertices)),
        )
        for configuration in CONFIGURATIONS
    ) == (
        (61, 26, 35, 125, 38, 87),
        (25, 18, 7, 47, 26, 21),
    )


def test_boundary_only_case_one_has_an_immediate_nonzero_coefficient() -> None:
    equations = bracket_equations(*boundary_supports(CASE_1))

    assert equations[(9, 16)] == (
        BracketTerm(
            p=Coefficient("P", (8, 16)),
            q=Coefficient("Q", (2, 1)),
            scalar=-24,
        ),
    )


def test_boundary_only_case_two_has_a_three_equation_hand_certificate() -> None:
    equations = bracket_equations(*boundary_supports(CASE_2))

    assert equations[(9, 15)] == (
        BracketTerm(
            p=Coefficient("P", (8, 15)),
            q=Coefficient("Q", (2, 1)),
            scalar=-22,
        ),
    )
    assert equations[(19, 35)] == (
        BracketTerm(
            p=Coefficient("P", (8, 14)),
            q=Coefficient("Q", (12, 22)),
            scalar=8,
        ),
        BracketTerm(
            p=Coefficient("P", (8, 15)),
            q=Coefficient("Q", (12, 21)),
            scalar=-12,
        ),
    )
    assert equations[(19, 37)] == (
        BracketTerm(
            p=Coefficient("P", (8, 14)),
            q=Coefficient("Q", (12, 24)),
            scalar=24,
        ),
        BracketTerm(
            p=Coefficient("P", (8, 15)),
            q=Coefficient("Q", (12, 23)),
            scalar=4,
        ),
        BracketTerm(
            p=Coefficient("P", (8, 16)),
            q=Coefficient("Q", (12, 22)),
            scalar=-16,
        ),
    )


def test_sympy_independently_recovers_the_boundary_coefficient_identities() -> None:
    x, y = sp.symbols("x y")

    for configuration, target_exponents in (
        (CASE_1, ((9, 16),)),
        (CASE_2, ((9, 15), (19, 35), (19, 37))),
    ):
        p_support, q_support = boundary_supports(configuration)
        p_coefficients = {
            point: sp.Symbol(f"p_{point[0]}_{point[1]}") for point in p_support
        }
        q_coefficients = {
            point: sp.Symbol(f"q_{point[0]}_{point[1]}") for point in q_support
        }
        p_polynomial = sum(
            p_coefficients[point] * x ** point[0] * y ** point[1] for point in p_support
        )
        q_polynomial = sum(
            q_coefficients[point] * x ** point[0] * y ** point[1] for point in q_support
        )
        bracket = sp.Poly(
            sp.diff(p_polynomial, x) * sp.diff(q_polynomial, y)
            - sp.diff(p_polynomial, y) * sp.diff(q_polynomial, x),
            x,
            y,
        )
        equations = bracket_equations(p_support, q_support)

        for exponent in target_exponents:
            expected = sum(
                term.scalar
                * p_coefficients[term.p.exponent]
                * q_coefficients[term.q.exponent]
                for term in equations[exponent]
            )
            assert (
                sp.expand(
                    bracket.coeff_monomial(x ** exponent[0] * y ** exponent[1])
                    - expected
                )
                == 0
            )


def test_case_one_needs_three_and_case_two_needs_four_interiors() -> None:
    summaries = (
        audit_sparse_interior_bound(CASE_1, maximum_interior_terms=2),
        audit_sparse_interior_bound(CASE_2, maximum_interior_terms=3),
    )

    assert tuple(
        (summary.support_patterns, summary.verified_certificates)
        for summary in summaries
    ) == ((7504, 7504), (3683, 3683))
    assert tuple(
        (summary.zero_product_certificates, summary.algebraic_certificates)
        for summary in summaries
    ) == ((7504, 0), (3678, 5))
    assert tuple(summary.maximum_forced_zero_steps for summary in summaries) == (
        21,
        22,
    )
    assert tuple(summary.certificate_sha256 for summary in summaries) == (
        "7c7a3824ac64f0f94989125079be5b4ae7a54817c7eaa4c3b3f2bf79b5bb1519",
        "ee64dd2ba8035ee29c4efd86ef7f80fc049f8ecd5d32aec638f719cfd31f83fb",
    )


def test_all_five_zero_product_exceptions_have_exact_algebraic_replays() -> None:
    possible = interior_coefficients(CASE_2)
    survivors = tuple(
        selected
        for selected in combinations(possible, 3)
        if find_zero_product_certificate(CASE_2, selected) is None
    )

    assert comb(len(possible), 3) == 3276
    assert survivors == tuple(
        certificate.selected_interior for certificate in CASE_2_ALGEBRAIC_CERTIFICATES
    )
    assert len(survivors) == 5
    assert all(
        verify_algebraic_support_certificate(CASE_2, certificate)
        for certificate in CASE_2_ALGEBRAIC_CERTIFICATES
    )


def test_full_polygons_are_not_claimed_and_four_terms_can_stop_the_method() -> None:
    assert all(
        find_relaxed_full_polygon_certificate(configuration) is None
        for configuration in CONFIGURATIONS
    )

    four_interior_terms = (
        Coefficient("P", (2, 3)),
        Coefficient("P", (3, 5)),
        Coefficient("Q", (3, 4)),
        Coefficient("Q", (4, 7)),
    )
    assert find_zero_product_certificate(CASE_2, four_interior_terms) is None


def test_fresh_replay_rejects_a_tampered_certificate() -> None:
    certificate = find_zero_product_certificate(CASE_2)
    assert certificate is not None
    assert verify_zero_product_certificate(CASE_2, (), certificate)

    tampered_contradiction = replace(
        certificate.contradiction,
        term=replace(certificate.contradiction.term, scalar=1),
    )
    tampered = ZeroProductCertificate(
        forced_zeros=certificate.forced_zeros,
        contradiction=tampered_contradiction,
    )
    assert not verify_zero_product_certificate(CASE_2, (), tampered)
