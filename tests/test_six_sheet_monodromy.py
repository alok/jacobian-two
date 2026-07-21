"""Regression tests for the exact six-sheet transitive-group certificate."""

import pytest

from scripts.six_sheet_monodromy import (
    A6_COLLISION_MERIDIAN,
    A6_CONTRACTED_SOURCE_ENDPOINTS,
    A6_MINIMAL_DICRITICAL_LEAF_LABEL_ENDPOINT,
    A6_EXCEPTIONAL_PARITY_ENDPOINT,
    A6_TORUS_2_5_LOCAL_GENERATORS,
    OREVKOV_FORBIDDEN_RAMIFICATION_INDEX,
    OREVKOV_ONE_DICRITICAL_TYPES,
    OREVKOV_SIX_SHEET_BUDGET,
    OrevkovBudgetProfile,
    DicriticalLeafLabelEndpoint,
    OneDicriticalContractedSourceEndpoint,
    S6_TWO_CURVE_COLLISION_GENERATORS,
    S6_CONTRACTED_SOURCE_ENDPOINTS,
    S6_MINIMAL_DICRITICAL_LEAF_LABEL_ENDPOINT,
    S6_SATURATED_CONTRACTED_SOURCE_ENDPOINTS,
    S6_TWO_CURVE_FIBER_PROFILES,
    TRANSITIVE_GROUPS,
    branch_inertia_realizations,
    candidate_ids,
    classify_six_sheet_groups,
    compose,
    conjugacy_class,
    conjugacy_class_analyses,
    full_generating_inertia_pairs,
    generated_group,
    normally_generated_overgroups,
    one_dicritical_passports,
    one_dicritical_s6_fiber_profiles,
    orevkov_budget_profiles,
    permutation_orbit_sizes,
    ramified_branch_profiles,
)


def test_catalogue_has_the_sixteen_gap_transitive_groups() -> None:
    rows = classify_six_sheet_groups()

    assert tuple(row.identifier for row in rows) == tuple(
        f"6T{index}" for index in range(1, 17)
    )
    assert tuple(row.order for row in rows) == (
        6,
        6,
        12,
        12,
        18,
        24,
        24,
        24,
        36,
        36,
        48,
        60,
        72,
        120,
        360,
        720,
    )


def test_only_regular_groups_are_galois() -> None:
    rows = classify_six_sheet_groups()
    galois = tuple(row.identifier for row in rows if row.is_galois)

    assert galois == ("6T1", "6T2")
    assert tuple(row.automorphism_group for row in rows[:2]) == ("C6", "S3")
    assert all(row.stabilizer_order == 1 for row in rows[:2])


def test_exact_deck_group_normalizer_quotients() -> None:
    rows = classify_six_sheet_groups()

    assert tuple(row.automorphism_order for row in rows) == (
        6,
        6,
        2,
        2,
        3,
        2,
        2,
        2,
        1,
        1,
        2,
        1,
        1,
        1,
        1,
        1,
    )
    assert candidate_ids(trivial_automorphisms=True) == (
        "6T9",
        "6T10",
        "6T12",
        "6T13",
        "6T14",
        "6T15",
        "6T16",
    )


def test_primitive_candidates_are_exactly_the_four_self_normalizing_cases() -> None:
    rows = classify_six_sheet_groups()
    primitive = tuple(row for row in rows if row.is_primitive)

    assert tuple(row.identifier for row in primitive) == (
        "6T12",
        "6T14",
        "6T15",
        "6T16",
    )
    assert all(row.automorphism_order == 1 for row in primitive)
    assert candidate_ids(primitive=True, trivial_automorphisms=True) == (
        "6T12",
        "6T14",
        "6T15",
        "6T16",
    )


def test_every_non_galois_group_with_deck_symmetry_is_imprimitive() -> None:
    rows = classify_six_sheet_groups()
    deck_symmetric = tuple(
        row
        for row in rows
        if not row.is_galois and row.automorphism_order > 1
    )

    assert tuple(row.identifier for row in deck_symmetric) == (
        "6T3",
        "6T4",
        "6T5",
        "6T6",
        "6T7",
        "6T8",
        "6T11",
    )
    assert all(not row.is_primitive for row in deck_symmetric)


def test_affine_and_refined_unrestricted_filters_are_exact() -> None:
    rows = classify_six_sheet_groups()

    assert tuple(row.fixed_point_normal_closure_order for row in rows) == (
        1,
        1,
        6,
        4,
        9,
        8,
        24,
        24,
        18,
        18,
        48,
        60,
        36,
        120,
        360,
        720,
    )
    assert candidate_ids(affine_inertia_compatible=True) == (
        "6T7",
        "6T8",
        "6T11",
        "6T12",
        "6T14",
        "6T15",
        "6T16",
    )
    assert candidate_ids(general_keller_compatible=True) == (
        "6T15",
        "6T16",
    )


def test_refined_primitive_branch_profiles_leave_only_a6_and_s6() -> None:
    profiles = {
        fixture.identifier: ramified_branch_profiles(fixture)
        for fixture in TRANSITIVE_GROUPS[11:]
        if fixture.identifier in {"6T12", "6T14", "6T15", "6T16"}
    }

    assert {
        identifier: tuple(
            (
                tuple(
                    (
                        analysis.cycle_type,
                        analysis.class_size,
                        analysis.normal_closure_order,
                        analysis.moved_sheet_count,
                    )
                    for analysis in profile.branch_classes
                ),
                profile.minimum_defect_cost,
                profile.residual_defect,
                profile.survives_irreducible_branch_obstruction,
            )
            for profile in group_profiles
        )
        for identifier, group_profiles in profiles.items()
    } == {
        "6T12": (
            ((((2, 2, 1, 1), 15, 60, 4),), 4, 1, False),
        ),
        "6T14": (
            ((((4, 1, 1), 30, 120, 4),), 4, 1, False),
        ),
        "6T15": (
            ((((3, 1, 1, 1), 40, 360, 3),), 3, 2, True),
            ((((2, 2, 1, 1), 45, 360, 4),), 4, 1, False),
        ),
        "6T16": (
            ((((2, 1, 1, 1, 1), 15, 720, 2),), 2, 3, True),
            ((((4, 1, 1), 90, 720, 4),), 4, 1, False),
            (
                (
                    ((2, 1, 1, 1, 1), 15, 720, 2),
                    ((2, 1, 1, 1, 1), 15, 720, 2),
                ),
                4,
                1,
                True,
            ),
            ((((3, 2, 1), 120, 720, 5),), 5, 0, False),
            (
                (
                    ((2, 1, 1, 1, 1), 15, 720, 2),
                    ((3, 1, 1, 1), 40, 360, 3),
                ),
                5,
                0,
                True,
            ),
        ),
    }

    assert all(
        OREVKOV_FORBIDDEN_RAMIFICATION_INDEX
        not in analysis.cycle_type
        for group_profiles in profiles.values()
        for profile in group_profiles
        for analysis in profile.branch_classes
    )


def test_a6_torus_2_5_local_monodromy_is_a_hostile_survivor() -> None:
    first, second = A6_TORUS_2_5_LOCAL_GENERATORS

    def word(*letters: tuple[int, ...]) -> tuple[int, ...]:
        result = tuple(range(6))
        for letter in letters:
            result = compose(result, letter)
        return result

    # Meridional Artin presentation of the T(2,5) knot group.
    assert word(first, second, first, second, first) == word(
        second,
        first,
        second,
        first,
        second,
    )
    local_group = generated_group(A6_TORUS_2_5_LOCAL_GENERATORS)
    global_group = generated_group(
        (*A6_TORUS_2_5_LOCAL_GENERATORS, A6_COLLISION_MERIDIAN)
    )

    assert len(local_group) == 60
    assert permutation_orbit_sizes(local_group) == (5, 1)
    assert len(global_group) == 360
    assert global_group == generated_group(TRANSITIVE_GROUPS[14].generators)


def test_a6_exceptional_source_parity_endpoint_is_inconsistent() -> None:
    endpoint = A6_EXCEPTIONAL_PARITY_ENDPOINT

    assert endpoint.lifted_degree == 10
    assert endpoint.forced_tangent_contact_order == 4
    assert endpoint.contradicts_odd_deck_invariance


def test_every_one_dicritical_contracted_source_endpoint_is_eliminated() -> None:
    assert tuple(
        (
            endpoint.local_degree,
            endpoint.jacobian_order,
            endpoint.forced_tangent_contact_order,
            endpoint.contradicts_odd_deck_invariance,
        )
        for endpoint in A6_CONTRACTED_SOURCE_ENDPOINTS
    ) == (
        (3, 2, 2, True),
        (5, 2, 4, True),
    )
    assert tuple(
        (
            endpoint.local_degree,
            endpoint.jacobian_order,
            endpoint.invariant_jacobian_lower_bound,
            endpoint.contradicts_invariant_order,
        )
        for endpoint in S6_CONTRACTED_SOURCE_ENDPOINTS
    ) == (
        (2, 1, 2, True),
        (3, 1, 2, True),
        (4, 1, 2, True),
        (5, 1, 2, True),
    )
    assert all(
        endpoint.eliminates_contracted_source
        for endpoint in (
            *A6_CONTRACTED_SOURCE_ENDPOINTS,
            *S6_CONTRACTED_SOURCE_ENDPOINTS,
        )
    )


def test_contracted_source_endpoint_preserves_its_exact_boundary() -> None:
    even_index_three = OneDicriticalContractedSourceEndpoint(
        ramification_index=3,
        local_degree=4,
    )

    assert even_index_three.forced_tangent_contact_order == 3
    assert not even_index_three.contradicts_invariant_order
    assert not even_index_three.contradicts_odd_deck_invariance
    assert not even_index_three.eliminates_contracted_source
    with pytest.raises(ValueError, match="ramification index 2 or 3"):
        OneDicriticalContractedSourceEndpoint(
            ramification_index=4,
            local_degree=5,
        )
    with pytest.raises(ValueError, match="local degree must be positive"):
        OneDicriticalContractedSourceEndpoint(
            ramification_index=3,
            local_degree=0,
        )


def test_saturated_s6_contracts_neither_boundary_chain() -> None:
    assert tuple(
        (
            endpoint.ramification_index,
            endpoint.local_degree,
            endpoint.eliminates_contracted_source,
        )
        for endpoint in S6_SATURATED_CONTRACTED_SOURCE_ENDPOINTS
    ) == (
        (2, 2, True),
        (3, 3, True),
    )


def test_joint_minimal_leaf_and_corner_blowup_orbits_are_exact() -> None:
    s6_endpoint = S6_MINIMAL_DICRITICAL_LEAF_LABEL_ENDPOINT
    a6_endpoint = A6_MINIMAL_DICRITICAL_LEAF_LABEL_ENDPOINT
    assert (s6_endpoint.dicritical_label, s6_endpoint.neighbor_label) == (2, 1)
    assert (a6_endpoint.dicritical_label, a6_endpoint.neighbor_label) == (3, 2)
    assert s6_endpoint.is_jointly_minimal
    assert a6_endpoint.is_jointly_minimal

    s6_orbit = [s6_endpoint]
    a6_orbit = [a6_endpoint]
    for _ in range(5):
        s6_orbit.append(s6_orbit[-1].corner_blowup())
        a6_orbit.append(a6_orbit[-1].corner_blowup())
    assert tuple(endpoint.neighbor_label for endpoint in s6_orbit) == (
        1,
        3,
        5,
        7,
        9,
        11,
    )
    assert tuple(endpoint.neighbor_label for endpoint in a6_orbit) == (
        2,
        5,
        8,
        11,
        14,
        17,
    )
    assert all(
        endpoint.dicritical_label == endpoint.ramification_index
        and endpoint.neighbor_label > 0
        and endpoint.has_forced_congruence
        for endpoint in (*s6_orbit, *a6_orbit)
    )
    assert not s6_orbit[1].is_jointly_minimal
    assert not a6_orbit[1].is_jointly_minimal
    with pytest.raises(ValueError, match="index must be at least two"):
        DicriticalLeafLabelEndpoint(
            ramification_index=1,
            negative_self_intersection=1,
        )
    with pytest.raises(ValueError, match="negative self-intersection"):
        DicriticalLeafLabelEndpoint(
            ramification_index=2,
            negative_self_intersection=0,
        )


def test_one_dicritical_s6_fiber_census_is_exhaustive() -> None:
    profiles = one_dicritical_s6_fiber_profiles()

    assert tuple(
        (
            profile.jumps,
            profile.boundary_block_sizes,
            profile.affine_fiber_size,
            profile.is_omitted,
        )
        for profile in profiles
    ) == (
        ((0,), (2,), 4, False),
        ((1,), (3,), 3, False),
        ((2,), (4,), 2, False),
        ((3,), (5,), 1, False),
        ((0, 0), (2, 2), 2, False),
        ((0, 1), (2, 3), 1, False),
        ((0, 2), (2, 4), 0, True),
        ((1, 1), (3, 3), 0, True),
        ((0, 0, 0), (2, 2, 2), 0, True),
    )


def test_saturated_two_curve_s6_collision_rows_are_exact() -> None:
    assert tuple(
        (
            profile.event,
            profile.boundary_block_sizes,
            profile.affine_fiber_size,
            profile.is_omitted,
        )
        for profile in S6_TWO_CURVE_FIBER_PROFILES
    ) == (
        ("three_cycle_self_collision", (3, 3), 0, True),
        ("transposition_double_self_collision", (2, 2), 2, False),
        ("transposition_triple_self_collision", (2, 2, 2), 0, True),
        ("cross_intersection", (2, 3), 1, False),
    )


def test_saturated_two_curve_collision_fixture_still_generates_s6() -> None:
    first_three_cycle, second_three_cycle, first_transposition, second_transposition = (
        S6_TWO_CURVE_COLLISION_GENERATORS
    )
    three_cycle_collision = generated_group(
        (first_three_cycle, second_three_cycle)
    )
    transposition_collision = generated_group(
        (first_transposition, second_transposition)
    )
    global_group = generated_group(S6_TWO_CURVE_COLLISION_GENERATORS)

    assert compose(first_three_cycle, second_three_cycle) == compose(
        second_three_cycle,
        first_three_cycle,
    )
    assert compose(first_transposition, second_transposition) == compose(
        second_transposition,
        first_transposition,
    )
    assert len(three_cycle_collision) == 9
    assert permutation_orbit_sizes(three_cycle_collision) == (3, 3)
    assert len(transposition_collision) == 4
    assert permutation_orbit_sizes(transposition_collision) == (2, 2, 1, 1)
    assert len(global_group) == 720
    assert global_group == generated_group(TRANSITIVE_GROUPS[15].generators)


def test_one_dicritical_filter_keeps_distinct_same_cycle_type_classes() -> None:
    rows = classify_six_sheet_groups()
    group_6t7 = generated_group(
        next(
            fixture.generators
            for fixture in TRANSITIVE_GROUPS
            if fixture.identifier == "6T7"
        )
    )
    double_transpositions = tuple(
        analysis
        for analysis in conjugacy_class_analyses(group_6t7)
        if analysis.cycle_type == (2, 2, 1, 1)
    )

    assert len(double_transpositions) == 2
    assert tuple(
        sorted(analysis.normal_closure_order for analysis in double_transpositions)
    ) == (4, 24)
    assert rows[6].one_dicritical_inertia_compatible


def test_orevkov_one_dicritical_cycle_types_and_survivors() -> None:
    rows = classify_six_sheet_groups()

    assert OREVKOV_ONE_DICRITICAL_TYPES == frozenset(
        {
            (2, 1, 1, 1, 1),
            (2, 2, 1, 1),
            (3, 1, 1, 1),
            (4, 1, 1),
        }
    )
    assert candidate_ids(one_dicritical_inertia_compatible=True) == (
        "6T7",
        "6T8",
        "6T12",
        "6T14",
        "6T15",
        "6T16",
    )
    assert candidate_ids(one_dicritical_keller_compatible=True) == (
        "6T15",
        "6T16",
    )
    assert candidate_ids(one_dicritical_deck_compatible=True) == (
        "6T12",
        "6T14",
        "6T15",
        "6T16",
    )
    assert not rows[10].one_dicritical_inertia_compatible
    assert rows[10].affine_inertia_compatible
    assert {
        row.identifier: tuple(
            (
                witness.cycle_type,
                witness.class_size,
                witness.normal_closure_order,
            )
            for witness in row.one_dicritical_witnesses
        )
        for row in rows
        if row.one_dicritical_witnesses
    } == {
        "6T7": (((2, 2, 1, 1), 6, 24),),
        "6T8": (((4, 1, 1), 6, 24),),
        "6T12": (((2, 2, 1, 1), 15, 60),),
        "6T14": (((4, 1, 1), 30, 120),),
        "6T15": (
            ((2, 2, 1, 1), 45, 360),
            ((3, 1, 1, 1), 40, 360),
        ),
        "6T16": (
            ((2, 1, 1, 1, 1), 15, 720),
            ((4, 1, 1), 90, 720),
        ),
    }


def test_one_dicritical_passports_preserve_all_three_filter_stages() -> None:
    passports = one_dicritical_passports()

    assert tuple(
        (
            passport.group_identifier,
            passport.ramification_index,
            passport.tangential_degree,
            passport.jump_defect,
        )
        for passport in passports
    ) == (
        ("6T7", 2, 2, 3),
        ("6T8", 4, 1, 1),
        ("6T12", 2, 2, 3),
        ("6T14", 4, 1, 1),
        ("6T15", 2, 2, 3),
        ("6T15", 3, 1, 2),
        ("6T16", 2, 1, 3),
        ("6T16", 4, 1, 1),
    )
    assert tuple(
        (
            passport.group_identifier,
            passport.ramification_index,
            passport.tangential_degree,
            passport.forced_tangential_jump,
            passport.residual_excess,
        )
        for passport in passports
    ) == (
        ("6T7", 2, 2, 2, 1),
        ("6T8", 4, 1, 0, 1),
        ("6T12", 2, 2, 2, 1),
        ("6T14", 4, 1, 0, 1),
        ("6T15", 2, 2, 2, 1),
        ("6T15", 3, 1, 0, 2),
        ("6T16", 2, 1, 0, 3),
        ("6T16", 4, 1, 0, 1),
    )
    assert tuple(
        (
            passport.group_identifier,
            passport.ramification_index,
            passport.tangential_degree,
        )
        for passport in passports
        if passport.survives_deck_rigidity
    ) == (
        ("6T12", 2, 2),
        ("6T14", 4, 1),
        ("6T15", 2, 2),
        ("6T15", 3, 1),
        ("6T16", 2, 1),
        ("6T16", 4, 1),
    )
    assert tuple(
        (
            passport.group_identifier,
            passport.ramification_index,
            passport.tangential_degree,
        )
        for passport in passports
        if passport.survives_contractible_branch_obstruction
    ) == (
        ("6T15", 3, 1),
        ("6T16", 2, 1),
    )


def test_6t11_needs_more_than_one_allowed_inertia_class() -> None:
    group_6t11 = generated_group(TRANSITIVE_GROUPS[10].generators)
    allowed_classes = tuple(
        analysis
        for analysis in conjugacy_class_analyses(group_6t11)
        if analysis.cycle_type in OREVKOV_ONE_DICRITICAL_TYPES
    )

    assert tuple(
        (
            analysis.cycle_type,
            analysis.class_size,
            analysis.normal_closure_order,
        )
        for analysis in allowed_classes
    ) == (
        ((2, 1, 1, 1, 1), 3, 8),
        ((2, 2, 1, 1), 3, 4),
        ((2, 2, 1, 1), 6, 24),
        ((4, 1, 1), 6, 24),
    )


def test_6t11_full_generating_pairs_have_exact_minimum_costs() -> None:
    pairs = full_generating_inertia_pairs(TRANSITIVE_GROUPS[10])

    assert tuple(
        (
            (
                pair.first.cycle_type,
                pair.first.class_size,
                pair.first.normal_closure_order,
            ),
            (
                pair.second.cycle_type,
                pair.second.class_size,
                pair.second.normal_closure_order,
            ),
            pair.joint_normal_closure_order,
            pair.minimum_ramification_cost,
            tuple(
                realization.components
                for realization in pair.minimum_realizations
            ),
        )
        for pair in pairs
    ) == (
        (
            ((2, 1, 1, 1, 1), 3, 8),
            ((2, 2, 1, 1), 6, 24),
            48,
            6,
            (((2, 1),), ((2, 2),)),
        ),
        (
            ((2, 1, 1, 1, 1), 3, 8),
            ((4, 1, 1), 6, 24),
            48,
            6,
            (((2, 1),), ((4, 1),)),
        ),
        (
            ((2, 2, 1, 1), 6, 24),
            ((4, 1, 1), 6, 24),
            48,
            8,
            (((2, 2),), ((4, 1),)),
        ),
    )


def test_refined_budget_eliminates_6t11_before_the_deck_component() -> None:
    pairs = full_generating_inertia_pairs(TRANSITIVE_GROUPS[10])
    budget_compatible = tuple(
        pair
        for pair in pairs
        if pair.minimum_ramification_cost + 1 <= OREVKOV_SIX_SHEET_BUDGET
    )

    # The cheapest class pair costs 2 + (2*2) = 6 after the forced
    # tangential ramification is counted.  The deck-forced index-one
    # component only strengthens the contradiction.
    assert budget_compatible == ()
    assert tuple(pair.minimum_ramification_cost for pair in pairs) == (6, 6, 8)


def test_double_transposition_component_decompositions_are_complete() -> None:
    realizations = branch_inertia_realizations((2, 2, 1, 1))

    assert tuple(
        (
            realization.components,
            realization.outer_sum_cost,
            realization.forced_tangential_jump,
            realization.ramification_cost,
        )
        for realization in realizations
    ) == (
        (((2, 2),), 2, 2, 4),
        (((2, 1), (2, 1)), 4, 0, 4),
    )


def test_6t7_pre_local_filter_and_6t8_obstruction_profiles_are_exact() -> None:
    group_6t7 = generated_group(TRANSITIVE_GROUPS[6].generators)
    full_6t7 = tuple(
        analysis
        for analysis in conjugacy_class_analyses(group_6t7)
        if analysis.cycle_type in OREVKOV_ONE_DICRITICAL_TYPES
        and analysis.normal_closure_order == len(group_6t7)
    )
    group_6t8 = generated_group(TRANSITIVE_GROUPS[7].generators)
    full_6t8 = tuple(
        analysis
        for analysis in conjugacy_class_analyses(group_6t8)
        if analysis.cycle_type in OREVKOV_ONE_DICRITICAL_TYPES
        and analysis.normal_closure_order == len(group_6t8)
    )

    assert tuple(analysis.cycle_type for analysis in full_6t7) == (
        (2, 2, 1, 1),
    )
    assert tuple(analysis.cycle_type for analysis in full_6t8) == (
        (4, 1, 1),
    )
    one_component_realization = orevkov_budget_profiles(
        required_indices=(2,),
        minimum_unramified_components=1,
    )
    assert tuple(
        profile for profile in one_component_realization if profile.jump_defect >= 2
    ) == (
        # This is the sole 6T7 case left before applying the local-orbit
        # obstruction checked in the next test.
        OrevkovBudgetProfile(
            dicritical_indices=(1, 2),
            jump_defect=2,
        ),
    )
    assert orevkov_budget_profiles(
        required_indices=(2, 2),
        minimum_unramified_components=1,
    ) == (
        # The alternative 6T7 realization has two (e,d)=(2,1)
        # components.  It is saturated and has no finite jump.
        OrevkovBudgetProfile(
            dicritical_indices=(1, 2, 2),
            jump_defect=0,
        ),
    )
    assert orevkov_budget_profiles(
        required_indices=(4,),
        minimum_unramified_components=1,
    ) == (
        # This is the saturated 6T8 profile used by the smooth-line
        # monodromy obstruction.
        OrevkovBudgetProfile(
            dicritical_indices=(1, 4),
            jump_defect=0,
        ),
    )


def test_6t7_residual_jump_has_no_required_local_degree_four_orbit() -> None:
    group_6t7 = generated_group(TRANSITIVE_GROUPS[6].generators)
    meridian = next(
        analysis.representative
        for analysis in conjugacy_class_analyses(group_6t7)
        if analysis.cycle_type == (2, 2, 1, 1)
        and analysis.normal_closure_order == len(group_6t7)
    )
    # Geometry supplies the premise used by this finite certificate: in the
    # saturated residual profile, purity and uniqueness of the ramified
    # divisor make the entire local discriminant the one unibranch germ.
    local_groups = normally_generated_overgroups(group_6t7, meridian)

    assert tuple(
        (analysis.order, analysis.orbit_sizes)
        for analysis in local_groups
    ) == (
        (2, (2, 2, 1, 1)),
        (6, (3, 3)),
        (6, (3, 3)),
        (24, (6,)),
    )
    assert all(4 not in analysis.orbit_sizes for analysis in local_groups)


def test_primitive_e4_jump_local_degree_five_orbit_is_group_theoretically_possible() -> None:
    group_6t14 = generated_group(TRANSITIVE_GROUPS[13].generators)
    meridian_6t14 = next(
        analysis.representative
        for analysis in conjugacy_class_analyses(group_6t14)
        if analysis.cycle_type == (4, 1, 1)
        and analysis.normal_closure_order == len(group_6t14)
    )
    local_groups_6t14 = normally_generated_overgroups(group_6t14, meridian_6t14)

    assert tuple(
        (analysis.order, analysis.orbit_sizes)
        for analysis in local_groups_6t14
        if 5 in analysis.orbit_sizes
    ) == (
        (20, (5, 1)),
        (20, (5, 1)),
    )

    # A direct F20 witness avoids enumerating the much larger S6 subgroup
    # lattice while proving that the same local orbit is possible in 6T16.
    group_6t16 = generated_group(TRANSITIVE_GROUPS[15].generators)
    meridian_6t16 = next(
        analysis.representative
        for analysis in conjugacy_class_analyses(group_6t16)
        if analysis.cycle_type == (4, 1, 1)
        and analysis.normal_closure_order == len(group_6t16)
    )
    five_cycle = (0, 2, 3, 5, 1, 4)
    local_f20 = generated_group((meridian_6t16, five_cycle))

    assert five_cycle in group_6t16
    assert len(local_f20) == 20
    assert permutation_orbit_sizes(local_f20) == (5, 1)
    assert generated_group(
        conjugacy_class(local_f20, meridian_6t16)
    ) == local_f20


def test_one_dicritical_e2d2_passport_forces_one_jump_of_three() -> None:
    one_component_profiles = tuple(
        profile
        for profile in orevkov_budget_profiles(required_indices=(2,))
        if len(profile.dicritical_indices) == 1
    )

    assert one_component_profiles == (
        OrevkovBudgetProfile(
            dicritical_indices=(2,),
            jump_defect=3,
        ),
    )

    # In 6T12=A5, the normally-meridian-generated local groups have no
    # four-sheet orbit.  The only five-sheet possibilities are D10.
    group_6t12 = generated_group(TRANSITIVE_GROUPS[11].generators)
    meridian_6t12 = next(
        analysis.representative
        for analysis in conjugacy_class_analyses(group_6t12)
        if analysis.cycle_type == (2, 2, 1, 1)
        and analysis.normal_closure_order == len(group_6t12)
    )
    local_groups_6t12 = normally_generated_overgroups(group_6t12, meridian_6t12)

    assert all(4 not in analysis.orbit_sizes for analysis in local_groups_6t12)
    assert tuple(
        (analysis.order, analysis.orbit_sizes)
        for analysis in local_groups_6t12
        if 5 in analysis.orbit_sizes
    ) == (
        (10, (5, 1)),
        (10, (5, 1)),
    )

    # The same five-sheet local passport remains abstractly possible in
    # 6T15=A6.  The geometric four-sheet exclusion uses that both 2-cycles
    # coalesce in one orbit; a bare A6 orbit table would be too coarse.
    group_6t15 = generated_group(TRANSITIVE_GROUPS[14].generators)
    meridian_6t15 = next(
        analysis.representative
        for analysis in conjugacy_class_analyses(group_6t15)
        if analysis.cycle_type == (2, 2, 1, 1)
        and analysis.normal_closure_order == len(group_6t15)
    )
    five_cycle = (0, 4, 5, 2, 3, 1)
    local_d10 = generated_group((meridian_6t15, five_cycle))

    assert five_cycle in group_6t15
    assert len(local_d10) == 10
    assert permutation_orbit_sizes(local_d10) == (5, 1)
    assert generated_group(
        conjugacy_class(local_d10, meridian_6t15)
    ) == local_d10

def test_imprimitive_factor_degrees_are_two_by_three_or_three_by_two() -> None:
    rows = classify_six_sheet_groups()
    non_galois = tuple(row for row in rows if not row.is_galois)

    assert {row.identifier: row.factor_degrees for row in non_galois} == {
        "6T3": ((2, 3), (3, 2)),
        "6T4": ((2, 3),),
        "6T5": ((3, 2),),
        "6T6": ((2, 3),),
        "6T7": ((2, 3),),
        "6T8": ((2, 3),),
        "6T9": ((3, 2),),
        "6T10": ((3, 2),),
        "6T11": ((2, 3),),
        "6T12": (),
        "6T13": ((3, 2),),
        "6T14": (),
        "6T15": (),
        "6T16": (),
    }


def test_group_closure_is_independent_of_repeated_generators() -> None:
    rows = classify_six_sheet_groups()
    first_generator = (1, 2, 3, 4, 5, 0)

    assert len(generated_group((first_generator, first_generator))) == rows[0].order
