"""Optional Sage/GAP cross-check for the pure-Python degree-six certificate.

Run from the repository root with:

    sage tools/check_six_sheet_gap.sage

Sage is intentionally not a project or CI dependency.  This script compares
the dependency-free permutation engine against GAP's transitive-groups library
and independently verifies every invariant used by the monodromy filters.
"""

from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from sage.libs.gap.libgap import libgap

from scripts.six_sheet_monodromy import (
    TRANSITIVE_GROUPS,
    classify_six_sheet_groups,
    conjugacy_class_analyses,
    generated_group,
    normally_generated_overgroups,
)


def gap_conjugacy_data(group):
    """Return exact ``(cycle type, class size, closure order)`` triples."""

    points = libgap([1, 2, 3, 4, 5, 6])
    result = []
    for conjugacy_class in group.ConjugacyClasses():
        representative = conjugacy_class.Representative()
        cycle_type = tuple(
            sorted(
                (int(length) for length in representative.CycleLengths(points)),
                reverse=True,
            )
        )
        cyclic_subgroup = libgap.Subgroup(group, [representative])
        closure = libgap.NormalClosure(group, cyclic_subgroup)
        result.append(
            (cycle_type, int(conjugacy_class.Size()), int(closure.Size()))
        )
    return sorted(result)


def custom_conjugacy_data(fixture):
    """Return the corresponding triples from the pure-Python engine."""

    group = generated_group(fixture.generators)
    return sorted(
        (
            analysis.cycle_type,
            analysis.class_size,
            analysis.normal_closure_order,
        )
        for analysis in conjugacy_class_analyses(group)
    )


def gap_local_overgroup_data(group, representative):
    """Return normally generated overgroup orders and sheet-orbit sizes."""

    cyclic_subgroup = libgap.Subgroup(group, [representative])
    intermediate = libgap.IntermediateSubgroups(group, cyclic_subgroup)
    overgroups = [cyclic_subgroup, *intermediate["subgroups"], group]
    points = libgap([1, 2, 3, 4, 5, 6])
    result = []
    for subgroup in overgroups:
        normal_closure = libgap.NormalClosure(
            subgroup,
            libgap.Subgroup(subgroup, [representative]),
        )
        if int(normal_closure.Size()) != int(subgroup.Size()):
            continue
        orbit_sizes = tuple(
            sorted(
                (len(orbit) for orbit in subgroup.Orbits(points).sage()),
                reverse=True,
            )
        )
        result.append((int(subgroup.Size()), orbit_sizes))
    return sorted(result)


def custom_6t7_local_overgroup_data():
    """Return the pure-Python local-monodromy obstruction data for 6T7."""

    fixture = TRANSITIVE_GROUPS[6]
    group = generated_group(fixture.generators)
    meridian = next(
        analysis.representative
        for analysis in conjugacy_class_analyses(group)
        if analysis.cycle_type == (2, 2, 1, 1)
        and analysis.normal_closure_order == len(group)
    )
    return sorted(
        (analysis.order, analysis.orbit_sizes)
        for analysis in normally_generated_overgroups(group, meridian)
    )


def main():
    rows = classify_six_sheet_groups()
    number_of_groups = int(libgap.NrTransitiveGroups(6))
    assert number_of_groups == len(rows) == len(TRANSITIVE_GROUPS) == 16

    for index, (row, fixture) in enumerate(
        zip(rows, TRANSITIVE_GROUPS, strict=True),
        start=1,
    ):
        group = libgap.TransitiveGroup(6, index)
        stabilizer = group.Stabilizer(1)
        normalizer = group.Normalizer(stabilizer)
        quotient = normalizer.FactorGroup(stabilizer)
        block_sizes = tuple(
            sorted({len(block) for block in group.AllBlocks().sage()})
        )
        fixed_point_elements = [
            element
            for element in group.Elements()
            if int(element.NrMovedPoints()) < 6
        ]
        fixed_point_subgroup = libgap.Subgroup(group, fixed_point_elements)
        fixed_point_closure = libgap.NormalClosure(group, fixed_point_subgroup)

        assert row.identifier == f"6T{index}"
        assert row.order == int(group.Size())
        assert row.structure == str(group.StructureDescription())
        assert row.stabilizer_order == int(stabilizer.Size())
        assert row.normalizer_order == int(normalizer.Size())
        assert row.automorphism_group == str(quotient.StructureDescription())
        assert row.is_galois == bool(group.IsNormal(stabilizer))
        assert row.is_primitive == bool(group.IsPrimitive())
        assert row.block_sizes == block_sizes
        assert row.fixed_point_normal_closure_order == int(
            fixed_point_closure.Size()
        )
        assert custom_conjugacy_data(fixture) == gap_conjugacy_data(group)

    group_6t7 = libgap.TransitiveGroup(6, 7)
    points = libgap([1, 2, 3, 4, 5, 6])
    meridian_6t7 = next(
        conjugacy_class.Representative()
        for conjugacy_class in group_6t7.ConjugacyClasses()
        if tuple(
            sorted(
                (
                    int(length)
                    for length in conjugacy_class.Representative().CycleLengths(
                        points
                    )
                ),
                reverse=True,
            )
        )
        == (2, 2, 1, 1)
        and int(
            libgap.NormalClosure(
                group_6t7,
                libgap.Subgroup(
                    group_6t7,
                    [conjugacy_class.Representative()],
                ),
            ).Size()
        )
        == 24
    )
    assert custom_6t7_local_overgroup_data() == gap_local_overgroup_data(
        group_6t7,
        meridian_6t7,
    )

    gap_version = libgap.eval("GAPInfo.Version")
    print(
        f"independent GAP {gap_version} cross-check: "
        "16/16 groups, every conjugacy class, and the 6T7 local "
        "overgroups match"
    )


if __name__ == "__main__":
    main()
