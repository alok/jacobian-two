"""Exact regression for the split double-contact overlap closure."""

from scripts.a6_delta_ten_split_c22_overlap_closure import (
    SAGE_ARC_LIFT_METADATA,
    exact_double_contact_overlap_closure_certificate,
)


def test_ordinary_two_contact_arc_dominates_the_split_overlap_incidence() -> None:
    certificate = exact_double_contact_overlap_closure_certificate()

    assert certificate.split_equation_residuals == (0, 0, 0)
    assert certificate.dominance_residuals == (0, 0, 0, 0)
    assert certificate.free_parameter_residual == 0
    assert certificate.arc_parameter_derivative != 0
    assert certificate.sage_arc_lift_metadata == SAGE_ARC_LIFT_METADATA
    assert certificate.algebraically_contained
    assert not certificate.topology_computed
    assert certificate.verified
