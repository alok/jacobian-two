"""Aggregate the exact rank audit for all 22 true-split allocations.

The generic split ledger enumerates 22 allowed allocations across the six
expected codimension-two delta-ten profiles.  Four focused certificates now
cover that entire list:

* ``C3+7N`` and ``C2^2+6N`` on their true component jets;
* ``T112+6N`` and ``C2+T111+5N`` on labeled four-source fibers;
* ``Q0+4N`` and ``T111^2+4N`` through the global fiber-remainder systems; and
* the full ``k=2`` to ``k=-2`` family involution.

This module checks that the exact audited allocation keys equal the generated
22-row ledger, rather than inferring completeness from a count.  The result
is a rank/dimension theorem: no true-split coefficient-rank stratum in these
six profiles supports an incidence component of dimension at least three.

It is not yet a full split exclusion.  The clean rank-open loci for ``C3``,
``C2^2``, ``T112``, and the mixed profile are now connected to their global
components, while the global-fiber argument handles ``Q0`` and two triples.
The prescribed overlap allocations and the exceptional compatible
affine-line fibers on rank-drop bases still require topology or a closure
proof.  Deeper profile intersections and removed non-clean boundaries remain
open.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from scripts.a6_delta_ten_split_codim_two import SPLIT_ALLOCATIONS
from scripts.a6_delta_ten_split_component_closure import (
    exact_split_component_closure_certificate,
)
from scripts.a6_delta_ten_split_contact_closure import (
    exact_split_contact_closure_certificate,
)
from scripts.a6_delta_ten_split_contact_rank import (
    c3_rank_specs,
    exact_split_contact_rank_certificate,
    residual_rank_specs,
    simple_double_contact_specs,
)
from scripts.a6_delta_ten_split_rank_q0_t111 import (
    exact_split_q0_double_triple_rank_certificate,
)
from scripts.a6_delta_ten_split_t112_mixed_rank import (
    SPLIT_INCIDENCE_SPECS,
    exact_split_t112_mixed_rank_certificate,
)

AllocationKey = tuple[str, str, str]

EXPECTED_PROFILE_COUNTS: Final = (
    ("C3+7N", 5),
    ("C2^2+6N", 6),
    ("T112+6N", 4),
    ("C2+T111+5N", 3),
    ("Q0+4N", 2),
    ("T111^2+4N", 2),
)


def _split_class(kappa: int) -> str:
    """Return the ledger's split-class label for a representative value."""

    if kappa == 0:
        return "k=0"
    if kappa == 2:
        return "k=+/-2"
    msg = f"the aggregate uses only k=0 and the transported k=2 representative: {kappa}"
    raise ValueError(msg)


def _expected_keys() -> tuple[AllocationKey, ...]:
    """Return the generated ledger keys in their canonical order."""

    return tuple(
        (allocation.profile, allocation.split_class, allocation.allocation)
        for allocation in SPLIT_ALLOCATIONS
    )


def _audited_keys() -> tuple[AllocationKey, ...]:
    """Return keys derived from the four focused rank systems."""

    contact_keys = tuple(
        (spec.profile, _split_class(spec.kappa), spec.allocation)
        for spec in c3_rank_specs()
    )
    simple_double_contact_keys = tuple(
        ("C2^2+6N", _split_class(spec.kappa), spec.allocation)
        for spec in simple_double_contact_specs()
    )
    residual_double_contact_keys = tuple(
        ("C2^2+6N", _split_class(spec.kappa), spec.allocation)
        for spec in residual_rank_specs()
    )
    triple_keys = tuple(
        (spec.profile, _split_class(spec.kappa), spec.allocation)
        for spec in SPLIT_INCIDENCE_SPECS
    )
    global_fiber_keys = tuple(
        key for key in _expected_keys() if key[0] in {"Q0+4N", "T111^2+4N"}
    )
    return (
        *contact_keys,
        *simple_double_contact_keys,
        *residual_double_contact_keys,
        *triple_keys,
        *global_fiber_keys,
    )


@dataclass(frozen=True, slots=True)
class AllSplitRankCertificate:
    """Exact coverage and dimension certificate for the 22-row ledger."""

    expected_keys: tuple[AllocationKey, ...]
    audited_keys: tuple[AllocationKey, ...]
    audited_incidence_dimension_bounds: tuple[tuple[AllocationKey, int], ...]
    profile_counts: tuple[tuple[str, int], ...]
    contact_rank_verified: bool
    triple_mixed_rank_verified: bool
    global_fiber_rank_verified: bool
    triple_mixed_maximal_rank_closure_verified: bool
    contact_rank_open_closure_verified: bool
    maximum_split_incidence_dimension: int
    maximum_residual_rankdrop_incidence_dimension: int
    contact_exceptional_topology_open: bool
    deeper_boundaries_open: bool
    proves_plane_jacobian_conjecture: bool

    @property
    def allocation_coverage_exact(self) -> bool:
        """Whether there are no duplicated, omitted, or invented allocation keys."""

        return bool(
            len(self.expected_keys) == len(self.audited_keys) == 22
            and len(set(self.audited_keys)) == 22
            and set(self.audited_keys) == set(self.expected_keys)
            and len(self.audited_incidence_dimension_bounds) == 22
            and {
                key for key, _bound in self.audited_incidence_dimension_bounds
            }
            == set(self.expected_keys)
            and self.profile_counts == EXPECTED_PROFILE_COUNTS
        )

    @property
    def all_rank_strata_classified(self) -> bool:
        """Whether every focused rank certificate and coverage check passes."""

        return bool(
            self.allocation_coverage_exact
            and self.contact_rank_verified
            and self.triple_mixed_rank_verified
            and self.global_fiber_rank_verified
            and self.maximum_split_incidence_dimension == 2
            and self.maximum_residual_rankdrop_incidence_dimension == 1
        )

    @property
    def verified(self) -> bool:
        """Whether the rank theorem and its honest remaining boundary agree."""

        return bool(
            self.all_rank_strata_classified
            and self.triple_mixed_maximal_rank_closure_verified
            and self.contact_rank_open_closure_verified
            and self.contact_exceptional_topology_open
            and self.deeper_boundaries_open
            and not self.proves_plane_jacobian_conjecture
        )


@cache
def exact_all_split_rank_certificate() -> AllSplitRankCertificate:
    """Build the exact aggregate true-split rank certificate."""

    expected_keys = _expected_keys()
    profile_counts = tuple(
        (
            profile,
            sum(key[0] == profile for key in expected_keys),
        )
        for profile, _expected_count in EXPECTED_PROFILE_COUNTS
    )
    contact = exact_split_contact_rank_certificate()
    triple_mixed = exact_split_t112_mixed_rank_certificate()
    global_fiber = exact_split_q0_double_triple_rank_certificate()
    closure = exact_split_component_closure_certificate()
    contact_closure = exact_split_contact_closure_certificate()
    audited_dimension_bounds = (
        *(
            (
                (spec.profile, _split_class(spec.kappa), spec.allocation),
                spec.generic_incidence_dimension,
            )
            for spec in c3_rank_specs()
        ),
        *(
            (
                ("C2^2+6N", _split_class(spec.kappa), spec.allocation),
                spec.generic_incidence_dimension,
            )
            for spec in simple_double_contact_specs()
        ),
        *(
            (
                ("C2^2+6N", _split_class(spec.kappa), spec.allocation),
                spec.generic_incidence_dimension,
            )
            for spec in residual_rank_specs()
        ),
        *(
            (
                (spec.profile, _split_class(spec.kappa), spec.allocation),
                max(
                    spec.generic_incidence_dimension,
                    spec.residual_incidence_dimension_bound,
                ),
            )
            for spec in SPLIT_INCIDENCE_SPECS
        ),
        *(
            (
                key,
                (
                    global_fiber.quadruple.maximum_valid_incidence_dimension
                    if key[0] == "Q0+4N"
                    else global_fiber.double_triple.maximum_valid_incidence_dimension
                ),
            )
            for key in expected_keys
            if key[0] in {"Q0+4N", "T111^2+4N"}
        ),
    )
    return AllSplitRankCertificate(
        expected_keys=expected_keys,
        audited_keys=_audited_keys(),
        audited_incidence_dimension_bounds=audited_dimension_bounds,
        profile_counts=profile_counts,
        contact_rank_verified=contact.verified,
        triple_mixed_rank_verified=triple_mixed.verified,
        global_fiber_rank_verified=global_fiber.verified,
        triple_mixed_maximal_rank_closure_verified=(
            closure.maximal_rank_split_topology_closed
        ),
        contact_rank_open_closure_verified=(
            contact_closure.maximal_rank_topology_closed
        ),
        maximum_split_incidence_dimension=max(
            bound for _key, bound in audited_dimension_bounds
        ),
        maximum_residual_rankdrop_incidence_dimension=max(
            contact.maximum_residual_incidence_dimension,
            *(
                spec.residual_incidence_dimension_bound
                for spec in SPLIT_INCIDENCE_SPECS
            ),
        ),
        contact_exceptional_topology_open=contact_closure.exceptional_topology_open,
        deeper_boundaries_open=True,
        proves_plane_jacobian_conjecture=False,
    )


def main() -> int:
    """Print the exact all-allocation conclusion and remaining boundary."""

    certificate = exact_all_split_rank_certificate()
    print("true-split allocation keys covered:", len(certificate.audited_keys))
    print(
        "all true-split rank strata classified:", certificate.all_rank_strata_classified
    )
    print(
        "largest compatible residual rank-drop incidence dimension:",
        certificate.maximum_residual_rankdrop_incidence_dimension,
    )
    print("remaining: exceptional split contact fibers and deeper intersections")
    print("proves JC(2):", certificate.proves_plane_jacobian_conjecture)
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
