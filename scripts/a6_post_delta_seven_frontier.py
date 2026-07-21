"""Exact one-pair large-link scan after collision delta seven.

Assuming the conditional one-genuine-pair and one-dicritical ``A6`` setting,
the genus identity leaves finitely many affine degree pairs at each collision
delta.  This module exhausts deltas eight through ten.  The centerless torus
link quotients eliminate every candidate at deltas eight and nine; at delta
ten, only affine degrees ``(4, 9)`` admit an ``A6``-generating representation
whose geometric meridian is a single three-cycle.

This is only the next coarse conditional frontier.  It neither constructs the
surviving curve nor derives the standing hypotheses for a general Keller map.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from scripts.a6_one_pair_infinity import (
    OnePairDegreeCandidate,
    TorusA6Census,
    meridian_exponents,
    one_pair_degree_candidates,
    permutation_from_cycles,
    permutation_power,
    torus_a6_census,
)
from scripts.six_sheet_monodromy import (
    IDENTITY,
    Permutation,
    compose,
    cycle_type,
    generated_group,
)

A6_T49_WITNESS: Final[tuple[Permutation, Permutation]] = (
    permutation_from_cycles((1, 6, 5, 3), (2, 4)),
    permutation_from_cycles((1, 3, 6), (2, 5, 4)),
)


@dataclass(frozen=True, slots=True)
class A6PostDeltaSevenFrontierCertificate:
    """Candidates, exact censuses, and a witness through delta ten."""

    delta_eight_candidates: tuple[OnePairDegreeCandidate, ...]
    delta_nine_candidates: tuple[OnePairDegreeCandidate, ...]
    delta_ten_candidates: tuple[OnePairDegreeCandidate, ...]
    censuses: tuple[TorusA6Census, ...]
    t49_witness: tuple[Permutation, Permutation]

    @property
    def verified(self) -> bool:
        """Whether the exact scan has the expected unique next survivor."""

        x, y = self.t49_witness
        u, v = meridian_exponents(4, 9)
        meridian = compose(
            permutation_power(x, u),
            permutation_power(y, v),
        )
        return bool(
            tuple(c.affine_degrees for c in self.delta_eight_candidates)
            == ((2, 21), (3, 11))
            and tuple(c.projective_pair for c in self.delta_eight_candidates)
            == ((19, 21), (8, 11))
            and tuple(c.affine_degrees for c in self.delta_nine_candidates)
            == ((2, 23),)
            and tuple(c.projective_pair for c in self.delta_nine_candidates)
            == ((21, 23),)
            and tuple(c.affine_degrees for c in self.delta_ten_candidates)
            == ((2, 25), (3, 13), (4, 9), (5, 7))
            and tuple(c.projective_pair for c in self.delta_ten_candidates)
            == ((23, 25), (10, 13), (5, 9), (2, 7))
            and tuple(c.pair for c in self.censuses)
            == (
                (2, 21),
                (3, 11),
                (2, 23),
                (2, 25),
                (3, 13),
                (4, 9),
                (5, 7),
            )
            and tuple(c.meridian_exponents for c in self.censuses)
            == (
                (1, -10),
                (-1, 4),
                (1, -11),
                (1, -12),
                (1, -4),
                (1, -2),
                (-2, 3),
            )
            and tuple((c.x_choices, c.y_choices) for c in self.censuses)
            == (
                (46, 81),
                (81, 1),
                (46, 1),
                (46, 145),
                (81, 1),
                (136, 81),
                (145, 1),
            )
            and tuple(c.generated_order_histogram for c in self.censuses)
            == (
                (
                    (1, 1),
                    (2, 45),
                    (3, 80),
                    (6, 720),
                    (12, 720),
                    (24, 720),
                    (60, 1440),
                ),
                ((1, 1), (3, 80)),
                ((1, 1), (2, 45)),
                (
                    (1, 1),
                    (2, 45),
                    (5, 144),
                    (10, 720),
                    (60, 2880),
                    (360, 2880),
                ),
                ((1, 1), (3, 80)),
                (
                    (1, 1),
                    (2, 45),
                    (3, 80),
                    (4, 90),
                    (6, 720),
                    (12, 720),
                    (24, 2160),
                    (36, 1440),
                    (60, 1440),
                    (360, 4320),
                ),
                ((1, 1), (5, 144)),
            )
            and tuple(c.a6_meridian_histogram for c in self.censuses)
            == (
                (),
                (),
                (),
                (((4, 2), 1440), ((5, 1), 1440)),
                (),
                (
                    ((3, 1, 1, 1), 720),
                    ((3, 3), 720),
                    ((5, 1), 2880),
                ),
                (),
            )
            and tuple(c.single_three_cycle_a6_pairs for c in self.censuses)
            == (0, 0, 0, 0, 0, 720, 0)
            and permutation_power(x, 4) == IDENTITY
            and permutation_power(y, 9) == IDENTITY
            and len(generated_group((x, y))) == 360
            and cycle_type(meridian) == (3, 1, 1, 1)
        )


@cache
def exact_post_delta_seven_frontier_certificate(
) -> A6PostDeltaSevenFrontierCertificate:
    """Construct the exact candidate and torus-quotient scan."""

    pairs = ((2, 21), (3, 11), (2, 23), (2, 25), (3, 13), (4, 9), (5, 7))
    return A6PostDeltaSevenFrontierCertificate(
        delta_eight_candidates=one_pair_degree_candidates(8),
        delta_nine_candidates=one_pair_degree_candidates(9),
        delta_ten_candidates=one_pair_degree_candidates(10),
        censuses=tuple(torus_a6_census(*pair) for pair in pairs),
        t49_witness=A6_T49_WITNESS,
    )


def main() -> None:
    """Print the next coarse conditional frontier."""

    certificate = exact_post_delta_seven_frontier_certificate()
    print("post-delta-seven frontier verified:", certificate.verified)
    print(
        "single-three-cycle A6 counts:",
        {
            census.pair: census.single_three_cycle_a6_pairs
            for census in certificate.censuses
        },
    )
    print("next survivor: Delta=10, affine (4,9), infinity (5,9)")


if __name__ == "__main__":
    main()
