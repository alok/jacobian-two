"""Exhaust the conditional residual ``A6`` delta-five family.

The generic and codimension-one analysis in
``scripts.a6_delta_five_family`` leaves, on the valid parameter set, the
union of the singular locus of the collision discriminant and its
intersection with the triple-collision divisor.  This module gives exact
rational normalizations of every irreducible residual curve, factors every
possible further degeneration, and replays the affine-complement
presentations on every generic and exceptional topology stratum.

The conclusion is conditional and computer-assisted.  It assumes the
one-genuine-pair and finite-singularity hypotheses of
``scripts.a6_one_pair_infinity`` and uses Sage-generated Zariski--van Kamp
presentations.  It does not prove those hypotheses for a Keller branch and
does not prove the plane Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final, cast

from sympy import (
    Expr,
    Poly,
    Rational,
    Symbol,
    cancel,
    diff,
    discriminant,
    expand,
    factor,
    gcd,
    sqrt,
)
from sympy.polys.domains import QQ

from scripts.a6_delta_five_family import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    CUSP_COLLISION_FACTOR,
    EXTRA_CRITICAL_FACTOR,
    GAMMA,
    S,
    TRIPLE_COLLISION_FACTOR,
)
from scripts.a6_one_pair_infinity import (
    OnePairDegreeCandidate,
    ThreeCyclePresentationCensus,
    TorusA6Census,
    one_pair_degree_candidates,
    three_cycle_presentation_census,
    torus_a6_census,
)

Relations = tuple[tuple[int, ...], ...]

H_PARAMETER: Final = Symbol("h")
M_PARAMETER: Final = Symbol("m")
N_PARAMETER: Final = Symbol("u")
E_PARAMETER: Final = Symbol("v")

DEGREE_SEVEN_DISCRIMINANT_FACTOR: Final = cancel(
    discriminant(COLLISION_POLYNOMIAL, S) / CUSP_COLLISION_FACTOR
)
SINGULAR_DISCRIMINANT_GENERATORS: Final = (
    DEGREE_SEVEN_DISCRIMINANT_FACTOR,
    diff(DEGREE_SEVEN_DISCRIMINANT_FACTOR, ALPHA),
    diff(DEGREE_SEVEN_DISCRIMINANT_FACTOR, BETA),
    diff(DEGREE_SEVEN_DISCRIMINANT_FACTOR, GAMMA),
)


def _parameters(alpha: Expr, beta: Expr, gamma: Expr) -> dict[Symbol, Expr]:
    """Return a substitution for the normalized three-parameter family."""

    return {ALPHA: alpha, BETA: beta, GAMMA: gamma}


def _valid_product(alpha: Expr, beta: Expr, gamma: Expr) -> Expr:
    """Evaluate ``alpha*L*C`` at one parameter point."""

    return expand(
        (
            ALPHA * EXTRA_CRITICAL_FACTOR * CUSP_COLLISION_FACTOR
        ).subs(_parameters(alpha, beta, gamma))
    )


def _rational_pullback(
    polynomial: Expr,
    parameter: Symbol,
    alpha: Expr,
    beta: Expr,
    gamma: Expr,
) -> Expr:
    """Evaluate a parameter polynomial in ``QQ(parameter)`` exactly.

    Using SymPy's rational-function domain avoids expensive general-purpose
    expression cancellation.  The returned expression is zero exactly when
    the rational parametrization lies on the polynomial's zero set.
    """

    field = QQ.frac_field(parameter)
    values = tuple(
        field.from_sympy(value) for value in (alpha, beta, gamma)
    )
    result = field.zero
    for exponents, coefficient in Poly(
        polynomial,
        ALPHA,
        BETA,
        GAMMA,
        domain=QQ,
    ).terms():
        term = field.convert(coefficient)
        for value, exponent in zip(values, exponents, strict=True):
            term *= value**exponent
        result += term
    return cast(Expr, result.as_expr())


# A: one triple root of H, hence one contact-three collision generically.
A_ALPHA: Final = cancel(
    -2
    * H_PARAMETER**3
    * (3 * H_PARAMETER**2 + 12 * H_PARAMETER + 13)
    / (3 * H_PARAMETER + 11)
)
A_BETA: Final = cancel(
    -(
        3 * H_PARAMETER**5
        + 9 * H_PARAMETER**4
        - 4 * H_PARAMETER**3
        - 24 * H_PARAMETER**2
        + 6 * H_PARAMETER
        - 4
    )
    / (3 * H_PARAMETER + 11)
)
A_GAMMA: Final = cancel(
    -6
    * (H_PARAMETER**2 + 3 * H_PARAMETER - 2)
    / (3 * H_PARAMETER + 11)
)
A_QUADRATIC_NUMERATOR: Final = (
    (3 * H_PARAMETER + 11) * S**2
    + (3 * H_PARAMETER**2 + 21 * H_PARAMETER + 34) * S
    + 6 * H_PARAMETER**2
    + 24 * H_PARAMETER
    + 26
)
A_FACTORIZATION: Final = cancel(
    (S - H_PARAMETER) ** 3
    * A_QUADRATIC_NUMERATOR
    / (3 * H_PARAMETER + 11)
)
A_L_FORMULA: Final = cancel(
    -2
    * (3 * H_PARAMETER + 1) ** 2
    * (3 * H_PARAMETER + 4) ** 3
    / (3 * H_PARAMETER + 11)
)
A_C_FORMULA: Final = cancel(
    -3 * (H_PARAMETER + 1) ** 5 / (3 * H_PARAMETER + 11)
)
A_T_FORMULA: Final = cancel(
    (H_PARAMETER + 3) ** 3
    * (3 * H_PARAMETER + 1) ** 3
    / (3 * H_PARAMETER + 11) ** 2
)
A_QUADRATIC_AT_TRIPLE_ROOT: Final = cancel(
    2
    * (H_PARAMETER + 1) ** 2
    * (3 * H_PARAMETER + 13)
    / (3 * H_PARAMETER + 11)
)
A_QUADRATIC_DISCRIMINANT: Final = (
    3
    * (H_PARAMETER + 1) ** 2
    * (3 * H_PARAMETER**2 + 12 * H_PARAMETER + 4)
)


# B: two double roots of H, hence two contact-two collisions generically.
B_P: Final = cancel(
    -2
    * (M_PARAMETER**2 - 4 * M_PARAMETER + 5)
    / (M_PARAMETER - 3)
)
B_Q: Final = cancel(
    -(
        2 * M_PARAMETER**3
        - 6 * M_PARAMETER**2
        + 3 * M_PARAMETER
        + 3
    )
    / (M_PARAMETER - 3)
)
B_R: Final = cancel(-2 * (M_PARAMETER - 2) / (M_PARAMETER - 3))
B_ALPHA: Final = cancel(
    2
    * (M_PARAMETER - 2)
    * (
        2 * M_PARAMETER**3
        - 6 * M_PARAMETER**2
        + 3 * M_PARAMETER
        + 3
    )
    ** 2
    / (M_PARAMETER - 3) ** 3
)
B_BETA: Final = cancel(
    (
        4 * M_PARAMETER**7
        - 28 * M_PARAMETER**6
        + 56 * M_PARAMETER**5
        + 36 * M_PARAMETER**4
        - 259 * M_PARAMETER**3
        + 247 * M_PARAMETER**2
        + 75 * M_PARAMETER
        - 155
    )
    / (M_PARAMETER - 3) ** 3
)
B_GAMMA: Final = cancel(
    -2
    * (2 * M_PARAMETER**2 - 8 * M_PARAMETER + 9)
    / (M_PARAMETER - 3)
)
B_QUADRATIC: Final = S**2 + B_P * S + B_Q
B_FACTORIZATION: Final = cancel(B_QUADRATIC**2 * (S - B_R))
B_C_FORMULA: Final = cancel(
    4
    * (M_PARAMETER - 2) ** 2
    * (M_PARAMETER - 1) ** 5
    / (M_PARAMETER - 3) ** 3
)
B_L_FORMULA: Final = cancel(
    2
    * M_PARAMETER
    * (3 * M_PARAMETER - 5) ** 2
    * (6 * M_PARAMETER**2 - 16 * M_PARAMETER + 9) ** 2
    / (M_PARAMETER - 3) ** 3
)
B_T_FORMULA: Final = cancel(
    4
    * M_PARAMETER
    * (M_PARAMETER - 2) ** 5
    * (3 * M_PARAMETER - 5) ** 2
    / (M_PARAMETER - 3) ** 4
)
B_QUADRATIC_DISCRIMINANT: Final = cancel(
    4
    * (M_PARAMETER - 2)
    * (M_PARAMETER - 1) ** 2
    * (3 * M_PARAMETER - 8)
    / (M_PARAMETER - 3) ** 2
)
B_SIMPLE_ROOT_VALUE: Final = cancel(
    -(M_PARAMETER - 1) ** 2
    * (2 * M_PARAMETER**2 - 12 * M_PARAMETER + 15)
    / (M_PARAMETER - 3) ** 2
)


# N: a nonordinary triple collision with one tangent pair, plus one node.
N_ALPHA: Final = cancel(
    (N_PARAMETER + 1) ** 2
    * (N_PARAMETER + 2) ** 2
    * (2 * N_PARAMETER + 1)
    / N_PARAMETER
)
N_BETA: Final = cancel(
    (N_PARAMETER**2 + 3 * N_PARAMETER + 1)
    * (
        N_PARAMETER**3
        + 4 * N_PARAMETER**2
        + 5 * N_PARAMETER
        + 3
    )
    / N_PARAMETER
)
N_GAMMA: Final = cancel(
    (N_PARAMETER**2 + 4 * N_PARAMETER + 1) / N_PARAMETER
)
N_DOUBLE_ROOT: Final = -N_PARAMETER - 2
N_NODE_ROOT: Final = cancel(-2 - 1 / N_PARAMETER)
N_QUADRATIC: Final = (
    S**2
    - S * N_PARAMETER
    + N_PARAMETER**2
    + 2 * N_PARAMETER
    + 1
)
N_FACTORIZATION: Final = cancel(
    (S + N_PARAMETER + 2) ** 2
    * (S * N_PARAMETER + 2 * N_PARAMETER + 1)
    * N_QUADRATIC
    / N_PARAMETER
)
N_C_FORMULA: Final = cancel(
    (N_PARAMETER + 1) ** 4 * (N_PARAMETER + 2) / N_PARAMETER
)
N_L_FORMULA: Final = cancel(
    (2 * N_PARAMETER + 3)
    * (3 * N_PARAMETER + 2) ** 2
    * (3 * N_PARAMETER + 5) ** 2
    / N_PARAMETER
)
N_QUADRATIC_AT_DOUBLE_ROOT: Final = (
    (N_PARAMETER + 1) * (3 * N_PARAMETER + 5)
)
N_QUADRATIC_AT_NODE_ROOT: Final = cancel(
    (N_PARAMETER + 1) ** 4 / N_PARAMETER**2
)
N_QUADRATIC_DISCRIMINANT: Final = (
    -(N_PARAMETER + 2) * (3 * N_PARAMETER + 2)
)
N_ROOT_COINCIDENCE_FORMULA: Final = cancel(
    (1 - N_PARAMETER**2) / N_PARAMETER
)


# E: an ordinary triple collision and a separate contact-two collision.
E_ALPHA: Final = Rational(9, 5) * E_PARAMETER - 9
E_BETA: Final = E_PARAMETER
E_GAMMA: Final = Rational(6)
E_TRIPLE_FACTOR: Final = (
    E_PARAMETER
    + 5 * S**3
    + 10 * S**2
    + 5 * S
    - 5
) / 5
E_FACTORIZATION: Final = (S + 3) ** 2 * E_TRIPLE_FACTOR
E_C_FORMULA: Final = 4 * (E_PARAMETER - 5) / 5
E_L_FORMULA: Final = 5 * (27 * E_PARAMETER - 155)
E_TRIPLE_DISCRIMINANT: Final = (
    -25 * (E_PARAMETER - 5) * (27 * E_PARAMETER - 155)
)
E_TRIPLE_AT_CONTACT_ROOT: Final = (E_PARAMETER - 65) / 5

# Factors at which a residual parametrization is invalid, respectively can
# change collision topology.  Removing the polynomial gcd leaves exactly the
# valid exceptional parameters audited below.
A_INVALID_FACTOR: Final = (
    (3 * H_PARAMETER + 11)
    * H_PARAMETER
    * (3 * H_PARAMETER**2 + 12 * H_PARAMETER + 13)
    * (3 * H_PARAMETER + 1)
    * (3 * H_PARAMETER + 4)
    * (H_PARAMETER + 1)
)
A_CHANGE_FACTOR: Final = (
    (H_PARAMETER + 3)
    * (3 * H_PARAMETER + 1)
    * (H_PARAMETER + 1)
    * (3 * H_PARAMETER + 13)
    * (3 * H_PARAMETER**2 + 12 * H_PARAMETER + 4)
)
A_EXCEPTION_FACTOR: Final = (
    (H_PARAMETER + 3)
    * (3 * H_PARAMETER + 13)
    * (3 * H_PARAMETER**2 + 12 * H_PARAMETER + 4)
)
B_INVALID_FACTOR: Final = (
    (M_PARAMETER - 3)
    * M_PARAMETER
    * (M_PARAMETER - 2)
    * (M_PARAMETER - 1)
    * (3 * M_PARAMETER - 5)
    * (6 * M_PARAMETER**2 - 16 * M_PARAMETER + 9)
    * (2 * M_PARAMETER**3 - 6 * M_PARAMETER**2 + 3 * M_PARAMETER + 3)
)
B_CHANGE_FACTOR: Final = (
    M_PARAMETER
    * (M_PARAMETER - 2)
    * (3 * M_PARAMETER - 5)
    * (M_PARAMETER - 1)
    * (3 * M_PARAMETER - 8)
    * (2 * M_PARAMETER**2 - 12 * M_PARAMETER + 15)
)
B_EXCEPTION_FACTOR: Final = (
    (3 * M_PARAMETER - 8)
    * (2 * M_PARAMETER**2 - 12 * M_PARAMETER + 15)
)
N_INVALID_FACTOR: Final = (
    N_PARAMETER
    * (N_PARAMETER + 1)
    * (N_PARAMETER + 2)
    * (2 * N_PARAMETER + 1)
    * (2 * N_PARAMETER + 3)
    * (3 * N_PARAMETER + 2)
    * (3 * N_PARAMETER + 5)
)
N_CHANGE_FACTOR: Final = (
    (N_PARAMETER - 1)
    * (N_PARAMETER + 1)
    * (3 * N_PARAMETER + 5)
    * (N_PARAMETER + 2)
    * (3 * N_PARAMETER + 2)
)
N_EXCEPTION_FACTOR: Final = N_PARAMETER - 1
E_INVALID_FACTOR: Final = (
    (E_PARAMETER - 5) * (27 * E_PARAMETER - 155)
)
E_CHANGE_FACTOR: Final = E_INVALID_FACTOR * (E_PARAMETER - 65)
E_EXCEPTION_FACTOR: Final = E_PARAMETER - 65


A_SAMPLE: Final = (
    Rational(16, 5),
    Rational(32, 5),
    Rational(24, 5),
)
B_SAMPLE: Final = (
    Rational(1, 6),
    Rational(127, 108),
    Rational(2),
)
N_SAMPLE: Final = (Rational(20, 3), Rational(1), Rational(2, 3))
E_SAMPLE: Final = (Rational(-9), Rational(0), Rational(6))

P4_PARAMETERS: Final = (
    Rational(-114244, 81),
    Rational(-63407, 81),
    Rational(34, 3),
)
PNT_PARAMETERS: Final = (Rational(108), Rational(65), Rational(6))

SQRT_SIX: Final = sqrt(6)
P32_PLUS_H: Final = -2 + Rational(2, 3) * SQRT_SIX
P32_MINUS_H: Final = -2 - Rational(2, 3) * SQRT_SIX


def _p32_parameters(h: Expr) -> tuple[Expr, Expr, Expr]:
    """Return the quadratic contact-three/contact-two parameter point."""

    return (
        expand(-24 * (109 * h + 40)),
        expand(-Rational(2, 3) * (2189 * h + 802)),
        expand(2 * (3 * h + 2)),
    )


P32_PLUS_PARAMETERS: Final = _p32_parameters(P32_PLUS_H)
P32_MINUS_PARAMETERS: Final = _p32_parameters(P32_MINUS_H)


# Raw Sage 10.8 affine Zariski--van Kamp relations.  Generators 1..3 are
# geometric fiber meridians.  The Sage reproducer regenerates these tuples.
A_RELATIONS: Final[Relations] = (
    (-3, -2, 3, 1, -3, 2, 3, 1, -3, 2, 3, 1, -3, -2, 3, -1, -3, -2, 3, -1),
    (-3, -2, 3, -1, -3, -2, -3, 2, 3, 1, -3, 2, 3, 1, -3, -2, 3, -1, -3, -2, 3, 2, 3, 1, -3, 2, 3, -1),
    (-3, -2, 3, -1, -3, -2, 3, 2, 3, 1),
    (-3, -2, -3, 2, 3, 1, -3, -2, 3, 2, 3, -1),
    (3, 2, 3, 2, 3, 2, -3, -2, -3, -2, -3, -2),
)
B_RELATIONS: Final[Relations] = (
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -1, 2, 3, 2, -3, -2, 1, 2, 3, -2),
    (-3, -2, 1, 2, 3, -2, 1, 2, 3, -2, -1, 2, -3, -2, -1, 2),
    (2, 1, 2, 1, -2, -1, -2, -1),
    (-3, 1),
)
N_RELATIONS: Final[Relations] = (
    (2, 3, 1, -3, -2, -1),
    (2, 3, 2, -1, -3, -2, -3, -2, 3, 1),
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -3, -2, 1, 2, 3, 2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
)
E_RELATIONS: Final[Relations] = (
    (2, 3, 1, -3, -2, -1),
    (2, -1, -3, -2, 3, 1),
    (-3, 1, 3, 1, 3, 1, -3, -1, -3, -1),
    (-3, -1, -3, 2, 3, 1),
    (3, 2, 3, 2, -3, -2, -3, -2),
)
P4_RELATIONS: Final[Relations] = (
    (-3, -2, 3, 1, -3, 2, 3, -1),
    (-3, -2, -3, -2, 3, 2, 3, 1, -3, -2, -3, 2, 3, 2, 3, 1, -3, -2, -3, 2, 3, 2, 3, 1, -3, -2, -3, -2, 3, 2, 3, -1, -3, -2, -3, -2, 3, 2, 3, -1),
    (-3, -2, -3, -2, 3, 2, 3, -1, -3, -2, -3, -2, 3, 2, 3, 2, 3, 1),
    (3, 2, 3, 2, 3, 2, 3, 2, -3, -2, -3, -2, -3, -2, -3, -2),
)
PNT_RELATIONS: Final[Relations] = (
    (-3, -2, 3, 1, -3, 2, 3, 1, -3, 2, 3, 1, -3, -2, 3, -1, -3, -2, 3, -1),
    (-3, -2, 3, -1, -3, -2, 3, 2, 3, 1),
    (2, 3, 1, -3, -2, -1),
    (2, 3, 2, 3, 2, -1, -3, -2, -3, -2, -3, -2, 3, 1),
)
P32_PLUS_RELATIONS: Final[Relations] = (
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -1, 2, 3, 2, -3, -2, 1, 2, 3, 2, -3, -2, 1, 2, 3, 2, -3, -2, 1, 2, 3, -2, -3, -2, -1, 2, 3, -2, -3, -2, -1, 2, 3, -2),
    (-3, -2, -1, 2, 3, -2, -3, -2, -1, 2, 3, -2, 1, 2, 3, 2),
    (-3, -2, 1, 2, 3, -2, 1, 2, 3, -2, -1, 2, -3, -2, -1, 2),
)
P32_MINUS_RELATIONS: Final[Relations] = (
    (2, 1, 2, 1, -2, -1, -2, -1),
    (-2, -1, -3, 1, 3, 1, 2, -1, -3, 1, 3, 1, 2, -1, -3, 1, 3, 1, -2, -1, -3, -1, 3, 1, -2, -1, -3, -1, 3, 1),
    (-2, -1, -3, -1, -3, 1, 3, 1, 2, -1, -3, 1, 3, 1),
    (-3, 1, 3, 1, 3, 1, 3, -1, -3, -1, -3, -1),
)


@dataclass(frozen=True, slots=True)
class PresentationCertificate:
    """A raw geometric-meridian presentation and its exact ``A6`` census."""

    name: str
    relations: Relations
    census: ThreeCyclePresentationCensus

    @property
    def verified(self) -> bool:
        """Whether no single-three-cycle assignment has image ``A6``."""

        return (
            self.census.assignments == 40**3
            and self.census.satisfying_assignments == 40
            and self.census.generated_order_histogram == ((3, 40),)
            and self.census.a6_assignments == 0
        )


def _presentation(name: str, relations: Relations) -> PresentationCertificate:
    """Replay one stored Sage presentation without a Sage dependency."""

    return PresentationCertificate(
        name=name,
        relations=relations,
        census=three_cycle_presentation_census(relations),
    )


@dataclass(frozen=True, slots=True)
class DeltaSixSevenCertificate:
    """Large-link corollary after eliminating collision delta five."""

    delta_six_candidates: tuple[OnePairDegreeCandidate, ...]
    delta_seven_candidates: tuple[OnePairDegreeCandidate, ...]
    t2_17: TorusA6Census
    t2_19: TorusA6Census
    t3_10: TorusA6Census
    t4_7: TorusA6Census

    @property
    def verified(self) -> bool:
        """Whether delta seven is the next possible large-link equality."""

        return (
            tuple(c.affine_degrees for c in self.delta_six_candidates)
            == ((2, 17),)
            and tuple(c.affine_degrees for c in self.delta_seven_candidates)
            == ((2, 19), (3, 10), (4, 7))
            and self.t2_17.single_three_cycle_a6_pairs == 0
            and self.t2_19.single_three_cycle_a6_pairs == 0
            and self.t3_10.single_three_cycle_a6_pairs == 720
            and self.t4_7.single_three_cycle_a6_pairs == 0
        )


@cache
def exact_delta_six_seven_certificate() -> DeltaSixSevenCertificate:
    """Build the exact genus and link census at deltas six and seven."""

    return DeltaSixSevenCertificate(
        delta_six_candidates=one_pair_degree_candidates(6),
        delta_seven_candidates=one_pair_degree_candidates(7),
        t2_17=torus_a6_census(2, 17),
        t2_19=torus_a6_census(2, 19),
        t3_10=torus_a6_census(3, 10),
        t4_7=torus_a6_census(4, 7),
    )


@dataclass(frozen=True, slots=True)
class A6DeltaFiveResidualCertificate:
    """Exact residual stratification and presentation exhaustion."""

    stratum_identities: tuple[Expr, ...]
    residual_ideal_identities: tuple[Expr, ...]
    exhaustion_identities: tuple[Expr, ...]
    exceptional_identities: tuple[Expr, ...]
    generic_valid_products: tuple[Expr, ...]
    exceptional_valid_products: tuple[Expr, ...]
    presentations: tuple[PresentationCertificate, ...]
    next_frontier: DeltaSixSevenCertificate

    @property
    def verified(self) -> bool:
        """Whether every residual and next-link checkpoint agrees."""

        return (
            all(identity == 0 for identity in self.stratum_identities)
            and all(
                identity == 0 for identity in self.residual_ideal_identities
            )
            and all(identity == 0 for identity in self.exhaustion_identities)
            and all(identity == 0 for identity in self.exceptional_identities)
            and self.generic_valid_products
            == (
                Rational(768, 5),
                Rational(1, 1296),
                Rational(250880, 9),
                Rational(-27900),
            )
            and all(value != 0 for value in self.exceptional_valid_products)
            and tuple(item.name for item in self.presentations)
            == (
                "contact three plus two nodes",
                "two contact-two points plus one node",
                "nonordinary triple plus one node",
                "ordinary triple plus separate contact two",
                "contact four plus one node",
                "higher nonordinary triple",
                "contact three plus contact two (positive embedding)",
                "contact three plus contact two (negative embedding)",
            )
            and all(item.verified for item in self.presentations)
            and self.next_frontier.verified
        )


@cache
def exact_a6_delta_five_residual_certificate(
) -> A6DeltaFiveResidualCertificate:
    """Construct the full conditional delta-five residual certificate."""

    a_substitution = _parameters(A_ALPHA, A_BETA, A_GAMMA)
    b_substitution = _parameters(B_ALPHA, B_BETA, B_GAMMA)
    n_substitution = _parameters(N_ALPHA, N_BETA, N_GAMMA)
    e_substitution = _parameters(E_ALPHA, E_BETA, E_GAMMA)

    p4_h = COLLISION_POLYNOMIAL.subs(_parameters(*P4_PARAMETERS))
    pnt_h = COLLISION_POLYNOMIAL.subs(_parameters(*PNT_PARAMETERS))
    p32_plus_h = COLLISION_POLYNOMIAL.subs(
        _parameters(*P32_PLUS_PARAMETERS)
    )
    p32_minus_h = COLLISION_POLYNOMIAL.subs(
        _parameters(*P32_MINUS_PARAMETERS)
    )
    p32_plus_other_root = -3 * (3 * P32_PLUS_H + 2) / 2
    p32_minus_other_root = -3 * (3 * P32_MINUS_H + 2) / 2

    residual_ideal_identities = tuple(
        _rational_pullback(
            polynomial,
            parameter,
            alpha,
            beta,
            gamma,
        )
        for parameter, alpha, beta, gamma, generators in (
            (
                H_PARAMETER,
                A_ALPHA,
                A_BETA,
                A_GAMMA,
                SINGULAR_DISCRIMINANT_GENERATORS,
            ),
            (
                M_PARAMETER,
                B_ALPHA,
                B_BETA,
                B_GAMMA,
                SINGULAR_DISCRIMINANT_GENERATORS,
            ),
            (
                N_PARAMETER,
                N_ALPHA,
                N_BETA,
                N_GAMMA,
                (
                    DEGREE_SEVEN_DISCRIMINANT_FACTOR,
                    TRIPLE_COLLISION_FACTOR,
                ),
            ),
            (
                E_PARAMETER,
                E_ALPHA,
                E_BETA,
                E_GAMMA,
                (
                    DEGREE_SEVEN_DISCRIMINANT_FACTOR,
                    TRIPLE_COLLISION_FACTOR,
                ),
            ),
        )
        for polynomial in generators
    )

    parameter_identifications = tuple(
        cancel(value.subs(parameter, point) - expected)
        for parameter, point, values, expected_values in (
            (
                H_PARAMETER,
                Rational(-13, 3),
                (A_ALPHA, A_BETA, A_GAMMA),
                P4_PARAMETERS,
            ),
            (
                M_PARAMETER,
                Rational(8, 3),
                (B_ALPHA, B_BETA, B_GAMMA),
                P4_PARAMETERS,
            ),
            (
                H_PARAMETER,
                P32_PLUS_H,
                (A_ALPHA, A_BETA, A_GAMMA),
                P32_PLUS_PARAMETERS,
            ),
            (
                H_PARAMETER,
                P32_MINUS_H,
                (A_ALPHA, A_BETA, A_GAMMA),
                P32_MINUS_PARAMETERS,
            ),
            (
                M_PARAMETER,
                3 - SQRT_SIX / 2,
                (B_ALPHA, B_BETA, B_GAMMA),
                P32_PLUS_PARAMETERS,
            ),
            (
                M_PARAMETER,
                3 + SQRT_SIX / 2,
                (B_ALPHA, B_BETA, B_GAMMA),
                P32_MINUS_PARAMETERS,
            ),
            (
                H_PARAMETER,
                Rational(-3),
                (A_ALPHA, A_BETA, A_GAMMA),
                PNT_PARAMETERS,
            ),
            (
                N_PARAMETER,
                Rational(1),
                (N_ALPHA, N_BETA, N_GAMMA),
                PNT_PARAMETERS,
            ),
            (
                E_PARAMETER,
                Rational(65),
                (E_ALPHA, E_BETA, E_GAMMA),
                PNT_PARAMETERS,
            ),
        )
        for value, expected in zip(values, expected_values, strict=True)
    )

    return A6DeltaFiveResidualCertificate(
        stratum_identities=(
            cancel(COLLISION_POLYNOMIAL.subs(a_substitution) - A_FACTORIZATION),
            cancel(EXTRA_CRITICAL_FACTOR.subs(a_substitution) - A_L_FORMULA),
            cancel(CUSP_COLLISION_FACTOR.subs(a_substitution) - A_C_FORMULA),
            cancel(TRIPLE_COLLISION_FACTOR.subs(a_substitution) - A_T_FORMULA),
            cancel(COLLISION_POLYNOMIAL.subs(b_substitution) - B_FACTORIZATION),
            cancel(EXTRA_CRITICAL_FACTOR.subs(b_substitution) - B_L_FORMULA),
            cancel(CUSP_COLLISION_FACTOR.subs(b_substitution) - B_C_FORMULA),
            cancel(TRIPLE_COLLISION_FACTOR.subs(b_substitution) - B_T_FORMULA),
            cancel(COLLISION_POLYNOMIAL.subs(n_substitution) - N_FACTORIZATION),
            cancel(EXTRA_CRITICAL_FACTOR.subs(n_substitution) - N_L_FORMULA),
            cancel(CUSP_COLLISION_FACTOR.subs(n_substitution) - N_C_FORMULA),
            cancel(TRIPLE_COLLISION_FACTOR.subs(n_substitution)),
            cancel(COLLISION_POLYNOMIAL.subs(e_substitution) - E_FACTORIZATION),
            cancel(EXTRA_CRITICAL_FACTOR.subs(e_substitution) - E_L_FORMULA),
            cancel(CUSP_COLLISION_FACTOR.subs(e_substitution) - E_C_FORMULA),
            cancel(TRIPLE_COLLISION_FACTOR.subs(e_substitution)),
        ),
        residual_ideal_identities=residual_ideal_identities,
        exhaustion_identities=(
            cancel(
                A_QUADRATIC_NUMERATOR.subs(S, H_PARAMETER)
                / (3 * H_PARAMETER + 11)
                - A_QUADRATIC_AT_TRIPLE_ROOT
            ),
            expand(discriminant(A_QUADRATIC_NUMERATOR, S) - A_QUADRATIC_DISCRIMINANT),
            cancel(discriminant(B_QUADRATIC, S) - B_QUADRATIC_DISCRIMINANT),
            cancel(B_QUADRATIC.subs(S, B_R) - B_SIMPLE_ROOT_VALUE),
            expand(N_QUADRATIC.subs(S, N_DOUBLE_ROOT) - N_QUADRATIC_AT_DOUBLE_ROOT),
            cancel(N_QUADRATIC.subs(S, N_NODE_ROOT) - N_QUADRATIC_AT_NODE_ROOT),
            expand(discriminant(N_QUADRATIC, S) - N_QUADRATIC_DISCRIMINANT),
            expand(
                discriminant(5 * E_TRIPLE_FACTOR, S)
                - E_TRIPLE_DISCRIMINANT
            ),
            expand(E_TRIPLE_FACTOR.subs(S, -3) - E_TRIPLE_AT_CONTACT_ROOT),
            cancel(
                N_DOUBLE_ROOT
                - N_NODE_ROOT
                - N_ROOT_COINCIDENCE_FORMULA
            ),
            expand(
                factor(
                    cancel(
                        A_CHANGE_FACTOR
                        / gcd(A_CHANGE_FACTOR, A_INVALID_FACTOR)
                    )
                )
                - A_EXCEPTION_FACTOR
            ),
            expand(
                factor(
                    cancel(
                        B_CHANGE_FACTOR
                        / gcd(B_CHANGE_FACTOR, B_INVALID_FACTOR)
                    )
                )
                - B_EXCEPTION_FACTOR
            ),
            expand(
                factor(
                    cancel(
                        N_CHANGE_FACTOR
                        / gcd(N_CHANGE_FACTOR, N_INVALID_FACTOR)
                    )
                )
                - N_EXCEPTION_FACTOR
            ),
            expand(
                factor(
                    cancel(
                        E_CHANGE_FACTOR
                        / gcd(E_CHANGE_FACTOR, E_INVALID_FACTOR)
                    )
                )
                - E_EXCEPTION_FACTOR
            ),
        ),
        exceptional_identities=(
            expand(
                p4_h
                - (S - 4) * (3 * S + 13) ** 4 / 81
            ),
            expand(pnt_h - (S + 3) ** 3 * (S**2 - S + 4)),
            expand(3 * P32_PLUS_H**2 + 12 * P32_PLUS_H + 4),
            expand(3 * P32_MINUS_H**2 + 12 * P32_MINUS_H + 4),
            expand(
                p32_plus_h
                - (S - P32_PLUS_H) ** 3
                * (S - p32_plus_other_root) ** 2
            ),
            expand(
                p32_minus_h
                - (S - P32_MINUS_H) ** 3
                * (S - p32_minus_other_root) ** 2
            ),
        )
        + parameter_identifications
        + (
            expand(N_C_FORMULA.subs(N_PARAMETER, -1)),
            expand(E_C_FORMULA.subs(E_PARAMETER, 5)),
            expand(
                E_L_FORMULA.subs(E_PARAMETER, Rational(155, 27))
            ),
        ),
        generic_valid_products=tuple(
            _valid_product(*sample)
            for sample in (A_SAMPLE, B_SAMPLE, N_SAMPLE, E_SAMPLE)
        ),
        exceptional_valid_products=(
            _valid_product(*P4_PARAMETERS),
            _valid_product(*PNT_PARAMETERS),
            _valid_product(*P32_PLUS_PARAMETERS),
            _valid_product(*P32_MINUS_PARAMETERS),
        ),
        presentations=(
            _presentation("contact three plus two nodes", A_RELATIONS),
            _presentation("two contact-two points plus one node", B_RELATIONS),
            _presentation("nonordinary triple plus one node", N_RELATIONS),
            _presentation("ordinary triple plus separate contact two", E_RELATIONS),
            _presentation("contact four plus one node", P4_RELATIONS),
            _presentation("higher nonordinary triple", PNT_RELATIONS),
            _presentation(
                "contact three plus contact two (positive embedding)",
                P32_PLUS_RELATIONS,
            ),
            _presentation(
                "contact three plus contact two (negative embedding)",
                P32_MINUS_RELATIONS,
            ),
        ),
        next_frontier=exact_delta_six_seven_certificate(),
    )


def main() -> int:
    """Print the residual exhaustion and fail on any exact regression."""

    certificate = exact_a6_delta_five_residual_certificate()
    print(
        "delta-five residual images:",
        {
            item.name: dict(item.census.generated_order_histogram)
            for item in certificate.presentations
        },
    )
    print(
        "next one-pair A6 frontier:",
        {
            "Delta=6": [
                candidate.affine_degrees
                for candidate in certificate.next_frontier.delta_six_candidates
            ],
            "Delta=7": [
                candidate.affine_degrees
                for candidate in certificate.next_frontier.delta_seven_candidates
            ],
            "surviving equality pair": (3, 10),
        },
    )
    print(f"conditional residual certificate verified: {certificate.verified}")
    print("conclusion under stated hypotheses: Delta >= 7")
    print("equality after the large-link test: affine (3,10), infinity (7,10)")
    print("claim boundary: conditional/computer-assisted; JC(2) remains open")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
