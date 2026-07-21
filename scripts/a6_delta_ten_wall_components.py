"""Certify the dominant delta-ten contact and ordinary-triple incidences.

The normalized degree-``(4, 9)`` family has coefficients

``P = t^2 + k*t^3 + t^4`` and
``Q = a*t^5 + b*t^6 + c*t^7 + d*t^8 + t^9``.

For a two-preimage collision, the collision decic ``H(s)`` and tangent
polynomial ``T(s)`` are affine-linear in ``(a,b,c,d)``.  This module checks
that their coefficient matrix has rank two away from exactly the forced
diagonal, pair-denominator, and split-incidence factors.  It also records the
exact restrictions on every rank-drop component and verifies that the known
contact-two member lies on the dominant component with one reduced incidence
root.

For a three-preimage collision, elementary symmetric coordinates ``(p,q,r)``
on the three source parameters satisfy the irreducible surface

``q^2 - p*r - q = 0``.

An exact Rabinowitsch-localized Groebner calculation proves that the two
``Q``-equality rows have rank two everywhere on the valid open of this
surface.  The calculation is dependency-free beyond SymPy; Sage is not
trusted as a runtime dependency.

These are component and connectedness certificates.  They do not classify
deeper intersections of the contact and triple-image walls, and they do not
prove the plane Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

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
    groebner,
)

from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    COLLISION_POLYNOMIAL,
    DELTA,
    GAMMA,
    KAPPA,
    PAIR_DENOMINATOR,
    PAIR_DIAGONAL_FACTOR,
    S,
    TANGENCY_POLYNOMIAL,
)
from scripts.a6_delta_ten_contact_wall import CONTACT_PARAMETERS

CONTACT_PARAMETERS_LINEAR: Final = (ALPHA, BETA, GAMMA, DELTA)
CONTACT_LOCALIZATION_VARIABLE: Final = Symbol("z_contact")

CONTACT_DIAGONAL_LINEAR_FORM: Final = (
    80 * ALPHA
    + 48 * BETA * S
    + 28 * GAMMA * S**2
    + 16 * DELTA * S**3
    + 9 * S**4
)
CONTACT_CUSP_VALUE_FORM: Final = (
    ALPHA + BETA * S + GAMMA * S**2 + DELTA * S**3 + S**4
)
CONTACT_CUSP_TANGENT_FORM: Final = (
    5 * ALPHA
    + 6 * BETA * S
    + 7 * GAMMA * S**2
    + 8 * DELTA * S**3
    + 9 * S**4
)

CONTACT_EXPECTED_RESIDUAL_BASIS: Final = (
    PAIR_DENOMINATOR**4,
    PAIR_DENOMINATOR**3 * (S**2 - 1),
    PAIR_DENOMINATOR**2 * (S**2 - 1) ** 2,
    (KAPPA + 2 * S**3) * (S**2 - 1) ** 3,
    S**2 * (S**2 - 1) ** 4,
)

TRIPLE_P: Final = Symbol("p")
TRIPLE_Q: Final = Symbol("q")
TRIPLE_R: Final = Symbol("r")
TRIPLE_T: Final = Symbol("triple_t")
TRIPLE_U: Final = Symbol("triple_u")
TRIPLE_V: Final = Symbol("triple_v")
TRIPLE_LOCALIZATION_VARIABLE: Final = Symbol("z_triple")

TRIPLE_BASE_POLYNOMIAL: Final = (
    TRIPLE_Q**2 - TRIPLE_P * TRIPLE_R - TRIPLE_Q
)
TRIPLE_BASE_DISCRIMINANT: Final = 4 * TRIPLE_P * TRIPLE_R + 1

TRIPLE_ELEMENTARY_P: Final = TRIPLE_T + TRIPLE_U + TRIPLE_V
TRIPLE_ELEMENTARY_Q: Final = (
    TRIPLE_T * TRIPLE_U
    + TRIPLE_T * TRIPLE_V
    + TRIPLE_U * TRIPLE_V
)
TRIPLE_ELEMENTARY_R: Final = TRIPLE_T * TRIPLE_U * TRIPLE_V
TRIPLE_BASE_ON_ROOTS: Final = expand(
    TRIPLE_BASE_POLYNOMIAL.subs(
        {
            TRIPLE_P: TRIPLE_ELEMENTARY_P,
            TRIPLE_Q: TRIPLE_ELEMENTARY_Q,
            TRIPLE_R: TRIPLE_ELEMENTARY_R,
        }
    )
)
TRIPLE_VALID_LOCALIZER: Final = expand(
    TRIPLE_T
    * TRIPLE_U
    * TRIPLE_V
    * (TRIPLE_T - TRIPLE_U)
    * (TRIPLE_T - TRIPLE_V)
    * (TRIPLE_U - TRIPLE_V)
    * TRIPLE_ELEMENTARY_Q
)
TRIPLE_SAMPLE: Final = (
    Rational(-3, 5),
    Rational(-2, 5),
    Rational(1, 5),
)


def _two_by_two_minors(
    first_row: tuple[Expr, ...],
    second_row: tuple[Expr, ...],
) -> tuple[Expr, ...]:
    """Return all ordered-column ``2 x 2`` minors of two equal-length rows."""

    if len(first_row) != len(second_row):
        msg = "minor rows must have equal length"
        raise ValueError(msg)
    return tuple(
        expand(
            first_row[left] * second_row[right]
            - first_row[right] * second_row[left]
        )
        for left in range(len(first_row))
        for right in range(left + 1, len(first_row))
    )


def _polynomial_gcd(polynomials: tuple[Expr, ...]) -> Expr:
    """Return the normalized iterated polynomial gcd of a nonempty tuple."""

    if not polynomials:
        msg = "a polynomial gcd needs at least one input"
        raise ValueError(msg)
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = gcd(common, polynomial)
    return common


def _groebner_expressions(
    polynomials: tuple[Expr, ...],
    *generators: Symbol,
    order: str,
) -> tuple[Expr, ...]:
    """Return a Groebner basis as ordinary expressions."""

    basis = groebner(polynomials, *generators, order=order)
    return tuple(polynomial.as_expr() for polynomial in basis.polys)


def _contact_coefficient_rows() -> tuple[tuple[Expr, ...], tuple[Expr, ...]]:
    """Return coefficient rows of ``H,T`` in ``(a,b,c,d)``."""

    return (
        tuple(
            diff(COLLISION_POLYNOMIAL, parameter)
            for parameter in CONTACT_PARAMETERS_LINEAR
        ),
        tuple(
            diff(TANGENCY_POLYNOMIAL, parameter)
            for parameter in CONTACT_PARAMETERS_LINEAR
        ),
    )


def _triple_equality_rows(maximum_power: int) -> tuple[tuple[Expr, ...], tuple[Expr, ...]]:
    """Return raw ``Q(t)=Q(u)=Q(v)`` rows through ``maximum_power``."""

    if maximum_power < 6:
        msg = "the triple equality matrix needs at least powers five and six"
        raise ValueError(msg)
    powers = range(5, maximum_power + 1)
    return (
        tuple(TRIPLE_T**power - TRIPLE_U**power for power in powers),
        tuple(TRIPLE_T**power - TRIPLE_V**power for power in powers),
    )


@dataclass(frozen=True, slots=True)
class DeltaTenWallComponentsCertificate:
    """Exact rank, component, and sample data for the first delta-ten walls."""

    contact_minor_gcd_identity: Expr
    contact_residual_groebner_basis: tuple[Expr, ...]
    contact_localized_groebner_basis: tuple[Expr, ...]
    contact_sample_ab_minor: Expr
    contact_sample_gcd: Expr
    contact_sample_second_derivative: Expr
    contact_sample_tangency_derivative: Expr
    contact_s_zero_identities: tuple[Expr, Expr]
    contact_denominator_identities: tuple[Expr, Expr]
    contact_diagonal_identities: tuple[Expr, Expr]
    contact_cusp_boundary_identities: tuple[Expr, Expr]
    triple_base_factorization: Expr
    triple_base_discriminant: Expr
    triple_base_discriminant_squarefree_part: Expr
    triple_base_discriminant_p_degree: int
    triple_coefficient_localized_groebner_basis: tuple[Expr, ...]
    triple_augmented_localized_groebner_basis: tuple[Expr, ...]
    triple_sample_base_value: Expr
    triple_sample_localizer: Expr
    triple_sample_kappa: Expr
    triple_sample_ab_minor: Expr

    @property
    def verified(self) -> bool:
        """Whether every exact component and rank check agrees."""

        expected_contact_basis = tuple(
            expand(polynomial) for polynomial in CONTACT_EXPECTED_RESIDUAL_BASIS
        )
        return bool(
            self.contact_minor_gcd_identity == 0
            and self.contact_residual_groebner_basis == expected_contact_basis
            and self.contact_localized_groebner_basis == (1,)
            and self.contact_sample_ab_minor == 131220
            and self.contact_sample_gcd == S + 2
            and self.contact_sample_second_derivative == Rational(2538, 5)
            and self.contact_sample_tangency_derivative == Rational(30456, 5)
            and self.contact_s_zero_identities == (0, 0)
            and self.contact_denominator_identities == (0, 0)
            and self.contact_diagonal_identities == (0, 0)
            and self.contact_cusp_boundary_identities == (0, 0)
            and self.triple_base_factorization == TRIPLE_BASE_POLYNOMIAL
            and self.triple_base_discriminant == TRIPLE_BASE_DISCRIMINANT
            and self.triple_base_discriminant_squarefree_part
            == TRIPLE_BASE_DISCRIMINANT
            and self.triple_base_discriminant_p_degree == 1
            and self.triple_coefficient_localized_groebner_basis == (1,)
            and self.triple_augmented_localized_groebner_basis == (1,)
            and self.triple_sample_base_value == 0
            and self.triple_sample_localizer != 0
            and self.triple_sample_kappa == 2
            and self.triple_sample_ab_minor == Rational(8652, 48828125)
        )


@cache
def exact_delta_ten_wall_components_certificate() -> DeltaTenWallComponentsCertificate:
    """Build the exact dominant-wall and exceptional-component certificate."""

    contact_rows = _contact_coefficient_rows()
    contact_minors = _two_by_two_minors(*contact_rows)
    contact_minor_gcd = _polynomial_gcd(contact_minors)
    expected_contact_gcd = (
        S * PAIR_DENOMINATOR**3 * PAIR_DIAGONAL_FACTOR
    )
    normalized_contact_minors = tuple(
        cancel(minor / contact_minor_gcd) for minor in contact_minors
    )
    contact_residual_basis = _groebner_expressions(
        normalized_contact_minors,
        KAPPA,
        S,
        order="lex",
    )
    contact_rank_localizer = (
        S * PAIR_DENOMINATOR * PAIR_DIAGONAL_FACTOR
    )
    contact_localized_basis = _groebner_expressions(
        (
            *normalized_contact_minors,
            1 - CONTACT_LOCALIZATION_VARIABLE * contact_rank_localizer,
        ),
        CONTACT_LOCALIZATION_VARIABLE,
        KAPPA,
        S,
        order="grevlex",
    )

    sample_collision = expand(COLLISION_POLYNOMIAL.subs(CONTACT_PARAMETERS))
    sample_tangency = expand(TANGENCY_POLYNOMIAL.subs(CONTACT_PARAMETERS))
    sample_rows = tuple(
        tuple(expand(entry.subs(CONTACT_PARAMETERS).subs(S, -2)) for entry in row)
        for row in contact_rows
    )
    sample_minors = _two_by_two_minors(sample_rows[0], sample_rows[1])
    sample_gcd = gcd(
        Poly(sample_collision, S),
        Poly(sample_tangency, S),
    ).as_expr()

    diagonal_kappa = -(2 * S**2 + 4) / (3 * S)
    cusp_kappa = -(S**2 + 1) / S
    contact_diagonal_identities = (
        cancel(
            81
            * S**2
            * COLLISION_POLYNOMIAL.subs(KAPPA, diagonal_kappa)
            - (S**2 - 1) ** 4 * CONTACT_DIAGONAL_LINEAR_FORM
        ),
        cancel(
            81
            * S**3
            * TANGENCY_POLYNOMIAL.subs(KAPPA, diagonal_kappa)
            - 4
            * (S**2 - 1) ** 4
            * (S**2 - 2)
            * CONTACT_DIAGONAL_LINEAR_FORM
        ),
    )
    contact_cusp_boundary_identities = (
        cancel(
            S**2 * COLLISION_POLYNOMIAL.subs(KAPPA, cusp_kappa)
            - (S**2 - 1) ** 4 * CONTACT_CUSP_VALUE_FORM
        ),
        cancel(
            S**3 * TANGENCY_POLYNOMIAL.subs(KAPPA, cusp_kappa)
            + 2 * (S**2 - 1) ** 4 * CONTACT_CUSP_TANGENT_FORM
        ),
    )

    triple_discriminant = discriminant(
        TRIPLE_BASE_POLYNOMIAL,
        TRIPLE_Q,
    )
    triple_discriminant_polynomial = Poly(
        triple_discriminant,
        TRIPLE_P,
        TRIPLE_R,
    )
    triple_coefficient_rows = _triple_equality_rows(8)
    triple_augmented_rows = _triple_equality_rows(9)
    triple_coefficient_minors = _two_by_two_minors(*triple_coefficient_rows)
    triple_augmented_minors = _two_by_two_minors(*triple_augmented_rows)
    triple_coefficient_localized_basis = _groebner_expressions(
        (
            TRIPLE_BASE_ON_ROOTS,
            *triple_coefficient_minors,
            1 - TRIPLE_LOCALIZATION_VARIABLE * TRIPLE_VALID_LOCALIZER,
        ),
        TRIPLE_LOCALIZATION_VARIABLE,
        TRIPLE_T,
        TRIPLE_U,
        TRIPLE_V,
        order="grevlex",
    )
    triple_augmented_localized_basis = _groebner_expressions(
        (
            TRIPLE_BASE_ON_ROOTS,
            *triple_augmented_minors,
            1 - TRIPLE_LOCALIZATION_VARIABLE * TRIPLE_VALID_LOCALIZER,
        ),
        TRIPLE_LOCALIZATION_VARIABLE,
        TRIPLE_T,
        TRIPLE_U,
        TRIPLE_V,
        order="grevlex",
    )
    triple_sample_substitution = dict(
        zip(
            (TRIPLE_T, TRIPLE_U, TRIPLE_V),
            TRIPLE_SAMPLE,
            strict=True,
        )
    )
    triple_sample_p = TRIPLE_ELEMENTARY_P.subs(triple_sample_substitution)
    triple_sample_q = TRIPLE_ELEMENTARY_Q.subs(triple_sample_substitution)
    triple_sample_r = TRIPLE_ELEMENTARY_R.subs(triple_sample_substitution)

    return DeltaTenWallComponentsCertificate(
        contact_minor_gcd_identity=expand(
            contact_minor_gcd - expected_contact_gcd
        ),
        contact_residual_groebner_basis=contact_residual_basis,
        contact_localized_groebner_basis=contact_localized_basis,
        contact_sample_ab_minor=sample_minors[0],
        contact_sample_gcd=sample_gcd,
        contact_sample_second_derivative=diff(
            sample_collision,
            S,
            2,
        ).subs(S, -2),
        contact_sample_tangency_derivative=diff(
            sample_tangency,
            S,
        ).subs(S, -2),
        contact_s_zero_identities=(
            expand(COLLISION_POLYNOMIAL.subs(S, 0) - ALPHA * KAPPA**2),
            expand(
                TANGENCY_POLYNOMIAL.subs(S, 0)
                - 10 * ALPHA * KAPPA**3
            ),
        ),
        contact_denominator_identities=(
            expand(
                COLLISION_POLYNOMIAL.subs(KAPPA, -2 * S)
                - S**2 * (S**2 - 1) ** 4
            ),
            expand(
                TANGENCY_POLYNOMIAL.subs(KAPPA, -2 * S)
                + 36 * S**3 * (S**2 - 1) ** 4
            ),
        ),
        contact_diagonal_identities=contact_diagonal_identities,
        contact_cusp_boundary_identities=contact_cusp_boundary_identities,
        triple_base_factorization=factor(TRIPLE_BASE_POLYNOMIAL),
        triple_base_discriminant=triple_discriminant,
        triple_base_discriminant_squarefree_part=(
            triple_discriminant_polynomial.sqf_part().as_expr()
        ),
        triple_base_discriminant_p_degree=Poly(
            triple_discriminant,
            TRIPLE_P,
        ).degree(),
        triple_coefficient_localized_groebner_basis=(
            triple_coefficient_localized_basis
        ),
        triple_augmented_localized_groebner_basis=(
            triple_augmented_localized_basis
        ),
        triple_sample_base_value=TRIPLE_BASE_POLYNOMIAL.subs(
            {
                TRIPLE_P: triple_sample_p,
                TRIPLE_Q: triple_sample_q,
                TRIPLE_R: triple_sample_r,
            }
        ),
        triple_sample_localizer=TRIPLE_VALID_LOCALIZER.subs(
            triple_sample_substitution
        ),
        triple_sample_kappa=cancel(
            triple_sample_r / triple_sample_q - triple_sample_p
        ),
        triple_sample_ab_minor=triple_coefficient_minors[0].subs(
            triple_sample_substitution
        ),
    )


def main() -> int:
    """Print the exact rank certificate and fail on any regression."""

    certificate = exact_delta_ten_wall_components_certificate()
    print("delta-ten wall component certificate:", certificate.verified)
    print(
        "contact rank-drop factors:",
        S * PAIR_DENOMINATOR * PAIR_DIAGONAL_FACTOR,
    )
    print(
        "triple base and sample AB minor:",
        TRIPLE_BASE_POLYNOMIAL,
        certificate.triple_sample_ab_minor,
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
