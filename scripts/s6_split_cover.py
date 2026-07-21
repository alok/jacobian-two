"""Exact degree-six ``S6`` cover with a globally split boundary preimage.

The finite flat map

``(x,z) |-> (x^2+z, x^3-x^4*z)``

extends to a finite degree-six morphism ``P1 x P1 -> P2`` whose sole target
boundary line has two noncontracted source components.  Its geometric
monodromy is ``S6``.  This is a hostile model for attempts to derive
Orevkov's condition (*) from degree, finiteness, monodromy, and boundary
intersection transport alone.

The map is not Keller: its affine Jacobian is nonconstant and it has affine
branching.  It is not a counterexample to the Jacobian conjecture and does not
realize the surviving Keller passport.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from sympy import (
    Expr,
    Matrix,
    Poly,
    QQ,
    Symbol,
    diff,
    discriminant,
    expand,
    gcd,
)

X: Final = Symbol("x")
Z: Final = Symbol("z")
P: Final = Symbol("P")
Q: Final = Symbol("Q")

X0: Final = Symbol("X0")
X1: Final = Symbol("X1")
Y0: Final = Symbol("Y0")
Y1: Final = Symbol("Y1")

AFFINE_P: Final = X**2 + Z
AFFINE_Q: Final = X**3 - X**4 * Z
FIBER_POLYNOMIAL: Final = X**6 - P * X**4 + X**3 - Q
AFFINE_JACOBIAN: Final = -2 * X**5 + 4 * X**3 * Z - 3 * X**2

BRANCH_FACTOR: Final = (
    1024 * P**6 * Q
    + 13824 * P**3 * Q**2
    - 8640 * P**3 * Q
    - 108 * P**3
    + 46656 * Q**3
    + 34992 * Q**2
    + 8748 * Q
    + 729
)


def projective_sections() -> tuple[Expr, Expr, Expr, Expr]:
    """Return ``D`` and the three basepoint-free ``O(3,1)`` sections."""

    diagonal = X0 * Y0 + X1 * Y1
    first = X0**2 * diagonal
    second = X0**3 * Y1 + X1**2 * diagonal
    third = X1**3 * Y0
    return diagonal, first, second, third


def split_boundary_matrix() -> Matrix:
    """Return the intersection matrix of ``C1=(X0=0)`` and ``C2=(D=0)``."""

    return Matrix([[0, 1], [1, 2]])


def blown_up_corridor_matrix() -> Matrix:
    """Return the boundary matrix after one corner and one free blowup."""

    return Matrix(
        [
            [-1, 1, 0, 0],
            [1, -2, 1, 0],
            [0, 1, -1, 1],
            [0, 0, 1, 1],
        ]
    )


def _entries(matrix: Matrix) -> tuple[Expr, ...]:
    """Return matrix entries in stable row-major order."""

    return tuple(expand(entry) for entry in matrix)


def _deleted_negative_determinant(matrix: Matrix, vertex: int) -> Expr:
    """Delete a vertex and take the determinant of the negative matrix."""

    return expand((-matrix).minor_submatrix(vertex, vertex).det())


@dataclass(frozen=True, slots=True)
class S6SplitCoverCertificate:
    """Exact algebraic, projective, boundary, and monodromy premises."""

    affine_relation: Expr
    affine_inverse_coordinate: Expr
    jacobian: Expr
    fiber_degree: int
    fiber_leading_coefficient: Expr
    x0_boundary_identities: tuple[Expr, Expr, Expr]
    diagonal_boundary_identities: tuple[Expr, Expr, Expr, Expr]
    affine_chart_identities: tuple[Expr, Expr, Expr, Expr]
    split_transport_defect: tuple[Expr, ...]
    split_degree: Expr
    split_pullback_square: Expr
    corridor_determinant: Expr
    corridor_labels: tuple[int, int, int, int]
    corridor_transport_defect: tuple[Expr, ...]
    corridor_degree: Expr
    corridor_pullback_square: Expr
    corridor_determinant_labels: tuple[Expr, Expr]
    discriminant_factor_identity: Expr
    branch_factor_gcd: Expr
    quadratic_after_cubic_residual: Expr
    cubic_after_quadratic_odd_term: Expr

    @property
    def finite_flat_premises_verified(self) -> bool:
        """Check the monic rank-six presentation and affine coordinate inverse."""

        return bool(
            self.affine_relation == 0
            and self.affine_inverse_coordinate == 0
            and self.fiber_degree == 6
            and self.fiber_leading_coefficient == 1
        )

    @property
    def projective_extension_verified(self) -> bool:
        """Check the exact chart identities used in basepoint-freeness."""

        return all(
            identity == 0
            for identity in (
                *self.x0_boundary_identities,
                *self.diagonal_boundary_identities,
                *self.affine_chart_identities,
            )
        )

    @property
    def split_boundary_verified(self) -> bool:
        """Check both the original split boundary and blown-up corridor."""

        return bool(
            all(entry == 0 for entry in self.split_transport_defect)
            and self.split_degree == 6
            and self.split_pullback_square == 6
            and self.corridor_determinant == -1
            and self.corridor_labels == (2, 1, 0, -1)
            and all(entry == 0 for entry in self.corridor_transport_defect)
            and self.corridor_degree == 6
            and self.corridor_pullback_square == 6
            and self.corridor_determinant_labels == (-3, -2)
        )

    @property
    def s6_monodromy_premises_verified(self) -> bool:
        """Check indecomposability and a reduced simple-branch factor.

        A degree-six polynomial has only decomposition patterns ``2 o 3`` and
        ``3 o 2``.  The two stored nonzero residuals are the final coefficient
        contradictions in those normalized patterns.  Indecomposability makes
        polynomial monodromy primitive; a primitive subgroup containing a
        transposition is ``S6``.
        """

        return bool(
            self.discriminant_factor_identity == 0
            and self.branch_factor_gcd == 1
            and self.quadratic_after_cubic_residual != 0
            and self.cubic_after_quadratic_odd_term == 1
        )

    @property
    def verified(self) -> bool:
        """Whether every exact premise and hostile boundary identity passes."""

        return bool(
            self.finite_flat_premises_verified
            and self.projective_extension_verified
            and self.split_boundary_verified
            and self.s6_monodromy_premises_verified
            and self.jacobian == AFFINE_JACOBIAN
            and self.jacobian != 1
        )


def exact_s6_split_cover_certificate() -> S6SplitCoverCertificate:
    """Build the complete exact stopping certificate."""

    relation = expand(
        FIBER_POLYNOMIAL.subs(
            {
                P: AFFINE_P,
                Q: AFFINE_Q,
            }
        )
    )
    inverse_coordinate = expand(Z - (AFFINE_P - X**2))
    jacobian = expand(
        diff(AFFINE_P, X) * diff(AFFINE_Q, Z)
        - diff(AFFINE_P, Z) * diff(AFFINE_Q, X)
    )
    fiber_poly = Poly(FIBER_POLYNOMIAL, X)

    diagonal, first, second, third = projective_sections()
    x0_substitution = {X0: 0}
    diagonal_parametrization = {Y0: -X1, Y1: X0}
    affine_substitution = {
        X0: 1,
        X1: X,
        Y0: 1 - X * Z,
        Y1: Z,
    }

    split_matrix = split_boundary_matrix()
    split_normal = Matrix([2, 1])
    split_tangential = Matrix([1, 4])

    corridor = blown_up_corridor_matrix()
    corridor_normal = Matrix([2, 2, 3, 1])
    corridor_tangential = Matrix([0, 1, 0, 4])

    computed_discriminant = discriminant(FIBER_POLYNOMIAL, X)
    branch_poly = Poly(BRANCH_FACTOR, Q, domain=QQ.frac_field(P))
    branch_derivative = Poly(diff(BRANCH_FACTOR, Q), Q, domain=QQ.frac_field(P))

    inner_cubic = X**3 - P * X / 2
    forced_quadratic_outer = inner_cubic**2 + inner_cubic
    polynomial_without_q = X**6 - P * X**4 + X**3

    return S6SplitCoverCertificate(
        affine_relation=relation,
        affine_inverse_coordinate=inverse_coordinate,
        jacobian=jacobian,
        fiber_degree=int(fiber_poly.degree()),
        fiber_leading_coefficient=expand(fiber_poly.LC()),
        x0_boundary_identities=(
            expand(first.subs(x0_substitution)),
            expand(second.subs(x0_substitution) - X1**3 * Y1),
            expand(third.subs(x0_substitution) - X1**3 * Y0),
        ),
        diagonal_boundary_identities=(
            expand(diagonal.subs(diagonal_parametrization)),
            expand(first.subs(diagonal_parametrization)),
            expand(second.subs(diagonal_parametrization) - X0**4),
            expand(third.subs(diagonal_parametrization) + X1**4),
        ),
        affine_chart_identities=(
            expand(diagonal.subs(affine_substitution) - 1),
            expand(first.subs(affine_substitution) - 1),
            expand(second.subs(affine_substitution) - AFFINE_P),
            expand(third.subs(affine_substitution) - AFFINE_Q),
        ),
        split_transport_defect=_entries(
            split_matrix * split_normal - split_tangential
        ),
        split_degree=expand((split_normal.T * split_tangential)[0]),
        split_pullback_square=expand(
            (split_normal.T * split_matrix * split_normal)[0]
        ),
        corridor_determinant=expand(corridor.det()),
        corridor_labels=(2, 1, 0, -1),
        corridor_transport_defect=_entries(
            corridor * corridor_normal - corridor_tangential
        ),
        corridor_degree=expand(
            (corridor_normal.T * corridor_tangential)[0]
        ),
        corridor_pullback_square=expand(
            (corridor_normal.T * corridor * corridor_normal)[0]
        ),
        corridor_determinant_labels=(
            _deleted_negative_determinant(corridor, 0),
            _deleted_negative_determinant(corridor, 1),
        ),
        discriminant_factor_identity=expand(
            computed_discriminant - Q**2 * BRANCH_FACTOR
        ),
        branch_factor_gcd=expand(gcd(branch_poly, branch_derivative).as_expr()),
        quadratic_after_cubic_residual=expand(
            forced_quadratic_outer - polynomial_without_q
        ),
        cubic_after_quadratic_odd_term=Poly(
            polynomial_without_q,
            X,
        ).coeff_monomial(X**3),
    )


def main() -> int:
    """Print the exact cover certificate and fail on any regression."""

    certificate = exact_s6_split_cover_certificate()
    print(
        "finite affine cover:",
        {
            "degree": certificate.fiber_degree,
            "Jacobian": certificate.jacobian,
            "finite-flat premises": certificate.finite_flat_premises_verified,
        },
    )
    print(
        "split boundary:",
        {
            "pullback": "2*C1+C2",
            "degrees": (1, 4),
            "corridor labels": certificate.corridor_labels,
            "verified": certificate.split_boundary_verified,
        },
    )
    print(
        "geometric monodromy premises:",
        {
            "indecomposable": certificate.quadratic_after_cubic_residual != 0,
            "simple branch factor": certificate.branch_factor_gcd == 1,
            "S6 premises": certificate.s6_monodromy_premises_verified,
        },
    )
    print(f"S6 split-cover certificate verified: {certificate.verified}")
    print("claim boundary: finite S6 cover, but not Keller")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
