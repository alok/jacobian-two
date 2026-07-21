"""Exact local multiplication fixtures for the hostile ``A6`` trace lattice.

The hostile trace matrix has the correct determinant, normalization cokernel,
and rank jumps, but those quadratic data do not include a multiplication law.
This module checks that its two forced singular fibers nevertheless support the
expected commutative associative algebras:

* ``C[z]/(z^5) x C`` at the degree-five cusp; and
* ``C[z]/(z^3) x C[w]/(w^3)`` at a ``3+3`` collision.

It also verifies the normalization presentation and a one-parameter hostile
family with the same ``3+1+1+1 -> 5+1`` and ``3+1+1+1 -> 3+3`` partition
jumps.  These are pointwise and normalization-level stopping fixtures.  They
do not construct a global rank-six algebra, an ``A6`` cover, or a Keller map.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from sympy import (
    Expr,
    I,
    Matrix,
    Poly,
    Symbol,
    diag,
    expand,
    factor,
    simplify,
    sqrt,
    zeros,
)

P: Final = Symbol("P")
Q: Final = Symbol("Q")
PARAMETER: Final = Symbol("t")
COLLISION_P: Final = Symbol("rho")
FIBER_VARIABLE: Final = Symbol("z")

BRANCH_EQUATION: Final = (
    -P**5 + 5 * P**3 * Q - 5 * P * Q**2 + Q**3 + Q**2
)
COLLISION_RELATION: Final = COLLISION_P**2 + COLLISION_P - 1
CYCLOTOMIC_FIVE: Final = (
    PARAMETER**4
    + PARAMETER**3
    + PARAMETER**2
    + PARAMETER
    + 1
)


def a6_presentation_matrix() -> Matrix:
    """Return the symmetric normalization presentation ``K``."""

    return Matrix(
        [
            [-P**2, P**2 - Q, P - Q],
            [P**2 - Q, -Q, -P],
            [P - Q, -P, P - 1],
        ]
    )


def a6_trace_data() -> tuple[Matrix, Matrix, Matrix]:
    """Return ``(Phi, H, T)`` with ``T = Phi.T * H * Phi``."""

    presentation = a6_presentation_matrix()
    identity = Matrix.eye(3)
    zero = zeros(3)
    phi = diag(presentation, identity)
    middle = zero.row_join(identity).col_join(identity.row_join(zero))
    return phi, middle, phi.T * middle * phi


def a6_trace_unit() -> Matrix:
    """Return the polynomial-basis coordinates of the trace unit."""

    return Matrix([0, 0, 1, 0, -3, -3])


def _entries(matrix: Matrix) -> tuple[Expr, ...]:
    """Return matrix entries in a stable flat tuple."""

    return tuple(entry for entry in matrix)


def _simplified_entries(matrix: Matrix) -> tuple[Expr, ...]:
    """Simplify all entries of a matrix."""

    return tuple(simplify(entry) for entry in matrix)


def _reduce_at_collision(expression: Expr) -> Expr:
    """Reduce an expression modulo ``rho^2 + rho - 1``."""

    polynomial = Poly(expand(expression), COLLISION_P)
    relation = Poly(COLLISION_RELATION, COLLISION_P)
    return factor(polynomial.rem(relation).as_expr())


def _collision_entries(matrix: Matrix) -> tuple[Expr, ...]:
    """Reduce every matrix entry in the collision residue field."""

    return tuple(_reduce_at_collision(entry) for entry in matrix)


def _cusp_basis() -> Matrix:
    """Return the basis realizing the ``5+1`` trace form at the cusp.

    Its columns correspond to ``(p, z, z^2, z^3, z^4, q)`` in
    ``C[z]/(z^5) x C``.  The last two nilpotents map to ``ker(Phi_0)``.
    """

    standard = Matrix.eye(6)
    unit = a6_trace_unit()
    delta = I * sqrt(5)
    coefficient_e3 = (5 + delta) / 6
    coefficient_f3 = (-5 + delta) / 2
    length_five_unit = (
        coefficient_e3 * standard[:, 2]
        + coefficient_f3 * standard[:, 5]
    )
    simple_unit = unit - length_five_unit
    return Matrix.hstack(
        length_five_unit,
        standard[:, 3],
        standard[:, 4],
        standard[:, 0],
        standard[:, 1],
        simple_unit,
    )


def _collision_basis() -> Matrix:
    """Return the basis realizing the ``3+3`` trace form at a node.

    Its columns correspond to ``(p_+, z, z^2, p_-, w, w^2)``.  The two
    top nilpotents map to ``ker(Phi_rho)`` modulo the collision relation.
    """

    standard = Matrix.eye(6)
    unit = a6_trace_unit()
    trace_orthogonal = standard[:, 2] - unit / 2
    first_unit = unit / 2 + I * trace_orthogonal
    second_unit = unit / 2 - I * trace_orthogonal
    first_linear = Matrix([0, 0, 0, 1, -COLLISION_P, 0])
    first_square = Matrix([1, -COLLISION_P, 0, 0, 0, 0])
    second_linear = Matrix([0, 0, 0, -1, 0, 1])
    second_square = Matrix([-1, 0, 1, 0, 0, 0])
    return Matrix.hstack(
        first_unit,
        first_linear,
        first_square,
        second_unit,
        second_linear,
        second_square,
    )


def _cusp_frobenius_data() -> tuple[Matrix, Matrix]:
    """Return the perfect middle form and multiplication by ``(z^2,1)``."""

    middle = zeros(6)
    for left in range(5):
        for right in range(5):
            if left + right == 4:
                middle[left, right] = 5
    middle[5, 5] = 1

    section = zeros(6)
    section[2, 0] = 1
    section[3, 1] = 1
    section[4, 2] = 1
    section[5, 5] = 1
    return middle, section


def _collision_frobenius_data() -> tuple[Matrix, Matrix]:
    """Return the perfect middle form and multiplication by ``(z,w)``."""

    middle = zeros(6)
    for offset in (0, 3):
        for left in range(3):
            for right in range(3):
                if left + right == 2:
                    middle[offset + left, offset + right] = 3

    section = zeros(6)
    section[1, 0] = 1
    section[2, 1] = 1
    section[4, 3] = 1
    section[5, 4] = 1
    return middle, section


def hostile_normalization_family() -> Expr:
    """Return a monic rank-six family over the normalization parameter.

    The family is deliberately only over ``C[t]``.  It is not claimed to
    descend through the singular branch curve or to realize the global trace
    matrix.
    """

    first_simple_root = 1 - CYCLOTOMIC_FIVE
    second_simple_root = 1 - (1 + PARAMETER) * CYCLOTOMIC_FIVE
    return expand(
        FIBER_VARIABLE**3
        * (FIBER_VARIABLE - first_simple_root)
        * (FIBER_VARIABLE - second_simple_root)
        * (FIBER_VARIABLE - 1)
    )


@dataclass(frozen=True)
class A6LocalMultiplicationCertificate:
    """Exact outputs for the local multiplication stopping fixtures."""

    presentation_determinant: Expr
    adjugate_identity: tuple[Expr, ...]
    parametrized_branch_equation: Expr
    normalization_kernel: tuple[Expr, ...]
    normalized_adjugate_factor_identity: tuple[Expr, ...]
    trace_factor_identity: tuple[Expr, ...]
    trace_unit_norm: Expr
    cusp_basis_determinant: Expr
    cusp_gram_identity: tuple[Expr, ...]
    cusp_unit_identity: tuple[Expr, ...]
    cusp_middle_determinant: Expr
    cusp_middle_trace_identity: tuple[Expr, ...]
    cusp_phi_kernel_columns: tuple[tuple[Expr, ...], tuple[Expr, ...]]
    collision_basis_determinant: Expr
    collision_gram_identity: tuple[Expr, ...]
    collision_unit_identity: tuple[Expr, ...]
    collision_middle_determinant: Expr
    collision_middle_trace_identity: tuple[Expr, ...]
    collision_phi_kernel_columns: tuple[tuple[Expr, ...], tuple[Expr, ...]]
    hostile_family_degree: int
    hostile_family_leading_coefficient: Expr
    hostile_family_cusp_identity: Expr
    hostile_family_collision_identity: Expr

    @property
    def verified(self) -> bool:
        """Whether every exact pointwise and normalization check passes."""

        zero_tuples = (
            self.adjugate_identity,
            self.normalization_kernel,
            self.normalized_adjugate_factor_identity,
            self.trace_factor_identity,
            self.cusp_gram_identity,
            self.cusp_unit_identity,
            self.cusp_middle_trace_identity,
            *self.cusp_phi_kernel_columns,
            self.collision_gram_identity,
            self.collision_unit_identity,
            self.collision_middle_trace_identity,
            *self.collision_phi_kernel_columns,
        )
        return bool(
            expand(self.presentation_determinant - BRANCH_EQUATION) == 0
            and self.parametrized_branch_equation == 0
            and all(
                entry == 0
                for entries in zero_tuples
                for entry in entries
            )
            and self.trace_unit_norm == 6
            and self.cusp_basis_determinant == -I * sqrt(5)
            and self.cusp_middle_determinant != 0
            and self.collision_basis_determinant == 3 * I
            and self.collision_middle_determinant != 0
            and self.hostile_family_degree == 6
            and self.hostile_family_leading_coefficient == 1
            and self.hostile_family_cusp_identity == 0
            and self.hostile_family_collision_identity == 0
        )


def exact_a6_local_multiplication_certificate(
) -> A6LocalMultiplicationCertificate:
    """Build all exact local multiplication and normalization checks."""

    presentation = a6_presentation_matrix()
    adjugate = presentation.adjugate()
    identity_three = Matrix.eye(3)
    parameter_p = PARAMETER**2 + PARAMETER**3
    parameter_q = PARAMETER**5
    parameter_substitution = {P: parameter_p, Q: parameter_q}
    parametrized_presentation = presentation.subs(parameter_substitution)
    normalization_generators = Matrix(
        [1, PARAMETER, PARAMETER**2]
    )
    expected_normalized_adjugate = (
        -PARAMETER**4
        * CYCLOTOMIC_FIVE
        * normalization_generators
        * normalization_generators.T
    )

    phi, middle, trace = a6_trace_data()
    expected_trace = zeros(3).row_join(presentation).col_join(
        presentation.row_join(zeros(3))
    )
    unit = a6_trace_unit()

    cusp_substitution = {P: 0, Q: 0}
    cusp_phi = phi.subs(cusp_substitution)
    cusp_trace = trace.subs(cusp_substitution)
    cusp_basis = _cusp_basis()
    cusp_expected_trace = diag(5, 0, 0, 0, 0, 1)
    cusp_unit_coordinates = Matrix([1, 0, 0, 0, 0, 1])
    cusp_middle, cusp_section = _cusp_frobenius_data()

    collision_substitution = {P: COLLISION_P, Q: 1}
    collision_phi = phi.subs(collision_substitution)
    collision_trace = trace.subs(collision_substitution)
    collision_basis = _collision_basis()
    collision_expected_trace = diag(3, 0, 0, 3, 0, 0)
    collision_unit_coordinates = Matrix([1, 0, 0, 1, 0, 0])
    collision_middle, collision_section = _collision_frobenius_data()

    family = hostile_normalization_family()
    family_polynomial = Poly(family, FIBER_VARIABLE)
    cusp_family = expand(family.subs(PARAMETER, 0))
    expected_cusp_family = FIBER_VARIABLE**5 * (FIBER_VARIABLE - 1)
    expected_collision_family = (
        FIBER_VARIABLE**3 * (FIBER_VARIABLE - 1) ** 3
    )
    collision_family_difference = Poly(
        expand(family - expected_collision_family),
        PARAMETER,
    ).rem(Poly(CYCLOTOMIC_FIVE, PARAMETER))

    return A6LocalMultiplicationCertificate(
        presentation_determinant=factor(presentation.det()),
        adjugate_identity=_entries(
            (presentation * adjugate - BRANCH_EQUATION * identity_three)
            .applyfunc(expand)
        ),
        parametrized_branch_equation=factor(
            BRANCH_EQUATION.subs(parameter_substitution)
        ),
        normalization_kernel=_entries(
            (parametrized_presentation * normalization_generators)
            .applyfunc(factor)
        ),
        normalized_adjugate_factor_identity=_entries(
            (
                adjugate.subs(parameter_substitution)
                - expected_normalized_adjugate
            ).applyfunc(factor)
        ),
        trace_factor_identity=_entries(
            (trace - phi.T * middle * phi).applyfunc(expand)
        )
        + _entries((trace - expected_trace).applyfunc(expand)),
        trace_unit_norm=factor((unit.T * trace * unit)[0]),
        cusp_basis_determinant=simplify(cusp_basis.det()),
        cusp_gram_identity=_simplified_entries(
            cusp_basis.T * cusp_trace * cusp_basis - cusp_expected_trace
        ),
        cusp_unit_identity=_simplified_entries(
            cusp_basis * cusp_unit_coordinates - unit
        ),
        cusp_middle_determinant=factor(cusp_middle.det()),
        cusp_middle_trace_identity=_simplified_entries(
            cusp_section.T * cusp_middle * cusp_section
            - cusp_expected_trace
        ),
        cusp_phi_kernel_columns=(
            _simplified_entries(cusp_phi * cusp_basis[:, 3]),
            _simplified_entries(cusp_phi * cusp_basis[:, 4]),
        ),
        collision_basis_determinant=_reduce_at_collision(
            collision_basis.det()
        ),
        collision_gram_identity=_collision_entries(
            collision_basis.T
            * collision_trace
            * collision_basis
            - collision_expected_trace
        ),
        collision_unit_identity=_collision_entries(
            collision_basis * collision_unit_coordinates - unit
        ),
        collision_middle_determinant=factor(collision_middle.det()),
        collision_middle_trace_identity=_entries(
            (
                collision_section.T
                * collision_middle
                * collision_section
                - collision_expected_trace
            ).applyfunc(expand)
        ),
        collision_phi_kernel_columns=(
            _collision_entries(collision_phi * collision_basis[:, 2]),
            _collision_entries(collision_phi * collision_basis[:, 5]),
        ),
        hostile_family_degree=int(family_polynomial.degree()),
        hostile_family_leading_coefficient=factor(family_polynomial.LC()),
        hostile_family_cusp_identity=factor(
            cusp_family - expected_cusp_family
        ),
        hostile_family_collision_identity=factor(
            collision_family_difference.as_expr()
        ),
    )


def main() -> int:
    """Print the exact stopping certificate and fail on any regression."""

    certificate = exact_a6_local_multiplication_certificate()
    print(
        "A6 local multiplication:",
        {
            "detK": certificate.presentation_determinant,
            "unit norm": certificate.trace_unit_norm,
            "cusp basis determinant": certificate.cusp_basis_determinant,
            "collision basis determinant": (
                certificate.collision_basis_determinant
            ),
            "hostile family degree": certificate.hostile_family_degree,
        },
    )
    print(f"certificate verified: {certificate.verified}")
    print("claim boundary: local fibers only; no global algebra or Keller map")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
