"""Reproduce the conditional delta-ten C2 + Q0 + 2N certificate.

This checker derives the compatibility quartic, performs exact full-localizer
rank and singular saturations over QQ, and regenerates the rational member's
implicit curve, singular scheme, and affine van Kamp presentation.  The
proper-isotopy propagation across the clean family remains a separate
theorem-level dependency.
"""

from itertools import combinations
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sage.schemes.curves.zariski_vankampen import fundamental_group

from scripts.a6_delta_ten_contact_quadruple import (
    CONTACT_QUADRUPLE_AUGMENTED_MATRIX,
    CONTACT_QUADRUPLE_COEFFICIENT_MATRIX,
    CONTACT_QUADRUPLE_COMPATIBILITY,
    CONTACT_QUADRUPLE_RELATIONS,
    CONTACT_QUADRUPLE_VALID_LOCALIZER,
    CONTACT_SUM,
    EXPECTED_AUGMENTED_DETERMINANT,
    SAMPLE_IMPLICIT,
    SAMPLE_NODE_X_POLYNOMIAL,
)
from scripts.a6_delta_ten_generic import KAPPA
from scripts.a6_delta_ten_quadruple import FIBER_VALUE


Base = PolynomialRing(QQ, ["k", "h", "u"], order="degrevlex")
k, h, u = Base.gens()


def as_base(expression):
    """Convert an exact SymPy expression into the Sage base ring."""

    return Base(str(expression).replace("**", "^"))


coefficient_matrix = matrix(
    Base,
    5,
    4,
    [as_base(entry) for entry in CONTACT_QUADRUPLE_COEFFICIENT_MATRIX],
)
augmented_matrix = matrix(
    Base,
    5,
    5,
    [as_base(entry) for entry in CONTACT_QUADRUPLE_AUGMENTED_MATRIX],
)
compatibility = as_base(CONTACT_QUADRUPLE_COMPATIBILITY)
valid_localizer = as_base(CONTACT_QUADRUPLE_VALID_LOCALIZER)
expected_determinant = as_base(EXPECTED_AUGMENTED_DETERMINANT)

assert augmented_matrix.det() == expected_determinant
compatibility_factorization = compatibility.factor()
assert len(compatibility_factorization) == 1
assert compatibility_factorization[0][0] == compatibility
assert compatibility_factorization[0][1] == 1

sample_base = {k: -4, h: 1, u: 1}
assert compatibility.subs(sample_base) == 0
assert tuple(
    compatibility.derivative(variable).subs(sample_base)
    for variable in (k, h, u)
) == (4, 0, 8)
assert valid_localizer.subs(sample_base) == -45563904
sample_coefficient = coefficient_matrix.subs(sample_base)
sample_augmented = augmented_matrix.subs(sample_base)
assert sample_coefficient.rank() == 4
assert sample_augmented.rank() == 4
assert tuple(
    sample_coefficient.matrix_from_rows(
        [row for row in range(5) if row != omitted]
    ).det()
    for omitted in range(5)
) == (-12288, 9728, -1536, -192, -32)


# The projective quartic has four smooth points at infinity and exactly two
# affine ordinary nodes, both on split fibers.  Thus p_a=3 and g=3-2=1.
ProjectiveRing.<K, U, W> = QQ[]
projective_compatibility = (
    K^3*U + 5*K^2*U^2 + 4*K*U^3
    + (3*K^2 + 16*K*U + 12*U^2)*W^2 + 4*W^4
)
assert projective_compatibility(W=1, K=k, U=u) == compatibility
infinity_points = ((0, 1, 0), (-1, 1, 0), (-4, 1, 0), (1, 0, 0))
for point in infinity_points:
    assert projective_compatibility(*point) == 0
    assert any(
        derivative(*point) != 0
        for derivative in (
            projective_compatibility.derivative(K),
            projective_compatibility.derivative(U),
            projective_compatibility.derivative(W),
        )
    )

Affine.<x_shift, y_shift> = QQ[]
for point in ((-2, 1), (2, -1)):
    shifted = Affine(
        compatibility(
            k=x_shift + point[0],
            h=0,
            u=y_shift + point[1],
        )
    )
    tangent_cone = sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            shifted.coefficients(),
            shifted.monomials(),
            strict=True,
        )
        if monomial.degree() == 2
    )
    assert tangent_cone == 2*x_shift^2 + 8*y_shift^2
    tangent_hessian = matrix([
        [
            tangent_cone.derivative(first).derivative(second)
            for second in (x_shift, y_shift)
        ]
        for first in (x_shift, y_shift)
    ])
    assert tangent_hessian.det() == 64

affine_singular = Base.ideal(
    [
        compatibility,
        compatibility.derivative(k),
        compatibility.derivative(u),
    ]
)
assert affine_singular.dimension() == 1
singular_radical = affine_singular.radical()
split_singular_lines = (
    Base.ideal([k + 2, u - 1]),
    Base.ideal([k - 2, u + 1]),
)
assert singular_radical == split_singular_lines[0].intersection(
    split_singular_lines[1]
)
assert (4 - 1) * (4 - 2) // 2 == 3
assert 3 - len(split_singular_lines) == 1


valid_factors = (
    k,
    k - 2,
    k + 2,
    h,
    u,
    k + 2*u,
    u^2 + k*u + 1,
    2*u^2 + 3*k*u + 4,
    256*h^2 + 27*h*k^4 - 144*h*k^2 + 128*h - 4*k^2 + 16,
    h*(k + 2*u)^2 + u*(u + k)*(u^2 + k*u + 1)^2,
)
assert prod(valid_factors) == valid_localizer


def strip_valid_units(polynomial):
    """Remove maximal valid-localizer powers and rational content."""

    reduced = Base(polynomial)
    if reduced == 0:
        return reduced
    for factor in valid_factors:
        while True:
            quotient, remainder = reduced.quo_rem(factor)
            if remainder != 0:
                break
            reduced = quotient
    reduced = Base(reduced / reduced.content())
    if reduced.leading_coefficient() < 0:
        reduced = -reduced
    return reduced


def normalized_minors(target, size):
    """Return distinct nonzero minors modulo valid-chart units."""

    normalized = []
    for rows in combinations(range(target.nrows()), size):
        for columns in combinations(range(target.ncols()), size):
            minor = strip_valid_units(
                target.matrix_from_rows_and_columns(rows, columns).det()
            )
            if minor != 0 and minor not in normalized:
                normalized.append(minor)
    return normalized


def exact_modular_saturation(ideal, factor):
    """Run Singular exact modular-reconstruction saturation over QQ."""

    singular.lib("modquotient.lib")
    answer = singular.modSat(
        singular(ideal),
        singular(Base.ideal([factor])),
    )
    return answer[1].sage(), answer[2].sage()


valid_singular, valid_singular_exponent = exact_modular_saturation(
    affine_singular,
    valid_localizer,
)
assert valid_singular.is_one()
assert valid_singular_exponent == 1
print("[1/4] compatibility smoothness saturation complete", flush=True)

coefficient_rank_three = Base.ideal(
    [compatibility, *normalized_minors(coefficient_matrix, 4)]
)
augmented_rank_three = Base.ideal(
    [compatibility, *normalized_minors(augmented_matrix, 4)]
)
coefficient_rank_two = Base.ideal(
    [compatibility, *normalized_minors(coefficient_matrix, 3)]
)
augmented_rank_two = Base.ideal(
    [compatibility, *normalized_minors(augmented_matrix, 3)]
)

coefficient_rank_three_valid, coefficient_rank_three_exponent = (
    exact_modular_saturation(coefficient_rank_three, valid_localizer)
)
augmented_rank_three_valid, augmented_rank_three_exponent = (
    exact_modular_saturation(augmented_rank_three, valid_localizer)
)
coefficient_rank_two_valid, coefficient_rank_two_exponent = (
    exact_modular_saturation(coefficient_rank_two, valid_localizer)
)
augmented_rank_two_valid, augmented_rank_two_exponent = (
    exact_modular_saturation(augmented_rank_two, valid_localizer)
)
assert coefficient_rank_three_valid.is_one()
assert augmented_rank_three_valid.is_one()
assert coefficient_rank_two_valid.is_one()
assert augmented_rank_two_valid.is_one()
assert (
    coefficient_rank_three_exponent,
    augmented_rank_three_exponent,
    coefficient_rank_two_exponent,
    augmented_rank_two_exponent,
) == (1, 1, 0, 0)
print("[2/4] full-valid rank saturations complete", flush=True)


# Regenerate the exact member over ZZ with X=P and Y=Q.
Source.<source_X, source_Y, t> = ZZ[]
source_x = t^4 - 4*t^3 + t^2
source_y = t^9 - 6*t^8 + 13*t^7 - 19*t^6 + 7*t^5
raw_implicit = (source_X - source_x).resultant(source_Y - source_y, t)
coefficient_content = gcd(
    abs(coefficient) for coefficient in raw_implicit.coefficients()
)
primitive_implicit = Source(raw_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(source_Y^4) < 0:
    primitive_implicit = -primitive_implicit

Plane.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], Plane)
regenerated_curve = source_to_plane(primitive_implicit)
curve = Plane(str(SAMPLE_IMPLICIT).replace("**", "^"))
assert coefficient_content == 1
assert regenerated_curve == curve
assert len(curve.factor()) == 1
ParameterRing.<z> = QQ[]
assert curve(
    X=z^4 - 4*z^3 + z^2,
    Y=z^9 - 6*z^8 + 13*z^7 - 19*z^6 + 7*z^5,
) == 0

jacobian = Plane.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 18
assert jacobian.radical().vector_space_dimension() == 5

node_ideal = Plane.ideal([
    X^2 + 564*X - 476,
    8491*X + 2*Y - 7154,
])
quadruple_ideal = Plane.ideal([X - 1, Y + 1])
contact_ideal = Plane.ideal([X - 3, Y + 7])
cusp_ideal = Plane.ideal([X, Y])
expected_radicals = (node_ideal, quadruple_ideal, contact_ideal, cusp_ideal)
component_data = []
for component in jacobian.primary_decomposition():
    radical = component.radical()
    lengths = (
        component.vector_space_dimension(),
        radical.vector_space_dimension(),
    )
    component_data.append((radical, lengths))
assert tuple(lengths for _, lengths in component_data) == (
    (2, 2),
    (9, 1),
    (3, 1),
    (4, 1),
)
assert tuple(radical for radical, _ in component_data) == expected_radicals
assert component_data[0][0] == component_data[0][0].radical()

NodeRing.<node_x> = QQ[]
node_x_polynomial = NodeRing(
    str(SAMPLE_NODE_X_POLYNOMIAL)
    .replace("X", "node_x")
    .replace("**", "^")
)
assert node_x_polynomial.is_irreducible()
assert node_x_polynomial.discriminant() == 320000
node_to_plane = NodeRing.hom([X], Plane)
assert node_ideal.reduce(node_to_plane(node_x_polynomial)) == 0
print("[3/4] sample curve and singular scheme complete", flush=True)


group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert len(actual_relations) == 9
assert actual_relations == CONTACT_QUADRUPLE_RELATIONS

isomorphism = group.simplification_isomorphism()
infinite_cyclic = isomorphism.codomain()
assert len(infinite_cyclic.gens()) == 1
assert len(infinite_cyclic.relations()) == 0
assert all(
    isomorphism(generator) == infinite_cyclic.gen(0)
    for generator in group.gens()
)
section = isomorphism.section()
assert section(infinite_cyclic.gen(0)) == group.gen(0)
assert isomorphism(section(infinite_cyclic.gen(0))) == infinite_cyclic.gen(0)
print("[4/4] affine complement reconstruction complete", flush=True)

print("PASS")
print("compatibility quartic irreducible over QQ: True")
print("valid compatibility singular saturation: unit, exponent", (
    valid_singular_exponent
))
print("coefficient/augmented rank <= 3 saturations: units, exponents", (
    coefficient_rank_three_exponent,
    augmented_rank_three_exponent,
))
print("coefficient/augmented rank <= 2 saturations: units, exponents", (
    coefficient_rank_two_exponent,
    augmented_rank_two_exponent,
))
print("Jacobian component lengths:", tuple(
    lengths for _, lengths in component_data
))
print("affine complement: Z")
