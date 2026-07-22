"""Reproduce the conditional delta-ten ``C3 + T111 + 4N`` certificate.

This checker independently derives the compatibility hypersurface from the
five incidence equations, performs exact full-localizer singular and rank
saturations over ``QQ``, verifies a rank-two coefficient-image differential,
and regenerates the rational member's singular scheme and affine van Kamp
presentation.  Propagation over the clean family remains a separate
theorem-level dependency.
"""

from itertools import combinations
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sage.schemes.curves.zariski_vankampen import fundamental_group

from scripts.a6_delta_ten_contact_three_triple import (
    CONTACT_THREE_SUM,
    CONTACT_THREE_TRIPLE_AUGMENTED_MATRIX,
    CONTACT_THREE_TRIPLE_COEFFICIENT_MATRIX,
    CONTACT_THREE_TRIPLE_COMPATIBILITY,
    CONTACT_THREE_TRIPLE_EQUATIONS,
    CONTACT_THREE_TRIPLE_IMPLICIT,
    CONTACT_THREE_TRIPLE_NODE_X_POLYNOMIAL,
    CONTACT_THREE_TRIPLE_RELATIONS,
    CONTACT_THREE_TRIPLE_VALID_LOCALIZER,
    CONTACT_TRIPLE_OVERLAP_B,
    CONTACT_TRIPLE_OVERLAP_O,
    EXPECTED_AUGMENTED_DETERMINANT,
    SAMPLE_BASE_TANGENTS,
    SAMPLE_COEFFICIENT_TANGENTS,
    SAMPLE_TRIPLE_SLOPE_POLYNOMIAL,
    TRIPLE_OMITTED_ROOT,
)


Base = PolynomialRing(QQ, ["k", "e", "w"], order="degrevlex")
k, e, w = Base.gens()


def as_base(expression):
    """Convert an exact SymPy expression into the Sage base ring."""

    return Base(str(expression).replace("**", "^"))


coefficient_matrix = matrix(
    Base,
    5,
    4,
    [as_base(entry) for entry in CONTACT_THREE_TRIPLE_COEFFICIENT_MATRIX],
)
augmented_matrix = matrix(
    Base,
    5,
    5,
    [as_base(entry) for entry in CONTACT_THREE_TRIPLE_AUGMENTED_MATRIX],
)
compatibility = as_base(CONTACT_THREE_TRIPLE_COMPATIBILITY)
valid_localizer = as_base(CONTACT_THREE_TRIPLE_VALID_LOCALIZER)
expected_determinant = as_base(EXPECTED_AUGMENTED_DETERMINANT)

assert augmented_matrix.det() == expected_determinant
known_factor = (
    2
    * (e^2 + k*e + 1)^2
    * (w^2 + k*w + 1)^3
    * as_base(CONTACT_TRIPLE_OVERLAP_O)^3
)
derived_residual, derived_remainder = augmented_matrix.det().quo_rem(known_factor)
assert derived_remainder == 0
assert derived_residual == compatibility
assert len(compatibility.monomials()) == 62
assert compatibility.degree() == 9
assert compatibility.degree(k) == 6
assert compatibility.degree(e) == 3
assert compatibility.degree(w) == 7
compatibility_factorization = compatibility.factor()
assert len(compatibility_factorization) == 1
assert compatibility_factorization[0][0] == compatibility
assert compatibility_factorization[0][1] == 1

sample_base = {k: -4, e: QQ(-1)/2, w: 1}
assert compatibility.subs(sample_base) == 0
assert tuple(
    compatibility.derivative(variable).subs(sample_base)
    for variable in (k, e, w)
) == (-84, -112, 0)
assert valid_localizer.subs(sample_base) == 7223580
sample_coefficient = coefficient_matrix.subs(sample_base)
sample_augmented = augmented_matrix.subs(sample_base)
assert sample_coefficient.rank() == 4
assert sample_augmented.rank() == 4
assert tuple(
    sample_coefficient.matrix_from_rows(rows).det()
    for rows in combinations(range(5), 4)
) == (-2504320, 3333120, 0, QQ(29575)/4, QQ(159705)/4)
sample_solution = sample_coefficient.solve_right(
    vector(QQ, [0, 48, 528, QQ(2409407)/256, QQ(-200427)/16])
)
assert tuple(sample_solution) == (
    QQ(39)/2,
    QQ(-409)/8,
    QQ(109)/4,
    QQ(-31)/4,
)


triple_discriminant = (
    16*e^4 + 8*e^3*k - 5*e^2*k^2 + 16*e^2
    + 3*e*k^3 - 12*e*k - k^2 + 4
)
overlap_b = as_base(CONTACT_TRIPLE_OVERLAP_B)
overlap_o = as_base(CONTACT_TRIPLE_OVERLAP_O)
valid_factors = (
    k,
    k - 2,
    k + 2,
    e,
    w,
    k + 2*w,
    w^2 + k*w + 1,
    2*w^2 + 3*k*w + 4,
    e^2 + k*e + 1,
    4*e^2 + 3*k*e + 2,
    triple_discriminant,
    overlap_b,
    overlap_o,
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
    """Run Singular exact modular-reconstruction saturation over ``QQ``."""

    singular.lib("modquotient.lib")
    answer = singular.modSat(
        singular(ideal),
        singular(Base.ideal([factor])),
    )
    return answer[1].sage(), answer[2].sage()


singular_ideal = Base.ideal([
    compatibility,
    compatibility.derivative(k),
    compatibility.derivative(e),
    compatibility.derivative(w),
])
valid_singular, valid_singular_exponent = exact_modular_saturation(
    singular_ideal,
    valid_localizer,
)
assert valid_singular.is_one()
assert valid_singular_exponent == 1
print("[1/5] compatibility smoothness saturation complete", flush=True)

coefficient_four_minors = normalized_minors(coefficient_matrix, 4)
augmented_four_minors = normalized_minors(augmented_matrix, 4)
coefficient_three_minors = normalized_minors(coefficient_matrix, 3)
augmented_three_minors = normalized_minors(augmented_matrix, 3)
assert tuple(map(len, (
    coefficient_four_minors,
    augmented_four_minors,
    coefficient_three_minors,
    augmented_three_minors,
))) == (5, 25, 40, 100)

coefficient_rank_three = Base.ideal([compatibility, *coefficient_four_minors])
augmented_rank_three = Base.ideal([compatibility, *augmented_four_minors])
coefficient_rank_two = Base.ideal([compatibility, *coefficient_three_minors])
augmented_rank_two = Base.ideal([compatibility, *augmented_three_minors])

coefficient_rank_three_valid, coefficient_rank_three_exponent = (
    exact_modular_saturation(coefficient_rank_three, valid_localizer)
)
assert not coefficient_rank_three_valid.is_one()
assert coefficient_rank_three_valid.dimension() == 1
assert coefficient_rank_three_exponent == 1
rank_three_hilbert = coefficient_rank_three_valid.homogenize().hilbert_polynomial()
assert rank_three_hilbert.degree() == 1
assert rank_three_hilbert[1] == 14
assert rank_three_hilbert[0] == -21
print("[2/5] coefficient rank-three curve complete", flush=True)

augmented_rank_three_valid, augmented_rank_three_exponent = (
    exact_modular_saturation(augmented_rank_three, valid_localizer)
)
assert augmented_rank_three_valid.is_one()
assert augmented_rank_three_exponent == 1
print("[3/5] augmented rank-three incompatibility complete", flush=True)

coefficient_rank_two_valid, coefficient_rank_two_exponent = (
    exact_modular_saturation(coefficient_rank_two, valid_localizer)
)
augmented_rank_two_valid, augmented_rank_two_exponent = (
    exact_modular_saturation(augmented_rank_two, valid_localizer)
)
assert coefficient_rank_two_valid.is_one()
assert augmented_rank_two_valid.is_one()
assert (coefficient_rank_two_exponent, augmented_rank_two_exponent) == (1, 1)
print("[4/5] redundant rank-two saturations complete", flush=True)


# Independently differentiate the five original equations at the sample.
Total = PolynomialRing(QQ, ["k", "e", "w", "a", "b", "c", "d"])
tk, te, tw, ta, tb, tc, td = Total.gens()


def as_total(expression):
    """Convert an exact SymPy expression into the total incidence ring."""

    return Total(str(expression).replace("**", "^"))


total_equations = [as_total(equation) for equation in CONTACT_THREE_TRIPLE_EQUATIONS]
sample_total = {
    tk: -4,
    te: QQ(-1)/2,
    tw: 1,
    ta: QQ(39)/2,
    tb: QQ(-409)/8,
    tc: QQ(109)/4,
    td: QQ(-31)/4,
}
base_jacobian = matrix(QQ, [
    [equation.derivative(variable).subs(sample_total) for variable in (tk, te, tw)]
    for equation in total_equations
])
parameter_jacobian = matrix(QQ, [
    [equation.derivative(variable).subs(sample_total) for variable in (ta, tb, tc, td)]
    for equation in total_equations
])
for base_tangent, coefficient_tangent in zip(
    SAMPLE_BASE_TANGENTS,
    SAMPLE_COEFFICIENT_TANGENTS,
    strict=True,
):
    sage_base_tangent = vector(QQ, [QQ(str(value)) for value in base_tangent])
    sage_coefficient_tangent = vector(
        QQ,
        [QQ(str(value)) for value in coefficient_tangent],
    )
    assert (
        base_jacobian * sage_base_tangent
        + parameter_jacobian * sage_coefficient_tangent
    ).is_zero()
image_tangents = matrix(QQ, [
    tuple(
        QQ(str(value))
        for value in (base_tangent[0], *coefficient_tangent)
    )
    for base_tangent, coefficient_tangent in zip(
        SAMPLE_BASE_TANGENTS,
        SAMPLE_COEFFICIENT_TANGENTS,
        strict=True,
    )
])
assert image_tangents.rank() == 2


# Regenerate the exact member over ZZ with X=P and Y=8Q.
Source.<source_X, source_Y, t> = ZZ[]
source_x = t^4 - 4*t^3 + t^2
source_y = 8*t^9 - 62*t^8 + 218*t^7 - 409*t^6 + 156*t^5
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
curve = Plane(str(CONTACT_THREE_TRIPLE_IMPLICIT).replace("**", "^"))
assert coefficient_content == 1
assert regenerated_curve == curve
assert len(curve.factor()) == 1
ParameterRing.<z, m> = QQ[]
assert curve(
    X=z^4 - 4*z^3 + z^2,
    Y=8*z^9 - 62*z^8 + 218*z^7 - 409*z^6 + 156*z^5,
) == 0
assert curve(X=QQ(13)/16, Y=QQ(-2197)/128) == 0
assert curve(X=3, Y=-199) == 0
assert curve(X=0, Y=0) == 0

# The ordinary-triple slopes are W(t)/(t-e), not merely W(t).  This
# corrected eliminant is deliberately checked independently here.
triple_cubic = z^3 - QQ(9)/2*z^2 + QQ(13)/4*z - QQ(13)/8
triple_quotient = (
    z^6 - QQ(13)/4*z^5 + QQ(75)/8*z^4 + QQ(13)/4*z^3
    - QQ(13)/8*z^2 - QQ(169)/64*z - QQ(169)/128
)
slope_eliminant = triple_cubic.resultant(
    triple_quotient - m*(z + QQ(1)/2),
    z,
)
stored_slope = ParameterRing(
    str(SAMPLE_TRIPLE_SLOPE_POLYNOMIAL)
    .replace("m_three_triple", "m")
    .replace("**", "^")
)
assert 262144*slope_eliminant + stored_slope == 0
assert stored_slope.discriminant(m) == -89062908728555597524369408000000

jacobian = Plane.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 17
assert jacobian.radical().vector_space_dimension() == 7

triple_ideal = Plane.ideal([16*X - 13, 128*Y + 2197])
contact_ideal = Plane.ideal([X - 3, Y + 199])
cusp_ideal = Plane.ideal([X, Y])
node_component = None
component_data = []
for component in jacobian.primary_decomposition():
    radical = component.radical()
    lengths = (
        component.vector_space_dimension(),
        radical.vector_space_dimension(),
    )
    component_data.append((radical, lengths))
    if radical not in (triple_ideal, contact_ideal, cusp_ideal):
        node_component = component
        assert lengths == (4, 4)
    elif radical == triple_ideal:
        assert lengths == (4, 1)
    elif radical == contact_ideal:
        assert lengths == (5, 1)
    else:
        assert radical == cusp_ideal
        assert lengths == (4, 1)
assert sorted(lengths for _, lengths in component_data) == sorted([
    (4, 4),
    (4, 1),
    (5, 1),
    (4, 1),
])
assert node_component is not None
assert node_component == node_component.radical()

NodeRing.<node_x> = QQ[]
node_x_polynomial = NodeRing(
    str(CONTACT_THREE_TRIPLE_NODE_X_POLYNOMIAL)
    .replace("X", "node_x")
    .replace("**", "^")
)
assert node_x_polynomial.is_irreducible()
assert node_x_polynomial.discriminant() == -4321892177659568
node_to_plane = NodeRing.hom([X], Plane)
assert node_component.reduce(node_to_plane(node_x_polynomial)) == 0

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert CONTACT_THREE_TRIPLE_RELATIONS
assert actual_relations == CONTACT_THREE_TRIPLE_RELATIONS

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
print("[5/5] sample singular scheme and complement complete", flush=True)

print("PASS")
print("compatibility residual terms/degrees:", (
    len(compatibility.monomials()),
    compatibility.degree(),
    compatibility.degree(k),
    compatibility.degree(e),
    compatibility.degree(w),
))
print("compatibility residual irreducible over QQ: True")
print("valid singular saturation: unit, exponent", valid_singular_exponent)
print("coefficient rank <= 3 curve: degree 14, exponent", (
    coefficient_rank_three_exponent
))
print("augmented rank <= 3 saturation: unit, exponent", (
    augmented_rank_three_exponent
))
print("rank <= 2 saturations: units, exponents", (
    coefficient_rank_two_exponent,
    augmented_rank_two_exponent,
))
print("coefficient-image tangent rank: 2")
print("Jacobian component lengths:", tuple(
    lengths for _, lengths in component_data
))
print("affine complement: Z")
