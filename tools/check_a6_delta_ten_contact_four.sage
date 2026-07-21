"""Reproduce the conditional delta-ten ``C4 + 6N`` certificate in Sage.

The first half is a determinantal saturation over the nonsplit unordered-pair
base.  It proves that the residual determinant curve cannot hide a
higher-dimensional coefficient fiber.  The second half independently
regenerates the rational member, its affine singular scheme, and its van Kamp
presentation.  No topology-propagation claim is made here.
"""

from itertools import combinations
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sage.schemes.curves.zariski_vankampen import fundamental_group

from scripts.a6_delta_ten_contact_four import (
    CONTACT_FOUR_AUGMENTED_MATRIX,
    CONTACT_FOUR_COEFFICIENT_MATRIX,
    CONTACT_FOUR_IMPLICIT,
    CONTACT_FOUR_NODE_X_POLYNOMIAL,
    CONTACT_FOUR_RELATIONS,
    CONTACT_FOUR_RESIDUAL_FACTOR,
    CONTACT_FOUR_SCALED_Q,
)
from scripts.a6_delta_ten_generic import (
    PAIR_DENOMINATOR,
    PAIR_DIAGONAL_FACTOR,
    PAIR_QUADRATIC,
)


Base = PolynomialRing(QQ, ["k", "s"], order="degrevlex")
k, s = Base.gens()


def as_base(expression):
    """Convert one exact SymPy expression into the Sage base ring."""

    return Base(str(expression).replace("**", "^"))


coefficient_matrix = matrix(
    Base,
    4,
    4,
    [as_base(entry) for entry in CONTACT_FOUR_COEFFICIENT_MATRIX],
)
augmented_matrix = matrix(
    Base,
    4,
    5,
    [as_base(entry) for entry in CONTACT_FOUR_AUGMENTED_MATRIX],
)
pair_denominator = as_base(PAIR_DENOMINATOR)
pair_quadratic = as_base(PAIR_QUADRATIC)
pair_diagonal = as_base(PAIR_DIAGONAL_FACTOR)
residual_factor = as_base(CONTACT_FOUR_RESIDUAL_FACTOR)

# Remove row factors that are nonzero on the valid chart.  The H row has a
# pair-denominator factor only in its coefficient entries, not in its
# inhomogeneous entry, so it is normalized only in the coefficient matrix.
for column in range(4):
    quotient, remainder = coefficient_matrix[0, column].quo_rem(pair_denominator)
    assert remainder == 0
    coefficient_matrix[0, column] = quotient

for target in (coefficient_matrix, augmented_matrix):
    for column in range(target.ncols()):
        for row, divisor in ((2, Base(2)), (3, Base(6))):
            quotient, remainder = target[row, column].quo_rem(divisor)
            assert remainder == 0
            target[row, column] = quotient

assert coefficient_matrix.det() == (
    2 * pair_denominator^3 * pair_quadratic^5 * residual_factor
)
assert len(residual_factor.factor()) == 1

valid_localizer = (
    k
    * (k^2 - 4)
    * s
    * pair_denominator
    * pair_quadratic
    * pair_diagonal
)
valid_ideal = Base.ideal([valid_localizer])

# Rank <= 2 means all sixteen 3 x 3 coefficient minors vanish.  The raw
# locus contains invalid chart curves; saturation removes all of it.
coefficient_rank_two_minors = [
    coefficient_matrix.matrix_from_rows_and_columns(rows, columns).det()
    for rows in combinations(range(4), 3)
    for columns in combinations(range(4), 3)
]
coefficient_rank_two_ideal = Base.ideal(coefficient_rank_two_minors)
coefficient_rank_two_valid, coefficient_rank_two_exponent = (
    coefficient_rank_two_ideal.saturation(valid_ideal)
)
assert coefficient_rank_two_ideal.dimension() == 1
assert coefficient_rank_two_valid.is_one()
assert coefficient_rank_two_exponent == 2

# On the residual determinant curve, consistency is equivalent to vanishing
# of the four augmented maximal minors that use the constant column.  The
# valid compatible scheme is finite of exact length ten.
compatibility_minors = [
    augmented_matrix.matrix_from_columns(columns).det()
    for columns in combinations(range(5), 4)
    if 4 in columns
]
assert all(gcd(residual_factor, minor) == 1 for minor in compatibility_minors)
compatibility_ideal = Base.ideal([residual_factor] + compatibility_minors)
compatibility_valid, compatibility_exponent = compatibility_ideal.saturation(
    valid_ideal
)
assert compatibility_ideal.dimension() == 0
assert not compatibility_valid.is_one()
assert compatibility_valid.dimension() == 0
assert compatibility_valid.vector_space_dimension() == 10
assert compatibility_exponent == 4

# These augmented rank strata are a deliberately redundant hostile check.
# In particular, there is no valid compatible base at which a coefficient
# fiber could jump to dimension at least two.
augmented_rank_two_minors = [
    augmented_matrix.matrix_from_rows_and_columns(rows, columns).det()
    for rows in combinations(range(4), 3)
    for columns in combinations(range(5), 3)
]
augmented_rank_two_ideal = Base.ideal(augmented_rank_two_minors)
augmented_rank_two_valid, augmented_rank_two_exponent = (
    augmented_rank_two_ideal.saturation(valid_ideal)
)
assert augmented_rank_two_ideal.dimension() == 1
assert augmented_rank_two_valid.is_one()
assert augmented_rank_two_exponent == 2

augmented_rank_one_minors = [
    augmented_matrix.matrix_from_rows_and_columns(rows, columns).det()
    for rows in combinations(range(4), 2)
    for columns in combinations(range(5), 2)
]
augmented_rank_one_ideal = Base.ideal(augmented_rank_one_minors)
augmented_rank_one_valid, augmented_rank_one_exponent = (
    augmented_rank_one_ideal.saturation(valid_ideal)
)
assert augmented_rank_one_ideal.dimension() == 0
assert augmented_rank_one_valid.is_one()
assert augmented_rank_one_exponent == 1


# Regenerate the exact rational member over ZZ after scaling Q by fourteen.
Source.<u, v, t> = ZZ[]
source_x = t^2 + t^3 + t^4
source_y = Source(str(CONTACT_FOUR_SCALED_Q).replace("**", "^"))
assert source_y == 18*t^5 + 27*t^6 + 42*t^7 + 24*t^8 + 14*t^9
raw_implicit = (u - source_x).resultant(v - source_y, t)
coefficient_content = gcd(
    abs(coefficient) for coefficient in raw_implicit.coefficients()
)
primitive_implicit = Source(raw_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(v^4) < 0:
    primitive_implicit = -primitive_implicit

R.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], R)
regenerated_curve = source_to_plane(primitive_implicit)
curve = R(str(CONTACT_FOUR_IMPLICIT).replace("**", "^"))
assert coefficient_content == 1
assert regenerated_curve == curve

ParameterRing.<z> = QQ[]
assert curve(
    X=z^2 + z^3 + z^4,
    Y=18*z^5 + 27*z^6 + 42*z^7 + 24*z^8 + 14*z^9,
) == 0
assert curve(X=0, Y=-1) == 0

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert len(actual_relations) == 10
assert actual_relations == CONTACT_FOUR_RELATIONS

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

jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 17
assert jacobian.radical().vector_space_dimension() == 8

origin_ideal = R.ideal([X, Y])
contact_ideal = R.ideal([X, Y + 1])
components = tuple(
    (component, component.radical())
    for component in jacobian.primary_decomposition()
)
origin_components = tuple(
    component for component, radical in components if radical == origin_ideal
)
contact_components = tuple(
    component for component, radical in components if radical == contact_ideal
)
node_components = tuple(
    component
    for component, radical in components
    if radical != origin_ideal and radical != contact_ideal
)
assert len(origin_components) == 1
assert len(contact_components) == 1
assert len(node_components) == 1
origin_component = origin_components[0]
contact_component = contact_components[0]
node_component = node_components[0]
assert origin_component.vector_space_dimension() == 4
assert origin_component.radical().vector_space_dimension() == 1
assert contact_component.vector_space_dimension() == 7
assert contact_component.radical().vector_space_dimension() == 1
assert node_component.vector_space_dimension() == 6
assert node_component == node_component.radical()

NodeRing.<x> = QQ[]
node_x_polynomial = NodeRing(
    str(CONTACT_FOUR_NODE_X_POLYNOMIAL).replace("X", "x").replace("**", "^")
)
assert node_x_polynomial.is_irreducible()
assert node_x_polynomial.discriminant() != 0
node_to_plane = NodeRing.hom([X], R)
assert node_component.reduce(node_to_plane(node_x_polynomial)) == 0

print("PASS")
print("residual factor irreducible:", residual_factor)
print(
    "valid coefficient rank <= 2 saturation: unit, exponent",
    coefficient_rank_two_exponent,
)
print(
    "valid compatibility scheme: dimension",
    compatibility_valid.dimension(),
    "length",
    compatibility_valid.vector_space_dimension(),
    "exponent",
    compatibility_exponent,
)
print("delta-ten C4 + 6N implicit equation:", curve)
print("delta-ten C4 + 6N simplified complement:", infinite_cyclic)
print(
    "delta-ten C4 + 6N affine singularities: one length-four cusp, "
    "one length-seven contact, and six reduced nodes"
)
