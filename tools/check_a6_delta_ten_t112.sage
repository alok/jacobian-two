"""Reproduce one conditional delta-ten T112 wall member."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


# Global labeled incidence component.  The valid-base product removes zero or
# repeated sources, the sigma_2 chart boundary, and coincidence with the
# fourth root.  On this localization the triple-root surface is smooth and
# every incidence matrix has rank three.
Base.<base_a, base_b, base_c> = QQ[]
base_sigma_1 = base_a + base_b + base_c
base_sigma_2 = base_a*base_b + base_a*base_c + base_b*base_c
base_sigma_3 = base_a*base_b*base_c
base_constraint = (
    base_sigma_2^2 - base_sigma_1*base_sigma_3 - base_sigma_2
)
base_factorization = base_constraint.factor()
assert len(base_factorization) == 1
assert base_factorization[0][1] == 1
BaseAlgebraic.<alg_a, alg_b, alg_c> = QQbar[]
absolute_factorization = BaseAlgebraic(base_constraint).factor()
assert len(absolute_factorization) == 1
assert absolute_factorization[0][1] == 1

valid_base_product = (
    base_a*base_b*base_c*base_sigma_2
    * (base_a - base_b)*(base_a - base_c)*(base_b - base_c)
    * (base_sigma_2 + base_b*base_c)
    * (base_sigma_2 + base_a*base_c)
    * (base_sigma_2 + base_a*base_b)
)
base_singular_ideal = Base.ideal(
    [
        base_constraint,
        base_constraint.derivative(base_a),
        base_constraint.derivative(base_b),
        base_constraint.derivative(base_c),
    ]
)
assert base_singular_ideal.radical() == Base.ideal([base_a, base_b, base_c])
base_smooth_saturation, base_smooth_exponent = base_singular_ideal.saturation(
    Base.ideal(valid_base_product)
)
assert base_smooth_exponent == 1
assert base_smooth_saturation.is_one()
assert base_singular_ideal.reduce(valid_base_product) == 0


def divided_power(degree, left, right):
    """Return (left^degree-right^degree)/(left-right)."""

    return sum(left^(degree - 1 - index)*right^index for index in range(degree))


row_ab = [
    divided_power(degree, base_a, base_b) for degree in range(5, 9)
]
row_ac = [
    divided_power(degree, base_a, base_c) for degree in range(5, 9)
]
kappa_numerator = base_sigma_3 - base_sigma_1*base_sigma_2


def p_derivative_numerator(parameter):
    """Return sigma_2 times P'(parameter) on the root chart."""

    return (
        base_sigma_2*(2*parameter + 4*parameter^3)
        + 3*kappa_numerator*parameter^2
    )


p_derivative_a = p_derivative_numerator(base_a)
p_derivative_b = p_derivative_numerator(base_b)
row_tangent = [
    degree
    * (
        base_a^(degree - 1)*p_derivative_b
        - base_b^(degree - 1)*p_derivative_a
    )
    for degree in range(5, 9)
]
incidence_matrix = matrix(Base, [row_ab, row_ac, row_tangent])
rank_drop_ideal = Base.ideal(
    [base_constraint]
    + [
        incidence_matrix.matrix_from_columns(columns).determinant()
        for columns in Subsets(range(4), 3)
    ]
)
assert rank_drop_ideal.dimension() == 1
rank_drop_saturation, rank_drop_exponent = rank_drop_ideal.saturation(
    Base.ideal(valid_base_product)
)
assert rank_drop_exponent == 3
assert rank_drop_saturation.is_one()
assert rank_drop_ideal.reduce(valid_base_product^3) == 0

# The exact clean member has only the two orientations of its tangent pair in
# the labeled incidence fiber.  This proves the coefficient projection is
# generically finite; on the clean locus it is the expected degree-two label
# cover.
sample_coefficients = (2089/11875, 3542/2375, 1788/475, 338/95)
sample_equations = [
    sum(row_ab[index]*sample_coefficients[index] for index in range(4))
    + divided_power(9, base_a, base_b),
    sum(row_ac[index]*sample_coefficients[index] for index in range(4))
    + divided_power(9, base_a, base_c),
    sum(row_tangent[index]*sample_coefficients[index] for index in range(4))
    + 9
    * (
        base_a^8*p_derivative_b
        - base_b^8*p_derivative_a
    ),
    kappa_numerator - 2*base_sigma_2,
]
sample_fiber_ideal = Base.ideal([base_constraint] + sample_equations)
sample_fiber_saturation, sample_fiber_exponent = sample_fiber_ideal.saturation(
    Base.ideal(valid_base_product)
)
expected_sample_fiber = Base.ideal(
    [
        base_b^2 + base_b + 6/25,
        base_a + base_b + 1,
        base_c - 1/5,
    ]
)
assert sample_fiber_exponent == 1
assert sample_fiber_saturation == expected_sample_fiber
assert sample_fiber_saturation.dimension() == 0
assert sample_fiber_saturation.vector_space_dimension() == 2
assert sample_fiber_saturation == sample_fiber_saturation.radical()


Source.<u, v, t> = QQ[]
source_x = t^2 + 2*t^3 + t^4
source_y = (
    2089/11875*t^5
    + 3542/2375*t^6
    + 1788/475*t^7
    + 338/95*t^8
    + t^9
)
raw_implicit = (u - source_x).resultant(v - source_y, t)
coefficient_denominator = lcm(
    coefficient.denominator() for coefficient in raw_implicit.coefficients()
)
integral_implicit = Source(coefficient_denominator * raw_implicit)
coefficient_content = gcd(
    abs(ZZ(coefficient)) for coefficient in integral_implicit.coefficients()
)
primitive_implicit = Source(integral_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(v^4) < 0:
    primitive_implicit = -primitive_implicit

R.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], R)
regenerated_curve = source_to_plane(primitive_implicit)
curve = (
    -19885406494140625*X^9
    + 98693762695312500*X^8
    - 12903111675000000*X^7
    - 23744619140625000*X^6*Y
    + 550166298120000*X^6
    - 100341934312500000*X^5*Y
    - 7329711534336*X^5
    + 133129766845703125*X^4*Y^2
    + 10731317108400000*X^4*Y
    - 540905323007812500*X^3*Y^2
    - 289807243200000*X^3*Y
    + 74936584472656250*X^2*Y^3
    + 142878836250000000*X^2*Y^2
    - 159150234375000000*X*Y^3
    - 10428071625000000*X*Y^2
    + 19885406494140625*Y^4
    - 4340460937500000*Y^3
    + 236852100000000*Y^2
)
assert coefficient_denominator == 19885406494140625
assert coefficient_content == 1
assert regenerated_curve == curve
factorization = curve.factor()
assert curve.total_degree() == 9
assert len(factorization) == 1
assert factorization[0][1] == 1

ParameterRing.<z> = QQ[]
parameter_x = z^2 + 2*z^3 + z^4
parameter_y = (
    2089/11875*z^5
    + 3542/2375*z^6
    + 1788/475*z^7
    + 338/95*z^8
    + z^9
)
assert curve(X=parameter_x, Y=parameter_y) == 0
triple_point = (36/625, 7776/37109375)
assert all(
    (parameter_x(parameter), parameter_y(parameter)) == triple_point
    for parameter in (-3/5, -2/5, 1/5)
)

# The Jacobian scheme consists of the intrinsic length-four cusp, a
# length-six T112 point, and reduced packets of two and four ordinary nodes.
jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 16
assert jacobian.radical().vector_space_dimension() == 8
primary_components = jacobian.primary_decomposition()
assert len(primary_components) == 4
origin_ideal = R.ideal([X, Y])
triple_ideal = R.ideal([X - triple_point[0], Y - triple_point[1]])
components = tuple((component, component.radical()) for component in primary_components)
origin_components = tuple(
    component for component, radical in components if radical == origin_ideal
)
triple_components = tuple(
    component for component, radical in components if radical == triple_ideal
)
node_components = tuple(
    component
    for component, radical in components
    if radical != origin_ideal and radical != triple_ideal
)
assert len(origin_components) == 1
assert len(triple_components) == 1
assert len(node_components) == 2
origin_component = origin_components[0]
triple_component = triple_components[0]
assert origin_component.vector_space_dimension() == 4
assert origin_component.radical().vector_space_dimension() == 1
assert origin_component != origin_component.radical()
assert triple_component.vector_space_dimension() == 6
assert triple_component.radical().vector_space_dimension() == 1
assert triple_component != triple_component.radical()
assert sorted(component.vector_space_dimension() for component in node_components) == [2, 4]
assert all(component == component.radical() for component in node_components)

vertical_node_x = 361*X^2 - 1692*X + 1296
graph_node_x = (
    50906640625*X^4
    - 462374666250*X^3
    + 2074717406016*X^2
    + 599188777920*X
    - 38011109760
)
vertical_node_x_univariate = vertical_node_x.univariate_polynomial()
graph_node_x_univariate = graph_node_x.univariate_polynomial()
assert vertical_node_x_univariate.is_irreducible()
assert graph_node_x_univariate.is_irreducible()
assert vertical_node_x_univariate.discriminant() != 0
assert graph_node_x_univariate.discriminant() != 0
vertical_components = tuple(
    component for component in node_components if component.reduce(vertical_node_x) == 0
)
graph_components = tuple(
    component for component in node_components if component.reduce(graph_node_x) == 0
)
assert len(vertical_components) == 1
assert len(graph_components) == 1
assert vertical_components[0].vector_space_dimension() == 2
assert graph_components[0].vector_space_dimension() == 4
assert vertical_components[0] != graph_components[0]
assert jacobian.radical().reduce(
    X * (625*X - 36) * vertical_node_x * graph_node_x
) == 0


# Sage 10.8's exact unsimplified affine Zariski--van Kamp presentation.
expected_relations = (
    (-3, -2, 1, 2, 3, -2, -1, 2),
    (
        -3, -2, -1, -2, 1, 2, 3, -2, 1, 2, -3, -2, -1, 2, 1, 2, 3,
        -2, -1, 2,
    ),
    (-4, 2, 4, 2, 4, 2, -4, -2, -4, -2),
    (-3, -2, -1, -4, -2, 1, 2, 4, 1, 2),
    (
        -4, -2, -4, -2, -4, -2, -4, -2, 1, 2, 4, 2, 4, 2, 4, 2, 4,
        -2, -4, -2, -4, -2, -4, -2, -1, 2, 4, 2, 4, 2, 4, 2,
    ),
    (-4, -2, -4, -2, -4, -2, 1, 2, 4, 2, 4, 2),
    (
        2, 4, 2, -4, -2, -4, -2, -4, -2, -1, 2, 4, 2, 4, -2, -4,
        -2, -4, -2, 1, 2, 4, 2, 4,
    ),
    (
        -4, -2, -1, -4, -2, -1, 2, 4, 1, 2, 4, 1, 2, 4, 2, 4, 2, 4,
        -2, -4, -4, -2, -4, -2,
    ),
    (
        -4, -2, -4, -2, -1, 2, 4, 2, 4, 2, -4, -2, -4, -2, 1, 2, 4,
        2, 4, -2,
    ),
    (-4, -3, 4, 2, -4, 3, 4, -2),
    (4, 3, -4, -3),
)

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert len(actual_relations) == 11
assert actual_relations == expected_relations
isomorphism = group.simplification_isomorphism()
simplified = isomorphism.codomain()
assert len(simplified.gens()) == 1
assert len(simplified.relations()) == 0
assert all(
    isomorphism(generator) == simplified.gen(0)
    for generator in group.gens()
)
section = isomorphism.section()
assert section(simplified.gen(0)) == group.gen(0)
assert isomorphism(section(simplified.gen(0))) == simplified.gen(0)

print(
    "delta-ten T112 global incidence:",
    "absolutely irreducible smooth valid base; rank-drop saturation = (1), exponent 3; "
    "sample projection fiber reduced of length 2",
)
print("delta-ten T112 implicit equation:", curve)
print("delta-ten T112 raw presentation:", group)
print("delta-ten T112 simplified complement:", simplified)
print(
    "delta-ten T112 affine singularities: one length-four cusp, one "
    "length-six T112 point, two reduced nodes, and four reduced nodes"
)
print("delta-ten T112 Sage certificate verified: PASS")
