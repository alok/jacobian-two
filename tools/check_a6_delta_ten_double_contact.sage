"""Reproduce the dominant delta-ten ``C2^2 + 6N`` member in Sage."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


Source.<u, v, t> = ZZ[]
source_x = t^2 + t^3 + t^4
source_y = 6*t^5 + 3*t^6 + 6*t^7 + 2*t^9
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
curve = (
    -16*X^9 - 72*X^8 - 1176*X^7 - 144*X^6*Y + 387*X^6
    + 882*X^5*Y + 60*X^4*Y^2 - 36*X^5 - 273*X^4*Y
    - 84*X^3*Y^2 + 18*X^2*Y^3 + 24*X^3*Y + 54*X^2*Y^2
    - 12*X*Y^3 + Y^4 - 12*X*Y^2 + 2*Y^3 + Y^2
)
assert coefficient_content == 1
assert regenerated_curve == curve

ParameterRing.<z> = QQ[]
assert curve(X=z^2 + z^3 + z^4, Y=6*z^5 + 3*z^6 + 6*z^7 + 2*z^9) == 0
assert curve(X=0, Y=-1) == 0
assert curve(X=-2, Y=7) == 0

# Sage 10.8's exact unsimplified affine Zariski--van Kamp presentation.
expected_relations = (
    (4, 3, -4, -3),
    (2, 1, 2, 1, -2, -1, -2, -1),
    (-2, -1, 3, 4, -3, -4, -3, 1, 2, 1),
    (
        -2, -1, 3, 4, -3, -4, -3, 1, -2, -1, 3, 4, 3, -4, -3, 1,
        2, -1, -2, -1, 3, 4, -3, 1, 2, 1,
    ),
    (
        -2, -1, 3, -4, -3, 1, 2, 1, -2, -1, 3, 4, -3, 1, 2, 1,
        -2, -1, 3, 4, -3, 1, 2, -1, -2, -1, 3, -4, -3, 1, 2, -1,
    ),
    (
        -2, -1, 3, 4, 3, -4, -3, 1, 2, -1, 3, 4, 3, -4, -3, 1,
        2, -1, 3, 4, 3, -4, -3, 1, -2, -1, 3, 4, -3, -4, -3, 1,
        -2, -1, 3, 4, -3, -4, -3, 1,
    ),
    (
        -2, -1, 3, 4, -3, -4, -3, 1, -2, -1, 3, 4, -3, 1, 2, 1,
        -2, -1, 3, -4, -3, 1, 2, -1, 3, 4, 3, -4, -3, 1, 2, -1,
        3, 4, -3, -4, -3, 1, -2, -1, 3, 4, -3, 1, 2, -1, -2, -1,
        3, -4, -3, 1, 2, -1, 3, 4, 3, -4, -3, 1,
    ),
    (
        -2, -1, 3, -4, -3, 1, 2, -1, 3, 4, -3, -4, -3, 1, -2, -1,
        3, 4, -3, 1, 2, 1, -2, -1, 3, -4, -3, 1, 2, -1, 3, 4,
        3, -4, -3, 1, -2, -1, 3, 4, -3, 1, 2, -1,
    ),
    (
        -2, -1, 3, -4, -3, 1, 2, -1, 3, 4, 3, -4, -3, 1, -2, -1,
        3, 4, -3, 1, 2, -1, 3, 4, -3, -4, -3, 1,
    ),
    (-2, -1, 3, 4, -3, 1, 2, -1, 3, -4, -3, 1),
    (-4, -3, 1, 3, 4, -3, -1, 3),
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

jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 16
assert jacobian.radical().vector_space_dimension() == 9
components = tuple(
    (component, component.radical())
    for component in jacobian.primary_decomposition()
)
assert len(components) == 4
origin_ideal = R.ideal([X, Y])
first_contact_ideal = R.ideal([X, Y + 1])
second_contact_ideal = R.ideal([X + 2, Y - 7])

def components_with_radical(expected):
    return tuple(
        component for component, radical in components if radical == expected
    )

origin_components = components_with_radical(origin_ideal)
first_contact_components = components_with_radical(first_contact_ideal)
second_contact_components = components_with_radical(second_contact_ideal)
special_radicals = (origin_ideal, first_contact_ideal, second_contact_ideal)
node_components = tuple(
    component for component, radical in components
    if radical not in special_radicals
)
assert len(origin_components) == 1
assert len(first_contact_components) == 1
assert len(second_contact_components) == 1
assert len(node_components) == 1
origin_component = origin_components[0]
first_contact_component = first_contact_components[0]
second_contact_component = second_contact_components[0]
node_component = node_components[0]
assert origin_component.vector_space_dimension() == 4
assert origin_component.radical().vector_space_dimension() == 1
assert first_contact_component.vector_space_dimension() == 3
assert first_contact_component.radical().vector_space_dimension() == 1
assert second_contact_component.vector_space_dimension() == 3
assert second_contact_component.radical().vector_space_dimension() == 1
assert node_component.vector_space_dimension() == 6
assert node_component == node_component.radical()

node_x_polynomial = (
    64*X^6 + 160*X^5 - 22352*X^4 + 14232*X^3
    - 3540*X^2 + 405*X - 18
)
node_x_univariate = node_x_polynomial.univariate_polynomial()
assert node_x_univariate.is_irreducible()
assert node_x_univariate.discriminant() != 0
assert node_component.reduce(node_x_polynomial) == 0
assert jacobian.radical().reduce(X * (X + 2) * node_x_polynomial) == 0

print("delta-ten C2^2 + 6N implicit equation:", curve)
print("delta-ten C2^2 + 6N raw presentation:", group)
print("delta-ten C2^2 + 6N simplified complement:", simplified)
print(
    "delta-ten C2^2 + 6N affine singularities: one length-four cusp, "
    "two length-three contacts, and six reduced nodes"
)
print("delta-ten C2^2 + 6N Sage certificate verified: PASS")
