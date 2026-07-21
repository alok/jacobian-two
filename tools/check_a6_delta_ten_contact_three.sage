"""Reproduce the conditional delta-ten ``C3 + 7N`` member in Sage."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


Source.<u, v, t> = ZZ[]
source_x = t^2 + t^3 + t^4
source_y = 48*t^5 + 8*t^7 - 21*t^8 + 16*t^9
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
    -65536*X^9 - 718463*X^8 - 8446608*X^7 + 264164*X^6*Y
    - 3767040*X^6 + 917752*X^5*Y - 628992*X^5 + 15430*X^4*Y^2
    + 397314*X^4*Y + 11048*X^3*Y^2 + 65520*X^3*Y + 228*X^2*Y^3
    + 5288*X^2*Y^2 + 16*X*Y^3 + 1544*X*Y^2 + Y^4 + 3*Y^3
    + 273*Y^2
)
assert coefficient_content == 1
assert regenerated_curve == curve

ParameterRing.<z> = QQ[]
assert curve(X=z^2 + z^3 + z^4, Y=48*z^5 + 8*z^7 - 21*z^8 + 16*z^9) == 0
assert curve(X=-2, Y=-464) == 0

expected_relations = (
    (2, 1, -2, -1),
    (4, 3, 4, 3, 4, 3, -4, -3, -4, -3, -4, -3),
    (-4, -3, -4, -3, -4, -3, 2, 3, 4, 3, 4, 3),
    (
        -4, -3, -4, -3, -1, 3, 4, 3, 4, 3, -4, -3, -4, -3, 1, 3,
        4, 3, 4, -3,
    ),
    (
        -4, -3, -4, -3, -2, 3, 4, 3, 4, 3, -4, -3, -4, -3, 2, 3,
        4, 3, 4, -3,
    ),
    (2, 1, -2, -1),
    (
        -4, -3, -4, -3, 1, 2, 1, -2, -1, 3, 4, 3, 4, -3, -4, -3,
        1, 2, -1, -2, -1, 3, 4, 3,
    ),
    (
        -4, -3, -4, -3, 1, 2, -1, 3, 4, 3, 4, -3, -4, -3, 1, -2,
        -1, 3, 4, 3,
    ),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, 3, 4, -3, -4, -3, 1, 2, 1),
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
assert jacobian.vector_space_dimension() == 16
assert jacobian.radical().vector_space_dimension() == 9

origin_ideal = R.ideal([X, Y])
contact_ideal = R.ideal([X + 2, Y + 464])
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
assert contact_component.vector_space_dimension() == 5
assert contact_component.radical().vector_space_dimension() == 1
assert node_component.vector_space_dimension() == 7
assert node_component == node_component.radical()

NodeRing.<x> = QQ[]
node_x_polynomial = (
    16777216*x^7 - 76283904*x^6 - 5288038400*x^5 - 6332232000*x^4
    - 3246794416*x^3 - 892022247*x^2 - 131282424*x - 8496306
)
assert node_x_polynomial.is_irreducible()
assert node_x_polynomial.discriminant() != 0
node_to_plane = NodeRing.hom([X], R)
assert node_component.reduce(node_to_plane(node_x_polynomial)) == 0

print("delta-ten C3 + 7N implicit equation:", curve)
print("delta-ten C3 + 7N raw presentation:", group)
print("delta-ten C3 + 7N simplified complement:", infinite_cyclic)
print(
    "delta-ten C3 + 7N affine singularities: one length-four cusp, "
    "one length-five contact, and seven reduced nodes"
)
print("delta-ten C3 + 7N Sage certificate verified: PASS")
