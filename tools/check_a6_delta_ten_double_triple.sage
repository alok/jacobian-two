"""Reproduce the delta-ten two-ordinary-triple member and complement."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


Source.<u, v, t> = ZZ[]
source_x = t^2 + 3*t^3 + t^4
source_y = 15*t^5 + 4*t^6 + 18*t^7 + 15*t^8 + 3*t^9
raw_implicit = (u - source_x).resultant(v - source_y, t)
coefficient_content = gcd(abs(coefficient) for coefficient in raw_implicit.coefficients())
primitive_implicit = Source(raw_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(v^4) < 0:
    primitive_implicit = -primitive_implicit

R.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], R)
regenerated_curve = source_to_plane(primitive_implicit)
curve = (
    -81*X^9 - 27*X^8 + 441*X^7 - 90*X^6*Y + 1495*X^6
    + 402*X^5*Y + 108*X^4*Y^2 + 1125*X^5 - 543*X^4*Y
    - 290*X^3*Y^2 + 21*X^2*Y^3 - 1085*X^3*Y - 240*X^2*Y^2
    - 33*X*Y^3 + Y^4 + 201*X*Y^2 - 40*Y^3 - 5*Y^2
)
assert coefficient_content == 1
assert regenerated_curve == curve
assert curve(X=source_x, Y=source_y) == 0
assert curve.total_degree() == 9
curve_factorization = curve.factor()
assert len(curve_factorization) == 1
assert curve_factorization[0][1] == 1

expected_relations = (
    (4, -3),
    (4, 3, -4, -3),
    (1, 3, 4, -3, -4, -3, -2, -1, 2, 3, 4, 3, -4, -3),
    (4, -3, -4, -3, 1, 2, 3, 4, 3, -4, -3, -2, -1, 3),
    (2, 1, -2, -1),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-2, -1, -2, -1, 3, -4, -3, 1, 2, 1, 2, 1),
    (-2, -1, 3, -4, -3, 1, 2, 1, -2, -1, 3, 4, -3, 1, 2, -1),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-4, 1, 3, 4, -3, -1),
    (-3, 4, 1, 3, -1, -4),
)

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert relations == expected_relations
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
assert jacobian.radical().vector_space_dimension() == 7
primary_components = jacobian.primary_decomposition()
assert len(primary_components) == 4

origin_ideal = R.ideal([X, Y])
first_triple_ideal = R.ideal([X + 1, Y - 1])
second_triple_ideal = R.ideal([X - 5, Y + 125])
components = tuple(
    (component, component.radical()) for component in primary_components
)
origin_components = tuple(
    component for component, radical in components if radical == origin_ideal
)
first_triple_components = tuple(
    component for component, radical in components if radical == first_triple_ideal
)
second_triple_components = tuple(
    component for component, radical in components if radical == second_triple_ideal
)
node_components = tuple(
    component
    for component, radical in components
    if radical not in (origin_ideal, first_triple_ideal, second_triple_ideal)
)
assert len(origin_components) == 1
assert len(first_triple_components) == 1
assert len(second_triple_components) == 1
assert len(node_components) == 1
assert origin_components[0].vector_space_dimension() == 4
assert origin_components[0].radical().vector_space_dimension() == 1
assert first_triple_components[0].vector_space_dimension() == 4
assert first_triple_components[0].radical().vector_space_dimension() == 1
assert second_triple_components[0].vector_space_dimension() == 4
assert second_triple_components[0].radical().vector_space_dimension() == 1
node_component = node_components[0]
assert node_component.vector_space_dimension() == 4
assert node_component == node_component.radical()

node_x_polynomial = 81*X^4 - 81*X^3 - 459*X^2 - 59*X + 6
node_x_univariate = node_x_polynomial.univariate_polynomial()
assert node_x_univariate.is_irreducible()
assert node_x_univariate.discriminant() == 402849850982400
assert node_component.reduce(node_x_polynomial) == 0
assert jacobian.radical().reduce(
    X * (X + 1) * (X - 5) * node_x_polynomial
) == 0

print("delta-ten double-triple implicit equation:", curve)
print("delta-ten double-triple raw presentation:", group)
print("delta-ten double-triple simplified complement:", infinite_cyclic)
print(
    "delta-ten double-triple affine singularities: one length-four cusp, "
    "two length-four ordinary triple points, and four reduced nodes"
)
print("delta-ten double-triple Sage certificate verified: PASS")
