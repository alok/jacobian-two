"""Reproduce the generic delta-ten affine curve and complement certificate."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


R.<X, Y> = QQ[]
curve = (
    -X^9 - 4*X^8 - 8*X^7 + 21*X^6*Y - 4*X^6 + 41*X^5*Y
    - X^5 + 27*X^4*Y^2 + 20*X^4*Y + 22*X^3*Y^2 + 5*X^3*Y
    + 9*X^2*Y^3 + 16*X^2*Y^2 + 5*X*Y^3 + 4*X*Y^2
    + Y^4 - Y^3 + Y^2
)

ParameterRing.<t> = QQ[]
assert curve(X=t^2 + t^3 + t^4, Y=t^5 + t^9) == 0

expected_relations = (
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-3, 1, 3, -1),
    (2, 1, -2, -1),
    (-3, 1, 3, -1),
    (-3, -2, -3, 2, 3, 1),
    (
        -3, -2, -3, -2, -3, -2, 3, 2, 3, 2, 3, -1, -3, -2, -1, 4,
        1, 2, 3, 1, -3, -2, -3, -2, -3, 2, 3, 2, 3, 2, 3, -1, -3,
        -2, -1, -4, 1, 2, 3, 1,
    ),
    (3, 2, -3, -2),
    (-3, -2, -1, -4, 1, 2, 3, 1, -3, -2, -1, 4, 1, 2, 3, -1),
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -1, 4, 1, 2),
    (-4, 1, 4, -1),
    (3, 2, -3, -2),
    (4, 3, -4, -3),
)

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert actual_relations == expected_relations

isomorphism = group.simplification_isomorphism()
infinite_cyclic = isomorphism.codomain()
assert len(infinite_cyclic.gens()) == 1
assert len(infinite_cyclic.relations()) == 0
assert all(
    isomorphism(generator) == infinite_cyclic.gen(0)
    for generator in group.gens()
)
assert isomorphism.section()(infinite_cyclic.gen(0)) == group.gen(0)

jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 14
assert jacobian.radical().vector_space_dimension() == 11

primary_components = jacobian.primary_decomposition()
assert len(primary_components) == 2
origin_ideal = R.ideal([X, Y])
components = tuple((component, component.radical()) for component in primary_components)
origin_components = tuple(
    component for component, radical in components if radical == origin_ideal
)
node_components = tuple(
    component for component, radical in components if radical != origin_ideal
)
assert len(origin_components) == 1
assert len(node_components) == 1
origin_component = origin_components[0]
node_component = node_components[0]
assert origin_component.vector_space_dimension() == 4
assert node_component.vector_space_dimension() == 10
assert node_component == node_component.radical()

node_x_polynomial = (
    X^10 - 10*X^9 + 3*X^8 + 191*X^7 - 226*X^6 - 712*X^5
    - 668*X^4 - 331*X^3 - 95*X^2 - 15*X - 1
)
assert node_component.reduce(node_x_polynomial) == 0
assert jacobian.radical().reduce(X * node_x_polynomial) == 0

print("delta-ten raw presentation:", group)
print("delta-ten simplified complement:", infinite_cyclic)
print("delta-ten affine singularities: one length-four cusp and ten reduced nodes")
print("generic conditional A6 delta-ten Sage certificate verified")
