"""Reproduce the delta-ten ordinary-quadruple member and complement."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


Source.<u, v, t> = ZZ[]
source_x = t^2 + t^3 + t^4
source_y = 2*t^9 - t^7 - 4*t^6 - 4*t^5
raw_implicit = (u - source_x).resultant(v - source_y, t)
coefficient_content = gcd(abs(coefficient) for coefficient in raw_implicit.coefficients())
primitive_implicit = Source(raw_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(v^4) < 0:
    primitive_implicit = -primitive_implicit

R.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], R)
regenerated_curve = source_to_plane(primitive_implicit)
curve = (
    -16*X^9 + 248*X^8 - 543*X^7 + 276*X^6*Y + 424*X^6
    - 583*X^5*Y + 116*X^4*Y^2 - 112*X^5 + 395*X^4*Y
    - 177*X^3*Y^2 + 18*X^2*Y^3 - 84*X^3*Y + 76*X^2*Y^2
    - 13*X*Y^3 + Y^4 - 16*X*Y^2 - Y^3 + 7*Y^2
)
assert coefficient_content == 1
assert regenerated_curve == curve
assert curve(X=t^2 + t^3 + t^4, Y=2*t^9 - t^7 - 4*t^6 - 4*t^5) == 0
assert curve.total_degree() == 9
curve_factorization = curve.factor()
assert len(curve_factorization) == 1
assert curve_factorization[0][1] == 1


# Sage 10.8's exact unsimplified affine Zariski--van Kamp presentation.  The
# signed integers name generators 1,...,4, with sign recording inversion.
expected_relations = (
    (-4, -3, -2, -1, 2, 3, 4, -3, 2, 3),
    (2, 1, -2, -1),
    (
        -4, -3, 1, 3, 4, -3, 1, 3, 4, -3, 1, 3, -4, -3, -1, 3,
        -4, -3, -1, 3,
    ),
    (2, 3, 4, 1, -3, -2, -1, -4),
    (1, -2, -1, -4, -3, 2, 3, 4),
    (1, -4, -3, -2, -1, 2, 3, 4),
    (-2, -1, 3, 1, 2, -1, -3, 1),
    (-4, -3, 1, 3, 4, -3, -1, 3),
    (-4, -3, -1, -3, 1, 3, 4, -3, 1, 3),
    (-3, 1, 3, -1),
)

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert len(actual_relations) == 10
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


# The Jacobian scheme consists of the intrinsic length-four T(2,5) cusp,
# the length-nine ordinary quadruple point at (1,-1), and four reduced nodes.
jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 17
assert jacobian.radical().vector_space_dimension() == 6
primary_components = jacobian.primary_decomposition()
assert len(primary_components) == 3
origin_ideal = R.ideal([X, Y])
quadruple_ideal = R.ideal([X - 1, Y + 1])
components = tuple((component, component.radical()) for component in primary_components)
origin_components = tuple(
    component for component, radical in components if radical == origin_ideal
)
quadruple_components = tuple(
    component for component, radical in components if radical == quadruple_ideal
)
node_components = tuple(
    component
    for component, radical in components
    if radical != origin_ideal and radical != quadruple_ideal
)
assert len(origin_components) == 1
assert len(quadruple_components) == 1
assert len(node_components) == 1
origin_component = origin_components[0]
quadruple_component = quadruple_components[0]
node_component = node_components[0]
assert origin_component.vector_space_dimension() == 4
assert origin_component.radical().vector_space_dimension() == 1
assert origin_component != origin_component.radical()
assert quadruple_component.vector_space_dimension() == 9
assert quadruple_component.radical().vector_space_dimension() == 1
assert quadruple_component != quadruple_component.radical()
assert node_component.vector_space_dimension() == 4
assert node_component == node_component.radical()

node_x_polynomial = 64*X^4 - 336*X^3 + 484*X^2 + 189*X - 588
node_x_univariate = node_x_polynomial.univariate_polynomial()
assert node_x_univariate.is_irreducible()
assert node_x_univariate.discriminant() == -265262438203392
assert node_component.reduce(node_x_polynomial) == 0
assert jacobian.radical().reduce(X * (X - 1) * node_x_polynomial) == 0

print("quadruple implicit equation:", curve)
print("quadruple raw presentation:", group)
print("quadruple simplified complement:", infinite_cyclic)
print(
    "quadruple affine singularities: one length-four cusp, one length-nine "
    "ordinary quadruple point, and four reduced nodes"
)
print("delta-ten ordinary-quadruple Sage certificate verified: PASS")
