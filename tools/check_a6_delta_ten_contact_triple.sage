"""Reproduce the conditional delta-ten ``C2 + T111 + 5N`` member in Sage."""

from pathlib import Path
import sys

from sage.schemes.curves.zariski_vankampen import fundamental_group

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.a6_delta_ten_contact_triple import (
    CONTACT_TRIPLE_RELATIONS as PYTHON_CONTACT_TRIPLE_RELATIONS,
)


Source.<u, v, t> = ZZ[]
source_x = 2*t^4 + 3*t^3 + 2*t^2
source_y = t^9 - 2*t^7 + t^6 - 2*t^5
raw_implicit = (u - source_x).resultant(v - source_y, t)
coefficient_content = gcd(
    abs(coefficient) for coefficient in raw_implicit.coefficients()
)
primitive_implicit = Source(raw_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(u^9) < 0:
    primitive_implicit = -primitive_implicit

R.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], R)
regenerated_curve = source_to_plane(primitive_implicit)
curve = (
    X^9 - 31*X^8 + 139*X^7 - 529*X^6*Y + 499*X^6
    + 4234*X^5*Y + 1376*X^5 - 1840*X^4*Y^2 + 4183*X^4*Y
    + 19799*X^3*Y^2 + 23392*X^3*Y - 1728*X^2*Y^3
    - 75323*X^2*Y^2 + 10344*X*Y^3 + 44168*X*Y^2
    - 512*Y^4 - 2503*Y^3 - 11008*Y^2
)
assert coefficient_content == 1
assert regenerated_curve == curve

ParameterRing.<z> = QQ[]
assert curve(X=2*z^4 + 3*z^3 + 2*z^2, Y=z^9 - 2*z^7 + z^6 - 2*z^5) == 0
assert curve(X=-5, Y=-2) == 0
assert curve(X=1, Y=4) == 0

expected_relations = (
    (-4, 2, 4, 2, 4, -2, -4, -2),
    (2, -1),
    (-4, 2, 4, 2, 4, 2, -4, -2, -4, -2),
    (-4, 2, 4, -2),
    (
        -4, -2, -4, -2, -4, -2, -1, 2, 4, 2, 4, 2, -4, -2, -4,
        -2, -4, -2, 1, 2, 4, 2, 4, 2, 4, 2, -4, -2, -4, -2, -4,
        -2, -1, 2, 4, 2, 4, 2, 4, -2, -4, -2, -4, -2, 1, 2, 4,
        2, 4, 2, 4, -2,
    ),
    (
        1, 2, 4, 2, 4, 2, -4, -2, -4, -2, -4, -2, -1, -4, -3, 4,
        2, 4, 2, 4, 2, 4, -2, -4, -2, -4, -2, -4, 3, 4,
    ),
    (
        4, 2, 4, -2, -4, -2, -4, -2, -4, 3, 4, 1, -4, -3, 4, 2,
        4, 2, 4, 2, -4, -2, -4, -2, -4, -2, -1, 2, 4, 2,
    ),
    (-4, -2, -4, -2, -4, 3, 4, 2, 4, 2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
    (-3, -2, -1, 2, 4, -2, 1, 2, 3, -2, -1, 2, -4, -2, 1, 2),
    (-4, -2, 1, 2, 4, -2, -1, 2),
)
assert PYTHON_CONTACT_TRIPLE_RELATIONS == expected_relations

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
simplified_group = isomorphism.codomain()
assert len(simplified_group.gens()) == 1
assert len(simplified_group.relations()) == 0
assert all(
    isomorphism(generator) == simplified_group.gen(0)
    for generator in group.gens()
)
section = isomorphism.section()
assert section(simplified_group.gen(0)) == group.gen(0)
assert isomorphism(section(simplified_group.gen(0))) == simplified_group.gen(0)

jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 16
assert jacobian.radical().vector_space_dimension() == 8

origin_ideal = R.ideal([X, Y])
contact_ideal = R.ideal([X + 5, Y + 2])
triple_ideal = R.ideal([X - 1, Y - 4])
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
triple_components = tuple(
    component for component, radical in components if radical == triple_ideal
)
node_components = tuple(
    component
    for component, radical in components
    if radical != origin_ideal
    and radical != contact_ideal
    and radical != triple_ideal
)
assert len(origin_components) == 1
assert len(contact_components) == 1
assert len(triple_components) == 1
assert len(node_components) == 1
origin_component = origin_components[0]
contact_component = contact_components[0]
triple_component = triple_components[0]
node_component = node_components[0]
assert origin_component.vector_space_dimension() == 4
assert origin_component.radical().vector_space_dimension() == 1
assert contact_component.vector_space_dimension() == 3
assert contact_component.radical().vector_space_dimension() == 1
assert triple_component.vector_space_dimension() == 4
assert triple_component.radical().vector_space_dimension() == 1
assert node_component.vector_space_dimension() == 5
assert node_component == node_component.radical()

NodeRing.<x> = QQ[]
node_x_polynomial = (
    512*x^5 - 62528*x^4 + 2718048*x^3 - 48220459*x^2
    + 255118656*x + 115495936
)
assert node_x_polynomial.is_irreducible()
assert node_x_polynomial.discriminant() == (
    -17004216373814107665135371132504044296002359787520
)
node_to_plane = NodeRing.hom([X], R)
assert node_component.reduce(node_to_plane(node_x_polynomial)) == 0

print("delta-ten C2 + T111 + 5N implicit equation:", curve)
print("delta-ten C2 + T111 + 5N raw presentation:", group)
print("delta-ten C2 + T111 + 5N simplified complement:", simplified_group)
print(
    "delta-ten C2 + T111 + 5N affine singularities: one length-four cusp, "
    "one length-three contact, one length-four triple point, and five nodes"
)
print("delta-ten C2 + T111 + 5N Sage certificate verified: PASS")
