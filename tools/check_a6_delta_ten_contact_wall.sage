"""Reproduce the delta-ten contact-two wall and its affine complement."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


# Regenerate the primitive implicit equation over the integers.  Keeping this
# calculation separate from the literal below catches transcription errors in
# both the parametrization and the stored certificate.
Source.<u, v, t> = ZZ[]
source_x = t^2 + t^3 + t^4
source_y = 5*t^9 - 16*t^6 - 12*t^5
raw_implicit = (u - source_x).resultant(v - source_y, t)
coefficient_content = gcd(abs(coefficient) for coefficient in raw_implicit.coefficients())
primitive_implicit = Source(raw_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(v^4) < 0:
    primitive_implicit = -primitive_implicit

R.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], R)
regenerated_curve = source_to_plane(primitive_implicit)
curve = (
    -625*X^9 + 12000*X^8 - 42600*X^7 + 4225*X^6*Y + 66496*X^6
    - 14715*X^5*Y + 675*X^4*Y^2 - 19152*X^5 + 20496*X^4*Y
    - 1797*X^3*Y^2 + 45*X^2*Y^3 - 3724*X^3*Y
    + 1410*X^2*Y^2 - 60*X*Y^3 + Y^4 - 48*X*Y^2 + 10*Y^3
    + 133*Y^2
)
assert coefficient_content == 1
assert regenerated_curve == curve

ParameterRing.<z> = QQ[]
assert curve(X=z^2 + z^3 + z^4, Y=5*z^9 - 16*z^6 - 12*z^5) == 0
assert curve(X=-2, Y=-128) == 0


# Sage 10.8's exact, unsimplified affine Zariski--van Kamp presentation.  The
# signed integers refer to generators 1,...,4, with sign recording inversion.
expected_relations = (
    (3, 2, -3, -2),
    (4, 3, -4, -3),
    (
        -4, -3, -2, -1, 2, 3, 4, 3, -4, -3, -2, 1, 2, 3, 4, 3,
        -4, -3, -2, 1, 2, 3, 4, -3, -4, -3, -2, -1, 2, 3, 4, -3,
    ),
    (-4, -3, 2, 3, 4, -3, -2, 3),
    (2, 1, -2, -1),
    (-4, -3, -2, -1, 3, -4, -3, 1, 2, 3, 4, 3),
    (
        -2, -1, 3, 4, -3, 1, 2, -1, 3, 4, -3, 1, 2, -1, 3, 4,
        -3, 1, -2, -1, 3, -4, -3, 1, -2, -1, 3, -4, -3, 1,
    ),
    (-2, -1, 3, -4, -3, 1, 3, 4, -3, 1),
    (-4, -3, -1, 3, 4, 3, -4, -3, 1, 3, 4, -3),
    (-4, -3, 1, 3, 4, -3, -1, 3),
    (-4, -3, -2, 3, 4, 3, -4, -3, 2, 3, 4, -3),
    (-4, -3, 2, 3, 4, -3, -2, 3),
)

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert len(actual_relations) == 12
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


# The Jacobian scheme separates into the intrinsic T(2,5) cusp at the origin,
# the contact-two collision at (-2,-128), and eight reduced ordinary nodes.
jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 15
assert jacobian.radical().vector_space_dimension() == 10

primary_components = jacobian.primary_decomposition()
assert len(primary_components) == 3
origin_ideal = R.ideal([X, Y])
contact_ideal = R.ideal([X + 2, Y + 128])
components = tuple((component, component.radical()) for component in primary_components)
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
assert contact_component.vector_space_dimension() == 3
assert contact_component.radical().vector_space_dimension() == 1
assert node_component.vector_space_dimension() == 8
assert node_component == node_component.radical()

node_x_polynomial = (
    15625*X^8 - 87500*X^7 + 332500*X^6 - 747150*X^5
    + 1299030*X^4 - 493779*X^3 - 58605*X^2 + 1202586*X - 636804
)
assert node_x_polynomial.degree() == 8
node_x_univariate = node_x_polynomial.univariate_polynomial()
assert node_x_univariate.is_irreducible()
assert node_x_univariate.discriminant() != 0
assert node_component.reduce(node_x_polynomial) == 0
assert jacobian.radical().reduce(X * (X + 2) * node_x_polynomial) == 0

print("delta-ten contact-wall implicit equation:", curve)
print("delta-ten contact-wall raw presentation:", group)
print("delta-ten contact-wall simplified complement:", infinite_cyclic)
print(
    "delta-ten contact-wall affine singularities: one length-four cusp, "
    "one length-three contact, and eight reduced nodes"
)
print("delta-ten contact-wall Sage certificate verified: PASS")
