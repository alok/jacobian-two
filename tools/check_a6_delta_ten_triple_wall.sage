"""Reproduce the delta-ten ordinary-triple wall and affine complement."""

from sage.schemes.curves.zariski_vankampen import fundamental_group


# Regenerate the primitive implicit equation directly from the rational
# parametrization.  This calculation is independent of the literal
# polynomial below and therefore catches transcription and scaling errors.
Source.<u, v, t> = QQ[]
source_x = t^2 + 2*t^3 + t^4
source_y = 114/625*t^5 - 1/5*t^6 - 3*t^7 - 3*t^8 + t^9
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
    -152587890625*X^9
    + 50791992187500*X^8
    - 8773818750000*X^7
    + 99914550781250*X^6*Y
    + 505313640000*X^6
    + 183958593750000*X^5*Y
    + 39825439453125*X^4*Y^2
    - 9701462016*X^5
    - 22186507500000*X^4*Y
    - 284648681640625*X^3*Y^2
    + 4577636718750*X^2*Y^3
    + 648518400000*X^3*Y
    + 106483359375000*X^2*Y^2
    + 14560546875000*X*Y^3
    + 152587890625*Y^4
    - 10259156250000*X*Y^2
    + 421875000000*Y^3
    + 291600000000*Y^2
)
assert coefficient_denominator == 152587890625
assert coefficient_content == 1
assert regenerated_curve == curve
assert curve.total_degree() == 9
curve_factorization = curve.factor()
assert len(curve_factorization) == 1
assert curve_factorization[0][1] == 1

ParameterRing.<z> = QQ[]
parameter_x = z^2 + 2*z^3 + z^4
parameter_y = 114/625*z^5 - 1/5*z^6 - 3*z^7 - 3*z^8 + z^9
assert curve(X=parameter_x, Y=parameter_y) == 0
assert all(
    (parameter_x(parameter), parameter_y(parameter)) == (36/625, 0)
    for parameter in (-3/5, -2/5, 1/5)
)


# Sage 10.8's exact unsimplified affine Zariski--van Kamp presentation.  The
# signed integers name generators 1,...,4, with sign recording inversion.
expected_relations = (
    (4, 3, -4, -3),
    (2, 1, -2, -1),
    (
        -4, -3, -2, -1, 3, -4, -3, 1, 2, 3, 4, 3, -4, -3, -2, -1,
        3, 4, -3, 1, 2, 3, 4, -3,
    ),
    (2, 1, -2, -1),
    (
        -4, -3, -2, -1, 3, -4, -3, -2, -1, 3, -4, -3, 1, 2, 3, 4,
        -3, 1, 2, 3, 4, 3,
    ),
    (2, 1, 2, 1, 2, -1, -2, -1, -2, -1),
    (-3, -2, -1, 3, -4, -3, 1, 2, 3, 4),
    (-2, -1, 2, 3, 4, -3, 1, 3, -4, -3),
    (-2, -1, -2, -1, -2, -1, 3, -4, -3, 1, 2, 1, 2, 1, 2, 1),
    (
        -2, -1, -2, -1, 3, -4, -3, 1, 2, 1, 2, 1, -2, -1, -2, -1,
        3, 4, -3, 1, 2, 1, 2, -1,
    ),
    (-2, -1, 3, -4, -3, 1, 2, 1, -2, -1, 3, 4, -3, 1, 2, -1),
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


# The Jacobian scheme has the intrinsic length-four T(2,5) cusp, a
# length-four ordinary triple point, and two reduced Galois packets containing
# four and three ordinary nodes respectively.
jacobian = R.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 15
assert jacobian.radical().vector_space_dimension() == 9

primary_components = jacobian.primary_decomposition()
assert len(primary_components) == 4
origin_ideal = R.ideal([X, Y])
triple_ideal = R.ideal([X - 36/625, Y])
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
assert triple_component.vector_space_dimension() == 4
assert triple_component.radical().vector_space_dimension() == 1
assert triple_component != triple_component.radical()
assert sorted(component.vector_space_dimension() for component in node_components) == [3, 4]
assert all(component == component.radical() for component in node_components)

cubic_node_x = 625*X^3 - 190861*X^2 + 203652*X - 20736
quartic_node_x = (
    390625*X^4
    + 332908125*X^3
    + 22873309446*X^2
    + 7957013400*X
    + 579156480
)
cubic_node_x_univariate = cubic_node_x.univariate_polynomial()
quartic_node_x_univariate = quartic_node_x.univariate_polynomial()
assert cubic_node_x_univariate.is_irreducible()
assert quartic_node_x_univariate.is_irreducible()
assert cubic_node_x_univariate.discriminant() != 0
assert quartic_node_x_univariate.discriminant() != 0
cubic_components = tuple(
    component for component in node_components if component.reduce(cubic_node_x) == 0
)
quartic_components = tuple(
    component for component in node_components if component.reduce(quartic_node_x) == 0
)
assert len(cubic_components) == 1
assert len(quartic_components) == 1
assert cubic_components[0].vector_space_dimension() == 3
assert quartic_components[0].vector_space_dimension() == 4
assert cubic_components[0] != quartic_components[0]
assert jacobian.radical().reduce(
    X * (625*X - 36) * cubic_node_x * quartic_node_x
) == 0

print("delta-ten triple-wall implicit equation:", curve)
print("delta-ten triple-wall raw presentation:", group)
print("delta-ten triple-wall simplified complement:", infinite_cyclic)
print(
    "delta-ten triple-wall affine singularities: one length-four cusp, "
    "one length-four triple point, four reduced nodes, and three reduced nodes"
)
print("delta-ten triple-wall Sage certificate verified: PASS")
