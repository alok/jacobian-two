"""Independently audit the P-critical triple boundary and two hostile members."""


# The removed fourth-root locus has three labeled components.  Saturation by
# the valid triple-source chart removes zero/repeated sources and the
# sigma_2 chart boundary.
Base.<root_a, root_b, root_c> = QQ[]
sigma_1 = root_a + root_b + root_c
sigma_2 = root_a*root_b + root_a*root_c + root_b*root_c
sigma_3 = root_a*root_b*root_c
base_constraint = sigma_2^2 - sigma_1*sigma_3 - sigma_2
fourth_factors = (
    sigma_2 + root_b*root_c,
    sigma_2 + root_a*root_c,
    sigma_2 + root_a*root_b,
)
chart_localizer = (
    root_a*root_b*root_c*sigma_2
    * (root_a - root_b)*(root_a - root_c)*(root_b - root_c)
)
boundary_ideal = Base.ideal(
    [base_constraint, prod(fourth_factors)]
)
boundary_saturation, boundary_exponent = boundary_ideal.saturation(
    Base.ideal(chart_localizer)
)

critical_a = Base.ideal(
    [
        3*root_a*(root_b + root_c) + 2*root_a^2 - 2,
        3*root_b*root_c - root_a^2 + 1,
    ]
)
critical_b = Base.ideal(
    [
        3*root_b*(root_a + root_c) + 2*root_b^2 - 2,
        3*root_a*root_c - root_b^2 + 1,
    ]
)
critical_c = Base.ideal(
    [
        3*root_c*(root_a + root_b) + 2*root_c^2 - 2,
        3*root_a*root_b - root_c^2 + 1,
    ]
)
critical_components = (critical_a, critical_b, critical_c)
expected_boundary = critical_a.intersection(critical_b).intersection(critical_c)
assert boundary_exponent >= 1
assert boundary_saturation == expected_boundary
assert all(component.dimension() == 1 for component in critical_components)
assert all(component.is_prime() for component in critical_components)


def primitive_implicit(source_x, source_y):
    """Return a primitive integral implicit equation."""

    Source.<u, v, t> = QQ[]
    parameter_to_source = source_x.parent().hom([t], Source)
    source_x_in_ring = parameter_to_source(source_x)
    source_y_in_ring = parameter_to_source(source_y)
    raw = (u - source_x_in_ring).resultant(v - source_y_in_ring, t)
    denominator = lcm(coefficient.denominator() for coefficient in raw.coefficients())
    integral = Source(denominator * raw)
    content = gcd(abs(ZZ(coefficient)) for coefficient in integral.coefficients())
    primitive = Source(integral / content)
    R.<X, Y> = QQ[]
    source_to_plane = Source.hom([X, Y, 0], R)
    return source_to_plane(primitive)


Parameter.<z> = QQ[]

# Hostile T112 fixture: the critical source z=2 is immersed vertically and
# the two roots of z^2+z+1 form the tangent pair.
t112_x = z^4 - 3*z^3 + z^2
t112_y = z^9 - 17/4*z^7 + 11/4*z^6 - 17/4*z^5
t112_curve = primitive_implicit(t112_x, t112_y)
T112Ring = t112_curve.parent()
X, Y = T112Ring.gens()
t112_jacobian = T112Ring.ideal(
    [t112_curve, t112_curve.derivative(X), t112_curve.derivative(Y)]
)
assert t112_jacobian.dimension() == 0
assert t112_jacobian.vector_space_dimension() == 16
assert t112_jacobian.radical().vector_space_dimension() == 8
t112_origin = T112Ring.ideal([X, Y])
t112_triple = T112Ring.ideal([X + 4, Y - 8])
t112_components = tuple(
    (component, component.radical())
    for component in t112_jacobian.primary_decomposition()
)
t112_origin_components = tuple(
    component for component, radical in t112_components if radical == t112_origin
)
t112_triple_components = tuple(
    component for component, radical in t112_components if radical == t112_triple
)
t112_node_components = tuple(
    component
    for component, radical in t112_components
    if radical != t112_origin and radical != t112_triple
)
assert len(t112_origin_components) == 1
assert len(t112_triple_components) == 1
assert sum(component.vector_space_dimension() for component in t112_node_components) == 6
assert t112_origin_components[0].vector_space_dimension() == 4
assert t112_triple_components[0].vector_space_dimension() == 6
assert all(component == component.radical() for component in t112_node_components)

# Hostile C2+T111 fixture: the same critical source is immersed vertically,
# while a separate pair over (2,89/5) has contact order two.
c2_x = z^4 - 3*z^3 + z^2
c2_y = z^9 - 173/10*z^7 + 361/10*z^6 - 173/10*z^5
c2_curve = primitive_implicit(c2_x, c2_y)
C2Ring = c2_curve.parent()
X2, Y2 = C2Ring.gens()
c2_jacobian = C2Ring.ideal(
    [c2_curve, c2_curve.derivative(X2), c2_curve.derivative(Y2)]
)
assert c2_jacobian.dimension() == 0
assert c2_jacobian.vector_space_dimension() == 16
assert c2_jacobian.radical().vector_space_dimension() == 8
c2_origin = C2Ring.ideal([X2, Y2])
c2_contact = C2Ring.ideal([X2 - 2, Y2 - 89/5])
c2_triple = C2Ring.ideal([X2 + 4, Y2 - 272/5])
c2_components = tuple(
    (component, component.radical())
    for component in c2_jacobian.primary_decomposition()
)
c2_origin_components = tuple(
    component for component, radical in c2_components if radical == c2_origin
)
c2_contact_components = tuple(
    component for component, radical in c2_components if radical == c2_contact
)
c2_triple_components = tuple(
    component for component, radical in c2_components if radical == c2_triple
)
c2_node_components = tuple(
    component
    for component, radical in c2_components
    if radical != c2_origin and radical != c2_contact and radical != c2_triple
)
assert len(c2_origin_components) == 1
assert len(c2_contact_components) == 1
assert len(c2_triple_components) == 1
assert sum(component.vector_space_dimension() for component in c2_node_components) == 5
assert c2_origin_components[0].vector_space_dimension() == 4
assert c2_contact_components[0].vector_space_dimension() == 3
assert c2_triple_components[0].vector_space_dimension() == 4
assert all(component == component.radical() for component in c2_node_components)

print("P-critical labeled boundary saturation exponent:", boundary_exponent)
print("P-critical labeled boundary components: three prime curves")
print("P-critical T112 hostile Jacobian lengths: 4 + 6 + 6")
print("P-critical C2+T111 hostile Jacobian lengths: 4 + 3 + 4 + 5")
print("P-critical triple boundary Sage certificate verified: PASS")
