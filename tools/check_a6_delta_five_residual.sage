"""Reproduce the residual A6 delta-five primary and complement certificates."""

import sys
from pathlib import Path

from sage.schemes.curves.zariski_vankampen import fundamental_group

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.a6_delta_five_residual import (
    A_RELATIONS,
    B_RELATIONS,
    E_RELATIONS,
    N_RELATIONS,
    P32_MINUS_RELATIONS,
    P32_PLUS_RELATIONS,
    P4_RELATIONS,
    PNT_RELATIONS,
)


ParameterRing.<alpha, beta, gamma> = QQ[]
CollisionRing.<s> = ParameterRing[]
H = (
    s^5 + (gamma+2)*s^4 + (4*gamma-2)*s^3
    + (-alpha+2*beta+3*gamma-4)*s^2
    + (-alpha+3*beta-gamma)*s + alpha
)
C = alpha-beta+gamma-1
L = 135*alpha-108*beta+84*gamma-64
T = (
    alpha*gamma-alpha-2*beta*gamma+3*beta
    + 2*gamma^2-5*gamma+3
)
G = ParameterRing(H.discriminant()/C)

singular_ideal = ParameterRing.ideal(
    [G, G.derivative(alpha), G.derivative(beta), G.derivative(gamma)]
).radical()
singular_components = singular_ideal.primary_decomposition()
assert len(singular_components) == 3
invalid_cusp_line = ParameterRing.ideal(
    [beta-2*gamma+3, alpha-gamma+2]
)
assert any(component == invalid_cusp_line for component in singular_components)
assert invalid_cusp_line.reduce(C) == 0

triple_tangency_components = ParameterRing.ideal([G, T]).radical().primary_decomposition()
assert len(triple_tangency_components) == 3
invalid_critical_line = ParameterRing.ideal(
    [27*beta-31*gamma+31, 27*alpha-8*gamma+12]
)
ordinary_triple_tangency_line = ParameterRing.ideal(
    [gamma-6, 5*alpha-9*beta+45]
)
assert any(
    component == invalid_critical_line
    for component in triple_tangency_components
)
assert any(
    component == ordinary_triple_tangency_line
    for component in triple_tangency_components
)
assert invalid_critical_line.reduce(L) == 0


def matching_components(components, values):
    substitution = {alpha: values[0], beta: values[1], gamma: values[2]}
    return [
        component
        for component in components
        if all(generator.subs(substitution) == 0 for generator in component.gens())
    ]


# Each nonconstant rational map below lands on exactly one one-dimensional
# primary component.  Hence its irreducible image is dense in that component.
ParameterField.<parameter> = FunctionField(QQ)
a_values = (
    -2*parameter^3*(3*parameter^2+12*parameter+13)/(3*parameter+11),
    -(3*parameter^5+9*parameter^4-4*parameter^3-24*parameter^2
      +6*parameter-4)/(3*parameter+11),
    -6*(parameter^2+3*parameter-2)/(3*parameter+11),
)
b_values = (
    2*(parameter-2)*(2*parameter^3-6*parameter^2+3*parameter+3)^2
      /(parameter-3)^3,
    (4*parameter^7-28*parameter^6+56*parameter^5+36*parameter^4
      -259*parameter^3+247*parameter^2+75*parameter-155)
      /(parameter-3)^3,
    -2*(2*parameter^2-8*parameter+9)/(parameter-3),
)
n_values = (
    (parameter+1)^2*(parameter+2)^2*(2*parameter+1)/parameter,
    (parameter^2+3*parameter+1)
      *(parameter^3+4*parameter^2+5*parameter+3)/parameter,
    (parameter^2+4*parameter+1)/parameter,
)
e_values = (9*parameter/5-9, parameter, 6)

a_matches = matching_components(singular_components, a_values)
b_matches = matching_components(singular_components, b_values)
n_matches = matching_components(triple_tangency_components, n_values)
e_matches = matching_components(triple_tangency_components, e_values)
assert all(len(matches) == 1 for matches in (a_matches, b_matches, n_matches, e_matches))
assert a_matches[0] != b_matches[0]
assert n_matches[0] != e_matches[0]
assert all(
    component.dimension() == 1
    for component in (a_matches[0], b_matches[0], n_matches[0], e_matches[0])
)
assert all(any(value.derivative() != 0 for value in values) for values in (
    a_values, b_values, n_values, e_values,
))


R.<X, Y> = QQ[]


def universal_implicit(ring, a, b, g):
    x, y = ring.gens()
    c = a-b+g-1
    return (
        -x^8 + (-3*a+3*b*g-2*b-g^3+g^2)*x^7
        + (-3*a^2+3*a*b*g-a*b-2*a*g^2+2*a*g-b^3+b^2*g-b^2)*x^6
        + (8-3*g)*x^5*y - a^2*c*x^5
        + (-3*a*g+13*a+3*b^2-8*b*g+b+7*g^2-9*g+2)*x^4*y
        + (5*a-2*b)*c*x^3*y + (-3*b+7*g-12)*x^2*y^2
        + (-5*a+6*b-7*g+8)*x*y^2 + y^3 + c*y^2
    )


rational_curves_and_relations = (
    (
        "A",
        (QQ(16)/5, QQ(32)/5, QQ(24)/5),
        5,
        # A: Q=16*t^5+32*t^6+24*t^7+5*t^8.
        -125*X^8 - 2224*X^7 + 2560*X^6 - 160*X^5*Y - 768*X^5
        - 22*X^4*Y + 48*X^3*Y + 12*X^2*Y^2 - 16*X*Y^2
        + Y^3 + 3*Y^2,
        A_RELATIONS,
    ),
    (
        "B",
        (QQ(1)/6, QQ(127)/108, QQ(2)),
        108,
        # B: Q=18*t^5+127*t^6+216*t^7+108*t^8.
        -1259712*X^8 + 256608*X^7 - 16795*X^6 + 23328*X^5*Y
        + 324*X^5 - 3777*X^4*Y + 164*X^3*Y - 165*X^2*Y^2
        + 24*X*Y^2 + Y^3 - Y^2,
        B_RELATIONS,
    ),
    (
        "N",
        (QQ(20)/3, QQ(1), QQ(2)/3),
        3,
        # N: Q=20*t^5+3*t^6+2*t^7+3*t^8.
        -27*X^8 - 536*X^7 - 3376*X^6 + 54*X^5*Y - 6400*X^5
        + 640*X^4*Y + 1504*X^3*Y - 31*X^2*Y^2 - 72*X*Y^2
        + Y^3 + 16*Y^2,
        N_RELATIONS,
    ),
    (
        "E",
        (QQ(-9), QQ(0), QQ(6)),
        1,
        # E: Q=-9*t^5+6*t^7+t^8.
        -X^8 - 153*X^7 + 297*X^6 - 10*X^5*Y + 324*X^5
        + 245*X^4*Y + 180*X^3*Y + 30*X^2*Y^2 + 11*X*Y^2
        + Y^3 - 4*Y^2,
        E_RELATIONS,
    ),
    (
        "P4",
        (QQ(-114244)/81, QQ(-63407)/81, QQ(34)/3),
        81,
        # P4: Q=-114244*t^5-63407*t^6+918*t^7+81*t^8.
        -531441*X^8 - 11769121800*X^7 + 274656502790000*X^6
        - 170586*X^5*Y + 652584576800000*X^5 + 12721440000*X^4*Y
        + 22220300000*X^3*Y + 195675*X^2*Y^2 + 185000*X*Y^2
        + Y^3 - 50000*Y^2,
        P4_RELATIONS,
    ),
    (
        "PNT",
        (QQ(108), QQ(65), QQ(6)),
        1,
        # PNT: Q=108*t^5+65*t^6+6*t^7+t^8.
        -X^8 + 536*X^7 - 175632*X^6 - 10*X^5*Y - 559872*X^5
        + 9280*X^4*Y + 19680*X^3*Y - 165*X^2*Y^2 - 184*X*Y^2
        + Y^3 + 48*Y^2,
        PNT_RELATIONS,
    ),
)

for name, parameters, y_scale, curve, _ in rational_curves_and_relations:
    normalized_curve = universal_implicit(R, *parameters)
    scaled_curve = R(y_scale^3 * normalized_curve(X, Y/y_scale))
    assert curve == scaled_curve
    print("implicit specialization verified:", name)


def signed_relations(group):
    return tuple(tuple(relation.Tietze()) for relation in group.relations())


K.<root_six> = QuadraticField(6)
NumberRing.<NX, NY> = K[]
quadratic_curves_and_relations = tuple(
    (
        sign,
        universal_implicit(
            NumberRing,
            -24*(109*(-2 + sign*2*root_six/3)+40),
            -2*(2189*(-2 + sign*2*root_six/3)+802)/3,
            2*(3*(-2 + sign*2*root_six/3)+2),
        ),
        expected_relations,
    )
    for sign, expected_relations in (
        (1, P32_PLUS_RELATIONS),
        (-1, P32_MINUS_RELATIONS),
    )
)
for sign, _, _ in quadratic_curves_and_relations:
    h = -2 + sign*2*root_six/3
    assert 3*h^2+12*h+4 == 0

if "--algebra-only" not in sys.argv:
    for name, _, _, curve, expected_relations in rational_curves_and_relations:
        group = fundamental_group(
            curve,
            simplified=False,
            projective=False,
            puiseux=True,
        )
        assert signed_relations(group) == expected_relations
        print("rational residual presentation:", name, group)

    for sign, curve, expected_relations in quadratic_curves_and_relations:
        group = fundamental_group(
            curve,
            simplified=False,
            projective=False,
            puiseux=True,
        )
        assert signed_relations(group) == expected_relations
        print("quadratic residual presentation:", sign, group)

    print("A6 delta-five residual Sage certificates verified")
else:
    print("A6 delta-five residual Sage algebra certificates verified")
