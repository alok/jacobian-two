"""Reproduce the three A6 delta-five affine complement presentations."""

R.<X, Y> = QQ[]

ParameterRing.<alpha, beta, gamma> = QQ[]
CollisionRing.<s> = ParameterRing[]
H = (
    s^5 + (gamma+2)*s^4 + (4*gamma-2)*s^3
    + (-alpha+2*beta+3*gamma-4)*s^2
    + (-alpha+3*beta-gamma)*s + alpha
)
discriminant_factors = list(H.discriminant().factor())
assert len(discriminant_factors) == 2
assert all(exponent == 1 for _, exponent in discriminant_factors)
assert sorted(factor.total_degree() for factor, _ in discriminant_factors) == [1, 7]
assert any(
    factor == alpha-beta+gamma-1 or factor == -(alpha-beta+gamma-1)
    for factor, _ in discriminant_factors
)
degree_seven_factor = next(
    factor for factor, _ in discriminant_factors if factor.total_degree() == 7
)
AbsoluteRing = PolynomialRing(QQbar, names=("A", "B", "Gamma"))
A, B, Gamma = AbsoluteRing.gens()
absolute_degree_seven = AbsoluteRing(
    degree_seven_factor(alpha=A, beta=B, gamma=Gamma)
)
absolute_factors = list(absolute_degree_seven.factor())
assert len(absolute_factors) == 1
assert absolute_factors[0][1] == 1
assert absolute_factors[0][0].total_degree() == 7

triple_divisor = (
    alpha*gamma-alpha-2*beta*gamma+3*beta
    + 2*gamma^2-5*gamma+3
)
assert len(list(triple_divisor.factor())) == 1

curves_and_relations = (
    (
        -X^8 - 5*X^7 - 6*X^6 + 8*X^5*Y + X^5 + 19*X^4*Y
        - 3*X^3*Y - 15*X^2*Y^2 + 9*X*Y^2 + Y^3 - Y^2,
        (
            (3, 2, -3, -2),
            (-3, -2, 3, 1),
            (-3, -2, -3, 2, 3, 1, -3, -2, 3, 2, 3, -1),
            (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
            (3, 2, -3, -2),
            (-3, -2, -3, -2, 1, 2, 3, 2, 3, -2, -3, -2, -1, 2, 3, 2),
            (-3, -2, 1, 2, 3, -2, -1, 2),
        ),
    ),
    (
        -X^8 - 9*X^7 - 27*X^6 + 8*X^5*Y - 18*X^5 + 41*X^4*Y
        + 30*X^3*Y - 12*X^2*Y^2 - 7*X*Y^2 + Y^3 + 2*Y^2,
        (
            (3, 2, -3, -2),
            (2, 3, 1, -3, -2, -1),
            (2, -1, -3, -2, 3, 1),
            (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
            (-3, -2, -3, -2, 1, 2, 3, 2),
            (-3, -2, 1, 2, 3, -2, -1, 2),
        ),
    ),
    (
        -X^8 + X^7 - 31*X^6 - X^5*Y + 48*X^5 + 2*X^4*Y
        + 66*X^3*Y + 6*X^2*Y^2 + 13*X*Y^2 + Y^3 - 3*Y^2,
        (
            (2, 1, -2, -1),
            (3, 2, -3, -2),
            (-3, -2, -1, 2, 3, 2, -3, -2, 1, 2, 3, 2, -3, -2, 1, 2, 3, 2, -3, -2, -1, 2, 3, -2, -3, -2, -1, 2, 3, -2),
            (-3, -2, -1, 2, 3, -2, -3, -2, -1, 2, 3, -2, 1, 2, 3, 2),
            (-3, -2, -1, 2, -3, -2, -1, 2, 3, -2, 1, 2, 3, 2, -3, -2, -1, 2, -3, -2, 1, 2, 3, -2, 1, 2, 3, -2),
            (-3, -2, 1, 2, 3, -2, 1, 2, 3, -2, -1, 2, -3, -2, -1, 2),
        ),
    ),
)

for curve, expected_relations in curves_and_relations:
    G = Curve(curve).fundamental_group(simplified=False)
    actual_relations = tuple(tuple(relation.Tietze()) for relation in G.relations())
    assert actual_relations == expected_relations
    isomorphism = G.simplification_isomorphism()
    Z = isomorphism.codomain()
    assert len(Z.gens()) == 1
    assert len(Z.relations()) == 0
    assert all(isomorphism(generator) == Z.gen(0) for generator in G.gens())
    assert isomorphism.section()(Z.gen(0)) == G.gen(0)
    print("presentation:", G)
    print("simplified group:", Z)

print("A6 delta-five Sage complement certificates verified")
