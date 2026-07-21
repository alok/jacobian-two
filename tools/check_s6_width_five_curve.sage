"""Reproduce the Sage complement group for the width-five S6 near-miss."""

from sage.schemes.curves.zariski_vankampen import fundamental_group

R.<p, q> = QQ[]
b = (
    p^7 - 2*p^6 + 9*p^5*q - 7*p^4*q^2 + 2*p^5 - 10*p^4*q
    + 28*p^3*q^2 - 34*p^2*q^3 + 12*p*q^4 - q^5 - 2*q^4
)

expected_relations = (
    (-5, 3, 5, -3),
    (-5, -3, -2, 1, 2, 3, 5, -3, -2, -1, 2, 3),
    (-5, -3, -2, -1, 2, 1, 2, 3, 5, -3, -2, -1, -2, 1, 2, 3),
    (-4, -3, -2, -1, -3, -2, -1, 2, 1, 2, 3, 1, 2, 3),
    (5, 1, 2, 3, 5, -3, -2, -1, -2, -1, 2, 1, 2, 3, -5, -3, -2, -1),
    (1, 2, 3, -5, -3, -2, -1, -5, -3, -2, -1, -3, -2, -1, 2, 1, 2, 3, 1, 2, 3, 5),
    (-3, -5, -3, -2, -1, -3, -2, -1, -2, 1, 2, 3, -2, 1, 2, -3, -2, -1, 2, 1, 2, 3, 1, 2, 3, 5),
    (-3, -2, -1, 2, 1, 2, 3, -2, -1, -2, 1, 2),
    (2, 1, -2, -1),
    (5, 4, -5, -4),
)

G = fundamental_group(b, simplified=False, projective=False, puiseux=True)
actual_relations = tuple(tuple(relation.Tietze()) for relation in G.relations())
assert actual_relations == expected_relations

isomorphism = G.simplification_isomorphism()
Z = isomorphism.codomain()
assert len(Z.gens()) == 1
assert len(Z.relations()) == 0
assert all(isomorphism(generator) == Z.gen(0) for generator in G.gens())
assert isomorphism.section()(Z.gen(0)) == G.gen(0)

print("width-five presentation:", G)
print("simplified group:", Z)
print("generator images:", [isomorphism(generator) for generator in G.gens()])
print("S6 width-five Sage complement certificate verified")
