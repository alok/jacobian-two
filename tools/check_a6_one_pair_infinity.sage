"""Reproduce the two Sage Zariski--van Kamp presentations used by the A6 audit."""

R.<X, Y> = QQ[]

generic_curve = -X^5 + 5*X^3*Y - 5*X*Y^2 + Y^3 + Y^2
contact_two_curve = (
    -X^5 - 100*X^4 - 10*X^3*Y + 40*X^2*Y
    + 15*X*Y^2 + Y^3 - 4*Y^2
)

generic_relations = (
    (3, 2, -3, -2),
    (3, 2, 3, 2, 3, -2, -3, -2, -3, -2),
    (-3, -2, -3, -2, 1, 2, 3, 2),
    (-3, -2, 1, 2, 3, -2, -1, 2),
)

contact_two_relations = (
    (
        -3, -2, -1, 2, 1, 2, 3, -2, -1, 2, 1, 2, 3, -2, -1,
        2, 1, 2, -3, -2, -1, -2, 1, 2, -3, -2, -1, -2, 1, 2,
    ),
    (-3, -2, -1, -2, -1, 2, 1, 2, 3, -2, -1, 2, 1, 2),
    (2, 1, 2, 1, -2, -1, -2, -1),
)


def signed_relations(group):
    return tuple(tuple(relation.Tietze()) for relation in group.relations())


G0 = Curve(generic_curve).fundamental_group(simplified=False)
G5 = Curve(contact_two_curve).fundamental_group(simplified=False)

assert signed_relations(G0) == generic_relations
assert signed_relations(G5) == contact_two_relations

G0_simplified = Curve(generic_curve).fundamental_group(simplified=True)
assert len(G0_simplified.gens()) == 1
assert len(G0_simplified.relations()) == 0

print("generic presentation:", G0)
print("generic simplified group:", G0_simplified)
print("contact-two presentation:", G5)
print("A6 one-pair Sage presentation certificate verified")
