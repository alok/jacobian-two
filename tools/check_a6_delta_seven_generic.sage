"""Reproduce the generic delta-seven affine-complement certificate."""

import sys
from pathlib import Path

from sage.schemes.curves.zariski_vankampen import fundamental_group

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.a6_delta_seven_generic import GENERIC_RELATIONS


CollisionRing.<s> = QQ[]
H = s^7 + 6*s^6 + 9*s^5 - 5*s^3 + 2*s^2 + 2*s - 2
D = 30*s^7 + 170*s^6 + 250*s^5 + 20*s^4 - 110*s^3 + 30*s^2 + 30*s - 20
assert H.discriminant() == 464000000
assert H.resultant(-s*(3*s+4)) == -4220
assert H.resultant(s^2+s) == 2
assert H.resultant(D) == -1958080000000

NodeRing.<NX> = QQ[]
node_coordinates = (
    NX^7 - 24*NX^6 - 31*NX^5 - 30*NX^4 - 65*NX^3
    + 28*NX^2 + 18*NX + 2
)
assert node_coordinates.discriminant() == 113281250000000000
assert node_coordinates(0) == 2

R.<X, Y> = QQ[]
curve = (
    -X^10 - 10*X^8 + 15*X^6*Y - 20*X^6 - 4*X^5*Y - 4*X^5
    + 50*X^4*Y + 10*X^3*Y^2 + 10*X^3*Y - 25*X^2*Y^2
    + Y^3 + Y^2
)

ParameterRing.<t> = QQ[]
assert curve(X=t^2+t^3, Y=2*t^5+t^10) == 0
assert (2*t^5+t^10).derivative()(t=-2/3) == 33760/19683
assert (2*t^5+t^10)(t=-1) == -1

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert actual_relations == GENERIC_RELATIONS

isomorphism = group.simplification_isomorphism()
Z = isomorphism.codomain()
assert len(Z.gens()) == 1
assert len(Z.relations()) == 0
assert all(isomorphism(generator) == Z.gen(0) for generator in group.gens())
assert isomorphism.section()(Z.gen(0)) == group.gen(0)

print("delta-seven raw presentation:", group)
print("delta-seven simplified complement:", Z)
print("generic conditional A6 delta-seven Sage certificate verified")
