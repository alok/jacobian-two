"""Reproduce the three excess-three delta-seven complement presentations."""

import sys
from pathlib import Path

from sage.schemes.curves.zariski_vankampen import fundamental_group

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.a6_delta_seven_deeper_wall import (
    CONTACT_FOUR_RELATIONS,
    MIXED_RELATIONS,
    THREE_DOUBLE_RELATIONS,
)


R.<X, Y> = QQ[]

contact_four = (
    64*X^10 - 146225*X^9 + 868720*X^8 - 833440*X^7
    - 36540*X^6*Y - 3091200*X^6 + 250320*X^5*Y - 1467648*X^5
    - 257600*X^4*Y - 2800*X^3*Y^2 - 349440*X^3*Y
    + 17360*X^2*Y^2 - 21120*X*Y^2 - 64*Y^3 + 3328*Y^2
)

mixed = (
    -7625597484987*X^10 + 3094601589800*X^9 - 472284263760*X^8
    + 33943539520*X^7 - 176547210*X^6*Y - 1158502400*X^6
    + 43414680*X^5*Y + 14992384*X^5 - 3202400*X^4*Y
    + 27600*X^3*Y^2 + 77440*X^3*Y - 8295*X^2*Y^2
    + 640*X*Y^2 + Y^3 - 16*Y^2
)

three_double = (
    -32768*X^10 + 114055*X^9 - 131320*X^8 + 55840*X^7
    - 2445*X^6*Y - 9840*X^6 + 4970*X^5*Y + 576*X^5
    - 1980*X^4*Y + 5*X^3*Y^2 + 240*X^3*Y
    - 90*X^2*Y^2 + 35*X*Y^2 + Y^3 - 4*Y^2
)

samples = (
    ("contact four", contact_four, CONTACT_FOUR_RELATIONS),
    ("contact three plus contact two", mixed, MIXED_RELATIONS),
    ("three contact twos", three_double, THREE_DOUBLE_RELATIONS),
)

for name, curve, expected_relations in samples:
    group = fundamental_group(
        curve,
        simplified=False,
        projective=False,
        puiseux=True,
    )
    actual_relations = tuple(
        tuple(relation.Tietze()) for relation in group.relations()
    )
    assert actual_relations == expected_relations

    isomorphism = group.simplification_isomorphism()
    simplified = isomorphism.codomain()
    assert len(simplified.gens()) == 1
    assert len(simplified.relations()) == 0
    assert all(isomorphism(generator) == simplified.gen(0) for generator in group.gens())
    print(name, "raw presentation:", group)
    print(name, "simplified complement:", simplified)

print("all excess-three delta-seven Sage certificates verified")
