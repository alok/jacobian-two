"""Independently verify split source bases and row transforms in Sage.

The Python certificates prove irreducibility by primitive-linear and flatness
arguments.  This checker asks Singular directly whether the total one-pair,
ordered two-pair, and mixed triple-plus-pair ideals are prime and confirms
their dimensions and relevant smooth loci.  It also saturates every split-row
transformation determinant by the complete valid-chart localizer.  It makes
no topology or Jacobian-conjecture claim.
"""


PairRing = PolynomialRing(QQ, ["k", "s", "p"], order="degrevlex")
k, s, p = PairRing.gens()
pair_equation = (2 * s + k) * p - s * (s^2 + k * s + 1)
pair_ideal = PairRing.ideal([pair_equation])

assert pair_ideal.is_prime()
assert pair_ideal.dimension() == 2


TwoPairRing = PolynomialRing(
    QQ,
    ["k", "s1", "p1", "s2", "p2"],
    order="degrevlex",
)
k, s1, p1, s2, p2 = TwoPairRing.gens()
first_pair = (2 * s1 + k) * p1 - s1 * (s1^2 + k * s1 + 1)
second_pair = (2 * s2 + k) * p2 - s2 * (s2^2 + k * s2 + 1)
two_pair_ideal = TwoPairRing.ideal([first_pair, second_pair])

assert two_pair_ideal.is_prime()
assert two_pair_ideal.dimension() == 3


MixedRing = PolynomialRing(QQ, ["q", "r", "s", "p"], order="degrevlex")
q, r, s, p = MixedRing.gens()
n = r^2 + q^2 - q^3
mixed_equation = (2 * q * r * s + n) * p - s * (q * r * (s^2 + 1) + n * s)
mixed_ideal = MixedRing.ideal([mixed_equation])

assert mixed_ideal.is_prime()
assert mixed_ideal.dimension() == 3

# The total mixed hypersurface is smooth on q*r != 0, which is the ordinary
# triple chart used by the global incidence.  This is stronger than merely
# proving that its defining polynomial is irreducible.
mixed_singular_ideal = MixedRing.ideal(
    [mixed_equation]
    + [mixed_equation.derivative(variable) for variable in MixedRing.gens()]
)
mixed_singular_valid, mixed_singular_exponent = mixed_singular_ideal.saturation(
    MixedRing.ideal([q * r])
)
assert mixed_singular_valid.is_one()
assert mixed_singular_exponent == 2

# The ordered two-pair complete intersection is not globally smooth.  Its
# entire singular scheme is supported at the three double-overlap fibers.
two_pair_jacobian = matrix(
    TwoPairRing,
    [
        [first_pair.derivative(variable) for variable in TwoPairRing.gens()],
        [second_pair.derivative(variable) for variable in TwoPairRing.gens()],
    ],
)
two_pair_singular_ideal = TwoPairRing.ideal(
    [first_pair, second_pair]
    + [
        two_pair_jacobian.matrix_from_columns(columns).det()
        for columns in Subsets(range(5), 2)
    ]
)
double_overlap_ideals = (
    TwoPairRing.ideal([k, s1, 2 * p1 - 1, s2, 2 * p2 - 1]),
    TwoPairRing.ideal([k - 2, s1 + 1, p1, s2 + 1, p2]),
    TwoPairRing.ideal([k + 2, s1 - 1, p1, s2 - 1, p2]),
)
expected_two_pair_singular_radical = double_overlap_ideals[0]
for point_ideal in double_overlap_ideals[1:]:
    expected_two_pair_singular_radical = (
        expected_two_pair_singular_radical.intersection(point_ideal)
    )
assert two_pair_singular_ideal.dimension() == 0
two_pair_singular_radical = two_pair_singular_ideal.radical()
assert two_pair_singular_radical == expected_two_pair_singular_radical
assert expected_two_pair_singular_radical.vector_space_dimension() == 3
assert two_pair_singular_radical.reduce(s1 - s2) == 0
assert two_pair_singular_radical.reduce(p1 - p2) == 0


print("PASS")
print("total pair base: prime, dimension 2")
print("ordered two-pair base: prime, dimension 3")
print("total mixed base: prime, dimension 3")
print(
    "total mixed q*r-localized singular ideal: unit, exponent",
    mixed_singular_exponent,
)
print("ordered two-pair singular support: three reduced double-overlap points")


# The Python row comparison records a nonzero transformation determinant for
# each split allocation.  Nonzero at one witness is not enough: the
# determinant must be a unit on the entire valid split localization.  Replay
# that assertion ideal-theoretically over the exact split base.  If
#
#   (<base equation>, <determinant>) : <valid localizer>^infinity = (1),
#
# then the determinant has no zero on the valid chart.
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.a6_delta_ten_split_component_closure import (
    _triple_total_coordinates,
    exact_split_component_embedding_certificates,
)
from scripts.a6_delta_ten_split_t112_mixed_rank import (
    SPLIT_INCIDENCE_SPECS,
    CONTACT_ROOT,
    THIRD_ROOT,
    TRIPLE_ROOT,
    split_triple_geometry,
    valid_base_localizer,
)


embedding_certificates = {
    certificate.name: certificate
    for certificate in exact_split_component_embedding_certificates()
}
transformation_unit_saturations = {}
mixed_qr_unit_saturations = {}
for spec in SPLIT_INCIDENCE_SPECS:
    variable_names = [str(TRIPLE_ROOT), str(THIRD_ROOT)]
    if spec.profile == "C2+T111+5N":
        variable_names.append(str(CONTACT_ROOT))
    SplitRing = PolynomialRing(QQ, variable_names, order="degrevlex")

    def as_split_polynomial(expression):
        return SplitRing(str(expression).replace("**", "^"))

    geometry = split_triple_geometry(spec.kappa)
    certificate = embedding_certificates[spec.name]
    determinant_zero_ideal = SplitRing.ideal(
        [
            as_split_polynomial(geometry.base_constraint),
            as_split_polynomial(certificate.transformation_determinant),
        ]
    )
    valid_ideal = SplitRing.ideal([as_split_polynomial(valid_base_localizer(spec))])
    saturated, exponent = determinant_zero_ideal.saturation(valid_ideal)
    assert saturated.is_one()
    transformation_unit_saturations[spec.name] = (exponent, True)

    if spec.profile == "C2+T111+5N":
        triple_pair_sum, triple_product = _triple_total_coordinates(spec)
        qr_zero_ideal = SplitRing.ideal(
            [
                as_split_polynomial(geometry.base_constraint),
                as_split_polynomial(triple_pair_sum * triple_product),
            ]
        )
        qr_saturated, qr_exponent = qr_zero_ideal.saturation(valid_ideal)
        assert qr_saturated.is_one()
        mixed_qr_unit_saturations[spec.name] = (qr_exponent, True)

assert transformation_unit_saturations == {
    "t112_k0_v": (10, True),
    "t112_k0_w": (4, True),
    "t112_k2_v": (2, True),
    "t112_k2_w": (4, True),
    "mixed_k0_w": (4, True),
    "mixed_k2_v": (2, True),
    "mixed_k2_w": (5, True),
}
assert mixed_qr_unit_saturations == {
    "mixed_k0_w": (4, True),
    "mixed_k2_v": (2, True),
    "mixed_k2_w": (2, True),
}

print("split transformation determinant saturations:", transformation_unit_saturations)
print("split mixed q*r unit saturations:", mixed_qr_unit_saturations)
print("topology computed: False")
print("proves JC(2): False")
