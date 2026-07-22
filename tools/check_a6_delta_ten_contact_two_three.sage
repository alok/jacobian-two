"""Reproduce the conditional delta-ten ``C2 + C3 + 5N`` certificate.

The first half derives and factors the 409-term compatibility hypersurface
directly from the five collision-jet equations, checks one smooth valid
rank-four point, and uses exact modular-reconstruction saturation over ``QQ``
to classify every lower-rank coefficient and augmented stratum.  The second
half independently regenerates the rational member, affine singular scheme,
and van Kamp presentation.  The checker certifies exact algebra and one
complement computation; propagation over a clean family is a separate
theorem-level dependency.
"""

from itertools import combinations
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sage.schemes.curves.zariski_vankampen import fundamental_group

from scripts.a6_delta_ten_contact_two_three import (
    CONTACT_THREE_SUM,
    CONTACT_TWO_SUM,
    CONTACT_TWO_THREE_AUGMENTED_MATRIX,
    CONTACT_TWO_THREE_COEFFICIENT_MATRIX,
    CONTACT_TWO_THREE_IMPLICIT,
    CONTACT_TWO_THREE_NODE_X_POLYNOMIAL,
    CONTACT_TWO_THREE_PAIR_DISJOINT_FACTOR,
    CONTACT_TWO_THREE_RELATIONS,
    CONTACT_TWO_THREE_VALID_LOCALIZER,
)
from scripts.a6_delta_ten_generic import KAPPA


Base = PolynomialRing(QQ, ["k", "u", "v"], order="degrevlex")
k, u, v = Base.gens()


def as_base(expression):
    """Convert an exact SymPy expression into the Sage base ring."""

    renamed = expression.subs(
        {
            KAPPA: KAPPA,
            CONTACT_TWO_SUM: CONTACT_TWO_SUM,
            CONTACT_THREE_SUM: CONTACT_THREE_SUM,
        }
    )
    return Base(str(renamed).replace("**", "^"))


coefficient_matrix = matrix(
    Base,
    5,
    4,
    [as_base(entry) for entry in CONTACT_TWO_THREE_COEFFICIENT_MATRIX],
)
augmented_matrix = matrix(
    Base,
    5,
    5,
    [as_base(entry) for entry in CONTACT_TWO_THREE_AUGMENTED_MATRIX],
)
valid_localizer = as_base(CONTACT_TWO_THREE_VALID_LOCALIZER)

pair_quadratic_u = u^2 + k*u + 1
pair_quadratic_v = v^2 + k*v + 1
known_boundary_factor = (
    -2 * (u-v)^6 * pair_quadratic_u * pair_quadratic_v^3
)
compatibility_residual, compatibility_remainder = (
    augmented_matrix.det().quo_rem(known_boundary_factor)
)
assert compatibility_remainder == 0
assert len(compatibility_residual.monomials()) == 409
assert compatibility_residual.degree() == 18
assert compatibility_residual.degree(k) == 13
assert compatibility_residual.degree(u) == 9
assert compatibility_residual.degree(v) == 10
compatibility_factorization = compatibility_residual.factor()
assert len(compatibility_factorization) == 1
assert compatibility_factorization[0][1] == 1
assert compatibility_factorization[0][0].degree() == compatibility_residual.degree()

sample_base = {k: QQ(2)/3, u: 1, v: -1}
assert compatibility_residual.subs(sample_base) == 0
sample_gradient = tuple(
    compatibility_residual.derivative(variable).subs(sample_base)
    for variable in (k, u, v)
)
assert sample_gradient == (
    QQ(1073741824)/59049,
    QQ(536870912)/59049,
    QQ(1073741824)/177147,
)
assert any(value != 0 for value in sample_gradient)
assert valid_localizer.subs(sample_base) != 0

sample_coefficient_matrix = coefficient_matrix.subs(sample_base)
sample_augmented_matrix = augmented_matrix.subs(sample_base)
assert sample_coefficient_matrix.rank() == 4
assert sample_augmented_matrix.rank() == 4
assert all(
    sample_coefficient_matrix.matrix_from_rows(rows).det() != 0
    for rows in combinations(range(5), 4)
)


# Strip only factors that are units on the explicitly declared valid chart.
# This is a localization, not a claim that any removed boundary is empty.
pair_disjoint_factor = as_base(CONTACT_TWO_THREE_PAIR_DISJOINT_FACTOR)
pair_diagonal_u = 2*u^2 + 3*k*u + 4
pair_diagonal_v = 2*v^2 + 3*k*v + 4
valid_factors = (
    k,
    k - 2,
    k + 2,
    u,
    v,
    u - v,
    k + 2*u,
    k + 2*v,
    pair_quadratic_u,
    pair_quadratic_v,
    pair_diagonal_u,
    pair_diagonal_v,
    k + u + v,
    pair_disjoint_factor,
)


def strip_valid_units(polynomial):
    """Remove every maximal valid-localizer power and rational content."""

    reduced = Base(polynomial)
    if reduced == 0:
        return reduced
    for factor in valid_factors:
        while True:
            quotient, remainder = reduced.quo_rem(factor)
            if remainder != 0:
                break
            reduced = quotient
    return Base(reduced / reduced.content())


def normalized_minors(target, size):
    """Return nonzero minors normalized by valid-chart units."""

    normalized = []
    for rows in combinations(range(target.nrows()), size):
        for columns in combinations(range(target.ncols()), size):
            minor = strip_valid_units(
                target.matrix_from_rows_and_columns(rows, columns).det()
            )
            if minor != 0:
                normalized.append(minor)
    return normalized


def exact_modular_saturation(ideal, factor):
    """Run Singular's exact modular-reconstruction saturation over ``QQ``."""

    singular.lib("modquotient.lib")
    answer = singular.modSat(
        singular(ideal),
        singular(Base.ideal([factor])),
    )
    return answer[1].sage(), answer[2].sage()


# The invertible linear change x=2u+k, y=2v+k cuts the residual from 409 to
# 299 terms and makes the exact characteristic-zero singular saturation
# tractable.  The transformed localizer differs from direct substitution
# only by nonzero rational units.
Shift = PolynomialRing(QQ, ["k_shift", "x_shift", "y_shift"], order="degrevlex")
k_shift, x_shift, y_shift = Shift.gens()
shifted_residual = Shift(
    compatibility_residual(
        k=k_shift,
        u=(x_shift-k_shift)/2,
        v=(y_shift-k_shift)/2,
    )
)
shifted_pair_quadratic_u = x_shift^2 - k_shift^2 + 4
shifted_pair_quadratic_v = y_shift^2 - k_shift^2 + 4
shifted_pair_diagonal_u = (
    x_shift^2 + k_shift*x_shift - 2*k_shift^2 + 8
)
shifted_pair_diagonal_v = (
    y_shift^2 + k_shift*y_shift - 2*k_shift^2 + 8
)


def shifted_target_numerator(value):
    """Return a cleared transformed collision-target numerator."""

    return -(value-k_shift)*(value+k_shift)*(value^2-k_shift^2+4)^2


def shifted_target_denominator(value):
    """Return a cleared transformed collision-target denominator."""

    return value^2


shifted_same_target, shifted_same_target_remainder = (
    (
        shifted_target_numerator(x_shift)
        * shifted_target_denominator(y_shift)
        - shifted_target_numerator(y_shift)
        * shifted_target_denominator(x_shift)
    ).quo_rem(x_shift-y_shift)
)
assert shifted_same_target_remainder == 0
shifted_valid_localizer = (
    k_shift
    * (k_shift^2-4)
    * (x_shift-k_shift)
    * (y_shift-k_shift)
    * (x_shift-y_shift)
    * x_shift
    * y_shift
    * shifted_pair_quadratic_u
    * shifted_pair_quadratic_v
    * shifted_pair_diagonal_u
    * shifted_pair_diagonal_v
    * shifted_same_target
)
shifted_singular_ideal = Shift.ideal(
    [shifted_residual]
    + [
        shifted_residual.derivative(variable)
        for variable in (k_shift, x_shift, y_shift)
    ]
)
singular.lib("modquotient.lib")
singular_answer = singular.modSat(
    singular(shifted_singular_ideal),
    singular(Shift.ideal([shifted_valid_localizer])),
)
valid_singular_saturation = singular_answer[1].sage()
valid_singular_exponent = singular_answer[2].sage()
assert valid_singular_saturation.is_one()
assert valid_singular_exponent == 2
print("[1/5] valid singular saturation complete", flush=True)


coefficient_rank_three_ideal = Base.ideal(
    normalized_minors(coefficient_matrix, 4)
)
augmented_rank_three_ideal = Base.ideal(
    normalized_minors(augmented_matrix, 4)
)
coefficient_rank_three_localizer = k*(k-2)*(k+2)*(u-v)
augmented_rank_three_localizer = (
    coefficient_rank_three_localizer * pair_diagonal_u * pair_diagonal_v
)
coefficient_rank_three_valid, coefficient_rank_three_exponent = (
    exact_modular_saturation(
        coefficient_rank_three_ideal,
        coefficient_rank_three_localizer,
    )
)
remaining_rank_three_localizer = prod(
    factor
    for index, factor in enumerate(valid_factors)
    if index not in (0, 1, 2, 5)
)
coefficient_rank_three_full_valid, remaining_rank_three_exponent = (
    exact_modular_saturation(
        coefficient_rank_three_valid,
        remaining_rank_three_localizer,
    )
)
augmented_rank_three_valid, augmented_rank_three_exponent = (
    exact_modular_saturation(
        augmented_rank_three_ideal,
        augmented_rank_three_localizer,
    )
)
assert not coefficient_rank_three_valid.is_one()
assert coefficient_rank_three_full_valid == coefficient_rank_three_valid
assert coefficient_rank_three_valid.dimension() == 1
assert coefficient_rank_three_exponent == 4
rank_three_hilbert = (
    coefficient_rank_three_valid.homogenize().hilbert_polynomial()
)
assert rank_three_hilbert.degree() == 1
assert rank_three_hilbert[1] == 30
assert rank_three_hilbert[0] == -51
print("[2/5] full-valid coefficient rank-three curve complete", flush=True)
assert augmented_rank_three_valid.is_one()
assert augmented_rank_three_exponent == 3
print("[3/5] augmented rank-three incompatibility complete", flush=True)

# Redundant lower-rank checks make the incompatibility conclusion hostile to
# accidental embedded rank-two strata.
coefficient_rank_two_ideal = Base.ideal(
    normalized_minors(coefficient_matrix, 3)
)
augmented_rank_two_ideal = Base.ideal(
    normalized_minors(augmented_matrix, 3)
)
rank_two_localizer = k*(k-2)*(k+2)
coefficient_rank_two_valid, coefficient_rank_two_exponent = (
    exact_modular_saturation(coefficient_rank_two_ideal, rank_two_localizer)
)
augmented_rank_two_valid, augmented_rank_two_exponent = (
    exact_modular_saturation(augmented_rank_two_ideal, rank_two_localizer)
)
assert coefficient_rank_two_valid.is_one()
assert coefficient_rank_two_exponent == 3
assert augmented_rank_two_valid.is_one()
assert augmented_rank_two_exponent == 3
print("[4/5] redundant rank-two saturations complete", flush=True)


# Regenerate the exact rational member over ZZ with X=3P and Y=Q.
Source.<source_X, source_Y, t> = ZZ[]
source_x = 3*t^2 + 2*t^3 + 3*t^4
source_y = 3*t^5 + t^6 + 3*t^7 + t^9
raw_implicit = (source_X - source_x).resultant(source_Y - source_y, t)
coefficient_content = gcd(
    abs(coefficient) for coefficient in raw_implicit.coefficients()
)
primitive_implicit = Source(raw_implicit / coefficient_content)
if primitive_implicit.monomial_coefficient(source_Y^4) < 0:
    primitive_implicit = -primitive_implicit

Plane.<X, Y> = QQ[]
source_to_plane = Source.hom([X, Y, 0], Plane)
regenerated_curve = source_to_plane(primitive_implicit)
curve = Plane(str(CONTACT_TWO_THREE_IMPLICIT).replace("**", "^"))
assert coefficient_content == 1
assert regenerated_curve == curve
assert len(curve.factor()) == 1

ParameterRing.<z> = QQ[]
assert curve(
    X=3*z^2 + 2*z^3 + 3*z^4,
    Y=3*z^5 + z^6 + 3*z^7 + z^9,
) == 0
assert curve(X=-5, Y=3) == 0
assert curve(X=-1, Y=-1) == 0
assert curve(X=0, Y=0) == 0

group = fundamental_group(
    curve,
    simplified=False,
    projective=False,
    puiseux=True,
)
actual_relations = tuple(tuple(relation.Tietze()) for relation in group.relations())
assert len(group.gens()) == 4
assert len(actual_relations) == 10
assert actual_relations == CONTACT_TWO_THREE_RELATIONS

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

jacobian = Plane.ideal([curve, curve.derivative(X), curve.derivative(Y)])
assert jacobian.dimension() == 0
assert jacobian.vector_space_dimension() == 17
assert jacobian.radical().vector_space_dimension() == 8

node_ideal = None
contact_three_ideal = Plane.ideal([X + 1, Y + 1])
contact_two_ideal = Plane.ideal([X + 5, Y - 3])
cusp_ideal = Plane.ideal([X, Y])
component_data = []
for component in jacobian.primary_decomposition():
    radical = component.radical()
    lengths = (
        component.vector_space_dimension(),
        radical.vector_space_dimension(),
    )
    component_data.append((component, radical, lengths))
    if radical not in (contact_three_ideal, contact_two_ideal, cusp_ideal):
        node_ideal = component

assert tuple(lengths for _, _, lengths in component_data) == (
    (5, 5),
    (5, 1),
    (3, 1),
    (4, 1),
)
assert component_data[1][1] == contact_three_ideal
assert component_data[2][1] == contact_two_ideal
assert component_data[3][1] == cusp_ideal
assert node_ideal is not None
assert node_ideal == node_ideal.radical()

NodeRing.<node_x> = QQ[]
node_x_polynomial = NodeRing(
    str(CONTACT_TWO_THREE_NODE_X_POLYNOMIAL)
    .replace("X", "node_x")
    .replace("**", "^")
)
assert node_x_polynomial.is_irreducible()
assert node_x_polynomial.discriminant() != 0
node_to_plane = NodeRing.hom([X], Plane)
assert node_ideal.reduce(node_to_plane(node_x_polynomial)) == 0
print("[5/5] sample curve and complement reconstruction complete", flush=True)

print("PASS")
print("compatibility residual terms:", len(compatibility_residual.monomials()))
print("compatibility residual degrees:", (
    compatibility_residual.degree(),
    compatibility_residual.degree(k),
    compatibility_residual.degree(u),
    compatibility_residual.degree(v),
))
print("compatibility residual irreducible over QQ: True")
print("smooth rational compatibility point:", sample_base)
print(
    "valid compatibility singular saturation: unit, exponent",
    valid_singular_exponent,
)
print(
    "coefficient rank <= 3 full-valid curve: degree 30, exponents",
    (coefficient_rank_three_exponent, remaining_rank_three_exponent),
)
print(
    "augmented rank <= 3 valid saturation: unit, exponent",
    augmented_rank_three_exponent,
)
print(
    "coefficient/augmented rank <= 2 saturations: units, exponents",
    (coefficient_rank_two_exponent, augmented_rank_two_exponent),
)
print("Jacobian component lengths:", tuple(
    lengths for _, _, lengths in component_data
))
print("affine complement: Z")
