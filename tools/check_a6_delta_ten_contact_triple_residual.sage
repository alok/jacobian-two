"""Certify the mixed contact--triple residual-rank boundary exactly.

This checker studies the coefficient matrix over the three-dimensional
``(q,r,w)`` base used by ``a6_delta_ten_contact_triple.py``.  It deliberately
does not make a topology claim: the only question here is whether the
residual determinant divisor can carry a three-dimensional incidence
component after the coefficient rank drops.
"""

from itertools import combinations
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.a6_delta_ten_contact_triple import (
    CONTACT_DENOMINATOR_FACTOR,
    CONTACT_DIAGONAL_FACTOR,
    CONTACT_PRODUCT_FACTOR,
    INCIDENCE_MATRIX,
    OVERLAP_FACTOR_C,
    OVERLAP_FACTOR_E,
    RESIDUAL_RANK_FACTOR,
    SPLIT_FACTOR,
    TRIPLE_DISCRIMINANT_NUMERATOR,
    TRIPLE_FOURTH_ROOT_FACTOR,
)


Base = PolynomialRing(
    QQ,
    ["q_contact_triple", "r_contact_triple", "w_contact_triple"],
    order="degrevlex",
)
q, r, w = Base.gens()


def as_base(expression):
    """Convert one exact SymPy expression into the Sage base ring."""

    return Base(str(expression).replace("**", "^"))


coefficient_matrix = matrix(
    Base,
    4,
    4,
    [as_base(entry) for entry in INCIDENCE_MATRIX],
)
contact_denominator = as_base(CONTACT_DENOMINATOR_FACTOR)

# The first two rows have the nonzero contact-denominator factor in common.
# Removing it does not alter rank on the valid contact chart and makes the
# determinantal calculations substantially smaller.
for row in (0, 1):
    for column in range(4):
        quotient, remainder = coefficient_matrix[row, column].quo_rem(
            contact_denominator
        )
        assert remainder == 0
        coefficient_matrix[row, column] = quotient

residual_rank_factor = as_base(RESIDUAL_RANK_FACTOR)
valid_localizer = prod(
    as_base(expression)
    for expression in (
        SPLIT_FACTOR,
        CONTACT_DENOMINATOR_FACTOR,
        CONTACT_DIAGONAL_FACTOR,
        CONTACT_PRODUCT_FACTOR,
        OVERLAP_FACTOR_C,
        OVERLAP_FACTOR_E,
        TRIPLE_DISCRIMINANT_NUMERATOR,
        TRIPLE_FOURTH_ROOT_FACTOR,
    )
) * q * r * w

# Rank at most two forces every 3 x 3 minor to vanish.  It is enough to use
# the eight minors containing the two triple rows: their common zero set is
# an over-approximation of the true rank-two locus.
rank_two_minors = []
for contact_row in (0, 1):
    rows = (contact_row, 2, 3)
    for columns in combinations(range(4), 3):
        minor = coefficient_matrix.matrix_from_rows_and_columns(rows, columns).det()
        rank_two_minors.append(minor)

rank_two_ideal = Base.ideal([residual_rank_factor] + rank_two_minors)
rank_two_valid, rank_two_saturation_exponent = rank_two_ideal.saturation(
    Base.ideal([valid_localizer])
)

# The two triple-equality rows already have rank two on the valid triple
# chart.  Their six 2 x 2 minors therefore suffice to rule out coefficient
# rank one.
triple_rows = coefficient_matrix.matrix_from_rows([2, 3])
rank_one_minors = [
    triple_rows.matrix_from_columns(columns).det()
    for columns in combinations(range(4), 2)
]
rank_one_ideal = Base.ideal(rank_one_minors)
rank_one_valid, rank_one_saturation_exponent = rank_one_ideal.saturation(
    Base.ideal(
        [
            q
            * r
            * as_base(TRIPLE_DISCRIMINANT_NUMERATOR)
            * as_base(TRIPLE_FOURTH_ROOT_FACTOR)
        ]
    )
)

assert rank_two_ideal.dimension() == 1
assert rank_two_valid.is_one()
assert rank_two_saturation_exponent == 9
assert rank_one_valid.is_one()
assert rank_one_saturation_exponent == 10

print("PASS")
print("rank <= 2 raw dimension:", rank_two_ideal.dimension())
print("rank <= 2 valid saturation: unit, exponent", rank_two_saturation_exponent)
print("rank <= 1 valid saturation: unit, exponent", rank_one_saturation_exponent)
