"""Replay the split ``C3``/``C2^2`` rank-stratum audit in Sage.

The Python certificate proves the localized ideal equalities by an exact
ideal sandwich.  This checker independently asks Singular for the actual
saturations, dimensions, lengths, radicals, and rank-two exclusions.
It makes no topology or plane-Jacobian-conjecture claim.
"""

from itertools import combinations
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.a6_delta_ten_split_contact_rank import (
    _linear_system,
    c3_rank_specs,
    residual_rank_specs,
    simple_double_contact_specs,
)


Base = PolynomialRing(QQ, ["u", "v"], order="lex")
u, v = Base.gens()


def as_base(expression):
    """Convert an exact SymPy expression to the Sage base ring."""

    return Base(
        str(expression)
        .replace("u_split", "u")
        .replace("v_split", "v")
        .replace("**", "^")
    )


def as_matrix(sympy_matrix):
    """Convert one exact SymPy matrix to a Sage matrix."""

    return matrix(
        Base,
        sympy_matrix.rows,
        sympy_matrix.cols,
        [as_base(entry) for entry in sympy_matrix],
    )


def minors(target, size):
    """Return all exact minors of a fixed size."""

    return [
        target.matrix_from_rows_and_columns(rows, columns).det()
        for rows in combinations(range(target.nrows()), size)
        for columns in combinations(range(target.ncols()), size)
    ]


# Every allowed C3 system has maximal rank on its exact clean root open.
for spec in c3_rank_specs():
    sympy_coefficient, _, _ = _linear_system(spec.equations)
    coefficient = as_matrix(sympy_coefficient)
    maximal_minors = minors(coefficient, coefficient.nrows())
    assert Base.ideal([gcd(maximal_minors)]) == Base.ideal(
        [as_base(spec.expected_maximal_minor_gcd)]
    )
    clean, _ = Base.ideal(maximal_minors).saturation(
        Base.ideal([as_base(spec.clean_localizer)])
    )
    assert clean.is_one()


# The three visible-factor C2^2 systems also have maximal rank on the clean
# open; no residual determinant remains in these allocations.
for spec in simple_double_contact_specs():
    sympy_coefficient, _, _ = _linear_system(spec.equations)
    coefficient = as_matrix(sympy_coefficient)
    maximal_minors = minors(coefficient, coefficient.nrows())
    assert Base.ideal([gcd(maximal_minors)]) == Base.ideal(
        [as_base(spec.expected_rank_drop_generator)]
    )
    clean, _ = Base.ideal(maximal_minors).saturation(
        Base.ideal([as_base(spec.clean_localizer)])
    )
    assert clean.is_one()


expected_compatibility_exponents = (2, 6, 1)
expected_rank_two_exponents = (2, 2, 1)
expected_lengths = (4, 6, 6)

for index, spec in enumerate(residual_rank_specs()):
    sympy_coefficient, _, sympy_augmented = _linear_system(spec.equations)
    coefficient = as_matrix(sympy_coefficient)
    augmented = as_matrix(sympy_augmented)
    visible = as_base(spec.determinant_visible_factor)
    residual, determinant_remainder = coefficient.det().quo_rem(visible)
    assert determinant_remainder == 0
    assert len(residual.factor()) == 1

    compatibility_minors = [
        augmented.matrix_from_columns(columns).det()
        for columns in combinations(range(augmented.ncols()), 4)
        if augmented.ncols() - 1 in columns
    ]
    compatibility_common = gcd(compatibility_minors)
    normalized_compatibility = []
    for compatibility_minor in compatibility_minors:
        quotient, remainder = compatibility_minor.quo_rem(compatibility_common)
        assert remainder == 0
        normalized_compatibility.append(quotient)

    raw_compatibility = Base.ideal([residual] + normalized_compatibility)
    clean_localizer = as_base(spec.clean_localizer)
    valid_compatibility, compatibility_exponent = raw_compatibility.saturation(
        Base.ideal([clean_localizer])
    )
    expected_compatibility = Base.ideal(
        [as_base(generator) for generator in spec.expected_compatibility_basis]
    )

    assert raw_compatibility.dimension() == 0
    assert compatibility_exponent == expected_compatibility_exponents[index]
    assert valid_compatibility == expected_compatibility
    assert valid_compatibility.dimension() == 0
    assert valid_compatibility.vector_space_dimension() == expected_lengths[index]
    assert valid_compatibility.radical() == valid_compatibility
    assert (valid_compatibility + Base.ideal([clean_localizer])).is_one()

    # On the residual determinant, rank <= 2 would require every 3 x 3
    # coefficient minor to vanish.  The localized ideal is the unit ideal,
    # so every compatible residual point has rank exactly three.
    coefficient_rank_two = Base.ideal(minors(coefficient, 3))
    coefficient_rank_two_valid, rank_two_exponent = (
        coefficient_rank_two.saturation(Base.ideal([clean_localizer]))
    )
    assert coefficient_rank_two_valid.is_one()
    assert rank_two_exponent == expected_rank_two_exponents[index]


print("PASS")
print("clean residual base lengths:", expected_lengths)
print("residual coefficient ranks: (3, 3, 3)")
print("coefficient fibers: affine lines")
print("maximum residual incidence dimension: 1")
print("topology computed: False")
print("proves JC(2): False")
