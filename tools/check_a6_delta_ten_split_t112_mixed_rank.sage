"""Certify the true-split ``T112`` and mixed rank strata exactly.

The systems are rebuilt independently from the ``V/W`` collision factors and
then cross-checked against the typed SymPy module.  Saturation removes only
repeated/critical sources, pair diagonals, component overlaps, and collisions
of the separate contact target with the triple target.  No determinant or
compatibility factor is inverted.
"""

from itertools import combinations
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.a6_delta_ten_split_t112_mixed_rank import (
    SPLIT_INCIDENCE_SPECS,
    split_incidence_matrices,
    valid_base_localizer,
)


Base = PolynomialRing(QQ, ["z", "w", "h"], order="degrevlex")
z, w, h = Base.gens()
Coefficient = PolynomialRing(Base, ["a", "b", "c", "d"])
a, b, c, d = Coefficient.gens()


def component_polynomials(kappa):
    """Return the true vertical and graph collision factors."""

    R = Coefficient["R"]
    r = R.gen()
    S = Coefficient["S"]
    s = S.gen()
    if kappa == 0:
        vertical = r^2 - c*r + a
        graph = (
            s^8
            + (8 - 2*c)*s^6
            + (-4*b + 8*d)*s^5
            + (-4*a + 2*c + 6)*s^4
            + (-8*b + 16*d)*s^3
            + (-16*a + 18*c - 16)*s^2
            + (12*b - 8*d)*s
            + 4*a - 2*c + 1
        )
    elif kappa == 2:
        vertical = (
            r^4 + (4*d - c - 10)*r^3
            + (a - 3*b - 10*d + 6*c + 15)*r^2
            + (-3*a + 4*b + 6*d - 5*c - 7)*r
            + a - b - d + c + 1
        )
        graph = (
            s^6 + 8*s^5 + (6 - 2*c + 8*d)*s^4
            + (-4*b + 16*d + 2*c - 16)*s^3
            + (-4*a - 8*b - 8*d + 18*c + 1)*s^2
            + (-16*a + 12*b - 2*c)*s + 4*a
        )
    else:
        raise ValueError("kappa must be 0 or 2")
    return R, vertical, S, graph


def coefficient_system(profile, kappa, allocation):
    """Build one true-component affine-linear incidence system."""

    R, vertical_in_r, S, graph_in_s = component_polynomials(kappa)
    r = R.gen()
    s = S.gen()
    if kappa == 0:
        base_constraint = z^2 + w^2 + 1
        vertical_root = -z^2
    else:
        base_constraint = z*(z + 1) + w*(w + 1)
        vertical_root = -z*(z + 1)
    graph_root = z + w
    vertical = Coefficient(vertical_in_r(vertical_root))
    graph = Coefficient(graph_in_s(graph_root))
    if profile == "T112":
        if allocation == "V":
            tangent = Coefficient(vertical_in_r.derivative(r)(vertical_root))
        else:
            tangent = Coefficient(graph_in_s.derivative(s)(graph_root))
        equations = (vertical, graph, tangent)
    elif profile == "mixed":
        if allocation == "V":
            contact_value = Coefficient(vertical_in_r(h))
            contact_derivative = Coefficient(vertical_in_r.derivative(r)(h))
        else:
            contact_value = Coefficient(graph_in_s(h))
            contact_derivative = Coefficient(graph_in_s.derivative(s)(h))
        equations = (vertical, graph, contact_value, contact_derivative)
    else:
        raise ValueError("unknown profile")

    matrix_rows = []
    right_hand_side = []
    zero = {a: 0, b: 0, c: 0, d: 0}
    for equation in equations:
        matrix_rows.append([Base(equation.derivative(parameter)) for parameter in (a,b,c,d)])
        right_hand_side.append(Base(-equation.subs(zero)))
    coefficient_matrix = matrix(Base, matrix_rows)
    augmented_matrix = coefficient_matrix.augment(vector(Base, right_hand_side))
    return base_constraint, coefficient_matrix, augmented_matrix


def minors(matrix_value, size):
    """Return all minors of a fixed size."""

    return [
        matrix_value.matrix_from_rows_and_columns(rows, columns).det()
        for rows in combinations(range(matrix_value.nrows()), size)
        for columns in combinations(range(matrix_value.ncols()), size)
    ]


def as_base(expression):
    """Convert a SymPy expression using this checker's short variable names."""

    text = str(expression).replace("**", "^")
    text = text.replace("z_split_t112_mixed", "z")
    text = text.replace("w_split_t112_mixed", "w")
    text = text.replace("h_split_t112_mixed", "h")
    return Base(text)


cases = {
    "t112_k0_v": ("T112", 0, "V"),
    "t112_k0_w": ("T112", 0, "W"),
    "t112_k2_v": ("T112", 2, "V"),
    "t112_k2_w": ("T112", 2, "W"),
    "mixed_k0_w": ("mixed", 0, "W"),
    "mixed_k2_v": ("mixed", 2, "V"),
    "mixed_k2_w": ("mixed", 2, "W"),
}

expected_t112 = {
    "t112_k0_v": (1, 1, 1),
    "t112_k0_w": (1, 1, 3),
    "t112_k2_v": (2, 2, 3),
    "t112_k2_w": (2, 2, 3),
}

expected_mixed = {
    "mixed_k0_w": {
        "prime_count": 1,
        "compatibility_unit": False,
        "compatibility_basis": (
            z^2 + 11/12,
            h^2 + 1/3,
            w - h/2,
        ),
    },
    "mixed_k2_v": {
        "prime_count": 1,
        "compatibility_unit": False,
        "compatibility_basis": (
            z^2 + z - h + 1/4,
            h^2 - h + 1/2,
            w - h + 3/2,
        ),
    },
    "mixed_k2_w": {
        "prime_count": 2,
        "compatibility_unit": True,
        "compatibility_basis": (Base(1),),
    },
}

python_specs = {spec.name: spec for spec in SPLIT_INCIDENCE_SPECS}

for name, (profile, kappa, allocation) in cases.items():
    base_constraint, coefficient_matrix, augmented_matrix = coefficient_system(
        profile, kappa, allocation
    )

    # Cross-check that the independent construction is exactly the system
    # exported by the typed Python certificate.
    python_matrix, python_augmented = split_incidence_matrices(
        python_specs[name]
    )
    assert coefficient_matrix == matrix(
        Base,
        python_matrix.rows,
        python_matrix.cols,
        [as_base(entry) for entry in python_matrix],
    )
    assert augmented_matrix == matrix(
        Base,
        python_augmented.rows,
        python_augmented.cols,
        [as_base(entry) for entry in python_augmented],
    )
    base_localizer = as_base(valid_base_localizer(python_specs[name]))
    valid = Base.ideal([base_localizer])
    generic_rank = coefficient_matrix.rank()
    assert generic_rank == python_specs[name].expected_rank

    rank_drop = Base.ideal(
        [base_constraint] + minors(coefficient_matrix, generic_rank)
    )
    rank_drop_valid, rank_drop_exponent = rank_drop.saturation(valid)
    compatibility = Base.ideal(
        [base_constraint] + minors(augmented_matrix, generic_rank)
    )
    compatibility_valid, compatibility_exponent = compatibility.saturation(
        valid
    )

    if profile == "T112":
        expected_rank_exponent, expected_compatibility_exponent, prime_count = (
            expected_t112[name]
        )
        assert rank_drop.dimension() == 1  # finite in (z,w), with unused h
        assert rank_drop_valid.is_one()
        assert rank_drop_exponent == expected_rank_exponent
        assert len(rank_drop.associated_primes()) == prime_count
        assert compatibility_valid.is_one()
        assert compatibility_exponent == expected_compatibility_exponent
        continue

    # A mixed rank-three divisor survives, but compatibility is much smaller.
    expected = expected_mixed[name]
    assert rank_drop_valid.dimension() == 1
    rank_drop_primes = tuple(rank_drop_valid.associated_primes())
    assert len(rank_drop_primes) == expected["prime_count"]
    assert all(prime.dimension() == 1 for prime in rank_drop_primes)
    assert rank_drop_exponent == 2
    assert compatibility_exponent == 2
    assert compatibility_valid.is_one() == expected["compatibility_unit"]
    assert tuple(compatibility_valid.groebner_basis()) == expected[
        "compatibility_basis"
    ]
    if not compatibility_valid.is_one():
        assert compatibility_valid.dimension() == 0
        assert compatibility_valid.vector_space_dimension() == 4
        assert compatibility_valid == compatibility_valid.radical()

    # A hidden threefold would require at least a curve of compatible rank-two
    # bases, or a compatible rank-one point.  Both rank-two ideals are empty
    # after the same valid localization.
    coefficient_rank_two = Base.ideal(
        [base_constraint] + minors(coefficient_matrix, 3)
    )
    coefficient_rank_two_valid, coefficient_rank_two_exponent = (
        coefficient_rank_two.saturation(valid)
    )
    augmented_rank_two = Base.ideal(
        [base_constraint] + minors(augmented_matrix, 3)
    )
    augmented_rank_two_valid, augmented_rank_two_exponent = (
        augmented_rank_two.saturation(valid)
    )
    assert coefficient_rank_two_valid.is_one()
    assert coefficient_rank_two_exponent == 1
    assert augmented_rank_two_valid.is_one()
    assert augmented_rank_two_exponent == 1


print("PASS")
print("split T112 valid rank drops: empty in four allocation charts")
print("mixed k=0/W compatibility: reduced length 4, rank 3, A1 fibers")
print("mixed k=2/V compatibility: reduced length 4, rank 3, A1 fibers")
print("mixed k=2/W compatibility: empty")
print("no hidden three-dimensional incidence component in the seven charts")
