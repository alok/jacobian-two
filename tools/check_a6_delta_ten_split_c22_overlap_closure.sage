"""Certify that the split ``C2^2`` overlap-plus-``W`` locus is a closure.

The ordinary ordered two-contact incidence uses the total unordered-pair
surface

``F(k,s,p) = (2*s+k)*p - s*(s^2+k*s+1)``.

Let ``C`` be the coefficient of ``t`` in ``Q`` modulo ``t^2-s*t+p`` and let
``D = F_p*d/ds - F_s*d/dp``.  For two ordinary pairs, the equations are
``C_1=D(C)_1=C_2=D(C)_2=0``.  This checker constructs an exact rational arc
in that global incidence whose first pair tends to the vertical/graph
intersection at ``k=0`` while its second pair tends to a clean ``W`` pair.

The limiting two-parameter family is compared with the complete affine-line
solution of the true split ``overlap+W`` system.  The comparison is dominant,
so the entire irreducible overlap incidence lies in the closure of the
ordinary two-contact component.  This is algebraic containment only; the
checker makes no complement-topology or Jacobian-conjecture claim.
"""

from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.a6_delta_ten_split_contact_rank import (
    simple_double_contact_specs,
)


# Build the universal affine-linear pair rows independently over QQ[s,p,k].
Universal = PolynomialRing(QQ, ["s", "p", "k"])
s, p, k = Universal.gens()

remainders = [Universal(0), Universal(1)]
for _degree in range(2, 10):
    remainders.append(s * remainders[-1] - p * remainders[-2])

# Entries are coefficients of (a,b,c,d,1) in C.
target_row = vector(
    Universal,
    [
        remainders[5],
        remainders[6],
        remainders[7],
        remainders[8],
        remainders[9],
    ],
)
fiber_derivative = vector(
    Universal,
    [
        (2 * s + k) * entry.derivative(s)
        - (2 * p - 3 * s^2 - 2 * k * s - 1) * entry.derivative(p)
        for entry in target_row
    ],
)


# The rational arc.  The first pair tends to (s,p)=(0,1/2); its pair
# equation determines k.  The second pair is the ordinary W graph pair.
ArcBase = PolynomialRing(QQ, ["epsilon", "lambda", "v"])
epsilon, arc_parameter, v = ArcBase.gens()
ArcField = ArcBase.fraction_field()
epsilon, arc_parameter, v = map(ArcField, (epsilon, arc_parameter, v))

first_sum = epsilon
first_product = QQ(1) / 2 + arc_parameter * epsilon
arc_kappa = first_sum * (
    first_sum^2 + 1 - 2 * first_product
) / (first_product - first_sum^2)
second_sum = v
second_product = second_sum * (
    second_sum^2 + arc_kappa * second_sum + 1
) / (2 * second_sum + arc_kappa)


def specialize_pair_row(row, pair_sum, pair_product):
    """Specialize one universal coefficient row to the arc field."""

    substitution = {s: pair_sum, p: pair_product, k: arc_kappa}
    return vector(ArcField, [ArcField(entry.subs(substitution)) for entry in row])


arc_rows = (
    specialize_pair_row(target_row, first_sum, first_product),
    specialize_pair_row(fiber_derivative, first_sum, first_product),
    specialize_pair_row(target_row, second_sum, second_product),
    specialize_pair_row(fiber_derivative, second_sum, second_product),
)
coefficient_matrix = matrix(ArcField, [list(row[:4]) for row in arc_rows])
right_side = vector(ArcField, [-row[4] for row in arc_rows])

assert coefficient_matrix.det() != 0
arc_solution = coefficient_matrix.solve_right(right_side)
assert all(
    row[:4] * arc_solution + row[4] == 0
    for row in arc_rows
)


def epsilon_valuation(polynomial):
    """Return the least epsilon exponent of a nonzero ArcBase polynomial."""

    assert polynomial != 0
    return min(monomial[0] for monomial in polynomial.dict())


LimitBase = PolynomialRing(QQ, ["lambda", "v"])
limit_parameter, limit_v = LimitBase.gens()
LimitField = LimitBase.fraction_field()


def arc_limit(expression):
    """Take the exact epsilon -> 0 limit of a regular rational function."""

    numerator = expression.numerator()
    denominator = expression.denominator()
    numerator_order = epsilon_valuation(numerator)
    denominator_order = epsilon_valuation(denominator)
    assert numerator_order == denominator_order
    numerator_lead = sum(
        coefficient * limit_parameter^monomial[1] * limit_v^monomial[2]
        for monomial, coefficient in numerator.dict().items()
        if monomial[0] == numerator_order
    )
    denominator_lead = sum(
        coefficient * limit_parameter^monomial[1] * limit_v^monomial[2]
        for monomial, coefficient in denominator.dict().items()
        if monomial[0] == denominator_order
    )
    assert denominator_lead != 0
    return LimitField(numerator_lead) / LimitField(denominator_lead)


limit_solution = vector(LimitField, [arc_limit(entry) for entry in arc_solution])
lam = LimitField(limit_parameter)
u = LimitField(limit_v)
limit_denominator = u^5 + 5 * u^3 - 8 * lam + 10 * u
expected_limit = vector(
    LimitField,
    [
        (
            3 * u^7 + 14 * u^5 + 21 * u^3 - 8 * lam + 14 * u
        )
        / (4 * limit_denominator),
        (
            -3 * lam * u^8
            + u^9
            - 13 * lam * u^6
            - u^7
            - 9 * lam * u^4
            - 21 * u^5
            + 19 * lam * u^2
            - 35 * u^3
            - 2 * lam
        )
        / (4 * u * limit_denominator),
        (
            3 * u^7
            + 15 * u^5
            + 26 * u^3
            - 16 * lam
            + 24 * u
        )
        / (2 * limit_denominator),
        (
            -3 * lam * u^8
            + 3 * u^9
            - 13 * lam * u^6
            - 3 * u^7
            + 5 * lam * u^4
            - 63 * u^5
            + 65 * lam * u^2
            - 105 * u^3
            - 6 * lam
        )
        / (8 * u * limit_denominator),
    ],
)
assert limit_solution == expected_limit


# Convert the independently generated true-split overlap equations to the
# limit field and verify that the arc limit satisfies all three.
overlap_spec = next(
    spec
    for spec in simple_double_contact_specs()
    if spec.allocation == "overlap+W"
)


def as_limit_expression(expression):
    """Convert a SymPy split equation to the Sage limit field."""

    rendered = str(expression)
    # Replace coefficient names as tokens, longest structural names first.
    rendered = rendered.replace("v_split", "v").replace("**", "^")
    evaluation = sage_eval(
        rendered,
        locals={
            "a": expected_limit[0],
            "b": expected_limit[1],
            "c": expected_limit[2],
            "d": expected_limit[3],
            "v": u,
        },
    )
    return LimitField(evaluation)


assert all(as_limit_expression(equation) == 0 for equation in overlap_spec.equations)


# The fourth limiting coefficient is a nonconstant fractional-linear
# function of lambda.  More strongly, solve it for lambda in terms of an
# arbitrary desired affine-line coordinate eta and recover the complete
# rank-three split solution.  This proves dominance, not merely one arc.
DominanceBase = PolynomialRing(QQ, ["eta", "v"])
eta, dominance_v = DominanceBase.gens()
DominanceField = DominanceBase.fraction_field()
eta = DominanceField(eta)
dominance_v = DominanceField(dominance_v)
lambda_for_eta = dominance_v^2 * (
    -8 * eta * dominance_v^4
    - 40 * eta * dominance_v^2
    - 80 * eta
    + 3 * dominance_v^7
    - 3 * dominance_v^5
    - 63 * dominance_v^3
    - 105 * dominance_v
) / (
    -64 * eta * dominance_v
    + 3 * dominance_v^8
    + 13 * dominance_v^6
    - 5 * dominance_v^4
    - 65 * dominance_v^2
    + 6
)


def to_dominance_field(expression):
    """Substitute lambda_for_eta and v into one limiting coefficient."""

    numerator = LimitBase(expression.numerator())
    denominator = LimitBase(expression.denominator())
    substitution = {
        limit_parameter: lambda_for_eta,
        limit_v: dominance_v,
    }
    return DominanceField(numerator.subs(substitution)) / DominanceField(
        denominator.subs(substitution)
    )


dominant_solution = vector(
    DominanceField,
    [to_dominance_field(entry) for entry in expected_limit],
)
assert dominant_solution[3] == eta

split_denominator = dominance_v^6 + 5 * dominance_v^4 - 5 * dominance_v^2 + 15
complete_split_solution = vector(
    DominanceField,
    [
        -16 * eta * dominance_v / split_denominator
        + (
            3 * dominance_v^8
            + 14 * dominance_v^6
            - 70 * dominance_v^2
            + 21
        )
        / (4 * split_denominator),
        eta
        * (
            2 * dominance_v^6
            + 10 * dominance_v^4
            + 10 * dominance_v^2
            + 10
        )
        / split_denominator
        + (
            -dominance_v^9
            + dominance_v^7
            + 21 * dominance_v^5
            + 35 * dominance_v^3
        )
        / (2 * split_denominator),
        -32 * eta * dominance_v / split_denominator
        + (
            3 * dominance_v^8
            + 15 * dominance_v^6
            + 5 * dominance_v^4
            - 75 * dominance_v^2
            + 36
        )
        / (2 * split_denominator),
        eta,
    ],
)
assert dominant_solution == complete_split_solution


print("PASS")
print("ordinary two-contact arc determinant nonzero: True")
print("limit equals complete overlap+W incidence: True")
print("dominant parameters: (v, eta)")
print("topology computed: False")
print("proves JC(2): False")
