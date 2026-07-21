"""Certify component completeness for the positive-dimensional delta-seven walls.

Run with::

    sage tools/check_a6_delta_seven_components.sage

The calculations are deliberately independent of the stored rational charts.
They start from unrestricted collision-factor ansaetze in the two affine
coefficient-slice equations, decompose or saturate the resulting ideals, and
then audit every locus on which a displayed chart denominator vanishes.
"""


def same_ideals(actual, expected):
    """Whether two short collections contain the same ideals, ignoring order."""

    if len(actual) != len(expected):
        return False
    unmatched = list(expected)
    for ideal in actual:
        for index, candidate in enumerate(unmatched):
            if ideal == candidate:
                unmatched.pop(index)
                break
        else:
            return False
    return not unmatched


def slice_equations(polynomial):
    """The two linear coefficient-slice equations for a monic septic."""

    coefficients = [polynomial[index] for index in range(6)]
    h0, h1, h2, h3, h4, h5 = coefficients
    return (
        h3 - 4 * h4 + 10 * h5 - 85,
        h0 - 2 * h1 + 3 * h2 - 11 * h4 + 34 * h5 - 306,
    )


def family_parameters(polynomial):
    """Recover ``alpha,beta,gamma,delta`` from a slice collision septic."""

    h0, h1, _, _, h4, h5 = [polynomial[index] for index in range(6)]
    alpha = -h0
    beta = h1 + h0
    delta = (h4 + h1 + h0 + 18 - 2 * h5) / 3
    gamma = h4 + h1 + h0 + 27 - 3 * h5
    return alpha, beta, gamma, delta


def wall_factors(parameters):
    """Return the cusp, extra-critical, and triple-collision factors."""

    alpha, beta, gamma, delta = parameters
    cusp = alpha + beta - gamma + delta - 1
    extra = 1215 * alpha + 756 * beta - 576 * gamma + 432 * delta - 320
    triple = (
        alpha^2
        - alpha * beta * gamma
        + 3 * alpha * beta * delta
        - alpha * beta
        + alpha * gamma^2
        - 3 * alpha * gamma * delta
        + 7 * alpha * delta
        - 9 * alpha
        - 2 * beta^2 * gamma
        + 6 * beta^2 * delta
        - 6 * beta^2
        + 5 * beta * gamma^2
        - 20 * beta * gamma * delta
        + 22 * beta * gamma
        + 15 * beta * delta^2
        - 32 * beta * delta
        + 17 * beta
        - 3 * gamma^3
        + 15 * gamma^2 * delta
        - 17 * gamma^2
        - 21 * gamma * delta^2
        + 46 * gamma * delta
        - 25 * gamma
        + 9 * delta^3
        - 29 * delta^2
        + 31 * delta
        - 11
    )
    return cusp, extra, triple


# ---------------------------------------------------------------------------
# Two unordered double roots: (2,2,1,1,1)
# ---------------------------------------------------------------------------

R2.<a, b, c1, c0> = PolynomialRing(QQ, order="degrevlex")
S2.<s> = PolynomialRing(R2)
quadratic = s^2 - a * s + b
residual_cubic = s^3 + (6 + 2 * a) * s^2 + c1 * s + c0
two_double = quadratic^2 * residual_cubic
two_f1, two_f2 = slice_equations(two_double)
two_ideal = R2.ideal([two_f1, two_f2])

two_denominator = (
    3 * a^3 + a^2 * b + 21 * a^2 + 10 * a * b
    + 42 * a + 2 * b^2 + 12 * b + 26
)
coefficient_matrix = matrix(
    R2,
    [
        [two_f1.derivative(c1), two_f1.derivative(c0)],
        [two_f2.derivative(c1), two_f2.derivative(c0)],
    ],
)
assert coefficient_matrix.det() == (a + b + 1) * two_denominator

# ``a+b+1=0`` means that the repeated quadratic has root ``-1`` and hence
# lies on the excluded cusp wall.  Saturation leaves the unique main surface.
two_main = two_ideal.saturation(R2.ideal([a + b + 1]))[0]
assert two_main.dimension() == 2
assert two_main.is_prime()
two_invalid = R2.ideal(
    [a + b + 1, 4 * b^2 + b * c1 - 13 * b - 3 * c1 + 2 * c0 + 11]
)
two_embedded = R2.ideal([a + 2, b - 1])
two_components = [primary.radical() for primary in two_ideal.primary_decomposition()]
assert same_ideals(two_components, [two_main, two_invalid, two_embedded])

# The main component's omitted denominator divisor has exactly three prime
# supports.  Two have the repeated quadratic ``(s+1)^2`` (so C=0); on the
# third, with roots -1/3 and -4/3, both L and T vanish identically.
two_boundary = (two_main + R2.ideal([two_denominator])).radical()
two_boundary_lt = R2.ideal([3 * a + 5, 9 * b - 4, c1 - 2 * c0 - 1])
two_boundary_c1 = R2.ideal([a + 2, b - 1, c1 - c0 - 1])
two_boundary_c2 = R2.ideal([a + 2, b - 1, c0])
two_boundary_components = [
    primary.radical() for primary in two_boundary.primary_decomposition()
]
assert same_ideals(
    two_boundary_components,
    [two_boundary_lt, two_boundary_c1, two_boundary_c2],
)
two_cusp, two_extra, two_triple = wall_factors(family_parameters(two_double))
assert two_extra.subs({a: -5/3, b: 4/9, c1: 2 * c0 + 1}) == 0
assert two_triple.subs({a: -5/3, b: 4/9, c1: 2 * c0 + 1}) == 0
assert two_cusp.subs({a: -2, b: 1, c1: c0 + 1}) == 0
assert two_cusp.subs({a: -2, b: 1, c0: 0}) == 0
print("two-double: one prime valid component; all determinant boundaries excluded")


# ---------------------------------------------------------------------------
# One quadruple root: (4,1,1,1)
# ---------------------------------------------------------------------------

R4.<v, d1, d0> = PolynomialRing(QQ, order="degrevlex")
S4.<s> = PolynomialRing(R4)
contact_four = (s - v)^4 * (s^3 + (6 + 4 * v) * s^2 + d1 * s + d0)
four_f1, four_f2 = slice_equations(contact_four)
four_ideal = R4.ideal([four_f1, four_f2])
four_main = four_ideal.saturation(R4.ideal([v + 1]))[0]
assert four_main.dimension() == 1
assert four_main.is_prime()
for prime in four_ideal.associated_primes():
    if prime != four_main:
        assert v + 1 in prime

# At the sole chart denominator, v=-13/3, the two slice equations are
# inconsistent.  Equivalently, the main prime plus that denominator is 1.
assert four_main + R4.ideal([3 * v + 13]) == R4.ideal([1])
four_cusp, _, _ = wall_factors(family_parameters(contact_four))
assert four_cusp.subs({v: -1}) == 0
print("contact-four: one prime valid curve; chart denominator has no point")


# ---------------------------------------------------------------------------
# One triple and one double root: (3,2,1,1)
# ---------------------------------------------------------------------------

RM.<h, k, q> = PolynomialRing(QQ, order="degrevlex")
SM.<s> = PolynomialRing(RM)
mixed_residual = s^2 + (6 + 3 * h + 2 * k) * s + q
mixed = (s - h)^3 * (s - k)^2 * mixed_residual
mixed_f1, mixed_f2 = slice_equations(mixed)
mixed_ideal = RM.ideal([mixed_f1, mixed_f2])

mixed_denominator = 3 * h^2 + 6 * h * k + 12 * h + k^2 + 8 * k + 10
mixed_numerator = (
    3 * h^4 + 18 * h^3 * k + 38 * h^3 + 18 * h^2 * k^2
    + 108 * h^2 * k + 132 * h^2 + 6 * h * k^3 + 66 * h * k^2
    + 204 * h * k + 180 * h + 8 * k^3 + 54 * k^2 + 120 * k + 85
)
mixed_parameter_curve = (
    3 * h^4 * k + 9 * h^4 + 18 * h^3 * k^2 + 83 * h^3 * k
    + 87 * h^3 + 18 * h^2 * k^3 + 168 * h^2 * k^2
    + 393 * h^2 * k + 273 * h^2 + 6 * h * k^4 + 108 * h * k^3
    + 444 * h * k^2 + 681 * h * k + 357 * h + 22 * k^4
    + 146 * k^3 + 366 * k^2 + 408 * k + 170
)
assert mixed_f1 == mixed_denominator * q - mixed_numerator
assert mixed_f1.resultant(mixed_f2, q) == (
    -(h + 1)^3 * (k + 1) * mixed_parameter_curve
)

# All components other than the main curve force h=-1 or k=-1, hence C=0.
mixed_main = mixed_ideal.saturation(RM.ideal([(h + 1) * (k + 1)]))[0]
assert mixed_main.dimension() == 1
assert mixed_main.is_prime()
for prime in mixed_ideal.associated_primes():
    if prime != mixed_main:
        assert (h + 1 in prime) or (k + 1 in prime)

# The plane parameter quintic is geometrically irreducible and its projective
# normalization has genus two.
AM.<hm, km> = PolynomialRing(QQ)
M = (
    3 * hm^4 * km + 9 * hm^4 + 18 * hm^3 * km^2 + 83 * hm^3 * km
    + 87 * hm^3 + 18 * hm^2 * km^3 + 168 * hm^2 * km^2
    + 393 * hm^2 * km + 273 * hm^2 + 6 * hm * km^4 + 108 * hm * km^3
    + 444 * hm * km^2 + 681 * hm * km + 357 * hm + 22 * km^4
    + 146 * km^3 + 366 * km^2 + 408 * km + 170
)
BM.<HM, KM> = PolynomialRing(QQbar)
assert len(BM(M).factor()) == 1
assert Curve(M).geometric_genus() == 2

# The dense rational chart q=N/D misses two genuine points, not a component.
# The third boundary point has h=k=-1 and is invalid.
mixed_boundary = (mixed_main + RM.ideal([mixed_denominator])).radical()
mixed_boundary_valid = RM.ideal(
    [22 * k - 43 * q + 46, 264 * h + 387 * q - 40,
     3741 * q^2 - 4704 * q + 1280]
)
mixed_boundary_invalid = RM.ideal([q, k + 1, h + 1])
assert mixed_boundary_valid.is_prime()
assert mixed_boundary == mixed_boundary_valid.intersection(mixed_boundary_invalid)

# Work in the quadratic residue field to prove that both conjugate points
# retain exactly the (3,2,1,1) partition and avoid C, L, and T.
U.<z> = PolynomialRing(QQ)
Kz = U.fraction_field()
z = Kz(z)
h_z = (40 - 387 * z) / 264
k_z = (43 * z - 46) / 22
boundary_polynomial = U(3741 * z^2 - 4704 * z + 1280)
assert boundary_polynomial.is_irreducible()
SZ.<s> = PolynomialRing(Kz)
residual_z = s^2 + (6 + 3 * h_z + 2 * k_z) * s + z
mixed_z = (s - h_z)^3 * (s - k_z)^2 * residual_z
cusp_z, extra_z, triple_z = wall_factors(family_parameters(mixed_z))
validity_values = (
    h_z - k_z,
    cusp_z,
    extra_z,
    triple_z,
    residual_z(h_z),
    residual_z(k_z),
    residual_z.discriminant(),
)
for value in validity_values:
    assert gcd(U(value.numerator()), boundary_polynomial) == 1
print("mixed: prime genus-two curve; dense chart misses two valid conjugate points")


# ---------------------------------------------------------------------------
# Three unordered double roots: (2,2,2,1)
# ---------------------------------------------------------------------------

R3.<e, f, g> = PolynomialRing(QQ, order="degrevlex")
S3.<s> = PolynomialRing(R3)
repeated_cubic = s^3 - e * s^2 + f * s + g
three_double = repeated_cubic^2 * (s + 6 + 2 * e)
three_f1, three_f2 = slice_equations(three_double)
three_ideal = R3.ideal([three_f1, three_f2])

three_good = R3.ideal(
    [9 * e + 3 * f - g + 17,
     36 * f^2 - 24 * f * g + 4 * g^2 - 105 * f + 8 * g + 85]
)
three_bad_cusp = R3.ideal(
    [e + f - g + 1, 4 * f^2 - 8 * f * g + 4 * g^2 - 13 * f + 12 * g + 11]
)
three_bad_triple = R3.ideal([e + 2, f - 1])
three_bad_cusp_degenerate = R3.ideal([e + g + 2, f - 2 * g - 1])
three_components = [
    primary.radical() for primary in three_ideal.primary_decomposition()
]
assert same_ideals(
    three_components,
    [three_good, three_bad_cusp, three_bad_triple, three_bad_cusp_degenerate],
)

# The sole component not contained in C=0 or T=0 is exactly the stored
# rational family.
stored_f = 4 * e^2 + 16 * e + 17
stored_g = 12 * e^2 + 57 * e + 68
assert three_good == R3.ideal([f - stored_f, g - stored_g])
assert three_bad_cusp == R3.ideal(
    [f - (4 * e^2 + 20 * e + 27), g - (4 * e^2 + 21 * e + 28)]
)
three_cusp, _, three_triple = wall_factors(family_parameters(three_double))
assert three_cusp.subs({f: stored_f, g: stored_g}) == -4 * (2 * e + 5)^5
assert three_triple.subs({f: stored_f, g: stored_g}) == (
    -4 * (e + 2)^8 * (12 * e + 29)^2
)
assert three_cusp.subs(
    {f: 4 * e^2 + 20 * e + 27, g: 4 * e^2 + 21 * e + 28}
) == 0
assert three_triple.subs({e: -2, f: 1}) == 0
assert three_cusp.subs({f: -2 * e - 3, g: -e - 2}) == 0
print("three-double: four prime supports; exactly one survives off C and T")

print("all positive-dimensional delta-seven component certificates verified")
