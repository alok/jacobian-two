"""Exact sparse-support certificates for the residual ``(72, 108)`` case.

Guccione--Guccione--Horruitiner--Valqui, Proposition 4.3, reduces the only
coordinate-degree pair below 125 that remained in their 2022 analysis to two
Newton-polygon configurations.  In their transformed Laurent coordinates the
putative pair satisfies ``[P, Q] = x^2`` and has one of the two pairs of
polygons stored below.

This module proves deliberately sparse-support statements.  Case 1 requires
at least three coefficients strictly inside the two Newton polygons, while
Case 2 requires at least four.  Boundary coefficients are otherwise
unrestricted: all boundary lattice points are allowed, and only the polygon
vertices are required to be nonzero.

The proof is an exact finite exhaustion.  A monomial pair

    p_ij x^i y^j, q_kl x^k y^l

contributes ``(i*l - j*k) p_ij q_kl`` to the coefficient of
``x^(i+k-1) y^(j+l-1)`` in ``[P,Q]``.  Whenever a coefficient equation has
only one product left and one factor is known nonzero, the other factor must
vanish.  Iterating these implications eventually leaves a nonzero product as
the sole term in a coefficient that must be zero.  Every certificate is
replayed from rebuilt coefficient equations before it is counted.

Nothing here excludes the full Proposition 4.3 configurations, proves the
Jacobian conjecture, or concerns the separate generic sheet degree.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from math import gcd
from typing import Final, Literal, TypeAlias

import sympy as sp

Point: TypeAlias = tuple[int, int]
Side: TypeAlias = Literal["P", "Q"]

TARGET_BRACKET_EXPONENT: Final[Point] = (2, 0)


@dataclass(frozen=True, order=True)
class Coefficient:
    """A coefficient of ``P`` or ``Q``, identified by its exponent."""

    side: Side
    exponent: Point

    def label(self) -> str:
        """Return a stable compact label used in audit digests."""

        x_degree, y_degree = self.exponent
        return f"{self.side}[{x_degree},{y_degree}]"


@dataclass(frozen=True, order=True)
class BracketTerm:
    """One nonzero coefficient product in a Jacobian coefficient equation."""

    p: Coefficient
    q: Coefficient
    scalar: int

    def factors(self) -> tuple[Coefficient, Coefficient]:
        """Return the two coefficient variables in the product."""

        return self.p, self.q


@dataclass(frozen=True)
class ForceZero:
    """A one-term equation forcing one coefficient to vanish."""

    output: Point
    term: BracketTerm
    known_nonzero: Coefficient
    forced_zero: Coefficient


@dataclass(frozen=True)
class Contradiction:
    """A one-term zero equation whose two factors are known nonzero."""

    output: Point
    term: BracketTerm


@dataclass(frozen=True)
class ZeroProductCertificate:
    """A replayable sequence of zero implications ending in contradiction."""

    forced_zeros: tuple[ForceZero, ...]
    contradiction: Contradiction


@dataclass(frozen=True)
class NewtonConfiguration:
    """One of the two residual Newton-polygon configurations."""

    name: str
    p_vertices: tuple[Point, ...]
    q_vertices: tuple[Point, ...]

    def vertices(self) -> frozenset[Coefficient]:
        """Return the coefficients forced nonzero by the exact polygons."""

        return frozenset(
            Coefficient("P", point) for point in self.p_vertices
        ) | frozenset(Coefficient("Q", point) for point in self.q_vertices)


CASE_1: Final = NewtonConfiguration(
    name="Proposition 4.3, case 1",
    p_vertices=((0, 0), (1, 0), (8, 14), (8, 16), (0, 8)),
    q_vertices=((0, 0), (2, 1), (12, 21), (12, 24), (0, 12)),
)

CASE_2: Final = NewtonConfiguration(
    name="Proposition 4.3, case 2",
    p_vertices=((0, 0), (1, 0), (8, 14), (8, 16)),
    q_vertices=((0, 0), (2, 1), (12, 21), (12, 24)),
)

CONFIGURATIONS: Final[tuple[NewtonConfiguration, ...]] = (CASE_1, CASE_2)


@dataclass(frozen=True)
class AuditSummary:
    """Reproducible summary of one exhaustive sparse-interior audit."""

    configuration: str
    p_lattice_points: int
    p_boundary_points: int
    p_interior_points: int
    q_lattice_points: int
    q_boundary_points: int
    q_interior_points: int
    maximum_interior_terms: int
    support_patterns: int
    verified_certificates: int
    zero_product_certificates: int
    algebraic_certificates: int
    maximum_forced_zero_steps: int
    certificate_sha256: str


@dataclass(frozen=True)
class AlgebraicSupportCertificate:
    """A small exact ideal certificate for a zero-product exception.

    Zero-product propagation is run first.  The listed bracket coefficients,
    after the propagated zeros and the normalization ``p_10 = q_21 = 1``,
    generate the unit ideal once inverses for ``required_nonzero`` are added.
    Thus the certificate is replayed from the full coefficient equations; its
    conclusion is not accepted as a hard-coded exception.
    """

    selected_interior: tuple[Coefficient, ...]
    equation_outputs: tuple[Point, ...]
    required_nonzero: tuple[Coefficient, ...]


CASE_2_ALGEBRAIC_CERTIFICATES: Final[tuple[AlgebraicSupportCertificate, ...]] = (
    AlgebraicSupportCertificate(
        selected_interior=(
            Coefficient("P", (1, 1)),
            Coefficient("Q", (2, 2)),
            Coefficient("Q", (2, 3)),
        ),
        equation_outputs=((2, 1), (2, 2), (2, 3), (2, 4)),
        required_nonzero=(Coefficient("Q", (2, 3)),),
    ),
    AlgebraicSupportCertificate(
        selected_interior=(
            Coefficient("P", (1, 1)),
            Coefficient("Q", (2, 2)),
            Coefficient("Q", (3, 5)),
        ),
        equation_outputs=((2, 1), (3, 4), (3, 5)),
        required_nonzero=(Coefficient("Q", (2, 2)), Coefficient("Q", (3, 5))),
    ),
    AlgebraicSupportCertificate(
        selected_interior=(
            Coefficient("P", (1, 1)),
            Coefficient("Q", (2, 2)),
            Coefficient("Q", (5, 9)),
        ),
        equation_outputs=((2, 1), (5, 8), (5, 9)),
        required_nonzero=(Coefficient("Q", (2, 2)), Coefficient("Q", (5, 9))),
    ),
    AlgebraicSupportCertificate(
        selected_interior=(
            Coefficient("P", (2, 3)),
            Coefficient("Q", (3, 4)),
            Coefficient("Q", (4, 7)),
        ),
        equation_outputs=(
            (3, 3),
            (4, 4),
            (5, 7),
            (5, 6),
            (6, 9),
            (7, 11),
            (6, 8),
            (7, 10),
            (8, 13),
            (8, 12),
            (14, 24),
            (9, 14),
            (16, 28),
        ),
        required_nonzero=(Coefficient("P", (2, 3)), Coefficient("P", (8, 14))),
    ),
    AlgebraicSupportCertificate(
        selected_interior=(
            Coefficient("P", (4, 7)),
            Coefficient("Q", (5, 8)),
            Coefficient("Q", (8, 15)),
        ),
        equation_outputs=((3, 2), (5, 7), (6, 9), (9, 16)),
        required_nonzero=(Coefficient("P", (4, 7)), Coefficient("P", (8, 16))),
    ),
)


EquationMap: TypeAlias = dict[Point, tuple[BracketTerm, ...]]


def _cross(start: Point, end: Point, point: Point) -> int:
    """Return the oriented cross product ``(end-start) x (point-start)``."""

    return (end[0] - start[0]) * (point[1] - start[1]) - (end[1] - start[1]) * (
        point[0] - start[0]
    )


def _signed_double_area(vertices: tuple[Point, ...]) -> int:
    """Return twice the signed area of a polygon."""

    return sum(
        x * next_y - y * next_x
        for (x, y), (next_x, next_y) in zip(vertices, vertices[1:] + vertices[:1])
    )


def polygon_lattice_points(vertices: tuple[Point, ...]) -> frozenset[Point]:
    """Enumerate all lattice points in or on a convex counterclockwise polygon."""

    if len(vertices) < 3 or _signed_double_area(vertices) <= 0:
        msg = "vertices must describe a nondegenerate counterclockwise polygon"
        raise ValueError(msg)

    edges = tuple(zip(vertices, vertices[1:] + vertices[:1]))
    minimum_x = min(point[0] for point in vertices)
    maximum_x = max(point[0] for point in vertices)
    minimum_y = min(point[1] for point in vertices)
    maximum_y = max(point[1] for point in vertices)
    return frozenset(
        (x_degree, y_degree)
        for x_degree in range(minimum_x, maximum_x + 1)
        for y_degree in range(minimum_y, maximum_y + 1)
        if all(_cross(start, end, (x_degree, y_degree)) >= 0 for start, end in edges)
    )


def polygon_boundary_lattice_points(
    vertices: tuple[Point, ...],
) -> frozenset[Point]:
    """Enumerate the lattice points on the polygon boundary exactly."""

    result: set[Point] = set()
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        delta_x = end[0] - start[0]
        delta_y = end[1] - start[1]
        steps = gcd(abs(delta_x), abs(delta_y))
        if steps == 0:
            msg = "consecutive polygon vertices must be distinct"
            raise ValueError(msg)
        result.update(
            (
                start[0] + delta_x * step // steps,
                start[1] + delta_y * step // steps,
            )
            for step in range(steps + 1)
        )
    return frozenset(result)


def polygon_interior_lattice_points(
    vertices: tuple[Point, ...],
) -> frozenset[Point]:
    """Enumerate the lattice points strictly inside the polygon."""

    return polygon_lattice_points(vertices) - polygon_boundary_lattice_points(vertices)


def interior_coefficients(
    configuration: NewtonConfiguration,
) -> tuple[Coefficient, ...]:
    """Return every possible strictly interior coefficient in stable order."""

    p_coefficients = (
        Coefficient("P", point)
        for point in polygon_interior_lattice_points(configuration.p_vertices)
    )
    q_coefficients = (
        Coefficient("Q", point)
        for point in polygon_interior_lattice_points(configuration.q_vertices)
    )
    return tuple(sorted((*p_coefficients, *q_coefficients)))


def _bracket_term(p_point: Point, q_point: Point) -> tuple[Point, BracketTerm] | None:
    """Return the Jacobian contribution of one monomial pair, if nonzero."""

    p_x, p_y = p_point
    q_x, q_y = q_point
    scalar = p_x * q_y - p_y * q_x
    if scalar == 0:
        return None
    output = (p_x + q_x - 1, p_y + q_y - 1)
    return output, BracketTerm(
        p=Coefficient("P", p_point),
        q=Coefficient("Q", q_point),
        scalar=scalar,
    )


def bracket_equations(
    p_support: Iterable[Point], q_support: Iterable[Point]
) -> EquationMap:
    """Build every nonempty coefficient equation in ``[P,Q]`` exactly."""

    equations: defaultdict[Point, list[BracketTerm]] = defaultdict(list)
    for p_point in sorted(set(p_support)):
        for q_point in sorted(set(q_support)):
            contribution = _bracket_term(p_point, q_point)
            if contribution is not None:
                output, term = contribution
                equations[output].append(term)
    return {output: tuple(sorted(terms)) for output, terms in sorted(equations.items())}


def boundary_supports(
    configuration: NewtonConfiguration,
) -> tuple[frozenset[Point], frozenset[Point]]:
    """Return all lattice points on the two polygon boundaries."""

    return (
        polygon_boundary_lattice_points(configuration.p_vertices),
        polygon_boundary_lattice_points(configuration.q_vertices),
    )


def _supports_for_selected_interior(
    configuration: NewtonConfiguration,
    selected_interior: Iterable[Coefficient],
) -> tuple[frozenset[Point], frozenset[Point], frozenset[Coefficient]]:
    """Build an exact support ansatz and validate its selected interior terms."""

    selected = frozenset(selected_interior)
    possible = frozenset(interior_coefficients(configuration))
    if not selected <= possible:
        invalid = ", ".join(
            coefficient.label() for coefficient in sorted(selected - possible)
        )
        msg = f"selected coefficients are not strictly interior: {invalid}"
        raise ValueError(msg)

    p_boundary, q_boundary = boundary_supports(configuration)
    p_support = p_boundary | frozenset(
        coefficient.exponent for coefficient in selected if coefficient.side == "P"
    )
    q_support = q_boundary | frozenset(
        coefficient.exponent for coefficient in selected if coefficient.side == "Q"
    )
    return p_support, q_support, selected


def _active_terms(
    terms: tuple[BracketTerm, ...],
    known_zero: frozenset[Coefficient] | set[Coefficient],
) -> tuple[BracketTerm, ...]:
    """Remove products already killed by a known-zero factor."""

    return tuple(
        term for term in terms if term.p not in known_zero and term.q not in known_zero
    )


def _find_certificate(
    equations: Mapping[Point, tuple[BracketTerm, ...]],
    known_nonzero: frozenset[Coefficient],
) -> ZeroProductCertificate | None:
    """Find a deterministic zero-product contradiction, if this method can."""

    forced_zeros, contradiction = _propagate_zero_products(equations, known_nonzero)
    if contradiction is None:
        return None
    return ZeroProductCertificate(
        forced_zeros=forced_zeros,
        contradiction=contradiction,
    )


def _propagate_zero_products(
    equations: Mapping[Point, tuple[BracketTerm, ...]],
    known_nonzero: frozenset[Coefficient],
) -> tuple[tuple[ForceZero, ...], Contradiction | None]:
    """Run the deterministic propagation and retain its terminal zero set."""

    known_zero: set[Coefficient] = set()
    forced_zeros: list[ForceZero] = []

    while True:
        next_step: ForceZero | None = None
        for output in sorted(equations):
            if output == TARGET_BRACKET_EXPONENT:
                continue
            active = _active_terms(equations[output], known_zero)
            if len(active) != 1:
                continue
            term = active[0]
            p_nonzero = term.p in known_nonzero
            q_nonzero = term.q in known_nonzero
            if p_nonzero and q_nonzero:
                return (
                    tuple(forced_zeros),
                    Contradiction(output=output, term=term),
                )
            if p_nonzero:
                next_step = ForceZero(
                    output=output,
                    term=term,
                    known_nonzero=term.p,
                    forced_zero=term.q,
                )
                break
            if q_nonzero:
                next_step = ForceZero(
                    output=output,
                    term=term,
                    known_nonzero=term.q,
                    forced_zero=term.p,
                )
                break

        if next_step is None:
            return tuple(forced_zeros), None
        known_zero.add(next_step.forced_zero)
        forced_zeros.append(next_step)


def find_zero_product_certificate(
    configuration: NewtonConfiguration,
    selected_interior: Iterable[Coefficient] = (),
) -> ZeroProductCertificate | None:
    """Certify impossibility for an exact choice of nonzero interior terms.

    Every boundary lattice coefficient is allowed, but may vanish.  Polygon
    vertices and exactly the coefficients in ``selected_interior`` are known
    nonzero; all other interior coefficients are absent.
    """

    p_support, q_support, selected = _supports_for_selected_interior(
        configuration, selected_interior
    )
    equations = bracket_equations(p_support, q_support)
    return _find_certificate(equations, configuration.vertices() | selected)


def find_relaxed_full_polygon_certificate(
    configuration: NewtonConfiguration,
) -> ZeroProductCertificate | None:
    """Run the same method with every polygon lattice coefficient allowed.

    Only vertices are declared nonzero.  Returning ``None`` is an important
    scope check: it means this zero-product method does not decide the full
    Proposition 4.3 coefficient system.
    """

    equations = bracket_equations(
        polygon_lattice_points(configuration.p_vertices),
        polygon_lattice_points(configuration.q_vertices),
    )
    return _find_certificate(equations, configuration.vertices())


def _verify_certificate_against_equations(
    equations: Mapping[Point, tuple[BracketTerm, ...]],
    known_nonzero: frozenset[Coefficient],
    certificate: ZeroProductCertificate,
) -> bool:
    """Replay a certificate without trusting the searcher's stored state."""

    known_zero: set[Coefficient] = set()
    for step in certificate.forced_zeros:
        if step.output == TARGET_BRACKET_EXPONENT:
            return False
        terms = equations.get(step.output)
        if terms is None:
            return False
        active = _active_terms(terms, known_zero)
        if active != (step.term,) or step.term.scalar == 0:
            return False
        if step.known_nonzero not in known_nonzero:
            return False
        if step.known_nonzero not in step.term.factors():
            return False
        if step.forced_zero not in step.term.factors():
            return False
        if step.known_nonzero == step.forced_zero:
            return False
        if step.forced_zero in known_nonzero or step.forced_zero in known_zero:
            return False
        known_zero.add(step.forced_zero)

    contradiction = certificate.contradiction
    if contradiction.output == TARGET_BRACKET_EXPONENT:
        return False
    terms = equations.get(contradiction.output)
    if terms is None:
        return False
    active = _active_terms(terms, known_zero)
    return (
        active == (contradiction.term,)
        and contradiction.term.scalar != 0
        and contradiction.term.p in known_nonzero
        and contradiction.term.q in known_nonzero
    )


def verify_zero_product_certificate(
    configuration: NewtonConfiguration,
    selected_interior: Iterable[Coefficient],
    certificate: ZeroProductCertificate,
) -> bool:
    """Independently replay a sparse-support contradiction certificate."""

    p_support, q_support, selected = _supports_for_selected_interior(
        configuration, selected_interior
    )
    equations = bracket_equations(p_support, q_support)
    return _verify_certificate_against_equations(
        equations, configuration.vertices() | selected, certificate
    )


def _symbol_for(coefficient: Coefficient) -> sp.Symbol:
    """Return the stable SymPy variable used in algebraic replays."""

    x_degree, y_degree = coefficient.exponent
    return sp.Symbol(f"{coefficient.side.lower()}_{x_degree}_{y_degree}")


def verify_algebraic_support_certificate(
    configuration: NewtonConfiguration,
    certificate: AlgebraicSupportCertificate,
) -> bool:
    """Replay one exceptional support as an exact unit-ideal calculation."""

    selected = tuple(sorted(certificate.selected_interior))
    if selected != certificate.selected_interior:
        return False
    p_support, q_support, selected_set = _supports_for_selected_interior(
        configuration, selected
    )
    known_nonzero = configuration.vertices() | selected_set
    if not set(certificate.required_nonzero) <= known_nonzero:
        return False
    if len(set(certificate.equation_outputs)) != len(certificate.equation_outputs):
        return False
    if TARGET_BRACKET_EXPONENT in certificate.equation_outputs:
        return False

    equations = bracket_equations(p_support, q_support)
    expected_target = (
        BracketTerm(
            p=Coefficient("P", (1, 0)),
            q=Coefficient("Q", (2, 1)),
            scalar=1,
        ),
    )
    if equations.get(TARGET_BRACKET_EXPONENT) != expected_target:
        return False

    forced_zeros, contradiction = _propagate_zero_products(equations, known_nonzero)
    if contradiction is not None:
        return False
    known_zero = frozenset(step.forced_zero for step in forced_zeros)
    normalization = {
        _symbol_for(Coefficient("P", (1, 0))): sp.Integer(1),
        _symbol_for(Coefficient("Q", (2, 1))): sp.Integer(1),
    }

    polynomials: list[sp.Expr] = []
    for output in certificate.equation_outputs:
        terms = equations.get(output)
        if terms is None:
            return False
        expression = sp.Add(
            *(
                sp.Integer(term.scalar) * _symbol_for(term.p) * _symbol_for(term.q)
                for term in _active_terms(terms, known_zero)
            )
        ).xreplace(normalization)
        polynomials.append(sp.expand(expression))

    for index, coefficient in enumerate(certificate.required_nonzero):
        inverse = sp.Symbol(f"inverse_{index}")
        polynomials.append(inverse * _symbol_for(coefficient) - 1)

    symbols = tuple(
        sorted(
            set().union(*(polynomial.free_symbols for polynomial in polynomials)),
            key=str,
        )
    )
    if not symbols:
        return False
    basis = sp.groebner(polynomials, *symbols, order="grevlex", domain=sp.QQ)
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


def _certificate_record(
    selected: tuple[Coefficient, ...], certificate: ZeroProductCertificate
) -> str:
    """Serialize one verified certificate deterministically for an audit digest."""

    support = ",".join(coefficient.label() for coefficient in selected) or "boundary"
    force_records = (
        f"{step.output}:{step.term.scalar}:{step.known_nonzero.label()}"
        f"->{step.forced_zero.label()}"
        for step in certificate.forced_zeros
    )
    contradiction = certificate.contradiction
    contradiction_record = (
        f"X{contradiction.output}:{contradiction.term.scalar}:"
        f"{contradiction.term.p.label()}*{contradiction.term.q.label()}"
    )
    return "|".join((support, *force_records, contradiction_record))


def _algebraic_certificate_record(
    certificate: AlgebraicSupportCertificate,
) -> str:
    """Serialize one replayed exceptional-support certificate."""

    support = ",".join(
        coefficient.label() for coefficient in certificate.selected_interior
    )
    outputs = ",".join(str(output) for output in certificate.equation_outputs)
    nonzero = ",".join(
        coefficient.label() for coefficient in certificate.required_nonzero
    )
    return f"{support}|algebraic:{outputs}|nonzero:{nonzero}|unit-ideal"


def audit_sparse_interior_bound(
    configuration: NewtonConfiguration, maximum_interior_terms: int = 2
) -> AuditSummary:
    """Exhaust and verify every support with at most the requested interiors."""

    if maximum_interior_terms < 0:
        msg = "maximum_interior_terms must be nonnegative"
        raise ValueError(msg)

    possible = interior_coefficients(configuration)
    digest = sha256()
    support_patterns = 0
    verified_certificates = 0
    zero_product_certificates = 0
    algebraic_certificates = 0
    maximum_steps = 0
    algebraic_by_support = (
        {
            certificate.selected_interior: certificate
            for certificate in CASE_2_ALGEBRAIC_CERTIFICATES
        }
        if configuration == CASE_2
        else {}
    )

    for count in range(maximum_interior_terms + 1):
        for selected in combinations(possible, count):
            support_patterns += 1
            certificate = find_zero_product_certificate(configuration, selected)
            if certificate is None:
                algebraic = algebraic_by_support.get(selected)
                if algebraic is None:
                    continue
                if not verify_algebraic_support_certificate(configuration, algebraic):
                    msg = (
                        "internally invalid algebraic certificate for "
                        f"{configuration.name}"
                    )
                    raise AssertionError(msg)
                algebraic_certificates += 1
                verified_certificates += 1
                record = _algebraic_certificate_record(algebraic)
            else:
                if not verify_zero_product_certificate(
                    configuration, selected, certificate
                ):
                    msg = (
                        "internally invalid zero-product certificate for "
                        f"{configuration.name}"
                    )
                    raise AssertionError(msg)
                zero_product_certificates += 1
                verified_certificates += 1
                maximum_steps = max(maximum_steps, len(certificate.forced_zeros))
                record = _certificate_record(selected, certificate)
            digest.update(record.encode())
            digest.update(b"\n")

    p_lattice = polygon_lattice_points(configuration.p_vertices)
    p_boundary = polygon_boundary_lattice_points(configuration.p_vertices)
    q_lattice = polygon_lattice_points(configuration.q_vertices)
    q_boundary = polygon_boundary_lattice_points(configuration.q_vertices)
    return AuditSummary(
        configuration=configuration.name,
        p_lattice_points=len(p_lattice),
        p_boundary_points=len(p_boundary),
        p_interior_points=len(p_lattice - p_boundary),
        q_lattice_points=len(q_lattice),
        q_boundary_points=len(q_boundary),
        q_interior_points=len(q_lattice - q_boundary),
        maximum_interior_terms=maximum_interior_terms,
        support_patterns=support_patterns,
        verified_certificates=verified_certificates,
        zero_product_certificates=zero_product_certificates,
        algebraic_certificates=algebraic_certificates,
        maximum_forced_zero_steps=maximum_steps,
        certificate_sha256=digest.hexdigest(),
    )


def main() -> None:
    """Print the two exhaustive audit summaries."""

    audits = ((CASE_1, 2), (CASE_2, 3))
    for configuration, maximum_interior_terms in audits:
        summary = audit_sparse_interior_bound(configuration, maximum_interior_terms)
        print(summary.configuration)
        print(
            "  P lattice/boundary/interior: "
            f"{summary.p_lattice_points}/"
            f"{summary.p_boundary_points}/"
            f"{summary.p_interior_points}"
        )
        print(
            "  Q lattice/boundary/interior: "
            f"{summary.q_lattice_points}/"
            f"{summary.q_boundary_points}/"
            f"{summary.q_interior_points}"
        )
        print(
            "  verified supports with <= "
            f"{summary.maximum_interior_terms} interior terms: "
            f"{summary.verified_certificates}/{summary.support_patterns}"
        )
        print(
            "  zero-product/algebraic certificates: "
            f"{summary.zero_product_certificates}/"
            f"{summary.algebraic_certificates}"
        )
        print(f"  maximum forced-zero steps: {summary.maximum_forced_zero_steps}")
        print(f"  certificate sha256: {summary.certificate_sha256}")


if __name__ == "__main__":
    main()
