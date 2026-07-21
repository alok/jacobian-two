"""Place the clean true-split strata inside the global incidence components.

The split-rank certificates work directly with the vertical and graph factors
of the collision polynomial.  That is the right coordinate system for rank
computations, but by itself it does not say whether a split surface is a new
irreducible component.  This module supplies that missing comparison.

For ``T112`` the comparison is with the existing labeled triple-source
incidence.  Exact invertible row transformations identify all four clean
split systems with restrictions of that global rank-three system.  Hence the
clean split surfaces are already contained in the connected incidence used
by the proper Whitney--Thom propagation argument.

For ``C2+T111`` retain both the sum ``s`` and product ``p`` of the separate
contact pair.  If ``q,r`` are the second and third elementary symmetric
functions of the triple fiber, put

``n = r^2 + q^2 - q^3`` and ``k = n/(q*r)``.

After clearing ``q*r``, equality of the two first coordinates is

``G = (2*q*r*s+n)*p - s*(q*r*(s^2+1)+n*s)``.

The two coefficients of this polynomial in ``p`` are coprime.  Thus ``G`` is
primitive and linear in ``p``, so the total mixed source base is irreducible.
It includes the denominator/vertical pair at ``k=2`` instead of deleting it.
Exact row transformations then place the rank-four parts of all three mixed
split allocations in the same irreducible Cramer graph as the nonsplit
cyclic sample.

This closes topology only on the clean maximal-rank split loci.  The two
compatible rank-three schemes found by the focused rank audit have affine
line coefficient fibers.  Those lower-dimensional exceptional fibers still
require a closure or direct complement calculation; this module deliberately
does not promote them to excluded.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final, Literal

from sympy import (
    Expr,
    Poly,
    Rational,
    Symbol,
    cancel,
    diff,
    expand,
    gcd,
    rem,
    together,
)

from scripts.a6_delta_ten_contact_triple import (
    KAPPA_NUMERATOR,
    KAPPA_PARAMETERIZATION,
    TRIPLE_PAIR_SUM,
    TRIPLE_PRODUCT,
)
from scripts.a6_delta_ten_generic import (
    FAMILY_P,
    FAMILY_Q,
    KAPPA,
    T,
)
from scripts.a6_delta_ten_split_codim_two import (
    K0_FIBER_ONE,
    K2_FIBER_ONE,
    exact_split_witness_certificates,
    split_witness_specs,
)
from scripts.a6_delta_ten_split_t112_mixed_rank import (
    CONTACT_ROOT,
    SPLIT_INCIDENCE_SPECS,
    THIRD_ROOT,
    TRIPLE_ROOT,
    SplitIncidenceSpec,
    split_incidence_equations,
    split_triple_geometry,
)
from scripts.a6_delta_ten_t112 import exact_delta_ten_t112_certificate

CONTACT_SUM_TOTAL: Final = Symbol("s_total_mixed")
CONTACT_PRODUCT_TOTAL: Final = Symbol("p_total_mixed")

TOTAL_PAIR_PRODUCT_COEFFICIENT: Final = expand(
    2 * TRIPLE_PAIR_SUM * TRIPLE_PRODUCT * CONTACT_SUM_TOTAL + KAPPA_NUMERATOR
)
TOTAL_PAIR_CONSTANT_COEFFICIENT: Final = expand(
    -CONTACT_SUM_TOTAL
    * (
        TRIPLE_PAIR_SUM * TRIPLE_PRODUCT * (CONTACT_SUM_TOTAL**2 + 1)
        + KAPPA_NUMERATOR * CONTACT_SUM_TOTAL
    )
)
TOTAL_MIXED_BASE_EQUATION: Final = expand(
    TOTAL_PAIR_PRODUCT_COEFFICIENT * CONTACT_PRODUCT_TOTAL
    + TOTAL_PAIR_CONSTANT_COEFFICIENT
)

Component = Literal["V", "W"]


def _reduce_on_split_base(expression: Expr, constraint: Expr) -> Expr:
    """Reduce a polynomial modulo the monic quadratic split-fiber equation."""

    return expand(
        rem(
            Poly(expand(expression), THIRD_ROOT),
            Poly(constraint, THIRD_ROOT),
        ).as_expr()
    )


def _pair_target_tangent_rows(
    kappa: int,
    pair_sum: Expr,
    pair_product: Expr,
) -> tuple[Expr, Expr]:
    """Return symmetric target-equality and slope-equality rows for a pair."""

    pair_polynomial = T**2 - CONTACT_SUM_TOTAL * T + CONTACT_PRODUCT_TOTAL
    q_remainder = rem(FAMILY_Q, pair_polynomial, T)
    target_row = Poly(q_remainder, T).coeff_monomial(T)

    p_derivative_remainder = rem(
        diff(FAMILY_P.subs(KAPPA, kappa), T),
        pair_polynomial,
        T,
    )
    q_derivative_remainder = rem(diff(FAMILY_Q, T), pair_polynomial, T)
    p_linear = Poly(p_derivative_remainder, T).coeff_monomial(T)
    p_constant = Poly(p_derivative_remainder, T).coeff_monomial(1)
    q_linear = Poly(q_derivative_remainder, T).coeff_monomial(T)
    q_constant = Poly(q_derivative_remainder, T).coeff_monomial(1)
    tangent_row = expand(q_linear * p_constant - q_constant * p_linear)
    substitution = {
        CONTACT_SUM_TOTAL: pair_sum,
        CONTACT_PRODUCT_TOTAL: pair_product,
    }
    return (
        expand(target_row.subs(substitution)),
        expand(tangent_row.subs(substitution)),
    )


def _component_row_transform(
    kappa: int,
    component: Component,
    root: Expr,
) -> tuple[Expr, Expr, Expr]:
    """Return ``target=a*H`` and ``tangent=b*H+c*H'`` coefficients."""

    if kappa == 0 and component == "V":
        return (
            root**2,
            10 * root**2 * (2 * root - 1),
            4 * root**3 * (2 * root - 1),
        )
    if kappa == 0 and component == "W":
        return (
            Rational(1, 16),
            -(root**2) / 8,
            -root * (root**2 + 2) / 8,
        )
    if kappa == 2 and component == "V":
        return (1, 4 * root, 2 * root * (4 * root - 1))
    if kappa == 2 and component == "W":
        return (
            root**2 / 16,
            -(root**2) * (root + 1) * (3 * root + 5) / 8,
            -(root**3) * (root + 1) * (root + 2) / 8,
        )
    msg = f"unsupported split component transform: k={kappa}, {component}"
    raise ValueError(msg)


def _mixed_contact_coordinates(spec: SplitIncidenceSpec) -> tuple[Expr, Expr]:
    """Return sum and product coordinates of the separate contact pair."""

    h = CONTACT_ROOT
    if spec.kappa == 0 and spec.allocation == "contact-W":
        return h, (h**2 + 1) / 2
    if spec.kappa == 2 and spec.allocation == "contact-V":
        return Rational(-1), h
    if spec.kappa == 2 and spec.allocation == "contact-W":
        return h, h * (h + 1) / 2
    msg = f"unsupported mixed allocation: {spec.name}"
    raise ValueError(msg)


def _triple_total_coordinates(spec: SplitIncidenceSpec) -> tuple[Expr, Expr]:
    """Return ``q,r`` for the first three sources in a split four-source fiber."""

    first, second, third, _fourth = split_triple_geometry(spec.kappa).sources
    pair_sum = expand(first * second + first * third + second * third)
    product = expand(first * second * third)
    return pair_sum, product


@dataclass(frozen=True, slots=True)
class TotalMixedBaseCertificate:
    """Primitive-linear irreducibility certificate for the total source base."""

    degree_in_contact_product: int
    coefficient_gcd: Expr
    nonsplit_sample_identity: Expr

    @property
    def irreducible(self) -> bool:
        """Whether Gauss's lemma applies to the total base equation."""

        return self.degree_in_contact_product == 1 and self.coefficient_gcd == 1

    @property
    def absolutely_irreducible(self) -> bool:
        """Whether the same primitive-linear proof works over ``QQbar``.

        A common factor over the algebraic closure would give, after taking
        its finite Galois norm, a nonconstant common factor over ``QQ``.
        Thus the exact rational gcd-one certificate is also geometric.
        """

        return self.irreducible

    @property
    def verified(self) -> bool:
        """Whether irreducibility and the ordinary pair chart agree exactly."""

        return bool(
            self.irreducible
            and self.absolutely_irreducible
            and self.nonsplit_sample_identity == 0
        )


@cache
def exact_total_mixed_base_certificate() -> TotalMixedBaseCertificate:
    """Build the exact primitive-linear total-base certificate."""

    solved_product = cancel(
        -TOTAL_PAIR_CONSTANT_COEFFICIENT / TOTAL_PAIR_PRODUCT_COEFFICIENT
    )
    generic_pair_product = cancel(
        CONTACT_SUM_TOTAL
        * (CONTACT_SUM_TOTAL**2 + KAPPA_PARAMETERIZATION * CONTACT_SUM_TOTAL + 1)
        / (2 * CONTACT_SUM_TOTAL + KAPPA_PARAMETERIZATION)
    )
    return TotalMixedBaseCertificate(
        degree_in_contact_product=Poly(
            TOTAL_MIXED_BASE_EQUATION,
            CONTACT_PRODUCT_TOTAL,
        ).degree(),
        coefficient_gcd=expand(
            gcd(
                TOTAL_PAIR_PRODUCT_COEFFICIENT,
                TOTAL_PAIR_CONSTANT_COEFFICIENT,
            )
        ),
        nonsplit_sample_identity=cancel(solved_product - generic_pair_product),
    )


@dataclass(frozen=True, slots=True)
class SplitComponentEmbeddingCertificate:
    """Exact restriction of one split system from a global incidence chart."""

    name: str
    profile: str
    expected_rank: int
    global_row_residuals: tuple[Expr, ...]
    total_base_identity: Expr
    kappa_identity: Expr
    transformation_determinant: Expr
    witness_transformation_determinant: Expr
    witness_base_gradient: tuple[Expr, Expr, Expr, Expr]
    clean_witness_verified: bool

    @property
    def transformation_is_generically_invertible(self) -> bool:
        """Whether the split and global row systems agree on a dense open."""

        return bool(
            self.transformation_determinant != 0
            and self.witness_transformation_determinant != 0
        )

    @property
    def verified(self) -> bool:
        """Whether the restriction, nonempty open, and clean witness all agree."""

        return bool(
            all(residual == 0 for residual in self.global_row_residuals)
            and self.total_base_identity == 0
            and self.kappa_identity == 0
            and self.transformation_is_generically_invertible
            and any(value != 0 for value in self.witness_base_gradient)
            and self.clean_witness_verified
        )


def _witness_root_substitution(spec: SplitIncidenceSpec) -> dict[Expr, Expr]:
    """Return the exact split-base point used by the clean witness ledger."""

    fiber = K0_FIBER_ONE if spec.kappa == 0 else K2_FIBER_ONE
    return {
        TRIPLE_ROOT: fiber[0],
        THIRD_ROOT: fiber[2],
        CONTACT_ROOT: spec.witness_contact_root,
    }


@cache
def exact_split_component_embedding_certificates() -> tuple[
    SplitComponentEmbeddingCertificate, ...
]:
    """Compare all four ``T112`` and three mixed systems with global rows."""

    witness_specs = {spec.name: spec for spec in split_witness_specs()}
    witness_certificates = {
        certificate.name: certificate
        for certificate in exact_split_witness_certificates()
    }
    results: list[SplitComponentEmbeddingCertificate] = []

    for spec in SPLIT_INCIDENCE_SPECS:
        geometry = split_triple_geometry(spec.kappa)
        constraint = geometry.base_constraint
        vertical_root = geometry.vertical_root
        graph_root = geometry.graph_roots[0]
        vertical_component, graph_component = split_incidence_equations(spec)[:2]

        vertical_sum = expand(geometry.sources[0] + geometry.sources[1])
        vertical_product = expand(geometry.sources[0] * geometry.sources[1])
        graph_sum = expand(geometry.sources[0] + geometry.sources[2])
        graph_product = expand(geometry.sources[0] * geometry.sources[2])
        vertical_target, vertical_tangent = _pair_target_tangent_rows(
            spec.kappa,
            vertical_sum,
            vertical_product,
        )
        graph_target, graph_tangent = _pair_target_tangent_rows(
            spec.kappa,
            graph_sum,
            graph_product,
        )
        a_vertical, b_vertical, c_vertical = _component_row_transform(
            spec.kappa,
            "V",
            vertical_root,
        )
        a_graph, b_graph, c_graph = _component_row_transform(
            spec.kappa,
            "W",
            graph_root,
        )

        canonical_rows: list[Expr] = [vertical_target, graph_target]
        transformed_rows: list[Expr] = [
            a_vertical * vertical_component,
            a_graph * graph_component,
        ]
        transform_determinant = a_vertical * a_graph

        if spec.profile == "T112+6N":
            selected_derivative = split_incidence_equations(spec)[2]
            if spec.allocation == "tangent-V":
                canonical_rows.append(vertical_tangent)
                transformed_rows.append(
                    b_vertical * vertical_component + c_vertical * selected_derivative
                )
                transform_determinant *= c_vertical
            else:
                canonical_rows.append(graph_tangent)
                transformed_rows.append(
                    b_graph * graph_component + c_graph * selected_derivative
                )
                transform_determinant *= c_graph
            total_base_identity: Expr = Rational(0)
            kappa_identity: Expr = Rational(0)
            gradient = (
                diff(constraint, TRIPLE_ROOT),
                diff(constraint, THIRD_ROOT),
                Rational(0),
                Rational(0),
            )
            witness_substitution = _witness_root_substitution(spec)
            gradient = tuple(
                expand(value.subs(witness_substitution)) for value in gradient
            )
        else:
            contact_component, contact_derivative = split_incidence_equations(spec)[2:]
            contact_sum, contact_product = _mixed_contact_coordinates(spec)
            contact_target, contact_tangent = _pair_target_tangent_rows(
                spec.kappa,
                contact_sum,
                contact_product,
            )
            contact_kind: Component = "V" if spec.allocation == "contact-V" else "W"
            a_contact, b_contact, c_contact = _component_row_transform(
                spec.kappa,
                contact_kind,
                CONTACT_ROOT,
            )
            canonical_rows.extend((contact_target, contact_tangent))
            transformed_rows.extend(
                (
                    a_contact * contact_component,
                    b_contact * contact_component + c_contact * contact_derivative,
                )
            )
            transform_determinant *= a_contact * c_contact

            triple_pair_sum, triple_product = _triple_total_coordinates(spec)
            total_substitution = {
                TRIPLE_PAIR_SUM: triple_pair_sum,
                TRIPLE_PRODUCT: triple_product,
                CONTACT_SUM_TOTAL: contact_sum,
                CONTACT_PRODUCT_TOTAL: contact_product,
            }
            total_base_identity = _reduce_on_split_base(
                TOTAL_MIXED_BASE_EQUATION.subs(total_substitution),
                constraint,
            )
            kappa_numerator, kappa_denominator = together(
                KAPPA_PARAMETERIZATION.subs(
                    {
                        TRIPLE_PAIR_SUM: triple_pair_sum,
                        TRIPLE_PRODUCT: triple_product,
                    }
                )
                - spec.kappa
            ).as_numer_denom()
            del kappa_denominator
            kappa_identity = _reduce_on_split_base(kappa_numerator, constraint)
            gradient_expressions = tuple(
                diff(TOTAL_MIXED_BASE_EQUATION, variable).subs(total_substitution)
                for variable in (
                    TRIPLE_PAIR_SUM,
                    TRIPLE_PRODUCT,
                    CONTACT_SUM_TOTAL,
                    CONTACT_PRODUCT_TOTAL,
                )
            )
            witness_substitution = _witness_root_substitution(spec)
            gradient = tuple(
                expand(value.subs(witness_substitution))
                for value in gradient_expressions
            )

        row_residuals = tuple(
            _reduce_on_split_base(canonical - transformed, constraint)
            for canonical, transformed in zip(
                canonical_rows,
                transformed_rows,
                strict=True,
            )
        )
        witness_substitution = _witness_root_substitution(spec)
        witness_spec = witness_specs[spec.witness_name]
        clean_witness = witness_certificates[spec.witness_name]
        results.append(
            SplitComponentEmbeddingCertificate(
                name=spec.name,
                profile=spec.profile,
                expected_rank=spec.expected_rank,
                global_row_residuals=row_residuals,
                total_base_identity=total_base_identity,
                kappa_identity=kappa_identity,
                transformation_determinant=_reduce_on_split_base(
                    transform_determinant,
                    constraint,
                ),
                witness_transformation_determinant=expand(
                    transform_determinant.subs(witness_substitution)
                ),
                witness_base_gradient=gradient,
                clean_witness_verified=clean_witness.verified(witness_spec),
            )
        )
    return tuple(results)


@dataclass(frozen=True, slots=True)
class SplitComponentClosureCertificate:
    """Aggregate component-containment and clean-topology conclusion."""

    total_mixed_base: TotalMixedBaseCertificate
    embeddings: tuple[SplitComponentEmbeddingCertificate, ...]
    t112_global_incidence_excluded: bool
    representative_exceptional_rank_three_base_length: int
    exceptional_rank_three_fiber_topology_open: bool

    @property
    def maximal_rank_split_topology_closed(self) -> bool:
        """Whether all seven clean maximal-rank split loci are already covered."""

        t112 = tuple(item for item in self.embeddings if item.profile == "T112+6N")
        mixed = tuple(item for item in self.embeddings if item.profile == "C2+T111+5N")
        return bool(
            self.total_mixed_base.verified
            and len(t112) == 4
            and len(mixed) == 3
            and all(item.verified for item in self.embeddings)
            and all(item.expected_rank == 3 for item in t112)
            and all(item.expected_rank == 4 for item in mixed)
            and self.t112_global_incidence_excluded
        )

    @property
    def verified(self) -> bool:
        """Whether the maximal-rank theorem and remaining boundary are exact."""

        return bool(
            self.maximal_rank_split_topology_closed
            and self.representative_exceptional_rank_three_base_length == 8
            and self.exceptional_rank_three_fiber_topology_open
        )


@cache
def exact_split_component_closure_certificate() -> SplitComponentClosureCertificate:
    """Build the exact clean-split component comparison."""

    t112 = exact_delta_ten_t112_certificate()
    return SplitComponentClosureCertificate(
        total_mixed_base=exact_total_mixed_base_certificate(),
        embeddings=exact_split_component_embedding_certificates(),
        t112_global_incidence_excluded=t112.conditional_generic_exclusion_supported,
        # The k=0 contact-W and k=2 contact-V schemes each have length four.
        representative_exceptional_rank_three_base_length=8,
        exceptional_rank_three_fiber_topology_open=True,
    )


def main() -> int:
    """Print the component-containment conclusion and its remaining boundary."""

    certificate = exact_split_component_closure_certificate()
    print("total mixed base irreducible:", certificate.total_mixed_base.irreducible)
    print(
        "clean maximal-rank split topology closed:",
        certificate.maximal_rank_split_topology_closed,
    )
    print(
        "remaining exceptional rank-three base length (k=0,+2 reps):",
        certificate.representative_exceptional_rank_three_base_length,
    )
    print("remaining: topology of the compatible rank-three affine-line fibers")
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
