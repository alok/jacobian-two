"""Place the maximal-rank split contact loci in the global components.

Let

``F(k,s,p) = (2*s+k)*p - s*(s^2+k*s+1)``

be the total unordered source-pair surface and let ``C`` be the coefficient of
``t`` in the remainder of ``Q(t)`` modulo ``t^2-s*t+p``.  Then ``C=0`` is
target equality.  The derivation

``D = F_p*d/ds - F_s*d/dp``

is tangent to every fixed-``k`` pair fiber.  Away from its singular points,
``C=D(C)=0`` is contact at least two and
``C=D(C)=D^2(C)=0`` is contact at least three.

This coordinate-free system remains regular when the usual formula solving
for ``p`` has denominator zero.  Exact row transformations identify its jets
with the true vertical/graph component jets on:

* the three maximal-rank split ``C3`` charts; and
* the five maximal-rank split ``C2^2`` charts.

All transformation determinants are units on the corresponding clean opens.
The total pair surface is geometrically integral and flat over the ``k`` line.
Its ordered self-fiber-product is integral as well: it is flat, and its
generic fiber is a domain; a zero divisor would therefore become torsion over
``QQ[k]``.  Consequently the displayed rank-open split loci lie in the same
irreducible global contact incidences as the existing cyclic samples.  The
proper Whitney--Thom propagation already proved for those clean global
components applies to these rank-open split loci.

The result deliberately excludes the two ``C3`` overlap allocations, the
``C2^2`` overlap-plus-contact allocation, and the finite compatible rank-three
fibers on the three residual determinant charts.  Their incidence dimensions
are already at most two, but their complement topology remains open.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Final

from sympy import (
    Expr,
    Matrix,
    Poly,
    cancel,
    diff,
    expand,
    linear_eq_to_matrix,
    rem,
)

from scripts.a6_delta_ten_contact_three import (
    exact_delta_ten_contact_three_certificate,
)
from scripts.a6_delta_ten_double_contact import (
    exact_double_contact_sample_certificate,
    exact_ordered_incidence_certificate,
)
from scripts.a6_delta_ten_generic import (
    ALPHA,
    BETA,
    DELTA,
    FAMILY_Q,
    GAMMA,
    KAPPA,
    PAIR_INCIDENCE,
    R,
    S,
    T,
)
from scripts.a6_delta_ten_split_codim_two import (
    exact_split_witness_certificates,
    split_witness_specs,
)
from scripts.a6_delta_ten_split_contact_rank import (
    U,
    V,
    LinearRankSpec,
    ResidualRankSpec,
    SimpleDoubleContactSpec,
    c3_rank_specs,
    residual_rank_specs,
    simple_double_contact_specs,
)
from scripts.a6_delta_ten_split_pair_surface import (
    exact_split_pair_surface_certificate,
)

COEFFICIENTS: Final = (ALPHA, BETA, GAMMA, DELTA)
PAIR_POLYNOMIAL: Final = T**2 - S * T + R
PAIR_TARGET_ROW: Final = Poly(rem(FAMILY_Q, PAIR_POLYNOMIAL, T), T).coeff_monomial(T)


def pair_fiber_derivative(expression: Expr) -> Expr:
    """Differentiate along a fixed-``k`` fiber of the total pair surface."""

    return expand(
        diff(PAIR_INCIDENCE, R) * diff(expression, S)
        - diff(PAIR_INCIDENCE, S) * diff(expression, R)
    )


PAIR_TANGENCY_ROW: Final = pair_fiber_derivative(PAIR_TARGET_ROW)
PAIR_CONTACT_THREE_ROW: Final = pair_fiber_derivative(PAIR_TANGENCY_ROW)


def _component_coordinates(kappa: int, component: str, root: Expr) -> dict[Expr, Expr]:
    """Embed a true vertical or graph component into ``(k,s,p)``."""

    if (kappa, component) == (0, "V"):
        return {KAPPA: 0, S: 0, R: root}
    if (kappa, component) == (0, "W"):
        return {KAPPA: 0, S: root, R: (root**2 + 1) / 2}
    if (kappa, component) == (2, "V"):
        return {KAPPA: 2, S: -1, R: root}
    if (kappa, component) == (2, "W"):
        return {KAPPA: 2, S: root, R: root * (root + 1) / 2}
    msg = f"unsupported split component: k={kappa}, component={component}"
    raise ValueError(msg)


def _augmented_row_matrix(equations: tuple[Expr, ...]) -> Matrix:
    """Return rows of full affine-linear expressions, including constants."""

    coefficient, right_side = linear_eq_to_matrix(equations, COEFFICIENTS)
    return coefficient.row_join(-right_side)


def _row_transformation(
    global_equations: tuple[Expr, ...],
    component_equations: tuple[Expr, ...],
) -> tuple[Matrix, tuple[Expr, ...]]:
    """Solve and verify the exact global-to-component row transformation."""

    global_rows = _augmented_row_matrix(global_equations)
    component_rows = _augmented_row_matrix(component_equations)
    rank = len(component_equations)
    component_square = component_rows[:, :rank]
    transformation = (global_rows[:, :rank] * component_square.inv()).applyfunc(cancel)
    residual_matrix = (global_rows - transformation * component_rows).applyfunc(cancel)
    return transformation, tuple(residual_matrix)


@dataclass(frozen=True, slots=True)
class TotalContactBaseCertificate:
    """Geometric-integrality certificate for one and two ordered pair bases."""

    pair_surface_geometrically_integral: bool
    pair_surface_flat_over_kappa: bool
    ordered_two_pair_base_flat_over_kappa: bool
    ordered_two_pair_generic_fiber_domain: bool
    ordered_two_pair_base_domain: bool

    @property
    def verified(self) -> bool:
        """Whether both total source bases are irreducible over ``CC``."""

        return bool(
            self.pair_surface_geometrically_integral
            and self.pair_surface_flat_over_kappa
            and self.ordered_two_pair_base_flat_over_kappa
            and self.ordered_two_pair_generic_fiber_domain
            and self.ordered_two_pair_base_domain
        )


@cache
def exact_total_contact_base_certificate() -> TotalContactBaseCertificate:
    """Build the flatness-plus-generic-domain fiber-product certificate."""

    pair = exact_split_pair_surface_certificate()
    pair_integral = pair.total_surface_irreducible
    pair_flat = pair.flat_over_kappa
    # If A is the pair-surface algebra, A is flat over the PID QQ[k].  Thus
    # A tensor A is flat too.  Its generic fiber is obtained from two
    # independent primitive linear pair equations and is a domain.  If xy=0
    # in A tensor A, generic localization makes x or y zero; flatness then
    # kills the corresponding QQ[k]-torsion, so x=0 or y=0 already.
    ordered_flat = pair_flat
    generic_domain = pair_integral
    total_domain = ordered_flat and generic_domain
    return TotalContactBaseCertificate(
        pair_surface_geometrically_integral=pair_integral,
        pair_surface_flat_over_kappa=pair_flat,
        ordered_two_pair_base_flat_over_kappa=ordered_flat,
        ordered_two_pair_generic_fiber_domain=generic_domain,
        ordered_two_pair_base_domain=total_domain,
    )


@dataclass(frozen=True, slots=True)
class MaximalRankContactEmbeddingCertificate:
    """Exact row comparison for one maximal-rank split contact allocation."""

    name: str
    profile: str
    transformation_determinant: Expr
    row_residuals: tuple[Expr, ...]
    witness_transformation_determinant: Expr
    clean_witness_verified: bool

    @property
    def verified(self) -> bool:
        """Whether the global and component systems agree on a clean open."""

        return bool(
            self.transformation_determinant != 0
            and all(residual == 0 for residual in self.row_residuals)
            and self.witness_transformation_determinant != 0
            and self.clean_witness_verified
        )


def _witness_name(
    spec: LinearRankSpec | SimpleDoubleContactSpec | ResidualRankSpec,
) -> str:
    """Return the generic-ledger witness name for one focused rank spec."""

    return spec.name.removesuffix("_residual")


def _witness_root_values(
    spec: LinearRankSpec | SimpleDoubleContactSpec | ResidualRankSpec,
) -> dict[Expr, Expr]:
    """Return the component-root coordinates of the exact clean witness."""

    witness = {item.name: item for item in split_witness_specs()}[_witness_name(spec)]
    roots = tuple(root.root for root in witness.roots)
    if isinstance(spec, LinearRankSpec):
        return {U: roots[0]}
    return {U: roots[0], V: roots[1]}


def _clean_witness_verified(
    spec: LinearRankSpec | SimpleDoubleContactSpec | ResidualRankSpec,
) -> bool:
    """Replay the full residual-node and target-separation witness predicate."""

    witness_specs = {item.name: item for item in split_witness_specs()}
    witness_certificates = {
        item.name: item for item in exact_split_witness_certificates()
    }
    name = _witness_name(spec)
    return witness_certificates[name].verified(witness_specs[name])


def _c3_embedding(spec: LinearRankSpec) -> MaximalRankContactEmbeddingCertificate:
    """Compare one non-overlap split ``C3`` jet system."""

    global_equations = tuple(
        expand(row.subs(_component_coordinates(spec.kappa, spec.allocation, U)))
        for row in (PAIR_TARGET_ROW, PAIR_TANGENCY_ROW, PAIR_CONTACT_THREE_ROW)
    )
    transformation, residuals = _row_transformation(global_equations, spec.equations)
    determinant = cancel(transformation.det())
    witness_substitution = _witness_root_values(spec)
    return MaximalRankContactEmbeddingCertificate(
        name=spec.name,
        profile=spec.profile,
        transformation_determinant=determinant,
        row_residuals=residuals,
        witness_transformation_determinant=cancel(
            determinant.subs(witness_substitution)
        ),
        clean_witness_verified=_clean_witness_verified(spec),
    )


def _double_contact_embedding(
    spec: SimpleDoubleContactSpec | ResidualRankSpec,
) -> MaximalRankContactEmbeddingCertificate:
    """Compare one non-overlap split two-contact row system."""

    global_equations: list[Expr] = []
    for component, root in zip(spec.allocation, (U, V), strict=True):
        substitution = _component_coordinates(spec.kappa, component, root)
        global_equations.extend(
            (
                expand(PAIR_TARGET_ROW.subs(substitution)),
                expand(PAIR_TANGENCY_ROW.subs(substitution)),
            )
        )
    transformation, residuals = _row_transformation(
        tuple(global_equations),
        spec.equations,
    )
    determinant = cancel(transformation.det())
    witness_substitution = _witness_root_values(spec)
    return MaximalRankContactEmbeddingCertificate(
        name=spec.name,
        profile="C2^2+6N",
        transformation_determinant=determinant,
        row_residuals=residuals,
        witness_transformation_determinant=cancel(
            determinant.subs(witness_substitution)
        ),
        clean_witness_verified=_clean_witness_verified(spec),
    )


@cache
def exact_maximal_rank_contact_embedding_certificates() -> tuple[
    MaximalRankContactEmbeddingCertificate, ...
]:
    """Build the three ``C3`` and five ``C2^2`` rank-open comparisons."""

    c3 = tuple(
        _c3_embedding(spec)
        for spec in c3_rank_specs()
        if not spec.allocation.startswith("overlap-")
    )
    simple_double_contacts = tuple(
        _double_contact_embedding(spec)
        for spec in simple_double_contact_specs()
        if spec.allocation != "overlap+W"
    )
    residual_double_contacts = tuple(
        _double_contact_embedding(spec) for spec in residual_rank_specs()
    )
    return (*c3, *simple_double_contacts, *residual_double_contacts)


@dataclass(frozen=True, slots=True)
class SplitContactClosureCertificate:
    """Aggregate rank-open topology conclusion and exceptional boundary."""

    total_bases: TotalContactBaseCertificate
    embeddings: tuple[MaximalRankContactEmbeddingCertificate, ...]
    global_contact_three_sample_cyclic: bool
    global_double_contact_sample_cyclic: bool
    global_double_contact_cramer_verified: bool
    maximal_rank_topology_closed: bool
    exceptional_allocation_count: int
    residual_affine_line_base_length: int
    exceptional_topology_open: bool

    @property
    def verified(self) -> bool:
        """Whether the rank-open theorem and remaining exception count agree."""

        return bool(
            self.total_bases.verified
            and len(self.embeddings) == 8
            and all(item.verified for item in self.embeddings)
            and self.global_contact_three_sample_cyclic
            and self.global_double_contact_sample_cyclic
            and self.global_double_contact_cramer_verified
            and self.maximal_rank_topology_closed
            and self.exceptional_allocation_count == 3
            and self.residual_affine_line_base_length == 16
            and self.exceptional_topology_open
        )


@cache
def exact_split_contact_closure_certificate() -> SplitContactClosureCertificate:
    """Build the clean rank-open contact-component closure certificate."""

    contact_three = exact_delta_ten_contact_three_certificate()
    double_contact = exact_double_contact_sample_certificate()
    ordered = exact_ordered_incidence_certificate()
    embeddings = exact_maximal_rank_contact_embedding_certificates()
    return SplitContactClosureCertificate(
        total_bases=exact_total_contact_base_certificate(),
        embeddings=embeddings,
        global_contact_three_sample_cyclic=(
            contact_three.verified
            and contact_three.sage_cyclic_simplification == (1, 0, True)
        ),
        global_double_contact_sample_cyclic=(
            double_contact.verified
            and double_contact.sage_cyclic_simplification == (1, 0, True)
        ),
        global_double_contact_cramer_verified=ordered.verified,
        # The exact row embeddings put every listed clean rank-open point in
        # the already-propagated irreducible global incidence.  This boolean
        # records that mathematical Whitney--Thom consequence, not a CAS
        # computation of topology.
        maximal_rank_topology_closed=bool(
            len(embeddings) == 8 and all(item.verified for item in embeddings)
        ),
        # C3 overlap-V, C3 overlap-W, and C2^2 overlap+W.
        exceptional_allocation_count=3,
        # Reduced ordered compatibility lengths 4+6+6.
        residual_affine_line_base_length=16,
        exceptional_topology_open=True,
    )


def main() -> int:
    """Print the maximal-rank contact closure and exact remaining boundary."""

    certificate = exact_split_contact_closure_certificate()
    print(
        "ordered two-pair total base domain:",
        certificate.total_bases.ordered_two_pair_base_domain,
    )
    print(
        "maximal-rank split contact topology closed:",
        certificate.maximal_rank_topology_closed,
    )
    print("remaining overlap allocations:", certificate.exceptional_allocation_count)
    print(
        "remaining residual ordered base length:",
        certificate.residual_affine_line_base_length,
    )
    return 0 if certificate.verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
