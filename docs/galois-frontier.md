# What Galois theory does—and does not—solve

## Verdict

The plane Jacobian conjecture is not solved by Galois theory.  The classical
result is conditional:

> If a characteristic-zero Keller map induces a Galois extension of rational
> function fields, then the map is a polynomial automorphism.

For a plane map `F=(P,Q)`, the relevant extension is

\[
  k(P,Q)\subset k(x,y).
\]

The nonzero constant Jacobian makes this a finite separable **field**
extension, so the associated polynomial morphism is dominant and generically
finite.  It need not be a finite morphism: that stronger statement would
require `k[x,y]` to be a finite `k[P,Q]`-module (equivalently, integral over
`k[P,Q]`) and would imply properness.  The Keller condition supplies none of
those global finiteness conclusions.  It also does not make the function-field
extension normal.  “Finite and separable” gives a finite Galois closure of the
function field; it does not say that the original extension is itself Galois
or by itself specify a global finite normal-cover model.  That missing
normality hypothesis is exactly why the classical theorem does not settle
`JC(2)`.

Campbell proved the complex Galois case, and Razar and Wright supplied
algebraic treatments.  The result is also recorded among the equivalent
conditions in the Bass–Connell–Wright reduction paper.

## The low-sheet frontier

Here “degree” means the generic number of points in a fiber, equivalently
`[k(x,y):k(P,Q)]`.  It is not the ordinary total degree of `P` or `Q`.

| Generic fiber degree | Status for a hypothetical plane Keller map |
|---:|---|
| 1 | Keller's birational theorem gives an automorphism |
| 2 | the extension is automatically Galois, so it is an automorphism |
| 3 | excluded by Orevkov's two-/three-sheet theorem |
| 4 | excluded in full generality by Domrina |
| 5 | included in Żołądek's theorem through degree 5 |
| 6 | first unresolved generic-fiber degree in the accepted literature |

Thus the “next dimensional case” is not dimension four.  The announced map
disproves `JC(n)` for every `n >= 3`, while `JC(2)` remains the sole surviving
fixed-dimensional case.  Inside `JC(2)`, generic degree six is the sharpest
small-sheet frontier located in this audit.

Alexander Borisov's framework paper independently states the lower bound six
and asks whether its birational-combinatorial methods can improve it.  That is
a real open problem, not something this repository claims to close.

The repository does now obtain a new necessary condition inside degree six.
Affine fixed-sheet inertia, Orevkov's defect identity, finite-normalization
deck symmetry, purity, and local monodromy eliminate every imprimitive
transitive action.  Thus a six-sheet counterexample would have to have one of
the four primitive actions `A5`, `S5`, `A6`, or `S6`.  This is a substantial
filter, but it does not exclude those four actions and therefore does not raise
the accepted lower bound.  See
[`six-sheet-monodromy.md`](six-sheet-monodromy.md) for the complete proof.

## Two tempting but unsafe shortcuts

### A claimed plane proof

Wolfgang Bartenwerfer's 2022 paper claimed to prove the plane conjecture.  The
journal's 2023 correction says that the main proof is wrong and describes the
correction as amounting to withdrawal of the whole paper.  It cannot be used
as evidence that `JC(2)` is solved.

### A claimed prime-degree theorem

A 2024 arXiv v1 preprint claims that no plane Keller map can have prime
function-field degree.  This audit found no peer-reviewed confirmation, and
the cited online discussion does not itself state the needed classification
theorem.  We therefore do not use the preprint as an accepted dependency.
This caution does not change the first open degree: six is composite.

## A known quadratic-in-one-variable class

Sabatini proves a positive theorem for real planar maps whose coordinates are
at most quadratic in one chosen variable.  In the polynomial
constant-Jacobian case, the map factors into triangular polynomial maps and
therefore has a polynomial inverse.

The Lean theorem in this repository is an independently proved algebraic
normal-form statement over an arbitrary characteristic-zero field.  It should
not be attributed verbatim to the real-analytic theorem: the relationship is
that Sabatini establishes prior mathematical coverage of the bounded-degree
class, while the repository supplies its own explicit field-algebra proof.

## What this repository attacks

There are four intentionally separate levels.

1. `JC(2)` itself, beginning at generic degree six, remains open.
2. At generic degree six, the repository proves the primitive-monodromy
   restriction above and sharper conditional one-dicritical passports.  These
   are necessary conditions, not existence or nonexistence theorems for all
   six-sheet maps.
3. The repository proves kernel-checked positive fragments for maps of bounded
   degree in one variable.  The quadratic-in-one-variable normal form is a
   Lean formalization and explicit reconstruction of a mathematically known
   special case; no novelty claim is made.
4. For the new three-dimensional counterexample, the repository determines
   the full fiber stratification, exact omitted curve, and complete
   nonproper-value set.  These are finite consequences of the explicit map and
   do not purport to solve `JC(2)`.

This is the honest boundary: Galois theory explains why degree two is safe and
why the three-dimensional example must be nonnormal, but the first unknown
plane degree requires geometry at infinity beyond the Galois-case theorem.

## Primary sources

- L. Andrew Campbell,
  [“A condition for a polynomial map to be invertible”](https://doi.org/10.1007/BF01349234),
  *Mathematische Annalen* 205 (1973), 243–248.
- Michael Razar,
  [“Polynomial maps with constant Jacobian”](https://doi.org/10.1007/BF02764906),
  *Israel Journal of Mathematics* 32 (1979), 97–106.
- David Wright,
  [“On the Jacobian conjecture”](https://projecteuclid.org/journals/illinois-journal-of-mathematics/volume-25/issue-3/On-the-Jacobian-conjecture/10.1215/ijm/1256047158.full),
  *Illinois Journal of Mathematics* 25 (1981), 423–440.
- Hyman Bass, Edwin Connell, and David Wright,
  [“The Jacobian conjecture: reduction of degree and formal expansion of the inverse”](https://doi.org/10.1090/S0273-0979-1982-15032-7),
  *Bulletin of the AMS* 7 (1982), 287–330.
- S. Yu. Orevkov,
  [“On three-sheeted polynomial mappings of `C^2`”](https://doi.org/10.1070/IM1987v029n03ABEH000984),
  *Mathematics of the USSR-Izvestiya* 29 (1987), 587–596.
- A. V. Domrina,
  [“On four-sheeted polynomial mappings of `C^2`. II. The general case”](https://doi.org/10.1070/im2000v064n01ABEH000273),
  *Izvestiya: Mathematics* 64 (2000), 39–76.
- Henryk Żołądek,
  [“An application of Newton–Puiseux charts to the Jacobian problem”](https://doi.org/10.1016/j.top.2008.04.001),
  *Topology* 47 (2008), 431–469.
- Alexander Borisov,
  [“Frameworks for two-dimensional Keller maps”](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v27i3p54),
  *Electronic Journal of Combinatorics* 27 (2020), P3.54.
- Marco Sabatini,
  [“Global injectivity of planar non-singular maps that are polynomial in one variable”](https://www.impan.pl/en/publishing-house/journals-and-series/colloquium-mathematicum/all/175/1/115518/global-injectivity-of-planar-non-singular-maps-that-are-polynomial-in-one-variable),
  *Colloquium Mathematicum* 175 (2024), 137–151.
- Ralph Chill and Gabriele Nebe,
  [“Correction to: The Cremona problem in dimension 2”](https://doi.org/10.1007/s00013-023-01863-0),
  *Archiv der Mathematik* 120 (2023), 665–666.
- Vered Moskowicz,
  [“There are no Keller maps having prime degree field extensions”](https://arxiv.org/abs/2407.13795),
  arXiv v1 (2024).  This is listed as an audited claim, not as an accepted
  dependency.
- Elia Bisi, Piotr Dyszewski, Nina Gantert, Samuel G. G. Johnston, Joscha
  Prochno, and Dominik Schmid,
  [“Random planar trees and the Jacobian conjecture”](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/jlms.70416),
  *Journal of the London Mathematical Society* 113 (2026), e70416, used here
  only as a current peer-reviewed status reference.
- David Rodríguez Díaz,
  [“On the origin of the Jacobian conjecture”](https://comptes-rendus.academie-sciences.fr/mathematique/articles/10.5802/crmath.831/),
  *Comptes Rendus Mathématique* 364 (2026), 363–370, used as a current account
  of the continuing geometry-at-infinity obstruction.
