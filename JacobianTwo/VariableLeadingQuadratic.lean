/-
Copyright (c) 2026 Alok Singh. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Alok Singh
-/
import JacobianTwo.AffineCoordinate
import Mathlib.RingTheory.UniqueFactorizationDomain.Multiplicity
import Mathlib.Tactic.Ring

/-!
# Keller maps with a variable-leading quadratic coordinate

This module begins the coefficient descent for a plane polynomial map whose
second coordinate has the form

`a(x) * y^2 + g(x) * y + f(x)`.

The coefficient immediately above the `y`-degree of the first coordinate is
independent of `g` and `f`.  Its vanishing gives the differential equation

`2 * a * p_n' - n * a' * p_n = 0`

for the leading `y`-coefficient `p_n` of the first coordinate.
-/

noncomputable section

open scoped Polynomial Polynomial.Bivariate

namespace JacobianTwo.VariableLeadingQuadratic

open Polynomial
open JacobianTwo.AffineCoordinate

variable {K : Type*} [Field K]

@[simp]
theorem C_natCast_inner (n : ℕ) : C (n : K[X]) = (n : K[X][Y]) := by
  rw [map_natCast]

@[simp]
theorem C_ofNat_inner (n : ℕ) [n.AtLeastTwo] :
    C (ofNat(n) : K[X]) = (ofNat(n) : K[X][Y]) := by
  rw [C_ofNat]

@[simp]
theorem CC_ofNat_scalar (n : ℕ) [n.AtLeastTwo] :
    C (C (ofNat(n) : K)) = (ofNat(n) : K[X][Y]) := by
  rw [C_ofNat, C_ofNat]

/-- A polynomial that is quadratic in `y`, allowing its leading coefficient
to vary with `x`. -/
def variableQuadraticCoordinate (a g f : K[X]) : K[X][Y] :=
  C a * Y ^ 2 + C g * Y + C f

/-- The exact UFD square-shape conclusion needed after the odd-degree
leading-coefficient identity.  Establishing this predicate from that identity
is a separate factorization step. -/
def IsNonzeroScalarTimesSquare (a : K[X]) : Prop :=
  ∃ eps : K, ∃ h : K[X], eps ≠ 0 ∧ a = C eps * h ^ 2

/-- In `K[x]`, if a square is a nonzero scalar times an odd power, then the
base is a nonzero scalar times a square.  The same square root also describes
the original square root up to a nonzero scalar. -/
theorem scalar_times_square_and_odd_root_of_sq_eq_scalar_mul_pow
    {p a : K[X]} {c : K} {n : ℕ}
    (hc : c ≠ 0) (ha : a ≠ 0) (hodd : Odd n)
    (hsq : p ^ 2 = C c * a ^ n) :
    ∃ eps : K, ∃ h : K[X], ∃ lam : K,
      eps ≠ 0 ∧ lam ≠ 0 ∧
      a = C eps * h ^ 2 ∧ p = C lam * h ^ n := by
  rcases hodd with ⟨m, hn⟩
  have hpowDvd : (a ^ m) ^ 2 ∣ p ^ 2 := by
    rw [hsq]
    refine ⟨C c * a, ?_⟩
    rw [hn, show 2 * m + 1 = m + m + 1 by omega, pow_succ, pow_add,
      pow_two]
    ring
  have hdvd : a ^ m ∣ p :=
    (UniqueFactorizationMonoid.pow_dvd_pow_iff_dvd (by omega : 2 ≠ 0)).mp
      hpowDvd
  obtain ⟨h, hpresentation⟩ := hdvd
  have hpresentation' : p = a ^ m * h := hpresentation
  have hcancel :
      a ^ (m + m) * h ^ 2 = a ^ (m + m) * (C c * a) := by
    calc
      a ^ (m + m) * h ^ 2 = (a ^ m * h) ^ 2 := by
        rw [pow_add, pow_two, pow_two]
        ring
      _ = p ^ 2 := by rw [← hpresentation']
      _ = C c * a ^ n := hsq
      _ = a ^ (m + m) * (C c * a) := by
        rw [hn, pow_succ, pow_add]
        ring
  have hahpow : a ^ (m + m) ≠ 0 := pow_ne_zero _ ha
  have hhsq : h ^ 2 = C c * a := mul_left_cancel₀ hahpow hcancel
  let eps : K := c⁻¹
  let lam : K := c⁻¹ ^ m
  have heps : eps ≠ 0 := by
    dsimp [eps]
    exact inv_ne_zero hc
  have hlam : lam ≠ 0 := by
    dsimp [lam]
    exact pow_ne_zero _ (inv_ne_zero hc)
  have haSquare : a = C eps * h ^ 2 := by
    rw [hhsq]
    symm
    calc
      C eps * (C c * a) = (C eps * C c) * a := by ring
      _ = C (eps * c) * a := by rw [C_mul]
      _ = a := by
        dsimp [eps]
        rw [inv_mul_cancel₀ hc, C_1, one_mul]
  have hpRoot : p = C lam * h ^ n := by
    rw [hpresentation', haSquare, hn]
    dsimp [lam, eps]
    rw [mul_pow, ← C_pow, pow_succ, pow_add]
    ring
  exact ⟨eps, h, lam, heps, hlam, haSquare, hpRoot⟩

/-- The square-shape projection of
`scalar_times_square_and_odd_root_of_sq_eq_scalar_mul_pow`. -/
theorem isNonzeroScalarTimesSquare_of_sq_eq_scalar_mul_odd_pow
    {p a : K[X]} {c : K} {n : ℕ}
    (hc : c ≠ 0) (ha : a ≠ 0) (hodd : Odd n)
    (hsq : p ^ 2 = C c * a ^ n) :
    IsNonzeroScalarTimesSquare a := by
  obtain ⟨eps, h, _lam, heps, _hlam, haSquare, _hpRoot⟩ :=
    scalar_times_square_and_odd_root_of_sq_eq_scalar_mul_pow hc ha hodd hsq
  exact ⟨eps, h, heps, haSquare⟩

@[simp]
theorem yDerivative_variableQuadraticCoordinate (a g f : K[X]) :
    yDerivative (variableQuadraticCoordinate a g f) =
      C (C (2 : K) * a) * Y + C g := by
  simp [variableQuadraticCoordinate, yDerivative_apply, derivative_pow]
  ring

@[simp]
theorem xDerivative_variableQuadraticCoordinate (a g f : K[X]) :
    xDerivative (variableQuadraticCoordinate a g f) =
      C (derivative a) * Y ^ 2 + C (derivative g) * Y + C (derivative f) := by
  simp [variableQuadraticCoordinate, Derivation.leibniz]

theorem natDegree_variableQuadraticCoordinate {a g f : K[X]} (ha : a ≠ 0) :
    (variableQuadraticCoordinate a g f).natDegree = 2 := by
  exact natDegree_quadratic ha

theorem leadingCoeff_variableQuadraticCoordinate {a g f : K[X]} (ha : a ≠ 0) :
    (variableQuadraticCoordinate a g f).leadingCoeff = a := by
  exact leadingCoeff_quadratic ha

/-- The coefficient of the Jacobian one step above the `y`-degree of `P`.
Only the leading coefficient of `P` and the quadratic coefficient `a` occur. -/
theorem coeff_natDegree_succ_jacobian_variableQuadraticCoordinate
    (P : K[X][Y]) (a g f : K[X]) :
    (jacobian P (variableQuadraticCoordinate a g f)).coeff
        (P.natDegree + 1) =
      C (2 : K) * derivative P.leadingCoeff * a -
        C (P.natDegree : K) * P.leadingCoeff * derivative a := by
  by_cases hn : P.natDegree = 0
  · rw [hn]
    simp only [Nat.zero_add, Nat.cast_zero, C_0, zero_mul]
    rw [jacobian, xDerivative_variableQuadraticCoordinate,
      yDerivative_variableQuadraticCoordinate, coeff_sub]
    simp only [mul_add, coeff_add, ← mul_assoc, coeff_mul_X, coeff_mul_C,
      coeff_xDerivative, yDerivative_apply, coeff_derivative]
    have hP1 : P.coeff 1 = 0 := by
      exact coeff_eq_zero_of_natDegree_lt (by omega)
    have hP2 : P.coeff 2 = 0 := by
      exact coeff_eq_zero_of_natDegree_lt (by omega)
    rw [hP1, hP2]
    simp only [map_zero, zero_mul, Nat.cast_one, add_zero]
    rw [show P.coeff 0 = P.leadingCoeff by simpa [hn] using P.coeff_natDegree]
    have hquadraticTerm :
        (derivative P * C (derivative a) * Y ^ 2).coeff 1 = 0 := by
      rw [pow_two, ← mul_assoc, coeff_mul_X, coeff_mul_X_zero]
    rw [hquadraticTerm]
    ring
  · obtain ⟨n, hnP⟩ := Nat.exists_eq_succ_of_ne_zero hn
    rw [hnP]
    rw [jacobian, xDerivative_variableQuadraticCoordinate,
      yDerivative_variableQuadraticCoordinate, coeff_sub]
    simp only [mul_add, coeff_add, ← mul_assoc, pow_two, coeff_mul_X,
      coeff_mul_C, coeff_xDerivative, yDerivative_apply, coeff_derivative]
    have hlead : P.coeff (n + 1) = P.leadingCoeff := by
      simpa [hnP] using P.coeff_natDegree
    have hnext : P.coeff (n + 1 + 1) = 0 := by
      simpa [hnP] using P.coeff_natDegree_succ_eq_zero
    have hnext2 : P.coeff (n + 1 + 1 + 1) = 0 := by
      exact coeff_eq_zero_of_natDegree_lt (by omega)
    rw [hlead, hnext, hnext2]
    simp only [map_zero, zero_mul, add_zero]
    push_cast
    rw [map_add]
    simp only [map_one, map_natCast]
    ring

/-- A constant Jacobian makes the top coefficient identity vanish. -/
theorem leadingCoeff_differential_eq_of_constant_jacobian
    (P : K[X][Y]) (a g f : K[X]) (k : K)
    (hjac : jacobian P (variableQuadraticCoordinate a g f) = CC k) :
    C (2 : K) * derivative P.leadingCoeff * a =
      C (P.natDegree : K) * P.leadingCoeff * derivative a := by
  have hcoeff := congrArg
    (fun R : K[X][Y] => R.coeff (P.natDegree + 1)) hjac
  rw [coeff_natDegree_succ_jacobian_variableQuadraticCoordinate] at hcoeff
  have hrhs : (CC k).coeff (P.natDegree + 1) = 0 := by
    simp
  exact sub_eq_zero.mp (hcoeff.trans hrhs)

/-- If the `y`-degree of `P` is a positive even number `2*r`, its leading
coefficient is a scalar multiple of `a^r`.  This is the algebraic input for
subtracting a scalar multiple of `Q^r` in the even-degree descent. -/
theorem leadingCoeff_eq_C_mul_pow_of_even_natDegree [CharZero K]
    {P : K[X][Y]} {a g f : K[X]} {k : K} {r : ℕ}
    (ha : a ≠ 0) (hr : r ≠ 0)
    (hdegree : P.natDegree = 2 * r)
    (hjac : jacobian P (variableQuadraticCoordinate a g f) = CC k) :
    ∃ lam : K, P.leadingCoeff = C lam * a ^ r := by
  have htop :=
    leadingCoeff_differential_eq_of_constant_jacobian P a g f k hjac
  rw [hdegree] at htop
  have htwo : (2 : K) ≠ 0 := by norm_num
  have htwoC : C (2 : K) ≠ (0 : K[X]) := C_ne_zero.mpr htwo
  have hscaled :
      C (2 : K) * (derivative P.leadingCoeff * a) =
        C (2 : K) *
          (C (r : K) * P.leadingCoeff * derivative a) := by
    calc
      C (2 : K) * (derivative P.leadingCoeff * a) =
          C ((2 * r : ℕ) : K) * P.leadingCoeff * derivative a := by
            simpa only [mul_assoc] using htop
      _ = C (2 : K) *
          (C (r : K) * P.leadingCoeff * derivative a) := by
            push_cast
            rw [map_mul]
            ring
  have hreduced : derivative P.leadingCoeff * a =
      C (r : K) * P.leadingCoeff * derivative a := by
    exact mul_left_cancel₀ htwoC hscaled
  exact eq_C_mul_pow_of_differential_eq ha hr hreduced

/-- At odd `y`-degree, the square of the leading coefficient is a nonzero
scalar multiple of the corresponding power of the quadratic coefficient.
This is the differential identity whose UFD factorization should imply
`IsNonzeroScalarTimesSquare a`. -/
theorem leadingCoeff_sq_eq_C_mul_pow_of_odd_natDegree [CharZero K]
    {P : K[X][Y]} {a g f : K[X]} {k : K}
    (ha : a ≠ 0) (hodd : Odd P.natDegree)
    (hjac : jacobian P (variableQuadraticCoordinate a g f) = CC k) :
    ∃ c : K, c ≠ 0 ∧
      P.leadingCoeff ^ 2 = C c * a ^ P.natDegree := by
  have hn : P.natDegree ≠ 0 := Nat.ne_of_gt hodd.pos
  have htop :=
    leadingCoeff_differential_eq_of_constant_jacobian P a g f k hjac
  have hsquareDifferential :
      derivative (P.leadingCoeff ^ 2) * a =
        C (P.natDegree : K) * P.leadingCoeff ^ 2 * derivative a := by
    rw [derivative_pow]
    simp only [Nat.reduceSubDiff, pow_one]
    calc
      C (2 : K) * P.leadingCoeff * derivative P.leadingCoeff * a =
          P.leadingCoeff *
            (C (2 : K) * derivative P.leadingCoeff * a) := by ring
      _ = P.leadingCoeff *
          (C (P.natDegree : K) * P.leadingCoeff * derivative a) := by
            rw [htop]
      _ = C (P.natDegree : K) * P.leadingCoeff ^ 2 * derivative a := by
            ring
  obtain ⟨c, hc⟩ :=
    eq_C_mul_pow_of_differential_eq ha hn hsquareDifferential
  have hPne : P ≠ 0 := by
    intro hzero
    rw [hzero, natDegree_zero] at hn
    exact hn rfl
  have hlead : P.leadingCoeff ≠ 0 := leadingCoeff_ne_zero.mpr hPne
  have hcne : c ≠ 0 := by
    intro hzero
    have hsquareZero : P.leadingCoeff ^ 2 = 0 := by
      simpa [hzero] using hc
    exact (pow_ne_zero 2 hlead) hsquareZero
  exact ⟨c, hcne, hc⟩

/-- The UFD consequence of the odd-degree leading-coefficient identity.  The
quadratic coefficient is a nonzero scalar times a square, and the leading
coefficient of `P` is a nonzero scalar times the matching odd power of the
same square root. -/
theorem odd_natDegree_quadraticCoeff_and_leadingCoeff_shape [CharZero K]
    {P : K[X][Y]} {a g f : K[X]} {k : K}
    (ha : a ≠ 0) (hodd : Odd P.natDegree)
    (hjac : jacobian P (variableQuadraticCoordinate a g f) = CC k) :
    ∃ eps : K, ∃ h : K[X], ∃ lam : K,
      eps ≠ 0 ∧ lam ≠ 0 ∧
      a = C eps * h ^ 2 ∧
      P.leadingCoeff = C lam * h ^ P.natDegree := by
  obtain ⟨c, hc, hsq⟩ :=
    leadingCoeff_sq_eq_C_mul_pow_of_odd_natDegree ha hodd hjac
  exact scalar_times_square_and_odd_root_of_sq_eq_scalar_mul_pow
    hc ha hodd hsq

/-- Repeated target shears remove every positive even `y`-degree from the
first coordinate.  A nonzero constant Jacobian excludes degree zero, so the
terminal residual has odd `y`-degree. -/
theorem exists_target_shear_add_odd_residual [CharZero K]
    {P : K[X][Y]} {a g f : K[X]} {k : K}
    (ha : a ≠ 0) (hk : k ≠ 0)
    (hjac : jacobian P (variableQuadraticCoordinate a g f) = CC k) :
    ∃ H : K[X], ∃ R : K[X][Y],
      P = (H.map C).comp (variableQuadraticCoordinate a g f) + R ∧
      jacobian R (variableQuadraticCoordinate a g f) = CC k ∧
      Odd R.natDegree := by
  let Q := variableQuadraticCoordinate a g f
  have descend : ∀ n : ℕ, ∀ R : K[X][Y], R.natDegree = n →
      jacobian R Q = CC k →
      ∃ H : K[X], ∃ S : K[X][Y],
        R = (H.map C).comp Q + S ∧ jacobian S Q = CC k ∧
        Odd S.natDegree := by
    intro n
    induction n using Nat.strong_induction_on with
    | h n ih =>
        intro R hdegree hRjac
        by_cases hn : n = 0
        · have hRconstant : R = C (R.coeff 0) := by
            apply eq_C_of_natDegree_eq_zero
            rw [hdegree, hn]
          have htop :=
            leadingCoeff_differential_eq_of_constant_jacobian R a g f k (by
              simpa only [Q] using hRjac)
          have htopzero :
              C (2 : K) * derivative R.leadingCoeff * a = 0 := by
            simpa [hdegree, hn] using htop
          have htwo : (2 : K) ≠ 0 := by norm_num
          have htwoC : C (2 : K) ≠ (0 : K[X]) := C_ne_zero.mpr htwo
          have hderivative : derivative R.leadingCoeff = 0 := by
            rcases mul_eq_zero.mp htopzero with hleft | haright
            · exact (mul_eq_zero.mp hleft).resolve_left htwoC
            · exact (ha haright).elim
          have hcoeff : R.coeff 0 = R.leadingCoeff := by
            simpa [hdegree, hn] using R.coeff_natDegree
          have hRjaczero : jacobian R Q = 0 := by
            rw [hRconstant]
            have hcoeffDerivative : derivative (R.coeff 0) = 0 := by
              rw [hcoeff]
              exact hderivative
            simp [jacobian, hcoeffDerivative]
          have hkCC : CC k ≠ (0 : K[X][Y]) :=
            C_ne_zero.mpr (C_ne_zero.mpr hk)
          exact (hkCC (hRjac.symm.trans hRjaczero)).elim
        · obtain heven | hodd := Nat.even_or_odd n
          · rcases heven with ⟨r, hneven⟩
            have hr : r ≠ 0 := by
              intro hrzero
              subst r
              simp at hneven
              exact hn hneven
            have hdegreeEven : R.natDegree = 2 * r := by
              rw [hdegree, hneven]
              omega
            obtain ⟨lam, hlamodel⟩ :=
              leadingCoeff_eq_C_mul_pow_of_even_natDegree ha hr hdegreeEven (by
                simpa only [Q] using hRjac)
            have hRne : R ≠ 0 := by
              intro hzero
              rw [hzero, natDegree_zero] at hdegree
              exact hn hdegree.symm
            have hlam : lam ≠ 0 := by
              intro hlamzero
              have hleadzero : R.leadingCoeff = 0 := by
                simp [hlamodel, hlamzero]
              exact (leadingCoeff_ne_zero.mpr hRne) hleadzero
            let model : K[X][Y] := CC lam * Q ^ r
            let tail : K[X][Y] := R - model
            have hQdegree : Q.natDegree = 2 := by
              exact natDegree_variableQuadraticCoordinate ha
            have hQne : Q ≠ 0 := by
              intro hzero
              rw [hzero, natDegree_zero] at hQdegree
              omega
            have hmodel_ne : model ≠ 0 := by
              dsimp [model]
              exact mul_ne_zero (C_ne_zero.mpr (C_ne_zero.mpr hlam))
                (pow_ne_zero _ hQne)
            have hmodel_natDegree : model.natDegree = R.natDegree := by
              dsimp [model]
              rw [natDegree_mul (C_ne_zero.mpr (C_ne_zero.mpr hlam))
                  (pow_ne_zero _ hQne), natDegree_C, zero_add, natDegree_pow,
                hQdegree, hdegreeEven]
              omega
            have hmodel_leadingCoeff : model.leadingCoeff = R.leadingCoeff := by
              dsimp [model]
              rw [leadingCoeff_mul, leadingCoeff_C, leadingCoeff_pow,
                show Q.leadingCoeff = a by
                  exact leadingCoeff_variableQuadraticCoordinate ha]
              exact hlamodel.symm
            have hmodel_degree : model.degree = R.degree := by
              rw [degree_eq_natDegree hmodel_ne, degree_eq_natDegree hRne,
                hmodel_natDegree]
            have htail_degree : tail.degree < R.degree := by
              dsimp [tail]
              exact degree_sub_lt hmodel_degree.symm hRne
                hmodel_leadingCoeff.symm
            have htail_natDegree : tail.natDegree < n := by
              by_cases htail : tail = 0
              · simp [htail, Nat.pos_of_ne_zero hn]
              · apply (natDegree_lt_iff_degree_lt htail).mpr
                rw [degree_eq_natDegree hRne, hdegree] at htail_degree
                exact htail_degree
            have htail_jac : jacobian tail Q = CC k := by
              dsimp [tail, model]
              rw [jacobian_sub_left, jacobian_CC_mul_pow_left, sub_zero]
              exact hRjac
            obtain ⟨H, S, htail_repr, hSjac, hSodd⟩ :=
              ih tail.natDegree htail_natDegree tail rfl htail_jac
            refine ⟨H + monomial r lam, S, ?_, hSjac, hSodd⟩
            calc
              R = tail + model := by simp [tail]
              _ = ((H.map C).comp Q + S) + CC lam * Q ^ r := by
                rw [htail_repr]
              _ = (((H + monomial r lam).map C).comp Q) + S := by
                rw [Polynomial.map_add, add_comp, map_monomial, monomial_comp]
                ring
          · refine ⟨0, R, ?_, hRjac, ?_⟩
            · simp
            · rw [hdegree]
              exact hodd
  exact descend P.natDegree P rfl hjac

end JacobianTwo.VariableLeadingQuadratic
