/-
Copyright (c) 2026 Alok Singh. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Alok Singh
-/
import JacobianTwo.AffineCoordinate
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

end JacobianTwo.VariableLeadingQuadratic
