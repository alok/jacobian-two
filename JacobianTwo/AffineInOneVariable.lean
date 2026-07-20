/-
Copyright (c) 2026 Alok Singh. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Alok Singh
-/
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.RingTheory.Polynomial.Wronskian
import Mathlib.Tactic.Ring

/-!
# A two-dimensional obstruction for maps affine in one variable

Consider a plane polynomial map

`G(x, y) = (a(x) * y + b(x), c(x) * y + d(x))`.

Its Jacobian determinant is affine in `y`; its two coefficients are

`a' * c - a * c'` and `b' * c - a * d'`.

This module proves that if the determinant is a nonzero constant, then both
`a'` and `c'` vanish. Over a characteristic-zero field, `a` and `c` are
therefore constants. Thus a two-dimensional counterexample cannot arise from
this direct analogue of the announced three-dimensional map's affine-in-`z`
shape.

## Main results

- `affineInY_slopes_derivative_eq_zero`: the `y`-slopes have zero derivative.
- `affineInY_slopes_constant`: in characteristic zero, those slopes are
  constant polynomials.
-/

noncomputable section

open scoped Polynomial

namespace JacobianTwo.AffineInOneVariable

open Polynomial

variable {K : Type*} [Field K]

/-- The Jacobian determinant of
`(a(x) * y + b(x), c(x) * y + d(x))`, represented as a polynomial in `y`
whose coefficients are polynomials in `x`. -/
def affineInYJacobian (a b c d : K[X]) : (K[X])[X] :=
  C (derivative a * c - a * derivative c) * X +
    C (derivative b * c - a * derivative d)

@[simp]
theorem coeff_one_affineInYJacobian (a b c d : K[X]) :
    (affineInYJacobian a b c d).coeff 1 = derivative a * c - a * derivative c := by
  simp [affineInYJacobian]

@[simp]
theorem coeff_zero_affineInYJacobian (a b c d : K[X]) :
    (affineInYJacobian a b c d).coeff 0 = derivative b * c - a * derivative d := by
  simp [affineInYJacobian]

/-- If an affine-in-`y` plane map has nonzero constant Jacobian determinant,
then the derivatives of both `y`-slope polynomials vanish. -/
theorem affineInY_slopes_derivative_eq_zero {a b c d : K[X]} {k : K}
    (hk : k ≠ 0)
    (hdet : affineInYJacobian a b c d = C (C k)) :
    derivative a = 0 ∧ derivative c = 0 := by
  have hy : derivative a * c - a * derivative c = 0 := by
    have h := congrArg (fun p : (K[X])[X] ↦ p.coeff 1) hdet
    simpa using h
  have hconst : derivative b * c - a * derivative d = C k := by
    have h := congrArg (fun p : (K[X])[X] ↦ p.coeff 0) hdet
    simpa using h
  have hac : IsCoprime a c := by
    refine ⟨-C (k⁻¹) * derivative d, C (k⁻¹) * derivative b, ?_⟩
    calc
      (-C (k⁻¹) * derivative d) * a + (C (k⁻¹) * derivative b) * c =
          C (k⁻¹) * (derivative b * c - a * derivative d) := by ring
      _ = C (k⁻¹) * C k := by rw [hconst]
      _ = 1 := by rw [← C_mul]; simp [hk]
  apply hac.wronskian_eq_zero_iff.mp
  calc
    wronskian a c = -(derivative a * c - a * derivative c) := by
      simp only [wronskian]
      ring
    _ = 0 := by rw [hy, neg_zero]

/-- Over a characteristic-zero field, the two `y`-slope polynomials in a
constant-Jacobian affine-in-`y` map are constant. -/
theorem affineInY_slopes_constant [CharZero K] {a b c d : K[X]} {k : K}
    (hk : k ≠ 0)
    (hdet : affineInYJacobian a b c d = C (C k)) :
    a = C (a.coeff 0) ∧ c = C (c.coeff 0) := by
  obtain ⟨ha, hc⟩ := affineInY_slopes_derivative_eq_zero hk hdet
  exact ⟨eq_C_of_derivative_eq_zero ha, eq_C_of_derivative_eq_zero hc⟩

end JacobianTwo.AffineInOneVariable
