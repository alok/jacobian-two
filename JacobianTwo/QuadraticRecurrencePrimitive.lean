/-
Copyright (c) 2026 Alok Singh. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Alok Singh
-/
import JacobianTwo.QuadraticRecurrence
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# A polynomial primitive for the quadratic recurrence

For a polynomial `S`, the downward quadratic recurrence asks for a polynomial
`T` satisfying

`T' = X * S' + (1/2) * S`.

In characteristic zero this file constructs `T` coefficientwise and proves
the identity.  The next descent layer evaluates this identity at the centered
residual `F`.
-/

noncomputable section

open Polynomial

namespace JacobianTwo.QuadraticRecurrencePrimitive

variable {K : Type*} [Field K] [CharZero K]

/-- The coefficientwise primitive of `X*S' + S/2`.  A monomial `a*X^n`
is sent to `((n+1/2)/(n+1))*a*X^(n+1)`. -/
def recurrencePrimitive (S : K[X]) : K[X] :=
  S.sum fun n a ↦
    monomial (n + 1)
      ((((n : K) + (1 / 2 : K)) / (n + 1 : K)) * a)

omit [CharZero K] in
theorem recurrencePrimitive_add (S T : K[X]) :
    recurrencePrimitive (S + T) =
      recurrencePrimitive S + recurrencePrimitive T := by
  rw [recurrencePrimitive, recurrencePrimitive, recurrencePrimitive]
  apply sum_add_index
  · intro n
    simp
  · intro n a b
    simp [mul_add]

omit [CharZero K] in
theorem recurrencePrimitive_monomial (n : ℕ) (a : K) :
    recurrencePrimitive (monomial n a) =
      monomial (n + 1)
        ((((n : K) + (1 / 2 : K)) / (n + 1 : K)) * a) := by
  rw [recurrencePrimitive]
  exact sum_monomial_index a _ (by simp)

/-- The constructed polynomial has the required derivative. -/
theorem derivative_recurrencePrimitive (S : K[X]) :
    derivative (recurrencePrimitive S) =
      X * derivative S + C (1 / 2 : K) * S := by
  induction S using Polynomial.induction_on' with
  | add P Q hP hQ =>
      rw [recurrencePrimitive_add, derivative_add, hP, hQ, derivative_add]
      ring
  | monomial n a =>
      cases n with
      | zero =>
          rw [recurrencePrimitive_monomial, derivative_monomial_succ,
            derivative_monomial]
          simp only [Nat.cast_zero, mul_zero, monomial_zero_right, mul_zero,
            zero_add, C_mul_monomial]
          ext m
          simp only [coeff_monomial]
          split_ifs <;> norm_num
      | succ n =>
          have hden : (n + 2 : K) ≠ 0 := by
            exact_mod_cast Nat.succ_ne_zero (n + 1)
          have hden' : (2 + n : K) ≠ 0 := by
            simpa [add_comm] using hden
          have hden'' : (n : K) + 1 + 1 ≠ 0 := by
            simpa only [show (2 : K) = 1 + 1 by norm_num, add_assoc]
              using hden
          rw [recurrencePrimitive_monomial, derivative_monomial_succ,
            derivative_monomial_succ]
          rw [X_mul_monomial, C_mul_monomial]
          ext m
          simp only [coeff_add, coeff_monomial]
          split_ifs
          · push_cast
            field_simp [hden, hden', hden'']
          · ring

end JacobianTwo.QuadraticRecurrencePrimitive
