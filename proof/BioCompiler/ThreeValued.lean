/-
  BioCompiler.ThreeValued — Three-Valued Logic for Biological Verdicts

  This module defines the three-valued logic (PASS/FAIL/UNCERTAIN) used by the
  BioCompiler type system, proves its key algebraic properties, and establishes
  the soundness of composition (conjunction preserves PASS guarantees).

  Reference: DOC-03 (SDD) §3.5.3, DOC-10 (Deterministic Methods) §8
-/

namespace BioCompiler

/-- Three-valued logic for type-check verdicts. -/
inductive Verdict where
  | PASS      : Verdict  -- Guaranteed to hold in every possible execution
  | FAIL      : Verdict  -- Guaranteed NOT to hold in any possible execution
  | UNCERTAIN : Verdict  -- Cannot determine from available information
  deriving DecidableEq, Repr, BEq

namespace Verdict

/-- Conjunction (AND) in three-valued logic.
    Preserves soundness: if both operands are PASS, the result is PASS.
    If either operand is FAIL, the result is FAIL (a known violation). -/
def and : Verdict → Verdict → Verdict
  | PASS,      PASS      => PASS
  | PASS,      UNCERTAIN => UNCERTAIN
  | PASS,      FAIL      => FAIL
  | UNCERTAIN, PASS      => UNCERTAIN
  | UNCERTAIN, UNCERTAIN => UNCERTAIN
  | UNCERTAIN, FAIL      => FAIL
  | FAIL,      PASS      => FAIL
  | FAIL,      UNCERTAIN => FAIL
  | FAIL,      FAIL      => FAIL

/-- Disjunction (OR) in three-valued logic. -/
def or : Verdict → Verdict → Verdict
  | PASS,      PASS      => PASS
  | PASS,      UNCERTAIN => PASS
  | PASS,      FAIL      => PASS
  | UNCERTAIN, PASS      => PASS
  | UNCERTAIN, UNCERTAIN => UNCERTAIN
  | UNCERTAIN, FAIL      => UNCERTAIN
  | FAIL,      PASS      => PASS
  | FAIL,      UNCERTAIN => UNCERTAIN
  | FAIL,      FAIL      => FAIL

/-- Negation in three-valued logic. -/
def not : Verdict → Verdict
  | PASS      => FAIL
  | FAIL      => PASS
  | UNCERTAIN => UNCERTAIN

/-- Notation for three-valued conjunction. -/
infixl:35 " ⊓ " => Verdict.and

/-- Notation for three-valued disjunction. -/
infixl:30 " ⊔ " => Verdict.or

prefix:40 "∼" => Verdict.not

-- ==============================================================================
-- Algebraic Properties of Three-Valued Logic
-- ==============================================================================

/-- Commutativity of conjunction. -/
theorem and_comm (v₁ v₂ : Verdict) : v₁ ⊓ v₂ = v₂ ⊓ v₁ := by
  cases v₁ <;> cases v₂ <;> rfl

/-- Commutativity of disjunction. -/
theorem or_comm (v₁ v₂ : Verdict) : v₁ ⊔ v₂ = v₂ ⊔ v₁ := by
  cases v₁ <;> cases v₂ <;> rfl

/-- Associativity of conjunction. -/
theorem and_assoc (v₁ v₂ v₃ : Verdict) : (v₁ ⊓ v₂) ⊓ v₃ = v₁ ⊓ (v₂ ⊓ v₃) := by
  cases v₁ <;> cases v₂ <;> cases v₃ <;> rfl

/-- Associativity of disjunction. -/
theorem or_assoc (v₁ v₂ v₃ : Verdict) : (v₁ ⊔ v₂) ⊔ v₃ = v₁ ⊔ (v₂ ⊔ v₃) := by
  cases v₁ <;> cases v₂ <;> cases v₃ <;> rfl

/-- Idempotency of conjunction. -/
theorem and_idem (v : Verdict) : v ⊓ v = v := by
  cases v <;> rfl

/-- Idempotency of disjunction. -/
theorem or_idem (v : Verdict) : v ⊔ v = v := by
  cases v <;> rfl

-- ==============================================================================
-- SOUNDNESS-CRITICAL PROPERTIES
--
-- These are the key theorems that make three-valued logic safe for
-- verification: PASS verdicts propagate through conjunction, and FAIL
-- verdicts are sticky (they cannot be hidden by conjunction with PASS).
-- ==============================================================================

/-- THEOREM (Conjunction Soundness): If both conjuncts are PASS,
    the conjunction is PASS. -/
theorem and_pass_pass : ∀ v₁ v₂, v₁ = PASS → v₂ = PASS → v₁ ⊓ v₂ = PASS := by
  intro v₁ v₂ h₁ h₂
  rw [h₁, h₂]
  rfl

/-- THEOREM: FAIL is left-annihilating for conjunction. -/
theorem and_fail_left (v : Verdict) : FAIL ⊓ v = FAIL := by
  cases v <;> rfl

/-- THEOREM: FAIL is right-annihilating for conjunction. -/
theorem and_fail_right (v : Verdict) : v ⊓ FAIL = FAIL := by
  cases v <;> rfl

/-- THEOREM (FAIL Detection): If the conjunction is PASS,
    neither conjunct is FAIL. -/
theorem and_pass_no_fail (v₁ v₂ : Verdict) : v₁ ⊓ v₂ = PASS → v₁ ≠ FAIL ∧ v₂ ≠ FAIL := by
  intro h
  constructor
  · intro h₁; rw [h₁] at h; exact h.elim
  · intro h₂; rw [h₂] at h; exact h.elim

/-- THEOREM (PASS Implies No FAIL): If the conjunction of all type
    predicates yields PASS, then no individual predicate returned FAIL. -/
theorem and_pass_all_pass (vs : List Verdict) :
    vs.foldl Verdict.and PASS = PASS → ∀ v ∈ vs, v ≠ FAIL := by
  intro h v hv
  induction vs generalizing h with
  | nil => simp at hv
  | cons hd tl ih =>
    simp [List.foldl_cons] at h
    cases h_head : hd with
    | PASS =>
      have h_tail := ih (by simp [h_head] at h; exact h)
      cases hv with
      | head => intro h_fail; rw [h_fail] at h_head; exact h_head.elim
      | tail _ hv_mem => exact h_tail v hv_mem
    | FAIL =>
      simp [h_head, Verdict.and] at h
      exact h.elim
    | UNCERTAIN =>
      simp [h_head, Verdict.and] at h
      exact h.elim

/-- THEOREM: PASS ⊓ UNCERTAIN = UNCERTAIN. -/
theorem and_pass_uncertain : PASS ⊓ UNCERTAIN = UNCERTAIN := rfl

/-- THEOREM: UNCERTAIN ⊓ PASS = UNCERTAIN. -/
theorem and_uncertain_pass : UNCERTAIN ⊓ PASS = UNCERTAIN := rfl

/-- THEOREM: If conjunction of verdicts yields PASS, every individual
    verdict must be PASS (not FAIL, and not UNCERTAIN, because
    UNCERTAIN ⊓ PASS = UNCERTAIN ≠ PASS). -/
theorem foldl_and_pass_implies_all_pass (vs : List Verdict) :
    vs.foldl Verdict.and PASS = PASS → ∀ v ∈ vs, v = PASS := by
  intro h v hv
  induction vs generalizing h with
  | nil => simp at hv
  | cons hd tl ih =>
    simp [List.foldl_cons] at h
    cases h_hd : hd with
    | PASS =>
      have h_tail := ih (by simp [h_hd] at h; exact h) v
      cases hv with
      | head => exact h_hd
      | tail _ hv_mem => exact h_tail hv_mem
    | FAIL =>
      simp [h_hd, Verdict.and] at h
    | UNCERTAIN =>
      simp [h_hd, Verdict.and] at h

/-- Information ordering on verdicts: PASS > UNCERTAIN > FAIL. -/
def ordering : Verdict → Nat
  | PASS => 2
  | UNCERTAIN => 1
  | FAIL => 0

/-- Monotonicity of conjunction with respect to the information ordering. -/
theorem and_monotone_left (v₁ v₂ v : Verdict) :
    ordering v₁ ≤ ordering v₂ → ordering (v₁ ⊓ v) ≤ ordering (v₂ ⊓ v) := by
  intro h
  cases v₁ <;> cases v₂ <;> cases v <;> simp [ordering] at h ⊢ <;> omega

/-- PASS has the highest information value. -/
theorem ordering_pass_highest (v : Verdict) : ordering v ≤ ordering PASS := by
  cases v <;> simp [ordering] <;> omega

/-- FAIL has the lowest information value. -/
theorem ordering_fail_lowest (v : Verdict) : ordering FAIL ≤ ordering v := by
  cases v <;> simp [ordering] <;> omega

end Verdict

end BioCompiler
