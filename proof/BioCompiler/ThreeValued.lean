/-
  BioCompiler.ThreeValued — Three-Valued Logic for Biological Verdicts
  Reference: DOC-03 (SDD) §3.5.3, DOC-10 (Deterministic Methods) §8
-/

namespace BioCompiler

inductive Verdict where
  | PASS      : Verdict
  | FAIL      : Verdict
  | UNCERTAIN : Verdict
  deriving DecidableEq, Repr, BEq

namespace Verdict

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

def not : Verdict → Verdict
  | PASS      => FAIL
  | FAIL      => PASS
  | UNCERTAIN => UNCERTAIN

-- Key reduction lemmas
@[simp] theorem and_PASS (v : Verdict) : Verdict.and PASS v = v := by cases v <;> rfl
@[simp] theorem and_FAIL_left (v : Verdict) : Verdict.and FAIL v = FAIL := by cases v <;> rfl
@[simp] theorem and_FAIL_right (v : Verdict) : Verdict.and v FAIL = FAIL := by cases v <;> rfl

-- Discrimination
private theorem FAIL_ne_PASS : FAIL ≠ PASS := by intro h; cases h
theorem UNCERTAIN_ne_PASS : UNCERTAIN ≠ PASS := by intro h; cases h

-- Verdict.and init hd = PASS ↔ init = PASS ∧ hd = PASS
theorem and_eq_PASS_iff {init hd : Verdict} :
    Verdict.and init hd = PASS ↔ init = PASS ∧ hd = PASS := by
  cases init <;> cases hd <;> simp [Verdict.and]

-- Algebraic properties
theorem and_comm (v₁ v₂ : Verdict) : Verdict.and v₁ v₂ = Verdict.and v₂ v₁ := by
  cases v₁ <;> cases v₂ <;> rfl

theorem or_comm (v₁ v₂ : Verdict) : Verdict.or v₁ v₂ = Verdict.or v₂ v₁ := by
  cases v₁ <;> cases v₂ <;> rfl

theorem and_assoc (v₁ v₂ v₃ : Verdict) : Verdict.and (Verdict.and v₁ v₂) v₃ = Verdict.and v₁ (Verdict.and v₂ v₃) := by
  cases v₁ <;> cases v₂ <;> cases v₃ <;> rfl

theorem or_assoc (v₁ v₂ v₃ : Verdict) : Verdict.or (Verdict.or v₁ v₂) v₃ = Verdict.or v₁ (Verdict.or v₂ v₃) := by
  cases v₁ <;> cases v₂ <;> cases v₃ <;> rfl

theorem and_idem (v : Verdict) : Verdict.and v v = v := by cases v <;> rfl
theorem or_idem (v : Verdict) : Verdict.or v v = v := by cases v <;> rfl

-- Soundness-critical properties
theorem and_pass_pass : ∀ v₁ v₂, v₁ = PASS → v₂ = PASS → Verdict.and v₁ v₂ = PASS := by
  intro v₁ v₂ h₁ h₂; rw [h₁, h₂]; rfl

theorem and_pass_uncertain : Verdict.and PASS UNCERTAIN = UNCERTAIN := rfl
theorem and_uncertain_pass : Verdict.and UNCERTAIN PASS = UNCERTAIN := rfl

/-- foldl starting from non-PASS init never reaches PASS. -/
theorem foldl_ne_pass_of_ne_pass (init : Verdict) (vs : List Verdict) :
    init ≠ PASS → (vs.foldl Verdict.and init) ≠ PASS := by
  intro h_init
  induction vs generalizing init with
  | nil => simp [List.foldl_nil]; exact h_init
  | cons hd tl ih =>
    simp only [List.foldl_cons]
    intro h_pass
    have h_and_ne : Verdict.and init hd ≠ PASS := by
      intro h; exact h_init (and_eq_PASS_iff.mp h).left
    exact ih (Verdict.and init hd) h_and_ne h_pass

/-- foldl starting from UNCERTAIN never reaches PASS. -/
theorem foldl_uncertain_ne_pass (vs : List Verdict) :
    (vs.foldl Verdict.and UNCERTAIN) ≠ PASS :=
  foldl_ne_pass_of_ne_pass UNCERTAIN vs UNCERTAIN_ne_PASS

/-- foldl starting from FAIL never reaches PASS. -/
theorem foldl_fail_ne_pass (vs : List Verdict) :
    (vs.foldl Verdict.and FAIL) ≠ PASS :=
  foldl_ne_pass_of_ne_pass FAIL vs FAIL_ne_PASS

/-- If foldl from PASS yields PASS, every element must be PASS. -/
theorem foldl_and_pass_implies_all_pass (vs : List Verdict) :
    (vs.foldl Verdict.and PASS) = PASS → ∀ v ∈ vs, v = PASS := by
  intro h v hv
  induction vs with
  | nil => simp at hv
  | cons hd tl ih =>
    simp only [List.foldl_cons, and_PASS] at h
    -- h : tl.foldl Verdict.and hd = PASS
    -- Prove hd = PASS: if hd ≠ PASS, foldl can't reach PASS
    have h_hd_pass : hd = PASS := by
      apply Decidable.byContradiction
      intro h_ne
      exact absurd h (foldl_ne_pass_of_ne_pass hd tl h_ne)
    -- h is already simplified to: tl.foldl Verdict.and PASS = PASS
    -- v ∈ hd :: tl
    rw [h_hd_pass] at h
    obtain h_eq | h_mem := List.mem_cons.mp hv
    · rw [h_eq]; exact h_hd_pass
    · exact ih h h_mem

/-- If foldl from PASS yields PASS, no element is FAIL. -/
theorem and_pass_all_pass (vs : List Verdict) :
    (vs.foldl Verdict.and PASS) = PASS → ∀ v ∈ vs, v ≠ FAIL := by
  intro h v hv
  have h_v_pass := foldl_and_pass_implies_all_pass vs h v hv
  rw [h_v_pass]; intro h; cases h

/-- Information ordering on verdicts. -/
def ordering : Verdict → Nat
  | PASS => 2
  | UNCERTAIN => 1
  | FAIL => 0

theorem and_monotone_left (v₁ v₂ v : Verdict) :
    ordering v₁ ≤ ordering v₂ → ordering (Verdict.and v₁ v) ≤ ordering (Verdict.and v₂ v) := by
  intro h
  cases v₁ <;> cases v₂ <;> cases v <;>
    simp only [Verdict.and, ordering] <;>
    first | rfl | omega | (nomatch h)

theorem ordering_pass_highest (v : Verdict) : ordering v ≤ ordering PASS := by
  cases v <;> simp [ordering] <;> decide

theorem ordering_fail_lowest (v : Verdict) : ordering FAIL ≤ ordering v := by
  cases v <;> simp [ordering] <;> decide

end Verdict

end BioCompiler
