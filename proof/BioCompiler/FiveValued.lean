/-
  BioCompiler.FiveValued — Five-Valued Logic for Biological Verdicts

  Extends the three-valued logic (PASS/UNCERTAIN/FAIL) with intermediate
  verdicts LIKELY_PASS and LIKELY_FAIL, providing finer-grained confidence
  distinctions while preserving the algebraic structure.

  The five values are ordered by confidence:
    PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL

  Key properties:
  1. five_valued_and refines three_valued_and (projection to 3 values)
  2. five_valued_and PASS PASS = PASS
  3. five_valued_and FAIL x = FAIL for all x
  4. Associativity and commutativity of five_valued_and

  REFERENCE: DOC-03 (SDD) §3.5.3, DOC-10 (Deterministic Methods) §8
-/

import BioCompiler.ThreeValued

namespace BioCompiler

-- ==============================================================================
-- Five-Valued Verdict Type
-- ==============================================================================

inductive FiveVerdict where
  | PASS         : FiveVerdict
  | LIKELY_PASS  : FiveVerdict
  | UNCERTAIN    : FiveVerdict
  | LIKELY_FAIL  : FiveVerdict
  | FAIL         : FiveVerdict
  deriving DecidableEq, Repr, BEq

namespace FiveVerdict

-- ==============================================================================
-- Conjunction (AND) — Kleene-style: takes the minimum
-- ==============================================================================

def and : FiveVerdict → FiveVerdict → FiveVerdict
  | PASS,         PASS         => PASS
  | PASS,         LIKELY_PASS  => LIKELY_PASS
  | PASS,         UNCERTAIN    => UNCERTAIN
  | PASS,         LIKELY_FAIL  => LIKELY_FAIL
  | PASS,         FAIL         => FAIL
  | LIKELY_PASS,  PASS         => LIKELY_PASS
  | LIKELY_PASS,  LIKELY_PASS  => LIKELY_PASS
  | LIKELY_PASS,  UNCERTAIN    => UNCERTAIN
  | LIKELY_PASS,  LIKELY_FAIL  => LIKELY_FAIL
  | LIKELY_PASS,  FAIL         => FAIL
  | UNCERTAIN,    PASS         => UNCERTAIN
  | UNCERTAIN,    LIKELY_PASS  => UNCERTAIN
  | UNCERTAIN,    UNCERTAIN    => UNCERTAIN
  | UNCERTAIN,    LIKELY_FAIL  => LIKELY_FAIL
  | UNCERTAIN,    FAIL         => FAIL
  | LIKELY_FAIL,  PASS         => LIKELY_FAIL
  | LIKELY_FAIL,  LIKELY_PASS  => LIKELY_FAIL
  | LIKELY_FAIL,  UNCERTAIN    => LIKELY_FAIL
  | LIKELY_FAIL,  LIKELY_FAIL  => LIKELY_FAIL
  | LIKELY_FAIL,  FAIL         => FAIL
  | FAIL,         PASS         => FAIL
  | FAIL,         LIKELY_PASS  => FAIL
  | FAIL,         UNCERTAIN    => FAIL
  | FAIL,         LIKELY_FAIL  => FAIL
  | FAIL,         FAIL         => FAIL

-- ==============================================================================
-- Disjunction (OR) — Kleene-style: takes the maximum
-- ==============================================================================

def or : FiveVerdict → FiveVerdict → FiveVerdict
  | PASS,         PASS         => PASS
  | PASS,         LIKELY_PASS  => PASS
  | PASS,         UNCERTAIN    => PASS
  | PASS,         LIKELY_FAIL  => PASS
  | PASS,         FAIL         => PASS
  | LIKELY_PASS,  PASS         => PASS
  | LIKELY_PASS,  LIKELY_PASS  => LIKELY_PASS
  | LIKELY_PASS,  UNCERTAIN    => LIKELY_PASS
  | LIKELY_PASS,  LIKELY_FAIL  => LIKELY_PASS
  | LIKELY_PASS,  FAIL         => LIKELY_PASS
  | UNCERTAIN,    PASS         => PASS
  | UNCERTAIN,    LIKELY_PASS  => LIKELY_PASS
  | UNCERTAIN,    UNCERTAIN    => UNCERTAIN
  | UNCERTAIN,    LIKELY_FAIL  => UNCERTAIN
  | UNCERTAIN,    FAIL         => UNCERTAIN
  | LIKELY_FAIL,  PASS         => PASS
  | LIKELY_FAIL,  LIKELY_PASS  => LIKELY_PASS
  | LIKELY_FAIL,  UNCERTAIN    => UNCERTAIN
  | LIKELY_FAIL,  LIKELY_FAIL  => LIKELY_FAIL
  | LIKELY_FAIL,  FAIL         => LIKELY_FAIL
  | FAIL,         PASS         => PASS
  | FAIL,         LIKELY_PASS  => LIKELY_PASS
  | FAIL,         UNCERTAIN    => UNCERTAIN
  | FAIL,         LIKELY_FAIL  => LIKELY_FAIL
  | FAIL,         FAIL         => FAIL

-- ==============================================================================
-- Negation
-- ==============================================================================

def not : FiveVerdict → FiveVerdict
  | PASS         => FAIL
  | LIKELY_PASS  => LIKELY_FAIL
  | UNCERTAIN    => UNCERTAIN
  | LIKELY_FAIL  => LIKELY_PASS
  | FAIL         => PASS

-- ==============================================================================
-- Projection from Five-Valued to Three-Valued
-- ==============================================================================

/-- Project a five-valued verdict to a three-valued verdict.

    The projection maps:
    - PASS → PASS (definite pass stays definite pass)
    - LIKELY_PASS → UNCERTAIN (likely pass is still uncertain in 3-valued)
    - UNCERTAIN → UNCERTAIN
    - LIKELY_FAIL → UNCERTAIN (likely fail is still uncertain in 3-valued)
    - FAIL → FAIL (definite fail stays definite fail)

    This is the abstraction function in the refinement relationship:
    the 5-valued logic refines the 3-valued logic by splitting UNCERTAIN
    into three sub-categories. -/
def project : FiveVerdict → Verdict
  | PASS         => Verdict.PASS
  | LIKELY_PASS  => Verdict.UNCERTAIN
  | UNCERTAIN    => Verdict.UNCERTAIN
  | LIKELY_FAIL  => Verdict.UNCERTAIN
  | FAIL         => Verdict.FAIL

-- ==============================================================================
-- Key Theorems
-- ==============================================================================

-- AND with PASS is identity (left and right)
@[simp] theorem and_PASS (v : FiveVerdict) : FiveVerdict.and PASS v = v := by
  cases v <;> rfl

@[simp] theorem and_PASS_right (v : FiveVerdict) : FiveVerdict.and v PASS = v := by
  cases v <;> rfl

-- AND with FAIL is always FAIL
@[simp] theorem and_FAIL_left (v : FiveVerdict) : FiveVerdict.and FAIL v = FAIL := by
  cases v <;> rfl

@[simp] theorem and_FAIL_right (v : FiveVerdict) : FiveVerdict.and v FAIL = FAIL := by
  cases v <;> rfl

-- Discrimination lemmas
private theorem FAIL_ne_PASS : FAIL ≠ PASS := by intro h; cases h
private theorem FAIL_ne_LIKELY_PASS : FAIL ≠ LIKELY_PASS := by intro h; cases h
private theorem FAIL_ne_LIKELY_FAIL : FAIL ≠ LIKELY_FAIL := by intro h; cases h
theorem UNCERTAIN_ne_PASS : UNCERTAIN ≠ PASS := by intro h; cases h
theorem UNCERTAIN_ne_FAIL : UNCERTAIN ≠ FAIL := by intro h; cases h
theorem LIKELY_PASS_ne_PASS : LIKELY_PASS ≠ PASS := by intro h; cases h
theorem LIKELY_FAIL_ne_FAIL : LIKELY_FAIL ≠ FAIL := by intro h; cases h

-- ==============================================================================
-- Algebraic Properties
-- ==============================================================================

theorem and_comm (v₁ v₂ : FiveVerdict) : FiveVerdict.and v₁ v₂ = FiveVerdict.and v₂ v₁ := by
  cases v₁ <;> cases v₂ <;> rfl

theorem or_comm (v₁ v₂ : FiveVerdict) : FiveVerdict.or v₁ v₂ = FiveVerdict.or v₂ v₁ := by
  cases v₁ <;> cases v₂ <;> rfl

theorem and_assoc (v₁ v₂ v₃ : FiveVerdict) :
    FiveVerdict.and (FiveVerdict.and v₁ v₂) v₃ = FiveVerdict.and v₁ (FiveVerdict.and v₂ v₃) := by
  cases v₁ <;> cases v₂ <;> cases v₃ <;> rfl

theorem or_assoc (v₁ v₂ v₃ : FiveVerdict) :
    FiveVerdict.or (FiveVerdict.or v₁ v₂) v₃ = FiveVerdict.or v₁ (FiveVerdict.or v₂ v₃) := by
  cases v₁ <;> cases v₂ <;> cases v₃ <;> rfl

theorem and_idem (v : FiveVerdict) : FiveVerdict.and v v = v := by cases v <;> rfl
theorem or_idem (v : FiveVerdict) : FiveVerdict.or v v = v := by cases v <;> rfl

-- ==============================================================================
-- Refinement: five_valued_and refines three_valued_and
-- ==============================================================================

/-- THEOREM (Refinement): The projection of five_valued_and equals
    three_valued_and of the projections. This means the five-valued
    conjunction refines the three-valued conjunction:

      project (five_valued_and a b) = three_valued_and (project a) (project b)

    This establishes that the 5-valued logic is a conservative refinement
    of the 3-valued logic: any conclusion drawn in the 3-valued world
    from 5-valued verdicts via projection is the same as if we had
    used 3-valued logic from the start. -/
theorem and_project_refines (a b : FiveVerdict) :
    project (FiveVerdict.and a b) = Verdict.and (project a) (project b) := by
  cases a <;> cases b <;> rfl

/-- THEOREM (OR Refinement): The projection of five_valued_or equals
    three_valued_or of the projections. -/
theorem or_project_refines (a b : FiveVerdict) :
    project (FiveVerdict.or a b) = Verdict.or (project a) (project b) := by
  cases a <;> cases b <;> rfl

-- ==============================================================================
-- Soundness-Critical Properties (mirroring ThreeValued.lean)
-- ==============================================================================

/-- five_valued_and PASS PASS = PASS -/
theorem and_pass_pass : FiveVerdict.and PASS PASS = PASS := rfl

/-- five_valued_and FAIL x = FAIL for all x -/
theorem and_fail_absorb (x : FiveVerdict) : FiveVerdict.and FAIL x = FAIL := by
  cases x <;> rfl

/-- five_valued_and x FAIL = FAIL for all x -/
theorem and_fail_absorb_right (x : FiveVerdict) : FiveVerdict.and x FAIL = FAIL := by
  cases x <;> rfl

/-- If five_valued_and a b = PASS, then a = PASS and b = PASS -/
theorem and_eq_PASS_iff {a b : FiveVerdict} :
    FiveVerdict.and a b = PASS ↔ a = PASS ∧ b = PASS := by
  cases a <;> cases b <;> simp [FiveVerdict.and]

/-- If five_valued_and a b ≠ FAIL, then a ≠ FAIL and b ≠ FAIL -/
theorem and_ne_fail_both {a b : FiveVerdict} :
    FiveVerdict.and a b ≠ FAIL → a ≠ FAIL ∧ b ≠ FAIL := by
  intro h
  constructor
  · intro ha; rw [ha, and_FAIL_left] at h; exact h rfl
  · intro hb; rw [hb, and_FAIL_right] at h; exact h rfl

-- ==============================================================================
-- Ordering and Monotonicity
-- ==============================================================================

/-- Information ordering on five-valued verdicts.
    Higher values are more positive. -/
def ordering : FiveVerdict → Nat
  | PASS         => 4
  | LIKELY_PASS  => 3
  | UNCERTAIN    => 2
  | LIKELY_FAIL  => 1
  | FAIL         => 0

/-- AND is monotonically decreasing in each argument:
    if v₁ ≤ v₂ in the ordering, then and v₁ v ≤ and v₂ v -/
theorem and_monotone_left (v₁ v₂ v : FiveVerdict) :
    ordering v₁ ≤ ordering v₂ → ordering (FiveVerdict.and v₁ v) ≤ ordering (FiveVerdict.and v₂ v) := by
  intro h
  cases v₁ <;> cases v₂ <;> cases v <;>
    simp only [FiveVerdict.and, ordering] <;>
    first | rfl | omega | (nomatch h)

theorem ordering_pass_highest (v : FiveVerdict) : ordering v ≤ ordering PASS := by
  cases v <;> simp [ordering] <;> decide

theorem ordering_fail_lowest (v : FiveVerdict) : ordering FAIL ≤ ordering v := by
  cases v <;> simp [ordering] <;> decide

-- ==============================================================================
-- Foldl Properties (mirroring ThreeValued.lean)
-- ==============================================================================

/-- foldl starting from non-PASS init never reaches PASS. -/
theorem foldl_ne_pass_of_ne_pass (init : FiveVerdict) (vs : List FiveVerdict) :
    init ≠ PASS → (vs.foldl FiveVerdict.and init) ≠ PASS := by
  intro h_init
  induction vs generalizing init with
  | nil => simp [List.foldl_nil]; exact h_init
  | cons hd tl ih =>
    simp only [List.foldl_cons]
    intro h_pass
    have h_and_ne : FiveVerdict.and init hd ≠ PASS := by
      intro h; exact h_init (and_eq_PASS_iff.mp h).left
    exact ih (FiveVerdict.and init hd) h_and_ne h_pass

/-- foldl starting from FAIL never reaches PASS. -/
theorem foldl_fail_ne_pass (vs : List FiveVerdict) :
    (vs.foldl FiveVerdict.and FAIL) ≠ PASS :=
  foldl_ne_pass_of_ne_pass FAIL vs FAIL_ne_PASS

/-- If foldl from PASS yields PASS, every element must be PASS. -/
theorem foldl_and_pass_implies_all_pass (vs : List FiveVerdict) :
    (vs.foldl FiveVerdict.and PASS) = PASS → ∀ v ∈ vs, v = PASS := by
  intro h v hv
  induction vs with
  | nil => simp at hv
  | cons hd tl ih =>
    simp only [List.foldl_cons, and_PASS] at h
    have h_hd_pass : hd = PASS := by
      apply Decidable.byContradiction
      intro h_ne
      exact absurd h (foldl_ne_pass_of_ne_pass hd tl h_ne)
    rw [h_hd_pass] at h
    obtain h_eq | h_mem := List.mem_cons.mp hv
    · rw [h_eq]; exact h_hd_pass
    · exact ih h h_mem

-- ==============================================================================
-- Projection Properties
-- ==============================================================================

/-- Projecting PASS gives PASS -/
theorem project_pass : project PASS = Verdict.PASS := rfl

/-- Projecting FAIL gives FAIL -/
theorem project_fail : project FAIL = Verdict.FAIL := rfl

/-- Projecting UNCERTAIN gives UNCERTAIN -/
theorem project_uncertain : project UNCERTAIN = Verdict.UNCERTAIN := rfl

/-- Projecting LIKELY_PASS gives UNCERTAIN -/
theorem project_likely_pass : project LIKELY_PASS = Verdict.UNCERTAIN := rfl

/-- Projecting LIKELY_FAIL gives UNCERTAIN -/
theorem project_likely_fail : project LIKELY_FAIL = Verdict.UNCERTAIN := rfl

/-- If project v = Verdict.PASS, then v = PASS -/
theorem project_pass_implies_pass {v : FiveVerdict} :
    project v = Verdict.PASS → v = PASS := by
  intro h; cases v <;> simp [project] at h <;> try rfl <;> cases h

/-- If project v = Verdict.FAIL, then v = FAIL -/
theorem project_fail_implies_fail {v : FiveVerdict} :
    project v = Verdict.FAIL → v = FAIL := by
  intro h; cases v <;> simp [project] at h <;> try rfl <;> cases h

-- ==============================================================================
-- De Morgan's Laws
-- ==============================================================================

/-- De Morgan: not (and a b) = or (not a) (not b) -/
theorem de_morgan_and (a b : FiveVerdict) :
    FiveVerdict.not (FiveVerdict.and a b) = FiveVerdict.or (FiveVerdict.not a) (FiveVerdict.not b) := by
  cases a <;> cases b <;> rfl

/-- De Morgan: not (or a b) = and (not a) (not b) -/
theorem de_morgan_or (a b : FiveVerdict) :
    FiveVerdict.not (FiveVerdict.or a b) = FiveVerdict.and (FiveVerdict.not a) (FiveVerdict.not b) := by
  cases a <;> cases b <;> rfl

-- ==============================================================================
-- Not Involutive
-- ==============================================================================

/-- not(not(v)) = v for all five-valued verdicts -/
theorem not_involutive (v : FiveVerdict) : FiveVerdict.not (FiveVerdict.not v) = v := by
  cases v <;> rfl

-- ==============================================================================
-- Absorption Laws
-- ==============================================================================

/-- AND absorbs OR: and a (or a b) = a -/
theorem and_absorbs_or (a b : FiveVerdict) :
    FiveVerdict.and a (FiveVerdict.or a b) = a := by
  cases a <;> cases b <;> rfl

/-- OR absorbs AND: or a (and a b) = a -/
theorem or_absorbs_and (a b : FiveVerdict) :
    FiveVerdict.or a (FiveVerdict.and a b) = a := by
  cases a <;> cases b <;> rfl

-- ==============================================================================
-- Distributivity Laws
-- ==============================================================================

/-- AND distributes over OR: and a (or b c) = or (and a b) (and a c) -/
theorem and_distrib_or (a b c : FiveVerdict) :
    FiveVerdict.and a (FiveVerdict.or b c) =
    FiveVerdict.or (FiveVerdict.and a b) (FiveVerdict.and a c) := by
  cases a <;> cases b <;> cases c <;> rfl

/-- OR distributes over AND: or a (and b c) = and (or a b) (or a c) -/
theorem or_distrib_and (a b c : FiveVerdict) :
    FiveVerdict.or a (FiveVerdict.and b c) =
    FiveVerdict.and (FiveVerdict.or a b) (FiveVerdict.or a c) := by
  cases a <;> cases b <;> cases c <;> rfl

-- ==============================================================================
-- De Morgan Refinement (projection consistency)
-- ==============================================================================

/-- De Morgan is consistent under projection: project(not v) = Verdict.not(project v) -/
theorem not_project_refines (v : FiveVerdict) :
    project (FiveVerdict.not v) = Verdict.not (project v) := by
  cases v <;> rfl

-- ==============================================================================
-- Monotonicity of OR (complements and_monotone_left for AND)
-- ==============================================================================

/-- OR is monotonically increasing in each argument -/
theorem or_monotone_left (v₁ v₂ v : FiveVerdict) :
    ordering v₁ ≤ ordering v₂ → ordering (FiveVerdict.or v₁ v) ≤ ordering (FiveVerdict.or v₂ v) := by
  intro h
  cases v₁ <;> cases v₂ <;> cases v <;>
    simp only [FiveVerdict.or, ordering] <;>
    first | rfl | omega | (nomatch h)

/-- AND is monotonically decreasing in the right argument -/
theorem and_monotone_right (v v₁ v₂ : FiveVerdict) :
    ordering v₁ ≤ ordering v₂ → ordering (FiveVerdict.and v v₁) ≤ ordering (FiveVerdict.and v v₂) := by
  intro h
  cases v <;> cases v₁ <;> cases v₂ <;>
    simp only [FiveVerdict.and, ordering] <;>
    first | rfl | omega | (nomatch h)

/-- OR is monotonically increasing in the right argument -/
theorem or_monotone_right (v v₁ v₂ : FiveVerdict) :
    ordering v₁ ≤ ordering v₂ → ordering (FiveVerdict.or v v₁) ≤ ordering (FiveVerdict.or v v₂) := by
  intro h
  cases v <;> cases v₁ <;> cases v₂ <;>
    simp only [FiveVerdict.or, ordering] <;>
    first | rfl | omega | (nomatch h)

-- ==============================================================================
-- Bounded Lattice Properties
-- ==============================================================================

/-- PASS is the top element for OR: or v PASS = PASS -/
theorem or_pass_top (v : FiveVerdict) : FiveVerdict.or v PASS = PASS := by
  cases v <;> rfl

/-- FAIL is the bottom element for OR: or v FAIL = v -/
@[simp] theorem or_FAIL_bottom (v : FiveVerdict) : FiveVerdict.or v FAIL = v := by
  cases v <;> rfl

/- Bounded lattice: and v PASS = v (identity) — already and_PASS_right -/
/- Bounded lattice: or v FAIL = v (identity) — already or_FAIL_bottom -/
/- Bounded lattice: and v FAIL = FAIL (absorbing) — already and_FAIL_right -/
/- Bounded lattice: or v PASS = PASS (absorbing) — already or_pass_top -/

end FiveVerdict

end BioCompiler
