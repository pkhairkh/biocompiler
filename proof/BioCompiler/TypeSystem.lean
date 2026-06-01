/-
  BioCompiler.TypeSystem — Type Predicates and Soundness Proof

  Central Theorem:

    ∀ (P : TypePredicate) (seq : Sequence) (ctx : CellularContext),
      evaluate P seq ctx = Verdict.PASS → propertyHolds P seq ctx

  "Well-typed genes don't go wrong." — after Milner, 1978

  PROOF STRATEGY:
  - Arithmetic predicates (CodonAdapted, GCInRange, InFrame): PASS condition = property.
    Proved via ite_fail_imp (contrapositive: ¬condition → FAIL = PASS, contradiction).
  - Scanner predicates (NoCrypticSplice, NoRestrictionSite, NoInstabilityMotif, NoCpGIsland):
    PASS → scanner = false → contrapositive of completeness → no match exists.
  - NoCrypticSplice uses dual-threshold (PASS/UNCERTAIN/FAIL):
    FAIL if cryptic site found, UNCERTAIN if borderline site found, PASS otherwise.
    PASS guarantee: all sites have score < uncertainLoThreshold.
  - NDFST predicate (SpliceCorrect): PASS → singleton output set → ctx matches.

  SORRY STATUS: 0 remaining. All proofs are sorry-free, including ndfstRun_complete
  (proved via ConsumesInput in NDFST.lean). NoCpGIsland added as 8th predicate.
  Dual-threshold NoCrypticSplice with PASS/UNCERTAIN/FAIL support.

  REFERENCE: DOC-03 (SDD) §3.5, DOC-10 (Deterministic Methods) §4
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners

namespace BioCompiler

open Verdict Sequence

-- ==============================================================================
-- Type Predicates
-- ==============================================================================

inductive TypePredicate where
  | SpliceCorrect (cellType : String) : TypePredicate
  | NoCrypticSplice : TypePredicate
  | CodonAdapted (organism : String) (threshold : Rat) : TypePredicate
  | GCInRange (lo hi : Rat) : TypePredicate
  | NoRestrictionSite (enzymeSites : List Sequence) : TypePredicate
  | InFrame (readingFrame : Nat) (exonBoundaries : List Nat) : TypePredicate
  | NoInstabilityMotif : TypePredicate
  | NoCpGIsland : TypePredicate
  deriving Repr

-- ==============================================================================
-- Auxiliary Lemmas for Proof Engineering
-- ==============================================================================

/-- If `if cond then PASS else FAIL = PASS`, then cond holds.
    Contrapositive: if ¬cond, the else branch gives FAIL ≠ PASS. -/
theorem ite_fail_imp {cond : Prop} [Decidable cond]
    (h : (if cond then PASS else FAIL) = PASS) : cond := by
  by_cases h_pos : cond
  · exact h_pos
  · rw [if_neg h_pos] at h; cases h

/-- If `if cond then FAIL else PASS = PASS`, then ¬cond holds.
    Contrapositive: if cond, the then branch gives FAIL ≠ PASS. -/
theorem ite_pass_imp_neg {cond : Prop} [Decidable cond]
    (h : (if cond then FAIL else PASS) = PASS) : ¬cond := by
  intro h_pos
  rw [if_pos h_pos] at h; cases h

/-- Bool: b ≠ true ↔ b = false. -/
theorem bool_ne_true_iff_false (b : Bool) : b ≠ true ↔ b = false := by
  cases b <;> simp

/-- Bool: (a || b) = false → a = false. -/
theorem Bool.or_false_left (a b : Bool) (h : (a || b) = false) : a = false := by
  cases a with
  | false => rfl
  | true => simp at h

/-- Bool: (a || b) = false → b = false. -/
theorem Bool.or_false_right (a b : Bool) (h : (a || b) = false) : b = false := by
  cases a with
  | false => simp at h; exact h
  | true => simp at h

/-- A singleton list has length 1. -/
theorem list_singleton_length {α : Type} {x : α} : ([x] : List α).length = 1 := by simp

/-- String equality from negation of inequality (using BEq). -/
theorem string_eq_of_not_ne (s₁ s₂ : String) (h : ¬(s₁ != s₂)) : s₁ = s₂ := by
  have h_beq : (s₁ == s₂) = true := by
    by_cases h_beq : (s₁ == s₂) = true
    · exact h_beq
    · exfalso
      have h_false : (s₁ == s₂) = false := by
        cases h_beq' : (s₁ == s₂) with
        | true => exact absurd h_beq' h_beq
        | false => rfl
      have : (s₁ != s₂) = true := by
        show (!(s₁ == s₂)) = true
        rw [h_false]
        rfl
      exact h this
  exact LawfulBEq.eq_of_beq h_beq

/-- Rat: ¬(a < b) ↔ b ≤ a. Proved using the underlying Bool-based definitions. -/
theorem Rat.not_lt_iff_le (a b : Rat) : ¬(a < b) ↔ b ≤ a := by
  constructor
  · intro h
    have : a.blt b = false := by
      cases h_blt : (a.blt b) with
      | true => exact absurd h_blt h
      | false => rfl
    exact this
  · intro h h_lt
    have h_blt_false : a.blt b = false := h
    have h_blt_true : a.blt b = true := h_lt
    rw [h_blt_false] at h_blt_true
    cases h_blt_true

-- ==============================================================================
-- Evaluation Function
-- ==============================================================================

def evaluate [SpliceSiteScanner] [CodonAdaptationIndex] [CpGIslandScanner]
    {State : Type} [DecidableEq State] [Inhabited State] [SplicingNDFST State] :
    TypePredicate → Sequence → CellularContext → Verdict
  | TypePredicate.SpliceCorrect cellType, seq, ctx =>
      if ctx.cellType != cellType then UNCERTAIN
      else
        match ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq with
        | [_] => PASS
        | _ => FAIL

  | TypePredicate.NoCrypticSplice, seq, _ =>
      if SpliceSiteScanner.hasCrypticSpliceSite seq = true then FAIL
      else if SpliceSiteScanner.hasBorderlineSpliceSite seq = true then UNCERTAIN
      else PASS

  | TypePredicate.CodonAdapted org threshold, seq, _ =>
      if CodonAdaptationIndex.computeCAI seq org ≥ threshold then PASS else FAIL

  | TypePredicate.GCInRange lo hi, seq, _ =>
      if lo ≤ gcContent seq ∧ gcContent seq ≤ hi then PASS else FAIL

  | TypePredicate.NoRestrictionSite enzymeSites, seq, _ =>
      if hasAnyRestrictionSite seq enzymeSites = true then FAIL else PASS

  | TypePredicate.InFrame rf boundaries, seq, _ =>
      if readingFrameConsistent boundaries rf = true ∧
         hasPrematureStop seq rf = false then PASS else FAIL

  | TypePredicate.NoInstabilityMotif, seq, _ =>
      if hasInstabilityMotif seq = true then FAIL else PASS

  | TypePredicate.NoCpGIsland, seq, _ =>
      if CpGIslandScanner.hasCpGIsland seq = true then FAIL else PASS

-- ==============================================================================
-- Property Semantics
-- ==============================================================================

def propertyHolds [SpliceSiteScanner] [CodonAdaptationIndex] [CpGIslandScanner]
    {State : Type} [DecidableEq State] [Inhabited State] [SplicingNDFST State] :
    TypePredicate → Sequence → CellularContext → Prop
  | TypePredicate.SpliceCorrect cellType, seq, ctx =>
      ctx.cellType = cellType ∧
        (ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq).length = 1
  | TypePredicate.NoCrypticSplice, seq, _ =>
      ∀ (pos : Nat) (site : SpliceSiteMatch),
        pos < seq.length → site.position = pos →
        site.score ≥ uncertainLoThreshold → False
  | TypePredicate.CodonAdapted org threshold, seq, _ =>
      CodonAdaptationIndex.computeCAI seq org ≥ threshold
  | TypePredicate.GCInRange lo hi, seq, _ =>
      lo ≤ gcContent seq ∧ gcContent seq ≤ hi
  | TypePredicate.NoRestrictionSite enzymeSites, seq, _ =>
      ∀ (site : Sequence) (pos : Nat),
        site ∈ enzymeSites → pos + site.length ≤ seq.length →
          (seq.drop pos).take site.length ≠ site
  | TypePredicate.InFrame rf boundaries, seq, _ =>
      readingFrameConsistent boundaries rf = true ∧ hasPrematureStop seq rf = false
  | TypePredicate.NoInstabilityMotif, seq, _ =>
      (∀ (pos : Nat), pos + atttaMotif.length ≤ seq.length →
        (seq.drop pos).take atttaMotif.length ≠ atttaMotif) ∧
      (∀ (pos : Nat), pos + uRichMotif.length ≤ seq.length →
        (seq.drop pos).take uRichMotif.length ≠ uRichMotif)
  | TypePredicate.NoCpGIsland, seq, _ =>
      ∀ (pos : Nat),
        pos + cpgIslandWindowSize ≤ seq.length →
          let window := (seq.drop pos).take cpgIslandWindowSize
          ((window.count Nucleotide.G + window.count Nucleotide.C : Rat) / window.length < cpgIslandGCThreshold) ∨
          (∃ (cpgCount : Nat),
            cpgCount = (List.zipWith (· == ·) window (window.drop 1)).count true ∧
            (cpgCount : Rat) * window.length <
              cpgIslandObsExpThreshold * (window.count Nucleotide.C) * (window.count Nucleotide.G))

-- ==============================================================================
-- Soundness Theorem
-- ==============================================================================

theorem type_soundness [SpliceSiteScanner] [CodonAdaptationIndex] [CpGIslandScanner]
    {State : Type} [DecidableEq State] [Inhabited State] [SplicingNDFST State]
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    @evaluate _ _ _ State _ _ _ P seq ctx = PASS →
    @propertyHolds _ _ _ State _ _ _ P seq ctx := by
  intro h_pass
  cases P with

  | SpliceCorrect cellType =>
    simp only [evaluate] at h_pass
    -- Case-split on the Bool condition
    cases h_cond : (ctx.cellType != cellType) with
    | true =>
      -- If ctx.cellType != cellType, evaluate returns UNCERTAIN ≠ PASS
      simp [h_cond] at h_pass
    | false =>
      -- ctx.cellType = cellType (from ¬(ctx.cellType != cellType))
      have h_cell_eq : ctx.cellType = cellType := by
        have : ¬(ctx.cellType != cellType) := by
          intro h; rw [h] at h_cond; cases h_cond
        exact string_eq_of_not_ne _ _ this
      simp [h_cond] at h_pass
      -- h_pass now concerns only the match on ndfstUniqueOutputSet
      have h_list_len : (ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq).length = 1 := by
        cases h_list : ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq with
        | nil => simp [h_list] at h_pass
        | cons hd tl =>
          cases tl with
          | nil => simp
          | cons hd' tl' => simp [h_list] at h_pass
      simp only [propertyHolds]
      exact ⟨h_cell_eq, h_list_len⟩

  | NoCrypticSplice =>
    simp only [evaluate] at h_pass
    have h_not_cryptic : SpliceSiteScanner.hasCrypticSpliceSite seq ≠ true := by
      intro h; rw [if_pos h] at h_pass; cases h_pass
    have h_false_cryptic : SpliceSiteScanner.hasCrypticSpliceSite seq = false :=
      (bool_ne_true_iff_false _).mp h_not_cryptic
    have h_not_borderline : SpliceSiteScanner.hasBorderlineSpliceSite seq ≠ true := by
      intro h
      have : (if SpliceSiteScanner.hasCrypticSpliceSite seq = true then FAIL
              else if SpliceSiteScanner.hasBorderlineSpliceSite seq = true then UNCERTAIN
              else PASS) = UNCERTAIN := by
        rw [if_neg h_not_cryptic, if_pos h]
      rw [this] at h_pass; cases h_pass
    have h_false_borderline : SpliceSiteScanner.hasBorderlineSpliceSite seq = false :=
      (bool_ne_true_iff_false _).mp h_not_borderline
    simp only [propertyHolds]
    intro pos site h_pos h_site_pos h_ge
    by_cases h_cryptic : site.score ≥ crypticThreshold
    · have h_absurd := SpliceSiteScanner.scanner_completeness seq pos site
                        h_pos h_site_pos h_cryptic h_false_cryptic
      exact h_absurd
    · have h_absurd := SpliceSiteScanner.borderline_completeness seq pos site
                        h_pos h_site_pos h_ge h_cryptic h_false_borderline
      exact h_absurd

  | CodonAdapted organism threshold =>
    simp only [evaluate, propertyHolds] at *
    exact ite_fail_imp h_pass

  | GCInRange lo hi =>
    simp only [evaluate, propertyHolds] at *
    exact ite_fail_imp h_pass

  | NoRestrictionSite enzymeSites =>
    simp only [evaluate] at h_pass
    have h_not_true : hasAnyRestrictionSite seq enzymeSites ≠ true := by
      intro h_true; rw [if_pos h_true] at h_pass; cases h_pass
    have h_false : hasAnyRestrictionSite seq enzymeSites = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    simp only [propertyHolds]
    intro site pos h_site_mem h_pos h_match
    have h_has_true := hasAnyRestrictionSite_complete seq enzymeSites site pos
                         h_site_mem h_pos h_match
    rw [h_false] at h_has_true
    cases h_has_true

  | InFrame rf boundaries =>
    simp only [evaluate, propertyHolds] at *
    exact ite_fail_imp h_pass

  | NoInstabilityMotif =>
    simp only [evaluate] at h_pass
    have h_not_true : hasInstabilityMotif seq ≠ true := by
      intro h_true; rw [if_pos h_true] at h_pass; cases h_pass
    have h_false : hasInstabilityMotif seq = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    simp only [propertyHolds]
    constructor
    · intro pos h_pos h_match
      have h_has_attta : hasPattern seq atttaMotif = true :=
        hasPattern_complete seq atttaMotif pos h_pos h_match
      unfold hasInstabilityMotif at h_false
      have : hasPattern seq atttaMotif = false :=
        Bool.or_false_left (hasPattern seq atttaMotif) (hasPattern seq uRichMotif) h_false
      rw [this] at h_has_attta
      cases h_has_attta
    · intro pos h_pos h_match
      have h_has_urich : hasPattern seq uRichMotif = true :=
        hasPattern_complete seq uRichMotif pos h_pos h_match
      unfold hasInstabilityMotif at h_false
      have : hasPattern seq uRichMotif = false :=
        Bool.or_false_right (hasPattern seq atttaMotif) (hasPattern seq uRichMotif) h_false
      rw [this] at h_has_urich
      cases h_has_urich

  | NoCpGIsland =>
    simp only [evaluate] at h_pass
    have h_not_true : CpGIslandScanner.hasCpGIsland seq ≠ true := by
      intro h_true; rw [if_pos h_true] at h_pass; cases h_pass
    have h_false : CpGIslandScanner.hasCpGIsland seq = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    simp only [propertyHolds]
    intro pos h_pos
    have h_em := Classical.em
      (((((seq.drop pos).take cpgIslandWindowSize).count Nucleotide.G +
         ((seq.drop pos).take cpgIslandWindowSize).count Nucleotide.C : Rat) /
        ((seq.drop pos).take cpgIslandWindowSize).length) < cpgIslandGCThreshold)
    cases h_em with
    | inl h_gc_lt => left; exact h_gc_lt
    | inr h_gc_not_lt =>
      right
      have h_gc_ge := Rat.not_lt_iff_le _ _ |>.mp h_gc_not_lt
      have h_absurd := CpGIslandScanner.scanner_completeness seq pos h_pos h_gc_ge h_false
      exact h_absurd.elim

end BioCompiler
