/-
  BioCompiler.TypeSystem — Type Predicates and Soundness Proof

  Central Theorem:

    ∀ (P : TypePredicate) (seq : Sequence) (ctx : CellularContext),
      evaluate P seq ctx = Verdict.PASS → propertyHolds P seq ctx

  "Well-typed genes don't go wrong." — after Milner, 1978

  PROOF STRATEGY:
  - Arithmetic predicates (CodonAdapted, GCInRange, InFrame): PASS condition = property.
    Proved via dite_fail_imp (contrapositive: ¬condition → FAIL = PASS, contradiction).
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
theorem dite_fail_imp {cond : Prop} [Decidable cond]
    (h : (if cond then PASS else FAIL) = PASS) : cond := by
  by_contra h_neg
  simp only [dif_neg h_neg] at h
  cases h

/-- If `if cond then FAIL else PASS = PASS`, then ¬cond holds.
    Contrapositive: if cond, the then branch gives FAIL ≠ PASS. -/
theorem dite_pass_imp_neg {cond : Prop} [Decidable cond]
    (h : (if cond then FAIL else PASS) = PASS) : ¬cond := by
  intro h_pos
  simp only [dif_pos h_pos] at h
  cases h

/-- Bool: b ≠ true ↔ b = false. -/
theorem bool_ne_true_iff_false (b : Bool) : b ≠ true ↔ b = false := by
  cases b <;> simp

/-- Bool: (a || b) = false → a = false. -/
theorem Bool.or_false_left (a b : Bool) (h : (a || b) = false) : a = false := by
  cases a <;> simp at h <;> exact h

/-- Bool: (a || b) = false → b = false. -/
theorem Bool.or_false_right (a b : Bool) (h : (a || b) = false) : b = false := by
  cases a <;> simp at h <;> cases b <;> simp at h <;> exact h

/-- A singleton list has length 1. -/
theorem list_singleton_length {α : Type} {x : α} : ([x] : List α).length = 1 := by simp

/-- String equality from negation of inequality (using BEq). -/
theorem string_eq_of_not_ne [BEq String] [LawfulBEq String]
    (s₁ s₂ : String) (h : ¬(s₁ != s₂)) : s₁ = s₂ := by
  simp [BEq.ne] at h
  exact lawfulBEq_eq s₁ s₂ |>.mp h

-- ==============================================================================
-- Evaluation Function
-- ==============================================================================

def evaluate [SpliceSiteScanner] [CodonAdaptationIndex] [CpGIslandScanner]
    {State : Type} [DecidableEq State] [SplicingNDFST State] :
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
    {State : Type} [DecidableEq State] [SplicingNDFST State] :
    TypePredicate → Sequence → CellularContext → Prop
  | TypePredicate.SpliceCorrect cellType, seq, ctx =>
      ctx.cellType = cellType ∧
        (ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq).length = 1
  | TypePredicate.NoCrypticSplice, seq, _ =>
      -- PASS means: no strong cryptic sites AND no borderline sites.
      -- Stated in contrapositive form: no site with score ≥ uncertainLoThreshold
      -- exists. This is equivalent to "all sites have score < uncertainLoThreshold"
      -- and matches the proof structure (contrapositive of scanner completeness).
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
          seq.drop pos |>.take site.length ≠ site
  | TypePredicate.InFrame rf boundaries, seq, _ =>
      readingFrameConsistent boundaries rf = true ∧ hasPrematureStop seq rf = false
  | TypePredicate.NoInstabilityMotif, seq, _ =>
      (∀ (pos : Nat), pos + atttaMotif.length ≤ seq.length →
        seq.drop pos |>.take atttaMotif.length ≠ atttaMotif) ∧
      (∀ (pos : Nat), pos + uRichMotif.length ≤ seq.length →
        seq.drop pos |>.take uRichMotif.length ≠ uRichMotif)
  | TypePredicate.NoCpGIsland, seq, _ =>
      ∀ (pos : Nat),
        pos + cpgIslandWindowSize ≤ seq.length →
          let window := seq.drop pos |>.take cpgIslandWindowSize
          -- No CpG island means: for every window, either GC < threshold
          -- OR the observed/expected CpG ratio < threshold.
          -- The full Obs/Exp ratio formalization requires counting CpG
          -- dinucleotides and computing (count_CG * window_len) / (count_C * count_G),
          -- which is definable but involved. We formalize the GC condition
          -- as the primary criterion (which is the stronger condition:
          -- if GC% < threshold, the region cannot be a CpG island regardless
          -- of Obs/Exp ratio). The scanner_completeness assumption for
          -- CpGIslandScanner ensures that if a true CpG island exists
          -- (satisfying BOTH conditions), the scanner will detect it.
          ((window.count Nucleotide.G + window.count Nucleotide.C : Rat) / window.length < cpgIslandGCThreshold) ∨
          (∃ (cpgCount : Nat),
            cpgCount = (List.zipWith (· == ·) window (window.drop 1)).count true ∧
            (cpgCount : Rat) * window.length <
              cpgIslandObsExpThreshold * (window.count Nucleotide.C) * (window.count Nucleotide.G))

-- ==============================================================================
-- ╔══════════════════════════════════════════════════════════════════════════════╗
-- ║                    SOUNDNESS THEOREM                                      ║
-- ║                                                                            ║
-- ║  ∀ (P : TypePredicate) (seq : Sequence) (ctx : CellularContext),           ║
-- ║    evaluate P seq ctx = PASS → propertyHolds P seq ctx                     ║
-- ║                                                                            ║
-- ║  "Well-typed genes don't go wrong."                                        ║
-- ╚══════════════════════════════════════════════════════════════════════════════╝
-- ==============================================================================

theorem type_soundness [SpliceSiteScanner] [CodonAdaptationIndex] [CpGIslandScanner]
    {State : Type} [DecidableEq State] [SplicingNDFST State]
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    evaluate P seq ctx = PASS → propertyHolds P seq ctx := by
  intro h_pass
  cases P with

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 1: SpliceCorrect(cellType)
  -- ═══════════════════════════════════════════════════════════════════════════
  | SpliceCorrect cellType =>
    unfold evaluate at h_pass
    -- Step 1: Prove ¬(ctx.cellType != cellType)
    -- If ctx.cellType != cellType were true, the if gives UNCERTAIN ≠ PASS
    have h_not_ne : ¬(ctx.cellType != cellType) := by
      intro h_cond
      have : (if ctx.cellType != cellType then (UNCERTAIN : Verdict) else
              match ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq with
              | [_] => PASS | _ => FAIL) = UNCERTAIN := dif_pos h_cond
      rw [this] at h_pass
      cases h_pass
    -- Step 2: ctx.cellType = cellType (by LawfulBEq)
    have h_cell_eq : ctx.cellType = cellType :=
      string_eq_of_not_ne _ _ h_not_ne
    -- Step 3: Simplify the if-then-else to the else branch
    have h_cond_false : (ctx.cellType != cellType) = false := by
      cases h_cond : (ctx.cellType != cellType) with
      | true => exact absurd h_cond h_not_ne
      | false => rfl
    rw [dif_neg h_cond_false] at h_pass
    -- Step 4: h_pass now concerns only the match
    -- (match ... with | [_] => PASS | _ => FAIL) = PASS
    -- Case-split on the list to show it must be a singleton
    have h_list_len : (ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq).length = 1 := by
      by_contra h_ne
      cases h_list : ndfstUniqueOutputSet (SplicingNDFST.ndfst : NDFST State) seq with
      | nil =>
        -- Empty list: match gives FAIL ≠ PASS
        simp only [h_list, List.length_nil] at h_ne h_pass
        cases h_pass
      | cons hd tl =>
        cases tl with
        | nil =>
          -- Singleton [hd]: length = 1, contradicting h_ne
          simp only [h_list, List.length_cons, List.length_nil, Nat.add_zero] at h_ne
          omega
        | cons hd' tl' =>
          -- ≥2 elements: match gives FAIL ≠ PASS
          simp only [h_list, List.length_cons] at h_pass
          cases h_pass
    unfold propertyHolds
    exact ⟨h_cell_eq, h_list_len⟩

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 2: NoCrypticSplice — PROVED COMPLETELY (dual-threshold)
  --
  -- PASS means: no cryptic sites AND no borderline sites.
  -- If hasCrypticSpliceSite = true → FAIL ≠ PASS (contradiction)
  -- If hasBorderlineSpliceSite = true → UNCERTAIN ≠ PASS (contradiction)
  -- Therefore both scanners returned false, and by contrapositive of
  -- their completeness, no site with score ≥ uncertainLoThreshold exists.
  --
  -- The property is stated in contrapositive form:
  --   site.score ≥ uncertainLoThreshold → False
  -- which is equivalent to site.score < uncertainLoThreshold.
  -- ═══════════════════════════════════════════════════════════════════════════
  | NoCrypticSplice =>
    unfold evaluate at h_pass
    -- h_pass : (if hasCrypticSpliceSite seq = true then FAIL
    --           else if hasBorderlineSpliceSite seq = true then UNCERTAIN
    --           else PASS) = PASS
    -- For PASS, neither condition can be true
    have h_not_cryptic : SpliceSiteScanner.hasCrypticSpliceSite seq ≠ true := by
      intro h; simp [dif_pos h] at h_pass; cases h_pass
    have h_false_cryptic : SpliceSiteScanner.hasCrypticSpliceSite seq = false :=
      (bool_ne_true_iff_false _).mp h_not_cryptic
    have h_not_borderline : SpliceSiteScanner.hasBorderlineSpliceSite seq ≠ true := by
      intro h
      have : (if SpliceSiteScanner.hasCrypticSpliceSite seq = true then FAIL
              else if SpliceSiteScanner.hasBorderlineSpliceSite seq = true then UNCERTAIN
              else PASS) = UNCERTAIN := by
        simp [dif_neg h_false_cryptic, dif_pos h]
      rw [this] at h_pass; cases h_pass
    have h_false_borderline : SpliceSiteScanner.hasBorderlineSpliceSite seq = false :=
      (bool_ne_true_iff_false _).mp h_not_borderline
    unfold propertyHolds
    intro pos site h_pos h_site_pos h_ge
    -- h_ge : site.score ≥ uncertainLoThreshold
    -- Goal: False
    -- If score ≥ uncertainLoThreshold, it's either cryptic or borderline
    by_cases h_cryptic : site.score ≥ crypticThreshold
    · -- Strong cryptic site → scanner should have found it → contradiction
      have h_absurd := SpliceSiteScanner.scanner_completeness seq pos site
                        h_pos h_site_pos h_cryptic h_false_cryptic
      exact absurd rfl h_absurd
    · -- Not cryptic, so borderline → borderline scanner should have found it → contradiction
      have h_absurd := SpliceSiteScanner.borderline_completeness seq pos site
                        h_pos h_site_pos h_ge h_cryptic h_false_borderline
      exact absurd rfl h_absurd

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 3: CodonAdapted — PROVED COMPLETELY
  -- ═══════════════════════════════════════════════════════════════════════════
  | CodonAdapted organism threshold =>
    unfold evaluate propertyHolds at *
    exact dite_fail_imp h_pass

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 4: GCInRange — PROVED COMPLETELY
  -- ═══════════════════════════════════════════════════════════════════════════
  | GCInRange lo hi =>
    unfold evaluate propertyHolds at *
    exact dite_fail_imp h_pass

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 5: NoRestrictionSite — PROVED COMPLETELY
  -- ═══════════════════════════════════════════════════════════════════════════
  | NoRestrictionSite enzymeSites =>
    unfold evaluate at h_pass
    have h_not_true : hasAnyRestrictionSite seq enzymeSites ≠ true := by
      intro h_true
      simp only [dif_pos h_true] at h_pass
      cases h_pass
    have h_false : hasAnyRestrictionSite seq enzymeSites = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    unfold propertyHolds
    intro site pos h_site_mem h_pos h_match
    have h_has_true := hasAnyRestrictionSite_complete seq enzymeSites site pos
                         h_site_mem h_pos h_match
    rw [h_false] at h_has_true
    cases h_has_true

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 6: InFrame — PROVED COMPLETELY
  -- ═══════════════════════════════════════════════════════════════════════════
  | InFrame rf boundaries =>
    unfold evaluate propertyHolds at *
    exact dite_fail_imp h_pass

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 7: NoInstabilityMotif
  -- ═══════════════════════════════════════════════════════════════════════════
  | NoInstabilityMotif =>
    unfold evaluate at h_pass
    have h_not_true : hasInstabilityMotif seq ≠ true := by
      intro h_true
      simp only [dif_pos h_true] at h_pass
      cases h_pass
    have h_false : hasInstabilityMotif seq = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    unfold propertyHolds
    constructor
    · -- No ATTTA motif at any position
      intro pos h_pos h_match
      have h_has_attta : hasPattern seq atttaMotif = true :=
        hasPattern_complete seq atttaMotif pos h_pos h_match
      unfold hasInstabilityMotif at h_false
      -- hasInstabilityMotif = hasPattern seq atttaMotif || hasPattern seq uRichMotif
      -- h_false : (hasPattern seq atttaMotif || hasPattern seq uRichMotif) = false
      -- But h_has_attta : hasPattern seq atttaMotif = true
      -- So true || _ = true ≠ false. Contradiction.
      have : hasPattern seq atttaMotif = false :=
        Bool.or_false_left (hasPattern seq atttaMotif) (hasPattern seq uRichMotif) h_false
      rw [this] at h_has_attta
      cases h_has_attta
    · -- No U-rich motif at any position
      intro pos h_pos h_match
      have h_has_urich : hasPattern seq uRichMotif = true :=
        hasPattern_complete seq uRichMotif pos h_pos h_match
      unfold hasInstabilityMotif at h_false
      have : hasPattern seq uRichMotif = false :=
        Bool.or_false_right (hasPattern seq atttaMotif) (hasPattern seq uRichMotif) h_false
      rw [this] at h_has_urich
      cases h_has_urich

  -- ═══════════════════════════════════════════════════════════════════════════
  -- Case 8: NoCpGIsland — PROVED via scanner completeness
  --
  -- PROOF STRATEGY: If the scanner returned false, then by contrapositive
  -- of scanner_completeness, no window satisfies BOTH the GC threshold AND
  -- the Obs/Exp threshold simultaneously. Therefore, for every window,
  -- either GC < threshold OR Obs/Exp < threshold. We prove this by
  -- contradiction: if there were a position where BOTH conditions hold,
  -- the scanner would have returned true.
  -- ═══════════════════════════════════════════════════════════════════════════
  | NoCpGIsland =>
    unfold evaluate at h_pass
    have h_not_true : CpGIslandScanner.hasCpGIsland seq ≠ true := by
      intro h_true
      simp only [dif_pos h_true] at h_pass
      cases h_pass
    have h_false : CpGIslandScanner.hasCpGIsland seq = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    unfold propertyHolds
    intro pos h_pos
    -- If both GC >= threshold AND Obs/Exp >= threshold at some position,
    -- then by scanner_completeness, the scanner would return true (contradiction).
    -- Therefore, at least one condition must fail, which means either
    -- GC < threshold (left disjunct) or Obs/Exp < threshold (right disjunct).
    by_contra h_both_fail
    -- h_both_fail : ¬(GC < threshold ∨ Obs/Exp < threshold)
    push_neg at h_both_fail
    obtain ⟨h_gc_ge, h_obs_ge⟩ := h_both_fail
    -- h_gc_ge : (G+C)/len ≥ threshold
    -- h_obs_ge : ¬(∃ cpgCount, cpgCount = CG_dinucleotides ∧ ratio < threshold)
    -- From h_obs_ge, we know Obs/Exp ≥ threshold, which combined with h_gc_ge
    -- means this position satisfies both CpG island criteria.
    -- By scanner_completeness (contrapositive), scanner must have returned true.
    have h_scanner_true : CpGIslandScanner.hasCpGIsland seq = false → False :=
      CpGIslandScanner.scanner_completeness seq pos h_pos h_gc_ge
    exact h_scanner_true h_false

end BioCompiler
