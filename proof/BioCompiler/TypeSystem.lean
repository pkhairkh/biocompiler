/-
  BioCompiler.TypeSystem — Type Predicates and Soundness Proof

  Central Theorem:

    ∀ (P : TypePredicate) (seq : Sequence) (ctx : CellularContext),
      evaluate P seq ctx = Verdict.PASS → propertyHolds P seq ctx

  "Well-typed genes don't go wrong." — after Milner, 1978

  PREDICATE ARCHITECTURE (32 constructors):
  - Core (13): fully evaluable from sequence + context, can return PASS/FAIL/UNCERTAIN
    SpliceCorrect, NoCrypticSplice, CodonAdapted, GCInRange, NoRestrictionSite,
    InFrame, NoInstabilityMotif, NoCpGIsland, NoGTDinucleotide, NoStopCodons,
    ValidCodingSeq, CodonOptimality, NoCrypticPromoter
  - SLOT-dependent (19): require external tools/FFI, evaluate ALWAYS returns UNCERTAIN
    ConservationScore, NoUnexpectedTMDomain, mRNASecondaryStructure,
    CoTranslationalFolding, StructureConfidence, NoMisfoldingRisk,
    CorrectFoldTopology, NoUnexpectedInteraction, StableFolding,
    NoDestabilizingMutation, DisulfideBondIntegrity, HydrophobicCoreQuality,
    SolubleExpression, NoAggregationProneRegion, ChargeComposition,
    NoLongHydrophobicStretch, LowImmunogenicity, NoStrongTCellEpitope,
    NoDominantBCellEpitope, PopulationCoverageSafe

  SORRY STATUS: 0 remaining. All proofs are sorry-free.

  REFERENCE: DOC-03 (SDD) §3.5, DOC-10 (Deterministic Methods) §4
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners

namespace BioCompiler

open Verdict Sequence

-- ==============================================================================
-- Type Predicates (32 constructors: 13 core + 19 SLOT-dependent)
-- ==============================================================================

inductive TypePredicate where
  -- Core predicates (fully evaluable from sequence + context)
  | SpliceCorrect (cellType : String) : TypePredicate
  | NoCrypticSplice : TypePredicate
  | CodonAdapted (organism : String) (threshold : Rat) : TypePredicate
  | GCInRange (lo hi : Rat) : TypePredicate
  | NoRestrictionSite (enzymeSites : List Sequence) : TypePredicate
  | InFrame (readingFrame : Nat) (exonBoundaries : List Nat) : TypePredicate
  | NoInstabilityMotif : TypePredicate
  | NoCpGIsland : TypePredicate
  | NoGTDinucleotide : TypePredicate
  | NoStopCodons : TypePredicate
  | ValidCodingSeq : TypePredicate
  | CodonOptimality (organism : String) (threshold : Rat) : TypePredicate
  | NoCrypticPromoter (organism : String) (threshold : Rat) : TypePredicate
  -- SLOT-dependent predicates (evaluate ALWAYS returns UNCERTAIN)
  | ConservationScore (minScore : Int) : TypePredicate
  | NoUnexpectedTMDomain (isCytosolic : Bool) (threshold : Rat) : TypePredicate
  | mRNASecondaryStructure (dgThreshold : Rat) : TypePredicate
  | CoTranslationalFolding (organism : String) : TypePredicate
  | StructureConfidence (threshold : Rat) : TypePredicate
  | NoMisfoldingRisk : TypePredicate
  | CorrectFoldTopology : TypePredicate
  | NoUnexpectedInteraction : TypePredicate
  | StableFolding (ddgThreshold : Rat) : TypePredicate
  | NoDestabilizingMutation (maxDDG : Rat) : TypePredicate
  | DisulfideBondIntegrity : TypePredicate
  | HydrophobicCoreQuality (threshold : Rat) : TypePredicate
  | SolubleExpression (minScore : Rat) : TypePredicate
  | NoAggregationProneRegion : TypePredicate
  | ChargeComposition (pILo pIHi : Rat) : TypePredicate
  | NoLongHydrophobicStretch (maxLen : Nat) : TypePredicate
  | LowImmunogenicity (maxScore : Rat) : TypePredicate
  | NoStrongTCellEpitope (ic50Threshold : Rat) : TypePredicate
  | NoDominantBCellEpitope (scoreThreshold : Rat) : TypePredicate
  | PopulationCoverageSafe (maxCoverage : Rat) : TypePredicate
  deriving Repr

-- ==============================================================================
-- SLOT Classification
-- ==============================================================================

/-- Classify a type predicate as SLOT-dependent (true) or core (false).
    SLOT-dependent predicates rely on FFI output and cannot produce PASS. -/
def isSLOT : TypePredicate → Bool
  | TypePredicate.SpliceCorrect _ => false
  | TypePredicate.NoCrypticSplice => false
  | TypePredicate.CodonAdapted _ _ => false
  | TypePredicate.GCInRange _ _ => false
  | TypePredicate.NoRestrictionSite _ => false
  | TypePredicate.InFrame _ _ => false
  | TypePredicate.NoInstabilityMotif => false
  | TypePredicate.NoCpGIsland => false
  | TypePredicate.NoGTDinucleotide => false
  | TypePredicate.NoStopCodons => false
  | TypePredicate.ValidCodingSeq => false
  | TypePredicate.CodonOptimality _ _ => false
  | TypePredicate.NoCrypticPromoter _ _ => false
  | TypePredicate.ConservationScore _ => true
  | TypePredicate.NoUnexpectedTMDomain _ _ => true
  | TypePredicate.mRNASecondaryStructure _ => true
  | TypePredicate.CoTranslationalFolding _ => true
  | TypePredicate.StructureConfidence _ => true
  | TypePredicate.NoMisfoldingRisk => true
  | TypePredicate.CorrectFoldTopology => true
  | TypePredicate.NoUnexpectedInteraction => true
  | TypePredicate.StableFolding _ => true
  | TypePredicate.NoDestabilizingMutation _ => true
  | TypePredicate.DisulfideBondIntegrity => true
  | TypePredicate.HydrophobicCoreQuality _ => true
  | TypePredicate.SolubleExpression _ => true
  | TypePredicate.NoAggregationProneRegion => true
  | TypePredicate.ChargeComposition _ _ => true
  | TypePredicate.NoLongHydrophobicStretch _ => true
  | TypePredicate.LowImmunogenicity _ => true
  | TypePredicate.NoStrongTCellEpitope _ => true
  | TypePredicate.NoDominantBCellEpitope _ => true
  | TypePredicate.PopulationCoverageSafe _ => true

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
    [PromoterScanner] [TMDomainScanner] [mRNAStructureOracle] [CoTranslationalFoldingOracle]
    {State : Type} [DecidableEq State] [Inhabited State] [SplicingNDFST State] :
    TypePredicate → Sequence → CellularContext → Verdict
  -- Core predicates
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

  | TypePredicate.NoGTDinucleotide, seq, _ =>
      if hasPattern seq spliceDonorConsensus = true then FAIL else PASS

  | TypePredicate.NoStopCodons, seq, _ =>
      if hasPrematureStop seq 0 = true then FAIL else PASS

  | TypePredicate.ValidCodingSeq, seq, _ =>
      if isValidCodingSeq seq = true then PASS else FAIL

  | TypePredicate.CodonOptimality org threshold, seq, _ =>
      if CodonAdaptationIndex.computeCAI seq org ≥ threshold then PASS else FAIL

  | TypePredicate.NoCrypticPromoter organism threshold, seq, _ =>
      if PromoterScanner.hasCrypticPromoter seq organism threshold = true then FAIL
      else if PromoterScanner.hasBorderlinePromoter seq organism threshold = true then UNCERTAIN
      else PASS

  -- SLOT-dependent predicates: evaluate ALWAYS returns UNCERTAIN
  | TypePredicate.ConservationScore _, _, _ => UNCERTAIN
  | TypePredicate.NoUnexpectedTMDomain _ _, _, _ => UNCERTAIN
  | TypePredicate.mRNASecondaryStructure _, _, _ => UNCERTAIN
  | TypePredicate.CoTranslationalFolding _, _, _ => UNCERTAIN
  | TypePredicate.StructureConfidence _, _, _ => UNCERTAIN
  | TypePredicate.NoMisfoldingRisk, _, _ => UNCERTAIN
  | TypePredicate.CorrectFoldTopology, _, _ => UNCERTAIN
  | TypePredicate.NoUnexpectedInteraction, _, _ => UNCERTAIN
  | TypePredicate.StableFolding _, _, _ => UNCERTAIN
  | TypePredicate.NoDestabilizingMutation _, _, _ => UNCERTAIN
  | TypePredicate.DisulfideBondIntegrity, _, _ => UNCERTAIN
  | TypePredicate.HydrophobicCoreQuality _, _, _ => UNCERTAIN
  | TypePredicate.SolubleExpression _, _, _ => UNCERTAIN
  | TypePredicate.NoAggregationProneRegion, _, _ => UNCERTAIN
  | TypePredicate.ChargeComposition _ _, _, _ => UNCERTAIN
  | TypePredicate.NoLongHydrophobicStretch _, _, _ => UNCERTAIN
  | TypePredicate.LowImmunogenicity _, _, _ => UNCERTAIN
  | TypePredicate.NoStrongTCellEpitope _, _, _ => UNCERTAIN
  | TypePredicate.NoDominantBCellEpitope _, _, _ => UNCERTAIN
  | TypePredicate.PopulationCoverageSafe _, _, _ => UNCERTAIN

-- ==============================================================================
-- Property Semantics
-- ==============================================================================

def propertyHolds [SpliceSiteScanner] [CodonAdaptationIndex] [CpGIslandScanner]
    [PromoterScanner] [TMDomainScanner] [mRNAStructureOracle] [CoTranslationalFoldingOracle]
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
  | TypePredicate.NoGTDinucleotide, seq, _ =>
      ∀ (pos : Nat), pos + spliceDonorConsensus.length ≤ seq.length →
        (seq.drop pos).take spliceDonorConsensus.length ≠ spliceDonorConsensus
  | TypePredicate.NoStopCodons, seq, _ =>
      hasPrematureStop seq 0 = false
  | TypePredicate.ValidCodingSeq, seq, _ =>
      isValidCodingSeq seq = true
  | TypePredicate.CodonOptimality org threshold, seq, _ =>
      CodonAdaptationIndex.computeCAI seq org ≥ threshold
  | TypePredicate.NoCrypticPromoter organism threshold, seq, _ =>
      ∀ (pm : PromoterMatch), pm.organism = organism →
        pm.score ≥ threshold * 8 / 10 → False
  -- SLOT-dependent predicates: propertyHolds is True (vacuously) since evaluate
  -- never returns PASS
  | TypePredicate.ConservationScore _, _, _ => True
  | TypePredicate.NoUnexpectedTMDomain _ _, _, _ => True
  | TypePredicate.mRNASecondaryStructure _, _, _ => True
  | TypePredicate.CoTranslationalFolding _, _, _ => True
  | TypePredicate.StructureConfidence _, _, _ => True
  | TypePredicate.NoMisfoldingRisk, _, _ => True
  | TypePredicate.CorrectFoldTopology, _, _ => True
  | TypePredicate.NoUnexpectedInteraction, _, _ => True
  | TypePredicate.StableFolding _, _, _ => True
  | TypePredicate.NoDestabilizingMutation _, _, _ => True
  | TypePredicate.DisulfideBondIntegrity, _, _ => True
  | TypePredicate.HydrophobicCoreQuality _, _, _ => True
  | TypePredicate.SolubleExpression _, _, _ => True
  | TypePredicate.NoAggregationProneRegion, _, _ => True
  | TypePredicate.ChargeComposition _ _, _, _ => True
  | TypePredicate.NoLongHydrophobicStretch _, _, _ => True
  | TypePredicate.LowImmunogenicity _, _, _ => True
  | TypePredicate.NoStrongTCellEpitope _, _, _ => True
  | TypePredicate.NoDominantBCellEpitope _, _, _ => True
  | TypePredicate.PopulationCoverageSafe _, _, _ => True

-- ==============================================================================
-- Soundness Theorem
-- ==============================================================================

theorem type_soundness [SpliceSiteScanner] [CodonAdaptationIndex] [CpGIslandScanner]
    [PromoterScanner] [TMDomainScanner] [mRNAStructureOracle] [CoTranslationalFoldingOracle]
    {State : Type} [DecidableEq State] [Inhabited State] [SplicingNDFST State]
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    @evaluate _ _ _ _ _ _ _ State _ _ _ P seq ctx = PASS →
    @propertyHolds _ _ _ _ _ _ _ State _ _ _ P seq ctx := by
  intro h_pass
  cases P with

  | SpliceCorrect cellType =>
    simp only [evaluate] at h_pass
    cases h_cond : (ctx.cellType != cellType) with
    | true =>
      simp [h_cond] at h_pass
    | false =>
      have h_cell_eq : ctx.cellType = cellType := by
        have : ¬(ctx.cellType != cellType) := by
          intro h; rw [h] at h_cond; cases h_cond
        exact string_eq_of_not_ne _ _ this
      simp [h_cond] at h_pass
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

  | NoGTDinucleotide =>
    simp only [evaluate] at h_pass
    have h_not_true : hasPattern seq spliceDonorConsensus ≠ true := by
      intro h; rw [if_pos h] at h_pass; cases h_pass
    have h_false : hasPattern seq spliceDonorConsensus = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    simp only [propertyHolds]
    intro pos h_pos h_match
    have h_has := hasPattern_complete seq spliceDonorConsensus pos h_pos h_match
    rw [h_false] at h_has
    cases h_has

  | NoStopCodons =>
    simp only [evaluate] at h_pass
    have h_not_true : hasPrematureStop seq 0 ≠ true := by
      intro h; rw [if_pos h] at h_pass; cases h_pass
    have h_false : hasPrematureStop seq 0 = false :=
      (bool_ne_true_iff_false _).mp h_not_true
    simp only [propertyHolds]
    exact h_false

  | ValidCodingSeq =>
    simp only [evaluate, propertyHolds] at *
    exact ite_fail_imp h_pass

  | CodonOptimality organism threshold =>
    simp only [evaluate, propertyHolds] at *
    exact ite_fail_imp h_pass

  | NoCrypticPromoter organism threshold =>
    simp only [evaluate] at h_pass
    have h_not_cryptic : PromoterScanner.hasCrypticPromoter seq organism threshold ≠ true := by
      intro h; rw [if_pos h] at h_pass; cases h_pass
    have h_false_cryptic : PromoterScanner.hasCrypticPromoter seq organism threshold = false :=
      (bool_ne_true_iff_false _).mp h_not_cryptic
    have h_not_borderline : PromoterScanner.hasBorderlinePromoter seq organism threshold ≠ true := by
      intro h
      have : (if PromoterScanner.hasCrypticPromoter seq organism threshold = true then FAIL
              else if PromoterScanner.hasBorderlinePromoter seq organism threshold = true then UNCERTAIN
              else PASS) = UNCERTAIN := by
        rw [if_neg h_not_cryptic, if_pos h]
      rw [this] at h_pass; cases h_pass
    have h_false_borderline : PromoterScanner.hasBorderlinePromoter seq organism threshold = false :=
      (bool_ne_true_iff_false _).mp h_not_borderline
    simp only [propertyHolds]
    intro pm h_org h_ge
    by_cases h_above : pm.score ≥ threshold
    · have h_absurd := PromoterScanner.scanner_completeness seq organism threshold pm h_org h_above h_false_cryptic
      exact h_absurd
    · have h_absurd := PromoterScanner.borderline_completeness seq organism threshold pm h_org h_ge h_above h_false_borderline
      exact h_absurd

  -- SLOT-dependent predicates: evaluate returns UNCERTAIN ≠ PASS, so vacuously true
  | ConservationScore _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoUnexpectedTMDomain _ _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | mRNASecondaryStructure _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | CoTranslationalFolding _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | StructureConfidence _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoMisfoldingRisk =>
    simp only [evaluate] at h_pass; cases h_pass
  | CorrectFoldTopology =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoUnexpectedInteraction =>
    simp only [evaluate] at h_pass; cases h_pass
  | StableFolding _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoDestabilizingMutation _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | DisulfideBondIntegrity =>
    simp only [evaluate] at h_pass; cases h_pass
  | HydrophobicCoreQuality _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | SolubleExpression _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoAggregationProneRegion =>
    simp only [evaluate] at h_pass; cases h_pass
  | ChargeComposition _ _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoLongHydrophobicStretch _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | LowImmunogenicity _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoStrongTCellEpitope _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | NoDominantBCellEpitope _ =>
    simp only [evaluate] at h_pass; cases h_pass
  | PopulationCoverageSafe _ =>
    simp only [evaluate] at h_pass; cases h_pass

end BioCompiler
