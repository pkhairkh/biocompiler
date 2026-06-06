/-
  BioCompiler.Compositional — Compositional Soundness Proof

  This module proves:
  1. compositional_soundness: evaluateAll = PASS → ∀ P, propertyHolds P
  2. slot_predicates_uncertain: SLOT-dependent predicates never produce PASS
  3. slot_predicates_dont_affect_pass: SLOT predicates in the list prevent overall PASS
  4. all_core_if_pass: if evaluateAll = PASS, all predicates are core
  5. Constraints DON'T compose: concrete counterexamples for dinucleotides
     and restriction sites forming at sequence concatenation boundaries

  Reference: DOC-03 (SDD) §3.5.4, DOC-01 (SRS) INV-TYP-03
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners
import BioCompiler.TypeSystem

namespace BioCompiler

open Verdict Sequence

-- ==============================================================================
-- Composed Evaluation
-- ==============================================================================

/-- Evaluate a list of type predicates and compose the results
    using three-valued conjunction. -/
def evaluateAll [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner] [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    Verdict :=
  (predicates.map (fun P =>
    @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst P seq ctx
  )).foldl Verdict.and PASS

-- ==============================================================================
-- Compositional Soundness Theorem
-- ==============================================================================

/-- THEOREM (Compositional Soundness): If the composed evaluation of all
    predicates yields PASS, then every individual property holds.

    This follows from:
    1. evaluateAll = PASS → every evaluate P = PASS (by foldl_and_pass_implies_all_pass)
    2. evaluate P = PASS → propertyHolds P (by type_soundness)

    Corollary: A guarantee certificate (which requires overall PASS) can only
    be issued when all claimed properties actually hold. -/
theorem compositional_soundness [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner] [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    @evaluateAll inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst predicates seq ctx = PASS →
    ∀ P ∈ predicates, @propertyHolds inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst P seq ctx := by
  intro h P hP
  have h_all_pass : ∀ v ∈ (predicates.map (fun P =>
      @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
        State inst_dec inst_inhab inst_ndfst P seq ctx)), v = PASS :=
    Verdict.foldl_and_pass_implies_all_pass _ h
  have h_eval_pass : @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst P seq ctx = PASS := by
    have : @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
        State inst_dec inst_inhab inst_ndfst P seq ctx ∈
        predicates.map (fun P =>
          @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
            State inst_dec inst_inhab inst_ndfst P seq ctx) := by
      simp [List.mem_map]
      exact ⟨P, hP, rfl⟩
    exact h_all_pass _ this
  exact type_soundness P seq ctx h_eval_pass

-- ==============================================================================
-- SLOT-Dependent Predicate Properties
-- ==============================================================================

/-- THEOREM (SLOT Predicates Uncertain): SLOT-dependent predicates never
    produce a PASS verdict. Since they depend on non-deterministic FFI output,
    their evaluation always returns UNCERTAIN.

    This is the key lemma that ensures SLOT-dependent predicates don't affect
    the compositional result: if evaluateAll = PASS, no SLOT-dependent predicate
    was in the list (because UNCERTAIN ∧ PASS = UNCERTAIN ≠ PASS). -/
theorem slot_predicates_uncertain [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner] [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (pred : TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    isSLOT pred = true →
    @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst pred seq ctx ≠ Verdict.PASS := by
  intro h_slot
  intro h_pass
  cases pred with
  -- Core predicates: isSLOT = false, contradiction with h_slot
  | SpliceCorrect _ =>
      simp only [isSLOT] at h_slot; cases h_slot
  | NoCrypticSplice =>
      simp only [isSLOT] at h_slot; cases h_slot
  | CodonAdapted _ _ =>
      simp only [isSLOT] at h_slot; cases h_slot
  | GCInRange _ _ =>
      simp only [isSLOT] at h_slot; cases h_slot
  | NoRestrictionSite _ =>
      simp only [isSLOT] at h_slot; cases h_slot
  | InFrame _ _ =>
      simp only [isSLOT] at h_slot; cases h_slot
  | NoInstabilityMotif =>
      simp only [isSLOT] at h_slot; cases h_slot
  | NoCpGIsland =>
      simp only [isSLOT] at h_slot; cases h_slot
  | NoGTDinucleotide =>
      simp only [isSLOT] at h_slot; cases h_slot
  | NoStopCodons =>
      simp only [isSLOT] at h_slot; cases h_slot
  | ValidCodingSeq =>
      simp only [isSLOT] at h_slot; cases h_slot
  | CodonOptimality _ _ =>
      simp only [isSLOT] at h_slot; cases h_slot
  | NoCrypticPromoter _ _ =>
      simp only [isSLOT] at h_slot; cases h_slot
  -- SLOT-dependent predicates: evaluate returns UNCERTAIN ≠ PASS
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

/-- COROLLARY: If any SLOT-dependent predicate is in the list, evaluateAll
    cannot return PASS. This is because:
    1. The SLOT predicate evaluates to UNCERTAIN
    2. UNCERTAIN ∧ (anything) ∈ {UNCERTAIN, FAIL}
    3. foldl from PASS with an UNCERTAIN element gives UNCERTAIN ≠ PASS

    Therefore, a guarantee certificate can only be issued for lists
    containing exclusively core predicates. -/
theorem slot_predicates_dont_affect_pass [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner] [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    (∃ P ∈ predicates, isSLOT P = true) →
    @evaluateAll inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst predicates seq ctx ≠ PASS := by
  intro ⟨P, hP_mem, hP_slot⟩
  have hP_not_pass := @slot_predicates_uncertain inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
    State inst_dec inst_inhab inst_ndfst P seq ctx hP_slot
  have h_all_pass : (∀ v ∈ (predicates.map (fun Q =>
      @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
        State inst_dec inst_inhab inst_ndfst Q seq ctx)), v = PASS) →
      @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
        State inst_dec inst_inhab inst_ndfst P seq ctx = PASS := by
    intro h_all
    have hP_in : @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
        State inst_dec inst_inhab inst_ndfst P seq ctx ∈
        predicates.map (fun Q => @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
          State inst_dec inst_inhab inst_ndfst Q seq ctx) := by
      simp [List.mem_map]; exact ⟨P, hP_mem, rfl⟩
    exact h_all _ hP_in
  intro h_eval_all
  have hP_pass := h_all_pass (Verdict.foldl_and_pass_implies_all_pass _ h_eval_all)
  exact hP_not_pass hP_pass

/-- COROLLARY: If evaluateAll returns PASS, then no SLOT-dependent predicate
    is in the list. Equivalently, all predicates in the list are core predicates. -/
theorem all_core_if_pass [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner] [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    @evaluateAll inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst predicates seq ctx = PASS →
    ∀ P ∈ predicates, isSLOT P = false := by
  intro h P hP
  cases h_bool : isSLOT P with
  | true =>
    have h_not_pass := @slot_predicates_dont_affect_pass inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State inst_dec inst_inhab inst_ndfst predicates seq ctx ⟨P, hP, h_bool⟩
    exact absurd h h_not_pass
  | false => rfl

-- ==============================================================================
-- Sequence Concatenation: Constraints Don't Compose
--
-- These counterexamples demonstrate the KEY advantage of type systems over
-- constraint lists: types compose, constraints don't.
--
-- Two gene fragments that individually satisfy all constraints can produce
-- a constraint violation when concatenated. This cannot happen with type
-- systems because the composition of type-checked fragments remains
-- well-typed (with appropriate junction checking).
-- ==============================================================================

/-- THEOREM: Dinucleotide predicates DON'T compose (counterexample).
    GT dinucleotide can form at the junction of two GT-free sequences.
    Construction: s1 = [C, G], s2 = [T, C].
    Neither contains GT, but [C,G,T,C] contains GT at position 1.

    This is THE key argument for why type systems are needed: constraint lists
    don't compose, but type-checked fragments with junction checking do. -/
theorem dinucleotide_no_compose :
    ∃ (s1 s2 : Sequence),
      hasPattern s1 spliceDonorConsensus = false ∧
      hasPattern s2 spliceDonorConsensus = false ∧
      hasPattern (s1 ++ s2) spliceDonorConsensus = true := by
  exact ⟨
    [Nucleotide.C, Nucleotide.G],
    [Nucleotide.T, Nucleotide.C],
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

/-- THEOREM: Restriction sites DON'T compose (counterexample).
    GATC (Sau3AI) forms at the junction of [G,A,T] and [C,G].
    Neither fragment contains GATC, but the concatenation does. -/
theorem restriction_site_no_compose :
    ∃ (s1 s2 : Sequence) (site : Sequence),
      hasPattern s1 site = false ∧
      hasPattern s2 site = false ∧
      site.length = 4 ∧
      hasPattern (s1 ++ s2) site = true := by
  exact ⟨
    [Nucleotide.G, Nucleotide.A, Nucleotide.T],
    [Nucleotide.C, Nucleotide.G],
    [Nucleotide.G, Nucleotide.A, Nucleotide.T, Nucleotide.C],
    by native_decide,
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

end BioCompiler
