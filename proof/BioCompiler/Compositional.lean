/-
  BioCompiler.Compositional — Compositional Soundness Proof

  This module proves:
  1. compositional_soundness: evaluateAll = PASS → ∀ P, propertyHolds P
  2. Constraints DON'T compose: concrete counterexamples for dinucleotides
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
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    Verdict :=
  (predicates.map (fun P =>
    @evaluate inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst P seq ctx
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
    {State : Type} [inst_dec : DecidableEq State] [inst_inhab : Inhabited State] [inst_ndfst : SplicingNDFST State]
    (predicates : List TypePredicate) (seq : Sequence) (ctx : CellularContext) :
    @evaluateAll inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst predicates seq ctx = PASS →
    ∀ P ∈ predicates, @propertyHolds inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst P seq ctx := by
  intro h P hP
  have h_all_pass : ∀ v ∈ (predicates.map (fun P =>
      @evaluate inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst P seq ctx)), v = PASS :=
    Verdict.foldl_and_pass_implies_all_pass _ h
  have h_eval_pass : @evaluate inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst P seq ctx = PASS := by
    have : @evaluate inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst P seq ctx ∈
        predicates.map (fun P =>
          @evaluate inst_splice inst_cai inst_cpg State inst_dec inst_inhab inst_ndfst P seq ctx) := by
      simp [List.mem_map]
      exact ⟨P, hP, rfl⟩
    exact h_all_pass _ this
  exact type_soundness P seq ctx h_eval_pass

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
