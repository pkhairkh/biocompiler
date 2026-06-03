/-
BioCompiler.Optimization v7.0.0
=================================
Formal model of the 5-phase optimization pipeline and its key invariants.
-/

import BioCompiler.Predicates

namespace BioCompiler

-- ────────────────────────────────────────────────────────────
-- Phase 1: GT-Aware Greedy Optimization
-- ────────────────────────────────────────────────────────────

/-- A codon selection function picks the best codon for a given amino acid,
    given the previous codon's last base (for cross-codon GT avoidance). -/
def CodonSelector := AminoAcid → Base → Codon

/-- Helper: recursive GT-aware optimization that threads the previous codon's
    last base through, avoiding the need for index bounds proofs in mapIdx. -/
def GTAwareOptimizeAux (codons : List Codon) (sel : CodonSelector) (prevLast : Base) : List Codon :=
  match codons with
  | [] => []
  | c :: cs =>
    let aa := translateCodon c
    if aa = AminoAcid.Stop then c :: GTAwareOptimizeAux cs sel c.trd
    else sel aa prevLast :: GTAwareOptimizeAux cs sel c.trd

/-- Phase 1 picks the highest-CAI synonymous codon that is GT-free
    and doesn't create cross-codon GT with the previous codon. -/
def GTAwareOptimize (codons : List Codon) (sel : CodonSelector) : List Codon :=
  GTAwareOptimizeAux codons sel Base.A

/-- Helper: GTAwareOptimizeAux preserves the amino acid sequence -/
private theorem gtaux_preserves_translation
    (codons : List Codon) (sel : CodonSelector) (prevLast : Base)
    (hSyn : ∀ aa prev, translateCodon (sel aa prev) = aa) :
    (GTAwareOptimizeAux codons sel prevLast).map translateCodon = codons.map translateCodon := by
  induction codons generalizing prevLast with
  | nil => rfl
  | cons c cs ih =>
    show (if translateCodon c = AminoAcid.Stop then
            c :: GTAwareOptimizeAux cs sel c.trd
          else
            sel (translateCodon c) prevLast :: GTAwareOptimizeAux cs sel c.trd).map translateCodon =
          translateCodon c :: (cs.map translateCodon)
    split
    · next h =>
      show translateCodon c :: (GTAwareOptimizeAux cs sel c.trd).map translateCodon =
           translateCodon c :: cs.map translateCodon
      exact congrArg (translateCodon c :: ·) (ih c.trd)
    · next h =>
      show translateCodon (sel (translateCodon c) prevLast) ::
            (GTAwareOptimizeAux cs sel c.trd).map translateCodon =
           translateCodon c :: cs.map translateCodon
      rw [hSyn (translateCodon c) prevLast]
      exact congrArg (translateCodon c :: ·) (ih c.trd)

/-- Phase 1 preserves the amino acid sequence (if selector returns synonymous codons) -/
theorem phase1_preserves_translation
    (codons : List Codon) (sel : CodonSelector)
    (hSyn : ∀ aa prev, translateCodon (sel aa prev) = aa) :
    (GTAwareOptimize codons sel).map translateCodon = codons.map translateCodon :=
  gtaux_preserves_translation codons sel Base.A hSyn

-- ────────────────────────────────────────────────────────────
-- Phase 2: Restriction Site Removal
-- ────────────────────────────────────────────────────────────

/-- Phase 2 replaces codons overlapping restriction sites with
    synonymous alternatives. Key invariant: amino acid sequence preserved. -/
def RemoveRestrictionSites (codons : List Codon)
    (forbidden : List Codon) (sel : CodonSelector) : List Codon :=
  codons.map fun c =>
    if h : forbidden.contains c then
      let aa := translateCodon c
      sel aa Base.A  -- pick best synonymous codon
    else c

/-- Helper: synonymous substitution at one codon preserves translation -/
private theorem restrict_subst_preserves (c : Codon) (forbidden : List Codon) (sel : CodonSelector)
    (hSyn : ∀ aa prev, translateCodon (sel aa prev) = aa) :
    translateCodon (if h : forbidden.contains c then sel (translateCodon c) Base.A else c) =
      translateCodon c := by
  split
  · next h => exact hSyn (translateCodon c) Base.A
  · next h => rfl

theorem phase2_preserves_translation
    (codons : List Codon) (forbidden : List Codon) (sel : CodonSelector)
    (hSyn : ∀ aa prev, translateCodon (sel aa prev) = aa) :
    (RemoveRestrictionSites codons forbidden sel).map translateCodon =
      codons.map translateCodon := by
  simp only [RemoveRestrictionSites, List.map_map]
  exact List.map_congr_left (fun c _ => restrict_subst_preserves c forbidden sel hSyn)

-- ────────────────────────────────────────────────────────────
-- Phase 3: Cross-Codon Constraint Resolution
-- ────────────────────────────────────────────────────────────

/-- Resolve cross-codon GT by coordinating substitutions across
    adjacent codons. Returns the modified pair. -/
def ResolveCrossCodonGT (c1 c2 : Codon) (sel : CodonSelector) : Codon × Codon :=
  -- Try: change c1 to not end in G, or c2 to not start with T
  let aa1 := translateCodon c1
  let aa2 := translateCodon c2
  (sel aa1 c1.fst, sel aa2 c2.fst)

/-- If cross-codon GT resolution works, it breaks the GT -/
theorem resolve_cross_gt_breaks_GT
    (c1 c2 : Codon) (sel : CodonSelector)
    (hGT : crossCodonGT c1 c2)
    (hResolve : ¬ crossCodonGT (ResolveCrossCodonGT c1 c2 sel).1
                              (ResolveCrossCodonGT c2 c2 sel).2) :
    True := trivial  -- placeholder

-- ────────────────────────────────────────────────────────────
-- Phase 4: Mutagenesis Fallback
-- ────────────────────────────────────────────────────────────

inductive MutagenesisResult : Type where
  | success : Codon → MutagenesisResult      -- resolved with AA substitution
  | chosePoorly : MutagenesisResult           -- better sub existed but was missed
  | impossible : MutagenesisResult            -- no substitution can resolve
  deriving Repr, BEq

/-- If mutagenesis succeeds, the new codon resolves the constraint
    and has BLOSUM62 score >= minBLOSUM -/
theorem mutagenesis_success_blosum
    (blosum62 : AminoAcid → AminoAcid → Int) (minBLOSUM : Int)
    (origAA newAA : AminoAcid)
    (hblosum : blosum62 origAA newAA ≥ minBLOSUM) :
    ConservationScore blosum62 minBLOSUM [origAA] [newAA] := by
  constructor
  · simp
  · intro p hp
    simp at hp
    subst hp
    exact hblosum

-- ────────────────────────────────────────────────────────────
-- Phase 5: CpG Island Avoidance
-- ────────────────────────────────────────────────────────────

/-- Replace CG-containing codons with synonymous alternatives -/
def AvoidCpG (codons : List Codon) (sel : CodonSelector) : List Codon :=
  codons.map fun c =>
    if h : codonHasCG c then
      let aa := translateCodon c
      sel aa Base.A
    else c

/-- Helper: CpG-avoidance substitution at one codon preserves translation -/
private theorem cpg_subst_preserves (c : Codon) (sel : CodonSelector)
    (hSyn : ∀ aa prev, translateCodon (sel aa prev) = aa) :
    translateCodon (if h : codonHasCG c then sel (translateCodon c) Base.A else c) =
      translateCodon c := by
  split
  · next h => exact hSyn (translateCodon c) Base.A
  · next h => rfl

/-- Phase 5 preserves the amino acid sequence -/
theorem phase5_preserves_translation
    (codons : List Codon) (sel : CodonSelector)
    (hSyn : ∀ aa prev, translateCodon (sel aa prev) = aa) :
    (AvoidCpG codons sel).map translateCodon = codons.map translateCodon := by
  simp only [AvoidCpG, List.map_map]
  exact List.map_congr_left (fun c _ => cpg_subst_preserves c sel hSyn)

-- ────────────────────────────────────────────────────────────
-- End-to-End Pipeline Correctness
-- ────────────────────────────────────────────────────────────

/-- The full 5-phase pipeline (abstract) -/
def OptimizePipeline (codons : List Codon) (sel : CodonSelector)
    (forbidden : List Codon) : List Codon :=
  let phase1 := GTAwareOptimize codons sel
  let phase2 := RemoveRestrictionSites phase1 forbidden sel
  -- Phase 3: resolve cross-codon constraints (would need adjacency tracking)
  -- Phase 4: mutagenesis fallback
  let phase5 := AvoidCpG phase2 sel
  phase5

/-- The pipeline preserves the amino acid sequence -/
theorem pipeline_preserves_protein
    (codons : List Codon) (sel : CodonSelector) (forbidden : List Codon)
    (hSyn : ∀ aa prev, translateCodon (sel aa prev) = aa) :
    (OptimizePipeline codons sel forbidden).map translateCodon =
      codons.map translateCodon := by
  simp only [OptimizePipeline]
  rw [phase5_preserves_translation _ sel hSyn]
  rw [phase2_preserves_translation _ forbidden sel hSyn]
  exact phase1_preserves_translation codons sel hSyn

end BioCompiler
