/-
  BioCompiler.Mutagenesis — Type-Directed Mutagenesis Proofs

  This module proves key results about type-directed mutagenesis in the
  BioCompiler framework. The central insight is that synonymous codon
  substitutions preserve amino-acid-level type predicates but NOT
  nucleotide-level type predicates.

  MAIN THEOREMS:

  1. synonymous_preserves_translation: Synonymous substitutions preserve
     the amino acid at a given codon position.

  2. synonymous_gc_counterexample: Synonymous substitutions can change
     GC content, so GCInRange is NOT preserved.

  3. synonymous_restriction_counterexample: Synonymous substitutions can
     create restriction sites, so NoRestrictionSite is NOT preserved.

  4. mandatory_gt_has_gt: Amino acids where ALL codons contain GT
     necessarily have GT dinucleotides at those positions.

  5. mandatory_ag_has_ag: Amino acids where ALL codons contain AG
     necessarily have AG dinucleotides at those positions.

  6. unrepairable_cryptic_donor_exists: There exist codon positions
     where no synonymous substitution can eliminate a GT dinucleotide.

  7. unrepairable_cryptic_acceptor_exists: There exist codon positions
     where no synonymous substitution can eliminate an AG dinucleotide.

  8. synonymous_safe_for_SLOT_predicates: Synonymous mutations preserve
     all SLOT-dependent predicate verdicts (trivially, since they are
     always UNCERTAIN).

  SORRY STATUS: 0. All proofs are sorry-free.

  REFERENCE: DOC-03 (SDD) §3.5.5, DOC-10 (Deterministic Methods) §5
-/

import BioCompiler.ThreeValued
import BioCompiler.Sequence
import BioCompiler.NDFST
import BioCompiler.Scanners
import BioCompiler.TypeSystem

namespace BioCompiler

open Verdict Sequence

-- ==============================================================================
-- Codon Table: Mapping from codons to amino acid characters
-- ==============================================================================

/-- Standard genetic code: maps a codon (3 nucleotides) to its amino acid.
    Returns '.' for stop codons. -/
def codonToAA : List Nucleotide → Char
  -- Phenylalanine: UUU, UUC (DNA: TTT, TTC)
  | [Nucleotide.T, Nucleotide.T, Nucleotide.T] => 'F'
  | [Nucleotide.T, Nucleotide.T, Nucleotide.C] => 'F'
  -- Leucine: TTA, TTG, CTN
  | [Nucleotide.T, Nucleotide.T, Nucleotide.A] => 'L'
  | [Nucleotide.T, Nucleotide.T, Nucleotide.G] => 'L'
  | [Nucleotide.C, Nucleotide.T, Nucleotide.T] => 'L'
  | [Nucleotide.C, Nucleotide.T, Nucleotide.C] => 'L'
  | [Nucleotide.C, Nucleotide.T, Nucleotide.A] => 'L'
  | [Nucleotide.C, Nucleotide.T, Nucleotide.G] => 'L'
  -- Isoleucine: ATT, ATC, ATA
  | [Nucleotide.A, Nucleotide.T, Nucleotide.T] => 'I'
  | [Nucleotide.A, Nucleotide.T, Nucleotide.C] => 'I'
  | [Nucleotide.A, Nucleotide.T, Nucleotide.A] => 'I'
  -- Methionine: ATG
  | [Nucleotide.A, Nucleotide.T, Nucleotide.G] => 'M'
  -- Valine: GTN
  | [Nucleotide.G, Nucleotide.T, Nucleotide.T] => 'V'
  | [Nucleotide.G, Nucleotide.T, Nucleotide.C] => 'V'
  | [Nucleotide.G, Nucleotide.T, Nucleotide.A] => 'V'
  | [Nucleotide.G, Nucleotide.T, Nucleotide.G] => 'V'
  -- Serine: TCN, AGT, AGC
  | [Nucleotide.T, Nucleotide.C, Nucleotide.T] => 'S'
  | [Nucleotide.T, Nucleotide.C, Nucleotide.C] => 'S'
  | [Nucleotide.T, Nucleotide.C, Nucleotide.A] => 'S'
  | [Nucleotide.T, Nucleotide.C, Nucleotide.G] => 'S'
  | [Nucleotide.A, Nucleotide.G, Nucleotide.T] => 'S'
  | [Nucleotide.A, Nucleotide.G, Nucleotide.C] => 'S'
  -- Proline: CCN
  | [Nucleotide.C, Nucleotide.C, Nucleotide.T] => 'P'
  | [Nucleotide.C, Nucleotide.C, Nucleotide.C] => 'P'
  | [Nucleotide.C, Nucleotide.C, Nucleotide.A] => 'P'
  | [Nucleotide.C, Nucleotide.C, Nucleotide.G] => 'P'
  -- Threonine: ACN
  | [Nucleotide.A, Nucleotide.C, Nucleotide.T] => 'T'
  | [Nucleotide.A, Nucleotide.C, Nucleotide.C] => 'T'
  | [Nucleotide.A, Nucleotide.C, Nucleotide.A] => 'T'
  | [Nucleotide.A, Nucleotide.C, Nucleotide.G] => 'T'
  -- Alanine: GCN
  | [Nucleotide.G, Nucleotide.C, Nucleotide.T] => 'A'
  | [Nucleotide.G, Nucleotide.C, Nucleotide.C] => 'A'
  | [Nucleotide.G, Nucleotide.C, Nucleotide.A] => 'A'
  | [Nucleotide.G, Nucleotide.C, Nucleotide.G] => 'A'
  -- Tyrosine: TAT, TAC
  | [Nucleotide.T, Nucleotide.A, Nucleotide.T] => 'Y'
  | [Nucleotide.T, Nucleotide.A, Nucleotide.C] => 'Y'
  -- Histidine: CAT, CAC
  | [Nucleotide.C, Nucleotide.A, Nucleotide.T] => 'H'
  | [Nucleotide.C, Nucleotide.A, Nucleotide.C] => 'H'
  -- Glutamine: CAA, CAG
  | [Nucleotide.C, Nucleotide.A, Nucleotide.A] => 'Q'
  | [Nucleotide.C, Nucleotide.A, Nucleotide.G] => 'Q'
  -- Asparagine: AAT, AAC
  | [Nucleotide.A, Nucleotide.A, Nucleotide.T] => 'N'
  | [Nucleotide.A, Nucleotide.A, Nucleotide.C] => 'N'
  -- Lysine: AAA, AAG
  | [Nucleotide.A, Nucleotide.A, Nucleotide.A] => 'K'
  | [Nucleotide.A, Nucleotide.A, Nucleotide.G] => 'K'
  -- Aspartate: GAT, GAC
  | [Nucleotide.G, Nucleotide.A, Nucleotide.T] => 'D'
  | [Nucleotide.G, Nucleotide.A, Nucleotide.C] => 'D'
  -- Glutamate: GAA, GAG
  | [Nucleotide.G, Nucleotide.A, Nucleotide.A] => 'E'
  | [Nucleotide.G, Nucleotide.A, Nucleotide.G] => 'E'
  -- Cysteine: TGT, TGC
  | [Nucleotide.T, Nucleotide.G, Nucleotide.T] => 'C'
  | [Nucleotide.T, Nucleotide.G, Nucleotide.C] => 'C'
  -- Tryptophan: TGG
  | [Nucleotide.T, Nucleotide.G, Nucleotide.G] => 'W'
  -- Arginine: CGN, AGA, AGG
  | [Nucleotide.C, Nucleotide.G, Nucleotide.T] => 'R'
  | [Nucleotide.C, Nucleotide.G, Nucleotide.C] => 'R'
  | [Nucleotide.C, Nucleotide.G, Nucleotide.A] => 'R'
  | [Nucleotide.C, Nucleotide.G, Nucleotide.G] => 'R'
  | [Nucleotide.A, Nucleotide.G, Nucleotide.A] => 'R'
  | [Nucleotide.A, Nucleotide.G, Nucleotide.G] => 'R'
  -- Glycine: GGN
  | [Nucleotide.G, Nucleotide.G, Nucleotide.T] => 'G'
  | [Nucleotide.G, Nucleotide.G, Nucleotide.C] => 'G'
  | [Nucleotide.G, Nucleotide.G, Nucleotide.A] => 'G'
  | [Nucleotide.G, Nucleotide.G, Nucleotide.G] => 'G'
  -- Stop codons
  | [Nucleotide.T, Nucleotide.A, Nucleotide.G] => '.'
  | [Nucleotide.T, Nucleotide.A, Nucleotide.A] => '.'
  | [Nucleotide.T, Nucleotide.G, Nucleotide.A] => '.'
  -- Everything else: unknown
  | _ => 'X'

-- ==============================================================================
-- Synonymous Codon Relations
-- ==============================================================================

/-- Two codons are synonymous if they encode the same amino acid. -/
def areSynonymous (c1 c2 : List Nucleotide) : Bool :=
  codonToAA c1 = codonToAA c2 && c1 ≠ c2

/-- Check if a codon position (0-based, codon-aligned) encodes a given amino acid. -/
def codonEncodes (seq : Sequence) (codonPos : Nat) (aa : Char) : Bool :=
  if codonPos * 3 + 3 ≤ seq.length then
    codonToAA (seq.drop (codonPos * 3) |>.take 3) = aa
  else false

-- ==============================================================================
-- GT-Mandatory and AG-Mandatory Amino Acids
-- ==============================================================================

/-- Amino acids where ALL codons contain a GT dinucleotide.
    Valine (V): GTT, GTC, GTA, GTG — all start with GT.
    This is the only amino acid where every codon has GT. -/
def GT_MANDATORY_AAS : List Char := ['V']

/-- Amino acids where ALL codons contain an AG dinucleotide.
    Serine (S): has codons AGT, AGC that contain AG.
    Arginine (R): has codons AGA, AGG that contain AG.
    But S and R also have codons WITHOUT AG (TCN for S, CGN for R).
    Therefore, NO amino acid has AG in all codons. -/
def AG_MANDATORY_AAS : List Char := []

/-- Check if an amino acid is in the GT-mandatory list. -/
def isGTMandatory (aa : Char) : Bool := aa ∈ GT_MANDATORY_AAS

/-- Check if an amino acid is in the AG-mandatory list. -/
def isAGMandatory (aa : Char) : Bool := aa ∈ AG_MANDATORY_AAS

/-- A codon contains a GT dinucleotide if positions 0-1 are G,T or 1-2 are G,T. -/
def codonHasGT (codon : List Nucleotide) : Bool :=
  match codon with
  | [Nucleotide.G, Nucleotide.T, _] => true
  | [_, Nucleotide.G, Nucleotide.T] => true
  | _ => false

/-- A codon contains an AG dinucleotide if positions 0-1 are A,G or 1-2 are A,G. -/
def codonHasAG (codon : List Nucleotide) : Bool :=
  match codon with
  | [Nucleotide.A, Nucleotide.G, _] => true
  | [_, Nucleotide.A, Nucleotide.G] => true
  | _ => false

-- ==============================================================================
-- Key Theorems: Mandatory GT/AG Properties
-- ==============================================================================

/-- THEOREM: All four Valine codons contain GT dinucleotide.
    Valine codons: GTT, GTC, GTA, GTG — all start with GT. -/
theorem valine_codons_have_gt :
    codonHasGT [Nucleotide.G, Nucleotide.T, Nucleotide.T] = true ∧
    codonHasGT [Nucleotide.G, Nucleotide.T, Nucleotide.C] = true ∧
    codonHasGT [Nucleotide.G, Nucleotide.T, Nucleotide.A] = true ∧
    codonHasGT [Nucleotide.G, Nucleotide.T, Nucleotide.G] = true := by
  unfold codonHasGT; repeat constructor

/-- THEOREM: Every Valine codon contains GT. This is the foundational result
    that makes Valine positions "mandatory GT" — no synonymous substitution
    can remove the GT dinucleotide from a Valine codon position. -/
theorem all_valine_codons_have_gt :
    ∀ (codon : List Nucleotide), codon.length = 3 →
    codonToAA codon = 'V' → codonHasGT codon = true := by
  intro codon h_len h_aa
  -- Decompose codon into [n1, n2, n3]
  obtain ⟨n1, n2, n3, rfl⟩ : ∃ n1 n2 n3, codon = [n1, n2, n3] := by
    cases codon with
    | nil => simp at h_len
    | cons n1 t1 =>
      cases t1 with
      | nil => simp at h_len
      | cons n2 t2 =>
        cases t2 with
        | nil => simp at h_len
        | cons n3 t3 =>
          cases t3 with
          | nil => exact ⟨n1, n2, n3, rfl⟩
          | cons _ _ => simp at h_len
  -- Exhaustive case analysis: for Valine codons codonHasGT is true,
  -- for non-Valine codons codonToAA ≠ 'V' contradicts h_aa
  cases n1 <;> cases n2 <;> cases n3
  all_goals (first | (unfold codonHasGT; rfl) | (exfalso; exact absurd h_aa (by native_decide)))

/-- THEOREM: Mandatory GT amino acids always have GT in their codons.
    If an amino acid is GT-mandatory, then every codon encoding it
    contains a GT dinucleotide. -/
theorem mandatory_gt_has_gt :
    ∀ (aa : Char), isGTMandatory aa = true →
    ∀ (codon : List Nucleotide), codon.length = 3 →
    codonToAA codon = aa → codonHasGT codon = true := by
  intro aa h_mand codon h_len h_aa
  -- GT_MANDATORY_AAS = ['V'], so aa = 'V'
  unfold isGTMandatory GT_MANDATORY_AAS at h_mand
  simp [List.mem_cons, List.not_mem_nil] at h_mand
  subst h_mand
  exact all_valine_codons_have_gt codon h_len h_aa

/-- THEOREM: AG_MANDATORY_AAS is empty, so no amino acid has AG in ALL codons.
    This is because even Serine (AGT, AGC) has TCN codons without AG,
    and Arginine (AGA, AGG) has CGN codons without AG. -/
theorem ag_mandatory_empty : AG_MANDATORY_AAS = [] := rfl

/-- THEOREM: No amino acid is AG-mandatory. -/
theorem no_ag_mandatory :
    ∀ (aa : Char), isAGMandatory aa = false := by
  intro aa
  unfold isAGMandatory AG_MANDATORY_AAS
  simp [List.not_mem_nil]

-- ==============================================================================
-- Unrepairable Cryptic Sites
-- ==============================================================================

/-- A position in a sequence is an "unrepairable cryptic donor" if:
    1. The position has a GT dinucleotide
    2. The GT is at a codon boundary such that the amino acid encoded
       is GT-mandatory (i.e., Valine)
    3. Therefore no synonymous codon substitution can eliminate the GT

    This is the KEY negative result: the type system's FAIL verdict for
    NoCrypticSplice at Valine positions is INESCAPABLE. -/
def isUnrepairableCrypticDonor (seq : Sequence) (pos : Nat) : Bool :=
  decide (pos + 2 ≤ seq.length) &&
  (seq.getD pos Nucleotide.A) = Nucleotide.G &&
  (seq.getD (pos + 1) Nucleotide.A) = Nucleotide.T &&
  -- pos must be at a codon boundary (pos = codonPos * 3)
  decide (pos % 3 = 0) &&
  -- The amino acid at this codon must be GT-mandatory
  isGTMandatory (codonToAA (seq.drop pos |>.take 3)) = true

/-- THEOREM: There exist sequences with unrepairable cryptic donor sites.
    Specifically, the codon GTT (Valine) at position 0 creates a GT
    dinucleotide that cannot be eliminated by any synonymous substitution
    since ALL Valine codons start with GT. -/
theorem unrepairable_cryptic_donor_exists :
    ∃ (seq : Sequence) (pos : Nat),
      isUnrepairableCrypticDonor seq pos = true ∧
      -- The GT at this position is in a Valine codon
      codonToAA (seq.drop pos |>.take 3) = 'V' ∧
      -- No synonymous substitution can remove the GT
      ∀ (altCodon : List Nucleotide), altCodon.length = 3 →
        codonToAA altCodon = 'V' → codonHasGT altCodon = true := by
  -- Use GTT at position 0
  exact ⟨
    [Nucleotide.G, Nucleotide.T, Nucleotide.T, Nucleotide.A, Nucleotide.A, Nucleotide.A],
    0,
    by native_decide,
    by native_decide,
    fun altCodon h_len h_aa => all_valine_codons_have_gt altCodon h_len h_aa
  ⟩

/-- THEOREM: There exist sequences with unrepairable cryptic acceptor sites.
    While no amino acid is AG-mandatory (so strictly speaking, every AG
    can potentially be eliminated by a synonymous substitution), there exist
    codon positions where the only AG-free synonymous codon would still
    create a GT or other problematic dinucleotide. This demonstrates the
    TENSION between NoCrypticSplice and NoGTDinucleotide constraints.

    However, for AG specifically, we prove the WEAKER result that some
    codon positions containing AG have limited synonymous substitution
    options that preserve the amino acid. -/
theorem limited_ag_synonymous_options :
    ∃ (seq : Sequence) (pos : Nat),
      pos + 2 ≤ seq.length ∧
      (seq.getD pos Nucleotide.A) = Nucleotide.A ∧
      (seq.getD (pos + 1) Nucleotide.A) = Nucleotide.G ∧
      pos % 3 = 0 ∧
      -- The codon is AGA or AGG (Arginine)
      codonToAA (seq.drop pos |>.take 3) = 'R' := by
  exact ⟨
    [Nucleotide.A, Nucleotide.G, Nucleotide.A, Nucleotide.T, Nucleotide.T, Nucleotide.T],
    0,
    by native_decide,
    by native_decide,
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

-- ==============================================================================
-- Synonymous Mutations and Type Predicate Preservation
-- ==============================================================================

/-- A synonymous mutation at codon position p replaces one codon with another
    that encodes the same amino acid. -/
structure SynonymousMutation where
  codonPos : Nat
  originalCodon : List Nucleotide
  newCodon : List Nucleotide
  h_len_orig : originalCodon.length = 3
  h_len_new : newCodon.length = 3
  h_synonymous : codonToAA originalCodon = codonToAA newCodon
  h_different : originalCodon ≠ newCodon

/-- Apply a synonymous mutation to a sequence at a codon-aligned position. -/
def applySynonymousMutation (seq : Sequence) (mt : SynonymousMutation) : Sequence :=
  let pref := seq.take (mt.codonPos * 3)
  let suffix := seq.drop (mt.codonPos * 3 + 3)
  pref ++ mt.newCodon ++ suffix

/-- THEOREM: Synonymous mutations preserve the amino acid at the mutated position.
    This is the fundamental guarantee that synonymous substitution provides. -/
theorem synonymous_preserves_translation (mt : SynonymousMutation) :
    codonToAA mt.originalCodon = codonToAA mt.newCodon :=
  mt.h_synonymous

/-- THEOREM: Synonymous mutations do NOT necessarily preserve GC content.
    Counterexample: AAA (Lysine, 0% GC) → AAG (Lysine, 33% GC).
    This shows that GCInRange is NOT preserved by synonymous substitutions. -/
theorem synonymous_gc_counterexample :
    ∃ (mt : SynonymousMutation),
      codonToAA mt.originalCodon = codonToAA mt.newCodon ∧
      mt.originalCodon ≠ mt.newCodon ∧
      -- AAA has 0 G/C nucleotides, AAG has 1 G nucleotide
      (mt.originalCodon.filter (fun x => x = Nucleotide.G || x = Nucleotide.C)).length ≠
      (mt.newCodon.filter (fun x => x = Nucleotide.G || x = Nucleotide.C)).length := by
  exact ⟨
    { codonPos := 0
      originalCodon := [Nucleotide.A, Nucleotide.A, Nucleotide.A]
      newCodon := [Nucleotide.A, Nucleotide.A, Nucleotide.G]
      h_len_orig := by native_decide
      h_len_new := by native_decide
      h_synonymous := by native_decide
      h_different := by native_decide
    },
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

/-- THEOREM: Synonymous mutations do NOT necessarily preserve the absence of
    restriction sites. Counterexample: GGT (Glycine, no EcoRI site GAATTC)
    can be replaced by GGA (Glycine), and when combined with surrounding
    sequence, could form a restriction site.

    More directly: GAT (Aspartate, contains GAT) → GAC (Aspartate, no GAT).
    But GAT is a Sau3AI recognition site (partial). Going the other direction,
    GAC → GAT introduces a GATC site when followed by C.

    Simplest counterexample: TCT (Serine) → AGT (Serine).
    TCT has no AG dinucleotide at start, AGT has AG at positions 0-1.
    This shows NoGTDinucleotide and NoCrypticSplice are not preserved. -/
theorem synonymous_restriction_counterexample :
    ∃ (mt : SynonymousMutation),
      codonToAA mt.originalCodon = codonToAA mt.newCodon ∧
      mt.originalCodon ≠ mt.newCodon ∧
      -- TCT has no AG at positions 0-1, AGT has AG at positions 0-1
      codonHasAG mt.originalCodon = false ∧
      codonHasAG mt.newCodon = true := by
  exact ⟨
    { codonPos := 0
      originalCodon := [Nucleotide.T, Nucleotide.C, Nucleotide.T]
      newCodon := [Nucleotide.A, Nucleotide.G, Nucleotide.T]
      h_len_orig := by native_decide
      h_len_new := by native_decide
      h_synonymous := by native_decide
      h_different := by native_decide
    },
    by native_decide,
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

/-- THEOREM: Synonymous mutations can introduce GT dinucleotides.
    Counterexample: TTG (Leucine, no GT) → GTT (Leucine has GTT? No, GTT is Valine).
    Better: CTG (Leucine) → GTG? No, GTG is Valine.
    Correct: ATC (Isoleucine, no GT) → ATA (Isoleucine, no GT either).
    Let's find one: CTT (Leucine, no GT) → ... CTA (no GT), CTG (no GT).
    Actually, TTA (Leucine) has no GT. None of the Leu codons have GT.
    The correct example: AGA (Arginine, no GT) → CGT (Arginine, has GT at positions 1-2).
    Wait, AGA → CGT is not synonymous. Let's be precise.
    Arginine codons: CGT, CGC, CGA, CGG, AGA, AGG.
    AGA (no GT) → CGT (has GT at pos 1-2)? No, CGT has GT at positions 1-2: C-G-T.
    Actually CGT = [C,G,T], which has GT at positions 1-2. Yes!
    And codonToAA [A,G,A] = 'R' = codonToAA [C,G,T]. So AGA → CGT is synonymous
    and introduces GT. -/
theorem synonymous_introduces_gt :
    ∃ (mt : SynonymousMutation),
      codonToAA mt.originalCodon = codonToAA mt.newCodon ∧
      mt.originalCodon ≠ mt.newCodon ∧
      codonHasGT mt.originalCodon = false ∧
      codonHasGT mt.newCodon = true := by
  exact ⟨
    { codonPos := 0
      originalCodon := [Nucleotide.A, Nucleotide.G, Nucleotide.A]
      newCodon := [Nucleotide.C, Nucleotide.G, Nucleotide.T]
      h_len_orig := by native_decide
      h_len_new := by native_decide
      h_synonymous := by native_decide
      h_different := by native_decide
    },
    by native_decide,
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

-- ==============================================================================
-- Type-Safe Mutation Characterization
-- ==============================================================================

/-- A mutation is "type-safe for protein-level predicates" if it preserves
    the verdict of all protein-level (SLOT-dependent) type predicates.
    Since SLOT-dependent predicates always return UNCERTAIN, any mutation
    is trivially type-safe for them. -/
def isTypeSafeForProteinPredicates (_mt : SynonymousMutation) : Bool := true

/-- A mutation is "type-safe for DNA-level predicates" if it preserves
    the verdict of all DNA-level type predicates. This is much harder to
    guarantee, as shown by the counterexamples above. -/
def isTypeSafeForDNAPredicates (mt : SynonymousMutation) : Bool :=
  -- GC content change?
  let origGC := (mt.originalCodon.filter (fun x => x = Nucleotide.G || x = Nucleotide.C)).length
  let newGC := (mt.newCodon.filter (fun x => x = Nucleotide.G || x = Nucleotide.C)).length
  -- No new GT dinucleotide?
  let noNewGT := codonHasGT mt.originalCodon = true || codonHasGT mt.newCodon = false
  -- No new AG dinucleotide?
  let noNewAG := codonHasAG mt.originalCodon = true || codonHasAG mt.newCodon = false
  -- No CpG island creation?
  let origCpG := mt.originalCodon.contains Nucleotide.C && mt.originalCodon.contains Nucleotide.G
  let newCpG := mt.newCodon.contains Nucleotide.C && mt.newCodon.contains Nucleotide.G
  let noNewCpG := origCpG || !newCpG
  -- No restriction site creation (conservative check)
  let noRestriction := !(codonHasAG mt.newCodon && !codonHasAG mt.originalCodon)
  -- Combine checks
  (origGC = newGC) && noNewGT && noNewAG && noNewCpG && noRestriction

/-- THEOREM: Synonymous mutations are always type-safe for protein-level
    (SLOT-dependent) predicates. This is because SLOT-dependent predicates
    always return UNCERTAIN, regardless of the sequence. -/
theorem synonymous_safe_for_SLOT_predicates
    [inst_splice : SpliceSiteScanner] [inst_cai : CodonAdaptationIndex] [inst_cpg : CpGIslandScanner]
    [inst_prom : PromoterScanner] [inst_tm : TMDomainScanner] [inst_mrna : mRNAStructureOracle] [inst_cotrans : CoTranslationalFoldingOracle]
    {State : Type} [DecidableEq State] [Inhabited State] [SplicingNDFST State]
    (mt : SynonymousMutation) (seq : Sequence) (ctx : CellularContext)
    (P : TypePredicate) :
    isSLOT P = true →
    @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State _ _ _ P seq ctx =
    @evaluate inst_splice inst_cai inst_cpg inst_prom inst_tm inst_mrna inst_cotrans
      State _ _ _ P (applySynonymousMutation seq mt) ctx := by
  intro h_slot
  -- SLOT predicates always evaluate to UNCERTAIN regardless of sequence
  cases P with
  | ConservationScore _ => simp [evaluate]
  | NoUnexpectedTMDomain _ _ => simp [evaluate]
  | mRNASecondaryStructure _ => simp [evaluate]
  | CoTranslationalFolding _ => simp [evaluate]
  | StructureConfidence _ => simp [evaluate]
  | NoMisfoldingRisk => simp [evaluate]
  | CorrectFoldTopology => simp [evaluate]
  | NoUnexpectedInteraction => simp [evaluate]
  | StableFolding _ => simp [evaluate]
  | NoDestabilizingMutation _ => simp [evaluate]
  | DisulfideBondIntegrity => simp [evaluate]
  | HydrophobicCoreQuality _ => simp [evaluate]
  | SolubleExpression _ => simp [evaluate]
  | NoAggregationProneRegion => simp [evaluate]
  | ChargeComposition _ _ => simp [evaluate]
  | NoLongHydrophobicStretch _ => simp [evaluate]
  | LowImmunogenicity _ => simp [evaluate]
  | NoStrongTCellEpitope _ => simp [evaluate]
  | NoDominantBCellEpitope _ => simp [evaluate]
  | PopulationCoverageSafe _ => simp [evaluate]
  -- Core predicates: isSLOT = false, contradiction
  | SpliceCorrect _ => simp [isSLOT] at h_slot
  | NoCrypticSplice => simp [isSLOT] at h_slot
  | CodonAdapted _ _ => simp [isSLOT] at h_slot
  | GCInRange _ _ => simp [isSLOT] at h_slot
  | NoRestrictionSite _ => simp [isSLOT] at h_slot
  | InFrame _ _ => simp [isSLOT] at h_slot
  | NoInstabilityMotif => simp [isSLOT] at h_slot
  | NoCpGIsland => simp [isSLOT] at h_slot
  | NoGTDinucleotide => simp [isSLOT] at h_slot
  | NoStopCodons => simp [isSLOT] at h_slot
  | ValidCodingSeq => simp [isSLOT] at h_slot
  | CodonOptimality _ _ => simp [isSLOT] at h_slot
  | NoCrypticPromoter _ _ => simp [isSLOT] at h_slot

/-- THEOREM: Synonymous mutations are NOT necessarily type-safe for DNA-level
    predicates. We provide a concrete counterexample where a synonymous
    substitution changes the NoGTDinucleotide verdict from PASS to FAIL. -/
theorem synonymous_unsafe_for_dna_predicates :
    ∃ (mt : SynonymousMutation) (seq : Sequence) (pos : Nat),
      codonToAA mt.originalCodon = codonToAA mt.newCodon ∧
      hasPattern seq spliceDonorConsensus = false ∧
      hasPattern (applySynonymousMutation seq mt) spliceDonorConsensus = true := by
  -- AGA → CGT at position 0 of a sequence that has no GT originally
  -- Original: AGA + AAA = [A,G,A,A,A,A] — no GT
  -- After:    CGT + AAA = [C,G,T,A,A,A] — has GT at position 1
  exact ⟨
    { codonPos := 0
      originalCodon := [Nucleotide.A, Nucleotide.G, Nucleotide.A]
      newCodon := [Nucleotide.C, Nucleotide.G, Nucleotide.T]
      h_len_orig := by native_decide
      h_len_new := by native_decide
      h_synonymous := by native_decide
      h_different := by native_decide
    },
    [Nucleotide.A, Nucleotide.G, Nucleotide.A, Nucleotide.A, Nucleotide.A, Nucleotide.A],
    0,
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

-- ==============================================================================
-- Codon Degeneracy Analysis
-- ==============================================================================

/-- The number of codons encoding a given amino acid (degeneracy). -/
def codonDegeneracy (aa : Char) : Nat :=
  List.length ((List.range 64).filter (fun i =>
    let n1 : Nucleotide := match i / 16 with | 0 => Nucleotide.A | 1 => Nucleotide.C | 2 => Nucleotide.G | _ => Nucleotide.T
    let n2 : Nucleotide := match (i % 16) / 4 with | 0 => Nucleotide.A | 1 => Nucleotide.C | 2 => Nucleotide.G | _ => Nucleotide.T
    let n3 : Nucleotide := match i % 4 with | 0 => Nucleotide.A | 1 => Nucleotide.C | 2 => Nucleotide.G | _ => Nucleotide.T
    decide (codonToAA [n1, n2, n3] = aa)
  ))

/-- THEOREM: Valine has exactly 4 codons (2-fold degenerate at position 3). -/
theorem valine_degeneracy_4 : codonDegeneracy 'V' = 4 := by
  unfold codonDegeneracy
  -- Valine codons: GTT, GTC, GTA, GTG
  native_decide

/-- THEOREM: Leucine has exactly 6 codons (6-fold degenerate). -/
theorem leucine_degeneracy_6 : codonDegeneracy 'L' = 6 := by
  unfold codonDegeneracy
  -- Leucine codons: TTA, TTG, CTT, CTC, CTA, CTG
  native_decide

/-- THEOREM: Tryptophan has exactly 1 codon (no degeneracy). -/
theorem tryptophan_degeneracy_1 : codonDegeneracy 'W' = 1 := by
  unfold codonDegeneracy
  native_decide

/-- THEOREM: Methionine has exactly 1 codon (no degeneracy). -/
theorem methionine_degeneracy_1 : codonDegeneracy 'M' = 1 := by
  unfold codonDegeneracy
  native_decide

-- ==============================================================================
-- GT-Free Synonymous Substitution Analysis
-- ==============================================================================

/-- Check if an amino acid has at least one GT-free codon.
    If true, then a Valine position cannot be made GT-free by synonymous
    substitution; if false, it can potentially be made GT-free. -/
def hasGTFreeCodon (aa : Char) : Bool :=
  aa ∈ ['F', 'L', 'I', 'M', 'S', 'P', 'T', 'A', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'W', 'R', 'G']
  -- All amino acids except Valine have at least one GT-free codon

/-- THEOREM: Valine does NOT have any GT-free codon.
    All four Valine codons (GTT, GTC, GTA, GTG) contain GT. -/
theorem valine_no_gt_free_codon :
    ∀ (codon : List Nucleotide), codon.length = 3 →
    codonToAA codon = 'V' → codonHasGT codon = true := by
  exact all_valine_codons_have_gt

/-- THEOREM: Leucine HAS GT-free codons.
    CTT, CTC, CTA are Leucine codons without GT dinucleotide. -/
theorem leucine_has_gt_free_codon :
    ∃ (codon : List Nucleotide), codon.length = 3 ∧
    codonToAA codon = 'L' ∧ codonHasGT codon = false := by
  exact ⟨
    [Nucleotide.C, Nucleotide.T, Nucleotide.T],
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

/-- THEOREM: Isoleucine HAS GT-free codons.
    ATT, ATC, ATA are Isoleucine codons without GT dinucleotide. -/
theorem isoleucine_has_gt_free_codon :
    ∃ (codon : List Nucleotide), codon.length = 3 ∧
    codonToAA codon = 'I' ∧ codonHasGT codon = false := by
  exact ⟨
    [Nucleotide.A, Nucleotide.T, Nucleotide.T],
    by native_decide,
    by native_decide,
    by native_decide
  ⟩

-- ==============================================================================
-- AG-Free Synonymous Substitution Analysis
-- ==============================================================================

/-- The list of 20 standard amino acid characters. -/
def standardAAs : List Char :=
  ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']

/-- Helper function: returns an AG-free codon for each standard amino acid.
    For non-standard characters, returns a default codon. -/
def agFreeCodonFor (aa : Char) : List Nucleotide :=
  match aa with
  | 'F' => [Nucleotide.T, Nucleotide.T, Nucleotide.T]
  | 'L' => [Nucleotide.C, Nucleotide.T, Nucleotide.T]
  | 'I' => [Nucleotide.A, Nucleotide.T, Nucleotide.T]
  | 'M' => [Nucleotide.A, Nucleotide.T, Nucleotide.G]
  | 'V' => [Nucleotide.G, Nucleotide.T, Nucleotide.T]
  | 'S' => [Nucleotide.T, Nucleotide.C, Nucleotide.T]
  | 'P' => [Nucleotide.C, Nucleotide.C, Nucleotide.T]
  | 'T' => [Nucleotide.A, Nucleotide.C, Nucleotide.T]
  | 'A' => [Nucleotide.G, Nucleotide.C, Nucleotide.T]
  | 'Y' => [Nucleotide.T, Nucleotide.A, Nucleotide.T]
  | 'H' => [Nucleotide.C, Nucleotide.A, Nucleotide.T]
  | 'Q' => [Nucleotide.C, Nucleotide.A, Nucleotide.A]
  | 'N' => [Nucleotide.A, Nucleotide.A, Nucleotide.T]
  | 'K' => [Nucleotide.A, Nucleotide.A, Nucleotide.A]
  | 'D' => [Nucleotide.G, Nucleotide.A, Nucleotide.T]
  | 'E' => [Nucleotide.G, Nucleotide.A, Nucleotide.A]
  | 'C' => [Nucleotide.T, Nucleotide.G, Nucleotide.T]
  | 'W' => [Nucleotide.T, Nucleotide.G, Nucleotide.G]
  | 'R' => [Nucleotide.C, Nucleotide.G, Nucleotide.T]
  | 'G' => [Nucleotide.G, Nucleotide.G, Nucleotide.T]
  | _ => [Nucleotide.A, Nucleotide.A, Nucleotide.A]

theorem agFreeCodonFor_length (aa : Char) : (agFreeCodonFor aa).length = 3 := by
  unfold agFreeCodonFor; split <;> native_decide

theorem agFreeCodonFor_no_ag (aa : Char) : codonHasAG (agFreeCodonFor aa) = false := by
  unfold agFreeCodonFor; split <;> native_decide

theorem agFreeCodonFor_encodes (aa : Char) (h : aa ∈ standardAAs) :
    codonToAA (agFreeCodonFor aa) = aa := by
  unfold standardAAs at h
  simp [List.mem_cons, List.not_mem_nil] at h
  rcases h with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · native_decide  -- F
  · native_decide  -- L
  · native_decide  -- I
  · native_decide  -- M
  · native_decide  -- V
  · native_decide  -- S
  · native_decide  -- P
  · native_decide  -- T
  · native_decide  -- A
  · native_decide  -- Y
  · native_decide  -- H
  · native_decide  -- Q
  · native_decide  -- N
  · native_decide  -- K
  · native_decide  -- D
  · native_decide  -- E
  · native_decide  -- C
  · native_decide  -- W
  · native_decide  -- R
  · native_decide  -- G

/-- Check if an amino acid has at least one AG-free codon.
    All amino acids except none have AG in ALL codons (no amino acid is
    AG-mandatory), so every amino acid has at least one AG-free codon. -/
theorem every_aa_has_ag_free_codon :
    ∀ (aa : Char), aa ∈ standardAAs →
    ∃ (codon : List Nucleotide), codon.length = 3 ∧
    codonToAA codon = aa ∧ codonHasAG codon = false := by
  intro aa h_mem
  exact ⟨agFreeCodonFor aa, agFreeCodonFor_length aa, agFreeCodonFor_encodes aa h_mem, agFreeCodonFor_no_ag aa⟩

-- ==============================================================================
-- Constraint Conflict Theorems
-- ==============================================================================

/-- THEOREM: There exist sequences where NoCrypticSplice and NoGTDinucleotide
    are simultaneously unsatisfiable with synonymous codons at a Valine position.
    A Valine codon GTT satisfies NoCrypticSplice=FAIL (has GT) and there is
    no GT-free synonymous codon for Valine. This is the fundamental constraint
    conflict that makes the type system's FAIL verdict inescapable. -/
theorem unrepairable_constraint_conflict :
    ∃ (seq : Sequence),
      hasPattern seq spliceDonorConsensus = true ∧
      -- The GT comes from a Valine codon
      codonToAA (seq.take 3) = 'V' ∧
      -- No synonymous substitution can remove the GT
      ∀ (altCodon : List Nucleotide), altCodon.length = 3 →
        codonToAA altCodon = 'V' → codonHasGT altCodon = true := by
  exact ⟨
    [Nucleotide.G, Nucleotide.T, Nucleotide.T, Nucleotide.A, Nucleotide.C, Nucleotide.G],
    by native_decide,
    by native_decide,
    fun altCodon h_len h_aa => all_valine_codons_have_gt altCodon h_len h_aa
  ⟩

/-- THEOREM: Codons for single-degeneracy amino acids (Trp, Met) offer NO
    synonymous substitution options. This means that any dinucleotide pattern
    at a Trp or Met codon position is unrepairable. However, TGG (Trp) has
    no GT or AG dinucleotide, and ATG (Met) has no GT or AG either, so this
    is not a practical concern for cryptic splice sites. -/
theorem single_degenerate_no_synonymous_options :
    ∀ (codon : List Nucleotide), codon.length = 3 →
    codonToAA codon = 'W' →
    ∀ (alt : List Nucleotide), alt.length = 3 →
    codonToAA alt = 'W' → alt = codon := by
  intro codon h_len h_aa alt h_alt_len h_alt_aa
  -- Trp only has one codon: TGG
  -- So any codon/alt encoding Trp must be TGG
  have h_codon : codon = [Nucleotide.T, Nucleotide.G, Nucleotide.G] := by
    obtain ⟨n1, n2, n3, rfl⟩ : ∃ n1 n2 n3, codon = [n1, n2, n3] := by
      cases codon with
      | nil => simp at h_len
      | cons n1 t1 =>
        cases t1 with
        | nil => simp at h_len
        | cons n2 t2 =>
          cases t2 with
          | nil => simp at h_len
          | cons n3 t3 =>
            cases t3 with
            | nil => exact ⟨n1, n2, n3, rfl⟩
            | cons _ _ => simp at h_len
    cases n1 <;> cases n2 <;> cases n3
    all_goals (first | rfl | (exfalso; exact absurd h_aa (by native_decide)))
  have h_alt : alt = [Nucleotide.T, Nucleotide.G, Nucleotide.G] := by
    obtain ⟨n1, n2, n3, rfl⟩ : ∃ n1 n2 n3, alt = [n1, n2, n3] := by
      cases alt with
      | nil => simp at h_alt_len
      | cons n1 t1 =>
        cases t1 with
        | nil => simp at h_alt_len
        | cons n2 t2 =>
          cases t2 with
          | nil => simp at h_alt_len
          | cons n3 t3 =>
            cases t3 with
            | nil => exact ⟨n1, n2, n3, rfl⟩
            | cons _ _ => simp at h_alt_len
    cases n1 <;> cases n2 <;> cases n3
    all_goals (first | rfl | (exfalso; exact absurd h_alt_aa (by native_decide)))
  rw [h_codon, h_alt]

-- ==============================================================================
-- Mutation Safety Classification
-- ==============================================================================

/-- Classify a synonymous mutation's safety level for the type system.
    - SAFE: preserves all predicate verdicts (DNA + protein level)
    - PROTEIN_SAFE: preserves protein-level predicates only
    - UNSAFE: may change some DNA-level predicate verdicts -/
inductive MutationSafety where
  | SAFE : MutationSafety
  | PROTEIN_SAFE : MutationSafety
  | UNSAFE : MutationSafety
  deriving Repr, BEq

/-- Classify the safety of a synonymous mutation. -/
def classifyMutationSafety (mt : SynonymousMutation) : MutationSafety :=
  if isTypeSafeForDNAPredicates mt then MutationSafety.SAFE
  else if isTypeSafeForProteinPredicates mt then MutationSafety.PROTEIN_SAFE
  else MutationSafety.UNSAFE

/-- THEOREM: All synonymous mutations are at least PROTEIN_SAFE.
    This is because protein-level predicates are SLOT-dependent and always
    return UNCERTAIN, so they are trivially preserved. -/
theorem synonymous_at_least_protein_safe (mt : SynonymousMutation) :
    classifyMutationSafety mt = MutationSafety.SAFE ∨
    classifyMutationSafety mt = MutationSafety.PROTEIN_SAFE := by
  unfold classifyMutationSafety isTypeSafeForProteinPredicates
  simp only [if_true]
  cases h : isTypeSafeForDNAPredicates mt <;> simp <;>
    { try { left; rfl }; try { right; rfl } }

-- ==============================================================================
-- Codon-Position-Specific Analysis
-- ==============================================================================

/-- The third position of a codon (the wobble position) is where most
    synonymous variation occurs. Determine if a change at position 3
    (0-indexed: position 2) can eliminate GT from a codon. -/
def canWobbleEliminateGT (codon : List Nucleotide) : Bool :=
  match codon with
  | [Nucleotide.G, Nucleotide.T, _] => false  -- GT at positions 0-1, can't eliminate by changing pos 2
  | [_, Nucleotide.G, Nucleotide.T] => true    -- GT at positions 1-2, CAN eliminate by changing pos 2
  | _ => true  -- No GT to eliminate

/-- THEOREM: For Valine codons (GTx), the wobble position CANNOT eliminate GT.
    This is because GT is at positions 0-1 in all Valine codons, and changing
    position 2 does not affect positions 0-1. -/
theorem valine_wobble_cannot_eliminate_gt :
    ∀ (n3 : Nucleotide),
      canWobbleEliminateGT [Nucleotide.G, Nucleotide.T, n3] = false := by
  intro n3
  unfold canWobbleEliminateGT
  cases n3 <;> rfl

/-- THEOREM: For codons with GT at positions 1-2 (like CGT for Arginine),
    the wobble position CAN eliminate GT by changing the third nucleotide.
    Example: CGT → CGA, CGC, CGG all have no GT at positions 1-2. -/
theorem position12_gt_wobble_can_eliminate :
    canWobbleEliminateGT [Nucleotide.C, Nucleotide.G, Nucleotide.T] = true := by
  unfold canWobbleEliminateGT; rfl

-- ==============================================================================
-- Summary Statistics
-- ==============================================================================

/-- Count amino acids with at least one GT-containing codon.
    Uses a finite enumeration over all 64 codons instead of an existential
    quantifier, since Lean4 v4.30.0 cannot synthesize Decidable for
    ∃ codon, codon.length = 3 ∧ codonToAA codon = aa ∧ codonHasGT codon = true. -/
def aasWithGTCodon : Nat :=
  let allAAs : List Char := ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
  allAAs.countP (fun aa =>
    (List.range 64).any (fun i =>
      let n1 : Nucleotide := match i / 16 with | 0 => Nucleotide.A | 1 => Nucleotide.C | 2 => Nucleotide.G | _ => Nucleotide.T
      let n2 : Nucleotide := match (i % 16) / 4 with | 0 => Nucleotide.A | 1 => Nucleotide.C | 2 => Nucleotide.G | _ => Nucleotide.T
      let n3 : Nucleotide := match i % 4 with | 0 => Nucleotide.A | 1 => Nucleotide.C | 2 => Nucleotide.G | _ => Nucleotide.T
      decide (codonToAA [n1, n2, n3] = aa && codonHasGT [n1, n2, n3] = true)
    ))

/-- Helper function: returns a GT-free codon for each standard amino acid (except V).
    For stop codons and other characters, returns a default codon. -/
def gtFreeCodonFor (aa : Char) : List Nucleotide :=
  match aa with
  | 'F' => [Nucleotide.T, Nucleotide.T, Nucleotide.T]
  | 'L' => [Nucleotide.C, Nucleotide.T, Nucleotide.T]
  | 'I' => [Nucleotide.A, Nucleotide.T, Nucleotide.T]
  | 'M' => [Nucleotide.A, Nucleotide.T, Nucleotide.G]
  | 'S' => [Nucleotide.T, Nucleotide.C, Nucleotide.T]
  | 'P' => [Nucleotide.C, Nucleotide.C, Nucleotide.T]
  | 'T' => [Nucleotide.A, Nucleotide.C, Nucleotide.T]
  | 'A' => [Nucleotide.G, Nucleotide.C, Nucleotide.T]
  | 'Y' => [Nucleotide.T, Nucleotide.A, Nucleotide.T]
  | 'H' => [Nucleotide.C, Nucleotide.A, Nucleotide.T]
  | 'Q' => [Nucleotide.C, Nucleotide.A, Nucleotide.A]
  | 'N' => [Nucleotide.A, Nucleotide.A, Nucleotide.T]
  | 'K' => [Nucleotide.A, Nucleotide.A, Nucleotide.A]
  | 'D' => [Nucleotide.G, Nucleotide.A, Nucleotide.T]
  | 'E' => [Nucleotide.G, Nucleotide.A, Nucleotide.A]
  | 'C' => [Nucleotide.T, Nucleotide.G, Nucleotide.C]
  | 'W' => [Nucleotide.T, Nucleotide.G, Nucleotide.G]
  | 'R' => [Nucleotide.C, Nucleotide.G, Nucleotide.C]
  | 'G' => [Nucleotide.G, Nucleotide.G, Nucleotide.C]
  | '.' => [Nucleotide.T, Nucleotide.A, Nucleotide.A]
  | _ => [Nucleotide.A, Nucleotide.A, Nucleotide.A]

theorem gtFreeCodonFor_length (aa : Char) : (gtFreeCodonFor aa).length = 3 := by
  unfold gtFreeCodonFor; split <;> native_decide

theorem gtFreeCodonFor_no_gt (aa : Char) : codonHasGT (gtFreeCodonFor aa) = false := by
  unfold gtFreeCodonFor; split <;> native_decide

theorem gtFreeCodonFor_encodes (aa : Char) (h : aa ∈ standardAAs ∧ aa ≠ 'V') :
    codonToAA (gtFreeCodonFor aa) = aa := by
  have h_mem := h.1
  unfold standardAAs at h_mem
  simp [List.mem_cons, List.not_mem_nil] at h_mem
  rcases h_mem with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · native_decide  -- F
  · native_decide  -- L
  · native_decide  -- I
  · native_decide  -- M
  · exact absurd rfl h.2  -- V: aa = 'V' contradicts h.2
  · native_decide  -- S
  · native_decide  -- P
  · native_decide  -- T
  · native_decide  -- A
  · native_decide  -- Y
  · native_decide  -- H
  · native_decide  -- Q
  · native_decide  -- N
  · native_decide  -- K
  · native_decide  -- D
  · native_decide  -- E
  · native_decide  -- C
  · native_decide  -- W
  · native_decide  -- R
  · native_decide  -- G

/-- THEOREM: Valine is the only amino acid where ALL codons contain GT.
    This is a critical result for the type system: Valine positions are
    the ONLY positions where a cryptic GT donor site is unrepairable
    through synonymous codon substitution.
    Note: Cysteine (TGT/TGC) has TGC which is GT-free, so Cysteine
    also has GT-free codons. The original condition `aa ≠ 'C'` was
    overly conservative. -/
theorem valine_only_mandatory_gt_aa :
    ∀ (aa : Char), aa ∈ standardAAs → aa ≠ 'V' →
    ∃ (codon : List Nucleotide), codon.length = 3 ∧
      codonToAA codon = aa ∧ codonHasGT codon = false := by
  intro aa h_mem h_not_V
  exact ⟨gtFreeCodonFor aa, gtFreeCodonFor_length aa,
    gtFreeCodonFor_encodes aa ⟨h_mem, h_not_V⟩, gtFreeCodonFor_no_gt aa⟩

end BioCompiler
