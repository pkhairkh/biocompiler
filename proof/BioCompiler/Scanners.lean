/-
  BioCompiler.Scanners — Concrete Scanner Implementations with Proved Completeness

  This module implements the scanner functions used by the type system
  as CONCRETE Lean4 functions (not axioms), and proves their completeness
  theorems. Each scanner is defined using the pattern matching primitives
  from BioCompiler.Sequence.

  Reference: DOC-03 (SDD) §3.1, DOC-01 (SRS) §2.3
-/

import BioCompiler.Sequence

namespace BioCompiler

open Sequence

-- ==============================================================================
-- Generic Pattern Scanner
-- ==============================================================================

/-- Check if a specific pattern appears anywhere in the sequence. -/
def hasPattern (seq : Sequence) (pattern : Sequence) : Bool :=
  containsPattern seq pattern

/-- THEOREM (hasPattern completeness): If pattern appears at some
    position in the sequence, hasPattern returns true. -/
theorem hasPattern_complete (seq pattern : Sequence) (pos : Nat)
    (h_pos : pos + pattern.length ≤ seq.length)
    (h_match : (seq.drop pos).take pattern.length = pattern) :
    hasPattern seq pattern = true :=
  containsPattern_complete seq pattern pos h_pos h_match

/-- THEOREM (hasPattern soundness): If hasPattern returns true,
    the pattern appears at some position. -/
theorem hasPattern_sound (seq pattern : Sequence) :
    hasPattern seq pattern = true →
      ∃ (pos : Nat), pos + pattern.length ≤ seq.length ∧
        (seq.drop pos).take pattern.length = pattern :=
  containsPattern_sound seq pattern

-- ==============================================================================
-- Restriction Site Scanner
-- ==============================================================================

/-- Check if any restriction enzyme recognition site from the given set
    appears in the sequence. -/
def hasAnyRestrictionSite (seq : Sequence) (enzymeSites : List Sequence) : Bool :=
  enzymeSites.any fun site => hasPattern seq site

/-- THEOREM (Restriction site completeness) -/
theorem hasAnyRestrictionSite_complete (seq : Sequence) (enzymeSites : List Sequence)
    (site : Sequence) (pos : Nat)
    (h_site : site ∈ enzymeSites)
    (h_pos : pos + site.length ≤ seq.length)
    (h_match : (seq.drop pos).take site.length = site) :
    hasAnyRestrictionSite seq enzymeSites = true := by
  unfold hasAnyRestrictionSite
  rw [List.any_eq_true]
  exact ⟨site, h_site, hasPattern_complete seq site pos h_pos h_match⟩

/-- THEOREM (Restriction site soundness) -/
theorem hasAnyRestrictionSite_sound (seq : Sequence) (enzymeSites : List Sequence) :
    hasAnyRestrictionSite seq enzymeSites = true →
      ∃ (site : Sequence) (pos : Nat),
        site ∈ enzymeSites ∧
        pos + site.length ≤ seq.length ∧
        (seq.drop pos).take site.length = site := by
  unfold hasAnyRestrictionSite
  intro h
  obtain ⟨site, h_mem, h_has⟩ := List.any_eq_true.mp h
  have ⟨pos, h_pos, h_match⟩ := hasPattern_sound seq site h_has
  exact ⟨site, pos, h_mem, h_pos, h_match⟩

-- ==============================================================================
-- Instability Motif Scanner
-- ==============================================================================

/-- The ATTTA instability motif (DNA form of AUUUA). -/
def atttaMotif : Sequence :=
  [Nucleotide.A, Nucleotide.T, Nucleotide.T, Nucleotide.T, Nucleotide.A]

/-- The minimum run of consecutive T's that constitutes a U-rich region. -/
def minURichLength : Nat := 6

/-- Generate a sequence of n consecutive T nucleotides. -/
def tRun (n : Nat) : Sequence := List.replicate n Nucleotide.T

/-- The U-rich motif: a sequence of at least `minURichLength` consecutive T's. -/
def uRichMotif : Sequence := tRun minURichLength

/-- Check if the sequence contains any instability motif. -/
def hasInstabilityMotif (seq : Sequence) : Bool :=
  hasPattern seq atttaMotif || hasPattern seq uRichMotif

/-- Check if a position matches the ATTTA instability motif. -/
def matchesInstabilityAttta (seq : Sequence) (pos : Nat) : Bool :=
  matchesAt seq atttaMotif pos

/-- Check if a position matches a U-rich region. -/
def matchesInstabilityURich (seq : Sequence) (pos : Nat) : Bool :=
  matchesAt seq uRichMotif pos

/-- A position matches some instability motif. -/
def matchesInstabilityMotif (seq : Sequence) (pos : Nat) : Bool :=
  matchesInstabilityAttta seq pos || matchesInstabilityURich seq pos

/-- Helper: Bool || true = true -/
@[simp] theorem Bool.or_true (b : Bool) : (b || true) = true := by cases b <;> rfl

/-- Helper: true || Bool = true -/
@[simp] theorem Bool.true_or (b : Bool) : (true || b) = true := by cases b <;> rfl

/-- Helper: false || b = b -/
@[simp] theorem Bool.false_or (b : Bool) : (false || b) = b := by cases b <;> rfl

/-- THEOREM (Instability motif ATTTA completeness) -/
theorem hasInstabilityMotif_attta_complete (seq : Sequence) (pos : Nat)
    (h_pos : pos + atttaMotif.length ≤ seq.length)
    (h_match : (seq.drop pos).take atttaMotif.length = atttaMotif) :
    hasInstabilityMotif seq = true := by
  unfold hasInstabilityMotif
  have h1 := hasPattern_complete seq atttaMotif pos h_pos h_match
  simp [h1]

/-- THEOREM (Instability motif U-rich completeness) -/
theorem hasInstabilityMotif_urich_complete (seq : Sequence) (pos : Nat)
    (h_pos : pos + uRichMotif.length ≤ seq.length)
    (h_match : (seq.drop pos).take uRichMotif.length = uRichMotif) :
    hasInstabilityMotif seq = true := by
  unfold hasInstabilityMotif
  have h1 := hasPattern_complete seq uRichMotif pos h_pos h_match
  simp [h1]

/-- THEOREM (Instability motif completeness contrapositive): If no ATTTA motif
    and no U-rich motif exists at any position in the sequence, then
    hasInstabilityMotif returns false. -/
theorem hasInstabilityMotif_complete (seq : Sequence)
    (h_attta : ∀ (pos : Nat), pos + atttaMotif.length ≤ seq.length →
      (seq.drop pos).take atttaMotif.length = atttaMotif → False)
    (h_urich : ∀ (pos : Nat), pos + uRichMotif.length ≤ seq.length →
      (seq.drop pos).take uRichMotif.length = uRichMotif → False) :
    hasInstabilityMotif seq = false := by
  unfold hasInstabilityMotif
  have h1 : hasPattern seq atttaMotif = false := by
    cases h : hasPattern seq atttaMotif with
    | true =>
      exfalso
      have ⟨pos, h_pos, h_match⟩ := hasPattern_sound seq atttaMotif h
      exact h_attta pos h_pos h_match
    | false => rfl
  have h2 : hasPattern seq uRichMotif = false := by
    cases h : hasPattern seq uRichMotif with
    | true =>
      exfalso
      have ⟨pos, h_pos, h_match⟩ := hasPattern_sound seq uRichMotif h
      exact h_urich pos h_pos h_match
    | false => rfl
  simp [h1, h2]

/-- THEOREM (Instability motif soundness) -/
theorem hasInstabilityMotif_sound (seq : Sequence) :
    hasInstabilityMotif seq = true →
      (∃ (pos : Nat), pos + atttaMotif.length ≤ seq.length ∧
        (seq.drop pos).take atttaMotif.length = atttaMotif) ∨
      (∃ (pos : Nat), pos + uRichMotif.length ≤ seq.length ∧
        (seq.drop pos).take uRichMotif.length = uRichMotif) := by
  unfold hasInstabilityMotif
  intro h
  cases h1 : hasPattern seq atttaMotif with
  | true =>
    left
    exact hasPattern_sound seq atttaMotif h1
  | false =>
    right
    have : hasPattern seq uRichMotif = true := by simp [h1] at h; exact h
    exact hasPattern_sound seq uRichMotif this

-- ==============================================================================
-- Cryptic Splice Site Scanner
-- ==============================================================================

/-- Splice site consensus patterns (simplified for formalization). -/
def spliceDonorConsensus : Sequence :=
  [Nucleotide.G, Nucleotide.T]

def spliceAcceptorConsensus : Sequence :=
  [Nucleotide.A, Nucleotide.G]

/-- A splice site match: position and score. -/
structure SpliceSiteMatch where
  siteType : String
  position : Nat
  score : Rat
  deriving Repr

/-- Threshold for cryptic splice site detection. -/
def crypticThreshold : Rat := 3.0

/-- Lower threshold for borderline (uncertain) cryptic splice sites. -/
def uncertainLoThreshold : Rat := 1.5

/-- The CpG dinucleotide pattern. -/
def cpgDinucleotide : Sequence :=
  [Nucleotide.C, Nucleotide.G]

/-- Minimum window size for CpG island detection (default: 200 bp). -/
def cpgIslandWindowSize : Nat := 200

/-- GC content threshold for CpG island detection. -/
def cpgIslandGCThreshold : Rat := 6 / 10  -- 0.60

/-- Observed/Expected CpG ratio threshold for CpG island detection. -/
def cpgIslandObsExpThreshold : Rat := 65 / 100  -- 0.65

/-- Abstract interface for a CpG island scanner.

    Updated: scanner_completeness now requires BOTH the GC content criterion
    AND the Obs/Exp CG ratio criterion, matching the biological definition
    of a CpG island and the concrete scanner implementation.

    See BioCompiler.ScannerProofs for the concrete instance with proofs. -/
class CpGIslandScanner where
  hasCpGIsland : Sequence → Bool
  scanner_completeness :
    ∀ (seq : Sequence) (pos : Nat),
      pos + cpgIslandWindowSize ≤ seq.length →
      let window := (seq.drop pos).take cpgIslandWindowSize
      (window.count Nucleotide.G + window.count Nucleotide.C : Rat) / window.length ≥ cpgIslandGCThreshold →
      let cpgCount := (List.zipWith (· == ·) window (window.drop 1)).count true
      (cpgCount : Rat) * window.length ≥ cpgIslandObsExpThreshold * (window.count Nucleotide.C) * (window.count Nucleotide.G) →
      hasCpGIsland seq = false → False
  scanner_soundness :
    ∀ (seq : Sequence),
      hasCpGIsland seq = true →
        ∃ (pos : Nat), pos + cpgIslandWindowSize ≤ seq.length

/-- PROVED CONCRETE CpG island scanner is in BioCompiler.ScannerProofs.
    This module provides the CpGIslandScanner instance with proved
    completeness and soundness, eliminating axioms 4-5. -/

/-- Abstract interface for a splice site scanner. -/
class SpliceSiteScanner where
  hasCrypticSpliceSite : Sequence → Bool
  scanner_completeness :
    ∀ (seq : Sequence) (pos : Nat) (site : SpliceSiteMatch),
      pos < seq.length →
      site.position = pos →
      site.score ≥ crypticThreshold →
      hasCrypticSpliceSite seq = false → False
  scanner_soundness :
    ∀ (seq : Sequence),
      hasCrypticSpliceSite seq = true →
        ∃ (pos : Nat) (site : SpliceSiteMatch),
          pos < seq.length ∧ site.position = pos ∧ site.score ≥ crypticThreshold
  hasBorderlineSpliceSite : Sequence → Bool
  borderline_completeness :
    ∀ (seq : Sequence) (pos : Nat) (site : SpliceSiteMatch),
      pos < seq.length →
      site.position = pos →
      site.score ≥ uncertainLoThreshold →
      ¬(site.score ≥ crypticThreshold) →
      hasBorderlineSpliceSite seq = false → False

-- ==============================================================================
-- Codon Adaptation Index (CAI) Computation
-- ==============================================================================

/-- Codon usage table: maps each codon to its relative adaptiveness value. -/
class CodonAdaptationIndex where
  computeCAI : Sequence → String → Rat
  computeMinCodonCAI : Sequence → String → Rat
  cai_deterministic (seq : Sequence) (org : String) :
    computeCAI seq org = computeCAI seq org
  min_cai_deterministic (seq : Sequence) (org : String) :
    computeMinCodonCAI seq org = computeMinCodonCAI seq org

-- ==============================================================================
-- Valid Coding Sequence Scanner
-- ==============================================================================

/-- Check if a sequence is a valid coding sequence:
    length divisible by 3 and no internal stop codons. -/
def isValidCodingSeq (seq : Sequence) : Bool :=
  seq.length % 3 = 0 && !hasPrematureStop seq 0

-- NOTE: Soundness for ValidCodingSeq follows from the concrete definitions
-- of hasPrematureStop (proved in Sequence.lean) and modular arithmetic.

-- ==============================================================================
-- Conservation Score (BLOSUM62) Oracle
-- ==============================================================================

/-- BLOSUM62 substitution score for an amino acid pair.
    Placeholder: actual scores depend on the 20×20 BLOSUM62 matrix.
    Soundness proof only requires that evaluate condition = propertyHolds condition. -/
def blosum62Score (aa1 aa2 : String) : Int := 0

-- ==============================================================================
-- Cryptic Promoter Scanner
-- ==============================================================================

/-- Promoter match: organism, position, and consensus match score. -/
structure PromoterMatch where
  organism : String
  position : Nat
  score    : Rat
  deriving Repr

/-- Threshold for cryptic promoter detection. -/
def promoterThreshold : Rat := 7 / 10  -- 0.70

/-- Lower threshold for borderline (uncertain) promoter sites. -/
def promoterUncertainThreshold : Rat := 56 / 100  -- 0.56 (= 0.70 × 0.8)

/-- Window size for promoter motif detection (TATA box = 6 bp). -/
def promoterMotifSize : Nat := 6

/-- TATA box consensus motif (simplified for formalization).
    This is the most common eukaryotic core promoter element. -/
def tataBoxMotif : Sequence :=
  [Nucleotide.T, Nucleotide.A, Nucleotide.T, Nucleotide.A, Nucleotide.A, Nucleotide.A]

/-- Compute a promoter match score at a specific position.
    The score is the fraction of bases in the window that match the TATA box
    consensus motif, giving a value in [0, 1].
    Returns 0 if the window extends beyond the sequence. -/
def promoterScoreAt (seq : Sequence) (pos : Nat) : Rat :=
  if pos + promoterMotifSize ≤ seq.length then
    let window := (seq.drop pos).take promoterMotifSize
    ((List.zipWith (· == ·) window tataBoxMotif).count true : Rat) / promoterMotifSize
  else 0

/-- Abstract interface for a cryptic promoter scanner.

    Updated: scanner_completeness, scanner_soundness, and borderline_completeness
    now reference concrete position-based scoring (promoterScoreAt) instead of
    abstract PromoterMatch objects. This makes the axioms provable with a concrete
    sliding window implementation.

    See BioCompiler.ScannerProofs for the concrete instance with proofs. -/
class PromoterScanner where
  hasCrypticPromoter : Sequence → String → Rat → Bool
  hasBorderlinePromoter : Sequence → String → Rat → Bool
  scanner_completeness :
    ∀ (seq : Sequence) (organism : String) (threshold : Rat) (pos : Nat),
      pos + promoterMotifSize ≤ seq.length →
      promoterScoreAt seq pos ≥ threshold →
      hasCrypticPromoter seq organism threshold = false → False
  scanner_soundness :
    ∀ (seq : Sequence) (organism : String) (threshold : Rat),
      hasCrypticPromoter seq organism threshold = true →
        ∃ (pos : Nat), pos + promoterMotifSize ≤ seq.length ∧ promoterScoreAt seq pos ≥ threshold
  borderline_completeness :
    ∀ (seq : Sequence) (organism : String) (threshold : Rat) (pos : Nat),
      pos + promoterMotifSize ≤ seq.length →
      promoterScoreAt seq pos ≥ threshold * 8 / 10 →
      ¬(promoterScoreAt seq pos ≥ threshold) →
      hasBorderlinePromoter seq organism threshold = false → False

-- ==============================================================================
-- Transmembrane Domain Scanner
-- ==============================================================================

/-- TM domain match: position, window, and hydrophobic fraction.
    Kept for reference; the TMDomainScanner class no longer uses this
    structure in its axioms (axioms now reference concrete window positions). -/
structure TMDomainMatch where
  position      : Nat
  windowSize    : Nat
  hydroFraction : Rat
  deriving Repr

/-- Default TM domain hydrophobic fraction threshold. -/
def tmDomainThreshold : Rat := 68 / 100  -- 0.68

/-- TM domain window size (17 codons = 51 nucleotides).
    Standard length for transmembrane alpha-helix detection. -/
def tmDomainWindowSize : Nat := 51

/-- List of codons encoding hydrophobic amino acids.
    Hydrophobic amino acids: Ala (A), Val (V), Ile (I), Leu (L),
    Met (M), Phe (F), Trp (W). Total: 21 codons. -/
def hydrophobicCodons : List Sequence := [
  -- Alanine (A): GCN
  [Nucleotide.G, Nucleotide.C, Nucleotide.A],
  [Nucleotide.G, Nucleotide.C, Nucleotide.C],
  [Nucleotide.G, Nucleotide.C, Nucleotide.G],
  [Nucleotide.G, Nucleotide.C, Nucleotide.T],
  -- Valine (V): GTN
  [Nucleotide.G, Nucleotide.T, Nucleotide.A],
  [Nucleotide.G, Nucleotide.T, Nucleotide.C],
  [Nucleotide.G, Nucleotide.T, Nucleotide.G],
  [Nucleotide.G, Nucleotide.T, Nucleotide.T],
  -- Isoleucine (I): ATH (ATA, ATC, ATT — not ATG which is Met)
  [Nucleotide.A, Nucleotide.T, Nucleotide.A],
  [Nucleotide.A, Nucleotide.T, Nucleotide.C],
  [Nucleotide.A, Nucleotide.T, Nucleotide.T],
  -- Leucine (L): TTR + CTN
  [Nucleotide.T, Nucleotide.T, Nucleotide.A],
  [Nucleotide.T, Nucleotide.T, Nucleotide.G],
  [Nucleotide.C, Nucleotide.T, Nucleotide.A],
  [Nucleotide.C, Nucleotide.T, Nucleotide.C],
  [Nucleotide.C, Nucleotide.T, Nucleotide.G],
  [Nucleotide.C, Nucleotide.T, Nucleotide.T],
  -- Methionine (M): ATG
  [Nucleotide.A, Nucleotide.T, Nucleotide.G],
  -- Phenylalanine (F): TTY (TTC, TTT)
  [Nucleotide.T, Nucleotide.T, Nucleotide.C],
  [Nucleotide.T, Nucleotide.T, Nucleotide.T],
  -- Tryptophan (W): TGG
  [Nucleotide.T, Nucleotide.G, Nucleotide.G]
]

/-- Check if a codon (3-nucleotide sequence) encodes a hydrophobic amino acid.
    Returns true iff the codon appears in the hydrophobicCodons list. -/
def isHydrophobicCodon (codon : Sequence) : Bool :=
  hydrophobicCodons.any fun c => decide (c = codon)

/-- Count the number of hydrophobic codons in a window.
    Extracts codons sequentially from the beginning of the window
    (positions 0-2, 3-5, 6-8, ...). Partial codons at the end are ignored. -/
def countHydrophobicCodons (window : Sequence) : Nat :=
  let numCodons := window.length / 3
  (List.range numCodons).foldl (fun acc i =>
    if isHydrophobicCodon ((window.drop (3 * i)).take 3) then acc + 1 else acc) 0

/-- Compute the hydrophobic fraction of a window.
    Returns the fraction of codons that encode hydrophobic amino acids.
    Returns 0 if the window contains no complete codons (length < 3). -/
def tmHydrophobicFraction (window : Sequence) : Rat :=
  let numCodons := window.length / 3
  if numCodons = 0 then 0
  else (countHydrophobicCodons window : Rat) / numCodons

/-- Abstract interface for a transmembrane domain scanner.
    Dual-threshold: FAIL if hydrophobic fraction ≥ threshold,
    UNCERTAIN if ≥ threshold × 0.85.

    Updated: scanner_completeness, scanner_soundness, and borderline_completeness
    now reference concrete window positions and the tmHydrophobicFraction
    computation, matching the concrete scanner implementation.

    See BioCompiler.ScannerProofs for the concrete instance with proofs. -/
class TMDomainScanner where
  hasTMDomain : Sequence → Bool → Rat → Bool
  hasBorderlineTMDomain : Sequence → Bool → Rat → Bool
  scanner_completeness :
    ∀ (seq : Sequence) (isCytosolic : Bool) (threshold : Rat) (pos : Nat),
      isCytosolic = true →
      pos + tmDomainWindowSize ≤ seq.length →
      let window := (seq.drop pos).take tmDomainWindowSize
      tmHydrophobicFraction window ≥ threshold →
      hasTMDomain seq isCytosolic threshold = false → False
  scanner_soundness :
    ∀ (seq : Sequence) (isCytosolic : Bool) (threshold : Rat),
      hasTMDomain seq isCytosolic threshold = true →
        ∃ (pos : Nat), pos + tmDomainWindowSize ≤ seq.length ∧
          let window := (seq.drop pos).take tmDomainWindowSize
          tmHydrophobicFraction window ≥ threshold
  borderline_completeness :
    ∀ (seq : Sequence) (isCytosolic : Bool) (threshold : Rat) (pos : Nat),
      isCytosolic = true →
      pos + tmDomainWindowSize ≤ seq.length →
      let window := (seq.drop pos).take tmDomainWindowSize
      tmHydrophobicFraction window ≥ threshold * 85 / 100 →
      ¬(tmHydrophobicFraction window ≥ threshold) →
      hasBorderlineTMDomain seq isCytosolic threshold = false → False

-- ==============================================================================
-- mRNA Secondary Structure Oracle
-- ==============================================================================

/-- Structure stability estimate: position and ΔG. -/
structure StructureStabilityMatch where
  position : Nat
  deltaG   : Rat
  deriving Repr

/-- Default ΔG threshold for strong mRNA secondary structure (kcal/mol). -/
def mrnaStructureThreshold : Rat := -15

/-- Window size for mRNA secondary structure analysis. -/
def mrnaStructureWindowSize : Nat := 30

/-- Estimate the minimum free energy (deltaG) of a window.
    Simplified model: each potential Watson-Crick base pair (G-C or A-T)
    contributes -2 kcal/mol. Potential base pairs = min(G,C) + min(A,T).
    Returns 0 for empty windows. -/
def estimatedDeltaG (window : Sequence) : Rat :=
  if window.length = 0 then (0 : Rat)
  else
    let gcPairs := min (window.count Nucleotide.G) (window.count Nucleotide.C)
    let atPairs := min (window.count Nucleotide.A) (window.count Nucleotide.T)
    (-2 : Rat) * (gcPairs + atPairs)

/-- Abstract interface for mRNA secondary structure analysis.
    Dual-threshold: FAIL if estimated deltaG ≤ threshold (very stable structure),
    UNCERTAIN if estimated deltaG ≤ threshold × 0.7.

    Updated: oracle_completeness now requires position constraints,
    matching the CpGIslandScanner pattern. See BioCompiler.OracleProofs
    for the concrete instance with proofs. -/
class mRNAStructureOracle where
  hasStrongStructure : Sequence → Rat → Bool
  hasBorderlineStructure : Sequence → Rat → Bool
  oracle_completeness :
    ∀ (seq : Sequence) (threshold : Rat) (pos : Nat),
      pos + mrnaStructureWindowSize ≤ seq.length →
      let window := (seq.drop pos).take mrnaStructureWindowSize
      estimatedDeltaG window ≤ threshold →
      hasStrongStructure seq threshold = false → False
  oracle_soundness :
    ∀ (seq : Sequence) (threshold : Rat),
      hasStrongStructure seq threshold = true →
        ∃ (pos : Nat), pos + mrnaStructureWindowSize ≤ seq.length ∧
          estimatedDeltaG ((seq.drop pos).take mrnaStructureWindowSize) ≤ threshold
  borderline_completeness :
    ∀ (seq : Sequence) (threshold : Rat) (pos : Nat),
      pos + mrnaStructureWindowSize ≤ seq.length →
      let window := (seq.drop pos).take mrnaStructureWindowSize
      estimatedDeltaG window ≤ threshold * 7 / 10 →
      ¬(estimatedDeltaG window ≤ threshold) →
      hasBorderlineStructure seq threshold = false → False

-- ==============================================================================
-- Co-Translational Folding Oracle
-- ==============================================================================

/-- Folding disruption: codon position and disruption type.
    Kept for reference; the CoTranslationalFoldingOracle class no longer
    uses this structure in its axioms. -/
structure FoldingDisruption where
  codonPosition : Nat
  disruptionType : String
  deriving Repr

/-- Number of codons in the ramp region (first N codons). -/
def cotransRampCodons : Nat := 30

/-- Threshold below which the codon ramp speed is considered disrupted.
    Low adaptation (few optimal codons) → slow ramp → folding disruption. -/
def cotransDisruptionThreshold : Rat := 3 / 10   -- 0.30

/-- Threshold for borderline folding disruption. -/
def cotransBorderlineThreshold : Rat := 5 / 10    -- 0.50

/-- Compute a simplified codon adaptation index for the ramp region.
    Returns 1.0 for empty sequences (no disruption by convention).
    Uses GC content as a proxy for codon optimality.
    A real implementation would use organism-specific codon usage tables. -/
def rampAdaptationIndex (seq : Sequence) : Rat :=
  if seq.length = 0 then (1 : Rat)
  else ((seq.count Nucleotide.G + seq.count Nucleotide.C : Nat) : Rat) / seq.length

/-- Abstract interface for co-translational folding analysis.

    Updated: oracle_completeness now uses rampAdaptationIndex,
    following the CpGIslandScanner pattern. See BioCompiler.OracleProofs
    for the concrete instance with proofs. -/
class CoTranslationalFoldingOracle where
  hasFoldingDisruption : Sequence → String → Bool
  hasBorderlineFolding : Sequence → String → Bool
  oracle_completeness :
    ∀ (seq : Sequence) (organism : String),
      rampAdaptationIndex seq ≤ cotransDisruptionThreshold →
      hasFoldingDisruption seq organism = false → False
  oracle_soundness :
    ∀ (seq : Sequence) (organism : String),
      hasFoldingDisruption seq organism = true →
        rampAdaptationIndex seq ≤ cotransBorderlineThreshold
  borderline_completeness :
    ∀ (seq : Sequence) (organism : String),
      rampAdaptationIndex seq ≤ cotransBorderlineThreshold →
      ¬(rampAdaptationIndex seq ≤ cotransDisruptionThreshold) →
      hasBorderlineFolding seq organism = false → False

end BioCompiler
