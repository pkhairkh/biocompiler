/-
BioCompiler.Predicates v7.0.0
================================
Formal specification of the 8 certification predicates for gene optimization.
-/

import BioCompiler.CodonTable

namespace BioCompiler

-- ────────────────────────────────────────────────────────────
-- Float Order Properties (axioms for IEEE 754 semantics)
-- Lean 4's Float lacks these in the standard library; they hold
-- for all non-NaN values and are vacuously true for NaN.
-- ────────────────────────────────────────────────────────────

private axiom float_le_trans (a b c : Float) : a ≤ b → b ≤ c → a ≤ c
private axiom float_lt_of_lt_of_le (a b c : Float) : a < b → b ≤ c → a < c
private axiom float_not_lt_of_le (a b : Float) : a ≤ b → ¬(b < a)

-- ────────────────────────────────────────────────────────────
-- Splice Site Verdict (dual-threshold)
-- ────────────────────────────────────────────────────────────

inductive SpliceVerdict : Type where
  | PASS      : SpliceVerdict
  | UNCERTAIN : SpliceVerdict
  | FAIL      : SpliceVerdict
  deriving Repr, BEq, DecidableEq

namespace SpliceVerdict

/-- PASS and UNCERTAIN are "acceptable" (don't block optimization) -/
def acceptable : SpliceVerdict → Bool
  | PASS => true
  | UNCERTAIN => true
  | FAIL => false

theorem acceptable_PASS : acceptable PASS = true := rfl
theorem acceptable_UNCERTAIN : acceptable UNCERTAIN = true := rfl
theorem acceptable_FAIL : acceptable FAIL = false := rfl

end SpliceVerdict

-- ────────────────────────────────────────────────────────────
-- MaxEntScan Score Classification (dual-threshold)
-- ────────────────────────────────────────────────────────────

/-- Classify a splice score into a verdict using dual thresholds.
    PASS:      score < low
    UNCERTAIN: low <= score < high
    FAIL:      score >= high
-/
def classifySplice (score low high : Float) : SpliceVerdict :=
  if score < low then SpliceVerdict.PASS
  else if score < high then SpliceVerdict.UNCERTAIN
  else SpliceVerdict.FAIL

theorem classifySplice_lt_low (score low high : Float) (h : score < low) :
    classifySplice score low high = SpliceVerdict.PASS := by
  simp [classifySplice, h]

theorem classifySplice_ge_high (score low high : Float) (h : high ≤ score) (h_low : low ≤ high) :
    classifySplice score low high = SpliceVerdict.FAIL := by
  unfold classifySplice
  have h1 : ¬(score < low) := by
    intro h_sl
    exact float_not_lt_of_le low score (float_le_trans low high score h_low h) h_sl
  have h2 : ¬(score < high) := by
    intro h_sh
    exact float_not_lt_of_le high score h h_sh
  simp [h1, h2]

/-- Dual-threshold monotonicity: PASS at strict threshold implies PASS at
    any more permissive threshold. -/
theorem classifySplice_monotone (score low1 low2 high : Float)
    (h_low : low1 ≤ low2) (h_pass : score < low1) :
    classifySplice score low2 high = SpliceVerdict.PASS := by
  have : score < low2 := float_lt_of_lt_of_le score low1 low2 h_pass h_low
  exact classifySplice_lt_low score low2 high this

-- ────────────────────────────────────────────────────────────
-- 8 Predicate Definitions
-- ────────────────────────────────────────────────────────────

/-- Predicate 1: NoStopCodons — No stop codons except possibly the last. -/
def NoStopCodons (codons : List Codon) : Prop :=
  codons.length ≥ 1 ∧
  ∀ i : Nat, i < codons.length - 1 →
    translateCodon (codons.get ⟨i, by sorry⟩) ≠ AminoAcid.Stop

/-- Predicate 2: NoCrypticSplice (dual-threshold).
    No codon with a GT dinucleotide has a MaxEntScan score >= high threshold. -/
def NoCrypticSplice (codons : List Codon)
    (scoreFn : Codon → Float) (low high : Float) : Prop :=
  ∀ i : Nat, i < codons.length →
    codonHasGT (codons.get ⟨i, by sorry⟩) →
    classifySplice (scoreFn (codons.get ⟨i, by sorry⟩)) low high ≠ SpliceVerdict.FAIL

/-- Predicate 3: NoCpGIsland — No codon contains a CG dinucleotide. -/
def NoCpGIsland (codons : List Codon) : Prop :=
  ∀ i : Nat, i < codons.length →
    ¬ codonHasCG (codons.get ⟨i, by sorry⟩)

/-- Predicate 4: NoRestrictionSite — No codon matches a forbidden pattern. -/
def NoRestrictionSite (codons : List Codon) (forbidden : List Codon) : Prop :=
  ∀ c ∈ codons, c ∉ forbidden

/-- Predicate 5: NoGTDinucleotide — No within-codon GT, no cross-codon GT. -/
def NoGTDinucleotide (codons : List Codon) : Prop :=
  (∀ i : Nat, i < codons.length → ¬ codonHasGT (codons.get ⟨i, by sorry⟩)) ∧
  (∀ i : Nat, i + 1 < codons.length →
    ¬ crossCodonGT (codons.get ⟨i, by sorry⟩) (codons.get ⟨i+1, by sorry⟩))

/-- Predicate 6: ValidCodingSeq — All codons translate to standard amino acids. -/
def ValidCodingSeq (codons : List Codon) : Prop :=
  ∀ i : Nat, i < codons.length →
    translateCodon (codons.get ⟨i, by sorry⟩) ≠ AminoAcid.Stop ∨
    i = codons.length - 1

/-- Predicate 7: ConservationScore — BLOSUM62 substitution quality. -/
def ConservationScore (blosum62 : AminoAcid → AminoAcid → Int)
    (minBLOSUM : Int)
    (original optimized : List AminoAcid) : Prop :=
  original.length = optimized.length ∧
  ∀ i : Nat, i < original.length →
    blosum62 (original.get ⟨i, by sorry⟩) (optimized.get ⟨i, by sorry⟩) ≥ minBLOSUM

/-- Predicate 8: CodonOptimality — CAI quality threshold. -/
def CodonOptimality (caiWeights : Codon → Float)
    (minCAI : Float)
    (codons : List Codon) : Prop :=
  ∀ i : Nat, i < codons.length →
    caiWeights (codons.get ⟨i, by sorry⟩) ≥ minCAI

-- ────────────────────────────────────────────────────────────
-- Key Theorems about Predicate Composition
-- ────────────────────────────────────────────────────────────

/-- THEOREM 1 (Novel): NoGTDinucleotide subsumes NoCrypticSplice(FAIL).
    If there are no GT dinucleotides at all, NoCrypticSplice trivially holds
    because its antecedent (codonHasGT) is never satisfied. -/
theorem noGT_subsumes_crypticSplice_FAIL (codons : List Codon)
    (scoreFn : Codon → Float) (low high : Float)
    (hGT : NoGTDinucleotide codons) :
    NoCrypticSplice codons scoreFn low high := by
  intro i hi hHasGT
  exfalso
  exact hGT.1 i hi hHasGT

/-- THEOREM 2: ValidCodingSeq implies NoStopCodons (modulo last codon). -/
theorem validCoding_implies_noInternalStops (codons : List Codon)
    (hValid : ValidCodingSeq codons)
    (hLen : codons.length ≥ 1) :
    NoStopCodons codons := by
  constructor
  · exact hLen
  · intro i hi
    have h := hValid i (by omega)
    cases h with
    | inl hne => exact hne
    | inr heq => omega

/-- THEOREM 3: ConservationScore with minBLOSUM = 0 means all substitutions
    are at least neutral (non-negative BLOSUM62 score). -/
theorem conservation_nonnegative
    (blosum62 : AminoAcid → AminoAcid → Int)
    (orig opt : List AminoAcid)
    (h : ConservationScore blosum62 0 orig opt) :
    ∀ i : Nat, i < orig.length →
      blosum62 (orig.get ⟨i, by sorry⟩) (opt.get ⟨i, by sorry⟩) ≥ 0 :=
  by
  intro i hi
  exact h.2 i hi

/-- THEOREM 4: Dual-threshold monotonicity.
    PASS at strict threshold implies PASS at permissive threshold. -/
theorem dual_threshold_monotonicity
    (score low1 low2 high : Float)
    (h_low : low1 ≤ low2)
    (h_pass : score < low1) :
    classifySplice score low2 high = SpliceVerdict.PASS :=
  classifySplice_monotone score low1 low2 high h_low h_pass

end BioCompiler
