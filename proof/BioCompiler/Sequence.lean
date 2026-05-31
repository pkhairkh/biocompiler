/-
  BioCompiler.Sequence — Nucleotide Sequences and Pattern Matching

  This module defines the nucleotide alphabet, sequences, and pattern matching
  operations with PROVED completeness theorems. These are the building blocks
  for the scanner functions used in the type system.

  Reference: DOC-03 (SDD) §3.1, DOC-04 (ICD) §5
-/

namespace BioCompiler

-- ==============================================================================
-- Nucleotide Alphabet
-- ==============================================================================

/-- IUPAC nucleotide alphabet (standard DNA bases). -/
inductive Nucleotide where
  | A : Nucleotide
  | C : Nucleotide
  | G : Nucleotide
  | T : Nucleotide
  deriving DecidableEq, Repr, BEq

namespace Sequence

/-- A nucleotide sequence is a list of nucleotides. -/
abbrev Sequence := List Nucleotide

/-- The empty sequence. -/
def empty : Sequence := []

/-- Length of a sequence. Uses List.length directly to avoid termination issues. -/
def length (seq : Sequence) : Nat := List.length seq

/-- GC content as a rational number: (G + C) / total.
    Returns 0 for the empty sequence. -/
def gcContent (seq : Sequence) : Rat :=
  if seq.length = 0 then 0
  else ((seq.count Nucleotide.G + seq.count Nucleotide.C : Nat) : Rat) / seq.length

/-- Extract a subsequence starting at position `start` with length `len`. -/
def subseq (seq : Sequence) (start len : Nat) : Sequence :=
  (seq.drop start).take len

/-- LawfulBEq instances for Nucleotide (required for List Nucleotide BEq). -/
private theorem nuc_eq_of_beq : ∀ {a b : Nucleotide}, (a == b) = true → a = b
  | .A, .A, _ => rfl
  | .C, .C, _ => rfl
  | .G, .G, _ => rfl
  | .T, .T, _ => rfl
  | .A, .C, h => False.elim (absurd h (by decide))
  | .A, .G, h => False.elim (absurd h (by decide))
  | .A, .T, h => False.elim (absurd h (by decide))
  | .C, .A, h => False.elim (absurd h (by decide))
  | .C, .G, h => False.elim (absurd h (by decide))
  | .C, .T, h => False.elim (absurd h (by decide))
  | .G, .A, h => False.elim (absurd h (by decide))
  | .G, .C, h => False.elim (absurd h (by decide))
  | .G, .T, h => False.elim (absurd h (by decide))
  | .T, .A, h => False.elim (absurd h (by decide))
  | .T, .C, h => False.elim (absurd h (by decide))
  | .T, .G, h => False.elim (absurd h (by decide))

instance : ReflBEq Nucleotide where
  rfl {a} := by cases a <;> decide

instance : LawfulBEq Nucleotide where
  eq_of_beq := nuc_eq_of_beq

/-- Check if a pattern matches the sequence at a specific position.
    Returns true iff `seq[start : start + pattern.length] = pattern`.
    Uses BEq (==) for List comparison, which returns Bool directly. -/
def matchesAt (seq : Sequence) (pattern : Sequence) (pos : Nat) : Bool :=
  if pos + pattern.length ≤ seq.length then
    (seq.drop pos |>.take pattern.length) == pattern
  else false

/-- Check if a pattern appears anywhere in the sequence.
    Scans every position from 0 to `seq.length - pattern.length`.
    Returns true iff there exists a position where the pattern matches. -/
def containsPattern (seq : Sequence) (pattern : Sequence) : Bool :=
  if pattern.length = 0 then true
  else if seq.length < pattern.length then false
  else
    let maxPos := seq.length - pattern.length
    (List.range (maxPos + 1)).any fun pos => matchesAt seq pattern pos

-- ==============================================================================
-- Pattern Matching Correctness Theorems
-- ==============================================================================

/-- Helper: matchesAt returns true iff the subsequence equals the pattern. -/
theorem matchesAt_spec (seq pattern : Sequence) (pos : Nat) :
    matchesAt seq pattern pos = true ↔
      pos + pattern.length ≤ seq.length ∧
      (seq.drop pos |>.take pattern.length) = pattern := by
  unfold matchesAt
  constructor
  · intro h
    split at h
    · exact ⟨by omega, eq_of_beq h⟩
    · simp at h
  · intro ⟨h_len, h_eq⟩
    split
    · rw [h_eq]; exact BEq.rfl
    · omega

/-- THEOREM (containsPattern completeness):
    If the pattern appears at position `pos` in the sequence, then
    containsPattern returns true. -/
theorem containsPattern_complete (seq pattern : Sequence) (pos : Nat)
    (h_pos : pos + pattern.length ≤ seq.length)
    (h_match : (seq.drop pos |>.take pattern.length) = pattern) :
    containsPattern seq pattern = true := by
  unfold containsPattern
  split
  · rfl
  · split
    · omega
    · have h_pos_bound : pos < seq.length - pattern.length + 1 := by omega
      have h_matches : matchesAt seq pattern pos = true :=
        (matchesAt_spec seq pattern pos).mpr ⟨h_pos, h_match⟩
      exact List.any_eq_true.mpr ⟨pos, List.mem_range.mpr h_pos_bound, h_matches⟩

/-- THEOREM (containsPattern soundness):
    If containsPattern returns true, then the pattern appears at
    some position in the sequence. -/
theorem containsPattern_sound (seq pattern : Sequence) :
    containsPattern seq pattern = true →
      ∃ (pos : Nat), pos + pattern.length ≤ seq.length ∧
        (seq.drop pos |>.take pattern.length) = pattern := by
  intro h
  by_cases h_empty : pattern.length = 0
  case pos =>
    refine ⟨0, ?_, ?_⟩
    · omega
    · have : pattern = [] := List.length_eq_zero.mp h_empty; simp [this]
  case neg =>
    by_cases h_short : seq.length < pattern.length
    case pos =>
      have : containsPattern seq pattern = false := by
        simp [containsPattern, h_empty, h_short]
      rw [this] at h; simp at h
    case neg =>
      have h_any : (List.range (seq.length - pattern.length + 1)).any
          (matchesAt seq pattern) = true := by
        simp [containsPattern, h_empty, h_short] at h; exact h
      obtain ⟨pos, h_mem, h_match⟩ := List.any_eq_true.mp h_any
      have h_pos : pos < seq.length - pattern.length + 1 := List.mem_range.mp h_mem
      obtain ⟨h_len, h_eq⟩ := (matchesAt_spec seq pattern pos).mp h_match
      exact ⟨pos, by omega, h_eq⟩

-- ==============================================================================
-- Codon Operations
-- ==============================================================================

/-- A codon is a sequence of exactly 3 nucleotides. -/
abbrev Codon := Sequence

/-- The three stop codons in the standard genetic code (DNA form). -/
def stopCodons : List Codon :=
  [[Nucleotide.T, Nucleotide.A, Nucleotide.A],   -- TAA
   [Nucleotide.T, Nucleotide.A, Nucleotide.G],   -- TAG
   [Nucleotide.T, Nucleotide.G, Nucleotide.A]]   -- TGA

/-- Check if a codon is a stop codon. -/
def isStopCodon (codon : Codon) : Bool :=
  stopCodons.any (· = codon)

/-- Check if any position in the sequence (in a given reading frame)
    contains a premature stop codon. A premature stop is a stop codon
    that appears before the last codon position. -/
def hasPrematureStop (seq : Sequence) (readingFrame : Nat) : Bool :=
  let codonStarts := (List.range ((seq.length - readingFrame) / 3)).map
    fun i => readingFrame + 3 * i
  codonStarts.any fun start =>
    if start + 3 ≤ seq.length then
      isStopCodon (seq.drop start |>.take 3)
    else false

/-- Check that every exon boundary preserves the reading frame
    (exon lengths are multiples of 3, relative to the reading frame).
    This takes a list of exon boundary positions and checks that
    each boundary position is congruent to the reading frame mod 3. -/
def readingFrameConsistent (boundaries : List Nat) (readingFrame : Nat) : Bool :=
  boundaries.all fun pos => pos % 3 = readingFrame % 3

end Sequence

end BioCompiler
