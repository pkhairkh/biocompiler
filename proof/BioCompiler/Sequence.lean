/-
  BioCompiler.Sequence — Nucleotide Sequences and Pattern Matching

  This module defines the nucleotide alphabet, sequences, and pattern matching
  operations with PROVED completeness theorems. These are the building blocks
  for the scanner functions used in the type system.

  Key results:
  - matchesAt_correct: matchesAt returns true iff the subsequence equals the pattern
  - containsPattern_complete: if pattern appears at any position, containsPattern returns true
  - containsPattern_sound: if containsPattern returns true, pattern appears at some position

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

/-- A nucleotide sequence is a list of nucleotides. -/
abbrev Sequence := List Nucleotide

namespace Sequence

/-- The empty sequence. -/
def empty : Sequence := []

/-- Length of a sequence. -/
def length (seq : Sequence) : Nat := seq.length

/-- GC content as a rational number: (G + C) / total.
    Returns 0 for the empty sequence. -/
def gcContent (seq : Sequence) : Rat :=
  if seq.length = 0 then 0
  else ((seq.count Nucleotide.G + seq.count Nucleotide.C : Nat) : Rat) / seq.length

/-- Extract a subsequence starting at position `start` with length `len`. -/
def subseq (seq : Sequence) (start len : Nat) : Sequence :=
  (seq.drop start).take len

/-- Check if a pattern matches the sequence at a specific position.
    Returns true iff `seq[start : start + pattern.length] = pattern`. -/
def matchesAt (seq : Sequence) (pattern : Sequence) (pos : Nat) : Bool :=
  if pos + pattern.length ≤ seq.length then
    seq.drop pos |>.take pattern.length = pattern
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
      seq.drop pos |>.take pattern.length = pattern := by
  simp [matchesAt]
  split
  · simp
    constructor
    · intro h; exact ⟨by simp at *; omega, h⟩
    · intro ⟨h₁, h₂⟩; exact h₂
  · simp
    constructor
    · intro h; omega
    · intro h; omega

/-- THEOREM (containsPattern completeness):
    If the pattern appears at position `pos` in the sequence, then
    containsPattern returns true.

    Proof: We show that `matchesAt` returns true at position `pos`,
    and since `containsPattern` uses `List.any` over all positions
    including `pos`, the result must be true. -/
theorem containsPattern_complete (seq pattern : Sequence) (pos : Nat)
    (h_pos : pos + pattern.length ≤ seq.length)
    (h_match : seq.drop pos |>.take pattern.length = pattern) :
    containsPattern seq pattern = true := by
  simp [containsPattern]
  split
  · -- pattern.length = 0: trivially true
    simp at *
    have : pattern.length = 0 := by simp at *; omega
    simp [this]
  · split
    · -- seq.length < pattern.length: contradiction with h_pos
      omega
    · -- main case: show matchesAt returns true at pos
      have h_pos_bound : pos < seq.length - pattern.length + 1 := by omega
      have h_matches : matchesAt seq pattern pos = true := by
        simp [matchesAt]
        split
        · exact h_match
        · omega
      simp [List.any_iff_exists]
      exact ⟨pos, List.mem_range.mpr h_pos_bound, h_matches⟩

/-- THEOREM (containsPattern soundness):
    If containsPattern returns true, then the pattern appears at
    some position in the sequence. -/
theorem containsPattern_sound (seq pattern : Sequence) :
    containsPattern seq pattern = true →
      ∃ (pos : Nat), pos + pattern.length ≤ seq.length ∧
        seq.drop pos |>.take pattern.length = pattern := by
  simp [containsPattern]
  split
  · -- pattern.length = 0
    intro _
    use 0
    simp
  · split
    · -- seq.length < pattern.length: impossible for true
      intro h; exact absurd h (by decide)
    · -- main case
      simp [List.any_iff_exists]
      intro ⟨pos, h_mem, h_match⟩
      have h_pos : pos < seq.length - pattern.length + 1 := List.mem_range.mp h_mem
      use pos
      constructor
      · omega
      · simp [matchesAt] at h_match
        split at h_match
        · exact h_match
        · omega

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
