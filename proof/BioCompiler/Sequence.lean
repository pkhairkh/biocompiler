/-
  BioCompiler.Sequence — Nucleotide Sequences and Pattern Matching

  Key results:
  - matchesAt_spec: matchesAt returns true iff the subsequence equals the pattern
  - containsPattern_complete: if pattern appears, containsPattern returns true
  - containsPattern_sound: if containsPattern returns true, pattern appears

  Reference: DOC-03 (SDD) §3.1, DOC-04 (ICD) §5
-/

namespace BioCompiler

inductive Nucleotide where
  | A : Nucleotide
  | C : Nucleotide
  | G : Nucleotide
  | T : Nucleotide
  deriving DecidableEq, Repr, BEq

abbrev Sequence := List Nucleotide

namespace Sequence

def empty : Sequence := []

def gcContent (seq : Sequence) : Rat :=
  if seq.length = 0 then 0
  else ((seq.count Nucleotide.G + seq.count Nucleotide.C : Nat) : Rat) / seq.length

/-- Check if a pattern matches the sequence at a specific position (Prop version). -/
def matchesAtProp (seq : Sequence) (pattern : Sequence) (pos : Nat) : Prop :=
  pos + pattern.length ≤ seq.length ∧ (seq.drop pos).take pattern.length = pattern

/-- Check if a pattern matches the sequence at a specific position (Bool version). -/
def matchesAt (seq : Sequence) (pattern : Sequence) (pos : Nat) : Bool :=
  if h : pos + pattern.length ≤ seq.length then
    decide ((seq.drop pos).take pattern.length = pattern)
  else false

theorem matchesAt_spec (seq pattern : Sequence) (pos : Nat) :
    matchesAt seq pattern pos = true ↔ matchesAtProp seq pattern pos := by
  unfold matchesAt matchesAtProp
  split
  · simp only [decide_eq_true_iff]
      -- Goal: decide ... = true ↔ pos + pattern.length ≤ seq.length ∧ ...
      -- But the iff for decide gives us the Prop directly
    constructor
    · intro h; exact ⟨by omega, h⟩
    · intro ⟨h₁, h₂⟩; exact h₂
  · -- else branch: pos + pattern.length > seq.length
    simp
    intro h₁ h₂; omega

def containsPattern (seq : Sequence) (pattern : Sequence) : Bool :=
  if pattern.length = 0 then true
  else if seq.length < pattern.length then false
  else
    let maxPos := seq.length - pattern.length
    (List.range (maxPos + 1)).any fun pos => matchesAt seq pattern pos

theorem containsPattern_complete (seq pattern : Sequence) (pos : Nat)
    (h_pos : pos + pattern.length ≤ seq.length)
    (h_match : (seq.drop pos).take pattern.length = pattern) :
    containsPattern seq pattern = true := by
  unfold containsPattern
  split
  · -- pattern.length = 0
    -- true branch: containsPattern = true regardless
    rfl
  · split
    · -- seq.length < pattern.length
      omega
    · -- main case
      have h_pos_bound : pos < seq.length - pattern.length + 1 := by omega
      have h_matches : matchesAt seq pattern pos = true :=
        (matchesAt_spec seq pattern pos).mpr ⟨by omega, h_match⟩
      rw [List.any_eq_true]
      exact ⟨pos, List.mem_range.mpr h_pos_bound, h_matches⟩

theorem containsPattern_sound (seq pattern : Sequence) :
    containsPattern seq pattern = true →
      ∃ (pos : Nat), pos + pattern.length ≤ seq.length ∧
        (seq.drop pos).take pattern.length = pattern := by
  intro h
  unfold containsPattern at h
  split at h
  · -- pattern.length = 0
    have h_nil : pattern = [] := by
      cases pattern with
      | nil => rfl
      | cons hd tl => simp at *
    refine ⟨0, by omega, by rw [h_nil]; simp⟩
  · split at h
    · -- seq.length < pattern.length: impossible
      -- containsPattern = false but h says true, contradiction
      simp at h
    · -- main case
      obtain ⟨pos, h_mem, h_match⟩ := List.any_eq_true.mp h
      have h_pos : pos < seq.length - pattern.length + 1 := List.mem_range.mp h_mem
      have ⟨h_len, h_eq⟩ := (matchesAt_spec seq pattern pos).mp h_match
      refine ⟨pos, by omega, h_eq⟩

-- ==============================================================================
-- Codon Operations
-- ==============================================================================

abbrev Codon := Sequence

def stopCodons : List Codon :=
  [[Nucleotide.T, Nucleotide.A, Nucleotide.A],
   [Nucleotide.T, Nucleotide.A, Nucleotide.G],
   [Nucleotide.T, Nucleotide.G, Nucleotide.A]]

def isStopCodon (codon : Codon) : Bool :=
  stopCodons.any fun c => decide (c = codon)

def hasPrematureStop (seq : Sequence) (readingFrame : Nat) : Bool :=
  let codonStarts := (List.range ((seq.length - readingFrame) / 3)).map
    fun i => readingFrame + 3 * i
  codonStarts.any fun start =>
    if start + 3 ≤ seq.length then
      isStopCodon ((seq.drop start).take 3)
    else false

def readingFrameConsistent (boundaries : List Nat) (readingFrame : Nat) : Bool :=
  boundaries.all fun pos => pos % 3 = readingFrame % 3

end Sequence

end BioCompiler
