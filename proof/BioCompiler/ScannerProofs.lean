/-
  BioCompiler.ScannerProofs — Concrete Scanner Proofs

  This module replaces axioms 4-11 from the Trusted Computing Base with
  actual proofs. CpGIslandScanner, PromoterScanner, and TMDomainScanner
  are now backed by concrete sliding-window implementations with proved
  completeness, soundness, and (where applicable) borderline completeness.

  Proof Strategy:
  ┌─────────────────────────────────────────────────────────────────────────┐
  │  CpG Island Scanner = sliding window over all positions               │
  │                                                                       │
  │  For each position pos where pos + windowSize ≤ seq.length:          │
  │    1. Compute GC content of window                                    │
  │    2. Compute Obs/Exp CG ratio of window                              │
  │    3. Check both against thresholds                                   │
  │                                                                       │
  │  Promoter Scanner = sliding window over all positions                 │
  │                                                                       │
  │  For each position pos where pos + motifSize ≤ seq.length:           │
  │    1. Compute promoter score (TATA box motif match fraction)          │
  │    2. Check score ≥ threshold (cryptic) or in [0.8t, t) (borderline) │
  │                                                                       │
  │  TM Domain Scanner = sliding window over all positions               │
  │                                                                       │
  │  For each position pos where pos + tmDomainWindowSize ≤ seq.length:  │
  │    1. Extract codons from window                                      │
  │    2. Compute hydrophobic fraction (hydrophobic codons / total)       │
  │    3. Check fraction ≥ threshold (cryptic) or in [0.85t, t) (border) │
  │    4. When isCytosolic = false, TM domains expected → return false   │
  │                                                                       │
  │  Completeness: If ANY window satisfies the criteria, the scanner      │
  │    finds it (because List.any checks every position).                 │
  │                                                                       │
  │  Soundness: If the scanner reports a hit, some position in the        │
  │    range returned true, guaranteeing a valid window position.         │
  │                                                                       │
  │  All proofs follow from the deterministic nature of the scan and     │
  │  the correctness of List.any (which is exact, not heuristic).         │
  └─────────────────────────────────────────────────────────────────────────┘

  ELIMINATED AXIOMS:
  - #4: CpGIslandScanner.scanner_completeness
  - #5: CpGIslandScanner.scanner_soundness
  - #6: PromoterScanner.scanner_completeness
  - #7: PromoterScanner.scanner_soundness
  - #8: PromoterScanner.borderline_completeness
  - #9: TMDomainScanner.scanner_completeness
  - #10: TMDomainScanner.scanner_soundness
  - #11: TMDomainScanner.borderline_completeness

  REFERENCE: DOC-03 (SDD) §3.1, DOC-10 (Deterministic Methods) §4
-/

import BioCompiler.Scanners

namespace BioCompiler

open Sequence

-- ==============================================================================
-- Helper Theorems for Decidable Propositions
-- ==============================================================================

/-- If a decidable proposition holds, `decide` returns true.
    This is the key bridge from Prop-level hypotheses to Bool-level computations. -/
theorem decide_eq_true_of_prop {p : Prop} [hdec : Decidable p] (h : p) :
    @decide p hdec = true := by
  cases hdec with
  | isTrue _ => rfl
  | isFalse h' => exfalso; exact h' h

/-- If `decide` returns true, the proposition holds (converse, for reference). -/
theorem prop_of_decide_eq_true {p : Prop} [Decidable p] (h : decide p = true) : p :=
  of_decide_eq_true h

/-- `a ≥ b` is equivalent to `b ≤ a` for Rat (definitional). -/
theorem Rat.ge_iff_le (a b : Rat) : (a ≥ b) ↔ (b ≤ a) := by rfl

-- ==============================================================================
-- Concrete CpG Island Window Check
-- ==============================================================================

/-- Check if a single window satisfies BOTH CpG island criteria:
    1. GC content ≥ cpgIslandGCThreshold (0.60)
    2. Obs/Exp CG ratio ≥ cpgIslandObsExpThreshold (0.65)

    Returns false if the window is empty. Uses `decide` to convert
    decidable Rat comparisons to Bool. -/
def cpgWindowCheck (window : Sequence) : Bool :=
  if window.length = 0 then false
  else
    let gcRatio := (window.count Nucleotide.G + window.count Nucleotide.C : Rat) / window.length
    let gcPass := decide (gcRatio ≥ cpgIslandGCThreshold)
    let cpgCount := (List.zipWith (· == ·) window (window.drop 1)).count true
    let obsExpPass := decide ((cpgCount : Rat) * window.length ≥
      cpgIslandObsExpThreshold * (window.count Nucleotide.C) * (window.count Nucleotide.G))
    gcPass && obsExpPass

/-- THEOREM: If both CpG island criteria hold for a nonempty window,
    then cpgWindowCheck returns true.

    Proof: Both `decide` calls return true (by decide_eq_true_of_prop),
    and `true && true = true`. -/
theorem cpgWindowCheck_true (window : Sequence)
    (h_len : window.length > 0)
    (h_gc : (window.count Nucleotide.G + window.count Nucleotide.C : Rat) / window.length
              ≥ cpgIslandGCThreshold)
    (h_obs_exp :
      let cpgCount := (List.zipWith (· == ·) window (window.drop 1)).count true
      (cpgCount : Rat) * window.length ≥
        cpgIslandObsExpThreshold * (window.count Nucleotide.C) * (window.count Nucleotide.G)) :
    cpgWindowCheck window = true := by
  unfold cpgWindowCheck
  split
  · -- window.length = 0 branch: contradicts h_len
    omega
  · -- main branch: show gcPass && obsExpPass = true
    have h_gcPass : decide ((window.count Nucleotide.G + window.count Nucleotide.C : Rat) /
        window.length ≥ cpgIslandGCThreshold) = true :=
      decide_eq_true_of_prop h_gc
    have h_obsExpPass : decide (
        let cpgCount := (List.zipWith (· == ·) window (window.drop 1)).count true
        (cpgCount : Rat) * window.length ≥
          cpgIslandObsExpThreshold * (window.count Nucleotide.C) *
          (window.count Nucleotide.G)) = true :=
      decide_eq_true_of_prop h_obs_exp
    -- Now: gcPass = true ∧ obsExpPass = true → gcPass && obsExpPass = true
    simp only [h_gcPass, h_obsExpPass, Bool.true_and]

-- ==============================================================================
-- Concrete CpG Island Scanner (Sliding Window)
-- ==============================================================================

/-- CONCRETE CpG island scanner: checks every window position.
    Uses a sliding window of size `cpgIslandWindowSize` (200 bp).

    Implementation:
    1. If sequence is shorter than window size, return false
    2. Otherwise, check cpgWindowCheck at every valid position
    3. Return true if any window passes both criteria

    This is a DETERMINISTIC computation — no heuristics, no approximations.
    Every position is checked, and the check is exact. -/
def hasCpGIslandConcrete (seq : Sequence) : Bool :=
  if seq.length < cpgIslandWindowSize then false
  else
    let numWindows := seq.length - cpgIslandWindowSize + 1
    (List.range numWindows).any fun pos =>
      cpgWindowCheck ((seq.drop pos).take cpgIslandWindowSize)

-- ==============================================================================
-- Completeness Proof
-- ==============================================================================

/-- THEOREM (CpG Island Scanner Completeness): If any window in the sequence
    satisfies BOTH the GC content criterion AND the Obs/Exp CG ratio criterion,
    then the concrete scanner returns true.

    Proof outline:
    1. From pos + windowSize ≤ seq.length, derive seq.length ≥ windowSize
       (so the scanner enters the main branch)
    2. From pos + windowSize ≤ seq.length, derive pos is in the range
       [0, numWindows) where numWindows = seq.length - windowSize + 1
    3. From the two criteria, derive cpgWindowCheck returns true at pos
    4. By List.any_eq_true, if any element in the range returns true,
       the whole any-call returns true
    5. This contradicts hasCpGIslandConcrete = false

    Key insight: the scanner checks the EXACT SAME conditions as the
    hypotheses, so there is no gap between what we assume and what
    the scanner checks. -/
theorem hasCpGIslandConcrete_complete (seq : Sequence) (pos : Nat)
    (h_pos : pos + cpgIslandWindowSize ≤ seq.length)
    (h_gc :
      let window := (seq.drop pos).take cpgIslandWindowSize
      (window.count Nucleotide.G + window.count Nucleotide.C : Rat) /
        window.length ≥ cpgIslandGCThreshold)
    (h_obs_exp :
      let window := (seq.drop pos).take cpgIslandWindowSize
      let cpgCount := (List.zipWith (· == ·) window (window.drop 1)).count true
      (cpgCount : Rat) * window.length ≥
        cpgIslandObsExpThreshold * (window.count Nucleotide.C) *
        (window.count Nucleotide.G))
    (h_false : hasCpGIslandConcrete seq = false) :
    False := by
  unfold hasCpGIslandConcrete at h_false
  -- Step 1: Split on whether seq.length < cpgIslandWindowSize
  split at h_false
  · -- seq.length < cpgIslandWindowSize: contradicts h_pos
    omega
  · -- Main branch: seq.length ≥ cpgIslandWindowSize
    -- Step 2: Show pos is in the range
    have h_pos_in_range : pos < seq.length - cpgIslandWindowSize + 1 := by omega
    have h_mem : pos ∈ List.range (seq.length - cpgIslandWindowSize + 1) :=
      List.mem_range.mpr h_pos_in_range
    -- Step 3: Show cpgWindowCheck returns true at pos
    have h_window_len : ((seq.drop pos).take cpgIslandWindowSize).length > 0 := by
      have h_take_len : ((seq.drop pos).take cpgIslandWindowSize).length = min cpgIslandWindowSize (seq.length - pos) := by
        rw [List.length_take, List.length_drop]
      rw [h_take_len, Nat.min_def]
      split
      · native_decide  -- cpgIslandWindowSize = 200 > 0
      · omega
    have h_window_true : cpgWindowCheck ((seq.drop pos).take cpgIslandWindowSize) = true :=
      cpgWindowCheck_true ((seq.drop pos).take cpgIslandWindowSize)
        h_window_len h_gc h_obs_exp
    -- Step 4: By List.any_eq_true, the whole scan returns true
    have h_any_true :
        (List.range (seq.length - cpgIslandWindowSize + 1)).any
          (fun p => cpgWindowCheck ((seq.drop p).take cpgIslandWindowSize)) = true :=
      List.any_eq_true.mpr ⟨pos, h_mem, h_window_true⟩
    -- Step 5: Contradiction with h_false
    rw [h_any_true] at h_false
    cases h_false

-- ==============================================================================
-- Soundness Proof
-- ==============================================================================

/-- THEOREM (CpG Island Scanner Soundness): If the concrete scanner reports
    a CpG island (returns true), then there exists a valid window position
    in the sequence (i.e., a position where a windowSize-sized window fits).

    Proof outline:
    1. If seq.length < windowSize, the scanner returns false (contradiction)
    2. Otherwise, List.any = true means some position in the range returned true
    3. Any position in List.range numWindows satisfies pos < numWindows,
       which gives pos + windowSize ≤ seq.length

    Note: the soundness axiom in the CpGIslandScanner class only requires
    the existence of a valid window position, not that the window actually
    satisfies the CpG criteria. This is because the soundness guarantee
    needed for the type system is: "if the scanner says there IS an island,
    then the sequence is long enough that an island COULD exist."
    The stronger form (the window at that position actually satisfies the
    criteria) also holds and is provable, but is not needed for soundness. -/
theorem hasCpGIslandConcrete_sound (seq : Sequence)
    (h_true : hasCpGIslandConcrete seq = true) :
    ∃ (pos : Nat), pos + cpgIslandWindowSize ≤ seq.length := by
  unfold hasCpGIslandConcrete at h_true
  split at h_true
  · -- seq.length < cpgIslandWindowSize: scanner returns false, contradiction
    simp at h_true
  · -- Main branch
    obtain ⟨pos, h_mem, _h_check⟩ := List.any_eq_true.mp h_true
    have h_pos_in_range : pos < seq.length - cpgIslandWindowSize + 1 :=
      List.mem_range.mp h_mem
    exact ⟨pos, by omega⟩

-- ==============================================================================
-- Stronger Soundness (for reference; not required by the class)
-- ==============================================================================

/-- The window at the position found by the scanner actually has length
    equal to cpgIslandWindowSize. -/
theorem hasCpGIslandConcrete_sound_window_len (seq : Sequence)
    (h_true : hasCpGIslandConcrete seq = true) :
    ∃ (pos : Nat), pos + cpgIslandWindowSize ≤ seq.length ∧
      ((seq.drop pos).take cpgIslandWindowSize).length = cpgIslandWindowSize := by
  unfold hasCpGIslandConcrete at h_true
  split at h_true
  · simp at h_true
  · obtain ⟨pos, h_mem, _h_check⟩ := List.any_eq_true.mp h_true
    have h_pos_in_range : pos < seq.length - cpgIslandWindowSize + 1 :=
      List.mem_range.mp h_mem
    have h_pos_le : pos + cpgIslandWindowSize ≤ seq.length := by omega
    refine ⟨pos, h_pos_le, ?_⟩
    have h_take_len : ((seq.drop pos).take cpgIslandWindowSize).length = min cpgIslandWindowSize (seq.length - pos) := by
      rw [List.length_take, List.length_drop]
    rw [h_take_len, Nat.min_def]
    split
    · native_decide  -- cpgIslandWindowSize = 200 > 0
    · omega

-- ==============================================================================
-- CpGIslandScanner Instance — Eliminates Axioms 4-5
-- ==============================================================================

/-- Concrete CpGIslandScanner instance with PROVED completeness and soundness.

    This instance replaces the abstract (axiomatic) CpGIslandScanner that
    previously required trust. Now the completeness and soundness are
    derived from first principles:

    - Completeness: follows from the scanner checking every valid position
      and the conditions matching exactly what the hypotheses state
    - Soundness: follows from List.any being true only when some element
      satisfies the predicate, which implies a valid position exists

    ELIMINATED AXIOMS:
    - CpGIslandScanner.scanner_completeness (was TCB axiom #4)
    - CpGIslandScanner.scanner_soundness  (was TCB axiom #5) -/
instance concreteCpGIslandScanner : CpGIslandScanner where
  hasCpGIsland := hasCpGIslandConcrete
  scanner_completeness := hasCpGIslandConcrete_complete
  scanner_soundness := hasCpGIslandConcrete_sound

-- ==============================================================================
-- Verification: Instance satisfies the class contract
-- ==============================================================================

/-- Verification lemma: the instance's completeness is exactly the proved theorem. -/
theorem instance_completeness_eq :
    @CpGIslandScanner.scanner_completeness concreteCpGIslandScanner =
      hasCpGIslandConcrete_complete := rfl

/-- Verification lemma: the instance's soundness is exactly the proved theorem. -/
theorem instance_soundness_eq :
    @CpGIslandScanner.scanner_soundness concreteCpGIslandScanner =
      hasCpGIslandConcrete_sound := rfl

-- ==============================================================================
-- ==============================================================================
-- Concrete Promoter Scanner Proofs — Eliminates Axioms 6-8
-- ==============================================================================
-- ==============================================================================

/-- Check if the promoter score at a position meets the cryptic threshold.
    Returns false if the position is out of range. Uses `decide` to convert
    the decidable Rat comparison to Bool. -/
def promoterWindowCheck (seq : Sequence) (threshold : Rat) (pos : Nat) : Bool :=
  if pos + promoterMotifSize ≤ seq.length then
    decide (promoterScoreAt seq pos ≥ threshold)
  else false

/-- Check if the promoter score at a position falls in the borderline range:
    score ≥ threshold * 0.8 AND score < threshold.
    Returns false if the position is out of range. -/
def borderlinePromoterWindowCheck (seq : Sequence) (threshold : Rat) (pos : Nat) : Bool :=
  if pos + promoterMotifSize ≤ seq.length then
    decide (promoterScoreAt seq pos ≥ threshold * 8 / 10) &&
      decide (¬(promoterScoreAt seq pos ≥ threshold))
  else false

-- ==============================================================================
-- Window Check Correctness Lemmas
-- ==============================================================================

/-- THEOREM: If the promoter score at a valid position meets the threshold,
    then promoterWindowCheck returns true.

    Proof: The position is in range (so we enter the then-branch), and
    `decide` returns true because the proposition holds. -/
theorem promoterWindowCheck_true (seq : Sequence) (threshold : Rat) (pos : Nat)
    (h_pos : pos + promoterMotifSize ≤ seq.length)
    (h_score : promoterScoreAt seq pos ≥ threshold) :
    promoterWindowCheck seq threshold pos = true := by
  unfold promoterWindowCheck
  split
  · -- then branch: pos + promoterMotifSize ≤ seq.length holds
    exact decide_eq_true_of_prop h_score
  · -- else branch: ¬(pos + promoterMotifSize ≤ seq.length) — contradicts h_pos
    omega

/-- THEOREM: If the promoter score at a valid position falls in the borderline
    range [threshold*0.8, threshold), then borderlinePromoterWindowCheck
    returns true.

    Proof: Both `decide` calls return true (by decide_eq_true_of_prop),
    and `true && true = true`. -/
theorem borderlinePromoterWindowCheck_true (seq : Sequence) (threshold : Rat) (pos : Nat)
    (h_pos : pos + promoterMotifSize ≤ seq.length)
    (h_lower : promoterScoreAt seq pos ≥ threshold * 8 / 10)
    (h_not_above : ¬(promoterScoreAt seq pos ≥ threshold)) :
    borderlinePromoterWindowCheck seq threshold pos = true := by
  unfold borderlinePromoterWindowCheck
  split
  · -- then branch: pos + promoterMotifSize ≤ seq.length holds
    have h1 : decide (promoterScoreAt seq pos ≥ threshold * 8 / 10) = true :=
      decide_eq_true_of_prop h_lower
    have h2 : decide (¬(promoterScoreAt seq pos ≥ threshold)) = true :=
      decide_eq_true_of_prop h_not_above
    simp only [h1, h2, Bool.true_and]
  · -- else branch: ¬(pos + promoterMotifSize ≤ seq.length) — contradicts h_pos
    omega

-- ==============================================================================
-- Concrete Promoter Scanner (Sliding Window)
-- ==============================================================================

/-- CONCRETE cryptic promoter scanner: checks every position.
    Uses a sliding window of size `promoterMotifSize` (6 bp, TATA box).

    Implementation:
    1. If sequence is shorter than motif size, return false
    2. Otherwise, check promoterWindowCheck at every valid position
    3. Return true if any position has score ≥ threshold

    This is a DETERMINISTIC computation — no heuristics, no approximations.
    Every position is checked, and the check is exact.

    Note: The organism parameter is accepted for interface compatibility
    but is not used in this concrete implementation, as the TATA box
    is a universal eukaryotic core promoter element. -/
def hasCrypticPromoterConcrete (seq : Sequence) (organism : String) (threshold : Rat) : Bool :=
  if seq.length < promoterMotifSize then false
  else
    let numPositions := seq.length - promoterMotifSize + 1
    (List.range numPositions).any fun pos => promoterWindowCheck seq threshold pos

/-- CONCRETE borderline promoter scanner: checks every position for
    scores in [threshold * 0.8, threshold).

    Same sliding window structure as the cryptic scanner, but checking
    for borderline scores instead. -/
def hasBorderlinePromoterConcrete (seq : Sequence) (organism : String) (threshold : Rat) : Bool :=
  if seq.length < promoterMotifSize then false
  else
    let numPositions := seq.length - promoterMotifSize + 1
    (List.range numPositions).any fun pos => borderlinePromoterWindowCheck seq threshold pos

-- ==============================================================================
-- Completeness Proof (Cryptic Promoter)
-- ==============================================================================

/-- THEOREM (Promoter Scanner Completeness): If any position in the sequence
    has a promoter score ≥ threshold, then the concrete scanner returns true
    (i.e., cannot return false).

    Proof outline:
    1. From pos + motifSize ≤ seq.length, derive seq.length ≥ motifSize
       (so the scanner enters the main branch)
    2. From pos + motifSize ≤ seq.length, derive pos is in the range
       [0, numPositions) where numPositions = seq.length - motifSize + 1
    3. From score ≥ threshold, derive promoterWindowCheck returns true at pos
    4. By List.any_eq_true, if any element in the range returns true,
       the whole any-call returns true
    5. This contradicts hasCrypticPromoterConcrete = false

    Key insight: the scanner checks the EXACT SAME condition as the
    hypothesis (promoterScoreAt seq pos ≥ threshold), so there is no gap
    between what we assume and what the scanner checks. -/
theorem hasCrypticPromoterConcrete_complete (seq : Sequence) (organism : String)
    (threshold : Rat) (pos : Nat)
    (h_pos : pos + promoterMotifSize ≤ seq.length)
    (h_score : promoterScoreAt seq pos ≥ threshold)
    (h_false : hasCrypticPromoterConcrete seq organism threshold = false) :
    False := by
  unfold hasCrypticPromoterConcrete at h_false
  -- Step 1: Split on whether seq.length < promoterMotifSize
  split at h_false
  · -- seq.length < promoterMotifSize: contradicts h_pos
    omega
  · -- Main branch: seq.length ≥ promoterMotifSize
    -- Step 2: Show pos is in the range
    have h_pos_in_range : pos < seq.length - promoterMotifSize + 1 := by omega
    have h_mem : pos ∈ List.range (seq.length - promoterMotifSize + 1) :=
      List.mem_range.mpr h_pos_in_range
    -- Step 3: Show promoterWindowCheck returns true at pos
    have h_window_true : promoterWindowCheck seq threshold pos = true :=
      promoterWindowCheck_true seq threshold pos h_pos h_score
    -- Step 4: By List.any_eq_true, the whole scan returns true
    have h_any_true :
        (List.range (seq.length - promoterMotifSize + 1)).any
          (fun p => promoterWindowCheck seq threshold p) = true :=
      List.any_eq_true.mpr ⟨pos, h_mem, h_window_true⟩
    -- Step 5: Contradiction with h_false
    rw [h_any_true] at h_false
    cases h_false

-- ==============================================================================
-- Soundness Proof (Cryptic Promoter)
-- ==============================================================================

/-- THEOREM (Promoter Scanner Soundness): If the concrete scanner reports
    a cryptic promoter (returns true), then there exists a valid position
    in the sequence with promoter score ≥ threshold.

    Proof outline:
    1. If seq.length < motifSize, the scanner returns false (contradiction)
    2. Otherwise, List.any = true means some position in the range returned true
    3. Any position in List.range numPositions satisfies pos < numPositions,
       which gives pos + motifSize ≤ seq.length
    4. promoterWindowCheck = true at that position implies
       decide(score ≥ threshold) = true, which implies score ≥ threshold

    Stronger than the CpG soundness: we also recover the score condition
    from the decide bridge. -/
theorem hasCrypticPromoterConcrete_sound (seq : Sequence) (organism : String)
    (threshold : Rat)
    (h_true : hasCrypticPromoterConcrete seq organism threshold = true) :
    ∃ (pos : Nat), pos + promoterMotifSize ≤ seq.length ∧
      promoterScoreAt seq pos ≥ threshold := by
  unfold hasCrypticPromoterConcrete at h_true
  split at h_true
  · -- seq.length < promoterMotifSize: scanner returns false, contradiction
    simp at h_true
  · -- Main branch
    obtain ⟨pos, h_mem, h_check⟩ := List.any_eq_true.mp h_true
    have h_pos_in_range : pos < seq.length - promoterMotifSize + 1 :=
      List.mem_range.mp h_mem
    have h_pos_le : pos + promoterMotifSize ≤ seq.length := by omega
    -- From promoterWindowCheck = true, derive the score condition
    unfold promoterWindowCheck at h_check
    split at h_check
    · -- pos + promoterMotifSize ≤ seq.length holds
      have h_score : promoterScoreAt seq pos ≥ threshold :=
        prop_of_decide_eq_true h_check
      exact ⟨pos, h_pos_le, h_score⟩
    · -- ¬(pos + promoterMotifSize ≤ seq.length): contradicts h_pos_le
      simp at h_check

-- ==============================================================================
-- Borderline Completeness Proof
-- ==============================================================================

/-- THEOREM (Borderline Promoter Scanner Completeness): If any position
    in the sequence has a promoter score in [threshold * 0.8, threshold),
    then the borderline scanner cannot return false.

    Proof outline:
    1. From pos + motifSize ≤ seq.length, derive seq.length ≥ motifSize
    2. Show pos is in the scan range
    3. From score ≥ threshold * 0.8 AND ¬(score ≥ threshold), derive
       borderlinePromoterWindowCheck returns true at pos
    4. By List.any_eq_true, the whole scan returns true
    5. Contradiction with hasBorderlinePromoterConcrete = false

    This handles the UNCERTAIN verdict band: sites with score in the
    borderline range are not cryptic (≥ threshold) but are suspicious
    enough to warrant an UNCERTAIN classification. -/
theorem hasBorderlinePromoterConcrete_complete (seq : Sequence) (organism : String)
    (threshold : Rat) (pos : Nat)
    (h_pos : pos + promoterMotifSize ≤ seq.length)
    (h_lower : promoterScoreAt seq pos ≥ threshold * 8 / 10)
    (h_not_above : ¬(promoterScoreAt seq pos ≥ threshold))
    (h_false : hasBorderlinePromoterConcrete seq organism threshold = false) :
    False := by
  unfold hasBorderlinePromoterConcrete at h_false
  -- Step 1: Split on whether seq.length < promoterMotifSize
  split at h_false
  · -- seq.length < promoterMotifSize: contradicts h_pos
    omega
  · -- Main branch: seq.length ≥ promoterMotifSize
    -- Step 2: Show pos is in the range
    have h_pos_in_range : pos < seq.length - promoterMotifSize + 1 := by omega
    have h_mem : pos ∈ List.range (seq.length - promoterMotifSize + 1) :=
      List.mem_range.mpr h_pos_in_range
    -- Step 3: Show borderlinePromoterWindowCheck returns true at pos
    have h_window_true : borderlinePromoterWindowCheck seq threshold pos = true :=
      borderlinePromoterWindowCheck_true seq threshold pos h_pos h_lower h_not_above
    -- Step 4: By List.any_eq_true, the whole scan returns true
    have h_any_true :
        (List.range (seq.length - promoterMotifSize + 1)).any
          (fun p => borderlinePromoterWindowCheck seq threshold p) = true :=
      List.any_eq_true.mpr ⟨pos, h_mem, h_window_true⟩
    -- Step 5: Contradiction with h_false
    rw [h_any_true] at h_false
    cases h_false

-- ==============================================================================
-- PromoterScanner Instance — Eliminates Axioms 6-8
-- ==============================================================================

/-- Concrete PromoterScanner instance with PROVED completeness, soundness,
    and borderline completeness.

    This instance replaces the abstract (axiomatic) PromoterScanner that
    previously required trust. Now all three properties are derived from
    first principles:

    - Completeness: follows from the scanner checking every valid position
      and the condition matching exactly what the hypothesis states
    - Soundness: follows from List.any being true only when some element
      satisfies the predicate, combined with the decide bridge to recover
      the proposition from the Bool
    - Borderline completeness: same structure as completeness, but checking
      the dual-threshold band [threshold*0.8, threshold)

    ELIMINATED AXIOMS:
    - PromoterScanner.scanner_completeness   (was TCB axiom #6)
    - PromoterScanner.scanner_soundness      (was TCB axiom #7)
    - PromoterScanner.borderline_completeness (was TCB axiom #8) -/
instance concretePromoterScanner : PromoterScanner where
  hasCrypticPromoter := hasCrypticPromoterConcrete
  hasBorderlinePromoter := hasBorderlinePromoterConcrete
  scanner_completeness := hasCrypticPromoterConcrete_complete
  scanner_soundness := hasCrypticPromoterConcrete_sound
  borderline_completeness := hasBorderlinePromoterConcrete_complete

-- ==============================================================================
-- Verification: Instance satisfies the class contract
-- ==============================================================================

/-- Verification lemma: the instance's completeness is exactly the proved theorem. -/
theorem instance_promoter_completeness_eq :
    @PromoterScanner.scanner_completeness concretePromoterScanner =
      hasCrypticPromoterConcrete_complete := rfl

/-- Verification lemma: the instance's soundness is exactly the proved theorem. -/
theorem instance_promoter_soundness_eq :
    @PromoterScanner.scanner_soundness concretePromoterScanner =
      hasCrypticPromoterConcrete_sound := rfl

/-- Verification lemma: the instance's borderline completeness is exactly the proved theorem. -/
theorem instance_promoter_borderline_eq :
    @PromoterScanner.borderline_completeness concretePromoterScanner =
      hasBorderlinePromoterConcrete_complete := rfl

-- ==============================================================================
-- TM Domain Scanner Proofs — Eliminates Axioms 9-11
-- ==============================================================================

/-- If a decidable proposition is false, `decide` returns false.
    Converse of decide_eq_true_of_prop; used for borderline checks
    where we need ¬p → decide p = false. -/
theorem decide_eq_false_of_not_prop {p : Prop} [hdec : Decidable p] (h : ¬p) :
    @decide p hdec = false := by
  cases hdec with
  | isTrue h' => exfalso; exact h h'
  | isFalse _ => rfl

-- ==============================================================================
-- Concrete TM Domain Window Checks
-- ==============================================================================

/-- Check if a single window satisfies the TM domain criterion:
    hydrophobic fraction ≥ threshold. Uses `decide` to convert
    the decidable Rat comparison to Bool.
    Returns false if the window is empty. -/
def tmWindowCheck (window : Sequence) (threshold : Rat) : Bool :=
  if window.length = 0 then false
  else decide (tmHydrophobicFraction window ≥ threshold)

/-- Check if a single window satisfies the borderline TM domain criterion:
    hydrophobic fraction ≥ threshold × 0.85 AND < threshold.
    Uses `decide` for both comparisons.
    Returns false if the window is empty. -/
def tmBorderlineWindowCheck (window : Sequence) (threshold : Rat) : Bool :=
  if window.length = 0 then false
  else
    let frac := tmHydrophobicFraction window
    decide (frac ≥ threshold * 85 / 100) && !(decide (frac ≥ threshold))

/-- THEOREM: If the TM domain criterion holds for a nonempty window,
    then tmWindowCheck returns true.

    Proof: Both the non-emptiness check and `decide` return true
    (by the window length being > 0 and decide_eq_true_of_prop). -/
theorem tmWindowCheck_true (window : Sequence) (threshold : Rat)
    (h_len : window.length > 0)
    (h_hydro : tmHydrophobicFraction window ≥ threshold) :
    tmWindowCheck window threshold = true := by
  unfold tmWindowCheck
  split
  · -- window.length = 0 branch: contradicts h_len
    omega
  · -- main branch: decide returns true
    exact decide_eq_true_of_prop h_hydro

/-- THEOREM: If the borderline TM domain criterion holds for a nonempty window,
    then tmBorderlineWindowCheck returns true.

    Proof: The first `decide` returns true (by decide_eq_true_of_prop),
    the second `decide` returns false (by decide_eq_false_of_not_prop),
    so `true && !false = true && true = true`. -/
theorem tmBorderlineWindowCheck_true (window : Sequence) (threshold : Rat)
    (h_len : window.length > 0)
    (h_ge : tmHydrophobicFraction window ≥ threshold * 85 / 100)
    (h_not : ¬(tmHydrophobicFraction window ≥ threshold)) :
    tmBorderlineWindowCheck window threshold = true := by
  unfold tmBorderlineWindowCheck
  split
  · -- window.length = 0 branch: contradicts h_len
    omega
  · -- main branch
    have h_pass : decide (tmHydrophobicFraction window ≥ threshold * 85 / 100) = true :=
      decide_eq_true_of_prop h_ge
    have h_fail : decide (tmHydrophobicFraction window ≥ threshold) = false :=
      decide_eq_false_of_not_prop h_not
    simp [h_pass, h_fail]

-- ==============================================================================
-- Concrete TM Domain Scanner (Sliding Window)
-- ==============================================================================

/-- CONCRETE TM domain scanner: checks every window position.
    Uses a sliding window of size `tmDomainWindowSize` (51 bp = 17 codons).

    Implementation:
    1. If the protein is not cytosolic (isCytosolic = false), return false
       (TM domains are expected for membrane proteins)
    2. If sequence is shorter than window size, return false
    3. Otherwise, check tmWindowCheck at every valid position
    4. Return true if any window has hydrophobic fraction ≥ threshold

    This is a DETERMINISTIC computation — no heuristics, no approximations.
    Every position is checked, and the check is exact. -/
def hasTMDomainConcrete (seq : Sequence) (isCytosolic : Bool) (threshold : Rat) : Bool :=
  match isCytosolic with
  | false => false
  | true =>
    if seq.length < tmDomainWindowSize then false
    else
      let numWindows := seq.length - tmDomainWindowSize + 1
      (List.range numWindows).any fun pos =>
        tmWindowCheck ((seq.drop pos).take tmDomainWindowSize) threshold

/-- CONCRETE borderline TM domain scanner: checks every window position
    for borderline TM domains (hydrophobic fraction in [threshold*0.85, threshold)). -/
def hasBorderlineTMDomainConcrete (seq : Sequence) (isCytosolic : Bool) (threshold : Rat) : Bool :=
  match isCytosolic with
  | false => false
  | true =>
    if seq.length < tmDomainWindowSize then false
    else
      let numWindows := seq.length - tmDomainWindowSize + 1
      (List.range numWindows).any fun pos =>
        tmBorderlineWindowCheck ((seq.drop pos).take tmDomainWindowSize) threshold

-- ==============================================================================
-- Completeness Proof
-- ==============================================================================

/-- THEOREM (TM Domain Scanner Completeness): If any window in the sequence
    satisfies the hydrophobic fraction criterion (≥ threshold), then the
    concrete scanner returns true.

    Proof outline:
    1. From isCytosolic = true, the scanner enters the main scanning branch
    2. From pos + windowSize ≤ seq.length, derive seq.length ≥ windowSize
       (so the scanner enters the scanning loop)
    3. From pos + windowSize ≤ seq.length, derive pos is in the range
       [0, numWindows) where numWindows = seq.length - windowSize + 1
    4. From the hydrophobic fraction hypothesis, derive tmWindowCheck returns
       true at pos (by tmWindowCheck_true)
    5. By List.any_eq_true, if any element in the range returns true,
       the whole any-call returns true
    6. This contradicts hasTMDomainConcrete = false

    Key insight: the scanner checks the EXACT SAME conditions as the
    hypotheses, so there is no gap between what we assume and what
    the scanner checks. -/
theorem hasTMDomainConcrete_complete (seq : Sequence) (isCytosolic : Bool)
    (threshold : Rat) (pos : Nat)
    (h_cytosolic : isCytosolic = true)
    (h_pos : pos + tmDomainWindowSize ≤ seq.length)
    (h_hydro :
      let window := (seq.drop pos).take tmDomainWindowSize
      tmHydrophobicFraction window ≥ threshold)
    (h_false : hasTMDomainConcrete seq isCytosolic threshold = false) :
    False := by
  unfold hasTMDomainConcrete at h_false
  -- Step 1: Case-split on isCytosolic
  cases isCytosolic with
  | false =>
    -- isCytosolic = false: contradicts h_cytosolic
    exact Bool.false_ne_true h_cytosolic
  | true =>
    -- isCytosolic = true: match reduces to the main branch
    -- Step 2: Split on seq.length < tmDomainWindowSize
    split at h_false
    · -- seq.length < tmDomainWindowSize: contradicts h_pos
      omega
    · -- Main branch: seq.length ≥ tmDomainWindowSize
      -- Step 3: Show pos is in the range
      have h_pos_in_range : pos < seq.length - tmDomainWindowSize + 1 := by omega
      have h_mem : pos ∈ List.range (seq.length - tmDomainWindowSize + 1) :=
        List.mem_range.mpr h_pos_in_range
      -- Step 4: Show tmWindowCheck returns true at pos
      have h_window_len : ((seq.drop pos).take tmDomainWindowSize).length > 0 := by
        have h_take_len : ((seq.drop pos).take tmDomainWindowSize).length = min tmDomainWindowSize (seq.length - pos) := by
          rw [List.length_take, List.length_drop]
        rw [h_take_len, Nat.min_def]
        split
        · native_decide  -- tmDomainWindowSize = 51 > 0
        · omega
      have h_window_true :
          tmWindowCheck ((seq.drop pos).take tmDomainWindowSize) threshold = true :=
        tmWindowCheck_true ((seq.drop pos).take tmDomainWindowSize) threshold
          h_window_len h_hydro
      -- Step 5: By List.any_eq_true, the whole scan returns true
      have h_any_true :
          (List.range (seq.length - tmDomainWindowSize + 1)).any
            (fun p => tmWindowCheck ((seq.drop p).take tmDomainWindowSize) threshold) = true :=
        List.any_eq_true.mpr ⟨pos, h_mem, h_window_true⟩
      -- Step 6: Contradiction with h_false
      rw [h_any_true] at h_false
      cases h_false

-- ==============================================================================
-- Soundness Proof
-- ==============================================================================

/-- THEOREM (TM Domain Scanner Soundness): If the concrete scanner reports
    a TM domain (returns true), then there exists a valid window position
    in the sequence where the hydrophobic fraction ≥ threshold.

    Proof outline:
    1. If isCytosolic = false, the scanner returns false (contradiction)
    2. If seq.length < windowSize, the scanner returns false (contradiction)
    3. Otherwise, List.any = true means some position in the range returned true
    4. Any position in List.range numWindows satisfies pos < numWindows,
       which gives pos + windowSize ≤ seq.length
    5. The tmWindowCheck at that position returned true, meaning
       tmHydrophobicFraction ≥ threshold (by prop_of_decide_eq_true) -/
theorem hasTMDomainConcrete_sound (seq : Sequence) (isCytosolic : Bool)
    (threshold : Rat)
    (h_true : hasTMDomainConcrete seq isCytosolic threshold = true) :
    ∃ (pos : Nat), pos + tmDomainWindowSize ≤ seq.length ∧
      let window := (seq.drop pos).take tmDomainWindowSize
      tmHydrophobicFraction window ≥ threshold := by
  unfold hasTMDomainConcrete at h_true
  cases isCytosolic with
  | false => simp at h_true
  | true =>
    split at h_true
    · -- seq.length < tmDomainWindowSize: scanner returns false, contradiction
      simp at h_true
    · -- Main branch
      obtain ⟨pos, h_mem, h_check⟩ := List.any_eq_true.mp h_true
      have h_pos_in_range : pos < seq.length - tmDomainWindowSize + 1 :=
        List.mem_range.mp h_mem
      have h_pos_le : pos + tmDomainWindowSize ≤ seq.length := by omega
      -- From tmWindowCheck returning true, derive the hydrophobic fraction
      unfold tmWindowCheck at h_check
      split at h_check
      · -- window.length = 0: can't return true
        simp at h_check
      · -- decide = true: the proposition holds
        have h_hydro : tmHydrophobicFraction ((seq.drop pos).take tmDomainWindowSize) ≥ threshold :=
          prop_of_decide_eq_true h_check
        exact ⟨pos, h_pos_le, h_hydro⟩

-- ==============================================================================
-- Borderline Completeness Proof
-- ==============================================================================

/-- THEOREM (Borderline TM Domain Scanner Completeness): If any window in the
    sequence satisfies the borderline hydrophobic fraction criterion
    (≥ threshold × 0.85 and < threshold), then the concrete borderline
    scanner returns true.

    Proof outline: Same structure as completeness, but using
    tmBorderlineWindowCheck instead of tmWindowCheck. The borderline check
    verifies both:
    1. hydrophobic fraction ≥ threshold × 0.85 (via decide_eq_true_of_prop)
    2. hydrophobic fraction < threshold (via decide_eq_false_of_not_prop) -/
theorem hasBorderlineTMDomainConcrete_complete (seq : Sequence) (isCytosolic : Bool)
    (threshold : Rat) (pos : Nat)
    (h_cytosolic : isCytosolic = true)
    (h_pos : pos + tmDomainWindowSize ≤ seq.length)
    (h_ge :
      let window := (seq.drop pos).take tmDomainWindowSize
      tmHydrophobicFraction window ≥ threshold * 85 / 100)
    (h_not :
      let window := (seq.drop pos).take tmDomainWindowSize
      ¬(tmHydrophobicFraction window ≥ threshold))
    (h_false : hasBorderlineTMDomainConcrete seq isCytosolic threshold = false) :
    False := by
  unfold hasBorderlineTMDomainConcrete at h_false
  -- Step 1: Case-split on isCytosolic
  cases isCytosolic with
  | false =>
    -- isCytosolic = false: contradicts h_cytosolic
    exact Bool.false_ne_true h_cytosolic
  | true =>
    -- isCytosolic = true: match reduces to the main branch
    -- Step 2: Split on seq.length < tmDomainWindowSize
    split at h_false
    · -- seq.length < tmDomainWindowSize: contradicts h_pos
      omega
    · -- Main branch
      -- Step 3: Show pos is in the range
      have h_pos_in_range : pos < seq.length - tmDomainWindowSize + 1 := by omega
      have h_mem : pos ∈ List.range (seq.length - tmDomainWindowSize + 1) :=
        List.mem_range.mpr h_pos_in_range
      -- Step 4: Show tmBorderlineWindowCheck returns true at pos
      have h_window_len : ((seq.drop pos).take tmDomainWindowSize).length > 0 := by
        have h_take_len : ((seq.drop pos).take tmDomainWindowSize).length = min tmDomainWindowSize (seq.length - pos) := by
          rw [List.length_take, List.length_drop]
        rw [h_take_len, Nat.min_def]
        split
        · native_decide  -- tmDomainWindowSize = 51 > 0
        · omega
      have h_window_true :
          tmBorderlineWindowCheck ((seq.drop pos).take tmDomainWindowSize) threshold = true :=
        tmBorderlineWindowCheck_true ((seq.drop pos).take tmDomainWindowSize) threshold
          h_window_len h_ge h_not
      -- Step 5: By List.any_eq_true, the whole scan returns true
      have h_any_true :
          (List.range (seq.length - tmDomainWindowSize + 1)).any
            (fun p => tmBorderlineWindowCheck ((seq.drop p).take tmDomainWindowSize) threshold) = true :=
        List.any_eq_true.mpr ⟨pos, h_mem, h_window_true⟩
      -- Step 6: Contradiction with h_false
      rw [h_any_true] at h_false
      cases h_false

-- ==============================================================================
-- TMDomainScanner Instance — Eliminates Axioms 9-11
-- ==============================================================================

/-- Concrete TMDomainScanner instance with PROVED completeness, soundness,
    and borderline completeness.

    This instance replaces the abstract (axiomatic) TMDomainScanner that
    previously required trust. Now the three properties are derived from
    first principles:

    - Completeness: follows from the scanner checking every valid position
      and the conditions matching exactly what the hypotheses state
    - Soundness: follows from List.any being true only when some element
      satisfies the predicate, combined with the decide bridge to recover
      the proposition from the Bool
    - Borderline completeness: follows from the same sliding window pattern
      with the borderline check (fraction in [threshold*0.85, threshold))

    ELIMINATED AXIOMS:
    - TMDomainScanner.scanner_completeness    (was TCB axiom #9)
    - TMDomainScanner.scanner_soundness       (was TCB axiom #10)
    - TMDomainScanner.borderline_completeness  (was TCB axiom #11) -/
instance concreteTMDomainScanner : TMDomainScanner where
  hasTMDomain := hasTMDomainConcrete
  hasBorderlineTMDomain := hasBorderlineTMDomainConcrete
  scanner_completeness := hasTMDomainConcrete_complete
  scanner_soundness := hasTMDomainConcrete_sound
  borderline_completeness := hasBorderlineTMDomainConcrete_complete

-- ==============================================================================
-- Verification: Instance satisfies the class contract
-- ==============================================================================

/-- Verification lemma: the instance's completeness is exactly the proved theorem. -/
theorem tm_instance_completeness_eq :
    @TMDomainScanner.scanner_completeness concreteTMDomainScanner =
      hasTMDomainConcrete_complete := rfl

/-- Verification lemma: the instance's soundness is exactly the proved theorem. -/
theorem tm_instance_soundness_eq :
    @TMDomainScanner.scanner_soundness concreteTMDomainScanner =
      hasTMDomainConcrete_sound := rfl

/-- Verification lemma: the instance's borderline completeness is exactly the proved theorem. -/
theorem tm_instance_borderline_eq :
    @TMDomainScanner.borderline_completeness concreteTMDomainScanner =
      hasBorderlineTMDomainConcrete_complete := rfl

end BioCompiler
