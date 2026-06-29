/-
  BioCompiler.BLOSUM62 — Standard BLOSUM62 Substitution Matrix

  This module formalizes the standard BLOSUM62 amino-acid substitution
  matrix (Henikoff & Henikoff, 1992, PNAS 89(22):10915-9) as a Lean 4
  function `blosum62 : AA → AA → Int`, where `AA` is the inductive type
  of the 20 standard amino acids.

  Mathematical properties proved from the matrix definition:
    1. `blosum62_symm`           : the matrix is symmetric
    2. `blosum62_diag_nonneg`    : every diagonal entry is non-negative
    3. `blosum62_diag_pos`       : every diagonal entry is strictly positive
    4. `blosum62_diag_min`       : every diagonal entry is ≥ 4 (smallest is 4)
    5. `blosum62_diag_max`       : every diagonal entry is ≤ 11 (largest is 11)
    6. `blosum62_min`            : every entry is ≥ -4 (smallest off-diagonal)
    7. `blosum62_max`            : every entry is ≤ 11 (largest diagonal)
    8. `blosum62_diag_self_sub_nonneg` : self-substitution scores are non-negative

  Constants:
    - `minDiagScore : Int := 4`   (smallest BLOSUM62 diagonal value)
    - `maxDiagScore : Int := 11`  (largest BLOSUM62 diagonal value)
    - `minScore     : Int := -4`  (smallest BLOSUM62 entry overall)
    - `maxScore     : Int := 11`  (largest BLOSUM62 entry overall)

  These theorems discharge the BLOSUM62-related proof obligations in
  `SLOTVerification.lean` (ConservationScore slotPropertySemantics),
  eliminating the former `blosum62_tool_sound` axiom and the
  `blosum62Score "" ""` placeholder.

  Reference: DOC-11 (Formal Soundness Proof) §SLOT-Verification.
-/

import BioCompiler.Sequence

namespace BioCompiler

/-- The 20 standard amino acids (excludes Stop).

    Constructor order follows the standard IUPAC single-letter code
    alphabetical order used by the canonical BLOSUM62 matrix layout:
    A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V. -/
inductive AA where
  | Ala | Arg | Asn | Asp | Cys
  | Gln | Glu | Gly | His | Ile
  | Leu | Lys | Met | Phe | Pro
  | Ser | Thr | Trp | Tyr | Val
  deriving Repr, DecidableEq, BEq

namespace BLOSUM62

/-- The 20 amino acids in canonical BLOSUM62 ordering. -/
def aaList : List AA :=
  [AA.Ala, AA.Arg, AA.Asn, AA.Asp, AA.Cys,
   AA.Gln, AA.Glu, AA.Gly, AA.His, AA.Ile,
   AA.Leu, AA.Lys, AA.Met, AA.Phe, AA.Pro,
   AA.Ser, AA.Thr, AA.Trp, AA.Tyr, AA.Val]

/-- Index of an amino acid in `aaList` (0..19). -/
def aaIndex (a : AA) : Nat :=
  match a with
  | AA.Ala => 0  | AA.Arg => 1  | AA.Asn => 2  | AA.Asp => 3  | AA.Cys => 4
  | AA.Gln => 5  | AA.Glu => 6  | AA.Gly => 7  | AA.His => 8  | AA.Ile => 9
  | AA.Leu => 10 | AA.Lys => 11 | AA.Met => 12 | AA.Phe => 13 | AA.Pro => 14
  | AA.Ser => 15 | AA.Thr => 16 | AA.Trp => 17 | AA.Tyr => 18 | AA.Val => 19

/-- Standard BLOSUM62 20×20 substitution matrix (Henikoff & Henikoff, 1992).

    Row and column order: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S,
    T, W, Y, V.  All entries are integers (log-odds substitution scores
    scaled by 1/2 bit). -/
def matrix : List (List Int) :=
  [ --  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0], -- A
    [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3], -- R
    [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3], -- N
    [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3], -- D
    [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], -- C
    [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2], -- Q
    [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2], -- E
    [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3], -- G
    [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3], -- H
    [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3], -- I
    [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1], -- L
    [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2], -- K
    [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1], -- M
    [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1], -- F
    [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2], -- P
    [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2], -- S
    [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0], -- T
    [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3], -- W
    [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1], -- Y
    [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4]  -- V
  ]

/-- The BLOSUM62 substitution score for an amino acid pair.

    Property: the matrix is symmetric (`blosum62_symm`), so the order of
    arguments is irrelevant.  All entries are integers in `[-4, 11]`
    (`blosum62_min`, `blosum62_max`). -/
def blosum62 (a b : AA) : Int :=
  (matrix.getD (aaIndex a) []).getD (aaIndex b) 0

-- ============================================================================
-- Matrix well-formedness constants
-- ============================================================================

/-- Smallest BLOSUM62 diagonal (self-substitution) value.

    The diagonal entries are: A=4, R=5, N=6, D=6, C=9, Q=5, E=5, G=6, H=8,
    I=4, L=4, K=5, M=5, F=6, P=7, S=4, T=5, W=11, Y=7, V=4.
    Minimum is 4 (Ala, Ile, Leu, Ser, Val). -/
def minDiagScore : Int := 4

/-- Largest BLOSUM62 diagonal (self-substitution) value.

    Maximum is 11 (Trp). -/
def maxDiagScore : Int := 11

/-- Smallest BLOSUM62 entry overall (most penalized substitution).

    Several off-diagonal entries are -4 (e.g., W-N, W-D, W-C, A-W). -/
def minScore : Int := -4

/-- Largest BLOSUM62 entry overall (Trp self-substitution). -/
def maxScore : Int := 11

-- ============================================================================
-- Theorems: matrix is well-formed
-- ============================================================================

/-- Helper: the matrix has 20 rows. -/
theorem matrix_length : matrix.length = 20 := by decide

/-- Helper: every row has 20 columns. -/
theorem matrix_row_length : ∀ (i : Nat), i < 20 → (matrix.getD i []).length = 20 := by decide

/-- THEOREM (BLOSUM62 Symmetry): the BLOSUM62 matrix is symmetric.

    `blosum62 a b = blosum62 b a` for all amino acid pairs.

    This is a mathematical fact of the standard BLOSUM62 matrix and is
    proved by case analysis (decidable computation) over all 400 ordered
    amino-acid pairs. -/
theorem blosum62_symm (a b : AA) : blosum62 a b = blosum62 b a := by
  cases a <;> cases b <;> rfl

/-- THEOREM (BLOSUM62 diagonal non-negativity): every diagonal entry is
    non-negative.

    `blosum62 a a ≥ 0` for all amino acids `a`. -/
theorem blosum62_diag_nonneg (a : AA) : blosum62 a a ≥ 0 := by
  cases a <;> decide

/-- THEOREM (BLOSUM62 diagonal positivity): every diagonal entry is
    strictly positive (≥ 4). -/
theorem blosum62_diag_pos (a : AA) : blosum62 a a > 0 := by
  cases a <;> decide

/-- THEOREM (BLOSUM62 diagonal lower bound): every diagonal entry is
    ≥ 4 (the minimum diagonal value, attained by Ala/Ile/Leu/Ser/Val). -/
theorem blosum62_diag_min (a : AA) : blosum62 a a ≥ minDiagScore := by
  cases a <;> decide

/-- THEOREM (BLOSUM62 diagonal upper bound): every diagonal entry is
    ≤ 11 (the maximum diagonal value, attained by Trp). -/
theorem blosum62_diag_max (a : AA) : blosum62 a a ≤ maxDiagScore := by
  cases a <;> decide

/-- THEOREM (BLOSUM62 minimum entry): every entry is ≥ -4. -/
theorem blosum62_min (a b : AA) : blosum62 a b ≥ minScore := by
  cases a <;> cases b <;> decide

/-- THEOREM (BLOSUM62 maximum entry): every entry is ≤ 11. -/
theorem blosum62_max (a b : AA) : blosum62 a b ≤ maxScore := by
  cases a <;> cases b <;> decide

/-- THEOREM (BLOSUM62 self-substitution non-negativity): the
    self-substitution score `blosum62 a a` is non-negative for every
    amino acid `a`.  This is a direct corollary of `blosum62_diag_nonneg`. -/
theorem blosum62_diag_self_sub_nonneg (a : AA) : blosum62 a a ≥ 0 :=
  blosum62_diag_nonneg a

/-- The proposition that the BLOSUM62 matrix is well-formed: it is
    symmetric and every diagonal entry is non-negative.  This is the
    necessary condition used by `slotPropertySemantics` for the
    `ConservationScore` SLOT predicate. -/
def wellFormed : Prop :=
  (∀ (a b : AA), blosum62 a b = blosum62 b a) ∧
  (∀ (a : AA), blosum62 a a ≥ 0)

/-- THEOREM (BLOSUM62 well-formedness): the matrix is symmetric and has
    non-negative diagonal entries. -/
theorem wellFormed_proof : wellFormed := by
  exact ⟨fun a b => blosum62_symm a b, fun a => blosum62_diag_nonneg a⟩

end BLOSUM62

end BioCompiler
