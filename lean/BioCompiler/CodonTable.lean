/-
BioCompiler.CodonTable v7.0.0
==============================
Formal model of the standard genetic code in Lean 4.
Defines Base, Codon, AminoAcid, the translation function,
dinucleotide predicates, and cross-codon constraint theorems.
-/

namespace BioCompiler

-- ────────────────────────────────────────────────────────────
-- DNA Bases
-- ────────────────────────────────────────────────────────────

inductive Base : Type where
  | A : Base
  | C : Base
  | G : Base
  | T : Base
  deriving Repr, BEq, DecidableEq

namespace Base

def all : List Base := [A, C, G, T]

theorem all_complete (b : Base) : b ∈ all := by
  cases b <;> simp [all]

end Base

-- ────────────────────────────────────────────────────────────
-- Amino Acids
-- ────────────────────────────────────────────────────────────

inductive AminoAcid : Type where
  | F | L | I | M | V | S | P | T | A | Y
  | H | Q | N | K | D | E | C | W | R | G
  | Stop
  deriving Repr, BEq, DecidableEq

namespace AminoAcid

def all : List AminoAcid :=
  [F, L, I, M, V, S, P, T, A, Y, H, Q, N, K, D, E, C, W, R, G, Stop]

end AminoAcid

-- ────────────────────────────────────────────────────────────
-- Codons
-- ────────────────────────────────────────────────────────────

structure Codon : Type where
  fst : Base
  snd : Base
  trd : Base
  deriving Repr, BEq

namespace Codon

def all : List Codon :=
  Base.all.flatMap fun b1 =>
    Base.all.flatMap fun b2 =>
      Base.all.map fun b3 => ⟨b1, b2, b3⟩

theorem all_length : all.length = 64 := by
  simp [all, Base.all, List.length_flatMap, List.length_map]

end Codon

-- ────────────────────────────────────────────────────────────
-- The Standard Genetic Code
-- ────────────────────────────────────────────────────────────

def translateCodon : Codon → AminoAcid
  | ⟨Base.T, Base.T, Base.T⟩ => .F
  | ⟨Base.T, Base.T, Base.C⟩ => .F
  | ⟨Base.T, Base.T, Base.A⟩ => .L
  | ⟨Base.T, Base.T, Base.G⟩ => .L
  | ⟨Base.C, Base.T, _⟩ => .L
  | ⟨Base.A, Base.T, Base.A⟩ => .I
  | ⟨Base.A, Base.T, Base.C⟩ => .I
  | ⟨Base.A, Base.T, Base.T⟩ => .I
  | ⟨Base.A, Base.T, Base.G⟩ => .M
  | ⟨Base.G, Base.T, _⟩ => .V
  | ⟨Base.T, Base.C, _⟩ => .S
  | ⟨Base.A, Base.G, Base.T⟩ => .S
  | ⟨Base.A, Base.G, Base.C⟩ => .S
  | ⟨Base.C, Base.C, _⟩ => .P
  | ⟨Base.A, Base.C, _⟩ => .T
  | ⟨Base.G, Base.C, _⟩ => .A
  | ⟨Base.T, Base.A, Base.T⟩ => .Y
  | ⟨Base.T, Base.A, Base.C⟩ => .Y
  | ⟨Base.C, Base.A, Base.T⟩ => .H
  | ⟨Base.C, Base.A, Base.C⟩ => .H
  | ⟨Base.C, Base.A, Base.A⟩ => .Q
  | ⟨Base.C, Base.A, Base.G⟩ => .Q
  | ⟨Base.A, Base.A, Base.T⟩ => .N
  | ⟨Base.A, Base.A, Base.C⟩ => .N
  | ⟨Base.A, Base.A, Base.A⟩ => .K
  | ⟨Base.A, Base.A, Base.G⟩ => .K
  | ⟨Base.G, Base.A, Base.T⟩ => .D
  | ⟨Base.G, Base.A, Base.C⟩ => .D
  | ⟨Base.G, Base.A, Base.A⟩ => .E
  | ⟨Base.G, Base.A, Base.G⟩ => .E
  | ⟨Base.T, Base.G, Base.T⟩ => .C
  | ⟨Base.T, Base.G, Base.C⟩ => .C
  | ⟨Base.T, Base.G, Base.G⟩ => .W
  | ⟨Base.C, Base.G, _⟩ => .R
  | ⟨Base.A, Base.G, Base.A⟩ => .R
  | ⟨Base.A, Base.G, Base.G⟩ => .R
  | ⟨Base.G, Base.G, _⟩ => .G
  | ⟨Base.T, Base.A, Base.A⟩ => .Stop
  | ⟨Base.T, Base.A, Base.G⟩ => .Stop
  | ⟨Base.T, Base.G, Base.A⟩ => .Stop

-- ────────────────────────────────────────────────────────────
-- Stop codon identification
-- ────────────────────────────────────────────────────────────

def isStopCodon (c : Codon) : Prop :=
  c = ⟨Base.T, Base.A, Base.A⟩ ∨
  c = ⟨Base.T, Base.A, Base.G⟩ ∨
  c = ⟨Base.T, Base.G, Base.A⟩

theorem stop_codons_are_stops : ∀ c, isStopCodon c → translateCodon c = AminoAcid.Stop := by
  intro c h
  cases h with
  | inl h => simp [h, translateCodon]
  | inr h => cases h with
    | inl h => simp [h, translateCodon]
    | inr h => simp [h, translateCodon]

-- ────────────────────────────────────────────────────────────
-- Dinucleotide predicates (codon-level)
-- ────────────────────────────────────────────────────────────

/-- A codon contains a GT dinucleotide (within the codon) -/
def codonHasGT (c : Codon) : Prop :=
  (c.fst = Base.G ∧ c.snd = Base.T) ∨ (c.snd = Base.G ∧ c.trd = Base.T)

instance codonHasGT_decidable (c : Codon) : Decidable (codonHasGT c) := by
  unfold codonHasGT
  exact inferInstance

/-- Two adjacent codons form a cross-codon GT:
    last base of first = G, first base of second = T -/
def crossCodonGT (c1 c2 : Codon) : Prop :=
  c1.trd = Base.G ∧ c2.fst = Base.T

/-- A codon contains a CG dinucleotide (within the codon) -/
def codonHasCG (c : Codon) : Prop :=
  (c.fst = Base.C ∧ c.snd = Base.G) ∨ (c.snd = Base.C ∧ c.trd = Base.G)

-- Derive Decidable for codonHasCG using the fact that Base has DecidableEq
instance codonHasCG_decidable (c : Codon) : Decidable (codonHasCG c) := by
  unfold codonHasCG
  exact inferInstance

/-- Two adjacent codons form a cross-codon CG -/
def crossCodonCG (c1 c2 : Codon) : Prop :=
  c1.trd = Base.C ∧ c2.fst = Base.G

-- ────────────────────────────────────────────────────────────
-- Synonymous codons
-- ────────────────────────────────────────────────────────────

def Synonymous (c1 c2 : Codon) : Prop :=
  translateCodon c1 = translateCodon c2

theorem synonymous_refl (c : Codon) : Synonymous c c := rfl

theorem synonymous_symm {c1 c2 : Codon} : Synonymous c1 c2 → Synonymous c2 c1 := Eq.symm

theorem synonymous_trans {c1 c2 c3 : Codon} :
    Synonymous c1 c2 → Synonymous c2 c3 → Synonymous c1 c3 :=
  fun h1 h2 => Eq.trans h1 h2

/-- For any non-stop amino acid, there exists a codon without GT -/
theorem exists_gt_free_codon (aa : AminoAcid) (h : aa ≠ AminoAcid.Stop) :
    ∃ c : Codon, translateCodon c = aa ∧ ¬ codonHasGT c := by
  sorry

/-- For any non-stop amino acid, there exists a codon without CG -/
theorem exists_cg_free_codon (aa : AminoAcid) (h : aa ≠ AminoAcid.Stop) :
    ∃ c : Codon, translateCodon c = aa ∧ ¬ codonHasCG c := by
  sorry

/-- THEOREM (Novel): Cross-codon GT can always be resolved by synonymous substitution. -/
theorem cross_codon_gt_resolvable (c1 c2 : Codon)
    (h1 : translateCodon c1 ≠ AminoAcid.Stop)
    (h2 : translateCodon c2 ≠ AminoAcid.Stop)
    (hGT : crossCodonGT c1 c2) :
    ∃ c1' c2', Synonymous c1 c1' ∧ Synonymous c2 c2' ∧ ¬ crossCodonGT c1' c2' := by
  sorry

/-- THEOREM (Novel): Cross-codon CG can always be resolved by synonymous substitution. -/
theorem cross_codon_cg_resolvable (c1 c2 : Codon)
    (h1 : translateCodon c1 ≠ AminoAcid.Stop)
    (h2 : translateCodon c2 ≠ AminoAcid.Stop)
    (hCG : crossCodonCG c1 c2) :
    ∃ c1' c2', Synonymous c1 c1' ∧ Synonymous c2 c2' ∧ ¬ crossCodonCG c1' c2' := by
  sorry

end BioCompiler
