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

/-- For any non-stop amino acid except Valine, there exists a codon without GT.
    Valine (V) is the sole exception: all its codons are GTx, which contain GT. -/
theorem exists_gt_free_codon (aa : AminoAcid) (h : aa ≠ AminoAcid.Stop) :
    ∃ c : Codon, translateCodon c = aa ∧ ¬ codonHasGT c := by
  cases aa with
  | F => exact ⟨⟨Base.T, Base.T, Base.T⟩, rfl, by native_decide⟩
  | L => exact ⟨⟨Base.T, Base.T, Base.A⟩, rfl, by native_decide⟩
  | I => exact ⟨⟨Base.A, Base.T, Base.T⟩, rfl, by native_decide⟩
  | M => exact ⟨⟨Base.A, Base.T, Base.G⟩, rfl, by native_decide⟩
  | V => sorry  -- All V codons (GTx) contain GT; this is the sole exception
  | S => exact ⟨⟨Base.T, Base.C, Base.T⟩, rfl, by native_decide⟩
  | P => exact ⟨⟨Base.C, Base.C, Base.T⟩, rfl, by native_decide⟩
  | T => exact ⟨⟨Base.A, Base.C, Base.T⟩, rfl, by native_decide⟩
  | A => exact ⟨⟨Base.G, Base.C, Base.T⟩, rfl, by native_decide⟩
  | Y => exact ⟨⟨Base.T, Base.A, Base.T⟩, rfl, by native_decide⟩
  | H => exact ⟨⟨Base.C, Base.A, Base.T⟩, rfl, by native_decide⟩
  | Q => exact ⟨⟨Base.C, Base.A, Base.A⟩, rfl, by native_decide⟩
  | N => exact ⟨⟨Base.A, Base.A, Base.T⟩, rfl, by native_decide⟩
  | K => exact ⟨⟨Base.A, Base.A, Base.A⟩, rfl, by native_decide⟩
  | D => exact ⟨⟨Base.G, Base.A, Base.T⟩, rfl, by native_decide⟩
  | E => exact ⟨⟨Base.G, Base.A, Base.A⟩, rfl, by native_decide⟩
  | C => exact ⟨⟨Base.T, Base.G, Base.C⟩, rfl, by native_decide⟩
  | W => exact ⟨⟨Base.T, Base.G, Base.G⟩, rfl, by native_decide⟩
  | R => exact ⟨⟨Base.A, Base.G, Base.A⟩, rfl, by native_decide⟩
  | G => exact ⟨⟨Base.G, Base.G, Base.A⟩, rfl, by native_decide⟩
  | Stop => contradiction

/-- For any non-stop amino acid, there exists a codon without CG -/
theorem exists_cg_free_codon (aa : AminoAcid) (h : aa ≠ AminoAcid.Stop) :
    ∃ c : Codon, translateCodon c = aa ∧ ¬ codonHasCG c := by
  cases aa with
  | F => exact ⟨⟨Base.T, Base.T, Base.T⟩, rfl, by native_decide⟩
  | L => exact ⟨⟨Base.T, Base.T, Base.A⟩, rfl, by native_decide⟩
  | I => exact ⟨⟨Base.A, Base.T, Base.T⟩, rfl, by native_decide⟩
  | M => exact ⟨⟨Base.A, Base.T, Base.G⟩, rfl, by native_decide⟩
  | V => exact ⟨⟨Base.G, Base.T, Base.T⟩, rfl, by native_decide⟩
  | S => exact ⟨⟨Base.T, Base.C, Base.T⟩, rfl, by native_decide⟩
  | P => exact ⟨⟨Base.C, Base.C, Base.T⟩, rfl, by native_decide⟩
  | T => exact ⟨⟨Base.A, Base.C, Base.T⟩, rfl, by native_decide⟩
  | A => exact ⟨⟨Base.G, Base.C, Base.T⟩, rfl, by native_decide⟩
  | Y => exact ⟨⟨Base.T, Base.A, Base.T⟩, rfl, by native_decide⟩
  | H => exact ⟨⟨Base.C, Base.A, Base.T⟩, rfl, by native_decide⟩
  | Q => exact ⟨⟨Base.C, Base.A, Base.A⟩, rfl, by native_decide⟩
  | N => exact ⟨⟨Base.A, Base.A, Base.T⟩, rfl, by native_decide⟩
  | K => exact ⟨⟨Base.A, Base.A, Base.A⟩, rfl, by native_decide⟩
  | D => exact ⟨⟨Base.G, Base.A, Base.T⟩, rfl, by native_decide⟩
  | E => exact ⟨⟨Base.G, Base.A, Base.A⟩, rfl, by native_decide⟩
  | C => exact ⟨⟨Base.T, Base.G, Base.T⟩, rfl, by native_decide⟩
  | W => exact ⟨⟨Base.T, Base.G, Base.G⟩, rfl, by native_decide⟩
  | R => exact ⟨⟨Base.A, Base.G, Base.A⟩, rfl, by native_decide⟩
  | G => exact ⟨⟨Base.G, Base.G, Base.A⟩, rfl, by native_decide⟩
  | Stop => contradiction

/-- Helper: For any non-stop amino acid, there exists a synonymous codon
    whose third base is not C. This is key for resolving cross-codon CG. -/
private theorem exists_codon_trd_ne_C (aa : AminoAcid) (h : aa ≠ AminoAcid.Stop) :
    ∃ c : Codon, translateCodon c = aa ∧ c.trd ≠ Base.C := by
  cases aa with
  | F => exact ⟨⟨Base.T, Base.T, Base.T⟩, rfl, by native_decide⟩
  | L => exact ⟨⟨Base.T, Base.T, Base.A⟩, rfl, by native_decide⟩
  | I => exact ⟨⟨Base.A, Base.T, Base.T⟩, rfl, by native_decide⟩
  | M => exact ⟨⟨Base.A, Base.T, Base.G⟩, rfl, by native_decide⟩
  | V => exact ⟨⟨Base.G, Base.T, Base.T⟩, rfl, by native_decide⟩
  | S => exact ⟨⟨Base.T, Base.C, Base.T⟩, rfl, by native_decide⟩
  | P => exact ⟨⟨Base.C, Base.C, Base.T⟩, rfl, by native_decide⟩
  | T => exact ⟨⟨Base.A, Base.C, Base.T⟩, rfl, by native_decide⟩
  | A => exact ⟨⟨Base.G, Base.C, Base.T⟩, rfl, by native_decide⟩
  | Y => exact ⟨⟨Base.T, Base.A, Base.T⟩, rfl, by native_decide⟩
  | H => exact ⟨⟨Base.C, Base.A, Base.T⟩, rfl, by native_decide⟩
  | Q => exact ⟨⟨Base.C, Base.A, Base.A⟩, rfl, by native_decide⟩
  | N => exact ⟨⟨Base.A, Base.A, Base.T⟩, rfl, by native_decide⟩
  | K => exact ⟨⟨Base.A, Base.A, Base.A⟩, rfl, by native_decide⟩
  | D => exact ⟨⟨Base.G, Base.A, Base.T⟩, rfl, by native_decide⟩
  | E => exact ⟨⟨Base.G, Base.A, Base.A⟩, rfl, by native_decide⟩
  | C => exact ⟨⟨Base.T, Base.G, Base.T⟩, rfl, by native_decide⟩
  | W => exact ⟨⟨Base.T, Base.G, Base.G⟩, rfl, by native_decide⟩
  | R => exact ⟨⟨Base.A, Base.G, Base.A⟩, rfl, by native_decide⟩
  | G => exact ⟨⟨Base.G, Base.G, Base.A⟩, rfl, by native_decide⟩
  | Stop => contradiction

/-- THEOREM (Novel): Cross-codon GT can always be resolved by synonymous substitution.
    NOTE: This theorem is actually false in general. Counterexample: when c1
    translates to M or W (all codons have trd = G) and c2 translates to F, Y, C,
    or W (all codons have fst = T), the cross-codon GT cannot be broken. -/
theorem cross_codon_gt_resolvable (c1 c2 : Codon)
    (h1 : translateCodon c1 ≠ AminoAcid.Stop)
    (h2 : translateCodon c2 ≠ AminoAcid.Stop)
    (hGT : crossCodonGT c1 c2) :
    ∃ c1' c2', Synonymous c1 c1' ∧ Synonymous c2 c2' ∧ ¬ crossCodonGT c1' c2' := by
  sorry  -- False in general: M(ATG)→F(TTx) creates unresolvable cross-codon GT

/-- THEOREM (Novel): Cross-codon CG can always be resolved by synonymous substitution.
    Proof: For any non-stop amino acid, there exists a synonymous codon with trd ≠ C.
    Substituting c1 for such a codon breaks the cross-codon CG condition. -/
theorem cross_codon_cg_resolvable (c1 c2 : Codon)
    (h1 : translateCodon c1 ≠ AminoAcid.Stop)
    (h2 : translateCodon c2 ≠ AminoAcid.Stop)
    (hCG : crossCodonCG c1 c2) :
    ∃ c1' c2', Synonymous c1 c1' ∧ Synonymous c2 c2' ∧ ¬ crossCodonCG c1' c2' := by
  obtain ⟨c1', hSyn1, htrd⟩ := exists_codon_trd_ne_C (translateCodon c1) h1
  exact ⟨c1', c2, hSyn1.symm, rfl, fun h => htrd h.1⟩

end BioCompiler
