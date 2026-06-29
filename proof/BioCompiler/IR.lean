/-
  BioCompiler.IR — Intermediate Representation Types and Lowering Correctness

  This module defines the five IR levels (IR-L0 through IR-L4) and proves
  that the lowering passes (transcription, splicing, translation) preserve
  the semantic content of the gene design.

  The central theorem is:

    compile_correctness:
      For any IR-L0 (genomic DNA), compiling to IR-L3 (polypeptide)
      via transcribe → splice → translate produces the correct protein,
      i.e., the amino acid sequence is the standard genetic code translation
      of the coding sequence in the genomic DNA.

  This is the "compiler correctness" theorem — analogous to CompCert's
  compilation preserves observable behavior.

  NOTE on proof status:
    All correctness theorems in this module are proved fully — no `sorry`
    remains. `codonToAA_stop_iff` is a finite case analysis over the 64
    codons (closed by `cases <;> decide`), and `translate_correctness`
    follows from it together with the helper `groupCodons_last?_spec`
    (which identifies the last codon of `groupCodons` with the last 3
    bases of the CDS). transcribe_correctness, splice_correctness, and
    compile_correctness (composition) are proved by reduction to the
    definitions, since `correctTranslation` is *defined* as the
    composition transcribe → extract exons → translate.

  Reference: DOC-03 (SDD) §3.2 (IR Pipeline), DOC-04 (ICD) §3
-/

import BioCompiler.Sequence
import BioCompiler.TypeSystem
import BioCompiler.ThreeValued

namespace BioCompiler

-- ==============================================================================
-- IR Level Enum
-- ==============================================================================

/--
The five IR levels of the BioCompiler pipeline.
  - L0: genomic DNA
  - L1: pre-mRNA (post-transcription, pre-splicing)
  - L2: mature mRNA (post-splicing, with explicit UTR/CDS structure)
  - L3: polypeptide (post-translation)
  - L4: folded protein (structure prediction delegated to SLOT oracles)
-/
inductive IRLevel where
  | L0 : IRLevel
  | L1 : IRLevel
  | L2 : IRLevel
  | L3 : IRLevel
  | L4 : IRLevel
  deriving Repr, DecidableEq

-- ==============================================================================
-- Gene Region (annotations on the genomic sequence)
-- ==============================================================================

/--
Type of an annotated region of a gene. Used by IR-L0 and IR-L1 to record
the biological role of each subsequence (exon, intron, UTR, CDS, regulatory).
-/
inductive RegionType where
  | exon : RegionType
  | intron : RegionType
  | fiveUTR : RegionType
  | threeUTR : RegionType
  | cds : RegionType
  | promoter : RegionType
  | terminator : RegionType
  deriving Repr, DecidableEq

/--
An annotated region: a half-open interval `[start, end_pos)` on the
sequence, together with its region type.
-/
structure GeneRegion where
  start : Nat
  end_pos : Nat
  region_type : RegionType
  deriving Repr

-- ==============================================================================
-- RNA Bases
--
-- The existing `Nucleotide` type (in `BioCompiler.Sequence`) models DNA bases
-- (A, C, G, T). RNA uses uracil (U) in place of thymine (T), so we define a
-- distinct `RNABase` inductive. This also enforces, at the type level, that
-- IR-L0 (DNA) and IR-L1/L2 (RNA) cannot be confused.
-- ==============================================================================

inductive RNABase where
  | A : RNABase
  | C : RNABase
  | G : RNABase
  | U : RNABase
  deriving DecidableEq, Repr, BEq, Inhabited

abbrev RNASequence := List RNABase

-- ==============================================================================
-- IR-L0: Genomic DNA
-- ==============================================================================

/--
IR-L0: the input to the compiler. A genomic DNA sequence together with
region annotations, the source organism, and an optional gene name.
-/
structure IR_L0 where
  sequence : Sequence           -- DNA sequence (List Nucleotide)
  regions : List GeneRegion
  organism : String
  gene_name : Option String
  deriving Repr

-- ==============================================================================
-- IR-L1: Pre-mRNA (post-transcription)
-- ==============================================================================

/--
IR-L1: the result of transcription. The DNA sequence has been converted
to RNA (T → U); region annotations are inherited from IR-L0.
-/
structure IR_L1 where
  sequence : RNASequence
  regions : List GeneRegion
  organism : String
  gene_name : Option String
  deriving Repr

-- ==============================================================================
-- IR-L2: Mature mRNA (post-splicing)
-- ==============================================================================

/--
IR-L2: the result of splicing. Introns have been removed; the mature mRNA
is split into a 5' UTR, a coding sequence (CDS), and a 3' UTR.
The CDS is required to start with AUG, end with a stop codon, and have
length divisible by 3 (these are well-formedness conditions, see
`wellFormedL2`).
-/
structure IR_L2 where
  five_utr : RNASequence
  cds : RNASequence
  three_utr : RNASequence
  organism : String
  gene_name : Option String
  deriving Repr

-- ==============================================================================
-- IR-L3: Polypeptide (post-translation)
-- ==============================================================================

/--
The 21 standard amino-acid values (20 amino acids + Stop). The order of
constructors is the standard IUPAC single-letter code alphabetical order.
-/
inductive AminoAcid where
  | Ala | Arg | Asn | Asp | Cys
  | Gln | Glu | Gly | His | Ile
  | Leu | Lys | Met | Phe | Pro
  | Ser | Thr | Trp | Tyr | Val
  | Stop  -- stop codon (terminates translation)
  deriving Repr, DecidableEq, BEq

/--
IR-L3: the result of translation. A list of amino acids (the polypeptide
chain), ending in `AminoAcid.Stop`.
-/
structure IR_L3 where
  sequence : List AminoAcid
  organism : String
  gene_name : Option String
  deriving Repr

-- ==============================================================================
-- IR-L4: Folded protein (placeholder)
--
-- The 3D structure of the protein is not modeled concretely in the formal IR;
-- it is delegated to SLOT-dependent oracles (see `OracleProofs.lean`:
-- StructureConfidence, CorrectFoldTopology, etc.). IR-L4 therefore carries
-- only the polypeptide sequence plus the organism/name metadata.
-- ==============================================================================

structure IR_L4 where
  sequence : List AminoAcid
  organism : String
  gene_name : Option String
  deriving Repr

-- ==============================================================================
-- Transcription: L0 → L1
-- ==============================================================================

/--
Transcription of a single DNA base to its RNA equivalent.
T → U; A, C, G are unchanged.
-/
def transcribeBase (b : Nucleotide) : RNABase :=
  match b with
  | Nucleotide.T => RNABase.U
  | Nucleotide.A => RNABase.A
  | Nucleotide.C => RNABase.C
  | Nucleotide.G => RNABase.G

/--
Transcription of an entire IR-L0 to IR-L1: map `transcribeBase` over the
DNA sequence, preserving region annotations and metadata.
-/
def transcribe (ir : IR_L0) : IR_L1 :=
  { sequence := ir.sequence.map transcribeBase
  , regions := ir.regions
  , organism := ir.organism
  , gene_name := ir.gene_name }

/--
Transcription correctness: the RNA sequence is the DNA sequence with each
base mapped via `transcribeBase` (T → U, A → A, C → C, G → G).
This is the L0 → L1 lowering correctness theorem.
-/
theorem transcribe_correctness (ir : IR_L0) :
  (transcribe ir).sequence = ir.sequence.map transcribeBase := by
  rfl

/--
Transcription preserves metadata and region annotations.
-/
theorem transcribe_preserves_metadata (ir : IR_L0) :
  (transcribe ir).regions = ir.regions ∧
  (transcribe ir).organism = ir.organism ∧
  (transcribe ir).gene_name = ir.gene_name := by
  simp [transcribe]

/--
Transcription preserves sequence length (the T → U substitution is bijective).
-/
theorem transcribe_preserves_length (ir : IR_L0) :
  (transcribe ir).sequence.length = ir.sequence.length := by
  rw [transcribe_correctness, List.length_map]

-- ==============================================================================
-- Splicing: L1 → L2
-- ==============================================================================

/--
An exon range: (start, length) of an exon within the pre-mRNA sequence.
-/
abbrev ExonRange := (Nat × Nat)

/--
Extract a subsequence from an RNA sequence given a (start, length) pair.
Out-of-range ranges yield the empty sequence.
-/
def extractSubseq (seq : RNASequence) (r : ExonRange) : RNASequence :=
  let (s, l) := r
  if _ : s + l ≤ seq.length then
    (seq.drop s).take l
  else
    []

/--
Splicing: L1 → L2.

Given a pre-mRNA and a list of exon ranges (in genomic order), produce
mature mRNA by concatenating the exon sub-sequences into the CDS.
The 5'/3' UTRs are modeled separately; here they are taken as inputs to
the splice function (typically extracted from `fiveUTR`/`threeUTR` regions).

The key correctness property is that introns are removed and exons are
concatenated in order — this is what `splice_correctness` states.
-/
def splice (ir : IR_L1) (exonRanges : List ExonRange)
    (five_utr three_utr : RNASequence) : IR_L2 :=
  { five_utr := five_utr
  , cds := exonRanges.flatMap (extractSubseq ir.sequence)
  , three_utr := three_utr
  , organism := ir.organism
  , gene_name := ir.gene_name }

/--
Splicing correctness: the CDS of the spliced mRNA is the concatenation
of the exon sub-sequences extracted from the pre-mRNA.
This is the L1 → L2 lowering correctness theorem.
-/
theorem splice_correctness (ir : IR_L1) (exonRanges : List ExonRange)
    (five_utr three_utr : RNASequence) :
  (splice ir exonRanges five_utr three_utr).cds =
    exonRanges.flatMap (extractSubseq ir.sequence) := by
  rfl

/--
Splicing preserves metadata.
-/
theorem splice_preserves_metadata (ir : IR_L1) (exonRanges : List ExonRange)
    (five_utr three_utr : RNASequence) :
  (splice ir exonRanges five_utr three_utr).organism = ir.organism ∧
  (splice ir exonRanges five_utr three_utr).gene_name = ir.gene_name := by
  simp [splice]

-- ==============================================================================
-- Translation: L2 → L3
-- ==============================================================================

/--
The standard genetic code: maps an RNA codon (three bases) to an amino acid.

This is the canonical NCBI translation table 1 (standard code).
The match is exhaustive over the 4^3 = 64 codons; the trailing wildcard
covers no reachable case but is included for totality robustness.
-/
def codonToAA (c1 c2 c3 : RNABase) : AminoAcid :=
  match c1, c2, c3 with
  -- U_U_
  | RNABase.U, RNABase.U, RNABase.U => AminoAcid.Phe
  | RNABase.U, RNABase.U, RNABase.C => AminoAcid.Phe
  | RNABase.U, RNABase.U, RNABase.A => AminoAcid.Leu
  | RNABase.U, RNABase.U, RNABase.G => AminoAcid.Leu
  -- U_C_ (Ser)
  | RNABase.U, RNABase.C, _ => AminoAcid.Ser
  -- U_A_
  | RNABase.U, RNABase.A, RNABase.U => AminoAcid.Tyr
  | RNABase.U, RNABase.A, RNABase.C => AminoAcid.Tyr
  | RNABase.U, RNABase.A, RNABase.A => AminoAcid.Stop
  | RNABase.U, RNABase.A, RNABase.G => AminoAcid.Stop
  -- U_G_
  | RNABase.U, RNABase.G, RNABase.U => AminoAcid.Cys
  | RNABase.U, RNABase.G, RNABase.C => AminoAcid.Cys
  | RNABase.U, RNABase.G, RNABase.A => AminoAcid.Stop
  | RNABase.U, RNABase.G, RNABase.G => AminoAcid.Trp
  -- C_U_ (Leu)
  | RNABase.C, RNABase.U, RNABase.U => AminoAcid.Leu
  | RNABase.C, RNABase.U, RNABase.C => AminoAcid.Leu
  | RNABase.C, RNABase.U, RNABase.A => AminoAcid.Leu
  | RNABase.C, RNABase.U, RNABase.G => AminoAcid.Leu
  -- C_C_ (Pro)
  | RNABase.C, RNABase.C, _ => AminoAcid.Pro
  -- C_A_
  | RNABase.C, RNABase.A, RNABase.U => AminoAcid.His
  | RNABase.C, RNABase.A, RNABase.C => AminoAcid.His
  | RNABase.C, RNABase.A, RNABase.A => AminoAcid.Gln
  | RNABase.C, RNABase.A, RNABase.G => AminoAcid.Gln
  -- C_G_ (Arg)
  | RNABase.C, RNABase.G, _ => AminoAcid.Arg
  -- A_U_
  | RNABase.A, RNABase.U, RNABase.U => AminoAcid.Ile
  | RNABase.A, RNABase.U, RNABase.C => AminoAcid.Ile
  | RNABase.A, RNABase.U, RNABase.A => AminoAcid.Ile
  | RNABase.A, RNABase.U, RNABase.G => AminoAcid.Met  -- start codon
  -- A_C_ (Thr)
  | RNABase.A, RNABase.C, _ => AminoAcid.Thr
  -- A_A_
  | RNABase.A, RNABase.A, RNABase.U => AminoAcid.Asn
  | RNABase.A, RNABase.A, RNABase.C => AminoAcid.Asn
  | RNABase.A, RNABase.A, RNABase.A => AminoAcid.Lys
  | RNABase.A, RNABase.A, RNABase.G => AminoAcid.Lys
  -- A_G_
  | RNABase.A, RNABase.G, RNABase.U => AminoAcid.Ser
  | RNABase.A, RNABase.G, RNABase.C => AminoAcid.Ser
  | RNABase.A, RNABase.G, RNABase.A => AminoAcid.Arg
  | RNABase.A, RNABase.G, RNABase.G => AminoAcid.Arg
  -- G_U_ (Val)
  | RNABase.G, RNABase.U, RNABase.U => AminoAcid.Val
  | RNABase.G, RNABase.U, RNABase.C => AminoAcid.Val
  | RNABase.G, RNABase.U, RNABase.A => AminoAcid.Val
  | RNABase.G, RNABase.U, RNABase.G => AminoAcid.Val
  -- G_C_ (Ala)
  | RNABase.G, RNABase.C, _ => AminoAcid.Ala
  -- G_A_
  | RNABase.G, RNABase.A, RNABase.U => AminoAcid.Asp
  | RNABase.G, RNABase.A, RNABase.C => AminoAcid.Asp
  | RNABase.G, RNABase.A, RNABase.A => AminoAcid.Glu
  | RNABase.G, RNABase.A, RNABase.G => AminoAcid.Glu
  -- G_G_ (Gly)
  | RNABase.G, RNABase.G, _ => AminoAcid.Gly
  -- (No catch-all: the patterns above are exhaustive over the 4^3 = 64 codons.
  --  Each wildcard `_` covers the remaining constructors of its position.)

/--
Group an RNA sequence into codons (triplets). Trailing 1–2 bases (if the
length is not divisible by 3) are silently dropped.
-/
def groupCodons : RNASequence → List (RNABase × RNABase × RNABase)
  | [] => []
  | [_] => []
  | [_a, _b] => []
  | a :: b :: c :: rest => (a, b, c) :: groupCodons rest

/--
Translate a CDS (list of RNA bases) to a list of amino acids by grouping
into codons and applying the standard genetic code.
-/
def translateCDS (cds : RNASequence) : List AminoAcid :=
  (groupCodons cds).map (fun (a, b, c) => codonToAA a b c)

/--
Translation of an entire IR-L2 to IR-L3: translate only the CDS (UTRs are
untranslated by definition).
-/
def translate (ir : IR_L2) : IR_L3 :=
  { sequence := translateCDS ir.cds
  , organism := ir.organism
  , gene_name := ir.gene_name }

/--
Translation preserves metadata.
-/
theorem translate_preserves_metadata (ir : IR_L2) :
  (translate ir).organism = ir.organism ∧
  (translate ir).gene_name = ir.gene_name := by
  simp [translate]

/--
The number of codons produced by `groupCodons` is `seq.length / 3`.
Proved by structural recursion mirroring `groupCodons` itself.
-/
theorem groupCodons_length : ∀ (seq : RNASequence),
    (groupCodons seq).length = seq.length / 3
  | [] => by simp [groupCodons]
  | [_] => by simp [groupCodons]
  | [_a, _b] => by simp [groupCodons]
  | a :: b :: c :: rest => by
    have ih := groupCodons_length rest
    simp only [groupCodons, List.length_cons, ih]
    -- (a :: b :: c :: rest).length = 3 + rest.length, and (3 + n)/3 = 1 + n/3
    omega

/--
`groupCodons_length` re-stated in terms of `translateCDS`: the polypeptide
length equals `cds.length / 3`.
-/
theorem translateCDS_length (cds : RNASequence) :
  (translateCDS cds).length = cds.length / 3 := by
  simp [translateCDS, groupCodons_length, List.length_map]

/--
Check whether an RNA codon is a stop codon (UAA, UAG, UGA).
-/
def isStopCodonRNA (c1 c2 c3 : RNABase) : Bool :=
  match c1, c2, c3 with
  | RNABase.U, RNABase.A, RNABase.A => true
  | RNABase.U, RNABase.A, RNABase.G => true
  | RNABase.U, RNABase.G, RNABase.A => true
  | _, _, _ => false

/--
Check whether an RNA sequence ends with a stop codon (i.e. its last 3 bases
form a stop codon). Length must be ≥ 3 for a meaningful answer.
-/
def endsWithStopCodon (cds : RNASequence) : Bool :=
  if _ : cds.length ≥ 3 then
    let last3 := (cds.drop (cds.length - 3)).take 3
    match last3 with
    | [a, b, c] => isStopCodonRNA a b c
    | _ => false
  else false

/--
`codonToAA` agrees with `isStopCodonRNA`: a codon is a stop codon iff
`codonToAA` returns `AminoAcid.Stop`.

This is a finite case analysis over 64 codons. The proof is mechanical:
we split into 64 concrete cases via `cases` and close each by `decide`.
-/
theorem codonToAA_stop_iff (c1 c2 c3 : RNABase) :
  codonToAA c1 c2 c3 = AminoAcid.Stop ↔ isStopCodonRNA c1 c2 c3 = true := by
  -- 64-way case analysis: enumerate all (c1, c2, c3) ∈ RNABase^3 (4^3 = 64
  -- codons) and close each branch by computation. `decide` cannot be used
  -- directly because the theorem statement is a Π-type over free variables,
  -- so we first split into 64 concrete cases via `cases` and then decide
  -- each (now ground) goal.
  cases c1 <;> cases c2 <;> cases c3 <;> decide

/--
Helper: when `seq.length` is divisible by 3 and positive, the last element of
`groupCodons seq` is the last 3 bases of `seq` (as a triple), and
`(seq.drop (seq.length - 3)).take 3` is that same 3-element list.

Proved by structural recursion mirroring `groupCodons` itself. The key
invariant is that `seq.length % 3 = 0` is preserved when peeling off the
first 3 bases (one codon).
-/
private theorem groupCodons_last?_spec : ∀ (seq : RNASequence),
    seq.length % 3 = 0 → 0 < seq.length →
    ∃ a b c, (groupCodons seq).getLast? = some (a, b, c) ∧
             (seq.drop (seq.length - 3)).take 3 = [a, b, c]
  | [] => fun _ h => by
    simp only [List.length_nil] at h
    exact absurd h (by decide)
  | [_] => fun h_len _ => by
    simp only [List.length_cons, List.length_nil, Nat.zero_add] at h_len
    exact absurd h_len (by decide)
  | [_, _] => fun h_len _ => by
    simp only [List.length_cons, List.length_nil, Nat.zero_add, Nat.add_assoc,
               Nat.add_one, Nat.one_add] at h_len
    exact absurd h_len (by decide)
  | a :: b :: c :: rest => fun h_len h_pos => by
    have h_seq_len : (a :: b :: c :: rest).length = 3 + rest.length := by
      have : (a :: b :: c :: rest).length = rest.length + 3 := by simp [List.length]
      omega
    have h_rest_mod : rest.length % 3 = 0 := by omega
    by_cases h_rest_empty : rest = []
    · -- rest = [], so seq = [a, b, c]; groupCodons = [(a,b,c)]; drop 0; take 3 = [a,b,c]
      subst h_rest_empty
      refine ⟨a, b, c, rfl, rfl⟩
    · -- rest ≠ [], so rest.length ≥ 3 (since rest.length % 3 = 0 and rest.length > 0)
      have h_rest_pos : 0 < rest.length := by
        cases rest with
        | nil => exact absurd rfl h_rest_empty
        | cons x rest1 =>
          have : (x :: rest1).length = rest1.length + 1 := rfl
          omega
      obtain ⟨a', b', c', h_ih_last, h_ih_drop⟩ :=
        groupCodons_last?_spec rest h_rest_mod h_rest_pos
      refine ⟨a', b', c', ?_, ?_⟩
      · -- (groupCodons (a::b::c::rest)).getLast? = some (a', b', c')
        have h_gc : groupCodons (a :: b :: c :: rest) = (a, b, c) :: groupCodons rest := rfl
        have h_gc_rest_len : (groupCodons rest).length = rest.length / 3 := groupCodons_length rest
        have h_gc_rest_pos : 0 < (groupCodons rest).length := by omega
        rw [h_gc]
        cases h_rc : (groupCodons rest) with
        | nil =>
          rw [h_rc, List.length_nil] at h_gc_rest_pos
          exact absurd h_gc_rest_pos (by decide)
        | cons x xs =>
          rw [h_rc] at h_ih_last
          -- In both sub-cases (xs = [] or xs = y :: ys), the last? of
          -- ((a,b,c) :: x :: xs) is definitionally the last? of (x :: xs).
          cases xs with
          | nil => exact h_ih_last
          | cons y ys => exact h_ih_last
      · -- ((a::b::c::rest).drop ((a::b::c::rest).length - 3)).take 3 = [a', b', c']
        rw [h_seq_len, Nat.add_sub_cancel_left]
        -- Now: ((a::b::c::rest).drop rest.length).take 3 = [a', b', c']
        have h_drop : (a :: b :: c :: rest).drop rest.length = rest.drop (rest.length - 3) := by
          have h1 : (a :: b :: c :: rest).drop 3 = rest := rfl
          have h2 : rest.length = 3 + (rest.length - 3) := by omega
          calc (a :: b :: c :: rest).drop rest.length
              = (a :: b :: c :: rest).drop (3 + (rest.length - 3)) := by rw [← h2]
            _ = ((a :: b :: c :: rest).drop 3).drop (rest.length - 3) := by rw [← List.drop_drop]
            _ = rest.drop (rest.length - 3) := by rw [h1]
        rw [h_drop, h_ih_drop]

/--
Translation correctness (strong form): if the CDS starts with AUG, ends with
a stop codon, and has length divisible by 3, then `translateCDS` produces a
polypeptide whose length is `cds.length / 3` and whose last amino acid is
`Stop`.

The length part follows from `translateCDS_length`. The "last is Stop" part
requires `codonToAA_stop_iff` (a finite case analysis) and the helper
`groupCodons_last?_spec` (which identifies the last codon of `groupCodons`
with the last 3 bases of the CDS).
-/
theorem translate_correctness (cds : RNASequence)
    (h_start : cds.take 3 = [RNABase.A, RNABase.U, RNABase.G])
    (h_stop : endsWithStopCodon cds = true)
    (h_len : cds.length % 3 = 0) :
    (translateCDS cds).length = cds.length / 3 ∧
    (translateCDS cds).getLast? = some AminoAcid.Stop := by
  refine ⟨translateCDS_length cds, ?_⟩
  -- 1. cds.length ≥ 3 (from h_start: cds.take 3 has length 3, so min 3 len = 3)
  have h_len3 : 3 ≤ cds.length := by
    have h_start_len : (cds.take 3).length = 3 := by rw [h_start]; rfl
    have h_take_len : (cds.take 3).length = min 3 cds.length := by rw [List.length_take]
    rw [h_take_len] at h_start_len
    omega
  -- 2. Extract the last codon (a, b, c) from groupCodons_last?_spec
  have h_pos : 0 < cds.length := by omega
  obtain ⟨a, b, c, h_gc_last, h_last3⟩ := groupCodons_last?_spec cds h_len h_pos
  -- h_gc_last : (groupCodons cds).getLast? = some (a, b, c)
  -- h_last3   : (cds.drop (cds.length - 3)).take 3 = [a, b, c]
  -- 3. From h_stop, conclude isStopCodonRNA a b c = true.
  --    Unfolding endsWithStopCodon: the match on the last 3 bases reduces
  --    to isStopCodonRNA a b c (by h_last3).
  have h_stop' : isStopCodonRNA a b c = true := by
    have h_eq : endsWithStopCodon cds = isStopCodonRNA a b c := by
      rw [endsWithStopCodon, dif_pos h_len3]
      simp only [h_last3]
    rw [h_eq] at h_stop
    exact h_stop
  -- 4. codonToAA a b c = AminoAcid.Stop (by codonToAA_stop_iff)
  have h_codon : codonToAA a b c = AminoAcid.Stop :=
    (codonToAA_stop_iff a b c).mpr h_stop'
  -- 5. (translateCDS cds).getLast? = some AminoAcid.Stop
  --    By List.getLast?_map and h_gc_last, this reduces to
  --    some (codonToAA a b c) = some AminoAcid.Stop, which follows from h_codon.
  rw [translateCDS, List.getLast?_map, h_gc_last]
  simp [h_codon]

-- ==============================================================================
-- Well-formedness predicates
-- ==============================================================================

/--
Well-formedness for IR-L0: a minimal sanity check that all regions lie
within the sequence bounds. A full well-formedness predicate would also
require region non-overlap, ordered boundaries, presence of a CDS region,
etc.; for the purpose of the compile-correctness theorem we only need
boundary soundness.
-/
def wellFormedL0 (ir : IR_L0) : Prop :=
  ir.regions.all (fun r => r.start ≤ r.end_pos ∧ r.end_pos ≤ ir.sequence.length)

/--
Well-formedness for IR-L2: the CDS has length divisible by 3, starts with
AUG, and ends with a stop codon.
-/
def wellFormedL2 (ir : IR_L2) : Prop :=
  ir.cds.length % 3 = 0 ∧
  ir.cds.take 3 = [RNABase.A, RNABase.U, RNABase.G] ∧
  endsWithStopCodon ir.cds = true

-- ==============================================================================
-- Full Compilation Correctness
-- ==============================================================================

/--
The reference (specification) translation of an IR-L0 design: transcribe,
extract the exon sub-sequences, translate via the standard genetic code.

`correctTranslation` is *defined* as the composition of the lowering
passes' specification, so `compile_correctness` below holds by reduction.
-/
def correctTranslation (ir : IR_L0) (exonRanges : List ExonRange)
    (_five_utr _three_utr : RNASequence) : List AminoAcid :=
  -- The 5'/3' UTRs are not translated; only the CDS (exon concatenation) is.
  -- The UTR parameters are present for structural symmetry with IR_L2.
  let rna := ir.sequence.map transcribeBase               -- = (transcribe ir).sequence
  let cds := exonRanges.flatMap (extractSubseq rna)        -- = (splice _ _ _ _).cds
  translateCDS cds                                          -- = (translate _).sequence

/--
The central theorem: compiling IR-L0 to IR-L3 via the lowering pipeline
`transcribe → splice → translate` produces a polypeptide whose sequence is
exactly the specification `correctTranslation`.

Proof: by reduction. Each lowering pass is defined to produce exactly the
components of `correctTranslation`, so the equality holds definitionally.
No additional inductive hypotheses are needed beyond well-formedness (which
ensures the CDS is a valid translatable sequence — but the *equality*
itself does not depend on well-formedness; well-formedness only ensures
the result is biologically meaningful).
-/
theorem compile_correctness (ir_l0 : IR_L0) (exonRanges : List ExonRange)
    (five_utr three_utr : RNASequence) :
    (translate (splice (transcribe ir_l0) exonRanges five_utr three_utr)).sequence =
      correctTranslation ir_l0 exonRanges five_utr three_utr := by
  -- Unfold each lowering pass and the specification:
  show translateCDS (splice (transcribe ir_l0) exonRanges five_utr three_utr).cds =
       translateCDS (exonRanges.flatMap (extractSubseq (ir_l0.sequence.map transcribeBase)))
  -- `(splice _ _ _ _).cds` reduces to `exonRanges.flatMap (extractSubseq (transcribe _).sequence)`:
  rw [splice_correctness]
  -- `(transcribe ir_l0).sequence` reduces to `ir_l0.sequence.map transcribeBase`;
  -- the two sides are then syntactically identical (rw closes with rfl):
  rw [transcribe_correctness]

/--
Corollary: the compiled polypeptide's length is `cds.length / 3`, where
`cds` is the spliced CDS — i.e., `correctTranslation` produces one amino
acid per codon.
-/
theorem compile_correctness_length (ir_l0 : IR_L0) (exonRanges : List ExonRange)
    (five_utr three_utr : RNASequence) :
    (translate (splice (transcribe ir_l0) exonRanges five_utr three_utr)).sequence.length =
      (exonRanges.flatMap (extractSubseq (ir_l0.sequence.map transcribeBase))).length / 3 := by
  rw [compile_correctness]
  -- `correctTranslation` unfolds (by kernel reduction of let-bindings) to:
  show (translateCDS (exonRanges.flatMap (extractSubseq (ir_l0.sequence.map transcribeBase)))).length =
       (exonRanges.flatMap (extractSubseq (ir_l0.sequence.map transcribeBase))).length / 3
  rw [translateCDS_length]

-- ==============================================================================
-- IR Level L3 → L4 (folding)
-- ==============================================================================

/--
Folding (L3 → L4): in the formal IR, this is the identity on the polypeptide
sequence. The actual 3D structure prediction is delegated to SLOT oracles
(see `OracleProofs.lean`: StructureConfidence, CorrectFoldTopology, etc.),
which are soundness-vacuous by design.
-/
def fold (ir : IR_L3) : IR_L4 :=
  { sequence := ir.sequence
  , organism := ir.organism
  , gene_name := ir.gene_name }

/--
Folding preserves the polypeptide sequence.
-/
theorem fold_preserves_sequence (ir : IR_L3) :
  (fold ir).sequence = ir.sequence := by
  rfl

/--
Folding preserves metadata.
-/
theorem fold_preserves_metadata (ir : IR_L3) :
  (fold ir).organism = ir.organism ∧
  (fold ir).gene_name = ir.gene_name := by
  simp [fold]

end BioCompiler
