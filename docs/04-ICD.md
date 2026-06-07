# DOC-04: Interface Control Document (ICD)

| Field | Value |
|---|---|
| **Document ID** | DOC-04 |
| **Version** | 12.0.0 |
| **Status** | Current |
| **Date** | 2026-06-07 |
| **Prepared By** | BioCompiler Project Team |
| **Reviewed By** | [TBD at baseline review] |
| **Approved By** | [TBD at baseline review] |
| **Standard** | MIL-STD-2549 / IEEE 1320 |

---

## 1. Introduction

### 1.1 Purpose

This Interface Control Document (ICD) defines the exact contracts for every Intermediate Representation (IR) schema, component Application Programming Interface (API), and Foreign Function Interface (FFI) boundary within the BioCompiler system. It is prepared in accordance with MIL-STD-2549 and IEEE 1320 (IEEE Standard for Functional Modeling Language — Syntax and Semantics for IDEF0) interface control practices.

The purpose of this document is to provide a single, authoritative reference that enables:

1. **Developers** to implement components independently against defined contracts, without requiring access to the source code of other components.
2. **Integrators** to compose components into pipelines with confidence that inputs and outputs are compatible.
3. **Verifiers** to test each interface in isolation using the defined schemas, error conditions, and invariants.
4. **External tool authors** to implement FFI adapters that conform to the adapter contract and correctly fill SLOT fields.

Every interface defined herein is **versioned**, **typed**, and **testable**. No component SHALL access another component's internal state; all inter-component communication is through the well-typed IR records published on the IR Bus.

### 1.2 Scope

This document covers:

- **IR Schema Definitions** (Section 2): The complete Protocol Buffer (protobuf) definitions for all four IR levels — IR-Seq (IF-DATA-01), IR-Peptide (IF-DATA-02), IR-Structure (IF-DATA-03), and IR-Circuit (IF-DATA-04) — including all messages, enums, required fields, optional fields, SLOT fields, semantic invariants, producers, and consumers.
- **Component Interface Specifications** (Section 3): The exact Python function signatures, input/output types, error conditions, determinism guarantees, and performance requirements for all ten components (COMP-01 through COMP-10), each assigned a unique interface identifier (IF-01 through IF-10).
- **FFI Adapter Contract** (Section 4): The abstract adapter interface, the folding adapter contract (AlphaFold/ColabFold), and the PTM adapter contract (NetPhos/PhosphoSitePlus), including timeout policies, validation rules, and fallback behavior.
- **CLI Interface** (Section 5): The command structure, option specifications, and exit code definitions for the BioCompiler command-line interface.

Out of scope:

- The internal algorithms and data structures of individual components (see DOC-03: Software Design Document).
- The architecture rationale for interface design decisions (see DOC-06: Design Rationale).
- Test cases for interface compliance (see DOC-05: Software Verification & Validation Plan).
- The API of external tools invoked via FFI (these are treated as black boxes with defined input/output contracts).

### 1.3 References

| Ref ID | Document / Source | Description |
|---|---|---|
| REF-01 | DOC-01: Software Requirements Specification | Defines all functional and non-functional requirements that these interfaces must satisfy. |
| REF-02 | DOC-02: Software Architecture Document | Defines the component decomposition, IR Bus, and data flow that these interfaces realize. |
| REF-03 | DOC-03: Software Design Document | Provides the detailed algorithms and data structures that implement these interfaces. |
| REF-04 | DOC-05: Software Verification & Validation Plan | Specifies how each interface contract is tested and verified. |
| REF-05 | DOC-06: Design Rationale | Justifies interface design decisions and documents alternatives considered. |
| REF-06 | DOC-08: Traceability Matrix | Maps each interface to its implementing component, satisfying requirement, and test case. |
| REF-07 | MIL-STD-2549 | Military Standard: Interface Control Document Preparation. |
| REF-08 | IEEE 1320-1998 | IEEE Standard for Functional Modeling Language — Syntax and Semantics for IDEF0. |
| REF-09 | Protocol Buffers Language Specification v3 | Google's language-neutral, platform-neutral extensible mechanism for serializing structured data. |
| REF-10 | Python 3.10+ Type Hint Specification | PEP 484, PEP 604, and related PEPs governing the type annotations used in interface signatures. |
| REF-11 | GENCODE v44 (GRCh38) | Source of splice site annotations referenced by IR-Seq. |
| REF-12 | Codon Usage Database (Kazusa) | Source of organism-specific codon usage tables referenced by IR-Peptide codon assignments. |
| REF-13 | REBASE | Source of restriction enzyme recognition site sequences referenced by IR-Seq. |
| REF-14 | AlphaFold Protein Structure Database | External folding oracle invoked via the Folding FFI Adapter. |
| REF-15 | NetPhos 3.1 Server | External PTM prediction oracle invoked via the PTM FFI Adapter. |

---

## 2. IR Schema Definitions

All IR schemas are defined in Protocol Buffers (proto3 syntax). Each schema is versioned independently following semantic versioning (MAJOR.MINOR.PATCH). The schema version is carried as a required field in every IR record; the pipeline rejects records with unknown major versions.

### 2.1 IR-Seq (IF-DATA-01)

**Schema File:** `ir/proto/ir_seq_v1.proto`
**Schema Version:** 1.0.0
**Interface ID:** IF-DATA-01
**Description:** Nucleotide sequence with annotated regions including exons, introns, splice sites, regulatory elements, UTRs, restriction sites, instability motifs, and variant calls. Each splice isoform from the Splicing Engine (COMP-02) is a distinct IR-Seq record with a unique isoform ID.

#### 2.1.1 Enum Definitions

```protobuf
// Coordinate system used to number positions in the sequence.
enum CoordinateSystem {
  COORDINATE_SYSTEM_UNSPECIFIED = 0;  // Must not appear in valid records
  ZERO_BASED = 1;                     // 0-based, half-open intervals [start, end)
  ONE_BASED   = 2;                    // 1-based, closed intervals [start, end]
  GENOMIC     = 3;                    // 1-based genomic coordinates with chromosome reference
}

// Strand of the nucleotide sequence.
enum Strand {
  STRAND_UNSPECIFIED = 0;  // Must not appear in valid records
  FORWARD = 1;             // Sense / coding strand (5' → 3')
  REVERSE = 2;             // Antisense / template strand (3' → 5')
}

// Reading frame offset within the first codon.
enum ReadingFrame {
  READING_FRAME_UNSPECIFIED = 0;  // Must not appear in valid records
  FRAME_0 = 1;  // First codon starts at position 0 (modulo 3)
  FRAME_1 = 2;  // First codon starts at position 1 (modulo 3)
  FRAME_2 = 3;  // First codon starts at position 2 (modulo 3)
}
```

#### 2.1.2 Message Definitions

```protobuf
// A contiguous range of positions in the nucleotide sequence.
message PositionRange {
  int32 start = 1;  // Inclusive start position (interpretation depends on CoordinateSystem)
  int32 end   = 2;  // Exclusive end for ZERO_BASED; inclusive end for ONE_BASED/GENOMIC
}

// A splice donor or acceptor site annotation.
message SpliceSite {
  enum SpliceSiteType {
    SPLICE_SITE_TYPE_UNSPECIFIED = 0;
    DONOR    = 1;   // 5' splice site (GT consensus)
    ACCEPTOR = 2;   // 3' splice site (AG consensus)
    BRANCH_POINT = 3;  // Branch point adenosine
  }

  SpliceSiteType site_type  = 1;   // Type of splice site
  PositionRange  position   = 2;   // Position range of the consensus motif
  string         consensus  = 3;   // Consensus sequence (e.g., "GT", "AG")
  float          score      = 4;   // Splice site strength score (e.g., MaxEntScan)
  string         intron_id  = 5;   // ID of the intron this site borders
}

// A regulatory element annotation (ESE, ESS, ISE, ISS).
message RegulatoryElement {
  enum RegulatoryElementType {
    REGULATORY_ELEMENT_TYPE_UNSPECIFIED = 0;
    ESE  = 1;  // Exonic Splicing Enhancer
    ESS  = 2;  // Exonic Splicing Silencer
    ISE  = 3;  // Intronic Splicing Enhancer
    ISS  = 4;  // Intronic Splicing Silencer
  }

  RegulatoryElementType element_type = 1;  // Type of regulatory element
  PositionRange         position      = 2;  // Position range of the element
  string                motif         = 3;  // Nucleotide motif sequence
  float                 score         = 4;  // Regulatory strength score
  string                binding_factor = 5; // Known binding factor (e.g., "SF2/ASF", "hnRNP A1")
}

// A restriction enzyme recognition site.
message RestrictionSite {
  string         enzyme_name = 1;  // Name of the restriction enzyme (e.g., "EcoRI")
  PositionRange  position    = 2;  // Position range of the recognition site
  string         recognition_sequence = 3;  // Recognition sequence (e.g., "GAATTC")
  bool           is_palindromic = 4;  // Whether the enzyme recognizes a palindromic sequence
}

// An mRNA instability motif.
message InstabilityMotif {
  enum InstabilityMotifType {
    INSTABILITY_MOTIF_TYPE_UNSPECIFIED = 0;
    AUUUA     = 1;  // AU-rich element (ARE) — single copy
    AUUUA_X3  = 2;  // Tandem AUUUA repeats (3+ copies)
    U_RICH    = 3;  // U-rich region (≥6 consecutive U residues)
  }

  InstabilityMotifType motif_type = 1;  // Type of instability motif
  PositionRange        position   = 2;  // Position range of the motif
  string               sequence   = 3;  // Exact nucleotide sequence of the motif
}

// A variant call annotation (for input sequences with known variants).
message VariantCall {
  enum VariantType {
    VARIANT_TYPE_UNSPECIFIED = 0;
    SNP       = 1;  // Single nucleotide polymorphism
    INSERTION = 2;  // Insertion
    DELETION  = 3;  // Deletion
    MNV       = 4;  // Multi-nucleotide variant
  }

  string      variant_id   = 1;  // Identifier (e.g., rsID from dbSNP)
  VariantType variant_type = 2;  // Type of variant
  int32       position     = 3;  // 0-based position of the variant
  string      ref_allele   = 4;  // Reference allele
  string      alt_allele   = 5;  // Alternate allele
  float       frequency    = 6;  // Allele frequency in source population (if known)
  string      source       = 7;  // Source database (e.g., "dbSNP", "COSMIC")
}

// Untranslated region annotations.
message UTRs {
  PositionRange five_prime_utr = 1;  // 5' UTR position range (before start codon)
  PositionRange three_prime_utr = 2;  // 3' UTR position range (after stop codon)
  bool          has_secis = 3;        // Whether the 3' UTR contains a SECIS element
  PositionRange secis_position = 4;   // Position of the SECIS element (if present)
}

// Provenance metadata for the IR record.
message Provenance {
  string tool_name    = 1;  // Name of the producing tool (e.g., "BioCompiler-Scanner")
  string tool_version = 2;  // Version of the producing tool (e.g., "1.0.0")
  string timestamp    = 3;  // ISO 8601 timestamp of record creation
  string run_id       = 4;  // UUID of the pipeline run
  string input_hash   = 5;  // SHA-256 hash of the input that produced this record
  map<string, string> parameters = 6;  // Key-value pairs of parameters used
}
```

#### 2.1.3 Main Message: IRSeq

```protobuf
// IR-Seq: Nucleotide-level intermediate representation.
// Produced by COMP-01 (Scanner) and COMP-02 (Splicing Engine).
// Consumed by COMP-02, COMP-03, COMP-05, COMP-09, COMP-10.
message IRSeq {
  // ── Schema metadata ──────────────────────────────────────────────────
  string schema_version = 1;  // Schema version (e.g., "1.0.0") — REQUIRED

  // ── Required fields ──────────────────────────────────────────────────
  string          sequence          = 2;   // Full nucleotide sequence (IUPAC: A/C/G/T/U/N) — REQUIRED
  CoordinateSystem coordinate_system = 3;  // Coordinate system for all positions — REQUIRED
  Strand           strand            = 4;  // Strand of the sequence — REQUIRED
  ReadingFrame     reading_frame     = 5;  // Reading frame offset — REQUIRED
  float            gc_content        = 6;  // GC content as a fraction [0.0, 1.0] — REQUIRED
  Provenance       provenance        = 7;  // Provenance metadata — REQUIRED

  // ── Sequence identifiers ─────────────────────────────────────────────
  string sequence_id  = 8;   // Unique identifier for this sequence record
  string isoform_id   = 9;   // Unique identifier for this splice isoform (set by COMP-02)
  string gene_name    = 10;  // Gene name (if known)
  string organism     = 11;  // Source organism (e.g., "Homo sapiens")

  // ── Optional annotation fields ───────────────────────────────────────
  repeated PositionRange      exon_boundaries     = 20;  // Exon position ranges (set by COMP-02)
  repeated PositionRange      intron_boundaries   = 21;  // Intron position ranges (set by COMP-02)
  repeated SpliceSite         splice_sites        = 22;  // Splice donor/acceptor/branch point sites
  repeated RegulatoryElement  regulatory_elements = 23;  // ESE/ESS/ISE/ISS elements
  UTRs                        utrs                = 24;  // 5' and 3' UTR annotations
  repeated RestrictionSite    restriction_sites   = 25;  // Restriction enzyme recognition sites
  repeated InstabilityMotif   instability_motifs  = 26;  // mRNA instability motifs
  repeated VariantCall        variant_calls       = 27;  // Known variant calls

  // ── Coding region annotation ─────────────────────────────────────────
  PositionRange coding_region    = 30;  // Position range of the coding sequence (CDS)
  PositionRange start_codon      = 31;  // Position of the start codon (ATG)
  PositionRange stop_codon       = 32;  // Position of the stop codon
  string        kozak_consensus  = 33;  // Kozak consensus sequence at initiation site

  // ── Computed properties ──────────────────────────────────────────────
  float         cai              = 40;  // Codon Adaptation Index [0.0, 1.0]
  string        cell_type        = 41;  // Cellular context (e.g., "HEK293T")
  int32         sequence_length  = 42;  // Length of the sequence in nucleotides
}
```

#### 2.1.4 Semantic Invariants

| ID | Invariant | Enforced By | Verification Method |
|---|---|---|---|
| INV-SEQ-01 | The sequence length in nucleotides is a positive integer. For coding regions, the coding region length is a positive multiple of 3. | COMP-01 | `len(sequence) > 0` and `len(coding_region) % 3 == 0` |
| INV-SEQ-02 | Every splice site annotation has a matching donor-acceptor pair within the same contig. Donor sites are followed (downstream) by a corresponding acceptor site. | COMP-02 | Each `SpliceSite` with `DONOR` type has a matching `SpliceSite` with `ACCEPTOR` type and the same `intron_id`. |
| INV-SEQ-03 | Isoform IDs are unique within a pipeline run. No two `IRSeq` records in the same `run_id` share the same `isoform_id`. | COMP-02 | Uniqueness check against all published IR-Seq records for the run. |
| INV-SEQ-04 | GC content is within [0.0, 1.0] and equals the computed GC fraction of the sequence: `gc_content == (G_count + C_count) / len(sequence)`. | COMP-01 | Recompute GC fraction and assert equality (within ±0.001 floating-point tolerance). |

#### 2.1.5 Producers and Consumers

| Role | Component | Description |
|---|---|---|
| **Producer** | COMP-01 (Scanner) | Produces the initial IR-Seq from raw nucleotide input, with tokenized annotations (splice sites, restriction sites, instability motifs, etc.). Does not set `exon_boundaries` or `intron_boundaries`. |
| **Producer** | COMP-02 (Splicing Engine) | Produces one or more IR-Seq records per input isoform, each with `exon_boundaries`, `intron_boundaries`, and `isoform_id` populated. May refine splice site annotations based on grammar parsing. |
| **Consumer** | COMP-02 (Splicing Engine) | Consumes IR-Seq from COMP-01 to run NDFST splicing parse. |
| **Consumer** | COMP-03 (Translation Engine) | Consumes IR-Seq (one per isoform) to translate coding regions. |
| **Consumer** | COMP-05 (Type System) | Consumes IR-Seq to evaluate sequence-level type predicates (SpliceCorrect, NoCrypticSplice, GCInRange, NoRestrictionSite, NoInstabilityMotif, InFrame). |
| **Consumer** | COMP-09 (Mutation Explorer) | Consumes IR-Seq to decompose the gene into grammar nonterminals for mutation enumeration. |
| **Consumer** | COMP-10 (ORF Analyzer) | Consumes IR-Seq to identify overlapping reading frames and compute shared constraint sets. |

---

### 2.2 IR-Peptide (IF-DATA-02)

**Schema File:** `ir/proto/ir_peptide_v1.proto`
**Schema Version:** 1.0.0
**Interface ID:** IF-DATA-02
**Description:** Amino acid chain with codon provenance, signal peptide annotations, domain boundaries, PTM site predictions (SLOT), secondary structure predictions (SLOT), and other peptide-level annotations.

#### 2.2.1 Message Definitions

```protobuf
// Mapping from a codon to its translated amino acid.
message CodonAssignment {
  int32  codon_position  = 1;  // 0-based position of the first nucleotide of the codon in the source IR-Seq
  string codon            = 2;  // Three-letter nucleotide codon (e.g., "AUG", "GCU")
  string amino_acid       = 3;  // Single-letter amino acid code (e.g., "M", "A", "U" for selenocysteine)
  float  codon_usage_freq = 4;  // Codon usage frequency for this amino acid in the target organism
}

// A protein domain boundary annotation.
message DomainBoundary {
  string domain_name   = 1;  // Domain identifier (e.g., "SH3", "kinase", "zf-C2H2")
  int32  start_residue = 2;  // Start residue position (0-based) in the amino acid sequence
  int32  end_residue   = 3;  // End residue position (0-based, inclusive)
  string source_db     = 4;  // Source database (e.g., "Pfam", "InterPro", "SMART")
  float  confidence    = 5;  // Domain prediction confidence [0.0, 1.0]
}

// A post-translational modification (PTM) site prediction.
// This is a SLOT field — filled by FFI adapters (NetPhos, PhosphoSitePlus, etc.).
message PTMSite {
  enum PTMType {
    PTM_TYPE_UNSPECIFIED = 0;
    PHOSPHORYLATION = 1;
    GLYCOSYLATION   = 2;  // N-linked or O-linked
    ACETYLATION     = 3;
    UBIQUITINATION  = 4;
    METHYLATION     = 5;
    SUMOYLATION     = 6;
    NITROSYLATION   = 7;
    HYDROXYLATION   = 8;
  }

  PTMType ptm_type        = 1;  // Type of post-translational modification
  int32   residue_position = 2;  // 0-based position of the modified residue
  string  residue          = 3;  // Single-letter amino acid code of the modified residue
  float   score            = 4;  // Prediction confidence score [0.0, 1.0]
  string  predictor        = 5;  // Name of the prediction tool (e.g., "NetPhos-3.1")
  string  evidence         = 6;  // Evidence type (e.g., "predicted", "experimentally_verified")
}

// Flag indicating selenocysteine insertion at a UGA codon.
message SelenocysteineFlag {
  int32 codon_position     = 1;  // 0-based position of the UGA codon
  int32 residue_position   = 2;  // 0-based position of the selenocysteine in the peptide
  bool  secis_verified     = 3;  // Whether a SECIS element was found in the 3' UTR
  string secis_type        = 4;  // Type of SECIS element (e.g., "Form 1", "Form 2")
}

// Warning about a potential programmed ribosomal frameshift.
message FrameshiftWarning {
  int32  position   = 1;  // 0-based position where frameshift may occur
  string motif      = 2;  // Frameshift motif sequence (e.g., "UUUAAAC" for -1 frameshift)
  int32  direction  = 3;  // Frameshift direction (-1 or +1)
  float  likelihood = 4;  // Not a probability — a deterministic binary annotation (1.0 = motif present, 0.0 = absent)
}

// Secondary structure prediction for a residue range.
// This is a SLOT field — filled by FFI adapters or internal heuristics.
message SecondaryStructurePred {
  enum StructureType {
    STRUCTURE_TYPE_UNSPECIFIED = 0;
    HELIX       = 1;  // Alpha-helix
    SHEET       = 2;  // Beta-sheet
    COIL        = 3;  // Random coil / loop
  }

  int32         start_residue = 1;  // Start residue position (0-based)
  int32         end_residue   = 2;  // End residue position (0-based, inclusive)
  StructureType structure_type = 3;  // Predicted secondary structure type
  float         confidence    = 4;  // Prediction confidence [0.0, 1.0]
  string        predictor     = 5;  // Name of the prediction method
}

// Solvent accessibility prediction for a residue.
// This is a SLOT field — filled by FFI adapters.
message SolventAccessibility {
  int32  residue_position = 1;  // 0-based residue position
  float  relative_asa    = 2;  // Relative accessible surface area [0.0, 1.0]
  bool   is_exposed      = 3;  // Binary: True if relative_asa >= 0.25 (exposed threshold)
  string predictor        = 4;  // Name of the prediction method
}

// Intrinsically disordered region prediction.
// This is a SLOT field — filled by FFI adapters.
message DisorderRegion {
  int32  start_residue = 1;  // Start residue position (0-based)
  int32  end_residue   = 2;  // End residue position (0-based, inclusive)
  float  score         = 3;  // Disorder propensity score [0.0, 1.0]
  string predictor     = 4;  // Name of the prediction method
}
```

#### 2.2.2 Main Message: IRPeptide

```protobuf
// IR-Peptide: Amino-acid-level intermediate representation.
// Produced by COMP-03 (Translation Engine).
// Consumed by COMP-04, COMP-05, COMP-09.
message IRPeptide {
  // ── Schema metadata ──────────────────────────────────────────────────
  string schema_version = 1;  // Schema version (e.g., "1.0.0") — REQUIRED

  // ── Required fields ──────────────────────────────────────────────────
  string amino_acid_sequence = 2;   // Full amino acid sequence (single-letter codes) — REQUIRED
  int32  peptide_length      = 3;   // Length of the amino acid sequence in residues — REQUIRED
  string source_isoform_id   = 4;   // ID of the source IR-Seq isoform — REQUIRED
  ReadingFrame reading_frame = 5;   // Reading frame used for translation — REQUIRED
  Provenance   provenance    = 6;   // Provenance metadata — REQUIRED

  // ── Codon provenance ─────────────────────────────────────────────────
  repeated CodonAssignment codon_assignments = 10;  // One per amino acid residue

  // ── Optional annotations ─────────────────────────────────────────────
  repeated DomainBoundary  domain_boundaries  = 20;  // Protein domain annotations
  bool    has_signal_peptide = 21;                    // Whether a signal peptide was detected
  int32   signal_peptide_end = 22;                    // Cleavage site position (0-based)

  // ── SLOT fields (initially empty; filled by FFI adapters) ────────────
  repeated SecondaryStructurePred secondary_structure_pred = 30;  // SLOT: Filled by folding adapter
  repeated float                 ss_confidence             = 31;  // SLOT: Per-residue confidence for SS prediction
  repeated PTMSite               ptm_sites                 = 32;  // SLOT: Filled by PTM adapter
  repeated SolventAccessibility  solvent_accessibility     = 33;  // SLOT: Filled by folding adapter
  repeated DisorderRegion        disorder_regions          = 34;  // SLOT: Filled by folding adapter

  // ── Flags and warnings ───────────────────────────────────────────────
  repeated SelenocysteineFlag selenocysteine_flags  = 40;  // Selenocysteine insertion events
  repeated FrameshiftWarning  frameshift_warnings   = 41;  // Potential programmed frameshifts
  bool   is_complete          = 42;  // True if translation reached a stop codon; False if premature or truncated
  bool   has_premature_stop   = 43;  // True if a stop codon was encountered before the expected end

  // ── Sequence identifiers ─────────────────────────────────────────────
  string peptide_id   = 50;  // Unique identifier for this peptide record
  string gene_name    = 51;  // Gene name (inherited from source IR-Seq)
  string organism     = 52;  // Organism (inherited from source IR-Seq)

  // ── Computed properties ──────────────────────────────────────────────
  float  molecular_weight_kda = 60;  // Estimated molecular weight in kilodaltons
  float  isoelectric_point    = 61;  // Estimated isoelectric point (pI)
  int32  num_cysteines        = 62;  // Count of cysteine residues
}
```

#### 2.2.3 Semantic Invariants

| ID | Invariant | Enforced By | Verification Method |
|---|---|---|---|
| INV-PEP-01 | Peptide length equals `(coding_region_length / 3) - 1` (stop codon excluded from the amino acid sequence). If the coding region is not divisible by 3, the peptide is truncated at the last complete codon. | COMP-03 | `peptide_length == (len(coding_region) // 3) - 1` |
| INV-PEP-02 | Every amino acid residue has a back-reference to its source codon in IR-Seq. The number of `codon_assignments` equals `peptide_length`. | COMP-03 | `len(codon_assignments) == peptide_length` and each `CodonAssignment` references a valid position in the source IR-Seq. |
| INV-PEP-03 | No duplicate isoform references within the same peptide record. A peptide is derived from exactly one source isoform. | COMP-03 | `source_isoform_id` is singular and non-empty. |

---

### 2.3 IR-Structure (IF-DATA-03)

**Schema File:** `ir/proto/ir_structure_v1.proto`
**Schema Version:** 1.0.0
**Interface ID:** IF-DATA-03
**Description:** Peptide with secondary/tertiary structure annotations from FFI oracles (AlphaFold, ColabFold, RoseTTAFold). Includes per-atom coordinates, confidence scores (pLDDT), predicted aligned error matrices, PTM site predictions, validation flags, and folding energy estimates.

#### 2.3.1 Message Definitions

```protobuf
// A single atom coordinate in the predicted structure.
message AtomCoordinate {
  string atom_name      = 1;  // PDB atom name (e.g., "CA", "N", "C", "O", "CB")
  string residue_name   = 2;  // Three-letter amino acid code (e.g., "ALA", "GLY")
  int32  residue_number = 3;  // Residue sequence number (0-based)
  double x              = 4;  // X coordinate in Angstroms
  double y              = 5;  // Y coordinate in Angstroms
  double z              = 6;  // Z coordinate in Angstroms
  float  occupancy      = 7;  // Occupancy value (typically 1.0 for predicted structures)
  float  b_factor       = 8;  // B-factor / temperature factor (often used for pLDDT)
}

// Confidence score for a predicted structure.
message ConfidenceScore {
  int32  residue_position = 1;  // 0-based residue position
  float  plddt            = 2;  // Predicted Local Distance Difference Test score [0.0, 100.0]
  float  pae              = 3;  // Predicted Aligned Error in Angstroms [0.0, +inf)
  float  global_confidence = 4; // Global confidence metric for the entire structure [0.0, 1.0]
}

// Predicted aligned error for a pair of residues.
message PredictedAlignedError {
  int32  residue_i = 1;  // First residue position (0-based)
  int32  residue_j = 2;  // Second residue position (0-based)
  float  pae_value = 3;  // Predicted aligned error in Angstroms
}

// A validation flag for a structural element.
message ValidationFlag {
  enum FlagType {
    FLAG_TYPE_UNSPECIFIED = 0;
    CLASH          = 1;  // Steric clash between atoms
    RAMACHANDRAN   = 2;  // Ramachandran outlier
    ROTAMER        = 3;  // Rotamer outlier
    BOND_LENGTH    = 4;  // Bond length outlier
    BOND_ANGLE     = 5;  // Bond angle outlier
    PEPTIDE_BOND   = 6;  // Non-planar peptide bond
  }

  FlagType flag_type      = 1;  // Type of validation flag
  int32    residue_position = 2; // Residue position (0-based)
  string   description     = 3;  // Human-readable description
  float    deviation       = 4;  // Deviation from ideal (in appropriate units)
}

// A folding energy estimate for the predicted structure.
message EnergyEstimate {
  enum EnergyType {
    ENERGY_TYPE_UNSPECIFIED = 0;
    TOTAL          = 1;  // Total energy
  VAN_DER_WAALS  = 2;  // Van der Waals energy
  ELECTROSTATIC  = 3;  // Electrostatic energy
  SOLVATION      = 4;  // Solvation energy
  DIHEDRAL       = 5;  // Dihedral angle energy
  }

  EnergyType energy_type = 1;  // Type of energy estimate
  double     value       = 2;  // Energy value in kcal/mol
  string     method      = 3;  // Estimation method (e.g., "Rosetta", "MMTF", "OpenMM")
  float      confidence  = 4;  // Confidence in the estimate [0.0, 1.0]
}
```

#### 2.3.2 Main Message: IRStructure

```protobuf
// IR-Structure: Structure-level intermediate representation.
// Produced by COMP-04 (FFI Manager) after invoking folding oracles.
// Consumed by COMP-05, COMP-06, COMP-07.
message IRStructure {
  // ── Schema metadata ──────────────────────────────────────────────────
  string schema_version = 1;  // Schema version (e.g., "1.0.0") — REQUIRED

  // ── Required fields ──────────────────────────────────────────────────
  string           source_peptide_id   = 2;  // ID of the parent IR-Peptide — REQUIRED
  string           amino_acid_sequence = 3;  // Amino acid sequence (inherited from IR-Peptide) — REQUIRED
  Provenance       provenance          = 4;  // Provenance metadata — REQUIRED
  float            mean_plddt          = 5;  // Mean pLDDT score across all residues — REQUIRED

  // ── SLOT fields (filled by folding FFI adapter) ──────────────────────
  repeated AtomCoordinate      atom_coordinates     = 20;  // SLOT: Per-atom 3D coordinates
  repeated ConfidenceScore     residue_confidences  = 21;  // SLOT: Per-residue confidence scores (pLDDT + PAE)
  repeated PredictedAlignedError pae_matrix         = 22;  // SLOT: Full PAE matrix (N×N residues)
  repeated ValidationFlag      validation_flags     = 23;  // SLOT: Structural validation flags
  repeated EnergyEstimate      energy_estimates     = 24;  // SLOT: Folding energy estimates

  // ── Secondary structure annotation ───────────────────────────────────
  repeated SecondaryStructurePred secondary_structure = 30;  // Secondary structure assignments (from folding)
  repeated SolventAccessibility  solvent_accessibility = 31;  // Per-residue solvent accessibility
  repeated DisorderRegion        disorder_regions      = 32;  // Disordered region predictions

  // ── PTM annotations (from PTM FFI adapter) ───────────────────────────
  repeated PTMSite ptm_sites = 40;  // Post-translational modification site predictions

  // ── Model metadata ───────────────────────────────────────────────────
  string model_identifier  = 50;  // Identifier for the predicted model (e.g., "AF2-RANKED_0")
  int32  model_rank        = 51;  // Rank of this model among predictions (1 = best)
  string folding_tool       = 52;  // Name of the folding tool (e.g., "AlphaFold2", "ColabFold")
  string folding_version   = 53;  // Version of the folding tool
  string msa_depth         = 54;  // Depth of the multiple sequence alignment used (if applicable)

  // ── Quality metrics ──────────────────────────────────────────────────
  float  iptm_score        = 60;  // Interface pTM score for multimers [0.0, 1.0]
  int32  num_residues       = 61;  // Number of residues in the predicted structure
  int32  num_atoms          = 62;  // Number of atoms in the predicted structure
}
```

#### 2.3.3 Semantic Invariants

| ID | Invariant | Enforced By | Verification Method |
|---|---|---|---|
| INV-STR-01 | Every structural element has a confidence score in [0.0, 1.0] (normalized from pLDDT scale: `pLDDT / 100.0`). `mean_plddt` is in [0.0, 100.0]. | COMP-04 | `0.0 <= mean_plddt <= 100.0` and for each `ConfidenceScore`, `0.0 <= plddt <= 100.0` and `0.0 <= pae`. |
| INV-STR-02 | PTM site predictions reference valid residue positions in the parent IR-Peptide. Each `PTMSite.residue_position` is in `[0, peptide_length)`. | COMP-04 | `0 <= ptm_site.residue_position < peptide_length` for all PTM sites. |

---

### 2.4 IR-Circuit (IF-DATA-04)

**Schema File:** `ir/proto/ir_circuit_v1.proto`
**Schema Version:** 1.0.0
**Interface ID:** IF-DATA-04
**Description:** Directed acyclic graph (DAG) of gene nodes with interaction edges. Each node references a per-gene certificate. Edges carry interaction constraints (promoter interference, metabolic burden, resource competition, splicing interference, RNA-RNA interaction). Produced by COMP-08 (Compositional Verifier) and consumed by COMP-05, COMP-07, COMP-08.

#### 2.4.1 Message Definitions

```protobuf
// A promoter element driving transcription of a gene.
message Promoter {
  string promoter_id    = 1;  // Unique identifier for the promoter
  string promoter_name  = 2;  // Human-readable name (e.g., "CMV", "EF1α", "T7")
  bool   is_constitutive = 3; // Whether the promoter is constitutive (always active)
  repeated string tf_activators  = 4;  // Transcription factors that activate this promoter
  repeated string tf_repressors  = 5;  // Transcription factors that repress this promoter
  float  strength       = 6;  // Relative promoter strength (arbitrary units)
  string organism       = 7;  // Source organism for the promoter element
}

// A terminator element ending transcription.
message Terminator {
  string terminator_id   = 1;  // Unique identifier for the terminator
  string terminator_name = 2;  // Human-readable name (e.g., "SV40 polyA", "bGH polyA")
  float  efficiency      = 3;  // Termination efficiency [0.0, 1.0]
  string terminator_type = 4;  // Type (e.g., "rho_independent", "rho_dependent", "polyA_signal")
}

// An overlapping open reading frame annotation.
message OverlappingORF {
  string orf_id         = 1;  // Unique identifier for the ORF
  int32  start_position = 2;  // Start position in the circuit sequence (0-based)
  int32  end_position   = 3;  // End position (exclusive for 0-based)
  int32  frame_offset   = 4;  // Reading frame offset (0, 1, or 2)
  string protein_name   = 5;  // Name of the encoded protein (if known)
  bool   is_annotated   = 6;  // Whether this ORF is an annotated (known) ORF
}

// A gene node in the circuit graph.
message GeneNode {
  string gene_id          = 1;  // Unique identifier for the gene
  string gene_name        = 2;  // Human-readable gene name
  Promoter promoter       = 3;  // Promoter driving this gene
  Terminator terminator   = 4;  // Terminator ending this gene's transcript
  string certificate_id   = 5;  // ID of the per-gene guarantee certificate from COMP-07
  string ir_seq_id        = 6;  // ID of the gene's IR-Seq record
  string ir_peptide_id    = 7;  // ID of the gene's IR-Peptide record
  string ir_structure_id  = 8;  // ID of the gene's IR-Structure record (if available)
  repeated OverlappingORF overlapping_orfs = 9;  // Overlapping ORFs in this gene's region
}

// An interaction edge between two genes.
message InteractionEdge {
  enum InteractionType {
    INTERACTION_TYPE_UNSPECIFIED = 0;
    PROMOTER_CONFLICT     = 1;  // Unintentional TF-mediated regulation
    RESOURCE_COMPETITION  = 2;  // Ribosome/resource demand exceeds capacity
    SPLICING_INTERFERENCE = 3;  // Cryptic splice site interference
    RNA_RNA_INTERACTION   = 4;  // Complementary transcript regions
    METABOLIC_BURDEN      = 5;  // Combined metabolic load
  }

  string source_gene_id     = 1;  // ID of the source gene
  string target_gene_id     = 2;  // ID of the target gene
  InteractionType interaction_type = 3;  // Type of interaction
  string verdict            = 4;  // Verdict: "PASS", "FAIL", or "UNCERTAIN"
  string evidence           = 5;  // Evidence description
  float  severity           = 6;  // Severity score [0.0, 1.0] (only meaningful for FAIL/UNCERTAIN)
}
```

#### 2.4.2 Main Message: IRCircuit

```protobuf
// IR-Circuit: Circuit-level intermediate representation.
// Produced by COMP-08 (Compositional Verifier).
// Consumed by COMP-05, COMP-07, COMP-08.
message IRCircuit {
  // ── Schema metadata ──────────────────────────────────────────────────
  string schema_version = 1;  // Schema version (e.g., "1.0.0") — REQUIRED

  // ── Required fields ──────────────────────────────────────────────────
  string circuit_id     = 2;  // Unique identifier for the circuit — REQUIRED
  string organism       = 3;  // Target organism — REQUIRED
  string cell_type      = 4;  // Target cellular context — REQUIRED
  Provenance provenance = 5;  // Provenance metadata — REQUIRED

  // ── Circuit topology ─────────────────────────────────────────────────
  string topology_type  = 10;  // "linear" or "circular"
  repeated GeneNode genes = 11;  // Gene nodes in the circuit
  repeated InteractionEdge interactions = 12;  // Interaction edges between genes

  // ── Circuit-level properties ─────────────────────────────────────────
  int32  num_genes            = 20;  // Number of genes in the circuit
  int32  total_sequence_length = 21; // Total nucleotide length of the circuit construct
  float  total_gc_content     = 22;  // Overall GC content of the circuit construct
  string overall_verdict      = 23;  // Composed verdict: "PASS", "FAIL", or "UNCERTAIN"

  // ── Certificate reference ────────────────────────────────────────────
  string circuit_certificate_id = 30;  // ID of the circuit-level guarantee certificate
}
```

#### 2.4.3 Semantic Invariants

| ID | Invariant | Enforced By | Verification Method |
|---|---|---|---|
| INV-CIR-01 | The circuit graph is acyclic — no circular promoter dependencies exist. A cycle would indicate mutual transcription factor regulation forming an infinite loop. | COMP-08 | Topological sort of the gene-interaction graph succeeds (no cycle detected). |
| INV-CIR-02 | Every gene node references a valid, non-expired certificate from COMP-07. The certificate's `design_id` matches the gene's `ir_seq_id` content hash. | COMP-08 | Certificate existence check + `sha256(sequence) == certificate.design_id` for each gene. |

---

## 3. Component Interface Specifications

This section defines the exact interface contract for each of the ten BioCompiler components. Each interface is assigned a unique identifier (IF-01 through IF-10). For each interface, the following are specified:

- **Component ID and name**
- **Interface ID**
- **Input specification** (types, constraints, required/optional)
- **Output specification** (types, invariants, error conditions)
- **Error conditions** (named exception types with triggering conditions)
- **Determinism guarantee** (conditions under which the interface is deterministic)
- **Performance requirement** (time and space bounds)
- **Python function signature** (exact type-annotated signature)

### 3.1 IF-01: Scanner Interface

**Component:** COMP-01 — Scanner
**Interface ID:** IF-01
**Description:** Performs DFA-based lexical analysis of nucleotide sequences, producing an IR-Seq record with tokenized annotations for start codons, stop codons, splice donor/acceptor sites, branch points, polypyrimidine tracts, Kozak consensus sequences, RNA instability motifs, restriction enzyme recognition sites, and other regulatory elements.

| Aspect | Specification |
|---|---|
| **Component** | COMP-01 (Scanner) |
| **Input** | `sequence` (str): Nucleotide sequence in IUPAC notation (A/C/G/T/U/N). `format` (Literal["fasta", "raw"]): Input format. `restriction_enzymes` (list[str] \| None): Restriction enzyme names to scan for (from REBASE). |
| **Output** | `IRSeq`: Populated with `sequence`, `coordinate_system`, `strand`, `reading_frame`, `gc_content`, `provenance`, `splice_sites`, `restriction_sites`, `instability_motifs`, `start_codon`, `stop_codon`, `kozak_consensus`, `coding_region`, `sequence_length`. |
| **Error Conditions** | `InvalidSequenceError`: Non-IUPAC character in input sequence (includes position and character). `FormatError`: Malformed FASTA header or multi-sequence FASTA (includes line number). `EmptySequenceError`: Zero-length input sequence. `UnknownEnzymeError`: Restriction enzyme name not found in REBASE. |
| **Determinism** | **Fully deterministic.** Given identical `sequence`, `format`, and `restriction_enzymes`, the scanner produces bit-identical `IRSeq` output. No randomness, no floating-point ambiguity (GC content is computed via integer count divided by sequence length). |
| **Performance** | Time: O(n × d × k) where n = sequence length, d = number of active DFAs, k = average DFA size. ≤ 1 second for n ≤ 10,000. Space: O(n + t) where t = number of tokens emitted. ≤ 1 MB for n ≤ 10,000. |

#### Python Signature

```python
from typing import Literal
from biocompiler.ir.proto.ir_seq_v1_pb2 import IRSeq

def scan(
    sequence: str,
    format: Literal["fasta", "raw"] = "raw",
    restriction_enzymes: list[str] | None = None,
) -> IRSeq:
    """Perform DFA-based lexical analysis of a nucleotide sequence.

    Args:
        sequence: Nucleotide sequence string in IUPAC notation (A/C/G/T/U/N).
            If format is "fasta", the sequence is extracted from the FASTA
            record (header line starting with '>' followed by sequence lines).
        format: Input format — "raw" for bare sequence, "fasta" for FASTA format.
        restriction_enzymes: Optional list of restriction enzyme names (from REBASE)
            to scan for. If None, no restriction site scanning is performed.

    Returns:
        IRSeq record with all tokenized annotations populated.

    Raises:
        InvalidSequenceError: If the sequence contains non-IUPAC characters.
        FormatError: If format is "fasta" and the input is malformed.
        EmptySequenceError: If the input sequence has zero length.
        UnknownEnzymeError: If a specified enzyme name is not found in REBASE.
    """
    ...
```

---

### 3.2 IF-02: Splicing Engine Interface

**Component:** COMP-02 — Splicing Engine
**Interface ID:** IF-02
**Description:** Implements an NDFST that parses the tokenized IR-Seq against the splicing grammar and computes the set of all possible splice isoforms. The cellular context parameter modulates regulatory element thresholds, producing different isoform sets for different cell types. The computation is deterministic (same isoform set for same input), but the output is set-valued (multiple isoforms possible).

| Aspect | Specification |
|---|---|
| **Component** | COMP-02 (Splicing Engine) |
| **Input** | `ir_seq` (IRSeq): Tokenized IR-Seq from COMP-01. `cell_type` (str): Target cellular context (e.g., "HEK293T"). `splicing_rules_path` (str): Path to the splicing grammar rules YAML file. |
| **Output** | `list[IRSeq]`: One IR-Seq per splice isoform, each with `exon_boundaries`, `intron_boundaries`, `isoform_id`, and `reading_frame` populated. Returns an empty list if no valid parse exists. |
| **Error Conditions** | `NoValidIsoformError`: The NDFST produces an empty isoform set — no parse path satisfies the grammar. Includes the specific grammar constraint(s) that prevent any valid parse. `UnknownCellTypeError`: The specified cell type has no entry in the cellular context configuration. `GrammarRuleError`: The splicing rules YAML is malformed or references an unknown motif. |
| **Determinism** | **Fully deterministic.** Given identical `ir_seq`, `cell_type`, and `splicing_rules_path`, the engine produces an identical set of `IRSeq` isoform records in the same order. Non-determinism of splicing biology is modeled as set-valued output; the set itself is computed deterministically. |
| **Performance** | Time: O(n × s) where n = sequence length, s = number of NDFST states. ≤ 5 seconds for n ≤ 10,000 with typical grammar. Space: O(n × i) where i = number of isoforms. ≤ 100 MB for n ≤ 10,000 with ≤ 64 isoforms. |

#### Python Signature

```python
from biocompiler.ir.proto.ir_seq_v1_pb2 import IRSeq

def parse_splicing(
    ir_seq: IRSeq,
    cell_type: str = "HEK293T",
    splicing_rules_path: str = "config/splicing_rules.yaml",
) -> list[IRSeq]:
    """Parse pre-mRNA sequence against the splicing grammar to compute all possible splice isoforms.

    Uses a Non-Deterministic Finite-State Transducer (NDFST) to explore all valid
    parse paths through the splicing grammar under the specified cellular context.
    The isoform set is complete (all valid parse paths are represented) and sound
    (every isoform satisfies the grammar).

    Args:
        ir_seq: Tokenized IR-Seq record from the Scanner (COMP-01).
        cell_type: Target cellular context. Determines regulatory element
            thresholds (ESE/ESS/ISE/ISS) and splicing factor concentrations.
            Defaults to "HEK293T".
        splicing_rules_path: Path to the YAML file containing splicing grammar
            rules, consensus sequences, and cell-type-specific thresholds.
            Defaults to "config/splicing_rules.yaml".

    Returns:
        List of IRSeq records, one per splice isoform. Each record has
        exon_boundaries, intron_boundaries, isoform_id, and reading_frame
        populated. Returns an empty list if no valid parse exists.

    Raises:
        NoValidIsoformError: If the NDFST produces no valid isoforms.
        UnknownCellTypeError: If cell_type has no entry in the cellular context config.
        GrammarRuleError: If the splicing rules YAML is malformed.
    """
    ...
```

---

### 3.3 IF-03: Translation Engine Interface

**Component:** COMP-03 — Translation Engine
**Interface ID:** IF-03
**Description:** Implements a deterministic finite-state transducer (FST) that maps each codon in a spliced mRNA sequence to its corresponding amino acid. Handles the standard genetic code, selenocysteine insertion (UGA recoding with SECIS element), pyrrolysine incorporation (UAG recoding in archaeal contexts), and detection of programmed ribosomal frameshifting motifs.

| Aspect | Specification |
|---|---|
| **Component** | COMP-03 (Translation Engine) |
| **Input** | `ir_seq` (IRSeq): A single splice isoform IR-Seq from COMP-02 (with `exon_boundaries` populated). |
| **Output** | `IRPeptide`: Populated with `amino_acid_sequence`, `peptide_length`, `source_isoform_id`, `reading_frame`, `provenance`, `codon_assignments`, `selenocysteine_flags`, `frameshift_warnings`, `is_complete`, `has_premature_stop`. |
| **Error Conditions** | `NoStartCodonError`: No start codon (ATG) found in the coding region. `FrameError`: The reading frame is inconsistent — exon boundary does not preserve the reading frame across the splice junction. `InvalidCodonError`: A codon contains non-standard nucleotides (should not occur if scanner validated input). |
| **Determinism** | **Fully deterministic.** Given an identical IR-Seq isoform, the engine produces a bit-identical IR-Peptide. The genetic code is a fixed mapping; selenocysteine and pyrrolysine handling are deterministic given the context. |
| **Performance** | Time: O(n) where n = coding region length. ≤ 10 ms for n ≤ 10,000. Space: O(n) for the codon assignment list. ≤ 500 KB for n ≤ 10,000. |

#### Python Signature

```python
from biocompiler.ir.proto.ir_seq_v1_pb2 import IRSeq
from biocompiler.ir.proto.ir_peptide_v1_pb2 import IRPeptide

def translate(
    ir_seq: IRSeq,
) -> IRPeptide:
    """Translate a spliced mRNA isoform to an amino acid sequence.

    Uses a deterministic finite-state transducer (FST) implementing the standard
    genetic code. Handles selenocysteine insertion (UGA recoding with SECIS element
    in 3' UTR), pyrrolysine incorporation (UAG recoding in archaeal contexts),
    and detection of programmed ribosomal frameshifting motifs.

    Args:
        ir_seq: A single splice isoform IR-Seq record from the Splicing Engine
            (COMP-02). Must have exon_boundaries, reading_frame, and
            coding_region populated.

    Returns:
        IRPeptide record with amino acid sequence, codon assignments,
        selenocysteine flags, frameshift warnings, and completion status.

    Raises:
        NoStartCodonError: If no start codon is found in the coding region.
        FrameError: If reading frame is inconsistent across splice junctions.
        InvalidCodonError: If a codon contains non-standard nucleotides.
    """
    ...
```

---

### 3.4 IF-04: FFI Manager Interface

**Component:** COMP-04 — FFI Manager
**Interface ID:** IF-04
**Description:** Orchestrates invocations of external tools (protein structure predictors, PTM prediction tools) through a common adapter interface. Responsible for input formatting, output parsing, invariant validation, SLOT field filling, provenance tracking, timeout handling, and fallback behavior. Treats every external tool as a non-deterministic black box.

| Aspect | Specification |
|---|---|
| **Component** | COMP-04 (FFI Manager) |
| **Input** | `ir_peptide` (IRPeptide): Peptide record with empty SLOT fields. `adapter_name` (str): Name of the registered adapter to invoke. `adapter_config` (dict): Adapter-specific configuration (timeout, parameters, model selection). |
| **Output** | For folding: `IRStructure` with SLOT fields filled. For PTM: `IRPeptide` with `ptm_sites` SLOT field filled. |
| **Error Conditions** | `AdapterNotFoundError`: No adapter registered with the given name. `ExternalToolError`: The external tool returned a non-zero exit code or exceeded the timeout. `OutputParseError`: The tool's output does not match the expected schema. `OutputValidationError`: Parsed output violates IR invariants (INV-STR-01, INV-STR-02, or INV-PEP-03). |
| **Determinism** | **Non-deterministic by design.** FFI invocations may produce different results across runs due to GPU non-determinism, MSA stochasticity, or tool version differences. Non-determinism is captured as UNCERTAIN verdicts in the type system. Provenance metadata records the exact tool version and parameters. |
| **Performance** | Folding: 10 minutes — 48 hours depending on sequence length and tool. PTM: ≤ 5 seconds for typical sequences. Timeout is configurable; default is 3600 seconds for folding and 300 seconds for PTM. |

#### Python Signatures

```python
from biocompiler.ir.proto.ir_peptide_v1_pb2 import IRPeptide
from biocompiler.ir.proto.ir_structure_v1_pb2 import IRStructure

def invoke_folding(
    ir_peptide: IRPeptide,
    adapter_name: str = "alphafold2",
    timeout: float = 3600.0,
    model_preset: str = "monomer",
    msa_depth: int | None = None,
    use_gpu: bool = True,
) -> IRStructure:
    """Invoke an external folding tool via FFI to predict protein structure.

    Formats the amino acid sequence from the IR-Peptide as input for the
    specified folding adapter (AlphaFold2, ColabFold, etc.), invokes the
    external tool, parses the output, validates invariants, fills SLOT
    fields in the IR-Structure, and records provenance.

    Args:
        ir_peptide: IR-Peptide record with amino acid sequence populated.
            SLOT fields may be empty (will not be used by folding adapter).
        adapter_name: Name of the registered folding adapter.
            Supported values: "alphafold2", "alphafold3", "colabfold",
            "rosettafold". Defaults to "alphafold2".
        timeout: Maximum wall-clock time for the external tool, in seconds.
            Defaults to 3600.0 (1 hour).
        model_preset: Folding model preset. "monomer" for single-chain,
            "multimer" for multi-chain. Defaults to "monomer".
        msa_depth: Maximum MSA depth to use (None = tool default).
        use_gpu: Whether to use GPU acceleration. Defaults to True.

    Returns:
        IRStructure record with atom coordinates, confidence scores (pLDDT),
        PAE matrix, validation flags, and energy estimates populated.

    Raises:
        AdapterNotFoundError: If adapter_name is not registered.
        ExternalToolError: If the tool times out or returns a non-zero exit code.
        OutputParseError: If the tool output cannot be parsed.
        OutputValidationError: If parsed output violates IR invariants.
    """
    ...


def invoke_ptm_prediction(
    ir_peptide: IRPeptide,
    adapter_name: str = "netphos",
    timeout: float = 300.0,
    ptm_types: list[str] | None = None,
    threshold: float = 0.5,
) -> IRPeptide:
    """Invoke an external PTM prediction tool via FFI.

    Formats the amino acid sequence from the IR-Peptide as input for the
    specified PTM adapter (NetPhos, PhosphoSitePlus, etc.), invokes the
    external tool, parses the output, validates invariants, fills the
    ptm_sites SLOT field in the IR-Peptide, and records provenance.

    Args:
        ir_peptide: IR-Peptide record with amino acid sequence populated.
        adapter_name: Name of the registered PTM adapter.
            Supported values: "netphos", "phosphositeplus", "dbptm",
            "musitedeep". Defaults to "netphos".
        timeout: Maximum wall-clock time for the external tool, in seconds.
            Defaults to 300.0 (5 minutes).
        ptm_types: Types of PTMs to predict. If None, all supported types
            are predicted. Supported values: "phosphorylation",
            "glycosylation", "acetylation", "ubiquitination",
            "methylation", "sumoylation".
        threshold: Minimum confidence score for a PTM site to be included
            in the output. Defaults to 0.5.

    Returns:
        The input IRPeptide with the ptm_sites SLOT field filled with
        PTMSite predictions meeting the threshold.

    Raises:
        AdapterNotFoundError: If adapter_name is not registered.
        ExternalToolError: If the tool times out or returns a non-zero exit code.
        OutputParseError: If the tool output cannot be parsed.
        OutputValidationError: If parsed output violates IR invariants.
    """
    ...
```

---

### 3.5 IF-05: Type System Interface

**Component:** COMP-05 — Type System
**Interface ID:** IF-05
**Description:** Performs static verification of biological correctness properties on mRNA sequences. Evaluates each declared type predicate against the mRNA and produces a three-valued verdict: PASS (guaranteed correct with derivation trace), FAIL (guaranteed incorrect with violation identification), or UNCERTAIN (cannot determine, with knowledge gap specification). This is the safety-critical component; its soundness (no false PASS) is the most important invariant of the entire system.

| Aspect | Specification |
|---|---|
| **Component** | COMP-05 (Type System) |
| **Input** | `ir_seq` (IRSeq): The mRNA sequence to check. `predicates` (list[TypePredicate]): The set of type predicates to evaluate. `target_isoform` (str \| None): The intended splice isoform ID (for SpliceCorrect predicate). |
| **Output** | `TypeCheckResult`: For each predicate, a verdict (PASS/FAIL/UNCERTAIN) with derivation trace, violation identification, or knowledge gap specification. |
| **Error Conditions** | `PredicateError`: A type predicate is malformed or references an unknown property. `MissingPrerequisiteError`: A predicate requires data that is not present in the IR record (e.g., SpliceCorrect requires isoform IDs). |
| **Determinism** | **Fully deterministic.** Given identical IR records and predicates, the type system produces identical verdicts. Three-valued composition preserves soundness (INV-TYP-01, INV-TYP-03). |
| **Performance** | Time: O(n × p) where n = sequence length, p = number of predicates. ≤ 2 seconds for n ≤ 10,000 with 13 core predicates; ≤ 10 seconds with all 33 predicates. Space: O(n) for token scanning within predicates. ≤ 2 MB. |

#### Data Structures

```python
from enum import Enum
from dataclasses import dataclass, field
from typing import Any

class Verdict(Enum):
    """Three-valued verdict for type predicate evaluation."""
    PASS = "PASS"          # Guaranteed correct (with derivation trace)
    FAIL = "FAIL"          # Guaranteed incorrect (with violation identification)
    UNCERTAIN = "UNCERTAIN"  # Cannot determine (with knowledge gap specification)


@dataclass(frozen=True)
class TypePredicate:
    """A biological correctness property to be checked against an IR record."""
    name: str                    # Predicate name (e.g., "SpliceCorrect", "NoCrypticSplice")
    parameters: dict[str, Any]   # Predicate parameters (e.g., {"cell_type": "HEK293T"})


@dataclass(frozen=True)
class DerivationStep:
    """A single step in a derivation trace for a PASS verdict."""
    step: int          # Step number (1-indexed)
    rule: str          # Rule applied (e.g., "NDFST_output_is_singleton")
    evidence: str      # Evidence for this step (e.g., "Isoform set = {isoform_1}")


@dataclass(frozen=True)
class Violation:
    """A violation identified by a FAIL verdict."""
    predicate: str       # Name of the violated predicate
    position: int | None # Position of the violation (0-based), if applicable
    rule: str            # Rule that was violated
    evidence: str        # Evidence for the violation


@dataclass(frozen=True)
class KnowledgeGap:
    """A knowledge gap that prevents resolution to PASS or FAIL."""
    predicate: str     # Name of the predicate
    description: str   # Description of what is needed to resolve
    required_data: str # Description of the data or knowledge required


@dataclass(frozen=True)
class TypeCheckResult:
    """Result of evaluating a single type predicate against an IR record."""
    predicate: TypePredicate          # The predicate that was evaluated
    verdict: Verdict                  # Three-valued verdict
    derivation: list[DerivationStep] | None = None  # Derivation trace (for PASS)
    violation: Violation | None = None              # Violation details (for FAIL)
    knowledge_gap: KnowledgeGap | None = None       # Knowledge gap (for UNCERTAIN)
```

#### Python Signature

```python
from biocompiler.ir.proto.ir_seq_v1_pb2 import IRSeq

def type_check(
    ir_seq: IRSeq,
    predicates: list[TypePredicate],
    target_isoform: str | None = None,
) -> list[TypeCheckResult]:
    """Evaluate type predicates against an mRNA sequence.

    Each predicate is evaluated independently, producing a three-valued verdict:
    PASS (with derivation trace), FAIL (with violation identification), or
    UNCERTAIN (with knowledge gap specification). Results are composed using
    three-valued AND logic: overall PASS only if all predicates PASS.

    Soundness guarantee (INV-TYP-01): If PASS is returned for a predicate,
    the property holds for the sequence under the specified conditions.

    Args:
        ir_seq: IR-Seq record to check. Must have sequence, coordinate_system,
            strand, reading_frame, gc_content, and provenance populated.
        predicates: List of TypePredicate objects to evaluate.
        target_isoform: The intended splice isoform ID. Required for the
            SpliceCorrect predicate. If None and SpliceCorrect is in
            predicates, raises MissingPrerequisiteError.

    Returns:
        List of TypeCheckResult objects, one per predicate, in the same
        order as the input predicates.

    Raises:
        PredicateError: If a predicate is malformed or references an
            unknown property.
        MissingPrerequisiteError: If a predicate requires data not present
            in the IR record.
    """
    ...
```

---

### 3.6 IF-06: Optimizer Interface

**Component:** COMP-06 — Optimizer
**Interface ID:** IF-06
**Description:** Solves a Constraint Satisfaction Problem (CSP) to find synonymous codon assignments that satisfy all hard constraints while maximizing a scalar objective (typically CAI). When the CSP is infeasible, it computes a Minimal Unsatisfiable Subset (MUS) to diagnose which constraints conflict. Uses AC-3 for constraint propagation and backtracking search with CAI-ordered domain values.

| Aspect | Specification |
|---|---|
| **Component** | COMP-06 (Optimizer) |
| **Input** | `ir_seq` (IRSeq): The mRNA sequence to optimize. `constraints` (list[Constraint]): Hard constraints that must be satisfied. `objective` (Objective): Scalar objective to maximize (typically CAI). |
| **Output** | `OptimizationResult`: Either `FeasibleResult` (optimized sequence + objective value + satisfied constraints) or `InfeasibleResult` (MUS identifying the minimal conflicting constraint set). |
| **Error Conditions** | `SolverTimeoutError`: The CSP solver exceeded the maximum allowed time. `ConstraintDefinitionError`: A constraint is malformed or references an unknown property. `VariableDomainError`: A codon position has an empty domain after constraint propagation (equivalent to INFEASIBLE). |
| **Determinism** | **Fully deterministic.** Given identical IR-Seq, constraints, and objective, the optimizer produces an identical result (either the same feasible assignment or the same INFEASIBLE report with the same MUS). Domain values are ordered by CAI (descending), ensuring deterministic search order. |
| **Performance** | Time: O(n × d^n) worst case (backtracking), but typically O(n × d) with effective constraint propagation. ≤ 60 seconds for n ≤ 10,000 codon positions. Space: O(n × d) for the CSP variable domains. ≤ 1 MB. |

#### Data Structures

```python
from dataclasses import dataclass, field
from typing import Any

@dataclass(frozen=True)
class Constraint:
    """A hard constraint for the CSP optimizer."""
    name: str                    # Constraint name (e.g., "NoCrypticSplice", "CAI_geq_0.8")
    constraint_type: str         # Type category (e.g., "splicing", "gc_content", "restriction_site")
    parameters: dict[str, Any]   # Constraint parameters (e.g., {"threshold": 0.8})


@dataclass(frozen=True)
class Objective:
    """Scalar objective to maximize during optimization."""
    name: str                    # Objective name (e.g., "CAI", "GC_content")
    direction: str               # "maximize" or "minimize"


@dataclass(frozen=True)
class CodonAssignment:
    """An assignment of a specific codon to a position."""
    position: int       # 0-based codon position
    original_codon: str # Original codon in the input sequence
    assigned_codon: str # Optimized codon assignment
    amino_acid: str     # Amino acid encoded (must be the same for synonymous substitution)


@dataclass(frozen=True)
class FeasibleResult:
    """Result when the CSP has a feasible solution."""
    assignments: list[CodonAssignment]   # Per-position codon assignments
    optimized_sequence: str              # Full optimized nucleotide sequence
    objective_value: float               # Value of the objective function
    constraints_satisfied: list[str]     # Names of all satisfied constraints
    cai: float                           # Codon Adaptation Index of the optimized sequence
    gc_content: float                    # GC content of the optimized sequence


@dataclass(frozen=True)
class InfeasibleResult:
    """Result when the CSP has no feasible solution."""
    mus: list[Constraint]                # Minimal Unsatisfiable Subset of constraints
    mus_description: str                 # Human-readable explanation of the conflict
    all_constraints: list[Constraint]    # Full constraint set (for reference)


@dataclass(frozen=True)
class OptimizationResult:
    """Union type: either a feasible result or an infeasible result."""
    feasible: FeasibleResult | None = None     # Populated if solution found
    infeasible: InfeasibleResult | None = None  # Populated if no solution exists
    is_feasible: bool = False                    # Convenience flag

    def __post_init__(self):
        assert (self.feasible is None) != (self.infeasible is None), \
            "Exactly one of feasible or infeasible must be populated"
```

#### Python Signature

```python
from biocompiler.ir.proto.ir_seq_v1_pb2 import IRSeq

def optimize(
    ir_seq: IRSeq,
    constraints: list[Constraint],
    objective: Objective = Objective(name="CAI", direction="maximize"),
    max_solver_time: float = 60.0,
) -> OptimizationResult:
    """Solve a constraint satisfaction problem for synonymous codon optimization.

    Searches for synonymous codon assignments that satisfy all hard constraints
    while maximizing the specified objective (typically CAI). Uses AC-3 constraint
    propagation followed by backtracking search with CAI-ordered domain values.

    When no feasible assignment exists, computes a Minimal Unsatisfiable Subset
    (MUS) — the smallest subset of constraints that remains unsatisfiable — to
    diagnose the conflict.

    Args:
        ir_seq: IR-Seq record to optimize. Must have sequence, reading_frame,
            and codon assignments populated.
        constraints: List of hard constraints that must be satisfied.
        objective: Scalar objective to maximize or minimize.
            Defaults to maximizing CAI.
        max_solver_time: Maximum time for the CSP solver, in seconds.
            Defaults to 60.0.

    Returns:
        OptimizationResult containing either a FeasibleResult (with optimized
        sequence, codon assignments, and objective value) or an InfeasibleResult
        (with MUS identifying the minimal conflicting constraint set).

    Raises:
        SolverTimeoutError: If the solver exceeds max_solver_time.
        ConstraintDefinitionError: If a constraint is malformed.
    """
    ...
```

---

### 3.7 IF-07: Certificate Generator Interface

**Component:** COMP-07 — Certificate Generator
**Interface ID:** IF-07
**Description:** Produces machine-checkable guarantee certificates in JSON format for designs passing all type checks. Certificates include the verified sequence, each type predicate with its verdict and derivation trace, the CSP constraint set and assignment, and provenance metadata. Circuit certificates additionally include individual gene certificates and composition check results. Certificates are independently verifiable by a separate checker program without access to the BioCompiler pipeline.

| Aspect | Specification |
|---|---|
| **Component** | COMP-07 (Certificate Generator) |
| **Input** | `ir_seq` (IRSeq): The verified sequence. `type_results` (list[TypeCheckResult]): Type check results for all predicates. `optimization_result` (OptimizationResult): Optimization result (feasible or infeasible). `composition_results` (list[CompositionCheck] \| None): Composition check results for circuit certificates. |
| **Output** | `dict`: JSON-serializable guarantee certificate conforming to the Certificate JSON Schema (§3.7.1). |
| **Error Conditions** | `CertificateError`: The certificate cannot be generated due to invalid input (e.g., a type result has no derivation trace for a PASS verdict). `VerificationError`: The generated certificate fails self-verification (internal consistency check). |
| **Determinism** | **Fully deterministic.** Given identical inputs, the generator produces an identical JSON certificate. Timestamps are derived from the pipeline run start time (recorded in provenance), not from the current clock. |
| **Performance** | Time: O(n) for certificate construction. ≤ 100 ms for n ≤ 10,000. Space: O(n) for the certificate JSON. ≤ 5 MB. |

#### Certificate JSON Schema

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://biocompiler.org/schemas/certificate/v1.0.0",
  "title": "BioCompiler Guarantee Certificate",
  "description": "Machine-checkable guarantee certificate for a verified gene design or circuit design. Independently verifiable without access to the BioCompiler pipeline.",
  "type": "object",
  "required": ["version", "design_id", "sequence", "types", "optimization", "provenance"],
  "properties": {
    "version": {
      "type": "string",
      "pattern": "^[0-9]+\\.[0-9]+\\.[0-9]+$",
      "description": "Certificate schema version (semantic versioning)."
    },
    "design_id": {
      "type": "string",
      "pattern": "^[a-f0-9]{64}$",
      "description": "SHA-256 hash of the verified sequence. Provides cryptographic integrity."
    },
    "sequence": {
      "type": "string",
      "pattern": "^[ACGTUNacgtun]+$",
      "description": "The verified nucleotide sequence."
    },
    "types": {
      "type": "array",
      "items": {
        "type": "object",
        "required": ["predicate", "verdict"],
        "properties": {
          "predicate": {
            "type": "string",
            "description": "Name of the type predicate (e.g., 'SpliceCorrect', 'NoCrypticSplice')."
          },
          "verdict": {
            "type": "string",
            "enum": ["PASS", "FAIL", "UNCERTAIN"],
            "description": "Three-valued verdict."
          },
          "derivation": {
            "type": "array",
            "items": {
              "type": "object",
              "required": ["step", "rule", "evidence"],
              "properties": {
                "step": {
                  "type": "integer",
                  "minimum": 1,
                  "description": "Step number in the derivation trace (1-indexed)."
                },
                "rule": {
                  "type": "string",
                  "description": "Rule applied at this step."
                },
                "evidence": {
                  "type": "string",
                  "description": "Evidence supporting this step."
                }
              }
            },
            "description": "Derivation trace for PASS verdicts."
          },
          "violation": {
            "type": "object",
            "properties": {
              "predicate": {
                "type": "string",
                "description": "Name of the violated predicate."
              },
              "position": {
                "type": ["integer", "null"],
                "minimum": 0,
                "description": "Position of the violation (0-based), if applicable."
              },
              "rule": {
                "type": "string",
                "description": "Rule that was violated."
              },
              "evidence": {
                "type": "string",
                "description": "Evidence for the violation."
              }
            },
            "description": "Violation details for FAIL verdicts."
          },
          "knowledge_gap": {
            "type": "object",
            "properties": {
              "predicate": {
                "type": "string",
                "description": "Name of the predicate with the knowledge gap."
              },
              "description": {
                "type": "string",
                "description": "Description of what is needed to resolve."
              },
              "required_data": {
                "type": "string",
                "description": "Description of the data or knowledge required."
              }
            },
            "description": "Knowledge gap for UNCERTAIN verdicts."
          }
        }
      },
      "description": "Type predicate evaluation results."
    },
    "optimization": {
      "type": "object",
      "required": ["objective", "value", "constraints_satisfied"],
      "properties": {
        "objective": {
          "type": "string",
          "description": "Name of the optimization objective (e.g., 'CAI')."
        },
        "value": {
          "type": "number",
          "description": "Objective function value achieved."
        },
        "constraints_satisfied": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "Names of all constraints satisfied by the optimized sequence."
        },
        "codon_assignments": {
          "type": "array",
          "items": {
            "type": "object",
            "required": ["position", "original_codon", "assigned_codon", "amino_acid"],
            "properties": {
              "position": {
                "type": "integer",
                "minimum": 0
              },
              "original_codon": {
                "type": "string",
                "pattern": "^[ACGTU]{3}$"
              },
              "assigned_codon": {
                "type": "string",
                "pattern": "^[ACGTU]{3}$"
              },
              "amino_acid": {
                "type": "string",
                "minLength": 1,
                "maxLength": 1
              }
            }
          },
          "description": "Per-position codon assignments (only for positions that were changed)."
        },
        "mus": {
          "type": "array",
          "items": {
            "type": "object",
            "required": ["name", "constraint_type", "parameters"],
            "properties": {
              "name": {
                "type": "string"
              },
              "constraint_type": {
                "type": "string"
              },
              "parameters": {
                "type": "object"
              }
            }
          },
          "description": "Minimal Unsatisfiable Subset (only present if optimization is infeasible)."
        }
      },
      "description": "Optimization result details."
    },
    "composition": {
      "type": "array",
      "items": {
        "type": "object",
        "required": ["check_type", "genes_involved", "verdict", "evidence"],
        "properties": {
          "check_type": {
            "type": "string",
            "enum": ["promoter_conflict", "resource_competition", "splicing_interference", "rna_interaction", "metabolic_burden"],
            "description": "Type of composition check."
          },
          "genes_involved": {
            "type": "array",
            "items": {
              "type": "string"
            },
            "description": "IDs of genes involved in this check."
          },
          "verdict": {
            "type": "string",
            "enum": ["PASS", "FAIL", "UNCERTAIN"],
            "description": "Composition check verdict."
          },
          "evidence": {
            "type": "string",
            "description": "Evidence for the verdict."
          }
        }
      },
      "description": "Composition check results (only present for circuit certificates)."
    },
    "gene_certificates": {
      "type": "array",
      "items": {
        "type": "string",
        "pattern": "^[a-f0-9]{64}$"
      },
      "description": "SHA-256 hashes of per-gene certificates (only present for circuit certificates)."
    },
    "provenance": {
      "type": "object",
      "required": ["tool", "version", "timestamp", "input_hash"],
      "properties": {
        "tool": {
          "type": "string",
          "description": "Name of the producing tool (e.g., 'BioCompiler')."
        },
        "version": {
          "type": "string",
          "pattern": "^[0-9]+\\.[0-9]+\\.[0-9]+$",
          "description": "Version of the producing tool."
        },
        "timestamp": {
          "type": "string",
          "format": "date-time",
          "description": "ISO 8601 timestamp of certificate generation."
        },
        "input_hash": {
          "type": "string",
          "pattern": "^[a-f0-9]{64}$",
          "description": "SHA-256 hash of the original input specification."
        },
        "parameters": {
          "type": "object",
          "additionalProperties": {
            "type": "string"
          },
          "description": "Key-value pairs of pipeline parameters used."
        },
        "run_id": {
          "type": "string",
          "format": "uuid",
          "description": "UUID of the pipeline run."
        }
      },
      "description": "Provenance metadata for the certificate."
    }
  }
}
```

#### Python Signature

```python
import json
from biocompiler.ir.proto.ir_seq_v1_pb2 import IRSeq

def generate_certificate(
    ir_seq: IRSeq,
    type_results: list[TypeCheckResult],
    optimization_result: OptimizationResult,
    composition_results: list[CompositionCheck] | None = None,
    gene_certificates: list[str] | None = None,
) -> dict:
    """Generate a machine-checkable guarantee certificate in JSON format.

    The certificate encodes the verified sequence, each type predicate with its
    verdict and derivation trace (for PASS) or violation (for FAIL) or knowledge
    gap (for UNCERTAIN), the CSP constraint set and assignment, and provenance
    metadata. Circuit certificates additionally include individual gene certificate
    hashes and composition check results.

    The certificate is independently verifiable by the standalone certificate
    verifier (biocompiler certificate verify) without access to the BioCompiler
    pipeline. The verifier re-evaluates each predicate against the sequence and
    checks the derivation trace.

    Args:
        ir_seq: The verified IR-Seq sequence.
        type_results: Type check results for all evaluated predicates.
        optimization_result: CSP optimization result (feasible or infeasible).
        composition_results: Composition check results. Required for circuit
            certificates; None for single-gene certificates.
        gene_certificates: SHA-256 hashes of per-gene certificates. Required
            for circuit certificates; None for single-gene certificates.

    Returns:
        JSON-serializable dictionary conforming to the Certificate JSON Schema.
        Can be written to a file with json.dump().

    Raises:
        CertificateError: If the certificate cannot be generated due to
            invalid input (e.g., a PASS verdict with no derivation trace).
        VerificationError: If the generated certificate fails self-verification
            (internal consistency check).
    """
    ...
```

---

### 3.8 IF-08: Compositional Verifier Interface

**Component:** COMP-08 — Compositional Verifier
**Interface ID:** IF-08
**Description:** Verifies cross-component constraints in multi-gene circuits, including promoter conflict detection, resource competition analysis, splicing interference detection, and RNA-RNA interaction screening. Composes per-gene certificates into circuit-level proofs without re-running the pipeline. Produces three-valued verdicts composed via the defined three-valued logic algebra.

| Aspect | Specification |
|---|---|
| **Component** | COMP-08 (Compositional Verifier) |
| **Input** | `ir_circuit` (IRCircuit): Circuit topology with gene nodes and interaction edges. `gene_certificates` (list[dict]): Per-gene guarantee certificates. |
| **Output** | `list[CompositionCheck]`: Results for each cross-gene constraint check, with verdict, evidence, and interacting genes. |
| **Error Conditions** | `CertificateExpiredError`: A per-gene certificate has expired or is invalid. `CircuitCycleError`: The circuit graph contains a cycle (violates INV-CIR-01). `MissingCertificateError`: A gene node references a certificate that is not provided. |
| **Determinism** | **Fully deterministic.** Given identical circuit topology and certificates, the verifier produces identical composition check results. |
| **Performance** | Time: O(g² × n) where g = number of genes, n = average gene length. ≤ 30 seconds for g ≤ 10 with n ≤ 10,000. Space: O(g² + g × n) for interaction checks. ≤ 50 MB. |

#### Data Structures

```python
from dataclasses import dataclass
from biocompiler.type_system import Verdict

@dataclass(frozen=True)
class CompositionCheck:
    """Result of a single composition check between genes in a circuit."""
    check_type: str                    # Type: "promoter_conflict", "resource_competition",
                                       #       "splicing_interference", "rna_interaction",
                                       #       "metabolic_burden"
    genes_involved: list[str]          # IDs of genes involved in this check
    verdict: Verdict                   # Three-valued verdict (PASS/FAIL/UNCERTAIN)
    evidence: str                      # Evidence description for the verdict
    severity: float = 0.0              # Severity score [0.0, 1.0] (only for FAIL/UNCERTAIN)
```

#### Python Signature

```python
from biocompiler.ir.proto.ir_circuit_v1_pb2 import IRCircuit

def verify_composition(
    ir_circuit: IRCircuit,
    gene_certificates: list[dict],
) -> list[CompositionCheck]:
    """Verify cross-component constraints in a multi-gene circuit.

    Performs four composition checks:
    1. Promoter conflict — detects unintentional TF-mediated regulation between genes.
    2. Resource competition — checks that total ribosome demand does not exceed
       estimated cellular ribosome capacity.
    3. Splicing interference — detects cryptic splice sites in one gene's transcript
       that could interfere with splicing of another.
    4. RNA-RNA interaction — screens for complementary transcript regions that
       could form dsRNA triggering silencing.

    Each check produces a three-valued verdict (PASS/FAIL/UNCERTAIN) with evidence.
    The overall circuit verdict is composed from individual check verdicts using
    three-valued AND logic.

    Args:
        ir_circuit: IRCircuit record with gene nodes and interaction edges.
            Must have circuit_id, organism, cell_type, genes, and interactions
            populated.
        gene_certificates: List of per-gene guarantee certificates (as dicts
            conforming to the Certificate JSON Schema). Each certificate must
            have a valid design_id matching the corresponding gene node's
            ir_seq_id content hash.

    Returns:
        List of CompositionCheck results, one per cross-gene constraint
        evaluated. The list may contain multiple checks of the same type
        (e.g., promoter_conflict for different gene pairs).

    Raises:
        CertificateExpiredError: If a per-gene certificate is invalid or expired.
        CircuitCycleError: If the circuit graph contains a cycle.
        MissingCertificateError: If a gene node has no matching certificate.
    """
    ...
```

---

### 3.9 IF-09: Mutation Explorer Interface

**Component:** COMP-09 — Mutation Explorer
**Interface ID:** IF-09
**Description:** Decomposes the mutation space of a gene into categories based on splicing grammar nonterminals (intra-exonic, splice site, regulatory element), enumerates legal multi-mutation combinations, exploits independence across exons, and reports constraint conflicts (mutations that are individually legal but jointly violate the splicing grammar).

| Aspect | Specification |
|---|---|
| **Component** | COMP-09 (Mutation Explorer) |
| **Input** | `ir_seq` (IRSeq): The mRNA sequence to explore mutations for. `max_mutations` (int): Maximum number of simultaneous mutations. `mutation_categories` (list[str]): Categories to explore ("intra_exonic", "splice_site", "regulatory"). `constraints` (list[Constraint]): Constraints that legal mutations must satisfy. |
| **Output** | `MutationReport`: Decomposed mutation space, legal combinations, constraint conflicts, and per-nonterminal mutation counts. |
| **Error Conditions** | `EmptyMutationSpaceError`: No legal mutations exist in any category. `ConstraintConflictError`: All mutations in a category are blocked by constraints (reported as part of the MutationReport, not as an exception). |
| **Determinism** | **Fully deterministic.** Given identical IR-Seq and parameters, the explorer produces an identical MutationReport. Mutation enumeration order is deterministic (sorted by position, then by lexicographic codon order). |
| **Performance** | Time: O(c × k^m) where c = number of nonterminals, k = mutations per nonterminal, m = max_mutations. Exponential in m; ≤ 300 seconds for m ≤ 5 with typical gene length. Space: O(r) where r = number of legal combinations found. ≤ 500 MB. |

#### Data Structures

```python
from dataclasses import dataclass, field

@dataclass(frozen=True)
class MutationCandidate:
    """A single mutation candidate at a specific position."""
    position: int           # 0-based position of the mutation
    original_nucleotide: str  # Original nucleotide
    mutated_nucleotide: str   # Mutated nucleotide
    category: str             # Mutation category ("intra_exonic", "splice_site", "regulatory")
    codon_original: str       # Original codon (3-letter)
    codon_mutated: str        # Mutated codon (3-letter)
    amino_acid_original: str  # Original amino acid (single letter)
    amino_acid_mutated: str   # Mutated amino acid (single letter, same for synonymous)
    is_synonymous: bool       # Whether the mutation is synonymous

@dataclass(frozen=True)
class MutationConflict:
    """A conflict between individually legal but jointly illegal mutations."""
    mutations: list[MutationCandidate]  # The conflicting mutations
    violated_constraint: str             # The constraint violated by the combination
    description: str                     # Human-readable conflict description

@dataclass(frozen=True)
class NonterminalDecomposition:
    """Decomposition of the mutation space by grammar nonterminal."""
    category: str                         # "intra_exonic", "splice_site", or "regulatory"
    region: tuple[int, int]               # (start, end) of the region
    candidates: list[MutationCandidate]   # Legal mutations in this region
    is_independent: bool                  # Whether this region is independent of others

@dataclass(frozen=True)
class MutationReport:
    """Complete report from mutation exploration."""
    nonterminals: list[NonterminalDecomposition]  # Per-region decomposition
    legal_combinations: int                         # Total number of legal multi-mutation combinations
    conflicts: list[MutationConflict]              # Constraint conflicts between mutations
    total_candidates: int                           # Total individual mutation candidates
    independent_groups: int                         # Number of independent mutation groups
    explored_combinations: int                      # Number of combinations actually explored
```

#### Python Signature

```python
from biocompiler.ir.proto.ir_seq_v1_pb2 import IRSeq

def explore_mutations(
    ir_seq: IRSeq,
    max_mutations: int = 3,
    mutation_categories: list[str] | None = None,
    constraints: list[Constraint] | None = None,
    max_explore_time: float = 300.0,
) -> MutationReport:
    """Enumerate grammar-guided legal mutation combinations for a gene design.

    Decomposes the gene's mutation space by splicing grammar nonterminals
    (intra-exonic, splice site, regulatory element), enumerates legal
    multi-mutation combinations, exploits independence across exons for
    efficient enumeration, and detects constraint conflicts between
    individually legal but jointly illegal mutations.

    Args:
        ir_seq: IR-Seq record to explore mutations for. Must have sequence,
            exon_boundaries, splice_sites, and regulatory_elements populated.
        max_mutations: Maximum number of simultaneous mutations in a
            combination. Defaults to 3.
        mutation_categories: Categories to explore. If None, all categories
            are explored. Supported values: "intra_exonic", "splice_site",
            "regulatory".
        constraints: Constraints that legal mutations must satisfy. If None,
            only the splicing grammar constraint is applied.
        max_explore_time: Maximum time for exploration, in seconds.
            Defaults to 300.0 (5 minutes).

    Returns:
        MutationReport with per-region decomposition, legal combination count,
        constraint conflicts, and exploration statistics.

    Raises:
        EmptyMutationSpaceError: If no legal mutations exist in any category.
    """
    ...
```

---

### 3.10 IF-10: ORF Analyzer Interface

**Component:** COMP-10 — ORF Analyzer
**Interface ID:** IF-10
**Description:** Computes shared constraint sets for overlapping reading frames in viral genomes and compact genomes, classifies nucleotide positions as high-coupling (affecting multiple proteins) or low-coupling (affecting one protein), and detects constraint conflicts between frames where the optimization target for one frame contradicts the target for another.

| Aspect | Specification |
|---|---|
| **Component** | COMP-10 (ORF Analyzer) |
| **Input** | `sequence` (str): Nucleotide sequence with multiple annotated reading frames. `reading_frames` (list[ReadingFrameSpec]): Specifications for each reading frame. |
| **Output** | `ORFAnalysisReport`: Shared constraint set, per-position coupling classification, and inter-frame constraint conflicts. |
| **Error Conditions** | `InvalidFrameSpecError`: A reading frame specification is out of bounds or has an invalid frame offset. `NoOverlappingORFsError`: The specified reading frames do not overlap (analysis is trivial). `SequenceTooShortError`: The sequence is too short to contain the specified reading frames. |
| **Determinism** | **Fully deterministic.** Given identical sequence and reading frame specifications, the analyzer produces an identical ORFAnalysisReport. |
| **Performance** | Time: O(f × n) where f = number of reading frames, n = sequence length. ≤ 10 seconds for f ≤ 6 with n ≤ 30,000. Space: O(f × n) for per-position frame membership. ≤ 50 MB. |

#### Data Structures

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class ReadingFrameSpec:
    """Specification for a reading frame to analyze."""
    name: str           # Frame name (e.g., "ORF1a", "ORF1b", "S", "E", "M", "N")
    start: int          # Start position (0-based, inclusive)
    end: int            # End position (0-based, exclusive)
    frame: int          # Reading frame offset (0, 1, or 2)
    strand: str         # "forward" or "reverse"

@dataclass(frozen=True)
class FrameConflict:
    """A constraint conflict between overlapping reading frames at a position."""
    position: int                     # 0-based position of the conflict
    frames: list[str]                 # Names of the conflicting frames
    optimal_codons: dict[str, str]    # Per-frame optimal codon at this position
    conflict_type: str                # "codon_mismatch" or "constraint_violation"
    description: str                  # Human-readable conflict description

@dataclass(frozen=True)
class PositionCoupling:
    """Coupling classification for a nucleotide position."""
    position: int               # 0-based position
    affecting_frames: list[str] # Names of frames affected by a mutation at this position
    coupling: str               # "high" (≥2 frames) or "low" (1 frame)

@dataclass(frozen=True)
class ORFAnalysisReport:
    """Complete report from overlapping reading frame analysis."""
    shared_constraint_positions: list[int]            # Positions affecting ≥2 frames
    position_couplings: list[PositionCoupling]         # Per-position coupling classification
    conflicts: list[FrameConflict]                     # Inter-frame constraint conflicts
    num_high_coupling: int                              # Number of high-coupling positions
    num_low_coupling: int                               # Number of low-coupling positions
    num_conflicts: int                                  # Number of constraint conflicts
    reading_frames: list[ReadingFrameSpec]              # Input reading frame specifications
    sequence_length: int                                # Total sequence length
```

#### Python Signature

```python
def analyze_overlapping_orfs(
    sequence: str,
    reading_frames: list[ReadingFrameSpec],
) -> ORFAnalysisReport:
    """Analyze overlapping reading frames in a nucleotide sequence.

    Computes the shared constraint set (positions where mutations affect multiple
    proteins simultaneously), classifies positions as high-coupling (≥2 frames
    affected) or low-coupling (1 frame affected), and detects constraint conflicts
    between frames where the optimization target for one frame contradicts the
    target for another.

    Args:
        sequence: Nucleotide sequence in IUPAC notation (A/C/G/T/U/N).
        reading_frames: List of ReadingFrameSpec objects, one per reading frame
            to analyze. Each spec defines the frame name, start, end, offset,
            and strand.

    Returns:
        ORFAnalysisReport with shared constraint positions, per-position
        coupling classification, and inter-frame constraint conflicts.

    Raises:
        InvalidFrameSpecError: If a frame spec is out of bounds or has an
            invalid frame offset.
        SequenceTooShortError: If the sequence is too short for the specified
            reading frames.
    """
    ...
```

---

## 4. FFI Adapter Contract

The FFI Adapter Contract defines the interface that all Foreign Function Interface (FFI) adapters must implement. Adapters are the mechanism by which the BioCompiler system invokes external tools (protein structure predictors, PTM prediction tools) through a defined boundary. The FFI boundary is strictly enforced: the deterministic core pipeline treats every FFI invocation as a non-deterministic black-box operation, and no guarantee certificate depends on FFI output for its core validity. FFI stages enrich the IR with optional SLOT fields that the type system may reference but does not require for its fundamental correctness guarantees.

### 4.1 Abstract Adapter Interface

Every FFI adapter MUST implement the following Python Abstract Base Class (ABC):

```python
from abc import ABC, abstractmethod
from biocompiler.ir.proto.ir_peptide_v1_pb2 import IRPeptide
from biocompiler.ir.proto.ir_structure_v1_pb2 import IRStructure

class FFIAdapter(ABC):
    """Abstract base class for all FFI adapters.

    Every adapter must implement name(), slot_fields(), invoke(),
    and validate_output(). The FFI Manager (COMP-04) uses these
    methods to orchestrate external tool invocations.

    Adapters are registered in the adapter registry by name. The
    FFI Manager selects adapters by name and delegates formatting,
    invocation, parsing, and validation to the adapter.
    """

    @abstractmethod
    def name(self) -> str:
        """Return the unique name of this adapter.

        The name is used to register and select the adapter in the
        FFI Manager's adapter registry. Must be lowercase with
        underscores (e.g., "alphafold2", "netphos").

        Returns:
            Adapter name string.
        """
        ...

    @abstractmethod
    def slot_fields(self) -> list[str]:
        """Return the list of IR SLOT fields that this adapter fills.

        SLOT fields are initially empty in the IR record and are
        filled by the adapter's output. This method allows the FFI
        Manager to verify that all expected SLOT fields are populated
        after a successful invocation.

        Returns:
            List of SLOT field names (e.g., ["atom_coordinates",
            "residue_confidences", "pae_matrix"]).
        """
        ...

    @abstractmethod
    def invoke(self, input_data: bytes, timeout: float) -> bytes:
        """Invoke the external tool with the formatted input data.

        This method is responsible for:
        1. Launching the external tool as a subprocess or API call.
        2. Passing the input data.
        3. Waiting for the tool to complete, up to the specified timeout.
        4. Returning the raw output bytes from the tool.

        The FFI Manager handles timeout enforcement. If the tool
        exceeds the timeout, the manager raises ExternalToolError.

        Args:
            input_data: Formatted input data for the external tool,
                as produced by format_input().
            timeout: Maximum wall-clock time for the external tool,
                in seconds.

        Returns:
            Raw output bytes from the external tool.

        Raises:
            ExternalToolError: If the tool times out or returns a
                non-zero exit code.
        """
        ...

    @abstractmethod
    def validate_output(self, output) -> bool:
        """Validate that the parsed output satisfies IR invariants.

        This method checks that the parsed output from the external
        tool satisfies the relevant IR invariants before the SLOT
        fields are filled. If validation fails, the FFI Manager
        raises OutputValidationError.

        Specific invariants to check depend on the adapter type:
        - Folding adapters: INV-STR-01 (confidence scores in [0.0, 1.0])
          and INV-STR-02 (PTM positions reference valid residues).
        - PTM adapters: INV-PEP-03 (no duplicate isoform references)
          and residue positions within bounds.

        Args:
            output: Parsed output dictionary from parse_output().

        Returns:
            True if the output passes all invariant checks; False otherwise.
        """
        ...

    # The following methods are concrete helpers that adapters MAY override
    # for custom formatting and parsing.

    def format_input(self, ir_peptide: IRPeptide, config: dict) -> bytes:
        """Format the IR-Peptide as input for the external tool.

        Default implementation produces a FASTA-formatted string with
        the amino acid sequence. Adapters may override for custom formats
        (e.g., A3M for AlphaFold, JSON for REST APIs).

        Args:
            ir_peptide: IR-Peptide record with amino acid sequence.
            config: Adapter-specific configuration dictionary.

        Returns:
            Formatted input data as bytes.
        """
        fasta = f">{ir_peptide.peptide_id or 'sequence'}\n{ir_peptide.amino_acid_sequence}\n"
        return fasta.encode("utf-8")

    def parse_output(self, raw_output: bytes) -> dict:
        """Parse the raw output bytes from the external tool.

        Default implementation attempts JSON parsing. Adapters MUST
        override this for tools that produce binary or custom-format
        output (e.g., PDB/mmCIF for folding tools).

        Args:
            raw_output: Raw output bytes from the external tool.

        Returns:
            Parsed output as a dictionary.

        Raises:
            OutputParseError: If the output cannot be parsed.
        """
        import json
        try:
            return json.loads(raw_output.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            raise OutputParseError(f"Failed to parse output as JSON: {e}")

    def fill_slots(self, ir_record, parsed: dict) -> None:
        """Fill SLOT fields in the IR record from parsed output.

        This method mutates the IR record in place, filling the SLOT
        fields that this adapter is responsible for (as declared by
        slot_fields()).

        Args:
            ir_record: IR record (IRPeptide or IRStructure) to fill.
            parsed: Parsed output dictionary from parse_output().
        """
        ...
```

### 4.2 Folding Adapter Contract (AlphaFold / ColabFold)

| Aspect | Specification |
|---|---|
| **Adapter Names** | `"alphafold2"`, `"alphafold3"`, `"colabfold"`, `"rosettafold"` |
| **SLOT Fields Filled** | `atom_coordinates`, `residue_confidences`, `pae_matrix`, `validation_flags`, `energy_estimates`, `secondary_structure`, `solvent_accessibility`, `disorder_regions` |
| **Input Format** | FASTA (amino acid sequence) for single-chain; A3M (MSA) for multi-chain. Optionally: MSA depth specification, template structures. |
| **Output Format** | PDB or mmCIF file (3D coordinates) + JSON confidence metrics (pLDDT per residue, PAE matrix). AlphaFold3 may output CIF directly. |
| **Timeout** | Default 3600 seconds (1 hour). Configurable via adapter_config. GPU-accelerated runs typically complete in 10–60 minutes for sequences ≤ 1000 residues. |
| **Validation Rules** | (1) INV-STR-01: All confidence scores in [0.0, 1.0] (pLDDT normalized from [0, 100] scale). (2) INV-STR-02: All residue positions reference valid indices in [0, peptide_length). (3) All energy estimates are finite IEEE 754 doubles (INV-STR-03). (4) Number of residues in output matches `peptide_length`. (5) No NaN or Inf values in coordinates or confidence scores. |
| **Fallback Behavior** | On timeout: Return IR-Structure with `mean_plddt = 0.0` and all SLOT fields empty; type system issues UNCERTAIN verdict. On parse error: Same fallback with diagnostic. On validation failure: Attempt partial fill of valid residues; flag invalid residues with `ValidationFlag`. |
| **Provenance Recorded** | `tool_name`, `tool_version`, `timestamp`, `run_id`, `input_hash`, `parameters` (including `model_preset`, `msa_depth`, `use_gpu`), `folding_tool`, `folding_version`, `model_identifier`, `model_rank`, `msa_depth` |
| **Non-determinism** | GPU floating-point operations may produce slightly different coordinates across runs. MSA construction may vary with database updates. These are captured as UNCERTAIN verdicts in the type system. |

### 4.3 PTM Adapter Contract (NetPhos / PhosphoSitePlus)

| Aspect | Specification |
|---|---|
| **Adapter Names** | `"netphos"`, `"phosphositeplus"`, `"dbptm"`, `"musitedeep"` |
| **SLOT Fields Filled** | `ptm_sites` (in IR-Peptide) |
| **Input Format** | FASTA (amino acid sequence). Some adapters accept additional context (organism, tissue type). |
| **Output Format** | JSON with list of PTM site predictions, each containing: residue position, residue type, PTM type, confidence score, predictor name, evidence type. NetPhos outputs a custom text format that is parsed. |
| **Timeout** | Default 300 seconds (5 minutes). Configurable via adapter_config. Typically completes in ≤ 60 seconds for sequences ≤ 2000 residues. |
| **Validation Rules** | (1) INV-STR-02: All PTM site positions are in [0, peptide_length). (2) All confidence scores are in [0.0, 1.0]. (3) All residue types match the amino acid at the specified position. (4) No duplicate PTM site entries (same position + same PTM type). |
| **Fallback Behavior** | On timeout: Return IR-Peptide with `ptm_sites` empty; type system issues UNCERTAIN verdict. On parse error: Same fallback with diagnostic. On validation failure: Discard invalid PTM sites; retain valid ones. |
| **Provenance Recorded** | `tool_name`, `tool_version`, `timestamp`, `run_id`, `input_hash`, `parameters` (including `ptm_types`, `threshold`) |
| **Non-determinism** | PTM prediction tools are generally deterministic for the same input and model version. Database-backed adapters (PhosphoSitePlus, dbPTM) may return different results as databases are updated. These are captured as UNCERTAIN verdicts. |

---

## 5. CLI Interface

The BioCompiler command-line interface (CLI) provides human-accessible entry points for all pipeline functions. The CLI is built using Python's `typer` library (which wraps `click`) and follows standard POSIX conventions for command structure, option naming, and exit codes.

### 5.1 Command Structure

The top-level command is `biocompiler`. Subcommands correspond to major pipeline functions:

```
biocompiler <subcommand> [options] [arguments]
```

| Subcommand | Interface ID | Description | Primary Component |
|---|---|---|---|
| `design` | IF-01 + IF-02 + IF-03 + IF-05 + IF-06 + IF-07 | Full pipeline: scan → splice → translate → type check → optimize → certify | COMP-01..07 |
| `verify` | IF-05 | Run type checking on an existing sequence without optimization | COMP-05 |
| `explore` | IF-09 | Explore legal mutation combinations | COMP-09 |
| `analyze-orf` | IF-10 | Analyze overlapping reading frames | COMP-10 |
| `verify-circuit` | IF-08 | Verify a multi-gene circuit | COMP-08 |
| `check-cert` | — | Standalone certificate verification (no pipeline) | COMP-07 (verifier) |

### 5.2 Design Command

The `design` subcommand executes the full single-gene pipeline:

```
biocompiler design [OPTIONS] <input>
```

Where `<input>` is the path to the input sequence file (FASTA format) or `-` for stdin.

| Option | Short | Type | Default | Description |
|---|---|---|---|---|
| `--organism` | `-o` | str | `"Homo sapiens"` | Target organism for codon optimization and splicing grammar. |
| `--cell-type` | `-c` | str | `"HEK293T"` | Target cellular context for splicing. |
| `--cai-threshold` | | float | `0.8` | Minimum Codon Adaptation Index threshold. |
| `--gc-min` | | float | `0.4` | Minimum GC content (as fraction). |
| `--gc-max` | | float | `0.6` | Maximum GC content (as fraction). |
| `--avoid-restriction` | `-r` | list[str] | `[]` | Restriction enzyme names to avoid (e.g., `-r EcoRI -r BamHI`). |
| `--splicing-rules` | | str | `"config/splicing_rules.yaml"` | Path to splicing grammar rules YAML. |
| `--folding-adapter` | | str | `None` | Folding FFI adapter to use. If not specified, no folding is performed. |
| `--ptm-adapter` | | str | `None` | PTM FFI adapter to use. If not specified, no PTM prediction is performed. |
| `--folding-timeout` | | float | `3600.0` | Timeout for folding FFI invocation (seconds). |
| `--ptm-timeout` | | float | `300.0` | Timeout for PTM FFI invocation (seconds). |
| `--max-solver-time` | | float | `60.0` | Maximum time for the CSP solver (seconds). |
| `--output` | | str | `None` | Output file path. If not specified, output is written to stdout. |
| `--certificate-output` | | str | `None` | Certificate output file path. If not specified, certificate is written alongside the output. |
| `--format` | | str | `"json"` | Output format: "json", "fasta", or "genbank". |
| `--verbose` | `-v` | flag | `False` | Enable verbose output with diagnostic details. |
| `--quiet` | `-q` | flag | `False` | Suppress all output except errors. |
| `--version` | `-V` | flag | `False` | Print version and exit. |
| `--help` | `-h` | flag | `False` | Print help and exit. |

**Example usage:**

```bash
biocompiler design gene_input.fasta \
    --organism "Homo sapiens" \
    --cell-type HEK293T \
    --cai-threshold 0.8 \
    --gc-min 0.4 --gc-max 0.6 \
    --avoid-restriction EcoRI BamHI XhoI \
    --folding-adapter alphafold2 \
    --output optimized_gene.fasta \
    --certificate-output certificate.json
```

### 5.3 Exit Codes

The BioCompiler CLI uses the following exit codes to communicate the pipeline result:

| Exit Code | Name | Meaning | Action for User |
|---|---|---|---|
| `0` | **Success** | The pipeline completed successfully. All type predicates returned PASS. A guarantee certificate was generated. | Use the output sequence and certificate. |
| `1` | **Infeasible** | The CSP optimizer found no feasible assignment. The constraint set is unsatisfiable. The MUS (Minimal Unsatisfiable Subset) is reported in the output. | Review the MUS to identify conflicting constraints. Relax one or more constraints and retry. |
| `2` | **Uncertain** | One or more type predicates returned UNCERTAIN. The design cannot be fully verified, but no FAIL verdicts were issued. The certificate includes UNCERTAIN verdicts with knowledge gap specifications. | Review the knowledge gaps. Provide additional data or accept the uncertainty. |
| `10` | **Input Error** | The input is invalid (malformed FASTA, non-IUPAC characters, missing required fields, unknown organism, etc.). | Fix the input and retry. |
| `11` | **Internal Error** | An unexpected internal error occurred (pipeline bug, invariant violation, out of memory). A diagnostic report is written to stderr. | Report the error to the development team with the diagnostic report. |

**Exit code composition rule:** If multiple error conditions apply, the most specific exit code is used. The priority order is: `0` (Success) > `1` (Infeasible) > `2` (Uncertain) > `10` (Input Error) > `11` (Internal Error). For example, if both an input error and an internal error occur, exit code `10` is used (the input error is the root cause).

---

## Appendix A: Glossary

| Term | Definition |
|---|---|
| **CAI** | Codon Adaptation Index — a deterministic metric measuring the similarity of codon usage in a gene to the codon usage of highly expressed genes in a target organism. Values range from 0 to 1, with higher values indicating better adaptation. |
| **CSP** | Constraint Satisfaction Problem — a mathematical problem defined by a set of variables, each with a domain of possible values, and a set of constraints that restrict the values the variables can simultaneously take. The solver finds assignments satisfying all constraints or proves that no such assignment exists. |
| **DFA** | Deterministic Finite Automaton — a finite-state machine that accepts or rejects a given string of symbols by running through a uniquely determined sequence of states. Used in the BioCompiler scanner for motif detection. |
| **FFI** | Foreign Function Interface — a mechanism by which the BioCompiler system invokes external tools (folding predictors, PTM predictors) through a defined adapter contract, treating the external tool as a black box with non-deterministic output. |
| **FST** | Finite-State Transducer — a finite automaton that produces an output string for each input string, implementing a mapping from input sequences to output sequences. The translation engine uses a deterministic FST. |
| **IR** | Intermediate Representation — a typed, schema-enforced data structure (defined in protocol buffers) that mediates data flow between pipeline stages. Each IR level (IR-Seq, IR-Peptide, IR-Structure, IR-Circuit) has defined invariants. |
| **MUS** | Minimal Unsatisfiable Subset — the smallest subset of constraints from an infeasible CSP that remains unsatisfiable when considered in isolation. Used to diagnose infeasibility by pinpointing the exact set of conflicting constraints. |
| **NDFST** | Non-Deterministic Finite-State Transducer — a finite-state transducer that maps each input string to a set of possible output strings (not a probability distribution). The set-valued computation is deterministic: the same input always produces the same set of outputs. Used for splicing grammar parsing. |
| **ORF** | Open Reading Frame — a contiguous sequence of codons from a start codon (AUG) to a stop codon (UAA/UAG/UGA) that potentially encodes a protein. |
| **pLDDT** | Predicted Local Distance Difference Test — a per-residue confidence score (0–100 scale) produced by AlphaFold and related structure prediction tools. Higher scores indicate higher confidence in the predicted structure. |
| **PTM** | Post-Translational Modification — a chemical modification applied to a protein after translation, such as phosphorylation, glycosylation, acetylation, or ubiquitination. Predicted by external tools via the FFI. |
| **SECIS** | Selenocysteine Insertion Sequence — a stem-loop structure in the 3' UTR of selenoprotein mRNAs that directs the ribosome to recode a UGA stop codon as selenocysteine (Sec, U) instead of terminating translation. |
| **SLOT** | A field in the Intermediate Representation (IR) that is initially empty and is filled by a specific FFI adapter. SLOT fields enable the separation of deterministic core computation from non-deterministic external tool output. The type system may reference SLOT fields but does not require them for its fundamental correctness guarantees. |
| **UTR** | Untranslated Region — the portions of mRNA before the start codon (5' UTR) and after the stop codon (3' UTR) that are not translated into protein but contain regulatory elements. |
| **ESE** | Exonic Splicing Enhancer — a short nucleotide motif within an exon that promotes the inclusion of that exon in the mature mRNA by recruiting splicing activator proteins (SR proteins). |
| **ESS** | Exonic Splicing Silencer — a short nucleotide motif within an exon that suppresses the inclusion of that exon in the mature mRNA by recruiting splicing repressor proteins (hnRNPs). |
| **ISE** | Intronic Splicing Enhancer — a short nucleotide motif within an intron that promotes the inclusion of the adjacent exon in the mature mRNA. |
| **ISS** | Intronic Splicing Silencer — a short nucleotide motif within an intron that suppresses the inclusion of the adjacent exon in the mature mRNA. |

---

*End of DOC-04: Interface Control Document (ICD) — Version 12.0.0 — Status: Current*
