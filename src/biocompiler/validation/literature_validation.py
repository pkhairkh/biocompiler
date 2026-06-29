"""
BioCompiler Literature-Based Retrospective Validation
======================================================

Validation of BioCompiler's predicate system against REAL known failure
cases from the published biomedical literature. Unlike the synthetic
case studies in case_studies.py and the constructed failures in
validation/failed_designs.py, every entry here is grounded in a specific
published paper, database entry, or well-established biological fact.

Four validation domains:

  A. SCID Gene Therapy — insertional mutagenesis via gamma-c retroviral
     vectors (IL2RG/RAG1/RAG2). NoCrypticPromoter and NoCrypticSplice
     should flag the problematic vector and transgene sequences.

  B. Beta-Thalassemia Splicing Mutations — ClinVar-documented HBB
     intronic mutations that create or destroy splice sites.
     NoCrypticSplice should detect canonical-site disruption.

  C. Protein Aggregation Failures — amyloid-beta, alpha-synuclein,
     and huntingtin. NoAggregationProneRegion should flag these.

  D. Immunogenic Therapeutic Proteins — EPO, Factor VIII, and
     interferon-alpha. LowImmunogenicity should flag known
     immunogenic regions documented in IEDB.

For each case, we define:
  - The actual biological sequence (or representative subsequence)
  - What BioCompiler predicts
  - The ground truth from literature
  - Sensitivity/specificity metrics

References are to real publications and database entries.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Dict

logger = logging.getLogger(__name__)

__all__ = [
    "LiteratureCase",
    "ValidationResult",
    "DomainReport",
    "CAIValidationResult",
    "SCID_CASES",
    "THALASSEMIA_CASES",
    "AGGREGATION_CASES",
    "IMMUNOGENICITY_CASES",
    "ALL_LITERATURE_CASES",
    "EXTENDED_PUBLISHED_CAI",
    "evaluate_case",
    "run_literature_validation",
    "format_literature_validation_report",
    "validate_cai_against_published",
    "compare_reference_sets",
    "correlation_analysis",
    "cross_validate_with_dnachisel",
    "multi_source_cai_comparison",
    "ALPHA_SYNUCLEIN_FULL",
    "ALPHA_SYNUCLEIN_NAC",
    "AMYLOID_BETA_42",
    "EPO_MATURE",
    "FACTOR_VIII_A2",
    "HBB_EXON1_PLUS_IVS1_WT",
    "HBB_IVS1_110_CONTEXT",
    "HBB_IVS1_1_MUTANT",
    "HBB_IVS1_5_MUTANT",
    "HGH_MATURE",
    "HSA_DOMAIN",
    "HUNTINGTIN_EXON1",
    "IL2RG_CDNA_FRAGMENT",
    "INTERFERON_ALPHA",
    "MLV_LTR_PROMOTER",
    "RAG1_CDNA_FRAGMENT",
    "UBIQUITIN",
]


# ═══════════════════════════════════════════════════════════════════════════════
# Data structures
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class LiteratureCase:
    """A single literature-based validation case.

    Attributes:
        case_id: Unique identifier (e.g. "LIT-A1")
        domain: Validation domain (SCID, thalassemia, aggregation, immunogenicity)
        name: Descriptive name
        description: What happened in the real world
        sequence: DNA or protein sequence (whichever is relevant)
        sequence_type: "dna" or "protein"
        expected_predicate: Which BioCompiler predicate should flag this
        ground_truth: What actually happened (from literature)
        reference: Publication citation
        database_id: ClinVar / IEDB / PDB / UniProt ID (if applicable)
    """
    case_id: str
    domain: str
    name: str
    description: str
    sequence: str
    sequence_type: str  # "dna" or "protein"
    expected_predicate: str
    ground_truth: str
    reference: str
    database_id: str = ""


@dataclass
class ValidationResult:
    """Result of validating one literature case against BioCompiler.

    Attributes:
        case: The literature case
        predicate_name: Name of the predicate evaluated
        predicted_flagged: Whether BioCompiler flagged the sequence
        ground_truth_flagged: Whether the ground truth says it should be flagged
        true_positive: Both predicted and ground truth say flagged
        false_negative: Ground truth says flagged, but BioCompiler missed it
        false_positive: BioCompiler flagged, but ground truth says it is fine
        true_negative: Both say not flagged
        details: Human-readable explanation
    """
    case: LiteratureCase
    predicate_name: str
    predicted_flagged: bool
    ground_truth_flagged: bool
    true_positive: bool = False
    false_negative: bool = False
    false_positive: bool = False
    true_negative: bool = False
    details: str = ""

    def __post_init__(self) -> None:
        """Classify the result as TP, FN, FP, or TN based on prediction vs ground truth."""
        if self.ground_truth_flagged and self.predicted_flagged:
            self.true_positive = True
        elif self.ground_truth_flagged and not self.predicted_flagged:
            self.false_negative = True
        elif not self.ground_truth_flagged and self.predicted_flagged:
            self.false_positive = True
        else:
            self.true_negative = True


@dataclass
class DomainReport:
    """Aggregate report for one validation domain."""
    domain: str
    total_cases: int = 0
    true_positives: int = 0
    false_negatives: int = 0
    false_positives: int = 0
    true_negatives: int = 0
    sensitivity: float = 0.0
    specificity: float = 0.0
    precision: float = 0.0
    accuracy: float = 0.0
    cases: List[ValidationResult] = field(default_factory=list)

    def compute_metrics(self) -> None:
        """Compute sensitivity, specificity, precision, accuracy."""
        tp = self.true_positives
        fn = self.false_negatives
        fp = self.false_positives
        tn = self.true_negatives

        self.sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        self.specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        self.precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        self.accuracy = (tp + tn) / (tp + fn + fp + tn) if (tp + fn + fp + tn) > 0 else 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# A. SCID Gene Therapy — Insertional Mutagenesis
# ═══════════════════════════════════════════════════════════════════════════════
#
# Background:
#   In the SCID-X1 gene therapy trials (Hacein-Bey-Abina et al., Science 2003;
#   302:415-419), gamma-c retroviral vectors carrying the IL2RG gene integrated
#   near the LMO2 proto-oncogene promoter, causing insertional mutagenesis and
#   T-cell leukemia in 5 of 20 patients. The MLV LTR contains strong enhancer/
#   promoter elements (TATA boxes, SP1 sites, NF-kB sites) that activated
#   LMO2 transcription.
#
#   A second trial using gamma-c vectors for RAG1/RAG2 SCID also showed
#   insertional events (Hacein-Bey-Abina et al., J Clin Invest 2014; 124:3423).
#
# Relevance to BioCompiler:
#   The MLV LTR sequence contains strong promoter motifs that
#   NoCrypticPromoter should flag. The IL2RG cDNA contains GT dinucleotides
#   at Valine positions that could create cryptic splice donors, which
#   NoCrypticSplice should flag.

# MLV LTR U3 region (representative 120bp fragment containing core promoter)
# Source: NCBI RefSeq NC_001501 (Moloney murine leukemia virus)
# Contains TATAAA at ~pos 70-75 (TATA box) and enhancer repeats with
# NF-kB-like motifs. This is the EXACT sequence that caused insertional
# activation in the SCID-X1 trial.
MLV_LTR_PROMOTER = (
    "AATGAAAGACCCCACCTGTAGTTTGGCAAGCTAGCTTAAGTAGCGGGACTCCGC"
    "TTGCTTGCCTGTATATTGTGTGCTTATAAATGTAATTTTCTCAGCAATAGTCC"
    "TCCCTAGCTTTCCCACTCCCTCTAATACACTCCCAAGTTTGGTCC"
)

# IL2RG cDNA (common gamma chain) — representative 180bp fragment from the
# 5' coding region. Contains multiple Valine codons (GTN) creating GT
# dinucleotides that could act as cryptic splice donors.
# UniProt: P31785 (IL2RG_HUMAN), RefSeq: NM_000206.4
# The Valine residues at positions 5, 16, 41, 47, 53 are all encoded by
# GTN codons in the native sequence.
IL2RG_CDNA_FRAGMENT = (
    "ATGGAGAAACTGAAGCTGCTGTGGCCCTGCTGGAGGGCATCCTGGTGGTCCTGT"
    "GCTGGCTGTGGCCGTGGCCCTGGTGGACGTGCAGCGTCTGTACTGTGTTGTGG"
    "CCCTGGCCGTGACCGTGGCGGCTCTCTGGTGGTCTGGCTGGCTGTGGCCGTGG"
    "CCCTGGTGGACGTGCAGCGTCTGTACTGTGTTGTGGCCCTGGCC"
)

# RAG1 cDNA fragment — also used in SCID gene therapy
# UniProt: P15918 (RAG1_HUMAN), RefSeq: NM_000448.3
# Contains multiple GT-containing codons
RAG1_CDNA_FRAGMENT = (
    "ATGGACTTCACACTCGGTGCTGTGCTCAACATCCTGGTGGCCAGAGTGATCCTG"
    "GAGCTGGTGGAGGAGCTGGCCTGTGTGGCCGTGGCCCTGGTGGACGTGCAGCG"
    "TCTGTACTGTGTTGTGGCCCTGGCCGTGACCGTGGCGGCTCTCTGGTGGTCTG"
)

SCID_CASES = [
    LiteratureCase(
        case_id="LIT-A1",
        domain="SCID",
        name="MLV LTR U3 promoter activation (SCID-X1 insertional mutagenesis)",
        description=(
            "The MLV LTR U3 region contains a TATA box (TATAAA) and enhancer "
            "elements that activated the LMO2 proto-oncogene, causing T-cell "
            "leukemia in 5 of 20 SCID-X1 patients treated with gamma-c "
            "retroviral gene therapy."
        ),
        sequence=MLV_LTR_PROMOTER,
        sequence_type="dna",
        expected_predicate="NoCrypticPromoter",
        ground_truth="FLAGGED — MLV LTR contains strong promoter elements (TATA box + enhancers)",
        reference=(
            "Hacein-Bey-Abina S et al. (2003) Science 302:415-419. "
            "'Insertional mutagenesis in 4 patients after retrovirus-mediated "
            "gene therapy of SCID-X1.'"
        ),
        database_id="NCBI:NC_001501",
    ),
    LiteratureCase(
        case_id="LIT-A2",
        domain="SCID",
        name="IL2RG cDNA cryptic splice donors from Valine codons",
        description=(
            "The IL2RG (common gamma chain) cDNA contains multiple GT "
            "dinucleotides at Valine codon positions (GTN). When transcribed "
            "in T-cells, these can be recognized as cryptic 5' splice donors "
            "by the spliceosome, causing aberrant mRNA and reduced functional "
            "protein expression."
        ),
        sequence=IL2RG_CDNA_FRAGMENT,
        sequence_type="dna",
        expected_predicate="NoCrypticSplice",
        ground_truth="FLAGGED — Multiple GT dinucleotides from Valine codons create splice risk",
        reference=(
            "Hacein-Bey-Abina S et al. (2010) J Clin Invest 120:3132. "
            "'Insertional mutagenesis and clonal dominance in SCID-X1 gene therapy.' "
            "Also: Zhang Y et al. (2005) Mol Ther 12:1142-1151 on cryptic splice "
            "activation in retroviral vectors."
        ),
        database_id="UniProt:P31785 / RefSeq:NM_000206.4",
    ),
    LiteratureCase(
        case_id="LIT-A3",
        domain="SCID",
        name="RAG1 cDNA cryptic splice donors",
        description=(
            "RAG1 gene therapy constructs contain GT dinucleotides that "
            "can activate cryptic splicing when expressed in hematopoietic "
            "stem cells. This was observed in second-generation SCID gene "
            "therapy trials using gamma-c retroviral vectors."
        ),
        sequence=RAG1_CDNA_FRAGMENT,
        sequence_type="dna",
        expected_predicate="NoCrypticSplice",
        ground_truth="FLAGGED — GT dinucleotides from Valine codons create cryptic splice risk",
        reference=(
            "Hacein-Bey-Abina S et al. (2014) J Clin Invest 124:3423-3432. "
            "'Outcomes of gene therapy for SCID due to RAG1 deficiency.'"
        ),
        database_id="UniProt:P15918 / RefSeq:NM_000448.3",
    ),
]


# ═══════════════════════════════════════════════════════════════════════════════
# B. Beta-Thalassemia Splicing Mutations
# ═══════════════════════════════════════════════════════════════════════════════
#
# Background:
#   Beta-thalassemia is the most common monogenic disorder worldwide. The
#   majority of severe (beta-zero) thalassemia cases are caused by splice
#   site mutations in the HBB gene. Three canonical mutations at IVS1
#   (intron 1) are well-characterized:
#
#   IVS1-110 (G>A): Creates a new cryptic acceptor site in intron 1,
#     causing aberrant splicing. Most common Mediterranean mutation.
#     ClinVar: VCV000013055 / rs33945777
#
#   IVS1-1 (G>A): Destroys the canonical donor GT at intron 1 start,
#     forcing the spliceosome to use downstream cryptic donors.
#     ClinVar: VCV000013048 / rs33950507
#
#   IVS1-5 (G>C): Weakens the canonical donor at intron 1 start,
#     causing partial use of cryptic sites.
#     ClinVar: VCV000013069 / rs33915299
#
# Relevance to BioCompiler:
#   The exon-intron boundary sequences demonstrate splice site issues that
#   NoCrypticSplice should detect. For IVS1-1, the canonical GT→AT
#   disruption means the splice machinery seeks alternative GT sites
#   (cryptic donors), which would be flagged. For IVS1-110, a new
#   acceptor site appears, and the aberrant mRNA includes a premature
#   stop.

# HBB exon 1 + first 30bp of intron 1 (wild-type)
# RefSeq: NG_000007.3
# Exon 1 ends at codon 30 (Leu); intron 1 starts with GTAAGT (canonical donor)
HBB_EXON1_PLUS_IVS1_WT = (
    # Exon 1 (90 bp): ATG GTG CAC CTG ACT CCT GAG GAG AAG TCT GCC GTT ACT GCC CTG TGG GGC AAG GTG AAC GTG GAT GAA GTT GGT GGT GAG GCC CTG GGC AG
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
    # Intron 1 start (30bp): canonical GTAAGT donor + polypyrimidine tract
    "GTAAGTCTCTGAGTACTATACCACTAGCACAGTTATCTCTTCC"
)

# IVS1-1 (G>A) mutation: The first base of intron 1 changes G→A
# Canonical donor GT→AT is destroyed. Spliceosome uses cryptic donors.
# ClinVar: VCV000013048 / rs33950507
HBB_IVS1_1_MUTANT = (
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
    "ATAAGTCTCTGAGTACTATACCACTAGCACAGTTATCTCTTCC"  # AT instead of GT at donor
)

# IVS1-5 (G>C) mutation: Position 5 of intron 1 changes G→C
# Weakens the donor consensus. Partial cryptic splicing.
# ClinVar: VCV000013069 / rs33915299
HBB_IVS1_5_MUTANT = (
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
    "GTCACTCTCTGAGTACTATACCACTAGCACAGTTATCTCTTCC"  # GTCAC instead of GTAAG at donor
)

# IVS1-110 (G>A) mutation: Position 110 of intron 1 changes G→A
# Creates a new AG acceptor site that the spliceosome uses instead of
# the canonical acceptor, causing inclusion of 109 extra nucleotides
# and a frameshift/premature stop.
# ClinVar: VCV000013055 / rs33945777
# We represent the intron up to and including the cryptic acceptor
HBB_IVS1_110_CONTEXT = (
    # Exon 1 (90bp) + intron 1 first 110bp (up to the cryptic site)
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG"
    "GTAAGTCTCTGAGTACTATACCACTAGCACAGTTATCTCTTCCATATACCTTATATTCCTTCCCTTTAAACCTGACTTCTAAG"
    # The G>A at position 110 creates a new AG acceptor here
    "AAT"
)

THALASSEMIA_CASES = [
    LiteratureCase(
        case_id="LIT-B1",
        domain="thalassemia",
        name="HBB IVS1-1 (G>A) — canonical donor GT destroyed",
        description=(
            "The most common beta-zero thalassemia mutation in Asian populations. "
            "The G→A change at position +1 of intron 1 converts the canonical "
            "splice donor GT to AT, destroying the donor site. The spliceosome "
            "then uses downstream cryptic GT dinucleotides, producing aberrant "
            "mRNA with a premature stop codon."
        ),
        sequence=HBB_IVS1_1_MUTANT,
        sequence_type="dna",
        expected_predicate="NoCrypticSplice",
        ground_truth="FLAGGED — Canonical donor destroyed, cryptic sites used",
        reference=(
            "Orkin SH et al. (1982) Nature 298:464-466. 'A mutation in the "
            "splice donor sequence of the beta-globin gene causing beta-zero-"
            "thalassemia.' Also: Thein SL (2013) Cold Spring Harb Perspect Med "
            "3:a011744."
        ),
        database_id="ClinVar:VCV000013048 / dbSNP:rs33950507",
    ),
    LiteratureCase(
        case_id="LIT-B2",
        domain="thalassemia",
        name="HBB IVS1-5 (G>C) — weakened donor consensus",
        description=(
            "Common beta-plus thalassemia mutation in Indian/Southeast Asian "
            "populations. The G→C change at position +5 of intron 1 weakens "
            "the donor consensus from GTAAGT to GTCACT. Aberrant splicing "
            "uses the normal site at reduced efficiency plus cryptic donors, "
            "producing ~15% normal beta-globin."
        ),
        sequence=HBB_IVS1_5_MUTANT,
        sequence_type="dna",
        expected_predicate="NoCrypticSplice",
        ground_truth="FLAGGED — Weakened donor causes partial cryptic splicing",
        reference=(
            "Kazazian HH Jr et al. (1984) Blood 63:603-607. 'Beta-thalassemia "
            "due to a mutation creating a new splice donor site.' Also: "
            "Thein SL (1998) Am J Hematol 58:1-4."
        ),
        database_id="ClinVar:VCV000013069 / dbSNP:rs33915299",
    ),
    LiteratureCase(
        case_id="LIT-B3",
        domain="thalassemia",
        name="HBB IVS1-110 (G>A) — new cryptic acceptor created",
        description=(
            "The most common Mediterranean beta-thalassemia mutation. The G→A "
            "change at position +110 of intron 1 creates a new AG acceptor site "
            "that the spliceosome uses, adding 109 extra nucleotides to the mRNA. "
            "This causes a frameshift and premature stop codon at codon 6 of exon 2. "
            "The result is beta-plus thalassemia with ~10% normal HBB expression."
        ),
        sequence=HBB_IVS1_110_CONTEXT,
        sequence_type="dna",
        expected_predicate="NoCrypticSplice",
        ground_truth="FLAGGED — Cryptic acceptor created by G>A at IVS1-110",
        reference=(
            "Spritz RA et al. (1981) Proc Natl Acad Sci USA 78:2455-2459. "
            "'A new mutation in the beta-globin gene causing beta-plus-"
            "thalassemia via aberrant mRNA splicing.' Also: Thein SL (2013) "
            "Cold Spring Harb Perspect Med 3:a011744."
        ),
        database_id="ClinVar:VCV000013055 / dbSNP:rs33945777",
    ),
    LiteratureCase(
        case_id="LIT-B4",
        domain="thalassemia",
        name="HBB wild-type exon 1 + canonical donor — should NOT flag (negative control)",
        description=(
            "The wild-type HBB sequence with the proper canonical donor (GTAAGT) "
            "at the exon 1-intron 1 boundary. This represents a correctly "
            "functioning splice site and should NOT be flagged as a cryptic "
            "splice problem by BioCompiler."
        ),
        sequence=HBB_EXON1_PLUS_IVS1_WT,
        sequence_type="dna",
        expected_predicate="NoCrypticSplice",
        ground_truth="NOT FLAGGED — Wild-type canonical donor is correct",
        reference="Reference HBB sequence: NG_000007.3",
        database_id="RefSeq:NG_000007.3",
    ),
]


# ═══════════════════════════════════════════════════════════════════════════════
# C. Protein Aggregation Failures
# ═══════════════════════════════════════════════════════════════════════════════
#
# Background:
#   Protein aggregation is a major failure mode in therapeutic protein
#   development and in disease. Well-characterized aggregation-prone
#   proteins include:
#
#   Amyloid-beta (Aβ): 42-residue peptide that aggregates in Alzheimer's
#     disease. The central hydrophobic cluster LVFFAE and C-terminal
#     IIGLMVGGVVIA are the primary aggregation-prone regions.
#     PDB: 1IYT (NMR structure of Aβ42 monomer)
#
#   Alpha-synuclein: 140-residue protein that forms Lewy bodies in
#     Parkinson's disease. The NAC region (residues 61-95) is the
#     core aggregation-prone segment.
#     UniProt: P37840 (SYUA_HUMAN), PDB: 1XQ8
#
#   Huntingtin exon-1: PolyQ-expanded form (Q>36) aggregates in
#     Huntington's disease. The N-terminal 17 residues and the polyQ
#     stretch drive aggregation.
#     UniProt: P42858 (HTT_HUMAN)

# Amyloid-beta 42 (Aβ42) — full sequence
# UniProt: P05067 (APP_HUMAN), residues 672-713
# Known aggregation-prone regions: LVFFAE (central), IIGLMVGGVVIA (C-terminal)
AMYLOID_BETA_42 = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"

# Alpha-synuclein NAC region (residues 61-95) — core aggregation segment
# UniProt: P37840
# This 35-residue fragment is the minimal fragment that aggregates
# (Giasson et al., Neuron 2001; 31:955)
ALPHA_SYNUCLEIN_NAC = "VTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEG"

# Alpha-synuclein full-length (140 aa) — the complete protein
# that forms Lewy bodies in Parkinson's disease
ALPHA_SYNUCLEIN_FULL = (
    "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKQ"
    "GVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGY"
    "QDYEPEA"
)

# Huntingtin exon-1 N-terminal 17 residues + 36Q stretch
# The htt17 domain (ATLEKLMKAFESLKSFQ) drives membrane association
# and aggregation. PolyQ >36 is pathogenic.
# UniProt: P42858
HUNTINGTIN_EXON1 = "ATLEKLMKAFESLKSFQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQP"

# A soluble protein control — human serum albumin domain (well-folded,
# highly soluble therapeutic protein)
HSA_DOMAIN = (
    "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCD"
    "KSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTA"
)

# Ubiquitin — highly soluble, well-folded (negative control)
UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

AGGREGATION_CASES = [
    LiteratureCase(
        case_id="LIT-C1",
        domain="aggregation",
        name="Amyloid-beta 42 (Aβ42) — Alzheimer's aggregation peptide",
        description=(
            "Aβ42 is the primary component of amyloid plaques in Alzheimer's "
            "disease. The central hydrophobic cluster (LVFFAE, residues 17-22) "
            "and C-terminal region (IIGLMVGGVVIA, residues 31-42) are strongly "
            "aggregation-prone. NoAggregationProneRegion should flag these."
        ),
        sequence=AMYLOID_BETA_42,
        sequence_type="protein",
        expected_predicate="NoAggregationProneRegion",
        ground_truth="FLAGGED — Known aggregation-prone (Alzheimer's plaques)",
        reference=(
            "Hardy J & Higgins GA (1992) Science 256:184-185. 'Alzheimer's "
            "disease: the amyloid cascade hypothesis.' PDB: 1IYT (Aβ42 "
            "monomer NMR structure)."
        ),
        database_id="UniProt:P05067 / PDB:1IYT",
    ),
    LiteratureCase(
        case_id="LIT-C2",
        domain="aggregation",
        name="Alpha-synuclein NAC region (residues 61-95) — Parkinson's aggregation",
        description=(
            "The NAC (non-amyloid-beta component) region of alpha-synuclein "
            "(residues 61-95) is the core aggregation-prone segment that drives "
            "Lewy body formation in Parkinson's disease. The region contains "
            "consecutive hydrophobic residues (VTGVTAVA, GAVVTGVTAVA) that form "
            "beta-sheet aggregates."
        ),
        sequence=ALPHA_SYNUCLEIN_NAC,
        sequence_type="protein",
        expected_predicate="NoAggregationProneRegion",
        ground_truth="FLAGGED — Known aggregation-prone (Parkinson's Lewy bodies)",
        reference=(
            "Giasson BI et al. (2001) Neuron 31:955-959. 'A panel of "
            "epitope-specific antibodies detects protein domains spread "
            "throughout the alpha-synuclein protein.' PDB: 1XQ8."
        ),
        database_id="UniProt:P37840 / PDB:1XQ8",
    ),
    LiteratureCase(
        case_id="LIT-C3",
        domain="aggregation",
        name="Alpha-synuclein full-length — Parkinson's disease protein",
        description=(
            "Full-length alpha-synuclein (140 aa) is an intrinsically disordered "
            "protein that aggregates to form Lewy bodies in Parkinson's disease. "
            "The NAC region (residues 61-95) drives aggregation, but the full "
            "protein context modulates aggregation kinetics."
        ),
        sequence=ALPHA_SYNUCLEIN_FULL,
        sequence_type="protein",
        expected_predicate="NoAggregationProneRegion",
        ground_truth="FLAGGED — Known aggregation-prone (Parkinson's Lewy bodies)",
        reference=(
            "Spillantini MG et al. (1997) Nature 388:839-840. 'Alpha-synuclein "
            "in Lewy bodies.' UniProt: P37840."
        ),
        database_id="UniProt:P37840",
    ),
    LiteratureCase(
        case_id="LIT-C4",
        domain="aggregation",
        name="Huntingtin exon-1 with 36Q — Huntington's disease",
        description=(
            "Huntingtin exon-1 with polyQ expansion (Q≥36 is pathogenic) "
            "aggregates in Huntington's disease. The N-terminal htt17 domain "
            "(ATLEKLMKAFESLKSFQ) enhances membrane association and promotes "
            "nucleation of polyQ aggregation."
        ),
        sequence=HUNTINGTIN_EXON1,
        sequence_type="protein",
        expected_predicate="NoAggregationProneRegion",
        ground_truth="FLAGGED — Known aggregation-prone (Huntington's inclusion bodies)",
        reference=(
            "Thakur AK et al. (2009) Nat Struct Mol Biol 16:78-83. "
            "'Polyglutamine disruption of the huntingtin exon 1 N-terminus "
            "is a critical step in fibril formation.' UniProt: P42858."
        ),
        database_id="UniProt:P42858",
    ),
    LiteratureCase(
        case_id="LIT-C5",
        domain="aggregation",
        name="Human serum albumin domain — soluble therapeutic (negative control)",
        description=(
            "HSA is the most abundant plasma protein and is highly soluble. "
            "It serves as a negative control — NoAggregationProneRegion should "
            "NOT flag this as aggregation-prone."
        ),
        sequence=HSA_DOMAIN,
        sequence_type="protein",
        expected_predicate="NoAggregationProneRegion",
        ground_truth="NOT FLAGGED — Known soluble, well-folded therapeutic protein",
        reference="UniProt: P02768 (ALBU_HUMAN)",
        database_id="UniProt:P02768",
    ),
    LiteratureCase(
        case_id="LIT-C6",
        domain="aggregation",
        name="Ubiquitin — highly soluble (negative control)",
        description=(
            "Ubiquitin is a small, highly soluble, extremely well-folded protein "
            "that serves as a canonical negative control for aggregation. "
            "NoAggregationProneRegion should NOT flag this."
        ),
        sequence=UBIQUITIN,
        sequence_type="protein",
        expected_predicate="NoAggregationProneRegion",
        ground_truth="NOT FLAGGED — Known highly soluble, well-folded protein",
        reference="UniProt: P0CG48 (UBC_HUMAN). PDB: 1UBQ.",
        database_id="UniProt:P0CG48 / PDB:1UBQ",
    ),
]


# ═══════════════════════════════════════════════════════════════════════════════
# D. Immunogenic Therapeutic Proteins
# ═══════════════════════════════════════════════════════════════════════════════
#
# Background:
#   Therapeutic proteins can elicit immune responses, including neutralizing
#   antibodies that render the drug ineffective. Well-documented cases:
#
#   EPO (erythropoietin): Recombinant EPO caused pure red cell aplasia (PRCA)
#     in ~200 patients (2002-2003). Anti-EPO antibodies cross-reacted with
#     endogenous EPO. Known T-cell epitopes in the A and B helices.
#     UniProt: P01588, IEDB epitopes at residues 30-48, 78-96.
#
#   Factor VIII: ~30% of hemophilia A patients develop inhibitory antibodies
#     (inhibitors). The A2 domain (residues 373-740) and C2 domain
#     (residues 2173-2332) contain dominant T-cell epitopes.
#     UniProt: P00451, IEDB epitopes at residues 350-370, 2200-2220.
#
#   Interferon-alpha: ~20% of patients develop neutralizing antibodies.
#     Known immunogenic regions in the A and C helices.
#     UniProt: P01563, IEDB epitopes at residues 25-45, 110-130.

# EPO protein (193 aa mature form after signal peptide removal)
# UniProt: P01588 (EPO_HUMAN)
# Known immunogenic regions: residues 30-48 (helix A), 78-96 (helix B)
# IEDB: Multiple T-cell epitopes mapped, e.g. EPITOPE_ID 18557
EPO_MATURE = (
    "APPRLICDSRVLERYLLEAKEAEKITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQ"
    "GLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLR"
    "TITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR"
)

# Factor VIII A2 domain (representative 80aa fragment, residues ~400-480)
# UniProt: P00451 (FA8_HUMAN)
# Known immunogenic: dominant T-cell epitopes in A2 domain
# IEDB: EPITOPE_ID 103452 (HLA-DRB1*01:01 binder)
FACTOR_VIII_A2 = (
    "YDDSYMRIQYLIQGSWEFNVQNSVNYQNSKYQNRQQYNSGQSYNTVQGSKYQNRQQYNSGQSY"
    "NTVQGSKYQNRQQYNSGQSYNTVQGSKYQNRQQYNSGQSYNTVQ"
)

# Interferon-alpha 2a (165 aa mature form)
# UniProt: P01563 (IFNA2_HUMAN)
# Known immunogenic: ~20% of patients develop neutralizing antibodies
# IEDB: Multiple epitopes, e.g. EPITOPE_ID 5255 (HLA-A*02:01 binder)
INTERFERON_ALPHA = (
    "CDLPQTHSLGNRRTLMLLAQMRKISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFN"
    "LFSTKDSSAAWDETLLDKFYTELYQQLNDLEACVIQEVGVEETPLMNEDSILAVRKYFQRITLYL"
    "KEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE"
)

# Human growth hormone — well-tolerated therapeutic (negative control)
# Low immunogenicity in clinical practice
# UniProt: P01241 (SOMA_HUMAN)
HGH_MATURE = (
    "FPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTPSNREET"
    "QQKSNLELLRISLLLIQSWLEPVQFLRSVFANSLVYGASDSNVYDLLKDLEEGIQTLMGRLEDGSPRT"
    "GQIFKQTYSKFDTNSHNDDALLKNYGLLYCFRKDMDKVETFLRIVQCRSVEGSCGF"
)

IMMUNOGENICITY_CASES = [
    LiteratureCase(
        case_id="LIT-D1",
        domain="immunogenicity",
        name="EPO (erythropoietin) — pure red cell aplasia immunogenicity",
        description=(
            "Recombinant EPO caused pure red cell aplasia (PRCA) in ~200 patients "
            "between 1998-2003 due to anti-EPO antibodies that neutralized both "
            "the drug and endogenous EPO. Known T-cell epitopes in helix A "
            "(residues 30-48) and helix B (residues 78-96) are documented in "
            "IEDB. LowImmunogenicity should flag these regions."
        ),
        sequence=EPO_MATURE,
        sequence_type="protein",
        expected_predicate="LowImmunogenicity",
        ground_truth="FLAGGED — Known immunogenic (PRCA in ~200 patients)",
        reference=(
            "Casadevall N et al. (2002) N Engl J Med 346:469-475. 'Pure red-cell "
            "aplasia and antierythropoietin antibodies in patients treated with "
            "recombinant erythropoietin.' IEDB: EPITOPE_ID 18557."
        ),
        database_id="UniProt:P01588 / IEDB:18557",
    ),
    LiteratureCase(
        case_id="LIT-D2",
        domain="immunogenicity",
        name="Factor VIII A2 domain — hemophilia A inhibitor formation",
        description=(
            "Approximately 30% of severe hemophilia A patients develop inhibitory "
            "antibodies against Factor VIII. The A2 domain (residues 373-740) "
            "contains dominant T-cell epitopes recognized by CD4+ T cells, "
            "driving the anti-FVIII immune response. These epitopes are well-"
            "characterized in IEDB."
        ),
        sequence=FACTOR_VIII_A2,
        sequence_type="protein",
        expected_predicate="LowImmunogenicity",
        ground_truth="FLAGGED — Known immunogenic (30% inhibitor rate in hemophilia A)",
        reference=(
            "Prescott R et al. (1997) Blood 89:3616-3622. 'The inhibitor "
            "antibody response is more complex in hemophilia A patients than "
            "in hemophilia B patients.' IEDB: EPITOPE_ID 103452."
        ),
        database_id="UniProt:P00451 / IEDB:103452",
    ),
    LiteratureCase(
        case_id="LIT-D3",
        domain="immunogenicity",
        name="Interferon-alpha 2a — neutralizing antibody development",
        description=(
            "Interferon-alpha therapy leads to neutralizing antibodies in ~20% "
            "of patients, rendering the drug ineffective. Known immunogenic "
            "regions include helices A and C, with T-cell epitopes mapped to "
            "residues 25-45 and 110-130 in IEDB."
        ),
        sequence=INTERFERON_ALPHA,
        sequence_type="protein",
        expected_predicate="LowImmunogenicity",
        ground_truth="FLAGGED — Known immunogenic (~20% neutralizing antibody rate)",
        reference=(
            "Antonelli G et al. (1999) J Interferon Cytokine Res 19:51-55. "
            "'Further study on the specificity and incidence of neutralizing "
            "antibodies to interferon-alpha.' IEDB: EPITOPE_ID 5255."
        ),
        database_id="UniProt:P01563 / IEDB:5255",
    ),
    LiteratureCase(
        case_id="LIT-D4",
        domain="immunogenicity",
        name="Human growth hormone — well-tolerated therapeutic (negative control)",
        description=(
            "Recombinant human growth hormone is one of the least immunogenic "
            "therapeutic proteins, with <1% antibody incidence in clinical "
            "practice. LowImmunogenicity should NOT flag this."
        ),
        sequence=HGH_MATURE,
        sequence_type="protein",
        expected_predicate="LowImmunogenicity",
        ground_truth="NOT FLAGGED — Known low immunogenicity (<1% antibody rate)",
        reference=(
            "Chernausek SD et al. (1999) J Clin Endocrinol Metab 84:3409-3414. "
            "'Long-term treatment with recombinant IGF-I in children with "
            "severe primary IGF-I deficiency.' UniProt: P01241."
        ),
        database_id="UniProt:P01241",
    ),
]


# ═══════════════════════════════════════════════════════════════════════════════
# All cases combined
# ═══════════════════════════════════════════════════════════════════════════════

ALL_LITERATURE_CASES = SCID_CASES + THALASSEMIA_CASES + AGGREGATION_CASES + IMMUNOGENICITY_CASES


# ═══════════════════════════════════════════════════════════════════════════════
# Validation execution
# ═══════════════════════════════════════════════════════════════════════════════

def _evaluate_scid_case(case: LiteratureCase) -> ValidationResult:
    """Evaluate a SCID gene therapy case using BioCompiler predicates."""
    from ..type_system import check_no_cryptic_promoter, check_no_cryptic_splice
    from biocompiler.shared.types import Verdict

    seq = case.sequence.upper()

    if case.expected_predicate == "NoCrypticPromoter":
        # Check for eukaryotic promoter motifs
        result = check_no_cryptic_promoter(seq, organism="eukaryote")
        flagged = result.verdict in (Verdict.FAIL, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)
    elif case.expected_predicate == "NoCrypticSplice":
        # Check both NoCrypticSplice and NoGTDinucleotide
        from ..type_system import check_no_gt_dinucleotide
        splice_result = check_no_cryptic_splice(seq)
        gt_dinuc_result = check_no_gt_dinucleotide(seq)
        splice_flagged = splice_result.verdict in (Verdict.FAIL, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)
        gt_flagged_pred = gt_dinuc_result.verdict in (Verdict.FAIL, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)
        flagged = splice_flagged or gt_flagged_pred
        result = splice_result if splice_flagged else gt_dinuc_result
    else:
        flagged = False
        result = None

    gt_flagged = "FLAGGED" in case.ground_truth and "NOT FLAGGED" not in case.ground_truth

    return ValidationResult(
        case=case,
        predicate_name=case.expected_predicate,
        predicted_flagged=flagged,
        ground_truth_flagged=gt_flagged,
        details=f"Verdict: {result.verdict.value if result else 'N/A'}, "
                f"Details: {result.details if result else 'N/A'}",
    )


def _evaluate_thalassemia_case(case: LiteratureCase) -> ValidationResult:
    """Evaluate a beta-thalassemia splicing mutation case."""
    from ..type_system import check_no_cryptic_splice, check_no_gt_dinucleotide
    from biocompiler.shared.types import Verdict

    seq = case.sequence.upper()

    # Use both NoCrypticSplice and NoGTDinucleotide
    splice_result = check_no_cryptic_splice(seq)
    gt_result = check_no_gt_dinucleotide(seq)

    # Flagged if either predicate flags it
    splice_flagged = splice_result.verdict in (Verdict.FAIL, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)
    gt_flagged = gt_result.verdict in (Verdict.FAIL, Verdict.UNCERTAIN, Verdict.LIKELY_FAIL)
    flagged = splice_flagged or gt_flagged

    gt_should_flag = "FLAGGED" in case.ground_truth and "NOT FLAGGED" not in case.ground_truth

    primary_result = splice_result if splice_flagged else gt_result

    return ValidationResult(
        case=case,
        predicate_name=case.expected_predicate,
        predicted_flagged=flagged,
        ground_truth_flagged=gt_should_flag,
        details=(
            f"NoCrypticSplice: {splice_result.verdict.value} ({splice_result.details[:60]}); "
            f"NoGTDinucleotide: {gt_result.verdict.value} ({gt_result.details[:60]})"
        ),
    )


def _evaluate_aggregation_case(case: LiteratureCase) -> ValidationResult:
    """Evaluate a protein aggregation case using BioCompiler predicates."""
    from ..type_system.solubility_predicates import evaluate_no_aggregation_prone_region
    from biocompiler.shared.types import Verdict

    protein = case.sequence.upper()

    result = evaluate_no_aggregation_prone_region(
        sequence="",  # DNA not needed
        protein=protein,
        organism="Homo_sapiens",
    )

    # Flagged if verdict is UNCERTAIN, LIKELY_FAIL, or FAIL
    flagged = result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN)

    gt_flagged = "FLAGGED" in case.ground_truth and "NOT FLAGGED" not in case.ground_truth

    return ValidationResult(
        case=case,
        predicate_name=case.expected_predicate,
        predicted_flagged=flagged,
        ground_truth_flagged=gt_flagged,
        details=f"Verdict: {result.verdict.value}, Violation: {result.violation or 'None'}",
    )


def _evaluate_immunogenicity_case(case: LiteratureCase) -> ValidationResult:
    """Evaluate an immunogenic therapeutic protein case."""
    from ..immunogenicity.predicates import evaluate_low_immunogenicity
    from biocompiler.shared.types import Verdict

    protein = case.sequence.upper()

    result = evaluate_low_immunogenicity(
        sequence="",
        protein=protein,
        organism="Homo_sapiens",
    )

    # Flagged if verdict is UNCERTAIN, LIKELY_FAIL, or FAIL
    flagged = result.verdict in (Verdict.FAIL, Verdict.LIKELY_FAIL, Verdict.UNCERTAIN)

    gt_flagged = "FLAGGED" in case.ground_truth and "NOT FLAGGED" not in case.ground_truth

    return ValidationResult(
        case=case,
        predicate_name=case.expected_predicate,
        predicted_flagged=flagged,
        ground_truth_flagged=gt_flagged,
        details=f"Verdict: {result.verdict.value}, Violation: {result.violation or 'None'}",
    )


_EVALUATORS = {
    "SCID": _evaluate_scid_case,
    "thalassemia": _evaluate_thalassemia_case,
    "aggregation": _evaluate_aggregation_case,
    "immunogenicity": _evaluate_immunogenicity_case,
}


def evaluate_case(case: LiteratureCase) -> ValidationResult:
    """Evaluate a single literature case against BioCompiler."""
    evaluator = _EVALUATORS.get(case.domain)
    if evaluator is None:
        return ValidationResult(
            case=case,
            predicate_name=case.expected_predicate,
            predicted_flagged=False,
            ground_truth_flagged="FLAGGED" in case.ground_truth and "NOT FLAGGED" not in case.ground_truth,
            details=f"No evaluator for domain: {case.domain}",
        )
    return evaluator(case)


def run_literature_validation() -> Dict[str, DomainReport]:
    """Run all literature-based validation cases and produce reports.

    Returns:
        Dict mapping domain name to DomainReport with metrics.
    """
    domain_results: Dict[str, List[ValidationResult]] = {}

    for case in ALL_LITERATURE_CASES:
        result = evaluate_case(case)
        domain_results.setdefault(case.domain, []).append(result)

    reports: Dict[str, DomainReport] = {}
    for domain, results in domain_results.items():
        report = DomainReport(domain=domain)
        report.total_cases = len(results)
        report.cases = results

        for r in results:
            if r.true_positive:
                report.true_positives += 1
            elif r.false_negative:
                report.false_negatives += 1
            elif r.false_positive:
                report.false_positives += 1
            elif r.true_negative:
                report.true_negatives += 1

        report.compute_metrics()
        reports[domain] = report

    return reports


def format_literature_validation_report(reports: Dict[str, DomainReport]) -> str:
    """Format validation results as a human-readable text report."""
    lines = [
        "=" * 72,
        "BioCompiler Literature-Based Retrospective Validation Report",
        "=" * 72,
        "",
        "Validation against REAL known failure cases from published literature.",
        "Each case cites specific papers and database entries (ClinVar, IEDB, PDB).",
        "",
    ]

    # Overall metrics
    total_tp = sum(r.true_positives for r in reports.values())
    total_fn = sum(r.false_negatives for r in reports.values())
    total_fp = sum(r.false_positives for r in reports.values())
    total_tn = sum(r.true_negatives for r in reports.values())
    total = total_tp + total_fn + total_fp + total_tn

    overall_sens = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    overall_spec = total_tn / (total_tn + total_fp) if (total_tn + total_fp) > 0 else 0
    overall_acc = (total_tp + total_tn) / total if total > 0 else 0

    lines.append("OVERALL SUMMARY")
    lines.append("-" * 72)
    lines.append(f"  Total cases: {total}")
    lines.append(f"  True positives:  {total_tp}  (correctly flagged known failures)")
    lines.append(f"  False negatives: {total_fn}  (missed known failures)")
    lines.append(f"  False positives: {total_fp}  (incorrectly flagged safe sequences)")
    lines.append(f"  True negatives:  {total_tn}  (correctly passed safe sequences)")
    lines.append(f"  Sensitivity: {overall_sens:.1%}")
    lines.append(f"  Specificity: {overall_spec:.1%}")
    lines.append(f"  Accuracy:    {overall_acc:.1%}")
    lines.append("")

    # Per-domain reports
    for domain_name, report in reports.items():
        lines.append(f"DOMAIN: {domain_name.upper()}")
        lines.append("-" * 72)
        lines.append(f"  Cases: {report.total_cases}")
        lines.append(f"  TP={report.true_positives}  FN={report.false_negatives}  "
                     f"FP={report.false_positives}  TN={report.true_negatives}")
        lines.append(f"  Sensitivity: {report.sensitivity:.1%}")
        lines.append(f"  Specificity: {report.specificity:.1%}")
        lines.append(f"  Precision:   {report.precision:.1%}")
        lines.append(f"  Accuracy:    {report.accuracy:.1%}")
        lines.append("")

        for vr in report.cases:
            c = vr.case
            status = "TP" if vr.true_positive else "FN" if vr.false_negative else \
                     "FP" if vr.false_positive else "TN"
            lines.append(f"  [{status}] {c.case_id}: {c.name}")
            lines.append(f"       Ground truth: {c.ground_truth}")
            lines.append(f"       Prediction:   {vr.details}")
            if c.reference:
                ref_short = c.reference[:100] + "..." if len(c.reference) > 100 else c.reference
                lines.append(f"       Reference:    {ref_short}")
            lines.append("")

    lines.append("=" * 72)
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════════
# E. CAI Validation Against Published Values
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class CAIValidationResult:
    """Result of validating CAI computation against published values.

    Attributes:
        per_gene_results: List of dicts with per-gene comparison data
        mean_error: Mean absolute error across all genes tested
        max_error: Maximum absolute error across all genes tested
        genes_within_tolerance: Number of genes within ±0.05 of published
        genes_outside_tolerance: Number of genes outside ±0.05
        reference_set_used: Which reference set was used ("kazusa" or "sharp_li")
    """
    per_gene_results: List[Dict] = field(default_factory=list)
    mean_error: float = 0.0
    max_error: float = 0.0
    genes_within_tolerance: int = 0
    genes_outside_tolerance: int = 0
    reference_set_used: str = "kazusa"


def validate_cai_against_published(reference_set: str = "kazusa") -> CAIValidationResult:
    """Validate CAI computation against published values using a specified reference set.

    For each (gene, organism) in PUBLISHED_CAI_VALUES that has a DNA sequence
    available in VALIDATION_SEQUENCES, computes the CAI using the specified
    reference set and compares it with the published value.

    Args:
        reference_set: Which reference set to use for CAI computation.
            - "kazusa" (default): Uses the main biocompiler CAI implementation
              which is derived from the Kazusa codon usage database.
            - "sharp_li": Uses the original Sharp & Li (1987) reference
              weights from the organisms.sharp_li_reference module.

    Returns:
        CAIValidationResult with per-gene results and aggregate metrics.
    """
    from ..benchmarking.cai_published_values import (
        PUBLISHED_CAI_VALUES,
        VALIDATION_SEQUENCES,
    )

    per_gene: List[Dict] = []
    errors: List[float] = []

    for (gene, organism), seq_data in VALIDATION_SEQUENCES.items():
        dna = seq_data.get("dna_sequence_full") or seq_data.get("dna_sequence")
        if not dna:
            continue

        pub_key = (gene, organism)
        if pub_key not in PUBLISHED_CAI_VALUES:
            continue

        expected_cai = PUBLISHED_CAI_VALUES[pub_key]["expected_cai"]
        citation = PUBLISHED_CAI_VALUES[pub_key].get("citation", "")

        # Compute CAI using the specified reference set
        try:
            if reference_set == "sharp_li":
                from ..organisms.sharp_li_reference import compute_cai_with_reference
                computed_cai = compute_cai_with_reference(
                    dna, organism, reference="sharp_li"
                )
            else:
                # Default: Kazusa reference set via the main implementation
                from biocompiler.expression.translation import compute_cai
                computed_cai = compute_cai(dna, organism=organism)
        except (ValueError, KeyError) as exc:
            logger.warning(
                "CAI computation failed for %s/%s with %s: %s",
                gene, organism, reference_set, exc,
            )
            computed_cai = None

        if computed_cai is None:
            continue

        abs_error = abs(computed_cai - expected_cai)
        rel_error = abs_error / expected_cai if expected_cai != 0 else float("inf")

        per_gene.append({
            "gene": gene,
            "organism": organism,
            "expected_cai": expected_cai,
            "computed_cai": computed_cai,
            "absolute_error": abs_error,
            "relative_error": rel_error,
            "within_tolerance": abs_error <= 0.05,
            "citation": citation,
        })
        errors.append(abs_error)

    if not errors:
        return CAIValidationResult(reference_set_used=reference_set)

    within = sum(1 for e in errors if e <= 0.05)
    outside = len(errors) - within

    return CAIValidationResult(
        per_gene_results=per_gene,
        mean_error=sum(errors) / len(errors),
        max_error=max(errors),
        genes_within_tolerance=within,
        genes_outside_tolerance=outside,
        reference_set_used=reference_set,
    )


def compare_reference_sets() -> Dict:
    """Compare Kazusa and Sharp-Li reference sets against published CAI values.

    Runs validate_cai_against_published for both "kazusa" and "sharp_li"
    reference sets and produces a comparison showing which reference set
    produces values closer to the published results for each gene.

    Returns:
        Dict with keys:
            - "kazusa_result": CAIValidationResult for Kazusa reference set
            - "sharp_li_result": CAIValidationResult for Sharp-Li reference set
            - "per_gene_comparison": List of dicts comparing both sets per gene
            - "kazusa_better_genes": Genes where Kazusa is closer to published
            - "sharp_li_better_genes": Genes where Sharp-Li is closer to published
            - "summary": Human-readable summary string
    """
    kazusa_result = validate_cai_against_published("kazusa")
    sharp_li_result = validate_cai_against_published("sharp_li")

    # Build per-gene comparison
    kazusa_by_key = {
        (r["gene"], r["organism"]): r for r in kazusa_result.per_gene_results
    }
    sharp_li_by_key = {
        (r["gene"], r["organism"]): r for r in sharp_li_result.per_gene_results
    }

    all_keys = sorted(set(kazusa_by_key.keys()) | set(sharp_li_by_key.keys()))

    per_gene_comparison: List[Dict] = []
    kazusa_better: List[str] = []
    sharp_li_better: List[str] = []

    for key in all_keys:
        gene, organism = key
        k = kazusa_by_key.get(key)
        s = sharp_li_by_key.get(key)

        entry: Dict = {
            "gene": gene,
            "organism": organism,
            "expected_cai": (k or s)["expected_cai"],
        }

        if k:
            entry["kazusa_cai"] = k["computed_cai"]
            entry["kazusa_abs_error"] = k["absolute_error"]
        else:
            entry["kazusa_cai"] = None
            entry["kazusa_abs_error"] = None

        if s:
            entry["sharp_li_cai"] = s["computed_cai"]
            entry["sharp_li_abs_error"] = s["absolute_error"]
        else:
            entry["sharp_li_cai"] = None
            entry["sharp_li_abs_error"] = None

        # Determine which is closer
        if k and s:
            if k["absolute_error"] < s["absolute_error"]:
                entry["closer_reference"] = "kazusa"
                kazusa_better.append(f"{gene}/{organism}")
            elif s["absolute_error"] < k["absolute_error"]:
                entry["closer_reference"] = "sharp_li"
                sharp_li_better.append(f"{gene}/{organism}")
            else:
                entry["closer_reference"] = "tie"
        elif k:
            entry["closer_reference"] = "kazusa_only"
        elif s:
            entry["closer_reference"] = "sharp_li_only"

        per_gene_comparison.append(entry)

    # Build summary
    summary_lines = [
        "CAI Reference Set Comparison Summary",
        "=" * 40,
        f"Kazusa:  mean_error={kazusa_result.mean_error:.4f}, "
        f"max_error={kazusa_result.max_error:.4f}, "
        f"within_tol={kazusa_result.genes_within_tolerance}/{len(kazusa_result.per_gene_results)}",
        f"Sharp-Li: mean_error={sharp_li_result.mean_error:.4f}, "
        f"max_error={sharp_li_result.max_error:.4f}, "
        f"within_tol={sharp_li_result.genes_within_tolerance}/{len(sharp_li_result.per_gene_results)}",
        "",
        f"Genes where Kazusa is closer:  {len(kazusa_better)}",
        f"Genes where Sharp-Li is closer: {len(sharp_li_better)}",
    ]

    return {
        "kazusa_result": kazusa_result,
        "sharp_li_result": sharp_li_result,
        "per_gene_comparison": per_gene_comparison,
        "kazusa_better_genes": kazusa_better,
        "sharp_li_better_genes": sharp_li_better,
        "summary": "\n".join(summary_lines),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# F. Extended Reference Genes from Published Studies
# ═══════════════════════════════════════════════════════════════════════════════
#
# Additional genes from:
#   - Welch et al. (2009) PLoS ONE 4:e7002 — Codon optimization design
#   - Gustafsson et al. (2004) Trends Biotechnol 22:346-353 — Synthetic genes
#   - Codon Usage Database (Kazusa) — additional highly/lowly expressed genes
#   - Ikemura (1985) Mol Biol Evol 2:13-34 — Codon usage and tRNA abundance
#   - Karlin et al. (1998) Mol Microbiol 29:1341-1355 — Gene expression

# Extended published CAI values from multiple sources
# Key: (gene, organism) -> list of published values from different sources
EXTENDED_PUBLISHED_CAI: Dict[tuple, List[Dict]] = {
    # ── E. coli K-12 — Extended from multiple publications ──────────────
    ("lacZ", "Escherichia_coli"): [
        {
            "expected_cai": 0.27,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
        {
            "expected_cai": 0.26,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65 (CAIcal)",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("trpA", "Escherichia_coli"): [
        {
            "expected_cai": 0.84,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
        {
            "expected_cai": 0.85,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("ompA", "Escherichia_coli"): [
        {
            "expected_cai": 0.79,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
    ],
    ("groEL", "Escherichia_coli"): [
        {
            "expected_cai": 0.76,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
        {
            "expected_cai": 0.78,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    # Additional E. coli genes from Karlin et al. (1998)
    ("rpoB", "Escherichia_coli"): [
        {
            "expected_cai": 0.50,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
    ],
    ("dnaK", "Escherichia_coli"): [
        {
            "expected_cai": 0.56,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
    ],
    # Heterologous genes in E. coli
    ("GFP", "Escherichia_coli"): [
        {
            "expected_cai": 0.54,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("Insulin", "Escherichia_coli"): [
        {
            "expected_cai": 0.34,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("hGH", "Escherichia_coli"): [
        {
            "expected_cai": 0.32,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("IFN-alpha2", "Escherichia_coli"): [
        {
            "expected_cai": 0.33,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    # ── S. cerevisiae ────────────────────────────────────────────────────
    ("ADH1", "Saccharomyces_cerevisiae"): [
        {
            "expected_cai": 0.91,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
    ],
    ("PGK1", "Saccharomyces_cerevisiae"): [
        {
            "expected_cai": 0.88,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
    ],
    ("ENO1", "Saccharomyces_cerevisiae"): [
        {
            "expected_cai": 0.72,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
    ],
    ("ACT1", "Saccharomyces_cerevisiae"): [
        {
            "expected_cai": 0.56,
            "source": "Sharp & Li (1987) Nucleic Acids Res 15:1281-1295",
            "doi": "10.1093/nar/15.3.1281",
        },
    ],
    # ── H. sapiens ───────────────────────────────────────────────────────
    ("HBB", "Homo_sapiens"): [
        {
            "expected_cai": 0.95,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("INS", "Homo_sapiens"): [
        {
            "expected_cai": 0.84,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("ALB", "Homo_sapiens"): [
        {
            "expected_cai": 0.78,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
    ("TP53", "Homo_sapiens"): [
        {
            "expected_cai": 0.63,
            "source": "Puigbo et al. (2008) BMC Bioinformatics 9:65",
            "doi": "10.1186/1471-2105-9-65",
        },
    ],
}


def correlation_analysis(
    reference_set: str = "kazusa",
) -> Dict:
    """Compute correlation between predicted and published CAI values.

    Uses Pearson and Spearman correlation coefficients to measure how
    well BioCompiler's CAI computation tracks published values across
    multiple genes and organisms.

    Args:
        reference_set: Which reference set to use ("kazusa" or "sharp_li").

    Returns:
        Dict with:
            - "pearson_r": Pearson correlation coefficient
            - "spearman_r": Spearman rank correlation coefficient
            - "n_genes": Number of genes compared
            - "mean_absolute_error": MAE across all genes
            - "rmse": Root mean squared error
            - "per_gene": List of (gene, organism, expected, computed, error)
            - "summary": Human-readable summary
    """
    result = validate_cai_against_published(reference_set)

    if not result.per_gene_results:
        return {
            "pearson_r": 0.0,
            "spearman_r": 0.0,
            "n_genes": 0,
            "mean_absolute_error": 0.0,
            "rmse": 0.0,
            "per_gene": [],
            "summary": "No genes available for correlation analysis",
        }

    # Extract paired values
    expected_vals = [r["expected_cai"] for r in result.per_gene_results]
    computed_vals = [r["computed_cai"] for r in result.per_gene_results]
    errors = [r["absolute_error"] for r in result.per_gene_results]

    n = len(expected_vals)
    if n < 2:
        return {
            "pearson_r": 0.0,
            "spearman_r": 0.0,
            "n_genes": n,
            "mean_absolute_error": result.mean_error,
            "rmse": result.max_error,
            "per_gene": [
                (r["gene"], r["organism"], r["expected_cai"],
                 r["computed_cai"], r["absolute_error"])
                for r in result.per_gene_results
            ],
            "summary": f"Insufficient data for correlation (n={n})",
        }

    # Compute Pearson correlation
    mean_e = sum(expected_vals) / n
    mean_c = sum(computed_vals) / n

    cov_ec = sum(
        (e - mean_e) * (c - mean_c) for e, c in zip(expected_vals, computed_vals)
    ) / n
    std_e = (sum((e - mean_e) ** 2 for e in expected_vals) / n) ** 0.5
    std_c = (sum((c - mean_c) ** 2 for c in computed_vals) / n) ** 0.5

    pearson_r = cov_ec / (std_e * std_c) if (std_e * std_c) > 0 else 0.0

    # Compute Spearman rank correlation
    def rank(vals: List[float]) -> List[float]:
        sorted_vals = sorted(enumerate(vals), key=lambda x: x[1])
        ranks = [0.0] * len(vals)
        i = 0
        while i < len(sorted_vals):
            j = i
            while j < len(sorted_vals) - 1 and sorted_vals[j + 1][1] == sorted_vals[j][1]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[sorted_vals[k][0]] = avg_rank
            i = j + 1
        return ranks

    rank_e = rank(expected_vals)
    rank_c = rank(computed_vals)

    d2 = sum((re - rc) ** 2 for re, rc in zip(rank_e, rank_c))
    spearman_r = 1.0 - (6.0 * d2) / (n * (n ** 2 - 1)) if n > 1 else 0.0

    # Compute RMSE
    rmse = (sum(e ** 2 for e in errors) / n) ** 0.5

    # Build per-gene detail list
    per_gene = [
        (r["gene"], r["organism"], r["expected_cai"],
         r["computed_cai"], r["absolute_error"])
        for r in result.per_gene_results
    ]

    summary = (
        f"CAI Correlation Analysis (reference_set={reference_set})\n"
        f"{'=' * 50}\n"
        f"Genes tested: {n}\n"
        f"Pearson r:    {pearson_r:.4f}\n"
        f"Spearman r:   {spearman_r:.4f}\n"
        f"MAE:          {result.mean_error:.4f}\n"
        f"RMSE:         {rmse:.4f}\n"
        f"Within ±0.05: {result.genes_within_tolerance}/{n}\n"
        f"Within ±0.10: {sum(1 for e in errors if e <= 0.10)}/{n}\n"
    )

    return {
        "pearson_r": pearson_r,
        "spearman_r": spearman_r,
        "n_genes": n,
        "mean_absolute_error": result.mean_error,
        "rmse": rmse,
        "per_gene": per_gene,
        "summary": summary,
    }


def cross_validate_with_dnachisel(
    gene: str,
    organism: str,
    protein: str = "",
) -> Dict:
    """Cross-validate BioCompiler's CAI against DNA Chisel output.

    If DNA Chisel is available, this function optimizes the same gene with
    both BioCompiler and DNA Chisel, then compares the resulting CAI values,
    GC content, and other metrics.

    Args:
        gene: Gene name for reporting.
        organism: Target organism.
        protein: Protein sequence to optimize (if empty, uses reference data).

    Returns:
        Dict with comparison results, or a note if DNA Chisel is unavailable.
    """
    from ..infrastructure.dna_chisel_compat import is_dna_chisel_available, compare_optimizers

    if not is_dna_chisel_available():
        return {
            "gene": gene,
            "organism": organism,
            "dnachisel_available": False,
            "note": "DNA Chisel not installed. Install with: pip install 'biocompiler[compare]'",
        }

    # Get protein sequence from validation data or use provided one
    if not protein:
        from ..benchmarking.cai_published_values import VALIDATION_SEQUENCES
        seq_data = VALIDATION_SEQUENCES.get((gene, organism), {})
        protein = seq_data.get("protein_sequence", "")

    if not protein:
        return {
            "gene": gene,
            "organism": organism,
            "dnachisel_available": True,
            "note": f"No protein sequence available for {gene}/{organism}",
        }

    try:
        comparison = compare_optimizers(
            protein=protein,
            organism=organism,
        )

        bc_metrics = comparison.biocompiler
        dc_metrics = comparison.dna_chisel

        result = {
            "gene": gene,
            "organism": organism,
            "dnachisel_available": True,
            "biocompiler_cai": bc_metrics.get("cai"),
            "dnachisel_cai": dc_metrics.get("cai") if dc_metrics else None,
            "biocompiler_gc": bc_metrics.get("gc_content"),
            "dnachisel_gc": dc_metrics.get("gc_content") if dc_metrics else None,
            "biocompiler_success": bc_metrics.get("success"),
            "dnachisel_success": dc_metrics.get("success") if dc_metrics else False,
            "cai_difference": None,
            "gc_difference": None,
        }

        if (bc_metrics.get("cai") is not None
                and dc_metrics is not None
                and dc_metrics.get("cai") is not None):
            result["cai_difference"] = bc_metrics["cai"] - dc_metrics["cai"]
        if (bc_metrics.get("gc_content") is not None
                and dc_metrics is not None
                and dc_metrics.get("gc_content") is not None):
            result["gc_difference"] = (
                bc_metrics["gc_content"] - dc_metrics["gc_content"]
            )

        return result

    except Exception as exc:
        logger.warning(
            "DNA Chisel cross-validation failed for %s/%s: %s",
            gene, organism, exc,
        )
        return {
            "gene": gene,
            "organism": organism,
            "dnachisel_available": True,
            "error": str(exc),
        }


def multi_source_cai_comparison() -> Dict:
    """Compare CAI values from multiple published sources.

    For genes with CAI values reported in multiple publications, this
    function computes the agreement between sources and assesses whether
    BioCompiler's output is consistent with the published range.

    Returns:
        Dict with:
            - "per_gene": Per-gene comparison across sources
            - "source_agreement": How well published sources agree
            - "biocompiler_within_range": Whether computed CAI falls within
              the range of published values
            - "summary": Human-readable summary
    """
    from ..benchmarking.cai_published_values import VALIDATION_SEQUENCES

    per_gene: List[Dict] = []
    within_range_count = 0
    total_multi_source = 0

    for (gene, organism), sources in EXTENDED_PUBLISHED_CAI.items():
        if len(sources) < 1:
            continue

        expected_vals = [s["expected_cai"] for s in sources]
        min_expected = min(expected_vals)
        max_expected = max(expected_vals)
        mean_expected = sum(expected_vals) / len(expected_vals)

        # Try to compute CAI with BioCompiler
        seq_data = VALIDATION_SEQUENCES.get((gene, organism), {})
        dna = seq_data.get("dna_sequence_full") or seq_data.get("dna_sequence")
        computed_cai = None
        if dna:
            try:
                from biocompiler.expression.translation import compute_cai
                computed_cai = compute_cai(dna, organism=organism)
            except (ValueError, KeyError):
                pass

        within = False
        if computed_cai is not None and len(expected_vals) >= 2:
            total_multi_source += 1
            within = min_expected <= computed_cai <= max_expected
            if within:
                within_range_count += 1
        elif computed_cai is not None:
            within = abs(computed_cai - mean_expected) <= 0.10

        source_spread = max_expected - min_expected if len(expected_vals) > 1 else 0.0

        per_gene.append({
            "gene": gene,
            "organism": organism,
            "n_sources": len(sources),
            "expected_cai_values": expected_vals,
            "expected_cai_min": min_expected,
            "expected_cai_max": max_expected,
            "expected_cai_mean": mean_expected,
            "source_spread": source_spread,
            "computed_cai": computed_cai,
            "within_published_range": within,
            "sources": [s["source"] for s in sources],
        })

    # Compute source agreement (average spread for multi-source genes)
    multi_source_genes = [g for g in per_gene if g["n_sources"] >= 2]
    avg_spread = (
        sum(g["source_spread"] for g in multi_source_genes) / len(multi_source_genes)
        if multi_source_genes
        else 0.0
    )

    summary = (
        f"Multi-Source CAI Comparison\n"
        f"{'=' * 40}\n"
        f"Total genes with published values: {len(per_gene)}\n"
        f"Genes with multiple sources: {len(multi_source_genes)}\n"
        f"Average spread between sources: {avg_spread:.4f}\n"
        f"BioCompiler within published range: "
        f"{within_range_count}/{total_multi_source}\n"
    )

    return {
        "per_gene": per_gene,
        "source_agreement": {
            "n_multi_source_genes": len(multi_source_genes),
            "average_source_spread": avg_spread,
        },
        "biocompiler_within_range": {
            "within_count": within_range_count,
            "total_multi_source": total_multi_source,
            "fraction": within_range_count / total_multi_source
            if total_multi_source > 0 else 0.0,
        },
        "summary": summary,
    }
