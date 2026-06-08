"""
BioCompiler Type System — Codon Tables and Core Data
====================================================
Defines the standard genetic code, amino-acid mappings, BLOSUM62 matrix,
predicate class names, promoter consensus sequences, and supporting
data structures shared across the type_system package.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Tuple

from ..types import Verdict

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# Module-level named constants
# ────────────────────────────────────────────────────────────
# BLOSUM62 score returned when an amino-acid pair is not found
_BLOSUM62_MISSING_SCORE: int = -10

# Simplified nearest-neighbor dG coefficients (kcal/mol per base pair)
_DG_GC_PAIR_KCAL: float = -1.5
_DG_AU_PAIR_KCAL: float = -0.5
_DG_GU_PAIR_KCAL: float = -0.3

# Threshold multipliers for predicate verdict levels
_PROMOTER_UNCERTAIN_RATIO: float = 0.8
_TM_BORDERLINE_RATIO: float = 0.85
_MRNA_MODERATE_DG_RATIO: float = 0.7

# Codon ramp (co-translational folding) constants
_CODON_RAMP_LENGTH: int = 30
_MIN_RAMP_FOR_WARNING: int = 10
_FAST_CODON_CAI_THRESHOLD: float = 0.7
_PAUSE_SITE_CAI_THRESHOLD: float = 0.3
_HIGH_AVG_CAI_THRESHOLD: float = 0.9

# Co-translational folding structure confidence thresholds
_COTRANS_HIGH_CONFIDENCE: float = 0.7
_COTRANS_LOW_CONFIDENCE: float = 0.5

# MRNAStability organism-specific thresholds
_MRNA_STABILITY_THRESHOLDS: Dict[str, float] = {
    "E_coli": 0.8,
    "Escherichia_coli": 0.8,
    "Homo_sapiens": 0.7,
    "Mus_musculus": 0.7,
    "CHO_K1": 0.7,
    "Saccharomyces_cerevisiae": 0.75,
}

# TM domain organism-specific minimum hydrophobic stretch lengths
_TM_EUKARYOTIC_MIN_STRETCH: int = 19
_TM_PROKARYOTIC_MIN_STRETCH: int = 17

# mRNA secondary structure organism-specific dG cutoffs
_MRNA_DG_PROKARYOTE_FAIL: float = -15.0
_MRNA_DG_EUKARYOTE_FAIL: float = -25.0

# NoRestrictionSite minimum site length
_RESTRICTION_SITE_MIN_LENGTH: int = 6

# NoCpGIsland GC-rich relaxation threshold
_CPG_GC_RICH_THRESHOLD: float = 0.60
_CPG_DENSITY_MULTIPLIER: float = 2.0

# MaxEntScan insufficient-context sentinel score
_MAXENT_INSUFFICIENT_CONTEXT_SCORE: float = -50.0

# Eukaryotic promoter initiator search offsets (bp downstream of TATA box)
_EUK_INITIATOR_OFFSET_MIN: int = 20
_EUK_INITIATOR_OFFSET_MAX: int = 40

# mRNA instability motif thresholds
_INSTABILITY_T_RUN_MIN: int = 6

# Eukaryotic GT tolerance
_EUKARYOTE_GT_PER_BP: float = 1.0 / 50.0  # 1 GT per 50 bp


# ────────────────────────────────────────────────────────────
# Standard Genetic Code — CODON_TABLE (fixed: no invalid entries)
# ────────────────────────────────────────────────────────────
CODON_TABLE: Dict[str, str] = {
    # Phenylalanine
    "TTT": "F", "TTC": "F",
    # Leucine
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine (start)
    "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine
    "TAT": "Y", "TAC": "Y",
    # Histidine
    "CAT": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAT": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartic acid
    "GAT": "D", "GAC": "D",
    # Glutamic acid
    "GAA": "E", "GAG": "E",
    # Cysteine
    "TGT": "C", "TGC": "C",
    # Tryptophan
    "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Serine (AG- group)
    "AGT": "S", "AGC": "S",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stop codons
    "TAA": "*", "TAG": "*", "TGA": "*",
}

# Reverse lookup: amino acid -> list of codons
AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    if _aa not in AA_TO_CODONS:
        AA_TO_CODONS[_aa] = []
    AA_TO_CODONS[_aa].append(_codon)

# Derived sets for convenience
START_CODONS: List[str] = ["ATG"]
STOP_CODONS: List[str] = ["TAA", "TAG", "TGA"]


# ────────────────────────────────────────────────────────────
# BLOSUM62 Substitution Matrix (20x20 standard)
# Tuple-key format for backward compatibility: BLOSUM62[(aa1, aa2)] = score
# Data sourced from constants.py canonical definition.
# ────────────────────────────────────────────────────────────
from ..constants import BLOSUM62 as _BLOSUM62_NESTED  # noqa: E402

BLOSUM62: Dict[Tuple[str, str], int] = {}
for _a1, _row in _BLOSUM62_NESTED.items():
    for _a2, _score in _row.items():
        BLOSUM62[(_a1, _a2)] = _score


# ────────────────────────────────────────────────────────────
# Certificate levels
# ────────────────────────────────────────────────────────────
class CertLevel(Enum):
    GOLD = "GOLD"
    SILVER = "SILVER"
    BRONZE = "BRONZE"


# ────────────────────────────────────────────────────────────
# Splice verdict (dual-threshold)
# ────────────────────────────────────────────────────────────
class SpliceVerdict(Enum):
    PASS = auto()        # MaxEnt score < low threshold
    UNCERTAIN = auto()   # low <= score < high  (warn but don't block)
    FAIL = auto()        # score >= high threshold


# ────────────────────────────────────────────────────────────
# 12 Predicate Classes for Certified Optimization
# ────────────────────────────────────────────────────────────
PREDICATE_NAMES = [
    "NoStopCodons",           # 1 — no internal stops
    "NoCrypticSplice",        # 2 — dual-threshold splice check
    "NoCpGIsland",            # 3 — CpG island avoidance
    "NoRestrictionSite",      # 4 — enzyme site removal
    "NoGTDinucleotide",       # 5 — GT dinucleotide avoidance (cross-codon aware)
    "ValidCodingSeq",         # 6 — in-frame, valid codons only
    "ConservationScore",      # 7 — BLOSUM62-based AA conservation
    "CodonOptimality",        # 8 — CAI-based codon quality
    "NoCrypticPromoter",      # 9 — cryptic promoter avoidance
    "NoUnexpectedTMDomain",   # 10 — unexpected transmembrane domain detection
    "mRNASecondaryStructure", # 11 — mRNA secondary structure around RBS
    "CoTranslationalFolding", # 12 — co-translational folding pause-site preservation
    # Structure predicates
    "StructureConfidence",    # 13 — ESMFold structure quality confidence
    "NoMisfoldingRisk",       # 14 — misfolding risk indicators
    "CorrectFoldTopology",    # 15 — fold topology validation
    "NoUnexpectedInteraction",# 16 — unwanted protein-protein interactions
    # Stability predicates
    "StableFolding",          # 17 — thermodynamic stability (dG)
    "NoDestabilizingMutation",# 18 — no high-ddG mutations
    "DisulfideBondIntegrity", # 19 — cysteine pairing check
    "HydrophobicCoreQuality", # 20 — hydrophobic core composition
    # Solubility predicates
    "SolubleExpression",      # 21 — CamSol solubility score
    "NoAggregationProneRegion",#22 — aggregation-prone region detection
    "ChargeComposition",      # 23 — charge balance and pI
    "NoLongHydrophobicStretch",#24 — long hydrophobic stretch detection
    # Immunogenicity predicates
    "LowImmunogenicity",      # 25 — overall immunogenicity score
    "NoStrongTCellEpitope",   # 26 — MHC binding epitope detection
    "NoDominantBCellEpitope", # 27 — B-cell epitope coverage
    "PopulationCoverageSafe", # 28 — MHC allele population coverage
]


# ────────────────────────────────────────────────────────────
# Promoter consensus sequences
# ────────────────────────────────────────────────────────────
PROMOTER_CONSENSUS: Dict[str, Dict[str, object]] = {
    "E_coli": {
        "type": "prokaryotic",
        "sigma": "sigma70",
        "-35_box": "TTGACA",
        "-10_box": "TATAAT",
        "spacer": 17,  # optimal spacer length between -35 and -10 boxes
    },
    "eukaryote": {
        "type": "eukaryotic",
        "TATA_box": "TATAAA",
        "Initiator": "YYANWYY",  # IUPAC: Y=C/T, W=A/T, N=any
    },
}

# IUPAC ambiguity code expansion for promoter matching
_IUPAC_DNA = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def _match_iupac(base: str, pattern_base: str) -> bool:
    """Check if a DNA base matches an IUPAC ambiguity code."""
    return base in _IUPAC_DNA.get(pattern_base.upper(), {pattern_base.upper()})


def _score_consensus(seq_region: str, consensus: str) -> float:
    """Score a sequence region against a consensus pattern.

    Returns a value between 0.0 (no match) and 1.0 (perfect match).
    """
    if len(seq_region) != len(consensus):
        return 0.0
    matches = sum(1 for s, c in zip(seq_region, consensus) if _match_iupac(s, c))
    return matches / len(consensus)


@dataclass
class PredicateResult:
    """Result of checking one predicate against a sequence."""
    predicate: str
    passed: bool
    verdict: Optional[Verdict] = None  # used by NoCrypticSplice and others
    details: str = ""
    positions: List[int] = field(default_factory=list)
    verification_evidence: Optional[Dict[str, Any]] = None  # SLOT verification evidence
    mutagenesis_applied: bool = False
    unavoidable_constraints: List[str] = field(default_factory=list)


# Mapping from standard organism names to SPECIES dict keys
_ORGANISM_TO_SPECIES_KEY: Dict[str, str] = {
    "Homo_sapiens": "human",
    "human": "human",
    "Escherichia_coli": "ecoli",
    "E_coli": "ecoli",
    "ecoli": "ecoli",
    "Mus_musculus": "human",  # fallback — closest available
    "CHO_K1": "human",       # fallback — closest available
    "Saccharomyces_cerevisiae": "ecoli",  # fallback — closest available
}
