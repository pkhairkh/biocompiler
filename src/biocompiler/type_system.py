"""
BioCompiler Type System v9.2.0
==============================
Defines the core types, codon tables, BLOSUM62 matrix, and 28 predicate classes
for certified gene optimization: 12 DNA-level + 4 structure + 4 stability +
4 solubility + 4 immunogenicity.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, List, Dict, Set, Optional, Tuple
from .types import Verdict, SLOTMode, TypeCheckResult
from .sliding_gc import evaluate_sliding_gc, check_sliding_gc, SlidingGCResult, WindowViolation

# ── NUMBA integration for dinucleotide counting ───────────────────────
try:
    from .numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        fast_dinucleotide_count as _numba_fast_dinuc_count,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False
    _numba_fast_dinuc_count = None  # type: ignore[assignment]
    _seq_to_bytes = None  # type: ignore[assignment]


def _count_dinucs_fast(seq: str, *dinucleotides: str) -> tuple[int, ...]:
    """Count multiple dinucleotides in a single pass using the NUMBA kernel.

    Falls back to pure-Python counting when NUMBA is unavailable.

    Args:
        seq: DNA sequence string (uppercase ACGT).
        *dinucleotides: One or more dinucleotide strings (e.g. "GT", "CG", "AG").

    Returns:
        Tuple of counts, one per dinucleotide, in the same order as input.
    """
    n_dinucs = len(dinucleotides)
    if n_dinucs == 0:
        return ()

    # Fast path: NUMBA kernel
    if _HAS_NUMBA and _numba_fast_dinuc_count is not None:
        import numpy as _np
        seq_bytes = _seq_to_bytes(seq)
        dinuc_keys = _np.array(
            [[ord(d[0]), ord(d[1])] for d in dinucleotides],
            dtype=_np.uint8,
        )
        counts = _numba_fast_dinuc_count(seq_bytes, dinuc_keys, n_dinucs)
        return tuple(int(c) for c in counts)

    # Pure-Python fallback
    results = []
    for di in dinucleotides:
        count = 0
        pos = 0
        while True:
            pos = seq.find(di, pos)
            if pos == -1:
                break
            count += 1
            pos += 1
        results.append(count)
    return tuple(results)


logger = logging.getLogger(__name__)


__all__ = [
    # Core data tables
    "CODON_TABLE", "AA_TO_CODONS", "BLOSUM62",
    # Enums
    "CertLevel", "SpliceVerdict",
    # Predicate names
    "PREDICATE_NAMES",
    # Promoter data
    "PROMOTER_CONSENSUS",
    # Data classes
    "PredicateResult",
    # Low-level predicate checks
    "check_no_stop_codons", "check_no_cryptic_splice", "check_no_cpg_island",
    "check_no_restriction_site", "check_no_gt_dinucleotide", "check_no_avoidable_gt",
    "check_no_gt_dinucleotide_soft",
    "check_valid_coding_seq", "check_conservation_score", "check_codon_optimality",
    "check_no_cryptic_promoter", "check_no_unexpected_tm_domain",
    "check_mrna_secondary_structure", "check_co_translational_folding",
    "check_mrna_stability", "evaluate_mrna_stability",
    # High-level evaluate API
    "evaluate_gc_in_range", "evaluate_no_cryptic_splice", "evaluate_splice_correct",
    "evaluate_codon_adapted", "evaluate_no_restriction_site", "evaluate_in_frame",
    "evaluate_no_instability_motif", "evaluate_no_unexpected_tm_domain",
    "evaluate_mrna_secondary_structure", "evaluate_no_cryptic_promoter",
    "evaluate_no_cpg_island", "analyze_codon_at_position",
    "evaluate_co_translational_folding", "evaluate_all_predicates",
    "evaluate_no_stop_codons", "evaluate_no_gt_dinucleotide",
    "evaluate_valid_coding_seq", "evaluate_conservation_score",
    "evaluate_codon_optimality",
    # Sliding-window GC
    "evaluate_sliding_gc", "check_sliding_gc", "SlidingGCResult", "WindowViolation",
    # Cross-codon helpers
    "find_cross_codon_gt", "find_cross_codon_cg", "find_cross_codon_restriction",
    # Organism-aware helpers
    "_is_prokaryotic_organism", "_compute_max_gt_count", "_EUKARYOTE_GT_PER_BP",
    # Registry
    "PredicateRegistry", "registry",
]


# ────────────────────────────────────────────────────────────
# Module-level named constants
# ────────────────────────────────────────────────────────────
# BLOSUM62 score returned when an amino-acid pair is not found
_BLOSUM62_MISSING_SCORE: int = -10

# Simplified nearest-neighbor ΔG coefficients (kcal/mol per base pair)
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

# mRNA secondary structure organism-specific ΔG cutoffs
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

# ────────────────────────────────────────────────────────────
# BLOSUM62 Substitution Matrix (20x20 standard)
# Tuple-key format for backward compatibility: BLOSUM62[(aa1, aa2)] = score
# Data sourced from constants.py canonical definition.
# ────────────────────────────────────────────────────────────
from .constants import BLOSUM62 as _BLOSUM62_NESTED

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
    "StableFolding",          # 17 — thermodynamic stability (ΔG)
    "NoDestabilizingMutation",# 18 — no high-ΔΔG mutations
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


# ────────────────────────────────────────────────────────────
# 12 Predicate check functions
# ────────────────────────────────────────────────────────────

def check_no_stop_codons(seq: str) -> PredicateResult:
    """Predicate 1: No internal stop codons.

    The last codon in the reading frame is allowed to be a stop
    (natural termination). Only stops that appear BEFORE the last
    codon are flagged as violations.
    """
    if len(seq) < 3:
        return PredicateResult("NoStopCodons", True, verdict=Verdict.PASS, details="Sequence too short")
    last_codon_start = len(seq) - 3
    violations = []
    for i in range(0, last_codon_start, 3):  # skip the last codon
        codon = seq[i:i+3]
        if codon in ("TAA", "TAG", "TGA"):
            violations.append(i)
    if violations:
        return PredicateResult("NoStopCodons", False, verdict=Verdict.FAIL, details="Internal stop codons found", positions=violations)
    return PredicateResult("NoStopCodons", True, verdict=Verdict.PASS, details="No internal stop codons")


def check_no_cryptic_promoter(seq: str, organism: str = "E_coli", threshold: float = 0.7) -> PredicateResult:
    """Predicate 9: No cryptic promoter sites.

    Scans for promoter motifs using position weight matrix scoring.
    For prokaryotes (E_coli), scans for -35 (TTGACA) and -10 (TATAAT)
    boxes separated by a 17bp spacer. For eukaryotes, scans for TATA box
    (TATAAA) and Initiator (YYANWYY), PLUS requires an additional promoter
    element (CAAT box or GC box) within 50bp for FAIL.

    Key improvement: single promoter-like motifs (e.g., a lone TATA box)
    are ubiquitous in coding sequences and should NOT trigger FAIL. Only
    FAIL when MULTIPLE promoter elements are found together (within 50bp
    for eukaryotes, or both -35 and -10 for prokaryotes) AND a TATA box
    is present.

    Scoring is based on match quality (how many positions match the consensus).
    - Multiple promoter elements + TATA box with score >= threshold: FAIL
    - TATA box only (no additional elements): PASS with warning
    - Borderline match: UNCERTAIN
    - Otherwise: PASS
    """
    # Additional eukaryotic promoter elements
    _CAAT_BOX = "CCAAT"
    _GC_BOX = "GGGCGG"
    _PROMOTER_ELEMENT_WINDOW = 50  # bp window for multi-element detection

    seq = seq.upper()
    if len(seq) < 6:
        return PredicateResult("NoCrypticPromoter", True, verdict=Verdict.PASS,
                               details="Sequence too short for promoter motifs")

    consensus_info = PROMOTER_CONSENSUS.get(organism, PROMOTER_CONSENSUS["E_coli"])
    worst_score = 0.0
    worst_pos = -1
    worst_verdict = Verdict.PASS
    promoter_positions: List[int] = []

    if consensus_info["type"] == "prokaryotic":
        box35 = consensus_info["-35_box"]
        box10 = consensus_info["-10_box"]
        spacer = consensus_info["spacer"]
        # Total promoter length: len(-35) + spacer + len(-10)
        promoter_len = len(box35) + spacer + len(box10)

        for i in range(len(seq) - promoter_len + 1):
            region_35 = seq[i:i + len(box35)]
            region_10 = seq[i + len(box35) + spacer:i + promoter_len]

            if len(region_10) < len(box10):
                continue

            score_35 = _score_consensus(region_35, box35)
            score_10 = _score_consensus(region_10, box10)
            # Combined score: average of both boxes
            # Both boxes must individually have reasonable scores for a real promoter
            combined = (score_35 + score_10) / 2.0

            if combined > worst_score:
                worst_score = combined
                worst_pos = i
                promoter_positions = [i, i + len(box35) + spacer]

    elif consensus_info["type"] == "eukaryotic":
        tata_box = consensus_info["TATA_box"]
        initiator = consensus_info["Initiator"]

        # Scan for TATA box
        for i in range(len(seq) - len(tata_box) + 1):
            score_tata = _score_consensus(seq[i:i + len(tata_box)], tata_box)

            # Look for initiator within ~25-35bp downstream of TATA box start
            for offset in range(_EUK_INITIATOR_OFFSET_MIN, _EUK_INITIATOR_OFFSET_MAX):
                ini_start = i + offset
                if ini_start + len(initiator) > len(seq):
                    break
                score_ini = _score_consensus(seq[ini_start:ini_start + len(initiator)], initiator)
                combined = (score_tata + score_ini) / 2.0

                if combined > worst_score:
                    worst_score = combined
                    worst_pos = i
                    promoter_positions = [i, ini_start]

    # Determine verdict based on worst score with multi-element requirement
    if worst_score >= threshold:
        # Check if this is truly a multi-element promoter (not just a single motif)
        has_tata = False
        has_additional_element = False

        if consensus_info["type"] == "eukaryotic":
            # For eukaryotes: require TATA box + additional element (CAAT/GC box)
            # within 50bp of the worst promoter match
            search_start = max(0, worst_pos - _PROMOTER_ELEMENT_WINDOW)
            search_end = min(len(seq), worst_pos + _PROMOTER_ELEMENT_WINDOW)
            search_region = seq[search_start:search_end]

            # Check for TATA box in the region
            for j in range(len(search_region) - len("TATAAA") + 1):
                if _score_consensus(search_region[j:j + 6], "TATAAA") >= 0.8:
                    has_tata = True
                    break

            # Check for CAAT box
            caat_pos = search_region.find(_CAAT_BOX)
            if caat_pos >= 0:
                has_additional_element = True

            # Check for GC box
            gc_pos = search_region.find(_GC_BOX)
            if gc_pos >= 0:
                has_additional_element = True

            if has_tata and has_additional_element:
                worst_verdict = Verdict.FAIL
            elif has_tata:
                # TATA box alone — not enough for a cryptic promoter
                worst_verdict = Verdict.PASS
            else:
                # No clear TATA box — not a eukaryotic promoter
                worst_verdict = Verdict.PASS
        else:
            # For prokaryotes: -35 and -10 combo already requires two elements
            # This is inherently multi-element, so FAIL is appropriate
            worst_verdict = Verdict.FAIL
    elif worst_score >= threshold * _PROMOTER_UNCERTAIN_RATIO:
        worst_verdict = Verdict.UNCERTAIN
    else:
        worst_verdict = Verdict.PASS

    passed = worst_verdict != Verdict.FAIL
    details = f"Worst promoter score {worst_score:.3f} at pos {worst_pos}"
    if worst_verdict == Verdict.PASS:
        if worst_score >= threshold and consensus_info["type"] == "eukaryotic":
            details = (f"Promoter-like motif found (score {worst_score:.3f}) but lacks "
                       f"multiple elements — likely false positive")
        else:
            details = f"No significant promoter motifs found (worst score {worst_score:.3f})"

    return PredicateResult(
        "NoCrypticPromoter", passed, verdict=worst_verdict,
        details=details,
        positions=promoter_positions,
    )


def check_no_cryptic_splice(seq: str, low_thresh: float = 3.0, high_thresh: float = 6.0, organism: str = "") -> PredicateResult:
    """Predicate 2: No cryptic splice sites (dual-threshold PASS/UNCERTAIN/FAIL).

    Uses the proper MaxEntScan log-odds scoring model (Yeo & Burge 2004) from
    the maxentscan module to evaluate both donor (GT) and acceptor (AG) splice
    sites.  Numeric scores are converted to SpliceVerdict thresholds:
      - score < low_thresh  → PASS
      - low_thresh ≤ score < high_thresh → UNCERTAIN
      - score ≥ high_thresh → FAIL

    Organism-specific thresholds:
    - Prokaryotes: auto-PASS (no splicing in prokaryotes)
    - Eukaryotes: high_thresh=8.0 (stricter than default of 6.0 to reduce
      false positives from common coding-sequence GT/AG dinucleotides)

    Skipped for prokaryotic organisms (splice sites are a eukaryote-specific
    concern).
    """
    # Skip cryptic splice check for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        logger.info(
            "Cryptic splice check skipped for prokaryotic organism '%s'",
            organism,
        )
        return PredicateResult(
            "NoCrypticSplice", True, verdict=Verdict.PASS,
            details=f"Cryptic splice check skipped for prokaryotic organism '{organism}'",
        )

    # For eukaryotes, use stricter high_thresh of 8.0 to reduce false positives
    effective_high_thresh = high_thresh
    if organism and not _is_prokaryotic_organism(organism):
        effective_high_thresh = max(high_thresh, 8.0)

    from .maxentscan import score_donor, score_acceptor

    seq = seq.upper()
    max_score = _MAXENT_INSUFFICIENT_CONTEXT_SCORE
    worst_pos = -1
    worst_verdict = Verdict.PASS

    # Scan donor sites (GT dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "GT":
            score = score_donor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0
            if score < low_thresh:
                v = Verdict.PASS
            elif score < effective_high_thresh:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.FAIL
            if score > max_score:
                max_score = score
                worst_pos = i
                worst_verdict = v

    # Scan acceptor sites (AG dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "AG":
            score = score_acceptor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0
            if score < low_thresh:
                v = Verdict.PASS
            elif score < effective_high_thresh:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.FAIL
            if score > max_score:
                max_score = score
                worst_pos = i
                worst_verdict = v

    if worst_pos < 0:
        return PredicateResult("NoCrypticSplice", True, verdict=Verdict.PASS,
                               details="No splice dinucleotides found")

    passed = worst_verdict != Verdict.FAIL
    return PredicateResult("NoCrypticSplice", passed, verdict=worst_verdict,
                           details=f"Worst splice score {max_score:.2f} at pos {worst_pos}",
                           positions=[worst_pos] if worst_pos >= 0 else [])


def _is_prokaryotic_organism(organism: str) -> bool:
    """Return True if the organism is prokaryotic.

    Uses :func:`biocompiler.organism_config.is_eukaryotic_organism` when
    available; falls back to a simple name-based heuristic for common
    prokaryotic identifiers.

    Args:
        organism: Organism name (e.g. ``"E_coli"``, ``"Homo_sapiens"``).

    Returns:
        True if the organism is prokaryotic, False otherwise.
    """
    if not organism:
        return False
    try:
        from .organism_config import is_eukaryotic_organism
        return not is_eukaryotic_organism(organism)
    except Exception:
        # Fallback: common prokaryotic identifiers
        prokaryotic_names = {
            "E_coli", "E_coli_K12", "E_coli_BL21",
            "Escherichia_coli", "ecoli",
            "Bacillus_subtilis", "bsub",
            "Pseudomonas_aeruginosa",
        }
        return organism in prokaryotic_names


def check_no_cpg_island(seq: str, window: int = 200, threshold: float = 0.6, organism: str = "") -> PredicateResult:
    """Predicate 3: No CpG islands (Obs/Exp CG ratio > threshold in any window).

    CpG island avoidance is primarily relevant for mammalian expression
    systems where CpG methylation can lead to gene silencing.  For
    prokaryotic organisms (e.g. *E. coli*), CpG islands have no known
    regulatory significance, so the check is skipped when a prokaryotic
    organism is specified.

    Optimized: Uses sliding window with O(1) updates per step instead
    of O(W) full window scan, reducing total complexity from O(N*W) to O(N).

    Args:
        seq: DNA sequence to evaluate.
        window: Sliding window size in nucleotides (default 200).
        threshold: Maximum allowed Obs/Exp CG ratio (default 0.6).
        organism: Target organism name.  If prokaryotic, the check is
            skipped and PASS is returned immediately.  If empty (default),
            the check runs as before for backward compatibility.

    Returns:
        PredicateResult with PASS/FAIL verdict.
    """
    # Skip CpG checking for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        logger.info(
            "CpG island check skipped for prokaryotic organism '%s'",
            organism,
        )
        return PredicateResult(
            "NoCpGIsland", True, verdict=Verdict.PASS,
            details=f"CpG island check skipped for prokaryotic organism '{organism}'",
        )

    n = len(seq)
    if n < window:
        return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                               details=f"Sequence length {n} < window size {window}")

    # Fast short-circuit: if total CG count is 0, no CpG island is possible
    total_cg = _count_dinucs_fast(seq, "CG")[0]
    if total_cg == 0:
        return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                               details="No CG dinucleotides found in sequence")

    # Pre-compute CG positions for fast lookup
    # cg_at[i] = 1 if seq[i:i+2] == "CG", else 0
    cg_at = [0] * (n - 1)
    for i in range(n - 1):
        if seq[i] == 'C' and seq[i + 1] == 'G':
            cg_at[i] = 1

    # Initialize first window
    c_count = seq[:window].count("C")
    g_count = seq[:window].count("G")
    cg_count = sum(cg_at[:window - 1])

    worst_ratio = 0.0
    worst_start = -1

    # Check first window
    expected = (c_count * g_count) / window if window > 0 else 0
    obs_exp = cg_count / expected if expected > 0 else 0.0
    if obs_exp > worst_ratio:
        worst_ratio = obs_exp
        worst_start = 0

    # Slide window — O(1) per step
    for start in range(1, n - window + 1):
        # Remove outgoing base at start-1, add incoming base at start+window-1
        outgoing = seq[start - 1]
        incoming = seq[start + window - 1]

        if outgoing == 'C':
            c_count -= 1
        elif outgoing == 'G':
            g_count -= 1
        if incoming == 'C':
            c_count += 1
        elif incoming == 'G':
            g_count += 1

        # Update CG count: remove cg_at[start-1], add cg_at[start+window-2]
        cg_count -= cg_at[start - 1]
        if start + window - 2 < n - 1:
            cg_count += cg_at[start + window - 2]

        expected = (c_count * g_count) / window if window > 0 else 0
        obs_exp = cg_count / expected if expected > 0 else 0.0
        if obs_exp > worst_ratio:
            worst_ratio = obs_exp
            worst_start = start

    if worst_ratio > threshold:
        # GC-content-aware relaxation: for GC-rich sequences, CpG density
        # is naturally higher. Only FAIL if CpG density is >2x the expected
        # density for the sequence's GC content.
        seq_gc = (seq.count("G") + seq.count("C")) / len(seq) if len(seq) > 0 else 0
        if seq_gc > _CPG_GC_RICH_THRESHOLD:
            # For GC-rich sequences: expected CpG density = GC_fraction^2
            # The obs_exp ratio already normalizes by C*G/window, but for
            # very GC-rich sequences, even the "expected" can be too low.
            # Use a relaxed threshold: 2x the normal threshold.
            relaxed_threshold = threshold * _CPG_DENSITY_MULTIPLIER
            if worst_ratio > relaxed_threshold:
                return PredicateResult("NoCpGIsland", False, verdict=Verdict.FAIL,
                                       details=(f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f} "
                                                f"> {relaxed_threshold:.3f} (GC-rich={seq_gc:.1%}, "
                                                f"relaxed threshold)"),
                                       positions=[worst_start])
            return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                                   details=(f"CpG Obs/Exp={worst_ratio:.3f} within relaxed threshold "
                                            f"{relaxed_threshold:.3f} for GC-rich sequence ({seq_gc:.1%})"))
        return PredicateResult("NoCpGIsland", False, verdict=Verdict.FAIL,
                               details=f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f} > {threshold}",
                               positions=[worst_start])
    return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                           details=f"Worst CpG Obs/Exp ratio {worst_ratio:.3f} <= {threshold}")


def check_no_restriction_site(seq: str, enzymes: List[str], min_site_length: int = _RESTRICTION_SITE_MIN_LENGTH) -> PredicateResult:
    """Predicate 4: No restriction enzyme recognition sites.

    Only checks restriction sites that are >= min_site_length bp (default 6bp).
    Short 4bp restriction sites are too common in coding sequences and
    cause excessive false positives.

    Also handles cross-codon sites that span 3+ codons properly.
    """
    from .restriction_sites import get_recognition_site
    violations = []
    for enzyme in enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        # Skip short restriction sites (< min_site_length bp)
        if len(site) < min_site_length:
            continue
        pos = seq.find(site)
        while pos != -1:
            violations.append(pos)
            pos = seq.find(site, pos + 1)
    if violations:
        return PredicateResult("NoRestrictionSite", False, verdict=Verdict.FAIL,
                               details=f"Restriction sites found at {violations}",
                               positions=violations)
    return PredicateResult("NoRestrictionSite", True, verdict=Verdict.PASS, details="No restriction sites found")


def check_no_gt_dinucleotide(seq: str) -> PredicateResult:
    """Predicate 5: No GT dinucleotides (5' splice donor mimic), including cross-codon.

    This is the STRICT version — any GT fails the predicate.

    Uses the NUMBA ``fast_dinucleotide_count`` kernel when available
    for fast count-based short-circuit; falls back to pure-Python otherwise.
    """
    # Fast short-circuit: if GT count is 0, no need to enumerate positions
    gt_count = _count_dinucs_fast(seq, "GT")[0]
    if gt_count == 0:
        return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS, details="No GT dinucleotides found")
    # Need positions for the result — enumerate only when count > 0
    positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
    return PredicateResult("NoGTDinucleotide", False, verdict=Verdict.FAIL,
                           details=f"GT dinucleotides at {positions}",
                           positions=positions)


def check_no_avoidable_gt(seq: str, organism: str = "") -> PredicateResult:
    """Predicate 5 (relaxed): No avoidable GT dinucleotides.

    A GT is "unavoidable" if ALL synonymous codons for that amino acid
    also contain GT or create a cross-codon GT.  This predicate PASSES
    if every remaining GT in the sequence is unavoidable — i.e., there
    is no synonymous substitution that could remove it.

    Specifically:
    - Within-codon GT: unavoidable if every synonymous codon for the AA
      also contains "GT" (e.g., Valine GTN where all 4 codons start with GT)
    - Cross-codon GT: unavoidable if no combination of synonymous codons
      for the two adjacent AAs eliminates the boundary GT

    Skipped for prokaryotic organisms (GT splice donor sites are a
    eukaryote-specific concern).
    """
    # Skip GT dinucleotide check for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        logger.info(
            "GT dinucleotide check skipped for prokaryotic organism '%s'",
            organism,
        )
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details=f"GT dinucleotide check skipped for prokaryotic organism '{organism}'",
        )

    # Fast short-circuit: if GT count is 0, no need to enumerate positions
    gt_count = _count_dinucs_fast(seq, "GT")[0]
    if gt_count == 0:
        return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS, details="No GT dinucleotides found")

    gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]

    avoidable_positions = []
    unavoidable_positions = []

    for pos in gt_positions:
        codon_idx = pos // 3  # which codon does position 'pos' fall in?
        codon_start = codon_idx * 3
        next_codon_start = codon_start + 3

        # Determine whether this GT is within a single codon or crosses a boundary
        if pos + 1 < next_codon_start:
            # Within-codon GT (both bases in the same codon)
            codon = seq[codon_start:codon_start + 3]
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                unavoidable_positions.append(pos)
                continue

            # Check if any synonymous codon avoids GT
            has_avoidable = False
            for alt in AA_TO_CODONS.get(aa, []):
                if "GT" not in alt:
                    # Also check this alt doesn't create cross-codon GT
                    # with the previous codon's last base
                    if codon_start > 0:
                        prev_base = seq[codon_start - 1]
                        if prev_base + alt[0] == "GT":
                            continue  # would create cross-codon GT
                    # And check it doesn't create GT with the next codon's first base
                    if next_codon_start + 3 <= len(seq):
                        next_base = seq[next_codon_start]
                        if alt[-1] + next_base == "GT":
                            continue  # would create cross-codon GT
                    has_avoidable = True
                    break

            if has_avoidable:
                avoidable_positions.append(pos)
            else:
                unavoidable_positions.append(pos)
        else:
            # Cross-codon GT (pos is last base of one codon, pos+1 is first of next)
            prev_codon_start = codon_start  # codon containing 'pos'
            curr_codon_start = next_codon_start  # codon containing 'pos+1'

            if curr_codon_start + 3 > len(seq):
                unavoidable_positions.append(pos)
                continue

            prev_codon = seq[prev_codon_start:prev_codon_start + 3]
            curr_codon = seq[curr_codon_start:curr_codon_start + 3]
            prev_aa = CODON_TABLE.get(prev_codon)
            curr_aa = CODON_TABLE.get(curr_codon)

            if prev_aa is None or curr_aa is None:
                unavoidable_positions.append(pos)
                continue

            # If one side is a stop codon, we can only try changing the other side
            if prev_aa == "*" and curr_aa == "*":
                unavoidable_positions.append(pos)
                continue

            has_avoidable = False

            if prev_aa == "*":
                # Can only change the current codon — check if any alt
                # avoids the boundary GT (don't require alt to be GT-free internally)
                for c_alt in AA_TO_CODONS.get(curr_aa, [curr_codon]):
                    if prev_codon[-1] + c_alt[0] != "GT":
                        has_avoidable = True
                        break
            elif curr_aa == "*":
                # Can only change the previous codon — check if any alt
                # avoids the boundary GT (don't require alt to be GT-free internally)
                for p_alt in AA_TO_CODONS.get(prev_aa, [prev_codon]):
                    if p_alt[-1] + curr_codon[0] != "GT":
                        has_avoidable = True
                        break
            else:
                # Both are regular AAs — try all combinations
                # For cross-codon GTs, only check if the specific boundary
                # GT can be eliminated; internal GTs in alternative codons
                # are separate positions that get checked independently
                prev_alts = AA_TO_CODONS.get(prev_aa, [prev_codon])
                curr_alts = AA_TO_CODONS.get(curr_aa, [curr_codon])

                for p_alt in prev_alts:
                    for c_alt in curr_alts:
                        if p_alt[-1] + c_alt[0] != "GT":
                            has_avoidable = True
                            break
                    if has_avoidable:
                        break

            if has_avoidable:
                avoidable_positions.append(pos)
            else:
                unavoidable_positions.append(pos)

    if avoidable_positions:
        return PredicateResult("NoGTDinucleotide", False, verdict=Verdict.FAIL,
                               details=(f"Avoidable GT dinucleotides at {avoidable_positions}; "
                                        f"unavoidable at {unavoidable_positions}"),
                               positions=avoidable_positions)
    return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS,
                           details=(f"All {len(unavoidable_positions)} GT dinucleotides are "
                                    f"unavoidable (no synonymous substitution can remove them)"),
                           positions=unavoidable_positions)


# ────────────────────────────────────────────────────────────
# Soft GT dinucleotide check for eukaryotes
# ────────────────────────────────────────────────────────────

# Default: eukaryotes can tolerate ~1 GT per 50 bp of coding sequence.
# This is biologically justified because in-codon GTs from optimal codons
# (GGT, TGT, GTT, etc.) are common in high-expression eukaryotic genes and
# only become problematic when they form strong cryptic splice donor sites
# (MaxEntScan score >= threshold).
_EUKARYOTE_GT_PER_BP: float = 1.0 / 50.0  # 1 GT per 50 bp


def _compute_max_gt_count(seq_len: int, organism: str = "") -> int:
    """Compute the maximum allowed GT dinucleotide count for a sequence.

    For prokaryotes: 0 (hard constraint — GT avoidance is irrelevant since
    prokaryotes have no spliceosome, but if a user explicitly requests
    GT checking on a prokaryote, any GT is a FAIL).

    For eukaryotes: ``max(1, int(seq_len * _EUKARYOTE_GT_PER_BP))``.
    This accounts for the biological reality that in-codon GTs from optimal
    codons are common and acceptable in eukaryotic CDS.

    Args:
        seq_len: Length of the DNA sequence in base pairs.
        organism: Target organism name. Used to determine prokaryotic vs
            eukaryotic classification.

    Returns:
        Maximum allowed GT count (0 for prokaryotes, length-scaled for
        eukaryotes).
    """
    if organism and _is_prokaryotic_organism(organism):
        return 0
    return max(1, int(seq_len * _EUKARYOTE_GT_PER_BP))


def check_no_gt_dinucleotide_soft(
    seq: str,
    organism: str = "",
    max_gt_count: int | None = None,
) -> PredicateResult:
    """Predicate 5 (soft): No GT dinucleotides with eukaryote-aware tolerance.

    This is the organism-aware soft-constraint version of the GT dinucleotide
    check, designed for eukaryotic gene optimization where destroying CAI to
    eliminate every GT is counter-productive.

    Evaluation semantics:

    - **Prokaryotes**: Hard constraint — any GT dinucleotide is a FAIL.
      (Prokaryotes have no spliceosome, so GT dinucleotides are
      biologically irrelevant; if a user explicitly checks GT for a
      prokaryote, we treat it as a hard constraint.)

    - **Eukaryotes**: Soft constraint — GTs are reported but the predicate
      uses ``LIKELY_FAIL`` (not ``FAIL``) when GTs exceed ``max_gt_count``,
      indicating a soft violation that should not block the optimization.
      The predicate PASSES (``PASS``) if GT count ≤ ``max_gt_count``.

    - ``max_gt_count``: If not provided, auto-computed from sequence length
      using :func:`_compute_max_gt_count` (default: 1 GT per 50 bp for
      eukaryotes, 0 for prokaryotes). Can be explicitly set to 0 for
      hard-constraint behavior on eukaryotes.

    The result distinguishes in-codon vs cross-codon GTs in the details,
    since in-codon GTs from optimal codons (GGT for Gly, TGT for Cys,
    GTT/GTC/GTA/GTG for Val) are biologically acceptable for eukaryotes.

    Args:
        seq: DNA sequence to evaluate.
        organism: Target organism name. If prokaryotic, any GT is FAIL.
        max_gt_count: Maximum GT count before triggering SOFT_FAIL for
            eukaryotes. If None, auto-computed from sequence length.

    Returns:
        PredicateResult with verdict:
        - PASS: No GT dinucleotides (or GT count ≤ max_gt_count for eukaryotes)
        - LIKELY_FAIL: GT count > max_gt_count for eukaryotes (soft fail)
        - FAIL: Any GT for prokaryotes (hard constraint)
    """
    seq = seq.upper()
    # Fast short-circuit: use NUMBA kernel for count, skip position enumeration if 0
    gt_count_fast = _count_dinucs_fast(seq, "GT")[0]
    if gt_count_fast == 0:
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details="No GT dinucleotides found",
        )

    gt_positions = [i for i in range(len(seq) - 1) if seq[i:i + 2] == "GT"]

    if not gt_positions:
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details="No GT dinucleotides found",
        )

    # Compute max_gt_count if not provided
    if max_gt_count is None:
        max_gt_count = _compute_max_gt_count(len(seq), organism)

    # Count in-codon vs cross-codon GTs for reporting
    in_codon_gt = []
    cross_codon_gt = []
    for pos in gt_positions:
        codon_of_g = pos // 3
        codon_of_t = (pos + 1) // 3
        if codon_of_g == codon_of_t:
            in_codon_gt.append(pos)
        else:
            cross_codon_gt.append(pos)

    gt_count = len(gt_positions)

    # Prokaryotes: hard constraint (FAIL for any GT)
    if organism and _is_prokaryotic_organism(organism):
        return PredicateResult(
            "NoGTDinucleotide", False, verdict=Verdict.FAIL,
            details=(
                f"GT dinucleotides: {gt_count} "
                f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
                f"Hard constraint for prokaryotes: max_gt_count=0."
            ),
            positions=gt_positions,
        )

    # Eukaryotes: soft constraint
    if gt_count <= max_gt_count:
        return PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details=(
                f"GT dinucleotides: {gt_count} ≤ max_gt_count={max_gt_count} "
                f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
                f"Acceptable for eukaryotes — in-codon GTs from optimal codons "
                f"are biologically common."
            ),
            positions=gt_positions,
        )

    # GT count exceeds tolerance — soft fail (LIKELY_FAIL, not FAIL)
    # This indicates a warning, not a hard block. For soft constraints,
    # passed=True so the predicate doesn't appear in failed_predicates
    # (which would trigger strict mode errors). The LIKELY_FAIL verdict
    # and details still convey the soft violation to users who check.
    return PredicateResult(
        "NoGTDinucleotide", True, verdict=Verdict.LIKELY_FAIL,
        details=(
            f"GT dinucleotides: {gt_count} > max_gt_count={max_gt_count} "
            f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
            f"Soft constraint for eukaryotes: in-codon GTs from optimal codons "
            f"are acceptable (CAI > GT avoidance). Consider if these GTs form "
            f"strong cryptic splice donors (MaxEntScan score ≥ threshold)."
        ),
        positions=gt_positions,
    )


def check_valid_coding_seq(seq: str) -> PredicateResult:
    """Predicate 6: Valid coding sequence (length divisible by 3, all valid codons)."""
    if len(seq) % 3 != 0:
        return PredicateResult("ValidCodingSeq", False, verdict=Verdict.FAIL,
                               details=f"Sequence length {len(seq)} not divisible by 3")
    invalid = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon not in CODON_TABLE:
            invalid.append((i, codon))
    if invalid:
        return PredicateResult("ValidCodingSeq", False, verdict=Verdict.FAIL,
                               details=f"Invalid codons: {invalid}")
    return PredicateResult("ValidCodingSeq", True, verdict=Verdict.PASS, details="All codons valid")


def check_conservation_score(dna, protein, min_score: int = 0) -> PredicateResult:
    """Predicate 7: BLOSUM62 conservation score between DNA-derived and target protein.

    Translates the DNA sequence to its amino-acid sequence, then compares each
    position against the target protein using the BLOSUM62 substitution matrix.
    After a correct optimization the two sequences are identical, so every
    diagonal score is positive and the predicate should PASS.

    Backward-compatible: if both *dna* and *protein* are single amino-acid
    characters, the old two-AA substitution check is performed instead.

    Args:
        dna: Optimized DNA coding sequence (translated internally).
            For backward compatibility, also accepts a single amino-acid
            character (paired with *protein* as a single amino-acid).
        protein: Target protein sequence (amino-acid string).
            For backward compatibility, also accepts a single amino-acid
            character (paired with *dna* as a single amino-acid).
        min_score: Minimum BLOSUM62 score per position for PASS (default 0).
    """
    # Backward compatibility: old callers passed two single AA characters
    # as (original_aa, new_aa). Detect this pattern and handle it.
    if isinstance(dna, str) and isinstance(protein, str) and len(dna) == 1 and len(protein) == 1:
        score = BLOSUM62.get((dna, protein), _BLOSUM62_MISSING_SCORE)
        passed = score >= min_score
        return PredicateResult(
            "ConservationScore", passed,
            verdict=Verdict.PASS if passed else Verdict.FAIL,
            details=f"BLOSUM62({dna},{protein})={score}, min={min_score}",
        )

    # Translate DNA → protein
    translated = ""
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            return PredicateResult(
                "ConservationScore", False,
                verdict=Verdict.FAIL,
                details=f"Invalid/stop codon '{codon}' at DNA position {i}",
            )
        translated += aa

    # Length mismatch is an automatic FAIL
    if len(translated) != len(protein):
        return PredicateResult(
            "ConservationScore", False,
            verdict=Verdict.FAIL,
            details=f"Length mismatch: translated {len(translated)} AA vs target {len(protein)} AA",
        )

    # Score each position with BLOSUM62
    total_score = 0
    min_found = 0
    for t_aa, p_aa in zip(translated, protein):
        s = BLOSUM62.get((t_aa, p_aa), _BLOSUM62_MISSING_SCORE)
        total_score += s
        if s < min_found:
            min_found = s

    passed = min_found >= min_score
    return PredicateResult(
        "ConservationScore", passed,
        verdict=Verdict.PASS if passed else Verdict.FAIL,
        details=f"BLOSUM62 total={total_score}, min_pos={min_found}, min={min_score} (translated={translated}, target={protein})",
    )


def check_codon_optimality(dna, organism, min_cai: float = 0.0) -> PredicateResult:
    """Predicate 8: Codon optimality (CAI score above threshold).

    Looks up the codon adaptiveness table for *organism* from
    CODON_ADAPTIVENESS_TABLES and computes the geometric-mean CAI across
    all codons in *dna*.

    Backward-compatible: if *organism* is a ``dict`` (old-style
    ``species_cai`` mapping), it is used directly as the CAI weight table
    and *dna* is treated as a single codon string (old calling convention).

    Args:
        dna: DNA coding sequence to evaluate.
            For backward compatibility, also accepts a single codon string
            when *organism* is a dict.
        organism: Target organism name (e.g. ``"e_coli"``, ``"Homo_sapiens"``).
            For backward compatibility, also accepts a ``Dict[str, float]``
            of codon→CAI weights (old ``species_cai`` parameter).
        min_cai: Minimum acceptable CAI for PASS (default 0.0).
    """
    # Backward compatibility: old callers passed (codon, species_cai_dict, min_cai)
    if isinstance(organism, dict):
        species_cai = organism
        cai = species_cai.get(dna, 0.0)
        passed = cai >= min_cai
        return PredicateResult(
            "CodonOptimality", passed,
            verdict=Verdict.PASS if passed else Verdict.FAIL,
            details=f"CAI({dna})={cai:.4f}, min={min_cai}",
        )

    # New-style: (dna_sequence, organism_name, min_cai)
    from .organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

    # Resolve organism name to its canonical key
    canonical = resolve_organism(organism)
    species_cai: Dict[str, float] = CODON_ADAPTIVENESS_TABLES.get(canonical, {})
    if not species_cai:
        # Fallback to E. coli if organism not found
        species_cai = CODON_ADAPTIVENESS_TABLES.get("Escherichia_coli", {})

    dna = dna.upper()
    num_codons = len(dna) // 3

    if num_codons == 0:
        return PredicateResult(
            "CodonOptimality", True, verdict=Verdict.PASS,
            details="Sequence too short for CAI computation",
        )

    # Compute geometric-mean CAI (Sharp & Li 1987)
    import math
    log_product = 0.0
    for i in range(num_codons):
        codon = dna[i * 3:(i + 1) * 3]
        w = species_cai.get(codon, 0.0)
        if w <= 0.0:
            log_product += math.log(1e-4)  # clamp to avoid log(0)
        else:
            log_product += math.log(w)
    cai = math.exp(log_product / num_codons)

    passed = cai >= min_cai
    return PredicateResult(
        "CodonOptimality", passed,
        verdict=Verdict.PASS if passed else Verdict.FAIL,
        details=f"CAI={cai:.4f}, min={min_cai}, organism={canonical}",
    )


def check_no_unexpected_tm_domain(
    seq: str,
    is_cytosolic: bool = True,
    window_size: int = 19,
    threshold: float = 0.68,
    organism: str = "",
) -> PredicateResult:
    """Predicate 10: No unexpected transmembrane (TM) domains after mutagenesis.

    If a cytosolic protein gains hydrophobic stretches from amino acid
    substitutions, that constitutes a FAIL. Transmembrane domains are
    detected by sliding a window of `window_size` amino acids and computing
    the fraction of hydrophobic residues (A, V, I, L, M, F, W, Y).

    Organism-aware window sizing:
    - Eukaryotes: minimum hydrophobic stretch of 19 aa (default)
    - Prokaryotes: minimum hydrophobic stretch of 17 aa

    Flanking charge check: Only flag as FAIL if the hydrophobic stretch
    also has appropriate flanking charges consistent with a true TM domain
    (positive-inside rule). Short hydrophobic stretches in soluble proteins
    that lack flanking positive charges are NOT flagged as FAIL.

    Verdict logic (only applies when is_cytosolic=True):
    - If any window exceeds `threshold` AND has TM-like flanking charges: FAIL
    - If any window exceeds `threshold` but lacks TM flanking charges: UNCERTAIN
    - If any window exceeds `threshold * _TM_BORDERLINE_RATIO`: UNCERTAIN
    - Otherwise: PASS

    If is_cytosolic=False (membrane protein), TM domains are expected: PASS.

    Args:
        seq: DNA sequence to evaluate.
        is_cytosolic: Whether the protein is cytosolic (default True).
        window_size: Sliding window size in amino acids (default 19).
            Overridden by organism-specific minimums if organism is provided.
        threshold: Hydrophobic fraction threshold for FAIL (default 0.68).
        organism: Target organism name. If provided, used to adjust window
            size (eukaryotes: 19, prokaryotes: 17).

    Returns:
        PredicateResult with PASS/UNCERTAIN/FAIL verdict.
    """
    seq = seq.upper()

    if not is_cytosolic:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.PASS,
            details="Membrane protein — TM domains are expected",
        )

    # Adjust window size based on organism
    if organism:
        if _is_prokaryotic_organism(organism):
            window_size = max(window_size, _TM_PROKARYOTIC_MIN_STRETCH)
        else:
            window_size = max(window_size, _TM_EUKARYOTIC_MIN_STRETCH)

    # Translate DNA to amino acids
    aa_seq = _translate_dna_to_aa(seq)

    if len(aa_seq) < window_size:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.PASS,
            details=f"Protein too short for TM domain scan ({len(aa_seq)} aa < {window_size} window)",
        )

    HYDROPHOBIC = set("AVILMFWY")
    POSITIVE = set("KR")
    worst_frac = 0.0
    worst_pos = -1
    borderline_positions: List[int] = []
    fail_positions: List[int] = []
    fail_no_flank_positions: List[int] = []

    for i in range(len(aa_seq) - window_size + 1):
        window = aa_seq[i:i + window_size]
        hydro_count = sum(1 for aa in window if aa in HYDROPHOBIC)
        frac = hydro_count / window_size
        if frac > worst_frac:
            worst_frac = frac
            worst_pos = i
        if frac > threshold * _TM_BORDERLINE_RATIO:
            borderline_positions.append(i)
        if frac > threshold:
            # Check flanking charges (positive-inside rule)
            # True TM domains have positive charges (K/R) on the cytoplasmic side
            # Look for positive residues in the 5 aa flanking each end
            flank_n = aa_seq[max(0, i - 5):i]
            flank_c = aa_seq[i + window_size:min(len(aa_seq), i + window_size + 5)]
            n_pos_count = sum(1 for aa in flank_n if aa in POSITIVE)
            c_pos_count = sum(1 for aa in flank_c if aa in POSITIVE)
            # A true TM domain typically has positive charges on one side
            # (cytoplasmic side) but NOT both sides. Signal peptides are
            # N-terminal with positive residues before the hydrophobic core.
            has_tm_flanking = (n_pos_count >= 1 or c_pos_count >= 1) and not (n_pos_count >= 2 and c_pos_count >= 2)
            if has_tm_flanking:
                fail_positions.append(i)
            else:
                # Hydrophobic stretch without TM-like flanking — likely
                # a soluble protein hydrophobic patch, not a true TM domain
                fail_no_flank_positions.append(i)

    if fail_positions:
        return PredicateResult(
            "NoUnexpectedTMDomain", False, verdict=Verdict.FAIL,
            details=(f"TM domain detected: worst hydrophobic fraction {worst_frac:.3f} "
                     f"at AA pos {worst_pos} exceeds threshold {threshold} "
                     f"with TM-like flanking charges "
                     f"({len(fail_positions)} window(s) failing)"),
            positions=fail_positions,
        )

    # If hydrophobic stretches exist but lack TM flanking charges → LIKELY_PASS
    # (not UNCERTAIN), since the absence of flanking charges is strong evidence
    # that this is NOT a true TM domain but rather a soluble protein hydrophobic patch.
    if fail_no_flank_positions:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.LIKELY_PASS,
            details=(f"Hydrophobic stretch without TM flanking charges: "
                     f"worst fraction {worst_frac:.3f} at AA pos {worst_pos} "
                     f"({len(fail_no_flank_positions)} window(s) — likely soluble protein patch, not TM domain)"),
            positions=fail_no_flank_positions,
        )

    # Borderline hydrophobic fraction (between threshold*0.85 and threshold)
    # without flanking charges is also likely a soluble patch → LIKELY_PASS
    if borderline_positions:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.LIKELY_PASS,
            details=(f"Borderline hydrophobic stretch: worst fraction {worst_frac:.3f} "
                     f"at AA pos {worst_pos} exceeds {threshold * _TM_BORDERLINE_RATIO:.3f} "
                     f"({len(borderline_positions)} window(s) — below TM threshold, likely soluble patch)"),
            positions=borderline_positions,
        )

    return PredicateResult(
        "NoUnexpectedTMDomain", True, verdict=Verdict.PASS,
        details=f"No TM domain detected (worst hydrophobic fraction {worst_frac:.3f})",
    )


def check_mrna_secondary_structure(
    seq: str,
    window_start: int = 0,
    window_end: int = 50,
    dg_threshold: float = -15.0,
    use_viennarna: bool = True,
    organism: str = "",
) -> PredicateResult:
    """Predicate 11: No strong mRNA secondary structure around RBS/start codon.

    Checks for stable secondary structure near the ribosome binding site
    that could block ribosome binding. When ViennaRNA or the Nussinov
    fallback is available, uses proper thermodynamic folding (Turner
    nearest-neighbor parameters). Otherwise falls back to the simplified
    hairpin model for backward compatibility.

    Organism-specific ΔG cutoffs:
    - Prokaryotes: ΔG < -15 kcal/mol is FAIL (RBS is critical for
      ribosome binding)
    - Eukaryotes: ΔG < -25 kcal/mol is FAIL (cap-dependent translation
      is less sensitive to RBS secondary structure)
    - If organism is not specified, uses the provided dg_threshold

    Args:
        seq: DNA sequence to evaluate.
        window_start: Start position of the analysis window (default 0).
        window_end: End position of the analysis window (default 50).
        dg_threshold: ΔG threshold for FAIL (default -15.0 kcal/mol).
            Overridden by organism-specific cutoffs if organism is provided.
        use_viennarna: If True (default), try ViennaRNA/Nussinov for
            accurate ΔG computation. If False, use the legacy toy model.
        organism: Target organism name. If provided, uses organism-specific
            ΔG cutoffs instead of dg_threshold.

    Returns:
        PredicateResult with PASS/UNCERTAIN/FAIL verdict.
    """
    seq = seq.upper()

    # Apply organism-specific ΔG threshold
    effective_threshold = dg_threshold
    if organism:
        if _is_prokaryotic_organism(organism):
            effective_threshold = _MRNA_DG_PROKARYOTE_FAIL
        else:
            effective_threshold = _MRNA_DG_EUKARYOTE_FAIL

    # Try ViennaRNA/Nussinov for accurate ΔG
    if use_viennarna:
        try:
            from .viennarna import is_viennarna_available, compute_5prime_dg
            if is_viennarna_available():
                dg = compute_5prime_dg(seq, window=window_end - window_start)
                if dg <= effective_threshold:
                    return PredicateResult(
                        "mRNASecondaryStructure", False, verdict=Verdict.FAIL,
                        details=(f"Strong mRNA secondary structure (ViennaRNA): "
                                 f"ΔG={dg:.1f} kcal/mol <= {effective_threshold}"),
                    )
                if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
                    return PredicateResult(
                        "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
                        details=(f"Moderate mRNA secondary structure (ViennaRNA): "
                                 f"ΔG={dg:.1f} kcal/mol <= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f}"),
                    )
                return PredicateResult(
                    "mRNASecondaryStructure", True, verdict=Verdict.PASS,
                    details=(f"Weak mRNA secondary structure (ViennaRNA): "
                             f"ΔG={dg:.1f} kcal/mol"),
                )
        except Exception:
            logger.debug("Falling through to Nussinov fallback", exc_info=True)

        try:
            from .viennarna_fallback import compute_approx_dg
            dg = compute_approx_dg(seq, region="5utr")
            if dg <= effective_threshold:
                return PredicateResult(
                    "mRNASecondaryStructure", False, verdict=Verdict.FAIL,
                    details=(f"Strong mRNA secondary structure (Nussinov fallback): "
                             f"ΔG≈{dg:.1f} kcal/mol <= {effective_threshold}"),
                )
            if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
                return PredicateResult(
                    "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
                    details=(f"Moderate mRNA secondary structure (Nussinov fallback): "
                             f"ΔG≈{dg:.1f} kcal/mol <= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f}"),
                )
            return PredicateResult(
                "mRNASecondaryStructure", True, verdict=Verdict.PASS,
                details=(f"Weak mRNA secondary structure (Nussinov fallback): "
                         f"ΔG≈{dg:.1f} kcal/mol"),
            )
        except Exception:
            logger.debug("Falling through to legacy toy model", exc_info=True)

    # Legacy toy model (original implementation, backward-compatible)
    # Uses a simplified hairpin approximation for ΔG
    seq = seq.upper()
    # Extract the window around the RBS/start codon
    effective_end = min(window_end, len(seq))
    window_seq = seq[window_start:effective_end]

    if len(window_seq) < 4:
        return PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.PASS,
            details="Sequence window too short for structure analysis",
        )

    # Convert DNA to RNA for pairing analysis (T -> U)
    rna = window_seq.replace("T", "U")

    # Count potential base pairs using a simplified complementary pairing
    # We look at each position and its potential complement in the sequence
    # A simplified approach: count pairs in the first half pairing with
    # the reversed second half (mimicking a hairpin stem)
    gc_pairs = 0
    au_pairs = 0
    gu_pairs = 0

    half = len(rna) // 2
    first_half = rna[:half]
    second_half = rna[half:2 * half]  # mirror region

    for i in range(min(len(first_half), len(second_half))):
        # Pair first_half[i] with second_half reversed
        j = len(second_half) - 1 - i
        if j < 0:
            break
        base_5 = first_half[i]
        base_3 = second_half[j]

        if (base_5 == "G" and base_3 == "C") or (base_5 == "C" and base_3 == "G"):
            gc_pairs += 1
        elif (base_5 == "A" and base_3 == "U") or (base_5 == "U" and base_3 == "A"):
            au_pairs += 1
        elif (base_5 == "G" and base_3 == "U") or (base_5 == "U" and base_3 == "G"):
            gu_pairs += 1

    # Simplified nearest-neighbor ΔG estimate
    dg = _DG_GC_PAIR_KCAL * gc_pairs + _DG_AU_PAIR_KCAL * au_pairs + _DG_GU_PAIR_KCAL * gu_pairs

    if dg <= effective_threshold:
        return PredicateResult(
            "mRNASecondaryStructure", False, verdict=Verdict.FAIL,
            details=(f"Strong mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
                     f"<= {effective_threshold} (GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
        )

    if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
        # Moderate ΔG: use GC content of the window to refine verdict.
        # AT-rich windows have weaker actual structures → LIKELY_PASS.
        # GC-rich windows can form stronger structures → LIKELY_FAIL.
        window_gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) if window_seq else 0.5
        if window_gc < 0.5:
            return PredicateResult(
                "mRNASecondaryStructure", True, verdict=Verdict.LIKELY_PASS,
                details=(f"Moderate mRNA secondary structure: ΔG={dg:.1f} kcal/mol, "
                         f"but AT-rich window (GC={window_gc:.0%}) weakens structure "
                         f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
            )
        return PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.LIKELY_FAIL,
            details=(f"Moderate mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
                     f"<= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f} "
                     f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
        )

    return PredicateResult(
        "mRNASecondaryStructure", True, verdict=Verdict.PASS,
        details=(f"Weak mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
                 f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
    )


# ────────────────────────────────────────────────────────────
# Predicate 12: Co-translational folding helpers
# ────────────────────────────────────────────────────────────

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


def _resolve_species_cai(key: str) -> Dict[str, float]:
    """Resolve an organism name or SPECIES key to a flat codon→CAI-weight dict.

    Uses CODON_ADAPTIVENESS_TABLES with resolve_organism() as the
    single source of truth for CAI weights.

    Args:
        key: Organism name, alias, or short species key
            (e.g. ``"ecoli"``, ``"Homo_sapiens"``, ``"human"``).

    Returns:
        Dict mapping codon strings to CAI adaptiveness values.
    """
    from .organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

    organism = resolve_organism(key)
    if organism in CODON_ADAPTIVENESS_TABLES:
        return dict(CODON_ADAPTIVENESS_TABLES[organism])
    # Fallback to ecoli for unknown organisms
    return dict(CODON_ADAPTIVENESS_TABLES["Escherichia_coli"])


def _compute_codon_ramp_score(seq: str, species_cai: Dict[str, float]) -> Dict[str, Any]:
    """Compute codon ramp score and identify pause sites / speed disruptions.

    The codon ramp is the first ~30 codons where slow codons are beneficial
    for proper ribosome loading. Outside the ramp, slow codons (CAI < 0.3)
    serve as pause sites that allow co-translational folding of protein
    domains before downstream sequence emerges from the ribosome.

    Args:
        seq: DNA coding sequence (uppercase, length divisible by 3).
        species_cai: Dict mapping codon strings to their CAI values
            (relative adaptiveness, 0.0–1.0).

    Returns:
        Dict with keys:
          - ramp_score: average CAI in the first 30 codons (lower = more
            ramp = better for ribosome loading).
          - pause_sites: list of (codon_position, CAI) tuples where
            CAI < 0.3 outside the ramp region (positions 30+).
          - speed_disruptions: list of (codon_position, original_slow_cai,
            new_fast_cai) tuples where a likely pause site was replaced by
            a fast codon during optimization (heuristic: CAI in 0.7–1.0
            where the slowest synonymous codon would have CAI < 0.3).
    """
    seq = seq.upper()
    num_codons = len(seq) // 3
    ramp_length = min(_CODON_RAMP_LENGTH, num_codons)  # first 30 codons = ramp region

    # Collect per-codon CAI values
    codon_cais: List[Tuple[int, float]] = []  # (codon_index, cai)
    for i in range(num_codons):
        codon = seq[i * 3:(i + 1) * 3]
        cai = species_cai.get(codon, 0.0)
        codon_cais.append((i, cai))

    # Ramp score: average CAI in first 30 codons
    ramp_cais = [cai for idx, cai in codon_cais[:ramp_length]]
    ramp_score = sum(ramp_cais) / len(ramp_cais) if ramp_cais else 0.0

    # Pause sites: slow codons outside the ramp
    pause_sites: List[Tuple[int, float]] = []
    for idx, cai in codon_cais[ramp_length:]:
        if cai < _PAUSE_SITE_CAI_THRESHOLD:
            pause_sites.append((idx, cai))

    # Speed disruptions: fast codons that likely replaced slow ones
    # A speed disruption is detected when the current codon is fast (CAI > 0.7)
    # but the slowest synonymous codon for the same AA has CAI < 0.3,
    # suggesting an optimization replaced a natural pause site.
    speed_disruptions: List[Tuple[int, float, float]] = []
    for idx, cai in codon_cais:
        if idx < ramp_length:
            continue  # ramp region — speed-ups are expected there
        if cai <= _FAST_CODON_CAI_THRESHOLD:
            continue  # not a fast codon
        codon = seq[idx * 3:(idx + 1) * 3]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            continue
        # Check if the slowest synonymous codon for this AA is a pause site
        syn_codons = AA_TO_CODONS.get(aa, [codon])
        slowest_cai = min(species_cai.get(c, 0.0) for c in syn_codons)
        if slowest_cai < _PAUSE_SITE_CAI_THRESHOLD:
            speed_disruptions.append((idx, slowest_cai, cai))

    return {
        "ramp_score": ramp_score,
        "pause_sites": pause_sites,
        "speed_disruptions": speed_disruptions,
    }


def check_co_translational_folding(
    seq: str,
    species_cai: Dict[str, float],
    domain_boundaries: List[int] | None = None,
    min_pause_cai: float = 0.3,
    structure_confidence: float | None = None,
    plddt_score: float | None = None,
) -> PredicateResult:
    """Predicate 12: Check co-translational folding preservation.

    Checks whether codon optimization has disrupted critical pause sites
    that are important for proper protein folding during translation.
    Slow codons (low tRNA adaptation) at domain boundaries allow the
    nascent chain to fold properly before downstream sequence emerges.

    When structure_confidence (ESMFold confidence or similar) is provided,
    UNCERTAIN verdicts can be resolved:
    - structure_confidence > 0.7 and pLDDT is good → resolve to PASS
    - structure_confidence < 0.5 → resolve to FAIL
    - Otherwise keep UNCERTAIN

    Args:
        seq: DNA coding sequence (uppercase).
        species_cai: Dict mapping codon strings to CAI values.
        domain_boundaries: Optional list of codon positions (0-based)
            where protein domains start/end. If provided, these positions
            are checked for speed-up (CAI > 0.7 where a pause is needed).
        min_pause_cai: CAI threshold below which a codon is considered
            a pause site (default 0.3).
        structure_confidence: Optional ESMFold confidence score (0–1).
            When provided, UNCERTAIN verdicts can be resolved using
            structural evidence.
        plddt_score: Optional pLDDT score from structure prediction.
            Used alongside structure_confidence for verdict resolution.

    Returns:
        PredicateResult with PASS/UNCERTAIN/FAIL verdict.
    """
    seq = seq.upper()
    num_codons = len(seq) // 3

    if num_codons == 0:
        return PredicateResult(
            "CoTranslationalFolding", True, verdict=Verdict.PASS,
            details="Sequence too short for co-translational folding analysis",
        )

    ramp_info = _compute_codon_ramp_score(seq, species_cai)
    ramp_score = ramp_info["ramp_score"]
    pause_sites = ramp_info["pause_sites"]
    speed_disruptions = ramp_info["speed_disruptions"]

    # Positions to report in result
    flagged_positions: List[int] = []
    details_parts: List[str] = []

    # --- Check ramp region ---
    ramp_length = min(_CODON_RAMP_LENGTH, num_codons)
    ramp_all_fast = all(
        species_cai.get(seq[i * 3:(i + 1) * 3], 0.0) > _FAST_CODON_CAI_THRESHOLD
        for i in range(ramp_length)
    )

    if ramp_all_fast and ramp_length >= _MIN_RAMP_FOR_WARNING:
        details_parts.append(
            f"Ramp region (first {ramp_length} codons) is entirely fast "
            f"(avg CAI={ramp_score:.3f}) — ribosome jam risk"
        )
        flagged_positions.extend(range(ramp_length))

    # --- Check domain boundaries ---
    domain_disrupted = 0
    if domain_boundaries:
        for boundary_pos in domain_boundaries:
            if boundary_pos < 0 or boundary_pos >= num_codons:
                continue
            codon = seq[boundary_pos * 3:(boundary_pos + 1) * 3]
            boundary_cai = species_cai.get(codon, 0.0)
            if boundary_cai > _FAST_CODON_CAI_THRESHOLD:
                # Fast codon at domain boundary — pause site was disrupted
                domain_disrupted += 1
                flagged_positions.append(boundary_pos)
                details_parts.append(
                    f"Domain boundary at codon {boundary_pos} has CAI={boundary_cai:.3f} "
                    f"(fast codon where pause is needed)"
                )
    else:
        # Use heuristics: proline-rich regions and hydrophobicity changes
        # suggest domain boundaries. Flag UNCERTAIN if the average CAI
        # is very high (>0.9) throughout (no pause sites at all).
        avg_cai = sum(
            species_cai.get(seq[i * 3:(i + 1) * 3], 0.0)
            for i in range(num_codons)
        ) / num_codons if num_codons > 0 else 0.0

        if avg_cai > _HIGH_AVG_CAI_THRESHOLD and not pause_sites:
            details_parts.append(
                f"Average CAI={avg_cai:.3f} with no pause sites detected — "
                f"co-translational folding may be disrupted (no domain boundaries provided)"
            )

    # --- Determine verdict ---
    if domain_boundaries and domain_disrupted > 0 and ramp_all_fast:
        # Ramp destroyed AND domain boundaries disrupted
        verdict = Verdict.FAIL
        passed = False
    elif domain_boundaries and domain_disrupted >= 2:
        # Multiple domain boundaries disrupted
        verdict = Verdict.LIKELY_FAIL
        passed = False
    elif domain_boundaries and domain_disrupted == 1:
        # Single domain boundary disrupted — LIKELY_FAIL (not UNCERTAIN)
        # since a single disrupted boundary is meaningful evidence
        verdict = Verdict.LIKELY_FAIL
        passed = True
    elif ramp_all_fast and ramp_length >= _MIN_RAMP_FOR_WARNING:
        # Ramp too fast — but if average CAI is moderate (not extremely high),
        # ribosome jam risk is reduced → LIKELY_PASS instead of UNCERTAIN.
        # Very high CAI (> 0.95) is more concerning → LIKELY_FAIL.
        if ramp_score > 0.95:
            verdict = Verdict.LIKELY_FAIL
        else:
            verdict = Verdict.LIKELY_PASS
        passed = True
    elif speed_disruptions:
        # Some pause sites may have been replaced by fast codons
        # Use disruption rate to determine severity
        num_codons_total = len(seq) // 3
        disruption_rate = len(speed_disruptions) / max(num_codons_total, 1)
        if len(speed_disruptions) <= 2 or disruption_rate < 0.05:
            verdict = Verdict.LIKELY_PASS
        else:
            verdict = Verdict.LIKELY_FAIL
        passed = True
        details_parts.append(
            f"{len(speed_disruptions)} potential pause site(s) replaced by fast codons"
        )
    else:
        verdict = Verdict.PASS
        passed = True
        if pause_sites:
            details_parts.append(
                f"Good ramp (avg CAI={ramp_score:.3f}) with {len(pause_sites)} "
                f"natural pause site(s) preserved"
            )
        else:
            details_parts.append(
                f"Ramp avg CAI={ramp_score:.3f}, no pause site concerns detected"
            )

    # --- Resolve UNCERTAIN using structural evidence ---
    if verdict == Verdict.UNCERTAIN and structure_confidence is not None:
        if structure_confidence > _COTRANS_HIGH_CONFIDENCE:
            # High ESMFold confidence and good pLDDT → structure is likely correct
            # despite codon changes; co-translational folding preserved
            verdict = Verdict.PASS
            passed = True
            details_parts.append(
                f"UNCERTAIN resolved to PASS: structure_confidence={structure_confidence:.3f} > {_COTRANS_HIGH_CONFIDENCE}"
            )
        elif structure_confidence < _COTRANS_LOW_CONFIDENCE:
            # Low confidence suggests the structure is wrong → folding disrupted
            verdict = Verdict.LIKELY_FAIL
            passed = False
            details_parts.append(
                f"UNCERTAIN resolved to LIKELY_FAIL: structure_confidence={structure_confidence:.3f} < {_COTRANS_LOW_CONFIDENCE}"
            )
        # else: keep UNCERTAIN

    details = "; ".join(details_parts) if details_parts else "Co-translational folding appears preserved"

    return PredicateResult(
        "CoTranslationalFolding",
        passed,
        verdict=verdict,
        details=details,
        positions=flagged_positions,
    )


# ────────────────────────────────────────────────────────────
# Predicate: MRNAStability — CAI-weighted codon stability
# ────────────────────────────────────────────────────────────

def check_mrna_stability(
    seq: str,
    organism: str = "Homo_sapiens",
    threshold: float | None = None,
) -> PredicateResult:
    """Check mRNA stability using CAI-weighted codon stability scores.

    Computes a composite stability score based on the mrna_stability module's
    motif scanning combined with CAI-weighted codon optimality. Uses
    organism-specific thresholds to determine PASS/UNCERTAIN/FAIL.

    Organism-specific stability thresholds:
    - E. coli: 0.8
    - Human: 0.7
    - Yeast: 0.75

    Args:
        seq: DNA coding sequence.
        organism: Target organism name.
        threshold: Override threshold for PASS. If None, uses organism default.

    Returns:
        PredicateResult with PASS/UNCERTAIN/FAIL verdict.
    """
    from .mrna_stability import score_mrna_stability

    seq = seq.upper()

    if len(seq) < 3:
        return PredicateResult(
            "MRNAStability", True, verdict=Verdict.PASS,
            details="Sequence too short for mRNA stability analysis",
        )

    # Determine threshold
    if threshold is None:
        threshold = _MRNA_STABILITY_THRESHOLDS.get(organism, 0.7)

    # Get motif-based stability score
    report = score_mrna_stability(seq, organism)
    stability_score = report.overall_score

    # Also compute CAI-weighted codon stability contribution
    species_cai = _resolve_species_cai(organism)
    num_codons = len(seq) // 3
    if num_codons > 0:
        cai_values = [
            species_cai.get(seq[i * 3:(i + 1) * 3], 0.0)
            for i in range(num_codons)
        ]
        avg_cai = sum(cai_values) / num_codons
        # Blend motif score with CAI: weighted average (70% motif, 30% CAI)
        combined_score = 0.7 * stability_score + 0.3 * avg_cai
    else:
        avg_cai = 0.0
        combined_score = stability_score

    # Determine verdict
    if combined_score >= threshold:
        verdict = Verdict.PASS
        passed = True
        details = (
            f"mRNA stability score {combined_score:.3f} >= {threshold:.3f} "
            f"(motif={stability_score:.3f}, avg_CAI={avg_cai:.3f}, "
            f"stabilizing={report.stabilizing_count}, "
            f"destabilizing={report.destabilizing_count})"
        )
    elif combined_score >= threshold * 0.85:
        # Borderline stability: use CAI quality to refine verdict.
        # If CAI is good (>= 0.5), the codon usage supports stability even
        # if motifs are borderline → LIKELY_PASS. Otherwise → LIKELY_FAIL.
        if avg_cai >= 0.5:
            verdict = Verdict.LIKELY_PASS
            passed = True
            details = (
                f"mRNA stability score {combined_score:.3f} borderline "
                f"(threshold={threshold:.3f}), but good CAI={avg_cai:.3f} supports stability "
                f"(motif={stability_score:.3f}, "
                f"stabilizing={report.stabilizing_count}, "
                f"destabilizing={report.destabilizing_count})"
            )
        else:
            verdict = Verdict.LIKELY_FAIL
            passed = True
            details = (
                f"mRNA stability score {combined_score:.3f} borderline "
                f"(threshold={threshold:.3f}) with weak CAI={avg_cai:.3f} "
                f"(motif={stability_score:.3f}, "
                f"stabilizing={report.stabilizing_count}, "
                f"destabilizing={report.destabilizing_count})"
            )
    else:
        verdict = Verdict.FAIL
        passed = False
        details = (
            f"mRNA stability score {combined_score:.3f} < {threshold:.3f} "
            f"(motif={stability_score:.3f}, avg_CAI={avg_cai:.3f}, "
            f"stabilizing={report.stabilizing_count}, "
            f"destabilizing={report.destabilizing_count})"
        )

    return PredicateResult(
        "MRNAStability", passed, verdict=verdict,
        details=details,
    )


def evaluate_mrna_stability(
    seq: str,
    organism: str = "Homo_sapiens",
    threshold: float | None = None,
) -> TypeCheckResult:
    """Evaluate mRNA stability for a coding sequence.

    Uses CAI-weighted codon stability scores with organism-specific
    thresholds. The stability score blends motif scanning (70%) with
    codon adaptation index (30%).

    Args:
        seq: DNA coding sequence.
        organism: Target organism name.
        threshold: Override threshold. If None, uses organism default
            (E. coli: 0.8, Human: 0.7, Yeast: 0.75).

    Returns:
        TypeCheckResult with PASS/UNCERTAIN/FAIL verdict.
    """
    result = check_mrna_stability(seq, organism, threshold)

    if result.verdict == Verdict.PASS:
        return TypeCheckResult(
            predicate=f"MRNAStability({organism})",
            verdict=Verdict.LIKELY_PASS if result.details.startswith("mRNA stability score") else Verdict.PASS,
        )
    elif result.verdict == Verdict.UNCERTAIN:
        return TypeCheckResult(
            predicate=f"MRNAStability({organism})",
            verdict=Verdict.UNCERTAIN,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate=f"MRNAStability({organism})",
            verdict=Verdict.FAIL,
            violation=result.details,
        )


# ────────────────────────────────────────────────────────────
# Helper: DNA-to-amino-acid translation
# ────────────────────────────────────────────────────────────

def _translate_dna_to_aa(seq: str) -> str:
    """Translate a DNA sequence to an amino acid string using the CODON_TABLE.

    Reads the sequence in-frame from position 0. Codons not found in
    CODON_TABLE are silently skipped.

    Args:
        seq: DNA sequence (uppercase or mixed case).

    Returns:
        Amino acid string (single-letter codes, stops as '*').
    """
    seq = seq.upper()
    aa_list: List[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is not None:
            aa_list.append(aa)
    return "".join(aa_list)


# ────────────────────────────────────────────────────────────
# Cross-codon constraint helpers
# ────────────────────────────────────────────────────────────

def find_cross_codon_gt(seq: str) -> List[int]:
    """Find GT dinucleotides that span codon boundaries (pos i-1,i where i%3==0)."""
    positions = []
    for i in range(3, len(seq) - 1):
        if i % 3 == 0 and seq[i-1] == "G" and seq[i] == "T":
            positions.append(i - 1)
    return positions


def find_cross_codon_cg(seq: str) -> List[int]:
    """Find CG dinucleotides that span codon boundaries."""
    positions = []
    for i in range(3, len(seq) - 1):
        if i % 3 == 0 and seq[i-1] == "C" and seq[i] == "G":
            positions.append(i - 1)
    return positions


def find_cross_codon_restriction(seq: str, site: str) -> List[int]:
    """Find restriction sites that span codon boundaries (3+ codons).

    A site spans 3+ codons when it covers more than 9 bp (3 codons),
    meaning it must cross at least 2 codon boundaries.
    """
    positions = []
    site_len = len(site)
    for i in range(len(seq) - site_len + 1):
        if seq[i:i+site_len] == site:
            codon_start_i = (i // 3) * 3
            codon_end_i = ((i + site_len - 1) // 3) * 3 + 3
            # A site spanning 3+ codons covers more than 9 bp
            if codon_end_i - codon_start_i > 9:
                positions.append(i)
    return positions


# ════════════════════════════════════════════════════════════════
# High-level evaluate_* API — returns TypeCheckResult objects
# ════════════════════════════════════════════════════════════════



def evaluate_gc_in_range(seq: str, gc_lo: float = 0.30, gc_hi: float = 0.70) -> TypeCheckResult:
    """Evaluate whether GC content falls within the specified range.

    Args:
        seq: DNA sequence to evaluate.
        gc_lo: Minimum acceptable GC fraction (inclusive).
        gc_hi: Maximum acceptable GC fraction (inclusive).

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    if not seq:
        return TypeCheckResult(
            predicate=f"GCInRange({gc_lo}, {gc_hi})",
            verdict=Verdict.UNCERTAIN,
            violation="Empty sequence",
        )
    seq = seq.upper()
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    if gc_lo <= gc <= gc_hi:
        return TypeCheckResult(
            predicate=f"GCInRange({gc_lo}, {gc_hi})",
            verdict=Verdict.PASS,
        )
    direction = "below" if gc < gc_lo else "above"
    return TypeCheckResult(
        predicate=f"GCInRange({gc_lo}, {gc_hi})",
        verdict=Verdict.FAIL,
        violation=f"GC content {gc:.3f} is {direction} range [{gc_lo}, {gc_hi}]",
    )


def evaluate_no_cryptic_splice(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    cryptic_threshold: float = 3.0,
    uncertain_lo: float = 0.0,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate whether the sequence contains cryptic splice sites.

    Uses dual-threshold MaxEntScan scoring: sites scoring >= cryptic_threshold
    are FAIL, sites scoring >= uncertain_lo (but < cryptic_threshold) are
    UNCERTAIN, and sites below uncertain_lo are PASS.

    When uncertain_lo=0 (default), only PASS/FAIL verdicts are produced,
    preserving backward compatibility.

    Organism-specific thresholds:
    - Prokaryotes: auto-PASS (no splicing in prokaryotes).
    - Eukaryotes: use cryptic_threshold=8.0 as the effective FAIL
      threshold (unless overridden by caller), matching the optimizer's
      splice-site elimination logic.

    Args:
        seq: DNA sequence to evaluate.
        boundaries: Exon boundaries (used for context; currently informational).
        cryptic_threshold: Score threshold above which a site is FAIL.
        uncertain_lo: Score threshold above which a site is UNCERTAIN.
            Set to 0 to disable UNCERTAIN zone (binary PASS/FAIL only).
        organism: Target organism.  If prokaryotic, auto-PASS; if eukaryotic
            and cryptic_threshold < 8.0, raises the effective threshold to 8.0
            to match the optimizer's behaviour.

    Returns:
        TypeCheckResult with PASS/UNCERTAIN/FAIL verdict.
    """
    # Skip for prokaryotic organisms
    if organism and _is_prokaryotic_organism(organism):
        return TypeCheckResult(
            predicate="NoCrypticSplice",
            verdict=Verdict.PASS,
        )

    # For eukaryotes, use a higher default threshold matching the optimizer
    effective_threshold = cryptic_threshold
    if organism and not _is_prokaryotic_organism(organism):
        effective_threshold = max(cryptic_threshold, 8.0)

    from .maxentscan import score_donor, score_acceptor

    seq = seq.upper()
    worst_score = _MAXENT_INSUFFICIENT_CONTEXT_SCORE
    worst_pos = -1
    worst_verdict = Verdict.PASS

    # Scan donor sites (GT dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "GT":
            score = score_donor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0

            if score >= effective_threshold:
                v = Verdict.FAIL
            elif uncertain_lo > 0 and score >= uncertain_lo:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.PASS

            if score > worst_score:
                worst_score = score
                worst_pos = i
                worst_verdict = v

    # Scan acceptor sites (AG dinucleotides)
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "AG":
            score = score_acceptor(seq, i)
            if score <= _MAXENT_INSUFFICIENT_CONTEXT_SCORE:
                score = 0.0

            if score >= effective_threshold:
                v = Verdict.FAIL
            elif uncertain_lo > 0 and score >= uncertain_lo:
                v = Verdict.UNCERTAIN
            else:
                v = Verdict.PASS

            if score > worst_score:
                worst_score = score
                worst_pos = i
                worst_verdict = v

    return TypeCheckResult(
        predicate="NoCrypticSplice",
        verdict=worst_verdict,
        violation=(
            f"Worst cryptic splice score {worst_score:.2f} at pos {worst_pos}"
            if worst_verdict != Verdict.PASS else None
        ),
    )


def evaluate_splice_correct(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    cellular_context: str = "HEK293T",
) -> TypeCheckResult:
    """Evaluate splice correctness — whether known exon boundaries are respected.

    This is an informational predicate that checks for canonical splice signals
    (GT..AG) at intron boundaries. It passes if all introns have canonical
    splice signals, and fails if any intron lacks them.

    Args:
        seq: Full pre-mRNA sequence.
        boundaries: List of (start, end) tuples for each exon.
            If ``None`` or an empty list, there are no splice boundaries to
            check and the predicate returns PASS immediately.
        cellular_context: Cell type context (currently informational).

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    seq = seq.upper()

    # Guard: no boundaries or empty boundaries → no splice correction needed
    if boundaries is None or not boundaries or len(boundaries) < 2:
        # Single-exon gene or no boundaries provided — nothing to check
        return TypeCheckResult(
            predicate="SpliceCorrect",
            verdict=Verdict.PASS,
        )

    # Check canonical splice signals at each intron boundary
    n = len(seq)
    for i in range(len(boundaries) - 1):
        try:
            intron_start = boundaries[i][1]
            intron_end = boundaries[i + 1][0]
        except (IndexError, TypeError):
            # Malformed boundary entry — skip
            continue

        if intron_start >= intron_end:
            continue

        # Bounds-check before indexing the sequence string
        # Check donor (GT) at intron start
        if intron_start + 2 <= n and intron_start >= 0:
            donor = seq[intron_start:intron_start + 2]
            if donor != "GT":
                return TypeCheckResult(
                    predicate="SpliceCorrect",
                    verdict=Verdict.FAIL,
                    violation=f"Non-canonical donor {donor} at pos {intron_start}",
                )

        # Check acceptor (AG) at intron end
        if intron_end - 2 >= 0 and intron_end <= n:
            acceptor = seq[intron_end - 2:intron_end]
            if acceptor != "AG":
                return TypeCheckResult(
                    predicate="SpliceCorrect",
                    verdict=Verdict.FAIL,
                    violation=f"Non-canonical acceptor {acceptor} at pos {intron_end - 2}",
                )

    return TypeCheckResult(
        predicate="SpliceCorrect",
        verdict=Verdict.PASS,
    )


def evaluate_codon_adapted(
    seq: str,
    organism: str = "Homo_sapiens",
    threshold: float = 0.5,
) -> TypeCheckResult:
    """Evaluate whether codon adaptation index meets the threshold.

    Args:
        seq: DNA coding sequence.
        organism: Target organism for CAI computation.
        threshold: Minimum CAI score for PASS.

    Returns:
        TypeCheckResult with PASS/FAIL verdict.

    Raises:
        UnsupportedOrganismError: if organism is not in the supported set.
    """
    from .translation import compute_cai

    cai = compute_cai(seq, organism)
    if cai >= threshold:
        return TypeCheckResult(
            predicate=f"CodonAdapted({organism}, {threshold})",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate=f"CodonAdapted({organism}, {threshold})",
        verdict=Verdict.FAIL,
        violation=f"CAI {cai:.4f} is below threshold {threshold}",
    )


def evaluate_no_restriction_site(
    seq: str,
    enzymes: List[str] | None = None,
    enzyme_set: List[str] | None = None,
    min_site_length: int = _RESTRICTION_SITE_MIN_LENGTH,
) -> TypeCheckResult:
    """Evaluate whether the sequence contains restriction enzyme recognition sites.

    Only checks sites >= min_site_length bp (default 6bp) to avoid
    false positives from short 4bp sites that are common in coding sequences.

    Args:
        seq: DNA sequence to evaluate.
        enzymes: List of enzyme names to check for.
        enzyme_set: Alias for enzymes (used by certificate verification).
        min_site_length: Minimum recognition site length to check (default 6).

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    from .restriction_sites import get_recognition_site

    effective_enzymes = enzymes or enzyme_set or []
    seq = seq.upper()
    violations = []

    for enzyme in effective_enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        # Skip short restriction sites (< min_site_length bp)
        if len(site) < min_site_length:
            continue
        # Check for IUPAC patterns
        has_iupac = any(b not in "ACGT" for b in site.upper())
        if has_iupac:
            # Use IUPAC matching
            for i in range(len(seq) - len(site) + 1):
                window = seq[i:i + len(site)]
                match = True
                for s_base, p_base in zip(window, site.upper()):
                    from .constants import IUPAC_EXPAND
                    allowed = IUPAC_EXPAND.get(p_base, p_base)
                    if s_base not in allowed:
                        match = False
                        break
                if match:
                    violations.append((i, enzyme))
        else:
            pos = seq.find(site)
            while pos != -1:
                violations.append((pos, enzyme))
                pos = seq.find(site, pos + 1)

    if violations:
        details = ", ".join(f"{e}@{p}" for p, e in violations)
        return TypeCheckResult(
            predicate="NoRestrictionSite",
            verdict=Verdict.FAIL,
            violation=f"Restriction sites found: {details}",
        )
    return TypeCheckResult(
        predicate="NoRestrictionSite",
        verdict=Verdict.PASS,
    )


def evaluate_in_frame(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    exon_boundaries: List[Tuple[int, int]] | None = None,
) -> TypeCheckResult:
    """Evaluate whether the coding sequence is in-frame (valid ORF).

    Checks that each exon's length is divisible by 3 and that all codons
    in the coding regions are valid.

    Args:
        seq: DNA sequence to evaluate.
        boundaries: List of (start, end) tuples for each exon.
        exon_boundaries: Alias for boundaries (used by certificate verification).

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    effective_boundaries = boundaries or exon_boundaries or [(0, len(seq))]
    seq = seq.upper()

    for start, end in effective_boundaries:
        exon_seq = seq[start:end]
        if len(exon_seq) % 3 != 0:
            return TypeCheckResult(
                predicate="InFrame",
                verdict=Verdict.FAIL,
                violation=f"Exon [{start}, {end}) length {len(exon_seq)} not divisible by 3",
            )
        # Check all codons are valid
        for i in range(0, len(exon_seq), 3):
            codon = exon_seq[i:i+3]
            if codon not in CODON_TABLE:
                return TypeCheckResult(
                    predicate="InFrame",
                    verdict=Verdict.FAIL,
                    violation=f"Invalid codon '{codon}' at position {start + i}",
                )

    return TypeCheckResult(
        predicate="InFrame",
        verdict=Verdict.PASS,
    )


def evaluate_no_instability_motif(seq: str) -> TypeCheckResult:
    """Evaluate whether the sequence contains mRNA instability motifs.

    Detects two classes of instability motifs:
    1. ATTTA — the canonical AUUUA mRNA destabilizing element (5 bases)
    2. U-rich regions — 6 or more consecutive T's in DNA (corresponding to
       poly-U in mRNA), which are known mRNA degradation signals.

    Five consecutive T's (TTTTT) is NOT flagged, as the threshold is 6.

    Args:
        seq: DNA sequence to evaluate.

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    seq = seq.upper()
    positions = []

    # Check for ATTTA (AUUUA in mRNA)
    for i in range(len(seq) - 4):
        if seq[i:i+5] == "ATTTA":
            positions.append(i)

    # Check for U-rich regions (6+ consecutive T's in DNA)
    i = 0
    while i < len(seq):
        if seq[i] == "T":
            run_start = i
            while i < len(seq) and seq[i] == "T":
                i += 1
            run_len = i - run_start
            if run_len >= _INSTABILITY_T_RUN_MIN:
                positions.append(run_start)
        else:
            i += 1

    if positions:
        return TypeCheckResult(
            predicate="NoInstabilityMotif",
            verdict=Verdict.FAIL,
            violation=f"Instability motifs at positions {positions}",
        )
    return TypeCheckResult(
        predicate="NoInstabilityMotif",
        verdict=Verdict.PASS,
    )


def evaluate_no_unexpected_tm_domain(
    seq: str,
    is_cytosolic: bool = True,
    window_size: int = 19,
    threshold: float = 0.68,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate whether a cytosolic protein has gained unexpected TM domains.

    Checks for hydrophobic stretches in the translated protein that could
    form transmembrane domains — a common side-effect of mutagenesis on
    cytosolic proteins.

    Organism-aware: prokaryotes use window_size=17, eukaryotes use 19.
    Also applies flanking charge check (positive-inside rule) to distinguish
    true TM domains from soluble protein hydrophobic patches.

    Verdicts use the five-valued logic:
    - LIKELY_PASS: No hydrophobic stretches detected and protein is cytosolic
    - UNCERTAIN: Hydrophobic stretch without TM flanking charges, or
      borderline fraction
    - FAIL: Clear TM domains with appropriate flanking charges in a
      cytosolic protein
    - PASS: If not cytosolic (membrane protein), TM domains are expected

    Args:
        seq: DNA sequence to evaluate.
        is_cytosolic: Whether the protein is cytosolic (default True).
        window_size: Sliding window size in amino acids (default 19).
        threshold: Hydrophobic fraction threshold for FAIL (default 0.68).
        organism: Target organism name for organism-specific window sizing.

    Returns:
        TypeCheckResult with LIKELY_PASS/UNCERTAIN/FAIL/PASS verdict.
    """
    # Delegate to the low-level check for consistent organism/flanking logic
    result = check_no_unexpected_tm_domain(
        seq, is_cytosolic=is_cytosolic, window_size=window_size,
        threshold=threshold, organism=organism,
    )

    if not is_cytosolic:
        return TypeCheckResult(
            predicate="NoUnexpectedTMDomain",
            verdict=Verdict.PASS,
        )

    # Map PredicateResult verdict to TypeCheckResult verdict
    if result.verdict == Verdict.FAIL:
        return TypeCheckResult(
            predicate=f"NoUnexpectedTMDomain({is_cytosolic}, {threshold})",
            verdict=Verdict.FAIL,
            violation=result.details,
        )
    elif result.verdict == Verdict.UNCERTAIN:
        return TypeCheckResult(
            predicate=f"NoUnexpectedTMDomain({is_cytosolic}, {threshold})",
            verdict=Verdict.UNCERTAIN,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate=f"NoUnexpectedTMDomain({is_cytosolic}, {threshold})",
            verdict=Verdict.LIKELY_PASS,
        )


def evaluate_mrna_secondary_structure(
    seq: str,
    window_start: int = 0,
    window_end: int = 50,
    dg_threshold: float = -15.0,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate mRNA secondary structure stability around the RBS/start codon.

    Strong secondary structure in the 5' UTR or around the start codon can
    block ribosome binding and reduce translation efficiency. This predicate
    uses a simplified nearest-neighbor ΔG approximation.

    Organism-specific ΔG cutoffs:
    - Prokaryotes: ΔG < -15 kcal/mol is FAIL
    - Eukaryotes: ΔG < -25 kcal/mol is FAIL (less relevant for cap-dependent)
    - If organism not specified, uses dg_threshold

    Verdicts use the five-valued logic:
    - LIKELY_PASS: Weak structure (ΔG close to 0, easy for ribosome access)
    - UNCERTAIN: Moderate structure (ΔG between effective_threshold*0.7 and effective_threshold)
    - FAIL: Very stable structure that blocks ribosome binding

    Args:
        seq: DNA sequence to evaluate.
        window_start: Start position of the analysis window (default 0).
        window_end: End position of the analysis window (default 50).
        dg_threshold: ΔG threshold for FAIL (default -15.0 kcal/mol).
        organism: Target organism name for organism-specific cutoffs.

    Returns:
        TypeCheckResult with LIKELY_PASS/UNCERTAIN/FAIL verdict.
    """
    # Apply organism-specific ΔG threshold
    effective_threshold = dg_threshold
    if organism:
        if _is_prokaryotic_organism(organism):
            effective_threshold = _MRNA_DG_PROKARYOTE_FAIL
        else:
            effective_threshold = _MRNA_DG_EUKARYOTE_FAIL

    seq = seq.upper()
    effective_end = min(window_end, len(seq))
    window_seq = seq[window_start:effective_end]

    if len(window_seq) < 4:
        return TypeCheckResult(
            predicate="mRNASecondaryStructure",
            verdict=Verdict.LIKELY_PASS,
        )

    # Convert DNA to RNA for pairing analysis (T -> U)
    rna = window_seq.replace("T", "U")

    # Count potential base pairs using simplified hairpin model
    gc_pairs = 0
    au_pairs = 0
    gu_pairs = 0

    half = len(rna) // 2
    first_half = rna[:half]
    second_half = rna[half:2 * half]

    for i in range(min(len(first_half), len(second_half))):
        j = len(second_half) - 1 - i
        if j < 0:
            break
        base_5 = first_half[i]
        base_3 = second_half[j]

        if (base_5 == "G" and base_3 == "C") or (base_5 == "C" and base_3 == "G"):
            gc_pairs += 1
        elif (base_5 == "A" and base_3 == "U") or (base_5 == "U" and base_3 == "A"):
            au_pairs += 1
        elif (base_5 == "G" and base_3 == "U") or (base_5 == "U" and base_3 == "G"):
            gu_pairs += 1

    # Simplified nearest-neighbor ΔG estimate
    dg = _DG_GC_PAIR_KCAL * gc_pairs + _DG_AU_PAIR_KCAL * au_pairs + _DG_GU_PAIR_KCAL * gu_pairs

    if dg <= effective_threshold:
        verdict = Verdict.FAIL
    elif dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_PASS

    violation = None
    if verdict == Verdict.FAIL:
        violation = (
            f"Strong mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
            f"<= {effective_threshold}"
        )
    elif verdict == Verdict.UNCERTAIN:
        # Refine using GC content of the window
        window_gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) if window_seq else 0.5
        if window_gc < 0.5:
            verdict = Verdict.LIKELY_PASS
            violation = (
                f"Moderate mRNA secondary structure: ΔG={dg:.1f} kcal/mol, "
                f"but AT-rich window (GC={window_gc:.0%}) weakens structure"
            )
        else:
            verdict = Verdict.LIKELY_FAIL
            violation = (
                f"Moderate mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
                f"<= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f}"
            )

    return TypeCheckResult(
        predicate=f"mRNASecondaryStructure({window_start}, {window_end}, {effective_threshold})",
        verdict=verdict,
        violation=violation,
    )


def evaluate_no_cryptic_promoter(
    seq: str,
    organism: str = "E_coli",
    threshold: float = 0.7,
) -> TypeCheckResult:
    """Evaluate whether the sequence contains cryptic promoter sites.

    Uses position weight matrix scoring to scan for promoter motifs.
    For prokaryotes (E_coli), scans for -35 and -10 boxes; for eukaryotes,
    scans for TATA box and Initiator with multi-element requirement.

    Key improvement: requires MULTIPLE promoter elements within 50bp AND
    a TATA box to classify as cryptic promoter for eukaryotes. Single
    motifs alone should not trigger FAIL.

    Verdicts use the five-valued logic:
    - LIKELY_PASS: No significant promoter motifs (best score well below threshold)
    - UNCERTAIN: Borderline match (score between threshold*0.8 and threshold)
    - FAIL: Strong cryptic promoter match with multiple elements + TATA

    Args:
        seq: DNA sequence to evaluate.
        organism: Organism whose promoter consensus to use (default: "E_coli").
        threshold: Score threshold above which a match is FAIL.

    Returns:
        TypeCheckResult with LIKELY_PASS/UNCERTAIN/FAIL verdict.
    """
    result = check_no_cryptic_promoter(seq, organism, threshold)

    if result.verdict == Verdict.FAIL:
        return TypeCheckResult(
            predicate=f"NoCrypticPromoter({organism}, {threshold})",
            verdict=Verdict.FAIL,
            violation=result.details,
        )
    elif result.verdict == Verdict.UNCERTAIN:
        return TypeCheckResult(
            predicate=f"NoCrypticPromoter({organism}, {threshold})",
            verdict=Verdict.UNCERTAIN,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate=f"NoCrypticPromoter({organism}, {threshold})",
            verdict=Verdict.LIKELY_PASS,
        )


def evaluate_no_cpg_island(
    seq: str,
    window: int = 200,
    threshold: float = 0.6,
    organism: str = "",
) -> TypeCheckResult:
    """Evaluate whether the sequence contains CpG islands.

    A CpG island is detected when the observed/expected CG ratio exceeds
    the threshold in any sliding window of the specified size.

    CpG island avoidance is primarily relevant for mammalian expression
    systems where CpG methylation can lead to gene silencing.  For
    prokaryotic organisms (e.g. *E. coli*), CpG islands have no known
    regulatory significance, so the check is skipped when a prokaryotic
    organism is specified.

    GC-content-aware: for GC-rich targets (>60% GC), the CpG threshold
    is relaxed to 2x the normal threshold since CpG density is naturally
    higher in GC-rich sequences.

    Args:
        seq: DNA sequence to evaluate.
        window: Window size for CpG island scanning.
        threshold: Maximum allowed Obs/Exp CG ratio.
        organism: Target organism name.  If prokaryotic, the check is
            skipped and PASS is returned immediately.  If empty (default),
            the check runs as before for backward compatibility.

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    # Delegate to the low-level check for consistent organism/GC-aware logic
    result = check_no_cpg_island(seq, window=window, threshold=threshold, organism=organism)

    if result.verdict == Verdict.FAIL:
        return TypeCheckResult(
            predicate="NoCpGIsland",
            verdict=Verdict.FAIL,
            violation=result.details,
        )
    return TypeCheckResult(
        predicate="NoCpGIsland",
        verdict=Verdict.PASS,
    )


def analyze_codon_at_position(
    seq: str,
    position: int,
    organism: str = "Homo_sapiens",
) -> Dict[str, Any]:
    """Analyze the codon at a given position for optimality and alternatives.

    Args:
        seq: DNA coding sequence.
        position: 0-based nucleotide position (will be rounded to codon start).
        organism: Target organism for CAI lookup.

    Returns:
        Dict with keys: codon, amino_acid, cai, alternatives, position.
    """
    codon_start = (position // 3) * 3
    if codon_start + 3 > len(seq):
        return {"codon": "N/A", "amino_acid": "N/A", "cai": 0.0, "alternatives": [], "position": codon_start}

    seq = seq.upper()
    codon = seq[codon_start:codon_start + 3]
    aa = CODON_TABLE.get(codon, "?")
    species_cai = _resolve_species_cai(organism)
    cai = species_cai.get(codon, 0.0)

    alternatives = []
    for alt in AA_TO_CODONS.get(aa, []):
        if alt != codon:
            alternatives.append({
                "codon": alt,
                "cai": species_cai.get(alt, 0.0),
            })

    return {
        "codon": codon,
        "amino_acid": aa,
        "cai": cai,
        "alternatives": sorted(alternatives, key=lambda x: x["cai"], reverse=True),
        "position": codon_start,
    }


def evaluate_co_translational_folding(
    seq: str,
    organism: str = "Homo_sapiens",
    domain_boundaries: List[int] | None = None,
    min_pause_cai: float = 0.3,
    structure_confidence: float | None = None,
    plddt_score: float | None = None,
) -> TypeCheckResult:
    """Evaluate co-translational folding preservation after codon optimization.

    Checks whether codon optimization has disrupted critical pause sites
    that are important for proper protein folding during translation. The
    key insight is that slow codons (low tRNA adaptation) at domain
    boundaries allow the nascent chain to fold properly before downstream
    sequence emerges from the ribosome.

    When structure_confidence is provided, UNCERTAIN verdicts can be
    resolved: confidence > 0.7 → PASS, confidence < 0.5 → LIKELY_FAIL.

    Verdicts use the five-valued logic:
    - PASS: Good ramp with appropriate pause sites preserved.
    - LIKELY_PASS: Acceptable ramp with minor concerns.
    - UNCERTAIN: Ramp too fast or domain boundaries speeded-up.
    - LIKELY_FAIL: Multiple domain boundaries disrupted.
    - FAIL: Ramp destroyed AND domain boundaries disrupted.

    Args:
        seq: DNA coding sequence.
        organism: Target organism for CAI computation (default: "Homo_sapiens").
        domain_boundaries: Optional list of codon positions (0-based) where
            protein domains start/end. If provided, these positions are
            checked for speed-up (CAI > 0.7 where a pause is needed).
        min_pause_cai: CAI threshold below which a codon is considered a
            pause site (default 0.3).
        structure_confidence: Optional ESMFold confidence score (0–1).
        plddt_score: Optional pLDDT score from structure prediction.

    Returns:
        TypeCheckResult with PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL verdict.
    """
    seq = seq.upper()

    # Resolve organism name to flat codon→CAI-weight dict
    species_cai = _resolve_species_cai(organism)

    # Run the low-level predicate check
    result = check_co_translational_folding(
        seq, species_cai, domain_boundaries, min_pause_cai,
        structure_confidence=structure_confidence,
        plddt_score=plddt_score,
    )

    # Map the PredicateResult verdict to a TypeCheckResult verdict using
    # the five-valued logic, incorporating additional context about the
    # severity of the findings.
    ramp_info = _compute_codon_ramp_score(seq, species_cai)
    ramp_score = ramp_info["ramp_score"]
    speed_disruptions = ramp_info["speed_disruptions"]
    pause_sites = ramp_info["pause_sites"]

    ramp_length = min(_CODON_RAMP_LENGTH, len(seq) // 3)
    ramp_all_fast = ramp_length >= _MIN_RAMP_FOR_WARNING and all(
        species_cai.get(seq[i * 3:(i + 1) * 3], 0.0) > _FAST_CODON_CAI_THRESHOLD
        for i in range(ramp_length)
    )

    # Count domain boundary disruptions
    domain_disrupted = 0
    if domain_boundaries:
        num_codons = len(seq) // 3
        for bp in domain_boundaries:
            if 0 <= bp < num_codons:
                codon = seq[bp * 3:(bp + 1) * 3]
                if species_cai.get(codon, 0.0) > _FAST_CODON_CAI_THRESHOLD:
                    domain_disrupted += 1

    # Determine the TypeCheckResult verdict
    if domain_boundaries and domain_disrupted > 0 and ramp_all_fast:
        # Ramp destroyed AND domain boundaries disrupted
        verdict = Verdict.FAIL
        violation = (
            f"Ramp destroyed (avg CAI={ramp_score:.3f}) and {domain_disrupted} "
            f"domain boundary(ies) disrupted by fast codons"
        )
    elif domain_boundaries and domain_disrupted >= 2:
        # Multiple domain boundaries disrupted
        verdict = Verdict.LIKELY_FAIL
        violation = (
            f"{domain_disrupted} domain boundary(ies) have fast codons "
            f"(CAI > 0.7) where pause sites are needed"
        )
    elif domain_boundaries and domain_disrupted == 1:
        # Single domain boundary disrupted
        verdict = Verdict.UNCERTAIN
        violation = (
            f"1 domain boundary has a fast codon (CAI > 0.7) where "
            f"a pause site may be needed"
        )
    elif ramp_all_fast and ramp_length >= _MIN_RAMP_FOR_WARNING:
        # Ramp too fast — ribosome jam risk
        verdict = Verdict.UNCERTAIN
        violation = (
            f"Codon ramp (first {ramp_length} codons) is entirely fast "
            f"(avg CAI={ramp_score:.3f}) — ribosome jam risk"
        )
    elif speed_disruptions:
        # Some pause sites may have been replaced by fast codons
        if len(speed_disruptions) <= 2:
            verdict = Verdict.LIKELY_PASS
            violation = (
                f"{len(speed_disruptions)} potential pause site(s) replaced by "
                f"fast codons — minor concern"
            )
        else:
            verdict = Verdict.UNCERTAIN
            violation = (
                f"{len(speed_disruptions)} potential pause site(s) replaced by "
                f"fast codons — may affect co-translational folding"
            )
    elif pause_sites:
        # Good — pause sites preserved
        verdict = Verdict.PASS
        violation = None
    else:
        # No pause sites, but no obvious issues either
        verdict = Verdict.LIKELY_PASS
        violation = None

    return TypeCheckResult(
        predicate=f"CoTranslationalFolding({organism}, {min_pause_cai})",
        verdict=verdict,
        violation=violation,
    )


def evaluate_no_stop_codons(seq: str) -> TypeCheckResult:
    """Evaluate whether the DNA sequence contains internal stop codons.

    The last codon in the reading frame is allowed to be a stop
    (natural termination). Only stops before the last codon are
    flagged as violations.

    Args:
        seq: DNA coding sequence to evaluate.

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    result = check_no_stop_codons(seq)
    if result.passed:
        return TypeCheckResult(
            predicate="NoStopCodons",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="NoStopCodons",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_no_gt_dinucleotide(
    seq: str,
    organism: str = "",
    max_gt_count: int | None = None,
) -> TypeCheckResult:
    """Evaluate whether the sequence contains avoidable GT dinucleotides.

    Uses the soft GT check (check_no_gt_dinucleotide_soft) which is
    organism-aware: for eukaryotes, GT avoidance is a soft constraint
    (LIKELY_FAIL when GT count exceeds max_gt_count); for prokaryotes,
    it's a hard constraint (FAIL for any GT).

    Args:
        seq: DNA coding sequence to evaluate.
        organism: Target organism.  If prokaryotic, auto-PASS (skipped).
        max_gt_count: Maximum GT count before soft-fail. If None,
            auto-computed from sequence length (1 per 50bp for eukaryotes).

    Returns:
        TypeCheckResult with verdict matching the soft check logic.
    """
    result = check_no_gt_dinucleotide_soft(
        seq, organism=organism, max_gt_count=max_gt_count,
    )

    if result.verdict == Verdict.PASS:
        return TypeCheckResult(
            predicate="NoGTDinucleotide",
            verdict=Verdict.PASS,
        )
    elif result.verdict == Verdict.LIKELY_FAIL:
        return TypeCheckResult(
            predicate="NoGTDinucleotide",
            verdict=Verdict.LIKELY_FAIL,
            violation=result.details,
        )
    else:
        return TypeCheckResult(
            predicate="NoGTDinucleotide",
            verdict=Verdict.FAIL,
            violation=result.details,
        )


def evaluate_valid_coding_seq(seq: str) -> TypeCheckResult:
    """Evaluate whether the DNA sequence is a valid coding sequence.

    Checks that the length is divisible by 3 and all codons are
    in the standard genetic code.

    Args:
        seq: DNA coding sequence to evaluate.

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    result = check_valid_coding_seq(seq)
    if result.passed:
        return TypeCheckResult(
            predicate="ValidCodingSeq",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="ValidCodingSeq",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_conservation_score(
    seq: str,
    protein: str = "",
    min_score: int = 0,
) -> TypeCheckResult:
    """Evaluate whether the BLOSUM62 conservation score meets the threshold.

    Translates the DNA sequence and compares each position against the
    target protein using the BLOSUM62 substitution matrix.  After correct
    optimization the two sequences are identical, so every diagonal score
    is positive and the predicate should PASS.

    Args:
        seq: DNA coding sequence to evaluate (translated internally).
        protein: Target protein sequence (amino-acid string).  If empty,
            the predicate auto-PASSes (no comparison target).
        min_score: Minimum BLOSUM62 score per position for PASS (default 0).

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    if not protein:
        return TypeCheckResult(
            predicate="ConservationScore",
            verdict=Verdict.PASS,
        )
    result = check_conservation_score(seq, protein, min_score=min_score)
    if result.passed:
        return TypeCheckResult(
            predicate="ConservationScore",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="ConservationScore",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_codon_optimality(
    seq: str,
    organism: str = "Homo_sapiens",
    threshold: float = 0.5,
) -> TypeCheckResult:
    """Evaluate whether the codon adaptation index (CAI) meets the threshold.

    Computes the geometric-mean CAI across all codons in the sequence
    using the organism-specific codon adaptiveness table.

    Args:
        seq: DNA coding sequence to evaluate.
        organism: Target organism for CAI computation.
        threshold: Minimum CAI score for PASS (default 0.5).

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    result = check_codon_optimality(seq, organism, min_cai=threshold)
    if result.passed:
        return TypeCheckResult(
            predicate="CodonOptimality",
            verdict=Verdict.PASS,
        )
    return TypeCheckResult(
        predicate="CodonOptimality",
        verdict=Verdict.FAIL,
        violation=result.details,
    )


def evaluate_all_predicates(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    known_exon_boundaries: List[Tuple[int, int]] | None = None,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: List[str] | None = None,
    cryptic_threshold: float = 3.0,
    uncertain_lo: float = 0.0,
    cai_threshold: float = 0.5,
    promoter_threshold: float = 0.7,
    is_cytosolic: bool = True,
    tm_threshold: float = 0.68,
    mrna_window: int = 50,
    mrna_dg_threshold: float = -15.0,
    folding_threshold: float = 0.3,
    domain_boundaries: List[int] | None = None,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> List[TypeCheckResult]:
    """Evaluate all 12 type predicates against a sequence.

    The 12 predicates are:
    1. NoCrypticSplice — no cryptic splice donors
    2. SpliceCorrect — canonical splice signals at intron boundaries
    3. GCInRange — GC content within acceptable range
    4. CodonAdapted — CAI above threshold
    5. NoRestrictionSite — no restriction enzyme sites
    6. InFrame — valid coding frame
    7. NoInstabilityMotif — no mRNA instability motifs
    8. NoCpGIsland — no CpG islands
    9. NoCrypticPromoter — no cryptic promoter sites
    10. NoUnexpectedTMDomain — no unexpected transmembrane domains
    11. mRNASecondaryStructure — no strong mRNA secondary structure
    12. CoTranslationalFolding — co-translational folding pause-site preservation

    SLOT Mode:
    Predicates 1, 4, 9, 10, 11, 12 are SLOT-dependent (rely on heuristic
    scanners or external tools). Their behavior depends on slot_mode:

    - CONSERVATIVE (default): SLOT predicates return UNCERTAIN, matching
      the Lean4 formal model exactly.
    - VERIFIED: SLOT predicates return PASS when verification conditions
      are met (tool available + result meets threshold).
    - PERMISSIVE: SLOT predicates return PASS with weaker evidence
      thresholds.

    Args:
        seq: DNA sequence to evaluate.
        boundaries: Exon boundary tuples [(start, end), ...].
        organism: Target organism for CAI and promoter computation.
        gc_lo: Minimum GC fraction.
        gc_hi: Maximum GC fraction.
        enzymes: Restriction enzymes to check. If None, uses a default set
            of common cloning enzymes: EcoRI, BamHI, XhoI, HindIII, NotI.
        cryptic_threshold: MaxEnt score threshold for cryptic splice FAIL.
        uncertain_lo: MaxEnt score threshold for UNCERTAIN.
        cai_threshold: Minimum CAI for CodonAdapted PASS.
        promoter_threshold: Score threshold for cryptic promoter FAIL.
        is_cytosolic: Whether the protein is cytosolic (for TM domain check).
        tm_threshold: Hydrophobic fraction threshold for TM domain FAIL.
        mrna_window: Window size for mRNA secondary structure analysis.
        mrna_dg_threshold: ΔG threshold for mRNA structure FAIL.
        folding_threshold: Minimum CAI for a codon to be considered a pause
            site (used by CoTranslationalFolding predicate).
        domain_boundaries: Codon positions where protein domains start/end
            (used by CoTranslationalFolding predicate).
        slot_mode: SLOT evaluation mode (default CONSERVATIVE).

    Returns:
        List of 12 TypeCheckResult objects.
    """
    # Backward compatibility: accept known_exon_boundaries as alias for boundaries
    if boundaries is None and known_exon_boundaries is not None:
        boundaries = known_exon_boundaries

    # Track whether the caller supplied real exon boundaries.
    # If neither boundaries nor known_exon_boundaries was given, the
    # sequence is treated as a single-exon gene and splice-correct
    # checking is skipped (no intron boundaries → nothing to check).
    has_real_boundaries = boundaries is not None and len(boundaries) >= 2

    if boundaries is None:
        boundaries = [(0, len(seq))]

    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
    # Map organism name to promoter organism key
    promoter_organism = "E_coli" if organism in ("E_coli", "ecoli", "Escherichia_coli") else "eukaryote"

    # Import slot_verification lazily to avoid circular imports
    from .slot_verification import is_slot_predicate, verify_slot_predicate

    # Core (non-SLOT) predicates: always evaluate normally
    # SpliceCorrect is only evaluated when real exon boundaries were provided;
    # otherwise we emit an auto-PASS (no introns → no splice correction needed).
    if has_real_boundaries:
        splice_result = evaluate_splice_correct(seq, boundaries)
    else:
        splice_result = TypeCheckResult(
            predicate="SpliceCorrect",
            verdict=Verdict.PASS,
        )

    results: List[TypeCheckResult] = [
        splice_result,
        evaluate_gc_in_range(seq, gc_lo, gc_hi),
        evaluate_no_restriction_site(seq, enzymes),
        evaluate_in_frame(seq, boundaries),
        evaluate_no_instability_motif(seq),
        evaluate_no_cpg_island(seq, organism=organism),
    ]

    # SLOT predicates: behavior depends on slot_mode
    slot_predicates = [
        ("NoCrypticSplice", lambda: evaluate_no_cryptic_splice(seq, boundaries, cryptic_threshold, uncertain_lo, organism=organism)),
        ("CodonAdapted", lambda: evaluate_codon_adapted(seq, organism, cai_threshold)),
        ("NoCrypticPromoter", lambda: evaluate_no_cryptic_promoter(seq, promoter_organism, promoter_threshold)),
        ("NoUnexpectedTMDomain", lambda: evaluate_no_unexpected_tm_domain(seq, is_cytosolic, 19, tm_threshold)),
        ("mRNASecondaryStructure", lambda: evaluate_mrna_secondary_structure(seq, 0, mrna_window, mrna_dg_threshold)),
        ("CoTranslationalFolding", lambda: evaluate_co_translational_folding(seq, organism, domain_boundaries, folding_threshold)),
    ]

    for pred_name, eval_fn in slot_predicates:
        if slot_mode == SLOTMode.CONSERVATIVE and is_slot_predicate(pred_name):
            # CONSERVATIVE: always UNCERTAIN for SLOT predicates
            result = eval_fn()
            result.verdict = Verdict.UNCERTAIN
            result.knowledge_gap = f"SLOT predicate: {pred_name} returns UNCERTAIN in CONSERVATIVE mode"
            results.append(result)
        elif slot_mode in (SLOTMode.VERIFIED, SLOTMode.PERMISSIVE) and is_slot_predicate(pred_name):
            # VERIFIED/PERMISSIVE: use slot verification conditions
            verdict, evidence = verify_slot_predicate(
                pred_name,
                slot_mode=slot_mode,
                seq=seq,
                low_thresh=uncertain_lo,
                high_thresh=cryptic_threshold,
                organism=promoter_organism,
                threshold=promoter_threshold,
                is_cytosolic=is_cytosolic,
                window_end=mrna_window,
                dg_threshold=mrna_dg_threshold,
                domain_boundaries=domain_boundaries,
                min_pause_cai=folding_threshold,
            )
            result = eval_fn()
            result.verdict = verdict
            result.knowledge_gap = f"SLOT predicate: {pred_name} evaluated in {slot_mode.value} mode"
            if evidence.verified:
                result.derivation = [{"evidence": evidence.to_dict()}]
            results.append(result)
        else:
            # Non-SLOT predicate: evaluate normally
            results.append(eval_fn())

    # Reorder to match the expected predicate order:
    # 1.NoCrypticSplice, 2.SpliceCorrect, 3.GCInRange, 4.CodonAdapted,
    # 5.NoRestrictionSite, 6.InFrame, 7.NoInstabilityMotif, 8.NoCpGIsland,
    # 9.NoCrypticPromoter, 10.NoUnexpectedTMDomain, 11.mRNASecondaryStructure,
    # 12.CoTranslationalFolding
    name_to_result = {r.predicate: r for r in results}
    # Handle parameterized names (e.g., CodonAdapted(Homo_sapiens, 0.5))
    ordered = []
    for canonical in [
        "NoCrypticSplice", "SpliceCorrect", "GCInRange", "CodonAdapted",
        "NoRestrictionSite", "InFrame", "NoInstabilityMotif", "NoCpGIsland",
        "NoCrypticPromoter", "NoUnexpectedTMDomain", "mRNASecondaryStructure",
        "CoTranslationalFolding",
    ]:
        if canonical in name_to_result:
            ordered.append(name_to_result[canonical])
        else:
            # Try prefix match for parameterized names
            for name, result in name_to_result.items():
                if name.startswith(canonical):
                    ordered.append(result)
                    break

    # If ordering didn't work, fall back to original
    if len(ordered) != 12:
        ordered = results

    return ordered


# ════════════════════════════════════════════════════════════════
# Predicate Registry — named dispatch for certificate verification
# ════════════════════════════════════════════════════════════════

from .exceptions import UnknownPredicateError


class PredicateRegistry:
    """Registry of named type predicates with evaluate() and verify() dispatch.

    The registry provides a single entry point for certificate generation
    and verification, mapping predicate names to their evaluate functions.
    It supports both evaluation (with default parameters) and verification
    (re-running a predicate with specific parameters from a certificate).
    """

    def __init__(self) -> None:
        self._predicates: Dict[str, callable] = {}
        self._verify_params: Dict[str, Dict[str, str]] = {}

    def register(
        self,
        name: str,
        fn: callable,
        verify_param_map: Dict[str, str] | None = None,
    ) -> None:
        """Register a predicate evaluation function.

        Args:
            name: Predicate name (e.g., 'NoCrypticSplice').
            fn: Callable that returns TypeCheckResult.
            verify_param_map: Optional mapping from certificate param names
                to function kwarg names.
        """
        self._predicates[name] = fn
        if verify_param_map:
            self._verify_params[name] = verify_param_map

    def names(self) -> List[str]:
        """Return sorted list of registered predicate names."""
        return sorted(self._predicates.keys())

    def evaluate(self, name: str, **kwargs) -> TypeCheckResult:
        """Evaluate a named predicate.

        Args:
            name: Predicate name.
            **kwargs: Arguments to pass to the predicate function.

        Returns:
            TypeCheckResult from the predicate.

        Raises:
            UnknownPredicateError: if name is not registered.
        """
        if name not in self._predicates:
            raise UnknownPredicateError(name)
        return self._predicates[name](**kwargs)

    def verify(self, name: str, **kwargs) -> TypeCheckResult:
        """Verify a predicate — same as evaluate but with certificate params.

        The verify method maps certificate parameter names to the function's
        expected kwargs before calling evaluate.

        Args:
            name: Predicate name.
            **kwargs: Certificate parameters, possibly with different names
                than the evaluate function expects.

        Returns:
            TypeCheckResult from re-evaluation.

        Raises:
            UnknownPredicateError: if name is not registered.
        """
        if name not in self._predicates:
            raise UnknownPredicateError(name)

        # Map certificate-style kwargs to evaluate-style kwargs
        param_map = self._verify_params.get(name, {})
        mapped_kwargs = {}
        for cert_key, fn_key in param_map.items():
            if cert_key in kwargs:
                mapped_kwargs[fn_key] = kwargs[cert_key]

        # Pass through any kwargs that match the function's signature directly
        fn = self._predicates[name]
        import inspect
        sig = inspect.signature(fn)
        for key, val in kwargs.items():
            if key in sig.parameters:
                mapped_kwargs[key] = val

        return fn(**mapped_kwargs)

    def __contains__(self, name: str) -> bool:
        return name in self._predicates


# ────────────────────────────────────────────────────────────
# Global registry instance with all 12 predicates registered
# ────────────────────────────────────────────────────────────
registry = PredicateRegistry()

registry.register(
    "NoCrypticSplice",
    evaluate_no_cryptic_splice,
    verify_param_map={
        "seq": "seq",
        "known_exon_boundaries": "boundaries",
        "organism": "organism",
        "cryptic_splice_threshold": "cryptic_threshold",
    },
)

registry.register(
    "SpliceCorrect",
    evaluate_splice_correct,
    verify_param_map={
        "seq": "seq",
        "known_exon_boundaries": "boundaries",
        "cellular_context": "cellular_context",
    },
)

registry.register(
    "GCInRange",
    evaluate_gc_in_range,
    verify_param_map={
        "seq": "seq",
        "gc_lo": "gc_lo",
        "gc_hi": "gc_hi",
    },
)

# SlidingGC predicate (local/sliding-window GC constraint)
from .sliding_gc import evaluate_sliding_gc as _evaluate_sliding_gc

registry.register(
    "SlidingGC",
    _evaluate_sliding_gc,
    verify_param_map={
        "seq": "seq",
        "window_size": "window_size",
        "gc_min": "gc_min",
        "gc_max": "gc_max",
    },
)

registry.register(
    "CodonAdapted",
    evaluate_codon_adapted,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "threshold": "threshold",
    },
)

registry.register(
    "NoRestrictionSite",
    evaluate_no_restriction_site,
    verify_param_map={
        "seq": "seq",
        "enzyme_set": "enzymes",
    },
)

registry.register(
    "InFrame",
    evaluate_in_frame,
    verify_param_map={
        "seq": "seq",
        "exon_boundaries": "boundaries",
    },
)

registry.register(
    "NoInstabilityMotif",
    evaluate_no_instability_motif,
    verify_param_map={
        "seq": "seq",
    },
)

registry.register(
    "NoCpGIsland",
    evaluate_no_cpg_island,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
    },
)

registry.register(
    "NoStopCodons",
    evaluate_no_stop_codons,
    verify_param_map={
        "seq": "seq",
    },
)

registry.register(
    "NoGTDinucleotide",
    evaluate_no_gt_dinucleotide,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
    },
)

registry.register(
    "ValidCodingSeq",
    evaluate_valid_coding_seq,
    verify_param_map={
        "seq": "seq",
    },
)

registry.register(
    "ConservationScore",
    evaluate_conservation_score,
    verify_param_map={
        "seq": "seq",
        "protein": "protein",
        "min_score": "min_score",
    },
)

registry.register(
    "CodonOptimality",
    evaluate_codon_optimality,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "threshold": "threshold",
    },
)

registry.register(
    "NoCrypticPromoter",
    evaluate_no_cryptic_promoter,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "threshold": "threshold",
    },
)

registry.register(
    "CoTranslationalFolding",
    evaluate_co_translational_folding,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "domain_boundaries": "domain_boundaries",
        "min_pause_cai": "min_pause_cai",
    },
)

registry.register(
    "NoUnexpectedTMDomain",
    evaluate_no_unexpected_tm_domain,
    verify_param_map={
        "seq": "seq",
        "is_cytosolic": "is_cytosolic",
        "threshold": "threshold",
    },
)

registry.register(
    "mRNASecondaryStructure",
    evaluate_mrna_secondary_structure,
    verify_param_map={
        "seq": "seq",
        "window_start": "window_start",
        "window_end": "window_end",
        "dg_threshold": "dg_threshold",
    },
)

# ────────────────────────────────────────────────────────────
# Protein-level predicates (stability, solubility,
# immunogenicity, structure)
# ────────────────────────────────────────────────────────────
from .stability_predicates import (
    evaluate_stable_folding,
    evaluate_no_destabilizing_mutation,
    evaluate_disulfide_bond_integrity,
    evaluate_hydrophobic_core_quality,
)
from .solubility_predicates import (
    evaluate_soluble_expression,
    evaluate_no_aggregation_prone_region,
    evaluate_charge_composition,
    evaluate_no_long_hydrophobic_stretch,
)
from .immuno_predicates import (
    evaluate_low_immunogenicity,
    evaluate_no_strong_t_cell_epitope,
    evaluate_no_dominant_b_cell_epitope,
    evaluate_population_coverage_safe,
)
from .structure.predicates import (
    evaluate_structure_confidence,
    evaluate_no_misfolding_risk,
    evaluate_correct_fold_topology,
    evaluate_no_unexpected_interaction,
)

# Structure predicates
registry.register(
    "StructureConfidence",
    evaluate_structure_confidence,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "NoMisfoldingRisk",
    evaluate_no_misfolding_risk,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "CorrectFoldTopology",
    evaluate_correct_fold_topology,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "NoUnexpectedInteraction",
    evaluate_no_unexpected_interaction,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)

# Stability predicates
registry.register(
    "StableFolding",
    evaluate_stable_folding,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "stability_threshold": "stability_threshold",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "NoDestabilizingMutation",
    evaluate_no_destabilizing_mutation,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "max_ddg": "max_ddg",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "DisulfideBondIntegrity",
    evaluate_disulfide_bond_integrity,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "HydrophobicCoreQuality",
    evaluate_hydrophobic_core_quality,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)

# Solubility predicates
registry.register(
    "SolubleExpression",
    evaluate_soluble_expression,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "min_solubility_score": "min_solubility_score",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "NoAggregationProneRegion",
    evaluate_no_aggregation_prone_region,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "ChargeComposition",
    evaluate_charge_composition,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)
registry.register(
    "NoLongHydrophobicStretch",
    evaluate_no_long_hydrophobic_stretch,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "pdb_string": "pdb_string",
    },
)

# Immunogenicity predicates
registry.register(
    "LowImmunogenicity",
    evaluate_low_immunogenicity,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "max_immunogenicity_score": "max_immunogenicity_score",
    },
)
registry.register(
    "NoStrongTCellEpitope",
    evaluate_no_strong_t_cell_epitope,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "mhc_alleles": "mhc_alleles",
    },
)
registry.register(
    "NoDominantBCellEpitope",
    evaluate_no_dominant_b_cell_epitope,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
    },
)
registry.register(
    "PopulationCoverageSafe",
    evaluate_population_coverage_safe,
    verify_param_map={
        "sequence": "sequence",
        "protein": "protein",
        "organism": "organism",
        "mhc_alleles": "mhc_alleles",
    },
)