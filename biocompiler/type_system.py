"""
BioCompiler Type System v7.0.0
==============================
Defines the core types, codon tables, BLOSUM62 matrix, and 8 predicate classes
for certified gene optimization.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List, Dict, Set, Optional, Tuple

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
# ────────────────────────────────────────────────────────────
_BLOSUM62_ROWS = [
    #  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],  # A
    [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
    [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],  # N
    [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],  # D
    [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
    [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],  # Q
    [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],  # E
    [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],  # G
    [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],  # H
    [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],  # I
    [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],  # L
    [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1, -1, -1, -3, -2, -2],  # K
    [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],  # M
    [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],  # F
    [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],  # P
    [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],  # S
    [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],  # T
    [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],  # W
    [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],  # Y
    [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],  # V
]

_BLOSUM_INDEX = list("ARNDCQEGHILKMFPSTWYV")

BLOSUM62: Dict[Tuple[str, str], int] = {}
for _i, _a1 in enumerate(_BLOSUM_INDEX):
    for _j, _a2 in enumerate(_BLOSUM_INDEX):
        BLOSUM62[(_a1, _a2)] = _BLOSUM62_ROWS[_i][_j]


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
# 8 Predicate Classes for Certified Optimization
# ────────────────────────────────────────────────────────────
PREDICATE_NAMES = [
    "NoStopCodons",        # 1 — no internal stops
    "NoCrypticSplice",     # 2 — dual-threshold splice check
    "NoCpGIsland",         # 3 — CpG island avoidance
    "NoRestrictionSite",   # 4 — enzyme site removal
    "NoGTDinucleotide",    # 5 — GT dinucleotide avoidance (cross-codon aware)
    "ValidCodingSeq",      # 6 — in-frame, valid codons only
    "ConservationScore",   # 7 — BLOSUM62-based AA conservation
    "CodonOptimality",     # 8 — CAI-based codon quality
]


@dataclass
class PredicateResult:
    """Result of checking one predicate against a sequence."""
    predicate: str
    passed: bool
    verdict: Optional[SpliceVerdict] = None  # used by NoCrypticSplice
    details: str = ""
    positions: List[int] = field(default_factory=list)


# ────────────────────────────────────────────────────────────
# 8 Predicate check functions
# ────────────────────────────────────────────────────────────

def check_no_stop_codons(seq: str) -> PredicateResult:
    """Predicate 1: No internal stop codons."""
    violations = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in ("TAA", "TAG", "TGA"):
            violations.append(i)
    if violations:
        return PredicateResult("NoStopCodons", False, details="Internal stop codons found", positions=violations)
    return PredicateResult("NoStopCodons", True, details="No internal stop codons")


def check_no_cryptic_splice(seq: str, low_thresh: float = 3.0, high_thresh: float = 6.0) -> PredicateResult:
    """Predicate 2: No cryptic splice sites (dual-threshold PASS/UNCERTAIN/FAIL)."""
    from .splice import maxent_score
    gt_positions = []
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "GT":
            gt_positions.append(i)
    if not gt_positions:
        return PredicateResult("NoCrypticSplice", True, verdict=SpliceVerdict.PASS,
                               details="No GT dinucleotides found")
    max_score = 0.0
    worst_pos = -1
    worst_verdict = SpliceVerdict.PASS
    for pos in gt_positions:
        context_start = max(0, pos - 3)
        context_end = min(len(seq), pos + 6)
        context = seq[context_start:context_end]
        score = maxent_score(context)
        if score < low_thresh:
            v = SpliceVerdict.PASS
        elif score < high_thresh:
            v = SpliceVerdict.UNCERTAIN
        else:
            v = SpliceVerdict.FAIL
        if score > max_score:
            max_score = score
            worst_pos = pos
            worst_verdict = v

    passed = worst_verdict != SpliceVerdict.FAIL
    return PredicateResult("NoCrypticSplice", passed, verdict=worst_verdict,
                           details=f"Worst splice score {max_score:.2f} at pos {worst_pos}",
                           positions=[worst_pos] if worst_pos >= 0 else [])


def check_no_cpg_island(seq: str, window: int = 200, threshold: float = 0.6) -> PredicateResult:
    """Predicate 3: No CpG islands (Obs/Exp CG ratio > threshold in any window)."""
    worst_ratio = 0.0
    worst_start = -1
    for start in range(0, len(seq) - window + 1):
        window_seq = seq[start:start + window]
        c_count = window_seq.count("C")
        g_count = window_seq.count("G")
        cg_count = sum(1 for i in range(len(window_seq) - 1) if window_seq[i:i+2] == "CG")
        expected = (c_count * g_count) / len(window_seq) if len(window_seq) > 0 else 0
        obs_exp = cg_count / expected if expected > 0 else 0.0
        if obs_exp > worst_ratio:
            worst_ratio = obs_exp
            worst_start = start
    if worst_ratio > threshold:
        return PredicateResult("NoCpGIsland", False,
                               details=f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f} > {threshold}",
                               positions=[worst_start])
    return PredicateResult("NoCpGIsland", True,
                           details=f"Worst CpG Obs/Exp ratio {worst_ratio:.3f} <= {threshold}")


def check_no_restriction_site(seq: str, enzymes: List[str]) -> PredicateResult:
    """Predicate 4: No restriction enzyme recognition sites."""
    from .restriction_sites import get_recognition_site
    violations = []
    for enzyme in enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        pos = seq.find(site)
        while pos != -1:
            violations.append(pos)
            pos = seq.find(site, pos + 1)
    if violations:
        return PredicateResult("NoRestrictionSite", False,
                               details=f"Restriction sites found at {violations}",
                               positions=violations)
    return PredicateResult("NoRestrictionSite", True, details="No restriction sites found")


def check_no_gt_dinucleotide(seq: str) -> PredicateResult:
    """Predicate 5: No GT dinucleotides (5' splice donor mimic), including cross-codon."""
    positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
    if positions:
        return PredicateResult("NoGTDinucleotide", False,
                               details=f"GT dinucleotides at {positions}",
                               positions=positions)
    return PredicateResult("NoGTDinucleotide", True, details="No GT dinucleotides found")


def check_valid_coding_seq(seq: str) -> PredicateResult:
    """Predicate 6: Valid coding sequence (length divisible by 3, all valid codons)."""
    if len(seq) % 3 != 0:
        return PredicateResult("ValidCodingSeq", False,
                               details=f"Sequence length {len(seq)} not divisible by 3")
    invalid = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon not in CODON_TABLE:
            invalid.append((i, codon))
    if invalid:
        return PredicateResult("ValidCodingSeq", False,
                               details=f"Invalid codons: {invalid}")
    return PredicateResult("ValidCodingSeq", True, details="All codons valid")


def check_conservation_score(original_aa: str, new_aa: str, min_score: int = 0) -> PredicateResult:
    """Predicate 7: BLOSUM62 conservation score for amino acid substitution."""
    score = BLOSUM62.get((original_aa, new_aa), -10)
    passed = score >= min_score
    return PredicateResult("ConservationScore", passed,
                           details=f"BLOSUM62({original_aa},{new_aa})={score}, min={min_score}")


def check_codon_optimality(codon: str, species_cai: Dict[str, float], min_cai: float = 0.0) -> PredicateResult:
    """Predicate 8: Codon optimality (CAI score above threshold)."""
    cai = species_cai.get(codon, 0.0)
    passed = cai >= min_cai
    return PredicateResult("CodonOptimality", passed,
                           details=f"CAI({codon})={cai:.4f}, min={min_cai}")


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
    """Find restriction sites that span codon boundaries."""
    positions = []
    site_len = len(site)
    for i in range(len(seq) - site_len + 1):
        if seq[i:i+site_len] == site:
            codon_start_i = (i // 3) * 3
            codon_end_i = ((i + site_len - 1) // 3) * 3 + 3
            if codon_end_i - codon_start_i > 3:
                positions.append(i)
    return positions
