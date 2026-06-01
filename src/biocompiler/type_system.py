"""
BioCompiler Type System v7.0.0
==============================
Defines the core types, codon tables, BLOSUM62 matrix, and 8 predicate classes
for certified gene optimization.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List, Dict, Set, Optional, Tuple
from .types import Verdict

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
    verdict: Optional[Verdict] = None  # used by NoCrypticSplice and others
    details: str = ""
    positions: List[int] = field(default_factory=list)


# ────────────────────────────────────────────────────────────
# 8 Predicate check functions
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


def check_no_cryptic_splice(seq: str, low_thresh: float = 3.0, high_thresh: float = 6.0) -> PredicateResult:
    """Predicate 2: No cryptic splice sites (dual-threshold PASS/UNCERTAIN/FAIL)."""
    from .splice import maxent_score
    gt_positions = []
    for i in range(len(seq) - 1):
        if seq[i:i+2] == "GT":
            gt_positions.append(i)
    if not gt_positions:
        return PredicateResult("NoCrypticSplice", True, verdict=Verdict.PASS,
                               details="No GT dinucleotides found")
    max_score = 0.0
    worst_pos = -1
    worst_verdict = Verdict.PASS
    for pos in gt_positions:
        context_start = max(0, pos - 3)
        context_end = min(len(seq), pos + 6)
        context = seq[context_start:context_end]
        score = maxent_score(context)
        if score < low_thresh:
            v = Verdict.PASS
        elif score < high_thresh:
            v = Verdict.UNCERTAIN
        else:
            v = Verdict.FAIL
        if score > max_score:
            max_score = score
            worst_pos = pos
            worst_verdict = v

    passed = worst_verdict != Verdict.FAIL
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
        return PredicateResult("NoCpGIsland", False, verdict=Verdict.FAIL,
                               details=f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f} > {threshold}",
                               positions=[worst_start])
    return PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
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
        return PredicateResult("NoRestrictionSite", False, verdict=Verdict.FAIL,
                               details=f"Restriction sites found at {violations}",
                               positions=violations)
    return PredicateResult("NoRestrictionSite", True, verdict=Verdict.PASS, details="No restriction sites found")


def check_no_gt_dinucleotide(seq: str) -> PredicateResult:
    """Predicate 5: No GT dinucleotides (5' splice donor mimic), including cross-codon.

    This is the STRICT version — any GT fails the predicate.
    """
    positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
    if positions:
        return PredicateResult("NoGTDinucleotide", False, verdict=Verdict.FAIL,
                               details=f"GT dinucleotides at {positions}",
                               positions=positions)
    return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS, details="No GT dinucleotides found")


def check_no_avoidable_gt(seq: str) -> PredicateResult:
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
    """
    gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]
    if not gt_positions:
        return PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS, details="No GT dinucleotides found")

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
                # Can only change the current codon
                for c_alt in AA_TO_CODONS.get(curr_aa, [curr_codon]):
                    if prev_codon[-1] + c_alt[0] != "GT" and "GT" not in c_alt:
                        has_avoidable = True
                        break
            elif curr_aa == "*":
                # Can only change the previous codon
                for p_alt in AA_TO_CODONS.get(prev_aa, [prev_codon]):
                    if p_alt[-1] + curr_codon[0] != "GT" and "GT" not in p_alt:
                        has_avoidable = True
                        break
            else:
                # Both are regular AAs — try all combinations
                prev_alts = AA_TO_CODONS.get(prev_aa, [prev_codon])
                curr_alts = AA_TO_CODONS.get(curr_aa, [curr_codon])

                for p_alt in prev_alts:
                    for c_alt in curr_alts:
                        if p_alt[-1] + c_alt[0] != "GT":
                            if "GT" not in p_alt and "GT" not in c_alt:
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


def check_conservation_score(original_aa: str, new_aa: str, min_score: int = 0) -> PredicateResult:
    """Predicate 7: BLOSUM62 conservation score for amino acid substitution."""
    score = BLOSUM62.get((original_aa, new_aa), -10)
    passed = score >= min_score
    return PredicateResult("ConservationScore", passed, verdict=Verdict.PASS if passed else SpliceVerdict.FAIL,
                           details=f"BLOSUM62({original_aa},{new_aa})={score}, min={min_score}")


def check_codon_optimality(codon: str, species_cai: Dict[str, float], min_cai: float = 0.0) -> PredicateResult:
    """Predicate 8: Codon optimality (CAI score above threshold)."""
    cai = species_cai.get(codon, 0.0)
    passed = cai >= min_cai
    return PredicateResult("CodonOptimality", passed, verdict=Verdict.PASS if passed else SpliceVerdict.FAIL,
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


# ════════════════════════════════════════════════════════════════
# High-level evaluate_* API — returns TypeCheckResult objects
# ════════════════════════════════════════════════════════════════

from .types import TypeCheckResult


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
) -> TypeCheckResult:
    """Evaluate whether the sequence contains cryptic splice sites.

    Uses dual-threshold MaxEntScan scoring: sites scoring >= cryptic_threshold
    are FAIL, sites scoring >= uncertain_lo (but < cryptic_threshold) are
    UNCERTAIN, and sites below uncertain_lo are PASS.

    When uncertain_lo=0 (default), only PASS/FAIL verdicts are produced,
    preserving backward compatibility.

    Args:
        seq: DNA sequence to evaluate.
        boundaries: Exon boundaries (used for context; currently informational).
        cryptic_threshold: Score threshold above which a site is FAIL.
        uncertain_lo: Score threshold above which a site is UNCERTAIN.
            Set to 0 to disable UNCERTAIN zone (binary PASS/FAIL only).

    Returns:
        TypeCheckResult with PASS/UNCERTAIN/FAIL verdict.
    """
    from .maxentscan import score_donor

    seq = seq.upper()
    gt_positions = [i for i in range(len(seq) - 1) if seq[i:i+2] == "GT"]

    if not gt_positions:
        return TypeCheckResult(
            predicate="NoCrypticSplice",
            verdict=Verdict.PASS,
        )

    worst_score = -50.0
    worst_pos = -1
    worst_verdict = Verdict.PASS

    for pos in gt_positions:
        score = score_donor(seq, pos)
        # If score_donor returns -50 (insufficient context), treat as 0
        if score <= -50.0:
            score = 0.0

        if score >= cryptic_threshold:
            v = Verdict.FAIL
        elif uncertain_lo > 0 and score >= uncertain_lo:
            v = Verdict.UNCERTAIN
        else:
            v = Verdict.PASS

        if score > worst_score:
            worst_score = score
            worst_pos = pos
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
        cellular_context: Cell type context (currently informational).

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    seq = seq.upper()

    if not boundaries or len(boundaries) < 2:
        # Single-exon gene or no boundaries provided — nothing to check
        return TypeCheckResult(
            predicate="SpliceCorrect",
            verdict=Verdict.PASS,
        )

    # Check canonical splice signals at each intron boundary
    for i in range(len(boundaries) - 1):
        intron_start = boundaries[i][1]
        intron_end = boundaries[i + 1][0]

        if intron_start >= intron_end:
            continue

        # Check donor (GT) at intron start
        if intron_start + 2 <= len(seq):
            donor = seq[intron_start:intron_start + 2]
            if donor != "GT":
                return TypeCheckResult(
                    predicate="SpliceCorrect",
                    verdict=Verdict.FAIL,
                    violation=f"Non-canonical donor {donor} at pos {intron_start}",
                )

        # Check acceptor (AG) at intron end
        if intron_end - 2 >= 0:
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
) -> TypeCheckResult:
    """Evaluate whether the sequence contains restriction enzyme recognition sites.

    Args:
        seq: DNA sequence to evaluate.
        enzymes: List of enzyme names to check for.
        enzyme_set: Alias for enzymes (used by certificate verification).

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
            if run_len >= 6:
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


def evaluate_no_cpg_island(
    seq: str,
    window: int = 200,
    threshold: float = 0.6,
) -> TypeCheckResult:
    """Evaluate whether the sequence contains CpG islands.

    A CpG island is detected when the observed/expected CG ratio exceeds
    the threshold in any sliding window of the specified size.

    Args:
        seq: DNA sequence to evaluate.
        window: Window size for CpG island scanning.
        threshold: Maximum allowed Obs/Exp CG ratio.

    Returns:
        TypeCheckResult with PASS/FAIL verdict.
    """
    seq = seq.upper()

    if len(seq) < window:
        # Sequence shorter than window — compute for the whole sequence
        c_count = seq.count("C")
        g_count = seq.count("G")
        cg_count = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")
        expected = (c_count * g_count) / len(seq) if len(seq) > 0 else 0
        obs_exp = cg_count / expected if expected > 0 else 0.0
        if obs_exp > threshold:
            return TypeCheckResult(
                predicate="NoCpGIsland",
                verdict=Verdict.FAIL,
                violation=f"CpG island Obs/Exp={obs_exp:.3f} > {threshold}",
            )
        return TypeCheckResult(
            predicate="NoCpGIsland",
            verdict=Verdict.PASS,
        )

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
        return TypeCheckResult(
            predicate="NoCpGIsland",
            verdict=Verdict.FAIL,
            violation=f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f} > {threshold}",
        )
    return TypeCheckResult(
        predicate="NoCpGIsland",
        verdict=Verdict.PASS,
    )


def analyze_codon_at_position(
    seq: str,
    position: int,
    organism: str = "Homo_sapiens",
) -> dict:
    """Analyze the codon at a given position for optimality and alternatives.

    Args:
        seq: DNA coding sequence.
        position: 0-based nucleotide position (will be rounded to codon start).
        organism: Target organism for CAI lookup.

    Returns:
        Dict with keys: codon, amino_acid, cai, alternatives, position.
    """
    from .species import SPECIES

    codon_start = (position // 3) * 3
    if codon_start + 3 > len(seq):
        return {"codon": "N/A", "amino_acid": "N/A", "cai": 0.0, "alternatives": [], "position": codon_start}

    seq = seq.upper()
    codon = seq[codon_start:codon_start + 3]
    aa = CODON_TABLE.get(codon, "?")
    species_cai = SPECIES.get(organism, SPECIES.get("ecoli", {}))
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


def evaluate_all_predicates(
    seq: str,
    boundaries: List[Tuple[int, int]] | None = None,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: List[str] | None = None,
    cryptic_threshold: float = 3.0,
    uncertain_lo: float = 0.0,
    cai_threshold: float = 0.5,
) -> List[TypeCheckResult]:
    """Evaluate all 8 type predicates against a sequence.

    The 8 predicates are:
    1. NoCrypticSplice — no cryptic splice donors
    2. SpliceCorrect — canonical splice signals at intron boundaries
    3. GCInRange — GC content within acceptable range
    4. CodonAdapted — CAI above threshold
    5. NoRestrictionSite — no restriction enzyme sites
    6. InFrame — valid coding frame
    7. NoInstabilityMotif — no mRNA instability motifs
    8. NoCpGIsland — no CpG islands

    Args:
        seq: DNA sequence to evaluate.
        boundaries: Exon boundary tuples [(start, end), ...].
        organism: Target organism for CAI computation.
        gc_lo: Minimum GC fraction.
        gc_hi: Maximum GC fraction.
        enzymes: Restriction enzymes to check. If None, uses a default set
            of common cloning enzymes: EcoRI, BamHI, XhoI, HindIII, NotI.
        cryptic_threshold: MaxEnt score threshold for cryptic splice FAIL.
        uncertain_lo: MaxEnt score threshold for UNCERTAIN.
        cai_threshold: Minimum CAI for CodonAdapted PASS.

    Returns:
        List of 8 TypeCheckResult objects.
    """
    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
    results: List[TypeCheckResult] = [
        evaluate_no_cryptic_splice(seq, boundaries, cryptic_threshold, uncertain_lo),
        evaluate_splice_correct(seq, boundaries),
        evaluate_gc_in_range(seq, gc_lo, gc_hi),
        evaluate_codon_adapted(seq, organism, cai_threshold),
        evaluate_no_restriction_site(seq, enzymes),
        evaluate_in_frame(seq, boundaries),
        evaluate_no_instability_motif(seq),
        evaluate_no_cpg_island(seq),
    ]
    return results


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
# Global registry instance with all 8 predicates registered
# ────────────────────────────────────────────────────────────
registry = PredicateRegistry()

registry.register(
    "NoCrypticSplice",
    evaluate_no_cryptic_splice,
    verify_param_map={
        "seq": "seq",
        "known_exon_boundaries": "boundaries",
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
    },
)
