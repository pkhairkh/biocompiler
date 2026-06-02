"""
BioCompiler Type System v7.3.0
==============================
Defines the core types, codon tables, BLOSUM62 matrix, and 28 predicate classes
for certified gene optimization: 12 DNA-level + 4 structure + 4 stability +
4 solubility + 4 immunogenicity.
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
PROMOTER_CONSENSUS: Dict[str, Dict] = {
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
    (TATAAA) and Initiator (YYANWYY).

    Scoring is based on match quality (how many positions match the consensus).
    - If match score >= threshold: FAIL (strong cryptic promoter)
    - If match score >= threshold * 0.8: UNCERTAIN (weak match)
    - Otherwise: PASS
    """
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
            for offset in range(20, 40):
                ini_start = i + offset
                if ini_start + len(initiator) > len(seq):
                    break
                score_ini = _score_consensus(seq[ini_start:ini_start + len(initiator)], initiator)
                combined = (score_tata + score_ini) / 2.0

                if combined > worst_score:
                    worst_score = combined
                    worst_pos = i
                    promoter_positions = [i, ini_start]

    # Determine verdict based on worst score
    if worst_score >= threshold:
        worst_verdict = Verdict.FAIL
    elif worst_score >= threshold * 0.8:
        worst_verdict = Verdict.UNCERTAIN
    else:
        worst_verdict = Verdict.PASS

    passed = worst_verdict != Verdict.FAIL
    details = f"Worst promoter score {worst_score:.3f} at pos {worst_pos}"
    if worst_verdict == Verdict.PASS:
        details = f"No significant promoter motifs found (worst score {worst_score:.3f})"

    return PredicateResult(
        "NoCrypticPromoter", passed, verdict=worst_verdict,
        details=details,
        positions=promoter_positions,
    )


def check_no_cryptic_splice(seq: str, low_thresh: float = 3.0, high_thresh: float = 6.0) -> PredicateResult:
    """Predicate 2: No cryptic splice sites (dual-threshold PASS/UNCERTAIN/FAIL)."""
    from .splicing import maxent_score
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


def check_no_unexpected_tm_domain(
    seq: str,
    is_cytosolic: bool = True,
    window_size: int = 19,
    threshold: float = 0.68,
) -> PredicateResult:
    """Predicate 10: No unexpected transmembrane (TM) domains after mutagenesis.

    If a cytosolic protein gains hydrophobic stretches from amino acid
    substitutions, that constitutes a FAIL. Transmembrane domains are
    detected by sliding a window of `window_size` amino acids and computing
    the fraction of hydrophobic residues (A, V, I, L, M, F, W, Y).

    Verdict logic (only applies when is_cytosolic=True):
    - If any window exceeds `threshold`: FAIL
    - If any window exceeds `threshold * 0.85`: UNCERTAIN
    - Otherwise: PASS

    If is_cytosolic=False (membrane protein), TM domains are expected: PASS.

    Args:
        seq: DNA sequence to evaluate.
        is_cytosolic: Whether the protein is cytosolic (default True).
        window_size: Sliding window size in amino acids (default 19).
        threshold: Hydrophobic fraction threshold for FAIL (default 0.68).

    Returns:
        PredicateResult with PASS/UNCERTAIN/FAIL verdict.
    """
    seq = seq.upper()

    if not is_cytosolic:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.PASS,
            details="Membrane protein — TM domains are expected",
        )

    # Translate DNA to amino acids
    aa_seq = _translate_dna_to_aa(seq)

    if len(aa_seq) < window_size:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.PASS,
            details=f"Protein too short for TM domain scan ({len(aa_seq)} aa < {window_size} window)",
        )

    HYDROPHOBIC = set("AVILMFWY")
    worst_frac = 0.0
    worst_pos = -1
    borderline_positions: List[int] = []
    fail_positions: List[int] = []

    for i in range(len(aa_seq) - window_size + 1):
        window = aa_seq[i:i + window_size]
        hydro_count = sum(1 for aa in window if aa in HYDROPHOBIC)
        frac = hydro_count / window_size
        if frac > worst_frac:
            worst_frac = frac
            worst_pos = i
        if frac > threshold * 0.85:
            borderline_positions.append(i)
        if frac > threshold:
            fail_positions.append(i)

    if fail_positions:
        return PredicateResult(
            "NoUnexpectedTMDomain", False, verdict=Verdict.FAIL,
            details=(f"TM domain detected: worst hydrophobic fraction {worst_frac:.3f} "
                     f"at AA pos {worst_pos} exceeds threshold {threshold} "
                     f"({len(fail_positions)} window(s) failing)"),
            positions=fail_positions,
        )

    if borderline_positions:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.UNCERTAIN,
            details=(f"Borderline TM domain: worst hydrophobic fraction {worst_frac:.3f} "
                     f"at AA pos {worst_pos} exceeds {threshold * 0.85:.3f} "
                     f"({len(borderline_positions)} window(s) borderline)"),
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
) -> PredicateResult:
    """Predicate 11: No strong mRNA secondary structure around RBS/start codon.

    Checks for stable secondary structure near the ribosome binding site
    that could block ribosome binding. Uses a simplified nearest-neighbor
    ΔG approximation based on counting potential base pairs in the window.

    Scoring (simplified ViennaRNA-style):
    - Count GC pairs (G-C and C-G) in the window
    - Count AU pairs (A-U and U-A) in the window
    - Count GU wobble pairs (G-U and U-G) in the window
    - ΔG ≈ -1.5 * gc_pairs - 0.5 * au_pairs - 0.3 * gu_pairs

    Verdict logic:
    - If ΔG <= dg_threshold (very stable structure): FAIL
    - If ΔG <= dg_threshold * 0.7: UNCERTAIN
    - Otherwise: PASS

    Args:
        seq: DNA sequence to evaluate.
        window_start: Start position of the analysis window (default 0).
        window_end: End position of the analysis window (default 50).
        dg_threshold: ΔG threshold for FAIL (default -15.0 kcal/mol).

    Returns:
        PredicateResult with PASS/UNCERTAIN/FAIL verdict.
    """
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
    dg = -1.5 * gc_pairs - 0.5 * au_pairs - 0.3 * gu_pairs

    if dg <= dg_threshold:
        return PredicateResult(
            "mRNASecondaryStructure", False, verdict=Verdict.FAIL,
            details=(f"Strong mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
                     f"<= {dg_threshold} (GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
        )

    if dg <= dg_threshold * 0.7:
        return PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
            details=(f"Moderate mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
                     f"<= {dg_threshold * 0.7:.1f} (GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
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


def _compute_codon_ramp_score(seq: str, species_cai: Dict[str, float]) -> Dict:
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
    ramp_length = min(30, num_codons)  # first 30 codons = ramp region

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
        if cai < 0.3:
            pause_sites.append((idx, cai))

    # Speed disruptions: fast codons that likely replaced slow ones
    # A speed disruption is detected when the current codon is fast (CAI > 0.7)
    # but the slowest synonymous codon for the same AA has CAI < 0.3,
    # suggesting an optimization replaced a natural pause site.
    speed_disruptions: List[Tuple[int, float, float]] = []
    for idx, cai in codon_cais:
        if idx < ramp_length:
            continue  # ramp region — speed-ups are expected there
        if cai <= 0.7:
            continue  # not a fast codon
        codon = seq[idx * 3:(idx + 1) * 3]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            continue
        # Check if the slowest synonymous codon for this AA is a pause site
        syn_codons = AA_TO_CODONS.get(aa, [codon])
        slowest_cai = min(species_cai.get(c, 0.0) for c in syn_codons)
        if slowest_cai < 0.3:
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
) -> PredicateResult:
    """Predicate 12: Check co-translational folding preservation.

    Checks whether codon optimization has disrupted critical pause sites
    that are important for proper protein folding during translation.
    Slow codons (low tRNA adaptation) at domain boundaries allow the
    nascent chain to fold properly before downstream sequence emerges.

    Args:
        seq: DNA coding sequence (uppercase).
        species_cai: Dict mapping codon strings to CAI values.
        domain_boundaries: Optional list of codon positions (0-based)
            where protein domains start/end. If provided, these positions
            are checked for speed-up (CAI > 0.7 where a pause is needed).
        min_pause_cai: CAI threshold below which a codon is considered
            a pause site (default 0.3).

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
    ramp_length = min(30, num_codons)
    ramp_all_fast = all(
        species_cai.get(seq[i * 3:(i + 1) * 3], 0.0) > 0.7
        for i in range(ramp_length)
    )

    if ramp_all_fast and ramp_length >= 10:
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
            if boundary_cai > 0.7:
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

        if avg_cai > 0.9 and not pause_sites:
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
        # Single domain boundary disrupted
        verdict = Verdict.UNCERTAIN
        passed = True
    elif ramp_all_fast and ramp_length >= 10:
        # Ramp too fast — ribosome jam risk
        verdict = Verdict.UNCERTAIN
        passed = True
    elif speed_disruptions:
        # Some pause sites may have been replaced by fast codons
        verdict = Verdict.UNCERTAIN
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

    details = "; ".join(details_parts) if details_parts else "Co-translational folding appears preserved"

    return PredicateResult(
        "CoTranslationalFolding",
        passed,
        verdict=verdict,
        details=details,
        positions=flagged_positions,
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


def evaluate_no_unexpected_tm_domain(
    seq: str,
    is_cytosolic: bool = True,
    window_size: int = 19,
    threshold: float = 0.68,
) -> TypeCheckResult:
    """Evaluate whether a cytosolic protein has gained unexpected TM domains.

    Checks for hydrophobic stretches in the translated protein that could
    form transmembrane domains — a common side-effect of mutagenesis on
    cytosolic proteins.

    Verdicts use the five-valued logic:
    - LIKELY_PASS: No hydrophobic stretches detected and protein is cytosolic
    - UNCERTAIN: Borderline hydrophobic stretches (fraction between
      threshold*0.85 and threshold)
    - FAIL: Clear TM domains in a cytosolic protein (fraction > threshold)
    - PASS: If not cytosolic (membrane protein), TM domains are expected

    Args:
        seq: DNA sequence to evaluate.
        is_cytosolic: Whether the protein is cytosolic (default True).
        window_size: Sliding window size in amino acids (default 19).
        threshold: Hydrophobic fraction threshold for FAIL (default 0.68).

    Returns:
        TypeCheckResult with LIKELY_PASS/UNCERTAIN/FAIL/PASS verdict.
    """
    seq = seq.upper()

    if not is_cytosolic:
        return TypeCheckResult(
            predicate="NoUnexpectedTMDomain",
            verdict=Verdict.PASS,
        )

    # Translate DNA to amino acids
    aa_seq = _translate_dna_to_aa(seq)

    if len(aa_seq) < window_size:
        return TypeCheckResult(
            predicate="NoUnexpectedTMDomain",
            verdict=Verdict.LIKELY_PASS,
        )

    HYDROPHOBIC = set("AVILMFWY")
    worst_frac = 0.0
    worst_pos = -1

    for i in range(len(aa_seq) - window_size + 1):
        window = aa_seq[i:i + window_size]
        hydro_count = sum(1 for aa in window if aa in HYDROPHOBIC)
        frac = hydro_count / window_size
        if frac > worst_frac:
            worst_frac = frac
            worst_pos = i

    # Determine verdict using five-valued logic
    if worst_frac > threshold:
        verdict = Verdict.FAIL
    elif worst_frac > threshold * 0.85:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_PASS

    violation = None
    if verdict == Verdict.FAIL:
        violation = (
            f"TM domain at AA pos {worst_pos}, hydrophobic fraction "
            f"{worst_frac:.3f} > {threshold}"
        )
    elif verdict == Verdict.UNCERTAIN:
        violation = (
            f"Borderline TM domain at AA pos {worst_pos}, hydrophobic fraction "
            f"{worst_frac:.3f} > {threshold * 0.85:.3f}"
        )

    return TypeCheckResult(
        predicate=f"NoUnexpectedTMDomain({is_cytosolic}, {threshold})",
        verdict=verdict,
        violation=violation,
    )


def evaluate_mrna_secondary_structure(
    seq: str,
    window_start: int = 0,
    window_end: int = 50,
    dg_threshold: float = -15.0,
) -> TypeCheckResult:
    """Evaluate mRNA secondary structure stability around the RBS/start codon.

    Strong secondary structure in the 5' UTR or around the start codon can
    block ribosome binding and reduce translation efficiency. This predicate
    uses a simplified nearest-neighbor ΔG approximation.

    Verdicts use the five-valued logic:
    - LIKELY_PASS: Weak structure (ΔG close to 0, easy for ribosome access)
    - UNCERTAIN: Moderate structure (ΔG between dg_threshold*0.7 and dg_threshold)
    - FAIL: Very stable structure that blocks ribosome binding
      (ΔG <= dg_threshold)

    Args:
        seq: DNA sequence to evaluate.
        window_start: Start position of the analysis window (default 0).
        window_end: End position of the analysis window (default 50).
        dg_threshold: ΔG threshold for FAIL (default -15.0 kcal/mol).

    Returns:
        TypeCheckResult with LIKELY_PASS/UNCERTAIN/FAIL verdict.
    """
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
    dg = -1.5 * gc_pairs - 0.5 * au_pairs - 0.3 * gu_pairs

    if dg <= dg_threshold:
        verdict = Verdict.FAIL
    elif dg <= dg_threshold * 0.7:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_PASS

    violation = None
    if verdict == Verdict.FAIL:
        violation = (
            f"Strong mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
            f"<= {dg_threshold}"
        )
    elif verdict == Verdict.UNCERTAIN:
        violation = (
            f"Moderate mRNA secondary structure: ΔG={dg:.1f} kcal/mol "
            f"<= {dg_threshold * 0.7:.1f}"
        )

    return TypeCheckResult(
        predicate=f"mRNASecondaryStructure({window_start}, {window_end}, {dg_threshold})",
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
    scans for TATA box and Initiator.

    Verdicts use the five-valued logic:
    - LIKELY_PASS: No significant promoter motifs (best score well below threshold)
    - UNCERTAIN: Borderline match (score between threshold*0.8 and threshold)
    - FAIL: Strong cryptic promoter match (score >= threshold)

    Args:
        seq: DNA sequence to evaluate.
        organism: Organism whose promoter consensus to use (default: "E_coli").
        threshold: Score threshold above which a match is FAIL.

    Returns:
        TypeCheckResult with LIKELY_PASS/UNCERTAIN/FAIL verdict.
    """
    seq = seq.upper()
    if len(seq) < 6:
        return TypeCheckResult(
            predicate="NoCrypticPromoter",
            verdict=Verdict.LIKELY_PASS,
        )

    consensus_info = PROMOTER_CONSENSUS.get(organism, PROMOTER_CONSENSUS["E_coli"])
    worst_score = 0.0
    worst_pos = -1

    if consensus_info["type"] == "prokaryotic":
        box35 = consensus_info["-35_box"]
        box10 = consensus_info["-10_box"]
        spacer = consensus_info["spacer"]
        promoter_len = len(box35) + spacer + len(box10)

        for i in range(len(seq) - promoter_len + 1):
            region_35 = seq[i:i + len(box35)]
            region_10 = seq[i + len(box35) + spacer:i + promoter_len]

            if len(region_10) < len(box10):
                continue

            score_35 = _score_consensus(region_35, box35)
            score_10 = _score_consensus(region_10, box10)
            combined = (score_35 + score_10) / 2.0

            if combined > worst_score:
                worst_score = combined
                worst_pos = i

    elif consensus_info["type"] == "eukaryotic":
        tata_box = consensus_info["TATA_box"]
        initiator = consensus_info["Initiator"]

        for i in range(len(seq) - len(tata_box) + 1):
            score_tata = _score_consensus(seq[i:i + len(tata_box)], tata_box)

            for offset in range(20, 40):
                ini_start = i + offset
                if ini_start + len(initiator) > len(seq):
                    break
                score_ini = _score_consensus(seq[ini_start:ini_start + len(initiator)], initiator)
                combined = (score_tata + score_ini) / 2.0

                if combined > worst_score:
                    worst_score = combined
                    worst_pos = i

    # Determine verdict using five-valued logic
    if worst_score >= threshold:
        verdict = Verdict.FAIL
    elif worst_score >= threshold * 0.8:
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.LIKELY_PASS

    violation = None
    if verdict == Verdict.FAIL:
        violation = f"Cryptic promoter at pos {worst_pos}, score={worst_score:.3f} >= {threshold}"
    elif verdict == Verdict.UNCERTAIN:
        violation = f"Borderline promoter at pos {worst_pos}, score={worst_score:.3f} >= {threshold * 0.8:.3f}"

    return TypeCheckResult(
        predicate=f"NoCrypticPromoter({organism}, {threshold})",
        verdict=verdict,
        violation=violation,
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
    from .organisms import SPECIES

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


def evaluate_co_translational_folding(
    seq: str,
    organism: str = "Homo_sapiens",
    domain_boundaries: List[int] | None = None,
    min_pause_cai: float = 0.3,
) -> TypeCheckResult:
    """Evaluate co-translational folding preservation after codon optimization.

    Checks whether codon optimization has disrupted critical pause sites
    that are important for proper protein folding during translation. The
    key insight is that slow codons (low tRNA adaptation) at domain
    boundaries allow the nascent chain to fold properly before downstream
    sequence emerges from the ribosome.

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

    Returns:
        TypeCheckResult with PASS/LIKELY_PASS/UNCERTAIN/LIKELY_FAIL/FAIL verdict.
    """
    from .organisms import SPECIES

    seq = seq.upper()

    # Resolve organism name to SPECIES dict key
    species_key = _ORGANISM_TO_SPECIES_KEY.get(organism, "human")
    species_cai = SPECIES.get(species_key, SPECIES.get("human", {}))

    # Run the low-level predicate check
    result = check_co_translational_folding(
        seq, species_cai, domain_boundaries, min_pause_cai
    )

    # Map the PredicateResult verdict to a TypeCheckResult verdict using
    # the five-valued logic, incorporating additional context about the
    # severity of the findings.
    ramp_info = _compute_codon_ramp_score(seq, species_cai)
    ramp_score = ramp_info["ramp_score"]
    speed_disruptions = ramp_info["speed_disruptions"]
    pause_sites = ramp_info["pause_sites"]

    ramp_length = min(30, len(seq) // 3)
    ramp_all_fast = ramp_length >= 10 and all(
        species_cai.get(seq[i * 3:(i + 1) * 3], 0.0) > 0.7
        for i in range(ramp_length)
    )

    # Count domain boundary disruptions
    domain_disrupted = 0
    if domain_boundaries:
        num_codons = len(seq) // 3
        for bp in domain_boundaries:
            if 0 <= bp < num_codons:
                codon = seq[bp * 3:(bp + 1) * 3]
                if species_cai.get(codon, 0.0) > 0.7:
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
    elif ramp_all_fast and ramp_length >= 10:
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

    Returns:
        List of 12 TypeCheckResult objects.
    """
    # Backward compatibility: accept known_exon_boundaries as alias for boundaries
    if boundaries is None and known_exon_boundaries is not None:
        boundaries = known_exon_boundaries
    if boundaries is None:
        boundaries = [(0, len(seq))]

    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]
    # Map organism name to promoter organism key
    promoter_organism = "E_coli" if organism in ("E_coli", "ecoli", "Escherichia_coli") else "eukaryote"
    results: List[TypeCheckResult] = [
        evaluate_no_cryptic_splice(seq, boundaries, cryptic_threshold, uncertain_lo),
        evaluate_splice_correct(seq, boundaries),
        evaluate_gc_in_range(seq, gc_lo, gc_hi),
        evaluate_codon_adapted(seq, organism, cai_threshold),
        evaluate_no_restriction_site(seq, enzymes),
        evaluate_in_frame(seq, boundaries),
        evaluate_no_instability_motif(seq),
        evaluate_no_cpg_island(seq),
        evaluate_no_cryptic_promoter(seq, promoter_organism, promoter_threshold),
        evaluate_no_unexpected_tm_domain(seq, is_cytosolic, 19, tm_threshold),
        evaluate_mrna_secondary_structure(seq, 0, mrna_window, mrna_dg_threshold),
        evaluate_co_translational_folding(seq, organism, domain_boundaries, folding_threshold),
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
# Global registry instance with all 12 predicates registered
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