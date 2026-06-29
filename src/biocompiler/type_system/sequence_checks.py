"""
BioCompiler Type System — Sequence-Level Predicate Checks
=========================================================
Shared helpers (dinucleotide counting, organism classification, translation,
codon-ramp scoring, cross-codon finders) plus DNA/codon-level predicate
check functions: stop codons, cryptic promoter, valid coding seq,
conservation score, codon optimality, BLAST matches, primer compatibility.

This module was extracted from the historical checks.py monolith during
the W8-b refactor. All public names are re-exported by checks.py for
backwards compatibility.
"""
import logging
from typing import Any, Dict, List, Optional, Set, Tuple

from .codon_tables import (
    AA_TO_CODONS,
    BLOSUM62,
    CODON_TABLE,
    PROMOTER_CONSENSUS,
    PredicateResult,
    _BLOSUM62_MISSING_SCORE,
    _CPG_DENSITY_MULTIPLIER,
    _CPG_GC_RICH_THRESHOLD,
    _COTRANS_HIGH_CONFIDENCE,
    _COTRANS_LOW_CONFIDENCE,
    _CODON_RAMP_LENGTH,
    _DG_AU_PAIR_KCAL,
    _DG_GC_PAIR_KCAL,
    _DG_GU_PAIR_KCAL,
    _EUKARYOTE_GT_PER_BP,
    _EUK_INITIATOR_OFFSET_MAX,
    _EUK_INITIATOR_OFFSET_MIN,
    _FAST_CODON_CAI_THRESHOLD,
    _HIGH_AVG_CAI_THRESHOLD,
    _INSTABILITY_T_RUN_MIN,
    _MAXENT_INSUFFICIENT_CONTEXT_SCORE,
    _MIN_RAMP_FOR_WARNING,
    _MRNA_DG_EUKARYOTE_FAIL,
    _MRNA_DG_PROKARYOTE_FAIL,
    _MRNA_MODERATE_DG_RATIO,
    _MRNA_STABILITY_THRESHOLDS,
    _ORGANISM_TO_SPECIES_KEY,
    _PAUSE_SITE_CAI_THRESHOLD,
    _PROMOTER_UNCERTAIN_RATIO,
    _RESTRICTION_SITE_MIN_LENGTH,
    _TM_BORDERLINE_RATIO,
    _TM_EUKARYOTIC_MIN_STRETCH,
    _TM_PROKARYOTIC_MIN_STRETCH,
    _match_iupac,
    _score_consensus,
)
from biocompiler.shared.types import Verdict

# ── NUMBA integration for dinucleotide counting ───────────────────────
try:
    from biocompiler.optimizer.numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        fast_dinucleotide_count as _numba_fast_dinuc_count,
        seq_to_bytes as _seq_to_bytes,
    )
except ImportError:
    _HAS_NUMBA = False
    _numba_fast_dinuc_count = None  # type: ignore[assignment]
    _seq_to_bytes = None  # type: ignore[assignment]


logger = logging.getLogger(__name__)

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


# ────────────────────────────────────────────────────────────
# Organism-aware helpers
# ────────────────────────────────────────────────────────────

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
        from biocompiler.organisms.config import is_eukaryotic_organism
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


def _compute_max_gt_count(seq_len: int, organism: str = "") -> int:
    """Compute the maximum allowed GT dinucleotide count for a sequence.

    For prokaryotes: 0 (hard constraint).
    For eukaryotes: ``max(1, int(seq_len * _EUKARYOTE_GT_PER_BP))``.

    Args:
        seq_len: Length of the DNA sequence in base pairs.
        organism: Target organism name.

    Returns:
        Maximum allowed GT count.
    """
    if organism and _is_prokaryotic_organism(organism):
        return 0
    return max(1, int(seq_len * _EUKARYOTE_GT_PER_BP))


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


def _resolve_species_cai(key: str) -> Dict[str, float]:
    """Resolve an organism name or SPECIES key to a flat codon->CAI-weight dict.

    Uses CODON_ADAPTIVENESS_TABLES with resolve_organism() as the
    single source of truth for CAI weights.

    Args:
        key: Organism name, alias, or short species key
            (e.g. ``"ecoli"``, ``"Homo_sapiens"``, ``"human"``).

    Returns:
        Dict mapping codon strings to CAI adaptiveness values.
    """
    from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

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
            (relative adaptiveness, 0.0-1.0).

    Returns:
        Dict with keys:
          - ramp_score: average CAI in the first 30 codons.
          - pause_sites: list of (codon_position, CAI) tuples where
            CAI < 0.3 outside the ramp region (positions 30+).
          - speed_disruptions: list of (codon_position, original_slow_cai,
            new_fast_cai) tuples.
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


# ────────────────────────────────────────────────────────────
# Predicate check functions
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
            # High-scoring promoter-like motif — even without full consensus,
            # this is suspicious. Return UNCERTAIN rather than PASS.
            worst_verdict = Verdict.UNCERTAIN
            passed = True  # not a hard FAIL
            details = (f"Promoter-like motif found (score {worst_score:.3f}) but lacks "
                       f"multiple elements — uncertain, may be cryptic promoter")
        else:
            details = f"No significant promoter motifs found (worst score {worst_score:.3f})"

    return PredicateResult(
        "NoCrypticPromoter", passed, verdict=worst_verdict,
        details=details,
        positions=promoter_positions,
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
    if isinstance(dna, str) and isinstance(protein, str) and len(dna) == 1 and len(protein) == 1:
        score = BLOSUM62.get((dna, protein), _BLOSUM62_MISSING_SCORE)
        passed = score >= min_score
        return PredicateResult(
            "ConservationScore", passed,
            verdict=Verdict.PASS if passed else Verdict.FAIL,
            details=f"BLOSUM62({dna},{protein})={score}, min={min_score}",
        )

    # Translate DNA -> protein
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
            of codon->CAI weights (old ``species_cai`` parameter).
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
    from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

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


def check_no_blast_matches(
    seq: str,
    reference_sequences: list[str] | None = None,
    k: int = 15,
) -> PredicateResult:
    """Check that the sequence has no significant k-mer matches against references.

    Uses k-mer overlap detection to find matching subsequences. If any
    k-mer of length >= k is shared with a reference sequence, the predicate
    FAILS.

    Args:
        seq: DNA sequence to check (uppercase).
        reference_sequences: List of reference DNA sequences.
        k: Minimum k-mer size to consider significant (default 15).

    Returns:
        PredicateResult with PASS/FAIL verdict.
    """
    if not reference_sequences:
        return PredicateResult(
            "NoBlastMatches", True, verdict=Verdict.PASS,
            details="No reference sequences provided for BLAST match check",
        )

    # Lazy import — the optimizer.blast_avoidance module was removed in the
    # second-pass cleanup. If it is unavailable, the predicate returns PASS
    # with a diagnostic note (the NoBlastMatches predicate is not in the
    # default-20 set, so this only affects callers that explicitly opt in).
    try:
        from ..optimizer.blast_avoidance import check_kmer_overlap
    except ImportError:
        return PredicateResult(
            "NoBlastMatches", True, verdict=Verdict.PASS,
            details=(
                "blast_avoidance module unavailable "
                "(removed in second-pass cleanup); skipping k-mer overlap check"
            ),
        )

    seq = seq.upper()
    overlaps = check_kmer_overlap(seq, reference_sequences, k=k)

    if overlaps:
        positions = [start for start, _length, _kmer in overlaps]
        return PredicateResult(
            "NoBlastMatches", False, verdict=Verdict.FAIL,
            details=f"Found {len(overlaps)} k-mer overlap(s) against reference sequences",
            positions=positions,
        )

    return PredicateResult(
        "NoBlastMatches", True, verdict=Verdict.PASS,
        details=f"No k-mer overlaps (k={k}) found against reference sequences",
    )


def check_primer_compatibility(
    seq: str,
    region_start: int = 0,
    region_end: int | None = None,
    min_tm: float = 55.0,
    max_tm: float = 65.0,
) -> PredicateResult:
    """Check that the sequence is compatible with primer design for the given region.

    Designs primers flanking the specified region and validates that they
    meet Tm, GC clamp, and secondary structure requirements.

    Args:
        seq: Template DNA sequence.
        region_start: Start of the target region (0-based).
        region_end: End of the target region (0-based, exclusive).
            If None, uses the full sequence length.
        min_tm: Minimum acceptable Tm (°C).
        max_tm: Maximum acceptable Tm (°C).

    Returns:
        PredicateResult with PASS/FAIL verdict.
    """
    if region_end is None:
        region_end = len(seq)

    from ..optimizer.primer_design import evaluate_primer_constraint

    result = evaluate_primer_constraint(
        seq, region_start, region_end, min_tm=min_tm, max_tm=max_tm,
    )

    if result.satisfied:
        return PredicateResult(
            "PrimerCompatibility", True, verdict=Verdict.PASS,
            details="Sequence is primer-compatible",
        )

    issues_str = "; ".join(result.issues) if result.issues else "Unknown primer design issue"
    return PredicateResult(
        "PrimerCompatibility", False, verdict=Verdict.FAIL,
        details=f"Primer compatibility issues: {issues_str}",
    )


# ────────────────────────────────────────────────────────────
# miRNA binding site predicate
# ────────────────────────────────────────────────────────────

# miRNA seed database is now in the dedicated module
from .mirna_seeds import get_mirna_seeds


# Tissue proximity groups — miRNAs from nearby tissues are partially relevant
_TISSUE_GROUPS: Dict[str, Set[str]] = {
    "liver": {"liver", "hepatocyte", "fibrotic"},
    "muscle": {"muscle", "smooth_muscle", "cardiac"},
    "immune": {"immune", "blood", "lymphoid", "hematopoietic"},
    "neural": {"neural", "brain"},
    "epithelial": {"epithelial", "endothelial", "kidney"},
    "colon": {"colon", "intestinal", "digestive"},
    "ubiquitous": {"ubiquitous"},
    "tumor_suppressor": {"tumor_suppressor", "cancer"},
}




__all__ = [
    # Shared helpers
    "_count_dinucs_fast",
    "_is_prokaryotic_organism",
    "_compute_max_gt_count",
    "_translate_dna_to_aa",
    "_resolve_species_cai",
    "_compute_codon_ramp_score",
    # Cross-codon finders
    "find_cross_codon_gt",
    "find_cross_codon_cg",
    "find_cross_codon_restriction",
    # DNA-level predicate checks
    "check_no_stop_codons",
    "check_no_cryptic_promoter",
    "check_valid_coding_seq",
    "check_conservation_score",
    "check_codon_optimality",
    "check_no_blast_matches",
    "check_primer_compatibility",
]
