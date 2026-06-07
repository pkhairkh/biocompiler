"""
BioCompiler Type System — Predicate Check Functions
====================================================
Low-level check_* functions that return PredicateResult objects.
Each check function implements the core logic for a single predicate.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

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
from ..types import Verdict

logger = logging.getLogger(__name__)


# ── NUMBA integration for dinucleotide counting ───────────────────────
try:
    from ..numba_kernels import (
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
        from ..organism_config import is_eukaryotic_organism
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
      - score < low_thresh  -> PASS
      - low_thresh <= score < high_thresh -> UNCERTAIN
      - score >= high_thresh -> FAIL

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

    from ..maxentscan import score_donor, score_acceptor

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
    from ..restriction_sites import get_recognition_site
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
                for c_alt in AA_TO_CODONS.get(curr_aa, [curr_codon]):
                    if prev_codon[-1] + c_alt[0] != "GT":
                        has_avoidable = True
                        break
            elif curr_aa == "*":
                for p_alt in AA_TO_CODONS.get(prev_aa, [prev_codon]):
                    if p_alt[-1] + curr_codon[0] != "GT":
                        has_avoidable = True
                        break
            else:
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

    - **Eukaryotes**: Soft constraint — GTs are reported but the predicate
      uses ``LIKELY_FAIL`` (not ``FAIL``) when GTs exceed ``max_gt_count``,
      indicating a soft violation that should not block the optimization.
      The predicate PASSES (``PASS``) if GT count <= ``max_gt_count``.

    - ``max_gt_count``: If not provided, auto-computed from sequence length
      using :func:`_compute_max_gt_count`.

    Args:
        seq: DNA sequence to evaluate.
        organism: Target organism name. If prokaryotic, any GT is FAIL.
        max_gt_count: Maximum GT count before triggering SOFT_FAIL for
            eukaryotes. If None, auto-computed from sequence length.

    Returns:
        PredicateResult with verdict.
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
                f"GT dinucleotides: {gt_count} <= max_gt_count={max_gt_count} "
                f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
                f"Acceptable for eukaryotes — in-codon GTs from optimal codons "
                f"are biologically common."
            ),
            positions=gt_positions,
        )

    # GT count exceeds tolerance — soft fail
    return PredicateResult(
        "NoGTDinucleotide", True, verdict=Verdict.LIKELY_FAIL,
        details=(
            f"GT dinucleotides: {gt_count} > max_gt_count={max_gt_count} "
            f"(in-codon: {len(in_codon_gt)}, cross-codon: {len(cross_codon_gt)}). "
            f"Soft constraint for eukaryotes: in-codon GTs from optimal codons "
            f"are acceptable (CAI > GT avoidance). Consider if these GTs form "
            f"strong cryptic splice donors (MaxEntScan score >= threshold)."
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
    (positive-inside rule).

    Verdict logic (only applies when is_cytosolic=True):
    - If any window exceeds `threshold` AND has TM-like flanking charges: FAIL
    - If any window exceeds `threshold` but lacks TM flanking charges: UNCERTAIN
    - If any window exceeds `threshold * _TM_BORDERLINE_RATIO`: UNCERTAIN
    - Otherwise: PASS

    If is_cytosolic=False (membrane protein), TM domains are expected: PASS.
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
            flank_n = aa_seq[max(0, i - 5):i]
            flank_c = aa_seq[i + window_size:min(len(aa_seq), i + window_size + 5)]
            n_pos_count = sum(1 for aa in flank_n if aa in POSITIVE)
            c_pos_count = sum(1 for aa in flank_c if aa in POSITIVE)
            has_tm_flanking = (n_pos_count >= 1 or c_pos_count >= 1) and not (n_pos_count >= 2 and c_pos_count >= 2)
            if has_tm_flanking:
                fail_positions.append(i)
            else:
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

    if fail_no_flank_positions:
        return PredicateResult(
            "NoUnexpectedTMDomain", True, verdict=Verdict.LIKELY_PASS,
            details=(f"Hydrophobic stretch without TM flanking charges: "
                     f"worst fraction {worst_frac:.3f} at AA pos {worst_pos} "
                     f"({len(fail_no_flank_positions)} window(s) — likely soluble protein patch, not TM domain)"),
            positions=fail_no_flank_positions,
        )

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

    Organism-specific dG cutoffs:
    - Prokaryotes: dG < -15 kcal/mol is FAIL
    - Eukaryotes: dG < -25 kcal/mol is FAIL
    - If organism is not specified, uses the provided dg_threshold
    """
    seq = seq.upper()

    # Apply organism-specific dG threshold
    effective_threshold = dg_threshold
    if organism:
        if _is_prokaryotic_organism(organism):
            effective_threshold = _MRNA_DG_PROKARYOTE_FAIL
        else:
            effective_threshold = _MRNA_DG_EUKARYOTE_FAIL

    # Try ViennaRNA/Nussinov for accurate dG
    if use_viennarna:
        try:
            from ..viennarna import is_viennarna_available, compute_5prime_dg
            if is_viennarna_available():
                dg = compute_5prime_dg(seq, window=window_end - window_start)
                if dg <= effective_threshold:
                    return PredicateResult(
                        "mRNASecondaryStructure", False, verdict=Verdict.FAIL,
                        details=(f"Strong mRNA secondary structure (ViennaRNA): "
                                 f"dG={dg:.1f} kcal/mol <= {effective_threshold}"),
                    )
                if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
                    return PredicateResult(
                        "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
                        details=(f"Moderate mRNA secondary structure (ViennaRNA): "
                                 f"dG={dg:.1f} kcal/mol <= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f}"),
                    )
                return PredicateResult(
                    "mRNASecondaryStructure", True, verdict=Verdict.PASS,
                    details=(f"Weak mRNA secondary structure (ViennaRNA): "
                             f"dG={dg:.1f} kcal/mol"),
                )
        except Exception:
            logger.debug("Falling through to Nussinov fallback", exc_info=True)

        try:
            from ..viennarna_fallback import compute_approx_dg
            dg = compute_approx_dg(seq, region="5utr")
            if dg <= effective_threshold:
                return PredicateResult(
                    "mRNASecondaryStructure", False, verdict=Verdict.FAIL,
                    details=(f"Strong mRNA secondary structure (Nussinov fallback): "
                             f"dG~{dg:.1f} kcal/mol <= {effective_threshold}"),
                )
            if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
                return PredicateResult(
                    "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
                    details=(f"Moderate mRNA secondary structure (Nussinov fallback): "
                             f"dG~{dg:.1f} kcal/mol <= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f}"),
                )
            return PredicateResult(
                "mRNASecondaryStructure", True, verdict=Verdict.PASS,
                details=(f"Weak mRNA secondary structure (Nussinov fallback): "
                         f"dG~{dg:.1f} kcal/mol"),
            )
        except Exception:
            logger.debug("Falling through to legacy toy model", exc_info=True)

    # Legacy toy model (original implementation, backward-compatible)
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

    # Simplified nearest-neighbor dG estimate
    dg = _DG_GC_PAIR_KCAL * gc_pairs + _DG_AU_PAIR_KCAL * au_pairs + _DG_GU_PAIR_KCAL * gu_pairs

    if dg <= effective_threshold:
        return PredicateResult(
            "mRNASecondaryStructure", False, verdict=Verdict.FAIL,
            details=(f"Strong mRNA secondary structure: dG={dg:.1f} kcal/mol "
                     f"<= {effective_threshold} (GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
        )

    if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
        window_gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) if window_seq else 0.5
        if window_gc < 0.5:
            return PredicateResult(
                "mRNASecondaryStructure", True, verdict=Verdict.LIKELY_PASS,
                details=(f"Moderate mRNA secondary structure: dG={dg:.1f} kcal/mol, "
                         f"but AT-rich window (GC={window_gc:.0%}) weakens structure "
                         f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
            )
        return PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.LIKELY_FAIL,
            details=(f"Moderate mRNA secondary structure: dG={dg:.1f} kcal/mol "
                     f"<= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f} "
                     f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
        )

    return PredicateResult(
        "mRNASecondaryStructure", True, verdict=Verdict.PASS,
        details=(f"Weak mRNA secondary structure: dG={dg:.1f} kcal/mol "
                 f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs)"),
    )


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

    When structure_confidence (ESMFold confidence or similar) is provided,
    UNCERTAIN verdicts can be resolved:
    - structure_confidence > 0.7 and pLDDT is good -> resolve to PASS
    - structure_confidence < 0.5 -> resolve to FAIL
    - Otherwise keep UNCERTAIN
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
                domain_disrupted += 1
                flagged_positions.append(boundary_pos)
                details_parts.append(
                    f"Domain boundary at codon {boundary_pos} has CAI={boundary_cai:.3f} "
                    f"(fast codon where pause is needed)"
                )
    else:
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
        verdict = Verdict.FAIL
        passed = False
    elif domain_boundaries and domain_disrupted >= 2:
        verdict = Verdict.LIKELY_FAIL
        passed = False
    elif domain_boundaries and domain_disrupted == 1:
        verdict = Verdict.LIKELY_FAIL
        passed = True
    elif ramp_all_fast and ramp_length >= _MIN_RAMP_FOR_WARNING:
        if ramp_score > 0.95:
            verdict = Verdict.LIKELY_FAIL
        else:
            verdict = Verdict.LIKELY_PASS
        passed = True
    elif speed_disruptions:
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
            verdict = Verdict.PASS
            passed = True
            details_parts.append(
                f"UNCERTAIN resolved to PASS: structure_confidence={structure_confidence:.3f} > {_COTRANS_HIGH_CONFIDENCE}"
            )
        elif structure_confidence < _COTRANS_LOW_CONFIDENCE:
            verdict = Verdict.LIKELY_FAIL
            passed = False
            details_parts.append(
                f"UNCERTAIN resolved to LIKELY_FAIL: structure_confidence={structure_confidence:.3f} < {_COTRANS_LOW_CONFIDENCE}"
            )

    details = "; ".join(details_parts) if details_parts else "Co-translational folding appears preserved"

    return PredicateResult(
        "CoTranslationalFolding",
        passed,
        verdict=verdict,
        details=details,
        positions=flagged_positions,
    )


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
    """
    from ..mrna_stability import score_mrna_stability

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


# ────────────────────────────────────────────────────────────
# BLAST match avoidance and primer compatibility checks
# ────────────────────────────────────────────────────────────

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

    from ..optimizer.blast_avoidance import check_kmer_overlap

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
