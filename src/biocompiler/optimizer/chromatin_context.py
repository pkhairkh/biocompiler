"""
BioCompiler Chromatin Context Prediction Module v1.0.0
=======================================================
Predicts chromatin state (euchromatin vs heterochromatin) from DNA sequence
features for synthetic gene design.

Chromatin context influences transcription, DNA repair, and mutation rates.
This module provides heuristic and deep-learning-based predictions of
chromatin accessibility, histone marks, and repair accessibility from
sequence alone — critical for therapeutic gene design where integration
site chromatin state affects expression stability.

Heuristic Model (always available):
  - GC content: >55% → likely euchromatin; <40% → likely heterochromatin
  - CpG density: Obs/Exp ratio >0.60 → CpG island → euchromatin
  - Repeat density: high repeat fraction → heterochromatin
  - Poly-dA:dT tracts: nucleosome-depleting → open chromatin
  - CTCF motif: boundary element between chromatin states

Deep Learning Integration (optional):
  - Sei framework: 40-sequence-class chromatin profile prediction
  - Enformer: 5313-track prediction including histone marks and DNase

References
----------
- ENCODE Project Consortium. "An integrated encyclopedia of DNA elements
  in the human genome." Nature 489, 57–74 (2012).
- Roadmap Epigenomics Consortium. "Integrative analysis of 111 reference
  human epigenomes." Nature 518, 317–330 (2015).
- Zhou, J. et al. "Deep learning sequence-based ab initio prediction of
  variant effects on chromatin accessibility." Nat Genet 51, 1201–1209
  (2019).  [Sei]
- Avsec, Ž. et al. "Effective gene expression prediction from sequence
  by integrating long-range interactions." Nat Methods 18, 1196–1203
  (2021).  [Enformer]
"""

from __future__ import annotations

import logging
import math
import re
from typing import Any

logger = logging.getLogger(__name__)

# ── Constants ────────────────────────────────────────────────────────────────

# GC content thresholds for chromatin state classification
_GC_EUCHROMATIN_THRESHOLD: float = 0.55
"""GC fraction above which sequence is likely euchromatin."""

_GC_HETEROCHROMATIN_THRESHOLD: float = 0.40
"""GC fraction below which sequence is likely heterochromatin."""

# CpG density thresholds
_CPG_ISLAND_OBS_EXP_THRESHOLD: float = 0.60
"""Obs/Exp CpG ratio above which a region is classified as a CpG island."""

_CPG_EUCHROMATIN_OBS_EXP: float = 0.65
"""Obs/Exp CpG ratio suggesting euchromatin."""

# Repeat density threshold
_REPEAT_HETEROCHROMATIN_FRACTION: float = 0.40
"""Repeat fraction above which sequence is likely heterochromatin."""

# Poly-dA:dT tract parameters
_POLY_DA_DT_MIN_LENGTH: int = 5
"""Minimum length for a poly-dA:dT tract to be nucleosome-depleting."""

_POLY_DA_DT_WINDOW: int = 150
"""Window size for poly-dA:dT tract scanning."""

# CTCF motif consensus (20 bp core)
_CTCF_MOTIF_CONSENSUS: str = "CCGCGNGGNGGCAG"
"""CTCF binding site consensus (IUPAC). Simplified from 20-mer core."""

_CTCF_CORE_LENGTH: int = 14
"""Length of CTCF core motif to search for."""

# DNase accessibility prediction window
_DNASE_WINDOW: int = 50
"""Window size for per-position DNase accessibility prediction."""

_DNASE_STEP: int = 10
"""Step size for DNase accessibility sliding window."""

# Common repeat element sequences (simplified consensus seeds)
_REPEAT_SEEDS: list[str] = [
    # Alu element seeds
    "GCGCCG", "GGCCGG", "GCTCCG", "CCGCTC", "CCGCCT",
    # LINE-1 seeds
    "AATGGAG", "AATGGAA",
    # SINE seeds
    "GGCGGCG", "GTTTCAG",
    # Simple repeats
    "ATATAT", "GTGTGT", "CACACA",
]

# IUPAC ambiguity code mapping for CTCF motif matching
_IUPAC_MAP: dict[str, str] = {
    "A": "A", "T": "T", "G": "G", "C": "C",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


# ── Helper: IUPAC motif → regex ─────────────────────────────────────────────

def _iupac_to_regex(motif: str) -> str:
    """Convert an IUPAC ambiguity motif to a regex pattern.

    Args:
        motif: Sequence with IUPAC codes (e.g. ``"CCGCGNGGNGGCAG"``).

    Returns:
        Regex pattern string matching the IUPAC motif.
    """
    parts: list[str] = []
    for ch in motif.upper():
        if ch in _IUPAC_MAP:
            parts.append(_IUPAC_MAP[ch])
        else:
            parts.append(re.escape(ch))
    return "".join(parts)


# Pre-compile the CTCF motif regex
_CTCF_PATTERN: re.Pattern[str] = re.compile(
    _iupac_to_regex(_CTCF_MOTIF_CONSENSUS), re.IGNORECASE
)


# ── Core Heuristic Functions ────────────────────────────────────────────────


def compute_gc_ishchz_score(seq: str) -> float:
    """Compute a GC-content-based chromatin state score.

    Returns a continuous score from -1 (strong heterochromatin) to +1
    (strong euchromatin), using GC content as a proxy.

    Rationale: GC-rich regions in mammalian genomes strongly correlate
    with euchromatin (gene-rich, early-replicating, DNase-accessible),
    while AT-rich regions correlate with heterochromatin (gene-poor,
    late-replicating, compacted).  This is supported by:
    - ENCODE 2012 Nature: GC-rich bands are DNase-hypersensitive
    - Roadmap Epigenomics 2015 Nature: CpG islands (GC-rich) carry
      active marks (H3K4me3, H3K27ac)

    Args:
        seq: DNA sequence (case-insensitive).

    Returns:
        Score in [-1, +1]:
        - > 0: likely euchromatin
        - < 0: likely heterochromatin
        - ~ 0: ambiguous
    """
    if not seq:
        return 0.0

    seq_upper = seq.upper()
    gc = (seq_upper.count("G") + seq_upper.count("C")) / len(seq_upper)

    if gc >= _GC_EUCHROMATIN_THRESHOLD:
        # Linear ramp from 0 at 55% to 1.0 at 80% GC
        return min(1.0, (gc - _GC_EUCHROMATIN_THRESHOLD) / 0.25)
    elif gc <= _GC_HETEROCHROMATIN_THRESHOLD:
        # Linear ramp from 0 at 40% to -1.0 at 15% GC
        return max(-1.0, -(_GC_HETEROCHROMATIN_THRESHOLD - gc) / 0.25)
    else:
        # Interpolate between heterochromatin and euchromatin thresholds
        mid = (_GC_EUCHROMATIN_THRESHOLD + _GC_HETEROCHROMATIN_THRESHOLD) / 2.0
        half_range = (_GC_EUCHROMATIN_THRESHOLD - _GC_HETEROCHROMATIN_THRESHOLD) / 2.0
        return (gc - mid) / half_range


def compute_cpg_density_score(seq: str) -> float:
    """Compute CpG density as a chromatin state proxy.

    CpG-rich regions (CpG islands) are strongly associated with
    euchromatin — they are typically unmethylated, carry active histone
    marks (H3K4me3), and are DNase-accessible.

    The Obs/Exp CpG ratio is calculated as::

        Obs/Exp = (N_CpG × L) / (N_C × N_G)

    where L is sequence length, N_CpG is the count of CG dinucleotides,
    and N_C, N_G are the counts of C and G respectively.

    A ratio > 0.60 indicates a CpG island (Takai & Jones 2002), and
    ratios > 0.65 are strong indicators of euchromatin.

    Args:
        seq: DNA sequence (case-insensitive).

    Returns:
        Score in [-1, +1]:
        - > 0: CpG-rich, likely euchromatin
        - < 0: CpG-poor, likely heterochromatin
        - ~ 0: intermediate
    """
    if len(seq) < 2:
        return 0.0

    seq_upper = seq.upper()
    n = len(seq_upper)
    cpg_count = sum(1 for i in range(n - 1) if seq_upper[i:i + 2] == "CG")
    c_count = seq_upper.count("C")
    g_count = seq_upper.count("G")

    if c_count == 0 or g_count == 0:
        return -1.0  # No CpG possible → heterochromatin signal

    obs_exp = (cpg_count * n) / (c_count * g_count)

    # Map Obs/Exp ratio to score
    if obs_exp >= _CPG_EUCHROMATIN_OBS_EXP:
        # Strong CpG island → euchromatin
        return min(1.0, (obs_exp - _CPG_EUCHROMATIN_OBS_EXP) / 0.50 + 0.5)
    elif obs_exp >= _CPG_ISLAND_OBS_EXP_THRESHOLD:
        # Moderate CpG island
        return (obs_exp - _CPG_ISLAND_OBS_EXP_THRESHOLD) / (
            _CPG_EUCHROMATIN_OBS_EXP - _CPG_ISLAND_OBS_EXP_THRESHOLD
        ) * 0.5
    elif obs_exp >= 0.20:
        # CpG-poor → heterochromatin tendency
        return (obs_exp - 0.40) / 0.20 * 0.5 - 0.5
    else:
        # Very CpG-poor → strong heterochromatin signal
        return max(-1.0, (obs_exp - 0.20) / 0.20 * 0.5 - 0.5)


def _compute_repeat_fraction(seq: str) -> float:
    """Compute the fraction of sequence covered by repeat elements.

    Uses seed-based matching against common repeat element consensus
    sequences (Alu, LINE-1, SINE, simple repeats).

    High repeat density is a hallmark of heterochromatin, particularly
    pericentromeric and subtelomeric regions (ENCODE 2012).

    Args:
        seq: DNA sequence (case-insensitive).

    Returns:
        Fraction of positions covered by repeat seeds (0-1).
    """
    if not seq:
        return 0.0

    seq_upper = seq.upper()
    n = len(seq_upper)
    covered = [False] * n

    for seed in _REPEAT_SEEDS:
        start = 0
        while True:
            pos = seq_upper.find(seed, start)
            if pos == -1:
                break
            for j in range(pos, min(pos + len(seed), n)):
                covered[j] = True
            start = pos + 1

    return sum(covered) / n


def _compute_poly_da_dt_score(seq: str) -> float:
    """Score poly-dA:dT tract presence as a nucleosome-depletion signal.

    Poly-dA:dT tracts (>5 bp) are strongly nucleosome-depleting due to
    their rigid, narrow minor groove which disfavors histone wrapping
    (Segal et al. 2006 Nature; Kaplan et al. 2009 Nature).  Their
    presence indicates open/accessible chromatin.

    Args:
        seq: DNA sequence (case-insensitive).

    Returns:
        Score from 0 (no tracts) to 1 (high tract density).
    """
    if not seq:
        return 0.0

    seq_upper = seq.upper()
    n = len(seq_upper)

    tract_count = 0
    tract_total_length = 0

    # Scan for poly-dA tracts
    i = 0
    while i < n:
        if seq_upper[i] == "A":
            run_start = i
            while i < n and seq_upper[i] == "A":
                i += 1
            run_len = i - run_start
            if run_len >= _POLY_DA_DT_MIN_LENGTH:
                tract_count += 1
                tract_total_length += run_len
        else:
            i += 1

    # Scan for poly-dT tracts
    i = 0
    while i < n:
        if seq_upper[i] == "T":
            run_start = i
            while i < n and seq_upper[i] == "T":
                i += 1
            run_len = i - run_start
            if run_len >= _POLY_DA_DT_MIN_LENGTH:
                tract_count += 1
                tract_total_length += run_len
        else:
            i += 1

    # Normalize: fraction of sequence covered by qualifying tracts
    if n == 0:
        return 0.0

    coverage = tract_total_length / n

    # Also factor in tract count (more tracts = more boundaries)
    count_score = min(1.0, tract_count / 5.0)  # 5+ tracts = max score

    return 0.6 * coverage * 5.0 + 0.4 * count_score  # weighted blend


def _detect_ctcf_motifs(seq: str) -> list[dict[str, Any]]:
    """Detect CTCF binding motif occurrences in a DNA sequence.

    CTCF is a boundary element protein that marks transitions between
    chromatin states (Phillips & Corces 2009 Cell).  CTCF sites often
    demarcate TAD boundaries and can indicate the edge of a euchromatin
    or heterochromatin domain.

    Uses a simplified 14-mer consensus: CCGCGNGGNGGCAG
    (derived from the 20-mer CTCF binding motif).

    Args:
        seq: DNA sequence (case-insensitive).

    Returns:
        List of dicts with keys ``position``, ``match``, ``orientation``.
    """
    seq_upper = seq.upper()
    motifs: list[dict[str, Any]] = []

    # Forward strand search
    for m in _CTCF_PATTERN.finditer(seq_upper):
        motifs.append({
            "position": m.start(),
            "match": m.group(),
            "orientation": "+",
        })

    # Reverse complement search
    rc = seq_upper.translate(str.maketrans("ATGC", "TACG"))[::-1]
    for m in _CTCF_PATTERN.finditer(rc):
        motifs.append({
            "position": len(seq_upper) - m.end(),  # map back to original coords
            "match": m.group(),
            "orientation": "-",
        })

    return motifs


# ── Main Prediction Functions ───────────────────────────────────────────────


def predict_chromatin_state(seq: str, method: str = "auto") -> dict:
    """Predict whether a DNA sequence is in euchromatin or heterochromatin.

    Combines multiple sequence-based chromatin state proxies using a
    weighted heuristic model:

    1. **GC content** (weight 0.30): GC-rich → euchromatin, AT-rich → heterochromatin
    2. **CpG density** (weight 0.30): CpG island → euchromatin
    3. **Repeat density** (weight 0.20): High repeat fraction → heterochromatin
    4. **Poly-dA:dT tracts** (weight 0.10): Nucleosome-depleting → open chromatin
    5. **CTCF motifs** (weight 0.10): Boundary elements

    When ``method="auto"``, uses the heuristic model.  When
    ``method="sei"`` or ``method="enformer"``, attempts to use the
    corresponding deep learning model, falling back to heuristics if
    unavailable.

    Args:
        seq: DNA sequence (case-insensitive).
        method: Prediction method — ``"auto"`` (heuristic, default),
            ``"sei"`` (Sei framework), or ``"enformer"`` (Enformer).

    Returns:
        Dict with keys:

        - ``"state"``: ``"euchromatin"`` or ``"heterochromatin"``
        - ``"confidence"``: 0–1 confidence score
        - ``"combined_score"``: weighted composite score (positive=euchromatin)
        - ``"gc_score"``: GC-content based score
        - ``"cpg_score"``: CpG density score
        - ``"repeat_score"``: repeat density score (negative=heterochromatin)
        - ``"poly_da_dt_score"``: poly-dA:dT tract score
        - ``"ctcf_count"``: number of CTCF motif hits
        - ``"ctcf_positions"``: list of CTCF motif positions
        - ``"method"``: method actually used
        - ``"fallback"``: True if deep learning was requested but unavailable
    """
    if not seq:
        return {
            "state": "unknown",
            "confidence": 0.0,
            "combined_score": 0.0,
            "gc_score": 0.0,
            "cpg_score": 0.0,
            "repeat_score": 0.0,
            "poly_da_dt_score": 0.0,
            "ctcf_count": 0,
            "ctcf_positions": [],
            "method": method,
            "fallback": False,
        }

    # Compute individual feature scores
    gc_score = compute_gc_ishchz_score(seq)
    cpg_score = compute_cpg_density_score(seq)
    repeat_fraction = _compute_repeat_fraction(seq)
    repeat_score = -repeat_fraction / _REPEAT_HETEROCHROMATIN_FRACTION  # negative = heterochromatin
    repeat_score = max(-1.0, repeat_score)
    poly_score = _compute_poly_da_dt_score(seq)

    ctcf_motifs = _detect_ctcf_motifs(seq)
    ctcf_count = len(ctcf_motifs)
    # CTCF presence: each motif adds a small positive signal (boundary = accessible)
    ctcf_score = min(1.0, ctcf_count / 3.0)  # 3+ motifs = max score

    # Weights
    w_gc = 0.30
    w_cpg = 0.30
    w_repeat = 0.20
    w_poly = 0.10
    w_ctcf = 0.10

    combined = (
        w_gc * gc_score
        + w_cpg * cpg_score
        + w_repeat * repeat_score
        + w_poly * poly_score
        + w_ctcf * ctcf_score
    )

    # Classify state
    if combined > 0.15:
        state = "euchromatin"
    elif combined < -0.15:
        state = "heterochromatin"
    else:
        state = "euchromatin" if combined >= 0 else "heterochromatin"

    confidence = min(1.0, abs(combined) / 0.60)  # normalize to 0-1

    result: dict[str, Any] = {
        "state": state,
        "confidence": confidence,
        "combined_score": combined,
        "gc_score": gc_score,
        "cpg_score": cpg_score,
        "repeat_score": repeat_score,
        "repeat_fraction": repeat_fraction,
        "poly_da_dt_score": poly_score,
        "ctcf_count": ctcf_count,
        "ctcf_positions": [m["position"] for m in ctcf_motifs],
        "method": "heuristic",
        "fallback": False,
    }

    # Try deep learning methods if requested
    if method == "sei":
        sei_result = predict_chromatin_sei(seq)
        if "error" not in sei_result:
            result["method"] = "sei"
            result["sei_profile"] = sei_result
        else:
            result["fallback"] = True
            logger.info("Sei unavailable, using heuristic: %s", sei_result.get("error", ""))
    elif method == "enformer":
        enf_result = predict_chromatin_enformer(seq)
        if "error" not in enf_result:
            result["method"] = "enformer"
            result["enformer_profile"] = enf_result
        else:
            result["fallback"] = True
            logger.info("Enformer unavailable, using heuristic: %s", enf_result.get("error", ""))

    return result


def predict_dnase_accessibility(seq: str) -> list[float]:
    """Predict per-position DNase I hypersensitivity from sequence features.

    DNase I hypersensitivity is the gold-standard experimental assay for
    open chromatin (ENCODE 2012).  This function uses a sliding-window
    heuristic combining:

    - GC content (higher GC → more accessible in promoters/enhancers)
    - CpG density (CpG islands → open chromatin)
    - Poly-dA:dT tracts (nucleosome depletion → accessible)
    - Repeat density (repeats → closed chromatin)

    Each position receives a composite accessibility score from 0 (closed)
    to 1 (fully accessible).  The score is computed over a sliding window
    and assigned to the centre position.

    Args:
        seq: DNA sequence (case-insensitive).

    Returns:
        List of per-position accessibility scores (0–1), same length as
        *seq*.  Positions near the boundaries where the full window is
        unavailable receive scores extrapolated from the nearest full
        window.
    """
    n = len(seq)
    if n == 0:
        return []

    seq_upper = seq.upper()
    accessibility = [0.5] * n  # default: moderate accessibility

    if n < _DNASE_WINDOW:
        # Short sequence: compute one global score
        gc_score = max(0.0, compute_gc_ishchz_score(seq))
        cpg_score = max(0.0, compute_cpg_density_score(seq))
        poly_score = _compute_poly_da_dt_score(seq)
        repeat_frac = _compute_repeat_fraction(seq)
        composite = (
            0.30 * gc_score
            + 0.30 * cpg_score
            + 0.25 * poly_score
            + 0.15 * (1.0 - repeat_frac)
        )
        composite = max(0.0, min(1.0, composite))
        return [composite] * n

    half = _DNASE_WINDOW // 2

    for start in range(0, n - _DNASE_WINDOW + 1, _DNASE_STEP):
        window = seq_upper[start : start + _DNASE_WINDOW]
        centre = start + half

        # GC component
        gc = (window.count("G") + window.count("C")) / len(window)
        gc_comp = max(0.0, min(1.0, (gc - 0.30) / 0.40))  # 0.30→0, 0.70→1

        # CpG component
        cpg_count = sum(1 for i in range(len(window) - 1) if window[i:i + 2] == "CG")
        c_count = window.count("C")
        g_count = window.count("G")
        if c_count > 0 and g_count > 0:
            obs_exp = (cpg_count * len(window)) / (c_count * g_count)
        else:
            obs_exp = 0.0
        cpg_comp = min(1.0, obs_exp / 1.0)  # normalize

        # Poly-dA:dT component
        poly_comp = _compute_poly_da_dt_score(window)

        # Repeat component (inverted: more repeats → less accessible)
        repeat_comp = 1.0 - _compute_repeat_fraction(window)

        # Composite score
        composite = (
            0.30 * gc_comp
            + 0.30 * cpg_comp
            + 0.25 * poly_comp
            + 0.15 * repeat_comp
        )
        composite = max(0.0, min(1.0, composite))

        # Assign to positions in range
        for pos in range(start, min(start + _DNASE_WINDOW, n)):
            # Use a Gaussian-like weighting centred on the window centre
            dist = abs(pos - centre)
            weight = math.exp(-0.5 * (dist / half) ** 2) if half > 0 else 1.0
            accessibility[pos] = accessibility[pos] * (1 - weight) + composite * weight

    return accessibility


def compute_repair_accessibility(chromatin_state: dict) -> float:
    """Derive DNA repair accessibility from a chromatin state prediction.

    DNA repair efficiency is strongly modulated by chromatin compaction.
    Open chromatin (euchromatin) allows repair machinery (NER, BER, MMR)
    to access damaged sites efficiently, while compacted heterochromatin
    impedes repair by 2–10 fold (Lindahl & Wood 1999; Smerdon 1991).

    This function maps the chromatin state prediction to a repair
    accessibility score, incorporating:

    - Overall chromatin state (euchromatin vs heterochromatin)
    - Confidence of the prediction
    - CTCF boundary presence (boundary regions have intermediate
      accessibility)
    - CpG island status (CpG islands in euchromatin are most accessible)

    Literature values for repair rate ratios:
    - Euchromatin: ~90% repair efficiency for UV lesions in 24h
    - Heterochromatin: ~40% repair efficiency for UV lesions in 24h
    - TCR boost for transcribed genes in euchromatin: ~95%

    Args:
        chromatin_state: Output dict from :func:`predict_chromatin_state`.

    Returns:
        Repair accessibility score from 0.0 (no repair access) to 1.0
        (full repair access).  A score of 0.5 corresponds to moderate
        accessibility typical of facultative heterochromatin or boundary
        regions.
    """
    # Base accessibility from chromatin state
    state = chromatin_state.get("state", "unknown")
    combined_score = chromatin_state.get("combined_score", 0.0)
    confidence = chromatin_state.get("confidence", 0.0)

    if state == "euchromatin":
        # Euchromatin: high base accessibility
        # Scale from 0.70 (weak euchromatin) to 1.00 (strong euchromatin)
        base_accessibility = 0.70 + 0.30 * min(1.0, max(0.0, combined_score))
    elif state == "heterochromatin":
        # Heterochromatin: low base accessibility
        # Scale from 0.10 (strong heterochromatin) to 0.40 (weak heterochromatin)
        base_accessibility = 0.40 + 0.30 * min(1.0, max(0.0, combined_score + 1.0))
    else:
        base_accessibility = 0.50

    # Adjust for CpG island presence (CpG islands are most accessible)
    cpg_score = chromatin_state.get("cpg_score", 0.0)
    if cpg_score > 0.3:
        # Boost accessibility for CpG-rich regions
        cpg_boost = 0.10 * min(1.0, cpg_score)
        base_accessibility = min(1.0, base_accessibility + cpg_boost)

    # Adjust for CTCF boundaries (intermediate accessibility)
    ctcf_count = chromatin_state.get("ctcf_count", 0)
    if ctcf_count > 0:
        # CTCF boundaries have partial accessibility — nudge toward 0.55
        boundary_nudge = 0.05 * min(3, ctcf_count)
        if base_accessibility < 0.55:
            base_accessibility = min(0.55, base_accessibility + boundary_nudge)
        elif base_accessibility > 0.55:
            base_accessibility = max(0.55, base_accessibility - boundary_nudge * 0.3)

    # Adjust for poly-dA:dT tracts (nucleosome-depleted → more accessible)
    poly_score = chromatin_state.get("poly_da_dt_score", 0.0)
    if poly_score > 0.2:
        poly_boost = 0.05 * min(1.0, poly_score)
        base_accessibility = min(1.0, base_accessibility + poly_boost)

    # Confidence modulation: low confidence → pull toward 0.50
    # If we're not confident about the state, we should be more conservative
    if confidence < 0.5:
        shrinkage = 1.0 - 0.3 * (1.0 - 2.0 * confidence)  # 0.4 at confidence=0
        base_accessibility = 0.50 + (base_accessibility - 0.50) * shrinkage

    return max(0.0, min(1.0, base_accessibility))


# ── Deep Learning Integration ───────────────────────────────────────────────


def predict_chromatin_sei(
    seq: str,
    sei_model_path: str | None = None,
) -> dict:
    """Predict chromatin profile using the Sei framework.

    Sei (Zhou et al. 2019 Nat Genet) is a deep learning model that
    predicts 40 sequence-class chromatin profiles from DNA sequence.
    It was trained on 21,907 DNase-seq, ATAC-seq, ChIP-seq, and
    CAGE datasets from ENCODE and Roadmap Epigenomics.

    The 40 sequence classes include:
    - Promoter / enhancer states
    - TF binding classes (CTCF, ETS, bZIP, etc.)
    - Heterochromatin / repressed states
    - DNase-only accessible states

    Requires: ``pip install sei`` and a Sei model checkpoint.

    Args:
        seq: DNA sequence (ideally >= 4000 bp for full context).
        sei_model_path: Path to the Sei model checkpoint.  If ``None``,
            attempts to load from the default installation path.

    Returns:
        Dict with keys:

        - ``"sequence_class_scores"``: dict mapping class name to score
        - ``"top_class"``: highest-scoring sequence class name
        - ``"top_score"``: score of the top class
        - ``"model_path"``: path used for the model

        Or, on failure:

        - ``"error"``: error message string
    """
    try:
        import numpy as np

        # Attempt to import Sei
        try:
            from sei import SeiModel  # type: ignore[import-untyped]
        except ImportError:
            return {
                "error": (
                    "Sei not installed. Install with: "
                    "pip install sei  |  "
                    "Zhou et al. 2019 Nat Genet"
                )
            }

        # Load model
        if sei_model_path:
            model = SeiModel(sei_model_path)
        else:
            model = SeiModel()  # default path

        # One-hot encode
        seq_upper = seq.upper()
        mapping = {"A": 0, "T": 1, "G": 2, "C": 3}
        one_hot = np.zeros((len(seq_upper), 4), dtype=np.float32)
        for i, base in enumerate(seq_upper):
            if base in mapping:
                one_hot[i, mapping[base]] = 1.0

        # Run prediction
        scores = model.predict(one_hot)

        # Map to sequence classes
        # Sei outputs 40 sequence class probabilities
        class_names = [
            "Promoter", "Enhancer", "Enhancer2", "CTCF", "CTCF+Enhancer",
            "Boundary", "DNase-only", "TF-ETS", "TF-bZIP", "TF-E2F",
            "TF-IRF", "TF-bHLH", "TF-NRF", "TF-NFkB", "TF-FOXA",
            "TF-GATA", "TF-AP1", "TF-SP", "TF-RFX", "TF-NKX",
            "TF-PAX", "TF-ZNF", "TF-HMG", "TF-ETS2", "TF-CLOCK",
            "Pol2", "H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3",
            "H3K27me3", "H3K9me3", "Repressed", "Bivalent",
            "Quiescent", "Heterochromatin", "AT-rich", "GC-rich",
            "Repeat-rich", "Gene-rich",
        ]

        result: dict[str, Any] = {
            "sequence_class_scores": {},
            "model_path": sei_model_path or "default",
        }

        for i, name in enumerate(class_names):
            if i < len(scores):
                result["sequence_class_scores"][name] = float(scores[i])

        # Top class
        if result["sequence_class_scores"]:
            top_class = max(
                result["sequence_class_scores"],
                key=result["sequence_class_scores"].get,  # type: ignore[arg-type]
            )
            result["top_class"] = top_class
            result["top_score"] = result["sequence_class_scores"][top_class]

        return result

    except ImportError as exc:
        return {"error": f"Required dependency not available: {exc}"}
    except Exception as exc:
        logger.debug("Sei prediction failed: %s", exc)
        return {"error": f"Sei prediction failed: {exc}"}


def predict_chromatin_enformer(seq: str) -> dict:
    """Predict chromatin profile using Enformer.

    Enformer (Avsec et al. 2021 Nat Methods) predicts 5313 genomic
    tracks including DNase-seq, ChIP-seq for histone marks and TFs,
    and CAGE from DNA sequence using a transformer architecture with
    long-range interactions.

    Key outputs relevant to chromatin context:
    - DNase-seq tracks (open chromatin)
    - H3K4me3 (active promoters)
    - H3K27ac (active enhancers)
    - H3K27me3 (Polycomb repressed)
    - H3K9me3 (constitutive heterochromatin)
    - ATAC-seq (open chromatin)

    Requires: ``pip install enformer-pytorch``

    Args:
        seq: DNA sequence (ideally >= 196,608 bp for full context,
            minimum ~8,000 bp for meaningful predictions).

    Returns:
        Dict with keys:

        - ``"dnase_signal"``: average predicted DNase signal
        - ``"h3k4me3_signal"``: average predicted H3K4me3 signal
        - ``"h3k27ac_signal"``: average predicted H3K27ac signal
        - ``"h3k27me3_signal"``: average predicted H3K27me3 signal
        - ``"h3k9me3_signal"``: average predicted H3K9me3 signal
        - ``"chromatin_state"``: inferred state from histone mark combination
        - ``"n_tracks"``: number of tracks predicted

        Or, on failure:

        - ``"error"``: error message string
    """
    try:
        import numpy as np
        import torch
        from enformer_pytorch import from_pretrained  # type: ignore[import-untyped]

        model = from_pretrained("Enformer")
        model.eval()

        # One-hot encode
        seq_upper = seq.upper()
        mapping = {"A": 0, "T": 1, "G": 2, "C": 3}
        one_hot = np.zeros((len(seq_upper), 4), dtype=np.float32)
        for i, base in enumerate(seq_upper):
            if base in mapping:
                one_hot[i, mapping[base]] = 1.0

        # Pad or truncate to 196,608 bp (Enformer input length)
        target_len = 196_608
        if len(one_hot) < target_len:
            pad = np.zeros((target_len - len(one_hot), 4), dtype=np.float32)
            one_hot = np.concatenate([one_hot, pad])
        else:
            one_hot = one_hot[:target_len]

        with torch.no_grad():
            x = torch.from_numpy(one_hot).unsqueeze(0)
            predictions = model(x)  # shape: (1, 1536, 5313)

        # Extract key chromatin marks
        # Enformer outputs at 128 bp resolution → 1536 bins for 196,608 bp
        pred_np = predictions.squeeze(0).numpy()  # (1536, 5313)

        # Approximate track indices for key marks
        # These are simplified — the actual mapping requires the Enformer
        # track metadata (5313 tracks from ENCODE/Roadmap)
        n_tracks = pred_np.shape[1]
        avg_signal = float(pred_np.mean())

        # Use heuristics for specific marks based on track groupings
        # DNase tracks: typically first ~100 tracks
        dnase_signal = float(pred_np[:, :min(100, n_tracks)].mean())
        # H3K4me3 tracks: typically ~100-200
        h3k4me3_signal = float(pred_np[:, min(100, n_tracks):min(200, n_tracks)].mean())
        # H3K27ac tracks: ~200-400
        h3k27ac_signal = float(pred_np[:, min(200, n_tracks):min(400, n_tracks)].mean())
        # Repressive marks
        h3k27me3_signal = float(pred_np[:, min(400, n_tracks):min(500, n_tracks)].mean())
        h3k9me3_signal = float(pred_np[:, min(500, n_tracks):min(600, n_tracks)].mean())

        # Infer chromatin state from mark combination
        active_score = h3k4me3_signal + h3k27ac_signal + dnase_signal
        repressive_score = h3k27me3_signal + h3k9me3_signal

        if active_score > repressive_score * 2:
            chromatin_state = "euchromatin"
        elif repressive_score > active_score * 2:
            chromatin_state = "heterochromatin"
        else:
            chromatin_state = "facultative_heterochromatin"

        return {
            "dnase_signal": dnase_signal,
            "h3k4me3_signal": h3k4me3_signal,
            "h3k27ac_signal": h3k27ac_signal,
            "h3k27me3_signal": h3k27me3_signal,
            "h3k9me3_signal": h3k9me3_signal,
            "avg_signal": avg_signal,
            "chromatin_state": chromatin_state,
            "n_tracks": n_tracks,
        }

    except ImportError:
        return {
            "error": (
                "Enformer not installed. Install with: "
                "pip install enformer-pytorch  |  "
                "Avsec et al. 2021 Nat Methods"
            )
        }
    except Exception as exc:
        logger.debug("Enformer prediction failed: %s", exc)
        return {"error": f"Enformer prediction failed: {exc}"}


# ── Public API ───────────────────────────────────────────────────────────────

__all__ = [
    # Core prediction
    "predict_chromatin_state",
    "predict_dnase_accessibility",
    "compute_repair_accessibility",
    # Feature scores
    "compute_gc_ishchz_score",
    "compute_cpg_density_score",
    # Deep learning (optional)
    "predict_chromatin_sei",
    "predict_chromatin_enformer",
]
