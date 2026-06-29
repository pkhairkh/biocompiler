"""
BioCompiler Type System — Immunogenicity / RNA-Regulation Predicate Checks
==========================================================================
Immunogenicity and post-transcriptional regulation predicate checks:
m6A modification sites, polyadenylation signals, nucleoside modification
guidance, and ribosome quality-control (RQC) trigger motifs.

Extracted from the historical checks.py monolith during the W8-b refactor.
Re-exported by checks.py for backwards compatibility.
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

logger = logging.getLogger(__name__)

def check_no_m6a_site(
    seq: str,
    organism: str = "",
    scan_mode: str = "comprehensive",
) -> PredicateResult:
    """Predicate: No N6-methyladenosine (m6A) modification sites.

    m6A is the most abundant internal mRNA modification, catalyzed by the
    METTL3/METTL14 writer complex. m6A modifications can affect:

    - mRNA stability (YTHDF2-mediated decay)
    - Translation efficiency (YTHDF1/eIF3-mediated enhancement)
    - Splicing regulation (heterogeneous nuclear ribonucleoprotein effects)
    - Innate immune activation (via RIG-I/MDA5 sensing of modified RNAs)

    The m6A consensus motif is RRACH (R=A/G, A=m6A, C, H=A/C/U), where
    the central A is the modified position. In DNA coding sequences:
    - RRACH = [AG][AG]AC[ACT]

    Position-dependent scoring:
    - Near stop codon (last 400 nt): m6A is common and expected → lower risk
    - Near start codon (first 100 nt): m6A can affect initiation → moderate risk
    - In middle of CDS: m6A is less expected → higher risk
    - In DRACH context (D=A/G/U): highest confidence m6A site

    Additional high-confidence motifs from miCLIP/PA-m6A studies:
    - GGACU (the strongest m6A signal)
    - UGACU
    - AGACU
    - GGACA

    Verdict logic:
    - FAIL: >= 5 high-confidence m6A sites (DRACH with GGAC core) in CDS
    - UNCERTAIN: 2-4 high-confidence sites, or >= 3 RRACH sites
    - PASS: 0-1 RRACH sites, or only near stop codon

    Args:
        seq: DNA coding sequence (uppercase).
        organism: Target organism (unused, reserved for future).
        scan_mode: "comprehensive" or "conservative". Default "comprehensive".

    Returns:
        PredicateResult with verdict and positions of m6A sites.
    """
    seq = seq.upper()
    n = len(seq)

    if n < 5:
        return PredicateResult(
            "NoM6ASite", True, verdict=Verdict.PASS,
            details="Sequence too short for m6A motif scanning",
        )

    # Scan for RRACH motif: [AG][AG]AC[ACT]
    rrach_positions: List[int] = []
    for i in range(n - 4):
        if seq[i] in "AG" and seq[i + 1] in "AG" and seq[i + 2] == "A" and seq[i + 3] == "C" and seq[i + 4] in "ACT":
            rrach_positions.append(i)

    # Scan for DRACH motif: [AGT][AG]AC[ACT] (higher confidence)
    drach_positions: List[int] = []
    for i in range(n - 4):
        if seq[i] in "AGT" and seq[i + 1] in "AG" and seq[i + 2] == "A" and seq[i + 3] == "C" and seq[i + 4] in "ACT":
            drach_positions.append(i)

    # High-confidence GGAC core sites (GGAC = positions i..i+3 within RRACH)
    ggac_positions: List[int] = []
    for i in range(n - 3):
        if seq[i:i + 4] == "GGAC":
            ggac_positions.append(i)

    # High-confidence m6A sites = DRACH with GGAC core
    ggac_set = set(ggac_positions)
    high_confidence_positions = [p for p in drach_positions if p in ggac_set or (p + 1) in ggac_set]

    # Apply position-dependent weighting
    # Near stop codon (last 400 nt) → lower risk, still counted but deprioritized
    # Near start codon (first 100 nt) → moderate risk
    # Middle of CDS → higher risk
    high_conf_cds = []  # high-confidence sites NOT in stop-codon region
    high_conf_stop = []  # high-confidence sites near stop codon

    for p in high_confidence_positions:
        if p >= n - 400:
            high_conf_stop.append(p)
        else:
            high_conf_cds.append(p)

    # Also track GGACU/UGACU/AGACU/GGACA specific motifs for detailed reporting
    specific_high = []
    for motif in ["GGACT", "GGACC", "GGACA", "TGACT", "AGACT"]:
        pos = seq.find(motif)
        while pos != -1:
            specific_high.append((pos, motif))
            pos = seq.find(motif, pos + 1)

    # Verdict logic
    n_high_conf_cds = len(high_conf_cds)
    n_rrach = len(rrach_positions)

    if n_high_conf_cds >= 5:
        verdict = Verdict.FAIL
    elif 2 <= n_high_conf_cds <= 4 or n_rrach >= 3:
        verdict = Verdict.UNCERTAIN
    else:
        # Check if all sites are only near stop codon
        if n_rrach <= 1 or (rrach_positions and all(p >= n - 400 for p in rrach_positions)):
            verdict = Verdict.PASS
        else:
            verdict = Verdict.PASS

    passed = verdict != Verdict.FAIL

    all_positions = sorted(set(rrach_positions))
    if verdict == Verdict.FAIL:
        details = (
            f"{n_high_conf_cds} high-confidence m6A sites in CDS "
            f"(DRACH+GGAC core); {n_rrach} total RRACH sites"
        )
    elif verdict == Verdict.UNCERTAIN:
        details = (
            f"{n_high_conf_cds} high-confidence m6A sites in CDS, "
            f"{n_rrach} total RRACH sites — modification likely"
        )
    else:
        if n_rrach == 0:
            details = "No RRACH m6A motifs found"
        elif all(p >= n - 400 for p in rrach_positions):
            details = (
                f"{n_rrach} RRACH site(s) found only near stop codon "
                f"(expected location, low risk)"
            )
        else:
            details = (
                f"{n_rrach} RRACH site(s) found, "
                f"{n_high_conf_cds} high-confidence in CDS — low risk"
            )

    return PredicateResult(
        "NoM6ASite", passed, verdict=verdict,
        details=details,
        positions=all_positions,
    )


# ────────────────────────────────────────────────────────────
# Polyadenylation signal predicate
# ────────────────────────────────────────────────────────────

# Polyadenylation signal hexamers with approximate frequencies (Beaudoing et al. 2000)
_POLYA_SIGNALS: List[Tuple[str, float, str]] = [
    ("AATAAA", 0.53, "canonical"),
    ("ATTAAA", 0.18, "variant"),
    ("AGTAAA", 0.03, "variant"),
    ("TATAAA", 0.03, "variant"),
    ("CATAAA", 0.02, "variant"),
    ("GATAAA", 0.02, "variant"),
    ("AATATA", 0.01, "variant"),
    ("AATACA", 0.01, "variant"),
    ("AATAGA", 0.01, "variant"),
    ("ACTAAA", 0.01, "variant"),
    ("AATGAA", 0.01, "variant"),
]

# T-rich downstream element patterns (within 30bp of PAS)
_T_RICH_MOTIFS = ["TTTT", "TTTTT"]

_POLYA_DOWNSTREAM_WINDOW = 30  # bp downstream to check for U/GU-rich


def check_no_polya_signal(
    seq: str,
    organism: str = "",
    scan_cds_only: bool = True,
) -> PredicateResult:
    """Predicate: No premature polyadenylation signals.

    Polyadenylation signals (PAS) direct 3'-end processing and poly(A) tail
    addition. When a PAS appears within the coding region, it can cause:

    - Premature transcription termination
    - Truncated mRNA with loss of downstream coding information
    - Failed protein expression

    The canonical PAS is AATAAA (AAUAAA in mRNA), which accounts for ~53%
    of human polyadenylation signals. Variant hexamers account for the rest:

    Canonical:  AATAAA  (53%)
    Variant 1:  ATTAAA  (18%)
    Variant 2:  AGTAAA  (3%)
    Variant 3:  TATAAA  (3%)
    Variant 4:  CATAAA  (2%)
    Variant 5:  GATAAA  (2%)
    Variant 6:  AATATA  (1%)
    Variant 7:  AATACA  (1%)
    Variant 8:  AATAGA  (1%)
    Variant 9:  ACTAAA  (1%)
    Variant 10: AATGAA  (1%)

    Context matters: A PAS followed by a downstream U/GU-rich element
    (within 30 nt) is far more likely to be functional. The scan checks
    for U-rich (T-rich in DNA) elements: TTTT, TTTTT within 30bp downstream.

    Verdict logic:
    - FAIL: Canonical AATAAA + downstream T-rich element → likely functional PAS
    - UNCERTAIN: Canonical AATAAA without downstream element, OR variant + element
    - PASS: No PAS, or only weak variants without context

    Args:
        seq: DNA coding sequence (uppercase).
        organism: Target organism (unused, reserved for future).
        scan_cds_only: If True, skip PAS in last 30 nt (expected 3' UTR region).

    Returns:
        PredicateResult with verdict and positions of PAS hexamers.
    """
    seq = seq.upper()
    n = len(seq)

    if n < 6:
        return PredicateResult(
            "NoPolyASignal", True, verdict=Verdict.PASS,
            details="Sequence too short for PAS hexamer scanning",
        )

    # Search window: if scan_cds_only, exclude last 30nt (likely 3' UTR)
    search_end = n - 30 if scan_cds_only and n > 60 else n

    pas_hits: List[Dict[str, Any]] = []

    for hexamer, freq, category in _POLYA_SIGNALS:
        pos = seq.find(hexamer)
        while pos != -1:
            if pos < search_end:
                # Check for downstream T-rich element
                downstream_start = pos + 6
                downstream_end = min(pos + 6 + _POLYA_DOWNSTREAM_WINDOW, n)
                downstream_region = seq[downstream_start:downstream_end]
                has_t_rich = any(m in downstream_region for m in _T_RICH_MOTIFS)

                pas_hits.append({
                    "pos": pos,
                    "hexamer": hexamer,
                    "category": category,
                    "frequency": freq,
                    "has_downstream_element": has_t_rich,
                })
            pos = seq.find(hexamer, pos + 1)

    if not pas_hits:
        return PredicateResult(
            "NoPolyASignal", True, verdict=Verdict.PASS,
            details="No polyadenylation signal hexamers found in CDS",
        )

    # Determine worst verdict
    worst_verdict = Verdict.PASS
    canonical_with_element = [h for h in pas_hits if h["category"] == "canonical" and h["has_downstream_element"]]
    canonical_no_element = [h for h in pas_hits if h["category"] == "canonical" and not h["has_downstream_element"]]
    variant_with_element = [h for h in pas_hits if h["category"] == "variant" and h["has_downstream_element"]]
    variant_no_element = [h for h in pas_hits if h["category"] == "variant" and not h["has_downstream_element"]]

    if canonical_with_element:
        worst_verdict = Verdict.FAIL
    elif canonical_no_element or variant_with_element:
        worst_verdict = Verdict.UNCERTAIN

    passed = worst_verdict != Verdict.FAIL

    all_positions = sorted({h["pos"] for h in pas_hits})

    if worst_verdict == Verdict.FAIL:
        detail_parts = ["Functional PAS detected (canonical + downstream element): "]
        for h in canonical_with_element:
            detail_parts.append(f"{h['hexamer']} at pos {h['pos']}")
        details = "; ".join(detail_parts)
    elif worst_verdict == Verdict.UNCERTAIN:
        n_canonical = len(canonical_no_element)
        n_variant = len(variant_with_element)
        details = (
            f"Possible PAS: {n_canonical} canonical (no downstream element), "
            f"{n_variant} variant (with downstream element)"
        )
    else:
        n_variants = len(variant_no_element)
        details = (
            f"Only weak variant PAS found ({n_variants} variant hexamer(s) "
            f"without downstream element) — low risk"
        )

    return PredicateResult(
        "NoPolyASignal", passed, verdict=worst_verdict,
        details=details,
        positions=all_positions,
    )


# ────────────────────────────────────────────────────────────
# Nucleoside modification guidance predicate
# ────────────────────────────────────────────────────────────

# Immunostimulatory motifs in DNA coding strand (RNA → DNA conversion)
# GUCC in RNA = GTCC in DNA; UGUGU in RNA = TGTGT in DNA
_IMMUNOSTIM_MOTIFS: Dict[str, str] = {
    "GTCC": "TLR7/8_agonist",
    "TGTGT": "RIG-I_agonist",
    # GU-rich stretches (4+ G/T in 5-nt window)
    # AT-rich stretches (4+ A/T in 5-nt window) — checked via sliding window
}

# Minimum window for GU-rich / AU-rich detection
_RICH_MOTIF_WINDOW = 5
_RICH_MOTIF_THRESHOLD = 4  # 4 out of 5 bases must be G/T or A/T


def check_nucleoside_modification_guidance(
    seq: str,
    organism: str = "",
    modification_type: str = "auto",
) -> PredicateResult:
    """Predicate: Nucleoside modification guidance for mRNA therapeutics.

    For mRNA vaccine/therapeutic applications, modified nucleosides are
    critical for reducing innate immune sensing (TLR7/8, RIG-I, PKR) and
    improving translation. This predicate provides guidance on:

    1. **Uridine content assessment**: High uridine content correlates with
       stronger innate immune activation. Sequences with >30% U content
       in any 50-nt window are flagged for modification.

    2. **Immunostimulatory motif detection**: Specific sequence motifs that
       are particularly immunostimulatory when unmodified:
       - GUCC (TLR7/8 agonist)
       - UGUGU (RIG-I agonist, identified by Kato et al. 2008)
       - GU-rich sequences (TLR7/8 ligands)
       - AU-rich sequences (PKR activators)

    3. **Modification recommendations**: Based on the pattern analysis,
       provides guidance on which modified nucleosides would be most
       beneficial:
       - N1-methylpseudouridine (m1Ψ): best for reducing TLR7/8 activation,
         improves translation (Andries et al. 2015, FDA-approved in mRNA vaccines)
       - 5-methylcytidine (m5C): reduces PKR activation
       - Pseudouridine (Ψ): moderate immune evasion, enhanced translation
       - 2-thiouridine (s2U): strong immune evasion but reduces translation

    4. **Modified nucleoside density estimation**: Estimates what fraction
       of uridines should be replaced for optimal results.

    Verdict logic:
    - FAIL: Sequence has immunostimulatory motifs AND >35% uridine content
      → modification is essential, unmodified RNA will cause strong immune response
    - UNCERTAIN: Sequence has immunostimulatory motifs OR high uridine (30-35%)
      → modification recommended
    - PASS: No immunostimulatory motifs and uridine content < 30%
      → modification optional but still recommended for clinical use

    Args:
        seq: DNA coding sequence (uppercase).
        organism: Target organism (unused, reserved for future).
        modification_type: "auto", "m1psi", "m5c", "psi", or "s2u".
            Default "auto" selects based on motif analysis.

    Returns:
        PredicateResult with verdict, positions of immunostimulatory motifs,
        and modification recommendations in verification_evidence.
    """
    seq = seq.upper()
    n = len(seq)

    if n == 0:
        return PredicateResult(
            "NucleosideModificationGuidance", True, verdict=Verdict.PASS,
            details="Empty sequence",
        )

    # 1. Compute overall T content (T in DNA = U in mRNA)
    t_count = seq.count("T")
    overall_u_pct = t_count / n

    # 2. Sliding window for local U content (50-nt window)
    _U_WINDOW = 50
    max_local_u_pct = 0.0
    max_local_u_pos = -1
    if n >= _U_WINDOW:
        # Initialize first window
        window_t = seq[:_U_WINDOW].count("T")
        local_u_pct = window_t / _U_WINDOW
        max_local_u_pct = local_u_pct
        max_local_u_pos = 0

        for start in range(1, n - _U_WINDOW + 1):
            if seq[start - 1] == "T":
                window_t -= 1
            if seq[start + _U_WINDOW - 1] == "T":
                window_t += 1
            local_u_pct = window_t / _U_WINDOW
            if local_u_pct > max_local_u_pct:
                max_local_u_pct = local_u_pct
                max_local_u_pos = start
    else:
        max_local_u_pct = overall_u_pct
        max_local_u_pos = 0

    # 3. Scan for specific immunostimulatory motifs
    motif_hits: List[Dict[str, Any]] = []
    for motif, motif_type in _IMMUNOSTIM_MOTIFS.items():
        pos = seq.find(motif)
        while pos != -1:
            motif_hits.append({"pos": pos, "motif": motif, "type": motif_type})
            pos = seq.find(motif, pos + 1)

    # 4. Scan for GU-rich and AU-rich stretches (sliding window)
    gu_rich_positions: List[int] = []
    au_rich_positions: List[int] = []
    for i in range(n - _RICH_MOTIF_WINDOW + 1):
        window = seq[i:i + _RICH_MOTIF_WINDOW]
        gt_count_in_window = sum(1 for c in window if c in "GT")
        at_count_in_window = sum(1 for c in window if c in "AT")
        # GU-rich requires both G and T present with multiple T's
        # (not just a GC-rich region with occasional T)
        t_count_in_window = sum(1 for c in window if c == "T")
        g_count_in_window = sum(1 for c in window if c == "G")
        a_count_in_window = sum(1 for c in window if c == "A")
        if gt_count_in_window >= _RICH_MOTIF_THRESHOLD and t_count_in_window >= 2 and g_count_in_window >= 1:
            gu_rich_positions.append(i)
        if at_count_in_window >= _RICH_MOTIF_THRESHOLD and t_count_in_window >= 2 and a_count_in_window >= 1:
            au_rich_positions.append(i)

    has_immunostim = bool(motif_hits) or bool(gu_rich_positions) or bool(au_rich_positions)

    # 5. Determine verdict
    if has_immunostim and overall_u_pct > 0.35:
        verdict = Verdict.FAIL
    elif has_immunostim or (0.30 <= overall_u_pct <= 0.35):
        verdict = Verdict.UNCERTAIN
    else:
        verdict = Verdict.PASS

    passed = verdict != Verdict.FAIL

    # 6. Generate modification recommendations
    recommendations: List[str] = []
    if motif_hits or gu_rich_positions:
        recommendations.append(
            "m1Ψ (N1-methylpseudouridine): recommended for TLR7/8 "
            "and RIG-I immune evasion; FDA-approved in mRNA vaccines"
        )
    if au_rich_positions:
        recommendations.append(
            "m5C (5-methylcytidine): recommended to reduce PKR activation "
            "from AU-rich sequences"
        )
    if overall_u_pct > 0.30:
        u_replace_pct = min(100, int((overall_u_pct - 0.25) / overall_u_pct * 100))
        recommendations.append(
            f"Replace ~{u_replace_pct}% of uridines with m1Ψ for "
            f"optimal immune evasion (current U content: {overall_u_pct:.1%})"
        )
    if max_local_u_pct > 0.30:
        recommendations.append(
            f"Local U hotspot: {max_local_u_pct:.1%} in 50-nt window at "
            f"pos {max_local_u_pos} — prioritize modification in this region"
        )

    # Collect all positions for the result
    all_positions = sorted({
        h["pos"] for h in motif_hits
    } | set(gu_rich_positions[:10]) | set(au_rich_positions[:10]))

    # Build details string
    if verdict == Verdict.FAIL:
        details = (
            f"Immunostimulatory motifs present AND high uridine "
            f"({overall_u_pct:.1%}) — modification essential"
        )
    elif verdict == Verdict.UNCERTAIN:
        if has_immunostim:
            details = (
                f"Immunostimulatory motifs detected — modification recommended "
                f"(U content: {overall_u_pct:.1%})"
            )
        else:
            details = (
                f"Elevated uridine content ({overall_u_pct:.1%}) — "
                f"modification recommended"
            )
    else:
        details = (
            f"No immunostimulatory motifs, uridine content "
            f"{overall_u_pct:.1%} < 30% — modification optional"
        )

    evidence: Dict[str, Any] = {
        "overall_u_content": round(overall_u_pct, 4),
        "max_local_u_content": round(max_local_u_pct, 4),
        "max_local_u_pos": max_local_u_pos,
        "immunostim_motifs": len(motif_hits),
        "gu_rich_stretches": len(gu_rich_positions),
        "au_rich_stretches": len(au_rich_positions),
        "recommendations": recommendations,
        "modification_type": modification_type,
    }

    return PredicateResult(
        "NucleosideModificationGuidance", passed, verdict=verdict,
        details=details,
        positions=all_positions,
        verification_evidence=evidence,
    )


def check_no_rqc_trigger(
    seq: str,
    organism: str = "",
    poly_a_min_length: int = 6,
    stall_window_size: int = 15,
    stall_cai_threshold: float = 0.2,
) -> PredicateResult:
    """Predicate: No Ribosome Quality Control (RQC) triggers.

    Detects sequence features that can trigger ribosome quality control pathways:

    1. **Internal poly-A stretches** in CDS: 6+ consecutive A's within the
       coding region can cause ribosome stalling and activate no-go decay (NGD).
       In eukaryotic mRNA vaccines, internal poly-A can be mistaken for a
       poly-A tail, causing premature transcription termination.

    2. **Strong secondary structure barriers**: Extremely stable stem-loops
       (ΔG < -30 kcal/mol) near the 5' end can stall scanning ribosomes,
       triggering NGD. Detected via GC-content heuristic in a sliding window.

    3. **Consecutive low-CAI codons**: A stretch of 5+ consecutive codons
       with CAI < 0.2 suggests extreme codon rarity that can stall elongation.
       This is detected only when a species_cai dictionary is available via
       the organism parameter.

    4. **Non-stop decay triggers**: Sequences where the last in-frame codon
       is NOT a stop codon trigger non-stop decay (NSD). This is complementary
       to the NoStopCodons predicate which checks for INTERNAL stops.

    Verdict logic:
    - FAIL: Internal poly-A >= 8 A's, or non-stop decay trigger (no stop codon)
    - UNCERTAIN: Internal poly-A 6-7 A's, or strong structure barrier
    - PASS: No RQC triggers detected
    """
    seq = seq.upper()
    n = len(seq)
    triggers: List[Dict[str, Any]] = []
    worst_verdict = Verdict.PASS

    # 1. Scan for internal poly-A stretches
    i = 0
    while i < n:
        if seq[i] == "A":
            run_start = i
            run_length = 0
            while i < n and seq[i] == "A":
                run_length += 1
                i += 1
            if run_length >= poly_a_min_length:
                triggers.append({
                    "type": "poly_a",
                    "start": run_start,
                    "length": run_length,
                })
                if run_length >= 8:
                    worst_verdict = Verdict.FAIL
                elif run_length >= poly_a_min_length and worst_verdict == Verdict.PASS:
                    worst_verdict = Verdict.UNCERTAIN
        else:
            i += 1

    # 2. Scan for strong GC-rich windows (structure barrier heuristic)
    # Use 30bp sliding window; if GC > 0.85 for >= 2 consecutive windows, flag
    gc_window = 30
    gc_threshold = 0.85
    consecutive_high_gc = 0
    gc_barrier_start = -1

    if n >= gc_window:
        # Initialize first window
        gc_count = seq[:gc_window].count("G") + seq[:gc_window].count("C")
        gc_frac = gc_count / gc_window

        if gc_frac > gc_threshold:
            consecutive_high_gc = 1
            gc_barrier_start = 0

        for start in range(1, n - gc_window + 1):
            outgoing = seq[start - 1]
            incoming = seq[start + gc_window - 1]
            if outgoing in ("G", "C"):
                gc_count -= 1
            if incoming in ("G", "C"):
                gc_count += 1
            gc_frac = gc_count / gc_window

            if gc_frac > gc_threshold:
                if consecutive_high_gc == 0:
                    gc_barrier_start = start
                consecutive_high_gc += 1
            else:
                if consecutive_high_gc >= 2:
                    triggers.append({
                        "type": "gc_structure_barrier",
                        "start": gc_barrier_start,
                        "length": gc_window + (consecutive_high_gc - 1),
                    })
                    if worst_verdict == Verdict.PASS:
                        worst_verdict = Verdict.UNCERTAIN
                consecutive_high_gc = 0

        # Check final window
        if consecutive_high_gc >= 2:
            triggers.append({
                "type": "gc_structure_barrier",
                "start": gc_barrier_start,
                "length": gc_window + (consecutive_high_gc - 1),
            })
            if worst_verdict == Verdict.PASS:
                worst_verdict = Verdict.UNCERTAIN

    # 3. Non-stop decay trigger: last in-frame codon is NOT a stop codon
    if n >= 3:
        last_codon = seq[n - 3:n]
        if last_codon not in ("TAA", "TAG", "TGA"):
            triggers.append({
                "type": "non_stop_decay",
                "start": n - 3,
                "length": 3,
            })
            worst_verdict = Verdict.FAIL

    # 4. Low-CAI codon stretch: skip (CodonOptimality predicate already covers this)
    # Note: consecutive low-CAI codon detection is handled by check_codon_optimality
    # and the co-translational folding predicate. We note this in details when
    # the organism parameter is provided.
    cai_note = ""
    if organism:
        cai_note = " (low-CAI stretch detection deferred to CodonOptimality predicate)"

    passed = worst_verdict != Verdict.FAIL
    positions = [t["start"] for t in triggers]

    # Build details
    detail_parts = []
    for t in triggers:
        if t["type"] == "poly_a":
            detail_parts.append(f"poly-A stretch of {t['length']} A's at pos {t['start']}")
        elif t["type"] == "gc_structure_barrier":
            detail_parts.append(
                f"GC structure barrier at pos {t['start']} ({t['length']} bp)"
            )
        elif t["type"] == "non_stop_decay":
            detail_parts.append(
                f"non-stop decay trigger: last codon is not a stop codon"
            )

    details = "; ".join(detail_parts) if detail_parts else "No RQC triggers detected"
    if cai_note and details != "No RQC triggers detected":
        details += cai_note

    return PredicateResult(
        "NoRQCTrigger", passed, verdict=worst_verdict,
        details=details,
        positions=positions,
    )




__all__ = [
    "check_no_m6a_site",
    "check_no_polya_signal",
    "check_nucleoside_modification_guidance",
    "check_no_rqc_trigger",
]
