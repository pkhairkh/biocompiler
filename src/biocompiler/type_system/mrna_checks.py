"""
BioCompiler Type System — mRNA-Level Predicate Checks
=====================================================
mRNA-structure and translation-dynamics predicate checks: unexpected
transmembrane domains, mRNA secondary structure, co-translational folding,
mRNA stability, and cryptic ORFs.

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

from .sequence_checks import (
    _is_prokaryotic_organism,
    _translate_dna_to_aa,
    _resolve_species_cai,
    _compute_codon_ramp_score,
)


logger = logging.getLogger(__name__)

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
            from biocompiler.engines.viennarna import is_viennarna_available, compute_5prime_dg
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
            from biocompiler.engines.viennarna_fallback import compute_approx_dg
            dg = compute_approx_dg(seq, region="5utr")
            if dg <= effective_threshold:
                # UNCERTAIN because Nussinov fallback has ~60-70% sensitivity
                # and cannot predict pseudoknots; FAIL would overstate confidence.
                return PredicateResult(
                    "mRNASecondaryStructure", False, verdict=Verdict.UNCERTAIN,
                    details=(f"Strong mRNA secondary structure (Nussinov fallback): "
                             f"dG~{dg:.1f} kcal/mol <= {effective_threshold} "
                             f"[UNCERTAIN: heuristic estimate, not ViennaRNA]"),
                )
            if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
                return PredicateResult(
                    "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
                    details=(f"Moderate mRNA secondary structure (Nussinov fallback): "
                             f"dG~{dg:.1f} kcal/mol <= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f} "
                             f"[UNCERTAIN: heuristic estimate, not ViennaRNA]"),
                )
            # UNCERTAIN because Nussinov fallback is a heuristic, not a
            # calibrated thermodynamic model like ViennaRNA.
            return PredicateResult(
                "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
                details=(f"Weak mRNA secondary structure (Nussinov fallback): "
                         f"dG~{dg:.1f} kcal/mol "
                         f"[UNCERTAIN: heuristic estimate, not ViennaRNA]"),
            )
        except Exception:
            logger.debug("Falling through to legacy toy model", exc_info=True)

    # Legacy toy model (original implementation, backward-compatible)
    # WARNING: This is an extremely simplified model; all verdicts are UNCERTAIN
    # because the toy model does not compute real thermodynamic energies.
    seq = seq.upper()
    # Extract the window around the RBS/start codon
    effective_end = min(window_end, len(seq))
    window_seq = seq[window_start:effective_end]

    if len(window_seq) < 4:
        return PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
            details="Sequence window too short for structure analysis "
                    "[UNCERTAIN: legacy toy model, not ViennaRNA]",
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
        # UNCERTAIN because the legacy toy model is an uncalibrated heuristic;
        # FAIL would overstate confidence in this rough dG estimate.
        return PredicateResult(
            "mRNASecondaryStructure", False, verdict=Verdict.UNCERTAIN,
            details=(f"Strong mRNA secondary structure: dG={dg:.1f} kcal/mol "
                     f"<= {effective_threshold} (GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs) "
                     f"[UNCERTAIN: legacy toy model, not ViennaRNA]"),
        )

    if dg <= effective_threshold * _MRNA_MODERATE_DG_RATIO:
        window_gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) if window_seq else 0.5
        # All verdicts from the legacy toy model are UNCERTAIN because it
        # uses a mirror-pairing heuristic, not real thermodynamic parameters.
        if window_gc < 0.5:
            return PredicateResult(
                "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
                details=(f"Moderate mRNA secondary structure: dG={dg:.1f} kcal/mol, "
                         f"but AT-rich window (GC={window_gc:.0%}) weakens structure "
                         f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs) "
                         f"[UNCERTAIN: legacy toy model, not ViennaRNA]"),
            )
        return PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
            details=(f"Moderate mRNA secondary structure: dG={dg:.1f} kcal/mol "
                     f"<= {effective_threshold * _MRNA_MODERATE_DG_RATIO:.1f} "
                     f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs) "
                     f"[UNCERTAIN: legacy toy model, not ViennaRNA]"),
        )

    # UNCERTAIN because the legacy toy model is an uncalibrated heuristic;
    # PASS would overstate confidence in this rough dG estimate.
    return PredicateResult(
        "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
        details=(f"Weak mRNA secondary structure: dG={dg:.1f} kcal/mol "
                 f"(GC={gc_pairs}, AU={au_pairs}, GU={gu_pairs} pairs) "
                 f"[UNCERTAIN: legacy toy model, not ViennaRNA]"),
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
    from biocompiler.expression.mrna_stability import score_mrna_stability

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

def check_no_cryptic_orf(
    seq: str,
    min_orf_length: int = 30,
    organism: str = "",
) -> PredicateResult:
    """Predicate: No cryptic open reading frames.

    Scans for out-of-frame ATG codons that form open reading frames
    (ATG to in-frame stop codon) in alternative reading frames (+1 and +2).
    Cryptic ORFs can compete for ribosomes (leaky scanning), produce
    immunogenic peptides via non-canonical translation, or trigger
    nonsense-mediated decay if the ORF is long enough.

    For prokaryotes, also scans for Shine-Dalgarno-like sequences (AGGAGG)
    upstream of out-of-frame start codons, which increase translation
    initiation likelihood.

    A cryptic ORF is flagged when:
    - Out-of-frame ATG followed by in-frame stop > min_orf_length codons away
    - FAIL: cryptic ORF >= 90 codons (270 nt) — strong translation risk
    - UNCERTAIN: cryptic ORF >= min_orf_length but < 90 codons
    - PASS: no cryptic ORFs above min_orf_length

    Args:
        seq: DNA coding sequence (uppercase).
        min_orf_length: Minimum ORF length (in codons) to flag. Default 30.
        organism: Organism name for prokaryotic SD-sequence detection.

    Returns:
        PredicateResult with verdict and positions of cryptic ORFs.
    """
    seq = seq.upper()
    stop_codons = {"TAA", "TAG", "TGA"}
    is_prokaryote = _is_prokaryotic_organism(organism) if organism else False

    cryptic_orfs: List[Dict[str, Any]] = []

    for frame in (1, 2):
        # Scan for ATG codons in this alternative frame
        pos = frame
        while pos + 3 <= len(seq):
            codon = seq[pos:pos + 3]
            if codon == "ATG":
                # Found out-of-frame ATG — look for in-frame stop codon
                orf_start = pos
                stop_pos = None
                scan_pos = pos + 3
                codon_count = 1  # count the ATG as the first codon

                while scan_pos + 3 <= len(seq):
                    next_codon = seq[scan_pos:scan_pos + 3]
                    if next_codon in stop_codons:
                        stop_pos = scan_pos
                        break
                    codon_count += 1
                    scan_pos += 3

                if codon_count >= min_orf_length:
                    # Check for Shine-Dalgarno in prokaryotes
                    has_sd = False
                    if is_prokaryote:
                        # Look for AGGAGG in a wider upstream window (up to
                        # 20 bp before the ATG), then verify the 3' end of
                        # the match is within 4-12 bp of the start codon.
                        sd_search_start = max(0, orf_start - 20)
                        upstream = seq[sd_search_start:orf_start]
                        sd_pos = upstream.find("AGGAGG")
                        while sd_pos >= 0:
                            # Distance from 3' end of SD to ATG start
                            sd_end_in_seq = sd_search_start + sd_pos + 6
                            distance = orf_start - sd_end_in_seq
                            if 4 <= distance <= 12:
                                has_sd = True
                                break
                            sd_pos = upstream.find("AGGAGG", sd_pos + 1)

                    cryptic_orfs.append({
                        "frame": frame,
                        "start": orf_start,
                        "stop": stop_pos,
                        "length_codons": codon_count,
                        "has_sd": has_sd,
                    })
            pos += 3

    if not cryptic_orfs:
        return PredicateResult(
            "NoCrypticORF", True, verdict=Verdict.PASS,
            details="No cryptic ORFs above minimum length threshold",
        )

    # Determine worst verdict
    worst_verdict = Verdict.PASS
    fail_orfs = []
    uncertain_orfs = []

    for orf in cryptic_orfs:
        length = orf["length_codons"]
        if length >= 90:
            # Long cryptic ORF — strong translation risk
            # Prokaryote with SD is even worse
            if orf["has_sd"]:
                worst_verdict = Verdict.FAIL
            else:
                worst_verdict = Verdict.FAIL if worst_verdict != Verdict.FAIL else Verdict.FAIL
            fail_orfs.append(orf)
        elif length >= min_orf_length:
            # Medium cryptic ORF — uncertain risk
            if orf["has_sd"]:
                # SD sequence increases risk — upgrade to FAIL
                fail_orfs.append(orf)
                worst_verdict = Verdict.FAIL
            else:
                uncertain_orfs.append(orf)
                if worst_verdict == Verdict.PASS:
                    worst_verdict = Verdict.UNCERTAIN

    positions = [orf["start"] for orf in cryptic_orfs]
    passed = worst_verdict != Verdict.FAIL

    # Build details string
    detail_parts = []
    if fail_orfs:
        detail_parts.append(
            f"FAIL-level cryptic ORFs: {[(o['start'], o['length_codons'], 'frame+' + str(o['frame'])) for o in fail_orfs]}"
        )
    if uncertain_orfs:
        detail_parts.append(
            f"UNCERTAIN-level cryptic ORFs: {[(o['start'], o['length_codons'], 'frame+' + str(o['frame'])) for o in uncertain_orfs]}"
        )
    details = "; ".join(detail_parts) if detail_parts else "No cryptic ORFs above minimum length threshold"

    return PredicateResult(
        "NoCrypticORF", passed, verdict=worst_verdict,
        details=details,
        positions=positions,
    )




__all__ = [
    "check_no_unexpected_tm_domain",
    "check_mrna_secondary_structure",
    "check_co_translational_folding",
    "check_mrna_stability",
    "check_no_cryptic_orf",
]
