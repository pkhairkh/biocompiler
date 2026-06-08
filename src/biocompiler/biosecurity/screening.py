"""
Main biosecurity screening API.

Contains the core screening function, the optimization gate hook,
risk classification, recommendation builder, and mode detection.
"""

from __future__ import annotations

import logging
import os
import warnings
from typing import Optional

from .fuzzy_matching import (
    _fuzzy_match_edit_distance,
    _fuzzy_match_hamming,
    reverse_complement,
)
from .hazard_signatures import (
    _DNA_SIGNATURES,
    _HAZARD_SIGNATURES,
    _PROTEIN_SIGNATURES,
)
from .kmer_similarity import (
    _RISK_ORDER,
    _compute_kmer_similarity,
    _extract_kmers,
    _max_risk,
)
from .pathogen_signatures import _MOTIF_TO_PATHOGEN, _PATHOGEN_SIGNATURES
from .types import (
    BiosecurityMode,
    BiosecurityReport,
    BiosecurityScreeningResult,
    HazardMatch,
    RiskLevel,
)

from ..exceptions import BiosecurityError

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Biosecurity mode
# ─────────────────────────────────────────────────────────────────────────────

def get_biosecurity_mode() -> BiosecurityMode:
    """Read the biosecurity screening mode from the environment.

    Returns one of ``"enforce"``, ``"warn"``, or ``"off"`` based on the
    ``BIOCOMPILER_BIOSECURITY_MODE`` environment variable.  Falls back
    to ``"enforce"`` when the variable is unset or set to an
    unrecognised value.  The value is case-insensitive.

    Returns
    -------
    BiosecurityMode
        The effective biosecurity mode.
    """
    raw = os.environ.get("BIOCOMPILER_BIOSECURITY_MODE", "").strip().lower()
    if raw in ("enforce", "warn", "off"):
        return raw  # type: ignore[return-value]
    return "enforce"


# ─────────────────────────────────────────────────────────────────────────────
# Core screening function
# ─────────────────────────────────────────────────────────────────────────────

def screen_hazardous_sequence(protein: str, dna: str = "") -> BiosecurityReport:
    """Screen a protein (and optional DNA) sequence against known hazardous
    signatures.

    Parameters
    ----------
    protein : str
        Amino acid sequence in single-letter code (uppercase).
    dna : str, optional
        Nucleotide sequence to screen for resistance marker signatures.

    Returns
    -------
    BiosecurityReport
        Screening result with risk level, matches, and recommendations.

    Notes
    -----
    - Protein screening uses sliding-window substring matching against
      short peptide motifs (8-12 aa).
    - Fuzzy matching via Hamming distance (1-2 substitutions) and
      Levenshtein edit distance (1 edit) is applied for motifs < 15 aa.
    - DNA screening uses substring matching against nucleotide patterns
      (15-21 nt) for antibiotic resistance markers, plus reverse
      complement screening.
    - Confidence scoring accounts for motif length (longer = higher
      confidence) and exact match position context.
    """
    protein = protein.upper().strip()
    dna = dna.upper().strip()

    matches: list[HazardMatch] = []
    flagged_categories: set[str] = set()
    # Track exact match positions per signature name for deduplication
    exact_positions: dict[str, set[int]] = {}

    # ── Screen protein against peptide motifs ────────────────────────────
    for sig in _PROTEIN_SIGNATURES:
        motif = sig["motif"].upper()
        motif_len = len(motif)
        pos = protein.find(motif)
        while pos != -1:
            # Adjust confidence by motif length (longer motifs are more
            # specific).  10-mer = base confidence, 8-mer = -0.05, 12-mer = +0.05
            length_adj = (motif_len - 10) * 0.025
            confidence = min(1.0, max(0.0, sig["confidence"] + length_adj))

            matches.append(HazardMatch(
                category=sig["category"],
                name=sig["name"],
                position=pos,
                matched_sequence=motif,
                confidence=round(confidence, 3),
                source=sig["source"],
                match_type="exact",
                distance=0,
                strand="forward",
                substitutions=[],
            ))
            flagged_categories.add(sig["category"])
            exact_positions.setdefault(sig["name"], set()).add(pos)

            # Look for additional occurrences
            pos = protein.find(motif, pos + 1)

        # ── Fuzzy matching (Hamming) for short protein motifs ────────────
        if motif_len < 15:
            hamming_results = _fuzzy_match_hamming(protein, motif, max_distance=2)
            for fpos, fdist, fsubs in hamming_results:
                # Skip fuzzy matches at positions already covered by exact matches
                if sig["name"] in exact_positions and fpos in exact_positions[sig["name"]]:
                    continue
                # Reduce confidence for fuzzy matches
                confidence = min(1.0, max(0.0, sig["confidence"] - 0.10 * fdist))
                length_adj = (motif_len - 10) * 0.025
                confidence = min(1.0, max(0.0, confidence + length_adj))

                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=fpos,
                    matched_sequence=protein[fpos:fpos + motif_len],
                    confidence=round(confidence, 3),
                    source=sig["source"],
                    match_type="fuzzy",
                    distance=fdist,
                    strand="forward",
                    substitutions=fsubs,
                ))
                flagged_categories.add(sig["category"])

            # ── Fuzzy matching (edit distance) for short protein motifs ───
            edit_results = _fuzzy_match_edit_distance(protein, motif, max_distance=1)
            for epos, edist in edit_results:
                # Skip if an exact match already covers this position
                if sig["name"] in exact_positions and epos in exact_positions[sig["name"]]:
                    continue
                # Skip if a Hamming fuzzy match already covers this position
                existing_fuzzy = [
                    m for m in matches
                    if m.name == sig["name"] and m.position == epos and m.match_type == "fuzzy"
                ]
                if existing_fuzzy:
                    continue
                confidence = min(1.0, max(0.0, sig["confidence"] - 0.10 * edist))
                length_adj = (motif_len - 10) * 0.025
                confidence = min(1.0, max(0.0, confidence + length_adj))

                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=epos,
                    matched_sequence=protein[epos:epos + motif_len],
                    confidence=round(confidence, 3),
                    source=sig["source"],
                    match_type="fuzzy",
                    distance=edist,
                    strand="forward",
                    substitutions=[],
                ))
                flagged_categories.add(sig["category"])

    # ── Screen DNA against nucleotide patterns ───────────────────────────
    if dna:
        # Forward strand
        for sig in _DNA_SIGNATURES:
            motif = sig["motif"].upper()
            pos = dna.find(motif)
            while pos != -1:
                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=pos,
                    matched_sequence=motif,
                    confidence=sig["confidence"],
                    source=sig["source"],
                    match_type="exact",
                    distance=0,
                    strand="forward",
                    substitutions=[],
                ))
                flagged_categories.add(sig["category"])

                pos = dna.find(motif, pos + 1)

        # Reverse complement strand
        rc_dna = reverse_complement(dna)
        for sig in _DNA_SIGNATURES:
            motif = sig["motif"].upper()
            pos = rc_dna.find(motif)
            while pos != -1:
                # Map position back to the original strand
                original_pos = len(rc_dna) - pos - len(motif)
                matches.append(HazardMatch(
                    category=sig["category"],
                    name=sig["name"],
                    position=original_pos,
                    matched_sequence=motif,
                    confidence=sig["confidence"],
                    source=sig["source"],
                    match_type="reverse_complement",
                    distance=0,
                    strand="reverse",
                    substitutions=[],
                ))
                flagged_categories.add(sig["category"])

                pos = rc_dna.find(motif, pos + 1)

    # ── Determine risk level ─────────────────────────────────────────────
    if not matches:
        risk_level: RiskLevel = "none"
    else:
        match_risks = [sig_risk_for_match(m) for m in matches]
        risk_level = _max_risk(*match_risks)

    is_hazardous = risk_level in ("medium", "high", "critical")

    # ── Build recommendations ────────────────────────────────────────────
    recommendations = _build_recommendations(risk_level, flagged_categories, matches)

    report = BiosecurityReport(
        is_hazardous=is_hazardous,
        risk_level=risk_level,
        flagged_categories=sorted(flagged_categories),
        matches=matches,
        recommendations=recommendations,
    )

    logger.info(
        "Biosecurity screening complete: risk=%s, matches=%d, categories=%s",
        risk_level, len(matches), flagged_categories,
    )

    return report


def sig_risk_for_match(match: HazardMatch) -> str:
    """Look up the risk level for a match's signature.

    Fuzzy matches have their risk downgraded:
      - distance 1 → one level below the signature risk
      - distance 2 → two levels below the signature risk
    Reverse complement matches keep the original risk level.
    """
    base_risk: str | None = None
    for sig in _HAZARD_SIGNATURES:
        if sig["name"] == match.name and sig["category"] == match.category:
            base_risk = sig["risk"]
            break

    if base_risk is None:
        # Fallback: infer from category
        _category_default_risk = {
            "select_agent": "critical",
            "viral_surface": "high",
            "antibiotic_resistance": "medium",
            "oncogene": "low",
        }
        base_risk = _category_default_risk.get(match.category, "low")

    if match.match_type == "reverse_complement":
        return base_risk

    if match.match_type == "fuzzy" and match.distance > 0:
        # Fuzzy matches get fixed downgraded risk levels:
        #   distance 1 → "medium"
        #   distance 2+ → "low"
        if match.distance == 1:
            return "medium"
        return "low"

    return base_risk


def _build_recommendations(
    risk_level: RiskLevel,
    flagged_categories: set[str],
    matches: list[HazardMatch],
) -> list[str]:
    """Build actionable recommendations based on screening results."""
    recs: list[str] = []

    if risk_level == "none":
        recs.append("No biosecurity concerns detected. Proceed with optimization.")
        return recs

    if "select_agent" in flagged_categories:
        toxin_names = sorted({m.name for m in matches if m.category == "select_agent"})
        recs.append(
            f"CRITICAL: Select agent toxin signature(s) detected: {', '.join(toxin_names)}. "
            f"Optimization is blocked. Contact your institutional biosafety officer (IBO) "
            f"and verify compliance with 42 CFR Part 73 before proceeding."
        )

    if "viral_surface" in flagged_categories:
        viral_names = sorted({m.name for m in matches if m.category == "viral_surface"})
        recs.append(
            f"Viral surface protein signature(s) detected: {', '.join(viral_names)}. "
            f"Verify the sequence is intended for legitimate vaccine/therapeutic development. "
            f"Review dual-use research of concern (DURC) policies."
        )

    if "antibiotic_resistance" in flagged_categories:
        ar_names = sorted({m.name for m in matches if m.category == "antibiotic_resistance"})
        recs.append(
            f"Antibiotic resistance marker(s) detected: {', '.join(ar_names)}. "
            f"Ensure appropriate containment and disposal procedures. "
            f"Verify compliance with NIH Guidelines for Research Involving Recombinant DNA."
        )

    if "oncogene" in flagged_categories:
        onco_names = sorted({m.name for m in matches if m.category == "oncogene"})
        recs.append(
            f"Oncogene/growth factor signature(s) detected: {', '.join(onco_names)}. "
            f"Verify intended use is for legitimate research. "
            f"Review institutional gene therapy oversight requirements."
        )

    # Mention fuzzy/homology matches
    fuzzy_matches = [m for m in matches if m.match_type == "fuzzy"]
    if fuzzy_matches:
        fuzzy_names = sorted({m.name for m in fuzzy_matches})
        recs.append(
            f"Fuzzy homology match(es) detected (substitution/indel tolerance): "
            f"{', '.join(fuzzy_names)}. "
            f"These are lower-confidence matches — review carefully before proceeding."
        )

    # Mention reverse complement matches
    rc_matches = [m for m in matches if m.match_type == "reverse_complement"]
    if rc_matches:
        rc_names = sorted({m.name for m in rc_matches})
        recs.append(
            f"Reverse complement / anti-sense match(es) detected on the opposite strand: "
            f"{', '.join(rc_names)}. "
            f"These indicate hazardous sequences on the reverse strand."
        )

    if risk_level in ("high", "critical"):
        recs.append(
            "Optimization BLOCKED due to high/critical biosecurity risk. "
            "Resolve all flagged issues or obtain explicit institutional approval."
        )
    elif risk_level == "medium":
        recs.append(
            "Optimization allowed with warning. Review flagged items and "
            "ensure compliance with institutional biosafety protocols."
        )
    elif risk_level == "low":
        recs.append(
            "Low-risk biosecurity flags detected. No action required but "
            "review recommended."
        )

    return recs


# ─────────────────────────────────────────────────────────────────────────────
# Integration hook
# ─────────────────────────────────────────────────────────────────────────────

def check_biosecurity_before_optimize(
    protein: str,
    organism: str = "",
    dna: str = "",
    biosecurity_mode: Optional[BiosecurityMode] = None,
    skip_biosecurity_check: bool = False,
) -> BiosecurityScreeningResult:
    """Biosecurity gate called at the start of optimization.

    This function screens the input sequence and enforces hard-stop or
    warning behavior depending on the risk level and the effective
    biosecurity mode:

    - ``enforce``: ``critical`` or ``high`` risk raises
      :class:`BiosecurityError`; ``medium`` emits :class:`UserWarning`.
    - ``warn``: all risk levels emit warnings but no
      :class:`BiosecurityError` is raised.
    - ``off``: screening is skipped entirely; returns a clean result
      noting that screening was skipped.

    If *biosecurity_mode* is ``None`` (the default), the mode is read
    from the ``BIOCOMPILER_BIOSECURITY_MODE`` environment variable via
    :func:`get_biosecurity_mode`.

    Parameters
    ----------
    protein : str
        Protein sequence to screen.
    organism : str, optional
        Target organism (for context in log messages).
    dna : str, optional
        DNA sequence to screen for resistance markers.
    biosecurity_mode : BiosecurityMode or None, optional
        Explicit override for the biosecurity mode.  When ``None``,
        the mode is read from the environment variable.
    skip_biosecurity_check : bool, optional
        If ``True``, skip biosecurity screening entirely and return a
        passed result.  Default is ``False``.  This should only be used
        for testing or when the user explicitly opts out.

    Returns
    -------
    BiosecurityScreeningResult
        The screening result with pass/fail status and details.

    Raises
    ------
    BiosecurityError
        If mode is ``enforce`` and risk_level is ``critical`` or ``high``.
    ValueError
        If *protein* is empty or whitespace-only.
    """
    if not protein or not protein.strip():
        raise ValueError("Protein sequence must not be empty")

    # Skip check: always return passed
    if skip_biosecurity_check:
        return BiosecurityScreeningResult(
            passed=True,
            screened_sequence_length=len(protein.strip()),
        )

    if biosecurity_mode is None:
        biosecurity_mode = get_biosecurity_mode()

    # Off mode: skip screening entirely
    if biosecurity_mode == "off":
        return BiosecurityScreeningResult(
            passed=True,
            screened_sequence_length=len(protein.strip()),
        )

    protein_upper = protein.upper().strip()

    # ── Check legacy pathogen signatures (exact substring match) ─────────
    flagged_pathogens: list[str] = []
    risk_levels: list[str] = []
    match_details: list[str] = []
    kmer_scores: dict[str, float] = {}

    for sig_seq, pathogen_name, risk_level, description in _PATHOGEN_SIGNATURES:
        if sig_seq.upper() in protein_upper:
            pos = protein_upper.find(sig_seq.upper())
            flagged_pathogens.append(pathogen_name)
            risk_levels.append(risk_level)
            match_details.append(
                f"{description}: matched at position {pos}"
            )
            kmer_scores[pathogen_name] = 1.0

    # ── Also run the motif-based screening ──────────────────────────────
    report = screen_hazardous_sequence(protein, dna)

    # Merge motif-based findings into the result
    for match in report.matches:
        # Map motif match names to pathogen names where possible
        pathogen_name = _MOTIF_TO_PATHOGEN.get(match.name, match.name)
        if pathogen_name not in flagged_pathogens:
            flagged_pathogens.append(pathogen_name)
            # Map category to risk level based on biosecurity classification
            _CATEGORY_RISK = {
                "SELECT_AGENT": "CRITICAL",
                "VIRAL_SURFACE": "HIGH",
                "ANTIBIOTIC_RESISTANCE": "MEDIUM",
                "ONCOGENE": "HIGH",
                "TOXIN": "CRITICAL",
            }
            cat_risk = _CATEGORY_RISK.get(match.category.upper(), "MEDIUM")
            risk_levels.append(cat_risk.lower() if cat_risk in ("CRITICAL", "HIGH") else "medium")
            match_details.append(
                f"{match.name}: {match.match_type} match at position {match.position}"
            )
        # Compute k-mer similarity for this pathogen
        if pathogen_name not in kmer_scores:
            query_kmers = _extract_kmers(protein_upper)
            sig_kmers = _extract_kmers(match.matched_sequence)
            sim = _compute_kmer_similarity(query_kmers, sig_kmers)
            kmer_scores[pathogen_name] = round(sim, 4)

    passed = len(flagged_pathogens) == 0

    # Collect categories and matches from the motif report for compatibility
    flagged_categories = list({
        match.category for match in report.matches
    }) if report.matches else []
    # Also include category from legacy pathogen signatures
    for sig_seq, pathogen_name, risk_level, description in _PATHOGEN_SIGNATURES:
        if sig_seq.upper() in protein_upper and pathogen_name not in flagged_categories:
            cat = "select_agent" if risk_level.lower() == "critical" else "hazardous"
            if cat not in flagged_categories:
                flagged_categories.append(cat)

    result = BiosecurityScreeningResult(
        passed=passed,
        flagged_pathogens=flagged_pathogens,
        risk_levels=risk_levels,
        match_details=match_details,
        kmer_scores=kmer_scores,
        screened_sequence_length=len(protein.strip()),
        flagged_categories=flagged_categories,
        matches=list(report.matches) if report.matches else [],
        recommendations=report.recommendations if hasattr(report, "recommendations") else [],
    )

    # Enforce mode: raise BiosecurityError for critical/high risk
    if biosecurity_mode == "enforce" and not passed:
        has_critical_or_high = any(
            rl.upper() in ("CRITICAL", "HIGH") for rl in risk_levels
        )
        if has_critical_or_high:
            logger.error(
                "Biosecurity gate BLOCKED optimization: pathogens=%s, "
                "organism=%s, protein_len=%d",
                flagged_pathogens, organism, len(protein),
            )
            err = BiosecurityError(
                protein=protein,
                flagged_pathogens=flagged_pathogens,
                risk_levels=risk_levels,
                match_details=match_details,
                flagged_categories=flagged_categories,
                matches=list(report.matches) if report.matches else [],
            )
            err.report = result  # attach screening result for introspection
            raise err

        # Medium risk in enforce mode: emit UserWarning (proceeds, but alerts)
        has_medium = any(
            rl.upper() == "MEDIUM" for rl in risk_levels
        )
        if has_medium:
            warnings.warn(
                f"Biosecurity screening detected medium-risk hazard(s): "
                f"{', '.join(flagged_pathogens)}. "
                f"Optimization is NOT blocked, but review is recommended.",
                UserWarning,
                stacklevel=2,
            )

    # Warn mode: emit warnings for any non-passed result
    if not passed and biosecurity_mode == "warn":
        warnings.warn(
            f"Biosecurity screening detected hazard(s): "
            f"{', '.join(flagged_pathogens)}. "
            f"Mode is 'warn' so optimization is NOT blocked, but review is strongly recommended.",
            UserWarning,
            stacklevel=2,
        )

    return result
