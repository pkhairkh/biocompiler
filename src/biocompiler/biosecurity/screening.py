"""
Main biosecurity screening API.

Contains the core screening function, the optimization gate hook,
risk classification, recommendation builder, and mode detection.

Screening pipeline:
  1. **Motif matching** (always): Sliding-window substring and fuzzy matching
     against built-in peptide and nucleotide hazard signatures.
  2. **BLAST homology search** (optional, when BLAST+ is available): If
     NCBI BLAST+ (blastn/blastp) is installed and a local database is
     configured via ``BIOCOMPILER_BLAST_DB_PATH``, a homology search is
     run as a second pass after motif matching.  BLAST hits that exceed
     the identity and e-value thresholds are merged into the report.
     If BLAST+ is not available, this step is skipped gracefully with
     no effect on the motif-based results.
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
from biocompiler.biosecurity.types import (
    BiosecurityMode,
    BiosecurityReport,
    BiosecurityScreeningResult,
    HazardMatch,
    RiskLevel,
)

from biocompiler.shared.exceptions import BiosecurityError

# Lazy import for BLAST screening — avoids hard dependency on blast_screening
# module and ensures graceful skip when BLAST+ is not installed.
try:
    from .blast_screening import BlastScreener as _BlastScreener
except ImportError:
    _BlastScreener = None  # type: ignore[assignment,misc]

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

def screen_hazardous_sequence(protein: str, dna: str = "", fast_mode: bool = False) -> BiosecurityReport:
    """Screen a protein (and optional DNA) sequence against known hazardous
    signatures.

    The screening pipeline has two passes:

    **Pass 1 — Motif matching** (always runs):
      Protein screening uses sliding-window substring matching against
      short peptide motifs (8-12 aa).  Fuzzy matching via Hamming distance
      (1-2 substitutions) and Levenshtein edit distance (1 edit) is applied
      for motifs < 15 aa.  DNA screening uses substring matching against
      nucleotide patterns (15-21 nt) for antibiotic resistance markers, plus
      reverse complement screening.

    **Pass 2 — BLAST homology search** (optional):
      If NCBI BLAST+ is installed and a local database is configured
      (via ``BIOCOMPILER_BLAST_DB_PATH`` environment variable), a homology
      search is run as a second pass.  BLAST hits that exceed identity and
      e-value thresholds are added to the report as additional matches.
      If BLAST+ is not available, this step is skipped gracefully with a
      debug log message — the motif-based results remain valid.

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
    - Confidence scoring accounts for motif length (longer = higher
      confidence) and exact match position context.
    - BLAST hits are attributed a confidence of 0.8 by default and
      classified as "blast_homology" match type.
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
        if motif_len < 15 and not fast_mode:
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
            if not fast_mode:
              edit_results = _fuzzy_match_edit_distance(protein, motif, max_distance=1)
            else:
              edit_results = []
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

    # ── Pass 2: BLAST homology search (optional) ────────────────────────
    # If BLAST+ is available and a database is configured, run a homology
    # search as a second pass after motif matching.  This catches hazards
    # that motif matching may miss (e.g., novel toxin variants with
    # sufficient sequence similarity to known hazardous sequences).
    # If BLAST+ is not available, this step is skipped gracefully.
    blast_hit_count = 0
    if _BlastScreener is not None:
        try:
            blast_screener = _BlastScreener()
            if blast_screener.is_blast_available():
                # Screen DNA if available, otherwise screen protein
                if dna:
                    blast_result = blast_screener.screen_sequence(dna)
                else:
                    blast_result = blast_screener.screen_protein(protein)

                for hit in blast_result.concerning_hits:
                    # Determine category based on pathogen/toxin flags
                    if hit.is_toxin and hit.is_pathogen:
                        blast_category = "select_agent"
                    elif hit.is_pathogen:
                        blast_category = "viral_surface"
                    elif hit.is_toxin:
                        blast_category = "select_agent"
                    else:
                        blast_category = "hazardous"

                    # BLAST hits use position 0 (alignment-level position
                    # is not directly comparable to motif position).
                    matches.append(HazardMatch(
                        category=blast_category,
                        name=hit.subject_id,
                        position=0,
                        matched_sequence=(
                            f"BLAST homology: {hit.subject_organism} "
                            f"id={hit.identity_percent:.1f}% "
                            f"e={hit.e_value:.1e}"
                        ),
                        confidence=0.8,
                        source=blast_result.screening_database,
                        match_type="blast_homology",
                        distance=0,
                        strand="forward",
                        substitutions=[],
                    ))
                    flagged_categories.add(blast_category)
                    blast_hit_count += 1

                logger.debug(
                    "BLAST screening second pass: %d concerning hits from %s",
                    blast_hit_count, blast_result.screening_database,
                )
            else:
                logger.debug(
                    "BLAST+ not available; skipping BLAST homology second pass."
                )
        except Exception as exc:
            # BLAST errors must never break the screening pipeline.
            logger.warning(
                "BLAST screening failed (non-fatal); motif-based results "
                "remain valid. Error: %s", exc,
            )
    else:
        logger.debug(
            "blast_screening module not importable; skipping BLAST second pass."
        )

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
        "Biosecurity screening complete: risk=%s, matches=%d (%d BLAST), categories=%s",
        risk_level, len(matches), blast_hit_count, flagged_categories,
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
        # Fallback: infer from category (Fix: # `australia_group` was previously omitted from this dict,
        # causing AG signatures (cholera toxin, saxitoxin, perfringens
        # epsilon toxin) to silently fall through to the "low" default
        # even though every AG entry in the DB has risk="high".)
        _category_default_risk = {
            "select_agent": "critical",
            "australia_group": "high",
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
        #
        # EXCEPTION — select_agent toxins (anthrax LF, botulinum, ricin,
        # abrin, shiga, diphtheria, tetanus, SEB, conotoxins, etc.):
        # These are CDC Select Agents (42 CFR Part 73) — the highest-
        # consequence hazards in the database.  A distance-2 fuzzy match
        # is still a strong near-miss signal (e.g. a single wobble-codon
        # variant or a sequencing error away from a known lethal toxin).
        # Downgrading to "low" would silently bypass is_hazardous (which
        # requires risk ∈ {medium, high, critical}) and let a real
        # select-agent near-miss through the gate.  Floor the risk at
        # "medium" so any fuzzy match on a select_agent toxin triggers
        # review.  See GAP-1 (anthrax LF distance-2 fix).
        if match.category == "select_agent":
            return "medium"
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

    # Fix: previously the australia_group branch was
    # missing from _build_recommendations, so AG hits (cholera toxin,
    # saxitoxin, perfringens epsilon toxin) produced no category-specific
    # guidance.  Australia Group covers dual-use biological agents under
    # export-control lists (e.g. AG Common Control List); risk is "high".
    if "australia_group" in flagged_categories:
        ag_names = sorted({m.name for m in matches if m.category == "australia_group"})
        recs.append(
            f"Australia Group export-controlled agent signature(s) detected: {', '.join(ag_names)}. "
            f"Verify compliance with the Australia Group Common Control List and "
            f"applicable national export-control regulations before synthesis."
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
    # NOTE (prior fix / C12): Previously this used fast_mode=True, which skips
    # fuzzy (approximate) matching and only does exact substring lookups.
    # That caused the anthrax LF reference to pass through the pre-opt gate
    # unflagged (README claimed 8/8 biosecurity detection, actual was 7/8).
    # The pre-optimization gate must use the SAME matching fidelity as the
    # post-optimization screen so that hazardous signatures are caught
    # BEFORE any optimization work is invested. fast_mode=False enables the
    # full fuzzy/k-mer matching path.
    report = screen_hazardous_sequence(protein, dna, fast_mode=False)

    # ── BLAST homology screening: SKIPPED in pre-optimization ─────────
    # BLAST is expensive (~0.18s) and the input protein is from a known source.
    # Full BLAST screening runs post-optimization on the optimized DNA.
    # (Commented out for performance — uncomment to enable pre-opt BLAST)
    # try:
    #     from .blast_integration import check_biosecurity_blast
    #     blast_report = check_biosecurity_blast(protein, dna=dna)
    #     if blast_report and blast_report.is_hazardous:
    #         report = BiosecurityReport(
    #             is_hazardous=True,
    #             risk_level=_max_risk(report.risk_level, blast_report.risk_level),
    #             flagged_categories=list(set(report.flagged_categories + blast_report.flagged_categories)),
    #             matches=report.matches + blast_report.matches,
    #             recommendations=report.recommendations + blast_report.recommendations,
    #         )
    # except (ImportError, AttributeError, RuntimeError, OSError) as blast_exc:
    #     logger.debug("BLAST screening skipped: %s", blast_exc)

    # Merge motif-based findings into the result
    for match in report.matches:
        # Map motif match names to pathogen names where possible
        pathogen_name = _MOTIF_TO_PATHOGEN.get(match.name, match.name)
        if pathogen_name not in flagged_pathogens:
            flagged_pathogens.append(pathogen_name)
            # prior fix / H16: Previously a hardcoded `_CATEGORY_RISK` dict
            # escalated ALL oncogene matches to "high", ignoring the
            # database's `risk="low"` field on oncogene signatures like
            # VEGF_receptor / EGFR_activation / p53_DNA_binding.  This
            # hard-blocked legitimate oncogene research proteins in
            # enforce mode.  Use the DB-aware `sig_risk_for_match()`
            # instead, which (1) looks up the signature's actual `risk`
            # field, (2) applies the fuzzy-distance downgrade
            # (distance-1 -> medium, distance-2 -> low) and (3) preserves
            # the GAP-1 select_agent distance-2 floor at "medium".
            cat_risk = sig_risk_for_match(match)
            risk_levels.append(cat_risk)
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
