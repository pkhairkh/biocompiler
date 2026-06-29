"""
BioCompiler Type System — Restriction-Site Predicate Checks
===========================================================
Restriction-site detection predicate check (exact IUPAC site match with
minimum length guard).

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

def check_no_restriction_site(seq: str, enzymes: List[str], min_site_length: int = _RESTRICTION_SITE_MIN_LENGTH) -> PredicateResult:
    """Predicate 4: No restriction enzyme recognition sites.

    Only checks restriction sites that are >= min_site_length bp (default 6bp).
    Short 4bp restriction sites are too common in coding sequences and
    cause excessive false positives.

    Also handles cross-codon sites that span 3+ codons properly.
    """
    from biocompiler.sequence.restriction_sites import get_recognition_site
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




__all__ = [
    "check_no_restriction_site",
]
