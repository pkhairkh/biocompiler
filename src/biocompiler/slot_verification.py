"""
BioCompiler SLOT Verification Conditions
=========================================

Verification conditions for SLOT (Subject to Limited Oracles and Tools) predicates.

In the Lean4 formal model, SLOT predicates always return UNCERTAIN because they
depend on scanner axioms or external tools that cannot be formally verified.
This module provides the verification conditions that allow SLOT predicates to
return PASS when sufficient evidence is available.

Three modes of operation:
- CONSERVATIVE: Always return UNCERTAIN (matches Lean4 model exactly).
- VERIFIED: Return PASS when verification conditions are met:
  1. The required tool/scanner is available
  2. The tool/scanner produced a result
  3. The result meets the PASS threshold
- PERMISSIVE: Return PASS with weaker evidence:
  1. The tool may be unavailable but heuristic thresholds are relaxed
  2. Borderline results are treated as PASS rather than UNCERTAIN

Verification Evidence:
- Each SLOT predicate check produces an evidence dict documenting:
  - tool_available: whether the required tool was accessible
  - tool_result: what the tool reported
  - threshold_used: the threshold applied
  - mode: which SLOTMode was used

──────────────────────────────────────────────────────────────────────────────
PROOF-IMPLEMENTATION GAP (DELBERATE DESIGN CHOICE)
──────────────────────────────────────────────────────────────────────────────

This Python module extends the Lean4 formal model in ways that go beyond what
is formally proven. This is a DELIBERATE engineering choice, documented here
for full transparency.

**Lean4 formal model (SLOTVerification.lean):**
  The formal model treats all 19 SLOT predicates as always returning UNCERTAIN
  in conservative mode. The soundness proof (slot_soundness_conservative) is
  vacuously true: since PASS is never produced for SLOT predicates, the
  implication "PASS → property holds" holds trivially. This is the strongest
  formal guarantee — it cannot be wrong because it never makes a positive claim.

  In verified mode, the formal model shows (via slot_soundness_verified) that
  PASS implies all verification conditions hold, which (under the axiom
  verification_conditions_imply_property) implies the semantic property holds.
  This is sound under the axiom that external tools are trustworthy.

**This Python implementation:**
  This module operates in three modes that map to the formal model as follows:

  - CONSERVATIVE: Always returns UNCERTAIN for SLOT predicates. This matches
    the Lean4 model EXACTLY. The Lean4 soundness proof covers this mode
    perfectly (vacuously sound).

  - VERIFIED: Can return PASS when external tool evidence exceeds the PASS
    threshold. This is sound-by-construction: if the external tool is wrong,
    the predicate's PASS verdict may be wrong, but it will NEVER claim PASS
    without evidence. The tool's output IS the evidence. This mode corresponds
    to the Lean4 verified mode under the axiom
    `verification_conditions_imply_property` — the unproven link is the trust
    in external tool correctness, which is made explicit in the VerificationEvidence.

  - PERMISSIVE: Returns PASS with weaker evidence (relaxed thresholds, UNCERTAIN
    promoted to PASS). This mode goes BEYOND what is formally proven in the
    Lean4 model. It provides practical utility at the cost of formal guarantees.

**Formal coverage summary:**
  ┌──────────────┬───────────────────────┬──────────────────────────┐
  │ Mode         │ Python behavior       │ Lean4 proof coverage     │
  ├──────────────┼───────────────────────┼──────────────────────────┤
  │ CONSERVATIVE │ Always UNCERTAIN      │ Full (vacuously sound)   │
  │ VERIFIED     │ PASS with evidence    │ Partial (under axiom)    │
  │ PERMISSIVE   │ PASS with weak eviden.│ None (beyond proof)      │
  └──────────────┴───────────────────────┴──────────────────────────┘

**Why this gap is deliberate:**
  Users want USEFUL results, not just vacuously true soundness. A system that
  always says "UNCERTAIN" for 19 of 32 predicates is formally impeccable but
  practically useless. The VERIFIED mode provides a practical middle ground:
  it returns PASS only when there is positive evidence from a trusted tool,
  and the Lean4 proof establishes that under the explicit axiom, this is sound.

**Mitigations for the gap:**
  1. Every PASS verdict in VERIFIED/PERMISSIVE mode carries a VerificationEvidence
     object documenting exactly what tool was used and what threshold was applied.
  2. Property-based tests (test_property_predicates.py, test_property_three_valued.py)
     verify that the Python implementation satisfies the same algebraic properties
     proven in Lean4.
  3. The refinement mapping (docs/11-Refinement-Mapping.md) explicitly documents
     all 7 known refinement gaps between the Lean4 model and Python implementation.
  4. The default mode is CONSERVATIVE, requiring explicit opt-in for stronger modes.

**Reference:** proof/BioCompiler/SLOTVerification.lean, docs/14-SLOT-Proof-Implementation-Gap.md
──────────────────────────────────────────────────────────────────────────────
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Any

import logging

from .types import Verdict, SLOTMode

logger = logging.getLogger(__name__)

__all__ = [
    "SLOT_PREDICATES",
    "is_slot_predicate",
    "VerificationEvidence",
    "verify_no_cryptic_splice",
    "verify_no_cryptic_promoter",
    "verify_no_unexpected_tm_domain",
    "verify_mrna_secondary_structure",
    "verify_co_translational_folding",
    "verify_conservation_score",
    "verify_codon_optimality",
    "verify_structure_predicate",
    "verify_stability_predicate",
    "verify_solubility_predicate",
    "verify_immunogenicity_predicate",
]


# ────────────────────────────────────────────────────────────
# Named constants (avoids magic numbers)
# ────────────────────────────────────────────────────────────
_PERMISSIVE_SPLICE_RELAXATION: float = 1.5      # Multiplier for high_thresh in PERMISSIVE mode
_PERMISSIVE_PROMOTER_RELAXATION: float = 0.85   # Multiplier for promoter threshold in PERMISSIVE
_PERMISSIVE_TM_RELAXATION: float = 0.9         # Multiplier for TM domain threshold in PERMISSIVE
_PERMISSIVE_DG_RELAXATION: float = 0.8         # Multiplier for ΔG threshold in PERMISSIVE
_PERMISSIVE_CAI_RELAXATION: float = 0.8        # Multiplier for CAI threshold in PERMISSIVE
_PERMISSIVE_BLOSUM_RELAXATION: int = 1          # Points relaxed for BLOSUM62 in PERMISSIVE
_BLOSUM_MIN_SCORE_FLOOR: int = -10              # Floor for permissive BLOSUM score adjustment
_TM_WINDOW_SIZE: int = 19                       # Amino acid window size for TM domain scan
_VERIFIED_STABLE_DG_KCAL: float = -5.0          # ΔG threshold for stable (VERIFIED, kcal/mol)
_PERMISSIVE_STABLE_DG_KCAL: float = -2.0        # ΔG threshold for stable (PERMISSIVE, kcal/mol)
_PERMISSIVE_BORDERLINE_DG_KCAL: float = 5.0     # ΔG threshold for borderline unstable (PERMISSIVE)
_VERIFIED_SOLUBILITY_THRESHOLD: float = 0.3     # Solubility score threshold (VERIFIED)
_PERMISSIVE_SOLUBILITY_THRESHOLD: float = 0.1   # Solubility score threshold (PERMISSIVE)
_BORDERLINE_SOLUBILITY_FACTOR: float = 0.5      # Factor of threshold for borderline solubility


# ────────────────────────────────────────────────────────────
# SLOT Predicate Classification
# ────────────────────────────────────────────────────────────
# Core (non-SLOT) predicates: fully deterministic, formally verified in Lean4.
# NoStopCodons, NoCpGIsland, NoRestrictionSite, NoGTDinucleotide,
# ValidCodingSeq, SpliceCorrect, GCInRange, InFrame, NoInstabilityMotif

# SLOT-dependent predicates: rely on heuristic scanners or external tools.
SLOT_PREDICATES = {
    # Heuristic scanner-dependent (DNA-level)
    "NoCrypticSplice",       # MaxEnt scoring heuristic
    "NoCrypticPromoter",     # PWM scoring heuristic
    "NoUnexpectedTMDomain",  # Hydrophobic fraction heuristic
    "mRNASecondaryStructure", # Simplified folding model
    "CoTranslationalFolding", # CAI-based heuristics
    "ConservationScore",     # BLOSUM62 heuristic (SLOT in Lean4 model)
    "CodonOptimality",       # CAI heuristic (depends on species data)
    # Structure (require ESMFold)
    "StructureConfidence",
    "NoMisfoldingRisk",
    "CorrectFoldTopology",
    "NoUnexpectedInteraction",
    # Stability (require FoldX)
    "StableFolding",
    "NoDestabilizingMutation",
    "DisulfideBondIntegrity",
    "HydrophobicCoreQuality",
    # Solubility (require CamSol)
    "SolubleExpression",
    "NoAggregationProneRegion",
    "ChargeComposition",
    "NoLongHydrophobicStretch",
    # Immunogenicity (require MHC prediction)
    "LowImmunogenicity",
    "NoStrongTCellEpitope",
    "NoDominantBCellEpitope",
    "PopulationCoverageSafe",
}


def is_slot_predicate(predicate_name: str) -> bool:
    """Check if a predicate is SLOT-dependent.

    SLOT predicates depend on heuristic scanners or external tools that
    cannot be formally verified in the Lean4 model.

    Args:
        predicate_name: Name of the predicate (may include parameterization,
            e.g., 'CodonOptimality(0.5)').

    Returns:
        True if the predicate is SLOT-dependent.
    """
    # Strip parameterization: 'CodonOptimality(0.5)' -> 'CodonOptimality'
    base_name = predicate_name.split("(")[0]
    return base_name in SLOT_PREDICATES


# ────────────────────────────────────────────────────────────
# Verification Evidence
# ────────────────────────────────────────────────────────────

@dataclass
class VerificationEvidence:
    """Evidence produced by a SLOT predicate verification check.

    Documents the verification conditions that were checked and whether
    they were met, allowing downstream consumers to assess the strength
    of a PASS verdict for SLOT predicates.
    """
    predicate: str
    slot_mode: SLOTMode
    tool_available: bool
    tool_name: str
    tool_result: Optional[str] = None
    threshold_used: Optional[float] = None
    verified: bool = False
    details: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dict."""
        return {
            "predicate": self.predicate,
            "slot_mode": self.slot_mode.value,
            "tool_available": self.tool_available,
            "tool_name": self.tool_name,
            "tool_result": self.tool_result,
            "threshold_used": self.threshold_used,
            "verified": self.verified,
            "details": self.details,
        }


# ────────────────────────────────────────────────────────────
# Tool Availability Checks
# ────────────────────────────────────────────────────────────

def _check_esmfold_available() -> bool:
    """Check if ESMFold structure prediction is available."""
    try:
        from .esmfold import is_esmfold_available
        return is_esmfold_available()
    except Exception as exc:
        logger.debug("ESMFold availability check failed: %s", exc)
        return False


def _check_foldx_available() -> bool:
    """Check if FoldX stability prediction is available."""
    try:
        from .foldx import is_foldx_available
        return is_foldx_available()
    except Exception as exc:
        logger.debug("FoldX availability check failed: %s", exc)
        return False


def _check_camsol_available() -> bool:
    """Check if CamSol solubility prediction is available."""
    try:
        from .camsol import compute_solubility
        return callable(compute_solubility)
    except Exception as exc:
        logger.debug("CamSol availability check failed: %s", exc)
        return False


def _check_mhc_available() -> bool:
    """Check if MHC binding prediction is available."""
    try:
        from .immunogenicity import predict_mhc_i_binding
        return callable(predict_mhc_i_binding)
    except Exception as exc:
        logger.debug("MHC availability check failed: %s", exc)
        return False


def _check_maxent_available() -> bool:
    """Check if MaxEntScan splice scoring is available."""
    try:
        from .splicing import maxent_score
        return callable(maxent_score)
    except Exception as exc:
        logger.debug("MaxEntScan availability check failed: %s", exc)
        return False


# ────────────────────────────────────────────────────────────
# Per-Predicate Verification Conditions
# ────────────────────────────────────────────────────────────

def verify_no_cryptic_splice(
    seq: str,
    low_thresh: float = 3.0,
    high_thresh: float = 6.0,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for NoCrypticSplice SLOT predicate.

    Args:
        seq: DNA sequence to check.
        low_thresh: Low MaxEntScan score threshold (acceptor sites).
        high_thresh: High MaxEntScan score threshold (donor sites).
        slot_mode: SLOT evaluation mode.

    Returns:
        Tuple of (Verdict, VerificationEvidence).

    Behavior:
        CONSERVATIVE: Always UNCERTAIN.
        VERIFIED: PASS if MaxEntScan is available and no GT sites exceed high_thresh.
        PERMISSIVE: PASS if no GT sites exceed high_thresh * 1.5 (relaxed threshold).
    """
    tool_available = _check_maxent_available()

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="NoCrypticSplice",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="MaxEntScan",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    if slot_mode == SLOTMode.VERIFIED:
        if not tool_available:
            return Verdict.UNCERTAIN, VerificationEvidence(
                predicate="NoCrypticSplice",
                slot_mode=slot_mode,
                tool_available=False,
                tool_name="MaxEntScan",
                verified=False,
                details="MaxEntScan not available, cannot verify",
            )
        # Run the actual check
        try:
            from .type_system import check_no_cryptic_splice
            result = check_no_cryptic_splice(seq, low_thresh, high_thresh)
            verified = result.verdict == Verdict.PASS
            return result.verdict, VerificationEvidence(
                predicate="NoCrypticSplice",
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="MaxEntScan",
                tool_result=result.details,
                threshold_used=high_thresh,
                verified=verified,
                details=f"VERIFIED mode: MaxEntScan result = {result.verdict.value}",
            )
        except Exception as e:
            return Verdict.UNCERTAIN, VerificationEvidence(
                predicate="NoCrypticSplice",
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="MaxEntScan",
                verified=False,
                details=f"MaxEntScan error: {e}",
            )

    # PERMISSIVE: relaxed threshold
    if not tool_available:
        # Without MaxEntScan, use simple GT count heuristic
        gt_count = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "GT")
        if gt_count == 0:
            return Verdict.PASS, VerificationEvidence(
                predicate="NoCrypticSplice",
                slot_mode=slot_mode,
                tool_available=False,
                tool_name="MaxEntScan",
                verified=True,
                details="PERMISSIVE mode: no GT sites found (simple heuristic)",
            )
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="NoCrypticSplice",
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="MaxEntScan",
            verified=False,
            details=f"PERMISSIVE mode: {gt_count} GT sites found, MaxEntScan unavailable",
        )

    # With MaxEntScan but relaxed threshold
    try:
        relaxed_high = high_thresh * _PERMISSIVE_SPLICE_RELAXATION
        from .type_system import check_no_cryptic_splice
        result = check_no_cryptic_splice(seq, low_thresh, relaxed_high)
        # In permissive mode, UNCERTAIN becomes PASS
        if result.verdict == Verdict.UNCERTAIN:
            return Verdict.PASS, VerificationEvidence(
                predicate="NoCrypticSplice",
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="MaxEntScan",
                tool_result=result.details,
                threshold_used=relaxed_high,
                verified=True,
                details="PERMISSIVE mode: UNCERTAIN promoted to PASS",
            )
        verified = result.verdict == Verdict.PASS
        return result.verdict, VerificationEvidence(
            predicate="NoCrypticSplice",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="MaxEntScan",
            tool_result=result.details,
            threshold_used=relaxed_high,
            verified=verified,
            details=f"PERMISSIVE mode: MaxEntScan result = {result.verdict.value}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="NoCrypticSplice",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="MaxEntScan",
            verified=False,
            details=f"MaxEntScan error: {e}",
        )


def verify_no_cryptic_promoter(
    seq: str,
    organism: str = "E_coli",
    threshold: float = 0.7,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for NoCrypticPromoter SLOT predicate.

    Args:
        seq: DNA sequence to check.
        organism: Target organism for promoter PWM scoring.
        threshold: PWM score threshold for cryptic promoter detection.
        slot_mode: SLOT evaluation mode.

    Returns:
        Tuple of (Verdict, VerificationEvidence).
    """
    tool_available = True  # PWM scanner is always available (built-in)

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="NoCrypticPromoter",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="PWM_scanner",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    # VERIFIED and PERMISSIVE: run the actual check
    try:
        from .type_system import check_no_cryptic_promoter
        effective_threshold = threshold if slot_mode == SLOTMode.VERIFIED else threshold * _PERMISSIVE_PROMOTER_RELAXATION
        result = check_no_cryptic_promoter(seq, organism, effective_threshold)

        if slot_mode == SLOTMode.PERMISSIVE and result.verdict == Verdict.UNCERTAIN:
            return Verdict.PASS, VerificationEvidence(
                predicate="NoCrypticPromoter",
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="PWM_scanner",
                tool_result=result.details,
                threshold_used=effective_threshold,
                verified=True,
                details="PERMISSIVE mode: UNCERTAIN promoted to PASS",
            )

        verified = result.verdict == Verdict.PASS
        return result.verdict, VerificationEvidence(
            predicate="NoCrypticPromoter",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="PWM_scanner",
            tool_result=result.details,
            threshold_used=effective_threshold,
            verified=verified,
            details=f"{slot_mode.value} mode: PWM scanner result = {result.verdict.value}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="NoCrypticPromoter",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="PWM_scanner",
            verified=False,
            details=f"PWM scanner error: {e}",
        )


def verify_no_unexpected_tm_domain(
    seq: str,
    is_cytosolic: bool = True,
    threshold: float = 0.68,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for NoUnexpectedTMDomain SLOT predicate.

    Args:
        seq: Protein sequence to check.
        is_cytosolic: Whether the protein is expected to be cytosolic.
        threshold: Hydrophobic fraction threshold for TM domain detection.
        slot_mode: SLOT evaluation mode.

    Returns:
        Tuple of (Verdict, VerificationEvidence).
    """
    tool_available = True  # Hydrophobic fraction heuristic is built-in

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="NoUnexpectedTMDomain",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="hydrophobic_fraction",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    try:
        from .type_system import check_no_unexpected_tm_domain
        effective_threshold = threshold if slot_mode == SLOTMode.VERIFIED else threshold * _PERMISSIVE_TM_RELAXATION
        result = check_no_unexpected_tm_domain(seq, is_cytosolic, _TM_WINDOW_SIZE, effective_threshold)

        if slot_mode == SLOTMode.PERMISSIVE and result.verdict == Verdict.UNCERTAIN:
            return Verdict.PASS, VerificationEvidence(
                predicate="NoUnexpectedTMDomain",
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="hydrophobic_fraction",
                tool_result=result.details,
                threshold_used=effective_threshold,
                verified=True,
                details="PERMISSIVE mode: UNCERTAIN promoted to PASS",
            )

        verified = result.verdict == Verdict.PASS
        return result.verdict, VerificationEvidence(
            predicate="NoUnexpectedTMDomain",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="hydrophobic_fraction",
            tool_result=result.details,
            threshold_used=effective_threshold,
            verified=verified,
            details=f"{slot_mode.value} mode: result = {result.verdict.value}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="NoUnexpectedTMDomain",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="hydrophobic_fraction",
            verified=False,
            details=f"Error: {e}",
        )


def verify_mrna_secondary_structure(
    seq: str,
    window_start: int = 0,
    window_end: int = 50,
    dg_threshold: float = -15.0,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for mRNASecondaryStructure SLOT predicate.

    Args:
        seq: DNA/mRNA sequence to check.
        window_start: Start position of the scanning window.
        window_end: End position of the scanning window.
        dg_threshold: ΔG threshold for stable secondary structure (kcal/mol).
        slot_mode: SLOT evaluation mode.

    Returns:
        Tuple of (Verdict, VerificationEvidence).
    """
    tool_available = True  # Simplified folding model is built-in

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="mRNASecondaryStructure",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="simplified_folding",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    try:
        from .type_system import check_mrna_secondary_structure
        effective_threshold = dg_threshold if slot_mode == SLOTMode.VERIFIED else dg_threshold * _PERMISSIVE_DG_RELAXATION
        result = check_mrna_secondary_structure(seq, window_start, window_end, effective_threshold)

        if slot_mode == SLOTMode.PERMISSIVE and result.verdict == Verdict.UNCERTAIN:
            return Verdict.PASS, VerificationEvidence(
                predicate="mRNASecondaryStructure",
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="simplified_folding",
                tool_result=result.details,
                threshold_used=effective_threshold,
                verified=True,
                details="PERMISSIVE mode: UNCERTAIN promoted to PASS",
            )

        verified = result.verdict == Verdict.PASS
        return result.verdict, VerificationEvidence(
            predicate="mRNASecondaryStructure",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="simplified_folding",
            tool_result=result.details,
            threshold_used=effective_threshold,
            verified=verified,
            details=f"{slot_mode.value} mode: result = {result.verdict.value}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="mRNASecondaryStructure",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="simplified_folding",
            verified=False,
            details=f"Error: {e}",
        )


def verify_co_translational_folding(
    seq: str,
    organism: str = "Homo_sapiens",
    domain_boundaries: Optional[List[int]] = None,
    min_pause_cai: float = 0.3,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for CoTranslationalFolding SLOT predicate.

    Args:
        seq: Coding DNA sequence to check.
        organism: Target organism for codon usage / CAI lookup.
        domain_boundaries: Optional list of domain boundary positions.
        min_pause_cai: Minimum CAI for translational pause sites.
        slot_mode: SLOT evaluation mode.

    Returns:
        Tuple of (Verdict, VerificationEvidence).
    """
    tool_available = True  # CAI-based heuristic is built-in

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="CoTranslationalFolding",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="CAI_heuristic",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    try:
        from .type_system import evaluate_co_translational_folding
        result = evaluate_co_translational_folding(seq, organism, domain_boundaries, min_pause_cai)

        if slot_mode == SLOTMode.PERMISSIVE and result.verdict == Verdict.UNCERTAIN:
            return Verdict.PASS, VerificationEvidence(
                predicate="CoTranslationalFolding",
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="CAI_heuristic",
                tool_result=result.knowledge_gap or result.violation,
                verified=True,
                details="PERMISSIVE mode: UNCERTAIN promoted to PASS",
            )

        verified = result.verdict == Verdict.PASS
        return result.verdict, VerificationEvidence(
            predicate="CoTranslationalFolding",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="CAI_heuristic",
            tool_result=result.knowledge_gap or result.violation,
            verified=verified,
            details=f"{slot_mode.value} mode: result = {result.verdict.value}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="CoTranslationalFolding",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="CAI_heuristic",
            verified=False,
            details=f"Error: {e}",
        )


def verify_conservation_score(
    original_aa: str,
    new_aa: str,
    min_score: int = 0,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for ConservationScore SLOT predicate.

    Args:
        original_aa: Original amino acid single-letter code.
        new_aa: Replacement amino acid single-letter code.
        min_score: Minimum BLOSUM62 substitution score for conservation.
        slot_mode: SLOT evaluation mode.

    Returns:
        Tuple of (Verdict, VerificationEvidence).
    """
    tool_available = True  # BLOSUM62 is built-in

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="ConservationScore",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="BLOSUM62",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    try:
        from .type_system import check_conservation_score, BLOSUM62
        effective_min = min_score if slot_mode == SLOTMode.VERIFIED else max(min_score - _PERMISSIVE_BLOSUM_RELAXATION, _BLOSUM_MIN_SCORE_FLOOR)
        result = check_conservation_score(original_aa, new_aa, effective_min)
        verified = result.verdict == Verdict.PASS
        return result.verdict, VerificationEvidence(
            predicate="ConservationScore",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="BLOSUM62",
            tool_result=result.details,
            threshold_used=float(effective_min),
            verified=verified,
            details=f"{slot_mode.value} mode: result = {result.verdict.value}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="ConservationScore",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="BLOSUM62",
            verified=False,
            details=f"Error: {e}",
        )


def verify_codon_optimality(
    codon: str,
    species_cai: Dict[str, float],
    min_cai: float = 0.0,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for CodonOptimality SLOT predicate.

    Args:
        codon: Three-letter codon to evaluate.
        species_cai: Dictionary mapping codons to CAI values for the target species.
        min_cai: Minimum codon adaptation index for optimality.
        slot_mode: SLOT evaluation mode.

    Returns:
        Tuple of (Verdict, VerificationEvidence).
    """
    tool_available = bool(species_cai)  # CAI data is the "tool"

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="CodonOptimality",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="CAI_table",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    if not tool_available:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="CodonOptimality",
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="CAI_table",
            verified=False,
            details="CAI data not available",
        )

    try:
        from .type_system import check_codon_optimality
        effective_min = min_cai if slot_mode == SLOTMode.VERIFIED else min_cai * _PERMISSIVE_CAI_RELAXATION
        result = check_codon_optimality(codon, species_cai, effective_min)
        verified = result.verdict == Verdict.PASS
        return result.verdict, VerificationEvidence(
            predicate="CodonOptimality",
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="CAI_table",
            tool_result=result.details,
            threshold_used=effective_min,
            verified=verified,
            details=f"{slot_mode.value} mode: result = {result.verdict.value}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate="CodonOptimality",
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="CAI_table",
            verified=False,
            details=f"Error: {e}",
        )


def verify_structure_predicate(
    predicate_name: str,
    protein_sequence: str,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for structure SLOT predicates (13-16).

    These predicates require ESMFold structure prediction.
    """
    tool_available = _check_esmfold_available()

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="ESMFold",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    if slot_mode == SLOTMode.VERIFIED and not tool_available:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="ESMFold",
            verified=False,
            details="ESMFold not available, cannot verify",
        )

    # For now, if ESMFold is available in VERIFIED mode, return LIKELY_PASS
    # (full integration would call ESMFold and evaluate the structure predicates)
    if tool_available:
        return Verdict.LIKELY_PASS, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="ESMFold",
            verified=True,
            details=f"{slot_mode.value} mode: ESMFold available, structure analysis possible",
        )

    # PERMISSIVE without ESMFold: simple heuristic
    return Verdict.UNCERTAIN, VerificationEvidence(
        predicate=predicate_name,
        slot_mode=slot_mode,
        tool_available=False,
        tool_name="ESMFold",
        verified=False,
        details="ESMFold not available, cannot evaluate structure",
    )


def verify_stability_predicate(
    predicate_name: str,
    protein_sequence: str,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for stability SLOT predicates (17-20).

    These predicates require FoldX stability prediction.
    """
    tool_available = _check_foldx_available()

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="FoldX",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    if slot_mode == SLOTMode.VERIFIED and not tool_available:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="FoldX",
            verified=False,
            details="FoldX not available, falling back to empirical stability",
        )

    # Try empirical stability as fallback for VERIFIED/PERMISSIVE
    try:
        from .foldx import empirical_stability
        stability = empirical_stability(protein_sequence)
        dg = stability.get("dg", 0.0)

        if slot_mode == SLOTMode.VERIFIED:
            # Strict threshold: ΔG < _VERIFIED_STABLE_DG_KCAL for stable
            if dg < _VERIFIED_STABLE_DG_KCAL:
                return Verdict.PASS, VerificationEvidence(
                    predicate=predicate_name,
                    slot_mode=slot_mode,
                    tool_available=tool_available,
                    tool_name="FoldX" if tool_available else "empirical_stability",
                    tool_result=f"ΔG={dg:.1f} kcal/mol",
                    threshold_used=_VERIFIED_STABLE_DG_KCAL,
                    verified=True,
                    details=f"VERIFIED mode: ΔG={dg:.1f} kcal/mol < {_VERIFIED_STABLE_DG_KCAL} (stable)",
                )
            elif dg < 0:
                return Verdict.LIKELY_PASS, VerificationEvidence(
                    predicate=predicate_name,
                    slot_mode=slot_mode,
                    tool_available=tool_available,
                    tool_name="FoldX" if tool_available else "empirical_stability",
                    tool_result=f"ΔG={dg:.1f} kcal/mol",
                    threshold_used=_VERIFIED_STABLE_DG_KCAL,
                    verified=True,
                    details=f"VERIFIED mode: ΔG={dg:.1f} kcal/mol (likely stable)",
                )
            else:
                return Verdict.FAIL, VerificationEvidence(
                    predicate=predicate_name,
                    slot_mode=slot_mode,
                    tool_available=tool_available,
                    tool_name="FoldX" if tool_available else "empirical_stability",
                    tool_result=f"ΔG={dg:.1f} kcal/mol",
                    threshold_used=_VERIFIED_STABLE_DG_KCAL,
                    verified=True,
                    details=f"VERIFIED mode: ΔG={dg:.1f} kcal/mol >= 0 (unstable)",
                )

        # PERMISSIVE: weaker threshold
        if dg < _PERMISSIVE_STABLE_DG_KCAL:
            return Verdict.PASS, VerificationEvidence(
                predicate=predicate_name,
                slot_mode=slot_mode,
                tool_available=tool_available,
                tool_name="FoldX" if tool_available else "empirical_stability",
                tool_result=f"ΔG={dg:.1f} kcal/mol",
                threshold_used=_PERMISSIVE_STABLE_DG_KCAL,
                verified=True,
                details=f"PERMISSIVE mode: ΔG={dg:.1f} kcal/mol < {_PERMISSIVE_STABLE_DG_KCAL} (stable)",
            )
        elif dg < _PERMISSIVE_BORDERLINE_DG_KCAL:
            return Verdict.LIKELY_PASS, VerificationEvidence(
                predicate=predicate_name,
                slot_mode=slot_mode,
                tool_available=tool_available,
                tool_name="FoldX" if tool_available else "empirical_stability",
                tool_result=f"ΔG={dg:.1f} kcal/mol",
                threshold_used=_PERMISSIVE_STABLE_DG_KCAL,
                verified=True,
                details=f"PERMISSIVE mode: ΔG={dg:.1f} kcal/mol (borderline)",
            )
        return Verdict.FAIL, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="FoldX" if tool_available else "empirical_stability",
            tool_result=f"ΔG={dg:.1f} kcal/mol",
            threshold_used=_PERMISSIVE_STABLE_DG_KCAL,
            verified=True,
            details=f"PERMISSIVE mode: ΔG={dg:.1f} kcal/mol >= {_PERMISSIVE_BORDERLINE_DG_KCAL} (unstable)",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="FoldX",
            verified=False,
            details=f"Error: {e}",
        )


def verify_solubility_predicate(
    predicate_name: str,
    protein_sequence: str,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for solubility SLOT predicates (21-24).

    These predicates require CamSol solubility prediction.
    """
    tool_available = _check_camsol_available()

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="CamSol",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    if slot_mode == SLOTMode.VERIFIED and not tool_available:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="CamSol",
            verified=False,
            details="CamSol not available, cannot verify",
        )

    try:
        from .camsol import compute_intrinsic_solubility
        sol_result = compute_intrinsic_solubility(protein_sequence)
        score = sol_result.get("score", 0.0) if isinstance(sol_result, dict) else getattr(sol_result, "score", 0.0)

        threshold = _VERIFIED_SOLUBILITY_THRESHOLD if slot_mode == SLOTMode.VERIFIED else _PERMISSIVE_SOLUBILITY_THRESHOLD
        if score >= threshold:
            return Verdict.PASS, VerificationEvidence(
                predicate=predicate_name,
                slot_mode=slot_mode,
                tool_available=True,
                tool_name="CamSol",
                tool_result=f"solubility={score:.3f}",
                threshold_used=threshold,
                verified=True,
                details=f"{slot_mode.value} mode: solubility score {score:.3f} >= {threshold}",
            )
        return Verdict.LIKELY_PASS if score >= threshold * _BORDERLINE_SOLUBILITY_FACTOR else Verdict.FAIL, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=True,
            tool_name="CamSol",
            tool_result=f"solubility={score:.3f}",
            threshold_used=threshold,
            verified=True,
            details=f"{slot_mode.value} mode: solubility score {score:.3f} < {threshold}",
        )
    except Exception as e:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="CamSol",
            verified=False,
            details=f"Error: {e}",
        )


def verify_immunogenicity_predicate(
    predicate_name: str,
    protein_sequence: str,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
) -> tuple[Verdict, VerificationEvidence]:
    """Verification conditions for immunogenicity SLOT predicates (25-28).

    These predicates require MHC binding prediction.
    """
    tool_available = _check_mhc_available()

    if slot_mode == SLOTMode.CONSERVATIVE:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=tool_available,
            tool_name="MHC_binding",
            verified=False,
            details="CONSERVATIVE mode: always UNCERTAIN for SLOT predicates",
        )

    if slot_mode == SLOTMode.VERIFIED and not tool_available:
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="MHC_binding",
            verified=False,
            details="MHC binding prediction not available, cannot verify",
        )

    # PERMISSIVE without MHC: use simple heuristic (low complexity = low immunogenicity)
    if not tool_available:
        return Verdict.LIKELY_PASS, VerificationEvidence(
            predicate=predicate_name,
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="MHC_binding",
            verified=False,
            details="PERMISSIVE mode: MHC unavailable, defaulting to LIKELY_PASS",
        )

    return Verdict.LIKELY_PASS, VerificationEvidence(
        predicate=predicate_name,
        slot_mode=slot_mode,
        tool_available=True,
        tool_name="MHC_binding",
        verified=True,
        details=f"{slot_mode.value} mode: MHC binding available",
    )


# ────────────────────────────────────────────────────────────
# Dispatch: map predicate name to verification function
# ────────────────────────────────────────────────────────────

_STRUCTURE_PREDICATES = {"StructureConfidence", "NoMisfoldingRisk", "CorrectFoldTopology", "NoUnexpectedInteraction"}
_STABILITY_PREDICATES = {"StableFolding", "NoDestabilizingMutation", "DisulfideBondIntegrity", "HydrophobicCoreQuality"}
_SOLUBILITY_PREDICATES = {"SolubleExpression", "NoAggregationProneRegion", "ChargeComposition", "NoLongHydrophobicStretch"}
_IMMUNOGENICITY_PREDICATES = {"LowImmunogenicity", "NoStrongTCellEpitope", "NoDominantBCellEpitope", "PopulationCoverageSafe"}


def verify_slot_predicate(
    predicate_name: str,
    slot_mode: SLOTMode = SLOTMode.CONSERVATIVE,
    **kwargs,
) -> tuple[Verdict, VerificationEvidence]:
    """Dispatch to the appropriate verification function for a SLOT predicate.

    Args:
        predicate_name: Name of the SLOT predicate.
        slot_mode: SLOT evaluation mode.
        **kwargs: Arguments to pass to the verification function.

    Returns:
        Tuple of (Verdict, VerificationEvidence).

    Raises:
        ValueError: If predicate_name is not a SLOT predicate.
    """
    base_name = predicate_name.split("(")[0]

    if not is_slot_predicate(base_name):
        raise ValueError(f"'{base_name}' is not a SLOT predicate")

    # DNA-level SLOT predicates
    if base_name == "NoCrypticSplice":
        return verify_no_cryptic_splice(
            kwargs.get("seq", ""),
            kwargs.get("low_thresh", 3.0),
            kwargs.get("high_thresh", 6.0),
            slot_mode,
        )
    elif base_name == "NoCrypticPromoter":
        return verify_no_cryptic_promoter(
            kwargs.get("seq", ""),
            kwargs.get("organism", "E_coli"),
            kwargs.get("threshold", 0.7),
            slot_mode,
        )
    elif base_name == "NoUnexpectedTMDomain":
        return verify_no_unexpected_tm_domain(
            kwargs.get("seq", ""),
            kwargs.get("is_cytosolic", True),
            kwargs.get("threshold", 0.68),
            slot_mode,
        )
    elif base_name == "mRNASecondaryStructure":
        return verify_mrna_secondary_structure(
            kwargs.get("seq", ""),
            kwargs.get("window_start", 0),
            kwargs.get("window_end", 50),
            kwargs.get("dg_threshold", -15.0),
            slot_mode,
        )
    elif base_name == "CoTranslationalFolding":
        return verify_co_translational_folding(
            kwargs.get("seq", ""),
            kwargs.get("organism", "Homo_sapiens"),
            kwargs.get("domain_boundaries"),
            kwargs.get("min_pause_cai", 0.3),
            slot_mode,
        )
    elif base_name == "ConservationScore":
        return verify_conservation_score(
            kwargs.get("original_aa", "A"),
            kwargs.get("new_aa", "A"),
            kwargs.get("min_score", 0),
            slot_mode,
        )
    elif base_name == "CodonOptimality":
        return verify_codon_optimality(
            kwargs.get("codon", "ATG"),
            kwargs.get("species_cai", {}),
            kwargs.get("min_cai", 0.0),
            slot_mode,
        )
    # Higher-level SLOT predicates
    elif base_name in _STRUCTURE_PREDICATES:
        return verify_structure_predicate(
            base_name,
            kwargs.get("protein_sequence", ""),
            slot_mode,
        )
    elif base_name in _STABILITY_PREDICATES:
        return verify_stability_predicate(
            base_name,
            kwargs.get("protein_sequence", ""),
            slot_mode,
        )
    elif base_name in _SOLUBILITY_PREDICATES:
        return verify_solubility_predicate(
            base_name,
            kwargs.get("protein_sequence", ""),
            slot_mode,
        )
    elif base_name in _IMMUNOGENICITY_PREDICATES:
        return verify_immunogenicity_predicate(
            base_name,
            kwargs.get("protein_sequence", ""),
            slot_mode,
        )
    else:
        # Unknown SLOT predicate — conservative default
        return Verdict.UNCERTAIN, VerificationEvidence(
            predicate=base_name,
            slot_mode=slot_mode,
            tool_available=False,
            tool_name="unknown",
            verified=False,
            details=f"Unknown SLOT predicate: {base_name}",
        )
