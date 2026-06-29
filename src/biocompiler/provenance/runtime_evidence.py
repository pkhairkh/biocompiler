"""
Runtime evidence checks for SLOT predicate tool-soundness axioms.

Each check corresponds to a narrowed axiom in
``proof/BioCompiler/SLOTVerification.lean`` (the 34 narrowed axioms left
after the W1-A5 refactoring).  These checks do NOT prove the axioms — they
verify that the external tool's reported output is *self-consistent* and
*within the expected range* so that obvious tool malfunctions are caught
at runtime rather than silently trusted.

The 34 narrowed axioms are grouped by tool:

  TMHMM (3)
    - tmhmm_window_size_contract
    - tmhmm_cytosolic_only_contract
    - tmhmm_threshold_sound_contract

  ViennaRNA (2)
    - viennarna_window_size_contract
    - viennarna_deltaG_sound_contract

  AlphaFold co-translational (3)
    - alphafold_ramp_window_contract
    - alphafold_cotrans_threshold_contract
    - alphafold_adaptation_index_sound_contract

  FoldX stable folding (3)
    - foldx_stability_meaning_contract
    - foldx_estimated_deltaG_proxy_contract
    - foldx_stable_folding_sound_contract

  FoldX stability margin (2)
    - foldx_ddg_threshold_meaningful_contract
    - foldx_stability_margin_sound_contract

  FoldX destabilizing mutation (2)
    - foldx_max_ddg_meaningful_contract
    - foldx_destabilizing_mutation_sound_contract

  FoldX hydrophobic core (2)
    - foldx_core_window_contract
    - foldx_core_quality_sound_contract

  ProteinSol (3)
    - proteinsol_score_range_contract
    - proteinsol_gc_proxy_contract
    - proteinsol_solubility_sound_contract

  Aggrescan (3)
    - aggrescan_window_size_contract
    - aggrescan_threshold_value_contract
    - aggrescan_no_aggregation_sound_contract

  ExPASy (3)
    - expasy_pi_range_contract
    - expasy_gc_proxy_contract
    - expasy_charge_composition_sound_contract

  NetMHC (2)
    - netmhc_score_nonneg_contract
    - netmhc_threshold_nonneg_contract

  NetMHCpan (2)
    - netmhcpan_ic50_positive_contract
    - netmhcpan_threshold_positive_contract

  BepiPred (2)
    - bepipred_score_nonneg_contract
    - bepipred_threshold_nonneg_contract

  IEDB (2)
    - iedb_coverage_range_contract
    - iedb_threshold_range_contract

Reference: ``proof/BioCompiler/SLOTVerification.lean`` (lines 610-998).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Optional

# ─────────────────────────────────────────────────────────────────────────────
# Canonical constants — mirror the Lean definitions in SLOTVerification.lean
# ─────────────────────────────────────────────────────────────────────────────
TM_DOMAIN_WINDOW_SIZE = 51          # nucleotides (tmDomainWindowSize)
MRNA_STRUCTURE_WINDOW_SIZE = 30     # nucleotides (mrnaStructureWindowSize)
COTRANS_RAMP_CODONS = 30            # codons (cotransRampCodons)
TM_DOMAIN_THRESHOLD = 0.68          # 68/100 (tmDomainThreshold, Aggrescan canonical)
COTRANS_DISRUPTION_THRESHOLD_DEFAULT = 0.30   # cotransDisruptionThreshold default


@dataclass
class EvidenceCheckResult:
    """Result of a single runtime evidence check.

    Attributes:
        axiom_name: the name of the narrowed Lean axiom being checked.
        passed: True iff the tool's output is self-consistent and within
            the expected range.
        message: human-readable explanation of the outcome.
    """

    axiom_name: str
    passed: bool
    message: str

    def __bool__(self) -> bool:  # convenience: `if result:` works as expected
        return self.passed


# ─────────────────────────────────────────────────────────────────────────────
# TMHMM (NoUnexpectedTMDomain): 3 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_tmhmm_window_size_contract(window_size: int) -> EvidenceCheckResult:
    """Check: TMHMM's scanning window = 51 nt (``tmDomainWindowSize``).

    Lean axiom: ``tmhmm_window_size_contract``.
    """
    passed = isinstance(window_size, int) and window_size == TM_DOMAIN_WINDOW_SIZE
    if passed:
        msg = (
            f"TMHMM window size = {window_size} nt "
            f"(matches tmDomainWindowSize = {TM_DOMAIN_WINDOW_SIZE})"
        )
    else:
        msg = (
            f"TMHMM window size = {window_size!r}, "
            f"expected {TM_DOMAIN_WINDOW_SIZE} nt"
        )
    return EvidenceCheckResult("tmhmm_window_size_contract", passed, msg)


def check_tmhmm_cytosolic_only_contract(is_cytosolic: Any) -> EvidenceCheckResult:
    """Check: TMHMM only checks cytosolic proteins; vacuous otherwise.

    Lean axiom: ``tmhmm_cytosolic_only_contract``.
    The Lean axiom's conclusion is ``isCytosolic = false → True`` (always
    true).  The runtime check verifies that the input cytosolic flag is a
    well-formed boolean (not ``None`` / missing), so the tool's filter is
    being applied correctly.
    """
    passed = isinstance(is_cytosolic, bool)
    if passed:
        if is_cytosolic:
            msg = "TMHMM cytosolic flag = True (TM-domain check applies)"
        else:
            msg = "TMHMM cytosolic flag = False (check vacuously satisfied)"
    else:
        msg = (
            "TMHMM cytosolic flag must be a bool, "
            f"got {type(is_cytosolic).__name__}"
        )
    return EvidenceCheckResult("tmhmm_cytosolic_only_contract", passed, msg)


def check_tmhmm_threshold_sound_contract(
    max_hydrophobic_fraction: float,
    threshold: float,
    is_cytosolic: bool = True,
) -> EvidenceCheckResult:
    """Check: when cytosolic, max window hydrophobic fraction < threshold.

    Lean axiom: ``tmhmm_threshold_sound_contract``.
    """
    if not is_cytosolic:
        return EvidenceCheckResult(
            "tmhmm_threshold_sound_contract",
            True,
            "Vacuously satisfied: isCytosolic=False",
        )
    try:
        passed = float(max_hydrophobic_fraction) < float(threshold)
    except (TypeError, ValueError) as exc:
        passed = False
        msg = f"Cannot compare max_hydrophobic_fraction={max_hydrophobic_fraction!r} "
        f"and threshold={threshold!r}: {exc}"
        return EvidenceCheckResult("tmhmm_threshold_sound_contract", passed, msg)
    if passed:
        msg = (
            f"max hydrophobic fraction = {max_hydrophobic_fraction} < "
            f"threshold {threshold}"
        )
    else:
        msg = (
            f"max hydrophobic fraction = {max_hydrophobic_fraction} >= "
            f"threshold {threshold}"
        )
    return EvidenceCheckResult("tmhmm_threshold_sound_contract", passed, msg)


# ─────────────────────────────────────────────────────────────────────────────
# ViennaRNA (mRNASecondaryStructure): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_viennarna_window_size_contract(window_size: int) -> EvidenceCheckResult:
    """Check: ViennaRNA scanning window = 30 nt (``mrnaStructureWindowSize``).

    Lean axiom: ``viennarna_window_size_contract``.
    """
    passed = isinstance(window_size, int) and window_size == MRNA_STRUCTURE_WINDOW_SIZE
    if passed:
        msg = (
            f"ViennaRNA window size = {window_size} nt "
            f"(matches mrnaStructureWindowSize = {MRNA_STRUCTURE_WINDOW_SIZE})"
        )
    else:
        msg = (
            f"ViennaRNA window size = {window_size!r}, "
            f"expected {MRNA_STRUCTURE_WINDOW_SIZE} nt"
        )
    return EvidenceCheckResult("viennarna_window_size_contract", passed, msg)


def check_viennarna_deltaG_sound_contract(
    min_estimated_delta_g: float,
    dg_threshold: float,
) -> EvidenceCheckResult:
    """Check: ViennaRNA reported ΔG above ``dgThreshold`` ⇒ minimum window
    ``estimatedDeltaG`` exceeds ``dgThreshold``.

    Lean axiom: ``viennarna_deltaG_sound_contract``.
    """
    try:
        passed = float(min_estimated_delta_g) > float(dg_threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "viennarna_deltaG_sound_contract",
            False,
            f"Cannot compare min_estimated_delta_g={min_estimated_delta_g!r} "
            f"and dg_threshold={dg_threshold!r}: {exc}",
        )
    if passed:
        msg = (
            f"min estimated ΔG = {min_estimated_delta_g} > "
            f"dgThreshold {dg_threshold}"
        )
    else:
        msg = (
            f"min estimated ΔG = {min_estimated_delta_g} <= "
            f"dgThreshold {dg_threshold}"
        )
    return EvidenceCheckResult("viennarna_deltaG_sound_contract", passed, msg)


# ─────────────────────────────────────────────────────────────────────────────
# AlphaFold co-translational (CoTranslationalFolding): 3 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_alphafold_ramp_window_contract(ramp_codons: int) -> EvidenceCheckResult:
    """Check: AlphaFold ramp region = 30 codons (``cotransRampCodons``).

    Lean axiom: ``alphafold_ramp_window_contract``.
    """
    passed = isinstance(ramp_codons, int) and ramp_codons == COTRANS_RAMP_CODONS
    if passed:
        msg = (
            f"AlphaFold ramp region = {ramp_codons} codons "
            f"(matches cotransRampCodons = {COTRANS_RAMP_CODONS})"
        )
    else:
        msg = (
            f"AlphaFold ramp region = {ramp_codons!r} codons, "
            f"expected {COTRANS_RAMP_CODONS}"
        )
    return EvidenceCheckResult("alphafold_ramp_window_contract", passed, msg)


def check_alphafold_cotrans_threshold_contract(threshold: float) -> EvidenceCheckResult:
    """Check: ``cotransDisruptionThreshold`` ∈ [0, 1].

    Lean axiom: ``alphafold_cotrans_threshold_contract``.
    """
    try:
        t = float(threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "alphafold_cotrans_threshold_contract",
            False,
            f"cotransDisruptionThreshold not numeric: {threshold!r} ({exc})",
        )
    passed = 0.0 <= t <= 1.0
    msg = (
        f"cotransDisruptionThreshold = {t} in [0, 1]"
        if passed
        else f"cotransDisruptionThreshold = {t} outside [0, 1]"
    )
    return EvidenceCheckResult("alphafold_cotrans_threshold_contract", passed, msg)


def check_alphafold_adaptation_index_sound_contract(
    ramp_adaptation_index: float,
    threshold: float,
) -> EvidenceCheckResult:
    """Check: ``rampAdaptationIndex`` > ``cotransDisruptionThreshold``.

    Lean axiom: ``alphafold_adaptation_index_sound_contract``.
    """
    try:
        passed = float(ramp_adaptation_index) > float(threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "alphafold_adaptation_index_sound_contract",
            False,
            f"Cannot compare ramp_adaptation_index={ramp_adaptation_index!r} "
            f"and threshold={threshold!r}: {exc}",
        )
    if passed:
        msg = (
            f"rampAdaptationIndex = {ramp_adaptation_index} > "
            f"threshold {threshold}"
        )
    else:
        msg = (
            f"rampAdaptationIndex = {ramp_adaptation_index} <= "
            f"threshold {threshold}"
        )
    return EvidenceCheckResult(
        "alphafold_adaptation_index_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# FoldX stable folding (NoMisfoldingRisk): 3 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_foldx_stability_meaning_contract(stable_threshold: Optional[float] = None) -> EvidenceCheckResult:
    """Check: FoldX "stableFolding" criterion is non-trivial.

    Lean axiom: ``foldx_stability_meaning_contract`` — the Lean statement
    is the tautology ``(0 : Rat) < (1 : Rat)`` (the criterion is
    non-trivial).  The runtime check verifies that the threshold value
    used by the tool is finite and well-defined (a meaningful stability
    criterion must have a numeric threshold).

    Passing ``None`` (the default) is interpreted as "the criterion is
    non-trivial by construction" (mirrors the Lean tautology).
    """
    if stable_threshold is None:
        return EvidenceCheckResult(
            "foldx_stability_meaning_contract",
            True,
            "FoldX stableFolding criterion is non-trivial (0 < 1)",
        )
    passed = isinstance(stable_threshold, (int, float)) and math.isfinite(float(stable_threshold))
    msg = (
        f"FoldX stability threshold = {stable_threshold} (finite, meaningful)"
        if passed
        else f"FoldX stability threshold = {stable_threshold!r} (non-finite or non-numeric)"
    )
    return EvidenceCheckResult("foldx_stability_meaning_contract", passed, msg)


def check_foldx_estimated_deltaG_proxy_contract(estimated_delta_g: float) -> EvidenceCheckResult:
    """Check: ``estimatedDeltaG`` is a valid (finite) proxy for FoldX's ΔG.

    Lean axiom: ``foldx_estimated_deltaG_proxy_contract`` — the Lean
    statement is tautological (a value's sign matches itself).  The
    runtime check verifies that the proxy is a finite number so the sign
    convention (negative = stable) is well-defined.
    """
    passed = (
        isinstance(estimated_delta_g, (int, float))
        and not isinstance(estimated_delta_g, bool)
        and math.isfinite(float(estimated_delta_g))
    )
    msg = (
        f"estimatedDeltaG = {estimated_delta_g} (finite, self-consistent sign)"
        if passed
        else f"estimatedDeltaG = {estimated_delta_g!r} (non-finite or non-numeric, proxy broken)"
    )
    return EvidenceCheckResult(
        "foldx_estimated_deltaG_proxy_contract", passed, msg
    )


def check_foldx_stable_folding_sound_contract(estimated_delta_g: float) -> EvidenceCheckResult:
    """Check: FoldX stable-folding report ⇒ ``estimatedDeltaG`` < 0.

    Lean axiom: ``foldx_stable_folding_sound_contract``.
    """
    try:
        passed = float(estimated_delta_g) < 0.0
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "foldx_stable_folding_sound_contract",
            False,
            f"estimatedDeltaG not numeric: {estimated_delta_g!r} ({exc})",
        )
    msg = (
        f"estimatedDeltaG = {estimated_delta_g} < 0 (stable folding confirmed)"
        if passed
        else f"estimatedDeltaG = {estimated_delta_g} >= 0 (not stable)"
    )
    return EvidenceCheckResult(
        "foldx_stable_folding_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# FoldX stability margin (StableFolding): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_foldx_ddg_threshold_meaningful_contract(ddg_threshold: float) -> EvidenceCheckResult:
    """Check: ``ddgThreshold`` ≥ 0 (meaningful stability margin).

    Lean axiom: ``foldx_ddg_threshold_meaningful_contract``.
    """
    try:
        passed = float(ddg_threshold) >= 0.0
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "foldx_ddg_threshold_meaningful_contract",
            False,
            f"ddgThreshold not numeric: {ddg_threshold!r} ({exc})",
        )
    msg = (
        f"ddgThreshold = {ddg_threshold} >= 0 (meaningful stability margin)"
        if passed
        else f"ddgThreshold = {ddg_threshold} < 0 (not meaningful)"
    )
    return EvidenceCheckResult(
        "foldx_ddg_threshold_meaningful_contract", passed, msg
    )


def check_foldx_stability_margin_sound_contract(
    estimated_delta_g: float,
    ddg_threshold: float,
) -> EvidenceCheckResult:
    """Check: ``estimatedDeltaG`` ≤ -``ddgThreshold`` (stability margin).

    Lean axiom: ``foldx_stability_margin_sound_contract``.
    """
    try:
        eg = float(estimated_delta_g)
        dt = float(ddg_threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "foldx_stability_margin_sound_contract",
            False,
            f"Cannot evaluate stability margin: {exc}",
        )
    passed = eg <= -dt
    msg = (
        f"estimatedDeltaG = {eg} <= -ddgThreshold = {-dt}"
        if passed
        else f"estimatedDeltaG = {eg} > -ddgThreshold = {-dt}"
    )
    return EvidenceCheckResult(
        "foldx_stability_margin_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# FoldX destabilizing mutation (NoDestabilizingMutation): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_foldx_max_ddg_meaningful_contract(max_ddg: float) -> EvidenceCheckResult:
    """Check: ``maxDDG`` ≥ 0 (meaningful destabilization bound).

    Lean axiom: ``foldx_max_ddg_meaningful_contract``.
    """
    try:
        passed = float(max_ddg) >= 0.0
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "foldx_max_ddg_meaningful_contract",
            False,
            f"maxDDG not numeric: {max_ddg!r} ({exc})",
        )
    msg = (
        f"maxDDG = {max_ddg} >= 0 (meaningful destabilization bound)"
        if passed
        else f"maxDDG = {max_ddg} < 0 (not meaningful)"
    )
    return EvidenceCheckResult(
        "foldx_max_ddg_meaningful_contract", passed, msg
    )


def check_foldx_destabilizing_mutation_sound_contract(
    estimated_delta_g: float,
    max_ddg: float,
) -> EvidenceCheckResult:
    """Check: ``estimatedDeltaG`` ≤ -``maxDDG`` (no destabilizing mutation).

    Lean axiom: ``foldx_destabilizing_mutation_sound_contract``.
    """
    try:
        eg = float(estimated_delta_g)
        md = float(max_ddg)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "foldx_destabilizing_mutation_sound_contract",
            False,
            f"Cannot evaluate destabilization bound: {exc}",
        )
    passed = eg <= -md
    msg = (
        f"estimatedDeltaG = {eg} <= -maxDDG = {-md}"
        if passed
        else f"estimatedDeltaG = {eg} > -maxDDG = {-md}"
    )
    return EvidenceCheckResult(
        "foldx_destabilizing_mutation_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# FoldX hydrophobic core (HydrophobicCoreQuality): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_foldx_core_window_contract(window_size: int) -> EvidenceCheckResult:
    """Check: FoldX core window = 51 nt (``tmDomainWindowSize``).

    Lean axiom: ``foldx_core_window_contract``.
    """
    passed = isinstance(window_size, int) and window_size == TM_DOMAIN_WINDOW_SIZE
    if passed:
        msg = (
            f"FoldX core window = {window_size} nt "
            f"(matches tmDomainWindowSize = {TM_DOMAIN_WINDOW_SIZE})"
        )
    else:
        msg = (
            f"FoldX core window = {window_size!r} nt, "
            f"expected {TM_DOMAIN_WINDOW_SIZE}"
        )
    return EvidenceCheckResult("foldx_core_window_contract", passed, msg)


def check_foldx_core_quality_sound_contract(
    max_hydrophobic_fraction: float,
    threshold: float,
) -> EvidenceCheckResult:
    """Check: FoldX core quality report ⇒ max hydrophobic fraction ≥ threshold.

    Lean axiom: ``foldx_core_quality_sound_contract``.
    """
    try:
        passed = float(max_hydrophobic_fraction) >= float(threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "foldx_core_quality_sound_contract",
            False,
            f"Cannot compare max_hydrophobic_fraction={max_hydrophobic_fraction!r} "
            f"and threshold={threshold!r}: {exc}",
        )
    msg = (
        f"max hydrophobic fraction = {max_hydrophobic_fraction} >= "
        f"threshold {threshold}"
        if passed
        else f"max hydrophobic fraction = {max_hydrophobic_fraction} < "
        f"threshold {threshold}"
    )
    return EvidenceCheckResult(
        "foldx_core_quality_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# ProteinSol (SolubleExpression): 3 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_proteinsol_score_range_contract(score: float) -> EvidenceCheckResult:
    """Check: ProteinSol solubility score ∈ [0, 1].

    Lean axiom: ``proteinsol_score_range_contract``.
    """
    try:
        s = float(score)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "proteinsol_score_range_contract",
            False,
            f"ProteinSol score not numeric: {score!r} ({exc})",
        )
    passed = 0.0 <= s <= 1.0
    msg = (
        f"ProteinSol score = {s} in [0, 1]"
        if passed
        else f"ProteinSol score = {s} outside [0, 1]"
    )
    return EvidenceCheckResult("proteinsol_score_range_contract", passed, msg)


def check_proteinsol_gc_proxy_contract(gc_content: float) -> EvidenceCheckResult:
    """Check: GC content ∈ [0, 1] (valid ProteinSol proxy).

    Lean axiom: ``proteinsol_gc_proxy_contract``.
    """
    try:
        g = float(gc_content)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "proteinsol_gc_proxy_contract",
            False,
            f"GC content not numeric: {gc_content!r} ({exc})",
        )
    passed = 0.0 <= g <= 1.0
    msg = (
        f"GC content = {g} in [0, 1] (valid ProteinSol proxy)"
        if passed
        else f"GC content = {g} outside [0, 1] (invalid proxy)"
    )
    return EvidenceCheckResult("proteinsol_gc_proxy_contract", passed, msg)


def check_proteinsol_solubility_sound_contract(
    gc_content: float,
    min_score: float,
) -> EvidenceCheckResult:
    """Check: ProteinSol reported solubility ≥ ``minScore`` ⇒ ``gcContent``
    ≥ ``minScore`` (in-repo proxy meets threshold).

    Lean axiom: ``proteinsol_solubility_sound_contract``.
    """
    try:
        g = float(gc_content)
        m = float(min_score)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "proteinsol_solubility_sound_contract",
            False,
            f"Cannot evaluate ProteinSol solubility: {exc}",
        )
    passed = g >= m
    msg = (
        f"gcContent = {g} >= minScore {m}"
        if passed
        else f"gcContent = {g} < minScore {m}"
    )
    return EvidenceCheckResult(
        "proteinsol_solubility_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# Aggrescan (NoAggregationProneRegion): 3 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_aggrescan_window_size_contract(window_size: int) -> EvidenceCheckResult:
    """Check: Aggrescan window = 51 nt (``tmDomainWindowSize``).

    Lean axiom: ``aggrescan_window_size_contract``.
    """
    passed = isinstance(window_size, int) and window_size == TM_DOMAIN_WINDOW_SIZE
    if passed:
        msg = (
            f"Aggrescan window = {window_size} nt "
            f"(matches tmDomainWindowSize = {TM_DOMAIN_WINDOW_SIZE})"
        )
    else:
        msg = (
            f"Aggrescan window = {window_size!r} nt, "
            f"expected {TM_DOMAIN_WINDOW_SIZE}"
        )
    return EvidenceCheckResult("aggrescan_window_size_contract", passed, msg)


def check_aggrescan_threshold_value_contract(threshold: float) -> EvidenceCheckResult:
    """Check: ``tmDomainThreshold`` = 0.68 (canonical Aggrescan threshold).

    Lean axiom: ``aggrescan_threshold_value_contract``.
    """
    try:
        t = float(threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "aggrescan_threshold_value_contract",
            False,
            f"tmDomainThreshold not numeric: {threshold!r} ({exc})",
        )
    passed = abs(t - TM_DOMAIN_THRESHOLD) < 1e-9
    msg = (
        f"tmDomainThreshold = {t} (matches canonical value {TM_DOMAIN_THRESHOLD})"
        if passed
        else f"tmDomainThreshold = {t}, expected {TM_DOMAIN_THRESHOLD}"
    )
    return EvidenceCheckResult("aggrescan_threshold_value_contract", passed, msg)


def check_aggrescan_no_aggregation_sound_contract(
    max_hydrophobic_fraction: float,
) -> EvidenceCheckResult:
    """Check: no window reaches hydrophobic fraction ≥ ``tmDomainThreshold``.

    Lean axiom: ``aggrescan_no_aggregation_sound_contract``.
    """
    try:
        m = float(max_hydrophobic_fraction)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "aggrescan_no_aggregation_sound_contract",
            False,
            f"max_hydrophobic_fraction not numeric: "
            f"{max_hydrophobic_fraction!r} ({exc})",
        )
    passed = m < TM_DOMAIN_THRESHOLD
    msg = (
        f"max hydrophobic fraction = {m} < "
        f"tmDomainThreshold {TM_DOMAIN_THRESHOLD}"
        if passed
        else f"max hydrophobic fraction = {m} >= "
        f"tmDomainThreshold {TM_DOMAIN_THRESHOLD}"
    )
    return EvidenceCheckResult(
        "aggrescan_no_aggregation_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# ExPASy (ChargeComposition): 3 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_expasy_pi_range_contract(pi: float) -> EvidenceCheckResult:
    """Check: ExPASy computed pI ∈ [0, 14].

    Lean axiom: ``expasy_pi_range_contract``.
    """
    try:
        p = float(pi)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "expasy_pi_range_contract",
            False,
            f"ExPASy pI not numeric: {pi!r} ({exc})",
        )
    passed = 0.0 <= p <= 14.0
    msg = (
        f"ExPASy pI = {p} in [0, 14]"
        if passed
        else f"ExPASy pI = {p} outside [0, 14]"
    )
    return EvidenceCheckResult("expasy_pi_range_contract", passed, msg)


def check_expasy_gc_proxy_contract(gc_content: float) -> EvidenceCheckResult:
    """Check: GC content ∈ [0, 1] (valid ExPASy pI proxy).

    Lean axiom: ``expasy_gc_proxy_contract``.
    """
    try:
        g = float(gc_content)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "expasy_gc_proxy_contract",
            False,
            f"GC content not numeric: {gc_content!r} ({exc})",
        )
    passed = 0.0 <= g <= 1.0
    msg = (
        f"GC content = {g} in [0, 1] (valid ExPASy proxy)"
        if passed
        else f"GC content = {g} outside [0, 1] (invalid proxy)"
    )
    return EvidenceCheckResult("expasy_gc_proxy_contract", passed, msg)


def check_expasy_charge_composition_sound_contract(
    gc_content: float,
    pi_lo: float,
    pi_hi: float,
) -> EvidenceCheckResult:
    """Check: ExPASy pI ∈ [pILo, pIHi] ⇒ ``gcContent`` ∈ [pILo, pIHi].

    Lean axiom: ``expasy_charge_composition_sound_contract``.
    """
    try:
        g = float(gc_content)
        lo = float(pi_lo)
        hi = float(pi_hi)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "expasy_charge_composition_sound_contract",
            False,
            f"Cannot evaluate ExPASy charge composition: {exc}",
        )
    if lo > hi:
        return EvidenceCheckResult(
            "expasy_charge_composition_sound_contract",
            False,
            f"Invalid range: pILo={lo} > pIHi={hi}",
        )
    passed = lo <= g <= hi
    msg = (
        f"gcContent = {g} in [pILo={lo}, pIHi={hi}]"
        if passed
        else f"gcContent = {g} outside [pILo={lo}, pIHi={hi}]"
    )
    return EvidenceCheckResult(
        "expasy_charge_composition_sound_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# NetMHC (LowImmunogenicity): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_netmhc_score_nonneg_contract(score: float) -> EvidenceCheckResult:
    """Check: NetMHC immunogenicity score ≥ 0.

    Lean axiom: ``netmhc_score_nonneg_contract``.
    """
    try:
        s = float(score)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "netmhc_score_nonneg_contract",
            False,
            f"NetMHC score not numeric: {score!r} ({exc})",
        )
    passed = s >= 0.0
    msg = (
        f"NetMHC score = {s} >= 0"
        if passed
        else f"NetMHC score = {s} < 0"
    )
    return EvidenceCheckResult("netmhc_score_nonneg_contract", passed, msg)


def check_netmhc_threshold_nonneg_contract(max_score: float) -> EvidenceCheckResult:
    """Check: ``maxScore`` ≥ 0 (meaningful immunogenicity bound).

    Lean axiom: ``netmhc_threshold_nonneg_contract``.
    """
    try:
        m = float(max_score)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "netmhc_threshold_nonneg_contract",
            False,
            f"maxScore not numeric: {max_score!r} ({exc})",
        )
    passed = m >= 0.0
    msg = (
        f"maxScore = {m} >= 0 (meaningful immunogenicity bound)"
        if passed
        else f"maxScore = {m} < 0 (not meaningful)"
    )
    return EvidenceCheckResult("netmhc_threshold_nonneg_contract", passed, msg)


# ─────────────────────────────────────────────────────────────────────────────
# NetMHCpan (NoStrongTCellEpitope): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_netmhcpan_ic50_positive_contract(ic50: float) -> EvidenceCheckResult:
    """Check: NetMHCpan binding IC50 > 0 (a concentration).

    Lean axiom: ``netmhcpan_ic50_positive_contract``.
    """
    try:
        v = float(ic50)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "netmhcpan_ic50_positive_contract",
            False,
            f"NetMHCpan IC50 not numeric: {ic50!r} ({exc})",
        )
    passed = v > 0.0
    msg = (
        f"NetMHCpan IC50 = {v} > 0"
        if passed
        else f"NetMHCpan IC50 = {v} <= 0"
    )
    return EvidenceCheckResult("netmhcpan_ic50_positive_contract", passed, msg)


def check_netmhcpan_threshold_positive_contract(ic50_threshold: float) -> EvidenceCheckResult:
    """Check: ``ic50Threshold`` > 0 (meaningful binding threshold).

    Lean axiom: ``netmhcpan_threshold_positive_contract``.
    """
    try:
        t = float(ic50_threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "netmhcpan_threshold_positive_contract",
            False,
            f"ic50Threshold not numeric: {ic50_threshold!r} ({exc})",
        )
    passed = t > 0.0
    msg = (
        f"ic50Threshold = {t} > 0 (meaningful binding threshold)"
        if passed
        else f"ic50Threshold = {t} <= 0 (not meaningful)"
    )
    return EvidenceCheckResult(
        "netmhcpan_threshold_positive_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# BepiPred (NoDominantBCellEpitope): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_bepipred_score_nonneg_contract(score: float) -> EvidenceCheckResult:
    """Check: BepiPred B-cell epitope score ≥ 0.

    Lean axiom: ``bepipred_score_nonneg_contract``.
    """
    try:
        s = float(score)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "bepipred_score_nonneg_contract",
            False,
            f"BepiPred score not numeric: {score!r} ({exc})",
        )
    passed = s >= 0.0
    msg = (
        f"BepiPred score = {s} >= 0"
        if passed
        else f"BepiPred score = {s} < 0"
    )
    return EvidenceCheckResult("bepipred_score_nonneg_contract", passed, msg)


def check_bepipred_threshold_nonneg_contract(score_threshold: float) -> EvidenceCheckResult:
    """Check: ``scoreThreshold`` ≥ 0.

    Lean axiom: ``bepipred_threshold_nonneg_contract``.
    """
    try:
        t = float(score_threshold)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "bepipred_threshold_nonneg_contract",
            False,
            f"scoreThreshold not numeric: {score_threshold!r} ({exc})",
        )
    passed = t >= 0.0
    msg = (
        f"scoreThreshold = {t} >= 0"
        if passed
        else f"scoreThreshold = {t} < 0"
    )
    return EvidenceCheckResult(
        "bepipred_threshold_nonneg_contract", passed, msg
    )


# ─────────────────────────────────────────────────────────────────────────────
# IEDB (PopulationCoverageSafe): 2 narrowed axioms
# ─────────────────────────────────────────────────────────────────────────────

def check_iedb_coverage_range_contract(coverage: float) -> EvidenceCheckResult:
    """Check: IEDB population-coverage ∈ [0, 1].

    Lean axiom: ``iedb_coverage_range_contract``.
    """
    try:
        c = float(coverage)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "iedb_coverage_range_contract",
            False,
            f"IEDB coverage not numeric: {coverage!r} ({exc})",
        )
    passed = 0.0 <= c <= 1.0
    msg = (
        f"IEDB coverage = {c} in [0, 1]"
        if passed
        else f"IEDB coverage = {c} outside [0, 1]"
    )
    return EvidenceCheckResult("iedb_coverage_range_contract", passed, msg)


def check_iedb_threshold_range_contract(max_coverage: float) -> EvidenceCheckResult:
    """Check: ``maxCoverage`` ∈ [0, 1] (meaningful coverage bound).

    Lean axiom: ``iedb_threshold_range_contract``.
    """
    try:
        m = float(max_coverage)
    except (TypeError, ValueError) as exc:
        return EvidenceCheckResult(
            "iedb_threshold_range_contract",
            False,
            f"maxCoverage not numeric: {max_coverage!r} ({exc})",
        )
    passed = 0.0 <= m <= 1.0
    msg = (
        f"maxCoverage = {m} in [0, 1]"
        if passed
        else f"maxCoverage = {m} outside [0, 1]"
    )
    return EvidenceCheckResult("iedb_threshold_range_contract", passed, msg)


# ─────────────────────────────────────────────────────────────────────────────
# Registry: maps tool-output dict keys → check functions
# ─────────────────────────────────────────────────────────────────────────────

# Mapping from a tool-output dict key to a callable that takes the tool_outputs
# dict and returns an EvidenceCheckResult.  Used by ``run_all_evidence_checks``.
def _evidence_check_registry() -> list:
    """Return a list of (predicate, callable) pairs.

    Each predicate takes the ``tool_outputs`` dict and returns ``True`` if
    the check is applicable (i.e., the dict has the keys needed by the
    callable).  The callable takes the dict and returns an
    ``EvidenceCheckResult``.
    """
    return [
        # ── TMHMM (3) ───────────────────────────────────────────────────────
        (
            lambda d: "tmhmm_window_size" in d,
            lambda d: check_tmhmm_window_size_contract(d["tmhmm_window_size"]),
        ),
        (
            lambda d: "tmhmm_is_cytosolic" in d,
            lambda d: check_tmhmm_cytosolic_only_contract(d["tmhmm_is_cytosolic"]),
        ),
        (
            lambda d: "tmhmm_max_hydrophobic_fraction" in d,
            lambda d: check_tmhmm_threshold_sound_contract(
                d["tmhmm_max_hydrophobic_fraction"],
                d.get("tmhmm_threshold", TM_DOMAIN_THRESHOLD),
                d.get("tmhmm_is_cytosolic", True),
            ),
        ),
        # ── ViennaRNA (2) ───────────────────────────────────────────────────
        (
            lambda d: "viennarna_window_size" in d,
            lambda d: check_viennarna_window_size_contract(d["viennarna_window_size"]),
        ),
        (
            lambda d: "viennarna_min_estimated_delta_g" in d,
            lambda d: check_viennarna_deltaG_sound_contract(
                d["viennarna_min_estimated_delta_g"],
                d.get("viennarna_dg_threshold", -10.0),
            ),
        ),
        # ── AlphaFold (3) ───────────────────────────────────────────────────
        (
            lambda d: "alphafold_ramp_codons" in d,
            lambda d: check_alphafold_ramp_window_contract(d["alphafold_ramp_codons"]),
        ),
        (
            lambda d: "alphafold_cotrans_threshold" in d,
            lambda d: check_alphafold_cotrans_threshold_contract(
                d["alphafold_cotrans_threshold"]
            ),
        ),
        (
            lambda d: "alphafold_ramp_adaptation_index" in d,
            lambda d: check_alphafold_adaptation_index_sound_contract(
                d["alphafold_ramp_adaptation_index"],
                d.get("alphafold_cotrans_threshold", COTRANS_DISRUPTION_THRESHOLD_DEFAULT),
            ),
        ),
        # ── FoldX stable folding (3) ────────────────────────────────────────
        (
            lambda d: "foldx_stability_criterion_present" in d,
            lambda d: check_foldx_stability_meaning_contract(
                d.get("foldx_stability_threshold")
            ),
        ),
        (
            lambda d: "foldx_estimated_delta_g" in d,
            lambda d: check_foldx_estimated_deltaG_proxy_contract(
                d["foldx_estimated_delta_g"]
            ),
        ),
        (
            lambda d: "foldx_stable_folding_passed" in d and d["foldx_stable_folding_passed"],
            lambda d: check_foldx_stable_folding_sound_contract(
                d.get("foldx_estimated_delta_g", 0.0)
            ),
        ),
        # ── FoldX stability margin (2) ──────────────────────────────────────
        (
            lambda d: "foldx_ddg_threshold" in d,
            lambda d: check_foldx_ddg_threshold_meaningful_contract(
                d["foldx_ddg_threshold"]
            ),
        ),
        (
            lambda d: "foldx_estimated_delta_g" in d and "foldx_ddg_threshold" in d,
            lambda d: check_foldx_stability_margin_sound_contract(
                d["foldx_estimated_delta_g"], d["foldx_ddg_threshold"]
            ),
        ),
        # ── FoldX destabilizing mutation (2) ────────────────────────────────
        (
            lambda d: "foldx_max_ddg" in d,
            lambda d: check_foldx_max_ddg_meaningful_contract(d["foldx_max_ddg"]),
        ),
        (
            lambda d: "foldx_estimated_delta_g" in d and "foldx_max_ddg" in d,
            lambda d: check_foldx_destabilizing_mutation_sound_contract(
                d["foldx_estimated_delta_g"], d["foldx_max_ddg"]
            ),
        ),
        # ── FoldX hydrophobic core (2) ──────────────────────────────────────
        (
            lambda d: "foldx_core_window_size" in d,
            lambda d: check_foldx_core_window_contract(d["foldx_core_window_size"]),
        ),
        (
            lambda d: "foldx_core_max_hydrophobic_fraction" in d,
            lambda d: check_foldx_core_quality_sound_contract(
                d["foldx_core_max_hydrophobic_fraction"],
                d.get("foldx_core_threshold", TM_DOMAIN_THRESHOLD),
            ),
        ),
        # ── ProteinSol (3) ──────────────────────────────────────────────────
        (
            lambda d: "proteinsol_score" in d,
            lambda d: check_proteinsol_score_range_contract(d["proteinsol_score"]),
        ),
        (
            lambda d: "proteinsol_gc_content" in d,
            lambda d: check_proteinsol_gc_proxy_contract(d["proteinsol_gc_content"]),
        ),
        (
            lambda d: "proteinsol_gc_content" in d and "proteinsol_min_score" in d,
            lambda d: check_proteinsol_solubility_sound_contract(
                d["proteinsol_gc_content"], d["proteinsol_min_score"]
            ),
        ),
        # ── Aggrescan (3) ───────────────────────────────────────────────────
        (
            lambda d: "aggrescan_window_size" in d,
            lambda d: check_aggrescan_window_size_contract(d["aggrescan_window_size"]),
        ),
        (
            lambda d: "aggrescan_threshold" in d,
            lambda d: check_aggrescan_threshold_value_contract(d["aggrescan_threshold"]),
        ),
        (
            lambda d: "aggrescan_max_hydrophobic_fraction" in d,
            lambda d: check_aggrescan_no_aggregation_sound_contract(
                d["aggrescan_max_hydrophobic_fraction"]
            ),
        ),
        # ── ExPASy (3) ──────────────────────────────────────────────────────
        (
            lambda d: "expasy_pi" in d,
            lambda d: check_expasy_pi_range_contract(d["expasy_pi"]),
        ),
        (
            lambda d: "expasy_gc_content" in d,
            lambda d: check_expasy_gc_proxy_contract(d["expasy_gc_content"]),
        ),
        (
            lambda d: "expasy_gc_content" in d
            and "expasy_pi_lo" in d
            and "expasy_pi_hi" in d,
            lambda d: check_expasy_charge_composition_sound_contract(
                d["expasy_gc_content"], d["expasy_pi_lo"], d["expasy_pi_hi"]
            ),
        ),
        # ── NetMHC (2) ──────────────────────────────────────────────────────
        (
            lambda d: "netmhc_score" in d,
            lambda d: check_netmhc_score_nonneg_contract(d["netmhc_score"]),
        ),
        (
            lambda d: "netmhc_max_score" in d,
            lambda d: check_netmhc_threshold_nonneg_contract(d["netmhc_max_score"]),
        ),
        # ── NetMHCpan (2) ───────────────────────────────────────────────────
        (
            lambda d: "netmhcpan_ic50" in d,
            lambda d: check_netmhcpan_ic50_positive_contract(d["netmhcpan_ic50"]),
        ),
        (
            lambda d: "netmhcpan_ic50_threshold" in d,
            lambda d: check_netmhcpan_threshold_positive_contract(
                d["netmhcpan_ic50_threshold"]
            ),
        ),
        # ── BepiPred (2) ────────────────────────────────────────────────────
        (
            lambda d: "bepipred_score" in d,
            lambda d: check_bepipred_score_nonneg_contract(d["bepipred_score"]),
        ),
        (
            lambda d: "bepipred_score_threshold" in d,
            lambda d: check_bepipred_threshold_nonneg_contract(
                d["bepipred_score_threshold"]
            ),
        ),
        # ── IEDB (2) ────────────────────────────────────────────────────────
        (
            lambda d: "iedb_coverage" in d,
            lambda d: check_iedb_coverage_range_contract(d["iedb_coverage"]),
        ),
        (
            lambda d: "iedb_max_coverage" in d,
            lambda d: check_iedb_threshold_range_contract(d["iedb_max_coverage"]),
        ),
    ]


def run_all_evidence_checks(tool_outputs: dict) -> list[EvidenceCheckResult]:
    """Run all applicable evidence checks on a dict of tool outputs.

    The dict keys are tool-output identifiers (e.g. ``"tmhmm_window_size"``,
    ``"viennarna_min_estimated_delta_g"``, ``"netmhcpan_ic50"``,
    ``"foldx_estimated_delta_g"`` …).  Missing keys result in skipped
    checks (no entry in the returned list for that axiom).

    Returns:
        A list of :class:`EvidenceCheckResult` instances — one per
        applicable check.  Each result's ``axiom_name`` matches a
        narrowed Lean axiom in ``SLOTVerification.lean``.
    """
    if not isinstance(tool_outputs, dict):
        return [
            EvidenceCheckResult(
                "run_all_evidence_checks",
                False,
                f"tool_outputs must be a dict, got {type(tool_outputs).__name__}",
            )
        ]

    results: list[EvidenceCheckResult] = []
    for predicate, run in _evidence_check_registry():
        try:
            if predicate(tool_outputs):
                results.append(run(tool_outputs))
        except Exception as exc:  # pragma: no cover — defensive
            results.append(
                EvidenceCheckResult(
                    "run_all_evidence_checks",
                    False,
                    f"Check raised unexpected exception: {type(exc).__name__}: {exc}",
                )
            )
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: the full list of 34 axiom names (one per narrowed axiom)
# ─────────────────────────────────────────────────────────────────────────────
ALL_AXIOM_NAMES: tuple[str, ...] = (
    "tmhmm_window_size_contract",
    "tmhmm_cytosolic_only_contract",
    "tmhmm_threshold_sound_contract",
    "viennarna_window_size_contract",
    "viennarna_deltaG_sound_contract",
    "alphafold_ramp_window_contract",
    "alphafold_cotrans_threshold_contract",
    "alphafold_adaptation_index_sound_contract",
    "foldx_stability_meaning_contract",
    "foldx_estimated_deltaG_proxy_contract",
    "foldx_stable_folding_sound_contract",
    "foldx_ddg_threshold_meaningful_contract",
    "foldx_stability_margin_sound_contract",
    "foldx_max_ddg_meaningful_contract",
    "foldx_destabilizing_mutation_sound_contract",
    "foldx_core_window_contract",
    "foldx_core_quality_sound_contract",
    "proteinsol_score_range_contract",
    "proteinsol_gc_proxy_contract",
    "proteinsol_solubility_sound_contract",
    "aggrescan_window_size_contract",
    "aggrescan_threshold_value_contract",
    "aggrescan_no_aggregation_sound_contract",
    "expasy_pi_range_contract",
    "expasy_gc_proxy_contract",
    "expasy_charge_composition_sound_contract",
    "netmhc_score_nonneg_contract",
    "netmhc_threshold_nonneg_contract",
    "netmhcpan_ic50_positive_contract",
    "netmhcpan_threshold_positive_contract",
    "bepipred_score_nonneg_contract",
    "bepipred_threshold_nonneg_contract",
    "iedb_coverage_range_contract",
    "iedb_threshold_range_contract",
)

assert len(ALL_AXIOM_NAMES) == 34, (
    f"Expected 34 narrowed axioms, got {len(ALL_AXIOM_NAMES)}"
)


__all__ = [
    # Data class
    "EvidenceCheckResult",
    # Constants
    "TM_DOMAIN_WINDOW_SIZE",
    "MRNA_STRUCTURE_WINDOW_SIZE",
    "COTRANS_RAMP_CODONS",
    "TM_DOMAIN_THRESHOLD",
    "COTRANS_DISRUPTION_THRESHOLD_DEFAULT",
    "ALL_AXIOM_NAMES",
    # Aggregator
    "run_all_evidence_checks",
    # TMHMM checks (3)
    "check_tmhmm_window_size_contract",
    "check_tmhmm_cytosolic_only_contract",
    "check_tmhmm_threshold_sound_contract",
    # ViennaRNA checks (2)
    "check_viennarna_window_size_contract",
    "check_viennarna_deltaG_sound_contract",
    # AlphaFold checks (3)
    "check_alphafold_ramp_window_contract",
    "check_alphafold_cotrans_threshold_contract",
    "check_alphafold_adaptation_index_sound_contract",
    # FoldX stable folding checks (3)
    "check_foldx_stability_meaning_contract",
    "check_foldx_estimated_deltaG_proxy_contract",
    "check_foldx_stable_folding_sound_contract",
    # FoldX stability margin checks (2)
    "check_foldx_ddg_threshold_meaningful_contract",
    "check_foldx_stability_margin_sound_contract",
    # FoldX destabilizing mutation checks (2)
    "check_foldx_max_ddg_meaningful_contract",
    "check_foldx_destabilizing_mutation_sound_contract",
    # FoldX hydrophobic core checks (2)
    "check_foldx_core_window_contract",
    "check_foldx_core_quality_sound_contract",
    # ProteinSol checks (3)
    "check_proteinsol_score_range_contract",
    "check_proteinsol_gc_proxy_contract",
    "check_proteinsol_solubility_sound_contract",
    # Aggrescan checks (3)
    "check_aggrescan_window_size_contract",
    "check_aggrescan_threshold_value_contract",
    "check_aggrescan_no_aggregation_sound_contract",
    # ExPASy checks (3)
    "check_expasy_pi_range_contract",
    "check_expasy_gc_proxy_contract",
    "check_expasy_charge_composition_sound_contract",
    # NetMHC checks (2)
    "check_netmhc_score_nonneg_contract",
    "check_netmhc_threshold_nonneg_contract",
    # NetMHCpan checks (2)
    "check_netmhcpan_ic50_positive_contract",
    "check_netmhcpan_threshold_positive_contract",
    # BepiPred checks (2)
    "check_bepipred_score_nonneg_contract",
    "check_bepipred_threshold_nonneg_contract",
    # IEDB checks (2)
    "check_iedb_coverage_range_contract",
    "check_iedb_threshold_range_contract",
]
