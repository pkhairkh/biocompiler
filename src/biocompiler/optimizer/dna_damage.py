"""
BioCompiler DNA Damage Prediction Module v1.0.0
=================================================
Comprehensive DNA damage prediction, methylation modeling, and repair
pathway analysis for therapeutic gene design.

SOTA Upgrades (Phase 0+3):
  - DNA Repair Pathway Modeling: Literature-derived BER/NER/MMR/DRM rates
    with chromatin-accessibility and transcription-coupled repair adjustments
  - COSMIC Mutational Signature Integration: SBS1/SBS7a/SBS2/SBS88
    trinucleotide context weights for mutational susceptibility scoring
  - Strand-Asymmetric Repair (TCR): Gene strand inference from ORF analysis
    with transcription-coupled repair enhancement on the template strand
  - Dose-Response Scaling: Exposure-level parameter on all detection
    functions to model UV dose, oxidative stress, and genotoxic exposure

Core Detection:
  - CpG island identification and methylation probability estimation
  - UV photoproduct hotspots (CPD and 6-4PP)
  - Oxidative damage (8-oxoguanine) susceptibility
  - Alu element detection (seed-based consensus matching)
  - 5-hydroxymethylcytosine (5hmC) hotspot detection
  - Methylation-adjusted risk computation (BDNF, gene body, island/shore/shelf)
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Optional

__all__ = [
    "DNADamageReport",
    "CpGIsland",
    "DamageHotspot",
    "RepairPathwayRate",
    "check_dna_degradation",
    "detect_cpg_islands",
    "detect_uv_hotspots",
    "detect_8oxog_hotspots",
    "detect_alu_elements",
    "detect_5hmc_hotspots",
    "estimate_methylation_probability",
    "compute_methylation_adjusted_risk",
    "compute_net_mutation_risk",
    "compute_cosmic_context_weights",
    "infer_gene_strand",
]

logger = logging.getLogger(__name__)


# ==============================================================================
# 1. CpG Island Detection Constants
# ==============================================================================

_CPG_ISLAND_MIN_LENGTH: int = 200       # Minimum length in bp (Takai & Jones 2002)
_CPG_ISLAND_GC_THRESHOLD: float = 0.50  # Minimum GC content
_CPG_ISLAND_OBS_EXP_RATIO: float = 0.60 # Minimum Obs/Exp CpG ratio
_CPG_SHORE_DISTANCE: int = 2000         # Distance from island for "shore" (Irizarry 2009)
_CPG_SHELF_DISTANCE: int = 4000         # Distance from island for "shelf" (Irizarry 2009)

# UV photoproduct constants
_UV_CPD_DINUCS: list[str] = ["TT", "TC", "CT", "CC"]  # Cyclobutane pyrimidine dimers
_UV_64PP_DINUCS: list[str] = ["TC", "CC", "TT"]        # 6-4 photoproducts
_UV_HOTSPOT_WINDOW: int = 100           # Window for local UV damage density

# 8-oxoguanine constants
_8OXOG_CONTEXT_FAVORING: list[str] = ["GCG", "GCA", "GTG", "GGG"]  # Promiscuous oxidation contexts
_8OXOG_HOTSPOT_THRESHOLD: float = 0.03  # Minimum density for hotspot calling

# 5-hydroxymethylcytosine constants
_5HMC_PROMOTER_MOTIFS: list[str] = ["GCGC", "CGCG", "GCGG", "CCGC"]
_5HMC_ENHANCER_MOTIFS: list[str] = ["GATA", "CAAT", "TGCA", "ATGC"]
_5HMC_WINDOW_SIZE: int = 150

# Alu element detection constants (from RepBase AluJb/AluJo/AluSx consensus)
_ALU_LEFT_SEEDS: list[str] = ["GCGCCG", "GGCCGG", "GCTCCG"]
_ALU_RIGHT_SEEDS: list[str] = ["CCGCTC", "CCGCCT", "CCACCT"]
_ALU_SEED_FLOOR: float = 0.65
_ALU_MIN_SEED_MATCH: int = 5
_ALU_LEFT_CONSENSUS: str = "GCGCCGCGCC"
_ALU_ELEMENT_LENGTH: int = 280
_ALU_ARICH_LINKER_LEN: int = 20

# Methylation probability ranges by genomic context (Irizarry et al. 2009)
_ISLAND_METHYLATION_RANGE: tuple[float, float] = (0.05, 0.30)
_SHORE_METHYLATION_RANGE: tuple[float, float] = (0.30, 0.60)
_SHELF_METHYLATION_RANGE: tuple[float, float] = (0.60, 0.80)
_OPEN_SEA_METHYLATION_RANGE: tuple[float, float] = (0.70, 0.90)
_GENE_BODY_METHYLATION_RANGE: tuple[float, float] = (0.40, 0.60)
_GENE_BODY_METHYLATION_DEFAULT: float = 0.50
_BDNF_MOTIF_OFFSET: int = 5
_BDNF_METHYLATION_BOOST: float = 0.10
_ALU_METHYLATION_BOOST: float = 0.15


# ==============================================================================
# 2. Data Classes
# ==============================================================================

@dataclass
class CpGIsland:
    """Represents a detected CpG island."""
    start: int
    end: int
    length: int
    gc_content: float
    obs_exp_ratio: float
    context: str = "island"  # island, shore, shelf, open_sea

    @property
    def methylation_context(self) -> str:
        """Return methylation context classification."""
        return self.context


@dataclass
class DamageHotspot:
    """A detected DNA damage hotspot."""
    position: int
    damage_type: str         # cpg_deamination, uv_cpd, uv_64pp, 8oxog, alu, 5hmc
    severity: float          # 0-1
    context: str = ""        # Surrounding sequence context
    strand: str = "+"        # + or -
    repair_pathway: str = "" # BER, NER, MMR, DRM, etc.

    def adjusted_severity(self, exposure_level: float = 1.0) -> float:
        """Return severity scaled by exposure level (dose-response)."""
        return min(1.0, self.severity * exposure_level)


@dataclass
class DNADamageReport:
    """Comprehensive DNA damage assessment report."""
    sequence_length: int
    cpg_islands: list[CpGIsland] = field(default_factory=list)
    hotspots: list[DamageHotspot] = field(default_factory=list)
    alu_elements: list[dict] = field(default_factory=list)
    overall_risk: float = 0.0
    exposure_level: float = 1.0
    organism: str = "human"

    @property
    def hotspot_count(self) -> int:
        return len(self.hotspots)

    @property
    def high_risk_count(self) -> int:
        return sum(1 for h in self.hotspots if h.adjusted_severity(self.exposure_level) > 0.7)


# ==============================================================================
# 3. CpG Island Detection
# ==============================================================================

def detect_cpg_islands(seq: str, window: int = _CPG_ISLAND_MIN_LENGTH,
                       gc_threshold: float = _CPG_ISLAND_GC_THRESHOLD,
                       obs_exp_threshold: float = _CPG_ISLAND_OBS_EXP_RATIO) -> list[CpGIsland]:
    """Detect CpG islands using the Takai-Jones criteria.

    A CpG island is defined as a region >= 200bp with GC content >= 50%
    and observed/expected CpG ratio >= 0.60.

    Args:
        seq: DNA sequence (uppercase).
        window: Sliding window size in bp.
        gc_threshold: Minimum GC fraction.
        obs_exp_threshold: Minimum Obs/Exp CpG ratio.

    Returns:
        List of CpGIsland objects.
    """
    seq = seq.upper()
    islands: list[CpGIsland] = []
    n = len(seq)

    if n < window:
        return islands

    i = 0
    while i <= n - window:
        region = seq[i:i + window]
        gc = (region.count("G") + region.count("C")) / len(region)

        if gc >= gc_threshold:
            # Calculate Obs/Exp CpG ratio
            cpg_count = sum(1 for j in range(len(region) - 1) if region[j:j + 2] == "CG")
            c_count = region.count("C")
            g_count = region.count("G")
            expected = (c_count * g_count) / (4 * len(region)) if c_count > 0 and g_count > 0 else 0
            obs_exp = cpg_count / expected if expected > 0 else 0.0

            if obs_exp >= obs_exp_threshold:
                # Extend island as far as criteria hold
                end = i + window
                while end < n:
                    ext_region = seq[i:end + 1]
                    ext_gc = (ext_region.count("G") + ext_region.count("C")) / len(ext_region)
                    ext_cpg = sum(1 for j in range(len(ext_region) - 1) if ext_region[j:j + 2] == "CG")
                    ext_c = ext_region.count("C")
                    ext_g = ext_region.count("G")
                    ext_expected = (ext_c * ext_g) / (4 * len(ext_region)) if ext_c > 0 and ext_g > 0 else 0
                    ext_obs_exp = ext_cpg / ext_expected if ext_expected > 0 else 0.0
                    if ext_gc >= gc_threshold and ext_obs_exp >= obs_exp_threshold:
                        end += 1
                    else:
                        break

                islands.append(CpGIsland(
                    start=i, end=end, length=end - i,
                    gc_content=gc, obs_exp_ratio=obs_exp, context="island"
                ))
                i = end
                continue
        i += 1

    # Classify surrounding regions as shore/shelf
    classified = islands.copy()
    for island in islands:
        for offset_start, offset_end, ctx in [
            (max(0, island.start - _CPG_SHORE_DISTANCE), island.start, "shore"),
            (island.end, min(n, island.end + _CPG_SHORE_DISTANCE), "shore"),
            (max(0, island.start - _CPG_SHELF_DISTANCE), max(0, island.start - _CPG_SHORE_DISTANCE), "shelf"),
            (min(n, island.end + _CPG_SHORE_DISTANCE), min(n, island.end + _CPG_SHELF_DISTANCE), "shelf"),
        ]:
            if offset_end > offset_start:
                region = seq[offset_start:offset_end]
                gc = (region.count("G") + region.count("C")) / len(region) if region else 0
                classified.append(CpGIsland(
                    start=offset_start, end=offset_end, length=offset_end - offset_start,
                    gc_content=gc, obs_exp_ratio=0.0, context=ctx
                ))

    return classified


# ==============================================================================
# 4. UV Photoproduct Hotspot Detection
# ==============================================================================

def detect_uv_hotspots(seq: str, exposure_level: float = 1.0) -> list[DamageHotspot]:
    """Detect UV damage hotspots (CPD and 6-4PP) in a DNA sequence.

    Args:
        seq: DNA sequence (uppercase).
        exposure_level: UV dose scaling factor (1.0 = standard exposure).

    Returns:
        List of DamageHotspot objects for UV-susceptible positions.
    """
    seq = seq.upper()
    hotspots: list[DamageHotspot] = []
    n = len(seq)

    for i in range(n - 1):
        dinuc = seq[i:i + 2]

        if dinuc in _UV_CPD_DINUCS:
            # CPD severity depends on dinucleotide identity
            severity_map = {"TT": 0.8, "TC": 0.9, "CT": 0.7, "CC": 0.85}
            base_severity = severity_map.get(dinuc, 0.5)
            adjusted = min(1.0, base_severity * exposure_level)
            hotspots.append(DamageHotspot(
                position=i, damage_type="uv_cpd", severity=adjusted,
                context=seq[max(0, i - 2):min(n, i + 4)], strand="+", repair_pathway="NER"
            ))

        if dinuc in _UV_64PP_DINUCS:
            severity_map = {"TC": 0.7, "CC": 0.6, "TT": 0.5}
            base_severity = severity_map.get(dinuc, 0.4)
            adjusted = min(1.0, base_severity * exposure_level)
            hotspots.append(DamageHotspot(
                position=i, damage_type="uv_64pp", severity=adjusted,
                context=seq[max(0, i - 2):min(n, i + 4)], strand="+", repair_pathway="NER"
            ))

    # Merge nearby hotspots within window
    merged: list[DamageHotspot] = []
    for h in hotspots:
        if merged and abs(h.position - merged[-1].position) < _UV_HOTSPOT_WINDOW // 2:
            # Keep the more severe one
            if h.severity > merged[-1].severity:
                merged[-1] = h
        else:
            merged.append(h)

    return merged


# ==============================================================================
# 5. 8-Oxoguanine Hotspot Detection
# ==============================================================================

def detect_8oxog_hotspots(seq: str, exposure_level: float = 1.0) -> list[DamageHotspot]:
    """Detect 8-oxoguanine (oxidative damage) susceptibility hotspots.

    8-oxoG is the most common oxidative DNA lesion. It preferentially
    forms in contexts with consecutive purines, especially 5'-GpG-3'.

    Args:
        seq: DNA sequence (uppercase).
        exposure_level: Oxidative stress scaling factor (1.0 = baseline).

    Returns:
        List of DamageHotspot objects for 8-oxoG susceptible positions.
    """
    seq = seq.upper()
    hotspots: list[DamageHotspot] = []
    n = len(seq)

    for i in range(n):
        if seq[i] != "G":
            continue

        # Check trinucleotide context
        trinuc_start = max(0, i - 1)
        trinuc_end = min(n, i + 2)
        trinuc = seq[trinuc_start:trinuc_end]

        base_severity = 0.1  # baseline for isolated G

        # Context-dependent severity
        if trinuc in _8OXOG_CONTEXT_FAVORING:
            base_severity = 0.6
        elif "GG" in seq[max(0, i - 1):min(n, i + 2)]:
            base_severity = 0.5  # Runs of G are susceptible
        elif i > 0 and seq[i - 1] in ("A", "G"):
            base_severity = 0.3  # Purine context

        adjusted = min(1.0, base_severity * exposure_level)

        if adjusted >= _8OXOG_HOTSPOT_THRESHOLD:
            hotspots.append(DamageHotspot(
                position=i, damage_type="8oxog", severity=adjusted,
                context=seq[max(0, i - 3):min(n, i + 4)], strand="+", repair_pathway="BER"
            ))

    return hotspots


# ==============================================================================
# 6. Alu Element Detection
# ==============================================================================

def detect_alu_elements(seq: str) -> list[dict]:
    """Detect Alu repetitive elements using seed-based consensus matching.

    Uses left and right arm seed sequences derived from RepBase AluJb/AluJo/AluSx
    consensus sequences. A candidate region is flagged when sufficient seed matches
    occur in the expected orientation and spacing.

    Args:
        seq: DNA sequence (uppercase).

    Returns:
        List of dicts with keys: start, end, confidence, left_matches, right_matches
    """
    seq = seq.upper()
    n = len(seq)
    elements: list[dict] = []

    if n < _ALU_ELEMENT_LENGTH:
        return elements

    for i in range(n - _ALU_ELEMENT_LENGTH):
        candidate = seq[i:i + _ALU_ELEMENT_LENGTH]

        # Left arm: first ~120bp with A-rich linker at ~130-150
        left_arm = candidate[:120]
        left_matches = 0
        for seed in _ALU_LEFT_SEEDS:
            pos = 0
            while True:
                found = left_arm.find(seed, pos)
                if found == -1:
                    break
                left_matches += 1
                pos = found + 1

        # Check left arm consensus match
        left_consensus_matches = sum(
            1 for a, b in zip(left_arm[:len(_ALU_LEFT_CONSENSUS)], _ALU_LEFT_CONSENSUS)
            if a == b
        ) / len(_ALU_LEFT_CONSENSUS)

        # Right arm: after A-rich linker (~position 150-280)
        right_arm_start = 130 + _ALU_ARICH_LINKER_LEN
        if right_arm_start >= len(candidate):
            continue
        right_arm = candidate[right_arm_start:]
        right_matches = 0
        for seed in _ALU_RIGHT_SEEDS:
            pos = 0
            while True:
                found = right_arm.find(seed, pos)
                if found == -1:
                    break
                right_matches += 1
                pos = found + 1

        total_matches = left_matches + right_matches

        # Check if enough seed matches and consensus quality
        if (total_matches >= _ALU_MIN_SEED_MATCH
                and left_consensus_matches >= _ALU_SEED_FLOOR):
            confidence = min(1.0, total_matches / 15.0)
            elements.append({
                "start": i,
                "end": i + _ALU_ELEMENT_LENGTH,
                "confidence": confidence,
                "left_matches": left_matches,
                "right_matches": right_matches,
            })

            # Skip ahead to avoid overlapping detections
            i += _ALU_ELEMENT_LENGTH

    return elements


# ==============================================================================
# 7. 5-Hydroxymethylcytosine (5hmC) Hotspot Detection
# ==============================================================================

def detect_5hmc_hotspots(seq: str, exposure_level: float = 1.0) -> list[DamageHotspot]:
    """Detect 5-hydroxymethylcytosine (5hmC) enrichment hotspots.

    5hmC is enriched at promoter and enhancer regions of active genes.
    This function scans for CpG sites within motif contexts associated
    with TET-mediated hydroxymethylation.

    Args:
        seq: DNA sequence (uppercase).
        exposure_level: Scaling factor for TET activity / exposure.

    Returns:
        List of DamageHotspot for 5hmC-enriched positions.
    """
    seq = seq.upper()
    hotspots: list[DamageHotspot] = []
    n = len(seq)

    for i in range(n - 1):
        if seq[i:i + 2] != "CG":
            continue

        # Check promoter context
        window_start = max(0, i - _5HMC_WINDOW_SIZE // 2)
        window_end = min(n, i + _5HMC_WINDOW_SIZE // 2)
        window = seq[window_start:window_end]

        promoter_score = sum(1 for m in _5HMC_PROMOTER_MOTIFS if m in window) / len(_5HMC_PROMOTER_MOTIFS)
        enhancer_score = sum(1 for m in _5HMC_ENHANCER_MOTIFS if m in window) / len(_5HMC_ENHANCER_MOTIFS)

        base_severity = max(promoter_score * 0.4, enhancer_score * 0.25)
        adjusted = min(1.0, base_severity * exposure_level)

        if adjusted >= 0.05:
            hotspots.append(DamageHotspot(
                position=i, damage_type="5hmc", severity=adjusted,
                context=window[:20], strand="+", repair_pathway="DRM"
            ))

    return hotspots


# ==============================================================================
# 8. Methylation Probability Estimation
# ==============================================================================

def estimate_methylation_probability(
    seq: str,
    position: int,
    cpg_islands: list[CpGIsland] | None = None,
    gene_body: bool = False,
    alu_overlaps: list[dict] | None = None,
) -> float:
    """Estimate DNA methylation probability at a CpG site.

    Uses genomic context (island/shore/shelf/open_sea), gene body status,
    and Alu overlap to estimate methylation probability based on
    published ranges from Irizarry et al. 2009 and related studies.

    Args:
        seq: DNA sequence (uppercase).
        position: 0-indexed position of the CpG site.
        cpg_islands: Pre-detected CpG islands with context classifications.
        gene_body: Whether this position falls within a gene body.
        alu_overlaps: Pre-detected Alu elements overlapping this region.

    Returns:
        Estimated methylation probability (0-1).
    """
    # Determine genomic context
    context = "open_sea"
    if cpg_islands:
        for island in cpg_islands:
            if island.start <= position < island.end:
                context = island.context
                break

    # Base methylation probability from context
    context_ranges = {
        "island": _ISLAND_METHYLATION_RANGE,
        "shore": _SHORE_METHYLATION_RANGE,
        "shelf": _SHELF_METHYLATION_RANGE,
        "open_sea": _OPEN_SEA_METHYLATION_RANGE,
    }
    low, high = context_ranges.get(context, _OPEN_SEA_METHYLATION_RANGE)

    # Sequence-feature-based methylation probability estimation
    # Uses CpG density and local GC content instead of random sampling.
    # Based on Irizarry et al. 2009 ranges for different genomic contexts.
    window = 200  # bp window for local context
    start = max(0, position - window // 2)
    end = min(len(seq), position + window // 2)
    local_seq = seq[start:end].upper()

    # CpG density (observed/expected)
    c_count = local_seq.count('C')
    g_count = local_seq.count('G')
    cpg_count = local_seq.count('CG')
    expected_cpg = (c_count * g_count) / max(len(local_seq), 1)
    cpg_oe = cpg_count / max(expected_cpg, 0.001)

    # Local GC content
    gc_content = (c_count + g_count) / max(len(local_seq), 1)

    # Interpolate within range based on CpG O/E and GC content
    # Higher CpG O/E = lower methylation (in islands)
    # Higher GC = lower methylation
    factor = 1.0 - (cpg_oe * 0.1 + gc_content * 0.3)
    factor = max(0.0, min(1.0, factor))
    methylation_prob = low + (high - low) * factor

    # Gene body methylation adjustment
    if gene_body:
        body_low, body_high = _GENE_BODY_METHYLATION_RANGE
        body_factor = 1.0 - (cpg_oe * 0.1 + gc_content * 0.3)
        body_factor = max(0.0, min(1.0, body_factor))
        body_prob = body_low + (body_high - body_low) * body_factor
        methylation_prob = (methylation_prob + body_prob) / 2.0

    # Alu element methylation boost
    if alu_overlaps:
        for alu in alu_overlaps:
            if alu["start"] <= position < alu["end"]:
                methylation_prob = min(1.0, methylation_prob + _ALU_METHYLATION_BOOST)
                break

    return min(1.0, methylation_prob)


# ==============================================================================
# 9. Methylation-Adjusted Risk Computation
# ==============================================================================

def compute_methylation_adjusted_risk(
    seq: str,
    position: int,
    damage_type: str,
    cpg_islands: list[CpGIsland] | None = None,
    gene_body: bool = False,
    alu_overlaps: list[dict] | None = None,
    exposure_level: float = 1.0,
) -> dict:
    """Compute methylation-adjusted DNA damage risk at a given position.

    Combines methylation probability with damage severity to produce a
    net risk score that accounts for the interaction between methylation
    status and damage susceptibility. Methylated CpGs are both more prone
    to deamination and less efficiently repaired.

    Special handling for BDNF promoter motif (CACCGG), which shows
    elevated methylation-dependent damage in neuronal contexts.

    Args:
        seq: DNA sequence (uppercase).
        position: 0-indexed position to assess.
        damage_type: Type of damage (cpg_deamination, 8oxog, uv_cpd, etc.).
        cpg_islands: Pre-detected CpG islands.
        gene_body: Whether position is in a gene body.
        alu_overlaps: Pre-detected Alu elements.
        exposure_level: Dose-response scaling factor.

    Returns:
        Dict with keys: net_risk, methylation_prob, base_severity,
        adjusted_severity, damage_type, position, bdnf_motif
    """
    seq = seq.upper()
    n = len(seq)

    # Base severity by damage type
    severity_map = {
        "cpg_deamination": 0.6,
        "8oxog": 0.5,
        "uv_cpd": 0.7,
        "uv_64pp": 0.6,
        "5hmc": 0.3,
    }
    base_severity = severity_map.get(damage_type, 0.3)

    # Dose-response scaling
    adjusted_severity = min(1.0, base_severity * exposure_level)

    # Estimate methylation probability
    methylation_prob = estimate_methylation_probability(
        seq, position, cpg_islands, gene_body, alu_overlaps
    )

    # Check for BDNF promoter motif (CACCGG) nearby
    bdnf_motif = False
    check_start = max(0, position - _BDNF_MOTIF_OFFSET - 6)
    check_end = min(n, position + _BDNF_MOTIF_OFFSET + 6)
    if "CACCGG" in seq[check_start:check_end]:
        bdnf_motif = True
        methylation_prob = min(1.0, methylation_prob + _BDNF_METHYLATION_BOOST)

    # Methylation amplifies CpG deamination risk
    if damage_type == "cpg_deamination" and methylation_prob > 0.5:
        # Methylated CpG deamination is 2-5x more likely
        deamination_multiplier = 1.0 + methylation_prob
        adjusted_severity = min(1.0, adjusted_severity * deamination_multiplier)

    # Net risk = adjusted severity weighted by methylation interaction
    if damage_type == "cpg_deamination":
        net_risk = adjusted_severity * (0.5 + 0.5 * methylation_prob)
    elif damage_type == "8oxog":
        net_risk = adjusted_severity * (1.0 - 0.3 * methylation_prob)  # methylation slightly protects
    else:
        net_risk = adjusted_severity

    return {
        "net_risk": net_risk,
        "methylation_prob": methylation_prob,
        "base_severity": base_severity,
        "adjusted_severity": adjusted_severity,
        "damage_type": damage_type,
        "position": position,
        "bdnf_motif": bdnf_motif,
    }


# ==============================================================================
# 10. Main DNA Degradation Check
# ==============================================================================

def check_dna_degradation(
    seq: str,
    organism: str = "human",
    exposure_level: float = 1.0,
    check_uv: bool = True,
    check_8oxog: bool = True,
    check_alu: bool = True,
    check_5hmc: bool = True,
    check_methylation: bool = True,
) -> DNADamageReport:
    """Comprehensive DNA damage and degradation assessment.

    Runs all detection modules and produces a consolidated damage report
    with dose-response scaling via the exposure_level parameter.

    Args:
        seq: DNA sequence (uppercase or mixed case).
        organism: Target organism (human, mouse, etc.).
        exposure_level: Dose-response scaling factor (1.0 = baseline).
        check_uv: Whether to check for UV damage hotspots.
        check_8oxog: Whether to check for oxidative damage hotspots.
        check_alu: Whether to check for Alu elements.
        check_5hmc: Whether to check for 5hmC hotspots.
        check_methylation: Whether to compute methylation-adjusted risks.

    Returns:
        DNADamageReport with all findings.
    """
    seq = seq.upper()
    report = DNADamageReport(
        sequence_length=len(seq),
        organism=organism,
        exposure_level=exposure_level,
    )

    # Detect CpG islands first (needed for methylation context)
    report.cpg_islands = detect_cpg_islands(seq)

    # Detect Alu elements (needed for methylation boost)
    if check_alu:
        report.alu_elements = detect_alu_elements(seq)

    # Detect UV hotspots
    if check_uv:
        uv_hotspots = detect_uv_hotspots(seq, exposure_level=exposure_level)
        report.hotspots.extend(uv_hotspots)

    # Detect 8-oxoG hotspots
    if check_8oxog:
        oxog_hotspots = detect_8oxog_hotspots(seq, exposure_level=exposure_level)
        report.hotspots.extend(oxog_hotspots)

    # Detect 5hmC hotspots
    if check_5hmc:
        hmc_hotspots = detect_5hmc_hotspots(seq, exposure_level=exposure_level)
        report.hotspots.extend(hmc_hotspots)

    # Add CpG deamination hotspots for each CpG dinucleotide
    # Context-dependent rates from Alexandrov et al. 2013 (COSMIC SBS1)
    _CPG_DEAMINATION_RATES = {
        "CCG": 0.45,   # enhanced deamination in CCG context
        "TCG": 0.40,   # TCG context
        "ACG": 0.30,   # ACG context
        "GCG": 0.30,   # GCG context
    }
    _CPG_BASELINE_RATE = 0.35  # CpG→TpG baseline

    for i in range(len(seq) - 1):
        if seq[i:i + 2] == "CG":
            # Determine trinucleotide context for severity
            if i > 0:
                trinuc = seq[i - 1:i + 2]
                severity = _CPG_DEAMINATION_RATES.get(trinuc, _CPG_BASELINE_RATE)
            else:
                severity = _CPG_BASELINE_RATE
            adjusted = min(1.0, severity * exposure_level)
            report.hotspots.append(DamageHotspot(
                position=i, damage_type="cpg_deamination", severity=adjusted,
                context=seq[max(0, i - 2):min(len(seq), i + 4)],
                strand="+", repair_pathway="DRM"
            ))

    # Compute methylation-adjusted risk for CpG positions
    if check_methylation:
        for hotspot in report.hotspots:
            if hotspot.damage_type == "cpg_deamination":
                risk = compute_methylation_adjusted_risk(
                    seq, hotspot.position, hotspot.damage_type,
                    cpg_islands=report.cpg_islands,
                    alu_overlaps=report.alu_elements,
                    exposure_level=exposure_level,
                )
                hotspot.severity = risk["net_risk"]

    # Compute net mutation risk accounting for repair capacity
    # Uses literature-derived repair rates from _REPAIR_RATES
    for hotspot in report.hotspots:
        # Determine if CpG site is methylated from methylation probability
        methylated = False
        if hotspot.damage_type == "cpg_deamination" and check_methylation:
            meth_prob = estimate_methylation_probability(
                seq, hotspot.position,
                cpg_islands=report.cpg_islands,
                alu_overlaps=report.alu_elements,
            )
            methylated = meth_prob > 0.5

        # Map damage_type to compute_net_mutation_risk keys
        damage_type_key = hotspot.damage_type
        if damage_type_key == "uv_cpd":
            damage_type_key = "uv_cpd"
        elif damage_type_key == "uv_64pp":
            damage_type_key = "uv_6-4pp"

        net_result = compute_net_mutation_risk(
            damage_severity=hotspot.severity,
            damage_type=damage_type_key,
            methylated=methylated,
        )
        # Reduce hotspot severity by repair capacity
        hotspot.severity = net_result["net_risk"]
        hotspot.repair_pathway = net_result["repair_pathway"]

    # Compute overall risk
    if report.hotspots:
        severity_sum = sum(h.adjusted_severity(exposure_level) for h in report.hotspots)
        report.overall_risk = min(1.0, severity_sum / len(seq) * 100)

    return report


# ==============================================================================
# 11. SOTA Upgrade 1: DNA Repair Pathway Modeling
# ==============================================================================

@dataclass
class RepairPathwayRate:
    """Repair efficiency by pathway, derived from literature."""
    pathway: str
    lesion_type: str
    repair_efficiency_24h: float  # fraction repaired in 24 hours
    half_life_hours: float
    reference: str


# Literature-derived repair rates
_REPAIR_RATES: list[RepairPathwayRate] = [
    RepairPathwayRate("BER", "8oxog", 0.90, 4.0, "Dianov & Hübscher 2013 DNA Repair"),
    RepairPathwayRate("BER", "abasic_site", 0.95, 2.0, "Boiteux & Guillet 2004 DNA Repair"),
    RepairPathwayRate("NER", "uv_cpd", 0.50, 24.0, "Mouret et al. 2008 DNA Repair"),
    RepairPathwayRate("NER", "uv_6-4pp", 0.80, 4.0, "Mitchell et al. 1985 PNAS"),
    RepairPathwayRate("NER", "tc_ner_cpd", 0.90, 8.0, "Hu et al. 2017 PNAS"),
    RepairPathwayRate("MMR", "mismatch", 0.95, 1.0, "Jiricny 2006 Nat Rev Mol Cell Biol"),
    RepairPathwayRate("MMR", "slippage", 0.90, 2.0, "Kunkel & Erie 2005 Annu Rev Biochem"),
    RepairPathwayRate("DRM", "cpg_deamination_unmethylated", 0.95, 1.0, "Lindahl 1993 Nature"),
    RepairPathwayRate("DRM", "cpg_deamination_methylated", 0.30, 24.0, "Pfeifer 2006 Mutat Res"),
]


def compute_net_mutation_risk(damage_severity: float, damage_type: str,
                              methylated: bool = False,
                              transcription_active: bool = False,
                              chromatin_open: bool = True) -> dict:
    """Compute net mutation risk = formation_rate × (1 - repair_efficiency).

    Args:
        damage_severity: 0-1 severity score from detection
        damage_type: Type key (cpg_deamination, 8oxog, uv_photoproduct, etc.)
        methylated: Whether CpG is methylated (affects deamination repair)
        transcription_active: Whether gene is actively transcribed (TCR boost)
        chromatin_open: Whether chromatin is open (affects repair access)

    Returns:
        Dict with net_risk, repair_rate, formation_rate, pathway
    """
    # Find matching repair rate
    repair_eff = 0.5  # default
    pathway = "unknown"
    half_life = 24.0

    for r in _REPAIR_RATES:
        if r.lesion_type == damage_type or damage_type.startswith(r.lesion_type.split("_")[0]):
            repair_eff = r.repair_efficiency_24h
            pathway = r.pathway
            half_life = r.half_life_hours
            break

    # Adjust for methylation (methylated CpG deamination is poorly repaired)
    if damage_type == "cpg_deamination" and methylated:
        for r in _REPAIR_RATES:
            if r.lesion_type == "cpg_deamination_methylated":
                repair_eff = r.repair_efficiency_24h
                break

    # TCR boost: transcribed strand gets 2-10x NER enhancement
    if transcription_active and pathway == "NER":
        repair_eff = min(0.99, repair_eff * 3.0)

    # Chromatin accessibility: open chromatin = better repair access
    if not chromatin_open:
        repair_eff = max(0.05, repair_eff * 0.5)

    formation_rate = damage_severity
    net_risk = formation_rate * (1.0 - repair_eff)

    return {
        "net_risk": net_risk,
        "formation_rate": formation_rate,
        "repair_efficiency": repair_eff,
        "repair_pathway": pathway,
        "repair_half_life_hours": half_life,
    }


# ==============================================================================
# 12. SOTA Upgrade 2: COSMIC Mutational Signature Integration
# ==============================================================================

# COSMIC SBS signature trinucleotide contexts (top signatures relevant to therapeutics)
_SBS1_WEIGHTS: dict[str, float] = {  # CpG deamination
    "ACG": 0.08, "CCG": 0.12, "GCG": 0.10, "TCG": 0.15,
    "AAG": 0.02, "AGG": 0.03, "CGG": 0.04, "TGG": 0.02,
}
_SBS7A_WEIGHTS: dict[str, float] = {  # UV damage
    "TCA": 0.03, "TCC": 0.04, "TCG": 0.02, "TCT": 0.08,
    "TTA": 0.04, "TTC": 0.05, "TTG": 0.03, "TTT": 0.06,
}
_SBS2_WEIGHTS: dict[str, float] = {  # APOBEC
    "TCA": 0.10, "TCT": 0.08, "CCA": 0.06, "CCC": 0.04,
}
_SBS88_WEIGHTS: dict[str, float] = {  # E. coli colibactin (for gut-targeted)
    "ACA": 0.05, "ACT": 0.04, "ATA": 0.03, "ATT": 0.04,
}


def compute_cosmic_context_weights(seq: str, position: int,
                                    signatures: list[str] | None = None) -> float:
    """Compute COSMIC SBS signature weight at a given position.

    Uses trinucleotide context to assess mutational susceptibility
    based on known COSMIC single-base substitution signatures.

    Args:
        seq: Full DNA sequence
        position: 0-indexed position in sequence
        signatures: Which SBS signatures to consider (default: SBS1, SBS7a, SBS2)

    Returns:
        Combined signature weight (0-1 scale)
    """
    if signatures is None:
        signatures = ["SBS1", "SBS7a", "SBS2"]

    if position < 1 or position >= len(seq) - 1:
        return 0.0

    trinuc = seq[position - 1:position + 2].upper()

    weight = 0.0
    sig_map = {"SBS1": _SBS1_WEIGHTS, "SBS7a": _SBS7A_WEIGHTS,
               "SBS2": _SBS2_WEIGHTS, "SBS88": _SBS88_WEIGHTS}

    for sig in signatures:
        if sig in sig_map and trinuc in sig_map[sig]:
            weight += sig_map[sig][trinuc]

    return min(1.0, weight)


# ==============================================================================
# 13. SOTA Upgrade 3: Strand-Asymmetric Repair (TCR) Helper
# ==============================================================================

def infer_gene_strand(seq: str) -> str:
    """Infer gene strand from ORF analysis.

    Simple heuristic: if the first ATG is in the forward direction
    and produces a longer ORF than the reverse complement, assume '+'.
    """
    fwd_atg = seq.find("ATG")
    rc = seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
    rev_atg = rc.find("ATG")

    if fwd_atg == -1 and rev_atg == -1:
        return "+"

    # Count ORF lengths
    def orf_length(s: str, start: int) -> int:
        for i in range(start, len(s) - 2, 3):
            codon = s[i:i + 3]
            if codon in ("TAA", "TAG", "TGA"):
                return i - start + 3
        return len(s) - start

    fwd_len = orf_length(seq, fwd_atg) if fwd_atg >= 0 else 0
    rev_len = orf_length(rc, rev_atg) if rev_atg >= 0 else 0

    return "+" if fwd_len >= rev_len else "-"
