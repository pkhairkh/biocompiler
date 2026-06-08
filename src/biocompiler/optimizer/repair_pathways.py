"""
BioCompiler DNA Repair Pathway Kinetics Module v1.0.0
=====================================================
Comprehensive DNA repair pathway modeling for therapeutic gene design,
covering five major repair pathways with literature-derived kinetic
parameters.

Repair Pathways Modeled:
  - BER (Base Excision Repair): 8-oxoG, uracil, abasic sites, alkylated bases
  - NER (Nucleotide Excision Repair): CPD, 6-4PP, bulky adducts, TCR-enhanced
  - MMR (Mismatch Repair): base-base mismatches, indel loops, slippage
  - HR (Homologous Recombination): DSB in S/G2, replication fork collapse
  - NHEJ (Non-Homologous End Joining): DSB in G1, alt-NHEJ/MMEJ

Key Features:
  - Literature-derived repair efficiencies and half-lives for 16 lesion types
  - Chromatin accessibility adjustment (open, heterochromatin, boundary)
  - Transcription-coupled repair (TCR) enhancement for bulky lesions
  - Cell cycle-dependent pathway selection (HR vs NHEJ for DSBs)
  - Methylation-aware repair prediction for CpG deamination

References:
  - Lindahl T (1993) Nature 362:709-715 — Instability and decay of DNA
  - Wood RD et al. (1996) Science 271:1663 — DNA repair in eukaryotes
  - Jiricny J (2006) Nat Rev Mol Cell Biol 7:335-346 — MMR
  - Dianov GL & Hubscher U (2013) DNA Repair 12:593-601 — BER
  - Sancar A (1996) Annu Rev Biochem 65:43-81 — NER
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Optional

__all__ = [
    "RepairPathway",
    "RepairPrediction",
    "predict_repair",
    "compute_net_mutation_risk",
    "get_pathway_for_lesion",
    "adjust_for_chromatin",
    "adjust_for_tcr",
    "REPAIR_PATHWAY_DATABASE",
]

logger = logging.getLogger(__name__)


# ==============================================================================
# 1. Data Classes
# ==============================================================================

@dataclass
class RepairPathway:
    """Describes a DNA repair pathway with its kinetic parameters.

    Attributes:
        pathway_name: Name of the repair pathway (BER, NER, MMR, HR, NHEJ).
        lesion_types: List of lesion type strings this pathway handles.
        repair_efficiency_24h: Fraction of lesions repaired within 24 hours (0-1).
        half_life_hours: Half-life of the lesion under this pathway (hours).
        reference: Literature reference for the kinetic parameters.
    """
    pathway_name: str
    lesion_types: list[str]
    repair_efficiency_24h: float
    half_life_hours: float
    reference: str

    def repair_fraction_at_time(self, t_hours: float) -> float:
        """Compute the fraction of lesions repaired at time *t_hours*.

        Uses first-order kinetics: f(t) = 1 - exp(-ln(2) * t / t_half).

        Args:
            t_hours: Time elapsed since lesion formation (hours).

        Returns:
            Fraction of lesions repaired at the given time (0-1).
        """
        k = math.log(2) / self.half_life_hours if self.half_life_hours > 0 else float("inf")
        fraction = 1.0 - math.exp(-k * t_hours)
        # Cap at the 24h efficiency if t_hours exceeds 24h and efficiency < 1
        if t_hours > 24.0 and self.repair_efficiency_24h < 1.0:
            fraction = min(fraction, self.repair_efficiency_24h)
        return min(fraction, 1.0)


@dataclass
class RepairPrediction:
    """Prediction result for a single lesion's repair outcome.

    Attributes:
        lesion_type: Type of DNA lesion (e.g. '8-oxoG', 'CPD').
        position: Genomic position of the lesion (0-indexed).
        formation_severity: Initial damage severity (0-1).
        repair_efficiency: Predicted repair efficiency (0-1) after adjustments.
        net_risk: Net mutation risk = formation_severity * (1 - repair_efficiency).
        pathway: Name of the repair pathway handling this lesion.
        adjusted_severity: Formation severity after context-specific adjustments.
    """
    lesion_type: str
    position: int
    formation_severity: float
    repair_efficiency: float
    net_risk: float
    pathway: str
    adjusted_severity: float


# ==============================================================================
# 2. Repair Pathway Database
# ==============================================================================

# --- BER (Base Excision Repair) ---
# References: Dianov & Hubscher 2013 DNA Repair; Lindahl 1993 Nature
_BER_OGOG = RepairPathway(
    pathway_name="BER",
    lesion_types=["8-oxoG", "8oxog", "8-oxoguanine"],
    repair_efficiency_24h=0.90,
    half_life_hours=4.0,
    reference="OGG1-initiated; Dianov & Hubscher 2013 DNA Repair; Lindahl 1993 Nature",
)

_BER_URACIL = RepairPathway(
    pathway_name="BER",
    lesion_types=["uracil", "U:G mismatch"],
    repair_efficiency_24h=0.95,
    half_life_hours=1.0,
    reference="UNG-initiated; Lindahl 1993 Nature; Dianov & Hubscher 2013 DNA Repair",
)

_BER_ABASIC = RepairPathway(
    pathway_name="BER",
    lesion_types=["abasic_site", "AP site", "abasic"],
    repair_efficiency_24h=0.95,
    half_life_hours=0.5,
    reference="APE1-initiated; Wood et al. 1996 Science; Dianov & Hubscher 2013 DNA Repair",
)

_BER_ALKYL = RepairPathway(
    pathway_name="BER",
    lesion_types=["alkylated_base", "3-methylA", "7-methylG", "alkyl"],
    repair_efficiency_24h=0.80,
    half_life_hours=8.0,
    reference="AAG/MPG-initiated; Dianov & Hubscher 2013 DNA Repair; Wood et al. 1996 Science",
)

_BER_CPG_DEAMINATION = RepairPathway(
    pathway_name="BER",
    lesion_types=["cpg_deamination", "cpg_deamination_unmethylated", "deamination"],
    repair_efficiency_24h=0.95,
    half_life_hours=1.0,
    reference="TDG/MBD4-initiated; Lindahl 1993 Nature; Dianov & Hubscher 2013 DNA Repair",
)

# --- NER (Nucleotide Excision Repair) ---
# References: Sancar 1996 Annu Rev Biochem; Wood et al. 1996 Science
_NER_CPD = RepairPathway(
    pathway_name="NER",
    lesion_types=["CPD", "cyclobutane_pyrimidine_dimer", "uv_cpd"],
    repair_efficiency_24h=0.50,
    half_life_hours=24.0,
    reference="XPC/TFIIH-initiated GGR; Sancar 1996 Annu Rev Biochem; Wood et al. 1996 Science",
)

_NER_64PP = RepairPathway(
    pathway_name="NER",
    lesion_types=["6-4PP", "6-4_photoproduct", "uv_64pp", "uv_6-4pp"],
    repair_efficiency_24h=0.80,
    half_life_hours=4.0,
    reference="XPC/TFIIH-initiated GGR; Sancar 1996 Annu Rev Biochem; Wood et al. 1996 Science",
)

_NER_BULKY = RepairPathway(
    pathway_name="NER",
    lesion_types=["bulky_adduct", "benzo[a]pyrene", "cisplatin", "aromatic_adduct"],
    repair_efficiency_24h=0.70,
    half_life_hours=12.0,
    reference="XPC/TFIIH-initiated GGR; Sancar 1996 Annu Rev Biochem; Wood et al. 1996 Science",
)

_NER_TCR = RepairPathway(
    pathway_name="NER",
    lesion_types=["TCR-enhanced_CPD", "TCR-enhanced_6-4PP", "tc_ner_cpd"],
    repair_efficiency_24h=0.90,
    half_life_hours=8.0,
    reference="CSB/UVSSA-initiated TCR on transcribed strand; Sancar 1996 Annu Rev Biochem",
)

# --- MMR (Mismatch Repair) ---
# References: Jiricny 2006 Nat Rev Mol Cell Biol
_MMR_MISMATCH = RepairPathway(
    pathway_name="MMR",
    lesion_types=["base-base_mismatch", "mismatch", "G:T mismatch", "A:C mismatch"],
    repair_efficiency_24h=0.95,
    half_life_hours=1.0,
    reference="MSH2/MSH6 (MutSalpha); Jiricny 2006 Nat Rev Mol Cell Biol",
)

_MMR_INDEL = RepairPathway(
    pathway_name="MMR",
    lesion_types=["insertion_deletion_loop", "IDL", "indel_loop", "+1 frameshift", "-1 frameshift"],
    repair_efficiency_24h=0.95,
    half_life_hours=1.0,
    reference="MSH2/MSH6 for 1-2nt loops; Jiricny 2006 Nat Rev Mol Cell Biol",
)

_MMR_SLIPPAGE = RepairPathway(
    pathway_name="MMR",
    lesion_types=["microsatellite_slippage", "slippage", "repeat_instability"],
    repair_efficiency_24h=0.90,
    half_life_hours=2.0,
    reference="MSH2/MSH6 + MLH1/PMS2; Jiricny 2006 Nat Rev Mol Cell Biol",
)

# --- HR (Homologous Recombination) ---
# References: Wood et al. 1996 Science; Jiricny 2006 (DSB context)
_HR_DSG = RepairPathway(
    pathway_name="HR",
    lesion_types=["DSB_SG2", "double_strand_break_SG2", "DSB_replicating"],
    repair_efficiency_24h=0.80,
    half_life_hours=4.0,
    reference="RAD51-mediated; S/G2 phase; Wood et al. 1996 Science",
)

_HR_FORK = RepairPathway(
    pathway_name="HR",
    lesion_types=["replication_fork_collapse", "fork_collapse", "stalled_fork"],
    repair_efficiency_24h=0.70,
    half_life_hours=6.0,
    reference="RAD51-mediated fork restart; Wood et al. 1996 Science",
)

# --- NHEJ (Non-Homologous End Joining) ---
# References: Wood et al. 1996 Science; Lieber 2010 Annu Rev Biochem
_NHEJ_DSG1 = RepairPathway(
    pathway_name="NHEJ",
    lesion_types=["DSB_G1", "double_strand_break_G1", "DSB_nonreplicating"],
    repair_efficiency_24h=0.60,
    half_life_hours=2.0,
    reference="KU70/80-initiated classical NHEJ; G1 phase; Wood et al. 1996 Science",
)

_NHEJ_ALT = RepairPathway(
    pathway_name="NHEJ",
    lesion_types=["alt-NHEJ", "MMEJ", "microhomology_mediated"],
    repair_efficiency_24h=0.50,
    half_life_hours=4.0,
    reference="POLQ-mediated alt-NHEJ/MMEJ; error-prone; Wood et al. 1996 Science",
)

# Master database: all pathway entries
REPAIR_PATHWAY_DATABASE: list[RepairPathway] = [
    # BER
    _BER_OGOG,
    _BER_URACIL,
    _BER_ABASIC,
    _BER_ALKYL,
    _BER_CPG_DEAMINATION,
    # NER
    _NER_CPD,
    _NER_64PP,
    _NER_BULKY,
    _NER_TCR,
    # MMR
    _MMR_MISMATCH,
    _MMR_INDEL,
    _MMR_SLIPPAGE,
    # HR
    _HR_DSG,
    _HR_FORK,
    # NHEJ
    _NHEJ_DSG1,
    _NHEJ_ALT,
]

# Reverse lookup: lesion_type -> RepairPathway
_LESION_LOOKUP: dict[str, RepairPathway] = {}
for _pw in REPAIR_PATHWAY_DATABASE:
    for _lt in _pw.lesion_types:
        _LESION_LOOKUP[_lt.lower()] = _pw


# ==============================================================================
# 3. Chromatin and TCR Adjustment Constants
# ==============================================================================

# Chromatin accessibility multipliers for repair efficiency
# Based on: Mao et al. 2004 J Biol Chem; Hu et al. 2017 PNAS
_CHROMATIN_MULTIPLIERS: dict[str, float] = {
    "open": 1.0,          # Euchromatin — full repair access
    "euchromatin": 1.0,   # Synonym for open
    "boundary": 0.75,     # Boundary/transition zone
    "heterochromatin": 0.5,  # Closed chromatin — impaired repair
    "closed": 0.5,        # Synonym for heterochromatin
    "constitutive_heterochromatin": 0.3,  # Pericentromeric — very impaired
}

# TCR enhancement factor for bulky lesions on transcribed strand
# Based on: Hu et al. 2017 PNAS; Sancar 1996 Annu Rev Biochem
# TCR can increase NER repair 2-10x for CPD on transcribed strand
_TCR_ENHANCEMENT_FACTOR: float = 1.8  # Multiplier applied to repair_efficiency
_TCR_BULKY_LESION_TYPES: set[str] = {
    "cpd", "cyclobutane_pyrimidine_dimer", "uv_cpd",
    "6-4pp", "6-4_photoproduct", "uv_64pp", "uv_6-4pp",
    "bulky_adduct", "benzo[a]pyrene", "cisplatin", "aromatic_adduct",
    "tcr-enhanced_cpd", "tcr-enhanced_6-4pp", "tc_ner_cpd",
}

# Cell cycle phase preferences for DSB repair
# HR is restricted to S/G2 (sister chromatid available)
# NHEJ operates throughout cell cycle but predominates in G1
_CELL_CYCLE_DSB_PATHWAY: dict[str, str] = {
    "G1": "NHEJ",
    "G1/S": "NHEJ",
    "S": "HR",
    "G2": "HR",
    "G2/M": "HR",
    "M": "NHEJ",
}


# ==============================================================================
# 4. Core Functions
# ==============================================================================

def get_pathway_for_lesion(lesion_type: str) -> str:
    """Look up which repair pathway handles a given lesion type.

    Args:
        lesion_type: String identifier for the lesion (e.g. '8-oxoG', 'CPD').

    Returns:
        Name of the repair pathway (BER, NER, MMR, HR, NHEJ),
        or 'unknown' if the lesion type is not in the database.
    """
    pw = _LESION_LOOKUP.get(lesion_type.lower())
    if pw is not None:
        return pw.pathway_name

    # Fuzzy matching: try partial match on common substrings
    lt_lower = lesion_type.lower()
    if "8-oxo" in lt_lower or "8oxo" in lt_lower or "oxog" in lt_lower:
        return "BER"
    if "cpd" in lt_lower or "cyclobutane" in lt_lower:
        return "NER"
    if "6-4" in lt_lower or "64pp" in lt_lower or "photoproduct" in lt_lower:
        return "NER"
    if "bulky" in lt_lower or "adduct" in lt_lower or "cisplatin" in lt_lower:
        return "NER"
    if "uracil" in lt_lower:
        return "BER"
    if "abasic" in lt_lower or "ap site" in lt_lower:
        return "BER"
    if "alkyl" in lt_lower or "methyl" in lt_lower:
        return "BER"
    if "mismatch" in lt_lower:
        return "MMR"
    if "indel" in lt_lower or "loop" in lt_lower or "frameshift" in lt_lower:
        return "MMR"
    if "slippage" in lt_lower or "repeat" in lt_lower or "microsatellite" in lt_lower:
        return "MMR"
    if "dsb" in lt_lower or "double_strand" in lt_lower:
        return "HR"  # Default to HR; cell cycle will refine
    if "fork" in lt_lower:
        return "HR"
    if "nhej" in lt_lower or "mmej" in lt_lower or "microhomology" in lt_lower:
        return "NHEJ"
    if "deamination" in lt_lower or "cpg" in lt_lower:
        return "BER"

    return "unknown"


def adjust_for_chromatin(repair_efficiency: float, chromatin_state: str) -> float:
    """Adjust repair rate for chromatin context.

    Closed/heterochromatin regions impede repair factor access,
    reducing observed repair efficiency. Open/euchromatin regions
    allow full repair access.

    Args:
        repair_efficiency: Baseline repair efficiency (0-1).
        chromatin_state: One of 'open', 'euchromatin', 'boundary',
            'heterochromatin', 'closed', 'constitutive_heterochromatin'.

    Returns:
        Adjusted repair efficiency (0-1).
    """
    multiplier = _CHROMATIN_MULTIPLIERS.get(chromatin_state.lower(), 0.75)
    adjusted = repair_efficiency * multiplier
    return max(0.0, min(1.0, adjusted))


def adjust_for_tcr(repair_efficiency: float, is_transcribed_strand: bool) -> float:
    """Apply transcription-coupled repair (TCR) enhancement for bulky lesions.

    TCR is a dedicated sub-pathway of NER that rapidly repairs
    lesions on the transcribed strand of active genes. It is
    particularly important for CPD, which is poorly repaired by
    global genome NER (GGR) alone.

    This function is intended to be applied to NER-pathway lesions.
    For non-NER lesions, TCR does not apply and the efficiency
    is returned unchanged.

    Args:
        repair_efficiency: Baseline repair efficiency (0-1).
        is_transcribed_strand: Whether the lesion is on the transcribed strand.

    Returns:
        Adjusted repair efficiency (0-1). Enhanced by TCR factor if
        on the transcribed strand, unchanged otherwise.
    """
    if not is_transcribed_strand:
        return repair_efficiency

    enhanced = repair_efficiency * _TCR_ENHANCEMENT_FACTOR
    return min(0.99, enhanced)  # Cap below 1.0 — no repair is perfectly efficient


def predict_repair(
    lesion_type: str,
    methylated: bool = False,
    transcribed_strand: bool = False,
    chromatin_open: bool = True,
    cell_cycle: str = "S",
) -> RepairPrediction:
    """Predict repair outcome for a specific DNA lesion.

    Combines pathway lookup with context-specific adjustments for
    chromatin accessibility, transcription-coupled repair, methylation
    status, and cell cycle phase to produce a comprehensive repair
    prediction.

    Args:
        lesion_type: Type of DNA lesion (e.g. '8-oxoG', 'CPD', 'DSB').
        methylated: Whether the CpG site is methylated (affects deamination repair).
        transcribed_strand: Whether the lesion is on the transcribed strand
            (enables TCR for NER-pathway lesions).
        chromatin_open: Whether the chromatin is in an open/euchromatic state.
            If False, heterochromatin is assumed.
        cell_cycle: Current cell cycle phase ('G1', 'S', 'G2', 'M', 'G1/S',
            'G2/M'). Affects DSB repair pathway selection (HR vs NHEJ).

    Returns:
        RepairPrediction with repair_efficiency, net_risk, pathway assignment,
        and adjusted severity.
    """
    lt_lower = lesion_type.lower()
    pathway_name = get_pathway_for_lesion(lesion_type)

    # --- Resolve DSB pathway based on cell cycle ---
    is_dsb = "dsb" in lt_lower or "double_strand" in lt_lower
    is_nhej_lesion = "nhej" in lt_lower or "mmej" in lt_lower or "microhomology" in lt_lower

    if is_dsb and not is_nhej_lesion:
        preferred_pathway = _CELL_CYCLE_DSB_PATHWAY.get(cell_cycle, "HR")
        # Override pathway if cell cycle dictates different pathway
        pathway_name = preferred_pathway

    # --- Look up base repair efficiency ---
    pw = _LESION_LOOKUP.get(lt_lower)
    if pw is not None:
        repair_eff = pw.repair_efficiency_24h
        half_life = pw.half_life_hours
    else:
        # Default conservative estimate for unknown lesions
        repair_eff = 0.50
        half_life = 12.0
        logger.warning(
            "Lesion type '%s' not in repair database; using default "
            "efficiency=0.50, half_life=12h", lesion_type
        )

    # --- DSB: adjust efficiency based on cell-cycle-appropriate pathway ---
    if is_dsb and not is_nhej_lesion:
        if pathway_name == "NHEJ":
            # NHEJ for G1-phase DSBs
            nhej_pw = _LESION_LOOKUP.get("dsb_g1")
            if nhej_pw is not None:
                repair_eff = nhej_pw.repair_efficiency_24h
                half_life = nhej_pw.half_life_hours
        elif pathway_name == "HR":
            # HR for S/G2-phase DSBs
            hr_pw = _LESION_LOOKUP.get("dsb_sg2")
            if hr_pw is not None:
                repair_eff = hr_pw.repair_efficiency_24h
                half_life = hr_pw.half_life_hours

    # --- Methylation adjustment ---
    # Methylated CpG deamination is poorly repaired (DRM pathway)
    formation_severity = 0.5  # default
    if "deamination" in lt_lower or "cpg" in lt_lower:
        formation_severity = 0.6
        if methylated:
            # Methylated CpG deamination repair is severely impaired
            # (Pfeifer 2006 Mutat Res)
            repair_eff = 0.30
            half_life = 24.0

    # --- TCR enhancement for NER lesions on transcribed strand ---
    if pathway_name == "NER" and transcribed_strand:
        if lt_lower in _TCR_BULKY_LESION_TYPES:
            repair_eff = adjust_for_tcr(repair_eff, transcribed_strand)
            # Use TCR half-life for CPD specifically
            if "cpd" in lt_lower:
                tcr_pw = _LESION_LOOKUP.get("tc_ner_cpd")
                if tcr_pw is not None:
                    half_life = tcr_pw.half_life_hours

    # --- Chromatin accessibility ---
    chromatin_state = "open" if chromatin_open else "heterochromatin"
    repair_eff = adjust_for_chromatin(repair_eff, chromatin_state)

    # --- Adjusted severity ---
    # Methylation can increase formation severity for CpG deamination
    adjusted_severity = formation_severity
    if methylated and ("deamination" in lt_lower or "cpg" in lt_lower):
        adjusted_severity = min(1.0, formation_severity * 1.8)

    # --- Net risk calculation ---
    net_risk = adjusted_severity * (1.0 - repair_eff)

    return RepairPrediction(
        lesion_type=lesion_type,
        position=-1,  # Position not specified in this call
        formation_severity=formation_severity,
        repair_efficiency=repair_eff,
        net_risk=net_risk,
        pathway=pathway_name,
        adjusted_severity=adjusted_severity,
    )


def compute_net_mutation_risk(
    formation_severity: float,
    lesion_type: str,
    context: dict | None = None,
    *,
    transcription_active: bool | None = None,
    chromatin_open: bool | None = None,
) -> dict:
    """Compute net mutation risk accounting for DNA repair.

    Net mutation risk = adjusted_formation_rate * (1 - repair_efficiency),
    where repair efficiency is adjusted for chromatin state, transcription,
    methylation, and cell cycle phase.

    This function provides a dictionary-based interface complementary to
    the dataclass-based ``predict_repair`` function, maintaining backward
    compatibility with the dna_damage module's ``compute_net_mutation_risk``.

    Args:
        formation_severity: Initial damage severity (0-1).
        lesion_type: Type of DNA lesion (e.g. '8-oxoG', 'CPD').
        context: Optional dict with keys:
            - 'methylated' (bool): CpG methylation status
            - 'transcribed_strand' (bool): On transcribed strand
            - 'chromatin_open' (bool): Open chromatin
            - 'chromatin_state' (str): 'open', 'heterochromatin', etc.
            - 'cell_cycle' (str): 'G1', 'S', 'G2', 'M'
            - 'position' (int): Genomic position
        transcription_active: Whether the gene is actively transcribed.
            Enables TCR boost for NER-pathway lesions. Takes precedence
            over context['transcribed_strand'] if provided.
        chromatin_open: Whether chromatin is in an open/euchromatic state.
            Takes precedence over context['chromatin_open'] if provided.

    Returns:
        Dict with keys:
            - net_risk: Net mutation risk (0-1)
            - formation_severity: Input formation severity
            - adjusted_severity: Severity after context adjustments
            - repair_efficiency: Final repair efficiency (0-1)
            - repair_pathway: Name of the primary repair pathway
            - repair_half_life_hours: Half-life of the lesion (hours)
            - lesion_type: Input lesion type
            - position: Genomic position (or -1)
            - context_applied: Dict of adjustments that were applied
    """
    if context is None:
        context = {}

    methylated: bool = context.get("methylated", False)
    # Direct keyword args take precedence over context dict
    if transcription_active is not None:
        transcribed_strand: bool = transcription_active
    else:
        transcribed_strand = context.get("transcribed_strand", False)
    if chromatin_open is not None:
        _chromatin_open: bool = chromatin_open
    else:
        _chromatin_open = context.get("chromatin_open", True)
    chromatin_state: str | None = context.get("chromatin_state", None)
    cell_cycle: str = context.get("cell_cycle", "S")
    position: int = context.get("position", -1)

    # Use predict_repair as the core engine
    prediction = predict_repair(
        lesion_type=lesion_type,
        methylated=methylated,
        transcribed_strand=transcribed_strand,
        chromatin_open=_chromatin_open,
        cell_cycle=cell_cycle,
    )

    # Override chromatin with explicit state if provided
    repair_eff = prediction.repair_efficiency
    if chromatin_state is not None:
        repair_eff = adjust_for_chromatin(repair_eff, chromatin_state)

    # Apply formation severity from input
    adjusted_severity = formation_severity
    if methylated and ("deamination" in lesion_type.lower() or "cpg" in lesion_type.lower()):
        adjusted_severity = min(1.0, formation_severity * 1.8)

    net_risk = adjusted_severity * (1.0 - repair_eff)

    # Track which adjustments were applied
    context_applied: dict[str, str | float | bool] = {}
    if methylated:
        context_applied["methylation_adjustment"] = True
    if transcribed_strand and prediction.pathway == "NER":
        context_applied["tcr_enhancement"] = True
    if chromatin_state is not None:
        context_applied["chromatin_state"] = chromatin_state
    elif not _chromatin_open:
        context_applied["chromatin_state"] = "heterochromatin"
    if "dsb" in lesion_type.lower() or "double_strand" in lesion_type.lower():
        context_applied["cell_cycle_pathway"] = prediction.pathway

    # Look up half-life from database
    pw = _LESION_LOOKUP.get(lesion_type.lower())
    half_life = pw.half_life_hours if pw is not None else 12.0

    return {
        "net_risk": net_risk,
        "formation_severity": formation_severity,
        "adjusted_severity": adjusted_severity,
        "repair_efficiency": repair_eff,
        "repair_pathway": prediction.pathway,
        "repair_half_life_hours": half_life,
        "lesion_type": lesion_type,
        "position": position,
        "context_applied": context_applied,
    }
