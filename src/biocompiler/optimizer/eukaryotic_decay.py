"""BioCompiler Eukaryotic mRNA Decay Pathway Modeling
=====================================================

Comprehensive modeling of eukaryotic mRNA decay pathways for therapeutic
mRNA design and heterologous expression optimization.

This module provides:

1. **DecayPathway** dataclass — Structured representation of individual
   decay pathways with rate modifiers, tissue specificity, and references.
2. **DecaySignal** dataclass — Detected activation signals for specific
   decay pathways with severity scoring.
3. **_DECAY_PATHWAYS** database — Eight curated eukaryotic decay pathways
   with literature-derived rate parameters.
4. **predict_decay_rate()** — Predict overall mRNA decay rate and half-life
   from sequence features, polyA tail length, and cap status.
5. **detect_decay_signals()** — Identify which decay pathways are activated
   for a given mRNA sequence and configuration.
6. **optimize_decay_resistance()** — Modify sequence codons to resist decay
   (reduce ARE creation, eliminate NMD triggers, etc.).
7. **estimate_halflife()** — Convert a first-order decay rate constant to
   half-life in hours.

Covered pathways:

- XRN1 5'→3' exonuclease
- CCR4-NOT deadenylation complex (ARE class-dependent rates)
- DCP1/DCP2 decapping (triggered by short polyA, DHH1/RCK, Pat/Lsm1-7)
- Exosome 3'→5' (DIS3/Rrp44, Rrp6)
- PARN deadenylase (polyA-specific)
- No-Go Decay (NGD) via Dom34/Hbs1
- Non-stop Decay (NSD) via Ski7
- Nonsense-Mediated Decay (NMD) via UPF1/2/3

Usage::

    from biocompiler.optimizer.eukaryotic_decay import (
        DecayPathway,
        DecaySignal,
        predict_decay_rate,
        detect_decay_signals,
        optimize_decay_resistance,
        estimate_halflife,
    )

    # Predict decay rate for a sequence
    result = predict_decay_rate(
        "ATGGCC...TAA",
        are_class="II",
        polyA_length=30,
        has_5cap=True,
    )
    print(f"Decay rate: {result['decay_rate']:.4f} /hr")
    print(f"Half-life:  {result['halflife_hours']:.2f} hr")

    # Detect activated decay signals
    signals = detect_decay_signals("ATGGCC...TAA", polyA_length=200)
    for sig in signals:
        print(f"{sig.signal_type}: {sig.description} (severity={sig.severity:.2f})")

    # Optimize sequence for decay resistance
    optimized = optimize_decay_resistance("ATGGCC...TAA")

References:
  Yamashita, A. et al. (2005). "Concerted action of the RNA helicase
  and the CCR4-NOT complex in deadenylation." *Cell* 121:527–539.
  Parker, R. (2012). "RNA degradation in Saccharomyces cerevisiae."
  *Annu Rev Biochem* 81:455–478.
  Franks, T.M. & Lykke-Andersen, J. (2010). "The control of mRNA
  decapping and P-body formation." *Mol Cell* 39:401–412.
  Doma, M.K. & Parker, R. (2006). "Endonucleolytic cleavage of
  eukaryotic mRNAs with stalls in translation elongation." *Cell*
  126:475–483.
  van Hoof, A. et al. (2002). "An exosome-dependent nuclease is
  involved in mRNA decay upon premature termination of translation."
  *Science* 296:2274–2277.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "DecayPathway",
    "DecaySignal",
    "predict_decay_rate",
    "detect_decay_signals",
    "optimize_decay_resistance",
    "estimate_halflife",
    "_DECAY_PATHWAYS",
]

# ────────────────────────────────────────────────────────────
# 1. Core Data Structures
# ────────────────────────────────────────────────────────────


@dataclass
class DecayPathway:
    """A single eukaryotic mRNA decay pathway.

    Attributes:
        name: Human-readable name of the pathway (e.g.,
            ``"XRN1 5'→3' Exonuclease"``).
        description: Detailed description of the mechanism and
            key molecular players.
        rate_modifier: Baseline rate modifier relative to default
            mRNA turnover (1.0 = no acceleration, >1.0 = faster
            decay, <1.0 = slower decay).
        tissue_specific: Whether the pathway has tissue-specific
            activity patterns.
        reference: Primary literature reference for the rate
            parameters.
    """

    name: str
    description: str
    rate_modifier: float
    tissue_specific: bool
    reference: str

    def __post_init__(self) -> None:
        if self.rate_modifier < 0:
            raise ValueError(
                f"rate_modifier must be non-negative, got {self.rate_modifier}"
            )


@dataclass
class DecaySignal:
    """A signal indicating activation of a specific decay pathway.

    Attributes:
        signal_type: Category of the decay signal (e.g.,
            ``"ccr4_not_deadenylation"``, ``"nmd_upf1"``).
        position: 0-based position in the sequence where the
            signal is detected.
        severity: Severity score in [0, 1].  Higher values
            indicate greater risk of accelerated decay.
        pathway: Decay pathway this signal activates (key in
            ``_DECAY_PATHWAYS``).
        description: Human-readable description of the signal.
    """

    signal_type: str
    position: int
    severity: float
    pathway: str
    description: str

    def __post_init__(self) -> None:
        self.severity = max(0.0, min(1.0, self.severity))


# ────────────────────────────────────────────────────────────
# 2. Comprehensive Pathway Database
# ────────────────────────────────────────────────────────────

_DECAY_PATHWAYS: dict[str, DecayPathway] = {
    "xrn1_5to3": DecayPathway(
        name="XRN1 5'→3' Exonuclease",
        description=(
            "Major 5'→3' exoribonuclease that degrades decapped mRNA. "
            "XRN1 (also known as Xrn1p/Kem1p) is the primary cytoplasmic "
            "5'→3' exoribonuclease in eukaryotes. It acts after DCP1/DCP2 "
            "removes the 5' cap structure, rapidly processively degrading "
            "the mRNA body. XRN1 is processive and highly efficient, "
            "making it the dominant decay route for most eukaryotic mRNAs."
        ),
        rate_modifier=2.0,
        tissue_specific=False,
        reference="Parker 2012 Annu Rev Biochem 81:455–478",
    ),
    "ccr4_not_deadenylation": DecayPathway(
        name="CCR4-NOT Deadenylation Complex",
        description=(
            "Major deadenylation complex that shortens polyA tails, "
            "typically the first step in mRNA decay. The CCR4-NOT complex "
            "contains two catalytic subunits: CCR4 (Pop2/Caf1) with "
            "3'→5' exonuclease activity. Deadenylation rate depends on "
            "ARE class: Class I AREs accelerate 1.5×, Class II AREs "
            "accelerate 3.0×, Class III AREs do not accelerate. "
            "Deadenylation is the rate-limiting step for most mRNA decay."
        ),
        rate_modifier=1.5,  # baseline; ARE class overrides below
        tissue_specific=False,
        reference="Yamashita 2005 Cell 121:527–539",
    ),
    "dcp1_dcp2_decapping": DecayPathway(
        name="DCP1/DCP2 Decapping",
        description=(
            "DCP2 is the catalytic subunit that hydrolyzes the 5' cap "
            "(m7GpppN), releasing m7GDP and a 5'-monophosphate mRNA that "
            "is a substrate for XRN1. DCP1 is an essential activator. "
            "Decapping is triggered by: (1) short polyA tail (<50 nt), "
            "(2) DHH1/RCK helicase recruitment after deadenylation, "
            "(3) Pat1b/Lsm1-7 complex binding to the 3' end of "
            "deadenylated mRNA. DCP1/DCP2 forms cytoplasmic foci (P-bodies) "
            "where decapping occurs."
        ),
        rate_modifier=1.0,
        tissue_specific=False,
        reference="Franks 2010 Mol Cell 39:401–412",
    ),
    "exosome_3to5": DecayPathway(
        name="Exosome 3'→5' Decay",
        description=(
            "The RNA exosome is a multi-subunit 3'→5' exoribonuclease "
            "complex that degrades RNA substrates. Cytoplasmic exosome "
            "provides a backup decay pathway when XRN1 is impaired. "
            "Catalytic subunits: DIS3/Rrp44 (processive exoribonuclease "
            "with endonuclease activity via PIN domain) and Rrp6 "
            "(nuclear exosome subunit, distributive). The exosome is "
            "slower than XRN1 but can handle structured RNA substrates "
            "that stall XRN1."
        ),
        rate_modifier=0.5,
        tissue_specific=False,
        reference="Parker 2012 Annu Rev Biochem 81:455–478",
    ),
    "parn_deadenylase": DecayPathway(
        name="PARN Deadenylase",
        description=(
            "PolyA-specific ribonuclease (PARN) is a 3'→5' exoribonuclease "
            "that specifically degrades polyA tails. PARN is cap-binding "
            "and its activity is stimulated by the presence of the 5' cap "
            "structure, creating a circular topology for deadenylation. "
            "PARN mutations cause dyskeratosis congenita via telomerase "
            "RNA misprocessing. PARN activity is particularly important "
            "in certain tissues (bone marrow, testis) and for specific "
            "RNA substrates (snoRNAs, telomerase RNA)."
        ),
        rate_modifier=0.8,
        tissue_specific=True,
        reference="Parker 2012 Annu Rev Biochem 81:455–478",
    ),
    "ngd_no_go_decay": DecayPathway(
        name="No-Go Decay (NGD)",
        description=(
            "No-Go Decay is triggered by ribosome stalling during "
            "translation elongation, typically caused by strong mRNA "
            "secondary structure, rare codons, or damaged mRNA. The "
            "stalled ribosome is recognized by Dom34 (Pelota in mammals) "
            "and Hbs1 (a GTPase related to eRF3), which mediate "
            "endonucleolytic cleavage of the mRNA near the stall site. "
            "The resulting fragments are degraded by XRN1 (5' fragment) "
            "and the exosome (3' fragment). NGD is a quality control "
            "mechanism to clear translationally stalled mRNAs."
        ),
        rate_modifier=1.5,
        tissue_specific=False,
        reference="Doma 2006 Cell 126:475–483",
    ),
    "nsd_non_stop_decay": DecayPathway(
        name="Non-stop Decay (NSD)",
        description=(
            "Non-stop Decay targets mRNAs that lack a stop codon, causing "
            "ribosomes to translate to the 3' end and stall at the polyA "
            "tail. Ski7 (a protein related to eRF3 with a GTPase domain) "
            "recognizes the stalled ribosome at the mRNA 3' end and "
            "recruits the exosome for 3'→5' degradation. In yeast, the "
            "Ski complex (Ski2, Ski3, Ski8) assists the exosome. "
            "NSD is distinct from NMD and NGD, and prevents accumulation "
            "of non-stop translation products (potentially toxic chimeric "
            "proteins from readthrough into the polyA tail)."
        ),
        rate_modifier=1.2,
        tissue_specific=False,
        reference="van Hoof 2002 Science 296:2274–2277",
    ),
    "nmd_nonsense_mediated": DecayPathway(
        name="Nonsense-Mediated Decay (NMD)",
        description=(
            "NMD targets mRNAs containing premature termination codons "
            "(PTCs) to prevent production of truncated, potentially "
            "dominant-negative proteins. The core NMD factors are UPF1 "
            "(a helicase/ATPase that is the central regulator), UPF2 "
            "(an adaptor protein), and UPF3/UPF3X (nuclear-cytoplasmic "
            "shuttling factor that bridges EJC to UPF2). Two triggers: "
            "(1) EJC-dependent: PTC ≥50-55 nt upstream of the last "
            "exon-exon junction; (2) Long 3'UTR: 3'UTR >2000 nt in "
            "intronless transcripts. UPF1 phosphorylation by SMG1 "
            "initiates decay via SMG5/6/7 recruitment and subsequent "
            "deadenylation, decapping, and exonucleolysis."
        ),
        rate_modifier=3.0,
        tissue_specific=False,
        reference="Parker 2012 Annu Rev Biochem 81:455–478",
    ),
}

# ARE class-dependent rate modifiers for the CCR4-NOT complex
# Following Yamashita 2005 Cell and Chen & Shyu 1995 classification
_ARE_CLASS_RATES: dict[str, float] = {
    "I": 1.5,   # Class I: scattered AUUUA in U-rich context — moderate acceleration
    "II": 3.0,  # Class II: overlapping AUUUA repeats — strong acceleration
    "III": 1.0, # Class III: U-rich without AUUUA — no acceleration
    "none": 1.0,
}

# ────────────────────────────────────────────────────────────
# 3. ARE Detection Patterns
# ────────────────────────────────────────────────────────────

# ARE classification follows Chen & Shyu (1995):
#   Class I: scattered AUUUA in U-rich context
#   Class II: overlapping AUUUA repeats (rapid decay)
#   Class III: U-rich without AUUUA

_ARE_PATTERNS: list[tuple[str, str, str, float]] = [
    # (pattern, are_class, description, base_severity)
    (r"ATTTA", "I", "ARE Class I (scattered AUUUA in U-rich context)", 0.3),
    (r"ATTTATTTA", "II", "ARE Class II (overlapping AUUUA repeats)", 0.7),
    (r"ATTTATTTATTTA", "II", "ARE Class II (triple AUUUA repeat — rapid decay)", 0.9),
    (r"T{6,}", "III", "ARE Class III (U-rich element without AUUUA)", 0.2),
]

# ────────────────────────────────────────────────────────────
# 4. Baseline decay parameters
# ────────────────────────────────────────────────────────────

# Baseline mRNA half-life in mammalian cells for a typical stable transcript
# with full cap and long polyA tail: ~10 hours (Schwanhausser 2011 Nature)
_BASELINE_HALFLIFE_HOURS: float = 10.0

# Baseline decay rate (first-order): k = ln(2) / t_half
_BASELINE_DECAY_RATE: float = math.log(2) / _BASELINE_HALFLIFE_HOURS

# Thresholds for pathway activation
_POLYA_SHORT_THRESHOLD: int = 50      # DCP1/DCP2 triggered below this
_POLYA_PARTIAL_THRESHOLD: int = 80    # XRN1 activated below this
_POLYA_PARN_OPTIMAL: int = 150        # PARN most active above this
_3UTR_LONG_NMD_THRESHOLD: int = 2000  # Long 3'UTR NMD trigger
_PTC_UPSTREAM_DISTANCE: int = 50      # PTC distance from last EJ for NMD

# NGD stall patterns: strong secondary structure indicators and poly-basic runs
_NGD_STALL_PATTERNS: list[tuple[str, str, float]] = [
    # (pattern, description, base_severity)
    (r"G{4,}", "G-quadruplex-forming sequence (stalls ribosomes)", 0.6),
    (r"C{4,}", "Poly-C stretch (strong secondary structure)", 0.4),
    (r"([AG]{3,})\1{2,}", "Triplet repeat (can stall ribosomes)", 0.5),
]

# Codon table for synonym-based optimization (standard genetic code)
_CODON_TABLE: dict[str, list[str]] = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "*": ["TAA", "TAG", "TGA"],
}

# Reverse lookup: codon → amino acid
_CODON_TO_AA: dict[str, str] = {}
for _aa, _codons in _CODON_TABLE.items():
    for _codon in _codons:
        _CODON_TO_AA[_codon] = _aa


# ────────────────────────────────────────────────────────────
# 5. Public API Functions
# ────────────────────────────────────────────────────────────


def estimate_halflife(decay_rate: float) -> float:
    """Convert a first-order decay rate constant to half-life in hours.

    Uses the relationship: t_1/2 = ln(2) / k

    Args:
        decay_rate: First-order decay rate constant (per hour).
            Must be non-negative.  A decay rate of 0 returns
            ``float('inf')`` (stable transcript).

    Returns:
        Half-life in hours.  Returns ``float('inf')`` for a
        decay rate of 0.

    Raises:
        ValueError: If decay_rate is negative.
    """
    if decay_rate < 0:
        raise ValueError(f"decay_rate must be non-negative, got {decay_rate}")
    if decay_rate == 0:
        return float("inf")
    return math.log(2) / decay_rate


def detect_decay_signals(
    seq: str,
    are_signals: list[DecaySignal] | None = None,
    polyA_length: int = 200,
) -> list[DecaySignal]:
    """Identify which decay pathways are activated for a given mRNA.

    Scans the sequence and evaluates cap/PolyA status to determine
    which eukaryotic decay pathways are active, and with what severity.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        are_signals: Pre-detected ARE signals.  If ``None``, ARE
            detection is run automatically.
        polyA_length: Length of polyA tail in nucleotides
            (default 200, typical for newly transcribed mRNA).

    Returns:
        List of :class:`DecaySignal` objects for activated pathways,
        sorted by severity (highest first).
    """
    signals: list[DecaySignal] = []
    seq_upper = seq.upper()
    seq_len = len(seq_upper)

    # ── ARE detection ──
    if are_signals is None:
        are_signals = _detect_are_signals_internal(seq_upper)

    # Determine dominant ARE class
    are_class = _determine_are_class(are_signals)

    # ── CCR4-NOT deadenylation ──
    ccr4_rate = _ARE_CLASS_RATES.get(are_class, 1.0)
    if ccr4_rate > 1.0:
        signals.append(DecaySignal(
            signal_type="ccr4_not_deadenylation",
            position=0,
            severity=min(1.0, ccr4_rate / 3.0),
            pathway="ccr4_not_deadenylation",
            description=(
                f"CCR4-NOT deadenylation accelerated {ccr4_rate:.1f}× "
                f"by ARE Class {are_class}"
            ),
        ))

    # ── PARN deadenylase ──
    if polyA_length > _POLYA_PARN_OPTIMAL:
        # PARN is active on long polyA tails with cap present
        parn_severity = 0.15  # low baseline — PARN is a minor pathway
        signals.append(DecaySignal(
            signal_type="parn_deadenylase",
            position=seq_len - 1,
            severity=parn_severity,
            pathway="parn_deadenylase",
            description=(
                f"PARN deadenylase active on polyA tail "
                f"({polyA_length} nt)"
            ),
        ))
    elif polyA_length < _POLYA_SHORT_THRESHOLD:
        # Short polyA reduces PARN activity but other pathways dominate
        pass

    # ── DCP1/DCP2 decapping ──
    # Triggered by: short polyA, DHH1/RCK, Pat/Lsm1-7
    if polyA_length < _POLYA_SHORT_THRESHOLD:
        cap_severity = max(0.0, 1.0 - polyA_length / _POLYA_SHORT_THRESHOLD)
        # Additional triggers: DHH1/RCK and Pat/Lsm1-7 are recruited
        # after deadenylation — model as additional severity when
        # CCR4-NOT has been active
        dhh1_boost = 0.1 if ccr4_rate > 1.0 else 0.0
        pat_lsm_boost = 0.15 if polyA_length < 20 else 0.0
        signals.append(DecaySignal(
            signal_type="dcp1_dcp2_decapping",
            position=0,
            severity=min(1.0, cap_severity + dhh1_boost + pat_lsm_boost),
            pathway="dcp1_dcp2_decapping",
            description=(
                f"DCP1/DCP2 decapping triggered by short polyA "
                f"({polyA_length} nt)"
                + (" + DHH1/RCK recruitment" if dhh1_boost > 0 else "")
                + (" + Pat/Lsm1-7 binding" if pat_lsm_boost > 0 else "")
            ),
        ))

    # ── XRN1 5'→3' exonuclease ──
    # Active on decapped mRNA; partially deadenylated mRNA more susceptible
    if polyA_length < _POLYA_PARTIAL_THRESHOLD:
        xrn1_severity = min(1.0, 0.3 + 0.7 * (1.0 - polyA_length / _POLYA_PARTIAL_THRESHOLD))
        # Boost if decapping has been triggered
        decapping_active = any(s.signal_type == "dcp1_dcp2_decapping" for s in signals)
        if decapping_active:
            xrn1_severity = min(1.0, xrn1_severity + 0.2)
        signals.append(DecaySignal(
            signal_type="xrn1_5to3",
            position=0,
            severity=xrn1_severity,
            pathway="xrn1_5to3",
            description=(
                f"XRN1 5'→3' decay active"
                + (" (decapping-dependent)" if decapping_active else "")
                + f" (polyA={polyA_length} nt)"
            ),
        ))

    # ── Exosome 3'→5' ──
    # Always present as backup; more active when XRN1 is impaired
    exosome_severity = 0.15  # baseline backup
    xrn1_active = any(s.signal_type == "xrn1_5to3" for s in signals)
    if not xrn1_active:
        exosome_severity = 0.3  # compensates when XRN1 not active
    # DIS3/Rrp44 and Rrp6 can degrade structured substrates
    signals.append(DecaySignal(
        signal_type="exosome_3to5",
        position=seq_len - 1,
        severity=exosome_severity,
        pathway="exosome_3to5",
        description=(
            f"Exosome 3'→5' decay (backup pathway)"
            + (" (compensating for low XRN1)" if not xrn1_active else "")
        ),
    ))

    # ── No-Go Decay (NGD) ──
    ngd_signals = _detect_ngd_signals(seq_upper)
    signals.extend(ngd_signals)

    # ── Non-stop Decay (NSD) ──
    nsd_signals = _detect_nsd_signals(seq_upper)
    signals.extend(nsd_signals)

    # ── Nonsense-Mediated Decay (NMD) ──
    nmd_signals = _detect_nmd_signals_internal(seq_upper)
    signals.extend(nmd_signals)

    # Sort by severity (highest first)
    signals.sort(key=lambda s: s.severity, reverse=True)

    return signals


def predict_decay_rate(
    seq: str,
    are_class: str,
    polyA_length: int,
    has_5cap: bool = True,
    is_transcribed: bool = True,
) -> dict[str, Any]:
    """Predict overall mRNA decay rate and half-life.

    Combines contributions from all active decay pathways into a
    composite first-order decay rate constant, then computes the
    corresponding half-life.

    The model uses additive rate contributions from each pathway,
    weighted by the pathway's rate modifier and severity of activation:

        k_total = k_baseline × (1 + Σ severity_i × (modifier_i - 1))

    where the sum is over all active pathways with modifier > 1.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        are_class: ARE class detected in the sequence
            (``"I"``, ``"II"``, ``"III"``, or ``"none"``).
        polyA_length: Length of polyA tail in nucleotides.
        has_5cap: Whether the mRNA has an intact 5' cap
            (default ``True``).  Uncapped mRNA is rapidly degraded
            by XRN1.
        is_transcribed: Whether the mRNA is being actively
            transcribed (default ``True``).  Affects transcription-
            coupled repair and nuclear surveillance.

    Returns:
        Dictionary with keys:

        - ``"decay_rate"`` — Composite first-order decay rate (per hour)
        - ``"halflife_hours"`` — Predicted half-life in hours
        - ``"pathway_contributions"`` — Dict mapping pathway key to
          its rate contribution
        - ``"dominant_pathway"`` — Key of the pathway with the
          largest contribution
        - ``"risk_level"`` — Overall risk level (``"low"``,
          ``"moderate"``, ``"high"``, or ``"very_high"``)
        - ``"signals"`` — List of :class:`DecaySignal` objects
          for activated pathways
    """
    seq_upper = seq.upper()

    # Detect signals
    signals = detect_decay_signals(seq_upper, polyA_length=polyA_length)

    # Determine ARE class from signals if not provided explicitly
    if are_class == "none":
        are_class = _determine_are_class(signals)

    # Compute pathway contributions
    contributions: dict[str, float] = {}

    # Start with baseline rate
    k_total = _BASELINE_DECAY_RATE
    contributions["baseline"] = _BASELINE_DECAY_RATE

    # Add contributions from each active pathway
    for signal in signals:
        pathway_key = signal.pathway
        pathway = _DECAY_PATHWAYS.get(pathway_key)
        if pathway is None:
            continue

        # CCR4-NOT rate depends on ARE class
        if pathway_key == "ccr4_not_deadenylation":
            effective_modifier = _ARE_CLASS_RATES.get(are_class, 1.0)
        else:
            effective_modifier = pathway.rate_modifier

        # Contribution = severity × (modifier - 1) × baseline_rate
        # Only pathways with modifier > 1 accelerate decay
        if effective_modifier > 1.0:
            contribution = signal.severity * (effective_modifier - 1.0) * _BASELINE_DECAY_RATE
        else:
            # Slowdown pathways reduce the total rate
            contribution = -signal.severity * (1.0 - effective_modifier) * _BASELINE_DECAY_RATE * 0.1

        contributions[pathway_key] = contribution
        k_total += contribution

    # ── No-5'cap penalty ──
    # Uncapped mRNA is rapidly degraded by XRN1 even without decapping
    if not has_5cap:
        cap_penalty = 2.0 * _BASELINE_DECAY_RATE
        contributions["no_5cap"] = cap_penalty
        k_total += cap_penalty

    # ── Not transcribed penalty (nuclear surveillance) ──
    # mRNA not being actively transcribed is subject to nuclear
    # surveillance and exosome degradation
    if not is_transcribed:
        transcription_penalty = 0.5 * _BASELINE_DECAY_RATE
        contributions["not_transcribed"] = transcription_penalty
        k_total += transcription_penalty

    # Ensure non-negative rate
    k_total = max(0.0, k_total)

    # Compute half-life
    halflife = estimate_halflife(k_total)

    # Determine dominant pathway
    non_baseline = {k: v for k, v in contributions.items() if k != "baseline"}
    dominant = max(non_baseline, key=lambda k: abs(non_baseline[k])) if non_baseline else "baseline"

    # Risk level
    if halflife > 8.0:
        risk_level = "low"
    elif halflife > 4.0:
        risk_level = "moderate"
    elif halflife > 2.0:
        risk_level = "high"
    else:
        risk_level = "very_high"

    return {
        "decay_rate": k_total,
        "halflife_hours": halflife,
        "pathway_contributions": contributions,
        "dominant_pathway": dominant,
        "risk_level": risk_level,
        "signals": signals,
    }


def optimize_decay_resistance(
    seq: str,
    are_signals: list[DecaySignal] | None = None,
) -> str:
    """Modify sequence to resist decay pathways.

    Applies codon-level optimizations to reduce the activation of
    decay pathways without changing the encoded protein:

    1. **ARE reduction** — Replace codons creating ATTTA motifs with
       synonymous codons that avoid AUUUA patterns.
    2. **NGD reduction** — Replace codons creating G-quadruplex-
       forming G-runs and strong stall signals.
    3. **NSD reduction** — Ensure a proper stop codon exists.
    4. **NMD reduction** — Reduce premature stop codons by codon
       substitution (when not at the intended ORF end).

    The algorithm iterates codon-by-codon and, for each problem
    motif, picks a synonymous codon that does not recreate the
    motif.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U). Must be
            a valid coding sequence (length divisible by 3).
        are_signals: Pre-detected ARE signals.  If ``None``,
            ARE detection is run automatically.

    Returns:
        Optimized mRNA sequence with reduced decay signals.

    Raises:
        ValueError: If the sequence length is not divisible by 3
            or contains invalid codons.
    """
    seq_upper = seq.upper()

    if len(seq_upper) % 3 != 0:
        raise ValueError(
            f"Sequence length {len(seq_upper)} is not divisible by 3. "
            f"optimize_decay_resistance requires a valid coding sequence."
        )

    # Validate codons
    codons = [seq_upper[i:i+3] for i in range(0, len(seq_upper), 3)]
    for idx, codon in enumerate(codons):
        if codon not in _CODON_TO_AA:
            raise ValueError(
                f"Invalid codon '{codon}' at position {idx * 3}"
            )

    # Detect ARE signals if not provided
    if are_signals is None:
        are_signals = _detect_are_signals_internal(seq_upper)

    # Collect positions of problem motifs
    problem_positions: set[int] = set()
    for sig in are_signals:
        # Mark the codons overlapping with this signal
        for pos in range(sig.position, min(sig.position + 8, len(seq_upper))):
            problem_positions.add(pos)

    # Also detect NGD stall signals
    for pattern, _desc, _sev in _NGD_STALL_PATTERNS:
        for match in re.finditer(pattern, seq_upper):
            for pos in range(match.start(), match.end()):
                problem_positions.add(pos)

    # Also detect NMD premature stops
    nmd_signals = _detect_nmd_signals_internal(seq_upper)
    for sig in nmd_signals:
        for pos in range(sig.position, min(sig.position + 3, len(seq_upper))):
            problem_positions.add(pos)

    # Build optimized sequence codon by codon
    optimized_codons = list(codons)

    for codon_idx in range(len(codons)):
        codon_start = codon_idx * 3
        codon_end = codon_start + 3
        # Check if this codon overlaps with any problem position
        overlaps_problem = any(
            codon_start <= p < codon_end for p in problem_positions
        )

        if not overlaps_problem:
            continue

        current_codon = optimized_codons[codon_idx]
        aa = _CODON_TO_AA[current_codon]

        # Skip methionine and tryptophan (no synonyms)
        if aa in ("M", "W"):
            continue

        # Skip the final stop codon (intentional)
        if aa == "*" and codon_idx == len(codons) - 1:
            continue

        # Try synonymous codons
        synonyms = _CODON_TABLE[aa]
        best_codon = current_codon
        best_score = _score_codon_context(
            optimized_codons, codon_idx, current_codon, problem_positions
        )

        for syn in synonyms:
            if syn == current_codon:
                continue
            score = _score_codon_context(
                optimized_codons, codon_idx, syn, problem_positions
            )
            if score < best_score:
                best_codon = syn
                best_score = score

        optimized_codons[codon_idx] = best_codon

    # Ensure the sequence ends with a stop codon (NSD resistance)
    if len(optimized_codons) > 0:
        last_aa = _CODON_TO_AA[optimized_codons[-1]]
        if last_aa != "*":
            # Append a stop codon
            optimized_codons.append("TAA")

    return "".join(optimized_codons)


# ────────────────────────────────────────────────────────────
# 6. Internal Helper Functions
# ────────────────────────────────────────────────────────────


def _detect_are_signals_internal(seq_upper: str) -> list[DecaySignal]:
    """Internal ARE detection returning DecaySignal objects."""
    signals: list[DecaySignal] = []

    for pattern, are_class, description, base_severity in _ARE_PATTERNS:
        for match in re.finditer(pattern, seq_upper):
            # Position-dependent severity: AREs in 3'UTR are more potent
            frac = match.start() / max(len(seq_upper), 1)
            if frac > 0.85:
                position_mult = 1.5  # 3'UTR
            elif frac > 0.75:
                position_mult = 1.2  # Near 3'UTR
            else:
                position_mult = 1.0  # CDS/5'UTR

            severity = min(1.0, base_severity * position_mult)
            signals.append(DecaySignal(
                signal_type=f"are_class_{are_class.lower()}",
                position=match.start(),
                severity=severity,
                pathway="ccr4_not_deadenylation",
                description=f"{description} at position {match.start()}: {match.group()}",
            ))

    return signals


def _determine_are_class(signals: list[DecaySignal]) -> str:
    """Determine the dominant ARE class from a list of decay signals."""
    if any(s.signal_type == "are_class_ii" for s in signals):
        return "II"
    elif any(s.signal_type == "are_class_i" for s in signals):
        return "I"
    elif any(s.signal_type == "are_class_iii" for s in signals):
        return "III"
    return "none"


def _detect_ngd_signals(seq_upper: str) -> list[DecaySignal]:
    """Detect No-Go Decay (NGD) activation signals.

    NGD is triggered by ribosome stalls from:
    - G-quadruplex forming sequences
    - Strong secondary structure (poly-C)
    - Triplet repeats
    """
    signals: list[DecaySignal] = []

    for pattern, description, base_severity in _NGD_STALL_PATTERNS:
        for match in re.finditer(pattern, seq_upper):
            signals.append(DecaySignal(
                signal_type="ngd_stall",
                position=match.start(),
                severity=base_severity,
                pathway="ngd_no_go_decay",
                description=(
                    f"NGD trigger: {description} at position "
                    f"{match.start()}: {match.group()}"
                ),
            ))

    return signals


def _detect_nsd_signals(seq_upper: str) -> list[DecaySignal]:
    """Detect Non-stop Decay (NSD) activation signals.

    NSD is triggered when an mRNA lacks a stop codon or when
    the stop codon is at the very end without adequate 3'UTR.
    """
    signals: list[DecaySignal] = []

    # Check if the sequence ends with a stop codon
    if len(seq_upper) >= 3:
        last_codon = seq_upper[-3:]
        if last_codon not in ("TAA", "TAG", "TGA"):
            # No stop codon at the end → NSD trigger
            signals.append(DecaySignal(
                signal_type="nsd_no_stop",
                position=len(seq_upper) - 3,
                severity=0.8,
                pathway="nsd_non_stop_decay",
                description=(
                    f"NSD trigger: no stop codon at 3' end "
                    f"(last codon: {last_codon})"
                ),
            ))

    # Check for stop codons followed by very short 3'UTR
    # (< 20 nt after last stop codon → borderline NSD risk)
    for i in range(len(seq_upper) - 3, -1, -3):
        codon = seq_upper[i:i+3]
        if codon in ("TAA", "TAG", "TGA"):
            trailing = len(seq_upper) - (i + 3)
            if 0 < trailing < 20:
                signals.append(DecaySignal(
                    signal_type="nsd_short_3utr",
                    position=i,
                    severity=0.3,
                    pathway="nsd_non_stop_decay",
                    description=(
                        f"NSD risk: very short 3'UTR ({trailing} nt) "
                        f"after stop codon at position {i}"
                    ),
                ))
            break

    return signals


def _detect_nmd_signals_internal(seq_upper: str) -> list[DecaySignal]:
    """Detect NMD activation signals.

    For synthetic/intronless mRNA:
    - Long 3'UTR (>2000 nt) triggers NMD via UPF1/2/3

    For mRNA with ORF analysis:
    - Premature stop codons (in-frame stops before the final stop)
    """
    signals: list[DecaySignal] = []

    # Find the first ORF (start with ATG, end with stop)
    orf_start = -1
    orf_end = -1
    for i in range(0, len(seq_upper) - 2, 3):
        if seq_upper[i:i+3] == "ATG":
            orf_start = i
            break

    if orf_start < 0:
        return signals

    # Find the first in-frame stop codon
    for i in range(orf_start, len(seq_upper) - 2, 3):
        codon = seq_upper[i:i+3]
        if codon in ("TAA", "TAG", "TGA"):
            orf_end = i + 3
            break

    if orf_end < 0:
        # No stop codon found → NSD, not NMD
        return signals

    # Check for premature stop codons (PTCs) before the expected end
    # Scan for additional in-frame stops before the last codon
    ptc_positions: list[int] = []
    for i in range(orf_start + 3, orf_end - 3, 3):
        codon = seq_upper[i:i+3]
        if codon in ("TAA", "TAG", "TGA"):
            ptc_positions.append(i)

    # If there are PTCs, it's likely a frameshift issue or real PTC
    # For NMD: check if PTC is >50 nt upstream of the final stop
    for ptc_pos in ptc_positions:
        distance_to_end = orf_end - ptc_pos
        if distance_to_end > _PTC_UPSTREAM_DISTANCE:
            severity = min(1.0, 0.4 + 0.6 * (distance_to_end / 1000))
            signals.append(DecaySignal(
                signal_type="nmd_ptc",
                position=ptc_pos,
                severity=severity,
                pathway="nmd_nonsense_mediated",
                description=(
                    f"NMD trigger: premature stop codon at position "
                    f"{ptc_pos} ({seq_upper[ptc_pos:ptc_pos+3]}), "
                    f"{distance_to_end} nt upstream of final stop"
                ),
            ))

    # Check 3'UTR length for long 3'UTR NMD
    utr3_length = len(seq_upper) - orf_end
    if utr3_length > _3UTR_LONG_NMD_THRESHOLD:
        severity = min(1.0, utr3_length / 5000)
        signals.append(DecaySignal(
            signal_type="nmd_long_3utr",
            position=orf_end,
            severity=severity,
            pathway="nmd_nonsense_mediated",
            description=(
                f"NMD trigger: long 3'UTR ({utr3_length} nt, "
                f"threshold: {_3UTR_LONG_NMD_THRESHOLD} nt)"
            ),
        ))

    return signals


def _score_codon_context(
    codons: list[str],
    codon_idx: int,
    candidate: str,
    problem_positions: set[int],
) -> float:
    """Score a candidate codon substitution by counting problem motifs.

    Lower score = better (fewer remaining problem motifs).

    Args:
        codons: Current list of codons (will be temporarily modified).
        codon_idx: Index of the codon being substituted.
        candidate: Candidate replacement codon.
        problem_positions: Set of positions with problem motifs.

    Returns:
        Score (number of problem positions still covered after
        substitution).  Lower is better.
    """
    # Temporarily substitute
    original = codons[codon_idx]
    codons[codon_idx] = candidate

    # Reconstruct sequence for the local window
    codon_start = codon_idx * 3
    # Check a window around the substitution for new problem motifs
    window_start = max(0, codon_start - 12)
    window_end = min(len(codons) * 3, codon_start + 15)
    window_seq = "".join(codons)[window_start:window_end]

    # Count remaining problem positions in the window
    score = 0
    for p in problem_positions:
        if window_start <= p < window_end:
            score += 1

    # Check for new ARE motifs introduced by the substitution
    for pattern, _, _, base_sev in _ARE_PATTERNS:
        for match in re.finditer(pattern, window_seq):
            abs_pos = window_start + match.start()
            if abs_pos not in problem_positions:
                score += 2  # Penalize new ARE creation

    # Check for new NGD stall motifs
    for pattern, _, base_sev in _NGD_STALL_PATTERNS:
        for match in re.finditer(pattern, window_seq):
            abs_pos = window_start + match.start()
            if abs_pos not in problem_positions:
                score += 2  # Penalize new stall creation

    # Check for new premature stop codons
    for i in range(0, len(window_seq) - 2, 3):
        codon = window_seq[i:i+3]
        if codon in ("TAA", "TAG", "TGA"):
            abs_pos = window_start + i
            if abs_pos not in problem_positions:
                score += 3  # Strongly penalize new stop codon creation

    # Restore original
    codons[codon_idx] = original

    return score
