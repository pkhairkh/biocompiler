"""BioCompiler Epitranscriptomic Modification Site Detection Module
==================================================================

Comprehensive detection and scoring of epitranscriptomic modification
sites in mRNA for therapeutic mRNA design and heterologous expression.

This module provides standalone detection of the six major mRNA
epitranscriptomic marks with enzyme-specific consensus motifs:

1. **m6A (N6-methyladenosine)** — DRACH and RRACH consensus motifs,
   catalysed by METTL3/METTL14 complex.
2. **m5C (5-methylcytosine)** — DNMT2, NSUN1-7 enzyme-specific motifs.
3. **Ψ (Pseudouridine)** — Pus1, Pus4, and Dyskerin/H/ACA snoRNA-guided
   motifs.
4. **m1A (N1-methyladenosine)** — TRMT6/61A consensus and tRNA T-loop
   motifs.
5. **2'-O-Me (2'-O-methylation)** — Fibrillarin/snoRNA-guided motifs.
6. **m6Am (N6,2'-O-dimethyladenosine)** — BCAA motif at transcription
   start, targeted by FTO demethylase.

Each detection function returns a list of :class:`EpitranscriptomicSite`
objects with modification type, position, sequence context, severity,
catalysing enzyme, and literature reference.

Impact assessment functions relate detected marks to mRNA stability,
translation, and immune evasion outcomes.

Usage::

    from biocompiler.optimizer.epitranscriptomics import (
        EpitranscriptomicSite,
        detect_m6a_sites,
        detect_m5c_sites,
        detect_pseudouridine_sites,
        detect_m1a_sites,
        detect_2om_sites,
        detect_m6am_sites,
        detect_all_epitranscriptomic_marks,
        assess_stability_impact,
        get_modification_function,
    )

    # Detect all marks
    sites = detect_all_epitranscriptomic_marks("ATGGCC...TAA")

    # Or select specific marks
    sites = detect_all_epitranscriptomic_marks("ATGGCC...TAA",
                                                marks=["m6a", "m5c"])

    # Assess stability impact
    impact = assess_stability_impact(sites)

References:
  Dominissini, D. et al. (2012). "Topology of the human and mouse m6A
  RNA methylomes revealed by m6A-seq." *Nature* 485:201-206.
  Meyer, K.D. et al. (2012). "Comprehensive analysis of mRNA
  methylation reveals enrichment in 3' UTRs and near stop codons."
  *Cell* 149:1635-1646.
  Squires, J.E. et al. (2012). "Widespread occurrence of 5-
  methylcytosine in human coding and non-coding RNA." *Nucleic Acids
  Research* 40:5023-5033.
  Carlile, T.M. et al. (2014). "Pseudouridine profiling reveals
  regulated mRNA pseudouridylation in yeast and human cells."
  *Nature* 515:143-146.
  Dominissini, D. et al. (2016). "The dynamic N1-methyladenosine
  methylome in eukaryotic messenger RNA." *Nature* 530:441-446.
  Wei, C.M. et al. (1976). "Methylated nucleotides block 5' terminus
  of HeLa cell messenger RNA." *Cell* 7:1-10. (m6Am)
  Mauer, J. et al. (2017). "Reversible methylation of m6Am in the
  5' cap controls mRNA stability." *Nature* 541:371-375. (FTO/m6Am)
"""

from __future__ import annotations

import re
import logging
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "EpitranscriptomicSite",
    "detect_m6a_sites",
    "detect_m5c_sites",
    "detect_pseudouridine_sites",
    "detect_m1a_sites",
    "detect_2om_sites",
    "detect_m6am_sites",
    "detect_all_epitranscriptomic_marks",
    "assess_stability_impact",
    "get_modification_function",
    "MODIFICATION_FUNCTIONS",
]

# ────────────────────────────────────────────────────────────
# 1. Core Data Structures
# ────────────────────────────────────────────────────────────


@dataclass
class EpitranscriptomicSite:
    """A single epitranscriptomic modification site detected in mRNA.

    Attributes:
        modification_type: Type of modification (e.g., ``"m6A"``,
            ``"m5C"``, ``"Psi"``, ``"m1A"``, ``"2OMe"``, ``"m6Am"``).
        position: 0-based position in the sequence where the
            modification is predicted.
        sequence_context: The matched subsequence surrounding the
            modification site.
        severity: Severity score in ``[0, 1]``.  Higher values indicate
            greater predicted functional impact on mRNA fate.
        enzyme: The writer enzyme or enzyme complex responsible for the
            modification (e.g., ``"METTL3/METTL14"``, ``"NSUN2"``,
            ``"Pus1"``).
        reference: Literature reference for the modification motif
            (e.g., ``"Dominissini et al. 2012 Nature"``).
    """

    modification_type: str
    position: int
    sequence_context: str
    severity: float
    enzyme: str
    reference: str

    def __post_init__(self) -> None:
        self.severity = max(0.0, min(1.0, self.severity))


# ────────────────────────────────────────────────────────────
# 2. Modification Functional Effects Database
# ────────────────────────────────────────────────────────────

MODIFICATION_FUNCTIONS: dict[str, dict[str, Any]] = {
    "m6A": {
        "effects": [
            "Promotes decay via YTHDF2 recruitment",
            "Enhances translation via YTHDF1 recruitment",
            "Facilitates nuclear export via YTHDC1",
        ],
        "stability_impact": -0.3,  # net destabilising
        "translation_impact": 0.2,
        "key_readers": ["YTHDF1", "YTHDF2", "YTHDC1"],
        "writer": "METTL3/METTL14",
        "eraser": "FTO, ALKBH5",
        "reference": "Dominissini et al. 2012 Nature; Meyer et al. 2012 Cell",
    },
    "m5C": {
        "effects": [
            "Promotes mRNA stability via ALYREF binding",
            "Facilitates nuclear export",
        ],
        "stability_impact": 0.3,  # net stabilising
        "translation_impact": 0.05,
        "key_readers": ["ALYREF"],
        "writer": "NSUN2, DNMT2",
        "eraser": "TET2 (potential)",
        "reference": "Squires et al. 2012 NAR",
    },
    "Psi": {
        "effects": [
            "Enhances mRNA stability",
            "Alters codon-anticodon pairing properties",
        ],
        "stability_impact": 0.4,  # net stabilising
        "translation_impact": 0.15,
        "key_readers": [],
        "writer": "Pus1, Pus4, Dyskerin/Cbf5",
        "eraser": "None known",
        "reference": "Carlile et al. 2014 Nature",
    },
    "m1A": {
        "effects": [
            "Promotes translation initiation",
            "Enhances ribosome interaction",
        ],
        "stability_impact": 0.1,
        "translation_impact": 0.35,
        "key_readers": ["YTHDF3 (potential)"],
        "writer": "TRMT6/61A",
        "eraser": "ALKBH3",
        "reference": "Dominissini et al. 2016 Nature",
    },
    "2OMe": {
        "effects": [
            "Enhances mRNA stability",
            "Innate immune evasion (prevents RIG-I/MDA5 sensing)",
        ],
        "stability_impact": 0.35,  # net stabilising
        "translation_impact": 0.1,
        "key_readers": [],
        "writer": "Fibrillarin/snoRNA",
        "eraser": "None known",
        "reference": "Karijolich et al. 2015 RNA",
    },
    "m6Am": {
        "effects": [
            "Promotes mRNA stability (FTO demethylase target)",
            "Cap-proximal modification that resists DCP2 decapping",
        ],
        "stability_impact": 0.5,  # strongly stabilising
        "translation_impact": 0.15,
        "key_readers": [],
        "writer": "PCIF1/CAPAM",
        "eraser": "FTO",
        "reference": "Mauer et al. 2017 Nature; Wei et al. 1976 Cell",
    },
}


def get_modification_function(mod_type: str) -> str:
    """Return the known functional effect description for a modification type.

    Args:
        mod_type: Modification type string (e.g., ``"m6A"``, ``"m5C"``,
            ``"Psi"``, ``"m1A"``, ``"2OMe"``, ``"m6Am"``).

    Returns:
        Human-readable description of the modification's functional
        effects, or a message indicating the type is unknown.
    """
    info = MODIFICATION_FUNCTIONS.get(mod_type)
    if info is None:
        return f"Unknown modification type: {mod_type}"
    effects = info["effects"]
    stability = info["stability_impact"]
    if stability > 0:
        stability_word = "stabilising"
    elif stability < 0:
        stability_word = "destabilising"
    else:
        stability_word = "neutral"
    return (
        f"{mod_type}: {'; '.join(effects)} "
        f"(net {stability_word}, impact={stability:+.1f}). "
        f"Writer: {info['writer']}; Eraser: {info['eraser']}. "
        f"Ref: {info['reference']}"
    )


# ────────────────────────────────────────────────────────────
# 3. m6A (N6-Methyladenosine) Detection
# ────────────────────────────────────────────────────────────

# DRACH consensus: D=A/G/U  R=A/G  H=A/C/U
# In DNA alphabet: [AG][AG]AC[ACT]
_M6A_DRACH = re.compile(r"[AG][AG]AC[ACT]")

# RRACH extended consensus: R=A/G (two positions) then A then C then H=A/C/U
# In DNA alphabet: [AG][AG]AC[ACT]  — same as DRACH because D includes G/A
# and the first R maps to [AG] already.  However, the *extended* RRACH
# consensus includes an additional purine 5' of the DRACH core:
# [AG][AG][AG]AC[ACT] — i.e. RRACH with a leading R.
_M6A_RRACH = re.compile(r"[AG][AG][AG]AC[ACT]")


def detect_m6a_sites(seq: str) -> list[EpitranscriptomicSite]:
    """Detect m6A (N6-methyladenosine) consensus motifs in mRNA.

    Scans for both the DRACH consensus (D=A/G/U, R=A/G, H=A/C/U)
    and the extended RRACH motif.  m6A modifications near stop codons
    and in 3'UTRs are particularly impactful for stability
    (Dominissini et al. 2012; Meyer et al. 2012).

    DRACH matches that are also part of an RRACH match receive a
    higher severity score due to stronger consensus support.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of :class:`EpitranscriptomicSite` objects for predicted
        m6A sites.
    """
    sites: list[EpitranscriptomicSite] = []
    seq_upper = seq.upper()
    n = len(seq_upper)

    # Collect RRACH positions for severity boosting
    rrach_positions: set[int] = set()
    for match in _M6A_RRACH.finditer(seq_upper):
        # The modified A is at position match.start() + 3 (0-based)
        rrach_positions.add(match.start() + 3)

    for match in _M6A_DRACH.finditer(seq_upper):
        # The modified A is at position match.start() + 2 (0-based)
        mod_pos = match.start() + 2
        frac = match.start() / max(n, 1)

        # Position-dependent severity:
        # m6A near stop codon / 3'UTR boundary is more impactful
        if 0.75 < frac < 0.90:
            base_severity = 0.6
        elif frac > 0.85:
            base_severity = 0.4
        else:
            base_severity = 0.2

        # Boost severity if the site also matches the extended RRACH
        if mod_pos in rrach_positions:
            base_severity = min(1.0, base_severity * 1.3)
            motif_label = "m6A RRACH"
        else:
            motif_label = "m6A DRACH"

        sites.append(EpitranscriptomicSite(
            modification_type="m6A",
            position=mod_pos,
            sequence_context=match.group(),
            severity=base_severity,
            enzyme="METTL3/METTL14",
            reference="Dominissini et al. 2012 Nature; Meyer et al. 2012 Cell",
        ))

    # Add RRACH-only sites (those not already covered by DRACH)
    drach_positions: set[int] = set()
    for match in _M6A_DRACH.finditer(seq_upper):
        drach_positions.add(match.start() + 2)

    for match in _M6A_RRACH.finditer(seq_upper):
        mod_pos = match.start() + 3
        if mod_pos not in drach_positions:
            frac = match.start() / max(n, 1)
            if 0.75 < frac < 0.90:
                base_severity = 0.55
            elif frac > 0.85:
                base_severity = 0.35
            else:
                base_severity = 0.18

            sites.append(EpitranscriptomicSite(
                modification_type="m6A",
                position=mod_pos,
                sequence_context=match.group(),
                severity=base_severity,
                enzyme="METTL3/METTL14",
                reference="Dominissini et al. 2012 Nature; Meyer et al. 2012 Cell",
            ))

    return sites


# ────────────────────────────────────────────────────────────
# 4. m5C (5-Methylcytosine) Detection
# ────────────────────────────────────────────────────────────

# m5C consensus motifs for different writer enzymes
# (Squires et al. 2012 NAR; Yang et al. 2017 Mol Cell)
_M5C_MOTIFS: list[tuple[re.Pattern[str], str, str, float]] = [
    # (compiled regex, enzyme name, motif description, base_severity)
    # DNMT2: recognises tRNA-like CCGG motif in mRNA
    (re.compile(r"CCGG"), "DNMT2", "DNMT2-type m5C methylation (CCGG)", 0.45),
    # NSUN1 (DNMT1 homolog): CpG context in mRNA
    (re.compile(r"CG"), "NSUN1", "NSUN1-type m5C methylation (CpG)", 0.25),
    # NSUN2: C-rich context with flanking purines
    (re.compile(r"[AT]C[AT]C[AT]"), "NSUN2", "NSUN2-type m5C methylation (NCNCN)", 0.35),
    # NSUN3: mitochondrial, CC motif near start
    (re.compile(r"CC[AT]C"), "NSUN3", "NSUN3-type m5C methylation (CCNC)", 0.30),
    # NSUN4: C in GC-rich context
    (re.compile(r"GC[GC]C"), "NSUN4", "NSUN4-type m5C methylation (GCSC)", 0.30),
    # NSUN6: CCUCC motif
    (re.compile(r"CC[ATGC]CC"), "NSUN6", "NSUN6-type m5C methylation (CCNCC)", 0.28),
    # NSUN7: GCUG context (DNA: GCTG)
    (re.compile(r"GCTG"), "NSUN7", "NSUN7-type m5C methylation (GCTG)", 0.25),
]


def detect_m5c_sites(seq: str) -> list[EpitranscriptomicSite]:
    """Detect m5C (5-methylcytosine) consensus motifs in mRNA.

    Scans for m5C modification sites using enzyme-specific consensus
    motifs for DNMT2 and NSUN1-7 writer enzymes.  m5C modifications
    promote mRNA stability and nuclear export via ALYREF binding
    (Squires et al. 2012; Yang et al. 2017).

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of :class:`EpitranscriptomicSite` objects for predicted
        m5C sites.
    """
    sites: list[EpitranscriptomicSite] = []
    seq_upper = seq.upper()
    n = len(seq_upper)

    seen: set[tuple[str, int]] = set()  # (enzyme, position) dedup

    for pattern, enzyme, description, base_severity in _M5C_MOTIFS:
        for match in pattern.finditer(seq_upper):
            # The modified C is typically the first C in the motif
            mod_pos = match.start() + match.group().index("C")
            key = (enzyme, mod_pos)
            if key in seen:
                continue
            seen.add(key)

            frac = mod_pos / max(n, 1)

            # m5C in CDS is more impactful than in UTRs
            if 0.05 < frac < 0.90:
                severity = base_severity
            else:
                severity = base_severity * 0.7  # UTR positions less severe

            sites.append(EpitranscriptomicSite(
                modification_type="m5C",
                position=mod_pos,
                sequence_context=match.group(),
                severity=severity,
                enzyme=enzyme,
                reference="Squires et al. 2012 NAR",
            ))

    return sites


# ────────────────────────────────────────────────────────────
# 5. Pseudouridine (Ψ) Detection
# ────────────────────────────────────────────────────────────

# Pseudouridylation motifs for different writer enzymes
# (Carlile et al. 2014 Nature; Schwartz et al. 2014 Cell)
_PSI_MOTIFS: list[tuple[re.Pattern[str], str, str, float]] = [
    # Pus1: position-specific isomerisation, GU-rich context
    # RNA: GURUC → DNA: G[AT]RTC  (R=A/G)
    (re.compile(r"G[AT][AG]TC"), "Pus1",
     "Pus1-mediated pseudouridylation (GURUC)", 0.45),
    # Pus1 alternative: AUUC context
    (re.compile(r"ATTC"), "Pus1",
     "Pus1-mediated pseudouridylation (AUUC)", 0.35),
    # Pus4: highly conserved U2 snRNA-like motif
    # RNA: GUNUC → DNA: G[AT]TC
    (re.compile(r"GA[AT]C"), "Pus4",
     "Pus4-mediated pseudouridylation (GUNUC)", 0.40),
    # Pus7: UGUAR context (R=A/G)
    # RNA: UGUAR → DNA: TGT[AG]
    (re.compile(r"TGT[AG]"), "Pus7",
     "Pus7-mediated pseudouridylation (UGUAR)", 0.38),
    # Dyskerin/Cbf5 (H/ACA snoRNA-guided):
    # H/ACA snoRNAs form guide-target duplexes.  The target U is
    # typically in a single-stranded region flanked by two short
    # helices.  A simplified consensus is NNNANN'N'N' where N'
    # pairs with the guide.  We use an AT-rich bulge motif.
    (re.compile(r"[AT]{2}[AG]C[AT]{2}"), "Dyskerin/Cbf5",
     "H/ACA snoRNA-guided pseudouridylation (AT-rich bulge)", 0.35),
    # rpS13-guided Ψ in 18S rRNA-like context: GUAA
    (re.compile(r"GTAA"), "Dyskerin/Cbf5",
     "H/ACA snoRNA-guided pseudouridylation (GUAA)", 0.30),
]


def detect_pseudouridine_sites(seq: str) -> list[EpitranscriptomicSite]:
    """Detect pseudouridine (Ψ) modification sites in mRNA.

    Scans for pseudouridylation sites using enzyme-specific consensus
    motifs for Pus1, Pus4, Pus7, and Dyskerin/Cbf5 (H/ACA snoRNA-
    guided).  Pseudouridine enhances mRNA stability and alters codon-
    anticodon pairing properties (Carlile et al. 2014; Schwartz et al.
    2014).

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of :class:`EpitranscriptomicSite` objects for predicted
        Ψ sites.
    """
    sites: list[EpitranscriptomicSite] = []
    seq_upper = seq.upper()
    n = len(seq_upper)

    seen: set[tuple[str, int]] = set()

    for pattern, enzyme, description, base_severity in _PSI_MOTIFS:
        for match in pattern.finditer(seq_upper):
            # The pseudouridylated U (T in DNA) position
            # Typically the first T in the motif
            t_pos = match.group().find("T")
            if t_pos == -1:
                mod_pos = match.start()
            else:
                mod_pos = match.start() + t_pos

            key = (enzyme, mod_pos)
            if key in seen:
                continue
            seen.add(key)

            frac = mod_pos / max(n, 1)

            # Ψ in CDS can cause recoding; in UTRs mostly structural
            if 0.05 < frac < 0.90:
                severity = base_severity
            else:
                severity = base_severity * 0.8

            sites.append(EpitranscriptomicSite(
                modification_type="Psi",
                position=mod_pos,
                sequence_context=match.group(),
                severity=severity,
                enzyme=enzyme,
                reference="Carlile et al. 2014 Nature",
            ))

    return sites


# ────────────────────────────────────────────────────────────
# 6. m1A (N1-Methyladenosine) Detection
# ────────────────────────────────────────────────────────────

# m1A consensus motifs (Dominissini et al. 2016 Nature;
# Li et al. 2017 Mol Cell)
_M1A_MOTIFS: list[tuple[re.Pattern[str], str, str, float]] = [
    # TRMT6/61A consensus: GU-rich motif surrounding the target A
    # RNA: G/A-G-m1A-G/A-G-C → DNA: [AG][AG]A[AG][AG]C
    (re.compile(r"[AG][AG]A[AG][AG]C"), "TRMT6/61A",
     "TRMT6/61A m1A consensus (G/AGAG/AC)", 0.50),
    # tRNA T-loop consensus: TΨC loop → DNA: T[AT]C
    # The m1A at position 58 in tRNA is in the TΨC loop context
    (re.compile(r"T[AT]C"), "TRMT61B",
     "tRNA T-loop m1A (TΨC context)", 0.40),
    # m1A near start codon (5'UTR/cap-proximal)
    # Dominissini 2016 found m1A enriched near TSS/GU-rich context
    (re.compile(r"GA[AG]C[AT]G"), "TRMT6/61A",
     "m1A in GU-rich mRNA context", 0.35),
    # GC-rich hairpin loop context (Safra et al. 2017)
    (re.compile(r"GCA[GC]G"), "TRMT6/61A",
     "m1A in GC-rich hairpin loop", 0.38),
]


def detect_m1a_sites(seq: str) -> list[EpitranscriptomicSite]:
    """Detect m1A (N1-methyladenosine) consensus motifs in mRNA.

    Scans for m1A modification sites using consensus motifs for the
    TRMT6/61A complex (cytoplasmic mRNA) and TRMT61B (mitochondrial).
    m1A promotes translation initiation and enhances ribosome
    interaction (Dominissini et al. 2016).

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of :class:`EpitranscriptomicSite` objects for predicted
        m1A sites.
    """
    sites: list[EpitranscriptomicSite] = []
    seq_upper = seq.upper()
    n = len(seq_upper)

    seen: set[tuple[str, int]] = set()

    for pattern, enzyme, description, base_severity in _M1A_MOTIFS:
        for match in pattern.finditer(seq_upper):
            # The modified A position in the motif
            a_pos = match.group().find("A")
            if a_pos == -1:
                mod_pos = match.start()
            else:
                mod_pos = match.start() + a_pos

            key = (enzyme, mod_pos)
            if key in seen:
                continue
            seen.add(key)

            frac = mod_pos / max(n, 1)

            # m1A near 5' end / start codon is more impactful for translation
            if frac < 0.10:
                severity = min(1.0, base_severity * 1.4)  # 5'UTR/cap-proximal
            elif 0.05 < frac < 0.30:
                severity = base_severity  # Near start codon
            else:
                severity = base_severity * 0.8  # CDS/3'UTR

            sites.append(EpitranscriptomicSite(
                modification_type="m1A",
                position=mod_pos,
                sequence_context=match.group(),
                severity=severity,
                enzyme=enzyme,
                reference="Dominissini et al. 2016 Nature",
            ))

    return sites


# ────────────────────────────────────────────────────────────
# 7. 2'-O-Methylation (Nm) Detection
# ────────────────────────────────────────────────────────────

# 2'-O-methylation motifs (Fibrillarin/snoRNA-guided;
# Dong et al. 2012 NAR; Karijolich et al. 2015 RNA)
_2OM_MOTIFS: list[tuple[re.Pattern[str], str, str, float]] = [
    # Fibrillarin/snoRNA-guided C/D box snoRNA targets:
    # The guide region pairs with substrate, the methylated residue
    # is exactly 5 nt upstream of the D' or D box.
    # In mRNA, common Am sites appear in GGACU-like context
    # DNA: GG[AT]C  (Am at position 3)
    (re.compile(r"GG[AT]C"), "Fibrillarin/snoRNA",
     "C/D box snoRNA-guided 2'-O-Me (GGAC/GGAU)", 0.35),
    # Um (2'-O-methyluridine) in GAUC context
    # DNA: GATC
    (re.compile(r"GATC"), "Fibrillarin/snoRNA",
     "C/D box snoRNA-guided 2'-O-Me (GAUC)", 0.32),
    # Cm in CGG context (common in rRNA-like motifs)
    (re.compile(r"CGG"), "Fibrillarin/snoRNA",
     "C/D box snoRNA-guided 2'-O-Me (CGG)", 0.28),
    # Gm in AGC context
    (re.compile(r"AGC"), "Fibrillarin/snoRNA",
     "C/D box snoRNA-guided 2'-O-Me (AGC)", 0.25),
    # Cap 2'-O-Me at position 1 (first transcribed nucleotide)
    # detected separately in m6Am detection
]


def detect_2om_sites(seq: str) -> list[EpitranscriptomicSite]:
    """Detect 2'-O-methylation (Nm) sites in mRNA.

    Scans for 2'-O-methylation sites guided by C/D box snoRNAs
    via Fibrillarin.  2'-O-methylation enhances mRNA stability and
    enables innate immune evasion by preventing RIG-I and MDA5 sensing
    (Karijolich et al. 2015; Dong et al. 2012).

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of :class:`EpitranscriptomicSite` objects for predicted
        2'-O-Me sites.
    """
    sites: list[EpitranscriptomicSite] = []
    seq_upper = seq.upper()
    n = len(seq_upper)

    seen: set[tuple[str, int]] = set()

    for pattern, enzyme, description, base_severity in _2OM_MOTIFS:
        for match in pattern.finditer(seq_upper):
            # The 2'-O-methylated nucleotide is at position 2
            # (0-based within match)
            mod_pos = match.start() + 2

            key = (enzyme, mod_pos)
            if key in seen:
                continue
            seen.add(key)

            frac = mod_pos / max(n, 1)

            # Cap-proximal 2'-O-Me is critical for immune evasion
            if frac < 0.05:
                severity = min(1.0, base_severity * 1.5)
            elif 0.05 < frac < 0.90:
                severity = base_severity
            else:
                severity = base_severity * 0.8

            sites.append(EpitranscriptomicSite(
                modification_type="2OMe",
                position=mod_pos,
                sequence_context=match.group(),
                severity=severity,
                enzyme=enzyme,
                reference="Karijolich et al. 2015 RNA; Dong et al. 2012 NAR",
            ))

    return sites


# ────────────────────────────────────────────────────────────
# 8. m6Am (N6,2'-O-Dimethyladenosine) Detection
# ────────────────────────────────────────────────────────────

# m6Am occurs at the first transcribed nucleotide (cap-proximal A)
# when it follows the 7-methylguanosine cap (m7Gppp).
# Writer: PCIF1/CAPAM; Eraser: FTO
# The BCAA motif (B=C/G/U, i.e. not A) at the +1 position indicates
# potential m6Am: the consensus is [CGT]A at positions +1/+2.
# In practice, m6Am requires the cap context (m7GpppN), so we detect
# the BCAA motif at the 5' end of the sequence.
# BCAA motif: B = C/G/T (not A), then C, A, A
# DNA: [CGT]CAA at position 0-3

_M6AM_BCAA = re.compile(r"^[CGT]CAA")


def detect_m6am_sites(seq: str) -> list[EpitranscriptomicSite]:
    """Detect m6Am (N6,2'-O-dimethyladenosine) sites at transcription start.

    m6Am occurs at the first transcribed nucleotide when it is an
    adenosine, following the m7Gppp cap.  The BCAA motif (B=C/G/U,
    C, A, A) at the transcription start site is a strong indicator
    of m6Am modification.  m6Am promotes mRNA stability and is a
    target of the FTO demethylase (Mauer et al. 2017).

    For synthetic mRNA, the cap is added enzymatically, so this
    detection identifies cap-proximal BCAA motifs that would be
    subject to m6Am modification in vivo.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U), starting from
            the first transcribed nucleotide.

    Returns:
        List of :class:`EpitranscriptomicSite` objects for predicted
        m6Am sites (at most one per sequence, at position 1).
    """
    sites: list[EpitranscriptomicSite] = []
    seq_upper = seq.upper()

    match = _M6AM_BCAA.match(seq_upper)
    if match:
        # m6Am is on the A at position 1 (second nucleotide, 0-based)
        sites.append(EpitranscriptomicSite(
            modification_type="m6Am",
            position=1,
            sequence_context=match.group(),
            severity=0.65,  # m6Am has strong stabilising effect
            enzyme="PCIF1/CAPAM",
            reference="Mauer et al. 2017 Nature; Wei et al. 1976 Cell",
        ))

    # Also check for m6Am with just an A at +1 without full BCAA
    # (weaker prediction but still relevant)
    if len(seq_upper) >= 2 and seq_upper[1] == "A" and not match:
        # A at +1 position without full BCAA motif — lower confidence
        sites.append(EpitranscriptomicSite(
            modification_type="m6Am",
            position=1,
            sequence_context=seq_upper[:4] if len(seq_upper) >= 4 else seq_upper,
            severity=0.30,  # weaker prediction
            enzyme="PCIF1/CAPAM",
            reference="Mauer et al. 2017 Nature",
        ))

    return sites


# ────────────────────────────────────────────────────────────
# 9. Composite Detection Function
# ────────────────────────────────────────────────────────────

# Map from short mark name to detector function
_DETECTOR_MAP: dict[str, tuple[Any, str]] = {
    "m6a": (detect_m6a_sites, "m6A (N6-methyladenosine)"),
    "m5c": (detect_m5c_sites, "m5C (5-methylcytosine)"),
    "pseudouridine": (detect_pseudouridine_sites, "Pseudouridine (Ψ)"),
    "psi": (detect_pseudouridine_sites, "Pseudouridine (Ψ)"),
    "m1a": (detect_m1a_sites, "m1A (N1-methyladenosine)"),
    "2om": (detect_2om_sites, "2'-O-methylation"),
    "2ome": (detect_2om_sites, "2'-O-methylation"),
    "2'o-me": (detect_2om_sites, "2'-O-methylation"),
    "m6am": (detect_m6am_sites, "m6Am (N6,2'-O-dimethyladenosine)"),
}

_ALL_MARKS: list[str] = ["m6a", "m5c", "pseudouridine", "m1a", "2om", "m6am"]


def detect_all_epitranscriptomic_marks(
    seq: str,
    marks: list[str] | None = None,
) -> list[EpitranscriptomicSite]:
    """Detect all requested epitranscriptomic modification marks in mRNA.

    Runs the specified modification detectors and returns combined
    results sorted by severity (highest first).  When *marks* is
    ``None``, all six modification types are detected.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        marks: Optional list of mark names to detect.  Accepted names
            (case-insensitive): ``"m6a"``, ``"m5c"``,
            ``"pseudouridine"`` (or ``"psi"``), ``"m1a"``, ``"2om"``
            (or ``"2ome"``, ``"2'o-me"``), ``"m6am"``.  Defaults to
            ``None`` (detect all).

    Returns:
        Combined list of :class:`EpitranscriptomicSite` objects from
        all requested detectors, sorted by severity (highest first).

    Raises:
        ValueError: If any mark name in *marks* is not recognised.
    """
    if marks is None:
        marks = _ALL_MARKS

    all_sites: list[EpitranscriptomicSite] = []

    for mark in marks:
        mark_lower = mark.lower()
        entry = _DETECTOR_MAP.get(mark_lower)
        if entry is None:
            valid = sorted(_DETECTOR_MAP.keys())
            raise ValueError(
                f"Unknown epitranscriptomic mark: {mark!r}. "
                f"Valid marks: {valid}"
            )
        detector_fn, description = entry
        try:
            sites = detector_fn(seq)
            all_sites.extend(sites)
        except Exception:
            logger.warning(
                "Detector for %s failed, skipping", description,
                exc_info=True,
            )

    # Sort by severity (highest first), then by position
    all_sites.sort(key=lambda s: (-s.severity, s.position))

    return all_sites


# ────────────────────────────────────────────────────────────
# 10. Impact Assessment
# ────────────────────────────────────────────────────────────


def assess_stability_impact(
    sites: list[EpitranscriptomicSite],
) -> dict[str, Any]:
    """Assess how detected epitranscriptomic marks affect mRNA stability.

    Aggregates stability, translation, and immune-evasion contributions
    from each detected modification site using the
    :data:`MODIFICATION_FUNCTIONS` database.

    Args:
        sites: List of detected :class:`EpitranscriptomicSite` objects.

    Returns:
        Dictionary with the following keys:

        - ``"total_sites"``: Total number of detected sites.
        - ``"sites_by_type"``: Count of sites per modification type.
        - ``"net_stability_impact"``: Weighted sum of stability
          contributions (positive = stabilising, negative =
          destabilising).
        - ``"net_translation_impact"``: Weighted sum of translation
          contributions.
        - ``"stabilising_marks"``: List of modification types that
          contribute positively to stability.
        - ``"destabilising_marks"``: List of modification types that
          contribute negatively to stability.
        - ``"immune_evasion_score"``: Score in [0, 1] representing
          innate immune evasion potential (primarily from 2'-O-Me).
        - ``"high_severity_sites"``: Sites with severity >= 0.5.
        - ``"recommendation"``: Human-readable recommendation string.
    """
    if not sites:
        return {
            "total_sites": 0,
            "sites_by_type": {},
            "net_stability_impact": 0.0,
            "net_translation_impact": 0.0,
            "stabilising_marks": [],
            "destabilising_marks": [],
            "immune_evasion_score": 0.0,
            "high_severity_sites": [],
            "recommendation": "No epitranscriptomic sites detected.",
        }

    sites_by_type: dict[str, int] = {}
    net_stability = 0.0
    net_translation = 0.0
    immune_evasion = 0.0
    stabilising: list[str] = []
    destabilising: list[str] = []
    high_severity: list[EpitranscriptomicSite] = []

    for site in sites:
        mod_type = site.modification_type
        sites_by_type[mod_type] = sites_by_type.get(mod_type, 0) + 1

        if site.severity >= 0.5:
            high_severity.append(site)

        func_info = MODIFICATION_FUNCTIONS.get(mod_type)
        if func_info is not None:
            # Weight the impact by severity
            weight = site.severity
            net_stability += func_info["stability_impact"] * weight
            net_translation += func_info["translation_impact"] * weight

            if func_info["stability_impact"] > 0 and mod_type not in stabilising:
                stabilising.append(mod_type)
            elif func_info["stability_impact"] < 0 and mod_type not in destabilising:
                destabilising.append(mod_type)

        # 2'-O-Me contributes to immune evasion
        if mod_type == "2OMe":
            immune_evasion += site.severity * 0.5

    # Cap immune evasion at 1.0
    immune_evasion = min(1.0, immune_evasion)

    # Build recommendation
    if net_stability > 0.5:
        rec = (
            "Net epitranscriptomic landscape is strongly stabilising. "
            "Consider retaining m5C and Ψ sites for therapeutic mRNA design."
        )
    elif net_stability > 0:
        rec = (
            "Net epitranscriptomic landscape is mildly stabilising. "
            "Most detected marks support mRNA stability."
        )
    elif net_stability > -0.3:
        rec = (
            "Net epitranscriptomic landscape is mildly destabilising. "
            "Consider reducing m6A sites near stop codons/3'UTR if "
            "stability is a priority."
        )
    else:
        rec = (
            "Net epitranscriptomic landscape is strongly destabilising. "
            "High m6A burden may accelerate decay via YTHDF2. Consider "
            "codon optimisation to reduce DRACH/RRACH motifs near 3' end."
        )

    if immune_evasion > 0.3:
        rec += " 2'-O-Me sites present may aid immune evasion for in vivo delivery."

    return {
        "total_sites": len(sites),
        "sites_by_type": sites_by_type,
        "net_stability_impact": round(net_stability, 4),
        "net_translation_impact": round(net_translation, 4),
        "stabilising_marks": stabilising,
        "destabilising_marks": destabilising,
        "immune_evasion_score": round(immune_evasion, 4),
        "high_severity_sites": high_severity,
        "recommendation": rec,
    }
