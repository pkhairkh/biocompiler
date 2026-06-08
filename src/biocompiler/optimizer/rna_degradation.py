"""BioCompiler RNA Degradation Analysis Module
==============================================

Comprehensive mRNA degradation signal detection and scoring for
therapeutic mRNA and heterologous expression design.

This module models the major eukaryotic and prokaryotic mRNA decay
pathways, including:

1. **miRNA-mediated silencing** — Seed-match detection against a curated
   miRNA database with ViennaRNA accessibility-aware severity scoring.
2. **AU-rich elements (AREs)** — Class I/II/III ARE detection with
   position-dependent weighting.
3. **m6A epitranscriptomic marks** — DRACH consensus motif detection
   affecting mRNA stability.
4. **Eukaryotic decay pathways** — XRN1 5'→3' decay, CCR4-NOT
   deadenylation, exosome 3'→5' decay, and DCP2 decapping modeling
   (replacing the bacterial RNase E-centric approach).
5. **Nonsense-mediated decay (NMD)** — Detection of premature
   termination codons (EJC-dependent and long 3'UTR pathways).
6. **Epitranscriptomic modifications** — m5C, pseudouridine (Ψ),
   m1A, and 2'-O-methylation motif detection.
7. **ViennaRNA accessibility** — Per-position mRNA accessibility
   computation using base-pairing probabilities (with GC-content
   heuristic fallback), integrated into miRNA binding site severity
   scoring.

Usage::

    from biocompiler.optimizer.rna_degradation import (
        DegradationSignal,
        detect_mirna_sites,
        detect_are_signals,
        detect_m6a_sites,
        detect_eukaryotic_decay_signals,
        detect_nmd_triggers,
        detect_epitranscriptomic_marks,
        compute_mrna_accessibility,
        analyze_rna_degradation,
    )

    # Full degradation analysis
    report = analyze_rna_degradation("ATGGCC...TAA", organism="Homo_sapiens")
    for signal in report:
        print(f"{signal.signal_type}: {signal.description} (severity={signal.severity:.2f})")

References:
  Bartel, D.P. (2009). "MicroRNAs: target recognition and regulatory
  functions." *Cell* 136:215–233.
  Shaw, G. & Kamen, R. (1986). "A conserved AU sequence from the 3'
  untranslated region of GM-CSF mRNA mediates selective mRNA
  degradation." *Cell* 46:659–667.
  Dominissini, D. et al. (2012). "Topology of the human and mouse m6A
  RNA methylomes revealed by m6A-seq." *Nature* 485:201–206.
  Jones, C.I. et al. (2012). "RNA-seq reveals transcriptome changes
  in nonsense-mediated mRNA decay-depleted cells." *Nature* (XRN1
  decay pathway).
  Yamashita, A. et al. (2005). "Concerted action of the RNA helicase
  and the CCR4-NOT complex in deadenylation." *Cell* (CCR4-NOT).
  Mitchell, P. et al. (1997). "The exosome: a conserved eukaryotic
  RNA processing complex." *PNAS* (Exosome 3'→5' decay).
  van Dijk, E. et al. (2002). "Human Dcp2: a catalytically active
  mRNA decapping enzyme." *EMBO J* (DCP2 decapping).
  Squires, J.E. et al. (2012). "Widespread occurrence of 5-
  methylcytosine in human coding and non-coding RNA." *Nucleic Acids
  Research* (m5C).
  Carlile, T.M. et al. (2014). "Pseudouridine profiling reveals
  regulated mRNA pseudouridylation in yeast and human cells."
  *Nature* (pseudouridine).
  Dominissini, D. et al. (2016). "The dynamic N1-methyladenosine
  methylome in eukaryotic messenger RNA." *Nature* (m1A).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "DegradationSignal",
    "detect_mirna_sites",
    "detect_are_signals",
    "detect_m6a_sites",
    "detect_eukaryotic_decay_signals",
    "detect_nmd_triggers",
    "detect_epitranscriptomic_marks",
    "compute_mrna_accessibility",
    "analyze_rna_degradation",
    "_MIRNA_DATABASE",
    "_EUKARYOTIC_DECAY_PATHWAYS",
]

# ────────────────────────────────────────────────────────────
# 1. Core Data Structures
# ────────────────────────────────────────────────────────────


@dataclass
class DegradationSignal:
    """A single mRNA degradation signal detected in a sequence.

    Attributes:
        signal_type: Category of the degradation signal (e.g.,
            ``"mirna"``, ``"are_class_ii"``, ``"m6a"``,
            ``"ccr4_not_deadenylation"``, ``"nmd_ejc_dependent"``).
        position: 0-based position in the sequence where the signal
            starts.
        sequence_context: The matched subsequence or a descriptive
            context string.
        severity: Severity score in [0, 1].  Higher values indicate
            greater risk of accelerated degradation.
        description: Human-readable description of the signal.
        pathway: Decay pathway this signal contributes to (e.g.,
            ``"mirna"``, ``"are"``, ``"epitranscriptomic_m6a"``,
            ``"deadenylation"``, ``"nmd"``).
    """

    signal_type: str
    position: int
    sequence_context: str
    severity: float
    description: str
    pathway: str = ""

    def __post_init__(self) -> None:
        self.severity = max(0.0, min(1.0, self.severity))


# ────────────────────────────────────────────────────────────
# 2. miRNA Database
# ────────────────────────────────────────────────────────────

# Curated miRNA seed-match database for therapeutically relevant miRNAs.
# Each entry contains:
#   - name: miRNA identifier
#   - mature_seq: Mature miRNA sequence (RNA, U not T)
#   - seed_2_8: Positions 2-8 of the mature sequence (DNA, T not U),
#     which is the canonical seed region for target recognition
#   - tissue: Tissues where this miRNA is highly expressed
#   - disease_relevance: Clinical relevance for mRNA therapeutic design
#   - conservation: Evolutionary conservation level

_MIRNA_DATABASE: dict[str, dict[str, Any]] = {
    "miR-1": {
        "name": "miR-1",
        "mature_seq": "UGGAAUGUAAAGAAGUAUGUAU",
        "seed_2_8": "GGAATGT",
        "tissue": ["heart", "skeletal_muscle"],
        "disease_relevance": "Cardiomyopathy",
        "conservation": "high",
    },
    "miR-122": {
        "name": "miR-122",
        "mature_seq": "UGGAGUGUGACAAUGGUGUUUG",
        "seed_2_8": "GGAGTGT",
        "tissue": ["liver"],
        "disease_relevance": "Hepatitis C, hepatocellular carcinoma",
        "conservation": "high",
    },
    "miR-143": {
        "name": "miR-143",
        "mature_seq": "UGAGAUGAAGCACUGUAGCUC",
        "seed_2_8": "GAGATGA",
        "tissue": ["colon", "prostate", "adipose"],
        "disease_relevance": "Colorectal cancer",
        "conservation": "medium",
    },
    "miR-145": {
        "name": "miR-145",
        "mature_seq": "GUCCAGUUUUCCCAGGAAUCCC",
        "seed_2_8": "UCCAGUU",
        "tissue": ["colon", "breast", "vascular_smooth_muscle"],
        "disease_relevance": "Colorectal cancer, vascular disease",
        "conservation": "high",
    },
    "miR-155": {
        "name": "miR-155",
        "mature_seq": "UUAAUGCUAAUCGUGAUAGGGGU",
        "seed_2_8": "TAATGCT",
        "tissue": ["immune_cells", "lymph_node", "spleen"],
        "disease_relevance": "Lymphoma, inflammation",
        "conservation": "high",
    },
    "miR-200b": {
        "name": "miR-200b",
        "mature_seq": "UAAUACUGCCUGGUAAUGAUG",
        "seed_2_8": "AATACTG",
        "tissue": ["kidney", "breast", "ovary"],
        "disease_relevance": "Epithelial-mesenchymal transition, cancer metastasis",
        "conservation": "high",
    },
    "miR-200c": {
        "name": "miR-200c",
        "mature_seq": "UAAUACUGCCUGGUAAUGAUG",
        "seed_2_8": "AATACTG",
        "tissue": ["kidney", "breast", "ovary"],
        "disease_relevance": "EMT, cancer metastasis",
        "conservation": "high",
    },
    "miR-21": {
        "name": "miR-21",
        "mature_seq": "UAGCUUAUCAGACUGAUGUUGA",
        "seed_2_8": "AGCTTAT",
        "tissue": ["ubiquitous", "lung", "breast", "colon"],
        "disease_relevance": "Glioblastoma, multiple cancers",
        "conservation": "high",
    },
    "miR-223": {
        "name": "miR-223",
        "mature_seq": "UGUCAGUUUGUCAAAUACCCCA",
        "seed_2_8": "GTCAGTT",
        "tissue": ["bone_marrow", "neutrophils", "immune_cells"],
        "disease_relevance": "Inflammation, leukemia",
        "conservation": "medium",
    },
    "miR-29a": {
        "name": "miR-29a",
        "mature_seq": "UAGCACCAUCUGAAAUCGGUUA",
        "seed_2_8": "AGCACCA",
        "tissue": ["lung", "kidney", "spleen"],
        "disease_relevance": "Pulmonary fibrosis, diabetes",
        "conservation": "high",
    },
    "miR-34a": {
        "name": "miR-34a",
        "mature_seq": "UGGCAGUGUCUUAGCUGGUUGU",
        "seed_2_8": "GGCAGTG",
        "tissue": ["brain", "testis", "spleen"],
        "disease_relevance": "p53 pathway, tumor suppressor",
        "conservation": "high",
    },
    "let-7a": {
        "name": "let-7a",
        "mature_seq": "UGAGGUAGUAGGUUGUAUAGUU",
        "seed_2_8": "GAGGTAG",
        "tissue": ["ubiquitous", "lung", "brain"],
        "disease_relevance": "Lung cancer, tumor suppressor",
        "conservation": "high",
    },
    "miR-142-3p": {
        "name": "miR-142-3p",
        "mature_seq": "UGUAGUGUUUCCUACUUUAUGGA",
        "seed_2_8": "GTAGTGT",
        "tissue": ["hematopoietic", "immune_cells"],
        "disease_relevance": "Hematological malignancies",
        "conservation": "medium",
    },
    "miR-146a": {
        "name": "miR-146a",
        "mature_seq": "UGAGAACUGAAUUCCAUGGGUU",
        "seed_2_8": "GAGAACT",
        "tissue": ["immune_cells", "thymus", "spleen"],
        "disease_relevance": "Inflammation, innate immunity",
        "conservation": "high",
    },
    "miR-16": {
        "name": "miR-16",
        "mature_seq": "UAGCAGCACGUAAAUAUUGGCG",
        "seed_2_8": "AGCAGCA",
        "tissue": ["ubiquitous", "lymph_node", "spleen"],
        "disease_relevance": "Chronic lymphocytic leukemia",
        "conservation": "high",
    },
    "miR-125b": {
        "name": "miR-125b",
        "mature_seq": "UCCCUGAGACCCUAACUUGUGA",
        "seed_2_8": "CCCTGAG",
        "tissue": ["brain", "breast", "hematopoietic"],
        "disease_relevance": "Neurodegeneration, breast cancer",
        "conservation": "high",
    },
}


# ────────────────────────────────────────────────────────────
# 3. AU-Rich Element (ARE) Detection
# ────────────────────────────────────────────────────────────

# ARE classification follows Chen & Shyu (1995):
#   Class I: scattered AUUUA in U-rich context
#   Class II: overlapping AUUUA repeats (rapid decay)
#   Class III: U-rich without AUUUA

_ARE_PATTERNS: list[tuple[str, str, str, float]] = [
    # (pattern, are_class, description, base_severity)
    (r"ATTTA", "I", "ARE Class I (scattered AUUUA in U-rich context)", 0.3),
    (r"ATTTATTTA", "II", "ARE Class II (overlapping AUUUA repeats - rapid decay)", 0.7),
    (r"ATTTATTTATTTA", "II", "ARE Class II (triple AUUUA repeat - very rapid decay)", 0.9),
    (r"T{6,}", "III", "ARE Class III (U-rich element without AUUUA)", 0.2),
]


def detect_are_signals(seq: str) -> list[DegradationSignal]:
    """Detect AU-rich element (ARE) degradation signals in an mRNA sequence.

    Scans for Class I, II, and III ARE motifs following the Chen & Shyu
    (1995) classification scheme.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of DegradationSignal objects for detected ARE motifs.
    """
    signals: list[DegradationSignal] = []
    seq_upper = seq.upper()

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
            signals.append(DegradationSignal(
                signal_type=f"are_class_{are_class.lower()}",
                position=match.start(),
                sequence_context=match.group(),
                severity=severity,
                description=f"{description} at position {match.start()}: {match.group()}",
                pathway="are",
            ))

    return signals


# ────────────────────────────────────────────────────────────
# 4. m6A (N6-Methyladenosine) Detection
# ────────────────────────────────────────────────────────────

# m6A consensus motif: DRACH (D=A/G/U, R=A/G, H=A/C/U)
# In DNA alphabet: [AG][AG]AC[ACT]
_M6A_CONSENSUS = r"[AG][AG]AC[ACT]"


def detect_m6a_sites(seq: str) -> list[DegradationSignal]:
    """Detect m6A (N6-methyladenosine) consensus motifs in mRNA.

    m6A modifications at DRACH motifs can affect mRNA stability,
    splicing, and translation.  The presence of m6A sites near stop
    codons and in 3'UTRs is particularly impactful for stability
    (Dominissini et al. 2012; Wang et al. 2014).

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of DegradationSignal objects for predicted m6A sites.
    """
    signals: list[DegradationSignal] = []
    seq_upper = seq.upper()

    for match in re.finditer(_M6A_CONSENSUS, seq_upper):
        frac = match.start() / max(len(seq_upper), 1)

        # m6A near stop codon / 3'UTR boundary is more impactful
        if 0.75 < frac < 0.90:
            severity = 0.6  # Near stop codon — strong effect
        elif frac > 0.85:
            severity = 0.4  # 3'UTR — moderate effect
        else:
            severity = 0.2  # CDS — weaker effect

        signals.append(DegradationSignal(
            signal_type="m6a_drach",
            position=match.start(),
            sequence_context=match.group(),
            severity=severity,
            description=f"m6A DRACH motif at position {match.start()}: {match.group()}",
            pathway="epitranscriptomic_m6a",
        ))

    return signals


# ────────────────────────────────────────────────────────────
# 5. Eukaryotic Decay Machinery (Replacing Bacterial RNase E)
# ────────────────────────────────────────────────────────────

# Eukaryotic mRNA decay pathways. RNase E is a bacterial endonuclease
# and not relevant for eukaryotic mRNA therapeutic design.  These
# pathways model the major eukaryotic decay mechanisms instead.

_EUKARYOTIC_DECAY_PATHWAYS = {
    "xrn1_5to3": {
        "name": "XRN1 5'→3' Decay",
        "description": "Major 5'→3' exoribonuclease, acts on decapped mRNA",
        "signal": "decapped_mRNA",
        "rate_modifier": 2.0,
        "tissue_specific": False,
        "reference": "Jones et al. 2012 Nature",
    },
    "ccr4_not_deadenylation": {
        "name": "CCR4-NOT Deadenylation",
        "description": "PolyA tail shortening, rate depends on ARE class",
        "rate_by_are_class": {
            "I": 1.5,   # Class I ARE: moderate deadenylation acceleration
            "II": 3.0,   # Class II ARE: strong deadenylation acceleration
            "III": 1.0,  # Class III ARE: no acceleration
        },
        "tissue_specific": False,
        "reference": "Yamashita et al. 2005 Cell",
    },
    "exosome_3to5": {
        "name": "Exosome 3'→5' Decay",
        "description": "3'→5' exoribonuclease complex, backup pathway",
        "rate_modifier": 0.5,  # Slower than XRN1
        "tissue_specific": False,
        "reference": "Mitchell et al. 1997 PNAS",
    },
    "decapping_dcp2": {
        "name": "DCP2 Decapping",
        "description": "Removes 5' cap, prerequisite for XRN1",
        "signal": "short_polyA",  # Triggered when polyA tail is short
        "rate_modifier": 1.0,
        "tissue_specific": False,
        "reference": "van Dijk et al. 2002 EMBO J",
    },
}


def detect_eukaryotic_decay_signals(seq: str, are_class: str = "none",
                                     polyA_tail_length: int = 200) -> list[DegradationSignal]:
    """Detect eukaryotic mRNA decay pathway activation signals.

    Models the major eukaryotic decay pathways:
    1. CCR4-NOT deadenylation (rate depends on ARE class)
    2. DCP2 decapping (triggered by short polyA tail)
    3. XRN1 5'→3' decay (requires decapping)
    4. Exosome 3'→5' decay (backup pathway)

    Args:
        seq: mRNA sequence (DNA alphabet, T not U)
        are_class: ARE class detected in sequence ("I", "II", "III", "none")
        polyA_tail_length: Length of polyA tail in nucleotides

    Returns:
        List of DegradationSignal objects for activated decay pathways
    """
    signals: list[DegradationSignal] = []
    seq_upper = seq.upper()

    # CCR4-NOT deadenylation rate
    ccr4_rate = _EUKARYOTIC_DECAY_PATHWAYS["ccr4_not_deadenylation"]["rate_by_are_class"]
    deadenylation_rate = ccr4_rate.get(are_class, 1.0)

    if deadenylation_rate > 1.0:
        signals.append(DegradationSignal(
            signal_type="ccr4_not_deadenylation",
            position=0,
            sequence_context=f"ARE_class_{are_class}",
            severity=min(1.0, deadenylation_rate / 3.0),
            description=f"CCR4-NOT deadenylation accelerated {deadenylation_rate:.1f}x by ARE Class {are_class}",
            pathway="deadenylation",
        ))

    # DCP2 decapping - triggered when polyA is short (<50 nt)
    if polyA_tail_length < 50:
        cap_severity = max(0, 1.0 - polyA_tail_length / 50.0)
        signals.append(DegradationSignal(
            signal_type="decapping_dcp2",
            position=0,
            sequence_context=f"polyA_{polyA_tail_length}nt",
            severity=cap_severity,
            description=f"DCP2 decapping triggered by short polyA tail ({polyA_tail_length} nt)",
            pathway="decapping",
        ))

    # XRN1 5'→3' decay - always active on decapped mRNA
    if polyA_tail_length < 80:  # Partially deadenylated
        xrn1_severity = min(1.0, deadenylation_rate * 0.5)
        signals.append(DegradationSignal(
            signal_type="xrn1_5to3",
            position=0,
            sequence_context="decapped",
            severity=xrn1_severity,
            description="XRN1 5'→3' decay active on decapped mRNA",
            pathway="5_to_3_decay",
        ))

    # Exosome 3'→5' decay - backup, always present at low level
    signals.append(DegradationSignal(
        signal_type="exosome_3to5",
        position=len(seq) - 1,
        sequence_context="3prime_end",
        severity=0.2,  # Low baseline
        description="Exosome 3'→5' decay (backup pathway, low rate)",
        pathway="3_to_5_decay",
    ))

    return signals


# ────────────────────────────────────────────────────────────
# 6. Nonsense-Mediated Decay (NMD) Detection
# ────────────────────────────────────────────────────────────


def detect_nmd_triggers(seq: str, orf_start: int = 0, orf_end: int | None = None,
                         has_introns: bool = False,
                         last_exon_junction: int | None = None) -> list[DegradationSignal]:
    """Detect nonsense-mediated decay triggers in mRNA.

    For synthetic mRNA (no introns):
    - Long 3'UTR (>2000 nt) triggers NMD via HNRNPDL/mTOR pathway

    For genes with introns:
    - PTC ≥50-55 nt upstream of last exon-exon junction triggers EJC-dependent NMD

    Args:
        seq: mRNA sequence (DNA alphabet)
        orf_start: Start position of main ORF (0-indexed)
        orf_end: End position of main ORF (0-indexed, exclusive)
        has_introns: Whether the gene has introns (EJC-dependent NMD)
        last_exon_junction: Position of last exon-exon junction (0-indexed)

    Returns:
        List of DegradationSignal objects for NMD triggers
    """
    signals: list[DegradationSignal] = []
    seq_upper = seq.upper()

    if orf_end is None:
        # Find the ORF end by scanning for stop codons
        for i in range(orf_start, len(seq_upper) - 2, 3):
            codon = seq_upper[i:i+3]
            if codon in ("TAA", "TAG", "TGA"):
                orf_end = i + 3
                break
        if orf_end is None:
            return signals

    # Check for premature stop codons (PTCs) in the ORF
    for i in range(orf_start + 3, orf_end - 2, 3):  # Skip initial ATG
        codon = seq_upper[i:i+3]
        if codon in ("TAA", "TAG", "TGA") and i < orf_end - 3:
            # This is a premature stop codon
            distance_to_end = orf_end - i

            if has_introns and last_exon_junction is not None:
                # EJC-dependent NMD: PTC ≥50 nt upstream of last EJ
                distance_to_ej = last_exon_junction - i
                if distance_to_ej >= 50:
                    severity = min(1.0, 0.3 + 0.7 * (distance_to_ej / 500))
                    signals.append(DegradationSignal(
                        signal_type="nmd_ejc_dependent",
                        position=i,
                        sequence_context=seq_upper[max(0,i-3):i+6],
                        severity=severity,
                        description=f"EJC-dependent NMD: PTC at {i}, {distance_to_ej}nt upstream of last EJ",
                        pathway="nmd",
                    ))
            else:
                # Intronless NMD: check 3'UTR length
                utr3_length = len(seq_upper) - orf_end
                if utr3_length > 2000:
                    severity = min(1.0, utr3_length / 5000)
                    signals.append(DegradationSignal(
                        signal_type="nmd_long_3utr",
                        position=orf_end,
                        sequence_context=f"3UTR_{utr3_length}nt",
                        severity=severity,
                        description=f"Long 3'UTR NMD: {utr3_length}nt 3'UTR (threshold: 2000nt)",
                        pathway="nmd",
                    ))

    return signals


# ────────────────────────────────────────────────────────────
# 7. Epitranscriptomic Modifications (beyond m6A)
# ────────────────────────────────────────────────────────────

# m5C (5-methylcytosine) consensus motifs (Squires et al. 2012)
_M5C_CONSENSUS = [
    (r"[CT]C[ACT]G", "m5C_DNMT2", 0.4, "DNMT2-type m5C methylation motif"),
    (r"CC[AT]C", "m5C_NSCUN1", 0.3, "NSUN1-type m5C methylation motif"),
    (r"[AT]C[AT]C[AT]", "m5C_NSCUN2", 0.35, "NSUN2-type m5C methylation motif"),
]

# Pseudouridine (Ψ) hotspots (Carlile et al. 2014)
_PSEUDOURIDINE_MOTIFS = [
    (r"GATCT", "Psi_GATCT", 0.5, "Pus1-mediated pseudouridylation motif"),
    (r"GACTGA", "Psi_GACTGA", 0.4, "Pus4-mediated pseudouridylation motif"),
    (r"RU[AC]U[AG]R", "Psi_dyskerin", 0.35, "Dyskerin/Cbf5 H/ACA snoRNA-guided motif"),
]

# m1A (N1-methyladenosine) motifs (Dominissini et al. 2016)
_M1A_MOTIFS = [
    (r"G[AG]G[AG]C", "m1A_tRNA_Tloop", 0.4, "m1A in tRNA T-loop consensus"),
    (r"[AT]C[AT]G", "m1A_mRNA", 0.3, "m1A in mRNA GU-rich motif"),
    (r"TCC", "m1A_TCC", 0.25, "m1A near TCC trinucleotide"),
]

# 2'-O-methylation (Nm) motifs
_2OM_MOTIFS = [
    (r"GG[AT]C", "Am_fibrillarin", 0.3, "Fibrillarin-guided 2'-O-methylation"),
]


def detect_epitranscriptomic_marks(seq: str) -> list[DegradationSignal]:
    """Detect epitranscriptomic modification motifs beyond m6A.

    Covers: m5C, pseudouridine (Ψ), m1A, and 2'-O-methylation.

    Args:
        seq: mRNA sequence (DNA alphabet)

    Returns:
        List of DegradationSignal objects for epitranscriptomic sites
    """
    signals: list[DegradationSignal] = []
    seq_upper = seq.upper()

    # m5C detection
    for pattern, name, severity, desc in _M5C_CONSENSUS:
        for match in re.finditer(pattern, seq_upper):
            signals.append(DegradationSignal(
                signal_type=f"m5C_{name}",
                position=match.start(),
                sequence_context=match.group(),
                severity=severity,
                description=f"m5C site: {desc}",
                pathway="epitranscriptomic_m5c",
            ))

    # Pseudouridine detection
    for pattern, name, severity, desc in _PSEUDOURIDINE_MOTIFS:
        expanded = pattern.replace("R", "[AG]").replace("Y", "[CT]").replace("U", "T")
        for match in re.finditer(expanded, seq_upper):
            signals.append(DegradationSignal(
                signal_type=f"psi_{name}",
                position=match.start(),
                sequence_context=match.group(),
                severity=severity,
                description=f"Ψ site: {desc}",
                pathway="epitranscriptomic_psi",
            ))

    # m1A detection
    for pattern, name, severity, desc in _M1A_MOTIFS:
        for match in re.finditer(pattern, seq_upper):
            signals.append(DegradationSignal(
                signal_type=f"m1A_{name}",
                position=match.start(),
                sequence_context=match.group(),
                severity=severity,
                description=f"m1A site: {desc}",
                pathway="epitranscriptomic_m1a",
            ))

    # 2'-O-methylation detection
    for pattern, name, severity, desc in _2OM_MOTIFS:
        for match in re.finditer(pattern, seq_upper):
            signals.append(DegradationSignal(
                signal_type=f"2OM_{name}",
                position=match.start(),
                sequence_context=match.group(),
                severity=severity,
                description=f"2'-O-Me site: {desc}",
                pathway="epitranscriptomic_2om",
            ))

    return signals


# ────────────────────────────────────────────────────────────
# 8. ViennaRNA Accessibility Integration
# ────────────────────────────────────────────────────────────


def compute_mrna_accessibility(seq: str, window: int = 80, span: int = 40,
                                unpaired_region: int = 20) -> list[float]:
    """Compute per-position mRNA accessibility (unpaired probability).

    Uses ViennaRNA RNAplfold-equivalent computation when available,
    falls back to GC-content heuristic.

    Args:
        seq: mRNA sequence (DNA alphabet)
        window: Sliding window size (default 80)
        span: Maximum base pair span (default 40)
        unpaired_region: Region length for unpaired probability (default 20)

    Returns:
        List of per-position accessibility scores (0-1, where 1 = fully accessible)
    """
    # Try ViennaRNA
    try:
        import RNA
        rna_seq = seq.upper().replace("T", "U")
        fc = RNA.fold_compound(rna_seq)

        # Compute partition function
        fc.pf()

        # Get base pairing probabilities
        bpp = fc.bpp()
        n = len(seq)
        unpaired = [1.0] * (n + 1)  # 1-indexed

        for i in range(1, n + 1):
            for j in range(i + 1, min(n + 1, i + span + 1)):
                if bpp[i][j] > 0:
                    unpaired[i] -= bpp[i][j]
                    unpaired[j] -= bpp[i][j]

        return [max(0.0, min(1.0, p)) for p in unpaired[1:]]

    except ImportError:
        # Fallback: GC-content heuristic for accessibility
        n = len(seq)
        accessibility = [0.5] * n
        for i in range(n):
            start = max(0, i - window // 2)
            end = min(n, i + window // 2)
            window_seq = seq[start:end].upper()
            gc = sum(1 for b in window_seq if b in "GC") / max(1, len(window_seq))
            # Higher GC = more structure = less accessible
            accessibility[i] = max(0.0, min(1.0, 1.0 - gc))

        return accessibility


# ────────────────────────────────────────────────────────────
# 9. miRNA Binding Site Detection (with accessibility)
# ────────────────────────────────────────────────────────────

def _reverse_complement_dna(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement.get(b, "N") for b in reversed(seq.upper()))


def detect_mirna_sites(seq: str, tissue: str | None = None,
                        accessibility: list[float] | None = None) -> list[DegradationSignal]:
    """Detect miRNA binding sites in an mRNA sequence.

    Scans for seed-match complementarity (positions 2-8 of the miRNA)
    against the curated miRNA database.  When accessibility data is
    provided, sites in regions with accessibility < 0.05 have their
    severity reduced by 80% (the site is buried in secondary structure
    and inaccessible).

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        tissue: Optional tissue filter — only report miRNAs expressed
            in the specified tissue.
        accessibility: Optional per-position accessibility scores from
            :func:`compute_mrna_accessibility`.  When provided,
            inaccessible sites have reduced severity.

    Returns:
        List of DegradationSignal objects for detected miRNA binding sites.
    """
    signals: list[DegradationSignal] = []
    seq_upper = seq.upper()

    for mirna_id, mirna_data in _MIRNA_DATABASE.items():
        # Tissue filter
        if tissue is not None:
            if tissue.lower() not in [t.lower() for t in mirna_data["tissue"]] \
                    and "ubiquitous" not in [t.lower() for t in mirna_data["tissue"]]:
                continue

        seed = mirna_data["seed_2_8"]
        if not seed or len(seed) != 7:
            logger.warning(
                "Skipping miRNA %s: seed_2_8 is %d nt (expected 7)",
                mirna_id, len(seed) if seed else 0,
            )
            continue

        # The miRNA seed (positions 2-8) binds the mRNA target.
        # We search for the reverse complement of the seed on the
        # mRNA coding strand (the seed in the miRNA pairs with the
        # target, so we look for the complement).
        seed_rc = _reverse_complement_dna(seed)

        # Search for seed match
        for match in re.finditer(seed_rc, seq_upper):
            # Base severity from conservation
            conservation = mirna_data.get("conservation", "medium")
            if conservation == "high":
                base_severity = 0.7
            elif conservation == "medium":
                base_severity = 0.5
            else:
                base_severity = 0.3

            # Position-dependent modifier: 3'UTR sites are more potent
            frac = match.start() / max(len(seq_upper), 1)
            if frac > 0.85:
                position_mult = 1.3  # 3'UTR
            elif frac > 0.75:
                position_mult = 1.1  # Near 3'UTR
            else:
                position_mult = 0.8  # CDS — less effective due to ribosomes

            severity = min(1.0, base_severity * position_mult)

            # Accessibility-based severity adjustment
            if accessibility is not None:
                # Check the accessibility of the binding site region
                site_start = match.start()
                site_end = match.end()
                if site_end <= len(accessibility):
                    site_accessibility = sum(
                        accessibility[site_start:site_end]
                    ) / max(1, site_end - site_start)

                    # If accessibility < 0.05, reduce severity by 80%
                    if site_accessibility < 0.05:
                        severity *= 0.2
                        accessibility_note = " (80% reduced: site buried in structure)"
                    else:
                        accessibility_note = f" (accessibility={site_accessibility:.2f})"
                else:
                    accessibility_note = ""
            else:
                accessibility_note = ""

            signals.append(DegradationSignal(
                signal_type="mirna",
                position=match.start(),
                sequence_context=match.group(),
                severity=severity,
                description=(
                    f"miRNA {mirna_data['name']} seed match at position "
                    f"{match.start()}: {match.group()} (seed RC of "
                    f"{seed}){accessibility_note}"
                ),
                pathway="mirna",
            ))

    return signals


# ────────────────────────────────────────────────────────────
# 10. Bacterial RNase E Detection (retained for prokaryotic targets)
# ────────────────────────────────────────────────────────────

# RNase E recognition motifs for prokaryotic mRNA design.
# These are AT-rich single-stranded regions preferred by RNase E.

_RNASE_E_PATTERNS: list[tuple[str, str, float]] = [
    # (pattern, description, base_severity)
    (r"[AT]{5,}", "RNase E recognition site (AT-rich single-stranded region)", 0.4),
    (r"A[AT]T[AT]A[AT]T", "RNase E cleavage consensus (AT-rich motif)", 0.5),
    (r"ATTA[AT]", "Endonuclease cleavage site (AUUA[AU] motif)", 0.35),
]


def detect_rnase_e_sites(seq: str) -> list[DegradationSignal]:
    """Detect RNase E recognition sites in a bacterial mRNA sequence.

    RNase E is the major bacterial endonuclease that initiates mRNA
    decay.  This function is retained for prokaryotic expression
    design; for eukaryotic targets, use
    :func:`detect_eukaryotic_decay_signals` instead.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).

    Returns:
        List of DegradationSignal objects for RNase E sites.
    """
    signals: list[DegradationSignal] = []
    seq_upper = seq.upper()

    for pattern, description, base_severity in _RNASE_E_PATTERNS:
        for match in re.finditer(pattern, seq_upper):
            signals.append(DegradationSignal(
                signal_type="rnase_e",
                position=match.start(),
                sequence_context=match.group(),
                severity=base_severity,
                description=f"{description} at position {match.start()}: {match.group()}",
                pathway="rnase_e",
            ))

    return signals


# ────────────────────────────────────────────────────────────
# 11. Comprehensive Degradation Analysis
# ────────────────────────────────────────────────────────────


def analyze_rna_degradation(
    seq: str,
    organism: str = "Homo_sapiens",
    tissue: str | None = None,
    are_class: str = "none",
    polyA_tail_length: int = 200,
    orf_start: int = 0,
    orf_end: int | None = None,
    has_introns: bool = False,
    last_exon_junction: int | None = None,
    use_accessibility: bool = True,
) -> list[DegradationSignal]:
    """Run comprehensive RNA degradation signal analysis.

    Combines all degradation detection modules into a single analysis
    pipeline.  For eukaryotic organisms, the analysis includes miRNA
    sites, ARE signals, m6A, eukaryotic decay pathways, NMD triggers,
    and epitranscriptomic modifications.  For prokaryotic organisms,
    RNase E detection is used instead of eukaryotic decay pathways.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        organism: Target organism name (e.g., ``"Homo_sapiens"``,
            ``"Escherichia_coli"``).
        tissue: Optional tissue filter for miRNA analysis.
        are_class: Pre-detected ARE class for eukaryotic decay modeling
            (``"I"``, ``"II"``, ``"III"``, or ``"none"``).
        polyA_tail_length: PolyA tail length for eukaryotic decay modeling.
        orf_start: Start position of main ORF (0-indexed).
        orf_end: End position of main ORF (0-indexed, exclusive).
        has_introns: Whether the gene has introns (affects NMD detection).
        last_exon_junction: Position of last exon-exon junction.
        use_accessibility: Whether to compute and use mRNA accessibility
            for miRNA severity adjustment (default True).

    Returns:
        Combined list of DegradationSignal objects from all detectors,
        sorted by severity (highest first).
    """
    all_signals: list[DegradationSignal] = []

    # Determine if organism is prokaryotic
    _PROKARYOTIC_ORGANISMS = {
        "Escherichia_coli", "E_coli", "Bacillus_subtilis",
        "Pseudomonas_aeruginosa", "Corynebacterium_glutamicum",
    }
    is_prokaryote = organism in _PROKARYOTIC_ORGANISMS

    # Compute accessibility if requested
    accessibility: list[float] | None = None
    if use_accessibility and not is_prokaryote:
        try:
            accessibility = compute_mrna_accessibility(seq)
        except Exception:
            logger.debug("Accessibility computation failed, proceeding without it",
                         exc_info=True)
            accessibility = None

    # ARE detection (all organisms — ARE-like motifs can affect
    # bacterial mRNA too, though via different mechanisms)
    all_signals.extend(detect_are_signals(seq))

    # Auto-detect ARE class from the detected signals
    if are_class == "none":
        are_signals = [s for s in all_signals if s.signal_type.startswith("are_class_")]
        if any(s.signal_type == "are_class_ii" for s in are_signals):
            are_class = "II"
        elif any(s.signal_type == "are_class_i" for s in are_signals):
            are_class = "I"
        elif any(s.signal_type == "are_class_iii" for s in are_signals):
            are_class = "III"

    # m6A detection
    all_signals.extend(detect_m6a_sites(seq))

    # Epitranscriptomic marks (m5C, Ψ, m1A, 2'-O-Me)
    all_signals.extend(detect_epitranscriptomic_marks(seq))

    if is_prokaryote:
        # Prokaryotic pathway: RNase E
        all_signals.extend(detect_rnase_e_sites(seq))
    else:
        # Eukaryotic pathways
        # miRNA detection with accessibility
        all_signals.extend(detect_mirna_sites(seq, tissue=tissue,
                                               accessibility=accessibility))

        # Eukaryotic decay pathways
        all_signals.extend(detect_eukaryotic_decay_signals(
            seq, are_class=are_class, polyA_tail_length=polyA_tail_length,
        ))

        # NMD detection
        all_signals.extend(detect_nmd_triggers(
            seq, orf_start=orf_start, orf_end=orf_end,
            has_introns=has_introns, last_exon_junction=last_exon_junction,
        ))

    # Sort by severity (highest first)
    all_signals.sort(key=lambda s: s.severity, reverse=True)

    return all_signals
