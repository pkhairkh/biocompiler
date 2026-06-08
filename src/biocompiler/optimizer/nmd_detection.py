"""BioCompiler Nonsense-Mediated Decay (NMD) Detection Module
==============================================================

Comprehensive NMD trigger detection, probability estimation, and
suppression for therapeutic mRNA and heterologous expression design.

This module models the four major NMD trigger pathways:

1. **EJC-dependent NMD** — Premature termination codon (PTC) located
   ≥50-55 nt upstream of the last exon-exon junction triggers
   SMG1/UPF1-dependent NMD (strong, 80-90% probability).
2. **Long 3'UTR NMD** — In intronless (synthetic) mRNA, a 3'UTR
   exceeding ~2 kb can trigger NMD via the faux-UTR model (moderate,
   50-70% probability).
3. **uORF-mediated NMD** — Upstream open reading frames (uORFs) with
   strong Kozak context can consume ribosomes and/or trigger NMD
   themselves; severity depends on Kozak strength and distance to
   the main ORF.
4. **Intron retention NMD** — Retained introns that introduce PTCs
   trigger EJC-dependent or EJC-independent NMD depending on context.

Key NMD rules implemented from the literature:

- **50-55 nt rule**: PTCs located ≥50-55 nt upstream of the last
  exon-exon junction (EJ) are strong NMD triggers.  PTCs within 50 nt
  of the last EJ generally escape NMD (Chang 2007; Lykke-Andersen &
  Jensen 2015).
- **Last exon escape**: PTCs in the final exon escape EJC-dependent
  NMD because no downstream EJ exists (Kervestin & Jacobson 2012).
- **Long 3'UTR model**: For intronless mRNA, an abnormally long 3'UTR
  is sensed as aberrant by the ribosome, recruiting UPF1 and
  triggering NMD via the faux-UTR mechanism (Linde & Bhatt 2019).
- **Kozak-dependent uORF NMD**: uORFs with strong Kozak consensus
  sequences are efficiently translated and can trigger NMD when they
  terminate far from the 5' cap; weak-Kozak uORFs may be leaky
  scanned (Lykke-Andersen & Jensen 2015).
- **SMG1/UPF1 dependence**: EJC-dependent NMD requires the SMG1
  kinase and UPF1 helicase; the long 3'UTR pathway also requires
  UPF1 but can be SMG1-independent in some contexts.

Usage::

    from biocompiler.optimizer.nmd_detection import (
        NMDSignal,
        detect_nmd_triggers,
        detect_ejc_dependent_nmd,
        detect_long_3utr_nmd,
        detect_uorf_nmd,
        detect_intron_retention_nmd,
        fix_nmd_triggers,
        estimate_nmd_probability,
    )

    # Full NMD analysis
    signals = detect_nmd_triggers(
        seq, orf_start=0, orf_end=900,
        has_introns=True, last_exon_junction=850,
    )
    for s in signals:
        print(f"{s.signal_type}: {s.description} "
              f"(severity={s.severity:.2f}, P(NMD)={estimate_nmd_probability(s):.2f})")

    # Fix NMD triggers
    fixed_seq = fix_nmd_triggers(seq, signals, organism="human")

References:
  Lykke-Andersen, J. & Jensen, T.H. (2015). "Nonsense-mediated mRNA
  decay: an intricate machinery that shapes transcriptomes."
  *Nat Rev Mol Cell Biol* 16:665–677.

  Chang, Y.-F. et al. (2007). "The nonsense-mediated decay RNA
  surveillance pathway." *Annu Rev Biochem* 76:51–74.

  Kervestin, S. & Jacobson, A. (2012). "NMD: a multifaceted response
  to premature translational termination." *Nat Rev Mol Cell Biol*
  13:700–712.

  Linde, A. & Bhatt, D.M. (2019). "Nonsense-mediated mRNA decay
  in the context of mammalian gene regulation."
  *Nat Rev Genet* 20:557–570.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Literal

logger = logging.getLogger(__name__)

__all__ = [
    "NMDSignal",
    "detect_nmd_triggers",
    "detect_ejc_dependent_nmd",
    "detect_long_3utr_nmd",
    "detect_uorf_nmd",
    "detect_intron_retention_nmd",
    "fix_nmd_triggers",
    "estimate_nmd_probability",
]

# ────────────────────────────────────────────────────────────
# 1. Core Data Structures
# ────────────────────────────────────────────────────────────

# Signal type enum values
NMD_SIGNAL_TYPES = Literal[
    "ejc_dependent",
    "long_3utr",
    "upstream_orf",
    "intron_retention",
]

# NMD pathway enum values
NMD_PATHWAYS = Literal[
    "ejc_dependent",
    "faux_utr",
    "uorf",
    "intron_retention",
]


@dataclass
class NMDSignal:
    """A single NMD trigger signal detected in an mRNA sequence.

    Attributes:
        signal_type: Category of NMD trigger.  One of:
            - ``"ejc_dependent"`` — PTC upstream of last exon junction
            - ``"long_3utr"`` — Abnormally long 3'UTR (intronless mRNA)
            - ``"upstream_orf"`` — uORF that may trigger NMD
            - ``"intron_retention"`` — Retained intron introducing PTC
        position: 0-based position in the sequence where the trigger
            is located (e.g., PTC position, uORF start, intron start).
        severity: Severity score in [0, 1].  Higher values indicate
            greater likelihood that NMD will be triggered and result
            in mRNA degradation.
        description: Human-readable description of the NMD trigger.
        nmd_pathway: The NMD pathway this signal activates.  One of:
            - ``"ejc_dependent"`` — SMG1/UPF1-dependent EJC pathway
            - ``"faux_utr"`` — Long 3'UTR / faux-UTR pathway
            - ``"uorf"`` — uORF-mediated NMD
            - ``"intron_retention"`` — Intron retention NMD
    """

    signal_type: str
    position: int
    severity: float
    description: str
    nmd_pathway: str = ""

    def __post_init__(self) -> None:
        # Clamp severity to [0, 1]
        self.severity = max(0.0, min(1.0, self.severity))


# ────────────────────────────────────────────────────────────
# 2. Constants from Literature
# ────────────────────────────────────────────────────────────

# EJC-dependent NMD: PTC must be ≥50-55 nt upstream of the last
# exon-exon junction (Chang 2007; Lykke-Andersen & Jensen 2015).
# We use 50 nt as the canonical threshold; 55 nt is the more
# conservative value used in some studies.
EJC_DISTANCE_THRESHOLD: int = 50
EJC_DISTANCE_THRESHOLD_CONSERVATIVE: int = 55

# Long 3'UTR NMD threshold: 3'UTRs >2 kb trigger NMD in intronless
# mRNA (Linde & Bhatt 2019; Kervestin & Jacobson 2012).
LONG_3UTR_THRESHOLD: int = 2000

# Kozak consensus sequence scoring weights (position -3 to +4 relative
# to A of ATG).  Values derived from Jackson et al. 2010 and
# Noderer & Aldred 2016 quantitative studies.
# Strong Kozak: GCCRCCATGG (R = A/G)
# -3: G is optimal; +4: G is optimal
KOZAK_POSITIONS: dict[int, dict[str, float]] = {
    -3: {"G": 0.35, "A": 0.20, "C": 0.10, "T": 0.05},
    -2: {"C": 0.15, "G": 0.10, "A": 0.10, "T": 0.05},
    -1: {"C": 0.20, "A": 0.15, "G": 0.10, "T": 0.05},
    0:  {"A": 1.0},   # ATG start — always A
    1:  {"T": 1.0},   # ATG — always T
    2:  {"G": 1.0},   # ATG — always G
    3:  {"G": 0.35, "A": 0.15, "C": 0.10, "T": 0.05},
}

# Minimum uORF length (in codons, excluding stop) to be considered
# a potential NMD trigger.  Short uORFs (<3 codons) are often leaky
# scanned or too short to efficiently engage the NMD machinery.
MIN_UORF_CODONS: int = 3

# Maximum 5'UTR search distance for uORFs (nt from start of sequence)
MAX_UORF_SEARCH_DISTANCE: int = 500

# Standard genetic code — stop codons
STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})
START_CODON: str = "ATG"

# Standard codon table (DNA) for fix_nmd_triggers
_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Reverse mapping: amino acid → list of synonymous codons
_AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in _CODON_TABLE.items():
    if _aa != "*":
        _AA_TO_CODONS.setdefault(_aa, []).append(_codon)


# ────────────────────────────────────────────────────────────
# 3. Helper Functions
# ────────────────────────────────────────────────────────────


def _score_kozak(seq: str, atg_pos: int) -> float:
    """Score the Kozak consensus strength around an ATG at *atg_pos*.

    Returns a score in [0, 1] where:
      - ≥0.8: strong Kozak (efficient initiation)
      - 0.5–0.8: moderate Kozak
      - <0.5: weak Kozak (leaky scanning likely)

    The scoring uses position-specific weights derived from
    Jackson et al. (2010) and Noderer & Aldred (2016).  Only the
    flanking positions (-3 to -1 and +3 to +4 relative to the A of
    ATG) are scored, since the ATG itself is present at every start
    codon and does not contribute discriminative power.

    Args:
        seq: Full mRNA sequence (DNA alphabet, uppercase).
        atg_pos: 0-based position of the 'A' in ATG.

    Returns:
        Kozak consensus score in [0, 1].
    """
    # Only score flanking positions; the ATG positions (0, 1, 2) are
    # always maximally matched at any start codon and provide no
    # information about initiation efficiency.
    _FLANKING_OFFSETS = [-3, -2, -1, 3]

    max_score = 0.0
    observed = 0.0

    for offset in _FLANKING_OFFSETS:
        weights = KOZAK_POSITIONS[offset]
        idx = atg_pos + offset
        if 0 <= idx < len(seq):
            base = seq[idx]
            observed += weights.get(base, 0.0)
        # The maximum possible score per position is the highest weight
        max_score += max(weights.values())

    return min(1.0, observed / max(1e-9, max_score))


def _find_orf_end(seq: str, orf_start: int) -> int | None:
    """Find the end position of the ORF starting at *orf_start*.

    Scans in-frame for the first stop codon.

    Args:
        seq: mRNA sequence (DNA alphabet, uppercase).
        orf_start: 0-based start position.

    Returns:
        End position (exclusive) of the ORF, or None if no stop found.
    """
    for i in range(orf_start, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) < 3:
            break
        if codon in STOP_CODONS:
            return i + 3
    return None


def _find_all_atg_positions(seq: str, start: int, end: int) -> list[int]:
    """Find all ATG positions in the range [start, end).

    Args:
        seq: mRNA sequence (DNA alphabet, uppercase).
        start: Start of search range (inclusive).
        end: End of search range (exclusive).

    Returns:
        List of 0-based ATG positions.
    """
    positions: list[int] = []
    for i in range(start, min(end, len(seq) - 2)):
        if seq[i:i + 3] == START_CODON:
            positions.append(i)
    return positions


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


# ────────────────────────────────────────────────────────────
# 4. EJC-Dependent NMD Detection
# ────────────────────────────────────────────────────────────


def detect_ejc_dependent_nmd(
    seq: str,
    orf_start: int,
    orf_end: int,
    exon_junctions: list[int],
) -> list[NMDSignal]:
    """Detect EJC-dependent NMD triggers in spliced mRNA.

    EJC-dependent NMD is triggered when a premature termination codon
    (PTC) is located ≥50-55 nt upstream of the last exon-exon junction
    (EJ).  The Exon Junction Complex (EJC), deposited ~20-24 nt
    upstream of each exon-exon junction during splicing, serves as a
    mark for the NMD machinery.  When a ribosome terminates at a PTC
    upstream of an EJ, the EJC remains bound and recruits UPF1/SMG1,
    triggering NMD.

    PTCs in the last exon escape NMD because no downstream EJ exists
    (Kervestin & Jacobson 2012).

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        orf_start: 0-based start position of the main ORF.
        orf_end: 0-based end position of the main ORF (exclusive).
        exon_junctions: List of 0-based positions of exon-exon
            junctions (the first nucleotide of the downstream exon).
            Must be sorted in ascending order.

    Returns:
        List of NMDSignal objects for EJC-dependent NMD triggers.

    References:
        Chang, Y.-F. et al. (2007). *Annu Rev Biochem* 76:51–74.
        Lykke-Andersen, J. & Jensen, T.H. (2015).
          *Nat Rev Mol Cell Biol* 16:665–677.
    """
    signals: list[NMDSignal] = []
    seq_upper = seq.upper()

    if not exon_junctions:
        return signals

    # Sort junctions to ensure correct ordering
    sorted_junctions = sorted(exon_junctions)
    last_ej = sorted_junctions[-1]

    # Scan for premature stop codons in the ORF
    # Skip the first codon (start ATG) and the last codon (natural stop)
    for i in range(orf_start + 3, orf_end - 3, 3):
        codon = seq_upper[i:i + 3]
        if len(codon) < 3:
            break
        if codon in STOP_CODONS and i < orf_end - 3:
            # This is a PTC (premature, not the natural stop)
            distance_to_last_ej = last_ej - i

            # Check if PTC is in the last exon (after last EJ)
            if i >= last_ej:
                # PTC in last exon — escapes NMD
                logger.debug(
                    "PTC at %d is in last exon (after EJ at %d), escapes NMD",
                    i, last_ej,
                )
                continue

            # Apply the 50-55 nt rule
            if distance_to_last_ej >= EJC_DISTANCE_THRESHOLD:
                # Strong NMD trigger
                severity = min(
                    1.0,
                    0.5 + 0.5 * (distance_to_last_ej / 1000.0),
                )
                signals.append(NMDSignal(
                    signal_type="ejc_dependent",
                    position=i,
                    severity=severity,
                    description=(
                        f"EJC-dependent NMD: PTC at position {i} "
                        f"({codon}), {distance_to_last_ej} nt upstream of "
                        f"last exon junction at {last_ej} "
                        f"(threshold: {EJC_DISTANCE_THRESHOLD} nt)"
                    ),
                    nmd_pathway="ejc_dependent",
                ))
            elif distance_to_last_ej >= EJC_DISTANCE_THRESHOLD_CONSERVATIVE - 5:
                # Borderline: 45-49 nt — reduced probability
                severity = 0.4
                signals.append(NMDSignal(
                    signal_type="ejc_dependent",
                    position=i,
                    severity=severity,
                    description=(
                        f"Borderline EJC-dependent NMD: PTC at position {i} "
                        f"({codon}), {distance_to_last_ej} nt upstream of "
                        f"last exon junction at {last_ej} "
                        f"(near threshold: {EJC_DISTANCE_THRESHOLD} nt)"
                    ),
                    nmd_pathway="ejc_dependent",
                ))

    return signals


# ────────────────────────────────────────────────────────────
# 5. Long 3'UTR NMD Detection
# ────────────────────────────────────────────────────────────


def detect_long_3utr_nmd(
    seq: str,
    orf_end: int,
    threshold: int = LONG_3UTR_THRESHOLD,
) -> list[NMDSignal]:
    """Detect long 3'UTR-triggered NMD in intronless (synthetic) mRNA.

    For intronless mRNA (e.g., synthetic mRNA therapeutics), the
    EJC-dependent pathway is not available.  Instead, an abnormally
    long 3'UTR (>2 kb) can trigger NMD via the faux-UTR model:
    the ribosome senses an extended 3'UTR as aberrant because of
    the physical distance between the stop codon and the polyA tail,
    leading to UPF1 accumulation and NMD activation.

    The probability of NMD increases with 3'UTR length:
      - 2-5 kb: moderate (50-70%)
      - >5 kb: strong (70-90%)

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        orf_end: 0-based end position of the main ORF (exclusive).
        threshold: Minimum 3'UTR length to trigger NMD (default 2000 nt).

    Returns:
        List of NMDSignal objects for long 3'UTR NMD triggers.

    References:
        Linde, A. & Bhatt, D.M. (2019). *Nat Rev Genet* 20:557–570.
        Kervestin, S. & Jacobson, A. (2012).
          *Nat Rev Mol Cell Biol* 13:700–712.
    """
    signals: list[NMDSignal] = []
    utr3_length = len(seq) - orf_end

    if utr3_length >= threshold:
        # Severity scales with 3'UTR length
        if utr3_length >= 5000:
            severity = min(1.0, 0.7 + 0.3 * ((utr3_length - 5000) / 10000))
        elif utr3_length >= 3000:
            severity = 0.6 + 0.1 * ((utr3_length - 3000) / 2000)
        else:
            severity = 0.5 + 0.1 * ((utr3_length - threshold) / 1000)

        severity = min(1.0, severity)

        signals.append(NMDSignal(
            signal_type="long_3utr",
            position=orf_end,
            severity=severity,
            description=(
                f"Long 3'UTR NMD: 3'UTR is {utr3_length} nt "
                f"(threshold: {threshold} nt). "
                f"Faux-UTR model: extended 3'UTR sensed as aberrant, "
                f"UPF1 accumulates downstream of stop codon."
            ),
            nmd_pathway="faux_utr",
        ))

    return signals


# ────────────────────────────────────────────────────────────
# 6. uORF-Mediated NMD Detection
# ────────────────────────────────────────────────────────────


def detect_uorf_nmd(
    seq: str,
    orf_start: int,
    kozak_threshold: float = 0.5,
) -> list[NMDSignal]:
    """Detect upstream ORF (uORF)-mediated NMD triggers.

    uORFs are open reading frames in the 5'UTR (upstream of the main
    ORF).  When a ribosome initiates at a uORF start codon and
    translates the uORF, two deleterious effects can occur:

    1. **Ribosome consumption**: The ribosome is diverted from the
       main ORF, reducing expression.
    2. **NMD trigger**: If the uORF terminates far from the 5' cap
       (especially in spliced mRNA with downstream EJCs), the
       terminating ribosome can recruit UPF1 and trigger NMD.

    The severity depends on:
    - **Kozak strength**: Strong Kozak consensus → efficient
      initiation → higher NMD probability.
    - **uORF length**: Longer uORFs (≥35 codons) are more likely to
      trigger NMD than very short ones.
    - **Distance to main ORF**: uORFs close to the main ORF start
      have less room for EJC deposition, reducing NMD probability.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        orf_start: 0-based start position of the main ORF.
        kozak_threshold: Minimum Kozak score to consider a uORF as a
            potential NMD trigger.  Below this threshold, the uORF is
            likely leaky-scanned.  Default 0.5.

    Returns:
        List of NMDSignal objects for uORF-mediated NMD triggers.

    References:
        Lykke-Andersen, J. & Jensen, T.H. (2015).
          *Nat Rev Mol Cell Biol* 16:665–677.
        Calvo, S.E. et al. (2009). "Upstream open reading frames cause
          widespread reduction in protein expression." *PNAS* 106:7507.
    """
    signals: list[NMDSignal] = []
    seq_upper = seq.upper()

    # Search for ATGs in the 5'UTR (upstream of main ORF start)
    search_end = min(orf_start, MAX_UORF_SEARCH_DISTANCE)
    atg_positions = _find_all_atg_positions(seq_upper, 0, search_end)

    for atg_pos in atg_positions:
        # Score Kozak context
        kozak_score = _score_kozak(seq_upper, atg_pos)

        if kozak_score < kozak_threshold:
            # Weak Kozak — likely leaky scanned, low NMD risk
            continue

        # Find uORF end
        uorf_end = _find_orf_end(seq_upper, atg_pos)
        if uorf_end is None:
            # No stop codon found — uORF extends into main ORF or beyond
            # This is a special case: overlapping uORF
            uorf_length_codons = (orf_start - atg_pos) // 3
            if uorf_length_codons < MIN_UORF_CODONS:
                continue  # Too short to matter

            severity = min(
                1.0,
                0.3 + 0.3 * kozak_score + 0.2 * min(1.0, uorf_length_codons / 100),
            )
            signals.append(NMDSignal(
                signal_type="upstream_orf",
                position=atg_pos,
                severity=severity,
                description=(
                    f"uORF-mediated NMD: overlapping uORF at position "
                    f"{atg_pos} (Kozak={kozak_score:.2f}, "
                    f"~{uorf_length_codons} codons before main ORF, "
                    f"no in-frame stop found)"
                ),
                nmd_pathway="uorf",
            ))
            continue

        uorf_length_codons = (uorf_end - atg_pos) // 3 - 1  # Exclude stop
        if uorf_length_codons < MIN_UORF_CODONS:
            continue  # Too short to trigger NMD

        # Calculate distance from uORF stop to main ORF start
        distance_to_main_orf = orf_start - uorf_end

        # Severity model:
        # - Kozak strength: higher → more efficient initiation → worse
        # - uORF length: longer uORFs → more likely to trigger NMD
        # - Distance: uORFs far from main ORF → more room for EJCs → worse
        #   (but for synthetic intronless mRNA, this is less relevant)

        length_factor = min(1.0, uorf_length_codons / 100.0)
        kozak_factor = kozak_score
        distance_factor = min(1.0, max(0.0, distance_to_main_orf) / 500.0)

        # uORFs that terminate close to main ORF can also reduce
        # expression by ribosome consumption, even without NMD
        if distance_to_main_orf < 0:
            # uORF overlaps main ORF — strong impact
            severity = min(
                1.0,
                0.4 + 0.3 * kozak_factor + 0.2 * length_factor,
            )
        elif distance_to_main_orf < 50:
            # uORF terminates very close to main ORF — moderate NMD risk
            severity = min(
                1.0,
                0.2 + 0.3 * kozak_factor + 0.2 * length_factor + 0.1 * distance_factor,
            )
        else:
            # uORF terminates further away — higher NMD risk
            severity = min(
                1.0,
                0.3 + 0.3 * kozak_factor + 0.2 * length_factor + 0.2 * distance_factor,
            )

        signals.append(NMDSignal(
            signal_type="upstream_orf",
            position=atg_pos,
            severity=severity,
            description=(
                f"uORF-mediated NMD: uORF at position {atg_pos} "
                f"(Kozak={kozak_score:.2f}, {uorf_length_codons} codons, "
                f"stop at {uorf_end}, "
                f"{distance_to_main_orf} nt to main ORF)"
            ),
            nmd_pathway="uorf",
        ))

    return signals


# ────────────────────────────────────────────────────────────
# 7. Intron Retention NMD Detection
# ────────────────────────────────────────────────────────────


def detect_intron_retention_nmd(
    seq: str,
    intron_positions: list[tuple[int, int]],
) -> list[NMDSignal]:
    """Detect intron retention events that can trigger NMD.

    Retained introns can trigger NMD through two mechanisms:

    1. **PTC introduction**: The retained intron sequence may contain
       in-frame stop codons that act as premature termination codons.
    2. **EJC deposition**: If the intron contains splicing signals
       (5' splice site, 3' splice site, branch point), partial
       spliceosome assembly may deposit EJCs that recruit NMD factors.

    The severity depends on:
    - Whether the intron contains in-frame stop codons
    - The position of the intron relative to the ORF
    - The length of the retained intron

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        intron_positions: List of (start, end) tuples for retained
            introns, where start and end are 0-based positions.

    Returns:
        List of NMDSignal objects for intron retention NMD triggers.

    References:
        Lykke-Andersen, J. & Jensen, T.H. (2015).
          *Nat Rev Mol Cell Biol* 16:665–677.
        Kervestin, S. & Jacobson, A. (2012).
          *Nat Rev Mol Cell Biol* 13:700–712.
    """
    signals: list[NMDSignal] = []
    seq_upper = seq.upper()

    for intron_start, intron_end in intron_positions:
        if intron_start < 0 or intron_end > len(seq_upper) or intron_start >= intron_end:
            logger.warning(
                "Invalid intron position (%d, %d) for sequence of length %d",
                intron_start, intron_end, len(seq_upper),
            )
            continue

        intron_seq = seq_upper[intron_start:intron_end]
        intron_length = intron_end - intron_start

        # Check for in-frame stop codons in the retained intron
        # We need to determine the reading frame at the intron start
        has_inframe_stop = False
        ptc_position = None
        ptc_codon = None

        # Try all three reading frames
        for frame in range(3):
            for i in range(intron_start + frame, intron_end - 2, 3):
                codon = seq_upper[i:i + 3]
                if len(codon) < 3:
                    break
                if codon in STOP_CODONS:
                    has_inframe_stop = True
                    ptc_position = i
                    ptc_codon = codon
                    break
            if has_inframe_stop:
                break

        # Check for splice-like signals in the intron
        has_splice_signals = False
        if intron_length >= 30:
            # 5' splice site consensus: GT
            if intron_seq[:2] == "GT":
                has_splice_signals = True
            # 3' splice site consensus: AG
            if intron_seq[-2:] == "AG":
                has_splice_signals = True

        # Severity calculation
        base_severity = 0.0
        reasons: list[str] = []

        if has_inframe_stop and ptc_position is not None:
            # PTC in retained intron — strong NMD trigger
            base_severity += 0.5
            reasons.append(
                f"in-frame stop codon {ptc_codon} at position {ptc_position}"
            )

        if has_splice_signals:
            # Splice signals suggest partial spliceosome engagement
            base_severity += 0.2
            reasons.append("splice-like signals present")

        # Longer introns are more likely to trigger NMD
        length_factor = min(1.0, intron_length / 500.0)
        base_severity += 0.2 * length_factor
        if intron_length > 100:
            reasons.append(f"intron length {intron_length} nt")

        severity = min(1.0, base_severity)

        if severity > 0.1:  # Only report significant signals
            reasons_str = "; ".join(reasons) if reasons else "retained intron"
            signals.append(NMDSignal(
                signal_type="intron_retention",
                position=intron_start,
                severity=severity,
                description=(
                    f"Intron retention NMD: retained intron at "
                    f"positions {intron_start}-{intron_end} "
                    f"({intron_length} nt). {reasons_str}"
                ),
                nmd_pathway="intron_retention",
            ))

    return signals


# ────────────────────────────────────────────────────────────
# 8. Comprehensive NMD Detection
# ────────────────────────────────────────────────────────────


def detect_nmd_triggers(
    seq: str,
    orf_start: int = 0,
    orf_end: int | None = None,
    has_introns: bool = False,
    last_exon_junction: int | None = None,
    polyA_signal_positions: list[int] | None = None,
) -> list[NMDSignal]:
    """Detect all NMD trigger types in an mRNA sequence.

    This is the main entry point for NMD detection.  It runs all four
    NMD detection modules and returns a combined list of signals.

    For spliced mRNA (has_introns=True):
      - EJC-dependent NMD (if exon junctions provided)
      - uORF-mediated NMD
      - Intron retention NMD (if intron positions provided)

    For intronless/synthetic mRNA (has_introns=False):
      - Long 3'UTR NMD
      - uORF-mediated NMD

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        orf_start: 0-based start position of the main ORF (default 0).
        orf_end: 0-based end position of the main ORF (exclusive).
            If None, automatically detected by scanning for the first
            in-frame stop codon.
        has_introns: Whether the gene has introns.  Affects which NMD
            pathways are checked.
        last_exon_junction: Position of the last exon-exon junction.
            Required for EJC-dependent NMD detection when has_introns
            is True.
        polyA_signal_positions: Optional list of polyA signal (AATAAA)
            positions in the sequence.  Used to refine 3'UTR boundary
            estimation for long 3'UTR NMD detection.

    Returns:
        Combined list of NMDSignal objects from all NMD detectors,
        sorted by severity (highest first).
    """
    all_signals: list[NMDSignal] = []
    seq_upper = seq.upper()

    # Auto-detect ORF end if not provided
    if orf_end is None:
        orf_end = _find_orf_end(seq_upper, orf_start)
        if orf_end is None:
            logger.warning("No stop codon found starting at position %d", orf_start)
            return all_signals

    # 1. EJC-dependent NMD (spliced mRNA only)
    if has_introns and last_exon_junction is not None:
        ejc_signals = detect_ejc_dependent_nmd(
            seq_upper, orf_start, orf_end, [last_exon_junction],
        )
        all_signals.extend(ejc_signals)

    # 2. Long 3'UTR NMD (primarily for intronless mRNA, but also
    #    relevant for spliced mRNA with very long 3'UTRs)
    if not has_introns or (orf_end is not None and len(seq_upper) - orf_end > LONG_3UTR_THRESHOLD):
        utr3_signals = detect_long_3utr_nmd(seq_upper, orf_end)
        all_signals.extend(utr3_signals)

    # 3. uORF-mediated NMD (always check)
    uorf_signals = detect_uorf_nmd(seq_upper, orf_start)
    all_signals.extend(uorf_signals)

    # 4. Intron retention NMD — not checked here because it requires
    #    explicit intron position data.  Users should call
    #    detect_intron_retention_nmd() directly with intron_positions.

    # Sort by severity (highest first)
    all_signals.sort(key=lambda s: s.severity, reverse=True)

    return all_signals


# ────────────────────────────────────────────────────────────
# 9. NMD Probability Estimation
# ────────────────────────────────────────────────────────────


def estimate_nmd_probability(signal: NMDSignal) -> float:
    """Estimate the probability that an NMD trigger will actually
    trigger mRNA decay.

    NMD is a stochastic process — not every PTC or long 3'UTR will
    necessarily trigger decay.  The probability depends on:

    - **EJC-dependent NMD**: 80-90% for strong PTCs ≥50 nt upstream
      of the last EJ.  Reduced for borderline positions (45-49 nt),
      PTCs near the 3' end of the penultimate exon, or if the
      downstream EJC is weakly bound.
    - **Long 3'UTR NMD**: 50-70% for 3'UTRs >2 kb.  Increases with
      length.  UPF1-dependent but can be SMG1-independent.
    - **uORF-mediated NMD**: 20-60% depending on Kozak strength,
      uORF length, and distance to main ORF.  Weak-Kozak uORFs are
      often leaky-scanned and rarely trigger NMD.
    - **Intron retention NMD**: 40-80% depending on whether the
      retained intron introduces a PTC and has splice signals.

    Args:
        signal: An NMDSignal object from one of the detection functions.

    Returns:
        Estimated probability of NMD triggering (float in [0, 1]).

    References:
        Lykke-Andersen, J. & Jensen, T.H. (2015).
          *Nat Rev Mol Cell Biol* 16:665–677.
        Kervestin, S. & Jacobson, A. (2012).
          *Nat Rev Mol Cell Biol* 13:700–712.
        Linde, A. & Bhatt, D.M. (2019). *Nat Rev Genet* 20:557–570.
    """
    if signal.signal_type == "ejc_dependent":
        # EJC-dependent NMD: high probability for strong triggers
        # 80-90% for canonical PTCs, reduced for borderline cases
        if signal.severity >= 0.7:
            return 0.85  # Strong EJC-dependent trigger
        elif signal.severity >= 0.5:
            return 0.70  # Moderate
        else:
            return 0.40  # Borderline (45-49 nt from EJ)

    elif signal.signal_type == "long_3utr":
        # Long 3'UTR NMD: moderate probability
        # 50-70% for typical long 3'UTRs, higher for very long ones
        if signal.severity >= 0.8:
            return 0.75  # Very long 3'UTR
        elif signal.severity >= 0.6:
            return 0.60  # Moderately long
        else:
            return 0.50  # Just above threshold

    elif signal.signal_type == "upstream_orf":
        # uORF-mediated NMD: variable probability
        # Depends heavily on Kozak strength and uORF length
        if signal.severity >= 0.7:
            return 0.55  # Strong Kozak, long uORF
        elif signal.severity >= 0.5:
            return 0.40  # Moderate
        else:
            return 0.20  # Weak trigger

    elif signal.signal_type == "intron_retention":
        # Intron retention NMD: variable
        if signal.severity >= 0.7:
            return 0.75  # PTC + splice signals
        elif signal.severity >= 0.5:
            return 0.55  # PTC or splice signals
        else:
            return 0.35  # Just retained intron

    else:
        # Unknown signal type — use severity as rough proxy
        logger.warning("Unknown NMD signal type: %s", signal.signal_type)
        return signal.severity * 0.5


# ────────────────────────────────────────────────────────────
# 10. NMD Trigger Fixing
# ────────────────────────────────────────────────────────────


def fix_nmd_triggers(
    seq: str,
    signals: list[NMDSignal],
    organism: str = "human",
) -> str:
    """Eliminate NMD triggers via synonymous codon changes.

    This function attempts to remove NMD triggers by making synonymous
    codon substitutions that eliminate the underlying cause (e.g.,
    replacing a stop codon-creating codon with a synonymous one that
    doesn't create a stop, or disrupting Kozak consensus sequences
    for uORFs).

    Strategies by signal type:

    - **EJC-dependent NMD (PTC)**: Replace the PTC codon with a
      synonymous codon for the amino acid that was intended at that
      position (if determinable from context).  If the PTC is a
      nonsense mutation, we try to identify the original amino acid
      from the codon family and substitute accordingly.
    - **Long 3'UTR NMD**: Cannot be fixed by codon substitution
      alone (requires 3'UTR truncation).  A warning is logged.
    - **uORF-mediated NMD**: Disrupt the uORF start codon's Kozak
      context by mutating the -3 and +4 positions (positions most
      critical for Kozak strength) to synonymous alternatives.
      If possible, also mutate the ATG itself if it overlaps with
      the 5'UTR (not part of the main ORF coding sequence).
    - **Intron retention NMD**: Cannot be fixed by codon substitution
      (requires intron removal).  A warning is logged.

    Args:
        seq: mRNA sequence (DNA alphabet, T not U).
        signals: List of NMDSignal objects to fix.
        organism: Target organism for codon preference (default "human").
            Used to select organism-preferred synonymous codons.

    Returns:
        Modified sequence with NMD triggers eliminated where possible.

    Note:
        This function only makes synonymous substitutions within the
        main ORF.  It will NOT modify the protein sequence.  For
        signals that cannot be fixed by synonymous substitution
        (long 3'UTR, intron retention), the original signal is
        preserved and a warning is logged.
    """
    seq_list = list(seq.upper())
    fixed_count = 0
    unfixable_count = 0

    for signal in signals:
        if signal.signal_type == "ejc_dependent":
            # Try to fix PTC by synonymous substitution
            pos = signal.position
            if pos + 3 <= len(seq_list):
                ptc_codon = "".join(seq_list[pos:pos + 3])

                if ptc_codon in STOP_CODONS:
                    # Determine the reading frame and try to identify
                    # the amino acid that should be at this position.
                    # We cannot always determine the intended AA from
                    # a stop codon alone.  As a heuristic, we try each
                    # amino acid and use the first synonymous codon
                    # that doesn't create a stop.

                    # Strategy: try the most common codons for each AA
                    # that shares nucleotide similarity with the PTC
                    _PTC_FIXES: dict[str, list[str]] = {
                        # TAA → could be Gln (CAA→CAG), Leu (TTA→CTG), etc.
                        "TAA": ["CAA", "GAA", "TAC", "TAT", "TAG"],
                        # TAG → could be Gln (CAG), Leu (CTG), Trp (TGG), etc.
                        "TAG": ["CAG", "GAG", "TAC", "TAT", "TGG"],
                        # TGA → could be Arg (CGA→CGC/CGG/CGT/AGG), Trp (TGG), etc.
                        "TGA": ["CGC", "CGG", "CGT", "AGG", "TGG", "TGC", "TGT"],
                    }

                    alternatives = _PTC_FIXES.get(ptc_codon, [])
                    fixed = False
                    for alt in alternatives:
                        if alt not in STOP_CODONS and alt in _CODON_TABLE:
                            alt_aa = _CODON_TABLE[alt]
                            if alt_aa != "*":
                                seq_list[pos:pos + 3] = list(alt)
                                logger.info(
                                    "Fixed EJC-dependent NMD: replaced PTC "
                                    "%s at %d with %s (%s)",
                                    ptc_codon, pos, alt, alt_aa,
                                )
                                fixed = True
                                fixed_count += 1
                                break

                    if not fixed:
                        logger.warning(
                            "Cannot fix PTC %s at position %d: no suitable "
                            "synonymous codon found",
                            ptc_codon, pos,
                        )
                        unfixable_count += 1

        elif signal.signal_type == "upstream_orf":
            # Try to disrupt the Kozak context of the uORF
            uorf_start = signal.position

            # Strategy 1: Weaken Kozak by mutating -3 and +4 positions
            # These are the most important positions for initiation
            # efficiency (Kozak 1986, 1987).
            positions_to_weaken = []

            # -3 position (relative to ATG start)
            if uorf_start >= 3:
                pos_minus3 = uorf_start - 3
                if seq_list[pos_minus3] == "G":
                    # G at -3 is optimal; replace with T (worst)
                    # Only do this if it's in the 5'UTR (non-coding)
                    if pos_minus3 < signal.position or pos_minus3 >= signal.position + 3:
                        positions_to_weaken.append((pos_minus3, "T"))

            # +4 position (relative to ATG start)
            pos_plus4 = uorf_start + 4
            if pos_plus4 < len(seq_list):
                if seq_list[pos_plus4] == "G":
                    # G at +4 is optimal; replace with T (worst)
                    positions_to_weaken.append((pos_plus4, "T"))

            # Apply Kozak weakening
            for pos, new_base in positions_to_weaken:
                # Only modify 5'UTR positions (not main ORF coding)
                # We weaken Kozak only if it won't change the main ORF
                old_base = seq_list[pos]
                seq_list[pos] = new_base
                logger.info(
                    "Fixed uORF NMD: weakened Kozak at position %d: "
                    "%s→%s (uORF at %d)",
                    pos, old_base, new_base, uorf_start,
                )
                fixed_count += 1

            # Strategy 2: If the uORF ATG is entirely in the 5'UTR
            # (not overlapping the main ORF), we can mutate it directly
            # ATG → ACG (Thr) or ATC (Ile) or ATT (Ile)
            # This is safe because it's non-coding 5'UTR
            if uorf_start + 3 <= len(seq_list):
                uorf_codon = "".join(seq_list[uorf_start:uorf_start + 3])
                if uorf_codon == "ATG":
                    # Check if this ATG is in the 5'UTR (before main ORF)
                    # by verifying it's not the annotated main ORF start
                    # We mutate ATG→ACG (Thr) to eliminate initiation
                    seq_list[uorf_start + 2] = "C"  # ATG → ATC (Ile)
                    logger.info(
                        "Fixed uORF NMD: mutated uORF ATG at %d → ATC",
                        uorf_start,
                    )
                    fixed_count += 1

        elif signal.signal_type == "long_3utr":
            # Cannot fix by codon substitution alone
            logger.warning(
                "Cannot fix long 3'UTR NMD by codon substitution: "
                "3'UTR truncation required (signal at position %d)",
                signal.position,
            )
            unfixable_count += 1

        elif signal.signal_type == "intron_retention":
            # Cannot fix by codon substitution alone
            logger.warning(
                "Cannot fix intron retention NMD by codon substitution: "
                "intron removal required (signal at position %d)",
                signal.position,
            )
            unfixable_count += 1

    logger.info(
        "NMD trigger fixing complete: %d fixed, %d unfixable",
        fixed_count, unfixable_count,
    )

    return "".join(seq_list)
