"""
BioCompiler UTR (Untranslated Region) Model
============================================

Organism-specific 5' and 3' UTR evaluation, scoring, and generation.

Each organism has distinct UTR architecture requirements that critically
affect translation efficiency and mRNA stability.  This module provides
a unified interface for scoring existing UTR sequences and suggesting
optimal ones, with organism-specific logic that accounts for the different
regulatory motifs and structural requirements across prokaryotes, yeast,
and mammals.

The scoring system is motif-based and weighted by biological importance:
prokaryotic 5' UTRs are scored primarily on Shine-Dalgarno presence and
spacing to the start codon, mammalian 5' UTRs on Kozak consensus match,
and 3' UTRs on polyadenylation signals and stability/instability motif
balance.  Suggestion functions return canonical optimal sequences that can
be directly prepended/appended to the coding sequence.

Five organisms are currently supported: *E. coli* (prokaryotic),
*S. cerevisiae* (yeast), *H. sapiens*, *M. musculus*, and CHO-K1
(mammalian).  Alias mappings (e.g., ``"human"`` → ``"Homo_sapiens"``)
are provided for convenience.

Usage::

    from biocompiler.utr_models import score_5utr, score_3utr, suggest_5utr, suggest_3utr

    # Score existing UTR sequences
    s5 = score_5utr("GCCACCATGG", "Homo_sapiens")   # Kozak consensus → ~1.0
    s3 = score_3utr("AATAAAAATAAA", "Homo_sapiens")  # PolyA signals → ~0.9+
    print(f"5' UTR score: {s5:.3f}, 3' UTR score: {s3:.3f}")

    # Get optimal UTR sequences for a given organism
    optimal_5utr = suggest_5utr("Escherichia_coli")  # "TAAGGAGGTAAAAAAAATG"
    optimal_3utr = suggest_3utr("Escherichia_coli")  # "GCGCCGCTTTTTTT"
    print(f"Suggested 5' UTR: {optimal_5utr}")
    print(f"Suggested 3' UTR: {optimal_3utr}")

Public API
----------
- :data:`ORGANISM_UTR_CONFIGS` — built-in registry of organism UTR parameters
- :func:`score_5utr` — score a 5' UTR sequence for a given organism (0–1)
- :func:`score_3utr` — score a 3' UTR sequence for a given organism (0–1)
- :func:`suggest_5utr` — generate an optimal 5' UTR for a given organism
- :func:`suggest_3utr` — generate an optimal 3' UTR for a given organism

References:
  Kozak, M. (1987). "An analysis of 5'-noncoding sequences from 699 vertebrate
  messenger RNAs." *Nucleic Acids Research* 15:8125–8148.
  Shine, J. & Dalgarno, L. (1974). "The 3'-terminal sequence of Escherichia coli
  16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding
  sites." *PNAS* 71:1342–1346.
  Salis, H.M., Mirsky, E.A., & Voigt, C.A. (2009). "Automated design of synthetic
  ribosome binding sites to control protein expression." *Nature Biotechnology*
  27:946–950.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

__all__: list[str] = [
    "UTRConfig",
    "ORGANISM_UTR_CONFIGS",
    "AVAILABLE_ORGANISMS",
    "score_5utr",
    "score_3utr",
    "suggest_5utr",
    "suggest_3utr",
]

logger = logging.getLogger(__name__)

# ─── Data model ────────────────────────────────────────────────────


@dataclass(frozen=True)
class UTRConfig:
    """Organism-specific UTR architecture parameters.

    Attributes:
        organism: Organism identifier (e.g. ``"Escherichia_coli"``).
        utr5_consensus: Consensus 5' UTR sequence for this organism
            (used as a reference; scoring is motif-based).
        utr3_consensus: Consensus 3' UTR sequence for this organism.
        kozak_sequence: Kozak consensus for eukaryotic translation
            initiation (e.g. ``"GCCACCATGG"``). ``None`` for prokaryotes.
        shine_dalgarno: Shine-Dalgarno ribosome binding site sequence
            (e.g. ``"AGGAGG"``). ``None`` for eukaryotes.
        polya_signal: Polyadenylation signal motif
            (e.g. ``"AATAAA"``). ``None`` for prokaryotes.
        stability_motifs: Sequence motifs that promote mRNA stability.
        instability_motifs: Sequence motifs that promote mRNA decay
            (e.g. AU-rich elements).
        splicing_signals: Whether splice-site awareness is needed
            when evaluating UTR sequences.
    """

    organism: str
    utr5_consensus: str
    utr3_consensus: str
    kozak_sequence: str | None
    shine_dalgarno: str | None
    polya_signal: str | None
    stability_motifs: list[str] = field(default_factory=list)
    instability_motifs: list[str] = field(default_factory=list)
    splicing_signals: bool = False


# ─── Built-in registry ────────────────────────────────────────────

ORGANISM_UTR_CONFIGS: dict[str, UTRConfig] = {
    # ── Escherichia coli ─────────────────────────────────────────
    "Escherichia_coli": UTRConfig(
        organism="Escherichia_coli",
        utr5_consensus="TAAGGAGGT",
        utr3_consensus="GCGCCGCTTTTTT",
        kozak_sequence=None,
        shine_dalgarno="AGGAGG",
        polya_signal=None,
        stability_motifs=[
            "GGAGG",   # Strong SD-like stabilising element
            "CGGCGG",  # Stem-forming GC-rich motif
        ],
        instability_motifs=[
            "ATTTA",    # AU-rich instability element (rare in E. coli but relevant)
            "RNaseE",   # Placeholder — single-stranded AU-rich regions
        ],
        splicing_signals=False,
    ),

    # ── Saccharomyces cerevisiae ─────────────────────────────────
    "Saccharomyces_cerevisiae": UTRConfig(
        organism="Saccharomyces_cerevisiae",
        utr5_consensus="TCATATAA",
        utr3_consensus="TAATAAATAA",
        kozak_sequence=None,  # Yeast uses a different mechanism (no Kozak)
        shine_dalgarno=None,
        polya_signal="AATAAA",  # PolyA-like signal
        stability_motifs=[
            "ATATTC",   # Positioning element for 3' end processing
            "AATAAA",   # PolyA signal also acts as stability element
        ],
        instability_motifs=[
            "ATTTA",     # AU-rich element
            "TTATTTAT",  # Extended ARE
        ],
        splicing_signals=True,  # Yeast has introns in ~5% of genes
    ),

    # ── Homo sapiens ─────────────────────────────────────────────
    "Homo_sapiens": UTRConfig(
        organism="Homo_sapiens",
        utr5_consensus="GCCACC",
        utr3_consensus="AATAAAAATAAA",
        kozak_sequence="GCCACCATGG",
        shine_dalgarno=None,
        polya_signal="AATAAA",
        stability_motifs=[
            "CUGB",       # C-rich stability element (simplified motif)
            "CCUCC",      # C-rich element
            "AATAAA",     # PolyA signal contributes to stability
        ],
        instability_motifs=[
            "ATTTA",       # Core AU-rich element (ARE)
            "TTATTTAT",    # Extended ARE (Class II)
            "TATTTAT",     # ARE variant
            "AUUUA",       # ARE pentamer (RNA form — scored as ATTTA on DNA)
        ],
        splicing_signals=True,
    ),

    # ── Mus musculus ─────────────────────────────────────────────
    "Mus_musculus": UTRConfig(
        organism="Mus_musculus",
        utr5_consensus="GCCACC",
        utr3_consensus="AATAAAAATAAA",
        kozak_sequence="GCCACCATGG",
        shine_dalgarno=None,
        polya_signal="AATAAA",
        stability_motifs=[
            "CUGB",
            "CCUCC",
            "AATAAA",
        ],
        instability_motifs=[
            "ATTTA",
            "TTATTTAT",
            "TATTTAT",
        ],
        splicing_signals=True,
    ),

    # ── CHO-K1 (Cricetulus griseus) ─────────────────────────────
    "CHO_K1": UTRConfig(
        organism="CHO_K1",
        utr5_consensus="GCCACC",
        utr3_consensus="AATAAAAATAAA",
        kozak_sequence="GCCACCATGG",
        shine_dalgarno=None,
        polya_signal="AATAAA",
        stability_motifs=[
            "CUGB",
            "CCUCC",
            "AATAAA",
        ],
        instability_motifs=[
            "ATTTA",
            "TTATTTAT",
            "TATTTAT",
        ],
        splicing_signals=True,
    ),
}

# ─── Alias mapping for flexible organism lookup ───────────────────

AVAILABLE_ORGANISMS: list[str] = list(ORGANISM_UTR_CONFIGS.keys())

_UTR_ALIASES: dict[str, str] = {
    "E_coli": "Escherichia_coli",
    "E_coli_K12": "Escherichia_coli",
    "ecoli": "Escherichia_coli",
    "human": "Homo_sapiens",
    "mouse": "Mus_musculus",
    "cho": "CHO_K1",
    "yeast": "Saccharomyces_cerevisiae",
}


def _resolve_organism(organism: str) -> str:
    """Resolve organism name to its canonical key in ORGANISM_UTR_CONFIGS.

    Args:
        organism: Organism identifier (canonical name or alias).

    Returns:
        Canonical organism key.

    Raises:
        ValueError: If the organism is not found.
    """
    if organism in ORGANISM_UTR_CONFIGS:
        return organism
    canonical = _UTR_ALIASES.get(organism)
    if canonical and canonical in ORGANISM_UTR_CONFIGS:
        return canonical
    available = list(ORGANISM_UTR_CONFIGS.keys())
    raise ValueError(
        f"Unknown organism {organism!r}. "
        f"Available: {available}"
    )


def _get_config(organism: str) -> UTRConfig:
    """Look up UTRConfig with alias resolution.

    Args:
        organism: Organism identifier.

    Returns:
        The matching UTRConfig.

    Raises:
        ValueError: If the organism is not found.
    """
    return ORGANISM_UTR_CONFIGS[_resolve_organism(organism)]


# ─── Scoring helpers ─────────────────────────────────────────────

# Shine-Dalgarno spacing: optimal is 5–9 bp from SD to start codon.
# Reference: Shine & Dalgarno 1974; Salis et al. 2009.
_SD_OPTIMAL_MIN: int = 5
_SD_OPTIMAL_MAX: int = 9
_SD_MAX_SPACING: int = 15  # Beyond this, score decays to ~0

# Kozak position weights: offsets relative to the A of ATG (start_pos + offset).
# Strong Kozak: GCCACCATGG — the -3 (A) and +3/+4 positions are most critical.
# Reference: Kozak 1987, 1991.
# In Kozak numbering: -3 is 3bp before A of ATG; +4 is the base after G of ATG.
# In 0-indexed sequence terms, +4 Kozak = offset +3 from start_pos.
_KOZAK_POSITIONS: dict[int, dict[str, float]] = {
    -6: {"G": 0.8, "A": 0.1, "C": 0.05, "T": 0.05},
    -5: {"C": 0.6, "G": 0.2, "A": 0.1, "T": 0.1},
    -4: {"C": 0.5, "G": 0.2, "A": 0.15, "T": 0.15},
    -3: {"A": 1.0, "G": 0.5, "C": 0.2, "T": 0.1},  # Most critical upstream position
    -2: {"C": 0.7, "A": 0.1, "G": 0.1, "T": 0.1},
    -1: {"C": 0.8, "A": 0.05, "G": 0.05, "T": 0.1},
    # offset 0 = A of ATG (always A, not scored)
    +3: {"G": 1.0, "A": 0.3, "C": 0.1, "T": 0.1},  # +4 in Kozak notation (2nd most critical)
}

# Rho-independent terminator parameters
_RHO_INDEP_MIN_GC_HAIRPIN: int = 4   # Minimum GC pairs in hairpin stem
_RHO_INDEP_MIN_U_RUN: int = 4        # Minimum U-run length after hairpin
_RHO_INDEP_HAIRPIN_MAX_DIST: int = 10  # Max distance from hairpin to U-run

# PolyA signal search window from stop codon (bases)
_POLYA_SEARCH_WINDOW: int = 50

# TATA box consensus pattern
_TATA_PATTERN: re.Pattern[str] = re.compile(r"TATA[AT]A[AT]")


# ─── 5' UTR scoring ──────────────────────────────────────────────


def _score_sd_presence(sequence: str, sd_seq: str) -> float:
    """Score the presence and strength of a Shine-Dalgarno sequence.

    Searches for the SD sequence (or close variants allowing 1-2 mismatches)
    in the 5' UTR. Returns a score from 0.0 (absent) to 1.0 (perfect match).

    Args:
        sequence: 5' UTR sequence (uppercased).
        sd_seq: Shine-Dalgarno consensus sequence.

    Returns:
        Score 0.0–1.0 reflecting SD match quality.
    """
    if len(sequence) < len(sd_seq):
        return 0.0

    best_score = 0.0
    sd_len = len(sd_seq)

    for i in range(len(sequence) - sd_len + 1):
        window = sequence[i : i + sd_len]
        matches = sum(1 for a, b in zip(window, sd_seq) if a == b)
        score = matches / sd_len
        if score > best_score:
            best_score = score

    return best_score


def _score_sd_spacing(sequence: str, sd_seq: str, gene_start: str = "ATG") -> float:
    """Score the spacing between Shine-Dalgarno and start codon.

    Optimal spacing is 5–9 bp. Scores decay for suboptimal spacing.

    Args:
        sequence: 5' UTR + start codon sequence.
        sd_seq: Shine-Dalgarno consensus.
        gene_start: Start codon sequence (default "ATG").

    Returns:
        Score 0.0–1.0 for spacing quality.
    """
    # Find the best SD match position
    sd_len = len(sd_seq)
    best_pos = -1
    best_score = 0.0

    upper_seq = sequence.upper()
    for i in range(len(upper_seq) - sd_len + 1):
        window = upper_seq[i : i + sd_len]
        matches = sum(1 for a, b in zip(window, sd_seq) if a == b)
        score = matches / sd_len
        if score > best_score:
            best_score = score
            best_pos = i

    if best_pos < 0 or best_score < 0.5:
        return 0.0

    # Find the start codon position
    start_pos = upper_seq.find(gene_start.upper())
    if start_pos < 0:
        return 0.0

    # Spacing = distance from end of SD to start of ATG
    spacing = start_pos - (best_pos + sd_len)
    if spacing < 0:
        return 0.0

    if _SD_OPTIMAL_MIN <= spacing <= _SD_OPTIMAL_MAX:
        return 1.0
    elif spacing < _SD_OPTIMAL_MIN:
        # Too close — linear penalty
        return spacing / _SD_OPTIMAL_MIN
    else:
        # Too far — exponential decay
        excess = spacing - _SD_OPTIMAL_MAX
        return max(0.0, (0.5 ** (excess / 3.0)))


def _score_kozak(sequence: str, gene_start: str = "ATG") -> float:
    """Score Kozak consensus match around the start codon.

    Evaluates positions -6 to +3 (offsets from A of ATG), with emphasis
    on the critical -3 (A) and +3/+4 Kozak (G) positions.

    Args:
        sequence: 5' UTR sequence including the start codon.
        gene_start: Start codon (default "ATG").

    Returns:
        Score 0.0–1.0 for Kozak match quality.
    """
    upper_seq = sequence.upper()
    start_pos = upper_seq.find(gene_start.upper())
    if start_pos < 0:
        return 0.0

    total_weight = 0.0
    max_weight = 0.0

    for offset, base_scores in _KOZAK_POSITIONS.items():
        seq_pos = start_pos + offset
        max_weight += max(base_scores.values())
        if 0 <= seq_pos < len(upper_seq):
            base = upper_seq[seq_pos]
            total_weight += base_scores.get(base, 0.0)
        # If position is outside the sequence, contribute 0

    if max_weight == 0:
        return 0.0

    return total_weight / max_weight


def _score_tata_box(sequence: str) -> float:
    """Score TATA box presence in the 5' UTR region.

    Args:
        sequence: 5' UTR sequence.

    Returns:
        Score 0.0–1.0 for TATA box presence (1.0 if found).
    """
    if _TATA_PATTERN.search(sequence.upper()):
        return 1.0
    # Partial match — look for TATA without the full consensus
    if "TATA" in sequence.upper():
        return 0.5
    return 0.0


def _score_yeast_initiator(sequence: str, gene_start: str = "ATG") -> float:
    """Score yeast-specific translation initiation context.

    Yeast lacks the mammalian Kozak consensus but prefers:
    - A at position -3 relative to ATG (critical in S. cerevisiae)
    - A-rich region immediately upstream of the start codon
    - Low secondary structure potential in the leader

    Args:
        sequence: 5' UTR sequence including the start codon.
        gene_start: Start codon (default "ATG").

    Returns:
        Score 0.0–1.0 for yeast initiator context quality.
    """
    upper_seq = sequence.upper()
    start_pos = upper_seq.find(gene_start.upper())
    if start_pos < 0:
        return 0.0

    score = 0.0
    max_score = 0.0

    # Position -3: A is strongly preferred (most important)
    max_score += 1.0
    if start_pos >= 3 and upper_seq[start_pos - 3] == "A":
        score += 1.0
    elif start_pos >= 3 and upper_seq[start_pos - 3] in "TG":
        score += 0.3

    # Position -1: A or C is preferred
    max_score += 0.5
    if start_pos >= 1 and upper_seq[start_pos - 1] in "AC":
        score += 0.5
    elif start_pos >= 1 and upper_seq[start_pos - 1] in "TG":
        score += 0.2

    # Position -2: any base, slight preference for A
    max_score += 0.2
    if start_pos >= 2 and upper_seq[start_pos - 2] == "A":
        score += 0.2

    # A-rich context in the 10 bp upstream of ATG
    max_score += 0.3
    if start_pos >= 10:
        upstream_10 = upper_seq[start_pos - 10 : start_pos]
        a_fraction = upstream_10.count("A") / len(upstream_10)
        score += 0.3 * min(1.0, a_fraction / 0.4)  # 40% A content = perfect
    elif start_pos > 0:
        upstream = upper_seq[:start_pos]
        a_fraction = upstream.count("A") / len(upstream) if upstream else 0
        score += 0.3 * min(1.0, a_fraction / 0.4)

    if max_score == 0:
        return 0.0
    return score / max_score


def score_5utr(sequence: str, organism: str) -> float:
    """Score a 5' UTR sequence for a given organism.

    Scoring logic by organism type:

    - **E. coli** (prokaryotic): Evaluates Shine-Dalgarno presence (40%)
      and spacing to start codon (60%). Optimal SD-ATG spacing is 5–9 bp.
    - **Yeast**: Evaluates TATA box presence (30%), yeast-specific
      initiator context (40%), and 5' UTR length (30%).
    - **Mammals** (human, mouse, CHO): Evaluates Kozak consensus match
      (70%) and TATA box presence (30%).

    The input sequence should include the start codon (ATG) at the 3' end
    for accurate spacing and Kozak scoring.

    Args:
        sequence: 5' UTR sequence (optionally including the start codon).
        organism: Organism identifier (e.g. ``"Homo_sapiens"`` or ``"human"``).

    Returns:
        Score between 0.0 and 1.0, where 1.0 is optimal.

    Examples::

        >>> score_5utr("TAAGGAGGTNNNNNATG", "E_coli")  # doctest: +SKIP
        0.9...
        >>> score_5utr("GCCACCATGG", "Homo_sapiens")  # doctest: +SKIP
        1.0
    """
    config = _get_config(organism)
    seq = sequence.upper()

    if not seq:
        return 0.0

    # ── Prokaryotic (E. coli): Shine-Dalgarno + spacing ──────────
    if config.shine_dalgarno is not None:
        sd_score = _score_sd_presence(seq, config.shine_dalgarno)
        spacing_score = _score_sd_spacing(seq, config.shine_dalgarno)
        score = 0.4 * sd_score + 0.6 * spacing_score
        logger.debug(
            "5' UTR score for %s: SD=%.3f, spacing=%.3f, total=%.3f",
            organism, sd_score, spacing_score, score,
        )
        return max(0.0, min(1.0, score))

    # ── Yeast: TATA box + initiator context ──────────────────────
    if config.organism == "Saccharomyces_cerevisiae":
        tata_score = _score_tata_box(seq)
        # Yeast uses a different initiation mechanism from the Kozak consensus.
        # Key features: A-rich context around ATG (especially A at -3),
        # and short leader length (typically 20–80 nt).
        initiator_score = _score_yeast_initiator(seq)
        # Yeast prefers 20–80 bp 5' UTRs; very short leaders (<10) are poor
        length_score = 1.0 if 20 <= len(seq) <= 100 else max(0.0, 0.5 - abs(len(seq) - 60) / 200)
        score = 0.3 * tata_score + 0.4 * initiator_score + 0.3 * length_score
        logger.debug(
            "5' UTR score for %s: TATA=%.3f, initiator=%.3f, length=%.3f, total=%.3f",
            organism, tata_score, initiator_score, length_score, score,
        )
        return max(0.0, min(1.0, score))

    # ── Mammalian (human, mouse, CHO): Kozak + TATA ─────────────
    kozak_score = _score_kozak(seq)
    tata_score = _score_tata_box(seq)
    score = 0.7 * kozak_score + 0.3 * tata_score
    logger.debug(
        "5' UTR score for %s: Kozak=%.3f, TATA=%.3f, total=%.3f",
        organism, kozak_score, tata_score, score,
    )
    return max(0.0, min(1.0, score))


# ─── 3' UTR scoring ──────────────────────────────────────────────


def _score_rho_independent_terminator(sequence: str) -> float:
    """Score a Rho-independent terminator in a 3' UTR.

    A Rho-independent terminator consists of:
    1. A GC-rich hairpin-forming region (≥4 GC pairs)
    2. A poly-U tract (≥4 U's) within 10 bp downstream

    Args:
        sequence: 3' UTR sequence (DNA, T's not U's).

    Returns:
        Score 0.0–1.0 for terminator quality.
    """
    seq = sequence.upper()
    best_score = 0.0

    # Scan for GC-rich regions that could form hairpin stems
    for i in range(len(seq)):
        # Look for a GC-rich region of at least 6 bases
        gc_run = 0
        j = i
        while j < len(seq) and seq[j] in "GC":
            gc_run += 1
            j += 1

        if gc_run >= _RHO_INDEP_MIN_GC_HAIRPIN:
            # Check for poly-T (U in RNA) tract downstream
            downstream_start = j
            downstream_end = min(len(seq), downstream_start + _RHO_INDEP_HAIRPIN_MAX_DIST)

            # Find the longest T-run in the downstream window
            t_run = 0
            max_t_run = 0
            for k in range(downstream_start, downstream_end):
                if seq[k] == "T":
                    t_run += 1
                    max_t_run = max(max_t_run, t_run)
                else:
                    t_run = 0

            if max_t_run >= _RHO_INDEP_MIN_U_RUN:
                # Score based on GC hairpin length and T-run length
                gc_score = min(1.0, gc_run / 8.0)  # 8+ GC bases = perfect
                t_score = min(1.0, max_t_run / 7.0)  # 7+ T's = perfect
                score = 0.5 * gc_score + 0.5 * t_score
                if score > best_score:
                    best_score = score

    return best_score


def _score_polya_signal(sequence: str, polya_signal: str) -> float:
    """Score polyadenylation signal presence and positioning in a 3' UTR.

    The canonical polyA signal (AATAAA) should appear within 10–30 bp
    downstream of the stop codon for optimal 3' end processing.

    Args:
        sequence: 3' UTR sequence.
        polya_signal: PolyA signal motif (e.g. ``"AATAAA"``).

    Returns:
        Score 0.0–1.0 for polyA signal quality.
    """
    seq = sequence.upper()
    signal_upper = polya_signal.upper()

    pos = seq.find(signal_upper)
    if pos < 0:
        # Check for variant polyA signals
        variants = ["ATTAAA", "AGTAAA", "TATAAA", "AATACA", "AATATA"]
        for variant in variants:
            pos = seq.find(variant)
            if pos >= 0:
                # Variant signals score lower
                return 0.6
        return 0.0

    # Score based on position — optimal is 10–30 bp from start of UTR
    # (i.e., from stop codon)
    if 10 <= pos <= 30:
        return 1.0
    elif pos < 10:
        return 0.7  # Too close — still functional but suboptimal
    elif pos <= _POLYA_SEARCH_WINDOW:
        # Decay with distance
        return max(0.3, 1.0 - (pos - 30) / 40.0)
    else:
        return 0.2  # Present but far


def _score_stability_motifs(sequence: str, config: UTRConfig) -> float:
    """Score stability and instability motif balance in a UTR.

    Args:
        sequence: UTR sequence.
        config: Organism UTR configuration with motif lists.

    Returns:
        Score 0.0–1.0 where higher = more stable mRNA expected.
    """
    seq = sequence.upper()

    stability_count = 0
    for motif in config.stability_motifs:
        stability_count += seq.count(motif.upper())

    instability_count = 0
    for motif in config.instability_motifs:
        instability_count += seq.count(motif.upper())

    # Base score from stability motifs (capped at 1.0)
    stability_score = min(1.0, stability_count * 0.3) if stability_count > 0 else 0.3

    # Penalty from instability motifs
    instability_penalty = min(1.0, instability_count * 0.25)

    return max(0.0, stability_score - instability_penalty * 0.5)


def score_3utr(sequence: str, organism: str) -> float:
    """Score a 3' UTR sequence for a given organism.

    Scoring logic by organism type:

    - **E. coli** (prokaryotic): Evaluates Rho-independent terminator
      quality (GC-rich hairpin + poly-U tract, 70%) and stability motif
      balance (30%).
    - **Yeast**: Evaluates polyA-like signal presence (50%) and
      stability motif balance (50%).
    - **Mammals** (human, mouse, CHO): Evaluates polyA signal presence
      and positioning (50%), and stability/instability motif balance (50%),
      including AU-rich element detection.

    Args:
        sequence: 3' UTR sequence.
        organism: Organism identifier (e.g. ``"Homo_sapiens"`` or ``"human"``).

    Returns:
        Score between 0.0 and 1.0, where 1.0 is optimal.

    Examples::

        >>> score_3utr("GCGCCGCTTTTTT", "E_coli")  # doctest: +SKIP
        0.8...
        >>> score_3utr("AATAAAAATAAA", "Homo_sapiens")  # doctest: +SKIP
        0.9...
    """
    config = _get_config(organism)
    seq = sequence.upper()

    if not seq:
        return 0.0

    # ── Prokaryotic (E. coli): Rho-independent terminator ────────
    if config.shine_dalgarno is not None:
        terminator_score = _score_rho_independent_terminator(seq)
        motif_score = _score_stability_motifs(seq, config)
        score = 0.7 * terminator_score + 0.3 * motif_score
        logger.debug(
            "3' UTR score for %s: terminator=%.3f, motifs=%.3f, total=%.3f",
            organism, terminator_score, motif_score, score,
        )
        return max(0.0, min(1.0, score))

    # ── Yeast: PolyA signal + stability ──────────────────────────
    if config.organism == "Saccharomyces_cerevisiae":
        if config.polya_signal:
            polya_score = _score_polya_signal(seq, config.polya_signal)
        else:
            polya_score = 0.0
        motif_score = _score_stability_motifs(seq, config)
        score = 0.5 * polya_score + 0.5 * motif_score
        logger.debug(
            "3' UTR score for %s: polyA=%.3f, motifs=%.3f, total=%.3f",
            organism, polya_score, motif_score, score,
        )
        return max(0.0, min(1.0, score))

    # ── Mammalian (human, mouse, CHO): PolyA + stability ────────
    if config.polya_signal:
        polya_score = _score_polya_signal(seq, config.polya_signal)
    else:
        polya_score = 0.0
    motif_score = _score_stability_motifs(seq, config)
    score = 0.5 * polya_score + 0.5 * motif_score
    logger.debug(
        "3' UTR score for %s: polyA=%.3f, motifs=%.3f, total=%.3f",
        organism, polya_score, motif_score, score,
    )
    return max(0.0, min(1.0, score))


# ─── UTR suggestion ──────────────────────────────────────────────


def suggest_5utr(organism: str, gene_start: str = "ATG") -> str:
    """Generate an optimal 5' UTR for the given organism.

    Constructs a UTR sequence incorporating the organism-specific
    initiation signals:

    - **E. coli**: ``"AGGAGGT"`` + 7 bp optimal A-rich spacing + start codon.
      The spacing region uses A/T-rich sequence to avoid secondary
      structure that could occlude the RBS.
    - **Yeast**: TATA box + A-rich leader (~25 bp) + A at -3 + ATG context.
    - **Mammals** (human, mouse, CHO): Kozak consensus
      ``"GCCACC"`` + ATG (the full Kozak context ``GCCACCATGG``).

    Args:
        organism: Organism identifier.
        gene_start: Start codon (default ``"ATG"``).

    Returns:
        Optimal 5' UTR sequence including the start codon.

    Examples::

        >>> suggest_5utr("E_coli")
        'AGGAGGTAAAAAAAATG'
        >>> suggest_5utr("Homo_sapiens")
        'GCCACCATGG'
    """
    config = _get_config(organism)
    start = gene_start.upper()

    # ── Prokaryotic (E. coli): SD + optimal spacing ──────────────
    if config.shine_dalgarno is not None:
        # Use AGGAGGT (includes the G that pairs with anti-SD on 16S rRNA)
        sd = "AGGAGGT"
        # Optimal spacing = 7 bp (middle of 5-9 range)
        # Use A-rich spacers to minimise secondary structure
        spacer = "AAAAAAA"
        utr = sd + spacer + start
        logger.debug("Suggested 5' UTR for %s: %s", organism, utr)
        return utr

    # ── Yeast: TATA box + leader + initiator ─────────────────────
    if config.organism == "Saccharomyces_cerevisiae":
        # Yeast 5' UTRs are typically 20-80 bp.
        # Include TATA box, A-rich leader (~30 bp), and A at -3 for
        # optimal yeast initiator context.
        leader = "AAAAA" * 5  # 25 bp A-rich leader (low 2ary structure)
        utr = "TATATAA" + leader + "AAA" + start
        logger.debug("Suggested 5' UTR for %s: %s", organism, utr)
        return utr

    # ── Mammalian: Kozak consensus ──────────────────────────────
    if config.kozak_sequence is not None:
        # Return the full Kozak consensus including ATG
        utr = config.kozak_sequence
        logger.debug("Suggested 5' UTR for %s: %s", organism, utr)
        return utr

    # Fallback: just the start codon
    logger.warning("No 5' UTR template for %s; returning start codon only", organism)
    return start


def suggest_3utr(organism: str) -> str:
    """Generate an optimal 3' UTR for the given organism.

    Constructs a UTR sequence incorporating the organism-specific
    termination/processing signals:

    - **E. coli**: GC-rich hairpin (``GCGCCGC``) + poly-U tract
      (``TTTTTTT``) — classic Rho-independent terminator.
    - **Yeast**: PolyA-like signal + short downstream element.
    - **Mammals** (human, mouse, CHO): ~15 bp spacer + polyadenylation
      signal (``AATAAA``) + downstream element + second polyA signal.

    Args:
        organism: Organism identifier.

    Returns:
        Optimal 3' UTR sequence.

    Examples::

        >>> suggest_3utr("E_coli")
        'GCGCCGCTTTTTTT'
        >>> suggest_3utr("Homo_sapiens")
        'AAAAAAAAAAAAAAAATAAAAATAAA'
    """
    config = _get_config(organism)

    # ── Prokaryotic (E. coli): Rho-independent terminator ────────
    if config.shine_dalgarno is not None:
        # GC-rich hairpin-forming region + poly-U tract
        utr = "GCGCCGC" + "TTTTTTT"
        logger.debug("Suggested 3' UTR for %s: %s", organism, utr)
        return utr

    # ── Yeast: PolyA-like signals ───────────────────────────────
    if config.organism == "Saccharomyces_cerevisiae":
        # Yeast positioning element (PE) + efficiency element (EE)
        # + polyA tract.  Add padding so polyA signal lands at ~15 bp.
        utr = "AAAAAAAAAAAAAATAATAAATAA"
        logger.debug("Suggested 3' UTR for %s: %s", organism, utr)
        return utr

    # ── Mammalian: Dual polyA signals ──────────────────────────
    if config.polya_signal is not None:
        # Optimal placement: polyA signal 10-30 bp from stop codon.
        # Add 15 bp spacer (A-rich, low 2ary structure), then AATAAA,
        # then downstream element and a second polyA signal.
        spacer = "A" * 15
        utr = spacer + config.polya_signal + "AATAAA"
        logger.debug("Suggested 3' UTR for %s: %s", organism, utr)
        return utr

    # Fallback: generic poly-T tract
    logger.warning("No 3' UTR template for %s; returning generic terminator", organism)
    return "TTTTTTT"
