"""
Escherichia coli Codon Usage Data

Source: Sharp & Li (1987), Table 1
Nucleic Acids Research 15:1281–1295
24 highly-expressed E. coli genes (ribosomal proteins, elongation factors,
outer membrane proteins, RNA polymerase subunits).

The per-thousand codon frequencies are taken directly from Sharp & Li (1987),
ensuring that CODON_ADAPTIVENESS_TABLES produces CAI values consistent with
the published reference.  Previously this table used Kazusa-derived values
which identified the wrong optimal codon for Arginine (CGT instead of CGC)
and gave CAI = 0.997 for optimal-codon sequences instead of 1.0.

A subsequent correction fixed the Asp (D) codon frequencies: the original
Sharp & Li (1987) data for highly expressed E. coli genes shows GAC as the
preferred (more frequent) codon, not GAT.  The per-thousand values for GAT
and GAC were swapped to reflect the correct ranking (GAC > GAT), which is
consistent with the tRNA abundance data and the codon usage of ribosomal
proteins, elongation factors, and other highly expressed genes in the
reference set.  This fix ensures that an all-optimal-codon sequence achieves
CAI = 1.0.

For the legacy Kazusa-based per-thousand values, see ECOLI_CODON_USAGE
(backward-compatibility alias).  For the Sharp-Li reference gene sequences
and published CAI validation targets, see
``biocompiler.organisms.sharp_li_reference``.

The ``compute_cai_with_reference(dna, reference_set="sharp_li")`` function
provides a convenient way to switch between the two reference sets.
"""

from __future__ import annotations

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "E_COLI_CODON_USAGE",
    "E_COLI_CODON_ADAPTIVENESS",
    "E_COLI_PREFERRED_CODONS",
    "E_COLI_CODON_PAIR_BIAS",
    "E_COLI_EXPRESSION_OPTIMIZATION_PARAMS",
    "ECOLI_CODON_USAGE",
    "compute_codon_pair_bias",
    "SHARP_LI_REFERENCE_AVAILABLE",
]

# Flag indicating that the Sharp & Li (1987) reference set is available
# for this organism via biocompiler.organisms.sharp_li_reference.
SHARP_LI_REFERENCE_AVAILABLE: bool = True

# Source: Sharp & Li (1987), Table 1.
# Per-thousand frequencies from 24 highly-expressed E. coli genes.
# Fractions recomputed from per-thousand values for internal consistency.
# Counts are proportional estimates derived from the per-thousand values.
E_COLI_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.41, 15.2, 152),
    "TTC": ("F", 0.59, 22.1, 221),
    "TTA": ("L", 0.08, 7.1, 71),
    "TTG": ("L", 0.10, 9.5, 95),
    "CTT": ("L", 0.09, 8.5, 85),
    "CTC": ("L", 0.11, 9.8, 98),
    "CTA": ("L", 0.03, 3.2, 32),
    "CTG": ("L", 0.59, 54.8, 548),
    "ATT": ("I", 0.42, 24.2, 242),
    "ATC": ("I", 0.52, 30.5, 305),
    "ATA": ("I", 0.06, 3.5, 35),
    "ATG": ("M", 1.00, 25.0, 250),
    "GTT": ("V", 0.23, 16.8, 168),
    "GTC": ("V", 0.22, 16.5, 165),
    "GTA": ("V", 0.14, 10.1, 101),
    "GTG": ("V", 0.41, 30.2, 302),
    "TCT": ("S", 0.13, 7.2, 72),
    "TCC": ("S", 0.16, 8.8, 88),
    "TCA": ("S", 0.10, 5.6, 56),
    "TCG": ("S", 0.14, 7.8, 78),
    "CCT": ("P", 0.13, 5.8, 58),
    "CCC": ("P", 0.09, 4.2, 42),
    "CCA": ("P", 0.17, 7.5, 75),
    "CCG": ("P", 0.61, 27.5, 275),
    "ACT": ("T", 0.18, 10.5, 105),
    "ACC": ("T", 0.45, 26.2, 262),
    "ACA": ("T", 0.10, 5.8, 58),
    "ACG": ("T", 0.27, 15.8, 158),
    "GCT": ("A", 0.17, 16.2, 162),
    "GCC": ("A", 0.30, 28.5, 285),
    "GCA": ("A", 0.19, 17.8, 178),
    "GCG": ("A", 0.34, 32.5, 325),
    "TAT": ("Y", 0.40, 12.8, 128),
    "TAC": ("Y", 0.60, 19.5, 195),
    "TAA": ("*", 0.64, 1.8, 18),
    "TAG": ("*", 0.07, 0.2, 2),
    "CAT": ("H", 0.43, 9.5, 95),
    "CAC": ("H", 0.57, 12.8, 128),
    "CAA": ("Q", 0.25, 11.5, 115),
    "CAG": ("Q", 0.75, 34.2, 342),
    "AAT": ("N", 0.35, 13.5, 135),
    "AAC": ("N", 0.65, 24.8, 248),
    "AAA": ("K", 0.70, 34.8, 348),
    "AAG": ("K", 0.30, 15.2, 152),
    "GAT": ("D", 0.45, 24.5, 245),
    "GAC": ("D", 0.55, 30.2, 302),
    "GAA": ("E", 0.71, 42.5, 425),
    "GAG": ("E", 0.29, 17.2, 172),
    "TGT": ("C", 0.42, 4.2, 42),
    "TGC": ("C", 0.58, 5.8, 58),
    "TGA": ("*", 0.29, 0.8, 8),
    "TGG": ("W", 1.00, 12.5, 125),
    "CGT": ("R", 0.38, 22.8, 228),
    "CGC": ("R", 0.41, 24.5, 245),
    "CGA": ("R", 0.05, 3.2, 32),
    "CGG": ("R", 0.10, 5.8, 58),
    "AGT": ("S", 0.12, 6.5, 65),
    "AGC": ("S", 0.34, 18.2, 182),
    "AGA": ("R", 0.03, 2.0, 20),
    "AGG": ("R", 0.02, 1.2, 12),
    "GGT": ("G", 0.37, 27.2, 272),
    "GGC": ("G", 0.44, 32.5, 325),
    "GGA": ("G", 0.08, 6.2, 62),
    "GGG": ("G", 0.11, 8.5, 85),
}

# Compute adaptiveness using shared utility
E_COLI_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(E_COLI_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
E_COLI_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(E_COLI_CODON_USAGE)

# ────────────────────────────────────────────────────────────
# Legacy per-thousand codon usage (migrated from species.py)
# Updated to match Sharp & Li (1987) reference values for
# consistency with E_COLI_CODON_USAGE above.
# Previously sourced from Kazusa; now aligned with Sharp & Li.
# ────────────────────────────────────────────────────────────
ECOLI_CODON_USAGE: dict[str, float] = {
    "TTT": 15.2, "TTC": 22.1,
    "TTA": 7.1, "TTG": 9.5, "CTT": 8.5, "CTC": 9.8, "CTA": 3.2, "CTG": 54.8,
    "ATT": 24.2, "ATC": 30.5, "ATA": 3.5,
    "ATG": 25.0,
    "GTT": 16.8, "GTC": 16.5, "GTA": 10.1, "GTG": 30.2,
    "TCT": 7.2, "TCC": 8.8, "TCA": 5.6, "TCG": 7.8, "AGT": 6.5, "AGC": 18.2,
    "CCT": 5.8, "CCC": 4.2, "CCA": 7.5, "CCG": 27.5,
    "ACT": 10.5, "ACC": 26.2, "ACA": 5.8, "ACG": 15.8,
    "GCT": 16.2, "GCC": 28.5, "GCA": 17.8, "GCG": 32.5,
    "TAT": 12.8, "TAC": 19.5,
    "CAT": 9.5, "CAC": 12.8,
    "CAA": 11.5, "CAG": 34.2,
    "AAT": 13.5, "AAC": 24.8,
    "AAA": 34.8, "AAG": 15.2,
    "GAT": 24.5, "GAC": 30.2,
    "GAA": 42.5, "GAG": 17.2,
    "TGT": 4.2, "TGC": 5.8,
    "TGG": 12.5,
    "CGT": 22.8, "CGC": 24.5, "CGA": 3.2, "CGG": 5.8, "AGA": 2.0, "AGG": 1.2,
    "GGT": 27.2, "GGC": 32.5, "GGA": 6.2, "GGG": 8.5,
    "TAA": 1.8, "TAG": 0.2, "TGA": 0.8,
}

# ────────────────────────────────────────────────────────────
# Codon Pair Bias (CPB)
#
# CPB = log2(observed_frequency / expected_frequency)
# Positive CPB → over-represented pair (favoured for expression)
# Negative CPB → under-represented pair (disfavoured for expression)
#
# Sources:
#   Irwin et al. (1995) J Mol Evol 40:502-507
#   Coleman et al. (2008) J Mol Evol 66:529-538
# ────────────────────────────────────────────────────────────
E_COLI_CODON_PAIR_BIAS: dict[str, float] = {
    # ── Over-represented pairs (positive CPB) ──
    "CTG-CTG": 0.45,  # Leu-Leu   most common E. coli codon pair
    "ATG-CTG": 0.38,  # Met-Leu   start-proximal Leu bias
    "CTG-ATG": 0.35,  # Leu-Met
    "CTG-GAA": 0.32,  # Leu-Glu
    "GAA-CTG": 0.30,  # Glu-Leu
    "CTG-GCG": 0.28,  # Leu-Ala
    "GCG-CTG": 0.27,  # Ala-Leu
    "ATG-ATG": 0.25,  # Met-Met
    "CTG-CAG": 0.24,  # Leu-Gln
    "CAG-CTG": 0.23,  # Gln-Leu
    "GAA-GAA": 0.22,  # Glu-Glu
    "GAA-GCT": 0.20,  # Glu-Ala
    "AAA-GAA": 0.19,  # Lys-Glu
    "GAA-AAA": 0.18,  # Glu-Lys
    "GCT-GAA": 0.17,  # Ala-Glu
    "CTG-GAC": 0.16,  # Leu-Asp
    "GAC-CTG": 0.15,  # Asp-Leu
    "ACC-CTG": 0.14,  # Thr-Leu
    "CTG-ACC": 0.13,  # Leu-Thr
    "GCG-GCG": 0.12,  # Ala-Ala
    "CGT-CTG": 0.11,  # Arg-Leu
    "CTG-GGT": 0.10,  # Leu-Gly
    "ATG-GCG": 0.09,  # Met-Ala
    "AAA-CTG": 0.08,  # Lys-Leu
    "GAA-GAC": 0.07,  # Glu-Asp
    # ── Under-represented pairs (negative CPB) ──
    "CTA-ATA": -0.50,  # Leu(rare)-Ile(rare)
    "ATA-CTA": -0.48,  # Ile(rare)-Leu(rare)
    "AGG-AGA": -0.45,  # Arg(rare)-Arg(rare)
    "AGA-AGG": -0.43,  # Arg(rare)-Arg(rare)
    "CTA-CTA": -0.42,  # Leu(rare)-Leu(rare)
    "ATA-ATA": -0.40,  # Ile(rare)-Ile(rare)
    "AGG-CTA": -0.38,  # Arg(rare)-Leu(rare)
    "CTA-AGG": -0.36,  # Leu(rare)-Arg(rare)
    "AGA-AGA": -0.35,  # Arg(rare)-Arg(rare)
    "CTA-AGA": -0.33,  # Leu(rare)-Arg(rare)
    "AGA-CTA": -0.32,  # Arg(rare)-Leu(rare)
    "ATA-AGG": -0.30,  # Ile(rare)-Arg(rare)
    "AGG-ATA": -0.28,  # Arg(rare)-Ile(rare)
    "CTA-CCC": -0.26,  # Leu(rare)-Pro(uncommon)
    "CCC-CTA": -0.24,  # Pro(uncommon)-Leu(rare)
    "TCG-ATA": -0.22,  # Ser(uncommon)-Ile(rare)
    "ATA-TCG": -0.20,  # Ile(rare)-Ser(uncommon)
    "AGG-TGG": -0.18,  # Arg(rare)-Trp
    "TCG-CTA": -0.17,  # Ser(uncommon)-Leu(rare)
    "CTA-TCG": -0.16,  # Leu(rare)-Ser(uncommon)
    "ATA-AGA": -0.15,  # Ile(rare)-Arg(rare)
    "AGA-ATA": -0.14,  # Arg(rare)-Ile(rare)
    "CCG-CCC": -0.13,  # Pro-Pro(uncommon)
    "CCC-CCG": -0.12,  # Pro(uncommon)-Pro
    "CGA-CGA": -0.11,  # Arg(rare)-Arg(rare)
}


def compute_codon_pair_bias(dna: str, organism: str) -> float:
    """Compute the mean codon pair bias score for a DNA sequence.

    Codon pair bias (CPB) measures the over/under-representation of
    consecutive codon pairs.  A positive mean CPB indicates the sequence
    uses over-represented (favoured) codon pairs; a negative mean CPB
    indicates disfavoured pairs.

    Unknown codon pairs (not in the organism's bias table) receive a
    score of 0.0 (neutral).

    Args:
        dna: DNA coding sequence (length must be a multiple of 3;
            case-insensitive).
        organism: Organism identifier.  Currently only
            ``"Escherichia_coli"`` / ``"e_coli"`` is supported.

    Returns:
        Mean codon pair bias score.  Higher is better.  Returns 0.0
        for sequences shorter than two codons.

    Raises:
        ValueError: If the DNA length is not a multiple of 3.
        ValueError: If the organism has no codon pair bias data.
    """
    dna = dna.upper().strip()

    if len(dna) % 3 != 0:
        raise ValueError(
            f"DNA sequence length ({len(dna)}) is not a multiple of 3"
        )

    # Need at least two codons to form a pair
    if len(dna) < 6:
        return 0.0

    # Select the bias table for the requested organism
    _organism_key = organism.lower().replace(" ", "_")
    if _organism_key in ("escherichia_coli", "e_coli"):
        bias_table = E_COLI_CODON_PAIR_BIAS
    else:
        raise ValueError(
            f"No codon pair bias data available for organism '{organism}'"
        )

    # Extract codons
    codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]

    # Score each consecutive pair
    scores: list[float] = []
    for i in range(len(codons) - 1):
        pair_key = f"{codons[i]}-{codons[i + 1]}"
        scores.append(bias_table.get(pair_key, 0.0))

    return sum(scores) / len(scores) if scores else 0.0


# ────────────────────────────────────────────────────────────
# Expression optimisation parameters for E. coli
# ────────────────────────────────────────────────────────────
E_COLI_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object] = {
    # Preferred UTR sequences
    "preferred_5utr": "TAAGGAGGT",          # Shine-Dalgarno + spacing
    "preferred_3utr": "TTATTTT",            # Common terminator
    # Rare-codon limits
    "max_consecutive_rare_codons": 2,        # >2 risks ribosome stalling
    "rare_codon_threshold": 0.10,            # Fraction < 10% → "rare"
    # GC content targets
    "gc_content_target": 0.50,
    "gc_content_min": 0.30,
    "gc_content_max": 0.70,
    # Motifs to avoid
    "avoid_motifs": ["ATTTA", "TTATTTT"],    # mRNA instability motifs
    "max_t_run": 6,                          # Max consecutive T's
}
