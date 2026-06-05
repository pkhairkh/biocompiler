"""
Mouse (Mus musculus) Codon Usage Data

Source: Kazusa Codon Usage Database
59,432 CDSs, 24,040,270 codons
Coding GC: 52.25%

Verified: all 61 sense codons present, preferred codons match
published mouse codon usage tables.  AGA is the preferred Arg
codon (consistent with mammalian pattern).
"""

from ._utils import CodonUsageTable, compute_codon_adaptiveness, compute_preferred_codons

__all__ = [
    "MOUSE_CODON_USAGE",
    "MOUSE_CODON_ADAPTIVENESS",
    "MOUSE_PREFERRED_CODONS",
    "MOUSE_CODON_PAIR_BIAS",
    "MOUSE_EXPRESSION_OPTIMIZATION_PARAMS",
    "MOUSE_UTR_MODELS",
]

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
MOUSE_CODON_USAGE: CodonUsageTable = {
    "TTT": ("F", 0.46, 17.5, 420837),
    "TTC": ("F", 0.54, 20.6, 494997),
    "TTA": ("L", 0.07, 7.6, 183106),
    "TTG": ("L", 0.13, 13.1, 315351),
    "CTT": ("L", 0.13, 12.6, 302711),
    "CTC": ("L", 0.22, 21.0, 504980),
    "CTA": ("L", 0.07, 7.2, 173960),
    "CTG": ("L", 0.38, 39.6, 952185),
    "ATT": ("I", 0.36, 16.0, 384673),
    "ATC": ("I", 0.48, 21.4, 514619),
    "ATA": ("I", 0.16, 7.3, 174852),
    "ATG": ("M", 1.00, 22.3, 535014),
    "GTT": ("V", 0.19, 11.7, 280603),
    "GTC": ("V", 0.24, 14.9, 358478),
    "GTA": ("V", 0.12, 7.2, 173960),
    "GTG": ("V", 0.46, 28.6, 687491),
    "TCT": ("S", 0.18, 14.6, 351289),
    "TCC": ("S", 0.23, 18.3, 440468),
    "TCA": ("S", 0.15, 12.1, 290522),
    "TCG": ("S", 0.05, 4.3, 103373),
    "CCT": ("P", 0.28, 17.4, 418454),
    "CCC": ("P", 0.33, 20.3, 488352),
    "CCA": ("P", 0.28, 16.9, 406531),
    "CCG": ("P", 0.11, 6.9, 165877),
    "ACT": ("T", 0.25, 13.2, 317071),
    "ACC": ("T", 0.37, 19.4, 467021),
    "ACA": ("T", 0.27, 14.4, 346178),
    "ACG": ("T", 0.10, 5.5, 131823),
    "GCT": ("A", 0.26, 18.4, 442344),
    "GCC": ("A", 0.41, 29.1, 699443),
    "GCA": ("A", 0.22, 15.8, 379664),
    "GCG": ("A", 0.11, 7.6, 182868),
    "TAT": ("Y", 0.45, 12.3, 295743),
    "TAC": ("Y", 0.55, 15.1, 362810),
    "TAA": ("*", 0.30, 1.0, 24040),
    "TAG": ("*", 0.23, 0.8, 19232),
    "CAT": ("H", 0.42, 10.8, 259636),
    "CAC": ("H", 0.58, 14.8, 356197),
    "CAA": ("Q", 0.26, 12.3, 295743),
    "CAG": ("Q", 0.74, 34.9, 839005),
    "AAT": ("N", 0.47, 17.5, 420717),
    "AAC": ("N", 0.53, 19.5, 469287),
    "AAA": ("K", 0.43, 24.2, 581773),
    "AAG": ("K", 0.57, 32.0, 769289),
    "GAT": ("D", 0.46, 22.0, 528886),
    "GAC": ("D", 0.54, 26.2, 629856),
    "GAA": ("E", 0.42, 29.4, 706784),
    "GAG": ("E", 0.58, 40.8, 980430),
    "TGT": ("C", 0.45, 10.5, 252422),
    "TGC": ("C", 0.55, 12.8, 307715),
    "TGA": ("*", 0.47, 1.7, 40868),
    "TGG": ("W", 1.00, 13.3, 319735),
    "CGT": ("R", 0.09, 4.9, 117798),
    "CGC": ("R", 0.19, 11.3, 271656),
    "CGA": ("R", 0.10, 6.3, 151454),
    "CGG": ("R", 0.20, 12.0, 288483),
    "AGT": ("S", 0.16, 12.3, 295743),
    "AGC": ("S", 0.24, 19.1, 459168),
    "AGA": ("R", 0.22, 12.1, 290522),
    "AGG": ("R", 0.21, 12.0, 288483),
    "GGT": ("G", 0.16, 11.0, 264443),
    "GGC": ("G", 0.35, 23.5, 564944),
    "GGA": ("G", 0.25, 16.8, 403876),
    "GGG": ("G", 0.24, 16.3, 391856),
}

# Compute relative adaptiveness using shared utility
MOUSE_CODON_ADAPTIVENESS: dict[str, float] = compute_codon_adaptiveness(MOUSE_CODON_USAGE)

# Preferred (highest-frequency) codon for each amino acid
MOUSE_PREFERRED_CODONS: dict[str, str] = compute_preferred_codons(MOUSE_CODON_USAGE)

# ────────────────────────────────────────────────────────────
# Codon Pair Bias (CPB)
#
# CPB = log2(observed_frequency / expected_frequency)
# Positive CPB → over-represented pair (favoured for expression)
# Negative CPB → under-represented pair (disfavoured for expression)
#
# Sources:
#   Coleman et al. (2008) J Mol Evol 66:529-538
#   Quax et al. (2015) Nature Reviews Genetics 16:322-330
#   Müllner et al. (2020) Nucleic Acids Res 48:5864-5876
#
# Mouse codon pair bias is similar to human but not identical.
# Key differences arise from slightly different tRNA gene copy numbers
# affecting optimal codon choices — e.g. mouse has fewer tRNA-Arg(UCU)
# genes making AGA slightly rarer, and mouse tRNA-Gly(UCC) is relatively
# more abundant, shifting GGC preference upward.
# ────────────────────────────────────────────────────────────
MOUSE_CODON_PAIR_BIAS: dict[str, float] = {
    # ── Over-represented pairs (positive CPB) ──
    "CTG-CTG": 0.42,  # Leu-Leu   most common mouse codon pair
    "CTG-CAG": 0.36,  # Leu-Gln   GC-rich pair favoured
    "CAG-CTG": 0.34,  # Gln-Leu
    "CTG-GAG": 0.32,  # Leu-Glu
    "GAG-CTG": 0.30,  # Glu-Leu
    "ATG-CTG": 0.28,  # Met-Leu   start-proximal Leu bias
    "CTG-ATG": 0.26,  # Leu-Met
    "GAG-CAG": 0.24,  # Glu-Gln
    "CAG-GAG": 0.23,  # Gln-Glu
    "CTG-GGC": 0.22,  # Leu-Gly   mouse favours GGC for Gly
    "GGC-CTG": 0.21,  # Gly-Leu
    "GAG-GAG": 0.20,  # Glu-Glu
    "CTG-ACC": 0.19,  # Leu-Thr
    "ACC-CTG": 0.18,  # Thr-Leu
    "CAG-CAG": 0.17,  # Gln-Gln
    "GAG-AAC": 0.16,  # Glu-Asn   AAC preferred over AAT in mouse
    "ATG-ATG": 0.15,  # Met-Met
    "CTG-GAC": 0.14,  # Leu-Asp   GAC preferred over GAT
    "GAC-CTG": 0.13,  # Asp-Leu
    "GAG-CTG": 0.12,  # Glu-Leu   (slightly lower than CTG-GAG)
    "GCC-CTG": 0.11,  # Ala-Leu
    "CTG-GCC": 0.10,  # Leu-Ala
    "TAC-CTG": 0.09,  # Tyr-Leu   TAC preferred over TAT
    "CTG-TAC": 0.08,  # Leu-Tyr
    # ── Under-represented pairs (negative CPB) ──
    "CTA-ATA": -0.48,  # Leu(rare)-Ile(rare)
    "ATA-CTA": -0.46,  # Ile(rare)-Leu(rare)
    "TTA-TTA": -0.44,  # Leu(rare)-Leu(rare)  AT-rich disfavoured
    "ATA-ATA": -0.42,  # Ile(rare)-Ile(rare)
    "AGA-AGG": -0.40,  # Arg(rare)-Arg(rare)  mouse has fewer tRNA-Arg(UCU)
    "AGG-AGA": -0.38,  # Arg(rare)-Arg(rare)
    "CTA-CTA": -0.36,  # Leu(rare)-Leu(rare)
    "TTA-ATA": -0.34,  # Leu(rare)-Ile(rare)
    "ATA-TTA": -0.32,  # Ile(rare)-Leu(rare)
    "CGA-CGA": -0.30,  # Arg(rare)-Arg(rare)
    "AGA-ATA": -0.28,  # Arg(rare)-Ile(rare)
    "ATA-AGA": -0.26,  # Ile(rare)-Arg(rare)
    "TTA-AGG": -0.24,  # Leu(rare)-Arg(rare)
    "AGG-TTA": -0.22,  # Arg(rare)-Leu(rare)
    "TCG-ATA": -0.20,  # Ser(uncommon)-Ile(rare)
    "ATA-TCG": -0.18,  # Ile(rare)-Ser(uncommon)
    "CTA-AGG": -0.17,  # Leu(rare)-Arg(rare)
    "AGG-CTA": -0.16,  # Arg(rare)-Leu(rare)
    "TTA-CGA": -0.15,  # Leu(rare)-Arg(rare)
}

# ────────────────────────────────────────────────────────────
# Expression optimisation parameters for Mus musculus
#
# Tuned for transgene expression in mouse models and in vivo
# gene therapy vectors.  Parameters account for:
#   - Kozak consensus slightly different from human (mouse often
#     has GCC at -3 position rather than just GCCRCC)
#   - Polyadenylation signal AATAAA (dominant in mouse)
#   - CpG island avoidance important for reducing silencing in
#     mouse embryonic and stem cell contexts
#   - Splice site awareness to prevent cryptic splicing
# ────────────────────────────────────────────────────────────
MOUSE_EXPRESSION_OPTIMIZATION_PARAMS: dict[str, object] = {
    # Preferred UTR sequences
    "preferred_5utr_kozak": "GCCACCATGG",            # Mouse Kozak consensus
    "preferred_3utr_polya_signal": "AATAAA",          # Dominant poly-A signal
    # Rare-codon limits
    "max_consecutive_rare_codons": 3,                 # >3 risks ribosome stalling
    "rare_codon_threshold": 0.10,                     # Fraction < 10% → "rare"
    # GC content targets (mouse coding GC ~52.25%)
    "gc_content_target": 0.52,
    "gc_content_min": 0.30,
    "gc_content_max": 0.70,
    # Motifs to avoid
    "avoid_motifs": ["TTATTTT"],                      # mRNA instability element
    # Splicing and epigenetic considerations
    "splice_site_awareness": True,                    # Avoid cryptic splice sites
    "max_t_run": 6,                                   # Max consecutive T's
    "cpg_island_avoidance": True,                     # Reduce transgene silencing
}

# ────────────────────────────────────────────────────────────
# UTR Models for Mus musculus
#
# 5' UTR: Context for gene therapy vectors and transgene
#   expression in mouse models.  Includes Kozak consensus
#   and 5' cap-proximal sequence for efficient translation
#   initiation in mouse cells.
#
# 3' UTR: Stability elements for sustained in vivo expression.
#   Uses the mouse β-globin 3' UTR scaffold (well-characterised
#   for transgene stability) with AU-rich element avoidance
#   and the canonical AATAAA polyadenylation signal.
# ────────────────────────────────────────────────────────────
MOUSE_UTR_MODELS: dict[str, object] = {
    "5utr": {
        "description": "Mouse 5' UTR context for gene therapy vectors",
        "kozak_consensus": "GCCACCATGG",
        "cap_proximal_sequence": "GCUUCACCUCA",
        "upstream_intron_3ss": "CAG",              # 3' splice acceptor
        "optimized_leader_length": (40, 80),        # Recommended leader length (nt)
        "avoid_upstream_aug": True,                 # Prevent upstream ORFs
        "secondary_structure_max_delta_g": -30.0,   # kcal/mol; avoid highly structured 5' UTR
    },
    "3utr": {
        "description": "Mouse 3' UTR with stability elements for in vivo expression",
        "scaffold": "mouse_beta_globin_3utr",       # Well-characterised stability scaffold
        "polya_signal": "AATAAA",                   # Canonical poly-A signal
        "polya_signal_position": -18,               # nt upstream of cleavage site
        "cleavage_site_motif": "CA",                # Cleavage typically after CA dinucleotide
        "downstream_u_content": True,               # U-rich downstream element
        "avoid_ares": True,                         # Avoid AU-rich destabilising elements
        "are_motifs_to_avoid": ["ATTTA", "TTATTTT"],
        "stability_elements": ["CDE-I", "CDE-II"],  # CstF-binding downstream elements
        "optimized_length": (200, 600),              # Recommended 3' UTR length (nt)
    },
}
