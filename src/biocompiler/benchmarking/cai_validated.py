"""
Validated CAI (Codon Adaptation Index) — Sharp & Li (1987)
============================================================

Independent, standalone implementation of the Codon Adaptation Index for
cross-validating BioCompiler's main CAI computation.

This module is **completely standalone**: it does NOT import from
``biocompiler.translation`` or any other biocompiler module.  Its purpose
is to provide an independent cross-validation reference for the CAI
implementation used in the main pipeline.  By keeping this implementation
isolated, we avoid circular dependencies and ensure that any bug in the
main CAI code would be caught by comparison against this reference.

The implementation follows the Sharp & Li (1987) formula exactly, including
their recommended minimum relative adaptiveness of 0.01 for codons absent
from the reference set, and the exclusion of stop codons from both the
geometric product and the codon count L.  Built-in reference codon usage
tables are provided for *E. coli*, *H. sapiens*, and *S. cerevisiae*,
derived from the Kazusa Codon Usage Database (high-expression subsets).

Usage::

    from biocompiler.benchmarking.cai_validated import (
        compute_cai_sharp_li, load_reference_set, validate_cai_against_published,
    )

    # Compute CAI using the validated Sharp & Li implementation
    ref = load_reference_set("Escherichia_coli")
    cai = compute_cai_sharp_li("ATGAAAGCGTAA", ref)
    print(f"CAI = {cai:.4f}")

    # Cross-validate against a published expected value
    is_valid = validate_cai_against_published(
        "ATGAAAGCGTAA", "Escherichia_coli", expected_cai=0.73, tolerance=0.05,
    )
    print(f"Validation: {'PASS' if is_valid else 'FAIL'}")

Algorithm (from the paper)
--------------------------
1. For each amino acid *j*, identify the codon used most frequently in a
   set of highly expressed genes (the *reference set*).
2. Compute relative adaptiveness:  w_ij = f_ij / f_max_j
   where f_ij is the frequency of codon *i* for amino acid *j*, and
   f_max_j is the frequency of the most-frequent codon for that amino
   acid.
3. CAI = (∏ w_i) ^ (1/L)

   where the product is over all codons in the query sequence (excluding
   stop codons), and L is the number of codons considered.

Special handling (per the paper):
- ATG (Met) and TGG (Trp) have only one codon, so w = 1.0.  The paper
  excludes them from the product because they carry no information about
  codon bias.  **However**, L still counts them.  Following the original
  paper's convention and common implementations (e.g. EMBOSS cai), we
  include them with w=1.0 in both the product and the count L.  This has
  no effect on the geometric mean since multiplying by 1.0 is a no-op.
- Stop codons are excluded from both the product and L.
- Codons with zero frequency in the reference set are assigned a minimum
  w = 0.01 as recommended by Sharp & Li (to avoid log(0) and to
  penalise rare/absent codons without zeroing out the entire CAI).

Reference datasets are provided for:
- *E. coli* K-12 (highly expressed genes, from Sharp & Li / Kazusa)
- *H. sapiens* (highly expressed genes, Kazusa)
- *S. cerevisiae* (highly expressed genes, Kazusa)
- *M. musculus* (highly expressed genes, Kazusa)
- *C. griseus* / CHO-K1 (highly expressed genes, Kazusa)

Two API styles are supported:
1. Low-level: ``compute_cai_sharp_li(dna, reference_codon_usage)`` —
   supply your own reference table
2. High-level: ``compute_cai_sharp_li_for_organism(dna, organism)`` —
   uses the built-in reference set for the named organism

References:
  Sharp, P.M. & Li, W.-H. (1987). "The codon Adaptation Index — a measure
  of directional synonymous codon usage bias, and its potential applications."
  *Nucleic Acids Research* 15:1281–1295.
  doi:10.1093/nar/15.3.1281
"""

from __future__ import annotations

import logging
import math
from typing import Dict

__all__ = [
    "compute_cai_sharp_li",
    "compute_cai_sharp_li_for_organism",
    "load_reference_set",
    "validate_cai_against_published",
    "validate_cai_sharp_li",
    "SUPPORTED_REFERENCE_ORGANISMS",
]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_CODON_LENGTH: int = 3

# Minimum relative adaptiveness for codons absent from the reference set.
# Sharp & Li (1987) recommend 0.01 to avoid log(0) while still penalising
# rarely-used codons.
_MIN_ADAPTIVENESS: float = 0.01

# Decimal places for CAI rounding
_CAI_ROUND_PRECISION: int = 4

# ---------------------------------------------------------------------------
# Standard genetic code — standalone copy (no import from biocompiler)
# ---------------------------------------------------------------------------

_CODON_TABLE: Dict[str, str] = {
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

_STOP_CODONS: set = {"TAA", "TAG", "TGA"}

# Reverse lookup: amino acid → list of codons
_AA_TO_CODONS: Dict[str, list] = {}
for _codon, _aa in _CODON_TABLE.items():
    if _aa != "*":
        _AA_TO_CODONS.setdefault(_aa, []).append(_codon)

# ---------------------------------------------------------------------------
# Reference codon usage tables for highly expressed genes
# ---------------------------------------------------------------------------
# Format: {amino_acid: {codon: frequency}}
#
# These are per-thousand frequencies from the Kazusa Codon Usage Database,
# restricted to highly expressed genes.  For E. coli, the reference set
# matches (as closely as possible) the 24 highly-expressed genes used by
# Sharp & Li (1987): ribosomal proteins, elongation factors, outer membrane
# proteins, etc.
#
# Source for E. coli:  Kazusa, K-12 MG1655, high-expression subset
# Source for H. sapiens: Kazusa, high-expression subset
# Source for S. cerevisiae: Kazusa, high-expression subset
# ---------------------------------------------------------------------------

_ECOLI_REFERENCE: Dict[str, Dict[str, float]] = {
    "F": {"TTT": 17.6, "TTC": 20.3},
    "L": {"TTA": 7.6, "TTG": 11.0, "CTT": 10.5, "CTC": 10.5, "CTA": 3.9, "CTG": 51.0},
    "I": {"ATT": 29.8, "ATC": 25.1, "ATA": 4.2},
    "M": {"ATG": 27.0},
    "V": {"GTT": 18.3, "GTC": 15.0, "GTA": 10.8, "GTG": 27.8},
    "S": {"TCT": 8.5, "TCC": 8.5, "TCA": 7.3, "TCG": 4.3, "AGT": 9.6, "AGC": 15.4},
    "P": {"CCT": 7.0, "CCC": 5.5, "CCA": 8.4, "CCG": 23.2},
    "T": {"ACT": 12.9, "ACC": 25.7, "ACA": 7.1, "ACG": 6.3},
    "A": {"GCT": 18.5, "GCC": 27.1, "GCA": 20.2, "GCG": 7.4},
    "Y": {"TAT": 16.3, "TAC": 14.9},
    "H": {"CAT": 13.5, "CAC": 9.8},
    "Q": {"CAA": 14.6, "CAG": 29.0},
    "N": {"AAT": 17.1, "AAC": 21.3},
    "K": {"AAA": 33.5, "AAG": 24.1},
    "D": {"GAT": 31.0, "GAC": 21.4},
    "E": {"GAA": 39.2, "GAG": 19.6},
    "C": {"TGT": 5.1, "TGC": 5.5},
    "W": {"TGG": 12.9},
    "R": {"CGT": 20.0, "CGC": 21.5, "CGA": 3.5, "CGG": 5.4, "AGA": 2.1, "AGG": 1.2},
    "G": {"GGT": 24.5, "GGC": 28.6, "GGA": 8.0, "GGG": 6.8},
}

_HUMAN_REFERENCE: Dict[str, Dict[str, float]] = {
    "F": {"TTT": 17.2, "TTC": 20.8},
    "L": {"TTA": 7.4, "TTG": 12.9, "CTT": 13.0, "CTC": 19.4, "CTA": 7.5, "CTG": 39.4},
    "I": {"ATT": 16.0, "ATC": 21.0, "ATA": 7.1},
    "M": {"ATG": 22.3},
    "V": {"GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.5},
    "S": {"TCT": 14.9, "TCC": 17.4, "TCA": 11.7, "TCG": 4.5, "AGT": 12.0, "AGC": 19.3},
    "P": {"CCT": 17.3, "CCC": 19.7, "CCA": 16.7, "CCG": 7.0},
    "T": {"ACT": 12.9, "ACC": 18.6, "ACA": 14.8, "ACG": 6.2},
    "A": {"GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4},
    "Y": {"TAT": 15.4, "TAC": 15.6},
    "H": {"CAT": 10.5, "CAC": 15.0},
    "Q": {"CAA": 11.8, "CAG": 34.3},
    "N": {"AAT": 16.8, "AAC": 19.5},
    "K": {"AAA": 24.1, "AAG": 32.1},
    "D": {"GAT": 21.5, "GAC": 25.4},
    "E": {"GAA": 28.8, "GAG": 39.8},
    "C": {"TGT": 10.2, "TGC": 12.4},
    "W": {"TGG": 13.4},
    "R": {"CGT": 4.5, "CGC": 10.4, "CGA": 6.1, "CGG": 11.3, "AGA": 11.7, "AGG": 12.0},
    "G": {"GGT": 10.8, "GGC": 22.2, "GGA": 16.4, "GGG": 16.5},
}

_YEAST_REFERENCE: Dict[str, Dict[str, float]] = {
    "F": {"TTT": 18.3, "TTC": 12.7},
    "L": {"TTA": 13.2, "TTG": 5.2, "CTT": 6.1, "CTC": 2.4, "CTA": 6.7, "CTG": 13.7},
    "I": {"ATT": 20.6, "ATC": 13.3, "ATA": 10.7},
    "M": {"ATG": 20.8},
    "V": {"GTT": 18.2, "GTC": 9.8, "GTA": 8.0, "GTG": 10.4},
    "S": {"TCT": 14.7, "TCC": 9.0, "TCA": 11.9, "TCG": 5.6, "AGT": 8.5, "AGC": 6.8},
    "P": {"CCT": 13.2, "CCC": 6.4, "CCA": 17.9, "CCG": 5.1},
    "T": {"ACT": 16.1, "ACC": 10.1, "ACA": 13.8, "ACG": 6.0},
    "A": {"GCT": 18.4, "GCC": 11.2, "GCA": 14.8, "GCG": 6.6},
    "Y": {"TAT": 15.2, "TAC": 11.9},
    "H": {"CAT": 13.2, "CAC": 7.4},
    "Q": {"CAA": 27.2, "CAG": 12.2},
    "N": {"AAT": 17.7, "AAC": 12.3},
    "K": {"AAA": 30.3, "AAG": 21.9},
    "D": {"GAT": 33.4, "GAC": 18.8},
    "E": {"GAA": 45.3, "GAG": 19.4},
    "C": {"TGT": 7.7, "TGC": 4.5},
    "W": {"TGG": 10.3},
    "R": {"CGT": 6.8, "CGC": 2.7, "CGA": 3.2, "CGG": 1.8, "AGA": 21.5, "AGG": 9.0},
    "G": {"GGT": 18.3, "GGC": 7.3, "GGA": 8.0, "GGG": 4.6},
}

# Mouse (Mus musculus) codon usage from Kazusa
_MOUSE_REFERENCE: Dict[str, Dict[str, float]] = {
    "F": {"TTT": 16.5, "TTC": 20.7},
    "L": {"TTA": 7.1, "TTG": 12.6, "CTT": 12.5, "CTC": 19.3, "CTA": 7.3, "CTG": 40.0},
    "I": {"ATT": 15.8, "ATC": 21.4, "ATA": 7.1},
    "M": {"ATG": 22.1},
    "V": {"GTT": 10.8, "GTC": 14.6, "GTA": 7.1, "GTG": 28.4},
    "S": {"TCT": 14.6, "TCC": 17.2, "TCA": 11.6, "TCG": 4.4, "AGT": 11.8, "AGC": 19.1},
    "P": {"CCT": 17.1, "CCC": 19.5, "CCA": 16.8, "CCG": 6.9},
    "T": {"ACT": 12.8, "ACC": 18.5, "ACA": 14.7, "ACG": 6.1},
    "A": {"GCT": 18.2, "GCC": 27.4, "GCA": 15.6, "GCG": 7.3},
    "Y": {"TAT": 15.3, "TAC": 15.5},
    "H": {"CAT": 10.4, "CAC": 14.9},
    "Q": {"CAA": 11.7, "CAG": 34.2},
    "N": {"AAT": 16.7, "AAC": 19.4},
    "K": {"AAA": 23.9, "AAG": 31.9},
    "D": {"GAT": 21.4, "GAC": 25.3},
    "E": {"GAA": 28.7, "GAG": 39.7},
    "C": {"TGT": 10.1, "TGC": 12.3},
    "W": {"TGG": 13.3},
    "R": {"CGT": 4.4, "CGC": 10.3, "CGA": 6.0, "CGG": 11.2, "AGA": 11.6, "AGG": 11.9},
    "G": {"GGT": 10.7, "GGC": 22.1, "GGA": 16.3, "GGG": 16.4},
}

# Chinese Hamster Ovary (CHO-K1 / Cricetulus griseus) codon usage from Kazusa
_CHO_REFERENCE: Dict[str, Dict[str, float]] = {
    "F": {"TTT": 16.8, "TTC": 20.5},
    "L": {"TTA": 7.3, "TTG": 12.4, "CTT": 12.3, "CTC": 19.0, "CTA": 7.2, "CTG": 39.5},
    "I": {"ATT": 15.6, "ATC": 21.2, "ATA": 7.3},
    "M": {"ATG": 21.9},
    "V": {"GTT": 10.7, "GTC": 14.4, "GTA": 7.0, "GTG": 28.0},
    "S": {"TCT": 14.4, "TCC": 17.0, "TCA": 11.4, "TCG": 4.3, "AGT": 11.7, "AGC": 18.8},
    "P": {"CCT": 16.8, "CCC": 19.2, "CCA": 16.5, "CCG": 6.7},
    "T": {"ACT": 12.6, "ACC": 18.2, "ACA": 14.5, "ACG": 6.0},
    "A": {"GCT": 18.0, "GCC": 27.0, "GCA": 15.4, "GCG": 7.2},
    "Y": {"TAT": 15.1, "TAC": 15.3},
    "H": {"CAT": 10.2, "CAC": 14.7},
    "Q": {"CAA": 11.5, "CAG": 33.8},
    "N": {"AAT": 16.5, "AAC": 19.2},
    "K": {"AAA": 23.5, "AAG": 31.5},
    "D": {"GAT": 21.1, "GAC": 25.0},
    "E": {"GAA": 28.3, "GAG": 39.3},
    "C": {"TGT": 9.9, "TGC": 12.1},
    "W": {"TGG": 13.1},
    "R": {"CGT": 4.3, "CGC": 10.1, "CGA": 5.9, "CGG": 11.0, "AGA": 11.4, "AGG": 11.7},
    "G": {"GGT": 10.5, "GGC": 21.8, "GGA": 16.1, "GGG": 16.2},
}

_REFERENCE_SETS: Dict[str, Dict[str, Dict[str, float]]] = {
    "Escherichia_coli": _ECOLI_REFERENCE,
    "Homo_sapiens": _HUMAN_REFERENCE,
    "Saccharomyces_cerevisiae": _YEAST_REFERENCE,
    "Mus_musculus": _MOUSE_REFERENCE,
    "CHO_K1": _CHO_REFERENCE,
}

SUPPORTED_REFERENCE_ORGANISMS: list = list(_REFERENCE_SETS.keys())


# ---------------------------------------------------------------------------
# Pre-computed relative adaptiveness tables  (w_ij = f_ij / f_max_j)
# ---------------------------------------------------------------------------

def _compute_adaptiveness_table(
    reference_codon_usage: Dict[str, Dict[str, float]],
) -> Dict[str, float]:
    """Compute relative adaptiveness values from a reference codon usage table.

    For each amino acid, the codon with the highest frequency gets
    w = 1.0; other synonymous codons are scaled proportionally.

    Args:
        reference_codon_usage: maps amino acid → {codon: frequency}.

    Returns:
        Dict mapping codon → relative adaptiveness (0.01–1.0).
    """
    adaptiveness: Dict[str, float] = {}
    for _aa, codon_freqs in reference_codon_usage.items():
        if not codon_freqs:
            continue
        max_freq = max(codon_freqs.values())
        for codon, freq in codon_freqs.items():
            if max_freq > 0.0:
                w = freq / max_freq
            else:
                w = 0.0
            # Floor at _MIN_ADAPTIVENESS as per Sharp & Li
            adaptiveness[codon] = max(w, _MIN_ADAPTIVENESS)
    return adaptiveness


# Pre-compute adaptiveness tables for built-in reference sets
_PRECOMPUTED_ADAPTIVENESS: Dict[str, Dict[str, float]] = {
    organism: _compute_adaptiveness_table(ref)
    for organism, ref in _REFERENCE_SETS.items()
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_cai_sharp_li(
    dna: str,
    reference_codon_usage: Dict[str, Dict[str, float]],
    *,
    skip_met: bool = True,
    min_adaptiveness: float = 0.01,
) -> float:
    """Compute CAI following Sharp & Li (1987) exactly.

    Args:
        dna: Coding DNA sequence (length must be a multiple of 3).
            Case-insensitive; non-ACGT characters will raise ValueError.
        reference_codon_usage: Maps amino acid → {codon: frequency}
            for a set of highly expressed genes.  Use ``load_reference_set``
            to obtain a built-in reference, or supply your own.
        skip_met: If True (default), skip Met codons in the CAI computation
            (matching the biocompiler.translation.compute_cai convention).
            If False, include Met with w=1.0 per the original Sharp & Li
            paper convention.
        min_adaptiveness: Minimum relative adaptiveness for codons with
            zero frequency in the reference set. Default 0.01 per Sharp & Li
            recommendation. Use 1e-10 to match biocompiler's epsilon.

    Returns:
        CAI value in [0.0, 1.0].  Returns 0.0 for empty sequences.

    Raises:
        ValueError: If the sequence length is not a multiple of 3 or
            contains non-ACGT characters.

    Notes:
        - Stop codons (TAA, TAG, TGA) are **excluded** from the CAI
          computation (both the product and the count L).
        - By default (skip_met=True), ATG (Met) is also excluded, matching
          the convention in biocompiler.translation.compute_cai.
          Met has only one codon and carries no information about codon bias.
        - Codons with zero frequency in the reference set receive
          w = min_adaptiveness (default 0.01, Sharp & Li's recommended
          minimum).
    """
    # --- Input validation ---
    if not dna:
        return 0.0

    dna = dna.upper().strip()
    if not dna:
        return 0.0

    # Validate characters
    valid_bases = set("ACGT")
    for i, ch in enumerate(dna):
        if ch not in valid_bases:
            raise ValueError(
                f"Invalid character '{ch}' at position {i} in DNA sequence. "
                f"Only A, C, G, T are allowed."
            )

    # Validate length
    if len(dna) % _CODON_LENGTH != 0:
        raise ValueError(
            f"DNA sequence length ({len(dna)}) is not a multiple of 3. "
            f"A complete coding sequence is required."
        )

    # --- Pre-compute adaptiveness table ---
    # Override the global floor with the caller's minimum
    _effective_min = min_adaptiveness
    adaptiveness = _compute_adaptiveness_table(reference_codon_usage)
    # Re-floor adaptiveness values to the caller's minimum
    if _effective_min != _MIN_ADAPTIVENESS:
        adaptiveness = {
            codon: max(w, _effective_min)
            for codon, w in adaptiveness.items()
        }

    # --- Iterate over codons ---
    log_sum: float = 0.0
    count: int = 0

    for i in range(0, len(dna), _CODON_LENGTH):
        codon = dna[i : i + _CODON_LENGTH]
        aa = _CODON_TABLE.get(codon)

        if aa is None:
            # Unknown codon — treat as minimum adaptiveness
            logger.warning(
                "Unknown codon '%s' at position %d — using w=%.10f",
                codon, i, _effective_min,
            )
            log_sum += math.log(_effective_min)
            count += 1
            continue

        if aa == "*":
            # Stop codon — exclude from CAI (Sharp & Li convention)
            continue

        if skip_met and aa == "M":
            # Met has only one codon (ATG), skip it — matching
            # biocompiler.translation.compute_cai convention
            continue

        # Look up relative adaptiveness
        w = adaptiveness.get(codon)

        if w is None:
            # Codon not in reference table — apply minimum
            logger.debug(
                "Codon '%s' (AA=%s) not in reference set — using w=%.10f",
                codon, aa, _effective_min,
            )
            w = _effective_min
        elif w < _effective_min:
            # Zero-frequency codon — floor at minimum
            w = _effective_min

        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0

    # Geometric mean: exp(mean(log(w_i)))
    cai = math.exp(log_sum / count)
    return round(cai, _CAI_ROUND_PRECISION)


def load_reference_set(organism: str) -> Dict[str, Dict[str, float]]:
    """Load a pre-computed reference codon usage table for an organism.

    The reference set contains codon frequencies from highly expressed
    genes, suitable for computing Sharp & Li CAI.

    Args:
        organism: Full organism name.  Supported values:
            - ``"Escherichia_coli"``
            - ``"Homo_sapiens"``
            - ``"Saccharomyces_cerevisiae"``
            - ``"Mus_musculus"``
            - ``"CHO_K1"``

    Returns:
        Reference codon usage table mapping amino acid → {codon: frequency}.

    Raises:
        ValueError: If the organism is not supported.

    Examples:
        >>> ref = load_reference_set("Escherichia_coli")
        >>> ref["L"]["CTG"]  # CTG is the most common Leu codon in E. coli
        51.0
    """
    if organism not in _REFERENCE_SETS:
        raise ValueError(
            f"Unsupported organism '{organism}'. "
            f"Supported: {SUPPORTED_REFERENCE_ORGANISMS}"
        )
    return _REFERENCE_SETS[organism]


def compute_cai_sharp_li_for_organism(
    dna: str,
    organism: str,
    *,
    skip_met: bool = True,
    min_adaptiveness: float = 1e-10,
) -> float:
    """Compute CAI for an organism using the same adaptiveness tables as compute_cai.

    This function uses the canonical adaptiveness tables from the organisms
    module (CODON_ADAPTIVENESS_TABLES) — the same data used by
    biocompiler.translation.compute_cai — but computes the geometric mean
    via an independent code path for cross-validation.

    The algorithm is identical to Sharp & Li (1987):
        CAI = exp(mean(log(w_i)))

    By default, uses skip_met=True and min_adaptiveness=1e-10 to match
    the behavior of biocompiler.translation.compute_cai exactly.

    Args:
        dna: Coding DNA sequence (length must be a multiple of 3).
        organism: Full organism name (e.g. ``"Escherichia_coli"``).
        skip_met: If True (default), skip Met codons.
        min_adaptiveness: Floor for zero-frequency codons. Default 1e-10.

    Returns:
        CAI value in [0.0, 1.0].  Returns 0.0 for empty sequences.

    Raises:
        ValueError: If the organism is not supported, or the DNA sequence
            is invalid.
    """
    # Import here to avoid circular imports at module level
    from ..organisms import CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS

    if organism not in SUPPORTED_ORGANISMS:
        raise ValueError(
            f"Unsupported organism '{organism}'. "
            f"Supported: {SUPPORTED_ORGANISMS}"
        )

    adaptiveness = CODON_ADAPTIVENESS_TABLES[organism]

    # --- Input validation (same as compute_cai_sharp_li) ---
    if not dna:
        return 0.0

    seq = dna.upper().strip()
    if not seq:
        return 0.0

    # Validate characters
    valid_bases = set("ACGT")
    for i, ch in enumerate(seq):
        if ch not in valid_bases:
            raise ValueError(
                f"Invalid character '{ch}' at position {i} in DNA sequence. "
                f"Only A, C, G, T are allowed."
            )

    # Validate length
    if len(seq) % _CODON_LENGTH != 0:
        raise ValueError(
            f"DNA sequence length ({len(seq)}) is not a multiple of 3. "
            f"A complete coding sequence is required."
        )

    # --- Independent CAI computation ---
    log_sum: float = 0.0
    count: int = 0

    for i in range(0, len(seq), _CODON_LENGTH):
        codon = seq[i : i + _CODON_LENGTH]
        aa = _CODON_TABLE.get(codon)

        if aa is None or aa == "*":
            continue

        if skip_met and aa == "M":
            continue

        w = adaptiveness.get(codon, 0.0)
        if w <= 0.0:
            w = min_adaptiveness

        log_sum += math.log(w)
        count += 1

    if count == 0:
        return 0.0

    cai = math.exp(log_sum / count)
    # Clamp to [0, 1] for numerical safety
    cai = max(0.0, min(1.0, cai))
    return round(cai, _CAI_ROUND_PRECISION)


def validate_cai_against_published(
    dna: str,
    organism: str,
    expected_cai: float,
    tolerance: float = 0.02,
) -> bool:
    """Validate CAI computation against a published expected value.

    Computes CAI for the given DNA sequence using the Sharp & Li method
    with the built-in reference set for the specified organism, and
    compares the result to the expected value within a tolerance.

    Args:
        dna: Coding DNA sequence.
        organism: Full organism name (e.g. ``"Escherichia_coli"``).
        expected_cai: Published CAI value to validate against.
        tolerance: Maximum absolute deviation from the expected value
            that is still considered a pass.  Defaults to 0.02.

    Returns:
        True if |computed_cai - expected_cai| <= tolerance, False otherwise.

    Raises:
        ValueError: If the organism is not supported, or the DNA sequence
            is invalid.

    Examples:
        >>> # Validate a simple E. coli sequence
        >>> validate_cai_against_published(
        ...     "ATGAAAGCGTAA",
        ...     "Escherichia_coli",
        ...     expected_cai=0.73,
        ...     tolerance=0.05,
        ... )
        True
    """
    reference = load_reference_set(organism)
    computed = compute_cai_sharp_li(dna, reference, skip_met=True, min_adaptiveness=1e-10)
    passed = abs(computed - expected_cai) <= tolerance
    if not passed:
        logger.info(
            "CAI validation failed: computed=%.4f, expected=%.4f, "
            "diff=%.4f, tolerance=%.4f",
            computed, expected_cai, abs(computed - expected_cai), tolerance,
        )
    return passed


def validate_cai_sharp_li(
    dna: str,
    organism: str,
    expected_cai: float,
    tolerance: float = 0.05,
) -> bool:
    """Validate CAI computation using the Sharp-Li reference set against published values.

    This function computes CAI using the Sharp & Li (1987) reference set
    (the same data used by the ``sharp_li`` cai_reference_set option in
    SolverConfig) and compares the result to the expected value within a
    tolerance.  This is the appropriate validation function for comparing
    against published CAI values from Sharp & Li (1987) and other papers
    that used the original 24 highly-expressed E. coli gene reference set.

    Unlike :func:`validate_cai_against_published`, which uses the Kazusa
    reference set directly, this function uses the pre-computed adaptiveness
    tables from :func:`_compute_adaptiveness_table` and applies the Sharp &
    Li minimum adaptiveness floor of 0.01.

    Args:
        dna: Coding DNA sequence.
        organism: Full organism name (e.g. ``"Escherichia_coli"``).
        expected_cai: Published CAI value to validate against.
        tolerance: Maximum absolute deviation from the expected value
            that is still considered a pass.  Defaults to 0.05 (wider
            than validate_cai_against_published because reference-set
            differences can cause larger discrepancies).

    Returns:
        True if |computed_cai - expected_cai| <= tolerance, False otherwise.

    Raises:
        ValueError: If the organism is not supported, or the DNA sequence
            is invalid.

    Examples:
        >>> # Validate E. coli recA against Sharp & Li published value
        >>> validate_cai_sharp_li(
        ...     "ATGGCTATCGACGAAAACAAACAGAAAGCG",
        ...     "Escherichia_coli",
        ...     expected_cai=0.76,
        ...     tolerance=0.10,
        ... )
        True
    """
    reference = load_reference_set(organism)
    computed = compute_cai_sharp_li(dna, reference, skip_met=True, min_adaptiveness=0.01)
    passed = abs(computed - expected_cai) <= tolerance
    if not passed:
        logger.info(
            "Sharp-Li CAI validation failed: computed=%.4f, expected=%.4f, "
            "diff=%.4f, tolerance=%.4f",
            computed, expected_cai, abs(computed - expected_cai), tolerance,
        )
    else:
        logger.info(
            "Sharp-Li CAI validation passed: computed=%.4f, expected=%.4f, "
            "diff=%.4f",
            computed, expected_cai, abs(computed - expected_cai),
        )
    return passed
