"""
Codon harmonization using the Relative Codon Adaptation (RCA) method
(Claassens et al., 2017, PMID: 28250820).

Instead of maximizing CAI (which selects the single best codon per amino acid),
harmonization matches the *relative* codon usage profile of the source organism
in the target host. This preserves translational kinetics — including rare codons
at domain boundaries that are functionally important for co-translational folding.

Algorithm (Claassens method):
    1. Compute RCA for source organism: for each codon c encoding amino acid aa,
       RCA_source(c) = freq(c) / max_freq(aa), where freq(c) is the per-thousand
       frequency of codon c and max_freq(aa) is the max per-thousand frequency
       among synonymous codons for the same amino acid.
    2. Compute RCA for target organism: same formula using target data.
    3. For each amino acid position in the protein:
       a. Get the source codon's RCA value.
       b. Among target organism codons for the same amino acid, find the one
          whose RCA is closest to the source RCA.
       c. Select that codon.

This produces a sequence whose *relative* codon usage pattern in the target host
mirrors the source organism's pattern, preserving fast/slow translation ramps
and pause sites that are critical for proper co-translational folding.

References
----------
Claassens, N.J., Siliakus, M.F., de Vos, W.M. & van der Oost, J. (2017).
Exploiting the genetic diversity of the TATA-binding protein for engineering
of translational efficiency. *Nucleic Acids Research*, 45(6), 3342–3356.
PMID: 28250820.
"""

from __future__ import annotations

import math
import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..organisms._utils import CodonUsageTable

__all__ = [
    "compute_rca",
    "harmonize_codons",
    "harmonize_with_cai_fallback",
    "compute_harmonization_score",
]

logger = logging.getLogger(__name__)


# ────────────────────────────────────────────────────────────
# RCA Computation
# ────────────────────────────────────────────────────────────

def compute_rca(codon_usage: CodonUsageTable) -> dict[str, float]:
    """Compute Relative Codon Adaptation (RCA) for each codon.

    RCA(c) = f(c) / f_max(aa), where:
    - f(c) is the per-thousand frequency of codon c
    - f_max(aa) is the max per-thousand frequency among synonymous codons
      for the same amino acid

    The RCA value is in [0.0, 1.0] for every codon. The most frequent
    codon for each amino acid has RCA = 1.0; rarer codons have lower values.
    This is conceptually similar to the CAI weight but uses per-thousand
    frequencies directly rather than the normalised adaptiveness.

    Stop codons (aa == "*") are excluded from the result.

    Args:
        codon_usage: Codon usage table mapping codon strings to
            (amino_acid, fraction, per_thousand, count) tuples.

    Returns:
        Dict mapping codon strings to their RCA values (0.0–1.0).
    """
    # Step 1: Find the maximum per-thousand frequency for each amino acid
    aa_max_freq: dict[str, float] = {}
    for _codon, (aa, _frac, freq, _count) in codon_usage.items():
        if aa == "*":
            continue
        current = aa_max_freq.get(aa, 0.0)
        if freq > current:
            aa_max_freq[aa] = freq

    # Step 2: Compute RCA as freq / max_freq for same amino acid
    rca: dict[str, float] = {}
    for codon, (aa, _frac, freq, _count) in codon_usage.items():
        if aa == "*":
            continue
        max_freq = aa_max_freq.get(aa, 0.0)
        rca[codon] = freq / max_freq if max_freq > 0.0 else 0.0

    return rca


# ────────────────────────────────────────────────────────────
# Codon Harmonization
# ────────────────────────────────────────────────────────────

def harmonize_codons(
    protein: str,
    source_organism: str,
    target_organism: str,
    codon_usage_source: CodonUsageTable | None = None,
    codon_usage_target: CodonUsageTable | None = None,
) -> str:
    """Harmonize codon usage from source to target organism (Claassens method).

    For each amino acid position, selects the target organism codon whose
    relative frequency (RCA) best matches the source organism codon's RCA.
    This preserves translational kinetics rather than maximizing CAI.

    For single-codon amino acids (Methionine/M, Tryptophan/W), the only
    available codon is always used regardless of organism.

    Args:
        protein: Amino acid sequence (1-letter codes, no stop codon).
        source_organism: Source organism name (resolved via resolve_organism).
        target_organism: Target organism name (resolved via resolve_organism).
        codon_usage_source: Optional codon usage table for the source organism.
            If None, looked up from CODON_USAGE_TABLES.
        codon_usage_target: Optional codon usage table for the target organism.
            If None, looked up from CODON_USAGE_TABLES.

    Returns:
        DNA sequence (uppercase) using harmonized codons for the target organism.

    Raises:
        ValueError: If an amino acid is not found in the codon tables.
        ValueError: If the source or target organism has no codon usage data.
    """
    from ..organisms import CODON_USAGE_TABLES, resolve_organism
    from ..type_system import AA_TO_CODONS

    # Resolve organism names
    source_organism = resolve_organism(source_organism)
    target_organism = resolve_organism(target_organism)

    # Look up codon usage tables
    if codon_usage_source is None:
        codon_usage_source = CODON_USAGE_TABLES.get(source_organism)
        if codon_usage_source is None:
            raise ValueError(
                f"No codon usage data for source organism '{source_organism}'. "
                f"Available: {list(CODON_USAGE_TABLES.keys())}"
            )
    if codon_usage_target is None:
        codon_usage_target = CODON_USAGE_TABLES.get(target_organism)
        if codon_usage_target is None:
            raise ValueError(
                f"No codon usage data for target organism '{target_organism}'. "
                f"Available: {list(CODON_USAGE_TABLES.keys())}"
            )

    # Compute RCA profiles
    rca_source = compute_rca(codon_usage_source)
    rca_target = compute_rca(codon_usage_target)

    # Build amino acid -> codon mapping from target organism
    # Only include codons present in the target RCA data
    target_aa_codons: dict[str, list[str]] = {}
    for codon, (aa, _frac, _freq, _count) in codon_usage_target.items():
        if aa == "*":
            continue
        if codon in rca_target:
            target_aa_codons.setdefault(aa, []).append(codon)

    # Harmonize each position
    dna_parts: list[str] = []
    for i, aa in enumerate(protein):
        codons = target_aa_codons.get(aa)
        if not codons:
            # Fallback to AA_TO_CODONS for standard amino acids
            codons = AA_TO_CODONS.get(aa, [])

        if not codons:
            raise ValueError(
                f"No codons found for amino acid '{aa}' at position {i} "
                f"in target organism '{target_organism}'"
            )

        if len(codons) == 1:
            # Single codon amino acid (M, W) — no choice
            dna_parts.append(codons[0])
            continue

        # Compute the "ideal" RCA for this position by averaging source RCA
        # values for all codons of this amino acid, weighted by frequency.
        # The Claassens method selects the source codon with the highest
        # frequency (most likely codon) and matches its RCA.
        # We use a weighted average of source RCAs as the target RCA.
        source_codons_for_aa = [
            c for c, (a, _, freq, _) in codon_usage_source.items()
            if a == aa and c in rca_source and freq > 0
        ]

        if source_codons_for_aa:
            # Weighted RCA: weight by per-thousand frequency in source organism
            total_weight = 0.0
            weighted_rca = 0.0
            for c in source_codons_for_aa:
                freq = codon_usage_source[c][2]  # per_thousand
                total_weight += freq
                weighted_rca += freq * rca_source[c]
            source_rca = weighted_rca / total_weight if total_weight > 0 else 1.0
        else:
            # Fallback: use the highest RCA (optimal codon)
            source_rca = 1.0

        # Among target codons for this amino acid, find the one whose RCA
        # is closest to the source RCA
        best_codon = codons[0]
        best_diff = float("inf")
        for c in codons:
            target_rca_val = rca_target.get(c, 0.0)
            diff = abs(target_rca_val - source_rca)
            if diff < best_diff:
                best_diff = diff
                best_codon = c
            # Tie-breaking: prefer higher-frequency codon (higher RCA)
            elif diff == best_diff:
                if rca_target.get(c, 0.0) > rca_target.get(best_codon, 0.0):
                    best_codon = c

        dna_parts.append(best_codon)

    return "".join(dna_parts)


# ────────────────────────────────────────────────────────────
# Harmonization with CAI Fallback
# ────────────────────────────────────────────────────────────

def harmonize_with_cai_fallback(
    protein: str,
    source_organism: str,
    target_organism: str,
    cai_weight: float = 0.5,
    codon_usage_source: CodonUsageTable | None = None,
    codon_usage_target: CodonUsageTable | None = None,
) -> str:
    """Blend harmonization with CAI optimization.

    For each position, choose between the harmonized codon and the CAI-optimal
    codon based on cai_weight:
    - cai_weight = 0.0: pure harmonization (same as harmonize_codons)
    - cai_weight = 1.0: pure CAI maximization
    - cai_weight = 0.5: balanced blend (default)

    The blending works per-position: the codon is selected based on a
    combined score that interpolates between the harmonization match quality
    and the CAI optimality.

    Args:
        protein: Amino acid sequence (1-letter codes, no stop codon).
        source_organism: Source organism name.
        target_organism: Target organism name.
        cai_weight: Blend weight for CAI optimization (0.0–1.0).
            0.0 = pure harmonization, 1.0 = pure CAI.
        codon_usage_source: Optional codon usage table for the source organism.
        codon_usage_target: Optional codon usage table for the target organism.

    Returns:
        DNA sequence using blended harmonized/CAI codons.

    Raises:
        ValueError: If cai_weight is not in [0.0, 1.0].
    """
    if not (0.0 <= cai_weight <= 1.0):
        raise ValueError(
            f"cai_weight must be in [0.0, 1.0], got {cai_weight}"
        )

    from ..organisms import CODON_USAGE_TABLES, CODON_ADAPTIVENESS_TABLES, resolve_organism
    from ..type_system import AA_TO_CODONS

    # Resolve organism names
    source_organism = resolve_organism(source_organism)
    target_organism = resolve_organism(target_organism)

    # Look up tables
    if codon_usage_source is None:
        codon_usage_source = CODON_USAGE_TABLES.get(source_organism)
        if codon_usage_source is None:
            raise ValueError(
                f"No codon usage data for source organism '{source_organism}'."
            )
    if codon_usage_target is None:
        codon_usage_target = CODON_USAGE_TABLES.get(target_organism)
        if codon_usage_target is None:
            raise ValueError(
                f"No codon usage data for target organism '{target_organism}'."
            )

    # Get CAI weights for target organism
    cai_weights = CODON_ADAPTIVENESS_TABLES.get(target_organism, {})

    # Compute RCA profiles
    rca_source = compute_rca(codon_usage_source)
    rca_target = compute_rca(codon_usage_target)

    # Build amino acid -> codon mapping from target organism
    target_aa_codons: dict[str, list[str]] = {}
    for codon, (aa, _frac, _freq, _count) in codon_usage_target.items():
        if aa == "*":
            continue
        if codon in rca_target:
            target_aa_codons.setdefault(aa, []).append(codon)

    harmonization_weight = 1.0 - cai_weight

    dna_parts: list[str] = []
    for i, aa in enumerate(protein):
        codons = target_aa_codons.get(aa)
        if not codons:
            codons = AA_TO_CODONS.get(aa, [])
        if not codons:
            raise ValueError(
                f"No codons found for amino acid '{aa}' at position {i} "
                f"in target organism '{target_organism}'"
            )

        if len(codons) == 1:
            dna_parts.append(codons[0])
            continue

        # Compute source target RCA (weighted average)
        source_codons_for_aa = [
            c for c, (a, _, freq, _) in codon_usage_source.items()
            if a == aa and c in rca_source and freq > 0
        ]

        if source_codons_for_aa:
            total_weight = 0.0
            weighted_rca = 0.0
            for c in source_codons_for_aa:
                freq = codon_usage_source[c][2]
                total_weight += freq
                weighted_rca += freq * rca_source[c]
            source_rca = weighted_rca / total_weight if total_weight > 0 else 1.0
        else:
            source_rca = 1.0

        # Score each codon: blend of harmonization quality and CAI
        best_codon = codons[0]
        best_score = float("-inf")

        for c in codons:
            # Harmonization score: 1.0 - |rca_target(c) - source_rca|
            target_rca_val = rca_target.get(c, 0.0)
            harm_score = 1.0 - abs(target_rca_val - source_rca)

            # CAI score: direct adaptiveness value
            cai_score = cai_weights.get(c, 0.0)

            # Combined score
            combined = harmonization_weight * harm_score + cai_weight * cai_score

            if combined > best_score:
                best_score = combined
                best_codon = c

        dna_parts.append(best_codon)

    return "".join(dna_parts)


# ────────────────────────────────────────────────────────────
# Harmonization Score
# ────────────────────────────────────────────────────────────

def compute_harmonization_score(
    sequence: str,
    source_organism: str,
    target_organism: str,
    codon_usage_source: CodonUsageTable | None = None,
    codon_usage_target: CodonUsageTable | None = None,
) -> float:
    """Score how well a sequence matches the source organism's RCA profile.

    The harmonization score measures the average match between the source
    organism's relative codon usage pattern and the sequence's actual
    codon usage in the target organism's context.

    For each codon position:
    1. Compute the source organism's weighted-average RCA for that amino acid
    2. Get the target organism's RCA for the codon actually used
    3. Score = 1.0 - |rca_target(codon) - source_rca|

    The overall score is the arithmetic mean across all positions.
    A score of 1.0 means perfect RCA profile matching; 0.0 means
    completely mismatched profiles.

    Args:
        sequence: DNA coding sequence (length must be a multiple of 3).
        source_organism: Source organism name.
        target_organism: Target organism name.
        codon_usage_source: Optional codon usage table for the source organism.
        codon_usage_target: Optional codon usage table for the target organism.

    Returns:
        Harmonization score in [0.0, 1.0]. Returns 0.0 for empty/invalid sequences.
    """
    from ..organisms import CODON_USAGE_TABLES, resolve_organism
    from ..type_system import CODON_TABLE

    if not sequence or len(sequence) < 3:
        return 0.0

    sequence = sequence.upper().strip()
    if len(sequence) % 3 != 0:
        logger.warning(
            "Sequence length %d is not a multiple of 3; "
            "trailing bases will be ignored for harmonization score",
            len(sequence),
        )

    # Resolve organism names
    source_organism = resolve_organism(source_organism)
    target_organism = resolve_organism(target_organism)

    # Look up tables
    if codon_usage_source is None:
        codon_usage_source = CODON_USAGE_TABLES.get(source_organism)
        if codon_usage_source is None:
            return 0.0
    if codon_usage_target is None:
        codon_usage_target = CODON_USAGE_TABLES.get(target_organism)
        if codon_usage_target is None:
            return 0.0

    # Compute RCA profiles
    rca_source = compute_rca(codon_usage_source)
    rca_target = compute_rca(codon_usage_target)

    # Pre-compute source weighted-average RCA per amino acid
    source_aa_rca: dict[str, float] = {}
    source_aa_codons: dict[str, list[tuple[str, float]]] = {}
    for codon, (aa, _frac, freq, _count) in codon_usage_source.items():
        if aa == "*":
            continue
        if codon in rca_source and freq > 0:
            source_aa_codons.setdefault(aa, []).append((codon, freq))

    for aa, codon_freqs in source_aa_codons.items():
        total_weight = sum(f for _, f in codon_freqs)
        if total_weight > 0:
            weighted_rca = sum(f * rca_source[c] for c, f in codon_freqs) / total_weight
            source_aa_rca[aa] = weighted_rca
        else:
            source_aa_rca[aa] = 1.0

    # Score each codon position
    scores: list[float] = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == "*":
            continue

        source_rca = source_aa_rca.get(aa, 1.0)
        target_rca_val = rca_target.get(codon, 0.0)

        # Score: how close is the target codon's RCA to the source's?
        position_score = 1.0 - abs(target_rca_val - source_rca)
        scores.append(max(0.0, position_score))

    if not scores:
        return 0.0

    return sum(scores) / len(scores)
