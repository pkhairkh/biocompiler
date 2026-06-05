"""
BioCompiler IUPAC Ambiguous Base Support
=========================================

Provides utilities for handling IUPAC ambiguous DNA bases in input sequences.
DNAchisel supports all IUPAC codes; this module brings that capability into
BioCompiler's input pipeline.

IUPAC DNA ambiguity codes:
  R = A/G (purine)        S = G/C (strong)
  Y = C/T (pyrimidine)    W = A/T (weak)
  K = G/T (keto)          M = A/C (amino)
  B = C/G/T (not A)       D = A/G/T (not C)
  H = A/C/T (not G)       V = A/C/G (not T)
  N = A/C/G/T (any)

Usage:
  from biocompiler.iupac import resolve_ambiguous, is_ambiguous, expand_ambiguous

  # Resolve ambiguous bases before optimization
  resolved = resolve_ambiguous("ATGRYSW", strategy="most_common")

  # Check if a base is ambiguous
  is_ambiguous("R")  # True
  is_ambiguous("A")  # False

  # Expand into all possible concrete sequences
  expand_ambiguous("ATR")  # ["ATA", "ATG"]
"""

from __future__ import annotations

import itertools
import logging
from typing import Callable

__all__ = [
    "IUPAC_DNA",
    "resolve_ambiguous",
    "is_ambiguous",
    "expand_ambiguous",
    "has_ambiguous",
    "validate_iupac_sequence",
]

logger = logging.getLogger(__name__)

# ==============================================================================
# IUPAC DNA ambiguity code table
# ==============================================================================

IUPAC_DNA: dict[str, set[str]] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"},   # purine
    "Y": {"C", "T"},   # pyrimidine
    "S": {"G", "C"},   # strong
    "W": {"A", "T"},   # weak
    "K": {"G", "T"},   # keto
    "M": {"A", "C"},   # amino
    "B": {"C", "G", "T"},  # not A
    "D": {"A", "G", "T"},  # not C
    "H": {"A", "C", "T"},  # not G
    "V": {"A", "C", "G"},  # not T
    "N": {"A", "C", "G", "T"},  # any
}

# Reverse lookup: set of concrete bases -> IUPAC code
_CONCRETE_TO_IUPAC: dict[frozenset[str], str] = {
    frozenset(v): k for k, v in IUPAC_DNA.items()
}

# All valid IUPAC DNA characters
VALID_IUPAC_BASES: set[str] = set(IUPAC_DNA.keys())

# Default base frequencies for "most_common" strategy (human genome averages)
# Source: human genome GC ~41%, so A≈T≈29.5%, G≈C≈20.5%
_DEFAULT_BASE_FREQ: dict[str, float] = {
    "A": 0.295, "C": 0.205, "G": 0.205, "T": 0.295,
}


# ==============================================================================
# Core functions
# ==============================================================================

def is_ambiguous(base: str) -> bool:
    """Check whether a single base is an ambiguous IUPAC code.

    Args:
        base: Single character (case-insensitive).

    Returns:
        True if the base represents more than one concrete nucleotide.
    """
    return base.upper() in IUPAC_DNA and len(IUPAC_DNA[base.upper()]) > 1


def has_ambiguous(dna: str) -> bool:
    """Check whether a DNA sequence contains any ambiguous IUPAC bases.

    Args:
        dna: DNA sequence string (case-insensitive).

    Returns:
        True if any base in the sequence is ambiguous.
    """
    return any(is_ambiguous(b) for b in dna.upper())


def validate_iupac_sequence(dna: str) -> str:
    """Validate that a DNA sequence contains only valid IUPAC characters.

    Args:
        dna: DNA sequence string (case-insensitive).

    Returns:
        Uppercased, validated DNA string.

    Raises:
        ValueError: If any character is not a valid IUPAC DNA code.
    """
    dna = dna.upper().strip()
    invalid = set(dna) - VALID_IUPAC_BASES
    if invalid:
        raise ValueError(
            f"Invalid DNA bases found in sequence: {invalid}. "
            f"Supported: A, C, G, T, N and IUPAC ambiguity codes "
            f"(R, Y, S, W, K, M, B, D, H, V)."
        )
    return dna


def resolve_ambiguous(
    dna: str,
    strategy: str = "most_common",
    *,
    base_freq: dict[str, float] | None = None,
    gc_target: float | None = None,
    cai_table: dict[str, float] | None = None,
) -> str:
    """Resolve ambiguous bases in a DNA sequence using a specified strategy.

    This function replaces each ambiguous IUPAC base with a single concrete
    base according to the chosen strategy. Non-ambiguous bases (A, C, G, T)
    are left unchanged.

    Strategies:
      - "most_common" (default): Pick the most frequent concrete base according
        to organism base frequency data.  If no ``base_freq`` dict is provided,
        human genome averages are used (A≈T≈29.5%, G≈C≈20.5%).
      - "cai_optimal": Pick the concrete base that maximizes codon-wise CAI
        when the ambiguous base falls within a codon.  Requires ``cai_table``
        (a codon→adaptiveness mapping).  For positions where CAI cannot be
        determined (e.g., partial codons), falls back to most_common.
      - "gc_balanced": Pick the concrete base that keeps GC content closest
        to ``gc_target``.  If no target is given, uses 0.50 (equal GC/AT).
      - "first": Simply pick the first concrete base in alphabetical order
        from the IUPAC expansion.  Deterministic but biologically naive.

    Args:
        dna: DNA sequence possibly containing IUPAC ambiguous bases.
        strategy: Resolution strategy name (see above).
        base_freq: Optional per-base frequency dict for "most_common".
        gc_target: Target GC fraction for "gc_balanced" (default 0.50).
        cai_table: Codon adaptiveness table for "cai_optimal".

    Returns:
        A concrete DNA string containing only A, C, G, T characters.

    Raises:
        ValueError: If the strategy is unknown or sequence has invalid chars.
    """
    dna = validate_iupac_sequence(dna)

    if not has_ambiguous(dna):
        return dna  # Already concrete — nothing to resolve

    resolvers: dict[str, Callable[..., str]] = {
        "most_common": _resolve_most_common,
        "cai_optimal": _resolve_cai_optimal,
        "gc_balanced": _resolve_gc_balanced,
        "first": _resolve_first,
    }

    if strategy not in resolvers:
        raise ValueError(
            f"Unknown resolution strategy '{strategy}'. "
            f"Supported: {sorted(resolvers.keys())}"
        )

    return resolvers[strategy](
        dna,
        base_freq=base_freq,
        gc_target=gc_target,
        cai_table=cai_table,
    )


def _resolve_most_common(
    dna: str,
    *,
    base_freq: dict[str, float] | None = None,
    gc_target: float | None = None,
    cai_table: dict[str, float] | None = None,
) -> str:
    """Resolve each ambiguous base to the most frequent concrete base."""
    freq = base_freq or _DEFAULT_BASE_FREQ
    result = []
    for base in dna:
        if not is_ambiguous(base):
            result.append(base)
        else:
            candidates = IUPAC_DNA[base]
            # Pick the candidate with highest frequency
            best = max(candidates, key=lambda b: freq.get(b, 0.0))
            result.append(best)
    return "".join(result)


def _resolve_cai_optimal(
    dna: str,
    *,
    base_freq: dict[str, float] | None = None,
    gc_target: float | None = None,
    cai_table: dict[str, float] | None = None,
) -> str:
    """Resolve ambiguous bases by maximizing CAI at the codon level.

    For each ambiguous base, check if it falls within a codon. If so,
    enumerate all possible concrete codons and pick the one with highest
    CAI. If the ambiguous base doesn't fall in a complete codon, fall
    back to most_common.
    """
    if cai_table is None:
        logger.warning(
            "No CAI table provided for cai_optimal strategy; "
            "falling back to most_common"
        )
        return _resolve_most_common(dna, base_freq=base_freq)

    dna_list = list(dna)
    n = len(dna_list)

    # Process codon by codon
    for codon_start in range(0, n - 2, 3):
        codon_end = codon_start + 3
        codon_bases = dna_list[codon_start:codon_end]
        # Check if this codon has any ambiguous bases
        ambiguous_positions = [
            i for i, b in enumerate(codon_bases) if is_ambiguous(b)
        ]
        if not ambiguous_positions:
            continue

        # Enumerate all possible concrete codons
        candidate_lists = []
        for b in codon_bases:
            if is_ambiguous(b):
                candidate_lists.append(sorted(IUPAC_DNA[b]))
            else:
                candidate_lists.append([b])

        best_codon = None
        best_cai = -1.0
        for combo in itertools.product(*candidate_lists):
            concrete_codon = "".join(combo)
            cai_val = cai_table.get(concrete_codon, 0.0)
            if cai_val > best_cai:
                best_cai = cai_val
                best_codon = concrete_codon

        if best_codon is not None:
            dna_list[codon_start:codon_end] = list(best_codon)

    # Handle remaining bases outside complete codons (fall back to most_common)
    for i in range((n // 3) * 3, n):
        if is_ambiguous(dna_list[i]):
            freq = base_freq or _DEFAULT_BASE_FREQ
            candidates = IUPAC_DNA[dna_list[i]]
            dna_list[i] = max(candidates, key=lambda b: freq.get(b, 0.0))

    return "".join(dna_list)


def _resolve_gc_balanced(
    dna: str,
    *,
    base_freq: dict[str, float] | None = None,
    gc_target: float | None = None,
    cai_table: dict[str, float] | None = None,
) -> str:
    """Resolve ambiguous bases to keep GC content closest to target.

    For each ambiguous base, pick the concrete base that brings the
    running GC count closest to the target GC fraction.
    """
    target = gc_target if gc_target is not None else 0.50
    result = []
    gc_count = 0
    total_resolved = 0

    # First pass: count concrete bases' GC
    for base in dna:
        if not is_ambiguous(base):
            if base in ("G", "C"):
                gc_count += 1
            total_resolved += 1

    # Second pass: resolve ambiguous bases greedily
    for base in dna:
        if not is_ambiguous(base):
            result.append(base)
        else:
            candidates = sorted(IUPAC_DNA[base])
            best = None
            best_diff = float("inf")
            for c in candidates:
                test_gc = gc_count + (1 if c in ("G", "C") else 0)
                test_total = total_resolved + 1
                test_frac = test_gc / test_total if test_total > 0 else 0.0
                diff = abs(test_frac - target)
                if diff < best_diff:
                    best_diff = diff
                    best = c
            result.append(best)
            if best in ("G", "C"):
                gc_count += 1
            total_resolved += 1

    return "".join(result)


def _resolve_first(
    dna: str,
    *,
    base_freq: dict[str, float] | None = None,
    gc_target: float | None = None,
    cai_table: dict[str, float] | None = None,
) -> str:
    """Resolve each ambiguous base to the first (alphabetical) concrete base."""
    result = []
    for base in dna:
        if not is_ambiguous(base):
            result.append(base)
        else:
            result.append(sorted(IUPAC_DNA[base])[0])
    return "".join(result)


# ==============================================================================
# Expansion
# ==============================================================================

# Cap to avoid combinatorial explosion
_EXPANSION_CAP: int = 4096


def expand_ambiguous(dna: str) -> list[str]:
    """Expand all ambiguous bases into all possible concrete sequences.

    For example, ``"ATR"`` expands to ``["ATA", "ATG"]`` and
    ``"RY"`` expands to ``["AC", "AT", "GC", "GT"]``.

    If the total number of combinations exceeds 4096, a warning is logged
    and an empty list is returned (to prevent memory explosion).

    Args:
        dna: DNA sequence possibly containing IUPAC ambiguous bases.

    Returns:
        List of all concrete DNA strings (containing only A, C, G, T).

    Raises:
        ValueError: If the sequence contains invalid characters.
    """
    dna = validate_iupac_sequence(dna)

    if not has_ambiguous(dna):
        return [dna]  # No ambiguity — return the single concrete sequence

    # Calculate total combinations
    total_combos = 1
    for base in dna:
        total_combos *= len(IUPAC_DNA[base])
        if total_combos > _EXPANSION_CAP:
            logger.warning(
                "IUPAC expansion of '%s' exceeds %d combinations; "
                "returning empty list to avoid memory explosion.",
                dna[:20] + ("..." if len(dna) > 20 else ""),
                _EXPANSION_CAP,
            )
            return []

    # Build all combinations
    results = [""]
    for base in dna:
        concrete = sorted(IUPAC_DNA[base])
        results = [r + c for r in results for c in concrete]

    return results
