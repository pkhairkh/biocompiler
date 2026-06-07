"""
BioCompiler Local GC Constraints — Region-Specific GC Content Control

Provides fine-grained GC content constraints for specific regions of a DNA
sequence, as opposed to the global GC range enforced by the main optimizer.

Use cases:
  - Ensure a specific exon has GC content within a therapeutic range
  - Avoid extremely high/low GC in regulatory regions
  - Enforce uniform GC distribution across a gene to prevent
    amplification bias

This module is independent of the DNA Chisel compatibility layer and works
with BioCompiler's own codon tables and translation engine.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Sequence

from .scanner import gc_content
from biocompiler.shared.constants import CODON_TABLE, AA_TO_CODONS

logger = logging.getLogger(__name__)

__all__ = [
    "LocalGCConstraint",
    "LocalGCResult",
    "check_local_gc",
    "optimize_local_gc",
]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LocalGCConstraint:
    """A GC content constraint for a specific region of a DNA sequence.

    Attributes:
        region_start: 0-based inclusive start of the region (in bp).
        region_end: 0-based exclusive end of the region (in bp).
        gc_min: Minimum allowed GC fraction for the region [0.0, 1.0].
        gc_max: Maximum allowed GC fraction for the region [0.0, 1.0].
    """

    region_start: int
    region_end: int
    gc_min: float
    gc_max: float

    def __post_init__(self) -> None:
        if self.region_start < 0:
            raise ValueError(
                f"region_start must be >= 0, got {self.region_start}"
            )
        if self.region_end <= self.region_start:
            raise ValueError(
                f"region_end ({self.region_end}) must be > region_start "
                f"({self.region_start})"
            )
        if not (0.0 <= self.gc_min <= 1.0):
            raise ValueError(
                f"gc_min must be in [0.0, 1.0], got {self.gc_min}"
            )
        if not (0.0 <= self.gc_max <= 1.0):
            raise ValueError(
                f"gc_max must be in [0.0, 1.0], got {self.gc_max}"
            )
        if self.gc_min > self.gc_max:
            raise ValueError(
                f"gc_min ({self.gc_min}) must be <= gc_max ({self.gc_max})"
            )


@dataclass
class LocalGCResult:
    """Result of checking or optimizing local GC constraints.

    Attributes:
        satisfied: True if all constraints are satisfied.
        violations: List of (constraint, actual_gc) tuples for violated constraints.
        sequence: The (possibly optimized) DNA sequence.
    """

    satisfied: bool
    violations: list[tuple[LocalGCConstraint, float]]
    sequence: str


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def check_local_gc(
    dna: str,
    constraints: Sequence[LocalGCConstraint],
) -> LocalGCResult:
    """Check whether a DNA sequence satisfies all local GC constraints.

    For each constraint, the GC content of the sub-sequence
    ``dna[region_start:region_end]`` is computed and compared against the
    allowed range.  Violations are recorded but do not stop the check — all
    constraints are evaluated.

    Args:
        dna: DNA sequence (upper-cased internally).
        constraints: Iterable of :class:`LocalGCConstraint` instances.

    Returns:
        A :class:`LocalGCResult` indicating satisfaction and listing any
        violations.
    """
    dna = dna.upper()
    violations: list[tuple[LocalGCConstraint, float]] = []

    for c in constraints:
        end = min(c.region_end, len(dna))
        if c.region_start >= len(dna):
            # Region entirely beyond sequence — skip
            logger.warning(
                "Constraint region [%d:%d) is beyond sequence length %d; skipping",
                c.region_start, c.region_end, len(dna),
            )
            continue

        region = dna[c.region_start:end]
        if not region:
            continue

        gc = gc_content(region)
        if gc < c.gc_min or gc > c.gc_max:
            violations.append((c, gc))

    return LocalGCResult(
        satisfied=len(violations) == 0,
        violations=violations,
        sequence=dna,
    )


def optimize_local_gc(
    dna: str,
    protein: str,
    constraints: Sequence[LocalGCConstraint],
    codon_table: dict[str, str] | None = None,
) -> LocalGCResult:
    """Optimize a DNA sequence to satisfy local GC constraints.

    The function iterates codon-by-codon.  For each codon whose position
    falls within a violated region, it tries synonymous codons (same amino
    acid) to move the local GC content toward the target range.  Codons
    that are not within any constrained region are left untouched.

    The algorithm is greedy and may require multiple passes.  It stops
    after a maximum of 20 iterations or when all constraints are satisfied,
    whichever comes first.

    Args:
        dna: Original DNA coding sequence.
        protein: Target protein sequence (single-letter codes).
        constraints: Iterable of :class:`LocalGCConstraint` instances.
        codon_table: Optional codon → amino acid mapping.  Defaults to the
            standard genetic code from :mod:`biocompiler.constants`.

    Returns:
        A :class:`LocalGCResult` with the optimized sequence and any
        remaining violations.
    """
    if codon_table is None:
        codon_table = CODON_TABLE

    dna = dna.upper()
    protein = protein.upper()

    # Quick check: if already satisfied, return immediately
    result = check_local_gc(dna, constraints)
    if result.satisfied:
        return result

    # Build reverse mapping: amino acid -> list of codons
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in codon_table.items():
        if aa != "*":
            aa_to_codons.setdefault(aa, []).append(codon)

    # Merge with AA_TO_CODONS for completeness
    for aa, codons in AA_TO_CODONS.items():
        if aa not in aa_to_codons:
            aa_to_codons[aa] = list(codons)

    # Work with a mutable list of codons
    codons = [dna[i:i + 3] for i in range(0, len(dna) - 2, 3)]
    seq_str = "".join(codons)

    MAX_ITERATIONS = 20
    for iteration in range(MAX_ITERATIONS):
        check = check_local_gc(seq_str, constraints)
        if check.satisfied:
            return LocalGCResult(
                satisfied=True,
                violations=[],
                sequence=seq_str,
            )

        # For each violated constraint, try to fix codons in that region
        for constraint, actual_gc in check.violations:
            # Determine which codon indices fall within this region
            codon_indices = [
                i for i in range(len(codons))
                if i * 3 + 3 > constraint.region_start
                and i * 3 < constraint.region_end
            ]

            for ci in codon_indices:
                if ci >= len(protein):
                    continue

                current_codon = codons[ci]
                aa = protein[ci]

                # Get synonymous codons for this amino acid
                synonyms = aa_to_codons.get(aa, [])
                if len(synonyms) <= 1:
                    continue  # No alternatives

                # Compute current regional GC
                region_start = constraint.region_start
                region_end = min(constraint.region_end, len(seq_str))
                region_seq = seq_str[region_start:region_end]
                current_region_gc = gc_content(region_seq)

                # Try each synonym and pick the one that moves GC toward target
                target_gc = (constraint.gc_min + constraint.gc_max) / 2.0
                best_codon = current_codon
                best_diff = abs(current_region_gc - target_gc)

                for syn in synonyms:
                    if syn == current_codon:
                        continue

                    # Temporarily swap and recompute
                    old = codons[ci]
                    codons[ci] = syn
                    trial_seq = "".join(codons)
                    trial_region = trial_seq[region_start:region_end]
                    trial_gc = gc_content(trial_region)
                    trial_diff = abs(trial_gc - target_gc)

                    if trial_diff < best_diff:
                        best_codon = syn
                        best_diff = trial_diff

                    # Restore
                    codons[ci] = old

                # Apply the best codon if it improved
                if best_codon != current_codon:
                    codons[ci] = best_codon

        seq_str = "".join(codons)

    # Final check after all iterations
    final_check = check_local_gc(seq_str, constraints)
    return LocalGCResult(
        satisfied=final_check.satisfied,
        violations=final_check.violations,
        sequence=seq_str,
    )
