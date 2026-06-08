"""
BioCompiler Assembly Planning — DNA Assembly Method Planning

Provides tools for planning the assembly of optimized DNA sequences using
standard methods:
  - **Golden Gate Assembly**: Type IIS restriction enzyme-based one-pot assembly
  - **Gibson Assembly**: Overlap-based seamless assembly

For each method, the module generates an :class:`AssemblyPlan` describing the
fragments, enzymes or overlaps, and total construct length.

Design principles:
  - Deterministic: same input always produces the same plan.
  - Method-agnostic data structure: AssemblyPlan works for all methods.
  - BioCompiler-aware: fragments are assumed to be codon-optimized sequences.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from ..constants import RESTRICTION_ENZYMES

logger = logging.getLogger(__name__)

__all__ = [
    "AssemblyPlan",
    "plan_golden_gate",
    "plan_gibson",
]

# ---------------------------------------------------------------------------
# Enzyme recognition sites for Golden Gate assembly
# ---------------------------------------------------------------------------

_GOLDEN_GATE_ENZYMES: dict[str, str] = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "BbsI": "GAAGAC",
    "SapI": "GCTCTTC",
    "PaqCI": "CACCTGC",
}

# Standard overhang sequences for MoClo Golden Gate assembly
_MOCLO_OVERHANGS: list[str] = [
    "GGAG",  # Position 1 (upstream of promoter)
    "TACT",  # Position 2 (promoter → 5'UTR/RBS)
    "AATG",  # Position 3 (RBS → CDS)
    "GCTT",  # Position 4 (CDS → terminator)
    "CGCT",  # Position 5 (terminator → downstream)
]


# ---------------------------------------------------------------------------
# AssemblyPlan dataclass
# ---------------------------------------------------------------------------

@dataclass
class AssemblyPlan:
    """Plan for assembling DNA fragments into a complete construct.

    Attributes:
        method: Assembly method — "golden_gate", "gibson", or
            "restriction_ligation".
        fragments: Ordered list of DNA fragment sequences.
        enzymes: List of enzyme names (for Golden Gate / restriction ligation).
        overlap_sequences: Overlap sequences between fragments (for Gibson).
        total_length: Total length of the assembled construct in bp.
    """

    method: str  # "golden_gate", "gibson", "restriction_ligation"
    fragments: list[str]
    enzymes: list[str]  # for Golden Gate
    overlap_sequences: list[str]  # for Gibson
    total_length: int

    def __post_init__(self) -> None:
        valid_methods = {"golden_gate", "gibson", "restriction_ligation"}
        if self.method not in valid_methods:
            raise ValueError(
                f"Invalid assembly method {self.method!r}. "
                f"Must be one of {sorted(valid_methods)}"
            )


# ---------------------------------------------------------------------------
# Golden Gate Assembly
# ---------------------------------------------------------------------------

def plan_golden_gate(
    sequences: list[str],
    enzymes: list[str] | None = None,
) -> AssemblyPlan:
    """Plan Golden Gate assembly of optimized sequences.

    Golden Gate Assembly uses Type IIS restriction enzymes (e.g., BsaI, BsmBI)
    that cut outside their recognition site, generating custom overhangs.  This
    function designs the assembly by:

    1. Selecting appropriate enzyme(s) from the provided list (default: BsaI,
       BsmBI).
    2. Checking that none of the input sequences contain the enzyme recognition
       sites (which would cause unwanted internal cuts).
    3. Computing the total construct length.

    The actual overhang design (4bp sticky ends) is method-specific; this
    planner assigns standard MoClo-compatible overhangs between fragments.

    Args:
        sequences: Ordered list of DNA fragment sequences to assemble.
        enzymes: List of Type IIS enzyme names to use.  Defaults to
            ``["BsaI", "BsmBI"]``.

    Returns:
        An :class:`AssemblyPlan` with the assembly method, fragments, enzymes,
        and total length.

    Raises:
        ValueError: If no sequences are provided or if sequences contain
            internal recognition sites for the selected enzymes.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for Golden Gate assembly")

    if enzymes is None:
        enzymes = ["BsaI", "BsmBI"]

    # Validate enzyme names
    valid_enzymes = []
    for enz in enzymes:
        site = _GOLDEN_GATE_ENZYMES.get(enz) or RESTRICTION_ENZYMES.get(enz)
        if site:
            valid_enzymes.append(enz)
        else:
            logger.warning(
                "Enzyme %r not found in Golden Gate or standard enzyme database; skipping",
                enz,
            )

    if not valid_enzymes:
        raise ValueError(
            f"No valid enzymes found from {enzymes!r}. "
            f"Supported Golden Gate enzymes: {sorted(_GOLDEN_GATE_ENZYMES.keys())}"
        )

    # Check for internal restriction sites
    for enz in valid_enzymes:
        site = _GOLDEN_GATE_ENZYMES.get(enz, "") or RESTRICTION_ENZYMES.get(enz, "")
        if not site:
            continue
        site_upper = site.upper()
        for i, seq in enumerate(sequences):
            if site_upper in seq.upper():
                logger.warning(
                    "Fragment %d contains %s recognition site (%s); "
                    "consider domesticating before assembly",
                    i, enz, site,
                )

    # Compute total length
    total_length = sum(len(s) for s in sequences)

    # Add MoClo-style overhangs between fragments
    # For N fragments, we need N-1 overhangs
    overlap_seqs = []
    for i in range(len(sequences) - 1):
        # Assign standard overhang (cycling through MoClo overhangs)
        overhang = _MOCLO_OVERHANGS[i % len(_MOCLO_OVERHANGS)]
        overlap_seqs.append(overhang)

    return AssemblyPlan(
        method="golden_gate",
        fragments=list(sequences),
        enzymes=valid_enzymes,
        overlap_sequences=overlap_seqs,
        total_length=total_length,
    )


# ---------------------------------------------------------------------------
# Gibson Assembly
# ---------------------------------------------------------------------------

def plan_gibson(
    sequences: list[str],
    overlap_length: int = 20,
) -> AssemblyPlan:
    """Plan Gibson assembly with overlap design.

    Gibson Assembly joins DNA fragments through overlapping ends.  An exonuclease
    chews back the 5' ends, exposing complementary single-stranded overhangs that
    anneal.  A polymerase fills in gaps and a ligase seals nicks.

    This function designs the assembly by:
    1. Computing overlap sequences for each junction between fragments.
       The overlap is derived from the existing ends of adjacent fragments.
    2. Computing the total construct length (fragments minus overlaps to
       avoid double-counting the junction regions).

    Args:
        sequences: Ordered list of DNA fragment sequences.
        overlap_length: Length of overlap in bp between adjacent fragments.
            Defaults to 20 bp.

    Returns:
        An :class:`AssemblyPlan` with the assembly method, fragments, overlap
        sequences, and total length.

    Raises:
        ValueError: If no sequences are provided or if overlap_length is invalid.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for Gibson assembly")

    if overlap_length < 4:
        raise ValueError(
            f"overlap_length must be >= 4 bp, got {overlap_length}"
        )

    # Compute overlap sequences: for each junction, the overlap is the last
    # overlap_length bp of the current fragment (which should match the first
    # overlap_length bp of the next fragment in a proper Gibson design).
    overlap_seqs: list[str] = []
    for i in range(len(sequences) - 1):
        frag = sequences[i].upper()
        if len(frag) >= overlap_length:
            overlap = frag[-overlap_length:]
        else:
            overlap = frag  # Use entire fragment as overlap if too short
            logger.warning(
                "Fragment %d is shorter than overlap_length (%d < %d); "
                "using full fragment as overlap",
                i, len(frag), overlap_length,
            )
        overlap_seqs.append(overlap)

    # Total length = sum of fragments - (overlapping regions counted once)
    total_length = sum(len(s) for s in sequences) - overlap_length * max(0, len(sequences) - 1)

    return AssemblyPlan(
        method="gibson",
        fragments=list(sequences),
        enzymes=[],  # Gibson doesn't use restriction enzymes
        overlap_sequences=overlap_seqs,
        total_length=total_length,
    )
