"""
Biosecurity type aliases and data classes.

Separated from the main module to avoid circular imports and keep
the data-layer dependency-free.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

# ─────────────────────────────────────────────────────────────────────────────
# Type aliases
# ─────────────────────────────────────────────────────────────────────────────

RiskLevel = Literal["none", "low", "medium", "high", "critical"]
BiosecurityMode = Literal["enforce", "warn", "off"]
MatchType = Literal["exact", "fuzzy", "reverse_complement"]
StrandType = Literal["forward", "reverse"]


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class HazardMatch:
    """A single match against a hazardous sequence signature."""

    category: str
    name: str
    position: int
    matched_sequence: str
    confidence: float
    source: str
    # Extended fields for fuzzy and reverse complement matching
    match_type: MatchType = "exact"
    distance: int = 0
    strand: StrandType = "forward"
    substitutions: list[tuple[int, str, str]] = field(default_factory=list)


@dataclass
class BiosecurityReport:
    """Result of biosecurity screening for a protein/DNA sequence."""

    is_hazardous: bool
    risk_level: RiskLevel
    flagged_categories: list[str]
    matches: list[HazardMatch]
    recommendations: list[str]


@dataclass
class BiosecurityScreeningResult:
    """User-facing result of biosecurity screening for a protein sequence.

    This is the high-level result type returned by
    :func:`check_biosecurity_before_optimize`.  It provides a simple
    ``passed`` / ``failed`` boolean along with details about any
    flagged pathogens, risk levels, k-mer similarity scores, and
    match positions.
    """

    passed: bool
    flagged_pathogens: list[str] = field(default_factory=list)
    risk_levels: list[str] = field(default_factory=list)
    match_details: list[str] = field(default_factory=list)
    kmer_scores: dict[str, float] = field(default_factory=dict)
    screened_sequence_length: int = 0
    # Compatibility fields matching BiosecurityReport interface
    flagged_categories: list[str] = field(default_factory=list)
    matches: list[Any] = field(default_factory=list)
    recommendations: list[str] = field(default_factory=list)

    @property
    def is_hazardous(self) -> bool:
        """True if any hazard was detected."""
        return not self.passed

    @property
    def risk_level(self) -> str:
        """Return the highest risk level, or 'none' if no hazards."""
        if not self.risk_levels:
            return "none"
        priority = {"CRITICAL": 4, "HIGH": 3, "MEDIUM": 2, "LOW": 1, "none": 0}
        return max(self.risk_levels, key=lambda r: priority.get(r.upper(), 0))

    def __str__(self) -> str:
        if self.passed:
            return "BiosecurityScreeningResult: PASSED"
        pathogens = ", ".join(self.flagged_pathogens)
        levels = ", ".join(self.risk_levels)
        return f"BiosecurityScreeningResult: FAILED (pathogens=[{pathogens}], risk=[{levels}])"
