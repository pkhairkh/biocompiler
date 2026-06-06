"""Type stubs for biocompiler.biosecurity — Hazardous sequence screening and biosecurity gate."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional


# ────────────────────────────────────────────────────────────
# Type aliases
# ────────────────────────────────────────────────────────────

RiskLevel = Literal["none", "low", "medium", "high", "critical"]
BiosecurityMode = Literal["enforce", "warn", "off"]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class HazardMatch:
    category: str
    name: str
    position: int
    matched_sequence: str
    confidence: float
    source: str


@dataclass
class BiosecurityReport:
    is_hazardous: bool
    risk_level: RiskLevel
    flagged_categories: list[str]
    matches: list[HazardMatch]
    recommendations: list[str]


# ────────────────────────────────────────────────────────────
# Core screening function
# ────────────────────────────────────────────────────────────

def screen_hazardous_sequence(protein: str, dna: str = ...) -> BiosecurityReport: ...


# ────────────────────────────────────────────────────────────
# Biosecurity mode
# ────────────────────────────────────────────────────────────

def get_biosecurity_mode() -> BiosecurityMode: ...


# ────────────────────────────────────────────────────────────
# Integration hook
# ────────────────────────────────────────────────────────────

def check_biosecurity_before_optimize(
    protein: str,
    organism: str = ...,
    dna: str = ...,
    biosecurity_mode: Optional[BiosecurityMode] = ...,
) -> BiosecurityReport: ...


# ────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────

def sig_risk_for_match(match: HazardMatch) -> str: ...
