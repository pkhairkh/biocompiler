"""Type stubs for biocompiler.types — core data structures."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any


# ────────────────────────────────────────────────────────────
# Enums
# ────────────────────────────────────────────────────────────

class SLOTMode(Enum):
    CONSERVATIVE: str
    VERIFIED: str
    PERMISSIVE: str


class Verdict(str, Enum):
    PASS: str
    LIKELY_PASS: str
    UNCERTAIN: str
    LIKELY_FAIL: str
    FAIL: str

    @property
    def confidence(self) -> float: ...

    @property
    def is_definite(self) -> bool: ...


# ────────────────────────────────────────────────────────────
# Five-valued logic
# ────────────────────────────────────────────────────────────

def five_valued_and(a: Verdict, b: Verdict) -> Verdict: ...
def five_valued_or(a: Verdict, b: Verdict) -> Verdict: ...
def three_valued_and(a: Verdict, b: Verdict) -> Verdict: ...
def three_valued_or(a: Verdict, b: Verdict) -> Verdict: ...
def combined_verdict(verdicts: list[Verdict]) -> Verdict: ...


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class PositionRange:
    start: int
    end: int

    def __len__(self) -> int: ...
    def overlaps(self, other: PositionRange) -> bool: ...
    def contains(self, position: int) -> bool: ...


@dataclass(frozen=True)
class Token:
    position: int
    element_type: str
    match_sequence: str
    score: float = 0.0
    frame: int | None = None
    strand: str = "+"

    @property
    def range(self) -> PositionRange: ...


@dataclass(frozen=True)
class SpliceIsoform:
    sequence: str
    exon_boundaries: list[tuple[int, int]]
    parse_path: list[str]
    score: float = 0.0


@dataclass
class TypeCheckResult:
    predicate: str
    verdict: Verdict
    derivation: list[dict] | None = None
    violation: str | None = None
    knowledge_gap: str | None = None

    @property
    def passed(self) -> bool: ...


@dataclass
class Certificate:
    version: str
    design_id: str
    sequence: str
    types: list[dict]
    provenance: dict
    hash_version: int = 2
    hash_algorithm: str = "sha256"

    def to_dict(self) -> dict: ...

    @classmethod
    def from_dict(cls, data: dict) -> Certificate: ...
