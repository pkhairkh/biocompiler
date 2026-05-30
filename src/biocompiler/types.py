"""
BioCompiler Core Types

Single canonical definition of all data structures.
No duplication — every other module imports from here.

Extended with:
- Frozen Token and SpliceIsoform for immutability guarantees
- Certificate validation in from_dict/to_dict
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class Verdict(str, Enum):
    """Three-valued logic for type-check verdicts (Kleene strong logic)."""
    PASS = "PASS"
    FAIL = "FAIL"
    UNCERTAIN = "UNCERTAIN"


def three_valued_and(a: Verdict, b: Verdict) -> Verdict:
    """Conjunction in three-valued logic (Kleene strong K3)."""
    if a == Verdict.FAIL or b == Verdict.FAIL:
        return Verdict.FAIL
    if a == Verdict.UNCERTAIN or b == Verdict.UNCERTAIN:
        return Verdict.UNCERTAIN
    return Verdict.PASS


def three_valued_or(a: Verdict, b: Verdict) -> Verdict:
    """Disjunction in three-valued logic (Kleene strong K3)."""
    if a == Verdict.PASS or b == Verdict.PASS:
        return Verdict.PASS
    if a == Verdict.UNCERTAIN or b == Verdict.UNCERTAIN:
        return Verdict.UNCERTAIN
    return Verdict.FAIL


def combined_verdict(verdicts: list[Verdict]) -> Verdict:
    """Compute the overall verdict from a list of verdicts using Kleene AND."""
    if not verdicts:
        return Verdict.UNCERTAIN
    result = verdicts[0]
    for v in verdicts[1:]:
        result = three_valued_and(result, v)
    return result


@dataclass(frozen=True)
class PositionRange:
    """Half-open interval [start, end) on a sequence."""
    start: int
    end: int

    def __len__(self) -> int:
        return self.end - self.start

    def overlaps(self, other: "PositionRange") -> bool:
        return self.start < other.end and other.start < self.end

    def contains(self, position: int) -> bool:
        return self.start <= position < self.end


@dataclass(frozen=True)
class Token:
    """An annotated region in a nucleotide sequence. Frozen for immutability."""
    position: int
    element_type: str
    match_sequence: str
    score: float = 0.0
    frame: int | None = None  # Reading frame (0, 1, 2), None if not frame-specific
    strand: str = "+"  # "+" or "-" for reverse complement matches

    @property
    def range(self) -> PositionRange:
        return PositionRange(self.position, self.position + len(self.match_sequence))


@dataclass
class SpliceIsoform:
    """A possible splice isoform computed by the NDFST."""
    sequence: str
    exon_boundaries: list[tuple[int, int]]
    parse_path: list[str]
    score: float = 0.0  # Path probability score


@dataclass
class TypeCheckResult:
    """Result of evaluating a type predicate against a sequence."""
    predicate: str
    verdict: Verdict
    derivation: Optional[list[dict]] = None
    violation: Optional[str] = None
    knowledge_gap: Optional[str] = None

    @property
    def passed(self) -> bool:
        return self.verdict == Verdict.PASS


_CERT_REQUIRED_KEYS = {"version", "design_id", "sequence", "types", "provenance"}


@dataclass
class Certificate:
    """A machine-checkable guarantee certificate."""
    version: str
    design_id: str
    sequence: str
    types: list[dict]
    provenance: dict

    def to_dict(self) -> dict:
        """Serialize to a plain dict (JSON-compatible)."""
        return {
            "version": self.version,
            "design_id": self.design_id,
            "sequence": self.sequence,
            "types": self.types,
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Certificate":
        """Deserialize from a plain dict with validation."""
        missing = _CERT_REQUIRED_KEYS - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize Certificate: missing keys {missing}"
            )
        return cls(
            version=data["version"],
            design_id=data["design_id"],
            sequence=data["sequence"],
            types=data["types"],
            provenance=data["provenance"],
        )
