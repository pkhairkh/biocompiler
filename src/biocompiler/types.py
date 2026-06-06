"""
BioCompiler Core Types

Single canonical definition of all data structures.
No duplication — every other module imports from here.

Extended with:
- Frozen Token and SpliceIsoform for immutability guarantees
- Certificate validation in from_dict/to_dict
"""

from dataclasses import dataclass
from enum import Enum


class SLOTMode(Enum):
    """Controls how SLOT-dependent predicates are evaluated.

    SLOT (Subject to Limited Oracles and Tools) predicates depend on heuristic
    scanners or external tools that cannot be formally verified in the Lean4 model.

    CONSERVATIVE: Always return UNCERTAIN for SLOT predicates (matches Lean4 model).
    VERIFIED: Return PASS when verification conditions are met (tool available + result OK).
    PERMISSIVE: Return PASS with weaker evidence thresholds.
    """
    CONSERVATIVE = "conservative"  # Always UNCERTAIN for SLOT predicates (current default)
    VERIFIED = "verified"          # PASS when verification conditions met
    PERMISSIVE = "permissive"      # PASS with weaker evidence


class Verdict(str, Enum):
    """Five-valued logic for type-check verdicts (Kleene-style).

    Ordering: PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL
    """
    PASS = "PASS"
    LIKELY_PASS = "LIKELY_PASS"
    UNCERTAIN = "UNCERTAIN"
    LIKELY_FAIL = "LIKELY_FAIL"
    FAIL = "FAIL"

    @property
    def confidence(self) -> float:
        """Return a confidence score: 1.0 for PASS, 0.75 for LIKELY_PASS,
        0.5 for UNCERTAIN, 0.25 for LIKELY_FAIL, 0.0 for FAIL."""
        _confidence_map = {
            Verdict.PASS: 1.0,
            Verdict.LIKELY_PASS: 0.75,
            Verdict.UNCERTAIN: 0.5,
            Verdict.LIKELY_FAIL: 0.25,
            Verdict.FAIL: 0.0,
        }
        return _confidence_map[self]

    @property
    def is_definite(self) -> bool:
        """True for PASS/FAIL (definite verdicts), False for LIKELY_PASS/UNCERTAIN/LIKELY_FAIL."""
        return self in (Verdict.PASS, Verdict.FAIL)


# Internal ordering: PASS > LIKELY_PASS > UNCERTAIN > LIKELY_FAIL > FAIL
_VERDICT_ORDER: dict[Verdict, int] = {
    Verdict.PASS: 4,
    Verdict.LIKELY_PASS: 3,
    Verdict.UNCERTAIN: 2,
    Verdict.LIKELY_FAIL: 1,
    Verdict.FAIL: 0,}


def five_valued_and(a: Verdict, b: Verdict) -> Verdict:
    """Conjunction in five-valued logic (Kleene-style). AND takes the minimum."""
    if _VERDICT_ORDER[a] <= _VERDICT_ORDER[b]:        return a
    return b


def five_valued_or(a: Verdict, b: Verdict) -> Verdict:
    """Disjunction in five-valued logic (Kleene-style). OR takes the maximum."""
    if _VERDICT_ORDER[a] >= _VERDICT_ORDER[b]:        return a
    return b


# Backward-compatible aliases
three_valued_and = five_valued_and
three_valued_or = five_valued_or


def combined_verdict(verdicts: list[Verdict]) -> Verdict:
    """Compute the overall verdict from a list of verdicts using five-valued AND.

    The combined verdict is the weakest link: if any predicate fails,
    the overall result fails. If all pass, the result passes.
    """
    if not verdicts:
        return Verdict.UNCERTAIN
    result = verdicts[0]
    for v in verdicts[1:]:
        result = five_valued_and(result, v)
    return result


__all__ = [
    "SLOTMode",
    "Verdict",
    "five_valued_and",
    "five_valued_or",
    "three_valued_and",
    "three_valued_or",
    "combined_verdict",
    "PositionRange",
    "Token",
    "SpliceIsoform",
    "TypeCheckResult",
    "Certificate",
]


@dataclass(frozen=True)
class PositionRange:
    """Half-open interval [start, end) on a sequence."""
    start: int
    end: int

    def __len__(self) -> int:
        """Return the length of the half-open interval."""
        return self.end - self.start

    def overlaps(self, other: "PositionRange") -> bool:
        """Return True if this range overlaps with *other*."""
        return self.start < other.end and other.start < self.end

    def contains(self, position: int) -> bool:
        """Return True if *position* falls within [start, end)."""
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
        """Return the PositionRange covered by this token."""
        return PositionRange(self.position, self.position + len(self.match_sequence))


@dataclass(frozen=True)
class SpliceIsoform:
    """A possible splice isoform computed by the NDFST. Frozen for immutability."""
    sequence: str
    exon_boundaries: list[tuple[int, int]]
    parse_path: list[str]
    score: float = 0.0  # Path probability score

    def __repr__(self) -> str:
        return (
            f"SpliceIsoform(len={len(self.sequence)}, "
            f"exons={len(self.exon_boundaries)}, "
            f"path={self.parse_path}, score={self.score:.2f})"
        )


@dataclass
class TypeCheckResult:
    """Result of evaluating a type predicate against a sequence."""
    predicate: str
    verdict: Verdict
    derivation: list[dict] | None = None
    violation: str | None = None
    knowledge_gap: str | None = None

    @property
    def passed(self) -> bool:
        """True if the verdict is PASS or LIKELY_PASS."""
        return self.verdict in (Verdict.PASS, Verdict.LIKELY_PASS)

    def __repr__(self) -> str:
        return f"TypeCheckResult({self.predicate}={self.verdict.value})"


_CERT_REQUIRED_KEYS: set[str] = {"version", "design_id", "sequence", "types", "provenance"}


@dataclass
class Certificate:
    """A machine-checkable guarantee certificate.

    Attributes:
        version: BioCompiler version that generated this certificate.
        design_id: Hash-based identifier for the certified design.
            For hash_version=2, this includes sequence + predicate results +
            optimization parameters.  For hash_version=1 (legacy), this was
            the SHA-256 of the sequence only.
        sequence: The DNA sequence being certified.
        types: List of predicate result dicts (predicate name, verdict, etc.).
        provenance: Dict with tool version, timestamp, parameters, input_hash, etc.
        hash_version: Version of the hash computation used for design_id.
            1 = legacy (sequence-only hash), 2 = full hash covering sequence +
            sorted predicate results + key optimization parameters.
            Defaults to 2 for new certificates; certificates deserialized from
            legacy data default to 1 for backward compatibility.
    """
    version: str
    design_id: str
    sequence: str
    types: list[dict]
    provenance: dict
    hash_version: int = 2

    def to_dict(self) -> dict:
        """Serialize to a plain dict (JSON-compatible)."""
        d = {
            "version": self.version,
            "design_id": self.design_id,
            "sequence": self.sequence,
            "types": self.types,
            "provenance": self.provenance,
        }
        # Always include hash_version for v2+ certificates.
        # Legacy v1 certificates may not have it.
        if self.hash_version != 1:
            d["hash_version"] = self.hash_version
        return d

    @classmethod
    def from_dict(cls, data: dict) -> "Certificate":
        """Deserialize from a plain dict with validation."""
        missing = _CERT_REQUIRED_KEYS - set(data.keys())
        if missing:
            raise ValueError(
                f"Cannot deserialize Certificate: missing keys {missing}"
            )
        # Legacy certificates without hash_version default to v1
        hash_version = data.get("hash_version", 1)
        return cls(
            version=data["version"],
            design_id=data["design_id"],
            sequence=data["sequence"],
            types=data["types"],
            provenance=data["provenance"],
            hash_version=hash_version,
        )
