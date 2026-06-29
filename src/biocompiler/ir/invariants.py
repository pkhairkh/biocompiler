"""
BioCompiler IR — Invariant Checkers
===================================

Each ``check_lN_invariants`` function validates that an IR object
satisfies its level's invariants — the semantic properties that
lowering passes must preserve.

These checks are deliberately **separate** from the dataclass
constructors: a malformed IR object can still be constructed (useful
for diagnostics), and then either pass or fail an explicit invariant
check.  This mirrors how LLVM lets you build malformed IR for testing
and then run the ``verify`` pass to detect it.

Each checker either returns ``True`` (if all invariants hold) or
raises :class:`IRError` with a descriptive message.
"""

from __future__ import annotations

from .types import (
    IR_L0_GenomicDNA,
    IR_L1_PreMRNA,
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    IR_L4_FoldedProtein,
    IRError,
    IRLevel,
)


# Canonical sets of valid symbols at each sequence level.
_VALID_DNA_BASES = frozenset("ACGTN")
_VALID_RNA_BASES = frozenset("ACGUN")
_VALID_AA = frozenset("ACDEFGHIKLMNPQRSTVWY*")
_RNA_STOP_CODONS = frozenset({"UAA", "UAG", "UGA"})


def check_l0_invariants(ir: IR_L0_GenomicDNA) -> bool:
    """L0 invariants.

    * Non-empty sequence.
    * Only DNA bases ``{A, C, G, T, N}`` (case-insensitive).
    * Every region is within ``[0, len(sequence)]`` and well-formed
      (``start >= 0``, ``end <= len``, ``start < end``).
    """
    if not ir.sequence:
        raise IRError(IRLevel.L0, "empty sequence")
    seq_upper = ir.sequence.upper()
    if not set(seq_upper) <= _VALID_DNA_BASES:
        bad = sorted(set(seq_upper) - _VALID_DNA_BASES)
        raise IRError(IRLevel.L0, f"invalid bases: {bad}")

    seq_len = len(ir.sequence)
    for r in ir.regions:
        if r.start < 0:
            raise IRError(IRLevel.L0, f"region start < 0: {r}")
        if r.end > seq_len:
            raise IRError(
                IRLevel.L0,
                f"region end {r.end} > sequence length {seq_len}: {r}",
            )
        if r.start >= r.end:
            raise IRError(
                IRLevel.L0,
                f"region start {r.start} >= end {r.end}: {r}",
            )
    return True


def check_l1_invariants(ir: IR_L1_PreMRNA) -> bool:
    """L1 invariants.

    * Non-empty sequence.
    * Only RNA bases ``{A, C, G, U, N}`` (case-insensitive).
    * Region annotations (if any) within bounds and well-formed.
    """
    if not ir.sequence:
        raise IRError(IRLevel.L1, "empty sequence")
    seq_upper = ir.sequence.upper()
    if not set(seq_upper) <= _VALID_RNA_BASES:
        bad = sorted(set(seq_upper) - _VALID_RNA_BASES)
        raise IRError(IRLevel.L1, f"invalid RNA bases: {bad}")

    seq_len = len(ir.sequence)
    for r in ir.regions:
        if r.start < 0 or r.end > seq_len or r.start >= r.end:
            raise IRError(IRLevel.L1, f"region out of bounds: {r}")
    return True


def check_l2_invariants(ir: IR_L2_MatureMRNA) -> bool:
    """L2 invariants.

    * CDS ends with a stop codon (start codon NOT required — back-translated genes may not start with M).
    * CDS ends with one of ``UAA``/``UAG``/``UGA`` (stop codon).
    * ``len(cds) % 3 == 0``.
    * ``sequence == five_utr + cds + three_utr`` (concatenation check).
    * All sub-sequences contain only RNA bases.
    """
    cds = ir.cds.upper()
    if not cds:
        raise IRError(IRLevel.L2, "empty CDS")
    # CDS does NOT need to start with AUG — back-translated genes may
    # start with any amino acid (e.g., antibody chains start with E = GAG).
    if cds[-3:] not in _RNA_STOP_CODONS:
        raise IRError(
            IRLevel.L2, f"CDS must end with stop codon: {cds[-3:]!r}"
        )
    if len(cds) % 3 != 0:
        raise IRError(
            IRLevel.L2, f"CDS length must be divisible by 3: {len(cds)}"
        )

    # Concatenation invariant.
    expected = ir.five_utr + ir.cds + ir.three_utr
    if ir.sequence != expected:
        raise IRError(
            IRLevel.L2,
            "sequence != five_utr + cds + three_utr",
        )

    # All sub-sequences must be valid RNA.
    for name, part in (
        ("five_utr", ir.five_utr),
        ("cds", ir.cds),
        ("three_utr", ir.three_utr),
    ):
        part_upper = part.upper()
        if not set(part_upper) <= _VALID_RNA_BASES:
            bad = sorted(set(part_upper) - _VALID_RNA_BASES)
            raise IRError(
                IRLevel.L2, f"invalid RNA bases in {name}: {bad}"
            )
    return True


def check_l3_invariants(ir: IR_L3_Polypeptide) -> bool:
    """L3 invariants.

    * Non-empty sequence.
    * Only valid single-letter amino acid codes (including ``*`` for
      stop and ``X`` for unknown).
    * Ends with ``*`` (stop).  Note: a polypeptide produced by
      :func:`biocompiler.ir.passes.translate` will always end with
      ``*`` because the CDS is required to end with a stop codon.
    """
    if not ir.sequence:
        raise IRError(IRLevel.L3, "empty polypeptide")
    seq_upper = ir.sequence.upper()
    invalid = set(seq_upper) - _VALID_AA
    if invalid:
        raise IRError(
            IRLevel.L3, f"invalid amino acids: {sorted(invalid)}"
        )
    if not seq_upper.endswith("*"):
        raise IRError(
            IRLevel.L3, "polypeptide must end with stop (*)"
        )
    return True


def check_l4_invariants(ir: IR_L4_FoldedProtein) -> bool:
    """L4 invariants.

    * Sequence (if non-empty) contains only valid amino acid codes.
    * Confidence, if present, is in ``[0, 100]`` (pLDDT range).
    * Coordinates, if present, must match the sequence length (the
      stub leaves this unchecked because the coordinate representation
      is not yet fixed — Phase 2 will pin this down).
    """
    if ir.sequence:
        seq_upper = ir.sequence.upper()
        invalid = set(seq_upper) - _VALID_AA
        if invalid:
            raise IRError(
                IRLevel.L4, f"invalid amino acids: {sorted(invalid)}"
            )

    if ir.confidence is not None:
        if not (0 <= ir.confidence <= 100):
            raise IRError(
                IRLevel.L4,
                f"confidence must be in [0, 100], got {ir.confidence}",
            )
    return True


# Convenience dispatcher — lets callers run the right checker based on
# the IR object's reported ``level`` property without a manual if/elif.
_CHECKERS = {
    IRLevel.L0: check_l0_invariants,
    IRLevel.L1: check_l1_invariants,
    IRLevel.L2: check_l2_invariants,
    IRLevel.L3: check_l3_invariants,
    IRLevel.L4: check_l4_invariants,
}


def check_invariants(ir: object) -> bool:
    """Dispatch to the appropriate level-specific checker.

    Uses the IR object's ``level`` property to pick the right checker.
    Raises :class:`IRError` (level L0 with a meta-message) if the
    object has no ``level`` attribute — this typically means a
    non-IR object was passed in.
    """
    level = getattr(ir, "level", None)
    if level not in _CHECKERS:
        raise IRError(
            IRLevel.L0,
            f"object has no recognised IR level: {type(ir).__name__}",
        )
    return _CHECKERS[level](ir)  # type: ignore[arg-type]


__all__ = [
    "check_l0_invariants",
    "check_l1_invariants",
    "check_l2_invariants",
    "check_l3_invariants",
    "check_l4_invariants",
    "check_invariants",
]
