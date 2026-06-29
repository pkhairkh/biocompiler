"""
BioCompiler IR — Core Data Structures
=====================================

Defines the five intermediate representation (IR) levels of BioCompiler,
mirroring the Central Dogma of molecular biology:

    IR-L0  GenomicDNA      →  raw DNA sequence with region annotations
    IR-L1  PreMRNA         →  transcribed RNA (T → U), regions preserved
    IR-L2  MatureMRNA      →  spliced mRNA: 5'UTR + CDS + 3'UTR
    IR-L3  Polypeptide     →  translated amino acid sequence
    IR-L4  FoldedProtein   →  3D structure (ESMFold/AlphaFold oracle — Phase 2)

Each level is a frozen-ish dataclass: fields are populated by a lowering
pass from the previous level (see :mod:`biocompiler.ir.passes`), and
every level exposes a ``level`` property returning the corresponding
:class:`IRLevel` enum member.

These dataclasses are the typed objects that flow through the compiler
pipeline.  They are intentionally minimal — semantic checks live in
:mod:`biocompiler.ir.invariants`, not in the dataclasses themselves —
so that a malformed IR can still be constructed and then rejected by an
explicit invariant check (which is more useful for diagnostics).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional


class IRLevel(Enum):
    """The five IR levels, ordered L0 → L4 (matches the Central Dogma)."""

    L0 = "genomic_dna"
    L1 = "pre_mrna"
    L2 = "mature_mrna"
    L3 = "polypeptide"
    L4 = "folded_protein"

    @property
    def order(self) -> int:
        """Numerical order of the level (0..4), independent of ``.value``.

        The string ``.value`` of an :class:`IRLevel` member is a
        human-readable label that is **not** alphabetically ordered, so
        comparing ``.value`` cannot be used to determine the pipeline
        stage.  Use this property (or compare enum members directly via
        :data:`IRLevel._ORDER`) instead.
        """
        return _LEVEL_ORDER.index(self)


# Canonical ordering of IR levels — used by ``IRLevel.order`` and by
# ``compile_gene`` in passes.py.  Defined explicitly (rather than via
# ``list(IRLevel)``) so that re-ordering the enum members above does not
# silently change pipeline semantics.
_LEVEL_ORDER: list[IRLevel] = [
    IRLevel.L0,
    IRLevel.L1,
    IRLevel.L2,
    IRLevel.L3,
    IRLevel.L4,
]


@dataclass
class GeneRegion:
    """An annotated region of a gene (exon, intron, UTR, promoter, ...).

    Coordinates are **0-indexed half-open** intervals: ``[start, end)``
    on the underlying sequence.  This is the same convention used by
    Python slicing and by the BED / GFF3 (1-based-start) formats after
    subtracting 1 from the start.

    Attributes
    ----------
    start:
        Inclusive start index (0-based).
    end:
        Exclusive end index (0-based).  Must satisfy ``end > start``.
    region_type:
        One of ``"exon"``, ``"intron"``, ``"5_utr"``, ``"3_utr"``,
        ``"promoter"``, ``"terminator"``, ``"cds"``.  Other strings are
        tolerated (treated as opaque annotations) but lowering passes
        only act on the canonical types.
    metadata:
        Free-form per-region metadata (e.g. ``{"frame": 0}``).
    """

    start: int
    end: int
    region_type: str
    metadata: dict = field(default_factory=dict)


@dataclass
class IR_L0_GenomicDNA:
    """IR-L0: Genomic DNA — the source-level representation.

    This is what a BioCompiler "frontend" produces: a raw DNA sequence
    plus region annotations (promoter, 5'UTR, exons, introns, 3'UTR,
    terminator).  For prokaryotic genes the region list is typically
    empty and the entire sequence is treated as the CDS.
    """

    sequence: str  # raw DNA sequence (A, C, G, T, N)
    regions: list[GeneRegion]
    organism: str
    gene_name: Optional[str] = None
    secis_positions: list[int] = field(default_factory=list)  # codon indices where TGA→U
    metadata: dict = field(default_factory=dict)

    @property
    def level(self) -> IRLevel:
        return IRLevel.L0


@dataclass
class IR_L1_PreMRNA:
    """IR-L1: Pre-mRNA — transcribed sequence.

    Identical to IR-L0 at the sequence level except that every ``T`` is
    replaced by ``U``.  Region annotations (exons, introns, UTRs) are
    preserved with the same coordinates, since transcription does not
    alter positions.
    """

    sequence: str  # RNA sequence (A, C, G, U)
    regions: list[GeneRegion]
    organism: str
    gene_name: Optional[str] = None
    secis_positions: list[int] = field(default_factory=list)  # codon indices where UGA→U
    metadata: dict = field(default_factory=dict)

    @property
    def level(self) -> IRLevel:
        return IRLevel.L1


@dataclass
class IR_L2_MatureMRNA:
    """IR-L2: Mature mRNA — spliced transcript.

    Introns have been removed and the three functional sub-regions of
    the mature transcript — 5'UTR, CDS, 3'UTR — are stored as separate
    strings.  The full ``sequence`` field is the concatenation
    ``five_utr + cds + three_utr``.

    Invariants enforced by :func:`biocompiler.ir.invariants.check_l2_invariants`:

    * ``cds`` starts with ``"AUG"`` (start codon).
    * ``cds`` ends with one of ``"UAA"``, ``"UAG"``, ``"UGA"`` (stop).
    * ``len(cds) % 3 == 0``.
    """

    sequence: str  # spliced RNA sequence = five_utr + cds + three_utr
    five_utr: str
    cds: str
    three_utr: str
    organism: str
    gene_name: Optional[str] = None
    secis_positions: list[int] = field(default_factory=list)  # CDS codon indices where UGA→U
    metadata: dict = field(default_factory=dict)

    @property
    def level(self) -> IRLevel:
        return IRLevel.L2


@dataclass
class IR_L3_Polypeptide:
    """IR-L3: Polypeptide — translated amino acid sequence.

    Single-letter amino acid codes with ``"*"`` for the stop codon.
    The sequence always ends with ``"*"`` if the CDS contained a
    proper stop codon (see IR-L2 invariants).
    """

    sequence: str
    organism: str
    gene_name: Optional[str] = None
    metadata: dict = field(default_factory=dict)

    @property
    def level(self) -> IRLevel:
        return IRLevel.L3


@dataclass
class IR_L4_FoldedProtein:
    """IR-L4: Folded protein — 3D structure + PTMs (oracle-dependent).

    Phase 1 (this code) leaves ``coordinates`` and ``confidence`` as
    ``None`` — folding is not yet implemented.  Phase 2 will populate
    these fields by calling ESMFold or AlphaFold as an oracle, with the
    resulting per-residue pLDDT stored in ``confidence``.

    Post-translational modifications (PTMs) are stored as a list of
    dicts, e.g. ``[{"position": 12, "mod": "phosphorylation", "residue": "S"}]``.
    """

    sequence: str
    coordinates: Optional[Any] = None  # 3D coords (from ESMFold/AlphaFold)
    confidence: Optional[float] = None  # pLDDT or similar, in [0, 100]
    ptms: list[dict] = field(default_factory=list)
    organism: str = ""
    gene_name: Optional[str] = None
    metadata: dict = field(default_factory=dict)

    @property
    def level(self) -> IRLevel:
        return IRLevel.L4


@dataclass
class IRError(Exception):
    """Raised when an IR invariant is violated.

    Carries the :class:`IRLevel` of the offending object so callers can
    report *where* in the pipeline the failure occurred.  Defined as a
    dataclass for uniform construction (``IRError(IRLevel.L2, "msg")``)
    while still being catchable as a normal ``Exception``.

    Note
    ----
    The dataclass-generated ``__init__`` does not call
    ``Exception.__init__``, so ``self.args`` is empty — but ``__str__``
    is overridden below, so ``str(error)`` and ``repr(error)`` behave
    intuitively.
    """

    level: IRLevel
    message: str

    def __str__(self) -> str:  # pragma: no cover - trivial formatter
        return f"IR-{self.level.name} Error: {self.message}"


__all__ = [
    "IRLevel",
    "GeneRegion",
    "IR_L0_GenomicDNA",
    "IR_L1_PreMRNA",
    "IR_L2_MatureMRNA",
    "IR_L3_Polypeptide",
    "IR_L4_FoldedProtein",
    "IRError",
]
