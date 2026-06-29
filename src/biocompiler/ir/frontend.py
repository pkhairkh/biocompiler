"""
BioCompiler IR Frontend — YAML Gene Specification Parser
========================================================

Users write gene designs in YAML:

    gene_name: HBB
    organism: human
    sequence: ATGGTGCATCTG...TAA
    regions:
      - type: cds
        start: 0
        end: 162

The frontend parses this into an :class:`IR_L0_GenomicDNA` object, which
can then be compiled through the IR pipeline (L0→L1→L2→L3→L4) via
:func:`biocompiler.ir.passes.compile_gene`.

This is the compiler's *source language*: a gene spec written in YAML is
to BioCompiler what a ``.c`` file is to gcc — the human-readable input
that the frontend turns into the IR that the rest of the compiler
operates on.

Usage
-----
::

    from biocompiler.ir.frontend import parse_spec, compile_from_spec
    from biocompiler.ir.types import IRLevel

    # Parse only — get the IR-L0 object:
    ir_l0 = parse_spec("gene.yaml")
    ir_l0 = parse_spec("gene_name: test\\norganism: e_coli\\nsequence: ATGGCTTAA\\n")

    # Parse + compile through the IR pipeline:
    ir_l3 = compile_from_spec("gene.yaml")                       # default: L3
    ir_l4 = compile_from_spec("gene.yaml", target_level=IRLevel.L4)

YAML schema
-----------
Required fields::

    sequence : str    — DNA sequence (A, C, G, T, N; case-insensitive)
    organism : str    — source organism (free-form, e.g. "human", "e_coli")

Optional fields::

    gene_name : str   — gene symbol (e.g. "HBB")
    regions   : list  — list of region annotations; each entry has:
        type          : str  — one of "exon", "intron", "5_utr", "3_utr",
                                 "promoter", "terminator", "cds"
                                 (other strings tolerated as opaque annotations)
        start         : int  — 0-based inclusive start
        end           : int  — 0-based exclusive end
        metadata      : dict — free-form per-region metadata (optional)
    metadata  : dict  — free-form top-level metadata

Coordinates follow the same 0-indexed half-open ``[start, end)``
convention as :class:`biocompiler.ir.types.GeneRegion` (and Python
slicing).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Union

import yaml

from .types import IRLevel, GeneRegion, IR_L0_GenomicDNA, IRError
from .passes import compile_gene, IRObject

# Set of valid DNA bases (uppercase).  Lowercase input is normalised by
# :func:`_validate_sequence`.  ``N`` is allowed (matches the IR-L0
# invariant — see ``biocompiler.ir.invariants.check_l0_invariants``).
_VALID_DNA_BASES = frozenset("ACGTN")

# Region types understood by the IR lowering passes (see the
# ``GeneRegion`` docstring).  Unknown strings are tolerated — they are
# treated as opaque annotations and ignored by ``splice`` — but listing
# them here lets us emit a helpful hint in error messages.
_CANONICAL_REGION_TYPES = frozenset({
    "exon", "intron", "5_utr", "3_utr",
    "promoter", "terminator", "cds",
})


class SpecError(ValueError):
    """Raised when a YAML gene spec is malformed.

    Subclasses ``ValueError`` so callers that catch ``ValueError`` (e.g.
    argparse type handlers) still work.  Kept distinct from
    :class:`biocompiler.ir.types.IRError`, which is raised by IR
    invariant checks *after* parsing — separating the two lets users
    distinguish "my YAML was wrong" from "the IR rejected the gene".
    """


# ────────────────────────────────────────────────────────────────────
# YAML loading
# ────────────────────────────────────────────────────────────────────
def _looks_like_path(s: str) -> bool:
    """Heuristic: does ``s`` look like a filesystem path to a YAML file?

    A path is assumed to be short, contain no newlines, and end in
    ``.yaml`` or ``.yml``.  This avoids mis-classifying a multi-line
    YAML document as a path (a YAML string almost always contains a
    newline and never ends in ``.yaml``).
    """
    if "\n" in s or "\r" in s:
        return False
    return s.endswith(".yaml") or s.endswith(".yml")


def _load_yaml(yaml_path_or_string: str) -> Any:
    """Load a YAML document from a file path or a YAML string.

    The input is treated as a file path iff it *looks* like one (see
    :func:`_looks_like_path`) AND that file exists on disk.  Otherwise
    it is parsed as a YAML string.  This dual-mode behaviour lets
    callers pass either ``"path/to/gene.yaml"`` or a literal YAML
    document without having to choose between two entry points.
    """
    if _looks_like_path(yaml_path_or_string):
        path = Path(yaml_path_or_string)
        if path.is_file():
            with path.open("r", encoding="utf-8") as fh:
                return yaml.safe_load(fh)
        # Looks like a path but doesn't exist — fall through and let
        # ``yaml.safe_load`` raise a clear parse error on the string.
    return yaml.safe_load(yaml_path_or_string)


# ────────────────────────────────────────────────────────────────────
# Field validators
# ────────────────────────────────────────────────────────────────────
def _validate_sequence(seq: Any) -> str:
    """Normalise and validate a DNA sequence.

    Returns the uppercased sequence.  Raises :class:`SpecError` if the
    sequence is empty, not a string, or contains characters outside
    ``{A, C, G, T, N}``.
    """
    if not isinstance(seq, str):
        raise SpecError(
            f"'sequence' must be a string, got {type(seq).__name__}"
        )
    if not seq:
        raise SpecError("'sequence' is required and must be non-empty")
    seq_upper = seq.upper()
    bad = set(seq_upper) - _VALID_DNA_BASES
    if bad:
        raise SpecError(
            f"invalid bases in sequence: {sorted(bad)} "
            f"(allowed: A, C, G, T, N)"
        )
    return seq_upper


def _validate_organism(organism: Any) -> str:
    """Validate the ``organism`` field."""
    if not isinstance(organism, str) or not organism.strip():
        raise SpecError(
            f"'organism' must be a non-empty string, got {organism!r}"
        )
    return organism.strip()


def _validate_gene_name(gene_name: Any) -> Union[str, None]:
    """Validate the optional ``gene_name`` field."""
    if gene_name is None:
        return None
    if not isinstance(gene_name, str):
        raise SpecError(
            f"'gene_name' must be a string, got {type(gene_name).__name__}"
        )
    return gene_name


def _validate_region(raw: Any, seq_len: int) -> GeneRegion:
    """Convert a raw YAML region dict into a :class:`GeneRegion`.

    Performs structural validation: required keys (``type``, ``start``,
    ``end``), integer types, coordinate ordering, and bounds.  Unknown
    region types are tolerated (treated as opaque annotations) to match
    the permissiveness of :class:`GeneRegion`.
    """
    if not isinstance(raw, dict):
        raise SpecError(
            f"region entry must be a mapping, got {type(raw).__name__}: {raw!r}"
        )

    # --- type ---
    rtype = raw.get("type")
    if rtype is None:
        raise SpecError(f"region is missing required 'type' field: {raw!r}")
    if not isinstance(rtype, str):
        raise SpecError(
            f"region 'type' must be a string, got {type(rtype).__name__}"
        )
    # Note: we deliberately do NOT reject unknown region types — the IR
    # dataclass treats them as opaque annotations (see GeneRegion
    # docstring).  We only check the canonical types when emitting
    # hint-bearing error messages elsewhere.

    # --- start / end ---
    start = raw.get("start")
    end = raw.get("end")
    if start is None:
        raise SpecError(f"region is missing required 'start' field: {raw!r}")
    if end is None:
        raise SpecError(f"region is missing required 'end' field: {raw!r}")
    # ``bool`` is a subclass of ``int`` in Python — guard explicitly so
    # that ``start: true`` doesn't silently become ``start: 1``.
    if not isinstance(start, int) or isinstance(start, bool):
        raise SpecError(
            f"region 'start' must be an integer, got {type(start).__name__}: {start!r}"
        )
    if not isinstance(end, int) or isinstance(end, bool):
        raise SpecError(
            f"region 'end' must be an integer, got {type(end).__name__}: {end!r}"
        )
    if end <= start:
        raise SpecError(
            f"region 'end' ({end}) must be > 'start' ({start}): {raw!r}"
        )
    if start < 0 or end > seq_len:
        raise SpecError(
            f"region [{start}, {end}) is out of bounds for sequence of "
            f"length {seq_len}: {raw!r}"
        )

    # --- metadata (optional) ---
    meta = raw.get("metadata", {})
    if meta is None:
        meta = {}
    if not isinstance(meta, dict):
        raise SpecError(
            f"region 'metadata' must be a mapping, got {type(meta).__name__}"
        )

    return GeneRegion(
        start=start, end=end, region_type=rtype, metadata=dict(meta),
    )


def _validate_regions(raw_regions: Any, seq_len: int) -> list[GeneRegion]:
    """Validate the optional ``regions`` list."""
    if raw_regions is None:
        return []
    if not isinstance(raw_regions, list):
        raise SpecError(
            f"'regions' must be a list, got {type(raw_regions).__name__}"
        )
    return [_validate_region(r, seq_len) for r in raw_regions]


def _validate_metadata(raw_meta: Any) -> dict:
    """Validate the optional top-level ``metadata`` mapping."""
    if raw_meta is None:
        return {}
    if not isinstance(raw_meta, dict):
        raise SpecError(
            f"'metadata' must be a mapping, got {type(raw_meta).__name__}"
        )
    return dict(raw_meta)


# ────────────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────────────
def parse_spec(yaml_path_or_string: str) -> IR_L0_GenomicDNA:
    """Parse a YAML gene specification into an :class:`IR_L0_GenomicDNA`.

    Parameters
    ----------
    yaml_path_or_string : str
        Either a path to a YAML file (must end in ``.yaml`` or ``.yml``)
        or a literal YAML string.

    Returns
    -------
    IR_L0_GenomicDNA
        The IR-L0 object ready to be fed into
        :func:`biocompiler.ir.passes.compile_gene`.

    Raises
    ------
    SpecError
        If the spec is malformed (missing required fields, invalid
        bases, out-of-bounds regions, ...).

    Notes
    -----
    The returned IR-L0 carries provenance in its ``metadata``:

    * ``metadata["pass"] == "frontend.parse_spec"``  — names the producer.
    * ``metadata["source_format"] == "yaml"``        — names the source
      language, so downstream passes know where the IR-L0 came from.

    Any user-supplied ``metadata`` in the YAML is preserved; the two
    provenance keys above are added only if not already present (so a
    user can override them if they really want to).
    """
    if not isinstance(yaml_path_or_string, str):
        raise SpecError(
            f"input must be a string (path or YAML), got "
            f"{type(yaml_path_or_string).__name__}"
        )
    if not yaml_path_or_string.strip():
        raise SpecError("input is empty")

    spec = _load_yaml(yaml_path_or_string)

    if spec is None:
        raise SpecError("YAML document is empty")
    if not isinstance(spec, dict):
        raise SpecError(
            f"YAML root must be a mapping, got {type(spec).__name__}"
        )

    # --- required fields ---
    if "sequence" not in spec:
        raise SpecError("spec is missing required field: 'sequence'")
    if "organism" not in spec:
        raise SpecError("spec is missing required field: 'organism'")

    sequence = _validate_sequence(spec["sequence"])
    organism = _validate_organism(spec["organism"])
    gene_name = _validate_gene_name(spec.get("gene_name"))

    # --- optional fields ---
    regions = _validate_regions(spec.get("regions"), len(sequence))
    metadata = _validate_metadata(spec.get("metadata"))

    # Stamp provenance: this IR-L0 was produced by the YAML frontend.
    # ``setdefault`` so a user can override these keys if they want.
    metadata.setdefault("pass", "frontend.parse_spec")
    metadata.setdefault("source_format", "yaml")

    return IR_L0_GenomicDNA(
        sequence=sequence,
        regions=regions,
        organism=organism,
        gene_name=gene_name,
        metadata=metadata,
    )


def compile_from_spec(
    yaml_path_or_string: str,
    target_level: IRLevel = IRLevel.L3,
) -> IRObject:
    """Parse a YAML gene spec and compile it through the IR pipeline.

    Convenience wrapper combining :func:`parse_spec` and
    :func:`biocompiler.ir.passes.compile_gene`::

        ir_l0 = parse_spec(yaml_path_or_string)
        return compile_gene(ir_l0, target_level)

    Parameters
    ----------
    yaml_path_or_string : str
        Path to a YAML file, or a literal YAML string.
    target_level : IRLevel, optional
        IR level to stop at.  Defaults to :attr:`IRLevel.L3`
        (polypeptide), the deepest level with a real lowering pass in
        Phase 1.

    Returns
    -------
    IRObject
        An IR object of the requested level (L0 through L4).

    Raises
    ------
    SpecError
        If the YAML spec is malformed (raised by :func:`parse_spec`).
    IRError
        If any lowering pass rejects its input (raised by
        :func:`compile_gene`).
    """
    ir_l0 = parse_spec(yaml_path_or_string)
    return compile_gene(ir_l0, target_level)


__all__ = [
    "parse_spec",
    "compile_from_spec",
    "SpecError",
]
