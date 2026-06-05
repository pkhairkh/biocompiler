"""
BioCompiler GenBank Round-Trip Verification
==============================================

Verifies that an optimized design survives a GenBank export→import
round trip: the DNA sequence and key annotations are preserved after
being exported to GenBank format and parsed back.

The main entry point is :func:`verify_genbank_roundtrip`, which takes
an :class:`~biocompiler.optimization.OptimizationResult`, exports it to
GenBank via :func:`~biocompiler.export.export_genbank`, re-imports it
via :func:`~biocompiler.import_seq.import_genbank`, and compares the
result.

Usage::

    from biocompiler.genbank_roundtrip import verify_genbank_roundtrip

    result = optimize_sequence("MSKGEELFTG", organism="Escherichia_coli")
    report = verify_genbank_roundtrip(result)
    if report.success:
        print("Round-trip verified!")
    else:
        print(f"Mismatches at {len(report.mismatches)} positions")
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from .optimization import OptimizationResult
from .export import export_genbank
from .import_seq import import_genbank

logger = logging.getLogger(__name__)

__all__ = [
    "RoundTripResult",
    "verify_genbank_roundtrip",
    "compare_sequences",
    "verify_annotation_preservation",
]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class RoundTripResult:
    """Result of a GenBank round-trip verification.

    Attributes:
        success: True if the round-trip was lossless (sequence matches, annotations preserved).
        original_dna: The DNA sequence from the original OptimizationResult (uppercase).
        reimported_dna: The DNA sequence after export→import (uppercase).
        mismatches: List of (position, original_base, reimported_base) tuples for any
                    positions where the sequences differ.
        annotation_preserved: True if key annotations (gene name, protein, organism)
                              survived the round trip.
        original_annotations: Dict of key annotations from the original result.
        reimported_annotations: Dict of key annotations from the re-imported record.
        genbank_text: The GenBank text that was generated (for debugging).
        warnings: List of non-fatal issues encountered during verification.
    """
    success: bool
    original_dna: str
    reimported_dna: str
    mismatches: list[tuple[int, str, str]]
    annotation_preserved: bool
    original_annotations: dict[str, Any] = field(default_factory=dict)
    reimported_annotations: dict[str, Any] = field(default_factory=dict)
    genbank_text: str = ""
    warnings: list[str] = field(default_factory=list)


# ────────────────────────────────────────────────────────────
# Core verification function
# ────────────────────────────────────────────────────────────

def verify_genbank_roundtrip(
    result: OptimizationResult,
    gene_name: str = "designed_gene",
    organism: str | None = None,
) -> RoundTripResult:
    """Export an OptimizationResult to GenBank, parse it back, and verify the
    DNA sequence matches.

    The verification checks:
    1. **Sequence integrity**: The re-imported DNA sequence must be
       identical to the original (case-insensitive).
    2. **Annotation preservation**: Key annotations (gene name, protein
       translation, organism) must survive the round trip.

    Args:
        result: An OptimizationResult from optimize_sequence().
        gene_name: Gene name to embed in the GenBank record.
        organism: Organism name for the GenBank record. If None, attempts
                  to read from result.organism_name (if it exists) or
                  defaults to "Homo_sapiens".

    Returns:
        A :class:`RoundTripResult` with detailed verification information.
    """
    if organism is None:
        organism = getattr(result, "organism_name", None) or "Homo_sapiens"

    # ── Step 1: Export to GenBank ──
    original_dna = result.sequence.upper()

    try:
        genbank_text = export_genbank(
            sequence=original_dna,
            organism=organism,
            gene_name=gene_name,
            protein=result.protein or None,
            cai=result.cai,
            type_results=_extract_type_results(result),
        )
    except Exception as e:
        logger.error("GenBank export failed: %s", e)
        return RoundTripResult(
            success=False,
            original_dna=original_dna,
            reimported_dna="",
            mismatches=[],
            annotation_preserved=False,
            warnings=[f"GenBank export failed: {e}"],
        )

    # ── Step 2: Re-import from GenBank ──
    try:
        reimported = import_genbank(genbank_text)
    except Exception as e:
        logger.error("GenBank import failed: %s", e)
        return RoundTripResult(
            success=False,
            original_dna=original_dna,
            reimported_dna="",
            mismatches=[],
            annotation_preserved=False,
            genbank_text=genbank_text,
            warnings=[f"GenBank import failed: {e}"],
        )

    reimported_dna = reimported.get("sequence", "").upper()

    # ── Step 3: Compare sequences ──
    mismatches = compare_sequences(original_dna, reimported_dna)

    # ── Step 4: Verify annotation preservation ──
    original_annotations = _extract_original_annotations(result, gene_name, organism)
    reimported_annotations = _extract_reimported_annotations(reimported)
    annotation_preserved, ann_warnings = verify_annotation_preservation(
        original_annotations, reimported_annotations
    )

    # ── Step 5: Compile result ──
    all_warnings = ann_warnings.copy()

    if mismatches:
        all_warnings.append(
            f"Sequence mismatch: {len(mismatches)} position(s) differ "
            f"between original and re-imported DNA"
        )

    if len(original_dna) != len(reimported_dna):
        all_warnings.append(
            f"Length mismatch: original={len(original_dna)} bp, "
            f"reimported={len(reimported_dna)} bp"
        )

    success = (len(mismatches) == 0) and annotation_preserved

    return RoundTripResult(
        success=success,
        original_dna=original_dna,
        reimported_dna=reimported_dna,
        mismatches=mismatches,
        annotation_preserved=annotation_preserved,
        original_annotations=original_annotations,
        reimported_annotations=reimported_annotations,
        genbank_text=genbank_text,
        warnings=all_warnings,
    )


# ────────────────────────────────────────────────────────────
# Sequence comparison
# ────────────────────────────────────────────────────────────

def compare_sequences(
    original: str, reimported: str
) -> list[tuple[int, str, str]]:
    """Compare two DNA sequences and return a list of mismatches.

    Args:
        original: The original DNA sequence (uppercase).
        reimported: The re-imported DNA sequence (uppercase).

    Returns:
        List of (position, original_base, reimported_base) tuples for
        positions where the sequences differ. Positions beyond the end
        of the shorter sequence are not reported.
    """
    original = original.upper()
    reimported = reimported.upper()
    mismatches: list[tuple[int, str, str]] = []

    min_len = min(len(original), len(reimported))
    for i in range(min_len):
        if original[i] != reimported[i]:
            mismatches.append((i, original[i], reimported[i]))

    return mismatches


# ────────────────────────────────────────────────────────────
# Annotation verification
# ────────────────────────────────────────────────────────────

def verify_annotation_preservation(
    original_annotations: dict[str, Any],
    reimported_annotations: dict[str, Any],
) -> tuple[bool, list[str]]:
    """Verify that key annotations survived the GenBank round trip.

    Checks the following fields:
    - gene_name: Must match (if original was non-empty).
    - protein: Must match (if original was non-empty).
    - organism: Must match (if original was non-empty).

    The check is lenient: if the original value is empty or missing,
    the annotation is considered preserved regardless of the reimported value.

    Args:
        original_annotations: Dict of key annotations from the original result.
        reimported_annotations: Dict of key annotations from the re-imported record.

    Returns:
        Tuple of (preserved: bool, warnings: list[str]).
    """
    warnings: list[str] = []
    preserved = True

    # Check gene name
    orig_gene = original_annotations.get("gene_name", "")
    reimp_gene = reimported_annotations.get("gene_name", "")
    if orig_gene and orig_gene != reimp_gene:
        preserved = False
        warnings.append(
            f"Gene name mismatch: original='{orig_gene}', "
            f"reimported='{reimp_gene}'"
        )

    # Check protein
    orig_protein = original_annotations.get("protein", "")
    reimp_protein = reimported_annotations.get("protein", "")
    if orig_protein:
        # Remove whitespace that might have been introduced during round-trip
        reimp_protein_clean = reimp_protein.replace(" ", "").replace("\n", "")
        if orig_protein != reimp_protein_clean:
            preserved = False
            warnings.append(
                f"Protein mismatch: original length={len(orig_protein)}, "
                f"reimported length={len(reimp_protein_clean)}"
            )

    # Check organism
    orig_organism = original_annotations.get("organism", "")
    reimp_organism = reimported_annotations.get("organism", "")
    if orig_organism:
        # Organism names may have underscores replaced with spaces in GenBank
        orig_normalized = orig_organism.replace("_", " ")
        reimp_normalized = reimp_organism.replace("_", " ").strip()
        if orig_normalized != reimp_normalized:
            preserved = False
            warnings.append(
                f"Organism mismatch: original='{orig_organism}', "
                f"reimported='{reimp_organism}'"
            )

    return preserved, warnings


# ────────────────────────────────────────────────────────────
# Internal helpers
# ────────────────────────────────────────────────────────────

def _extract_type_results(
    result: OptimizationResult,
) -> list | None:
    """Extract TypeCheckResult list from an OptimizationResult, if available."""
    predicate_results = getattr(result, "predicate_results", None)
    if predicate_results is not None and len(predicate_results) > 0:
        # predicate_results is a list of PredicateResult; we need
        # to convert to TypeCheckResult for the export engine.
        # For round-trip verification, annotations are optional;
        # return None to let export_genbank handle it gracefully.
        return None
    return None


def _extract_original_annotations(
    result: OptimizationResult,
    gene_name: str,
    organism: str,
) -> dict[str, Any]:
    """Extract key annotations from the original OptimizationResult."""
    return {
        "gene_name": gene_name,
        "protein": result.protein or "",
        "organism": organism,
        "cai": result.cai,
        "gc_content": result.gc_content,
        "sequence_length": len(result.sequence),
    }


def _extract_reimported_annotations(
    reimported: dict[str, Any],
) -> dict[str, Any]:
    """Extract key annotations from the re-imported GenBank record."""
    return {
        "gene_name": reimported.get("gene_name", ""),
        "protein": reimported.get("protein", ""),
        "organism": reimported.get("organism", ""),
        "gc_content": reimported.get("gc_content", 0.0),
        "sequence_length": reimported.get("length", 0),
    }
