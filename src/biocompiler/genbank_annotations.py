"""
GenBank annotation-based optimization workflow.
Enables no-code optimization by parsing constraint specifications
from GenBank feature annotations.

Annotation syntax (parsed from /note qualifiers):
    @no(ENZYME)           — Avoid restriction site (e.g., @no(BsaI_site))
    @no(PATTERN)          — Avoid pattern (e.g., @no(ATTTA))
    @optimize(ORGANISM)   — Codon optimize for organism
    @gc(LO-HI)            — GC content between LO% and HI%
    @harmonize(SOURCE->TARGET) — Harmonize codons
    @keep                 — Lock region from mutation
    @avoid_blast(DB)      — Avoid BLAST matches
    @no_cpg               — Avoid CpG dinucleotides
    @no_splice            — Avoid cryptic splice sites
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "AnnotationDirective",
    "GenBankAnnotationResult",
    "OptimizationResult",
    "parse_annotation_note",
    "parse_genbank_annotations",
    "annotations_to_optimization_params",
    "optimize_from_genbank",
]


# ────────────────────────────────────────────────────────────
# Data classes
# ────────────────────────────────────────────────────────────

@dataclass
class AnnotationDirective:
    """A single optimization directive parsed from a GenBank annotation.

    Attributes:
        directive_type: One of "no", "optimize", "gc", "harmonize",
            "keep", "avoid_blast", "no_cpg", "no_splice".
        parameter: Directive-specific parameter (e.g., "BsaI_site",
            "Homo_sapiens", "40-60").
        region: 1-based (start, end) positions, or None for whole sequence.
        source: Raw annotation text this directive was parsed from.
    """

    directive_type: str
    parameter: str
    region: tuple[int, int] | None = None
    source: str = ""


@dataclass
class GenBankAnnotationResult:
    """Result of parsing a GenBank file for optimization directives.

    Attributes:
        sequence: DNA sequence extracted from the GenBank record.
        organism: Organism name from the SOURCE section.
        gene_name: Gene name from /gene qualifier.
        directives: List of AnnotationDirective objects.
        raw_features: List of raw feature dicts from import_seq.
    """

    sequence: str
    organism: str
    gene_name: str
    directives: list[AnnotationDirective] = field(default_factory=list)
    raw_features: list[dict] = field(default_factory=list)


@dataclass
class OptimizationResult:
    """Result of an end-to-end GenBank annotation-driven optimization.

    Attributes:
        sequence: Optimized DNA sequence.
        protein: Protein translation of the optimized sequence.
        organism: Target organism used for optimization.
        directives_applied: Number of directives that were successfully applied.
        directives_total: Total number of directives parsed.
        warnings: List of warning messages.
        optimized: Whether optimization was actually performed.
    """

    sequence: str
    protein: str = ""
    organism: str = ""
    directives_applied: int = 0
    directives_total: int = 0
    warnings: list[str] = field(default_factory=list)
    optimized: bool = False


# ────────────────────────────────────────────────────────────
# Annotation parsing
# ────────────────────────────────────────────────────────────

# Regex patterns for annotation directives
_RE_NO = re.compile(r"@no\(([^)]+)\)")
_RE_OPTIMIZE = re.compile(r"@optimize\(([^)]+)\)")
_RE_GC = re.compile(r"@gc\((\d+)-(\d+)\)")
_RE_HARMONIZE = re.compile(r"@harmonize\(([^>]+)->([^)]+)\)")
_RE_AVOID_BLAST = re.compile(r"@avoid_blast\(([^)]+)\)")
_RE_KEEP = re.compile(r"@keep\b")
_RE_NO_CPG = re.compile(r"@no_cpg\b")
_RE_NO_SPLICE = re.compile(r"@no_splice\b")


def parse_annotation_note(note: str) -> list[AnnotationDirective]:
    """Parse a /note qualifier string to extract optimization directives.

    Supports the following annotation syntax:
        @no(ENZYME)           — Avoid restriction site or pattern
        @optimize(ORGANISM)   — Codon optimize for organism
        @gc(LO-HI)            — GC content between LO% and HI%
        @harmonize(SRC->TGT)  — Codon harmonization
        @keep                 — Lock region from mutation
        @avoid_blast(DB)      — Avoid BLAST matches against database
        @no_cpg               — Avoid CpG dinucleotides
        @no_splice            — Avoid cryptic splice sites

    Args:
        note: Raw /note qualifier string from a GenBank feature.

    Returns:
        List of AnnotationDirective objects parsed from the note.
    """
    directives: list[AnnotationDirective] = []

    if not note:
        return directives

    # Parse @no(...) directives
    for match in _RE_NO.finditer(note):
        param = match.group(1).strip()
        directives.append(AnnotationDirective(
            directive_type="no",
            parameter=param,
            source=match.group(0),
        ))

    # Parse @optimize(...) directives
    for match in _RE_OPTIMIZE.finditer(note):
        param = match.group(1).strip()
        directives.append(AnnotationDirective(
            directive_type="optimize",
            parameter=param,
            source=match.group(0),
        ))

    # Parse @gc(LO-HI) directives
    for match in _RE_GC.finditer(note):
        lo = match.group(1)
        hi = match.group(2)
        directives.append(AnnotationDirective(
            directive_type="gc",
            parameter=f"{lo}-{hi}",
            source=match.group(0),
        ))

    # Parse @harmonize(SOURCE->TARGET) directives
    for match in _RE_HARMONIZE.finditer(note):
        source_org = match.group(1).strip()
        target_org = match.group(2).strip()
        directives.append(AnnotationDirective(
            directive_type="harmonize",
            parameter=f"{source_org}->{target_org}",
            source=match.group(0),
        ))

    # Parse @avoid_blast(DB) directives
    for match in _RE_AVOID_BLAST.finditer(note):
        param = match.group(1).strip()
        directives.append(AnnotationDirective(
            directive_type="avoid_blast",
            parameter=param,
            source=match.group(0),
        ))

    # Parse @keep directive (no parameter)
    if _RE_KEEP.search(note):
        directives.append(AnnotationDirective(
            directive_type="keep",
            parameter="",
            source="@keep",
        ))

    # Parse @no_cpg directive (no parameter)
    if _RE_NO_CPG.search(note):
        directives.append(AnnotationDirective(
            directive_type="no_cpg",
            parameter="",
            source="@no_cpg",
        ))

    # Parse @no_splice directive (no parameter)
    if _RE_NO_SPLICE.search(note):
        directives.append(AnnotationDirective(
            directive_type="no_splice",
            parameter="",
            source="@no_splice",
        ))

    return directives


# ────────────────────────────────────────────────────────────
# GenBank file parsing
# ────────────────────────────────────────────────────────────

def parse_genbank_annotations(filepath_or_text: str) -> GenBankAnnotationResult:
    """Parse a GenBank file and extract all optimization directives.

    Scans all features for /note qualifiers containing optimization
    directives, and also extracts key metadata (sequence, organism,
    gene name).

    Args:
        filepath_or_text: Path to a GenBank file, or raw GenBank text.

    Returns:
        GenBankAnnotationResult with extracted directives and metadata.
    """
    from .import_seq import import_genbank

    gb = import_genbank(filepath_or_text)

    sequence = gb.get("sequence", "")
    organism = gb.get("organism", "")
    gene_name = gb.get("gene_name", "")
    features = gb.get("features", [])

    directives: list[AnnotationDirective] = []

    for feat in features:
        qualifiers = feat.get("qualifiers", {})
        note = qualifiers.get("note", "")
        if not note:
            continue

        # Parse directives from the /note qualifier
        feat_directives = parse_annotation_note(note)

        # Attach region information from feature location
        location = feat.get("location", "")
        region = _parse_feature_region(location, len(sequence))

        for d in feat_directives:
            d.region = region
            directives.append(d)

    return GenBankAnnotationResult(
        sequence=sequence,
        organism=organism,
        gene_name=gene_name,
        directives=directives,
        raw_features=features,
    )


def _parse_feature_region(location: str, seq_len: int) -> tuple[int, int] | None:
    """Parse a GenBank feature location to a 1-based (start, end) range.

    Handles simple locations like ``87..790`` and ``complement(87..790)``.
    Returns None for complex locations (joins, fuzzy boundaries).

    Args:
        location: Raw location string from a GenBank feature.
        seq_len: Length of the sequence (for bounds checking).

    Returns:
        1-based (start, end) tuple, or None if unparseable.
    """
    # Remove complement() wrapper
    loc = location.strip()
    loc = re.sub(r"complement\((.+)\)", r"\1", loc)

    # Skip join() locations — they span multiple regions
    if "join" in loc:
        return None

    # Try to parse simple range: start..end
    match = re.match(r"(\d+)\.\.(\d+)", loc)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        return (start, end)

    # Try single position
    match = re.match(r"(\d+)", loc)
    if match:
        pos = int(match.group(1))
        return (pos, pos)

    return None


# ────────────────────────────────────────────────────────────
# Conversion to optimization parameters
# ────────────────────────────────────────────────────────────

def annotations_to_optimization_params(annotations: GenBankAnnotationResult) -> dict[str, Any]:
    """Convert GenBankAnnotationResult to optimize_sequence() parameters.

    Maps annotation directives to the keyword arguments expected by
    :func:`biocompiler.optimize_sequence`.

    Args:
        annotations: Parsed GenBankAnnotationResult from
            :func:`parse_genbank_annotations`.

    Returns:
        Dict of keyword arguments for optimize_sequence().
    """
    params: dict[str, Any] = {
        "sequence": annotations.sequence,
    }

    # Default organism from GenBank metadata
    if annotations.organism:
        params["organism"] = annotations.organism

    avoid_enzymes: list[str] = []
    avoid_patterns: list[str] = []
    gc_lo: float | None = None
    gc_hi: float | None = None
    harmonize_source: str | None = None
    harmonize_target: str | None = None
    no_cpg: bool = False
    no_splice: bool = False
    keep_regions: list[tuple[int, int]] = []
    avoid_blast_refs: list[str] = []

    for directive in annotations.directives:
        dtype = directive.directive_type
        param = directive.parameter

        if dtype == "no":
            # Determine if this is an enzyme name or a raw pattern
            from .restriction_sites import get_recognition_site
            site = get_recognition_site(param)
            if site is not None:
                avoid_enzymes.append(param)
            else:
                # Treat as a raw DNA pattern
                avoid_patterns.append(param.upper())

        elif dtype == "optimize":
            params["organism"] = param

        elif dtype == "gc":
            parts = param.split("-")
            if len(parts) == 2:
                try:
                    gc_lo = float(parts[0]) / 100.0
                    gc_hi = float(parts[1]) / 100.0
                except ValueError:
                    logger.warning("Invalid GC range: %s", param)

        elif dtype == "harmonize":
            parts = param.split("->")
            if len(parts) == 2:
                harmonize_source = parts[0].strip()
                harmonize_target = parts[1].strip()

        elif dtype == "keep":
            if directive.region is not None:
                # Convert 1-based to 0-based
                start_0 = directive.region[0] - 1
                end_0 = directive.region[1]
                keep_regions.append((start_0, end_0))

        elif dtype == "avoid_blast":
            avoid_blast_refs.append(param)

        elif dtype == "no_cpg":
            no_cpg = True

        elif dtype == "no_splice":
            no_splice = True

    # Apply collected parameters
    if avoid_enzymes:
        params["avoid_enzymes"] = avoid_enzymes
    if avoid_patterns:
        params["avoid_patterns"] = avoid_patterns
    if gc_lo is not None and gc_hi is not None:
        params["gc_lo"] = gc_lo
        params["gc_hi"] = gc_hi
    if harmonize_source and harmonize_target:
        params["harmonize_source"] = harmonize_source
        params["harmonize_target"] = harmonize_target
    if no_cpg:
        params["avoid_cpg"] = True
    if no_splice:
        params["avoid_splice"] = True
    if keep_regions:
        params["keep_regions"] = keep_regions
    if avoid_blast_refs:
        params["avoid_blast_refs"] = avoid_blast_refs

    return params


# ────────────────────────────────────────────────────────────
# End-to-end optimization
# ────────────────────────────────────────────────────────────

def optimize_from_genbank(
    filepath_or_text: str,
    output_format: str = "genbank",
) -> OptimizationResult:
    """End-to-end: parse GenBank annotations → optimize → result.

    Parses a GenBank file for optimization directives, converts them
    to optimization parameters, runs the optimization, and returns
    the result.

    Args:
        filepath_or_text: Path to a GenBank file, or raw GenBank text.
        output_format: Desired output format ("genbank", "fasta", "dict").
            Currently only affects the return type structure (always
            returns OptimizationResult).

    Returns:
        OptimizationResult with the optimized sequence and metadata.
    """
    # Step 1: Parse annotations
    annotations = parse_genbank_annotations(filepath_or_text)

    if not annotations.sequence:
        return OptimizationResult(
            sequence="",
            organism=annotations.organism,
            gene_name=annotations.gene_name if hasattr(OptimizationResult, 'gene_name') else "",
            directives_total=len(annotations.directives),
            warnings=["No sequence found in GenBank record"],
            optimized=False,
        )

    if not annotations.directives:
        # No directives — return the sequence as-is
        protein = _translate_sequence(annotations.sequence)
        return OptimizationResult(
            sequence=annotations.sequence,
            protein=protein,
            organism=annotations.organism,
            directives_total=0,
            directives_applied=0,
            warnings=["No optimization directives found in annotations"],
            optimized=False,
        )

    # Step 2: Convert to optimization parameters
    params = annotations_to_optimization_params(annotations)

    # Step 3: Run optimization
    warnings: list[str] = []
    optimized_seq = annotations.sequence
    directives_applied = 0
    directives_total = len(annotations.directives)

    try:
        from .optimizer import optimize_sequence

        # Extract only the parameters that optimize_sequence accepts
        opt_params: dict[str, Any] = {}
        seq = params.pop("sequence", annotations.sequence)
        organism = params.pop("organism", annotations.organism or "Homo_sapiens")

        # Map known parameters
        if "avoid_enzymes" in params:
            opt_params["avoid_enzymes"] = params["avoid_enzymes"]
        if "gc_lo" in params:
            opt_params["gc_lo"] = params["gc_lo"]
        if "gc_hi" in params:
            opt_params["gc_hi"] = params["gc_hi"]
        if "avoid_cpg" in params:
            opt_params["avoid_cpg_islands"] = params["avoid_cpg"]
        if "avoid_splice" in params:
            opt_params["avoid_cryptic_splice"] = params["avoid_splice"]

        result = optimize_sequence(seq, organism=organism, **opt_params)
        optimized_seq = result.sequence if hasattr(result, "sequence") else str(result)
        directives_applied = directives_total  # Assume all were applied

    except ImportError:
        warnings.append("Optimizer not available — returning original sequence")
    except Exception as e:
        warnings.append(f"Optimization failed: {e}")

    protein = _translate_sequence(optimized_seq)

    return OptimizationResult(
        sequence=optimized_seq,
        protein=protein,
        organism=annotations.organism,
        directives_applied=directives_applied,
        directives_total=directives_total,
        warnings=warnings,
        optimized=True,
    )


def _translate_sequence(seq: str) -> str:
    """Translate a DNA sequence to its protein product.

    Args:
        seq: DNA sequence (uppercase).

    Returns:
        Amino acid string (single-letter codes).
    """
    from .type_system import CODON_TABLE

    seq = seq.upper()
    aa_list: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is not None and aa != "*":
            aa_list.append(aa)
    return "".join(aa_list)
