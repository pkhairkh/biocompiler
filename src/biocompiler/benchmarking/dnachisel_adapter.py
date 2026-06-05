"""
BioCompiler Benchmarking — DNAchisel Adapter
=============================================

Wraps DNAchisel's constraint-solving API to produce results directly comparable
with BioCompiler, enabling fair head-to-head benchmarking.

This adapter translates BioCompiler's constraint model into DNAchisel
specifications, runs the optimization, and converts the results back into
a common ``OptimizationResult`` format.  Because DNAchisel lacks a native
CAI optimization constraint, the adapter seeds the initial sequence with
BioCompiler's highest-CAI codons and computes CAI post-hoc using BioCompiler's
own evaluator — ensuring metric consistency across tools.

The adapter gracefully handles the case where DNAchisel is not installed,
raising ``ImportError`` with install instructions only when the adapter is
actually used (not at import time for the rest of the benchmarking package).

Usage::

    from biocompiler.benchmarking.dnachisel_adapter import DNAchiselAdapter, is_dnachisel_available

    if is_dnachisel_available():
        adapter = DNAchiselAdapter()
        result = adapter.optimize(
            protein="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
            organism="Homo_sapiens",
            constraints=[
                {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                {"type": "avoid_restriction", "enzymes": ["EcoRI", "BamHI"]},
            ],
        )
        print(f"CAI={result.cai:.4f}, GC={result.gc_content:.3f}, "
              f"Sites={result.restriction_site_count}, "
              f"Time={result.execution_time_s:.2f}s")
    else:
        print("DNAchisel not installed — pip install dnachisel")

DNAchisel Constraint Mapping:
    biocompiler constraint       -> DNAchisel specification
    ----------------------------    ---------------------------
    GC range (gc_lo, gc_hi)      -> EnforceGCContent(mini, maxi, window)
    Restriction site avoidance   -> AvoidPattern(site)
    Amino acid identity          -> EnforceTranslation(translation=protein)
    Codon optimization           -> CodonOptimize(species=species) [objective]
                                    (CAI recomputed with compute_cai_validated)

Note: DNAchisel's CodonOptimize objective is used to drive codon selection,
but CAI is always recomputed post-hoc using ``compute_cai_validated`` from
the metrics module — DNAchisel's own CAI output is NOT trusted. This ensures
fair comparison with BioCompiler using a single, consistent CAI evaluator.

References:
  Zulkower, V., Rosas, A., & Pujos, P. (2020). "DNA Chisel: a versatile
  sequence optimizer." *Bioinformatics* 36(16):4512–4519.
"""

import logging
import time
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "OptimizationResult",
    "DNAchiselAdapter",
    "is_dnachisel_available",
]

# ─── DNAchisel availability check ────────────────────────────────────

_DNACHISEL_AVAILABLE: bool = False
_DNACHISEL_ERROR: str = ""

try:
    from dnachisel import (
        DnaOptimizationProblem,
        AvoidPattern,
        CodonOptimize,
        EnforceGCContent,
        EnforceTranslation,
    )
    _DNACHISEL_AVAILABLE = True
    logger.debug("DNAchisel is available for benchmarking")
except ImportError as exc:
    _DNACHISEL_ERROR = str(exc)
    logger.debug(
        "DNAchisel not installed — benchmarking adapter will raise on use. "
        "Install with: pip install dnachisel"
    )


def is_dnachisel_available() -> bool:
    """Return True if DNAchisel is installed and importable."""
    return _DNACHISEL_AVAILABLE


# ─── Named constants ─────────────────────────────────────────────────

GC_ENFORCEMENT_WINDOW = 50
"""Sliding window size (bp) for local GC content enforcement in DNAchisel."""

MAX_RESTRICTION_ENZYMES = 10
"""Maximum number of restriction enzymes to pass to DNAchisel."""

# ─── Optimization result ─────────────────────────────────────────────


@dataclass
class OptimizationResult:
    """Result of DNAchisel optimization, comparable with BioCompiler output.

    Attributes:
        sequence: Optimized DNA sequence.
        protein: Target protein sequence (single-letter codes).
        cai: Codon Adaptation Index computed post-hoc.
        gc_content: GC content fraction of the optimized sequence.
        restriction_site_count: Number of restriction enzyme sites found.
        execution_time_s: Wall-clock time for optimization (seconds).
        success: Whether optimization completed without error.
        error: Error message if optimization failed, else None.
        constraints_applied: List of DNAchisel specification class names used.
    """
    sequence: str
    protein: str
    cai: float
    gc_content: float
    restriction_site_count: int
    execution_time_s: float
    success: bool
    error: str | None = None
    constraints_applied: list[str] = field(default_factory=list)


# ─── Amino acid / organism helpers ───────────────────────────────────


def _build_initial_sequence(protein: str, organism: str = "Homo_sapiens") -> str:
    """Build an initial DNA sequence from a protein using highest-CAI codons.

    DNAchisel requires an initial sequence to optimize from. We seed it
    with BioCompiler's best-codon-per-position sequence, giving DNAchisel
    a strong starting point.

    Args:
        protein: Amino acid sequence (single-letter codes).
        organism: Target organism for codon usage.

    Returns:
        DNA sequence using preferred codons for each amino acid.

    Raises:
        ImportError: If biocompiler internals are unavailable.
        ValueError: If an unknown amino acid is encountered.
    """
    from ..organisms import CODON_ADAPTIVENESS_TABLES, SUPPORTED_ORGANISMS
    from ..constants import AA_TO_CODONS

    if organism not in SUPPORTED_ORGANISMS:
        organism = "Homo_sapiens"

    usage = CODON_ADAPTIVENESS_TABLES.get(
        organism, CODON_ADAPTIVENESS_TABLES["Homo_sapiens"]
    )

    sorted_codons: dict[str, list[str]] = {}
    for aa in set(protein):
        codons = AA_TO_CODONS.get(aa, [])
        codons_sorted = sorted(
            codons, key=lambda c: usage.get(c, 0.0), reverse=True
        )
        sorted_codons[aa] = codons_sorted

    sequence_chars: list[str] = []
    for aa in protein:
        codons = sorted_codons.get(aa, [])
        if codons:
            sequence_chars.append(codons[0])
        else:
            raise ValueError(
                f"Unknown amino acid '{aa}' at position "
                f"{len(sequence_chars) + 1}. Only standard single-letter "
                f"IUPAC codes are supported."
            )

    return "".join(sequence_chars)


# Mapping from BioCompiler organism keys to DNAchisel species names
_DNACHISEL_SPECIES_MAP: dict[str, str] = {
    "Escherichia_coli": "e_coli",
    "Homo_sapiens": "h_sapiens",
    "Saccharomyces_cerevisiae": "s_cerevisiae",
    "Mus_musculus": "m_musculus",
}


def _build_dnachisel_spec(
    protein: str,
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    restriction_enzymes: list[str] | None = None,
) -> tuple[list[Any], list[Any]]:
    """Convert biocompiler constraints to DNAchisel specification format.

    Creates lists of DNAchisel constraint and objective objects that
    correspond to biocompiler's type predicates:

    Constraints:
    - EnforceTranslation: ensure the sequence encodes the target protein
    - EnforceGCContent: keep GC content within bounds
    - AvoidPattern: avoid restriction enzyme recognition sites

    Objectives:
    - CodonOptimize: maximize codon adaptation for the target species

    CAI is also evaluated post-hoc using ``compute_cai_validated`` to
    ensure metric consistency with BioCompiler.

    Args:
        protein: Target protein sequence (single-letter codes).
        organism: Target organism for codon optimization.
        gc_lo: Minimum GC content fraction.
        gc_hi: Maximum GC content fraction.
        restriction_enzymes: List of enzyme names to avoid.

    Returns:
        Tuple of (constraints, objectives) — each a list of DNAchisel
        Specification objects.

    Raises:
        ImportError: If DNAchisel is not installed.
    """
    if not _DNACHISEL_AVAILABLE:
        raise ImportError(
            f"DNAchisel is required for specification building: "
            f"{_DNACHISEL_ERROR}. "
            "Install with: pip install dnachisel"
        )

    constraints: list[Any] = []
    objectives: list[Any] = []

    # Enforce that the sequence translates to the target protein
    constraints.append(EnforceTranslation(translation=protein))

    # Enforce GC content range
    constraints.append(
        EnforceGCContent(
            mini=gc_lo,
            maxi=gc_hi,
            window=GC_ENFORCEMENT_WINDOW,
        )
    )

    # Avoid restriction enzyme sites
    if restriction_enzymes:
        from ..constants import RESTRICTION_ENZYMES

        for enz_name in restriction_enzymes:
            site = RESTRICTION_ENZYMES.get(enz_name)
            if site:
                # Skip IUPAC-ambiguous sites (DNAchisel handles these
                # differently)
                if any(b not in "ACGT" for b in site.upper()):
                    continue
                constraints.append(AvoidPattern(site))

    # Add CodonOptimize objective for the target species
    species = _DNACHISEL_SPECIES_MAP.get(organism)
    if species:
        objectives.append(CodonOptimize(species=species))
        logger.debug("Added CodonOptimize objective for species='%s'", species)
    else:
        logger.debug(
            "No DNAchisel species mapping for organism='%s'; "
            "skipping CodonOptimize objective",
            organism,
        )

    return constraints, objectives


def _count_restriction_sites(
    sequence: str,
    restriction_enzymes: list[str] | None = None,
) -> int:
    """Count restriction enzyme recognition sites in a sequence.

    Checks both forward and reverse complement strands.

    Args:
        sequence: DNA sequence to check.
        restriction_enzymes: List of enzyme names to check for.

    Returns:
        Total count of restriction sites found.
    """
    from ..constants import RESTRICTION_ENZYMES, reverse_complement

    if not restriction_enzymes:
        restriction_enzymes = list(RESTRICTION_ENZYMES.keys())

    count = 0
    seq_upper = sequence.upper()
    for enz_name in restriction_enzymes:
        site = RESTRICTION_ENZYMES.get(enz_name, "")
        if not site:
            continue
        site_upper = site.upper()
        # Only check concrete sites (no IUPAC ambiguity)
        if any(b not in "ACGT" for b in site_upper):
            continue
        # Forward strand
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
        # Reverse complement strand
        site_rc = reverse_complement(site_upper)
        if site_rc != site_upper:  # Avoid double-counting palindromes
            start = 0
            while True:
                pos = seq_upper.find(site_rc, start)
                if pos == -1:
                    break
                count += 1
                start = pos + 1

    return count


# ─── Adapter class ───────────────────────────────────────────────────


class DNAchiselAdapter:
    """Wraps DNAchisel's API to produce results comparable with BioCompiler.

    This adapter provides a uniform ``optimize()`` interface that mirrors
    biocompiler's optimization API, translating biocompiler-style
    constraints into DNAchisel specifications and converting results
    back into a common format.

    Usage::

        adapter = DNAchiselAdapter()
        result = adapter.optimize(
            protein="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
            organism="Homo_sapiens",
            constraints=[
                {"type": "gc_range", "gc_lo": 0.30, "gc_hi": 0.70},
                {"type": "avoid_restriction", "enzymes": ["EcoRI", "BamHI"]},
            ],
        )
        print(result.cai, result.gc_content, result.restriction_site_count)

    Raises:
        ImportError: If DNAchisel is not installed, with install instructions.
    """

    def __init__(self) -> None:
        """Initialize the adapter, checking DNAchisel availability."""
        if not _DNACHISEL_AVAILABLE:
            raise ImportError(
                "DNAchisel is not installed. "
                "Install with: pip install dnachisel"
            )

    def optimize(
        self,
        protein: str,
        organism: str = "Homo_sapiens",
        constraints: list[dict[str, Any]] | None = None,
    ) -> OptimizationResult:
        """Optimize a protein sequence using DNAchisel's constraint solver.

        This method translates biocompiler-style constraints into DNAchisel
        specifications, runs the optimization, and returns results in a
        format comparable with BioCompiler's output.

        Args:
            protein: Target protein sequence (single-letter amino acid codes).
            organism: Target organism for codon usage.
            constraints: List of constraint dicts. Each dict must have a
                ``"type"`` key. Supported types:

                - ``"gc_range"``: ``{"type": "gc_range", "gc_lo": 0.3, "gc_hi": 0.7}``
                - ``"avoid_restriction"``: ``{"type": "avoid_restriction", "enzymes": ["EcoRI", ...]}``
                - ``"preserve_translation"``: ``{"type": "preserve_translation"}``
                  (always added automatically)

                If None, defaults to GC range [0.30, 0.70] with the
                default enzyme set.

        Returns:
            OptimizationResult with optimized sequence and metrics.

        Raises:
            ImportError: If DNAchisel is not installed.
            ValueError: If protein contains invalid amino acid codes.
        """
        # Parse constraints into DNAchisel specs
        gc_lo = 0.30
        gc_hi = 0.70
        enzyme_names: list[str] = []

        if constraints:
            for c in constraints:
                ctype = c.get("type", "")
                if ctype == "gc_range":
                    gc_lo = c.get("gc_lo", 0.30)
                    gc_hi = c.get("gc_hi", 0.70)
                elif ctype == "avoid_restriction":
                    enzyme_names = c.get("enzymes", [])
                elif ctype == "preserve_translation":
                    pass  # Always enforced automatically

        # Limit enzyme count for performance
        if not enzyme_names:
            from ..constants import RESTRICTION_ENZYMES

            enzyme_names = list(RESTRICTION_ENZYMES.keys())[:MAX_RESTRICTION_ENZYMES]
        else:
            enzyme_names = enzyme_names[:MAX_RESTRICTION_ENZYMES]

        t0 = time.perf_counter()
        try:
            # Build initial sequence from preferred codons
            initial_seq = _build_initial_sequence(protein, organism)

            # Build constraint and objective specifications
            constraints, objectives = _build_dnachisel_spec(
                protein, organism, gc_lo, gc_hi, enzyme_names
            )

            # Record constraint class names for provenance
            constraint_names = [type(s).__name__ for s in constraints]
            constraint_names += [type(o).__name__ for o in objectives]

            # Create and solve the optimization problem
            problem = DnaOptimizationProblem(
                sequence=initial_seq,
                constraints=constraints,
                objectives=objectives,
            )

            # Resolve constraints first, then optimize objectives
            problem.resolve_constraints()
            if objectives:
                problem.optimize()

            # Extract the optimized sequence
            optimized = str(problem.sequence)

            elapsed = time.perf_counter() - t0

            # Compute CAI using validated evaluator for fairness —
            # DNAchisel's own CAI output is NOT trusted.
            from .metrics import compute_cai_validated
            from ..scanner import gc_content

            cai = compute_cai_validated(optimized, organism)
            gc = gc_content(optimized)
            rs_count = _count_restriction_sites(optimized, enzyme_names)

            return OptimizationResult(
                sequence=optimized,
                protein=protein,
                cai=cai,
                gc_content=gc,
                restriction_site_count=rs_count,
                execution_time_s=elapsed,
                success=True,
                constraints_applied=constraint_names,
            )

        except ImportError as exc:
            raise ImportError(
                "DNAchisel is not installed. "
                "Install with: pip install dnachisel"
            ) from exc
        except Exception as exc:
            elapsed = time.perf_counter() - t0
            logger.error(
                "DNAchisel optimization failed: %s", exc, exc_info=True
            )
            return OptimizationResult(
                sequence="",
                protein=protein,
                cai=0.0,
                gc_content=0.0,
                restriction_site_count=0,
                execution_time_s=elapsed,
                success=False,
                error=str(exc),
            )
