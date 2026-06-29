"""
Batch optimization function for the optimizer pipeline.

Contains the batch_optimize() function for optimizing multiple proteins
efficiently using a shared HybridOptimizer instance.

Decomposition: Extracted from pipeline_core.py (Task pipeline-decompose).
"""

from typing import Any, Dict, List, Optional

import logging

from ..type_system import AA_TO_CODONS, PredicateResult
from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
from biocompiler.provenance.certificate import format_certificate
from biocompiler.shared.exceptions import InvalidProteinError
from .constraints import _organism_to_species_key
from .utils import OptimizationResult

logger = logging.getLogger(__name__)


def batch_optimize(
    proteins: list[str],
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: list | None = None,
    cai_threshold: float = 0.5,
    strategy: str = "hybrid",
    full_postprocessing: bool = True,
    skip_biosecurity_check: bool = False,
    **kwargs,
) -> list[OptimizationResult]:
    """Optimize multiple proteins in batch.

    Each result is processed through the full :func:`optimize_sequence`
    pipeline (mRNA stability, codon pair bias, CAI recovery, UTR
    suggestions, provenance tracking). This delegates to
    :func:`optimize_sequence` per protein, providing full feature parity
    with the single-protein path.

    Note: A streamlined ``HybridOptimizer``-based batch path was removed
    in v0.9.0 when ``hybrid_optimizer.py`` was deleted as part of the
    legacy slow-path cleanup. ``full_postprocessing`` is retained as a
    keyword argument (default ``True``) for backward compatibility but
    the streamlined ``False`` branch is no longer supported — passing
    ``False`` now also routes through :func:`optimize_sequence`.

    Args:
        proteins: List of amino acid sequences (single-letter codes, no stop).
        organism: Target organism. Accepts canonical binomials, short keys,
            abbreviated binomials, or display names (same as
            :func:`optimize_sequence`).
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        enzymes: List of restriction enzyme names to avoid.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        strategy: Optimization strategy (default 'hybrid').
        full_postprocessing: Kept for backward compatibility. Always
            effectively ``True``; the streamlined ``False`` path was
            removed with ``HybridOptimizer`` in v0.9.0.
        **kwargs: Additional arguments passed to :func:`optimize_sequence`.

    Returns:
        List of OptimizationResult, one per input protein, in the same order.

    Raises:
        InvalidProteinError: if any protein contains invalid amino acid codes.
        UnsupportedOrganismError: if the organism is not supported.

    Examples:
        >>> from biocompiler.optimizer import batch_optimize
        >>> proteins = ['MSKGEELFTG', 'MALWMRLLPL', 'MVHLTPEEKS']
        >>> results = batch_optimize(proteins, organism='Escherichia_coli')
        >>> assert len(results) == 3
        >>> for r in results:
        ...     assert r.cai > 0.5
    """
    # Import here to avoid circular imports
    from .pipeline_core import optimize_sequence

    if not proteins:
        return []

    # Resolve organism name once
    organism = resolve_organism(organism, strict=False)

    # Validate all proteins upfront
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    for protein in proteins:
        protein_upper = protein.strip().upper()
        invalid = set(protein_upper) - valid_aas
        if invalid:
            raise InvalidProteinError(protein, invalid)

    # ── Per-protein path ───────────────────────────────────────────
    # Delegate to optimize_sequence per protein. This is the only path
    # since the streamlined HybridOptimizer branch was removed in v0.9.0
    # (hybrid_optimizer.py deleted as part of legacy slow-path cleanup).
    # ``full_postprocessing`` is retained as a kwarg for backward
    # compatibility but no longer triggers a different code path.
    results: list[OptimizationResult] = []
    for protein in proteins:
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
            enzymes=enzymes,
            strategy=strategy,
            skip_biosecurity_check=skip_biosecurity_check,
            **kwargs,
        )
        results.append(result)
    return results
