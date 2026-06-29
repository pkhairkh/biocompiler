"""
BioCompiler Application — CLI Service Functions
=================================================
Service functions extracted from the CLI command handlers (Wave 4b SoC refactoring).

These orchestrate business logic for the CLI layer:
- run_optimization: single protein optimization
- run_batch_optimization: batch protein optimization
- run_check_predicates: predicate checking
- resolve_organism_arg / resolve_source_organism_arg: organism resolution
- clear_engine_caches: engine cache management
- format_optimization_json / format_batch_json: JSON result formatting

This is separate from the API-layer services in optimization_service.py,
assessment_service.py, etc. which serve the REST API routes.
"""

from __future__ import annotations

import json
import logging
import random
from dataclasses import dataclass, field
from typing import List, Optional

logger = logging.getLogger(__name__)


# ── Data classes for structured results ──────────────────────────────────────

@dataclass
class OptimizationResult:
    """Structured result from a single optimization run."""
    optimized: str
    pred_results: list = field(default_factory=list)
    cert_text: str = ""
    cai: Optional[float] = None
    gc_content: Optional[float] = None
    codon_pair_bias: Optional[float] = None
    strategy: str = "hybrid"
    organism: str = ""
    protein: Optional[str] = None
    input_seq: Optional[str] = None


@dataclass
class BatchOptimizationResult:
    """Structured result from a batch optimization run."""
    results: list = field(default_factory=list)  # List of dicts
    organism: str = ""
    strategy: str = "hybrid"


@dataclass
class CheckResult:
    """Structured result from a predicate check run."""
    type_results: list = field(default_factory=list)
    cert_text: str = ""
    cert_level: str = ""
    predicate_results: list = field(default_factory=list)


# ── Organism resolution ──────────────────────────────────────────────────────

def resolve_organism_arg(raw: Optional[str], species: Optional[str] = None) -> str:
    """Resolve the organism from --organism or --species (alias).

    --organism takes precedence; --species is kept as a backward-compatible
    alias.  The value is normalised through ``resolve_organism`` so that
    shorthand names like 'ecoli' or 'human' are accepted.
    """
    from ..organisms import resolve_organism as _resolve
    value = raw or species or "Homo_sapiens"
    try:
        return _resolve(value)
    except Exception:
        return value


def resolve_source_organism_arg(raw: Optional[str]) -> Optional[str]:
    """Resolve the source organism from --source-organism.

    Returns None if --source-organism is not specified.
    """
    if raw is None:
        return None
    from ..organisms import resolve_organism as _resolve
    try:
        return _resolve(raw)
    except Exception:
        return raw


# ── Engine cache management ─────────────────────────────────────────────────

def clear_engine_caches() -> None:
    """Clear engine caches for a fresh optimization run."""
    try:
        from biocompiler.engines.foldx import clear_cache as foldx_clear_cache
        foldx_clear_cache()
    except ImportError:
        logger.debug("foldx module not available; skipping cache clear")
    try:
        from biocompiler.engines.camsol import clear_cache as camsol_clear_cache
        camsol_clear_cache()
    except ImportError:
        logger.debug("camsol module not available; skipping cache clear")
    try:
        from biocompiler.immunogenicity.core import clear_cache as immunogenicity_clear_cache
        immunogenicity_clear_cache()
    except ImportError:
        logger.debug("immunogenicity module not available; skipping cache clear")


# ── Optimization orchestration ───────────────────────────────────────────────

def run_optimization(
    *,
    protein: Optional[str] = None,
    input_seq: Optional[str] = None,
    organism: str = "Homo_sapiens",
    strategy: str = "hybrid",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    enzymes: Optional[List[str]] = None,
    no_splice_check: bool = False,
    use_codon_pair_bias: bool = False,
    source_organism: Optional[str] = None,
    therapeutic: bool = False,
    organism_domain_raw: str = "auto",
    splice_low: float = 3.0,
    splice_high: float = 6.0,
    avoid_gt: bool = True,
    seed: Optional[int] = None,
) -> OptimizationResult:
    """Run a single gene optimization.

    This is the core business logic extracted from cmd_optimize.
    It orchestrates the choice of optimizer backend and returns
    a structured result without any I/O or formatting.
    """
    if seed is not None:
        random.seed(seed)
        logger.info("Random seed set to %d for reproducible optimization", seed)

    # Clear caches
    clear_engine_caches()

    # Resolve organism domain
    if no_splice_check and organism_domain_raw == "auto":
        organism_domain_raw = "prokaryote"

    from ..api import resolve_organism_domain
    resolved_domain = resolve_organism_domain(organism, organism_domain_raw)

    enzymes = enzymes or []

    # ── Strategy dispatch ──────────────────────────────────────────────────
    if strategy == "csp":
        from ..solver.dispatch import csp_optimize, is_solver_available
        if not is_solver_available():
            raise RuntimeError("CSP solver not available. Install ortools or z3-solver.")

        opt_result = csp_optimize(
            protein if protein else "",
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
        )
        optimized = opt_result.sequence if hasattr(opt_result, "sequence") else str(opt_result)
        pred_results = []
        cert_text = ""
        cai_val = getattr(opt_result, "cai", None)

    elif strategy in ("constraint_first", "hybrid"):
        if protein:
            from biocompiler.optimizer import optimize_sequence
            opt_result = optimize_sequence(
                protein, organism=organism,
                gc_lo=gc_lo, gc_hi=gc_hi,
                consider_codon_pair_bias=use_codon_pair_bias,
                source_organism=source_organism,
                therapeutic=therapeutic,
            )
            optimized = opt_result.sequence
            pred_results = opt_result.predicate_results
            cert_text = opt_result.certificate_text
            cai_val = getattr(opt_result, "cai", None)
        else:
            # Legacy FASTA-input path
            from biocompiler.optimizer import BioOptimizer
            opt = BioOptimizer(
                species=organism,
                enzymes=enzymes,
                splice_low=splice_low,
                splice_high=splice_high,
                avoid_gt=avoid_gt,
                organism_domain=resolved_domain,
            )
            optimized, pred_results, cert_text = opt.optimize(input_seq)
            cai_val = None
    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    # ── Codon-pair bias scoring ────────────────────────────────────────────
    cpb_score: Optional[float] = None
    if use_codon_pair_bias:
        try:
            from biocompiler.expression.codon_pair_scoring import compute_cpb
            cpb_score = compute_cpb(optimized)
        except (ImportError, Exception):
            cpb_score = None

    # ── GC content ─────────────────────────────────────────────────────────
    from biocompiler.sequence.scanner import gc_content
    gc_val = gc_content(optimized) if optimized else 0.0

    return OptimizationResult(
        optimized=optimized,
        pred_results=pred_results,
        cert_text=cert_text,
        cai=cai_val,
        gc_content=gc_val,
        codon_pair_bias=cpb_score,
        strategy=strategy,
        organism=organism,
        protein=protein,
        input_seq=input_seq,
    )


def run_batch_optimization(
    *,
    proteins: List[tuple],  # List of (name, sequence)
    organism: str = "Homo_sapiens",
    strategy: str = "hybrid",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    use_codon_pair_bias: bool = False,
    source_organism: Optional[str] = None,
    therapeutic: bool = False,
) -> BatchOptimizationResult:
    """Run batch optimization for multiple proteins."""
    from biocompiler.optimizer import optimize_sequence

    results: List[dict] = []
    for name, prot_seq in proteins:
        try:
            opt_result = optimize_sequence(
                prot_seq, organism=organism,
                gc_lo=gc_lo, gc_hi=gc_hi,
                consider_codon_pair_bias=use_codon_pair_bias,
                source_organism=source_organism,
                therapeutic=therapeutic,
            )
            results.append({
                "name": name,
                "status": "ok",
                "sequence": opt_result.sequence,
                "gc_content": opt_result.gc_content,
                "cai": opt_result.cai,
                "codon_pair_bias": opt_result.codon_pair_bias if use_codon_pair_bias else None,
            })
        except Exception as exc:
            results.append({
                "name": name,
                "status": "error",
                "error": str(exc),
            })

    return BatchOptimizationResult(
        results=results,
        organism=organism,
        strategy=strategy,
    )


# ── Predicate checking ───────────────────────────────────────────────────────

def run_check_predicates(
    *,
    seq: str,
    organism: str,
    enzymes: Optional[List[str]] = None,
    cryptic_threshold: float = 3.0,
    uncertain_lo: float = 3.0,
    predicate_filter: Optional[List[str]] = None,
    species: str = "human",
) -> CheckResult:
    """Run predicate checks on a DNA sequence.

    This is the core business logic extracted from cmd_check.
    """
    from ..type_system import (
        evaluate_all_predicates,
        registry as predicate_registry,
        PredicateResult,
    )
    from biocompiler.provenance.certificate import format_certificate, compute_certificate

    type_results = evaluate_all_predicates(
        seq=seq,
        organism=organism,
        enzymes=enzymes if enzymes else None,
        cryptic_threshold=cryptic_threshold,
        uncertain_lo=uncertain_lo,
    )

    # Filter by predicate names if specified
    if predicate_filter:
        registered = set(predicate_registry.names())
        unknown = [p for p in predicate_filter if p not in registered]
        if unknown:
            raise ValueError(f"Unknown predicates: {', '.join(unknown)}")
        type_results = [
            r for r in type_results
            if any(r.predicate == p or r.predicate.startswith(p + "(") for p in predicate_filter)
        ]

    # Convert TypeCheckResult → PredicateResult for certificate
    results: List[PredicateResult] = []
    for r in type_results:
        results.append(PredicateResult(
            predicate=r.predicate,
            passed=r.passed,
            verdict=r.verdict,
            details=r.violation or "",
        ))

    cert_text = format_certificate(results, seq, species)
    cert_level = compute_certificate(results)

    return CheckResult(
        type_results=type_results,
        cert_text=cert_text,
        cert_level=cert_level.value if cert_level else "",
        predicate_results=results,
    )


# ── JSON formatting ─────────────────────────────────────────────────────────

def format_optimization_json(
    result: OptimizationResult,
    *,
    no_splice_check: bool = False,
    source_organism: Optional[str] = None,
    therapeutic: bool = False,
) -> str:
    """Format an OptimizationResult as a JSON string."""
    from .. import __version__

    result_dict: dict = {
        "version": __version__,
        "organism": result.organism,
        "strategy": result.strategy,
        "gc_content": round(result.gc_content, 4) if result.gc_content else None,
        "sequence_length": len(result.optimized) if result.optimized else 0,
        "sequence": result.optimized,
        "no_splice_check": no_splice_check,
        "codon_pair_bias": result.codon_pair_bias,
        "source_organism": source_organism,
        "therapeutic": therapeutic,
    }
    if result.pred_results:
        result_dict["predicate_results"] = [
            {"name": p.predicate, "passed": p.passed, "details": p.details}
            for p in result.pred_results
        ]
    if result.cert_text:
        result_dict["certificate_text"] = result.cert_text
    return json.dumps(result_dict, indent=2)


def format_batch_json(
    result: BatchOptimizationResult,
    *,
    no_splice_check: bool = False,
    total_proteins: int = 0,
    source_organism: Optional[str] = None,
    therapeutic: bool = False,
) -> str:
    """Format a BatchOptimizationResult as a JSON string."""
    from .. import __version__

    output_dict = {
        "version": __version__,
        "organism": result.organism,
        "strategy": result.strategy,
        "no_splice_check": no_splice_check,
        "total_proteins": total_proteins,
        "source_organism": source_organism,
        "therapeutic": therapeutic,
        "results": result.results,
    }
    return json.dumps(output_dict, indent=2)
