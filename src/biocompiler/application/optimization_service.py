"""
BioCompiler Application — Optimization service.

Orchestrates type-checking and optimization, calling the optimizer,
type_system, certificate generator, etc.
"""

import logging
from typing import Any, Optional

from biocompiler.sequence.scanner import gc_content
from biocompiler.expression.translation import translate
from ..type_system import evaluate_all_predicates
from biocompiler.provenance.certificate import generate_certificate
from biocompiler.optimizer import optimize_sequence
from biocompiler.shared.types import Verdict
from biocompiler.shared.constants import RESTRICTION_ENZYMES
from biocompiler.shared.exceptions import CertificateGenerationError
from biocompiler.organisms.config import is_eukaryotic_organism

logger = logging.getLogger(__name__)


def resolve_organism_domain(organism: str, organism_domain: str = "auto") -> str:
    """Resolve the effective organism domain from user input.

    When *organism_domain* is ``"auto"``, the domain is detected from
    the organism name using :func:`is_eukaryotic_organism`.  When
    explicitly set to ``"eukaryote"`` or ``"prokaryote"``, that value
    is used regardless of the organism name.

    Args:
        organism: Organism name (e.g. ``"Escherichia_coli"``).
        organism_domain: One of ``"auto"``, ``"eukaryote"``, or
            ``"prokaryote"``.  Defaults to ``"auto"``.

    Returns:
        The resolved domain string: ``"eukaryote"`` or ``"prokaryote"``.

    Raises:
        ValueError: If *organism_domain* is not one of the valid choices.
    """
    if organism_domain not in ("auto", "eukaryote", "prokaryote"):
        raise ValueError(
            f"Invalid organism_domain: {organism_domain!r}. "
            "Must be one of: auto, eukaryote, prokaryote"
        )

    if organism_domain == "auto":
        return "eukaryote" if is_eukaryotic_organism(organism) else "prokaryote"

    return organism_domain


def type_check_sequence(
    seq: str,
    organism: str,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: Optional[list[str]] = None,
    cellular_context: str = "HEK293T",
) -> dict[str, Any]:
    """Run type-check predicates on a DNA sequence.

    Returns a dict matching the TypeCheckResponse structure:
    {
        "sequence_length": int,
        "gc_content": float,
        "protein": str,
        "results": list[dict],
        "overall_verdict": str,
        "certificate": Optional[dict],
    }
    """
    seq = seq.upper()
    exon_boundaries = exon_boundaries or [(0, len(seq))]

    # Type check
    results = evaluate_all_predicates(
        seq=seq,
        boundaries=exon_boundaries,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        enzymes=enzymes,
    )

    # Overall verdict
    verdicts = [r.verdict for r in results]
    overall = Verdict.PASS
    for v in verdicts:
        if overall == Verdict.FAIL or v == Verdict.FAIL:
            overall = Verdict.FAIL
        elif overall == Verdict.UNCERTAIN or v == Verdict.UNCERTAIN:
            overall = Verdict.UNCERTAIN

    # Build response
    result_dicts = [
        {
            "predicate": r.predicate,
            "verdict": r.verdict.value,
            "violation": r.violation,
            "knowledge_gap": r.knowledge_gap,
        }
        for r in results
    ]

    # Protein translation
    coding_seq = "".join(seq[start:end] for start, end in exon_boundaries)
    protein = translate(coding_seq)

    # Generate certificate if all PASS
    cert_dict = None
    if overall == Verdict.PASS:
        try:
            cert = generate_certificate(
                coding_seq, results,
                {
                    "organism": organism,
                    "exon_boundaries": exon_boundaries,
                    "gc_lo": gc_lo,
                    "gc_hi": gc_hi,
                    "cai_threshold": cai_threshold,
                    "enzymes": enzymes or list(RESTRICTION_ENZYMES.keys()),
                    "cell_type": cellular_context,
                },
            )
            cert_dict = cert.to_dict()
        except CertificateGenerationError as exc:
            logger.warning("Certificate generation failed: %s", exc)

    return {
        "sequence_length": len(seq),
        "gc_content": gc_content(seq),
        "protein": protein[:100] + ("..." if len(protein) > 100 else ""),
        "results": result_dicts,
        "overall_verdict": overall.value,
        "certificate": cert_dict,
    }


def optimize_protein(
    protein: str,
    organism: str,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.2,
    restriction_sites: Optional[list[str]] = None,
    cryptic_splice_threshold: float = 3.0,
    organism_domain: str = "auto",
    source_organism: Optional[str] = None,
    therapeutic: bool = False,
    self_protein: Optional[bool] = None,
    track_provenance: bool = True,
    strict_mode: bool = True,
) -> dict[str, Any]:
    """Run optimization for a target protein.

    Returns a dict matching the OptimizeResponse structure:
    {
        "sequence": str,
        "protein": str,
        "cai": float,
        "gc_content": float,
        "satisfied_predicates": list[str],
        "failed_predicates": list[str],
        "fallback_used": bool,
        "organism_domain": str,
        "source_organism": Optional[str],
        "therapeutic": bool,
        "self_protein": Optional[bool],
        "decision_trail": Optional[OptimizationDecisionTrail],
    }
    """
    resolved_domain = resolve_organism_domain(organism, organism_domain)

    result = optimize_sequence(
        target_protein=protein,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        restriction_sites=restriction_sites,
        cryptic_splice_threshold=cryptic_splice_threshold,
        organism_domain=resolved_domain,
        source_organism=source_organism,
        therapeutic=therapeutic,
        self_protein=self_protein,
        track_provenance=track_provenance,
        strict_mode=strict_mode,
    )

    return {
        "sequence": result.sequence,
        "protein": result.protein,
        "cai": result.cai,
        "gc_content": result.gc_content,
        "satisfied_predicates": result.satisfied_predicates,
        "failed_predicates": result.failed_predicates,
        "fallback_used": result.fallback_used,
        "organism_domain": resolved_domain,
        "source_organism": source_organism,
        "therapeutic": therapeutic,
        "self_protein": self_protein,
        "decision_trail": result.decision_trail,
    }


def type_check_batch_item(
    sequence: str,
    organism: str,
    exon_boundaries: Optional[list[tuple[int, int]]] = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: Optional[list[str]] = None,
    cellular_context: str = "HEK293T",
) -> dict[str, Any]:
    """Process a single type-check item for batch processing."""
    return type_check_sequence(
        seq=sequence,
        organism=organism,
        exon_boundaries=exon_boundaries,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        enzymes=enzymes,
        cellular_context=cellular_context,
    )


def optimize_batch_item(
    protein: str,
    organism: str,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.2,
    enzymes: Optional[list[str]] = None,
    cryptic_splice_threshold: float = 3.0,
    organism_domain: str = "auto",
    source_organism: Optional[str] = None,
    therapeutic: bool = False,
    self_protein: Optional[bool] = None,
) -> dict[str, Any]:
    """Process a single optimization item for batch processing."""
    resolved_domain = resolve_organism_domain(organism, organism_domain)

    result = optimize_sequence(
        target_protein=protein,
        organism=organism,
        gc_lo=gc_lo,
        gc_hi=gc_hi,
        cai_threshold=cai_threshold,
        restriction_sites=enzymes,
        cryptic_splice_threshold=cryptic_splice_threshold,
        organism_domain=resolved_domain,
        source_organism=source_organism,
        therapeutic=therapeutic,
        self_protein=self_protein,
    )

    return {
        "sequence": result.sequence,
        "protein": result.protein,
        "cai": result.cai,
        "gc_content": result.gc_content,
        "satisfied_predicates": result.satisfied_predicates,
        "failed_predicates": result.failed_predicates,
        "fallback_used": result.fallback_used,
        "organism_domain": resolved_domain,
        "source_organism": source_organism,
        "therapeutic": therapeutic,
        "self_protein": self_protein,
    }
