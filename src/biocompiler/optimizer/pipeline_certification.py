"""Post-optimization certification for the integrated fast path.

evaluate all 20 predicates on the fast-path's optimized
DNA, generate the certificate, populate provenance, and enforce ``strict_mode``.
Mirrors the slow path's predicate-evaluation block at
``pipeline_core.py:1017-1509`` so the :class:`OptimizationResult` field
contract is identical across paths.

This module exists as a dedicated, unit-testable helper (per W2-c spec §3.2)
rather than inlining ~150 lines of certification logic into the fast-path
``try:`` block in ``pipeline_core.py``.

Public API:
    - :func:`certify_fast_path_result` — populate certificate / predicate /
      provenance fields on an :class:`OptimizationResult` in place, and
      enforce ``strict_mode`` by raising :class:`OptimizationConstraintError`
      when any predicate fails.
"""
from __future__ import annotations

import logging
import time as _time
from datetime import datetime, timezone
from typing import Any

from ..type_system import (
    PredicateResult,
    Verdict,
    check_no_stop_codons,
    check_valid_coding_seq,
    check_no_restriction_site,
    check_no_cryptic_splice,
    check_no_cpg_island,
    check_no_gt_dinucleotide_soft,
)
from .pipeline_paths import evaluate_extended_predicates
from .utils import OptimizationResult
from ..provenance import OptimizationRecord, _get_biocompiler_version
from ..provenance.certificate import format_certificate
from ..provenance.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    DecisionProvenanceCollector,
    OptimizationDecisionTrail,
)
from ..sequence.sliding_gc import check_sliding_gc
from ..shared.exceptions import OptimizationConstraintError


_DEFAULT_ENZYMES = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]


def certify_fast_path_result(
    result: OptimizationResult,
    target_protein: str,
    organism: str,
    species_key: str,
    gc_lo: float,
    gc_hi: float,
    cai_threshold: float,
    enzymes: list[str],
    effective_splice_low: float,
    effective_avoid_gt: bool,
    is_prokaryote: bool,
    gc_window_size: int,
    gc_window_min: float | None,
    gc_window_max: float | None,
    cpg_mode: str,
    strict_mode: bool,
    track_provenance: bool,
    include_utr: bool,
    seed: int | None,
    start_time: float,
    tissue: str = "",
    enable_mutagenesis: bool = False,
    logger: logging.Logger | None = None,
) -> None:
    """Populate certificate / predicate / provenance fields on ``result`` in place.

    Mirrors the slow path's predicate evaluation + certificate generation +
    provenance population block (``pipeline_core.py:1017-1509``) so the
    :class:`OptimizationResult` field contract is identical across the
    integrated fast path and the sequential slow path.

    The integrated optimizer (``integrated_optimize``) is a greedy
    forward-pass + multi-pass cleanup; it never applies amino-acid
    mutagenesis, so ``aa_substitutions=[]``, ``mutagenesis_applied=False``,
    and ``fallback_used=False`` unconditionally. mRNA stability, codon-pair
    bias, and custom-objective passes are NOT run on the fast path (callers
    wanting those features should pass ``use_integrated=False``); the
    corresponding result fields stay at their dataclass defaults.

    Args:
        result: The :class:`OptimizationResult` produced by the fast path
            (must already have ``sequence``, ``cai``, ``gc_content``, and
            ``protein`` populated). Mutated in place.
        target_protein: The original input amino-acid sequence.
        organism: Target organism identifier (e.g. ``"h_sapiens"``).
        species_key: Short species key for certificate rendering
            (e.g. ``"human"``); produced by ``_organism_to_species_key``.
        gc_lo, gc_hi: GC content acceptance range (0–1).
        cai_threshold: Minimum CAI for the ``CodonOptimality`` predicate.
        enzymes: Restriction enzymes to avoid (``None`` → default 5).
        effective_splice_low: MaxEntScan donor-score threshold below which
            a splice site is considered weak. ``999.0`` disables the check.
        effective_avoid_gt: Whether GT dinucleotide avoidance is requested.
        is_prokaryote: Skip eukaryote-specific predicates (NoCrypticSplice,
            NoCpGIsland, NoGTDinucleotide) when ``True``.
        gc_window_size, gc_window_min, gc_window_max: Sliding-window GC
            parameters (``window_size <= 0`` disables the check).
        cpg_mode: CpG elimination mode (``"aggressive"`` or ``"off"``).
            Currently informational on this path — the integrated optimizer
            already ran its own CpG sweep; we just record the mode.
        strict_mode: If ``True``, raise :class:`OptimizationConstraintError`
            with ``partial_result=result`` when any predicate fails.
        track_provenance: If ``True``, build an
            :class:`OptimizationDecisionTrail` via
            :class:`DecisionProvenanceCollector`. The
            :class:`OptimizationRecord` provenance is always populated.
        include_utr: If ``True``, populate ``suggested_5utr`` /
            ``suggested_3utr`` / ``utr_score_5`` / ``utr_score_3``.
        seed: RNG seed used for this run (``None`` if none provided).
        start_time: monotonic-clock start time of the optimization
            (used for the ``solve_time`` field of the provenance record).
        tissue: Optional tissue context for tissue-aware extended
            predicates (e.g. miRNA binding).
        enable_mutagenesis: Currently informational — the integrated
            optimizer does not apply amino-acid mutagenesis today.
        logger: Logger instance. If ``None``, a module logger is used.

    Raises:
        OptimizationConstraintError: If ``strict_mode=True`` and any
            predicate fails. The ``partial_result`` attribute is the
            populated ``result`` (with certificate + predicate_results
            already filled in, so callers can inspect what went wrong).
    """
    seq = result.sequence
    if logger is None:
        logger = logging.getLogger(__name__)

    # ── Predicates 1–10 (mirror slow-path pipeline_core.py:1043-1138) ──
    pred_results: list[PredicateResult] = []

    # 1. NoStopCodons
    pred_results.append(check_no_stop_codons(seq))

    # 2. NoCrypticSplice (skip for prokaryotes — no spliceosome)
    if is_prokaryote:
        pred_results.append(PredicateResult(
            "NoCrypticSplice", True,
            details="Skipped for prokaryotic organism",
        ))
    else:
        pred_results.append(check_no_cryptic_splice(
            seq, low_thresh=effective_splice_low, organism=organism,
        ))

    # 3. NoCpGIsland (skip for prokaryotes)
    if is_prokaryote:
        pred_results.append(PredicateResult(
            "NoCpGIsland", True,
            details="Skipped for prokaryotic organism",
        ))
    else:
        pred_results.append(check_no_cpg_island(seq, organism=organism))

    # 4. NoRestrictionSite
    pred_results.append(check_no_restriction_site(
        seq, enzymes or _DEFAULT_ENZYMES,
    ))

    # 5. NoGTDinucleotide (skip for prokaryotes; soft for eukaryotes)
    if is_prokaryote:
        pred_results.append(PredicateResult(
            "NoGTDinucleotide", True,
            details="Skipped for prokaryotic organism",
        ))
    elif effective_avoid_gt:
        pred_results.append(check_no_gt_dinucleotide_soft(
            seq, organism=organism,
        ))
    else:
        pred_results.append(PredicateResult(
            "NoGTDinucleotide", True,
            details="GT avoidance not requested",
        ))

    # 6. ValidCodingSeq
    pred_results.append(check_valid_coding_seq(seq))

    # 7. ConservationScore (auto-pass: integrated optimizer never mutagenizes)
    pred_results.append(PredicateResult(
        "ConservationScore", True,
        details="All AA conservation scores >= 0",
    ))

    # 8. CodonOptimality (use CAI already computed by the fast path)
    cai_val = result.cai
    pred_results.append(PredicateResult(
        "CodonOptimality", cai_val >= cai_threshold,
        details=f"CAI={cai_val:.4f}, min={cai_threshold}",
    ))

    # 9. GCInRange (use GC already computed by the fast path)
    gc_val = result.gc_content
    pred_results.append(PredicateResult(
        "GCInRange", gc_lo <= gc_val <= gc_hi,
        details=f"GC content: {gc_val:.3f} (range [{gc_lo}, {gc_hi}])",
    ))

    # 10. SlidingGC
    _eff_w_min = gc_window_min if gc_window_min is not None else gc_lo
    _eff_w_max = gc_window_max if gc_window_max is not None else gc_hi
    if gc_window_size > 0 and len(seq) >= gc_window_size:
        _sgc = check_sliding_gc(
            seq, window_size=gc_window_size,
            gc_min=_eff_w_min, gc_max=_eff_w_max,
        )
        pred_results.append(PredicateResult(
            "SlidingGC", _sgc.passed,
            details=(
                f"Window={gc_window_size}, range=[{_eff_w_min:.2f}, {_eff_w_max:.2f}], "
                f"min_gc={_sgc.min_gc:.3f}, max_gc={_sgc.max_gc:.3f}, "
                f"violations={len(_sgc.violations)}"
            ),
        ))
    else:
        pred_results.append(PredicateResult(
            "SlidingGC", True,
            details="Sliding-window GC check skipped (window_size=0 or sequence too short)",
        ))

    # ── Predicates 11–20 (extended diagnostic) ──
    # evaluate_extended_predicates appends NoInstabilityMotif,
    # NoCrypticPromoter, NoUnexpectedTMDomain, NoCrypticORF, NoRQCTrigger,
    # NoAluRepeat, NoMiRNABindingSite, NoM6ASite, NoPolyASignal,
    # mRNASecondaryStructure. Each individual predicate is wrapped in
    # try/except inside the helper, so evaluation cannot propagate
    # exceptions (it appends a synthetic "Unavailable" pass on failure).
    #
    # NOTE: skipped for prokaryotes to mirror the slow path's
    # run_prokaryote_hybrid_path (pipeline_paths.py:414-488), which does
    # NOT call evaluate_extended_predicates. This preserves field parity
    # (prokaryote result has 10 predicate_results on both paths) and
    # avoids surfacing eukaryote-only failures (MRNAStability,
    # NoRQCTrigger) that would make strict_mode=True raise on prokaryotic
    # inputs the slow path accepts.
    if not is_prokaryote:
        pred_results = evaluate_extended_predicates(
            seq, organism, pred_results, tissue=tissue,
        )

    # ── Certificate text (mirror slow-path pipeline_core.py:1148) ──
    cert_text = format_certificate(pred_results, seq, species_key)

    # ── Uncertain-predicate count + warning (mirror _emit_uncertainty_warnings) ──
    uncertain_count = sum(
        1 for r in pred_results
        if getattr(r, 'verdict', None) == Verdict.UNCERTAIN
    )
    if uncertain_count:
        _uncertain_names = [
            r.predicate for r in pred_results
            if getattr(r, 'verdict', None) == Verdict.UNCERTAIN
        ]
        logger.warning(
            "%d predicate(s) returned UNCERTAIN verdict: %s. "
            "Certificate level is capped. Install external tools "
            "(FoldX, ViennaRNA, ESMFold) for higher-confidence verdicts.",
            uncertain_count,
            ", ".join(_uncertain_names[:5])
            + ("..." if len(_uncertain_names) > 5 else ""),
        )

    failed = [r.predicate for r in pred_results if not r.passed]
    satisfied = [r.predicate for r in pred_results if r.passed]

    # ── Provenance record (always populated — cheap; mirror slow-path
    # pipeline_core.py:1398-1412) ──
    provenance_record = OptimizationRecord(
        input_sequence=target_protein,
        output_sequence=seq,
        organism=organism,
        constraints_applied=sorted({r.predicate for r in pred_results}),
        mutations_made=[],  # integrated optimizer never applies mutagenesis today
        solver_backend="integrated",
        solve_time=round(_time.monotonic() - start_time, 6),
        seed_used=seed,
        timestamp=datetime.now(timezone.utc).isoformat(),
        biocompiler_version=_get_biocompiler_version(),
    )

    # ── Decision trail (when track_provenance=True; mirror slow-path
    # pipeline_core.py:1434-1479) ──
    decision_trail: OptimizationDecisionTrail | None = None
    if track_provenance:
        try:
            _collector = DecisionProvenanceCollector()
            _collector.start_optimization(
                protein=target_protein,
                organism=organism,
                constraints=[r.predicate for r in pred_results],
                solver_backend="integrated",
                seed=seed,
            )
            # ── Record one CodonDecision per codon position ──
            # The integrated optimizer is a greedy CAI-first forward pass;
            # it DOES make a codon choice at each position, so we record
            # one CodonDecision per codon. The slow path's prokaryote
            # hybrid records these inside HybridProkaryoteOptimizer
            # (hybrid_prokaryote.py:176); the slow path's eukaryote hybrid
            # does NOT record them (a known slow-path gap). Recording them
            # here for both organism types gives the fast path a strictly
            # richer decision trail and satisfies the existing
            # test_optimize_decision_trail_has_codon_decisions contract
            # (expects N codon decisions for an N-AA protein on ecoli).
            #
            # alternatives_considered is reconstructed from the organism's
            # codon-adaptiveness table (not exposed by integrated_optimize
            # directly). Each alternative dict matches the slow path's
            # format (hybrid_prokaryote.py:155-161): codon, cai_contribution,
            # gc_contribution, violates_constraints, rejected_because.
            from biocompiler.shared.constants import AA_TO_CODONS
            from biocompiler.organisms import (
                CODON_ADAPTIVENESS_TABLES, resolve_organism,
            )
            _cai_table = CODON_ADAPTIVENESS_TABLES.get(organism, {})
            _n_codons = len(seq) // 3
            for _pos in range(_n_codons):
                _codon = seq[_pos * 3:_pos * 3 + 3]
                _aa = target_protein[_pos] if _pos < len(target_protein) else ""
                _chosen_cai = _cai_table.get(_codon, 0.0)
                # Build alternatives_considered from the synonymous codon set
                _alts: list[dict[str, Any]] = []
                for _alt_codon in AA_TO_CODONS.get(_aa, []):
                    _alt_cai = _cai_table.get(_alt_codon, 0.0)
                    _alt_gc = (
                        (_alt_codon.count("G") + _alt_codon.count("C")) / 3.0
                    )
                    if _alt_codon == _codon:
                        _rejected = ""
                        _violates = False
                    else:
                        # Best-effort: the integrated optimizer is CAI-first,
                        # so non-chosen alternatives are most often rejected
                        # for lower CAI. Constraint-violation detection per
                        # alternative would require substituting and
                        # re-checking — too expensive for the fast path.
                        _rejected = (
                            "lower CAI" if _alt_cai < _chosen_cai
                            else "constraint conflict"
                        )
                        _violates = False
                    _alts.append({
                        "codon": _alt_codon,
                        "cai_contribution": round(_alt_cai, 4),
                        "gc_contribution": round(_alt_gc, 2),
                        "violates_constraints": _violates,
                        "rejected_because": _rejected,
                    })
                # cai_impact = chosen_cai - max_alt_cai (negative = CAI lost
                # to constraints; 0.0 = best-CAI codon chosen)
                _alt_cais = [
                    _cai_table.get(a["codon"], 0.0)
                    for a in _alts if a["codon"] != _codon
                ]
                _best_alt_cai = max(_alt_cais) if _alt_cais else _chosen_cai
                _cai_impact = round(_chosen_cai - _best_alt_cai, 6)
                _collector.record_codon_decision(CodonDecision(
                    position=_pos,
                    amino_acid=_aa,
                    original_codon=None,  # de novo optimization
                    chosen_codon=_codon,
                    alternatives_considered=_alts,
                    constraint_reason="maximize_cai",
                    confidence=1.0,
                    cai_impact=_cai_impact,
                ))
            # ── Record one ConstraintDecision per predicate (mirror slow
            # path pipeline_core.py:1437-1465) ──
            for _pr in pred_results:
                _action = "satisfied" if _pr.passed else "conflicted"
                _details = _pr.details or f"Predicate {_pr.predicate}"
                _collector.record_constraint_decision(ConstraintDecision(
                    constraint_name=_pr.predicate,
                    constraint_type="hard",
                    action_taken=_action,
                    positions_affected=[],
                    tradeoff_description=_details,
                    impact_on_cai=0.0,
                ))
            decision_trail = _collector.finalize(
                output_dna=seq, cai=cai_val, gc=gc_val,
            )
        except Exception as _prov_exc:
            logger.debug(
                "Decision-trail finalization failed on fast path: %s: %s",
                type(_prov_exc).__name__, _prov_exc, exc_info=True,
            )
            decision_trail = None

    # ── UTR suggestions (when include_utr=True; mirror slow-path
    # pipeline_core.py:1420-1431) ──
    utr5_seq: str | None = None
    utr3_seq: str | None = None
    utr_score5: float | None = None
    utr_score3: float | None = None
    if include_utr:
        from biocompiler.expression.utr_models import (
            suggest_5utr, suggest_3utr, score_5utr, score_3utr,
        )
        try:
            utr5_seq = suggest_5utr(organism)
            utr_score5 = score_5utr(utr5_seq, organism)
        except ValueError:
            logger.debug("No 5' UTR suggestion for organism '%s'", organism)
        try:
            utr3_seq = suggest_3utr(organism)
            utr_score3 = score_3utr(utr3_seq, organism)
        except ValueError:
            logger.debug("No 3' UTR suggestion for organism '%s'", organism)

    # ── Populate result fields in place ──
    result.predicate_results = pred_results
    result.certificate_text = cert_text
    result.failed_predicates = failed
    result.satisfied_predicates = satisfied
    result.provenance = provenance_record
    result.decision_trail = decision_trail
    result.uncertain_predicate_count = uncertain_count
    result.fallback_used = False             # integrated optimizer never mutagenizes today
    result.aa_substitutions = []
    result.mutagenesis_applied = False
    result.suggested_5utr = utr5_seq
    result.suggested_3utr = utr3_seq
    result.utr_score_5 = utr_score5
    result.utr_score_3 = utr_score3
    # mrna_stability_score / destabilizing_motifs_removed / stability_improvement
    # remain None / 0 / None — the integrated optimizer doesn't run the
    # mRNA-stability pass (callers wanting it should pass use_integrated=False).
    # codon_pair_bias and objective_score remain None for the same reason.

    # ── Strict-mode enforcement (mirror slow-path pipeline_core.py:1538-1542) ──
    # IMPORTANT: must run AFTER all fields are populated so the
    # ``partial_result`` carried by OptimizationConstraintError has a
    # populated certificate + predicate_results (per W2-c spec §4.1).
    if strict_mode and result.failed_predicates:
        raise OptimizationConstraintError(
            failed_predicates=result.failed_predicates,
            partial_result=result,
        )
