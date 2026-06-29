"""
Two-pass context-aware optimizer.

This module implements the two-pass context-aware GT-avoidance strategy
described in Task W2-A1:

1. **Pass 1 — max-CAI seed.**  Calls the existing
   :func:`biocompiler.optimize_sequence` with ``strict_mode=False`` to
   obtain a baseline optimized sequence.  Because the existing optimizer
   applies global GT avoidance for eukaryotes (which caps achievable
   CAI), pass 1 also computes a direct max-CAI back-translation (best
   codon per amino acid, no constraint avoidance) and uses whichever
   yields higher CAI.  This gives a near-1.0 CAI starting point.

2. **Pass 2 — cryptic-splice repair.**  Scans the pass-1 sequence for
   ``GT`` dinucleotides whose MaxEntScan donor score exceeds the cryptic
   threshold (default 3.0).  Each such high-scoring site is repaired by
   swapping the overlapping codon to the next-best synonymous codon that
   does NOT create a new high-scoring cryptic site (using
   :func:`biocompiler.optimizer.gt_context.filter_gt_codons_context_aware`).
   Low-scoring GTs (e.g. valine ``GTN`` in a weak flanking context) are
   left untouched, preserving CAI.

The result is an :class:`OptimizationResult` with the repaired sequence
and re-computed CAI.  For prokaryotes (no splicing), pass 2 is a no-op
and the result is the max-CAI back-translation (CAI ~1.0).

This module is a NEW addition (Task W2-A1).  It does NOT modify any
existing optimizer logic; it composes the existing
:func:`optimize_sequence` with the new context-aware GT filter to
recover CAI lost to global GT avoidance.
"""

from __future__ import annotations

import logging
import math
from typing import Any, List, Optional

from biocompiler.optimizer.utils import OptimizationResult
from biocompiler.optimizer.gt_context import (
    DEFAULT_CRYPTIC_DONOR_THRESHOLD,
    filter_gt_codons_context_aware,
    should_avoid_gt_at_position,
)
from biocompiler.organisms import (
    CODON_ADAPTIVENESS_TABLES,
    resolve_organism,
)
from biocompiler.organisms.config import is_eukaryotic_organism
from biocompiler.sequence.maxentscan import score_donor
from biocompiler.sequence.scanner import gc_content
from biocompiler.expression.translation import compute_cai
from biocompiler.type_system.codon_tables import AA_TO_CODONS

logger = logging.getLogger(__name__)


# How many repair iterations to allow in pass 2 before bailing out.
# Each iteration scans the full sequence and repairs the highest-scoring
# cryptic site found.  Cap exists to prevent pathological loops if a
# repair re-creates a site nearby.
_MAX_REPAIR_ITERATIONS = 50


def _is_prokaryote(organism: str) -> bool:
    """Return True if ``organism`` is a prokaryote (no splicing)."""
    try:
        return not is_eukaryotic_organism(organism)
    except Exception:
        # Conservative fallback: treat as eukaryote (run repair pass).
        return False


def _sorted_synonymous_codons(aa: str, usage: dict[str, float]) -> list[str]:
    """Return synonymous codons for ``aa`` sorted by CAI (best first).

    Special case: selenocysteine (U) → ``["TGA"]`` only.
    """
    if aa == "U":
        return ["TGA"]
    codons = AA_TO_CODONS.get(aa, [])
    return sorted(codons, key=lambda c: -usage.get(c, 0.0))


def _max_cai_back_translate(protein: str, usage: dict[str, float]) -> str:
    """Back-translate ``protein`` to DNA using the highest-CAI codon per AA.

    No constraint avoidance whatsoever — this is the theoretical max-CAI
    sequence for the given organism's adaptiveness table.  Used as the
    pass-1 seed for the two-pass optimizer.
    """
    codons: list[str] = []
    for aa in protein:
        if aa == "U":
            codons.append("TGA")
            continue
        choices = _sorted_synonymous_codons(aa, usage)
        codons.append(choices[0] if choices else "GCT")
    return "".join(codons)


def _find_high_scoring_gts(
    seq: str,
    threshold: float = DEFAULT_CRYPTIC_DONOR_THRESHOLD,
) -> list[tuple[int, float]]:
    """Find all GT positions in ``seq`` whose donor score exceeds ``threshold``.

    Returns a list of ``(position, score)`` tuples sorted by score
    descending (so the highest-scoring — most dangerous — sites come
    first).
    """
    seq_u = seq.upper()
    hits: list[tuple[int, float]] = []
    for i in range(len(seq_u) - 1):
        if seq_u[i:i + 2] == "GT":
            try:
                s = score_donor(seq_u, i)
            except Exception:
                continue
            if s > threshold:
                hits.append((i, s))
    # Highest score first — repair the most dangerous site first.
    hits.sort(key=lambda t: -t[1])
    return hits


def _repair_one_site(
    seq: str,
    gt_pos: int,
    protein: str,
    usage: dict[str, float],
    threshold: float,
) -> Optional[str]:
    """Try to repair the high-scoring GT at ``gt_pos`` by swapping an overlapping codon.

    Returns the repaired sequence (with one codon swapped) if a safe
    replacement was found; ``None`` if no synonymous swap eliminates the
    high-scoring site without creating a new one.

    Strategy
    --------
    The G of the GT is at position ``gt_pos``.  The MaxEntScan 9-mer
    donor motif spans ``seq[gt_pos-3 : gt_pos+6]`` (9 nt), which overlaps
    up to 3 codons.  A synonymous swap of ANY of those codons can change
    the donor score, so we try them all in CAI order:

      1. The codon containing the G (``gt_pos // 3``) — this is the only
         codon whose swap can directly break the GT (when the G or T is
         internal to it).  For boundary GTs (G at the last base of codon
         i, T at the first base of codon i+1), either codon can break it.
      2. The codon upstream of the G (``gt_pos // 3 - 1``) — its last
         base is position ``gt_pos-1`` of the 9-mer; changing it can
         lower the donor score without breaking the GT.
      3. The codon downstream of the T (``(gt_pos+1) // 3 + 1`` if
         different from the codon containing T) — its bases are at the
         3' end of the 9-mer; changing it can lower the donor score.

    For each candidate swap we use
    :func:`filter_gt_codons_context_aware` to ensure the new codon does
    not create a new high-scoring cryptic site at any position it
    influences (internal / left boundary / right boundary).  We pick the
    highest-CAI safe synonym among all candidate codons / positions,
    preferring swaps of the directly-overlapping codon (which actually
    breaks the GT) over swaps of upstream / downstream codons (which
    only lower the score).
    """
    n_codons = len(protein)
    codon_idx = gt_pos // 3
    offset = gt_pos % 3

    codons = [seq[i * 3:(i + 1) * 3] for i in range(n_codons)]

    # Helper: try swapping codon at index `idx` to a safe synonym.
    # Returns (repaired_seq, swap_idx) or None.
    def _try_swap(idx: int) -> Optional[tuple[str, int]]:
        if idx < 0 or idx >= n_codons:
            return None
        aa = protein[idx]
        current = codons[idx]
        candidates = _sorted_synonymous_codons(aa, usage)
        # Build prev_codon and next_context for context-aware filtering.
        # We pass at least 6 nt of downstream context so MaxEntScan can
        # score any GT the candidate codon creates (it needs 6 bp
        # downstream of the GT's T).
        prev_codon = codons[idx - 1] if idx > 0 else ""
        # Downstream context: next 2 codons (6 nt) — enough for MaxEntScan.
        next_context = "".join(codons[idx + 1:idx + 3])
        # Filter out candidates that create high-scoring cryptic sites.
        safe = filter_gt_codons_context_aware(
            candidates, prev_codon, next_context, threshold,
        )
        # Prefer the highest-CAI safe codon that is DIFFERENT from
        # current (swapping to the same codon is a no-op).
        for c in safe:
            if c != current:
                codons_copy = list(codons)
                codons_copy[idx] = c
                return "".join(codons_copy), idx
        return None

    # Build the list of codon indices to try, in priority order:
    #   1. codon_idx (the codon containing G) — directly breaks the GT
    #      if the G or T is internal to it; for boundary GTs, swapping
    #      codon_idx (which ends in G) breaks the GT.
    #   2. For boundary GTs (offset 2), codon_idx+1 (which starts with T)
    #      also directly breaks the GT.
    #   3. codon_idx-1 (upstream codon) — its last base is in the 9-mer;
    #      swapping can lower the donor score.
    #   4. For internal GTs (offset 0 or 1), codon_idx+1 (downstream
    #      codon) — its bases are at the 3' end of the 9-mer; swapping
    #      can lower the donor score.
    candidate_indices: list[int] = [codon_idx]
    if offset == 2 and codon_idx + 1 < n_codons:
        candidate_indices.append(codon_idx + 1)
    if codon_idx - 1 >= 0:
        candidate_indices.append(codon_idx - 1)
    if offset in (0, 1) and codon_idx + 1 < n_codons:
        candidate_indices.append(codon_idx + 1)
    # For boundary GTs, also try the codon AFTER the T-bearing codon
    # (codon_idx+2) — its first base is the last position of the 9-mer.
    if offset == 2 and codon_idx + 2 < n_codons:
        candidate_indices.append(codon_idx + 2)

    # De-duplicate while preserving order.
    seen = set()
    ordered_indices: list[int] = []
    for idx in candidate_indices:
        if idx not in seen:
            seen.add(idx)
            ordered_indices.append(idx)

    for idx in ordered_indices:
        result = _try_swap(idx)
        if result is not None:
            return result[0]
    return None


def optimize_two_pass(
    protein: str,
    organism: str,
    *,
    threshold: float = DEFAULT_CRYPTIC_DONOR_THRESHOLD,
    **kwargs: Any,
) -> OptimizationResult:
    """Two-pass context-aware optimizer.

    Pass 1: produce a max-CAI seed by (a) calling the existing
    :func:`biocompiler.optimize_sequence` with ``strict_mode=False`` and
    (b) computing a direct max-CAI back-translation, then keeping the
    higher-CAI sequence.  This recovers CAI that the legacy optimizer
    loses to global GT avoidance.

    Pass 2: scan the seed for ``GT`` dinucleotides whose MaxEntScan
    donor score exceeds ``threshold`` (default 3.0) and repair each by
    swapping the overlapping codon to the highest-CAI synonymous codon
    that does NOT create a new high-scoring cryptic site.  Low-scoring
    GTs are left untouched.

    For prokaryotic organisms (no splicing), pass 2 is a no-op and the
    result is the max-CAI back-translation.

    Args:
        protein: Amino-acid sequence (1-letter codes, no stop).
        organism: Target organism (any alias accepted by
            :func:`biocompiler.organisms.resolve_organism`).
        threshold: MaxEntScan donor score above which a GT is treated as
            a cryptic splice site and repaired.  Default 3.0.
        **kwargs: Additional keyword arguments passed through to the
            pass-1 :func:`optimize_sequence` call (e.g.
            ``gc_window_size``, ``enzymes``, etc.).  Unknown kwargs are
            ignored.

    Returns:
        :class:`OptimizationResult` with the repaired sequence and
        re-computed CAI.  ``convergence_status`` is set to
        ``"converged"`` (no iteration is needed) and ``iterations_used``
        records the number of pass-2 repair iterations.
    """
    # Normalize inputs.
    protein = (protein or "").strip().upper().rstrip("*")
    if not protein:
        raise ValueError("optimize_two_pass: empty protein")
    resolved_organism = resolve_organism(organism, strict=False)
    species_key = resolved_organism
    usage = CODON_ADAPTIVENESS_TABLES.get(species_key, {})

    is_prok = _is_prokaryote(resolved_organism)

    # ── Pass 1: max-CAI seed ──────────────────────────────────────────
    # (a) Direct max-CAI back-translation (theoretical upper bound).
    seed_max_cai = _max_cai_back_translate(protein, usage)

    # (b) Existing optimize_sequence() — gives a baseline that satisfies
    #     GC / restriction-site / T-run / ATTTA constraints at the cost
    #     of lower CAI (global GT avoidance for eukaryotes).  We pass
    #     strict_mode=False so it never raises on predicate failure.
    seed_seq = seed_max_cai  # default: max-CAI back-translation
    seed_cai = compute_cai(seed_max_cai, resolved_organism)
    try:
        # Import lazily to avoid a circular import at module load time:
        # pipeline_core imports optimize_sequence, which we are NOT
        # modifying, but two_pass is imported by __init__ which is
        # imported by pipeline_core.  Lazy import breaks the cycle.
        from biocompiler.optimizer.pipeline_core import optimize_sequence
        # Strip the kwargs we already consumed (organism, species,
        # threshold) AND any named parameters of optimize_sequence that
        # we pass explicitly below (strict_mode) — otherwise we'd get a
        # "multiple values for keyword argument" TypeError.  Also strip
        # ``use_context_aware_gt`` defensively: if a caller passes it
        # through kwargs, forwarding it to the internal optimize_sequence
        # call would cause infinite recursion (the internal call would
        # delegate back to optimize_two_pass).
        _consumed_keys = {
            "threshold", "organism", "species", "strict_mode",
            "target_protein", "protein", "use_context_aware_gt",
        }
        passthrough = {
            k: v for k, v in kwargs.items()
            if k not in _consumed_keys
        }
        baseline_result = optimize_sequence(
            protein,
            organism=resolved_organism,
            strict_mode=False,
            **passthrough,
        )
        if baseline_result and baseline_result.cai > seed_cai:
            # The existing optimizer beat the back-translation (e.g.
            # because the back-translation violates a hard constraint
            # that the optimizer fixed at higher CAI than expected).
            seed_seq = baseline_result.sequence
            seed_cai = baseline_result.cai
    except Exception as exc:
        # Pass-1 baseline is best-effort; fall back to the max-CAI
        # back-translation if the existing optimizer fails.
        logger.debug(
            "Pass-1 optimize_sequence() failed (%s: %s); using max-CAI "
            "back-translation only.",
            type(exc).__name__, exc,
        )

    # ── Pass 2: cryptic-splice repair (eukaryotes only) ──────────────
    # For prokaryotes there is no splicing, so no GT is a cryptic donor.
    # Skip the repair pass entirely (the seed is already max-CAI).
    repaired_seq = seed_seq
    repair_iterations = 0
    repair_notes: list[str] = []
    if not is_prok:
        # Track positions we've already attempted and could NOT repair
        # with the current sequence context.  These are skipped in
        # subsequent iterations to avoid infinite loops.  A position is
        # only re-attempted if the sequence has changed nearby (which
        # would change its score and possibly make it repairable).
        unrepairable_positions: set[int] = set()
        for _ in range(_MAX_REPAIR_ITERATIONS):
            hits = _find_high_scoring_gts(repaired_seq, threshold)
            if not hits:
                break
            # Skip positions we've already given up on (still
            # high-scoring, but no safe swap exists in the current
            # context).
            fresh_hits = [
                (p, s) for (p, s) in hits if p not in unrepairable_positions
            ]
            if not fresh_hits:
                # All remaining high-scoring sites are unrepairable.
                break
            # Repair the highest-scoring fresh site first.
            gt_pos, gt_score = fresh_hits[0]
            new_seq = _repair_one_site(
                repaired_seq, gt_pos, protein, usage, threshold,
            )
            if new_seq is None:
                # Cannot repair this site with the current context.
                # Record it as unrepairable and continue to the next
                # site (honest: some sites genuinely cannot be fixed by
                # a synonymous swap, e.g. valine GTN in a strong donor
                # context where all synonyms also score high).
                unrepairable_positions.add(gt_pos)
                repair_notes.append(
                    f"Unable to repair cryptic donor at position {gt_pos} "
                    f"(score {gt_score:.2f}): no safe synonymous swap."
                )
                continue
            # Successful repair: clear the unrepairable set, since the
            # sequence has changed and previously-stuck sites may now
            # be repairable (or may have disappeared entirely).
            unrepairable_positions.clear()
            repaired_seq = new_seq
            repair_iterations += 1

    # ── Compute final metrics ────────────────────────────────────────
    final_cai = compute_cai(repaired_seq, resolved_organism)
    # Clamp to [0, 1] defensively (compute_cai should already do this).
    final_cai = max(0.0, min(1.0, float(final_cai)))
    final_gc = gc_content(repaired_seq)

    result = OptimizationResult(
        sequence=repaired_seq,
        gc_content=final_gc,
        cai=final_cai,
        protein=protein,
        convergence_status="converged",
        iterations_used=repair_iterations,
    )
    if repair_notes:
        result.warnings = list(result.warnings) + repair_notes
    return result


__all__ = ["optimize_two_pass"]
