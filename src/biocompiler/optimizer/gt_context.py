"""
Context-aware GT avoidance for the BioCompiler optimizer.

This module implements *context-aware* GT dinucleotide avoidance: instead
of avoiding every GT dinucleotide globally (which forces the optimizer to
discard otherwise-optimal codons like valine ``GTN``), it only avoids GTs
that would form a high-scoring cryptic 5' splice donor site as scored by
MaxEntScan (Yeo & Burge 2004).

Background
----------
The legacy optimizer globally avoids the ``GT`` dinucleotide in eukaryotic
sequences because ``GT`` is the canonical 5' splice donor consensus.  This
is overly conservative: the vast majority of internal GTs (e.g. valine
codons ``GTT/GTC/GTA/GTG``, cysteine ``TGT/TGC``) are NOT functional
splice donors because their flanking sequence context does not match the
extended 9-mer donor motif (positions -3..+6 around the GT) that
MaxEntScan models.  Avoiding them anyway forces the optimizer to use
lower-CAI synonymous codons, tanking CAI on human genes.

The fix is *context-aware*: only avoid a GT if its MaxEntScan donor score
exceeds the cryptic splice threshold (default 3.0).  Low-scoring GTs are
safe to keep, preserving CAI.

This module is a NEW addition (Task W2-A1).  It does NOT modify any
existing optimizer logic; it provides standalone helpers that the
two-pass optimizer (``optimizer.two_pass``) uses.

Public API
----------
- :func:`should_avoid_gt_at_position` — predicate: does the GT at ``pos``
  score above ``threshold``?
- :func:`filter_gt_codons_context_aware` — given a list of candidate
  codons, return those that do NOT create a high-scoring cryptic donor.
"""

from __future__ import annotations

from typing import List

from biocompiler.sequence.maxentscan import score_donor

# Default threshold for cryptic splice donor detection.
# MaxEntScan donor scores are log2-scale; published work (Yeo & Burge 2004)
# uses ~3.0 as the cutoff between real and cryptic/weak donors.  We keep
# this as the default but expose it as a parameter for callers that want
# stricter or more permissive filtering.
DEFAULT_CRYPTIC_DONOR_THRESHOLD: float = 3.0


def should_avoid_gt_at_position(
    seq: str,
    pos: int,
    threshold: float = DEFAULT_CRYPTIC_DONOR_THRESHOLD,
) -> bool:
    """Return ``True`` if a GT at ``pos`` would form a high-scoring cryptic donor.

    The check uses MaxEntScan's :func:`score_donor` (Yeo & Burge 2004) on
    the 9-mer ``seq[pos-3:pos+6]`` centered on the GT.  A score above
    ``threshold`` (default 3.0) means the GT looks like a functional 5'
    splice donor and should be avoided; a score at or below threshold
    means the GT is biologically inert and safe to keep.

    Conservative semantics
    ----------------------
    - If ``pos`` does not point at the ``G`` of a ``GT`` dinucleotide
      (e.g. ``seq[pos:pos+2] != "GT"``), returns ``False`` (no GT to
      avoid).
    - If MaxEntScan cannot score the GT (insufficient flanking sequence),
      :func:`score_donor` returns ``-20.0``, which is below any
      reasonable threshold — we return ``False`` (treat as safe).  This
      is conservative: we do NOT filter codons we cannot confidently
      score, so we never over-avoid.
    - If ``pos`` is out of range, returns ``False``.

    Args:
        seq: DNA sequence (upper or lower case; the function uppercases
            it internally).
        pos: Index of the ``G`` in a putative ``GT`` dinucleotide.
        threshold: MaxEntScan donor score above which a GT is considered
            a cryptic splice site.  Default 3.0 (matches
            ``CRYPTIC_SPLICE_THRESHOLD`` in :mod:`biocompiler.sequence.maxentscan`).

    Returns:
        ``True`` if the GT at ``pos`` should be avoided (high donor
        score); ``False`` otherwise (no GT, or low / unscoreable GT).

    Examples
    --------
    >>> # Internal GT in a valine codon with weak flanking → safe
    >>> should_avoid_gt_at_position("ATGGTGCACCTGACCCCCCC", 3)
    False
    >>> # GT in a strong donor context → avoid
    >>> strong = "CAGGTAAGT"  # canonical strong donor 9-mer
    >>> should_avoid_gt_at_position(strong, 3)
    True
    """
    if not seq or pos < 0:
        return False
    seq_u = seq.upper()
    if pos + 1 >= len(seq_u):
        return False
    if seq_u[pos:pos + 2] != "GT":
        return False
    try:
        score = score_donor(seq_u, pos)
    except Exception:
        # Be conservative on any MaxEntScan failure: do not filter.
        return False
    return score > threshold


def filter_gt_codons_context_aware(
    codon_choices: list[str],
    prev_codon: str,
    next_context: str,
    threshold: float = DEFAULT_CRYPTIC_DONOR_THRESHOLD,
) -> list[str]:
    """Filter candidate codons, keeping only those that don't create a high-scoring cryptic donor.

    For each candidate codon, build the local context
    ``prev_codon + codon + next_context`` and check every ``GT`` that the
    candidate codon *creates* (internally, at the left boundary with
    ``prev_codon``, or at the right boundary with ``next_context``).
    Codons that create a GT whose MaxEntScan donor score exceeds
    ``threshold`` are filtered out.  Codons that create only low-scoring
    GTs (or no GTs at all) are kept — this is the whole point of
    context-aware avoidance.

    The order of ``codon_choices`` is preserved in the output: callers
    are expected to pass codons already sorted by CAI (best first), so
    the first kept codon is the highest-CAI safe choice.

    Args:
        codon_choices: Candidate synonymous codons (typically CAI-sorted,
            best first).  May be empty.
        prev_codon: The codon placed immediately before the position
            being filled (upper or lower case; may be empty string if
            this is the first codon).
        next_context: The downstream nucleotide context starting
            immediately after the candidate codon (upper or lower case;
            may be empty).  At least 6 bp of downstream context is
            needed for MaxEntScan to score right-boundary and internal
            GTs; if less is provided, such GTs are treated as
            unscoreable (safe) — see :func:`should_avoid_gt_at_position`.
        threshold: MaxEntScan donor score above which a GT is considered
            a cryptic splice site.  Default 3.0.

    Returns:
        A new list containing the codons from ``codon_choices`` that do
        not create a high-scoring cryptic donor.  Order is preserved.
        If all codons create high-scoring sites, the list may be empty —
        callers should fall back to the highest-CAI candidate in that
        case (honest: cannot always fix).

    Examples
    --------
    >>> # Valine codons — all contain GT, but only those that create
    >>> # high-scoring donors are filtered out.
    >>> codons = ["GTG", "GTT", "GTC", "GTA"]
    >>> # In a weak donor context, all valine codons are kept.
    >>> safe = filter_gt_codons_context_aware(codons, "CCC", "CCC")
    >>> set(safe) == set(codons)
    True
    """
    if not codon_choices:
        return []

    prev_u = prev_codon.upper() if prev_codon else ""
    next_u = next_context.upper() if next_context else ""
    codon_len = 3  # canonical codon length; codon_choices are 3-nt strings
    codon_start = len(prev_u)

    safe: list[str] = []
    for codon in codon_choices:
        codon_u = codon.upper()
        if len(codon_u) != codon_len:
            # Not a canonical codon — pass through unchanged (don't guess).
            safe.append(codon)
            continue
        # Build the local context for MaxEntScan scoring.
        local = prev_u + codon_u + next_u
        # Candidate codon occupies [codon_start, codon_start+3) in `local`.
        codon_end = codon_start + codon_len

        # Positions where a GT could be created or influenced by this codon:
        #   - Internal positions: [codon_start, codon_end - 1) — the G of
        #     an internal GT must lie at codon_start or codon_start+1.
        #   - Left boundary: codon_start - 1 — the G is the last base of
        #     prev_codon, T is the first base of this codon.
        #   - Right boundary: codon_end - 1 — the G is the last base of
        #     this codon, T is the first base of next_context.
        candidate_gt_positions: list[int] = []
        # Internal GTs (G at codon_start, codon_start+1)
        for offset in (0, 1):
            p = codon_start + offset
            if p + 1 < codon_end and local[p:p + 2] == "GT":
                candidate_gt_positions.append(p)
        # Left boundary GT (G = prev_codon[-1], T = codon[0])
        if codon_start > 0:
            p = codon_start - 1
            if p + 1 < len(local) and local[p:p + 2] == "GT":
                candidate_gt_positions.append(p)
        # Right boundary GT (G = codon[-1], T = next_context[0])
        if codon_end > 0 and codon_end < len(local):
            p = codon_end - 1
            if p + 1 < len(local) and local[p:p + 2] == "GT":
                candidate_gt_positions.append(p)

        # Check each GT position; if any has a high donor score, filter out.
        is_safe = True
        for p in candidate_gt_positions:
            if should_avoid_gt_at_position(local, p, threshold):
                is_safe = False
                break
        if is_safe:
            safe.append(codon)

    return safe


__all__ = [
    "should_avoid_gt_at_position",
    "filter_gt_codons_context_aware",
    "DEFAULT_CRYPTIC_DONOR_THRESHOLD",
]
