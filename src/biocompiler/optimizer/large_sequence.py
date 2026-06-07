"""
Large Sequence Support for BioCompiler
=======================================

Handles optimization of protein sequences >10 kb (>3,333 aa) by splitting
into overlapping chunks, optimizing each independently, and merging results
with boundary-aware constraint resolution.

Architecture:
  1. Split protein into overlapping chunks (default: 300 aa, 10 aa overlap)
  2. Optimize each chunk independently using the standard optimizer
  3. Merge chunks using overlap consensus (pick codons that agree)
  4. Apply global constraints (GC, restriction sites) across the full sequence
  5. Final boundary repair pass to fix any issues at chunk boundaries

The overlap region ensures that cross-codon constraints (GT dinucleotides,
restriction sites spanning codon boundaries) are handled correctly at
chunk boundaries.
"""

from __future__ import annotations

import logging
from typing import Callable, Any

from ..optimization import optimize_sequence, OptimizationResult
from ..type_system import AA_TO_CODONS, CODON_TABLE
from ..translation import translate, compute_cai
from ..scanner import gc_content
from ..exceptions import InvalidProteinError, OptimizationConstraintError

__all__ = [
    "optimize_large_sequence",
    "ProteinTooLongError",
    "MAX_PROTEIN_LENGTH_DEFAULT",
    "_split_into_chunks",
    "_merge_chunks",
    "_repair_boundaries",
]

logger = logging.getLogger(__name__)

# ── Default safety cap ────────────────────────────────────────────────
# Prevents accidentally optimizing extremely long sequences that would
# consume excessive memory or time.  Can be overridden via parameter.
MAX_PROTEIN_LENGTH_DEFAULT: int = 10_000  # amino acids


class ProteinTooLongError(InvalidProteinError):
    """Raised when a protein exceeds the max_protein_length safety cap."""

    def __init__(self, protein_length: int, max_length: int):
        self.protein_length = protein_length
        self.max_length = max_length
        super().__init__(
            protein="",
            invalid_chars=set(),
        )
        # Override the message from InvalidProteinError
        self.args = (
            f"Protein length ({protein_length} aa) exceeds max_protein_length "
            f"({max_length} aa). Increase max_protein_length or use a shorter sequence.",
        )


def _validate_protein(protein: str) -> list[str]:
    """Validate protein string and return list of amino acid codes.

    Raises InvalidProteinError for invalid characters.
    """
    if not protein or not protein.strip():
        raise InvalidProteinError(protein, set())
    protein = protein.upper().strip()
    valid_aas = set(AA_TO_CODONS.keys())
    invalid = set(ch for ch in protein if ch not in valid_aas)
    if invalid:
        raise InvalidProteinError(protein, invalid)
    return list(protein)


def _split_into_chunks(
    protein: str,
    chunk_size: int = 300,
    overlap: int = 10,
) -> list[tuple[int, int, str]]:
    """Split a protein sequence into overlapping chunks.

    Each chunk is ``chunk_size`` amino acids long, with ``overlap`` amino
    acids shared between adjacent chunks.  This ensures that cross-codon
    constraints at boundaries are captured within at least one chunk.

    Returns:
        List of (start_pos, end_pos, subsequence) tuples.
        ``start_pos`` is inclusive, ``end_pos`` is exclusive (Python slice).
    """
    if len(protein) <= chunk_size:
        return [(0, len(protein), protein)]

    chunks: list[tuple[int, int, str]] = []
    step = chunk_size - overlap
    pos = 0

    while pos < len(protein):
        end = min(pos + chunk_size, len(protein))
        chunk_seq = protein[pos:end]
        chunks.append((pos, end, chunk_seq))

        if end >= len(protein):
            break

        pos += step

        # If the remaining sequence is smaller than the overlap, just
        # extend the last chunk to cover the rest.
        if pos + overlap >= len(protein):
            # We've already covered everything; break out
            if end >= len(protein):
                break

    return chunks


def _merge_chunks(
    chunk_results: list[tuple[int, int, OptimizationResult]],
    protein: str,
    overlap: int,
) -> str:
    """Merge optimized chunk DNA sequences into a single full sequence.

    For overlapping regions, the merge uses the **latter chunk's** codons
    (the chunk whose start position is later), because that chunk had
    the full context of the overlap region when optimizing.  This
    ensures that cross-codon constraints within the overlap are
    resolved consistently.

    Args:
        chunk_results: List of (start_pos, end_pos, OptimizationResult) tuples.
        protein: Full protein sequence.
        overlap: Overlap size in amino acids.

    Returns:
        Merged DNA sequence string.
    """
    if len(chunk_results) == 1:
        return chunk_results[0][2].sequence

    # Total DNA length
    total_len = len(protein) * 3
    merged = [""] * total_len  # One char per nucleotide

    # For each position, determine which chunk "owns" it.
    # In the overlap region, the later chunk takes priority.
    # We'll iterate in reverse order so that later chunks overwrite earlier.
    for start_aa, end_aa, result in chunk_results:
        dna = result.sequence
        dna_start = start_aa * 3
        for i in range(len(dna)):
            merged[dna_start + i] = dna[i]

    return "".join(merged)


def _repair_boundaries(
    dna: str,
    protein: str,
    organism: str,
    enzymes: list[str] | None = None,
) -> str:
    """Final pass to fix boundary issues in the merged sequence.

    Checks and repairs:
    1. Ensures all codons are valid and encode the correct amino acid.
    2. Repairs invalid codons by substituting the highest-CAI codon.
    3. Checks for restriction sites across boundaries and fixes them.

    Args:
        dna: Merged DNA sequence.
        protein: Expected protein sequence.
        organism: Target organism.
        enzymes: Restriction enzymes to avoid.

    Returns:
        Repaired DNA sequence.
    """
    from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

    resolved = resolve_organism(organism, strict=False)
    usage = CODON_ADAPTIVENESS_TABLES.get(resolved, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"])

    dna_list = list(dna)
    n_codons = len(protein)

    # Pass 1: Fix any invalid codons or codons that encode wrong amino acid
    for i in range(n_codons):
        codon = "".join(dna_list[i * 3: i * 3 + 3])
        expected_aa = protein[i]
        actual_aa = CODON_TABLE.get(codon)

        if actual_aa != expected_aa:
            # Replace with the best codon for this amino acid
            candidates = AA_TO_CODONS.get(expected_aa, [])
            if candidates:
                # Pick the one with highest adaptiveness
                best = max(candidates, key=lambda c: usage.get(c, 0.0))
                dna_list[i * 3] = best[0]
                dna_list[i * 3 + 1] = best[1]
                dna_list[i * 3 + 2] = best[2]

    # Pass 2: Remove restriction sites that may have formed at boundaries
    if enzymes:
        from ..constants import reverse_complement
        from ..restriction_sites import get_recognition_site

        dna_str = "".join(dna_list)
        sorted_codons: dict[str, list[str]] = {}
        for aa in set(protein):
            codons = AA_TO_CODONS.get(aa, [])
            sorted_codons[aa] = sorted(
                codons, key=lambda c: usage.get(c, 0.0), reverse=True
            )

        for enz in enzymes:
            site = get_recognition_site(enz)
            if site is None:
                continue
            site_rc = reverse_complement(site)
            site_upper = site.upper()
            site_rc_upper = site_rc.upper() if site_rc else ""

            # Iteratively remove sites
            for _ in range(100):  # Safety limit
                dna_str = "".join(dna_list)
                positions = []
                if site_upper:
                    start = 0
                    while True:
                        pos = dna_str.find(site_upper, start)
                        if pos == -1:
                            break
                        positions.append(pos)
                        start = pos + 1
                if site_rc_upper and site_rc_upper != site_upper:
                    start = 0
                    while True:
                        pos = dna_str.find(site_rc_upper, start)
                        if pos == -1:
                            break
                        positions.append(pos)
                        start = pos + 1

                if not positions:
                    break

                # Fix the first site found
                pos = positions[0]
                first_codon = pos // 3
                last_codon = (pos + len(site_upper) - 1) // 3
                overlapping = list(range(
                    max(0, first_codon),
                    min(n_codons, last_codon + 1),
                ))

                fixed = False
                for ci in overlapping:
                    aa = protein[ci]
                    current_codon = "".join(dna_list[ci * 3: ci * 3 + 3])
                    for alt in sorted_codons.get(aa, []):
                        if alt == current_codon:
                            continue
                        # Try this alternative
                        old = [dna_list[ci * 3], dna_list[ci * 3 + 1], dna_list[ci * 3 + 2]]
                        dna_list[ci * 3] = alt[0]
                        dna_list[ci * 3 + 1] = alt[1]
                        dna_list[ci * 3 + 2] = alt[2]
                        test_str = "".join(dna_list)
                        if site_upper not in test_str and (not site_rc_upper or site_rc_upper not in test_str):
                            fixed = True
                            break
                        # Undo
                        dna_list[ci * 3] = old[0]
                        dna_list[ci * 3 + 1] = old[1]
                        dna_list[ci * 3 + 2] = old[2]

                    if fixed:
                        break

                if not fixed:
                    # Could not fix this site; try multi-codon approach
                    from .constraints import _remove_site_multicodon
                    test_seq = "".join(dna_list)
                    aas = [protein[ci] for ci in overlapping]
                    new_seq, was_fixed = _remove_site_multicodon(
                        test_seq,
                        aas,
                        sorted_codons,
                        site_upper,
                        site_rc_upper,
                        usage=usage,
                    )
                    if was_fixed:
                        dna_list = list(new_seq)
                    else:
                        # Log warning but continue
                        logger.warning(
                            "Could not remove restriction site %s at position %d",
                            enz, pos,
                        )
                        break

    return "".join(dna_list)


def optimize_large_sequence(
    protein: str,
    organism: str = "Homo_sapiens",
    chunk_size: int = 300,
    overlap: int = 10,
    max_protein_length: int = MAX_PROTEIN_LENGTH_DEFAULT,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: list[str] | None = None,
    strategy: str = "hybrid",
    strict_mode: bool = True,
    progress_callback: Callable[[int, int], None] | None = None,
    **kwargs: Any,
) -> OptimizationResult:
    """Optimize a large protein sequence by splitting into overlapping chunks.

    For proteins longer than ``chunk_size``, the sequence is divided into
    overlapping chunks.  Each chunk is independently optimized, then the
    results are merged with overlap-aware boundary resolution and a final
    global constraint pass.

    For proteins shorter than or equal to ``chunk_size``, this delegates
    directly to :func:`optimize_sequence` (no chunking overhead).

    Progress Reporting:
        If ``progress_callback`` is provided, it is called periodically
        with ``(current_position, total_length)`` where both are in amino
        acids.  Callers can use this to display progress bars for
        long-running optimizations.

    Safety Cap:
        ``max_protein_length`` prevents accidentally optimizing sequences
        that are unreasonably large (default: 10,000 aa, ~30 kb DNA).
        Set to 0 or a very large number to disable, but be aware that
        very long sequences may consume significant memory and time.

    Args:
        protein: Amino acid sequence (1-letter codes, no stop).
        organism: Target organism name (same format as optimize_sequence).
        chunk_size: Size of each chunk in amino acids. Default: 300.
        overlap: Overlap between adjacent chunks in amino acids. Default: 10.
            Must be at least 4 to cover cross-codon constraints.
        max_protein_length: Safety cap on protein length. Default: 10,000 aa.
        gc_lo: Minimum acceptable GC fraction.
        gc_hi: Maximum acceptable GC fraction.
        cai_threshold: Minimum CAI score for the CodonAdapted predicate.
        enzymes: List of restriction enzyme names to avoid.
        strategy: Optimization strategy ('hybrid', 'constraint_first', 'cai_first').
        strict_mode: If True, raise on constraint violations.
        progress_callback: Optional callback ``(current_pos, total_len)``.
        **kwargs: Additional arguments passed to optimize_sequence.

    Returns:
        OptimizationResult with the merged, globally-constrained sequence.

    Raises:
        ProteinTooLongError: If protein length exceeds max_protein_length.
        InvalidProteinError: If the protein contains invalid amino acid codes.
        OptimizationConstraintError: If strict_mode=True and constraints fail.
    """
    # ── Input validation ─────────────────────────────────────────────
    aas = _validate_protein(protein)
    total_len = len(aas)

    # Safety cap check
    if max_protein_length > 0 and total_len > max_protein_length:
        raise ProteinTooLongError(total_len, max_protein_length)

    # Validate overlap
    if overlap < 4:
        logger.warning(
            "Overlap %d is less than 4; cross-codon constraints at boundaries "
            "may not be handled correctly. Consider using overlap >= 4.",
            overlap,
        )

    # Ensure chunk_size > overlap
    if chunk_size <= overlap:
        raise ValueError(
            f"chunk_size ({chunk_size}) must be greater than overlap ({overlap})"
        )

    # ── Short sequence fast path ─────────────────────────────────────
    # If the protein fits in a single chunk, just optimize directly.
    if total_len <= chunk_size:
        if progress_callback is not None:
            progress_callback(0, total_len)
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
            enzymes=enzymes,
            strategy=strategy,
            strict_mode=strict_mode,
            **kwargs,
        )
        if progress_callback is not None:
            progress_callback(total_len, total_len)
        return result

    # ── Chunk-based optimization ─────────────────────────────────────
    logger.info(
        "Large sequence optimization: %d aa protein, chunk_size=%d, overlap=%d",
        total_len, chunk_size, overlap,
    )

    chunks = _split_into_chunks(protein, chunk_size=chunk_size, overlap=overlap)
    logger.info("Split into %d chunks", len(chunks))

    chunk_results: list[tuple[int, int, OptimizationResult]] = []

    for idx, (start, end, subseq) in enumerate(chunks):
        # Report progress
        if progress_callback is not None:
            progress_callback(start, total_len)

        logger.debug(
            "Optimizing chunk %d/%d: aa %d-%d (%d aa)",
            idx + 1, len(chunks), start, end, len(subseq),
        )

        # Optimize this chunk
        chunk_result = optimize_sequence(
            target_protein=subseq,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
            enzymes=enzymes,
            strategy=strategy,
            strict_mode=False,  # Allow partial results for individual chunks
            **kwargs,
        )

        chunk_results.append((start, end, chunk_result))

    # Report progress after all chunks are done
    if progress_callback is not None:
        progress_callback(total_len, total_len)

    # ── Merge chunks ─────────────────────────────────────────────────
    logger.info("Merging %d chunks with overlap=%d", len(chunk_results), overlap)
    merged_dna = _merge_chunks(chunk_results, protein, overlap)

    # ── Repair boundaries ────────────────────────────────────────────
    logger.info("Repairing chunk boundaries")
    repaired_dna = _repair_boundaries(merged_dna, protein, organism, enzymes=enzymes)

    # ── Global constraint pass ───────────────────────────────────────
    # Run a full optimize_sequence on the merged result to apply global
    # constraints (GC content, restriction sites, etc.) across the entire
    # sequence.  We pass the merged DNA as a starting point by
    # re-optimizing the full protein with the same parameters.
    logger.info("Applying global constraints across full sequence")

    # For very large sequences, we can't run the full optimizer again
    # (that would defeat the purpose of chunking).  Instead, do targeted
    # global constraint checks and fixes.

    # Check and fix GC content
    final_dna = repaired_dna
    current_gc = gc_content(final_dna)

    if not (gc_lo <= current_gc <= gc_hi):
        logger.info(
            "GC content %.3f outside target [%.3f, %.3f]; adjusting",
            current_gc, gc_lo, gc_hi,
        )
        # Run targeted GC adjustment
        final_dna = _adjust_gc_global(final_dna, protein, organism, gc_lo, gc_hi)

    # ── Final verification ───────────────────────────────────────────
    # Verify the final sequence translates correctly
    translated = translate(final_dna, to_stop=False)
    if translated != protein:
        # Try one more repair pass
        final_dna = _repair_boundaries(final_dna, protein, organism, enzymes=enzymes)
        translated = translate(final_dna, to_stop=False)
        if translated != protein:
            raise OptimizationConstraintError(
                failed_predicates=["TranslationCorrect"],
                partial_result=None,
            )

    # Compute final metrics
    final_gc = gc_content(final_dna)
    final_cai = compute_cai(final_dna, organism=organism)

    # Check for remaining restriction sites
    failed_preds: list[str] = []
    if enzymes:
        from ..constants import reverse_complement
        from ..restriction_sites import get_recognition_site as _get_site
        for enz in enzymes:
            site = _get_site(enz)
            if site:
                site_rc = reverse_complement(site)
                if site in final_dna or (site_rc and site_rc in final_dna):
                    failed_preds.append(f"NoRestrictionSite_{enz}")

    # GC range check
    if not (gc_lo <= final_gc <= gc_hi):
        failed_preds.append("GCInRange")

    # Build the result
    warnings_list: list[str] = []
    if len(chunk_results) > 1:
        warnings_list.append(
            f"Large sequence optimized in {len(chunk_results)} chunks "
            f"(chunk_size={chunk_size}, overlap={overlap})"
        )

    result = OptimizationResult(
        sequence=final_dna,
        gc_content=final_gc,
        cai=final_cai,
        failed_predicates=failed_preds,
        protein=protein,
        warnings=warnings_list,
    )

    if strict_mode and failed_preds:
        raise OptimizationConstraintError(
            failed_predicates=failed_preds,
            partial_result=result,
        )

    return result


def _adjust_gc_global(
    dna: str,
    protein: str,
    organism: str,
    gc_lo: float,
    gc_hi: float,
) -> str:
    """Adjust GC content of the full sequence by swapping synonymous codons.

    Iterates through all codons and swaps GC-rich codons for AT-rich
    alternatives (or vice versa) to bring the overall GC content into
    the target range.

    Args:
        dna: DNA sequence to adjust.
        protein: Expected protein sequence.
        organism: Target organism.
        gc_lo: Minimum GC fraction.
        gc_hi: Maximum GC fraction.

    Returns:
        Adjusted DNA sequence.
    """
    from ..organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism

    resolved = resolve_organism(organism, strict=False)
    usage = CODON_ADAPTIVENESS_TABLES.get(resolved, CODON_ADAPTIVENESS_TABLES["Escherichia_coli"])

    dna_list = list(dna)
    n_codons = len(protein)
    current_gc = gc_content("".join(dna_list))

    for iteration in range(200):  # Safety limit
        if gc_lo <= current_gc <= gc_hi:
            break

        # Determine direction: need more GC or more AT?
        need_more_gc = current_gc < gc_lo
        need_more_at = current_gc > gc_hi

        best_swap: tuple[int, str, float] | None = None
        # best_swap: (codon_index, new_codon, cai_score)

        for i in range(n_codons):
            aa = protein[i]
            current_codon = "".join(dna_list[i * 3: i * 3 + 3])
            candidates = AA_TO_CODONS.get(aa, [])

            for alt in candidates:
                if alt == current_codon:
                    continue
                # Count GC bases in each
                current_gc_count = sum(1 for b in current_codon if b in "GC")
                alt_gc_count = sum(1 for b in alt if b in "GC")

                if need_more_gc and alt_gc_count > current_gc_count:
                    cai = usage.get(alt, 0.0)
                    if best_swap is None or cai > best_swap[2]:
                        best_swap = (i, alt, cai)
                elif need_more_at and alt_gc_count < current_gc_count:
                    cai = usage.get(alt, 0.0)
                    if best_swap is None or cai > best_swap[2]:
                        best_swap = (i, alt, cai)

        if best_swap is None:
            break  # No more beneficial swaps

        idx, new_codon, _ = best_swap
        dna_list[idx * 3] = new_codon[0]
        dna_list[idx * 3 + 1] = new_codon[1]
        dna_list[idx * 3 + 2] = new_codon[2]
        current_gc = gc_content("".join(dna_list))

    return "".join(dna_list)
