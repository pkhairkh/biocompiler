"""
BLAST match avoidance for orthogonal sequence design.

Unlike biosecurity screening (which checks for hazardous pathogen/toxin sequences),
BLAST avoidance ensures sequence orthogonality — that the optimized DNA does not
have significant matches against a reference genome or database. This is critical
for designing non-interfering genetic circuits and synthetic biology parts.

Supports two modes:
1. Local BLAST+ (if installed): blastn-short for precise alignment checking
2. K-mer based: Fast approximate matching using k-mer overlap with a reference
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass

from ..type_system import CODON_TABLE, AA_TO_CODONS

logger = logging.getLogger(__name__)

__all__ = [
    "BLASTMatch",
    "check_blast_matches",
    "check_blast_matches_cli",
    "avoid_blast_matches",
    "check_kmer_overlap",
    "eliminate_kmer_overlaps",
]


@dataclass
class BLASTMatch:
    """A single BLAST match against a reference sequence or database."""

    query_start: int
    query_end: int
    reference_start: int
    reference_end: int
    identity: float
    e_value: float
    reference_id: str = ""


def check_blast_matches(
    sequence: str,
    reference_sequences: list[str] | None = None,
    reference_file: str | None = None,
    word_size: int = 15,
    e_value: float = 1.0,
    min_identity: float = 80.0,
) -> list[BLASTMatch]:
    """Check for BLAST matches between the sequence and references.

    If ``reference_sequences`` is provided, use k-mer matching (always
    available, no external tools required).  If ``reference_file`` is
    provided, try BLAST+ CLI for precise alignment checking.

    Args:
        sequence: DNA sequence to check (uppercase).
        reference_sequences: List of reference DNA sequences for k-mer matching.
        reference_file: Path to a FASTA file for BLAST+ search.
        word_size: K-mer size for k-mer matching, or BLAST word size.
        e_value: E-value threshold for BLAST+.
        min_identity: Minimum percent identity to report (0-100).

    Returns:
        List of BLASTMatch objects for significant matches found.
    """
    sequence = sequence.upper()
    matches: list[BLASTMatch] = []

    # K-mer based matching (always available)
    if reference_sequences:
        overlaps = check_kmer_overlap(sequence, reference_sequences, k=word_size)
        for start, length, kmer in overlaps:
            matches.append(BLASTMatch(
                query_start=start,
                query_end=start + length,
                reference_start=0,  # k-mer mode doesn't track ref position precisely
                reference_end=length,
                identity=100.0,  # exact k-mer match
                e_value=0.0,
                reference_id="kmer_overlap",
            ))

    # BLAST+ CLI mode (if available and reference_file provided)
    if reference_file:
        blast_matches = _run_blastn_short(sequence, reference_file, word_size, e_value, min_identity)
        matches.extend(blast_matches)

    return matches


def avoid_blast_matches(
    sequence: str,
    protein: str,
    organism: str,
    reference_sequences: list[str] | None = None,
    reference_file: str | None = None,
    word_size: int = 15,
    max_iterations: int = 100,
) -> str:
    """Iteratively replace codons at match positions to eliminate BLAST matches.

    Preserves the protein sequence by only using synonymous codon substitutions.

    Args:
        sequence: Current DNA sequence (uppercase).
        protein: Target protein sequence (amino acid string).
        organism: Target organism for codon preference.
        reference_sequences: Reference sequences for k-mer matching.
        reference_file: Path to a FASTA file for BLAST+ search.
        word_size: K-mer size for overlap detection.
        max_iterations: Maximum optimization iterations.

    Returns:
        Optimized DNA sequence with reduced BLAST matches.
    """
    sequence = sequence.upper()

    for iteration in range(max_iterations):
        # Check for overlaps
        overlaps: list[tuple[int, int, str]] = []
        if reference_sequences:
            overlaps = check_kmer_overlap(sequence, reference_sequences, k=word_size)

        # Also check BLAST+ matches if file provided
        blast_matches: list[BLASTMatch] = []
        if reference_file:
            blast_matches = check_blast_matches(
                sequence, reference_file=reference_file,
                word_size=word_size, min_identity=80.0,
            )

        if not overlaps and not blast_matches:
            break  # No more matches — done!

        # Collect positions that need fixing
        positions_to_fix: set[int] = set()
        for start, length, _ in overlaps:
            for pos in range(start, min(start + length, len(sequence))):
                positions_to_fix.add(pos)
        for m in blast_matches:
            for pos in range(m.query_start, min(m.query_end, len(sequence))):
                positions_to_fix.add(pos)

        # Try to fix by replacing codons at overlap positions
        fixed_any = False
        for pos in sorted(positions_to_fix):
            codon_idx = pos // 3
            codon_start = codon_idx * 3
            if codon_start + 3 > len(sequence):
                continue

            current_codon = sequence[codon_start:codon_start + 3]
            aa = CODON_TABLE.get(current_codon)
            if aa is None or aa == "*":
                continue

            # Try synonymous alternatives
            alternatives = AA_TO_CODONS.get(aa, [])
            for alt in alternatives:
                if alt == current_codon:
                    continue
                # Build trial sequence
                trial = sequence[:codon_start] + alt + sequence[codon_start + 3:]

                # Check if this substitution reduces overlaps
                new_overlaps = []
                if reference_sequences:
                    new_overlaps = check_kmer_overlap(trial, reference_sequences, k=word_size)

                new_blast = []
                if reference_file:
                    # Skip BLAST+ in inner loop for performance; rely on k-mer
                    pass

                if len(new_overlaps) < len(overlaps):
                    sequence = trial
                    fixed_any = True
                    break  # Accept this substitution and re-check

        if not fixed_any:
            # Could not reduce overlaps further — try eliminate_kmer_overlaps
            if reference_sequences:
                sequence = eliminate_kmer_overlaps(
                    sequence, protein, organism, overlaps, max_iterations=1,
                )
            break

    return sequence


def check_kmer_overlap(
    sequence: str,
    reference_sequences: list[str],
    k: int = 15,
) -> list[tuple[int, int, str]]:
    """Find k-mer overlaps between the sequence and references.

    Args:
        sequence: DNA sequence to check (uppercase).
        reference_sequences: List of reference DNA sequences.
        k: K-mer size (default 15).

    Returns:
        List of (start_pos, length, matching_kmer) tuples for each overlap found.
    """
    sequence = sequence.upper()
    overlaps: list[tuple[int, int, str]] = []

    if k < 1 or len(sequence) < k:
        return overlaps

    # Build k-mer index for the query sequence
    query_kmers: dict[str, int] = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if all(b in "ACGT" for b in kmer):
            query_kmers[kmer] = i

    # Check against each reference
    for ref_idx, ref in enumerate(reference_sequences):
        ref = ref.upper()
        if len(ref) < k:
            continue

        for j in range(len(ref) - k + 1):
            ref_kmer = ref[j:j + k]
            if ref_kmer in query_kmers:
                start = query_kmers[ref_kmer]
                # Extend the match as far as possible
                length = k
                q_pos = start + k
                r_pos = j + k
                while (q_pos < len(sequence) and r_pos < len(ref)
                       and sequence[q_pos] == ref[r_pos]):
                    length += 1
                    q_pos += 1
                    r_pos += 1

                overlaps.append((start, length, sequence[start:start + length]))

    # Deduplicate: keep longest overlap at each position
    if overlaps:
        overlaps = _deduplicate_overlaps(overlaps)

    return overlaps


def eliminate_kmer_overlaps(
    sequence: str,
    protein: str,
    organism: str,
    overlaps: list[tuple],
    max_iterations: int = 100,
) -> str:
    """Replace codons at overlap positions to eliminate k-mer overlaps.

    Args:
        sequence: Current DNA sequence.
        protein: Target protein sequence.
        organism: Target organism for codon preference.
        overlaps: List of (start, length, matching_kmer) tuples.
        max_iterations: Maximum number of codon substitutions to attempt.

    Returns:
        Modified DNA sequence with reduced overlaps.
    """
    sequence = sequence.upper()
    iterations = 0

    for start, length, _kmer in overlaps:
        if iterations >= max_iterations:
            break

        # Find codon indices that overlap with this region
        for pos in range(start, min(start + length, len(sequence))):
            codon_idx = pos // 3
            codon_start = codon_idx * 3
            if codon_start + 3 > len(sequence):
                continue

            current_codon = sequence[codon_start:codon_start + 3]
            aa = CODON_TABLE.get(current_codon)
            if aa is None or aa == "*":
                continue

            # Try synonymous alternatives that break the overlap
            alternatives = AA_TO_CODONS.get(aa, [])
            for alt in alternatives:
                if alt == current_codon:
                    continue
                # Check if this alt would break the k-mer at this position
                trial = sequence[:codon_start] + alt + sequence[codon_start + 3:]
                # Check if the overlap region still contains the kmer
                region = trial[max(0, start - 3):min(start + length + 3, len(trial))]
                if _kmer not in trial[max(0, start):min(start + length, len(trial))]:
                    sequence = trial
                    iterations += 1
                    break  # Move to next overlap

            if iterations >= max_iterations:
                break

    return sequence


# ── Internal helpers ─────────────────────────────────────────────────


def _run_blastn_short(
    sequence: str,
    reference_file: str,
    word_size: int = 15,
    e_value: float = 1.0,
    min_identity: float = 80.0,
) -> list[BLASTMatch]:
    """Run blastn-short against a reference FASTA file.

    Returns empty list if BLAST+ is not installed.
    """
    blastn_path = shutil.which("blastn")
    if blastn_path is None:
        logger.debug("blastn not found — skipping BLAST+ check")
        return []

    if not os.path.isfile(reference_file):
        logger.warning("Reference file not found: %s", reference_file)
        return []

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(f">query\n{sequence}\n")
        tmp_path = tmp.name

    try:
        cmd = [
            blastn_path,
            "-query", tmp_path,
            "-subject", reference_file,
            "-task", "blastn-short",
            "-word_size", str(word_size),
            "-evalue", str(e_value),
            "-outfmt", "6 qstart qend sstart send pident evalue sseqid",
            "-perc_identity", str(min_identity),
        ]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            logger.warning("blastn failed: %s", result.stderr[:200])
            return []

        matches: list[BLASTMatch] = []
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue
            try:
                matches.append(BLASTMatch(
                    query_start=int(fields[0]) - 1,  # 1-based to 0-based
                    query_end=int(fields[1]),
                    reference_start=int(fields[2]) - 1,
                    reference_end=int(fields[3]),
                    identity=float(fields[4]),
                    e_value=float(fields[5]),
                    reference_id=fields[6],
                ))
            except (ValueError, IndexError):
                continue

        return matches
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        logger.warning("blastn execution error: %s", e)
        return []
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


def check_blast_matches_cli(
    sequence: str,
    reference_file: str,
    word_size: int = 15,
    e_value: float = 1.0,
) -> list[dict]:
    """Use BLAST+ CLI to check for matches, returning empty list if not installed.

    This is a convenience wrapper around the internal ``_run_blastn_short``
    that returns plain dicts instead of :class:`BLASTMatch` objects, suitable
    for CLI / pipeline integration.

    Args:
        sequence: DNA sequence to check (uppercase).
        reference_file: Path to a FASTA file for BLAST+ search.
        word_size: BLAST word size.
        e_value: E-value threshold for BLAST+.

    Returns:
        List of dicts with keys: query_start, query_end, reference_start,
        reference_end, identity, e_value, reference_id. Returns empty list
        if BLAST+ is not installed or the reference file is not found.
    """
    matches = _run_blastn_short(sequence, reference_file, word_size, e_value)
    return [
        {
            "query_start": m.query_start,
            "query_end": m.query_end,
            "reference_start": m.reference_start,
            "reference_end": m.reference_end,
            "identity": m.identity,
            "e_value": m.e_value,
            "reference_id": m.reference_id,
        }
        for m in matches
    ]


def _deduplicate_overlaps(
    overlaps: list[tuple[int, int, str]],
) -> list[tuple[int, int, str]]:
    """Deduplicate overlaps, keeping the longest at each position range."""
    # Sort by start position, then by length descending
    overlaps.sort(key=lambda x: (x[0], -x[1]))

    deduped: list[tuple[int, int, str]] = []
    for overlap in overlaps:
        start, length, kmer = overlap
        # Check if this overlap is contained within an existing one
        contained = False
        for existing in deduped:
            e_start, e_length, _ = existing
            if start >= e_start and start + length <= e_start + e_length:
                contained = True
                break
        if not contained:
            deduped.append(overlap)

    return deduped
