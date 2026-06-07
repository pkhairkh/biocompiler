"""
Fuzzy matching utilities for biosecurity screening.

Includes reverse complement, Hamming distance, Levenshtein distance,
and fuzzy match search functions.
"""

from __future__ import annotations

# ─────────────────────────────────────────────────────────────────────────────
# Sequence utilities
# ─────────────────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    dna : str
        DNA sequence consisting of A, C, G, T (uppercase expected).

    Returns
    -------
    str
        The reverse complement.
    """
    return dna.translate(_COMPLEMENT)[::-1]


def _hamming_distance(s1: str, s2: str) -> int:
    """Compute the Hamming distance between two equal-length strings.

    Raises ``ValueError`` if the strings differ in length.
    """
    if len(s1) != len(s2):
        raise ValueError(
            f"Hamming distance requires equal-length strings, "
            f"got {len(s1)} and {len(s2)}"
        )
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def _levenshtein_distance(s1: str, s2: str) -> int:
    """Compute the Levenshtein (edit) distance between two strings.

    Uses dynamic programming with O(min(len(s1), len(s2))) memory.
    """
    if len(s1) < len(s2):
        s1, s2 = s2, s1

    # Now len(s1) >= len(s2)
    previous = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (c1 != c2)
            current.append(min(insertions, deletions, substitutions))
        previous = current

    return previous[-1]


def _fuzzy_match_hamming(
    sequence: str,
    motif: str,
    max_distance: int = 2,
) -> list[tuple[int, int, list[tuple[int, str, str]]]]:
    """Find fuzzy matches of *motif* in *sequence* using Hamming distance.

    Only returns matches with distance >= 1 (i.e. excludes exact matches)
    and distance <= *max_distance*.

    Parameters
    ----------
    sequence : str
        The sequence to search in.
    motif : str
        The motif to search for.
    max_distance : int
        Maximum Hamming distance to report (default 2).

    Returns
    -------
    list of (position, distance, substitutions)
        Each element is a tuple ``(pos, dist, subs)`` where *subs* is a
        list of ``(position_in_window, original_char, replacement_char)``
        tuples describing the substitutions.
    """
    mlen = len(motif)
    slen = len(sequence)
    if mlen == 0 or slen < mlen:
        return []

    results: list[tuple[int, int, list[tuple[int, str, str]]]] = []
    for i in range(slen - mlen + 1):
        window = sequence[i : i + mlen]
        dist = _hamming_distance(window, motif)
        if 1 <= dist <= max_distance:
            subs = [
                (j, motif[j], window[j])
                for j in range(mlen)
                if motif[j] != window[j]
            ]
            results.append((i, dist, subs))

    return results


def _fuzzy_match_edit_distance(
    sequence: str,
    motif: str,
    max_distance: int = 1,
) -> list[tuple[int, int]]:
    """Find fuzzy matches of *motif* in *sequence* using Levenshtein distance.

    Uses a sliding-window approach with windows of varying length around
    the motif length to catch insertions and deletions.

    Only returns matches with distance >= 1 and <= *max_distance*.

    Parameters
    ----------
    sequence : str
        The sequence to search in.
    motif : str
        The motif to search for.
    max_distance : int
        Maximum edit distance to report (default 1).

    Returns
    -------
    list of (position, distance)
        Each element is a tuple ``(pos, dist)``.
    """
    mlen = len(motif)
    slen = len(sequence)
    if mlen == 0 or slen == 0:
        return []

    results: list[tuple[int, int]] = []
    # Check windows of length mlen-1 through mlen+max_distance
    for window_len in range(max(1, mlen - max_distance), mlen + max_distance + 1):
        if window_len > slen:
            continue
        for i in range(slen - window_len + 1):
            window = sequence[i : i + window_len]
            dist = _levenshtein_distance(window, motif)
            if 1 <= dist <= max_distance:
                results.append((i, dist))

    # Deduplicate: keep the best distance per position
    best: dict[int, int] = {}
    for pos, dist in results:
        if pos not in best or dist < best[pos]:
            best[pos] = dist

    return sorted(best.items())
